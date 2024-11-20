#include "kinematics_solver.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
//#include <QuadProg++/QuadProg++.hh>
#undef solve
#include "eigquadprog.hpp"

namespace py = pybind11;

namespace kincpp
{
KinematicsSolver::KinematicsSolver(const MatX& M, const MatX& S, const VecX& lower_joint_limits,
                                   const VecX& upper_joint_limits)
    : M(M), S(S), lower_joint_limits(lower_joint_limits.array()),
      upper_joint_limits(upper_joint_limits.array()) {
}

KinematicsSolver::~KinematicsSolver() = default;

/* Function: Compute end effector frame (used for current spatial position calculation)
 * Inputs: Home configuration (position and orientation) of end-effector
 *		   The joint screw axes in the space frame when the manipulator
 *             is at the home position
 * 		   A list of joint coordinates.
 * Returns: Transformation matrix representing the end-effector frame when the joints are
 *				at the specified coordinates
 * Notes: FK means Forward Kinematics
 */
MatX KinematicsSolver::ForwardKinematics(const VecX& q) {
    MatX T = M;
    for (int i = (q.size() - 1); i > -1; i--) {
        T = MatrixExp6(VecToSE3(S.col(i) * q(i))) * T;
    }
    return T;
}

/* Function: Gives the geometric Jacobian using velocity twist.
 * Inputs: Joint configuration.
 * Returns: 6xn geometric Jacobian using velocity twist.
 */
MatX KinematicsSolver::VelocityTwistJacobian(const VecX& q) {
    MatX Jvt = S;
    Mat4 T = Mat4::Identity();
    VecX sListTemp(S.col(0).size());
    for (int i = 1; i < q.size(); i++) {
        sListTemp << S.col(i - 1) * q(i - 1);
        T = T * MatrixExp6(VecToSE3(sListTemp));
        Jvt.col(i) = Adjoint(T) * S.col(i);
    }

    return Jvt;
}

/* Function: Gives the geometric Jacobian using spatial velocity.
 * Inputs: Joint configuration.
 * Returns: 6xn geometric Jacobian using spatial velocity.
 */
MatX KinematicsSolver::SpatialVelocityJacobian(const VecX& q) {
    MatX Jvt = VelocityTwistJacobian(q);
    // The top three rows correspond to angular and the bottom three rows are translational.
    // Here, we must swap them.
    Jvt.topRows(3).swap(Jvt.bottomRows(3));
    MatX velocity_tf = Eigen::Matrix<double, 6, 6>::Identity();
    Vec3 ee_pos = ForwardKinematics(q).block<3, 1>(0, 3);
    velocity_tf.block<3, 3>(0, 3) = -VecToSO3(ee_pos);
    return velocity_tf * Jvt;
}

std::pair<bool, VecX> KinematicsSolver::InverseKinematics(const IKParams& params) {
    switch (params.solver_type) {
        case NEWTON:
            return IK_NM(params.common_params, params.newton_params);
        case QP:
            return IK_QP(params.common_params, params.qp_params);
        default:
            throw std::runtime_error("Invalid solver type");
    }
}

std::pair<bool, VecX> KinematicsSolver::IK_NM(const CommonParams& common_params,
                                              const NewtonParams& newton_params) {
    Mat4 desired_ee_tf = common_params.desired_ee_tf;
    VecX q_guess = common_params.q_guess;
    int max_iterations = common_params.iters;
    double position_tolerance = newton_params.position_tolerance;
    double orientation_tolerance = newton_params.orientation_tolerance;
    bool project_to_joint_limits = newton_params.project_to_joint_limits;
    bool use_pseudo_inverse = newton_params.use_pseudo_inverse;

    int i = 0;
    Mat4 Tfk = ForwardKinematics(q_guess);
    Mat4 Tdiff = TransInv(Tfk) * desired_ee_tf;
    VecX Vs = Adjoint(Tfk) * SE3ToVec(MatrixLog6(Tdiff));
    Vec3 angular(Vs(0), Vs(1), Vs(2));
    Vec3 linear(Vs(3), Vs(4), Vs(5));

    bool err = (angular.norm() > orientation_tolerance || linear.norm() > position_tolerance);
    MatX Jvt;
    VecX curr_q = q_guess;
    while (err && i++ < max_iterations) {
        Jvt = VelocityTwistJacobian(curr_q);
        if (use_pseudo_inverse) {
            curr_q += Jvt.completeOrthogonalDecomposition().pseudoInverse() * Vs;
        }
        else {
            curr_q += Jvt.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vs);
        }
        if (project_to_joint_limits) {
            curr_q = curr_q.array().max(lower_joint_limits).min(upper_joint_limits).matrix();
        }
        Tfk = ForwardKinematics(curr_q);
        Tdiff = TransInv(Tfk) * desired_ee_tf;
        Vs = Adjoint(Tfk) * SE3ToVec(MatrixLog6(Tdiff));

        // Compute error
        angular = Vec3(Vs(0), Vs(1), Vs(2));
        linear = Vec3(Vs(3), Vs(4), Vs(5));
        err = (angular.norm() > orientation_tolerance || linear.norm() > position_tolerance);
    }

    // If we weren't projecting to the joint limits during the Newton updates, let's post clip
    // the results by wrapping joint positions to [-pi, pi) and then clipping to the limits.
    if (!project_to_joint_limits) {
        // First, wrap joint positions to [-pi, pi)
        curr_q = (curr_q.array()).unaryExpr([](double angle) {
            return std::fmod(angle + PI, 2 * PI) - PI;
        });
        // Clamp joint positions to limits
        curr_q = curr_q.array().max(lower_joint_limits).min(upper_joint_limits).matrix();
    }

    return std::make_pair(!err, curr_q);
}

std::pair<bool, VecX> KinematicsSolver::IK_QP(const CommonParams& common_params,
                                              const QPParams& qp_params) {
    Mat4 desired_ee_tf = common_params.desired_ee_tf;
    VecX curr_q = common_params.q_guess;
    int iters = common_params.iters;

    double kj = 0.01;
    double ks = 1.00;
    double kq = 0.00;
    double km = 0.00;
    double ps = 0.00;
    double pi = 0.30;

    VecX lb(12);
    lb << -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -1e10, -1e10, -1e10, -1e10, -1e10, -1e10;
    VecX ub = -lb;

    for (int i = 0; i < iters; i++) {
        Mat4 curr_ee_tf = ForwardKinematics(curr_q);
        Vec6 angle_axis_error = AngleAxisCpp(curr_ee_tf, desired_ee_tf);
        //    double quad_error = 0.5 * angle_axis_error *
        MatX Js = SpatialVelocityJacobian(curr_q);

        int nj = curr_q.size();
        int nm = nj + 6;
        // Quadratic component of objective function
        MatX Q = MatX::Identity(nm, nm);

        // Joint velocity component of Q
        Q.block(0, 0, nj, nj) *= kj;

        // Slack component of Q
        Q.block(nj, nj, nj, nj) =
            ks * (1. / angle_axis_error.array().abs().sum()) * MatX::Identity(nj, nj);

        // Equality constraints
        MatX Aeq(Js.rows(), 2 * Js.cols());
        Aeq << Js, MatX::Identity(6, 6);
        Vec6 beq = angle_axis_error;

        std::optional<MatX> Ain_tmp;
        std::optional<VecX> bin_tmp;
        if (kq > 0.0) {
            Ain_tmp = MatX::Zero(nj + 6, nj + 6);
            bin_tmp = VecX::Zero(nj + 6);
            // Add joint limit velocity damper
            for (int j = 0; j < nj; ++j) {
                double ql0 = lower_joint_limits(j);  // Lower joint limit
                double ql1 = upper_joint_limits(j);  // Upper joint limit

                if (ql1 - curr_q[j] <= pi) {
                    bin_tmp.value()(j) = ((ql1 - curr_q[j]) - ps) / (pi - ps);
                    Ain_tmp.value()(j, j) = 1;
                }

                if (curr_q[j] - ql0 <= pi) {
                    bin_tmp.value()(j) = -((ql0 - curr_q[j]) + ps) / (pi - ps);
                    Ain_tmp.value()(j, j) = -1;
                }
            }
            bin_tmp.value().head(nj) *= (1.0 / kq);
        }
        // TODO: add manipulability maximization
        VecX c = VecX::Zero(nj + 6);

        MatX Ain;
        VecX bin;
        std::tie(Ain, bin) = AddBoundConstraints(Ain_tmp, bin_tmp, lb, ub);


//        quadprogpp::Vector<double> sol(12);
//        quadprogpp::Matrix<double> Q_qp = ConvertEigenMatToQPMat(Q);
//        quadprogpp::Vector<double> c_qp = ConvertEigenVecToQPVec(c);
//        quadprogpp::Matrix<double> Aeq_qp = ConvertEigenMatToQPMat(-Aeq.transpose());
//        quadprogpp::Vector<double> beq_qp = ConvertEigenVecToQPVec(-beq);
//        quadprogpp::Matrix<double> Ain_qp = ConvertEigenMatToQPMat(-Ain.transpose());
//        quadprogpp::Vector<double> bin_qp = ConvertEigenVecToQPVec(-bin);
//        double v = quadprogpp::solve_quadprog(Q_qp, c_qp, Aeq_qp, beq_qp, Ain_qp, bin_qp, sol);
//        std::cout << "quad matrix:\n" << Q << std::endl;
//        std::cout << "quad vector:\n" << c << std::endl;
//        std::cout << "inequality matrix:\n" << -Ain.transpose() << std::endl;
//        std::cout << "inequality vector:\n" << -bin << std::endl;
//        std::cout << "equality matrix:\n" << -Aeq.transpose() << std::endl;
//        std::cout << "equality vector:\n" << -beq << std::endl;

        // NOTE:
        // https://www.labri.fr/perso/guenneba/code/QuadProg/eiquadprog.hpp

        VecX sol(12);
        double v = Eigen::solve_quadprog(Q, c, -Aeq.transpose(), beq, -Ain.transpose(), bin, sol);
        for (int k = 0; k < 6; ++k) {
            std::cout << sol[k] << " ";

        }
        std::cout << std::endl;
        std::cout << v << std::endl;
        exit(0);

//        curr_q += sol
    }
        return std::make_pair(false, Eigen::Vector<double, 6>::Zero());
    }
}  // namespace kincpp

// Pybind11 binding function
PYBIND11_MODULE(kincpp, m) {
    // Binding for SolverType enum
    /**
     * Enum representing the type of solver to use for inverse kinematics.
     *
     * Values:
     * - NEWTON (0): Use Newton's method for solving IK.
     * - QP (1): Use Quadratic Programming (QP) for solving IK.
     */
    py::enum_<kincpp::SolverType>(m, "SolverType")
        .value("NEWTON", kincpp::SolverType::NEWTON, "Newton's method for inverse kinematics.")
        .value("QP", kincpp::SolverType::QP, "Quadratic Programming (QP) for inverse kinematics.")
        .export_values();

    /**
     * kincpp::CommonParams struct.
     *
     * Parameters:
     *     desired_ee_tf (np.ndarray[4, 4]): Desired end-effector transformation matrix.
     *     q_guess (np.ndarray[J, 1]): Initial guess for joint positions.
     *     iters (int): Number of iterations for solving the IK problem.
     */
    py::class_<kincpp::CommonParams>(m, "CommonParams")
        .def(py::init<>(),
             R"doc(
                Common parameters for inverse kinematics.

                Parameters:
                    desired_ee_tf (np.ndarray[4, 4]): Desired end-effector transformation matrix.
                    q_guess (np.ndarray[J, 1]): Initial guess for joint positions.
                    iters (int): Number of iterations for solving the IK problem.
             )doc")
        .def_readwrite("desired_ee_tf", &kincpp::CommonParams::desired_ee_tf,
                       "Desired end-effector transformation matrix.")
        .def_readwrite("q_guess", &kincpp::CommonParams::q_guess,
                       "Initial guess for joint positions.")
        .def_readwrite("iters", &kincpp::CommonParams::iters,
                       "Number of iterations for solving the IK problem.");

    /**
     * kincpp::NewtonParams struct.
     *
     * Parameters:
     *     position_tolerance (double): Tolerance for the end-effector Cartesian position.
     *     orientation_tolerance (double): Tolerance for the end-effector orientation.
     *     project_to_joint_limits (bool): Whether to respect joint limits during IK solve.
     *         If True, joint updates are projected to limits (default: False).
     *     use_pseudo_inverse (bool): Whether to use pseudo-inverse for Jacobian inversion.
     *         If True, pseudo-inverse is used; ignored for non-Newton solvers (default: False).
     */
    py::class_<kincpp::NewtonParams>(m, "NewtonParams")
        .def(py::init<>(),
             R"doc(
                Parameters specific to Newton's method for inverse kinematics.

                Parameters:
                    position_tolerance (double): Tolerance for the end-effector Cartesian position.
                    orientation_tolerance (double): Tolerance for the end-effector orientation.
                    project_to_joint_limits (bool): Whether to respect joint limits during IK solve.
                        If True, joint updates are projected to limits (default: False).
                    use_pseudo_inverse (bool): Whether to use pseudo-inverse for Jacobian inversion.
                        If True, pseudo-inverse is used; ignored for non-Newton solvers (default: False).
             )doc")
        .def_readwrite("position_tolerance", &kincpp::NewtonParams::position_tolerance,
                       "Tolerance for the end-effector Cartesian position.")
        .def_readwrite("orientation_tolerance", &kincpp::NewtonParams::orientation_tolerance,
                       "Tolerance for the end-effector orientation.")
        .def_readwrite("project_to_joint_limits", &kincpp::NewtonParams::project_to_joint_limits,
                       "Whether to respect joint limits during IK solve.")
        .def_readwrite("use_pseudo_inverse", &kincpp::NewtonParams::use_pseudo_inverse,
                       "Whether to use pseudo-inverse for Jacobian inversion.");

    /**
     * kincpp::IKParams struct.
     *
     * Parameters:
     *     common_params (kincpp::CommonParams): Common parameters for inverse kinematics.
     *     newton_params (kincpp::NewtonParams): Parameters specific to Newton's method.
     *     solver_type (kincpp::SolverType): The solver type to use for IK (default: NEWTON).
     */
    py::class_<kincpp::IKParams>(m, "IKParams")
        .def(py::init<>(),
             R"doc(
                Parameters for inverse kinematics, including solver selection.

                Parameters:
                    common_params (kincpp::CommonParams): Common parameters for inverse kinematics.
                    newton_params (kincpp::NewtonParams): Parameters specific to Newton's method.
                    solver_type (kincpp::SolverType): The solver type to use for IK (default: NEWTON).
             )doc")
        .def_readwrite("common_params", &kincpp::IKParams::common_params,
                       "Common parameters for inverse kinematics.")
        .def_readwrite("newton_params", &kincpp::IKParams::newton_params,
                       "Parameters specific to Newton's method.")
        .def_readwrite("solver_type", &kincpp::IKParams::solver_type,
                       "The solver type to use for IK.");
    // Binding for KinematicsSolver class
    /**
     * Class representing a kinematics solver for robotics.
     *
     * This solver supports both forward and inverse kinematics computations.
     */
    py::class_<kincpp::KinematicsSolver>(m, "KinematicsSolver")
        .def(py::init<const kincpp::MatX&, const kincpp::MatX&, const kincpp::VecX&,
                      const kincpp::VecX&>(),
             py::arg("M"), py::arg("S"), py::arg("lower_joint_limits"),
             py::arg("upper_joint_limits"),
             R"doc(
                Constructor for the KinematicsSolver class.

                Parameters:
                    M (np.ndarray[4, 4]): Home configuration of end effector.
                    S (np.ndarray[6, J]): Screw axes of joints when in home configuration.
                    lower_joint_limits (np.ndarray[J, 1]): Lower joint limits.
                    upper_joint_limits (np.ndarray[J, 1]): Upper joint limits.
            )doc")
        .def("forward_kinematics", &kincpp::KinematicsSolver::ForwardKinematics, py::arg("q"),
             R"doc(
                Forward kinematics function.

                Parameters:
                    q (np.ndarray[J, 1]): Joint positions for which to solve FK.

                Returns:
                    np.ndarray[4, 4]: The end effector transformation matrix.
            )doc")
        .def("inverse_kinematics", &kincpp::KinematicsSolver::InverseKinematics, py::arg("params"),
             R"doc(
                Inverse kinematics function.

                Parameters:
                    params: IK params.

                Returns:
                    tuple:
                        - bool: Whether the IK solver succeeded or not.
                        - np.ndarray[J, 1]: A joint angle solution corresponding to the desired end-effector pose.
                          If the IK solver fails, this result is undefined.
            )doc")
        .def(
            "velocity_twist_jacobian",
            [](kincpp::KinematicsSolver& self, const kincpp::VecX& q) {
                // Call the original VelocityTwistJacobian function
                kincpp::MatX Js = self.VelocityTwistJacobian(q);

                // Perform the swap of top three rows with bottom three rows
                Js.topRows(3).swap(Js.bottomRows(3));

                return Js;  // Return the swapped matrix
            },
            py::arg("q"),
            R"doc(
                Compute the geometric Jacobian using velocity twist.

                Parameters:
                    q (np.ndarray[J, 1]): Joint configuration.

                Returns:
                    np.ndarray[6, J]: Geometric Jacobian using velocity twist.
             )doc")
        .def("spatial_velocity_jacobian", &kincpp::KinematicsSolver::SpatialVelocityJacobian,
             py::arg("q"),
             R"doc(
                Compute the geometric Jacobian using spatial velocity.

                Parameters:
                    q (np.ndarray[J, 1]): Joint configuration.

                Returns:
                    np.ndarray[6, J]: Geometric Jacobian using spatial velocity.
             )doc");

    /**
     * Computes the 6D pose error between two transformation matrices.
     *
     * Parameters:
     *     T (np.ndarray[4, 4]): Current transformation matrix.
     *     Td (np.ndarray[4, 4]): Desired transformation matrix.
     *
     * Returns:
     *     np.ndarray[6, 1]: A 6D vector where:
     *         - The first three elements represent translational error.
     *         - The last three elements represent rotational error in axis-angle form.
     */
    m.def("angle_axis_cpp", &kincpp::AngleAxisCpp,
          R"doc(
            Computes the 6D pose error between two transformation matrices.

            Parameters:
                T (np.ndarray[4, 4]): Current transformation matrix.
                Td (np.ndarray[4, 4]): Desired transformation matrix.

            Returns:
                np.ndarray[6, 1]: A 6D vector where:
                    - The first three elements represent translational error.
                    - The last three elements represent rotational error in axis-angle form.
          )doc");
}
