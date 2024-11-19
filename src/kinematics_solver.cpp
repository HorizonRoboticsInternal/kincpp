#include "kinematics_solver.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace kincpp
{
KinematicsSolver::KinematicsSolver(const MatX& M, const MatX& S, const VecX& lower_joint_limits,
                                   const VecX& upper_joint_limits, const SolverType& solver_type)
    : M(M), S(S), lower_joint_limits(lower_joint_limits.array()),
      upper_joint_limits(upper_joint_limits.array()), solver_type(solver_type) {
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
    MatX Js = S;
    MatX T = MatX::Identity(4, 4);
    VecX sListTemp(S.col(0).size());
    for (int i = 1; i < q.size(); i++) {
        sListTemp << S.col(i - 1) * q(i - 1);
        T = T * MatrixExp6(VecToSE3(sListTemp));
        Js.col(i) = Adjoint(T) * S.col(i);
    }

    return Js;
}

/* Function: Gives the geometric Jacobian using spatial velocity.
 * Inputs: Joint configuration.
 * Returns: 6xn geometric Jacobian using spatial velocity.
 */
MatX KinematicsSolver::SpatialVelocityJacobian(const VecX& q) {
    MatX Js = VelocityTwistJacobian(q);
    // The top three rows correspond to angular and the bottom three rows are translational.
    // Here, we must swap them.
    Js.topRows(3).swap(Js.bottomRows(3));
    MatX velocity_tf = Eigen::Matrix<double, 6, 6>::Identity();
    Vec3 ee_pos = ForwardKinematics(q).block<3, 1>(0, 3);
    velocity_tf.block<3, 3>(0, 3) = -VecToSO3(ee_pos);
    return velocity_tf * Js;
}

std::pair<bool, VecX> KinematicsSolver::InverseKinematics(const MatX& desired_ee_tf, VecX& q_guess,
                                                          double position_tolerance,
                                                          double orientation_tolerance,
                                                          bool project_to_joint_limits,
                                                          bool use_pseudo_inverse,
                                                          int max_iterations) {
    switch (solver_type) {
        case NEWTON:
            return IK_Newton(desired_ee_tf, q_guess, position_tolerance, orientation_tolerance,
                             project_to_joint_limits, use_pseudo_inverse, max_iterations);
        case QP:
            return IK_QP(desired_ee_tf, q_guess, position_tolerance, orientation_tolerance,
                         max_iterations);
        default:
            throw std::runtime_error("Invalid solver type");
    }
}

std::pair<bool, VecX> KinematicsSolver::IK_Newton(const MatX& desired_ee_tf, VecX& q_guess,
                                                  double position_tolerance,
                                                  double orientation_tolerance,
                                                  bool project_to_joint_limits,
                                                  bool use_pseudo_inverse, int max_iterations) {
    int i = 0;
    MatX Tfk = ForwardKinematics(q_guess);
    MatX Tdiff = TransInv(Tfk) * desired_ee_tf;
    VecX Vs = Adjoint(Tfk) * SE3ToVec(MatrixLog6(Tdiff));
    Vec3 angular(Vs(0), Vs(1), Vs(2));
    Vec3 linear(Vs(3), Vs(4), Vs(5));

    bool err = (angular.norm() > orientation_tolerance || linear.norm() > position_tolerance);
    MatX Js;
    VecX curr_q = q_guess;
    while (err && i++ < max_iterations) {
        Js = VelocityTwistJacobian(curr_q);
        if (use_pseudo_inverse) {
            curr_q += Js.completeOrthogonalDecomposition().pseudoInverse() * Vs;
        }
        else {
            curr_q += Js.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vs);
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

std::pair<bool, VecX> KinematicsSolver::IK_QP(const MatX& desired_ee_tf, VecX& q_guess,
                                              double position_tolerance,
                                              double orientation_tolerance, int max_iterations) {
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

    // Binding for KinematicsSolver class
    /**
     * Class representing a kinematics solver for robotics.
     *
     * This solver supports both forward and inverse kinematics computations.
     */
    py::class_<kincpp::KinematicsSolver>(m, "KinematicsSolver")
        .def(py::init<const kincpp::MatX&, const kincpp::MatX&, const kincpp::VecX&,
                      const kincpp::VecX&, const kincpp::SolverType&>(),
             py::arg("M"), py::arg("S"), py::arg("lower_joint_limits"),
             py::arg("upper_joint_limits"), py::arg("solver_type") = kincpp::SolverType::QP,
             R"doc(
                Constructor for the KinematicsSolver class.

                Parameters:
                    M (np.ndarray[4, 4]): Home configuration of end effector.
                    S (np.ndarray[6, J]): Screw axes of joints when in home configuration.
                    lower_joint_limits (np.ndarray[J, 1]): Lower joint limits.
                    upper_joint_limits (np.ndarray[J, 1]): Upper joint limits.
                    solver_type (SolverType): Solver type to use for IK (default: QP).
            )doc")
        .def("forward_kinematics", &kincpp::KinematicsSolver::ForwardKinematics, py::arg("q"),
             R"doc(
                Forward kinematics function.

                Parameters:
                    q (np.ndarray[J, 1]): Joint positions for which to solve FK.

                Returns:
                    np.ndarray[4, 4]: The end effector transformation matrix.
            )doc")
        .def("inverse_kinematics", &kincpp::KinematicsSolver::InverseKinematics,
             py::arg("desired_ee_tf"), py::arg("q_guess"), py::arg("position_tolerance") = 1e-3,
             py::arg("orientation_tolerance") = 1e-3, py::arg("project_to_joint_limits") = true,
             py::arg("use_pseudo_inverse") = false, py::arg("max_iterations") = 20,
             R"doc(
                Inverse kinematics function.

                Parameters:
                    desired_ee_tf (np.ndarray[4, 4]): Desired end effector position and orientation.
                    q_guess (np.ndarray[J, 1]): Initial guess for joint positions.
                    position_tolerance (double): The end effector Cartesian position tolerance.
                    orientation_tolerance (double): The end effector orientation tolerance.
                    project_to_joint_limits (bool): Whether to respect joint limits during IK solve.
                        If using NEWTON, each update will be projected to the limits.
                        If using QP, this argument is ignored.
                        (default: True)
                    use_pseudo_inverse (bool): Whether to use pseudo-inverse for Jacobian inversion.
                        If using NEWTON, uses psuedo-inverse for Jacobian inversion.
                        If using QP, this argument is ignored.
                        (default: False)
                    max_iterations (int): Maximum number of iterations before solver quits (default: 20).

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
}
