/*
 * Adapted from https://github.com/Le0nX/ModernRoboticsCpp
 */
#include <Eigen/Dense>
#include <cmath>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <vector>

namespace py = pybind11;

namespace kincpp
{

/* Function: Find if the value is negligible enough to consider 0
 * Inputs: value to be checked as a double
 * Returns: Boolean of true-ignore or false-can't ignore
 */
bool NearZero(const double val) {
    return (std::abs(val) < .000001);
}

/* Function: Returns the skew symmetric matrix representation of an angular velocity vector
 * Input: Eigen::Vector3d 3x1 angular velocity vector
 * Returns: Eigen::MatrixXd 3x3 skew symmetric matrix
 */
Eigen::Matrix3d VecToSO3(const Eigen::Vector3d& omg) {
    Eigen::Matrix3d m_ret;
    m_ret << 0, -omg(2), omg(1), omg(2), 0, -omg(0), -omg(1), omg(0), 0;
    return m_ret;
}

/* Function: Returns a normalized version of the input vector
 * Input: Eigen::MatrixXd
 * Output: Eigen::MatrixXd
 * Note: MatrixXd is used instead of VectorXd for the case of row vectors
 * 		Requires a copy
 *		Useful because of the MatrixXd casting
 */
Eigen::MatrixXd Normalize(Eigen::MatrixXd V) {
    V.normalize();
    return V;
}

/* Function: Returns angular velocity vector represented by the skew symmetric matrix
 * Inputs: Eigen::MatrixXd 3x3 skew symmetric matrix
 * Returns: Eigen::Vector3d 3x1 angular velocity
 */
Eigen::Vector3d SO3ToVec(const Eigen::MatrixXd& so3mat) {
    Eigen::Vector3d v_ret;
    v_ret << so3mat(2, 1), so3mat(0, 2), so3mat(1, 0);
    return v_ret;
}

/* Function: Translates an exponential rotation into it's individual components
 * Inputs: Exponential rotation (rotation matrix in terms of a rotation axis
 *				and the angle of rotation)
 * Returns: The axis and angle of rotation as [x, y, z, theta]
 */
Eigen::Vector4d AxisAng3(const Eigen::Vector3d& expc3) {
    Eigen::Vector4d v_ret;
    v_ret << Normalize(expc3), expc3.norm();
    return v_ret;
}

/* Function: Translates an exponential rotation into a rotation matrix
 * Inputs: exponenential representation of a rotation
 * Returns: Rotation matrix
 */
Eigen::Matrix3d MatrixExp3(const Eigen::Matrix3d& so3mat) {
    Eigen::Vector3d omgtheta = SO3ToVec(so3mat);

    Eigen::Matrix3d m_ret = Eigen::Matrix3d::Identity();
    if (NearZero(so3mat.norm())) {
        return m_ret;
    }
    else {
        double theta = (AxisAng3(omgtheta))(3);
        Eigen::Matrix3d omgmat = so3mat * (1 / theta);
        return m_ret + std::sin(theta) * omgmat + ((1 - std::cos(theta)) * (omgmat * omgmat));
    }
}

/* Function: Computes the matrix logarithm of a rotation matrix
 * Inputs: Rotation matrix
 * Returns: matrix logarithm of a rotation
 */
Eigen::Matrix3d MatrixLog3(const Eigen::Matrix3d& R) {
    double acosinput = (R.trace() - 1) / 2.0;
    Eigen::MatrixXd m_ret = Eigen::MatrixXd::Zero(3, 3);
    if (acosinput >= 1)
        return m_ret;
    else if (acosinput <= -1) {
        Eigen::Vector3d omg;
        if (!NearZero(1 + R(2, 2)))
            omg = (1.0 / std::sqrt(2 * (1 + R(2, 2)))) *
                  Eigen::Vector3d(R(0, 2), R(1, 2), 1 + R(2, 2));
        else if (!NearZero(1 + R(1, 1)))
            omg = (1.0 / std::sqrt(2 * (1 + R(1, 1)))) *
                  Eigen::Vector3d(R(0, 1), 1 + R(1, 1), R(2, 1));
        else
            omg = (1.0 / std::sqrt(2 * (1 + R(0, 0)))) *
                  Eigen::Vector3d(1 + R(0, 0), R(1, 0), R(2, 0));
        m_ret = VecToSO3(M_PI * omg);
        return m_ret;
    }
    else {
        double theta = std::acos(acosinput);
        m_ret = theta / 2.0 / sin(theta) * (R - R.transpose());
        return m_ret;
    }
}

/* Function: Combines a rotation matrix and position vector into a single
 * 				Special Euclidian Group (SE3) homogeneous transformation matrix
 * Inputs: Rotation Matrix (R), Position Vector (p)
 * Returns: Matrix of T = [ [R, p],
 *						    [0, 1] ]
 */
Eigen::MatrixXd RpToTrans(const Eigen::Matrix3d& R, const Eigen::Vector3d& p) {
    Eigen::MatrixXd m_ret(4, 4);
    m_ret << R, p, 0, 0, 0, 1;
    return m_ret;
}

/* Function: Separates the rotation matrix and position vector from
 *				the transfomation matrix representation
 * Inputs: Homogeneous transformation matrix
 * Returns: std::vector of [rotation matrix, position vector]
 */
std::vector<Eigen::MatrixXd> TransToRp(const Eigen::MatrixXd& T) {
    std::vector<Eigen::MatrixXd> Rp_ret;
    Eigen::Matrix3d R_ret;
    // Get top left 3x3 corner
    R_ret = T.block<3, 3>(0, 0);

    Eigen::Vector3d p_ret(T(0, 3), T(1, 3), T(2, 3));

    Rp_ret.emplace_back(R_ret);
    Rp_ret.emplace_back(p_ret);

    return Rp_ret;
}

/* Function: Translates a spatial velocity vector into a transformation matrix
 * Inputs: Spatial velocity vector [angular velocity, linear velocity]
 * Returns: Transformation matrix
 */
Eigen::MatrixXd VecToSE3(const Eigen::VectorXd& V) {
    // Separate angular (exponential representation) and linear velocities
    Eigen::Vector3d exp(V(0), V(1), V(2));
    Eigen::Vector3d linear(V(3), V(4), V(5));

    // Fill in values to the appropriate parts of the transformation matrix
    Eigen::MatrixXd m_ret(4, 4);
    m_ret << VecToSO3(exp), linear, 0, 0, 0, 0;

    return m_ret;
}

/* Function: Translates a transformation matrix into a spatial velocity vector
 * Inputs: Transformation matrix
 * Returns: Spatial velocity vector [angular velocity, linear velocity]
 */
Eigen::VectorXd SE3ToVec(const Eigen::MatrixXd& T) {
    Eigen::VectorXd m_ret(6);
    m_ret << T(2, 1), T(0, 2), T(1, 0), T(0, 3), T(1, 3), T(2, 3);

    return m_ret;
}

/* Function: Provides the adjoint representation of a transformation matrix
 *			 Used to change the frame of reference for spatial velocity vectors
 * Inputs: 4x4 Transformation matrix SE(3)
 * Returns: 6x6 Adjoint Representation of the matrix
 */
Eigen::MatrixXd Adjoint(const Eigen::MatrixXd& T) {
    std::vector<Eigen::MatrixXd> R = TransToRp(T);
    Eigen::MatrixXd ad_ret(6, 6);
    ad_ret = Eigen::MatrixXd::Zero(6, 6);
    Eigen::MatrixXd zeroes = Eigen::MatrixXd::Zero(3, 3);
    ad_ret << R[0], zeroes, VecToSO3(R[1]) * R[0], R[0];
    return ad_ret;
}

/* Function: Rotation expanded for screw axis
 * Inputs: se3 matrix representation of exponential coordinates (transformation matrix)
 * Returns: 6x6 Matrix representing the rotation
 */
Eigen::MatrixXd MatrixExp6(const Eigen::MatrixXd& se3mat) {
    // Extract the angular velocity vector from the transformation matrix
    Eigen::Matrix3d se3mat_cut = se3mat.block<3, 3>(0, 0);
    Eigen::Vector3d omgtheta = SO3ToVec(se3mat_cut);

    Eigen::MatrixXd m_ret(4, 4);

    if (NearZero(omgtheta.norm())) {
        // Reuse previous variables that have our required size
        se3mat_cut = Eigen::MatrixXd::Identity(3, 3);
        omgtheta << se3mat(0, 3), se3mat(1, 3), se3mat(2, 3);
        m_ret << se3mat_cut, omgtheta, 0, 0, 0, 1;
        return m_ret;
    }
    else {
        double theta = (AxisAng3(omgtheta))(3);
        Eigen::Matrix3d omgmat = se3mat.block<3, 3>(0, 0) / theta;
        Eigen::Matrix3d expExpand = Eigen::MatrixXd::Identity(3, 3) * theta +
                                    (1 - std::cos(theta)) * omgmat +
                                    ((theta - std::sin(theta)) * (omgmat * omgmat));
        Eigen::Vector3d linear(se3mat(0, 3), se3mat(1, 3), se3mat(2, 3));
        Eigen::Vector3d GThetaV = (expExpand * linear) / theta;
        m_ret << MatrixExp3(se3mat_cut), GThetaV, 0, 0, 0, 1;
        return m_ret;
    }
}

Eigen::MatrixXd MatrixLog6(const Eigen::MatrixXd& T) {
    Eigen::MatrixXd m_ret(4, 4);
    auto rp = kincpp::TransToRp(T);
    Eigen::Matrix3d omgmat = MatrixLog3(rp.at(0));
    Eigen::Matrix3d zeros3d = Eigen::Matrix3d::Zero(3, 3);
    if (NearZero(omgmat.norm())) {
        m_ret << zeros3d, rp.at(1), 0, 0, 0, 0;
    }
    else {
        double theta = std::acos((rp.at(0).trace() - 1) / 2.0);
        Eigen::Matrix3d logExpand1 = Eigen::MatrixXd::Identity(3, 3) - omgmat / 2.0;
        Eigen::Matrix3d logExpand2 =
            (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2) * omgmat * omgmat / theta;
        Eigen::Matrix3d logExpand = logExpand1 + logExpand2;
        m_ret << omgmat, logExpand * rp.at(1), 0, 0, 0, 0;
    }
    return m_ret;
}

/* Function: Compute end effector frame (used for current spatial position calculation)
 * Inputs: Home configuration (position and orientation) of end-effector
 *		   The joint screw axes in the space frame when the manipulator
 *             is at the home position
 * 		   A list of joint coordinates.
 * Returns: Transformation matrix representing the end-effector frame when the joints are
 *				at the specified coordinates
 * Notes: FK means Forward Kinematics
 */
Eigen::MatrixXd FKinSpace(const Eigen::MatrixXd& M, const Eigen::MatrixXd& Slist,
                          const Eigen::VectorXd& thetaList) {
    Eigen::MatrixXd T = M;
    for (int i = (thetaList.size() - 1); i > -1; i--) {
        T = MatrixExp6(VecToSE3(Slist.col(i) * thetaList(i))) * T;
    }
    return T;
}

/* Function: Gives the space Jacobian
 * Inputs: Screw axis in home position, joint configuration
 * Returns: 6xn Spatial Jacobian
 */
Eigen::MatrixXd JacobianSpace(const Eigen::MatrixXd& Slist, const Eigen::MatrixXd& thetaList) {
    Eigen::MatrixXd Js = Slist;
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(4, 4);
    Eigen::VectorXd sListTemp(Slist.col(0).size());
    for (int i = 1; i < thetaList.size(); i++) {
        sListTemp << Slist.col(i - 1) * thetaList(i - 1);
        T = T * MatrixExp6(VecToSE3(sListTemp));
        Js.col(i) = Adjoint(T) * Slist.col(i);
    }

    return Js;
}

/*
 * Function: Gives the body Jacobian
 * Inputs: Screw axis in BODY position, joint configuration
 * Returns: 6xn Bobdy Jacobian
 */
Eigen::MatrixXd JacobianBody(const Eigen::MatrixXd& Blist, const Eigen::MatrixXd& thetaList) {
    Eigen::MatrixXd Jb = Blist;
    Eigen::MatrixXd T = Eigen::MatrixXd::Identity(4, 4);
    Eigen::VectorXd bListTemp(Blist.col(0).size());
    for (int i = thetaList.size() - 2; i >= 0; i--) {
        bListTemp << Blist.col(i + 1) * thetaList(i + 1);
        T = T * MatrixExp6(VecToSE3(-1 * bListTemp));
        Jb.col(i) = Adjoint(T) * Blist.col(i);
    }
    return Jb;
}

Eigen::MatrixXd TransInv(const Eigen::MatrixXd& transform) {
    auto rp = kincpp::TransToRp(transform);
    auto Rt = rp.at(0).transpose();
    auto t = -(Rt * rp.at(1));
    Eigen::MatrixXd inv(4, 4);
    inv = Eigen::MatrixXd::Zero(4, 4);
    inv.block(0, 0, 3, 3) = Rt;
    inv.block(0, 3, 3, 1) = t;
    inv(3, 3) = 1;
    return inv;
}

Eigen::MatrixXd RotInv(const Eigen::MatrixXd& rotMatrix) {
    return rotMatrix.transpose();
}

Eigen::VectorXd ScrewToAxis(Eigen::Vector3d q, Eigen::Vector3d s, double h) {
    Eigen::VectorXd axis(6);
    axis.segment(0, 3) = s;
    axis.segment(3, 3) = q.cross(s) + (h * s);
    return axis;
}

Eigen::VectorXd AxisAng6(const Eigen::VectorXd& expc6) {
    Eigen::VectorXd v_ret(7);
    double theta = Eigen::Vector3d(expc6(0), expc6(1), expc6(2)).norm();
    if (NearZero(theta))
        theta = Eigen::Vector3d(expc6(3), expc6(4), expc6(5)).norm();
    v_ret << expc6 / theta, theta;
    return v_ret;
}

Eigen::MatrixXd ProjectToSO3(const Eigen::MatrixXd& M) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();
    if (R.determinant() < 0)
        // In this case the result may be far from M; reverse sign of 3rd column
        R.col(2) *= -1;
    return R;
}

Eigen::MatrixXd ProjectToSE3(const Eigen::MatrixXd& M) {
    Eigen::Matrix3d R = M.block<3, 3>(0, 0);
    Eigen::Vector3d t = M.block<3, 1>(0, 3);
    Eigen::MatrixXd T = RpToTrans(ProjectToSO3(R), t);
    return T;
}

double DistanceToSO3(const Eigen::Matrix3d& M) {
    if (M.determinant() > 0)
        return (M.transpose() * M - Eigen::Matrix3d::Identity()).norm();
    else
        return 1.0e9;
}

double DistanceToSE3(const Eigen::Matrix4d& T) {
    Eigen::Matrix3d matR = T.block<3, 3>(0, 0);
    if (matR.determinant() > 0) {
        Eigen::Matrix4d m_ret;
        m_ret << matR.transpose() * matR, Eigen::Vector3d::Zero(3), T.row(3);
        m_ret = m_ret - Eigen::Matrix4d::Identity();
        return m_ret.norm();
    }
    else
        return 1.0e9;
}

bool TestIfSO3(const Eigen::Matrix3d& M) {
    return std::abs(DistanceToSO3(M)) < 1e-3;
}

bool TestIfSE3(const Eigen::Matrix4d& T) {
    return std::abs(DistanceToSE3(T)) < 1e-3;
}

std::pair<bool, Eigen::VectorXd>
IKinSpace(const Eigen::MatrixXd& M, const Eigen::MatrixXd& Slist, const Eigen::MatrixXd& T,
          Eigen::VectorXd& thetalist, const Eigen::ArrayXd& lower_limits,
          const Eigen::ArrayXd& upper_limits, bool project_to_joint_limits, bool use_pseudo_inverse,
          double position_tolerance, double orientation_tolerance, int max_iterations = 20) {
    int i = 0;
    Eigen::MatrixXd Tfk = FKinSpace(M, Slist, thetalist);
    Eigen::MatrixXd Tdiff = TransInv(Tfk) * T;
    Eigen::VectorXd Vs = Adjoint(Tfk) * SE3ToVec(MatrixLog6(Tdiff));
    Eigen::Vector3d angular(Vs(0), Vs(1), Vs(2));
    Eigen::Vector3d linear(Vs(3), Vs(4), Vs(5));

    bool err = (angular.norm() > orientation_tolerance || linear.norm() > position_tolerance);
    Eigen::MatrixXd Js;

    while (err && i++ < max_iterations) {
        Js = JacobianSpace(Slist, thetalist);
        if (use_pseudo_inverse) {
            thetalist += Js.completeOrthogonalDecomposition().pseudoInverse() * Vs;
        }
        else {
            thetalist += Js.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vs);
        }
        if (project_to_joint_limits) {
            thetalist = thetalist.array().max(lower_limits).min(upper_limits).matrix();
        }
        Tfk = FKinSpace(M, Slist, thetalist);
        Tdiff = TransInv(Tfk) * T;
        Vs = Adjoint(Tfk) * SE3ToVec(MatrixLog6(Tdiff));

        // Compute error
        angular = Eigen::Vector3d(Vs(0), Vs(1), Vs(2));
        linear = Eigen::Vector3d(Vs(3), Vs(4), Vs(5));
        err = (angular.norm() > orientation_tolerance || linear.norm() > position_tolerance);
    }

    return std::make_pair(!err, thetalist);
}
}  // namespace kincpp

// Pybind11 binding function
PYBIND11_MODULE(kincpp, m) {
    m.def("forward", &kincpp::FKinSpace,
          "Forward kinematics function. \n\n"
          "Parameters:\n"
          "    M (np.ndarray[4, 4]): Home configuration of end effector.\n"
          "    S (np.ndarray[6, J]): Screw axes of joints when in home configuration.\n"
          "    joint_positions (np.ndarray[J, 1]): Joint positions for which to solve FK for.\n"
          "Returns:\n"
          "    np.ndarray[4, 4]: The end effector transformation matrix.",
          py::arg("M"), py::arg("S"), py::arg("joint_positions"));

    m.def("inverse", &kincpp::IKinSpace,
          "Inverse kinematics function. \n\n"
          "Parameters:\n"
          "    M (np.ndarray[4, 4]): Home configuration of end effector.\n"
          "    S (np.ndarray[6, J]): Screw axes of joints when in home configuration.\n"
          "    T (np.ndarray[4, 4]): Desired end effector position and orientation.\n"
          "    joint_position_guess (np.ndarray[J, 1]): Joint positions initial guess for IK "
          "solver.\n"
          "    joint_lower_limits (np.ndarray[J, 1]): Joint lower limits. Only used if "
          "project_to_joint_limits is True.\n"
          "    joint_upper_limits (np.ndarray[J, 1]): Joint upper limits. Only used if "
          "project_to_joint_limits is True.\n"
          "    project_to_joint_limits (bool): Whether to respect joint limits during IK solve.\n"
          "    use_pseudo_inverse (bool): Whether to use pseudo inverse as the inverse Jacobian "
          "during IK solve.\n"
          "    position_tolerance (double): The end effector Cartesian position tolerance.\n"
          "    orientation_tolerance (double): The end effector orientation tolerance.\n"
          "    max_iterations (int) = 20: Maximum number of iterations before solver quits.\n"
          "Returns:\n"
          "    bool: Whether the IK solver succeeded or not\n"
          "    np.ndarray[J, 1]: A joint angle solution corresponding to T.\n"
          "                      If IK solver failed, result is undefined. \n",
          py::arg("M"), py::arg("S"), py::arg("T"), py::arg("joint_position_guess"),
          py::arg("joint_lower_limits"), py::arg("joint_upper_limits"),
          py::arg("project_to_joint_limits"), py::arg("use_pseudo_inverse"),
          py::arg("position_tolerance"), py::arg("orientation_tolerance"),
          py::arg("max_iterations") = 20);
}
