#include "utils.h"

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
 * Input: Vec3 3x1 angular velocity vector
 * Returns: MatX 3x3 skew symmetric matrix
 */
Mat3 VecToSO3(const Vec3& omg) {
    Mat3 m_ret;
    m_ret << 0, -omg(2), omg(1), omg(2), 0, -omg(0), -omg(1), omg(0), 0;
    return m_ret;
}

/* Function: Returns a normalized version of the input vector
 * Input: MatX
 * Output: MatX
 * Note: MatrixXd is used instead of VectorXd for the case of row vectors
 * 		Requires a copy
 *		Useful because of the MatrixXd casting
 */
MatX Normalize(MatX V) {
    V.normalize();
    return V;
}

/* Function: Returns angular velocity vector represented by the skew symmetric matrix
 * Inputs: MatX 3x3 skew symmetric matrix
 * Returns: Vec3 3x1 angular velocity
 */
Vec3 SO3ToVec(const MatX& so3mat) {
    Vec3 v_ret;
    v_ret << so3mat(2, 1), so3mat(0, 2), so3mat(1, 0);
    return v_ret;
}

/* Function: Translates an exponential rotation into it's individual components
 * Inputs: Exponential rotation (rotation matrix in terms of a rotation axis
 *				and the angle of rotation)
 * Returns: The axis and angle of rotation as [x, y, z, theta]
 */
Eigen::Vector4d AxisAng3(const Vec3& expc3) {
    Eigen::Vector4d v_ret;
    v_ret << Normalize(expc3), expc3.norm();
    return v_ret;
}

/* Function: Translates an exponential rotation into a rotation matrix
 * Inputs: exponenential representation of a rotation
 * Returns: Rotation matrix
 */
Mat3 MatrixExp3(const Mat3& so3mat) {
    Vec3 omgtheta = SO3ToVec(so3mat);

    Mat3 m_ret = Mat3::Identity();
    if (NearZero(so3mat.norm())) {
        return m_ret;
    }
    else {
        double theta = (AxisAng3(omgtheta))(3);
        Mat3 omgmat = so3mat * (1 / theta);
        return m_ret + std::sin(theta) * omgmat + ((1 - std::cos(theta)) * (omgmat * omgmat));
    }
}

/* Function: Computes the matrix logarithm of a rotation matrix
 * Inputs: Rotation matrix
 * Returns: matrix logarithm of a rotation
 */
Mat3 MatrixLog3(const Mat3& R) {
    double acosinput = (R.trace() - 1) / 2.0;
    MatX m_ret = MatX::Zero(3, 3);
    if (acosinput >= 1)
        return m_ret;
    else if (acosinput <= -1) {
        Vec3 omg;
        if (!NearZero(1 + R(2, 2)))
            omg = (1.0 / std::sqrt(2 * (1 + R(2, 2)))) * Vec3(R(0, 2), R(1, 2), 1 + R(2, 2));
        else if (!NearZero(1 + R(1, 1)))
            omg = (1.0 / std::sqrt(2 * (1 + R(1, 1)))) * Vec3(R(0, 1), 1 + R(1, 1), R(2, 1));
        else
            omg = (1.0 / std::sqrt(2 * (1 + R(0, 0)))) * Vec3(1 + R(0, 0), R(1, 0), R(2, 0));
        m_ret = VecToSO3(M_PI * omg);
        return m_ret;
    }
    else {
        double theta = std::acos(acosinput);
        m_ret = theta / 2.0 / sin(theta) * (R - R.transpose());
        return m_ret;
    }
}

/* Function: Separates the rotation matrix and position vector from
 *				the transfomation matrix representation
 * Inputs: Homogeneous transformation matrix
 * Returns: std::vector of [rotation matrix, position vector]
 */
std::vector<MatX> TransToRp(const MatX& T) {
    std::vector<MatX> Rp_ret;
    Mat3 R_ret;
    // Get top left 3x3 corner
    R_ret = T.block<3, 3>(0, 0);

    Vec3 p_ret(T(0, 3), T(1, 3), T(2, 3));

    Rp_ret.emplace_back(R_ret);
    Rp_ret.emplace_back(p_ret);

    return Rp_ret;
}

/* Function: Translates a spatial velocity vector into a transformation matrix
 * Inputs: Spatial velocity vector [angular velocity, linear velocity]
 * Returns: Transformation matrix
 */
MatX VecToSE3(const VecX& V) {
    // Separate angular (exponential representation) and linear velocities
    Vec3 exp(V(0), V(1), V(2));
    Vec3 linear(V(3), V(4), V(5));

    // Fill in values to the appropriate parts of the transformation matrix
    MatX m_ret(4, 4);
    m_ret << VecToSO3(exp), linear, 0, 0, 0, 0;

    return m_ret;
}

/* Function: Translates a transformation matrix into a spatial velocity vector
 * Inputs: Transformation matrix
 * Returns: Spatial velocity vector [angular velocity, linear velocity]
 */
VecX SE3ToVec(const MatX& T) {
    VecX m_ret(6);
    m_ret << T(2, 1), T(0, 2), T(1, 0), T(0, 3), T(1, 3), T(2, 3);

    return m_ret;
}

/* Function: Provides the adjoint representation of a transformation matrix
 *			 Used to change the frame of reference for spatial velocity vectors
 * Inputs: 4x4 Transformation matrix SE(3)
 * Returns: 6x6 Adjoint Representation of the matrix
 */
MatX Adjoint(const MatX& T) {
    std::vector<MatX> R = TransToRp(T);
    MatX ad_ret(6, 6);
    ad_ret = MatX::Zero(6, 6);
    MatX zeroes = MatX::Zero(3, 3);
    ad_ret << R[0], zeroes, VecToSO3(R[1]) * R[0], R[0];
    return ad_ret;
}

/* Function: Rotation expanded for screw axis
 * Inputs: se3 matrix representation of exponential coordinates (transformation matrix)
 * Returns: 6x6 Matrix representing the rotation
 */
MatX MatrixExp6(const MatX& se3mat) {
    // Extract the angular velocity vector from the transformation matrix
    Mat3 se3mat_cut = se3mat.block<3, 3>(0, 0);
    Vec3 omgtheta = SO3ToVec(se3mat_cut);

    MatX m_ret(4, 4);

    if (NearZero(omgtheta.norm())) {
        // Reuse previous variables that have our required size
        se3mat_cut = MatX::Identity(3, 3);
        omgtheta << se3mat(0, 3), se3mat(1, 3), se3mat(2, 3);
        m_ret << se3mat_cut, omgtheta, 0, 0, 0, 1;
        return m_ret;
    }
    else {
        double theta = (AxisAng3(omgtheta))(3);
        Mat3 omgmat = se3mat.block<3, 3>(0, 0) / theta;
        Mat3 expExpand = MatX::Identity(3, 3) * theta + (1 - std::cos(theta)) * omgmat +
                         ((theta - std::sin(theta)) * (omgmat * omgmat));
        Vec3 linear(se3mat(0, 3), se3mat(1, 3), se3mat(2, 3));
        Vec3 GThetaV = (expExpand * linear) / theta;
        m_ret << MatrixExp3(se3mat_cut), GThetaV, 0, 0, 0, 1;
        return m_ret;
    }
}

MatX MatrixLog6(const MatX& T) {
    MatX m_ret(4, 4);
    auto rp = kincpp::TransToRp(T);
    Mat3 omgmat = MatrixLog3(rp.at(0));
    Mat3 zeros3d = Mat3::Zero(3, 3);
    if (NearZero(omgmat.norm())) {
        m_ret << zeros3d, rp.at(1), 0, 0, 0, 0;
    }
    else {
        double theta = std::acos((rp.at(0).trace() - 1) / 2.0);
        Mat3 logExpand1 = MatX::Identity(3, 3) - omgmat / 2.0;
        Mat3 logExpand2 = (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2) * omgmat * omgmat / theta;
        Mat3 logExpand = logExpand1 + logExpand2;
        m_ret << omgmat, logExpand * rp.at(1), 0, 0, 0, 0;
    }
    return m_ret;
}

MatX TransInv(const MatX& transform) {
    auto rp = kincpp::TransToRp(transform);
    auto Rt = rp.at(0).transpose();
    auto t = -(Rt * rp.at(1));
    MatX inv(4, 4);
    inv = MatX::Zero(4, 4);
    inv.block(0, 0, 3, 3) = Rt;
    inv.block(0, 3, 3, 1) = t;
    inv(3, 3) = 1;
    return inv;
}

}  // namespace kincpp
