#include "kincpp/utils.h"

namespace kincpp
{

/* Function: Check if a value is approximately zero
 * Inputs:
 *   - val: The value to check (double).
 * Returns:
 *   - true if |val| < 1e-6, false otherwise.
 */
bool NearZero(const double val) {
    return (std::abs(val) < 1e-6);
}

/* Function: Compute the skew-symmetric matrix of an angular velocity vector
 * Inputs:
 *   - omg: 3x1 angular velocity vector (Vec3).
 * Returns:
 *   - 3x3 skew-symmetric matrix (Mat3).
 */
Mat3 VecToSO3(const Vec3& omg) {
    return (Mat3() << 0, -omg(2), omg(1), omg(2), 0, -omg(0), -omg(1), omg(0), 0).finished();
}

/* Function: Convert skew-symmetric matrix to angular velocity vector
 * Inputs:
 *   - so3mat: 3x3 skew-symmetric matrix (Mat3).
 * Returns:
 *   - 3x1 angular velocity vector (Vec3).
 */
Vec3 SO3ToVec(const Mat3& so3mat) {
    return Vec3(so3mat(2, 1), so3mat(0, 2), so3mat(1, 0));
}

/* Function: Compute the rotation matrix from a skew-symmetric matrix
 * Inputs:
 *   - so3mat: Skew-symmetric matrix (Mat3).
 * Returns:
 *   - 3x3 rotation matrix (Mat3).
 */
Mat3 MatrixExp3(const Mat3& so3mat) {
    const Vec3 omgtheta = SO3ToVec(so3mat);
    const double theta = omgtheta.norm();

    if (NearZero(theta)) {
        return Mat3::Identity();
    }

    const Mat3 omgmat = so3mat / theta;
    return Mat3::Identity() + std::sin(theta) * omgmat + (1 - std::cos(theta)) * (omgmat * omgmat);
}

/* Function: Compute the matrix logarithm of a rotation matrix
 * Inputs:
 *   - R: 3x3 rotation matrix (Mat3).
 * Returns:
 *   - 3x3 skew-symmetric matrix (Mat3).
 */
Mat3 MatrixLog3(const Mat3& R) {
    const double acosinput = (R.trace() - 1) / 2.0;

    if (acosinput >= 1) {
        return Mat3::Zero();  // No rotation
    }
    if (acosinput <= -1) {
        // Handle 180-degree rotation
        Vec3 omg;
        if (!NearZero(1 + R(2, 2))) {
            omg = (1.0 / std::sqrt(2 * (1 + R(2, 2)))) * Vec3(R(0, 2), R(1, 2), 1 + R(2, 2));
        }
        else if (!NearZero(1 + R(1, 1))) {
            omg = (1.0 / std::sqrt(2 * (1 + R(1, 1)))) * Vec3(R(0, 1), 1 + R(1, 1), R(2, 1));
        }
        else {
            omg = (1.0 / std::sqrt(2 * (1 + R(0, 0)))) * Vec3(1 + R(0, 0), R(1, 0), R(2, 0));
        }
        return VecToSO3(PI * omg);
    }

    // General case
    const double theta = std::acos(acosinput);
    return theta / (2.0 * std::sin(theta)) * (R - R.transpose());
}

/* Function: Extract rotation matrix and position vector from a transformation matrix
 * Inputs:
 *   - T: 4x4 transformation matrix (Mat4).
 * Returns:
 *   - Vector containing [rotation matrix, position vector].
 */
std::vector<MatX> TransToRp(const Mat4& T) {
    return {T.block<3, 3>(0, 0), T.block<3, 1>(0, 3)};
}

/* Function: Construct a transformation matrix from a spatial velocity vector
 * Inputs:
 *   - V: Spatial velocity vector (Vec6).
 * Returns:
 *   - 4x4 transformation matrix (Mat4).
 */
Mat4 VecToSE3(const Vec6& V) {
    Mat4 result = Mat4::Zero();
    result.block<3, 3>(0, 0) = VecToSO3(V.head<3>());
    result.block<3, 1>(0, 3) = V.tail<3>();
    return result;
}

/* Function: Convert transformation matrix to spatial velocity vector
 * Inputs:
 *   - T: 4x4 transformation matrix (Mat4).
 * Returns:
 *   - Spatial velocity vector (Vec6).
 */
Vec6 SE3ToVec(const Mat4& T) {
    Vec6 result;
    result.head<3>() = SO3ToVec(T.block<3, 3>(0, 0));
    result.tail<3>() = T.block<3, 1>(0, 3);
    return result;
}

/* Function: Compute the adjoint representation of a transformation matrix
 * Inputs:
 *   - T: 4x4 transformation matrix (Mat4).
 * Returns:
 *   - 6x6 adjoint matrix (Mat6).
 */
Mat6 Adjoint(const Mat4& T) {
    const auto rp = TransToRp(T);
    Mat6 result = Mat6::Zero();
    result.block<3, 3>(0, 0) = rp[0];
    result.block<3, 3>(3, 0) = VecToSO3(rp[1]) * rp[0];
    result.block<3, 3>(3, 3) = rp[0];
    return result;
}

/* Function: Compute the exponential map for a screw axis (SE(3))
 * Inputs:
 *   - se3mat: 4x4 matrix representation of exponential coordinates (Mat4).
 * Returns:
 *   - 4x4 transformation matrix (Mat4).
 */
Mat4 MatrixExp6(const Mat4& se3mat) {
    // Extract the angular velocity (rotation part) from the transformation matrix
    Mat3 se3mat_cut = se3mat.block<3, 3>(0, 0);
    Vec3 omgtheta = SO3ToVec(se3mat_cut);
    double theta = omgtheta.norm();

    // Initialize the return matrix
    Mat4 m_ret = Mat4::Identity();

    // Case 1: No rotation (theta near zero)
    if (NearZero(theta)) {
        m_ret.block<3, 1>(0, 3) = se3mat.block<3, 1>(0, 3);  // Only translation
        return m_ret;
    }

    // Case 2: Non-zero rotation
    Mat3 omgmat = se3mat_cut / theta;        // Normalize the skew-symmetric matrix
    Mat3 rotation = MatrixExp3(se3mat_cut);  // Compute the rotational part

    // Compute the translation part using the closed-form expansion
    Mat3 expExpand = Mat3::Identity() * theta + (1 - std::cos(theta)) * omgmat +
                     ((theta - std::sin(theta)) * (omgmat * omgmat));
    Vec3 linear = se3mat.block<3, 1>(0, 3);
    Vec3 GThetaV = expExpand * linear / theta;

    // Combine rotation and translation
    m_ret.block<3, 3>(0, 0) = rotation;
    m_ret.block<3, 1>(0, 3) = GThetaV;

    return m_ret;
}

/* Function: Compute the matrix logarithm for a transformation matrix (SE(3))
 * Inputs:
 *   - T: 4x4 transformation matrix (Mat4).
 * Returns:
 *   - 4x4 matrix logarithm (Mat4).
 */
Mat4 MatrixLog6(const Mat4& T) {
    // Decompose the transformation matrix into rotation and position
    auto rp = kincpp::TransToRp(T);
    Mat3 R = rp[0];
    Vec3 p = rp[1];

    // Compute the logarithm of the rotation matrix
    Mat3 omgmat = MatrixLog3(R);
    double theta = omgmat.norm();

    // Initialize the return matrix
    Mat4 m_ret = Mat4::Zero();

    // Case 1: No rotation (theta near zero)
    if (NearZero(theta)) {
        m_ret.block<3, 1>(0, 3) = p;  // Only translation
        return m_ret;
    }

    // Case 2: Non-zero rotation
    Mat3 logExpand1 = Mat3::Identity() - omgmat / 2.0;
    Mat3 logExpand2 =
        ((1.0 / theta) - (1.0 / std::tan(theta / 2.0)) / 2.0) * (omgmat * omgmat) / theta;
    Mat3 logExpand = logExpand1 + logExpand2;

    // Combine rotation and translation components
    m_ret.block<3, 3>(0, 0) = omgmat;
    m_ret.block<3, 1>(0, 3) = logExpand * p;

    return m_ret;
}

/* Function: Compute the inverse of a transformation matrix
 * Inputs:
 *   - transform: 4x4 transformation matrix (Mat4).
 * Returns:
 *   - Inverse transformation matrix (Mat4).
 */
Mat4 TransInv(const Mat4& transform) {
    const auto rp = TransToRp(transform);
    Mat4 result = Mat4::Identity();
    result.block<3, 3>(0, 0) = rp[0].transpose();
    result.block<3, 1>(0, 3) = -rp[0].transpose() * rp[1];
    return result;
}

Vec6 AngleAxisDiff(const Mat4& T, const Mat4& Td) {
    Vec6 e;

    // Extract rotation matrices and position vectors
    const auto T_rp = TransToRp(T);
    const auto Td_rp = TransToRp(Td);

    const Mat3 R = Td_rp.at(0) * T_rp.at(0).transpose();
    const Vec3 p = Td_rp.at(1) - T_rp.at(1);

    // Set translational error
    e.head<3>() = p;

    // Compute rotational error
    const Vec3 li = SO3ToVec(R - R.transpose());
    if (NearZero(li.norm())) {
        // Handle diagonal case
        if (R.trace() > 0) {
            e.tail<3>().setZero();  // No rotation
        }
        else {
            e.tail<3>() = (PI / 2.0) * (R.diagonal().array() + 1.0).matrix();
        }
    }
    else {
        // Non-diagonal case
        const double ln = li.norm();
        e.tail<3>() = std::atan2(ln, R.trace() - 1.0) * li / ln;
    }

    return e;
}

/* Function: Add bounds as linear constraints to an inequality system
 * Inputs:
 *   - Ain: Optional input constraint matrix (MatX).
 *   - bin: Optional input constraint vector (VecX).
 *   - lb: Lower bounds (VecX).
 *   - ub: Upper bounds (VecX).
 * Returns:
 *   - Pair of updated constraint matrix (MatX) and vector (VecX).
 */
std::pair<MatX, VecX> AddBoundConstraints(const std::optional<MatX>& Ain,
                                          const std::optional<VecX>& bin, const VecX& lb,
                                          const VecX& ub) {
    // Initialize A and b with optional inputs or default empty values
    MatX A = Ain.value_or(MatX(0, 0));
    VecX b = bin.value_or(VecX(0));

    // Helper function to add constraints (lower or upper)
    auto add_constraints = [](MatX& A, VecX& b, const MatX& bound_A, const VecX& bound_b) {
        if (A.rows() == 0) {
            // Initialize with the bound constraints
            A = bound_A;
            b = bound_b;
        }
        else {
            // Append the bound constraints
            A.conservativeResize(A.rows() + bound_A.rows(), Eigen::NoChange);
            b.conservativeResize(b.size() + bound_b.size());
            A.bottomRows(bound_A.rows()) = bound_A;
            b.tail(bound_b.size()) = bound_b;
        }
    };

    // Add lower bound constraints
    if (lb.size() > 0) {
        MatX lower_bound_A = -MatX::Identity(lb.size(), lb.size());
        VecX lower_bound_b = -lb;
        add_constraints(A, b, lower_bound_A, lower_bound_b);
    }

    // Add upper bound constraints
    if (ub.size() > 0) {
        MatX upper_bound_A = MatX::Identity(ub.size(), ub.size());
        const VecX& upper_bound_b = ub;
        add_constraints(A, b, upper_bound_A, upper_bound_b);
    }

    // Return the updated constraint matrix and vector
    return {A, b};
}
}  // namespace kincpp
