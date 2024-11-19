#pragma once

#include "global_definitions.h"

namespace kincpp
{

/* Function: Find if the value is negligible enough to consider 0 */
bool NearZero(const double val);

/* Function: Returns the skew symmetric matrix representation of an angular velocity vector */
Mat3 VecToSO3(const Vec3& omg);

/* Function: Returns angular velocity vector represented by the skew symmetric matrix */
Vec3 SO3ToVec(const Mat3& so3mat);

/* Function: Translates an exponential rotation into a rotation matrix */
Mat3 MatrixExp3(const Mat3& so3mat);

/* Function: Computes the matrix logarithm of a rotation matrix */
Mat3 MatrixLog3(const Mat3& R);

/* Function: Separates the rotation matrix and position vector from the transformation matrix
 * representation */
std::vector<MatX> TransToRp(const Mat4& T);

/* Function: Translates a spatial velocity vector into a transformation matrix */
Mat4 VecToSE3(const Vec6& V);

/* Function: Translates a transformation matrix into a spatial velocity vector */
Vec6 SE3ToVec(const Mat4& T);

/* Function: Provides the adjoint representation of a transformation matrix */
Mat6 Adjoint(const Mat4& T);

/* Function: Computes the inverse of a homogeneous transformation matrix */
Mat4 TransInv(const Mat4& transform);

/* Function: Rotation expanded for screw axis */
Mat4 MatrixExp6(const Mat4& se3mat);

/* Function: Computes the matrix logarithm of a transformation matrix */
Mat4 MatrixLog6(const Mat4& T);

/* Function: Compute the difference as an angle axis between two transformations */
Vec6 AngleAxisDiff(const Mat4& T, const Mat4& Td);

/* Function: Add bound constraints to an inequality matrix */
std::pair<MatX, VecX> AddBoundConstraints(const std::optional<MatX>& Ain,
                                          const std::optional<VecX>& bin, const VecX& lb,
                                          const VecX& ub);
}  // namespace kincpp
