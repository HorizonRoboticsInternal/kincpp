#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <optional>
#include <QuadProg++/Array.hh>

namespace kincpp
{
constexpr double PI = 3.141592653589793;

using ArrX = Eigen::ArrayXd;
using Vec3 = Eigen::Vector3d;
using Vec6 = Eigen::Vector<double, 6>;
using VecX = Eigen::VectorXd;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;
using MatX = Eigen::MatrixXd;

/* Function: Find if the value is negligible enough to consider 0 */
bool NearZero(const double val);

/* Function: Returns the skew symmetric matrix representation of an angular velocity vector */
Mat3 VecToSO3(const Vec3& omg);

/* Function: Returns a normalized version of the input vector */
MatX Normalize(MatX V);

/* Function: Returns angular velocity vector represented by the skew symmetric matrix */
Vec3 SO3ToVec(const MatX& so3mat);

/* Function: Translates an exponential rotation into its individual components */
Eigen::Vector4d AxisAng3(const Vec3& expc3);

/* Function: Translates an exponential rotation into a rotation matrix */
Mat3 MatrixExp3(const Mat3& so3mat);

/* Function: Computes the matrix logarithm of a rotation matrix */
Mat3 MatrixLog3(const Mat3& R);

/* Function: Separates the rotation matrix and position vector from the transformation matrix
 * representation */
std::vector<MatX> TransToRp(const MatX& T);

/* Function: Translates a spatial velocity vector into a transformation matrix */
MatX VecToSE3(const VecX& V);

/* Function: Translates a transformation matrix into a spatial velocity vector */
VecX SE3ToVec(const MatX& T);

/* Function: Provides the adjoint representation of a transformation matrix */
MatX Adjoint(const MatX& T);

/* Function: Rotation expanded for screw axis */
MatX MatrixExp6(const MatX& se3mat);

/* Function: Computes the matrix logarithm of a transformation matrix */
MatX MatrixLog6(const MatX& T);

/* Function: Computes the inverse of a homogeneous transformation matrix */
MatX TransInv(const MatX& transform);

VecX AngleAxisCpp(const MatX& T, const MatX& Td);

std::pair<MatX, VecX> AddBoundConstraints(const std::optional<MatX>& Ain,
                                          const std::optional<VecX>& bin,
                                          const VecX& lb,
                                          const VecX& ub);

quadprogpp::Vector<double> ConvertEigenVecToQPVec(const VecX& eigen_vec);

quadprogpp::Matrix<double> ConvertEigenMatToQPMat(const MatX& eigen_mat);

}  // namespace kincpp
