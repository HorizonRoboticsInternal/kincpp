#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <optional>
#include <vector>

namespace kincpp
{
constexpr double PI = 3.141592653589793;
constexpr double INF = std::numeric_limits<double>::infinity();

using ArrX = Eigen::ArrayXd;
using Vec3 = Eigen::Vector3d;
using Vec6 = Eigen::Vector<double, 6>;
using VecX = Eigen::VectorXd;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;
using Mat6 = Eigen::Matrix<double, 6, 6>;
using MatX = Eigen::MatrixXd;

}  // namespace kincpp
