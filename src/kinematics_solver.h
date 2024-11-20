#pragma once

#include "utils.h"

namespace kincpp
{

typedef enum {
    NEWTON = 0,
    QP = 1
} SolverType;

struct CommonParams
{
    Mat4 desired_ee_tf;
    VecX q_guess;
    int iters{};
};

struct NewtonParams
{
    double position_tolerance = 1e-3;
    double orientation_tolerance = 1e-3;
    bool project_to_joint_limits = false;
    bool use_pseudo_inverse = false;
};

struct QPParams
{};

struct IKParams
{
    CommonParams common_params;
    NewtonParams newton_params;
    QPParams qp_params;
    SolverType solver_type = NEWTON;
};

class KinematicsSolver
{
  public:
    KinematicsSolver(const MatX& M, const MatX& S, const VecX& lower_joint_limits,
                     const VecX& upper_joint_limits);
    ~KinematicsSolver();

    MatX ForwardKinematics(const VecX& q);
    std::pair<bool, VecX> InverseKinematics(const IKParams& params);

    MatX VelocityTwistJacobian(const VecX& q);
    MatX SpatialVelocityJacobian(const VecX& q);

  protected:
    std::pair<bool, VecX> IK_NM(const CommonParams& common_params,
                                const NewtonParams& newton_params);
    std::pair<bool, VecX> IK_QP(const CommonParams& common_params, const QPParams& qp_params);

    MatX M;
    MatX S;
    ArrX lower_joint_limits;
    ArrX upper_joint_limits;
};

}  // namespace kincpp
