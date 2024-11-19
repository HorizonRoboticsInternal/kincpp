#pragma once

#include "utils.h"

namespace kincpp
{

typedef enum {
    NEWTON = 0,
    QP = 1
} SolverType;

class KinematicsSolver
{
  public:
    KinematicsSolver(const MatX& M, const MatX& S, const VecX& lower_joint_limits,
                     const VecX& upper_joint_limits, const SolverType& solver_type = QP);
    ~KinematicsSolver();

    MatX ForwardKinematics(const VecX& q);
    std::pair<bool, VecX>
    InverseKinematics(const MatX& desired_ee_tf, VecX& q_guess, double position_tolerance = 1e-3,
                      double orientation_tolerance = 1e-3, bool project_to_joint_limits = true,
                      bool use_pseudo_inverse = false, int max_iterations = 20);

    MatX VelocityTwistJacobian(const VecX& q);
    MatX SpatialVelocityJacobian(const VecX& q);

  protected:
    std::pair<bool, VecX> IK_NM(const MatX& desired_ee_tf, VecX& q_guess, double position_tolerance,
                                double orientation_tolerance, bool project_to_joint_limits = true,
                                bool use_pseudo_inverse = false, int max_iterations = 20);
    std::pair<bool, VecX> IK_QP(const MatX& desired_ee_tf, VecX& q_guess, double position_tolerance,
                                double orientation_tolerance, int max_iterations = 20);

    SolverType solver_type;
    MatX M;
    MatX S;
    ArrX lower_joint_limits;
    ArrX upper_joint_limits;
};

}  // namespace kincpp
