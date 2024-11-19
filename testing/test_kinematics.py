"""
    Forward and inverse kinematics unit tests.
"""
import numpy as np
import kincpp

M = np.array([[-1, 0, 0, 0],
              [0, 1, 0, 6],
              [0, 0, -1, 2],
              [0, 0, 0, 1]])

S = np.array([[0, 0, 0],
              [0, 0, 0],
              [1, 0, -1],
              [4, 0, -6],
              [0, 1, 0],
              [0, 0, -0.1]])


def fk_test():
    """
    Tests forward kinematics function.
    """
    solver = kincpp.KinematicsSolver(M, S,
                                     lower_joint_limits=-np.ones(3) * np.inf,
                                     upper_joint_limits=np.ones(3) * np.inf)

    q_guess = np.array([np.pi / 2.0, 3, np.pi])

    expected_res = np.array([[0, 1, 0, -5],
                             [1, 0, 0, 4],
                             [0, 0, -1, 1.68584073],
                             [0, 0, 0, 1]])

    fk_res = solver.forward_kinematics(q_guess)

    np.testing.assert_allclose(fk_res, expected_res, atol=1e-4)

    print("Forward Kinematics Test Passed.")


def ik_test():
    """
    Tests inverse kinematics function.
    """
    desired_ee_tf = np.array([[0, 1, 0, -5],
                              [1, 0, 0, 4],
                              [0, 0, -1, 1.6858],
                              [0, 0, 0, 1]])

    solver = kincpp.KinematicsSolver(M, S,
                                     lower_joint_limits=-np.ones(3) * np.inf,
                                     upper_joint_limits=np.ones(3) * np.inf)

    q_guess = np.array([1.5, 2.5, 3])

    position_tolerance = 1e-3
    orientation_tolerance = 1e-3

    expected_res = np.array([1.57073783, 2.99966384, 3.1415342], dtype=np.float64)

    ik_params = kincpp.IKParams()
    ik_params.common_params.desired_ee_tf = desired_ee_tf
    ik_params.common_params.q_guess = q_guess
    ik_params.common_params.iters = 20
    ik_params.newton_params.project_to_joint_limits = False
    ik_params.newton_params.use_pseudo_inverse = False
    ik_params.solver_type = kincpp.NEWTON

    success, ik_res = solver.inverse_kinematics(ik_params)

    assert success == True
    assert np.allclose(ik_res, expected_res, atol=1e-4)

    print("Inverse Kinematics Test Passed.")


def ja_test():
    solver = kincpp.KinematicsSolver(M, S,
                                     lower_joint_limits=-np.ones(3) * np.inf,
                                     upper_joint_limits=np.ones(3) * np.inf)

    curr_q = np.array([1.5, 2.5, 3])

    vt_jac = solver.velocity_twist_jacobian(curr_q)
    sv_jac = solver.spatial_velocity_jacobian(curr_q)

    ee_pos = solver.forward_kinematics(curr_q)[:3, -1]
    # skew-symmetric matrix
    ee_pos_mat = np.array([
        [0, -ee_pos[2], ee_pos[1]],
        [ee_pos[2], 0, -ee_pos[0]],
        [-ee_pos[1], ee_pos[0], 0]
    ])

    vel_tf = np.eye(6)
    vel_tf[:3, 3:] = -ee_pos_mat

    assert np.allclose(sv_jac, vel_tf @ vt_jac, atol=1e-6)
    print("Jacobian Test Passed.")


if __name__ == '__main__':
    fk_test()
    ik_test()
    ja_test()
