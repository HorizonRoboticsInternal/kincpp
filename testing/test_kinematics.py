"""
    Forward and inverse kinematics unit tests.
"""
import numpy as np
import kincpp


def fk_test():
    """
    Tests forward kinematics function.
    """
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

    joint_positions = np.array([np.pi / 2.0, 3, np.pi])

    expected_res = np.array([[0, 1, 0, -5],
                             [1, 0, 0, 4],
                             [0, 0, -1, 1.68584073],
                             [0, 0, 0, 1]])

    fk_res = kincpp.forward(M, S, joint_positions)

    np.testing.assert_allclose(fk_res, expected_res, atol=1e-4)

    print("Forward Kinematics Test Passed.")


def ik_test():
    """
    Tests inverse kinematics function.
    """
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

    T = np.array([[0, 1, 0, -5],
                  [1, 0, 0, 4],
                  [0, 0, -1, 1.6858],
                  [0, 0, 0, 1]])

    joint_position_guess = np.array([1.5, 2.5, 3])

    position_tolerance = 1e-3
    orientation_tolerance = 1e-3

    expected_res = np.array([1.57073783, 2.99966384, 3.1415342], dtype=np.float64)

    success, ik_res = kincpp.inverse(M, S, T, joint_position_guess, position_tolerance, orientation_tolerance)

    assert success == True
    assert np.allclose(ik_res, expected_res, atol=1e-4)

    print("Inverse Kinematics Test Passed.")


if __name__ == '__main__':
    fk_test()
    ik_test()
