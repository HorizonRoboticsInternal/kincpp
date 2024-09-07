import numpy as np
from loguru import logger
from typing import List, Optional
from scipy.spatial.transform import Rotation as R
from abc import ABC, abstractmethod
import lib.kincpp as kincpp


class KinematicsSolver(ABC):
    """
    Abstract class for solving manipulator forward and inverse kinematics.

    Users can implement a solver for a particular arm by inheriting this class
    and specifying the ee_home_config and joint_screw_axes params.
    """

    @property
    @abstractmethod
    def ee_home_config(self) -> np.ndarray:
        """
        Getter for the end effector home configuration as a 4x4 TF matrix.

        Returns:
            The end effector home configuration as a 4x4 TF matrix.

        """
        pass

    @property
    @abstractmethod
    def joint_screw_axes(self) -> np.ndarray:
        """
        Getter for the joint screw axes.

        Returns:
            The joint screw axes as a 6xN matrix, where N is the number of joints.

        """
        pass

    @property
    @abstractmethod
    def lower_joint_limits(self) -> np.ndarray:
        """
        Getter for the joint lower limits.

        Returns:
            The joint lower limits as an N size vector, where N is the number of joints.

        """
        pass

    @property
    @abstractmethod
    def upper_joint_limits(self) -> np.ndarray:
        """
        Getter for the joint upper limits.

        Returns:
            The joint upper limits as an N size vector, where N is the number of joints.

        """
        pass

    @property
    def rev(self) -> float:
        return 2 * np.pi

    def fwd_kin(self, joint_positions: np.ndarray) -> np.ndarray:
        """
        Kincpp forward kinematics wrapper function.

        Args:
            joint_positions: The current joint angles.

        Returns:
            The current end effector TF.
        """
        return _kincpp.forward(self.ee_home_config, self.joint_screw_axes,
                               joint_positions)

    def inv_kin(self,
                desired_ee_tf: np.ndarray,
                joint_position_guess: np.ndarray,
                position_tolerance: float = 1e-3,
                orientation_tolerance: float = 1e-3,
                max_iterations: int = 20) -> (bool, bool, np.ndarray):
        """
        Kincpp inverse kinematics wrapper function.

        Args:
            desired_ee_tf: The desired end effector TF.
            joint_position_guess: The joint position initial guess for the IK solver.
                If the ee delta is small enough, a safe guess is simply the current
                joint position.
            position_tolerance: The end effector Cartesian position tolerance.
            orientation_tolerance: The end effector orientation tolerance.
            max_iterations: The number of iterations before IK solver quits.

        Returns:
            A tuple containing three things:
                1) A bool for whether IK succeeded
                    Note if IK failed, the joint angle results are undefined.
                2) A bool for whether any joint limits were triggered
                    If this is true, then the output position will not
                    match the desired_ee_tf
                3) The joint angles.
        """
        success, joint_positions = _kincpp.inverse(
            self.ee_home_config, self.joint_screw_axes, desired_ee_tf,
            joint_position_guess, position_tolerance, orientation_tolerance,
            max_iterations)

        joint_limit_violation, joint_positions = self._wrap_joint_positions(
            joint_positions)

        return success, joint_limit_violation, joint_positions

    def create_cartesian_traj(
        self,
        initial_joint_position: np.ndarray,
        desired_xyz: Optional[np.ndarray] = None,
        desired_rpy: Optional[np.ndarray] = None,
        resolution: int = 40,
    ) -> List[np.ndarray]:
        """
        Generates a Cartesian trajectory based on the desired end-effector position and orientation,
        the initial joint positions, and the desired resolution.

        Args:
            initial_joint_position: The initial joint positions of the arm. This is an
                initial seed point for the Newton solver, so changing can result in different
                converged solutions.
            desired_xyz: The desired end-effector position in Cartesian coordinates.
                If None, the current end-effector position is used.
            desired_rpy: The desired end-effector orientation as Roll-Pitch-Yaw angles.
                If None, the current end-effector orientation is used.
            resolution: The number of steps to interpolate the trajectory.

        Returns:
            A list of joint positions corresponding to the interpolated Cartesian trajectory.
        """
        starting_ee_tf = self.fwd_kin(initial_joint_position)
        curr_xyz = starting_ee_tf[:3, -1]

        if desired_xyz is None:
            delta_xyz = np.zeros(3, dtype=np.float32)
        else:
            delta_xyz = (desired_xyz - curr_xyz) / resolution

        if desired_rpy is None:
            delta_rotation_matrix = np.eye(3, dtype=np.float32)
        else:
            # Calculate the desired rotation matrix using scipy
            desired_rotation_matrix = R.from_euler('xyz',
                                                   desired_rpy).as_matrix()

            # Calculate the current rotation matrix from the initial transformation
            current_rotation_matrix = starting_ee_tf[:3, :3]

            # Calculate the incremental rotation matrix
            delta_rotation_matrix = R.from_matrix(
                desired_rotation_matrix @ np.linalg.inv(
                    current_rotation_matrix)).as_rotvec() / resolution
            delta_rotation_matrix = R.from_rotvec(
                delta_rotation_matrix).as_matrix()

        joint_trajectory: List[np.ndarray] = []
        curr_ee_tf = starting_ee_tf
        guess = initial_joint_position
        for _ in range(resolution):
            curr_ee_tf[:3, -1] += delta_xyz
            curr_ee_tf[:3, :3] @= delta_rotation_matrix

            success, joint_limit_violation, joint_positions = self.inv_kin(
                curr_ee_tf, guess)

            if not success:
                raise RuntimeError(
                    "Inverse kinematics has failed to find a solution.")
            if joint_limit_violation:
                logger.warning("Joint limit violation detected.")

            # Deployment requires the gripper position as well, so we append a zero
            joint_trajectory.append(joint_positions)

            # Update the guess for the next target.
            guess = joint_positions

        return joint_trajectory

    def _wrap_joint_positions(
            self, joint_positions: np.ndarray) -> (bool, np.ndarray):
        """
        Wrap an array of joint commands to [-pi, pi) and between the joint limits.

        Args:
            joint_positions: joint positions to wrap

        Returns:
            A tuple containing a bool indicating whether the joint limits were met
            and an array of joint positions wrapped between [-pi, pi).
        """
        joint_positions = (joint_positions + np.pi) % self.rev - np.pi

        under_limit = joint_positions < self.lower_joint_limits
        over_limit = joint_positions > self.upper_joint_limits

        # Check if any joints were under or over the limits
        joint_limit_violation = np.any(under_limit) or np.any(over_limit)

        # Apply corrections and set to limits
        joint_positions[under_limit] = self.lower_joint_limits[under_limit]
        joint_positions[over_limit] = self.upper_joint_limits[over_limit]

        return joint_limit_violation, joint_positions
