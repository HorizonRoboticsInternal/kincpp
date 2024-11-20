from kincpp import angle_axis_cpp
import numpy as np
import roboticstoolbox as rtb
# from roboticstoolbox.tools import angle_axis_python


# Define M as a NumPy array
M = np.array([
    [1.0, 0.0, 0.0, 0.0872],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.1605],
    [0.0, 0.0, 0.0, 1.0]
])

# Define Slist as a NumPy array
Slist = np.array([
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    [0.0, 1.0, 1.0, 1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    [-0.0, -0.1035, -0.1035, -0.1605, -0.0, -0.0],
    [-0.0, -0.0, 0.0, 0.0, 0.062, 0.1605],
    [-0.0, -0.0, -0.35, -0.132, -0.0, -0.0]
])

# Define lower_lim as a NumPy array
lower_lim = np.array([
    -2.6179938779914944, 0.0, -2.885592653589793, 
    -1.5184364492350666, -1.3439035240356338, -2.792526803190927
])

# Define upper_lim as a NumPy array
upper_lim = np.array([
    2.6179938779914944, 2.8797932657906435, 0.0, 
    1.5184364492350666, 1.3439035240356338, 2.792526803190927
])

tf_matrix1 = np.array([
    [1, 0, 0, 3],  # Translate by 3 along X
    [0, 1, 0, 4],  # Translate by 4 along Y
    [0, 0, 1, 0.5],  # Translate by 5 along Z
    [0, 0, 0, 1]   # Homogeneous coordinate
])

theta = np.pi / 4  # 45 degrees
cos_theta, sin_theta = np.cos(theta), np.sin(theta)

tf_matrix2 = np.array([
    [cos_theta, -sin_theta, 0, 0],  # Rotation around Z-axis
    [sin_theta, cos_theta,  0, 0],
    [0,         0,          1, 0],
    [0,         0,          0, 1]
])

theta = np.pi / 2  # 90 degrees
cos_theta, sin_theta = np.cos(theta), np.sin(theta)

tf_matrix3 = np.array([
    [cos_theta,  0, sin_theta, 1],  # Rotation + Translation
    [0,          1, 0,         2],
    [-sin_theta, 0, cos_theta, 3],
    [0,          0, 0,         1]
])

# print(angle_axis_cpp(tf_matrix1, tf_matrix2))
# print(angle_axis_python(tf_matrix1, tf_matrix2))
# print('====')
# print(angle_axis_cpp(tf_matrix2, tf_matrix3))
# print(angle_axis_python(tf_matrix2, tf_matrix3))
# print('====')
# print(angle_axis_cpp(tf_matrix3, tf_matrix1))
# print(angle_axis_python(tf_matrix3, tf_matrix1))
# print('====')
# print(angle_axis_cpp(tf_matrix1, tf_matrix1))
# print(angle_axis_python(tf_matrix1, tf_matrix1))

z1 = rtb.models.URDF.Z1()
solver = rtb.IK_QP(slimit=1, ilimit=5)

q0 = np.zeros(6)
delta_limits = np.ones((2, 6))
delta_limits[0] *= -0.05
delta_limits[1] *= 0.05

try:
    print(solver.solve(z1.ets(), tf_matrix1.copy(), q0, delta_limits).q)
except RuntimeError:
    pass


import kincpp

solver = kincpp.KinematicsSolver(M, Slist, lower_lim, upper_lim)

ik_param = kincpp.IKParams()
ik_param.common_params.desired_ee_tf = tf_matrix1
ik_param.common_params.q_guess = q0
ik_param.common_params.iters = 5
ik_param.solver_type = kincpp.QP

print(solver.inverse_kinematics(ik_param))
# solver.solve(tf_matrix1, tf_matrix2, tf_matrix3)
