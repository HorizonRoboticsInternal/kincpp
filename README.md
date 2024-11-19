# kincpp
kincpp is a C++ library for computing both forward and inverse kinematics for arbitrary robot
manipulators. It is written entirely in C++ for speed and has pybind bindings for python.
Any manipulator can be modelled by simply providing its home configuration end-effector transform `M` and the joint screw axes `S`.
A tool to compute this for modelled robots can be found [here](https://github.com/Interbotix/kinematics_from_description).

Core logic to compute FK and Newton-based IK was largely extracted / refactored from
[ModernRoboticsCpp](https://github.com/Le0nX/ModernRoboticsCpp)
(a C++ implementation of [ModernRobotics](https://github.com/NxRLab/ModernRobotics)).
Constraint-based IK using a quadratic programming approach is also available and was ported over
to C++ from [robotics-toolbox-python](https://github.com/petercorke/robotics-toolbox-python).
Tentative profiling results show an order of magnitude speed increase over the original Python
version.


### Build from Source
To build from source, first install the necessary dependencies `Eigen` and `pybind`.
```bash
sudo apt install libeigen3-dev pybind11-dev
```
Then, simply build the shared library file
```bash
mkdir build && cd build && cmake ..
make -j4
```
This will produce a `kincpp.cpython-{YOUR_PYTHON_VERSION}-{YOUR_ARCH}.so` file that you can use directly for your python applications. You can also install it as python package by running
```bash
pip install .
```
A static library is also created which can be linked to in cmake using the target `kincpp_lib`.

### Testing
Once `kincpp` is installed in some way, you can test that everything is working by running the provided unit test.
```bash
python testing/test_kinematics.py
python testing/ik_qp_test.py
```
These files also contain useful examples for instantiating a `KinematicsSolver` and calling
the forward and inverse kinematics functions.

### Troubleshooting
If your system cannot find `Eigen`, you may need to create symlinks if your `Eigen` package is actually named `eigen3/Eigen`.
Simply locate your `Eigen` installation location, `locate eigen`, and then create the necessary symlinks.
For example, if `Eigen` is found in `/usr/include`, then you would do the following:
```bash
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
```
