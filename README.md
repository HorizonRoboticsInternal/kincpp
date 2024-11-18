# kincpp
This repo contains the tools to create the `kincpp` pip package, which contains `pybind` wrapper functions for forward and inverse kinematics written in C++.

The core C++ logic is almost entirely from [ModernRoboticsCpp](https://github.com/Le0nX/ModernRoboticsCpp) (a C++ implementation of [ModernRobotics](https://github.com/NxRLab/ModernRobotics)) with slight modifications.

### Pip installation

A pip package currently exists but binary wheels exist solely for Python3.10 on Linux. If compatible, you may install from pip as follows
```bash
python3.10 -m pip install kincpp
```

### Build from Source
To build from source, first install the necessary dependencies `Eigen` and `pybind`.
```bash
sudo apt install libeigen3-dev
pip install pybind11
```
Then, simply build the shared library file
```bash
mkdir build && cd build && cmake ..
make -j4
```
This will produce a `kincpp.cpython-{YOUR_PYTHON_VERSION}-{YOUR_ARCH}.so` file that you can use directly for your python applications. You can also install it as python package by running
```bash
python setup.py bdist_wheel && cd dist
pip install kincpp.0.1-py3-none-any.whl
```

### Update pip package

To rebuild the python wheel and update the package, run
```bash
python setup.py bdist_wheel
twine upload dist/{NAME_OF_WHEEL}
```

### Testing
Once `kincpp` is installed in some way, you can test that everything is working by running the provided unit test.
```bash
python testing/test_kinematics.py
```

### Troubleshooting
If your system cannot find `Eigen`, you may need to create symlinks if your `Eigen` package is actually named `eigen3/Eigen`.
Simply locate your `Eigen` installation location, `locate eigen`, and then create the necessary symlinks.
For example, if `Eigen` is found in `/usr/include`, then you would do the following:
```bash
cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
```
