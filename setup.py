from setuptools import setup


module_name = "kincpp"

setup(
    name=module_name,
    version='0.1',
    author='Andrew Choi',
    author_email='andrew.choi@horizon.cc',
    description='A pybind C++ kinematics module for manipulators.',
    packages=[module_name],
    package_dir={module_name: module_name},
    package_data={module_name: ['*.so']},
    python_requires='==3.10',
)
