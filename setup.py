from setuptools import setup, find_packages

setup(
    name="kincpp",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'kincpp': ['*.so', '*.py'],
    },
    description="A pybind C++ kinematics module for manipulators.",
    author="Andrew Choi",
    author_email="neffneff4@gmail.com",
)
