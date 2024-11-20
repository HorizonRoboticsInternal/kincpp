cmake_minimum_required(VERSION 3.0)
project(kincpp)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

set(LIBRARY_OUTPUT_PATH "../kincpp")

pybind11_add_module(kincpp
        src/kinematics_solver.cpp
        src/utils.cpp)

target_link_libraries(kincpp PRIVATE
    Eigen3::Eigen
    quadprog
)
