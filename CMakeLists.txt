cmake_minimum_required(VERSION 3.10)
project(kincpp)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")
# cmake_policy(SET CMP0148 NEW)


find_package(Python3 3.9 EXACT REQUIRED COMPONENTS Interpreter Development)
# set(PYBIND11_PYTHON_VERSION 3.9)
set(PYBIND11_FINDPYTHON OFF)
# find_package(Python3 3.9 EXACT REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

# find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})
find_package(Eigen3 REQUIRED)

set(LIBRARY_OUTPUT_PATH "../kincpp")

pybind11_add_module(kincpp
        src/kinematics_solver.cpp
        src/utils.cpp)

target_link_libraries(kincpp PRIVATE
    Eigen3::Eigen)

add_library(kincpp_lib STATIC
        src/kinematics_solver.cpp
        src/utils.cpp)

set_target_properties(kincpp_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(kincpp_lib PUBLIC ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(kincpp_lib PRIVATE
    Eigen3::Eigen)
