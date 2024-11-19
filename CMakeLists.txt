cmake_minimum_required(VERSION 3.10)
project(kincpp)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)
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
