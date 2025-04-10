cmake_minimum_required(VERSION 3.10)
project(kincpp)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

message(STATUS "Using compiler ${CMAKE_CXX_COMPILER_ID} ...")

find_package(pybind11 REQUIRED PATHS ${CMAKE_PREFIX_PATH} NO_DEFAULT_PATH)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(Eigen3 REQUIRED)

include_directories(${Python3_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")

set(LIBRARY_OUTPUT_PATH "../kincpp")

add_library(kincpp_lib STATIC
        src/kinematics_solver.cpp
        src/utils.cpp)
target_link_libraries(kincpp_lib PRIVATE
    pybind11::pybind11 Eigen3::Eigen)

set_target_properties(kincpp_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

pybind11_add_module(kincpp
    src/kinematics_solver.cpp
    src/utils.cpp)
