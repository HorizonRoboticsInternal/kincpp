cmake_minimum_required(VERSION 3.10)
project(kincpp)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations -fPIC")
set(CMAKE_CXX_STANDARD 17)

# Use Python from venv if one is active
if(DEFINED ENV{VIRTUAL_ENV})
    set(Python3_EXECUTABLE $ENV{VIRTUAL_ENV}/bin/python)
    message(STATUS "Using Python from: $ENV{VIRTUAL_ENV}")
endif()

# Find required packages
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

set(SOURCES
    src/kinematics_solver.cpp
    src/utils.cpp
)

# Static library
add_library(kincpp_lib STATIC ${SOURCES})
target_include_directories(kincpp_lib PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(kincpp_lib PUBLIC Eigen3::Eigen)

# Pybind module
add_library(kincpp MODULE ${SOURCES})
target_include_directories(kincpp PRIVATE
    ${Python3_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
)
target_link_libraries(kincpp PRIVATE
    Eigen3::Eigen
    pybind11::headers
)

# Set module properties
set_target_properties(kincpp PROPERTIES
    PREFIX ""
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/kincpp"
)
