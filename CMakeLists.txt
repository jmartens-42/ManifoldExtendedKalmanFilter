cmake_minimum_required(VERSION 3.10)
SET(PROJ_NAME ManifoldExtendedKalmanFilter)
PROJECT(${PROJ_NAME})

set(CMAKE_CXX_COMPILER /usr/bin/g++-10)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3 ")
set(CMAKE_VERBOSE_MAKEFILE on)


find_package(Eigen3 REQUIRED)

include_directories(include/
                    ${EIGEN3_INCLUDE_DIRS})


SET(SRC_FILES src/attitude_kalman_9dof.cpp
              src/tools.cpp)


ADD_LIBRARY(ManifoldExtendedKalmanFilter ${SRC_FILES})