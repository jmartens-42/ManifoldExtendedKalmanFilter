cmake_minimum_required(VERSION 3.10)

find_package(Eigen3 REQUIRED)

SET(SRC_FILES src/attitude_kalman_9dof.cpp
              src/tools.cpp)


ADD_LIBRARY(ManifoldExtendedKalmanFilter ${SRC_FILES})

target_include_directories(ManifoldExtendedKalmanFilter 
                            PUBLIC
                                include/)

target_link_libraries(ManifoldExtendedKalmanFilter Eigen3::Eigen)
