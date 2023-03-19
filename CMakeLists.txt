cmake_minimum_required(VERSION 3.10)

find_package(Eigen3 REQUIRED)

include_directories(include/)


SET(SRC_FILES src/attitude_kalman_9dof.cpp
              src/tools.cpp)


ADD_LIBRARY(ManifoldExtendedKalmanFilter ${SRC_FILES})
target_link_libraries(ManifoldExtendedKalmanFilter Eigen3::Eigen)