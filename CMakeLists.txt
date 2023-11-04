cmake_minimum_required(VERSION 3.16)

find_package(Eigen3 REQUIRED)

SET(SRC_FILES src/attitude_kalman_9dof.cpp
              src/tools.cpp)


ADD_LIBRARY(ManifoldExtendedKalmanFilter ${SRC_FILES})

target_compile_definitions(ManifoldExtendedKalmanFilter
    PUBLIC 
        EIGEN_NO_MALLOC
        EIGEN_DONT_PARALLELIZE
        EIGEN_NO_DEBUG
        EIGEN_NO_AUTOMATIC_RESIZING
        EIGEN_NO_IO
        EIGEN_DONT_VECTORIZE

)

target_include_directories(ManifoldExtendedKalmanFilter 
                            PUBLIC
                                include/)

target_link_libraries(ManifoldExtendedKalmanFilter Eigen3::Eigen)
