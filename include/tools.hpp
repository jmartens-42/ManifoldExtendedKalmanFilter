#pragma once
#include <eigen3/Eigen/Dense>

Eigen::Matrix<float, 3, 3> crossProdMat(Eigen::Vector<float, 3> v);

Eigen::Quaternionf angularVelocityToRotation(Eigen::Vector<float, 3> w, float dt);