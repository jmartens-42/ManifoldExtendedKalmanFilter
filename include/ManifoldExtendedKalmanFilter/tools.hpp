#pragma once

#include <Eigen/Dense>

Eigen::Matrix<float, 3, 3> crossProdMat(const Eigen::Vector<float, 3>& v);

Eigen::Quaternionf angularVelocityToRotation(const Eigen::Vector<float, 3>& w, float dt);
