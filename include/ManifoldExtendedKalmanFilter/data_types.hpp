#pragma once
#include <Eigen/Dense>

struct RotationalState{
    // attitude with respect to world frame
    Eigen::Quaternionf attitude = Eigen::Quaternionf(1, 0, 0, 0);

    // angular velocity in ego frame
    Eigen::Vector<float, 3> angular_velocity = Eigen::Vector<float, 3>(0,0,0);

};

struct Parameter{

    uint16_t parm_id;
    float parm_val;
    bool acked = false;
};