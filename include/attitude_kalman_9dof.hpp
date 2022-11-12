#pragma once
#include "data_types.hpp"
#include <eigen3/Eigen/Dense>


class AttitudeKalman9Dof{

public:

    AttitudeKalman9Dof(){

    }

    ~AttitudeKalman9Dof(){

    }

    // returns a new estimate of the current attitude and angular velocity
    // measaurement: Eigen::Vector of accelerometer, magnetometer, gyro data
    // expected_worldframe_acceleration_deviation, a prediction of how much external acceleration is contributing to the acceleration measurement (in world frame!)
    // dt: seconds between this step and the last
    RotationalState step(Eigen::Vector<float, 9> measurement, Eigen::Vector<float, 3> expected_worldframe_acceleration_deviation, float dt);


private:

    // Inverse chart from Euclidian R^3 -> Quaternion manifold S^3
    Eigen::Quaternionf inverseChart(Eigen::Vector<float, 3> e);

    // Rodrigues Parameters chart from Quaternion manifold S^3 -> Euclidian R^3
    Eigen::Vector<float, 3> chart(Eigen::Quaternionf q);

    void PredictState(float dt);
    void PredictMeasurement(Eigen::Vector<float, 3> expected_acc_disturbance);
    void updateEstimate(Eigen::Vector<float, 9> measurement);


    Eigen::Matrix<float, 6, 6> P_ = Eigen::Matrix<float, 6, 6>::Identity()*10; //state covariance matrix
    Eigen::Matrix<float, 9, 6> H_ = Eigen::Matrix<float, 9, 6>::Zero(); // state -> measurement transformation matrix
    Eigen::Matrix<float, 9, 9> S_ = Eigen::Matrix<float, 9, 9>::Zero(); // measurement covariance matrix
    Eigen::Matrix<float, 6, 9> K_ = Eigen::Matrix<float, 6, 9>::Zero(); // Kalman gain matrix 

    Eigen::Matrix<float, 3, 3> Q_w_ = Eigen::Matrix<float, 3, 3>::Identity()*0.6; //process covariance
    Eigen::Matrix<float, 3, 3> Q_a_m_ = Eigen::Matrix<float, 3, 3>::Identity()*0.1; //covariance representing external disturbances
    Eigen::Matrix<float, 3, 3> Q_m_m_ = Eigen::Matrix<float, 3, 3>::Identity()*0.1; //covariance representing external disturbances

    Eigen::Matrix<float, 3, 3> R_w_ = Eigen::Matrix<float, 3, 3>::Identity()*0.1; //gyro measurement noise covariance
    Eigen::Matrix<float, 3, 3> R_a_m_ = Eigen::Matrix<float, 3, 3>::Identity()*0.1; // accelerometer measurement noise covariance
    Eigen::Matrix<float, 3, 3> R_m_m_ = Eigen::Matrix<float, 3, 3>::Identity()*0.1; // magnet measurement noise covariance

    const Eigen::Vector<float, 3> world_frame_acc_ = Eigen::Vector<float, 3>({0, 0, 9.81}); // this is what measurement predictions are made off of
    
    Eigen::Vector<float, 3> expected_mag_disturbance_ = Eigen::Vector<float, 3>::Zero(); // just zero for now. Not sure how this would get implemented dynamically
    const Eigen::Vector<float, 3> world_frame_mag_ = Eigen::Vector<float, 3>(0, 0.707, -0.707); // this is what measurement predictions are made off of
    
    RotationalState current_prediction_; // used in the step logic
    RotationalState current_estimate_; // actual output

    Eigen::Vector<float, 9> expected_measurement_ = Eigen::Vector<float, 9>::Zero(); // storage for the expected measurement in the predictMeasurement step

    Eigen::Vector<float, 6> euclidian_estimate_ = Eigen::Vector<float, 6>::Zero(); // used as the R^3 analog of our attitude and angular velocity measurement


};