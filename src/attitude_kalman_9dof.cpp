#include "attitude_kalman_9dof.hpp"
#include <iostream>
#include "tools.hpp"


Eigen::Quaternionf AttitudeKalman9Dof::inverseChart(Eigen::Vector<float, 3> e){

    auto multiplier = 1/(sqrtf32(4 + e.squaredNorm()));
    auto scaled_e = multiplier*e;
    return Eigen::Quaternionf(multiplier*2, scaled_e.x(), scaled_e.y(), scaled_e.z());
}


Eigen::Vector<float, 3> AttitudeKalman9Dof::chart(Eigen::Quaternionf q){

    return 2*q.vec()/q.w();
}


RotationalState AttitudeKalman9Dof::step(Eigen::Vector<float, 9> measurement, Eigen::Vector<float, 3> expected_worldframe_acceleration_disturbance, float dt){    

    PredictState(dt);

    PredictMeasurement(expected_worldframe_acceleration_disturbance);

    updateEstimate(measurement);

    return current_estimate_;

}


void AttitudeKalman9Dof::PredictState(float dt){

    // just step forward the last angular velocity
    current_prediction_.angular_velocity = current_estimate_.angular_velocity;

    // estimate how much we rotated during the time step
    auto rotation_estimate = angularVelocityToRotation(current_prediction_.angular_velocity, dt);

    // update the current attitude estimate
    current_prediction_.attitude = current_estimate_.attitude * rotation_estimate;

    // F rotates P into current prediction frame
    Eigen::Matrix<float, 6, 6> F;
    F << rotation_estimate.toRotationMatrix().transpose(), (Eigen::Matrix<float, 3, 3>::Identity()*dt), Eigen::Matrix<float, 3, 3>::Zero(), Eigen::Matrix<float, 3, 3>::Identity();

    // Qn is process noise (how much do we trust the prediction?)
    Eigen::Matrix<float, 6, 6> Qn;
    Qn << Q_w_ * powf32(dt, 3) / 3.0, -1*Q_w_ * powf32(dt, 3) / 3.0, -1 * Q_w_ * powf32(dt, 2) / 2.0, Q_w_ * dt;

    // Rotate P into new predicted frame
    P_ = F * (P_ + Qn) * F.transpose();

}


void AttitudeKalman9Dof::PredictMeasurement(Eigen::Vector<float, 3> expected_acc_disturbance){

    // this is used a lot, just grab it now
    const auto predicted_rotation_matrix = current_prediction_.attitude.toRotationMatrix();

    // predict the measurements based on our current predicted attitude and normalize them
    const Eigen::Vector<float, 3> expected_gravity_vector = Eigen::Vector<float, 3>((expected_acc_disturbance + world_frame_acc_).transpose() * predicted_rotation_matrix).normalized()*10;
    const Eigen::Vector<float, 3> expected_mag_data = Eigen::Vector<float, 3>((expected_mag_disturbance_ + world_frame_mag_).transpose() * predicted_rotation_matrix).normalized()*10;
    const auto expected_gyro_data = current_prediction_.angular_velocity;

    expected_measurement_ << expected_gravity_vector, expected_mag_data, expected_gyro_data;

    // H transforms State -> measurement frame
    H_ << crossProdMat(expected_gravity_vector),  Eigen::Matrix<float, 3, 3>::Zero(), crossProdMat(expected_mag_data), Eigen::Matrix<float, 3, 3>::Zero(),  Eigen::Matrix<float, 3, 3>::Zero(), Eigen::Matrix<float, 3, 3>::Identity();

    Eigen::Matrix<float, 9, 9> N;

    // N is matrix representing sensor disturbances and noise (how much do we trust the sensors?)
    N << predicted_rotation_matrix.transpose() * Q_a_m_ * predicted_rotation_matrix + R_a_m_, Eigen::Matrix<float, 3, 3>::Zero(), Eigen::Matrix<float, 3, 3>::Zero(), 
            Eigen::Matrix<float, 3, 3>::Zero(), predicted_rotation_matrix.transpose() * Q_m_m_ * predicted_rotation_matrix + R_m_m_, Eigen::Matrix<float, 3, 3>::Zero(),
            Eigen::Matrix<float, 3, 3>::Zero(), Eigen::Matrix<float, 3, 3>::Zero(),  R_w_;

    // S is the measurement covariance (in measurement frame)
    S_ = H_*P_*H_.transpose() + N;

}

void AttitudeKalman9Dof::updateEstimate(Eigen::Vector<float, 9> measurement){

    // generate kalman gain matrix
    K_ = P_ * H_.transpose() * S_.inverse();

    // set up the euclidian estimate "centered at our current attitude" i.e (0, w)
    euclidian_estimate_ << Eigen::Vector<float, 3>::Zero(), current_prediction_.angular_velocity;

    // normalize the accelerometer measurement
    measurement.head(3) = (measurement.head(3)).normalized()*10.0;
    // normalize the mag measurement 
    measurement.block<3, 1>(3, 0) = measurement.block<3, 1>(3, 0).normalized()*10.0;

    // update the euclidian state based on measurement residual
    euclidian_estimate_ = euclidian_estimate_ + K_*(measurement - expected_measurement_);

    // update estimate covariance
    P_ = (Eigen::Matrix<float, 6, 6>::Identity() - K_ * H_) * P_;

    // update attitude estimate from euclidian state estimate
    current_estimate_.attitude = current_prediction_.attitude * inverseChart(euclidian_estimate_.head(3));
    current_estimate_.angular_velocity = euclidian_estimate_.tail<3>();

}


