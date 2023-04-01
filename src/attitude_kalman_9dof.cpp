#include "attitude_kalman_9dof.hpp"
#include "tools.hpp"


Eigen::Quaternionf AttitudeKalman9Dof::inverseChart(const Eigen::Vector<float, 3>& e){

    auto multiplier = 1/(sqrtf(4 + e.squaredNorm()));
    Eigen::Vector<float, 3> scaled_e = multiplier*e;
    return Eigen::Quaternionf(multiplier*2, scaled_e.x(), scaled_e.y(), scaled_e.z()).normalized();
}


Eigen::Vector<float, 3> AttitudeKalman9Dof::chart(const Eigen::Quaternionf& q){

    return 2*q.vec()/q.w();
}


RotationalState AttitudeKalman9Dof::step(const Eigen::Vector<float, 9>& measurement, const Eigen::Vector<float, 3>& expected_worldframe_acceleration_disturbance, float dt){    

    if(calibrated_){

        PredictState(dt);

        PredictMeasurement(expected_worldframe_acceleration_disturbance);

        updateEstimate(measurement);

    }
    else{
        calibrationStep(measurement);
    }

    return current_estimate_;

}

void AttitudeKalman9Dof::calibrationStep(const Eigen::Vector<float, 9>& measurement){

    static uint32_t sample_count = 0;
    static Eigen::Vector<float, 3> accel_temp(0,0,0);
    static Eigen::Vector<float, 3> gyro_temp(0,0,0);
    static Eigen::Vector<float, 3> mag_temp(0,0,0);

    switch(cal_step_){

        case CalibrationStep::SensorBaseline:

            accel_temp += measurement.block<3, 1>(0,0);
            gyro_temp += measurement.block<3, 1>(6,0);
            mag_temp += measurement.block<3, 1>(3,0);

            sample_count++;
            if(sample_count > 50){
                world_frame_acc_ = accel_temp / 50.0;
                gyro_bias_ = gyro_temp / 50.0;
                world_frame_mag_ = mag_temp / 50.0;
                cal_step_ = CalibrationStep::Complete;
            }

            break;
        case CalibrationStep::Complete:
            calibrated_ = true;
            break;
        default:
            break;
    }


}


void AttitudeKalman9Dof::PredictState(float dt){

    // just step forward the last angular velocity
    current_prediction_.angular_velocity = current_estimate_.angular_velocity;

    // estimate how much we rotated during the time step
    Eigen::Quaternionf rotation_estimate = angularVelocityToRotation(current_prediction_.angular_velocity, dt);

    // update the current attitude estimate
    current_prediction_.attitude = current_estimate_.attitude * rotation_estimate;

    // F rotates P into current prediction frame
    Eigen::Matrix<float, 6, 6> F;
    F << rotation_estimate.toRotationMatrix().transpose(), (Eigen::Matrix<float, 3, 3>::Identity()*dt), Eigen::Matrix<float, 3, 3>::Zero(), Eigen::Matrix<float, 3, 3>::Identity();

    // Qn is process noise (how much do we trust the prediction?)
    Eigen::Matrix<float, 6, 6> Qn;
    Qn << Q_w_ * powf(dt, 3) / 3.0, -1*Q_w_ * powf(dt, 3) / 3.0, -1 * Q_w_ * powf(dt, 2) / 2.0, Q_w_ * dt;

    // Rotate P into new predicted frame
    P_ = F * (P_ + Qn) * F.transpose();
}


void AttitudeKalman9Dof::PredictMeasurement(const Eigen::Vector<float, 3>& expected_acc_disturbance){

    // this is used a lot, just grab it now
    const Eigen::Matrix<float, 3, 3> predicted_rotation_matrix = current_prediction_.attitude.toRotationMatrix();

    // predict the measurements based on our current predicted attitude and normalize them
    expected_measurement_.block<3, 1>(0,0) = Eigen::Vector<float, 3>((expected_acc_disturbance + world_frame_acc_).transpose() * predicted_rotation_matrix).normalized();
    expected_measurement_.block<3, 1>(3,0) = Eigen::Vector<float, 3>((expected_mag_disturbance_ + world_frame_mag_).transpose() * predicted_rotation_matrix ).normalized();
    expected_measurement_.block<3, 1>(6,0) = current_prediction_.angular_velocity;

    // H transforms State -> measurement frame
    H_.block<3, 3>(0,0) = crossProdMat(expected_measurement_.block<3, 1>(0,0));
    H_.block<3, 3>(3,0) = crossProdMat(expected_measurement_.block<3, 1>(3,0));
    H_.block<3, 3>(6,3) = Eigen::Matrix<float, 3, 3>::Identity();

    // N is matrix representing sensor disturbances and noise (how much do we trust the sensors?)
    Eigen::Matrix<float, 9, 9> N = Eigen::Matrix<float, 9, 9>::Zero();
    N.block<3, 3>(0, 0) = predicted_rotation_matrix.transpose() * Q_a_m_ * predicted_rotation_matrix + R_a_m_;
    N.block<3, 3>(3, 3) = predicted_rotation_matrix.transpose() * Q_m_m_ * predicted_rotation_matrix + R_m_m_;
    N.block<3, 3>(6, 6) = R_w_;

    // S is the measurement covariance (in measurement frame)
    S_ = H_*P_*H_.transpose() + N;
}

void AttitudeKalman9Dof::updateEstimate(const Eigen::Vector<float, 9>& measurement){

    // generate kalman gain matrix
    K_ = (P_ * H_.transpose()) * S_.inverse();

    // set up the euclidian estimate "centered at our current attitude" i.e (0, w)
    euclidian_estimate_ << Eigen::Vector<float, 3>::Zero(), current_prediction_.angular_velocity;

    Eigen::Vector<float, 9> normalized_measurement;
    // normalize the accelerometer measurement
    normalized_measurement_.block<3, 1>(0, 0) = measurement.block<3, 1>(0, 0).normalized();

    // normalize the mag measurement 
    normalized_measurement_.block<3, 1>(3, 0) = measurement.block<3, 1>(3,0).normalized();

    // de-bias the gyro measurement
    normalized_measurement_.block<3,1>(6, 0) =  measurement.block<3, 1>(6, 0) - gyro_bias_;

    // update the euclidian state based on measurement residual
    euclidian_estimate_ += K_*(normalized_measurement_ - expected_measurement_);

    // update estimate covariance
    P_ = (Eigen::Matrix<float, 6, 6>::Identity() - K_ * H_) * P_;

    // update attitude estimate from euclidian state estimate
    current_estimate_.attitude = (current_prediction_.attitude * inverseChart(euclidian_estimate_.head(3))).normalized();
    current_estimate_.angular_velocity = euclidian_estimate_.tail<3>(); //normalized_measurement.block<3,1>(6, 0); //euclidian_estimate_.tail<3>(); 

}


