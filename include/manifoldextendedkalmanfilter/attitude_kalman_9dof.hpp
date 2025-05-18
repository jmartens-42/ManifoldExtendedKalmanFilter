#pragma once
#include <manifoldextendedkalmanfilter/data_types.hpp>
#include <Eigen/Dense>


class AttitudeKalman9Dof{

public:

    AttitudeKalman9Dof():
                          expected_measurement_(Eigen::Vector<float, 9>::Zero()),
                          P_(Eigen::Matrix<float, 6, 6>::Identity()*10),
                          H_(Eigen::Matrix<float, 9, 6>::Zero()),
                          S_(Eigen::Matrix<float, 9, 9>::Zero()),
                          K_(Eigen::Matrix<float, 6, 9>::Zero()),
                          Q_w_(Eigen::Matrix<float, 3, 3>::Identity()*0.05),
                          Q_a_m_(Eigen::Matrix<float, 3, 3>::Identity()*0.0001),
                          Q_m_m_(Eigen::Matrix<float, 3, 3>::Identity()*0.0001),
                          R_w_(Eigen::Matrix<float, 3, 3>::Identity()*0.0001),
                          R_a_m_(Eigen::Matrix<float, 3, 3>::Identity()*0.0001),
                          R_m_m_(Eigen::Matrix<float, 3, 3>::Identity()*0.0001),
                          world_frame_acc_(Eigen::Vector<float, 3>({0, 0, 9.81})),
                          gyro_bias_(Eigen::Vector<float, 3>::Zero()),
                          world_frame_mag_(28.71, 0.3, -53.83),
                          expected_mag_disturbance_(Eigen::Vector<float, 3>::Zero()),
                          euclidian_estimate_(Eigen::Vector<float, 6>::Zero()){

    }

    ~AttitudeKalman9Dof(){

    }

    void rezero(){
        calibrated_ = false;
    }

    std::string getDebugString(){
        return std::string() + "expected: " + std::to_string(expected_measurement_[0]) + " " + std::to_string(expected_measurement_[1]) + " " + std::to_string(expected_measurement_[2]) + " " + std::to_string(expected_measurement_.block<3,1>(0,0).norm()) +
                                " actual: " + std::to_string(normalized_measurement_.block<3, 1>(0, 0)[0]) + " " + std::to_string(normalized_measurement_.block<3, 1>(0, 0)[1]) + " " + std::to_string(normalized_measurement_.block<3, 1>(0, 0)[2]) + " " + std::to_string(normalized_measurement_.block<3, 1>(0, 0).norm()) +
                                " expected: " + std::to_string(expected_measurement_[3]) + " " + std::to_string(expected_measurement_[4]) + " " + std::to_string(expected_measurement_[5]) + " " + std::to_string(expected_measurement_.block<3,1>(3,0).norm())+
                                " actual: " + std::to_string(normalized_measurement_.block<3,1>(3, 0)[0]) + " " + std::to_string(normalized_measurement_.block<3,1>(3, 0)[1]) + " " + std::to_string(normalized_measurement_.block<3,1>(3, 0)[2]) + " " + std::to_string(normalized_measurement_.block<3,1>(3, 0).norm()) + "\n\r";

    }

    void setParameter(const Parameter& p){

        switch(p.parm_id){

            case 1:
                Q_w_ = Eigen::Matrix<float, 3, 3>::Identity()*p.parm_val;
                break;
            case 2:
                Q_a_m_ = Eigen::Matrix<float, 3, 3>::Identity()*p.parm_val;
                break;
            case 3:
                Q_m_m_ = Eigen::Matrix<float, 3, 3>::Identity()*p.parm_val;
                break;
            case 4:
                R_w_ = Eigen::Matrix<float, 3, 3>::Identity()*p.parm_val;
                break;
            case 5:
                R_a_m_ = Eigen::Matrix<float, 3, 3>::Identity()*p.parm_val;
                break;
            case 6:
                R_m_m_ = Eigen::Matrix<float, 3, 3>::Identity()*p.parm_val;
                break;
            default:
                break;
        }
    }

    // returns a new estimate of the current attitude and angular velocity
    // measaurement: Eigen::Vector of accelerometer, magnetometer, gyro data
    // expected_worldframe_acceleration_deviation, a prediction of how much external acceleration is contributing to the acceleration measurement (in world frame!)
    // dt: seconds between this step and the last
    RotationalState step(const Eigen::Vector<float, 9>& measurement, const Eigen::Vector<float, 3>& expected_worldframe_acceleration_deviation, float dt);
    RotationalState step(float dt);
    Eigen::Vector<float, 9> expected_measurement_; // storage for the expected measurement in the predictMeasurement step

private:

    // Inverse chart from Euclidian R^3 -> Quaternion manifold S^3
    Eigen::Quaternionf inverseChart(const Eigen::Vector<float, 3>& e);

    // Rodrigues Parameters chart from Quaternion manifold S^3 -> Euclidian R^3
    Eigen::Vector<float, 3> chart(const Eigen::Quaternionf& q);

    void PredictState(float dt);
    void PredictMeasurement(const Eigen::Vector<float, 3>& expected_acc_disturbance);
    void updateEstimate(const Eigen::Vector<float, 9>& measurement);

    void calibrationStep(const Eigen::Vector<float, 9>& measurement);


    Eigen::Matrix<float, 6, 6> P_; //state covariance matrix
    Eigen::Matrix<float, 9, 6> H_; // state -> measurement transformation matrix
    Eigen::Matrix<float, 9, 9> S_; // measurement covariance matrix
    Eigen::Matrix<float, 6, 9> K_; // Kalman gain matrix

    Eigen::Matrix<float, 3, 3> Q_w_;   //process covariance - how much do we trust the dynamical model
    Eigen::Matrix<float, 3, 3> Q_a_m_; //covariance representing external disturbances
    Eigen::Matrix<float, 3, 3> Q_m_m_; //covariance representing external disturbances

    Eigen::Matrix<float, 3, 3> R_w_; //gyro measurement noise covariance
    Eigen::Matrix<float, 3, 3> R_a_m_; // accelerometer measurement noise covariance
    Eigen::Matrix<float, 3, 3> R_m_m_; // magnet measurement noise covariance

    Eigen::Vector<float, 3> world_frame_acc_; // this is what measurement predictions are made off of
    Eigen::Vector<float, 3> gyro_bias_;
    Eigen::Vector<float, 3> world_frame_mag_; // this is what measurement predictions are made off of

    const Eigen::Vector<float, 3> expected_mag_disturbance_; // just zero for now. Not sure how this would get implemented dynamically

    RotationalState current_prediction_; // used in the step logic
    RotationalState current_estimate_; // actual output

    Eigen::Vector<float, 9> normalized_measurement_;

    bool calibrated_ = false;

    Eigen::Vector<float, 6> euclidian_estimate_; // used as the R^3 analog of our attitude and angular velocity measurement

    enum class CalibrationStep{
        SensorBaseline,
        Complete
    };
    CalibrationStep cal_step_ = CalibrationStep::SensorBaseline;


};
