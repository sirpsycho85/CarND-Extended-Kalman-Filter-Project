#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // state transition
  // will be modified later with dt based on measurement
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  // observation model - laser
  H_laser_ = MatrixXd(2, 4);
  H_laser_ <<  1, 0, 0, 0,
              0, 1, 0, 0;

  // observation noise covariance - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ <<  0.0225, 0,
              0, 0.0225;

  // observation model - radar
  // Jacobian calculated with measurement
  Hj_ = MatrixXd(3, 4);

  // observation noise covariance - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ <<  0.0225, 0, 0,
               0, 0.0225, 0,
               0, 0, 0.0225;

  noise_ax = 5;
  noise_ay = 5;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "Kalman Filter Initialization " << endl;
    cout << "EKF: " << endl;

    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    P_ = MatrixXd(4, 4);
    P_ <<  1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1000, 0,
           0, 0, 0, 1000;

    ekf_.P_ = P_;
    ekf_.F_ = F_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "INITIALIZE RADAR" << endl;

      float ro;
      float theta;
      float ro_dot;

      ro = measurement_pack.raw_measurements_[0];
      theta = measurement_pack.raw_measurements_[1];
      ro_dot = measurement_pack.raw_measurements_[2];
      
      float px;
      float py;
      float vx;
      float vy;

      px = ro * cos(theta);
      py = ro * sin(theta);

      // TODO: fix the vx and vy initialization
      vx = ro_dot * cos(theta);
      vy = ro_dot * cos(theta);

      ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "INITIALIZE LASER" << endl;
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     dt seems really tiny so Q has almost no effect...
     Q seems to have no effect anyway...can << be used to update or not?
   */
  float dt;
  cout << "dt = " << dt << endl;
  dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // process noise covariance
  ekf_.Q_ = MatrixXd(4, 4);

  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Q_ << 100,0,1,0,
             0,1,0,1,
             1,0,100,0,
             0,1,0,1;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.

     Use Jacobian for the radar. Currently just returning

   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    return;
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
