#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  time_us_ = 0;
  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  
  x_ = VectorXd(n_x_);
  x_.fill(0.0);
  
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);
  
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  nis_radar = 0.0;
  nis_lidar = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
 if (!is_initialized_) {
   double px = meas_pack.raw_measurements_[0];
   double py = meas_pack.raw_measurements_[1];
   x_ << px, py, 0, 0, 0;

   P_ << 0.15, 0, 0, 0, 0,
   		 0, 0.15, 0, 0, 0,
   		 0, 0, 1, 0, 0,
   		 0, 0, 0, 1, 0,
   		 0, 0, 0, 0, 1;
   
   time_us_ = meas_pack.timestamp_;
   is_initialized_ = true;
   return;
 }

 double delta_t = (meas_pack.timestamp_ - time_us_) / 1000000.0;
 time_us_ = meas_pack.timestamp_;
  
 Prediction(delta_t);
  
 if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
   UpdateRadar(meas_pack);
 } else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
   UpdateLidar(meas_pack);
 }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Create augmented mean state
  VectorXd x_aug(7);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;
  
  // Create augmented covariance matrix
  MatrixXd P_aug(7, 7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  // Square root matrix
  MatrixXd L = P_aug.llt().matrixL();
    
  // Create augmented sigma points
  MatrixXd Xsig_aug(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(n_aug_ + i + 1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  
  // Predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p = 0.0;
    double py_p = 0.0;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p += nu_a * delta_t;

    yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
    
  // Predicted state mean
  VectorXd x_pred(n_x_);
  x_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_pred += weights_(i) * Xsig_pred.col(i); 
  }
  
  // Predicted state covariance matrix
  MatrixXd P_pred(n_x_, n_x_);
  P_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // Calculate the difference vector from predicted
    // sigma points to predicted mean
    VectorXd diff = Xsig_pred.col(i) - x_pred;
    // Normalization
    if (diff(3) > 100 || diff(3) < -100) {
      printf("--x_(%d)--\n", i);
      std::cout << x_ << std::endl;
    }
    while (diff(3) > M_PI) diff(3) -= 2.0 * M_PI;
    while (diff(3) < -M_PI) diff(3) += 2.0 * M_PI;
    
    P_pred += weights_(i) * diff * diff.transpose();
  }
    
  x_ = x_pred;
  P_ = P_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  MatrixXd Zsig(2, 2 * n_aug_ + 1);
  for(int i = 0; i < Zsig.cols(); i++) {
    Zsig(0, i) = Xsig_pred(0, i); 
    Zsig(1, i) = Xsig_pred(1, i); 
  }
  
  VectorXd z_pred(2);
  z_pred.fill(0.0);
  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred += weights_(i) * Zsig.col(i); 
  }
  
  MatrixXd S(2, 2);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd diff = Zsig.col(i) - z_pred;
    S += weights_(i) * diff * diff.transpose();
  }
  
  // lidar measurement noise
  MatrixXd R(2, 2);
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  S += R;
  
  int n_z_ = 2;
  
  VectorXd z = meas_package.raw_measurements_;
  
  MatrixXd Tc(n_x_, n_z_);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z - z_pred;
  
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
  
  nis_lidar = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  MatrixXd Zsig(3, 2 * n_aug_ + 1);
  // Convert sigma points to measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // extract values for better readability
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;
    
    if (fabs(p_x) < 0.0001 || fabs(p_x * p_x + p_y * p_y) < 0.0001) {
      cout << "updateRadar () - Error - Division by Zero" << endl;
      continue;
    }
    
    // measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                       	// r
    Zsig(1, i) = atan2(p_y, p_x);                                		// phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);	// r_dot
  }
  
  VectorXd z_pred(3);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  MatrixXd S(3, 3);
  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // residual
    VectorXd diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (diff(1) > M_PI) diff(1) -= 2.0 * M_PI;
    while (diff(1) < -M_PI) diff(1) += 2.0 * M_PI;

    S += weights_(i) * diff * diff.transpose();
  }
  
  // add measurement noise covariance matrix
  MatrixXd R(3,3);
  R <<  std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;
  S += R;
  
  const int n_z_ = 3;
  
  // Create z from raw_measurements
  double ro  = meas_pack.raw_measurements_[0];
  double theta = meas_pack.raw_measurements_[1];
  double ro_dot = meas_pack.raw_measurements_[2];
  VectorXd z(n_z_);
  z << ro, theta, ro_dot;
  
  // matrix for cross correlation Tc
  MatrixXd Tc(n_x_, n_z_);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.0 * M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2.0 * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.0 * M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
  
  // NIS for radar
  nis_radar = z_diff.transpose() * S.inverse() * z_diff;
}
