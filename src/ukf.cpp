#include "ukf.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // guessed inital value
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // guessed initial value
  std_yawdd_ = 1.0;

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

  //as given from lecture
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  Q_ = MatrixXd(n_aug_ - n_x_, n_aug_ - n_x_);
  Q_ << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_(0, 0) = std_radr_ * std_radr_;
  R_radar_(1, 1) = std_radphi_ * std_radphi_;
  R_radar_(2, 2) = std_radrd_ * std_radrd_;

  R_laser_ = MatrixXd(2, 2);
  R_laser_(0, 0) = std_laspx_ * std_laspx_;
  R_laser_(1, 1) = std_laspy_ * std_laspy_;

  sig_points_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  sig_points_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_[0] = lambda_ / (lambda_ + n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement) {

  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (measurement.sensor_type_ == MeasurementPackage::LASER) {
      std::cout <<"LASER" <<std::endl;
  } else {
    std::cout <<"RADAR" <<std::endl;
  }

  if (!is_initialized_) {

    if (measurement.sensor_type_ == MeasurementPackage::LASER) {
      x_ << measurement.raw_measurements_[0], measurement.raw_measurements_[1], 0, 0, 0;
      P_ = 0.15 * MatrixXd::Identity(5, 5);
    } else if (measurement.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement.raw_measurements_[0];
      double theta = measurement.raw_measurements_[1];
      double x = rho * cos(theta);
      double y = rho * sin(theta);
      x_ << x, y, 0, 0, 0;
      P_ = 0.03 * MatrixXd::Identity(5, 5);
    }

    is_initialized_ = true;
    time_us_ = measurement.timestamp_;

    return;
  }

  double delta_t = (measurement.timestamp_ - time_us_) / 1000000.0;
  cout << "delta " << delta_t <<"\n\n";
  Prediction(delta_t);
  time_us_ = measurement.timestamp_;

  if (measurement.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(measurement);
  }

  if (measurement.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(measurement);
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

  GenerateSigmaPoints();

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double v = sig_points_aug_(2, i);
    double yaw = sig_points_aug_(3, i);
    double yawd = sig_points_aug_(4, i);
    double nu_a = sig_points_aug_(5, i);
    double nu_yawdd = sig_points_aug_(6, i);

    VectorXd noise = VectorXd(n_x_);
    noise << .5 * delta_t * delta_t * cos(yaw) * nu_a,
             .5 * delta_t * delta_t * sin(yaw) * nu_a,
             delta_t * nu_a,
             .5 * delta_t * delta_t * nu_yawdd,
             delta_t * nu_yawdd;
    VectorXd pred = VectorXd(n_x_);

    if (fabs(yawd) < 0.0001) {
      pred << v * cos(yaw) * delta_t,
              v * sin(yaw) * delta_t,
              0,
              0,
              0;
    } else {
      pred << (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw)),
              (v / yawd) * (-cos(yaw + yawd * delta_t) + cos(yaw)),
              0,
              yawd * delta_t,
              0;
    }

    sig_points_.col(i) = sig_points_aug_.col(i).head(n_x_) + pred + noise;
  }

  x_.fill(0.0);
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ += weights_[i] * sig_points_.col(i);
  }

  //predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd diff = (sig_points_.col(i) - x_);
    diff[3] = atan2(sin(diff[3]), cos(diff[3]));
    P_ += weights_[i] * diff * diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement) {

  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = measurement.raw_measurements_;

  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px  = sig_points_(0, i);
    double py  = sig_points_(1, i);
    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }
  VectorXd Zpred = VectorXd(2);
  Zpred.fill(0);

  //calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Zpred += weights_[i] * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  MatrixXd S = R_laser_.replicate(1, 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd diff = Zsig.col(i) - Zpred;
    S += weights_[i] * diff * diff.transpose();
  }

  MatrixXd Tc = MatrixXd(n_x_, 2);
  Tc.fill(0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = (sig_points_.col(i) - x_);
    VectorXd z_diff = Zsig.col(i) - Zpred;

    x_diff[3] = atan2(sin(x_diff[3]), cos(x_diff[3]));
    z_diff[1] = atan2(sin(z_diff[1]), cos(z_diff[1]));

    Tc += weights_[i] * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - Zpred;
  z_diff[1] = atan2(sin(z_diff[1]), cos(z_diff[1]));
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement) {

  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z = measurement.raw_measurements_;

  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px  = sig_points_aug_(0, i);
    double py  = sig_points_aug_(1, i);
    double v   = sig_points_aug_(2, i);
    double yaw = sig_points_aug_(3, i);
    double ps  = sqrt((px * px) + (py * py));
    Zsig(0, i) = ps;
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / ps;
  }
  VectorXd Zpred = VectorXd(3);
  Zpred.fill(0);

  //calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Zpred += weights_[i] * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  MatrixXd S = R_radar_.replicate(1, 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd diff = Zsig.col(i) - Zpred;
    S += weights_[i] * diff * diff.transpose();
  }

  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = (sig_points_.col(i) - x_);
    VectorXd z_diff = Zsig.col(i) - Zpred;

    x_diff[3] = atan2(sin(x_diff[3]), cos(x_diff[3]));
    z_diff[1] = atan2(sin(z_diff[1]), cos(z_diff[1]));
        
    Tc += weights_[i] * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - Zpred;
  z_diff[1] = atan2(sin(z_diff[1]), cos(z_diff[1]));
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

void UKF::GenerateSigmaPoints() {

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
  P_aug.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;

  sig_points_aug_.col(0) = x_aug;
  MatrixXd A = P_aug.llt().matrixL();
  double sq_lambda = sqrt(lambda_ + n_aug_);

  for (int i = 1; i < n_aug_ + 1; i++) {
    sig_points_aug_.col(i) = x_aug + sq_lambda * A.col(i - 1);
    sig_points_aug_.col(i + n_aug_) = x_aug - sq_lambda * A.col(i - 1);
  }
}
