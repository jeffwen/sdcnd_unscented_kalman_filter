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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  /**DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.**/
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
  /**DO NOT MODIFY**/
  
  // set the initialized flag for the measurement
  is_initialized_ = false;
  
  // set the initial process covariance matrix as expected difference between true and initialized state
  VectorXd diag = VectorXd(5);
  diag << 0.15, 0.15, 1.0, 1.0, 1.0;
  P_ = diag.asDiagonal();
  
  ///* State dimension
  n_x_ = 5;
  
  ///* Augmented state dimension
  n_aug_ = 7;
  
  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;
  
  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  
  //create weight matrix
  weights = VectorXd(2*n_aug_+1);
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Make sure you switch between lidar and radar measurements.
  **/
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    
    // first timestamp needs to be the first in the measurement NOT 1970...
    time_us_ = meas_package.timestamp_;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      
      // convert from radar measurements (polar) to cartesian and initialize state
      x_ << rho*cos(phi), rho*sin(phi), 3.0, 0.1, 0.1;
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      
      // initialize the state
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 3.0, 0.1, 0.1;
      
    }
    
    // make sure that there are no 0s in the state vector
    if (fabs(x_(0)) < 0.001 and fabs(x_(1)) < 0.001) {
      x_(0) = 0.001;
      x_(1) = 0.001;
    }
    
    is_initialized_ = true;
    return;
  }
  
  // predict
  // compute the time elapsed between the current and previous measurements
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;    //delta_t - expressed in seconds
  time_us_ = meas_package.timestamp_;
  
  Prediction(delta_t);
  
  // update
  if ((use_radar_) && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
    UpdateRadar(meas_package);
  } else if ((use_laser_) && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  **/
  
  /***SIGMA POINT GENERATION***/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  //augmented sigma point matrix with process noise
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  
  //create augmented covariance matrix
  P_aug.topLeftCorner(P_.rows(),P_.cols()) = P_;
  
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_*std_a_, 0,
  0, std_yawdd_*std_yawdd_;
  
  P_aug.bottomRightCorner(n_aug_-n_x_,n_aug_-n_x_) = Q;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0,1,n_aug_,n_aug_) = (sqrt(lambda_ + n_aug_) * A).colwise() + x_aug;
  Xsig_aug.block(0,n_aug_+1,n_aug_,n_aug_) = (-1 * (sqrt(lambda_ + n_aug_) * A)).colwise() + x_aug;
  
  /***SIGMA POINT PREDICTION***/
  for (int i=0; i<2*n_aug_+1; i++){
    //extract the x state vector without noise
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double psi = Xsig_aug(3,i);
    double psi_dot = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_psi_a = Xsig_aug(6,i);
    
    VectorXd temp_vect = VectorXd(5);
    
    //predict sigma points and check for zeros
    if(fabs(psi_dot) < 0.001){
      temp_vect << px + v*cos(psi)*delta_t + 0.5*(delta_t*delta_t)*cos(psi)*nu_a,
        py + v*sin(psi)*delta_t + 0.5*(delta_t*delta_t)*sin(psi)*nu_a,
        v + 0 + (delta_t*nu_a),
        psi + (psi_dot*delta_t) + 0.5*(delta_t*delta_t)*nu_psi_a,
        psi_dot + 0 + (delta_t*nu_psi_a);
      
    } else {
      temp_vect << px + (v/psi_dot)*(sin(psi+psi_dot*delta_t)-sin(psi)) + 0.5*(delta_t*delta_t)*cos(psi)*nu_a,
        py + (v/psi_dot)*(-1*cos(psi+psi_dot*delta_t)+cos(psi)) + 0.5*(delta_t*delta_t)*sin(psi)*nu_a,
        v + 0 + (delta_t*nu_a),
        psi + (psi_dot*delta_t) + 0.5*(delta_t*delta_t)*nu_psi_a,
        psi_dot + 0 + (delta_t*nu_psi_a);
    }

    Xsig_pred_.col(i) = temp_vect;
  }
  
  /***MEAN STATE AND COVARIANCE MATRIX PREDICTION***/
  // predicted mean and covariance matrix calculation
  x_.fill(0.0);
  
  for (int i = 0; i<2*n_aug_+1; i++){
    //set weights
    if (i==0){
      weights(i) = lambda_/(lambda_+n_aug_);
    } else {
      weights(i) = 1/(2*(lambda_+n_aug_));
    }
    //predict state mean
    x_ += weights(i)*Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i<2*n_aug_+1; i++){
    //predict state covariance matrix
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) {x_diff(3)-=2.*M_PI;};
    while (x_diff(3)<-M_PI) {x_diff(3)+=2.*M_PI;};
    P_ += weights(i) * x_diff * x_diff.transpose();
  }
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  **/
  
  /***PREDICT MEASUREMENT***/
  //set measurement dimension, lidar can measure px and py
  int n_z_ = 2;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  
  // measurement noise matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R.fill(0.0);
  
  VectorXd temp_vect = VectorXd(n_z_);
  
  for (int i=0; i<2*n_aug_+1; i++){
    float p_x = Xsig_pred_(0,i);
    float p_y = Xsig_pred_(1,i);
    
    //transform sigma points into measurement space
    temp_vect << p_x, p_y;
    Zsig.col(i) = temp_vect;
    
    //calculate mean predicted measurement
    z_pred += weights(i)*Zsig.col(i);
    
  }

  for (int i=0; i<2*n_aug_+1; i++){
    
    //calculate innovation covariance matrix S
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    S += weights(i) * z_diff * z_diff.transpose();
  }
  
  R(0,0) = pow(std_laspx_,2);
  R(1,1) = pow(std_laspy_,2);
  S += R;
  
  /***UPDATE STATE***/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
    
    //calculate cross correlation matrix
    VectorXd x_diff = Xsig_pred_.col(i)-x_;
    VectorXd z_diff = Zsig.col(i)-z_pred;
    
    Tc += weights(i)*(x_diff)*(z_diff).transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  
  // calculate NIS for radar
  float nis_lidar = (meas_package.raw_measurements_-z_pred).transpose()*S.inverse()*(meas_package.raw_measurements_-z_pred);
  lidar_nis_vect.push_back(nis_lidar);
  
  // calculate percentage of NIS over the Chi-squred value we care about (only 5% should be greater than this balue given 2 degrees of freedom px and py)
  float nis_lidar_percentage = std::count_if(lidar_nis_vect.begin(), lidar_nis_vect.end(), [](float i){return i > 5.991;});
  std::cout<<"Lidar NIS: "<<nis_lidar_percentage/lidar_nis_vect.size()<<"\n";

  //file writer
//  ofstream file;
//  file.open("/Users/jwen/Courses/sdcnd/sdcnd_unscented_kalman_filter/etc/nis_measurements.txt", std::ofstream::out | std::ofstream::app);
//  file<<"lidar,"<<nis_lidar<<"\n";
//  file.close();
  
  //update state mean and covariance matrix
  x_ += K*(meas_package.raw_measurements_-z_pred);
  P_ -= K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  **/

  /***PREDICT MEASUREMENT***/
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  
  // measurement noise matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R.fill(0.0);
  
  VectorXd temp_vect = VectorXd(n_z_);
  
  for (int i=0; i<2*n_aug_+1; i++){
    float p_x = Xsig_pred_(0,i);
    float p_y = Xsig_pred_(1,i);
    float nu = Xsig_pred_(2,i);
    float psi = Xsig_pred_(3,i);
    float psi_dot = Xsig_pred_(4,i);
    float sqrt_px_py = sqrt(pow(p_x,2)+pow(p_y,2));
    
    //check no zeros
    if(fabs(sqrt_px_py) < 0.0001){
      sqrt_px_py = 0.0001;
    }
    
    //transform sigma points into measurement space
    temp_vect << sqrt_px_py, atan2(p_y, p_x), (p_x*cos(psi)*nu+p_y*sin(psi)*nu)/sqrt_px_py;
    Zsig.col(i) = temp_vect;
    
    //calculate mean predicted measurement
    z_pred += weights(i)*Zsig.col(i);
  
  }
  
  for (int i=0; i<2*n_aug_+1; i++){
    
    //calculate innovation covariance matrix S
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    S += weights(i) * z_diff * z_diff.transpose();
  }
  
  R(0,0) = pow(std_radr_,2);
  R(1,1) = pow(std_radphi_,2);
  R(2,2) = pow(std_radrd_,2);
  S += R;
  
  /***UPDATE STATE***/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
    
    //calculate cross correlation matrix
    VectorXd x_diff = Xsig_pred_.col(i)-x_;
    VectorXd z_diff = Zsig.col(i)-z_pred;
    
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //check no zeros
    if(fabs(z_diff(1)) < 0.0001){
      z_diff(1) = 0.0001;
    }
    
    Tc += weights(i)*(x_diff)*(z_diff).transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  
  // calculate NIS for radar
  float nis_radar = (meas_package.raw_measurements_-z_pred).transpose()*S.inverse()*(meas_package.raw_measurements_-z_pred);
  radar_nis_vect.push_back(nis_radar);
  
  // calculate percentage of NIS over the Chi-squred value we care about (only 5% should be greater than this balue given 3 degrees of freedom px and py)
  float nis_radar_percentage = std::count_if(radar_nis_vect.begin(), radar_nis_vect.end(), [](float i){return i > 7.815;});
  std::cout<<"Radar NIS: "<<nis_radar_percentage/radar_nis_vect.size()<<"\n";

  //file writer
//  ofstream file;
//  file.open("/Users/jwen/Courses/sdcnd/sdcnd_unscented_kalman_filter/etc/nis_measurements.txt", std::ofstream::out | std::ofstream::app);
//  file<<"radar,"<<nis_radar<<"\n";
//  file.close();
  
  //update state mean and covariance matrix
  x_ += K*(meas_package.raw_measurements_-z_pred);
  P_ -= K*S*K.transpose();
}
