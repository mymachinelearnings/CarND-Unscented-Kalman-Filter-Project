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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
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
  n_x_ = 5;
  n_aug_ = 7;
  
  is_initialized_ = false;
  Xsig_pred_ = MatrixXd(7, 15);
  P_ = MatrixXd(5,5);
  weights_ = VectorXd(7);

  P_ << 1000, 0, 0, 0, 0, 
        0, 1000, 0, 0, 0, 
        0, 0, 1000, 0, 0, 
        0, 0, 0, 1000, 0, 
        0, 0, 0, 0, 1000;

  //Process noise covariance matrix
  MatrixXd Q_ = MatrixXd(2,2);
  Q_ << std_a_*std_a_, 0, 
        0, std_yawdd_*std_yawdd_;

  //Measurement noise covariance matrix
  Matrix R_ = MatrixXd(3,3);
  R_ << std_radr_*std_radr_, 0, 0, 
        0, std_radphi_*std_radphi_, 0, 
        0, 0, std_radrd_*std_radrd_;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_) {
    time_us_ = meas_package.timestamp_;



    if(meas_package.sensor_type_ == MeasurementPackage::'LASER') {
      x_ << meas_package.raw_measurements_[1], meas_package.raw_measurements_[2], 0.0, 0.0, 0.0;
    }
    if(meas_package.raw_measurements_[0] == MeasurementPackage::'RADAR') {
      float ro     = measurement_pack.raw_measurements_[0];
      float theta  = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      float px_r = ro * cos(theta);
      float py_r = ro * sin(theta);
      float vx_r = ro_dot * cos(theta);
      float vy_r = ro_dot * sin(theta);
      float v = sqrt((vx_r * vx_r) + (vy_r * vy_r));

      x_ << px_r, py_r, v, 0.0, 0.0;
    }
    
    is_initialized_ = true;

    cout << "Initialization Complete!!!";

  }


  double delta_t = meas_package.timestamp_ - time_us_;
  UKF::Prediction(&delta_t);

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

  /*
  * This step involves the following
  * 1. Augmentation : Generating augmented sigma points & others : x, P_aug, Xsig_aug
  * 2. Predicting Sigma points from t=k(Xsig_aug) to t=k+1
  * 3. Calculating mean and covariance of this predicted sigma points
  */


  /*****************************STEP 1 : GENERATING SIGMA POINTS****************************/

  //Initializations
  int mat_len = 2 * n_aug_ + 1;
  lambda_ = 3 - n_aug_;
  
  //Initializing Augmentation Matrices
  VectorXd x_aug    = VectorXd(7);
  MatrixXd P_aug    = MatrixXd(7,7);
  MatrixXd Xsig_aug = MatrixXd(7, 15);

  // AUGMENTING X (STAE VECTOR)
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // GENERATING PROCESS COVARIANCE MATRIX
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.topLeftCorner(n_x, n_x) = P;
  P_aug.bottomRightCorner(2, 2) << Q_;

  //Caclulating squareroot
  MatrixXd L = P_aug.llt().matrixL();
  
  // GENERATING SIGMA POINTS
  Xsig_aug.col(0) = x_;
  for(int i=0; i<mat_len; i++) {
    Xsig_aug.col(i+1)             = x + sqrt(lambda + n_aug_) * A.col(i);
    Xsig_aug.col(n_aug_ + i + +1) = x - sqrt(lambda + n_aug_) * A.col(i);
  }


  /*****************************STEP 2 : PREDICTING SIGMA POINTS****************************/
  /*
   * This step involves predicting sigma points from Xsig_aug at t=k to t=k+1
   * This involves passing each of the sima points through the process function
   * Note that the output matrix of this step will be Xsig_pred_ = [5 * 15] but not [7 * 15]
   * We consider noise while generating the sigma points, and the process function will 
   * use those to generate the predicted sigma points which are of the same dimension as the state vector
   * Noise parameters won't be considered after this step
   * the resultant predicted sigma points will later be used while transforming into measurement space aswell
   */
  for (int i = 0; i< 2*n_aug+1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //yawd = 0 ==> vechicle is driving in a straight line
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  cout << "Sigma Points Prediction Complete!!!";


  /*****************************STEP 3 : CALCULATING MEAN & COVARIANCE OF PREDICTING SIGMA POINTS****************************/
  /*
   * This step involves calculating the mean and covariance
   * This involves calculating the weights first
   */

  weights_ = VectorXd(2*n_x_+1);
  for(int i=0; i<2*n_aug_+1; i++) {
      double ind_weight = 0.0;
      if(i==0) {
          ind_weight = lambda / (lambda + n_aug_);
      } else {
          ind_weight = 1 / (2 * (lambda + n_aug_));
      }
      
      weights(i) = ind_weight;
      
  }

  // Calculating Mean
  VectorXd Xpred_mean = VectorXd(n_x_);
  MatrixXd Ppred = MatrixXd(n_x_, n_x_);
  for(int i=0; i<2*n_aug+1; i++) {
    Xpred_mean += weights(i) * Xsig_pred_.col(i);
  }
  
  //Calculating predicted state covariance matrix
  for(int i=0; i<2*n_aug+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    
    //Normalizing yaw angle(psi)
    while (x_diff(3)> M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2*M_PI;
    
    Ppred += weights(i) * x_diff * x_diff.transpose();
  }

  if(meas_package.sensor_type_ == MeasurementPackage::'LASER') {
    UKF::UpdateLidar(&meas_package);
  } else if(meas_package.sensor_type_ == MeasurementPackage::'RADAR') {
    UKF::UpdateRadar(&meas_package);
  }




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
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
