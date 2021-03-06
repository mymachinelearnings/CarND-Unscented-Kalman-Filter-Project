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
  use_laser_ = true;
  use_radar_ = true;
  x_ = VectorXd(5);
  P_ = MatrixXd(5, 5);

  std_a_ = 1.5;
  std_yawdd_ = 0.5;

  //Measurement Noise - given by Manufacturer - DO NOT CHANGE
  std_laspx_ = 0.15;
  std_laspy_ = 0.15;
  std_radr_ = 0.3;
  std_radphi_ = 0.03;
  std_radrd_ = 0.3;

  n_x_ = 5;
  n_aug_ = 7;

  is_initialized_ = false;
  Xsig_pred_ = MatrixXd(5, 15);
  weights_ = VectorXd(15);

  //Augmented variables will be used within the function, hence they are declared locally
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_) {
    time_us_ = meas_package.timestamp_;
    double cov = 0.5;
    P_ << cov, 0, 0, 0, 0,
          0, cov, 0, 0, 0,
          0, 0, cov, 0, 0,
          0, 0, 0, cov, 0,
          0, 0, 0, 0, cov;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0.0,
            0.0,
            0.0;
    } else if(meas_package.raw_measurements_[0] == MeasurementPackage::RADAR) {
      float ro     = meas_package.raw_measurements_[0];
      float theta  = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[2];

      float px_r = ro * cos(theta);
      float py_r = ro * sin(theta);
      float vx_r = ro_dot * cos(theta);
      float vy_r = ro_dot * sin(theta);
      float v = sqrt((vx_r * vx_r) + (vy_r * vy_r));

      x_ << px_r,
            py_r,
            v,
            0.0,
            0.0;
    }

    is_initialized_ = true;
    cout << "Initialization Complete!!!\n";
    return;
  }


  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
    //cout << "Lidar Update complete!\n";

  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
    //cout << "Radar Update complete!\n";
}

}

/*
* This step involves the following
* 1. Augmentation : Generating augmented sigma points & others : x, P_aug, Xsig_aug
* 2. Predicting Sigma points from t=k(Xsig_aug) to t=k+1
* 3. Calculating mean and covariance of this predicted sigma points
*/

void UKF::Prediction(double delta_t) {


  /*****************************STEP 1 : GENERATING SIGMA POINTS****************************/

  //Initializations
  lambda_ = 3 - n_aug_;

  //Initializing Augmentation Matrices
  VectorXd x_aug    = VectorXd(7);
  MatrixXd P_aug    = MatrixXd(7,7);
  MatrixXd Xsig_aug = MatrixXd(7, 15);
  x_aug.fill(0.0);
  P_aug.fill(0.0);
  Xsig_aug.fill(0.0);


  // Augmenting x (state vector with extra noise elements initialized to 0.0 each)

  //Process noise covariance matrix
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_*std_a_, 0,
        0, std_yawdd_*std_yawdd_;

  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // Augmenting Process covariance matrix (P)
  P_aug.topLeftCorner(5, 5) = P_;
  //P_aug.bottomRightCorner(2, 2) = Q;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //cout << "Augmentation Complete \n----------------------\n" << x_aug << "\n----\n" << P_aug << "\n-----------------------\n";

  //Caclulating squareroot
  MatrixXd P_sqrt = P_aug.llt().matrixL();
  // GENERATING SIGMA POINTS
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_; i++) {
    Xsig_aug.col(i+1)            = x_aug + sqrt(lambda_ + n_aug_) * P_sqrt.col(i);
    Xsig_aug.col(n_aug_ + 1 + i) = x_aug - sqrt(lambda_ + n_aug_) * P_sqrt.col(i);
    //cout << "Augmented Sigma Points at step --------- " << i << "\n" << Xsig_aug << "\n ---------------- \n" ;
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
  for (int i = 0; i< 2*n_aug_+1; i++) {
    //extract values for better readability
    double p_x  = Xsig_aug(0,i);
    double p_y  = Xsig_aug(1,i);
    double v    = Xsig_aug(2,i);
    double yaw  = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //vechicle is driving in a straight line ==> yawd=0
    if (fabs(yawd) > 0.001) {
        px_p = p_x + (v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw)));
        py_p = p_y + (v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) ));
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    //v, yawd came from the original x_, dont confuse with 0 in the formula
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p   += 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p   += 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p    += nu_a*delta_t;
    yaw_p  += 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p += nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;


  }
  //cout << "Predicted Sigma Points\n--------------\n" << Xsig_pred_ << "\n ---------------- \n" ;




  /*****************************STEP 3 : CALCULATING MEAN & COVARIANCE OF PREDICTING SIGMA POINTS****************************/
  /*
   * This step involves calculating the mean and covariance
   * This involves calculating the weights first
   */

  //weights_ = VectorXd(2*n_x_+1); This is already defined
  //weights_.fill(0.0); //REMOVE THIS if no change
  for(int i=0; i<2*n_aug_+1; i++) {
      double ind_weight = 0.0;
      if(i==0) {
          ind_weight = lambda_ / (lambda_ + n_aug_);
      } else {
          ind_weight = 1 / (2 * (lambda_ + n_aug_));
      }

      weights_(i) = ind_weight;
  }

  //cout << "Weights \n" << weights_ << endl;

  // Calculating Mean
  //VectorXd Xpred_mean = VectorXd(n_x_);
  //MatrixXd Ppred = MatrixXd(n_x_, n_x_);

  //temp x_ & P_
  VectorXd x_temp = VectorXd(5);
  MatrixXd P_temp = MatrixXd(5, 5);
  x_temp.fill(0.0);
  P_temp.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++) {
	  x_temp += weights_(i) * Xsig_pred_.col(i);
  }
  //x_temp = weights_ * Xsig_pred_;

  //Calculating predicted state covariance matrix
  for(int i=0; i<2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_temp;

    //Normalizing yaw angle(psi)
    while (x_diff(3)> M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2*M_PI;

    P_temp += weights_(i) * x_diff * x_diff.transpose();
  }

  x_ = x_temp;
  P_ = P_temp;

  //cout << "Predicted Mean & Covar\n-----------------\n" << x_ << "\n------\n" << P_ << "\n-----------------\n";
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /*****************************STEP 1 : TRANSFORMING PREDITECTED MEAN & COVARIANCE AT TIME k+1 TO MEASUREMENT SPACE****************************/
  /*
   * This step involves transforming the mean x [5] vector & covariance P [5 * 15]
   * into measurement space with mean z [3] vector & covariance S [3 * 21] matrix
   */
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_mean_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  //Measurement noise covariance matrix
  MatrixXd R = MatrixXd(2,2);
  R << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  //Zsig is for eventually calculating the predicted covariance in measurement space S
 Zsig.fill(0.0);
 for(int i=0; i<2*n_aug_+1; i++) {
     Zsig(0, i) = Xsig_pred_(0, i);
     Zsig(1, i) = Xsig_pred_(1, i);
 }

 //cout << "After conversion into measurement space \n----------------\n" << Zsig << "\n------------\n";


//calculate mean predicted measurement
z_mean_pred.fill(0.0);
for(int i=0; i<2*n_aug_+1; i++) {
    z_mean_pred += weights_(i) * Zsig.col(i);
}

//calculate covariance matrix S
S.fill(0.0);
for(int i=0; i<2*n_aug_+1; i++) {
    VectorXd Z_diff = Zsig.col(i) - z_mean_pred;

    S += weights_(i) * Z_diff * Z_diff.transpose();
  }

  S = S + R;

  /*****************************STEP 2 : UPDATING THE BELIEF****************************/
  /*
   * This step involves
   *
   */
  //Getting actual measurement at time k+1
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1];

  //cout << "LIDAR Prediction vs Measurement\n---------------------\n" << z_mean_pred << "\n-----------\n" << z << "\n---------------------\n";

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_mean_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_mean_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

/**
* In measurement space, the vecor dimension is 3 as radar provides 3 dimensions only
*
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /*****************************STEP 1 : TRANSFORMING PREDITECTED MEAN & COVARIANCE AT TIME k+1 TO MEASUREMENT SPACE****************************/
  /*
   * This step involves transforming the mean x [5] vector & covariance P [5 * 15]
   * into measurement space with mean z [3] vector & covariance S [3 * 21] matrix
   */
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_mean_pred = VectorXd(n_z);


  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  //Measurement noise covariance matrix
  MatrixXd R = MatrixXd(3,3);
  R << std_radr_*std_radr_, 0, 						 0,
       0, 					std_radphi_*std_radphi_, 0,
       0, 					0, 						 std_radrd_*std_radrd_;

  //Zsig is for eventually calculating the predicted covariance in measurement space S
 double px, py, v, psi, vx, vy;
 for(int i=0; i<2*n_aug_+1; i++) {
     px = Xsig_pred_(0, i);
     py = Xsig_pred_(1, i);
     v = Xsig_pred_(2, i);
     psi = Xsig_pred_(3, i); //you don't require Xsig_pred_.col(4) as we don't use it for transformation into measurement space

     vx = v * cos(psi);
     vy = v * sin(psi);

     Zsig(0, i) = sqrt(px*px + py*py);
     Zsig(1, i) = atan2(py, px);
     Zsig(2, i) = (px*vx + py*vy) / sqrt(px*px + py*py);
 }
 //cout << "After conversion into measurement space \n----------------\n" << Zsig << "\n------------\n";

//calculate mean predicted measurement
 z_mean_pred.fill(0.0);
 for(int i=0; i<2*n_aug_+1; i++) {
    z_mean_pred += weights_(i) * Zsig.col(i);
 }

//calculate covariance matrix S
S.fill(0.0);
for(int i=0; i<2*n_aug_+1; i++) {
    VectorXd Z_diff = Zsig.col(i) - z_mean_pred;

    while(Z_diff(1) > M_PI) Z_diff(1) -= 2*M_PI;
    while(Z_diff(1) < -M_PI) Z_diff(1) += 2*M_PI;


    S += weights_(i) * Z_diff * Z_diff.transpose();
  }

  S = S + R;

  /*****************************STEP 2 : UPDATING THE BELIEF****************************/
  /*
   * For that, first calculate the cross correlation between the predicted state and the measurement space
   * Use this to calculate K & uddate belief
   */
  //Getting actual measurement at time k+1
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
       meas_package.raw_measurements_[1],
       meas_package.raw_measurements_[2];
  //cout << "RADAR Prediction vs Measurement\n---------------------\n" << z_mean_pred << "\n-----------\n" << z << "\n---------------------\n";

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_mean_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_mean_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

}
