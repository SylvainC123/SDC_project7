#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <cmath>
#include <iostream>
#define minVal 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
 
UKF::UKF() {
  /// initially set to false, set to true in first call of ProcessMeasurement	
  is_initialized_ = false;
  
  /// if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  /// if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  /// initial state vector
  /// note: x_ is [px, py, vel, ang, ang_rate]
  x_ = VectorXd(5);

  /// initial covariance matrix
  P_ = MatrixXd(5, 5);

  /// Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  /// Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;

  /// Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  /// Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  /// Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  /// Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  /// Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /// State dimension
  n_x_=x_.size();
  
  /// Augmented state dimension
  n_aug_= n_x_ +2;
    
  /// Set the predicted sigma points matrix dimentions
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ +1);
  
  /// Sigma point spreading parameter
  lambda_= 3 - n_aug_;
  
  /// Laser covariance matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << pow(std_laspx_,2), 0,
              0, pow(std_laspy_,2);
  
  /// Radar covariance 
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << pow(std_radr_,2), 0, 0,
              0, pow(std_radphi_,2), 0,
              0, 0, pow(std_radrd_,2);
              
  /// Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  /// init P_ covariance matrix
  P_<< 1,0,0,0,0,
		0,1,0,0,0,
		0,0,1,0,0,
		0,0,0,1,0,
		0,0,0,0,1;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
 
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  ///  Initialization
  if (!is_initialized_) {
	//cout<<P_<<endl;
	/// init timestamp for dt calculation
    time_us_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      
      /// RADAR init
      float ro     = meas_package.raw_measurements_[0];
      float phi    = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[2];
      /// from polar to cartesian
      float x =ro * cos(phi);
	  float y =ro * sin(phi);
	  float vx =ro_dot * cos(phi);
	  float vy =ro_dot * sin(phi);
	  float v  = sqrt(vx * vx + vy * vy);
      x_ << x, y, v, 0.0, 0.0;
    }        
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /// LIDAR init
      float x = meas_package.raw_measurements_[0];
      float y = meas_package.raw_measurements_[1];
      /// to avoid zero radius and div per zero error
      if (fabs(x) < minVal and fabs(y) < minVal){
		x = minVal;
		y = minVal;
	    }
      x_ << x,y,0.0,0.0,0.0;
    } 
    ///set weights
    for (int i=1;i<weights_.size();i++)
	  {
	  weights_(i)=0.5/(lambda_+n_aug_);
	}
    weights_(0)=lambda_/(lambda_+n_aug_);
    /// init done
    is_initialized_ = true; 
    //cout << "initialisation performed" <<endl;
    //cout<< "x_:" << x_ << endl;
    return;
  }
  /// Measurement call section
  /// timestep update
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  /// Radar case
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
	  //cout << "Radar " << meas_package.raw_measurements_[0] << " " << meas_package.raw_measurements_[1] << endl;
      UpdateRadar(meas_package);
  }

  /// Lidar case
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
	  //cout << "Lidar " << meas_package.raw_measurements_[0] << " " << meas_package.raw_measurements_[1] << endl;
      UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /// Generate sigma points
  
  ///	create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);  
  ///	create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  ///	create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  /// Predict sigma points
  ///	fill augmented mean state
  x_aug.fill(0.0);
  //  for (int i=0; n_x_;i++){
  //	  	 x_aug(i)=x_(i);
  //		}
  x_aug.head(n_x_) = x_;
  ///	fill augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(x_.size(),x_.size()) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  //cout << std_yawdd_<< endl;
  //cout << P_aug << endl;
  ///	Square root of P matrix
  MatrixXd L = P_aug.llt().matrixL();
  ///create augmented sigma points
  ///	fill central point
  Xsig_aug.col(0)= x_aug;
  ///	fill pair of points
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  //cout << Xsig_aug << endl;
  ///predict sigma points
  
  ///	predicted state values
  double px_p,py_p;

  for (int j=0; j< 2 * n_aug_ + 1;j++)
    {
    double p_x = Xsig_aug(0,j);
    double p_y = Xsig_aug(1,j);
    double v = Xsig_aug(2,j);
    double yaw = Xsig_aug(3,j);
    double yawd = Xsig_aug(4,j);
    double nu_a = Xsig_aug(5,j);
    double nu_yawdd = Xsig_aug(6,j);
 
    ///	  avoid division by zero
    if (fabs(yawd) > minVal) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
        }
        
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    double v_p = v;
  
    ///	  add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
  
    ///	  write predicted sigma points into right column
    Xsig_pred_(0,j) = px_p;
    Xsig_pred_(1,j) = py_p;
    Xsig_pred_(2,j) = v_p;
    Xsig_pred_(3,j) = yaw_p;
    Xsig_pred_(4,j) = yawd_p;
    }

  ///predicted state mean 
  x_ = Xsig_pred_ * weights_;
 
  /// Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    /// iterate over sigma points
    /// State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    /// Angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
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
  /// Set measurement dimension, radar can measure px, py
  int n_z = 2;
  /// Create matrix for sigma points in measurement space
  /// Transform sigma points into measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);

  ///mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred  = Zsig * weights_;

  ///measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    ///residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    ///angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  ///add measurement noise covariance matrix
  S = S + R_lidar_;

  /// Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  ///calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    ///residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    ///angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    /// state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    ///angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
  /// Measurements
  VectorXd z = meas_package.raw_measurements_;
  ///Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  ///residual
  VectorXd z_diff = z - z_pred;

  ///update state mean and covariance matrix
  x_= x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  /// Calculate NIS
  NIS_lidar_ = z.transpose() * S.inverse() * z;
   
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  ///set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  ///create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);  
    
  ///transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    /// extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    /// measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  ///mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
      }

  ///measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    ///residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    ///angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  ///add measurement noise covariance matrix
  S = S + R_radar_;

  /// Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  ///calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    ///residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    ///angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    /// state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    ///angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
  /// Measurements
  VectorXd z = meas_package.raw_measurements_;
  ///Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  ///residual
  VectorXd z_diff = z - z_pred;

  ///angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  ///update state mean and covariance matrix
  x_= x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  /// Calculate NIS
  NIS_radar_ = z.transpose() * S.inverse() * z;
 
}
