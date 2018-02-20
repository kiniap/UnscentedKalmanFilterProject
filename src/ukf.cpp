#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  // set is_initialized to false
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //set state dimension
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3-n_x_;
  n_z_=3; // initially set it to the number of radar measurements

  // Sigma prediction matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_+1);

  // previous timestamp
  previous_timestamp_ = 0.0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // **assuming max accel of a bicycle to be 3m/s^2 (2*std_a)
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // **assuming max yaw accel of slightly more than pi/4 rad/s^2, std_yawdd of greater than pi/8 seems to work pretty well
  std_yawdd_ = 0.5;
  
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

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix to identity matrix
  // **Identity seems to produce the best rmse and NIS
  P_ = MatrixXd::Identity(n_x_, n_x_);
  //state covariance matrix P
  // initialize the P matrix using the variances from the measurement noise
  //  P_ << std_laspx_*std_laspx_, 0, 0, 0, 0, // using std_laspx to initialize noise in px
  //		 0, std_laspy_*std_laspy_, 0, 0, 0, // using std_laspy to initialize noise in py
  //		 0, 0, 1, 0, 0,
  //		 0, 0, 0, 1, 0,
  //		 0, 0, 0, 0, 1;

  //laser measurement matrix
  H_laser_ = MatrixXd(2,5);
  H_laser_ << 1, 0, 0, 0, 0,
			  0, 1, 0, 0, 0;

  //add measurement noise covariance matrix
  R_laser_ = MatrixXd(2,2);
  R_laser_ <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

  //radar measurement noise covariance matrix
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
          	  0, std_radphi_*std_radphi_, 0,
			  0, 0,std_radrd_*std_radrd_;


  // initialize NIS values for radar and lidar
  epsilon_NIS_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if(!is_initialized_){

	 // first measurement
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	  /**
	  Convert radar from polar to cartesian coordinates and initialize state.
	  */
	  const double rho = meas_package.raw_measurements_(0);
	  const double phi = meas_package.raw_measurements_(1);
	  const double rho_dot = meas_package.raw_measurements_(2);

	  const double px = rho*cos(phi);
	  const double py = rho*sin(phi);

	  // Assuming rho_dot is also at angle psi initially
	  const double v = rho_dot;

	  // Assuming psi = phi initially
	  const double psi = phi;

	  // assuming psi_dot = 0 initially
	  // TODO: should I be assuming something else here?
	  const double psi_dot = 0;
	  x_ << px, py, v, psi, psi_dot;

	}
	else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {

		const double px = meas_package.raw_measurements_(0);
		const double py = meas_package.raw_measurements_(1);

		// Assuming yaw is the atan2(px,py)
		const double psi = atan2(px, py);

		// Assuming zero velocity and zero yaw rate
		x_ << px, py, 0, psi, 0;
	}

    // done initializing, no need to predict or update
	previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;
  }

  //compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  /*****************************************************************************
   *  Prediction Step
   ****************************************************************************/

  Prediction(dt);

  /*****************************************************************************
   *  Update step
   ****************************************************************************/
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
	  UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
	  UpdateLidar(meas_package);
  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
	// Generate sigma points
  //MatrixXd Xsig = MatrixXd::Zero(n_x_, 2*n_x_+1);
  //GenerateSigmaPoints(Xsig);

  // Augment the sigma points to account for process noise
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2*n_aug_+1);
  AugmentedSigmaPoints(Xsig_aug);

  // Sigma point predictions
  SigmaPointPrediction(Xsig_aug, delta_t);

  // Update weights
  VectorXd weights = VectorXd::Zero(2*n_aug_+1);
  UpdateWeights(weights);

  // Predict mean and covariance
  PredictMeanAndCovariance(weights);
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {

	// set measurement dimension. Laser can measure px and py
	n_z_ = 2;

	const VectorXd z = meas_package.raw_measurements_;
	const VectorXd z_pred = H_laser_*x_;
	// Measurement update
	const VectorXd y = z - z_pred;

	const MatrixXd PHt = P_ * H_laser_.transpose();
	const MatrixXd S = H_laser_ * PHt + R_laser_;
	const MatrixXd K = PHt * S.inverse();

	x_ += K * y;
	P_ -= K * H_laser_ * P_;

	// Compute NIS
	epsilon_NIS_ = NIS(meas_package.raw_measurements_, z_pred, S);
	std::cout << "Lidar NIS: " << epsilon_NIS_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidarUKF(MeasurementPackage meas_package)
{
  // set measurement dimension. Laser can measure px and py
  n_z_ = 2;

  // Create vector to store predicted Z values
  VectorXd z_pred = VectorXd::Zero(n_z_);
  // Created matrix to store predicted covariance matrix
  MatrixXd S_pred = MatrixXd::Zero(n_z_,n_z_);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //set vector for weights
  VectorXd weights = VectorXd::Zero(2*n_aug_+1);
  UpdateWeights(weights);

  // Predict lidar measurements
  PredictLidarMeasurement(Zsig, z_pred, S_pred, weights);

  //Update state vector
  UpdateState(Zsig, z_pred, meas_package, S_pred, weights);

  // Compute NIS
  epsilon_NIS_ = NIS(meas_package.raw_measurements_, z_pred, S_pred);
  std::cout << "Lidar NIS: " << epsilon_NIS_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;

  // Create vector to store predicted Z values
  VectorXd z_pred = VectorXd::Zero(n_z_);
  // Created matrix to store predicted covariance matrix
  MatrixXd S_pred = MatrixXd::Zero(n_z_,n_z_);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z_, 2 * n_aug_ + 1);

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //set vector for weights
  VectorXd weights = VectorXd::Zero(2*n_aug_+1);

  // Update weights
  UpdateWeights(weights);

  // Predict radar measurement
  PredictRadarMeasurement(Zsig, z_pred, S_pred, weights);

  // Predict state
  UpdateState(Zsig, z_pred, meas_package, S_pred, weights);

  // Compute NIS
  epsilon_NIS_ = NIS(meas_package.raw_measurements_, z_pred, S_pred);
  std::cout << "Radar NIS: " << epsilon_NIS_ << std::endl;
}

/**
 *  Generate Sigma Points
 */
void UKF::GenerateSigmaPoints(MatrixXd& Xsig)
{
  //define spreading parameter
  lambda_ = 3 - n_x_;

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig
  Xsig.col(0) = x_;
  Xsig.block<5,5>(0,1) = x_.replicate(1,5) + std::sqrt(lambda_+n_x_)*A;
  Xsig.block<5,5>(0,6) = x_.replicate(1,5) - std::sqrt(lambda_+n_x_)*A;
}

/**
 *  Augmented Sigma Points
 */
void UKF::AugmentedSigmaPoints(MatrixXd& Xsig_aug)
{
  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create augmented mean state
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_)=P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //calculate square root of P_aug
  MatrixXd B = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block<7,7>(0,1) = x_aug.replicate(1,7) + std::sqrt(lambda_+n_aug_)*B;
  Xsig_aug.block<7,7>(0,8) = x_aug.replicate(1,7) - std::sqrt(lambda_+n_aug_)*B;
}

/**
 *  Sigma point prediction
 */
void UKF::SigmaPointPrediction(const MatrixXd& Xsig_aug, double delta_t)
{
  //predict sigma points
  for(unsigned i=0; i< 2*n_aug_+1; ++i )
  {
	  const double px = Xsig_aug(0,i);
	  const double py = Xsig_aug(1,i);
	  const double v = Xsig_aug(2,i);
	  const double yaw = Xsig_aug(3,i);
	  const double yawd = Xsig_aug(4,i);
	  const double nu_a = Xsig_aug(5,i);
	  const double nu_yawdd = Xsig_aug(6,i);

	  //predicted state
	  double px_pred, py_pred, v_pred, yaw_pred, yawd_pred;

	  // avoid division by zero
	  if (fabs(yawd) > 0.001)
	  {
		  px_pred = px + (v/yawd)*(sin(yaw+yawd*delta_t)-sin(yaw));
		  py_pred = py + (v/yawd)*(-cos(yaw+yawd*delta_t)+cos(yaw));
	  }else
	  {
		  px_pred = px + v*cos(yaw)*delta_t;
		  py_pred = py + v*sin(yaw)*delta_t;
	  }

	  v_pred = v;
	  yaw_pred = yaw+yawd*delta_t;
	  yawd_pred = yawd;

	  // add noise
	  px_pred += 0.5*delta_t*delta_t*cos(yaw)*nu_a;
	  py_pred += 0.5*delta_t*delta_t*sin(yaw)*nu_a;
	  v_pred += delta_t*nu_a;
	  yaw_pred += 0.5*delta_t*delta_t*nu_yawdd;
	  yawd_pred += delta_t*nu_yawdd;

	  //write predicted sigma points into right column
	  Xsig_pred_(0,i) = px_pred;
	  Xsig_pred_(1,i) = py_pred;
	  Xsig_pred_(2,i) = v_pred;
	  Xsig_pred_(3,i) = yaw_pred;
	  Xsig_pred_(4,i) = yawd_pred;
	}
}

/**
*  Predict mean and covariance
*/
void UKF::PredictMeanAndCovariance(const VectorXd& weights) {
//  //create vector for predicted state
//  VectorXd x = VectorXd::Zero(n_x_);
//

//
//  //predict state mean
//  for (unsigned i=0;i<2*n_aug_+1;++i)
//  {
//	  x = x + weights(i)*Xsig_pred_.col(i);
//  }
//
  x_ = Xsig_pred_* weights;

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd::Zero(n_x_, n_x_);

  //predict state covariance matrix
  for (unsigned i=0;i<2*n_aug_+1;++i)
  {
	  P = P + weights(i)*(Xsig_pred_.col(i) - x_)*(Xsig_pred_.col(i) - x_).transpose();
  }

  P_ = P;
  //MatrixXd Xdiff = Xsig_pred_ - x_.replicate(1,15);
  //P_ = (Xdiff * Xdiff.transpose()) * weights.replicate(1,n_x_);

//
//  //write result
//  x_ = x;
}

/**
 *  Update weights
 */
void UKF::UpdateWeights(VectorXd& weights)
{
	weights(0) = lambda_/(lambda_+n_aug_);
	for (int i=1; i<2*n_aug_+1; i++)
	{
		weights(i) = 0.5/(n_aug_+lambda_);
	}
}

/*
 * Predict radar measurement
 */
void UKF::PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S_pred, const VectorXd& weights) {

  //transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    const double p_x = Xsig_pred_(0,i);
    const double p_y = Xsig_pred_(1,i);
    const double v  = Xsig_pred_(2,i);
    const double yaw = Xsig_pred_(3,i);

    const double v1 = cos(yaw)*v;
    const double v2 = sin(yaw)*v;
    const double eps = 0.001;

    // measurement model
    const double rho = sqrt(p_x*p_x + p_y*p_y);
    double phi = 0.0;
    if (p_x != 0 || p_y != 0) // both can't be zero at the same time
    	phi = atan2(p_y,p_x);

    const double rho_dot = (p_x*v1 + p_y*v2 ) / std::max(eps, rho);

    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_dot;
  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  S_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // normalize angles
    //z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
    NormalizeAngle(z_diff(1));

    S_pred = S_pred + weights(i) * z_diff * z_diff.transpose();
  }

  S_pred = S_pred + R_radar_;
}

/*
 * Predict lidar measurement
 */

void UKF::PredictLidarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S_pred, const VectorXd& weights)
{
  //transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  S_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S_pred = S_pred + weights(i) * z_diff * z_diff.transpose();
  }

  S_pred = S_pred + R_laser_;
}

/*
 * Update state
 */
void UKF::UpdateState(const MatrixXd& Zsig, const VectorXd& z_pred, const MeasurementPackage& meas_package, const MatrixXd& S_pred, const VectorXd& weights)
{
  // create the measurements vector
  VectorXd z = meas_package.raw_measurements_;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for(unsigned i=0;i < 2*n_aug_+1;++i)
  {
    // calculate residuals
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i)-z_pred;

    // normalize angles for radar type sensor
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    	//z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
    	NormalizeAngle(z_diff(1));
    }

    Tc = Tc + weights(i)*x_diff*z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_z_,n_z_);
  K = Tc*S_pred.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K*(z-z_pred);
  P_ = P_ - K*S_pred*K.transpose();
}

/*
 * Compute NIS to check consistency of the filter
*/
double UKF::NIS(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S_pred)
{
	double epsilon = (z - z_pred).transpose()*S_pred.inverse()*(z - z_pred);
	return epsilon;
}

void UKF::NormalizeAngle(double& phi)
{
	phi = atan2(sin(phi), cos(phi));
}
