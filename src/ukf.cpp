#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include "Eigen/src/StlSupport/StdVector.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {


	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;
	//set state dimension
	n_x_ = 5;
	
	// initial state vector
	//Changed this to match the practice code
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 2;// 0.2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = .55;// 0.2;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;// 0.0175;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = .3;// 0.1;

	/**
	  TEST: Complete the initialization. See ukf.h for other member properties.

	  Hint: one or more values initialized above might be wildly off...
	  */
	
	// time when the state is true, in us
	time_us_ = 0;

	//set augmented dimension
	n_aug_ = 7;

	// Weights of sigma points
	weights_ = VectorXd(2 * n_aug_ + 1);

	//define spreading parameter
	lambda_ = 3 - n_aug_;
	// the current NIS for radar
	NIS_radar_ = 0.0;

	// the current NIS for laser
	NIS_laser_ = 0.0;

	// predicted sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	is_initialized_ = false;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  DONE:  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
		(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_))
	{
		if (!is_initialized_)
		{
			//Initialize Matrix x
			x_ << 1, 1, 1, 1, .1;

			//Initialize the covariance matrix
			//P_ << 0, 0, 0, 0, 0,
			//	0, 0, 0, 0, 0,
			//	0, 0, 0, 0, 0,
			//	0, 0, 0, 0, 0,
			//	0, 0, 0, 0, 0;
			P_ <<
				0.0054342, -0.002405, 0.0034157, -0.0034819, -0.00299378,
				-0.002405, 0.01084, 0.001492, 0.0098018, 0.00791091,
				0.0034157, 0.001492, 0.0058012, 0.00077863, 0.000792973,
				-0.0034819, 0.0098018, 0.00077863, 0.011923, 0.0112491,
				-0.0029937, 0.0079109, 0.00079297, 0.011249, 0.0126972;
			

			//set the time to the TimeStamp in Measurement package
			time_us_ = meas_package.timestamp_;

			//Determine if it is radar or lidar
			if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
			{
				/**
				Convert radar from polar to cartesian coordinates and initialize state.
				*/
				float ro = meas_package.raw_measurements_(0);
				float phi = meas_package.raw_measurements_(1);
				float ro_dot = meas_package.raw_measurements_(2);
				x_(0) = ro * cos(phi);
				x_(1) = ro * sin(phi);
				x_(2) = ro_dot * cos(phi);
				x_(3) = ro_dot * sin(phi);
			}
			else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
			{
				/**
				Initialize state.
				*/
				x_(0) = meas_package.raw_measurements_(0);
				x_(1) = meas_package.raw_measurements_(1);
			}

			//variables initialized
			is_initialized_ = true;
			//Return as there is nothing to predict or update. 
			return;
		}

		//Once everything is initialized, then move on to prediction.
		//get the time difference between the current and the last measurement
		float timeDifference = (meas_package.timestamp_ - time_us_) / 1000000.0;
		time_us_ = meas_package.timestamp_;

		Prediction(timeDifference);


		//Once prediction is comeplete, update the proper sensor's data
		if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
		{
			UpdateLidar(meas_package);
		}
		else
		{
			UpdateRadar(meas_package);
		}
	}


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TEST: Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

	//DONE check this, seems like we are generating the sigma points twice here...
	////Generate the Sigma points.
	//MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
	//MatrixXd A = P_.llt().matrixL();
	//double lambda1 = 3 - n_x_;
	//Xsig.col(0) = x_;
	//for(int i = 0; i < n_x_; i++)
	//{
	//	Xsig.col(i + 1) = x_ + sqrt(lambda1 + n_x_) * A.col(i);
	//	Xsig.col(1 + 1 + n_x_) = x_ - sqrt(lambda1 + n_x_) * A.col(i);
	//}

	//Augment the Sigma points. 
	VectorXd x_aug = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//Create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	//create the square root matrix
	MatrixXd sqRtmatrix = P_aug.llt().matrixL();

	//Create the augmented Sigma points.
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i<n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * sqRtmatrix.col(i);
		Xsig_aug.col(i + 1+ n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * sqRtmatrix.col(i);
	}

	//Now perform the prediction for the augmented Sigma Points

	//predict Sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double posX = Xsig_aug(0, i);
		double posY = Xsig_aug(1, i);
		double velocity = Xsig_aug(2, i);
		double yawAngle = Xsig_aug(3, i);
		double yawRate = Xsig_aug(4, i);
		double longitAccel = Xsig_aug(5, i);
		double yawNoise = Xsig_aug(6, i);

		//these will hold the predicted values
		double predPosX, predPosY;

		//guard for divide by zero
		if(fabs(yawRate) > 0.001)
		{
			predPosX = posX + velocity / yawRate*(sin(yawAngle + yawRate * delta_t) - sin(yawAngle));
			predPosY = posY + velocity / yawRate*(cos(yawAngle) - cos(yawAngle + yawRate*delta_t));
		}
		else
		{
			predPosX = posX + velocity*delta_t*cos(yawAngle);
			predPosY = posY + velocity*delta_t*sin(yawAngle);
		}

		double predVelocity = velocity;
		double predYawAngle = yawAngle + yawRate*delta_t;
		double predYawRate = yawRate;

		//factor in the noise
		predPosX = predPosX + .5*longitAccel*delta_t*delta_t*cos(yawAngle);
		predPosY = predPosY + .5*longitAccel*delta_t*delta_t*sin(yawAngle);
		predVelocity = predVelocity + longitAccel*delta_t;
		predYawAngle = predYawAngle + .5*yawNoise*delta_t*delta_t;
		predYawRate = predYawRate + yawNoise*delta_t;

		//store the predicted values back in to the matrix
		Xsig_pred_(0, i) = predPosX;
		Xsig_pred_(1, i) = predPosY;
		Xsig_pred_(2, i) = predVelocity;
		Xsig_pred_(3, i) = predYawAngle;
		Xsig_pred_(4, i) = predYawRate;
	}

	//now we need to predict the mean and covariance
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for(int i = 1; i<2*n_aug_+1; i++)
	{
		double weight = .5 / (lambda_ + n_aug_);
		weights_(i) = weight;		
	}
	x_.fill(0.0);
	P_.fill(0.0);

	//predict the state mean
	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		x_ = x_ + weights_(i)*Xsig_pred_.col(i);
	}
	//predict state covariance matrix
	for (int i = 0; i<2 * n_aug_ + 1; i++)
	{
		VectorXd diffX = Xsig_pred_.col(i) - x_;

		//normalize the angles.
		while (diffX(3) > M_PI)diffX(3) -= 2.*M_PI;
		
		while (diffX(3) < -M_PI)diffX(3) += 2.*M_PI;		

		P_ = P_ + weights_(i) * diffX * diffX.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TEST: Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

	VectorXd z = meas_package.raw_measurements_;
	int n_z = 2;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S = MatrixXd(n_z, n_z);
	

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{  												
		double posX = Xsig_pred_(0, i);
		double posY = Xsig_pred_(1, i);

		Zsig(0, i) = posX;
		Zsig(1, i) = posY;
	}

	//calculate mean predicted measurement
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//calculate measurement covariance matrix S
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{ 												
		VectorXd z_diff = Zsig.col(i) - z_pred;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add the noise
	MatrixXd noise = MatrixXd(n_z, n_z);
	noise << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;
	S = S + noise;

	//Perform the state update
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	Tc.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{
		VectorXd diffZ = Zsig.col(i) - z_pred;
		//set the state difference
		VectorXd diffX = Xsig_pred_.col(i) - x_;

		Tc = Tc + weights_(i) * diffX * diffZ.transpose();
	}

	// calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd diffZ = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * diffZ;
	P_ = P_ - K*S*K.transpose();
	//calculate NIS
	NIS_laser_ = diffZ.transpose() * S.inverse() * diffZ;


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TEST: Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	VectorXd z = meas_package.raw_measurements_;
	int n_z = 3;
    //create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	//std::cout << "ZSig Initialized" << std::endl;
	//std::cout << Zsig << std::endl;
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S = MatrixXd(n_z, n_z);
    //transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double posX = Xsig_pred_(0, i);
		double posY = Xsig_pred_(1, i);
		double velocity = Xsig_pred_(2, i);
		double yawAngle = Xsig_pred_(3, i);

		double v1 = cos(yawAngle) * velocity;
		double v2 = sin(yawAngle) * velocity;

		//measurement model
		Zsig(0, i) = sqrt(posX*posX + posY*posY);
		Zsig(1, i) = atan2(posY, posX);
		//std::cout<<sqrt(posX*posX + posY*posY)<<std::endl;
		Zsig(2, i) = (posX*v1 + posY*v2) / sqrt(posX*posX + posY*posY);
		//std::cout << i << std::endl;
		//std::cout << "ZSig After" << std::endl;
		//std::cout << Zsig << std::endl;
		//calculate mean predicted measurement
		
		//std::cout << "weights_" << std::endl;
		//std::cout << weights_ << std::endl;
	}
	z_pred.fill(0.0);
	//mean predicted measurement
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		auto weight1 = weights_(i);
		auto ZsigCol = Zsig.col(i);
		z_pred = z_pred + weights_(i) * Zsig.col(i);
		//std::cout << "z_pred After" + i << std::endl;
		//std::cout << z_pred << std::endl;
	}
		
	//calculate measurement covariance matrix S
	S.fill(0.0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		
		VectorXd diffZ = Zsig.col(i) - z_pred;
		//Always normalize the angles
		//std::cout << diffZ << std::endl;
		while (diffZ(1) > M_PI)diffZ(1) -= 2 * M_PI;
		
		while (diffZ(1) <-M_PI)diffZ(1) += 2 * M_PI;

		S = S + weights_(i) * diffZ * diffZ.transpose();
	}

	//add the noise
	MatrixXd noise = MatrixXd(n_z, n_z);
	noise << std_radr_*std_radr_, 0, 0, 
			0, std_radphi_*std_radphi_, 0,
			0, 0, std_radrd_*std_radrd_;
	S = S + noise;

	//Perform the state update
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	Tc.fill(0.0);

	for(int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		VectorXd diffZ = Zsig.col(i) - z_pred;
		//always do angle normalization
		while (diffZ(1)> M_PI) diffZ(1) -= 2.*M_PI;
		while (diffZ(1)<-M_PI) diffZ(1) += 2.*M_PI;

		//set the state difference
		VectorXd diffX = Xsig_pred_.col(i) - x_;
		//always do angle normalization
		while (diffX(3)> M_PI) diffX(3) -= 2.*M_PI;
		while (diffX(3)<-M_PI) diffX(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * diffX * diffZ.transpose();
	}
	// calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();
	VectorXd diffZ = z - z_pred;
	while (diffZ(1)> M_PI) diffZ(1) -= 2.*M_PI;
	while (diffZ(1)<-M_PI) diffZ(1) += 2.*M_PI;
	//update state mean and covariance matrix
	x_ = x_ + K * diffZ;
	P_ = P_ - K * S * K.transpose();

	//Calculate the NIS
	NIS_radar_ = diffZ.transpose() * S.inverse() * diffZ;
	//}
}
