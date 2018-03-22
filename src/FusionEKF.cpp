#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);
  
  //state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  //the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

  //the initial measurement function for laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

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
    previous_timestamp_ = measurement_pack.timestamp_; //Updating timestamp to previous_timestamp var
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_[0];
      float phie = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      float pos_x = ro * cos(phie);
      //check division by zero
      if(pos_x<0.00001){pos_x=0.00001;}

      float pos_y = ro * sin(phie);
      
      //check division by zero
      if(pos_y<0.00001){pos_y=0.00001;}
      
      float vel_x = ro_dot * cos(phie);

      float vel_y = ro_dot * sin(phie);

      ekf_.x_ << pos_x , pos_y , vel_x , vel_y;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_<<measurement_pack.raw_measurements_[0] , measurement_pack.raw_measurements_[1], 0 , 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //Calculating dt which is the time difference between measurements
  float delt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  //Updating the previous_timestamp to reflect the current timestamp
  previous_timestamp_ =  measurement_pack.timestamp_;

  //Define some useful/repetetive powers for delta t
  float delt_2 = delt * delt;
  float delt_3 = delt_2 * delt;
  float delt_4 = delt_3 * delt;

  //As suggested in skeleton comments set the noise values accordingly
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  //Modify the state transition matrix  as per the updated timestamp
  ekf_.F_(0,2) = delt;
  ekf_.F_(1,3) = delt;

  //set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ <<  delt_4/4*noise_ax, 0, delt_3/2*noise_ax, 0,
			   0, delt_4/4*noise_ay, 0, delt_3/2*noise_ay,
			   delt_3/2*noise_ax, 0, delt_2*noise_ax, 0,
			   0, delt_3/2*noise_ay, 0, delt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tool; //Create Tool object
    //Obtain the Jacobian by calling the calculateJacobian from tools class
    Hj_ = tool.CalculateJacobian(ekf_.x_);
    //Set it to the measurement function
    ekf_.H_ = Hj_; 
    //Set the measurement noise to the one for Radar
    ekf_.R_ = R_radar_;
    //Call the update function
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    //Set the Measurement function
    ekf_.H_ = H_laser_;
    //Set the measurement noise to the one for Lidar
    ekf_.R_ = R_laser_;
    //Call the update function
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
