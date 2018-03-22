#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd F_trans = F_.transpose();
  P_ = F_ * P_ * F_trans + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd H_trans = H_.transpose();
  
  VectorXd z_p = H_ * x_;

  VectorXd y = z - z_p;

  MatrixXd S = H_ * P_ * H_trans + R_;

  MatrixXd S_inv = S.inverse();

  MatrixXd P_Ht_Prod = P_ * H_trans;

  MatrixXd K = P_Ht_Prod * S_inv;

  //Update X

  x_ = x_ +(K * y);
  int size = x_.size();
  MatrixXd Iden = MatrixXd::Identity(size, size);
  P_ = (Iden - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float ro = sqrt(px*px +py*py);
  
  //check division by zero
        if(ro < 0.00001){
                ro = 0.00001;
        }

  float theta = atan2(py, px);
  float ro_dot = (px * vx + py * vy) / ro;

  VectorXd hp = VectorXd(3);
  hp << ro , theta, ro_dot;

  VectorXd y = z - hp;

  //Normalizing Angles

  if( y[1] > M_PI ){ y[1] -= 2*M_PI;}
  else if( y[1] < -M_PI ){ y[1] += 2*M_PI;}

  //Update the state using our steps as discussed in the course
  MatrixXd H_trans = H_.transpose();

  MatrixXd S = H_ * P_ * H_trans + R_;

  MatrixXd S_inv = S.inverse();

  MatrixXd K = P_ * H_trans *S_inv;

  //Update X

  x_ = x_ +(K * y);
  int size = x_.size();
  MatrixXd Iden = MatrixXd::Identity(size, size);
  P_ = (Iden - K * H_) * P_;
}
