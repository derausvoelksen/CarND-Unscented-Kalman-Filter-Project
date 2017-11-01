#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
  * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse.fill(0);
  double max_value = 1 << 30;

  if (estimations.size() != ground_truth.size()) {
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse /= estimations.size();
  //calculate the squared root
  rmse = rmse.array().sqrt();

  //Allow huge and nan rmse to go through
  for (int i = 0; i < rmse.size(); ++i) {
    if (rmse(i) != rmse(i)) {
  	  rmse(i) = max_value;
    }
    rmse(i) = min(rmse(i), max_value);
  }

  return rmse;
}