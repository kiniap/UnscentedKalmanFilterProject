#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// Ensure that the estimation vector size is non-zero
	if (estimations.size() == 0){
		cout << "Error! Size of estimations is zero!" << endl;
		return rmse;
	}

	// Ensure that the size of the estimations vector is the same as the size of the ground-truth vector
	if (estimations.size() != ground_truth.size()){
		cout << "Error! Size of estimations does not match size of ground truth!" << endl;
		return rmse;
	}

	// Loop through the estimations and accumulate residuals
	for(int i=0;i<estimations.size();++i){
		VectorXd residual = ground_truth[i] - estimations[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse/estimations.size();

	// Calculate the square root
	rmse = rmse.array().sqrt();
	cout << "RMSE = " << rmse << endl;
	return rmse;
}
