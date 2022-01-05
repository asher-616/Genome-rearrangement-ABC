#include "ReABC_params.h"

Re_params::Re_params(string inputTreeFileName)
{
	_inputTree = tree(inputTreeFileName);
}

void Re_params::init_params(double inv_rate, double trans_rate, double a_param) {
	// M0 parameters initiation
	_inv_rate = inv_rate;
	_trans_rate = trans_rate;
	_a_param = a_param;

	// M1 parameters are initialized to 0 for m0 runs
	_fis_rate = 0.0;
	_fus_rate = 0.0;

	// M2 parameters are initialized to 0 for lower models
	_dup_rate = 0.0;
	_loss_rate = 0.0;
	_root_a_param = 0.0;
}

void Re_params::init_params(double inv_rate, double trans_rate, double a_param, double fis_rate, double fus_rate) {
	// M0 parameters initiation
	_inv_rate = inv_rate;
	_trans_rate = trans_rate;
	_a_param = a_param;

	// M1 parameters
	_fis_rate = fis_rate;
	_fus_rate = fus_rate;

	// M2 parameters are initialized to 0 for lower models
	_dup_rate = 0.0;
	_loss_rate = 0.0;
	_root_a_param = 0.0;
}
Re_params::~Re_params()
{
}
