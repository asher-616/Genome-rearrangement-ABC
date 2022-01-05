#pragma once
// new parameter object. 20.12.21
// The new C++ code will perform a single simulation 
// therefore the old parameter object (GRABC_options) is irrelevant

#include <string>
#include "tree.h"
// include stuff here as needed

using namespace std;

class Re_params
{
public:
	Re_params(string inputTreeFileName);
	~Re_params();

	void init_params(double inv_rate, double trans_rate, double a_param);
	void init_params(double inv_rate, double trans_rate, double a_param, double fis_rate, double fus_rate);

	// M0 parameters
	static tree _inputTree;
	static double _inv_rate;
	static double _trans_rate;
	static double _a_param;
	
	// M1 parameters
	static double _fis_rate;
	static double _fus_rate;

	// M2 parameters (for future versions)
	static double _dup_rate;
	static double _loss_rate;
	static double _root_a_param;
private:

};

