#pragma once
// new parameter object. 20.12.21
// The new C++ code will perform a single simulation 
// therefore the old parameter object (GRABC_options) is irrelevant

#include <string>
// include stuff here as needed

using namespace std;

class Re_params
{
public:
	Re_params();
	~Re_params();

	// M0 parameters
	static string _inputTreeFileName;
	static double inv_rate;
	static double trans_rate;
	static double a_param;
	
	// M1 parameters
	static double fis_rate;
	static double fus_rate;

private:

};

Re_params::Re_params()
{
}

Re_params::~Re_params()
{
}
