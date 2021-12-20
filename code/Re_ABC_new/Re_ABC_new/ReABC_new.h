#pragma once

#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <vector>


#include "ReABC_params.h"



vector<vector<double> > perform_simulation(vector<int> ChromosomeLengths, double a_param, double inv_rate, double trans_rate, 
	double fus_rate, double fis_rate, double dup_rate, double loss_rate, double root_a_param);
	//Note that double dup_rate, double loss_rate, double root_a_param are not yet in used. to be used in M2. A.M 20.12.21