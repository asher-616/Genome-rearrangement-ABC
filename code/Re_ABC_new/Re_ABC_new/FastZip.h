#pragma once
#include <vector>
#include <iostream>
//#include "param_priors.h"
#include "RandomGenerators.h"
#include <cmath>
#include <algorithm>

using namespace std;

class FastZip
{
public:
	FastZip(double aParam, int max);
	~FastZip();
	int drawZip();

private:
	int numBins;
	vector<int> lowVector;
	vector<int> highVector;
	vector<double> lowHighRateVector;
};
