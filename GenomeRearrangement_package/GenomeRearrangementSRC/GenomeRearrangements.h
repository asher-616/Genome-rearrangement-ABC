#pragma once
// all imports commented. I want to make sure only the needed ones remain
//#include "someUtil.h"
#include <list> // used to perform events. check if can replaced with vectors
//#include <cstdlib>
#include "RandomGenerators.h"
//#include "param_priors.h"
//#include "GRABC_options.h" // should be changed to #include "ReABC_params.h"
#include "FastZip.h"
#include <vector>
//#include <algorithm>
#include <iostream>

using namespace std;

// 21.12.21 to count events on branch. 31.1.22 moved here from simulator
extern int inv_counter, trans_counter, fis_counter, fus_counter;

void drawRandomLocationOld(vector<int>& location, const vector<vector<int> > & genome, const int genomelength);
// drawRandomLocationOld was used pre-Itsik. Code left for reference
bool drawRandomLocation(vector<int>& location, const vector<vector<int> > & genome, const int genomelength);
void SimulateEvent(vector<vector<int> > & genome, int genomeLength, double transpositionRate, double inversionRate, FastZip blockLengthparameter,
	double fusionRate, double fissionRate, double duplicationRate, double deletionRate);
void Inversion(vector<vector<int> > & genome, vector<int> & eventLocation, int eventSize);
void Transposition(vector<vector<int> > & genome, vector<int> eventStartPoint, int eventSize, vector<int> eventdestination);
void Fusion(vector<vector<int> > & genome, int chromosome1, int chromosome2);
void Fission(vector<vector<int> > & genome, int genomeLength);
//to be added in the future
void GeneDuplication(vector<vector<int> > & genome, vector<int> eventStartPoint, vector<int> eventdestination);
void GeneLoss(vector<vector<int> > & genome, vector<int> eventStartPoint);

