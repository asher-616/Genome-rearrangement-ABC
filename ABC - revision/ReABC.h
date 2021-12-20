#pragma once
#ifndef _GRABC
#define _GRABC
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>

#include "GRABC_options.h"
//#include "summary_stats_defs.h"
#include "euclidean_distance.h"
//#include "param_priors.h"
#include "summ_stats_wrapper.h"
#include "genomes.h"
#include "simulator.h"



using namespace std;

vector<vector<double>> simulateRealData();
vector<vector<double> > getStatVec(genomes & currMSA);
//vector<double> getWeightsVec();
vector<string> readIndelibleTemplateControlFile(string indelibleTemplateCtrlFile);
vector<double> tranformSumStatVector(vector<vector<double> >&, int);
vector<vector<double>> tranformSumStatVector2(vector<vector<double> >&, int);
vector<int> generate_chromosome_sizes(int chromosome_number, int genome_size);
vector<vector<double> > simulateSequencesAndReturnSummaryStatistics(vector<int> randChromosomeLength, double randAParam, double randInvertRatio,
	double randTranslocateRatio, double randFusionRate, double randFissionRate, double randDuplicationRate, double randLossRate, double randRootAparam);
#endif