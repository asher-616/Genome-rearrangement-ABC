#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cfloat>
#include "errorMsg.h"
#include "Parameters.h"

using namespace std;

class GRABC_options {

public:
	static void initOptions(const string& paramFileName);
	virtual ~GRABC_options();
	static void initDefault();
	static void getParamsFromFile(const string& paramFileName);
	//static void setRLPrior(int msa_min_len, int msa_max_len);

public:
	//input real genomes
	static string _inputGenomesFile;

	//input tree file name
	static string _inputTreeFileName;

	//output collected samples
	static string _outputParamsFile;

	//parameters prior options
	static string _priorDistTypeRL;
	static string _priorDistTypeA;
	static string _priorDistTypeIR;
	static string _priorDistTypeFuR; //fusion rate distribution type
	static string _priorDistTypeFiR; //fission rate distribution type
	
	static string _priorDistTypeDR; //gene duplication rate distribution type for the future
	static string _priorDistTypeLR; //gene loss rate distribution type for the future


	static int _minGeneNumVal; //range of genes per chromosome
	static int _maxGeneNumVal;
	static int _minChNum; //range of chromosomes number
	static int _maxChNum;

	static double _minAVal; // range of a param for zipfian
	static double _maxAVal;

	static double _minInvVal; //inversion rate range
	static double _maxInvVal;

	static double _minTransVal;//tranposition rate range
	static double _maxTransVal;

	static double _minFissionVal;//fission rate range
	static double _maxFissionVal;

	static double _minFusionVal;//fusion rate range
	static double _maxFusionVal;

	static double _minDuplicationVal;//gene duplication rate range for the future
	static double _maxDuplicationVal;

	static double _minLossVal;//gene loss rate range for the future
	static double _maxLossVal;

	static double _minRootA; //used for gene family size in the root
	static double _maxRootA;
	static int _rootMaxFamilySize;

	static int _maxBlock; //should insert to param file used for summStat vector 

	
	//ABC Reject params
	static double _distanceCutOff; //this parameter isn't used in SpartaABC and was replaced bt a local variable
	static int _numberOfSamplesToKeep; //A.M 8.4 added for now because it is used in main and can't figure how to remove without messing algorithm

	// Testing fusion fission rates.
	// here we add counters for fusion and fission events to see if those are predicted correctly
	static int _fus_counter;
	static int _fis_counter;

	// for revision. counting inversions and translocations
	static int _inv_counter;
	static int _trans_counter;
};


