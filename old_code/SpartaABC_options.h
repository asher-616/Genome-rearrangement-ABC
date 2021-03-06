#ifndef _SPARTAABC_OPTIONS
#define _SPARTAABC_OPTIONS

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cfloat>
#include "errorMsg.h"
#include "Parameters.h"

using namespace std;

class SpartaABC_options {
		
public:
    static void initOptions(const string& paramFileName);
	virtual ~SpartaABC_options();	
	static void initDefault();
	static void getParamsFromFile(const string& paramFileName);
	static void setRLPrior(int msa_min_len, int msa_max_len);

public:
	//summary statistics weigths (0 = not used, 1 for all = unweighted, any positive double = given weight, -1 = estimate [not implmented yet])
	static double _wAvgGapSize; 
	static double _wAvgUniqueGapSize;
	static double _wMSALen;
	static double _wMSAMax;
	static double _wMSAMin; 
	static double _wTotNumGaps; // total number of gaps
	static double _wNumGapsLenOne; // number of gaps of length one
	static double _wNumGapsLenTwo; // number of gaps of length two
	static double _wNumGapsLenThree; // number of gaps of length three
	static double _wNumGapsLenAtLeastFour;
	static double _wTotNumUniqueGaps; // total number of unique gaps
	static double _wNumberOfMSA_position_with_0_gaps;
	static double _wNumberOfMSA_position_with_1_gaps;
	static double _wNumberOfMSA_position_with_2_gaps;
	static double _wNumberOfMSA_position_with_n_minus_1_gaps;

	//template control
	static string _indelibleTemplateControlFile;
	static string _dawgTemplateControlFile;

	//input real MSA
	static string _inputRealMSAFile;

	//input tree file name
	static string _inputTreeFileName;

	//output collected samples
	static string _outputGoodParamsFile;

	//parameters prior options
	static string _priorDistTypeRL;
	static string _priorDistTypeA;
	static string _priorDistTypeIR;

	static int _minRLVal;
	static int _maxRLVal;

	static double _minAVal;
	static double _maxAVal;

	static double _minIRVal;
	static double _maxIRVal;

	

	//ABC Reject params
	static int _numberOfSamplesToKeep;
	static double _distanceCutOff;

	//pairwise or multiple
	static int _alignmentMode;

	//similarity - simple / nucleotide / aa
	static int _similarityMode;

	// should dawg simulator be used (defaults to 'false')
	static bool _dawgSimulator;

};


#endif