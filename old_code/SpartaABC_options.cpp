#include "SpartaABC_options.h"


double SpartaABC_options::_wAvgGapSize;
double SpartaABC_options::_wAvgUniqueGapSize;
double SpartaABC_options::_wMSALen; 
double SpartaABC_options::_wMSAMax; 
double SpartaABC_options::_wMSAMin; 
double SpartaABC_options::_wTotNumGaps;
double SpartaABC_options::_wTotNumUniqueGaps;
double SpartaABC_options::_wNumGapsLenOne;
double SpartaABC_options::_wNumGapsLenTwo;
double SpartaABC_options::_wNumGapsLenThree;
double SpartaABC_options::_wNumGapsLenAtLeastFour;

double SpartaABC_options::_wNumberOfMSA_position_with_0_gaps;
double SpartaABC_options::_wNumberOfMSA_position_with_1_gaps;
double SpartaABC_options::_wNumberOfMSA_position_with_2_gaps;
double SpartaABC_options::_wNumberOfMSA_position_with_n_minus_1_gaps;

string SpartaABC_options::_indelibleTemplateControlFile;
string SpartaABC_options::_dawgTemplateControlFile;
string SpartaABC_options::_inputRealMSAFile;
string SpartaABC_options::_outputGoodParamsFile;
string SpartaABC_options::_priorDistTypeRL;
string SpartaABC_options::_priorDistTypeA;
string SpartaABC_options::_priorDistTypeIR;
string SpartaABC_options::_inputTreeFileName;
int SpartaABC_options::_minRLVal;
int SpartaABC_options::_maxRLVal;
double SpartaABC_options::_minAVal;
double SpartaABC_options::_maxAVal;
double SpartaABC_options::_minIRVal;
double SpartaABC_options::_maxIRVal;
int SpartaABC_options::_numberOfSamplesToKeep;
double SpartaABC_options::_distanceCutOff;
int SpartaABC_options::_alignmentMode;
int SpartaABC_options::_similarityMode;
bool SpartaABC_options::_dawgSimulator;


SpartaABC_options::~SpartaABC_options() 
{}

void SpartaABC_options::setRLPrior(int msa_min_len, int msa_max_len) {
	_minRLVal = 0.8 * msa_min_len;
	_maxRLVal = 1.1 * msa_max_len;
}


void SpartaABC_options::initOptions(const string& paramFileName){
	initDefault();
	getParamsFromFile(paramFileName);
}

// This function initiates the diffult values for all parameters
void SpartaABC_options::initDefault() {
	// default values for the weight of each summary statistics.
	// the weights are computed within the "distance" function
	_wAvgGapSize = 1.0; // weight of the average gap size
	_wAvgUniqueGapSize = 1.0; // weight of the average gap size for unique gaps
	_wMSALen = 1.0; // weight of the total MSA length
	_wMSAMax = 1.0; // weight of the longest sequence within the MSA (without gaps)
	_wMSAMin = 1.0; // weight of the shortest sequence within the MSA (without gaps)
	_wTotNumGaps = 1.0;
	_wTotNumUniqueGaps = 1.0;
	_wNumGapsLenOne = 1.0;
	_wNumGapsLenTwo = 1.0;
	_wNumGapsLenThree = 1.0;
	_wNumGapsLenAtLeastFour = 1.0;

	_wNumberOfMSA_position_with_0_gaps = 1.0;;
	_wNumberOfMSA_position_with_1_gaps = 1.0;;
	_wNumberOfMSA_position_with_2_gaps = 1.0;;
	_wNumberOfMSA_position_with_n_minus_1_gaps = 1.0;


	// control files and tree file and real MSA file
	_indelibleTemplateControlFile = "";
	_inputTreeFileName = "";
	_dawgTemplateControlFile = "";
    _inputRealMSAFile = "";

	_outputGoodParamsFile = "";

	// assumptions about the prior distribution
	_priorDistTypeRL = "uniform";
	_priorDistTypeA = "uniform";
	_priorDistTypeIR = "uniform";

	// The minimum and maximum value for the prior distributions
	_minRLVal = 190;
	_maxRLVal = 300;
	_minAVal = 1.001;
	_maxAVal = 2.0;
	_minIRVal = 0.01;
	_maxIRVal = 0.15;

	_numberOfSamplesToKeep = 10;
	_distanceCutOff = DBL_MAX;

	_alignmentMode = 0; //multiple sequence alignment
	_similarityMode = 0; //simple mode (1-nuc, 2-aa)

	_dawgSimulator = false;
	
	Parameters::addParameter("_wAvgGapSize",_wAvgGapSize);
	Parameters::addParameter("_wAvgUniqueGapSize",_wAvgUniqueGapSize);
	Parameters::addParameter("_wMSALen",_wMSALen);
	Parameters::addParameter("_wMSAMax",_wMSAMax);
	Parameters::addParameter("_wMSAMin",_wMSAMin);
	Parameters::addParameter("_wTotNumGaps",_wTotNumGaps);
	Parameters::addParameter("_wTotNumUniqueGaps",_wTotNumUniqueGaps);
	Parameters::addParameter("_wNumGapsLenOne",_wNumGapsLenOne);
	Parameters::addParameter("_wNumGapsLenTwo",_wNumGapsLenTwo);
	Parameters::addParameter("_wNumGapsLenThree",_wNumGapsLenThree);
	Parameters::addParameter("_wNumGapsLenAtLeastFour",_wNumGapsLenAtLeastFour);

	Parameters::addParameter("_wNnumberOfMSA_position_with_0_gaps", _wNumberOfMSA_position_with_0_gaps);
	Parameters::addParameter("_wNnumberOfMSA_position_with_1_gaps", _wNumberOfMSA_position_with_1_gaps);
	Parameters::addParameter("_wNnumberOfMSA_position_with_2_gaps", _wNumberOfMSA_position_with_2_gaps);
	Parameters::addParameter("_wNumberOfMSA_position_with_n_minus_1_gaps", _wNumberOfMSA_position_with_n_minus_1_gaps);



	Parameters::addParameter("_indelibleTemplateControlFile", _indelibleTemplateControlFile);
	Parameters::addParameter("_dawgTemplateControlFile", _dawgTemplateControlFile);
	Parameters::addParameter("_inputRealMSAFile", _inputRealMSAFile);
	Parameters::addParameter("_inputTreeFileName", _inputTreeFileName);
	Parameters::addParameter("_outputGoodParamsFile", _outputGoodParamsFile);

	Parameters::addParameter("_priorDistTypeRL", _priorDistTypeRL);
	Parameters::addParameter("_priorDistTypeA", _priorDistTypeA);
	Parameters::addParameter("_priorDistTypeIR", _priorDistTypeIR);

	Parameters::addParameter("_minRLVal", _minRLVal);
	Parameters::addParameter("_maxRLVal", _maxRLVal);
	Parameters::addParameter("_minAVal", _minAVal);
	Parameters::addParameter("_maxAVal", _maxAVal);
	Parameters::addParameter("_minIRVal", _minIRVal);
	Parameters::addParameter("_maxIRVal", _maxIRVal);

	Parameters::addParameter("_numberOfSamplesToKeep", _numberOfSamplesToKeep);
	Parameters::addParameter("_distanceCutOff", _distanceCutOff);

	Parameters::addParameter("_alignmentMode", _alignmentMode);
	Parameters::addParameter("_similarityMode", _similarityMode);

	Parameters::addParameter("_dawgSimulator", (_dawgSimulator == true) ? 1 : 0); // if set to 1 then _dawgSimulator is true
}


void SpartaABC_options::getParamsFromFile(const string& paramFileName)
{
	ifstream params(paramFileName.c_str());
	if(params.fail())
        errorMsg::reportError("cannot open parameter file: " + paramFileName);
	if(params.good())
        Parameters::readParameters(params);
	params.close();

	
	_wAvgGapSize = Parameters::getFloat("_wAvgGapSize");
	_wAvgUniqueGapSize = Parameters::getFloat("_wAvgUniqueGapSize");
	_wMSALen = Parameters::getFloat("_wMSALen");
	_wMSAMax = Parameters::getFloat("_wMSAMax");
	_wMSAMin = Parameters::getFloat("_wMSAMin");
	_wTotNumGaps = Parameters::getFloat("_wTotNumGaps");
	_wTotNumUniqueGaps = Parameters::getFloat("_wTotNumUniqueGaps");
	_wNumGapsLenOne = Parameters::getFloat("_wNumGapsLenOne");
	_wNumGapsLenTwo = Parameters::getFloat("_wNumGapsLenTwo");
	_wNumGapsLenThree = Parameters::getFloat("_wNumGapsLenThree");
	_wNumGapsLenAtLeastFour = Parameters::getFloat("_wNumGapsLenAtLeastFour");
	_wNumberOfMSA_position_with_0_gaps = Parameters::getFloat("_wNumberOfMSA_position_with_0_gaps");
	_wNumberOfMSA_position_with_1_gaps = Parameters::getFloat("_wNumberOfMSA_position_with_1_gaps");
	_wNumberOfMSA_position_with_2_gaps = Parameters::getFloat("_wNumberOfMSA_position_with_2_gaps");
	_wNumberOfMSA_position_with_n_minus_1_gaps = Parameters::getFloat("_wNumberOfMSA_position_with_n_minus_1_gaps");

	_indelibleTemplateControlFile = Parameters::getString("_indelibleTemplateControlFile");
	_dawgTemplateControlFile = Parameters::getString("_dawgTemplateControlFile");
	_inputRealMSAFile = Parameters::getString("_inputRealMSAFile");
	_outputGoodParamsFile = Parameters::getString("_outputGoodParamsFile");
	_inputTreeFileName = Parameters::getString("_inputTreeFileName");

	_priorDistTypeRL = Parameters::getString("_priorDistTypeRL");
	_priorDistTypeA = Parameters::getString("_priorDistTypeA");
	_priorDistTypeIR = Parameters::getString("_priorDistTypeIR");

	_minRLVal = Parameters::getInt("_minRLVal");
	_maxRLVal = Parameters::getInt("_maxRLVal");
	_minAVal = Parameters::getFloat("_minAVal");
	_maxAVal = Parameters::getFloat("_maxAVal");
	_minIRVal = Parameters::getFloat("_minIRVal");
	_maxIRVal = Parameters::getFloat("_maxIRVal");

	_numberOfSamplesToKeep = Parameters::getInt("_numberOfSamplesToKeep");
	_distanceCutOff = Parameters::getFloat("_distanceCutOff");

	_alignmentMode = Parameters::getInt("_alignmentMode");
	_similarityMode = Parameters::getInt("_similarityMode");

	_dawgSimulator = (Parameters::getInt("_dawgSimulator") == 1) ? true : false;

}
