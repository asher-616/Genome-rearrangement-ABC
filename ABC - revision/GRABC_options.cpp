#include "GRABC_options.h"



string GRABC_options::_inputGenomesFile;
string GRABC_options::_outputParamsFile;
string GRABC_options::_priorDistTypeRL;
string GRABC_options::_priorDistTypeA;
string GRABC_options::_priorDistTypeIR;
string GRABC_options::_inputTreeFileName;
int GRABC_options::_minGeneNumVal;
int GRABC_options::_maxGeneNumVal;
int GRABC_options::_minChNum;
int GRABC_options::_maxChNum;
int GRABC_options::_maxBlock;
double GRABC_options::_minAVal;
double GRABC_options::_maxAVal;
double GRABC_options::_minInvVal;
double GRABC_options::_maxInvVal;
double GRABC_options::_minTransVal;
double GRABC_options::_maxTransVal;
double GRABC_options::_distanceCutOff;
int GRABC_options::_numberOfSamplesToKeep; //A.M 8.4 added for now because it is used in main and can't figure how to remove without messing algorithm

string GRABC_options::_priorDistTypeFuR; //fusion rate distribution type
string GRABC_options::_priorDistTypeFiR; //fission rate distribution type

double GRABC_options::_minFissionVal;//fission rate range
double GRABC_options::_maxFissionVal;

double GRABC_options::_minFusionVal;//fusion rate range
double GRABC_options::_maxFusionVal;

string GRABC_options::_priorDistTypeDR; //gene duplication rate distribution type for the future
string GRABC_options::_priorDistTypeLR; //gene loss rate distribution type for the future

double GRABC_options::_minDuplicationVal;//gene duplication rate range for the future
double GRABC_options::_maxDuplicationVal;

double GRABC_options::_minLossVal;//gene loss rate range for the future
double GRABC_options::_maxLossVal;

double GRABC_options::_minRootA; //used for gene family size in the root
double GRABC_options::_maxRootA;
int GRABC_options::_rootMaxFamilySize;

// for testing fission and fusion numbers
int GRABC_options::_fus_counter = 0;
int GRABC_options::_fis_counter = 0;

GRABC_options::~GRABC_options()
{}

/* for now we get number of chromosomes and number of genes ranges as parameters. might decide to change that later
void GRABC_options::setRLPrior(int msa_min_len, int msa_max_len) {
	_minRLVal = 0.8 * msa_min_len;
	_maxRLVal = 1.1 * msa_max_len;
}
*/

void GRABC_options::initOptions(const string& paramFileName) {
	initDefault();
	getParamsFromFile(paramFileName);
}

// This function initiates the diffult values for all parameters
void GRABC_options::initDefault() {
	// default values for the weight of each summary statistics.
	// the weights are computed within the "distance" function
	//tree file and real genomes file
	_inputTreeFileName = "";
	_inputGenomesFile = "";

	_outputParamsFile = "";

	// assumptions about the prior distribution
	_priorDistTypeRL = "uniform";
	_priorDistTypeA = "uniform";
	_priorDistTypeIR = "uniform";//A.M maybe change name later

	_priorDistTypeFuR = "uniform";
	_priorDistTypeFiR = "uniform";

	_priorDistTypeDR = "uniform"; //gene duplication rate distribution type for the future
	_priorDistTypeLR = "uniform"; //gene loss rate distribution type for the future

	// The minimum and maximum value for the prior distributions
	_minGeneNumVal = 190;
	_maxGeneNumVal = 300;
	_minChNum = 4;
	_maxChNum = 18;
	_minAVal = 1.001; //notice that A value min/max are not taken as user parameters. that is because small changes in A val range can greatly increase computation time
	_maxAVal = 2.0;
	_minInvVal = 0.0;
	_maxInvVal = 1.5;
	_minTransVal = 0.0;
	_maxTransVal = 1.5;

	_minFissionVal = 0.00;//fission rate range
	_maxFissionVal = 0.03;
	_minFusionVal = 0.00;//fusion rate range
	_maxFusionVal = 0.8;

	//for the future
	_minDuplicationVal = 0.0;//gene duplication rate range for the future
	_maxDuplicationVal = 0.0;

	_minLossVal = 0.0;//gene loss rate range for the future
	_maxLossVal = 0.0;

	_maxBlock = 10;//chosen arbitrary defualt value for now
	_distanceCutOff = DBL_MAX;
	_numberOfSamplesToKeep = 50;

	//a param used for gene families sizes in the root . NOTE: root generation WILL be changed for a more accurate version later on,
	_minRootA = 1.001;
	_maxRootA = 2.0;
	_rootMaxFamilySize = 5; //chosen randomly

	Parameters::addParameter("_inputGenomesFile", _inputGenomesFile);
	Parameters::addParameter("_inputTreeFileName", _inputTreeFileName);
	Parameters::addParameter("_outputParamsFile", _outputParamsFile);

	Parameters::addParameter("_priorDistTypeRL", _priorDistTypeRL); //consider changing to something more genome oriented
	Parameters::addParameter("_priorDistTypeA", _priorDistTypeA);
	Parameters::addParameter("_priorDistTypeIR", _priorDistTypeIR); //consider to change to inversion/transposiotion

	Parameters::addParameter("_priorDistTypeFuR", _priorDistTypeFuR);
	Parameters::addParameter("_priorDistTypeFiR", _priorDistTypeFiR);

	Parameters::addParameter("_priorDistTypeDR", _priorDistTypeDR);//for the future
	Parameters::addParameter("_priorDistTypeLR", _priorDistTypeLR);





	Parameters::addParameter("_minGeneNumVal", _minGeneNumVal);
	Parameters::addParameter("_maxGeneNumVal", _maxGeneNumVal);
	Parameters::addParameter("_minChNum", _minChNum);
	Parameters::addParameter("_maxChNum", _maxChNum);
	Parameters::addParameter("_minAVal", _minAVal);
	Parameters::addParameter("_maxAVal", _maxAVal);
	Parameters::addParameter("_minInvVal", _minInvVal);
	Parameters::addParameter("_maxInvVal", _maxInvVal);
	Parameters::addParameter("_minTransVal", _minTransVal);
	Parameters::addParameter("_maxTransVal", _maxTransVal);
	Parameters::addParameter("_maxBlock", _maxBlock);

	Parameters::addParameter("_minFissionVal", _minFissionVal);
	Parameters::addParameter("_maxFissionVal", _maxFissionVal);
	Parameters::addParameter("_minFusionVal",  _minFusionVal);
	Parameters::addParameter("_maxFusionVal",  _maxFusionVal);


	Parameters::addParameter("_minDuplicationVal", _minDuplicationVal);// for the future
	Parameters::addParameter("_maxDuplicationVal", _maxDuplicationVal);
	Parameters::addParameter("_minLossVal",        _minLossVal);
	Parameters::addParameter("_maxLossVal",        _maxLossVal);
	
	Parameters::addParameter("_distanceCutOff", _distanceCutOff);
	Parameters::addParameter("_numberOfSamplesToKeep", _numberOfSamplesToKeep);

}


void GRABC_options::getParamsFromFile(const string& paramFileName)
{
	ifstream params(paramFileName.c_str());
	if (params.fail())
		errorMsg::reportError("cannot open parameter file: " + paramFileName);
	if (params.good())
		Parameters::readParameters(params);
	params.close();


	_inputGenomesFile = Parameters::getString("_inputGenomesFile");
	_outputParamsFile = Parameters::getString("_outputParamsFile");
	_inputTreeFileName = Parameters::getString("_inputTreeFileName");

	_priorDistTypeRL = Parameters::getString("_priorDistTypeRL");
	_priorDistTypeA = Parameters::getString("_priorDistTypeA");
	_priorDistTypeIR = Parameters::getString("_priorDistTypeIR");

	_minGeneNumVal = Parameters::getInt("_minGeneNumVal");
	_maxGeneNumVal = Parameters::getInt("_maxGeneNumVal");
	_minChNum = Parameters::getInt("_minChNum"); 
	_maxChNum = Parameters::getInt("_maxChNum");
	_minAVal = Parameters::getFloat("_minAVal");
	_maxAVal = Parameters::getFloat("_maxAVal");
	_minInvVal = Parameters::getFloat("_minInvVal");
	_maxInvVal = Parameters::getFloat("_maxInvVal");
	_minTransVal = Parameters::getFloat("_minTransVal");
	_maxTransVal = Parameters::getFloat("_maxTransVal");

	_minFissionVal = Parameters::getFloat("_minFissionVal");
	_maxFissionVal = Parameters::getFloat("_maxFissionVal");
	_minFusionVal = Parameters::getFloat("_minFusionVal");
	_maxFusionVal = Parameters::getFloat("_maxFusionVal");

	_minDuplicationVal = Parameters::getFloat("_minDuplicationVal");//for the future
	_maxDuplicationVal = Parameters::getFloat("_maxDuplicationVal");
	_minLossVal = Parameters::getFloat("_minLossVal");
	_maxLossVal = Parameters::getFloat("_maxLossVal");


	_distanceCutOff = Parameters::getFloat("_distanceCutOff");
	_numberOfSamplesToKeep = Parameters::getInt("_numberOfSamplesToKeep");

	_maxBlock = Parameters::getInt("_maxBlock");
}
