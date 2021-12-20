
#include "summ_stats_wrapper.h"
#include "param_priors.h"




vector<vector<double> > simulateSequencesAndReturnSummaryStatistics(vector<int> randChromosomeLength, double randAParam, double randInvertRatio,
	double randTranslocateRatio, double randFusionRate, double randFissionRate, double randDuplicationRate, double randLossRate, double randRootAparam);
vector<double> tranformSumStatVector(vector<vector<int> > &sumStatOld, int maxBin);

vector<vector<double> > getStatVec(genomes &currGenomes) {
	return currGenomes.get_summary_stats_vector();
}

/*
vector<vector<int> > getStatVec() {
	vector<vector<int> > statVals;//does nothing. used as filler until I write a function that can extract sum stats from genomes file
	return statVals;
}*/

vector<vector<double> > getStatVec(string inputFile) { // this function return the summary statistics vector from an input MSA read from a file
	// curr_MSA(inputFile);
	//return getStatVec(curr_MSA);
	genomes realGenomes(inputFile);
	return getStatVec(realGenomes);
}

string NiceHeader() {
	//string niceHeader = "DISTANCE\tRL\tA\tIR\tAVG_GAP_SIZE\tMSA_LEN\tMSA_MAX_LEN\tMSA_MIN_LEN\tTOT_NUM_GAPS\tNUM_GAPS_LEN_ONE\tNUM_GAPS_LEN_TWO\tNUM_GAPS_LEN_THREE\tNUM_GAPS_LEN_AT_LEAST_FOUR\tAVG_UNIQUE_GAP_SIZE\tTOT_NUM_UNIQE_GAPS";
	string niceHeader = "DISTANCE\tRL\tA\tIR\tDR\tAVG_GAP_SIZE\tMSA_LEN\tMSA_MAX_LEN\tMSA_MIN_LEN\tTOT_NUM_GAPS\tNUM_GAPS_LEN_ONE\tNUM_GAPS_LEN_TWO\tNUM_GAPS_LEN_THREE\tNUM_GAPS_LEN_AT_LEAST_FOUR\tAVG_UNIQUE_GAP_SIZE\tTOT_NUM_UNIQE_GAPS";
	return niceHeader;
}

vector<double> getWeightsVectorUsingSimulations() {
	// simulate datasets from the prior
	const size_t numberOfSimulations = 5000;
	vector<double> summaryStatisticsSum;
	vector<double> summaryStatisticsSquareSum;
	vector<double> summStatistics;
	vector<vector<double> > preSummStatistics;
	for (size_t i = 0; i < numberOfSimulations; ++i) {
		//cout << '\n';
		cout << "weight simulation number " << i << endl;

		int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
		//A.M should be switched to chromosome number parameters when they are added

		int totalGenes = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal);
		vector<int> randChromosomeLength = generate_chromosome_sizes(randChromosomeNum, totalGenes);



		/*
		int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
		vector<int> randChromosomeLength;
		for (size_t i = 0; i < randChromosomeNum; i++)
		{
			randChromosomeLength.push_back(getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal));
		}
		int totalGenes = 0; //used to sum total genes to be added to parameters vector
		for (auto& n : randChromosomeLength)
			totalGenes += n;
		*/



		double randAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minAVal, GRABC_options::_maxAVal);
		double randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
		double randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
		//A.M inversion/translocation ratio
		double randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
		double randFissionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

		//for the future
		double randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
		double randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

		double randRootAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minRootA, GRABC_options::_maxRootA);

		while ((randInvertRatio == 0 && randTranslocateRatio == 0) && randFissionRate == 0 && randDuplicationRate == 0 && randLossRate == 0) {
			//A.M we want to avoid both being zero since we divide by their sum (also simulation will be pointless)
			//we do not check fusion rate because if every other rate is zero, we will get only fusion events
			//until only one chromosome left and we will get exponent error (exponent function doesn't like zero as input)
			randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
			randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
			randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
			randFissionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

			randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
			randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

		}
		cout << "\tChN " << randChromosomeNum << "\tParams are:\tRL " << totalGenes << "\tAparam " << randAParam << "\tInvR " << randInvertRatio << "\tDelR " << randTranslocateRatio << "\tfusionRate "<< randFusionRate<< "\tfissionRate "<< randFissionRate <<endl;
		preSummStatistics = simulateSequencesAndReturnSummaryStatistics(randChromosomeLength, randAParam, randInvertRatio, randTranslocateRatio,
			randFusionRate, randFissionRate, randDuplicationRate, randLossRate, randRootAParam);
		
		summStatistics = tranformSumStatVector(preSummStatistics,  GRABC_options::_maxBlock);
		
		if (i == 0) {
			summaryStatisticsSum.resize(summStatistics.size());
			summaryStatisticsSquareSum.resize(summStatistics.size());
		}

		for (size_t j = 0; j < summStatistics.size(); ++j) {
			summaryStatisticsSum[j] += summStatistics[j];
			summaryStatisticsSquareSum[j] += (summStatistics[j] * summStatistics[j]);
		}
	}// end of simulations
	for (size_t j = 0; j < summStatistics.size(); ++j) {
		summaryStatisticsSum[j] /= numberOfSimulations;
		summaryStatisticsSquareSum[j] /= numberOfSimulations;
	}
	vector<double> variance(summaryStatisticsSum.size());
	for (size_t j = 0; j < variance.size(); ++j) {
		variance[j] = summaryStatisticsSquareSum[j] - summaryStatisticsSum[j] * summaryStatisticsSum[j];
		variance[j] = sqrt(variance[j]); // now it is standard deviation
		variance[j] = 1.0 / variance[j]; // now it is weights.
	}
	return variance; // now the weights are updated.
}

/* A.M not used. we are using simulations only for generating weights
vector<double> getWeightsVector() {
	double wAvgGapSize = SpartaABC_options::_wAvgGapSize;
	double wMSALen = SpartaABC_options::_wMSALen;
	double wMSAMax = SpartaABC_options::_wMSAMax;
	double wMSAMin = SpartaABC_options::_wMSAMin;
	double wTotNumGaps = SpartaABC_options::_wTotNumGaps;
	double wNumGapsLenOne = SpartaABC_options::_wNumGapsLenOne;
	double wNumGapsLenTwo = SpartaABC_options::_wNumGapsLenTwo;
	double wNumGapsLenThree = SpartaABC_options::_wNumGapsLenThree;
	double wNumGapsLenAtLeastFour = SpartaABC_options::_wNumGapsLenAtLeastFour;
	double wAvgUniqueGapSize = SpartaABC_options::_wAvgUniqueGapSize;
	double wTotNumUniqueGaps = SpartaABC_options::_wTotNumUniqueGaps;
	double wNumberOfMSA_position_with_0_gaps = SpartaABC_options::_wNumberOfMSA_position_with_0_gaps;
	double wNumberOfMSA_position_with_1_gaps = SpartaABC_options::_wNumberOfMSA_position_with_1_gaps;
	double wNumberOfMSA_position_with_2_gaps = SpartaABC_options::_wNumberOfMSA_position_with_2_gaps;
	double wNumberOfMSA_position_with_n_minus_1_gaps = SpartaABC_options::_wNumberOfMSA_position_with_n_minus_1_gaps;

	if ((wAvgGapSize == -1.0) && (wMSALen == -1.0) && (wMSAMax == -1.0) && (wMSAMin == -1.0) && (wTotNumGaps == -1.0)
		&& (wNumGapsLenOne == -1.0) && (wNumGapsLenTwo == -1.0) && (wNumGapsLenThree == -1.0) && (wNumGapsLenAtLeastFour == -1.0)
		&& (wAvgUniqueGapSize == -1.0) && (wTotNumUniqueGaps == -1.0) && (wNumberOfMSA_position_with_0_gaps == -1.0)
		&& (wNumberOfMSA_position_with_1_gaps == -1.0) && (wNumberOfMSA_position_with_2_gaps == -1.0) && (wNumberOfMSA_position_with_n_minus_1_gaps == -1.0))
	{
		return getWeightsVectorUsingSimulations();
	}
	else if ((wAvgGapSize == -1.0) || (wMSALen == -1.0) || (wMSAMax == -1.0) || (wMSAMin == -1.0) || (wTotNumGaps == -1.0)
		|| (wNumGapsLenOne == -1.0) || (wNumGapsLenTwo == -1.0) || (wNumGapsLenThree == -1.0) || (wNumGapsLenAtLeastFour == -1.0)
		|| (wAvgUniqueGapSize == -1.0) || (wTotNumUniqueGaps == -1.0) || (wNumberOfMSA_position_with_0_gaps == -1.0)
		|| (wNumberOfMSA_position_with_1_gaps == -1.0) || (wNumberOfMSA_position_with_2_gaps == -1.0) || (wNumberOfMSA_position_with_n_minus_1_gaps == -1.0))
	{
		cout << " error. Either all weights are given or all of them not. You cannot give only some of the weights";
		exit(4);
	}
	else { // read from the input file
		vector<double> _summStatWeights;
		_summStatWeights.push_back(wAvgGapSize);
		_summStatWeights.push_back(wMSALen);
		_summStatWeights.push_back(wMSAMax);
		_summStatWeights.push_back(wMSAMin);
		_summStatWeights.push_back(wTotNumGaps);
		_summStatWeights.push_back(wNumGapsLenOne);
		_summStatWeights.push_back(wNumGapsLenTwo);
		_summStatWeights.push_back(wNumGapsLenThree);
		_summStatWeights.push_back(wNumGapsLenAtLeastFour);
		_summStatWeights.push_back(wAvgUniqueGapSize);
		_summStatWeights.push_back(wTotNumUniqueGaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_0_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_1_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_2_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_n_minus_1_gaps);
		return _summStatWeights;
	}

}*/
// ----------------------------- old version
/*
vector<double> simulateSequencesAndReturnSummaryStatistics(double randRootLength,
	double randAParam,
	double randInsertRatio,
	double randDeletionRatio);

void summ_stats_wrapper::computeAvgOnEachMsa()
{
	int numberOfMsas = _msas.size();
	_avgSummStatVals = getStatVec(_msas[0]);
	if(numberOfMsas > 1)
	{
		for(int i=1; i<numberOfMsas; i++)
		{
			vector<double> currSummStatsVals = getStatVec(_msas[i]);
			for(int j=0; j<_avgSummStatVals.size(); j++)
			{
				_avgSummStatVals[j] += currSummStatsVals[j];
			}
		}
		for(int j=0; j<_avgSummStatVals.size(); j++)
		{
			_avgSummStatVals[j] = _avgSummStatVals[j] / numberOfMsas;
		}
	}
}
*/
/*
vector<double> summ_stats_wrapper::getStatVec (MSA &currMSA)
{
	vector<double> statVals;
	
	statVals.push_back(currMSA.getStatValByType(AVG_GAP_SIZE));
	statVals.push_back(currMSA.getStatValByType(MSA_LEN));
	statVals.push_back(currMSA.getStatValByType(MSA_MAX_LEN));
	statVals.push_back(currMSA.getStatValByType(MSA_MIN_LEN));
	statVals.push_back(currMSA.getStatValByType(TOT_NUM_GAPS));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_ONE));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_TWO));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_THREE));
	statVals.push_back(currMSA.getStatValByType(NUM_GAPS_LEN_AT_LEAST_FOUR));
	statVals.push_back(currMSA.getStatValByType(AVG_UNIQUE_GAP_SIZE));
	statVals.push_back(currMSA.getStatValByType(TOT_NUM_UNIQUE_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_0_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_1_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_2_GAPS));
	statVals.push_back(currMSA.getStatValByType(MSA_POSITION_WITH_N_MINUS_1_GAPS));

	return statVals;
}

*/

/*
void summ_stats_wrapper::fillWeightsVec()
{

	double wAvgGapSize = SpartaABC_options::_wAvgGapSize;
	double wMSALen = SpartaABC_options::_wMSALen; 
	double wMSAMax = SpartaABC_options::_wMSAMax; 
	double wMSAMin = SpartaABC_options::_wMSAMin; 
	double wTotNumGaps = SpartaABC_options::_wTotNumGaps;
	double wNumGapsLenOne = SpartaABC_options::_wNumGapsLenOne;
	double wNumGapsLenTwo = SpartaABC_options::_wNumGapsLenTwo;
	double wNumGapsLenThree = SpartaABC_options::_wNumGapsLenThree;
	double wNumGapsLenAtLeastFour = SpartaABC_options::_wNumGapsLenAtLeastFour;
	double wAvgUniqueGapSize = SpartaABC_options::_wAvgUniqueGapSize;
	double wTotNumUniqueGaps = SpartaABC_options::_wTotNumUniqueGaps;
	double wNumberOfMSA_position_with_0_gaps = SpartaABC_options::_wNumberOfMSA_position_with_0_gaps;
	double wNumberOfMSA_position_with_1_gaps = SpartaABC_options::_wNumberOfMSA_position_with_1_gaps;
	double wNumberOfMSA_position_with_2_gaps = SpartaABC_options::_wNumberOfMSA_position_with_2_gaps;
	double wNumberOfMSA_position_with_n_minus_1_gaps = SpartaABC_options::_wNumberOfMSA_position_with_n_minus_1_gaps;
	
	if ((wAvgGapSize == -1.0) && (wMSALen == -1.0) && (wMSAMax == -1.0) && (wMSAMin == -1.0) && (wTotNumGaps == -1.0)
		&& (wNumGapsLenOne == -1.0) && (wNumGapsLenTwo == -1.0) && (wNumGapsLenThree == -1.0) && (wNumGapsLenAtLeastFour == -1.0)
		&& (wAvgUniqueGapSize == -1.0) && (wTotNumUniqueGaps == -1.0) && (wNumberOfMSA_position_with_0_gaps == -1.0)
		&& (wNumberOfMSA_position_with_1_gaps == -1.0) && (wNumberOfMSA_position_with_2_gaps == -1.0) && (wNumberOfMSA_position_with_n_minus_1_gaps == -1.0))
	{
		//eval...
		// simulate 100,000 datasets from the prior
		const size_t numberOfSimulations = 1000;
		vector<double> summaryStatisticsSum;
		vector<double> summaryStatisticsSquareSum;
		vector<double> summStatistics;
		for (size_t i = 0; i < numberOfSimulations; ++i) {
			cout << "simulation number " << i << endl;
			int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);
			double randAParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minAVal, SpartaABC_options::_maxAVal);
			double randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
			double randDeletionRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
			while (randInsertRatio == 0 && randDeletionRatio == 0) {
				randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
				randDeletionRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
			}

			summStatistics = simulateSequencesAndReturnSummaryStatistics(randRootLength,
				randAParam,
				randInsertRatio,
				randDeletionRatio);
			if (i == 0) {
				summaryStatisticsSum.resize(summStatistics.size());
				summaryStatisticsSquareSum.resize(summStatistics.size());
			}
			for (size_t j = 0; j < summStatistics.size(); ++j) {
				summaryStatisticsSum[j] += summStatistics[j];
				summaryStatisticsSquareSum[j] += (summStatistics[j]* summStatistics[j]);
			}
		}// end of simulations
		for (size_t j = 0; j < summStatistics.size(); ++j) {
			summaryStatisticsSum[j] /= numberOfSimulations;
			summaryStatisticsSquareSum[j] /= numberOfSimulations;
		}
		vector<double> variance(summaryStatisticsSum.size());
		for (size_t j = 0; j < variance.size(); ++j) {
			variance[j] = summaryStatisticsSum[j]* summaryStatisticsSum[j] - summaryStatisticsSquareSum[j];
			variance[j] = sqrt(variance[j]); // now it is standard deviation
			variance[j] = 1.0 / variance[j]; // now it is weights.
			_summStatWeights.push_back(variance[j]); // now the weights are updated.
		}
	}
	else if ((wAvgGapSize == -1.0) || (wMSALen == -1.0) || (wMSAMax == -1.0) || (wMSAMin == -1.0) || (wTotNumGaps == -1.0)
		|| (wNumGapsLenOne == -1.0) || (wNumGapsLenTwo == -1.0) || (wNumGapsLenThree == -1.0) || (wNumGapsLenAtLeastFour == -1.0)
		|| (wAvgUniqueGapSize == -1.0) || (wTotNumUniqueGaps == -1.0) || (wNumberOfMSA_position_with_0_gaps == -1.0)
		|| (wNumberOfMSA_position_with_1_gaps == -1.0) || (wNumberOfMSA_position_with_2_gaps == -1.0) || (wNumberOfMSA_position_with_n_minus_1_gaps == -1.0))
	{
		cout << " error. Either all weights are given or all of them not. You cannot give only some of the weights";
		exit(4);
	} 
	else {
		_summStatWeights.push_back(wAvgGapSize);
		_summStatWeights.push_back(wMSALen);
		_summStatWeights.push_back(wMSAMax);
		_summStatWeights.push_back(wMSAMin);
		_summStatWeights.push_back(wTotNumGaps);
		_summStatWeights.push_back(wNumGapsLenOne);
		_summStatWeights.push_back(wNumGapsLenTwo);
		_summStatWeights.push_back(wNumGapsLenThree);
		_summStatWeights.push_back(wNumGapsLenAtLeastFour);
		_summStatWeights.push_back(wAvgUniqueGapSize);
		_summStatWeights.push_back(wTotNumUniqueGaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_0_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_1_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_2_gaps);
		_summStatWeights.push_back(wNumberOfMSA_position_with_n_minus_1_gaps);
	}
}
*/

/*
string summ_stats_wrapper::getNiceHeader()
{
	string niceHeader = "DISTANCE\tRL\tA\tIR\tAVG_GAP_SIZE\tMSA_LEN\tMSA_MAX_LEN\tMSA_MIN_LEN\tTOT_NUM_GAPS\tNUM_GAPS_LEN_ONE\tNUM_GAPS_LEN_TWO\tNUM_GAPS_LEN_THREE\tNUM_GAPS_LEN_AT_LEAST_FOUR\tAVG_UNIQUE_GAP_SIZE\tTOT_NUM_UNIQE_GAPS";
	return niceHeader;
}
*/