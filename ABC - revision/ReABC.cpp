//
//  main.cpp
//  SpartaABC
//	Integrate ABC concepts to the SPARTA algorithm
//	

#include "ReABC.h"



//using namespace std;


int simNo = 0;
size_t simulationCounter = 0; //counts simulation performed so far 

//vector<string> mainSimualteAlongAtree();

vector<vector<double> > simulateSequencesAndReturnSummaryStatistics(vector<int> randChromosomeLength, double randAParam, double randInvertRatio,
	double randTranslocateRatio, double randFusionRate, double randFissionRate, double randDuplicationRate, double randLossRate, double randRootAparam) {//A.M insert and deletion should be changed to inversion and transposition
	int numberOfGeneFamilies = 0;
	for (auto& n : randChromosomeLength)
		numberOfGeneFamilies += n;
	Simulator sim(randChromosomeLength, randAParam, randInvertRatio, randTranslocateRatio, randFusionRate,
		randFissionRate, randDuplicationRate, randLossRate, randRootAparam);
	//cout << "start sim" << endl;
	vector<genomeType> simulatedGenomes = sim.simulateBasedOnTree(GRABC_options::_inputTreeFileName);

	//summ_stats_wrapper simulatedSummStats(simulated_MSA, SpartaABC_options::_alignmentMode, SpartaABC_options::_similarityMode);
	//vector<double> summStatsSim = simulatedSummStats.getSummStatsValsVector();	
	//A.M here I should create the 'genomes' class object
	/*this is used to comment/uncomment the print block. to uncomment add  * /  at the end of line 
	ofstream resFile;
	resFile.open(GRABC_options::_outputParamsFile);
	if (!resFile)
	{
		cout << "couldn't open output file: " << GRABC_options::_outputParamsFile << endl;
		exit(1);
	}
	for (size_t i = 0; i < simulatedGenomes.size(); i++) //used to print genomes for testing purposes
	{
		for (size_t j = 0; j < simulatedGenomes[i].size(); j++)
		{
			cout << simulatedGenomes[i][j].size() << ' ';
			for (size_t k = 0; k < simulatedGenomes[i][j].size(); k++)
			{
				resFile << ' ' << simulatedGenomes[i][j][k] ;
			}
			resFile << endl;
		}
		resFile << endl;
		cout << endl;
	}
	/**/
	//cout << "start SS" << endl;
	genomes genomesObj(simulatedGenomes);
	vector<vector<double> > summStatsSim = getStatVec(genomesObj);//A.M need to adjust function to project
	
	return summStatsSim;
}

/**/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@---OLD MAIN---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
int oldmain(int argc, const char * argv[])
{

	srand((unsigned)(time(0))); // for random number generation. Do not delete.

	if (argc < 2) {// argc is 1 (the name of the program) plus the number of arguments
		cout << "GR-ABC requires a command file as an argument and no such file was provided." << endl;
		exit(1);
	}

	string GRABCConfigFile = argv[1]; //A.M param file adjusted to genome
	GRABC_options::initOptions(GRABCConfigFile);//A.M adjusted to genome rearrangement

	//summ_stats_wrapper realSummStats(SpartaABC_options::_inputRealMSAFile,SpartaABC_options::_alignmentMode,SpartaABC_options::_similarityMode);
	//vector<double> weightStats = realSummStats.getSummStatsWeightsVector();
	//vector<double> summStatsReal = realSummStats.getSummStatsValsVector();

	vector<vector<double> > summStatsReal = getStatVec(GRABC_options::_inputGenomesFile); //A.M summary statistics for real data. adjusted to GR

	//SpartaABC_options::setRLPrior(summStatsReal[3], summStatsReal[2]); //A.M current version get root parameters in the param file. may be changed later 



	//open result file where the parameter combinations that yielded small enough distance are kept //A.M I don't understand this
	ofstream resFile;
	resFile.open(GRABC_options::_outputParamsFile);
	if (!resFile)
	{
		cout << "couldn't open output file: " << GRABC_options::_outputParamsFile << endl;
		exit(1);
	}

	// print params and summ stats header
	string niceHeader = NiceHeader(); //A.M change header
	resFile << niceHeader << endl;

	//print statistics for input MSA
	resFile << "input_msa\t?\t?\t?\t?";
	for (int i = 0; i < summStatsReal.size(); i++)
	{
		for (size_t j = 0; j < summStatsReal[i].size(); j++)
		{
			resFile << "\t" << summStatsReal[i][j];
		}
		resFile << "\n";
	}
	resFile << endl;
	vector<double> summStatRealModified = tranformSumStatVector(summStatsReal, GRABC_options::_maxBlock);

	int numberOfCollectedSamples = 0;
	int numberSimulations = 0;
	vector<vector<double> > allSimulationsResults;

	vector<double> weightStats = getWeightsVectorUsingSimulations(); //A.M changed it so that weight are generated only using simulations
	for (int i = 0; i < weightStats.size(); i++)
	{
		resFile << "\t" << weightStats[i];
	}
	resFile << endl;

	euclidean_distance euclideanDistObj;
	euclideanDistObj.setWeightsVector(weightStats); //A.M weights are used to calculate the (weighted) distance
	double currentCutOffValue = 100.0; // very large //A.M why are these set here and not as const in the h file?
	size_t totalNumberOfTests = 250000;
	//SpartaABC_options::_numberOfSamplesToKeep
	while (simulationCounter < totalNumberOfTests)
	{
		simulationCounter++;
		int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
		//A.M should be switched to chromosome number parameters when they are added
		
		int totalGenes = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal);
		vector<int> randChromosomeLength;
		int chromosome_size = totalGenes / randChromosomeNum;
		for (size_t i = 0; i < randChromosomeNum-1; i++)
		{
			randChromosomeLength.push_back(chromosome_size);

		}
		int last_chromosome = totalGenes - chromosome_size * (randChromosomeNum - 1); //putting the remainder in case totalGenes is not divisible by number of chromosomes
		randChromosomeLength.push_back(last_chromosome);

		//int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);//obsolete 
		//A.M generate number of chromosome (range from param file) and for each chromosome generates number of genes (again from param file)
		double randAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minAVal, GRABC_options::_maxAVal);
		//A.M a parameter is used for block length distribution (zipfian)
		//double randIndelRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		double randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
		double randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
		//A.M inversion/translocation ratio
		double randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
		double randFissionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

		//for the future
		double randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
		double randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

		double randRootAparam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minRootA, GRABC_options::_maxRootA);

		while (randInvertRatio == 0 && randTranslocateRatio == 0 &&  randFissionRate == 0 && randDuplicationRate == 0 && randLossRate == 0) {
			//A.M we want to avoid both being zero since we divide by their sum (also simulation will be pointless)
			randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
			randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
			randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
			randFissionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

			randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
			randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

		}
		cout << "Simulation number:\t" << numberSimulations << /*endl;//*/ "\tKept thus far:\t" << numberOfCollectedSamples << "\tChN " << randChromosomeNum << "\tParams are:\tRL " << totalGenes << "\tA " << randAParam << "\tInsR " << randInvertRatio << "\tDelR " << randTranslocateRatio << endl;

		//MSA simulated_MSA = simulateSingleMSA(randRootLength, randAParam, randIndelRatio, templateCtrl);

		// our own version to get an MSA
		// 1. simulate along the tree

		vector<vector<double> > summStatsSim = simulateSequencesAndReturnSummaryStatistics(randChromosomeLength, randAParam, randInvertRatio,
			randTranslocateRatio, randFusionRate, randFissionRate, randDuplicationRate, randLossRate, randRootAparam); //A.M in this file, above the main
		vector<double> summStatSimModified = tranformSumStatVector(summStatsSim, GRABC_options::_maxBlock);

		double euclideanDist = euclideanDistObj.computeWeightedDist(summStatRealModified, summStatSimModified);
		cout << euclideanDist << endl;

		if (euclideanDist < currentCutOffValue)
		{
			vector<double> tmp;
			tmp.push_back(randChromosomeNum); 
			tmp.push_back(totalGenes); 
			tmp.push_back(randAParam);
			//tmp.push_back(randIndelRatio);
			tmp.push_back(randInvertRatio);
			tmp.push_back(randTranslocateRatio);
			tmp.push_back(randFusionRate);
			tmp.push_back(randFissionRate);
			tmp.push_back(euclideanDist); // not a parameter
			size_t k = 0;
			for (; k < allSimulationsResults.size(); ++k) {
				if (allSimulationsResults[k][tmp.size() - 1] > euclideanDist) break;
			}
			allSimulationsResults.insert(allSimulationsResults.begin() + k, tmp);
			if (allSimulationsResults.size() > GRABC_options::_numberOfSamplesToKeep) {
				allSimulationsResults.pop_back();
				//allSimulationsResults.push_back(tmp);
				currentCutOffValue = allSimulationsResults[allSimulationsResults.size() - 1][tmp.size() - 1];
			}
			//resFile<<euclideanDist<<"\t"<<randRootLength<<"\t"<<randAParam<<"\t"<< randInsertRatio<<"\t"<< randDeletionRatio;
			//for(int i=0; i<summStatsSim.size(); i++)
			//{
			//	resFile<<"\t"<<summStatsSim[i];
			//}
			//resFile<<endl;

			//for (size_t p = 0; p < allSimulationsResults.size(); ++p) {
			//	cout << allSimulationsResults[p][tmp.size() - 1] << " ";
			//} 
			//cout << endl;
			numberOfCollectedSamples++;
		}
		numberSimulations++;
	}
	resFile << "root\tAparam\tinversion\ttransp\tfusion\tfission\tdistance";
	for (size_t i = 0; i < allSimulationsResults.size(); ++i) {
		resFile << endl;
		resFile << allSimulationsResults[i][0] << allSimulationsResults[i][1] << '\t' << allSimulationsResults[i][2] << '\t' << allSimulationsResults[i][3] << '\t' << allSimulationsResults[i][4] << '\t' << allSimulationsResults[i][5] << '\t' << allSimulationsResults[i][6] << '\t' << allSimulationsResults[i][7] << endl;
		 
	}

	// here we will find the posterior expectation for each parameter
	double posteriorRL = 0.0;
	double posteriorChN = 0.0;
	//double posteriorIR=0.0;
	double posteriorInvR = 0.0;
	double posteriorTransR = 0.0;
	double posteriorA = 0.0;
	double posteriorFus = 0.0;
	double posteriorFis = 0.0;
	for (size_t i = 0; i < allSimulationsResults.size(); ++i) {
		posteriorChN += allSimulationsResults[i][0];
		posteriorRL += allSimulationsResults[i][1];
		posteriorA += allSimulationsResults[i][2];
		//posteriorIR += allSimulationsResults[i][2];
		posteriorInvR += allSimulationsResults[i][3];
		posteriorTransR += allSimulationsResults[i][4];
		posteriorFus += allSimulationsResults[i][5];
		posteriorFis += allSimulationsResults[i][6];
	}
	posteriorChN /= allSimulationsResults.size();
	posteriorRL /= allSimulationsResults.size();
	//posteriorIR /= allSimulationsResults.size();
	posteriorInvR /= allSimulationsResults.size();
	posteriorTransR /= allSimulationsResults.size();
	posteriorA /= allSimulationsResults.size();
	posteriorFus /= allSimulationsResults.size();
	posteriorFis /= allSimulationsResults.size();

	resFile << endl;
	resFile << "The posterior expectation for the chromosome number (ChN) is: " << posteriorChN << endl;
	resFile << "The posterior expectation for the root length (RL) is: " << posteriorRL << endl;
	//resFile << "The posterior expectation for the indel rate (IR) is: " << posteriorIR << endl;
	resFile << "The posterior expectation for the inversion rate (InvR) is: " << posteriorInvR << endl;
	resFile << "The posterior expectation for the translocation rate (TransR) is: " << posteriorTransR << endl;
	resFile << "The posterior expectation for the fusion rate (FusR) is: " << posteriorFus << endl;
	resFile << "The posterior expectation for the fission rate (FisR) is: " << posteriorFis << endl;
	resFile << "The posterior expectation for the A parameter is: " << posteriorA << endl;
	resFile << "Total number of simulations: " << simulationCounter << endl;
	resFile.close();

	return 0;
}

vector<string> readIndelibleTemplateControlFile(string indelibleTemplateCtrlFile)
{
	vector<string> templateInstructionString;
	ifstream InFile;
	InFile.open(indelibleTemplateCtrlFile);
	if (!InFile.is_open())
	{
		cout << "can't open control file " << indelibleTemplateCtrlFile << endl;
		exit(1);
	}
	string line;
	while (!InFile.eof())
	{
		getline(InFile, line);
		templateInstructionString.push_back(line);
	}
	return templateInstructionString;
}




vector<double> tranformSumStatVector(vector<vector<double> >& sumStatOld, int maxBin) {
	//input is the summary statistics vector of vectors of ints (one for all blocks and the other for reversed only)
	//output is a vector of doubles in which all blocks of size maxBin and above are summed together making it a vector of size 2*(maxBin+1)
	vector<double> sumStatNew((maxBin)*(sumStatOld.size() - 1), 0);
	
	for (size_t i = 0; i < sumStatOld.size() - 1; i++)
	{
		int sumTail = 0; //for summing all bins bigger
		for (size_t j = 0; j < sumStatOld[i].size(); j++)
		{
			if (j < maxBin - 1)
				sumStatNew[maxBin*i + j ] = sumStatOld[i][j];
			else
			{
				sumTail += sumStatOld[i][j];
			}
			
		}
		sumStatNew[maxBin*(i + 1) - 1] = sumTail;
	}
	//last vector in sumStatOld is all values that are not distribution and can be inserted as is and inserted last not to intefere with distribution stats that have a tail that needs treatment
	size_t last_vect = sumStatOld.size() - 1; //last vect is not unique blocks but features for fusion-fission
	for (size_t j = 0; j < sumStatOld[last_vect].size(); j++)
	{
		sumStatNew.push_back(sumStatOld[last_vect][j]);
	}
	return sumStatNew;
}

vector < vector<double> > tranformSumStatVector2(vector<vector<double> >& sumStatOld, int maxBin) {
	//input is the summary statistics vector of vectors of ints (one for all blocks and the other for reversed only)
	//output is a vector of doubles in which all blocks of size maxBin and above are summed together making it a vector of size 2*(maxBin+1)
	vector<double> sumStatOrder((maxBin)*(sumStatOld.size() - 1), 0);

	for (size_t i = 0; i < sumStatOld.size() - 1; i++)
	{
		int sumTail = 0; //for summing all bins bigger
		for (size_t j = 0; j < sumStatOld[i].size(); j++)
		{
			if (j < maxBin - 1)
				sumStatOrder[maxBin*i + j] = sumStatOld[i][j];
			else
			{
				sumTail += sumStatOld[i][j];
			}

		}
		sumStatOrder[maxBin*(i + 1) - 1] = sumTail;
	}
	//last vector in sumStatOld is all values that are not distribution and can be inserted as is and inserted last not to intefere with distribution stats that have a tail that needs treatment
	size_t last_vect = sumStatOld.size() - 1; //last vect is not unique blocks but features for fusion-fission
	vector<double> sumStatweight;
	for (size_t j = 0; j < sumStatOld[last_vect].size(); j++)
	{
		sumStatweight.push_back(sumStatOld[last_vect][j]);
	}
	vector <vector <double> > output_vect{ sumStatOrder, sumStatweight };
	return output_vect;
}




vector<vector<double>> simulateRealData() {

	ofstream resFile;
	resFile.open(GRABC_options::_outputParamsFile);
	if (!resFile)
	{
		cout << "couldn't open output file: " << GRABC_options::_outputParamsFile << endl;
		exit(1);
	}



	//this function simulate a single data set.
	// it will generate random parameters based on the ranges (and distribution) supplied in the param file
	//it outputs the parameters used to file. maybe add option to output the entire genome and the summary stats


	int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
	//A.M should be switched to chromosome number parameters when they are added

	int totalGenes = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal);
	
	//generates chromosome sizes using (naive implementation of) multinomial distribution (with minimal size of 1 for each chromosome). 
	vector<int> randChromosomeLength = generate_chromosome_sizes(randChromosomeNum, totalGenes);

	/* old way to generate chromosome sizes by generating them one by one
	int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
	//A.M should be switched to chromosome number parameters when they are added
	vector<int> randChromosomeLength;
	for (size_t i = 0; i < randChromosomeNum; i++)
	{
		randChromosomeLength.push_back(getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal));
	}
	int totalGenes = 0; //used to sum total genes to be added to parameters vector
	for (auto& n : randChromosomeLength)
		totalGenes += n;
	*/

	//int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);//obsolete 
	//A.M generate number of chromosome (range from param file) and for each chromosome generates number of genes (again from param file)
	double randAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minAVal, GRABC_options::_maxAVal);
	//A.M a parameter is used for block length distribution (zipfian)
	//double randIndelRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
	double randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
	double randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
	//A.M inversion/translocation ratio
	double randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
	double randFissionRate =  getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

	//for the future
	double randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
	double randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

	double randRootAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minRootA, GRABC_options::_maxRootA);

	while (randInvertRatio == 0 && randTranslocateRatio == 0 && randFissionRate == 0 && randDuplicationRate == 0 && randLossRate == 0) {
		//A.M we want to avoid both being zero since we divide by their sum (also simulation will be pointless)
		randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
		randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
		randFusionRate =  getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
		randFissionRate =  getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

		randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
		randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

	}

	//resFile << "simulation parameters" << endl;
	//resFile << "chromosomeNum\troot\tAparam\tinversion\ttransp\tfusion\tfission"<< endl;
	resFile << "simulation number\tnumber of Chr\trootlen\ta param\tinvertion rate\ttranslocation rate\tfusion rate\tfission rate\tfusion num\tfission num" << endl;
	//resFile << "real\t" << randChromosomeNum << "\t" << totalGenes << "\t" << randAParam << "\t" << randInvertRatio << "\t" << randTranslocateRatio << "\t" << randFusionRate << "\t" << randFissionRate <<endl;

	cout << randChromosomeNum << "\t" << totalGenes << "\t" << randAParam << "\t" << randInvertRatio << "\t" << randTranslocateRatio << "\t" << randFusionRate << "\t" << randFissionRate << endl;

	vector<vector<double> > summStatsSim = simulateSequencesAndReturnSummaryStatistics(randChromosomeLength, randAParam, randInvertRatio,
		randTranslocateRatio, randFusionRate, randFissionRate, randDuplicationRate, randLossRate, randRootAParam); //A.M in this file, above the main
  	
	resFile << "real\t" << randChromosomeNum << "\t" << totalGenes << "\t" << randAParam << "\t" <<
		randInvertRatio << "\t" << randTranslocateRatio << "\t" << randFusionRate << "\t" << randFissionRate << "\t" <<
		GRABC_options::_fus_counter << "\t" << GRABC_options::_fis_counter << endl;

	vector<vector<double>> summStatSimModified = tranformSumStatVector2(summStatsSim, GRABC_options::_maxBlock);
	for (size_t i = 0; i < summStatSimModified[0].size(); i++)
	{
		resFile  << summStatSimModified[0][i] << "\t";
	}
	resFile << endl;
	for (size_t i = 0; i < summStatSimModified[1].size(); i++)
	{
		resFile  << summStatSimModified[1][i] << "\t";
	}
	resFile << endl;
	return summStatSimModified;

}



void RunABC(vector<vector<double>> & dataSummStat) {

	ofstream resFile;
	resFile.open(GRABC_options::_outputParamsFile, ios_base::app);
	if (!resFile)
	{
		cout << "couldn't open output file: " << GRABC_options::_outputParamsFile << endl;
		exit(1);
	}

	int numberOfCollectedSamples = 0;
	int numberSimulations = 0;
	vector<vector<double> > allSimulationsResults;

	//vector<double> weightStats = getWeightsVectorUsingSimulations(); //A.M changed it so that weight are generated only using simulations
	//for (int i = 0; i < weightStats.size(); i++)
	//{
	//	if (weightStats[i] == INFINITY) 
	//	{
	//		//for cases where some summ stats stay the same (mainly for cases with 1 chromosome and no fis/fus events)
	//		weightStats[i] = 0;
	//	}
	//	resFile << "\t" << weightStats[i];
	//}
	//resFile << endl;

	//euclidean_distance euclideanDistObj;
	//euclideanDistObj.setWeightsVector(weightStats); //A.M weights are used to calculate the (weighted) distance
	double currentCutOffValue = 100.0; // very large //A.M why are these set here and not as const in the h file?
	size_t totalNumberOfTests = 200000;
	//SpartaABC_options::_numberOfSamplesToKeep
	resFile << "simulation number\tnumber of Chr\trootlen\ta param\tinvertion rate\ttranslocation rate\tfusion rate\tfission rate"<< endl;
	while (simulationCounter < totalNumberOfTests)
	{
		simulationCounter++;
		//here we assume all chromosomes have the same size (last chromosome will recieve the reminder of the division)
		int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
		//A.M should be switched to chromosome number parameters when they are added

		int totalGenes = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal);
		vector<int> randChromosomeLength = generate_chromosome_sizes(randChromosomeNum, totalGenes);
		


		/*
		int randChromosomeNum = getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minChNum, GRABC_options::_maxChNum);
		//A.M should be switched to chromosome number parameters when they are added
		vector<int> randChromosomeLength;
		for (size_t i = 0; i < randChromosomeNum; i++)
		{
			randChromosomeLength.push_back(getRandIntParamVal(GRABC_options::_priorDistTypeRL, GRABC_options::_minGeneNumVal, GRABC_options::_maxGeneNumVal));
		}
		int totalGenes = 0; //used to sum total genes to be added to parameters vector
		for (auto& n : randChromosomeLength)
			totalGenes += n;

		*/
		//int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);//obsolete 
		//A.M generate number of chromosome (range from param file) and for each chromosome generates number of genes (again from param file)
		double randAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minAVal, GRABC_options::_maxAVal);
		//A.M a parameter is used for block length distribution (zipfian)
		//double randIndelRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
		double randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
		double randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
		//A.M inversion/translocation ratio
		double randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
		double randFissionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

		//for the future
		double randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
		double randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

		double randRootAParam = getRandDoubleParamVal(GRABC_options::_priorDistTypeA, GRABC_options::_minRootA, GRABC_options::_maxRootA);

		while (randInvertRatio == 0 && randTranslocateRatio == 0 && randFissionRate == 0 && randDuplicationRate == 0 && randLossRate == 0) {
			//A.M we want to avoid both being zero since we divide by their sum (also simulation will be pointless)
			randInvertRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minInvVal, GRABC_options::_maxInvVal);
			randTranslocateRatio = getRandDoubleParamVal(GRABC_options::_priorDistTypeIR, GRABC_options::_minTransVal, GRABC_options::_maxTransVal);
			randFusionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFuR, GRABC_options::_minFusionVal, GRABC_options::_maxFusionVal);
			randFissionRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeFiR, GRABC_options::_minFissionVal, GRABC_options::_maxFissionVal);

			randDuplicationRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeDR, GRABC_options::_minDuplicationVal, GRABC_options::_maxDuplicationVal);
			randLossRate = getRandDoubleParamVal(GRABC_options::_priorDistTypeLR, GRABC_options::_minLossVal, GRABC_options::_maxLossVal);

		}
		cout << "Simulation number:\t" << numberSimulations << /*endl;//*/ "\tKept thus far:\t" << numberOfCollectedSamples << "\tChN " << randChromosomeNum << "\tParams are:\tRL " << totalGenes << "\tA " << randAParam << "\tInsR " << randInvertRatio << "\tDelR " << randTranslocateRatio << endl;
		
		//MSA simulated_MSA = simulateSingleMSA(randRootLength, randAParam, randIndelRatio, templateCtrl);

		// our own version to get an MSA
		// 1. simulate along the tree

		vector<vector<double> > summStatsSim = simulateSequencesAndReturnSummaryStatistics(randChromosomeLength, randAParam, randInvertRatio,
			randTranslocateRatio, randFusionRate, randFissionRate, randDuplicationRate, randLossRate, randRootAParam); //A.M in this file, above the main
		
		resFile << numberSimulations << "\t" << randChromosomeNum << "\t" << totalGenes << "\t" << randAParam << "\t" <<
			randInvertRatio << "\t" << randTranslocateRatio << "\t" << randFusionRate << "\t" << randFissionRate << "\t" <<
			GRABC_options::_fus_counter << "\t" << GRABC_options::_fis_counter << endl;

		vector < vector<double> > summStatSimModified = tranformSumStatVector2(summStatsSim, GRABC_options::_maxBlock);
		for (size_t i = 0; i < summStatSimModified[0].size(); i++)
		{
			resFile << summStatSimModified[0][i] << "\t";
		}
		resFile << endl;
		for (size_t i = 0; i < summStatSimModified[1].size(); i++)
		{
			resFile << summStatSimModified[1][i] << "\t";
		}
		resFile << endl;
		//double euclideanDist = euclideanDistObj.computeWeightedDist(dataSummStat, summStatSimModified);
		//cout << euclideanDist << endl;

	//	if (euclideanDist < currentCutOffValue)
	//	{
	//		vector<double> tmp;
	//		tmp.push_back(randChromosomeNum);
	//		tmp.push_back(totalGenes);
	//		tmp.push_back(randAParam);
	//		//tmp.push_back(randIndelRatio);
	//		tmp.push_back(randInvertRatio);
	//		tmp.push_back(randTranslocateRatio);
	//		tmp.push_back(randFusionRate);
	//		tmp.push_back(randFissionRate);
	//		tmp.push_back(euclideanDist); // not a parameter
	//		int distance_index = 7; //used to fix the Dana bug
	//		size_t k = 0;
	//		for (; k < allSimulationsResults.size(); ++k) {
	//			if (allSimulationsResults[k][distance_index] > euclideanDist) break;
	//		}
	//		allSimulationsResults.insert(allSimulationsResults.begin() + k, tmp);
	//		if (allSimulationsResults.size() > GRABC_options::_numberOfSamplesToKeep) {
	//			allSimulationsResults.pop_back();
	//			//allSimulationsResults.push_back(tmp);
	//			currentCutOffValue = allSimulationsResults[allSimulationsResults.size() - 1][distance_index];
	//		}
	//		//resFile<<euclideanDist<<"\t"<<randRootLength<<"\t"<<randAParam<<"\t"<< randInsertRatio<<"\t"<< randDeletionRatio;
	//		//for(int i=0; i<summStatsSim.size(); i++)
	//		//{
	//		//	resFile<<"\t"<<summStatsSim[i];
	//		//}
	//		//resFile<<endl;

	//		//for (size_t p = 0; p < allSimulationsResults.size(); ++p) {
	//		//	cout << allSimulationsResults[p][tmp.size() - 1] << " ";
	//		//} 
	//		//cout << endl;
	//		numberOfCollectedSamples++;
	//	}
		numberSimulations++;
	//}
	//resFile << "root\tAparam\tinversion\ttransp\tfusion\tfission\tdistance";
	//for (size_t i = 0; i < allSimulationsResults.size(); ++i) {
	//	resFile << endl;
	//	resFile << allSimulationsResults[i][0] << '\t' << allSimulationsResults[i][1] << '\t' << allSimulationsResults[i][2] << '\t' << allSimulationsResults[i][3] << '\t' << allSimulationsResults[i][4] << '\t' << allSimulationsResults[i][5] << '\t' << allSimulationsResults[i][6] << '\t' << allSimulationsResults[i][7] << endl;

	}

	// here we will find the posterior expectation for each parameter
	//double posteriorRL = 0.0;
	//double posteriorChN = 0.0;
	////double posteriorIR=0.0;
	//double posteriorInvR = 0.0;
	//double posteriorTransR = 0.0;
	//double posteriorA = 0.0;
	//double posteriorFus = 0.0;
	//double posteriorFis = 0.0;
	//for (size_t i = 0; i < allSimulationsResults.size(); ++i) {
	//	posteriorChN += allSimulationsResults[i][0];
	//	posteriorRL += allSimulationsResults[i][1];
	//	posteriorA += allSimulationsResults[i][2];
	//	//posteriorIR += allSimulationsResults[i][2];
	//	posteriorInvR += allSimulationsResults[i][3];
	//	posteriorTransR += allSimulationsResults[i][4];
	//	posteriorFus += allSimulationsResults[i][5];
	//	posteriorFis += allSimulationsResults[i][6];
	//}
	//posteriorChN /= allSimulationsResults.size();
	//posteriorRL /= allSimulationsResults.size();
	////posteriorIR /= allSimulationsResults.size();
	//posteriorInvR /= allSimulationsResults.size();
	//posteriorTransR /= allSimulationsResults.size();
	//posteriorA /= allSimulationsResults.size();
	//posteriorFus /= allSimulationsResults.size();
	//posteriorFis /= allSimulationsResults.size();

	//resFile << endl;
	//resFile << "The posterior expectation for the chromosome number (ChN) is: " << posteriorChN << endl;
	//resFile << "The posterior expectation for the root length (RL) is: " << posteriorRL << endl;
	////resFile << "The posterior expectation for the indel rate (IR) is: " << posteriorIR << endl;
	//resFile << "The posterior expectation for the inversion rate (InvR) is: " << posteriorInvR << endl;
	//resFile << "The posterior expectation for the translocation rate (TransR) is: " << posteriorTransR << endl;
	//resFile << "The posterior expectation for the fusion rate (FusR) is: " << posteriorFus << endl;
	//resFile << "The posterior expectation for the fission rate (FisR) is: " << posteriorFis << endl;
	//resFile << "The posterior expectation for the A parameter is: " << posteriorA << endl;
	//resFile << "Total number of simulations: " << simulationCounter << endl;
	resFile.close();

	return;

}

vector<int> generate_chromosome_sizes(int chromosome_number, int genome_size) {
	// super naive way to create multinomial distribution
	// each chromosome start with size 1 (to prevent empty chromosomes) and the rest of the genome_size - chromosome_number genes
	// are each assigned a chromosome uniformly
	vector<int> chromosomes_vector(chromosome_number, 1);
	for (size_t i = 0; i < genome_size- chromosome_number; i++)
	{
		int rand_chromome = getRandIntParamVal(GRABC_options::_priorDistTypeRL, 0, chromosome_number-1);
		// note that range is 0-(chr_num-1) as vectors start at 0
		chromosomes_vector[rand_chromome] ++;
	}
	return chromosomes_vector;
}




//#####################################################
//################### NEW MAIN ######################## no duplication with all SS
//#####################################################
int main(int argc, const char * argv[]) {

 	srand((unsigned)(time(0))); // for random number generation. Do not delete.

	if (argc < 2) {// argc is 1 (the name of the program) plus the number of arguments. therfore 2 because we need the param file
		cout << "GR-ABC requires a parameter file as an argument and no such file was provided." << endl;
		exit(1);
	}

	string GRABCConfigFile = argv[1]; //A.M param file was adjusted to genome
	GRABC_options::initOptions(GRABCConfigFile);//A.M adjusted to genome rearrangement
	
	vector<vector<double>> summStatRealModified;

	if (GRABC_options::_inputGenomesFile != "") //real data input
	{
		ofstream resFile;
		resFile.open(GRABC_options::_outputParamsFile);
		if (!resFile)
		{
			cout << "couldn't open output file: " << GRABC_options::_outputParamsFile << endl;
			exit(1);
		}
		resFile << "ABC calculated using genomes data from: " << GRABC_options::_inputGenomesFile << endl;
		vector<vector<double> > summStatsReal = getStatVec(GRABC_options::_inputGenomesFile); //A.M summary statistics for real data. adjusted to GR
		summStatRealModified = tranformSumStatVector2(summStatsReal, GRABC_options::_maxBlock);
		resFile << "Real data" << endl;
		for (size_t i = 0; i < summStatRealModified[0].size(); i++)
		{
			resFile << summStatRealModified[0][i] << "\t";
		}
		resFile << endl;
		for (size_t i = 0; i < summStatRealModified[1].size(); i++)
		{
			resFile << summStatRealModified[1][i] << "\t";
		}
		resFile << endl;
	}
	else //if no input genome is supplied. a simulated genome will be used and parameters will be added in the output file (to compare with the results)
	{
		summStatRealModified = simulateRealData();
	}
	RunABC(summStatRealModified);



}