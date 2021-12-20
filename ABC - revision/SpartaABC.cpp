////
////  main.cpp
////  SpartaABC
////	Integrate ABC concepts to the SPARTA algorithm
////	
//#include "simulator.h"
//#include "SpartaABC.h"
//
//
//int simNo=0;
//size_t simulationCounter = 0;
//
////vector<string> mainSimualteAlongAtree();
//
//vector<double> simulateSequencesAndReturnSummaryStatistics(double randRootLength,
//	double randAParam,
//	double randInsertRatio,
//	double randDeletionRatio) {
//	Simulator sim(randRootLength, randAParam, randInsertRatio, randDeletionRatio);
//	vector<string> simulatedSeqs = sim.simulateBasedOnTree(SpartaABC_options::_inputTreeFileName);
//	MSA simulated_MSA(simulatedSeqs);
//	//summ_stats_wrapper simulatedSummStats(simulated_MSA, SpartaABC_options::_alignmentMode, SpartaABC_options::_similarityMode);
//	//vector<double> summStatsSim = simulatedSummStats.getSummStatsValsVector();
//	vector<double> summStatsSim = getStatVec(simulated_MSA);
//	return summStatsSim;
//}
//
//int main(int argc, const char * argv[])
//{
//
//	srand((unsigned)(time(0))); // for random number generation. Do not delete.
//
//	if (argc < 2) {// argc is 1 (the name of the program) plus the number of arguments
//		cout<<"SpartaABC requires a command file as an argument and no such file was provided."<<endl;
//		exit(1);
//	}
//
//	string spartaABCConfigFile = argv[1];
//	SpartaABC_options::initOptions(spartaABCConfigFile);
//	
//	//summ_stats_wrapper realSummStats(SpartaABC_options::_inputRealMSAFile,SpartaABC_options::_alignmentMode,SpartaABC_options::_similarityMode);
//	//vector<double> weightStats = realSummStats.getSummStatsWeightsVector();
//	//vector<double> summStatsReal = realSummStats.getSummStatsValsVector();
//	
//	vector<double> summStatsReal = getStatVec(SpartaABC_options::_inputRealMSAFile);
//
//	SpartaABC_options::setRLPrior(summStatsReal[3],summStatsReal[2]);
//
//	
//
//	//open result file where the parameter combinations that yielded small enough distance are kept
//	ofstream resFile;
//	resFile.open(SpartaABC_options::_outputGoodParamsFile);
//	if (!resFile) 
//	{
//       cout<<"couldn't open output file: "<<SpartaABC_options::_outputGoodParamsFile<<endl;
//       exit(1);
//	}
//	
//	// print params and summ stats header
//	string niceHeader = NiceHeader();
//	resFile<<niceHeader<<endl;
//	
//	//print statistics for input MSA
//	resFile<<"input_msa\t?\t?\t?\t?";
//	for(int i=0; i<summStatsReal.size(); i++)
//	{
//		resFile<<"\t"<<summStatsReal[i];
//	}
//	resFile<<endl;
//
//	//read template control of Indelible/Dawg into vector<string>
//	vector<string> templateCtrl;
//	if (SpartaABC_options::_dawgSimulator)
//	{
//		templateCtrl = readIndelibleTemplateControlFile(SpartaABC_options::_dawgTemplateControlFile);
//	}
//	else
//	{
//		templateCtrl = readIndelibleTemplateControlFile(SpartaABC_options::_indelibleTemplateControlFile);
//	}
//
//	int numberOfCollectedSamples = 0;
//	int numberSimulations = 0;
//	vector<vector<double> > allSimulationsResults;
//
//	vector<double> weightStats = getWeightsVector();
//	for (int i = 0; i<weightStats.size(); i++)
//	{
//		resFile << "\t" << weightStats[i];
//	}
//	resFile << endl;
//
//	euclidean_distance euclideanDistObj;
//	euclideanDistObj.setWeightsVector(weightStats);
//	double currentCutOffValue = 100.0; // very large
//	size_t totalNumberOfSpartaTests = 100000;
//	//SpartaABC_options::_numberOfSamplesToKeep
//	while(simulationCounter < totalNumberOfSpartaTests)
//	{
//		simulationCounter++;
//		int randRootLength = getRandIntParamVal(SpartaABC_options::_priorDistTypeRL, SpartaABC_options::_minRLVal, SpartaABC_options::_maxRLVal);
//		double randAParam = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeA, SpartaABC_options::_minAVal, SpartaABC_options::_maxAVal);
//		//double randIndelRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
//		double randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
//		double randDeletionRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
//		while (randInsertRatio == 0 && randDeletionRatio == 0) {
//			 randInsertRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
//			 randDeletionRatio = getRandDoubleParamVal(SpartaABC_options::_priorDistTypeIR, SpartaABC_options::_minIRVal, SpartaABC_options::_maxIRVal);
//		}
//		cout << "Simulation number:\t" << numberSimulations << endl;// "\tKept thus far:\t" << numberOfCollectedSamples << "\tParams are:\tRL " << randRootLength << "\tA " << randAParam << "\tInsR " << randInsertRatio << "\tDelR " << randDeletionRatio << endl;
//
//		//MSA simulated_MSA = simulateSingleMSA(randRootLength, randAParam, randIndelRatio, templateCtrl);
//
//		// our own version to get an MSA
//		// 1. simulate along the tree
//
//		vector<double> summStatsSim = simulateSequencesAndReturnSummaryStatistics(randRootLength,randAParam,
//			 randInsertRatio,randDeletionRatio);
//
//		double euclideanDist = euclideanDistObj.computeWeightedDist(summStatsReal, summStatsSim);
//		
//		if(euclideanDist < currentCutOffValue)
//		{
//			vector<double> tmp;
//			tmp.push_back(randRootLength);
//			tmp.push_back(randAParam);
//			//tmp.push_back(randIndelRatio);
//			tmp.push_back(randInsertRatio);
//			tmp.push_back(randDeletionRatio);
//			tmp.push_back(euclideanDist); // not a parameter
//			size_t k = 0;
//			for (; k < allSimulationsResults.size(); ++k) {
//				if (allSimulationsResults[k][tmp.size() - 1] > euclideanDist) break;
//			}
//			allSimulationsResults.insert(allSimulationsResults.begin() + k, tmp);
//			if (allSimulationsResults.size() > SpartaABC_options::_numberOfSamplesToKeep) allSimulationsResults.pop_back();
//			//allSimulationsResults.push_back(tmp);
//			currentCutOffValue = allSimulationsResults[allSimulationsResults.size() - 1][tmp.size() - 1];
//			//resFile<<euclideanDist<<"\t"<<randRootLength<<"\t"<<randAParam<<"\t"<< randInsertRatio<<"\t"<< randDeletionRatio;
//			//for(int i=0; i<summStatsSim.size(); i++)
//			//{
//			//	resFile<<"\t"<<summStatsSim[i];
//			//}
//			//resFile<<endl;
//
//			//for (size_t p = 0; p < allSimulationsResults.size(); ++p) {
//			//	cout << allSimulationsResults[p][tmp.size() - 1] << " ";
//			//} 
//			//cout << endl;
//			//numberOfCollectedSamples++;
//		}
//		numberSimulations++;
//	}	
//
//	// here we will find the posterior expectation for each parameter
//	double posteriorRL=0.0;
//	//double posteriorIR=0.0;
//	double posteriorInsR=0.0;
//	double posteriorDelR=0.0;
//	double posteriorA=0.0;
//	for (size_t i = 0; i < allSimulationsResults.size(); ++i) {
//		posteriorRL+= allSimulationsResults[i][0];
//		posteriorA += allSimulationsResults[i][1];
//		//posteriorIR += allSimulationsResults[i][2];
//		posteriorInsR+=allSimulationsResults[i][2];
//		posteriorDelR += allSimulationsResults[i][3];
//		
//	}
//	posteriorRL /= allSimulationsResults.size();
//	//posteriorIR /= allSimulationsResults.size();
//	posteriorInsR /= allSimulationsResults.size();
//	posteriorDelR /= allSimulationsResults.size();
//	posteriorA /= allSimulationsResults.size();
//	resFile << endl;
//	resFile << "The posterior expectation for the root length (RL) is: " << posteriorRL << endl;
//	//resFile << "The posterior expectation for the indel rate (IR) is: " << posteriorIR << endl;
//	resFile << "The posterior expectation for the indel rate (InsR) is: " << posteriorInsR << endl;
//	resFile << "The posterior expectation for the indel rate (DelR) is: " << posteriorDelR << endl;
//	resFile << "The posterior expectation for the A parameter is: " << posteriorA << endl;
//	resFile << "Total number of simulations: " << simulationCounter << endl;
//	resFile.close();
//	
//    return 0;
//}
//
//vector<string> readIndelibleTemplateControlFile(string indelibleTemplateCtrlFile)
//{
//	vector<string> templateInstructionString;
//	ifstream InFile;
//	InFile.open(indelibleTemplateCtrlFile); 
//	if(!InFile.is_open()) 
//	{
//		cout<<"can't open control file "<<indelibleTemplateCtrlFile<<endl;
//		exit(1);
//	}
//	string line;
//	while(! InFile.eof()) 
//	{
//  		getline(InFile,line);
//		templateInstructionString.push_back(line);
//	}
//	return templateInstructionString;
//}
//
//
//MSA simulateSingleMSA(int rootLength, double a_param, double indel_rate_ratio, vector<string> templateCtrl)
//{
//
//	Simulation sim(indel_rate_ratio, rootLength, a_param);
//	// if compilation is with Dawg (still user can choose between INDELible and Dawg)
//	#ifdef WITH_DAWG
//	if (SpartaABC_options::_dawgSimulator)
//	{
//		sim.simulateDawgMSA(1, templateCtrl);
//	}
//	else
//	{
//		sim.simulateMSA(1, templateCtrl);
//	}
//	#elif WITH_INDELIBLE
//
//	// compilation is without Dawg
//	sim.simulateMSA(1, templateCtrl);
//	#else // USING OUR OWN SIMULATOR
//	sim.simulateMSA(1, templateCtrl);
//	#endif
//
//	MSA simulatedMSA = sim.msaVec[0];
//	return simulatedMSA;
//}
