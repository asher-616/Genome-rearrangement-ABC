#include <iostream>
#include "tree.h"
#include "treeIt.h"
#include "simulator.h"

using namespace std;

const int MaxIndelSize = 50; // maximum length of insertion or deletion
//int currIdToInsert; //A.M used in SpartaABC for continuity of inserstion IDs may want something similar when adding duplication/loss models
vector<nucID> superSequence;
int tmparray[MaxIndelSize];


vector<genomeType> Simulator::simulateBasedOnTree(string & treeFileName) {
	// In the beginning, each genome is stored as a vector 
	// of vectors of integeres, i.e., vector<vector <int> > for example
	// a root genome of length six in two chromosomes would be V = {{0,1,2},{3,4,5}}
	vector<genomeType> simulatedLeavesGenomes; 
	 
	//added to test prediction of fusion and fission events. counting the number of each event and trying to predict them
	fis_counter = 0;
	fus_counter = 0;

	// Added for revision (15.12.21). Used to count the number of inversions and translocations in a tree.
	// Note- Tal prefers to count the events per branch but it demands some major changes to the code.
	// 21.12.21 changed to local variables (insread of static ones stored in the parameter object. 
	//Added vectors to save for each branch. initiation should be moved inside (for each branch instead of once per simulation)
	inv_counter = 0;
	trans_counter = 0;

	//currIdToInsert = _rootLength; //A.M for insertion IDs see definition for details
	// Initating the root sequence
	genomeType ancestralGenome = generateRootGenomeNoDupWLOG(); //no duplications!! starting with all positive genes (as suggested by Itsik)
//A.M uses the chromosomelegths vector to create the genome
	//A.M this function results in a genome which is a vector of vectors like this:
	//A.M ancestralSequence = {{0,1,2,3,...},{...},...,{...,N}}
	//superSequence = ancestralSequence; //A.M I don't need this.
	tree t(treeFileName);
	simulateAlongTree(t.getRoot(), ancestralGenome, simulatedLeavesGenomes);
	//A.M genomes are transferred to genomeClass objects after we have the leaves vector.
	
	return simulatedLeavesGenomes;
}

void Simulator::test() //used to test various class objects that cannot be accessed directly
{
	genomeType genome = generateRootGenome();
	genomeType genomeOutput;
	double branchLength = 1.3;
	simualteEventsAlongAspecificBranch(genome, branchLength, genomeOutput);
}


genomeType Simulator::generateRootGenomeOld() { //CHECH AND CHANGE XXXXXXXXXXXXXX
	// Uses predetermined chromosome sizes (drawn from multinomial distribution somewhere else)
	//A.M NOTE: a newer version exist that allows duplicates of the same gene (may be used for M2 model)
	genomeType genome;
	int GeneNum = 1; //used for continuity of gene numbers throughout the genome
	for (size_t i = 0; i < _rootChromosomes.size(); i++)
	{
		vector<int> tempChromosome;
		for (int j = 0; j < _rootChromosomes[i]; j++)
		{
			double orient = uniform();
			if (orient < 0.5)
				tempChromosome.push_back(GeneNum*(-1));
			else
				tempChromosome.push_back(GeneNum);
			GeneNum++;
		}
		genome.push_back(tempChromosome);
	}
	return genome;
}

genomeType Simulator::generateRootGenomeNoDupWLOG() { 
	// Uses predetermined chromosome sizes (drawn from multinomial distribution somewhere else)
	//A.M NOTE: a newer version exist that allows duplicates of the same gene (may be used for M2 model)
	// As suggested by Itsik Pe'er we start with positive only genome (without loss of generality)
	genomeType genome;
	int GeneNum = 1; //used for continuity of gene numbers throughout the genome
	for (size_t i = 0; i < _rootChromosomes.size(); i++)
	{
		vector<int> tempChromosome;
		for (int j = 0; j < _rootChromosomes[i]; j++)
		{
			tempChromosome.push_back(GeneNum);
			GeneNum++;
		}
		genome.push_back(tempChromosome);
	}
	return genome;
}


genomeType Simulator::generateRootGenome() {

	//will use zip to generate duplicates. will loop over genes and for each gene generate number of zip distribution {1,2...maxCopies}
	//for each gene that have >1 copies we will need to insert its duplicates randomly. this is done after generating the base genome
	genomeType genome;
	int GeneNum = 1; //used for continuity of gene numbers throughout the genome
	for (size_t i = 0; i < _rootChromosomes.size(); i++)
	{
		vector<int> tempChromosome;
		for (int j = 0; j < _rootChromosomes[i]; j++)
		{
			double orient = uniform();
			if (orient < 0.5)
				tempChromosome.push_back(GeneNum*(-1));
			else
				tempChromosome.push_back(GeneNum);
			GeneNum++;
		}
		genome.push_back(tempChromosome);
	} //generated base genome
	//now adding duplicates
	int genomeLength = GeneNum - 1;
	for (int i = 0; i < GeneNum; i++)
	{
		int duplicates = rootZip.drawZip();
		while (duplicates > 1)
		{
			vector<int> loc;
			drawRandomLocation(loc,genome, genomeLength);
			int chr = loc[0];
			int pos = loc[1];
			vector<int>::iterator it = genome[chr].begin();
			int gene = i;
			double orient = uniform();
			if (orient < 0.5)
				gene*= -1;
			genome[chr].insert(it + pos, gene);
			genomeLength++;
			duplicates--;
		}
	}
	return genome;
}


void printSequences(const vector<vector<nucID> >& v) {
	for (size_t i = 0; i < v.size(); ++i) {
		for (size_t j = 0; j < v[i].size(); ++j) {
			cout << v[i][j] << " ";
		}
		cout << endl;
	}
}

void Simulator::simualteEventsAlongAspecificBranch_rec(const genomeType& ancestralGenome,
													   double branchLength,
													   genomeType& outputGenome){ // the output
	//cout << "entering simulate along a branch" << endl;
	// The son sequence is initiated to be the father, and then, later, rearrangements are introduced
	outputGenome = ancestralGenome; //A.M is this a hard copy? why does it occur inside the function?

	// The waiting time for an event depends on the insertion rate and the deletion rate
	// The waiting time is for the entire sequence and not for a specific position
	// The _IR param is for position, so to get the sequence-leven insertion rate
	// one has to multiply the _IR parameter with the current length of the sequence.
	// Same for the DR.
	// However, assume a sequence of size 3 xyz
	// deletions can start before x, before y and before z -> 3 places.
	// insertions can start before x, before y, before z or after z -> 4 places
	unsigned int genomeLength = 0;// calculate genome length
	for (size_t i = 0; i < ancestralGenome.size(); i++)
	{
		genomeLength += ancestralGenome[i].size();
	}
	double sequenceWiseInvertionRate = 1.0*_InvR*(genomeLength); //A.M IR to be replaced by inversion rate (can be still named IR)
	double sequenceWiseTranspositionRate = 1.0*_TrR *genomeLength; //A.M DR to be replaced by transposition rate (maybe TR?)
	double fusionRateOverall = 1.0 * _FuR * (ancestralGenome.size()-1); //for one chromosome, fusion rate is always zero
	double fissionRateOverall = 1.0 * _FiR * genomeLength;
	double duplicationRateOverall = 1.0 * _DR * genomeLength;
	double lossRateOverall = 1.0 * _DR * genomeLength;
	double totalEventRate = sequenceWiseInvertionRate + sequenceWiseTranspositionRate + fusionRateOverall + fissionRateOverall + duplicationRateOverall + lossRateOverall;

	double waitingTime = drawExp(totalEventRate);
	if (waitingTime > branchLength)
		return;
	// generate one events, i.e., an insertion or a deletion
	//A.M need to figure best way to hard-copy genome
	SimulateEvent(outputGenome, genomeLength, sequenceWiseTranspositionRate,sequenceWiseInvertionRate, zip, fusionRateOverall, fissionRateOverall, 
		duplicationRateOverall, lossRateOverall);
	/*for (size_t i = 0; i < outputGenome.size(); i++) //for testing. used to observe the chromosome size after each event
	{
		cout << outputGenome[i].size() << " ";
	}
	cout << endl;*/
	/*
	bool isDeletion = true;
	
	if (uniform() < sequenceWiseInsertionRate / (sequenceWiseInsertionRate + sequenceWiseDeletionRate)) isDeletion = false;
	if (isDeletion) {
		//cout << "entering isDeletion ";
		int theStartingPoint = uniform(0, ancestralSequence.size() - 1);
		int deletionLength = drawZipf(_A_param, MaxIndelSize);
		if (theStartingPoint == 0 && deletionLength >= ancestralSequence.size()) {
			deletionLength = ancestralSequence.size()-1; // we do not want the entire sequence to be deleted...
		}
	//	cout << "deletion length=" << deletionLength<< endl;
		//int deletionLength = powerLaw(powerLawParam, N);
		size_t positionInWhichDeletionEnds = theStartingPoint + deletionLength;
		if (positionInWhichDeletionEnds > outputSeq.size()) positionInWhichDeletionEnds = outputSeq.size();
		outputSeq.erase(outputSeq.begin() + theStartingPoint, outputSeq.begin() + positionInWhichDeletionEnds);
		simualteWithIndelsAlongAspecificBranch(outputSeq, branchLength - waitingTime, outputSeq);
	//	cout << "exiting deletion" << endl;
		return;
	}
	else { // insertions
		//cout << "entering insertion ";
		int theStartingPoint = uniform(0, ancestralSequence.size()+1);
		int insertionLength = drawZipf(_A_param, MaxIndelSize);
		//cout << "insertion length=" << insertionLength << endl;
		//int insertionLength = powerLaw(powerLawParam, N);
		vector<nucID>::iterator it;
		if (theStartingPoint > ancestralSequence.size()) {
			it = outputSeq.end();
		}
		else {
			it = outputSeq.begin()+theStartingPoint;
		}
		//Initiating the array of new IDs to insert
		for (size_t m = 0; m < insertionLength; m++) {
			tmparray[m] = currIdToInsert;
			currIdToInsert++;
		}

		// this part updates the super sequence
		int numberBeforeInsertion = -1; // this indicates an insertion before 0.
		if (theStartingPoint > outputSeq.size()) {
			numberBeforeInsertion = outputSeq[outputSeq.size() - 1];
		}
		else {
			if (theStartingPoint >= 1) numberBeforeInsertion = outputSeq[theStartingPoint - 1];
		}
		int j = 0; 
		if (numberBeforeInsertion == -1) j = -1;
		else {
			for (; j < superSequence.size(); ++j) {
				if (superSequence[j] == numberBeforeInsertion) break;
			}
		}
		vector<nucID>::iterator it2 = superSequence.begin()+(j +1);
		superSequence.insert(it2, tmparray, tmparray+ insertionLength);

		// updating the sequence itself
		//it2 = it2 + insertionLength;
		outputSeq.insert(it, tmparray, tmparray + insertionLength);
		
		simualteEventsAlongAspecificBranch(outputGenome, branchLength - waitingTime, outputGenome);

	}*/
	simualteEventsAlongAspecificBranch(outputGenome, branchLength - waitingTime, outputGenome);

}

//new iterative version of "simulate along a specific branch"
void Simulator::simualteEventsAlongAspecificBranch(const genomeType& ancestralGenome, double branchLength,
	genomeType& outputGenome) {
	outputGenome = ancestralGenome;
	
	// initialize all counters
	inv_counter = 0;
	trans_counter = 0;
	fis_counter = 0;
	fus_counter = 0;
	
	do
	{
		unsigned int genomeLength = 0;// calculate genome length
		for (size_t i = 0; i < outputGenome.size(); i++)
		{
			genomeLength += outputGenome[i].size();
		}
		double sequenceWiseInvertionRate = 1.0*_InvR*(genomeLength); //A.M IR to be replaced by inversion rate (can be still named IR)
		double sequenceWiseTranspositionRate = 1.0*_TrR *genomeLength; //A.M DR to be replaced by transposition rate (maybe TR?)
		double fusionRateOverall = _FuR * outputGenome.size()*(outputGenome.size() - 1) * 2;
		//number of chromosomes choose 2 and multiplied by 4 because there are 4 fusion combinations (s-s, s-e, e-s, e-e) 
		double fissionRateOverall = 1.0 * _FiR * genomeLength;
		double duplicationRateOverall = 1.0 * _DR * genomeLength;
		double lossRateOverall = 1.0 * _DR * genomeLength;
		double totalEventRate = sequenceWiseInvertionRate + sequenceWiseTranspositionRate + fusionRateOverall + fissionRateOverall + duplicationRateOverall + lossRateOverall;
		double waitingTime = drawExp(totalEventRate);
		if (waitingTime > branchLength)
			break;
		branchLength -= waitingTime;
		SimulateEvent(outputGenome, genomeLength, sequenceWiseTranspositionRate, sequenceWiseInvertionRate, zip, fusionRateOverall, fissionRateOverall,
			duplicationRateOverall, lossRateOverall);




	} while (true);
	inv_counter_vec.push_back(inv_counter);
	trans_counter_vec.push_back(trans_counter);
	fis_counter_vec.push_back(fis_counter);
	fus_counter_vec.push_back(fus_counter);
	return;

}
// In each step, this function simulate the sons of node t and if it is a leaf
// it pushes the results to the simGenomes vector
void Simulator::simulateAlongTree(tree::nodeP t,
								  const genomeType & nodeGenome,
								  vector<genomeType> & simGenomes) {
	genomeType sonGenome;
	if (t->isLeaf()) { 
		simGenomes.push_back(nodeGenome);
		return;
	}
	else {// internal node}
		for (size_t i = 0; i < t->getNumberOfSons(); ++i) {
			// first we simulate the sequence of son i,
			// the results is stored in "outSequence"
			simualteEventsAlongAspecificBranch(nodeGenome, t->getSon(i)->dis2father(), sonGenome);
			simulateAlongTree(t->getSon(i), sonGenome, simGenomes);
		}
	}
}

/* A.M Not relevant to my project
void outputInFastaFormat(const vector<nucID> & leafSequence) {
	size_t j = 0;
	for (size_t i = 0; i < superSequence.size(); ++i) {
		if (j == leafSequence.size()) cout << "- ";
		else if (superSequence[i] == leafSequence[j]) {
			cout << superSequence[i]<<" ";
			++j;
		}
		else {
			cout << "- ";
		}
	}
}

void outputInFastaFormat(const vector<vector<nucID> > & simulatedLeavesSequences) {
	for (size_t i = 0; i < simulatedLeavesSequences.size(); ++i) {
		cout << ">S" << i << endl;
		outputInFastaFormat(simulatedLeavesSequences[i]);
		cout << endl;
	}
}

void simulationToMSA(const vector<vector<nucID> > & simulatedLeavesSequences, vector<string> & msa) {
	for (size_t l = 0; l < simulatedLeavesSequences.size(); ++l) {//go over all the leaves
		string tmp;												  // reset string for the current sequence
		size_t j = 0;											  // reset index for the current sequence
		for (size_t i = 0; i < superSequence.size(); ++i) {		  // go over the superSequence
			if (j == simulatedLeavesSequences[l].size()) {		  // if currenr sequence ended, add "-"
				tmp.append("-");
			}
			else if (superSequence[i] == simulatedLeavesSequences[l][j]) { // if current sequence = superSequence
				tmp.append("A");											// add arbitrary character
				++j;
			}
			else {														  // else add "-"
				tmp.append("-");
			}
		}
		msa.push_back(tmp);												// add current string to MSA
	}
}
*/



//vector<string> Simulator::mainSimualteAlongAtree()
//{
//	int length = 3000;
//	vector<vector<nucID> > simulatedLeavesSequences;
//	vector<nucID> ancestralSequence;
//	int currIdToInsert = length;
//	for (int i = 0; i < length; ++i) {
//		ancestralSequence.push_back(i);
//	}
//	superSequence = ancestralSequence;
//	string treeFile = "C:\\Users\\Dana Rapoport\\Downloads\\SpartaABC_with_INDELible_20170320_bundle\\SpartaABC_with_INDELible_20170320_bundle\\run_example\\small_tree.txt";
//	tree t(treeFile);
//	simulateAlongTree(t.getRoot(), ancestralSequence, simulatedLeavesSequences, 1.3, 0.02, 0.02);
//	vector<string> msa;
//	simulationToMSA(simulatedLeavesSequences, msa);
//	return msa;
//	//for (int k = 0; k < msa.size(); ++k) { cout << msa[k] << endl; }
//	/*printSequences(simulatedLeavesSequences);
//	// printing the super sequence:
//	cout << "\nthe super sequence is: " << endl;
//	for (size_t k = 0; k < superSequence.size(); ++k) {
//		cout << superSequence[k] << " ";
//	}
//	cout << endl;
//	cout << " the final leaf sequences in Fasta format\n";
//	outputInFastaFormat(simulatedLeavesSequences);
//	*/
//	//return 0;
//}

/*
int main2() {

	size_t rootLength = 50;
	double A_param = 1.2;
	double IR = 0.01;
	double DR = 0.03;
	//Simulator sim(rootLength, A_param, IR, DR);//changed builder parameters. should adjust if testing needed
	//sim.test();
}
*/