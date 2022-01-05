
#include "simulator.h"

using namespace std;

const int MaxIndelSize = 50; // maximum length of insertion or deletion
//int currIdToInsert; //A.M used in SpartaABC for continuity of inserstion IDs may want something similar when adding duplication/loss models
vector<nucID> superSequence;
int tmparray[MaxIndelSize];

/*
explicit Simulator(vector<int> ChromosomeLength, double AParam, double InvertRate, double TranslocateRatio, double FusionRate,
	double FissionRate, double DuplicationRate, double LossRate, double randRootAparam, tree & t) :
	_rootChromosomes(ChromosomeLength), _A_param(AParam), _InvR(InvertRate), _TrR(TranslocateRatio), _FuR(FusionRate), _FiR(FissionRate),
	_DR(DuplicationRate), _LR(LossRate), zip(AParam, Re_params::_maxBlock), rootZip(randRootAparam, Re_params::_rootMaxFamilySize), _t(t) {};
	*/
Simulator::Simulator(vector<int> ChromosomeLength, double AParam, double InvertRate, double TranslocateRatio, double FusionRate,
	double FissionRate, double DuplicationRate, double LossRate, double randRootAparam, tree & t) {
	_rootChromosomes = ChromosomeLength;
	_A_param = AParam;
	_InvR = InvertRate;
	_TrR = TranslocateRatio;
	_FuR = FusionRate;
	_FiR = FissionRate;
	_DR = DuplicationRate;
	_LR = LossRate;
	zip = FastZip(AParam, Re_params::_maxBlock); 
	rootZip = FastZip(randRootAparam, Re_params::_rootMaxFamilySize);
}

vector<genomeType> Simulator::simulateBasedOnTree() {
	// In the beginning, each genome is stored as a vector 
	// of vectors of integeres, i.e., vector<vector <int> > for example
	// a root genome of length six in two chromosomes could be V = {{0,1,2},{3,4,5}}
	vector<genomeType> simulatedLeavesGenomes; 
	 
	//currIdToInsert = _rootLength; //A.M for insertion IDs see definition for details
	// Initating the root sequence
	genomeType ancestralGenome = generateRootGenomeNoDupWLOG(); //no duplications!! starting with all positive genes (as suggested by Itsik)
	//A.M uses the chromosomelegths vector to create the genome
	//A.M this function results in a genome which is a vector of vectors like this:
	//A.M ancestralSequence = {{0,1,2,3,...},{...},...,{...,N}}
	//tree t(treeFileName); // changed function so we get the tree object as input, thus reducing the times we read from file.
	simulateAlongTree(_t.getRoot(), ancestralGenome, simulatedLeavesGenomes);
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

// create a vector of all event vectors and return it. 22.12.21. 
// Might be changed later according to the relations between the C++ and python
vector<vector<int>> Simulator::get_event_counter_vectors() {
	vector<vector<int>> output_vector;
	output_vector.push_back(inv_counter_vec);
	output_vector.push_back(trans_counter_vec);
	output_vector.push_back(fis_counter_vec);
	output_vector.push_back(fus_counter_vec);
	return output_vector;
}


