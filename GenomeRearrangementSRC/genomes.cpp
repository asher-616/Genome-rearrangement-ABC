#include "genomes.h"

#include <ctime>







void genomes::init(const vector<genomeType> & genomeVect) {
	for (vector<genomeType>::const_iterator it = genomeVect.begin(); it != genomeVect.end(); ++it)
	{
		//need to change each genome to genomeClass and add to vector
		GenomeClass tempGenome(*it);
		genomesVect.push_back(tempGenome);
	}
	//calcSumStatAllVSAll(); // old version
	calcSumStatBrothersOnly(); // new version
	
}

genomes::genomes(const vector<genomeType> &genomeVect, tree sim_tree){
	//A.M here I receive a vector of genomestype objects and create a 'genomes' class object by creating a 'genomeClass' object from
	//each genometype object
	_t = sim_tree;
	init(genomeVect);
}

genomes::genomes(string inputFile, string treeFile) {
	// get an input file (minimalistic format) read it and extracts summary statistics
	ifstream file(inputFile);
	string line;
	vector<genomeType> genomesVector;
	genomeType genome;
	chromosomeType chromosome;

	while (getline(file, line))
	{
		cout << line << endl;
		if (line.empty())
		{
			continue; //ignoring empty lines
		}

		if (line[0] == '>')
		{
			//reached a new genome. time to stow the previous one away
			if (!chromosome.empty())
			{
				genome.push_back(chromosome);
				chromosome.clear();
			}
			if (!genome.empty())
			{
				genomesVector.push_back(genome);
				genome.clear();
			}

			continue;//currently we don't need the genome name. may change later
		}

		if (line[0] == 'c' || line[0] == 'C')
		{
			//reached a new chromosome. time to stow the previous one away
			if (!chromosome.empty())
				genome.push_back(chromosome);
			chromosome.clear();
			int pos2remove = line.find_first_of(" \t");
			line.erase(0, pos2remove + 1); //clearing "chr# from beginning of line so I can add genes to the (empty) chromosome
		}
		//add current line to chromosome
		size_t pos = 0;
		string token;
		while ((pos = line.find("\t")) != string::npos) {
			token = line.substr(0, pos);
			//cout << genomesVector.size() + 1 << genome.size() + 1 << token << std::endl; //testing to make sure it looks ok
			line.erase(0, pos + 1);
			int gene = stoi(token, nullptr);
			chromosome.push_back(gene);
		}
		int gene = stoi(line, nullptr);
		chromosome.push_back(gene);
	}
	if (!chromosome.empty())
	{
		genome.push_back(chromosome);
		chromosome.clear();
	}
	if (!genome.empty())
	{
		genomesVector.push_back(genome);
		genome.clear();
	}

	_t = tree(treeFile);
	init(genomesVector);

}

vector<vector<double> > genomes::get_summary_stats_vector() {
	// returns a vector of vectors that contain all summary statistics
	vector<vector<double> > sumStat;
	vector<double> importantStats; //this will be the last values vector in the vector of vectors and will contain all stand-alone values 
	//unlike the rest which are all of the same type (showing block sizes distribution) the values in this array will be inserted ignoring the maxBin value used for the other vectors
	importantStats.push_back(minChromosomes);
	importantStats.push_back(maxChromosomes);
	importantStats.push_back(averageChromosomes);
	importantStats.push_back(varCromosomes);
	importantStats.push_back(totalTips);
	importantStats.push_back(totalChrs);
	importantStats.push_back(uniqueTipCounter);
	importantStats.push_back(FisFusCounter);
	vector<double>	uniqueBlocksVecDouble(uniqueBlocksVec.begin(), uniqueBlocksVec.end()); //changing to double to add to the same vector
	vector<double>	uniqueRevBlocksVecdouble(uniqueRevBlocksVec.begin(), uniqueRevBlocksVec.end());
	sumStat.push_back(uniqueBlocksVecDouble);
	sumStat.push_back(uniqueRevBlocksVecdouble);
	//sumStat.push_back(importantStats);
	sumStat.push_back(importantStats);
	cout << "min\tmax\tavg" << endl; // for testing
	cout << minChromosomes << "\t" << maxChromosomes << "\t" << averageChromosomes << endl; // for testing
	return sumStat;
}






void genomes::calcSumStatAllVSAll() {
	//going over every pair of genomes in genomesVect, extract data and summarize it in repective internal variables
	for (vector<GenomeClass>::iterator it = genomesVect.begin(); it != genomesVect.end()-1; ++it)
	{
		for (vector<GenomeClass>::iterator it2 = it + 1; it2 != genomesVect.end(); ++it2)
		{
			//print_genome(*it);
			//print_genome(*it2);
			//do pair calculations here
			compareBlocksGenomes(*it, *it2);
		}
	}
	size_t sumChromosomes = 0;

	for (vector<GenomeClass>::iterator it = genomesVect.begin(); it != genomesVect.end(); ++it)
	{
		//here I will calculate additional summary statistics used for fission and fusion
		size_t numberOfChromosomes = (*it).genome.size();
		all_chromosomes_vect.push_back(numberOfChromosomes);
		if (numberOfChromosomes < minChromosomes)
		{
			minChromosomes = numberOfChromosomes;
		}
		if (numberOfChromosomes > maxChromosomes)
		{
			maxChromosomes = numberOfChromosomes;
		}
		sumChromosomes += numberOfChromosomes;
	}
	averageChromosomes = 1.0 * sumChromosomes / genomesVect.size();
	for (vector<GenomeClass>::iterator it = genomesVect.begin(); it != genomesVect.end(); ++it)
	{
		//here I will calculate additional summary statistics used for fission and fusion
		size_t numberOfChromosomes = (*it).genome.size();
		varCromosomes += pow(numberOfChromosomes - averageChromosomes, 2);
	}
	varCromosomes /= genomesVect.size();
	resultType FisFusPairRes = genomes::calcFisFusBasedOnTree();
	uniqueTipCounter = FisFusPairRes.second.first.size();
	FisFusCounter = FisFusPairRes.second.second;
	totalTips = FisFusPairRes.first.first.size();
	totalChrs = FisFusPairRes.first.second.size();
}



void genomes::calcSumStatBrothersOnly() {
	//this function DOES NOT go over every pair of genomes in genomesVect to calculate unique-blocks.
	//this function only calculates unique blocks between bother leaves

	// to be replaced with a brother comparison. First I need to figure the following: WHY IS IT STILL ALL VS ALL?
	// 1. what to do with a leaf for which its brother is an internal node? use closest leaf descendent
	// 2. how to compare a leaf that has more than one brother? compare the one closest to the parent to all others
	// 3. what about events that occured in on a branch leading to an internal node? (won't show when comparing brother leaves only) - we use closest leaf to represent internal nodes


	calcUniqueBlocksBasedOnTree(); // needs testing
	size_t sumChromosomes = 0;


	for (vector<GenomeClass>::iterator it = genomesVect.begin(); it != genomesVect.end(); ++it)
	{
		//here I will calculate additional summary statistics used for fission and fusion
		size_t numberOfChromosomes = (*it).genome.size();
		all_chromosomes_vect.push_back(numberOfChromosomes);
		if (numberOfChromosomes < minChromosomes)
		{
			minChromosomes = numberOfChromosomes;
		}
		if (numberOfChromosomes > maxChromosomes)
		{
			maxChromosomes = numberOfChromosomes;
		}
		sumChromosomes += numberOfChromosomes;
	}
	averageChromosomes = 1.0 * sumChromosomes / genomesVect.size();
	for (vector<GenomeClass>::iterator it = genomesVect.begin(); it != genomesVect.end(); ++it)
	{
		//here I will calculate additional summary statistics used for fission and fusion
		size_t numberOfChromosomes = (*it).genome.size();
		varCromosomes += pow(numberOfChromosomes - averageChromosomes, 2);
	}
	varCromosomes /= genomesVect.size();
	resultType FisFusPairRes = genomes::calcFisFusBasedOnTree();
	uniqueTipCounter = FisFusPairRes.second.first.size();
	FisFusCounter = FisFusPairRes.second.second;
	totalTips = FisFusPairRes.first.first.size();
	totalChrs = FisFusPairRes.first.second.size();
}

void genomes::calcUniqueBlocksBasedOnTree() { //need testing
	//tree t(Re_params::_inputTreeFileName); // replaced to reduce file reading
	tree t = _t;
	vector<GenomeClass>::iterator genomeIt = genomesVect.begin();
	genomeAndLengthPair TreeResSet = calcUniqueBlocksAlongTree(t.getRoot(), genomeIt); // not final!
	// might consider add chromosome set and tip set sizes to the SS
	return ;
}
//genomeAndLengthPair genomes::calcUniqueBlocksAlongTree(tree::nodeP t, vector<GenomeClass>::iterator &genomeIt) 
// should go here. was moved temporarily to compare to fisfus along a tree (for reference)


vector<int> genomes::compareBlocksGenomes(GenomeClass& genome1, GenomeClass& genome2) {//need to 
	// gets two genomes and returns a vector of all common blocks sizes
	size_t chromosomes = genome1.genome.size();
	vector<int> blocksVec;


	for (size_t i = 0; i < chromosomes; i++)
	{
		int block = 1;
		int pre1 = 0; //will contain gene number before current block for inversion test
		int pre2 = genome2.prevGene(genome1.genome[i][0]);
		vector<int> blockseq; //contains all block sizes for the genome pair. right now we don't use it. either use or remove
		size_t chromosomeSize = genome1.genome[i].size();
		for (size_t j = 0; j < chromosomeSize; j++)
		{
			//cout << j << ' '; //for testing
			blockseq.push_back(genome1.genome[i][j]);
			if (j + 1 == chromosomeSize)  //reached end of chromosome
			{
				blocksVec.push_back(block); 


				updateUniqueBlocksHash(blockseq);


				if (((isReversed(genome1.genome[i][j],genome2)) && ((pre1 == genome2.nextGene(genome1.genome[i][j]) ))) 
					|| (!(isReversed(genome1.genome[i][j], genome2)) && (pre1 == -1* genome2.nextGene(genome1.genome[i][j]))))
				{
					updateUniqueRevBlocksHash(blockseq);
				}
			}
			else
			{
				int nextGene1 = genome1.genome[i][j + 1]; //next gene in genome 1
				int nextGene2 = genome2.nextGene(genome1.genome[i][j]);
				if (abs(nextGene1) ==abs(nextGene2))
				{
					block += 1;
				}
				else
				{
					blocksVec.push_back(block);
					updateUniqueBlocksHash(blockseq);
					

					bool a = isReversed(genome1.genome[i][j], genome2); //testing for following if conditions
					bool b = pre1 == nextGene2;
					bool c = pre2 == nextGene1;
					if ((isReversed(genome1.genome[i][j],genome2) && (pre1 == nextGene2) || pre2 == nextGene1) ||
						(isReversed(genome1.genome[i][j], genome2) && (pre1 == nextGene2) || pre2 == nextGene1))
					{
						updateUniqueRevBlocksHash(blockseq);
					}
					blockseq.clear();
					block = 1;
					pre1 = genome1.genome[i][j];
					pre2 = genome2.prevGene(genome1.genome[i][j+1]);
					
				}
			}
		}
	}
	return blocksVec;
}



bool genomes::isReversed(int gene, GenomeClass & genome2) {
	//gets a gene and a gnome and checks if gene appears in genome in revesed orientation
	if (gene * genome2.genomeArrangement[abs(gene) - 1][1]<0)
	{
		return true;
	}
	return false; 
}

/* Not sure if used. commented for now
void genomes::updateUniqueBlocksTrees(vector<int> &seq) {
	if (seq.size()< _minUniqueBlockSize || seq.size() > _maxUniqueBlockSize)
		return;
	if (abs(seq[0])> abs(seq.back()))
		reverseSubSeq(seq);
	updateTree(uniqueBlocksTree, seq);
}
void genomes::updateUniqueRevBlocksTrees(vector<int> &seq) {
	if (seq.size() < _minUniqueRevBlockSize || seq.size() > _maxUniqueRevBlockSize)
		return;
	if (seq[0] < seq.back())
		reverseSubSeq(seq);
	updateTree(uniqueRevBlocksTree, seq);
}
*/
/* regular revBlocks feature was cancelled, so probable revBlock is stored in revBlock now
void genomes::updateUniqueProbRevBlocksTrees(vector<int> seq) {
	if (seq.size() < _minUniqueProbRevBlockSize || seq.size() > _maxUniqueProbRevBlockSize)
		return;
	if (seq[0] < seq.back())
		reverseSubSeq(seq);
	updateTree(uniqueProbRevBlocksTree, seq);
}
*/
/*
void genomes::updateTree(Trie & seqTree, vector<int> &seq) {
	//get a tree and a sequence and add sequence to the tree
	//Note: sequence should arrive so that is lower edge number is on the beginning 
	//Note: might want to add a variable that counts number of leaves in a tree if so here is a good place to update it

		seqTree.insert(seq);
}
*/
void genomes::reverseSubSeq(vector<int> & seq) {
	//get a gene sequence (by ref) and reverse it (gene orientation reversed too)
	//e.g. input [1,-2,3] output [-3,2,-1]
	reverse(seq.begin(), seq.end());
	for (vector<int>::iterator it = seq.begin(); it != seq.end(); ++it)//note that I use ++it here (stackoverflow) if test wont work may try and change it to it++
	{
		*it = *it * (-1);
	}
}


//@@@@@changing to hash instead of trees@@@@@@@
void genomes::updateUniqueBlocksHash(vector<int> &seq) {
	if (seq.size() < _minUniqueBlockSize || seq.size() > _maxUniqueBlockSize)
		return;
	if (abs(seq[0]) > abs(seq.back()))
		reverseSubSeq(seq);
	updateHashAndVector(seq, uniqueBlocksVec, uniqueBlocksMap);
}

void genomes::updateUniqueRevBlocksHash(vector<int> &seq) {
	if (seq.size() < _minUniqueBlockSize || seq.size() > _maxUniqueBlockSize)
		return;
	if (abs(seq[0]) > abs(seq.back()))
		reverseSubSeq(seq);
	updateHashAndVector(seq, uniqueRevBlocksVec, uniqueRevBlocksMap);
}
void genomes::updateHashAndVector(vector<int>& seq, vector<int>& blockVec, map<vector<int>, int>& hashMap){

	if (hashMap.find(seq) != hashMap.end())
	{
		hashMap[seq]++;
		return;
	}
	hashMap[seq] = 1;
	if (blockVec.size() < seq.size())
	{
		for (size_t i = blockVec.size(); i < seq.size(); i++)
		{
			blockVec.push_back(0);
		}
	}
	blockVec[seq.size() - 1]++;
}

resultType genomes::calcFisFusBasedOnTree() {
	// tree t(Re_params::_inputTreeFileName); // replaced to reduce file handling (costly)
	tree t = _t;
	vector<GenomeClass>::iterator genomeIt = genomesVect.begin();
	resultType TreeResSet = calcFisFusAlongTree(t.getRoot(), genomeIt);
	// might consider add chromosome set and tip set sizes to the SS
	return TreeResSet;
}


genomeAndLengthPair genomes::calcUniqueBlocksAlongTree(tree::nodeP t, vector<GenomeClass>::iterator &genomeIt) { //not finished
 /* short description of algorithm
 we wish to compare unique blocks only between brother nodes (as opposed to all VS all that was previously used
 for each internal node, we look at his sons, for each son, if it's a leaf, we use its genome,
 otherwise we first call this function on it (which returns the closest leaf genome as output, as well as the distance to it)
 and use the returned genome to compare.
 To deal with multifurcating nodes, we find the closest son and compare it to the rest
 (in bifurcating nodes this is a identical to simply comparing one son to the other (NOTE: might want to check number of sons
to reduce the need to find closest son as most nodes should be bifurcating)
comparison between two genomes is done using the previously used compareBlocksGenomes method*/

// Note: One of the reviewer complained about a bias in choosing representatives for inner nodes and wish for it to be random. a new version should be implemented 24.2.22 AM

	if (t->isLeaf()) {
		// nothing to do. just return genome and its distance to the parent
		genomeAndLengthPair returnVal(&(*genomeIt), t->dis2father()); //somewhat weird way to change iterator to pointer
		genomeIt++;
		return returnVal;
	}
	assert(t->getNumberOfSons() >= 2); 
	// We assume that each node is either a leaf or has at least 2 sons. It is not biological relevant to have nodes with 1 child
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (t->getNumberOfSons() == 2) {
		genomeAndLengthPair son0 = calcUniqueBlocksAlongTree(t->getSon(0), genomeIt); //should make sure genomeIt is sent by reference and thus is always updated
		genomeAndLengthPair son1 = calcUniqueBlocksAlongTree(t->getSon(1), genomeIt);
		compareBlocksGenomes(*(son0.first), *(son1.first)); // updates the class attributes directly, thus does not need to return anything
		if (son0.second < son1.second)
		{
			son0.second += t->dis2father();
			return son0; // son0 is closer and would be used for its parent
		}
		son1.second += t->dis2father();
		return son1;
	}
	//if we got this far, we have > 2 children for this node. we need to find the closest and compare it to the rest
	int closestSon;
	double closestDist = 100000.0; //very very big compared to any branch length of trees
	genomeAndLengthPair son2return; //will hold the closest son to return upwards
	vector<GenomeClass*> sonsVector; // will hold the results for all sons so we can find the closest to compare it to the rest
	for (size_t i = 0; i < t->getNumberOfSons(); ++i) {
		genomeAndLengthPair son = calcUniqueBlocksAlongTree(t->getSon(i), genomeIt);
		sonsVector.push_back(son.first);
		if (son.second < closestDist)
		{
			son2return = son;
			closestDist = son.second;
			closestSon = i;
		}
		// found the closest son. now compare to others
	}
	for (size_t i = 0; i < sonsVector.size(); i++)
	{
		if (i == closestSon)
			continue; // We compare all sons to the closest. No need to compare it to itself
		compareBlocksGenomes(*(sonsVector[closestSon]), *(sonsVector[i]));
	}
	return son2return;
}





resultType calcFisFusAlongTree(tree::nodeP t, vector<GenomeClass>::iterator &genomeIt) {
	/* this is a recursive function that go over the tree and for each node it calculate its sons fisfus values
	(which are number of unique tips and number of probable fisfus see more info in the specific functions)
	it also construct the node sets as a union of its sons sets)
	if current node is a leaf it return its sets and (0,0) as fisfus values
	else it claculates the fisfus of it sons and return the sum of these values and its sets */
	if (t->isLeaf()) {
		pair< set <chromosomeTipType>, set <minimizedChromosomeType> > setsPair((*genomeIt).genomeTipsSet, (*genomeIt).chromosomesSet);
		set <chromosomeTipType> emptySet;
		pair<set <chromosomeTipType>, int> SSpair(emptySet, 0);
		resultType outPair(setsPair, SSpair);
		genomeIt++;
		return outPair;
	}
	set <chromosomeTipType> mergedTipsSet;
	set <minimizedChromosomeType> mergedChrSet;
	set <chromosomeTipType> uniqueTipsSet;
	int FisFusCounter = 0;
	map<int, int> tipCounter; // counting tips. used to count new tips (i.e. tips that appear only in one son)
	for (size_t i = 0; i < t->getNumberOfSons(); ++i) {
		resultType	sonRes = calcFisFusAlongTree(t->getSon(i), genomeIt); //recursive call
		//divide the res, merge sets and continue
		uniqueTipsSet.insert(sonRes.second.first.begin(), sonRes.second.first.end());
		FisFusCounter += sonRes.second.second;
		set <chromosomeTipType> sonTipsSet = sonRes.first.first;
		set <minimizedChromosomeType> sonChrSet = sonRes.first.second;

		// check chromosomes for FisFus events. we test current son against all previous ( by comparing against the merged sets)
		FisFusCounter += countChromosomesFromAFissionedInB(sonChrSet, mergedChrSet, mergedTipsSet);
		FisFusCounter += countChromosomesFromAFissionedInB(mergedChrSet, sonChrSet, sonTipsSet);
		/*NOTE: There are some cases in which we count the same event multiple times. one
		one case is when we have multiple sons (more than two) and fusion occured in the first son,
		thus for each of the follwing sons we will count a FisFus event.
		Another more important case is when a fusion occur in a node, and for each of his parent nodes we will count a FisFus event.
		Tal suggested we ignore it for now, as this is only a summary statistic. we might want to adjust the code later if results
		are not sufficient */

		for (auto elem : sonTipsSet) {
			// here we merge tip sets and count the appearances of tips for uniqueTipsCount assesment
			if (tipCounter.find(elem) != tipCounter.end())
				tipCounter[elem] ++;
			else { //haven't encounter it yet. adding to map and merged set
				tipCounter[elem] = 1;
				mergedTipsSet.insert(elem);
			}
		}

		//merge chr sets 
		for (auto elem : sonChrSet) {
			bool flag = false;
			for (auto elem2 : mergedChrSet)
			{
				if (isMinChrEquals(elem, elem2))
				{
					flag = true;
					break;
				}
			}
			if (flag)
				continue;
			mergedChrSet.insert(elem);
		}
	}
	for (auto elemPair : tipCounter)
		if (elemPair.second == 1) // is good so we count it
			uniqueTipsSet.insert(elemPair.first);
	pair< set <chromosomeTipType>, set <minimizedChromosomeType> > setsPair(mergedTipsSet, mergedChrSet);
	pair<set <chromosomeTipType>, int> SSpair(uniqueTipsSet, FisFusCounter);
	resultType outPair(setsPair, SSpair);
	return outPair;
}
bool isMinChrEquals(minimizedChromosomeType chr1, minimizedChromosomeType chr2) {
	// input: two sets of minimized chromosomes. output: true if they are equal in one of the orientations
	return (chr1.first == chr2.first && chr1.second == chr2.second) || (chr1.first == -chr2.second  && chr1.second == -chr2.first);
}
pair<int, int> calcFisfusPair(set <chromosomeTipType> genomeTipsSet1, set <minimizedChromosomeType> chromosomesSet1, set <chromosomeTipType> genomeTipsSet2, set <minimizedChromosomeType> chromosomesSet2) {
	// used to compare tip sets for pairwise comparison (somewhat more efficient than what we used before).
	// not used for now as we can't guarantee for nodes to have two sons. might be used later if current procedure is too slow
	// can be used for nodes that have two sons (most of them) after varifyication
	int max_size = genomeTipsSet1.size() + genomeTipsSet2.size();
	int xorTipsCount = 0;
	set <chromosomeTipType>::iterator it1 = genomeTipsSet1.begin();
	set <chromosomeTipType>::iterator it2 = genomeTipsSet2.begin();

	while (it1 != genomeTipsSet1.end() || it2 != genomeTipsSet2.end()) {
		// calculating xor size between tip sets
		if (it1 == genomeTipsSet1.end())
		{
			//finished tipSet1 need to finish tipSet2
			xorTipsCount++; //all other items in tipSet2 are unique
			it2++;
			continue;
		}
		if (it2 == genomeTipsSet2.end())
		{
			//finished tipSet2 need to finish tipSet1
			xorTipsCount++; //all other items in tipSet1 are unique
			it1++;
			continue;
		}
		if ((*it1) == (*it2)) {
			it1++;
			it2++;

		}
		if ((*it1) < (*it2)) {
			//tipsSet1 has a tip not in tipsSet2
			xorTipsCount++;
			it1++;
		}
		if ((*it1) > (*it2)) {
			//tipsSet2 has a tip not in tipsSet1
			xorTipsCount++;
			it2++;
		}
	}
	int FisFusCounter = 0;
	// find probable fis-fus events
	FisFusCounter += countChromosomesFromAFissionedInB(chromosomesSet1, chromosomesSet2, genomeTipsSet2);
	FisFusCounter += countChromosomesFromAFissionedInB(chromosomesSet2, chromosomesSet1, genomeTipsSet1);
	pair<int, int> outputPair(xorTipsCount, FisFusCounter);
	return outputPair;
}

int countChromosomesFromAFissionedInB(set <minimizedChromosomeType> genomeA, set <minimizedChromosomeType> genomeB, set <chromosomeTipType> chrTipsB) {
	// goes over first genome and for each chromosome look if fissioned in the other
	// checking a chromosome (in minimized format which is a pair of tips) is fissioned if both tips can be found in 
	// the tips set but the (minimized) chromosome is not in the chromosome set (in either orientation)
	int FisFusCounter = 0;
	// find probable fis-fus events
	if (genomeA.size() == 0)
		return FisFusCounter;
	set <minimizedChromosomeType>::iterator chrIt = genomeA.begin();
	while (chrIt != genomeA.end()) {
		minimizedChromosomeType chr = *chrIt;
		minimizedChromosomeType revChr = reverseChromosome(chr);
		bool firstIn = chrTipsB.find(abs(chr.first)) != chrTipsB.end();
		bool secondIn = chrTipsB.find(abs(chr.second)) != chrTipsB.end();
		bool chrIn = genomeB.find(chr) != genomeB.end();
		bool revChrIn = genomeB.find(revChr) != genomeB.end();
		if (firstIn && secondIn && !chrIn && !revChrIn)
		{
			//we first make sure both chromosome tips are in the other tip set and then check they're not together
			//then we know this chromosome underwent a FisFus event
			FisFusCounter++;
		}
		chrIt++;
	}
	return FisFusCounter;
}
minimizedChromosomeType reverseChromosome(minimizedChromosomeType chr) {
	minimizedChromosomeType revChr;
	revChr.first = -chr.second;
	revChr.second = -chr.first;
	return revChr;
}

genomes::~genomes()
{
	//destructor
}