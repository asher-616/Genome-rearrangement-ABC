#pragma once

#include <vector>
#include "genome.h"
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include "tree.h"
#include <assert.h>
#include  <math.h>


using namespace std;

typedef pair< pair< set <chromosomeTipType>, set <minimizedChromosomeType> >, pair<set <chromosomeTipType>, int>> resultType;

typedef pair< GenomeClass*, double> genomeAndLengthPair; 
// used to return closest leaf genome from each son for unique blocks along tree

class genomes
{
public:
	genomes(const vector<genomeType> &, tree sim_tree);
	genomes(string, string);
	~genomes();


	vector<int> compareBlocksGenomes(GenomeClass & genome1, GenomeClass & genome2); //why is it here? A.M
	vector<vector<double> > get_summary_stats_vector(); //get summary statistic as a vector
	void saveGenomeToFile(string filePath);

private:
	void init(const vector<genomeType>&);//used in constructors so that both constructors will share the same code
	vector<GenomeClass> genomesVect;
	double _avgMaxBlock;
	// define max and min blocks size, both normal and reversed.
	const size_t _minUniqueBlockSize = 1; //min size for unique blocks we care about
	const size_t _maxUniqueBlockSize = 10; //same but max
	// const size_t _minUniqueRevBlockSize = 1; //min size for unique reversed blocks we care about
	// const size_t _maxUniqueRevBlockSize = 100; //same but max
	
	size_t minChromosomes = 100000; //initiate with extremly big number
	size_t maxChromosomes = 0; //initiate with zero
	vector<int> all_chromosomes_vect;
	double averageChromosomes;
	double varCromosomes = 0;
	//@@@ changing to hash
	map< vector<int>, int> uniqueBlocksMap;
	map< vector<int>, int> uniqueRevBlocksMap;
	vector<int>	uniqueBlocksVec;
	vector<int>	uniqueRevBlocksVec;


	void calcSumStatBrothersOnly();
	/* calculate summary statistics. in this version we don't count unique blocks in a all Vs all fashion but we compare only 
	brother nodes. for leaves whose brother is not a leaf, I need to consult with Tal and see how to treat them. TO BE UPDATED*/

	void calcUniqueBlocksBasedOnTree(); // calculates unique blocks on tree comparing for sons each internal node.

	genomeAndLengthPair calcUniqueBlocksAlongTree(tree::nodeP t, vector<GenomeClass>::iterator &genomeIt);

	void calcSumStatAllVSAll(); // this function calcualte unique block SS in all VS all fashion (all genome pairs in tree)
	/* this is the original calcSumStat function that was used pre-Itsik. as Itsik suggestions required re-run of all data, 
	 we also implemented (and tested) a newer (and hopefully faster) version of this calculation in which we only check and count 
	 unique blocks between two leaves that share a parent (need to discuss and add how to treat cases in which one son is a leaf 
	 but the other is an internal node */

	void reverseSubSeq(vector<int>& ); //get sequence and reverse it (and its gene orientation)
	void updateUniqueBlocksTrees(vector<int>&);
	void updateUniqueRevBlocksTrees(vector<int>&);
	//void updateUniqueProbRevBlocksTrees(vector<int>);
	//void updateTree(Trie&, vector<int>&);
	bool isReversed(int, GenomeClass &);

	void updateUniqueBlocksHash(vector<int>& seq);
	void updateUniqueRevBlocksHash(vector<int>& seq);
	void updateHashAndVector(vector<int>& seq, vector<int>& blockVec, map< vector<int>, int>& hashMap);

	//summary statistics
	/* all vectors are included in their respective Trie
	double avg_max_block;
	vector<double> unique_common_blocks_sum_vec; 
	//a vector that contain for k in range _minUniqueBlockSize to _maxUniqueBlockSize sum of all unique blocks of the relevant size
	vector<double> unique_rev_common_blocks_sum_vec; //similar to the above just for reversed blocks
	vector<double> unique_prob_rev_common_blocks_sum_vec; //similar to the above just for reversed blocks
	*/

	/* fis-fus SS. these summary statistics are number of new tips (in pairwise comparisons)
	and number of apperent fis-fus event (cases in which one genome has the (a,z) chromosome
	and the other has the (a,b), (y,z) chromosomes which is a good evidence for a fis-fus event
	*/
	//pair<int, int> calcFisFusAll(); //this calculate all vs all fisfus SS
	resultType calcFisFusBasedOnTree(); //this calculate tree-based fisfus SS
	int uniqueTipCounter;
	int FisFusCounter;
	int totalTips; //size of root tip set
	int totalChrs; // size of root chr set

	tree _t;
};


pair<int, int> calcFisfusPair(set <chromosomeTipType> genomeTipsSet1, set <minimizedChromosomeType> chromosomesSet1, set <chromosomeTipType> genomeTipsSet2, set <minimizedChromosomeType> chromosomesSet2);
int countChromosomesFromAFissionedInB(set <minimizedChromosomeType> genomeA, set <minimizedChromosomeType> genomeB, set <chromosomeTipType> chrTipsB);
minimizedChromosomeType reverseChromosome(minimizedChromosomeType chr);
bool isMinChrEquals(minimizedChromosomeType chr1, minimizedChromosomeType chr2);
resultType calcFisFusAlongTree(tree::nodeP t, vector<GenomeClass>::iterator &genomeIt);

