#pragma once


#include <vector>
#include <string>
#include "tree.h"
#include "GenomeRearrangements.h"
#include "FastZip.h"
#include "RandomGenerators.h"
#include <iostream>

using namespace std;

typedef unsigned short nucID;
typedef std::vector<int> chromosomeType; //A.M
typedef std::vector<chromosomeType> genomeType; //A.M

class Simulator {
	public:
		/*explicit Simulator(vector<int> ChromosomeLength, double AParam, double InvertRate, double TranslocateRatio, double FusionRate,
			double FissionRate, double DuplicationRate, double LossRate, double randRootAparam, tree & t) :
			_rootChromosomes(ChromosomeLength), _A_param(AParam), _InvR(InvertRate), _TrR(TranslocateRatio), _FuR(FusionRate), _FiR(FissionRate),
			_DR(DuplicationRate), _LR(LossRate), zip(AParam, Re_params::_maxBlock), rootZip(randRootAparam, Re_params::_rootMaxFamilySize), _t(t) {};*/
		Simulator(const string& treePath);
		
		void initSim(vector<int> ChromosomeLength, double AParam, size_t maxBlockSize, double InvertRate, double TranslocateRatio, double FusionRate,
			double FissionRate, double DuplicationRate, double LossRate, double randRootAparam);

		void setSeed(int seedNum);


		vector<genomeType> simulateBasedOnTree();
		void test();

		
		vector<vector<int>> get_event_counter_vectors();

		tree getSimTree();
	private:
		void simulateAlongTree(tree::nodeP t,
			const genomeType& fatherSeq,
			vector<genomeType> & simSequences);

		void simualteEventsAlongAspecificBranch(const genomeType& ancestralSequence,
			double branchLength, genomeType& outputSeq);
		void simualteEventsAlongAspecificBranch_rec(const genomeType& ancestralGenome,
			double branchLength, genomeType& outputGenome); //an old recursive version

		genomeType generateRootGenome();//new root method. used for duplications (for M2. not used in M0 or M1)
		//genomeType generateRootGenomeOld(); // no duplications, each gene direction decided arbitrarily REMOVED 03.01.21
		genomeType generateRootGenomeNoDupWLOG(); // no duplications, all genes start as positive (WLOG)
		size_t _rootLength;
		vector<int> _rootChromosomes; //vector of chromosome lengths in root
		double _A_param; // this parameter dictates the Zippfian distribution of indel sizes
		double _InvR; // note that this is the invertion rate per position
		double _TrR; // note that this is the translocation rate per position
		double _FuR; //fusion rate per chromosome number (calculated as 2*(chrNum -1))
		double _FiR; //fission rate per position
		double _DR; // gene duplication rate per gene. For the future
		double _LR; // gene loss rate per gene. for the furute

		tree _t;

		FastZip zip;
		FastZip rootZip;

		vector<int> inv_counter_vec, trans_counter_vec, fis_counter_vec, fus_counter_vec; // 21.12.21 to save events for all branches
		//21.12.21 need to write a function that extract this



};
