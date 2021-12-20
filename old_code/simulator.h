#ifndef _SIMULATOR_H
#define _SIMULATOR_H

#include <vector>
#include <string>
#include "tree.h"
#include "GenomeRearrangements.h"
#include "FastZip.h"
#include "danaRandomGenerators.h"
#include "GRABC_options.h"

using namespace std;

typedef unsigned short nucID;
typedef std::vector<int> chromosomeType; //A.M
typedef std::vector<chromosomeType> genomeType; //A.M

class Simulator {
	public:
		explicit Simulator(vector<int> ChromosomeLength, double AParam, double InvertRate, double TranslocateRatio, double FusionRate,
			double FissionRate, double DuplicationRate, double LossRate, double randRootAparam) :
			_rootChromosomes(ChromosomeLength), _A_param(AParam), _InvR(InvertRate), _TrR(TranslocateRatio), _FuR(FusionRate), _FiR(FissionRate),
			_DR(DuplicationRate), _LR(LossRate), zip(AParam, GRABC_options::_maxBlock), rootZip(randRootAparam,GRABC_options::_rootMaxFamilySize) {};
		vector<genomeType> simulateBasedOnTree(string & treeFileName);
		void test();
	private:
		void simulateAlongTree(tree::nodeP t,
			const genomeType& fatherSeq,
			vector<genomeType> & simSequences);

		void simualteEventsAlongAspecificBranch(const genomeType& ancestralSequence,
			double branchLength, genomeType& outputSeq);
		void simualteEventsAlongAspecificBranch_rec(const genomeType& ancestralGenome,
			double branchLength, genomeType& outputGenome); //an old recursive version

		genomeType generateRootGenome();//new root method. used for duplications (for M2. not used in M0 or M1)
		genomeType generateRootGenomeOld(); // no duplications, each gene direction decided arbitrarily
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

		FastZip zip;
		FastZip rootZip;



};
// The simulator simulates sequences along a tree without substitutions
// The output is a list of strings, each corresponding to a single sequence
// The name of each sequence is not given in this output.
// For example, if we simulate along the tree ((s1:0.3; s2: 0.4), s3:0.5);
// We will get as output a vector of string, say, V, so that:
// v[0] = "ACGAAG-"
// v[1] = "A---AG-"
// v[2] = "-CGAAGA"
// In truth, we do not simulate substitutions, so that vectors will actually look like:
// v[0] = "AAAAAA-"
// v[1] = "A---AA-"
// v[2] = "-AAAAAA"
#endif