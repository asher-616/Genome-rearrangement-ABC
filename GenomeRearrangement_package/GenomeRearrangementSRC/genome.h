#pragma once


#include <vector>
#include <cstdlib>
#include <iostream>
#include <set>
using namespace std;


typedef vector<int> chromosomeType; //A.M
typedef vector<chromosomeType> genomeType; //A.M
typedef vector<int> locationType;
typedef vector<locationType> locationVectorType;
//typedef vector<locationVectorType> genomearrangeType; //is used when there option for gene families of size >1
typedef int chromosomeTipType; //this may be changed if we decide that a tip is more than one gene
typedef pair <chromosomeTipType, chromosomeTipType> minimizedChromosomeType;


class GenomeClass
{
public:
	GenomeClass(const genomeType &);
	~GenomeClass();
	genomeType genome;
	locationVectorType genomeArrangement;
	//int geneFamilies;
	int nextGene(int gene);
	int prevGene(int gene);


	//new data structures for fis-fus summary statistics
	set <chromosomeTipType> genomeTipsSet;
	set <minimizedChromosomeType> chromosomesSet;


private:
	void getarrangement();
	unsigned int genomeSize;
};






