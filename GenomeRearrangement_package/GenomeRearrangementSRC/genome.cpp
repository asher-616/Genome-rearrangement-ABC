#include "genome.h"

GenomeClass::GenomeClass(const genomeType & newgenome) {
	//geneFamilies = numberOfGeneFamilies;
	genome = newgenome;
	unsigned int genomeLength = 0; //temp for calculating genome size
	for (size_t i = 0; i < newgenome.size(); i++)
	{
		genomeLength += newgenome[i].size();
		// follwing are used for fis-fus SS. for now I assume tips are last gene only. may be changed later
		int first = newgenome[i].front();
		int last = newgenome[i].back();
		this->genomeTipsSet.insert(abs(first));
		this->genomeTipsSet.insert(abs(last));
		pair <int, int> chrMinimized(first, last);
		this->chromosomesSet.insert(chrMinimized);
	}
	genomeSize = genomeLength;
	getarrangement();




}

void GenomeClass::getarrangement() {
	//get a genome as input and returns a vector of locations for each gene
	//NOTE: this assumes that all genes are numbered 1...n. in other cases you should adjust your gene list
	
	locationType loc_initiator(2,0);
	locationVectorType genomeArrangementTemp(genomeSize, loc_initiator);
	unsigned int chromosomes = genome.size();
	for (int i = 0; i < chromosomes; i++)
	{
		unsigned int chromosome_size = genome[i].size();
		for (size_t j = 0; j < chromosome_size; j++)
		{
			/*int signedlocation = j + 1; //location starts at 1 to allow direction
			if (genome[i][j] < 0)
			{
				signedlocation *= -1;
			}
			locationType location{ i,signedlocation };
			genomeArrangementTemp[abs(genome[i][j]) - 1].push_back(location);
			*/

			/* old code- used before introducing duplication and loss */
			genomeArrangementTemp[abs(genome[i][j])-1][0] = i;
			int signedlocation = j + 1; //location starts at 1 to allow direction
			if (genome[i][j]<0)
			{
				signedlocation = -signedlocation; 
			}
			genomeArrangementTemp[abs(genome[i][j])-1][1] = signedlocation;
		}
	}
	genomeArrangement = genomeArrangementTemp;
	return;

}


int GenomeClass::nextGene(int gene) {
	//gets a genome and a gene number and return the next gene for the gene orientation
	//(forward if both same sign, backword otherwise)
	//used for compareBlocksGenomes
	int gCh = genomeArrangement[abs(gene)-1][0]; //gene chromosome
	int gLoc = genomeArrangement[abs(gene)-1][1]; //gene location (on chromosome)
	if (gene * gLoc > 0)
	{//going forward since both genes are in same orientation
		if (abs(gLoc) < genome[gCh].size()) //to avoid going over the edge
		{
			return genome[gCh][abs(gLoc)]; //remember that indeces in genomeArrangement start at 1
		}
		else //got to edge of chromosome
		{
			return 0; 
		}
	}
	else
	{
		if (abs(gLoc) > 1)
		{
			return genome[gCh][abs(gLoc) - 2];
		}
		else //got to edge of chromosome
		{
			return 0;
		}
	}
	

}


int GenomeClass::prevGene(int gene) {
	//gets a genome and a gene number and return the next gene for the gene orientation
	//(forward if both same sign, backword otherwise)
	//used for compareBlocksGenomes
	int index = abs(gene) - 1;
	int gCh = genomeArrangement[index][0]; //gene chromosome
	int gLoc = genomeArrangement[index][1]; //gene location (on chromosome)
	if (gene * gLoc < 0)
	{//going forward since genes are in opposite orientation
		if (abs(gLoc) < genome[gCh].size()) //to avoid going over the edge
		{
			return genome[gCh][abs(gLoc) ];
		}
		else //got to edge of chromosome
		{
			return 0;
		}
	}
	else
	{
		if (abs(gLoc) > 1)
		{
			return genome[gCh][abs(gLoc) - 2];
		}
		else //got to edge of chromosome
		{
			return 0;
		}
	}
}















GenomeClass::~GenomeClass(){ //destructor
}