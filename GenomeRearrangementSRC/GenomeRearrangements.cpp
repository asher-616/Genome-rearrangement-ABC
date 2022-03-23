#include "GenomeRearrangements.h"


int inv_counter, trans_counter, fis_counter, fus_counter;

void drawRandomLocationOld(vector<int> & location, const vector<vector<int> > & genome, const int genomelength) {
	//receives a genome (and length) and returns a randomly chosen (uniform distribution) gene location (vector of chromosome number and gene number)
	unsigned int GeneNum = RandomGenerators::uniform(0,genomelength); //need to draw a random number in range (1-genomelength)
	unsigned int chromosome = 0;
	location.clear();
	while (GeneNum > genome[chromosome].size())
	{
		GeneNum -= genome[chromosome].size();
		chromosome++;
	}
	if (chromosome < genome.size() - 1 && GeneNum == genome[chromosome].size())
	{
		double someUniformRandomNum = RandomGenerators::uniform(); //random in range [0,1]
		if (someUniformRandomNum < 0.5)
		{
			chromosome++;
			GeneNum = 0;
		}
	}
	location.push_back(chromosome);
	location.push_back(GeneNum);
}

 bool drawRandomLocation(vector<int> & location, const vector<vector<int> > & genome, const int genomelength) {
	 // input: location vector - to update with drawn location,
	 // genome vector - for reference to depict exact location (chromosome and position within)
	 // genomelength - to draw global position in the genome (save the need to figure genome size)
	 // process: draw global position uniformly, find exact position, update it into location,
	 // draw event direction - if location is one of the genome ends, direction is determined inwards (only available choice)
	 // otherwise direction is determined randomly with 1:1 ratio
	 // output: true if direction is forward, false if direction is backwards (direction is relevant for inversion and transposition only)
	 unsigned int GeneNum = RandomGenerators::uniform(0, genomelength - 1); //need to draw a random number in range (0-genomelength-1)
	 unsigned int chromosome = 0;
	 location.clear();
	 if (GeneNum == 0) //first gene in first chromosome
	 {
		 // genome beginning 
		 location.push_back(0);
		 location.push_back(0);
		 return true;
	 }
	 if (GeneNum == genomelength - 1) // last gene in last chromosome
	 {
		 location.push_back(genome.size() - 1); //last chromosome
		 location.push_back(genome[genome.size() - 1].size() - 1); //last gene in last chromosome
		 return false;
	 }
	 while (GeneNum > genome[chromosome].size())
	 {
		 GeneNum -= genome[chromosome].size();
		 chromosome++;
	 }
	 bool forward = RandomGenerators::uniform() > 0.5;
	 if ((GeneNum == genome[chromosome].size()))
	 {
		 GeneNum = 0;
		 chromosome++;
	 }
	 location.push_back(chromosome);
	 location.push_back(GeneNum);
	 return forward;
 }

void SimulateEvent(vector<vector<int> > & genome, int genomeLength, double transpositionRate, double inversionRate, FastZip blockLengthparameter,
					double fusionRate, double fissionRate, double duplicationRate, double deletionRate) {
	/* NOTE: this function occur only AFTER an event occured. 
	This function receives a genome (vector of vector of int), and rates for all possible events and an a-parameter to choose block size (for inversion and tranposition).
	the function chooses block size according to a distribution (zipfian for now) and decides which of the events occured according to their rates.
	for the inversion event, a single break position is drawn from a distribution (uniform)
	for a tranposition, two break points are drawn (one for the original block position and another for its destination)
	for fusion, two chromosomes are chosen
	for fission a single break point is chosen
	these parameters are sent to the relevant function to perform the event
	*/
	for (size_t i = 0; i < genome.size(); i++)
	{
		if (genome[i].size() == 0){
			cout << "BUG!!!!!!!" << endl;
			exit(0);
		}
	}
	double someUniformRandomNum = RandomGenerators::uniform(); //random in range [0,1] 
	vector<int> location;
	/* //old code for location drawing (pre-Itsik). no direction was drawn (location was starting point)
	do
	{
		//cout << "draw location ";
		drawRandomLocation(location, genome, genomeLength);
		//cout << location[0] << '\t' << location[1] << endl;
	} while (location[1] == genome[location[0]].size());
	*/
	bool forward = drawRandomLocation(location, genome, genomeLength);
	//cout <<"start location "<< location[0] << '\t' << location[1] << "\t" << forward << endl;
	unsigned int blockSize;
	/* //old code for block-size (pre-Itsik) not sure why comparing to genomeLength was needed
	do
	{
		blockSize = blockLengthparameter.drawZip();
		//draw from zip distribution. Note max parameter is set to 50 (arbitrarily)
		//cout << "block size " << blockSize << endl;
		//cout << "genome size " << genomeLength << endl;
	} while (blockSize > genome[location[0]].size() - location[1] || blockSize >= genomeLength);
	/**/
	bool forwardOK = false;
	bool backwardOK = false;
	//int counter = 0; // for testing
	while (!(forwardOK) || !(backwardOK))
	{
		blockSize = blockLengthparameter.drawZip();
		//draw from zip distribution. Note max parameter is set to 50 (arbitrarily)
		//cout << "block size " << blockSize << endl;
		//cout << "genome size " << genomeLength << endl;
		/* used for testing that we don't have an event that go before the start 
		if (!forward)
			cout << "start " << location[1] << " size " << blockSize << endl;
		*/
		forwardOK = (!forward) || (blockSize + location[1] <= genome[location[0]].size());
		backwardOK = (forward) || (location[1] >= blockSize - 1);
		/*
		cout << !forwardOK << !backwardOK << endl;
		counter++;
		if (counter > 100)
			cout << "bug";
		*/
	}
	//cout << "block size " << blockSize <<" ch " << genome[location[0]].size() << endl;
	/* //old. had problem with backword size 
	do
	{
		blockSize = blockLengthparameter.drawZip();
		//draw from zip distribution. Note max parameter is set to 50 (arbitrarily)
		//cout << "block size " << blockSize << endl;
		//cout << "genome size " << genomeLength << endl;
		if (!forward)
			cout << "start " << location[1] << " size " << blockSize << endl;
		bool forwardOK = (!forward) || (blockSize + location[1] <= genome[location[0]].size());
		bool backwardOK = (forward) || (location[1] >= blockSize);
	} while (((blockSize + location[1] > genome[location[0]].size()) && forward)||( ( location[1] - blockSize < 0) && !forward ));
	/**/
	// checking that the block doesn't transpass the relevant edge (in the relevant direction)

	double totalRates = transpositionRate + inversionRate + fusionRate + fissionRate + duplicationRate + deletionRate;
	if (someUniformRandomNum < inversionRate/totalRates)
	{
		if (!forward)
			location[1] -= (blockSize - 1) ; //change location to event start point. the -1 explained in the transposition event
		inv_counter++;

		Inversion(genome, location, blockSize);
		//cout << "inversion " << location[0] << " " << location[1] << endl;
	}
	else if (someUniformRandomNum < (inversionRate + transpositionRate) / totalRates)

	{
		//cout << "transposition ";
		if (!forward)
			location[1] -= (blockSize - 1); //change location to event start point. the -1 is to include the start point
		// example: if we have a backword block of size 2 in position 3, 3-2 = 1 and will result in block of position 1,2 
		// and not the wanted 2,3 block
		if (location[1] <0)
		{
			cout << "ERROR" << endl;
		}
		vector<int> destination;
		do
		{
			drawRandomLocation(destination, genome, genomeLength);
			//cout << "destination " << destination[0] <<" "<< destination[1] << endl; //for testing
		} while (destination[0] == location[0] && destination[1] >= location[1] && destination[1]<= location[1]+blockSize);
		//cout << "destination " << destination[0] << " " << destination[1] << "ch size " << genome[destination[0]].size() << endl;
		trans_counter++;
		Transposition(genome, location, blockSize, destination);
		//cout << "xxxxx" << endl;
	}
	else if (someUniformRandomNum < (inversionRate + transpositionRate + fusionRate) / totalRates)
	{
		//do fusion
		int chromosome1 = RandomGenerators::uniform(0, genome.size()-1);
		int chromosome2;
		do
		{
			chromosome2 = RandomGenerators::uniform(0, genome.size()-1);
		} while (chromosome1 == chromosome2);
		//cout << "fusion " << chromosome1 << "+ " << chromosome2 << endl;
		fus_counter++;
		Fusion(genome, chromosome1, chromosome2);
	}
	else if (someUniformRandomNum < (inversionRate + transpositionRate + fusionRate + fissionRate) / totalRates) //for when duplications and deletions are added
	{
		//do fission
		fis_counter++;
		Fission(genome, genomeLength);
	}
	
	else if (someUniformRandomNum < (inversionRate + transpositionRate + fusionRate + fissionRate + duplicationRate) / totalRates)
	{
		//do duplication
	}
	else
	{
		//do deletion
	}
	
}





void Inversion(vector<vector<int> > & genome, vector<int> & eventLocation, int eventSize){
	//function receives a genome, event location (a vector of chromosome and location) and block size (number of genes to be inverted) invert genes while changing annotation
	//eventLocation contain two int: (chromosomeNumber,genePosition)
	//cout << "invesion" << endl;
	std::list<int> EventBlock;
	//cout << "inversion\n" << endl;
	for (size_t i = eventLocation[1]; i < eventLocation[1] + eventSize; i++){
		//copy block to a list
		EventBlock.push_front(genome[eventLocation[0]][i]*(-1));
	}
	int i = 0;
	for (list<int> :: iterator it = EventBlock.begin(); it != EventBlock.end(); it++)
	{
		genome[eventLocation[0]][eventLocation[1] + i] = *it;
		i++;
	}
	//cout << "out of inversion" << endl;
}





void Transposition(vector<vector<int> > & genome, vector<int> eventStartPoint, int eventSize, vector<int> eventdestination){
	//takes a block (of eventSize genes beginning at eventStartPoint and moves it to eventDestination
	double someUniformRandomNum = RandomGenerators::uniform(); //random in range [0,1] 
	//cout << "transposition" << endl;
	if (someUniformRandomNum < 0.5) // used to allow block to be reinserted in both orientations (1:1 ratio).
	{
		Inversion(genome, eventStartPoint, eventSize);
	}
	std::vector<int>::iterator itStart, itDest;
	//itStart = genome[eventStartPoint[0]].begin();
	//itDest = genome[eventdestination[0]].begin();
	//pushing block into destination 

	genome[eventdestination[0]].insert(genome[eventdestination[0]].begin() + eventdestination[1], genome[eventStartPoint[0]].begin() + eventStartPoint[1], genome[eventStartPoint[0]].begin() + eventStartPoint[1] + eventSize);
	//removing block from old position
	if (eventdestination[0] == eventStartPoint[0] && eventdestination[1] < eventStartPoint[1])
	{
		eventStartPoint[1] += eventSize;
	}
	genome[eventStartPoint[0]].erase(genome[eventStartPoint[0]].begin() + eventStartPoint[1], genome[eventStartPoint[0]].begin() + eventStartPoint[1] + eventSize );
	if (genome[eventStartPoint[0]].size() == 0 )
	{
		genome.erase(genome.begin() + eventStartPoint[0]);
	}
}






void Fusion(vector<vector<int>>& genome, int chromosome1, int chromosome2)
{
	//get a genome and two (different) chromosome numbers and merge them into one by choosing what orientation the second is attached to the first
	//cout << "fusion" << endl;
	double someUniformRandomNum = RandomGenerators::uniform(); //random in range [0,1] 
	if (someUniformRandomNum < 0.5) // used to allow block to be reinserted in both orientations (1:1 ratio).
	{
		reverse(genome[chromosome2].begin(), genome[chromosome2].end());
		for (size_t i = 0; i < genome[chromosome2].size(); i++)
		{
			genome[chromosome2][i] *= -1;
		}
	}
	someUniformRandomNum = RandomGenerators::uniform(); //random in range [0,1] 
	if (someUniformRandomNum < 0.5) // used to allow block to be reinserted in both orientations (1:1 ratio).
	{
		reverse(genome[chromosome1].begin(), genome[chromosome1].end());
		for (size_t i = 0; i < genome[chromosome1].size(); i++)
		{
			genome[chromosome1][i] *= -1;
		}
	}
	genome[chromosome1].insert(genome[chromosome1].end(), genome[chromosome2].begin(), genome[chromosome2].end());
	genome.erase(genome.begin() + chromosome2);
}






void Fission(vector<vector<int>>& genome, int genomeLength)
{
	vector<int> eventLocation;
	do
	{
		drawRandomLocation(eventLocation, genome, genomeLength);
	} while (eventLocation[1] == 0 || eventLocation[1] >= genome[eventLocation[0]].size() - 1);
	//cout << "fission " << eventLocation[1] <<' ' << genome[eventLocation[0]].size() - eventLocation[1] << endl;
	int chromosomeToSplit = eventLocation[0];
	//vector<int>::const_iterator breakPoint = genome[chromosomeToSplit].begin() + eventLocation[1];
	//vector<int>::const_iterator endPoint = genome[chromosomeToSplit].end();
	vector<int>::iterator breakPoint = genome[chromosomeToSplit].begin() + eventLocation[1];
	vector<int>::iterator endPoint = genome[chromosomeToSplit].end();
	vector<int> newChromosome(breakPoint, endPoint );
	genome[chromosomeToSplit].erase(breakPoint, endPoint);
	genome.push_back(newChromosome);
}

void GeneDuplication(vector<vector<int>>& genome, vector<int> eventStartPoint, vector<int> eventdestination)
{
	vector<int>::iterator it = genome[eventdestination[0]].begin();
	int gene = genome[eventStartPoint[0]][eventStartPoint[1]];
	double someUniformRandomNum = RandomGenerators::uniform(); //random in range [0,1] 
	if (someUniformRandomNum < 0.5)
	{
		gene *= -1;
	}
	genome[eventdestination[0]].insert(it + eventdestination[1], gene);
}

void GeneLoss(vector<vector<int>>& genome, vector<int> eventStartPoint)
{
	vector<int>::iterator it = genome[eventStartPoint[0]].begin();
	genome[eventStartPoint[0]].erase(it + eventStartPoint[1]);
}







/* Cancelled for now
void reciprocalTranslocation(vector<vector<int> > & genome,int chromosomeBlock1, int block1Start, int block1End, int chromosomeBlock2, int block2Start, int block2End) {
	//recieve a genome and 2 blocks (represented by 2 ints: chromosome number, block start and end) and replaces them
	//NOTICE: start to end switch isn't handled yet. do we want to take centromers and telomers into account?
	vector<int>::const_iterator first1 = genome[chromosomeBlock1].begin() + block1Start;
	vector<int>::const_iterator last1 = genome[chromosomeBlock1].begin() + block1End;
	vector<int> block1(first1, last1);
	vector<int>::const_iterator first2 = genome[chromosomeBlock2].begin() + block2Start;
	vector<int>::const_iterator last2 = genome[chromosomeBlock2].begin() + block2End;
	vector<int> block2(first2, last2);
	genome[chromosomeBlock1].erase(first1, last1);
	genome[chromosomeBlock1].erase(first2, last2);
	genome[chromosomeBlock1].insert(first1, block2.begin(), block2.end());
	genome[chromosomeBlock1].insert(first2, block1.begin(), block1.end());
}
*/
/* used to test rearrangement events
int main() {
	vector<vector<int> > genome;
	int chromosomeNum = 3;
	int numOfGenes = 20;
	int GeneNum = 1; //used for continuity of gene numbers throughout the genome
	for (size_t i = 0; i < chromosomeNum; i++)
	{
		vector<int> tempChromosome;
		for (int i = 0; i < numOfGenes; i++)
		{
			double orient = RandomGenerators::uniform();
			if (orient < 0.5)
				tempChromosome.push_back(GeneNum*(-1));
			else
				tempChromosome.push_back(GeneNum);
			GeneNum++;
		}
		genome.push_back(tempChromosome);
	}
	cout << "original genome" << endl;
	for (size_t i = 0; i < genome.size(); i++)
	{
		for (std::vector<int>::const_iterator j = genome[i].begin(); j != genome[i].end(); ++j)
			std::cout << *j << ' ';
		cout << endl;
	}
	SimulateEvent(genome, GeneNum, 0.5, 0.5, 1);
	for (size_t i = 0; i < genome.size(); i++)
	{
		for (std::vector<int>::const_iterator j = genome[i].begin(); j != genome[i].end(); ++j)
			std::cout << *j << ' ';
		cout << endl;
	}

}
*/
void printGenome(vector<vector<int> > genome) { //used to print genome to screen, for testing purposes
	for (size_t i = 0; i < genome.size(); i++)
	{
		for (std::vector<int>::const_iterator j = genome[i].begin(); j != genome[i].end(); ++j)
			std::cout << *j << ' ';
		cout << endl;
	}
}
