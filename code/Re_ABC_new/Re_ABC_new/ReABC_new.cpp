// new main 20.12.21
// We are aiming to remove the ABC part from the C++ code to python and clean the C++ as much as possible
// Updated by: Asher Moshe & Elya Wygoda
#include "ReABC_new.h"

vector<vector<double> > perform_simulation(vector<int> ChromosomeLengths, double a_param, double inv_rate, double trans_rate,
	double fus_rate, double fis_rate, double dup_rate, double loss_rate, double root_a_param) {
	
}





// ######################################################################
// ######                           Main                           ######
// ######################################################################
int main(int argc, const char * argv[]) {
	srand((unsigned)(time(0))); // for random number generation. Do not delete.
	// We should get all relevant parameters from the user (python code) and just run a single simulation here
	
	// Read parameters into parameter object (need to discuss with Elya how we receive those)

	// We need two modes: 
	// 1. Simulation mode, in which we get parameters (and tree file) and create a simulation. 
	//	  We need to return either genome data (in our format), summary statistics, or both.
	// 2. Summary statistics extraction mode, we receive genomes and tree file and return the summary statistics.


	// If mode 1 
	// Simulate genomes

	// Extract summary statistics (unless not requested)
	// If genome data is requested, create it now.

	// If mode 2
	// Read input genomes
	// Extract summary statistics

	// Return Summary statistics (unless mode 1 and not requested).

}