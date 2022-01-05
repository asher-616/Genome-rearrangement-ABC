// new main 20.12.21
// We are aiming to remove the ABC part from the C++ code to python and clean the C++ as much as possible
// Updated by: Asher Moshe & Elya Wygoda
#include "ReABC_new.h"

vector<vector<double> > perform_simulation(vector<int> ChromosomeLengths) {

	Simulator sim(ChromosomeLengths, Re_params::_a_param, Re_params::_inv_rate, Re_params::_trans_rate, Re_params::_fus_rate,
		Re_params::_fis_rate, Re_params::_dup_rate, Re_params::_loss_rate, Re_params::_root_a_param, Re_params::_inputTree); // create simulation object

	vector<genomeType> simulatedGenomes = sim.simulateBasedOnTree();
	// here I need to copy the event counter vectors from the sim object
	vector<vector<int>> counters_vec = sim.get_event_counter_vectors(); // order of counter vectors is inv, trans, fis, fus
	// what next? need to return it to the python
	
}
// pipeline scheme: (03.01.22)
// step 1 - set parameters into param object
// Can be set using param object. param object is initialized with a tree and assigned parameters using init_params 
// (overloaded for both M0 and M1 parameters)
// IMPORTANT: currently the parameters are static inside the param class.
// Can this work along with the python or should we change it?

// step 2 - simulate and return simulation and counters

// step 3 - calculate summary statistics



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