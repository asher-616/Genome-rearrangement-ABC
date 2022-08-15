#include "RandomGenerators.h"
// was danaRandomGenerators.h

using namespace std;

int RandomGenerators::seed;
default_random_engine RandomGenerators::generator;
mt19937 RandomGenerators::mt_rand;


void RandomGenerators::initRandomGenerator(int seed){
	if (seed < 0){
		int chrono_seed = static_cast<int>(chrono::system_clock::now().time_since_epoch().count());
		generator.seed(chrono_seed);
		mt_rand.seed(chrono_seed);
	} else {
		generator.seed(seed);
		mt_rand.seed(seed);
	}
}





double RandomGenerators::drawExp(double lambda) {
	
	exponential_distribution<double> distribution(lambda);
	double number = distribution(generator);
	return number;
}

double RandomGenerators::uniform() { // uniform between 0 and 1
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return number;
}

int RandomGenerators::uniform(int a, int b) { // a random number between a and b, including a and b.
	uniform_int_distribution<int> distribution(a, b);
	int number = distribution(generator);
	return number;
}


RandomGenerators::~RandomGenerators()
{
}
