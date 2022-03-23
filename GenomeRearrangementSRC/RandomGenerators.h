#pragma once
// was danaRandomGenerators.h
#include <algorithm>
#include <cmath>
#include <chrono>
#include <random>
#include <cassert>
using namespace std;


class RandomGenerators
{
private:
	static int seed;
	static default_random_engine generator;
	static mt19937 mt_rand;

public:

	static void initRandomGenerator(int seed);

	static double drawExp(double lambda);
	static double uniform();
	static int uniform(int a, int b);

	virtual ~RandomGenerators();
};



