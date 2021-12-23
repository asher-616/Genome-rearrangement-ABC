#include "FastZip.h"





FastZip::FastZip(double aParam, int max)
{
	vector<int> valueVector;
	vector<double> ratesVector;
	double sumRates = 0;
	for (size_t i = 1; i <= max; i++)
	{
		valueVector.push_back(i);
		ratesVector.push_back(1 / pow(i, aParam));
		sumRates += ratesVector.back();
	}
	for (size_t i = 0; i < ratesVector.size(); i++)
	{
		ratesVector[i] /= sumRates;
		//cout << ratesVector[i]<< "\t";
	}
	cout << endl;
	for (size_t i = 0; i < ratesVector.size(); i++)
	{
		ratesVector[i] *= ratesVector.size();
	}
	size_t j = 0;
	for (; j < ratesVector.size(); j++)
		if (ratesVector[j] < 1) break;
	vector<int> lowVals;
	vector<double> lowRates;
	vector<int> hiVals;
	vector<double> hiRates;
	

	for (size_t k = 0; k < j; k++)
	{
		hiVals.push_back(valueVector[k]);
		hiRates.push_back(ratesVector[k]);
	}
	for (size_t k = j; k < ratesVector.size(); k++)
	{
		lowVals.push_back(valueVector[k]);
		lowRates.push_back(ratesVector[k]);
	}


	reverse(hiVals.begin(), hiVals.end()); //so I can pull values high to low using pop_back
	reverse(hiRates.begin(), hiRates.end());

	int lowVal = lowVals.back();
	lowVals.pop_back();
	double lowRate = lowRates.back();
	lowRates.pop_back();

	int hiVal = hiVals.back();
	hiVals.pop_back();
	double hiRate = hiRates.back();
	hiRates.pop_back();

	while (lowRates.size() + hiRates.size() > 0) {
		lowVector.push_back(lowVal);
		highVector.push_back(hiVal);
		lowHighRateVector.push_back(lowRate);
		hiRate -= 1 - lowRate;
		if (hiRate < 1)
		{
			lowRate = hiRate;
			lowVal = hiVal;
			hiVal = hiVals.back();
			hiVals.pop_back();
			hiRate = hiRates.back();
			hiRates.pop_back();
		}
		else
		{
			lowVal = lowVals.back();
			lowVals.pop_back();
			lowRate = lowRates.back();
			lowRates.pop_back();
		}
	}
	lowVector.push_back(lowVal);
	highVector.push_back(hiVal);
	lowHighRateVector.push_back(lowRate);
	lowVector.push_back(hiVal);
	highVector.push_back(hiVal);
	lowHighRateVector.push_back(1.0);
	numBins = highVector.size();

}

int FastZip::drawZip() {
	int bin = uniform(0, numBins - 1);
	double randomNum = uniform();
	if (randomNum < lowHighRateVector[bin])
		return lowVector[bin];
	else
		return highVector[bin];
}

FastZip::~FastZip()
{
}


/*
int test() {
	double a_param = 1.3;
	int max_block = 10;
	FastZip zip(a_param, max_block);

	int tries = 1000000;
	vector<int> counts(max_block, 0);
	for (size_t i = 0; i < tries; i++)
	{
		int rand = zip.drawZip();
		counts[rand - 1]++;
	}
	for (size_t i = 0; i < max_block; i++)
	{
		cout << 1.0 * counts[i] / tries << "\t"; 
	}
	return 0;
}

*/