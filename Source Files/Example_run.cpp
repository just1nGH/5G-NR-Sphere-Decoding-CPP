// sphereDecoder.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "nrSphereDecorder.h"
#include <random>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <Eigen/dense>
#include <Eigen/QR>
using namespace std;
using chrono::system_clock;

void example_run();
int main()
{
	std::cout << "Hello World!\n";
	example_run();
}

void example_run()
{
	double N0 = 0.1; //noise
	int nTxAnts = 3;
	int nRxAnts = 4;

	vector<string> moduTypes{ "16qam","16qam","64qam" };

	SphereDecoder sd(nRxAnts, nTxAnts, moduTypes);

	// random engines
	default_random_engine random_engine;
	bernoulli_distribution  bern_dist;
	normal_distribution<double> norm_dist(0, 1);

	int K = sd.getTotNumBits();
	Veci msg; msg.reserve(K);// generate msg
	for (unsigned j = 0; j < K; j++)
		msg.push_back(bern_dist(random_engine));

	Veci offsets = sd.getOffsets();
	Veci Ks = sd.getKs();

	ComplexVec txSymbs; txSymbs.reserve(nTxAnts);
	Veci bitlabels;
	for (auto i = 0; i < nTxAnts; i++) {
		bitlabels = Veci(msg.begin() + offsets[i], msg.begin() + offsets[i] + Ks[i]);
		txSymbs.push_back(nrModuMapper(bitlabels, moduTypes[i])[0]);
	}

	int nRels = 10;
	int nErrsSoft = 0;

	for (int irel = 0; irel < nRels; irel++) {
		//flat channel
		random_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
		ComplexMatrix2D H(nRxAnts);
		for (auto i = 0; i < nRxAnts; i++) {
			H[i].reserve(nTxAnts);
			for (auto j = 0; j < nTxAnts; j++) {
				H[i].push_back(1.0 / sqrt(2.0 * nTxAnts) * complex<double>(norm_dist(random_engine), norm_dist(random_engine)));
			}
		}

		// rxSymbls
		ComplexVec rxSymbs; rxSymbs.reserve(nRxAnts);
		for (auto i = 0; i < nRxAnts; i++) {
			rxSymbs.push_back(inner_product(H[i].begin(), H[i].end(), txSymbs.begin(), complex<double>(0.0)));
		}

		//noise
		random_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
		for (auto& e : rxSymbs)
			e += sqrt(N0 / 2) * complex<double>(norm_dist(random_engine), norm_dist(random_engine));

		// decoded bits
		Vecd softBits = sd(H, rxSymbs);
		Veci decBits; decBits.reserve(K);
		for (auto i = 0; i < K; i++) decBits.push_back((softBits[i] < 0) ? 1 : 0);
		for (auto i = 0; i < K; i++) nErrsSoft += (decBits[i] == msg[i]) ? 0 : 1;
	}
	cout << "nErrsSoft = " << nErrsSoft << " of " << K * nRels << endl;
}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started:
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file