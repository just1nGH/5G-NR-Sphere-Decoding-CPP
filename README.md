# 5G-NR-Sphere-Decoding-CPP

C++ Implenmentation of 5G NR MIMO Sphere Decoder
- Tree traverse stratigies: single tree search 
- Support  different modulation schemes for Tx antennas

The implementation sphere decoding algorithms(single tree seach)for seeking the maximum-likelihood solution 
for a set of symbols transmitted over the MIMO channel. 


Note: the algorithm relys on QR decompostion which performed by using Eigen package, in order to use the code, the package must be installed.


```
//-----------------------------------------------------------------------------------------------
// softbits = nrSphereDecoder(H, rxSymbs) uses sphere decoding algorithms (single tree seach)
// for seeking the maximum - likelihood solution for a set of symbols transmitted over the MIMO channel.
// ***Input***
//	* H - complex channel matrix Nr-by-Nt, Nr is the number of rx antennas, Nt is the number of Tx antennas
//	* rxSymbs - complex received symbol vector  Nr-by-1
// ***Output*** softbits(llr) - soft bits (not sacled by 1/N0)
% Author: Dr J Mao
% Email: juquan.justin.mao@gmail.com
% 2021 Nov
%--------------------------------------------------------------------------
```

### Example run simlation results
Start by try `example_run.m`
```
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
```
