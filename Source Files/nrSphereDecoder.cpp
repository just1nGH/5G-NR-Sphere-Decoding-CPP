//--------------------------------------------------------------------------------------------------------------
//                                       ** sphere decoder **
//						Author : Dr J Mao  Email : juquan.justin.mao@gmail.com  2021 Dec
//          source paper :  "Soft-Output Sphere Decoding: Performance and Implementation Aspects"
//-------------------------------------------------------------------------------------------------------------
#include <iostream>
#include "nrSphereDecorder.h"
#include <numeric>
#include <Algorithm>
#include <array>
#include <Eigen/dense>
#include <Eigen/QR>

using namespace std;

SphereDecoder::SphereDecoder(int nRxAnts, int nTxAnts, vector<string>& moduTypes, double searchRadius)
{
	NumRxAnts = nRxAnts; // number of rx antennas
	NumTxAnts = nTxAnts; // number of transmit antennas

	// Search Radius controls tree pruning for counter hypothesis, 0: hard decision,inf: max-log-map
	// values in between is the tradeoff between complexityand performance, 0.2 is a thrumb of rule value.
	SearchRadius = searchRadius;

	//modulation handling
	if (moduTypes.size() != NumTxAnts) {
		if (moduTypes.size() == 1) {
			//string tmp = moduTypes[0];
			moduTypes.insert(moduTypes.end(), NumTxAnts - 1, moduTypes[0]);
		}
		else {
			std::cerr << "The number of modulations must be 1 or equal to number of Tx antennas!" << std::endl;
			abort();
		}
	}

	// constellations for each Tx Ants
	Cstltns.reserve(NumTxAnts);
	for (auto i = 0; i < NumTxAnts; i++) Cstltns.push_back(Constellation{ moduTypes[i] });

	// record number of bits per decoded symbol
	Ks.reserve(nTxAnts);
	for (auto i = 0; i < nTxAnts; i++) Ks.push_back(Cstltns[i].K);

	// total number of bits corresponding the decoded symbols
	Offsets.reserve(nTxAnts); // record the starting position for bit labels of each decoded symbols
	NumTotBits = 0;
	for (auto i = 0; i < nTxAnts; i++) {
		Offsets.push_back(NumTotBits);
		NumTotBits += Ks[i];
	}
}
SphereDecoder::~SphereDecoder()
{
}

void SphereDecoder::nodeInit(TreeNode& node)
{
	// initilize to root node;
	// add the first symbol of the root costellation
	node.level = NumTxAnts - 1;
	Constellation& cstltn = Cstltns[node.level];
	node.psv.push_front(cstltn.ConstlSymbs[0]);
	node.symIndices.push_front(0);
	node.bitLabels.insert(node.bitLabels.end(), cstltn.BitLabels[0].begin(), cstltn.BitLabels[0].end());
}

Vecd SphereDecoder::operator()(ComplexMatrix2D& H, ComplexVec& rxSymbs)
{
	//-----------------------------------------------------------------------------------------------
	// softbits = nrSphereDecoder(H, rxSymbs) uses sphere decoding algorithms (single tree seach)
	// for seeking the maximum - likelihood solution for a set of symbols transmitted over the MIMO channel.
	// ***Input***
	//	* H - complex channel matrix Nr-by-Nt, Nr is the number of rx antennas, Nt is the number of Tx antennas
	//	* rxSymbs - complex received symbol vector  Nr-by-1
	// ***Output*** softbits(llr) - soft bits (not sacled by 1/N0)
	//--------------------------------------------------------------------------------------------------

	// hypothesis maximum likelihood (ML) bits
	Veci bitLabelsML; bitLabelsML.reserve(NumTotBits);
	for (auto i = 0; i < NumTxAnts; i++)
		bitLabelsML.insert(bitLabelsML.end(), Cstltns[i].BitLabels[0].begin(), Cstltns[i].BitLabels[0].end());

	// QR Decomposition
	ComplexMatrix2D Qh, R; qrDecomposeEigen(H, Qh, R);

	// perform Q'*y, y is rx symbol vector
	ComplexVec Qhy; Qhy.reserve(NumTxAnts);
	for (auto i = 0; i < NumTxAnts; i++) {
		Qhy.push_back(inner_product(Qh[i].begin(), Qh[i].end(), rxSymbs.begin(), complex<double>(0.0)));
	}

	// initalize node to root
	TreeNode node; nodeInit(node);

	double di; //distance increment
	Vecd ped(NumTxAnts, 0.0);// partial euclidean distance
	//hypothesis minimum distances for ML bitsand its counterpart
	double lambdaML = std::numeric_limits<double>::infinity(); // minimum distance
	Vecd lambdaMLbar(NumTotBits, std::numeric_limits<double>::infinity());

	// varibles related to tree traveral
	bool isTravelDone = false;
	NumNodesVisted = 0; // only count leaf nodes

	while (!isTravelDone)
	{
		//----------------- disply node ----------------------------------------------------
		//cout << "node ["; for (auto e : node.symIndices) std::cout << e << " "; cout << "\b]" << endl;
		//
		// compute distance increment of current node
		di = norm(Qhy[node.level] -
			inner_product(node.psv.begin(), node.psv.end(), R[node.level].begin() + node.level, complex<double>(0.0)));

		// compute partial Euclidean distance, none root node distance increamented di from its parent node.
		ped[node.level] = (node.level == NumTxAnts - 1) ? di : ped[node.level + 1] + di;

		if (node.level == 0) { // leaf node
			NumNodesVisted += 1;
			if (ped.front() < lambdaML) { // smaller euclidean distance found
				//update the counter hypotheses
				for (auto i = 0; i < NumTotBits; i++) {
					if (node.bitLabels[i] != bitLabelsML[i]) lambdaMLbar[i] = lambdaML;
				}

				//update the hypotheses
				lambdaML = ped.front(); bitLabelsML = Veci(node.bitLabels.begin(), node.bitLabels.end());

				//LLR clipping for counter hypothese
				for (auto& e : lambdaMLbar) e = min(e, SearchRadius + lambdaML);
			}
			else { // No smaller euclidean distance found in the leaf node
				// update the counter hypotheses only
				for (auto i = 0; i < NumTotBits; i++) {
					if ((node.bitLabels[i] != bitLabelsML[i]) & (lambdaMLbar[i] > ped.front())) lambdaMLbar[i] = ped.front();
				}
			}
			isTravelDone = moveRight(node);
		}
		else { // none leaf node
			if (ped[node.level] < lambdaML)
				moveDown(node);
			else { // for hypothesis there is no need to go down
				//Check if there is a smaller Euclidean distance for the couther hypothesis in the sub - tree

				//  find the maximum lambadaMLBar for those bits not yet traversed
				int nBitsNotVisted = accumulate(Ks.begin(), Ks.begin() + node.level, 0);
				double lambaMax = *max_element(lambdaMLbar.begin(), lambdaMLbar.begin() + nBitsNotVisted);

				// find the maximum lambadaMLBar for those bits already traversed
				for (auto i = 0; i < NumTotBits - nBitsNotVisted; i++) {
					if (node.bitLabels[i] != bitLabelsML[nBitsNotVisted + i])
						if (lambdaMLbar[nBitsNotVisted + i] > lambaMax)
							lambaMax = lambdaMLbar[nBitsNotVisted + i];
				}

				if (ped[node.level] < lambaMax) // not punning as counter hypothesis  has a chance to update
					moveDown(node);
				else
					isTravelDone = moveRight(node); // prunning
			}
		}
	}
	// comput soft bits
	Vecd llr; llr.reserve(NumTotBits);
	for (auto i = 0; i < NumTotBits; i++) {
		llr.push_back((bitLabelsML[i] == 1) ? lambdaML - lambdaMLbar[i] : lambdaMLbar[i] - lambdaML);
	}
	return llr;
}

Veci SphereDecoder::hardSphereDecode(ComplexMatrix2D& H, ComplexVec& rxSymbs)
{
	//-----------------------------------------------------------------------------------------------
	// ***Input***
	//	* H - complex channel matrix Nr-by-Nt, Nr is the number of rx antennas, Nt is the number of Tx antennas
	//	* rxSymbs - complex received symbol vector  Nr-by-1
	// ***Output***
	//  * hard detected bits
	//--------------------------------------------------------------------------------------------------

	// hypothesis maximum likelihood (ML) bits
	Veci bitLabelsML; bitLabelsML.reserve(NumTotBits);
	for (auto i = 0; i < NumTxAnts; i++)
		bitLabelsML.insert(bitLabelsML.end(), Cstltns[i].BitLabels[0].begin(), Cstltns[i].BitLabels[0].end());

	// QR Decomposition
	ComplexMatrix2D Qh, R; qrDecomposeEigen(H, Qh, R);

	// perform Q'*y, y is rx symbol vector
	ComplexVec Qhy; Qhy.reserve(NumTxAnts);
	for (auto i = 0; i < NumTxAnts; i++) {
		Qhy.push_back(inner_product(Qh[i].begin(), Qh[i].end(), rxSymbs.begin(), complex<double>(0.0)));
	}

	// initalize node to root
	TreeNode node;
	// add the first symbol of the costellation for the last antenna
	node.level = NumTxAnts - 1;
	Constellation& cstltn = Cstltns[node.level];
	node.psv.push_front(cstltn.ConstlSymbs[0]);
	node.symIndices.push_front(0);
	node.bitLabels.insert(node.bitLabels.end(), cstltn.BitLabels[0].begin(), cstltn.BitLabels[0].end());

	Vecd ped(NumTxAnts, 0.0); // partial euclidean distance
	double di; //distance increment
	double lambdaML = std::numeric_limits<double>::infinity();

	// varibles related to tree traveral
	bool isTravelDone = false;
	NumNodesVisted = 0; // only count leaf nodes

	while (!isTravelDone)
	{
		//----------------- disply node ----------------------------------------------------
		//cout << "node ["; for (auto e : node.symIndices) std::cout << e << " "; cout << "\b]" << endl;

		// compute distance increment
		di = norm(Qhy[node.level] -
			inner_product(node.psv.begin(), node.psv.end(), R[node.level].begin() + node.level, complex<double>(0.0)));

		// compute partial Euclidean distance, none root node distance increamented di from its parent node.
		ped[node.level] = (node.level == NumTxAnts - 1) ? di : ped[node.level + 1] + di;

		if (node.level == 0) { // leaf node
			NumNodesVisted++;
			if (ped.front() < lambdaML) { // smaller euclidean distance found
				//update the hypotheses
				lambdaML = ped.front(); bitLabelsML = Veci(node.bitLabels.begin(), node.bitLabels.end());
			}

			isTravelDone = moveRight(node);
		}
		else { // none leaf node
			if (ped[node.level] < lambdaML)
				moveDown(node);
			else
				isTravelDone = moveRight(node); // trunc
		}
	}

	return bitLabelsML;
}

void SphereDecoder::qrDecomposeEigen(const ComplexMatrix2D& H, ComplexMatrix2D& Qh, ComplexMatrix2D& R)
{
	//produces an economy-size decomposition, computes only the first N columns of Q and the first N rows of R.
	// Qh is hermitan tranpose of Q sized N-by-M;

	using namespace Eigen;

	Qh = ComplexMatrix2D(NumTxAnts); //hermitan tranpose of Q sized N-by-M
	R = ComplexMatrix2D(NumTxAnts);

	// Eigen domain perform QR decompostion
	typedef Matrix<complex<double>, Dynamic, Dynamic> cmat;
	cmat EigenH(NumRxAnts, NumTxAnts);
	for (auto i = 0; i < NumRxAnts; i++)
		for (auto j = 0; j < NumTxAnts; j++)
			EigenH(i, j) = H[i][j];

	HouseholderQR<cmat> qr(EigenH.rows(), EigenH.cols());
	qr.compute(EigenH);
	cmat EigenQ = qr.householderQ() * cmat::Identity(EigenH.rows(), EigenH.cols());
	cmat temp = qr.matrixQR().triangularView<Upper>();
	cmat EigenR = temp.topRows(EigenH.cols());

	// hermitan transpose
	cmat EigenQh = EigenQ.adjoint();

	//back to std namespace
	for (auto i = 0; i < NumTxAnts; i++) {
		Qh[i].reserve(NumRxAnts);
		for (auto j = 0; j < NumRxAnts; j++) {
			Qh[i].push_back(EigenQh(i, j));
		}
	}

	for (auto i = 0; i < NumTxAnts; i++) {
		R[i].reserve(NumTxAnts);
		for (auto j = 0; j < NumTxAnts; j++) {
			R[i].push_back(EigenR(i, j));
		}
	}
}

void SphereDecoder::moveDown(TreeNode& node)
{
	// node move done a level, always add the first symbol to the psv in the constellation of that level
	node.level--;
	Constellation& cstl = Cstltns[node.level];
	node.symIndices.push_front(0);
	node.psv.push_front(cstl.ConstlSymbs[0]);
	node.bitLabels.insert(node.bitLabels.begin(), cstl.BitLabels[0].begin(), cstl.BitLabels[0].end());
};

bool SphereDecoder::moveRight(TreeNode& node)
{
	// try to move right, if already the rightest node of a subtree, then move up a level and go right
	while (true) {
		Constellation& cstl = Cstltns[node.level];
		if (node.symIndices[0] < cstl.M - 1) { // not yet to the far right,
			//level no change, the fist symbol of the psv will be changed to the next symbol in the constellation
			int symbIdx = node.symIndices[0] + 1;
			node.symIndices[0] = symbIdx;
			node.psv[0] = cstl.ConstlSymbs[symbIdx];
			for (auto i = 0; i < cstl.K; i++) node.bitLabels[i] = cstl.BitLabels[symbIdx][i];
			return false;
		}
		else { // reach the rightest child of the parent node, UP
			if (node.level == NumTxAnts - 1) // root node
				return true; // traveral completed
			else { // go up a level
				node.psv.pop_front();
				node.symIndices.pop_front();
				for (auto i = 0; i < cstl.K; i++) node.bitLabels.pop_front();
				node.level++;
			}
		}
	}
};

Veci de2bi(int decNum, int NOut, string msb)
{
	// This function Convert decimal numbers to binary numbers.
	int tmp = decNum;
	Veci out(NOut, 0);
	if (msb == "right-msb") {
		int i = 0;
		while ((tmp > 0) & (i < NOut)) {
			out[i++] = tmp % 2;
			tmp /= 2;
		}
	}
	else {
		int i = NOut - 1;
		while ((tmp > 0) & (i >= 0)) {
			out[i--] = tmp % 2;
			tmp /= 2;
		}
	}

	return out;
}