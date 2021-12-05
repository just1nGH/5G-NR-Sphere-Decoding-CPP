#include "nrModulation.h"
using namespace std;
//===================================== NR Modulation Mapper ===================================
ComplexVec nrModuMapper(Veci& bitsIn, string moduType)
{
	//----------------------------------------------------------------------------------------
	// This function maps the bit sequence into 
	// complex modulation symbols using modulation scheme specified in TS 38.211 Section 5.1
	// The modulation scheme, 'moduType' must be one of  'BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'.
	// Author Dr. J Mao 2021 Nov
	//----------------------------------------------------------------------------------------------
	Veci& b = bitsIn;
	std::transform(moduType.begin(), moduType.end(), moduType.begin(),
		[](unsigned char c) { return std::tolower(c); });

	if (moduType == "bpsk") {
		ComplexVec symbOut(b.size());
		for (auto i = 0; i < b.size(); i++) {
			symbOut[i] = 1.0 / sqrt(2) * complex<double>(1.0 - 2.0 * b[i], 1.0 - 2.0 * b[i]);
		}
		return symbOut;
	}
	else if (moduType == "qpsk") {
		assert(b.size() % 2 == 0);
		int nOutSymbs = b.size() / 2;
		ComplexVec symbOut(nOutSymbs);
		for (auto i = 0; i < nOutSymbs; i++) {
			symbOut[i] = 1.0 / sqrt(2) * complex<double>(1.0 - 2.0 * b[2 * i], 1.0 - 2.0 * b[2 * i + 1]);
		}
		return symbOut;

	}
	else if (moduType == "16qam") {
		assert(b.size() % 4 == 0);
		int nOutSymbs = b.size() / 4;
		ComplexVec symbOut(nOutSymbs);
		for (auto i = 0; i < nOutSymbs; i++) {
			symbOut[i] = 1.0 / sqrt(10) * complex<double>((1.0 - 2.0 * b[4 * i]) * (2.0 - (1 - 2.0 * b[4 * i + 2])),
				(1.0 - 2.0 * b[4 * i + 1]) * (2.0 - (1 - 2.0 * b[4 * i + 3])));
		}
		return symbOut;
	}

	else if (moduType == "64qam") {
		assert(b.size() % 6 == 0);
		int nOutSymbs = b.size() / 6;
		ComplexVec symbOut(nOutSymbs);
		for (auto i = 0; i < nOutSymbs; i++) {
			symbOut[i] = 1.0 / sqrt(42) * complex<double>((1.0 - 2.0 * b[6 * i]) * (4.0 - (1.0 - 2.0 * b[6 * i + 2]) * (2.0 - (1.0 - 2.0 * b[6 * i + 4]))),
				(1.0 - 2.0 * b[6 * i + 1]) * (4.0 - (1.0 - 2.0 * b[6 * i + 3]) * (2.0 - (1.0 - 2.0 * b[6 * i + 5]))));
		}
		return symbOut;
	}
	else if (moduType == "256qam") {
		assert(b.size() % 8 == 0);
		int nOutSymbs = b.size() / 8;
		ComplexVec symbOut(nOutSymbs);
		for (auto i = 0; i < nOutSymbs; i++) {
			symbOut[i] = 1.0 / sqrt(170) * complex<double>((1.0 - 2.0 * b[8 * i]) * (8.0 - (1.0 - 2.0 * b[8 * i + 2]) * (4.0 - (1.0 - 2.0 * b[8 * i + 4]) * (2.0 - (1.0 - 2.0 * b[8 * i + 6])))), //real part
				(1.0 - 2.0 * b[8 * i + 1]) * (8.0 - (1.0 - 2.0 * b[8 * i + 3]) * (4.0 - (1.0 - 2.0 * b[8 * i + 5]) * (2.0 - (1.0 - 2.0 * b[8 * i + 7]))))); // imag part

		}
		return symbOut;
	}
	else {
		return ComplexVec();
	}

}


//===================================== NR soft Modulation DeMapper ===================================
Vecd nrSoftModuDemapper(ComplexVec symbsIn, string moduType, double N0, string method)
{
	//----------------------------------------------------------------------------------------------------------
	//  nrSoftModuDemapper demodulates symbols to soft bits
	//   softBits = nrSoftModuDemapper(symbsIn,moduType,N0,method) demodulates the complex symbols
	//   symbsIn using soft decision. The modulation scheme, moduType must be one
	//   of  'BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'. symbsIn must be
	//   a column vector.
	//   'mothod'   -   Specified as 'max-log-map' or 'approx'.
	//
	//   Demodulation is performed according to the constellations given in
	//   TS 38.211 section 5.1 including the power normalization factors
	//   specified.
	//
	//   Source Paper:   Mao, Juquan, et al. "A low complexity 256QAM soft demapper for 5G mobile system." 
	//	 Author Dr. J Mao NOv 2021
	//---------------------------------------------------------------------------------------------------------------
	std::transform(moduType.begin(), moduType.end(), moduType.begin(),
		[](unsigned char c) { return std::tolower(c); });

	// avoid unbounded softbits 
	N0 = (N0 < 0.001) ? 0.001 : N0;

	if (moduType == "bpsk") {
		return bpskSoftDemodulation(symbsIn, N0);
	}
	else if (moduType == "qpsk") {
		return qpskSoftDemodulation(symbsIn, N0);
	}
	else if (moduType == "16qam") {
		return qam16SoftDemodulation(symbsIn, N0,  method);
	}

	else if (moduType == "64qam") {
		return qam64SoftDemodulation(symbsIn, N0, method);

	}
	else if (moduType == "256qam") {
		return qam256SoftDemodulation(symbsIn, N0, method);

	}
	else {
		return Vecd();
	}
}

//------------------------------------------------- bpsk----------------------------------------------
Vecd bpskSoftDemodulation(ComplexVec symbsIn, double N0)
{
	auto nSymbs = symbsIn.size();
	Vecd softBits(nSymbs);
	double A = 4.0 / N0 * 1 / sqrt(2);
	for (auto i = 0; i < nSymbs; i++) {
		softBits[i] = A * (complex<double>(1, -1) * symbsIn[i]).real();
	}
	return softBits;
}

//-------------------------------------------------qpsk----------------------------------------------
Vecd qpskSoftDemodulation(ComplexVec symbsIn, double N0)
{
	auto nSymbs = symbsIn.size();
	Vecd softBits(nSymbs * 2);
	double A = 4.0 / N0 * 1 / sqrt(2);
	for (auto i = 0; i < nSymbs; i++) {
		softBits[2 * i] = A * symbsIn[i].real();
		softBits[2 * i + 1] = A * symbsIn[i].imag();
	}
	return softBits;
}

//------------------------------------------------- 16QAM----------------------------------------------
Vecd qam16SoftDemodulation(ComplexVec symbsIn, double N0, string method)
{

	auto nSymbs = symbsIn.size();

	Vecd softBits(nSymbs * 4);
	double A =  1 / sqrt(10);

	double reX,imX;

	std::transform(method.begin(), method.end(), method.begin(),
		[](unsigned char c) { return std::tolower(c); });

	if (method == "approx") {
		for (auto i = 0; i < nSymbs; i++) {
			reX = symbsIn[i].real();
			imX = symbsIn[i].imag();
			softBits[4 * i    ] = 4 * A / N0 * reX;
			softBits[4 * i + 1] = 4 * A / N0 * imX;
			softBits[4 * i + 2] = 4 * A / N0 * (-abs(reX) + 2.0 * A);
			softBits[4 * i + 3] = 4 * A / N0 * (-abs(imX) + 2.0 * A);
		}
		return softBits;
	}
	else {

		for (auto i = 0; i < nSymbs; i++) {
			reX = symbsIn[i].real();
			imX = symbsIn[i].imag();
			softBits[4 * i] = maxLogMapQam16SoftDemodulation_0_1(reX, N0);
			softBits[4 * i + 1] = maxLogMapQam16SoftDemodulation_0_1(imX, N0);
			softBits[4 * i + 2] = 4 * A / N0 * (-abs(reX) + 2.0 * A);
			softBits[4 * i + 3] = 4 * A / N0 * (-abs(imX) + 2.0 * A);
		}
		return softBits;

	}

}
double maxLogMapQam16SoftDemodulation_0_1(double realSymbIn, double N0)
{
	double A = 1 / sqrt(10);
	double X = realSymbIn;

	if (X < -2*A)
		return 8 * A / N0 * (X + A);
	else if (X < 2 * A)
		return 4 * A / N0 * X;
	else
		return 8 * A / N0 * (X - A);
}

//------------------------------------------------- 64QAM----------------------------------------------
Vecd qam64SoftDemodulation(ComplexVec symbsIn, double N0, string method)
{

	auto nSymbs = symbsIn.size();

	Vecd softBits(nSymbs * 6);
	double A = 1 / sqrt(42);

	double reX, imX;

	std::transform(method.begin(), method.end(), method.begin(),
		[](unsigned char c) { return std::tolower(c); });

	if (method == "approx") {
		for (auto i = 0; i < nSymbs; i++) {
			reX = symbsIn[i].real();
			imX = symbsIn[i].imag();
			softBits[6 * i] = 4 * A / N0 * reX;
			softBits[6 * i + 1] = 4 * A / N0 * imX;
			softBits[6 * i + 2] = 4 * A / N0 * (-abs(reX) + 4.0 * A);
			softBits[6 * i + 3] = 4 * A / N0 * (-abs(imX) + 4.0 * A);
			softBits[6 * i + 4] = 4 * A / N0 * (-abs(-abs(reX) + 4 * A) + 2 * A);
			softBits[6 * i + 5] = 4 * A / N0 * (-abs(-abs(imX) + 4 * A) + 2 * A);
		}
		return softBits;
	}
	else {

		for (auto i = 0; i < nSymbs; i++) {
			reX = symbsIn[i].real();
			imX = symbsIn[i].imag();
			softBits[6 * i] = maxLogMapQam64SoftDemodulation_0_1(reX, N0);
			softBits[6 * i + 1] = maxLogMapQam64SoftDemodulation_0_1(imX, N0);
			softBits[6 * i + 2] = maxLogMapQam64SoftDemodulation_2_3(reX, N0);
			softBits[6 * i + 3] = maxLogMapQam64SoftDemodulation_2_3(imX, N0);
			softBits[6 * i + 4] = 4 * A / N0 * (-abs(-abs(reX) + 4 * A) + 2 * A);
			softBits[6 * i + 5] = 4 * A / N0 * (-abs(-abs(imX) + 4 * A) + 2 * A);
		}
		return softBits;

	}

}
double maxLogMapQam64SoftDemodulation_0_1(double realSymbIn, double N0)
{
	double A = 1 / sqrt(42);
	double X = realSymbIn;

	if (X < -6 * A)
		return 16 * A / N0 * (X + 3 * A);
	else if (X < -4 * A)
		return 12 * A / N0 * (X + 2 * A);
	else if (X < -2 * A)
		return 8 * A / N0 * (X + A);
	else if (X < 2 * A)
		return 4 * A / N0 * X;
	else if (X < 4 * A)
		return 8 * A / N0 * (X - A);
	else if (X < 6 * A)
		return 12 * A / N0 * (X - 2 * A);
	else
		return 16 * A / N0 * (X - 3 * A);
}
double maxLogMapQam64SoftDemodulation_2_3(double realSymbIn, double N0)
{
	double A = 1 / sqrt(42);
	double X = realSymbIn;

	if (X < -6 * A)
		return  8 * A / N0 * (X + 5 * A);
	else if (X < -2 * A)
		return 4 * A / N0 * (X + 4 * A);
	else if (X < 0)
		return 8 * A / N0 * (X + 3 * A);
	else if (X < 2 * A)
		return 8 * A / N0 * (-X + 3 * A);
	else if (X < 6 * A)
		return 4 * A / N0 * (-X + 4 * A);
	else
		return 8 * A / N0 * (-X + 5 * A);

}


//------------------------------------------------- 256QAM----------------------------------------------
Vecd qam256SoftDemodulation(ComplexVec symbsIn, double N0, string method)
{

	auto nSymbs = symbsIn.size();

	Vecd softBits(nSymbs * 8);
	double A = 1 / sqrt(170);

	double reX, imX;

	std::transform(method.begin(), method.end(), method.begin(),
		[](unsigned char c) { return std::tolower(c); });

	if (method == "approx") {
		for (auto i = 0; i < nSymbs; i++) {
			reX = symbsIn[i].real();
			imX = symbsIn[i].imag();
			softBits[8 * i] = 4 * A / N0 * reX;
			softBits[8 * i + 1] = 4 * A / N0 * imX;
			softBits[8 * i + 2] = 4 * A / N0 * (-abs(reX) + 8.0 * A);
			softBits[8 * i + 3] = 4 * A / N0 * (-abs(imX) + 8.0 * A);
			softBits[8 * i + 4] = 4 * A / N0 * (-abs(-abs(reX) + 8 * A) + 4 * A);
			softBits[8 * i + 5] = 4 * A / N0 * (-abs(-abs(imX) + 8 * A) + 4 * A);
			softBits[8 * i + 6] = 4 * A / N0 * (-abs(-abs(-abs(reX) + 8 * A) + 4 * A) + 2 * A);
			softBits[8 * i + 7] = 4 * A / N0 * (-abs(-abs(-abs(imX) + 8 * A) + 4 * A) + 2 * A);
		}
		return softBits;
	}
	else {

		for (auto i = 0; i < nSymbs; i++) {
			reX = symbsIn[i].real();
			imX = symbsIn[i].imag();
			softBits[8 * i]		= maxLogMapQam256SoftDemodulation_0_1(reX, N0);
			softBits[8 * i + 1] = maxLogMapQam256SoftDemodulation_0_1(imX, N0);
			softBits[8 * i + 2] = maxLogMapQam256SoftDemodulation_2_3(reX, N0);
			softBits[8 * i + 3] = maxLogMapQam256SoftDemodulation_2_3(imX, N0);
			softBits[8 * i + 4] = maxLogMapQam256SoftDemodulation_4_5(reX, N0);
			softBits[8 * i + 5] = maxLogMapQam256SoftDemodulation_4_5(imX, N0);
			softBits[8 * i + 6] = 4 * A / N0 * (-abs(-abs(-abs(reX) + 8 * A) + 4 * A) + 2 * A);
			softBits[8 * i + 7] = 4 * A / N0 * (-abs(-abs(-abs(imX) + 8 * A) + 4 * A) + 2 * A);
		}
		return softBits;

	}

}
double maxLogMapQam256SoftDemodulation_0_1(double realSymbIn, double N0)
{
	double A = 1 / sqrt(170);
	double X = realSymbIn;

	if (X < -14 * A)
		return 32 * A / N0 * (X + 7 * A);
	else if(X < -12 * A)
		return 28 * A / N0 * (X + 6 * A);
	else if(X < -10 * A)
		return 24 * A / N0 * (X + 5 * A);
	else if(X < -8 * A)
		return 20 * A / N0 * (X + 4 * A);
	else if(X < -6 * A)
		return 16 * A / N0 * (X + 3 * A);
	else if(X < -4 * A)
		return 12 * A / N0 * (X + 2 * A);
	else if(X < -2 * A)
		return 8 * A / N0 * (X + A);
	else if(X < 2 * A)
		return 4 * A / N0 * X;
	else if(X < 4 * A)
		return 8 * A / N0 * (X - A);
	else if(X < 6 * A)
		return 12 * A / N0 * (X - 2 * A);
	else if(X < 8 * A)
		return 16 * A / N0 * (X - 3 * A);
	else if(X < 10 * A)
		return 20 * A / N0 * (X - 4 * A);
	else if(X < 12 * A)
		return 24 * A / N0 * (X - 5 * A);
	else if(X < 14 * A)
		return 28 * A / N0 * (X - 6 * A);
	else
		return 32 * A / N0 * (X - 7 * A);

}
double maxLogMapQam256SoftDemodulation_2_3(double realSymbIn, double N0)
{
	double A = 1 / sqrt(170);
	double X = realSymbIn;

	if (X < -14 * A)
		return 16 * A / N0 * (X + 11 * A);
	else if (X < -12 * A)
		return 12 * A / N0 * (X + 10 * A);
	else if (X < -10 * A)
		return 8 * A / N0 * (X + 9 * A);
	else if (X < -6 * A)
		return 4 * A / N0 * (X + 8 * A);
	else if (X < -4 * A)
		return 8 * A / N0 * (X + 7 * A);
	else if (X < -2 * A)
		return 12 * A / N0 * (X + 6 * A);
	else if (X < 0)
		return 16 * A / N0 * (X + 5 * A);
	else if (X < 2 * A)
		return 16 * A / N0 * (-X + 5 * A);
	else if (X < 4 * A)
		return 12 * A / N0 * (-X + 6 * A);
	else if (X < 6 * A)
		return 8 * A / N0 * (-X + 7 * A);
	else if (X < 10 * A)
		return 4 * A / N0 * (-X + 8 * A);
	else if (X < 12 * A)
		return 8 * A / N0 * (-X + 9 * A);
	else if (X < 14 * A)
		return 12 * A / N0 * (-X + 10 * A);
	else
		return 16 * A / N0 * (-X + 11 * A);

}
double maxLogMapQam256SoftDemodulation_4_5(double realSymbIn, double N0)
{
	double A = 1 / sqrt(170);
	double X = realSymbIn;

	if (X < -14 * A)
		return 8 * A / N0 * (X + 13 * A);
	else if(X < -10 * A)
		return 4 * A / N0 * (X + 12 * A);
	else if(X < -8 * A)
		return 8 * A / N0 * (X + 11 * A);
	else if(X < -6 * A)
		return 8 * A / N0 * (-X - 5 * A);
	else if(X < -2 * A)
		return 4 * A / N0 * (-X - 4 * A);
	else if(X < 0)
		return 8 * A / N0 * (-X - 3 * A);
	else if(X < 2 * A)
		return 8 * A / N0 * (X - 3 * A);
	else if(X < 6 * A)
		return 4 * A / N0 * (X - 4 * A);
	else if(X < 8 * A)
		return 8 * A / N0 * (X - 5 * A);
	else if(X < 10 * A)
		return 8 * A / N0 * (-X + 11 * A);
	else if(X < 14 * A)
		return 4 * A / N0 * (-X + 12 * A);
	else
		return 8 * A / N0 * (-X + 13 * A);

}

