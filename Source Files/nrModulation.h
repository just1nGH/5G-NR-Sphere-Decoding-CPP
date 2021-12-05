#ifndef NR_MODULATION_H
#define NR_MODULATION_H

#include <vector>
#include <string>
#include <complex>
#include <assert.h>
#include <algorithm>

#include <vector>
#include <complex>

typedef std::vector<std::complex<double>>  ComplexVec;
typedef std::vector<ComplexVec>  ComplexMatrix2D;
typedef std::vector<ComplexMatrix2D>  ComplexMatrix3D;
typedef std::vector<ComplexMatrix3D>  ComplexMatrix4D;

typedef std::vector<int> Veci;
typedef std::vector<Veci>  Matrix2Di;
typedef std::vector<Matrix2Di>  Matrix3Di;
typedef std::vector<Matrix3Di>  Matrix4Di;

typedef std::vector<double>  Vecd;
typedef std::vector<Vecd>  Matrix2Dd;
typedef std::vector<Matrix2Dd>  Matrix3Dd;
typedef std::vector<Matrix3Dd>  Matrix4Dd;

//===================================== Modulation Information ===================================
struct ModuInfo {
	std::string modulation;
	int M;  // number of points in the constellation
	int Qm; // number of bits per complex symbol
	ModuInfo(std::string moduType) {
		std::transform(moduType.begin(), moduType.end(), moduType.begin(),
			[](unsigned char c) { return std::tolower(c); });
		modulation = moduType;
		if (moduType == "bpsk") {
			M = 2;
			Qm = 1;
		}
		else if (moduType == "qpsk") {
			M = 4;
			Qm = 2;
		}
		else if (moduType == "16qam") {
			M = 16;
			Qm = 4;
		}
		else if (moduType == "64qam") {
			M = 64;
			Qm = 6;
		}
		else if (moduType == "256qam") {
			M = 256;
			Qm = 8;
		}
		else {
		}
	}
};

//===================================== NR Modulation Mapper ===================================
// ComplexMatrix2D nrModuMapper(Veci& bitsIn, std::string moduType)
// This function maps the bit sequence into
// complex modulation symbols using modulation scheme specified in TS 38.211 Section 5.1
// The modulation scheme, 'moduType' must be one of  'BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'.
// Author Dr. J Mao 2021 Nov
//----------------------------------------------------------------------------------------------
ComplexVec nrModuMapper(Veci& bitsIn, std::string moduType);

//===================================== NR soft Modulation DeMapper ===================================
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
Vecd nrSoftModuDemapper(ComplexVec symbsIn, std::string moduType, double N0, std::string method = "max-log-map");
Vecd bpskSoftDemodulation(ComplexVec symbsIn, double N0);
Vecd qpskSoftDemodulation(ComplexVec symbsIn, double N0);
Vecd qam16SoftDemodulation(ComplexVec symbsIn, double N0, std::string method = "max-log-map");
double maxLogMapQam16SoftDemodulation_0_1(double realSymbIn, double N0);
Vecd qam64SoftDemodulation(ComplexVec symbsIn, double N0, std::string method = "max-log-map");
double maxLogMapQam64SoftDemodulation_0_1(double realSymbIn, double N0);
double maxLogMapQam64SoftDemodulation_2_3(double realSymbIn, double N0);
Vecd qam256SoftDemodulation(ComplexVec symbsIn, double N0, std::string method = "max-log-map");
double maxLogMapQam256SoftDemodulation_0_1(double realSymbIn, double N0);
double maxLogMapQam256SoftDemodulation_2_3(double realSymbIn, double N0);
double maxLogMapQam256SoftDemodulation_4_5(double realSymbIn, double N0);

#endif
