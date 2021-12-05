# 5G-NR-Sphere-Decoding-CPP

C++ Implenmentation of 5G NR MIMO Sphere Decoder
- Tree traverse stratigies: single tree search 
- Support  different modulation schemes for Tx antennas

The implementation sphere decoding algorithms(single tree seach)for seeking the maximum-likelihood solution 
for a set of symbols transmitted over the MIMO channel. 
```
	//-----------------------------------------------------------------------------------------------
	// softbits = nrSphereDecoder(H, rxSymbs) uses sphere decoding algorithms (single tree seach)
	// for seeking the maximum - likelihood solution for a set of symbols transmitted over the MIMO channel.
	// ***Input***
	//	* H - complex channel matrix Nr-by-Nt, Nr is the number of rx antennas, Nt is the number of Tx antennas
	//	* rxSymbs - complex received symbol vector  Nr-by-1
	// ***Output*** softbits(llr) - soft bits (not sacled by 1/N0)
	//--------------------------------------------------------------------------------------------------
% Author: Dr J Mao
% Email: juquan.justin.mao@gmail.com
% 2021 Nov
%--------------------------------------------------------------------------
```

### Example run simlation results
Start by try `example_run.m`
