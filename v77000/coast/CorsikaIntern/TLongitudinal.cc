/* $Id: TLongitudinal.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/TLongitudinal.h>
using namespace crs;

#include <iostream>


TLongitudinal::TLongitudinal (const float *data) {
    
  fData = data;
}


void TLongitudinal::Dump () {

  std::cout << "long:> " 
	    << " X/(gcm^{-2}): " << GetSlantDepth ()
	    << " N_g " << GetNGammas ()
	    << " N_e- " << GetNElectrons ()
	    << " N_e+ " << GetNPositrons ()
	    << " N_mu- " << GetNMuons ()
	    << " N_mu+ " << GetNAntiMuons ()
	    << " N_h " << GetNHadrons ()
	    << " N_ch " << GetNCharged ()
	    << " N_nuc " << GetNNuclei ()
	    << " N_ckov " << GetNCherenkov ()
	    << std::endl;
}
