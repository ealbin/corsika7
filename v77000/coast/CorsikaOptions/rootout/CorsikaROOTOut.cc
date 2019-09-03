/* $Id: CorsikaROOTOut.cc 6252 2017-09-01 07:39:34Z winchen $   */

#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>

#include <crs2r/TC2R.h>

#include <sstream>
#include <iostream>



///  Global object to handle all ROOT I/O operations and data conversions
crs2r::TC2R gC2R;



/*
  Data is one CORSIKA data-block constining of 21 SubBlocks.
  A SubBlock can be:
  - thinned mode:     39 (Particles) * 8 (ENTRIES) * 4 (BYTES) 
  - not-thinned mode: 39 (Particles) * 7 (ENTRIES) * 4 (BYTES) 
*/
void
wrida_ (const CREAL *Data) 
{    
  crs::TSubBlock subBlock (Data, gC2R.IsThinned ());    
  gC2R.Write (subBlock);
}



void
inida_ (const char* filename,
        const int& thinning,
        const int& /*curved*/,
        const int& /*slant*/,
        const int& /*stackinput*/,
        const int& /*preshower*/,
        int str_length) 
{  
  std::ostringstream FileName;
  for (int i=0; i<str_length; i++) {

    if (filename [i]==' ')
      break;

    FileName << filename [i];
  }

  FileName << ".root";

  std::cout << " ######################################################### "
	    << std::endl
	    << " open file: " << FileName.str ()
	    << " thinned: " << thinning
	    << std::endl;

  gC2R.Open (FileName.str (), thinning);
}

void
cloda_ () 
{
  std::cout << " ######################################################### "
	    << std::endl
	    << " closing "
	    << std::endl;

  gC2R.Close ();
}


void 
interaction_ (const crs::CInteraction& interaction)
{
}


void
track_ (const crs::CParticle& pre, const crs::CParticle& post)
{
} 


void tabularizedatmosphere_(const int &nPoints, const double* height, const double* refractiveIndex)
{}
