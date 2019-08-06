/* $Id: TBlock.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/TBlock.h>
using namespace crs;

#include <iostream>
#include <sstream>
#include <string>

const int TBlock::fgNSubBlocks = 21;
const int TBlock::fgNParticles = 39;
const int TBlock::fgNLongitudinals = 26;
const int TBlock::fgNEntriesNotThinned = 7;
const int TBlock::fgNEntriesThinned = 8;


TBlock::TBlock(const CREAL *data, bool thinned) 
{    
  fThinned = thinned;

  // FORTRAN:  data (MAXBUF,NSUBBL)
  // -> data [nsubbl] [maxbuf]

  const int EntriesPerParticle 
    = (fThinned ? fgNEntriesThinned : fgNEntriesNotThinned);
    
  for (int iSubBlock=0;
       iSubBlock < fgNSubBlocks;
       ++iSubBlock) {

    int index = (iSubBlock * 
		 fgNParticles * 
		 EntriesPerParticle);


    TSubBlock testSubBlock (data+index, thinned); 
	
    fSubBlocks.push_back (testSubBlock);
	
  }
    
}




