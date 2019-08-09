/* $Id: MLongitudinalBlock.cc,v 1.1.1.1 2007-07-31 07:00:51 rulrich Exp $   */

#include <crs/TBlock.h>
#include <crs/MLongitudinalBlock.h>
using namespace crs;

const int MLongitudinalBlock::fgEntriesPerLongInfo = 10;
const int MLongitudinalBlock::fgLongInfoOffset = 13;

MLongitudinalBlock::MLongitudinalBlock (const TSubBlock &right) {
    
  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;  

  ScanLongBlock ();
}


void MLongitudinalBlock::ScanLongBlock () {
    
  for (int iLongitudinal = 0; 
       iLongitudinal < TBlock::fgNLongitudinals;
       ++iLongitudinal) {
	
    int CurrentStep = 
      (GetBlockNumber () - 1) * TBlock::fgNLongitudinals + iLongitudinal;

    if (CurrentStep>=GetNSteps ()) {
      break;
    }

    int index = fgLongInfoOffset + fgEntriesPerLongInfo * iLongitudinal;
      
    TLongitudinal longData (fSubBlockData+index);
    fLongitudinals.push_back (longData);
	
  }

}
