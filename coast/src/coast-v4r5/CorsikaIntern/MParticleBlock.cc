/* $Id: MParticleBlock.cc,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $   */

#include <crs/TBlock.h>
#include <crs/MParticleBlock.h>
#include <crs/TParticleBlockEntry.h>
#include <crs/MCherenkov.h>
#include <crs/MParticle.h>
#include <crs/MMuonProductionInfo.h>
using namespace crs;


MParticleBlock::MParticleBlock (const TSubBlock &right) {
    
  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
    
    
  // find particles  
  int EntriesPerParticle = (fThinned ?
			    TBlock::fgNEntriesThinned :
			    TBlock::fgNEntriesNotThinned);
    

  for (int iParticle=0; 
       iParticle<TBlock::fgNParticles;
       iParticle++) {
	
    int index = EntriesPerParticle * iParticle;
	
    TParticleBlockEntry testParticle (fSubBlockData+index, fThinned);
	

    if (testParticle.IsEmpty ()) {
      break;
    }
	
    if (testParticle.IsCherenkov ()) {
      fParticles.push_back ((MCherenkov)testParticle);
    }

    if (testParticle.IsParticle () || 
	testParticle.IsNucleus ()) {
      fParticles.push_back ((MParticle)testParticle);
    }

    if (testParticle.IsMuonProductionInfo ()) {
      fParticles.push_back ((MMuonProductionInfo)testParticle);
    }

	
  }
    
}
