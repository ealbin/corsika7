/* $Id: MParticleBlock.cc 5924 2017-03-23 21:47:23Z rulrich $   */

#include <crs/TBlock.h>
#include <crs/MParticleBlock.h>
#include <crs/TParticleBlockEntry.h>
#include <crs/MCherenkov.h>
#include <crs/MParticle.h>
#include <crs/MMuonProductionInfo.h>

#include <iostream>
using namespace std;
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
  
  
  int lastEmpty = TBlock::fgNParticles - 1;
  for (int iParticle=0;
       iParticle<TBlock::fgNParticles;
       iParticle++) {
    
    int index = EntriesPerParticle * iParticle;
    
    TParticleBlockEntry testParticle (fSubBlockData+index, fThinned);
    
    if (testParticle.IsEmpty ()) {
      lastEmpty = iParticle;
      continue;
    }
    
    else if (testParticle.IsCherenkov ()) {
      fParticles.push_back ((MCherenkov)testParticle);
    }
    
    else if (testParticle.IsParticle () ||
	testParticle.IsNucleus ()) {
      fParticles.push_back ((MParticle)testParticle);
    }
    
    else if (testParticle.IsMuonProductionInfo ()) {
      fParticles.push_back ((MMuonProductionInfo)testParticle);
    }
  }

  /*
  if (lastEmpty != TBlock::fgNParticles - 1)
    cerr << " MParticleBlock::MParticleBlock() - WARNING:\n"
         << "\t ==> empty particle block within TSubBlock! \n"
         << "\t ==> (last empty = " << lastEmpty << ", particles in block: "
         << TBlock::fgNParticles << ")" << endl;
  */
}
