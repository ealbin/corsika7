/* $Id: MMuonProductionInfo.cc,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $   */

#include <crs/MMuonProductionInfo.h>
using namespace crs;

#include <iostream>
#include <string>


MMuonProductionInfo::MMuonProductionInfo (const float *data, bool thinned) :
TParticleBlockEntry (data, thinned) {
}


MMuonProductionInfo::MMuonProductionInfo (const TParticleBlockEntry &p) :
TParticleBlockEntry (p) {
}

std::string MMuonProductionInfo::GetParticleName () const {

  return std::string ("muon prod. info");
}


void MMuonProductionInfo::Dump () const {
    
  std::cout << *this 
	    << std::endl;
}



