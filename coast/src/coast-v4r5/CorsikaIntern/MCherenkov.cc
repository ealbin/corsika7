/* $Id: MCherenkov.cc,v 1.1.1.1 2007-07-31 07:00:52 rulrich Exp $   */

#include <crs/MCherenkov.h>
#include <crs/TParticleBlockEntry.h>
using namespace crs;

#include <iostream>
#include <string>



MCherenkov::MCherenkov (const float *data, bool thinned) :
TParticleBlockEntry (data,  thinned ) {
}

MCherenkov::MCherenkov (const TParticleBlockEntry &p) :
TParticleBlockEntry (p) {
}

void MCherenkov::Dump () const {
    
  std::cout << *this
	    << std::endl;
}

 
std::string MCherenkov::GetParticleName () const {

  return std::string ("cherenkov");
}

