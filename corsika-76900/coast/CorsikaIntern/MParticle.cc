/* $Id: MParticle.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/CorsikaConsts.h>
#include <crs/MParticle.h>
using namespace crs;

#include <iostream>
#include <string>


MParticle::MParticle (const float *data, bool thinned) 
: TParticleBlockEntry (data, thinned) {
}


MParticle::MParticle (const TParticleBlockEntry &p) 
: TParticleBlockEntry (p) {
}

std::string MParticle::GetParticleName () const {
    
  int i = GetParticleId ();
  if (i>0 && i<=gNParticles)
    return gParticleName [i-1];

  return std::string ("undefinied");
}


void MParticle::Dump () const {
    
  std::cout << *this << std::endl;
}

