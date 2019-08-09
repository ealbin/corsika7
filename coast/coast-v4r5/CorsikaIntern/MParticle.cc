/* $Id: MParticle.cc,v 1.2 2007-10-19 11:13:40 rulrich Exp $   */

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
    
  int i = GetParticleID ();
  if (i>0 && i<=gNParticles)
    return gParticleName [i-1];

  return std::string ("undefinied");
}


void MParticle::Dump () const {
    
  std::cout << *this << std::endl;
}

