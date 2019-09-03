/* $Id: MCherenkov.cc 5849 2017-03-06 21:49:35Z rulrich $   */

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

CREAL MCherenkov::GetPhotonsInBunch () const
{
  const double v = fData[0];
  if (v >= 99e5) {
    return (v-99e5-1.0)/10.;
  }
  return v;
}

void MCherenkov::Dump () const {
  std::cout << *this << std::endl;
}

 
std::string MCherenkov::GetParticleName () const {
  return std::string ("cherenkov");
}

