/* $Id: PersistentParticle.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/TBlock.h>
#include <crs/PersistentParticle.h>
using namespace crs;

#include <iostream>

PersistentParticle::PersistentParticle(const float* data, const bool thinned) :
  TParticleBlockEntry(0, thinned)
{

  const int n = fThinned ?
    TBlock::fgNEntriesThinned : TBlock::fgNEntriesNotThinned;
  fPersistentData.resize(n);
  copy(data, data + n, fPersistentData.begin());
  fData = &fPersistentData.front();

}


PersistentParticle::PersistentParticle(const TParticleBlockEntry& p) :
  TParticleBlockEntry(0, p.fThinned)
{

  const int n = fThinned ?
    TBlock::fgNEntriesThinned : TBlock::fgNEntriesNotThinned;
  fPersistentData.resize(n);
  copy(p.fData, p.fData + n, fPersistentData.begin());
  fData = &fPersistentData.front();
}

PersistentParticle::PersistentParticle(const PersistentParticle& p) :
  TParticleBlockEntry(0, p.fThinned),
  fPersistentData(p.fPersistentData)
{
  fData = &fPersistentData.front();
}

PersistentParticle& PersistentParticle::operator=(const PersistentParticle& p)
{
  if (&p != this) {
    fPersistentData = p.fPersistentData;
    fThinned = p.fThinned;
    fData = &fPersistentData.front();
  }
  return *this;
}



