#ifndef _INCLUDE_CORSIKA_PERSISTENTPARTICLE_
#define _INCLUDE_CORSIKA_PERSISTENTPARTICLE_

#include <crs/CorsikaTypes.h>
#include <crs/TParticleBlockEntry.h>

#include <vector>


namespace crs {

  /**
      \class PersitentParticle
      \brief CORSIKA particle data

      One CORSIKA particle. In contrast to MParticle or TParticleBlockEntry,
      this class retains a copy of the particle data that does not get invalid
      while file is read.
      Use this class if for some reason you want to keep a copy of a particle
      during the particle loop.

      \author M. Unger
      \date Mon Jun  4 22:34:39 CEST 2012
      \version $Id: PersistentParticle.h 5116 2016-01-04 19:09:04Z darko $
  */

  class PersistentParticle : public TParticleBlockEntry {


  public:

    PersistentParticle(const CREAL* data, const bool thinned);
    PersistentParticle(const TParticleBlockEntry& p);
    PersistentParticle(const PersistentParticle& p);
    PersistentParticle& operator=(const PersistentParticle& p);

  private:

    std::vector<CREAL> fPersistentData;

  };
}

#endif
