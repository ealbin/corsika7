/* $Id: TParticleBlockEntry.h 5849 2017-03-06 21:49:35Z rulrich $   */

#ifndef _INCLUDE_TPARTCLEDATA_
#define _INCLUDE_TPARTCLEDATA_

#include <crs/CorsikaTypes.h>

#include <string>
#include <iostream>

namespace crs {

  class MParticle;
  class MCherenkov;
  class MMuonProductionInfo;
  class PersistentParticle;

  typedef enum {
    eParticle,
    eNucleus,
    eCherenkov,
    eMuonProductionInfo,
    eEmpty,
    eUnknown
  } ParticleType;


  class TParticleBlockEntry {

    friend class TParticle;
    friend class TCherenkov;
    friend class TMuonProductionInfo;
    friend class PersistentParticle;

    friend std::ostream &operator<< (std::ostream &o, const TParticleBlockEntry &p);
    
  public:
    TParticleBlockEntry (const CREAL *data, bool thinn);
    TParticleBlockEntry (const TParticleBlockEntry& p);
    virtual ~TParticleBlockEntry () {}

    inline int GetParticleID () const {return (int)fData [0]/1000;}
    inline int GetParticleId () const {return (int)fData [0]/1000;}

    virtual void Dump () const;

    bool IsParticle () const;
    bool IsNucleus () const;
    bool IsCherenkov () const;
    bool IsMuonProductionInfo () const;
    bool IsEmpty () const;
    ParticleType GetType () const;

    virtual std::string GetParticleName () const;

  protected:
    const CREAL *fData;
    bool fThinned;

  };

  // std::ostream &operator<< (std::ostream &o, const TParticle &p);
  inline std::ostream &operator<< (std::ostream &o, const TParticleBlockEntry &p) {
    o << " entry: " << p.fData[0];    
    for (int i=1; i<7; ++i)
      o << " " << p.fData[i];
    if (p.fThinned)
      o << " " << (int)p.fData[7]; 
    return o;


  }

};


#endif
