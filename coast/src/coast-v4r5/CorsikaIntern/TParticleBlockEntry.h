/* $Id: TParticleBlockEntry.h,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $   */

#ifndef _INCLUDE_TPARTCLEDATA_
#define _INCLUDE_TPARTCLEDATA_

#include <crs/CorsikaTypes.h>

#include <string>


namespace crs {

  class MParticle;
  class MCherenkov;
  class MMuonProductionInfo;


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

  public:
    TParticleBlockEntry (const CREAL *data, bool thinn);
    TParticleBlockEntry (const TParticleBlockEntry& p);
    virtual ~TParticleBlockEntry () {}

    inline int GetParticleID () const {return (int)fData [0]/1000;}
	
    virtual void Dump () const {}

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

};


#endif
