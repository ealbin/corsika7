#ifndef _INCLUDE_CORSIKAIO_TPARTICLE_
#define _INCLUDE_CORSIKAIO_TPARTICLE_

#include <TObject.h>

namespace crs {
  class MParticle;
  class GroundParticle;
  class MEventHeader;
}

namespace crsIO {

  /** 
      \class TParticle
      \brief One Particle-info.

      One Particle. This class also provides simple conversion algorithms. 

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TParticle.h,v 1.4 2007-10-19 08:00:45 rulrich Exp $
  */

  class TParticle : public TObject {
    
  public:
    TParticle (const crs::MParticle &right);
    TParticle (const crs::MEventHeader &header , 
	       const crs::GroundParticle &right);
    TParticle ();

    void Dump() const;

    bool IsMuonAdditionalInfo() const;
    double GetMuonProductionHeight () const; /// muon production height in cm
    
    int CorsikaID;
    int ParticleID;
    int ObservationLevel;
    int HadronicGeneration;
    
    double Px;
    double Py;
    double Pz;
    
    double x;
    double y;
    
    double Time;
    double Weight;
    
    ClassDef (TParticle, 5);
  };
}


#endif
