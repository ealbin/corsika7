#ifndef _INCLUDE_CORSIKA_MPARTICLE_
#define _INCLUDE_CORSIKA_MPARTICLE_

#include <crs/CorsikaTypes.h>
#include <crs/TParticleBlockEntry.h>
#include <crs/IParticleReadout.h>

#include <cmath>
#include <ostream>
#include <string>


namespace crs {

  /** 
      \class TParticle
      \brief CORSIKA particle data
      
      One CORSIKA particle.
      
      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: MParticle.h,v 1.1.1.1 2007-07-31 07:00:51 rulrich Exp $
  */
  
  class MParticle : public TParticleBlockEntry, public IParticleReadout {
    
    
  public:
    
    MParticle (const CREAL *data, bool thinned);
    MParticle (const TParticleBlockEntry &p);
    
    virtual void Dump () const;
    
    
  public:
    
    CREAL GetWeight () const {return (fThinned ? ValueAt (7) : 1);}
    virtual std::string GetParticleName () const;
    
    
  private:
    
    inline CREAL ValueAt (int i) const {return fData [i];}
    
  };
  
  
  
  //std::ostream &operator<< (std::ostream &o, const TParticle &p);
  inline std::ostream &operator<< (std::ostream &o, const MParticle &p) {
    
    float P = p.GetPx ()*p.GetPx ();
    P += p.GetPy ()*p.GetPy ();
    P += p.GetPz ()*p.GetPz ();
    P = std::sqrt (P);
    
    o << " part: " << p.GetParticleID() 
      << " \'" << p.GetParticleName () << "\'"
      << " obs.level: " << p.GetObservationLevel ()
      << " gener.: " << p.GetHadronicGeneration ()
      << " E/GeV: " << p.GetKinEnergy ()
      << " theta/deg: " << std::acos (p.GetPz ()/P)*180./3.141
      << " (x|y)/m: (" << p.GetX ()/100.
      << "|" << p.GetY ()/100. << ")"
      << " T/ns: " << p.GetTime ()
      << " w: " << p.GetWeight ();
    
    return o;
  }
  
};    



#endif
