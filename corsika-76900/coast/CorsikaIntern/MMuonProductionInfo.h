#ifndef _INCLUDE_CORSIKA_MMUONPRODUCTIONINFO_
#define _INCLUDE_CORSIKA_MMUONPRODUCTIONINFO_

#include <crs/CorsikaTypes.h>
#include <crs/TParticleBlockEntry.h>
#include <crs/IParticleReadout.h>

#include <cmath>
#include <ostream>
#include <string>

/*
  #define IsParticle(p) (((p) > 0) && ((p) < 100000))
  #define IsNucleus(p) (((p) >= 100000) && ((p) < 9900000))
  #define IsCherenkov(p) ((p) >= 9900000)
*/


namespace crs {

  /** 
      \class MMuonProductionInfo
      \brief CORSIKA muonproductioninfo data

      One CORSIKA muonproductioninfo.

      \author Ralf Ulrich
      \date Sat Jun 18 00:14:28 CEST 2005
      \version $Id: MMuonProductionInfo.h 5116 2016-01-04 19:09:04Z darko $
  */

  class MMuonProductionInfo : public TParticleBlockEntry, public IParticleReadout {
	
	
  public:
    MMuonProductionInfo (const CREAL *data, bool thinned);
    MMuonProductionInfo (const TParticleBlockEntry &p);
	
    virtual void Dump () const;
	
  public:
    CREAL GetZ () const {return ValueAt (6);}	
    CREAL GetWeight () const {return (fThinned ? ValueAt (7) : 1);}
    virtual std::string GetParticleName () const;	

  private:
    inline CREAL ValueAt (int i) const {return fData [i];}

  };


  inline std::ostream &operator<< (std::ostream &o, 
				   const MMuonProductionInfo &p) {

    float P = p.GetPx ()*p.GetPx ();
    P += p.GetPy ()*p.GetPy ();
    P += p.GetPz ()*p.GetPz ();
    P = std::sqrt (P);

    o << " Muon produced: " /*<< p.GetParticleID () 
			      << " \'" << p.GetParticleName () << "\'"
			      << " obs.level: " << p.GetObservationLevel ()
			      << " gener.: " << p.GetHadronicGeneration ()*/
      << " E/GeV: " << p.GetKinEnergy ()
      << " theta/deg: " << std::acos (p.GetPz ()/P)*180./3.141
      << " at (x|y|z)/m: (" << p.GetX ()/100.
      << "|" << p.GetY ()/100.
      << "|" << p.GetZ ()/100. << ")"
      << " w: " << p.GetWeight ();

    return o;
	
  }

};    



#endif
