#ifndef _INCLUDE_CORSIKA_MCHERENKOV_
#define _INCLUDE_CORSIKA_MCHERENKOV_


#include <crs/CorsikaTypes.h>
#include <crs/TParticleBlockEntry.h>

#include <ostream>


namespace crs {

  /** 
      \class MCherenkov
      \brief CORSIKA Cherenkov data

      One CORSIKA Cherenkov-bunch data object.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: MCherenkov.h 5849 2017-03-06 21:49:35Z rulrich $
  */

  class MCherenkov : public TParticleBlockEntry {

  public:
    MCherenkov (const CREAL *data, bool fThinned);
    MCherenkov (const TParticleBlockEntry &p);

    virtual void Dump () const;
	
  public:
    CREAL GetPhotonsInBunch () const;
	
    CREAL GetX () const {return fData [1];} 
    CREAL GetY () const {return fData [2];}

    CREAL GetU () const {return fData [3];}
    CREAL GetV () const {return fData [4];}

    CREAL GetTime () const {return fData [5];}
	
    CREAL GetProductionHeight () const {return fData [6];}
	
    CREAL GetWeight () const {return (fThinned ? fData [7] : 1);}


  public:
    virtual std::string GetParticleName () const;
  };


  inline std::ostream &operator<< (std::ostream &o, const MCherenkov &p) {

    o << " ckov: " << p.GetPhotonsInBunch ()
      << " x: " << p.GetX ()
      << " y: " << p.GetY ()
      << " u: " << p.GetU () 
      << " v: " << p.GetV ()
      << " t: " << p.GetTime () 
      << " w: " << p.GetWeight ();
	
    return o;
	
  }

};    

#endif
