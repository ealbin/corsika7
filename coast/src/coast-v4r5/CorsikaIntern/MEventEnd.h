#ifndef _INCLUDE_CORSIKA_MEVENTEND_
#define _INCLUDE_CORSIKA_MEVENTEND_

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>

namespace crs {

  /** 
      \class MEventEnd
      \brief CORSIKA event end (EVTE) sub-block

      This class knows how to access the data stored in this sub-block type.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: MEventEnd.h,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $
  */

  class MEventEnd : public TSubBlock {
	
  public:
    MEventEnd () {}
    MEventEnd (const TSubBlock &right);
    virtual ~MEventEnd () {}
	
  public:
	
    CREAL GetEventNumber () const {return fSubBlockData [1];}
	
    CREAL GetPhotons () const {return fSubBlockData [2];}
    CREAL GetElectrons () const {return fSubBlockData [3];}
    CREAL GetHadrons () const {return fSubBlockData [4];}
    CREAL GetMuons () const {return fSubBlockData [5];}
    CREAL GetParticles () const {return fSubBlockData [6];}
	

    // NKG output
    /*
      CREAL GetLateral1X[21] () const {return fSubBlockData [7];}
      CREAL GetLateral1Y[21] () const {return fSubBlockData [];} // cm^-2
      CREAL GetLateral1XY[21] () const {return fSubBlockData [];}
      CREAL GetLateral1YX[21] () const {return fSubBlockData [];} // cm^-2

      CREAL GetLateral2X[21] () const {return fSubBlockData [];}
      CREAL GetLateral2Y[21] () const {return fSubBlockData [];} // cm^-2
      CREAL GetLateral2XY[21] () const {return fSubBlockData [];}
      CREAL GetLateral2YX[21] () const {return fSubBlockData [];} // cm^-2

      CREAL GetElectronNumber[10] () const {return fSubBlockData [];} // in steps of 100 g cm^-2
      CREAL GetAge[10] () const {return fSubBlockData [];}            // in steps of 100 g cm^-2
      CREAL GetDistances[10] () const {return fSubBlockData [];} // cm
      CREAL GetLocalAge1[10] () const {return fSubBlockData [];}

      CREAL GetLevelHeightMass[10] () const {return fSubBlockData [];} // g cm^-2
      CREAL GetLevelHeightDistance[10] () const {return fSubBlockData [];} // cm
      CREAL GetDistanceBinsAge[10] () const {return fSubBlockData [];} // cm
      CREAL GetLocalAge2[10] () const {return fSubBlockData [];}
    */

    /// Longitudinal distribution (see manual p. 64)
    CREAL GetNmax () const {return fSubBlockData [255];}
    CREAL GetX0 () const {return fSubBlockData [256];}
    CREAL GetXmax () const {return fSubBlockData [257];}
    CREAL GetLongi_a () const {return fSubBlockData [258];}
    CREAL GetLongi_b () const {return fSubBlockData [259];}
    CREAL GetLongi_c () const {return fSubBlockData [260];}
    CREAL GetLongi_Chi2 () const {return fSubBlockData [261];}
    
    /// Added according to the CORSIKA manual
    CREAL GetWeightedPhotons () const {return fSubBlockData [262];}
    CREAL GetWeightedElectrons () const {return fSubBlockData [263];}
    CREAL GetWeightedHadrons () const {return fSubBlockData [264];} 
    CREAL GetWeightedMuons () const {return fSubBlockData [265];}
    CREAL GetNPreshower () const {return fSubBlockData [266];}

  };
};

#endif
