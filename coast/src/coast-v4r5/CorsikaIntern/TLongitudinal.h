#ifndef _INCLUDE_CORSIKA_TLONGITUDINAL_
#define _INCLUDE_CORSIKA_TLONGITUDINAL_


#include <crs/CorsikaTypes.h>


namespace crs {

  /** 
      \class TLongitudinal
      \brief CORSIKA longitudinal step data

      One CORSIKA longitudinal slice, containing particle numbers
      for a certain atmospheric depth.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TLongitudinal.h,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $
  */

  class TLongitudinal {

    
  public:
    TLongitudinal (const CREAL *data);
    
    void Dump ();
	
  public:
    CREAL GetSlantDepth () const {return fData [0];}

    unsigned long long GetNGammas () const {return (unsigned long long) fData [1];}
    unsigned long long GetNPositrons () const {return (unsigned long long) fData [2];}
    unsigned long long GetNElectrons () const {return (unsigned long long) fData [3];}
    unsigned long long GetNAntiMuons () const {return (unsigned long long) fData [4];}
    unsigned long long GetNMuons () const {return (unsigned long long) fData [5];}
    unsigned long long GetNHadrons () const {return (unsigned long long) fData [6];}
    unsigned long long GetNCharged () const {return (unsigned long long) fData [7];}
    unsigned long long GetNNuclei () const {return (unsigned long long) fData [8];}
    unsigned long long GetNCherenkov () const {return (unsigned long long) fData [9];}       

  private:
    //TLongitudinal (const TLongitudinal&);

  private:
    const CREAL *fData;
	
  };
};    



#endif
