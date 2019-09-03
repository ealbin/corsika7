#ifndef _INCLUDE_CORSIKAIO_TLONGITUDINAL_
#define _INCLUDE_CORSIKAIO_TLONGITUDINAL_

#include <TObject.h>

namespace crs {
  class TLongitudinal;
};

namespace crsIO {

  /** 
      \class TLongitudinal
      \brief One longitudinal step info.
      
      One longitudinal step.
      
      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TLongitudinal.h 5116 2016-01-04 19:09:04Z darko $
  */
  
  class TLongitudinal : public TObject {
    
  public:
    TLongitudinal (const crs::TLongitudinal &right);
    TLongitudinal ();
    
    float Depth;
    
    unsigned long long nGammas;
    unsigned long long nElectrons;
    unsigned long long nPositrons;
    unsigned long long nMuons;
    unsigned long long nAntiMuons;
    unsigned long long nHadrons;
    unsigned long long nCharged;
    unsigned long long nNuclei;
    unsigned long long nCherenkov;
    
    ClassDef (TLongitudinal, 2);
  };
};


#endif
