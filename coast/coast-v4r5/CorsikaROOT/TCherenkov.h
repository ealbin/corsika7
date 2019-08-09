#ifndef _INCLUDE_CORSIAKIO_TCHERENKOV_
#define _INCLUDE_CORSIAKIO_TCHERENKOV_

#include <TObject.h>

namespace crs {
  class MCherenkov;
};


namespace crsIO {

  /** 
      \class TCherenkov
      \brief One Cherenkov bunch info.

      One Cherenkov bunch. 

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TCherenkov.h,v 1.1.1.1 2007-07-31 07:00:32 rulrich Exp $
  */

  class TCherenkov : public TObject {

  public:
    TCherenkov (const crs::MCherenkov &ckov);
    TCherenkov ();

    float nPhotons; // Number of photons in this bunch
	
    float x; // cm
    float y; // cm
	
    float u; // direction cosine to x axis
    float v; // direction cosine to y axis
	
    float Time; // time since first interction   
	
    float ProductionHeight; // production height

    float Weight; // bunch weight
	
    ClassDef (TCherenkov, 1)
  };
};

#endif

