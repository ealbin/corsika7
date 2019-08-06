/* $Id: CInteraction.h 5116 2016-01-04 19:09:04Z darko $   */


#ifndef _INCLUDE_CInteraction_H__
#define _INCLUDE_CInteraction_H__

#include <crs/CorsikaTypes.h>


namespace crs {
  
  // ***************************************************
  struct CInteraction {	
  public:
  CInteraction() :
      x(0),
      y(0),
      z(0),
      etot(0),
      sigma(0),
      kela(0),
      projId(0),
      targetId(0),
      time(0) {
    }
    void Dump () const;

  public:

    CDOUBLEPRECISION x;
    CDOUBLEPRECISION y;
    CDOUBLEPRECISION z;
    CDOUBLEPRECISION etot;
    CDOUBLEPRECISION sigma;
    CDOUBLEPRECISION kela;
    CINT projId;
    CINT targetId;	
    CDOUBLEPRECISION time;
    
    int GetProjectileId() const { return projId; }
    int GetTargetId() const { return targetId; }    
    
  };
};


#endif
