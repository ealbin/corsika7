#ifndef _INCLUDE_CORSIKA_TSUBBLOCH_H_
#define _INCLUDE_CORSIKA_TSUBBLOCH_H_

#include <crs/CorsikaTypes.h>

#include <string>

namespace crsRead {
  class TSubBlock;
}


namespace crs {

  class MRunHeader;
  class MRunEnd;
  class MEventHeader;
  class MEventEnd;
  class MParticleBlock;
  class MLongitudinalBlock;

  /** 
      \class TSubBlock
      \brief One CORSIKA sub-block, the container for all data

      All data is comming in sub-blocks of different types, that are
      encoded in the first REAL*4 of a sub-block. A sub-block is 
      able to contain 39 particles. The number of REAL*4 per particle is
      7 for normal CORSIKA but 8 for thinned CORSIKA.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TSubBlock.h,v 1.2 2007-10-19 11:13:40 rulrich Exp $
  */

  class TSubBlock {

    //friend class crsRead::TSubBlock;
	
    friend class MRunHeader;
    friend class MRunEnd;
    friend class MEventHeader;
    friend class MEventEnd;
    friend class MParticleBlock;
    friend class MLongitudinalBlock;
	
  public:

    typedef enum { 
      eRUNH, 	  // 0 Run Header
      eEVTH, 	  // 1 Event Header
      ePARTDATA,  // 2 Particle Block
      eLONG,      // 3 Longitudinal Block
      eEVTE, 	  // 4 Event End
      eRUNE, 	  // 5 Run End
      eNODATA,    // 6 empty SubBlock
      eUNKNOWN    // 7 not checked
    } SubBlockType;
	
	
    TSubBlock () : fSubBlockData (0), fType (eUNKNOWN), fThinned (false) {}
    TSubBlock (const TSubBlock &);
    TSubBlock (const CREAL *subblockdata, bool thinned);
    virtual ~TSubBlock () {}

    SubBlockType GetBlockType () const {return fType;}
    std::string GetBlockTypeName () const;
    bool IsThinned () const {return fThinned;}

    const CREAL* GetData () const {return fSubBlockData;}

  private:
    //void operator= (const TSubBlock &);

  public:
    void FindBlockType ();
    
  protected:
    const CREAL *fSubBlockData;
    SubBlockType fType;
    bool fThinned;
    
    static const std::string fgRunStart;
    static const std::string fgRunEnd;
    static const std::string fgEventStart;
    static const std::string fgEventEnd;
    static const std::string fgLongitudinal;
	
    static const int fgMaxObsLevels;
  };

};

#endif
