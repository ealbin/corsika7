#ifndef _INCLUDE_CORSIKAIO_TC2R_
#define _INCLUDE_CORSIKAIO_TC2R_

#include <vector>
#include <string>

#include <crs/CorsikaTypes.h>

class TFile;
class TTree;
class TBranch;
class TClonesArray;

namespace crsIO {
  class TRun;
  class TShower;
  class TParticle;
};

namespace crs {
  class TBlock;
  class TSubBlock;
  class GroundParticle;
  class MParticleBlock;
  class MLongitudinalBlock;
  class MEventHeader;
};

/**
   \defgroup CONVERTER Corsika to ROOT Converter
   
   Provides the functionality to convert CORSIKA (block/sub-block) data
   structure into ROOT objects.

   \author Ralf Ulrich
   \date Sat Feb  5 17:28:09 CET 2005
   \version $Id: TC2R.h,v 1.3 2007-10-18 14:33:47 pierog Exp $
*/


/** 
    \namespace crs2r
    \brief tools for CORSIKA to ROOT conversion
    \ingroup CONVERTER
    
    \author Ralf Ulrich
    \date Wed Jun 15 15:41:54 CEST 2005
    \version $Id: TC2R.h,v 1.3 2007-10-18 14:33:47 pierog Exp $
*/

namespace crs2r {
    
  /** 
      \class TC2R
      \brief CORSIKA to ROOT converter
      \ingroup CONVERTER

      This class is able to convert a standard CORSIKA block/sub-block data
      structure into an object-oriented ROOT structure. All data is saved
      into TTrees.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TC2R.h,v 1.3 2007-10-18 14:33:47 pierog Exp $
  */

  class TC2R {


  public:
    TC2R ();
    ~TC2R ();
    
    void Open (std::string filename,
	       bool thinned);

    void Write (const CREAL *data);
    void Write (const crs::TBlock &Block);
    void Write (const crs::TSubBlock &SubBlock);
    void Write (const crs::MEventHeader &header , const crs::GroundParticle& particle);
	
    void Close ();

    bool IsThinned () const {return fThinned;}
    
  protected:
    void AddParticleBlock (const crs::MParticleBlock &SubBlock);
    void AddLongitudinalBlock (const crs::MLongitudinalBlock &SubBlock);
	
	
  private:
    
    bool fThinned;
    
    TFile *fFile;
    
    TTree *fShower;
    TTree *fRun;
    
    crsIO::TRun *fCurrentRun;
    crsIO::TShower *fCurrentShower;

    TClonesArray *fParticles;
    TClonesArray *fCherenkov;
    TClonesArray *fLongitudinal;
	
    int fNParticles;
    int fNLong;
    int fNCherenkov;
	
  };

};
    
#endif

