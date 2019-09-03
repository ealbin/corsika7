#ifndef _M_CORSIKAREADER_H
#define _M_CORSIKAREADER_H

#include <crsRead/TCorsikaReader.h>

#include <crs/TSubBlock.h>

namespace crs {
  class TBlock;
  class TSubBlock;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
  class MRunEnd;
  class MParticleBlock;
  class MLongitudinalBlock;
  class MCherenkovBlock;
  class MParticle;
  class MCherenkov;
  class MMuonProductionInfo;
  class TLongitudinal;
};


#include <string>



namespace crsRead {
  
  class TSubBlock;
  
  /** 
      \class MCorsikaReader
      \brief Most straightforward implementation of a CORSIKA file reader
      
      An easy to use class for CORSIKA binary file reading.
      
      \author Ralf Ulrich
      \date Wed Jun 15 16:00:41 CEST 2005
      \version $Id: MCorsikaReader.h 5116 2016-01-04 19:09:04Z darko $ 
  */
  
  class MCorsikaReader : public TCorsikaReader {
    
  public:
    MCorsikaReader (const std::string &filename, int verbose=0);
    
    bool GetRun (crs::MRunHeader &run);
    bool GetShower (crs::MEventHeader &shower);
    bool GetData (crs::TSubBlock &data);
    bool GetShowerSummary (crs::MEventEnd &summary);
    
  private:
    virtual void HandleSubBlock (const crs::TSubBlock &CrsSubBlock);
    
  };
  
};

#endif
