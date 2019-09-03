#include <crsRead/TCorsikaReader.h>

namespace crs {
  class TSubBlock;
};

namespace crs2r {

  class TC2R;

  /** 
      \class MCorsikaReader
      \brief Implementation of a CORSIKA file reader writing ROOT files
      \ingroup CONVERTER

      This is an implementation of the crsRead::TCorsikaReader class 
      to read CORSIKA binary data files and write them as ROOT object
      files. It is an example of a CORSIKA to ROOT conversion programm.

      \author Ralf Ulrich
      \date Wed Jun 15 15:44:59 CEST 2005
      \version $Id: MCorsikaReader.h 5116 2016-01-04 19:09:04Z darko $
  */

  class MCorsikaReader : public crsRead::TCorsikaReader {
      
  public:
    MCorsikaReader(const std::string &fname, int verbose=0);
    
    virtual void Init();
    virtual void HandleSubBlock(const crs::TSubBlock &sb);
    virtual void Exit();

  private:
    TC2R *fC2R;
	
  };

};
