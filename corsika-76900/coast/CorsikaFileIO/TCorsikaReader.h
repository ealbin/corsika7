/**
   \file TCorsikaReader.h

   An interface to read CORSIKA files
   technicalities of Corsika file reading 

   \author Ralf Ulrich 
   \date Fri Jun 18 14:06:43 EDT 2004   
   \version $Id: TCorsikaReader.h 5116 2016-01-04 19:09:04Z darko $ 
*/


#ifndef TT_CORSIKAREADER_H
#define TT_CORSIKAREADER_H



#include <crsRead/TSubBlockIO.h>

#include <crs/TSubBlock.h>

#include <string>
#include <fstream>
#include <memory>



namespace crs {
  class TBlock;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
  class MRunEnd;
  class MParticleBlock;
  class MLongitudinalBlock;
  class MParticle;
  class MCherenkov;
  class MMuonProductionInfo;
  class TLongitudinal;
};


/**
   \namespace crsRead
   \brief Technicalisties for CORSIKA data file reading
   
   Inside this namespace all functionality for CORSIKA data file reading is 
   bundled. The class TCorsikaReader provides auto-detecion of the binary 
   file-format (thinned / not thinned, addidtional record-length info bytes).
   An example how to use TCorsikaReader is given in C2R/MCorsikaReader, that
   is used as an executable CORSIKA to ROOT converter. Using all the given 
   features of this package, the actual converter just needs a couple of lines
   of source code.
   
   \author Ralf Ulrich
   \date Sat Feb  5 16:54:58 CET 2005
   \version $Id: TCorsikaReader.h 5116 2016-01-04 19:09:04Z darko $
*/

namespace crsRead {


  /**
     \class TCorsikaReader
     \brief This class defines the interface to read CORSIKA files 

     This call of the Read() method will read the actual CORSIKA file
     and invoke the user methods Init, EventStart, EventEnd, RunStart,
     RunEnd, Particle, LongitudinalBlock, MuonAddInfo, CherenkovBunch and Exit
     to provide an interface for the handling of the CORSIKA data.

     \author Ralf Ulrich
     \date 05.2002
     \version $Id: TCorsikaReader.h 5116 2016-01-04 19:09:04Z darko $ 
  */

  class TCorsikaReader {
	
    typedef std::auto_ptr <std::ifstream> CorsikaFilePtr;

  public:
    TCorsikaReader ();
    TCorsikaReader (const std::string &filename, int verbose=-1);
    virtual ~TCorsikaReader ();

    // Do all the reading work
    bool Read ();
    //void KeepReading ();
    void InterruptReader ();
    void Rewind  ();



    // virtual Data Sub-Block handlers
    virtual void HandleSubBlock (const crs::TSubBlock &sb);
    virtual void HandleParticleBlock (const crs::MParticleBlock &p);
    virtual void HandleLongitudinalBlock(const crs::MLongitudinalBlock &s);
    virtual void HandleEventStart (const crs::MEventHeader& /*sb*/) {}
    virtual void HandleEventEnd (const crs::MEventEnd& /*sb*/) {}
    virtual void HandleRunStart (const crs::MRunHeader& /*sb*/) {}
    virtual void HandleRunEnd (const crs::MRunEnd& /*sb*/) {}
    virtual void HandleParticle (const crs::MParticle &p);
    virtual void HandleCherenkov (const crs::MCherenkov &c);
    virtual void HandleMuonProductionInfo (const crs::MMuonProductionInfo &m);
	
    // virtual user-event handlers
    virtual void Init () {} 
    virtual void Exit () {}


		
    // Setter inline functions
    inline void SetThinning (bool t) {fThinning=t;}
    inline void SetNBlockSizeInfo (int a) {fNBlockSizeInfo=a;}
    inline void SetVerboseLevel (int v) {fVerboseLevel=v;}
    inline void SetInputFileName (const std::string &fname) {Reset (); fInputFileName = fname;}
	
    // getter inline functions
    inline bool GetThinning () const {return fThinning;}
    inline int  GetNBlockSizeInfo () const {return fNBlockSizeInfo;}
    inline int  GetVerboseLevel () const {return fVerboseLevel;}
    inline const std::string & GetInputFileName () const {return fInputFileName;}
	
    //inline int  GetBlockCounter () const {return fBlockCounter;}
    //inline int  GetShowerCounter () const {return fShowerCounter;}
	
  private:
    void Reset ();
    CorsikaFilePtr TestFileBinaryStructure ();
    bool RetrieveData ();
	
  protected:
    /// flag that defines if the file uses thinning
    bool  fThinning;
	
    /// length of FORTRAN block size info
    int   fNBlockSizeInfo; 
	
    /// Verbose Level of the text ouput  (0-9)
    int   fVerboseLevel;     
	
    /// Name of the Input file
    std::string fInputFileName;  
	
    int   fLastShowerSubBlockCounter;  ///< Counter of SubBlocks
    int   fShowerCounter;              ///< Counter of Showers
    //crs::MParticleBlock::ParticleListConstIterator fParticleStart;
	
    bool fInterruption;
	
    int  fLastShowerPointer; // reference to the last shower read from the file

  protected:
    TSubBlockIO fSubBlock;
    CorsikaFilePtr fCorsFile;

  };

};

#endif
