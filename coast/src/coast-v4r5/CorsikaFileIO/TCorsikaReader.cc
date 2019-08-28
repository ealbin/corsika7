/********************************************************************/
/* TCorsikaReader.C                                                */
/* technicalities of Corsika file reading                           */
/* Ralf Ulrich Fri Jun 18 14:06:43 EDT 2004                         */
/* $Id: TCorsikaReader.cc,v 1.1.1.1 2007-07-31 06:58:57 rulrich Exp $   */
/********************************************************************/

#include <crsRead/TCorsikaReader.h>
using namespace crsRead;

#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MRunEnd.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>
#include <crs/MCherenkov.h>
#include <crs/MMuonProductionInfo.h>

#include <string>
#include <iostream>
#include <fstream>
using namespace std;


TCorsikaReader::TCorsikaReader() 
  : fVerboseLevel(1) {

  fInputFileName = "";  
  Reset();
}




/// fname File name of the CORSIKA data file
TCorsikaReader::TCorsikaReader(const std::string &DataFileName, 
			       int verbose) 
{
  fVerboseLevel = verbose;

  if(fVerboseLevel>10) {
    cout << " entering TCorsikaReader::TCorsikaReader " << endl;
  }
  
  if(fVerboseLevel>0) {
    cout << " reading CORSIKA file: \"" << DataFileName << "\"" << endl;
  }
  
  Reset();

  fInputFileName = DataFileName;
}




TCorsikaReader::~TCorsikaReader() 
{
}




void 
TCorsikaReader::Reset() 
{
  if(fVerboseLevel>10) {
    cout << " entering TCorsikaReader::Reset " << endl;
  }

  // reset counters
  fShowerCounter = 0;              
  fLastShowerSubBlockCounter = 0;  
 
  // reset pointer, references
  //fParticleStart = 0;
  fLastShowerPointer = 0;
 
  // reset thinning mode
  fThinning = false;

  // reset data structure
  fNBlockSizeInfo = 4;        // bytes(at least for LINUX file structure)
 
  fInterruption = false;

  // reset smart-pointers
  fCorsFile.reset ();
}





void 
TCorsikaReader::Rewind () 
{    
  if (fVerboseLevel>10) {
    cout << " entering TCorsikaReader::Rewind " << endl;
  }

  if (fCorsFile.get()) {
    
    // rewind file pointer
    fCorsFile->seekg(fLastShowerPointer);  
	
    // rewind sub block counter
    fSubBlock.SetSubBlockCounter (fLastShowerSubBlockCounter);
	
    // rewind shower counter
    fShowerCounter--;
	
    // reset particle data reader
    //fParticleStart = 0;
    fInterruption = false;
  }
    
}



/**
   \fn int TCorsikaReader::Read()
   Reads a CORSIKA data file and invokes the USER methods to handle the 
   different SubBlocks
   First shower is fStartShower and last shower is fEndShower
*/
bool 
TCorsikaReader::Read() 
{    
  if (fVerboseLevel>10) {
    cout << " entering TCorsikaReader::Read " << endl;
  }

  fInterruption = false;
    
  // look if file is already open
  if (NULL == fCorsFile.get()) {
    
    Reset();

    // check the CORSIKA file (thinned ... )
    fCorsFile = TestFileBinaryStructure();
    
    //
    // Call User Init Function
    //
    Init();
    if (fVerboseLevel > 3) {
      std::cout<<" USER INIT"<<std::endl;
    }
  }
    
  // now read the actual CORSIKA file, or go on, when interupted
  bool status = RetrieveData(); 

  if (false==fInterruption) {

    //
    // Call user Exit function
    //
    Exit();
    if (fVerboseLevel > 3) {
      std::cout<<" USER EXIT"<<std::endl;
    }

  }

  return status;
}




void 
TCorsikaReader::InterruptReader() 
{  
  if(fVerboseLevel>10) {
    cout << " entering TCorsikaReader::InterruptReader " << endl;
  }

  fInterruption = true;
}




/// \fn int TCorsikaReader::TestFileBinaryStructure()
/// Is automatically invoked after opening a CORSIKA data file
/// checks if thinned or not thinned
TCorsikaReader::CorsikaFilePtr 
TCorsikaReader::TestFileBinaryStructure() 
{
  if(fVerboseLevel>10) {
    cout << " entering TCorsikaReader::TestFileBinaryStructure " << endl;
  }
  
  //
  // OPEN INPUT FILE 
  //
  CorsikaFilePtr CorsFile(new std::ifstream(fInputFileName.c_str(), 
					     std::ios::in|std::ios::binary));
    
  if (!(CorsFile->good())) {
    std::cerr << "COULD NOT OPEN FILE " << fInputFileName << std::endl;
    return CorsikaFilePtr();
  }
    
    
  // 
  // Read first Sub-block and check for heading-bytes -> skip
  //
    
  //CorsFile->seekg(fSkipFirst , std::ios::beg);
  TSubBlockIO RUNHBlock(false, 0); // this is not a 'real' sub-block !!!!!
 (*CorsFile) >> RUNHBlock;
  fNBlockSizeInfo = RUNHBlock.FindRUNH();
    
  if(-1==fNBlockSizeInfo) {
    std::cerr << "COULD NOT FIND RUNH in " << fInputFileName << std::endl;
    return CorsikaFilePtr();
  }
    
  //
  // CHECK for THINNING
  //
  if(fVerboseLevel > 1) {
    std::cout << " TRY TO DETECT CORSIKA THINNING MODE" << std::endl;
  }
    
  // use blocksize without thinning
  TSubBlockIO SmallTestBlock(false, fNBlockSizeInfo);
    
  CorsFile->seekg(0, std::ios::beg);
 (*CorsFile) >> SmallTestBlock;
 (*CorsFile) >> SmallTestBlock;
    
  bool SmallEVTH = 
    SmallTestBlock.GetBlockType() == crs::TSubBlock::eEVTH;
    
  if(SmallEVTH) {
	
    fThinning = false;
	
  } else {
	
    // use blocksize with thinning
    TSubBlockIO LargeTestBlock(true, fNBlockSizeInfo);
	
    CorsFile->seekg(0 , std::ios::beg);
   (*CorsFile) >> LargeTestBlock;
   (*CorsFile) >> LargeTestBlock;
	
    bool LargeEVTH = 
      LargeTestBlock.GetBlockType() == crs::TSubBlock::eEVTH;
	
    if(LargeEVTH) {
	    
      fThinning = true;
	    
    } else {
	    
      std::cerr << "COULD NOT AUTOMATICALLY DETECT IF THE "
		<< " CORSIKA FILE "
		<< fInputFileName << " "
		<< "WAS THINNED OR NOT !" << std::endl;
      return CorsikaFilePtr();
    }
  }

  if (fVerboseLevel>1) {
    std::cout << " AUTO-DETECT of thinning mode result: "
	      <<(fThinning ? "THINNED" : "NOT THINNED") << std::endl;
  }


  // spool to the beginning of the CORSIKA file
  CorsFile->seekg(0, std::ios::beg);  
  
  // prepare a proper corsika subblock, suitable to store all 
  // incomming data
  fSubBlock.InitSubBlockSize(fThinning, fNBlockSizeInfo);
    
  // and return it as an open stream
  return CorsFile;
}
 

 

/**
   \fn int TCorsikaReader::RetrieveData()
   Fetches the data of one single shower
*/
bool 
TCorsikaReader::RetrieveData() 
{    
  if(fVerboseLevel>10) {
    cout << " entering TCorsikaReader::RetrieveData " << endl;
  }
  
  if (NULL==fCorsFile.get()) {
    std::cerr << " file pointer NULL while reading from " 
	      << fInputFileName 
	      << std::endl;
    return false;
  } 
    
  if (!(fCorsFile->good())) {
    // reached file end
    return false;
  }
    
    
  //
  // receiving SubBlocks until shower is complete
  //
  while (fCorsFile->good() &&
	  !fInterruption) {
	
   (*fCorsFile) >> fSubBlock;    
	
    HandleSubBlock(fSubBlock);
	
  } // while
    
  if (!(fCorsFile->good())) {
    // reached file end
    return false;
  }
    
  return true;
}




void 
TCorsikaReader::HandleSubBlock(const crs::TSubBlock &CrsSubBlock) 
{
  if (fVerboseLevel>10) {
    cout << " entering TCorsikaReader::HandleSubBlock " << CrsSubBlock.GetBlockTypeName() << endl;
  }

  int CurrentBlockPointer = fCorsFile->tellg();
  int LastShowerSubBlockCounter = fSubBlock.GetSubBlockCounter ();

  switch (CrsSubBlock.GetBlockType()) {
	
  case crs::TSubBlock::eRUNH:                             // RUN HEADER

    HandleRunStart(CrsSubBlock);
    if (fVerboseLevel>3) {
      std::cout << " USER RunStart" << std::endl;
    }
    break;
		
		

  case crs::TSubBlock::eRUNE:       	                // Run End

    HandleRunEnd(CrsSubBlock);
    if (fVerboseLevel>3) {
      std::cout << " USER RunEnd" << std::endl;
    }	       
    break;
	    
		
	    
  case crs::TSubBlock::eEVTH:      	                // Event Header
	    
    fShowerCounter++;
	    
    fLastShowerPointer = CurrentBlockPointer;
    fLastShowerSubBlockCounter = LastShowerSubBlockCounter;
	    
    HandleEventStart(CrsSubBlock);
    if (fVerboseLevel>3) {
      std::cout<<" USER EventStart"<<std::endl;
    }              
    break;
	    
	    
	    
  case crs::TSubBlock::eEVTE:	                        // Event End
	    
    HandleEventEnd(CrsSubBlock);
    if (fVerboseLevel>3) {
      std::cout<<" USER EventEnd"<<std::endl;
    }	       
    break;
	    
	    

	    
  case crs::TSubBlock::ePARTDATA: 	              // PARTICLE Block
	    
    HandleParticleBlock(CrsSubBlock);
    if (fVerboseLevel>3) {
      std::cout<<" USER particle block"<<std::endl;
    }	       
    break;
		
	    
	    
  case crs::TSubBlock::eLONG:	              // Longitudinal Sub Block

    HandleLongitudinalBlock(CrsSubBlock);
    if (fVerboseLevel>3) {
      std::cout << " USER Longitudinal" << std::endl;
    }	       
    break;


  case crs::TSubBlock::eNODATA:
    if (fVerboseLevel>3) {
      std::cout << " empty data block !" << std::endl;
    }	       
    break;
	    
  case crs::TSubBlock::eUNKNOWN: // to prevent compiler warning
    break;

  } // switch
    
}



void TCorsikaReader::HandleLongitudinalBlock(const crs::MLongitudinalBlock& /*b*/){
}






void TCorsikaReader::HandleParticleBlock(const crs::MParticleBlock &Data) {
	    
  //
  // 39 particles are coming through buffer,
  //
    
  crs::MParticleBlock::ParticleListConstIterator iEntry;
  for(iEntry = Data.FirstParticle();
       iEntry != Data.LastParticle();
       ++iEntry) {
	
    switch(iEntry->GetType()) {

    case crs::eParticle:
    case crs::eNucleus:
      HandleParticle(*iEntry);
      break;
		
    case crs::eCherenkov:
      HandleCherenkov(*iEntry);
      break;
		
    case crs::eMuonProductionInfo:
      HandleMuonProductionInfo(*iEntry);
      break;

    case crs::eEmpty:
      break;

    case crs::eUnknown:
      break;

    }
	
  } // for-loop
    
    // start with fresh particle block, this one is empty
    /*
      if (iParticle == Data.LastParticle()) {
      fNewParticleBlock = true;
      }
    */

}




void 
TCorsikaReader::HandleParticle(const crs::MParticle &p) 
{
  if (fVerboseLevel>4) {	   
    std::cout << p << std::endl;
  }    
}

void 
TCorsikaReader::HandleCherenkov(const crs::MCherenkov& /*c*/) 
{
}

void
TCorsikaReader::HandleMuonProductionInfo(const crs::MMuonProductionInfo& /*m*/) 
{
}

