/* $Id: TSubBlockIO.cc,v 1.2 2007-10-19 11:13:40 rulrich Exp $   */


#include <crsRead/TSubBlockIO.h>
using namespace crsRead; 

#include <crs/TBlock.h>
#include <crs/TSubBlock.h>

#include <sstream>
#include <fstream>




TSubBlockIO::TSubBlockIO () 
: TSubBlock (), //0,false) {
  fBlockSize (0),
  fNBlockSizeInfo (0),
  fSubBlockCounter (0),
  fBlockSizeInfo (0),  
  fOwnBuffer (false) {
}


/** 
    AdditionalCharacters the number of extra CHARs after each 
    CORSIKA Block (21 SubBlocks)

*/
TSubBlockIO::TSubBlockIO (bool thinned, int AdditionalCharacters) 
: fNBlockSizeInfo (AdditionalCharacters),
  fSubBlockCounter (0),
  fBlockSizeInfo (new char [AdditionalCharacters]),
  fOwnBuffer (true) 
{  
  fBlockSize = (thinned ? 
		crs::TBlock::fgNParticles * crs::TBlock::fgNEntriesThinned :
		crs::TBlock::fgNParticles * crs::TBlock::fgNEntriesNotThinned);

  // and also init the TSubBlock
  fSubBlockData = new CREAL [fBlockSize];
  fThinned = thinned;
  //TSubBlock (new CREAL [fBlockSize], thinned);
}


TSubBlockIO::~TSubBlockIO () 
{
  delete [] fBlockSizeInfo;
  if (fSubBlockData && fOwnBuffer)
    delete [] fSubBlockData;
}


void 
TSubBlockIO::SetBuffer(const CREAL* buffer) 
{    
  // get rid of the allocated buffer
  if (fSubBlockData && fOwnBuffer)
    delete [] fSubBlockData;
    
  // class does not own the buffer, so it should not delete it
  fOwnBuffer = false;
  fSubBlockData = buffer;
}

void 
TSubBlockIO::SetBuffer(const crs::TSubBlock& sb)
{
  // get rid of the allocated buffer
  if (fSubBlockData && fOwnBuffer)
    delete [] fSubBlockData;
  
  // class does not own the buffer, so it should not delete it
  fOwnBuffer = false;
  fSubBlockData = sb.GetData();
}


void 
TSubBlockIO::InitSubBlockSize (bool thinning, int AdditionalCharacters) 
{   
  if (fOwnBuffer)
    delete [] fSubBlockData;

  delete [] fBlockSizeInfo;

  fThinned = thinning;

  fNBlockSizeInfo = AdditionalCharacters;
  fBlockSizeInfo = new char [AdditionalCharacters];
  fSubBlockCounter = 0;
  fBlockSize = (thinning ? 
		crs::TBlock::fgNParticles * crs::TBlock::fgNEntriesThinned :
		crs::TBlock::fgNParticles * crs::TBlock::fgNEntriesNotThinned);
  fOwnBuffer = true;
  fSubBlockData = new CREAL [fBlockSize];


  /*
    int BlockSize = 21 * fBlockSize * 4; 
    
    for (int i=0; i<n; i++) {	

    if (i<4) {
    ((char*)fBlockSizeInfo) [i] = ((char*)&BlockSize) [i];
    } else {
    ((char*)fBlockSizeInfo) [i] = 0;
    }

    }
  */
}



int 
TSubBlockIO::FindRUNH () 
{    
  std::string strRUNH = "RUNH";

  int posRUNH = -1;
    
  char *TestChar = (char*) fSubBlockData;
  for (int i=0; i<fBlockSize; i++) {

    std::ostringstream TestStr;
    TestStr << TestChar [i+0]
	    << TestChar [i+1]
	    << TestChar [i+2]
	    << TestChar [i+3];
	
    if (TestStr.str () == strRUNH) {
      posRUNH = i;
      break;
    }

  }
    
  return posRUNH;
}

int
TSubBlockIO::FinishBlock(std::ofstream& s) 
{
  // and also init the TSubBlock
  if (fOwnBuffer)
    delete [] fSubBlockData;

  fOwnBuffer = true;
  CREAL* subBlockData = new CREAL[fBlockSize];
  for (int i=0; i<fBlockSize; ++i) 
    subBlockData[i] = 0;
  fSubBlockData = subBlockData;
  
  int c = 0;
  while (fSubBlockCounter < crs::TBlock::fgNSubBlocks) {
    s << *this;
    ++c;
  }
  return c;
}
