/* $Id: TSubBlockIO.cc 5116 2016-01-04 19:09:04Z darko $   */


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
TSubBlockIO::TSubBlockIO (bool thinned, unsigned int AdditionalCharacters) 
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

  InitBlockSizeInfo();
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

/* to overtake foreign data buffer
void
TSubBlockIO::Grab() 
{
  if (fOwnBuffer)
    return;
  
  CREAL* data = fSubBlockData;
  
  fOwnBuffer = true;
  fSubBlockData = new CREAL[fBlockSize];
  
  for (int i=0; i<fBlockSize; ++i) 
    subBlockData[i] = data[i];  
}
*/

void 
TSubBlockIO::InitSubBlockSize (bool thinning, unsigned int AdditionalCharacters) 
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

  InitBlockSizeInfo();
}


void
TSubBlockIO::InitBlockSizeInfo()
{
  union {
    int BlockSize;
    float BlockSizeFloat;
  } record_length_val;
  
  record_length_val.BlockSize = 21 * fBlockSize * 4;
    
  
  for (unsigned int i=0; i<fNBlockSizeInfo; i++) {	
    
    if (i<sizeof(float)) {
      ((char*)fBlockSizeInfo) [i] = ((char*)&record_length_val) [i];
    } else {
      ((char*)fBlockSizeInfo) [i] = 0;
    }    
  }
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
  while ((fSubBlockCounter % crs::TBlock::fgNSubBlocks)) {
    s << *this;
    ++c;
  }
  return c;
}
