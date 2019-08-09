/* $Id: TSubBlock.cc,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $   */

#include <crs/TSubBlock.h>
using namespace crs;

#include <string>
#include <sstream>
#include <iostream>

const std::string TSubBlock::fgRunStart   = "RUNH";
const std::string TSubBlock::fgRunEnd     = "RUNE";
const std::string TSubBlock::fgEventStart = "EVTH";
const std::string TSubBlock::fgEventEnd   = "EVTE";
const std::string TSubBlock::fgLongitudinal = "LONG";

// CORSIKA LIMITATION
const int TSubBlock::fgMaxObsLevels = 10;


TSubBlock::TSubBlock(const TSubBlock &sb) {

  fSubBlockData = sb.fSubBlockData;
  fType = sb.fType;
  fThinned = sb.fThinned;
}



TSubBlock::TSubBlock(const CREAL *data, bool thinned) {

  fSubBlockData = data;
  fThinned = thinned;
    
  FindBlockType();
}



void 
TSubBlock::FindBlockType() 
{    
  CCHAR *cBuf =(CCHAR*) fSubBlockData;
    
  std::ostringstream BlockOStr;
  BlockOStr << char(cBuf [0]) 
	    << char(cBuf [1]) 
	    << char(cBuf [2]) 
	    << char(cBuf [3]);
    
  std::string BlockStr = BlockOStr.str();

  int FirstEntry = int(fSubBlockData [0]);

  fType = eNODATA;

  if(BlockStr == fgRunStart) {

    fType = eRUNH;

  } else if(BlockStr == fgRunEnd) {

    fType = eRUNE;
	
  } else if(BlockStr == fgEventStart) {

    fType = eEVTH;
	
  } else if(BlockStr == fgEventEnd ) {
	
    fType = eEVTE;
	
  } else if(BlockStr == fgLongitudinal ) {
	
    fType = eLONG;
	
  } else if(FirstEntry > 0) {

    fType = ePARTDATA;

  } 

}


std::string 
TSubBlock::GetBlockTypeName() const 
{  
  std::string str = "";
  
  switch(fType) {
      case eUNKNOWN: str = "UNKOWN"; break;
      case eNODATA: str = "NODATA/CORRUPT"; break;
      case eRUNE: str = fgRunEnd; break;
      case eEVTE: str = fgEventEnd; break;
      case eLONG: str = fgLongitudinal; break;
      case ePARTDATA: str = "PARTICLEDATA"; break;
      case eEVTH: str = fgEventStart; break;
      case eRUNH: str = fgRunStart; break;
      default: str = "UNKNOWN"; break;
  }

  return str;

}
