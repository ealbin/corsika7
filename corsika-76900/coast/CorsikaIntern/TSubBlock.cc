/* $Id: TSubBlock.cc 5924 2017-03-23 21:47:23Z rulrich $   */

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

  const double FirstEntry = double(fSubBlockData[0]);

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

  } else if(FirstEntry != 0) {

    fType = ePARTDATA;

  }

}


std::string
TSubBlock::GetBlockTypeName() const
{
  switch (fType) {
      case eUNKNOWN: return "UNKOWN"; break;
      case eNODATA: return "NODATA/CORRUPT"; break;
      case eRUNE: return fgRunEnd; break;
      case eEVTE: return fgEventEnd; break;
      case eLONG: return fgLongitudinal; break;
      case ePARTDATA: return "PARTICLEDATA"; break;
      case eEVTH: return fgEventStart; break;
      case eRUNH: return fgRunStart; break;
      default: break;
  }
  return "UNKNOWN";
}

std::string
TSubBlock::GetBlockTypeName(const SubBlockType type)
{
  switch (type) {
      case eUNKNOWN: return "UNKOWN"; break;
      case eNODATA: return "NODATA/CORRUPT"; break;
      case eRUNE: return fgRunEnd; break;
      case eEVTE: return fgEventEnd; break;
      case eLONG: return fgLongitudinal; break;
      case ePARTDATA: return "PARTICLEDATA"; break;
      case eEVTH: return fgEventStart; break;
      case eRUNH: return fgRunStart; break;
      default: break;
  }
  return "UNKNOWN";
}
