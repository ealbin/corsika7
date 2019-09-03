/* $Id: MEventHeader.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/MEventHeader.h>
using namespace crs;

MEventHeader::MEventHeader (const TSubBlock &right) {

  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
}
