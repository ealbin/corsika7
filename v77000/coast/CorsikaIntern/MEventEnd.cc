/* $Id: MEventEnd.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/MEventEnd.h>
using namespace crs;

MEventEnd::MEventEnd (const TSubBlock &right) {

  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
}

