/* $Id: MRunEnd.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/MRunEnd.h>
using namespace crs;

MRunEnd::MRunEnd (const TSubBlock &right) {

  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
}
