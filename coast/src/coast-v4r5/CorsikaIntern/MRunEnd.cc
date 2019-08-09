/* $Id: MRunEnd.cc,v 1.1.1.1 2007-07-31 07:00:51 rulrich Exp $   */

#include <crs/MRunEnd.h>
using namespace crs;

MRunEnd::MRunEnd (const TSubBlock &right) {

  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
}
