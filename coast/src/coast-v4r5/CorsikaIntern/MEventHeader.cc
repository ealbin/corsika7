/* $Id: MEventHeader.cc,v 1.1.1.1 2007-07-31 07:00:52 rulrich Exp $   */

#include <crs/MEventHeader.h>
using namespace crs;

MEventHeader::MEventHeader (const TSubBlock &right) {

  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
}
