/* $Id: MRunHeader.cc,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $   */

#include <crs/MRunHeader.h>
using namespace crs;

MRunHeader::MRunHeader (const TSubBlock &right) {

    // ctor
    fSubBlockData = right.fSubBlockData;
    fType = right.fType;
    fThinned = right.fThinned;
}
