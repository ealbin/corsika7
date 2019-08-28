/* $Id: CInteraction.cc 1177 2009-09-24 20:32:08Z rulrich $   */

#include <crs/CInteraction.h>
#include <crs/CorsikaConsts.h>

#include <iostream>
#include <cstdlib>

using namespace crs;
using namespace std;


void CInteraction::Dump () const {
  cout << "Interaction projectileId: " << projId << ", targetId: " << targetId
       << ", energy: " << etot
       << ", x: " << x 
       << ", y: " << y 
       << ", z: " << z
       << ", kela:" << kela
       << endl;
}


