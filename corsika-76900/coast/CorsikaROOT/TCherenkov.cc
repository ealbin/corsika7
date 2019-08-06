/* $Id: TCherenkov.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crsIO/TCherenkov.h>
using namespace crsIO;

#include <crs/MCherenkov.h>


// -----------------
ClassImp (TCherenkov)

TCherenkov::TCherenkov (const crs::MCherenkov &ckov) {

  nPhotons = ckov.GetPhotonsInBunch ();
  x = ckov.GetX ();
  y = ckov.GetY ();

  u = ckov.GetU ();
  v = ckov.GetV ();
  Time = ckov.GetTime ();
  ProductionHeight = ckov.GetProductionHeight ();

  Weight = ckov.GetWeight ();
}


TCherenkov::TCherenkov () :
nPhotons (0),
  x (0),
  y (0),
  u (0),
  v (0),
  Time (0),
  ProductionHeight (0),
  Weight (0) {
}

