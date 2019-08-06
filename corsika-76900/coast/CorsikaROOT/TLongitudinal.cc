/* $Id: TLongitudinal.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crsIO/TLongitudinal.h>
using namespace crsIO;

#include <crs/TLongitudinal.h>


// -----------------
ClassImp (TLongitudinal)


TLongitudinal::TLongitudinal (const crs::TLongitudinal &right) {

  // CONVERT here ------------------
  Depth = right.GetSlantDepth ();
    
  nGammas = right.GetNGammas ();
  nElectrons = right.GetNElectrons ();
  nPositrons = right.GetNPositrons ();
  nMuons = right.GetNMuons ();
  nAntiMuons = right.GetNAntiMuons ();
  nHadrons = right.GetNHadrons ();
  nNuclei = right.GetNNuclei ();
  nCherenkov = right.GetNCherenkov ();

}


TLongitudinal::TLongitudinal () :
    
Depth (0),
  nGammas (0),
  nElectrons (0),
  nPositrons (0),
  nMuons (0),
  nAntiMuons (0),
  nHadrons (0),
  nNuclei (0),
  nCherenkov (0) {
}

