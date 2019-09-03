/* $Id: TShower.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <TTree.h>

#include <crsIO/TShower.h>
using namespace crsIO;

#include <crs/MRunHeader.h>
#include <crs/MRunEnd.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>

// -----------------
ClassImp (TShower)

TShower::TShower () {

  Clear ();
}

TShower::~TShower () {

  delete fParticles;
  delete fCherenkov;
}

void TShower::Clear () {

  nPhotons = 0; 
  nElectrons = 0;
  nHadrons = 0;
  nMuons = 0;
  nPreshower = 0;
  nPhotonsWritten = 0;
  nElectronsWritten = 0;
  nHadronsWritten = 0;
  nMuonsWritten = 0;
  CPUtime = 0;

  fParticles = 0;
  fCherenkov = 0;
}

void TShower::AddRunHeader (const crs::MRunHeader &RunHeader) {
}


void TShower::AddRunEnd (const crs::MRunEnd &RunEnd) {
}

void TShower::AddEventHeader (const crs::MEventHeader &EventHeader) {

  EventID = (Int_t) EventHeader.GetEventNumber ();
  // EventHeader.GetParticleId (); 
  Energy = EventHeader.GetEnergy (); 
  StartingAltitude = EventHeader.GetStartingAltitude (); 
  FirstTarget = (Int_t) EventHeader.GetFirstTarget (); 
  FirstHeight = EventHeader.GetZFirst (); 
    
  //EventHeader.GetPx (); 
  //EventHeader.GetPy (); 
  //EventHeader.GetPz (); 
    
  Theta = EventHeader.GetTheta ();
  Phi = EventHeader.GetPhi (); 

  /*
    EventHeader.GetNRandomSequences ();
    
    int index;
    EventHeader.GetSeed (index); 
    EventHeader.GetInitialCallsMod (index);
    EventHeader.GetInitialCallsDiv (index); 
  */
}


void TShower::AddEventEnd (const crs::MEventEnd &EventEnd) {

  
  nPhotonsWritten = (int)EventEnd.GetPhotons (); 
  nElectronsWritten = (int)EventEnd.GetElectrons (); 
  nHadronsWritten = (int)EventEnd.GetHadrons (); 
  nMuonsWritten = (int)EventEnd.GetMuons (); 
  nParticlesWritten = (int)EventEnd.GetParticles (); 

  /*
  // NKG output
  EventEnd.GetLateral1X[21] (); 
  EventEnd.GetLateral1Y[21] (); 
  EventEnd.GetLateral1XY[21] ();
  EventEnd.GetLateral1YX[21] ();

  EventEnd.GetLateral2X[21] (); 
  EventEnd.GetLateral2Y[21] (); 
  EventEnd.GetLateral2XY[21] (); 
  EventEnd.GetLateral2YX[21] (); 

  EventEnd.GetElectronNumber[10] (); 
  EventEnd.GetAge[10] (); 
  EventEnd.GetDistances[10] (); 
  EventEnd.GetLocalAge1[10] (); 

  EventEnd.GetLevelHeightMass[10] (); 
  EventEnd.GetLevelHeightDistance[10] ();
  EventEnd.GetDistanceBinsAge[10] ();
  EventEnd.GetLocalAge2[10] (); 
  */

  GH_Nmax = EventEnd.GetNmax ();
  GH_t0 = EventEnd. GetX0 ();
  GH_tmax = EventEnd.GetXmax ();
  GH_a = EventEnd.GetLongi_a ();
  GH_b = EventEnd.GetLongi_b ();
  GH_c = EventEnd.GetLongi_c ();
  GH_Chi2 = EventEnd.GetLongi_Chi2 ();



  nPhotons = EventEnd.GetWeightedPhotons (); 
  nElectrons = EventEnd.GetWeightedElectrons ();
  nHadrons = EventEnd.GetWeightedHadrons ();
  nMuons = EventEnd.GetWeightedMuons ();

}


