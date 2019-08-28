#include <crsIO/TParticle.h>
using namespace crsIO;

#include <crs/MParticle.h>
#include <crs/CParticle.h>
#include <crs/CorsikaConsts.h>
#include <crs/MEventHeader.h>

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

// -----------------
ClassImp(TParticle);


TParticle::TParticle (const crs::MParticle &right) 
{
  // CONVERT here ------------------
  ParticleID = right.GetParticleID ();
  HadronicGeneration = right.GetHadronicGeneration ();
  ObservationLevel = right.GetObservationLevel ();
  CorsikaID = ParticleID*1000+(HadronicGeneration%100)*10+ObservationLevel;
  
  Px = right.GetPx ();
  Py = right.GetPy ();
  Pz = right.GetPz ();
  
  x = right.GetX ();
  y = right.GetY ();
  Time = right.GetTime ();

  Weight = right.GetWeight ();
}

TParticle::TParticle(const crs::MEventHeader& header , 
		     const crs::GroundParticle& p) 
{
  HadronicGeneration = p.GetHadronicGeneration();
  ObservationLevel = 0; //p.GetObservationLevel();
  CorsikaID = p.GetCorsikaId();
  ParticleID = p.Id;

  const double ARRANR = header.GetArrayRotation();
  const double COSANG = cos(ARRANR);
  const double SINANG = sin(ARRANR);   	     
  
  Px = p.Px * COSANG + p.Py * SINANG;
  Py = p.Py * COSANG - p.Px * SINANG;
  Pz = p.Pz;
  x =  p.x * COSANG + p.y * SINANG;
  y =  p.y * COSANG - p.x * SINANG;

  Time = p.time;
  Weight = p.weight;
  
  /*  if (header.GetFlagExtraMuonInformation()) {
      cout << "TParticle with additionnal muon information not yet implemented !! " << endl;
      }*/
}

TParticle::TParticle () :    
  CorsikaID (0),
  ParticleID (0),
  ObservationLevel (0),
  HadronicGeneration (0),
  Px (0),
  Py (0),
  Pz (0),
  x (0),
  y (0),
  Time (0),
  Weight (0) {
}

bool
TParticle::IsMuonAdditionalInfo()
  const
{
  return (ParticleID == 75 || ParticleID == 76);
}

double 
TParticle::GetMuonProductionHeight() 
  const
{
  if (IsMuonAdditionalInfo()) {
    return Time; // special corsika mapping ...
    // use only if IsMuonAdditioanlInfo() is true
  } 
  cerr << "\nTParticle::GetMuonProductionHeight() called for particle-id=" << ParticleID << ", which is NOT a muon-additional-info! Please check with \"IsMuonAdditioanlInfo()\"...\n" << endl;
  return 0;
}

void
TParticle::Dump()
const
{
  cout << " CorsikaID=" << CorsikaID;
  cout << " ParticleID=" << ParticleID;
  cout << " ObservationLevel=" << ObservationLevel;
  cout << " HadronicGeneration=" << HadronicGeneration;
    
  cout << " Px=" << Px;
  cout << " Py=" << Py;
  cout << " Pz=" << Pz;
  
  cout << " x=" << x;
  cout << " y=" << y;
    
  cout << " Time=" << Time;
  cout << " Weight=" << Weight;
  cout << endl;

}
