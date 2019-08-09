/* $Id: TRun.cc,v 1.1.1.1 2007-07-31 07:00:32 rulrich Exp $   */

#include <crsIO/TRun.h>
using namespace crsIO;

#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>

#include <TTree.h>

// -----------------
ClassImp(TRun);


TRun::TRun() {

  Clear();
}


TRun::~TRun() {

  ObservationLevel.clear();
  //RandomSeed.clear();
}


void TRun::Clear() {
    
  RunID = 0;
  Date = 0;
  Version = 0;
    
  ObservationLevel.clear();

  ParticleID = 0;

  EnergySlope = 0;
  EnergyMin = 0;
  EnergyMax = 0;

  AzimuthMin = 0;
  AzimuthMax = 0;

  ZenithMin = 0;
  ZenithMax = 0;

  CutHadrons = 0;
  CutMuons = 0;
  CutElectrons = 0;
  CutPhotons = 0;

  // atmosphere
  AtmosphereA.clear();
  AtmosphereB.clear();
  AtmosphereC.clear();

  BFieldX = 0;
  BFieldZ = 0;

  // Flags
  EGS4 = false;
  NKG = false;
  Cherenkov = false;
  Neutrino = false;
  Curved = false;
  MuonAdditionalInfo = false;
    
  MuonMultScatteringMoliere = false;

  RadialRangeNKG = 0;

  // Hadronic interactions
  TransitionEnergy = 0; // HILOW
  LowEHadModel = 0;
  HighEHadModel = 0;
  VersionSIBYLL_interaction = 0;  // 1=v1.6, 2=v2.1
  VersionSIBYLL_crosssection = 0;
  VersionQGSJET_interaction = 0;  // 1=old, 2=QGSJET01C
  VersionQGSJET_crosssection = 0;
  VersionDPMJET_interaction = 0;  // 1=DPMJET
  VersionDPMJET_crosssection = 0;
  VersionVENUSNEXUS_crosssection = 0; // 1=venus, 2=nexus

  // HDPM
  NFLAIN = 0;
  NFLDIF = 0;
  NFLPI0 = 0;
  NFLPIF = 0;
  NFLCHE = 0;
  NFRAGM = 0;

  
  // thinning
  HadronicThinningFraction = 0;  // EFRCTHN
  EMThinningFraction = 0;        // EFRACTN * THINRAT
  HadronicThinningtLimit = 0;       // WMAX
  EMThinningLimit = 0;             //WMAX * WEITRAT
  RadialThiningRMax = 0; // cm

  ViewConeMin = 0; // deg
  ViewConeMax = 0; // deg

  // Cherenkov
  nCherenkovDetectorsX = 0;
  nCherenkovDetectorsY = 0;
  GridCherenkovDetectorX = 0; // cm
  GridCherenkovDetectorY = 0; // cm
  LengthCherenkovDetectorX = 0; // cm
  LengthCherenkovDetectorY = 0; // cm
  CherenkovOutputSeparate = 0;
    
  CherenkovBandwidthMin = 0; // nm
  CherenkovBandwidthMax = 0; // nm
    
  nUseCherenkovEvent = 0;
  CherenkovCoreX.clear();
  CherenkovCoreY.clear();
    
  OrientationArray = 0; /// rad between array-xand magn.orth
    
  StepLengthFactorMultiScatter = 0; // EGS4


  Computer = 0;

    
  /*
    int i;
    for(i=0; i<50; i++)
    C [i] = 0;

    for(i=0; i<40; i++)
    CKA [i]=0;

    for(i=0; i<5; i++) {
    CETA [i] = 0;
    AATM [i] = 0;
    BATM [i] = 0;
    CATM [i] = 0;
    }

    for(i=0; i<11; i++)
    CSTRBA [i] = 0;
  */
}

void TRun::AddRunHeader(const crs::MRunHeader &RunHeader) {
    
  RunID =(int) RunHeader.GetRunID();
  Date = RunHeader.GetDateStart();
  Version = RunHeader.GetVersion();
    
  int nObsLevel =(int) RunHeader.GetNObservationLevels();
  ObservationLevel.resize(nObsLevel);
  for(int iObsLevel = 0; 
       iObsLevel<nObsLevel;
       ++iObsLevel) {
	
    ObservationLevel [iObsLevel] =
      RunHeader.GetObservationHeight(iObsLevel);
  }

  EnergySlope = RunHeader.GetSpectralSlope();
  EnergyMin = RunHeader.GetEMin();
  EnergyMax = RunHeader.GetEMax();

  CutHadrons = RunHeader.GetCutoffHadrons();
  CutMuons = RunHeader.GetCutoffMuons();
  CutElectrons = RunHeader.GetCutoffElectrons();
  CutPhotons = RunHeader.GetCutoffPhotons();

  // atmosphere
  AtmosphereA.resize(5);
  AtmosphereB.resize(5);
  AtmosphereC.resize(5);
  int iAtm;
  for(iAtm=0; iAtm<5; iAtm++) {
    AtmosphereA [iAtm] = RunHeader.GetAtmosphereA(iAtm);
    AtmosphereB [iAtm] = RunHeader.GetAtmosphereB(iAtm);
    AtmosphereC [iAtm] = RunHeader.GetAtmosphereC(iAtm);
  }

  // Flags
  EGS4 = RunHeader.GetFlagEGS4();
  NKG = RunHeader.GetFlagNKG();

  // HDPM
  NFLAIN = RunHeader.GetConstNFLAIN();
  NFLDIF = RunHeader.GetConstNFLDIF();
  NFLPI0 = RunHeader.GetConstNFLPI0();
  NFLPIF = RunHeader.GetConstNFLPIF();
  NFLCHE = RunHeader.GetConstNFLCHE();
  NFRAGM = RunHeader.GetConstNFRAGM();

}




void 
TRun::AddEventHeader(const crs::MEventHeader &EventHeader) 
{
  ParticleID =(int) EventHeader.GetParticleId();

  AzimuthMin = EventHeader.GetPhiMin();
  AzimuthMax = EventHeader.GetPhiMax();

  ZenithMin = EventHeader.GetThetaMin();
  ZenithMax = EventHeader.GetThetaMax();

  BFieldX = EventHeader.GetBx();
  BFieldZ = EventHeader.GetBz();


  Cherenkov = EventHeader.GetFlagCherenkov();
  Neutrino = EventHeader.GetFlagNeutrino();
  Curved = EventHeader.GetFlagCurved();
  MuonAdditionalInfo = EventHeader.GetFlagExtraMuonInformation();
    
  MuonMultScatteringMoliere = EventHeader.GetFlagMuonMultiple();
    
  RadialRangeNKG = EventHeader.GetNKGRadialRange();

  // Hadronic interactions
  //TransitionEnergy = EventHeader.Get();  //HILOW
  LowEHadModel =(int)EventHeader.GetHadronicLowEModell();
  HighEHadModel =(int)EventHeader.GetHadronicHighEModell();
  VersionSIBYLL_interaction =(int) EventHeader.GetFlagSIBYLL();// 1=v1.6, 2=v2.1
  VersionSIBYLL_crosssection =(int) EventHeader.GetFlagSIBYLLCross();
  VersionQGSJET_interaction =(int) EventHeader.GetFlagQGSJET();//1=old,2=01C
  VersionQGSJET_crosssection =(int) EventHeader.GetFlagQGSJETCross();
  VersionDPMJET_interaction =(int) EventHeader.GetFlagDPMJET();  // 1=DPMJET
  VersionDPMJET_crosssection =(int) EventHeader.GetFlagDPMJETCross();
  VersionVENUSNEXUS_crosssection =(int) EventHeader.GetFlagVENUSCross(); // 1=venus, 2=nexus


  
  // thinning
  HadronicThinningFraction = EventHeader.GetEFractionThinningH(); // EFRCTHN
  EMThinningFraction = EventHeader.GetEFractionThinningEM();//EFRACTN *THINRAT
  HadronicThinningtLimit =(int) EventHeader.GetWMaxHadronic();       // WMAX
  EMThinningLimit =(int) EventHeader.GetWMaxEM();             //WMAX * WEITRAT
  RadialThiningRMax = EventHeader.GetRMaxThinning(); // cm

  ViewConeMin = EventHeader.GetInnerAngle(); // deg
  ViewConeMax = EventHeader.GetOuterAngle(); // deg

  // Cherenkov
  nCherenkovDetectorsX =(int) EventHeader.GetCherenkovNumberX();
  nCherenkovDetectorsY =(int) EventHeader.GetCherenkovNumberY();
  GridCherenkovDetectorX = EventHeader.GetCherenkovGridX(); // cm
  GridCherenkovDetectorY = EventHeader.GetCherenkovGridY(); // cm
  LengthCherenkovDetectorX = EventHeader.GetCherenkovDetectorX(); // cm
  LengthCherenkovDetectorY = EventHeader.GetCherenkovDetectorY(); // cm
  CherenkovOutputSeparate = EventHeader.GetCherenkovOutputFlag();
    
  CherenkovBandwidthMin = EventHeader.GetCherenkovBandwidthMin(); // nm
  CherenkovBandwidthMax = EventHeader.GetCherenkovBandwidthMax(); // nm

  int nCherEvents =(int) EventHeader.GetNUsesOfEvent();
  nUseCherenkovEvent = nCherEvents;
  CherenkovCoreX.resize(nCherEvents);
  CherenkovCoreY.resize(nCherEvents);
  for(int iCherEvents = 0;
       iCherEvents<nCherEvents;
       iCherEvents++) {
	
    CherenkovCoreX [iCherEvents] = EventHeader.GetCherenkovCoreX(iCherEvents);
    CherenkovCoreY [iCherEvents] = EventHeader.GetCherenkovCoreY(iCherEvents);
  }

  OrientationArray = EventHeader.GetArrayRotation(); /// rad between array-xand magn.orth
    
  StepLengthFactorMultiScatter = EventHeader.GetMultipleScatteringStep(); // EGS4


  Computer =(int) EventHeader.GetFlagComputer();
}


void 
TRun::AddEventEnd(const crs::MEventEnd& /*SubBlock*/) 
{
}

