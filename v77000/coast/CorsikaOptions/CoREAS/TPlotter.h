/****************************************************************************
 *                                                                          *
 *  Copyright and any other appropriate legal protection of these           *
 *  computer programs and associated documentation are reserved in          *
 *  all countries of the world.                                             *
 *                                                                          *
 *  These programs or documentation may not be reproduced by any method     *
 *  without prior written consent of Karlsruhe Institute of Technology (KIT)*
 *  ot its delegate. Commercial and military use are explicitly forbidden.  *
 *                                                                          *
 *  The Karlsruhe Institute of Technology welcomes comments concerning the  *
 *  COAST code but undertakes no obligation for maintenance of the programs,*
 *  nor responsibility for their correctness, and accepts no liability      *
 *  whatsoever resulting from the use of its programs.                      *
 *                                                                          *
 ****************************************************************************/

#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <map>
#include <vector>

#include "threevector.h"

#include "config.h" // for access to CORSIKA option flags

class TCorsika;
class GroundArea;
class GroundElement;
class ScenarioParams;

namespace crs {
  class CParticle;
  class CInteraction;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
};

class TPlotter {

public:
  TPlotter ();
  ~TPlotter ();

  void Welcome() const;
  void InitializeRadioSimulation();
  void CommunicateNodeResults();
  void SetDirectory (std::string dir) {fDirectory=dir;}
  void SetFileName (std::string file) {fFileName=file;}
  void SetThreadName (std::string name) {fThreadName=name;}
  void SetThinning (bool thinn) {fThinning=thinn;}
  void SetSlant (bool slant) {fSlant=slant;}
  void SetCurved (bool flag) {fCurved=flag;}
  void SetStackInput(bool flag) {fStackInput=flag;}
  void SetPreshower(bool flag) {fPreshower=flag;}
    
  void SetRunHeader (const crs::MRunHeader &header);
  void SetShowerHeader (const crs::MEventHeader &header);
  void SetShowerTrailer (const crs::MEventEnd &trailer);
    
  void SetEvent (int no) {fEventNo=no;}
  void SetRun (int no) {fRunNo=no;}
  void SetShowerZenith (float zenith);
  void SetShowerAzimuth (float azimuth);

  void SetFirstInteraction(const double x, const double y, const double z, 
			   const double X, const double t);
  void SetShowerAxis(const double zenith, const double azimuth);
  
  bool IsThinned() const {return fThinning;}
    
  void Init();
  void Write();

  void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);
  void AddTrackZHS (const crs::CParticle &pre, const crs::CParticle &post);
  void AddInteraction (const crs::CInteraction& interaction);
  void FillRefractivityTablesFromGDAS(int nPoints, const double* height, const double* refractiveIndex);

  bool IsMainThread() const {return ((fThreadName == "") || fThreadName.substr(fThreadName.size()-9) == "000000001");}

private:

  class Exception { };  // exception class

  double GetRefractiveIndexAtPosition(const ThreeVector& p1) const;
  double GetEffectiveRefractiveIndexBetween(GroundElement* const observer, const ThreeVector& p2) const;
  double GetNumericallyIntegratedRefractiveIndexBetween(const ThreeVector& p1, const double trackDistance, const ThreeVector& travelDirection) const;
  
  double GetHeightAtPosition(const ThreeVector& p1) const;
  int GetIntersectionsLineCircle(double r, double b, double c, double& ax, double& ay, double& bx, double& by);
	
  // vectors and functionality for tabulated atmosphere
	std::vector<double> fRefractivityTable;
	std::vector<double> fIntegratedRefractivityTable;
  std::vector<double> fHeightTable; // in principle this table is not needed, but would have to remove its filling from core CORSIKA part
	double InterpolateInPlanarAtmosphere(double h, const std::vector<double> &vec) const;
  double fDistanceIncrement;

  // vectors and functionality for tabulated atmosphere in curved geometry
  double fTabulationThresholdZenithAngle;
  unsigned int fNumTanThetaBins;
  double fTanThetaBinSpacing;
  std::vector<std::vector<float> > fIntegratedRefractivityInCurvedAtmosphereTables;
  std::vector<double> fDistanceTableIncrements;
  std::vector<double> fDistanceTableMinDistances; // all the same at the moment, could remove table
  double InterpolateInCurvedAtmosphere(int tanThetaBin, double distance, const std::vector<std::vector<float> >& table) const;

  //~ void Rotate(double x, double y, double z,
	      //~ double &sx, double &sy, double &sz,
	      //~ int inverse) const;
  
  int fVerbosityLevel;

  std::string fDirectory;
  std::string fFileName;
  std::string fThreadName;
  
  bool fThinning;
  bool fSlant;
  bool fCurved;
  bool fStackInput;
  bool fPreshower;
  
  // to figure out the primary track
  bool  fPrimaryTrack;
  
  // point of the first interaction
  double fFirstInteractionX;
  double fFirstInteractionY;
  double fFirstInteractionZ;
  double fFirstInteractionTime;
  double fFirstInteractionDist;
  double fFirstInteractionSlantDepth; // slant depth

  // intermediate storage of magnetic field data
  double fBx;
  double fBz;

  // catesian shower axis
  //~ double fAxisX;
  //~ double fAxisY;
  //~ double fAxisZ;
  
  // for rotations
  double fCosZenith;
  double fSinZenith;
  double fCosAzimuth;
  double fSinAzimuth;
  
  // information from CORSIKA
  double fHeightFirstInt; 
  float fXmax; // careful, this is not in g/cm^2
  float fSkimmingAltitude;
  bool fSkimming;
  int fEventNo;
  int fRunNo;
  double fObservationLevel;
  double fMaxShowerRange;
  float fZenith;
  float fAzimuth;
  float fShowerEnergy;
  int fShowerPrimary;
  
  struct ParticleDef {
      std::string name;
      int color;
  };
  
  std::map<int, ParticleDef> fParticles;
  
  // to keep track of versions
  std::string fCorsikaVersion;
  std::string fCoastVersion;

  TCorsika* fCORSIKA;

  // for directRadio
  GroundArea* fGroundArea;
  double fCoreHitTime;
  ThreeVector fRadioCore;
  double fSeaLevelRefractivity;
  ScenarioParams* itsScenarioParams;
  ThreeVector itsShowerAxis; // set in InitializeRadioSimulation()

/*** Histogramming
  
  double itsHeightBinsMeters;
  double itsLengthBinsMeters;
  long itsNumHeightBins;
  long itsNumLengthBins;
  std::vector<long> itsHeightHistograms;

*/
 
};



#endif
