#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <list>
#include <vector>
#include <fstream>
#include <map>
#include <crs/MEventHeader.h>
//#include <COASTconfig.h>
#include "Plane.h"

#define PLANE_VERSION "0.1 alpha"
/*
class TH3D;
class TH2D;
class TH1D;
class TH1;

namespace crs2r {
    class TC2R;
};
*/
namespace crsRead {
  class TSubBlockIO;
  class MParticleSubBlockOutput;
};

namespace crs {
  class CParticle;
  class GroundParticle;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
  class TSubBlock;
};




class TPlotter {

public:
  enum OutputOption {
    eNONE,
    //    eROOT,
    eCORSIKA
  };

public:
  TPlotter (OutputOption option=eNONE);
  ~TPlotter ();

  void Welcome() const;
  void Init();
  void Close();
  
  void SetOutputOption(OutputOption option) {fOutputOption=option;}
  void SetDirectory (std::string dir) {fDirectory=dir;}
  void SetFileName (std::string file) {fFileName=file;}
  void SetThinning (bool thinn) {fThinning=thinn;}
  void SetSlant (bool slant) {fSlant=slant;}
  void SetCurved (bool flag) {fCurved=flag;}
  void SetStackInput(bool flag) {fStackInput=flag;}
  // void SetCollectingPlane()                       // TODO: Define and Implement 
    
  void SetRunHeader (const crs::MRunHeader &header);
  void SetShowerHeader (const crs::MEventHeader &header);
  void SetShowerTrailer (const crs::MEventEnd &trailer);
    
  void SetEvent (int no) {fEventNo=no;}
  void SetRun (int no) {fRunNo=no;}
  void SetShowerZenith (float zenith);
  void SetShowerAzimuth (float azimuth);

  void SetShowerAxis (double zenith, double azimuth);
  
  bool IsThinned() const {return fThinning;}
  
  void Write ();
  void Write (const crs::TSubBlock &SubBlock);
  void Write (const crs::GroundParticle& p);

  void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);

private:
  void Clear ();
  
  void Rotate (double x, double y, double z,
	       double &sx, double &sy, double &sz,
	       int inverse);
  
  
  
private:    
  
  int fVerbosityLevel;
  std::string fDirectory;
  std::string fFileName;
  
  bool fThinning;
  bool fSlant;
  bool fCurved;
  bool fStackInput;
  
  // catesian shower axis
  double fAxisX;
  double fAxisY;
  double fAxisZ;
  
  // for rotations
  double fCosZenith;
  double fSinZenith;
  double fCosAzimuth;
  double fSinAzimuth;
  
  // information from CORSIKA
  float  fSkimmingAltitude;
  bool   fSkimming;
  int    fEventNo;
  int    fRunNo;
  double fMaxShowerRange;
  float  fZenith;
  float  fAzimuth;
  float  fShowerEnergy;
  int    fShowerPrimary;
  double fXOff;
  double fYOff;
  double fXPlane;
  double fYPlane;
  double fZPlane;
  double fXNorm;
  double fYNorm;
  double fZNorm;
  
  // to keep track of versions
  std::string fCorsikaVersion;
  std::string fCoastVersion;
  std::string fHistVersion;

  // output object options
  OutputOption fOutputOption;
  //  crs2r::TC2R* fC2R;
  
  std::ofstream* fOutputStream;
  crsRead::TSubBlockIO* fSubBlockIO;
  crsRead::MParticleSubBlockOutput* fParticleSubBlockIO;
  static const int fgNBlockSizeInfo;

  
  crs::MEventHeader fEventHeader;

  crs::Plane fPlane;
};



#endif
