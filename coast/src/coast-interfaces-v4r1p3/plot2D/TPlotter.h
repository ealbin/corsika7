#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <list>
#include <vector>
#include <map>

namespace crs {
  class CParticle;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
};

extern "C" double thick_ (const double &height/* cm */);
extern "C" double heigh_ (const double &depth /* g/cm2 */);
extern "C" double thcksi_ (const double &slantd);
extern "C" double heightd_ (const double &slantd, const double &rImpact);
extern "C" double distand_ (const double &dist);

class TPlotter {

  struct ParticleType {
    std::string name;
    int id;
    int red;
    int green;
    int blue;
  };

public:
  TPlotter ();
  ~TPlotter ();
  
  void SetFromGroundParticleFile(const bool from2d) {fFromGround=true;}
  
  void SetFileName (std::string file) {fFileName=file;}
  void SetThinning (bool thinn) {fThinning=thinn;}

  bool IsThinned() const {return fThinning;}
    
  void SetRunHeader (const crs::MRunHeader &header);
  void SetShowerHeader (const crs::MEventHeader &header);
  void SetShowerTrailer (const crs::MEventEnd &trailer);
    
  void SetEvent (int no) {fEventNo=no;}
  void SetRun (int no) {fRunNo=no;}
  void SetShowerZenith (float zenith);
  void SetShowerAzimuth (float azimuth);

  void SetFirstInteraction(double x, double y, double z, double t);
  void SetShowerAxis (double zenith, double azimuth);
  
  void SetObservationLevel (float level) {fObservationLevel=level;}
  void SetHeightFirstInt (float z) {fHeightFirstInt=z;}
  void SetXmax (float xmax);

  bool ParseConfig(const std::string& configName);
  void ReadInit();
  void Init();
  void Write();

  void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);

  void DrawLine(const double red, const double green, const double blue,
                const double w,
                const double x1, const double y1,
                const double x2, const double y2);

private:
  void Clear ();
  
  void Rotate (double x, double y, double z,
	       double &sx, double &sy, double &sz,
	       int inverse);

private:    

  int fVerbosityLevel;

  bool fFromGround;
  
  bool fPlotXmaxLevel;
  bool fPlotObservationLevels;
  
  bool fPrimaryTrack;
  bool fPrimaryTrackFirst;
  float fPrimaryTrackEnergy;
  int fPrimaryTrackId;
  
  double fFirstParticleX;
  double fFirstParticleY;
  double fFirstParticleZ;
  double fFirstParticleTime;

  // catesian
  double fAxisX;
  double fAxisY;
  double fAxisZ;
  
  // for rotations
  double fCosZenith;
  double fSinZenith;
  double fCosAzimuth;
  double fSinAzimuth;
  
  double fHeightFirstInt;
  double fDepthFirstInteraction;
  
  std::string fFileName;

  bool fThinning;
  
  double fTopOfAtmosphere;
  int fEventNo;
  int fRunNo;
  double fObservationLevel;

    
  double fXmax;
  
  double *fRed;
  double *fGreen;
  double *fBlue;
  
  double fMaxRed;
  double fMaxGreen;
  double fMaxBlue;

  bool fAutoRange;
  double fPixelPerKm;
  double fXstretch;
  double fXstart;
  double fXstop;
  double fYstart;
  double fYstop;
  int fNX;
  int fNY;
  int fBgRed;
  int fBgGreen;
  int fBgBlue;
  double fWeightBoostMuons;
  double fWeightBoostHadrons;
  std::string fExtension;

  std::map<int, ParticleType> fParticleType;  
};



#endif
