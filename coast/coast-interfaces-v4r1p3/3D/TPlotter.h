#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <vector>
#include <map>

class TFile;
class TCanvas;
class TObjArray;

class TCorsika;

namespace crs {
  class CParticle;
  class CInteraction;
  class MEventHeader;
  class MRunHeader;
  class MEventEnd;
};



class TPlotter {
  
  struct ParticleType {
    ParticleType() : lines(0) {}
    std::string name;
    int color;
    int maxDraw;
    double minEnergy;
    int count;
    TObjArray* lines;
  };

public:
    
public:
  TPlotter ();
  ~TPlotter ();

  void SetFileName (std::string file) {fFileName=file;}
  void Set(const bool thinn, 
	   const bool curved, 
	   const bool slant, 
	   const bool stackinput, 
	   const bool preshower) {
    fThinning = thinn;
    fCurved = curved;
    fSlant = slant;
    fStackInput = stackinput;
    fPreshower = preshower;
  }
    
  void SetRunHeader (const crs::MRunHeader &header);
  void SetShowerHeader (const crs::MEventHeader &header);
  void SetShowerTrailer (const crs::MEventEnd &trailer);
    
  void SetEvent (int no) {fEventNo=no;}
  void SetShowerZenith (float zenith);
  void SetShowerAzimuth (float azimuth);
  void SetShowerEnergy (float energy) {fEnergy0=energy;}
  void SetObservationLevel (float level) {fObservationLevel=level;}
  void SetHeightFirstInt (float z) {fHeightFirstInt=z;}
  //void SetXmax (float xmax) {fX_max=xmax;}
  void SetPrimary (int id) {fPrimary=id;}

  bool IsThinned () const {return fThinning;}
    
  void Init ();
  void Write ();
  void Close ();
    
  void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);
  void AddInteraction(const crs::CInteraction &inter);

private:

  bool ParseConfig(const std::string& config);
  void Clear();

  double MoliereRadius(double Density, double Temperature);

  double Temperature(double height);

  void Rotate(const double x, const double y, const double z,
	      double &sx, double &sy, double &sz,
	      const int inverse);

  void SetFirstInteraction(const double, const double, const double);			   


private:    

  TCorsika *fCORSIKA;

  bool fDebug;
  bool fDrawInShowerCS;

  bool fPrimaryTrack;
  bool fDoNotDrawPrimaryTrack;

  bool fFirstInteraction;
  double fFirstInteractionX;
  double fFirstInteractionY;
  double fFirstInteractionZ;
  double fFirstInteractionDist;

  std::string fFileName;
  bool fThinning;
  bool fCurved;
  bool fSlant;
  bool fStackInput;
  bool fPreshower;
  
  int fEventNo;
  double fObservationLevel;
  double fHeightFirstInt;
  float fZenith;
  float fAzimuth;
  float fEnergy0;
  int fPrimary;

  double fCosZenith;
  double fSinZenith;
  double fCosAzimuth;
  double fSinAzimuth;
    
  TFile *fFile;
  TCanvas *fCanvas;
  double fXmin;
  double fXmax;
  double fYmin;
  double fYmax;
  double fZmin;
  double fZmax;
  bool fFirst;
  
  std::map<int, ParticleType> fParticleType;
  /*
    int fNmaxElec;
    int fNmaxMuon;
    int fNmaxHadr;
    
    double fEminElec;
    double fEminMuon;
    double fEminHadr;
    
    TObjArray *fLinesElec;
    TObjArray *fLinesMuon;
    TObjArray *fLinesHadr;
    
    int fNElec;
    int fNMuon;
    int fNHadr;
  */
  /*
  static const int fgColorHadron;
  static const int fgColorElectron;
  static const int fgColorMuon;
  static const int fgColcxHadron;
  static const int fgColcxElectron;
  static const int fgColcxMuon;
  */
};



#endif
