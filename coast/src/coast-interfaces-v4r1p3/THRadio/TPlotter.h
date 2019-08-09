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
#include <list>
#include <vector>
#include <map>

#define HISTOGRAMING_VERSION "v4r1"

class TH3D;
class TH2D;
class TH1D;
class TH1;
class TFile;

class TCorsika;

namespace crs {
  class CParticle;
  class CInteraction;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
};



struct Center {
  Center():x(0),y(0),w(0),x2(0),y2(0),w2(0) {}
  double x;
  double y;
  double w;
  double x2;
  double y2;
  double w2;
};

// This keeps the histrogramms for one (depth) layer for a specific 
// type of particles
class Hists {
public:
  Hists() : Time(0), Angle(0), Energy(0) {}
  Hists(TH3D *t, TH3D *a, TH1D* e) : Time(t), Angle(a), Energy(e) {}
  //~Hists () {delete Time; delete Angle;}
  //   void Clear () {delete Time; delete Angle;}
  
  std::map<int,Center> center; // vector sum above all particles in plane

  TH3D* Time;
  TH3D* Angle;
  TH1D* Energy;
  //TH2D* Energy;
  //TH2D* Lateral;
};

// this is just for typing convinience
typedef std::map<int, Hists>  HistLayer;




class TPlotter {

public:
  TPlotter ();
  ~TPlotter ();

  void Welcome() const;
  void SetDirectory (std::string dir) {fDirectory=dir;}
  void SetFileName (std::string file) {fFileName=file;}
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
  
  //void SetZfirst (float start) {fZstart=start;}
  //void SetZlast (float stop) {fZstop=stop;}

  bool IsThinned() const {return fThinning;}
    
  void Init (const int nBins, const double maxRange=2500); // [g/cm2]
  void Write ();

  void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);
  void AddInteraction (const crs::CInteraction& interaction);

private:
  void Clear();
  void CalculatePlanes();
  
  void Rotate(double x, double y, double z,
	      double &sx, double &sy, double &sz,
	      int inverse) const;
  
  TH1* Interpolate3D(TH1* h1, double age1,
		     TH1* h2, double age2,
		     double age) const;
  
  
  
private:    
  
  int fVerbosityLevel;
  int fMaxInvalids;
  double fLateralScaleRadius;

  std::string fDirectory;
  std::string fFileName;
  TFile *fFile;
  
  bool fThinning;
  bool fSlant;
  bool fCurved;
  bool fStackInput;
  bool fPreshower;
  double fTimeErrorMargin;
  
  // to figure out the primary track
  bool  fPrimaryTrack;
  bool  fFirstInteraction;
  
  // point of the first interaction
  double fFirstInteractionX;
  double fFirstInteractionY;
  double fFirstInteractionZ;
  double fFirstInteractionTime;
  double fFirstInteractionDist;
  double fFirstInteractionSlantDepth; // slant depth

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
  double fHeightFirstInt; 
  float fXmax;
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
  
  // for the histograms
  double fTimeUnderflow;
  double fRadiusUnderflow;
  double fEnergyUnderflow;
  int fCountInvalid;
  double fdepthStart;
  double fdeltaSlice;

//double itssavedemweight;
//double itssavedhdweight;
  
  struct ParticleDef {
      std::string name;
      int color;
  };
  std::map<int, ParticleDef> fParticles;
  
  std::vector<HistLayer> fHists;
  std::vector<double> fLayerDepth;
  std::vector<double> fLayerHeight;
  std::vector<double> fLayerDist;
  int fNhist;
  
  // to keep track of versions
  std::string fCorsikaVersion;
  std::string fCoastVersion;
  std::string fHistVersion;

  TCorsika* fCORSIKA;

};



#endif
