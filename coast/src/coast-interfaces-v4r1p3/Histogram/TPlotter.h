/***************************************************************************
 *                                                                         *
 *  Copyright and any other appropriate legal protection of these          *
 *  computer programs and associated documentation are reserved in         *
 *  all countries of the world.                                            *
 *                                                                         *
 *  These programs or documentation may not be reproduced by any method    *
 *  without prior written consent of Karlsruhe Institute of Technology     *
 *  delegate. Commercial and military use are explicitly forbidden.        *
 *                                                                         *
 *  Karlsruhe Institute of Technology welcomes comments concerning the     *
 *  code but undertakes no obligation for maintenance of the programs,     *
 *  nor responsibility for their correctness, and accepts no liability     *
 *  whatsoever resulting from the use of its programs.                     *
 *                                                                         *
 ***************************************************************************/

#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <list>
#include <vector>
#include <map>

#define HISTOGRAMING_VERSION "v4r0"

class TH1;
class TFile;
class TObject;

class TCorsika;

namespace crs {
  class CParticle;
  class CInteraction;
  class MEventHeader;
  class MEventEnd;
  class MRunHeader;
};




class TPlotter {

  // To keeps the histrogramms for one (depth) layer for a specific type of particles
  typedef std::map<std::string, TH1*> HistDef;
  typedef HistDef::iterator HistDefIterator;
  typedef HistDef::const_iterator HistDefConstIterator;
  
public:
  TPlotter();
  ~TPlotter();

  void Welcome() const;
  
  void SetDirectory(std::string dir) {fDirectory=dir;}
  void SetFileName(std::string file) {fFileName=file;}
  void SetThinning(bool thinn) {fThinning=thinn;}
  void SetSlant(bool slant) {fSlant=slant;}
  void SetCurved(bool flag) {fCurved=flag;}
  void SetStackInput(bool flag) {fStackInput=flag;}
  void SetPreshower(bool flag) {fPreshower=flag;}
    
  void SetRunHeader(const crs::MRunHeader& header);
  void SetShowerHeader(const crs::MEventHeader& header);
  void SetShowerTrailer(const crs::MEventEnd& trailer);
    
  void SetEvent(int no) {fEventNo=no;}
  void SetRun(int no) {fRunNo=no;}
  void SetShowerZenith(float zenith);
  void SetShowerAzimuth(float azimuth);

  void SetFirstInteraction(double x, double y, double z, double X, double t);
  void SetShowerAxis(double zenith, double azimuth);
  
  bool IsThinned() const {return fThinning;}
    
  void Init(int nBins, double maxRange=-1);
  void Write();

  void AddTrack(const crs::CParticle& pre, const crs::CParticle& post);
  void AddInteraction(const crs::CInteraction& interaction);

private:
  void Clear();
  void CalculatePlanes();
  
  void Rotate(double x, double y, double z,
	      double& sx, double& sy, double& sz,
	      int inverse);  
  
  // --------- USER FUNCTIONS ----------------
  void InitParticles();
  void InitHistograms(HistDef& hists);
  void FillHistograms(HistDef& hists,
		      const double& x, const double& y, const double& time, const double& r, const double& rm,
		      const double& energy, const double& theta, const double& phi, const double& weight);
  // --------------------------------------
  
  
private:    
  
  int fVerbosityLevel;

  std::string fDirectory;
  std::string fFileName;
  TFile* fFile;
  
  // general information from CORSIKA
  bool fThinning;
  bool fSlant;
  bool fCurved;
  bool fStackInput;
  bool fPreshower;

  // to figure out the primary track
  bool fPrimaryTrack;
  bool fFirstInteraction;

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

  // event information from CORSIKA
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
  
  // for the histogramming
  double fTimeUnderflow;
  double fTimeErrorMargin;
  int fCountInvalid;
  int fCountErrors;
  int fMaxInvalids;
  
  // -------------- particle definitions ---------------
  struct ParticleDef {
    ParticleDef() : name("undefined"), color(0) {}
    ParticleDef(const std::string& theName, int theC)
      : name(theName), color(theC) {}
    std::string name;
    int color;
  };
  std::map<int, ParticleDef> fParticles;
  
  // -------------- Histogram definitions ---------------
  typedef std::map<int, HistDef> HistLayer; // this is just for typing convinience
  
  std::vector<HistLayer> fHists;
  std::vector<double> fLayerDepth;
  std::vector<double> fLayerHeight;
  std::vector<double> fZPlane;
  int fNhist;
  
  // to keep track of versions
  std::string fCorsikaVersion;
  std::string fCoastVersion;
  std::string fHistVersion;
  
  TCorsika* fCORSIKA;
};



#endif
