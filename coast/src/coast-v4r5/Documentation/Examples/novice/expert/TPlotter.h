#ifndef __INCLUDE_TPLOTTER_H__
#define __INCLUDE_TPLOTTER_H__

#include <string>
#include <vector>

class TH1D;
class TH2D;
class TH3D;
class TFile;

namespace crs {
    class CParticle;
    class MEventHeader;
    class MEventEnd;
};


// corsika routine:
// DOUBLE PRECISION FUNCTION THICK( ARG )
extern "C" double thick_ (const double &height/*cm*/); 
extern "C" double heigh_ (const double &depth/*g/cm*/);


class TPlotter {

 public:
    
    enum PlotMode {	
	eHeight,
	eDepth/*,
		eTime*/
    };
    

 public:
    TPlotter ();
    ~TPlotter ();

    void SetFileName (std::string file) {fFileName=file;}
    void SetThinning (bool thinn) {fThinning=thinn;}
    
    void SetShowerHeader (const crs::MEventHeader &header);
    void SetShowerTrailer (const crs::MEventEnd &trailer);
    
    void SetEvent (int no) {fEventNo=no;}
    void SetShowerZenith (float zenith);
    void SetShowerAzimuth (float azimuth);
    void SetShowerEnergy (float energy) {fEnergy0=energy;}
    void SetObservationLevel (float level) {fObservationLevel=level;}
    void SetHeightFirstInt (float z) {fHeightFirstInt=z;}
    void SetXmax (float xmax) {fXmax=xmax;}
    void SetPrimary (int id) {fPrimary=id;}

    void SetZfirst (float start) {fZstart=start;}
    void SetZlast (float stop) {fZstop=stop;}
    void SetNHistograms (int i) {fNhist=i;}
    void SetPlotMode (PlotMode m) {fMode=m;}

    bool IsThinned () const {return fThinning;}
    
    void Init (int nBins, PlotMode mode);
    void Write ();
    
    void AddTrack (const crs::CParticle &pre, const crs::CParticle &post);
    

 private:

    void Clear ();

    /*
    double MoliereRadius (double Temperature, 
			  double Pressure);
    */
    double MoliereRadius (double Density, double Temperature);

    double Temperature (double height);


    void Rotate (double x, double y, double z,
		 double &sx, double &sy, double &sz,
		 int inverse);


 private:    
    bool fPrimaryTrack;

    std::string fFileName;
    bool fThinning;
    float fXmax;
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
    std::vector<TH3D*> fHistTGamma;
    std::vector<TH3D*> fHistTElectron;
    //std::vector<TH3D*> fHistTPositron;
    std::vector<TH3D*> fHistTMuon;
    std::vector<TH3D*> fHistTHadron;

    std::vector<TH2D*> fHistEGamma;
    std::vector<TH2D*> fHistEElectron;
    //std::vector<TH2D*> fHistEPositron;
    std::vector<TH2D*> fHistEMuon;
    std::vector<TH2D*> fHistEHadron;

    std::vector<TH2D*> fHistLGamma;
    std::vector<TH2D*> fHistLElectron;
    //std::vector<TH2D*> fHistLPositron;
    std::vector<TH2D*> fHistLMuon;
    std::vector<TH2D*> fHistLHadron;

    std::vector<TH3D*> fHistAGamma;
    std::vector<TH3D*> fHistAElectron;
    //std::vector<TH3D*> fHistAPositron;
    std::vector<TH3D*> fHistAMuon;
    std::vector<TH3D*> fHistAHadron;

    PlotMode fMode;
    std::vector<double> fZPlane;
    double fZstart;
    double fZstop;
    int fNhist;
};



#endif
