#ifndef _INCLUDE_CORSIKAIO_TSHOWER_
#define _INCLUDE_CORSIKAIO_TSHOWER_

#include <TObject.h>

class TTree;

namespace crs {

  class MRunHeader;
  class MRunEnd;
  class MEventHeader;
  class MEventEnd;
  class MParticleBlock;
};

namespace crs2r {
  class TC2R;
};

namespace crsIO {

  /** 
      \class TShower
      \brief All available info for a CORSIKA shower.

      All available info for one simulated CORSIKA event.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TShower.h,v 1.1.1.1 2007-07-31 07:00:25 rulrich Exp $
  */

  class TShower : public TObject {
    
    friend class crs2r::TC2R;

  public:
    TShower ();
    ~TShower ();
    
    TTree *GetParticles (const int /* ObservationLevel */) const 
    {return fParticles;}
	
  protected:
    void Clear ();

    void AddRunHeader (const crs::MRunHeader &SubBlock);
    void AddRunEnd (const crs::MRunEnd &RunEnd);
    void AddEventHeader (const crs::MEventHeader &SubBlock);
    void AddEventEnd (const crs::MEventEnd &SubBlock);
    void AddLongitudinal (const crs::MParticleBlock &SubBlock);
 
  private:
    TTree *fParticles;
    TTree *fCherenkov;
    
  public:
    int EventID;

    float Energy;            // GeV
    float StartingAltitude;  // g/cm2
    int FirstTarget;
    float FirstHeight;        // cm
    float Theta;               // rad
    float Phi;             // rad


    int RandomSeed [10]; // #max = 10
    int RandomOffset [10]; // #max = 10
    
    // statistics
    float nPhotons;        // weighted
    float nElectrons;      // weighted
    float nHadrons;
    float nMuons;
    int nParticlesWritten; // +cherenkov -muonaddinfos
    int nPhotonsWritten;
    int nElectronsWritten;
    int nHadronsWritten;
    int nMuonsWritten;


    // NKG info
    //std::vector <float> LateralDistribution;


    // Gaisser Hillas charged particles
    float GH_Nmax;
    float GH_t0;
    float GH_tmax;
    float GH_a;
    float GH_b;
    float GH_c;
    float GH_Chi2;
	

	    



    int nPreshower;
	
    //int nPhotonsWritten;
    //int nElectronsWritten;
    //int nHadronsWritten;
    //int nMuonsWritten;
	

	
    // ... NKG
	
    // -------
    int CPUtime; // [s]
	
	
    ClassDef (TShower, 2)  
  };

};

#endif
