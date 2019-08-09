#ifndef _INCLUDE_CORSIKAIO_TRUN_
#define _INCLUDE_CORSIKAIO_TRUN_

#include <TObject.h>

#include <vector>

class TTree;

namespace crs {
  class MRunHeader;
  class MEventHeader;
  class MEventEnd;
};

namespace crs2r {
  class TC2R;
};

/**
   \namespace crsIO
   \brief Representation of CORSIKA data in ROOT-objects

   This namespace contains all classes inherited from ROOT TObject
   that are streamable using the flexible ROOT I/O.
     
   \author Ralf Ulrich
   \date Thu Feb  3 13:01:59 CET 2005
*/

namespace crsIO {

  /** 
      \class TRun
      \brief All available info for a CORSIKA run.

      All available info for a CORSIKA run.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TRun.h,v 1.1.1.1 2007-07-31 07:00:32 rulrich Exp $
  */

  class TRun : public TObject {
    
    friend class crs2r::TC2R;

  public:
    TRun();
    ~TRun();

  protected:
    void Clear();

    void AddRunHeader(const crs::MRunHeader &SubBlock);
    void AddEventHeader(const crs::MEventHeader &SubBlock);
    void AddEventEnd(const crs::MEventEnd &SubBlock);


  public:
    int RunID;
    float Date;
    float Version;
    
    std::vector <float> ObservationLevel; // in cm, #max = 10

    // primary particle
    int ParticleID;

    // energy spectrum range
    float EnergySlope;
    float EnergyMin;
    float EnergyMax;
    
    // azimuth range
    float AzimuthMin;
    float AzimuthMax;

    // zenith range
    float ZenithMin;
    float ZenithMax;

    float CutHadrons; // GeV
    float CutMuons; // GeV
    float CutElectrons; // GeV
    float CutPhotons; // GeV
    
    //float C [50];          constants inside corsika
    //float CKA [40];        
    //float CETA [5];
    //float CSTRBA [11];

    // atmosphere
    std::vector <float> AtmosphereA;
    std::vector <float> AtmosphereB;
    std::vector <float> AtmosphereC;

    // Magnetic fields
    float BFieldX; // muT
    float BFieldZ; // muT
    
    // Flags
    bool EGS4;
    bool NKG;    
    bool Cherenkov;
    bool Neutrino;
    bool Curved;
    bool MuonAdditionalInfo; 

    bool MuonMultScatteringMoliere;

    float RadialRangeNKG; // cm

    // Hadronic interactions
    float TransitionEnergy; // HILOW
    int LowEHadModel;
    int HighEHadModel;
    int VersionSIBYLL_interaction;  // 1=v1.6, 2=v2.1
    int VersionSIBYLL_crosssection;
    int VersionQGSJET_interaction;  // 1=old, 2=QGSJET01C
    int VersionQGSJET_crosssection;
    int VersionDPMJET_interaction;  // 1=DPMJET
    int VersionDPMJET_crosssection;
    int VersionVENUSNEXUS_crosssection; // 1=venus, 2=nexus
	
    // HDPM
    float NFLAIN;
    float NFLDIF;
    float NFLPI0;
    float NFLPIF;
    float NFLCHE;
    float NFRAGM;
	
    // thinning
    float HadronicThinningFraction;  // EFRCTHN
    float EMThinningFraction;        // EFRACTN * THINRAT
    int HadronicThinningtLimit;       // WMAX
    int EMThinningLimit;             //WMAX * WEITRAT
    float RadialThiningRMax; // cm

    float ViewConeMin; // deg
    float ViewConeMax; // deg

    // Cherenkov
    int nCherenkovDetectorsX;
    int nCherenkovDetectorsY;
    float GridCherenkovDetectorX; // cm
    float GridCherenkovDetectorY; // cm
    float LengthCherenkovDetectorX; // cm
    float LengthCherenkovDetectorY; // cm
    bool CherenkovOutputSeparate;

    float CherenkovBandwidthMin; // nm
    float CherenkovBandwidthMax; // nm

    int nUseCherenkovEvent;
    std::vector <float> CherenkovCoreX;
    std::vector <float> CherenkovCoreY;

    float OrientationArray; /// rad between array-xand magn.orth

    float StepLengthFactorMultiScatter; // EGS4


    // used computer/architecture 
    int Computer;
    
    // ...

    ClassDef(TRun, 2);
  };

};

#endif
