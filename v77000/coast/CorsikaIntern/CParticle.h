/* $Id: CParticle.h 5900 2017-03-20 14:40:18Z rulrich $   */


#ifndef _INCLUDE_CPARTICLE_H__
#define _INCLUDE_CPARTICLE_H__

#include <crs/CorsikaTypes.h>
// #include <crs/MParticle.h> does not work so good

/*
  This is the CORSIKA internal particle definition:

  C  curpar(0)   = phix,    DIRECTION COSINE X-DIRECTION
  C  CURPAR(1)   = PARTICLE TYPE
  C  CURPAR(2)   = GAMMA,  LORENTZ FACTOR IN LAB
  C  CURPAR(3)   = COSTHE, COSINE THETA
  C  CURPAR(4)   = PHI,    PHI
  CCCC  curpar(4)   = phiy,   DIRECTION COSINE Y-DIRECTION     ????????????????
  C  CURPAR(5)   = H,      HEIGHT (TRUE HEIGHT)
  C  CURPAR(6)   = T,      ACCUMULATED TIME IN SEC
  C  CURPAR(7)   = X,      X-POSITION
  C  CURPAR(8)   = Y,      Y-POSITION
  C  CURPAR(9)   = CHI,    PENETRATED MATERIAL UNTIL DECAY OR REACTION
  C                (G/CM**2)  CALCULATED IN BOX2
  C  CURPAR(10)  = BETA,   V/C, CALCULATED IN BOX2
  C  CURPAR(11)  = GCM,    GAMMA  IN CM, CALCULATED IN NUCINT
  C  CURPAR(12)  = ECM,    ENERGY IN CM, CALCULATED IN NUCINT
  +IF, THIN.
  C  CURPAR(13)  = WEIGHT, WEIGHT FOR THINNING
  +ENDIF.
  +IF, CURVED.
  C  CURPAR(14)  = HAPP    APPARENT HEIGHT  IN CARTESIAN COORDINATE SYSTEM
  C  CURPAR(15)  = COSTAP  APPARENT ZENITH ANGLE IN CART.COORDINATE SYSTEM
  C  CURPAR(16)  = COSTEA  ANGLE PARTICLE TO MID DETECT AT CENTER EARTH
  +ENDIF.
  +IF, INTTEST.
  C  CURPAR(17)  = TRANSVERSE MOMENTUM


  In fact the array always has a length of 16 only in the INTTEST it is 
  increased to 17. But the data is only partly filled according to the CORSIKA
  program settings.
  
*/

namespace crs {

  class MEventHeader;
  
  // ***************************************************
  // all the data a CORSIKA ground particle needs
  class GroundParticle {
  public:
    GroundParticle(int id=0, int gen=0, 
		   double theTime=0, double w=0,
		   double thepx=0, double thepy=0, double thepz=0, 
		   double thex=0, double they=0, double thez=0) 
      : Id(id), generation(gen), time(theTime), weight(w),
	Px(thepx), Py(thepy), Pz(thepz), 
	x(thex), y(they), z(thez) {}
    int GetCorsikaId() const { return Id * 1000 + (generation%100) * 10 + 0;}
    int GetHadronicGeneration() const { return generation; }
    int Id;
    int generation;
    double time;
    double weight;
    double Px;
    double Py;
    double Pz;
    double x;
    double y;
    double z;    
  };
  
  // ***************************************************
  struct CParticle {	

  public:
    enum EParticleId {
      eGamma = 1,
      ePositron = 2,
      eElectron = 3,
      eMuonBar = 5,
      eMuon = 6,
      ePi0 = 7,
      ePiP = 8,
      ePiM = 9,
      eKlong = 10,
      eKP = 11,
      eKM = 12,
      eNeutron = 13,
      eProton = 14,
      eProtonBar = 15,
      eKshort = 16,
      eEta = 17,
      eNeutronBar = 25,
    };

  public:
  CParticle() 
    : x(0),
      y(0),
      z(0),
      depth(0),
      time(0),
      energy(0),
      weight(0),
      particleId(0),
      hadronicGeneration(0) {}
    
    void Dump () const;
    // MParticle TransformParticle(const crs::MEventHeader& header) const; does not work so well
    
  public:
	   
    CDOUBLEPRECISION x;
    CDOUBLEPRECISION y;
    CDOUBLEPRECISION z;
    CDOUBLEPRECISION depth;
    CDOUBLEPRECISION time;
    CDOUBLEPRECISION energy;
    CDOUBLEPRECISION weight;
    CINT particleId;
    CINT hadronicGeneration;
    
    
    int GetParticleId() const { return particleId; }
    int GetHadronicGeneration() const { return hadronicGeneration; }
    int GetObservationLevel() const { return 1; }
    int GetCorsikaId() const { return GetParticleId() * 1000 + (GetHadronicGeneration()%100) * 10 + GetObservationLevel();}
    double GetTime() const { return time * 1.e9; } // convert s to ns
    double GetWeight() const { return weight; } 
    // GroundCoordinates TransformToGroundCoordinates(const crs::MEventHeader& header) const;
    

    /*
      CDOUBLEPRECISION cosPhiX;
      CDOUBLEPRECISION particleID;
      CDOUBLEPRECISION gamma;
      CDOUBLEPRECISION cosTheta;
      CDOUBLEPRECISION phi;
      //CDOUBLEPRECISION cosPhiY;        ~~~~~~~~~
      CDOUBLEPRECISION height;
      CDOUBLEPRECISION time;
      CDOUBLEPRECISION x;
      CDOUBLEPRECISION y;
      CDOUBLEPRECISION chi;
      CDOUBLEPRECISION beta;
      CDOUBLEPRECISION gammaCM;
      CDOUBLEPRECISION energyCM;
      CDOUBLEPRECISION weight;
      CDOUBLEPRECISION heightApp;
      CDOUBLEPRECISION cosThetaApp;
      CDOUBLEPRECISION cosToCenterEarth;
      // CDOUBLEPRECISION p_t;   // for INTTEST
      */
  };


  /**
     Definition of helper class to get CParticle data via function
     call to CORSIKA/fortran

     This class does not "own" data. It is just a pointer to some
     common-block data!
   */
  typedef struct __CParticleFortranPtr {
    __CParticleFortranPtr() : fInternal(0) {}
    CParticle* operator->() { return fInternal; }
    operator CParticle*() { return fInternal; }
    operator CParticle**() { return &fInternal; }
  private:
    CParticle* fInternal;
  } CParticleFortranPtr;
  
}; // end namespace crs


#endif
