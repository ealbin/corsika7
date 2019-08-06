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

#ifndef __include_TCorsika_h__
#define __include_TCorsika_h__


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Access to CORSIKA fortran COMMON blocks

   
      INTEGER          LNGMAX
      PARAMETER        (LNGMAX = 1825)
#if __SLANT__
      COMMON /CRSLANT/ RLONG,THCKRL,CTH,STHCPH,STHSPH,RLOFF
      DOUBLE PRECISION RLONG(0:LNGMAX),THCKRL(0:LNGMAX),
     *                 CTH,STHCPH,STHSPH,RLOFF
#endif


      INTEGER          LNGMAX
      PARAMETER        (LNGMAX = 1825)
      COMMON /CRLONGI/ ADLONG,AELONG,APLONG,DLONG,ELONG,HLONG,PLONG,
     *                 SDLONG,SELONG,SPLONG,THSTEP,THSTPI,
     *                 LHEIGH,NSTEP,
     *                 LLONGI,FLGFIT
#if __ANAHIST__
      DOUBLE PRECISION ADLONG(0:LNGMAX,19),AELONG(0:LNGMAX,10),
     *                 APLONG(0:LNGMAX,20),DLONG(0:LNGMAX,19),
     *                 ELONG(0:LNGMAX,10),
     *                 HLONG(0:LNGMAX),PLONG(0:LNGMAX,20),
     *                 SDLONG(0:LNGMAX,19),SELONG(0:LNGMAX,10),
     *                 SPLONG(0:LNGMAX,20),THSTEP,THSTPI
#else
      DOUBLE PRECISION ADLONG(0:LNGMAX,19),AELONG(0:LNGMAX,10),
     *                 APLONG(0:LNGMAX,10),DLONG(0:LNGMAX,19),
     *                 ELONG(0:LNGMAX,10),
     *                 HLONG(0:LNGMAX),PLONG(0:LNGMAX,10),
     *                 SDLONG(0:LNGMAX,19),SELONG(0:LNGMAX,10),
     *                 SPLONG(0:LNGMAX,10),THSTEP,THSTPI
#endif
      INTEGER          LHEIGH,NSTEP
      LOGICAL          LLONGI,FLGFIT

*/

//static const int LNGMAX = 1825+1; // HARDCODED VALUE FROM CORSIKA
static const int LNGMAX = 15000+1;    // HARDCODED VALUE FROM CORSIKA __AFTER__ 2007-10-19 18:40:2 (r287) !!!
extern "C" {
    extern struct CRSLANT {
	double RLONG[LNGMAX];
	double THCKRL[LNGMAX];
	double CTH;
	double STHCPH;
	double STHSPH;
	double RLOFF;
    } crslant_;
}
extern "C" {

    extern struct CRLONGI {
	double ADLONG[19][LNGMAX];
	double AELONG[10][LNGMAX];
	double APLONG[10][LNGMAX];
	double DLONG[19][LNGMAX];
	double ELONG[10][LNGMAX];
	double HLONG[LNGMAX];
	double PLONG[10][LNGMAX];
	double SDLONG[19][LNGMAX];
	double SELONG[10][LNGMAX];
	double SPLONG[10][LNGMAX];
	double THSTEP;
	double THSTPI;
	int    LHEIGH;
	int    NSTEP;
	bool   LLONGI;
	bool   FLGFIT;
    } crlongi_; 

  /*
    extern struct CRLONGI {
	double ADLONG[19][LNGMAX];
	double AELONG[10][LNGMAX];
	double APLONG[11][LNGMAX];
	double DLONG[19][LNGMAX];
	double ELONG[10][LNGMAX];
	double HLONG[LNGMAX];
	double PLONG[10][LNGMAX];
	double AMLONG[LNGMAX];
	double MLONG[LNGMAX];
	double SMLONG[LNGMAX];
	double SDLONG[19][LNGMAX];
	double SELONG[10][LNGMAX];
	double SPLONG[10][LNGMAX];
	double THSTEP;
	double THSTPI;
	int    LHEIGH;
	int    NSTEP;
	bool   LLONGI;
	bool   FLGFIT;
    } crlongi_; 
  */
}

// CORSIKA routine:
// DOUBLE PRECISION FUNCTION THICK( ARG )
extern "C" double thick_ (const double &height/* cm */);
extern "C" double heigh_ (const double &depth /* g/cm2 */);
extern "C" double thcksi_ (const double &slantd);
extern "C" double heightd_ (const double &slantd, const double &rImpact);
extern "C" double distand_ (const double &dist);
extern "C" double rhof_ (const double &height/* cm */);


/*
  class TCorsika provides access to internal CORSIKA 
  data structures (COMMONS)
  and functions
*/

class TCorsika {

 public:
  // constructor for flat+curved geometries
  TCorsika(const double zenith,
	   const double obsLevel,
	   const bool slant,
	   const bool curved);

  // constructor for skimming geometries
  TCorsika(const double skimmingHeight,
	   const double obsLevel);
  
  void SetHeightOfFirstInteraction(const double z) {fHeightFirstInt=z;}
  void SetVerbosityLevel(const int i) {fVerbosityLevel=i;}


  //double GetTopOfAtmosphere() const {return cTopOfAtmosphere;}
  
  double GetSlantDepthOfDistance(double dist) const;
  double GetDistanceOfSlantDepth(double slantDepth) const;
  double GetHeightOfSlantDepth(double slantdepth) const;
  double GetVerticalHeightOfDistance(double dist) const;
  double GetVerticalDepthOfHeight(double height) const;
  double GetHeightOfPoint(const double x, const double y, const double z) const;
  double GetDistanceOfHeight(const double height) const;
  double GetDistanceOfPoint(const double x, const double y, const double z) const;
  
  void DumpAtmosphereTable() const;
  
  double GetLastLevelSlantDepth() const;
  
 private:
  // functions for interpolation of tabulated values in curved showers
  long GetBin(double Value, double* Table) const;
  long GetBinReverse(double Value, double* Table) const;
  double GetTabulatedSlantDepthOfDistance(double distance) const;
  double GetTabulatedDistanceOfSlantDepth(double slantdepth) const;
  double GetTabulatedHeightOfDistance(double distance) const;
  double GetTabulatedHeightOfSlantDepth (double slantdepth) const;
  double GetTabulatedSlantDepthOfHeight(double height) const;
  double GetTabulatedDistanceOfHeight(double height) const;
  
 private:
  int fVerbosityLevel;
  
  double fCosZenith;
  double fSinZenith;
  double fZenith;
  double fSkimmingAltitude;
  //double fHeightFirstInteraction;
  double fObservationLevel;
  double fHeightFirstInt;
  bool fSlant;
  bool fCurved;
  bool fSkimming;
  
  double* heightTable;	// pointer to the CORSIKA-provided height table
  double* depthTable;		// pointer to the CORSIKA-provided slant depth table
  double* distanceTable;	// pointer to the CORSIKA-provided slant distance table
  long TableLength;		// number of entries in the tables
  
  bool fDecreasingHeightTable;
};




#endif
