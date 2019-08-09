#ifndef _INCLUDE_CORSIKA_MRUNHEADER_
#define _INCLUDE_CORSIKA_MRUNHEADER_

#include <crs/TSubBlock.h>
#include <crs/CorsikaTypes.h>

namespace crs {

/** 
    \class MRunHeader
    \brief CORSIKA run header (RUNH) sub-block

    This class knows how to access the data stored in this sub-block type.

    \author Ralf Ulrich
    \date Thu Feb  3 13:04:50 CET 2005
    \version $Id: MRunHeader.h,v 1.2 2007-10-24 10:13:31 pierog Exp $
*/

    class MRunHeader : public TSubBlock {
	
    public:
	MRunHeader () {}
	MRunHeader (const TSubBlock &right);
	virtual ~MRunHeader () {}

    public:
	
	CREAL GetRunID () const {return fSubBlockData [1];}
	CINT GetDateStart () const {return (int)fSubBlockData [2];} ///< yymmdd
	CREAL GetVersion () const {return fSubBlockData [3];}

	CINT GetNObservationLevels () const {return (int)fSubBlockData [4];}
	CREAL GetObservationHeight (int index) const 
	    {return (index<fgMaxObsLevels ? fSubBlockData [5+index] : -1);}//cm

	CREAL GetSpectralSlope () const {return fSubBlockData [15];}
	CREAL GetEMin () const {return fSubBlockData [16];}  // in GeV
	CREAL GetEMax () const {return fSubBlockData [17];}  // in GeV
	
	CREAL GetFlagEGS4 () const {return fSubBlockData [18];}
	CREAL GetFlagNKG () const {return fSubBlockData [19];}
	
	CREAL GetCutoffHadrons () const {return fSubBlockData [20];} // in GeV
	CREAL GetCutoffMuons () const {return fSubBlockData [21];}   // in GeV
	CREAL GetCutoffElectrons () const {return fSubBlockData [22];}// in GeV
	CREAL GetCutoffPhotons () const {return fSubBlockData [23];} // in GeV

	CREAL GetSamplingPlanePointX () const {return fSubBlockData [74];} // in cm
	CREAL GetSamplingPlanePointY () const {return fSubBlockData [75];} // in cm
	CREAL GetSamplingPlanePointZ () const {return fSubBlockData [76];} // in cm
	CREAL GetSamplingPlaneNormalX () const {return fSubBlockData [77];}
	CREAL GetSamplingPlaneNormalY () const {return fSubBlockData [78];}
	CREAL GetSamplingPlaneNormalZ () const {return fSubBlockData [79];}
	
	/*
	  CREAL[50] GetConstC () const {return fSubBlockData + 24;}
	  CREAL[20] GetConstCC () const {return fSubBlockData + 74;} 
	  CREAL[40] GetConstCKA () const {return fSubBlockData + 94;} 
	  CREAL[5] GetConstCETA () const {return fSubBlockData + 134;}
	  CREAL[11] GetConstCSTRBA () const {return fSubBlockData + 139;}
	  CREAL[4] GetConstUNUSED () const {return fSubBlockData + 150;}
	  CREAL[50] GetConstCAN () const {return fSubBlockData + 154;}
	  CREAL[50] GetConstCANN () const {return fSubBlockData + 204;}
	*/

	CREAL GetAtmosphereLayerBoundary (int index) const
	    {return (index<5 ? fSubBlockData [249+index] : -1);}

	CREAL GetAtmosphereA (int index) const 
	    {return (index<5 ? fSubBlockData [254+index] : -1);}
	CREAL GetAtmosphereB (int index) const 
	    {return (index<5 ? fSubBlockData [259+index] : -1);}
	CREAL GetAtmosphereC (int index) const 
	    {return (index<5 ? fSubBlockData [264+index] : -1);}


	CREAL GetConstNFLAIN () const {return fSubBlockData [269];}  ///< redundant
	CREAL GetConstNFLDIF () const {return fSubBlockData [270];}///< redundant
	CREAL GetConstNFLPI0 () const {return (int)fSubBlockData [271]%100;}///< redundant
	CREAL GetConstNFLPIF () const {return (int)fSubBlockData [271]/100;}///< redundant
	CREAL GetConstNFLCHE () const {return (int)fSubBlockData [272]%100;}///< redundant
	CREAL GetConstNFRAGM () const {return (int)fSubBlockData [272]/100;}///< redundant
    };
    
};

#endif
    
