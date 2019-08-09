#ifndef _INCLUDE_CORSIKA_MEVENTHEADER_
#define _INCLUDE_CORSIKA_MEVENTHEADER_

#include <crs/TSubBlock.h>
#include <crs/CorsikaTypes.h>

namespace crs {

  /** 
      \class MEventHeader
      \brief CORSIKA event header (EVTH) sub-block

      This class knows how to access the data stored in this sub-block type.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: MEventHeader.h,v 1.2 2007-10-24 10:13:31 pierog Exp $
  */

  class MEventHeader : public TSubBlock {

  public:
    MEventHeader () {}
    MEventHeader (const TSubBlock &right);
    virtual ~MEventHeader () {}

  public:
    CINT GetEventNumber () const {return (int)fSubBlockData [1];}
    CREAL GetParticleId () const {return fSubBlockData [2];}
    CREAL GetEnergy () const {return fSubBlockData [3];}
    CREAL GetStartingAltitude () const {return fSubBlockData [4];}///< gcm^-2
    CREAL GetFirstTarget () const {return fSubBlockData [5];}
    CREAL GetZFirst () const {return fSubBlockData [6];}          ///< cm
    CREAL GetPx () const {return fSubBlockData [7];}    ///< GeV
    CREAL GetPy () const {return fSubBlockData [8];}    ///< GeV
    CREAL GetPz () const {return fSubBlockData [9];}    ///< GeV
    CREAL GetTheta () const {return fSubBlockData [10];}///< zenith in rad
    CREAL GetPhi () const {return fSubBlockData [11];}  ///< azimuth in rad
	
    CINT GetNRandomSequences () const {return (int)fSubBlockData [12];}
	
    CINT GetSeed (int index) const 
    {return (index<fSubBlockData[12] ? (int)fSubBlockData[13+index*3]:-1);}
    CINT GetInitialCallsMod (int index) const 
    {return (index<fSubBlockData[12] ? (int)fSubBlockData[14+index*3]:-1);}
    CINT GetInitialCallsDiv (int index) const 
    {return (index<fSubBlockData[12] ? (int)fSubBlockData[15+index*3]:-1);}
	
    CREAL GetRunNumber () const {return fSubBlockData [43];}
    CINT GetDateStart () const {return (int)fSubBlockData [44];} ///< yymmdd
    CREAL GetVersion () const {return fSubBlockData [45];}
	
    CINT GetNObservationLevels () const {return (int)fSubBlockData [46];}
    CREAL GetObservationHeight (int index) const 
    {return (index<fgMaxObsLevels ? fSubBlockData [47+index]:-1);}//cm
	
    CREAL GetSpectralSlope () const {return fSubBlockData [57];}
    CREAL GetEMin () const {return fSubBlockData [58];} ///< in GeV
    CREAL GetEMax () const {return fSubBlockData [59];} ///< in GeV
	
    CREAL GetCutoffHadrons () const {return fSubBlockData [60];}///< in GeV
    CREAL GetCutoffMuons () const {return fSubBlockData [61];}  ///< in GeV
    CREAL GetCutoffElectrons () const {return fSubBlockData [62];}///< in GeV
    CREAL GetCutoffPhotons () const {return fSubBlockData [63];}///< in GeV
	
    CREAL GetNFLAIN () const {return fSubBlockData [64];}
    CREAL GetNFLDIF () const {return fSubBlockData [65];}
    CREAL GetNFLPI0 () const {return fSubBlockData [66];}
    CREAL GetNFLPIF () const {return fSubBlockData [67];}
    CREAL GetNFLCHE () const {return fSubBlockData [68];}
    CREAL GetNFRAGM () const {return fSubBlockData [69];}
	
    CREAL GetBx () const {return fSubBlockData [70];}///< mag field in mu T
    CREAL GetBz () const {return fSubBlockData [71];}///< mag field in mu T
	
    CREAL GetFlagEGS4 () const {return fSubBlockData [72];}
    CREAL GetFlagNKG () const {return fSubBlockData [73];}
	
    CREAL GetHadronicLowEModell () const {return fSubBlockData [74];}
    CREAL GetHadronicHighEModell () const {return fSubBlockData [75];}

    CREAL GetFlagCherenkov () const {return fSubBlockData [76];}
    CREAL GetFlagNeutrino () const {return fSubBlockData [77];}
    CREAL GetFlagCurved () const {return fSubBlockData [78];}

    ///1: IBM, 2: Transputer, 3: DEC/UNIX, 4: Mac, 5: VAX/VMS, 6: GNU/Linux
    CINT GetFlagComputer () const {return (int)fSubBlockData [79];} 
	
    CREAL GetThetaMin () const {return fSubBlockData [80];}   ///< degrees
    CREAL GetThetaMax () const {return fSubBlockData [81];}   ///< degrees
    CREAL GetPhiMin () const {return fSubBlockData [82];}     ///< degrees
    CREAL GetPhiMax () const {return fSubBlockData [83];}     ///< degrees
	
    CREAL GetCherenkovBunch () const {return fSubBlockData [84];}
    CREAL GetCherenkovNumberX () const {return fSubBlockData [85];}
    CREAL GetCherenkovNumberY () const {return fSubBlockData [86];}
    CREAL GetCherenkovGridX () const {return fSubBlockData [87];}
    CREAL GetCherenkovGridY () const {return fSubBlockData [88];}    ///< cm
    CREAL GetCherenkovDetectorX () const {return fSubBlockData [89];}
    CREAL GetCherenkovDetectorY () const {return fSubBlockData [90];}///< cm
    CREAL GetCherenkovOutputFlag () const {return fSubBlockData [91];} ///< 0: PARTOUT 1:CKOUT
	
    CREAL GetArrayRotation () const {return fSubBlockData [92];}
    CREAL GetFlagExtraMuonInformation () const {return fSubBlockData [93];}
	
    CREAL GetMultipleScatteringStep () const {return fSubBlockData [94];}
    CREAL GetCherenkovBandwidthMin () const {return fSubBlockData [95];}
    CREAL GetCherenkovBandwidthMax () const {return fSubBlockData [96];}//nm
    CINT GetNUsesOfEvent () const {return (int)fSubBlockData [97];}
	
    CREAL GetCherenkovCoreX (int index) const 
    {return (index<20 ? fSubBlockData [98+index]  : -1);} // cm
    CREAL GetCherenkovCoreY (int index) const 
    {return (index<20 ? fSubBlockData [118+index] : -1);} // cm

    CREAL GetFlagSIBYLL () const {return fSubBlockData [138];} // cm
    CREAL GetFlagSIBYLLCross () const {return fSubBlockData [139];}
    CREAL GetFlagQGSJET () const {return fSubBlockData [140];} // cm
    CREAL GetFlagQGSJETCross () const {return fSubBlockData [141];}
    CREAL GetFlagDPMJET () const {return fSubBlockData [142];} // cm
    CREAL GetFlagDPMJETCross () const {return fSubBlockData [143];}
    CREAL GetFlagVENUSCross () const {return fSubBlockData [144];}

    CREAL GetFlagMuonMultiple () const {return fSubBlockData [145];} // 0: Gauss, 1: Moilere
    CREAL GetNKGRadialRange () const {return fSubBlockData [146];} // cm

    CREAL GetEFractionThinningH () const {return fSubBlockData [147];} // Energy fraction of thinning level hadronic
    // These are in the CORSIKA manual but not in Lukas's original code
    CREAL GetEFractionThinningEM () const {return fSubBlockData [148];} // Energy fraction of thinning level EM
    CREAL GetWMaxHadronic () const {return fSubBlockData [149];} // cm
    CREAL GetWMaxEM () const {return fSubBlockData [150];}
    CREAL GetRMaxThinning () const {return fSubBlockData [151];}

    CREAL GetInnerAngle () const {return fSubBlockData [152];} // cm
    CREAL GetOuterAngle () const {return fSubBlockData [153];}

    CREAL GetTransitionEnergy() const {return fSubBlockData [154];}
    CREAL GetSkimmingIncidence() const {return fSubBlockData [155];}
    CREAL GetSkimmingAltitude() const {return fSubBlockData [156];}
    CREAL GetStartingHeight() const {return fSubBlockData [157];}
  };
};

#endif
