#ifndef _T_SIMPLECORSIKAREADER_H
#define _T_SIMPLECORSIKAREADER_H

#include <crsRead/TCorsikaReader.h>

#include <crs/TSubBlock.h>

#include <string>

namespace crs {
    class TBlock;
    class MEventHeader;
    class MEventEnd;
    class MRunHeader;
    class MRunEnd;
    class MParticleBlock;
    class MLongitudinalBlock;
    class MParticle;
    class MCherenkov;
    class MMuonProductionInfo;
    class TLongitudinal;
};

/** 
 \class TReaderExample
 \brief shows all the possibilites of the TCorsikaReader interface

 Does nothing, but stream in a CORSIKA binary data file and demonstrates
 how the blocks and subblocks are handed through the TCorsikaReader class.

 \author Ralf Ulrich
 \date Thu Jun 16 10:11:23 CEST 2005
 \version $Id: TReaderExample.h,v 1.1.1.1 2007-07-31 07:00:43 rulrich Exp $ 
*/

class TReaderExample : public crsRead::TCorsikaReader {
	
 public:
    TReaderExample (const std::string &filename, int verbose=0);
    ~TReaderExample ();
    

 private:
    // overload the virtual handler methods
    virtual void HandleSubBlock (const crs::TSubBlock &sb);
    virtual void HandleParticleBlock (const crs::MParticleBlock &p);
    virtual void HandleLongitudinalBlock(const crs::MLongitudinalBlock &s);
    virtual void HandleEventStart (const crs::MEventHeader &sb);
    virtual void HandleEventEnd (const crs::MEventEnd &sb);
    virtual void HandleRunStart (const crs::MRunHeader &sb);
    virtual void HandleRunEnd (const crs::MRunEnd &sb);
    virtual void HandleParticle (const crs::MParticle &p);
    virtual void HandleCherenkov (const crs::MCherenkov &p);
    virtual void HandleMuonProductionInfo (const crs::MMuonProductionInfo &p);

	
    // virtual user-event handlers
    virtual void Init ();
    virtual void Exit ();
	

 private:
};


#endif
