#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>

#include <crs2r/TC2R.h>

#include <TPlotter.h>

#include <sstream>
#include <iostream>
using namespace std;


///  Global object to handle all ROOT I/O operations and data conversions
TPlotter gPlotter;




void wrida_ (const CREAL *Data) {

    crs::TSubBlock subBlock (Data, gPlotter.IsThinned ());
    
    if (subBlock.GetBlockType ()!=crs::TSubBlock::ePARTDATA)
	cout << "RU: SubBlock: " << subBlock.GetBlockType () << endl;

    switch (subBlock.GetBlockType ()) {
	
	case crs::TSubBlock::eEVTH:
	    gPlotter.SetShowerHeader (subBlock);
	    gPlotter.Init (50, TPlotter::eDepth);
	    break;
	    
	case crs::TSubBlock::eEVTE:
	    gPlotter.SetShowerTrailer (subBlock);
	    gPlotter.Write ();
	    break;
	    
	default:
	    break;
    }
    
}



void inida_ (const char *filename,
	     const bool &thinning,
	     const bool &curved,
	     const bool &slant,
             const bool &stackinput,
	     int str_length) {

    if (str_length==0)
	str_length = 9;

    std::ostringstream FileName;
    for (unsigned int i=0; i<str_length; i++) {

	if (filename [i]==' ')
	    break;

	FileName << filename [i];
    }
    
    gPlotter.SetFileName (FileName.str ());
    gPlotter.SetThinning (thinning);
}



void cloda_ () {
    
}



void track_ (const crs::CParticle *pre, const crs::CParticle *post) {
    
    gPlotter.AddTrack (*pre, *post);
}




