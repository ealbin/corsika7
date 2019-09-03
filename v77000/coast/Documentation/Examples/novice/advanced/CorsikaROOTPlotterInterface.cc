#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TBlock.h>
#include <crs/TSubBlock.h>

#include <crs2r/TC2R.h>

#include <TCorsikaPlotter.h>

#include <sstream>
#include <iostream>



///  Global object to handle all ROOT I/O operations and data conversions
TCorsikaPlotter gPlotter;
bool gThinned;

/*
  Data is one CORSIKA data-subblock. A CORIKA block consists of 21 SubBlocks.
  A SubBlock can be:
   - thinned mode:     39 (Particles) * 8 (ENTRIES) * 4 (BYTES) 
   - not-thinned mode: 39 (Particles) * 7 (ENTRIES) * 4 (BYTES) 
 */
void wrida_ (const CREAL *Data) {
    
    crs::TSubBlock subBlock (Data, gThinned);
    gPlotter.PlotSubBlock (subBlock);
}



void inida_ (const char *filename,
	     const bool &thinning,
	     const bool &curved,
	     const bool &slant,
             const bool &stackinput,
	     int str_length) {
    
    std::ostringstream FileName;
    for (unsigned int i=0; i<str_length; i++) {

	if (filename [i]==' ')
	    break;

	FileName << filename [i];
    }

    gThinned = thinning;
    gPlotter = TCorsikaPlotter (FileName.str ());
}



void cloda_ () {
}




