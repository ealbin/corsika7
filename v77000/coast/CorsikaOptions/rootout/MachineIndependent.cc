#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>

#include <crsRead/TSubBlockIO.h>

#include <sstream>
#include <fstream>
#include <iostream>






// definition of CORSIKA binary data structure
const int gNBlockSizeInfo = 4;



// global variables for i/o operations
std::ofstream gBinaryFile;
crsRead::TSubBlockIO gSubBlock;




/*
  Data is one CORSIKA data-block constining of 21 SubBlocks.
  A SubBlock can be:
   - thinned mode:     39 (Particles) * 8 (ENTRIES) * 4 (BYTES) 
   - not-thinned mode: 39 (Particles) * 7 (ENTRIES) * 4 (BYTES) 
 */
void
wrida_ (const CREAL* Data) 
{    
    gSubBlock.SetBuffer (Data);
    gBinaryFile << gSubBlock;
}





void
inida_ (const char* filename,
        const int& thinning,
        const int& /*curved*/,
        const int& /*slant*/,
        const int& /*stackinput*/,
        const int& /*preshower*/,
        int str_length) 
{    
    std::ostringstream FileName;
    for (int i=0; i<str_length; i++) {

	if (filename [i]==' ')
	    break;

	FileName << filename [i];
    }
    
    std::cout << " ######################################################### "
	      << std::endl
	      << " Write machine independent file: " << FileName.str ()
	      << " thinned: " << thinning
	      << std::endl;
    
    gBinaryFile.open (FileName.str ().c_str (),
                      std::ios::binary | std::ios::out | std::ios::trunc);


    // create the CORSIKA SubBlock, that will handle all further streaming
    gSubBlock.InitSubBlockSize(thinning, gNBlockSizeInfo);
}





void
cloda_ () 
{
    gSubBlock.FinishBlock(gBinaryFile); // fill up last block
    gBinaryFile.flush ();
    gBinaryFile.close ();
}



void 
interaction_ (const crs::CInteraction& /*interaction*/)
{
}

void
track_ (const crs::CParticle& /*pre*/, const crs::CParticle& /*post*/)
{
} 
