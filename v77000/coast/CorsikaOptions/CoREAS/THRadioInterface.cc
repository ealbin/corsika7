/***************************************************************************
 *                                                                         *
 *  Copyright and any other appropriate legal protection of these          *
 *  computer programs and associated documentation are reserved in         *
 *  all countries of the world.                                            *
 *                                                                         *
 *  These programs or documentation may not be reproduced by any method    *
 *  without prior written consent of Karlsruhe Institute of Technology     *
 *  delegate. Commercial and military use are explicitly forbidden.        *
 *                                                                         *
 *  Karlsruhe Institute of Technology welcomes comments concerning the     *
 *  code but undertakes no obligation for maintenance of the programs,     *
 *  nor responsibility for their correctness, and accepts no liability     *
 *  whatsoever resulting from the use of its programs.                     *
 *                                                                         *
 ***************************************************************************/

#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>

#include <TPlotter.h>

#include <sstream>
#include <iostream>
using namespace std;


///  Global object to handle all ROOT I/O operations and data conversions
TPlotter gPlotter;

void 
wrida_ (const CREAL *Data) 
{
  crs::TSubBlock subBlock(Data, gPlotter.IsThinned());
  
  switch (subBlock.GetBlockType()) {
    
      case crs::TSubBlock::eRUNH:	   
        gPlotter.SetRunHeader (subBlock);
        break;
        
      case crs::TSubBlock::eEVTH:
        gPlotter.SetShowerHeader (subBlock);
        gPlotter.Init();
        break;
        
      case crs::TSubBlock::eEVTE:
        gPlotter.SetShowerTrailer (subBlock);
        gPlotter.Write();
        break;
        
      case crs::TSubBlock::eRUNE:
        break;
        
      default:
        break;
  }
  
}


void 
inida_ (const char *filename,
	const int &thinning,
	const int &curved,
	const int &slant,
	const int &stackinput,
	const int &preshower,
	int str_length) 
{  
  if (str_length == 0) 
    str_length = 9;
  
  std::ostringstream FileName;
  for (int i=0; i<str_length; i++) {
    char ch = filename [i];
    if (ch==' ')
      break;
    FileName << ch;
  }
  std::ostringstream directory;
  if (FileName.str().rfind('/') != std::string::npos) {
    directory << FileName.str().substr(0,FileName.str().rfind('/'));
    directory << "/";
  }
  std::string fileName;
  if (FileName.str().rfind('/') != std::string::npos) {
    fileName =
      FileName.str().substr(FileName.str().rfind('/')+1,
			    FileName.str().size()-FileName.str().rfind('/')-1);
  } else {
    fileName = FileName.str();
  }
  
  gPlotter.SetDirectory(directory.str());
  gPlotter.SetFileName(fileName);
  gPlotter.SetThreadName(fileName.size()>9 ? fileName.substr(9,fileName.size()-9) : ""); // save name of parallel thread for debug output
  gPlotter.SetThinning(thinning);
  gPlotter.SetSlant(slant);
  gPlotter.SetCurved(curved);
  gPlotter.SetStackInput(stackinput);
  gPlotter.SetPreshower(preshower);
  
  gPlotter.Welcome();
  
}


void cloda_ () 
{
}


void 
track_ (const crs::CParticle& pre, const crs::CParticle& post) 
{
  gPlotter.AddTrack (pre, post);
}


void 
interaction_ (const crs::CInteraction& interaction) 
{
  gPlotter.AddInteraction(interaction);
}

void mpinodedone_ ()
{
  gPlotter.CommunicateNodeResults();	
}



void tabularizedatmosphere_(const int &nPoints, const double* height, const double* refractiveIndex)
{
	// copy the atmosphere to the coreas part - this wastes memory
	gPlotter.FillRefractivityTablesFromGDAS(nPoints, height, refractiveIndex);
}
