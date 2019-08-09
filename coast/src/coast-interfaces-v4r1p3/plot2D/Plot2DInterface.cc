#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>

#include <crs2r/TC2R.h>

#include <TPlotter.h>

#include <sstream>
#include <iostream>
using namespace std;


///  Global object to handle all ROOT I/O operations and data conversions
TPlotter gPlotter;

void
wrida_ (const CREAL *Data) 
{
  crs::TSubBlock subBlock (Data, gPlotter.IsThinned());
  
  switch (subBlock.GetBlockType ()) {
    
      case crs::TSubBlock::eRUNH:	   
        cout << "COAST> ****** RUNH ******** " << endl;
        gPlotter.SetRunHeader(subBlock);
        break;
        
      case crs::TSubBlock::eEVTH:
        cout << "COAST> ****** EVTH ******** " << endl;
        gPlotter.SetShowerHeader(subBlock);
        break;
        
      case crs::TSubBlock::eEVTE:
        cout << "COAST> ****** EVTE ******** " << endl;
        gPlotter.SetShowerTrailer(subBlock);
        gPlotter.Write();
        break;
        
      case crs::TSubBlock::eRUNE:
        cout << "COAST> ****** RUNE ******** " << endl;
        break;
        
      default:
        break;
  }
  
}



void 
inida_ (const char* filename,
	const bool& thinning,
	const bool& /*curved*/,
	const bool& /*slant*/,
	const bool& /*stackinput*/,
	const bool& /*preshower*/,
	int str_length) 
{
  if (str_length==0) 
    str_length = 9;
  
  std::ostringstream FileName;
  for (unsigned int i=0; i<str_length; i++) {
    char ch = filename [i];
    if (ch==' ')
      break;
    FileName << ch;
  }
  
  gPlotter.SetThinning(thinning);
}



void cloda_ () 
{    
}



void
track_ (const crs::CParticle& pre, const crs::CParticle& post) 
{
  gPlotter.AddTrack(pre, post);
}

void interaction_ (const crs::CInteraction& /*interaction*/) 
{
}



