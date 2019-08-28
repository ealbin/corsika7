#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MRunHeader.h>
#include <crs/MEventEnd.h>

#include <crs/CParticle.h>
#include <crs/CInteraction.h>

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
  crs::TSubBlock subBlock (Data, gPlotter.IsThinned ());
  
  switch (subBlock.GetBlockType ()) {
    
  case crs::TSubBlock::eRUNH:	   
    gPlotter.SetRunHeader (subBlock);
    break;
    
  case crs::TSubBlock::eEVTH:
    gPlotter.SetShowerHeader (subBlock);
    gPlotter.Init ();
    break;
    
  case crs::TSubBlock::eEVTE:
    gPlotter.SetShowerTrailer (subBlock);
    gPlotter.Write ();
    break;
    
  default:
    break;
  }
}



void 
inida_ (const char *filename,
	const bool &thinning,
	const bool &curved,
	const bool &slant,
	const bool &stackinput,
	const bool &preshower,
	int str_length) 
{
  gPlotter.Set(thinning, curved, slant, stackinput, preshower);
}



void 
cloda_ () 
{
  gPlotter.Close ();
}



void 
track_(const crs::CParticle& pre, const crs::CParticle& post) 
{    
  gPlotter.AddTrack (pre, post);
}


void 
interaction_ (const crs::CInteraction& interaction) 
{
  gPlotter.AddInteraction(interaction);
}

