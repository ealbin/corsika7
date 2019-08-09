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
  crs::TSubBlock subBlock(Data, gPlotter.IsThinned());
  
  switch (subBlock.GetBlockType()) {
    
      case crs::TSubBlock::eRUNH:	   
        gPlotter.SetRunHeader (subBlock);    // New Run
        gPlotter.Write(subBlock);
        break;
        
      case crs::TSubBlock::eEVTH:
        gPlotter.SetShowerHeader (subBlock); // New Event
        gPlotter.Write(subBlock);
        break;
        
      case crs::TSubBlock::eEVTE:
        gPlotter.SetShowerTrailer (subBlock);// Event End
        gPlotter.Write(subBlock);
        break;
        
      case crs::TSubBlock::eRUNE:
      case crs::TSubBlock::eLONG:
        gPlotter.Write(subBlock);
        break;
        
      default:
        break;
  }
  
}



void 
inida_ (const char* filename,
	const bool& thinning,
	const bool& curved,
	const bool& slant,
	const bool& stackinput,
	const bool& /*preshower*/,
	int str_length) 
{
  if (str_length==0) 
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
  
  gPlotter.SetOutputOption(TPlotter::eCORSIKA);
  //gPlotter.SetOutputOption(TPlotter::eROOT);
  gPlotter.SetDirectory(directory.str());
  gPlotter.SetFileName(fileName);
  gPlotter.SetThinning(thinning);
  gPlotter.SetSlant(slant);
  gPlotter.SetCurved(curved);
  gPlotter.SetStackInput(stackinput);
  // gPlotter.SetCollectingPlane();   // TODO: DEFINE PLANE FOR PARTICLE COLLECTION AND PASS IT TO "TPlotter"
  
  gPlotter.Init();
  gPlotter.Welcome();
}



void 
cloda_ () 
{
  gPlotter.Close();
}



void 
track_ (const crs::CParticle& pre, const crs::CParticle& post) 
{
  gPlotter.AddTrack(pre, post);
}

void 
interaction_ (const crs::CInteraction& /*interaction*/) 
{
}





