#include <TPlotter.h>

#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>
using namespace crs;

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TASImage.h>
#include <TH3D.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;



typedef UInt_t CARD32;
typedef CARD32 ARGB32;

#if 1
#define MAKE_ARGB32(a,r,g,b)    ((( (CARD32)a)        <<24)|((((CARD32)r)&0x00FF)<<16)| \
                                     ((((CARD32)g)&0x00FF)<<8 )|(( (CARD32)b)&0x00FF))
#else
#define MAKE_ARGB32(a,r,g,b)    ((((a)&0x00FF)<<24)|(((r)&0x00FF)<<16)| \
                                     (((g)&0x00FF)<<8)|((b)&0x00FF))
#endif

#define MAKE_ARGB32_GREY(a,l)   (((a)<<24)|(((l)&0x00FF)<<16)|          \
                                     (((l)&0x00FF)<<8)|((l)&0x00FF))
#define ARGB32_ALPHA8(c)        (((c)>>24)&0x00FF)
#define ARGB32_RED8(c)          (((c)>>16)&0x00FF)
#define ARGB32_GREEN8(c)        (((c)>>8 )&0x00FF)
#define ARGB32_BLUE8(c)         ( (c)     &0x00FF)
#define ARGB32_CHAN8(c,i)       (((c)>>((i)<<3))&0x00FF)
#define MAKE_ARGB32_CHAN8(v,i)  (((v)&0x0000FF)<<((i)<<3))

#define ARGB32_ALPHA16(c)       ((((c)>>16)&0x00FF00)|0x00FF)
#define ARGB32_RED16(c)         ((((c)>>8)&0x00FF00)|0x00FF)
#define ARGB32_GREEN16(c)       (( (c)    &0x00FF00)|0x00FF)
#define ARGB32_BLUE16(c)        ((((c)<<8)&0x00FF00)|0x00FF)
#define ARGB32_CHAN16(c,i)      ((ARGB32_CHAN8(c,i)<<8)|0x00FF)
#define MAKE_ARGB32_CHAN16(v,i) ((((v)&0x00FF00)>>8)<<((i)<<3))


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   - first find the point of the first interaction, and define the range
     where the shower should get histogramed

   - rotate all particles in shower coordinate system, with (0,0,0) at 
     the first interaction

   - define planes in equal slant depth bins 

   - fill histograms 
   
 */

TPlotter::TPlotter() 
{
  fVerbosityLevel = 0;
  fFromGround = false;
  fPlotXmaxLevel = false;
  fPlotObservationLevels = false;
  fTopOfAtmosphere = 112.8 * km;  // from CORSIKA
  ReadInit();
  Clear ();
}



TPlotter::~TPlotter () 
{  
}



void
TPlotter::Clear () 
{
  fPrimaryTrack = true;
  fPrimaryTrackFirst = true;
  fFirstParticleX = 0;
  fFirstParticleY = 0;
  fFirstParticleZ = 0;

  fAxisX = 0;
  fAxisY = 0;
  fAxisZ = 0;

  fCosZenith = 0;
  fSinZenith = 0;
  fCosAzimuth = 0;
  fSinAzimuth = 0;
  
  fEventNo = 0;
  fObservationLevel = 0;
  fHeightFirstInt = 0;
}

void 
TPlotter::ReadInit() 
{ 
  const string COAST_USER_LIB((gSystem->Getenv("COAST_USER_LIB")?gSystem->Getenv("COAST_USER_LIB"):"none"));
  const string COAST_DIR((gSystem->Getenv("COAST_DIR")?gSystem->Getenv("COAST_DIR"):"none"));
  const string PWD(gSystem->Getenv("PWD"));

  if (!ParseConfig((PWD+"/COAST2DConfig.config").c_str())) {
    if (COAST_USER_LIB == "" ||
	!ParseConfig((COAST_USER_LIB+"/COAST2DConfig.config").c_str())) {
      if (COAST_DIR == "" ||
	  !ParseConfig((COAST_DIR+"/config/COAST2DConfig.config").c_str())) {
	cerr << "\n ***********************************************\n"
	     << " * COAST-2D: config error (COAST2DConfig.config not found) !\n";
	if (COAST_USER_LIB=="")
	  cerr << " *\n * Define COAST_USER_LIB environment variable ! \n";
	if (COAST_DIR=="")
	  cerr << " *\n * Define COAST_DIR environment variable ! \n";
	cerr << " *************************************************"
	     << endl;
	exit(2);
      } else {
	cout << "COAST> using config file: " << (COAST_DIR+"/config/COAST2DConfig.config") << endl;
      }
    } else {
      cout << "COAST> using config file: " << (COAST_USER_LIB+"/COAST2DConfig.config") << endl;
    }
  } else {
    cout << "COAST> using config file: " << (PWD+"/COAST2DConfig.config") << endl;
  }    
  Clear();
}


bool 
TPlotter::ParseConfig(const string& configName) 
{
  ifstream config(configName.c_str());
  if (!config.good())
    return false;

  fXstretch = 1;
  fAutoRange = false;
  fWeightBoostMuons = 1;
  fWeightBoostHadrons = 1;
  fExtension = "png";

  map<string,bool> checker;

  bool error = false;
  while(config.good()) {
    
    string word;
    config >> word;
    if (word.length()<1 ||
	word == "" ||
	word[0] == '#') {
      config.ignore(99999, '\n');
      continue;
    }
    
    if (word == "auto-range") {
      int flag;
      config >> flag;
      fAutoRange = flag > 0;
    } else if (word == "xmin") {
      config >> fXstart;
      fXstart *= km;
    } else if (word == "xmax") {
      config >> fXstop;
      fXstop *= km;
    } else if (word == "ymin") {
      config >> fYstart;
      fYstart *= km;
    } else if (word == "ymax") {
      config >> fYstop;
      fYstop *= km;
    } else if (word == "stretch_x") {
      config >> fXstretch;
    } else if (word == "pxl_per_km") {
      config >> fPixelPerKm;
      fPixelPerKm /= km;
    } else if (word == "weight-boost-muons") {
      config >> fWeightBoostMuons;
    } else if (word == "weight-boost-hadrons") {
      config >> fWeightBoostHadrons;
    } else if (word == "bg-color") {
      double readRed=0, readGreen=0, readBlue=0;
      config >> readRed >> readGreen >> readBlue;
      fBgRed = readRed * 255;
      fBgGreen = readGreen * 255;
      fBgBlue = readBlue * 255;
      if (fBgRed<0 || fBgRed>255 ||
          fBgGreen<0 || fBgGreen>255 ||
          fBgBlue<0 || fBgBlue>255) {
        cout << "COAST> config: bg-color out of range: [0, 1]" << endl;
        error = true;
      }
    } else if (word == "particle") {
      string name;
      int id;
      double red=0, green=0, blue=0;
      config >> name >> id >> red >> green >> blue;
      if (fParticleType.count(id)) {
        cerr << " TPlotter::ParseConfig> Multiple definitions of particle Id=" << id
             << endl;
        return false;
      }      
      ParticleType type;
      type.name = name;
      type.id = id;
      type.red = red;
      type.green = green;
      type.blue = blue;
      fParticleType[id] = type;
    } else if (word == "output-extension") {
      config >> fExtension;
    }
    
    checker[word] = true;

    if (!config.good()) {
      error = true;
    }    
    config.ignore(99999, '\n');
  }
  
  if (!checker.count("xmin") ||
      !checker.count("xmax") ||
      !checker.count("ymin") ||
      !checker.count("ymax") ||
      !checker.count("pxl_per_km") ||
      !checker.count("bg-color")) {
    cout << "COAST> Missing data in config file" << endl;
    for (map<string,bool>::const_iterator i = checker.begin(); i != checker.end(); ++i) {
      cout << "COAST> \"" << i->first << "\"" << endl;
    }
    error = true;
    return false;
  }

  config.close();  
  return !error;
}

/*
  fXstart = -1.5*km;
  fXstop  =  1.5*km;
  
  fYstart = -0.1*km;
  fYstop  = 33.*km;
  double width  = fXstop-fXstart;
  double height = fYstop-fYstart;
  double pixelPerKm = 30.;
  fNX = int(width/km*pixelPerKm);
  fNY = int(height/km*pixelPerKm);
  fMaxRed   = 0;
  fMaxGreen = 0;
  fMaxBlue  = 0;
  cout << "COAST> w=" << width << " h=" << height << " nx=" << fNX << "ny=" << fNY << endl; 
*/

void
TPlotter::SetShowerZenith(float zenith) 
{
  cout << "COAST> zenith/deg: " << zenith/deg << endl;
  fCosZenith = cos(zenith);
  fSinZenith = sin(zenith);
}



void
TPlotter::SetShowerAzimuth(float azimuth) 
{
  cout << "COAST> azimuth/deg: " << azimuth/deg << endl;
  fCosAzimuth = cos(azimuth);
  fSinAzimuth = sin(azimuth);
}


void
TPlotter::SetRunHeader(const crs::MRunHeader &header) 
{
  SetRun ((int)header.GetRunID ());
  // create filename
  ostringstream fname;
  fname << "DAT" 
	<< setw (6) << setfill('0') << fRunNo;
  SetFileName(fname.str ());
}


void 
TPlotter::SetShowerHeader(const crs::MEventHeader &header) 
{
  // to flag the track of the primary particle
  fPrimaryTrack = true;
  fPrimaryTrackFirst = true;
  fFirstParticleX = 0;
  fFirstParticleY = 0;
  fFirstParticleZ = 0;
  
  SetEvent (header.GetEventNumber ());
  SetShowerAxis(header.GetTheta ()*rad, header.GetPhi ()*rad);
  SetHeightFirstInt (header.GetZFirst ()*cm);
  SetObservationLevel (header.GetObservationHeight (header.GetNObservationLevels ()-1)*cm);
  
  if (fHeightFirstInt<0) {
    if (fVerbosityLevel>=1) {
      cout << "COAST> height of first interaction is NEGATIVE \n";
    }
    fHeightFirstInt *= -1;
  }
  cout  << "COAST> firstint " << fHeightFirstInt/m << "\n";
  cout  << "COAST> obslev " << fObservationLevel/m << "\n";
  if (fAutoRange) {
    fYstart = fObservationLevel;
    fYstop = fHeightFirstInt;
    const double height = (fYstop - fYstart) * 0.02;
    fYstart -= height;
    fYstop += height;
  }
  Init();
}


void
TPlotter::SetShowerTrailer(const crs::MEventEnd &trailer) 
{
  fXmax = trailer.GetXmax();
  fPrimaryTrack = true;
  fPrimaryTrackFirst = true;
}


void 
TPlotter::Init() 
{
  cout << "COAST> **************** Init ****************** " << endl;

  // init bitmap and arrays
  fXstop *= fXstretch;
  fXstart *= fXstretch;
  const double width  = fXstop - fXstart;
  const double height = fYstop - fYstart;
  fNX = int(width * fPixelPerKm);
  fNY = int(height * fPixelPerKm);
  fMaxRed   = 0;
  fMaxGreen = 0;
  fMaxBlue  = 0;

  const int size = fNX*fNY;
  fRed   = new double[size];
  fGreen = new double[size];
  fBlue  = new double[size];
  for (int i=0; i<size; i++) {
    fRed[i] = 0;
    fBlue[i] = 0;
    fGreen[i] = 0;
  }

  cout << "COAST> width=" << width/km << "km,x "
       << "height=" << height/km << "km, "
       << "nx=" << fNX << "px, "
       << "ny=" << fNY << "px, "
       << "px/km=" << fPixelPerKm*km << endl; 
}




void 
TPlotter::Write() 
{
  cout << "COAST> **************** Write ****************** " << endl;

  if (fVerbosityLevel>=2) {
    cout << "COAST> Write " << endl;
  }
  
 
  if (!fFromGround) {
    
    // const int bgColor = MAKE_ARGB32(255, fBgRed, fBgGreen, fBgBlue);
    
    cout << "COAST> maxRed="   << fMaxRed
	 << " maxGreen=" << fMaxGreen
	 << " maxBlue="  << fMaxBlue
	 << endl;
    
    double maxRed   = fMaxRed;
    double maxGreen = fMaxGreen;
    double maxBlue  = fMaxBlue;
    
    double maxTot = max(max(maxRed, maxGreen), maxBlue);
    //maxTot *= 0.001; // 1e10
    //maxTot *= 0.009; // 1e8
    //maxTot *= 0.04; // 1e6
    maxTot *= 0.2; // 1e4
    //maxTot *= 0.3; // 1e2
    maxRed   = maxTot;
    maxGreen = maxTot;
    maxBlue  = maxTot;

    double minRed   = 0;
    double minGreen = 0;
    double minBlue  = 0;
    
    // check color range
    for(int iX=0; iX<fNX; ++iX) {
      for(int iY=0; iY<fNY; ++iY) {
	const int bin = iY*fNX + iX;
	if ((fRed[bin] > 0 && fRed[bin] < minRed) || minRed==0) {
          minRed = fRed[bin];
        }
	if ((fGreen[bin] > 0 && fGreen[bin] < minGreen) || minGreen==0) {
          minGreen = fGreen[bin];
        }
	if ((fBlue[bin] > 0 && fBlue[bin] < minBlue) || minBlue==0) {
          minBlue = fBlue[bin];
        }
      }
    }
    
    cout << "COAST> minRed=" << minRed
	 << " minGreen=" << minGreen
	 << " minBlue="  << minBlue
	 << endl;
        
    double minTot = min(min(minRed, minGreen), minBlue);
    minTot = 0;
    minRed   = minTot;
    minGreen = minTot;
    minBlue  = minTot;
    
    TASImage* img = new TASImage(fNX, fNY);
    UInt_t* data = img->GetArgbArray();
    for(int iX=0; iX<fNX; ++iX) {
      for(int iY=0; iY<fNY; ++iY) {
	
	const int bin = (fNY-1-iY)*fNX + iX;
	const int binOut = iY*fNX + iX;

	const double red   = (maxRed   ? (fRed[bin]-minRed)/maxRed : 0.);
	const double green = (maxGreen ? (fGreen[bin]-minGreen)/maxGreen : 0.);
        const double blue  = (maxBlue  ? (fBlue[bin]-minBlue)/maxBlue : 0.);

        /*
          // works, but looks funny
        const double red   = (log(fRed[bin])-log(minRed)) / (log(maxRed)-log(minRed));
        const double green = (log(fGreen[bin])-log(minGreen)) / (log(maxGreen)-log(minGreen));
        const double blue  = (log(fBlue[bin])-log(minBlue)) / (log(maxBlue)-log(minBlue));
        */

	int redRGB   = (2*fBgRed   + int(red*255) - int(green*255) - int(blue*255)) / 2;
	int greenRGB = (2*fBgGreen + int(green*255) - int(red*255) - int(blue*255)) / 2;
	int blueRGB  = (2*fBgBlue  + int(blue*255) - int(green*255) - int(red*255)) / 2;

	if (redRGB<0)     redRGB = 0; 
	if (redRGB>255)   redRGB = 255;
	if (blueRGB<0)    blueRGB = 0; 
	if (blueRGB>255)  blueRGB = 255;
	if (greenRGB<0)   greenRGB = 0; 
	if (greenRGB>255) greenRGB = 255;
	
	const int color = MAKE_ARGB32(255, redRGB, greenRGB, blueRGB);
	data[binOut] = color;
      
      }
    }
  
    const double DY = (fYstop - fYstart) / fNY;
    
    if (fPlotXmaxLevel) {
      
      // -------------------------------------------------------------------
      // plot xmax level
      int yXmax  = int((heigh_(fXmax)-fYstart)/DY);
      cout << "COAST> yXmax=" << yXmax << endl;
      for (int iX=0; iX<fNX; ++iX) {
      int binOut = (fNY-1-yXmax)*fNX + iX;
      data[binOut] = MAKE_ARGB32(255, 255, 50, 50);
      binOut = (fNY-yXmax)*fNX + iX;
      data[binOut] = MAKE_ARGB32(255, 255, 50, 50);
      }
      
    }
    
    
    if (fPlotObservationLevels) {
      
      // -------------------------------------------------------------------
      // plot observation level
      int yObs0  = int((0.*m-fYstart)/DY); // sea level
      int yObs1  = int((110.*m-fYstart)/DY); // kascade
      int yObs2  = int((1400.*m-fYstart)/DY); // auger
      int yObs3  = int((3300.*m-fYstart)/DY); // south pole
      cout << "COAST> yObs=" << yObs1 << " " << fNY-1-yObs1 << endl;
      cout << "COAST> yObs=" << yObs2 << " " << fNY-1-yObs2 << endl;
      for (int iX=0;iX<fNX; ++iX) {
	int binOut = (fNY-1-yObs0)*fNX + iX;
	data[binOut] = MAKE_ARGB32(255, 0, 0, 0);
	int binOut1 = (fNY-1-yObs1)*fNX + iX;
	data[binOut1] = MAKE_ARGB32(255, 0, 0, 0);
	int binOut2 = (fNY-1-yObs2)*fNX + iX;
	data[binOut2] = MAKE_ARGB32(255, 0, 0, 0);
	int binOut3 = (fNY-1-yObs3)*fNX + iX;
	data[binOut3] = MAKE_ARGB32(255, 0, 0, 0);
      }
    }
    
    
  
    TCanvas* c = new TCanvas("c", "c", fNX, fNY);
    c->Range(fXstart, fYstart, fXstop, fYstop);
    c->SetFillColor(0);
    /*
      TH2F* frame = new TH2F("frame", "frame", 2, fXstart, fXstop, 2, fYstart, fYstop);
      frame->SetXTitle("X     [cm]");
      frame->GetXaxis()->SetTitleOffset(1.5);
      frame->GetXaxis()->SetTicks("+-");
      frame->SetYTitle("Y     [cm]");
      frame->GetYaxis()->SetTitleOffset(1.5);
      frame->GetYaxis()->SetTicks("+-");
      frame->Draw();
    */
    img->Draw();
    /*
      TLine *line = new TLine(fXstart, heigh_(fXmax), fXstop, heigh_(fXmax));
      line->SetLineColor(103);
      line->SetLineWidth(2);
      line->Draw("same");
    */
    //frame->Draw("same");
    
    
    ostringstream outName;
    outName << fFileName << "_" << fEventNo << "." << fExtension;
    c->Print(outName.str().c_str());
  }
  
  Clear (); // clean up
}




//#define Rotate(x, y, sx, sy, cosa, sina) { \
inline
void TPlotter::Rotate (double x, double y, double z,
		       double &sx, double &sy, double &sz,
		       int inverse) {

  // rotate around z by azimuth
  double sx_ =             x*fCosAzimuth + inverse * y*fSinAzimuth;
  double sy_ = - inverse * x*fSinAzimuth +           y*fCosAzimuth; 
  double sz_ =   z;
  
  // rotate around y` by zenith
  double sx__ =             sx_*fCosZenith - inverse * sz_*fSinZenith;
  double sy__ = sy_;
  double sz__ = + inverse * sx_*fSinZenith +           sz_*fCosZenith; 

  // rotate back around z`` by -azimuth 
  sx =             sx__*fCosAzimuth - inverse * sy__*fSinAzimuth;
  sy = + inverse * sx__*fSinAzimuth +           sy__*fCosAzimuth; 
  sz =   sz__;

}


void 
TPlotter::AddTrack(const crs::CParticle &pre, 
                   const crs::CParticle &post) 
{  
  int particleId0 = (int)pre.particleId;
    
  if (fVerbosityLevel>=9) {
    cout << "COAST> PRE> "
	 << " id=" << (int)pre.particleId
	 << " E=" << pre.energy
	 << " w=" << pre.weight
	 << " x=" << pre.x << "cm"
	 << " y=" << pre.y << "cm"
	 << " z=" << pre.z << "cm"
	 << " t=" << pre.time*s/ns << "ns"
	 << " X=" << pre.depth << "g/cm^2"
	 << "\n";
  }
  
  if (fVerbosityLevel>=10) {
    cout << "COAST> POST> "
	 << " id=" << (int)post.particleId
	 << " E=" << post.energy
	 << " w=" << post.weight
	 << " x=" << post.x << "cm"
	 << " y=" << post.y << "cm"
	 << " z=" << post.z << "cm"
	 << " t=" << post.time*s/ns << "ns"
	 << " X=" << post.depth << "g/cm^2"
	 << "\n";
    
  }

  const double pre_x = pre.x * fXstretch;
  const double pre_y = pre.y * fXstretch;
  const double post_x = post.x * fXstretch;
  const double post_y = post.y * fXstretch;
  
  /*
    Skip the track of the primary particle, which is in a different 
    reference system as the shower !!!
  */
  if (fPrimaryTrack) {
    
    if (fPrimaryTrackFirst) {
      fPrimaryTrackFirst = false;
      fPrimaryTrackEnergy = pre.energy;
      fPrimaryTrackId = (int)pre.particleId;
    }
    
    if (fPrimaryTrackId==(int)pre.particleId) {
      return;
    } 
    
    SetFirstInteraction(pre_x, pre_y, pre.z, pre.time);
    
    /*
      if (fHeightFirstInt<0) {
      fHeightFirstInt = pre.z;
      // hier weiter
      }
    */    
    
    if (fVerbosityLevel>=2) {
      cout << "COAST> Primary track end-pos: x=" << fFirstParticleX/m 
	   << " y=" << fFirstParticleY/m 
	   << " z=" << fFirstParticleZ/m << "\n";
    } else if (fVerbosityLevel>=1) {
      cout << "COAST> Primary track " << "\n";
    }
    
    fPrimaryTrack = false;
  }
   
  const double dX = post_x - pre_x; 
  const double dY = post_y - pre_y; 
  const double dZ = post.z - pre.z; 
  const double length = sqrt (dX*dX + dY*dY + dZ*dZ);
    
  if (length==0) {
    /*   ------------------------ 
	 This happens for particles that are NOT tracked by Corsika
	 e.g. pi0 -> decay immediatly
	 SKIP them
    */
    return;
  }	
  
  
  const int particleId = (int)pre.particleId;
  const map<int,ParticleType>::const_iterator it = fParticleType.find(particleId);
  if (it == fParticleType.end())
    return;

  if (fVerbosityLevel>=2) {
    cout << "COAST> *********** track"  << " w=" << pre.weight 
         << " x1=" << pre_x*cm/km << " z1=" << pre.z*cm/km
         << " x2=" << post_x*cm/km << " z2=" << post.z*cm/km
         << endl;
  }
  
  const double red = it->second.red;
  const double green = it->second.green;
  const double blue = it->second.blue;

  double weight = pre.weight;
  if (particleId>6)
    weight *= fWeightBoostHadrons;
  else if (particleId>3)
    weight *= fWeightBoostMuons;

  DrawLine(red, green, blue, 
           weight, pre_x*cm, pre.z*cm, post_x*cm, post.z*cm);  
}

void 
TPlotter::DrawLine(const double red, const double green, const double blue,
                   const double w,
                   const double x1, const double y1, 
                   const double x2, const double y2) 
{  
  const double DX = (fXstop - fXstart) / fNX;
  const double DY = (fYstop - fYstart) / fNY;
  
  // pixel address
  const double x1pix = (x1 - fXstart)  / DX;
  const double y1pix = (y1 - fYstart)  / DY;
  const double x2pix = (x2 - fXstart)  / DX;
  const double y2pix = (y2 - fYstart)  / DY;
  
  // current pixel
  int xCurr = int(x1pix);
  int yCurr = int(y1pix);
  
  double dx = x2pix - x1pix;
  double dy = y2pix - y1pix;
  const double dTot = sqrt(dx*dx + dy*dy);
  double dCurr = 0;

  int dxBin = (dx>0 ? +1 : -1);
  int dyBin = (dy>0 ? +1 : -1);
  if (dx==0) dxBin = 0;
  if (dy==0) dyBin = 0;
  const double d = sqrt(dx*dx + dy*dy);
  if (d==0) 
    return;
  dx /= d;
  dy /= d;
  
  //cout << " dx=" << dx << " dy=" << dy << endl;
  //cout << " NX=" << fNX << " NY=" << fNY << endl;
  
  
  while(true) {
    
    double xTest = xCurr + ( dx>0 ? 1 : 0);
    double yTest = yCurr + ( dy>0 ? 1 : 0);
    
    double faraway = 999999999;
    double d1 = faraway;
    double d2 = faraway;
    
    if (dx!=0) d1 = (xTest-x1pix) / dx;
    if (dy!=0) d2 = (yTest-y1pix) / dy;
    
    if (d1<0) d1 = faraway;
    if (d2<0) d2 = faraway;
    double dMin = min(d1,d2);
    dMin = min(dMin, dTot);
    double dPix = dMin - dCurr;

    /*
    cout << " xCurr=" << xCurr << " yCurr=" << yCurr 
         << " dTot="<< dTot << " d1=" << d1 << " d2=" << d2
         << " dCurr=" << dMin;
    */
  
    if (xCurr>=0 && yCurr>=0 &&
        xCurr<fNX && yCurr<fNY) {

      if (dPix*w<0) {
        cout << "COAST> smaler zero " << dPix << " " << w << endl;
      }
        
      const int bin = (yCurr*fNX + xCurr);
      
      //cout << " bin=" << bin;
      
      fRed[bin] += red*dPix*w;
      fGreen[bin] += green*dPix*w;
      fBlue[bin] += blue*dPix*w;
      if (fRed[bin]>fMaxRed) fMaxRed = fRed[bin];
      if (fGreen[bin]>fMaxGreen) fMaxGreen = fGreen[bin];
      if (fBlue[bin]>fMaxBlue) fMaxBlue = fBlue[bin];
    }

    //cout << endl;
    
    if (d1>dTot && d2>dTot) { // finished      
      break;
    }
    
    if (d1>d2) {
      if (dy>0) yCurr++;
      else yCurr--;
    } else {
      if (dx>0) xCurr++;
      else xCurr--;
    }
    
    dCurr = dMin;
    
  }
}


void
TPlotter::SetFirstInteraction(double x, double y, double z, double t) 
{
  fFirstParticleX = x;
  fFirstParticleY = y;
  fFirstParticleZ = z;
  fFirstParticleTime = t;
}


void
TPlotter::SetShowerAxis(double zenith, double azimuth) 
{  
  fCosZenith = cos(zenith);
  fSinZenith = sin(zenith);
  fCosAzimuth = cos(azimuth);
  fSinAzimuth = sin(azimuth);
  fAxisZ = fCosZenith;
  fAxisX = fSinZenith * fCosAzimuth;;
  fAxisZ = fSinZenith * fSinAzimuth;
}



