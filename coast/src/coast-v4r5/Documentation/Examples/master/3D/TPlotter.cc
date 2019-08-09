#include <TPlotter.h>

#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TSystem.h>
#include <TObjArray.h>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>

using namespace crs;
using namespace std;


TPlotter::TPlotter() :
  fDebug(false),  
  fDrawInShowerCS(false),
  fDoNotDrawPrimaryTrack(true) {
  
  const string COAST_USER_LIB(gSystem->Getenv("COAST_USER_LIB"));
  const string COAST_DIR(gSystem->Getenv("COAST_DIR"));
  const string PWD(gSystem->Getenv("PWD"));
  
  if (!ParseConfig((PWD+"/COAST3DConfig.config").c_str())) {
    if (COAST_USER_LIB == "" ||
	!ParseConfig((COAST_USER_LIB+"/COAST3DConfig.config").c_str())) {
      if (COAST_DIR == "" ||
	  !ParseConfig((COAST_DIR+"/config/COAST3DConfig.config").c_str())) {
	cerr << "\n ***********************************************\n"
	     << " * COAST-3D: config error (COAST3DConfig.config not found) !\n";
	if (COAST_USER_LIB=="")
	  cerr << " *\n * Define COAST_USER_LIB environment variable ! \n";
	if (COAST_DIR=="")
	  cerr << " *\n * Define COAST_DIR environment variable ! \n";
	cerr << " *************************************************"
	     << endl;
	exit(2);
      } else {
	cout << "COAST> using config file: " << (COAST_DIR+"/config/COAST3DConfig.config") << endl;
      }
    } else {
      cout << "COAST> using config file: " << (COAST_USER_LIB+"/COAST3DConfig.config") << endl;
    }
  } else {
    cout << "COAST> using config file: " << (PWD+"/COAST3DConfig.config") << endl;
  }    
  Clear();
}

TPlotter::~TPlotter() 
{
  Clear();
  delete fFile;
}

bool 
TPlotter::ParseConfig(const string& configName) 
{
  ifstream config(configName.c_str());
  if (!config.good())
    return false;
  
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
    
    const string particleName = word;
    int id;
    config >> id;
    int color;
    config >> color;
    double minE;
    config >> minE;
    int maxDraw;
    config >> maxDraw;
    
    if (!config.good()) {
      error = true;
    }

    
    config.ignore(99999, '\n');

    if (fParticleType.count(id)) {
      cerr << " TPlotter::ParseConfig> Multiple definitions of particle Id=" << id
	   << endl;
      return false;
    }
    
    ParticleType type;
    type.name = particleName;
    type.color = color;
    type.minEnergy = minE;
    type.maxDraw = maxDraw;
    fParticleType[id] = type;
    fParticleType[id].lines = new TObjArray();
    fParticleType[id].lines->SetOwner(kFALSE);
    fParticleType[id].count = 0;
  }

  config.close();
  return true;
}


void 
TPlotter::Clear() 
{
  fPrimaryTrack = true;
  fFirstInteraction = false;
  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  fFirstInteractionDist = 0;

  fFirst = true;
  
  fEventNo = 0;
  fObservationLevel = 0;
  fHeightFirstInt = 0;
  fZenith = 0;
  fAzimuth = 0;
  fEnergy0 = 0;
  fPrimary = 0;
  
  fCosZenith = 0;
  fSinZenith = 0;
  fCosAzimuth = 0;
  fSinAzimuth = 0;
  
  for (map<int,ParticleType>::iterator i = fParticleType.begin(); 
       i != fParticleType.end(); ++i) {
    i->second.lines->Delete(); 
    i->second.lines = new TObjArray();
    i->second.lines->SetOwner(kFALSE);
    i->second.count = 0;     
  }
}



void 
TPlotter::SetShowerZenith (float zenith) 
{
  fZenith = zenith*rad;
  fCosZenith = cos (zenith);
  fSinZenith = sin (zenith);
}



void 
TPlotter::SetShowerAzimuth (float azimuth) 
{
  fAzimuth = azimuth*rad;
  fCosAzimuth = cos (azimuth);
  fSinAzimuth = sin (azimuth);
}


void 
TPlotter::SetRunHeader (const crs::MRunHeader &header) 
{
  int run = (int)header.GetRunID ();
  // create filename
  ostringstream fname;
  fname << "DAT" 
	<< setw (6) << setfill('0') << run;
  SetFileName (fname.str ());
}


void 
TPlotter::SetShowerHeader(const crs::MEventHeader &header) 
{
  // to flag the track of the primary particle
  fPrimaryTrack = true;
  if (fStackInput || fPreshower) {
    fPrimaryTrack = false;
    fFirstInteraction = true;
  }
  
  SetEvent (header.GetEventNumber ());
  SetPrimary ((int)header.GetParticleId ());
  SetShowerZenith (header.GetTheta ()*rad);
  SetShowerAzimuth (header.GetPhi ()*rad);    
  SetShowerEnergy (header.GetEnergy ()*GeV);
  SetHeightFirstInt (header.GetZFirst ()*cm);
  SetObservationLevel (header.GetObservationHeight(header.GetNObservationLevels ()-1)*cm);

  if (fDebug) {
    cout << " COAST> zenith=" << fZenith/deg
	 << " azimuth=" << fAzimuth/deg << endl;
  }

}
    
    
void 
TPlotter::SetShowerTrailer (const crs::MEventEnd &trailer) 
{
}


void 
TPlotter::Init() 
{
  ostringstream file_name;
  file_name << fFileName 
	    << "_"
	    << fEventNo
	    << "_3d.root";
    
  fFile = TFile::Open (file_name.str ().c_str (), "RECREATE");
  fFile->cd ();    

  ostringstream cname;
  cname << "ev_" << fEventNo; 

  fCanvas = new TCanvas(cname.str().c_str());
  fCanvas->SetView(new TView3D());
}


void 
TPlotter::Close () 
{
}


void 
TPlotter::Write () 
{
  fCanvas->cd();
  
  cout << " COAST> +++++++++++++++++++++++++++++++++++\n";
  for (map<int,ParticleType>::iterator i=fParticleType.begin(); 
       i != fParticleType.end(); ++i) {
    cout << " COAST> plotted " << setw(8) << i->second.count << " tracks of type \'" << i->second.name << "\'\n";
    if ( i->second.count>0 ) { 
      i->second.lines->Draw();
    }
  }
  cout << " COAST> +++++++++++++++++++++++++++++++++++ " << endl;
  
  if (fDebug) {
    cout << "COAST> set range: " << fXmin << " " << fYmin << " "<< fZmin
	 << " " << fXmax << " " << fYmax << " " << fZmax << endl;
  }

  fCanvas->GetView ()->SetRange (fXmin, fYmin, fZmin, fXmax, fYmax, fZmax);
  fCanvas->Write ("", TFile::kOverwrite);
  fFile->Close ();
    
  Clear();
}





// #define Rotate(x, y, sx, sy, cosa, sina) { 
inline
void 
TPlotter::Rotate(const double x, const double y, const double z,
		 double &sx, double &sy, double &sz,
		 const int inverse) 
{
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
TPlotter::AddInteraction(const crs::CInteraction &interaction) 
{
  if (fDebug) {
    cout << " COAST interaction : " << endl;
    interaction.Dump();
  }
  
  if (fPrimaryTrack) 
    fFirstInteraction = true;
  else 
    fFirstInteraction = false;
  
  fPrimaryTrack = false;
}


void 
TPlotter::AddTrack(const crs::CParticle &pre, 
		   const crs::CParticle &post) 
{
  if (fDebug) {
    cout << " COAST track pre(id/x/y/z) : "
	 << " " << pre.particleId
	 << " " << pre.x
	 << " " << pre.y 
	 << " " << pre.z 
	 << " post(id/x/y/z) : " << post.particleId
	 << " " << post.x
	 << " " << post.y
	 << " " << post.z
	 << endl;
  }


  /*
    Skip the track of the primary particle, which is in a different 
    reference system as the shower !!!
  */
  if (fPrimaryTrack && fDoNotDrawPrimaryTrack) {
    return;
  }

  if (fPrimaryTrack && fDrawInShowerCS)
    return;
  
  if (fFirstInteraction) {
    
    fFirstInteraction = false;
    SetFirstInteraction(pre.x, pre.y, pre.z);
       
    if (fDebug) {
      cout << " COAST: Primary track end-pos: x=" << fFirstInteractionX/m 
	   << " y=" << fFirstInteractionY/m 
	   << " z=" << fFirstInteractionZ/m << endl;
    }
  }


  double xPre = pre.x/km;
  double yPre = pre.y/km;
  double zPre = pre.z/km;
  // const double zPre = -pre.time*1e9;

  double xPost = post.x/km;
  double yPost = post.y/km;
  double zPost = post.z/km;
  //  const double zPost = -post.time*1e9;

  if (fDrawInShowerCS) {
    
    /*   ------------------------
	 convert to shower frame
    */
    Rotate (xPre - fFirstInteractionX, yPre - fFirstInteractionY, zPre - fFirstInteractionZ,
	    xPre, yPre, zPre, +1);
    zPre += fFirstInteractionDist;
    Rotate (xPost - fFirstInteractionX, yPost - fFirstInteractionY, zPost - fFirstInteractionZ,
	    xPost, yPost, zPost, +1);
    zPost += fFirstInteractionDist;
  }
  
  const double e = post.energy/GeV;
  int particleId = (int)pre.particleId;
  
  const bool fconex = (particleId < 0);
  if (fconex) {
    particleId = -particleId;
  }
  
  map<int,ParticleType>::iterator it = fParticleType.find(particleId);
  if (it == fParticleType.end())
    return;
  
  if (it->second.minEnergy>e)
    return;
  
  /*  
    cout << " RU: "
    << " " << particleId
    << " " << xPre 
    << " " << yPre 
    << " " << zPre 
    << " " << particleId
    << " " << xPost 
    << " " << yPost 
    << " " << zPost 
    << " " << fconex
    << endl;
*/  
  

  /*   ------------------------
       direction of particle in shower frame
  */
  const double dX = xPost - xPre; 
  const double dY = yPost - yPre; 
  const double dZ = zPost - zPre; 
  const double length = sqrt (dX*dX + dY*dY + dZ*dZ);
  
  if (length==0) {
    /*   ------------------------ 
	 This happens for particles that are NOT tracked by Corsika
	 e.g. pi0 -> decay immediatly
	 SKIP them
    */
    return;
  }	
  
  /*  double speed=length/(post.time-pre.time)/29.9792458E4 ;
  if (speed>1.0001) cout << " speed " << speed << " p " << particleId << " " <<length << " " << post.z << " " << fconex << endl  ;
  */
  if (it->second.count >= it->second.maxDraw)
    return;
  
  int color = it->second.color;
  if (fconex) 
    color = color + 5 ;

  TPolyLine3D* line = new TPolyLine3D(2);
  line->SetPoint(0, xPre, yPre, zPre);
  line->SetPoint(1, xPost, yPost, zPost);
  line->SetLineColor(color);
  it->second.lines->Add(line);
  it->second.count++;
  //  if (particleId<4) cout << " conex  " << fconex << endl  ;
  
  /*
    // add track end/start points
  TPolyMarker3D* step = new TPolyMarker3D(1, 20);
  step->SetPoint(0, xPre, yPre, zPre);
  step->SetMarkerColor(1);
  step->SetMarkerSize(0.3);
  it->second.lines->Add(step);

  TPolyMarker3D* step2 = new TPolyMarker3D(1, 20);
  step2->SetPoint(0, xPost, yPost, zPost);
  step2->SetMarkerColor(2);
  step2->SetMarkerSize(0.3);
  it->second.lines->Add(step2);
  */
  
  if (fFirst) {

    fXmin = std::min (xPre, xPost);
    fYmin = std::min (yPre, yPost);
    fZmin = std::min (zPre, zPost);

    fXmax = std::max (xPre, xPost);
    fYmax = std::max (yPre, yPost);
    fZmax = std::max (zPre, zPost);

    fFirst = false;

  } else {

    if (std::min (xPre, xPost)<fXmin)
      fXmin = std::min (xPre, xPost);
    if (std::min (yPre, yPost)<fYmin)
      fYmin = std::min (yPre, yPost);
    if (std::min (zPre, zPost)<fZmin)
      fZmin = std::min (zPre, zPost);

    if (std::max (xPre, xPost)>fXmax)
      fXmax = std::max (xPre, xPost);
    if (std::max (yPre, yPost)>fYmax)
      fYmax = std::max (yPre, yPost);
    if (std::max (zPre, zPost)>fZmax)
      fZmax = std::max (zPre, zPost);

  }    
}


void 
TPlotter::SetFirstInteraction(const double x, const double y, const double z)
{
  fFirstInteractionX = x;
  fFirstInteractionY = y;
  fFirstInteractionZ = z;  
  cout << " COAST-SetFirstInteraction: x,y,z=" << x << " " << y << " " << z
       << endl;
}
