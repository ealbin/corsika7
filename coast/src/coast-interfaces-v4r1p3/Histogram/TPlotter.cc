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

/*
  RU, Do 23. Aug 12:42:47 CEST 2007
  TODO/CHECK
  - photon primaries are first tracked to their first interaction, then the
    Event-Header is passed to COAST -> COAST follows the first secondary to 
    its first interaction
  -> nuclear primaries first pass their Event-Header to COAST and are then 
     tracked to their first interaction.
  
  RU, Sa 11. Aug 19:06:01 CEST 2007 
  - changed security offset of last layer above the observation level
    from 1m to 1cm (1mm is too close and leads to a loss of particles already)
  - changed security offset for STACKIN from 1mm to 0.5mm (0,1mm is not enough)
  
  RU, Do 9. Aug 11:54:36 CEST 2007
  - fixed missing weight initialisation for weighted mean position calculations
  - fixed streaming of characters (version numbers + particle types)
  - added support for STACKIN
  - generally lift last histogram-layer 1m above the CORSIKA observation level, to
    assure getting the full shower particle data
  - since times, energies and radii are histogrammed lograithmically:
    added special treatment for values below zero 
    (->unphysical, log impossible, numerical problem)
    all data goes into the underflow-bins
  - started general cleanup of TPlotter data members and functions
  - added variances of mean <x> and <y> position
  - added division by zero check in calculation of mean <x> and <y>

*/

#include <TPlotter.h>
#include <TCorsika.h>

#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>
using namespace crs;

#include <COASTconfig.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdlib>
using namespace std;





/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   - first find the point of the first interaction, and define the range
     where the shower should get histogramed

   - rotate all particles in shower coordinate system, with (0,0,0) at 
     the first interaction

   - define planes in equal slant depth bins 

   - fill histograms 
   
 */

TPlotter::TPlotter() {
  
  //fVerbosityLevel = 10; // debugging
  fVerbosityLevel = 0;
  
  fMaxInvalids = 10000;
  fTimeErrorMargin = 0.0;
  
  fPrimaryTrack = true;
  fFirstInteraction = false;
  
  fSlant = false;
  fCurved = false;
  fStackInput = false;
  fPreshower = false;

  fTimeUnderflow = 1.e-16; // very close to zero time
  fCountInvalid = 0;
  fCountErrors = 0;
  
  fCORSIKA = 0;
  
  InitParticles();
  
  Clear();
  fDirectory = "";
}


TPlotter::~TPlotter() {
  fParticles.clear();
  Clear();
  delete fFile;
}


void 
TPlotter::Welcome() const 
{
  cout << "\n"
       << "COAST-Histogram> *************************************************\n"
       << "COAST-Histogram> **                                               \n"
       << "COAST-Histogram> ** You are using COAST-histogramming THRadio      \n"
       << "COAST-Histogram> **    COAST version: " << VERSION << " \n"
       << "COAST-Histogram> **  THRadio version: " << HISTOGRAMING_VERSION << " \n"
       << "COAST-Histogram> **                                               \n"
       << "COAST-Histogram> **  THIN/CURVED/SLANT/STACKIN: " << fThinning 
       << "/" << fCurved << "/" << fSlant << "/" << fStackInput << "\n"
       << "COAST-Histogram> **                                               \n"
       << "COAST-Histogram> *************************************************\n\n"
       << endl;
}




void 
TPlotter::Clear() 
{  
  fPrimaryTrack = true;
  fFirstInteraction = false;
  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  
  fAxisX = 0;
  fAxisY = 0;
  fAxisZ = 0;
  
  fCosZenith = 0;
  fSinZenith = 0;
  fCosAzimuth = 0;
  fSinAzimuth = 0;
  
  fXmax = 0;
  fEventNo = 0;
  fObservationLevel = 0;
  fHeightFirstInt = 0;
  fZenith = 0;
  fAzimuth = 0;
  fShowerEnergy = 0;
  fShowerPrimary = 0;
  
  fCountInvalid = 0;
  fCountErrors = 0;
  
  fHists.clear();
  fLayerDepth.clear();
  fLayerHeight.clear();
  fZPlane.clear();
  fNhist = 0;
  
  fCORSIKA = 0;
}


void
TPlotter::SetRunHeader(const crs::MRunHeader& header) 
{  
  ostringstream ssTmp; 
  ssTmp.str(""); ssTmp << header.GetVersion();  fCorsikaVersion = ssTmp.str();
  ssTmp.str(""); ssTmp << VERSION;              fCoastVersion = ssTmp.str();
  ssTmp.str(""); ssTmp << HISTOGRAMING_VERSION; fHistVersion = ssTmp.str();
  
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
  if (header.GetSkimmingIncidence()) {
    fSkimming = true;
    fSkimmingAltitude = header.GetSkimmingAltitude() * cm;
    cout << " COAST-Histogram> Detected skimming geometry with impact=" 
         << fSkimmingAltitude/m
         << " m" << endl;
  }
  
  // to flag the track of the primary particle
  fPrimaryTrack = true;
  fFirstInteraction = false;
  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  
  SetEvent(header.GetEventNumber());
  fShowerPrimary = (int)header.GetParticleId();
  
  SetShowerAxis(header.GetTheta()*rad, header.GetPhi()*rad);
  fShowerEnergy = header.GetEnergy()*GeV;
  
  fHeightFirstInt = header.GetZFirst()*cm;
  fObservationLevel = header.GetObservationHeight(header.GetNObservationLevels()-1)*cm;
  fObservationLevel += 10.*cm; // particles below the observation level are lost,
                               // therefore build in a 10cm safety margin to assure 
                               // getting the full particle data
  
  if (fHeightFirstInt<0) {
    if (fVerbosityLevel>=0) {
      cout << " COAST-Histogram> height of first interaction is NEGATIVE, corrected." 
           << endl;
    }
    fHeightFirstInt = fabs(fHeightFirstInt);
  }
  
  cout  << " COAST-Histogram> firstint " << fHeightFirstInt/m << "m\n"
        << " COAST-Histogram> obslev " << fObservationLevel/m 
        << "m, including safety offset of 10cm" 
        << endl;
}


void
TPlotter::SetShowerTrailer(const crs::MEventEnd &trailer) 
{
  fXmax = trailer.GetXmax()*g/cm/cm;
  if (!fSlant && fCosZenith!=0) fXmax /= fCosZenith;
  cout << " COAST-Histogram>  --------- shower end ------ \n"
       << " COAST-Histogram> Xmax: " << fXmax/g*cm*cm << "\n";
  if (fCountInvalid>0) 
    cout << " COAST-Histogram> invalid-counter: " << fCountInvalid << "\n";
  if (fCountErrors>0) 
    cout << " COAST-Histogram> error-counter: " << fCountErrors << "\n"
	 << "                  CHECK YOUR LOG FILE! \n";
  cout << flush;
}


void 
TPlotter::Init(int nBins, double maxRange) 
{
   if (!fSkimming) 
    fCORSIKA = new TCorsika(fZenith, 
			    fObservationLevel,
			    fSlant, fCurved);
  else 
    fCORSIKA = new TCorsika(fSkimmingAltitude,
			    fObservationLevel);
  
  if (fVerbosityLevel > 0) 
    fCORSIKA->SetVerbosityLevel(fVerbosityLevel);

  fNhist = nBins;
  fMaxShowerRange = maxRange; // for skmimming geometries
  
  fLayerDepth.clear();
  fLayerHeight.clear();
  fZPlane.clear();

  fCORSIKA->SetHeightOfFirstInteraction(fHeightFirstInt);

  if (fCurved && (not fSlant)) { // curved without slant is not supported in REAS
    cout << "\n\n"
         << " ###################################\n"
         << "   CURVED without SLANT is not supported   \n"
         << "    please switch on the SLANT option \n"
         << " ###################################\n\n"
         << endl;
    exit (11);
  }
  if (fStackInput || fPreshower) { // STACKIN OPTION NEEDS SPECIAL TREATMENT !!! 

    // define safety margin to prevent invalid times (in front of shower front)
    // caused by rounding problems
    fTimeErrorMargin = 3.3e-3*ns;	// ignore time errors corresponding to ~1 mm before shower front

    fPrimaryTrack = false;
    fFirstInteraction = true;
  }
  
  if (fVerbosityLevel>=10 && fCurved) {
    fCORSIKA->DumpAtmosphereTable();
  }
}


void 
TPlotter::CalculatePlanes() 
{  
  if (fVerbosityLevel>=10) {
    cout << "TPlotter::CalculatePlanes" << endl; 
  }
  
  // tell TCorsika the real starting point for histogramm layers
  //fCORSIKA->SetHeightOfFirstInteraction(fFirstInteractionZ);
  //fCORSIKA->SetHeightOfFirstInteraction(fCORSIKA->GetHeightOfSlantDepth(fFirstInteractionSlantDepth));

  double zFirst = fFirstInteractionSlantDepth; 
  double zLast = 0;
  if (fSkimming) { zLast = zFirst += fMaxShowerRange*g/cm/cm;  }
  else {           zLast = fCORSIKA->GetLastLevelSlantDepth(); }
    
  double zUnit = g/cm/cm;
  string ModeStr = "slant depth = ";
  string UnitStr = "g/cm^2";
  
  cout << " COAST-Histogram> start shower=" << zFirst/zUnit << " " << UnitStr
       << ", last histogram plane=" << zLast/zUnit << " " << UnitStr
       << endl;
  
  // --------------------------------------------------------------
  // decide where/what to plot
  
  
  ostringstream file_name;
  file_name << fDirectory
            << fFileName 
	    << "_"
	    << fEventNo
	    << ".hist.root";
  
  fFile = TFile::Open(file_name.str().c_str(), "RECREATE");
  fFile->cd ();
  
  const int nBins = fNhist;
  const double deltaSlice = (zLast-zFirst)/nBins;
  
  // --------------------------------------------------------------
  // create histogramms for each layer for this shower
  for (int iLayer=0; iLayer<nBins; ++iLayer) {
    
    const double zSlice = zFirst + deltaSlice*(iLayer+1);
    
    /* ----------------
       height of planes in shower frame, with (0/0/0) at core (OBS level)
    */
    double zPlane = 0;
    //case eHeight:
    //zPlane = (zFirst-zLast) - (zFirst-zLast)/(nBins-1)*i;
    //break;
      
    zPlane = fCORSIKA->GetDistanceOfSlantDepth(zSlice); 
    
    fLayerDepth.push_back(fCORSIKA->GetSlantDepthOfDistance(zPlane));
    fLayerHeight.push_back(fCORSIKA->GetVerticalHeightOfDistance(zPlane));
    fZPlane.push_back(zPlane);
    
    if (fVerbosityLevel>=2) {
      cout << " COAST-Histogram> plane " 
	   << " zSlice=" << zSlice/zUnit << UnitStr
	   << ", dist=" << zPlane/m << "m"
           << ", h=" << heigh_(zSlice/g*cm*cm*fCosZenith)*cm/km << "km"
	   << endl;
    }
    
    ostringstream h_title;
    h_title << ModeStr
	    << zSlice/zUnit
	    << " " << UnitStr;
    
    
    HistLayer histlayer; // this contains all histogramms for one layer
    
    // ------------------------------------------------------------------
    // loop over all particle types
    for (map<int,ParticleDef>::const_iterator iParticle = fParticles.begin();
	 iParticle != fParticles.end();
	 ++iParticle) {
      
      const string particleName = iParticle->second.name;
      const int color = iParticle->second.color;
      const int id = iParticle->first; // used as color
      
      // --------------------------------------------------
      // generate histograms
      HistDef hists;
      InitHistograms(hists);
      for (HistDefIterator iHist = hists.begin(); iHist != hists.end(); ++iHist) {
	ostringstream hName;
	hName << iHist->second->GetName() << "_" << particleName << "_" << iLayer;
	iHist->second->SetName(hName.str().c_str());
	iHist->second->SetLineColor(color);
	iHist->second->SetMarkerColor(color);
      }
      
      // -------------------------------------------------
      // create data object for one layer
      histlayer[id] = hists;
      
    } // end loop particle types
    
    // -------------------------------------------------
    // collect all layers
    fHists.push_back(histlayer); 
    
  } // end loop layers
  
}







void
TPlotter::Write() 
{  
  if (fVerbosityLevel>=2) {
    cout << " COAST-Histogram> Write " << endl;
  }

  fXmax /= g*cm*cm;

  // create the particle-string
  ostringstream particles_str;
  for (map<int,ParticleDef>::iterator iParticle = fParticles.begin();
       iParticle != fParticles.end();
       ++iParticle) {
    particles_str << iParticle->second.name << " ";
  }
  string particle_list = particles_str.str();
  char* particles = const_cast<char*>(particle_list.c_str());
  
  // create the histogram-string
  ostringstream hist_str;
  for (HistDefIterator iHist = fHists.begin()->begin()->second.begin();
       iHist != fHists.begin()->begin()->second.end(); ++iHist) {
    const string hName = string(iHist->second->GetName());
    hist_str << hName.substr(0, hName.find('_')) << " ";
  }
  string histogram_list = hist_str.str();
  char* histograms = const_cast<char*>(histogram_list.c_str());

  
  fFile->cd ();
  TTree *run = new TTree ("shower", "shower");
  run->Branch ("xmax"    , &fXmax            , "xmax/F");
  run->Branch ("evntno"  , &fEventNo         , "eventno/I");
  run->Branch ("obslevel", &fObservationLevel, "obslevel/D");
  run->Branch ("firstint", &fHeightFirstInt  , "firstint/D");
  run->Branch ("zenith"  , &fZenith          , "zenith/F");
  run->Branch ("azimuth" , &fAzimuth         , "azimuth/F");
  run->Branch ("energy"  , &fShowerEnergy    , "energy/F");
  run->Branch ("primary" , &fShowerPrimary   , "primary/I");
  run->Branch ("particles" , particles       , "particles/C");
  run->Branch ("histograms" , histograms     , "histograms/C");
  char *versionCorsika = const_cast<char*>(fCorsikaVersion.c_str());
  char *versionCoast   = const_cast<char*>(fCoastVersion.c_str());
  char *versionHist    = const_cast<char*>(fHistVersion.c_str());
  run->Branch ("corsika_version" , versionCorsika  , "corsika_version/C");
  run->Branch ("coast_version"   , versionCoast    , "coast_version/C");
  run->Branch ("hist_version"    , versionHist     , "hist_version/C");
  run->Fill();
    
    
  Float_t time, slantDepth, dist, age, rm, temp, density, pressure;
  Float_t height;
  TTree *tree = new TTree ("data_layer", "General information of the layer");
  tree->Branch ("time"  , &time,   "time/F");
  tree->Branch ("depth" , &slantDepth,  "depth/F");
  tree->Branch ("dist", &dist, "dist/F");
  tree->Branch ("height", &height, "height/F");
  tree->Branch ("age"   , &age,    "age/F");
  tree->Branch ("rm", &rm, "rm/F");
  tree->Branch ("temp", &temp, "temp/F");
  tree->Branch ("pressure", &pressure, "pressure/F");
  tree->Branch ("density", &density, "density/F");
  
  
  // create hist trees for each particle type
  map<int,TTree *> trees;
  for (map<int,ParticleDef>::iterator iParticle = fParticles.begin();
       iParticle != fParticles.end();
       ++iParticle) {
	
    int id = iParticle->first;
    string particleName = iParticle->second.name;

    ostringstream tree_name, tree_info;
    tree_name << "data_" << particleName;
    tree_info << "All the histogramms for " << particleName;
    
    fFile->cd();
    TTree *newTree = new TTree (tree_name.str().c_str(), 
				tree_info.str().c_str());
    
    trees[id] = newTree;
  }

  bool firstFill = true;
  for (int iLayer=0; iLayer<fNhist; ++iLayer) {

    time = (fZPlane[iLayer]-fZPlane[0]) / crs::cSpeedOfLight;
    dist = fZPlane[iLayer];  
    slantDepth = fLayerDepth[iLayer]; 
    height = fLayerHeight[iLayer];
    
    const double T = fCORSIKA->Temperature(height);	
    rm = fCORSIKA->MoliereRadius(slantDepth*fCosZenith, T)/m;
    
    static double g_earth = 9.81 * m/(s*s);
    double Pressure = slantDepth * fCosZenith * g_earth;
    static double MeanAirDensity = 28.95 * g/mol;
    static double R_gas = 8.314 * joule/mol/kelvin;
    double Density = Pressure / (R_gas*T/MeanAirDensity);
    
    temp = T;
    pressure = Pressure/bar;
    density = Density/g*cm*cm*cm;
    
    time /= ns;
    height /= km;
    dist /= km;
    slantDepth /= g*cm*cm;
    
    age = 3.0 / (1.0 + 2.0*fXmax/slantDepth); 
        
    tree->Fill();

    
    for (map<int,ParticleDef>::iterator iParticle = fParticles.begin();
	 iParticle != fParticles.end();
	 ++iParticle) {
      
      const int id = iParticle->first;
      const string particleName = iParticle->second.name;
      
      map<string, TH1*> treeMap; 
      for (HistDefIterator iHist = fHists[iLayer][id].begin(); 
	   iHist != fHists[iLayer][id].end(); ++iHist) {
	
	const string hname(string(iHist->second->GetName()));
	
	ostringstream bname;
	bname << hname.substr(0, hname.rfind('_'));
	
	if (firstFill) {
	  trees[id]->Branch(bname.str().c_str(), iHist->second->ClassName(), &(iHist->second));
	} else {
	  trees[id]->SetBranchAddress(bname.str().c_str(), &(iHist->second));
	}
      }
      
      trees[id]->Fill();
    }

    firstFill = false;
    
    // delete hists
    for (map<int,ParticleDef>::iterator iParticle = fParticles.begin();
	 iParticle != fParticles.end();
	 ++iParticle) {
      const int id = iParticle->first;
      for (HistDefIterator iHist = fHists[iLayer][id].begin(); 
	   iHist != fHists[iLayer][id].end(); ++iHist) {
	delete iHist->second;
      }
    }
  } // loop layers
  
  run->Write("", TFile::kOverwrite);
  tree->Write("", TFile::kOverwrite);
  for (map<int,TTree*>::iterator iTree = trees.begin();
       iTree!=trees.end(); ++iTree) {
    iTree->second->Write("", TFile::kOverwrite);
  }
  fFile->Close();
  
  // clean up
  Clear();
}




inline
void TPlotter::Rotate(double x, double y, double z,
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
TPlotter::AddTrack(const crs::CParticle& pre, 
                   const crs::CParticle& post) 
{  
  if (fVerbosityLevel>=9) {
    cout << " COAST-Histogram>   PRE> "
	 << " id=" << (int)pre.particleId
	 << " E=" << pre.energy
	 << " w=" << pre.weight
	 << " x=" << pre.x << "cm"
	 << " y=" << pre.y << "cm"
	 << " z=" << pre.z << "cm"
	 << " t=" << pre.time*s/ns << "ns"
	 << " X=" << pre.depth << "g/cm^2"
	 << endl;
  }
  
  if (fVerbosityLevel>=10) {
    cout << " COAST-Histogram>  POST> "
	 << " id=" << (int)post.particleId
	 << " E=" << post.energy
	 << " w=" << post.weight
	 << " x=" << post.x << "cm"
	 << " y=" << post.y << "cm"
	 << " z=" << post.z << "cm"
	 << " t=" << post.time*s/ns << "ns"
	 << " X=" << post.depth << "g/cm^2"
	 << endl;
  }
    
  /*
    Skip the track of the primary particle, which is in a different 
    reference system as the shower !!!
  */
  if (fPrimaryTrack) {
    if (fVerbosityLevel >= 2) {
      cout << " COAST: Primary track  " << endl;
    }    
    return;
  }

  if (fFirstInteraction) {

    fFirstInteraction = false;

    SetFirstInteraction(pre.x, pre.y, pre.z, pre.depth*g/cm/cm, pre.time);
    CalculatePlanes();
    
    if (fVerbosityLevel >= 2) {
      cout << " COAST: Primary track end-pos: x=" << fFirstInteractionX/m 
	   << " y=" << fFirstInteractionY/m 
	   << " z=" << fFirstInteractionZ/m << endl;
    }
  }
  
  const int particleId = abs((int)pre.particleId);
  if (!fParticles.count(particleId)) 
    return;
  
  /*   ------------------------
       convert to shower frame
       rotate around point of first interaction
  */
  double pre_shX, pre_shY, pre_shZ;
  Rotate(pre.x-fFirstInteractionX, pre.y-fFirstInteractionY, fFirstInteractionZ-pre.z, 
	 pre_shX, pre_shY, pre_shZ, +1);
  
  double post_shX, post_shY, post_shZ;
  Rotate(post.x-fFirstInteractionX, post.y-fFirstInteractionY, fFirstInteractionZ-post.z, 
	 post_shX, post_shY, post_shZ, +1);
  
  pre_shZ  += fFirstInteractionDist;
  post_shZ += fFirstInteractionDist;
  
  /*   ------------------------
       direction of particle in shower frame
  */
  const double dX = post_shX - pre_shX; 
  const double dY = post_shY - pre_shY; 
  const double dZ = post_shZ - pre_shZ; 
  const double length = sqrt (dX*dX + dY*dY + dZ*dZ);
    
  if (length==0) {
    /*   ------------------------ 
	 This happens for particles that are NOT tracked by Corsika
	 e.g. pi0 -> decay immediatly
	 SKIP them
    */
    return;
  }	
  
  const double theta   = acos(dZ/length) * rad;
  const double norm2Da = sqrt(pre_shX*pre_shX + pre_shY*pre_shY);
  const double norm2Db = sqrt(dX*dX + dY*dY);
  double phi;
  if (norm2Da==0 && norm2Db!=0) {
    phi = 0.*deg; // this particles starts at the axis and moves outward
  } else if (norm2Db==0) {
    phi = -1.*deg; // this particles moves parallel to the axis
  } else {
    double arg = ((pre_shX*dX + pre_shY*dY) / (norm2Da * norm2Db));
    if (arg > 1.0) { arg = 1.0; } else if (arg < -1.0) { arg = -1.0; } // necessary to compensate for rounding errors
    phi = acos(arg) * rad;
  }
  
  
  
  // check for crossing of planes
  for (int iLayer=0; iLayer<fNhist; ++iLayer) {
	
    const double zPlane = fZPlane[iLayer];
	
    if (fVerbosityLevel>=12) {
      cout << " COAST-Histogram> plane=" << zPlane << "cm, pre=" << pre_shZ << "cm"
	   << " post=" << post_shZ << "cm";
    }
    
    if (!((pre_shZ<zPlane && post_shZ>zPlane) || 
	  (post_shZ<zPlane && pre_shZ>zPlane))) {
      if (fVerbosityLevel>=12) {
	cout << " no hit " << endl;
      }
      continue;
    }
    
    if (fVerbosityLevel>=12) {
      cout << " intersect plane=" << iLayer << endl;
    }
    
    
    const double delta = (pre_shZ-zPlane) / (pre_shZ-post_shZ);	   
    
    const double e_post = post.energy - gParticleMass[(int)post.particleId-1];
    const double e_pre = pre.energy - gParticleMass[(int)pre.particleId-1];
    
    const double dT = post.time - pre.time;
    const double dE = e_post - e_pre;
    const double dW = post.weight - pre.weight;
    
    const double energy = e_pre + dE * delta;
    const double weight = pre.weight;
    const double x = pre_shX + dX * delta;
    const double y = pre_shY + dY * delta;
    // const double z = pre_shZ + dZ * delta; 

    double time = (pre.time + dT * delta); // * (s/ns);
    
    if (fVerbosityLevel>=10) {
      cout << " COAST-Histogram> PLANE> "
           << " e_pre=" << e_pre
           << " e_post=" << e_post
           << " dT=" << dT
           << " dE=" << dE
           << " dW=" << dW
           << " delta=" << delta
           << "\n";    
    }
    
    
    const double r = sqrt (pow (x, 2) + pow (y, 2));
    // double l = sqrt (pow (r, 2) + pow (z, 2));
	
    
    // times relative to what
    time -= fFirstInteractionTime;
    time *= s/ns; // convert from seconds to nanoseconds
    time -= ((zPlane-fFirstInteractionDist)/crs::cSpeedOfLight);  // PLANE SHOWER FRONT
    
    if (fCountInvalid>=fMaxInvalids) {
	cout << " COAST-Histogram> too many invalid histogram entries. Aborting simulation. Check your log file." 
	     << endl;
	exit(1);
    }
    
    if (time<=-fTimeErrorMargin) {
      cout << " COAST-Histogram> invalid TIME: " 
           << time/ns << " ns (filling in underflow-bin)" << endl;
      cout << " Track-pre>  "; pre.Dump();
      cout << " Track-post> "; post.Dump();
      ++fCountInvalid;
    }
    
    if (energy<=0) {
      cout << " COAST-Histogram>  invalid ENERGY: " 
           << energy << " GeV (filling in underflow-bin)" << endl;
      cout << " Track-pre>  "; pre.Dump();
      cout << " Track-post> "; post.Dump();
      ++fCountErrors;
      return;
    }
    
    if (r<=0) {
      cout << " COAST-Histogram>  invalid RADIUS: " 
           << r << " cm (filling in underflow-bin)" << endl;
      cout << " Track-pre>  "; pre.Dump();
      cout << " Track-post> "; post.Dump();
      ++fCountErrors;
      return;
    }
    
    
    const double height = fCORSIKA->GetVerticalHeightOfDistance(zPlane);
    const double T = fCORSIKA->Temperature(height);
    const double depth = fCORSIKA->GetVerticalDepthOfHeight(height);
    const double rm = fCORSIKA->MoliereRadius(depth, T);
    
    const int particleId = (int)pre.particleId;
    if (fParticles.count(particleId)) {
      
      // calculate logarithm, and put to underflow bin if not possible
      const double timeValue = (time>0 ? time : pow(10., fTimeUnderflow));
      
      FillHistograms(fHists[iLayer][particleId], 
		     x, y, timeValue, r, rm, energy, theta, phi, weight);
      
    }
  } 
}






void
TPlotter::SetFirstInteraction(double x, double y, double z, double X, double t) 
{
  fFirstInteractionX = x;
  fFirstInteractionY = y;
  fFirstInteractionZ = z;
  fFirstInteractionTime = t;
  fFirstInteractionSlantDepth = X;
  fFirstInteractionDist = fCORSIKA->GetDistanceOfSlantDepth(fFirstInteractionSlantDepth);
  //  fFirstInteractionDist = fCORSIKA->GetDistanceOfHeight(fFirstInteractionZ);

  cout << " COAST-Histogram> fFirstInteractionDist=" << fFirstInteractionDist/m << "m "
       << ", fFirstInteractionZ=" << fFirstInteractionZ/m << "m "
       << ", fFirstInteractionSlantDepth=" << fFirstInteractionSlantDepth/g*cm*cm << "g/cm2 "
       << ", fFirstInteractionTime=" << fFirstInteractionTime << "s "
       << endl;
}


void 
TPlotter::SetShowerAxis(double zenith, double azimuth) 
{  
  cout << " COAST-Histogram> zenith=" << zenith/deg << "deg"	
       << "   azimuth=" << azimuth/deg << "deg"	
       << endl;
  fZenith = zenith;
  fAzimuth = azimuth;
  fCosZenith = cos(zenith);
  fSinZenith = sin(zenith);
  fCosAzimuth = cos(azimuth);
  fSinAzimuth = sin(azimuth);
  fAxisZ = fCosZenith;
  fAxisX = fSinZenith * fCosAzimuth;;
  fAxisZ = fSinZenith * fSinAzimuth;
}

void 
TPlotter::AddInteraction(const crs::CInteraction& interaction) 
{
  if (fVerbosityLevel >= 10) {
    cout << " COAST: ";
    interaction.Dump();
  }
  
  if (fPrimaryTrack) 
    fFirstInteraction = true;
  else 
    fFirstInteraction = false;
  
  fPrimaryTrack = false;
}

