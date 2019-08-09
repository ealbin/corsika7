/****************************************************************************
 *                                                                          *
 *  Copyright and any other appropriate legal protection of these           *
 *  computer programs and associated documentation are reserved in          *
 *  all countries of the world.                                             *
 *                                                                          *
 *  These programs or documentation may not be reproduced by any method     *
 *  without prior written consent of Karlsruhe Institute of Technology (KIT)*
 *  ot its delegate. Commercial and military use are explicitly forbidden.  *
 *                                                                          *
 *  The Karlsruhe Institute of Technology welcomes comments concerning the  *
 *  COAST code but undertakes no obligation for maintenance of the programs,*
 *  nor responsibility for their correctness, and accepts no liability      *
 *  whatsoever resulting from the use of its programs.                      *
 *                                                                          *
 ****************************************************************************/

/*
  TP, Tue August  16 18:04:41 EDT 2010
  Fix problems related to CONEX
  - Correct first interaction distance taking into account the height
    of the observation level with CURVED option
  - use absolute value of ParticleId  

  RU, Mon May  3 11:04:41 EDT 2010
  Fixed problems related to STACKIN + PRESHOWER  

  RU, Wed Dec 16 17:38:59 EST 2009
  Fixed several problems concering CURVED + SKIMMING
  - Significant development of TCorsika
  - Avoiding the use of tabulated atmosphere as far as possible
  - Fixed the issue of rare empty histograms (caused by first 
    interaction point extremely high in the atmosphere <5g/cm2)
  - Improved the output of histprofile

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



/* -------------------------------------------------------------------
   Define the particle names/types and CORSIKA ID here !!!!!!!!!!!!!!!!
*/
const unsigned int nParticles = 4;
//string particle_name [3] = {"EM", "Muon", "Hadron"};
const string gParticle_name [nParticles] = {"electron", "positron"};//,"mu-","mu+"};
const int gParticle_id      [nParticles] = {3, 2};//, 6, 5};
const int gParticle_color   [nParticles] = {kBlue, kRed};//, 6, 5};
// -------------------------------------------------------------------



/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   - first find the point of the first interaction, and define the range
     where the shower should get histogramed

   - rotate all particles in shower coordinate system, with (0,0,0) at 
     the first interaction

   - define planes in equal slant depth bins 

   - fill histograms    
 */

TPlotter::TPlotter() {
  
  // do all lateral histogramming with respect to this scale radius
  fLateralScaleRadius = 100.*m;

  //fVerbosityLevel = 100; // all includning TCorsika
  //fVerbosityLevel = 40; // all from TPlotter
  fVerbosityLevel = 0; // none
  
  fMaxInvalids = 10000;
  fTimeErrorMargin = 0.0;
  
  fPrimaryTrack = true;
  fFirstInteraction = false;
  
  fSlant = false;
  fCurved = false;
  fStackInput = false;
  fPreshower = false;
  
  fTimeUnderflow   = 0;
  fRadiusUnderflow = 0;
  fEnergyUnderflow = 0;
  fCountInvalid  = 0;
  
//itssavedemweight = 0;
//itssavedhdweight = 0;
  
  fCORSIKA = 0;
  
  Clear ();
  fDirectory = "";
    
  for (unsigned int j=0; j<nParticles; j++) {
    fParticles[gParticle_id[j]].name = gParticle_name[j];
    fParticles[gParticle_id[j]].color = gParticle_color[j];
  }
}



TPlotter::~TPlotter () {
  
  fParticles.clear();
  Clear();
  delete fFile;
}


void 
TPlotter::Welcome() 
	 const 
{
  cout << "\n"
       << "COAST: *************************************************\n"
       << "COAST: **                                               \n"
       << "COAST: ** You are using COAST-histogramming THRadio      \n"
       << "COAST: **    COAST version: " << VERSION << " \n"
       << "COAST: **  THRadio version: " << HISTOGRAMING_VERSION << " \n"
       << "COAST: **   Histograms for: REAS3 	                   \n"
       << "COAST: **                                               \n"
       << "COAST: **  THIN/CURVED/SLANT/STACKIN/PRESHOWER: "
       << fThinning << "/" << fCurved << "/" << fSlant << "/" << fStackInput << "/" << fPreshower << "\n"
       << "COAST: **                                               \n"
       << "COAST: *************************************************\n\n"
       << endl;
}




void 
TPlotter::Clear () 
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
  
  fHists.clear();
  fLayerDepth.clear();
  fLayerHeight.clear();
  fLayerDist.clear();
  fNhist = 0;
  fCORSIKA = 0;
}


void TPlotter::SetRunHeader(const crs::MRunHeader &header) 
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
    cout << " COAST: Detected skimming geometry with impact=" 
         << fSkimmingAltitude/km
         << " km" << endl;
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
      cout << " COAST: height of first interaction is NEGATIVE, corrected." 
           << endl;
    }
    fHeightFirstInt = fabs(fHeightFirstInt);
  }
  
  cout  << " COAST: firstint " << fHeightFirstInt/m << "m\n"
        << " COAST: obslev " << fObservationLevel/m 
        << "m, including safety offset of 10cm" 
        << endl;
}


void 
TPlotter::SetShowerTrailer(const crs::MEventEnd &trailer) 
{
  fXmax = trailer.GetXmax()*g/cm/cm;
  if (!fSlant && fCosZenith!=0) 
    fXmax /= fCosZenith;
  
  cout << " COAST: --------- shower end ------ \n"
       << " COAST: xmax: " << fXmax/g*cm*cm << "\n";
  if (fCountInvalid>0) 
    cout << " COAST: invalid-counter: " << fCountInvalid ;
  cout << endl;
}


void 
TPlotter::Init(const int nBins, const double maxRange) 
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
  fLayerDist.clear();

  if( fCORSIKA == 0 )
    cout << "Warning fCORSIKA not defined" << endl;
  else
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
  
  if (fStackInput || fPreshower) { // STACKIN OPTION NEEDS SPECIAL TREATMENT !!!!
    // define safety margin to prevent invalid times (in front of shower front)
    // caused by rounding problems
    // #warning is this error margin still required?
    fTimeErrorMargin = 3.3e-3*ns;	// ignore time errors corresponding to ~1 mm before shower front
    
    fPrimaryTrack = false;
    fFirstInteraction = true;
  }
  
  if (fVerbosityLevel > 1 && fCurved) {
    fCORSIKA->DumpAtmosphereTable();
  }
}


void 
TPlotter::CalculatePlanes() 
{
  if (fVerbosityLevel >= 10) 
    cout << "TPlotter::CalculatePlanes" << endl; 
  
  // tell TCorsika the real starting point for histogramm layers
  //fCORSIKA->SetHeightOfFirstInteraction(fFirstInteractionZ);
  //fCORSIKA->SetHeightOfFirstInteraction(fCORSIKA->GetHeightOfSlantDepth(fFirstInteractionSlantDepth));
  
  const double depthStart = fFirstInteractionSlantDepth; 
  const double depthEnd = (fSkimming ? 
			   depthStart + fMaxShowerRange*g/cm/cm :
			   min(depthStart + fMaxShowerRange*g/cm/cm, 
			       fCORSIKA->GetLastLevelSlantDepth()));
  
  const double depthUnit = g/cm/cm;
  const string depthStr = "slant depth = ";
  const string unitStr = "g/cm^2";
  
  cout << " COAST: start shower=" << depthStart/depthUnit << " " << unitStr
       << ", last histogram plane=" << depthEnd/depthUnit << " " << unitStr
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
  const double deltaSlice = (depthEnd - depthStart) / nBins;

  fdeltaSlice = deltaSlice / depthUnit ;
  fdepthStart = depthStart / depthUnit;
  
  // --------------------------------------------------------------
  // create histogramms for each layer for this shower
  for (int i=0; i<nBins; ++i) {
    
    const double depthPlane = depthStart + deltaSlice * (i + 1);
    
    /* ----------------
       height of planes in shower frame
    */
    const double distPlane = fCORSIKA->GetDistanceOfSlantDepth(depthPlane); 
    const double heightPlane = fCORSIKA->GetVerticalHeightOfDistance(distPlane);
    
    fLayerDepth.push_back(fCORSIKA->GetSlantDepthOfDistance(distPlane));
    fLayerHeight.push_back(heightPlane);
    fLayerDist.push_back(distPlane);
    
    if (fVerbosityLevel >= 2) {
      cout << " COAST: plane " 
	   << " depthPlane=" << depthPlane/depthUnit << unitStr
	   << ", dist=" << distPlane/km << "km"
           << ", h=" << heightPlane/km << "km"
	   << endl;
    }
    
    ostringstream h_title;
    h_title << depthStr
	    << depthPlane/depthUnit
	    << " " << unitStr;
    
    
    HistLayer hists; // this contains all histogramms for one layer

    // ------------------------------------------------------------------
    // loop over all particle types
    for (map<int,ParticleDef>::const_iterator iParticle = fParticles.begin();
	 iParticle != fParticles.end();
	 ++iParticle) {
      
      const string name = iParticle->second.name;
      const int color = iParticle->second.color;
      const int id = iParticle->first;
      
      ostringstream hT_name, hA_name, hE_name;
      hT_name << "hT_" << i << "_" << name;
      hA_name << "hA_" << i << "_" << name;
      hE_name << "hE_" << i << "_" << name;
      
      
      // --------------------------------------------------
      // generate histograms

      const int nBinsTime  = 60;
      const double minLogTime = -3.0; // log10 (time/ns)
      const double maxLogTime = 4.5;  // log10 (time/ns)
      
      const int nBinsR  = 50;
      const double minLogR = -3.; // log10 (r/100m) 
      const double maxLogR = 2.;  // log10 (r/100m) 
      
      const int nBinsE  = 40;
      const double minLogE = -4.; // log10 (e/GeV)
      const double maxLogE = 2.;  // log10 (e/GeV)
      
      const int nBinsZenith  = 50;
      const double minZenith = 0;     // [deg]
      const double maxZenith = 90;    // [deg]
      
      const int    nBinsAzimuth  = 72;   //36 in original
      const double minAzimuth = 0;       // [deg]
      const double maxAzimuth = 360;     // [deg] // 180 in original
      
      TH3D *histT = new TH3D (hT_name.str ().c_str (),
			      h_title.str ().c_str (),
			      nBinsTime, minLogTime, maxLogTime, 
			      nBinsR, minLogR, maxLogR,   
			      nBinsE, minLogE, maxLogE);  
      
      TH3D *histA = new TH3D (hA_name.str ().c_str (),
			      h_title.str ().c_str (),
			      nBinsZenith, minZenith, maxZenith, 
			      nBinsAzimuth, minAzimuth, maxAzimuth,
                              nBinsE, minLogE, maxLogE);  
      
      TH1D *histE = 0;
      //hE TH1D *histE = new TH1D (hE_name.str ().c_str (),
      //hE		      h_title.str ().c_str (),
      //hE		      200, -3.4, -2.7); // log10 (e/GeV)
      
      // --------------------------------------------------
      // for particles giving un-physical values because of
      // numerical problems, use these values to fill into
      // the underflow bin of the histogram
      fTimeUnderflow   = histT->GetXaxis()->GetXmin() - 1.0;
      fRadiusUnderflow = histT->GetYaxis()->GetXmin() - 1.0;
      fEnergyUnderflow = std::min(histT->GetZaxis()->GetXmin() - 1.0,
                                  histA->GetZaxis()->GetXmin() - 1.0);
      
      // --------------------------------------------------
      // set histogram properties
      histT->SetLineColor(color);
      histT->SetXTitle("log_{10}(time/ns)");
      histT->SetYTitle("log_{10}(r/100m)");
      histT->SetZTitle("log_{10}(E_{kin}/GeV)");

      histA->SetLineColor(color);
      histA->SetXTitle("theta/deg");
      //histA->SetYTitle ("log_{10}(r/100m)");
      histA->SetYTitle("phi/deg");
      histA->SetZTitle("log_{10}(E_{kin}/GeV)");
      
      //hE histE->SetLineColor (particle_id);
      //hE histE->SetXTitle ("log_{10}(E_{kin}/GeV)");
      
      // -------------------------------------------------
      // create data object for one layer
      hists[id] = Hists(histT, histA, histE);
      
      for (int iEnergyCenter = 0; iEnergyCenter<=nBinsE; ++iEnergyCenter) {
        Center center;
        center.x = 0;
        center.y = 0;
        center.w = 0;
        center.x2= 0;
        center.y2= 0;
        center.w2= 0;
        hists[id].center[iEnergyCenter] = center;
      }
    } // end loop particle types
    
    // -------------------------------------------------
    // collect all layers
    fHists.push_back(hists); 
    
  } // end loop layers
  
}





void 
TPlotter::Write () 
{
  if (fVerbosityLevel >= 2) {
    cout << " COAST: Write " << endl;
    /*
    ostringstream gif_dir;
    gif_dir << fFileName 
    << "_"
	    << fEventNo;
    */
  }

  fXmax /= g*cm*cm;

  // create the particle-string
  ostringstream particles_ss;
  for (map<int,ParticleDef>::iterator iParticle = fParticles.begin ();
       iParticle != fParticles.end ();
       ++iParticle) {
    particles_ss << iParticle->second.name << " ";
  }
  const string particles_str(particles_ss.str());
  
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
  char *particles = const_cast<char*>(particles_str.c_str());
  run->Branch ("particles" , particles      , "particles/C");
  char *versionCorsika = const_cast<char*>(fCorsikaVersion.c_str());
  char *versionCoast   = const_cast<char*>(fCoastVersion.c_str());
  char *versionHist    = const_cast<char*>(fHistVersion.c_str());
  run->Branch ("corsika_version" , versionCorsika  , "corsika_version/C");
  run->Branch ("coast_version"   , versionCoast    , "coast_version/C");
  run->Branch ("thradio_version" , versionHist     , "thradio_version/C");
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
  for (map<int,ParticleDef>::iterator iParticle = fParticles.begin ();
       iParticle != fParticles.end ();
       ++iParticle) {
	
    int id = iParticle->first;
    string name = iParticle->second.name;

    ostringstream tree_name, tree_info;
    tree_name << "data_" << name;
    tree_info << "All the histograms for " << name;
    
    TTree *newTree = new TTree (tree_name.str ().c_str (), 
				tree_info.str ().c_str ());
    
    trees[id] = newTree;
  }

  bool firstFill = true;
  for (int i=0; i<fNhist; ++i) {

    for (map<int,ParticleDef>::iterator iParticle = fParticles.begin ();
	 iParticle != fParticles.end ();
	 ++iParticle) {
	
      int id = iParticle->first;
      string name = iParticle->second.name;

      TH3D *hT = 0;
      TH3D *hA = 0;
      //hE TH1D *hE = 0;
      
      hT = fHists[i][id].Time;
      hA = fHists[i][id].Angle;
      //hE hE = fHists [i][iParticle->second].Energy;
      
      vector<double> vCenter_E, vCenter_E_err;
      vector<double> vCenter_x, vCenter_x_err;
      vector<double> vCenter_y, vCenter_y_err;
      for (map<int,Center>::const_iterator iCenterBin 
             = fHists[i][id].center.begin();
           iCenterBin != fHists[i][id].center.end();
           ++iCenterBin) {
        
        int centerEnergyBin = iCenterBin->first;
        const Center& center = iCenterBin->second;
        
        double binEnergy = hT->GetZaxis()->GetBinCenter(centerEnergyBin);
        double binEnergyErr = 0;
        double center_x = 0;
        double center_y = 0;
        double center_x_err = 0;
        double center_y_err = 0;
        
        if (center.w>0) {
          
          center_x = center.x / center.w;
          center_y = center.y / center.w;
          
          double neff = pow(center.w, 2) / center.w2;
          
          double x_var = center.x2 / center.w - center_x*center_x;
	  if (x_var < 0.0) { x_var = 0.0; }	// correct rounding errors
          double x_rms = sqrt(x_var);
          center_x_err = x_rms / sqrt(neff);
          
          double y_var = center.y2 / center.w - center_y*center_y;
	  if (y_var < 0.0) { y_var = 0.0; }	// correct rounding errors
          double y_rms = sqrt(y_var);
          center_y_err = y_rms / sqrt(neff);

          vCenter_E.push_back(binEnergy);
          vCenter_E_err.push_back(binEnergyErr);
          vCenter_x.push_back(center_x);
          vCenter_y.push_back(center_y);
          vCenter_x_err.push_back(center_x_err);
          vCenter_y_err.push_back(center_y_err);
        }
      } // loop energy bins
      
      ostringstream bnameT, bnameA, bnameOx, bnameOy;
      
      bnameT << string("hT"+ iParticle->second.name).c_str();
      bnameA << string("hA"+ iParticle->second.name).c_str();
      bnameOx << "offsetX";
      bnameOy << "offsetY";
      
      TGraphErrors* hOx 
        = new TGraphErrors(vCenter_E.size(),
                           &vCenter_E.front(), &vCenter_x.front(),
                           &vCenter_E_err.front(), &vCenter_x_err.front());
      TGraphErrors* hOy 
        = new TGraphErrors(vCenter_E.size(),
                           &vCenter_E.front(), &vCenter_y.front(),
                           &vCenter_E_err.front(), &vCenter_y_err.front());
      
      
      if (firstFill) {
        
	trees[id]->Branch(bnameT.str().c_str(),
                          "TH3D", &hT, 32000, 0);
	trees[id]->Branch(bnameA.str().c_str(),
                          "TH3D", &hA, 32000, 0);
	trees[id]->Branch(bnameOx.str().c_str(),
                          "TGraphErrors", &hOx, 32000, 0);
	trees[id]->Branch(bnameOy.str().c_str(),
                          "TGraphErrors", &hOy, 32000, 0);

      } else {

	trees[id]->SetBranchAddress(bnameT.str().c_str(), &hT);
	trees[id]->SetBranchAddress(bnameA.str().c_str(), &hA);
	trees[id]->SetBranchAddress(bnameOx.str().c_str(), &hOx);
	trees[id]->SetBranchAddress(bnameOy.str().c_str(), &hOy);
      }
      
      time = (fLayerDist[i] - fLayerDist[0]) / crs::cSpeedOfLight;
      dist = fLayerDist[i];  
      slantDepth = fLayerDepth[i]; 
      height = fLayerHeight[i];
      
      double T = fCORSIKA->Temperature(height);	
      //rm = fCORSIKA->MoliereRadius(slantDepth*fCosZenith, T)/m;
      rm = fLateralScaleRadius/m; // for consistency with histogramming, has to be written into the root branch in the units of meter!
      
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
      
      trees[id]->Fill ();
      
      //hT->Print ();
      //hA->Print ();
    }
	
    tree->Fill ();
    firstFill = false;
	
    // delete hists
    for (map<int,ParticleDef>::iterator iParticle = fParticles.begin ();
	 iParticle != fParticles.end ();
	 ++iParticle) {
      int id = iParticle->first;
      delete (fHists [i][id].Time);
      delete (fHists [i][id].Angle);
      delete (fHists [i][id].Energy);

    }

  }
    
  run->Write ("", TFile::kOverwrite);
  tree->Write ("", TFile::kOverwrite);
  for (map<int,TTree*>::iterator iTree = trees.begin();
       iTree!=trees.end(); ++iTree) {
    iTree->second->Write("", TFile::kOverwrite);
  }
  fFile->Close ();
    
  // clean up
  Clear ();
}




// #define Rotate(x, y, sx, sy, cosa, sina) { 
inline
void 
TPlotter::Rotate(double x, double y, double z,
		 double &sx, double &sy, double &sz,
		 int inverse) const
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
TPlotter::AddTrack(const crs::CParticle &pre, 
		   const crs::CParticle &post) 
{  
  if (fVerbosityLevel >= 9) {
    /*  double dX = pre.x-fFirstInteractionX; 
  double dY = pre.y-fFirstInteractionY; 
  double dZ = pre.z-fFirstInteractionZ; 
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);*/
    cout << " COAST:   PRE> "
	 << " id=" << (int)pre.particleId
	 << " E=" << pre.energy
	 << " w=" << pre.weight
	 << " x=" << pre.x/km << "km"
	 << " y=" << pre.y/km << "km"
	 << " z=" << pre.z/km << "km"
	 << " t=" << pre.time*s/ns << "ns"
	 << " X=" << pre.depth << "g/cm^2"
      //         << " v= " <<  length/(pre.time*s/ns-fFirstInteractionTime*s/ns)/crs::cSpeedOfLight 
	 << endl;
  }
  
  if (fVerbosityLevel >= 10) {
    /*  double dX = post.x-fFirstInteractionX; 
  double dY = post.y-fFirstInteractionY; 
  double dZ = post.z-fFirstInteractionZ; 
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);*/
    cout << " COAST:  POST> "
	 << " id=" << (int)post.particleId
	 << " E=" << post.energy
	 << " w=" << post.weight
	 << " x=" << post.x/km << "km"
	 << " y=" << post.y/km << "km"
	 << " z=" << post.z/km << "km"
	 << " t=" << post.time*s/ns << "ns"
	 << " X=" << post.depth << "g/cm^2"
      //        << " v= " <<  length/(post.time*s/ns-fFirstInteractionTime*s/ns)/crs::cSpeedOfLight
	 << endl;
  }
    
  /*
    Skip the track of the primary particle, which is in a different 
    reference system as the shower !!!
  */
  //static double cosZenith_track = fCosZenith;
  //static double azimuth_track = fAzimuth;
  if (fPrimaryTrack) {
    if (fVerbosityLevel >= 2) {
      cout << " COAST: Primary track  " << endl;
    }    
    /*
      const double dX_track = post.x - pre.x; 
      const double dY_track = post.y - pre.y; 
      const double dZ_track = post.z - pre.z; 
      const double length_track = sqrt(pow(dX_track, 2)+ 
      pow(dY_track, 2)+ 
      pow(dZ_track, 2));    
      if (fabs(length_track) < 1e-7) {
      return; // pi0 could be the first secondary 
      }
    */
    /*
    // Unfortunately this is too unstable ... so we have to protect it
    const double new_cosZenith_track = -dZ_track/length_track;
    const double new_azimuth_track = atan2(dY_track, dX_track);
    if (std::fabs(new_cosZenith_track-cosZenith_track) > 0.1*deg ||
    std::fabs(new_azimuth_track-azimuth_track) > 0.1*deg) {
    if (fVerbosityLevel>=5) {
    cout << " COAST-debug: primary tracksegment: too much fluctuating! "
    << " cosZenith=" << new_cosZenith_track
    << " zenith=" << acos(new_cosZenith_track)/deg << " deg"
    << " azimuth=" << new_azimuth_track/deg << " deg"	     
    << endl;
    }      
    } else {
    cosZenith_track = -dZ_track/length_track;
    azimuth_track = atan2(dY_track, dX_track);
    if (fVerbosityLevel>=5) {
    cout << " COAST-debug: primary tracksegment: cosZenith=" << cosZenith_track
    << " zenith=" << acos(cosZenith_track)/deg << " deg"
    << " azimuth=" << azimuth_track/deg << " deg"
    <<  endl;
    }
    } 
    */
    return;
  }
  
  if (fFirstInteraction) {
    
    fFirstInteraction = false;
    
    /*
      // this was not leading to improvements ...
      if (fCurved) { 
      SetShowerAxis(acos(cosZenith_track), azimuth_track); 
      } 
    */   
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
  Rotate(pre.x - fFirstInteractionX, pre.y - fFirstInteractionY, fFirstInteractionZ - pre.z, 
	 pre_shX, pre_shY, pre_shZ, +1);
  
  double post_shX, post_shY, post_shZ;
  Rotate(post.x - fFirstInteractionX, post.y - fFirstInteractionY, fFirstInteractionZ - post.z, 
	 post_shX, post_shY, post_shZ, +1);
  
  pre_shZ  += fFirstInteractionDist;
  post_shZ += fFirstInteractionDist;
  
  /*   ------------------------
       direction of particle in shower frame
  */
  double dX = post_shX - pre_shX; 
  double dY = post_shY - pre_shY; 
  double dZ = post_shZ - pre_shZ; 
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);
  //  cout << "speed " <<  length/(post.time*s/ns-pre.time*s/ns)/crs::cSpeedOfLight << endl;

  if (length==0) {
    /*   ------------------------ 
	 This happens for particles that are NOT tracked by Corsika
	 e.g. pi0 -> decay immediatly
	 SKIP them
    */
    return;
  }	
  
  double theta   = acos(dZ/length) * rad;
//double norm2Da = sqrt(pre_shX*pre_shX + pre_shY*pre_shY);
  double norm2Db = sqrt(dX*dX + dY*dY);
  double phi;

/* REAS2 histograms
  if (norm2Da==0 && norm2Db!=0) {
    phi = 0.*deg; // this particles starts at the axis and moves outward
  } else if (norm2Db==0) {
    phi = -1.*deg; // this particles moves parallel to the axis
  } else {
    double arg = ((pre_shX*dX + pre_shY*dY) / (norm2Da * norm2Db));
    if (arg > 1.0) { arg = 1.0; } else if (arg < -1.0) { arg = -1.0; } // necessary to compensate for rounding errors
    phi = acos(arg) * rad;
  } */

  // REAS3 histograms
  if (norm2Db == 0.0) 
  {
    phi = -1.*deg; // this particles moves parallel to the axis
  } 
  else 
  {
    phi = atan2(dY, dX) * rad;
    if (phi < 0.0)
    { phi += 360.*deg; }
  }
  
  int iStartLayer = max(0,int((pre.depth - fdepthStart) / fdeltaSlice ) - 1 );
  int iEndLayer = min(fNhist,int((post.depth - fdepthStart) / fdeltaSlice) + 1);
  // check for crossing of planes
    //  for (int iLayer=0; iLayer<fNhist; ++iLayer) {
  for (int iLayer=iStartLayer; iLayer<iEndLayer; ++iLayer) {
	
    const double distPlane = fLayerDist[iLayer];
    
    if (fVerbosityLevel >= 12) {
      cout << " COAST: plane=" << distPlane/km << "km, pre=" << pre_shZ/km << "km"
	   << " post=" << post_shZ/km << "km " << iStartLayer << " " << iLayer << " " << iEndLayer - 2 ;
    }
    
    if (!((pre_shZ < distPlane && post_shZ > distPlane) || 
	  (post_shZ < distPlane && pre_shZ > distPlane))) {
      if (fVerbosityLevel >= 12) 
	cout << "  ----> no hit " << endl;
      continue;
    }
    
    if (fVerbosityLevel >= 12) {
      cout << " intersect plane=" << iLayer << endl;
    }
    
    
    double delta = (pre_shZ-distPlane) / (pre_shZ-post_shZ);	   
    
    double e_post = post.energy - gParticleMass[abs((int)post.particleId)-1];
    double e_pre = pre.energy - gParticleMass[abs((int)pre.particleId)-1];
    
    double dT = post.time - pre.time;
    double dE = e_post - e_pre;
    // double dW = post.weight - pre.weight;

    
    double time = (pre.time + dT * delta); // * (s/ns);
    double energy = e_pre + dE * delta;
    double weight = pre.weight;
    double x = pre_shX + dX * delta;
    double y = pre_shY + dY * delta;
    //double z = pre_shZ + dZ * delta; 

    if (fVerbosityLevel >= 10) {
      cout << " COAST: PLANE> "
	   << " planeId: " << iLayer
	   << " pId: " << pre.particleId
           << " e_pre=" << e_pre
           << " e_post=" << e_post
           << " dT=" << dT
           << " dE=" << dE
           << " delta=" << delta
           << "\n";    
    }
    
    
    double r = sqrt(pow (x, 2) + pow (y, 2));
    // double l = sqrt (pow (r, 2) + pow (z, 2));
    
    
    // times relative to what
    ///// //time -= l / crs::cSpeedOfLight;  // CURVED SHOWER FRONT
    time -= fFirstInteractionTime;
    time *= s/ns; // convert from seconds to nanoseconds
    time -= ((distPlane - fFirstInteractionDist)/crs::cSpeedOfLight);  // PLANE SHOWER FRONT
    // time -= (distPlane/crs::cSpeedOfLight);  // PLANE SHOWER FRONT
    
    if (fCountInvalid >= fMaxInvalids) {
	cout << " COAST: too many invalid histogram entries. Aborting simulation. Check your log file." 
	     << endl;
	exit(1);
    }

    if (time <= -fTimeErrorMargin) {
      cout << " COAST: invalid TIME: " 
           << time/ns << " ns (filling in underflow-bin)" << endl;
      cout << " Track-pre>  "; pre.Dump();
      cout << " Track-post> "; post.Dump();
      fCountInvalid++;
      //	exit(1);
    }
    
    if (energy <= 0) {
      cout << " invalid ENERGY: " 
           << energy << " GeV (filling in underflow-bin)" << endl;
      cout << " Track-pre>  "; pre.Dump();
      cout << " Track-post> "; post.Dump();
      fCountInvalid++;
    }
    
    if (r <= 0) {
      cout << " invalid RADIUS: " 
           << r/km << " km (filling in underflow-bin)" << endl;
      cout << " Track-pre>  "; pre.Dump();
      cout << " Track-post> "; post.Dump();
      fCountInvalid++;
    }
    
    
    /*cout << " time/ns = " << time/ns << " offset=" << ((distPlane-fFirstInteractionDist)/crs::cSpeedOfLight)/ns
	 << " t_1=" << fFirstInteractionTime * s/ns
	 << "\n";
    */
    
    if (fParticles.count(particleId)) {
      
      string particle_name = fParticles[particleId].name;
      
      // calculate logarithm, and put to underflow bin if not possible
      double timeValue   = (time>0   ? log10(time/ns) : fTimeUnderflow);
      double radiusValue = (r>0      ? log10(r/fLateralScaleRadius)  : fRadiusUnderflow);
      double energyValue = (energy>0 ? log10(energy)  : fEnergyUnderflow);

      int energyBin 
        = fHists [iLayer][particleId].Time->GetZaxis()->FindBin(energyValue);
      
      fHists[iLayer][particleId].Time
	->Fill(timeValue, radiusValue, energyValue, weight);
      
      fHists[iLayer][particleId].Angle
	->Fill(theta/deg, phi/deg, energyValue, weight);
      
      /*
	fHists [iLayer][particle_name].Energy
	->Fill (energy,
	weight);
      */
            
      fHists[iLayer][particleId].center[energyBin].x += x * weight;
      fHists[iLayer][particleId].center[energyBin].y += y * weight;
      fHists[iLayer][particleId].center[energyBin].w += weight;
      fHists[iLayer][particleId].center[energyBin].x2+= x*x * weight;
      fHists[iLayer][particleId].center[energyBin].y2+= y*y * weight;
      fHists[iLayer][particleId].center[energyBin].w2+= weight*weight;
      
    }
  } 
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




TH1*
TPlotter::Interpolate3D (TH1* h1, double age1,
			 TH1* h2, double age2,
			 double age) 
  const
{
  if (!((age<age1 && age>age2) ||
	(age<age2 && age>age1))) {
	
    cerr << " CANNOT interpolate because age is outside [age1,age2]." 
	 << " age=" << age << " age1="<<age1 << " age2=" << age2 
	 << "\n";
    return new TH3D ();	
  }
  
  ostringstream hName, hTitle;
  hName << h1->GetName ()
	<< "_interpol";
  hTitle << h1->GetTitle ()
	 << ", age=" << age;
    
  TH1 *h = (TH1*) h2->Clone (hName.str ().c_str ());
  
  h->Reset ("ICE");
  
  h->Add (h1, h2, std::fabs(age-age1), std::fabs(age-age2));
  h->SetTitle (hTitle.str ().c_str());
  
  return h;
}



void 
TPlotter::SetFirstInteraction(const double x, const double y, const double z, 
			      const double X, const double t) 
{
  fFirstInteractionX = x;
  fFirstInteractionY = y;
  fFirstInteractionZ = z;
  fFirstInteractionTime = t;
  fFirstInteractionSlantDepth = X;
  if (X<0) {
    cout << " COAST: WARNING FirstInteractionDepth smaller than 0 (" << X/g*cm*cm << " g/cm2) " << endl;
    fFirstInteractionSlantDepth = 0;
  }
  //fFirstInteractionDist = fCORSIKA->GetDistanceOfSlantDepth(fFirstInteractionSlantDepth);
  fFirstInteractionDist = fCORSIKA->GetDistanceOfPoint(x, y, z); 
  
  cout << " COAST-SetFirstInteraction: Dist=" << fFirstInteractionDist/km << "km "
       << ", x=" << fFirstInteractionX/km << "km "
       << ", y=" << fFirstInteractionY/km << "km "
       << ", z=" << fFirstInteractionZ/km << "km "
       << ", SlantDepth=" << fFirstInteractionSlantDepth/g*cm*cm << "g/cm2 "
       << ", Time=" << fFirstInteractionTime << "s "
       << endl;
}


void 
TPlotter::SetShowerAxis(const double zenith, const double azimuth) 
{  
  cout << " COAST-SetShowerAxis: zenith=" << zenith/deg << "deg"	
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
  fAxisY = fSinZenith * fSinAzimuth;
}


