#include <TPlotter.h>

#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
using namespace crs;

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <string>
#include <iostream>
#include <sstream>
using namespace std;



TPlotter::TPlotter () {

    Clear ();
}



TPlotter::~TPlotter () {

    Clear ();
    delete fFile;
}


void TPlotter::Clear () {

    fPrimaryTrack = true;

    fXmax = 0;
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
    
    for (vector<TH3D*>::iterator i = fHistTGamma.begin ();
	 i!=fHistTGamma.end ();++i)
	delete *i;
    fHistTGamma.clear ();
    for (vector<TH3D*>::iterator i = fHistTElectron.begin ();
	 i!=fHistTElectron.end ();++i)
	delete *i;
    fHistTElectron.clear ();
    /*for (vector<TH3D*>::iterator i = fHistTPositron.begin ();
	 i!=fHistTPositron.end ();++i)
	delete *i;
	fHistTPositron.clear ();*/
    for (vector<TH3D*>::iterator i = fHistTMuon.begin ();
	 i!=fHistTMuon.end ();++i)
	delete *i;
    fHistTMuon.clear ();
    for (vector<TH3D*>::iterator i = fHistTHadron.begin ();
	 i!=fHistTHadron.end ();++i)
	delete *i;
    fHistTHadron.clear ();



    for (vector<TH2D*>::iterator i = fHistEGamma.begin ();
	 i!=fHistEGamma.end ();++i)
	delete *i;
    fHistEGamma.clear ();
    for (vector<TH2D*>::iterator i = fHistEElectron.begin ();
	 i!=fHistEElectron.end ();++i)
	delete *i;
    fHistEElectron.clear ();
    /*for (vector<TH2D*>::iterator i = fHistEPositron.begin ();
	 i!=fHistEPositron.end ();++i)
	delete *i;
	fHistEPositron.clear ();*/
    for (vector<TH2D*>::iterator i = fHistEMuon.begin ();
	 i!=fHistEMuon.end ();++i)
	delete *i;
    fHistEMuon.clear ();
    for (vector<TH2D*>::iterator i = fHistEHadron.begin ();
	 i!=fHistEHadron.end ();++i)
	delete *i;
    fHistEHadron.clear ();



    for (vector<TH2D*>::iterator i = fHistLGamma.begin ();
	 i!=fHistLGamma.end ();++i)
	delete *i;
	fHistLGamma.clear ();
    for (vector<TH2D*>::iterator i = fHistLElectron.begin ();
	 i!=fHistLElectron.end ();++i)
	delete *i;
    fHistLElectron.clear ();
    /*
    for (vector<TH2D*>::iterator i = fHistLPositron.begin ();
	 i!=fHistLPositron.end ();++i)
	delete *i;
	fHistLPositron.clear ();*/
    for (vector<TH2D*>::iterator i = fHistLMuon.begin ();
	 i!=fHistLMuon.end ();++i)
	delete *i;
    fHistLMuon.clear ();
    for (vector<TH2D*>::iterator i = fHistLHadron.begin ();
	 i!=fHistLHadron.end ();++i)
	delete *i;
    fHistLHadron.clear ();



    for (vector<TH3D*>::iterator i = fHistAGamma.begin ();
	 i!=fHistAGamma.end ();++i)
	delete *i;
	fHistAGamma.clear ();
    for (vector<TH3D*>::iterator i = fHistAElectron.begin ();
	 i!=fHistAElectron.end ();++i)
	delete *i;
    fHistAElectron.clear ();
    /*for (vector<TH3D*>::iterator i = fHistAPositron.begin ();
	 i!=fHistAPositron.end ();++i)
	delete *i;
	fHistAPositron.clear ();*/
    for (vector<TH3D*>::iterator i = fHistAMuon.begin ();
	 i!=fHistAMuon.end ();++i)
	delete *i;
    fHistAMuon.clear ();
    for (vector<TH3D*>::iterator i = fHistAHadron.begin ();
	 i!=fHistAHadron.end ();++i)
	delete *i;
    fHistAHadron.clear ();

    fZPlane.clear ();
    fZstart = 0;
    fZstop = 0;
    fNhist = 0;
}



void TPlotter::SetShowerZenith (float zenith) {

    cout << " RU: zenith/deg: " << zenith/deg << endl;

    fZenith = zenith*rad;
    fCosZenith = cos (zenith);
    fSinZenith = sin (zenith);
}



 void TPlotter::SetShowerAzimuth (float azimuth) {

     cout << " RU: azimuth/deg: " << azimuth/deg << endl;

     fAzimuth = azimuth*rad;
     fCosAzimuth = cos (azimuth);
     fSinAzimuth = sin (azimuth);
}



void TPlotter::SetShowerHeader (const crs::MEventHeader &header) {

    // to flag the track of the primary particle
    fPrimaryTrack = true;

    SetEvent (header.GetEventNumber ());
    SetPrimary ((int)header.GetParticleId ());
    SetShowerZenith (header.GetTheta ()*rad);
    SetShowerAzimuth (header.GetPhi ()*rad);    
    SetShowerEnergy (header.GetEnergy ()*GeV);
    SetHeightFirstInt (header.GetZFirst ()*cm);
    SetObservationLevel (header.GetObservationHeight (
			     header.GetNObservationLevels ()-1)*cm);

}
    
    
void TPlotter::SetShowerTrailer (const crs::MEventEnd &trailer) {
    
    //cout << " RU: TRAILER " << endl;
    cout << " RU: xmax: " << trailer.GetXmax ()*g/cm/cm << endl;
    SetXmax (trailer.GetXmax ()*g/cm/cm);
}


void TPlotter::Init (int nBins, PlotMode mode) {
    
    //cout << " RU: INIT " << endl;

    double zFirst = 0, zLast = 0;
    double zUnit = 1;
    string ModeStr;
    string UnitStr;

    fZPlane.clear ();

    switch (mode) {

	case eHeight:
	    zFirst = fHeightFirstInt;
	    zLast = fObservationLevel;
	    zUnit = km;
	    ModeStr = "height z = ";
	    UnitStr = "km";
	    break;
	    
	case eDepth:
	    zFirst = thick_ (fHeightFirstInt/cm)*g/cm/cm / fCosZenith;
	    zLast = thick_ (fObservationLevel/cm)*g/cm/cm / fCosZenith;
	    zUnit = g*cm*cm;
	    ModeStr = "slant depth = ";
	    UnitStr = "g/cm^2";
	    break;
	    
	    /*case eTime:
	    zFirst = 0 * ns;
	    zLast = (fHeightFirstInt-fObservationLevel) / fCosZenith
		/ crs::cSpeedOfLight;
	    zUnit = ns;
	    ModeStr = "time t = ";
	    UnitStr = "ns";
	    break;*/
	    
    }
    
    // decide where/what to plot
    SetZfirst (zFirst);
    SetZlast (zLast);
    SetNHistograms (nBins);
    SetPlotMode (mode);
    

    ostringstream file_name;
    file_name << fFileName 
	      << "_"
	      << fEventNo
	      << ".root";
    
    fFile = TFile::Open (file_name.str ().c_str (), "RECREATE");
    fFile->cd ();
    
    // create histogramms
    for (int i=0; i<fNhist; ++i) {
	
	double zSlice = (zFirst - (zFirst-zLast)/(nBins-1)*i);

	
	/* ----------------
	   height of planes in shower frame, with (0/0/0) at core (OBS level)
	*/
	double zPlane = 0;
	switch (mode) {
	    
	    case eHeight:
		zPlane = (zFirst-zLast) 
		    - (zFirst-zLast)/(nBins-1)*i;
		break;

	    case eDepth:
		zPlane = fHeightFirstInt 
		    - heigh_ (zSlice/g*cm*cm * fCosZenith);
		break;

		/*case eTime:
		zPlane = (zFirst-zLast) 
		    - zSlice * crs::cSpeedOfLight * fCosZenith;
		    break;		*/
	}
	zPlane /= fCosZenith;
	//cout << " RU: " << zPlane/m << " " << zSlice/zUnit << endl;
	fZPlane.push_back (zPlane);
	
	
	
	ostringstream h_title;
	h_title << ModeStr
		<< zSlice/zUnit
		<< " " << UnitStr;
	
	
	///cout << " RU: " << i << " " << h_title.str () << endl;
	
	
	string paricle_name [4] = {"Gamma", "EM", "Muon", "Hadron"};
	//string paricle_name [2] = {"Electron", "Positron"};
	for (int j=0; j<4; j++) {
	
	    ostringstream hT_name, hE_name, hL_name, hA_name;
	    hT_name << "hT_"
		    << i
		    << "_"
		    << paricle_name [j];
	    hE_name << "hE_"
		    << i
		    << "_"
		    << paricle_name [j];
	    hL_name << "hL_"
		    << i
		    << "_"
		    << paricle_name [j];
	    hA_name << "hA_"
		    << i
		    << "_"
		    << paricle_name [j];
	    
	    TH3D *histT = new TH3D (hT_name.str ().c_str (),
				    h_title.str ().c_str (),
				    50, -2., 5., // log10 (time/ns)
				    20, -2., 2.,  // log10 (r/rm) 
				    50, -4., 2.);// log10 (e/GeV)

	    TH2D *histE = new TH2D (hE_name.str ().c_str (),
				    h_title.str ().c_str (),
				    50, -4., 2.,  // log10 (e/GeV)
				    20, 0., 2.);  // log10 (r/rm) 
				    
	    TH2D *histL = new TH2D (hL_name.str ().c_str (),
				    h_title.str ().c_str (),
				    50, -2., 2.,  // log10 (r/rm)
				    10, -4., 2.);// log10 (e/GeV)

	    TH3D *histA = new TH3D (hA_name.str ().c_str (),
				    h_title.str ().c_str (),
				    50, 0., 90., // [deg]
				    20, -2., 2.,  // log10 (r/rm) 
				    50, -4., 2.);// log10 (e/GeV)
				    
	    
	    histT->SetLineColor (j+2);
	    histT->SetXTitle ("log_{10}(time/ns)");
	    histT->SetYTitle ("counts");

	    histE->SetLineColor (j+2);
	    histE->SetXTitle ("log_{10}(E/GeV)");
	    histE->SetYTitle ("counts");

	    histL->SetLineColor (j+2);
	    histL->SetXTitle ("log_{10}(dist/m)");
	    histL->SetYTitle ("counts");

	    histA->SetLineColor (j+2);
	    histA->SetXTitle ("cos(theta)");
	    histA->SetYTitle ("counts");
	    
	    switch (j) {
		case 0:
		    fHistTGamma.push_back (histT);
		    fHistEGamma.push_back (histE);
		    fHistLGamma.push_back (histL);
		    fHistAGamma.push_back (histA);
		    break;

		case 1:
		    fHistTElectron.push_back (histT);
		    fHistEElectron.push_back (histE);
		    fHistLElectron.push_back (histL);
		    fHistAElectron.push_back (histA);
		    break;
/*
		case 1:
		    fHistTPositron.push_back (histT);
		    fHistEPositron.push_back (histE);
		    fHistLPositron.push_back (histL);
		    fHistAPositron.push_back (histA);
		    break;
*/

		case 3:
		    fHistTMuon.push_back (histT);
		    fHistEMuon.push_back (histE);
		    fHistLMuon.push_back (histL);
		    fHistAMuon.push_back (histA);
		    break;
		case 4:
		    fHistTHadron.push_back (histT);
		    fHistEHadron.push_back (histE);
		    fHistLHadron.push_back (histL);
		    fHistAHadron.push_back (histA);
		    break;

	    }
	}
    }


    //cout << "x/F:y/F:z/F:h/F:d/F:t/F:rm/F:i/I:str/C   RURU " << endl;
    //cout << "x/F:y/F:z/F:i/I:str/C   RURU " << endl;

}



void TPlotter::Write () {

    //cout << " RU: Write " << endl;
    
    ostringstream gif_dir;
    gif_dir << fFileName 
	    << "_"
	    << fEventNo;
    

    fXmax /= g*cm*cm;
    fXmax /= fCosZenith;

    fFile->cd ();
    TTree *run = new TTree ("shower", "shower");
    run->Branch ("xmax"    , &fXmax            , "xmax/F");
    run->Branch ("evntno"  , &fEventNo         , "eventno/I");
    run->Branch ("obslevel", &fObservationLevel, "obslevel/D");
    run->Branch ("firstint", &fHeightFirstInt  , "firstint/D");
    run->Branch ("zenith"  , &fZenith          , "zenith/F");
    run->Branch ("azimuth" , &fAzimuth         , "azimuth/F");
    run->Branch ("energy"  , &fEnergy0         , "energy/F");
    run->Branch ("primary" , &fPrimary         , "primary/I");
    run->Fill ();


    Float_t time, depth, height, age, rm, temp, density, pressure;
    TTree *tree = new TTree ("data", "data");
    tree->Branch ("time"  , &time,   "time/F");
    tree->Branch ("depth" , &depth,  "depth/F");
    tree->Branch ("height", &height, "height/F");
    tree->Branch ("age"   , &age,    "age/F");
    tree->Branch ("rm", &rm, "rm/F");
    tree->Branch ("temp", &temp, "temp/F");
    tree->Branch ("pressure", &pressure, "pressure/F");
    tree->Branch ("density", &density, "density/F");
    /*
    TH3D *hTelectron = new TH3D ();
    TH2D *hEelectron = new TH2D ();
    TH2D *hLelectron = new TH2D ();
    TH3D *hAelectron = new TH3D ();
    TH3D *hTPositron = new TH3D ();
    TH2D *hEPositron = new TH2D ();
    TH2D *hLPositron = new TH2D ();
    TH3D *hAPositron = new TH3D ();
    */
    TH3D *hTgamma = new TH3D ();
    TH2D *hEgamma = new TH2D ();
    TH2D *hLgamma = new TH2D ();
    TH3D *hAgamma = new TH3D ();
    TH3D *hTelectron = new TH3D ();
    TH2D *hEelectron = new TH2D ();
    TH2D *hLelectron = new TH2D ();
    TH3D *hAelectron = new TH3D ();
    TH3D *hTmuon = new TH3D ();
    TH2D *hEmuon = new TH2D ();
    TH2D *hLmuon = new TH2D ();
    TH3D *hAmuon = new TH3D ();
    TH2D *hEhadron = new TH2D ();
    TH3D *hThadron = new TH3D ();
    TH2D *hLhadron = new TH2D ();
    TH3D *hAhadron = new TH3D ();
    
    // the histogramms
    /*
    tree->Branch ("hTelectron"     , "TH3D",   &hTelectron, 32000, 0);
    tree->Branch ("hEelectron"     , "TH2D",   &hEelectron, 32000, 0);
    tree->Branch ("hLelectron"     , "TH2D",   &hLelectron, 32000, 0);
    tree->Branch ("hAelectron"     , "TH3D",   &hAelectron, 32000, 0);

    tree->Branch ("hTPositron"     , "TH3D",   &hTPositron, 32000, 0);
    tree->Branch ("hEPositron"     , "TH2D",   &hEPositron, 32000, 0);
    tree->Branch ("hLPositron"     , "TH2D",   &hLPositron, 32000, 0);
    tree->Branch ("hAPositron"     , "TH3D",   &hAPositron, 32000, 0);
    */
    tree->Branch ("hTgamma"     , "TH3D",   &hTgamma, 32000, 0);
    tree->Branch ("hEgamma"     , "TH2D",   &hEgamma, 32000, 0);
    tree->Branch ("hLgamma"     , "TH2D",   &hLgamma, 32000, 0);
    tree->Branch ("hAgamma"     , "TH3D",   &hAgamma, 32000, 0);

    tree->Branch ("hTelectron"     , "TH3D",   &hTelectron, 32000, 0);
    tree->Branch ("hEelectron"     , "TH2D",   &hEelectron, 32000, 0);
    tree->Branch ("hLelectron"     , "TH2D",   &hLelectron, 32000, 0);
    tree->Branch ("hAelectron"     , "TH3D",   &hAelectron, 32000, 0);

    tree->Branch ("hTmuon"   , "TH3D",   &hTmuon, 32000, 0);
    tree->Branch ("hEmuon"   , "TH2D",   &hEmuon, 32000, 0);
    tree->Branch ("hLmuon"   , "TH2D",   &hLmuon, 32000, 0);
    tree->Branch ("hAmuon"   , "TH3D",   &hAmuon, 32000, 0);

    tree->Branch ("hThadron" , "TH3D",   &hThadron, 32000, 0);
    tree->Branch ("hEhadron" , "TH2D",   &hEhadron, 32000, 0);
    tree->Branch ("hLhadron" , "TH2D",   &hLhadron, 32000, 0);
    tree->Branch ("hAhadron" , "TH3D",   &hAhadron, 32000, 0);
    
    

    double firstZ = fHeightFirstInt/km;
    double lastZ = fObservationLevel/km;
    double dZ = (lastZ-firstZ)/fNhist;
    
    double firstX = thick_ (fHeightFirstInt/cm)*g/cm/cm / fCosZenith;
    double lastX = thick_ (fObservationLevel/cm)*g/cm/cm / fCosZenith;
    double dX = (lastX-firstX)/fNhist;
    
    double firstT = 0*ns;
    double lastT = (fHeightFirstInt-fObservationLevel) / fCosZenith
	/ crs::cSpeedOfLight;
    double dT = (lastT-firstT)/fNhist;
    
    
    for (int i=0; i<fNhist; ++i) {
	
	time   = firstT + dT * i;
	height = firstZ + dZ * i;
	depth = firstX + dX * i;

	switch (fMode) {
	    case eHeight:
		depth  = thick_ (height/cm)*g/cm/cm / fCosZenith;
		break;
		
	    case eDepth:
		height = heigh_ (depth*fCosZenith/g*cm*cm)*cm;
		time = firstT + (fHeightFirstInt-height) / fCosZenith
		    / crs::cSpeedOfLight;
		break;
		
		/*case eTime:
		depth = thick_ (height/cm)*g/cm/cm / fCosZenith;
		break;*/
	}

	double T = Temperature (height);	
	rm = MoliereRadius (depth*fCosZenith, T)/m;

	static double g_earth = 9.81 * m/(s*s);
	double Pressure = depth * fCosZenith * g_earth;
	static double MeanAirDensity = 28.95 * g/mol;
	static double R_gas = 8.314 * joule/mol/kelvin;
	double Density = Pressure / (R_gas*T/MeanAirDensity);

	temp = T;
	pressure = Pressure/bar;
	density = Density/g*cm*cm*cm;

	time /= ns;
	height /= km;
	depth /= g*cm*cm;

	age = 3.0 / (1.0 + 2.0*fXmax/depth); 
	
/*
 	*hTem = *(fHistTEM [i]);
 	*hEem = *(fHistEEM [i]);
 	*hLem = *(fHistLEM [i]);
 	*hAem = *(fHistAEM [i]);

 	*hTem = *(fHistTEM [i]);
 	*hEem = *(fHistEEM [i]);
 	*hLem = *(fHistLEM [i]);
 	*hAem = *(fHistAEM [i]);
	*/

 	*hTgamma = *(fHistTGamma [i]); 
 	*hEgamma = *(fHistEGamma [i]);
 	*hLgamma = *(fHistLGamma [i]);
 	*hAgamma = *(fHistAGamma [i]);

 	*hTelectron = *(fHistTElectron [i]); 
 	*hEelectron = *(fHistEElectron [i]);
 	*hLelectron = *(fHistLElectron [i]);
 	*hAelectron = *(fHistAElectron [i]);

	*hTmuon = *(fHistTMuon [i]);
	*hEmuon = *(fHistEMuon [i]);
	*hLmuon = *(fHistLMuon [i]);
	*hAmuon = *(fHistAMuon [i]);

	*hEhadron = *(fHistEHadron [i]);
	*hThadron = *(fHistTHadron [i]);
	*hLhadron = *(fHistLHadron [i]);
	*hAhadron = *(fHistAHadron [i]);

	
	tree->Fill ();

	// delete hists
	delete &(*fHistTGamma [i]);
	delete &(*fHistEGamma [i]);
	delete &(*fHistLGamma [i]);
	delete &(*fHistAGamma [i]);

	delete &(*fHistTElectron [i]);
	delete &(*fHistEElectron [i]);
	delete &(*fHistLElectron [i]);
	delete &(*fHistAElectron [i]);
	/*
	delete &(*fHistTPositron [i]);
	delete &(*fHistEPositron [i]);
	delete &(*fHistLPositron [i]);
	delete &(*fHistAPositron [i]);
	*/

	delete &(*fHistTHadron [i]);
	delete &(*fHistEHadron [i]);
	delete &(*fHistLHadron [i]);
	delete &(*fHistAHadron [i]);

	delete &(*fHistEMuon [i]);
	delete &(*fHistTMuon [i]);
	delete &(*fHistLMuon [i]);
	delete &(*fHistAMuon [i]);
	
    }
    

    run->Write ();
    tree->Write ();
    fFile->Close ();

    
    // clean up
    fHistTGamma.clear ();    
    fHistEGamma.clear ();    
    fHistLGamma.clear ();    
    fHistAGamma.clear ();    

    fHistTElectron.clear ();    
    fHistEElectron.clear ();    
    fHistLElectron.clear ();    
    fHistAElectron.clear ();    

    /*
    fHistTPositron.clear ();    
    fHistEPositron.clear ();    
    fHistLPositron.clear ();    
    fHistAPositron.clear ();    
    */
 
    fHistTMuon.clear ();    
    fHistEMuon.clear ();    
    fHistLMuon.clear ();    
    fHistAMuon.clear ();    

    fHistTHadron.clear ();    
    fHistEHadron.clear ();    
    fHistLHadron.clear ();    
    fHistAHadron.clear ();    


    Clear ();
}




//#define Rotate(x, y, sx, sy, cosa, sina) { \
inline
void TPlotter::Rotate (double x, double y, double z,
		       double &sx, double &sy, double &sz,
		       int inverse) {
    

    double sx_ =             x*fCosAzimuth + inverse * y*fSinAzimuth;
    double sy_ = - inverse * x*fSinAzimuth +           y*fCosAzimuth; 
    double sz_ =   z;


    sx =             sx_*fCosZenith - inverse * sz_*fSinZenith;
    sy = sy_;
    sz = + inverse * sx_*fSinZenith +           sz_*fCosZenith; 

}


void TPlotter::AddTrack (const crs::CParticle &pre, 
			 const crs::CParticle &post) {
    

    /*
      Skip the track of the primary particle, which is in a different 
      reference system as the shower !!!
     */
    if (fPrimaryTrack) {
	cout << " RU: Primary track " << endl;
	fPrimaryTrack = false;
	return;
    }

    /*
    cout << post.x/m << " " 
	 << post.y/m << " " 
	 << post.z/m << " " 
	 << (int)pre.particleId << " " 
	 << " RURU "
	 << endl;
    return;
    */
  
    // direction of particle
    //double cosTheta = pre.cosTheta;
    //double theta = acos (cosTheta) * rad;


    /*   ------------------------
      convert to shower frame
    */
    double pre_shX, pre_shY, pre_shZ;
    Rotate (pre.x, pre.y, fHeightFirstInt-pre.z, 
	    pre_shX, pre_shY, pre_shZ, +1);


    double post_shX, post_shY, post_shZ;
    Rotate (post.x, post.y, fHeightFirstInt-post.z, 
	    post_shX, post_shY, post_shZ, +1);

/*
    cout << post_shX/m << " " 
	 << post_shY/m << " " 
	 << post_shZ/m << " " 
	 << (int)pre.particleId << " " 
	 << " RURU "
	 << endl;
    return;
*/  

    

    /*   ------------------------
      direction of particle in shower frame
    */
    double dX = post_shX - pre_shX; 
    double dY = post_shY - pre_shY; 
    double dZ = post_shZ - pre_shZ; 
    double length = sqrt (dX*dX + dY*dY + dZ*dZ);

    if (length==0) {
	/*   ------------------------ 
	  This happens for particles that are NOT tracked by Corsika
	  e.g. pi0 -> decay immediatly
	  SKIP them
	 */
	return;
    }	
    double theta = acos (dZ/length) * rad;



    /*
    cout << " RU: theta: " << theta/deg 
	 << " " << ThetaInShwFrame/deg
	 << endl;
    */
    

    /*   ------------------------
      what to look for (planes in time, depth, height)
    */
    /*
    double zPre = 0;
    double zPost = 0;
    switch (fMode) {
	case eHeight:
	    zPre = pre_shZ * cm/km;
	    zPost = post_shZ * cm/km;
	    break;

	case eDepth:
	    zPre = pre_depth;
	    zPost = post_depth;
	    break;

	case eTime:
	    zPre = pre.time * s/ns;
	    zPost = post.time * s/ns;
	    break;
    }
    */
    
    
    
    // check for crossing of planes
    for (int i=0; i<fNhist; ++i) {
	
	double zPlane = fZPlane [i]; // / fCosZenith;
	
	// downward going particle
	//if (pre_shZ  post_shZ) {
	    
	    //cout << " RU:" << zPlane << " " << zPre << " " << zPost << endl;

	    //if (!((zPre<zPlane && zPost>zPlane) || 
	    //(zPost<zPlane && zPre>zPlane)))
	    if (!((pre_shZ<zPlane && post_shZ>zPlane) || 
		  (post_shZ<zPlane && pre_shZ>zPlane)))
		continue;
	    
	    
	    double delta = (pre_shZ-zPlane) / (pre_shZ-post_shZ);	   
	    
	    double e_post = post.energy 
		- gParticleMass [(int)post.particleId-1];
	    double e_pre = pre.energy 
		- gParticleMass [(int)pre.particleId-1];
	    
	    double dT = post.time - pre.time;
	    double dE = e_post - e_pre;
	    double dW = post.weight - pre.weight;
	    
	    double time = (pre.time + dT * delta) * (s/ns);
	    double energy = e_pre + dE * delta;
	    double weight = 1;
	    double x = pre_shX + dX * delta;
	    double y = pre_shY + dY * delta;
	    double z = pre_shZ + dZ * delta;   // NOT NEEDED
	    if (fThinning) {
		weight = pre.weight + dW * delta;
	    }
	    
	    double height = fHeightFirstInt - zPlane * fCosZenith;
	    double depth = thick_ (height/cm) * (g/cm/cm);
	    
	    double r = sqrt (pow (x, 2) + pow (y, 2));
	    double l = sqrt (pow (r, 2) + pow (z, 2));
	    
	    //cout << "   RU:  d:" << delta << " t:" << time;
	    ///cout << " x:" << x << " y:" << y << " z:" << z;
	    ///cout << " l:" << l << " r:" << r;	    
	    
	    time -= l / crs::cSpeedOfLight;
	    //time -= zPlane/crs::cSpeedOfLight;
	    
	    
	    if (time<=0) {
		/*
		  cout << " invalid TIME: " 
		  << time << " ns " << endl;
		  pre.Dump ();
		  post.Dump ();
		*/
		continue;
	    }
	    
	    if (energy<=0) {
		/*
		  cout << " invalid ENERGY: " 
		  << energy << " GeV" << endl;
		  pre.Dump ();
		  post.Dump ();
		*/
		continue;
	    }

	    if (r<=0) {
		/*
		  cout << " invalid RADIUS: " 
		  << r << " cm " << endl;
		  pre.Dump ();
		  post.Dump ();
		*/
		continue;
	    }
	    
	    
	    double T = Temperature (height);
	    double rm = MoliereRadius (depth, T);

	    /*
	    cout << height/m << " "
		 << depth/g*cm*cm << " "
		 << T/kelvin << " " 
		 << rm/m << " "
		 << endl;
	    */
	    
	    /*
	    cout << x/m << " " 
		 << y/m << " " 
		 << z/m << " " 
		 << height/m << " "
		 << depth/g*cm*cm << " "
		 << time/ns << " " 
		 << rm/m << " "
		 << (int)pre.particleId << " " 
		 << " RURU "
		 << endl;
	    */
	    

	    int particleId = (int)pre.particleId;
	    switch (particleId) {

		case 1:
		    fHistTGamma [i]->Fill (log10 (time), 
					   log10 (r/rm),
					   log10 (energy),
					   weight);

		    fHistEGamma [i]->Fill (log10 (energy), log10 (r/rm),
					   weight);
		    
		    fHistLGamma [i]->Fill (log10 (r/rm), log10 (energy),
					   weight);
		    
		    fHistAGamma [i]->Fill (theta/deg,
					   log10 (r/rm),
					   log10 (energy),
					   weight);
		    break;

		case 2:
		case 3:
		    fHistTElectron [i]->Fill (log10 (time), 
					 log10 (r/rm),
					 log10 (energy),
					 weight);

		    fHistEElectron [i]->Fill (log10 (energy), log10 (r/rm),
					 weight);
		    
		    fHistLElectron [i]->Fill (log10 (r/rm), log10 (energy),
					 weight);
		    
		    fHistAElectron [i]->Fill (theta/deg,
					 log10 (r/rm),
					 log10 (energy),
					 weight);
		    break;

/*		    
		case 3:
		    fHistTPositron [i]->Fill (log10 (time), 
					log10 (r/rm),
					log10 (energy),
					weight);

		    fHistEPositron [i]->Fill (log10 (energy), log10 (r/rm),
					weight);

		    fHistLPositron [i]->Fill (log10 (r/rm), log10 (energy),
					weight);

		    fHistAPositron [i]->Fill (theta/deg,
					log10 (r/rm),
					log10 (energy),
					weight);
		    break;
*/	    
	
		case 5:
		case 6:
		    fHistTMuon [i]->Fill (log10 (time), 
					  log10 (r/rm),
					  log10 (energy),
					  weight);

		    fHistEMuon [i]->Fill (log10 (energy), log10 (r/rm),
					  weight);

		    fHistLMuon [i]->Fill (log10 (r/rm), log10 (energy),
					  weight);

		    fHistAMuon [i]->Fill (theta/deg,
					  log10 (r/rm),
					  log10 (energy),
					  weight);
		    break;

		default:

		    if (particleId>7 && particleId<54) {
			fHistTHadron [i]->Fill (log10 (time), 
						log10 (r/rm),
						log10 (energy),
						weight);

			fHistEHadron [i]->Fill (log10 (energy), log10 (r/rm),
						weight);

			fHistLHadron [i]->Fill (log10 (r/rm), log10 (energy),
						weight);

			fHistAHadron [i]->Fill (theta/deg,
						log10 (r/rm),
						log10 (energy),
						weight);
		     }

		    break;
	    }
	    
    }

}



double TPlotter::MoliereRadius (double VDepth,
				double Temperature) {

    static double E_scale = 21.*MeV;
    static double X_rad_length = 37*g/cm/cm;
    static double e_crit = 81.*MeV;
    static double R_m = E_scale * X_rad_length / e_crit;   // [g/cm2]

    // correct for 2 radiation length
    VDepth -= 2. * X_rad_length * fCosZenith;
    if (VDepth<=0) VDepth = 1.e-8;

    static double MeanAirDensity = 28.95 * g/mol;
    static double g_earth = 9.81 * m/(s*s);
    static double R_gas = 8.314 * joule/mol/kelvin;
    static double Na = 6.022e23 / mol;
    double Pressure = VDepth * g_earth;
    double Density = Pressure / (R_gas*Temperature/MeanAirDensity);

    return R_m / Density;
}


/*
double TPlotter::MoliereRadius (double Temperature, 
				double Pressure) {
    
    static double exponent = 1./5.25588;
    
    Temperature /= kelvin;
    Pressure /= millibar;
    
    double CorrectedPressure = Pressure - 73.94*fCosZenith;
    if (CorrectedPressure<=0) {
	//cout << " RU: CorrectedPressure<0: " << CorrectedPressure << endl;
	//cout << " T: " << Temperature << " P: " << Pressure << endl;
	CorrectedPressure = 0.01;
	Pressure = CorrectedPressure + 73.94*fCosZenith;
    }

    double Rm = 272.5 * Temperature * 
	pow (CorrectedPressure/Pressure, exponent) / CorrectedPressure;
    
    return Rm * m;
}
*/

double TPlotter::Temperature (double height) {
/*
  h1(km)   h2(km) 	dT/dh (K/km)
  0 	11 	-6.5
  11 	20 	0.0
  20 	32 	1.0
  32 	47 	2.8
  47 	51 	0.0
  51 	71 	-2.8
  71 	84.852 	-2.0

  Note: 84.852 km geopotential=86 km geometric

  These data along with the sea level standard values of
  Sea level pressure = 101325 N/m2
  Sea level temperature = 288.15 K
  Hydrostatic constant = 34.1631947 kelvin/km
  define the atmosphere. The sea level density of 1.225 kg/m3 is 
  derived from the fundamental quantities above
*/

    double temp = 288.15*kelvin;

    if (height < 11*km) {

	temp += -6.5*(kelvin/km) * height;
	
    } else if (height < 20*km) {

	temp += -6.5*(kelvin/km) * 11*km;

    } else if (height < 32*km) {

	temp += -6.5*(kelvin/km) * 11*km
	    + 1.0*(kelvin/km) * (height-20*km);

    } else if (height < 47*km) {

	temp += -6.5*(kelvin/km) * 11*km
	    + 1.0*(kelvin/km) * (32*km-20*km)
	    + 2.8*(kelvin/km) * (height-32*km);

    } else if (height < 51*km) {

	temp += -6.5*(kelvin/km) * 11*km
	    + 1.0*(kelvin/km) * (32*km-20*km)
	    + 2.8*(kelvin/km) * (47*km-32*km);

    } else if (height < 71*km) {

	temp += -6.5*(kelvin/km) * 11*km
	    + 1.0*(kelvin/km) * (32*km-20*km)
	    + 2.8*(kelvin/km) * (47*km-32*km)
	    - 2.8*(kelvin/km) * (height-51*km);

    } else if (height < 84.852*km) {

	temp += -6.5*(kelvin/km) * 11*km
	    + 1.0*(kelvin/km) * (32*km-20*km)
	    + 2.8*(kelvin/km) * (47*km-32*km)
	    - 2.8*(kelvin/km) * (71*km-51*km)
	    - 2.0*(kelvin/km) * (height-71*km);

    } else {

	temp += -6.5*(kelvin/km) * 11*km
	    + 1.0*(kelvin/km) * (32*km-20*km)
	    + 2.8*(kelvin/km) * (47*km-32*km)
	    - 2.8*(kelvin/km) * (71*km-51*km)
	    - 2.0*(kelvin/km) * (84.852*km-71*km);

    }

    return temp;
	
}
