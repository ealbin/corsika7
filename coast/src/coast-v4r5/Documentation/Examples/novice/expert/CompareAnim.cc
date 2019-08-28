#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLine.h>
#include <TStyle.h>

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
using namespace std;




void CompareAnim (const string &fileName, 
		  bool Normalize, 
		  double Age, 
		  bool FirstDraw, int color, int marker,
		  bool PlotFit,
		  TLegend *legend);


TH1 *Interpolate3D (TH1* h1, double age1,
		    TH1* h2, double age2,
		    double age);

Double_t NKG (Double_t *x, Double_t *p);



int main (int argc, char **argv) {

    gStyle->SetPalette (1);
    gStyle->SetPadRightMargin (.05);
    gStyle->SetOptStat (0);
    gStyle->SetOptTitle (0);

    if (argc<2) {
	cerr << "      please specify data file (s) " << endl;
	return 1;
    }

    vector <string> files;
    for (int i=1; i<argc; i++) {
	
	string fname (argv [i]);
	files.push_back (fname);
	cout << " File: " << fname << endl;
    }


    bool Normalize = true; 

    /*
      if (argc==3) {
	
      istringstream sMode (argv [2]);
      sMode >> Normalize;	
      }
    */
    
    if (Normalize) {
      cout << " --> Normalize all spectra " << endl;
    } else {
      cout << " --> Do not normalize spectra " << endl;
    }
    


    
    /*
      prepare commands
    */

    ostringstream gif_dir;
    gif_dir << "TMPDIF";//fileName.substr (0, fileName.rfind (".root")); 

    string DATRunEvent = "age";/*fileName.substr (fileName.rfind ('/')+1,
					  fileName.rfind (".root") -
					  fileName.rfind ('/') - 1);*/

    cout << " RU: Make Directory: " << gif_dir.str () << endl;
    gSystem->MakeDirectory (gif_dir.str ().c_str ());
    



    /*
      Interpolate series
    */
    int nAge = 3;
    double Age1 = .8;
    double dAge = .2; 

    int delay = 100/nAge * 15; // total time: 15*(s/100)*100ms = 15s
    
    
    
    TCanvas *cCanvas = new TCanvas ("cCanvas", "one slice", 800, 700);
    //cCanvas->SetLogy (1);
    TLegend *legend = new TLegend (.3, .6, .999, .999);

    for (int iAge=0; iAge<nAge; iAge++) {
	
	double Age = Age1 + dAge * iAge;
	
	/*
	  create canvas
	*/
		
	ostringstream sAge;
	sAge << " age = " << setw (6) << setfill ('0') 
	     << fixed << setprecision (4) << left << Age;
	//TPaveLabel *AgeText = new TPaveLabel (.4, .473, .6, .527, 
	//sAge.str ().c_str (), "");
	//AgeText->SetFillColor (70);
	//AgeText->Draw ();


	/*
	TPaveLabel *Title = new TPaveLabel (.0, .95, 1.0, 1.0, sAge.str ().c_str ());
	//sTitle.str ().c_str ());
	Title->SetFillColor (188);
	Title->Draw ();
	TPad *cTmp = new TPad ("cTmp", "cTmp", .0, .0, 1.0, .95, -1, 0, 0);
	cTmp->SetLogy (1);
	cTmp->Draw ();
	cTmp->cd ();
	*/

	static bool FirstDraw = true;
	int color = 100;
	for (vector <string>::iterator iFile = files.begin ();
	     iFile != files.end ();
	     ++iFile) {
	    
	    int marker = 20 + iAge;
	    bool PlotFit = (iFile==(files.end ()-1));
	    
	    CompareAnim ((*iFile), Normalize, Age, 
			 FirstDraw, color, marker, PlotFit, legend);
	    FirstDraw = false;
	    color -= (49/6);
	}

	/*
	  plot parameterisation
	*/
	/*
	TF1 *Para = new TF1 ("para", 
	   "pow(10,x+3)/((pow(10,x+3)+[0])*(pow(10,x+3)+[1]))", -10, 10);  
	Para->SetParameter (0, 6.425-1.532*Age);
	Para->SetParameter (1, 168.13-42.14*Age);
	*/
	//Para->Draw ("lsame");
	/*
	TF1 *Para = new TF1 ("para", NKG, -10, 10, 2);
	Para->SetParNames ("norm", "age");
	Para->SetParameter (0, 1);
	Para->FixParameter (1, Age);
	Para->DrawCopy ("lsame");
	*/

	legend->Draw ();

	// create and save gif
	/*
	ostringstream gif_name, eps_name;
	gif_name << gif_dir.str () << "/h_";
	gif_name.width (3);
	gif_name.fill ('0');
	gif_name << right << iAge << ".gif";
	
	eps_name << gif_dir.str () << "/h_";
	eps_name.width (3);
	eps_name.fill ('0');
	eps_name << right << iAge << ".eps";
    
	ostringstream eps2gif_cmd1;
	eps2gif_cmd1 << "pstopnm -ppm -xborder 0 -yborder 0 "
		     << "-xsize 900 -ysize 700 -portrait"
		     << " cTmp.eps";
    
	ostringstream eps2gif_cmd2;
	eps2gif_cmd2 << "ppmtogif cTmp.eps001.ppm > "
		     << gif_name.str ();
	
	ostringstream root_save;
	root_save << gif_dir.str () << "/h_"
		  << setw (3) << setfill ('0')
		  << right << iAge << ".C";
	
	cout << " RU: save gif " << eps2gif_cmd1.str () 
	     << " " << eps2gif_cmd2.str () << endl;
	cCanvas->Print ("cTmp.eps");
	cCanvas->Print (root_save.str ().c_str ());
	gSystem->Exec (eps2gif_cmd1.str ().c_str ());
	gSystem->Exec (eps2gif_cmd2.str ().c_str ());
	
	delete cCanvas;
	*/
    }

    cCanvas->Print ("canvas.C");
    cCanvas->Print ("canvas.root");
    cCanvas->Print ("canvas.eps");

    /*
    ostringstream makeanim;
    makeanim << "gifsicle --delay=" << delay << " `ls "
	     << gif_dir.str () << "/h*.gif` > "
	     << DATRunEvent << "_comp.gif";
    
    ostringstream rm_dir;
    rm_dir << " rm -rf " << gif_dir.str ().c_str ();
    

    cout << makeanim.str () << endl;


    gSystem->Exec (makeanim.str ().c_str ());
    //gSystem->Exec (rm_dir.str ().c_str ());
    */

    return 0;
}





void CompareAnim (const string &fileName, 
		  bool Normalize, 
		  double Age, 
		  bool FirstDraw, int color, int marker,
		  bool PlotFit, TLegend *legend) {
    
    
    int cut1 = 0;
    if (fileName.rfind ('/')!=string::npos)
	cut1 = fileName.rfind ('/');
    
    int cut2 = fileName.length () - cut1;
    if (fileName.rfind (".root")!=string::npos)
	cut2 = fileName.rfind (".root") - cut1;
    
    string DATRunEvent = fileName.substr (cut1, cut2);
	
    TFile *file = TFile::Open (fileName.c_str (), "READ");

    if (file==NULL || file->IsZombie ())
	return;
    
    file->cd ();
    TTree *run = (TTree*) file->Get ("shower");
    TTree *tree = (TTree*) file->Get ("data");


    Float_t Xmax, Zenith, Azimuth, Energy0;
    Double_t ObservationLevel, HeightFirstInt;
    Int_t EventNo, Primary;

    run->SetBranchAddress ("xmax"    , &Xmax  );
    run->SetBranchAddress ("evntno"  , &EventNo);
    run->SetBranchAddress ("obslevel", &ObservationLevel);
    run->SetBranchAddress ("firstint", &HeightFirstInt );
    run->SetBranchAddress ("zenith"  , &Zenith        );
    run->SetBranchAddress ("azimuth" , &Azimuth        );
    run->SetBranchAddress ("energy"  , &Energy0    );
    run->SetBranchAddress ("primary" , &Primary   );
    run->GetEntry (0);
    



    Float_t time, depth, height, age, rm;

    tree->SetBranchAddress ("time"  , &time );
    tree->SetBranchAddress ("depth" , &depth );
    tree->SetBranchAddress ("height", &height);
    tree->SetBranchAddress ("age"   , &age );
    tree->SetBranchAddress ("rm"    , &rm );
 
    TH3D *hist3D = 0;
    /*
    TH3D *hTmuon = 0;
    TH3D *hThadron = 0;
    TH2D *hEem = 0;
    TH2D *hEmuon = 0;
    TH2D *hEhadron = 0;
    */
    TH2D *hLem = 0;
    TH2D *hLmuon = 0;
    /*
    TH2D *hLhadron = 0;
    TH3D *hAem = 0;
    TH3D *hAmuon = 0;
    TH3D *hAhadron = 0;
    */
    
    string histUnit = "log_{10}(r/Rm)";
    //string histUnit = "log_{10}(E/GeV)";
    //string histUnit = "log_{10}(t/ns)";
    //string histUnit = "theta/deg";
    
    // the histogramms
    //tree->SetBranchAddress ("hTem"     , &hist3D);
    tree->SetBranchAddress ("hAem"     , &hist3D);

    /*
    tree->SetBranchAddress ("hTmuon"   , &hTmuon);
    tree->SetBranchAddress ("hThadron" , &hThadron);
    tree->SetBranchAddress ("hEem"     , &hEem);
    tree->SetBranchAddress ("hEmuon"   , &hEmuon);
    tree->SetBranchAddress ("hEhadron" , &hEhadron);
    */
    tree->SetBranchAddress ("hLem"     , &hLem);
    tree->SetBranchAddress ("hLmuon"   , &hLmuon);
    /*
    tree->SetBranchAddress ("hLhadron" , &hLhadron);
    tree->SetBranchAddress ("hAem"     , &hAem);
    tree->SetBranchAddress ("hAmuon"   , &hAmuon);
    tree->SetBranchAddress ("hAhadron" , &hAhadron);
    */

    int Nhist = tree->GetEntries ();


    /*
      float AgeMin = 9999, AgeMax = 0;
      for (int j=1; j<Nhist; j++) {

      tree->GetEntry (j);
      
      if (age>AgeMax)
      AgeMax = age;
      
      if (age<AgeMin)
      AgeMin = age;
      
      //cout << " AGE: " << age << endl;
      }  
    */



    /*
      Normalizing function
    */

    // dA = 2*pi*((r+dr)^2-r^2) = 2*pi* (2rdr + dr^2)  =propto=  r
    TF1 *dA = new TF1 ("dA", "pow(10,x)", -10, 10);  

    // dOmega =propto= sin(theta)
    TF1 *dO = new TF1 ("dO", "sin(x*3.141/180.)", 0, 360);
    //TF1 *dO = new TF1 ("dO", "cos(x*3.141/180.)", 0, 360);
    
    //TF1 *Norm = new TF1 ("norm", "1", -1e10, 1e10);
    TF1 *Norm = dA;
    

    /*
      Fit functions
    */
    TF1 *Para = new TF1 ("para", NKG, -10, 10, 3);
    Para->SetParNames ("norm", "age", "rm");
    Para->SetParameter (0, 1);
    Para->SetParameter (1, Age);
    Para->SetParameter (2, 1);
    

    /*
      interpolate in AGE
     */
    int i;
    for (i=1; i<Nhist; i++) {
	
	tree->GetEntry (i);
	    
	if (age>Age)
	    break;
    }
    
    TH3D *h1 = (TH3D*) hist3D->Clone ("h1");
    TH2D *h2 = (TH2D*) hLem->Clone ("h2");
    TH2D *h3 = (TH2D*) hLmuon->Clone ("h3");
    double hAge2 = age;
    
    tree->GetEntry (i-1);
    TH3D *hInter = (TH3D*) Interpolate3D (h1, hAge2, hist3D, age, Age);
    TH2D *hInterLem = (TH2D*) Interpolate3D (h2, hAge2, hLem, age, Age);
    TH2D *hInterLmuon = (TH2D*) Interpolate3D (h3, hAge2, hLmuon, age, Age);
    delete h1;
    delete h2;
    delete h3;
    
    
    
    int NBinsZ = 0;
    int nProj = 1;
    double maxT[nProj];
    double minT[nProj];
    for (int i=0; i<nProj; i++) {
	maxT [i] = 0;
	minT [i] = 1e10;
    }
    
    
    
    if (Normalize) {
	
	for (int iProj=0; iProj<nProj; iProj++) {
	    
	    maxT [iProj] = 1;
		
	    // make projection:
	    NBinsZ = hInter->GetYaxis ()->GetNbins ();
	    int bin1 = 1 + NBinsZ*iProj/nProj;
	    int bin2 = 0 + NBinsZ*(iProj+1)/nProj;
	    //hInter->GetYaxis()->SetRange (bin1, bin2);
	    //TH2D *hist = (TH2D*) hInter->Project3D ("zx");
	    TH1D *hist = (TH1D*) hInter->Project3D ("z");
	    hist->Divide (Norm);

	    double iT = hist->Integral ();
	    if (iT) {
		double histT_min = hist->GetMinimum (.9)/iT;
		double histT_max = hist->GetMaximum ()/iT;
		if (histT_min<minT [iProj])
		    minT [iProj] = histT_min;
		if (histT_max>maxT [iProj])
		    maxT [iProj] = histT_max;
	    }
	}

    } else {

	
	for (int iProj=0; iProj<nProj; iProj++) {
	    
	    minT [iProj] = 1;
	    
	    // make projection:
	    NBinsZ =hInter->GetYaxis ()->GetNbins ();
	    int bin1 = 1 + NBinsZ*iProj/nProj;
	    int bin2 = 0 + NBinsZ*(iProj+1)/nProj;
	    //hInter->GetYaxis()->SetRange (bin1, bin2);
	    //TH2D *hist = (TH2D*) hInter->Project3D ("zx");
	    TH1D *hist = (TH1D*) hInter->Project3D ("z");
	    hist->Divide (Norm);
	    
	    double histT_max = hist->GetMaximum ();
	    if (histT_max>maxT [iProj])
		maxT [iProj] = histT_max;
	    
	}

	for (int iProj=0; iProj<nProj; iProj++)
	    maxT [iProj] *= 1.1;	   

    }




    /*
      do the legend entry
    */
    ostringstream sTitle;
    switch (Primary) {
	case 1:
	    sTitle << "Gamma";
	    break;
	case 5626:
	    sTitle << "Iron";
	    break;
	case 14:
	    sTitle << "Proton";
	    break;
	default:
	    sTitle << "Unknown";
    }

    sTitle << " " 
	   << ", "/*", E_{0}="*/ << scientific << setprecision(0) << Energy0 << "GeV "
	   << ", "/*", #theta= "*/ << fixed << setprecision (0) << Zenith/3.141*180. << "deg "
	   << ", "/*", X_{max}="*/ << fixed << setprecision (0) << Xmax << "g/cm^{2} ";

		
	
    nProj = 1;
    for (int iProj=0; iProj<nProj; iProj++) {
	    
	ostringstream sProj, sProj2, sProj3;
	//sProj << "zx"
	sProj << "z"
	      << "_" << iProj+1 
	      << "_" << int(Age*100)
	      << "_" << DATRunEvent;
	
	sProj2 << "L2"
	      << "_" << iProj+1 
	      << "_" << int(Age*100)
	      << "_" << DATRunEvent;
	
	sProj3 << "L3"
	      << "_" << iProj+1 
	      << "_" << int(Age*100)
	      << "_" << DATRunEvent;
	
	// make projection:
	NBinsZ =hInter->GetZaxis ()->GetNbins ();
	int bin1 = 1 + NBinsZ*iProj/nProj;
	int bin2 = 0 + NBinsZ*(iProj+1)/nProj;
	cout << " b1/2: " << bin1 << " " << bin2 << endl;
	hInter->GetZaxis()->SetRange (bin1, bin2);
	//TH2D *hist =(TH2D*) hInter->Project3D (sProj.str ().c_str ());
	TH1D *hist =(TH1D*) hInter->Project3D (sProj.str ().c_str ());
	TH1D *hLem_proj = hInterLem->ProjectionX (sProj2.str ().c_str ());
	TH1D *hLmuon_proj = hInterLmuon->ProjectionX (sProj3.str ().c_str ());
	hist->Divide (Norm);
	hLem_proj->Divide (Norm);
	hLmuon_proj->Divide (Norm);
	
	
	if (Normalize) {
	    
	    if (hist->Integral ()) 
		hist->Scale (1./hist->Integral ());

	    if (hLmuon_proj->Integral ()) 
		hLmuon_proj->Scale (1./hLmuon_proj->Integral ());

	    if (hLem_proj->Integral ()) 
		hLem_proj->Scale (1./hLem_proj->Integral ());
	} 
	
	// set maxima/minima
	//hist->SetMaximum (maxT [iProj]);
	//hist->SetMinimum (minT [iProj]);
	hist->SetMaximum (.5);
	hist->SetMinimum (0);
	//hist->SetMaximum (1);
	//hist->SetMinimum (1e-4);
	
	
	// set axis titles
	hist->SetXTitle (histUnit.c_str ());
	//hist->SetYTitle ("log_{10}(E/GeV)");
	hist->SetYTitle ("dN");

	hLem_proj->SetXTitle (histUnit.c_str ());
	hLem_proj->SetYTitle ("1/norm * dN/dA");
	
	
	/*
 	hist->SetLineColor (color);
	hist->SetLineWidth (2);
	//hist->SetLineStyle (color);
	hist->SetDirectory (0);
	hist->Draw ((FirstDraw &&iProj==0) ? "" : "same");
	if (iProj==0)
	    legend->AddEntry (hist, sTitle.str ().c_str ());
	*/

	hLem_proj->SetLineColor (color);
	hLem_proj->SetMarkerColor (color);
	hLem_proj->SetMarkerStyle (marker);
	hLem_proj->SetMarkerSize (1.1);
	hLem_proj->SetLineWidth (2);
	hLem_proj->SetLineStyle (1);
	hLem_proj->SetDirectory (0);
	hLem_proj->GetXaxis ()->SetRangeUser (-2., 1.2);

	hLmuon_proj->SetLineColor (color);
	hLmuon_proj->SetMarkerColor (color);
	hLmuon_proj->SetMarkerStyle (marker);
	hLmuon_proj->SetMarkerSize (1.1);
	hLmuon_proj->SetLineWidth (2);
	hLmuon_proj->SetLineStyle (3);
	hLmuon_proj->SetDirectory (0);

	//hLem_proj->Draw ((FirstDraw &&iProj==0) ? "" : "same");
	hLem_proj->Draw ((FirstDraw &&iProj==0) ? "p" : "samep");
	hLem_proj->Fit (Para, "0Q1I");
	Para->SetLineColor (1);
	Para->SetLineWidth (1);
	Para->SetLineStyle (marker-19);
	
	if (iProj==0 && marker==20) {
	    legend->AddEntry (hLem_proj, sTitle.str ().c_str (), "p");
	    //legend->AddEntry (hist, sTitle.str ().c_str ());
	}

	if (PlotFit) {
	    TH1 *f = Para->GetHistogram ()->DrawCopy ("samel");
	    ostringstream sParaTitle;
	    sParaTitle << "NKG Funktion, Age=" << setprecision (2) << Age;
	    legend->AddEntry (f, sParaTitle.str ().c_str (), "l");
	}

	
	//hLmuon_proj->Draw ("same");
	//hLmuon_proj->Fit (Para, "0Q");
	//Para->GetHistogram ()->DrawCopy ("samel");
	
	
    }


    
	

    file->Close ();
}






TH1 *Interpolate3D (TH1* h1, double age1,
		     TH1* h2, double age2,
		    double age) {


    if (!((age<age1 && age>age2) ||
	  (age<age2 && age>age1))) {
	
	cerr << " CANNOT interpolate because age is outside [age1,age2]." 
	     << " age=" << age << " age1="<<age1 << " age2=" << age2 
	     << endl;
	return new TH3D ();	
    }

    ostringstream hName;
    hName << h1->GetName ()
	  << "_interpol";
    
    TH1 *h = (TH1*) h2->Clone (hName.str ().c_str ());
    
    h->Reset ("ICE");

    h->Add (h1, h2, std::fabs(age-age1), std::fabs(age-age2));

    return h;
}


Double_t NKG (Double_t *x, Double_t *p) {

    double R = pow (10, x[0]);//*p[1];
    double s = p [1];

    double rho = p [0];
    rho *= pow (R, s-2.);
    rho *= pow (1.+R, s-4.5);
    return rho;
}

/*
double SdSimpleSim::NKG (double N, double Rmoliere, double R, double s) {
    
    // normalize to Rmoliere
    R /= Rmoliere;

    double C = TMath::Gamma (4.5-s) /
	(TMath::Gamma (s) * TMath::Gamma (4.5-2*s));

    double rho = C * N/(2.*kPi*Rmoliere*Rmoliere);

    rho *= pow (R, s-2.);
    rho *= pow (1.+R, s-4.5);
    

    return rho;
}

*/
