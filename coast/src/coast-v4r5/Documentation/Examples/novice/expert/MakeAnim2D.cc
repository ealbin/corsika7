#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TStyle.h>

#include <iostream>
#include <string>
#include <sstream>
using namespace std;


#include "extractTH2D.C"


void MakeAnim (const string &fileName, bool Normalize);




int main (int argc, char **argv) {

    if (argc<2) {
	cerr << "      please specify data file. " << endl;
	return 1;
    }

    bool Normalize = false; 
    if (argc==3) {
	
	istringstream sMode (argv [2]);
	sMode >> Normalize;	
    }

    if (Normalize) {
	cout << " --> Normalize all spectra " << endl;
    } else {
	cout << " --> Do not normalize spectra " << endl;
    }

    string fname (argv [1]);
    MakeAnim (fname, Normalize);
    return 0;
}




void MakeAnim (const string &fileName, bool Normalize) {
    
    gStyle->SetPalette (1);
    gStyle->SetPadRightMargin (.15);
    gStyle->SetOptStat (0);

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
    TH2D *hLem = 0;
    TH2D *hLmuon = 0;
    TH2D *hLhadron = 0;
    TH3D *hAem = 0;
    TH3D *hAmuon = 0;
    TH3D *hAhadron = 0;
    */
    
    string histUnit = "log_{10}(t/ns)";
    //string histUnit = "theta/deg";
    
    // the histogramms
    tree->SetBranchAddress ("hTem"     , &hist3D);
    //tree->SetBranchAddress ("hAem"     , &hist3D);

    /*
    tree->SetBranchAddress ("hTmuon"   , &hTmuon);
    tree->SetBranchAddress ("hThadron" , &hThadron);
    tree->SetBranchAddress ("hEem"     , &hEem);
    tree->SetBranchAddress ("hEmuon"   , &hEmuon);
    tree->SetBranchAddress ("hEhadron" , &hEhadron);
    tree->SetBranchAddress ("hLem"     , &hLem);
    tree->SetBranchAddress ("hLmuon"   , &hLmuon);
    tree->SetBranchAddress ("hLhadron" , &hLhadron);
    tree->SetBranchAddress ("hAem"     , &hAem);
    tree->SetBranchAddress ("hAmuon"   , &hAmuon);
    tree->SetBranchAddress ("hAhadron" , &hAhadron);
    */

    int Nhist = tree->GetEntries ();
  




    ostringstream gif_dir;
    gif_dir << fileName.substr (0, fileName.rfind (".root")); 

    string DATRunEvent = fileName.substr (fileName.rfind ('/')+1,
					  fileName.rfind (".root") -
					  fileName.rfind ('/') - 1);

    cout << " RU: Make Directory: " << gif_dir.str () << endl;
    gSystem->MakeDirectory (gif_dir.str ().c_str ());
    
    ostringstream tmp_name;
    tmp_name << gif_dir.str () << "/cTmp.eps";       

    int NBinsZ = 0;


    static const int nProj = 4;

    double maxT[nProj] = {0,0,0,0};
    double minT[nProj] = {1e10,1e10,1e10,1e10};

    if (Normalize) {
	
	for (int i=0; i<Nhist; ++i) {

	    tree->GetEntry (i);

	    for (int iProj=0; iProj<nProj; iProj++) {
		    
		maxT [iProj] = 1;


		// make projection:
		NBinsZ = hist3D->GetYaxis ()->GetNbins ();
		int bin1 = 0 + NBinsZ*iProj/nProj;
		int bin2 = 1 + NBinsZ*(iProj+1)/nProj;
		//hist3D->GetYaxis()->SetRange (bin1, bin2);
		//TH2D *hist = (TH2D*) hist3D->Project3D ("zx");
		TH2D *hist = extractY (hist3D, "zx", bin1, 0);

		double iT = hist->Integral ();
		if (iT) {
		    double histT_min = hist->GetMinimum (.9)/iT;
		    double histT_max = hist->GetMaximum ()/iT;
		    if (histT_min<minT [iProj])
			minT [iProj] = histT_min;
		    if (histT_max>maxT [iProj])
			maxT [iProj] = histT_max;
		}
		
		delete hist;
	    }
	}

    } else {

	
	for (int i=0; i<Nhist; ++i) {

	    tree->GetEntry (i);
	
	    for (int iProj=0; iProj<nProj; iProj++) {
		
		minT [iProj] = 1;

		// make projection:
		NBinsZ = hist3D->GetYaxis ()->GetNbins ();
		int bin1 = 0 + NBinsZ*iProj/nProj;
		int bin2 = 1 + NBinsZ*(iProj+1)/nProj;
		//hist3D->GetYaxis()->SetRange (bin1, bin2);
		//TH2D *hist = (TH2D*) hist3D->Project3D ("zx");
		TH2D *hist = extractY (hist3D, "zx", bin1, 0);

				
		double histT_max = hist->GetMaximum ();
		if (histT_max>maxT [iProj])
		    maxT [iProj] = histT_max;
		
		delete hist;
	    }	    

	}
	
	for (int iProj=0; iProj<nProj; iProj++)
	    maxT [iProj] *= 1.1;	   

    }



    int delay = 100/Nhist * 15; // total time: 15*(s/100)*100ms = 15s

	
    for (int i=0; i<Nhist; ++i) {
	
	tree->GetEntry (i);


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
	       << ", E_{0}=" << Energy0 << "GeV "
	       << ", #theta= " << Zenith/3.141*180. << "deg "
	       << ", X_{max}=" << Xmax << "g/cm^{2} ";
	
	
	TCanvas *cCanvas = new TCanvas ("cCanvas", "one slice", 800, 700);
	TPaveLabel *Title = new TPaveLabel (.0, .95, 1.0, 1.0, 
					    sTitle.str ().c_str ());
	Title->SetFillColor (188);
	Title->Draw ();
	TPad *cTmp = new TPad ("cTmp", "cTmp", .0, .0, 1.0, .95, -1, 0, 0);
	cTmp->Draw ();
	cTmp->Divide (2, 2);
	
	
	for (int iProj=0; iProj<nProj; iProj++) {
	    
	    ostringstream sProj;
	    sProj << "zx"
		  << "_" << i+1 << "_" << iProj+1;

	    // make projection:
	    NBinsZ = hist3D->GetYaxis ()->GetNbins ();
	    int bin1 = 0 + NBinsZ*iProj/nProj;
	    int bin2 = 1 + NBinsZ*(iProj+1)/nProj;
	    cout << bin1 << " " << bin2 << endl;
	    //hist3D->GetYaxis()->SetRange (bin1, bin2);
	    //TH2D *hist = (TH2D*) hist3D->Project3D (sProj.str ().c_str ());
	    TH2D *hist = extractY (hist3D, sProj.str ().c_str (),bin1, 0);
	    
	    if (Normalize) {
		
		if (hist->Integral ()) 
		    hist->Scale (1./hist->Integral ());
	    } 

	    // set maxima/minima
	    hist->SetMaximum (maxT [iProj]);
	    hist->SetMinimum (minT [iProj]);
	    

	    // set axis titles
	    hist->SetXTitle (histUnit.c_str ());
	    hist->SetYTitle ("log_{10}(E/GeV)");
	    hist->SetZTitle ("dN");
	    
	    
	    cTmp->cd (iProj+1);
	    gPad->SetLogz (1);
	    hist->Draw ("colz");
	    //gPad->SetLogy (1);
	    //hist->Draw ("");
	    hist->Dump ();
	    
    
	    cTmp->cd ();
	    ostringstream sAge;
	    sAge << " age = "; 
	    sAge.width (6);
	    sAge.fill('0');
	    sAge.precision (4);
	    sAge << left << age;
	    TPaveLabel *AgeText = new TPaveLabel (.4, .473, .6, .527, 
						  sAge.str ().c_str (), "");
	    AgeText->SetFillColor (70);
	    AgeText->Draw ();
	    
	}
	
	// create and save gif
	ostringstream gif_name, eps_name;
	gif_name << gif_dir.str () << "/h_";
	gif_name.width (3);
	gif_name.fill ('0');
	gif_name << right << i << ".gif";
	
	eps_name << gif_dir.str () << "/h_";
	eps_name.width (3);
	eps_name.fill ('0');
	eps_name << right << i << ".eps";
	
	ostringstream eps2gif_cmd1;
	eps2gif_cmd1 << "pstopnm -ppm -xborder 0 -yborder 0 "
		     << "-xsize 900 -ysize 700 -portrait -stdout "
		     << tmp_name.str ()
		     << " > " << gif_dir.str () << "/cTmp.ppm";
	
	ostringstream eps2gif_cmd2;
	eps2gif_cmd2 << "ppmtogif "
		     << gif_dir.str () << "/cTmp.ppm "
		     << " > " << gif_name.str ();
	
	

	cout << " RU: save gif " << eps2gif_cmd1.str () 
	     << " " << eps2gif_cmd2.str () << endl;
	cCanvas->Print (tmp_name.str ().c_str ());
	gSystem->Exec (eps2gif_cmd1.str ().c_str ());
	gSystem->Exec (eps2gif_cmd2.str ().c_str ());
	
	delete cCanvas;
    }

    
    ostringstream makeanim;
    makeanim << "gifsicle --delay=" << delay << " `ls "
	     << gif_dir.str () << "/h*.gif` > "
	     << DATRunEvent << "_2D.gif";
    

    ostringstream rm_dir;
    rm_dir << " rm -rf " << gif_dir.str ().c_str ();


    cout << makeanim.str () << endl;


    gSystem->Exec (makeanim.str ().c_str ());
    gSystem->Exec (rm_dir.str ().c_str ());

    file->Close ();
    
}
