#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TH1D.h>

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;




void MakeAnim (const string &fileName);




int main (int argc, char **argv) {

    if (argc<2) {
	cerr << "      please specify data file. " << endl;
	return 1;
    }

    string fname (argv [1]);
    MakeAnim (fname);
    return 0;
}




void MakeAnim (const string &fileName) {

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
    



    Float_t time, depth, height, age;

    tree->SetBranchAddress ("time"  , &time );
    tree->SetBranchAddress ("depth" , &depth );
    tree->SetBranchAddress ("height", &height);
    tree->SetBranchAddress ("age"   , &age );
 
    TH1D *hTem = 0;
    TH1D *hTmuon = 0;
    TH1D *hThadron = 0;
    TH1D *hEem = 0;
    TH1D *hEmuon = 0;
    TH1D *hEhadron = 0;
    TH1D *hLem = 0;
    TH1D *hLmuon = 0;
    TH1D *hLhadron = 0;
    TH1D *hAem = 0;
    TH1D *hAmuon = 0;
    TH1D *hAhadron = 0;

    // the histogramms
    tree->SetBranchAddress ("hTem"     , &hTem);
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

    int Nhist = tree->GetEntries ();
  




    ostringstream gif_dir;
    gif_dir << fileName.substr (0, fileName.rfind (".root")); 

    string DATRunEvent = fileName.substr (fileName.rfind ('/')+1,
					  fileName.rfind (".root") -
					  fileName.rfind ('/') - 1);

    cout << " RU: Make Directory: " << gif_dir.str () << endl;
    gSystem->MakeDirectory (gif_dir.str ().c_str ());
    


    for (int i=0; i<Nhist; ++i) {
	
	tree->GetEntry (i);


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
		     << "-xsize 900 -ysize 700 -portrait"
		     << " cTmp.eps";
	
	ostringstream eps2gif_cmd2;
	eps2gif_cmd2 << "ppmtogif cTmp.eps001.ppm > "
		     << gif_name.str ();
	
	TCanvas *cTmp = new TCanvas ("cTmp", "one slice", 800, 600);
	cTmp->Divide (2, 2);
	
	cTmp->cd (1);
	gPad->SetLogy (1);
	hTem->Draw ();
	hTmuon->Draw ("same");
	hThadron->Draw ("same");
	
	cTmp->cd (2);
	gPad->SetLogy (1);
	hEem->Draw ();
	hEmuon->Draw ("same");
	hEhadron->Draw ("same");
	
	cTmp->cd (3);
	gPad->SetLogy (1);
	hLem->Draw ();
	hLmuon->Draw ("same");
	hLhadron->Draw ("same");
	
	cTmp->cd (4);
	gPad->SetLogy (1);
	hAem->Draw ();
	hAmuon->Draw ("same");
	hAhadron->Draw ("same");
	
	cout << " RU: save gif " << eps2gif_cmd1.str () 
	     << " " << eps2gif_cmd2.str () << endl;
	cTmp->Print ("cTmp.eps");	
	gSystem->Exec (eps2gif_cmd1.str ().c_str ());
	gSystem->Exec (eps2gif_cmd2.str ().c_str ());
	
	delete cTmp;
    }

    
    ostringstream makeanim;
    makeanim << "gifsicle --delay=10 `ls "
	     << gif_dir.str () << "/h*.gif` > "
	     << DATRunEvent << ".gif";
    

    cout << makeanim.str () << endl;


    gSystem->Exec (makeanim.str ().c_str ());
    
    file->Close ();
    
}
