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

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH3D.h>

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
using namespace std;




void MakeAnim (const string &fileName,
	       const char *plot_cmd);




int main (int argc, char **argv) {

    if (argc<2) {
	cerr << "      please specify data file (and cmd). " << endl;
	cerr << "       ./MakeAnim <file> <cmd=\"x\"> " << endl; 
	return 1;
    }

    string fname (argv [1]);
    
    const char * plot_cmd = "x";
    if (argc>=3)
	plot_cmd = argv [2];


    MakeAnim (fname, plot_cmd);
    return 0;
}




void MakeAnim (const string &fileName,
	       const char *plot_cmd) {
    


    
    
    map <int,string> Particles;
    
    /* -------------------------------------------------------------------
      Define the particle names/types and CORSIKA ID here !!!!!!!!!!!!!!!!
    */
    const unsigned int nParticles = 2;
    //string particle_name [3] = {"EM", "Muon", "Hadron"};
    string particle_name [nParticles] = {"electron", "positron"};
    int particle_id      [nParticles] = {3, 2};
    // -------------------------------------------------------------------
    
    for (unsigned int j=0; j<nParticles; j++) {
	Particles [particle_id [j]] = particle_name [j];
    }
    
    
    
    
    
    
    
    TFile *file = TFile::Open (fileName.c_str (), "READ");

    if (file==NULL || file->IsZombie ())
	return;
    
    file->cd ();
    TTree *run = (TTree*) file->Get ("shower");
    TTree *tree = (TTree*) file->Get ("data_layer");

    if (run == 0) {
	cout << "run tree is missing" << endl;
	return;
    }
    if (tree == 0) {
	cout << "data_layer tree is missing" << endl;
	return;
    }
    

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
    
    
    
    
    
    
    // Float_t time, depth, height;
    Float_t age;
    /*
    tree->SetBranchAddress ("time"  , &time );
    tree->SetBranchAddress ("depth" , &depth );
    tree->SetBranchAddress ("height", &height);
    */
    tree->SetBranchAddress ("age"   , &age );


    // create hist trees for each particle type
    map<string,TTree *> trees;
    for (map<int,string>::iterator iParticle = Particles.begin ();
	 iParticle != Particles.end ();
	 ++iParticle) {
	
	ostringstream tree_name, tree_info;
	tree_name << "data_" << iParticle->second;
	tree_info << "All the histogramms for " << iParticle->second;
	
	TTree *newTree = (TTree*) file->Get (tree_name.str ().c_str ());
	
	if (newTree == 0) {
	    cout << tree_name.str () << " tree is missing " << endl;
	    return;
	}

	trees [iParticle->second] = newTree;
    }

 


    int Nhist = static_cast<int>(tree->GetEntries ());
  


    ostringstream gif_dir;
    gif_dir << fileName.substr (0, fileName.rfind (".root")); 

    string DATRunEvent = fileName.substr (fileName.rfind ('/')+1,
					  fileName.rfind (".root") -
					  fileName.rfind ('/') - 1);

    cout << " RU: Make Directory: " << gif_dir.str () << endl;
    gSystem->MakeDirectory (gif_dir.str ().c_str ());
    

    double maxT=0, maxA=0, minT=0, minA=0;
    bool firstM = true;
    for (int i=0; i<Nhist; ++i) {
	
	TCanvas *cTmp = new TCanvas ("cTmp", "one slice", 800, 600);
	cTmp->Divide (2);
	
	tree->GetEntry (i);

	for (map<int,string>::const_iterator iParticle = Particles.begin ();
	     iParticle != Particles.end ();
	     ++iParticle) {

	    TH3D *hT = 0;
	    TH3D *hA = 0;
	    
	    ostringstream htname;
	    htname << "hT" << iParticle->second;
	    ostringstream haname;
	    haname << "hA" << iParticle->second;

	    trees[iParticle->second]->SetBranchAddress (htname.str ().c_str (),
							&hT);
	    trees[iParticle->second]->SetBranchAddress (haname.str ().c_str (),
							&hA);
	    trees[iParticle->second]->GetEntry (i);
	    

	    double MaxT = hT->Project3D (plot_cmd)->GetMaximum ();
	    double MinT = hT->Project3D (plot_cmd)->GetMinimum ();
	    double MaxA = hA->Project3D (plot_cmd)->GetMaximum ();
	    double MinA = hA->Project3D (plot_cmd)->GetMinimum ();


	    if (firstM) {
		
		maxT = MaxT;
		minT = MinT;
		maxA = MaxA;
		minA = MinA;
		firstM = false;

	    } else {

		if (maxT<MaxT)
		    maxT = MaxT;
		
		if (minT>MinT)
		    minT = MinT;

		if (maxA<MaxA)
		    maxA = MaxA;
		
		if (minA>MinA)
		    minA = MinA;

	    }

	}

    }


    for (int i=0; i<Nhist; ++i) {
	

	TCanvas *cTmp = new TCanvas ("cTmp", "one slice", 800, 600);
	cTmp->Divide (2);
	
	tree->GetEntry (i);
	
	bool first = true;
	for (map<int,string>::const_iterator iParticle = Particles.begin ();
	     iParticle != Particles.end ();
	     ++iParticle) {
	    
	    TH3D *hT = 0;
	    TH3D *hA = 0;
	    
	    ostringstream htname;
	    htname << "hT" << iParticle->second;
	    ostringstream haname;
	    haname << "hA" << iParticle->second;
	    
	    trees[iParticle->second]->SetBranchAddress (htname.str ().c_str (),
							&hT);
	    trees[iParticle->second]->SetBranchAddress (haname.str ().c_str (),
							&hA);
	    trees[iParticle->second]->GetEntry (i);
	    

	    TH1* hTproj = (TH1*)  hT->Project3D (plot_cmd);
	    TH1* hAproj = (TH1*)  hA->Project3D (plot_cmd);

	    hTproj->SetMaximum (maxT*1.1);
	    hAproj->SetMaximum (maxA*1.1);
	    hTproj->SetLineColor (iParticle->first);
	    hAproj->SetLineColor (iParticle->first);
	    //hTproj->SetLineStyle (iParticle->first);
	    //hAproj->SetLineStyle (iParticle->first);
	    
	    cTmp->cd (1);
	    //gPad->SetLogy (1);
	    if (first) hTproj->Draw ();
	    else       hTproj->Draw ("same");

	
	    cTmp->cd (2);
	    //gPad->SetLogy (1);
	    if (first) hAproj->Draw ();
	    else       hAproj->Draw ("same");

	    first = false;
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
		     << "-xsize 900 -ysize 700 -portrait"
		     << " cTmp.eps";
	
	ostringstream eps2gif_cmd2;
	eps2gif_cmd2 << "ppmtogif cTmp.eps001.ppm > "
		     << gif_name.str ();
	
	
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
