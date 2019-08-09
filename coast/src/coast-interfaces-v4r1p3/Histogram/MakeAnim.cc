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
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
using namespace std;




void MakeAnim(const string &fileName);




int main(int argc, char **argv) {
  
  if (argc<2) {
    cerr << "      please specify data file (and cmd). " << endl;
    cerr << "       ./MakeAnim <file>" << endl; 
    return 1;
  }
  string fname(argv [1]);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  
  MakeAnim(fname);
  return 0;
}





void MakeAnim(const string& fileName) {
  
  TFile* file = TFile::Open (fileName.c_str (), "READ");
  
  if (file==NULL || file->IsZombie ()) {
    cerr << " Cannot open file " << fileName << endl;
    return;
  }
  
  file->cd ();
  TTree* run = (TTree*) file->Get("shower");
  TTree* tree = (TTree*) file->Get("data_layer");
  
  if (run == 0) {
    cout << "shower tree is missing" << endl;
    return;
  }
  if (tree == 0) {
    cout << "data_layer tree is missing" << endl;
    return;
  }
  
  
  Float_t Xmax, Zenith, Azimuth, Energy0;
  Double_t ObservationLevel, HeightFirstInt;
  Int_t EventNo, Primary;
  Char_t* particles = new Char_t[5000];
  Char_t* histograms = new Char_t[5000];

  run->SetBranchAddress("xmax"    , &Xmax  );
  run->SetBranchAddress("evntno"  , &EventNo);
  run->SetBranchAddress("obslevel", &ObservationLevel);
  run->SetBranchAddress("firstint", &HeightFirstInt );
  run->SetBranchAddress("zenith"  , &Zenith        );
  run->SetBranchAddress("azimuth" , &Azimuth        );
  run->SetBranchAddress("energy"  , &Energy0    );
  run->SetBranchAddress("primary" , &Primary   );
  run->SetBranchAddress("particles" , particles   );
  run->SetBranchAddress("histograms" , histograms   );
  run->GetEntry(0);
  
  
  vector<string> Particles;
  int start = 0;
  while (string(particles).find(' ', start) != string::npos) {
    const int pos = string(particles).find(" ", start);
    Particles.push_back(string(particles).substr(start, pos-start));
    start = pos+1;
  }
  if (string(particles).substr(start, string(particles).size()).size()>2) {
    Particles.push_back(string(particles).substr(start, string(particles).size()));
  }
  
  vector<string> Histograms;
  start = 0;
  while (string(histograms).find(' ', start) != string::npos) {
    const int pos = string(histograms).find(" ", start);
    Histograms.push_back(string(histograms).substr(start, pos-start));
    start = pos+1;
  }
  if (string(histograms).substr(start, string(histograms).size()).size()>2) {
    Histograms.push_back(string(histograms).substr(start, string(histograms).size()));
  }
  
  
  
  // Float_t time, depth, height;
  Float_t age;
  Float_t depth;
  //tree->SetBranchAddress ("time"  , &time );
  //tree->SetBranchAddress ("height", &height);
  tree->SetBranchAddress ("depth" , &depth );
  tree->SetBranchAddress ("age"   , &age );
  
  // create hist trees for each particle type
  map<string, TTree *> trees;
  for (vector<string>::iterator iParticle = Particles.begin();
       iParticle != Particles.end();
       ++iParticle) {
    
    ostringstream tree_name, tree_info;
    tree_name << "data_" << *iParticle;
    tree_info << "All the histogramms for " << *iParticle;
    
    TTree *newTree = (TTree*) file->Get(tree_name.str().c_str());
    if (newTree == 0) {
      cout << tree_name.str () << " tree is missing " << endl;
      return;
    }
    trees[*iParticle] = newTree;
  }
  
  
  
  
  const int Nhist = static_cast<int>(tree->GetEntries ());
  
  
  
  //ostringstream gif_dir;
  //gif_dir << fileName.substr(0, fileName.rfind (".root")); 
  
  const string DATRunEvent = fileName.substr(fileName.rfind ('/')+1,
					     fileName.rfind (".root") -
					     fileName.rfind ('/') - 1);
  
  const string gif_dir = DATRunEvent + "_anim";
  cout << " COAST: Make Directory: " << gif_dir << endl;
  gSystem->MakeDirectory(gif_dir.c_str());
  
  
  map<string, double> maxV, minV;
    
  map<string, map<string, TH1*> > obj_map;
  for (vector<string>::const_iterator iParticle = Particles.begin();
       iParticle != Particles.end();
       ++iParticle) {
    
    for (vector<string>::const_iterator iHist = Histograms.begin();
	 iHist != Histograms.end();
	 ++iHist) {
      
      ostringstream hnamess;
      hnamess << *iHist << "_" << *iParticle;
      const string hname = hnamess.str();
      
      obj_map[*iParticle][*iHist] = 0;	
      trees[*iParticle]->SetBranchAddress(hname.c_str(), &(obj_map[*iParticle][*iHist]));
      
      for (int iLayer=0; iLayer<Nhist; ++iLayer) {
	
	tree->GetEntry(iLayer);
	trees[*iParticle]->GetEntry(iLayer);	
	TH1* h = obj_map[*iParticle][*iHist];
	
	if (minV.count(*iHist)==0) {
	  minV[*iHist] = h->GetMinimum() ;
	  maxV[*iHist] = h->GetMaximum() ;
	} else {
	  if (minV[*iHist] > h->GetMinimum())
	    minV[*iHist] = h->GetMinimum();
	  if (maxV[*iHist] < h->GetMaximum()) 
	    maxV[*iHist] = h->GetMaximum();
	}
      }
    }
  } // end loop layers
  
  
    //bool first = true;
  for (vector<string>::const_iterator iParticle = Particles.begin();
       iParticle != Particles.end();
       ++iParticle) {
    
    for (vector<string>::const_iterator iHist = Histograms.begin();
	 iHist != Histograms.end();
	 ++iHist) {
      
      ostringstream hnamess;
      hnamess << *iHist << "_" << *iParticle;
      const string hname = hnamess.str();
      
      obj_map[*iParticle][*iHist] = 0;	
      trees[*iParticle]->SetBranchAddress(hname.c_str(), &(obj_map[*iParticle][*iHist]));
	
      for (int iLayer=0; iLayer<Nhist; ++iLayer) {
	//for (int iLayer=15; iLayer<18; ++iLayer) {
	
	tree->GetEntry(iLayer);
	trees[*iParticle]->GetEntry(iLayer);
	TH1* h = obj_map[*iParticle][*iHist];
	
	//const double delta = maxV[*iHist] - minV[*iHist];
	h->SetMaximum(maxV[*iHist]);
	h->SetMinimum(max(1.,minV[*iHist])); // for log scales
	
	TCanvas* cTmp = new TCanvas("cTmp", "one slice", 800, 600);
	if (string(h->ClassName())=="TH3D") {
	  cTmp->SetLogy(1);
	  ((TH3*)h)->Project3D("x")->Draw();
	} else if (string(h->ClassName())=="TH2D") {
	  cTmp->SetLogz(1);
	  //h->Draw("surf3");
	  h->Draw("cont0");
	} else if (string(h->ClassName())=="TH1D") {
	  cTmp->SetLogy(1);
	  h->Draw("l");
	} else if (string(h->ClassName())=="TProfile3D") {
	  cTmp->SetLogy(1);
	  ((TH3*)h)->Project3D("x")->Draw();
	} else if (string(h->ClassName())=="TProfile2D") {
	  cTmp->SetLogz(1);
	  h->Draw();
	} else if (string(h->ClassName())=="TProfile") {
	  cTmp->SetLogy(1);
	  h->Draw();
	} 
	
	// create and save gif
	ostringstream gif_name, eps_name;
	gif_name << gif_dir << "/" << *iHist << "_" << *iParticle << "_";
	gif_name.width(3);
	gif_name.fill('0');
	gif_name << right << iLayer << ".gif";
	
	eps_name << gif_dir << "/" << *iHist << "_" << *iParticle << "_";
	eps_name.width(3);
	eps_name.fill('0');
	eps_name << right << iLayer << ".eps";
	
	ostringstream eps2gif_cmd1;
	eps2gif_cmd1 << "pstopnm -ppm -xborder 0 -yborder 0 "
		     << "-xsize 900 -ysize 700 -portrait"
		     << " " << eps_name.str().c_str();
	
	ostringstream eps2gif_cmd2, cleanup;
	eps2gif_cmd2 << "ppmtogif " << eps_name.str().c_str() << "001.ppm > "
		     << gif_name.str();
	cleanup << "rm " << eps_name.str().c_str() << "001.ppm";
	
	cout << " COAST> save eps " << eps2gif_cmd1.str() 
	     << " " << eps2gif_cmd2.str() << endl;
	cout << " COAST> save gif " << eps2gif_cmd1.str() 
	     << " " << eps2gif_cmd2.str() << endl;
	cTmp->Print(eps_name.str().c_str());
	gSystem->Exec(eps2gif_cmd1.str().c_str());
	gSystem->Exec(eps2gif_cmd2.str().c_str());
	gSystem->Exec(cleanup.str().c_str());
	
	delete cTmp;
	
      }
    }
  }
  
  for (vector<string>::const_iterator iParticle = Particles.begin();
       iParticle != Particles.end();
       ++iParticle) {
    
    for (vector<string>::const_iterator iHist = Histograms.begin();
	 iHist != Histograms.end();
	 ++iHist) {
      
      ostringstream makeanim;
      makeanim << "gifsicle --delay=10 `ls "
	       << gif_dir << "/" << *iHist << "_" << *iParticle << "_" << "*.gif` > "
	       << DATRunEvent << "_" << *iHist << "_" << *iParticle << ".gif";      
      cout << makeanim.str() << endl;
      gSystem->Exec(makeanim.str ().c_str ());
    }
  }
  
  file->Close ();
  
}
