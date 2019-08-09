/***************************************************************************
 *   Copyright (C) 2006 by Tim Huege                                       *
 *   tim.huege@ik.fzk.de, s.lafebre@astro.ru.nl                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <TFile.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TLegend.h>

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace std;


void readLongFile(const string& fname, const int nProf,
		  TGraph& gEl, TGraph& gPo);

void 
histprofile(const string& filename, const string& histtype,
	    const bool hasLongFile, const string& longName)
{
  if ((histtype != "hA") && (histtype != "hT"))
    {
      cout << "\nError: specify either 'hA' for angle histograms or 'hT' for time histograms! Aborting ...\n\n";
      exit(10);
    }
  
  TGraph* gElLong = 0;
  TGraph* gPoLong = 0;  
  if (hasLongFile) {
    gElLong = new TGraph();
    gElLong->SetLineColor(kBlue);
    gElLong->SetLineWidth(3);
    gElLong->SetLineStyle(2);
    gPoLong = new TGraph();
    gPoLong->SetLineColor(kRed);
    gPoLong->SetLineWidth(3);
    gPoLong->SetLineStyle(2);
    readLongFile(longName, 0, *gElLong, *gPoLong);
  }
  
  TFile *file = TFile::Open (filename.c_str(), "READ");
    
  TTree *dataTree = (TTree*)file->Get("data_layer");
  TTree *electree = (TTree*)file->Get("data_electron");
  TTree *positree = (TTree*)file->Get("data_positron");

  const long layers = static_cast<long>(dataTree->GetEntries());

  Float_t rm, time, depth;
  dataTree->SetBranchAddress("rm", &rm);
  dataTree->SetBranchAddress("time", &time);
  dataTree->SetBranchAddress("depth", &depth);

  // map empty TH3D into TTree branch
  TH3D *pelechisto = 0;
  TH3D *pposihisto = 0;

  if (histtype == "hT")
    {
      electree->SetBranchAddress ("hTelectron", &pelechisto);
      positree->SetBranchAddress ("hTpositron", &pposihisto);
    }
  else
    {
      electree->SetBranchAddress ("hAelectron", &pelechisto);
      positree->SetBranchAddress ("hApositron", &pposihisto);
    }

  TGraph* gElHist = 0;
  TGraph* gPoHist = 0;  
  if (hasLongFile) {
    gElHist = new TGraph();
    gElHist->SetLineColor(kBlue);
    gElHist->SetLineWidth(3);
    gElHist->SetLineStyle(1);
    gPoHist = new TGraph();
    gPoHist->SetLineColor(kRed);
    gPoHist->SetLineWidth(3);
    gPoHist->SetLineStyle(1);
  }

  cout << "\n"
       << setw(11) << "#depth "
       << setw(12) << "electrons "
       << setw(12) << "positrons "
       << setw(11) << "%ufx "
       << setw(11) << "%ofx "
       << setw(11) << "%ufy "
       << setw(11) << "%ofy "
       << setw(11) << "%ufz "
       << setw(11) << "%ofz \n\n";

  double elecsum = 0.0, posisum = 0.0;
  for (int i = 0; i<layers; ++i)
    {
      dataTree->GetEntry (i);
      electree->GetEntry (i);
      positree->GetEntry (i);
      double elecintegral=0, posiintegral=0, ufx=0, ofx=0, ufy=0, ofy=0, ufz=0, ofz=0;
      for (long ii=0; ii<=pelechisto->GetNbinsX()+1; ++ii)
	{
	  for (long j=0; j<=pelechisto->GetNbinsY()+1; ++j)
	    {
	      for (long k=0; k<=pelechisto->GetNbinsZ()+1; ++k)
		{
		  double elecs = pelechisto->GetBinContent(pelechisto->GetBin(ii, j, k));
		  double posis = pposihisto->GetBinContent(pposihisto->GetBin(ii, j, k));
		  if (ii==0) { ufx += elecs+posis; }
		  if (ii==pelechisto->GetNbinsX()+1) { ofx += elecs+posis; }
		  if (j==0) { ufy += elecs+posis; }
		  if (j==pelechisto->GetNbinsY()+1) { ofy += elecs+posis; }
		  if (k==0) { ufz += elecs+posis; }
		  if (k==pelechisto->GetNbinsZ()+1) { ofz += elecs+posis; }
		  elecintegral += elecs;
		  posiintegral += posis;
		}
	    }
	}
      elecsum += elecintegral;
      posisum += posiintegral;
      const double total = elecintegral + posiintegral;
      if (total!=0) {
	ufx /= total;
	ofx /= total;
	ufy /= total;
	ofy /= total;
	ufz /= total;
	ofz /= total;
	cout << setw(10) << setprecision(4) << fixed << depth << " "
	     << setw(11) << (elecintegral<1e5 ? fixed : scientific) 
	     << elecintegral << " "
	     << setw(11) << (posiintegral<1e5 ? fixed : scientific) 
	     << posiintegral  << " "
	     << scientific << setprecision(3) << setw(10) << ufx << " "
	     << setw(10) << ofx << " "
	     << setw(10) << ufy << " "
	     << setw(10) << ofy  << " "
	     << setw(10) << ufz  << " "
	     << setw(10) << ofz << fixed << '\n';
      } else {
	cout << setw(10) << setprecision(4) << fixed << depth << " "
	     << setw(11) << (elecintegral<1e5 ? fixed : scientific) 
	     << elecintegral << " "
	     << setw(11) <<  (posiintegral<1e5 ? fixed : scientific) 
	     << posiintegral  << " "
	     << scientific << setprecision(3) << setw(10) << "- "
	     << setw(10) << "- "
	     << setw(10) << "- "
	     << setw(10) << "- " 
	     << setw(10) << "- " 
	     << setw(10) << "- " << '\n';
	
      }

      if (hasLongFile) {
	gElHist->SetPoint(i, depth, elecintegral);
	gPoHist->SetPoint(i, depth, posiintegral);
      }

    } // loop layers
  
  cout << "\n#elecsum: " << elecsum 
       << "\tposisum: " << posisum 
       << "\telec/posi: " << elecsum/posisum << "\n";
  
  if (hasLongFile) {
    TCanvas* canvas = new TCanvas("profiles compare");
    TLegend* legend = new TLegend(0.15, 0.7, 0.475, 0.89);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    gElLong->Draw("al");
    gPoLong->Draw("l");
    gElHist->Draw("l");
    gPoHist->Draw("l");
    legend->AddEntry(gElLong, "electrons longfile", "l");
    legend->AddEntry(gPoLong, "positrons longfile", "l");
    legend->AddEntry(gElHist, "electrons histfile", "l");
    legend->AddEntry(gPoHist, "positrons histfile", "l");
    legend->Draw();
    canvas->SaveAs("histprofile.eps");
  }
  
}


int 
main (int argc, char** argv) 
{
  if (argc < 3) {
    cout << "\n Tool to compute the longitudinal profile for electron/positrons from hist files\n"
	 << " Usage: histprofile <histfile> <hT/hA> [long-file]\n" << endl;
    return 0;
  }

  const bool hasLongFile = argc >= 4;
  const string longFile = hasLongFile ? argv[3] : "none";
  
  gROOT->SetStyle("Plain");

  histprofile(argv[1], argv[2], hasLongFile, longFile);
  return 0;
}

// --------------- HELPER FUNCTIONS ----------------------
void 
readLongFile(const string& fname, const int nProf,
	     TGraph& gEl, TGraph& gPo)
{
  ifstream longFile(fname.c_str());
  
  int n = 0;
  int iProf = 0;
  while(longFile.good()) {
    char linech[1000];
    longFile.getline(linech, 1000);
    string line(linech);
    cout << " ======> " << line << endl;
    if (line.find("LONGITUDINAL DISTRIBUTION IN") == 1) {
      string dummy;
      istringstream ss(line);
      ss >> dummy >> dummy >> dummy >> n;
      cout << " found prof with length " << n << endl;
      if (iProf == nProf)
	break;
      ++iProf;
      break;
    }
    if (!longFile.good()) {
      cout << "Corsika long file has bad format. Parsing stopped." << endl;
      longFile.close();
      return;
    }
  }
  longFile.ignore(99999, '\n');
  gEl.Set(n);
  gPo.Set(n);
  for (int i=0; i<n; ++i) {
    double depth, nel, npo, dummy;
    longFile >> depth >> dummy >> npo >> nel;
    longFile.ignore(99999, '\n');  
    gEl.SetPoint(i, depth, nel);
    gPo.SetPoint(i, depth, npo);
  }  
  
  longFile.close();
}
