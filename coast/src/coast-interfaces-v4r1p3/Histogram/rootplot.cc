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
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraphErrors.h>

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

void DrawFrame(const TH1*, const TGraph*);
void rootplot(string filename, string slicenumber, string plottype, TCanvas& c1);

void rootplot(string filename, string slicenumber, string plottype, TCanvas& c1)
{

    long layer;
    stringstream ss;
    ss << slicenumber;
    ss >> layer;
	
    TFile *file = TFile::Open (filename.c_str(), "READ");

    TTree *dataTree = (TTree*)file->Get("data_layer");
    TTree *electree = (TTree*)file->Get("data_electron");
    TTree *positree = (TTree*)file->Get("data_positron");

    Float_t rm, time, depth;
    dataTree->SetBranchAddress("rm", &rm);
    dataTree->SetBranchAddress("time", &time);
    dataTree->SetBranchAddress("depth", &depth);

    // map empty TH3D into TTree branch
    TH3D *histAelectron = 0;
    TH3D *histTelectron = 0;
    TH3D *histApositron = 0;
    TH3D *histTpositron = 0;
    TGraphErrors *histOffsetXelectron = 0;
    TGraphErrors *histOffsetYelectron = 0;
    TGraphErrors *histOffsetXpositron = 0;
    TGraphErrors *histOffsetYpositron = 0;
    electree->SetBranchAddress ("hTelectron", &histTelectron);
    positree->SetBranchAddress ("hTpositron", &histTpositron);
    electree->SetBranchAddress ("hAelectron", &histAelectron);
    positree->SetBranchAddress ("hApositron", &histApositron);
    electree->SetBranchAddress ("offsetX_electron", &histOffsetXelectron);
    positree->SetBranchAddress ("offsetX_positron", &histOffsetXpositron);
    electree->SetBranchAddress ("offsetY_electron", &histOffsetYelectron);
    positree->SetBranchAddress ("offsetY_positron", &histOffsetYpositron);

    dataTree->GetEntry (layer);
    electree->GetEntry (layer);
    positree->GetEntry (layer);

      // now create plots

    TH3D* histT = histTelectron;
    TH3D* histA = histAelectron;
    TGraphErrors* histOffsetX = histOffsetXelectron;
    TGraphErrors* histOffsetY = histOffsetYelectron;
    if (plottype == "posi")
    { 
      histT = histTpositron; 
      histA = histApositron; 
      histOffsetX = histOffsetXpositron;
      histOffsetY = histOffsetYpositron;
    }
    
    if ((histOffsetX != 0) && (histOffsetY !=0))
    {
      histOffsetX->SetTitle("Offset in north direction");
      histOffsetY->SetTitle("Offset in west direction");
      histOffsetX->GetHistogram()->SetXTitle("log_{10}(E/eV)");
      histOffsetY->GetHistogram()->SetXTitle("log_{10}(E/eV)");
      histOffsetX->GetHistogram()->SetYTitle("Offset [cm]");
      histOffsetY->GetHistogram()->SetYTitle("Offset [cm]");
      histOffsetX->SetMarkerStyle(25);
      histOffsetY->SetMarkerStyle(25);
      histOffsetX->SetMarkerSize(0.5);
      histOffsetY->SetMarkerSize(0.5);
      histOffsetX->SetMarkerColor(histT->GetLineColor());
      histOffsetY->SetMarkerColor(histT->GetLineColor());
      histOffsetX->SetLineColor(histT->GetLineColor());
      histOffsetY->SetLineColor(histT->GetLineColor());
    }

    c1.cd(1);
    histT->Project3D(("x"+filename).c_str())->Draw("same");
    c1.cd(3);
    histT->Project3D(("y"+filename).c_str())->Draw("same");
    c1.cd(5);
    histT->Project3D(("z"+filename).c_str())->Draw("same");
    c1.cd(2);
    histA->Project3D(("x"+filename).c_str())->Draw("same");
    c1.cd(4);
    histA->Project3D(("y"+filename).c_str())->Draw("same");
    c1.cd(6);
    histA->Project3D(("z"+filename).c_str())->Draw("same");
    if ((histOffsetX != 0) && (histOffsetY !=0))
    {
      c1.cd(7);
      DrawFrame(histT->Project3D("z"), histOffsetX);
      histOffsetX->Draw("lp,same");
      c1.cd(8);
      DrawFrame(histT->Project3D("z"), histOffsetY);
      histOffsetY->Draw("lp,same");
    }
    c1.Update();

    delete histAelectron;
    delete histApositron;
    delete histTelectron;
    delete histTpositron;
    
}

int main (int argc, char** argv) {

	if ((argc != 4) && (argc != 7)) {
		cout << "\nUsage: rootplot <histfile1> <sliceno> <elec/posi> [<histfile2> <sliceno> <elec/posi>]\n" 
                     << endl;
		return 0;
	}

	TCanvas c1("c1",argv[1],800,600);
	c1.SetFillColor(kWhite);
	gStyle->SetOptStat(0);
	gStyle->SetMarkerSize(0.5);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(2);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetCanvasBorderSize(2);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetPadBorderSize(2);
	gStyle->SetPadBorderMode(0);
	c1.Divide(2,4);

	rootplot(argv[1], argv[2], argv[3], c1);
	if (argc == 7) { rootplot(argv[4], argv[5], argv[6], c1); }

	c1.Print("rootplot.eps");

	return 0;
}


void DrawFrame(const TH1* h, const TGraph* g) {

  double xMin = h->GetXaxis()->GetXmin();
  double xMax = h->GetXaxis()->GetXmax();
  double yMin = 0;
  double yMax = 0;
  bool first = true;
  for (int i=0; i<g->GetN(); ++i) {
    if (first) {
      first = false;
      yMin = g->GetY()[i];
      yMax = g->GetY()[i];
    } else {
      if (yMin>g->GetY()[i]) yMin = g->GetY()[i];
      if (yMax<g->GetY()[i]) yMax = g->GetY()[i];
    }
  }
  
  yMin -= (yMax-yMin)*0.15;
  yMax += (yMax-yMin)*0.15;
  yMax *= 1.5;
  
  ostringstream frameName;
  frameName << "frame_" << g->GetTitle() << "_" << h->GetName();
  
//  TH2F *frame = new TH2F(frameName.str().c_str(), g->GetTitle(), 2, xMin, xMax, 2, yMin, yMax);  
  TH2F *frame = new TH2F(frameName.str().c_str(), g->GetTitle(), 2, xMin, xMax, 2, -yMax, yMax);  
  frame->SetXTitle(g->GetHistogram()->GetXaxis()->GetTitle());
  frame->SetYTitle(g->GetHistogram()->GetYaxis()->GetTitle());
  frame->Draw("same");
}
