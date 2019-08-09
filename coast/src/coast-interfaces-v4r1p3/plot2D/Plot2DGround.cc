#include <crs/CorsikaTypes.h>
#include <crs/CorsikaConsts.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <crs2r/TC2R.h>
#include <crsRead/MCorsikaReader.h>

#include <TH2D.h>
#include <TLine.h>
#include <TFile.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>

#include <TCanvas.h>

#include <sstream>
#include <iostream>
using namespace std;
using namespace crs;

#define PHOTONS 1

#define ELECTRONS 2:\
case 3 

#define MUONS 5: \
case 6

#define PIZERO 7

#define CHARGEDPI 8: \
case 9

#define NEUTRON 13

#define PROTON 14:\
case 15 


const double GetCerenkovThrEnergy(int id){
  
  double mass=0; 
  switch (id){
    case MUONS: 
      mass=0.105;
      break; 
    case CHARGEDPI:
      mass= 0.140; 
      break; 
    default: 
      mass=0; 
  };

  double refIndex= 1.49; 
  double gamma= refIndex/sqrt(refIndex*refIndex-1.);
  
  double energy = mass* gamma;
  if (energy>0)
    //    cout<<"CherTher:"<<energy<<" "<<gamma<<endl;
  return energy;   
}



int main(int argc, char** argv) {

  if (argc<2) {
    cout << "  please specify Corsika file " << endl;
    return 0;
  }
    
  gROOT->ForceStyle("plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  

  string fname(argv[1]);
  crsRead::MCorsikaReader cr(fname, 3);
  
  string inputFile;
  if (fname.rfind('/')==string::npos) {
    inputFile = fname;
  } else {
    inputFile = fname.substr(fname.rfind('/')+1, 
                             fname.size()-fname.rfind('/'));
  }
  double xaxisup=  4; //log10(2.);
  double xaxislow= -4; //log10(0.1);
  int nbins = 50;
  TH2D histoEl("histoEl", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoEl.GetXaxis()->SetTitle("x (electrons) [km]");
  histoEl.GetYaxis()->SetTitle("y [km]");
  
  TH2D histoMu("histoMu", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoMu.GetXaxis()->SetTitle("x (muons) [km]");
  histoMu.GetYaxis()->SetTitle("y [km]");
 
  TH2D histoGa("histoGa", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoGa.GetXaxis()->SetTitle("x (photons) [km]");
  histoGa.GetYaxis()->SetTitle("y [km]");
 

  TH2D histoElT("histoElT", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoElT.GetXaxis()->SetTitle("x (electromagnetic) [km]");
  histoElT.GetYaxis()->SetTitle("y [km]");

  TH2D histoOt("histoOt", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoOt.GetXaxis()->SetTitle("x (others) [km]");
  histoOt.GetYaxis()->SetTitle("y [km]");

  TH2D histoPi("histoPi", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPi.GetXaxis()->SetTitle("x (pions) [km]");
  histoPi.GetYaxis()->SetTitle("y [km]");


  TH2D histoPiZero("histoPiZero", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPiZero.GetXaxis()->SetTitle("x (pions) [km]");
  histoPiZero.GetYaxis()->SetTitle("y [km]");



  TH2D histoPEl("histoPEl", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPEl.GetXaxis()->SetTitle("x [km]");
  histoPEl.GetYaxis()->SetTitle("y [km]");
  
  TH2D histoPMu("histoPMu", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPMu.GetXaxis()->SetTitle("x [km]");
  histoPMu.GetYaxis()->SetTitle("y [km]");
 
  TH2D histoPGa("histoPGa", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPGa.GetXaxis()->SetTitle("x [km]");
  histoPGa.GetYaxis()->SetTitle("y [km]");
 

  TH2D histoPElT("histoPElT", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPElT.GetXaxis()->SetTitle("x [km]");
  histoPElT.GetYaxis()->SetTitle("y [km]");

  TH2D histoPTot("histoPTot", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPTot.GetXaxis()->SetTitle("x [km]");
  histoPTot.GetYaxis()->SetTitle("y [km]");
  
  TH2D histoTot("histoTot", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoTot.GetXaxis()->SetTitle("x [km]");
  histoTot.GetYaxis()->SetTitle("y [km]");
 
  TH2D histoPOt("histoPOt", "", nbins, xaxislow, xaxisup, nbins, xaxislow, xaxisup);
  histoPOt.GetXaxis()->SetTitle("x [km]");
  histoPOt.GetYaxis()->SetTitle("y [km]");
  int nbinsldf=10;

  xaxislow=log10(0.200);
  TH1D ldfTot("ldfTot", "",nbinsldf , xaxislow, xaxisup);
  ldfTot.GetXaxis()->SetTitle("lg(r/km)");
  ldfTot.GetYaxis()->SetTitle("Particles ");

  TH1D ldfEl("ldfEl", "",nbinsldf , xaxislow, xaxisup);
  ldfEl.GetXaxis()->SetTitle("lg(r/km)");
  ldfEl.GetYaxis()->SetTitle("Particles ");
  TH1D ldfGa("ldfGa", "",nbinsldf , xaxislow, xaxisup);
  ldfGa.GetXaxis()->SetTitle("lg(r/km)");
  ldfGa.GetYaxis()->SetTitle("Particles ");

  TH1D ldfMu("ldfMu", "",nbinsldf , xaxislow, xaxisup);
  ldfMu.GetXaxis()->SetTitle("lg(r/km)");
  ldfMu.GetYaxis()->SetTitle("Particles ");


  TH1D ldfPi("ldfPi", "",nbinsldf , xaxislow, xaxisup);
  ldfPi.GetXaxis()->SetTitle("lg(r/km)");
  ldfPi.GetYaxis()->SetTitle("Particles ");


  TH1D pPi("pPi", "", 200 , -1, -1);
  pPi.GetXaxis()->SetTitle("p (GeV)");
  pPi.GetYaxis()->SetTitle("N ");

  TH1D pMu("pMu", "", 200 , -1, -1);
  pMu.GetXaxis()->SetTitle("p (GeV)");
  pMu.GetYaxis()->SetTitle("N ");
  
  TH1D pEl("pEl", "", 200 , -1, -1);
  pEl.GetXaxis()->SetTitle("p (GeV)");
  pEl.GetYaxis()->SetTitle("N ");

  TH1D ldfOt("ldfOt", "",nbinsldf ,0, log10(xaxisup));
  ldfOt.GetXaxis()->SetTitle("lg(r/km)");
  ldfOt.GetYaxis()->SetTitle("Particles ");
 
  int ShowerCounter = 0;
  int RunCounter = 0;
  
  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {
    
    RunCounter++;
    ShowerCounter = 0;
    
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      
      ShowerCounter ++;
      
      
      crs::TSubBlock Data;
      while (cr.GetData (Data)) {
        
        switch (Data.GetBlockType ()) {
	    

            case crs::TSubBlock::ePARTDATA:
            {
              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.FirstParticle();
                   iEntry != ParticleData.LastParticle();
                   ++iEntry) {
                
                if (iEntry->IsParticle()) {
                  
                  crs::MParticle iPart(*iEntry);
                   
                  int id    = iPart.GetParticleID();
                  
                  int level = iPart.GetObservationLevel();
                  double w  = iPart.GetWeight();
                  double e  = iPart.GetKinEnergy();
		  
		  double momentum= sqrt(iPart.GetPx()* iPart.GetPx()+
					iPart.GetPy()* iPart.GetPy()+  
					iPart.GetPz()* iPart.GetPz());
		  
		  //                  double x  = iPart.GetX()/10000.;
		  //                  double y  = iPart.GetY()/10000.;
		  
		  double x = iPart.GetX()/10000; //log10(fabs(iPart.GetX())/10000.);
		  double y = iPart.GetY()/10000; //log10(fabs(iPart.GetY())/10000.);
		  double cherthr= GetCerenkovThrEnergy(id); 
		  
// 		  if(cherthr >0) 
// 		    cout<<iPart.GetKinEnergy()<<" "<<GetCerenkovThrEnergy(id)<<endl; 

// 		  if (iPart.GetKinEnergy()< GetCerenkovThrEnergy(id) ){
// 		    continue;
// 		  } 
		  double r  = sqrt(iPart.GetX()*iPart.GetX()
				   +iPart.GetY()*iPart.GetY());
		  r/=10000.;
		  switch (id){

		    case ELECTRONS:
		      
		      histoEl.Fill(x, y, w);
		      histoPEl.Fill(x, y, w*iPart.GetKinEnergy());
		      ldfEl.Fill(log10(r), w);
		      pEl.Fill(momentum, w);

		      break; 
		      
		    case MUONS:
		      histoMu.Fill(x, y, w);
		      histoPMu.Fill(x, y, w*iPart.GetKinEnergy());
		      ldfMu.Fill(log10(r), w);
		      pMu.Fill(momentum, w);
		      break; 
		    case PHOTONS: 
		      histoGa.Fill(x, y, w);
		      histoPGa.Fill(x, y, w*iPart.GetKinEnergy());
		      ldfGa.Fill(log10(r), w);
		      break; 
		    case CHARGEDPI: 
		      histoPi.Fill(x, y, w);
		      ldfPi.Fill(log10(r), w);
		      pPi.Fill(momentum, w);
		      break; 

		    case PIZERO: 
		      //if(iPart.GetKinEnergy()> GetCerenkovThrEnergy(iPart.GetMass()))
		      histoPiZero.Fill(x, y, w);
		      break; 

		    default: break;  
		  };
		  if(id<65){
		    histoOt.Fill(x, y, w);
		    histoPOt.Fill(x, y, w*iPart.GetKinEnergy());
		    ldfOt.Fill(log10(r), w);
		  }
		  if(id<4){
		    histoElT.Fill(x, y, w);
		    histoPElT.Fill(x, y, w*iPart.GetKinEnergy());
		  }
		  
		  if(id<65){
		      ldfTot.Fill(log10(r), w);
		      histoTot.Fill(x, y, w);
		      histoPTot.Fill(x, y, w*iPart.GetKinEnergy());
		  }
                  // TODO: add particle to images

                }
                
              } // end particle loop
              
              break;
            }
            
            case crs::TSubBlock::eLONG:
              break;
              
            default:
              break;
        } // end data block
        
      } // loop data
      
      crs::MEventEnd ShowerSummary;
      crs::MEventHeader ShowerHeader;

      cr.GetShowerSummary(ShowerSummary);
      
      double Xmax = ShowerSummary.GetXmax();
      cout<<" "<<90.-Shower.GetTheta()*180./3.14<<" "<<
	  Shower.GetPhi()*180./3.14<<endl;

      double phi= Shower.GetTheta();
      cout << "---------------------------------\n"
           << " Shower info:\n"
           << "  Xmax = " << Xmax << "\n";
      
      TCanvas c2d("c2d","c2d",700,700);
      c2d.SetFillColor(0);
      c2d.SetLogz(1);
      c2d.SetTopMargin(0.15);
      c2d.SetRightMargin(0.15);
      c2d.SetBottomMargin(0.12);
      c2d.SetLeftMargin(0.12);
      histoTot.SetMinimum(100);
      histoTot.Draw("contz");
      c2d.Print("shower2d.eps");      
      
      TFile newfile("testfile.root", "RECREATE");
      TLine line(0.,0. ,-5.*cos(phi), 5*sin(phi));
      line.SetLineColor(kRed);
      
      TCanvas c;
      c.SetLogz();
      c.cd();
      c.SetFillColor(0);
      c.SetBorderMode(0);
      c.SetFrameBorderMode(0);
      c.SetCanvasSize(700, 700);
      histoElT.Draw("cont0");
      line.Draw("same");
      TCanvas c2;
      c2.SetFillColor(0);
      c2.SetBorderMode(0);
      c2.SetFrameBorderMode(0);
      c2.cd();
      histoMu.Draw("cont0z");
      line.Draw("same");
      c2.Write();
      c.Write();
      histoPi.Write();
      histoPiZero.Write();
      pPi.Write();
      pMu.Write();
      pEl.Write();

      histoEl.Write();
      histoElT.Write();
      histoMu.GetZaxis()->SetRangeUser(0, histoMu.GetMaximum());
      histoPi.GetZaxis()->SetRangeUser(0, histoMu.GetMaximum());
      histoMu.Write();
      histoGa.Write();
      histoOt.Write();
      histoTot.Write();
      histoPEl.Write();
      histoPElT.Write();
      histoPMu.Write();
      histoPGa.Write();
      histoPOt.Write();
      histoPTot.Write();
     
      ldfTot.Write();
      ldfEl.Write();
      ldfGa.Write();
      ldfMu.Write();      
      ldfPi.Write();
      ldfOt.Write();
      line.Write();

      newfile.Close();
    } // loop shower
    
  } // loop runs (usually just 1)    
  
  cout << " Read " << ShowerCounter << " showers from file " << endl;
  
  return 1;
}


extern "C" double heigh_(const double &f) {return 0;}




