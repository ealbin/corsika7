#include <TPlotter.h>

#include <crs/CorsikaConsts.h>
using namespace crs;

#include <TProfile.h>
#include <TProfile2D.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>

#include <iostream>
using namespace std;

/** ***********************************************************************
  User function to define particle type that should be histogrammed
**/
void
TPlotter::InitParticles() 
{  
  /*
    Syntax:  fParticles[<corsika-particle-id>] = ParticleDef("<name>", <color-code>);    
    Warning: <corsika-particle-id> needs to be unique
  */
  
  //fParticles[2] = ParticleDef("positron", kRed);
  fParticles[3] = ParticleDef("electron", kBlue);
  //fParticles[5] = ParticleDef("mu+", kGreen+1);
  fParticles[6] = ParticleDef("mu-", kGreen+1);  
}




/** ************************************************************************
  User function to define histograms.
  These histograms are filled for each depth layer for all particle types
    -> Be carfull with the dimensions ...
**/
void
TPlotter::InitHistograms(HistDef& hists) 
{  
  /*
    Syntax: hists[<hist-id-string>] = <historam>
    Info:   All histograms are automatically generated for all particle
            types and depth layers.
            The histogram name is automatically appended with appropriate 
            additions.
   */
  
  // --------------------------------------------------
  // generate histograms

  /*  
  TH2D* hist1 = new TH2D("hPlane", "plane",
			 50, -1, 1,  // [km]
			 50, -1, 1); // [km]
  hist1->SetXTitle("X   [km]");
  hist1->SetYTitle("Y   [km]");  
  hists["1"] = hist1; 
  */
  TProfile* hist2 = new TProfile("hAngle", "angle",
				 6, -2.5, 2.5, "s");  // [log(r/rm)]
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerSize(0.8);
  hist2->SetLineWidth(2);
  hist2->SetXTitle("log_{10}(r/r_{m})");
  hist2->SetYTitle("Theta   [deg]");
  hists["2"] = hist2; 
}




/** ************************************************************************
  User function to fill histograms
 */
void
TPlotter::FillHistograms(HistDef& hists,
                         const double& /*x*/, const double& /*y*/, const double& /*time*/, 
                         const double& r, const double& rm,
                         const double& /*energy*/, const double& theta, const double& /*phi*/, 
                         const double& weight) 
{  
  // ((TH2D*)hists["1"])->Fill(x/km, y/km, weight);
  ((TProfile*)hists["2"])->Fill(log10(r/rm), theta/deg, weight);
}      
