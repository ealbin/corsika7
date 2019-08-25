
#ifndef tools_results_h
#define tools_results_h

#include <iostream>
#include <iomanip>
#include <math.h>

#include <TROOT.h>
#include <TObject.h>
#include <TAxis.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include "settings.h"
#include "run.h"
#include "sim.h"

class Results {

public:

    Run *fRun;
    Sim *fSim;
    Bool_t fAll_good; 

    Colors fColors;
    DensityOpts  fDensityOpts;
    SliceOpts    fSliceOpts;
    SpectrumOpts fSpectrumOpts;
    ContentOpts  fContentOpts;
    ImpactOpts   fImpactOpts;

    TList *fDensity;  // TH1D histograms for longitudinal density
    TList *fSlice;    // TH2D histograms for x-y slice
    TList *fSpectrum; // TH1D histograms for energy spectrum
    TH1D  *fContent;  // TH1D histogram  for particle content
    TH2D  *fImpact;   // TH2D histogram  for first impact

    void logBinning(TAxis *vAxis);
    void initializeHistograms();
    void printStdError();    

    // Assumes there is a file open, locates sim and run ttrees if they exist
    Results();

    // Uses supplied TChain or TTree
    //Results(TTree *vRun=NULL, TTree *vSim=NULL);

    ~Results() {
        cout << "delete!\n";
        if (gDirectory->IsOnHeap()) {
            if (fRun != NULL )
                delete fRun;
            if (fSim != NULL )
                delete fSim;

            if (fDensity != NULL) {
                fDensity->Delete();
                delete fDensity;
            }
            if (fSlice != NULL) {
                fSlice->Delete();
                delete fSlice;
            }
            if (fSpectrum != NULL) {
                fSpectrum->Delete();
                delete fSpectrum;
            }
            if (fContent != NULL)
                delete fContent;
            if (fImpact != NULL)
                delete fImpact;
        }
    }

    // Get run ttree if it was found
    Run* GetRun() { return fRun; }

    // Get sim ttree if it was found
    Sim* GetSim() { return fSim; }

    // Print observation level choices
    void ObservationLevels();

    // Plot longitudinal density for one or all observation levels (vObsLevel=-1)
    void PlotDensity(Int_t vObsLevel=-1);

};

#endif

