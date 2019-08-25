
#include "results.h"

// start root by loading the data file:
// $ root -l DATxxxxxx.root  #(aka, assumes your data is _file0)
// [] .L results.cc
// [] Results r;
// [] r.blablabla

Results::Results(TTree *vRun, TTree *vSim) {
   
    fRun = NULL;
    fSim = NULL; 
    
    fDensity  = NULL;
    fSlice    = NULL;
    fSpectrum = NULL;
    fContent  = NULL;
    fImpact   = NULL;

    if (vRun == NULL || vSim == NULL ||
        TString(vRun->GetName()) != "run" || 
        TString(vSim->GetName()) != "sim") {
        cout << " [ERROR] Input is invalid (null pointers or wrong TTrees)\n";
        return;
    }

    fRun = new Run(vRun);
    fSim = new Sim(vSim);
    fRun->GetEntry(0);
    fSim->GetEntry(0);
    fAll_good = true;
    cout << " [info] saul goodman\n";
    initializeHistograms();
}

// Locates sim and run ttrees if they exist
//------------------------------------------------------------------------
Results::Results() {

    fRun = NULL;
    fSim = NULL; 
    
    fDensity  = NULL;
    fSlice    = NULL;
    fSpectrum = NULL;
    fContent  = NULL;
    fImpact   = NULL;

    TString *filename = new TString( gDirectory->GetName() );
    Bool_t is_data_file = filename->EndsWith(".root");
    
    fAll_good = false;
    Bool_t has_run_tree = false;
    Bool_t has_sim_tree = false;

    if (is_data_file) {
        int n_keys = gDirectory->GetNkeys();
        TList *key_list = gDirectory->GetListOfKeys();

        TString *key_name; 
        for (int i = 0; i < n_keys; i++) {
            key_name = new TString(key_list->At(i)->GetName());

            if (*key_name == "run") {
                fRun = new Run( (TTree*)gDirectory->Get("run") );
                fRun->GetEntry(0);
                has_run_tree = true;
                cout << " [info] found run tree\n";
            }
            else if (*key_name == "sim") {
                fSim = new Sim( (TTree*)gDirectory->Get("sim") );
                fSim->GetEntry(0);
                has_sim_tree = true;
                cout << " [info] found sim tree\n";
            } 
        
            delete key_name;
        }
    }
    
    delete filename;

    fAll_good = has_run_tree && has_sim_tree;

    if (fAll_good) {
        cout << " [info] saul goodman\n";
        initializeHistograms();
    }
    else {
        cout << " [ERROR] snot goodman\n";
        if (!is_data_file)
            cout << "         - no root file!\n";
        else {
            if (!has_run_tree)
                cout << "         - no run tree!\n";
            if (!has_sim_tree)
                cout << "         - no sim tree!\n";
        }
    }
}


//////////////////////////////////////////////////////////////////////////////

void Results::ObservationLevels() {
    if (!fAll_good) {
        printStdError();
        return;
    }

    cout << endl;
    cout << "index     alititude [km]\n";
    cout << "------------------------\n";
    
    Int_t nlevels = fRun->run_ObservationLevel.size();
    for (Int_t i = 0; i < nlevels && i < 9; i++)
        cout << setw(3)  << (i+1) << "      " 
             << setw(7) << fRun->run_ObservationLevel.at(i)/1e5 << endl;
    if (nlevels == 10)
        cout << setw(3)  << '0' << "      " 
             << setw(7) << fRun->run_ObservationLevel.at(9)/1e5 << endl;
    cout << endl;
}


//////////////////////////////////////////////////////////////////////////////

void Results::PlotDensity(Int_t vObsLevel) {
    if (!fAll_good) {
        printStdError();
        return;
    }
    
    TString name = fDensityOpts.name;
    TString expression = "sqrt(particle.x**2 + particle.y**2)";
    
    TString title;
    TString selection = "";
    Int_t end   = fRun->run_ObservationLevel.size();
    if (vObsLevel >= 0 && vObsLevel < end) {
        title = Form("Observation Level = %f [km]", 
                     fRun->run_ObservationLevel.at(vObsLevel)/1e5);
        selection = Form("(particle.ObservationLevel == %d) && ", vObsLevel);
    }
    else if (vObsLevel == -1)
        title = "Composit of All Observation Levels";
    else {
        cout << " [ERROR] Invalid observation level index\n";
        ObservationLevels();
        return;
    }

    Int_t entries = fDensity->GetEntries();
    for (Int_t i = 0; i < entries; i++) {
        ((TH1D*)fDensity->At(i))->SetTitle(title);
        ((TH1D*)fDensity->At(i))->Reset("ICEM");
    }

    //TTree *sim = fSim->fChain;
    TString id = "particle.ParticleID == ";
    TString select;

    select = Form("%s(%s%d)", selection.Data(), id.Data(), 1);
    fSim->fChain->Project(name + "_photon", expression, select);
/*
    select = Form("%s(%s%d || %s%d)", selection.Data(), id.Data(), 2, id.Data(), 3);
    sim->Project(name + "_electron", expression, select);

    select = Form("%s(%s%d || %s%d)", selection.Data(), id.Data(), 5, id.Data(), 6);
    sim->Project(name + "_muon", expression, select);

    select = Form("%s(%s%d || %s%d)", selection.Data(), id.Data(), 14, id.Data(), 15);
    sim->Project(name + "_proton", expression, select);

    select = Form("%s(%s%d || %s%d)", selection.Data(), id.Data(), 13, id.Data(), 25);
    sim->Project(name + "_neutron", expression, select);

    TCanvas *density_plot = new TCanvas(name, title, 
                                        fDensityOpts.width, fDensityOpts.height);
    for (Int_t i = 0; i < entries; i++) {
        ((TH1D*)fDensity->At(i))->Draw("lep same"); 
    } 
*/
}


//////////////////////////////////////////////////////////////////////////////

void Results::logBinning(TAxis *vAxis) {
    Int_t nbins = vAxis->GetNbins();
    Axis_t amin = vAxis->GetXmin();
    Axis_t amax = vAxis->GetXmax();
    Axis_t width = (amax - amin) / (Double_t) nbins;

    Axis_t *logbins = new Axis_t[nbins + 1];
    for (int i = 0; i <= nbins; i++) 
        logbins[i] = TMath::Power(10., amin + i * width);
        
    vAxis->Set(nbins, logbins);
    delete logbins;    
}


//////////////////////////////////////////////////////////////////////////////

void Results::initializeHistograms() {
    TString  name;
    TString  title;
    TString  xtitle;
    TString  ytitle;
    Bool_t   stats;
    Int_t    xbins;
    Int_t    ybins;
    Double_t xmin;
    Double_t ymin;
    Double_t xmax;
    Double_t ymax;
    Int_t    entries;


    // Longitudinal Density
    //------------------------------------------------------------------------
    name   = fDensityOpts.name; 
    title  = fDensityOpts.title;
    xtitle = fDensityOpts.xtitle;
    ytitle = fDensityOpts.ytitle;
    title  = title + ";" + xtitle + ";" + ytitle;
    stats  = fDensityOpts.stats;
    xbins  = fDensityOpts.xbins;
    xmin   = fDensityOpts.xmin;
    xmax   = fDensityOpts.xmax;
    TH1D *den_photon   = new TH1D(name + "photon",   title, xbins, xmin, xmax);
    TH1D *den_electron = new TH1D(name + "electron", title, xbins, xmin, xmax);
    TH1D *den_muon     = new TH1D(name + "muon",     title, xbins, xmin, xmax);
    TH1D *den_proton   = new TH1D(name + "proton",   title, xbins, xmin, xmax);
    TH1D *den_neutron  = new TH1D(name + "neutron",  title, xbins, xmin, xmax);
    TH1D *den_ocharged = new TH1D(name + "ocharged", title, xbins, xmin, xmax);
    TH1D *den_oneutral = new TH1D(name + "oneutral", title, xbins, xmin, xmax);
    TH1D *den_nuclei   = new TH1D(name + "nuclei",   title, xbins, xmin, xmax);
    den_photon  ->SetLineColor(fColors.photon);
    den_electron->SetLineColor(fColors.electron);
    den_muon    ->SetLineColor(fColors.muon);
    den_proton  ->SetLineColor(fColors.proton);
    den_neutron ->SetLineColor(fColors.neutron);
    den_ocharged->SetLineColor(fColors.ocharged);
    den_oneutral->SetLineColor(fColors.oneutral);
    den_nuclei  ->SetLineColor(fColors.nuclei);
    fDensity = new TList();
    fDensity->Add(den_photon);
    fDensity->Add(den_electron);
    fDensity->Add(den_muon);
    fDensity->Add(den_proton);
    fDensity->Add(den_neutron);
    fDensity->Add(den_ocharged);
    fDensity->Add(den_oneutral);
    fDensity->Add(den_nuclei);
    entries = fDensity->GetEntries();
    for (Int_t i = 0; i < entries; i++)
        ((TH1D*)fDensity->At(i))->SetStats(stats);


    // X-Y Slice
    //------------------------------------------------------------------------
    name   = fSliceOpts.name;
    title  = fSliceOpts.title;
    xtitle = fSliceOpts.xtitle;
    ytitle = fSliceOpts.ytitle;
    title  = title + ";" + xtitle + ";" + ytitle;
    stats  = fSliceOpts.stats;
    xbins  = fSliceOpts.xbins;
    ybins  = fSliceOpts.ybins;
    xmin   = fSliceOpts.xmin;
    ymin   = fSliceOpts.ymin;
    xmax   = fSliceOpts.xmax;
    ymax   = fSliceOpts.ymax;
    TH2D *sli_photon   = new TH2D(name + "photon",   "photon"   + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_electron = new TH2D(name + "electron", "electron" + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_muon     = new TH2D(name + "muon",     "muon"     + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_proton   = new TH2D(name + "proton",   "proton"   + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_neutron  = new TH2D(name + "neutron",  "neutron"  + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_ocharged = new TH2D(name + "ocharged", "ocharged" + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_oneutral = new TH2D(name + "oneutral", "oneutral" + title, xbins, xmin, xmax, ybins, ymin, ymax);
    TH2D *sli_nuclei   = new TH2D(name + "nuclei",   "nuclei"   + title, xbins, xmin, xmax, ybins, ymin, ymax);
    fSlice = new TList();
    fSlice->Add(sli_photon);
    fSlice->Add(sli_electron);
    fSlice->Add(sli_muon);
    fSlice->Add(sli_proton);
    fSlice->Add(sli_neutron);
    fSlice->Add(sli_ocharged);
    fSlice->Add(sli_oneutral);
    fSlice->Add(sli_nuclei);
    entries = fDensity->GetEntries();
    for (Int_t i = 0; i < entries; i++)
        ((TH2D*)fDensity->At(i))->SetStats(stats);


    // Energy Spectrum
    //------------------------------------------------------------------------
    name   = fSpectrumOpts.name;
    title  = fSpectrumOpts.title;
    xtitle = fSpectrumOpts.xtitle;
    ytitle = fSpectrumOpts.ytitle;
    title  = title + ";" + xtitle + ";" + ytitle;
    stats  = fSpectrumOpts.stats;
    xbins  = fSpectrumOpts.xbins;
    xmin   = fSpectrumOpts.xmin;
    xmax   = fSpectrumOpts.xmax;
    TH1D *spe_photon   = new TH1D(name + "photon",   title, xbins, xmin, xmax);
    TH1D *spe_electron = new TH1D(name + "electron", title, xbins, xmin, xmax);
    TH1D *spe_muon     = new TH1D(name + "muon",     title, xbins, xmin, xmax);
    TH1D *spe_proton   = new TH1D(name + "proton",   title, xbins, xmin, xmax);
    TH1D *spe_neutron  = new TH1D(name + "neutron",  title, xbins, xmin, xmax);
    TH1D *spe_ocharged = new TH1D(name + "ocharged", title, xbins, xmin, xmax);
    TH1D *spe_oneutral = new TH1D(name + "oneutral", title, xbins, xmin, xmax);
    TH1D *spe_nuclei   = new TH1D(name + "nuclei",   title, xbins, xmin, xmax);
    spe_photon  ->SetLineColor(fColors.photon);
    spe_electron->SetLineColor(fColors.electron);
    spe_muon    ->SetLineColor(fColors.muon);
    spe_proton  ->SetLineColor(fColors.proton);
    spe_neutron ->SetLineColor(fColors.neutron);
    spe_ocharged->SetLineColor(fColors.ocharged);
    spe_oneutral->SetLineColor(fColors.oneutral);
    spe_nuclei  ->SetLineColor(fColors.nuclei);
    fSpectrum = new TList();
    fSpectrum->Add(spe_photon);
    fSpectrum->Add(spe_electron);
    fSpectrum->Add(spe_muon);
    fSpectrum->Add(spe_proton);
    fSpectrum->Add(spe_neutron);
    fSpectrum->Add(spe_ocharged);
    fSpectrum->Add(spe_oneutral);
    fSpectrum->Add(spe_nuclei);
    entries = fSpectrum->GetEntries();
    for (Int_t i = 0; i < entries; i++) {
        logBinning( ((TH1D*)fSpectrum->At(i))->GetXaxis() );
        ((TH1D*)fSpectrum->At(i))->SetStats(stats);
    }

    
    // Particle Content
    //------------------------------------------------------------------------
    name   = fContentOpts.name;
    title  = fContentOpts.title;
    ytitle = fContentOpts.ytitle;
    title  = title + ";;" + ytitle;
    stats  = fContentOpts.stats;
    xbins  = fContentOpts.xbins; 
    fContent = new TH1D(name, title, xbins, 1, xbins+1);
    for (Int_t i = 1; i <= xbins; i++)
        fContent->GetXaxis()->SetBinLabel(i, fContentOpts.GetLabel(i));


    // First Impact
    //------------------------------------------------------------------------
    name   = fImpactOpts.name;
    title  = fImpactOpts.title;
    ytitle = fImpactOpts.ytitle;
    title  = title + ";;" + ytitle;
    stats  = fImpactOpts.stats;
    xbins  = fImpactOpts.xbins; 
    ybins  = fImpactOpts.ybins;
    xmin   = fImpactOpts.xmin;
    ymin   = fImpactOpts.ymin;
    xmax   = fImpactOpts.xmax;
    ymax   = fImpactOpts.ymax;
    fImpact = new TH2D(name, title, xbins, xmin, xmax, ybins, ymin, ymax);
    logBinning(fImpact->GetYaxis());
    for (Int_t i = 1; i <= xbins; i++)
        fImpact->GetXaxis()->SetBinLabel(i, fImpactOpts.GetLabel(i));

}


//////////////////////////////////////////////////////////////////////////////

void Results::printStdErr() {
    cout << " [ERROR] Sorry, data file was not successfully loaded.\n";
    cout << "         You'll need to try again from scratch.\n";
}

