//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug 24 14:02:13 2019 by ROOT version 5.34/38
// from TTree run/simulation run
// found on file: DAT000001.root
// Then fucked with by Eric
//////////////////////////////////////////////////////////

#ifndef run_h
#define run_h

#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxrun = 1;

class Run {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //crsIO::TRun     *run_;
   UInt_t          run_TObject_fUniqueID;
   UInt_t          run_TObject_fBits;
   Int_t           run_RunID;
   Float_t         run_Date;
   Float_t         run_Version;
   vector<float>   run_ObservationLevel;
   Int_t           run_ParticleID;
   Float_t         run_EnergySlope;
   Float_t         run_EnergyMin;
   Float_t         run_EnergyMax;
   Float_t         run_AzimuthMin;
   Float_t         run_AzimuthMax;
   Float_t         run_ZenithMin;
   Float_t         run_ZenithMax;
   Float_t         run_CutHadrons;
   Float_t         run_CutMuons;
   Float_t         run_CutElectrons;
   Float_t         run_CutPhotons;
   vector<float>   run_AtmosphereA;
   vector<float>   run_AtmosphereB;
   vector<float>   run_AtmosphereC;
   Float_t         run_BFieldX;
   Float_t         run_BFieldZ;
   Bool_t          run_EGS4;
   Bool_t          run_NKG;
   Bool_t          run_Cherenkov;
   Bool_t          run_Neutrino;
   Bool_t          run_Curved;
   Bool_t          run_MuonAdditionalInfo;
   Bool_t          run_MuonMultScatteringMoliere;
   Float_t         run_RadialRangeNKG;
   Float_t         run_TransitionEnergy;
   Int_t           run_LowEHadModel;
   Int_t           run_HighEHadModel;
   Int_t           run_VersionSIBYLL_interaction;
   Int_t           run_VersionSIBYLL_crosssection;
   Int_t           run_VersionQGSJET_interaction;
   Int_t           run_VersionQGSJET_crosssection;
   Int_t           run_VersionDPMJET_interaction;
   Int_t           run_VersionDPMJET_crosssection;
   Int_t           run_VersionVENUSNEXUS_crosssection;
   Float_t         run_NFLAIN;
   Float_t         run_NFLDIF;
   Float_t         run_NFLPI0;
   Float_t         run_NFLPIF;
   Float_t         run_NFLCHE;
   Float_t         run_NFRAGM;
   Float_t         run_HadronicThinningFraction;
   Float_t         run_EMThinningFraction;
   Int_t           run_HadronicThinningtLimit;
   Int_t           run_EMThinningLimit;
   Float_t         run_RadialThiningRMax;
   Float_t         run_ViewConeMin;
   Float_t         run_ViewConeMax;
   Int_t           run_nCherenkovDetectorsX;
   Int_t           run_nCherenkovDetectorsY;
   Float_t         run_GridCherenkovDetectorX;
   Float_t         run_GridCherenkovDetectorY;
   Float_t         run_LengthCherenkovDetectorX;
   Float_t         run_LengthCherenkovDetectorY;
   Bool_t          run_CherenkovOutputSeparate;
   Float_t         run_CherenkovBandwidthMin;
   Float_t         run_CherenkovBandwidthMax;
   Int_t           run_nUseCherenkovEvent;
   vector<float>   run_CherenkovCoreX;
   vector<float>   run_CherenkovCoreY;
   Float_t         run_OrientationArray;
   Float_t         run_StepLengthFactorMultiScatter;
   Int_t           run_Computer;

   // List of branches
   TBranch        *b_run_TObject_fUniqueID;   //!
   TBranch        *b_run_TObject_fBits;   //!
   TBranch        *b_run_RunID;   //!
   TBranch        *b_run_Date;   //!
   TBranch        *b_run_Version;   //!
   TBranch        *b_run_ObservationLevel;   //!
   TBranch        *b_run_ParticleID;   //!
   TBranch        *b_run_EnergySlope;   //!
   TBranch        *b_run_EnergyMin;   //!
   TBranch        *b_run_EnergyMax;   //!
   TBranch        *b_run_AzimuthMin;   //!
   TBranch        *b_run_AzimuthMax;   //!
   TBranch        *b_run_ZenithMin;   //!
   TBranch        *b_run_ZenithMax;   //!
   TBranch        *b_run_CutHadrons;   //!
   TBranch        *b_run_CutMuons;   //!
   TBranch        *b_run_CutElectrons;   //!
   TBranch        *b_run_CutPhotons;   //!
   TBranch        *b_run_AtmosphereA;   //!
   TBranch        *b_run_AtmosphereB;   //!
   TBranch        *b_run_AtmosphereC;   //!
   TBranch        *b_run_BFieldX;   //!
   TBranch        *b_run_BFieldZ;   //!
   TBranch        *b_run_EGS4;   //!
   TBranch        *b_run_NKG;   //!
   TBranch        *b_run_Cherenkov;   //!
   TBranch        *b_run_Neutrino;   //!
   TBranch        *b_run_Curved;   //!
   TBranch        *b_run_MuonAdditionalInfo;   //!
   TBranch        *b_run_MuonMultScatteringMoliere;   //!
   TBranch        *b_run_RadialRangeNKG;   //!
   TBranch        *b_run_TransitionEnergy;   //!
   TBranch        *b_run_LowEHadModel;   //!
   TBranch        *b_run_HighEHadModel;   //!
   TBranch        *b_run_VersionSIBYLL_interaction;   //!
   TBranch        *b_run_VersionSIBYLL_crosssection;   //!
   TBranch        *b_run_VersionQGSJET_interaction;   //!
   TBranch        *b_run_VersionQGSJET_crosssection;   //!
   TBranch        *b_run_VersionDPMJET_interaction;   //!
   TBranch        *b_run_VersionDPMJET_crosssection;   //!
   TBranch        *b_run_VersionVENUSNEXUS_crosssection;   //!
   TBranch        *b_run_NFLAIN;   //!
   TBranch        *b_run_NFLDIF;   //!
   TBranch        *b_run_NFLPI0;   //!
   TBranch        *b_run_NFLPIF;   //!
   TBranch        *b_run_NFLCHE;   //!
   TBranch        *b_run_NFRAGM;   //!
   TBranch        *b_run_HadronicThinningFraction;   //!
   TBranch        *b_run_EMThinningFraction;   //!
   TBranch        *b_run_HadronicThinningtLimit;   //!
   TBranch        *b_run_EMThinningLimit;   //!
   TBranch        *b_run_RadialThiningRMax;   //!
   TBranch        *b_run_ViewConeMin;   //!
   TBranch        *b_run_ViewConeMax;   //!
   TBranch        *b_run_nCherenkovDetectorsX;   //!
   TBranch        *b_run_nCherenkovDetectorsY;   //!
   TBranch        *b_run_GridCherenkovDetectorX;   //!
   TBranch        *b_run_GridCherenkovDetectorY;   //!
   TBranch        *b_run_LengthCherenkovDetectorX;   //!
   TBranch        *b_run_LengthCherenkovDetectorY;   //!
   TBranch        *b_run_CherenkovOutputSeparate;   //!
   TBranch        *b_run_CherenkovBandwidthMin;   //!
   TBranch        *b_run_CherenkovBandwidthMax;   //!
   TBranch        *b_run_nUseCherenkovEvent;   //!
   TBranch        *b_run_CherenkovCoreX;   //!
   TBranch        *b_run_CherenkovCoreY;   //!
   TBranch        *b_run_OrientationArray;   //!
   TBranch        *b_run_StepLengthFactorMultiScatter;   //!
   TBranch        *b_run_Computer;   //!

   Run(TTree *tree=0);
   virtual ~Run();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
};

#endif

Run::Run(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), die.
   if (tree == 0) {
       std::cout << " [ERROR] No TTree or TChain\n";
      return;
   }
   Init(tree);
}

Run::~Run()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Run::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void Run::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run.TObject.fUniqueID", &run_TObject_fUniqueID, &b_run_TObject_fUniqueID);
   fChain->SetBranchAddress("run.TObject.fBits", &run_TObject_fBits, &b_run_TObject_fBits);
   fChain->SetBranchAddress("run.RunID", &run_RunID, &b_run_RunID);
   fChain->SetBranchAddress("run.Date", &run_Date, &b_run_Date);
   fChain->SetBranchAddress("run.Version", &run_Version, &b_run_Version);
   fChain->SetBranchAddress("run.ObservationLevel", &run_ObservationLevel, &b_run_ObservationLevel);
   fChain->SetBranchAddress("run.ParticleID", &run_ParticleID, &b_run_ParticleID);
   fChain->SetBranchAddress("run.EnergySlope", &run_EnergySlope, &b_run_EnergySlope);
   fChain->SetBranchAddress("run.EnergyMin", &run_EnergyMin, &b_run_EnergyMin);
   fChain->SetBranchAddress("run.EnergyMax", &run_EnergyMax, &b_run_EnergyMax);
   fChain->SetBranchAddress("run.AzimuthMin", &run_AzimuthMin, &b_run_AzimuthMin);
   fChain->SetBranchAddress("run.AzimuthMax", &run_AzimuthMax, &b_run_AzimuthMax);
   fChain->SetBranchAddress("run.ZenithMin", &run_ZenithMin, &b_run_ZenithMin);
   fChain->SetBranchAddress("run.ZenithMax", &run_ZenithMax, &b_run_ZenithMax);
   fChain->SetBranchAddress("run.CutHadrons", &run_CutHadrons, &b_run_CutHadrons);
   fChain->SetBranchAddress("run.CutMuons", &run_CutMuons, &b_run_CutMuons);
   fChain->SetBranchAddress("run.CutElectrons", &run_CutElectrons, &b_run_CutElectrons);
   fChain->SetBranchAddress("run.CutPhotons", &run_CutPhotons, &b_run_CutPhotons);
   fChain->SetBranchAddress("run.AtmosphereA", &run_AtmosphereA, &b_run_AtmosphereA);
   fChain->SetBranchAddress("run.AtmosphereB", &run_AtmosphereB, &b_run_AtmosphereB);
   fChain->SetBranchAddress("run.AtmosphereC", &run_AtmosphereC, &b_run_AtmosphereC);
   fChain->SetBranchAddress("run.BFieldX", &run_BFieldX, &b_run_BFieldX);
   fChain->SetBranchAddress("run.BFieldZ", &run_BFieldZ, &b_run_BFieldZ);
   fChain->SetBranchAddress("run.EGS4", &run_EGS4, &b_run_EGS4);
   fChain->SetBranchAddress("run.NKG", &run_NKG, &b_run_NKG);
   fChain->SetBranchAddress("run.Cherenkov", &run_Cherenkov, &b_run_Cherenkov);
   fChain->SetBranchAddress("run.Neutrino", &run_Neutrino, &b_run_Neutrino);
   fChain->SetBranchAddress("run.Curved", &run_Curved, &b_run_Curved);
   fChain->SetBranchAddress("run.MuonAdditionalInfo", &run_MuonAdditionalInfo, &b_run_MuonAdditionalInfo);
   fChain->SetBranchAddress("run.MuonMultScatteringMoliere", &run_MuonMultScatteringMoliere, &b_run_MuonMultScatteringMoliere);
   fChain->SetBranchAddress("run.RadialRangeNKG", &run_RadialRangeNKG, &b_run_RadialRangeNKG);
   fChain->SetBranchAddress("run.TransitionEnergy", &run_TransitionEnergy, &b_run_TransitionEnergy);
   fChain->SetBranchAddress("run.LowEHadModel", &run_LowEHadModel, &b_run_LowEHadModel);
   fChain->SetBranchAddress("run.HighEHadModel", &run_HighEHadModel, &b_run_HighEHadModel);
   fChain->SetBranchAddress("run.VersionSIBYLL_interaction", &run_VersionSIBYLL_interaction, &b_run_VersionSIBYLL_interaction);
   fChain->SetBranchAddress("run.VersionSIBYLL_crosssection", &run_VersionSIBYLL_crosssection, &b_run_VersionSIBYLL_crosssection);
   fChain->SetBranchAddress("run.VersionQGSJET_interaction", &run_VersionQGSJET_interaction, &b_run_VersionQGSJET_interaction);
   fChain->SetBranchAddress("run.VersionQGSJET_crosssection", &run_VersionQGSJET_crosssection, &b_run_VersionQGSJET_crosssection);
   fChain->SetBranchAddress("run.VersionDPMJET_interaction", &run_VersionDPMJET_interaction, &b_run_VersionDPMJET_interaction);
   fChain->SetBranchAddress("run.VersionDPMJET_crosssection", &run_VersionDPMJET_crosssection, &b_run_VersionDPMJET_crosssection);
   fChain->SetBranchAddress("run.VersionVENUSNEXUS_crosssection", &run_VersionVENUSNEXUS_crosssection, &b_run_VersionVENUSNEXUS_crosssection);
   fChain->SetBranchAddress("run.NFLAIN", &run_NFLAIN, &b_run_NFLAIN);
   fChain->SetBranchAddress("run.NFLDIF", &run_NFLDIF, &b_run_NFLDIF);
   fChain->SetBranchAddress("run.NFLPI0", &run_NFLPI0, &b_run_NFLPI0);
   fChain->SetBranchAddress("run.NFLPIF", &run_NFLPIF, &b_run_NFLPIF);
   fChain->SetBranchAddress("run.NFLCHE", &run_NFLCHE, &b_run_NFLCHE);
   fChain->SetBranchAddress("run.NFRAGM", &run_NFRAGM, &b_run_NFRAGM);
   fChain->SetBranchAddress("run.HadronicThinningFraction", &run_HadronicThinningFraction, &b_run_HadronicThinningFraction);
   fChain->SetBranchAddress("run.EMThinningFraction", &run_EMThinningFraction, &b_run_EMThinningFraction);
   fChain->SetBranchAddress("run.HadronicThinningtLimit", &run_HadronicThinningtLimit, &b_run_HadronicThinningtLimit);
   fChain->SetBranchAddress("run.EMThinningLimit", &run_EMThinningLimit, &b_run_EMThinningLimit);
   fChain->SetBranchAddress("run.RadialThiningRMax", &run_RadialThiningRMax, &b_run_RadialThiningRMax);
   fChain->SetBranchAddress("run.ViewConeMin", &run_ViewConeMin, &b_run_ViewConeMin);
   fChain->SetBranchAddress("run.ViewConeMax", &run_ViewConeMax, &b_run_ViewConeMax);
   fChain->SetBranchAddress("run.nCherenkovDetectorsX", &run_nCherenkovDetectorsX, &b_run_nCherenkovDetectorsX);
   fChain->SetBranchAddress("run.nCherenkovDetectorsY", &run_nCherenkovDetectorsY, &b_run_nCherenkovDetectorsY);
   fChain->SetBranchAddress("run.GridCherenkovDetectorX", &run_GridCherenkovDetectorX, &b_run_GridCherenkovDetectorX);
   fChain->SetBranchAddress("run.GridCherenkovDetectorY", &run_GridCherenkovDetectorY, &b_run_GridCherenkovDetectorY);
   fChain->SetBranchAddress("run.LengthCherenkovDetectorX", &run_LengthCherenkovDetectorX, &b_run_LengthCherenkovDetectorX);
   fChain->SetBranchAddress("run.LengthCherenkovDetectorY", &run_LengthCherenkovDetectorY, &b_run_LengthCherenkovDetectorY);
   fChain->SetBranchAddress("run.CherenkovOutputSeparate", &run_CherenkovOutputSeparate, &b_run_CherenkovOutputSeparate);
   fChain->SetBranchAddress("run.CherenkovBandwidthMin", &run_CherenkovBandwidthMin, &b_run_CherenkovBandwidthMin);
   fChain->SetBranchAddress("run.CherenkovBandwidthMax", &run_CherenkovBandwidthMax, &b_run_CherenkovBandwidthMax);
   fChain->SetBranchAddress("run.nUseCherenkovEvent", &run_nUseCherenkovEvent, &b_run_nUseCherenkovEvent);
   fChain->SetBranchAddress("run.CherenkovCoreX", &run_CherenkovCoreX, &b_run_CherenkovCoreX);
   fChain->SetBranchAddress("run.CherenkovCoreY", &run_CherenkovCoreY, &b_run_CherenkovCoreY);
   fChain->SetBranchAddress("run.OrientationArray", &run_OrientationArray, &b_run_OrientationArray);
   fChain->SetBranchAddress("run.StepLengthFactorMultiScatter", &run_StepLengthFactorMultiScatter, &b_run_StepLengthFactorMultiScatter);
   fChain->SetBranchAddress("run.Computer", &run_Computer, &b_run_Computer);
}

