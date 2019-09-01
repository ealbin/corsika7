//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug 24 14:02:19 2019 by ROOT version 5.34/38
// from TTree sim/simulated showers
// found on file: DAT000001.root
// Then fucked with by Eric
//////////////////////////////////////////////////////////

#ifndef sim_h
#define sim_h

#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObject.h>
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxshower = 1;
   const Int_t kMaxlong_ = 1;
   const Int_t kMaxcherenkov_ = 1;

  Int_t kMaxparticle_ = 0; 

class Sim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
 //crsIO::TShower  *shower_;
   UInt_t          shower_TObject_fUniqueID;
   UInt_t          shower_TObject_fBits;
   Int_t           shower_EventID;
   Float_t         shower_Energy;
   Float_t         shower_StartingAltitude;
   Int_t           shower_FirstTarget;
   Float_t         shower_FirstHeight;
   Float_t         shower_Theta;
   Float_t         shower_Phi;
   Int_t           shower_RandomSeed[10];
   Int_t           shower_RandomOffset[10];
   Float_t         shower_nPhotons;
   Float_t         shower_nElectrons;
   Float_t         shower_nHadrons;
   Float_t         shower_nMuons;
   Int_t           shower_nParticlesWritten;
   Int_t           shower_nPhotonsWritten;
   Int_t           shower_nElectronsWritten;
   Int_t           shower_nHadronsWritten;
   Int_t           shower_nMuonsWritten;
   Float_t         shower_GH_Nmax;
   Float_t         shower_GH_t0;
   Float_t         shower_GH_tmax;
   Float_t         shower_GH_a;
   Float_t         shower_GH_b;
   Float_t         shower_GH_c;
   Float_t         shower_GH_Chi2;
   Int_t           shower_nPreshower;
   Int_t           shower_CPUtime;
   Int_t           particle__;
   UInt_t          *particle__fUniqueID;//[kMaxparticle_];   //[particle._]
   UInt_t          *particle__fBits;//[kMaxparticle_];   //[particle._]
   Int_t           *particle__CorsikaID;//[kMaxparticle_];   //[particle._]
   Int_t           *particle__ParticleID;//[kMaxparticle_];   //[particle._]
   Int_t           *particle__ObservationLevel;//[kMaxparticle_];   //[particle._]
   Int_t           *particle__HadronicGeneration;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__Px;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__Py;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__Pz;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__x;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__y;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__Time;//[kMaxparticle_];   //[particle._]
   Double_t        *particle__Weight;//[kMaxparticle_];   //[particle._]
   Int_t           long__;
   UInt_t          long__fUniqueID[kMaxlong_];   //[long._]
   UInt_t          long__fBits[kMaxlong_];   //[long._]
   Float_t         long__Depth[kMaxlong_];   //[long._]
   ULong64_t       long__nGammas[kMaxlong_];   //[long._]
   ULong64_t       long__nElectrons[kMaxlong_];   //[long._]
   ULong64_t       long__nPositrons[kMaxlong_];   //[long._]
   ULong64_t       long__nMuons[kMaxlong_];   //[long._]
   ULong64_t       long__nAntiMuons[kMaxlong_];   //[long._]
   ULong64_t       long__nHadrons[kMaxlong_];   //[long._]
   ULong64_t       long__nCharged[kMaxlong_];   //[long._]
   ULong64_t       long__nNuclei[kMaxlong_];   //[long._]
   ULong64_t       long__nCherenkov[kMaxlong_];   //[long._]
   Int_t           cherenkov__;
   UInt_t          cherenkov__fUniqueID[kMaxcherenkov_];   //[cherenkov._]
   UInt_t          cherenkov__fBits[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__nPhotons[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__x[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__y[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__u[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__v[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__Time[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__ProductionHeight[kMaxcherenkov_];   //[cherenkov._]
   Float_t         cherenkov__Weight[kMaxcherenkov_];   //[cherenkov._]

   // List of branches
   TBranch        *b_shower_TObject_fUniqueID;   //!
   TBranch        *b_shower_TObject_fBits;   //!
   TBranch        *b_shower_EventID;   //!
   TBranch        *b_shower_Energy;   //!
   TBranch        *b_shower_StartingAltitude;   //!
   TBranch        *b_shower_FirstTarget;   //!
   TBranch        *b_shower_FirstHeight;   //!
   TBranch        *b_shower_Theta;   //!
   TBranch        *b_shower_Phi;   //!
   TBranch        *b_shower_RandomSeed;   //!
   TBranch        *b_shower_RandomOffset;   //!
   TBranch        *b_shower_nPhotons;   //!
   TBranch        *b_shower_nElectrons;   //!
   TBranch        *b_shower_nHadrons;   //!
   TBranch        *b_shower_nMuons;   //!
   TBranch        *b_shower_nParticlesWritten;   //!
   TBranch        *b_shower_nPhotonsWritten;   //!
   TBranch        *b_shower_nElectronsWritten;   //!
   TBranch        *b_shower_nHadronsWritten;   //!
   TBranch        *b_shower_nMuonsWritten;   //!
   TBranch        *b_shower_GH_Nmax;   //!
   TBranch        *b_shower_GH_t0;   //!
   TBranch        *b_shower_GH_tmax;   //!
   TBranch        *b_shower_GH_a;   //!
   TBranch        *b_shower_GH_b;   //!
   TBranch        *b_shower_GH_c;   //!
   TBranch        *b_shower_GH_Chi2;   //!
   TBranch        *b_shower_nPreshower;   //!
   TBranch        *b_shower_CPUtime;   //!
   TBranch        *b_particle__;   //!
   TBranch        *b_particle__fUniqueID;   //!
   TBranch        *b_particle__fBits;   //!
   TBranch        *b_particle__CorsikaID;   //!
   TBranch        *b_particle__ParticleID;   //!
   TBranch        *b_particle__ObservationLevel;   //!
   TBranch        *b_particle__HadronicGeneration;   //!
   TBranch        *b_particle__Px;   //!
   TBranch        *b_particle__Py;   //!
   TBranch        *b_particle__Pz;   //!
   TBranch        *b_particle__x;   //!
   TBranch        *b_particle__y;   //!
   TBranch        *b_particle__Time;   //!
   TBranch        *b_particle__Weight;   //!
   TBranch        *b_long__;   //!
   TBranch        *b_long__fUniqueID;   //!
   TBranch        *b_long__fBits;   //!
   TBranch        *b_long__Depth;   //!
   TBranch        *b_long__nGammas;   //!
   TBranch        *b_long__nElectrons;   //!
   TBranch        *b_long__nPositrons;   //!
   TBranch        *b_long__nMuons;   //!
   TBranch        *b_long__nAntiMuons;   //!
   TBranch        *b_long__nHadrons;   //!
   TBranch        *b_long__nCharged;   //!
   TBranch        *b_long__nNuclei;   //!
   TBranch        *b_long__nCherenkov;   //!
   TBranch        *b_cherenkov__;   //!
   TBranch        *b_cherenkov__fUniqueID;   //!
   TBranch        *b_cherenkov__fBits;   //!
   TBranch        *b_cherenkov__nPhotons;   //!
   TBranch        *b_cherenkov__x;   //!
   TBranch        *b_cherenkov__y;   //!
   TBranch        *b_cherenkov__u;   //!
   TBranch        *b_cherenkov__v;   //!
   TBranch        *b_cherenkov__Time;   //!
   TBranch        *b_cherenkov__ProductionHeight;   //!
   TBranch        *b_cherenkov__Weight;   //!

   Sim(TTree *tree=0);
   virtual Long64_t GetEntries();
   virtual Int_t    GetEntry(Long64_t entry);

private:
   virtual void     Init();
   virtual void     Allocate();
   virtual void     Link();

   virtual         ~Sim();
   virtual void     Free();
};

#endif

Sim::Sim(TTree *tree) : fChain(tree) 
{
// if parameter tree is not specified (or zero), die. 
   if (fChain == 0) {
        std::cout << " [ERROR] No TTree or TChain\n";
        return;
   }
   if (fChain->GetEntries() == 0) {
        std::cout << " [ERROR] TTree or TChain has 0 entries\n";
        return;
   }
   
   Init();
}

Long64_t Sim::GetEntries()
{
   return fChain->GetEntries();
}

Int_t Sim::GetEntry(Long64_t entry)
{
    return fChain->GetEntry(entry);
}


void Sim::Init()
{
   fChain->SetMakeClass(1);
   fChain->SetBranchStatus("*", 0);
   fChain->SetBranchStatus("particle.", 1);
   fChain->SetBranchAddress("particle.", &particle__, &b_particle__);
   Long64_t entries = fChain->GetEntries();
   
   for (Long64_t i = 0; i < entries; i++) {
        fChain->GetEntry(i);
        if (particle__ > kMaxparticle_)
            kMaxparticle_ = particle__;
   }
   fChain->SetBranchStatus("*", 1);
   Allocate();
   Link();
}

void Sim::Allocate()
{
   std::cout << "Allocating for " << kMaxparticle_ << " (maximal) particles\n";
   particle__fUniqueID          = new UInt_t[kMaxparticle_];
   particle__fBits              = new UInt_t[kMaxparticle_];
   particle__CorsikaID          = new Int_t[kMaxparticle_];
   particle__ParticleID         = new Int_t[kMaxparticle_];
   particle__ObservationLevel   = new Int_t[kMaxparticle_];
   particle__HadronicGeneration = new Int_t[kMaxparticle_];
   particle__Px                 = new Double_t[kMaxparticle_];
   particle__Py                 = new Double_t[kMaxparticle_];
   particle__Pz                 = new Double_t[kMaxparticle_];
   particle__x                  = new Double_t[kMaxparticle_];
   particle__y                  = new Double_t[kMaxparticle_];
   particle__Time               = new Double_t[kMaxparticle_];
   particle__Weight             = new Double_t[kMaxparticle_];

   for (Int_t i = 0; i < kMaxparticle_; i++) {
       particle__fUniqueID[i]          = 0;
       particle__fBits[i]              = 0;
       particle__CorsikaID[i]          = 0;
       particle__ParticleID[i]         = 0;
       particle__ObservationLevel[i]   = 0;
       particle__HadronicGeneration[i] = 0;
       particle__Px[i]                 = 0.;
       particle__Py[i]                 = 0.;
       particle__Pz[i]                 = 0.;
       particle__x[i]                  = 0.;
       particle__y[i]                  = 0.;
       particle__Time[i]               = 0.;
       particle__Weight[i]             = 0.;
   }
}

void Sim::Link()
{
   // Set branch addresses and branch pointers
   fChain->SetBranchAddress("shower.TObject.fUniqueID", &shower_TObject_fUniqueID, &b_shower_TObject_fUniqueID);
   fChain->SetBranchAddress("shower.TObject.fBits", &shower_TObject_fBits, &b_shower_TObject_fBits);
   fChain->SetBranchAddress("shower.EventID", &shower_EventID, &b_shower_EventID);
   fChain->SetBranchAddress("shower.Energy", &shower_Energy, &b_shower_Energy);
   fChain->SetBranchAddress("shower.StartingAltitude", &shower_StartingAltitude, &b_shower_StartingAltitude);
   fChain->SetBranchAddress("shower.FirstTarget", &shower_FirstTarget, &b_shower_FirstTarget);
   fChain->SetBranchAddress("shower.FirstHeight", &shower_FirstHeight, &b_shower_FirstHeight);
   fChain->SetBranchAddress("shower.Theta", &shower_Theta, &b_shower_Theta);
   fChain->SetBranchAddress("shower.Phi", &shower_Phi, &b_shower_Phi);
   fChain->SetBranchAddress("shower.RandomSeed[10]", shower_RandomSeed, &b_shower_RandomSeed);
   fChain->SetBranchAddress("shower.RandomOffset[10]", shower_RandomOffset, &b_shower_RandomOffset);
   fChain->SetBranchAddress("shower.nPhotons", &shower_nPhotons, &b_shower_nPhotons);
   fChain->SetBranchAddress("shower.nElectrons", &shower_nElectrons, &b_shower_nElectrons);
   fChain->SetBranchAddress("shower.nHadrons", &shower_nHadrons, &b_shower_nHadrons);
   fChain->SetBranchAddress("shower.nMuons", &shower_nMuons, &b_shower_nMuons);
   fChain->SetBranchAddress("shower.nParticlesWritten", &shower_nParticlesWritten, &b_shower_nParticlesWritten);
   fChain->SetBranchAddress("shower.nPhotonsWritten", &shower_nPhotonsWritten, &b_shower_nPhotonsWritten);
   fChain->SetBranchAddress("shower.nElectronsWritten", &shower_nElectronsWritten, &b_shower_nElectronsWritten);
   fChain->SetBranchAddress("shower.nHadronsWritten", &shower_nHadronsWritten, &b_shower_nHadronsWritten);
   fChain->SetBranchAddress("shower.nMuonsWritten", &shower_nMuonsWritten, &b_shower_nMuonsWritten);
   fChain->SetBranchAddress("shower.GH_Nmax", &shower_GH_Nmax, &b_shower_GH_Nmax);
   fChain->SetBranchAddress("shower.GH_t0", &shower_GH_t0, &b_shower_GH_t0);
   fChain->SetBranchAddress("shower.GH_tmax", &shower_GH_tmax, &b_shower_GH_tmax);
   fChain->SetBranchAddress("shower.GH_a", &shower_GH_a, &b_shower_GH_a);
   fChain->SetBranchAddress("shower.GH_b", &shower_GH_b, &b_shower_GH_b);
   fChain->SetBranchAddress("shower.GH_c", &shower_GH_c, &b_shower_GH_c);
   fChain->SetBranchAddress("shower.GH_Chi2", &shower_GH_Chi2, &b_shower_GH_Chi2);
   fChain->SetBranchAddress("shower.nPreshower", &shower_nPreshower, &b_shower_nPreshower);
   fChain->SetBranchAddress("shower.CPUtime", &shower_CPUtime, &b_shower_CPUtime);
   fChain->SetBranchAddress("particle.", &particle__, &b_particle__);
   fChain->SetBranchAddress("particle..fUniqueID", particle__fUniqueID, &b_particle__fUniqueID);
   fChain->SetBranchAddress("particle..fBits", particle__fBits, &b_particle__fBits);
   fChain->SetBranchAddress("particle..CorsikaID", particle__CorsikaID, &b_particle__CorsikaID);
   fChain->SetBranchAddress("particle..ParticleID", particle__ParticleID, &b_particle__ParticleID);
   fChain->SetBranchAddress("particle..ObservationLevel", particle__ObservationLevel, &b_particle__ObservationLevel);
   fChain->SetBranchAddress("particle..HadronicGeneration", particle__HadronicGeneration, &b_particle__HadronicGeneration);
   fChain->SetBranchAddress("particle..Px", particle__Px, &b_particle__Px);
   fChain->SetBranchAddress("particle..Py", particle__Py, &b_particle__Py);
   fChain->SetBranchAddress("particle..Pz", particle__Pz, &b_particle__Pz);
   fChain->SetBranchAddress("particle..x", particle__x, &b_particle__x);
   fChain->SetBranchAddress("particle..y", particle__y, &b_particle__y);
   fChain->SetBranchAddress("particle..Time", particle__Time, &b_particle__Time);
   fChain->SetBranchAddress("particle..Weight", particle__Weight, &b_particle__Weight);
   fChain->SetBranchAddress("long.", &long__, &b_long__);
   fChain->SetBranchAddress("long..fUniqueID", &long__fUniqueID, &b_long__fUniqueID);
   fChain->SetBranchAddress("long..fBits", &long__fBits, &b_long__fBits);
   fChain->SetBranchAddress("long..Depth", &long__Depth, &b_long__Depth);
   fChain->SetBranchAddress("long..nGammas", &long__nGammas, &b_long__nGammas);
   fChain->SetBranchAddress("long..nElectrons", &long__nElectrons, &b_long__nElectrons);
   fChain->SetBranchAddress("long..nPositrons", &long__nPositrons, &b_long__nPositrons);
   fChain->SetBranchAddress("long..nMuons", &long__nMuons, &b_long__nMuons);
   fChain->SetBranchAddress("long..nAntiMuons", &long__nAntiMuons, &b_long__nAntiMuons);
   fChain->SetBranchAddress("long..nHadrons", &long__nHadrons, &b_long__nHadrons);
   fChain->SetBranchAddress("long..nCharged", &long__nCharged, &b_long__nCharged);
   fChain->SetBranchAddress("long..nNuclei", &long__nNuclei, &b_long__nNuclei);
   fChain->SetBranchAddress("long..nCherenkov", &long__nCherenkov, &b_long__nCherenkov);
   fChain->SetBranchAddress("cherenkov.", &cherenkov__, &b_cherenkov__);
   fChain->SetBranchAddress("cherenkov..fUniqueID", &cherenkov__fUniqueID, &b_cherenkov__fUniqueID);
   fChain->SetBranchAddress("cherenkov..fBits", &cherenkov__fBits, &b_cherenkov__fBits);
   fChain->SetBranchAddress("cherenkov..nPhotons", &cherenkov__nPhotons, &b_cherenkov__nPhotons);
   fChain->SetBranchAddress("cherenkov..x", &cherenkov__x, &b_cherenkov__x);
   fChain->SetBranchAddress("cherenkov..y", &cherenkov__y, &b_cherenkov__y);
   fChain->SetBranchAddress("cherenkov..u", &cherenkov__u, &b_cherenkov__u);
   fChain->SetBranchAddress("cherenkov..v", &cherenkov__v, &b_cherenkov__v);
   fChain->SetBranchAddress("cherenkov..Time", &cherenkov__Time, &b_cherenkov__Time);
   fChain->SetBranchAddress("cherenkov..ProductionHeight", &cherenkov__ProductionHeight, &b_cherenkov__ProductionHeight);
   fChain->SetBranchAddress("cherenkov..Weight", &cherenkov__Weight, &b_cherenkov__Weight);
}

Sim::~Sim()
{
   if (!fChain) return;
   Free();
   delete fChain->GetCurrentFile();
}

void Sim::Free()
{
   if (particle__fUniqueID)          delete[] particle__fUniqueID;
   if (particle__fBits)              delete[] particle__fBits;
   if (particle__CorsikaID)          delete[] particle__CorsikaID;
   if (particle__ParticleID)         delete[] particle__ParticleID;
   if (particle__ObservationLevel)   delete[] particle__ObservationLevel;
   if (particle__HadronicGeneration) delete[] particle__HadronicGeneration;
   if (particle__Px)                 delete[] particle__Px;
   if (particle__Py)                 delete[] particle__Py;
   if (particle__Pz)                 delete[] particle__Pz;
   if (particle__x)                  delete[] particle__x;
   if (particle__y)                  delete[] particle__y;
   if (particle__Time)               delete[] particle__Time;
   if (particle__Weight)             delete[] particle__Weight;
   
   particle__fUniqueID          = 0;
   particle__fBits              = 0;
   particle__CorsikaID          = 0;
   particle__ParticleID         = 0;
   particle__ObservationLevel   = 0;
   particle__HadronicGeneration = 0;
   particle__Px                 = 0;
   particle__Py                 = 0;
   particle__Pz                 = 0;
   particle__x                  = 0;
   particle__y                  = 0;
   particle__Time               = 0;
   particle__Weight             = 0;
}

