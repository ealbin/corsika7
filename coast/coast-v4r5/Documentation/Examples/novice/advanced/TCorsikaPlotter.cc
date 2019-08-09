#include <TCorsikaPlotter.h>

#include <crs/TBlock.h>
#include <crs/TSubBlock.h>
#include <crs/MParticleBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MParticle.h>
#include <crs/MCherenkov.h>

#include <TFile.h>
#include <TH1D.h>
#include <TDirectory.h>

#include <string>
#include <iostream>
#include <sstream>


TCorsikaPlotter::TCorsikaPlotter ()
: fFileName ("") {
}


TCorsikaPlotter::TCorsikaPlotter (const std::string& filename)
: fFileName (filename) {
}


TCorsikaPlotter::~TCorsikaPlotter () {
}


void TCorsikaPlotter::PlotSubBlock (const crs::TSubBlock &SubBlock) {
  
  std::ostringstream oFileName;
  
  switch (SubBlock.GetBlockType ()) {
    
  case crs::TSubBlock::eRUNH:
    break;
    
  case crs::TSubBlock::eEVTH:
    oFileName << fFileName
	      << "_"
	      << ((crs::MEventHeader)SubBlock).GetEventNumber ()
	      << ".root";
    fFile = TFile::Open (oFileName.str ().c_str (), "RECREATE");
    std::cout << " write to: " << oFileName.str () << std::endl;
    fEnergy = new TH1D ("energy", "log(e) for EM", 200, -3.2, -2.); 
    break;
    
  case crs::TSubBlock::ePARTDATA:
    PlotParticles (SubBlock);
    break;
    
  case crs::TSubBlock::eLONG:
    break;
    
  case crs::TSubBlock::eEVTE:
    fFile->cd ();
    fEnergy->Write ();
    fFile->Close ();	
    delete fFile;
    break;
    
  case crs::TSubBlock::eRUNE:
    break;
    
  case crs::TSubBlock::eNODATA:
    break;
    
  }
}



void TCorsikaPlotter::PlotParticles (const crs::MParticleBlock &Data) {
  
  crs::MParticleBlock::ParticleListConstIterator iEntry;
  
  for (iEntry = Data.FirstParticle ();
       iEntry != Data.LastParticle ();
       ++iEntry) {
    
    if (iEntry->IsParticle ()) {
      
      crs::MParticle Particle (*iEntry);
      
      //	    if (Particle.GetParticleID ()<4) {
      if (Particle.GetParticleID ()==3 || Particle.GetParticleID ()==2) {
	fEnergy->Fill (log10 (Particle.GetKinEnergy ()), 
		       Particle.GetWeight ());
      }
    }
  }
}
