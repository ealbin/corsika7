/* $Id: TC2R.cc,v 1.3 2007-10-18 14:33:47 pierog Exp $   */

#include <crs2r/TC2R.h>
using namespace crs2r;

#include <crs/TBlock.h>
#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/CParticle.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MEventEnd.h>
#include <crs/MRunEnd.h>
#include <crs/MParticle.h>
#include <crs/MCherenkov.h>
#include <crs/MMuonProductionInfo.h>

#include <crsIO/TRun.h>
#include <crsIO/TShower.h>
#include <crsIO/TParticle.h>
#include <crsIO/TCherenkov.h>
#include <crsIO/TLongitudinal.h>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObjectTable.h>

#include <iostream>



/*
 */
TC2R::TC2R() :

fFile(0),
  fShower(0),
  fRun(0),
  fCurrentRun(0),
  fCurrentShower(0),
  fParticles(0),
  fCherenkov(0),
  fLongitudinal(0),
  fNParticles(0),
  fNLong(0),
  fNCherenkov(0) {
}


TC2R::~TC2R() {

  delete fCurrentRun;
  delete fCurrentShower;
  delete fFile;
}


void TC2R::Open(std::string filename,
		 bool thinned) {

  fFile = TFile::Open(filename.c_str(), "RECREATE");
  fFile->cd();
  fThinned = thinned;    

  fCurrentShower = new crsIO::TShower();
  fShower = new TTree("sim", "simulated showers");
  fShower->Branch("shower.", "crsIO::TShower", &fCurrentShower);

  fCurrentRun = new crsIO::TRun();
  fRun = new TTree("run", "simulation run");
  fRun->Branch("run.", "crsIO::TRun", &fCurrentRun);
    
  fParticles = new TClonesArray("crsIO::TParticle", 1000);
  fShower->Branch("particle.", &fParticles);
  fNParticles = 0;
    
  fLongitudinal = new TClonesArray("crsIO::TLongitudinal", 100000);
  fShower->Branch("long.", &fLongitudinal);
  fNLong = 0;
    
  fCherenkov = new TClonesArray("crsIO::TCherenkov", 100000);
  fShower->Branch("cherenkov.", &fCherenkov);
  fNCherenkov = 0;
}


void TC2R::Write(const CREAL *data) {
    
  crs::TBlock block(data, fThinned);
  Write(block);

}

void TC2R::Write(const crs::TBlock &Block) {

  crs::TBlock::SubBlockConstIterator iSubBlock ;

  for(iSubBlock = Block.GetFirstSubBlock(); 
       iSubBlock != Block.GetLastSubBlock(); 
       ++iSubBlock) {

    Write(*iSubBlock);
  }
}

void TC2R::Write(const crs::TSubBlock &SubBlock) {

  if(fFile) {

    switch(SubBlock.GetBlockType()) {
	    
    case crs::TSubBlock::eRUNH:
      fCurrentRun->Clear();
      fCurrentShower->Clear();
      fCurrentRun->AddRunHeader(SubBlock);
      fCurrentShower->AddRunHeader(SubBlock);
      break;
		
    case crs::TSubBlock::eEVTH:
      fCurrentRun->AddEventHeader(SubBlock);
      fCurrentShower->AddEventHeader(SubBlock);
      break;
		
    case crs::TSubBlock::ePARTDATA:
      AddParticleBlock(SubBlock);
      break;

    case crs::TSubBlock::eLONG:
      AddLongitudinalBlock(SubBlock);
      break;
		
    case crs::TSubBlock::eEVTE:
      fCurrentShower->AddEventEnd(SubBlock);
		
      fShower->Fill();
		
      fNParticles = 0;
      fParticles->Clear();
				
      fNLong = 0;
      fLongitudinal->Clear();
		
      fNCherenkov = 0;
      fCherenkov->Clear();

      //gObjectTable->Dump();
      break;
		
    case crs::TSubBlock::eRUNE:
      fRun->Fill();
      //fCurrentRun->Write();
      //fCurrentShower->Wr
      // do nothing
      break;

    case crs::TSubBlock::eNODATA:
      break;
    }

  }

}

/*
  For the inclined shower front sampling
 */
void TC2R::Write(const crs::MEventHeader &header, 
		 const crs::GroundParticle& particle) {
  new((*fParticles)[fNParticles++]) crsIO::TParticle(header, particle);
}

void TC2R::AddParticleBlock(const crs::MParticleBlock &Data) {
    
  crs::MParticleBlock::ParticleListConstIterator iEntry;
  for(iEntry = Data.FirstParticle();
       iEntry != Data.LastParticle();
       ++iEntry) {
	
    switch(iEntry->GetType()) {

    case crs::eParticle:
    case crs::eNucleus:
      new((*fParticles) [fNParticles++]) 
	crsIO::TParticle(*iEntry);
      break;
		
    case crs::eCherenkov:
      new((*fCherenkov) [fNCherenkov++]) 
	crsIO::TCherenkov(*iEntry);
      break;
		
    case crs::eMuonProductionInfo:
      new((*fParticles) [fNParticles++]) 
	crsIO::TParticle(*iEntry);
      break;
		
    case crs::eEmpty:
    case crs::eUnknown:
      break;
		
    }
  }
}



void TC2R::AddLongitudinalBlock(const crs::MLongitudinalBlock &Data) {
    
  crs::MLongitudinalBlock::LongitudinalListConstIterator iLong;

  for(iLong = Data.FirstLongitudinal();
       iLong != Data.LastLongitudinal();
       ++iLong) {
	
    new((*fLongitudinal) [fNLong++]) crsIO::TLongitudinal(*iLong);
	
  }
}



void TC2R::Close() {

  // fShower->Write();
  // fRun->Write();

  if(fFile) {
    fFile->Write();
    fFile->Close();
  }
  delete fFile;
  fFile = 0;
}





