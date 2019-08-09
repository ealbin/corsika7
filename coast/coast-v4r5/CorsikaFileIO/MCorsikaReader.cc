#include <crsRead/MCorsikaReader.h>
#include <crsRead/TSubBlockIO.h>
using namespace crsRead;

#include <string>
#include <iostream>
#include <fstream>
#include <memory>

#include <crs/TBlock.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>
#include <crs/MRunEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>
#include <crs/MMuonProductionInfo.h>
#include <crs/MCherenkov.h>
#include <crs/TLongitudinal.h>


MCorsikaReader::MCorsikaReader (const std::string &filename, 
				int verbose) 
  : TCorsikaReader (filename, verbose) {
}



bool MCorsikaReader::GetRun (crs::MRunHeader &run) {
  
  do {
    
    if (!Read ()) {
      return false;
    }
    
    run = fSubBlock;
    
  } while (run.GetBlockType ()!=crs::TSubBlock::eRUNH);
  
  return true;
}



bool MCorsikaReader::GetShower (crs::MEventHeader &shower) {
  
  do {
    
    if (!Read ()) {
      return false;
    }
    
    shower = fSubBlock;
    
  } while (shower.GetBlockType ()!=crs::TSubBlock::eEVTH);
  
  return true;
}



bool MCorsikaReader::GetData (crs::TSubBlock &data) {
  
  if (!Read ()) {
    return false;
  }
  
  data = fSubBlock;
  
  if (data.GetBlockType ()==crs::TSubBlock::eEVTE) 
    return false;
  
  return true;
}



bool MCorsikaReader::GetShowerSummary (crs::MEventEnd &summary) {
  
  summary = (fSubBlock);
  
  while (summary.GetBlockType ()!=crs::TSubBlock::eEVTE) {
    
    if (!Read ()) {
      return false;
    }
    
    summary = fSubBlock;
    
  }
  
  return true;
}



// this is needed to prevent the Read function from reading the full file
// at once.
void 
MCorsikaReader::HandleSubBlock (const crs::TSubBlock& /*CrsSubBlock*/) 
{  
  InterruptReader ();
}
