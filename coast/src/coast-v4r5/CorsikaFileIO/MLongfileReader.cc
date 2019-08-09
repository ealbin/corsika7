#include <crsRead/MLongfileReader.h>
#include <crsRead/TSubBlockIO.h>
using namespace crsRead;

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <fstream>
using namespace std;

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


MLongfileReader::MLongfileReader(const std::string &filename, 
                                 int verbose) 
  : fFileName(filename),
    fVerbose(verbose) {
  fFile = new std::ifstream(fFileName.c_str(), std::ios::in);
}


bool MLongfileReader::ReadShower(int n) {
  
  // return to start
  fFile->seekg(0 , std::ios::beg);

  if (!FindShower(n)) return false;
  if (!ParseParticleTable(n)) return false;
  return ParseEdepTable(n);
}


bool MLongfileReader::FindShower(int n) {
  
  while (fFile->good()) {
    
    string LONGITUDINAL, DISTRIBUTION;
    (*fFile) >> LONGITUDINAL;
    if (LONGITUDINAL != "LONGITUDINAL") {
      fFile->ignore(9999, '\n'); // next line
      continue;
    }
    
    (*fFile) >> DISTRIBUTION;
    if (DISTRIBUTION != "DISTRIBUTION") {
      fFile->ignore(9999, '\n'); // next line
      continue;
    }
    
    string IN, STEPS, OF, UNITS, FOR, SHOWER;
    string mode;
    int nSteps, Id;
    double deltaX;
    (*fFile) >> /*LONGITUDINAL >> DISTRIBUTION >>*/ IN >> nSteps >> mode
          >> STEPS >> OF >> deltaX >> UNITS >> FOR >> SHOWER >> Id;
    fFile->ignore(99999, '\n');  // skip first line
    
    if (Id != n) { continue; }
    
    fNBins = nSteps;
    fDeltaX = deltaX;
    
    if (mode == "VERTICAL") {
      fIsSlant = false;
    } else if (mode == "SLANT") {
      fIsSlant = true;
    } else {
      cerr << " MLongfileReader> Shower mode \"" 
           << mode
           << "\" unknown in long file! "
           << endl;
      return false;
    }
    
    return true;
    
  }
  
  return false;
}



bool MLongfileReader::ParseParticleTable(int /*n*/) {
  
  //DEPTH     GAMMAS   POSITRONS   ELECTRONS         MU+         MU-     HADRONS     CHARGED      NUCLEI   CHERENKOV
  fFile->ignore(99999, '\n');  // skip line
  
  fParticleEntry.clear();

  // loop data
  for (int i=0; i<fNBins; ++i) {
    double depth, gammas, positrons, electrons, mu, muBar, 
      had, charged, nuclei, ckov;
    
    (*fFile) >> depth >> gammas >> positrons >> electrons >> mu >> muBar 
          >> had >> charged >> nuclei >> ckov;
    fFile->ignore(999, '\n');
    
    fParticleEntry.push_back(ParticleEntry(depth,
                                           gammas,
                                           positrons,
                                           electrons,
                                           mu,
                                           muBar ,
                                           had,
                                           charged,
                                           nuclei,
                                           ckov));
    
  }
  return true;
}


bool MLongfileReader::ParseEdepTable(int n) {
  
  string LONGITUDINAL, ENERGY, DEPOSIT, IN, STEPS, OF, UNIT, 
    FOR, SHOWER;
  string mode;
  int nSteps, id;
  double deltaX;
  
  (*fFile) >> LONGITUDINAL >> ENERGY >> DEPOSIT >> IN >> nSteps >> mode >> STEPS 
        >> OF >> deltaX >> UNIT >> FOR >> SHOWER >> id;
  fFile->ignore(99999, '\n');  // skip rest of line
  
  //DEPTH       GAMMA    EM IONIZ     EM CUT    MU IONIZ      MU CUT  HADR IONIZ    HADR CUT   NEUTRINO        SUM
  fFile->ignore(99999, '\n');  // skip line
  
  if (id != n) {
    cerr << "MLongfileReader: Error 1" << endl;
    return false;
  }

  fNBinsEdep = nSteps;
  fDeltaXEdep = deltaX;
  
  if (mode == "VERTICAL") {
    fIsSlantEdep = false;
  } else if (mode == "SLANT") {
    fIsSlantEdep = true;
  } else {
    cerr << " MLongfileReader> Shower mode \"" 
         << mode
         << "\" unknown in long file (edep)! "
         << endl;
    return false;
  }
  
  fEdepEntry.clear();
  
  // loop data
  for (int i=0; i<nSteps; ++i) {
    
    double depth,  gamma, emIoniz, emCut, muIoniz, muCut, hadIoniz, hadCut,
      neutrino, sum;
    
    
    (*fFile) >> depth >> gamma >> emIoniz >> emCut >> muIoniz >> muCut 
          >> hadIoniz >> hadCut >> neutrino >> sum;
    fFile->ignore(999, '\n');

    fEdepEntry.push_back(EdepEntry(depth,
                                   gamma,
                                   emIoniz,
                                   emCut,
                                   muIoniz,
                                   muCut,
                                   hadIoniz,
                                   hadCut,
                                   neutrino,
                                   sum));
    
  }
  return true;
}

