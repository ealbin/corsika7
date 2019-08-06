#include <TReaderExample.h>

#include <crs/TBlock.h>
#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MRunEnd.h>
#include <crs/MEventHeader.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>
#include <crs/MCherenkov.h>
#include <crs/MMuonProductionInfo.h>

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
using namespace std;


TReaderExample::TReaderExample (const std::string &filename, 
                                int verbose) :
    TCorsikaReader (filename, verbose) {

    fCountSubBlock = 0;
    fCountParticleBlock = 0;
    fCountLongitudinalBlock = 0;
    fCountParticle = 0;
    fCountCherenkov = 0;
    fCountMuonProduction = 0;
}


TReaderExample::~TReaderExample () {
}




void 
TReaderExample::HandleSubBlock (const crs::TSubBlock &sb) 
{
  if (fCountSubBlockIds.count(sb.GetBlockType())==0)
    fCountSubBlockIds[sb.GetBlockType()] = 0;
  fCountSubBlockIds[sb.GetBlockType()]++;
  
    if (fCountSubBlock<10) {
	cout << "TReaderExample::HandleSubBlock # " 
	     << fCountSubBlock
	     << " type: " 
	     << sb.GetBlockType ()
             << " ("
	     << sb.GetBlockTypeName ()
             << ")"
	     << endl;
    } if (fCountSubBlock == 10 ) {
	cout << "TReaderExample::HandleSubBlock -> "
	     << "stop showing individual SubBlocks ... too many"
	     << endl;
    }

    fCountSubBlock ++;
    
    TCorsikaReader::HandleSubBlock (sb);
}


void
TReaderExample::HandleParticleBlock (const crs::MParticleBlock &p) 
{
    if (fCountParticleBlock<10) {
	cout << "TReaderExample::HandleParticleBlock # " 
	     << fCountParticleBlock
	     << endl;
    } if (fCountParticleBlock == 10) {
	cout << "TReaderExample::HandleParticleBlock -> "
	     << "stop showing individual ParticleBlocks ... too many"
	     << endl;
    }

    fCountParticleBlock ++;
    
    TCorsikaReader::HandleParticleBlock (p);
}



void
TReaderExample::HandleLongitudinalBlock(const crs::MLongitudinalBlock &s)
{
    if (fCountLongitudinalBlock<10) {
	cout << "TReaderExample::HandleLongitudinalBlock # " 
	     << fCountLongitudinalBlock
	     << endl;
    } if (fCountLongitudinalBlock == 10) {
	cout << "TReaderExample::HandleLongitudinalBlock -> "
	     << "stop showing individual LongitudinalBlocks ... too many"
	     << endl;
    }

    fCountLongitudinalBlock ++;

    TCorsikaReader::HandleLongitudinalBlock (s);
}


void
TReaderExample::HandleEventStart (const crs::MEventHeader &sb) 
{
    cout << "TReaderExample::EventStart EventNo "
	 << sb.GetEventNumber ()
	 << endl;
}


void
TReaderExample::HandleEventEnd (const crs::MEventEnd &sb) 
{
    cout << "TReaderExample::EventEnd"
	 << endl;
}

void
TReaderExample::HandleRunStart (const crs::MRunHeader &sb) 
{
    cout << "TReaderExample::RunStart"
	 << endl;
}

void
TReaderExample::HandleRunEnd (const crs::MRunEnd &sb) 
{
    cout << "TReaderExample::RunEnd"
	 << endl;
}

void
TReaderExample::HandleParticle (const crs::MParticle &p) 
{
    if (fCountParticle<20) {
	cout << "TReaderExample::Particle # "
	     << fCountParticle
	     << " > " << p
	     << endl;
    } else if (fCountParticle == 20 ) {
	cout << "TReaderExample::Particle -> "
	     << "stop showing individual Particles ... too many"
	     << endl;
    }

    const int pId = p.GetParticleId();
    if (fCountParticleIds.count(pId)==0)
      fCountParticleIds[pId] = 0;
    fCountParticleIds[pId]++;
    
    fCountParticle++;
}


void 
TReaderExample::HandleCherenkov (const crs::MCherenkov &c) 
{
    if (fCountCherenkov<20) {
	cout << "TReaderExample::Cherenkov # "
	     << fCountCherenkov
	     << " > " << c
	     << endl;
    } else if (fCountCherenkov == 20 ) {
	cout << "TReaderExample::Cherenkov -> "
	     << "stop showing individual Cherenkovs ... too many"
	     << endl;
    }
    
    fCountCherenkov ++;
}

void
TReaderExample::HandleMuonProductionInfo (const crs::MMuonProductionInfo &m) 
{
    if (fCountMuonProduction<20) {
	cout << "TReaderExample::MuonProduction # "
	     << fCountMuonProduction
	     << " > " << m
	     << endl;
    } else if (fCountMuonProduction == 20 ) {
	cout << "TReaderExample::MuonProduction -> "
	     << "stop showing individual MuonProductions ... too many"
	     << endl;
    }
    
    fCountMuonProduction ++;
}
	
// virtual user-event handlers
void
TReaderExample::Init () 
{
    cout << "TReaderExample::Init"
	 << endl;
}


void
TReaderExample::Exit () 
{
    cout << "TReaderExample::Exit"
	 << endl;
}


void 
TReaderExample::Summary() const
{
  cout << "\n===================================================\n"
       << " subblocks:    " << fCountSubBlock << "\n"
       << " partcl-blocks:" << fCountParticleBlock << "\n";
  for (map<int,int>::const_iterator i=fCountSubBlockIds.begin(); i!=fCountSubBlockIds.end(); ++i)
    cout << "       id=" << setw(5) << i->first << " (" << setw(14) << crs::TSubBlock::GetBlockTypeName((crs::TSubBlock::SubBlockType)i->first) << ") counted " << i->second << "\n"; 
  cout << " longi-blocks:  " << fCountLongitudinalBlock << "\n"
       << " particles:     " << fCountParticle << "\n";
  for (map<int,int>::const_iterator i=fCountParticleIds.begin(); i!=fCountParticleIds.end(); ++i)
    cout << "       id=" << setw(6) << i->first << " counted " << i->second << "\n"; 
  cout << " cher-photons:  " << fCountCherenkov << "\n"
       << " muon-prods:    " << fCountMuonProduction << "\n"
       << endl;
}

