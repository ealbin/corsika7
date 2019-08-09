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
}


TReaderExample::~TReaderExample () {
}




void TReaderExample::HandleSubBlock (const crs::TSubBlock &sb) {

    static int CountSubBlock = 0;

    if (CountSubBlock<10) {
	cout << "TReaderExample::HandleSubBlock # " 
	     << CountSubBlock
	     << " type: " 
	     << sb.GetBlockType ()
	     << endl;
    } if (CountSubBlock == 10 ) {
	cout << "TReaderExample::HandleSubBlock -> "
	     << "stop showing individual SubBlocks ... too many"
	     << endl;
    }

    CountSubBlock ++;
    
    TCorsikaReader::HandleSubBlock (sb);
}


void TReaderExample::HandleParticleBlock (const crs::MParticleBlock &p) {

    static int CountParticleBlock = 0;

    if (CountParticleBlock<10) {
	cout << "TReaderExample::HandleParticleBlock # " 
	     << CountParticleBlock
	     << endl;
    } if (CountParticleBlock == 10) {
	cout << "TReaderExample::HandleParticleBlock -> "
	     << "stop showing individual ParticleBlocks ... too many"
	     << endl;
    }

    CountParticleBlock ++;
    
    TCorsikaReader::HandleParticleBlock (p);
}



void TReaderExample::HandleLongitudinalBlock(const crs::MLongitudinalBlock &s){

    static int CountLongitudinalBlock = 0;

    if (CountLongitudinalBlock<10) {
	cout << "TReaderExample::HandleLongitudinalBlock # " 
	     << CountLongitudinalBlock
	     << endl;
    } if (CountLongitudinalBlock == 10) {
	cout << "TReaderExample::HandleLongitudinalBlock -> "
	     << "stop showing individual LongitudinalBlocks ... too many"
	     << endl;
    }

    CountLongitudinalBlock ++;

    TCorsikaReader::HandleLongitudinalBlock (s);
}


void TReaderExample::HandleEventStart (const crs::MEventHeader &sb) {

    cout << "TReaderExample::EventStart EventNo "
	 << sb.GetEventNumber ()
	 << endl;
}


void TReaderExample::HandleEventEnd (const crs::MEventEnd &sb) {

    cout << "TReaderExample::EventEnd"
	 << endl;
}

void TReaderExample::HandleRunStart (const crs::MRunHeader &sb) {

    cout << "TReaderExample::RunStart"
	 << endl;
}

void TReaderExample::HandleRunEnd (const crs::MRunEnd &sb) {

    cout << "TReaderExample::RunEnd"
	 << endl;
}

void TReaderExample::HandleParticle (const crs::MParticle &p) {

    static int CountParticle = 0;

    if (CountParticle<20) {
	cout << "TReaderExample::Particle # "
	     << CountParticle
	     << " > " << p
	     << endl;
    } else if (CountParticle == 20 ) {
	cout << "TReaderExample::Particle -> "
	     << "stop showing individual Particles ... too many"
	     << endl;
    }
    
    CountParticle ++;
}


void TReaderExample::HandleCherenkov (const crs::MCherenkov &c) {

    static int CountCherenkov = 0;

    if (CountCherenkov<20) {
	cout << "TReaderExample::Cherenkov # "
	     << CountCherenkov
	     << " > " << c
	     << endl;
    } else if (CountCherenkov == 20 ) {
	cout << "TReaderExample::Cherenkov -> "
	     << "stop showing individual Cherenkovs ... too many"
	     << endl;
    }
    
    CountCherenkov ++;
}

void TReaderExample::HandleMuonProductionInfo (const crs::MMuonProductionInfo &m) {

    static int CountMuonProduction = 0;

    if (CountMuonProduction<20) {
	cout << "TReaderExample::MuonProduction # "
	     << CountMuonProduction
	     << " > " << m
	     << endl;
    } else if (CountMuonProduction == 20 ) {
	cout << "TReaderExample::MuonProduction -> "
	     << "stop showing individual MuonProductions ... too many"
	     << endl;
    }
    
    CountMuonProduction ++;
}
	
// virtual user-event handlers
void TReaderExample::Init () {

    cout << "TReaderExample::Init"
	 << endl;
}


void TReaderExample::Exit () {

    cout << "TReaderExample::Exit"
	 << endl;
}


