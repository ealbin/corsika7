/* $Id: coast.cc 6960 2018-12-20 16:38:09Z pierog $   */

#include <interface/CorsikaInterface.h>

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>

#include <sstream>
#include <iostream>

using namespace std;

// this is just a trick for this example. It can be done better in 
// real user code
bool gIsThinned__ = false;


/*
  Data is one CORSIKA data-block constining of 21 SubBlocks.
  A SubBlock can be:
  - thinned mode:     39 (Particles) * 8 (ENTRIES) * 4 (BYTES) 
  - not-thinned mode: 39 (Particles) * 7 (ENTRIES) * 4 (BYTES) 
*/
void
wrida_ (const CREAL *Data) 
{    
  crs::CParticleFortranPtr p;
  const bool isF = prminfo_(p); 
  std::cout << "coast wrida_ isF=" << isF << std::endl;
  p->Dump();
  std::cout << std::endl;


  crs::TSubBlock subBlock (Data, gIsThinned__);    
}




void
inida_ (const char* filename,
        const int& thinning,
        const int& /*curved*/,
        const int& /*slant*/,
        const int& /*stackinput*/,
        const int& /*preshower*/,
        int str_length) 
{
  gIsThinned__ = thinning;
  std::cout << " ######################################################### "
	    << std::endl
	    << " This is the COAST example init function. \n"
	    << " CORSIKA output: " << filename
	    << " thinned: " << thinning
	    << std::endl;
    
  crs::CParticleFortranPtr pptr;
  const bool isF = prminfo_(pptr);
  std::cout << "coast inida_ isF=" << isF << std::endl;
  pptr->Dump();
  std::cout << std::endl;
}

void
cloda_ () 
{
  std::cout << " ######################################################### "
	    << std::endl
	    << " COAST example closing "
	    << std::endl;

  crs::CParticleFortranPtr pptr;
  const bool isF = prminfo_(pptr); 
  std::cout << "coast cloda_ isF=" << isF << std::endl;
  pptr->Dump();
  std::cout << std::endl;
}


void 
interaction_ (const crs::CInteraction& interaction)
{
  /*
    all interactions in the shower are available in this function ! 
    
    the information availabel in the CInteraction class are:
    
    double x;
    double y;
    double z;      
    double etot;      // lab energy
    double sigma;     // cross-section of process
    double kela;      // elasticity
    int    projId;    // projectile
    int    targetId;  // target
    double time;
  */
}


void
track_ (const crs::CParticle& pre, const crs::CParticle& post)
{
  /*
    all particles in the shower are available in this function !
    
    The pre and post objecte are the two endpoints for one single track 
    in the shower, where the information available in CParticle is:

    double x;
    double y;
    double z;
    double depth;
    double time;
    double energy;
    double weight;
    int    particleId;
    int    hadronicGeneration;

  */
} 


void tabularizedatmosphere_(const int &nPoints, const double* height, const double* refractiveIndex)
{
	// for special use only but should be defined because it is delcared in CORSIKA.F
}
