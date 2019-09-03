/* $Id: CParticle.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/CParticle.h>
#include <crs/CorsikaConsts.h>
#include <crs/MEventHeader.h>

#include <iostream>
#include <cstdlib>

using namespace crs;
using namespace std;


void 
CParticle::Dump () 
  const 
{
  cout << " CParticle> Id: " << particleId 
       << ", hadGen: " << hadronicGeneration
       << ", E: " << energy
       << ", x: " << x 
       << ", y: " << y 
       << ", z: " << z
       << ", X: " << depth
       << ", t: " << time  
       << ", w: " << weight
       << endl;
}


/*  nice try ...
MParticle TransformParticle(const crs::MEventHeader& header) const {
  
  CREAL* newParticleData = 
  
  return particle(newParticleData, header.IsThinned());
}
*/


/*
crs::GroundCoordinates 
CParticle::TransformToGroundCoordinates(const double cosTheta,
					const double cosPhiX,
					const double cosPhiY,
					const double arrayRotation) 
  const 
{
  const double ETOT = energy;
  const double PAMA = crs::gParticleMass[GetParticleId()];
  const double PTOT = sqrt( (ETOT-PAMA)*(ETOT+PAMA) ) ;
  const double STT  = sqrt( (1-cosTheta)*(1+cosTheta) );
  double PHIPAR = 0;
  if ( cosPhiY != 0  ||  cosPhiX != 0 ) {
    PHIPAR = atan2( cosPhiY, cosPhiX ); 
  } 
  
  const double COSANG = cos(arrayRotation);
  const double SINANG = sin(arrayRotation);
  
  return crs::GroundCoordinates(PTOT * STT * cos( PHIPAR + arrayRotation ), // Px
				PTOT * STT * sin( PHIPAR + arrayRotation ), // Py
				PTOT * cosTheta,                     // Pz
				x * COSANG + y * SINANG,             // x 
				y * COSANG - x * SINANG);            // y
}
*/
