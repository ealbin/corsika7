/* $Id: MParticleSubBlockOutput.cc,v 1.1 2007-10-19 11:13:40 rulrich Exp $   */


#include <crs/MEventHeader.h>

#include <crsRead/MParticleSubBlockOutput.h>

using namespace crsRead; 

/*
  MParticleSubBlockOutput::MParticleSubBlockOutput () 
  : TSubBlockIO(),
  fParticleWriteIndex(0) {
  }
*/


MParticleSubBlockOutput::MParticleSubBlockOutput (bool thinned, int AdditionalCharacters, 
						  const crs::MEventHeader& header) 
: TSubBlockIO(thinned, AdditionalCharacters),
  fParticleWriteIndex(0),
  fArrayRotation(header.GetArrayRotation())        // <- this is for particle coordinate transformations
{                
  for (int i=0; i<fBlockSize; ++i) 
    ((CREAL*)fSubBlockData)[i] = 0;   
}

void
MParticleSubBlockOutput::Clear() 
{
  for (int i=0; i<fBlockSize; ++i) 
    ((CREAL*)fSubBlockData)[i] = 0; 
  fParticleWriteIndex = 0;
}
