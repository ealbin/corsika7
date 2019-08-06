/**
   \file MParticleSubBlockOutput.h                                      

   Enables the writing of individual particles into CORSIKA subblocks
   
   \author Ralf Ulrich 
   \date Fri Oct 19 10:43:17 CEST 2007
   \version $Id: MParticleSubBlockOutput.h 5116 2016-01-04 19:09:04Z darko $ 

*/


#ifndef _INCLUDE_MParticleSubBlockOutput_H_
#define _INCLUDE_MParticleSubBlockOutput_H_

#include <crsRead/TSubBlockIO.h>

#include <crs/CParticle.h>
#include <crs/CorsikaTypes.h>
#include <crs/TBlock.h>

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

namespace crs {
  class GroundParticle;
  class MEventHeader;
}


namespace crsRead {

  /**
     
     \class MParticleSubBlockOutput
     \brief Class to write individual particles into subblocks
     
  */


  class MParticleSubBlockOutput : public TSubBlockIO {
    
    // friend inline std::ofstream& operator<<(std::ofstream& ,MParticleSubBlockOutput&);
    friend inline MParticleSubBlockOutput& operator<<(MParticleSubBlockOutput&, const crs::GroundParticle& p);
    
  public:
    // MParticleSubBlockOutput();                                ///< Default constructor
    MParticleSubBlockOutput(bool Thinned, int NBlockSizeInfo, const crs::MEventHeader& header);
    
    bool IsFull() const {return fParticleWriteIndex>=crs::TBlock::fgNParticles;}
    void Clear();
    
  private:
    MParticleSubBlockOutput(const MParticleSubBlockOutput&);                ///< copy constructor
    MParticleSubBlockOutput& operator=(const MParticleSubBlockOutput &sb);
    
  protected:
    int fParticleWriteIndex;
    
    const double fArrayRotation;
  };
  
  
  
  /**
     \fn inline ofstream& operator<<( ofstream& s, MParticleSubBlockOutput& b)
     \brief Corsika data file streaming operator

     This function defines the operator<< to stream out CORSIKA/COAST CParticles into 
     binary particle sub-blocks.
       
     \author Ralf Ulrich
     \date Tue Jun 14 17:45:10 CEST 2005
     \version $Id: MParticleSubBlockOutput.h 5116 2016-01-04 19:09:04Z darko $
  */
  inline MParticleSubBlockOutput& 
  operator<<(MParticleSubBlockOutput& b, const crs::GroundParticle& p) 
  {  
    if (b.IsFull()) {
      std::cerr << "COAST> MParticleSubBlockOutput::ERROR trying to add more than "
		<< crs::TBlock::fgNParticles 
		<< " to particle subblock! Check your code! " 
		<< std::endl;
      return b;
    }
    
    const int writeIndex = b.fParticleWriteIndex 
      * (b.fThinned ? crs::TBlock::fgNEntriesThinned : crs::TBlock::fgNEntriesNotThinned );
    
    CREAL* particleData = &(((CREAL*)b.fSubBlockData)[writeIndex]);
    
    const double ARRANR = b.fArrayRotation;
    const double COSANG = std::cos(ARRANR);
    const double SINANG = std::sin(ARRANR);   	     
    
    particleData[0] = p.GetCorsikaId();
    particleData[1] = p.Px * COSANG + p.Py * SINANG;
    particleData[2] = p.Py * COSANG - p.Px * SINANG;
    particleData[3] = p.Pz;
    particleData[4] = p.x * COSANG + p.y * SINANG;
    particleData[5] = p.y * COSANG - p.x * SINANG;
    particleData[6] = p.time;
    if (b.fThinned) {
      particleData[7] = p.weight;
    }

    ++b.fParticleWriteIndex;
    
    return b;
  }
    
};


#endif
