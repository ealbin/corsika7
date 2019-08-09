/**
   \file TSubBlockIO.h                                      

   Corsika File Format utilities.                           
   Basic input of CORSIKA data from files.

   \author Ralf Ulrich 
   \date Fri Jun 18 14:06:43 EDT 2004                 
   \version $Id: TSubBlockIO.h,v 1.2 2007-10-19 11:13:40 rulrich Exp $ 

*/


#ifndef TCORSIKASUBBLOCKIO_H
#define TCORSIKASUBBLOCKIO_H

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/TBlock.h>

#include <fstream>
#include <iostream>
#include <iomanip>


namespace crsRead {

  /**

     \class TSubBlockIO
     \brief This class defines one subblock of a CORSIKA binary data file
 
     The purpose of this class is to count SubBlocks, know about the Block 
     structure of the bianry file and the FORTRAN steering bytes added here
     and there. There are streaming operators >> and << implemented for 
     convinient read and write operations.
     The idea is to hide the CORSIKA(fortran) block structure from the user and 
     concentrate on the data containing subblocks.  
 
     One SubBlock contains 
     - 39*7*(REAL*4) = 1092 Bytes in the not-thinned mode and 
     - 39*8*(REAL*4) = 1248 Bytes in thinned mode.

     One Block contains 21 SubBlocks
     - 21*1092 = 22932 Bytes in the not-thinned mode and 
     - 21*1248 = 26208 Bytes in thinned mode.
 
 
     before each Block fortran puts the block length in bytes, BUT in a machine 
     dependent format 
  */
  
  class MParticleSubBlockOutput;

  class TSubBlockIO : public crs::TSubBlock {

    friend inline std::ifstream& operator>>(std::ifstream& ,TSubBlockIO&);
    friend inline std::ofstream& operator<<(std::ofstream& ,TSubBlockIO&);

  public:
    TSubBlockIO();                                ///< Default constructor
    TSubBlockIO(bool Thinned, int NBlockSizeInfo);
    ~TSubBlockIO(); 
    
    void InitSubBlockSize(bool thinn, int NBlockSizeInfo);
    
    void SetSubBlockCounter(int n) {fSubBlockCounter = n;}
    int GetSubBlockCounter() const {return fSubBlockCounter;}
    
    int FinishBlock(std::ofstream& s);

    int FindRUNH();
    void SetBuffer(const CREAL* buffer);
    void SetBuffer(const crs::TSubBlock& sb);
	
  private:
    TSubBlockIO(const TSubBlockIO&);                ///< copy constructor
    TSubBlockIO& operator=(const TSubBlockIO &sb);
    void InitBuffer();
	
  protected:
    int    fBlockSize;            ///< Actual size of one CORSIKA SubBlock
    int    fNBlockSizeInfo;       ///< Number of Additional bytess after 
    ///                                each Block(21 SubBlocks)
    int    fSubBlockCounter;      ///< Counter of SubBlocks read from file
    char  *fBlockSizeInfo;        ///< FORTRAN block size info
	
    bool fOwnBuffer;      ///< flag for ownership over buffer
  };
  
    
  
  /**
     \fn inline ifstream& operator>>( ifstream& s, TSubBlockIO& b)
     \brief Corsika data file streaming operator

     This function defines the operator>> to stream in CORSIKA data into 
     Corsika sub-blocks
       
     \author Ralf Ulrich
     \date 05.2002
     \version $Id: TSubBlockIO.h,v 1.2 2007-10-19 11:13:40 rulrich Exp $
  */
  inline 
  std::ifstream& 
  operator>>(std::ifstream& s, TSubBlockIO& b) 
  {
    // BLOCK begin
    if( !((b.fSubBlockCounter)%crs::TBlock::fgNSubBlocks) ) {
      s.read(b.fBlockSizeInfo, b.fNBlockSizeInfo);
    }
    
    
    s.read((char*)(b.fSubBlockData), b.fBlockSize*4);
    b.FindBlockType();
		

    // BLOCK end
    if( !((b.fSubBlockCounter+1)%crs::TBlock::fgNSubBlocks) ) {
      s.read(b.fBlockSizeInfo, b.fNBlockSizeInfo);
    }
	
    ++b.fSubBlockCounter;
	
    return s;
  }




  /**
     \fn inline ofstream& operator<<( ofstream& s, TSubBlockIO& b)
     \brief Corsika data file streaming operator

     This function defines the operator<< to stream out CORSIKA sub-blocks
     into binary files.
       
     \author Ralf Ulrich
     \date Tue Jun 14 17:45:10 CEST 2005
     \version $Id: TSubBlockIO.h,v 1.2 2007-10-19 11:13:40 rulrich Exp $
  */
  inline 
  std::ofstream& 
  operator<<(std::ofstream& s, TSubBlockIO& b) 
  {
    // BLOCK begin
    if( !((b.fSubBlockCounter)%crs::TBlock::fgNSubBlocks) ) {
      s.write(b.fBlockSizeInfo, b.fNBlockSizeInfo);
      // -------- debug -------------
      //std::cout << "COAST> ~~~~~~~ writing " << b.fNBlockSizeInfo << std::endl;
      // -------- debug -------------
    }

    
    s.write((char*)(b.fSubBlockData), b.fBlockSize*4);
	
    // -------- debug -------------
    //for (int i=0; i<b.fBlockSize*4; ++i) {
    //  std::cout << " sb-io: " << std::setw(5) << i << " \'"
    //<< *(((char*)b.fSubBlockData)+i) << "\' " 
    //<< std::setw(4) << (int)(*(((char*)b.fSubBlockData)+i));
    //if (!(i%4))
    //std::cout << " " << (*(((CREAL*)b.fSubBlockData)+i/4));
    //std::cout << std::endl;
    //}
    // -------- debug -------------


    // BLOCK ends
    if( !((b.fSubBlockCounter+1)%crs::TBlock::fgNSubBlocks) ) {
      s.write(b.fBlockSizeInfo, b.fNBlockSizeInfo);
      // -------- debug -------------
      //std::cout << "COAST> ~~~~~~~ writing " << b.fNBlockSizeInfo << std::endl;
      // -------- debug -------------
    }
	
    ++b.fSubBlockCounter;
	
    return s;
  }
    
};


#endif
