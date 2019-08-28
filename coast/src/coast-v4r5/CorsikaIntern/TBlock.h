#ifndef _INCLUDE_CORSIKA_TBLOCK_
#define _INCLUDE_CORSIKA_TBLOCK_

#include <crs/TSubBlock.h>
#include <crs/CorsikaTypes.h>

#include <vector>

/**
   \namespace crs      
   \brief Representation of internal CORSIKA data-structures

   This namespace contains all C++ wrapper objects for the internal CORSIKA
   data (block/ subblock). This is realized by a set of classes, 
   knowing how to access the binary CORSIKA data.

   \author Ralf Ulrich
   \date Thu Feb  3 13:01:59 CET 2005
*/

namespace crs {

  /** 
      \class TBlock
      \brief One CORSIKA data-block, containing 21 sub-blocks.

      This is the main CORSIKA data object. It mirrors the binary 
      block of fortran REAL*4 data
      with a variable size for thinned and not-thinned CORSIKA runs.
      During creation of a TBlock the underlying data is scanned for 
      valid TSubBlock. It is possible to access these TSubBlock using 
      iterators.

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: TBlock.h,v 1.1.1.1 2007-07-31 07:00:49 rulrich Exp $
  */

  class TBlock {
	
  public:
    typedef std::vector<TSubBlock> SubBlocks; 
    typedef SubBlocks::iterator SubBlockIterator;
    typedef SubBlocks::const_iterator SubBlockConstIterator;
	
    TBlock (const CREAL *data, bool thinned);
    ~TBlock () {}
	
    SubBlockConstIterator GetFirstSubBlock () const 
    {return fSubBlocks.begin ();}
    SubBlockConstIterator GetLastSubBlock () const 
    {return fSubBlocks.end ();}
	
  private:
    SubBlocks fSubBlocks;
    bool fThinned;
	
  public:
    static const int fgNSubBlocks;
    static const int fgNParticles;
    static const int fgNLongitudinals;
    static const int fgNEntriesThinned;
    static const int fgNEntriesNotThinned;
  };

};

#endif
