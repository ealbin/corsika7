#ifndef _INCLUDE_CORSIKA_MLONGITUDINAL_BLOCK
#define _INCLUDE_CORSIKA_MLONGITUDINAL_BLOCK

#include <crs/TSubBlock.h>
#include <crs/CorsikaTypes.h>
#include <crs/TLongitudinal.h>

#include <vector>

namespace crs {

  /** 
      \class MLongitudinalBlock
      \brief CORSIKA longitudinal sub-block, constaining 26 steps

      While converting a TSubBlock into a MLongitudinalBlock the underlying data 
      of the sub-block is scanned for valid longitudinal steps, and a list of 
      these 
      steps is internally stored. This list of steps is accessible 
      using iterators. 

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: MLongitudinalBlock.h 5116 2016-01-04 19:09:04Z darko $
  */

  class MLongitudinalBlock : public TSubBlock {
	
  public:
    typedef std::vector <TLongitudinal> LongitudinalList;
    typedef LongitudinalList::iterator LongitudinalListIterator;
    typedef LongitudinalList::const_iterator LongitudinalListConstIterator;
	
    MLongitudinalBlock () {}
    MLongitudinalBlock (const TSubBlock &right);
    ~MLongitudinalBlock () {}
	
    int GetNSteps () const {return (int) fSubBlockData [4]/100;} // total #
    int GetNLongBlocks () const {return (int) fSubBlockData [4]%100;}
    int GetBlockNumber () const {return (int) fSubBlockData [5];}


    LongitudinalListConstIterator FirstLongitudinal () const 
    {return fLongitudinals.begin ();}
    LongitudinalListConstIterator LastLongitudinal () const 
    {return fLongitudinals.end ();}
	
  private:
    void ScanLongBlock ();

    static const int fgEntriesPerLongInfo;
    static const int fgLongInfoOffset;

  private:	
    LongitudinalList fLongitudinals;

  };
};

#endif
