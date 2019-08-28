#ifndef _INCLUDE_CORSIKA_MPARTICLE_BLOCK
#define _INCLUDE_CORSIKA_MPARTICLE_BLOCK

#include <crs/CorsikaTypes.h>
#include <crs/TSubBlock.h>
#include <crs/TParticleBlockEntry.h>

#include <vector>

namespace crs {

  /** 
      \class MParticleBlock
      \brief CORSIKA particle sub-block, constains 39 particles

      While converting a TSubBlock into a MParticleBlock the underlying data 
      of the sub-block is scanned for valid particles, and a list of these 
      particles is internally stored. This list of particles is accessible 
      using iterators. 

      \author Ralf Ulrich
      \date Thu Feb  3 13:04:50 CET 2005
      \version $Id: MParticleBlock.h,v 1.1.1.1 2007-07-31 07:00:48 rulrich Exp $
  */

  class MParticleBlock : public TSubBlock {
	
  public:
    typedef std::vector <TParticleBlockEntry> ParticleList;
    typedef ParticleList::iterator ParticleListIterator;
    typedef ParticleList::const_iterator ParticleListConstIterator;
	
    MParticleBlock () {}
    MParticleBlock (const TSubBlock &right);
    virtual ~MParticleBlock () {}
	
    ParticleListConstIterator FirstParticle () const 
    {return fParticles.begin ();}
    ParticleListConstIterator LastParticle () const 
    {return fParticles.end ();}
	
	
  private:	
    ParticleList fParticles;

  };

};

#endif
