#ifndef _INCLUDE_CORSIKA_IPARTICLEREADOUT_
#define _INCLUDE_CORSIKA_IPARTICLEREADOUT_

#include <crs/CorsikaTypes.h>
#include <cmath>

namespace crs {

  /**
      \class IParticleReadout
      \brief CORSIKA common particle data readout interface

      One CORSIKA particle definition.

      \author Ralf Ulrich
      \date Sun Jun 19 12:06:45 CEST 2005
      \version $Id: IParticleReadout.h 5116 2016-01-04 19:09:04Z darko $
  */

  class IParticleReadout {

  public:
    IParticleReadout () {}
    virtual ~IParticleReadout () {}


    virtual int GetObservationLevel () const
    {return int(std::abs(ValueAt(0)))%10;}
    virtual int GetHadronicGeneration () const
    {return int(std::abs(ValueAt(0)))%1000/10;}


    virtual CREAL GetPx () const {return ValueAt (1);}
    virtual CREAL GetPy () const {return ValueAt (2);}
    virtual CREAL GetPz () const {return ValueAt (3);}

    virtual CREAL GetX () const {return ValueAt (4);}
    virtual CREAL GetY () const {return ValueAt (5);}

    virtual CREAL GetTime () const {return ValueAt(6); }
    virtual CREAL GetWeight () const {return ValueAt (7);}

  public:
    double GetMass () const;         ///< mass in GeV
    int GetPDGCode () const;
    double GetKinEnergy () const;    ///< kin. energy in GeV
    double GetTheta () const;        ///< zenith angle in rad
    double GetPSquared() const;      ///< squared of momentum in GeV
    double GetE() const { return std::sqrt(GetPSquared() + std::pow(GetMass(),2)); }
    /// wrt. beam particle (e.g. mother)
    double GetFeynmanX(const IParticleReadout& beam) const;
  private:
    virtual CREAL ValueAt (int i) const = 0;

  };



};



#endif
