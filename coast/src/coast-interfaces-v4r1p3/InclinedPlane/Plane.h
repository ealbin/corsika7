#ifndef _crs_Plane_h_
#define _crs_Plane_h_

#include <crs/CParticle.h>
#include <crs/CorsikaConsts.h> // for particle masses

#include <cmath>

namespace crs {

  class Plane {
  public:
    Plane() : fD(0)
    { fPoint[0] = fPoint[1] = fPoint[2] = fNormal[0] = fNormal[1] = fNormal[2] = 0; }

    Plane(const double point[3], const double normal[3])
    {
      for (int i = 0; i < 3; ++i) {
        fPoint[i] = point[i];
        fNormal[i] = normal[i];
      }
      fD = Scalar(fPoint, fNormal);
    }

    Plane(const double x,  const double y,  const double z,
          const double nx, const double ny, const double nz)
    {
      fPoint[0] = x;
      fPoint[1] = y;
      fPoint[2] = z;
      fNormal[0] = nx;
      fNormal[1] = ny;
      fNormal[2] = nz;
      fD = Scalar(fPoint, fNormal);
    }

    double GetZ(const double x, const double y) const
    { return (fD - x * fNormal[0] - y * fNormal[1]) / fNormal[2]; }

    const double* GetPoint() const { return fPoint; }

    const double* GetNormal() const { return fNormal; }

    // note that p1 has to be on the normal side, p2 on the opposite
    // for efficiency p2 should be final point and normal oriented towards
    // incoming shower
    bool IsIntersecting(const CParticle& p1, const CParticle& p2) const
    { return Scalar(&p2.x, fNormal) <= fD && fD <= Scalar(&p1.x, fNormal); }

    // check for intersection and interpolate particle if true
    bool
    IsIntersecting(const CParticle& p1, const CParticle& p2, crs::GroundParticle& interp)
      const
    {
      const double p2n = Scalar(&p2.x, fNormal);
      if (p2n <= fD) {
        const double p1n = Scalar(&p1.x, fNormal);
        if (fD <= p1n) {
          interp.Id = p1.particleId;
          interp.generation = p1.hadronicGeneration;
          interp.weight = p1.weight;
          const double d1 = p1n - fD;
          const double d2 = p2n - fD;
          const double eps = d1 / (d1 - d2);
          const Interpolator interpolate(eps);
          // interpolate
          interp.time = interpolate(p1.time, p2.time);
          interp.x = interpolate(p1.x, p2.x);
          interp.y = interpolate(p1.y, p2.y);
          interp.z = interpolate(p1.z, p2.z);

          const double eTot = interpolate(p1.energy, p2.energy);
	  const double mass = crs::gParticleMass[interp.Id];
	  const double pTot = std::sqrt( (eTot-mass)*(eTot+mass) ) ;
	  
	  const double n[3] = {p2.x-p1.x, p2.y-p1.y, p2.z-p1.z};
	  const double length = std::sqrt(std::pow(n[0], 2) +
					  std::pow(n[1], 2) +
					  std::pow(n[2], 2));
	  if (length) {
	    interp.Px = n[0] / length * pTot;
	    interp.Py = n[1] / length * pTot;
	    interp.Pz = n[2] / length * pTot;
	  } else {
	    interp.Px = 0;	    
	    interp.Py = 0;	    
	    interp.Pz = 0;	    
	  }
          return true;
        }
        return false;
      }
      return false;
    }

  private:
    // this should be a template since CDOUBLEPRECISION can change...
    template<typename T, typename U>
    static
    double
    Scalar(const T* x, const U* y)
    {
      double sum = 0;
      for (int i = 0; i < 3; ++i)
        sum += *(x++) * (*(y++));
      return sum;
    }

    class Interpolator {
    public:
      Interpolator(const double eps) : fEps1(1 - eps), fEps2(eps) { }
      double operator()(const double x1, const double x2) const
	{ return fEps1 * x1 + fEps2 * x2; }
    private:
      double fEps1;
      double fEps2;
    };

    double fPoint[3];
    double fNormal[3];
    double fD;
  };

}

#endif
