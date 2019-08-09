/* $Id: IParticleReadout.cc,v 1.2 2007-08-09 15:46:59 rulrich Exp $   */

#include <crs/CorsikaConsts.h>
#include <crs/IParticleReadout.h>

#include <cmath>

using namespace crs;
using namespace std;



double IParticleReadout::GetMass () const {
  const int id = int(abs(ValueAt(0))/1000);
  if (id>0 && id<=gNParticles) 
    return gParticleMass [id-1];
  return -1;	    
}



int IParticleReadout::GetPDGCode () const {
  const int id = int(abs(ValueAt(0)/1000));
  if (id>0 && id<=gNParticles)
    return gParticleCodes [id-1];	
  return -1;
}



double IParticleReadout::GetKinEnergy () const {
  double Px = GetPx();
  double Py = GetPy();
  double Pz = GetPz();
  double P2 = Px*Px + Py*Py + Pz*Pz;
  double mass = GetMass();
  return std::sqrt(P2+mass*mass) - mass;
}



double IParticleReadout::GetTheta () const {
  double Px = GetPx();
  double Py = GetPy();
  double Pz = GetPz();
  double P2 = Px*Px + Py*Py + Pz*Pz;
  return acos(Pz/sqrt(P2));
}
