/* $Id: IParticleReadout.cc 5116 2016-01-04 19:09:04Z darko $   */

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
  double P2 = GetPSquared();
  double mass = GetMass();
  return std::sqrt(P2+mass*mass) - mass;
}


double
IParticleReadout::GetPSquared()
  const
{
  return pow(GetPx(), 2) + pow(GetPy(), 2) +  pow(GetPz(), 2);
}



double IParticleReadout::GetTheta () const {
  double Pz = GetPz();
  double P2 = GetPSquared();
  return acos(Pz/sqrt(P2));
}

double
IParticleReadout::GetFeynmanX(const IParticleReadout& beam)
  const
{

  // (m_p + m_n)/2
  static const double targetMass =
    0.5 * (gParticleMass[12] + gParticleMass[13]);
  static const double targetMass2 = targetMass * targetMass;
  const double beamMass2 = pow(beam.GetMass(), 2);

  const double pBeam2 = beam.GetPSquared();
  const double beamEnergy = sqrt(pBeam2 + beamMass2);
  const double s2 =  beamMass2 + targetMass2 + 2 * beamEnergy * targetMass;
  const double sqrtS = sqrt(s2);
  const double gammaCMS = (beamEnergy + targetMass)/sqrtS;
  const double betaCMS = sqrt(pBeam2) / (beamEnergy + targetMass);

  const double p2 = GetPSquared();
  const double labEnergy = sqrt(p2 + pow(GetMass(), 2));
  const double pZ =
    (GetPx() * beam.GetPx() +
     GetPy() * beam.GetPy() +
     GetPz() * beam.GetPz()) / sqrt(pBeam2);
  const double pZCMS =
    -gammaCMS * betaCMS * labEnergy + gammaCMS * pZ;

  return 2.*pZCMS/sqrtS;

}

