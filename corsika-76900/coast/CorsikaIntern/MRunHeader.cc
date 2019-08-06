/* $Id: MRunHeader.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs/CorsikaConsts.h>
#include <crs/MRunHeader.h>

#include <cmath>

using namespace crs;

MRunHeader::MRunHeader (const TSubBlock &right)
{
  // ctor
  fSubBlockData = right.fSubBlockData;
  fType = right.fType;
  fThinned = right.fThinned;
}

double
MRunHeader::GetSamplingPlaneNormalX () const
{
  return std::sin(GetSamplingPlaneTheta()*deg) * std::cos(GetSamplingPlanePhi()*deg);
}

double
MRunHeader::GetSamplingPlaneNormalY () const
{
  return std::sin(GetSamplingPlaneTheta()*deg) * std::sin(GetSamplingPlanePhi()*deg);
}

double
MRunHeader::GetSamplingPlaneNormalZ () const
{
  return std::cos(GetSamplingPlaneTheta()*deg);
}

CREAL
MRunHeader::GetVerticalDepth(const CREAL height)
  const
{

  const int n = 5;
  int atmosphereSection = 0;
  if (height < GetAtmosphereLayerBoundary(0))
    return GetAtmosphereA(0);
  else if (height >= GetAtmosphereLayerBoundary(n-1))
    atmosphereSection = n-1;
  else {
    for (int l=0; l<n-1; ++l) {
      if (height >= GetAtmosphereLayerBoundary(l) &&
          height < GetAtmosphereLayerBoundary(l+1)) {
        atmosphereSection = l;
        break;
      }
    }
  }

  const double A = GetAtmosphereA(atmosphereSection);
  const double B = GetAtmosphereB(atmosphereSection);
  const double C = GetAtmosphereC(atmosphereSection);

  return atmosphereSection < 4 ?
    A + B * exp(-height/C) :
    A - B * height/C;
}

