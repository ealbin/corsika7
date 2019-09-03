#ifndef __RESPONSETABLE_H__
#define __RESPONSETABLE_H__

#include <string>
#include <vector>

//#include <TH2D.h>
//#include <TMath.h>


//#if ROOT_VERSION_CODE < ROOT_VERSION(5,26,0)
//#error "Need ROOT >5.25 for TH2::Interpolate()."
//#endif // if ROOT_VERSION_CODE < ROOT_VERSION(5,26,0)


class ResponseTable
{
 public:
  // Constructor from file
  ResponseTable(std::string fname);
  ResponseTable() {};
  // Returns effective antenna area for given direction
  double GetEffectiveAntennaArea(double theta, double phi);

  // Returns copy of internal histogram
//TH2D GetHistogram() const {
//  return fEffectiveAreaHist;
//};


 private:
  // Checks whether direction is contained in histogram
  bool IsWithinLimits(double theta, double phi) const;

  // Shifts angle into (pi, pi] range
  double NormaliseAngle(double phi) const;

  // Generates bin edges for equidistant bins
  std::vector<double> GenerateBinEdges(double xLow,
                                       double xHigh,
                                       size_t nBins) const;

  // Angular limits
  double fPhiMin, fPhiMax;
  double fThetaMin, fThetaMax;

  // Histogram for effective antenna area
//TH2D fEffectiveAreaHist;
};

#endif
