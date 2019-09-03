#include <fstream>
#include <stdexcept>
#include "ResponseTable.h"

using namespace std;


/*
  Constructs a response table from a file.

  Parameters
  ----------
  fname : string
    Path to the data file.
 */
ResponseTable::ResponseTable(string fname)
  : fPhiMin(0.0),
    fPhiMax(0.0),
    fThetaMin(0.0),
    fThetaMax(0.0) //,
//  fEffectiveAreaHist()
{
  ifstream input(fname.c_str());

  // Read angular limits and number of steps
  unsigned int phiSteps, thetaSteps;

  input >> fThetaMin >> fThetaMax >> thetaSteps;
  input >> fPhiMin >> fPhiMax >> phiSteps;

  // Convert angles to radians
//fThetaMin *= TMath::Pi() / 180.0;
//fThetaMax *= TMath::Pi() / 180.0;
//fPhiMin *= TMath::Pi() / 180.0;
//fPhiMax *= TMath::Pi() / 180.0;

  // Generate the edges of equidistant bins covering the range
  // [min - delta / 2, max + delta / 2], so that each bin will be
  // centered at the data point
  const double dPhi = (fPhiMax - fPhiMin) / (phiSteps - 1);
  const vector<double> phiBins = GenerateBinEdges((fPhiMin - dPhi / 2.0),
                                                  (fPhiMax + dPhi / 2.0),
                                                  phiSteps);

  const double dTheta = (fThetaMax - fThetaMin) / (thetaSteps - 1);
  const vector<double> thetaBins = GenerateBinEdges((fThetaMin - dTheta / 2.0),
                                                    (fThetaMax + dTheta / 2.0),
                                                    thetaSteps);

  const string histTitle = "hEffectiveArea_" + fname;
//fEffectiveAreaHist = TH2D(histTitle.c_str(), "",
//                          thetaSteps, &thetaBins[0],
//                          phiSteps, &phiBins[0]);

  // Read the table of effective areas and fill the histogram
  for (unsigned int row = 1; row < thetaSteps + 1; row++) {
    for (unsigned int col = 1; col < phiSteps + 1; col++) {
      double effectiveArea;
      input >> effectiveArea;
//    fEffectiveAreaHist.SetBinContent(row, col, effectiveArea);
    }
  }
}


/*
  Returns the effective area for the given direction.

  Parameters
  ----------
  theta : double
    Zenith angle (rad).

  phi : double
    Azimuth angle counting ccw from east (rad).
 */
double ResponseTable::GetEffectiveAntennaArea(double theta, double phi)
{
  if (!IsWithinLimits(theta, phi)) {
    return 0.0;
  }

  return 1.0; // dummy TH
//return fEffectiveAreaHist.Interpolate(theta, NormaliseAngle(phi));

}


/*
  Checks whether the given direction is contained in the histogram.

  Parameters
  ----------
  theta : double
    Zenith angle (rad).

  phi : double
    Azimuth angle counting ccw from east (rad).
 */
bool ResponseTable::IsWithinLimits(double theta, double phi) const
{
  phi = NormaliseAngle(phi);

  return ((theta >= fThetaMin) && (theta <= fThetaMax) &&
          (phi >= fPhiMin) && (phi <= fPhiMax));
}


/*
  Shifts the given angle into the (-pi, pi] range.

  Parameters
  ----------
  phi : double
    Angle to shift (rad).
 */
double ResponseTable::NormaliseAngle(double phi) const
{
//  while (phi <= -TMath::Pi())
//    phi += 2 * TMath::Pi();

//  while (phi > TMath::Pi())
//    phi -= 2 * TMath::Pi();

  return phi;
}


/*
  Returns the edges of nBins equidistant bins between xLow and xHigh.

  Parameters
  ----------
  xLow, xHigh : double
    Lower (upper) edge of the first (last) bin.

  nBins : size_t
    Number of bins to generate in between.
*/
std::vector<double> ResponseTable::GenerateBinEdges(double xLow,
                                                    double xHigh,
                                                    size_t nBins) const
{
  if (nBins == 0) {
    throw runtime_error("ResponseTable::GenerateBinEdges(): "
                        "Need at least one bin.");
  }

  if (xHigh <= xLow) {
    throw runtime_error("ResponseTable::GenerateBinEdges(): "
                        "xHigh must be larger than xLow.");
  }

  const double dx = (xHigh - xLow) / nBins;
  std::vector<double> edges(nBins + 1);
  edges[0] = xLow;
  edges[nBins] = xHigh;
  for (size_t i = 1; i < nBins; i++) {
    edges[i] = xLow + i * dx;
  }

  return edges;
}
