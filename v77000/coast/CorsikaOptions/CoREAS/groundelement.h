/***************************************************************************
                          groundelement.h  -  description
                             -------------------
    begin                : Tue May 27 2003
    copyright            : (C) 2003 by Tim Huege
    email                : tim.huege@kit.edu  
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *  Copyright and any other appropriate legal protection of these          *
 *  computer programs and associated documentation are reserved in         *
 *  all countries of the world.                                            *
 *                                                                         *
 *  These programs or documentation may not be reproduced by any method    *
 *  without prior written consent of Karlsruhe Institute of Technology or its    *
 *  delegate. Commercial and military use are explicitly forbidden.        *
 *                                                                         *
 *  Karlsruhe Institute of Technology welcomes comments concerning the REAS      *
 *  code but undertakes no obligation for maintenance of the programs,     *
 *  nor responsibility for their correctness, and accepts no liability     *
 *  whatsoever resulting from the use of its programs.                     *
 *                                                                         *
 ***************************************************************************/

#ifndef GROUNDELEMENT_H
#define GROUNDELEMENT_H

#include <vector>
#include <algorithm>
#include "datapoint.h"
#include "definitions.h"
#include "grid.h"
#include "threevector.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "ResponseTable.h"

using std::string;
using std::ostream;
using std::vector;
using std::sort;
using std::setprecision;
using std::cout;
using std::find;
using std::stringstream;

/**An actual bin on the ground.
  *@author Tim Huege
  */

enum State { Active, Passive, Inactive };

enum Mode { Normal, Pattern, Gamma, Height, SlantDepth, Distance };


class GroundElement {
public: 
  GroundElement(ThreeVector parPosition, double parTimeResolution, double parTimeLowerBoundary,
                double parTimeUpperBoundary, string parObserverName, string parModeString, string parModeData1, string parModeData2);
  ~GroundElement();

  void Activate(vector<GroundElement*>& parActiveGroundElementPointers);

  double GetRadius() const { return sqrt(itsPosition.GetX()*itsPosition.GetX()+itsPosition.GetY()*itsPosition.GetY()); }
  double GetAzimuth() const { return atan3(itsPosition.GetX(),itsPosition.GetY()); }
  double GetX() const { return itsPosition.GetX(); }
  double GetY() const { return itsPosition.GetY(); }
  double GetZ() const { return itsPosition.GetZ(); }
  bool IsActive() const { return (itsState == Active); }

  double GetHeight() const { return itsHeight; }
  void SetHeight(double parHeight) { itsHeight = parHeight; }

  Mode GetMode() const { return itsMode; }
  double GetModeDataLow() const { return itsModeDataLow; }
  double GetModeDataHigh() const { return itsModeDataHigh; }  
  
  ThreeVector GetPosition() const { return itsPosition; }

  double GetEffectiveAntennaArea(double theta, double phi);
  void CollectValuePair(const EFieldDataPoint& start, const EFieldDataPoint& end) {  pitsTimeGrid->IncorporateValuePair(start, end); }

  double GetTimeResolution() const { return itsTimeResolution; }
  double GetTimeResolutionReductionFactor() const { return itsTimeResolutionReductionFactor; }
  double GetTimeLowerBoundary() const { return itsTimeLowerBoundary; }
  double GetTimeUpperBoundary() const { return itsTimeUpperBoundary; }
  string GetObserverName() const { return itsObserverName; }
  
  long GetNumSamples() const { return pitsTimeGrid->GetNumEntries(); }
  void CopySample(long num, double& x, double& y, double& z) const;
  void EnlargeSampleLengthTo(unsigned long numSamples);

  // needed to overwrite time boundaries with optimised values
  void SetTimeBoundaries(double parLowerBoundary, double parUpperBoundary) { itsTimeLowerBoundary=parLowerBoundary; itsTimeUpperBoundary=parUpperBoundary; }
  // needed to adjust time resolution to optimised values
  void DecreaseTimeResolution(double parFactor){ itsTimeResolutionReductionFactor*=parFactor; itsTimeResolution*=parFactor; }

  // allow direct access to member tables, no encapsulation but practical for this purpose
  std::vector<float>& GetEarthIntersectionDistances() { return fEarthIntersectionDistances; }
  std::vector<unsigned short>& GetEarthIntersectionAngleBins() { return fTanEarthIntersectionAngleBins; }
      
  friend ostream& operator<<(ostream& os, const GroundElement& rhs);
  
private:

  void CreateTimeGrid();

  ThreeVector itsPosition;
  double itsHeight;

  bool itsUseResponseTable;
  ResponseTable itsResponseTable;

  State itsState;
  Mode itsMode;
  
  double itsModeDataLow;
  double itsModeDataHigh;
  
  double itsTimeResolution;
  double itsTimeResolutionReductionFactor;
  double itsTimeLowerBoundary;
  double itsTimeUpperBoundary;

  bool itsTimeResolutionReduction;
  double itsLastShowerFractionFinished;
  string itsObserverName;
  
  Grid<timestamp, EFieldVector>* pitsTimeGrid;                      // the current grid (during the particle block)

  std::vector<float> fEarthIntersectionDistances;                 // tables for curved refractivity
  std::vector<unsigned short> fTanEarthIntersectionAngleBins;

};


//is a friend
inline
ostream& operator<<(ostream& os, const GroundElement& rhs)
{
  if (rhs.itsState == Active)
  {
    os << setprecision(14);            // set high precision
    os << (*rhs.pitsTimeGrid);
  }
  return os;
}


#endif
