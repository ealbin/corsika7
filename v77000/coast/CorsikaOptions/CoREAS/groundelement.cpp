/***************************************************************************
                          groundelement.cpp  -  description
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

#include "groundelement.h"
#include <sstream>
#include <stdexcept>
#include <fstream>

GroundElement::GroundElement(ThreeVector parPosition, double parTimeResolution, double parTimeLowerBoundary,
          double parTimeUpperBoundary, string parObserverName, string parModeString, string parModeData1, string parModeData2)
: itsPosition(parPosition), itsUseResponseTable(false), itsState(Inactive), itsMode(Normal),
  itsModeDataLow(0.0), itsModeDataHigh(0.0),
  itsTimeResolution(parTimeResolution), itsTimeResolutionReductionFactor(1.0),
  itsTimeLowerBoundary(parTimeLowerBoundary), itsTimeUpperBoundary(parTimeUpperBoundary),
  itsLastShowerFractionFinished(0.0), itsObserverName(parObserverName)
{
  // see if a special mode is foreseen for this groundbin and set it up appropriately

/*  if (parModeString == "pattern") {
    itsMode = Pattern;
    itsUseResponseTable = true;
    itsResponseTable = ResponseTable(parModeData1);
  } else */
  if (parModeString == "gamma")
  {
    itsMode = Gamma;
    stringstream ss;
    ss << parModeData1 << "\t" << parModeData2;
    ss >> itsModeDataLow >> itsModeDataHigh;
    if (itsModeDataLow >= itsModeDataHigh)
      swapvars(itsModeDataLow, itsModeDataHigh);
  }
  else if (parModeString == "height")
  {
    itsMode = Height;
    stringstream ss;
    ss << parModeData1 << "\t" << parModeData2;
    ss >> itsModeDataLow >> itsModeDataHigh;
    if (itsModeDataLow >= itsModeDataHigh)
      swapvars(itsModeDataLow, itsModeDataHigh);
  }
  else if (parModeString == "slantdepth")
  {
    itsMode = SlantDepth;
    stringstream ss;
    ss << parModeData1 << "\t" << parModeData2;
    ss >> itsModeDataLow >> itsModeDataHigh;
    if (itsModeDataLow >= itsModeDataHigh)
      swapvars(itsModeDataLow, itsModeDataHigh);
  }
  else if (parModeString == "distance")
  {
    itsMode = Distance;
    stringstream ss;
    ss << parModeData1 << "\t" << parModeData2;
    ss >> itsModeDataLow >> itsModeDataHigh;
    if (itsModeDataLow >= itsModeDataHigh)
      swapvars(itsModeDataLow, itsModeDataHigh);
  }
}


GroundElement::~GroundElement()
{
  if ( itsState == Active )  // ground elements are always reactivated at the end of the calculation, there are no passive ground elements
  {
    delete pitsTimeGrid;
  }
  // destruction of the groundelements, no need to clean up further (e.g., keep track of active groundelements)
}



void GroundElement::Activate(vector<GroundElement*>& parActiveGroundElementPointers)
{
  if (itsState == Inactive)
  {
    itsState = Active;
    CreateTimeGrid();                  // initialize the inner time grid

    // should check for out of memory!!!!!!
    parActiveGroundElementPointers.push_back(this);
  }
}


void GroundElement::CreateTimeGrid()
{
  pitsTimeGrid = new Grid<timestamp, EFieldVector>(itsTimeResolution, itsTimeLowerBoundary, itsTimeUpperBoundary);
}


double GroundElement::GetEffectiveAntennaArea(double theta, double phi)
{
  if (!itsUseResponseTable)
    return 1.0;
  else
    return itsResponseTable.GetEffectiveAntennaArea(theta, phi);
}

void GroundElement::CopySample(long num, double& x, double& y, double& z) const
{
  ThreeVector tempData;
  double t;
  pitsTimeGrid->ReadSample(num, t, tempData);
  x = tempData.GetX();
  y = tempData.GetY();
  z = tempData.GetZ(); 
}

void GroundElement::EnlargeSampleLengthTo(unsigned long numSamples)
{
  pitsTimeGrid->EnlargeNumberOfEntriesBy(numSamples-pitsTimeGrid->GetNumEntries());  
}

