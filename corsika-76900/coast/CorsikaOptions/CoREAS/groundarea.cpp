/***************************************************************************
                          groundarea.cpp  -  description
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

#include "groundarea.h"
#include <sstream>      // for conversion of numerics to strings
#include <sys/stat.h>   // for mkdir
#include <sys/types.h>  // for mkdir
#include <errno.h>      // for mkdir


GroundArea::GroundArea( double parTimeResolution, double parLowerTimeBoundary, double parUpperTimeBoundary,
                        double parAutomaticTimeBoundaries, double parResolutionReductionScale, const ThreeVector& parPositionOfFirstInteraction, double parTimeOfFirstInteraction,
                        vector<AntennaPosition> parAntennaPositions, ThreeVector parRadioCore, ThreeVector parAxis)
:    itsLowestHeight(1.e12), itsRadioCore(parRadioCore), itsAxis(parAxis)
{
  // create and activate the desired antenna positions
  for (vector<AntennaPosition>::const_iterator i=parAntennaPositions.begin(); i!=parAntennaPositions.end(); ++i)
  {
    CreateBin(i->GetPosition()+itsRadioCore, parTimeResolution, parLowerTimeBoundary, parUpperTimeBoundary, parResolutionReductionScale,
              parAutomaticTimeBoundaries, parPositionOfFirstInteraction, parTimeOfFirstInteraction, i->GetName(),
              i->GetModeString(), i->GetModeData1(), i->GetModeData2());
    // scan lowest observer height of all observer positions
    if ((i->GetPosition()+itsRadioCore).GetZ() < itsLowestHeight) { itsLowestHeight = (i->GetPosition()+itsRadioCore).GetZ(); }     
  }
}

GroundArea::~GroundArea()
{
  for (vector<GroundElement*>::iterator i=itsActiveGroundElementPointers.begin(); i!=itsActiveGroundElementPointers.end(); ++i)
  {
    delete (*i);
  }
}


void GroundArea::CreateBin(ThreeVector parPosition, double parTimeResolution, double parTimeLowerBoundary,
               double parTimeUpperBoundary, double parResolutionReductionScale, double parAutomaticTimeBoundaries,
               const ThreeVector& parPositionOfFirstInteraction, double parTimeOfFirstInteraction, string parObserverName,
               string parModeString, string parModeData1, string parModeData2)
{
  // calculate the axis Distance of the observer
  double axisDistance=(parPosition-itsRadioCore).CrossedWith(itsAxis).GetLength();

  // first set up ground element with user-defined time boundaries and resolutions
  GroundElement* ptmpGroundElement = new GroundElement(parPosition,
                                     parTimeResolution,
                                     parTimeLowerBoundary,
                                     parTimeUpperBoundary,
                                     parObserverName,
                                     parModeString,
                                     parModeData1,
                                     parModeData2);
            
  // then decrease time resolution if applicable, overwriting time boundaries if necessary
  if (parResolutionReductionScale >0)
  {
    int numReductions = static_cast<int>(axisDistance/parResolutionReductionScale);
    if (fabs(axisDistance-(numReductions+1)*parResolutionReductionScale) < 1.e-4*parResolutionReductionScale)
      ++numReductions; // compensate for rounding errors

    ptmpGroundElement->DecreaseTimeResolution(numReductions+1);
    double newResolution = parTimeResolution*(numReductions+1);            
    double newLowerBoundary = floor(parTimeLowerBoundary/newResolution)*newResolution;
    double newUpperBoundary = ceil(parTimeUpperBoundary/newResolution)*newResolution;
    ptmpGroundElement->SetTimeBoundaries(newLowerBoundary, newUpperBoundary);
  }
            
  // now, if applicable, overwrite time boundaries with optimised values
  if (parAutomaticTimeBoundaries > 0)
  {
    ThreeVector binPosition = ptmpGroundElement->GetPosition();
    double theRes=ptmpGroundElement->GetTimeResolution();
    double interactionTime = parTimeOfFirstInteraction;
    double lightTravelTime = (parPositionOfFirstInteraction-binPosition).GetLength()/SpeedOfLight;
    double startTime = floor((interactionTime+lightTravelTime)/theRes)*theRes-50.0*theRes;    // round to time resolution and put in 50 bins as safety
    double endTime = ceil((interactionTime+lightTravelTime+(parAutomaticTimeBoundaries*theRes/parTimeResolution))/theRes)*theRes+50.0*theRes; // round to time resolution and put in 50 bins as safety
    ptmpGroundElement->SetTimeBoundaries(startTime, endTime);
  }
            
  ptmpGroundElement->Activate(itsActiveGroundElementPointers);
}

// only create directory and write out .bins file
bool GroundArea::WriteGenericFiles(string parFullPath, string parRelativePath)
{
  errno = 0;
  int result = mkdir(parFullPath.c_str(), 0777);
  if ( (result == 0) || (errno == EEXIST) )
  {
    ofstream bout( (parFullPath+".bins").c_str() ); // open file to list the active bins
    bool allok = static_cast<bool>(bout);                                // opening ok?
    for (vector<GroundElement*>::const_iterator i = itsActiveGroundElementPointers.begin(); i != itsActiveGroundElementPointers.end(); ++i)
    {
      stringstream ss;
      //use relative path
      ss << parRelativePath << "/raw_" << (*i)->GetObserverName() << ".dat";
      // for write-out calculate again antenna positions relative to the radio core, except for z-component and also axis distance
      bout << setprecision(12) << ss.str() << "\t" << (*i)->GetPosition()-ThreeVector(itsRadioCore.GetX(),itsRadioCore.GetY(),0.0) << "\t0.0\t" << ((*i)->GetPosition()-itsRadioCore).CrossedWith(itsAxis).GetLength() << "\n";        // list the active bin
    }
    bout.close();
    return allok;
  }
  else { return false; }        // error
}

// write the results out to files
bool GroundArea::WriteResults(string parFullPath, string parRelativePath)
{
  errno = 0;
  int result = mkdir(parFullPath.c_str(), 0777);
  if ( (result == 0) || (errno == EEXIST) )
  {
    ofstream bout( (parFullPath+".bins").c_str() ); // open file to list the active bins
    bool allok = static_cast<bool>(bout);                                // opening ok?
    for (vector<GroundElement*>::const_iterator i = itsActiveGroundElementPointers.begin(); i != itsActiveGroundElementPointers.end(); ++i)
    {
      stringstream ss;
      //use relative path
      ss << parRelativePath << "/raw_" << (*i)->GetObserverName() << ".dat";
      // for write-out calculate again antenna positions relative to the radio core, except for z-component and also axis distance
      bout << setprecision(12) << ss.str() << "\t" << (*i)->GetPosition()-ThreeVector(itsRadioCore.GetX(),itsRadioCore.GetY(),0.0) << "\t0.0\t" << ((*i)->GetPosition()-itsRadioCore).CrossedWith(itsAxis).GetLength() << "\n";        // list the active bin
      stringstream rawss;
      rawss << parFullPath << "/raw_" << (*i)->GetObserverName() << ".dat";
            
      ofstream fout(rawss.str().c_str());
      if (fout)
      {
        fout << *(*i);                                                                                                // write out the timeseries data
        fout.close();
      }
      else { allok = false; }
    }
    bout.close();
    return allok;
  }
  else { return false; }        // error
}
