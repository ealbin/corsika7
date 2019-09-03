/***************************************************************************
                          groundarea.h  -  description
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

#ifndef GROUNDAREA_H
#define GROUNDAREA_H

#include "definitions.h"
#include "groundelement.h"
#include "antennaposition.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>

using std::ofstream;
using std::vector;
using std::stringstream;
using std::setprecision;
using std::string;

/**A section on the ground.
  *@author Tim Huege
  */

class GroundArea {
public: 
  GroundArea( double parTimeResolution, double parLowerTimeBoundary, double parUpperTimeBoundary,
              double parAutomaticTimeBoundaries, double parResolutionReductionScale, const ThreeVector& parPositionOfFirstInteraction, double parTimeOfFirstInteraction,
              vector<AntennaPosition> parAntennaPositions, ThreeVector parRadioCore, ThreeVector parRadioAxis);

  ~GroundArea();

  bool WriteGenericFiles(string parFullPath, string parRelativePath);
  bool WriteResults(string parFullPath, string parRelativePath);

  void CreateBin(  ThreeVector parPosition, double parTimeResolution, double parTimeLowerBoundary,
    double parTimeUpperBoundary, double parResolutionReductionScale, double parAutomaticTimeBoundaries,
    const ThreeVector& parPositionOfFirstInteraction, double parTimeOfFirstInteraction, string parObserverName, 
    string parModeString, string parModeData1, string parModeData2 );

  // retrieving the GroundElements to be calculated
  const vector<GroundElement*>& GetActiveGroundElementPointers() { return itsActiveGroundElementPointers; }
    
  // retrieving other data
  double GetLowestHeight() const { return itsLowestHeight; }
  long GetNumActiveGroundElements() const { return itsActiveGroundElementPointers.size(); }
    
  class Exception { };  // exception class
  
private:

  // data
  double itsLowestHeight;
  ThreeVector itsRadioCore;
  ThreeVector itsAxis;

  vector<GroundElement*> itsActiveGroundElementPointers;
};


#endif
