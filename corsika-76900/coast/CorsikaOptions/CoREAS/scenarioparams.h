/***************************************************************************
                          scenarioparams.h  -  description
                             -------------------
    begin                : Wed Oct 29 2003
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

#ifndef SCENARIOPARAMS_H
#define SCENARIOPARAMS_H

#include <string>
#include <vector>
#include "threevector.h"
#include "definitions.h"
#include <fstream>
#include <iostream>
#include "antennaposition.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

class Scenario;
class Shower;
class Observer;
class GroundBinList;

/**
  *Class handling the input and output as well as storing of parameters for a scenario.
  *Provides parsing and version checking of parameter files, makes sure parameters are consistent and writes back parameters used at end of simulation. Contains all parameters of the the scenario in one concentrated class.
  *@author Tim Huege
  */

class ScenarioParams {
public: 
  ScenarioParams(const string& parName, const string& parOutName, const string& parGroundFileName); ///< constructor creating scenario parName with associated ground file parGroundFileName
       
  // do conversion and sanity check of imported parameters
  void ConversionAndSanityCheck();  ///< converts angles into vectors and does some sanity checks, has to be called at end of parsing input files

  // input and output from/to parameter file
  bool ReadFromFiles();  ///< reads in the parameter file
  bool WriteToFiles() const;  ///< writes the parameter file to disk

  ~ScenarioParams();

  friend ostream& operator<< (ostream& os, const ScenarioParams& rhs);

private:
  // define friends that need to access the data
  friend class TPlotter;
  
  // function to import the groundfile
  bool ImportAntennaFile(); ///< import the antenna list file
  bool ImportCorsikaParameterFile();  ///< read in the corsika input parameter file
  
  // member variables
  string itsName;
  string itsOutName;
  string itsGroundFileName;
  double itsPrimaryParticleEnergy;
  double itsDepthOfShowerMaximum;
  double itsDistanceOfShowerMaximum;
  double itsShowerZenithAngle;
  double itsShowerAzimuthAngle;
  double itsMagneticFieldStrength;
  double itsMagneticFieldInclinationAngle;
  double itsCoreCoordinateNorth;
  double itsCoreCoordinateWest;
  double itsCoreCoordinateVertical;
  double itsTimeResolution;
  double itsTimeLowerBoundary;
  double itsTimeUpperBoundary;

  bool itsResponseTableUsed;

  vector<AntennaPosition> itsAntennaPositions;
  
  double itsAutomaticTimeBoundaries;
  double itsResolutionReductionScale;
  
  mutable string itsCorsikaParameterFileName;
  double itsGroundLevelRefractiveIndex;
  bool itsGenerateHistograms;
    
  long itsEventNumber;
  long itsRunNumber;
  long itsGPSSecs;
  double itsGPSNanoSecs;
  long itsPrimaryParticleType;
  double itsCoreEastingOffline;
  double itsCoreNorthingOffline;
  double itsCoreVerticalOffline;
  string itsOfflineCoordinateSystem;
  double itsRotationAngleForMagfieldDeclination;
  string itsComment;
         
  // readin routine needs to be able to write to member variables
  friend istream& operator>> (istream& is, ScenarioParams& rhs);

};


#endif
