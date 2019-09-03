/***************************************************************************
                          scenarioparams.cpp  -  description
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

#include "scenarioparams.h"
#include <iomanip>
#include <sstream>

using std::setprecision;
using std::stringstream;


ScenarioParams::ScenarioParams(const string& parName, const string& parOutName, const string& parGroundFileName)
: itsName(parName),
  itsOutName(parOutName),
  itsGroundFileName(parGroundFileName),
  itsPrimaryParticleEnergy(0.0),
  itsDepthOfShowerMaximum(0.0),                           // careful: slant depth!
  itsDistanceOfShowerMaximum(0.0),
  itsShowerZenithAngle(0.0),
  itsShowerAzimuthAngle(0.0),
  itsMagneticFieldStrength(0.3),
  itsMagneticFieldInclinationAngle(0.0),
  itsCoreCoordinateNorth(0.0),
  itsCoreCoordinateWest(0.0),
  itsCoreCoordinateVertical(0.0),
  itsTimeResolution(1.0e-9),
  itsTimeLowerBoundary(-5.e-5),
  itsTimeUpperBoundary(5.e-5),
  itsAutomaticTimeBoundaries(4.e-7),
  itsResolutionReductionScale(5000),
  itsCorsikaParameterFileName(""),
  itsGroundLevelRefractiveIndex(1.000292),
  itsGenerateHistograms(false),
  itsEventNumber(-1),
  itsRunNumber(-1),
  itsGPSSecs(0),
  itsGPSNanoSecs(0.0),
  itsPrimaryParticleType(0),
  itsCoreEastingOffline(0.0),
  itsCoreNorthingOffline(0.0),
  itsCoreVerticalOffline(0.0),
  itsOfflineCoordinateSystem("Reference"),
  itsRotationAngleForMagfieldDeclination(0.0),
  itsComment("")
{
}

ScenarioParams::~ScenarioParams()
{
}

bool ScenarioParams::ReadFromFiles()
{
  // define flag for successful read operations
  bool Ok;

  // read in and parse reas parameter file
  ifstream fin((itsName+".reas").c_str());
  Ok = fin.is_open();
  fin >> *this;
  fin.close();
  if (Ok)
  {  cout << " CoREAS: Parameter file " << itsName << ".reas successfully imported!\n";  }
  else
  {  cout << " CoREAS: Error reading parameter file " << itsName << ".reas!\n\n" << endl;  return false; }

  // read in and interpret antenna list file
  Ok = ImportAntennaFile();
  if (Ok)
  {  cout << " CoREAS: Antenna list file " << itsGroundFileName << ".list successfully imported!\n";  }
  else
  { cout << " CoREAS: Error reading antenna list file " << itsGroundFileName << ".list!\n\n" << endl; return false; }
  
  return true;  // all files read successfully
}

bool ScenarioParams::WriteToFiles() const
{
  // export the scenario parameters to the out file
  bool FirstOk;
  ofstream fout((itsOutName+".reas").c_str());
  fout << *this;
  FirstOk = fout.is_open();
  fout.close();

  if (not FirstOk) { cout << "Error writing parameter file " << itsName << ".reas!" << endl; }

  return FirstOk;
}




// friend function
ostream& operator<< (ostream& os, const ScenarioParams& rhs)
{
  os << "# CoREAS V" << setprecision(2) << ProgramVersion << " by Tim Huege <tim.huege@kit.edu> with contributions by Marianne Ludwig and Clancy James - parameter file" << setprecision(10) << "\n";
  os << "\n";
  os << "# parameters setting up the spatial observer configuration:\n";  
  os << "\n";
  os << "CoreCoordinateNorth = " << rhs.itsCoreCoordinateNorth << "\t\t\t; in cm\n";
  os << "CoreCoordinateWest = " << rhs.itsCoreCoordinateWest << "\t\t\t; in cm\n";
  os << "CoreCoordinateVertical = " << rhs.itsCoreCoordinateVertical << "\t\t\t; in cm\n";
  os << "\n";
  os << "# parameters setting up the temporal observer configuration:\n";  
  os << "\n";
  os << "TimeResolution = " << rhs.itsTimeResolution << "\t\t\t\t; in s" << "\n";
  os << "AutomaticTimeBoundaries = " << rhs.itsAutomaticTimeBoundaries << "\t\t\t; 0: off, x: automatic boundaries with width x in s" << "\n";
  os << "TimeLowerBoundary = " << rhs.itsTimeLowerBoundary << "\t\t\t\t; in s, only if AutomaticTimeBoundaries set to 0" << "\n";
  os << "TimeUpperBoundary = " << rhs.itsTimeUpperBoundary << "\t\t\t\t; in s, only if AutomaticTimeBoundaries set to 0" << "\n";
  os << "ResolutionReductionScale = " << rhs.itsResolutionReductionScale << "\t\t\t; 0: off, x: decrease time resolution linearly every x cm in radius" << "\n";  
  os << "\n";
  os << "# parameters setting up the simulation functionality:\n";
  os << "GroundLevelRefractiveIndex = " << rhs.itsGroundLevelRefractiveIndex << "\t\t; specify refractive index at 0 m asl\n";
  os << "\n";
  os << "# event information for Offline simulations:\n";  
  os << "\n";
  os << "EventNumber = " << rhs.itsEventNumber << "\n";
  os << "RunNumber = " << rhs.itsRunNumber << "\n";
  os << "GPSSecs = " << rhs.itsGPSSecs << "\n";
  os << "GPSNanoSecs = " << rhs.itsGPSNanoSecs << "\n";
  os << "CoreEastingOffline = " << rhs.itsCoreEastingOffline << "\t\t\t\t; in meters\n";
  os << "CoreNorthingOffline = " << rhs.itsCoreNorthingOffline << "\t\t\t\t; in meters\n";
  os << "CoreVerticalOffline = " << rhs.itsCoreVerticalOffline << "\t\t\t\t; in meters\n";
  os << "OfflineCoordinateSystem = " << rhs.itsOfflineCoordinateSystem << "\n";
  os << "RotationAngleForMagfieldDeclination = " << rhs.itsRotationAngleForMagfieldDeclination << "\t\t; in degrees\n";
  os << "Comment =" << rhs.itsComment << "\n";
  os << "\n";
  os << "# event information for your convenience and backwards compatibility with other software, these values are not used as input parameters for the simulation:\n";  
  os << "\n";
  os << "ShowerZenithAngle = " << rhs.itsShowerZenithAngle << "\t\t\t; in degrees" << "\n";
  os << "ShowerAzimuthAngle = " << rhs.itsShowerAzimuthAngle << "\t\t; in degrees, 0: shower propagates to north, 90: to west" << "\n";
  os << "PrimaryParticleEnergy = " << rhs.itsPrimaryParticleEnergy << "\t\t\t; in eV" << "\n";
  os << "PrimaryParticleType = " << rhs.itsPrimaryParticleType << "\t\t\t; as defined in CORSIKA" << "\n";
  os << "DepthOfShowerMaximum = " << rhs.itsDepthOfShowerMaximum << "\t\t; slant depth in g/cm^2" << "\n";
  os << "DistanceOfShowerMaximum = " << rhs.itsDistanceOfShowerMaximum << "\t\t; geometrical distance of shower maximum from core in cm" << "\n";
  os << "MagneticFieldStrength = " << rhs.itsMagneticFieldStrength << "\t\t\t; in Gauss" << "\n";
  os << "MagneticFieldInclinationAngle = " << rhs.itsMagneticFieldInclinationAngle << "\t\t; in degrees, >0: in northern hemisphere, <0: in southern hemisphere" << "\n";
  os << "GeomagneticAngle = " << acos(cos(rhs.itsShowerAzimuthAngle*M_PI/180.0)*sin(rhs.itsShowerZenithAngle*M_PI/180.0)*cos(rhs.itsMagneticFieldInclinationAngle*M_PI/180.0)+cos(rhs.itsShowerZenithAngle*M_PI/180.0)*sin(rhs.itsMagneticFieldInclinationAngle*M_PI/180.0))*180.0/M_PI << "\t\t\t; in degrees" << "\n";
  os << "CorsikaFilePath = ./\n";
  os << "CorsikaParameterFile = " << rhs.itsCorsikaParameterFileName << "\n";
  return os;
}


istream& operator>> (istream& is, ScenarioParams& rhs)
{
   char lineBuf[1024];
   while ( is.getline(lineBuf,1024) )
   {
    string currLine(lineBuf);
    if (currLine.substr(0,1) == string("#")) { continue; }  // commentary line
    unsigned int pos = currLine.find("=",1);                // search for "=" starting with second character
    if ( pos < currLine.size() )
    {
      string type = currLine.substr(0,pos);              // everything before the "="

      string value;
      unsigned int pose = currLine.find(";",pos);        // search for ";" starting with the position of the "=" found before
      if ( pose < currLine.size() )
      {
        // comment found, only use part before the ";"
        value =  currLine.substr(pos+1,pose-pos-1);    // only part before the ";"
      }
      else
      {
        value = currLine.substr(pos+1);                // everything after the "="
      }

      unsigned int tend = type.size();
      stringstream ss;
      double dd; long ll; ThreeVector tv; string sr; // temporary variables

      ss << value;  // feed the value into the stringstream for conversion

      // now check for the individual types
      if (type.find("CoreCoordinateNorth") < tend) { ss >> dd; rhs.itsCoreCoordinateNorth=dd; continue; }
      if (type.find("CoreCoordinateWest") < tend) { ss >> dd; rhs.itsCoreCoordinateWest=dd; continue; }
      if (type.find("CoreCoordinateVertical") < tend) { ss >> dd; rhs.itsCoreCoordinateVertical=dd; continue; }
      if (type.find("TimeResolution") < tend) { ss >> dd; rhs.itsTimeResolution=dd; continue; }
      if (type.find("TimeLowerBoundary") < tend) { ss >> dd; rhs.itsTimeLowerBoundary=dd; continue; }
      if (type.find("TimeUpperBoundary") < tend) { ss >> dd; rhs.itsTimeUpperBoundary=dd; continue; }
      if (type.find("GroundLevelRefractiveIndex") < tend) { ss >> dd; rhs.itsGroundLevelRefractiveIndex=dd; continue; }
      if (type.find("AutomaticTimeBoundaries") < tend) { ss >> dd; rhs.itsAutomaticTimeBoundaries=dd; continue; }
      if (type.find("ResolutionReductionScale") < tend) { ss >> dd; rhs.itsResolutionReductionScale=dd; continue; }
     
      if (type.find("EventNumber") < tend) { ss >> ll; rhs.itsEventNumber=ll; continue; }
      if (type.find("RunNumber") < tend) { ss >> ll; rhs.itsRunNumber=ll; continue; }
      if (type.find("GPSSecs") < tend) { ss >> ll; rhs.itsGPSSecs=ll; continue; }
      if (type.find("GPSNanoSecs") < tend) { ss >> dd; rhs.itsGPSNanoSecs=dd; continue; }
      if ((type.find("CoreEastingOffline") < tend) || (type.find("CoreEastingPampaAmarilla") < tend)) { ss >> dd; rhs.itsCoreEastingOffline=dd; continue;}
      if ((type.find("CoreNorthingOffline") < tend) || (type.find("CoreNorthingPampaAmarilla") < tend)) { ss >> dd; rhs.itsCoreNorthingOffline=dd; continue;}
      if ((type.find("CoreVerticalOffline") < tend) || (type.find("CoreVerticalPampaAmarilla") < tend)) { ss >> dd; rhs.itsCoreVerticalOffline=dd; continue;}
      if (type.find("OfflineCoordinateSystem") < tend) { rhs.itsOfflineCoordinateSystem=ss.str(); continue; }
      if (type.find("RotationAngleForMagfieldDeclination") < tend) { ss >> dd; rhs.itsRotationAngleForMagfieldDeclination=dd; continue;}
      if (type.find("Comment") < tend) { rhs.itsComment=ss.str(); continue; }
      if (type.find("CorsikaParameterFile") < tend) { ss >> sr; rhs.itsCorsikaParameterFileName=sr; continue; }
    }
  }

  rhs.ConversionAndSanityCheck();
  return is;        // actually not of much use, since file was read to the end anyway
}


void ScenarioParams::ConversionAndSanityCheck()
{
  // construct axes from user-specified angles

  if (itsGroundLevelRefractiveIndex < 1.0) { itsGroundLevelRefractiveIndex = 1.0; } // cannot be smaller than unity
  if (itsResolutionReductionScale < 0.0) {itsResolutionReductionScale = 0.0; }
  
  // overrides (e.g., because a function is not correctly implemented yet)
}


bool ScenarioParams::ImportAntennaFile()
{
  bool successfulAntennaFileImport = false;
  // actually antenna list setup should be capsulated in own class, imported through operator<< and exported through operator>>
  ifstream gin( (itsGroundFileName+".list").c_str() );
  if (gin)
  {
     char lineBuf[1024];
     while ( gin.getline(lineBuf,1024) )
     {
      string currLine(lineBuf);
      if (currLine.substr(0,1) == string("#")) { continue; }  // commentary line
      unsigned int pos = currLine.find("=",1);              // search for "=" starting with second character
      if ( pos < currLine.size() )
      {
        string type = currLine.substr(0,pos);              // everything before the "="

        string value;
        unsigned int pose = currLine.find(";",pos);        // search for ";" starting with the position of the "=" found before
        if ( pose < currLine.size() )
        {
          // comment found, only use part before the ";"
          value =  currLine.substr(pos+1,pose-pos-1);    // only part before the ";"
        }
        else
        {
          value = currLine.substr(pos+1);                // everything after the "="
        }

        unsigned int tend = type.size();
        stringstream ss;

        ss << value;  // feed the value into the stringstream for conversion

        if (type.find("AntennaPosition") < tend)
        {
          double dx, dy, dz;
          string id, s1, s2, s3;
          ss >> dx >> dy >> dz >> id >> s1 >> s2 >> s3;
//        if (s1 == "pattern") {
//          cout << " CoREAS: Antenna response pattern file name for observer " << id << " is " << s2 << endl;
//        }
          // calculate the antenna coordinates relative to the radio core (caution, handling of z component was different in REAS3)
          itsAntennaPositions.push_back(AntennaPosition(ThreeVector(dx-itsCoreCoordinateNorth,dy-itsCoreCoordinateWest,dz-itsCoreCoordinateVertical), id, s1, s2, s3));
          continue;
        }
      }
    }
    successfulAntennaFileImport = true;
  }
  else
  {
    successfulAntennaFileImport = false;
  }
  gin.close();

  return successfulAntennaFileImport;
}

