/****************************************************************************
 *                                                                          *
 *  Copyright and any other appropriate legal protection of these           *
 *  computer programs and associated documentation are reserved in          *
 *  all countries of the world.                                             *
 *                                                                          *
 *  These programs or documentation may not be reproduced by any method     *
 *  without prior written consent of Karlsruhe Institute of Technology (KIT)*
 *  ot its delegate. Commercial and military use are explicitly forbidden.  *
 *                                                                          *
 *  The Karlsruhe Institute of Technology welcomes comments concerning the  *
 *  COAST code but undertakes no obligation for maintenance of the programs,*
 *  nor responsibility for their correctness, and accepts no liability      *
 *  whatsoever resulting from the use of its programs.                      *
 *                                                                          *
 ****************************************************************************/

#include <TPlotter.h>
#include <TCorsika.h>

#include <interface/CorsikaInterface.h>
#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/CInteraction.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>
using namespace crs;

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include "groundelement.h"
#include "groundarea.h"
#include "groundelement.h"
#include "scenarioparams.h"
using namespace std;

#if __PARALLELIB__
#include <mpi.h>
#include <stdio.h>
#include <cstring>
#endif



/* -------------------------------------------------------------------
   Define the particle names/types and CORSIKA ID here !!!!!!!!!!!!!!!!
*/
const unsigned int nParticles = 2;
const char* gParticle_name [nParticles] = {"electron", "positron"};
const int gParticle_id     [nParticles] = {3, 2};
// -------------------------------------------------------------------


TPlotter::TPlotter()
: fTabulationThresholdZenithAngle(68.0*M_PI/180.),
  fNumTanThetaBins(10000),
  fTanThetaBinSpacing(0.006),
  fCORSIKA(0),
  fGroundArea(0),
  itsScenarioParams(0),
  itsShowerAxis(0,0,0)
{
  fVerbosityLevel = 0; // no debug output (set to 40 for all from TPlotter, 100 for all including TCORSIKA)

  // set up such that first particle is considered primary
  fPrimaryTrack = true;

  // set up safe initialization values
  fThinning = false;
  fSlant = false;
  fCurved = false;
  fStackInput = false;
  fPreshower = false;

  fCorsikaVersion = "";
  fCoastVersion = "";
  
  fBx = 0;
  fBz = 0;

  fDirectory = "";
  fFileName = "";
  fThreadName = "";

  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  fFirstInteractionTime = 0;
  fFirstInteractionDist = 0;
  fFirstInteractionSlantDepth = 0;
  
  //~ fAxisX = 0;
  //~ fAxisY = 0;
  //~ fAxisZ = 0;
  
  fCosZenith = 0;
  fSinZenith = 0;
  fCosAzimuth = 0;
  fSinAzimuth = 0;
  
  fEventNo = 0;
  fRunNo = 0;
  fObservationLevel = 0;
  fMaxShowerRange = 0;
  fHeightFirstInt = 0;
  fZenith = 0;
  fAzimuth = 0;
  fShowerEnergy = 0;
  fShowerPrimary = 0;
  fXmax = 0;

  fSkimming = false;
  fSkimmingAltitude = 0;

  fCoreHitTime = 0.0;
  fRadioCore = ThreeVector(0,0,0);
  fSeaLevelRefractivity = 0;
   
  // set up table of particle names
  for (unsigned int j=0; j<nParticles; j++) {
    fParticles[gParticle_id[j]].name = gParticle_name[j];
  }
  
/*** Histogramming
  itsHeightBinsMeters = 1000.*meter; // one bin every km
  itsLengthBinsMeters = 0.1*meter; // one bin every 10 cm
  itsNumHeightBins = 100;
  itsNumLengthBins = 500;

  // initialize the std::vector
  itsHeightHistograms.reserve(itsNumHeightBins*itsNumLengthBins);
  for (int i=0; i<itsNumHeightBins*itsNumLengthBins; ++i)
    itsHeightHistograms.push_back(0);
*/

}



TPlotter::~TPlotter () {

  fParticles.clear();

  if (fGroundArea) {
    delete fGroundArea;
    fGroundArea = 0;
  }

  if (itsScenarioParams) {
    delete itsScenarioParams;   
    itsScenarioParams = 0;
  }

  if (fCORSIKA) {
    delete fCORSIKA;
    fCORSIKA = 0;
  }

}


void 
TPlotter::Welcome() 
   const 
{
  if (IsMainThread()) {
    cout << "\n"
         << " *******************************************************\n"
         << " **                                               \n"
         << " ** You are using the CoREAS V" << setprecision(2) << ProgramVersion << " radio simulation code\n"
         << " **                                               \n"
         << " ** PERFORMANCE OPTIMIZED VERSION v1.3           \n"
         << " **                                               \n"
         << " **  THIN/CURVED/SLANT/STACKIN/PRESHOWER: "
         << fThinning << "/" << fCurved << "/" << fSlant << "/" << fStackInput << "/" << fPreshower << "\n"
         << " **                                               \n"
         << " ** Please cite the following article:            \n"
         << " ** T. Huege, M. Ludwig, C.W. James, AIP Conf. Proc.\n"
         << " ** 1535, 128-132 (2013), doi:10.1063/1.4807534     \n"
         << " **                                               \n"
         << " *******************************************************\n\n"
         << endl;
  }
}

void
TPlotter::InitializeRadioSimulation()
{
  ostringstream fname;
  fname << "SIM" << setw (6) << setfill('0') << fRunNo; // read-in name constructed from run number
  ostringstream outfname;
  outfname << "SIM" << fFileName.substr(3,fFileName.size()-3); // write-out file name coming from CORSIKA (for parallel option)
  itsScenarioParams = new ScenarioParams(fDirectory+fname.str(), fDirectory+outfname.str(), fDirectory+fname.str());
  
  bool readok = itsScenarioParams->ReadFromFiles();

  if (readok)
  {
    // overwrite values with the correct ones from CORSIKA
    itsScenarioParams->itsShowerZenithAngle = fZenith*180./Pi;
    itsScenarioParams->itsShowerAzimuthAngle = fAzimuth*180./Pi;
    itsScenarioParams->itsPrimaryParticleType = fShowerPrimary;
    itsScenarioParams->itsPrimaryParticleEnergy = fShowerEnergy*1.e9;
    itsScenarioParams->itsMagneticFieldStrength = sqrt(fBx*fBx+fBz*fBz)/100.0;       // convert from microTesla to Gauss before writing out
    itsScenarioParams->itsMagneticFieldInclinationAngle = atan2(fBz,fBx)*180.0/M_PI; // convert from radians to degrees before writing out
    // temporarily delete values for depth of shower maximum and distance of shower maximum
    itsScenarioParams->itsDepthOfShowerMaximum = -1.0;
    itsScenarioParams->itsDistanceOfShowerMaximum = -1.0;

    if (IsMainThread()) {
      cout << " CoREAS" << fThreadName <<": The following parameters are used for the radio calculation\n\n";
      cout << (*itsScenarioParams); // write used parameter values into CORSIKA output log
      cout << "\n";
    }

    // now calculate the (radio) core position at the height of itsCoreCoordinateVertical
    itsShowerAxis = ThreeVector(cos(fAzimuth)*sin(fZenith),sin(fAzimuth)*sin(fZenith),-cos(fZenith));
    const double p = (itsScenarioParams->itsCoreCoordinateVertical-fObservationLevel)/cos(fZenith);
    fRadioCore = ThreeVector(0.0,0.0,fObservationLevel)-p*itsShowerAxis;
    
    //~ if (IsMainThread())
      //~ cout << " CoREAS radioCore: " << fRadioCore << "\n" << endl;
     
    const ThreeVector firstInteraction(fFirstInteractionX, fFirstInteractionY, fFirstInteractionZ);
    const double distFirstIntToCore = (firstInteraction-fRadioCore).GetLength();
    fCoreHitTime = fFirstInteractionTime+distFirstIntToCore/SpeedOfLight; // SpeedOfLight is in cgs units

    // check deviation of local zenith angle if OBSLEV is not equal to CoreCoordinateVertical within 10 meters
    if (fabs(fObservationLevel-itsScenarioParams->itsCoreCoordinateVertical) > 1000.)
    {
      if (IsMainThread())
      {
        const ThreeVector r1 = ThreeVector(0,0,-cRearth)-fRadioCore;
        const ThreeVector r2 = itsShowerAxis;
        const double Theta = acos(r1.GetDirection().DottedWith(r2.GetDirection()));
        cout << " CoREAS WARNING: OBSLEV is at " << fObservationLevel/m << " m while radio core is at " << itsScenarioParams->itsCoreCoordinateVertical/m << " m.\n";
        cout << " The local zenith angle at the radio core is thus " << Theta*180./M_PI << " degrees rather than the requested " << fZenith*180./M_PI << "degrees.\n";
        cout << " To avoid this mismatch, set OBSLEV the same as CoreCoordinateVertical!" << "\n\n";
      }
    }

    // initialize ground area
    fGroundArea = new GroundArea(itsScenarioParams->itsTimeResolution, 
      itsScenarioParams->itsTimeLowerBoundary, 
      itsScenarioParams->itsTimeUpperBoundary, 
      itsScenarioParams->itsAutomaticTimeBoundaries, 
      itsScenarioParams->itsResolutionReductionScale, 
      firstInteraction, 
      fFirstInteractionTime-fCoreHitTime,
      itsScenarioParams->itsAntennaPositions,
      fRadioCore, itsShowerAxis);

    long maxSampleLength = 0;
    // check that all antennas have the same number of samples, sometimes rounding errors lead to off-by-one problems which we will fix on the go, needs two passes
    for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); it != fGroundArea->GetActiveGroundElementPointers().end(); ++it)
      maxSampleLength = max(maxSampleLength, (*it)->GetNumSamples());
    for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); it != fGroundArea->GetActiveGroundElementPointers().end(); ++it)
      if ((*it)->GetNumSamples()!= maxSampleLength)
        (*it)->EnlargeSampleLengthTo(maxSampleLength);

    if (fObservationLevel > fGroundArea->GetLowestHeight()) {
      if (IsMainThread())
        cout << " CoREAS warning: lowest observation level is at " << fObservationLevel << " cm and thus above lowest antenna height at " << fGroundArea->GetLowestHeight() << " cm!\n\n";
    }

    fSeaLevelRefractivity = itsScenarioParams->itsGroundLevelRefractiveIndex-1.0;

  }
  else
  {    
    // error message is being printed out by ScenarioParams class
    exit (11);
  }

  // set up the internal refractive index model for vertical geometry, either using the Gladstone-Dale law or a GDAS-based refractivity
  // if GDAS-based, the fRefractivityTable-vector has been pre-filled at this stage
  size_t numberOfPoints = fRefractivityTable.size();
  double  minimumHeight, maximumHeight;
  const double safetyMargin = 1.0 * km; // provide a safety margin because particles can in fact propagate beyond the observation level
  if (numberOfPoints == 0)
  {
    // set up Gladstone-Dale law tables of refractivity versus height

    if (IsMainThread())
      cout << "Initializing Gladstone-Dale lookup table for refractivity.\n";

    minimumHeight = -safetyMargin;
    if (fObservationLevel-safetyMargin < minimumHeight )
      minimumHeight = fObservationLevel-safetyMargin;

    if (fGroundArea->GetLowestHeight()-safetyMargin < minimumHeight )
      minimumHeight = fGroundArea->GetLowestHeight()-safetyMargin;

    maximumHeight = 112.82920 * km; // ToDo: get value from CORSIKA
    fDistanceIncrement = 1.0 * meter;
    numberOfPoints = (maximumHeight-minimumHeight) / fDistanceIncrement + 1;

    fRefractivityTable.resize(numberOfPoints);
    fHeightTable.resize(numberOfPoints);

    double rho_0 = rhof_(0.);
    for (size_t i = 0; i < numberOfPoints; i++)
    {
      fHeightTable[i] = minimumHeight + i * fDistanceIncrement;
      double rho = rhof_(fHeightTable[i]);
      fRefractivityTable[i] = fSeaLevelRefractivity * rho / rho_0;
    }
  }
  else
  {
    // use refractivity table pre-filled from GDAS data

    if (IsMainThread())
      cout << "Using refractivity lookup table from GDAS-based file.\n";

    minimumHeight = fHeightTable.front();
    maximumHeight = fHeightTable.back();
    fDistanceIncrement = (maximumHeight - minimumHeight) / (fHeightTable.size()-1);

    if ((minimumHeight > (fGroundArea->GetLowestHeight() - safetyMargin))
    ||  (minimumHeight > (fObservationLevel - safetyMargin)))
    {
      if (IsMainThread()) {
        cout << "\n\n"
             << " ###################################\n"
             << " Minimum height in GDAS-based refractivity profile is: " << minimumHeight / km << " km.\n"
             << " The minimum height has to be " << safetyMargin / km << " km lower than both the lowest antenna and the observation level, because particles can propagate beyond the observation level.\n"
             << " ###################################\n\n"
             << endl;
      }
      exit (11);
    };
  }

  // tabulate the integrated refractivity as a function of height, possibly half-bin shifted
  fIntegratedRefractivityTable.resize(fRefractivityTable.size());
  fIntegratedRefractivityTable[0] = (fRefractivityTable[0]) * fDistanceIncrement;
  for (size_t i = 1; i < fIntegratedRefractivityTable.size(); i++)
  {
    fIntegratedRefractivityTable[i] = fIntegratedRefractivityTable[i-1] + (fRefractivityTable[i]) * fDistanceIncrement;
  }

  if (IsMainThread())
  {
    cout << " Tabulated atmosphere inside CoREAS has following parameters" << std::endl
         << "     Minimum height: " << minimumHeight / m << " m" <<  std::endl
         << "     Maximum height: " << maximumHeight / m << " m" << std::endl
         << "         Resolution: " << fDistanceIncrement / m << " m" << std::endl
         << "   Number of points: " << numberOfPoints << "\n" << std::endl;
  }

  // set up refractivity tables for zenith angles up to 89 degrees corresponding to tan(theta) = 0..60
  if ((fCurved) && (fZenith>=fTabulationThresholdZenithAngle))
  {
    if (IsMainThread())
    {
    cout << " Setting up tabulated integrated refractivities for curved geometries." << std::endl
         << " Number of tan(theta) bins: " << fNumTanThetaBins << std::endl
         << " Spacing of tan(theta) bins: " << fTanThetaBinSpacing << std::endl;

    }

    // reserve memory for tanThetaTables
    fIntegratedRefractivityInCurvedAtmosphereTables.resize(fNumTanThetaBins);
    fDistanceTableIncrements.reserve(fNumTanThetaBins);
    fDistanceTableMinDistances.reserve(fNumTanThetaBins);
    for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); it != fGroundArea->GetActiveGroundElementPointers().end(); ++it)
    {
      (*it)->GetEarthIntersectionDistances().reserve(fNumTanThetaBins);
      (*it)->GetEarthIntersectionAngleBins().reserve(fNumTanThetaBins);
      (*it)->SetHeight(GetHeightAtPosition((*it)->GetPosition()));
    }
    
    const unsigned long numDistanceBins = fRefractivityTable.size()/4; // define number of bins used for atmospheric tables, quarter of original table length
    
    for (unsigned int i = 0; i<fNumTanThetaBins; ++i)
    {
      const double tanTheta = i*fTanThetaBinSpacing;
      const double theta = atan(tanTheta);
      fIntegratedRefractivityInCurvedAtmosphereTables.at(i).resize(numDistanceBins);

      // calculate distance corresponding to edge of atmosphere at this zenith angle
      double x1, y1, x2, y2, maxDistanceAnalytic;
      double y0 = cRearth;
      const double deduction = 100.; // lower maxHeight by 1 m as to not cause exceptions
      int numPoints = GetIntersectionsLineCircle(cRearth+maximumHeight-deduction, -tanTheta, y0*tanTheta, x1, y1, x2, y2);
      if (numPoints == 2)
      {
        double d1 = sqrt((y1-y0)*(y1-y0)+x1*x1);
        double d2 = sqrt((y2-y0)*(y2-y0)+x2*x2);
        //ThreeVector edgePosition(x2,0,y2-cRearth);
        maxDistanceAnalytic = min(d1,d2);
        //std::cout << x2 << "\t" << y2 << "\t" << maxDistanceAnalytic << "\t" << edgePosition << "\t" << GetHeightAtPosition(edgePosition) << "\n";
      }
      else
      {
        std::cout << "ERROR! numPoints = " << numPoints << " in calculation of edge of atmosphere.\n";
        exit(11);
      }

      // calculate distance to Earth and impact zenith angle for heights of our given antennas at this zenith angle
      for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); it != fGroundArea->GetActiveGroundElementPointers().end(); ++it)
      {
        double currentHeight = (*it)->GetHeight();
        double edgeDistance;
        ThreeVector edgePosition;
        y0 = cRearth+currentHeight;
        numPoints = GetIntersectionsLineCircle(cRearth, -tanTheta, y0*tanTheta, x1, y1, x2, y2);
        if (numPoints == 2)
        {
          double d1 = sqrt((y1-y0)*(y1-y0)+x1*x1);
          double d2 = sqrt((y2-y0)*(y2-y0)+x2*x2);
          if (d2 < d1)
          {
            edgePosition = ThreeVector(x2,0,y2-cRearth);
            edgeDistance = d2;
          }
          else 
          {
            edgePosition = ThreeVector(x1,0,y1-cRearth);
            edgeDistance = d1;
          }
          // now calculate and note down Earth intersaction angle
          const ThreeVector r1 = ThreeVector(0,0,-cRearth)-edgePosition;
          const ThreeVector r2 = edgePosition-ThreeVector(0,0,currentHeight);
          const double cosTheta = r1.GetDirection().DottedWith(r2.GetDirection());
          const double tanEdgeAngle = sqrt(1.-cosTheta*cosTheta)/cosTheta; // see https://www.rapidtables.com/math/trigonometry/arccos/tan-of-arccos.html
          (*it)->GetEarthIntersectionAngleBins().push_back(static_cast<unsigned int>(tanEdgeAngle/fTanThetaBinSpacing+0.5));
        }
        else
        {
          // tangent at Earth surface or not intersecting with Earth surface, so intersection angle is 90° and thus definitely outside tabulation
          edgeDistance = 0;
          (*it)->GetEarthIntersectionAngleBins().push_back(fNumTanThetaBins+1); // i.e. not tabulated
        }
        (*it)->GetEarthIntersectionDistances().push_back(edgeDistance);
      }      
      
      const ThreeVector travelDirection(sin(theta),0,cos(theta));
      const double increment = (maxDistanceAnalytic-minimumHeight)/(numDistanceBins-1);
      fDistanceTableIncrements.push_back(increment);
      fDistanceTableMinDistances.push_back(minimumHeight);

      double integratedRefractivity = 0;
      for (unsigned long j = 0; j<numDistanceBins; ++j)
      {
        double distance = minimumHeight+j*increment;
        double localHeight = GetHeightAtPosition(distance*travelDirection);
        double localRefractivity = InterpolateInPlanarAtmosphere(localHeight, fRefractivityTable);
        integratedRefractivity += localRefractivity*increment;
        fIntegratedRefractivityInCurvedAtmosphereTables.at(i)[j] = integratedRefractivity;
      }
      
    }
    cout << "\n";
  }

}


void TPlotter::SetRunHeader(const crs::MRunHeader &header) 
{
  ostringstream ssTmp; 
  ssTmp.str(""); ssTmp << header.GetVersion(); fCorsikaVersion = ssTmp.str();
  ssTmp.str(""); ssTmp << setprecision(2) << ProgramVersion; fCoastVersion = ssTmp.str();
  
  SetRun ((int)header.GetRunID ());

  if ((int)header.GetNumberOfShowers() != 1) {
    cout << "\n\n\nCoREAS" << fThreadName <<": Error, NSHOW is set to a value of " << (int)header.GetNumberOfShowers() <<". CoREAS only supports NSHOW = 1.\n\n\n" << endl;
    exit(11);
  }
  
}

void 
TPlotter::SetShowerHeader(const crs::MEventHeader &header) 
{
  if (header.GetSkimmingIncidence()) {
    fSkimming = true;
    fSkimmingAltitude = header.GetSkimmingAltitude() * cm;
    
    if (IsMainThread()) {
      cout << " CoREAS" << fThreadName <<": Detected skimming geometry with impact=" 
          << fSkimmingAltitude/km
          << " km" << endl;
    }
  }
  
  // to flag the track of the primary particle
  fPrimaryTrack = true;
  fFirstInteractionX = 0;
  fFirstInteractionY = 0;
  fFirstInteractionZ = 0;
  
  SetEvent(header.GetEventNumber());
  fShowerPrimary = (int)header.GetParticleId();
  
  SetShowerAxis(header.GetTheta()*rad, header.GetPhi()*rad);
  fShowerEnergy = header.GetEnergy()*GeV;
  
  fHeightFirstInt = header.GetZFirst()*cm;
  fObservationLevel = header.GetObservationHeight(header.GetNObservationLevels()-1)*cm;
  
  if (fHeightFirstInt<0) {
    if (fVerbosityLevel>=10) {
      cout << " CoREAS" << fThreadName <<": height of first interaction is NEGATIVE, corrected." << endl;
    }
    fHeightFirstInt = fabs(fHeightFirstInt);
  }

  // store magnetic field configuration for later access at writeout of .reas file
  fBx = header.GetBx();
  fBz = header.GetBz();

  if (IsMainThread()) {
    cout  << " CoREAS" << fThreadName <<": height of first interaction " << fHeightFirstInt/m << "m\n"
          << " CoREAS" << fThreadName <<": observation level " << fObservationLevel/m  << "m\n"
          << endl;
  }

}


void 
TPlotter::SetShowerTrailer(const crs::MEventEnd &trailer) 
{
  fXmax = trailer.GetXmax()*g/cm/cm;
  if (!fSlant && fCosZenith!=0) 
    fXmax /= fCosZenith;
  
  if (IsMainThread()) {
    cout << " CoREAS" << fThreadName <<": --------- shower end ------ \n"
         << " CoREAS" << fThreadName <<": xmax: " << fXmax/g*cm*cm << "\n"
         << endl;
  }
}


void 
TPlotter::Init() 
{
  if (!fSkimming) 
    fCORSIKA = new TCorsika(fZenith, 
          fObservationLevel,
          fSlant, fCurved);
  else 
    fCORSIKA = new TCorsika(fSkimmingAltitude,
          fObservationLevel);
  
  if (fVerbosityLevel > 0) 
    fCORSIKA->SetVerbosityLevel(fVerbosityLevel);

  if( fCORSIKA == 0 )
    cout << "Warning fCORSIKA not defined" << endl;
  else
    fCORSIKA->SetHeightOfFirstInteraction(fHeightFirstInt);

  if (fCurved && (not fSlant)) { // curved without slant is not supported in CoREAS
    cout << "\n\n"
         << " ###################################\n"
         << "   CURVED without SLANT is not supported   \n"
         << "    please switch on the SLANT option \n"
         << " ###################################\n\n"
         << endl;
    exit (11);
  }
  
  //~ if (fVerbosityLevel > 1 && fCurved) {
    //~ if (IsMainThread())
      //~ fCORSIKA->DumpAtmosphereTable();
  //~ }
}

void 
TPlotter::Write() 
{
  if (fVerbosityLevel >= 2) {
    cout << " CoREAS" << fThreadName <<": Write " << endl;
  }

  // collect results differently for parallelized and unparallelized CoREAS
  
  #if __PARALLELIB__
  // do nothing, TPlotter object is recycled for next sub-shower until CommunicateNodeResults() is called
  #else
  // only write out results if radio simulation was initialized
  if (itsScenarioParams)
  {
    // Write some values into the ScenarioParams and then output the file
    itsScenarioParams->itsDepthOfShowerMaximum = fXmax/g*cm*cm;
    
    // use CORSIKA table to calculate distance to shower maximum
    //~ if ((fXmax > 0) && (fXmax < fCORSIKA->GetMaximumTabulatedDepth()))
      //~ itsScenarioParams->itsDistanceOfShowerMaximum = (fCORSIKA->GetDistanceOfPoint(fRadioCore.GetX(), fRadioCore.GetY(), fRadioCore.GetZ())-fCORSIKA->GetDistanceOfSlantDepth(fXmax))/cm;

    // do our own integration to calculate distance to shower maximum
    const double stepSize = 1.0 * m;
    const ThreeVector arrivalAxis = -1 * itsShowerAxis;
    double integratedGrammage = 0;
    double localHeight = 0;
    long i = 0;
    
    // first calculate total grammage of atmosphere, might be available elsewhere but make sure it is precise enough by calculating it ourselves
    while (localHeight < 112.82920 * km) // ToDo: get value from CORSIKA
    { 
      localHeight = GetHeightAtPosition(fRadioCore+stepSize*(i+0.5)*arrivalAxis);
      integratedGrammage += stepSize*rhof_(localHeight);
      ++i;
    }
    
    // now count distance to Xmax
    const double xMaxGrammage = integratedGrammage-fXmax/g*cm*cm;
    integratedGrammage = 0;
    i = 0;
    while (integratedGrammage < xMaxGrammage) {
      localHeight = GetHeightAtPosition(fRadioCore+stepSize*(i+0.5)*arrivalAxis);
      integratedGrammage += stepSize*rhof_(localHeight);
      ++i;
    }
    double distanceToXmax = i*stepSize;
    itsScenarioParams->itsDistanceOfShowerMaximum = distanceToXmax/cm;

    // write out data to .reas file
    itsScenarioParams->WriteToFiles();                                                                                             

    // output direct radio results
    ostringstream fullPath;
    fullPath << fDirectory << "SIM" << fFileName.substr(3,fFileName.size()-3) << "_coreas";
    ostringstream relativePath;
    relativePath << "SIM" << fFileName.substr(3,fFileName.size()-3) << "_coreas";
    fGroundArea->WriteResults(fullPath.str(), relativePath.str());
  }
  #endif

/*** Histogramming
  ostringstream histofilename;
  histofilename << fDirectory << "SIM" << fFileName.substr(3,fFileName.size()-3) << ".histo";
  ofstream histofile(histofilename.str().c_str());
  for (long i=0; i<itsNumHeightBins; ++i) {
    histofile << "# height bin " << i << " from " << i*itsHeightBinsMeters/meter << " to " << (i+1)*itsHeightBinsMeters/meter << " meters\n";
    for (long j=0; j<itsNumLengthBins; ++j) {
      histofile << j*itsLengthBinsMeters/meter << " " << itsHeightHistograms.at(i*itsNumLengthBins+j) << "\n";
    }
    histofile << "\n\n"; // two blank lines for gnuplot-index
  }
  histofile.close();
*/
 
}

extern int dbg_infos;
void
TPlotter::CommunicateNodeResults()
{
  // check if the simulation has actually run, otherwise do nothing
  if (!itsScenarioParams)
    return;

#if __PARALLELIB__
#define MASTER 0
  int rank, num_proc, mpierr;

  // determine the paths for writeout - chopping off the thread name at the end!!!
  ostringstream fullPath;
  fullPath << fDirectory << "SIM" << fFileName.substr(3,6) << "_coreas";
  ostringstream relativePath;
  relativePath << "SIM" << fFileName.substr(3,6) << "_coreas";

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_proc);

  if (rank == 1)
  {
    // for one of the nodes only write out the generic files (.reas, .bins) and create the subdirectory for the antenna data files
    itsScenarioParams->itsDepthOfShowerMaximum = -1.;    // not known in case of parallel simulation
    itsScenarioParams->itsDistanceOfShowerMaximum = -1.; // not known in case of parallel simulation
    itsScenarioParams->WriteToFiles();
    fGroundArea->WriteGenericFiles(fullPath.str(),relativePath.str());
  }

  MPI_Barrier(MPI_COMM_WORLD); // Barrier: to let slave with rank 1 create the forder for Storing CoREAS simulation data

  // determine length of the antenna traces - number of entries has to be the same for all antennas!
  long numAntennas = fGroundArea->GetActiveGroundElementPointers().size();
  long length = fGroundArea->GetActiveGroundElementPointers().front()->GetNumSamples();

  // cross-check that trace lengths are same for all antennas BEFORE doing any time-consuming communication
  for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); it != fGroundArea->GetActiveGroundElementPointers().end(); ++it)
  {
    if ((*it)->GetNumSamples() != length)
    {
      cout << " CoREAS" << fThreadName <<": Error, length of traces not identical for all antennas! Aborting ..." << endl;
      MPI_Abort(MPI_COMM_WORLD,111);
    }
  }
  long bufferlength = length;

  if (rank == 1)
  {
    long AntParams[2] = {numAntennas, length};
    mpierr=MPI_Send(AntParams,2,MPI_LONG,MASTER,6666,MPI_COMM_WORLD);
    if (mpierr!= MPI_SUCCESS)
    {
        cout << " CoREAS" << fThreadName <<": Error, In MPI_Send communication about number of Antennas\n";
        MPI_Abort(MPI_COMM_WORLD,mpierr);
    }
  }

  // distribute antennas among slaves and master
  int minAntPerProc = numAntennas / num_proc;
  int remAnt = numAntennas % num_proc;
  int myAntNum = minAntPerProc + ((rank < remAnt)?1:0);
  int masterAntNum = minAntPerProc + ((remAnt > 0)?1:0);

  char statfile[255],numstr[9];
  FILE* file;
  if (dbg_infos)
  {  //Record status of Master to Slave communication in statistic file
      strcpy(statfile,fullPath.str().c_str());//statdir
      strcat(statfile,"/mpiid-");
      sprintf(numstr, "%d", rank);
      strcat(statfile,numstr);
      strcat(statfile,"-antenna.txt");
      file = fopen(statfile,"a");
      fprintf(file,"\n[%f] FINALIZE: %ld antennas with length of the first =%ld",MPI_Wtime(),numAntennas,length);
      fprintf(file,"\n\n[%f] MY NUMBER OF ANTENNAS TO COLLECT IS: %d \n",MPI_Wtime(),myAntNum);
      fclose(file);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //allocate memory for result array after knowing msgsize ---
  double * antenna = (double *)calloc(bufferlength*3, sizeof(double));

  // now send SOME antenna data to MPI master
  long antnum = 1;
  if(rank == 1)
  {
      for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); (antnum <= masterAntNum) && ( it != fGroundArea->GetActiveGroundElementPointers().end() ); it+=num_proc, ++antnum)
      {
          char filename[300];

          // get the file name for the given antenna ??? antnum must be less or equal to numAntennas
          ostringstream rawss;
          rawss << fullPath.str() << "/raw_" << (*it)->GetObserverName() << ".dat";

          strcpy(filename, rawss.str().c_str());
          if (dbg_infos)
          {
              file =fopen(statfile,"a");
              fprintf(file,"\n[%f]: File Name of antenna #%ld is %s",MPI_Wtime(),1+(antnum-1)*num_proc,filename);
              fclose(file);
          }

          double timeparams[2];
          timeparams[0] = (*it)->GetTimeLowerBoundary();
          timeparams[1] = (*it)->GetTimeResolution();

          // transfer data to the MPI master
          MPI_Send(filename, 300, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD);  // File name for antenna data antnum as TAG only when submitting as a joint structure - be careful when editing length of file name
          if (dbg_infos)
          {
              file =fopen(statfile,"a");
              fprintf(file,"\n[%f]: Successfully SEND File Name of antenna #%ld",MPI_Wtime(), antnum);
              fclose(file);
          }
          MPI_Send(timeparams, 2, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);  // File name for antenna data antnum as TAG only when submitting as a joint structure - be careful when editing length of file name
          if (dbg_infos)
          {
              file =fopen(statfile,"a");
              fprintf(file,"\n[%f]: Successfully SEND timing parameters of antenna #%ld",MPI_Wtime(), antnum);
              fclose(file);
          }
      }
  }
  if (dbg_infos)
  {
      file =fopen(statfile,"a");
      fprintf(file,"\n[%f]: BARRIER waiting to slave #1 to sent antenna filenames and master finish receiving",MPI_Wtime());
      fclose(file);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  typedef struct antenna_s
      {
      char filename[300];//text for file name of antenna-identifier
      double time0;//time lower bound
      double dt; // time step
      double *x;//X coordinate buffer
      double *y;//Y coordinate buffer
      double *z;//Z coordinate buffer
      } antennastruc;
  antennastruc * antennas = 0;

  if(myAntNum > 0)
  {
      antennas=(antennastruc *)malloc((myAntNum)*sizeof(antennastruc));//array of structures
  }

  for(int j=0; j<myAntNum; j++)
  {
      antennas[j].x = (double*)malloc(length*sizeof(double));
      antennas[j].y = (double*)malloc(length*sizeof(double));
      antennas[j].z = (double*)malloc(length*sizeof(double));
  }
  antnum = 1;
  double * recv_buffer = (double *)malloc(bufferlength*3*sizeof(double));

  int localnum=0;
  for (vector<GroundElement*>::const_iterator it = fGroundArea->GetActiveGroundElementPointers().begin(); it != fGroundArea->GetActiveGroundElementPointers().end(); ++it, ++antnum)
  {
    // copy the data into the common buffer
    for (int j=0;j<length;j++) // x,y,z arrays of length size
    {
        (*it)->CopySample(j, antenna[j], antenna[j+bufferlength], antenna[j+bufferlength*2]);
    }
    int gather_rank = (antnum-1)%num_proc;
    mpierr=MPI_Reduce(antenna, recv_buffer, bufferlength*3, MPI_DOUBLE, MPI_SUM, gather_rank, MPI_COMM_WORLD);

    if (dbg_infos)
    {
      file =fopen(statfile,"a");
      fprintf(file,"\n[%f]: Successfully sent data of antenna #%ld to RANK %d",MPI_Wtime(), antnum, gather_rank);
      fclose(file);
    }
    if(rank == gather_rank)
    {
        if (dbg_infos)
        {
            file =fopen(statfile,"a");
            fprintf(file,"\n[%f]: RECEIVED antenna #%d (global number %ld)",MPI_Wtime(),localnum+1,antnum);
            fclose(file);
        }
        // get the file name for the given antenna ??? antnum must be less or equal to numAntennas
        ostringstream rawss;
        rawss << fullPath.str() << "/raw_" << (*it)->GetObserverName() << ".dat";

        strcpy(antennas[localnum].filename, rawss.str().c_str());

        if (dbg_infos)
        {
            file =fopen(statfile,"a");
            fprintf(file,"\n[%f]: File Name of MY antenna  #%d (global number %ld) is %s",MPI_Wtime(),localnum+1,antnum,antennas[localnum].filename);
            fclose(file);
        }

        antennas[localnum].time0 = (*it)->GetTimeLowerBoundary();
        antennas[localnum].dt = (*it)->GetTimeResolution();

        for(int k=0;k<length;k++)
        {
            antennas[localnum].x[k] =  recv_buffer[k];
            antennas[localnum].y[k] =  recv_buffer[k+bufferlength];
            antennas[localnum].z[k] =  recv_buffer[k+bufferlength*2];
        }
        localnum++;
    }
  }
  free(recv_buffer);

  for(int j=0; j<myAntNum; j++)
  {
      FILE *antfile;
      antfile=fopen(antennas[j].filename,"w");
      for(int k=0;k<length;k++)
      {
          fprintf(antfile,"%.12e\t%.12e\t%.12e\t%.12e\n",antennas[j].time0 + antennas[j].dt*k,antennas[j].x[k],antennas[j].y[k],antennas[j].z[k]);
      }
      fclose(antfile);
  }

  for(int j=0; j<myAntNum; j++)
  {
      free(antennas[j].x);
      free(antennas[j].y);
      free(antennas[j].z);
  }

  free(antenna);
#endif
}


//~ inline
//~ void 
//~ TPlotter::Rotate(double x, double y, double z,
     //~ double &sx, double &sy, double &sz,
     //~ int inverse) const
//~ {

  //~ // rotate around z by azimuth
  //~ double sx_ =             x*fCosAzimuth + inverse * y*fSinAzimuth;
  //~ double sy_ = - inverse * x*fSinAzimuth +           y*fCosAzimuth; 
  //~ double sz_ =   z;
  
  //~ // rotate around y` by zenith
  //~ double sx__ =             sx_*fCosZenith - inverse * sz_*fSinZenith;
  //~ double sy__ = sy_;
  //~ double sz__ = + inverse * sx_*fSinZenith +           sz_*fCosZenith; 

  //~ // rotate back around z`` by -azimuth 
  //~ sx =             sx__*fCosAzimuth - inverse * sy__*fSinAzimuth;
  //~ sy = + inverse * sx__*fSinAzimuth +           sy__*fCosAzimuth; 
  //~ sz =   sz__;

//~ }


void 
TPlotter::AddTrack(const crs::CParticle &pre, 
                   const crs::CParticle &post) 
{  

  /*
    Check whether we are dealing with primary particle. If so,
    skip the track of the primary particle, which is in a different 
    reference system as the shower. Otherwise set first interaction.
  */

  if (fPrimaryTrack) {
    crs::CParticleFortranPtr pptr;
    fPrimaryTrack = !prminfo_(pptr);          // the flag is false if the primary info is not yet available because we are still following the primary
    if (fPrimaryTrack) {
      if (fVerbosityLevel >= 2) {
        cout << " CoREAS" << fThreadName <<": Primary track, skipping  " << endl;
      }    
      return;
    }
    else {
      SetFirstInteraction(pptr->x, pptr->y, pptr->z, pptr->depth*g/cm/cm, pptr->time);
      if (fVerbosityLevel >= 2) {
        cout << " CoREAS" << fThreadName <<": Primary track end-pos: x=" << fFirstInteractionX/m 
             << " y=" << fFirstInteractionY/m 
             << " z=" << fFirstInteractionZ/m << endl;
      }
    }
  }

  if (fVerbosityLevel >= 9) {
  double dX = post.x-pre.x; 
  double dY = post.y-pre.y; 
  double dZ = post.z-pre.z;
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);
  double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
  double postGamma = post.energy/gParticleMass[abs((int)post.particleId)-1];
  double preBetaValue = sqrt(1.-1./(preGamma*preGamma));
  double postBetaValue = sqrt(1.-1./(postGamma*postGamma));  
  cout << " CoREAS" << fThreadName <<":   PRE> "
   << " id=" << (int)pre.particleId
   << " E=" << pre.energy
   << " w=" << pre.weight
   << " x=" << pre.x/km << "km"
   << " y=" << pre.y/km << "km"
   << " z=" << pre.z/km << "km"
   << " t=" << pre.time*s/ns << "ns"
   << " X=" << pre.depth << "g/cm^2"
/* << " Px=" << pre.Px << ""
   << " Py=" << pre.Py << ""
   << " Pz=" << pre.Pz << "" */
   << " v= " <<  length/(post.time*s/ns-pre.time*s/ns)/crs::cSpeedOfLight 
   << " prebeta = " << preBetaValue
   << " postbeta = " << postBetaValue
   << " ratio = " << length/(post.time*s/ns-pre.time*s/ns)/crs::cSpeedOfLight/preBetaValue
   << endl;
  }
  
  if (fVerbosityLevel >= 10) {
/*double dX = post.x-fFirstInteractionX; 
  double dY = post.y-fFirstInteractionY; 
  double dZ = post.z-fFirstInteractionZ; 
  double length = sqrt (dX*dX + dY*dY + dZ*dZ);*/
  cout << " CoREAS" << fThreadName <<":  POST> "
   << " id=" << (int)post.particleId
   << " E=" << post.energy
   << " w=" << post.weight
   << " x=" << post.x/km << "km"
   << " y=" << post.y/km << "km"
   << " z=" << post.z/km << "km"
   << " t=" << post.time*s/ns << "ns"
   << " X=" << post.depth << "g/cm^2"
// << " v= " <<  length/(post.time*s/ns-fFirstInteractionTime*s/ns)/crs::cSpeedOfLight
   << endl;
  }
    
  /* Skip the particle if pre and post times are identical */
  if (pre.time-post.time == 0.0)
    return;  
 
  const int particleId = abs((int)pre.particleId);
  
  // check if particles to be counted at all
  if (!fParticles.count(particleId)) 
    return;

  // calculate direct radio
  // find a better way to specify charge signs
  double particleCharge = gParticleCharge[abs((int)pre.particleId)-1];
      
  ThreeVector prePosition(pre.x, pre.y, pre.z);
  ThreeVector postPosition(post.x, post.y, post.z);

/*** Histogramming

  double stepLength = (postPosition-prePosition).GetLength();
  double stepHeight = pre.z;  // assume flat geometry
  if ((stepLength < itsNumLengthBins*itsLengthBinsMeters)
     && (stepLength >= 0)
     && (stepHeight < itsNumHeightBins*itsHeightBinsMeters)
     && (stepHeight >= 0))
  {
    // calculate bin - height bins are outer loop
    long heightBin = stepHeight/itsHeightBinsMeters;
    long stepBin = stepLength/itsLengthBinsMeters;
    itsHeightHistograms.at(heightBin*itsNumLengthBins+stepBin)++;
  }
*/

  double preIndex = GetRefractiveIndexAtPosition(prePosition);
  double postIndex = GetRefractiveIndexAtPosition(postPosition);
      
  double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
  double postGamma = post.energy/gParticleMass[abs((int)post.particleId)-1];
      
  double preBetaValue = sqrt(1.-1./(preGamma*preGamma));
  double postBetaValue = sqrt(1.-1./(postGamma*postGamma));
    
  ThreeVector currDirection(post.x-pre.x, post.y-pre.y, post.z-pre.z);
  currDirection*=1./currDirection.GetLength();
      
//ThreeVector preBeta = preBetaValue*currDirection;     // neglects curvature along track from pre to post!
//ThreeVector postBeta = postBetaValue*currDirection;   // neglects curvature along track from pre to post!

  double corrBetaValue = (postPosition-prePosition).GetLength()/(SpeedOfLight*(post.time-pre.time));
  ThreeVector corrBeta = corrBetaValue * currDirection;
  const double velocityWeight = 0.5*(postBetaValue+preBetaValue)/corrBetaValue;

/*cout << "--------------------\n";
  cout << "pre.particleId: " << pre.particleId << "\n";
  cout << "post.particleId: " << post.particleId << "\n";
  cout << "particleCharge: " << particleCharge << "\n";
  cout << "pre.energy: " << pre.energy << "\n";
  cout << "post.energy: " << post.energy << "\n";
  cout << "preGamma: " << preGamma << "\n";
  cout << "postGamma: " << postGamma << "\n";
  cout << "preBetaValue: " << preBetaValue << "\n";
  cout << "postBetaValue: " << postBetaValue << "\n";
  cout << "prePosition: " << prePosition.GetX() << ", " << prePosition.GetY() << ", " << prePosition.GetZ() << "\n";     
  cout << "postPosition: " << postPosition.GetX() << ", " << postPosition.GetY() << ", " << postPosition.GetZ() << "\n";      
  cout << "preBeta: " << preBeta.GetX() << ", " << preBeta.GetY() << ", " << preBeta.GetZ() << "\n";     
  cout << "postBeta: " << postBeta.GetX() << ", " << postBeta.GetY() << ", " << postBeta.GetZ() << "\n";      
  cout << "corrBeta: " << corrBeta.GetX() << ", " << corrBeta.GetY() << ", " << corrBeta.GetZ() << "\n";      
*/

  const double approxthreshold = 1.0e-3; // set threshold for application of ZHS-like approximation
  ThreeVector startE;
  ThreeVector endE;
  
  const vector<GroundElement*>& groundElements = fGroundArea->GetActiveGroundElementPointers();
  for (vector<GroundElement*>::const_iterator currGroundElement = groundElements.begin();
      currGroundElement != groundElements.end(); ++currGroundElement)  
  {
    ThreeVector obsPosition = (*currGroundElement)->GetPosition();

    // check if particle contribution should be skipped due to special mode of observer bin
    if ((*currGroundElement)->GetMode() != Normal)
    {
      // check whether to apply cuts in particle gamma
      if ((*currGroundElement)->GetMode() == Gamma)
      {
        double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
        if ((preGamma < (*currGroundElement)->GetModeDataLow()) || (preGamma >= (*currGroundElement)->GetModeDataHigh())) {
          continue;
        }
      }
      // check whether to apply cuts in particle height asl
      else if ((*currGroundElement)->GetMode() == Height)
      {
        double height = GetHeightAtPosition(ThreeVector(pre.x, pre.y, pre.z));
        if ((height < (*currGroundElement)->GetModeDataLow())
           || (height >= (*currGroundElement)->GetModeDataHigh())) {
          continue;
        }
      }
      // check whether to apply cuts in particle slant depth
      else if ((*currGroundElement)->GetMode() == SlantDepth)
      {
        if ((pre.depth < (*currGroundElement)->GetModeDataLow())
           || (pre.depth >= (*currGroundElement)->GetModeDataHigh())) {
          continue;
        }
      }
      // check whether to apply cuts in distance from shower plane through the core
      else if ((*currGroundElement)->GetMode() == Distance)
      {
        //~ double distance = fCORSIKA->GetDistanceOfPoint(fRadioCore.GetX(), fRadioCore.GetY(), fRadioCore.GetZ())-fCORSIKA->GetDistanceOfSlantDepth(pre.depth*g/cm/cm);
        double distance = (fRadioCore-prePosition).DottedWith(itsShowerAxis);
        if ((distance < (*currGroundElement)->GetModeDataLow())
           || (distance >= (*currGroundElement)->GetModeDataHigh())) {
          continue;
        }
      }
    }    
        
    ThreeVector preN = obsPosition-prePosition;
    double preDistance = preN.GetLength();
    preN *= 1./preDistance;
    double preT = pre.time-fCoreHitTime;
    double preResT = preT+GetEffectiveRefractiveIndexBetween(*currGroundElement,prePosition)*preDistance/SpeedOfLight;
        
    ThreeVector postN = obsPosition-postPosition;
    double postDistance = postN.GetLength();
    postN *= 1./postDistance;
    double postT = post.time-fCoreHitTime;
    double postResT = postT+GetEffectiveRefractiveIndexBetween(*currGroundElement,postPosition)*postDistance/SpeedOfLight;

    // skip particle if both start- and end-contributions lie outside time-window to be recorded
    if (((preResT < (*currGroundElement)->GetTimeLowerBoundary()) || (preResT >= (*currGroundElement)->GetTimeUpperBoundary())) &&
       ((postResT < (*currGroundElement)->GetTimeLowerBoundary()) || (postResT >= (*currGroundElement)->GetTimeUpperBoundary())))
      continue;
        
    double preDoppler = (1.0 - preIndex*corrBeta.DottedWith(preN));
    double postDoppler = (1.0 - postIndex*corrBeta.DottedWith(postN));

    // apply special treatment for refractive index, this is motivated by a finite observation time resolution and a ZHS-like limit calculation
    if ( (fSeaLevelRefractivity > 0.0) && ((fabs(preDoppler)<approxthreshold) || (fabs(postDoppler)<approxthreshold)) )
    {
      ThreeVector midPosition = 0.5*(prePosition+postPosition);
      ThreeVector midN = obsPosition-midPosition;
      double midDistance = midN.GetLength();
      midN *= 1./midDistance;
      double midIndex = GetRefractiveIndexAtPosition(midPosition);
      double midDoppler = (1.0 - midIndex*corrBeta.DottedWith(midN));

      // check if midDoppler has become zero because of numerical limitations
      if (midDoppler == 0)
      {
        // redo calculation with higher precision
        long double indexL = midIndex;
        long double betaX = corrBeta.GetX();
        long double betaY = corrBeta.GetY();
        long double betaZ = corrBeta.GetZ();
        long double nX = midN.GetX();
        long double nY = midN.GetY();
        long double nZ = midN.GetZ();
        long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
        midDoppler = doppler;
      }

      startE = pre.weight * particleCharge * UnitCharge * ( midN.CrossedWith( midN.CrossedWith(corrBeta))) /
               (SpeedOfLight * midDistance * midDoppler);      
                               
      endE = -1.0 * startE;
        
      // force deltaT to be proportional to the Doppler factor (i.e. contribution = constant * track length), analogous to ZHS
      const double midResT = 0.5*(preT+postT)+GetEffectiveRefractiveIndexBetween(*currGroundElement,midPosition)*midDistance/SpeedOfLight;
      const double segmentlength = (postPosition-prePosition).GetLength();
      double deltaT = segmentlength/(SpeedOfLight*corrBetaValue)*fabs(midDoppler);
        
      if (preResT < postResT)
      { // startE arrives earlier
        preResT = midResT-0.5*deltaT;
        postResT = midResT+0.5*deltaT;
      }
      else
      { // endE arrives earlier
        postResT = midResT-0.5*deltaT;
        preResT = midResT+0.5*deltaT;
      }

      const long double gridresolution = (*currGroundElement)->GetTimeResolution();
      deltaT = postResT-preResT; //have to recalculate, because deltaT can be negative
        
      // redistribute contributions over time scale defined by the observation time resolution
      if (fabs(deltaT) < gridresolution)
      {
        startE *= fabs(deltaT/gridresolution);
        endE *= fabs(deltaT/gridresolution);
          
        const long startBin = static_cast<long>(floor(preResT/gridresolution+0.5l));
        const long endBin = static_cast<long>(floor(postResT/gridresolution+0.5l));
        const double startBinFraction = (preResT/gridresolution)-floor(preResT/gridresolution);
        const double endBinFraction = (postResT/gridresolution)-floor(postResT/gridresolution);

        // only do timing modification if contributions would land in same bin
        if (startBin == endBin)
        {
          // if startE arrives before endE
          if (deltaT >= 0)
          { 
            if ((startBinFraction >= 0.5) && (endBinFraction >= 0.5)) // both points left of bin center
            {
              preResT -= gridresolution; // shift startE to previous gridpoint
            }
            else if ((startBinFraction < 0.5) && (endBinFraction < 0.5)) // both points right of bin center
            {
              postResT += gridresolution; // shift endE to next gridpoint
            }
            else // points on both sides of bin center
            {
              const double leftDist = 1.0-startBinFraction;
              const double rightDist = endBinFraction;
              // check if asymmetry to right or left
              if (rightDist >= leftDist)
              {
                postResT += gridresolution; // shift endE to next gridpoint
              }
              else
              {
                preResT -= gridresolution;  // shift startE to previous gridpoint
              }
            }          
          }
          else // if endE arrives before startE
          {
            if ((startBinFraction >= 0.5) && (endBinFraction >= 0.5)) // both points left of bin center
            {
              postResT -= gridresolution; // shift endE to previous gridpoint
            }
            else if ((startBinFraction < 0.5) && (endBinFraction < 0.5)) // both points right of bin center
            {
              preResT += gridresolution; // shift startE to next gridpoint
            }
            else // points on both sides of bin center
            {
              const double leftDist = 1.0-endBinFraction;
              const double rightDist = startBinFraction;
              // check if asymmetry to right or left
              if (rightDist >= leftDist)
              {
                preResT += gridresolution; // shift startE to next gridpoint
              }
              else
              {
                postResT -= gridresolution;  // shift startE to previous gridpoint
              }
            }
          }
        } 
      }
    }
    else
    {
      // refractive index is unity or not near Cherenkov angle

      // check if preDoppler has become zero in case of refractive index of unity because of numerical limitations
      if (preDoppler == 0)
      {
        // redo calculation with higher precision
        long double indexL = preIndex;
        long double betaX = corrBeta.GetX();
        long double betaY = corrBeta.GetY();
        long double betaZ = corrBeta.GetZ();
        long double nX = preN.GetX();
        long double nY = preN.GetY();
        long double nZ = preN.GetZ();
        long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
        preDoppler = doppler;
      }

      startE = pre.weight * particleCharge * UnitCharge * ( preN.CrossedWith( preN.CrossedWith(corrBeta))) / 
               (SpeedOfLight * preDistance * preDoppler);

      // check if postDoppler has become zero in case of refractive index of unity because of numerical limitations
      if (postDoppler == 0)
      {
        // redo calculation with higher precision
        long double indexL = postIndex;
        long double betaX = corrBeta.GetX();
        long double betaY = corrBeta.GetY();
        long double betaZ = corrBeta.GetZ();
        long double nX = postN.GetX();
        long double nY = postN.GetY();
        long double nZ = postN.GetZ();
        long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
        postDoppler = doppler;
      }
      
      endE = -1.0 * post.weight * particleCharge * UnitCharge * ( postN.CrossedWith( postN.CrossedWith(corrBeta))) / 
             (SpeedOfLight * postDistance * postDoppler);

      // if preDoppler or postDoppler are below a certain threshold, redistribute contributions over two consecutive bins (take into account finite detector time resolution)
      if ((preDoppler<1.e-9) || (postDoppler<1.e-9))
      {
        const long double gridresolution = (*currGroundElement)->GetTimeResolution();
        double deltaT = postResT-preResT;
        if (fabs(deltaT) < gridresolution)
        {
          startE *= fabs(deltaT/gridresolution);
          endE *= fabs(deltaT/gridresolution);
          
          const long startBin = static_cast<long>(floor(preResT/gridresolution+0.5l));
          const long endBin = static_cast<long>(floor(postResT/gridresolution+0.5l));
          const double startBinFraction = (preResT/gridresolution)-floor(preResT/gridresolution);
          const double endBinFraction = (postResT/gridresolution)-floor(postResT/gridresolution);

          // only do timing modification if contributions would land in same bin
          if (startBin == endBin)
          {
            if ((startBinFraction >= 0.5) && (endBinFraction >= 0.5)) // both points left of bin center
            {
              preResT -= gridresolution; // shift startE to previous gridpoint
            }
            else if ((startBinFraction < 0.5) && (endBinFraction < 0.5)) // both points right of bin center
            {
              postResT += gridresolution; // shift endE to next gridpoint
            }
            else // points on both sides of bin center
            {
              const double leftDist = 1.0-startBinFraction;
              const double rightDist = endBinFraction;
              // check if asymmetry to right or left
              if (rightDist >= leftDist)
              {
                postResT += gridresolution; // shift endE to next gridpoint
              }
              else
              {
                preResT -= gridresolution;  // shift startE to previous gridpoint
              }
            }
          }
        }          
      }
    }
  
    if ((*currGroundElement)->GetMode() == Pattern)
    {
      //get angles in ResponseTable interface coordinate system (CCW from east, CORSIKA: CCW from North)
      const double startTheta = atan2(sqrt(preN.GetX()*preN.GetX()+preN.GetY()*preN.GetY()), -preN.GetZ() );
      const double startPhi = atan2(preN.GetY(), preN.GetX()) + Pi/2.;
      const double endTheta = atan2(sqrt(postN.GetX()*postN.GetX()+postN.GetY()*postN.GetY()), -postN.GetZ() );
      const double endPhi = atan2(postN.GetY(), postN.GetX()) + Pi/2.;

      //get weight factors from antenna response pattern 
      //if no response pattern is used (*currGroundElement)->GetEffectiveAntennaArea() returns 1
      const double weightStartE = sqrt((*currGroundElement)->GetEffectiveAntennaArea(startTheta, startPhi));
      const double weightEndE = sqrt((*currGroundElement)->GetEffectiveAntennaArea(endTheta, endPhi));
      startE *= weightStartE;
      endE *= weightEndE;
    }

    (*currGroundElement)->CollectValuePair(EFieldDataPoint(preResT,velocityWeight*startE),EFieldDataPoint(postResT,velocityWeight*endE));
  }
}


//~ void 
//~ TPlotter::AddTrackZHS(const crs::CParticle &pre, 
                   //~ const crs::CParticle &post) 
//~ {  

  //~ /*
    //~ Check whether we are dealing with primary particle. If so,
    //~ skip the track of the primary particle, which is in a different 
    //~ reference system as the shower. Otherwise set first interaction.
  //~ */

  //~ if (fPrimaryTrack) {
    //~ crs::CParticleFortranPtr pptr;
    //~ fPrimaryTrack = !prminfo_(pptr);          // the flag is false if the primary info is not yet available because we are still following the primary
    //~ if (fPrimaryTrack) {
      //~ if (fVerbosityLevel >= 2) {
        //~ cout << " CoREAS" << fThreadName <<": Primary track, skipping  " << endl;
      //~ }    
      //~ return;
    //~ }
    //~ else {
      //~ SetFirstInteraction(pptr->x, pptr->y, pptr->z, pptr->depth*g/cm/cm, pptr->time);
      //~ if (fVerbosityLevel >= 2) {
        //~ cout << " CoREAS" << fThreadName <<": Primary track end-pos: x=" << fFirstInteractionX/m 
             //~ << " y=" << fFirstInteractionY/m 
             //~ << " z=" << fFirstInteractionZ/m << endl;
      //~ }
    //~ }
  //~ }

  //~ if (fVerbosityLevel >= 9) {
  //~ double dX = post.x-pre.x; 
  //~ double dY = post.y-pre.y; 
  //~ double dZ = post.z-pre.z;
  //~ double length = sqrt (dX*dX + dY*dY + dZ*dZ);
  //~ double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
  //~ double postGamma = post.energy/gParticleMass[abs((int)post.particleId)-1];
  //~ double preBetaValue = sqrt(1.-1./(preGamma*preGamma));
  //~ double postBetaValue = sqrt(1.-1./(postGamma*postGamma));  
  //~ cout << " CoREAS" << fThreadName <<":   PRE> "
   //~ << " id=" << (int)pre.particleId
   //~ << " E=" << pre.energy
   //~ << " w=" << pre.weight
   //~ << " x=" << pre.x/km << "km"
   //~ << " y=" << pre.y/km << "km"
   //~ << " z=" << pre.z/km << "km"
   //~ << " t=" << pre.time*s/ns << "ns"
   //~ << " X=" << pre.depth << "g/cm^2"
//~ /* << " Px=" << pre.Px << ""
   //~ << " Py=" << pre.Py << ""
   //~ << " Pz=" << pre.Pz << "" */
   //~ << " v= " <<  length/(post.time*s/ns-pre.time*s/ns)/crs::cSpeedOfLight 
   //~ << " prebeta = " << preBetaValue
   //~ << " postbeta = " << postBetaValue
   //~ << " ratio = " << length/(post.time*s/ns-pre.time*s/ns)/crs::cSpeedOfLight/preBetaValue
   //~ << endl;
  //~ }
  
  //~ if (fVerbosityLevel >= 10) {
//~ /*double dX = post.x-fFirstInteractionX; 
  //~ double dY = post.y-fFirstInteractionY; 
  //~ double dZ = post.z-fFirstInteractionZ; 
  //~ double length = sqrt (dX*dX + dY*dY + dZ*dZ);*/
  //~ cout << " CoREAS" << fThreadName <<":  POST> "
   //~ << " id=" << (int)post.particleId
   //~ << " E=" << post.energy
   //~ << " w=" << post.weight
   //~ << " x=" << post.x/km << "km"
   //~ << " y=" << post.y/km << "km"
   //~ << " z=" << post.z/km << "km"
   //~ << " t=" << post.time*s/ns << "ns"
   //~ << " X=" << post.depth << "g/cm^2"
//~ // << " v= " <<  length/(post.time*s/ns-fFirstInteractionTime*s/ns)/crs::cSpeedOfLight
   //~ << endl;
  //~ }
    
  //~ /* Skip the particle if pre and post times are identical */
  //~ if (pre.time-post.time == 0.0)
    //~ return;  
 
  //~ const int particleId = abs((int)pre.particleId);
  
  //~ // check if particles to be counted at all
  //~ if (!fParticles.count(particleId)) 
    //~ return;

  //~ // calculate direct radio
  //~ // find a better way to specify charge signs
  //~ double particleCharge = gParticleCharge[abs((int)pre.particleId)-1];
      
  //~ ThreeVector prePosition(pre.x, pre.y, pre.z);
  //~ ThreeVector postPosition(post.x, post.y, post.z);

//~ /*** Histogramming

  //~ double stepLength = (postPosition-prePosition).GetLength();
  //~ double stepHeight = pre.z;  // assume flat geometry
  //~ if ((stepLength < itsNumLengthBins*itsLengthBinsMeters)
     //~ && (stepLength >= 0)
     //~ && (stepHeight < itsNumHeightBins*itsHeightBinsMeters)
     //~ && (stepHeight >= 0))
  //~ {
    //~ // calculate bin - height bins are outer loop
    //~ long heightBin = stepHeight/itsHeightBinsMeters;
    //~ long stepBin = stepLength/itsLengthBinsMeters;
    //~ itsHeightHistograms.at(heightBin*itsNumLengthBins+stepBin)++;
  //~ }
//~ */

  //~ double preIndex = GetRefractiveIndexAtPosition(prePosition);
  //~ double postIndex = GetRefractiveIndexAtPosition(postPosition);
      
  //~ double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
  //~ double postGamma = post.energy/gParticleMass[abs((int)post.particleId)-1];
      
  //~ double preBetaValue = sqrt(1.-1./(preGamma*preGamma));
  //~ double postBetaValue = sqrt(1.-1./(postGamma*postGamma));
    
  //~ ThreeVector currDirection(post.x-pre.x, post.y-pre.y, post.z-pre.z);
  //~ currDirection*=1./currDirection.GetLength();
      
//~ //ThreeVector preBeta = preBetaValue*currDirection;     // neglects curvature along track from pre to post!
//~ //ThreeVector postBeta = postBetaValue*currDirection;   // neglects curvature along track from pre to post!

  //~ double corrBetaValue = (postPosition-prePosition).GetLength()/(SpeedOfLight*(post.time-pre.time));
  //~ ThreeVector corrBeta = corrBetaValue * currDirection;
  //~ const double velocityWeight = 0.5*(postBetaValue+preBetaValue)/corrBetaValue;

//~ /*cout << "--------------------\n";
  //~ cout << "pre.particleId: " << pre.particleId << "\n";
  //~ cout << "post.particleId: " << post.particleId << "\n";
  //~ cout << "particleCharge: " << particleCharge << "\n";
  //~ cout << "pre.energy: " << pre.energy << "\n";
  //~ cout << "post.energy: " << post.energy << "\n";
  //~ cout << "preGamma: " << preGamma << "\n";
  //~ cout << "postGamma: " << postGamma << "\n";
  //~ cout << "preBetaValue: " << preBetaValue << "\n";
  //~ cout << "postBetaValue: " << postBetaValue << "\n";
  //~ cout << "prePosition: " << prePosition.GetX() << ", " << prePosition.GetY() << ", " << prePosition.GetZ() << "\n";     
  //~ cout << "postPosition: " << postPosition.GetX() << ", " << postPosition.GetY() << ", " << postPosition.GetZ() << "\n";      
  //~ cout << "preBeta: " << preBeta.GetX() << ", " << preBeta.GetY() << ", " << preBeta.GetZ() << "\n";     
  //~ cout << "postBeta: " << postBeta.GetX() << ", " << postBeta.GetY() << ", " << postBeta.GetZ() << "\n";      
  //~ cout << "corrBeta: " << corrBeta.GetX() << ", " << corrBeta.GetY() << ", " << corrBeta.GetZ() << "\n";      
//~ */

  //~ ThreeVector startE;
  //~ ThreeVector endE;
  
  //~ const vector<GroundElement*>& groundElements = fGroundArea->GetActiveGroundElementPointers();
  //~ for (vector<GroundElement*>::const_iterator currGroundElement = groundElements.begin();
      //~ currGroundElement != groundElements.end(); ++currGroundElement)  
  //~ {
    //~ ThreeVector obsPosition = (*currGroundElement)->GetPosition();

    //~ // check if particle contribution should be skipped due to special mode of observer bin
    //~ if ((*currGroundElement)->GetMode() != Normal)
    //~ {
      //~ // check whether to apply cuts in particle gamma
      //~ if ((*currGroundElement)->GetMode() == Gamma)
      //~ {
        //~ double preGamma = pre.energy/gParticleMass[abs((int)pre.particleId)-1];
        //~ if ((preGamma < (*currGroundElement)->GetModeDataLow()) || (preGamma >= (*currGroundElement)->GetModeDataHigh())) {
          //~ continue;
        //~ }
      //~ }
      //~ // check whether to apply cuts in particle height asl
      //~ else if ((*currGroundElement)->GetMode() == Height)
      //~ {
        //~ double height = GetHeightAtPosition(ThreeVector(pre.x, pre.y, pre.z));
        //~ if ((height < (*currGroundElement)->GetModeDataLow())
           //~ || (height >= (*currGroundElement)->GetModeDataHigh())) {
          //~ continue;
        //~ }
      //~ }
      //~ // check whether to apply cuts in particle slant depth
      //~ else if ((*currGroundElement)->GetMode() == SlantDepth)
      //~ {
        //~ if ((pre.depth < (*currGroundElement)->GetModeDataLow())
           //~ || (pre.depth >= (*currGroundElement)->GetModeDataHigh())) {
          //~ continue;
        //~ }
      //~ }
      //~ // check whether to apply cuts in distance from shower plane through the core
      //~ else if ((*currGroundElement)->GetMode() == Distance)
      //~ {
        //~ double distance = fCORSIKA->GetDistanceOfPoint(fRadioCore.GetX(), fRadioCore.GetY(), fRadioCore.GetZ())-fCORSIKA->GetDistanceOfSlantDepth(pre.depth*g/cm/cm);
        //~ if ((distance < (*currGroundElement)->GetModeDataLow())
           //~ || (distance >= (*currGroundElement)->GetModeDataHigh())) {
          //~ continue;
        //~ }
      //~ }
    //~ }    
        
    //~ ThreeVector preN = obsPosition-prePosition;
    //~ double preDistance = preN.GetLength();
    //~ preN *= 1./preDistance;
    //~ double preT = pre.time-fCoreHitTime;
    //~ double preResT = preT+GetEffectiveRefractiveIndexBetween(obsPosition,prePosition)*preDistance/SpeedOfLight;
        
    //~ ThreeVector postN = obsPosition-postPosition;
    //~ double postDistance = postN.GetLength();
    //~ postN *= 1./postDistance;
    //~ double postT = post.time-fCoreHitTime;
    //~ double postResT = postT+GetEffectiveRefractiveIndexBetween(obsPosition,postPosition)*postDistance/SpeedOfLight;

    //~ double dTp = postT-preT;


    //~ // skip particle if both start- and end-contributions lie outside time-window to be recorded
    //~ if (((preResT < (*currGroundElement)->GetTimeLowerBoundary()) || (preResT >= (*currGroundElement)->GetTimeUpperBoundary())) &&
       //~ ((postResT < (*currGroundElement)->GetTimeLowerBoundary()) || (postResT >= (*currGroundElement)->GetTimeUpperBoundary())))
      //~ continue;
        
    //~ double preDoppler = (1.0 - preIndex*corrBeta.DottedWith(preN));
    //~ double postDoppler = (1.0 - postIndex*corrBeta.DottedWith(postN));

    //~ const long double gridresolution = (*currGroundElement)->GetTimeResolution();

//~ //****************************

    //~ ThreeVector midPosition = 0.5*(prePosition+postPosition);
    //~ ThreeVector U = obsPosition-midPosition;
    //~ double midDistance = U.GetLength();
    //~ ThreeVector midN = 1./midDistance*U;
    //~ double midIndex = GetRefractiveIndexAtPosition(midPosition);
    //~ double midDoppler = (1.0 - midIndex*corrBeta.DottedWith(midN));

    //~ // check if midDoppler has become zero because of numerical limitations
    //~ if (midDoppler == 0)
    //~ {
      //~ // redo calculation with higher precision
      //~ long double indexL = midIndex;
      //~ long double betaX = corrBeta.GetX();
      //~ long double betaY = corrBeta.GetY();
      //~ long double betaZ = corrBeta.GetZ();
      //~ long double nX = midN.GetX();
      //~ long double nY = midN.GetY();
      //~ long double nZ = midN.GetZ();
      //~ long double doppler = 1.0l - indexL * (betaX*nX+betaY*nY+betaZ*nZ);
      //~ midDoppler = doppler;
    //~ }


    //~ double fDtDTpLimit=1.e-15;

    //~ //calculate vector potentials

    //~ ThreeVector V=SpeedOfLight*corrBeta;
    //~ double VdotU= V.DottedWith(U);
    //~ double constA=particleCharge*UnitCharge/midDistance; //mu*q/4/pi/R in volts * meter^-2 * seconds^2   -> CGS in SI // Now the vector potential is in V*s/m.
    //~ ThreeVector A=constA*(V-VdotU*U);

    //~ //two ways to calculate dT/dTp (Doppler factor), need to test which is best
    //~ double deltaT = postResT-preResT;
    //~ //double dT_dTp=fabs(deltaT)/dTp; //brute ratio way
    //~ //G4double dT_dTp=fabs(1.-(n/c)*(VdotU));  //formal derivative way, usually used
    //~ double dT_dTp=midDoppler;
    
    //~ double factor=0.;

    //~ long startBin = static_cast<long>(floor(preResT/gridresolution+0.5l));
    //~ long endBin = static_cast<long>(floor(postResT/gridresolution+0.5l));
    //~ const double startBinFraction = (preResT/gridresolution)-floor(preResT/gridresolution);
    //~ const double endBinFraction = (postResT/gridresolution)-floor(postResT/gridresolution);


    //~ // The binning and weighting to each bin depends on T1,T2, and dT/dTp 
    //~ //B. THE TRACK LANDS WITHIN ONE TIME BIN
    //~ if (startBin == endBin) {          
      //~ //   a. dt_Ant is not too close to zero
      //~ if (dT_dTp > fDtDTpLimit) {
        //~ //(T2-T1)/DT_BIN corrects for the ammount of time actually projected into the bin
        //~ factor=(fabs(deltaT)/gridresolution)/dT_dTp;
      //~ }
      //~ //  b. dt_Ant is too close to zero and the Cherenkov angle approximation is needed
      //~ else if (dT_dTp <= fDtDTpLimit) { //cherenkov limit, apply nascent delta function
        //~ factor=(dTp/gridresolution); // dTp is always positive.
      //~ }
      //~ (*currGroundElement)->CollectValuePair(EFieldDataPoint(preResT,factor*gridresolution*A),EFieldDataPoint(postResT,0.0*startE));
    //~ }
    //~ else  
    //~ //C. THE TRACK MAPS TO MULTIPLE BINS
    //~ {
      //~ // check correct ordering of bins
      //~ if (startBin>endBin) {
        //~ long tmp = startBin;
        //~ startBin = endBin;
        //~ endBin = tmp;
      //~ }
      //~ for (long iBin=startBin; iBin<=endBin; iBin++){
        //~ factor=1./dT_dTp;
        //~ //  a. Bin edges need correction factors so as not to over count
        //~ if (preResT <= postResT) {
          //~ if (iBin == startBin)
            //~ factor=(double(startBin+1)*gridresolution-preResT)/gridresolution/dT_dTp;
          //~ else if (iBin == endBin)
            //~ factor=((postResT-double(endBin)*gridresolution)/gridresolution)/dT_dTp;
        //~ }
        //~ else
        //~ {
          //~ if (iBin == startBin)
            //~ factor=(double(startBin+1)*gridresolution-postResT)/gridresolution/dT_dTp;
          //~ else if (iBin == endBin)
            //~ factor=((preResT-double(endBin)*gridresolution)/gridresolution)/dT_dTp;
        //~ }
      //~ (*currGroundElement)->CollectValuePair(EFieldDataPoint(preResT+gridresolution*(iBin-startBin),factor*gridresolution*A),EFieldDataPoint(postResT,0.0*startE));
      //~ }
    //~ }
  

//~ //****************************

  //~ }
//~ }




void 
TPlotter::AddInteraction(const crs::CInteraction& interaction) 
{
  if (fVerbosityLevel >= 10) {
    cout << " CoREAS" << fThreadName <<": ";
    interaction.Dump();
  }
}


void 
TPlotter::SetFirstInteraction(const double x, const double y, const double z, 
            const double X, const double t) 
{
  fFirstInteractionX = x;
  fFirstInteractionY = y;
  fFirstInteractionZ = z;
  fFirstInteractionTime = t;
  fFirstInteractionSlantDepth = X;
  if (X<0) {
    if (IsMainThread())
      cout << " CoREAS" << fThreadName <<": WARNING FirstInteractionDepth smaller than 0 (" << X/g*cm*cm << " g/cm2) " << endl;
    fFirstInteractionSlantDepth = 0;
  }
  //~ fFirstInteractionDist = fCORSIKA->GetDistanceOfPoint(x, y, z);
  
  // now that first interaction point is known, set up radio simulation, but do not re-initialize if already initialized (in particular in case of parallel version)
  if (!fGroundArea)
  {
    if (IsMainThread()) {
      cout << " CoREAS" << fThreadName // << ": SetFirstInteraction Dist=" << fFirstInteractionDist/km << "km "
           << ": x=" << fFirstInteractionX/km << "km "
           << ", y=" << fFirstInteractionY/km << "km "
           << ", z=" << fFirstInteractionZ/km << "km "
           << ", SlantDepth=" << fFirstInteractionSlantDepth/g*cm*cm << "g/cm2 "
           << ", Time=" << fFirstInteractionTime << "s "
           << "\n\n";
    }
    InitializeRadioSimulation();                                   
  }
}


void 
TPlotter::SetShowerAxis(const double zenith, const double azimuth) 
{  
  if (IsMainThread()) {
    cout << " \n\n CoREAS" << fThreadName <<": SetShowerAxis zenith=" << zenith/deg << "deg"  
        << "   azimuth=" << azimuth/deg << "deg"  
        << endl;
  }
  fZenith = zenith;
  fAzimuth = azimuth;
  fCosZenith = cos(zenith);
  fSinZenith = sin(zenith);
  fCosAzimuth = cos(azimuth);
  fSinAzimuth = sin(azimuth);
  //~ fAxisZ = fCosZenith;  //TH: I think this should be -fCosZenith (downward-pointing)
  //~ fAxisX = fSinZenith * fCosAzimuth;;
  //~ fAxisY = fSinZenith * fSinAzimuth;
}

double TPlotter::GetRefractiveIndexAtPosition(const ThreeVector& p1) const
{
  if (fSeaLevelRefractivity == 0.0)
    return 1.0;

  if (fCurved)
    return 1. + InterpolateInPlanarAtmosphere(GetHeightAtPosition(p1), fRefractivityTable); // this is correct in any geometry
  else
    return 1 + InterpolateInPlanarAtmosphere(p1.GetZ(), fRefractivityTable);
  
}

double TPlotter::GetEffectiveRefractiveIndexBetween(GroundElement* const observer, const ThreeVector& p2) const
{
  if (fSeaLevelRefractivity == 0.0)
    return 1.0; // more efficient than averaging (although same result)

  if ((!fCurved) || (fZenith < fTabulationThresholdZenithAngle))
  {
    // calculate in planar atmosphere
    double h1 = observer->GetZ();
    double h2 = p2.GetZ();

    double r1 = InterpolateInPlanarAtmosphere(h1, fIntegratedRefractivityTable);
    double r2 = InterpolateInPlanarAtmosphere(h2, fIntegratedRefractivityTable);

    double res = 1. + (r2 - r1) / (h2 - h1); // integrated refractivity is tabulated with growing values to greater heights
    return res;
  }
  else
  {
    // calculate in curved atmosphere, if possible using tabulated refractivities

    // calculate local zenith angle
    const double p1Height = observer->GetHeight();
    const ThreeVector p1 = observer->GetPosition();
    const ThreeVector p1p2 = p2-p1;
    const double trackDistance = p1p2.GetLength();
    const ThreeVector travelDirection = 1./trackDistance*p1p2; // some of this stuff could be passed to the function as parameter!

    if (p1Height < 0)
      return GetNumericallyIntegratedRefractiveIndexBetween(p1, trackDistance, travelDirection); // for heights below sea level, revert to numerical integration

    const ThreeVector radialDirection = (p1-ThreeVector(0,0,-cRearth)).GetDirection();
    const double cosTheta = travelDirection.DottedWith(radialDirection);

    if (cosTheta <= 0)
      return GetNumericallyIntegratedRefractiveIndexBetween(p1, trackDistance, travelDirection); // for exactly 90° or upgoing geometry revert to numerical integration

    const double tanTheta = sqrt(1.-cosTheta*cosTheta)/cosTheta; // see https://www.rapidtables.com/math/trigonometry/arccos/tan-of-arccos.html

    // look up zenith angle at Earth intersection from local zenith angle and antenna height
    const unsigned int tanThetaBinRaw = static_cast<unsigned int>(tanTheta/fTanThetaBinSpacing+0.5); // nearest neighbour interpolation

    if (tanThetaBinRaw >= fNumTanThetaBins)
      return GetNumericallyIntegratedRefractiveIndexBetween(p1, trackDistance, travelDirection); // for non-tabulated zenith angles, approximately >88.5°, revert to numerical integration

    // get Earth intersection angle from observer
    const unsigned int tanThetaBinAtEarth = observer->GetEarthIntersectionAngleBins().at(tanThetaBinRaw);
    
    if (tanThetaBinAtEarth >= fNumTanThetaBins)
      return GetNumericallyIntegratedRefractiveIndexBetween(p1, trackDistance, travelDirection); // for non-tabulated zenith angle at Earth intersection, approximately >88.5°, revert to numerical integration

    // calculate from integrated refractivity tables for curved geometries
    
    double d1 = observer->GetEarthIntersectionDistances().at(tanThetaBinRaw);
    double d2 = d1+trackDistance;
    double r1 = InterpolateInCurvedAtmosphere(tanThetaBinAtEarth, d1, fIntegratedRefractivityInCurvedAtmosphereTables);
    double r2 = InterpolateInCurvedAtmosphere(tanThetaBinAtEarth, d2, fIntegratedRefractivityInCurvedAtmosphereTables);
    double res = 1. + (r2 - r1) / (d2 - d1); // integrated refractivity is tabulated with growing values to greater distances

    return res;
  }
}

double TPlotter::GetNumericallyIntegratedRefractiveIndexBetween(const ThreeVector& p1, const double trackDistance, const ThreeVector& travelDirection) const
{
  //~ std::cout << ".";
  // numerically integrate the density of the atmosphere along the line of sight between p1 and p2 to calculate averaged refractive index
  int numSteps = static_cast<int>(trackDistance/500000.)+1; // one step every 5 km, relevant for zenith angles >85 degrees
  if (numSteps < 25)
    numSteps = 25; // do at least 25 steps, relevant for zenith angles <85 degrees
  const double stepSize = trackDistance/numSteps;
  double traversedN = 0.0;
  for (int jj = 0; jj<numSteps; ++jj)
  {
    double localHeight = GetHeightAtPosition(p1+stepSize*(jj+0.5)*travelDirection);
    traversedN += InterpolateInPlanarAtmosphere(localHeight, fRefractivityTable);
  }
  return 1. + traversedN / numSteps;
}


double TPlotter::GetHeightAtPosition(const ThreeVector& p1) const
{
  //warning: unify TPlotter::GetHeightAtPosition() and TCorsika::GetHeightOfPoint(), the latter of which seems to be wrong
  if (fCurved)
    return sqrt(p1.GetX()*p1.GetX() + p1.GetY()*p1.GetY() + p1.GetZ()*p1.GetZ() + 2*cRearth*p1.GetZ() + cRearth*cRearth) - cRearth;
  else
    return p1.GetZ();
}

// This function is called from within THRadioInterface to provide the GDAS table that has been read in by core CORSIKA to the radio part
void TPlotter::FillRefractivityTablesFromGDAS(int nPoints, const double* height, const double* refractiveIndex)
{
  fHeightTable.reserve(nPoints);
  fRefractivityTable.reserve(nPoints);
  for (int i=0; i<nPoints; i++)
  {
    fHeightTable.push_back(height[i] * crs::meter);
    fRefractivityTable.push_back(refractiveIndex[i] - 1.);
  }
}

/// Assumes equally spaced height values!
double TPlotter::InterpolateInPlanarAtmosphere(double h, const std::vector<double>& vec) const
{
  double fidx = (h - fHeightTable.front()) / fDistanceIncrement;
  size_t idx = size_t(fidx);
  double f = fidx - idx;

  // outside lower edge of table? then interpolate linearly from evolution in first 10 bins
  if (fidx <= 0) {
    const double slope10 = (vec[10]-vec[0]) / 10.;
    return vec[0]+slope10*fidx;
  }

  const unsigned long lastElement = vec.size()-1;
  if (fidx >= lastElement) {
    const double slope10 = (vec[lastElement]-vec[lastElement-10]) / 10.;
    return vec[lastElement]+slope10*(fidx-lastElement);
  }

  return  (1-f) * vec[idx] + f * vec[idx+1];
}

/// Assumes equally spaced distance values within each given table!
double TPlotter::InterpolateInCurvedAtmosphere(int tanThetaBin, double distance, const std::vector<std::vector<float> >& table) const
{
  // not nice that this function diretcly uses member variables not passed as arguments
  
  double minDistance = fDistanceTableMinDistances.at(tanThetaBin);
  double delta = fDistanceTableIncrements.at(tanThetaBin);

  double fidx = (distance - minDistance) / delta;
  size_t idx = size_t(fidx);
  double f = fidx - idx;

  // outside lower edge of table? then interpolate linearly from evolution in first 10 bins
  if (fidx <= 0) {
    const double slope10 = (table.at(tanThetaBin).at(10)-table.at(tanThetaBin).at(0)) / 10.;
    return table.at(tanThetaBin).at(0)+slope10*fidx;
  }

  // outside upper edge of table? then interpolate linearly from evolution in last 10 bins
  const unsigned long lastElement = table.at(tanThetaBin).size()-1;
  if (fidx >= lastElement) {
    const double slope10 = (table.at(tanThetaBin).at(lastElement)-table.at(tanThetaBin).at(lastElement-10)) / 10.;
    return table.at(tanThetaBin).at(lastElement)+slope10*(fidx-lastElement);
  }

  return  (1-f) * table.at(tanThetaBin)[idx] + f * table.at(tanThetaBin)[idx+1];
}

int TPlotter::GetIntersectionsLineCircle(double r, double b, double c, double& ax, double& ay, double& bx, double& by)
{
  // solution from: https://cp-algorithms.com/geometry/circle-line-intersection.html
  const double a = 1.; // without loss of generality set a to 1
  const double EPS = 1.e-6;
  double x0 = -a*c/(a*a+b*b);
  double y0 = -b*c/(a*a+b*b);
  if (c*c > r*r*(a*a+b*b)+EPS) {
    ax = x0;
    ay = y0;
    return 0;
  }
  else if (fabs(c*c - r*r*(a*a+b*b)) < EPS) {
    ax = x0;
    ay = y0;
    return 1;
  }
  else {
    double d = r*r - c*c/(a*a+b*b);
    double mult = sqrt (d / (a*a+b*b));
    ax = x0 + b * mult;
    bx = x0 - b * mult;
    ay = y0 - a * mult;
    by = y0 + a * mult;
    return 2;
  }
}
  

