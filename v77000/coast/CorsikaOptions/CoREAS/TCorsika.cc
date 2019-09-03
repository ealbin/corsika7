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

#include <TCorsika.h>
#include <crs/CorsikaConsts.h>

#include <iostream>
#include <iomanip>
#include <cstdlib>

using namespace crs;
using namespace std;


// constructor for flat+curved geometries
TCorsika::TCorsika(const double zenith, 
                   const double obsLevel,
                   const bool slant, 
                   const bool curved) :
  fVerbosityLevel(0),
  fCosZenith(cos(zenith)),
  fSinZenith(sin(zenith)),
  fZenith(zenith),
  fSkimmingAltitude(0),
  fObservationLevel(obsLevel),
  fHeightFirstInt(0.),
  fSlant(slant),
  fCurved(curved),
  fSkimming(false),
  heightTable(crlongi_.HLONG),
  depthTable(crslant_.THCKRL),
  distanceTable(crslant_.RLONG),
  TableLength(crlongi_.NSTEP+1)
  //~ fDecreasingHeightTable(heightTable[TableLength-1]<heightTable[0]) // set to true if height is decreasing with increasing slant depth
{
  if (curved && !slant) {
    cerr << " COAST: *************************** \n"
         << " COAST: ** \n"
         << " COAST: ** impossible to run \n"
         << " COAST: ** CURVED or SKIMMING geometries \n"
         << " COAST: ** without the SLANT option! \n"
         << " COAST: *************************** \n"
         << " COAST: aborting ... " << endl;
    exit(11);
  }

  // Check CRSLANT_ for existence
  //warning maybe remove CRSLANT_ test 
  if (curved) {
    const char* test = (char*) &crslant_;
    if (test[0]=='W' &&
        test[1]=='R' &&
        test[2]=='O' &&
        test[3]=='N' &&
        test[4]=='G') {
      cerr << " COAST: *************************** \n"
           << " COAST: ** \n"
           << " COAST: ** impossible to run \n"
           << " COAST: ** CURVED or SKIMMING geometries \n"
           << " COAST: ** without the SLANT option! \n"
           << " COAST: *************************** \n"
           << " COAST: aborting ... " << endl;
      exit(10);
    }
  }
}


// constructor for skimming geometries
TCorsika::TCorsika(const double skimmingHeight,
                   const double obsLevel) :
  fVerbosityLevel(0),
  fCosZenith(0),
  fSinZenith(0),
  fZenith(0),
  fSkimmingAltitude(skimmingHeight),
  fObservationLevel(obsLevel),
  fSlant(true),
  fCurved(true),
  fSkimming(true),
  heightTable(crlongi_.HLONG),
  depthTable(crslant_.THCKRL),
  distanceTable(crslant_.RLONG),
  TableLength(crlongi_.NSTEP+1)
  //~ fDecreasingHeightTable(heightTable[TableLength-1]<heightTable[0]) // set to true if height is decreasing with increasing slant depth
{
  // disallow skimming geometries for now
  cerr << " COAST: *************************** \n"
       << " COAST: ** \n"
       << " COAST: ** SKIMMING geometries are not supported yet\n"
       << " COAST: ** (implemented but untested)! \n"
       << " COAST: *************************** \n"
       << " COAST: aborting ... " << endl;
  exit(11);

  if (!fSlant) {
    cerr << " COAST: *************************** \n"
         << " COAST: ** \n"
         << " COAST: ** impossible to run \n"
         << " COAST: ** CURVED or SKIMMING geometries \n"
         << " COAST: ** without the SLANT option! \n"
         << " COAST: *************************** \n"
         << " COAST: aborting ... " << endl;
    exit(11);
  }
  
  // Check CRSLANT_ for existence
  //warning maybe remove CRSLANT_ test
  const char* test = (char*) &crslant_;
  if (test[0]=='W' &&
      test[1]=='R' &&
      test[2]=='O' &&
      test[3]=='N' &&
      test[4]=='G') {
    cerr << " COAST: *************************** \n"
	 << " COAST: ** \n"
	 << " COAST: ** impossible to run \n"
	 << " COAST: ** CURVED or SKIMMING geometries \n"
	 << " COAST: ** without the SLANT option! \n"
	 << " COAST: *************************** \n"
	 << " COAST: aborting ... " << endl;
    exit(10);
  }
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//                slant depth profle interpolation
//
// slantDepth, distance: starting at top of atmosphere, 
//                       increasing in direction of motion of the shower
// vertical height: vertical height above sea level
// vertical depth:  vertical depth in gcm^-2
// 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// universal irrespective of shower geometry
//~ double 
//~ TCorsika::GetVerticalDepthOfHeight(double height)
  //~ const 
//~ {
  //~ // ---------------------------------------------------
  //~ // pure vertical CORSIKA
  //~ return thick_(height/cm)*g/cm/cm;
//~ }


/*
  Get slant depth in a distance of the first interactionpoint
  in case of skimming incidence take tabulated atmosphere
*/
//~ double 
//~ TCorsika::GetSlantDepthOfDistance(double distance) 
  //~ const 
//~ {
  //~ // ---------------------------------------------------
  //~ if (not fCurved)
  //~ {
    //~ // pure vertical CORSIKA
    //~ double height = cTopOfAtmosphere - distance*fCosZenith;
    //~ double vertDepth = thick_(height/cm) * g/cm/cm;
    //~ return vertDepth / fCosZenith;
  //~ }

  //~ // take tabulated atmosphere
  //~ return GetTabulatedSlantDepthOfDistance(distance/cm) * g/cm/cm;
//~ }


//~ double 
//~ TCorsika::GetDistanceOfSlantDepth(double slantDepth) 
  //~ const 
//~ {
  //~ // ---------------------------------------------------
  //~ if (not fCurved)
  //~ {
    //~ // pure vertical CORSIKA
    //~ return (cTopOfAtmosphere - heigh_(slantDepth/g*cm*cm*fCosZenith)*cm) / fCosZenith; 
  //~ }

  //~ // take tabulated atmosphere
  //~ return GetTabulatedDistanceOfSlantDepth(slantDepth/g*cm*cm) * cm;
//~ }


//~ double 
//~ TCorsika::GetVerticalHeightOfDistance(double distance) 
  //~ const 
//~ {
  //~ // ---------------------------------------------------
  //~ if (not fCurved)
  //~ {
    //~ // pure vertical CORSIKA
    //~ return cTopOfAtmosphere - distance*fCosZenith;
  //~ }

  //~ // take tabulated atmosphere
  //~ return GetTabulatedHeightOfDistance(distance/cm) * cm;
//~ }


//~ void
//~ TCorsika::DumpAtmosphereTable()
  //~ const 
//~ {
  //~ // set up some plausibility tests to identify problems reading out common block
  //~ double stepSize = (depthTable[crlongi_.NSTEP]-depthTable[0])/(TableLength-1);
  //~ double previousDistance = -1.0;

  //~ cout << " COAST: Atmosphere definition table for curved atmosphere: \n";
  //~ for (int i=0; i<TableLength; ++i) {
    //~ cout << " COAST:" 
         //~ << " i=" << setw(4) << i
         //~ << " h=" << setw(9) << heightTable[i]/100000 << " km, "
         //~ << " d=" << setw(9) << distanceTable[i]/100000 << " km, "
         //~ << " x=" << setw(5) << depthTable[i] << " g/cm2"
         //~ << endl;

    //~ // check if RLONG-table is strictly monotonously rising
    //~ if (distanceTable[i] <= previousDistance)
    //~ {
      //~ cout << "Error: RLONG-table used in COAST is not strictly monotonous! "
           //~ << "Aborting ...\n";
      //~ exit(10);
    //~ }
    //~ else
    //~ {
      //~ previousDistance = distanceTable[i];
    //~ }

    //~ // check if THCKRL-table is equidistant
    //~ if (fabs(depthTable[i]-i*stepSize)>(1.0e-5*stepSize))
    //~ {
      //~ cout << "Error: THCKRL-table used in COAST is not equidistantly spaced! "
           //~ << "Aborting ...\n";
      //~ exit(10);
    //~ }
  //~ }
//~ }

/* deprecated and not generally well-define

double 
TCorsika::GetSlantDepthOfHeight(double height) 
  const 
{
  if (fCurved) { 
    return GetTabulatedSlantDepthOfHeight(height/cm)*g/cm/cm; 
  }
  return GetSlantDepthOfDistance((cTopOfAtmosphere-height)/fCosZenith);
}
*/


/*
double TCorsika::GetFirstLevelSlantDepth() const {
  if (fCurved) { 
    return GetTabulatedSlantDepthOfHeight(fHeightFirstInteraction/cm)*g/cm/cm; 
  }
  return GetSlantDepthOfDistance((cTopOfAtmosphere-fHeightFirstInteraction)/fCosZenith);
}*/


//~ double 
//~ TCorsika::GetLastLevelSlantDepth() 
  //~ const 
//~ {
  //~ if (fCurved) { 
    //~ return GetTabulatedSlantDepthOfHeight(fObservationLevel/cm)*g/cm/cm; 
  //~ }
  //~ return GetSlantDepthOfDistance((cTopOfAtmosphere-fObservationLevel)/fCosZenith);
//~ }

//~ double 
//~ TCorsika::GetHeightOfSlantDepth(double slantdepth) 
  //~ const
//~ {
  //~ if (fCurved) { 
    //~ return GetTabulatedHeightOfSlantDepth(slantdepth/g*cm*cm)*cm; 
  //~ }
  //~ return heigh_(slantdepth/g*cm*cm*fCosZenith)*cm;
//~ }

//warning TH: this function gives different results than GetHeightAtPosition(), it is not currently used in CoREAS
/*
double
TCorsika::GetHeightOfPoint(const double x, const double y, const double z) 
  const 
{
  if (fSkimming) {

    //warning check that coordinate origin for skimming showers is at Rp
    const double r = sqrt(pow(x, 2) + pow(y, 2));
    
    return sqrt(pow(cRearth+fSkimmingAltitude, 2) + pow(r, 2)) - cRearth;
  }
  if (fCurved) {
    //const double r = sqrt(pow(x, 2) + pow(y, 2));
    const double alpha = 90*deg - fZenith;    
    const double r = z / tan(alpha);
    const double S = pow(cRearth + z, 2) + pow(r, 2);
    if (S<0) {
      cout << "\n FATAL-ERROR: GetHeightOfZ has S=" << S << "<0 !!! \n" << endl;
      exit(10);
    }
    const double h = sqrt(S) - cRearth;
    if (fVerbosityLevel>60) {
      cout << " TCorsika::GetHeightOfZ: " << sqrt(S)/km << " " << cRearth/km << " h=" << h/km << " z=" << z/km << endl;
    }
    return h;
  }
  return z; // in the flat case z=height
}
*/

//~ double 
//~ TCorsika::GetDistanceOfPoint(const double x, const double y, const double z)
  //~ const
//~ {
  //~ if (fSkimming) {
    
    //~ // distance to margin of atmosphere
    //~ const double D = sqrt(pow(cRearth+cTopOfAtmosphere, 2) - pow(cRearth+fSkimmingAltitude, 2));

    //~ //warning check that coordinate origin for skimming showers is at Rp
    //~ const double r = sqrt(pow(x, 2) + pow(y, 2));
    
    //~ return D-r;
  //~ }
  //~ if (fCurved) {
    
    //~ // distance of point from origin
    //~ const double r_h = sqrt(pow(x, 2) + pow(y, 2));
    //~ const double r = sqrt(pow(r_h, 2) + pow(z-fObservationLevel, 2));
    
    //~ const double alpha = 180*deg - fZenith;
    //~ const double cosAlpha = cos(alpha);
    
    //~ // distance to margin of the atmosphere
    //~ const double D = (cRearth+fObservationLevel)*cosAlpha + sqrt(pow((cRearth+fObservationLevel)*cosAlpha, 2) - pow(cRearth+fObservationLevel, 2) + pow(cRearth+cTopOfAtmosphere, 2));

    //~ return D-r;
  //~ }
  //~ return (cTopOfAtmosphere - z) / fCosZenith;  
//~ }

//~ double 
//~ TCorsika::GetDistanceOfHeight(const double height) 
  //~ const
//~ {
  //~ if (fSkimming) {
    //~ cout << "\n\n"
	 //~ << "###############################################\n"
	 //~ << "# must not call GetDistanceOfHeight           #\n"
	 //~ << "# for skimming showers. It is not well defined#\n"
	 //~ << "###############################################\n\n"
	 //~ << endl;
    //~ exit (11);
  //~ }
 
  //~ double startingPoint ;
  //~ double finalSign ;
  //~ if (fZenith>90*deg) {  //upward going showers
    //~ startingPoint = fHeightFirstInt  ;
    //~ finalSign = -1. ;
  //~ } else {  //downward going showers
    //~ startingPoint = cTopOfAtmosphere ;
    //~ finalSign = 1 ;
  //~ }

  //~ if (fCurved) {
    //~ // return GetTabulatedDistanceOfHeight(height); can get really crude ... 
    //~ const double alpha = 180*deg - fZenith;
    //~ const double sinAlpha = sin(alpha);
    //~ const double beta_prime = asin(sinAlpha * (cRearth+fObservationLevel) / (cRearth + startingPoint));
    //~ const double beta = asin(sinAlpha * (cRearth+fObservationLevel) / (cRearth + height));
    //~ const double gamma_prime_prime = beta - beta_prime;
    //~ if (fVerbosityLevel>60) {
      //~ cout << " GetDistanceOfHeight: " << sin(gamma_prime_prime) * (cRearth + startingPoint) / sin(180*deg - beta) << endl;
    //~ }
    //~ return finalSign * sin(gamma_prime_prime) * (cRearth + startingPoint) / sin(180*deg - beta);
  //~ }  
  //~ return finalSign * (startingPoint - height) / fCosZenith;
//~ }


//~ long 
//~ TCorsika::GetBin(double Value, double* Table) 
  //~ const
//~ {
  //~ // outside table?
  //~ if ((Value < Table[0]) || (Value > Table[TableLength-1]))
  //~ {
    //~ std::cout << "\nError (GetBin): Out of range in AtmoTable!\n";
    //~ std::cout << "Asked for value " << Value << " when table covers values from " 
              //~ << Table[0] << " to " << Table[TableLength-1] << ".\n";
    //~ std::cout << "Aborting ...\n";
    //~ exit(10);
  //~ }

  //~ long lastLowerBin = 0;
  //~ long lastUpperBin = TableLength-1;
  //~ long currentBin = TableLength/2;
  //~ bool finished;
  //~ do	// interval search for the correct bin
  //~ {
    //~ finished = true;
    //~ if (Value < Table[currentBin])
    //~ {
      //~ lastUpperBin = currentBin;
      //~ currentBin = (lastLowerBin+lastUpperBin)/2;
      //~ finished = false;
    //~ }
    //~ else
    //~ {
      //~ if (Value > Table[currentBin+1])
      //~ {
        //~ lastLowerBin = currentBin;
        //~ currentBin = (lastLowerBin+lastUpperBin)/2;        // cannot (and must not) point further than to second last bin (AtmoTableLength-2)
        //~ finished = false;
      //~ }
    //~ }
  //~ } while (not finished);
  //~ return currentBin;  // return index of next-lower bin
//~ }


//~ long 
//~ TCorsika::GetBinReverse(double Value, double* Table) 
  //~ const
//~ {
  //~ // outside table?
  //~ if ((Value > Table[0]) || (Value < Table[TableLength-1]))
  //~ {
    //~ std::cout << "\nError (GetBinReverse): Out of range in AtmoTable!\n";
    //~ std::cout << "Asked for value " << Value << " when table covers values from " 
              //~ << Table[0] << " to " << Table[TableLength-1] << ".\n";
    //~ std::cout << "Aborting ...\n";
    //~ exit(10);
  //~ }

  //~ long lastLowerBin = 0;
  //~ long lastUpperBin = TableLength-1;
  //~ long currentBin = TableLength/2;
  //~ bool finished;
  //~ do	// interval search for the correct bin
  //~ {
    //~ finished = true;
    //~ if (Value > Table[currentBin])
    //~ {
      //~ lastUpperBin = currentBin;
      //~ currentBin = (lastLowerBin+lastUpperBin)/2;
      //~ finished = false;
    //~ }
    //~ else
    //~ {
      //~ if (Value < Table[currentBin+1])
      //~ {
        //~ lastLowerBin = currentBin;
        //~ currentBin = (lastLowerBin+lastUpperBin)/2;        // cannot (and must not) point further than to second last bin (AtmoTableLength-2)
        //~ finished = false;
      //~ }
    //~ }
  //~ } while  (not finished);
  //~ return currentBin;  // return index of next-lower bin
//~ }



//~ double 
//~ TCorsika::GetTabulatedSlantDepthOfDistance (double distance) 
  //~ const
//~ {
  //~ if (fVerbosityLevel > 50) {
    //~ cout << "TCorsika::GetTabulatedSlantDepthOfDistance d=" 
         //~ << distance/m << " m" << endl;
  //~ }
  //~ long currentBin = GetBin(distance, distanceTable);
  //~ double lowerD = distanceTable[currentBin];
  //~ double upperD = distanceTable[currentBin+1];
  //~ double lowerX = depthTable[currentBin];
  //~ double upperX = depthTable[currentBin+1];

  //~ if (currentBin==0) {
    //~ cout << "\nWARNING: (GetTabulatedSlantDepthOfDistance) Interpolation within first bin of atmosphere table is VERY INACCURATE !!!\n" << endl;
  //~ }

  //~ // linear interpolation
  //~ return (lowerX+(distance-lowerD)*(upperX-lowerX)/(upperD-lowerD));
//~ }

//~ double 
//~ TCorsika::GetTabulatedDistanceOfSlantDepth (double slantdepth) 
  //~ const
//~ {
  //~ if (fVerbosityLevel > 50) {
    //~ cout << "TCorsika::GetTabulatedDistanceOfSlantDepth x="
         //~ << slantdepth << " g/cm2" << endl;
  //~ }
  //~ long currentBin = GetBin(slantdepth, depthTable);
  //~ double lowerX = depthTable[currentBin];
  //~ double upperX = depthTable[currentBin+1];
  //~ double lowerD = distanceTable[currentBin];
  //~ double upperD = distanceTable[currentBin+1];

  //~ if (currentBin==0) {
    //~ cout << "\nWARNING: (GetTabulatedDistanceOfSlantDepth) Interpolation within first bin of atmosphere table is VERY INACCURATE !!!\n" << endl;
  //~ }

  //~ // linear interpolation
  //~ return (lowerD+(slantdepth-lowerX)*(upperD-lowerD)/(upperX-lowerX));
//~ }


//~ double 
//~ TCorsika::GetTabulatedHeightOfDistance (double distance) 
  //~ const
//~ {
  //~ if (fVerbosityLevel > 50) {
    //~ cout << "TCorsika::GetTabulatedHeightOfDistance d="
         //~ << distance/m << " m" << endl;
  //~ }
  //~ long currentBin = GetBin(distance, distanceTable);
  //~ double lowerD = distanceTable[currentBin];
  //~ double upperD = distanceTable[currentBin+1];
  //~ double lowerH = heightTable[currentBin];
  //~ double upperH = heightTable[currentBin+1];

  //~ if (currentBin==0) {
    //~ cout << "\nWARNING: (GetTabulatedHeightOfDistance) Interpolation within first bin of atmosphere table is VERY INACCURATE !!!\n" << endl;
  //~ }

  //~ // linear interpolation
  //~ return (lowerH+(distance-lowerD)*(upperH-lowerH)/(upperD-lowerD));
//~ }

//~ double 
//~ TCorsika::GetTabulatedHeightOfSlantDepth (double slantdepth) 
  //~ const
//~ {
  //~ if (fVerbosityLevel > 50) {
    //~ cout << "TCorsika::GetTabulatedHeightOfSlantDepth x="
         //~ << slantdepth/g*cm*cm << " g/cm^2" << endl;
  //~ }
  //~ long currentBin = GetBin(slantdepth, depthTable);
  //~ double lowerX = depthTable[currentBin];
  //~ double upperX = depthTable[currentBin+1];
  //~ double lowerH = heightTable[currentBin];
  //~ double upperH = heightTable[currentBin+1];

  //~ if (currentBin==0) {
    //~ cout << "\nWARNING: (GetTabulatedHeightOfSlantDepth) Interpolation within first bin of atmosphere table is VERY INACCURATE !!!\n" << endl;
  //~ }

  //~ // linear interpolation
  //~ return (lowerH+(slantdepth-lowerX)*(upperH-lowerH)/(upperX-lowerX));
//~ }

//~ // This function might need to be removed, since for skimming showers it is not well defined
//~ double 
//~ TCorsika::GetTabulatedSlantDepthOfHeight(double height) 
  //~ const
//~ {
  //~ if (fVerbosityLevel > 50) {
    //~ cout << "TCorsika::GetTabulatedSlantDepthOfHeight h="
         //~ << height/m << " m" << endl;
  //~ }
  //~ if (fSkimming) {
    //~ cout << "\n\n"
	 //~ << "##########################################################\n"
	 //~ << "# must not call TCorsika::GetTabulatedSlantDepthOfHeight #\n"
	 //~ << "# for skimming showers. It is not well defined           #\n"
	 //~ << "##########################################################\n\n"
	 //~ << endl;
    //~ exit (11);
  //~ }
  //~ long currentBin = 0;
  //~ if (fDecreasingHeightTable)
  //~ {
    //~ // "downward" going shower
    //~ currentBin = GetBinReverse(height, heightTable);
  //~ }
  //~ else
  //~ {
    //~ // "upward" going shower
    //~ currentBin = GetBin(height, heightTable);
  //~ }

  //~ if (currentBin==0) {
    //~ cout << "\nWARNING: (GetTabulatedSlantDepthOfHeight) Interpolation within first bin of atmosphere table is VERY INACCURATE !!!\n" << endl;
  //~ }

  //~ double lowerH = heightTable[currentBin];
  //~ double upperH = heightTable[currentBin+1];
  //~ double lowerX = depthTable[currentBin];
  //~ double upperX = depthTable[currentBin+1];
  //~ // linear interpolation
  //~ return (lowerX+(height-lowerH)*(upperX-lowerX)/(upperH-lowerH));
//~ }

//~ // This function might need to be removed, since for skimming showers it is not well defined
//~ double 
//~ TCorsika::GetTabulatedDistanceOfHeight(double height)
  //~ const
//~ {
  //~ if (fVerbosityLevel > 50) {
    //~ cout << "TCorsika::GetTabulatedDistanceOfHeight h="
         //~ << height/m << " m" << endl;
  //~ }
  //~ if (fSkimming) {
    //~ cout << "\n\n"
	 //~ << "##########################################################\n"
	 //~ << "# must not call TCorsika::GetTabulatedDistanceOfHeight   #\n"
	 //~ << "# for skimming showers. It is not well defined           #\n"
	 //~ << "##########################################################\n\n"
	 //~ << endl;
    //~ exit (11);
  //~ }
  //~ long currentBin = 0;
  //~ if (fDecreasingHeightTable)
    //~ {
    //~ // "downward" going shower
    //~ currentBin = GetBinReverse(height, heightTable);
  //~ }
  //~ else
  //~ {
    //~ // "upward" going shower
    //~ currentBin = GetBin(height, heightTable);
  //~ }

  //~ if (currentBin==0) {
    //~ cout << "\nWARNING: (GetTabulatedDistanceOfHeight) Interpolation within first bin of atmosphere table is VERY INACCURATE !!!\n" << endl;
  //~ }

  //~ double lowerH = heightTable[currentBin];
  //~ double upperH = heightTable[currentBin+1];
  //~ double lowerD = distanceTable[currentBin];
  //~ double upperD = distanceTable[currentBin+1];
  
  //~ // linear interpolation
  //~ return (lowerD+(height-lowerH)*(upperD-lowerD)/(upperH-lowerH));
//~ }

//~ double TCorsika::GetMaximumTabulatedDepth() const
//~ {
  //~ return depthTable[TableLength-1]*g/cm/cm;
//~ }
