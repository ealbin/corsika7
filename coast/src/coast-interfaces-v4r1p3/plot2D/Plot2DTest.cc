#include <crs/CorsikaTypes.h>
#include <crs/CorsikaConsts.h>
#include <crs/TSubBlock.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>

#include <crs2r/TC2R.h>

#include <TPlotter.h>

#include <sstream>
#include <iostream>
using namespace std;
using namespace crs;


///  Global object to handle all ROOT I/O operations and data conversions
TPlotter gPlotter;


int main() {

  gPlotter.Init();
  gPlotter.SetThinning(true);
  gPlotter.SetRun (0);
  gPlotter.SetFileName("test");

  gPlotter.SetEvent (0);
  
  const int weight = 1;
  gPlotter.DrawLine(0, 0, 0, weight, -.5*km, -4.1*km, 3*km, 5*km );
  gPlotter.Write();

  return 1;
}


extern "C" double heigh_(const double &f) {return 0;}




