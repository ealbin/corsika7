#include <crsRead/MCorsikaReader.h>

#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <iostream>
#include <sstream>
#include <map>
using namespace std;


// to hold data of one observation level
struct ObsLevel {
  ObsLevel() : nPart(0), x(0), y(0), x2(0), y2(0), w(0) {}
  int nPart;
  double x;
  double y;
  double x2;
  double y2;
  double w;
};




int
main (int argc, char **argv)
{
  if (argc<2) {
    cout << "  please specify Corsika binary input file " << endl;
    return 1;
  }

  string fname(argv[1]);
  crsRead::MCorsikaReader cr(fname, 3);

  string inputFile;
  if (fname.rfind('/')==string::npos) {
    inputFile = fname;
  } else {
    inputFile = fname.substr(fname.rfind('/')+1,
                             fname.size()-fname.rfind('/'));
  }

  cout << "Analysing input file: " << inputFile << endl;

  int ShowerCounter = 0;

  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {

    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {

      ++ShowerCounter;

      const int nObsLevel = Shower.GetNObservationLevels();

      cout << "Shower Event #" << Shower.GetEventNumber() << " with "
	   << nObsLevel << " observation levels"
	   << endl;

      map<int, ObsLevel> obsLevel;

      for (int iObsLevel=1; iObsLevel<=nObsLevel; ++iObsLevel) {

        const double height = Shower.GetObservationHeight(iObsLevel-1);
        cout << " found observation-level " << iObsLevel << " at h=" << height << endl;
        ObsLevel emptyLevel;
        obsLevel[iObsLevel] = emptyLevel;

      } // end loop observation levels


      crs::TSubBlock Data;
      while (cr.GetData (Data)) {

        switch (Data.GetBlockType ()) {

            case crs::TSubBlock::ePARTDATA:
            {
              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.ParticlesBegin();
                   iEntry != ParticleData.ParticlesEnd();
                   ++iEntry) {

                if (iEntry->IsParticle()) {

                  crs::MParticle iPart(*iEntry);

                  const int level = iPart.GetObservationLevel();
                  const double w  = iPart.GetWeight();
                  const double x  = iPart.GetX();
                  const double y  = iPart.GetY();

		  if (obsLevel.count(level)==0) {
		    cout << " detected new obs-level " << level
			 << ". Possibly INCLIN-Option " << endl;
		    ObsLevel emptyLevel;
		    obsLevel[level] = emptyLevel;
		  }

                  obsLevel[level].nPart++;
                  obsLevel[level].x  += x*w;
                  obsLevel[level].y  += y*w;
                  obsLevel[level].x2 += x*x*w;
                  obsLevel[level].y2 += y*y*w;
                  obsLevel[level].w  += w;

                }

              } // end particle loop

              break;
            }

            case crs::TSubBlock::eLONG:
              break;

            default:
              break;

        } // end data block

      } // loop data

      crs::MEventEnd ShowerSummary;
      cr.GetShowerSummary(ShowerSummary);
      const double Xmax = ShowerSummary.GetXmax();

      cout << "---------------------------------\n"
           << " Shower info:\n"
           << "  Xmax = " << Xmax << "\n";

      for (map<int, ObsLevel>::iterator iLevel = obsLevel.begin();
	   iLevel != obsLevel.end();
	   ++iLevel) {

	cout << "----------------------------------------------\n";
	cout << "   observation level: " << iLevel->first
	     << " with " << iLevel->second.nPart << " particles "
	     << "\n";


	/*
	  cout << "    particle id: " << iPId->first
	  << " <x>=" << iPId->second.x/iPId->second.w
	  << " <y>=" << iPId->second.y/iPId->second.w
	  << " <x_rms>=" << sqrt(iPId->second.x2/iPId->second.w-pow(iPId->second.x/iPId->second.w,2))
	  << " <y_rms>=" << sqrt(iPId->second.y2/iPId->second.w-pow(iPId->second.y/iPId->second.w,2))
	  << "\n";
	*/

      } // loop observation levels

    } // loop shower

  } // loop runs (usually just 1)

  cout << " Read " << ShowerCounter << " showers from file " << endl;

  return 0;
}



