#include <crsRead/MCorsikaReader.h>

#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <TFile.h>
#include <TH1D.h>

#include <iostream>
#include <sstream>
using namespace std;



void PlotParticles (const crs::MParticleBlock &Data, TH1D &hist);


int main (int argc, char **argv) {

    if (argc<2) {
	cout << "  please specify Corsika file " << endl;
	return 0;
    }

    string fname (argv [1]);
    crsRead::MCorsikaReader cr (fname, 3);

    int ShowerCounter = 0;

    crs::MRunHeader Run;
    while (cr.GetRun (Run)) {

	crs::MEventHeader Shower;
	while (cr.GetShower (Shower)) {

	    ShowerCounter ++;

	    ostringstream oFileName;
	    oFileName << "Shower"
		      << Shower.GetEventNumber () << ".root";
	    TFile oFile (oFileName.str ().c_str (), "RECREATE");
	    TH1D hEnergy ("energy", "energy spectra of EM", 50, -3, 1);

	    crs::TSubBlock Data;
	    while (cr.GetData (Data)) {

		switch (Data.GetBlockType ()) {

		    case crs::TSubBlock::ePARTDATA:
			PlotParticles (Data, hEnergy);
			break;

		    case crs::TSubBlock::eLONG:
			break;

		    default:
			break;
		}


	    } // loop data

	    crs::MEventEnd ShowerSummary;
	    cr.GetShowerSummary (ShowerSummary);

	    oFile.cd ();
	    hEnergy.Write ();
	    oFile.Close ();

	} // loop shower

    } // loop runs (usually just 1)

    cout << " Read " << ShowerCounter << " showers from file " << endl;

    return 1;
}



void PlotParticles (const crs::MParticleBlock &Data, TH1D &hist) {

    crs::MParticleBlock::ParticleListConstIterator iEntry;
    for (iEntry = Data.ParticlesEnd();
	 iEntry != Data.ParticlesBegin();
	 ++iEntry) {

	if (iEntry->IsParticle ()) {

	    crs::MParticle iPart (*iEntry);

	    if (iPart.GetParticleID ()<4) {

		hist.Fill (log10 (iPart.GetKinEnergy ()),
			   iPart.GetWeight ());

	    }

	}

    }

}
