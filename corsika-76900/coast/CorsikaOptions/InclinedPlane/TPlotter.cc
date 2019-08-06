/*
  RU, Di 9. Okt 18:14:43 CEST 2007
  
  - started new code for inclined particle sampling in CORSIKA

*/

#include <TPlotter.h>

#include <Plane.h>

#include <crs/CorsikaConsts.h>
#include <crs/CParticle.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MRunHeader.h>

//#include <crs2r/TC2R.h>

#include <crsRead/TSubBlockIO.h>
#include <crsRead/MParticleSubBlockOutput.h>
/*
#include <COASTconfig.h>

#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TSystem.h>
*/
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace crs;
//using namespace crs2r;
using namespace crsRead;
using namespace std;


const int TPlotter::fgNBlockSizeInfo = 4;  // magic fortran number


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

 */

TPlotter::TPlotter(OutputOption option) 
{  
  //fVerbosityLevel = 10;
  fVerbosityLevel = 0;
  
  fSlant = false;
  fCurved = false;
  fStackInput = false;
  
  fOutputOption = option;
  //  fC2R = 0;
  fOutputStream = 0;
  fSubBlockIO = 0;
  fParticleSubBlockIO = 0;
  
  Clear ();
  fDirectory = "";
}



TPlotter::~TPlotter () 
{
  //  delete fC2R;
  delete fOutputStream;
  delete fSubBlockIO;
  delete fParticleSubBlockIO;
}


void 
TPlotter::Welcome() 
	 const 
{
  cout << "\n"
       << "COAST: *************************************************\n"
       << "COAST: **                                               \n"
       << "COAST: ** You are using COAST inclined particle collecting  \n"
    //       << "COAST: **     COAST version: " << VERSION << " \n"
       << "COAST: ** Collector version: " << PLANE_VERSION << " \n"
       << "COAST: **                                               \n"
       << "COAST: **  THIN/CURVED/SLANT/STACKIN: " << fThinning 
       << "/" << fCurved << "/" << fSlant << "/" << fStackInput << "\n"
       << "COAST: **                                               \n"
       << "COAST: *************************************************\n\n"
       << endl;
}




void 
TPlotter::Clear() 
{  
  fAxisX = 0;
  fAxisY = 0;
  fAxisZ = 0;
  
  fCosZenith = 0;
  fSinZenith = 0;
  fCosAzimuth = 0;
  fSinAzimuth = 0;
  
  fEventNo = 0;
  fZenith = 0;
  fAzimuth = 0;
  fShowerEnergy = 0;
  fShowerPrimary = 0;
}


void
TPlotter::SetRunHeader(const crs::MRunHeader &header) 
{  
  ostringstream ssTmp; 
  ssTmp.str(""); ssTmp << header.GetVersion();  fCorsikaVersion = ssTmp.str();
  //  ssTmp.str(""); ssTmp << VERSION;              fCoastVersion = ssTmp.str();
  ssTmp.str(""); ssTmp << PLANE_VERSION;        fHistVersion = ssTmp.str();
  
  SetRun ((int)header.GetRunID ());
  
  //Define Obervation plane coordinates
  fXPlane=header.GetSamplingPlanePointX();
  fYPlane=header.GetSamplingPlanePointY();
  fZPlane=header.GetSamplingPlanePointZ();
  fXNorm=header.GetSamplingPlaneNormalX();
  fYNorm=header.GetSamplingPlaneNormalY();
  fZNorm=header.GetSamplingPlaneNormalZ();
}


void
TPlotter::SetShowerHeader(const crs::MEventHeader &header) 
{  
  if (header.GetSkimmingIncidence()) {
    fSkimming = true;
    fSkimmingAltitude = header.GetSkimmingAltitude() * cm;
    cout << " COAST: Detected skimming geometry with impact=" 
         << fSkimmingAltitude/m
         << " m" << endl;
  }
  
  // to flag the track of the primary particle
  
  fEventHeader = header;
  SetEvent(header.GetEventNumber());
  fShowerPrimary = (int)header.GetParticleId();

  fParticleSubBlockIO = new crsRead::MParticleSubBlockOutput(fThinning, fgNBlockSizeInfo, fEventHeader);
  
  SetShowerAxis(header.GetTheta()*rad, header.GetPhi()*rad);
  fShowerEnergy = header.GetEnergy()*GeV;

  if (!fCurved) {
    
    if ( header.GetZFirst() > 0. ) {   // if track start at first interaction point
      
      const double DXY = tan(header.GetTheta())*(header.GetZFirst() - header.GetObservationHeight(1));
      fXOff = cos(header.GetPhi()) * DXY;
      fYOff = sin(header.GetPhi()) * DXY;
      
    } else {  // if track start at first interaction point
      
      const double DXY = tan(header.GetTheta())*(header.GetStartingHeight() - header.GetObservationHeight(1));
      fXOff = cos(header.GetPhi()) * DXY;
      fYOff = sin(header.GetPhi()) * DXY;
    }
    
    fPlane = Plane(fXPlane+fXOff, fYPlane+fYOff, fZPlane, fXNorm, fYNorm, fZNorm);
    
  } else {  // No offset if CURVED option
    
    fPlane = Plane(fXPlane, fYPlane, fZPlane, fXNorm, fYNorm, fZNorm);
    
  }
}


void 
TPlotter::SetShowerTrailer (const crs::MEventEnd& /*trailer*/) 
{  
  // for binary output we have to write the remaining particle subblock
  if (fOutputOption==eCORSIKA && !fParticleSubBlockIO->IsFull()) {
    // write to disk
    fSubBlockIO->SetBuffer(fParticleSubBlockIO->GetData());
    (*fOutputStream) << (*fSubBlockIO);
    fParticleSubBlockIO->Clear();
  }
}


void 
TPlotter::Init() 
{  
  ostringstream filename;
  filename << fDirectory << "/" << fFileName << ".inclined";
  
  switch(fOutputOption) {
  case eNONE:
    cerr << " COAST: TPlotter::Init() No output option specified. EXIT " << endl;
    exit(1);
    break;
    /*  case eROOT:
    filename << ".root";
    cout << " COAST: Open file: \"" << filename.str() << "\" for the output of inclined plane sampling.\n" 
	 << "        Output-mode is set to ROOT."
	 << endl;
    fC2R = new crs2r::TC2R();  
    fC2R->Open(filename.str(), fThinning);
    break;*/
  case eCORSIKA:
    filename << "";
    cout << " COAST: Open file: \"" << filename.str() << "\" for the output of inclined plane sampling.\n"
	 << "        Output-mode is set to CORSIKA."
	 << endl;
    fOutputStream = new std::ofstream(filename.str().c_str(), 
				      std::ios::binary | std::ios::out | std::ios::trunc);
    fSubBlockIO = new crsRead::TSubBlockIO(fThinning, fgNBlockSizeInfo);
    fParticleSubBlockIO = 0; 
    break;
  }
}

void
TPlotter::Close() 
{  
  switch(fOutputOption) {
    /*  case eROOT:
    fC2R->Close();
    break;*/
  case eCORSIKA:
    // we have to fill up the remaining binary BLOCK with zeros, before closing the file
    while ( ((fSubBlockIO->GetSubBlockCounter())%crs::TBlock::fgNSubBlocks) ) {

      const int subBlockSize = crs::TBlock::fgNParticles *
	(fThinning ? 
	 crs::TBlock::fgNEntriesThinned :
	 crs::TBlock::fgNEntriesNotThinned);
      
      // and also init the TSubBlock
      CREAL* zeroBlock = new CREAL[subBlockSize];
      for (int i=0; i<subBlockSize; ++i) { zeroBlock[i] = 0; }
      fSubBlockIO->SetBuffer(zeroBlock);
      (*fOutputStream) << (*fSubBlockIO);
    }
    fOutputStream->flush();
    fOutputStream->close();
    break;
  case eNONE:
    break;
  }
}


void
TPlotter::Write (const crs::TSubBlock &SubBlock) 
{  
  switch(fOutputOption) {
    /*case eROOT:
    fC2R->Write(SubBlock);
    break;*/
  case eCORSIKA:
    fSubBlockIO->SetBuffer(SubBlock);
    (*fOutputStream) << (*fSubBlockIO);
    break;
  case eNONE:
    break;
  }
}

void
TPlotter::Write(const crs::GroundParticle& p) 
{
  switch(fOutputOption) {
    /*case eROOT:

    fC2R->Write(fEventHeader, p);
    break;
    */

  case eCORSIKA:

    (*fParticleSubBlockIO) << p;

    if (fParticleSubBlockIO->IsFull()) { // write to disk
      
      fSubBlockIO->SetBuffer(fParticleSubBlockIO->GetData());
      (*fOutputStream) << (*fSubBlockIO);
      fParticleSubBlockIO->Clear();
    }
    break;

  case eNONE:
    break;
  }
}


void 
TPlotter::Write () 
{
  if (fVerbosityLevel>=2) {
    cout << " COAST: Write " << endl;
  }
    
  // clean up
  Clear ();
}




// #define Rotate(x, y, sx, sy, cosa, sina) { 
inline
void 
TPlotter::Rotate (double x, double y, double z,
		  double &sx, double &sy, double &sz,
		  int inverse) 
{
  // rotate around z by azimuth
  double sx_ =             x*fCosAzimuth + inverse * y*fSinAzimuth;
  double sy_ = - inverse * x*fSinAzimuth +           y*fCosAzimuth; 
  double sz_ =   z;
  
  // rotate around y` by zenith
  double sx__ =             sx_*fCosZenith - inverse * sz_*fSinZenith;
  double sy__ = sy_;
  double sz__ = + inverse * sx_*fSinZenith +           sz_*fCosZenith; 

  // rotate back around z`` by -azimuth 
  sx =             sx__*fCosAzimuth - inverse * sy__*fSinAzimuth;
  sy = + inverse * sx__*fSinAzimuth +           sy__*fCosAzimuth; 
  sz =   sz__;

}


void 
TPlotter::AddTrack(const crs::CParticle &pre, 
		   const crs::CParticle &post) 
{  
  if (fVerbosityLevel>=9) {
    cout << " COAST:   PRE> "
	 << " id=" << (int)pre.particleId
	 << " E=" << pre.energy
	 << " w=" << pre.weight
	 << " x=" << pre.x << "cm"
	 << " y=" << pre.y << "cm"
	 << " z=" << pre.z << "cm"
	 << " t=" << pre.time*s/ns << "ns"
	 << " X=" << pre.depth << "g/cm^2"
	 << endl;
  }
  
  if (fVerbosityLevel>=10) {
    cout << " COAST:  POST> "
	 << " id=" << (int)post.particleId
	 << " E=" << post.energy
	 << " w=" << post.weight
	 << " x=" << post.x << "cm"
	 << " y=" << post.y << "cm"
	 << " z=" << post.z << "cm"
	 << " t=" << post.time*s/ns << "ns"
	 << " X=" << post.depth << "g/cm^2"
	 << endl;
  }
  
  
  
  /*   ------------------------
       direction of particle in shower frame
  */
  const double dX = post.x - pre.x; 
  const double dY = post.y - pre.y; 
  const double dZ = post.z - pre.z; 
  const double length = sqrt (dX*dX + dY*dY + dZ*dZ);
    
  if (length==0) {
    /*   ------------------------ 
	 This happens for particles that are NOT tracked by Corsika
	 e.g. pi0 -> decay immediatly
	 SKIP them
    */
    return;
  }	
  
  
  // check for crossing of plane

  crs::GroundParticle interpolatedParticle;
  if (fPlane.IsIntersecting(pre, post, interpolatedParticle)) {

    if (!fCurved) {
      interpolatedParticle.x -= fXOff + fXPlane;
      interpolatedParticle.y -= fYOff + fYPlane; 
      //cout << fXOff << ' ' << fYOff << ' ' << fXPlane << endl ; 
    }

    Write(interpolatedParticle);
  }
}







void 
TPlotter::SetShowerAxis(double zenith, double azimuth) 
{  
  fZenith = zenith;
  fAzimuth = azimuth;
  fCosZenith = cos(zenith);
  fSinZenith = sin(zenith);
  fCosAzimuth = cos(azimuth);
  fSinAzimuth = sin(azimuth);
  fAxisZ = fCosZenith;
  fAxisX = fSinZenith * fCosAzimuth;;
  fAxisZ = fSinZenith * fSinAzimuth;
}


