/* $Id: MCorsikaReader.cc 5116 2016-01-04 19:09:04Z darko $   */

#include <crs2r/MCorsikaReader.h>
#include <crs2r/TC2R.h>
using namespace crs2r;

#include <crs/TSubBlock.h>

#include <string>
#include <iostream>
using namespace std;


MCorsikaReader::MCorsikaReader(const std::string &fname, int verbose)
: crsRead::TCorsikaReader(fname, verbose) {
}


/*
  This is called by TCorsikaReader 
*/
void MCorsikaReader::Init() {

  int iStart = 0;
  std::string::size_type iPos = GetInputFileName().rfind('/');
  if(iPos != std::string::npos)
    iStart = iPos+1;
    
  int iEnd = GetInputFileName().size();
    
  std::string OutputFileName =
    GetInputFileName().substr(iStart, iEnd-iStart);
  OutputFileName += ".root";

  fC2R = new TC2R();
  fC2R->Open(OutputFileName, GetThinning());
}	



/*
  This is called by TCorsikaReader 
*/
void MCorsikaReader::HandleSubBlock(const crs::TSubBlock &sb) {

  if(fVerboseLevel>10) {
    cout << " entering MCorsikaReader::HandleSubBlock " << sb.GetBlockTypeName() << endl;
  }

  fC2R->Write(sb);
}



/*
  This is called by TCorsikaReader 
*/
void MCorsikaReader::Exit() {

  fC2R->Close();
  delete fC2R;
}	
