#include <TSimpleReader.h>
#include <TCorsikaPlotter.h>

#include <crs/TSubBlock.h>

#include <string>
#include <iostream>
#include <fstream>
#include <memory>



TSimpleReader::TSimpleReader (const std::string &filename, 
			      int verbose)
: TCorsikaReader (filename, verbose) ,
  fPlotter (new TCorsikaPlotter ("Shower")) {
}


TSimpleReader::~TSimpleReader () {
  delete fPlotter;
}


void TSimpleReader::HandleSubBlock (const crs::TSubBlock &CrsSubBlock) {
  
  fPlotter->PlotSubBlock (CrsSubBlock);
}
