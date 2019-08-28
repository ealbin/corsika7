#ifndef _T_SIMPLECORSIKAREADER_H
#define _T_SIMPLECORSIKAREADER_H

#include <crsRead/TCorsikaReader.h>

#include <crs/TSubBlock.h>

#include <string>

class TCorsikaPlotter;


/** 
 \class TSimpleReader
 \brief simples possible CORSIKA binary file reader

 Does nothing, but stream in TSubBlocks and give them to a plotter
 class.

 \author Ralf Ulrich
 \date Wed Jun 15 17:00:08 CEST 2005
 \version $Id: TSimpleReader.h,v 1.1.1.1 2007-07-31 07:00:43 rulrich Exp $ 
*/

class TSimpleReader : public crsRead::TCorsikaReader {
  
public:
  TSimpleReader (const std::string &filename, int verbose=0);
  ~TSimpleReader ();
  
  
private:
  // overload the HandleSubBlock method to catch all the TSubBlocks
  virtual void HandleSubBlock (const crs::TSubBlock &CrsSubBlock);
  
  
private:
  TCorsikaPlotter *fPlotter;
  
};


#endif
