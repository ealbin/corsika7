/***************************************************************************
                          grid.h  -  description
                             -------------------
    begin                : Tue Jun 3 2003
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

#ifndef GRID_H
#define GRID_H

#include <map>
#include <vector>
#include <iostream>
#include "datapoint.h"
#include "definitions.h"

//#include <fstream>

#include "threevector.h"

using std::ostream;
using std::endl;
using std::vector;
using std::cout;
using std::endl;

/**
  *@author Tim Huege
  */

//forward declaration necessary to let the following forward declarations know about the template arguments of Grid
template <typename K, typename T>
class Grid;

//forward declaration necessary to make operator<< a friend template function
template <typename K, typename T>
ostream& operator<< (ostream& os, const Grid<K,T>& rhs);


template <typename K, typename T>
class Grid {
public: 
//  Grid();
  Grid(const K& askedBinSize, const K& askedLowerBoundary, const K& askedUpperBoundary);
  ~Grid();
  void IncorporateValuePair(const DataPoint<K,T>& start, const DataPoint<K,T>& end);
  ostream& Print(ostream& os) const;  // member function to be called by operator<<
  long GetNumEntries() const { return itsNumEntries-TimeWindowSafetyCutBins; }
  void EnlargeNumberOfEntriesBy(unsigned long numSamples);

  void ReadSample(long num, K& key, T& value) const;

  
protected:
  K itsBinSize;
  K itsLowerBoundary;
  K itsUpperBoundary;
  vector<T>  itsGridVector;
  long itsLowerBoundaryIndex;
  long itsNumEntries;

  long GetGridIndex(K pointPosition) { return static_cast<long>(floor(pointPosition/this->itsBinSize+0.5l))-itsLowerBoundaryIndex; }
};

#endif
