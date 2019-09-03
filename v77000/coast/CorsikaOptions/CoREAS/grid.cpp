/***************************************************************************
                          grid.cpp  -  description
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

#include "grid.h"


using std::make_pair;


template <typename K, typename T>
Grid<K,T>::Grid(const K& askedBinSize, const K& askedLowerBoundary, const K& askedUpperBoundary)
: itsBinSize(askedBinSize),itsLowerBoundary(askedLowerBoundary),itsUpperBoundary(askedUpperBoundary),
  itsLowerBoundaryIndex(static_cast<long>(floor(askedLowerBoundary/askedBinSize+0.5l))),
  itsNumEntries(static_cast<long>(floor((askedUpperBoundary-askedLowerBoundary)/askedBinSize+1.5l)))
{
  itsGridVector.resize(itsNumEntries);  // reserve space with empty entries
}

template <typename K, typename T>
Grid<K,T>::~Grid()
{
}

template <typename K, typename T>
void Grid<K,T>::IncorporateValuePair(const DataPoint<K,T>& start, const DataPoint<K,T>& end)
{
  const K startK = start.itsKey;
  if ((startK < itsUpperBoundary) && (startK >= itsLowerBoundary))
  {
    const long startKIndex = GetGridIndex(startK);
    itsGridVector[startKIndex]+=start.itsData;  // have to divide by binsize on printout!
  }

  const K endK = end.itsKey;
  if ((endK < itsUpperBoundary) && (endK >= itsLowerBoundary))
  {
    const long endKIndex = GetGridIndex(endK);
    itsGridVector[endKIndex]+=end.itsData;  // have to divide by binsize on printout!
  }
}

template <typename K, typename T>
void Grid<K,T>::EnlargeNumberOfEntriesBy(unsigned long numSamples)
{
  itsGridVector.resize(itsNumEntries+numSamples);
  itsUpperBoundary += itsBinSize*numSamples;
}


template <typename K, typename T>
ostream& Grid<K,T>::Print(ostream& os) const
{
  for (long i=0; i<(itsNumEntries-TimeWindowSafetyCutBins); ++i)
  {
    os << (itsLowerBoundaryIndex+i)*this->itsBinSize << "\t" << itsGridVector.at(i)/this->itsBinSize << "\n"; // have to divide by binsize on printout!
  }
  return os;
}

template <typename K, typename T>
void Grid<K,T>::ReadSample(long num, K& key, T& value) const
{ 
  key = (itsLowerBoundaryIndex+num)*this->itsBinSize;
  value = itsGridVector.at(num)/this->itsBinSize; // have to divide by binsize on printout!
}

template <typename K, typename T>
ostream& operator<<(ostream& os, const Grid<K,T>& rhs)
{
  rhs.Print(os);  // revert to implementation of grid
  return os;
}

// necessary to avoid linker errors
// see http://www.fmi.uni-konstanz.de/~kuehl/c++-faq/containers-and-templates.html#faq-34.14
#include "threevector.h"
template class Grid<double,ThreeVector>;
template ostream& operator<< <double,ThreeVector> (ostream& os, const Grid<double,ThreeVector>& rhs);
