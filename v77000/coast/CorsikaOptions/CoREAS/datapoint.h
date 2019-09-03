/***************************************************************************
                          datapoint.h  -  description
                             -------------------
    begin                : Wed Jun 4 2003
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

#ifndef DATAPOINT_H
#define DATAPOINT_H

#include <iostream>

using std::ostream;
using std::istream;
using std::endl;

//for issues with friends and templates see http://gcc.gnu.org/fom_serv/cache/32.html

//forward declaration necessary to let the following forward declarations know about the template arguments of DataPoint
template <typename K, typename T>
class DataPoint;

//forward declaration necessary to make operator<< a friend template function
template <typename K, typename T>
ostream& operator<< (ostream& os, const DataPoint<K,T>& rhs);

//forward declaration necessary to make operator>> a friend template function
template <typename K, typename T>
istream& operator>> (istream& is, DataPoint<K,T>& rhs);

//forward declaration necessary to make SimpleGrid<K,T> a friend template class
template <typename K, typename T>
class Grid;


/*leave out any <K,T> for parameters and return values in the original class
 *declaration,but include them in all definitions and declarations outside
 *original class declaration (although not necessary everywhere)
 */

/**
  *A datapoint consisting of an absolute tag and associated data.
  *K is the key (e.g., a time stamp) and T the associated data (e.g., a vectorial quantity).
  *@author Tim Huege
  */
  
template <typename K, typename T>
class DataPoint {
  friend class Grid<K,T>; // possible due to forward declaration
public:
  DataPoint(); ///< standard constructor
  DataPoint(const K& parKey, const T& parData); ///< constructor with explicit data
  DataPoint(const DataPoint& rhs); ///< copy-constructor
  ~DataPoint();
  const K& GetKey() const { return itsKey; } ///< get key of data point
  const T& GetData() const { return itsData; } ///< get data of data point
  void SetKey(const double& newKey) { itsKey = newKey; } ///< set key of data point
  DataPoint& operator=(const DataPoint& rhs);
  bool operator<(const DataPoint& rhs) const; ///< compare keys of *this and rhs
  const DataPoint Interpolate(const K& parKey, const DataPoint& otherPoint) const; ///< interpolate data at key parKey from *this and otherPoint
  friend ostream& operator<< <> (ostream& os, const DataPoint<K,T>& rhs); // possible due to forward declaration
  friend istream& operator>> <> (istream& is, DataPoint<K,T>& rhs); // possible due to forward declaration
private:
  K itsKey;
  T itsData;
};

// constructors
template <typename K, typename T>
inline
DataPoint<K,T>::DataPoint()
: itsKey(0.0), itsData(T())
{}

template <typename K, typename T>
inline
DataPoint<K,T>::DataPoint(const K& parKey, const T& parData)
: itsKey(parKey), itsData(parData)
{}

template <typename K, typename T>
inline
DataPoint<K,T>::DataPoint(const DataPoint<K,T>& rhs)
: itsKey(rhs.itsKey), itsData(rhs.itsData)
{}


// defined here because it's inline
template <typename K, typename T>
inline
bool DataPoint<K,T>::operator<(const DataPoint<K,T>& rhs) const
{
  return (itsKey < rhs.itsKey);
}

template <typename K, typename T>
inline
const DataPoint<K,T> DataPoint<K,T>::Interpolate(const K& parKey, const DataPoint<K,T>& otherPoint) const
{
  return DataPoint<K,T>(parKey, itsData + ((parKey - itsKey)/(otherPoint.itsKey - itsKey)) * (otherPoint.itsData - itsData) );
}

#endif
