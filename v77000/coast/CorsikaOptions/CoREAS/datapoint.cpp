/***************************************************************************
                          datapoint.cpp  -  description
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

#include "datapoint.h"

template <typename K, typename T>
DataPoint<K,T>::~DataPoint()
{}

template <typename K, typename T>
DataPoint<K,T>& DataPoint<K,T>::operator=(const DataPoint<K,T>& rhs)
{
  itsKey = rhs.itsKey;
  itsData = rhs.itsData;
  return *this;
}

//is a friend
template <typename K, typename T>
ostream& operator<<(ostream& os, const DataPoint<K,T>& rhs)
{
  os << rhs.itsKey << "\t" << rhs.itsData;
  return os;
}

template <typename K, typename T>
istream& operator>>(istream& is, DataPoint<K,T>& rhs)
{
  is >> rhs.itsKey >> rhs.itsData;
  return is;
}


// necessary to avoid linker errors
// see http://www.fmi.uni-konstanz.de/~kuehl/c++-faq/containers-and-templates.html#faq-34.14
#include "threevector.h"
template class DataPoint<double,ThreeVector>;
template class DataPoint<double,double>;
template ostream& operator<< <double,ThreeVector> (ostream& os, const DataPoint<double,ThreeVector>& rhs);
template istream& operator>> <double,ThreeVector> (istream& is, DataPoint<double,ThreeVector>& rhs);
//template DataPoint<double,ThreeVector>::~DataPoint();
