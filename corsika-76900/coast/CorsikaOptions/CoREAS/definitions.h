 /***************************************************************************
                          definitions.h  -  description
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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <cmath>
#include <vector>
#include "threevector.h"
#include "datapoint.h"

using std::vector;


// version number
const double ProgramVersion = 1.3;

// constants
const double Pi = M_PI;
const double Pi_double = 2.0 * Pi;
const double Pi_half = Pi/2.0;
const double SpeedOfLight = 29979245800.0;					// in cm/s
const double UnitCharge =  4.803204196927238e-10;	// in cgs units
const long TimeWindowSafetyCutBins = 20;	// cut off so many bins at the end when writing out time series data

// typedefs
typedef double timestamp;
typedef ThreeVector EFieldVector;
typedef DataPoint<timestamp,EFieldVector> EFieldDataPoint;
typedef vector<EFieldDataPoint> EFieldTimeSeries;
typedef DataPoint<double,double> ValuePair;

//handy functions

/*
template<class T>
inline T min(T a, T b)
{
	if (a < b) return a;
	return b;
}

template<class T>
inline T max(T a, T b)
{
	if (a > b) return a;
	return b;
}
*/

template<class T>
inline void swapvars(T& a, T& b)
{
	T c = a;
	a = b;
	b = c;
}

template<class T>
inline bool fequal(T a, T b, double c)
{
	// special case: one or both parameters are exactly 0
	if (b == 0)
	{
		if (a == 0)
		{	return true; }
		else
		{ return false; }
	}
	// usual case
	return ((a/b-1 < c) && (a/b-1 > -c));
}


/*inline double min(double a, double b)
{
	if (a < b) return a;
	return b;
}

inline int min(int a, int b)
{
	if (a < b) return a;
	return b;
}

inline double max(double a, double b)
{
	if (a > b) return a;
	return b;
}

inline int max(int a, int b)
{
	if (a > b) return a;
	return b;
}


inline void swap(double& a, double& b)
{
	double c = a;
	a = b;
	b = c;
} */

inline double norm_phi(double phi)	// renormalizes an angle to be in the range [0, 2Pi)
{
	double tmpPhi = fmod(phi, Pi_double);
	if (tmpPhi < 0.0 ) { tmpPhi += Pi_double; }
	return tmpPhi;
}

inline double beta(double gamma)
{
	return sqrt(1.0 - 1.0 / (gamma * gamma));
}

inline double atan3(double x, double y)
{
	return atan2(y, x);
}

// pure helper class, not worth extra file
/// a complex number denoted by real and imaginary part
class Complex {
	public:
	Complex(double parReal, double parImag) : itsReal(parReal), itsImag(parImag) { }
	~Complex() { }
	double itsReal;
	double itsImag;
	double GetAmplitude() const { return sqrt(itsReal*itsReal+itsImag*itsImag); }
};

#endif
