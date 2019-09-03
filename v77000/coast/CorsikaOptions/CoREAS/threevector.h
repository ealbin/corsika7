/***************************************************************************
                          threevector.h  -  description
                             -------------------
    begin                : Tue May 27 2003
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

#ifndef THREEVECTOR_H
#define THREEVECTOR_H

#include <cmath>
#include <iostream>

using std::ostream;
using std::istream;

/// enum to specify rotation axes
enum Axis { XAxis, YAxis, ZAxis };
 
/**
  *A generic threedimensional vector class.
  *@author Tim Huege
  */

class ThreeVector {
public: 
  ThreeVector(); ///< constructor creating zero-vector
  ThreeVector(const double& px, const double& py, const double& pz); ///< constructor with given vector coordinates
  ~ThreeVector();

  const double& GetX() const;  ///< get cartesian x coordinate
  const double& GetY() const; ///< get cartesian y coordinate
  const double& GetZ() const; ///< get cartesian z coordinate

  ThreeVector& operator*= (const double& factor);
  ThreeVector& operator+= (const ThreeVector& rhs);
  ThreeVector& operator-= (const ThreeVector& rhs);
  const ThreeVector operator- (const ThreeVector& rhs) const;
  const ThreeVector operator+ (const ThreeVector& rhs) const;
  bool operator== (const ThreeVector& rhs) const;

  double GetLength() const; ///< get vector length
  double DottedWith(const ThreeVector& rhs) const; ///< calc dot product of *this with rhs
  const ThreeVector CrossedWith(const ThreeVector& rhs) const; ///< calc cross product of *this with rhs
  const ThreeVector GetDirection() const; ///< get unit vector pointing in direction of *this
  double GetAngleTo(const ThreeVector& rhs) const; ///< get angle of *this to rhs
  ThreeVector& Rotate(const double& parAngle, const Axis parAxis); ///< rotate *this by parAngle (radians) around parAxis
  
  // friends
  friend const ThreeVector operator* (const double& factor, const ThreeVector& rhs);
  friend const ThreeVector operator/ (const ThreeVector& rhs, const double& factor);
  friend ostream& operator<< (ostream& os, const ThreeVector& rhs);
  friend istream& operator>> (istream& is, ThreeVector& rhs);

  class Exception { };  // exception class
  
protected:
  double x;
  double y;
  double z;
};

// constructors
inline
ThreeVector::ThreeVector()
: x(0.0), y(0.0), z(0.0)
{}

inline
ThreeVector::ThreeVector(const double& px, const double& py, const double& pz)
: x(px), y(py), z(pz)
{}





// access functions

inline
const double& ThreeVector::GetX() const
{
  return x;
}

inline
const double& ThreeVector::GetY() const
{
  return y;
}

inline
const double& ThreeVector::GetZ() const
{
  return z;
}




// operators

inline
ThreeVector& ThreeVector::operator*= (const double& factor)
{
  x*=factor;
  y*=factor;
  z*=factor;
  return *this;
}

inline
ThreeVector& ThreeVector::operator+= (const ThreeVector& rhs)
{
  x+=rhs.x;
  y+=rhs.y;
  z+=rhs.z;
  return *this;
}

inline
ThreeVector& ThreeVector::operator-= (const ThreeVector& rhs)
{
  x-=rhs.x;
  y-=rhs.y;
  z-=rhs.z;
  return *this;
}

inline
const ThreeVector ThreeVector::operator- (const ThreeVector& rhs) const
{
  return ThreeVector(x - rhs.x, y - rhs.y, z - rhs.z);
}

inline
const ThreeVector ThreeVector::operator+ (const ThreeVector& rhs) const
{
  return ThreeVector(x + rhs.x, y + rhs.y, z + rhs.z);
}

inline
bool ThreeVector::operator== (const ThreeVector& rhs) const
{
  return ((rhs.x == this->x) && (rhs.y == this->y) && (rhs.z == this->z));
}



// friend operators

inline
const ThreeVector operator* (const double& factor, const ThreeVector& rhs)
{
  return ThreeVector(rhs.x*factor, rhs.y*factor, rhs.z*factor);
}

inline
const ThreeVector operator/ (const ThreeVector& rhs, const double& factor)
{
  return ThreeVector(rhs.x/factor, rhs.y/factor, rhs.z/factor);
}




// others

inline
double ThreeVector::DottedWith(const ThreeVector& rhs) const
{
  return (this->x * rhs.x + this->y * rhs.y + this->z * rhs.z);
}

inline
const ThreeVector ThreeVector::CrossedWith ( const ThreeVector& rhs) const
{
  return ThreeVector( this->y * rhs.z - rhs.y * this->z, this->z * rhs.x - rhs.z * this->x, this->x * rhs.y - rhs.x * this->y);
}

inline
const ThreeVector ThreeVector::GetDirection() const
{
  return (1.0 / this->GetLength() * (*this));
}


inline
double ThreeVector::GetAngleTo(const ThreeVector& rhs) const
{
  double arg = this->DottedWith(rhs)/(this->GetLength()*rhs.GetLength()); 
  if (arg>1.0) { return 0.0; } else { return acos(arg); }
}


inline
double ThreeVector::GetLength() const
{
  return sqrt(x*x + y*y + z*z);
}





// friends

inline
ostream& operator<<(ostream& os, const ThreeVector& rhs)
{
  os << rhs.x << "\t" << rhs.y << "\t" << rhs.z;
  return os;
}

inline
istream& operator>>(istream& is, ThreeVector& rhs)
{
  is >> rhs.x >> rhs.y >> rhs.z;
  return is;
}



#endif
