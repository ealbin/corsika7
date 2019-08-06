/***************************************************************************
                          threevector.cpp  -  description
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

#include "threevector.h"


ThreeVector::~ThreeVector()
{}

ThreeVector& ThreeVector::Rotate(const double& parAngle, const Axis parAxis)
{
  ThreeVector srcVec(*this);        // needed for upcoming calculations
  double cosAngle = cos(parAngle);  // for performance
  double sinAngle = sin(parAngle);
  switch (parAxis)
  {
    case  XAxis:
          {
            y=cosAngle*srcVec.GetY()-sinAngle*srcVec.GetZ();
            z=sinAngle*srcVec.GetY()+cosAngle*srcVec.GetZ();
          }
          break;
    case  YAxis:
          {
            x=cosAngle*srcVec.GetX()+sinAngle*srcVec.GetZ();
            z=-sinAngle*srcVec.GetX()+cosAngle*srcVec.GetZ();
          }
          break;
    case  ZAxis:
          {
            x=cosAngle*srcVec.GetX()-sinAngle*srcVec.GetY();
            y=sinAngle*srcVec.GetX()+cosAngle*srcVec.GetY();
          }
          break;
    default: throw Exception();
          break;
  }

  return *this;
}
