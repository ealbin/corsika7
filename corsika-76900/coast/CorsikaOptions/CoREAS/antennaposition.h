/***************************************************************************
 *   Copyright (C) 2005 by Tim Huege                                       *
 *   tim.huege@kit.edu                                                     *
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

#ifndef ANTENNAPOSITION_H
#define ANTENNAPOSITION_H

#include <string>
#include "threevector.h"

using std::string;

/**
@author Tim Huege
*/
class AntennaPosition{
public:
  AntennaPosition(ThreeVector parPosition, string parName, string parMode = "", string parModeData1 = "", string parModeData2 = "");
  ~AntennaPosition();

  ThreeVector GetPosition() const { return itsPosition; }
  string GetName() const { return itsName; }

  // access to special data fields
  string GetModeString() const { return itsMode; }
  string GetModeData1() const { return itsModeData1; }
  string GetModeData2() const { return itsModeData2; }

private:
  ThreeVector itsPosition;
  string itsName;
  string itsMode;
  string itsModeData1;
  string itsModeData2;
};

#endif
