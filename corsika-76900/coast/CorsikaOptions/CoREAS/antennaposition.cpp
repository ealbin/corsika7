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


#include "antennaposition.h"

AntennaPosition::AntennaPosition(ThreeVector parPosition, string parName, string parMode, string parModeData1, string parModeData2)
  : itsPosition(parPosition), itsName(parName), itsMode(parMode), itsModeData1(parModeData1), itsModeData2(parModeData2)
{
}


AntennaPosition::~AntennaPosition()
{
}


