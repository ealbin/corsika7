/*-- Author :    D. Horn & M. Raabe, DESY        13/01/97*/
/*====================================================================*/

/*    subroutine timer( iseco )                                       */

/*--------------------------------------------------------------------*/
/*  c-routine to read out system time in sec.                         */
/*  To be used in combination with CORSIKA                            */
/*  for old g77 compiler (version <0.5.21, with gcc 2.7.2.3).         */
/*  Not needed with newer g77 compilers.                              */
/*  To be compiled with g77 -funderscoring                            */
/*  This subroutine is called from AAMAIN and SEKDAT                  */
/*  Argument:                                                         */
/*   iseco  =   seconds since 00:00:00 GMT on Jan. 1th 1970           */
/*--------------------------------------------------------------------*/

#include <sys/types.h>
#include <sys/timeb.h>

void timer_(long *UnixSysTime)
{
  int dummy = 0;
  struct timeb Time_Struct;
  dummy=ftime(&Time_Struct);
  *UnixSysTime=Time_Struct.time;
}
/*--------------------------------------------------------------------*/
