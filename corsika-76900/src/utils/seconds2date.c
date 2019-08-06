/* seconds2date.c:
**************************************************************************
   calculates the date of the total number of seconds since 01.01.1970.
**************************************************************************
   g++ -fbounds-check seconds2date.c -o seconds2date
*************************************************************************/
#define _ISOC99_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

int main(int argcnt, char *argvec[]) 
{
   double dsecs;
   int jtage[13]={0,31,28,31,30,31,30,31,31,30,31,30,31}, jm,
       ijahr, imont, ileap, itage, ihour, imint, isecs, msecs, mdays, mtage;

// - - - - - - - - less than 2 arguments or too many arguments: 
   if ( argcnt < 2 || argcnt > 2 ) 
      printf("\n        usage:  ./seconds2date  seconds \n \n");

// - - - - - - - - program + 3 numbers (= 4 arguments):
   else if ( argcnt == 2 ) {
   dsecs = atof( argvec[1] );
   if ( 0. <= dsecs && dsecs < 2147483648. ) {
      msecs = int( dsecs );
      mdays = msecs;
      isecs = mdays%60; // printf("           isecs = %d \n",isecs);
      mdays = mdays - isecs;
      mdays = mdays / 60;
      imint = mdays%60; // printf("           imint = %d \n",imint);
      mdays = mdays - imint; 
      mdays = mdays / 60;
      ihour = mdays%24; // printf("           ihour = %d \n",ihour);
      mdays = mdays - ihour;
      mdays = mdays / 24; // printf("           mdays = %d \n",mdays);
      ileap = int( (671.+mdays) / 1461. ); 
      ijahr = 1970 + int( mdays / 365.25 );
      // printf("           ijahr = %d \n",ijahr);
      itage = mdays - (ijahr-1970) * 365 - ileap + 1;
      if ( itage > 365 ) {
         itage = 1;
         ijahr = ijahr + 1;
      }
      jtage[2] = 28;
      if ( ijahr%4 == 0 ) {
         itage++;
         jtage[2] = jtage[2] + 1;
      }
      if ( ijahr%4 == 0 ) itage++;
      if ( ijahr%4 == 0 ) jtage[2] = jtage[2] + 1; 
      // printf("           itage = %d \n",itage);
      imont = 0;
      mtage = 0;
      while( mtage < itage ) {
         imont++;
         mtage = mtage + jtage[imont]; 
      }
      for( jm=1; jm<imont; jm++ ) itage = itage - jtage[jm];
      printf("\n        %.2d.%.2d.%4d  %.2d:%.2d:%.2d    totalseconds"
         "  %d \n \n", itage, imont, ijahr, ihour, imint, isecs, msecs);
   }
   else {
      if ( dsecs >= 2147483648. ) 
         printf("\n        warning:  seconds invalid ...>=2^31... \n \n");
      else if ( dsecs < 0. )
         printf("\n        warning:  seconds invalid ...<0... \n \n");
   }
   } // end-of argcnt == 2.

// - - - - - - - - end-of program.
   return 0;
}
 
