
/* pllselectecuts.c:
**************************************************************************
complink:
       gcc -O0 -fbounds-check pllselectecuts.c -o pllselectecuts -lm
usage:
       ./pllselectecuts [Energy]
*************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>

void energyratio( double Engy, double Rati, double *Erat);

/* main program to display parallel corsika submit parameters:
*************************************************************************/

int main(int argc, char *argv[])
{
   double Energy=7.5e7, Eratio, Ratio, elog10, cfactor=2.6123;
   double qselect[30][4]={
          {  0.00,     0.,    0.,   0.},
          { 14.00,    56.,    6.,   0.}, { 14.25,    56.,    7.,   0.},
          { 14.50,    56.,    8.,   0.}, { 14.75,    56.,    9.,   0.},
          { 15.00,    56.,   13.,   0.}, { 15.25,    56.,   17.,   0.},
          { 15.50,    56.,   20.,   0.}, { 15.75,    56.,   30.,   0.},
          { 16.00,    56.,   60.,   0.}, { 16.25,    56.,   70.,   0.},
          { 16.50,    56.,   80.,   0.}, { 16.75,    56.,  100.,   0.},
          { 17.00,    84.,  120.,   0.}, { 17.25,    84.,  150.,   0.},
          { 17.50,    84.,  270.,   0.}, { 17.75,    84.,  340.,   0.},
          { 18.00,    84.,  370.,   0.}, { 18.25,   140.,  700.,   0.},
          { 18.50,   196.,  800.,   0.}, { 18.75,   252.,  900.,   0.},
          { 19.00,   308., 1100.,   0.}, { 19.25,   308., 1700.,   0.},
          { 19.50,   476., 2700.,   0.}, { 19.75,   616., 3200.,   0.},
          { 20.00,   924., 3500.,   0.}, { 20.25,   924., 4000.,   0.},
          { 20.50,   924., 4300.,   0.}, { 20.75,   980., 4310.,   0.},
          { 21.00,   980., 4320.,   0.}};
   int ieng, icnt, iprc, irat, ifiles, mfiles, iengy;

   // check execution argument (or use default value):
   if ( 1 <= argc && argc <= 3 ) {
   if ( argc >= 2 ) Energy = atof(argv[1]);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// for( ieng=1; ieng<38; ieng++ ) { Energy = 15.875 + 0.125001 * ieng;

   if ( Energy <  23.456789 ) Energy = pow(10.,Energy); // value to eV.
   if ( Energy >= 1.2345e12 ) Energy = 1.e-9 * Energy; // value to GeV.
   elog10 = 9.+log10(Energy); // log10 of GeV.
   // - - - - energy < 1.00e5 GeV:
   if ( Energy < 1.0e5 ) {
      printf("  \n");
      printf("          Energy = %10.4e GeV, log10(Energy/eV) = %7.4f\n",
         Energy,elog10);
      printf("               ==> run shower simulation as single"
             " sequential job \n");
      printf("                   ./corsika74567_..... < inp321023"
             " > DAT321023.lst & \n");
   }
   // - - - - energy >= 1.00e5 GeV:
   else {
      printf("                                                  "
             "          subshow \n");
      printf("  log10(E)    Energy    nproc  ecutmin     ecutmax"
             "  cpumin  nfiles    Ratio\n \n");
      icnt = (0.4484 + elog10 - 14.) / 0.25;
      iengy = elog10;
      ifiles = 1; 
      /* * * * * calculate submit quantities of original energy */
      iprc = qselect[icnt][1];
      // for( irat=-1; irat<=1; irat++ ) {  
      for( irat=-1; irat<=7; irat+=2 ) {  
         Eratio = cfactor/4.5678 * Energy / iprc;
         Ratio = Energy / Eratio * (1.4567+0.4567*irat);
         energyratio( Energy, Ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.5e7; 
         ifiles =  cfactor*Energy/Eratio; 
         if ( cfactor*Energy/Eratio < 34567. && mfiles != ifiles ) 
           if ( Eratio > 1.e4 )  
             printf(
             " %8.4f %12.4e %5d. %7.0f. %11.1e %6.0f. %6.0f. %8.1f \n",
               elog10, Energy, iprc, Eratio*1.e-3, Eratio,
               qselect[icnt][2], cfactor*Energy/Eratio, Energy/Eratio);
           else 
             printf(
             " %8.4f %12.4e %5d. %8.1f %11.1e %6.0f. %6.0f. %8.1f \n",
               elog10, Energy, iprc, Eratio*1.e-3, Eratio,
               qselect[icnt][2], cfactor*Energy/Eratio, Energy/Eratio); 
         printf("                                                     "
                "      %6.0f.\n",2.39*Ratio+5.234);
         mfiles = ifiles;
      }
   }

   // printf("  \n");
   // printf("                                                  "
   //        "          subshow \n");
   // printf("  log10(E)    Energy    nproc    ectcut      ectmax"
   //        "    cpumin  nfiles    Ratio\n \n");

// } // end-of optional for-loop for testing selections.

   printf("  \n");

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   }

   return 0;
}

/* energyratio.c:
**************************************************************************
   function to calculate Ratio=energy/ectmax
**************************************************************************
      double qvalue[13]={1.0,1.5,2.0,2.5,3.0,3.5,4.,5.,6.,7.,8.,9.,10.};*/

void energyratio( double Energy, double Ratio, double *Eratio) {

   double qlog10[13]={0.0, 0.17609125905568124, 0.30102999566398120,     
      0.39794000867203760, 0.47712125471966244, 0.54406804435027567,     
      0.60205999132796240, 0.69897000433601886, 0.77815125038364363,     
      0.84509804001425681, 0.90308998699194354, 0.95424250943932487}; 
   double eratlog, fratlog, dlog;

   int i, iratlog, ilog;

      eratlog = log10( Energy / Ratio );
      iratlog = eratlog;
      fratlog = eratlog - iratlog;
      dlog = 0.5 * qlog10[1];
      for( i=1; i<=12; i++ ) {
         if ( fabs(qlog10[i]-fratlog) < dlog ) {
            dlog = fabs(qlog10[i]-fratlog);
            ilog = i;
         } 
      }
      eratlog = qlog10[ilog] + iratlog;

   *Eratio = pow(10.,eratlog);

}

/****************************************************************************
  log10(E)    Energy    nproc    ectcut      ectmax    cpumin  nfiles   Ratio
 
  16.0000   1.0000e+07    16.     250.0    2.500e+05      90.     98.    40.0
  16.0000   1.0000e+07    16.     200.0    2.000e+05      90.    123.    50.0
 
  16.5000   3.1623e+07    24.     500.0    5.000e+05     120.    155.    63.2
  16.5000   3.1623e+07    24.     350.0    3.500e+05     120.    222.    90.4
 
  17.0000   1.0000e+08    32.    1000.0    1.000e+06     150.    246.   100.0
  17.0000   1.0000e+08    32.     900.0    9.000e+05     150.    273.   111.1
 
  17.5000   3.1623e+08    64.    2000.0    2.000e+06     360.    388.   158.1
  17.5000   3.1623e+08    64.    1500.0    1.500e+06     360.    518.   210.8
 
  18.0000   1.0000e+09    96.    5000.0    5.000e+06     500.    491.   200.0
  18.0000   1.0000e+09    96.    3500.0    3.500e+06     500.    702.   285.7

  18.5000   3.1623e+09   120.   15000.0    1.500e+07    1400.    518.   210.8
  18.5000   3.1623e+09   160.   10000.0    1.000e+07    1400.    777.   316.2

  19.0000   1.0000e+10   240.   25000.0    2.500e+07    2000.    983.   400.0
  19.0000   1.0000e+10   240.   15000.0    1.500e+07    2000.   1638.   666.7

  log10(E)    Energy    nproc    ectcut      ectmax    cpumin  nfiles   Ratio
****************************************************************************/
