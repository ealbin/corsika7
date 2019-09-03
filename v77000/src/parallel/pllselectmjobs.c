
/* pllselectmjobs.c:
*  -----------------
*     for given [mjobs] (i.e. number of subshowers per processor) as
*     argument, examples of the important parallel parameters `dectcut`
*     and `dectmax` and others to run a parallel corsika simulation to
*     get a well estimated number of resulting particle data files,
*     directly connected to the ratio `energy/dectmax`.
************************************************************************
* compilation:
      gcc -fbounds-check pllselectmjobs.c -o pllselectmjobs -lm
* execution:
      ./pllselectmjobs  [mjobs]
************************************************************************
   log10(E)  nproc    dectcut   dectmax  cputime  nfiles  Ratio 
     17.67      32    3.0E+04   3.0E+06   159.7     377   154.5
     17.67      32    2.0E+04   2.0E+06   154.4     540   231.7
     17.67      32    2.0E+04   1.0E+06   146.6    1057   463.4
     17.67      64    2.0E+04   1.0E+06    73.8    1057   463.4
     17.67      96    2.0E+04   1.0E+06    50.9    1057   463.4
     17.67     128    2.0E+04   1.0E+06    39.4    1057   463.4
***********************************************************************/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>

void energyratio( double Engy, double Rati, double *Erat);

/* main program to display simulation parameters:
*************************************************************************/

int main(int argc, char *argv[])
{
   double Energy=7.50e7, Eratio, Ratio, ratio, elog10, efiles,
          efact=1./2.4321, efactor, eratlog;
   double qeratio[8]={0.,111.,144.,222.,433.,677.,1133.,1322.};
   double qprtime[8]={0.,100.,240.,480.,1500.,2400.,3600.,4300.};
   int i, ilog, iprc, iratlog, iprocs, mjobs=6;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - call with 1 or 2 arguments:
   if ( 1 <= argc && argc <= 3 ) {
      if ( argc == 2 ) {
         Energy = atof(argv[1]);
         if ( Energy < 1234. ) mjobs = Energy;
      }
      if ( argc == 3 ) {
         Energy = atof(argv[1]);
         mjobs = atoi(argv[2]);
      }     
      for( i=2; i<=40; i++ ) {
      Energy = 15.75002 + 0.125 * i;
      if ( Energy <  23.456789 ) Energy = pow(10.,Energy); // to eV.
      if ( Energy >= 1.2345e12 ) Energy = 1.e-9 * Energy; // to GeV.
      elog10 = 9.+log10(Energy); // log10 of GeV.
      efactor = 16.+8.*(elog10-15.99)*2.;
      efactor = efactor - 2.*((i-2)%4);
      iprc = efactor; // now number of processors.
      if ( Energy < 1.0e7 ) {
         printf("  \n");
         printf("          Energy = %10.4e GeV   log10(Energy/eV) = %7.4f \n",
            Energy,elog10);
         printf("               ==> run shower simulation as single"
                " sequential job: \n");
         printf("                   ./corsika72456_..... < inp310023"
                " > DAT310023.lst \n \n");
      }
      else {
         efiles = 1. * mjobs * iprc; // 6 jobs per processor.
         ratio = efact * efiles;
         if ( elog10 < 17.5 ) ilog = 1;
         else { if ( elog10 < 18.0 ) ilog = 2;
            else { if ( elog10 < 18.5 ) ilog = 3;
               else { if ( elog10 < 19.0 ) ilog = 4; 
                  else { if ( elog10 < 19.5 ) ilog = 5; 
                     else { if ( elog10 < 20.0 ) ilog = 6;
                        else { ilog = 7;
         }  }  }  }  }  } 
         printf("  \n");
         printf("  log10(E)    Energy    nproc   dectcut     dectmax"
                "   cputime  nfiles   Ratio \n \n");
      // - - - - - first choise:
         eratlog = log10( Energy / ratio );
         Eratio = pow(10.,eratlog);
         if ( elog10 < 18. )
         printf(" %8.4f   %10.4e   %3d  %9.1f  %11.1f    %4.0f   %5.0f  %7.1f \n",
                  elog10, Energy, iprc, Eratio*1.e-3, Eratio,
                  qprtime[ilog], 2.4567*Energy/Eratio, Energy/Eratio);
         else
         printf(" %8.4f   %10.4e   %3d  %9.1f  %11.3e    %4.0f   %5.0f  %7.1f \n",
                  elog10, Energy, iprc, Eratio*1.e-3, Eratio,
                  qprtime[ilog], 2.4567*Energy/Eratio, Energy/Eratio);
      // - - - - - second choise:
         energyratio( Energy, ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.4e7; 
         if ( 1./efact*Energy/Eratio < 5000. ) 
         if ( elog10 < 18. )
         printf(" %8.4f   %10.4e   %3d  %9.1f  %11.1f    %4.0f   %5.0f  %7.1f \n",
                  elog10, Energy, iprc, Eratio*1.e-3, Eratio,
                  qprtime[ilog], 2.4567*Energy/Eratio, Energy/Eratio);
         else
         printf(" %8.4f   %10.4e   %3d  %9.1f  %11.3e    %4.0f   %5.0f  %7.1f \n",
                  elog10, Energy, iprc, Eratio*1.e-3, Eratio,
                  qprtime[ilog], 2.4567*Energy/Eratio, Energy/Eratio);
      // - - - - - third choise:
         Ratio = qeratio[ilog];
         energyratio( Energy, Ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.4e7; 
         if ( 1./efact*Energy/Eratio < 6000. ) 
      // - - - - - fourth choise:
         ratio = 0.5*(Ratio+ratio);
         energyratio( Energy, ratio, &Eratio );
         if ( elog10 > 19.5 ) Eratio = Eratio - 1.4e7; 
         if ( 1./efact*Energy/Eratio < 6000. ) 
         if ( elog10 < 18. )
         printf(" %8.4f   %10.4e   %3d  %9.1f  %11.1f    %4.0f   %5.0f  %7.1f \n",
                  elog10, Energy, iprc, Eratio*1.e-3, Eratio,
                  qprtime[ilog], 2.4567*Energy/Eratio, Energy/Eratio);
         else
         printf(" %8.4f   %10.4e   %3d  %9.1f  %11.3e    %4.0f   %5.0f  %7.1f \n",
                  elog10, Energy, iprc, Eratio*1.e-3, Eratio,
                  qprtime[ilog], 2.4567*Energy/Eratio, Energy/Eratio);
      }
      } // end-of for loop for testing selection.
   }

}

/* energyratio.c:
** function to calculate ratio=energy/dectmax
*************************************************************************/

void energyratio( double Energy, double ratio, double *Eratio) {

   double qvalue[13]={1.0,1.5,2.0,2.5,3.0,3.5,4.,5.,6.,7.,8.,9.,10.};

   double qlog10[13]={0.0, 0.17609125905568124, 0.30102999566398120,     
      0.39794000867203760, 0.47712125471966244, 0.54406804435027567,     
      0.60205999132796240, 0.69897000433601886, 0.77815125038364363,     
      0.84509804001425681, 0.90308998699194354, 0.95424250943932487}; 

   double eratlog, fratlog, dlog;

   int i, iratlog, ilog;

      eratlog = log10( Energy / ratio );
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

