
// bcreinpcont.cpp:
// =======================================================================
//    g++ -O0 -fbounds-check bcreinpcont.cpp -o bcreinpcont -lm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    Automatic creation of successive steering files to run corsika
//    simulations including corresponding shell script files for
//    submission to IKP computing cluster in Bldg. 425 of KIT_CN;
//    monthly atmospheres, M. Roth, A. Bridgeman.
// usage:
//       ./bcr.sh
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    #!/bin/bash
//    # creation script `./bcr.sh`
//    g++ -O0 -fbounds-check bcreinpcont.cpp -o bcreinpcont -lm
//    ./bcreinpcont > bcreinpcont.tabout
//    cat bcreinpcont.tabout
//    chmod 750 saug*
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                                      juergen.oehlschlaeger@kit.edu
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cfloat>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
using namespace std;

const double cpi=3.141592653589793, cpi180=cpi/180., c180pi=180./cpi;

const double velight = 29.9792458; // velocity of light in cm/nsec.

const double AATM[6] = { 0.e0,       // param for US standard atmosph.
                   -186.5562e0, -94.919e0, 0.61289e0, 0.e0, 1.128292e-2 };
const double BATM[6] = { 0.e0,       // param for US standard atmosph.
                1222.6562e0, 1144.9069e0, 1305.5948e0, 540.1778e0, 0.e0 };
const double CATM[6] = { 0.e0,       // param for US standard atmosph.
              994186.38e0, 878153.55e0, 636143.04e0, 772170.16e0, 1.e-9 };

int *RandomSequenceInit(int);

float RandomUniform(float randmax, int *mrandseq);

double heightcm(double g); // calculate height (cm) for a given thickness.

double thickgr(double h); // thickness of atmosphere depending on height.

double ap[13][17]={ 
    {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
    {0., -129.86987635, -13.912578027, 1.137905219,-4.5501979957e-4,
         1170.0778448, 1310.6961298, 1490.6966654, 503.61356785,
       971950.03798, 682326.92144, 615751.05849, 795110.76,
      1070000., 1460000., 3660000., 10000000.},
    {0., -140.02703037, -32.154019677,1.3324218584,-3.4759783013e-4,
         1177.9948918, 1226.57938, 1544.3942757, 529.86754375,
       990841.49847, 746882.89502, 608546.43122, 787929.8796,
       980000., 1430000., 3550000., 10000000.},
    {0., -140.79368856, -36.408818366,1.4042651261,-2.4583883619e-4,
         1177.4425334, 1189.02147, 1556.0070888, 539.33390781,
       984028.35984, 766300.21844, 604631.06014, 782862.57672,
       900000., 1450000., 3500000., 10000000.},
    {0., -133.89496569, -47.089899137,0.83278306486,-2.105436733e-4,
         1174.1385323, 1215.7355043, 1425.6918596, 503.93169671,
       967164.95864, 768155.10252, 621852.25389, 785601.62044,
       990000., 1210000., 3800000., 10000000.},
    {0., -152.40981237,-12.095554799,0.85136995339,-4.5928490669e-5,
         1190.3622315, 1265.7764387, 1437.9375049, 524.79662686,
       985619.00774, 681017.97599, 617989.22267, 776008.70875,
       950000., 1500000., 3750000., 10000000.},
    {0., -187.01662033, -25.280960312, 0.89432403304, 3.47276315e-5,
         1227.6177212, 1178.8761669, 1445.4602411, 535.28753109,
       1012299.9934, 732521.03651, 614783.70338, 771077.62332,
       750000., 1450000., 3700000., 10000000.},
    {0., -136.34855518, -20.914318979, 0.87280369409,3.732759787e-5,
         1182.905474, 1213.2533581, 1413.9801198, 547.58686206,
       947431.21788, 710893.29507, 619412.05748, 769605.79586,
       890000., 1400000., 3700000., 10000000.},
    {0., -139.77721439,-25.948986634,0.70618984827,-9.6699826826e-6,
         1184.5463072, 1205.4324091, 1366.057468, 542.40761127,
       953640.44565, 723050.32076, 629158.2012, 772372.2386,
       880000., 1280000., 3830000., 10000000.},
    {0., -158.54644384,-95.110327312,0.76408408672,-1.5130518952e-4,
         1200.3725926, 1165.8522743, 1373.2445748, 521.58648708,
       986103.78389, 881190.80912, 628689.50007, 781005.59304,
       690000., 1100000., 3820000., 10000000.},
    {0., -150.1396913, -58.140415154,0.87315735073,-3.6221952628e-4,
         1189.9866918, 1193.1675932, 1422.7442647, 490.3837468,
       990579.96621, 800597.35355, 623333.46715, 793328.53093,
       900000., 1200000., 3800000., 10000000.},
    {0., -194.83705256, -68.148895317,0.91697530761,-5.928015954e-4,
         1231.0108397, 1164.4301252, 1412.7462127, 460.92737873,
      1042817.5709, 832823.86466, 624149.1426, 805671.77828,
       700000., 1200000., 3810000., 10000000.},
    {0., -152.92783842, -37.51920314, 1.5800882678,-5.9051575509e-4,
         1190.2311741, 1189.4576063, 1575.714709, 484.1030213,
       997989.25368, 766606.85869, 599559.45014, 802421.45312,
       870000., 1450000., 3480000., 10000000.} };

// - - - - - - set number of files, showers, energies - - - - - - - - - -
double theta[10]={ 0., 0., 12., 22., 38., 65., 0., 0., 0., 0.},
       qengy[10]={ 0., 1.00e8, 3.1623e8, 1.00e9, 3.1623e9,
                      1.00e10, 3.1623e10, 1.e11, 3.1623e11, 0.};        

int mprim[10]={0, 14, 5626, 1206, 2814, 402, 1, 1, 1, 1};

// = = = = = = = = = = = = = = = = = main  = = = = = = = = = = = = = = = =
//     write steering files, submit scripts, and submit command file.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int main(int argcnt, char *argvec[])
{
   FILE *coutext, *coutsub;

   char choutext[256], choutname[256], chexecute[256], chsubmpath[256],
        cskinput[16], cskscript[16], chdatfile[16], chdatlist[16],
        chsubmit[24], chmodel[8];

   double obsgrm, obslev, engya, engyb, eslope, engyval, thadron, tweight,
       tradius, themin, themax, rdnr, ELL, EUL, SLEX;
  
   int iset, iegy, iprm, ithe, iatm, isnr, iprim, ishow, irun, isnmax,
       msrun, mstrt0, mevcnt=5, mtext, iranseq=23, *nranseq;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if ( argcnt >= 2 ) mevcnt = atoi(argvec[1]);
   nranseq = RandomSequenceInit( iranseq );

// - - - - - - some initialisations for corsika files - - - - - - - - - -
   ishow = 1;
   eslope = -1.;
   obslev = 1452.e2;
   obsgrm = thickgr(obslev); 
   isnmax = mevcnt; // default at 5 (for test cases).

// - - - - - - define (file-)names:
   strcpy(chsubmpath, "/cr/auger02/joe/corsika.trunk/run"); // length 32.
   strcpy(chexecute, "corsika76400_stnd_QGSII4_gheisha"); // length 32.
   strcpy(chsubmit, "saugsubmit.sh");
   coutsub = fopen(chsubmit, "w");

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - set of showers for each primary, energy, angle, atmos:
   for( iset=1; iset<=1; iset++ ) {

// - - - - - - primary particle energy loop - - - - - - - - - - - - - - - 
   for( iegy=1; iegy<=1; iegy++ ) {
      engya = qengy[iegy];
      engyb = qengy[iegy];
      printf("\n - - - - - - - - - - - - - - - - - - - \n");
      if ( engya < engyb )
         printf("   energy=%9.2e ...%9.2e \n", engya+6.,engyb+6.);
      else 
         printf("   energy=%9.2e \n", engya+6.);

// - - - - - - primary particle type loop - - - - - - - - - - - - - - - -
   for( iprm=1; iprm<=2; iprm++ ) {
      iprim = mprim[iprm];
      printf("          prim= %d \n",mprim[iprm]);

// - - - - - - angle theta loop - - - - - - - - - - - - - - - - - - - - -
   for( ithe=1; ithe<=1; ithe++ ) {
      themin = theta[ithe];
      themax = themin;
      printf("               theta=%6.2f \n",themin);
      mstrt0 = 437300 + isnmax*(iprm-1)+ 2*(ithe-1) + 10*(iegy-1);
      iatm = 0;

// - - - - - - shower steering file creation loop - - - - - - - - - - - -
   for( isnr=1; isnr<=isnmax; isnr++ ) {
      msrun = mstrt0 + isnr - 1;
      if ( engya < engyb ) { 
         // - - - - - - - calculate random energy value:
         // call rmmard( rdnr(1),1,1 ) ! calculate random energy value.
         rdnr = RandomUniform( 1., nranseq ); 
         if ( eslope != -1. ) {
            ELL = pow( engya, eslope+1.);
            EUL = pow( engyb, eslope+1.);
            SLEX = 1. / (eslope+1.);
            engyval = pow( rdnr*EUL+(1.-rdnr)*ELL, SLEX);
         }
         else 
            engyval = engya * pow( engyb/engya, rdnr);
      }
      else
         engyval = engya;
      // energy value now copied to quantity engyval.
      thadron = 100.;
      tweight = 1.e-6 * engyval;
      tradius = 1.;
      if ( engyval >= 1.00e+09 ) tradius = 10.;
      if ( engyval >= 3.16e+09 ) tradius = 32.; 
      if ( engyval >= 1.00e+10 ) tradius = 100.; 
      if ( engyval >= 3.16e+10 ) tradius = 120.; 
      if ( engyval >= 1.00e+11 ) tradius = 140.; 
      if ( isnr%10 == 1 ) iatm = iatm + 1; // change atm after 10 sh.
      if ( isnr == 1 || isnr == isnmax ) 
         printf("                     msrun= %d \n",msrun);
      irun = msrun;
      if ( irun > 999999 ) irun = irun%1000000;
      sprintf(cskinput, "aug%.6d", irun);
      sprintf(chdatfile, "DAT%.6d", irun);      
      sprintf(chdatlist, "DAT%.6d.lst", irun); 

   // - - - - - - - create corsika steering file:
      coutext = fopen(cskinput, "w");
      fprintf(coutext, "RUNNR %10d\n", irun);
      fprintf(coutext, "EVTNR %10d\n", 1);
      fprintf(coutext, "SEED  %10d         0         0\n", msrun*3);
      fprintf(coutext, "SEED  %10d         0         0\n", msrun*3+1);
      fprintf(coutext, "SEED  %10d         0         0\n", msrun*3+2);
   // fprintf(coutext, "SEED  %10d         0         0\n", msrun*3+3);
   // fprintf(coutext, "SEED  %10d         0         0\n", msrun*3+4);
   // fprintf(coutext, "SEED  %10d         0         0\n", msrun*3+5);
      fprintf(coutext, "PRMPAR%10d\n", iprim);
      if ( engya < engyb ) fprintf(coutext, "ESLOPE   %9.2f\n",eslope);
      fprintf(coutext, "ERANGE    %14.4e %14.4e\n", engyval,engyval);
      fprintf(coutext, "THETAP   %9.2f %9.2f\n", themin,themax);
      if ( themin > 0.1234 )
         fprintf(coutext, "PHIP     %9.2f %9.2f\n", 0.,360.);
      else
         fprintf(coutext, "PHIP          0.00      0.00\n");
      fprintf(coutext, "NSHOW %10d\n", ishow);
      fprintf(coutext, "OBSLEV %11.2fE2 %13.3f  g/cm^2\n",
              obslev*1.e-2, obsgrm);      
      fprintf(coutext, "ECUTS       %6.2f %6.2f %8.1e %8.1e\n",
              0.1,0.01,5.e-5,5.e-5);
   // - - - - - - - tabular of FIXHEI quantities:
   // fprintf(coutext, "FIXHEI    36233.e2       0      (5.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    32765.e2       0      (8.18 gr/cm^2)\n");
   // fprintf(coutext, "FIXHEI    33394.15e2     0     (10.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    25254.68e2     0     (25.25468 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    23457.35e2     0     (33.3 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    22346.845e2    0     (39.534 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    20831.93e2     0     (50.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    19876.e2       0     (58.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    18765.43e2     0     (68.956 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    18669.00e2     0     (70.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    16383.17e2     0    (100.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    13790.77e2     0    (150.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI    11954.18e2     0    (200.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI     9347.04e2     0    (300.0 g/cm^2)\n");
   // fprintf(coutext, "FIXHEI     8300.1876e2   0    (350.0 g/cm^2)\n");
      fprintf(coutext, "MAXPRT         1\n");
      fprintf(coutext, "ECTMAP     1.E11\n");
      fprintf(coutext, "RADNKG    200.E2\n"); 
      fprintf(coutext, "HADFLG    0    0    0    0    0    2\n");
      fprintf(coutext, "ELMFLG         T         T\n");
      if ( strstr(chexecute, "EPOS") ) { 
         fprintf(coutext, "EPOPAR input ../epos/epos.param\n");
         fprintf(coutext, "EPOPAR fname inics ../epos/epos.inics\n");
         fprintf(coutext, "EPOPAR fname iniev ../epos/epos.iniev\n");
         fprintf(coutext, "EPOPAR fname initl ../epos/epos.initl\n");
         fprintf(coutext, "EPOPAR fname inirj ../epos/epos.inirj\n");
         fprintf(coutext, "EPOPAR fname inihy ../epos/epos.ini1b\n");
         fprintf(coutext, "EPOPAR fname check none\n");
         fprintf(coutext, "EPOPAR fname histo none\n");
         fprintf(coutext, "EPOPAR fname data  none\n");
         fprintf(coutext, "EPOPAR fname copy  none\n");
      }
      fprintf(coutext, "MUMULT         T\n");
      fprintf(coutext, "MUADDI         T\n");
   // fprintf(coutext, "EMADDI         T\n");
      fprintf(coutext, "STEPFC        1.\n");
      fprintf(coutext, "HILOW       100.\n");
      fprintf(coutext, "THIN      1.e-06%10.0f.%7.0f.e2\n",
              tweight, tradius);
      fprintf(coutext, "THINH         1.%10.0f.\n", thadron);
      fprintf(coutext, "ATMOD   0            monthly-atmos-%.2d\n",iatm);
      fprintf(coutext, "ATMA     %15.8e %15.8e %15.8e %15.8e\n",
              ap[iatm][1], ap[iatm][2], ap[iatm][3], ap[iatm][4]);
      fprintf(coutext, "ATMB     %15.8e %15.8e %15.8e %15.8e\n",
              ap[iatm][5], ap[iatm][6], ap[iatm][7], ap[iatm][8]);
      fprintf(coutext, "ATMC     %15.8e %15.8e %15.8e %15.8e\n",
              ap[iatm][9], ap[iatm][10], ap[iatm][11], ap[iatm][12]);
      fprintf(coutext, "ATMLAY   %15.8e %15.8e %15.8e %15.8e\n",
              ap[iatm][13], ap[iatm][14], ap[iatm][15], ap[iatm][16]);
      fprintf(coutext, "LONGI          T      5.     T      T\n");
      fprintf(coutext, "MAGNET       19.46     -14.16\n");
      if ( engya <= 1.012345e8 )
         fprintf(coutext, "CORECUT       0.001e2\n");
      else
         fprintf(coutext, "CORECUT %11.2fe2\n",tradius);
      fprintf(coutext, "DIRECT %s\n", "./"); 
      fprintf(coutext, "HOST   iklx282\n");
      fprintf(coutext, "USER   you\n");
      fprintf(coutext, "EXIT\n");
   // - - - - - - - close corsika steering file.
      fclose(coutext);
 
   // - - - - - - - create corsika run script file (IKP computing cluster):
      sprintf(cskscript, "saug%.6d", msrun);
      strcpy(choutext, cskscript);
      coutext = fopen(choutext, "w");
      fprintf(coutext, "#!/bin/bash\n");
      fprintf(coutext, "#\n");
      fprintf(coutext, "#$ -cwd\n");
      fprintf(coutext, "#$ -j y\n");
      fprintf(coutext, "#$ -o %s\n",chsubmpath);
      fprintf(coutext, "#\n");
      fprintf(coutext, "cd %s\n",chsubmpath);
      fprintf(coutext, "./%s < %s > %s\n",chexecute,cskinput,chdatlist);

   // - - - - - - - close corsika run script file.
      fclose(coutext);

   // - - - - - - - write submit command to submit file:
      fprintf(coutsub, "qsub %s\n",cskscript); 
      fprintf(coutsub, "sleep 12\n");      

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   } // end-of loop isnr.

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   } // end-of loop ithe.

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   } // end-of loop iprm.

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   } // end-of loop iegy.

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   } // end-of loop iset.

   fclose(coutsub); // close submit script after all run numbers.
   printf("\n");

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   return 0;

}

// = = = = = = = = = = = = = Init Random Sequence  = = = = = = = = = = = =
// Initialize a random number sequence for the various random number
// generators. Instead of the original array in a common block for the
// seed, which is addressed via a sequence number, here only a pointer
// for the dynamically allocated array is tranfered to each call.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int *RandomSequenceInit(int iseed) {

  float uu[98], *u, *c, s, t, twom24;
  int i, j, k, l, m, ii, ij, jj, kl, *mrandseq;

  mrandseq = (int *)malloc(4*104); /* 104 ints for the mrandseq */
  u = (float *) &(mrandseq[3]);
  c = (float *) &(mrandseq[101]);         
  twom24 = (float) pow(2.,-24.); /* define as constant */
  mrandseq[1] = iseed;
  mrandseq[2] = 0L;
  mrandseq[3] = 0L;
  ij = mrandseq[1] / 30082;
  kl = mrandseq[1] - 30082 * ij;   
  i  = (ij/177) % 177 + 2;   
  j  = ij % 177 + 2;      
  k  = (kl/169) % 178 + 1;   
  l  = kl % 169;      

  for (ii=1; ii<=97; ii++) {
    s = (float)0.0;
    t = (float)0.5;
    for (jj=1; jj<=24; jj++) {
      m = ( ((i*j) % 179)*k) % 179;   
      i = j;            
      j = k;            
      k = m;            
      l = (53*l+1) % 169;      
      if ( ((l*m)%64) >= 32) s=s+t;   
      t = (float)0.5 * t;         
    }               
    uu[ii] = s;            
  }               
  mrandseq[102] = 97;
  mrandseq[103] = 33;
  c[0] = (float)362436.0 * twom24;
  for (jj=1; jj<=97; jj++) u[jj] = uu[jj];      

  return(mrandseq); /* return a sequence for the generator */

}

// = = = = = = = = = = = = Uniform Random Dist = = = = = = = = = = = = = =
// Uniform Random Distribution is copied from Cern RMMAR
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float RandomUniform(float randmax, int *mrandseq) {

  float uni,*u,*c,cd,cm,twom24;

  *( (int *) (&twom24) ) = 0x33800000L; /* twom24 = 2 ** -24 */
  *( (int *) (&cd) ) = 0x3ee99762L;     /* cd = 7654321. * 2 ** -24 */
  *( (int *) (&cm) ) = 0x3f7ffffdL;     /* cm = 16777213. * 2 ** -24 */

  u = (float *) &(mrandseq[3]);         
  c = (float *) &(mrandseq[101]);         
  uni = u[mrandseq[102]] - u[mrandseq[103]];   

  if (uni < 0.) uni = uni + (float)1.0;      
  u[mrandseq[102]] = uni;               

  mrandseq[102]--;            
  if (mrandseq[102] == 0) mrandseq[102] = 97;   

  mrandseq[103]--;            
  if (mrandseq[103] == 0) mrandseq[103] = 97;   

  c[0] = c[0] - cd;         
  if ( c[0] < 0. ) c[0] = c[0] + cm;   

  uni = uni - c[0];         
  if ( uni < 0. ) uni = uni + (float)1.0;   

  if (uni == 0.) {
    uni = twom24 * u[2];      
    if ( uni == 0.) uni = twom24 * twom24;   /* Only to be safe */
  }                    

  return(uni*randmax);
                                                                        
}

// = = = = = = = = = = = = = = = heightcm  = = = = = = = = = = = = = = = =
//     calculate height (cm) a.s.l. for a given thickness (gramms/cm**2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double heightcm(double g)
{
   double height;
   if (g > 631.1e0)
      { height = CATM[1] * log( BATM[1] / (g-AATM[1]) ); }
   else if (g > 271.7e0)
      { height = CATM[2] * log( BATM[2] / (g-AATM[2]) ); }
   else if (g > 3.0395e0)
      { height = CATM[3] * log( BATM[3] / (g-AATM[3]) ); }
   else if (g > 1.28292e-3)
      { height = CATM[4] * log( BATM[4] / (g-AATM[4]) ); }
   else
      { height = (AATM[5] - g) / CATM[5]; }
   return height;
}

// = = = = = = = = = = = = = = = = thickgr = = = = = = = = = = = = = = = =
//       thickgr (gramms/cm**2) of atmosphere depending on height (cm)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double thickgr(double h)
{
   double thickn;
   if (h < 4.e5)
      { thickn = AATM[1] + BATM[1] * exp(-h/CATM[1]); }
   else if (h < 1.e6)
      { thickn = AATM[2] + BATM[2] * exp(-h/CATM[2]); }
   else if (h < 4.e6)
      { thickn = AATM[3] + BATM[3] * exp(-h/CATM[3]); }
   else if (h < 1.e7)
      { thickn = AATM[4] + BATM[4] * exp(-h/CATM[4]); }
   else
      { thickn = AATM[5] - h*CATM[5]; }
   return thickn;
}

