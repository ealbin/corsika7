// =======================================================================
// 
//  c o r s 2 i n p u t . c p p
//  ===========================
//    creates a (possible) corsika steering (i.e. input) file using the
//    existing particle data output file, and also already available for 
//    corsika simulations of pure 64 bit architecture and a particle
//    data file of a parallel simulation (MPI used); data files of the
//    former `simprod` simulation sets are also correctly read;
//    but: some steering quantities may be defined additionally and
//    individually;
// -----------------------------------------------------------------------
// CompLink:
//      g++ -O0 -fbounds-check cors2input.cpp -o cors2inputcpp -lm
// Runprog:
//      ./cors2inputcpp irunnr
//      ./cors2inputcpp DATirunnr
//      ./cors2inputcpp /cr/auger02/joe/DAT111888
//      ./cors2inputcpp prE16M562_41
// -----------------------------------------------------------------------
// "standard" Corsika:
//      output format for particle output (blocklength = 22932+8 fixed)
//      each record consists of 21 subblocks of 273 words (5733).
// "thinning" Corsika:
//      output format for particle output (blocklength = 26208+8 fixed)
//      each record consists of 21 subblocks of 312 words (6552).
// -----------------------------------------------------------------------
//      Transfer indices of thinning elements:
//          evth(149) = thin(1); evth(151) = thin(2); evth(152) = thin(3);
//          evth(148) = thin(1)/thinh(1); evth(150) = thin(2)/thinh(2);
//      Reading DAT-file: if (pdata(l) .gt. 9.9e6) bs = mod(pdata(l),1.e5)
//      Reading CER-file: bs = pdata(l) ! bunchsize
// -----------------------------------------------------------------------
//            RUNH = 211285.2812500000000000;
//            EVTH = 217433.0781250000000000;
//            LONG =  52815.2968750000000000;
//            EVTE =   3397.3918457031250000;
//            RUNE =   3301.3325195312500000;
// -----------------------------------------------------------------------
//                                      juergen.oehlschlaeger@kit.edu
// -----------------------------------------------------------------------
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
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
using namespace std;

// int mposubstring( char *chbigstr, char *chsubstr) // see down.

double heightcm(double g); // calculate height (cm) for a given thickness.

double thickgr(double h); // thickness of atmosphere depending on height.

double thicksouth(double h); // thickness of atmosphere at South Pole.

void pamafill(); // fill particle masses up to 59_Ni.

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = global constants  = = = = = = = = = = = = =

   const int nrecstd = 22940;           // "standard corsika" record size,
                                        //    or  pdata[0] = 3.213456e-41.
   const int nrecthi = 26216;           // "thinning corsika" record size,
                                        //    or  pdata[0] = 3.672523e-41.
   const int ndim = 20000000;
   const int nsubblo = 21;
   const int nsblstd = 273;
   const int nsblthi = 312;
   const int nmaxrec = 2180000;              // valid up to 50 GBytes.
   const int numbstd = nrecstd / 4;          // =  5735.
   const int numbthi = nrecthi / 4;          // =  6554.

   float pdata[numbthi]; // to read a single corsika record thinned data.
   float rdata[numbthi]; // to keep first corsika record thinned data.
   float sdata[numbstd]; // to read a single corsika record standard data.
   float hdata[819]; // to read end of first corsika record thinned data.
   float zdata[2]; // to read record length information on 64bit machines. 

   const double AATM[6] = { 0.e0,       // param for US standard atmosph.
                   -186.5562e0, -94.919e0, 0.61289e0, 0.e0, 1.128292e-2 };
   const double BATM[6] = { 0.e0,       // param for US standard atmosph.
                1222.6562e0, 1144.9069e0, 1305.5948e0, 540.1778e0, 0.e0 };
   const double CATM[6] = { 0.e0,       // param for US standard atmosph.
              994186.38e0, 878153.55e0, 636143.04e0, 772170.16e0, 1.e-9 };

   const double velight = 29.9792458,   // velocity of light in cm/nsec.
         cpi=4.*atan(1.), cpi180=cpi/180., c180pi=180./cpi;

   const double
      aatmax[6]={ 0.,-129.86987,-13.912578,1.1379052,-4.5502e-4,1.12829e-2},
      batmax[6]={ 0., 1170.07784, 1310.69613, 1490.6966, 503.613568, 1.0},
      catmax[6]={ 0., 971950.04, 682326.92, 615751.06, 795110.76, 1.e-9},
      atmlay[6]={ 0., 10.7e5, 14.6e5, 36.6e5, 100.e5, 112.8292e5};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - number of showers of simprod simulation series:
//     1.0000E+14....1.7783E+14      28117 sh  =  18 * 1500 sh  + 1117 sh
//     1.7783E+14....3.1623E+14      15811 sh  =  15 * 1000 sh  +  811 sh
//     3.1623E+14....5.6234E+14       8891 sh  =  22 *  400 sh  +   91 sh
//     5.6234E+14....1.0000E+15       5000 sh  =  25 *  200 sh
//     1.0000E+15....1.7783E+15       2812 sh  =  28 *  100 sh  +   12 sh
//     1.7783E+15....3.1623E+15       1581 sh  =  26 *   60 sh  +   21 sh
//     3.1623E+15....5.6234E+15        889 sh  =  29 *   30 sh  +   19 sh
//     5.6234E+15....1.0000E+16        500 sh  =  25 *   20 sh
//     1.0000E+16....1.7783E+16        281 sh  =  25 *   11 sh  +    6 sh
//     1.7783E+16....3.1623E+16        158 sh  =  31 *    5 sh  +    3 sh
//     3.1623E+16....5.6234E+16         89 sh  =  44 *    2 sh  +    1 sh
//     5.6234E+16....1.0000E+17         50 sh  =  50 *    1 sh
//     1.0000E+17....1.7783E+17         28 sh  =  28 *    1 sh
//     1.7783E+17....3.1623E+17         16 sh  =  16 *    1 sh
//     3.1623E+17....5.6234E+17          9 sh  =   9 *    1 sh
//     5.6234E+17....1.0000E+18          5 sh  =   5 *    1 sh
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   const int nevent[20]={ 0, 1500, 1000, 400, 200, 100,
          60, 30, 20, 11, 5, 2, 1, 1, 1, 1, 1, 0, 0, 0};

   const int ncount[20]={ 0, 18, 15, 22, 25, 28,
          26, 29, 25, 25, 31, 44, 50, 28,16, 9, 5, 0, 0, 0};

   const int nclast[20]={ 0, 1117, 811, 91, 0, 12, 21, 19, 0,
          6, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0};

   const int nbit[4]={32,32,64,64};

   const int ipartyp[202] = {
                 0,1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
           21,22,23,24,25,26,27,28,29,30,31,32,17*0,50,51,52,53,54,55,
         56,57,58,59,60,61,62,63,64,65,66,67,68,69,0,71,72,73,74,41*0,
        116,117,118,119,120,121,122,123,124,125,126,127,128,0,130,131,
         132,133,134,2*0,137,138,139,140,141,142,143,144,145,3*0,149,
           150,151,152,153,154,155,156,157,3*0,161,162,163,7*0,171,
                 172,173,25*0,199,0,0 };

   const double parmas[202] = { 0.0e0,               // particle masses.
         0.e0     ,5.10998902e-4,5.10998902e-4, 0.e0     ,0.105658357e0,
         0.105658357e0,0.1349766e0,0.13957018e0,0.13957018e0,0.497672e0,
         0.493677e0, 0.493677e0, 0.93956533e0,0.93827200e0,0.93827200e0,
         0.497672e0 , 0.54730e0  , 1.115683e0 , 1.18937e0 , 1.192642e0 ,
         1.197449e0 , 1.31483e0  , 1.32131e0  , 1.67245e0 ,0.93956533e0,
         1.115683e0 , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31483e0  ,
         1.32131e0  , 1.67245e0  , 0.e0       ,  16*0.e0  , 0.78257e0  ,
         0.7690e0   , 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
         1.2331e0   , 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   ,
         1.2349e0   , 0.89610e0  , 0.89166e0  , 0.89166e0 , 0.89610e0  ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.54730e0  , 0.54730e0  , 0.54730e0  , 0.54730e0 , 0.105658e0 ,
         0.105658e0 , 0.0e0      , 0.0e0      , 0.0e0     ,  36*0.0e0  ,
         1864.5e0   , 1869.3e0   , 1869.3e0   , 1864.5e0  , 1968.6e0   ,
         1968.5e0   , 2979.7e0   , 2006.7e0   , 2010.0e0  , 2010.0e0   ,
         2006.7e0   , 2112.4e0   , 2112.4e0   , 3510.51e0 , 3096.87e0  ,
         1776.99e0  , 1776.99e0  , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 2284.9e0   , 2466.3e0   , 2471.8e0  , 2452.6e0   ,
         2451.3e0   , 2452.2e0   , 2574.1e0   , 2578.8e0  , 2697.5e0   ,
         0.e0       , 0.e0       , 0.e0       , 2284.9e0  , 2466.3e0   ,
         2471.8e0   , 2452.6e0   , 2451.3e0   , 2452.2e0  , 2574.1e0   ,
         2578.8e0   , 2697.5e0   , 0.e0       , 0.e0      ,     0.e0   ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      ,   6*0.e0   ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      ,  27*0.e0  }; 

   const char qnuclei[30][4]={" ",
         "H ","He","Li","Be","B ","C ","N ","O ","F ","Ne",
         "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca",
         "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu"};

   const double chargs[202] = { 0.,                 // particle charges.
           0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,
          +1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  ,+1.e0, 0.  ,
          -1.e0, 0.  ,-1.e0,-1.e0, 0.  , 0.  ,-1.e0, 0.  ,+1.e0, 0.  ,
          +1.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,11*0.,
           0.  ,+1.e0,-1.e0,+2.e0,+1.e0, 0.  ,-1.e0,-2.e0,-1.e0, 0.  ,
          +1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,41*0.,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,
          -1.e0, 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  ,
          +1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  ,+1.e0,+1.e0, 0.  ,+2.e0,
          +1.e0, 0.  ,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,-1.e0,-1.e0,
           0.  ,-2.e0,-1.e0, 0.  ,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,
          +2.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
          -2.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,22*0.}; 

   const char chptext[202][20] = {"                   ", 
      " gamma             "," positron          "," electron          ",
      "                   "," muon +            "," muon -            ",
      " pion 0            "," pion +            "," pion -            ",
      " Kaon 0 long       "," Kaon +            "," Kaon -            ",
      " neutron           "," proton            "," anti proton       ",
      " Kaon 0 short      ","                   "," Lambda            ",
      " Sigma +           "," Sigma 0           "," Sigma -           ",
      " Xi 0              ","  Xi -             "," Omega -           ",
      " anti neutron      "," anti Lambda       "," anti Sigma -      ",
      " anti Sigma 0      "," anti Sigma +      "," anti Xi 0         ",
                            // 0 to 30.
      " anti Xi +         "," anti Omega +      ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   "," omega             "," rho 0             ",
      " rho +             "," rho -             "," Delta ++          ",
      " Delta +           "," Delta 0           "," Delta -           ",
      " anti Delta --     "," anti Delta -      "," anti Delta 0      ",
                            // 31 to 60.
      " anti Delta +      "," Kaon * 0          "," Kaon * +          ",
      " Kaon * -          "," anti Kaon * 0     "," electron neutrino ",
      " anti elec neutrino"," muon neutrino     "," anti muon neutrino",
      "                   "," eta=>2*gamma      "," eta=>3*pi0        ",
      " eta=>pi+ pi- pi0  "," eta=>pi+ pi- gamma","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
                            // 61 to 90.
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   "," D 0               "," D +               ",
      " anti D -          "," anti D 0          "," D s +             ",
                            // 91 to 120.
      " anti D s -        "," eta c             "," D * 0             ",
      " D * +             "," anti D * -        "," anti D * 0        ",
      " D * s +           "," anti D * s -      ","                   ",
      " J/psi             "," tau +             "," tau -             ",
      " tau neutrino      "," anti tau neutrino ","                   ",
      "                   "," Lambda c +        "," Xi c +            ",
      " Xi c 0            "," Sigma c ++        "," Sigma c +         ",
      " Sigma c 0         "," Xi c prime +      "," Xi c prime 0      ",
      " Omega c 0         ","                   ","                   ",
      "                   "," anti Lambda c -   "," anti Xi c -       ",
                            // 121 to 150.
      " anti Xi c 0       "," anti Sigma c --   "," anti Sigma c -    ",
      " anti Sigma c 0    "," anti Xi c prime - "," anti Xi c prime 0 ",
      " anti Omega c 0    ","                   ","                   ",
      "                   "," Sigma c * ++      "," Sigma c * +       ",
      " Sigma c * 0       ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   "," anti Sigma c * -- ",
      " anti Sigma c * -  "," anti Sigma c * 0  ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
                            // 151 to 180.
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      "                   ","                   ","                   ",
      " Cherenkov photon  ","                   ","                   "};
                            // 181 to 201.

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

int mposubstring( char *chbigstr, char *chsubstr) {
   // find the position of a substring inside a bigger string:
   // ipos = mposubstring( satzc, chsub); // position counted from 0.

   int icnt=0, ifnd=0, isub=0, iret=0;

   if ( strlen(chbigstr) >= strlen(chsubstr) ) {
      iret = -1;
      for( icnt=strlen(chbigstr)-strlen(chsubstr); icnt>=0; icnt-- ) {
         ifnd = 1;
         for( isub=0; isub<int(strlen(chsubstr)); isub++ ) {
            if ( chbigstr[icnt+isub] != chsubstr[isub] ) {
               ifnd = 0;
               break;
            }
         }
         if ( ifnd == 1 ) break;
      }
      iret = icnt;
   }

   return iret;

} // end-of-function int mposubstring.

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   double drunh[275], devth[275], qevte[275], qdata[101], obslev, obsgrm,
          pmass[6002], psign[6002], gramlong, gramms, hilow, tradius;
 
   int ndatpar, nreclen, nsblock, lbit, iobs, irunnr, ilast,
       nxpix=750, nypix=650, nprimry, lobslev, jobslev, ifil, iruntyp=2,
       i, ii, is, ibit, iche, idat, isla, isubr, isho, iprm, itpa, iusc,
       lthi, iday, imon, iyer, lename, mpos, ienglow, ienghig, ltype8,
       lemadd, multhi;

   ifstream inputcors;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = main program= = = = = = = = = = = = = = = =

int main(int argcnt, char *argvec[])
{
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - character arrays and current date and time:
   FILE *coutext;
   char chfile[240], chfname[240], chrunnr[240], cskinput[240],
        chusnam[16], cfilcher[8], cfildatx[8], cfilpart[8], cfislash[8],
        cfipusla[8], cfiuscor[8], chflag[2][8]={"      F","      T"};
   time_t tseconds;
   tm *datetime;
   if ( argcnt > 1 ) irunnr = atoi(argvec[1]);
   tseconds = time(&tseconds);
   datetime = localtime(&tseconds);
   printf("\n        _date_of_analysis_ ");
   printf(" %.2d. %.2d. %.4d.  %.2d.%.2d%.2d  (hh.mmss %ld) \n",
      datetime->tm_mday, datetime->tm_mon+1, datetime->tm_year+1900,
      datetime->tm_hour, datetime->tm_min, int(tseconds)%60, tseconds);
   pamafill();  // fill particle masses up to 59_Ni.
   strcpy(chusnam, "joe");
   strcpy(cfilpart,".part");
   strcpy(cfildatx,"DAT");
   strcpy(cfipusla,"./");
   strcpy(cfislash,"/");
   strcpy(cfiuscor,"_");
 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - read next name of a corsika data file (without quotes):
   iusc = 0;
   if ( argcnt > 1 ) {
      ii = strlen(argvec[1]);
      if ( ii > 6 ) {
         strcpy(chfile,argvec[1]);
         strcpy(chrunnr, &chfile[idat+3]);
         strcpy(&chrunnr[6], "\0");
         irunnr = atoi(chrunnr);
      }
      else {
         irunnr = atoi(argvec[1]);
         sprintf(chfile, "DAT%.6d", irunnr);
      }
   }
   else {
      irunnr = 0; // no argument given:
      sprintf(chfile, "DAT%.6d", irunnr);
      printf("\n         usage: ./cors2inputcpp 99435\n");
      printf("                ./cors2inputcpp DAT099435\n");
      printf("                ./cors2inputcpp /cr/auger02/joe/DAT111888\n");
      printf("                ./cors2inputcpp prE16M562_41\n");
   }
   lename = strlen(chfile);
   idat = mposubstring(chfile, cfildatx); // idat=0: first char of chfile.
   isla = mposubstring(chfile, cfislash); // last slash (counted from 0).  
   itpa = mposubstring(chfile, cfilpart);
   iusc = mposubstring(chfile, cfiuscor);
   iche = mposubstring(chfile, cfilcher);
   if ( iusc > 9 ) printf("        idat=%d  iusc=%d  isla=%d  itpa=%d\n",
      idat,iusc,isla,itpa);
   // - - - - - - ignore `./` before particle data file name:
   if ( idat == 0 ) ilast = 0;
   else if ( idat == -1 ) ilast = 0;
   if ( isla >= 0 ) ilast = isla++;
   if ( ilast == 0 ) strcpy(chfname, chfile);
   else if ( ilast > 0 ) {
      ilast++;
      strncpy(chfname, &chfile[ilast], lename-ilast+1);
   }
   if ( idat+1 > 0 ) sprintf(cskinput,"%s.inp",&chfile[ilast]);
   else if ( iche+1 > 0 ) sprintf(cskinput,"%s.inp",&chfile[ilast]);
   if ( iusc+1 > 9 ) { // found `coe16m158_01`, `prE16M562_41` a.s.o.
      strcpy(chusnam, "simprod");
      sprintf(cskinput,"%s.inp",&chfile[ilast]);
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - check binary corsika/cherenkov particle data file:
   inputcors.open(chfile,ios::in); // allocate file.
   if ( inputcors.peek() == EOF ) { printf(
      "\n         file `%s` or irunnr=%d do not exist at this path!\n",
         chfile,irunnr);
      inputcors.clear();
      inputcors.close();
      goto continue29;
   }
   else { // read one record:
      inputcors.read((char*)&sdata,sizeof(sdata));
      if ( 217433.0 < sdata[nsblstd+1] && sdata[nsblstd+1] < 217433.2 ) {
         lbit = 32;
         lthi = 0;
         for( ii=numbstd-1; ii<=numbthi; ii++ ) rdata[ii] = 0.;
         for( ii=0; ii<=numbstd; ii++ ) rdata[ii] = sdata[ii];
      }
      else if ( 217433.0 < sdata[nsblstd+2] && sdata[nsblstd+2] < 217433.2 ) {
         lbit = 64;
         lthi = 0;
         inputcors.read((char*)&zdata,sizeof(zdata));
         for( ii=numbstd-1; ii<=numbthi; ii++ ) rdata[ii] = 0.;
         for( ii=0; ii<=numbstd; ii++ ) rdata[ii] = sdata[ii+1];
      }
      else if ( 217433.0 < sdata[nsblthi+1] && sdata[nsblthi+1] < 217433.2 ) {
         lbit = 32;
         lthi = 1;
         inputcors.read((char*)&hdata,sizeof(hdata));
         for( ii=0; ii<numbstd; ii++ ) rdata[ii] = sdata[ii];
         for( ii=numbstd; ii<numbthi; ii++ ) rdata[ii] = hdata[ii-numbstd];
      }
      else if ( 217433.0 < sdata[nsblthi+2] && sdata[nsblthi+2] < 217433.2 ) {
         lbit = 64;
         lthi = 1;
         inputcors.read((char*)&hdata,sizeof(hdata));
         inputcors.read((char*)&zdata,sizeof(zdata));
         for( ii=0; ii<numbstd; ii++ ) rdata[ii] = sdata[ii+1];
         for( ii=numbstd; ii<numbthi; ii++ ) rdata[ii] = hdata[ii-numbstd+1];
      }
      lbit = (lbit-32) / 32;
      nsblock = nsblstd + lthi*39;
      ndatpar = nsblock / 39;
      inputcors.clear();
      inputcors.close();
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// check 32 or 64 bit and standard or thinning simulation:
   if ( sdata[1] >= 211284. && sdata[1] <= 211286. ) { 
      ibit = 0;
      if ( sdata[313] >= 217432. && sdata[313] <= 217434. ) isubr=312;
      else if (sdata[274] >= 217432. && sdata[274] <= 217434.) isubr=273;
   }
   else if ( abs(sdata[1]) < 1.e-6 ) {
      ibit = 1;
      if ( sdata[2] >= 211284. && sdata[2] <= 211286. ) {
         if ( sdata[314] >= 217432. && sdata[314] <= 217434.) isubr=312;
         else if (sdata[275] >= 217432. && sdata[275] <= 217434.) isubr=273;
      }
   }
   else if ( sdata[3] >= 2.0202 && sdata[3] < 9.9999 ) {
      ibit = -1;
      if ( sdata[312] >= 217432. && sdata[312] <= 217434.) isubr=312;
      else if (sdata[273] >= 217432. && sdata[273] <= 217434.) isubr=273;
   }
   else {
      printf("     sdata[1] = %f case should not occur!\n",sdata[1]);
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - copy and use quantities of runh end evth subblocks:
   isubr = nsblock;
   drunh[1] = 211285.2812500000000000;
   for( i=2; i<=isubr; i++ ) drunh[i] = double(sdata[ibit+i]);
   if ( 22. < drunh[3]/10000. && drunh[3]/10000. < 99. )
      drunh[3] = 110101.; // set default 01.Jan.2011.
   devth[1] = 217433.0781250000000000;
   for( i=isubr+2; i<=isubr*2; i++ )
      devth[i-isubr] = double(sdata[ibit+i]);
   if ( devth[7] < 0. ) devth[7] = -devth[7]; // i.e. curved version.
   // check subblocks on multi-thinning identification:
   lemadd = 0;
   ltype8 = 0;
   for( i=1+2*isubr; i<5*isubr; i+=(isubr/39) ) {
      is = int(1.e-3*rdata[i]*1.0000001);
      if ( is == 8888 ) ltype8 = ltype8 + 1;
      if ( is == 85 || is == 86 || is == 95 || is == 96 ) lemadd++;
   }
   iday = int(drunh[3]) % 100;
   imon = (int(drunh[3]) % 10000) / 100;
   iyer = int(drunh[3]) / 10000 + 1900;
   if ( iyer < 1988 ) iyer = iyer + 100;
   if ( drunh[4] > 1.e6 ) drunh[4] = drunh[4] * 1.e-6;
   printf("\n            %s \n",chfile);
   printf("\n    runh         runnumb      date         versprog"
         "     nobslev      obslev1      obslev2       eslope\n");
   printf(" %12.5e    %.6d. %14.5e %12.5e %12.5e %12.5e %12.5e %10.4f\n",
      drunh[1],int(drunh[2]),drunh[3],drunh[4],drunh[5],drunh[6],
      drunh[7],devth[58]);
   printf("\n    engymin      engymax      flagEGS4     flagNKG"
          "      ecut.hadr    ecut.muon    ecut.elec    ecut.phot\n");
   printf("%14.7e %14.7e %8.4f %11.4f %14.5e %12.5e %12.5e %12.5e\n",
      drunh[17],drunh[18],drunh[19],drunh[20],
      drunh[21],drunh[22],drunh[23],drunh[24]);  
   printf("\n    evth         evtnumb      particle     energy"
      "       altit/gr     nrtarget     height/cm\n");
   printf(" %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %15.8e\n",
      devth[1],devth[2],devth[3],devth[4],devth[5],devth[6],devth[7]);
   printf("\n    theta/deg    phi/deg       hilow       nr.seeds"
      "       seed1\n");
   printf(" %12.7f %13.7f %10.4f %11.4f %17.8e\n",devth[11]*c180pi,
      devth[12]*c180pi,devth[155],devth[13],double(devth[14]));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - print type/kind of corsika simulation:
   printf("\n");
   if ( isubr == 312 ) {
      if ( nbit[1+ibit] == 32 ) printf("                   thinning"
           " simulation. lsubrec = %3d\n",isubr);
      else printf("            %2d bit thinning"
           " simulation. lsubrec = %3d\n",nbit[1+ibit],isubr);
   }
   else if ( isubr == 273 ) {
      if ( nbit[1+ibit] == 32 ) {
         if ( ltype8 == 0 ) printf("                   standard"
            " simulation. lsubrec = %3d\n",isubr);
         else printf("                   multi-thinning"
            " simu. lsubrec = %3d\n",isubr);
      }
      else printf("            %2d bit standard"
           " simulation. lsubrec = %3d\n",nbit[1+ibit],isubr);
   }
   printf("\n            %s \n",cskinput);
   iprm = int(devth[3]);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - runnr, parallel option (later check `.lst`) - - - - - - -
   coutext = fopen(cskinput, "w");
   fprintf(coutext,
      "RUNNR   %9d    version%9.6f    rundate %.2d.%.2d.%4d \n",      
         int(drunh[2]),drunh[4],iday,imon,iyer);
   if ( drunh[2] < 10000. && devth[13] >= 6. ) { 
      if ( devth[4] > 1.e7 ) {
         drunh[80] = 4.5432e5 * pow(1.95583,
                     int(2.*log10(drunh[17]*1.000001)-14.));
         fprintf(coutext, "* PARALLEL %8.0f. %10.0f.   1   F \n",
            drunh[80]/1.e3,drunh[80]);
      }
   }
   fprintf(coutext, "EVTNR   %9d \n",int(devth[2]));
   is = int(devth[13]);
   for( ii=1; ii<=is; ii++ ) 
      fprintf(coutext, "SEED    %9d %9d %9d \n",int(devth[13+3*(ii-1)+1]),
              int(devth[13+3*(ii-1)+2]),int(devth[13+3*(ii-1)+3]));
   if ( iprm == 4 ) printf("* should not occur: iprm == 4 (neutrino)\n");
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - primary, energy, angles - - - - - - - - - - - - - - - - -
   if ( ( 1 <= iprm && iprm <= 74 ) || ( 116 <= iprm && iprm <= 194 ) ) {
      fprintf(coutext, "PRMPAR  %9d                       %s \n",
         iprm,chptext[iprm]);
   }
   else if ( 201 <= iprm && iprm <= 5629 ) {
      fprintf(coutext, "PRMPAR  %9d                      %4d-%s \n",
         iprm,int(iprm/100),qnuclei[(iprm%100)]);
   }
   else if ( iprm%100 <= 0 || 29 < iprm%100 || drunh[3] < 100. ) {
      printf("             Warning: invalid runh[3] = %f\n\n",drunh[3]);
      goto continue29;
   }
   ienglow = int(log10(devth[59]+4.678e-34));
   ienghig = int(log10(devth[60]+4.678e-34));
   fprintf(coutext, "ERANGE  %9.4fe%+.2d   %9.4fe%+.2d \n",
      devth[59]/pow(10.,ienglow),ienglow,
      devth[60]/pow(10.,ienghig),ienghig);
   fprintf(coutext, "ESLOPE %13.3f \n",devth[58]);
   fprintf(coutext, "THETAP %13.3f %10.3f \n",devth[81],devth[82]);
   fprintf(coutext, "PHIP   %13.3f %10.3f \n",devth[83],devth[84]);
   iobs = int(devth[47]);
   isho = 1; // int(drunh[94]);
   if ( isho > 1 && strcmp(&chusnam[0],"simprod") == 0 ) {
      // user `simprod` special selection of shower numbers:.
      printf("         strcmp( chusnam[0],`simprod` == 0 ) \n"); 
      printf("        runh(94)=%f   runh(95)=%f   runh(96)=%f \n",
         drunh[94],drunh[95],drunh[96]);
      printf("NSHOW   %9d \n",isho);
      is = 4. * log10(devth[4]) - 19.;
      printf("         simprod is=%d \n",is);
      isho = nevent[is];
      // if ( lename-iusc == 2 && msimu > ncount[is] ) isho = nclast[is];
      // else isho = 1;
   }
   fprintf(coutext, "NSHOW   %9d \n",isho);
   iobs = 1;
   for( i=1; i<=iobs; i++ ) {
      gramms = thickgr(double(devth[47+i]));
      if ( 2838 <= int(devth[47+i]*1.00001234e-2) && 
                   int(devth[47+i]*1.00001234e-2) <= 2859 ) 
          gramms = thicksouth(double(devth[47+i]));
      if ( gramms < 900. ) gramms = gramms + 0.00123456;
      fprintf(coutext, "OBSLEV %12.2fe2 %16.7f g/cm^2 \n",
           devth[47+i]/100.,gramms);
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - cuts, prints, flags - - - - - - - - - - - - - - - - - - -
   fprintf(coutext, "ECUTS     %10.3f %10.3f %10.2e %10.2e \n",
      devth[61],devth[62],devth[63],devth[64]);
   fprintf(coutext, "ECTMAP         1.e11 \n");
   fprintf(coutext, "RADNKG %12.2fe2 \n",devth[147]/100.);
   fprintf(coutext, "MAXPRT          1 \n");
   fprintf(coutext, "HADFLG  %4d %4d %4d %4d %4d %4d \n",
      int(devth[65]),int(devth[66]),int(devth[67]),
      int(devth[68]),int(devth[69]),int(devth[70]));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - model quantities - - - - - - - - - - - - - - - - - - - - -
   fprintf(coutext, "ELMFLG    %s   %s \n",
      chflag[int(devth[73])],chflag[int(devth[74])]);
   if ( devth[76] == 0. )
      fprintf(coutext, "* HDPM            T         0 \n");
   if ( devth[76] == 1. ) {
      fprintf(coutext, "VENUS           T         0 \n");
      if ( devth[145] == 1. )
         fprintf(coutext, "VENSIG          T \n");
      else if ( devth[145] == 2. )
         fprintf(coutext, "NEXSIG          T \n");
      else if ( devth[145] == 3. )
         fprintf(coutext, "NEXSIG          T \n");
      else if ( devth[145] == 4. )
         fprintf(coutext, "EPOSIG          T \n");
      else if ( devth[145] >= 0 )
         fprintf(coutext, "VENSIG          T \n");
      fprintf(coutext, "* VENPAR   \'      \'         0. \n"); // VENPAR parcha parval
   }
   if ( devth[76] == 2. ) {
      fprintf(coutext, "SIBYLL          T         0 \n");
      if ( devth[140] >= 1. )
         fprintf(coutext, "SIBSIG          T \n");
      else printf("          check evth-elements 139,140 \n");
   }
   if ( devth[76] == 3. ) {
      fprintf(coutext, "QGSJET          T         0 \n");
      if ( devth[142] >= 1. )
         fprintf(coutext, "QGSSIG          T \n");
      else printf("       check evth-elements 141,142 \n");
   }
   if ( devth[76] == 4. ) {
      fprintf(coutext, "DPMJET          T         0 \n");
      if ( devth[144] >= 1. )
         fprintf(coutext, "DPJSIG          T \n");
      else printf("       check evth-elements 143,144 \n");
   }
   if ( devth[76] == 5. ) {
      fprintf(coutext, "NEXUS           T         0 \n");
      if ( devth[145] == 1. ) {
         fprintf(coutext, "VENSIG          T \n");
      }
      else if ( devth[145] == 2. ) {
         fprintf(coutext, "NEXSIG          T \n");
         // init input files for nexus use:
         fprintf(coutext, "NEXPAR fname inics nexus/nexus.inics\n");
         fprintf(coutext, "NEXPAR fname iniev nexus/nexus.iniev\n");
         fprintf(coutext, "NEXPAR fname initl nexus/nexus.initl\n");
         fprintf(coutext, "NEXPAR fname inirj nexus/nexus.inirj\n");
         // dummy out files for nexus (debug case):
         fprintf(coutext, "NEXPAR fname check none\n");
         fprintf(coutext, "NEXPAR fname histo none\n");
         fprintf(coutext, "NEXPAR fname data  none\n");
         fprintf(coutext, "NEXPAR fname copy  none\n");
      }
   }
   if ( devth[76] == 6. ) {
      fprintf(coutext, "EPOS            T         0 \n");
      // init input files for epos use:
      fprintf(coutext, "EPOPAR input epos/epos.param\n");
      fprintf(coutext, "EPOPAR fname inics epos/epos.inics\n");
      fprintf(coutext, "EPOPAR fname iniev epos/epos.iniev\n");
      fprintf(coutext, "EPOPAR fname initl epos/epos.initl\n");
      fprintf(coutext, "EPOPAR fname inirj epos/epos.inirj\n");
      fprintf(coutext, "EPOPAR fname inihy epos/epos.ini1b\n");
      // dummy out files for epos (debug case):
      fprintf(coutext, "EPOPAR fname check none\n");
      fprintf(coutext, "EPOPAR fname histo none\n");
      fprintf(coutext, "EPOPAR fname data  none\n");
      fprintf(coutext, "EPOPAR fname copy  none\n");
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - logicals, hilow, longi, magnet - - - - - - - - - - - - -
   tradius = 0.; // [cm]
   if ( devth[4] >= 1.0e+08 ) tradius = 1.e0;
   if ( devth[4] >= 3.2e+08 ) tradius = 3.e0;
   if ( devth[4] >= 1.0e+09 ) tradius = 10.e0;
   if ( devth[4] >= 3.2e+09 ) tradius = 30.e0;
   if ( devth[4] >= 1.0e+10 ) tradius = 100.e0;
   if ( devth[4] >= 3.2e+10 ) tradius = 120.e0;
   if ( devth[4] >= 1.0e+11 ) tradius = 140.e0; // [m]
   fprintf(coutext, "MUMULT          T \n");
   if ( int(devth[94]) == 1 )
      fprintf(coutext, "MUADDI    %s \n",chflag[int(devth[94])]);
   if ( lemadd > 0 ) fprintf(coutext, "EMADDI          T \n");
   if ( devth[152] > 0.9876 ) 
      fprintf(coutext, "CORECUT %12.3f  cm \n",devth[152]);
   if ( devth[95] < 1. ) devth[95] = 1.;
   fprintf(coutext, "STEPFC  %12.3f \n",devth[95]);
   if ( devth[93] != 0. ) fprintf(coutext, "ARRANG  %12.2f \n",devth[93]);
   gramms = thickgr(double(devth[48]));
   if ( 2838 <= int(devth[48]*1.00001234e-2) && 
                int(devth[48]*1.00001234e-2) <= 2859 ) 
       gramms = thicksouth(double(devth[7]));
   fprintf(coutext, "FIXHEI %15.5fe2 %7d %15.7f g/cm^2 \n",
      devth[7]/100.,int(devth[6]),gramms);
   fprintf(coutext, "* ERANGE %17.7e %15.7e \n",devth[4],devth[4]);
   fprintf(coutext, "* THETAP %17.7f %15.7f \n", 
      devth[11]*c180pi,devth[11]*c180pi);
   fprintf(coutext, "* PHIP   %17.7f %15.7f \n",
      devth[12]*c180pi,devth[12]*c180pi);
   fprintf(coutext, "* OUTFILE DAT%.6d.firstint\n",int(drunh[2]));
   if ( devth[75] >= 3. ) {
      hilow = 200.;
      if ( devth[155] >= 80. ) hilow = devth[155];
      fprintf(coutext, "* low energy model fluka used. \n");
   }
   else if ( devth[75] == 2. ) {
      hilow = 80.;
      if ( devth[155] >= 80. ) hilow = devth[155];
      fprintf(coutext, "URQMD             T \n");
   }
   else if ( devth[75] <= 1. ) {
      hilow = 80.;
      if ( devth[155] >= 80. ) hilow = devth[155];
   }
   fprintf(coutext, "HILOW  %12.2f \n",hilow);
   gramlong = 5.0; // possibly check `.long` format(30x,i5,19x,f5.0).
   fprintf(coutext, "LONGI           T %8.1f      T      T \n",gramlong);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - check quantities of magnetic field - - - - - - - - - -
   if ( ( 17. < devth[71] && devth[71] <  21. ) &&
       ( -16. < devth[72] && devth[72] < -12. ) ) {
      fprintf(coutext, "MAGNET  %11.2f %11.2f     Auger \n",
         devth[71],devth[72]);
   }
   else if ( ( 18. < devth[71] && devth[71] < 22. ) &&
             ( 41. < devth[72] && devth[72] < 45. ) ) {
      fprintf(coutext, "MAGNET  %11.2f %11.2f     Karlsruhe \n",
         devth[71],devth[72]);
   }
   else if ( fabs(devth[71]) < 1.e-2 && fabs(devth[72]) < 1.e-2 ) {
      fprintf(coutext, "MAGNET  %12.3e %12.3e     NoMag \n",
         devth[71],devth[72]);
   }
   else if ( ( 26. < devth[71] && devth[71] <  29. ) &&
            ( -46. < devth[72] && devth[72] < -50. ) ) {
      fprintf(coutext, "MAGNET  %11.2f %11.2f     SKA (West-Australia)\n",
         devth[71],devth[72]);
   }
   else if ( ( 14. < devth[71] && devth[71] <  18. ) &&
            ( -55. < devth[72] && devth[72] < -51. ) ) {
      fprintf(coutext, "MAGNET  %11.2f %11.2f      South Pole \n",
         devth[71],devth[72]);
   }
   else {
      fprintf(coutext, "MAGNET   %11.2f %11.2f \n",devth[71],devth[72]);
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - thinning quantities - - - - - - - - - - - - - - - - - - -
   if ( devth[150] > 0. ) {
      devth[149] = 1.00001 * devth[149];
      fprintf(coutext, "THIN      %11.4e %11.4e %11.4e \n",
         devth[149],devth[151],devth[152]);
      fprintf(coutext, "THINH  %10.0f. %10.0f \n",1.,100.);
         // devth(149)/devth(148),devth(151)/devth(150));
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - cherenkov quantities  - - - - - - - - - - - - - - - - - -
   if ( devth[77] > 0. ) {
      if ( devth[85] > 0. ) {
         fprintf(coutext, "CERSIZ   %8.0f. \n",devth[85]);
         fprintf(coutext, "CERFIL   %8d \n",int(devth[92])); 
         // telescope coord. x,y,z, and radius r (IACT or not):
         fprintf(coutext, "* TELESCOPE  -3.17411e+06  1.50956e+06"
               "  21054.8    15.e2   Heat    `IACT`\n");
      }
      if ( devth[96] > 0. )
         fprintf(coutext, "CWAVLG  %8.0f. %8.0f. \n",devth[96],devth[97]);
      if ( devth[86] > 0. || devth[87] > 0. )
         fprintf(coutext, "CERARY  %4d %4d %10.2e %10.2e %10.2e %10.2e "
            " `not IACT`\n",int(devth[86]),int(devth[87]),devth[88],
            devth[89],devth[90],devth[91]);
      if ( devth[98] > 1. ) 
         fprintf(coutext, "CSCAT    %8d  %7.0f.e2%7.0f.e2 \n",
            int(devth[98]),drunh[248]/100.,drunh[249]/100.); 
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - core positions  - - - - - - - - - - - - - - - - - - - - -
   if ( devth[98] > 0. ) {
      mpos = int(devth[98]);
      for( is=1; is<=mpos; is++ )
         printf("COREPOS    %13.6e %13.6e        (pos.%.2d) \n",
            devth[98+is]+11.,devth[118+is]+12.,is); 
         fprintf(coutext, "COREPOS    %13.6e %13.6e        (pos.%.2d) \n",
            devth[98+is]+11.,devth[118+is]+12.,is);
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - multithin parameters  - - - - - - - - - - - - - - - - - -
   if ( devth[177] > 0 ) {
      multhi = int(devth[177]);
      fprintf(coutext, "MTHINR  %13.3f  cm \n",devth[152]);
      for( is=1; is<=multhi; is++ ) {
         fprintf(coutext, "MTHINH   %9.1e %9.1e %9.1e %9.1e \n",
            devth[177+is],devth[183+is],devth[189+is]/devth[177+is],
            devth[195+is]/devth[183+is]);
         fprintf(coutext, "MSEED    %9d %9d %9d \n",
            int(devth[201+is]),int(devth[207+is]),int(devth[213+is])); 
      }
   }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - direct, host, user  - - - - - - - - - - - - - - - - - - -
   if ( 2838 <= int(devth[48]*1.00001234e-2) &&
                int(devth[48]*1.00001234e-2) <= 2859 ) 
      fprintf(coutext, "ATMOD         13      South Pole October\n");
   else if ( strcmp(&chusnam[0],"maximo") == 0 ) {
      fprintf(coutext, "ATMOD       0 \n");
      fprintf(coutext, "ATMA   %14.7e %14.7e %14.7e %14.7e %12.5e\n",
         aatmax[1],aatmax[2],aatmax[3],aatmax[4],aatmax[5]);
      fprintf(coutext, "ATMB   %14.7e %14.7e %14.7e %14.7e %12.5e\n",
         batmax[1],batmax[2],batmax[3],batmax[4],batmax[5]);
      fprintf(coutext, "ATMC   %14.7e %14.7e %14.7e %14.7e %12.5e\n",
         catmax[1],catmax[2],catmax[3],catmax[4],catmax[5]);
      fprintf(coutext, "ATMLAY %8.3fe5 %8.3fe5 %8.3fe5 %8.3fe5 %8.3fe5\n", 
         atmlay[1]*1.e-5,atmlay[2]*1.e-5,atmlay[3]*1.e-5,
         atmlay[4]*1.e-5,atmlay[5]*1.e-5);
      fprintf(coutext, "DIRECT /cr/data02/joe/corsika-run/\n");
   }
   else {
      fprintf(coutext, "* ATMOD       0            monthly-atmosphere-04"
             "    (optional)\n");
      fprintf(coutext, "* ATMA     -1.33894966E+02 -4.70898991E+01 "
             " 8.32783065E-01 -2.10543673E-04\n");
      fprintf(coutext, "* ATMB      1.17413853E+03  1.21573550E+03 "
             " 1.42569186E+03  5.03931697E+02\n");
      fprintf(coutext, "* ATMC      9.67164959E+05  7.68155103E+05 "
             " 6.21852254E+05  7.85601620E+05\n");
      fprintf(coutext, "* ATMLAY    9.90000000E+05  1.21000000E+06 "
             " 3.80000000E+06  1.00000000E+07\n");
      /* fprintf(coutext, "* ATMOSPHERE    20    F     (needs file "
             "`atmprof20.dat`)\n"); */
   }
   if ( devth[60] > 3.e6 && devth[13] >= 6. ) {
      fprintf(coutext, "* OUTFILE csk%.6d/DAT%.6d.firstint\n",
                             int(drunh[2]),int(drunh[2]));
      fprintf(coutext, "DIRECT csk%.6d/\n",int(drunh[2]));
      fprintf(coutext, "HOST   uc1n996\n");
      fprintf(coutext, "USER   jl5949\n");
   }
   else {
      fprintf(coutext, "DIRECT ./\n");
      fprintf(coutext, "HOST   iklx288\n");
      fprintf(coutext, "USER   joe\n");
   }
   fprintf(coutext, "EXIT\n");

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - end-of program cors2input - - - - - - - - - - - - - - - - -
 continue29:
   printf("  \n");
   return 0;

}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = pamafill= = = = = = = = = = = = = = = = =
//                    fill particle masses up to 59_Ni.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void pamafill()
{
   double az,eb;
   int ia,in,ip,is;

// - - - - - - calculate particle (nuclei) masses:   
   for( ia=0; ia<=76; ia++ )
      { pmass[ia] = parmas[ia]; psign[ia] = chargs[ia]; }
   for( ia=77; ia<=6000; ia++ )
      { pmass[ia] = 0.; psign[ia] = 0.; }
   for( ia=1; ia<60; ia++ ) {
      for( ip=1; ip<=ia; ip++ ) {
         in = ia - ip;
         is = ia * 100 + ip;
         // - - without binding energy effects.
         pmass[is] = parmas[13] * in + parmas[14] * ip;
         // - - nuclei are assumed to be fully ionized.
         psign[is] = ip;
         // - - applying binding energy of nuclei.
         az = ia;
         eb = 14.1 * az - 0.595 * ip*ip * pow(az,-1./3.)
            - 13. * pow(az,2./3.) - 19. * (ip-in)*(ip-in) / az;
         if ( ip%2 == 0 && in%2 == 0 )
            { eb = eb + 33.5 * pow(az,-0.75); }
         else
            if ( ip%2 == 1 && in%2 == 1 )
            { eb = eb - 33.5 * pow(az,-0.75); }
         if ( eb > 0. ) eb=1.e-3*eb;
         else eb=0.; 
         pmass[is] = pmass[is]-eb;
      }
   }
   // masses of multi-neutron clusters.
   for( in=1; in<60; in++ ) {
      is = in * 100;
      pmass[is] = parmas[13] * in;
      psign[is] = 0.;
   }
}

 
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = heightcm  = = = = = = = = = = = = = = = =
//     calculate height (cm) a.s.l. for a given thickness (gramms/cm**2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double heightcm(double g)
{
   double height;
   if (g > 631.1)
      { height = CATM[1] * log( BATM[1] / (g-AATM[1]) ); }
   else if (g > 271.7)
      { height = CATM[2] * log( BATM[2] / (g-AATM[2]) ); }
   else if (g > 3.0395)
      { height = CATM[3] * log( BATM[3] / (g-AATM[3]) ); }
   else if (g > 1.28292e-3)
      { height = CATM[4] * log( BATM[4] / (g-AATM[4]) ); }
   else
      { height = (AATM[5] - g) / CATM[5]; }
   return height;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
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


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = thicksouth  = = = = = = = = = = = = = = =
//     thicksouth (gramms/cm**2) of atmosphere depending on height (cm)
//     for the South Pole atmosphere of October 01 of the year.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double thicksouth(double h)
{
      double thicks,
      spaatm[6]={0.,-142.801,-70.1538,1.14855,9.10269e-4,1.52236e-3},
      spbatm[6]={0.,1177.19,1125.11,1304.77,433.823,1.},
      spcatm[6]={0.,861745.,765925.,581351.,775155.,7.4095699},
      sphlay[6]={0.,-5.e5, 2.66667e5, 5.33333e5, 8.e5, 9.8765e6};
      if ( h < sphlay[2] ) {
         thicks = spaatm[1] + spbatm[1] * exp( -h/spcatm[1] ); }
      else if ( h < sphlay[3] ) {
         thicks = spaatm[2] + spbatm[2] * exp( -h/spcatm[2] ); }
      else if ( h < sphlay[4] ) {
         thicks = spaatm[3] + spbatm[3] * exp( -h/spcatm[3] ); }
      else if ( h < sphlay[5] ) {
         thicks = spaatm[4] + spbatm[4] * exp( -h/spcatm[4] ); }
      else {
         thicks = spaatm[5] -h/spcatm[5]; }
      return thicks;
}

