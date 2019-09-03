//========================================================================
// 
//  r e a d c o r s i k a . c p p
//  =============================
//    This program is able to read all binary corsika particle data files 
//    of "standard", "multithinning", and "thinning" simulations;
//    it writes a protocol file to standard output containing infos of 
//        primary, energy, runnumb, simnrsh, #hadrs, #muons, #gammas,
//            #elecs, #nkgelecs, obslev, theta, phi, h1km, h1gr.
// - - - - - - - - - - - - - - CompLinkRun - - - - - - - - - - - - - - - -
// CompLink:
//           g++ -O0 -fbounds-check readcorsika.cpp -o readcorsika -lm
// RunProgr:
//           ./readcorsika DAT370375 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// input-file: current particle data file 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         RUNH = 211285.2812500000000000;
//         EVTH = 217433.0781250000000000;
//         LONG =  52815.2968750000000000;
//         EVTE =   3397.3918457031250000;
//         RUNE =   3301.3325195312500000;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         standard simulation 22932 == 3.2134576383896705E-41
//         thinning simulation 26208 == 3.6725230153024806E-41
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//   naming conventions of corsika particles: 
//    1   gamma           24   Omega           64   K* -
//    2   positron        25   anti neutron    65   anti K* 0
//    3   electron        26   anti Lambda     66   electron neutrino
//    4                   27   anti Sigma -    67   electron anti neutrino
//    5   muon +          28   anti Sigma 0    68   muon neutrino
//    6   muon -          29   anti Sigma +    69   muon anti neutrino
//    7   pion 0          30   anti Xi 0       71   eta-> 2*gam        
//    8   pion +          31   anti Xi +       72   eta-> 3*pi0
//    9   pion -          32   anti Omega      73   eta-> pi+ + pi- + pi0
//   10   Kaon 0 long     50   omega           74   eta-> pi+ + pi- + gam
//   11   Kaon            51   rho 0           201   Deuteron  
//   12   Kaon -          52   rho +           301   Tritium
//   13   neutron         53   rho -           402   alpha
//   14   proton          54   Delta ++       1206   Carbon
//   15   anti proton     55   Delta +        1407   Nitrogen
//   16   Kaon 0 short    56   Delta 0        1608   Oxygen  
//   17   eta (71..74)    57   Delta -        2713   Aluminium  
//   18   Lambda          58   anti Delta --  2814   Silicon 
//   19   Sigma +         59   anti Delta -   3216   Sulfur
//   20   Sigma 0         60   anti Delta 0   5626   Iron    
//   21   Sigma -         61   anti Delta +   9900   Cherenkov photons   
//   22   Xi 0            62   K* 0       
//   23   Xi -            63   K* +     
//  116   D 0            131   tau +           150   anti Xi c -
//  117   D +            132   tau -           151   anti Xi c 0
//  118   anti D -       133   tau neutrino    152   anti Sigma c --
//  119   anti D 0       134   anti tau neutr  153   anti Sigma c -
//  120   D s +          137   Lambda c +      154   anti Sigma c 0
//  121   anti D s -     138   Xi c +          155   anti Xi c prime -
//  122   eta c          139   Xi c 0          156   anti Xi c prime 0
//  123   D*0            140   Sigma c ++      157   anti Omega c 0
//  124   D*+            141   Sigma c +       161   Sigma c * ++
//  125   anti D*-       142   Sigma c 0       162   Sigma c * +
//  126   anti D*0       143   Xi c prime +    163   Sigma c * 0
//  127   D* s +         144   Xi c prime 0    171   anti Sigma c * --
//  128   anti D* s -    145   Omega c 0       172   anti Sigma c * -
//  130   J/psi          149   anti Lambda c-  173   anti Sigma c * 0
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

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = function definitions  = = = = = = = = = = = = =

void analyzeEvte(); // analyze subblock EVTE.

void analyzeEvth(); // analyze subblock EVTH.

void analyzeLong(); // analyze subblock LONG.

void analyzePart(); // analyze particle data subblock.

void analyzeRune(); // analyze subblock RUNE.

void analyzeRunh(); // analyze subblock RUNH.
 
double heightcm(double g); // calculate height (cm) for a given thickness.

void pamafill(); // fill particle masses up to 59_Ni.

void pamasscalc(); // calculate particle masses up to 59_Ni.

void printable(int m); // formatted output table of shower quantities.

void printpartic(int m); // print particles of the first subblock.

double thickgr(double h); // thickness of atmosphere depending on height.

double thicksouth(double h); // thickness of atmosphere at South Pole.

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = global constants  = = = = = = = = = = = = =

   const int nrecstd = 22940;           // "standard corsika" record size,
                                        //        pdata[0] = 6.08270e-311.
                                        //    or  pdata[0] = 3.213456e-41.
   const int nrecthi = 26216;           // "thinning corsika" record size,
                                        //        pdata[0] = 6.95166e-311.
                                        //    or  pdata[0] = 3.672523e-41.
   const int ndim = 20000000;
   const int nsubblo = 21;
   const int nsblstd = 273;
   const int nsblthi = 312;
   const int nmaxrec = 2180000;              // valid up to 50 GBytes.
   const int numbstd = nrecstd / 4;          // =  5735.
   const int numbthi = nrecthi / 4;          // =  6554.
   int ndatpar, nreclen, nsblock, lbit,
       nxpix=750, nypix=650, nprimry, lobslev, jobslev, ifil, iruntyp=2;

   float pdata[numbthi]; // to read a single corsika record thinned data.
   float rdata[numbthi]; // to keep first corsika record thinned data.
   float sdata[numbstd]; // to read a single corsika record standard data.
   float hdata[819]; // to read end of first corsika record thinned data.
   float zdata[2]; // to read record length information on 64bit machines. 

   double qrunh[275], qevth[275], qevte[275], qdata[101], obslev, obsgrm;

   const double velight = 29.9792458; // velocity of light in cm/nsec.

   const double AATM[6] = { 0.e0,       // param for US standard atmosph.
                   -186.5562e0, -94.919e0, 0.61289e0, 0.e0, 1.128292e-2 };
   const double BATM[6] = { 0.e0,       // param for US standard atmosph.
                1222.6562e0, 1144.9069e0, 1305.5948e0, 540.1778e0, 0.e0 };
   const double CATM[6] = { 0.e0,       // param for US standard atmosph.
              994186.38e0, 878153.55e0, 636143.04e0, 772170.16e0, 1.e-9 };
 
   const double parmas[202] = { 0.0e0,          // particle masses in GeV.
         0.e0     ,5.10998902e-4,5.10998902e-4, 0.e0     ,0.105658357e0,
         0.105658357e0,0.1349766e0,0.13957018e0,0.13957018e0,0.497672e0,
         0.493677e0, 0.493677e0, 0.93956533e0,0.93827200e0,0.93827200e0,
         0.497672e0 , 0.54730e0  , 1.115683e0 , 1.18937e0 , 1.192642e0 ,
         1.197449e0 , 1.31483e0  , 1.32131e0  , 1.67245e0 ,0.93956533e0,
         1.115683e0 , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31483e0  ,
         1.32131e0  , 1.67245e0  , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.78257e0  ,
         0.7690e0   , 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
         1.2331e0   , 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   ,
         1.2349e0   , 0.89610e0  , 0.89166e0  , 0.89166e0 , 0.89610e0  ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.54730e0  , 0.54730e0  , 0.54730e0  , 0.54730e0 , 0.105658e0 ,
         0.105658e0 , 0.0e0      , 0.0e0      , 0.0e0     , 0.0e0      ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         1864.5e0   , 1869.3e0   , 1869.3e0   , 1864.5e0  , 1968.6e0   ,
         1968.5e0   , 2979.7e0   , 2006.7e0   , 2010.0e0  , 2010.0e0   ,
         2006.7e0   , 2112.4e0   , 2112.4e0   , 3510.51e0 , 3096.87e0  ,
         1776.99e0  , 1776.99e0  , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 2284.9e0   , 2466.3e0   , 2471.8e0  , 2452.6e0   ,
         2451.3e0   , 2452.2e0   , 2574.1e0   , 2578.8e0  , 2697.5e0   ,
         0.e0       , 0.e0       , 0.e0       , 2284.9e0  , 2466.3e0   ,
         2471.8e0   , 2452.6e0   , 2451.3e0   , 2452.2e0  , 2574.1e0   ,
         2578.8e0   , 2697.5e0   , 0.e0       , 0.e0      , 0.e0       ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         2519.4e0   , 2515.9e0   , 2517.5e0   , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0       ,
         0.e0       , 0.e0       , 0.e0       , 0.e0      , 0.e0, 0.e0 };

   const double chargs[202] = { 0.,                 // particle charges.
           0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  ,
          +1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  ,+1.e0, 0.  ,
          -1.e0, 0.  ,-1.e0,-1.e0, 0.  , 0.  ,-1.e0, 0.  ,+1.e0, 0.  ,
          +1.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  ,+1.e0,-1.e0,+2.e0,+1.e0, 0.  ,-1.e0,-2.e0,-1.e0, 0.  ,
          +1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,
          -1.e0, 0.  , 0.  ,+1.e0,-1.e0, 0.  ,+1.e0,-1.e0, 0.  , 0.  ,
          +1.e0,-1.e0, 0.  , 0.  , 0.  , 0.  ,+1.e0,+1.e0, 0.  ,+2.e0,
          +1.e0, 0.  ,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,-1.e0,-1.e0,
           0.  ,-2.e0,-1.e0, 0.  ,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  ,
          +2.e0,+1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
          -2.e0,-1.e0, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
           0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  };

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = global quantities = = = = = = = = = = = = =
  
   char corsfile[21][240];
  
   int mprimary, ibytes, mpinull, mpiplus, mpiminus, mcher=0, mbsize=0,
       nshowf[21], nfils, nshow, nshof, ishow, ishof, isho, irec, itext,
       idpa, iobs, igen, icod, mutr, ii, la, lb, lz, lp, lthi=0;
  
   double epart, eplog, equad, busize; // extra memory for quantities.
  
   double pmass[6001], psign[6001], prest[6001], rx, ry;
  
   double slopest, engylow, enghigh, thetlow, thehigh, obslev1;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = = = main= = = = = = = = = = = = = = = = =
//          
//=== II. - read file name as argument - - - - -
//=== III.- loop over all file names read in
//==== III.a.- - test next corsika file
//==== III.b.- - loop over all data records
//==== III.b.1.- - - read a single record from corsika file
//==== III.b.2.- - - analyze subblocks of the record
//=== IV. - after all showers and particles
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int main(int argcnt, char *argvec[])
{
   char chtext[240];
   time_t tseconds;
   tm *datetime;

// - - - - - copy execution argument (i.e. file name) - - - - - - - - - -
   printf(" _ _ _ _ _ _ readcorsika.cpp _ _ _ _ _ _"
          " _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n");
   if ( argcnt > 1 ) {
      strcpy(chtext,argvec[1]);
      itext = strlen(chtext);
      chtext[itext] = '\0'; 
      strcpy(corsfile[1],chtext);
      nshowf[ifil] = 1;
      nshow = 1;
      nfils = 1;
      isho = 1;
   }
   else {
      printf("\n             no particle data to process.\n");
      printf("             usage: ./readcorsika DAT123456\n");
      goto donothing;
   } 

// - - - - - print title and date of analysis - - - - - - - - - - - - - -
   tseconds = time(&tseconds);
   datetime = localtime(&tseconds);
   cout << "   " << endl;
   printf("        _date_of_analysis_ ");
   printf(" %.2d. %.2d. %.4d.  %.2d.%.2d%.2d hms  (%ld sec) \n",
      datetime->tm_mday, datetime->tm_mon+1, datetime->tm_year+1900,
      datetime->tm_hour, datetime->tm_min, int(tseconds)%60, tseconds); 

// - - - - - fill particle masses up to 5928_Ni; valid <= 5656 - - - - - -
   pamafill();

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== III. - loop over all file names read in- - - - - - - - - - - - - - -
   nshof = 0;
   ishow = 0;
   for( ifil=1; ifil<=nfils; ifil++ )
   {
   nshof = nshof + nshowf[ifil]; // total showers till end of current file
   ishof = 0;                    // counting showers of current file.
   cout << "   " << endl;
   printf("        %s       %d shower\n",corsfile[ifil],isho);

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.1.- - read/check first record from current corsika file - - -
   ifstream inputcors(corsfile[ifil],ios::in); // allocate file.
   if ( inputcors.peek() == EOF ) {
      printf("         file ==> does not exist at this path! Stop!\n");
      inputcors.clear();
      inputcors.close();
      break;
   }
   else {
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
      // check particle data as of cherenkov bunches incl. bunch size:
      mcher = int(sdata[nsblock+85]) * (2*int(sdata[nsblock+92])-1);
      // mcher > 0 for separate cherenkov bunches data file, 
      // mcher < 0 as of cherenkov bunches in particle data file.
   }

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b - - loop over all data records of current file  - - - - - - -
   for( irec=1; irec<nmaxrec; irec++ ) {
   // if ( irec%10 == 1 ) printf("        irec= %d \n",irec);

// _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.1.- - - get record of corsika particle data file: - - - - - -
   if ( irec > 1 ) {
   // - - - - - - - - read data irec>1 of known corsika record length:
      if ( lthi == 1 ) { // - - - - - - - long record (thinning).
         inputcors.read((char*)&pdata,sizeof(pdata));
         if ( pdata[0] < 1.21345763838967e-41 ) { // EOF found.
            break;
         }
         // - - - - - - - not yet at end of file:
         if ( lbit == 1 ) inputcors.read((char*)&zdata,sizeof(zdata));
         pdata[0] = 0.21345763838967e-41;
         if ( lbit == 1 ) {
            for( ii=1; ii<numbthi; ii++ ) pdata[ii] = pdata[ii+lbit];
            pdata[numbthi] = zdata[0];
         }
         for( ii=1; ii<=numbstd; ii++ ) sdata[ii] = pdata[ii];
      }
      else if ( lthi == 0 ) { // - - - - - - - short record (standard).
         inputcors.read((char*)&sdata,sizeof(sdata));
         if ( sdata[0] < 1.21345763838967e-41 ) { // EOF found.
            break;
         }
         // - - - - - - - not yet at end of file:
         if ( lbit == 1 ) inputcors.read((char*)&zdata,sizeof(zdata));
         sdata[0] = 0.21345763838967e-41;
         for( ii=numbstd; ii<=numbthi; ii++ ) pdata[ii] = 0.;
         for( ii=1; ii<=numbstd; ii++ ) pdata[ii] = sdata[ii+lbit];
      }
   }
   // - - - - - - - - take data of first record from read/check:
   else if ( irec == 1 ) {
      if ( iruntyp > 1 || ifil == 1 ) {
         if ( mcher != 0 ) printf("                        "
            "bunchsize= %d  (Cherenkov)\n",int(fabs(1.*mcher))); 
         printf("                          nsblock= %d ",nsblock);
         if ( lbit == 1 ) printf("     64bit file");
         printf(" \n");
      }
      for( ii=0; ii<numbthi ; ii++ ) pdata[ii] = rdata[ii];
      for( lz=1; lz<=nsblstd; lz++ ) qrunh[lz] = rdata[lz];
      for( lz=1; lz<=nsblstd; lz++ ) qevth[lz] = rdata[lz+nsblock];
      // - - - - - - - - optional test print to verify contents:
      if ( lbit > 2 ) {
         printf("        rdata[0-1]   %18.15e  %18.15e \n",rdata[0],rdata[1]);
         printf("        rdata[2-6]   %f  %f  %f  %f  %f \n",rdata[2],rdata[3],
             rdata[4],rdata[5],rdata[6]);
         printf("       rdata[7-11]   %f  %f  %f  %f  %f \n",rdata[7],rdata[8],
             rdata[9],rdata[10],rdata[11]);
         if ( nsblock == 273 ) {
            printf("  rdata[5725-5729]   %f  %f  %f  %f  %f \n",rdata[5725],
             rdata[5726],rdata[5727],rdata[5728],rdata[5729]);
            printf("  rdata[5730-5734]   %f  %f  %f  %f  %f \n",rdata[5730],
             rdata[5731],rdata[5732],rdata[5733],rdata[5734]);
            printf("  rdata[5735-5736]   %18.15e  %18.15e \n",
             rdata[5735],rdata[5736]);
         }
      }
      // - - - - - - - - optional print out of quantities of first record:
      if ( lbit < 2 ) {
         if ( iruntyp > 1 ) { 
            cout << "   " << endl;
            lp = 1;
            printf("            RUNH          runnr        date_of_run"
               "  vers_of_prog    nobslev       obslevcm\n");  
            printf(
               " subl%.2d:%14.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
               int((lp+nsblock)/nsblock),rdata[lp],rdata[lp+1],
               rdata[lp+2],rdata[lp+3],rdata[lp+4],rdata[lp+5]);
            lp = nsblock+1; 
            printf("            EVTH          event_nr     particle_id"
               "    energy      start_altit.    target_1st   h1st_intact\n");
            printf(
               " subl%.2d:%14.6e %13.6e %13.6e %13.6e %13.6e %13.6e"
               " %13.6e\n",int((lp+nsblock)/nsblock),rdata[lp],rdata[lp+1],
               rdata[lp+2],rdata[lp+3],rdata[lp+4],rdata[lp+5],rdata[lp+6]);
            lp = 2*nsblock+1; 
            printf("            particle      px_moment.    py_moment."
               "    pz_moment.    x_coord       y_coord        time nsec\n");
            printf(
               " subl%.2d:%14.6e %13.6e %13.6e %13.6e %13.6e %13.6e"
               " %13.6e\n",int((lp+nsblock)/nsblock),rdata[lp],rdata[lp+1],
               rdata[lp+2],rdata[lp+3],rdata[lp+4],rdata[lp+5],rdata[lp+6]);
         }
         cout << "   " << endl;
      }
   }

// _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.2.- - - analyze subblocks of the record - - - - - - - - - - -
      for( lp=1; lp<nsblock*21; lp+=nsblock )
      {

//========== RUNH subblock found =========================================
      if ( 211285.2 < pdata[lp] && pdata[lp] < 211285.3 ) {
         if ( ifil < 0 ) cout << "        RUNH subblock found. " << endl;
         analyzeRunh(); 
      }

//========== EVTH subblock found =========================================
      else if ( 217433.0 < pdata[lp] && pdata[lp] < 217433.1 ) {
         if ( ifil < 0 ) cout << "        EVTH subblock found. " << endl;      
         analyzeEvth();
      }

//========== LONG subblock found =========================================
      else if ( 52815.2 < pdata[lp] && pdata[lp] < 52815.3 ) {
         if ( ifil < 0 ) cout << "        LONG subblock ignored." << endl;
         analyzeLong();
      }

//========== EVTE subblock found =========================================
      else if ( 3397.3 < pdata[lp] && pdata[lp] < 3397.4 ) {
         if ( ifil < 0 ) cout << "        EVTE subblock found. " << endl;
         analyzeEvte();
         if ( ishow == nshow ) break; // all showers analyzed. 
      }  

//========== RUNE subblock found =========================================
      else if ( 3301.3 < pdata[lp] && pdata[lp] < 3301.4 ) {
         if ( ifil < 0 ) cout << "        RUNE subblock found. " << endl;
         analyzeRune();
      }

//========== particle data subblock ======================================
      else { // if ( 0. < pdata[lp] ) {
             // check particles of option EMADDI etc.
         analyzePart();
      }

      } // end_of_loop of analyzing subblocks - - - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.c.- - - continue with next data record or close file. 
 
      if ( ishow == nshow ) break; // all showers analyzed.

    } // endfor_loop over all data records - - - - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

    inputcors.close(); // close corsika file also after break at EOF;
 
    if ( ishow == nshow ) break; // all showers analyzed.

   } // endfor_loop over all file names read in  - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== IV. - after all showers and particles- - - - - - - - - - - - - - - -
   cout << "   " << endl;
   cout << "        _end_of_data"; 
   if ( ishow == 1 )
      cout << "       1 shower       nrec= " << irec;                
   else if ( ishow > 1 )
      cout << "           " << ishow << " showers.";
   else 
      cout << "           WARNING, case should not occur! ";
   cout << endl;
 donothing:
   cout << endl;
   
   return 0;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeEvte = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock event end
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeEvte()
{
   ishow++;
   if ( ishow <= nshow ) {
      // - - test-print of some quantities of evte subblock - -
      cout.precision(5);
      cout.setf(ios::scientific,ios::fixed);
   }
   if (pdata[lp+4]+pdata[lp+5] > 0.) {
      // - - fill electron numbers also at NKG simulations - -
      if ( pdata[lp+184] < 1. ) {
         la = 1;
         pdata[lp+184] = pdata[lp+184-la];
         while( pdata[lp+184] < 1. )
            { la++; pdata[lp+184] = pdata[lp+184-la]; } 
      }
      if ( pdata[lp+3] < 1. ) pdata[lp+3] = pdata[lp+184];
      if ( pdata[lp+194] < 0.25 ) {
         la = 1;
         pdata[lp+194] = pdata[lp+194-la];
         while( pdata[lp+194] < 0.25)
            { la++; pdata[lp+194] = pdata[lp+194-la]; } 
      }
   }
//- - - - - - - - - - number of hadrons at obslev - - - - - - - - - - - -
   qdata[5] = pdata[lp+4]; 
//- - - - - - - - - - number of muons at obslev - - - - - - - - - - - - -
   qdata[6] = pdata[lp+5]; 
//- - - - - - - - - - number of photons at obslev - - - - - - - - - - - -
   qdata[7] = pdata[lp+2]; 
//- - - - - - - - - - number of electrons at obslev - - - - - - - - - - -
   qdata[8] = pdata[lp+3]; 
//- - - - - - - - - - NKG number of electrons at obslev - - - - - - - - -
   qdata[9] = pdata[lp+184]; 
//- - - - - - - - - - shower age at obslev- - - - - - - - - - - - - - - -
   qdata[10] = pdata[lp+194]; 
   printable(ishow);           // print shower information to file.
   return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeEvth = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock event header
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeEvth()
{
//- - - - - - - - - - number of shower in current file - - - - - - - - - -
   ishof++;
//- - - - - - - - - - set counters to zero - - - - - - - - - - - - - - - -
   mutr = 0; 
   mpinull = 0;
   mpiplus = 0;
   mpiminus = 0;
//- - - - - - - - - - number of observation level - - - - - - - - - - - -
   lobslev = int(pdata[lp+46]);
//- - - - - - - - - - primary particle code - - - - - - - - - - - - - - -
   qdata[1] = pdata[lp+2];
//- - - - - - - - - - primary particle energy in GeV  - - - - - - - - - -
   qdata[2] = pdata[lp+3];
//- - - - - - - - - - - run number  - - - - - - - - - - - - - - - - - - -
   qdata[3] = pdata[lp+43];
//- - - - - - - - - - simulated shower number - - - - - - - - - - - - - -
   qdata[4] = pdata[lp+1];
//- - - - - - - - - - copy some parameters and all observation levels - -
   for( la=44; la<61; la++ ) { qdata[la] = pdata[lp+la-1]; }
//- - - - - - - - - - - observation level (meter) - - - - - - - - - - - -
   qdata[11] = 1.e-2 * pdata[lp+47];
//- - - - - - - - - - - theta in degrees- - - - - - - - - - - - - - - - -
   qdata[12] = pdata[lp+10] * 57.295779513e0;
//- - - - - - - - - - - phi in degrees- - - - - - - - - - - - - - - - - -
   qdata[13] = pdata[lp+11] * 57.295779513e0;
//- - - - - - - - - - height of first interaction (km)  - - - - - - - - -
   if ( pdata[lp+6] < 0. ) pdata[lp+6] = fabs(pdata[lp+6]);
   qdata[14] = 1.e-5 * pdata[lp+6];
//- - - - - - - - - - height of first interaction (gramm/cm2) - - - - - -
   qdata[15] = thickgr(pdata[lp+6]);
   obslev = double(pdata[lp+6]);
   obsgrm = thickgr(obslev);
   if ( 2838 <= int(obslev*1.00001234e-2) &&  
        int(obslev*1.00001234e-2) <= 2844 ) {
      obsgrm = thicksouth(obslev);
      qdata[15] = obsgrm;
   }
//- - - - - - - - - - runtime of light from first intact to obslev- - - -
   qdata[16] = ( pdata[lp+6] - pdata[lp+46+lobslev] )
             / ( velight * cos(pdata[lp+10]) );
   return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeLong = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock of longitudinal tables
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeLong()
{
   // analyze longitudinal tables of 
   // vertical depth of step n, number of gammas, positrons,
   // electrons, muons(+), muons(-), hadrons, all charged,
   // nuclei, and cherenkov photons at step n, test fit parameters.    
   // see separate program code `rcorsik1.cpp`.
   // cout << "        LONG subblock ignored. " << endl;
   return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzePart = = = = = = = = = = = = = = =
//        analyze or check corsika particle data subblock
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzePart()
{
  for( lz=lp; lz<lp+nsblock; lz+=ndatpar ) {

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
    if ( mcher <= 0 ) {
      // - - - - - - cherenkov bunch data within particle data - - - - - -
      if ( 5657000. <= pdata[lz] ) { // && pdata[lz] < 1.95e7 ) {
      // ( 9900000. // pdata[lz] = 99.e5 + nint(photons in bunch)*10 + 1;
         busize = (pdata[lz]-99.e5-1.)/10.;
         mbsize = int(busize)+1;
         iobs = 1;
         if ( irec == 1 && lz < 1234 )
            printf("      cdata[%d]= %f    %f   %f\n",
               lz,busize,pdata[lz+1],pdata[lz+2]);
      } 
      else if ( 1001. <= pdata[lz] ) { // < 5657000. ) {
         // - - - - - valid particle id found > 1000 - - - - - - - - - - -
         icod = int(pdata[lz]);
         iobs = icod % 10;    
         if ( iobs == 0 ) iobs = 10; 
         igen = (icod % 1000) / 10; // - generation of interaction.
         idpa = icod / 1000; // - particle id in corsika notation.
         // - - - - - - - - calculate total energy also for nuclei:
         if ( idpa < 200 ) { 
            equad = parmas[idpa]*parmas[idpa] + pdata[lz+1]*pdata[lz+1]
                    + pdata[lz+2]*pdata[lz+2] + pdata[lz+3]*pdata[lz+3];
         }
         else if ( idpa == 9900 ) {
            equad = parmas[1]*parmas[1] + pdata[lz+1]*pdata[lz+1]
                    + pdata[lz+2]*pdata[lz+2] + pdata[lz+3]*pdata[lz+3];
         }
         else if ( idpa < 6000 ) {
            la = idpa%100; // number of protons in nucleus.
            lb = int(double(idpa/100)); // number of neutrons in nucleus.
            equad = parmas[14]*double(la) + parmas[13]*double(lb-la);
            equad = equad*equad + pdata[lz+1]*pdata[lz+1]
                  + pdata[lz+2]*pdata[lz+2] + pdata[lz+3]*pdata[lz+3];
         }
         epart = sqrt(equad);
         eplog = log10(epart);
         // - - - - - - - - switch coordinates to meters instead of cm. 
         rx = 1.e-2 * pdata[lz+4];
         ry = 1.e-2 * pdata[lz+5];
         // =========== muons =========== 
         if ( idpa == 5 || idpa == 6 ) {
            // - - - - - - - - - - - count truncated muons:
            if (1600. <= rx*rx+ry*ry && rx*rx+ry*ry <= 40000.) mutr++;
         } 
         // ========== hadrons ========== 
         if ( 7 <= idpa && idpa <= 65 ) { 
            // ========== pion(0) ========== 
            if ( idpa == 7 ) { mpinull++; } 
            // ========== pion(+) ========== 
            if ( idpa == 8 )  { mpiplus++; } 
            // ========== pion(-) ========== 
            if ( idpa == 9 ) { mpiminus++; } 
         } // end_of_hadrons. 
      }
    } // end-of ( mcher <= 0 ).
    // - - - - - - separate cherenkov bunch data - - - - - - - - - - - - -
    else if ( mcher > 0 ) { // i.e. bunchsize.   
      busize = pdata[lz];
      mbsize = int(busize)+1;
      iobs = 1;
      if ( 0.9997*fabs(1.*mcher) < pdata[lz] && irec < 10 )
         printf("      pdata[%d]= %f    %f   %f\n",
            lz,pdata[lz],pdata[lz+1],pdata[lz+2]);
    } // end-of ( mcher > 0 ).

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  } // end-of for( lz=lp; 
   
  return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeRune = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock run end
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeRune()
{
   // cout << "      RUNE " << endl;
   return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = analyzeRunh = = = = = = = = = = = = = = =
//        analyze or check corsika data subblock run header
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void analyzeRunh()
{
   // cout << "      RUNH " << endl;
   return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
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


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = pamafill= = = = = = = = = = = = = = = = =
//                    fill particle masses up to 5928_Ni; valid <= 5656.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void pamafill()
{
   double az,eb;
   int ia,in,ip,is;

// - - - - - - initialize particle masses, charges, names:
   /* const int ipartyp[202] = {
                 0,1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
           21,22,23,24,25,26,27,28,29,30,31,32,17*0,50,51,52,53,54,55,
         56,57,58,59,60,61,62,63,64,65,66,67,68,69,0,71,72,73,74,41*0,
        116,117,118,119,120,121,122,123,124,125,126,127,128,0,130,131,
         132,133,134,2*0,137,138,139,140,141,142,143,144,145,3*0,149,
           150,151,152,153,154,155,156,157,3*0,161,162,163,7*0,171,
                 172,173,25*0,199,0,0 }; */
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
   /* char chptext[202][20] = {"                   ", 
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
      */
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
// = = = = = = = = = = = = = = = printable = = = = = = = = = = = = = = = =
//        formatted output table of shower quantities (to txt-file):
//    primary, energy, runnumb, simnrsh, #hadrs, #muons, #gammas, #elecs,
//            #nkgelecs, obslev, theta, phi, h1km, h1gr. 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void printable(int ish)
{
   if ( ish == 1 ) cout << " primary     energy     runnumb  simnrsh    "
      << "   #hadrs       #muons      #gammas       #elecs    #nkgelecs  " 
      << "  obslev   theta    phi     h1km    h1gr" << endl; 
   cout.precision(1);
   cout.setf(ios::fixed,ios::floatfield);
   cout << setw(8) << qdata[1];
   cout.precision(4);
   cout.setf(ios::scientific,ios::floatfield);
   cout << setw(13) << qdata[2];
   cout.precision(1);
   cout.setf(ios::fixed,ios::floatfield);
   cout << setw(10) << qdata[3];
   cout << setw(9)  << qdata[4];
   cout.precision(5);
   cout.setf(ios::scientific,ios::floatfield);
   cout << setw(13) << qdata[5];
   cout << setw(13) << qdata[6];
   cout << setw(13) << qdata[7];
   cout << setw(13) << qdata[8];
   cout << setw(13) << qdata[9];
   cout.precision(1);
   cout.setf(ios::fixed,ios::floatfield);
   cout << setw(10) << qdata[11];
   cout.precision(2);
   cout << setw(8) << qdata[12];
   cout << setw(8) << qdata[13];
   cout << setw(8) << qdata[14];
   cout << setw(8) << qdata[15];
   cout << endl; 
   return;
}


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = printpartic = = = = = = = = = = = = = = =
//         print particle quantities of the first particle data subblock
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void printpartic(int m)
{
   cout << setw(5) << int(pdata[m]/1000.) << "." << int(pdata[m])%10;
   cout.precision(5);
   for( int i=1; i<ndatpar; i++ ) {
      cout.setf(ios::scientific,ios::floatfield);
      cout << setw(13) << pdata[i+m];
   }
   cout.precision();
   cout << endl;
   return;
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


// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = pamasscalc= = = = = = = = = = = = = = = =
//                 calculate particle masses up to 59_Ni.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void pamasscalc()
{
   double az,eb;
   int ia,in,ip,is;
   for( ia=0; ia<=76; ia++ ) {
      pmass[ia] = parmas[ia];
      psign[ia] = chargs[ia];
   }
   for( ia=77; ia<=6000; ia++ ) {
      pmass[ia] = 0.;
      psign[ia] = 0.;
   }
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
         if (ip%2 == 0 && in%2 == 0)
            eb = eb + 33.5 * pow(az,-0.75);
         else
            if (ip%2 == 1 && in%2 == 1)
               eb = eb - 33.5 * pow(az,-0.75);
         if (eb > 0.) eb=1.e-3*eb;
         else eb=0.; 
      }
   }
   // - - - masses of multineutron clusters.
   for( in=1; in<60; in++ ) {
      is = in * 100;
      pmass[is] = parmas[13] * in;
      psign[is] = 0.;
   }
   return;
}

