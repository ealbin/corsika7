//========================================================================
// 
//  r c o r s i k 2 b e o k . c p p
//  ===============================
//    This program is able to read all binary corsika particle data files 
//    of "standard" and "thinning" simulations to test on completeness;
// - - - - - - - - - - - - - - CompLinkRun - - - - - - - - - - - - - - - -
// CompLink:
//    g++ -O0 -fbounds-check rcorsik2beok.cpp -o rcorsik2beok -lm
// RunProgr:
//    ./rcorsik2beok < rcorsik2beok.iname
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// input-file rcorsik2beok.iname (file names without quotes):
//    DAT000333 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         RUNH = 211285.2812500000000000;
//         EVTH = 217433.0781250000000000;
//         LONG =  52815.2968750000000000;
//         EVTE =   3397.3918457031250000;
//         RUNE =   3301.3325195312500000;
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
// = = = = = = = = = = = = = = global constants  = = = = = = = = = = = = =

   const int nrecstd = 22940;           // "standard corsika" record size,
   const int nrecthi = 26216;           // "thinning corsika" record size,
   const int ndim = 20000000;
   const int nsubblo = 21;
   const int nsblstd = 273;
   const int nsblthi = 312;
   const int nmaxrec = 2180000;              // valid up to 50 GBytes.
   const int numbstd = nrecstd / 4;          // =  5735.
   const int numbthi = nrecthi / 4;          // =  6554.
   int ndatpar, nreclen, nsblock, lbit, lsub,
       nxpix=750, nypix=650, nprimry, lobslev, jobslev, ifil, iruntyp=0;

   float pdata[numbthi]; // to read a single corsika record thinned data.
   float rdata[numbthi]; // to keep first corsika record thinned data.
   float sdata[numbstd]; // to read a single corsika record standard data.
   float hdata[819]; // to read end of first corsika record thinned data.
   float zdata[2]; // to read record length information on 64bit machines. 

   double qrunh[275], qevth[275], qevte[275], qdata[101], obslev, obsgrm;

   const double velight = 29.9792458; // velocity of light in cm/nsec.
 
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = global quantities = = = = = = = = = = = = =
  
   char corsfile[2000][120];
  
   int mprimary, ibytes, mpinull, mpiplus, mpiminus,
       nshowf[21], nfils, nshow, nshof, ishow, ishof, isho, irec, itext,
       idpa, iobs, igen, icod, mutr, ii, la, lb, lz, lp, lend=0, lthi=1;
  
   double slopest, engylow, enghigh, thetlow, thehigh, obslev1,
          epart, eplog, equad; // extra memory for quantities.
 
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// = = = = = = = = = = = = = = = = = main= = = = = = = = = = = = = = = = =

int main(int argcnt, char *argvec[])
{
   char chtext[120];
   time_t tseconds;
   tm *datetime;
   if ( argcnt > 1 ) iruntyp = atoi(argvec[1]);
   tseconds = time(&tseconds);
   datetime = localtime(&tseconds);
   cout << "   " << endl;
   printf("        date ");
   printf(" %.2d. %.2d. %.4d   %.2d.%.2d%.2d   (%ld sec) \n",
      datetime->tm_mday, datetime->tm_mon+1, datetime->tm_year+1900,
      datetime->tm_hour, datetime->tm_min, int(tseconds)%60, tseconds); 

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== II. - read file names - - - - -
   for( ifil=1; ifil<2000; ifil++ ) {
      if ( cin.peek() == EOF ) {
         ifil--;
         cin.clear();
         break;
      }
      cin.getline(chtext,120);
      itext = strlen(chtext);
      chtext[itext+1] = '\0';
      strcpy(corsfile[ifil],&chtext[0]);
   }
   nfils = ifil--;

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== III. - loop over all file names read in- - - - - - - - - - - - - - -
   nshof = 0;
   nshow = 1;
   ishow = 0;
   for( ifil=1; ifil<=nfils; ifil++ )
   {
   nshof = nshof + nshowf[ifil]; // total showers till end of current file
   ishof = 0;                    // counting showers of current file.
   printf("\n             %s \n",corsfile[ifil]);

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
   }
   ishow++;

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b - - loop over all data records of current file  - - - - - - -
   for( irec=1; irec<nmaxrec; irec++ ) {
   if ( inputcors.peek() == EOF && lend == 0 ) {
      cout << "             end of particle data at EOF:" << endl;
      goto endofdata;
   }

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
         printf("             nsblock: %d ",nsblock);
         if ( lbit == 1 ) printf("   64bit file");
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
   }

// _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.b.2.- - - analyze subblocks of the record - - - - - - - - - - -
      for( lp=1; lp<nsblock*21; lp+=nsblock ) {
      lsub = 1 + lp/nsblock;
      if ( 211285.2 < pdata[lp] && pdata[lp] < 211285.3 ) {
         if ( ifil < 0 ) cout << "        RUNH subblock found. " << endl;
      }
      else if ( 217433.0 < pdata[lp] && pdata[lp] < 217433.1 ) {
         if ( ifil < 0 ) cout << "        EVTH subblock found. " << endl;      
      }
      else if ( 3397.3 < pdata[lp] && pdata[lp] < 3397.4 ) {
         if ( ifil > 0 )
            cout << "        EVTE subblock found: " << lsub << endl;
      }  
      else if ( 3301.3 < pdata[lp] && pdata[lp] < 3301.4 ) {
         if ( ifil > 0 )
            cout << "        RUNE subblock found: " << lsub << endl;
         lend = lsub;
         if ( lp > nsblock ) {
            if ( pdata[lp-nsblock] < 3397.3 || 3397.4 < pdata[lp-nsblock] ) {
               cout << "             incomplete particle data!" << endl;
               goto endofdata;
            }                  
         }  
      }
      else if ( ( 1000.0 < pdata[lp] && pdata[lp] < 3301.3 ) ||
                ( 3301.4 < pdata[lp] && pdata[lp] < 3397.3 ) ||
                  3397.4 < pdata[lp] ) {
         if ( lp == 20*nsblock && inputcors.peek() == EOF ) { 
            cout << "             end of particle data at EOF!" << endl;
            goto endofdata;
         }
      }
      } // end_of_loop of analyzing subblocks - - - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//==== III.c.- - - continue with next data record or close file. 
    } // endfor_loop over all data records - - - - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
 endofdata:
    inputcors.close(); // close corsika file also after break at EOF;
   } // endfor_loop over all file names read in  - - - - - - - - - - - - -

// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
//=== IV. - after all showers and particles- - - - - - - - - - - - - - - -
   if ( ishow > 1 )
      cout << "             " << ishow << " showers processed.";
   else if ( ishow != 1 )
      cout << " WARNING, case should not occur! ";
   cout << endl << "  " << endl;
   return 0;
}
