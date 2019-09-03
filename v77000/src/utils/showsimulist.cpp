
// showsimulist.cpp
// ======================================================================
// CompLink:
//    g++ -O0 -fbounds-check showsimulist.cpp -o showsimulist -lm
// RunProg:
//    ls -l DAT?????? | ./showsimulist [ > showsimulist.ca2-showtest ]
//    ls -l CER?????? | ./showsimulist > showsimulist.ca2-cherenkov
//    ls -l DAT??????.cher | ./showsimulist
//    ls -l DAT??????.cher-tel??? | ./showsimulist
//    ./showsimulist < showsimulist.datnames
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    reading corsika particle data files from a (one name per line)
//    tabular each simulation ([DAT,CER]iiiiii) will be printed out
//    in one line the following quantities:
//        Primary, lg(E), theta, phi, nsh,  runnr, size,
//            obslvme, h1stme, thilev, thiwmax, lg(thirad),
//                verspgm, models, rundate, Xmagn, Zmagn;
//    Cherenkov files will be identified by `_ce` after the primary,
//    particle data files must be available as [DAT,CER]iiiiii,
//    protocol files as ...lst or additional conditions must be added.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//          RUNH = 211285.2812500000000000;
//          EVTH = 217433.0781250000000000;
//          LONG =  52815.2968750000000000;
//          EVTE =   3397.3918457031250000;
//          RUNE =   3301.3325195312500000;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//                             author: juergen.oehlschlaeger@kit.edu
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//         Primary   lg(E)  theta    phi  nsh  runnr   sizeM  obslvme  h1stme   .....
// Iron       5626   16.00    0.0    0.0   1  199079     3.9  1413.82   18350.  .....
// _stck_in_     4   15.09    0.0    0.0   1  169051     2.4   194.00   18765.  .....
// Fluorine   1909   16.00    0.0    0.0   1  199070     3.9  1416.51   22222.  .....
// proton       14   16.00   30.0   -3.3   1  199080     3.4  1429.25   22224.  .....
// Manganese  5525   16.00   30.0   -3.3   1  199082    10.1  1428.13   22222.  .....
// proton _ce   14   17.50   37.9 -138.0   1  000044    90.8  -500.00  -24880.  .....
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
using namespace std;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

const double cpi=4.*atan(1.), cpi180=cpi/180., c180pi=180./cpi;

const double velight = 29.9792458; // velocity of light in cm/nsec.

const int nfmx=50000, nchx=160;

const double AATM[6] = { 0.e0,       // param of US standard atmosph.
                   -186.5562e0, -94.919e0, 0.61289e0, 0.e0, 1.128292e-2 };
const double BATM[6] = { 0.e0,       // param of US standard atmosph.
                1222.6562e0, 1144.9069e0, 1305.5948e0, 540.1778e0, 0.e0 };
const double CATM[6] = { 0.e0,       // param of US standard atmosph.
              994186.38e0, 878153.55e0, 636143.04e0, 772170.16e0, 1.e-9 };

const char chemical[102][16]={       " unknown    ",
       " Hydrogen   "," Helium     "," Lithium    "," Beryllium  ",
       " Boron      "," Carbon     "," Nitrogen   "," Oxygen     ",
       " Fluorine   "," Neon       "," Sodium     "," Magnesium  ",
       " Aluminium  "," Silicon    "," Phosphorus "," Sulfur     ",
       " Chlorine   "," Argon      "," Potassium  "," Calcium    ",
       " Scandium   "," Titanium   "," Vanadium   "," Chromium   ",
       " Manganese  "," Iron       "," Cobalt     "," Nickel     ",
       " Copper     "," Zinc       "," Gallium    "," Germanium  ",
       " Arsenic    "," Selenium   "," Bromine    "," Krypton    ",
       " Rubidium   "," Strontium  "," Yttrium    "," Zirconium  ",
       " Niobium    "," Molybdenum "," Technetium "," Ruthenium  ",
       " Rhodium    "," Palladium  "," Silver     "," Cadmium    ",
       " Indium     "," Tin        "," Antimony   "," Tellurium  ",
       " Iodine     "," Xenon      "," Caesium    "," Barium     ",
       " Lanthanum  "," Cerium     "," Praseodym. "," Neodymium  ",
       " Promethium "," Samarium   "," Europium   "," Gadolinium ",
       " Terbium    "," Dysprosium "," Holmium    "," Erbium     ",
       " Thulium    "," Ytterbium  "," Lutetium   "," Hafnium    ",
       " Tantalum   "," Tungsten   "," Rhenium    "," Osmium     ",
       " Iridium    "," Platinum   "," Gold       "," Mercury    ",
       " Thallium   "," Lead       "," Bismuth    "," Polonium   ",
       " Astatine   "," Radon      "," Francium   "," Radium     ",
       " Actinium   "," Thorium    "," Protactin. "," Uranium    ",
       " Neptunium  "," Plutonium  "," Americium  "," Curium     ",
       " Berkelium  "," Californium"," Einsteinium","            ",
       "            "};

char qpatext[202][20]={"                   ",
       " gamma             "," positron          "," electron          ",
       " _stck_in_         "," muon+             "," muon-             ",
       " pi0               "," pi+               "," pi-               ",
       " K0long            "," K+                "," K-                ",
       " neutron           "," proton            "," anti proton       ",
       " K0short           ","                   "," Lambda            ",
       " Sigma+            "," Sigma0            "," Sigma-            ",
       " Xi0               "," Xi-               "," Omega-            ",
       " anti neutron      "," anti Lambda       "," anti Sigma-       ",
       " anti Sigma0       "," anti Sigma+       "," anti Xi0          ",
       " anti Xi+          "," anti Omega+       ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   "," omega             "," rho0              ",
       " rho+              "," rho-              "," Delta++           ",
       " Delta+            "," Delta0            "," Delta-            ",
       " anti Delta--      "," anti Delta-       "," anti Delta0       ",
       " anti Delta+       "," K*0               "," K*+               ",
       " K*-               "," anti K*0          "," electron neutrino ",
       " anti elec neutrino"," muon neutrino     "," anti muon neutrino",
       "                   "," eta=>2*gamma      "," eta=>3*pi0        ",
       " eta=>pi+pi-pi0    "," eta=>pi+pi-gamma  ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   "," D0                "," D+                ",
       " anti D-           "," anti D0           "," Ds+               ",
       " anti Ds-          "," eta c             "," D*0               ",
       " D*+               "," anti D*-          "," anti D*0          ",
       " D*s+              "," anti D*s-         ","                   ",
       " J/psi             "," tau +             "," tau -             ",
       " tau neutrino      "," anti tau neutrino ","                   ",
       "                   "," Lambda c +        "," Xi c +            ",
       " Xi c 0            "," Sigma c ++        "," Sigma c +         ",
       " Sigma c 0         "," Xi c prime +      "," Xi c prime 0      ",
       " Omega c 0         ","                   ","                   ",
       "                   "," anti Lambda c -   "," anti Xi c -       ",
       " anti Xi c 0       "," anti Sigma c --   "," anti Sigma c -    ",
       " anti Sigma c 0    "," anti Xi c prime - "," anti Xi c prime 0 ",
       " anti Omega c 0    ","                   ","                   ",
       "                   "," Sigma c * ++      "," Sigma c * +       ",
       " Sigma c * 0       ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   "," anti Sigma c * -- ",
       " anti Sigma c * -  "," anti Sigma c * 0  ","                   ",
       " Stnd              ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   ","                   ",
       "                   ","                   "," Thin              ",
       "                   ","                   ","                   ",
       " Cherenkov photon  ","                   ","                   "};

   const int nrecstd = 22940;           // "standard corsika" record size,
                                        //        pdata[0] = 6.08270e-311.
                                        //    or  pdata[0] = 3.213456e-41.
   const int nrecthi = 26216;           // "thinning corsika" record size,
                                        //        pdata[0] = 6.95166e-311.
                                        //    or  pdata[0] = 3.672523e-41.
   const int nsblstd = 273;
   const int nsblthi = 312;
   const int numbstd = nrecstd / 4;     // =  5735.
   const int numbthi = nrecthi / 4;     // =  6554.

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

int mposubstring( char *chbigstr, char *chsubstr) {
   // find the position of a substring inside a bigger string:
   // ipos = mposubstring( satzc, chsub); // usage.

   int icnt=0, ifnd=0, isub=0, iret=0;

   if ( strlen(chbigstr) >= strlen(chsubstr) ) {
      iret = -1;  
      for( icnt=strlen(chbigstr)-strlen(chsubstr); icnt>=0; icnt-- ) {
         ifnd=1;
         for( isub=0; isub<int(strlen(chsubstr)); isub++ ) {
            if (chbigstr[icnt+isub] != chsubstr[isub]) {
               ifnd = 0; break;
            }
         }
         if (ifnd == 1) break;
      }
      iret = icnt;
   }

   return iret;

} // end-of-function int mposubstring.

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

float pdata[numbthi]; // to read a single corsika record of thinned data.

double fsize[nfmx], qrunh[315], qevth[315];

char corsfile[nfmx][160];

int nfsho[nfmx], nfrun[nfmx], nflen[nfmx];

// = = = = = = = = = = = = = = = = = main  = = = = = = = = = = = = = = = =

int main(int argcnt, char *argvec[])
{
   double obslev, engyprev, thetprev, phiaprev, obslprev,
       energy, theta, phia, rdate;

   int ifil, idate, itext, iprmtxt, nfilstop, nplotsh, mcodprev,
       isla, lthi, im, lthiprev, irunnr, imont, ijahr, nfil, nsblock, 
       mcode, models, mjahr, mmon, mtag, mhour, mmint;

   char chlsline[256], chtext1[14][256], chrunnr[16], cprimary[16],  
        chmonth[18][8]={"   ","Jan","Feb","Mar","Apr","May","Jun","Jul",
            "Aug","Sep","Oct","Nov","Dec","Mär","Mai","Okt","Dez","   "};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - initializations:
   iprmtxt = 14;
   nfilstop = 0;
   obslev = 0.1;
   nplotsh = 0;
   lthi = 0;

// --read first file name and fix position of file size-------------------
// ------------------------------------------------------
      ifil = 1;
      cin.getline(chlsline,250);
      itext = strlen(chlsline);
      chlsline[itext+1] = '\0';
      if ( strcmp(&chlsline[0],"-rw") == 0 &&
           strcmp(&chlsline[0],"-r-") == 0 ) {
         printf("              Use of `ls -1 ...` not supported, but"
                " `ls -l ...` will be ok;\n");
         printf("              or file access rights do not begin"
                " as `-rw` or `-r-`; first use `chmod u+w *`;\n\n");
         goto continue445;
      }
   // check and copy first particle data file name:
      sscanf(chlsline,"%s %s %s %s %s %s %s %s %s",
         chtext1[1],chtext1[2],chtext1[3],chtext1[4],
         chtext1[5],chtext1[6],chtext1[7],chtext1[8],chtext1[9]);
      strcpy(corsfile[ifil], chtext1[9]);
   // check and fix file size of particle data file:   
      fsize[ifil] = 1.e-6 * atof(chtext1[5]); // now in Mbytes.
      if ( fsize[ifil] < 0.1 ) fsize[ifil] = 0.1;
   // copy current run number to integer element:
      strcpy(chrunnr, &corsfile[ifil][3]);
      strcpy(&chrunnr[6], "\0");
      irunnr = atoi(chrunnr);
      nflen[ifil] = strlen(chtext1[9]); // length of data file name.
      nfrun[ifil] = irunnr; 
      nfsho[ifil] = 1;
   // convert date of copied particle data file to integers:
      mtag = atoi(chtext1[7]);
      mmon = 0; mmint = 0; mhour = 0;
      for( im=1; im<17; im++ ) {
         mmon = mposubstring( chlsline, chmonth[im]);
         if ( mmon > 25 ) {
            mmon = im;
            if ( mmon == 13 ) mmon = 3;  // "Mär".
            if ( mmon == 14 ) mmon = 5;  // "Mai".
            if ( mmon == 15 ) mmon = 10; // "Okt".
            if ( mmon == 16 ) mmon = 12; // "Dez".
            break;
         }
      }
      mjahr = 2017;
      if ( strcmp(chtext1[8],":") != 0 ) mjahr = atoi(chtext1[8]);
      else {
         mmint = atoi(&chtext1[8][3]);
         strcpy(&chtext1[8][2], " ");
         mhour = atoi(chtext1[8]);
      }
      rdate = 10000.*mjahr + 100.*mmon + mtag + 1.e-2*mhour + 1.e-4*mmint;
      if ( rdate < 0. ) printf("        rdate= %18.4f \n",rdate);   
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// --read all run parameters including file names (more than 1 file)------
   for( ifil=2; ifil<nfmx; ifil++ ) {
      if ( cin.peek() == EOF ) { cin.clear(); break; }
      cin.getline(chlsline,250);
      itext = strlen(chlsline);
      chlsline[itext+1] = '\0';
   // check and copy ifil-th particle data file name:
      sscanf(chlsline,"%s %s %s %s %s %s %s %s %s",
         chtext1[1],chtext1[2],chtext1[3],chtext1[4],
         chtext1[5],chtext1[6],chtext1[7],chtext1[8],chtext1[9]);
      strcpy(corsfile[ifil], chtext1[9]);
   // check and fix file size of particle data file:   
      fsize[ifil] = 1.e-6 * atof(chtext1[5]); // now in Mbytes.
      if ( fsize[ifil] < 0.1 ) fsize[ifil] = 0.1;
   // copy current run number to integer element:
      strcpy(chrunnr, &corsfile[ifil][3]);
      strcpy(&chrunnr[6], "\0");
      irunnr = atoi(chrunnr);
      nflen[ifil] = strlen(chtext1[9]); // length of data file name.
      nfrun[ifil] = irunnr;
      nfsho[ifil] = 1;
   }
   nfil = ifil--;

// --print title line-----------------------------------------------------
// --------------------
   printf("  \n");
   printf("         primary   lg(E)  theta   phi   nsh  "
      " runnr   sizeM  obslvme  h1stme  thilev  wmax  thirad"
      "  verspgm    models   rundate  Xmagn  Zmagn\n");
   mcodprev = 0;
   engyprev = 0.;
   thetprev = 0.;
   phiaprev = 0.;
   obslprev = 0.;
   lthiprev = 0;

// --loop over all particle data files of `ls -l` command-----------------
// -----------------------------------------------------------------------
   printf("  \n");
   nfilstop = 0;
   for( ifil=1; ifil<nfil; ifil++ ) {
      irunnr = nfrun[ifil];
   // - - - - - - - - - - - - - - - - - - - - - - -
      ifstream inputcors(corsfile[ifil],ios::in); // allocate file.
      if ( inputcors.peek() == EOF ) {
         printf("         file ==> does not exist at this path! Stop!\n");
         goto continue443;
      }
      else {
         inputcors.read((char*)&pdata,sizeof(pdata));
         lthi = -1;
         if (217433.0 < pdata[nsblstd+1] && pdata[nsblstd+1] < 217433.1)
            lthi = 0;
         if (217433.0 < pdata[nsblthi+1] && pdata[nsblthi+1] < 217433.1)
            lthi = 1;
         nsblock = nsblstd + lthi*39; // ndatpar = nsblock / 39;
         for( im=0; im<=nsblock; im++ )
              qrunh[im] = pdata[im];
         for( im=nsblock+1; im<=nsblock+nsblstd; im++ )
              qevth[im-nsblock] = pdata[im];
      }
      if ( fabs(pdata[1]) < 1.e-6 ) goto continue443; 
      nfilstop++;
   // - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - check on primary particle:
      if ( 0 < qevth[3] && qevth[3] < 200. ) {
         iprmtxt = int(qevth[3]);
      }
      else if ( qevth[3] < 6000. ) {
         iprmtxt = 200;
         im = int(qevth[3])%100;
         if ( im > 1 ) { 
            strcpy(qpatext[iprmtxt], chemical[im]);
         }
         else {
            if ( int(qevth[3]) == 201 ) 
               strcpy(qpatext[iprmtxt], " Deuteron          ");
            else if ( int(qevth[3]) == 301 )
               strcpy(qpatext[iprmtxt], " Tritium           ");
            else
               strcpy(qpatext[iprmtxt], chemical[im]);
         }
         strcpy(cprimary, qpatext[iprmtxt]);
         strcpy(&cprimary[11], "\0");
      }
      else {
         printf("       found invalid particle id %f \n",qevth[3]);
         // should not - but may - occur of special simulations.
      }
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - check models and date of simulation - - - - -
      if ( int(qevth[77]) > 7 ) qevth[77] = 7.0123;
      models = 0;
      for( im=73; im<=80; im++ ) {
         models = models + pow(10,(80-im))*int(qevth[im]);
      }
      imont = (int(qevth[45])) % 10000;
      ijahr = int(qevth[45]) / 10000;
      if ( ijahr < 30 ) ijahr = ijahr + 2000;
      if ( ijahr < 100) ijahr = ijahr + 1900;
      idate = 10000. * ijahr + imont;
      if ( idate%100 > 31 ) idate = 31 + 100*int(idate/100);
      if ( idate%100 == 0 ) idate =  1 + idate;
      if ( qevth[148] == 0. ) qevth[148] = 1.;
      strcpy(cprimary, qpatext[iprmtxt]);
      strcpy(&cprimary[11], "\0");
   // - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - check on thinning parameters - - - - -
      if ( qevth[151] < 1. ) {
         isla = log10(qevth[4]) + log10(qevth[148]+4.67735141e-34);
         qevth[151] = 0.;
         for( im=1; im<=isla; im++ ) {
            qevth[151] = qevth[151] + pow(10.,(isla-im))*double(10-im);
         }
      }
      if ( qevth[152] > 1. ) {
         qevth[152] = 0.01234 + log10(qevth[152]*1.e-2);
      }
      else {
         qevth[152] = 0.;
         qevth[151] = 0.;
      } 
   // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - check previous parameters of shower - - - - -
      if ( iprmtxt > 0 ) {
        mcode = int(qevth[3]);
        // - - - - - - mcode=0 for stackin simulations:
        if ( mcode == 0 ) mcode = 4;
        energy = qevth[4];
        if ( qevth[4] > 0. ) energy = 9. + log10(qevth[4]);
        theta = c180pi*qevth[11];
        phia = c180pi*qevth[12];
        if ( mcode != mcodprev || energy != engyprev ||
             theta != thetprev || phia != phiaprev ||
            obslev != obslprev || lthi != lthiprev ) {
            // print/check extra blank line.
            // if ( ifil > 1 ) printf("  \n");
        }
   // - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - print out tabular quantities- - - - -
        if ( fsize[ifil] < 10000. ) {
          printf("%s%5d  %6.2f %5.1f %6.1f %5d  %.6d"
          " %8.1f %7.1f %7.0f. %5.1f %8.1e %5.1f %8.5f  %.8d"
          " %9d %6.2f %6.2f \n", cprimary, int(qevth[3]), energy,
          c180pi*qevth[11], c180pi*qevth[12], nfsho[ifil],
          int(qevth[44]), fsize[ifil], 1.e-2*qevth[48],
          1.e-2*qevth[7], log10(qevth[148]+4.67735141e-34),
          qevth[151], qevth[152]+0.001, 1.e-5+qevth[46], models,
          idate, qevth[71], qevth[72]); // "7.0f." ("%8.1e")
        }
        else {
          printf("%s%5d  %6.2f %5.1f %6.1f %5d  %.6d"
          " %7.0f. %7.1f %7.0f. %5.1f %8.1e %5.1f %8.5f  %.8d"
          " %9d %6.2f %6.2f \n", cprimary, int(qevth[3]), energy,
          c180pi*qevth[11], c180pi*qevth[12], nfsho[ifil],
          int(qevth[44]), fsize[ifil], 1.e-2*qevth[48],
          1.e-2*qevth[7], log10(qevth[148]+4.67735141e-34),
          qevth[151], qevth[152]+0.001, 1.e-5+qevth[46], models,
          idate, qevth[71], qevth[72]); // "7.0f." ("%8.1e")
        } 
      }
   // - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - copy current quantities to `old` variables:
      mcodprev = mcode;
      engyprev = energy;
      thetprev = theta;
      phiaprev = phia;
      obslprev = obslev;
      lthiprev = lthi; 
continue443:
      inputcors.close();
      nplotsh = 0;
   } // --end-of loop ifil=1...nfil.--------------------------------------

// --print closing comment lines------------------------------------------
// -------------------------------
   nplotsh = 0;
   printf("  \n");
   if ( nplotsh > 0 ) {
      printf("         primary   lg(E)  theta   phi   nsh  "
         " runnr   sizeM  obslvme  h1stme  thilev  wmax  thirad"
         "  verspgm    models   rundate  Xmagn  Zmagn   "
         "hadron-  muon-  electr- photoncut \n");
   }
   else {
      printf("         primary   lg(E)  theta   phi   nsh  "
         " runnr   sizeM  obslvme  h1stme  thilev  wmax  thirad"
         "  verspgm    models   rundate  Xmagn  Zmagn \n");
   }
   printf("  \n");

// --print explanation of model digits------------------------------------
// -------------------------------------
continue445:
   printf("              Total number of files: %7d \n", nfilstop);
   printf("              Appendix `_ce`: Cherenkov output read,"
      " otherwise no such extra info (i.e. only blanks).\n");
   printf("              Explanation of column `models`: "
      " (10^7): EGS flag;   (10^6): NKG flag;   (10^5):"
      " lowEnergy flag, 1=gheisha, 2=urqmd,\n");
   printf("              3=fluka; (10^4): highEnergy, 0=hdpm,"
      " 1=venus, 2=sibyll, 3=qgsjet, 4=dpmjet, 5=nexus, 6=epos;"
      " (10^3): Cherenkov flag;\n");
   printf("              (10^2): Neutrino flag; (10^1): Curved"
      " flag, 0=standard, 2=curved; (10^0):"
      " Computer, 3=unix, 4=macintosh.\n");
   printf("              thirad: radial thinning to 10^thirad meter.\n");
   printf("  \n  \n  \n");

   return 0;

}
  
// = = = = = = = = = = = = = = = heightcm  = = = = = = = = = = = = = = = =
//     calculate height (cm) a.s.l. for a given thickness (gramms/cm**2)
//     US standard atmosphere
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
//     thickgr (gramms/cm**2) of atmosphere depending on height (cm)
//     US standard atmosphere
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

// = = = = = = = = = = = = = = = thicksouth = = = = = = = = = = = = = = =
//     thicksouth (gramms/cm**2) of atmosphere depending on height (cm)
//     for the South Pole atmosphere of October 01.
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
  
