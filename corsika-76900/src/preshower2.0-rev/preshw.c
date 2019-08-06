/*====================================================================*/
/*      preshw.c                                                      */
/*--------------------------------------------------------------------*/
/* HEADER DECLARATIONS & DEFINITIONS                                  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "veto.h"
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/*--------------------------------------------------------------------*/
/*------------------- PRESHOWER WITH VETO ALGORITHM-------------------*/
/*--------------------------------------------------------------------*/
/*                                                                    */
/* In PRESHOWER 2.0 the veto algorithm replaces a standard method of  */
/* propagation with 2 loops - external over trajectory with a constant*/
/* step and the inner loop over particles. In the veto algorithm the  */
/* step is being computed dynamically (see the publication for the    */
/* references)                                                        */
/* DESCRIPTION of the procedure:                                      */
/* Propagation of ultra high energy photon before its entering the    */
/* earth's atmosphere is simulated, taking into account gamma         */
/* conversion into e+/- pair and subsequent bremsstrahlung of the     */
/* electrons in the geomagnetic field. As a result we obtain a bunch  */
/* (a preshower) of particles, mainly photons and a few electrons,    */
/* instead of the primary gamma. The information about all the        */
/* particles in the preshower is returned to CORSIKA or saved (if in  */
/* stand-alone mode).                                                 */
/*                                                                    */
/* Piotr Homola <Piotr.Homola@ifj.edu.pl>                             */
/*--------------------------------------------------------------------*/
/* Note:                                                              */
/* in C language: float: 32 bits,  range 3.4*1e-38  to 3.4*1e+38      */
/*                double: 64 bits, range 1.7*1e-308 to 1.7*1e+308     */
/*--------------------------------------------------------------------*/
/*                                                                    */
/* INPUT:                                                             */
/*  int *id: particle id (always id=1 (gamma))                        */
/*  double pE_gamma: initial energy of particle in GeV                */
/*  float *pthe_loc, *pphi_loc: zenith and azimuth angles of the      */
/*             shower (in radians) in the local coordinate system     */
/*             (as defined in CORSIKA)                                */
/*  float *ptoa: top of the atmosphere in cm                          */
/*  double *pglong: longitude position in deg. of experiment          */
/*             (Greenwich = 0., eastward is positive)                 */
/*  double *pglat:  latitude position in deg. of experiment           */
/*             (Northpole = +90., Southpole = -90.)                   */
/*  float *pigrf_year: the year in which the magnetic field is to be  */
/*             computed (this parameter is an input for the igrf      */
/*             procedure)                                             */
/*  int *iiprint: print flag = 0: suppress printing from preshw       */
/*                           = 1: enable   printing from preshw       */
/*  int *pcorsika:  1: source compatible with CORSIKA - random        */
/*                     generator from CORSIKA is used (rndm);         */
/*                  0: random generator from NR (ran2) ok for         */
/*                     stand-alone version                            */
/*  int *pnruns:    number of runs                                    */
/*  config *pconf: a structure which can hold different config        */ 
/*                 parameters for preshw_veto. One can redefine       */
/*                 config strucuture implementation in utils.h        */
/*                                                                    */
/* OUTPUT:                                                            */
/*  part_out: double array 100000x8, n: from 0 to the last secondary, */
/*             max 99999:                                             */
/*  part_out[n][0]: particle id: 1 gamma, 2 e+, 3 e-                  */
/*  part_out[n][1]: Energy [GeV]                                      */
/*  part_out[n][2]: cos(pi_the_loc) for particle trajectory (not used */
/*             since all preshower particles are assumed to have      */
/*             trajectories the same as the primary particle)         */
/*  part_out[n][3]: pi_phi_loc as above                               */
/*  part_out[n][4]: t as above, relative to t=0 for prim. part.       */
/*             passing top of atm (not used since all the particles   */
/*             are assumed to enter the atmosphere at same time)      */
/*  part_out[n][5]: x_loc: coordinates of the particle in cm in the   */
/*             local system (not used - all preshower particles are   */
/*             assumed to have trajectories the same as the primary   */
/*             particle, thus x_loc=0 and y_loc=0)                    */
/*  part_out[n][6]: y_loc, related to the shower axis at the top of   */
/*             atmosphere (not used, see above)                       */
/*  part_out[n][7]: it holds actual location of 'n' particle          */
/*  float *r_zero: [km] height above sea level for the simulation     */
/*             starting point                                         */
/*  float *r_fpp: [km] height a.s.l. for the first pair production    */
/*             (0 for surviving photons)                              */
/*  int *N_part_out: number of preshower particles at the top of the  */
/*             atmosphere (number of entries in part_out array)       */
/*--------------------------------------------------------------------*/
/* FRAMES of REFERENCE:                                               */
/* local: centered at the observatory site (theta,phi,R_earth)        */
/*            x-axis towards geogr. N, y-axis westward, z-axis upward */
/*            phi_local=0 for y=0, x>0                                */
/* local CORSIKA: centered at observatory site (theta,phi,R_earth)    */
/*             x towards magnetic North (along the projection of the  */
/*             local B line on the xy plane), y - westward, z-upward; */
/*             'local CORSIKA' is 'local' rotated around z-axis by    */
/*             DECLINATION angle (the angle between geographical      */
/*             North and the local direction to the magnetic North)   */
/*             phi_local=0 for y=0, x<0 (difference of 180 deg        */
/*             comparing to 'local')                                  */
/* global: centered in the Earth's center,                            */
/*             z - Earth's axis towards North, x - towards Greenwich, */
/*             y - eastward from x; phi=0 for y=0,x>0                 */
/*             (at Greenwich), theta=0 at geogr. North                */
/*--------------------------------------------------------------------*/

void preshwveto_(int *id, double *pE_gamma, float *pthe_loc,
                 float *pphi_loc, float *ptoa, double *pglong,
                 double *pglat, float *pigrf_year,
                 int *iiprint,int *pcorsika, int *pnruns,
                 double part_out[100000][8], float *r_zero, float *r_fpp,
                 int *N_part_out)
//                 int *N_part_out, config *pconf)
{
    float order,x,B_tr,r_0,theta,phi,r_curr,r_curr_0,
    r_step,E_slice,rand_b,the_loc,phi_loc,glong,
    glat,toa,tr_par,sing,cosg,the_deg,phi_deg,surv_percent,ppf;
    float wektor[4],vz[4],sites[4],sitec[4],r_loc[4], r_glob_0[4],
    tn[4],B_site_f[4],B_site_loc[4],vy[4],vy_C[4]   ;
    double B_site[4];

    double alpha_gamma,unit_dist,*i,*k,*di,*dk,i_ok,k_ok,di_ok,
    dk_ok,E_gamma,E_conv_tresh,chi_limit,E_prev,E_e,E_sg,E_sg0,
    I_brem,dE,E_int,E_ch,frac,Etot_g,Etot_e,Ee_max,Eg_max,E_0,
    prob_factor,prob_brem_total;
    int file_part_out,multirun_file,ii,jj,mm,nrun,conversion,last_ii,
    first_conv,n_surv,n_phot,n_eplus,n_eminus,corsika,
    igrf_year,part_out_size,iprint;
    long *seed1,seed_ok;
    FILE *file_1,*file_3;
    char fname_1[30],fname_3[30];
    srand(  time (NULL)  );

    /*-----------ADVANCED PARAMETERS, compile after changing -------------*/
    nrun=*pnruns;         /*number of runs                              */
    E_sg0=(1e+12)*eV2erg; /* preshower particles below this energy      */
    /* level are neglected                        */
    r_0=5.0*R_e;
    /* starting distance in Earth's radii (the    */
    /* length of trajectory from the beginning to */
    /* the Earth's surface)                       */
    multirun_file=0;      /* 1: if in multirun mode write output        */
    /*    (n_surv, fpp,...) to file multirun.dat; */
    /* 0: no writing                              */
    chi_limit=0.1;        /* a parameter important in computing         */
    /* pair prod. function T                      */
    igrf_year=*pigrf_year;/* input parameter for the field routine      */
    /* restrictions (limitations of the IGRF      */
    /* routine): 1965 <= igrf_year <= 2005;       */
    /* if one sets igrf_year outside this interval*/
    /* the igr_year value is reset automatically  */
    /* to the relevant limit, e.g. user sets      */
    /* igrf_year=1960 -> the procedure resets     */
    /* igrf_year=1965                             */
    iprint=*iiprint;      /* 1: printing enabled                        */
    /* 0: printing disabled (but errors printed)  */
    corsika=*pcorsika;    /* 1: source compatible with CORSIKA - random */
    /*    generator from CORSIKA is used (rndm);  */
    /* 0: random generator from NR (ran2) ok for  */
    /*    stand-alone version                     */
    file_part_out=0;      /* 1: part_out written to file part_out.dat;  */
    /* 0: no writing                              */
    part_out_size=99999;  /* upper limit for the size of the preshower  */
    /* array (part_out[][]); exit the program if  */
    /* there are more particles than              */
    /* part_out_size                              */
    E_conv_tresh=1e+16;   /* (eV) below this threshold photons are not  */
    /* checked for conversion                     */
    r_step=1000.0; /* step in km in the initial loop over trajectory to */
    /* calculate the maximum probabilities of conversion and            */
    /* bremsstrahlung along the track                                   */
    

    /*--------------------------------------------------------------------*/

    /*----------------- OPENING FILES ------------------------------------*/
    /* open multirun file: a line contains the summary of a run           */
    /*   (stand alone mode only)                                          */
    /*--------------------------------------------------------------------*/
    if (multirun_file==1)
    {
        sprintf(fname_3,"multirun.dat");
        if ((file_3=fopen(fname_3,"w")) == NULL)
        {
            printf(" Can not open file<br>\n");
        }
    }
    /*--------------------------------------------------------------------*/
    /* a file to save data on each particle in the preshower at the top   */
    /*  of atmosphere (toa) level:                                        */
    /*--------------------------------------------------------------------*/
    if (file_part_out==1)
    {
        sprintf(fname_1,"part_out.dat");
        if ((file_1=fopen(fname_1,"w")) == NULL)
        {
            printf(" Can not open file<br>\n");
        }
    }

    /*------------- read and preprocess input parameters -----------------*/
    /* initial gamma energy in eV                                         */
    E_gamma=*pE_gamma/eV2GeV;     /* now E-gamma in eV                  */
    the_loc=*pthe_loc;            /* shower trajectory in radians       */
    phi_loc=*pphi_loc+PI;
    /* in CORSIKA phi=0 means shower coming from magnetic SOUTH,          */
    /* in preshower all calculations are made for phi=0 -> shower from    */
    /* mag. NORTH, have to convert phi before starting simulations        */
    /* in the stand-alone mode this shift is also needed since there      */
    /* is no rotation by g.                                               */
    the_deg=the_loc*rad2deg;      /* shower trajectory in degrees       */
    phi_deg=phi_loc*rad2deg;
    /* more suitable format of azimuth in the printout                    */
    if (corsika==0)
    {
        phi_deg=phi_deg-180;
    }
    else
    {
        phi_deg=phi_deg-360;
    }
    toa=*ptoa;                    /* top of atmosphere                  */
    /* geographical position of the observatory site in global spherical  */
    /* coord. theta and phi                                               */
    glong=*pglong;                /* site longitude                     */
    glat=*pglat;                  /* site latitude                      */
    theta=(90.0-glat)*deg2rad;    /* site spherical theta, 0 at North   */
    phi=glong*deg2rad;            /* site spherical phi, 0 at Greenwich,*/
    /* positive eastwards                 */
    sites[1]=R_e;                 /* sites: global spherical coords. of */
    sites[2]=theta;               /* the observatory                    */
    sites[3]=phi;
    sph2car(sites,sitec);         /* get sitec: global cart. coords. of */
    /* the observatory                    */
    /*------------ computing g (= DECLINATION - 180 deg) -----------------*/
    /* declination of the local magnetic field is the angle between the   */
    /* local field and the local direction to the geographic north        */
    /* (positive declination means that the local field points eastward   */
    /* with respect to the geographic north (detailed definition at       */
    /* www.ngdc.noaa.gov). Declination is needed to transform local       */
    /* coordinates (for local frame see the CORSIKA manual) to the global */
    /* ones.                                                              */

    /* 1. get local B (B_site) vector in global coord.                    */
    b_igrf(igrf_year,sites[1],sites[2],sites[3],B_site);
    /* input float: sites: site global, spherical coordinates             */
    /* output double: B_site: local field vector in global cartesian      */
    /* coordinates                                                        */

    /* 2. convert B_site[] from double to float (B_site_f[]) to suit the  */
    /*    format of glob2loc:                                             */
    B_site_f[1]=B_site[1];
    B_site_f[2]=B_site[2];
    B_site_f[3]=B_site[3];

    /* 3. get B in local coords (B_site_loc)                              */
    glob2loc(B_site_f,B_site_loc,theta,phi);

    /* 4. get normalized local CORSIKA y-axis (vy_C)                      */
    /* local z axis: z axis in local cart. coords                         */
    vz[1]=0.0;
    vz[2]=0.0;
    vz[3]=1.0;
    cross(vz,B_site_loc,wektor);
    norm(wektor,vy_C);

    /* 5. get the g angle: the angle between vy and vy_C - using          */
    /*    cross product of vy and vy_C                                    */
    /* local y axis: y axis in local cart. coords                         */
    vy[1]=0.0;
    vy[2]=1.0;
    vy[3]=0.0;
    cross(vy,vy_C,wektor);
    /* check the sign of the angle g                                      */
    if (wektor[3]>0.0)
    {
        sing=normv(wektor);
    }
    else
    {
        sing=-normv(wektor);
    }
    cosg=dot(vy,vy_C);
    /* in the stand-alone mode no need of rotating by g - the azimuth is  */
    /* initialy given in local frame (counterclockwise from geogr. North  */
    /* so we set g=0                                                      */
    if (corsika==0)
    {
        sing=0;
        cosg=1;
    }
    /*------------ end of computing angle g ------------------------------*/

    /*--------------------------------------------------------------------*/
    /* getting the trajectory unit vector (tn[]) and length (r_zero):     */
    /* and starting altitude                                              */
    /* r_loc: starting point of primary photon trajectory in local coords.*/
    r_loc[1]=r_0;
    r_loc[2]=the_loc;
    r_loc[3]=phi_loc;
    /* transform. r_loc to global coords -> r_glob_0                      */
    locsphC2glob(r_loc,r_glob_0,theta,phi,sing,cosg,sitec);
    /* r_glob (starting point of the primary photon trajectory in global  */
    /* coord.) is used, together with sitec - geographical position of    */
    /* the observatory - to find the unit vector tn[] of the photon (and  */
    /* subsequently preshower) trajectory                                 */
    /*                                                                    */
    /* unit vector of the preshower trajectory:                           */
    tn[1]=sitec[1]-r_glob_0[1];
    tn[2]=sitec[2]-r_glob_0[2];
    tn[3]=sitec[3]-r_glob_0[3];
    norm(tn,tn);
    r_curr_0=normv(r_glob_0);     /* the distance from the simulation   */
    /* starting point to Earth's center   */
    *r_zero=(r_curr_0-R_e)/km2cm; /* simulation starting distance above */
    /* sea level (in km)                  */

    /*--------------------------------------------------------------------*/
    /*                printing out preshower info                         */
    /*--------------------------------------------------------------------*/
    /* standard output */
    if (iprint==1) 
    {
	printf(" ---------------- PRESHOWER 2.0 ----------------- \n");
	printf(" B model: IGRF-11, year = %d \n", igrf_year);
	printf(" B[Gauss] in local cart. coord.: ");
	printf(" Bx = %5.3f, By = %5.3f, Bz = %5.3f \n",
               -B_site_loc[1],-B_site_loc[2],B_site_loc[3]);
	printf(" B[microTesla] in CORSIKA format: BX = %5.3f, BZ = %5.3f \n",
               100.*sqrt(B_site_loc[1]*B_site_loc[1] 
	       +B_site_loc[2]*B_site_loc[2]),
               -100.*B_site_loc[3]);
	if (corsika==1) 
	{
    	    printf(" g = 180 deg - declination: sin(g) = %f, cos(g) = %f\n",
                   sing,cosg);
	}
	if (corsika==1) 
	{
	    printf(" random generator: CORSIKA \n");
	} 
	else 
	{
    	    printf(" random generator: Numerical Recipes \n");
	}
	printf(" initial trajectory in local frame [deg]: ");
	printf(" zenith = %5.1f azimuth = %5.1f \n", the_deg, phi_deg);
	printf(" initial trajectory azimuth in CORSIKA [deg]: %5.1f \n",
                phi_deg+180);
	printf(" top of atmosphere from CORSIKA: %6.2f km \n", toa);
	printf(" site: longitude: %6.2f; latitude %6.2f \n", glong,glat);
	printf(" simulation starting altitude: %5.0f km a.s.l. \n",
    	       *r_zero);
	printf(" primary photon energy: E_gamma = %8.2e eV \n",E_gamma);
    }
    /* part_out.dat header */
    if (file_part_out==1) 
    {
	fprintf(file_1," ---------------- PRESHOWER 2.0 ----------------- \n");
	fprintf(file_1," B model: IGRF-11, year = %d \n", igrf_year);
	if (corsika==1) 
	{
	    fprintf(file_1," random generator: CORSIKA \n");
	} 
	else 
	{
	    fprintf(file_1," random generator: Numerical Recipes \n");
	}
	fprintf(file_1," B[Gauss] in local cart. coord.: ");
	fprintf(file_1," Bx = %5.3f, By = %5.3f, Bz = %5.3f \n",
    	       -B_site_loc[1],-B_site_loc[2],B_site_loc[3]);
	fprintf(file_1,
    	     " B[microTesla] in CORSIKA format: BX = %5.3f, BZ = %5.3f \n",
	    100.*sqrt(B_site_loc[1]*B_site_loc[1]+B_site_loc[2]*B_site_loc[2]),
            -100.*B_site_loc[3]);
	if (corsika==1) 
	{
	    fprintf(file_1,
            " g = 180 deg - declination: sin(g) = %f, cos(g) = %f\n",sing,cosg);
	}
	fprintf(file_1," initial trajectory in local frame [deg]: ");
	fprintf(file_1,
         " zenith = %5.1f azimuth = %5.1f \n", the_deg, phi_deg);
	fprintf(file_1,
            " initial trajectory azimuth in CORSIKA format [deg]: %5.1f \n",
            phi_deg+180);
	fprintf(file_1," top of atmosphere: %6.2f km \n", toa);
	fprintf(file_1," site: longitude: %6.2f; latitude %6.2f \n",
    	     glong,glat);
	fprintf(file_1," simulation starting altitude: %5.0f km a.s.l. \n",
    	     *r_zero);
	fprintf(file_1,
    	     " primary photon energy: E_gamma = %8.2e eV \n",E_gamma);
	fprintf(file_1," ------------------------------------------------ \n");
	fprintf(file_1," columns: \n");
	fprintf(file_1," 1. particle id (1-photon/2-positron/3-electron) \n");
	fprintf(file_1," 2. particle energy [eV] \n");
	fprintf(file_1," ------------------------------------------------ \n");
    }

    /* multirun.dat header */
    if (multirun_file==1)
    {
        fprintf(file_3," ---------------- PRESHOWER 2.0 ----------------- \n");
        fprintf(file_3," B model: IGRF-11, year = %d \n", igrf_year);
        if (corsika==1)
        {
            fprintf(file_3," random generator: CORSIKA \n");
        }
        else
        {
            fprintf(file_3," random generator: Numerical Recipes \n");
        }
        fprintf(file_3," B[Gauss] in local cart. coord.: ");
        fprintf(file_3," Bx = %5.3f, By = %5.3f, Bz = %5.3f \n",
                -B_site_loc[1],-B_site_loc[2],B_site_loc[3]);
        fprintf(file_3,
                " B[microTesla] in CORSIKA format: BX = %5.3f, BZ = %5.3f \n",
                100.*sqrt(B_site_loc[1]*B_site_loc[1]+B_site_loc[2]*B_site_loc[2]),
                -100.*B_site_loc[3]);
        if (corsika==1)
        {
            fprintf(file_3," g = 180 deg - declination: sin(g) = %f, cos(g) = %f\n",
                    sing,cosg);
        }
        fprintf(file_3," initial trajectory in local frame [deg]: ");
        fprintf(file_3," zenith = %5.1f azimuth = %5.1f \n",
                the_deg, phi_deg);
        fprintf(file_3,
                " initial trajectory azimuth in CORSIKA format [deg]: %5.1f \n",
                phi_deg+180);
        fprintf(file_3," top of atmosphere: %6.2f km \n", toa);
        fprintf(file_3," site: longitude: %6.2f; latitude %6.2f \n",
                glong,glat);
        fprintf(file_3," simulation starting altitude: %5.0f km a.s.l. \n",
                *r_zero);
        fprintf(file_3," primary photon energy: E_gamma = %8.2e eV \n",
                E_gamma);
        fprintf(file_3," ------------------------------------------------ \n");
        fprintf(file_3," columns: \n");
        fprintf(file_3," 1. run number \n");
        fprintf(file_3,
                " 2. altitude of the first gamma conversion in km a.s.l. \n");
        fprintf(file_3,
                " 3. total number of particles at the top of the atmosphere \n");
        fprintf(file_3," 4. number of photons \n");
        fprintf(file_3," 5. number of electrons \n");
        fprintf(file_3," 6. maximum electron energy \n");
        fprintf(file_3," 7. maximum photon energy \n");
        fprintf(file_3," 8. fraction of energy carried by electrons \n");
        fprintf(file_3,
                " 9. total energy carried by preshower particles = primary energy \n");
        fprintf(file_3," 10. total energy carried by electrons \n");
        fprintf(file_3," 11. total energy carried by photons \n");
        fprintf(file_3," ------------------------------------------------ \n");
    }

    /*------------------ INITIALIZATION ----------------------------------*/
    first_conv=0;
    i=&i_ok;
    k=&k_ok;
    di=&di_ok;
    dk=&dk_ok;
    x=0.0;
    order=0.0;
    n_surv=0;
    n_phot=1;
    n_eplus=0;
    n_eminus=0;
    seed_ok=-1000000;
    seed1=&seed_ok;

	/** @author:agata
	 * log_flag defines if some logging will be made 
	 * while counting max_alpha variables
	 */
    int log_flag =0;
    tr_par=0.0;

    E_e=E_gamma*eV2erg;
    /**
    *      @author: agata :preparing max values for veto algorithm
    */
    double alpha_max = 0.0;
    double alpha_max_brem = 0.0;

    alpha_max =  max_photon_alpha_rough( E_gamma,   igrf_year,  toa,
                                      r_step,  tn,  sitec,  chi_limit,  
				      log_flag);
    alpha_max*=1.1;
    
    alpha_max_brem =  max_electron_alpha_rough( E_e,   E_sg0,  igrf_year, toa,
                      r_step,  tn,  sitec,
                      tr_par, chi_limit, log_flag)   ;
    alpha_max_brem*=1.1;
    /*--------------------------------------------------------------------*/
    /*             PRESHOWER CREATION - main loop over # of runs          */
    /*--------------------------------------------------------------------*/
    for (mm=0; mm<nrun; mm++)
    {

        /* initialize                                                       */
        tr_par=0.0;                   /* trajectory parameter, current path */
        /* length travelled by the preshower  */
        n_phot=1;
        n_eplus=0;
        n_eminus=0;
        conversion=0;
        first_conv=0;
        jj=0;
        last_ii=0;                    /* highest not empty field in part_out */
        /* have to reload:                                                   */
        E_gamma=*pE_gamma/eV2GeV;     /* now E-gamma in eV                   */
        the_loc=*pthe_loc;
        phi_loc=*pphi_loc+PI;
        /* in corsika phi=0 means shower coming from magnetic SOUTH, in       */
        /* preshower all calculations are made for                            */
        /* phi=0 -> shower from mag. NORTH,                                   */
        /* have to convert phi from the CORSIKA format before starting      */
        /* simulations                                                         */
        /*--------------------------------------------------------------------*/
        if (iprint==1)
        {
            printf(" ------------------------------------------------ \n");
            if (nrun >= 2)
            {
                printf(" preshower run: %d \n", mm+1);
            }
        }
        /* initialization of preshower array part_out[];                      */
        /* first row: the initial photon (the input id is always 1)           */
        part_out[0][0]=*id;           /* particle id: 1 gamma, 2 e+, 3 e-   */
        part_out[0][1]=E_gamma;       /* Energy [eV]                        */
        part_out[0][2]=cos(the_loc);  /* particle zen. angle in local syst. */
        part_out[0][3]=phi_loc-PI;    /* particle azi. angle in local syst. */
        /* PI: have to return phi in the convention of CORSIKA                */
        part_out[0][4]=0.0;           /* t as above                         */
        part_out[0][5]=0.0;           /* x_loc: coordinates of the particle */
        /* in cm in the local system          */
        part_out[0][6]=0.0;           /* y_loc, related to the shower axis  */
        part_out[0][7]=tr_par;/**  @author agata: distance to the start point  */
        /* at the top of atm.                 */
        /* next rows: zeros for the beginning                                 */
        for (ii=1;ii<100000;ii++)
        {
            for (jj=0;jj<8;jj++)
            {
                part_out[ii][jj]=0.0;
            }
        }
        unit_dist=0.00000001;
        r_curr=r_curr_0;
        
        /*--------------------------------------------------------------------
        * PRESHOWER propagation: LOOP 1:  over current cascade particles 
        * repetition until none of the particles converts (photons) or
        * radiate (electrons).
        --------------------------------------------------------------------*/
        ii=0; 
        int veto_test_conver;
        while (part_out[ii][0]>0.0)
        {
            /* switch to different particle types                                 */
            switch ((int)part_out[ii][0])
            {
            case 1:
                /*--------------------------------------------------------------------
                *-----------------         PHOTON       -----------------------------
                *--------------------------------------------------------------------*/
                E_gamma=part_out[ii][1];  /* photon energy in eV                */

                if (E_gamma>E_conv_tresh )
                {
                    float dist = part_out[ii][7]; // @agata : getting actual particle distance

                    veto_res photon_res = photon_veto(E_gamma, igrf_year, toa, alpha_max,
                                                      tn, r_glob_0, r_curr_0, dist, chi_limit, seed1, corsika);
                    // @agata: counting new conversion distance
                    float r_veto =  photon_res.distance;

                    if(r_veto>0.0)
                    {
                        // B_tr and alpha value for point in which conversion occured
                        B_tr=photon_res.B_tr;
                        alpha_gamma=photon_res.alpha;

                        conversion=1;
                        ppf=ppfrac(seed1, corsika, B_tr, E_gamma*eV2erg);
                        /* update part_out                    */
                        part_out[ii][0]=2.0;   /* e+ replaces the photon             */
                        E_0=part_out[ii][1];   /* save primary energy to print later */
                        part_out[ii][1]*=ppf;  /* energy of e+                       */
                        part_out[ii][7]=r_veto; //@agata: updating current track distance
                        last_ii++;             /* one more particle appears, last_ii */

                        error_array_size(last_ii, part_out_size);

                        part_out[last_ii][0]=3.0; /* e- added at the end             */
                        part_out[last_ii][2]=cos(the_loc);
                        /* particle direction: zenith          */
                        part_out[last_ii][3]=phi_loc-PI;
                        /* particle direction: azimuth         */

                        //@agata: updating current track distance
                        // for second electron
                        part_out[last_ii][7]=r_veto;

                        part_out[last_ii][1]=part_out[ii][1]*(1/ppf-1);
                        /* particle energy                     */
                        if (iprint==1)
                        {
                            r_curr =tr_par2_r_curr(tn, r_glob_0, r_veto);
                            printf(" conv.: %5.0Lf km a.s.l: ",(r_curr-R_e)/km2cm);
                            printf(" Egamma = %8.2e -> Ee+ = %8.2e, Ee- = %8.2e \n",
                                   E_0, part_out[ii][1],part_out[last_ii][1]);
                        }
                        /* particle counters                   */
                        n_phot--;
                        n_eplus++;
                        n_eminus++;
                        if (first_conv==0)
                        {  /* saving distance for first pair prod */
                            r_curr =tr_par2_r_curr(tn, r_glob_0, r_veto);
                            *r_fpp=(r_curr-R_e)/km2cm;
                            first_conv=1;
                            /* for all  the next particles         */
                        }
                    }   // i : r veto found
                    else
                    { // if conversion doesn't occur the algorithm goes to the
                        // next particle
                        ii++;
                    }
                }                        /* end conv_thresh condition / 4 /     */
                else
                { // if energy isn't high enough the algorithm goes to the
                    // next particle
                    ii++;
                }
                break;
            case 2:
                /* no break: the same bremsstrahlung treatment done for cases 2 & 3   */
            case 3:
                /*--------------------------------------------------------------------*/
                /*-------       ELECTRON  (Sokolov, Ternov, 1986)     ----------------*/
                /*--------------------------------------------------------------------*/
                veto_test_conver=0;
                if(veto_test_conver!=1)
                {
                    E_e=part_out[ii][1]*eV2erg;
                    /* for calculations need E in ergs     */
                    double dist = part_out[ii][7];
                    veto_el_res el_res= electron_veto(E_e, E_sg0, igrf_year, toa,
                                                      tn, r_glob_0, r_curr_0, dist, chi_limit,
                                                      seed1, corsika, alpha_max_brem);
                    float r_veto =  el_res.distance;

                    if (r_veto>0.0)
                    {
                      // ctp120321 put that here to avoid copy of strange value when r_veto < 0
                      B_tr = el_res.B_tr;
                      prob_brem_total = el_res.prob_brem_total*unit_dist;

                        /* if a photon was emitted its energy is determined by the spectrum   */
                        /* I_brem:                                                            */
                        /*  the probability density function, I_brem*unit_dist*prob_factor,   */
                        /*  is integrated until randomly chosen fraction is reached,          */
                        /*  the energy E_ch at which integration end
                         parameter                 is assigned to newly    */
                        /*  emitted photon - most lines are the same as above in the calc. of */
                        /*  prob_brem_total                                                   */
                        n_phot++;              /* increase photon counter             */
                        rand_b=getrand(seed1,corsika);
                        /* get random number                   */
                        E_int=0.0;             /* initialize integrated energy        */
                        E_sg=E_sg0;
                        E_prev=E_sg;
                        E_slice=0.01;
                        frac=0.0;
                        while (frac<rand_b)
                        {
                            dE=E_sg*pow(10.0,E_slice)-E_sg;
                            E_sg=E_sg*pow(10.0,E_slice);
                            I_brem=brem(B_tr,E_e,E_sg);
                            prob_factor=dE/E_prev;
                            I_brem=I_brem*unit_dist*prob_factor;
                            E_int=E_int+I_brem;  /* integrated probab. dens. function  */
                            frac=E_int/prob_brem_total;
                            /* normalized probability             */
                            E_prev=E_sg;
                        }
                        E_ch=E_sg-dE;

                        /* the energy of the emitted particle: the last energy for which      */
                        /* integral E_int < rand_b*prob_brem_total                            */
                        /*      updating preshower particle array part_out[][]                */
                        part_out[ii][1]-=E_ch/eV2erg;
                        /* update primary e energy             */
                        part_out[ii][7]=r_veto;//@agata updating track distance
                        last_ii++;             /* one more particle appears           */
                        /* SECURITY: last_ii cannot exceed the output array size part_out_size*/
                        error_array_size(last_ii, part_out_size);

                        part_out[last_ii][0]=1.0;
                        /* a photon added at the end           */

                        part_out[last_ii][2]=cos(the_loc);
                        /* particle direct.: cosine            */
                        part_out[last_ii][3]=phi_loc-PI;
                        /* particle direct.: local azimuth     */
                        part_out[last_ii][1]=E_ch/eV2erg;
                        /* photon energy                       */


                        part_out[last_ii][7]=r_veto;//@agata updating track distance
                       
                        if (iprint==19)
                        { // @agata additional info about bremsstrahlung
                            r_curr =tr_par2_r_curr(tn, r_glob_0, r_veto);
                            printf(" brem.: %5.0Lf km a.s.l: ",(r_curr-R_e)/km2cm);
                            printf(" E_el = %8.2e -> Eel_new = %8.2e, E_gamma_new = %8.2e \n",
                                   E_e/(double)eV2erg, part_out[ii][1],part_out[last_ii][1]);

                        }

                    }
                    else
                    { //if bremsstrahlung doesn't occur the algorithm
                        //  goes to the next particle
                        ii++;
                    }                        /* end 'if a photon was emitted'       */
                    break;
                }
                else
                { //if bremsstrahlung doesn't occur the algorithm
                    //  goes to the next particle
                    ii++;
                }  // if test_conversion condition
            };
            /* end 2 swich     ************************************                    */

        }                          /* end 1 while for each particle       */
        if ((conversion==0))
        {
            /* photon didn't convert before entering the atmosphere */
            n_surv++;          /* count surviving photons (in         */
            /* multirun  mode)                    */
            if (iprint==1)
            {
                printf(" photon SURVIVED ! \n");
            }
            *r_fpp=toa;        /* first interaction at the top of     */
            /* atmosphere (toa)                    */
        }

        *N_part_out=last_ii+1;     /* # of entries in part_out array      */
        Ee_max=0.0;                  /* reset maximum electron energy in run*/
        Eg_max=0.0;                  /* reset maximum gamma energy in run   */
        Etot_g=0.0;                  /* reset total energy in gammas        */
        Etot_e=0.0;                  /* reset total energy in e+/-          */
        jj=0;
        /* store preshower particles data in a file, compute energy sums      */
        /* (gammas and electrons) and maximum energies                        */
        while (part_out[jj][0]>0.0)
        {
            for (ii=0;ii<2;ii++)
            {
                if (file_part_out==1)
                {
                    fprintf(file_1,"%8.2e ",part_out[jj][ii]);
                }
                if (ii==1)
                {             /* energy                              */
                    if (part_out[jj][0]==1)
                    {
                        Etot_g+=part_out[jj][ii];
                        /* increase total energy in gammas     */
                        if (Eg_max<part_out[jj][ii])
                        {
                            Eg_max=part_out[jj][ii];
                            /* update maximum gamma energy         */
                        }
                    }
                    if (part_out[jj][0]>1)
                    {
                        Etot_e+=part_out[jj][ii];
                        /* increase total energy in electrons  */
                        if (Ee_max<part_out[jj][ii])
                        {
                            Ee_max=part_out[jj][ii];
                            //  printf("%g  , e %g\n", Ee_max=part_out[jj][ii], Ee_max);
                            /* update maximum electron energy      */
                        }
                    }
                }
            }
            if (file_part_out==1)
            {
                fprintf(file_1,"\n");
            }
            part_out[jj][1]*=eV2GeV;   /* output energy in GeVs               */
            jj++;
        }
        if (iprint==1)
        {
            printf(" summary (energies in eV): \n");
            printf(" Etot_g = %8.2e, Etot_e = %8.2e,", Etot_g,Etot_e);
            printf(" Eg_max = %8.2e, Ee_max = %8.2e\n", Eg_max,Ee_max);
        }
        /* print to multi-run file                                            */
        if (multirun_file==1)
        {
            fprintf(file_3,
                    "%3d %8.2f %4d %4d %2d %8.2e %8.2e %6.4f %8.2e %8.2e %8.2e \n",
                    mm+1,*r_fpp,*N_part_out,n_phot,n_eplus+n_eminus,Ee_max,
                    Eg_max,Etot_e/(Etot_g+Etot_e),Etot_g+Etot_e,Etot_e,Etot_g);
        }
        if (iprint==1)
        {
            printf(" n_part = %d; n_phot = %d; n_e+ = %d; n_e- = %d \n",
                   *N_part_out,n_phot,n_eplus,n_eminus);
        }


    }                              /* end run loop               */
    /*--------------------------------------------------------------------*/
    /*      PRESHOWER CREATION - end of main loop over # of runs          */
    /*--------------------------------------------------------------------*/
    /* in multi-run mode determine percentage of surviving (those which   */
    /* didn't undergo pair-production) photons                            */
    if (iprint==1)
    {
        if (nrun >= 2)
        {
            surv_percent=(float)n_surv/(float)nrun;
            printf(" ------------------------------------------------ \n");
            printf(" # of surviving photons = %d; frac. = %f \n",
                   n_surv,surv_percent);
        }
        printf(" ------------ END of PRESHOWER output ----------- \n");
    }

    if (file_part_out==1)
    {
        fclose(file_1)
        ;
    }
    if (multirun_file==1)
    {
        fclose(file_3)
        ;
    }
    fflush(stdout);
    fflush(stderr);

}
