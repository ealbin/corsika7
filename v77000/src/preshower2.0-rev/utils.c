#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include <sys/types.h>

/**
 * @file utils.c
 * @author agata
 * @brief package with function and structures definitions which are used in program
 * but cannot be qualified to any special package
 */

////////////////////////

/**
 * Alpha function counts alpha value from some given arguments.
 * It returns  also B_tr value which is used
 * in pararrel to the alpha value so it doesn't need to be counted twice.
 * this function is used in photon_veto() and countMaxAlpha() functions.
 */
alpha_res alpha(  float r_glob[4],float igrf_year,float tn[4], double E_gamma, double chi_limit)
{
    float r_glob_sph[4];
    double Bglob[4], cernbess[101];
    double T_Y, kk;
    alpha_res res;

    car2sph(r_glob,r_glob_sph);
    b_igrf(igrf_year,r_glob_sph[1],r_glob_sph[2],r_glob_sph[3],Bglob);
    double B_tr=B_transverse(Bglob,tn);
    // double B_tr= b_tr_temp;

    double Y_gamma=E_gamma*eV2erg*B_tr*0.5/(m_e*c_light*c_light*B_cr);
    if (Y_gamma<chi_limit)
    {
        T_Y=0.46*exp(-4/(3*Y_gamma)); /* Erber 3.3c                  */
    }
    else
    {
        double bess_arg=2/(3*Y_gamma);
        int bn1=1;
        int bd=3;
        int bnl=0;
        dbska_(&bess_arg,&bn1,&bd,&bnl,cernbess);
        kk=cernbess[0];

        T_Y=0.16*kk*kk*(1/Y_gamma);   /* pair production function,   */
    }
    double alpha_gamma=0.5*alpha_fs*m_e*c_light*B_tr*T_Y/(h_bar*B_cr);
    res.alpha = alpha_gamma;
    res.B_tr = B_tr;

    return res;
}

brem_res prob_brem_total(double E_sg0, double E_e, float r_glob[4],
                         float igrf_year, float tn[4])
{
    float r_glob_sph[4];
    double Bglob[4];
    // counting B_tr
    car2sph(r_glob,r_glob_sph);
    b_igrf(igrf_year,r_glob_sph[1],r_glob_sph[2],r_glob_sph[3],Bglob);
    double B_tr=B_transverse(Bglob,tn);
    // double B_tr= b_tr_temp;

    brem_res result;

    double E_sg=E_sg0;              /* sg = secondary gamma, particles     */
    /* below energy E_sg0 are neglected    */
    double E_prev=E_sg;             /* initial 'previous' energy           */
    /* calculate the total probability of emitting a bremsstrahl. photon  */
    /*  (sg) - prob_brem_total. - do a loop in log steps changing brem.   */
    /*  photon energy                                                     */
    double prob_brem_total=0.0;     /* initialize total probability        */
    float E_slice=0.01;            /* logarithmic step of integration the */
    /* probability function                */
    // printf(" E_sg %g, E_slice %g, E_e %g \n", E_sg, E_slice,E_e );

    while (E_sg*pow(10.0,E_slice)<E_e)
    {
        /* loop over brem. photon energies     */
        E_sg=E_sg*pow(10.0,E_slice);

        /* current brem. phot. energy          */
        /* spectral distribution of radiated energy: I_brem is the energy     */
        /*   emitted into photons of energies between E_sg and E_sg+dE,       */
        /*   per unit distance:                                               */
        /*   I_brem=(E_prev*dN)/(dE*unit_dist); dN, being the number of       */
        /*   emitted photons, below is used as probability of emitting a      */
        /*   photon - prob_g, small simulation step provides dN << 1          */
        double I_brem=brem(B_tr,E_e,E_sg);
        double dE=E_sg-E_prev;        /* in ergs                             */

        double prob_factor=dE/E_prev; /* auxilliary                          */
        double prob_g=I_brem*prob_factor;
        /* probability of emitting photon of   */
        /*energy in interval (E_sg, E_sg+dE)   */
        E_prev=E_sg;           /* save E_sg for the next round        */
        prob_brem_total=prob_brem_total+prob_g;
        /* total probability emitting a photon */
    }
    result.prob_brem_total= prob_brem_total;
    result.B_tr=B_tr;

    return result;
}



double brem(float B_tr, double E_e, double E_sg)

/* bremsstrahlung spectrum according to Sokolov                       */
/* B_tr: transverse magnetic field                                    */
/* E_e: electron energy                                               */
/* E_sg: energy of the emitted bremsstrahlung photon                  */
{
    double *i,*k,*di,*dk,i_ok,k_ok,di_ok,dk_ok,k_2_3_y,dy_dhv,
    kap_y,I_brem,ksi,y,fy,p3,Wcl;
    double cernbess[101];
    int bn2,bd,bnl;

    i=&i_ok;
    k=&k_ok;
    di=&di_ok;
    dk=&dk_ok;
    ksi=1.5*(B_tr/B_cr)*(E_e*m_e_erg_inv);
    y=E_sg/(E_e*ksi*(1.0-E_sg/E_e)); /* argument of Bessel function     */
    //  if (y>100) {                     /* assymtotic expansion            */
    //                                   /* (from Numerical Recipes)        */
    //    b1=PI/sqrt(2*PI*y);
    //    b2=b1*exp(-y);
    //    k_2_3_y=b2;
    //  } else {

    /* DBSKA: CERNLIB Bessel function of order 2/3                        */
    bn2=2;
    bd=3;
    bnl=0;
    dbska_(&y,&bn2,&bd,&bnl,cernbess);
    k_2_3_y=cernbess[0];

    //    bessik(y,0.66666667,&i_ok,&k_ok,&di_ok,&dk_ok);
    //    k_2_3_y=*k;
    //  }
    p3=1.0+y*ksi;
    kap_y=kappa(y);               /* Erber's kappa without multip. by y */
    p3=pow(p3,3.0);
    fy=(kap_y+k_2_3_y*ksi*y*ksi*y/(1.+ksi*y))*(9.*sqrt(3.)/(8*PI))*y/p3;
    Wcl=B_tr*B_tr*E_e*E_e*r0*r0*m_e_erg_inv*m_e_erg_inv*2.0/3.0;
    dy_dhv=E_e/(ksi*(E_e-E_sg)*(E_e-E_sg));
    I_brem=fy*Wcl*dy_dhv;
    return I_brem;
}


float getrand(long *seed1, int corsika)
/*--------------------------------------------------------------------*/
/*  get random number from (0,1)                                      */
/* choose random number generator depending on the working mode:      */
/* PRESHOWER-CORSIKA: rnmd(0) - the same generator as in CORSIKA      */
/* PRESHOWER stand-alone - generator from Numerical Recipes           */
/*--------------------------------------------------------------------*/
{
    float rand;
    rand = rndm(0);             // in CORSIKA                       
/*    if (corsika==1)
    {              // random generator: the one used   
        rand = rndm(0);             // in CORSIKA
    }
    else
    {                       // random generator from Numerical Recipes
        rand=ran2(seed1);
    }
*/
    return rand;
}
/*--------------------------------------------------------------------*/

float kappa(float x)

/*--------------------------------------------------------------------*/
/* Returns bremsstrahlung auxiliary function according to Erber       */
/* numerical values  (Table I, p. 632) + 0.003; 0.005; 3.0 from       */
/* calculations of Piotr Homola, the function is slightly modified in */
/* order to be used in computation of bremsstrahlung expression given */
/* by Sokolov.                                                        */
/* interpolation between points: linear                               */
/*--------------------------------------------------------------------*/
{
    int i;
    float kapp,x1,x2,y1,y2;
    float k[31][31]={
                        {0.001,0.003,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
                         0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
                         1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0},
                        {0.213,0.308,0.361,0.445,0.547,0.614,0.663,0.701,0.733,0.760,
                         0.782,0.801,0.818,0.903,0.924,0.905,0.871,0.831, 0.788,0.742,
                         0.696,0.651,0.3,0.129,0.053,0.0214,0.00845,0.00338,0.00129,
                         0.00049,0.00019}
                    };

    if (x<k[0][0])
    {               /* k[0][0]=.001,approximation: x<<1 */
        kapp=2.1495*pow(x,1.0/3.0)-1.8138*x;
        /* (Erber A27)                      */
    }
    else if (x>=k[0][0] && x<=k[0][30])
    {
        /* k[0][30]=10.0                    */
        for (i=0;i<30;i++)
        {
            if ((x>=k[0][i]) && (x<=k[0][i+1]))
            {
                /* x between points i, i+1, linear  */
                /* interpolation                    */
                x1=k[0][i];
                x2=k[0][i+1];
                y1=k[1][i];
                y2=k[1][i+1];
                kapp=y1+(x-x1)*(y2-y1)/(x2-x1);
            }
        }
    }
    else if (x>k[0][30])
    {      /* k[0][30]=10.0, approximation: x>>1*/
        kapp=1.2533*sqrt(x)*exp(-x);
        /* Erber A30                         */
    };
    /* the below is different from kappa used for Erber bremsstrahlung    */
    /* expression, in order to have the Sokolov's one we have to divide   */
    /* original Erber's  kappa by its argument                            */
    kapp/=x;
    return kapp;
}
/*--------------------------------------------------------------------*/

float B_transverse(double B[4], float tn[4])

/* compute B component perpendicular to the preshower trajectory      */
{
    float B_mod,B_ii,B_tr;
    B_mod=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]); /* B value              */
    B_ii=tn[1]*B[1]+tn[2]*B[2]+tn[3]*B[3];     /* B paralell (tn - unit*/
    /*   trajectory vector) */
    B_tr=sqrt(B_mod*B_mod-B_ii*B_ii);          /* B transverse         */
    return B_tr;
}
/*--------------------------------------------------------------------*/

void cross(float a[4], float b[4], float c[4])

/*    calculates vector cross product c = a x b                       */
{
    c[1]=a[2]*b[3]-a[3]*b[2];
    c[2]=a[3]*b[1]-a[1]*b[3];
    c[3]=a[1]*b[2]-a[2]*b[1];
}
/*--------------------------------------------------------------------*/

float dot(float a[4], float b[4])

/* returns dot product of two vectors                                 */
{
    return a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
}
/*--------------------------------------------------------------------*/

void norm(float a[4], float n[4])

/* returns normalized vector                                          */
{
    float av;
    av=sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3]);
    n[1]=a[1]/av;
    n[2]=a[2]/av;
    n[3]=a[3]/av;
}
/*--------------------------------------------------------------------*/

float normv(float a[4])

/* returns norm of a vector                                           */
{
    return sqrt(a[1]*a[1]+a[2]*a[2]+a[3]*a[3]);
}
/*--------------------------------------------------------------------*/

void glob2loc(float glob[4], float loc[4], float theta, float phi)

/* converts global cartesian coords. to local cartesian               */
/* 3dim rotation                                                      */
{
    loc[1]=glob[1]*cos(theta)*cos(phi)
           +glob[2]*cos(theta)*sin(phi)
           -glob[3]*sin(theta);
    loc[2]=glob[2]*cos(phi)-glob[1]*sin(phi);
    loc[3]=glob[1]*sin(theta)*cos(phi)
           +glob[2]*sin(theta)*sin(phi)
           +glob[3]*cos(theta);
}
/*--------------------------------------------------------------------*/

void locsphC2glob(float locsph[4], float glob[4], float theta,
                  float phi, float sing, float cosg, float y[4])

/* local spherical (CORSIKA frame) to global cart. coords.            */
{
    float w[4],loc[4];
    sph2car(locsph,loc);              /* local spher. to local cart.   */
    w[1]=cosg*loc[1]-sing*loc[2];     /* rotation by -g                */
    w[2]=cosg*loc[2]+sing*loc[1];
    w[3]=loc[3];
    /* rotation by the geographic position angles                         */
    glob[1]=w[1]*cos(theta)*cos(phi)
            -w[2]*sin(phi)
            +w[3]*sin(theta)*cos(phi);
    glob[2]=w[1]*cos(theta)*sin(phi)
            +w[2]*cos(phi)
            +w[3]*sin(theta)*sin(phi);
    glob[3]=w[3]*cos(theta)-w[1]*sin(theta);
    glob[1]+=y[1];                    /* y[] are the observatory site  */
    glob[2]+=y[2];                    /* global carthesian coords.     */
    glob[3]+=y[3];
}
/*--------------------------------------------------------------------*/

void sph2car(float r_loc_sph[], float r_loc_cart[])

/* converts spherical to cartesian coords.                            */
{
    r_loc_cart[1]=r_loc_sph[1]*sin(r_loc_sph[2])*cos(r_loc_sph[3]);
    r_loc_cart[2]=r_loc_sph[1]*sin(r_loc_sph[2])*sin(r_loc_sph[3]);
    r_loc_cart[3]=r_loc_sph[1]*cos(r_loc_sph[2]);
}

/*--------------------------------------------------------------------*/

void car2sph(float car[4], float sph[4])

/* inverse to sph2car                                                 */
{
    double sq;
    sq=car[1]*car[1]+car[2]*car[2];
    sph[1]=(float)sqrt(sq+car[3]*car[3]);
    if (sq==0.0)
    {
        sph[3]=0.0;
        if (car[3]<0.0)
        {
            sph[2]=PI;
        }
        else
        {
            sph[2]=0.0;
        }
        return;
    }
    sq=sqrt(sq);
    sph[3]=(float)atan2(car[2],car[1]);
    sph[2]=(float)atan2(sq,car[3]);
    if (sph[3]<0.0)
    {
        sph[3]+=2*PI;
    }
    return;
}
/*--------------------------------------------------------------------*/

void bsph2car(float colat, float longit, float bsph[4], double bcar[4])

/*  calculates cartesian field components from spherical ones         */
/*  INPUT:   colat,longit - spherical angles of the point in radians  */
/*           bsph[4] -  spherical components of the field             */
/*           bcar[4] - cartesian components of the field              */
/* used in b_igrf only                                                */
{
    float s,c,sf,cf,be;
    s=sin(colat);
    c=cos(colat);
    sf=sin(longit);
    cf=cos(longit);
    be=bsph[1]*s+bsph[2]*c;
    bcar[1]=be*cf-bsph[3]*sf;
    bcar[2]=be*sf+bsph[3]*cf;
    bcar[3]=bsph[1]*c-bsph[2]*s;
    return;
}
/*--------------------------------------------------------------------*/

void b_igrf(int igrf_year, float igrf_r, float igrf_theta,
            float igrf_phi, double bcar[4])

/* calling a fortran procedure IGRF by                                */
/* N.A.Tsyganenko, USRA, Greenbelt, USA                               */
/* IGRF - geomagnetic field model, www.ngdc.noaa.gov                  */
{
    float br,btheta,bphi;
    float bsph[4];
    int n,iy;
    iy=igrf_year;                 /* year of the simulation            */
    n=13;                         /* highest order of spher. harmonics */
    igrf_r=igrf_r/R_e;            /* r in units of R_e                 */
    /* call igrf routine (external, in fortran)                           */
    igrf_(&iy,&n,&igrf_r,&igrf_theta,&igrf_phi,&br,&btheta,&bphi);
    bsph[1]=(1e-5)*br;            /* change units to G                 */
    bsph[2]=(1e-5)*btheta;
    bsph[3]=(1e-5)*bphi;
    /* need to transform from  B in spherical coords. to cartesian        */
    bsph2car(igrf_theta,igrf_phi,bsph,bcar);
    /* returns bcar[]: magnetic field in global cart. coords.             */
}
/*--------------------------------------------------------------------*/

double ppp(float efrac, float B_tr, float E_g0)

/* returns pair production probability according to Daugherty'83      */
/* (equivalent to Erber'66 as shown in Tsai, Erber'74, this function  */
/* is used only  within ppfrac (computation of efrac) while in the    */
/* main program the probability is computed using expression Erber'66 */
/* E_g0 is the primary gamma energy                                   */
/* B_tr is the transverse B component                                 */
/* efrac is the fractional energy of one electron                     */
{
    double *i,*k,*di,*dk,i_ok,k_ok,di_ok,dk_ok,k_2_3_y,ppprob;
    double y,chi,cernbess[101];
    int bn2,bd,bnl;

    i=&i_ok;
    k=&k_ok;
    di=&di_ok;
    dk=&dk_ok;
    chi=0.5*(B_tr/B_cr)*(E_g0*m_e_erg_inv);
    /* same as in Erber          */
    y=1.0/(3.0*chi*efrac*(1.0-efrac));     /* bessel function argument  */
    //  if (y>100) {                           /* assymtotic expansion      */
    //                                         /* (Numerical Recipes)       */
    //    b1=PI/sqrt(2*PI*y);
    //    b2=b1*exp(-y);
    //    k_2_3_y=b2;
    //  } else {

    /* DBSKA: CERNLIB Bessel function of order 2/3                        */
    bn2=2;
    bd=3;
    bnl=0;
    dbska_(&y,&bn2,&bd,&bnl,cernbess);
    k_2_3_y=cernbess[0];

    //    bessik(y,0.66666667,&i_ok,&k_ok,&di_ok,&dk_ok);
    //    k_2_3_y=*k;
    //  }

    ppprob = k_2_3_y * m_e * c_light * alpha_fs * B_tr * sqrt(3.0)
             * (2.0+efrac*(1.0-efrac))
             /(h_bar*B_cr*9.0*PI*chi*efrac*(1.0-efrac));
    /* Eq. 19                    */
    return ppprob;
}
/*--------------------------------------------------------------------*/

float ppfrac(long *s, int corsika, float B_tr, float E_g0)

// determines the fractional energy of a pair member,                 */
/* Daugherty, Astroph. J. 273, 1983                                   */
{
    double ppprob,ppprob_max;
    float estep,efrac_i,r1,r2;
    int nestep,ii;
    estep=0.001;                   /* step in electron/positron energy  */
    nestep=ceil(0.5/estep)+1;      /* # of steps                        */
    efrac_i=0.0;
    ppprob_max=0;

    for(ii=1;ii<=nestep;ii++)
    {
        efrac_i=ii*estep;
        ppprob=ppp(efrac_i,B_tr,E_g0);
        if (ppprob_max<ppprob)
        {
            ppprob_max=ppprob;
        }
    }                              /* end of seeking maximum - needed   */
    /* for normalization                 */
    r1=getrand(s,corsika);
    r2=getrand(s,corsika);
    ppprob=ppp(r1,B_tr,E_g0)/ppprob_max;

    /* get fractional energy of a pair member according to distribution   */
    while (r2>ppprob)
    {
        r1=getrand(s,corsika);
        r2=getrand(s,corsika);
        ppprob=ppp(r1,B_tr,E_g0)/ppprob_max;
    }
    return r1;
}

extern void rmmard_(double *, int *, int *);
/* external use of the RMMARD random number generator of CORSIKA      */
/* used only in CORSIKA mode                                          */                                       

static double rndm()
/* a function to call external rmmard_;                               */
/* returns random number from (0,1);  used only in CORSIKA mode       */
{
    static double rtmp[10];
    int num = 1;
    int seq = 2;                /* Use sequence number 2 of RMMARD     */
    rmmard_(rtmp,&num,&seq);    /* Call Fortran subroutine             */
    return rtmp[0];
}


/**
 * error_array_size is a helpful function that checks if number of particles
 * in part_out array doesn't reach the permitted number (size)
 * If it is it show some info on the screen en exits program
 * It is used in photon and electron part of simulation
 */
void error_array_size(int last_ii, int part_out_size)
{
    if (last_ii>part_out_size)
    {
        printf(" -----------------------------------\n");
        printf(" PRESHOWER ABORT: # of particles (%d) >", last_ii);
        printf("  the output array size \n");
        exit(1);             /* exit the program                   */
    }
}

