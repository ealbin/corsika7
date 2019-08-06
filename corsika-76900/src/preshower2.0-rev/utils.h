#ifndef UTILS_H_
#define UTILS_H_
/**
 * @file utils.h
 * @author agata
 * @brief package with function and structures definitions which are used in program
 * but cannot be qualified to any special package
 */
 
/////////////////////////////////////////////
// CONSTANT DEFINITIONS
/////////////////////////////////////////////
#define R_e 637131500.0l       // earth radius in cm
#define km2cm 1.0e+5l           /* km to cm                           */
#define PI 3.141592653589793l
/*--------------------------------------------------------------------*/

#define deg2rad PI/180.0;       /* degrees to radians                 */
#define rad2deg 180.0/PI;       /* radians to degrees                 */
#define c_light 29979245800.0l  /* speed of light in vacuum in cm/s   */
#define m_e 9.1093897e-28l      /* electron rest mass in g            */
#define m_e_erg_inv 1.22142e+6l /* inverse electron rest mass in ergs */
#define alpha_fs 0.00729735308l /* fine structure constant            */
/* (dimensionless)                    */
#define h_bar 1.054557266e-27l  /* Planck constant x 1/2PI (h-bar)    */
/* in erg*s (1J=10^7erg)              */
#define e_SI 1.60217733e-19l    /* elementary charge in C             */
#define e_CGS 4.803e-10l        /* elementary charge in CGS units;    */
/* 1C=2.998e+9 CGS units of charge    */
#define B_cr 4.414e+13l         /* critical magnetic field acc. to    */
/* Erber, Rev.Mod.Phys. 38 (1966) 626 */
/* B_cr=m_e^2*c^3/(e*h_bar) in G      */
#define eV2erg 1.60219e-12l     /* electronvolts to ergs              */
#define eV2GeV 1.0e-9l          /* eV to GeV                          */
#define GeV2erg 1.60219e-3l     /* GeV to ergs                        */
#define r0 2.8179409238e-13l    /* classical electron radius [cm]     */
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
///////////////////////////////////////////////
// STRUCTURES DEFINITIONS
////////////////////////////////////////////

/**
 * config structure which is used while calling preshw()  function
 * it may contain some parameters which will be used in main program
 * it can be easily changed  by adding some elements without changing 
 * the preshw() argument list
 */
typedef struct settings
{
    FILE *out;
    double phi_temp;
    double b_tr_temp;
    int step_convert;
    int step_not_convert;
    int n_surv;
    double p_conv;
    double param;

}
config;

/**
 * veto_res structure is a container for variables that are counted in photon_veto()
 * function.
 * distance - is an interval that was found by veto algorithm at which conversion occurs.
 *           If it is equal to -1  then interval wasn't found.
 * B_tr  - is the transverse magnetic field counted for the found distance
 * alpha - alpha function value computed for the found distance
 */
typedef struct veto_result
{
    float distance;
    double B_tr;
    double alpha ;
}
veto_res;

/**
 * alpha_res structure is a container for variables that are computed in alpha() function.
 */
typedef struct alpha_result
{
	/**
	 * B_tr  - is the transverse magnetic field computed for distance given in alpha() 
	 */
    double B_tr;
	/**
	 * alpha function value computed for the distance given as an argument in alpha().
	 */
    double alpha ;
}
alpha_res;
////////////////////////////////////////////////

/**
 * veto_res structure is a container for variables that are computed in electron_veto()
 * function.
 * distance - is an interval that was found by veto algorithm at which bremsstrahlung occurs.
 *           If it is equal to -1  then interval wasn't found.
 * B_tr  - is the magnetic field computed for the distance found
 * alpha - alpha function value computed for the distance found
 */
typedef struct veto_electron_result
{
    float distance;
    double B_tr;
    double prob_brem_total;
}
veto_el_res;

/**
 * brem_res structure is a container for variables that are computed in 
 * prob_brem_total() function.
 * B_tr  - is the magnetic field computed for the distance given in alpha() 
 * prob_brem_total - alpha function value computed for the distance given as an argument in alpha.
 * IMPORTANT : prob_brem_total is not the same as in the original version of program
 * it must be multiplied by unit_dist to have the same role in main program.
 */
typedef struct brem_result
{
    double B_tr;
    double prob_brem_total ;
}
brem_res;


double alpha_from_tr_par(float tr_par,  float r_glob_0[4],float igrf_year,float tn[4],
                         double E_gamma, double chi_limit);


double prob_brem_total_from_tr_par(double tr_par, float r_glob_0[4], double E_sg0,
                                   double E_e,  float igrf_year, float tn[4]);

alpha_res alpha(  float r_glob[4],float igrf_year,float tn[4],
                  double E_gamma, double chi_limit);

brem_res prob_brem_total(double E_sg0, double E_e, float r_glob[4],
                         float igrf_year, float tn[4]);

void error_array_size(int last_ii, int part_out_size);

///////////////////////////////////////////////////////////////////
// FUNCTION DEFINITIONS
////////////////////////////////////////////

//extern void bessik(float x, float xnu, double *ri, double *rk,
//                   double *rip, double *rkp);
//                      /* cumputes modified Bessel functions */
//                                /* of fractional order (Num. Recipes) */
float getrand(long *, int);     /* get random number from (0,1)       */
/* choosing an appropriate generator  */
float kappa(float);             /* auxiliary function used in brems-  */
/* strahlung calculation procedure    */
double brem(float, double, double);
/* computes bremsstrahlung radiation  */
void conv(double part_out[50000][8],int id);
/* computes gamma conversion probabil.*/
float B_transverse(double B[4], float tn[4]);
/* getting B component perpendicular  */
/* to the gamma trajectory            */
extern float ran2(long *);      /* random number generator used only  */
/* in stand-alone mode (Num. Recip.)  */
void cross(float a[4], float b[4], float c[4]);
/* cross product of vectors           */
float dot(float a[4], float b[4]);
/* dot product of vectros             */
void norm(float a[4], float n[4]);
/* returns normalized vector          */
float normv(float a[4]);        /* returns vector's length            */
void glob2loc(float glob[4], float loc[4], float theta, float phi);
/* conv. global cart. coord. to local */
void sph2car(float r_s[], float r_c[]);
/* converts spherical to cart. coords.*/
void car2sph(float car[4], float sph[4]);
/* inverse to sph2car                 */
void locsphC2glob(float locsph[4], float glob[4], float theta, float phi,
                  float sing, float cosg, float sitec[4]);
/* local spherical (CORSIKA frame) to */
/* global cart. coords.               */
double ppp(float efrac, float B_tr, float E_g0);
/* returns gamma convers. probability,*/
/*  used within ppfrac                */
float ppfrac(long *s, int corsika, float B_tr, float E_g0);
/* returns gamma convers. probability */
/* using ppp()                        */
extern void igrf_(int *,int *,float *,float *,float *,float *,float *,
                      float *);     /* external procedure by Tsyganenko   */
/* computes magnetic field according  */
/* to the IGRF model, used within     */
/* b_igrf                             */
extern void dbska_(double *,int *,int *,int *,double cernbess[101]);
/* external CERNlib bessel function  */
void b_igrf(int igrf_year, float igrf_r, float igrf_theta,
            float igrf_phi, double bcar[4]);
/* returns magn. field in sph. coords.*/
/*  accord. IGRF model using igrf_()  */
void bsph2car(float colat, float longit, float bsph[4], double bcar[4]);
/* converts b in spherical coords.    */
/* into B in cart. coords.,           */
/* used in b_igrf only                */
static double rndm();           /* calls random generator of CORSIKA  */


extern void rmmard_(double *, int *, int *);
/* external use of the RMMARD random number generator of CORSIKA      */
/* used only in CORSIKA mode*/

#endif /*VETO_H_*/
