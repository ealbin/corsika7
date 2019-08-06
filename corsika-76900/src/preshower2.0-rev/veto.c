#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "math.h"
#include "veto.h"

/**
 * @file veto.c
 * @author agata
 * @brief package for functions which are used by veto algorithm 
 * 
 */
 
double b_tr_temp;
/**
 * photon_veto function is used to find a distance at which photon
 * conversion occurs. Function uses veto algoritm (described e.g. here:
 * http://cepa.fnal.gov/psm/simulation/mcgen/lund/v6_214/pythia6206/node31.html ).
 * 
 * When the algorithm finds the distance of conversion it returns  
 * also alpha value and B_tr value that were previously computed at 
 * this distance. They are then easily used in the main preshower simulation.
 * 
 * If the algoritm doesn't find the distance it returns the same structure
 * but the distance value is set to -1.0 which can be recognized from the upper 
 * level of program.
 * 
 * function can write some output to the file called logVeto.dat, but to achieve
 * this, line with fprintf(... must be uncommented and all recompiled.
 */
veto_res photon_veto(double E_gamma,  float igrf_year,float toa,  float alpha_max,
                     float tn[4], float r_glob_0[4],
                     float r_curr_0, float tr_par, double chi_limit,
                     long *seed1, int corsika)
{
    float tr, rc, tr_prev, rnd, last_step;
    float r_glob[4];
    int   doLoop, i;
    veto_res result;
//    FILE *logfile;
//    logfile=fopen("logVeto.dat", "w");
    i=0; // this variable is used only in log mode - it is written to file

    doLoop=1;

    tr = tr_par; // variable tr will be used as tr_par in following code

    // initialize rc
    r_glob[1]=r_glob_0[1] + tr*tn[1];
    r_glob[2]=r_glob_0[2] + tr*tn[2];
    r_glob[3]=r_glob_0[3] + tr*tn[3];
    /* current distance of the simulation point to the Earth's center     */
    rc=normv(r_glob); // rc is the same as r_curr in older algorithm

    while(doLoop==1)
    {
        tr_prev=tr;

//        fprintf(logfile, "%d %g %g\n",i,tr, tr-tr_prev );
        i++;
        rnd = getrand(seed1,corsika);    /* get random number           */

        tr = (alpha_max*tr_prev - log(rnd) ) /alpha_max;

        r_glob[1]=r_glob_0[1] + tr*tn[1];
        r_glob[2]=r_glob_0[2] + tr*tn[2];
        r_glob[3]=r_glob_0[3] + tr*tn[3];
        /* current distance of the simulation point to the Earth's center     */
        rc=normv(r_glob);
        last_step = tr-tr_prev;

        // dummy way to prevent going on the other site of the Earth
        if (r_curr_0<tr)
        {
            result.distance=-1.0;
            doLoop=0;
        }
        else if(rc>R_e+toa*km2cm)
        {
            rnd = getrand(seed1,corsika);
            alpha_res ares = alpha(r_glob, igrf_year, tn, E_gamma, chi_limit);
            double alfa_r = ares.alpha;
  
            if(alfa_r > alpha_max )
            {
                printf("Have to come back: on distance: %g alpha has value: %g \
                but actual alpha_max value is %g \n", 
                tr, alfa_r, alpha_max); // make some debug
                alpha_max = alfa_r; // assign the new alpha value
                tr = tr_prev; // go back and make step one more time ..
                continue; // in the next loop pass
            }
            
            double test_r = alfa_r/alpha_max;

            if(test_r <=rnd )
            {
                doLoop=1;
            }
            else
            {
                result.distance=tr;
                result.B_tr = ares.B_tr;
                result.alpha = alfa_r;
                doLoop=0;
            }
        }
        else
        {
            result.distance=-1.0;
            doLoop=0;
        }

    }// while
//    fclose(logfile);
    return result;
}

/**
 * tr_par2_r_curr function used to convert from tr_par trajectory distance
 * to r_curr distance to the earth center.
 * 
*/
float tr_par2_r_curr( float tn[4], float r_glob_0[4], float tr_par )
{
    float r_glob[4];
    /* parametric equation of the trajectory                              */
    r_glob[1]=r_glob_0[1]+tr_par*tn[1];
    r_glob[2]=r_glob_0[2]+tr_par*tn[2];
    r_glob[3]=r_glob_0[3]+tr_par*tn[3];
    /* current distance of the simulation point to the Earth's center     */

    float r_curr=normv(r_glob);
    return r_curr;
}


/** max_photon_alpha_rough function counts the value of alpha at the end of the
 * track (toa: the point with strongest B) and for a few further locations.
 * The maximum of all these alpha values is assumed to be the absolute maximum
 * of alpha along the track - max_alpha. This works well enough only for 
 * dipole-like magnetic field models
 * @param log_flag is an integer number which indicates if log file is created
 *        If it is set to 1 log file will be created with some values counted
 *        during the track   
  */

double max_photon_alpha_rough(double E_gamma,  float igrf_year, float toa,
                        float r_step, float tn[4], float sitec[4], 
			double chi_limit, int log_flag)
{
    float r_glob[4];
    double alpha_max = 0.0;
    double tr_par = 0;
    double r_curr =  (double)tr_par2_r_curr(tn, sitec, tr_par);
    double r_curr_0 = r_curr;
    double myr=0.0;
    long i=0;
    while (r_curr<=2*r_curr_0)
    {
        i++;
        /* parametric equation of the trajectory                              */
        r_glob[1]=sitec[1]-tr_par*tn[1];
        r_glob[2]=sitec[2]-tr_par*tn[2];
        r_glob[3]=sitec[3]-tr_par*tn[3];
        /* current distance of the simulation point to the Earth's center     */
        tr_par=r_step*km2cm*i;
        r_curr=normv(r_glob);
        alpha_res ares=  alpha(r_glob, igrf_year, tn, E_gamma,
                               chi_limit);
        double alpha_gamma = ares.alpha;
        if(alpha_gamma>alpha_max)
        {
            alpha_max=alpha_gamma;
            myr=r_curr;
        }
    }
    return alpha_max;
}

/** maxPhotonProb function counts maximal value of alpha function for distance
 * from the simulation start point to the end of the track (toa)
 * @param log_flag is an integer number which indicates if log file is created
 *        If it is set to 1 log file will be created with some values computed
 *        along the track   
 * @return max_alpha which is the max_value of alfa multiplied by 1.1 which ensueres
 * this value will never be less than alpha function.
  */


/** max_electron_alpha_rough function computes the maximum value of the 
 * bremsstralung probability function
 * at a distance from the simulation start point to the end of the track (toa)
 * @param log_flag is an integer number which indicates if log file is created
 *        If it is set to 1 log file will be created with some values computed
 *        during the track   
  */

double max_electron_alpha_rough(double E_e,  double E_sg0, float igrf_year, float toa,
                          float r_step, float tn[4], float r_glob_0[4],
                          float tr_par2, double chi_limit, int log_flag)
{
    double p_max=0.0;
    /////////////////////////////
    float r_glob[4];
    FILE *out;
    if(log_flag==1)
    {
        out = fopen("max_electron_prob_run.dat", "w");
    }
    double tr_par = 0;
    float myr=0.0;
    float mE_e=0.0;
    double r_curr =  (double)tr_par2_r_curr(tn, r_glob_0, tr_par);
    double r_curr_0 = r_curr;
    long i=0;
    while (E_e>E_sg0)
    {
      while (r_curr<=2*r_curr_0)
      {  
        i++;
        /* parametric equation of the trajectory                              */
        r_glob[1]=r_glob_0[1]-tr_par*tn[1];
        r_glob[2]=r_glob_0[2]-tr_par*tn[2];
        r_glob[3]=r_glob_0[3]-tr_par*tn[3];
        /* current distance of the simulation point to the Earth's center     */
        tr_par = r_step*km2cm*i;
        r_curr=normv(r_glob);
        //////////////////////////////////////////////
        //if(i%1000==0) printf(" %d | r: %g  | rc :%g\n", i, tr_par, r_curr);

        double prob_brem = prob_brem_total(E_sg0, E_e, r_glob, igrf_year, tn).prob_brem_total;


        //////////////////////////////////////////////////////
        if(log_flag==1)
        {
            fprintf(out, "%g %g\n", r_curr, prob_brem);
        }

//        printf("%e %g %g\n", E_e, r_curr, prob_brem);
        if(prob_brem>p_max)
        {
            p_max=prob_brem;
            myr=r_curr;
            mE_e=E_e;	    
        }
      }
      r_curr=r_curr_0;
      tr_par=0;
      i=0;
      E_e/=10.;
    }

    if(log_flag==1)
    {
        fclose(out);
    }
    printf("E_e %e, r %g, electron_alpha_max %g\n", mE_e, myr, p_max);

    ///////////////////////////////
    return p_max;
}



/** maxElectronProb function computes the maximum value of the bremsstralung probability 
 * function
 * for a distance from simulation start point to the end of track (toa)
 * @param log_flag is an integer number which indicates if log file is created
 *        If it is set to 1 log file will be created with some values computed
 *        along the track   
  */

/**
 * electron_veto function is used to find a distance at which electron bremsstrahlung
 * occurs. The function uses a veto algoritm (described e.g. here:
 * http://cepa.fnal.gov/psm/simulation/mcgen/lund/v6_214/pythia6206/node31.html ).
 * 
 * When algorithm finds a distance of bremsstrahlung it returns
 * also prob_brem_total value and B_tr value that were previously counted for 
 * this distance. They can be easily used in the main preshower simulation.
 * IMPORTANT: prob_brem_total is not the same value as used previously in the original
 * Preshower program. So to use it on a higher level it must be multiplayed by unit_dist
 * 
 * If the algorithm doesn't find the distance it returns the same structure
 * but the distance value is set to -1.0 which can be recognized from the upper 
 * level of program.
 * 
 * function can write some output to the file called logVeto_electron.dat, but to achieve
 * this, line with fprintf(...) must be uncommented and all recompilled.
 */
veto_el_res electron_veto(double E_e, double E_sg0, float igrf_year, float toa,
                          float tn[4], float r_glob_0[4],
                          float r_curr_0, float tr_par, double chi_limit,
                          long *seed1, int corsika, float electron_prob_max)
{
    float tr, rc,r_prev, rnd;
    float r_glob[4];
    int  doLoop, i;
    veto_el_res result;
//    FILE *logfile;
//    logfile=fopen("logVeto_electron.dat", "w");

    i=0;
    doLoop=1;
    tr = tr_par;

    // initialize rc
    r_glob[1]=r_glob_0[1] + tr*tn[1];
    r_glob[2]=r_glob_0[2] + tr*tn[2];
    r_glob[3]=r_glob_0[3] + tr*tn[3];
    rc=normv(r_glob);

    while(doLoop==1)
    {
        r_prev=tr;
        //    fprintf(logfile, "%d %g %g\n",i,r, r-r_prev );

        i++;
        rnd = getrand(seed1,corsika);
        tr = (electron_prob_max*r_prev - log(rnd) ) /electron_prob_max;

        // initializing r_glob table, which will be used to count rc
        // and to count prob_brem
        r_glob[1]=r_glob_0[1] + tr*tn[1];
        r_glob[2]=r_glob_0[2] + tr*tn[2];
        r_glob[3]=r_glob_0[3] + tr*tn[3];
        rc=normv(r_glob);

        // dummy way to prevent going on the other earth site
        if (r_curr_0<tr)
        {
            result.distance=-1.0;
            doLoop=0;
        }
        else if(rc>R_e+toa*km2cm)
        {
            rnd = getrand(seed1,corsika);
            brem_res bres =  prob_brem_total(E_sg0, E_e, r_glob, igrf_year, tn);
            double prob_brem_r = bres.prob_brem_total;
            
            if(prob_brem_r > electron_prob_max )
            {
                printf("Have to come back: on distance: %g alpha has value: %g \
                but actual alpha_max value is %g \n", 
                tr, prob_brem_r, electron_prob_max); // make some debug
                electron_prob_max = prob_brem_r; // assign the new alpha value
                tr = r_prev; // go back and make step one more time ..
                continue; // in the next loop pass
            }
           
            double test_r = prob_brem_r/electron_prob_max;

            if(test_r <=rnd )
            {
                doLoop=1;
            }
            else
            {
                result.distance=tr;
                result.B_tr = bres.B_tr;
                result.prob_brem_total = prob_brem_r;
                doLoop=0;
            }
        } // if atmosphere is reached
        else
        {
            result.distance=-1.0;
            doLoop=0;
        }

    }// while
//    fclose(logfile);
    return result;
}

