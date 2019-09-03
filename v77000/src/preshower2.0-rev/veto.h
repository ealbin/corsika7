#ifndef VETO_H_
#define VETO_H_
/**
 * @file veto.h
 * @author agata
 * @brief package for functions which are used by veto algorithm 
 * 
 */
  
veto_res photon_veto(double E_gamma,  float igrf_year,float toa,  float alpha_max,
                     float tn[4], float r_glob_0[4],
                     float r_curr_0, float tr_par, double chi_limit,
                     long *seed1, int corsika)  ;


///////////////////////////////////////////


veto_el_res electron_veto(double E_e, double E_sg0, float igrf_year, float toa,
                          float tn[4], float r_glob_0[4],
                          float r_curr, float tr_par, double chi_limit,
                          long *seed1, int corsika, float electron_prob_max);

///////////////////////////////////////////

float tr_par2_r_curr( float tn[4], float r_glob_0[4], float tr_par );


double max_photon_alpha_rough(double E_gamma,  float igrf_year, float toa,
                           float r_step, float tn[4], float r_glob_0[4],
                           double chi_limit, int log_flag);

///////////////////////////////


double max_electron_alpha_rough(double E_e,  double E_sg0, float igrf_year, float toa,
                          float r_step, float tn[4], float r_glob_0[4],
                          float tr_par2, double chi_limit, int log_flag);

#endif /*VETO_H_*/
