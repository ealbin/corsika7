/* $Id: CorsikaConsts.h 5116 2016-01-04 19:09:04Z darko $   */

#ifndef _INCLUDE_CORSIKA_CONSTS_H_
#define _INCLUDE_CORSIKA_CONSTS_H_

#include <string>
#include <cmath>


namespace crs {

  static const double eSI   = 1.602176462e-19; // positron charge in coulomb


  static const double ns = 1.0;          // CORSIKA unit
  static const double nanosecond = ns;
  static const double microsecond = 1.e3*ns;
  static const double millisecond = 1.e6*ns;
  static const double second = 1.e9*ns;
  static const double s = 1.e9*ns;
  
  static const double cm = 1.0;          // CORSIKA unit
  static const double mm = 1.e-1*cm; 
  static const double meter = 1.e2*cm;
  static const double m = meter;
  static const double km = 1.e3*m;

  static const double m2 = m*m;

  static const double GeV = 1;           // CORSIKA unit
  static const double MeV = 1.e-3*GeV;
  static const double keV = 1.e-6*GeV;
  static const double eV = 1.e-9*GeV;
  static const double TeV = 1.e3*GeV;
  static const double PeV = 1.e6*GeV;
  static const double EeV = 1.e9*GeV;
  static const double ZeV = 1.e12*GeV;
  static const double joule = eV/eSI; // joule = 6.24150 e+12 * MeV

  static const double mol = 1;

  static const double rad = 1.;          // CORSIKA unit
  static const double deg = M_PI/180.;

  static const double kelvin = 1.;

  static const double newton = joule/m; // newton = 6.24150 e+9 * MeV/mm

  static const double  kilogram = joule*second*second/(m*m);   
  static const double      gram = 1.e-3*kilogram;
  static const double milligram = 1.e-3*gram;

  static const double g = gram;     
  static const double kg = kilogram;

  static const double pascal     = newton/m2; // pascal=6.24150 e+3 * MeV/mm3
  static const double bar        = 100000*pascal;// bar=6.24150 e+8 * MeV/mm3
  static const double atmosphere = 101325*pascal;// atm=6.32420 e+8 * MeV/mm3
  static const double millibar   = 1e-3 * bar;
  
  
  
  static const float cSpeedOfLight = 29.9792458*cm/ns;	// c in [cm/ns]
  static const double cRearth = 6371315.00*m;      // m REarth = 637131500cm
  static const double cTopOfAtmosphere = 112.8292 * km;



  static const int gNParticles = 200;

  // particle masses (REVISED SEPT 2000 BY D. HECK)
  // MASSES REVISED          NOV 2004 BY D. HECK
  // CHARMED PARTICLES ADDED     2008 BY D. HECK
  // BOTTOM PARTICLES ADDED JULY 2012 BY A. GASCON
  static const float gParticleMass [gNParticles] = {
     0.e0	,.51099892e-3,.51099892e-3,  0.e0     ,.105658369e0,
     .105658369e0, .1349766e0, .13957018e0,.13957018e0, 0.497648e0 ,//10
     0.493677e0 , 0.493677e0 ,.93956536e0 ,.93827203e0,.93827203e0 ,
     0.497648e0 , 0.54775e0  , 1.115683e0 , 1.18937e0 , 1.192642e0 ,//20
     1.197449e0 , 1.31483e0  , 1.32131e0  , 1.67245e0 ,.93956536e0 ,
     1.115683e0 , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31483e0  ,//30
     1.32131e0  , 1.67245e0  , 0.e0	  , 0.e0      , 0.e0	   ,
     0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//40
     1.e9	, 580.e0     , 1.e5	  , 0.e0      , 0.e0	   ,
     0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.78259e0  ,//50
     0.7690e0	, 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
     1.2331e0	, 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   ,//60
     1.2349e0	, 0.89610e0  , 0.89166e0  , 0.89166e0 , 0.89610e0  ,
     0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//70
     0.54775e0  , 0.54775e0  , 0.54775e0  , 0.54775e0 ,.105658369e0,
     .105658369e0 , 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//80
     0.e0	, 0.e0       , 0.e0	  , 0.e0      ,.105658369e0,
     .105658369e0, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//90
     0.e0	, 0.e0       , 0.e0	  , 0.e0      ,.105658369e0,
     .105658369e0, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//100
     0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,
     0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//110
     0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,
     1.8645e0	, 1.8697e0   , 1.8697e0   , 1.8645e0   , 1.9682e0   ,//120
     1.9682e0	, 2.9804e0   , 2.0067e0   , 2.0100e0   , 2.0100e0   ,
     2.0067e0	, 2.1121e0   , 2.1121e0   , 0.0e0      , 3.096916e0 ,//130
     1.77699e0  , 1.77699e0  , 0.e0	  , 0.e0       , 0.e0	    ,
     0.e0	, 2.28646e0  , 2.4679e0   , 2.4710e0   , 2.45402e0  ,//140
     2.4529e0	, 2.45376e0  , 2.5757e0   , 2.5780e0   , 2.6975e0   ,
     0.e0	, 0.e0       , 0.e0	  , 2.28646e0  , 2.4679e0   ,//150
     2.4710e0	, 2.45402e0  , 2.4529e0   , 2.45376e0  , 2.5757e0   ,
     2.5780e0	, 2.6975e0   , 0.e0	  , 0.e0       , 0.e0	    ,//160
     2.5184e0	, 2.5175e0   , 2.5180e0   , 0.e0       , 0.e0	    ,
     0.e0	, 0.e0       , 0.e0	  , 0.e0       , 0.e0	    ,//170
     2.5184e0	, 2.5175e0   , 2.5180e0   , 0.e0       , 0.e0	    ,
     5.27958e0  , 5.27925e0  , 5.27925e0  , 5.27958e0  , 5.36677e0  ,//180
     5.36677e0  , 6.277e0    , 6.277e0    , 5.6194e0   , 5.8155e0   ,
     5.8113e0	, 5.788e0    , 5.7911e0   , 6.071e0    , 5.6194e0   ,//190
     5.8155e0	, 5.8113e0   , 5.788e0    , 5.7911e0   , 6.071e0    ,
     0.e0	, 0.e0       , 0.e0	  , 0.e0       , 0.e0	    }; //200


  // particle charges
  static const float gParticleCharge [gNParticles] = {
     0.0,+1.0,-1.0, 0.0,+1.0,-1.0, 0.0,+1.0,-1.0, 0.0, //10
    +1.0,-1.0, 0.0,+1.0,-1.0, 0.0, 0.0, 0.0,+1.0, 0.0, //20
    -1.0, 0.0,-1.0,-1.0, 0.0, 0.0,-1.0, 0.0,+1.0, 0.0, //30
    +1.0,+1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //40
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //50
     0.0,+1.0,-1.0,+2.0,+1.0, 0.0,-1.0,-2.0,-1.0, 0.0, //60
    +1.0, 0.0,+1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //70
     0.0, 0.0, 0.0, 0.0,+1.0,-1.0, 0.0, 0.0, 0.0, 0.0, //80
     0.0, 0.0, 0.0, 0.0,+1.0,-1.0, 0.0, 0.0, 0.0, 0.0, //90
     0.0, 0.0, 0.0, 0.0,+1.0,-1.0, 0.0, 0.0, 0.0, 0.0, //100
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //110
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,+1.0,-1.0, 0.0,+1.0, //120
    -1.0, 0.0, 0.0,+1.0,-1.0, 0.0,+1.0,-1.0, 0.0, 0.0, //130
    -1.0,+1.0, 0.0, 0.0, 0.0, 0.0,+1.0,+1.0, 0.0,+2.0, //140
    +1.0, 0.0,+1.0, 0.0, 0.0, 0.0, 0.0, 0.0,-1.0,-1.0, //150
     0.0,-2.0,-1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, //160
    +2.0,+1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,+2.0, //170
    -2.0,-1.0, 0.0, 0.0, 0.0, 0.0,+1.0,-1.0, 0.0, 0.0, //180
     0.0,+1.0,-1.0, 0.0,-1.0,+1.0, 0.0,-1.0,-1.0, 0.0, //190
    +1.0,-1.0, 0.0,+1.0,+1.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //200


  // particle names (might be incomplete)
  static const std::string gParticleName [gNParticles] = {
    "gamma", "e+", "e-", "", "mu+", 
    "mu-", "pi0", "pi+", "pi-", "K0_L", // 10
    "K+", "K-", "n", "p", "p_bar", 
    "K0_S", "eta", "Lambda", "Sigma+", "Sigma0", // 20
    "Sigma-", "Xi0", "Xi-", "Omega-", "n_bar", 
    "Lambda_bar", "Sigma-_bar", "Sigma0_bar", "Sigma+_bar", "Xi0_bar",// 30
    "Xi+_bar", "Omega+_bar", "", "", "",
    "", "", "", "", "", // 40
    "", "", "", "", "", 
    "", "", "", "", "omega", // 50
    "rho0", "rho+", "rho-", "Delta++", "Delta+", 
    "Delta0", "Delta-", "Delta--_bar", "Delta-_bar", "Delta0_bar", // 60
    "Delta+_bar", "K*0", "K*+", "K*-", "K*0_bar", 
    "nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "", // 70
    "eta->2gamma", "eta->3pi0", "eta->pi+pi-pi0", "eta->pi+pi-gamma", 
    "mu+ prod. info", 
    "mu- prod. info", "", "", "", "", // 80
    "", "", "", "", "dky mu+ start info","dky mu- start info",
    "", "", "", "",  // 90
    "", "", "", "", "dky mu+ end info","dky mu- end info",
    "", "", "", "",  // 100
    "", "", "", "", "", 
    "", "", "", "", "", // 110
    "", "", "", "", "", 
    "D0", "D+", "D-_bar", "D0_bar", "D_a+", // 120
    "D_a-_bar", "eta_c", "D*0", "D*+", "D*-_bar", 
    "D*0_bar", "D*_a+", "D*_a-_bar", "J/Psi", "", // 130
    "tau+", "tau-", "nu_tau", "nu_tau_bar", "",
    "", "Lambda+_c", "Xi+_c", "Xi0_c", "Sigma++_c", // 140
    "Sigma+_c", "Sigma0_c", "Xi`+_c", "Xi`0_c", "Omega0_c", 
    "", "", "", "Lambda-_c", "Xi-_c_bar", // 150
    "Xi0_c_bar","Sigma--_c_bar","Sigma-_c_bar","Sigma0_c_bar","Xi`-_c_bar",
    "Xi`0_c_bar", "Omega0_c_bar", "", "", "", // 160
    "Sigma*++_c", "Sigma*+_c", "Sigma*0_c", "", "", 
    "", "", "", "", "", // 170
    "Sigma*--_c_bar", "Sigma*-_c_bar", "Sigma*0_c_bar", "", "", 
    "", "", "", "", "", // 180
    "", "", "", "", "", 
    "", "", "", "", "", // 190
    "", "", "", "", "", 
    "", "", "", "", "" // 200
  };
    
    
    
  // particle codes
  static const int gParticleCodes [gNParticles] = {
    22, -11, 11, 0, -13, 
    13, 111, 211, -211, 130,              // 10
    321, -321, 2112, 2212, -2212,
    310, 221, 3122, 3222, 3212,        // 20
    3112, 3322, 3312, 3334, -2112, 
    -3122, -3112, -3212, -3222, -3322, // 30
    -3312, -3334, 0, 0, 0, 
    0, 0, 0, 0, 0,                     // 40
    0, 0, 0, 0, 0,                  
    0, 0, 0, 0, 223,                     // 50
    113, 213, -213, 2224, 2214,
    2114, 1114, -2224, -2214, -2114,     // 60
    -2214, 313, 323, -323, -313,
    12, -12, 14, -14, 0,                 // 70
    0, 0, 0, 0, -13,
    13, 0, 0, 0, 0,                        // 80
    0, 0, 0, 0, -13,
    13, 0, 0, 0, 0,                        // 90
    0, 0, 0, 0, -13,
    13, 0, 0, 0, 0,                        // 100
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,                        // 110
    0, 0, 0, 0, 0,
    421, 411, -411, -421, 431,           // 120
    -431, 441, 423, 413, -413,
    -423, 433, -433, 0, 443,              // 130
    -15, 15, 16, -16, 0,
    0, 4122, 4232, 4132, 4222,           // 140
    4212, 4112, 4322, 4312, 4332,
    0,0,0, -4122, -4323,                 // 150
    -4132, -4222, -4212, -4112, -4322,
    -4312, -4332, 0, 0, 0,                // 160
    4224, 4214, 4114, 0, 0,
    0, 0, 0, 0, 0,                        // 170
    -4224, -4214, -4114, 0, 0,
    0, 0, 0, 0, 0,                        // 180
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,                        // 190
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0                        // 200
  };
    
    
};

#endif
