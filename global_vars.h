
/********************************************************************************/

/**! \file global_vars.h                                                          
 *                                                
 *  \brief defines common variables and constants.
 */

/********************************************************************************/

double ztime,omegam,omegal,omegab,h100,box100,Xh;

#ifdef SELF_SHIELD
  double n0_z;
#endif

#if defined(TEST_KERNEL) && !defined(NO_PECVEL)
#error "TEST_KERNEL requires NO_PECVEL"
#endif


/* Numbers */
#define  PI    3.14159265358979323846
#define  GAMMA (5.0/3.0)

/* Physical constants (cgs units) */
/* See http://physics.nist.gov/cuu/Constants/index.html */ 
#define  GRAVITY      6.67384e-8
#define  BOLTZMANN    1.3806488e-16
#define  C            2.99792458e10
#define  AMU          1.66053886e-24 /* 1 a.m.u or Dalton */
#define  MPC          3.08568025e24
#define  KPC          3.08568025e21
#define  SIGMA_T      6.652458734e-25 
#define  SOLAR_MASS   1.989e33
#define  ELECTRONVOLT 1.602176565e-12

/* Atomic data. See Table 2, Morton 2003, ApJS, 149, 205 and
   https://ciaaw.org/index.htm */

/* Hydrogen */
#define  HMASS          1.00794  /* Hydrogen atomic weight in a.m.u. */
#define  LAMBDA_LYA_H1  1215.6701e-8 /* cm */
#define  FOSC_LYA_H1    0.416400
#define  GAMMA_LYA_H1   6.265e8  /* s^-1 */

/* Helium */
#define  HEMASS         4.002602 /* Helium-4 mass in a.m.u. */
#define  LAMBDA_LYA_HE2 303.7822e-8
#define  FOSC_LYA_HE2   0.416400
#define  GAMMA_LYA_HE2  6.265e8  /* s^-1 */

/* Silicon */
#define  SIMASS          28.085  /* Silicon atomic weight in a.m.u */
#define  SI_SOLAR        3.236e-5  /* Asplund et al (2009), Table 1,  10^(7.60-12) */
#define  LAMBDA_1190_Si2 1190.4158e-8
#define  FOSC_1190_Si2   0.2920
#define  LAMBDA_1193_Si2 1193.2897e-8
#define  FOSC_1193_Si2   0.5820
#define  LAMBDA_1260_Si2 1260.4221e-8
#define  FOSC_1260_Si2   1.180
#define  LAMBDA_1207_Si3 1206.50e-8
#define  FOSC_1207_Si3   1.63000
 

/* Gaussian profile table */
#define GXMAX  100.0
#define NGXTAB 1000000

