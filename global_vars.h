
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
#define  AMU          1.66053886e-24 /* 1 a.m.u */
#define  MPC          3.08568025e24
#define  KPC          3.08568025e21
#define  SIGMA_T      6.652458734e-25 
#define  SOLAR_MASS   1.989e33
#define  ELECTRONVOLT 1.602176565e-12

/* Atomic data (from VPFIT) */
#define  LAMBDA_LYA_H1  1215.6701e-8 /* cm */
#define  LAMBDA_LYA_HE2 303.7822e-8
#define  FOSC_LYA       0.416400
#define  GAMMA_LYA      6.265e8  /* s^-1 */
#define  HMASS          1.00794  /* Hydrogen mass in a.m.u. */
#define  HEMASS         4.002602 /* Helium-4 mass in a.m.u. */

/* Gaussian profile table */
#define GXMAX  100.0
#define NGXTAB 1000000

