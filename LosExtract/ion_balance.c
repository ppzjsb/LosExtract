
/********************************************************************************/

/*! \file ion_balance.c                                                                                
 *                                                                                                         
 *  \brief Routine for recomputing the HI fraction in pixels that are self-shielded
 */

/********************************************************************************/

#ifdef SELF_SHIELD

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global_vars.h"
#include "parameters.h"
#include "proto.h"

double gJH0, gJHe0, gJHep, gJH0_SS, J_UV;
double aHp, aHep, aHepp, ad, geH0, geHe0, geHep;

double alpha1_z, alpha2_z, beta_z, fval_z;

static double deltaT;
static double *AlphaHp, *AlphaHep, *AlphaHepp, *Alphad;
static double *GammaeH0, *GammaeHe0, *GammaeHep;
static double logTmin, logTmax;


#define SMALLNUM 1.0e-60
#define TMIN 1.0
#define TMAX 1.0e9
#define MAXITER  150
#define NCOOLTAB 2000

/********************************************************************************/

/* \brief Routine to (re-)compute the neutral hydrogen fraction using
 *  the self-shielding correction of Chardin et al. 2018, MNRAS, 478,
 *  1065 and Rahmati et al. 2013, MNRAS, 430, 2427 
 *
 * \param nHcgs The hydrogen number density in cgs units [cm^-3]
 * \param logT  log10 of the gas temperature [K]
 *
 * \return nH0 The neutral hydrogen fraction, nHI/nH [dimensionless]
 */

/********************************************************************************/

double self_shield(double nHcgs, double logT)
{
  double fnH_SS = nHcgs/n0_z; 
  
  /* Eq. (1), Chardin et al. 2018 and Eq. (A1) Rahmati et al. 2013  */
  gJH0_SS = gJH0 * ((1.0-fval_z) * pow(1.0 + pow(fnH_SS, beta_z), alpha1_z) 
		    + fval_z * pow(1.0 + fnH_SS, alpha2_z));
  
  double nH0 = ion_balance(nHcgs, logT);
  
  return nH0; /* HI/H */
}
/********************************************************************************/

/* \brief This function computes the equilibrium abundance ratios 
 *
 * \param nHcgs The hydrogen number density in cgs units [cm^-3]
 * \param logT  log10 of the gas temperature [K]
 *
 * \return nH0 The neutral hydrogen fraction, nHI/nH [dimensionless]
 */

/********************************************************************************/

double ion_balance(double nHcgs, double logT)
{
  
  double gJH0ne  = 0.0;
  double gJHe0ne = 0.0;
  double gJHepne = 0.0;
  double yhelium = (1.0 - Xh) / (4.0 * Xh);

  /* Initialised values only */
  double nHp   = 0.0;
  double nHep  = 0.0;
  double nHepp = 0.0;
  double nH0   = 0.0;
  double nHe0  = 0.0;
  double ne    = 1.0;
  
  if(logT <= logTmin) /* everything neutral */
    {
      nH0 = 1.0;
      return nH0;
    }
  
  if(logT >= logTmax) /* everything is ionised */
    {
      nH0   = 0.0;
      return nH0;
    }
  
  double t    = (logT - logTmin) / deltaT;
  int j       = (int) t;
  double fhi  = t - j;
  double flow = 1 - fhi;
  
  double necgs = ne * nHcgs;
  double neold = ne;
  
  int niter = 0;
 
  /* Evaluate number densities iteratively in units of nH */
  do
    {
      niter++;
          
      /* Recombination rates*/
      aHp   = flow * AlphaHp[j]   + fhi * AlphaHp[j + 1];
      aHep  = flow * AlphaHep[j]  + fhi * AlphaHep[j + 1];
      aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
      ad    = flow * Alphad[j]    + fhi * Alphad[j + 1];

      /* Collisional ionisation rates */
      geH0  = flow * GammaeH0[j]  + fhi * GammaeH0[j + 1];
      geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
      geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];

    
      if(necgs <= 1.e-25 || J_UV == 0 )
	{
	  gJH0ne  = 0.0;
	  gJHe0ne = 0.0;
	  gJHepne = 0.0;
	}
      else
	{
	  gJH0ne  = gJH0_SS / necgs; /* NOTE: uses RT correction */
	  gJHe0ne = gJHe0   / necgs;
	  gJHepne = gJHep   / necgs;
	}
      
      /* Follows Katz, Weinberg & Hernquist, 1996, ApJS, 105, 19, eq. 33-38 */
      nH0 = aHp / (aHp + geH0 + gJH0ne);	
      nHp = 1.0 - nH0;		
      
      if((gJHe0ne + geHe0) <= SMALLNUM)	/* no ionisation at all */
	{
	  nHep = 0.0;
	  nHepp = 0.0;
	  nHe0 = yhelium;
	}
      else
	{
	  nHep  = yhelium / (1.0 + (aHep+ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	
	  nHe0  = nHep * (aHep+ad) / (geHe0 + gJHe0ne);	
	  nHepp = nHep * (geHep + gJHepne) / aHepp;	
	}
      
      neold = ne;
      ne    = nHp + nHep + 2.0 * nHepp;	
      necgs = ne * nHcgs;
      
      
      if(J_UV == 0)
	break; 
      
      double nenew = 0.5 * (ne + neold);
      
      ne    = nenew;
      necgs = ne * nHcgs;
      
      
      if(fabs(ne - neold) < 1.0e-4)
	break;
      
      if(niter > (MAXITER - 10))
	printf("ne= %g  niter=%d\n", ne, niter);
    }
  while(niter < MAXITER);
  
  
  if(niter >= MAXITER)
    {
      printf("no convergence reached in ion_balance()\n");
      exit(0);
    }

  return nH0;
  
}

/********************************************************************************/

/* \brief Initialise cooling and ionisation table at start of run
 */

/********************************************************************************/

void InitCool(void)
{
  char fname[400];
  
  sprintf(fname,"./treecool/%s",UVBFILE);

  InitCoolMemory();
  
  MakeCoolingTable();
  
  ReadIonizeParams(fname);

  IonizeParamsTable();
  
  SelfShieldFit();
}

/******************************************************************************/

/* \brief Make the look-up table for the collisional ionisation rates
 * and recombination coefficients
 */

/******************************************************************************/

void MakeCoolingTable(void)  
{
  
  /* Verner & Ferland case-A recombination rate fit coefficients */
  double a_vf96[3] = {7.982e-11, 9.356e-10, 1.891e-10};
  double b_vf96[3] = {0.7480   , 0.7892   , 0.7524};
  double T0_vf96[3]= {3.148e0  , 4.266e-2 , 9.370e0};
  double T1_vf96[3]= {7.036e5  , 4.677e6  , 2.7674e6};
  
  /* Voronov collisional ionisation rate fit coefficients */
  double dE_vor97[3] = {13.6    , 24.6    , 54.4};
  double P_vor97[3]  = {0       , 0       , 1};
  double A_vor97[3]  = {0.291e-7, 0.175e-7, 0.205e-8};
  double X_vor97[3]  = {0.232   , 0.180   , 0.265};
  double K_vor97[3]  = {0.39    , 0.35    , 0.25};
  
  
  logTmin = log10(TMIN);
  logTmax = log10(TMAX);
  deltaT  = (logTmax - logTmin) / NCOOLTAB;
 
  for(int i = 0; i <= NCOOLTAB; i++)
    {
      double T      = pow(10.0, logTmin + deltaT * i);
      double T_eV   = T * BOLTZMANN/ELECTRONVOLT;
  
      
      /* Case-A recombination rates, [cm^3 s^-1]
	 Verner & Ferland, 1996, ApJS, 103, 467 */
      AlphaHp[i]  = a_vf96[0]/(sqrt(T/T0_vf96[0])*pow(1.0+sqrt(T/T0_vf96[0]),1.0-b_vf96[0])
			       *pow(1.0+sqrt(T/T1_vf96[0]),1.0+b_vf96[0]));
      
      AlphaHep[i] = a_vf96[1]/(sqrt(T/T0_vf96[1])*pow(1.0+sqrt(T/T0_vf96[1]),1.0-b_vf96[1])
			       *pow(1.0+sqrt(T/T1_vf96[1]),1.0+b_vf96[1]));
	
      AlphaHepp[i]= a_vf96[2]/(sqrt(T/T0_vf96[2])*pow(1.0+sqrt(T/T0_vf96[2]),1.0-b_vf96[2])
			       *pow(1.0+sqrt(T/T1_vf96[2]),1.0+b_vf96[2]));
      
      
      /* He+ dielectronic recombination rate [cm^3 s^-1]
	 Aldrovandi & Pequignot, 1973, A&A, 25, 137 */
      if(4.7e5/T < 70) 
        Alphad[i] = 1.9e-3*(1.0 + 0.3*exp(-9.4e4/T))*exp(-4.7e5/T)*pow(T,-1.5);

      
      /* Collisional ionisation rates [cm^3 s^-1]
	 Voronov 1997, ADNDT, 65, 1 */
      if(dE_vor97[0]/T_eV < 70)
	GammaeH0[i] = A_vor97[0]*pow(dE_vor97[0]/T_eV, K_vor97[0])*exp(-dE_vor97[0]/T_eV) 
	  *(1.0+P_vor97[0]*sqrt(dE_vor97[0]/T_eV)) / (X_vor97[0] + dE_vor97[0]/T_eV);
      
      if(dE_vor97[1]/T_eV < 70)
        GammaeHe0[i] = A_vor97[1]*pow(dE_vor97[1]/T_eV, K_vor97[1])*exp(-dE_vor97[1]/T_eV) 
	  *(1.0+P_vor97[1]*sqrt(dE_vor97[1]/T_eV)) / (X_vor97[1] + dE_vor97[1]/T_eV);
      
      if(dE_vor97[2]/T_eV < 70)
        GammaeHep[i] = A_vor97[2]*pow(dE_vor97[2]/T_eV, K_vor97[2])*exp(-dE_vor97[2]/T_eV) 
	  *(1.0+P_vor97[2]*sqrt(dE_vor97[2]/T_eV)) / (X_vor97[2] + dE_vor97[2]/T_eV);
      
    }
  
}

/******************************************************************************/

/* \brief Reads in the information from TREECOOL at start of run 
 */

/******************************************************************************/

#define TABLESIZE 500

static float inlogz[TABLESIZE];
static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
static int nheattab;		/* length of table */


void ReadIonizeParams(char *fname)
{
  FILE *fdcool;
  
  if(!(fdcool = fopen(fname, "r")))
    {
      printf(" Cannot read ionization table in file `%s'\n", fname);
      exit(0);
    }

  for(int i = 0; i < TABLESIZE; i++)
    gH0[i] = 0;
  
  for(int i = 0; i < TABLESIZE; i++)
    if(fscanf(fdcool, "%g %g %g %g %g %g %g",
	      &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF)
      break;

  fclose(fdcool);

  /*  nheattab is the number of entries in the table */
  for(int i = 0, nheattab = 0; i < TABLESIZE; i++)
    if(gH0[i] != 0.0)
      nheattab++;
    else
      break;
  
}

/******************************************************************************/

/* \brief  Computes the photo-ionisation and photo-heating rates from the
 *  external TREECOOL file
 */

/******************************************************************************/

void IonizeParamsTable()
{
  
  float logz = log10(ztime + 1.0);
  int ilow   = 0;
  
  for(int i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
	ilow = i;
      else
	break;
    }
  
  float dzlow = logz - inlogz[ilow];
  float dzhi = inlogz[ilow + 1] - logz;
  
  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      gJHe0 = gJHep = gJH0 = 0.0;
      J_UV  = 0.0;
      return;
    }
  else
    J_UV = 1.e-21;		/* irrelevant as long as it's not 0 */
  
  
  gJH0   = pow(10., (dzhi * log10(gH0[ilow]) + dzlow  * log10(gH0[ilow + 1]))  / (dzlow + dzhi));
  gJHe0  = pow(10., (dzhi * log10(gHe[ilow]) + dzlow  * log10(gHe[ilow + 1]))  / (dzlow + dzhi));
  gJHep  = pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ilow + 1])) / (dzlow + dzhi));
  
  return;
}
/******************************************************************************/

/* \brief Performs a linear interpolation to obtain the parameters for
 *  the self-shielding correction
 *
 *  The tabulated fits for the self-shielding correction to the HI
 *  photo-ionisation rate use Rahmati et al. (2013) from z=0 to z=2,
 *  and Chardin et al. (2017) from z=3 to z=10.  The redshift bins are
 *  separated by dz=1.0.  See Table A1 and Eq. A1 in Rahmati et
 *  al. (2013), and Table A1 in Chardin et al. (2017)
 */

/******************************************************************************/

void SelfShieldFit()
{
  
  double n0[11]     = {1.148e-3, 5.129e-3, 8.710e-3, 9.0e-3, 9.3e-3, 1.03e-2,
		       7.0e-3, 2.7e-3, 4.0e-3, 4.6e-3, 4.7e-3};
  
  double alpha1[11] = {-3.98, -2.94, -2.22, -1.12, -0.95, -1.29, 
		       -0.94, -0.86, -0.74, -0.64, -0.39};
  
  double alpha2[11] = {-1.09, -0.90, -1.09, -1.65, -1.50, -1.60, 
		       -1.51, -1.27, -1.40, -1.21, -0.86};
  
  double beta[11]   = {1.29, 1.21, 1.75, 5.32, 5.87, 5.06,
		       6.11, 7.08, 7.12, 9.99, 12.94};
  
  double fval[11]   = {0.01, 0.03, 0.03, 0.018, 0.015, 0.024, 
		       0.029, 0.041, 0.041, 0.029, 0.006};
    
  double zmin   = 0.0;
  double deltaz = 1.0; 
  
  double t    = (ztime - zmin) / deltaz;
  int j       = (int)t;
  double fhi  = t - j;
  double flow = 1.0 - fhi;
  
  if(j < 10)
    {
      n0_z     = flow * n0[j]     + fhi * n0[j+1];
      alpha1_z = flow * alpha1[j] + fhi * alpha1[j+1];
      alpha2_z = flow * alpha2[j] + fhi * alpha2[j+1];
      beta_z   = flow * beta[j]   + fhi * beta[j+1];
      fval_z   = flow * fval[j]   + fhi * fval[j+1];
    }
  else
    {
      n0_z     = n0[10];
      alpha1_z = alpha1[10];
      alpha2_z = alpha2[10];
      beta_z   = beta[10];
      fval_z   = fval[10];
    }
    
  return;
}

/******************************************************************************/

/* \brief  Allocate the arrays for the cooling function table 
 */

/******************************************************************************/

void InitCoolMemory(void)
{
  AlphaHp = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==AlphaHp)
    {
      free(AlphaHp);
      printf("Memory allocation failed for AlphaHp.\n");
      exit(0);
    }

  AlphaHep = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==AlphaHep)
    {
      free(AlphaHep);
      printf("Memory allocation failed for AlphaHep.\n");
      exit(0);
    }
  
  AlphaHepp = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==AlphaHepp)
    {
      free(AlphaHepp);
      printf("Memory allocation failed for AlphaHepp.\n");
      exit(0);
    }
  
  Alphad = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==Alphad)
    {
      free(Alphad);
      printf("Memory allocation failed for Alphad.\n");
      exit(0);
    }

  GammaeH0 = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==GammaeH0)
    {
      free(GammaeH0);
      printf("Memory allocation failed for GammaeH0.\n");
      exit(0);
    }

  GammaeHe0 = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==GammaeHe0)
    {
      free(GammaeHe0);
      printf("Memory allocation failed for GammaeHe0.\n");
      exit(0);
    }

  GammaeHep = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==GammaeHep)
    {
      free(GammaeHep);
      printf("Memory allocation failed for GammaeHep.\n");
      exit(0);
    }
}

#endif

/******************************************************************************/
