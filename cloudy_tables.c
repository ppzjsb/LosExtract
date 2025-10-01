
/********************************************************************************/

/**! \file cloudy_tables.c                                                                                
 *                                                                                                         
 *  \brief Routine for reading in pre-computed cloudy tables with
 *  metal ion abundances.  These can then be interpolated with respect
 *  to log10(nH/cm^-3) and log10(T/K) 
 */

/********************************************************************************/

#ifdef SILICON

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global_vars.h"
#include "parameters.h"
#include "proto.h"

double *Si2frac, *Si3frac, *Si4frac;
double *Si2frac_tab, *Si3frac_tab, *Si4frac_tab;
double ztab_min,ztab_max,ttab_min,ttab_max,dtab_min,dtab_max;

int n_ztab,n_ttab,n_dtab;


/** \brief Gets the ion fractions from the cloudy external table for a
 *   given redshift, temperature log10(T/K) and hydrogen number density
 *   log10(nH/cm^-3).  Linear interpolation is used.
 *
 *  \param logT   The input gas temperature, log10(T/K)
 *  \param lognH  The input hydrogen number density, log10(nH/cm^-3)
 *  \param *Si2frac   Pointer containing the requested Si+/Si (i.e. Si-II)
 *  \param *Si3frac   Pointer containing the requested Si2+/Si (i.e. Si-III)
 *  \param *Si4frac   Pointer containing the requested Si3+/Si (i.e. Si-IV)
 */

void get_ionfrac(double logT, double lognH, double *Si2frac, double *Si3frac, double *Si4frac)
{ 
  double *Si2frac_zint = (double *)calloc(n_ttab*n_dtab, sizeof(double));
  if(NULL==Si2frac_zint){free(Si2frac_zint); printf("Memory allocation failed.\n"); exit(0);}
  
  double *Si3frac_zint = (double *)calloc(n_ttab*n_dtab, sizeof(double));
  if(NULL==Si3frac_zint){free(Si3frac_zint); printf("Memory allocation failed.\n"); exit(0);}
  
  double *Si4frac_zint = (double *)calloc(n_ttab*n_dtab, sizeof(double));
  if(NULL==Si4frac_zint){free(Si4frac_zint); printf("Memory allocation failed.\n"); exit(0);}
  
  
  /* Redshift interpolation */
  double ztab_bin_inv = ((double)n_ztab - 1.0) / (ztab_max - ztab_min);  
  
  double t    = (ztime - ztab_min) * ztab_bin_inv;
  int j       = floor(t);
  double fhi  = t - j;
  double flow = 1 - fhi;
  
  if(j > n_ztab - 2 || j < 0) /* protect boundaries */
    {
      printf("z=%.3f is outside range of cloudy tables! Stopping...\n",ztime);
      printf("z_min=%.3f z_max<%.3f\n\n",ztab_min,ztab_max);
      exit(0);
    }
  
  int it, id;
  long long zoff = n_ttab * n_dtab; /* offset between redshift grid points in 1D array */
  for(id=0; id<n_dtab; id++)
    for(it=0; it<n_ttab; it++)
      {
	long long ind = it + id*n_ttab;
	Si2frac_zint[ind] = flow * Si2frac_tab[ind + j*zoff] + fhi * Si2frac_tab[ind + (j+1)*zoff];
	Si3frac_zint[ind] = flow * Si3frac_tab[ind + j*zoff] + fhi * Si3frac_tab[ind + (j+1)*zoff];
	Si4frac_zint[ind] = flow * Si4frac_tab[ind + j*zoff] + fhi * Si4frac_tab[ind + (j+1)*zoff];
      }
  
  
 
  double *Si2frac_dint = (double *)calloc(n_ttab, sizeof(double));
  if(NULL==Si2frac_dint){free(Si2frac_dint); printf("Memory allocation failed.\n"); exit(0);}

  double *Si3frac_dint = (double *)calloc(n_ttab, sizeof(double));
  if(NULL==Si3frac_dint){free(Si3frac_dint); printf("Memory allocation failed.\n"); exit(0);}

  double *Si4frac_dint = (double *)calloc(n_ttab, sizeof(double));
  if(NULL==Si4frac_dint){free(Si4frac_dint); printf("Memory allocation failed.\n"); exit(0);}
  
  
  /* log(nH/cm^-3) interpolation */
  lognH = dmin(dmax(lognH, dtab_min), dtab_max-1.0e-6); /* protect boundaries */
  double dtab_bin_inv = ((double)n_dtab - 1.0) / (dtab_max - dtab_min);
  
  t    = (lognH - dtab_min ) * dtab_bin_inv;
  j    = floor(t);
  fhi  = t - j;
  flow = 1 - fhi;
  
  long long doff = n_ttab; /* offset between density grid points in 1D array */
  for(it=0; it<n_ttab; it++)
    {
      Si2frac_dint[it] = flow * Si2frac_zint[it + j*doff] + fhi * Si2frac_zint[it + (j+1)*doff];
      Si3frac_dint[it] = flow * Si3frac_zint[it + j*doff] + fhi * Si3frac_zint[it + (j+1)*doff];
      Si4frac_dint[it] = flow * Si4frac_zint[it + j*doff] + fhi * Si4frac_zint[it + (j+1)*doff];
    }
  
  /* log(T/K) interpolation */
  logT = dmin(dmax(logT, ttab_min), ttab_max-1.0e-6); /* protect boundaries */
  double ttab_bin_inv = ((double)n_ttab - 1.0) / (ttab_max - ttab_min);
  
  t    = (logT - ttab_min ) * ttab_bin_inv;
  j    = floor(t);
  fhi  = t - j;
  flow = 1 - fhi;
  
  *Si2frac = pow(10.0, flow * Si2frac_dint[j] + fhi * Si2frac_dint[j+1]);
  *Si3frac = pow(10.0, flow * Si3frac_dint[j] + fhi * Si3frac_dint[j+1]);
  *Si4frac = pow(10.0, flow * Si4frac_dint[j] + fhi * Si4frac_dint[j+1]);
  
  free(Si2frac_zint);
  free(Si3frac_zint);
  free(Si4frac_zint);
  
  free(Si2frac_dint);
  free(Si3frac_dint);
  free(Si4frac_dint);
}


/** \brief Initialises the external cloudy table at start up and gets
    the table boundaries (redshift, logT, log nH).
 */

void InitCloudy()
{
  char *fname;
  FILE *input;
  
  fname = "./cloudy_tables/tablesize_p19.dat";
  if(!(input=fopen(fname,"rb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  
  fread(&ztab_min,sizeof(double),1,input);
  fread(&ztab_max,sizeof(double),1,input);
  fread(&n_ztab,sizeof(int),1,input);
  
  fread(&ttab_min,sizeof(double),1,input);
  fread(&ttab_max,sizeof(double),1,input);
  fread(&n_ttab,sizeof(int),1,input);
  
  fread(&dtab_min,sizeof(double),1,input);
  fread(&dtab_max,sizeof(double),1,input);
  fread(&n_dtab,sizeof(int),1,input);
    
  fclose(input);
  
  InitIonTableMemory();
  
  read_iontable("./cloudy_tables/cloudytable_p19.dat");  
}


/** \brief Reads in the external cloudy table with the ion abundances
 *   
 *  \param fname The file name
 */

void read_iontable(char *fname)
{
  FILE *input; 

  if(!(input=fopen(fname,"rb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  fread(Si2frac_tab,sizeof(double),n_dtab*n_ttab*n_ztab,input);
  fread(Si3frac_tab,sizeof(double),n_dtab*n_ttab*n_ztab,input);
  fread(Si4frac_tab,sizeof(double),n_dtab*n_ttab*n_ztab,input);
  fclose(input);
}


/** \brief Initialises the memory for the external cloudy table
 */

void InitIonTableMemory(void)
{  
  Si2frac_tab = (double *)calloc(n_ttab*n_dtab*n_ztab, sizeof(double));
  if(NULL==Si2frac_tab){free(Si2frac_tab); printf("Memory allocation failed.\n"); exit(0);}
  
  Si3frac_tab = (double *)calloc(n_ttab*n_dtab*n_ztab, sizeof(double));
  if(NULL==Si3frac_tab){free(Si3frac_tab); printf("Memory allocation failed.\n"); exit(0);}

  Si4frac_tab = (double *)calloc(n_ttab*n_dtab*n_ztab, sizeof(double));
  if(NULL==Si4frac_tab){free(Si4frac_tab); printf("Memory allocation failed.\n"); exit(0);}
}

#endif
