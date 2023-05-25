
/********************************************************************************/

/**! \file extract_los.c 
 *
 * \brief Reads in an on-the-fly LOS output and computes the optical depths
 *
 */

/********************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef OPENMP
#include <omp.h>
#endif

#include "global_vars.h"
#include "parameters.h"
#include "proto.h"

char *path;

double *rho_wk_H, *posaxis, *velaxis;
double *rho_wk_H1, *vel_wk_H1, *temp_wk_H1, *tau_H1;
double *xlos, *ylos, *zlos;

double ztime_file;

int *ilos, nbins, nlos;

#ifdef TAU_WEIGHT
double *rho_tau_H1, *temp_tau_H1;
#endif

#ifdef QUICK_GAUSS
static double *gauss;
static double dgx_inv;
#endif

/********************************************************************************/

/** \brief Calculate the Lyman-alpha optical depths for each sight-line  
 *
 * \param argc number of command line arguments passed
 * \param argv pointer array to each argument passed
 *
 * \return (0)                                       
 */

/********************************************************************************/


int main(int argc, char **argv)
{

#ifdef OPENMP
  double start;
#pragma omp parallel 
  start = omp_get_wtime();
#else
  double start = clock();
#endif

  
  if(argc != 2)
    {
      fprintf(stderr, "\n\nIncorrect argument(s).  Specify:\n\n");
      fprintf(stderr, "<path>     (path to output files)\n");
      fprintf(stderr, "Usage: ./LosExtract <path>\n");
      fprintf(stderr, "\n\n");
      exit(0);
    }
  
  path = argv[1];

  ztime_file = ZI_LOS;

#ifdef QUICK_GAUSS  
  GaussianProfileTable();
#endif

#ifdef OPENMP  
  int nthreads;    
#pragma omp parallel 
  nthreads = omp_get_num_threads();
  printf("\nRunning with OpenMP on %d threads\n",nthreads);
#endif


  int i = 0;
  for(i=0; i < NLOSFILES; i++)
    {
      printf("\nReading los output at z=%.3f:\n",ztime_file);
      read_los();

#ifdef SELF_SHIELD  
      InitCool();
#endif
      
      printf("\nComputing optical depths...\n");
      compute_absorption();
      printf("Done.\n\n");
      
      write_tau();
      
      free_los_memory();

      ztime_file += DZ_LOS;
      
    }

#ifdef OPENMP
  double end;
#pragma omp parallel
  end = omp_get_wtime();
  printf("\nTotal wall time %lf s\n",end - start);
#else
  double end = clock();
  printf("\nTotal CPU time %lf s\n", (end - start)/CLOCKS_PER_SEC);
#endif
  
  return(0);
}

/********************************************************************************/

/** \brief Calculate the Lyman-alpha optical depths for each sight-line                                                 
 */

/********************************************************************************/

void compute_absorption()
{
  
  double u_H1[nbins], b_H1_inv[nbins], b2_H1_inv[nbins], tau_H1_dR[nbins];
  
#ifdef TAU_WEIGHT
  double rho_wk_H_tauw[nbins], temp_wk_H1_tauw[nbins];
#endif
  
  double atime  = 1.0/(1.0+ztime);
  double rscale = (KPC*atime)/h100;     /* comoving kpc/h to cm */
  double escale = 1.0e10;               /* (km s^-1)^2 to (cm s^-1)^2 */
  double drbin  = box100/(double)nbins; /* comoving kpc/h */
  
  double Hz = 100.0*h100*sqrt(omegam/(atime*atime*atime)+omegal); /* km s^-1 Mpc^-1 */
  double H0 = 1.0e7 / MPC; /* 100 km s^-1 Mpc^-1 in cgs */ 
  
  double vmax  = box100 * Hz * rscale / MPC; /* box size km s^-1 */
  double vmax2 = 0.5 * vmax;
 
  double rhoc  = 3.0 * (H0*h100) * (H0*h100) / (8.0 * PI * GRAVITY); /* g cm^-3 */
  double critH = rhoc * omegab * Xh / (atime * atime * atime); /* g cm^-3*/

  double sigma_Lya = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_H1 * FOSC_LYA; /* cm^2 */
  double k1_conv   = 2.0 * BOLTZMANN / (HMASS * AMU);
  double k2_conv   = sigma_Lya * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
  
#ifdef VOIGT
  double aa_H1[nbins];
  double k_voigt   = GAMMA_LYA_H1 * LAMBDA_LYA_H1 / (4.0 * PI * sqrt(PI));  /* cm s^-1 */
#endif  

  double pccount = 0.0;
  
  int iproc = 0;
  
  for(iproc=0; iproc<nlos; iproc++)
    {
      
      int j = 0;

#ifdef OPENMP      
#pragma omp parallel for
#endif
      for(j=0; j<nbins; j++)
  	{
  	  int convol_index =  j + nbins*iproc;
	  
#ifdef SELF_SHIELD
	  
	  double nHcgs = critH * rho_wk_H[convol_index] / (HMASS * AMU); /* cm^-3 */
	  
	  if (nHcgs >= n0_z)
	    rho_wk_H1[convol_index] = self_shield(nHcgs, log10(temp_wk_H1[convol_index]));
#endif

	  
#ifdef TEST_KERNEL
	  rho_wk_H1[convol_index]  = 1.0e-5;
	  temp_wk_H1[convol_index] = 1.0e3;
#endif

	  
#ifdef NO_PECVEL
  	  u_H1[j] = velaxis[j]; 
#else
	  u_H1[j] = velaxis[j] + vel_wk_H1[convol_index]; /* km s^-1 */
#endif
	  
  	  b_H1_inv[j]  = 1.0 / sqrt(k1_conv*temp_wk_H1[convol_index]); /* (cm s^-1)^-1 */
	  b2_H1_inv[j] = escale * b_H1_inv[j] * b_H1_inv[j];  /* (km s^-1)^-2 */
	  tau_H1_dR[j] = k2_conv * rho_wk_H[convol_index] * rho_wk_H1[convol_index] * b_H1_inv[j];
	  
#ifdef VOIGT
	  aa_H1[j] = k_voigt * b_H1_inv[j];
#endif	  

#ifdef TAU_WEIGHT
	  rho_wk_H_tauw[j]   = rho_wk_H[convol_index];
	  temp_wk_H1_tauw[j] = temp_wk_H1[convol_index];
#endif
	  
  	}
      
      int i = 0;
      
      for(i=0; i<nbins; i++)
  	{
	  int pixel_index = i + nbins*iproc;
	  
	  double tau_H1_sum = 0.0;
	  
#ifdef TAU_WEIGHT
	  double rho_tau_H1_sum  = 0.0;
	  double temp_tau_H1_sum = 0.0;
#endif

	  int j = 0;
	  
#ifdef OPENMP
#ifdef TAU_WEIGHT
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum)
#else
#pragma omp parallel for reduction(+:tau_H1_sum) 
#endif
#endif
	  for(j=0; j<nbins; j++)
  	    {
	      
	      double vdiff_H1 = fabs(velaxis[i] - u_H1[j]);
	      if (vdiff_H1 > vmax2) vdiff_H1 = vmax - vdiff_H1;
	      
  	      double VH0  = vdiff_H1 * vdiff_H1 * b2_H1_inv[j]; 

#ifdef QUICK_GAUSS
	      double t    = VH0 * dgx_inv;
	      int tint    = (int)t;
	      double fhi  = t - tint;
	      double flow = 1 - fhi;
	      
	      double VH1  = VH0 < GXMAX ? flow*gauss[tint]+fhi*gauss[tint+1] : 0.0;
#else
	      double VH1  = exp(-VH0);
#endif
	      
#ifdef VOIGT
	      /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025.
		 See their footnote 4 in the 2007 erratum.  Note below
		 we still use the Gaussian profile close to the line
		 centre, and the factor of 1/sqrt(PI) has been
		 absorbed into aa_H1 */
	      
	      double VH2 = 1.5 / VH0;
	      
	      double profile_H1 = (VH0 < 1.0e-6) 
		? VH1 : VH1 - aa_H1[j] / VH0 * (VH1*VH1 * (4.0*VH0*VH0 + 7.0*VH0 + 4.0 + VH2) - VH2 - 1.0);
	      
#else
	      double profile_H1 = VH1;
#endif
	      tau_H1_sum += tau_H1_dR[j] * profile_H1;
	      
#ifdef TAU_WEIGHT
	      /* HI optical depth weighted quantities */
	      rho_tau_H1_sum  += rho_wk_H_tauw[j]   * tau_H1_dR[j] * profile_H1;
	      temp_tau_H1_sum += temp_wk_H1_tauw[j] * tau_H1_dR[j] * profile_H1;
#endif
	    }
	  
	  tau_H1[pixel_index] = tau_H1_sum;
	  
	  
#ifdef TAU_WEIGHT	  
	  rho_tau_H1[pixel_index]  = rho_tau_H1_sum / tau_H1[pixel_index];
	  temp_tau_H1[pixel_index] = temp_tau_H1_sum / tau_H1[pixel_index];
#endif	  
  	}
      
      double pcdone = 100.0*(double)iproc/((double)nlos-1.0);
      if(pcdone >= pccount)
  	{
  	  printf("%3.2f%%\n",pcdone);
  	  pccount += 10.0;
  	}
    } 
}

/********************************************************************************/

/** \brief Read in the LOS binary file
 */

/********************************************************************************/

void read_los()
{
  FILE *input;
  
  char fname[400];
  
  sprintf(fname, "%s/%s_z%.3f.dat",path,LOSBASE,ztime_file);
  if(!(input=fopen(fname,"rb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  
  fread(&ztime,sizeof(double),1,input);
  fread(&omegam,sizeof(double),1,input);
  fread(&omegal,sizeof(double),1,input);
  fread(&omegab,sizeof(double),1,input);
  fread(&h100,sizeof(double),1,input);
  fread(&box100,sizeof(double),1,input);
  fread(&Xh,sizeof(double),1,input); /* hydrogen fraction by mass */
  fread(&nbins,sizeof(int),1,input);
  fread(&nlos,sizeof(int),1,input);
  
  printf("z      = %f\n",ztime);
  printf("omegaM = %f\n",omegam);
  printf("omegaL = %f\n",omegal);
  printf("omegab = %f\n",omegab);
  printf("h100   = %f\n",h100);
  printf("box100 = %f\n",box100);
  printf("Xh     = %f\n",Xh);
  printf("nbins  = %d\n",nbins);
  printf("nlos   = %d\n",nlos);
    
  if(fabs(ztime-ztime_file) > 1.0e-3)
    {
      printf("\nRedshift in the header does not match the file name\n");
      printf("Header: z=%3.4f, File: z=%3.4f\n",ztime,ztime_file);
      printf("Stopping...\n\n");
      exit(0);
    }
  
  allocate_los_memory();
  
  fread(ilos,sizeof(int),nlos,input); /* LOS axis */
  fread(xlos,sizeof(double),nlos,input); /* LOS positions, comoving kpc/h */
  fread(ylos,sizeof(double),nlos,input);
  fread(zlos,sizeof(double),nlos,input);
  fread(posaxis,sizeof(double),nbins,input); /* pixel positions, comoving kpc/h */
  fread(velaxis,sizeof(double),nbins,input); /* pixel positions, km s^-1 */
  fread(rho_wk_H,sizeof(double),nbins*nlos,input); /* gas overdensity, Delta=rho/rho_crit */
  
  fread(rho_wk_H1,sizeof(double),nbins*nlos,input); /* n_HI/n_H */
  fread(temp_wk_H1,sizeof(double),nbins*nlos,input); /* T [K], HI weighted */
  fread(vel_wk_H1,sizeof(double),nbins*nlos,input); /* v_pec [km s^-1], HI weighted */

  fclose(input);
}

/********************************************************************************/

/** \brief Write the optical depths to a binary file
 */

/********************************************************************************/

void write_tau()
{
  char fname[400];
  FILE *output;
   
  sprintf(fname, "%s/tau%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
 
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }

  fwrite(tau_H1,sizeof(double),nbins*nlos,output);

#ifdef SELF_SHIELD
  fwrite(rho_wk_H1,sizeof(double),nbins*nlos,output);
#endif
  
  fclose(output);
  
#ifdef TAU_WEIGHT
  sprintf(fname, "%s/tauw%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  fwrite(rho_tau_H1,sizeof(double),nbins*nlos,output);
  fwrite(temp_tau_H1,sizeof(double),nbins*nlos,output);
  fclose(output);
#endif  


}

/********************************************************************************/

/** \brief Make a look-up table for the Gaussian line profile 
 */

/********************************************************************************/

#ifdef QUICK_GAUSS

void GaussianProfileTable()
{
  
  gauss = (double *)calloc(NGXTAB+1, sizeof(double));

  if(NULL==gauss)
    {
      free(gauss);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  double dgx =  GXMAX / NGXTAB;
  
  dgx_inv = 1.0 / dgx;
  
  int i = 0;

#ifdef OPENMP      
#pragma omp parallel for
#endif
  for(i = 0; i <= NGXTAB; i++)
    {
      double gx = dgx * (double)i;
      gauss[i] = exp(-gx);
    }
}

#endif

/********************************************************************************/

/** \brief Allocate LOS memory
 */

/********************************************************************************/

void allocate_los_memory()
{
  
  ilos = (int *)calloc(nlos, sizeof(int));
  if(NULL==ilos)
    {
      free(ilos);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  xlos = (double *)calloc(nlos, sizeof(double));
  if(NULL==xlos)
    {
      free(xlos);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }

  ylos = (double *)calloc(nlos, sizeof(double));
  if(NULL==ylos)
    {
      free(ylos);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }

  zlos = (double *)calloc(nlos, sizeof(double));
  if(NULL==zlos)
    {
      free(zlos);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  posaxis = (double *)calloc(nbins, sizeof(double));
  if(NULL==posaxis)
    {
      free(posaxis);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }

  velaxis = (double *)calloc(nbins, sizeof(double));
  if(NULL==velaxis)
    {
      free(velaxis);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  rho_wk_H = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_H)
    {
      free(rho_wk_H);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  rho_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_H1)
    {
      free(rho_wk_H1);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  vel_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==vel_wk_H1)
    {
      free(vel_wk_H1);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  temp_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_wk_H1)
    {
      free(temp_wk_H1);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_H1)
    {
      free(tau_H1);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }

#ifdef TAU_WEIGHT
  
  rho_tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_tau_H1)
    {
      free(rho_tau_H1);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
  temp_tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_tau_H1)
    {
      free(temp_tau_H1);
      printf("Memory allocation failed in extract_spectra.c\n");
      exit(0);
    }
  
#endif

}

/********************************************************************************/

/** \brief Free LOS memory
 */

/********************************************************************************/

void free_los_memory()
{
  
  free(ilos);
  free(xlos);
  free(ylos);
  free(zlos);

  free(posaxis);
  free(velaxis);
  
  free(rho_wk_H);
  free(rho_wk_H1);
  free(vel_wk_H1);
  free(temp_wk_H1);

  free(tau_H1);
  
#ifdef TAU_WEIGHT
  free(rho_tau_H1);
  free(temp_tau_H1);
#endif

}

/********************************************************************************/
