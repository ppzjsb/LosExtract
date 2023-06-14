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
double *rho_wk_H_hires, *velaxis_hires;
double *rho_wk_H1_hires, *vel_wk_H1_hires, *temp_wk_H1_hires;
#ifdef HE2LYA
double *rho_wk_He2, *vel_wk_He2, *temp_wk_He2,*tau_He2;
double *rho_wk_He2_hires, *vel_wk_He2_hires, *temp_wk_He2_hires;
#endif
double *xlos, *ylos, *zlos;

double ztime_file;

int *ilos, nbins, nlos;
int hires;

#ifdef TAU_WEIGHT
double *rho_tau_H1, *temp_tau_H1;
#ifdef HE2LYA
double *rho_tau_He2, *temp_tau_He2;
#endif
#endif

#ifdef QUICK_GAUSS
static double *gauss;
static double dgx_inv;
#endif




/** \brief Calculate the Lyman-alpha optical depths for each sight-line  
 *
 * \param argc number of command line arguments passed
 * \param argv pointer array to each argument passed
 *
 * \return (0)                                       
 */

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
      free_los_hires_memory();

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




/** \brief Calculate the Lyman-alpha optical depths for each sight-line
 */

void compute_absorption()
{

#ifdef REBIN_KERNEL
  double vth = 1.0e-5 * sqrt(BOLTZMANN * 1.0e4 / (HMASS * AMU)); /* in km/s, assumes T=10^4 K */
#ifdef HE2LYA
  double vth_He2 = 1.0e-5 * sqrt(BOLTZMANN * 1.0e4 /(HEMASS * AMU));
  vth = dmin(vth, vth_He2);  
#endif
  double dvbin = velaxis[1] - velaxis[0];
  if(dvbin > vth)
    {
      printf("\nCurrent pixel size is %f km/s\n",dvbin);
      printf("The thermal lines widths are likely to be unresolved!\n");
      hires = ceil(dvbin/vth);
      printf("Increasing sample rate by a factor %d, to nbins=%d\n",hires,hires*nbins);
      printf("New pixel size is %f km/s\n\n",dvbin/(double)hires);
    }
  else
    hires = 1;  /* if hires=1, there is no resampling */
#else
  hires = 1; 
#endif
  
  
  double u_H1[nbins*hires], b_H1_inv[nbins*hires], b2_H1_inv[nbins*hires], tau_H1_dR[nbins*hires];
#ifdef TAU_WEIGHT
  double rho_wk_H_tauw_H1[nbins*hires], temp_wk_H1_tauw[nbins*hires];
#endif
  
#ifdef HE2LYA
  double u_He2[nbins*hires], b_He2_inv[nbins*hires], b2_He2_inv[nbins*hires], tau_He2_dR[nbins*hires];
#ifdef TAU_WEIGHT
  double rho_wk_H_tauw_He2[nbins*hires], temp_wk_He2_tauw[nbins*hires];
#endif
#endif
  
  double atime  = 1.0/(1.0+ztime);
  double rscale = (KPC*atime)/h100;     /* comoving kpc/h to cm */
  double escale = 1.0e10;               /* (km s^-1)^2 to (cm s^-1)^2 */
  double drbin  = box100/(double)(nbins*hires); /* comoving kpc/h */
  
  double Hz = 100.0*h100*sqrt(omegam/(atime*atime*atime)+omegal); /* km s^-1 Mpc^-1 */
  double H0 = 1.0e7 / MPC; /* 100 km s^-1 Mpc^-1 in cgs */ 
  
  double vmax  = box100 * Hz * rscale / MPC; /* box size km s^-1 */
  double vmax2 = 0.5 * vmax;
  
  double rhoc  = 3.0 * (H0*h100) * (H0*h100) / (8.0 * PI * GRAVITY); /* g cm^-3 */
  double critH = rhoc * omegab * Xh / (atime * atime * atime); /* g cm^-3*/
  
  double sigma_Lya_H1 = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_H1 * FOSC_LYA; /* cm^2 */
  double k1_conv      = 2.0 * BOLTZMANN / (HMASS * AMU);
  double k2_conv      = sigma_Lya_H1 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);

#ifdef VOIGT
  double aa_H1[nbins*hires];
  double k_voigt = GAMMA_LYA * LAMBDA_LYA_H1 / (4.0 * PI * sqrt(PI));  /* cm s^-1 */
#endif  
  
#ifdef HE2LYA
  double sigma_Lya_He2 = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_HE2 * FOSC_LYA; /* cm^2 */
  double k1_conv_He2   = 2.0 * BOLTZMANN / (HEMASS * AMU);  
  double k2_conv_He2   = sigma_Lya_He2 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
#ifdef VOIGT
  double aa_He2[nbins*hires];
  double k_voigt_He2 = GAMMA_LYA * LAMBDA_LYA_HE2 / (4.0 * PI * sqrt(PI));  /* cm s^-1 */
#endif
#endif

  double pccount = 0.0;
  
  
  allocate_los_hires_memory();
  resample_los(vmax);

  
  int ilos = 0;
  
  for(ilos=0; ilos<nlos; ilos++)
    {
      
      int iconv = 0;
      
#ifdef OPENMP      
#pragma omp parallel for
#endif
      for(iconv=0; iconv<nbins*hires; iconv++)
  	{
  	  long long convol_index =  iconv + nbins*hires*ilos;
	  
#ifdef SELF_SHIELD
	  
	  double nHcgs = critH * rho_wk_H_hires[convol_index] / (HMASS * AMU); /* cm^-3 */
	  
	  if (nHcgs >= n0_z)
	    rho_wk_H1_hires[convol_index] = self_shield(nHcgs, log10(temp_wk_H1_hires[convol_index]));

#endif

	  
#ifdef TEST_KERNEL
	  rho_wk_H1_hires[convol_index]  = 1.0e-5;
	  temp_wk_H1_hires[convol_index] = 1.0e3;
#ifdef HE2LYA
	  rho_wk_He2_hires[convol_index]  = 1.0e-4;
	  temp_wk_He2_hires[convol_index] = 1.0e3;
#endif
#endif

	  
#ifdef NO_PECVEL
  	  u_H1[iconv] = velaxis_hires[iconv]; 
#else
	  u_H1[iconv] = velaxis_hires[iconv] + vel_wk_H1_hires[convol_index]; /* km s^-1 */
#endif
	  
  	  b_H1_inv[iconv]  = 1.0 / sqrt(k1_conv * temp_wk_H1_hires[convol_index]); /* (cm s^-1)^-1 */
	  b2_H1_inv[iconv] = escale * b_H1_inv[iconv] * b_H1_inv[iconv];  /* (km s^-1)^-2 */
	  tau_H1_dR[iconv] = k2_conv * rho_wk_H_hires[convol_index] * rho_wk_H1_hires[convol_index] * b_H1_inv[iconv];
	  
#ifdef VOIGT
	  aa_H1[iconv] = k_voigt * b_H1_inv[iconv];
#endif	  
	  
#ifdef TAU_WEIGHT
	  rho_wk_H_tauw_H1[iconv] = rho_wk_H_hires[convol_index];
	  temp_wk_H1_tauw[iconv]  = temp_wk_H1_hires[convol_index];
#endif
	  
	  
	  
#ifdef HE2LYA	  
#ifdef NO_PECVEL
  	  u_He2[iconv] = velaxis_hires[iconv]; 
#else
	  u_He2[iconv] = velaxis_hires[iconv] + vel_wk_He2_hires[convol_index]; /* km s^-1 */
#endif
	  
  	  b_He2_inv[iconv]  = 1.0 / sqrt(k1_conv_He2 * temp_wk_He2_hires[convol_index]); /* (cm s^-1)^-1 */
	  b2_He2_inv[iconv] = escale * b_He2_inv[iconv] * b_He2_inv[iconv];  /* (km s^-1)^-2 */
	  tau_He2_dR[iconv] = k2_conv_He2 * rho_wk_H_hires[convol_index] * rho_wk_He2_hires[convol_index] * b_He2_inv[iconv];
	  
#ifdef VOIGT
	  aa_He2[iconv] = k_voigt_He2 * b_He2_inv[iconv];
#endif	  

#ifdef TAU_WEIGHT
	  rho_wk_H_tauw_He2[iconv] = rho_wk_H_hires[convol_index];
	  temp_wk_He2_tauw[iconv]  = temp_wk_He2_hires[convol_index];
#endif
#endif
	  
	}
      
      int ipix = 0;
      
      for(ipix=0; ipix<nbins*hires; ipix++)
  	{
	  long long pixel_index = ipix + nbins*hires*ilos;
	  
	  double tau_H1_sum = 0.0;
#ifdef TAU_WEIGHT
	  double rho_tau_H1_sum  = 0.0;
	  double temp_tau_H1_sum = 0.0;	  
#endif
	  
#ifdef HE2LYA
	  double tau_He2_sum = 0.0;
#ifdef TAU_WEIGHT	  
	  double rho_tau_He2_sum  = 0.0;
	  double temp_tau_He2_sum = 0.0;
#endif
#endif
	  
	  int iconv = 0;
	  
#ifdef OPENMP
#ifdef TAU_WEIGHT
	  
#ifdef HE2LYA
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum, tau_He2_sum, rho_tau_He2_sum, temp_tau_He2_sum)
#else
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum)	  
#endif
	  
#else
	  
#ifdef HE2LYA
#pragma omp parallel for reduction(+:tau_H1_sum, tau_He2_sum)
#else
#pragma omp parallel for reduction(+:tau_H1_sum)
#endif
	  
#endif	  
#endif
	  for(iconv=0; iconv<nbins*hires; iconv++)
  	    {
	      
	      double vdiff_H1 = fabs(velaxis_hires[ipix] - u_H1[iconv]);
	      if (vdiff_H1 > vmax2) vdiff_H1 = vmax - vdiff_H1;
	      
  	      double VH1_0  = vdiff_H1 * vdiff_H1 * b2_H1_inv[iconv]; 
	      
#ifdef QUICK_GAUSS
	      double t    = VH1_0 * dgx_inv;
	      int tint    = (int)t;
	      double fhi  = t - tint;
	      double flow = 1 - fhi;
	      
	      double VH1_1  = VH1_0 < GXMAX ? flow*gauss[tint]+fhi*gauss[tint+1] : 0.0;
#else
	      double VH1_1  = exp(-VH1_0);
#endif
	      
#ifdef VOIGT
	      /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025.
		 See their footnote 4 in the 2007 erratum.  We keep
		 the Gaussian profile close to the line centre, and
		 note the factor of 1/sqrt(PI) has been absorbed into
		 aa_H1. */
	      
	      double VH1_2 = 1.5 / VH1_0;
	      
	      double profile_H1 = (VH1_0 < 1.0e-6) 
		? VH1_1 : VH1_1 - aa_H1[iconv] / VH1_0 * (VH1_1*VH1_1 * (4.0*VH1_0*VH1_0 + 7.0*VH1_0 + 4.0 + VH1_2) - VH1_2 - 1.0);
	      
#else
	      double profile_H1 = VH1_1;
#endif
	      tau_H1_sum += tau_H1_dR[iconv] * profile_H1;
	      
#ifdef TAU_WEIGHT
	      rho_tau_H1_sum  += rho_wk_H_tauw_H1[iconv] * tau_H1_dR[iconv] * profile_H1;
	      temp_tau_H1_sum += temp_wk_H1_tauw[iconv] * tau_H1_dR[iconv] * profile_H1;
#endif

#ifdef HE2LYA	      
	      double vdiff_He2 = fabs(velaxis_hires[ipix] - u_He2[iconv]);
	      if (vdiff_He2 > vmax2) vdiff_He2 = vmax - vdiff_He2;
	      
  	      double VHe2_0  = vdiff_He2 * vdiff_He2 * b2_He2_inv[iconv]; 
	      
#ifdef QUICK_GAUSS
	      double t    = VHe2_0 * dgx_inv;
	      int tint    = (int)t;
	      double fhi  = t - tint;
	      double flow = 1 - fhi;
	      
	      double VHe2_1  = VHe2_0 < GXMAX ? flow*gauss[tint]+fhi*gauss[tint+1] : 0.0;
#else
	      double VHe2_1  = exp(-VHe2_0);
#endif
	      
#ifdef VOIGT
	      double VHe2_2 = 1.5 / VHe2_0;
	      
	      double profile_He2 = (VHe2_0 < 1.0e-6) 
		? VHe2_1 : VHe2_1 - aa_He2[iconv] / VHe2_0 * (VHe2_1*VHe2_1 * (4.0*VHe2_0*VHe2_0 + 7.0*VHe2_0 + 4.0 + VHe2_2) - VHe2_2 - 1.0);
	      
#else
	      double profile_He2 = VHe2_1;
#endif
	      tau_He2_sum += tau_He2_dR[iconv] * profile_He2;
	      
#ifdef TAU_WEIGHT
	      rho_tau_He2_sum  += rho_wk_H_tauw_He2[iconv] * tau_He2_dR[iconv] * profile_He2;
	      temp_tau_He2_sum += temp_wk_He2_tauw[iconv]  * tau_He2_dR[iconv] * profile_He2;
#endif
#endif

	    }

	  if(pixel_index % hires == 0)
	    {
	      int base_index = floor(pixel_index / hires);
	   
	      tau_H1[base_index] = tau_H1_sum;
	      
#ifdef TAU_WEIGHT	  
	      rho_tau_H1[base_index]  = rho_tau_H1_sum / tau_H1[base_index];
	      temp_tau_H1[base_index] = temp_tau_H1_sum / tau_H1[base_index];
#endif
	      
#ifdef HE2LYA
	      tau_He2[base_index] = tau_He2_sum;
	      
#ifdef TAU_WEIGHT	  
	      rho_tau_He2[base_index]  = rho_tau_He2_sum / tau_He2[base_index];
	      temp_tau_He2[base_index] = temp_tau_He2_sum / tau_He2[base_index];
#endif  
#endif
	    }
	  
  	}
      
      double pcdone = 100.0*(double)ilos/((double)nlos-1.0);
      if(pcdone >= pccount)
  	{
  	  printf("%3.2f%%\n",pcdone);
  	  pccount += 10.0;
  	}
    } 
}





/** \brief Resamples the line of sight data if the thermal broadening
 *  kernel is judged to be under-resolved at the native resolution.
 *   Requires the RESAMPLE_KERNEL flag to do anything, otherwise it
 *   just defaults to the native sample rate.
 *   
 *  \param vmax The box size in km/s
 *
 */

void resample_los(double vmax)
{

  int ilos = 0, ipix = 0;
  
  if(hires > 1)
    {
      double dvbin_hires = vmax/(double)(nbins*hires);

      for(ipix=0; ipix<nbins*hires - 1; ipix++)
	velaxis_hires[ipix+1] = velaxis_hires[ipix] + dvbin_hires;
      
      for(ilos=0; ilos<nlos; ilos++)
	for(ipix=0; ipix<nbins*hires; ipix++)
	  {
	    long long pixel_index = ipix + nbins*hires*ilos;
	    
	    long long base_index_lo = floor(pixel_index / hires);
	    long long base_index_hi = base_index_lo < (ilos+1)*nbins-1 ?  base_index_lo+1 : base_index_lo-nbins+1;
	    
	    int base_vel_index_lo = base_index_lo - ilos*nbins;
	    int base_vel_index_hi = base_vel_index_lo < nbins-1 ?  base_vel_index_lo+1 : base_vel_index_lo-nbins+1;
	    
	    
	    rho_wk_H_hires[pixel_index] = lerp(velaxis_hires[ipix],
					       velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
					       rho_wk_H[base_index_lo], rho_wk_H[base_index_hi]);
	    
	    rho_wk_H1_hires[pixel_index] = lerp(velaxis_hires[ipix],
						velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
						rho_wk_H1[base_index_lo], rho_wk_H1[base_index_hi]);
	    
	    vel_wk_H1_hires[pixel_index] = lerp(velaxis_hires[ipix],
						velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
						vel_wk_H1[base_index_lo], vel_wk_H1[base_index_hi]);
	    
	    temp_wk_H1_hires[pixel_index] = lerp(velaxis_hires[ipix],
						 velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
						 temp_wk_H1[base_index_lo], temp_wk_H1[base_index_hi]);
	    
#ifdef HE2LYA
	    rho_wk_He2_hires[pixel_index] = lerp(velaxis_hires[ipix],
						 velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
						 rho_wk_He2[base_index_lo], rho_wk_He2[base_index_hi]);
	    
	    vel_wk_He2_hires[pixel_index] = lerp(velaxis_hires[ipix],
						 velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
						 vel_wk_He2[base_index_lo], vel_wk_He2[base_index_hi]);
	    
	    temp_wk_He2_hires[pixel_index] = lerp(velaxis_hires[ipix],
						  velaxis[base_vel_index_lo], velaxis[base_vel_index_hi],
						  temp_wk_He2[base_index_lo], temp_wk_He2[base_index_hi]);
#endif
	  }
    }
  else
    {
      for(ilos=0; ilos<nlos; ilos++)
	for(ipix=0; ipix < nbins*hires; ipix++)
	  {
	    
	    long long pixel_index = ipix + nbins*hires*ilos;
	    
	    rho_wk_H_hires[pixel_index]   = rho_wk_H[pixel_index];
	    
	    rho_wk_H1_hires[pixel_index]  = rho_wk_H1[pixel_index];
	    vel_wk_H1_hires[pixel_index]  = vel_wk_H1[pixel_index];
	    temp_wk_H1_hires[pixel_index] = temp_wk_H1[pixel_index];
	    
#ifdef HE2LYA
	    rho_wk_He2_hires[pixel_index]  = rho_wk_He2[pixel_index];
	    vel_wk_He2_hires[pixel_index]  = vel_wk_He2[pixel_index];
	    temp_wk_He2_hires[pixel_index] = temp_wk_He2[pixel_index];
#endif
	    if(ilos == 0)
	      velaxis_hires[ipix] = velaxis[ipix];
	    
	  }
    }
}




/** \brief Read in the LOS binary file.
 */

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
  fread(&Xh,sizeof(double),1,input); 
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
  
#ifdef HE2LYA
  fread(rho_wk_He2,sizeof(double),nbins*nlos,input); /* n_HeII/n_H */
  fread(temp_wk_He2,sizeof(double),nbins*nlos,input); /* T [K], HeII weighted */
  fread(vel_wk_He2,sizeof(double),nbins*nlos,input); /* v_pec [km s^-1], HeII weighted */
#endif
  
  fclose(input);
}




/** \brief Write the optical depths to a binary file.
 */

void write_tau()
{
  char fname[400];
  FILE *output;
   
  sprintf(fname, "%s/tauH1_x%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
 
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
  sprintf(fname, "%s/tauwH1_x%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  fwrite(rho_tau_H1,sizeof(double),nbins*nlos,output);
  fwrite(temp_tau_H1,sizeof(double),nbins*nlos,output);
  fclose(output);
#endif  

  
#ifdef HE2LYA  
  sprintf(fname, "%s/tauHe2_x%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  fwrite(tau_He2,sizeof(double),nbins*nlos,output);
  fclose(output);

  
#ifdef TAU_WEIGHT
  sprintf(fname, "%s/tauwHe2_x%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(0);
    }
  fwrite(rho_tau_He2,sizeof(double),nbins*nlos,output);
  fwrite(temp_tau_He2,sizeof(double),nbins*nlos,output);
  fclose(output);
#endif
  
#endif
  

}




/** \brief Make a look-up table for the Gaussian line profile.
 */

#ifdef QUICK_GAUSS

void GaussianProfileTable()
{
  
  gauss = (double *)calloc(NGXTAB+1, sizeof(double));

  if(NULL==gauss){free(gauss); printf("Memory allocation failed.\n"); exit(0);}
  
  double dgx =  GXMAX / NGXTAB;
  
  dgx_inv = 1.0 / dgx;
  
  int i = 0;

  for(i = 0; i <= NGXTAB; i++)
    {
      double gx = dgx * (double)i;
      gauss[i] = exp(-gx);
    }
}

#endif




/** \brief Allocate LOS memory.  Would be neater to use a structure.
 */

void allocate_los_memory()
{
  
  ilos = (int *)calloc(nlos, sizeof(int));
  if(NULL==ilos){free(ilos); printf("Memory allocation failed.\n"); exit(0);}
  
  xlos = (double *)calloc(nlos, sizeof(double));
  if(NULL==xlos){free(xlos); printf("Memory allocation failed.\n"); exit(0);}

  ylos = (double *)calloc(nlos, sizeof(double));
  if(NULL==ylos){free(ylos); printf("Memory allocation failed.\n"); exit(0);}

  zlos = (double *)calloc(nlos, sizeof(double));
  if(NULL==zlos){free(zlos); printf("Memory allocation failed.\n"); exit(0);}
  
  posaxis = (double *)calloc(nbins, sizeof(double));
  if(NULL==posaxis){free(posaxis); printf("Memory allocation failed.\n"); exit(0);}

  velaxis = (double *)calloc(nbins, sizeof(double));
  if(NULL==velaxis){free(velaxis); printf("Memory allocation failed.\n"); exit(0);}
  
  rho_wk_H = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_H){free(rho_wk_H); printf("Memory allocation failed.\n"); exit(0);}
  
  rho_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_H1){free(rho_wk_H1); printf("Memory allocation failed.\n"); exit(0);}
  
  vel_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==vel_wk_H1){free(vel_wk_H1); printf("Memory allocation failed.\n"); exit(0);}
  
  temp_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_wk_H1){free(temp_wk_H1); printf("Memory allocation failed.\n"); exit(0);}
  
  tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_H1){free(tau_H1); printf("Memory allocation failed.\n"); exit(0);}
  
  
#ifdef TAU_WEIGHT
  rho_tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_tau_H1){free(rho_tau_H1); printf("Memory allocation failed.\n"); exit(0);}
  
  temp_tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_tau_H1){free(temp_tau_H1); printf("Memory allocation failed.\n"); exit(0);}
#endif  

  
#ifdef HE2LYA
  rho_wk_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_He2){free(rho_wk_He2); printf("Memory allocation failed.\n"); exit(0);}
  
  vel_wk_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==vel_wk_He2){free(vel_wk_He2); printf("Memory allocation failed.\n"); exit(0);}
  
  temp_wk_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_wk_He2){free(temp_wk_He2); printf("Memory allocation failed.\n"); exit(0);}
  
  tau_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_He2){free(tau_He2); printf("Memory allocation failed.\n"); exit(0);}
  
#ifdef TAU_WEIGHT
  rho_tau_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_tau_He2){free(rho_tau_He2); printf("Memory allocation failed.\n"); exit(0);}
  
  temp_tau_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_tau_He2){free(temp_tau_He2); printf("Memory allocation failed.\n"); exit(0);}
#endif
#endif

}




/** \brief Free LOS memory.
 */

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
  
#ifdef HE2LYA
  free(rho_wk_He2);
  free(vel_wk_He2);
  free(temp_wk_He2);
  free(tau_He2);
#ifdef TAU_WEIGHT
  free(rho_tau_He2);
  free(temp_tau_He2);
#endif
#endif
}




/** \brief Allocate LOS memory for increased sample rate.  Would be
    neater to use a structure.
 */

void allocate_los_hires_memory()
{
  velaxis_hires = (double *)calloc(hires*nbins, sizeof(double));
  if(NULL==velaxis_hires){free(velaxis_hires); printf("Memory allocation failed.\n"); exit(0);}
  
  rho_wk_H_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==rho_wk_H_hires){free(rho_wk_H_hires); printf("Memory allocation failed.\n"); exit(0);}
  
  rho_wk_H1_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==rho_wk_H1_hires){free(rho_wk_H1_hires); printf("Memory allocation failed.\n"); exit(0);}
  
  vel_wk_H1_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==vel_wk_H1_hires){free(vel_wk_H1_hires); printf("Memory allocation failed.\n"); exit(0);}
  
  temp_wk_H1_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==temp_wk_H1_hires){free(temp_wk_H1_hires); printf("Memory allocation failed.\n"); exit(0);}
  
#ifdef HE2LYA
  rho_wk_He2_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==rho_wk_He2_hires){free(rho_wk_He2_hires); printf("Memory allocation failed.\n"); exit(0);}
  
  vel_wk_He2_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==vel_wk_He2_hires){free(vel_wk_He2_hires); printf("Memory allocation failed.\n"); exit(0);}
  
  temp_wk_He2_hires = (double *)calloc(nlos*hires*nbins, sizeof(double));
  if(NULL==temp_wk_He2_hires){free(temp_wk_He2_hires); printf("Memory allocation failed.\n"); exit(0);}
#endif
}



/** \brief Free LOS memory for increased sample rate.
 */

void free_los_hires_memory()
{
  free(velaxis_hires);
  
  free(rho_wk_H_hires);

  free(rho_wk_H1_hires);
  free(vel_wk_H1_hires);
  free(temp_wk_H1_hires);

#ifdef HE2LYA
  free(rho_wk_He2_hires);
  free(vel_wk_He2_hires);
  free(temp_wk_He2_hires);
#endif
}


