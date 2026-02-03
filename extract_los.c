
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

double *posaxis, *velaxis;
double *rho_wk_H;
double *rho_wk_H1, *vel_wk_H1, *temp_wk_H1, *tau_H1;

#ifdef HE2LYA
double *rho_wk_He2, *vel_wk_He2, *temp_wk_He2,*tau_He2;
#endif

#ifdef TAU_WEIGHT
double *rho_tau_H1, *temp_tau_H1;
#ifdef HE2LYA
double *rho_tau_He2, *temp_tau_He2;
#endif
#endif

#ifdef SILICON
double *tau_1190_Si2, *tau_1193_Si2, *tau_1260_Si2, *tau_1207_Si3, *tau_1394_Si4, *tau_1403_Si4;
#endif

double *xlos, *ylos, *zlos;
double ztime_file;

int *ilos, nbins, nlos;
int nhires;

#ifdef QUICK_LINE
static double *exp_negx;
static double dx_inv;
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
      exit(1);
    }
  
  path = argv[1];

  ztime_file = ZI_LOS;

#ifdef QUICK_LINE  
  LineProfileTable();
#endif

#ifdef SILICON
  InitCloudy();
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




/** \brief Calculate the Lyman-alpha optical depths for each sight-line
 */

void compute_absorption()
{

#ifdef RESAMPLE

  double Tth = 1.0e3;
  double vth = 1.0e-5 * sqrt(BOLTZMANN * Tth / (HMASS * AMU)); /* in km/s, assumes T=10^4 K */
  
#ifdef HE2LYA
  double vth_He2 = 1.0e-5 * sqrt(BOLTZMANN * Tth /(HEMASS * AMU));
  vth = dmin(vth, vth_He2);  
#endif
  
#ifdef SILICON
  double vth_Si = 1.0e-5 * sqrt(BOLTZMANN * Tth /(SIMASS * AMU));
  vth = dmin(vth, vth_Si);
#endif  

  double dvbin = velaxis[1] - velaxis[0];
  
  nhires = (dvbin > vth) ? ceil(dvbin / vth) : 1;
  
  printf("Resample %dx: n1=%d n2=%d dv1=%f dv2=%f\n\n", nhires, nbins, nbins*nhires, dvbin, dvbin/nhires);
  
#else
  nhires = 1; 
#endif
  
    
  double atime  = 1.0/(1.0+ztime);
  double rscale = (KPC*atime)/h100;     /* comoving kpc/h to cm */
  double escale = 1.0e10;               /* (km s^-1)^2 to (cm s^-1)^2 */
  double drbin  = box100/(double)(nbins*nhires); /* comoving kpc/h */
  
  double Hz    = 100.0*h100*sqrt(omegam/(atime*atime*atime)+omegal); /* km s^-1 Mpc^-1 */
  double H0    = 1.0e7 / MPC; /* 100 km s^-1 Mpc^-1 in cgs */   
  double vmax  = box100 * Hz * rscale / MPC; /* box size km s^-1 */
  double vmax2 = 0.5 * vmax;
  double rhoc  = 3.0 * (H0*h100) * (H0*h100) / (8.0 * PI * GRAVITY); /* g cm^-3 */
  double critH = rhoc * omegab * Xh / (atime * atime * atime); /* g cm^-3*/

    
    
  
  /* HI Lyman-alpha (1216) */
  double u_H1[nbins*nhires], b_H1, b2_H1[nbins*nhires], tau_H1_dR[nbins*nhires];
#ifdef TAU_WEIGHT
  double rho_wk_H_tauw_H1[nbins*nhires], temp_wk_H1_tauw[nbins*nhires];
#endif
  double sigma_Lya_H1 = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_H1 * FOSC_LYA_H1; /* cm^2 */
  double k1_conv      = 2.0 * BOLTZMANN / (HMASS * AMU);
  double k2_conv      = sigma_Lya_H1 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
#ifdef VOIGT
  double aa_H1[nbins*nhires];
  double k_voigt = GAMMA_LYA_H1 * LAMBDA_LYA_H1 / (4.0 * PI * sqrt(PI));  /* cm s^-1 */
#endif  

  
  /* HeII Lyman-alpha (304) */
#ifdef HE2LYA
  double u_He2[nbins*nhires], b_He2, b2_He2[nbins*nhires], tau_He2_dR[nbins*nhires];
#ifdef TAU_WEIGHT
  double rho_wk_H_tauw_He2[nbins*nhires], temp_wk_He2_tauw[nbins*nhires];
#endif
  double sigma_Lya_He2 = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_HE2 * FOSC_LYA_HE2; /* cm^2 */
  double k1_conv_He2   = 2.0 * BOLTZMANN / (HEMASS * AMU);  
  double k2_conv_He2   = sigma_Lya_He2 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
#ifdef VOIGT
  double aa_He2[nbins*nhires];
  double k_voigt_He2 = GAMMA_LYA_HE2 * LAMBDA_LYA_HE2 / (4.0 * PI * sqrt(PI));  /* cm s^-1 */
#endif
#endif
  
  
  /* SiII (1190,1193,1260), SiIII (1207) and SiIV (1394,1403) */
#ifdef SILICON
  double b_Si, b2_Si[nbins*nhires];

  double tau_1190_Si2_dR[nbins*nhires], tau_1193_Si2_dR[nbins*nhires];
  double tau_1260_Si2_dR[nbins*nhires], tau_1207_Si3_dR[nbins*nhires];
  double tau_1394_Si4_dR[nbins*nhires], tau_1403_Si4_dR[nbins*nhires];
  
  double Si2frac, Si3frac, Si4frac;

  /* Note: these should be in functions */
  double sigma_1190_Si2   = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_1190_Si2 * FOSC_1190_Si2;
  double sigma_1193_Si2   = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_1193_Si2 * FOSC_1193_Si2;
  double sigma_1260_Si2   = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_1260_Si2 * FOSC_1260_Si2;
  double sigma_1207_Si3   = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_1207_Si3 * FOSC_1207_Si3;
  double sigma_1394_Si4   = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_1394_Si4 * FOSC_1394_Si4;
  double sigma_1403_Si4   = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_1403_Si4 * FOSC_1403_Si4;

  double k1_conv_Si       = 2.0 * BOLTZMANN / (SIMASS * AMU);

  double k2_conv_1190_Si2 = sigma_1190_Si2 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
  double k2_conv_1193_Si2 = sigma_1193_Si2 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
  double k2_conv_1260_Si2 = sigma_1260_Si2 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
  double k2_conv_1207_Si3 = sigma_1207_Si3 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
  double k2_conv_1394_Si4 = sigma_1394_Si4 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
  double k2_conv_1403_Si4 = sigma_1403_Si4 * C * rscale * drbin * critH / (sqrt(PI) * HMASS * AMU);
#endif


  double velaxis_hires[nbins*nhires];
  double dvbin_hires = vmax / (nbins*nhires);
  int i = 0;
  for(i=0; i < nbins*nhires; i++)
    velaxis_hires[i] = i * dvbin_hires;

  
  double pccount = 0.0;
  int ilos = 0;
  for(ilos=0; ilos<nlos; ilos++)
    {
      
      int iconv = 0;
#ifdef OPENMP      
#pragma omp parallel for
#endif
      for(iconv = 0; iconv < nbins*nhires; iconv++)
  	{
	  
	  double rho_wk_H_hires    = resample(ilos, iconv, vmax, rho_wk_H);
	  double rho_wk_H1_hires   = resample(ilos, iconv, vmax, rho_wk_H1);
	  double temp_wk_H1_hires  = resample(ilos, iconv, vmax, temp_wk_H1);
#ifdef HE2LYA
	  double rho_wk_He2_hires  = resample(ilos, iconv, vmax, rho_wk_He2);
	  double temp_wk_He2_hires = resample(ilos, iconv, vmax, temp_wk_He2);
#endif 

	  
#ifdef SELF_SHIELD	  
	  double nHcgs = critH * rho_wk_H_hires / (HMASS * AMU); /* cm^-3 */
	  if (nHcgs >= n0_z)
	    rho_wk_H1_hires = self_shield(nHcgs, log10(temp_wk_H1_hires));
#endif 
	  
	  
#ifdef TEST_KERNEL
	  rho_wk_H1_hires  = 1.0e-5;
	  temp_wk_H1_hires = 1.0e3;
#ifdef HE2LYA
	  rho_wk_He2_hires  = 1.0e-4;
	  temp_wk_He2_hires = 1.0e3;
#endif
	  
#endif // End TEST_KERNEL
	  
	  
#ifdef NO_PECVEL
  	  u_H1[iconv] = velaxis_hires[iconv]; 
#else
	  double vel_wk_H1_hires  = resample(ilos, iconv, vmax, vel_wk_H1);
	  u_H1[iconv] = velaxis_hires[iconv] + vel_wk_H1_hires; /* km s^-1 */
#endif
	  
  	  b_H1             = sqrt(k1_conv * temp_wk_H1_hires); /* cm s^-1 */
	  b2_H1[iconv]     = b_H1 * b_H1 / escale;  /* (km s^-1)^2 */
	  tau_H1_dR[iconv] = k2_conv * rho_wk_H_hires * rho_wk_H1_hires / b_H1;

	  
#ifdef VOIGT
	  aa_H1[iconv] = k_voigt / b_H1;
#endif	  

	  
#ifdef TAU_WEIGHT
	  rho_wk_H_tauw_H1[iconv] = rho_wk_H_hires;
	  temp_wk_H1_tauw[iconv]  = temp_wk_H1_hires;
#endif
	  
	  
#ifdef HE2LYA	  
#ifdef NO_PECVEL
  	  u_He2[iconv] = velaxis_hires[iconv]; 
#else
	  double vel_wk_He2_hires  = resample(ilos, iconv, vmax, vel_wk_He2);
	  u_He2[iconv] = velaxis_hires[iconv] + vel_wk_He2_hires; 
#endif
	  
  	  b_He2             = sqrt(k1_conv_He2 * temp_wk_He2_hires); 
	  b2_He2[iconv]     = b_He2 * b_He2 / escale; 
	  tau_He2_dR[iconv] = k2_conv_He2 * rho_wk_H_nhires * rho_wk_He2_hires / b_He2;
	  
#ifdef VOIGT
	  aa_He2[iconv] = k_voigt_He2 / b_He2;
#endif	  

#ifdef TAU_WEIGHT
	  rho_wk_H_tauw_He2[iconv] = rho_wk_H_hires;
	  temp_wk_He2_tauw[iconv]  = temp_wk_He2_hires;
#endif
#endif // End HE2LYA
	  
	  
#ifdef SILICON
	  /* Here we use the HI weighted temperature and peculiar
	     velocity, but for a hydro simulation that
	     self-consistently includes metals we may rather use the
	     appropriate ion weighted quantity.  To compute relative
	     metallicities (relative to hydrogen by number) we use the
	     standard metallicity definition:
	     
	     [X/Y] = log10(X/Y)_obs - log10(X/Y)_sol 
	     
	     so CIV/H = CIV/C * 10^[C/H] * (C/H)_sol or
	     HeII/H = (HeII/He) * (yhelium/(He/H)_sol) * (He/H)_sol */
	     
#ifdef SILICON_POD
	  /* For [C/H], we use Schaye et al. 2003, ApJ, 596, 768, eq. (8)
	     and to get [Si/H] we use [Si/C]=0.77 in Table 2 of Aguirre
	     et al. 2004, ApJ, 602, 38.  Hence [Si/H] = [Si/C] + [C/H] */  
	  double ZSi_rel = pow(10.0, -2.70 + 0.08 * (ztime - 3.0) + 0.65*(log10(rho_wk_H_hires) - 0.5)); /* 10^[Si/H] */
#else
	  double ZSi_rel = pow(10.0, Z_SI);
#endif
	  double logT  = log10(temp_wk_H1_hires); /* K */
	  double lognH = log10(critH * rho_wk_H_hires / (HMASS * AMU)); /* cm^-3 */
	  
	  get_ionfrac(logT, lognH, &Si2frac, &Si3frac, &Si4frac); 
	  
	  double Si2_Hfrac = Si2frac * ZSi_rel * SI_SOLAR; /* SiII/H */
	  double Si3_Hfrac = Si3frac * ZSi_rel * SI_SOLAR; /* SiIII/H */
	  double Si4_Hfrac = Si4frac * ZSi_rel * SI_SOLAR; /* SiIV/H */
	  
	  b_Si         = sqrt(k1_conv_Si * temp_wk_H1_hires); /* cm s^-1 */
	  b2_Si[iconv] = b_Si * b_Si / escale;  /* (km s^-1)^2 */
	  
	  tau_1190_Si2_dR[iconv] = k2_conv_1190_Si2 * rho_wk_H_hires * Si2_Hfrac / b_Si;
	  tau_1193_Si2_dR[iconv] = (k2_conv_1193_Si2 / k2_conv_1190_Si2) * tau_1190_Si2_dR[iconv];
	  tau_1260_Si2_dR[iconv] = (k2_conv_1260_Si2 / k2_conv_1190_Si2) * tau_1190_Si2_dR[iconv]; 
	  tau_1207_Si3_dR[iconv] = k2_conv_1207_Si3 * rho_wk_H_hires * Si3_Hfrac / b_Si;
	  tau_1394_Si4_dR[iconv] = k2_conv_1394_Si4 * rho_wk_H_hires * Si4_Hfrac / b_Si;
	  tau_1403_Si4_dR[iconv] = (k2_conv_1403_Si4 / k2_conv_1394_Si4) * tau_1394_Si4_dR[iconv];
	  
#endif	  // End SILICON
	}



      
      int ipix = 0;
      for(ipix=0; ipix<nbins*nhires; ipix++)
  	{

	  long long pixel_index = ipix + nbins*nhires*ilos;
	  
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
	  
#ifdef SILICON 
	  double tau_1190_Si2_sum = 0.0;
	  double tau_1193_Si2_sum = 0.0;
	  double tau_1260_Si2_sum = 0.0;
	  double tau_1207_Si3_sum = 0.0;
	  double tau_1394_Si4_sum = 0.0;
	  double tau_1403_Si4_sum = 0.0;
#endif
	  
	  int iconv = 0;




	  /* This could be cleaned up a bit */
#ifdef OPENMP
#ifdef TAU_WEIGHT
	  
#if defined(HE2LYA) && defined(SILICON)
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum, tau_He2_sum, rho_tau_He2_sum, temp_tau_He2_sum, tau_1190_Si2_sum, tau_1193_Si2_sum, tau_1260_Si2_sum, tau_1207_Si3_sum, tau_1394_Si4_sum, tau_1403_Si4_sum)
#endif
	 	  
#if defined(HE2LYA) && !defined(SILICON)
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum, tau_He2_sum, rho_tau_He2_sum, temp_tau_He2_sum)
#endif
	  
#if !defined(HE2LYA) && defined(SILICON)
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum, tau_1190_Si2_sum, tau_1193_Si2_sum, tau_1260_Si2_sum, tau_1207_Si3_sum, tau_1394_Si4_sum, tau_1403_Si4_sum)
#endif

#if !defined(HE2LYA) && !defined(SILICON)	  
#pragma omp parallel for reduction(+:tau_H1_sum, rho_tau_H1_sum, temp_tau_H1_sum)	  
#endif
	  
#else // Else TAUWEIGHT

#if defined(HE2LYA) && defined(SILICON)
#pragma omp parallel for reduction(+:tau_H1_sum, tau_He2_sum, tau_1190_Si2_sum, tau_1193_Si2_sum, tau_1260_Si2_sum, tau_1207_Si3_sum, tau_1394_Si4_sum, tau_1403_Si4_sum)
#endif
	 	  
#if defined(HE2LYA) && !defined(SILICON)
#pragma omp parallel for reduction(+:tau_H1_sum, tau_He2_sum)
#endif
	  
#if !defined(HE2LYA) && defined(SILICON)
#pragma omp parallel for reduction(+:tau_H1_sum, tau_1190_Si2_sum, tau_1193_Si2_sum, tau_1260_Si2_sum, tau_1207_Si3_sum, tau_1394_Si4_sum, tau_1403_Si4_sum)
#endif
	  
#if !defined(HE2LYA) && !defined(SILICON)	  
#pragma omp parallel for reduction(+:tau_H1_sum)	  
#endif	  

#endif // End TAUWEIGHT
#endif // End OPENMP
	  
	  for(iconv=0; iconv<nbins*nhires; iconv++)
  	    {
	      
	      double vdiff_H1 = fabs(velaxis_hires[ipix] - u_H1[iconv]);

	      if (vdiff_H1 > vmax2) vdiff_H1 = vmax - vdiff_H1;
	      
  	      double VH1_0 = vdiff_H1 * vdiff_H1 / b2_H1[iconv]; 
	      
#ifdef QUICK_LINE

	      double VH1_1;

	      if(VH1_0 < 1.0e-4)
		VH1_1 = 1.0 - VH1_0;
	      else
		{
		  double t     = VH1_0 * dx_inv;
		  int tint     = (int)t;
		  double fhi   = t - tint;
		  double flow  = 1 - fhi;
		  
		  VH1_1 = VH1_0 < XMAX ? flow*exp_negx[tint] + fhi*exp_negx[tint+1] : 0.0;
		}
#else
	      double VH1_1 = exp(-VH1_0);
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


	      
	      /* Optionally compute HeII Lyman-alpha absorption */
#ifdef HE2LYA	      
	      double vdiff_He2 = fabs(velaxis_hires[ipix] - u_He2[iconv]);
	      
	      if (vdiff_He2 > vmax2) vdiff_He2 = vmax - vdiff_He2;
	      
  	      double VHe2_0  = vdiff_He2 * vdiff_He2 / b2_He2[iconv]; 
	      
#ifdef QUICK_LINE
	      double VHe2_1;

	      if(VHe2_0 < 1.0e-4)
		VHe2_1 = 1.0 - VHe2_0;
	      else
		{
		  double t_He2    = VHe2_0 * dx_inv;
		  int tint_He2    = (int)t_He2;
		  double fhi_He2  = t_He2 - tint_He2;
		  double flow_He2 = 1 - fhi_He2;
		  
		  VHe2_1   = VHe2_0 < XMAX ? flow_He2*exp_negx[tint_He2] + fhi_He2*exp_negx[tint_He2+1] : 0.0;
		}
#else
	      VHe2_1   = exp(-VHe2_0);
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
	      
	      
	      /* Optionally compute silicon absorption */
#ifdef SILICON 
	      double VSi_0  = vdiff_H1 * vdiff_H1 / b2_Si[iconv]; 
	      
#ifdef QUICK_LINE

	      double VSi_1;
	      
	      if(VSi_0 < 1.0e-4)
		VSi_1 = 1.0 - VSi_0;
	      else
		{
		  double t_Si    = VSi_0 * dx_inv;
		  int tint_Si    = (int)t_Si;
		  double fhi_Si  = t_Si - tint_Si;
		  double flow_Si = 1 - fhi_Si;
	      
		  VSi_1  = VSi_0 < XMAX ? flow_Si*exp_negx[tint_Si] + fhi_Si*exp_negx[tint_Si+1] : 0.0;
		}
#else
	      double VSi_1 = exp(-VSi_0);
#endif
	      
	      tau_1190_Si2_sum += tau_1190_Si2_dR[iconv] * VSi_1;
	      tau_1193_Si2_sum += tau_1193_Si2_dR[iconv] * VSi_1;
	      tau_1260_Si2_sum += tau_1260_Si2_dR[iconv] * VSi_1;
	      tau_1207_Si3_sum += tau_1207_Si3_dR[iconv] * VSi_1;
	      tau_1394_Si4_sum += tau_1394_Si4_dR[iconv] * VSi_1;
	      tau_1403_Si4_sum += tau_1403_Si4_dR[iconv] * VSi_1;
	      
#endif 
	      
	    }

	  if(pixel_index % nhires == 0)
	    {
	      int base_index = floor(pixel_index / nhires);
	      
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

#ifdef SILICON 
	      tau_1190_Si2[base_index] = tau_1190_Si2_sum;
	      tau_1193_Si2[base_index] = tau_1193_Si2_sum;
	      tau_1260_Si2[base_index] = tau_1260_Si2_sum;
	      tau_1207_Si3[base_index] = tau_1207_Si3_sum;
	      tau_1394_Si4[base_index] = tau_1394_Si4_sum;
	      tau_1403_Si4[base_index] = tau_1403_Si4_sum;
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




/** \brief Resamples the line of sight data during the line
 *   convolution if the thermal broadening kernel is judged to be
 *   under-resolved at the native resolution. This is necessary to
 *   ensure that optical depth is conserved, regardless of the mass
 *   resolution of the simulation.  Uses linear interpolation.
 *
 *  \param ilos  The current line of sight number (max is nnumlos-1)
 *  \param iconv The current line of sight pixel number (max is nbins*nhires-1)
 *  \param vmax  The box size in km/s
 *  \param field  The field quantity that will be resampled
 *
 *  \return field_nhires The resampled field quantity.
 *
 */

double resample(int ilos, int iconv, double vmax, double *field)
{
  double field_hires;

  int iv = floor(iconv / nhires);
  int pixel_index = iv + nbins*ilos;
  
  if(nhires > 1)
    {
      double dvbin = vmax / nbins;
      double velhub = iv * dvbin;

      double dvbin_hires = dvbin / nhires;
      double velhub_hires = iconv * dvbin_hires;

      double velhub_next = (velhub_hires < velhub) ? velhub - dvbin : velhub + dvbin;

      int pixel_index_min = nbins*ilos;
      int pixel_index_max = pixel_index_min + nbins - 1;
      
      int pixel_index_next = (velhub_hires < velhub) ? pixel_index - 1 : pixel_index + 1;

      while (pixel_index_next < pixel_index_min)
        pixel_index_next = pixel_index_max;
      while (pixel_index_next > pixel_index_max)
        pixel_index_next = pixel_index_min;

      field_hires = lerp(velhub_hires, velhub, velhub_next, field[pixel_index], field[pixel_index_next]);
    }
  else
    field_hires = field[pixel_index];

  return field_hires;
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
      exit(1);
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
      exit(1);
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

}




/** \brief Write the optical depths to a binary file.
 */

void write_tau()
{
  char fname[400];
  FILE *output;
   
  sprintf(fname, "%s/tauH1_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  
  
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }

  fwrite(tau_H1,sizeof(double),nbins*nlos,output);
#ifdef SELF_SHIELD
  fwrite(rho_wk_H1,sizeof(double),nbins*nlos,output);
#endif
  fclose(output);
  
  
#ifdef TAU_WEIGHT
  sprintf(fname, "%s/tauwH1_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(rho_tau_H1,sizeof(double),nbins*nlos,output);
  fwrite(temp_tau_H1,sizeof(double),nbins*nlos,output);
  fclose(output);
#endif  

  
#ifdef HE2LYA  
  sprintf(fname, "%s/tauHe2r_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_He2,sizeof(double),nbins*nlos,output);
  fclose(output);

  
#ifdef TAU_WEIGHT
  sprintf(fname, "%s/tauwHe2_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(rho_tau_He2,sizeof(double),nbins*nlos,output);
  fwrite(temp_tau_He2,sizeof(double),nbins*nlos,output);
  fclose(output);
#endif
  
#endif

  
#ifdef SILICON 

  sprintf(fname, "%s/tauSi2_1190_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_1190_Si2,sizeof(double),nbins*nlos,output);
  fclose(output);

  sprintf(fname, "%s/tauSi2_1193_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_1193_Si2,sizeof(double),nbins*nlos,output);
  fclose(output);

  sprintf(fname, "%s/tauSi2_1260_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_1260_Si2,sizeof(double),nbins*nlos,output);
  fclose(output);
  
  sprintf(fname, "%s/tauSi3_1207_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_1207_Si3,sizeof(double),nbins*nlos,output);
  fclose(output);

  sprintf(fname, "%s/tauSi4_1394_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_1394_Si4,sizeof(double),nbins*nlos,output);
  fclose(output);
  
  sprintf(fname, "%s/tauSi4_1403_v%d_n%d_z%.3f.dat",path,nbins,nlos,ztime_file);
  if(!(output=fopen(fname,"wb")))
    {
      printf("can't open file `%s`\n\n",fname);
      exit(1);
    }
  fwrite(tau_1403_Si4,sizeof(double),nbins*nlos,output);
  fclose(output);

#endif

}




/** \brief Make a look-up table for the negative exponential in the
    line profile.
 */

#ifdef QUICK_LINE

void LineProfileTable()
{
  
  exp_negx = (double *)calloc(NXTAB+1, sizeof(double));
  
  if(NULL==exp_negx){free(exp_negx); printf("Memory allocation failed.\n"); exit(1);}
  
  double dx =  XMAX / NXTAB;
  
  dx_inv = 1.0 / dx;
  
  int i = 0;
  
  for(i = 0; i <= NXTAB; i++)
    {
      double x = dx * (double)i;
      exp_negx[i] = exp(-x);
    }
}

#endif




/** \brief Allocate LOS memory.  Would be neater to use a structure.
 */

void allocate_los_memory()
{
  
  ilos = (int *)calloc(nlos, sizeof(int));
  if(NULL==ilos){free(ilos); printf("Memory allocation failed.\n"); exit(1);}
  
  xlos = (double *)calloc(nlos, sizeof(double));
  if(NULL==xlos){free(xlos); printf("Memory allocation failed.\n"); exit(1);}

  ylos = (double *)calloc(nlos, sizeof(double));
  if(NULL==ylos){free(ylos); printf("Memory allocation failed.\n"); exit(1);}

  zlos = (double *)calloc(nlos, sizeof(double));
  if(NULL==zlos){free(zlos); printf("Memory allocation failed.\n"); exit(1);}
  
  posaxis = (double *)calloc(nbins, sizeof(double));
  if(NULL==posaxis){free(posaxis); printf("Memory allocation failed.\n"); exit(1);}

  velaxis = (double *)calloc(nbins, sizeof(double));
  if(NULL==velaxis){free(velaxis); printf("Memory allocation failed.\n"); exit(1);}
  
  rho_wk_H = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_H){free(rho_wk_H); printf("Memory allocation failed.\n"); exit(1);}
  
  rho_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_H1){free(rho_wk_H1); printf("Memory allocation failed.\n"); exit(1);}
  
  vel_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==vel_wk_H1){free(vel_wk_H1); printf("Memory allocation failed.\n"); exit(1);}
  
  temp_wk_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_wk_H1){free(temp_wk_H1); printf("Memory allocation failed.\n"); exit(1);}
  
  tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_H1){free(tau_H1); printf("Memory allocation failed.\n"); exit(1);}
  
  
#ifdef TAU_WEIGHT
  rho_tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_tau_H1){free(rho_tau_H1); printf("Memory allocation failed.\n"); exit(1);}
  
  temp_tau_H1 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_tau_H1){free(temp_tau_H1); printf("Memory allocation failed.\n"); exit(1);}
#endif  

  
#ifdef HE2LYA
  rho_wk_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_wk_He2){free(rho_wk_He2); printf("Memory allocation failed.\n"); exit(1);}
  
  vel_wk_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==vel_wk_He2){free(vel_wk_He2); printf("Memory allocation failed.\n"); exit(1);}
  
  temp_wk_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_wk_He2){free(temp_wk_He2); printf("Memory allocation failed.\n"); exit(1);}
  
  tau_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_He2){free(tau_He2); printf("Memory allocation failed.\n"); exit(1);}
  
#ifdef TAU_WEIGHT
  rho_tau_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==rho_tau_He2){free(rho_tau_He2); printf("Memory allocation failed.\n"); exit(1);}
  
  temp_tau_He2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==temp_tau_He2){free(temp_tau_He2); printf("Memory allocation failed.\n"); exit(1);}
#endif
#endif


#ifdef SILICON 
  tau_1190_Si2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_1190_Si2){free(tau_1190_Si2); printf("Memory allocation failed.\n"); exit(1);}
  
  tau_1193_Si2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_1193_Si2){free(tau_1193_Si2); printf("Memory allocation failed.\n"); exit(1);}
  
  tau_1260_Si2 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_1260_Si2){free(tau_1260_Si2); printf("Memory allocation failed.\n"); exit(1);}
  
  tau_1207_Si3 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_1207_Si3){free(tau_1207_Si3); printf("Memory allocation failed.\n"); exit(1);}

  tau_1394_Si4 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_1394_Si4){free(tau_1394_Si4); printf("Memory allocation failed.\n"); exit(1);}
  
  tau_1403_Si4 = (double *)calloc(nlos*nbins, sizeof(double));
  if(NULL==tau_1403_Si4){free(tau_1403_Si4); printf("Memory allocation failed.\n"); exit(1);}
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

#ifdef SILICON 
  free(tau_1190_Si2);
  free(tau_1193_Si2);
  free(tau_1260_Si2);
  free(tau_1207_Si3);
  free(tau_1394_Si4);
  free(tau_1403_Si4);
#endif
  
}


