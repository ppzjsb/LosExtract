
/********************************************************************************/

/**! \file proto.h                                                          
 *                                                
 *  \brief header for defining function prototypes.
 */

/********************************************************************************/

/* extract_los.c */
void read_los();
void resample_los(double vmax);
void write_tau();
void compute_absorption();
void allocate_los_memory();
void allocate_los_hires_memory();
void free_los_memory();
void free_los_hires_memory();
#ifdef QUICK_GAUSS
void GaussianProfileTable();
#endif

/* ion_balance.c */
#ifdef SELF_SHIELD
double self_shield(double nHcgs, double logT);
double ion_balance(double nHcgs, double logT);
void InitCool();
void ReadIonizeParams(char *fname);
void IonizeParamsTable();
void MakeCoolingTable();
void InitCoolMemory();
void SelfShieldFit();
#endif

/* utils.c */
double dmax(double x, double y);
double dmin(double x, double y);
double lerp(double x, double x0, double x1, double y0, double y1);
