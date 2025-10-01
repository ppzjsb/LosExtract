
/********************************************************************************/

/**! \file proto.h                                                          
 *                                                
 *  \brief header for defining function prototypes.
 */

/********************************************************************************/

/* extract_los.c */
void read_los(void);
double resample(int ilos, int iconv, double vmax, double *field);
void write_tau(void);
void compute_absorption(void);
void allocate_los_memory(void);
void free_los_memory(void);
#ifdef QUICK_LINE
void LineProfileTable(void);
#endif

/* ion_balance.c */
#ifdef SELF_SHIELD
double self_shield(double nHcgs, double logT);
double ion_balance(double nHcgs, double logT);
void InitCool(void);
void ReadIonizeParams(char *fname);
void IonizeParamsTable(void);
void MakeCoolingTable(void);
void InitCoolMemory(void);
void SelfShieldFit(void);
#endif

/* utils.c */
double dmax(double x, double y);
double dmin(double x, double y);
double lerp(double x, double x0, double x1, double y0, double y1);

/* cloudy_tables.c */
#ifdef SILICON
void InitCloudy(void);
void get_ionfrac(double logT, double lognH, double *Si2frac, double *Si3frac, double *Si4frac);
void read_iontable(char *fname);
void InitIonTableMemory(void);
#endif
