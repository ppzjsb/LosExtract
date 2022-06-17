
/* Function prototypes */

/* extract_los.c */
void read_los();
void write_tau();
void compute_absorption();
void allocate_los_memory();
void free_los_memory();
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
