
#########################################################################
#									#
#  LosExtract can be run on a Gadget-3/4 LOS file to extract		#
#  Lyman-alpha forest absorption spectra.  Look at end of file for a	#
#  brief guide to the compile-time options. 				#
#									#
#########################################################################

#--------------------------------------- Compile-time options

OPTS +=-DOPENMP
#OPTS +=-DTAU_WEIGHT
#OPTS +=-DHE2LYA
OPTS +=-DSILICON
#OPTS +=-DSILICON_POD
#OPTS +=-DSELF_SHIELD
#OPTS +=-DNO_PECVEL
OPTS +=-DVOIGT
OPTS +=-DRESAMPLE
OPTS +=-DQUICK_LINE
#OPTS +=-DTEST_KERNEL

#--------------------------------------- Select system
SYSTYPE="pppzjsb02"
#SYSTYPE="brahan"
#SYSTYPE="macbook13"


ifeq ($(SYSTYPE),"pppzjsb02")
CC = gcc	
OPTIMIZE = -O3
CFLAGS   =  $(OPTIMIZE) -Wall 
OMPINCL  = -fopenmp
OMPLIB   = -lgomp
endif

ifeq ($(SYSTYPE),"brahan")
CC = gcc	
OPTIMIZE = -O3
CFLAGS   =  $(OPTIMIZE) -Wall -fcommon -Wno-unused-result
OMPINCL  = -fopenmp
OMPLIB   = -lgomp
endif

ifeq ($(SYSTYPE),"macbook13")
CC = gcc	
OPTIMIZE = -O3
CFLAGS   =  $(OPTIMIZE) -Wall
OMPINCL  = -I/opt/local/include/libomp/ -Xclang -fopenmp
OMPLIB   = -L/opt/local/lib/libomp/ -lomp
endif


CFLAGS += $(OPTS) $(OMPINCL)
LIBS = -lm $(OMPLIB)

EXEC = LosExtract
OBJS = extract_los.o ion_balance.o utils.o cloudy_tables.o
INCL = proto.h global_vars.h


$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS): $(INCL)

.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: tidy
tidy:
	rm -f $(OBJS) $(EXEC) *~
	cd ./idl_powerspec/; rm -f *~
	cd ./idl_readdata/; rm -f *~
	cd ./idl_utils/; rm -f *~
	cd ./treecool/; rm -f *~ 
	cd ./cloudy_tables/; rm -f *~

##############################################################################
#
# Compile-time options:
#
#	- OPENMP	Use OpenMP.  Requires library installation
#
#	- TAUWEIGHT 	Computes the Lyman-alpha optical depth weighted
#			density and temperature and saves to an output
#			file.
#
#	- HE2LYA	Optionally computes the HeII Lyman-alpha optical depths. 
#
#	- SILICON 	Optionally computes SiII (1190, 1193, 1260), SiIII
#			(1207) and SiIV (1394, 1403) optical depths.
#			Always uses a Gaussian profile.  The density
#			independent [Si/H] is set in parameters.h
#
#	- SILICON_POD 	Density dependent pixel optical depth value for
#			[Si/H], using [C/H] from Schaye et al. 2003,
#			ApJ, 596, 768, eq. (8) and [Si/C]=0.77 from
#			Table 2 of Aguirre et al. 2004, ApJ, 602, 38.
#			Overides the default [Si/H] model.  Note the
#			POD metallicities are derived assuming a
#			Haardt & Madau UVB model.
#
#	- SELF_SHIELD 	Adds a post-processed correction for the
#			self-shielding of neutral hydrogen.  Does not
#			make much difference for quick-Lya runs due to
#			missing dense gas.
#
#	- NO_PECVEL 	Computes the spectra ignoring gas peculiar velocities.
#
#	- VOIGT 	Computes the Lyman-alpha spectra with a Voigt profile
#			instead of a Gaussian, following Tepper-Garcia
#			(2006, MNRAS, 369, 2025).  Slower than the
#			Gaussian and can introduce some noise on small
#			scales in the power spectrum.  This should be
#			well below observable scales, but this needs
#			to be verified by the user.
#
#	- RESAMPLE 	Performs a check on the sampling rate of the
#			thermal broadening kernel.  If the native
#			resolution of the input files are
#			undersampling the kernel, this performs a
#			resample using linear interpolation. It is
#			*strongly recommended* this option is always
#			left on to ensure optical depth convergence.
#
#	- QUICK_LINE 	Uses a combination of a Taylor expansion and a
#			look-up table for the negative exponential
#			term in the line profile.  This is faster than
#			the default and just as accurate for the cost
#			of an extra 7.6 MiB memory.
#
#	- TEST_KERNEL 	Test option.  Can be used to assess if the LOS
#			pixel scale properly resolves the thermal
#			broadening kernel.
#
##############################################################################
