#########################################################################
#									#
#  LosExtract can be run on a Gadget-3/4  LOS file to extract		#
#  Lyman-alpha forest absorption spectra.  Look at end of file for a	#
#  brief guide to the compile-time options. 				#
#									#
#########################################################################

#--------------------------------------- Compile-time options

#OPTS +=-DOPENMP
#OPTS +=-DTAU_WEIGHT
#OPTS +=-DSELF_SHIELD
#OPTS +=-DNO_PECVEL
#OPTS +=-DVOIGT
#OPTS +=-DQUICK_GAUSS
#OPTS +=-DTEST_KERNEL

#--------------------------------------- Select system
#SYSTYPE="pppzjsb02"
#SYSTYPE="brahan"
SYSTYPE="macbook13"


ifeq ($(SYSTYPE),"pppzjsb02")
CC = gcc	
OPTIMIZE = -O3
CFLAGS   =  $(OPTIMIZE) 
OMPINCL  = -fopenmp
OMPLIB   = -lgomp
endif

ifeq ($(SYSTYPE),"brahan")
CC = gcc	
OPTIMIZE = -O3
CFLAGS   =  $(OPTIMIZE) 
OMPINCL  = -fopenmp
OMPLIB   = -lgomp
endif

ifeq ($(SYSTYPE),"macbook13")
CC = gcc	
OPTIMIZE = -O3
CFLAGS   =  $(OPTIMIZE) 
OMPINCL  = -I/opt/local/include/libomp/ -Xclang -fopenmp
OMPLIB   = -L/opt/local/lib/libomp/ -lomp
endif


CFLAGS += -Wall $(OPTS) $(OMPINCL)
LIBS = -lm $(OMPLIB)

EXEC = LosExtract
OBJS = extract_los.o ion_balance.o
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
	cd ./treecool/; rm -f *~

##############################################################################
#
# Compile-time options:
#
#	- OPENMP	Use OpenMP.  Requires library installation
#
#	- TAUWEIGHT 	Computes the HI Lyman-alpha optical depth weighted
#			density and temperature and saves to an output
#			file.
#
#	- SELF_SHIELD 	Adds a post-processed correction for the
#			self-shielding of neutral hydrogen. Requires
#			selection of a UV background file in
#			parameters.h
#
#	- NO_PECVEL 	Computes the spectra ignoring gas peculiar velocities.
#
#	- VOIGT 	Computes the spectra with a Voigt profile instead of a
#			Gaussian, following Tepper-Garcia (2006,
#			MNRAS, 369, 2025).  Slower than the Gaussian
#			and introduces some shot noise on small scales
#			in the power spectrum due to precision errors.
#			This should occur well below observable scales,
#			but this needs to be verified by the user.
#			Use with care.
#
#	- QUICK_GAUSS 	Uses a look-up table for the Gaussian profile.
#			Faster but slightly less accurate far from the
#			line centre.  Can introduce some shot noise on
#			small scales due to precision errors.
#
#	- TEST_KERNEL 	Test flag.  Runs a model with a constant HI
#			fraction and temperature, with vpec=0.  Used
#			to assess if the LOS pixel scale properly
#			resolves the thermal broadening kernel.
#
##############################################################################
