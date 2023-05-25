---------------------------------------------------------------

LosExtract

James S. Bolton

---------------------------------------------------------------

LosExtract can be run on a Gadget-3/4 Sherwood LOS file to extract
Lyman-alpha optical depths. 

RUNNING THE CODE:

Code requirements are a C compiler with (optionally) OpenMP. See e.g.

https://www.openmp.org/resources/openmp-compilers-tools/

To compile the code, at the command line type:
make

To run the code, at the command line type:
./LosExtract &lt path &rt

where:
&lt path &rt     (path to LOS files)

After running a binary output file will appear in the directory
containing the LOS files.

Parameters that can be varied in the code are contained in the file
parameters.h.

Compile time options are included in the Makefile, along with descriptions.

If changing the compile time options one should generate a new binary.
To do this, at the command line type:
make clean
make


ADDITIONAL FILES:

treecool/ contains a selection of ASCII tables with the (spatially
uniform) photoionisation and photoheating rate tables used within
Gadget-4.  These can be selected within parameters.h.  Currently these
are mainly variations on the Puchwein et al. 2019 model used for the
Sherwood-Relics project.  See e.g.

Puchwein et al. 2019, MNRAS, 485, 47.
Puchwein et al. 2022, MNRAS, 519, 6162.

Some other models from these papers are also included in case useful

Haardt &amp Madau, 2012, ApJ, 746, 125.
Khaire &amp Srianand, 2019, MNRAS, 484, 4174.
Faucher-Giguere 2020, MNRAS, 493, 1614.

idl_readdata/ contains a simple IDL code that reads in the output from
SpecExtract.

idl_powerspec/ contains an IDL code that computes the 1D Lyman-alpha
transmission power spectrum.

idl_utils/ contains some useful IDL code for common tasks (e.g. making the
outputs.txt list for a Gadget-3/4 run)

---------------------------------------------------------------

