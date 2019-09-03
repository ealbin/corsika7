PRESHOWER 2.0; 13.03.2012
--------------------------------
1. Authors: P. Homola(1), R. Engel (2), A. Pysz (3) and H. Wilczynski (1)
   (1): Institute of Nuclear Physics, Krakow, Poland
   (2): Forshungszentrum Karlsruhe, Institut fuer Kernphysik, Germany	      
   (3): AGH University of Science and Technology, Krakow, Poland
   
2. Corresponding author: Piotr Homola (Piotr.Homola@ifj.edu.pl)

3. Special credits 

- N. A. Tsyganenko (author of the subroutine computing the geomagnetic field)
- C. Bleve (author of the update of the geomagnetic field model)

4. Publications

- P. Homola et al., Comp. Phys. Comm. 173 (2005) 71
- P. Homola et al., Comp. Phys. Comm. 184 (2013) 1468

5. DESCRIPTION of the program

Propagation of ultra high energy photon before its entering the Earth's
atmosphere is simulated, taking into account gamma conversion into 
e+/- pair and subsequent bremsstrahlung of the electrons in the geomagnetic
field. As a result we obtain a bunch (a preshower) of particles, mainly 
photons and a few electrons, instead of the primary gamma. The information
about all the particles in the preshower is returned to CORSIKA or saved 
(if in stand-alone mode).   

6. Summary of changes comparing to PRESHOWER 1.0

- A veto algorithm was introduced in the gamma conversion and bremsstrahlung 
  tracking procedures. The length of tracking step is now variable along 
  the track and depends on the probability of the process expected to occur. 
  The new algorithm reduces significantly the number of tracking steps and 
  speeds up the execution of the program by 5-10 times. 

- The geomagnetic field model has been updated from IGRF-8 to IGRF-11, 
  allowing for interpolations up to the year 2015. 

- One minor bug has been found and fixed in the auxiliary function kappa(x).
  used for calculation of bremsstrahlung probability. The interpolation performed 
  in this function failed for the rare case of x=10.0. This happened because 
  of a faulty definition of the last interval where the interpolation was done. 
  As a result of this bug the input value x=10.0 was excluded from the computations. 
  
- the size of the array part_out, which stores the output particle data, has been 
  increased to 100000 entries.
   
- Numerical Recipes procedure to calculate modified Bessel functions have been
  replaced with an open source CERN routine DBSKA

7. Content of the packet

-IGRF-11.f: external procedure for calculation of the geomagnetic field acc. to 
         the IGRF-11 model (author: N. Tsyganenko, update of coefficients by 
	 C. Bleve), it is called from within PRESHOWER 2.0 	   
-Makefile: compilation instructions (type 'make') to compile PRESHOWER 2.0
           in a strand-alone mode
-preshw.c: source of the PRESHOWER 2.0 procedure
-prog.c: program that calls PRESHOWER 2.0 procedure
-README.txt: this file
-rmmard.f: dummy routine that prevents core dump while linking the stand_alone
         version. The real rmmard (a part of CORSIKA package) is used in 
	 the PRESHOWER-CORSIKA mode.
-utils.c: auxiliary functions 
-utils.h: a header file for utils.c
-veto.c: veto subroutines for gamma conversion and bremsstrahlung of electrons + related
         auxiliary functions
-veto.h: a header file for veto.c

examples of test input and output files:
-INPUT: an example of input file
-multirun.dat: typical summary output file, obtained with INPUT
-part_out.dat: typical detailed particle data file, obtained with INPUT
-out.txt: typical messages printed to the standard output

8. Compilation and linking on Linux system (stand-alone version)

- create a file nr_fun.c in the PRESHOWER 2.0 home directory and place
  the header:
  
#include <stdio.h>
#include <math.h>
  
  at the beginning of this file;
    
- below the header, the source code of the Numerical Recipes procedures
  ran2 and nrerror should be pasted;
	
- in the PRESHOWER 2.0 home directory type 'make' for compilation and
  linking the program

- troublesooting:
If "make" fails with a message like "make: f77: Command not found" you most likely
your default Fortran compiler is most likely gfortran. Adding a line:

FC= gfortran

in Makefile should fix the problem.

9. Running PRESHOWER 2.0

Prepare the input file and the type './preshower' to run.

10. Using PRESHOWER 2.0 within other air shower simulation codes

While linking PRESHOWER 2.0 to air shower simulation codes other than CORSIKA
one has to take care about proper format of the input parameters and adequate
treatment of the PRESHOWER 2.0 output by the external code. In case of 
difficulties do not hesitate to e-mail Piotr.Homola@ifj.edu.pl
