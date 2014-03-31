This directory contains a set of functions for performing a time 
series hydrodynamic simulation of a heaving buoy. Also provided 
are multiple sets of buoy data files, generated using WAMIT for 
a variety of buoy shapes and sizes, providing hydrodynamic 
coefficients and excitation forces for each buoy for a range
of frequencies. 

The main functions of interest for users are the following:


buoysimsetup - a function for setting up a 


seasetup.m -  a function for generating appropriate sets of sinusoids 
              for the simulation of sea waves by the function buoyodesim.m. 
              Currently capable of generating single sinusoids, multiple 
              sinusoids or a Pierson-Moscovitch sea spectrum.


buoysetup.m - a function for locating the simulation files for a 
              particular buoy from the buoy library, and loading 
              the physical parameters of the buoy from the appropriate
              data file. The function buoysimsetup uses the file 
              locations returned by buoysetup.m to set up the buoy 
              simulation.
