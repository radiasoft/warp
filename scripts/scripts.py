scripts_version = "$Id: scripts.py,v 1.1.1.1 2000/10/16 18:34:19 dave Exp $"
print """
These are some of the scripts available for use:
(Note that some are not yet available though in python)
------------------------------------------------------------------------------
appendablearray: Creates a subclass of arrays that can be appended to.

------------------------------------------------------------------------------
ctl: Provides the control commands generate and step.

------------------------------------------------------------------------------
egun_like: Sets the code up to do a steady-state calculation.

------------------------------------------------------------------------------
emittance: Calculates the emittance after removing particles in regions with
	   a density below a cutoff value. (incomplete)

------------------------------------------------------------------------------
env_match: Provides various matching routines using the envelope code.

------------------------------------------------------------------------------
envtuner: Provides primitive mouse driven tuning of the beam envelope.

------------------------------------------------------------------------------
drifts: Creates drift lattice elements.

------------------------------------------------------------------------------
find_sorparam: Uses a simple search algorithm to find a good value for the
               over-relaxation parameter for the SOR field solver.
find_mgparam: Uses a simple search algorithm to find a good value for the
              over-relaxation parameter for the multigrid field solver.

------------------------------------------------------------------------------
fringedquads: Converts hard-edaged quads into fringed quads.

------------------------------------------------------------------------------
fixwxy: Adjusts particles transverse positions and velocities to exactly match
        the input beam parameters (should only be used immediately after a
	wxy generate and before any time steps).

------------------------------------------------------------------------------
getzmom: Contains the calls needed to perform the z moments calculation
         from the basis interpreter.

------------------------------------------------------------------------------
grid_1d: Loads 1-D set of data onto a 1-D grid.
grid_2d: Loads 1-D set of data onto a 2-D grid.

------------------------------------------------------------------------------
hibeam:
hibeamlattice:
hibeamdefaults:
  A series of scripts which attempt to use hibeam input files to run WARP.
  (incomplete)

------------------------------------------------------------------------------
lattice: Provides a class structure for describing an accelerator lattice.
         Contains subroutines to convert the input from a MAD lattice into the
	 WARP lattice.

------------------------------------------------------------------------------
lbl_diagnostic: Makes plots of x and y phase space like the diagnostic plots
                produced on LBL experiments.

------------------------------------------------------------------------------
namelist: Reads fortrans namelist files (incomplete).

------------------------------------------------------------------------------
plane_restore: Saves particle and phi data at a plane (for restarts).
plane_save: Restores the particle and phi data saved by plane_save.
plane_reduce: Specialized script to reduce the size of the grid data saved
              using the plane_save routine.
plane_filter: Specialized script to filter particle data to model transport
              through a grid at the plane saved (via the plane_save script)
              and to reduce the size of the grid data saved.

------------------------------------------------------------------------------
plot_conductor: Functions to make plots of conductors and fields in each plane.

------------------------------------------------------------------------------
plot_extrapolated_particles: Functions to plot particle data extrapolated in
                             the moments calculation to the z plane specified.

------------------------------------------------------------------------------
plot_lab_moments: Functions to make plots of the lab window moments.

------------------------------------------------------------------------------
pyBasis: Imports pybasis module and replicates many functions of basis.
	 Note that this is automatically imported with warp.

------------------------------------------------------------------------------
realboundaries: Implements the capacity matrix solver in the slice code. Can
                create quadrupole rods and round pipes based on the apertures
		of elements.

------------------------------------------------------------------------------
residual: Calculates eps0*(del**2 phi) - rho.

------------------------------------------------------------------------------
rfalsi1: Regula falsi root finder in 1 dimension.
rfalsi2: Regula falsi root finder in 2 dimension.
rfalsi3: Regula falsi root finder in 3 dimension.

------------------------------------------------------------------------------
scripts: This file.

------------------------------------------------------------------------------
sig0_matrix: Calculates sigma0 by tracking many particles over one lattice 
             period, calculating the sigma0 for each pair of particles, and
             then averaging over all of the values of sigma0.

------------------------------------------------------------------------------
singleparticle: Sets WARP up for single particle calculations.

------------------------------------------------------------------------------
warp: The primary WARP module which imports the various physics packages.

------------------------------------------------------------------------------
warpplots: Contains various convenient plotting routines such as particle
           plots.
	   Note that this is automatically imported with warp.

------------------------------------------------------------------------------
wxy_match: Matching routines which work with the slice code.

------------------------------------------------------------------------------
"""
