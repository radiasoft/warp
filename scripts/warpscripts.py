warpscripts_version = "$Id: warpscripts.py,v 1.12 2008/04/03 18:28:19 dave Exp $"

def warpscriptsdoc():
  print """
Contains the command warpscripts() which prints a description of
available scripts
"""

def warpscripts():
  print """
For more info on any scripts, import the script and type the command
scriptnamedoc()


adjustmesh3d.py: provides routines for changing the mesh size
appendablearray.py: declares an array class which can be appended to
b_fields.py: Functions to generate and plot 'bgrd' arrays for solenoids
cir_match.py: routines for matching beam using the circe module
cirplots.py: routines for plotting circe ouput
colorbar.py: Obsolete
ctl.py: step and generate commands (automatically imported)
drawlattice.py: routine for drawing the lattice (automatically imported)
drifts.py: routine to create drift elements in blank spots in the lattice
egun_like.py: routine for steady-state calculations (like diode problems)
eguntowarp.py: inject particles into a simulation, based on egun results
env_match.py: routines for matching beam using the env module
env_tools.py: matchint tools for env code - still in development
envtuner.py: experimental script for mouse driven tuning of the beam envelope
errorcheck.py: makes numerous consistency and error checks
extpart.py: sets up extrapolated particle windows
find_mgparam.py: routine for finding optimal multigrid relaxation parameter
find_sorparam.py: routine for finding optimal sor relaxation parameter
Fitting.py: a general function to do a least squares fit
fixwxy.py: routine to shift particles to satisfy average quantities (wxy only)
fringedquads.py: routine which takes a hard-edged elements and generates
                 elements with fringe fields
generateconductors.py: classes for generating a broad class of geometric
                       objects for boundary conditions for the fieldsolver
generatelattice.py: Obsolete
getzmom.py: python function for calculating particle moments
gistdummy.py: stub routines when gist is not available
grid_1d.py: routine to bin one-dimensional data onto a grid
grid_2d.py: routine to bin two-dimensional data onto a grid
hermestools.py: provides tools for using Hermes
histplot.py: routine which makes a standard set of history plots
             (uses histplots)
histplots.py: routines to make various history plots (automatically imported)
lattice.py: script for MAD-like lattice definition
latticegenerator.py: generates lattice based on design code parameters
loadbalance.py: provides load balancing operations for the parallel version
matchenv.py: envelope matching routine
monitor.py: allows remote monitoring of a run
mplot.py: routines for mountain range style plots
mphoto.py: generates greyscale "photos" of beam density using 'rho' arra
noparens.py: allows functions to be called without the parenthesis '()'
optimizer.py: implementation various minimization algorithms
orbitrack.py: suite for self-consistent particle tracking and multispecies work
ParaKV.py: functions to impose a parabolic temperature distribution
parallel.py: functions for the parallel version
particlescraper.py: class to handle particle scraping on geometric objects
particles.py: function for dealing with particles (automatically imported)
PepperGrid.py: provides several masks for pepper-pot type diagnostics
photo_processing.py: combines tif photos into a montage or a gif animation
plot_conductor.py: routines for ploting internal conductors (automatically
                   imported)
plane_restore.py: restores data that was save at a plane (using plane_save)
plane_save.py: saves data as beam passes through a plane so it can be restored
               in subsequent simulations (using plane_restore)
plarr3d.py: routines to plot 3-D array of points with or without connecting
            lines
postrack.py: suite for postprocessing results where 'orbitrack' is used
printarraysizes.py: for debugging, prints sizes of all allocated arrays
printparameters.py: prints basic run parameters (automatically imported)
printparameters3d.py: prints basic run parameters (automatically imported)
printparametersrz.py: prints basic run parameters (automatically imported)
PRpyt.py: interface to read hdf files using PyTables
ptob.py: routine to save multi-dimensional array in a basis readable pdb file
PWpyt.py: interface to write hdf files using PyTables
Opyndx.py: defines convenience routines for using the OpenDX wrapper
pzplots.py: plots of z moments of particles (automatically imported)
rami_match.py: Slice code matching routine w/starting point in between quads
rami_scripts.py: convenience functions and alternate history plotting routines
realboundaries.py: routines allowing automatic boundaries for wxy
residual.py: calculates residual of Poisson's equation, del**2 phi - rho
runcounter.py: implements a counter for series of simulations
setupvalidation.py: convenience function for creating code validation decks
singleparticle.py: sets WARP up for single particle calculations
slitscanner.py: emulates slit scanners, both parallel and crosses slits
sortwxy.py: sorts particles based on z location
subcycle.py: functions for setting up subcycling of particles
timedependentvoltage.py: functions for setting up time dependent voltages
warp.py: fundamental warp script (automatically imported)
warpfortran.py: documentation of some fortran routines available
warphelp.py: brief list of some available commands
warp_objects.py: generate surface of revolution conductors from simple
                 description
warpparallel.py: functions for parallel version (automatically imported)
warpplots.py: various convenient plotting routines such as particle plots
              (automatically imported)
warpscripts.py: lists all available scripts for warp
wxy_match.py: routines for matching beam using the wxy module

============ UMD SCripts ==================
ParaKV.py: functions to impose a parabolic temperature distribution
mphoto.py: generates greyscale "photos" of beam density using 'rho' arra
photo_processing.py: combines tif photos into a montage or a gif animation
Fitting.py: a general function to do a least squares fit
orbitrack.py: suite for self-consistent particle tracking and multispecies work
postrack.py: suite for postprocessing results where 'orbitrack' is used
rami_scripts.py: convenience functions and alternate history plotting routines
rami_match.py: Slice code matching routine w/starting point in between quads
b_fields.py: Functions to generate and plot 'bgrd' arrays for solenoids

=========== Other files =================
work.gs: style file describing labels and axis for gist plots
warpstyle.gs: same as work.gs

Gist color palettes:
cool.gp: black box radiation colors in reverse order
earth.gp: elevation colors from map, blue to green to white (default)
gray.gp: gray scale
grey.gp: grey scale
heat.gp: black box radiation colors
ncar.gp: color table from NCAR package, chosen to enhance contrast
rainbowaf.gp: rainbow colors with white at the bottom
rainbow.gp: rainbow colors
rainbowinv.gp: rainbow colors in reverse order
stern.gp: color table from IDL (c) Research Systems, Inc.
yarg.gp: gray scale in reverse order


############################################################################
Basic scripts:
  warp
  warputils

Particles:
  extpart
  fixwxy
  ParaKV
  particlescraper
  particles
  postrack
  orbitrack
  singleparticle
  sortwxy

Envelope and matching:
  cir_match
  cirplots
  env_match
  env_tools
  envtuner
  hermestools
  matchenv
  rami_match
  wxy_match

Conductors and Field solver:
  adjustmesh3d
  find_mgparam
  find_sorparam
  generateconductors
  plot_conductor
  realboundaries
  timedependentvoltage
  warp_objects

Lattice:
  b_fields
  drawlattice
  drifts
  fringedquads
  generatelattice
  lattice
  latticegenerator

Run Modes:
  ctl
  egun_like
  eguntowarp
  plane_restore
  plane_save
  subcycle

Parallel:
  loadbalance
  parallel
  warpparallel

Plotting and Diagnostics:
  colorbar
  errorcheck
  getzmom
  gistdummy
  histplot
  histplots
  mplot
  mphoto
  PepperGrid
  photo_processing
  plot_conductor
  plarr3d
  printarraysizes
  printparameters
  printparameters3d
  printparametersrz
  Opyndx
  pzplots
  rami_scripts
  residual
  slitscanner
  warpplots

Documentation:
  warpfortran
  warphelp
  warpscripts

General tools:
  appendablearray
  Fitting
  grid_1d
  grid_2d
  monitor
  noparens
  optimizer
  PRpyt
  ptob
  PWpyt
  runcounter
  setupvalidation

"""
