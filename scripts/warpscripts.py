warpscripts_version = "$Id: warpscripts.py,v 1.3 2001/01/11 22:51:46 dave Exp $"

def warpscriptsdoc():
  print """
Contains the command warpscripts() which prints a description of
available scripts
"""

def warpscripts():
  print """
For more info on any scripts, import the script and type the command
scriptnamedoc()

appendablearray.py: declares an array class which can be appended to
cir_match.py: routines for matching beam using the circe module
cirplots.py: routines for plotting circe ouput
ctl.py: step and generate commands (automatically imported)
drifts.py: routine to create drift elements in blank spots in the lattice
egun_like.py: routine for steady-state calculations (like diode problems)
env_match.py: routines for matching beam using the env module
envtuner.py: experimental script for mouse driven tuning of the beam envelope
errorcheck.py: makes numerous consistency and error checks
expt_diagnostic.py: emulates experimental slit scanner plots of phase space
find_mgparam.py: routine for finding optimal multigrid relaxation parameter
find_sorparam.py: routine for finding optimal sor relaxation parameter
fixwxy.py: routine to shift particles to satisfy average quantities (wxy only)
fringedquads.py: routine which takes a hard-edged elements and generates
                 elements with fringe fields
grid_1d.py: routine to bin one-dimensional data onto a grid
grid_2d.py: routine to bin two-dimensional data onto a grid
histplot.py: routine which makes a standard set of history plots
             (uses histplots)
histplots.py: routines to make various history plots
lattice.py: script for MAD-like lattice definition
matchenv.py: envelope matching routine
monitor.py: allows remote monitoring of a run
mplot.py: routines for mountain range style plots
noparens.py: allows functions to be called without the parenthesis '()'
optimizer.py: implementation various minimization algorithms
plot_conductor.py: routines for ploting internal conductors
ptob.py: routine to save multi-dimensional array in a basis readable pdb file
pyBasis.py: defines routines to emulate basis functionality
            (automatically imported)
realboundaries.py: routines allowing automatic boundaries for wxy
residual.py: calculates residual of Poisson's equation, del**2 phi - rho
runcounter.py: implements a counter for series of simulations
singleparticle.py: sets WARP up for single particle calculations
warp.py: fundamental warp script
warphelp.py: brief list of some available commands
warpplots.py: various convenient plotting routines such as particle plots
              (automatically imported)
warpstyle.gs: style file for gist graphics
warpscripts.py: routine to descrive available scripts
wxy_match.py: routines for matching beam using the wxy module
"""

