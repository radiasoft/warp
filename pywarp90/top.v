top
#@(#) File TOP.V, version $Revision: 3.103 $, $Date: 2003/08/25 23:04:29 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package TOP of code WARP
# TOP - all variables and code which are needed by more than one package.
#       --> These common blocks are available to all packages <--
# Alex Friedman,  LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194

{
# --- Macro for setting nohazard depending on compiler ---
%define([NOHAZARD],[Immediate([c])])
%ifelse(COMPILER,CFT77,[define([NOHAZARD],[Immediate([cdir$$ ivdep])])])
# --- Macro for setting novector depending on compiler ---
%define([NOVECTOR],[Immediate([c])])
%ifelse(COMPILER,CFT77,[define([NOVECTOR],[Immediate([cdir$$ novector])])])
%ifelse(COMPILER,CIVIC,[define([NOVECTOR],[Immediate([cdir$$ nextscalar])])])
# --- Macro for recursive command, needed for F90 but not understood in F77
%ifelse(COMPILER,F77|PGF77,[define([RECURSIVE],[])])
# --- Define large and small for the machine ---
MACHEPS = 1.0e-14
#%ifelse(WORDSIZE,64,[define([LARGEPOS],[1.0e+99])],
#%                   [define([LARGEPOS],[1.0e+36])])
#%ifelse(WORDSIZE,64,[define([SMALLPOS],[1.0e-99])],
#%                   [define([SMALLPOS],[1.0e-36])])
LARGEPOS = 1.0e+36
SMALLPOS = 1.0e-36
# --- Protect divides ---
%define dvnz sign(abs($1)+SMALLPOS,$1) 
NCONTROL = 50              # Length of diagnostic control arrays
NELEMENT = 100             # Default length of arrays describing lattice
NPARPGRP = 256             # Number of particle per group
NSUBSETS = 3               # Max number of ptcl "subsets" for scatter plots
NWINDOWS = 9               # Max number of diagnostic windows
NUMZMMNT = 28              # Number of z moments
TNWINM  = 2*NWINDOWS               # Used only for data statements
NWPNSP1 = NWINDOWS + NSUBSETS + 1  # Used only for data statements
NEVER   = 0
SELDOM  = 1
ALWAYS  = 2
# Parallel message tags
M_PARTICLES  = 100
M_GET_PART   = M_PARTICLES + 1
M_GET_ZMMNTS = M_GET_PART + 1
M_GET_HIST   = M_GET_ZMMNTS + 1
M_MAX        = M_GET_HIST + 1
M_PHI        = M_MAX + 1
M_RHO        = M_PHI + 1
M_TRANSPOSE  = M_RHO + 1
M_SETEEARS   = M_TRANSPOSE + 1
M_GET_PHI    = M_SETEEARS + 1
M_GET_RHO    = M_GET_PHI + 1
M_SUM        = M_GET_RHO + 1
}

*********** Code_version:
# -------- CVS updates this when a new major "release" occurs ------------
codeid   character*8  /"warp r2"/     # Name of code, and major version

*********** TOPversion:
# Version control for global commons
verstop character*19 /"$Revision: 3.103 $"/ # Global common version, set by CVS

*********** Machine_param:
wordsize integer /64/ # Wordsize on current machine--used in bas.wrp
largepos real    /LARGEPOS/ # Large positive number
smallpos real    /SMALLPOS/ # Small positive number

*********** GlobalVars:
nparpgrp  integer /NPARPGRP/ # Number of particles per group. Effects size
                             # of temporaries and cache use.
dirichlet integer /0/        # constant for boundary condition (constant potential)
neumann   integer /1/        # constant for boundary condition (derivative = 0)
periodic  integer /2/        # constant for boundary condition 
absorb    integer /0/        # constant for particle absorption at boundaries
reflect   integer /1/        # constant for particle reflection at boundaries

*********** Timers dump parallel:
starttime   real /0./ # CPU start time (in seconds)
gentime     real /0./ # CPU for generate (in seconds)
steptime    real /0./ # Total CPU run time minus gentime (in seconds)
plottime    real /0./ # Time making automatic plots (in seconds)
momentstime real /0./ # Time to calculate the moments (in seconds)
fstime      real /0./ # Time to do the field solve, including beforefs and
                      # afterfs (in seconds)
latticetime real /0./ # Time to apply the fields from the lattice
dumptime    real /0./ # Time to do the data dumps, using the dump command
                      # (in seconds)
temperaturestime real /0./ # Time to calculate the temperatures (in seconds)

*********** Beam_acc dump:
# Beam and Accelerator variables needed by more than one package
aion          real /0./               [1] # A of ion (Carbon)
a0            real /0./               [m] # Initial beam width in x
ap0           real /0./               [1] # Initial beam envelope vx/vz
b0            real /0./               [m] # Initial beam width in y
bp0           real /0./               [1] # Initial beam envelope vy/vz
x0            real /0./               [m] # Initial beam centroid in x
xp0           real /0./               [1] # Initial beam centroid vx/vz
y0            real /0./               [m] # Initial beam centroid in y
yp0           real /0./               [1] # Initial beam centroid vy/vz
tunelen       real /0./               [m] # Length for tune calc (lattice per.)
dedr          real /0./             [E/m] # Uniform focusing 
                                          #   radial Electric field gradient
dexdx         real /0./             [E/m] # Uniform focusing X-Electric field 
                                          #   gradient (dE_x/dx) 
deydy         real /0./             [E/m] # Unifrom focusing Y-Electric field 
                                          #   gradient (dE_y/dy) 
dbdr          real /0./             [T/m] # Uniform focusing B-field gradient
ekin          real /0./              [eV] # Input beam kinetic energy
emit          real /0./           [m-rad] # Perp Emittance of beam (rms-edge)
emitx         real /0./           [m-rad] #    X-Emittance of beam (rms-edge)
                                          #    auto-set to emit, if zero 
emity         real /0./           [m-rad] #    Y-Emittance of beam (rms-edge)
                                          #    auto-set to emit, if zero 
emitn         real /0./           [m-rad] # Normalized   emittance of beam 
emitnx        real /0./           [m-rad] # Normalized X-emittance of beam 
                                          #    auto-set to emitn, if zero 
emitny        real /0./           [m-rad] # Normalized Y-emittance of beam
                                          #    auto-set to emitn, if zero  
ibeam         real /0./               [A] # Input beam current
zion          real /0./               [1] # Z of ion 
straight      real /0./               [1] # Percent of beam that isn't cigar
emitlong      real /0./           [m-rad] # Longitudinal emittance
eears         real /1./             [E/m] # Cigar Eears multiplier / switch
gfactor       real /0./               [1] # Geometric factor (is set if 0)
rwall         real                    [m] # Effective wall radius
gammabar      real                    [1] # Relativistic gamma factor
vbeam         real /0./             [m/s] # Beam speed: use 0 if ekin sets it
lrelativ      logical /.false./           # Flag for relativity

*********** Constant:
# Physical constants, in one place so they stay consistent
pi            real /3.14159265358979323/  # Pi
amu           real /1.66053873e-27/  [kg] # Atomic Mass Unit
clight        real /2.99792458e+8/  [m/s] # Speed of light in vacuum (exact)
echarge       real /1.602176462e-19/  [C] # Proton charge
emass         real /9.10938188e-31/  [kg] # Electron mass
eps0          real /0.0/     [F/m]        # Permittivity of free space
                                          # set from clight, mu0 in derivqty
euler         real /0.57721566490153386/  # Euler-Masceroni constant
jperev        real /0.0/    [J/eV]        # Conversion factor, Joules per eV
                                          # Set to echarge in derivqty
mu0           real /0.0/    [H/m]         # Permeability of free space
boltzmann     real /1.3806503e-23/ [J/K]  # Boltzmann's constant
fuz           real /3.e-5/                # for integer-real comparisons

**** Ch_var dump:
# Mostly character variables associated with run identification
# Top label pline3 is loaded with (it, time, zbeam) when pic packages active
# Middle lines pline2, pline1 are available to the user
# Bottom label line is loaded with (codeid, runid, date, time, runmaker)
numframe integer /0/ -dump    # current frame number for frame index
frameti  character*(11)       # ascii "Stepnnnnn" for frame index
framett  character*(26)       # ascii top    title for frame index
frametr  character*(30)       # ascii right  title for frame index
frametb  character*(8)        # ascii bottom title for frame index
frametl  character*(8)        # ascii left   title for frame index
pline3   character*(59)       # 3 above bottom (big, loaded by W6 w/ IT etc.)
pline2   character*(80)       # 2 above bottom (usually main label)
pline1   character*(80)       # 1 above bottom (usually aux. label)
runid    character*(40) /"warp"/ # Four-character run name
rundate  character*(8)  /" "/    # Run date
runtime  character*(8)  /" "/    # Run time
runmaker character*(35) /" "/    # Name of person running code / special notes

*********** Lattice dump:
# Arrays describing the focusing and bending lattice
# Element starts must fall in [0,zlatperi), except perhaps for element 0.
# Element ends must fall in (0,zlatperi], except for last element.
# Element 0 is set using periodicity if user doesn't set it.
zlatperi  real    /0./    [m] # Periodicity length of lattice
zlatstrt  real    /0./    [m] # Z of lattice start (added to element z's)
zlatbuffer real   /0./    [m] # Buffer added to element lengths so nearby
                              # elements are considered overlapping.
acclzstt  real    /-1.e9/ [m] # Z where accl gaps start
acclbeamframe real /0./   [m] # Location in beam frame where gap accelerates
                              # the beam frame.
lacclzl   logical /.false./   # When true, accelerating gaps are zero length
dipotype character*8 /"Userset"/ # Use "box" to autoset to box dipoles
iqerr     integer /0/         # index of first quad for which a position error
                              # is stored in qoffx,y
qerrxrms  real    /0./    [m] # rms x position error of quads that have errors
qerryrms  real    /0./    [m] # rms y position error of quads that have errors
ndrft     integer /NELEMENT/  # No. of drift elements in lattice
nbend     integer /NELEMENT/  # No. of bend elements in lattice
ndipo     integer /NELEMENT/  # No. of dipole elements in lattice
nquad     integer /NELEMENT/  # No. of quadrupole elements in lattice
ntquad    integer /0/         # No. data points for time dependent quad fields
nsext     integer /NELEMENT/  # No. of sextupole elements in lattice
nhele     integer /NELEMENT/  # No. of hard-edge multipole elements in lattice
nhmlt     integer /25/        # No. of h.e. multipoles in each h.e. element
nqerr     integer /1000/      # No. of quad elements for which error
                              # data is stored
naccl     integer /NELEMENT/  # No. of accelerator elements in lattice
ntaccl    integer /0/   # No. times for which the gap field is stored in acclet
nemlt     integer /NELEMENT/  # No. of electric multipole elements in lattice
nmmlt     integer /NELEMENT/  # No. of magnetic multipole elements in lattice
neerr     integer /1000/      # No. of electric mult elements for which error
                              # data is stored
nmerr     integer /1000/      # No. of magnetic mult elements for which error
                              # data is stored
nbgrd     integer /NELEMENT/  # No. of elements with B field on a 3-D grid
npgrd     integer /NELEMENT/  # No. of elements with potential on a 3-D grid
drftzs(0:ndrft)   _real [m]   # Z's of drift starts
drftze(0:ndrft)   _real [m]   # Z's of drift ends
drftap(0:ndrft)   _real [m]   # Aperture in drifts
drftox(0:ndrft)   _real [m]   # X-offsets of drifts
drftoy(0:ndrft)   _real [m]   # Y-offsets of drifts
drftol(0:ndrft)   _integer    # Overlap level of the element (autoset)
bendzs(0:nbend)   _real [m]   # Z's of bend starts
bendze(0:nbend)   _real [m]   # Z's of bend ends
bendrc(0:nbend)   _real [m]   # Radii of curvature of bends
bendap(0:nbend)   _real [m]   # Aperture in bends
bendox(0:nbend)   _real [m]   # X-offsets of bends (for aperture only)
bendoy(0:nbend)   _real [m]   # Y-offsets of bends (for aperture only)
bendol(0:nbend)   _integer    # Overlap level of the element (autoset)
dipozs(0:ndipo)   _real [m]   # Z's of dipo starts (set from bendzs if =dipoze)
dipoze(0:ndipo)   _real [m]   # Z's of dipo ends   (set from bendze if =dipozs)
dipoby(0:ndipo)   _real [T]   # By's of dipos (set from bendrc if 0 & dipoex=0)
dipobx(0:ndipo)   _real [T]   # Bx's of dipos
dipota(0:ndipo)   _real [1]   # Tan of dipo entry face angle (auto-set if 0)
dipotb(0:ndipo)   _real [1]   # Tan of dipo exit  face angle (auto-set if 0)
dipoex(0:ndipo)   _real [V/m] # Ex's of dipos
dipoey(0:ndipo)   _real [V/m] # Ey's of dipos
dipoap(0:ndipo)   _real [m]   # Aperture in dipoles
dipox1(0:ndipo)   _real [m]   # X location of first dipole plates
dipox2(0:ndipo)   _real [m]   # X location of second dipole plates
dipov1(0:ndipo)   _real [V]   # Voltage of first dipole plates
dipov2(0:ndipo)   _real [V]   # Voltage of second dipole plates
dipol1(0:ndipo)   _real [m]   # Length of first dipole plates
dipol2(0:ndipo)   _real [m]   # Length of second dipole plates
dipow1(0:ndipo)   _real [m]   # Width of first dipole plates
dipow2(0:ndipo)   _real [m]   # Width of second dipole plates
dipool(0:ndipo)   _integer    # Overlap level of the element (autoset)
quadzs(0:nquad)   _real [m]   # Z's of quad starts (hard-edge measure) 
quadze(0:nquad)   _real [m]   # Z's of quad ends   (hard-edge measure)
quaddb(0:nquad)   _real [T/m] # Magnetic quad strengths (field gradients)
quadde(0:nquad)   _real [V/m**2] # Electric quad strengths (field gradients)
quadet(0:ntquad,0:nquad) _real [V/m**2] # Electric quad strengths as a
                                        # function of time
quadbt(0:ntquad,0:nquad) _real [V/m**2] # Magnetic quad strengths as a
                                        # function of time
quadts(0:nquad)   _real [t]   # Time of start of quad field in quadet, quadbt
quaddt(0:nquad)   _real [t]   # Delta t for quad field data in quadet, quadbt
quadvx(0:nquad)   _real [V]   # Voltage of electric quad on x axis
quadvy(0:nquad)   _real [V]   # Voltage of electric quad on y axis
quadap(0:nquad)   _real [m]   # Aperture of quad
quadrr(0:nquad)   _real [m]   # Radius of electrostatic quadrupole rod
quadrl(0:nquad)   _real [m]   # Length of electrostatic quadrupole rod
quadgl(0:nquad)   _real [m]   # Length of electrostatic quadrupole gap
quadgp(0:nquad)   _real [ ]   # Gap position of ESQ, only sign is used
quadpw(0:nquad)   _real [m]   # End plate width of electrostatic quadrupole
quadpa(0:nquad)   _real [m]   # End plate aperture of electrostatic quadrupole
quadpr(0:nquad)   _real [m]   # End plate max radius
quadsl(0:nquad)   _real [m]   # Slant of rod which makes it a cone
quadol(0:nquad)   _integer    # Overlap level of the element (autoset)
quaddo(0:nquad)   _real [1]   # Relative strength of dodecopole at pole tip
qdelglx(0:nquad)  _real [m]   # Change in gap length on x axis
qdelgly(0:nquad)  _real [m]   # Change in gap length on y axis
qdelaxp(0:nquad)  _real [m]   # Change in aperture of rod on plus  x axis
qdelaxm(0:nquad)  _real [m]   # Change in aperture of rod on minus x axis
qdelayp(0:nquad)  _real [m]   # Change in aperture of rod on plus  y axis
qdelaym(0:nquad)  _real [m]   # Change in aperture of rod on minus y axis
qdelrxp(0:nquad)  _real [m]   # Change in radius of rod on plus  x axis
qdelrxm(0:nquad)  _real [m]   # Change in radius of rod on minus x axis
qdelryp(0:nquad)  _real [m]   # Change in radius of rod on plus  y axis
qdelrym(0:nquad)  _real [m]   # Change in radius of rod on minus y axis
qdelvxp(0:nquad)  _real [m]   # Change in voltage of rod on plus  x axis
qdelvxm(0:nquad)  _real [m]   # Change in voltage of rod on minus x axis
qdelvyp(0:nquad)  _real [m]   # Change in voltage of rod on plus  y axis
qdelvym(0:nquad)  _real [m]   # Change in voltage of rod on minus y axis
qdeloxp(0:nquad)  _real [m]   # Perpendicular offset of rod on plus  x axis
qdeloxm(0:nquad)  _real [m]   # Perpendicular offset of rod on minus x axis
qdeloyp(0:nquad)  _real [m]   # Perpendicular offset of rod on plus  y axis
qdeloym(0:nquad)  _real [m]   # Perpendicular offset of rod on minus y axis
qdelpwl(0:nquad)  _real [m]   # Change on left  plate width
qdelpwr(0:nquad)  _real [m]   # Change on right plate width
qdelpal(0:nquad)  _real [m]   # Change on left  plate aperture
qdelpar(0:nquad)  _real [m]   # Change on right plate aperture
qdelprl(0:nquad)  _real [m]   # Change on left  plate max radius
qdelprr(0:nquad)  _real [m]   # Change on right plate max radius
qoffx(0:nqerr)    _real [m]   # Quad offsets in x (note dim. nqerr) (auto. set)
qoffy(0:nqerr)    _real [m]   # Quad offsets in y (note dim. nqerr) (auto. set)
sextzs(0:nsext)   _real [m]   # Z's of sextupole starts
sextze(0:nsext)   _real [m]   # Z's of sextupole ends
sextdb(0:nsext)   _real [ ]   # d^2 B/dx^2 field of sextupole (6*V33)
sextde(0:nsext)   _real [ ]   # d^2 E/dx^2 field of sextupole (6*V33)
sextol(0:nsext)   _integer    # Overlap level of the element (autoset)
helezs(0:nhele)   _real   [m] # Z's of hard-edge (h.e.) multipole starts
heleze(0:nhele)   _real   [m] # Z's of hard-edge (h.e.) multipole ends
heleap(0:nhele)   _real [m]   # Aperture in hard-edge elements
heleae(1:nhmlt,0:nhele)   _real    [var]
                              # amplitude of h.e. electric multipoles
heleam(1:nhmlt,0:nhele)   _real    [var]
                              # amplitude of h.e. magnetic multipoles
heleep(1:nhmlt,0:nhele)   _real    [var]
                              # axial derivative of amplitude of h.e.
                              # electric multipoles
helemp(1:nhmlt,0:nhele)   _real    [var]
                              # axial derivative of amplitude of h.e.
                              # magnetic multipoles
hele_n(1:nhmlt,0:nhele)   _real    [1]
                              # Harmonic number of h.e. multipoles
hele_v(1:nhmlt,0:nhele)   _real    [1]
                              # Order of pseudomultipole of h.e. multipoles
helepe(1:nhmlt,0:nhele)   _real    [rad]
                              # phase angle of h.e. electric multipoles
helepm(1:nhmlt,0:nhele)   _real    [rad]
                              # phase angle of h.e. magnetic multipoles
helene(0:nhele)   _integer [1] # no. elec multipoles in h.e. element (auto-set)
helenm(0:nhele)   _integer [1] # no. mag  multipoles in h.e. element (auto-set)
heleox(0:nhele)   _real    [m] # X-offsets of h.e. multipole centers
heleoy(0:nhele)   _real    [m] # Y-offsets of h.e. multipole centers
helerr(0:nhele)   _real [m]   # Radius of electrostatic quadrupole rod
helerl(0:nhele)   _real [m]   # Length of electrostatic quadrupole rod
helegl(0:nhele)   _real [m]   # Length of electrostatic quadrupole gap
helegp(0:nhele)   _real [ ]   # Gap position of ESQ, only sign is used
helepw(0:nhele)   _real [m]   # End plate width of electrostatic quadrupole
helepa(0:nhele)   _real [m]   # End plate aperture of electrostatic quadrupole
heleol(0:nhele)   _integer    # Overlap level of the element (autoset)
emltzs(0:nemlt)   _real [m]   # Z's of electric multipole element starts
emltze(0:nemlt)   _real [m]   # Z's of electric multipole element ends
emltap(0:nemlt)   _real [m]   # Aperture in electric multipole elements
emltph(0:nemlt)   _real [rad] # Phase angle of electric multipole element field
emltsf(0:nemlt)   _real [1] /0./ # Scale factor for electric multipole element
                                 # Field is scaled by (emltsc+emltsf)
emltsc(0:nemlt)   _real [1] /1./ # Scale factor for electric multipole element
                                 # Field is scaled by (emltsc+emltsf)
emltid(0:nemlt)   _integer    # Index of electric multipole dataset 
emltox(0:neerr)   _real [m]   # Offset in x of electric multipole centers
emltoy(0:neerr)   _real [m]   # Offset in y of electric multipole centers
emltrr(0:nemlt)   _real [m]   # Radius of electrostatic quadrupole rod
emltrl(0:nemlt)   _real [m]   # Length of electrostatic quadrupole rod
emltgl(0:nemlt)   _real [m]   # Length of electrostatic quadrupole gap
emltgp(0:nemlt)   _real [ ]   # Gap position of ESQ, only sign is used
emltpw(0:nemlt)   _real [m]   # End plate width of electrostatic quadrupole
emltpa(0:nemlt)   _real [m]   # End plate aperture of electrostatic quadrupole
emltol(0:nemlt)   _integer    # Overlap level of the element (autoset)
mmltzs(0:nmmlt)   _real [m]   # Z's of magnetic multipole element starts
mmltze(0:nmmlt)   _real [m]   # Z's of magnetic multipole element ends
mmltap(0:nmmlt)   _real [m]   # Aperture in magnetic multipole elements
mmltph(0:nmmlt)   _real [rad] # Phase angle of magnetic multipole element field
mmltsf(0:nmmlt)   _real [1] /0./ # Scale factor for magnetic multipole element
                                 # Field is scaled by (mmltsc+mmltsf)
mmltsc(0:nmmlt)   _real [1] /1./ # Scale factor for magnetic multipole element
                                 # Field is scaled by (mmltsc+mmltsf)
mmltid(0:nmmlt)   _integer    # Index of magnetic multipole dataset
mmltox(0:nmerr)   _real [m]   # Offset in x of magnetic multipole centers
mmltoy(0:nmerr)   _real [m]   # Offset in y of magnetic multipole centers
mmltol(0:nmmlt)   _integer    # Overlap level of the element (autoset)
acclzs(0:naccl)   _real [m]   # Z's of acceleration gap starts
acclze(0:naccl)   _real [m]   # Z's of acceleration gap ends
acclez(0:naccl)   _real [V/m] # Ez's of acceleration gaps, constant part
acclap(0:naccl)   _real [m]   # Aperture in acceleration gaps
acclox(0:naccl)   _real [m]   # Offset in x of accl elements
accloy(0:naccl)   _real [m]   # Offset in y of accl elements
acclxw(0:naccl) _real [V/m^2] # Weights for linear x components of gap Ez's
# Acceleration gap Ez's are of the form, E_z = acclez + acclxw*x.
acclsw(0:naccl)   _real [1]   # Switch for grid accel by gaps (0 = on,1 = off)
acclet(0:ntaccl,0:naccl) _real [V/m] # Ez's of acceleration gaps as a function
                                     # of time.
acclts(0:naccl)   _real [t]   # Time of start of gap field in acclet
accldt(0:naccl)   _real [t]   # Delta t for gap field data in acclet
acclol(0:naccl)   _integer    # Overlap level of the element (autoset)
bgrdzs(0:nbgrd)   _real [m]   # Z starts of 3-D grid of B field data (BGRDdata)
bgrdze(0:nbgrd)   _real [m]   # Z ends of 3-D grid of B field data (BGRDdata)
bgrdxs(0:nbgrd)   _real [m]   # X starts of 3-D grid of B field data (BGRDdata)
bgrdys(0:nbgrd)   _real [m]   # Y starts of 3-D grid of B field data (BGRDdata)
bgrdap(0:nbgrd)   _real [m]   # Aperture in bgrd elements
bgrdox(0:nbgrd)   _real [m]   # Offset in x of bgrd elements
bgrdoy(0:nbgrd)   _real [m]   # Offset in y of bgrd elements
bgrdph(0:nbgrd)   _real [rad] # Phase angle of bgrd elements 
bgrdsp(0:nbgrd)   _real [1]   # Sine   of bgrdph (auto-set) 
bgrdcp(0:nbgrd)   _real [1]   # Cosine of bgrdph (auto-set) 
bgrdid(0:nbgrd)   _integer    # Index of to 3-D B field data sets (BGRDdata)
bgrdsf(0:nbgrd) _real [1] /0./ # Scale factor to multiply 3-D B field data set
                               # BGRDdata. Field is scaled by (bgrdsc+bgrdsf)
bgrdsc(0:nbgrd) _real [1] /1./ # Scale factor to multiply 3-D B field data set
                               # BGRDdata. Field is scaled by (bgrdsc+bgrdsf)
bgrdsy(0:nbgrd) _integer /0/   # Level of symmetry in the bgrd data.
                               # (0, no symmetry; 2, quadrupole)
                               # Defaul is no symmetry.
bgrdol(0:nbgrd)   _integer    # Overlap level of the element (autoset)
pgrdzs(0:npgrd)   _real [m]   # Z starts of 3-D grid of potential data(PGRDdata)
pgrdze(0:npgrd)   _real [m]   # Z ends of 3-D grid of potential data(PGRDdata)
pgrdxs(0:npgrd)   _real [m]   # X starts of 3-D grid of potential data(PGRDdata)
pgrdys(0:npgrd)   _real [m]   # Y starts of 3-D grid of potential data(PGRDdata)
pgrdap(0:npgrd)   _real [m]   # Aperture in pgrd elements
pgrdox(0:npgrd)   _real [m]   # Offset in x of pgrd elements
pgrdoy(0:npgrd)   _real [m]   # Offset in y of pgrd elements
pgrdph(0:npgrd)   _real [rad] # Phase angle of pgrd elements
pgrdsp(0:npgrd)   _real [1]   # Sine   of pgrdph (auto-set) 
pgrdcp(0:npgrd)   _real [1]   # Cosine of pgrdph (auto-set)  
pgrdid(0:npgrd)   _integer    # Index of to 3-D potential data sets (PGRDdata)
pgrdsf(0:npgrd) _real [1] /0./ # Scale factor to multiply 3-D potential data set
                               # PGRDdata. Field is scaled by (pgrdsc+pgrdsf)
pgrdsc(0:npgrd) _real [1] /1./ # Scale factor to multiply 3-D potential data set
                               # PGRDdata. Field is scaled by (pgrdsc+pgrdsf)
pgrdrr(0:npgrd)   _real [m]   # Radius of electrostatic quadrupole rod
pgrdrl(0:npgrd)   _real [m]   # Length of electrostatic quadrupole rod
pgrdgl(0:npgrd)   _real [m]   # Length of electrostatic quadrupole gap
pgrdgp(0:npgrd)   _real [ ]   # Gap position of ESQ, only sign is used
pgrdpw(0:npgrd)   _real [m]   # End plate width of electrostatic quadrupole
pgrdpa(0:npgrd)   _real [m]   # End plate aperture of electrostatic quadrupole
pgrdol(0:npgrd)   _integer    # Overlap level of the element (autoset)
drfts     logical             # Flag for existence of drfts (auto set)
bends     logical             # Flag for existence of bends (auto set)
dipos     logical             # Flag for existence of dipos (auto set)
quads     logical             # Flag for existence of quads (auto set)
sexts     logical             # Flag for existence of sexts (auto set)
heles     logical             # Flag for exist. of hard-edge mults (auto set)
accls     logical             # Flag for existence of accel (auto set)
emlts     logical             # Flag for existence of emlts (auto set)
mmlts     logical             # Flag for existence of mmlts (auto set)
bgrds     logical             # Flag for existence of bgrds (auto set)
pgrds     logical             # Flag for existence of pgrds (auto set)
diposet   logical  /.true./   # Auto-set dipoles from bend locations and radii 

******** Mult_data dump:
nemltsets          integer /0/  # Number of different electrostatic data sets
nesmult            integer /0/  # Number of electrostatic multipoles
nzemltmax          integer /0/  # Maximum number of multipole Z data points
nzemlt(nemltsets) _integer      # Number of multipole Z data points
dzemlt(nemltsets) _real         # Grid cell size for multipole data
esemlt(0:nzemltmax,nesmult,nemltsets) _real
                               # Electrostatic multipoles as a function of Z
esemltp(0:nzemltmax,nesmult,nemltsets) _real
                               # Axial derivatives of electrostatic multipoles
                               # as a function of Z.  If not set by user,
                               # autoset via finite difference of esemlt.
esemltph(0:nzemltmax,nesmult,nemltsets) _real
                               # Phase angle of electric multipoles
esemltphp(0:nzemltmax,nesmult,nemltsets) _real
                               # Axial derivatives of phase angle of electric
                               # multipoles.  If not set by user,
                               # autoset via finite difference of esemltph.
emlt_n(nesmult)   _real [1]    # Harmonic number of electric multipoles
emlt_v(nesmult)   _real [1]    # Order of pseudomultipole of electric multipoles

nmmltsets          integer /0/  # Number of magnetostatic multipole data sets
nmsmult            integer /0/  # Number of magnetostatic multipoles
nzmmltmax          integer /0/  # Maximum number of multipole Z data points
nzmmlt(nmmltsets) _integer      # Number of multipole Z data points
dzmmlt(nmmltsets) _real         # Grid cell size for multipole data
msmmlt(0:nzmmltmax,nmsmult,nmmltsets) _real
                               # Magnetostatic multipoles as a function of Z
msmmltp(0:nzmmltmax,nmsmult,nmmltsets) _real
                               # Axial derivatives of magnetostatic multipoles
                               # as a function of Z.  If not set by user,
                               # autoset via finite difference of msmmlt.
msmmltph(0:nzmmltmax,nmsmult,nmmltsets) _real
                               # Phase angle of magnetic multipoles
msmmltphp(0:nzmmltmax,nmsmult,nmmltsets) _real
                               # Axial derivatives of phase angle of magnetic
                               # multipoles.  If not set by user,
                               # autoset via finite difference of msmmltph.
mmlt_n(nmsmult)   _real [1]    # Harmonic number of magnetic multipoles
mmlt_v(nmsmult)   _real [1]    # Order of pseudomultipole of magnetic multipoles

******** BGRDdata dump:
# Data for the 3-D B field lattice element
bgrdnx integer /0/ # Number of X cells
bgrdny integer /0/ # Number of Y cells
bgrdnz integer /0/ # Number of Z cells
bgrdns integer /0/ # Number of data sets
bgrddx(bgrdns)                      _real [m]   # X cell size
bgrddy(bgrdns)                      _real [m]   # Y cell size
bgrddz(bgrdns)                      _real [m]   # Z cell size
bgrddxi(bgrdns)                     _real [1/m] # 1 over X cell size (autoset)
bgrddyi(bgrdns)                     _real [1/m] # 1 over Y cell size (autoset)
bgrddzi(bgrdns)                     _real [1/m] # 1 over Z cell size (autoset)
bgrdbx(0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns) _real [T] # Bx
bgrdby(0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns) _real [T] # By
bgrdbz(0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns) _real [T] # Bz

******** PGRDdata dump:
# Data for the 3-D potential lattice element
pgrdnx integer /0/ # Number of X cells
pgrdny integer /0/ # Number of Y cells
pgrdnz integer /0/ # Number of Z cells
pgrdns integer /0/ # Number of data sets
pgrddx(pgrdns)                      _real [m]   # X cell size
pgrddy(pgrdns)                      _real [m]   # Y cell size
pgrddz(pgrdns)                      _real [m]   # Z cell size
pgrddxi(pgrdns)                     _real [1/m] # 1 over X cell size (autoset)
pgrddyi(pgrdns)                     _real [1/m] # 1 over Y cell size (autoset)
pgrddzi(pgrdns)                     _real [1/m] # 1 over Z cell size (autoset)
pgrd(0:pgrdnx,0:pgrdny,-1:pgrdnz+1,pgrdns) _real [V] # Potential

*********** LatticeInternal dump:
# Internal lattice arrays, all derived from Lattice data
# nzl is set to nz by pkg w3d (etc.) at generation
dzl                   real [m]     # LatticeInternal mesh grid cell size
dzli                  real [m]     # LatticeInternal grid cell size inverse
zlframe               real [m]     # Location of LatticeInternal frame
zltime                real [m]     # Time of LatticeInternal frame
zlmin                 real [m]     # LatticeInternal mesh maximum in z
zlmax                 real [m]     # LatticeInternal mesh minimum in z
nzl                integer /0/ [1] # Number of LatticeInternal points
nzlmax             integer /0/ [1] # Length of LatticeInternal arrays
zlmesh(0:nzlmax)     _real [m]     # LatticeInternal Z mesh
ndrftol            integer /1/ # Maximum level of overlapping drft elements
odrftoi(0:ndrft)   _integer     # Overlap indices for drft elements
odrftio(0:ndrft)   _integer     # Overlap indices for drft elements
odrftii(ndrftol)   _integer     # Overlap indices for drft elements
odrftnn(ndrftol)   _integer     # Number of drft elements in overlap levels
nbendol            integer /1/ # Maximum level of overlapping bend elements
obendoi(0:nbend)   _integer     # Overlap indices for bend elements
obendio(0:nbend)   _integer     # Overlap indices for bend elements
obendii(nbendol)   _integer     # Overlap indices for bend elements
obendnn(nbendol)   _integer     # Number of bend elements in overlap levels
ndipool            integer /1/ # Maximum level of overlapping dipo elements
odipooi(0:ndipo)   _integer     # Overlap indices for dipo elements
odipoio(0:ndipo)   _integer     # Overlap indices for dipo elements
odipoii(ndipool)   _integer     # Overlap indices for dipo elements
odiponn(ndipool)   _integer     # Number of dipo elements in overlap levels
nquadol            integer /1/ # Maximum level of overlapping quad elements
oquadoi(0:nquad)   _integer     # Overlap indices for quad elements
oquadio(0:nquad)   _integer     # Overlap indices for quad elements
oquadii(nquadol)   _integer     # Overlap indices for quad elements
oquadnn(nquadol)   _integer     # Number of quad elements in overlap levels
nsextol            integer /1/ # Maximum level of overlapping sext elements
osextoi(0:nsext)   _integer     # Overlap indices for sext elements
osextio(0:nsext)   _integer     # Overlap indices for sext elements
osextii(nsextol)   _integer     # Overlap indices for sext elements
osextnn(nsextol)   _integer     # Number of sext elements in overlap levels
nheleol            integer /1/ # Maximum level of overlapping hele elements
oheleoi(0:nhele)   _integer     # Overlap indices for hele elements
oheleio(0:nhele)   _integer     # Overlap indices for hele elements
oheleii(nheleol)   _integer     # Overlap indices for hele elements
ohelenn(nheleol)   _integer     # Number of hele elements in overlap levels
nemltol            integer /1/ # Maximum level of overlapping emlt elements
oemltoi(0:nemlt)   _integer     # Overlap indices for emlt elements
oemltio(0:nemlt)   _integer     # Overlap indices for emlt elements
oemltii(nemltol)   _integer     # Overlap indices for emlt elements
oemltnn(nemltol)   _integer     # Number of emlt elements in overlap levels
nmmltol            integer /1/ # Maximum level of overlapping mmlt elements
ommltoi(0:nmmlt)   _integer     # Overlap indices for mmlt elements
ommltio(0:nmmlt)   _integer     # Overlap indices for mmlt elements
ommltii(nmmltol)   _integer     # Overlap indices for mmlt elements
ommltnn(nmmltol)   _integer     # Number of mmlt elements in overlap levels
nacclol            integer /1/ # Maximum level of overlapping accl elements
oaccloi(0:naccl)   _integer     # Overlap indices for accl elements
oacclio(0:naccl)   _integer     # Overlap indices for accl elements
oacclii(nacclol)   _integer     # Overlap indices for accl elements
oacclnn(nacclol)   _integer     # Number of accl elements in overlap levels
nbgrdol            integer /1/ # Maximum level of overlapping bgrd elements
obgrdoi(0:nbgrd)   _integer     # Overlap indices for bgrd elements
obgrdio(0:nbgrd)   _integer     # Overlap indices for bgrd elements
obgrdii(nbgrdol)   _integer     # Overlap indices for bgrd elements
obgrdnn(nbgrdol)   _integer     # Number of bgrd elements in overlap levels
npgrdol            integer /1/ # Maximum level of overlapping pgrd elements
opgrdoi(0:npgrd)   _integer     # Overlap indices for pgrd elements
opgrdio(0:npgrd)   _integer     # Overlap indices for pgrd elements
opgrdii(npgrdol)   _integer     # Overlap indices for pgrd elements
opgrdnn(npgrdol)   _integer     # Number of pgrd elements in overlap levels
cdrftzs(0:nzlmax,ndrftol)    _real [m]     # by z, Z's of drft starts
cdrftze(0:nzlmax,ndrftol)    _real [m]     # by z, Z's of drft ends
cdrftid(0:nzlmax,ndrftol) _integer         # by z, Index to drft arrays
cbendzs(0:nzlmax)            _real [m]     # by z, Z's of bend starts
cbendze(0:nzlmax)            _real [m]     # by z, Z's of bend ends
cbendrc(0:nzlmax)            _real [m]     # by z, Radii of curvature of bends
cbendid(0:nzlmax)         _integer         # by z, Index to bend arrays
cdipozs(0:nzlmax,ndipool)    _real [m]     # by z, Z's of dipo starts 
cdipoze(0:nzlmax,ndipool)    _real [m]     # by z, Z's of dipo ends  
cdipoby(0:nzlmax,ndipool)    _real [T]     # by z, By's of dipos
cdipobx(0:nzlmax,ndipool)    _real [T]     # by z, Bx's of dipos
cdipota(0:nzlmax,ndipool)    _real [1]     # by z, tan's of dipos entry angles
cdipotb(0:nzlmax,ndipool)    _real [1]     # by z, tan's of dipos exit angles
cdipoex(0:nzlmax,ndipool)    _real [V/m]   # by z, Ex's of dipos
cdipoey(0:nzlmax,ndipool)    _real [V/m]   # by z, Ey's of dipos
cdipoid(0:nzlmax,ndipool) _integer         # by z, Index to dipo arrays
cquadzs(0:nzlmax,nquadol)    _real [m]     # by z, Z's of quad starts
cquadze(0:nzlmax,nquadol)    _real [m]     # by z, Z's of quad ends
cquaddb(0:nzlmax,nquadol)    _real [ ]     # by z, Magnetic quad fld gradients
cquadde(0:nzlmax,nquadol)    _real [ ]     # by z, Electric quad fld gradients
cquadvx(0:nzlmax,nquadol)    _real [V]     # by z, Voltage on x axis rods
cquadvy(0:nzlmax,nquadol)    _real [V]     # by z, Voltage on y axis rods
cquadid(0:nzlmax,nquadol) _integer         # by z, Index to quad arrays
cqoffx(0:nzlmax,nquadol)     _real [m]     # by z, quad offset in x
cqoffy(0:nzlmax,nquadol)     _real [m]     # by z, quad offset in y
csextzs(0:nzlmax,nsextol)    _real [m]     # by z, Z's of sextupole starts
csextze(0:nzlmax,nsextol)    _real [m]     # by z, Z's of sextupole ends
csextdb(0:nzlmax,nsextol)    _real [ ]     # by z, d^2 B/dx^2 of sext (6*V33)
csextde(0:nzlmax,nsextol)    _real [ ]     # by z, d^2 E/dx^2 of sext (6*V33)
csextid(0:nzlmax,nsextol) _integer         # by z, Index to sext arrays
cheleid(0:nzlmax,nheleol) _integer         # by z, Index to hele arrays
chelezs(0:nzlmax,nheleol)    _real [m]     # by z, Z's of h.e. element starts
cheleze(0:nzlmax,nheleol)    _real [m]     # by z, Z's of h.e. element ends
cemltzs(0:nzlmax,nemltol)    _real [m]     # by z, Z's of electric mult starts
cemltze(0:nzlmax,nemltol)    _real [m]     # by z, Z's of electric mult ends
cemltph(0:nzlmax,nemltol)    _real [rad]   # by z, Phase of electric mult
cemltsf(0:nzlmax,nemltol)    _real [1]     # by z, Scale factor of e.s. mult
cemltsc(0:nzlmax,nemltol)    _real [1]     # by z, Scale factor of e.s. mult
cemltid(0:nzlmax,nemltol) _integer         # by z, Index of emlt arrays
cemltim(0:nzlmax,nemltol) _integer         # by z, Index of multipole dataset
cemltox(0:nzlmax,nemltol)    _real [m]     # by z, Offset in x of mult centers
cemltoy(0:nzlmax,nemltol)    _real [m]     # by z, Offset in y of mult centers
cmmltzs(0:nzlmax,nmmltol)    _real [m]     # by z, Z's of magnetic mult starts
cmmltze(0:nzlmax,nmmltol)    _real [m]     # by z, Z's of magnetic mult ends
cmmltph(0:nzlmax,nmmltol)    _real [rad]   # by z, Phase of magnetic mult
cmmltsf(0:nzlmax,nmmltol)    _real [1]     # by z, Scale factor of m.s. mult
cmmltsc(0:nzlmax,nmmltol)    _real [1]     # by z, Scale factor of m.s. mult
cmmltid(0:nzlmax,nmmltol) _integer         # by z, Index of mmlt arrays
cmmltim(0:nzlmax,nmmltol) _integer         # by z, Index of multipole dataset
cmmltox(0:nzlmax,nmmltol)    _real [m]     # by z, Offset in x of mult centers
cmmltoy(0:nzlmax,nmmltol)    _real [m]     # by z, Offset in y of mult centers
cacclzs(0:nzlmax,nacclol)    _real [m]     # by z, Z's of accelerator gap start
cacclze(0:nzlmax,nacclol)    _real [m]     # by z, Z's of accelerator gap ends
cacclez(0:nzlmax,nacclol)    _real [V/m]   # by z, const Ez's of accl gaps
cacclxw(0:nzlmax,nacclol)    _real [V/m^2] # by z, x-weights for accl gap Ez's 
cacclsw(0:nzlmax,nacclol)    _real [1]     # by z, switch for grid accel by gap
cacclid(0:nzlmax,nacclol) _integer         # by z, Index of accl arrays
cbgrdzs(0:nzlmax,nbgrdol)    _real [m]     # by z, Z's of 3-D B field start
cbgrdze(0:nzlmax,nbgrdol)    _real [m]     # by z, Z's of 3-D B field end
cbgrdid(0:nzlmax,nbgrdol) _integer [1]     # by z, Index of bgrd arrays
cpgrdzs(0:nzlmax,npgrdol)    _real [m]     # by z, Z's of 3-D potential start
cpgrdze(0:nzlmax,npgrdol)    _real [m]     # by z, Z's of 3-D potential end
cpgrdid(0:nzlmax,npgrdol) _integer [1]     # by z, Index of pgrd arrays
lindrft(0:ndrftol)        _logical         # Flag for when drft element in mesh
linbend                    logical         # Flag for when bend element in mesh
lindipo(0:ndipool)        _logical         # Flag for when dipo element in mesh
linquad(0:nquadol)        _logical         # Flag for when quad element in mesh
linsext(0:nsextol)        _logical         # Flag for when sext element in mesh
linhele(0:nheleol)        _logical         # Flag for when hele element in mesh
linemlt(0:nemltol)        _logical         # Flag for when emlt element in mesh
linmmlt(0:nmmltol)        _logical         # Flag for when mmlt element in mesh
linaccl(0:nacclol)        _logical         # Flag for when accl element in mesh
linbgrd(0:nbgrdol)        _logical         # Flag for when bgrd element in mesh
linpgrd(0:npgrdol)        _logical         # Flag for when pgrd element in mesh

*********** Ctl_to_pic:
# Communication between CTL and pic packages.  In TOP since it's "global"
maxcalls integer    # Total number of calls requested within this STEP or RUN
ncall    integer    # Number of this call within the current STEP or RUN

*********** InDiag dump:
# General diagnostic specifications (input qtys)
rwindows(2,0:NWINDOWS)    real  [m]  /0.,99.,TNWINM*0./
   # radial "window" limits for z-vz phase space plots
xwindows(2,0:NWINDOWS)    real  [m]  /-99.,99.,TNWINM*0./
   # "window" limits for y-z phase space plots
ywindows(2,0:NWINDOWS)    real  [m]  /-99.,99.,TNWINM*0./
   # "window" limits for x-z phase space plots
zwindows(2,0:NWINDOWS)    real  [m] /-99.,99.,TNWINM*0./
   # "window" limits for x-y phase space plots, emittance calcs
   #  Window 0 is set to full mesh at generation
xplmin   real  [m]          /0./ # Minimum x for plots; if 0, pkgs should set
xplmax   real  [m]          /0./ # Maximum x for plots; if 0, pkgs should set
yplmin   real  [m]          /0./ # Minimum y for plots; if 0, pkgs should set
yplmax   real  [m]          /0./ # Maximum y for plots; if 0, pkgs should set
zplmin   real  [m]          /0./ # Minimum z for plots; if 0, pkgs should set
zplmax   real  [m]          /0./ # Maximum z for plots; if 0, pkgs should set
xpplmin  real  [1]        /-.04/ # Minimum x' for plots
xpplmax  real  [1]        / .04/ # Maximum x' for plots
ypplmin  real  [1]        /-.04/ # Minimum y' for plots
ypplmax  real  [1]        / .04/ # Maximum y' for plots
vzrng    real  [1]          /0./ # Percent range of vz for plots; if 0, autoset
vzshift  real  [m/s]        /0./ # Shift of vz for plots
eplmin   real  [E]          /0./ # Minimun E, for E field plots
eplmax   real  [E]          /0./ # Maximun E, for E field plots
rhoplmin real  [Coul/m**3]  /0./ # Minimun rho, for rho plots
rhoplmax real  [Coul/m**3]  /0./ # Maximun rho, for rho plots
epsplfac real  [1]        /0.05/ # Factor used in plot limits of epsx and y
xxpslope real  [1/m]        /0./ # Slope removed from phase space plots
phiplmin real  [V]          /0./ # Min phi for on axis plots;both zero-autoscale
phiplmax real  [V]          /0./ # Max phi for on axis plots;both zero-autoscale

ifzmmnt                   integer /2/
   # specifies z moments calculation (0:none, 1:global moments only,
   # 2:full z moments(with extrapolation and linear weighting))
laccumulate_zmoments      logical /.false./
   # When true, zmoments are accumulated over multiple time steps.  Note that
   # then the routines which initially zero the arrays and which finish the
   # calculation must be called by hand.
itmomnts(NCONTROL)        integer /NCONTROL*0/
   # time steps to do calculation of z moments and print one-liner of info;
   # first 3 are a do loop
itplps(NCONTROL)          integer /NCONTROL*0/
   # time steps to do full set of phase space plots; first 3 are a do loop
itplfreq(NCONTROL)        integer /NCONTROL*0/
   # time steps to do "frequent" phase space plots; first 3 are a do loop
itplseldom(NCONTROL)      integer /NCONTROL*0/
   # time steps to do "seldom" plots; first 3 are a do loop
itplalways(NCONTROL)      integer /NCONTROL*0/
   # time steps to do "always" plots; first 3 are a do loop
zzmomnts(NCONTROL)        real    /NCONTROL*0/
   # z locations to do calculation of z moments and print one-liner of info;
   # first 3 are a do loop
zzplps(NCONTROL)          real    /NCONTROL*0/
   # z locations to do full set of phase space plots; first 3 are a do loop
zzplfreq(NCONTROL)        real    /NCONTROL*0/
   # z locations to do "frequent" phase space plots; first 3 are a do loop
zzplseldom(NCONTROL)          real    /NCONTROL*0/
   # z locations to do "seldom" plots; first 3 are a do loop
zzplalways(NCONTROL)        real    /NCONTROL*0/
   # z locations to do "always" plots; first 3 are a do loop
nhist                     integer /5/
   # Interval between timesteps at which history data is saved
npplot(NSUBSETS)          integer /50000,20000,5000/
   # Nominal no. of particles in each subset for phase space plots
   # Actual no. plotted is float(ntopick)/inclump - np
inclump(NSUBSETS)         integer /51,51,51/
   # clump size for choosing which ptcls to plot-MUST DIVIDE npsplt(in Pspwork)
charsbig integer /46/ # Number of characters per line for one-to-a-page plots
charssma integer /60/ # Number of characters per line for four-to-a-page plots
labels   logical /.false./ # Determines whether or not contours are labelled
nskipcol integer /5/       # number of particles to skip in color plots
nstepcol integer /10/      # number of sweeps to plots colot plots
ncolor   integer /5/       # number of colors in color plots
iscolor  integer /2/       # starting color number in color plots
lframadv logical /.true./  # When true, advance frame after each fortran plot
ppmark   character*8 /"dot"/ # Marker used in particle plots.
lprntpara logical /.true./ # When true print output parameters at start
loneliner logical /.false./ # When true, calls interpreter function 'oneliner'
                            # otherwise calls fortran function oneliner.
lpsplots logical /.false./ # When true, the compiled psplots calls the
                           # interpreter function 'psplots' with the
                           # appropriate arguments and then returns.
lonedplts logical /.false./ # When true, the compiled onedplts calls the
                            # interpreter function 'onedplts' with the
                            # appropriate arguments and then returns.

*********** InPltCtl dump:
# Controls for when the various plots are made
# Elements -1, -2, ... control "subset" plots
# Elements 0,1,2,3,... control plots for corresponding window
never    integer /NEVER/  # Allows NEVER to be used in the interpreter
seldom   integer /SELDOM/ # Allows SELDOM to be used in the interpreter
always   integer /ALWAYS/ # Allows ALWAYS to be used in the interpreter
ipzx(-NSUBSETS:NWINDOWS)   integer /NWPNSP1*NEVER/ # X vs Z ywindows
ipzy(-NSUBSETS:NWINDOWS)   integer /NWPNSP1*NEVER/ # Y vs Z xwindows
ipzxp(-NSUBSETS:NWINDOWS)  integer /NWPNSP1*NEVER/ # XP vs Z ywindows
ipzyp(-NSUBSETS:NWINDOWS)  integer /NWPNSP1*NEVER/ # YP vs Z xwindows
ipzvz(-NSUBSETS:NWINDOWS)  integer /NWPNSP1*NEVER/ # VZ vs Z rwindows
ipzxy(-NSUBSETS:0)         integer /NSUBSETS*NEVER,NEVER/
    # X vs Z and Y vs Z on one page ... subsets and all-particles only
ipxy(-NSUBSETS:NWINDOWS)   integer /NWPNSP1*NEVER/   # Y vs X zwindows
ipxxp(-NSUBSETS:NWINDOWS)  integer /NWPNSP1*NEVER/ # XP vs X zwindows
ipyyp(-NSUBSETS:NWINDOWS)  integer /NWPNSP1*NEVER/ # YP vs Y zwindows
ipxpyp(-NSUBSETS:NWINDOWS) integer /NWPNSP1*NEVER/ # YP vs XP zwindows
iptrace(0:NWINDOWS)        integer /NEVER,NWINDOWS*NEVER/
ipzx4                      integer /NEVER/ # X vs Z ywindows
ipzy4                      integer /NEVER/ # Y vs Z xwindows
ipzxp4                     integer /NEVER/ # XP vs Z ywindows
ipzyp4                     integer /NEVER/ # YP vs Z xwindows
ipzvz4                     integer /NEVER/ # VZ vs Z rwindows
ipxy4                      integer /NEVER/ # Y vs X zwindows
ipxxp4                     integer /NEVER/ # XP vs X zwindows
ipyyp4                     integer /NEVER/ # YP vs Y zwindows
ipxpyp4                    integer /NEVER/ # YP vs XP zwindow
ipcurr                     integer /NEVER/ # Current vs. z
ipegap                     integer /NEVER/ # NEVER Ez (smeared model) vs. z
iplchg                     integer /NEVER/ # Line charge vs. z
ipvzofz                    integer /NEVER/ # Axial velocity vs. z
ipezax                     integer /NEVER/ # Ez on axis vs. z
ipphiax                    integer /NEVER/ # Potential on axis vs. z
iprhoax                    integer /NEVER/ # Charge density on axis vs. z
ipzxyco                    integer /NEVER/
    # X vs Z and Y vs Z on one page ... in color, skip NSKIPCOL particles
ipzvzco                    integer /NEVER/
    # VZ vs Z in color, skip NSKIPCOL particles
ipxxpco(0:NWINDOWS)        integer /NEVER,NWINDOWS*NEVER/
    # X vs XP in color, skip NSKIPCOL particles
ipyypco(0:NWINDOWS)        integer /NEVER,NWINDOWS*NEVER/
    # Y vs YP in color, skip NSKIPCOL particles

*********** InGaps dump:
# Parameters associated with accelerating gap models
# (so far, just "smeared" gaps)
cgap    real [Farads]  /0./      # Capacitance of an individual gap
lgap    real [Henrys]  /1.e20/   # Inductance of an individual gap
rgap    real [Ohms]    /1.e-20/  # Resistance of an individual gap
zsgap   real [m]       /1./      # Axial separation of gap centers
ifgap   logical        /.false./ # turns gaps on or off

*********** InGen dump:
# General parameters which control the mechanics of the run (input qtys)
dt                        real  [s]  /0./
   # Timestep size
tstop                     real  [s]  /LARGEPOS/
   # Time at which run will stop
zstop                     real  [m]  /LARGEPOS/
   # Z location at which run will stop
bz0                       real [T]   /0./
   # Uniform axial "guide" B field
ez0                       real /0./
   # Linear Ez field for bunching force
vbeamfrm                  real [M/S]
   # Velocity of the beam frame.  Normally the same as vbeam.
allspecl                  logical  /.false./
   # flag for making all time steps "special" ones
depos                     character*8 /"vector"/
   # Specifies charge deposition algorithm, "scalar" or "vector"
laccumulate_rho           logical /.false./
   # When true, rho is accumulated over multiple time steps. (i.e. it is not
   # zeroed each time step before the deposition.)
gamadv                    character*8 /"stndrd"/
   # Specifies type of gamma advance, "stndrd", "fast 1", or "fast 2"
efetch                    integer /1/
   # Specifies type of fetch used to set internal E field
   # 1: direct calculation from phi
   # 2: direct calculation from phi, vectorization over 8 grid cell corners
   # 3: indirect calculation from pre-calculated E (finite differences of phi)
   # Method 3 is generally fastest but requires lots of extra storage space.
   # Next best is 1.
periinz                   logical /.true./
   # Specifies whether or not there is periodicity in z
stickyz                   logical /.false./
   # Specifies whether or not there are sticky walls in z
stickyxy                  logical /.true./
   # Specifies whether or not there are sticky walls in x and y
pbound0                   integer /2/
   # boundary condition for particles at iz = 0 (absorb/dirichlet: absorption, 
   # reflect/neumann: reflection, periodic: periodicity)
pboundnz                  integer /2/
   # boundary condition for particles at iz = nz (absorb/dirichlet: absorption, 
   # reflect/neumann: reflection, periodic: periodicity)
pboundxy                  integer /0/
   # boundary condition for particles at x and y (absorb/dirichlet: absorption, 
   # reflect/neumann: reflection, periodic: periodicity)
ibpush                    integer /1/
   # Specifies type of B-advance: 0-none, 1-fast, 2-with tan(alpha)/alpha
ifeears                   integer /0/
   # Specifies type of Eears: 0-none, 1-linear, 2-Ezax
fstype                    integer /0/
   # Specifies type of field solve (-1:none, 0:normal FFT, 1:pipe, 
   #                                2:cap, 3:psor, 4:2d FFT + tridiag in z)
itdump                    integer    /500/
   # Restart dump interval
nrestart                  character*8  /" "/
   # Name of run from which this one is restarting
nt                        integer  /20000000/
   # Number of timesteps to run
lgridqnt                  logical    /.false./
   # When true, zgrid is always an integer number of dz's nearest zbeam,
   # zgrid = int(zbeam/dz+.5)*dz, and so zgrid quantized number of dz's.
   # When false, zgrid is always equal to zbeam.
lbeamcom                  logical    /.false./
   # When true, zbeam follows the beam center of mass.
relativity                integer /0/
   # Level of relativitistic corrections.
   #  1: scale transverse self E-field by 1/gamma**2
clearlostpart             integer /1/
   # When 0, do not clear lost particles, when 1, swap lost particles with
   # ones at the end of the array (much faster), when 2, shift particles to
   # fill in the gaps (much slower)

*********** InPart dump:
# Particle input quantities (input qtys)
ns              integer         /1/
   # Number of species
npmax_s(0:ns)  _integer [1]     /0/ +parallel
   # Maximum number of particles of each species
   # Index of 0 is used as a guard and always has a value of zero.
np_s(ns)       _integer [1]     /0/
   # Number of particles by species: if zero, set to npmax; if npmax is zero,
   # it is set to sum(np_s).
a0_s(ns)       _real    [m]     /0./
   # Initial beam width in x by species: if zero, set to a0; if a0 is zero,
   # it is set to a0_s(1)
ap0_s(ns)      _real    [rad]   /0./
   # Initial beam envelope vx/vz by species: if zero, set to ap0;if ap0 is zero,
   # it is set to ap0_s(1)
b0_s(ns)       _real    [m]     /0./
   # Initial beam width in y by species: if zero, set to b0; if b0 is zero,
   # it is set to b0_s(1)
bp0_s(ns)      _real    [rad]   /0./
   # Initial beam envelope vy/vz by species: if zero, set to bp0;if bp0 is zero,
   # it is set to bp0_s(1)
x0_s(ns)       _real    [m]     /0./
   # Initial beam centroid in x by species: if zero, set to x0; if x0 is zero,
   # it is set to x0_s(1)
xp0_s(ns)      _real    [rad]   /0./
   # Initial beam centroid vx/vz by species: if zero, set to xp0;if xp0 is zero,
   # it is set to xp0_s(1)
y0_s(ns)       _real    [m]     /0./
   # Initial beam centroid in y by species: if zero, set to y0; if y0 is zero,
   # it is set to y0_s(1)
yp0_s(ns)      _real    [rad]   /0./
   # Initial beam centroid vy/vz by species: if zero, set to yp0;if yp0 is zero,
   # it is set to yp0_s(1)
aion_s(ns)     _real    [1]     /0./
   # Atomic number by species: if zero, set to aion; if aion is zero,
   # it is set to aion_s(1)
ekin_s(ns)     _real    [V]     /0./
   # Particle energy by species: if zero, set to ekin; if ekin is zero,
   # it is set to ekin_s(1)
emit_s(ns)     _real    [m-rad] /0./
   # Emittance by species: if zero, set to emit; if emit is zero,
   # it is set to emit_s(1)
emitx_s(ns)    _real    [m-rad] /0./
   # X-Emittance by species: if zero, set to xemit; if emitx is zero,
   # it is set to emitx_s(1)
emity_s(ns)    _real    [m-rad] /0./
   # Y-Emittance by species: if zero, set to yemit; if emit is zero,
   # it is set to emity_s(1)
emitn_s(ns)    _real    [m-rad] /0./
   # Normalized emittance by species: if zero, set to emitn; if emitn is zero,
   # it is set to emitn_s(1)
emitnx_s(ns)   _real    [m-rad] /0./
   # Normalized X-emittance by species: if zero, set to emitnx; 
   # if emitnx is zero, it is set to emitnx_s(1)
emitny_s(ns)   _real    [m-rad] /0./
   # Normalized Y-emittance by species: if zero, set to emitny; 
   # if emitny is zero, it is set to emitny_s(1)
ibeam_s(ns)    _real    [A]     /0./
   # Beam current by species: if zero, set to ibeam; if ibeam is zero,
   # it is set to ibeam_s(1)
zion_s(ns)     _real    [1]     /0./
   # Charge state by species: if zero, set to zion; if zion is zero,
   # it is set to zion_s(1)
straight_s(ns) _real    [1]     /0./
   # Percent of beam that isn't cigar: if zero, set to straight; if straight
   # is zero, it is set to straight_s
vbeam_s(ns)    _real    [m/s]   /0./
   # Particle axial velocity by species: if zero, set to vbeam; if vbeam is
   # zero, it is set to vbeam_s(1)
vtilt_s(ns)    _real    [1]     /0./
   # Axial velocity tilt by species: if zero, set to vtilt; if vtilt is zero,
   # it is set to vtilt_s(1)
vthperp_s(ns)  _real    [m/s]   /0./
   # Transverse thermal spread by species: if zero, set to vthperp; if vthperp
   # is zero, it is set to vthperp_s(1)
vthz_s(ns)     _real    [m/s]   /0./
   # Axial thermal spread by species: if zero, set to vthz; if vthz is zero,
   # it is set to vthz_s(1)
zimin_s(ns)    _real    [m]     /0./
   # Minimum initial z of beam by species: if zero, set to zimin; if zimin is
   # zero, it is set to zimin_s(1)
zimax_s(ns)    _real    [m]     /0./
   # Maximum initial z of beam by species: if zero, set to zimax; if zimax is
   # zero, it is set to zimax_s(1)
sp_fract(ns)   _real    [1]     /1./
   # For each species, fraction of current that particles will make up.
xcent_s(ns)    _real    [m]     /0./
   # Center in X of each species
ycent_s(ns)    _real    [m]     /0./
   # Center in Y of each species
xpcent_s(ns)   _real    [m]     /0./
   # Center X' of each species
ypcent_s(ns)   _real    [m]     /0./
   # Center Y' of each species

vtilt           real    [1]   /0./
   # Velocity tilt, vz = vbeam*[1+vtilt*(zmid-z)/(zimax-zimin)]
   # thus deltav/vbeam = vtilt
vthperp         real    [m/s] /0./
   # Input transverse thermal velocity (rms).
vthz            real    [m/s] /0./
   # Input z thermal velocity (rms). Fix old decks with vthz=vtz/2 if cigar.
   # Old variable vtz used to be vz_max-vz_bar for cigars which use flat f(vz).
zimax           real    [m]   /0./
   # Maximum initial z of beam
zimin           real    [m]   /0./
   # Minimum initial z of beam
vzperlam        real    [m]   /0./
   # Wave length of sinusoidal perturbation in vz
vzperamp        real    [1]   /0./
   # Amplitude of sinusoidal perturbation in vz
vzperphs        real    [rad] /0./
   # Phase of sinusoidal perturbation in vz
prwall          real    [m]   /0./
   # Optional radius of cylindrical wall that absorbs particles
   # When zero, uses largest cylinder that fits in grid
prwallx         real    [m]   /0./
   # X location of the center of the cylindrical wall that absorbs particles
prwally         real    [m]   /0./
   # Y location of the center of the cylindrical wall that absorbs particles
np_pert         integer [1]   /0/
   # Number of particles in a perturbation to the main distribution.
pertz_den       real    [1]   /0./
   # Ratio of peak perturbation density to peak density of main distribution.
pertz_ctr       real    [m]   /0./
   # Center of perturbation.
pertz_len       real    [m]   /0.2/
   # Length of perturbation.
pertz_strgt     real    [m]   /0.4/
   # Fraction of perturbation in straight section beween parabolic end caps.
randoffset      integer [1]   /0/
   # Offset added to quiet start random seeds.

*********** InjectVars dump:
inject    integer    /0/   # Type of injection, (0: turned off,
                           # 1: constant current,
                           # 2: space-charge limited (Child-Langmuir),
                           # 3: space-charge limited (Gauss's law))
inj_param real       /1./  # Relaxation parameter for inject.  Mainly used
                           # for Egun iterative mode - set to 1 for time
                           # dependent injection, 0 for steady-state injection.
injpid    integer    /0/   # pid index for injection information for particles
ninject   integer    /1/   # Number of injection sources
leninjct  real       /0./  # Length of region into which partcls are injected
zinject(ninject)  _real /0./  # Z Start of injection in lab frame
rinject(ninject)  _real /LARGEPOS/ # Radius of source curvature
ainject(ninject)  _real    # Width of injection in x
binject(ninject)  _real    # Width of injection in y
ainjmin(ninject)  _real    # Minimum of injection in x (hole in source)
binjmin(ninject)  _real    # Minimum of injection in y (hole in source)
apinject(ninject) _real    # Convergence angle of injection in x
bpinject(ninject) _real    # Convergence angle of injection in y
xinject(ninject)  _real    # X location of injection source
yinject(ninject)  _real    # Y location of injection source
xpinject(ninject) _real    # X angle of injection source
ypinject(ninject) _real    # Y angle of injection source
vzinject(ninject) _real /0./ [m/s] # Starting velocity.
                                   # For inject == 1, autoset to vbeam
finject(ninject,ns) _real  # Species fraction for each source
inj_zstart(ninject) _real /0./ [m] # Starting location relative to the emitting
                                   # surface location.
inj_d(ninject)      _real /1./ # Distance from surface where phi is fetched.
                               # In units of dz.
inj_f(ninject)      _real /1./ # Scaling factor on the number of particles
                               # injected.
inj_dtemp(ninject)  _real /0./ # Distance from surface where temperature is added.
                               # In units of dz.
lvinject logical /.false./ # Sets whether injection source is included
                           # in field solve
vinject(ninject)  _real    # Voltage on the injection source
npinject     integer /0/   # Number of particles injected each step
npinje_s(ns) _integer      # Number of particles injected each step
                           # for each species
npinjtmp(ns,ninject) _integer # Temporary for saving actual number of
                           # particles injected for inject=1. Only meaningful
                           # in parallel version.
injctspc  integer    /0/   # Extra space added to particle arrays
injctcnt  integer          # Count for random number generators
jmininj(ninject) _real /0/        # Minimum current density emited from the
                                  # source.
jmaxinj(ninject) _real /LARGEPOS/ # Maximum current density emittable from
                                  # the source.
inj_nsmooth integer /0/    # Number of smoothing iterations to do on the
                           # E field in front of the emitting surface.
linj_spherical logical /.true./ # Flags whether or not to include the
                           # spherical correction term in the Child-Langmuir
                           #current.
linj_enormcl logical /.true./ # Flags whether or not the normal E field
                           # near the emitting surface is varies as z**(1/3)
                           # as in the Child-Langmuir solution or is uniform.
linj_eperp logical /.false./ # Flags whether or not to include the tangential
                           # electric field near the emitting surface.
linj_rectangle logical /.false./ # Flags whether injection is over a rectangular region and not 
# elliptical  

ntinj          integer /0/ # Number of transverse emitting surfaces.
nztinjmn(ntinj)   _integer # Minimum z extent of transverse injection.
nztinjmx(ntinj)   _integer # Maximum z extent of transverse injection.
nztmax         integer /0/ # Maximum length (in z grid cells) of
                           # transverse injection surfaces
rtinject(0:nztmax,ntinj) _real
                           # Radius as a function of z of transverse
                           # emitting surface.
nttinjmax      integer /0/ # Max number of trasnversely injected particles
nttinj(ntinj)     _integer # Number of azimuthal points that transverse
                           # injection is broken down into.
                           # Automatically calculated.
vtinject(ntinj)    _real   # Voltage on transverse emitting surface.
ftinject(ns,ntinj) _real   # Species fractions for transverse emitting surface
xtinject(ntinj)    _real   # X location of center of transverse emitting surface
ytinject(ntinj)    _real   # Y location of center of transverse emitting surface
vztinject(ntinj) _real /0./ [m/s] # Starting velocity.
                                  # For inject == 1, autoset to vbeam
ntinject(ns)      _integer # Number of particles injected off of the
                           # transverse injection sources.
jmaxtinj(ntinj)    _real   # Maximum injected current for transverse injection.
tinjprev(0:nttinjmax-1,0:nztmax,ntinj) _real
   # Grid holding number of particles injected on transverse surface
   # on previous time step.  Only used when inject=2


*********** OutParams dump:
# Various output quantities associated with the initial particle load 
currdens          real [amps/m**2] # Current density at center of beam
chrgdens          real [Coul/m**3] # Charge density at center of beam
numbdens          real [m**-3]     # Number density of real particles
omegap            real [s**-1]     # Plasma frequency
omegapdt          real [1]         # Plasma frequency times dt
omegaptq          real [1]         # Plasma frequency times quad period
taup              real [s]         # Plasma period (2pi/omegap)
vthx              real [m/s]       # Transverse X-thermal velocity
vthy              real [m/s]       # Transverse Y-thermal velocity
vthxdt            real [m/s]       # Transverse X-thermal velocity times dt
vthydt            real [m/s]       # Transverse Y-thermal velocity times dt
vthxdtodx         real [m/s]       # Transverse X-thermal velocity times dt/dx
vthydtody         real [m/s]       # Transverse Y-thermal velocity times dt/dy
lamdebx           real [m]         # Transverse X-Debye wavelength
lamdeby           real [m]         # Transverse Y-Debye wavelength
lamdebxodx        real [1]         # Transverse X-Debye wavelength over dx
lamdebyody        real [1]         # Transverse Y-Debye wavelength over dy
vthzdt            real [m]         # Longitudinal thermal velocity times dt
vthzdtodz         real [1]         # Longitudinal thermal velocity times dt/dz
lamdebz           real [m]         # Longitudinal Debye wavelength
lamdebzodz        real [1]         # Longitudinal Debye wavelength over dz
vbeamoc           real [1]         # Normalized beam velocity
npreal            real [1]         # Number of real particles
totmass           real [m]         # Total mass
totchrg           real [Coul]      # Total charge
genperv           real [1]         # Generalized perveance
charcurr          real [amps]      # Characteristic current
budker            real [1]         # Budker parameter
vwave             real [m/s]       # Space charge wave velocity
femtxofscx        real [1]         # Ratio X-emittance to X-space charge forces
femtyofscy        real [1]         # Ratio X-emittance to X-space charge forces
sigmax            real [degrees]   # Depressed   particle X-phase advance
sigma0x           real [degrees]   # Undepressed particle X-phase advance 
sigmay            real [degrees]   # Depressed   particle Y-phase advance
sigma0y           real [degrees]   # Undepressed particle Y-phase advance 
omegabx           real [s**-1]     # Depressed   X-betatron frequency
omegab0x          real [s**-1]     # Undepressed X-betatron frequency
omegaby           real [s**-1]     # Depressed   Y-betatron frequency
omegab0y          real [s**-1]     # Undepressed Y-betatron frequency
taubx             real [s]         # Depressed   X-betatron period
taub0x            real [s]         # Undepressed X-betatron period
tauby             real [s]         # Depressed   Y-betatron period
taub0y            real [s]         # Undepressed Y-betatron period
lambdabx          real [m]         # Depressed   X-betatron wavelength
lambdab0x         real [m]         # Undepressed X-betatron wavelength
lambdaby          real [m]         # Depressed   Y-betatron wavelength
lambdab0y         real [m]         # Undepressed Y-betatron wavelength

*********** Z_arrays dump:
# 1d arrays used by >1 pkg
# nzzarr is usually set to nz by pkg w3d (etc.) at generation
dzz                     real   [m]     # Z_arrays mesh grid cell size
dzzi                    real   [m]     # Z_arrays mesh grid cell size inverse
zzmin                   real   [m]     # Z_arrays mesh maximum in z
zzmax                   real   [m]     # Z_arrays mesh minimum in z
nzzarr            integer /0/  [1]     # Length of arrays in group Z_arrays
zplmesh(0:nzzarr)      _real   [m]     # Z mesh used for qtys in group Z_arrays
curr(0:nzzarr)         _real   [A]     # Beam current
egap(0:nzzarr)         _real   [V/m]   # Gap electric field (smeared in z)
linechg(0:nzzarr)      _real   [C/m]   # Line charge density
vzofz(0:nzzarr)        _real   [m/s]   # Mean axial speed vs z
rhoax(0:nzzarr)        _real   [C/m^3] # charge density on axis
phiax(0:nzzarr)        _real   [V]     # potential on axis
ezax(0:nzzarr)         _real   [V/m]   # space charge E field on axis
eearsofz(0:nzzarr)     _real   [V/m]   # confining Eears, as a function of z
prwallz(0:nzzarr)      _real   [m]     # Radius at which particles are absorbed
prwallxz(0:nzzarr)     _real   [m]     # X of center of cylindrical wall
prwallyz(0:nzzarr)     _real   [m]     # Y of center of cylindrical wall
prwelips(0:nzzarr)     _real   [1]     # ellipticity of cylindrical wall
                                       # (ay = prwelips*ax)
lostpars(0:nzzarr)  _integer           # number of lost particles by zcells
lamkreal(0:nzzarr)     _real   [C/m]   # Real part of FFT of lambda
lamkimag(0:nzzarr)     _real   [C/m]   # Imaginary part of FFT of lambda


xmaxz(0:nzzarr)        _real   [m]   /+LARGEPOS/ # z-dependent locations used 
# to remove particles outside a rectangular region. 
xminz(0:nzzarr)        _real   [m]   /-LARGEPOS/ # z-dependent locations used
# to remove particles outside a rectangular region.
ymaxz(0:nzzarr)        _real   [m]   /+LARGEPOS/ # z-dependent locations used
# to remove particles outside a rectangular region.
yminz(0:nzzarr)        _real   [m]   /-LARGEPOS/ # z-dependent locations used 
# to remove particles outside a rectangular region.

*********** Io dump:
# Quantities associated with input / output
warpout                   Filedes  /-1/ -dump
   # Output unit number
nfiche                    integer  /-1/
   # Number of the current fiche
verbosity                 integer /2/
   # Specifies the amount of output during a simulation.  Currently only
   # affects whether the one line diagnostic is printed on special steps.
   # Line is not printed when less than 2.

*********** Win_Moments dump:
# Particle and field moment data (including emittances) at current timestep
nzwind       integer # Actual number of z windows used for emittance calcs
pnum(0:nzwind)          _real [1]     # Total no. of (physical) ions in window
xbar(0:nzwind)          _real [m]     # Mean X coordinate in window
ybar(0:nzwind)          _real [m]     # Mean Y coordinate in window
zbar(0:nzwind)          _real [m]     # Mean axial location in window
xpbar(0:nzwind)         _real [1]     # Mean X' in window
ypbar(0:nzwind)         _real [1]     # Mean Y' in window
vxbar(0:nzwind)         _real [m/s]   # Mean Vx in window
vybar(0:nzwind)         _real [m/s]   # Mean Vy in window
vzbar(0:nzwind)         _real [m/s]   # Mean Vz in window
xybar(0:nzwind)         _real [m**2]  # Mean product of X and Y in window
xypbar(0:nzwind)        _real [m]     # Mean product of X  and Y'
yxpbar(0:nzwind)        _real [m]     # Mean product of Y  and X'
xpypbar(0:nzwind)       _real [1]     # Mean product of X' and Y'
xsqbar(0:nzwind)        _real [m**2]  # Mean X-squared in window
ysqbar(0:nzwind)        _real [m**2]  # Mean Y-squared in window
zsqbar(0:nzwind)        _real [m**2]  # Mean Z-squared in window
xpsqbar(0:nzwind)       _real [1]     # Mean X' squared in window
ypsqbar(0:nzwind)       _real [1]     # Mean Y' squared in window
vxsqbar(0:nzwind)       _real [m/s]   # Mean Vx squared in window
vysqbar(0:nzwind)       _real [m/s]   # Mean Vy squared in window
vzsqbar(0:nzwind)       _real [m/s]   # Mean Vz squared in window
xxpbar(0:nzwind)        _real [m]     # Mean product of X and X' in window
yypbar(0:nzwind)        _real [m]     # Mean product of Y and Y' in window
zvzbar(0:nzwind)        _real [m]     # Mean product of Z and Vz in window
xvzbar(0:nzwind)        _real [m]     # Mean product of X and Vz in window
yvzbar(0:nzwind)        _real [m]     # Mean product of X and Vz in window
vxvzbar(0:nzwind)       _real [m]     # Mean product of Vx and Vz in window
vyvzbar(0:nzwind)       _real [m]     # Mean product of Vy and Vz in window
xrms(0:nzwind)          _real [m]     # RMS X in window
yrms(0:nzwind)          _real [m]     # RMS Y in window
zrms(0:nzwind)          _real [m]     # RMS Z in window
rrms(0:nzwind)          _real [m]     # RMS R in window
xprms(0:nzwind)         _real [m]     # RMS X' in window
yprms(0:nzwind)         _real [m]     # RMS Y' in window
epsx(0:nzwind)          _real [m-rad] # X-X' emittance
epsy(0:nzwind)          _real [m-rad] # Y-Y' emittance
epsz(0:nzwind)          _real [m-rad] # Z-Z' emittance
epsnx(0:nzwind)         _real [mm-mrad] # X-X' normalized emittance
epsny(0:nzwind)         _real [mm-mrad] # Y-Y' normalized emittance
epsnz(0:nzwind)         _real [mm-mrad] # Z-Z' normalized emittance
epsg(0:nzwind)          _real [m-rad] # Generalized emittance
epsh(0:nzwind)          _real [m-rad] # Generalized emittance
epsng(0:nzwind)         _real [mm-mrad] # Generalized normalized emittance
epsnh(0:nzwind)         _real [mm-mrad] # Generalized normalized emittance
vxrms(0:nzwind)         _real [m/s]   # True RMS Vx in window
vyrms(0:nzwind)         _real [m/s]   # True RMS Vy in window
vzrms(0:nzwind)         _real [m/s]   # True RMS Vz in window
rhomid(0:nzwind)        _real [C/m^3] # Charge dens. on axis at ctr of window
rhomax(0:nzwind)        _real [C/m^3] # Charge dens. max-over-X,Y at ctr of win.

*********** Z_Moments dump:
# Particle and field moment data (including emittances) at current timestep
# as a function of Z 
zmmntmax             real         # Moments grid maximum in Z
zmmntmin             real         # Moments grid minimum in Z
nzmmnt               integer /0/  # Number of points in z moments grid
dzm                  real         # Moments grid cell size
dzmi                 real         # Moments grid cell size inverse
numzmmnt             integer /NUMZMMNT/ # Number of moments calculated
zmntmesh(0:nzmmnt)  _real [m]     # Z mesh associated with Z moments
pnumz(0:nzmmnt)     _real [1]     # No. of (physical) ions at grid point
xbarz(0:nzmmnt)     _real [m]     # Mean X coordinate at grid point
ybarz(0:nzmmnt)     _real [m]     # Mean Y coordinate at grid point
zbarz(0:nzmmnt)     _real [m]     # Mean axial location at grid point
xpbarz(0:nzmmnt)    _real [1]     # Mean X' at grid point
ypbarz(0:nzmmnt)    _real [1]     # Mean Y' at grid point
vxbarz(0:nzmmnt)    _real [m/s]   # Mean Vx at grid point
vybarz(0:nzmmnt)    _real [m/s]   # Mean Vy at grid point
vzbarz(0:nzmmnt)    _real [m/s]   # Mean Vz at grid point
xybarz(0:nzmmnt)    _real [m**2]  # Mean product of X  and Y  at grid point
xypbarz(0:nzmmnt)   _real [m]     # Mean product of X  and Y' at grid point
yxpbarz(0:nzmmnt)   _real [m]     # Mean product of Y  and X' at grid point
xpypbarz(0:nzmmnt)  _real [1]     # Mean product of X' and Y' at grid point 
xsqbarz(0:nzmmnt)   _real [m**2]  # Mean X-squared at grid point
ysqbarz(0:nzmmnt)   _real [m**2]  # Mean Y-squared at grid point
zsqbarz(0:nzmmnt)   _real [m**2]  # Mean Z-squared at grid point
xpsqbarz(0:nzmmnt)  _real [1]     # Mean X' squared at grid point
ypsqbarz(0:nzmmnt)  _real [1]     # Mean Y' squared at grid point
vxsqbarz(0:nzmmnt)  _real [m/s]   # Mean Vx squared at grid point
vysqbarz(0:nzmmnt)  _real [m/s]   # Mean Vy squared at grid point
vzsqbarz(0:nzmmnt)  _real [m/s]   # Mean Vz squared at grid point
xxpbarz(0:nzmmnt)   _real [m]     # Mean product of X and X' at grid point
yypbarz(0:nzmmnt)   _real [m]     # Mean product of Y and Y' at grid point
zvzbarz(0:nzmmnt)   _real [m]     # Mean product of Z and Vz at grid point
xvzbarz(0:nzmmnt)   _real [m]     # Mean product of X and Vz at grid point
yvzbarz(0:nzmmnt)   _real [m]     # Mean product of X and Vz at grid point
vxvzbarz(0:nzmmnt)  _real [m]     # Mean product of Vx and Vz at grid point
vyvzbarz(0:nzmmnt)  _real [m]     # Mean product of Vy and Vz at grid point
xrmsz(0:nzmmnt)     _real [m]     # RMS X at grid point
yrmsz(0:nzmmnt)     _real [m]     # RMS Y at grid point
zrmsz(0:nzmmnt)     _real [m]     # RMS Z at grid point
rrmsz(0:nzmmnt)     _real [m]     # RMS R at grid point
xprmsz(0:nzmmnt)    _real [m]     # RMS X' at grid point
yprmsz(0:nzmmnt)    _real [m]     # RMS Y' at grid point
epsxz(0:nzmmnt)     _real [m-rad] # X-X' emittance at grid point
epsyz(0:nzmmnt)     _real [m-rad] # Y-Y' emittance at grid point
epszz(0:nzmmnt)     _real [m-rad] # Z-Z' emittance at grid point
epsnxz(0:nzmmnt)    _real [mm-mrad] # X-X' normalized emittance at grid point
epsnyz(0:nzmmnt)    _real [mm-mrad] # Y-Y' normalized emittance at grid point
epsnzz(0:nzmmnt)    _real [mm-mrad] # Z-Z' normalized emittance at grid point
epsgz(0:nzmmnt)     _real [m-rad]   # Generalized emittance on grid
epshz(0:nzmmnt)     _real [m-rad]   # Generalized emittance on grid
epsngz(0:nzmmnt)    _real [mm-mrad] # Generalized normalized emittance on grid
epsnhz(0:nzmmnt)    _real [mm-mrad] # Generalized normalized emittance on grid
vxrmsz(0:nzmmnt)    _real [m/s]   # True RMS Vx at grid point
vyrmsz(0:nzmmnt)    _real [m/s]   # True RMS Vy at grid point
vzrmsz(0:nzmmnt)    _real [m/s]   # True RMS Vz at grid point
rhomidz(0:nzmmnt)   _real [C/m^3] # Charge dens. on axis at grid point
rhomaxz(0:nzmmnt)   _real [C/m^3] # Charge dens. max-over-X,Y at grid point
tempmaxp(6)                    real # Temporary work array
tempminp(6)                    real # Temporary work array
tempzmmnts0(NUMZMMNT)          real # Temporary work array
tempzmmnts(0:nzmmnt,NUMZMMNT) _real # Temporary work array

********** Lab_Moments dump:
# Particle moment data as a function of time at locations in the lab frame
nlabwn  integer   /50/ # number of lab windows
zlw(nlabwn) _real [m] /0./ # z for lab windows
iflabwn integer /1/ # turns on lab window moments (0 off; 1 on)
itlabwn integer /0/ # Sets how often the lab moments are calculated
ntlabwn integer     # Maximum number of times lab frame moments are calculated
ilabwn(nlabwn) _integer # Number of times lab frame moments have been calculated
timelw(ntlabwn,nlabwn)  _real # Time in lab frame
pnumlw(ntlabwn,nlabwn)  _real # Number of particles in lab frame
xbarlw(ntlabwn,nlabwn)  _real # X bar in lab frame
ybarlw(ntlabwn,nlabwn)  _real # Y bar in lab frame
vzbarlw(ntlabwn,nlabwn) _real # Vz bar in lab frame
epsxlw(ntlabwn,nlabwn)  _real # X emittance in lab frame
epsylw(ntlabwn,nlabwn)  _real # Y emittance in lab frame
epszlw(ntlabwn,nlabwn)  _real # Z emittance in lab frame
vxrmslw(ntlabwn,nlabwn) _real # Vx RMS in lab frame
vyrmslw(ntlabwn,nlabwn) _real # Vy RMS in lab frame
vzrmslw(ntlabwn,nlabwn) _real # Vz RMS in lab frame
xrmslw(ntlabwn,nlabwn)  _real # X RMS in lab frame
yrmslw(ntlabwn,nlabwn)  _real # Y RMS in lab frame
rrmslw(ntlabwn,nlabwn)  _real # R RMS in lab frame
currlw(ntlabwn,nlabwn)  _real # Current in lab frame
linechglw(ntlabwn,nlabwn)  _real # Line-charge in lab frame
lostparslw(ntlabwn,nlabwn)  _real # Number of lost particles in lab frame

*********** Moments dump:
# Scalar moments of general interest
ek      real [J]      # Total Kinetic energy
ekzmbe  real [J]      # Total Z Kinetic energy minus beam energy
                      # 1/2m(vzsqbar - vbeam**2)
ekzbeam real [J]      # Z Kinetic energy in the beam frame
                      # 1/2m ave[(vz - vbeam)**2]
ekperp  real [J]      # Perp Kinetic energy
ese     real [J]      # Electrostatic energy
pz      real [kg-m/s] # Total axial momentum (subtracting out Vbeam)
xmaxp   real [m]      # Maximum X  over particles (set intermittently)
xminp   real [m]      # Minimum X  over particles (set intermittently)
ymaxp   real [m]      # Maximum Y  over particles (set intermittently)
yminp   real [m]      # Minimum Y  over particles (set intermittently)
zmaxp   real [m]      # Maximum Z  over particles (set intermittently)
zminp   real [m]      # Minimum Z  over particles (set intermittently)
vxmaxp  real [m/s]    # Maximum Vx over particles (set intermittently)
vxminp  real [m/s]    # Minimum Vx over particles (set intermittently)
vymaxp  real [m/s]    # Maximum Vy over particles (set intermittently)
vyminp  real [m/s]    # Minimum Vy over particles (set intermittently)
vzmaxp  real [m/s]    # Maximum Vz over particles (set intermittently)
vzminp  real [m/s]    # Minimum Vz over particles (set intermittently)

*********** Hist dump history:
# History data
jhist                          integer  /-1/
   # pointer to current entry in history arrays
lenhist                        integer  /0/
   # length of the history arrays
thist(0:lenhist)              _real [s]     limited (0:jhist)
   # Times at which data is saved
hzbeam(0:lenhist)             _real [m]     limited (0:jhist)
   # Beam frame location in lab frame
hvbeam(0:lenhist)             _real [m/s]   limited (0:jhist)
   # Beam frame velocity
hbmlen(0:lenhist)             _real [m]     limited (0:jhist)
   # RMS beam length
hefld(0:lenhist)              _real [J]     limited (0:jhist)
   # Field energy
hekzmbe(0:lenhist)            _real [J]     limited (0:jhist)
   # Total Z Kinetic energy minus beam energy
hekzbeam(0:lenhist)           _real [J]     limited (0:jhist)
   # Z Kinetic energy in the beam frame
hekperp(0:lenhist)            _real [J]     limited (0:jhist)
   # Perp Kinetic energy
hxmaxp(0:lenhist)             _real [m]     limited (0:jhist)
   # History of maximum X over particles
hxminp(0:lenhist)             _real [m]     limited (0:jhist)
   # History of minimum X over particles
hymaxp(0:lenhist)             _real [m]     limited (0:jhist)
   # History of maximum Y over particles
hyminp(0:lenhist)             _real [m]     limited (0:jhist)
   # History of minimum Y over particles
hzmaxp(0:lenhist)             _real [m]     limited (0:jhist)
   # History of maximum Z over particles
hzminp(0:lenhist)             _real [m]     limited (0:jhist)
   # History of minimum Z over particles
hvxmaxp(0:lenhist)            _real [m/s]   limited (0:jhist)
   # History of maximum Vx over particles
hvxminp(0:lenhist)            _real [m/s]   limited (0:jhist)
   # History of minimum Vx over particles
hvymaxp(0:lenhist)            _real [m/s]   limited (0:jhist)
   # History of maximum Vy over particles
hvyminp(0:lenhist)            _real [m/s]   limited (0:jhist)
   # History of minimum Vy over particles
hvzmaxp(0:lenhist)            _real [m/s]   limited (0:jhist)
   # History of maximum Vz over particles
hvzminp(0:lenhist)            _real [m/s]   limited (0:jhist)
   # History of minimum Vz over particles
hepsx(0:nzwind,0:lenhist)     _real [m-r]   limited (0:nzwind,0:jhist) +winhist
   # X-X' emittance by window as a function of time
hepsy(0:nzwind,0:lenhist)     _real [m-r]   limited (0:nzwind,0:jhist) +winhist
   # Y-Y' emittance by window as a function of time
hepsz(0:nzwind,0:lenhist)     _real [m-r]   limited (0:nzwind,0:jhist) +winhist
   # Z-Z' emittance by window as a function of time
hepsnx(0:nzwind,0:lenhist)    _real [mm-mr] limited (0:nzwind,0:jhist) +winhist
   # X-X' normalized emittance by window as a function of time
hepsny(0:nzwind,0:lenhist)    _real [mm-mr] limited (0:nzwind,0:jhist) +winhist
   # Y-Y' normalized emittance by window as a function of time
hepsnz(0:nzwind,0:lenhist)    _real [mm-mr] limited (0:nzwind,0:jhist) +winhist
   # Z-Z' normalized emittance by window as a function of time
hepsg(0:nzwind,0:lenhist)     _real [m-r]   limited (0:nzwind,0:jhist) +winhist
   # Generalized emittance by window as a function of time
hepsh(0:nzwind,0:lenhist)     _real [m-r]   limited (0:nzwind,0:jhist) +winhist
   # Generalized emittance by window as a function of time
hepsng(0:nzwind,0:lenhist)    _real [mm-mr] limited (0:nzwind,0:jhist) +winhist
   # Generalized normalized emittance by window as a function of time
hepsnh(0:nzwind,0:lenhist)    _real [mm-mr] limited (0:nzwind,0:jhist) +winhist
   # Generalized normalized emittance by window as a function of time
hpnum(0:nzwind,0:lenhist)     _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Number of particles in each z window (species 1)
hrhomid(0:nzwind,0:lenhist)   _real [C/m^3] limited (0:nzwind,0:jhist) +winhist
   # Charge density on axis at center of z window as a fcn of time
hrhomax(0:nzwind,0:lenhist)   _real [C/m^3] limited (0:nzwind,0:jhist) +winhist
   # Charge dens. max-over-x,y at ctr of z window as a fcn of time
hxbar(0:nzwind,0:lenhist)     _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True mean x by window as a function of time
hybar(0:nzwind,0:lenhist)     _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True mean y by window as a function of time
hxybar(0:nzwind,0:lenhist)    _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True mean xy by window as a function of time
hxrms(0:nzwind,0:lenhist)     _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True RMS x by window as a function of time
hyrms(0:nzwind,0:lenhist)     _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True RMS y by window as a function of time
hrrms(0:nzwind,0:lenhist)     _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True RMS r by window as a function of time
hxprms(0:nzwind,0:lenhist)    _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True RMS x' by window as a function of time
hyprms(0:nzwind,0:lenhist)    _real [m]     limited (0:nzwind,0:jhist) +winhist
   # True RMS y' by window as a function of time
hxsqbar(0:nzwind,0:lenhist)   _real [m^2]   limited (0:nzwind,0:jhist) +winhist
   # x squared bar by window as a function of time
hysqbar(0:nzwind,0:lenhist)   _real [m^2]   limited (0:nzwind,0:jhist) +winhist
   # y squared bar by window as a function of time
hvxbar(0:nzwind,0:lenhist)    _real [m/s]   limited (0:nzwind,0:jhist) +winhist
   # Mean vx by window as a function of time
hvybar(0:nzwind,0:lenhist)    _real [m/s]   limited (0:nzwind,0:jhist) +winhist
   # Mean vy by window as a function of time
hvzbar(0:nzwind,0:lenhist)    _real [m/s]   limited (0:nzwind,0:jhist) +winhist
   # Mean vz by window as a function of time
hxpbar(0:nzwind,0:lenhist)    _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean x' by window as a function of time
hypbar(0:nzwind,0:lenhist)    _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean y' by window as a function of time
hvxrms(0:nzwind,0:lenhist)    _real [m/s]   limited (0:nzwind,0:jhist) +winhist
   # True RMS vx by window as a function of time
hvyrms(0:nzwind,0:lenhist)    _real [m/s]   limited (0:nzwind,0:jhist) +winhist
   # True RMS vy by window as a function of time
hvzrms(0:nzwind,0:lenhist)    _real [m/s]   limited (0:nzwind,0:jhist) +winhist
   # True RMS vz by window as a function of time
hxpsqbar(0:nzwind,0:lenhist)  _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean x' squared by window as a function of time
hypsqbar(0:nzwind,0:lenhist)  _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean y' squared by window as a function of time
hxxpbar(0:nzwind,0:lenhist)   _real [m]     limited (0:nzwind,0:jhist) +winhist
   # Mean x * x' by window as a function of time
hyypbar(0:nzwind,0:lenhist)   _real [m]     limited (0:nzwind,0:jhist) +winhist
   # Mean y  * y' by window as a function of time
hxypbar(0:nzwind,0:lenhist)   _real [m]     limited (0:nzwind,0:jhist) +winhist
   # Mean x  * y' by window as a function of time 
hyxpbar(0:nzwind,0:lenhist)   _real [m]     limited (0:nzwind,0:jhist) +winhist
   # Mean y  * x' by window as a function of time 
hxpypbar(0:nzwind,0:lenhist)  _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean x' * y' by window as a function of time
hxvzbar(0:nzwind,0:lenhist)   _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean x * vz by window as a function of time
hyvzbar(0:nzwind,0:lenhist)   _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean y * vz by window as a function of time
hvxvzbar(0:nzwind,0:lenhist)  _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean vx * vz by window as a function of time
hvyvzbar(0:nzwind,0:lenhist)  _real [1]     limited (0:nzwind,0:jhist) +winhist
   # Mean vy * vz by window as a function of time
lhlinechg logical /.true./   # Turns on history of line charge
ihlinechg integer /1/        # Multiplier for hlinechg memory size (autoset)
hlinechg(0:nzzarr*ihlinechg,0:lenhist) _real [C/m] limited (0:nzzarr,0:jhist)
            +zhist           # Line charge density vs. space and time
lhvzofz logical /.true./     # Turns on history of vz
ihvzofz integer /1/          # Multiplier for hvzofz memory size (autoset)
hvzofz(0:nzzarr*ihvzofz,0:lenhist)  _real [m/s] limited (0:nzzarr,0:jhist)
            +zhist           # Vz versus space and time
lhcurrz logical /.false./    # Turns on history of current
ihcurrz integer /1/          # Multiplier for hcurrz memory size (autoset)
hcurrz(0:nzzarr*ihcurrz,0:lenhist)  _real [m/s] limited (0:nzzarr,0:jhist)
            +zhist           # Current versus space and time
lhepsxz logical /.false./    # Turns on history of X emittance
ihepsxz integer /0 /         # Multiplier for hepsxz memory size (autoset)
hepsxz(0:nzmmnt*ihepsxz,0:lenhist)  _real [m-r] limited (0:nzmmnt,0:jhist)
            +zhist           # X emittance versus space and time
lhepsyz logical /.false./    # Turns on history of Y emittance
ihepsyz integer /0 /         # Multiplier for hepsyz memory size (autoset)
hepsyz(0:nzmmnt*ihepsyz,0:lenhist)  _real [m-r] limited (0:nzmmnt,0:jhist)
            +zhist           # Y emittance versus space and time
lhepsnxz logical /.false./   # Turns on history of X normalized emittance
ihepsnxz integer /0 /        # Multiplier for hepsnxz memory size (autoset)
hepsnxz(0:nzmmnt*ihepsnxz,0:lenhist)  _real [mm-mrad] limited (0:nzmmnt,0:jhist)
            +zhist           # X normalized emittance versus space and time
lhepsnyz logical /.false./   # Turns on history of Y normalized emittance
ihepsnyz integer /0 /        # Multiplier for hepsnyz memory size (autoset)
hepsnyz(0:nzmmnt*ihepsnyz,0:lenhist)  _real [mm-mrad] limited (0:nzmmnt,0:jhist)
            +zhist           # Y normalized emittance versus space and time
lhepsgz logical /.false./    # Turns on history of Generalized emittance
ihepsgz integer /0 /         # Multiplier for hepsgz memory size (autoset)
hepsgz(0:nzmmnt*ihepsgz,0:lenhist)  _real [m-rad] limited (0:nzmmnt,0:jhist)
            +zhist           # Generalized emittance versus space and time
lhepshz logical /.false./    # Turns on history of Generalized emittance
ihepshz integer /0 /         # Multiplier for hepshz memory size (autoset)
hepshz(0:nzmmnt*ihepshz,0:lenhist)  _real [m-rad] limited (0:nzmmnt,0:jhist)
            +zhist           # Generalized emittance versus space and time
lhepsngz logical /.false./   # Turns on history of Generalized nrmlzd emittance
ihepsngz integer /0 /        # Multiplier for hepsngz memory size (autoset)
hepsngz(0:nzmmnt*ihepsngz,0:lenhist)  _real [mm-mrad] limited (0:nzmmnt,0:jhist)
            +zhist           # Generalized nrmlzd emittance versus space andtime
lhepsnhz logical /.false./   # Turns on history of Generalized nrmlzd emittance
ihepsnhz integer /0 /        # Multiplier for hepsnhz memory size (autoset)
hepsnhz(0:nzmmnt*ihepsnhz,0:lenhist)  _real [mm-mrad] limited (0:nzmmnt,0:jhist)
            +zhist           # Generalized nrmlzd emittance versus space andtime
lhxbarz logical /.false./    # Turns on history of X bar
ihxbarz integer /0 /         # Multiplier for hxbarz memory size (autoset)
hxbarz(0:nzmmnt*ihxbarz,0:lenhist)  _real [m] limited (0:nzmmnt,0:jhist)
            +zhist           # X bar versus space and time
lhybarz logical /.false./    # Turns on history of Y bar
ihybarz integer /0 /         # Multiplier for hybarz memory size (autoset)
hybarz(0:nzmmnt*ihybarz,0:lenhist)  _real [m] limited (0:nzmmnt,0:jhist)
            +zhist           # Y bar versus space and time
lhxybarz logical /.false./   # Turns on history of XY bar
ihxybarz integer /0 /        # Multiplier for hxybarz memory size (autoset)
hxybarz(0:nzmmnt*ihxybarz,0:lenhist)  _real [m**2] limited (0:nzmmnt,0:jhist)
            +zhist           # XY bar versus space and time
lhxrmsz logical /.false./    # Turns on history of X rms
ihxrmsz integer /0 /         # Multiplier for hxrmsz memory size (autoset)
hxrmsz(0:nzmmnt*ihxrmsz,0:lenhist)  _real [m] limited (0:nzmmnt,0:jhist)
            +zhist           # X rms versus space and time
lhyrmsz logical /.false./    # Turns on history of Y rms
ihyrmsz integer /0 /         # Multiplier for hyrmsz memory size (autoset)
hyrmsz(0:nzmmnt*ihyrmsz,0:lenhist)  _real [m] limited (0:nzmmnt,0:jhist)
            +zhist           # Y rms versus space and time
lhrrmsz logical /.false./    # Turns on history of X rms
ihrrmsz integer /0 /         # Multiplier for hrrmsz memory size (autoset)
hrrmsz(0:nzmmnt*ihrrmsz,0:lenhist)  _real [m] limited (0:nzmmnt,0:jhist)
            +zhist           # X rms versus space and time
lhxprmsz logical /.false./   # Turns on history of X' rms
ihxprmsz integer /0 /        # Multiplier for hxprmsz memory size (autoset)
hxprmsz(0:nzmmnt*ihxprmsz,0:lenhist)  _real [rad] limited (0:nzmmnt,0:jhist)
            +zhist           # X' rms versus space and time
lhyprmsz logical /.false./   # Turns on history of Y' rms
ihyprmsz integer /0 /        # Multiplier for hyprmsz memory size (autoset)
hyprmsz(0:nzmmnt*ihyprmsz,0:lenhist)  _real [rad] limited (0:nzmmnt,0:jhist)
            +zhist           # Y' rms versus space and time
lhxsqbarz logical /.false./  # Turns on history of X**2 bar
ihxsqbarz integer /0 /       # Multiplier for hxsqbarz memory size (autoset)
hxsqbarz(0:nzmmnt*ihxsqbarz,0:lenhist)  _real [m**2] limited (0:nzmmnt,0:jhist)
            +zhist           # X**2 bar versus space and time
lhysqbarz logical /.false./  # Turns on history of Y**2 bar
ihysqbarz integer /0 /       # Multiplier for hysqbarz memory size (autoset)
hysqbarz(0:nzmmnt*ihysqbarz,0:lenhist)  _real [m**2] limited (0:nzmmnt,0:jhist)
            +zhist           # Y**2 bar versus space and time
lhvxbarz logical /.false./   # Turns on history of Vx bar
ihvxbarz integer /0 /        # Multiplier for hvxbarz memory size (autoset)
hvxbarz(0:nzmmnt*ihvxbarz,0:lenhist)  _real [m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # Vx bar versus space and time
lhvybarz logical /.false./   # Turns on history of Vy bar
ihvybarz integer /0 /        # Multiplier for hvybarz memory size (autoset)
hvybarz(0:nzmmnt*ihvybarz,0:lenhist)  _real [m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # Vy bar versus space and time
lhvzbarz logical /.false./   # Turns on history of Vz bar
ihvzbarz integer /0 /        # Multiplier for hvzbarz memory size (autoset)
hvzbarz(0:nzmmnt*ihvzbarz,0:lenhist)  _real [m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # Vz bar versus space and time
lhxpbarz logical /.false./   # Turns on history of X' bar
ihxpbarz integer /0 /        # Multiplier for hxpbarz memory size (autoset)
hxpbarz(0:nzmmnt*ihxpbarz,0:lenhist)  _real [rad] limited (0:nzmmnt,0:jhist)
            +zhist           # X' bar versus space and time
lhypbarz logical /.false./   # Turns on history of Y' bar
ihypbarz integer /0 /        # Multiplier for hypbarz memory size (autoset)
hypbarz(0:nzmmnt*ihypbarz,0:lenhist)  _real [rad] limited (0:nzmmnt,0:jhist)
            +zhist           # Y' bar versus space and time
lhvxrmsz logical /.false./   # Turns on history of Vx rms
ihvxrmsz integer /0 /        # Multiplier for hvxrmsz memory size (autoset)
hvxrmsz(0:nzmmnt*ihvxrmsz,0:lenhist)  _real [m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # Vx rms versus space and time
lhvyrmsz logical /.false./   # Turns on history of Vy rms
ihvyrmsz integer /0 /        # Multiplier for hvyrmsz memory size (autoset)
hvyrmsz(0:nzmmnt*ihvyrmsz,0:lenhist)  _real [m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # Vy rms versus space and time
lhvzrmsz logical /.false./   # Turns on history of Vz rms
ihvzrmsz integer /0 /        # Multiplier for hvzrmsz memory size (autoset)
hvzrmsz(0:nzmmnt*ihvzrmsz,0:lenhist)  _real [m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # Vz rms versus space and time
lhxpsqbarz logical /.false./ # Turns on history of X'**2 bar
ihxpsqbarz integer /0 /      # Multiplier for hxpsqbarz memory size (autoset)
hxpsqbarz(0:nzmmnt*ihxpsqbarz,0:lenhist)  _real [rad**2] limited (0:nzmmnt,0:jhist)
            +zhist           # X'**2 bar versus space and time
lhypsqbarz logical /.false./ # Turns on history of Y'**2 bar
ihypsqbarz integer /0 /      # Multiplier for hypsqbarz memory size (autoset)
hypsqbarz(0:nzmmnt*ihypsqbarz,0:lenhist)  _real [rad**2] limited (0:nzmmnt,0:jhist)
            +zhist           # Y'**2 bar versus space and time
lhxxpbarz logical /.false./  # Turns on history of XX' bar
ihxxpbarz integer /0 /       # Multiplier for hxxpbarz memory size (autoset)
hxxpbarz(0:nzmmnt*ihxxpbarz,0:lenhist)  _real [m-rad] limited (0:nzmmnt,0:jhist)
            +zhist           # XX' bar versus space and time
lhyypbarz logical /.false./  # Turns on history of YY' bar
ihyypbarz integer /0 /       # Multiplier for hyypbarz memory size (autoset)
hyypbarz(0:nzmmnt*ihyypbarz,0:lenhist)  _real [m-rad] limited (0:nzmmnt,0:jhist)
            +zhist           # YY' bar versus space and time
lhxypbarz logical /.false./  # Turns on history of XY' bar
ihxypbarz integer /0 /       # Multiplier for hxypbarz memory size (autoset)
hxypbarz(0:nzmmnt*ihxypbarz,0:lenhist)  _real [m-rad] limited (0:nzmmnt,0:jhist)
            +zhist           # XY' bar versus space and time
lhyxpbarz logical /.false./  # Turns on history of YX' bar
ihyxpbarz integer /0 /       # Multiplier for hyxpbarz memory size (autoset)
hyxpbarz(0:nzmmnt*ihyxpbarz,0:lenhist)  _real [m-rad] limited (0:nzmmnt,0:jhist)
            +zhist           # YX' bar versus space and time
lhxpypbarz logical /.false./ # Turns on history of X'Y' bar
ihxpypbarz integer /0 /      # Multiplier for hxpypbarz memory size (autoset)
hxpypbarz(0:nzmmnt*ihxpypbarz,0:lenhist)  _real [rad**2] limited (0:nzmmnt,0:jhist)
            +zhist           # X'Y' bar versus space and time
lhxvzbarz logical /.false./  # Turns on history of XVz bar
ihxvzbarz integer /0 /       # Multiplier for hxvzbarz memory size (autoset)
hxvzbarz(0:nzmmnt*ihxvzbarz,0:lenhist)  _real [m*m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # XVz bar versus space and time
lhyvzbarz logical /.false./  # Turns on history of YVz bar
ihyvzbarz integer /0 /       # Multiplier for hyvzbarz memory size (autoset)
hyvzbarz(0:nzmmnt*ihyvzbarz,0:lenhist)  _real [m*m/s] limited (0:nzmmnt,0:jhist)
            +zhist           # YVz bar versus space and time
lhvxvzbarz logical /.false./ # Turns on history of VxVz bar
ihvxvzbarz integer /0 /      # Multiplier for hvxvzbarz memory size (autoset)
hvxvzbarz(0:nzmmnt*ihvxvzbarz,0:lenhist)  _real [(m/s)**2] limited (0:nzmmnt,0:jhist)
            +zhist           # VxVz bar versus space and time
lhvyvzbarz logical /.false./ # Turns on history of VyVz bar
ihvyvzbarz integer /0 /      # Multiplier for hvyvzbarz memory size (autoset)
hvyvzbarz(0:nzmmnt*ihvyvzbarz,0:lenhist)  _real [(m/s)**2] limited (0:nzmmnt,0:jhist)
            +zhist           # VyVz bar versus space and time

*********** Particles dump parallel:
# Dynamic particle arrays, and related data
np     integer    /0/  # Total no. of particles (including lost ones).
nplive integer    /0/  # No. of "live" particles
npmax  integer    /0/  # Maximum no. of particles
                       # (user input for some loadings)
npmaxb integer    /0/  # Maximum no. of particles for xp, yp, uxp, uyp
npidmax integer   /1/  # Maximum number of columns for pid.
                       # This is used so that the pid array is always allocated
npid   integer    /0/  # number of columns for pid.
npmaxi integer    /1/  # Maximum no. of particles for pid.
wpid   integer    /0/  # position of particle weights in array pid (FORTRAN indexed: based 1)
tpid   integer    /0/  # position of particle creation time in array pid (FORTRAN indexed: based 1)
rpid   integer    /0/  # position of particle initial radius in array pid (FORTRAN indexed: based 1)
sm(ns) _real [kg] /0./ -parallel # Species mass
sq(ns) _real [C]  /0./ -parallel # Species charge
sw(ns) _real [1]  /0./ -parallel # Species weight
                       # (real particles per simulation particles)
ins(ns)  _integer /1/  # Index of first particle in species
nps(ns)  _integer /0/  # Number of particles in species
gaminv(npmax)  _real  [1]   /1./ # inverse relativistic gamma factor
xp(npmaxb)     _real  [m]        # X-positions of particles
yp(npmaxb)     _real  [m]        # Y-positions of particles
zp(npmax)      _real  [m]        # Z-positions of particles
uxp(npmaxb)    _real  [m/s]      # gamma * X-velocities of particles
uyp(npmaxb)    _real  [m/s]      # gamma * Y-velocities of particles
uzp(npmax)     _real  [m/s]      # gamma * Z-velocities of particles
pid(npmaxi,npidmax) _real [1]    # Particle index - user for various purposes

*********** Scraped_Particles dump parallel:
# Arrays for scraped particles
scr_np             integer  /0/   # Total no. of scraped particles. 
scr_npmax          integer  /0/   # Max. no. of scraped particles.
scr_npbunch        integer  /10000/ # Number of particles in a bunch
scr_ns             integer  /1/   # number of species
scr_ins(scr_ns)        _integer /1/   # Index of first particle in species
scr_nps(scr_ns)        _integer /0/   # Number of particles in species
scr_xp(scr_npmax)  _real    [m]   # X-position
scr_yp(scr_npmax)  _real    [m]   # Y-position
scr_zp(scr_npmax)  _real    [m]   # Z-position
scr_uxp(scr_npmax) _real    [m/s] # gamma * X-velocities of particles
scr_uyp(scr_npmax) _real    [m/s] # gamma * Y-velocities of particles
scr_uzp(scr_npmax) _real    [m/s] # gamma * Z-velocities of particles

*********** LostParticles dump parallel:
lsavelostpart logical /.false./ # Flag setting whether lost particles are saved
npmaxlost           integer /0/ # Size of lost particle arrays
lostpartchunksize   integer /1000/
inslost(ns)        _integer /0/ # Index of first lost particles of species
npslost(ns)        _integer /0/ # Number of lost particles in species
npmaxlost_s(0:ns)  _integer /0/ # Max index of lost particles for species
xplost(npmaxlost)  _real [m]    # X-positions of lost particles
yplost(npmaxlost)  _real [m]    # Y-positions of lost particles
zplost(npmaxlost)  _real [m]    # Z-positions of lost particles
uxplost(npmaxlost) _real [m/s]  # gamma * X-velocities of lost particles
uyplost(npmaxlost) _real [m/s]  # gamma * Y-velocities of lost particles
uzplost(npmaxlost) _real [m/s]  # gamma * Z-velocities of lost particles
gaminvlost(npmaxlost) _real [1] # gamma inverse of lost particles
tplost(npmaxlost)  _real [s]    # time particles were lost
pidlost(npmaxlost,npidmax) _real [1] # Particle index of lost particles

*********** Picglb dump:
# Globally useful quantities for PIC simulation
time                      real  /0./
   # Problem time
zbeam                     real  [m]
   # Distance the center of the beam (actually, the "moving window") has moved
   # Advanced only when particles are advanced
zgrid                     real  [m]
   # Location of the center of the grid (actually, the "moving window") in the
   # lab frame.  Advanced only when IT is advanced
zgridprv                  real  [m]
   # Prvious location of the grid.  Needed for fetch of E-field from grid.
it                        integer           # Timestep counter
ldump                     logical
   # Flag set when this step will end in a restart dump
lfirst                    logical
   # Flag set when this is first step of a STEP or RUN command
llast                     logical
   # Flag set when this is last step of a STEP or RUN command
lfldplt                   logical
   # Flag set when fields will be plotted this step
lfinishd                  logical
   # Flag set when this step is the last of the run
lalways                   logical
   # Flag set when "always" plots are to be made
lseldom                   logical
   # Flag set when "seldom" plots are to be made
lmoments                  logical
   # Flag set when the moments are to be calculated
llabwn                    logical
   # Flag set when lab moments are to be calculated
lhist                     logical
   # Flag set when histories are to be save (so moments will be calculated)
lspecial                  logical
   # Flag set when this is a "special" timestep

*********** ExtPart dump:
# Arrays that hold particle data that is extrapolated to grid cell centers
# in the GETZMMNT routine.
nepwin                  integer /0/ # Number of grid locations (windows)
izepwin(nepwin)        _integer /-1/ # List of grid locations (indx of zmntmesh)
zzepwin(nepwin)        _real        # List of lab frame locations
wzepwin(nepwin)        _real        # List of lab frame widths
nepmax                  integer /0/ # Maximum number of particles
nep(nepwin,ns)         _integer     # Number of particles in each grid location
tep(nepmax,nepwin,ns)  _real        # time of particles at grid cell centers
xep(nepmax,nepwin,ns)  _real        # X coordinates at grid cell centers
yep(nepmax,nepwin,ns)  _real        # Y coordinates at grid cell centers
uxep(nepmax,nepwin,ns) _real        # X velocities at grid cell centers
uyep(nepmax,nepwin,ns) _real        # Y velocities at grid cell centers
uzep(nepmax,nepwin,ns) _real        # Z velocities at grid cell centers

*********** Pspwork:
# Scratch arrays for phase space plots
# (perhaps we should overlay some of this with other scratch space)
npsplt              integer /255/ # Number of particles to process at a time
                                  # 255 gives better "subsets" plots than 256
isubset(npsplt,NSUBSETS) _integer # Holds particle subsets for plotting
ntopick(NSUBSETS)         integer # Number of particles to pick from each clump
xxpsp(npsplt)               _real # Scratch for abscissa values of particles
yypsp(npsplt)               _real # Scratch for ordinate values of particles
xppsp(npsplt)               _real # Scratch for abscissa values of particles
yppsp(npsplt)               _real # Scratch for ordinate values of particles
xxpspw(npsplt,4)            _real # Scratch for abscissa values of particles
yypspw(npsplt,4)            _real # Scratch for ordinate values of particles

*********** TopPhys:
# "Physics" subroutines at top level
cigar (np:real,zunifrm:real,zpunifrm:real,z:real,zp:real,perpscal:real,
       straight:real,scrtch1:real,scrtch2:real,scrtch3:real,scrtch4:real)
            subroutine # Adjusts z and vz to make a finite, cigar beam
derivqty()  subroutine # Calculates global derived qtys.
getzmmnt(np,xp:real,yp:real,zp:real,uxp:real,uyp:real,uzp:real,gaminv:real,
         q:real,m:real,w:real,dt:real,itask,nplive,
         uxpo:real,uypo:real,uzpo:real,is:integer,ns:integer,
         maxp:real,minp:real,zmmnts0:real,zmmnts:real)
            subroutine # Sets moments as a function of z for species 1
getzmmnt_weights(np,xp:real,yp:real,zp:real,uxp:real,uyp:real,uzp:real,gaminv:real,
         wp:real,q:real,m:real,w:real,dt:real,itask,nplive,
         uxpo:real,uypo:real,uzpo:real,is:integer,ns:integer)
            subroutine # Sets moments as a function of z for species 1 with variables weights
periz(np,zp:real,zgrid:real,zmmax:real,zmmin:real)
            subroutine # Imposes periodicity on z
griddedparticlescraper(is:integer,distance:real,
                       nx:integer,ny:integer,nz:integer,
                       dx:real,dy:real,dz:real,xmin:real,ymin:real,zmin:real,
                       zbeam:real,l2symtry:logical,l4symtry:logical) subroutine
                       # General particle scraper which allows scraping
                       # in complex geometries.
setgamma(lrelativ)
            subroutine # Converts v to u, sets gammainv for all ptcls
gammaadv(np,gaminv:real,uxp:real,uyp:real,uzp:real,gamadv:real,lrelativ)
            subroutine # Advances gamma
resetlat()  subroutine # Resizes lattice arrays to their true lengths
setlatt()   subroutine # Sets lattice pointers for the current beam location
setlattzt(zbeam:real,time:real,fstype:integer)
            subroutine # Sets lattice pointers at zbeam and time
getelemid(z:real,offset:real,nelem:integer,elemzs:real,elemze:real,
          oelemnn:integer,oelemoi:integer,oelemio:integer,id:integer)
            subroutine # Gets id of element located nearest z
species()   subroutine # Sets species related arrays.
stckyz(np,zp:real,zmmax:real,zmmin:real,dz:real,uxp:real,uyp:real,uzp:real,
       gaminv:real,zgrid:real)
            subroutine # Enforces sticky z walls

*********** TopDiag:
# Subroutines in package TOP
setgrid1d(np:integer,x:real,nx:integer,grid:real,xmin:real,xmax:real)
        subroutine
        # Deposits data onto a 1-D grid.
deposgrid1d(itask:integer,np:integer,x:real,z:real,nx:integer,
            grid:real,gridcount,xmin:real,xmax:real)
        subroutine
        # Deposits data onto a 1-D grid.
getgrid1d(np:integer,x:real,z:real,nx:integer,grid:real,xmin:real,xmax:real)
        subroutine
        # Gathers data from a 1-D grid.
setgrid2d(np:integer,x:real,y:real,nx:integer,ny:integer,grid:real,
          xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Deposits uniform data onto a 2-D grid.
deposgrid2d(itask:integer,np:integer,x:real,y:real,z:real,nx:integer,ny:integer,
            grid:real,gridcount,xmin:real,xmax:real,ymin:real,ymax:real)
        subroutine
        # Deposits data onto a 2-D grid.
setgrid2dw(np:integer,x:real,y:real,w:real,nx:integer,ny:integer,grid:real,
          xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Deposits uniform data onto a 2-D grid.
deposgrid2dw(itask:integer,np:integer,x:real,y:real,z:real,w:real,nx:integer,ny:integer,
            grid:real,gridcount,xmin:real,xmax:real,ymin:real,ymax:real)
        subroutine
        # Deposits data onto a 2-D grid.
getgrid2d(np:integer,x:real,y:real,z:real,nx:integer,ny:integer,grid:real,
          xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Gathers data from a 2-D grid.
getgridngp2d(np:integer,x:real,y:real,z:real,nx:integer,ny:integer,grid:real,
          xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Gathers data from a 2-D grid using nearest grid point
setgrid3d(np:integer,x:real,y:real,z:real,nx:integer,ny:integer,nz:integer,
          grid:real,xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real)        subroutine
        # Deposits uniform data onto a 3-D grid.
deposgrid3d(itask:integer,np:integer,x:real,y:real,z:real,q:real,
            nx:integer,ny:integer,nz:integer,grid:real,gridcount,
            xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real)
        subroutine
        # Deposits data onto a 3-D grid.
getgrid3d(np:integer,x:real,y:real,z:real,f:real,
          nx:integer,ny:integer,nz:integer,grid:real,
          xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real,
          l2symtry:logical,l4symtry:logical) subroutine
        # Gathers data from a 3-D grid.
getgridngp3d(np:integer,x:real,y:real,z:real,f:real,
             nx:integer,ny:integer,nz:integer,grid:real,
             xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real,
             l2symtry:logical,l4symtry:logical) subroutine
        # Gathers data from a 3-D grid using nearest grid point.
grid2grid(unew:real,nxnew:integer,nynew:integer,
          xminnew:real,xmaxnew:real,yminnew:real,ymaxnew:real,
          uold:real,nxold:integer,nyold:integer,
          xminold:real,xmaxold:real,yminold:real,ymaxold:real) subroutine
        # project field from one grid to another
getpsgrd(np,xp:real,uxp:real,nw,nh,psgrd:real,wmin:real,wmax:real,hmin:real,
         hmax:real,zl:real,zr:real,zp:real,uzp:real,slope:real)
              subroutine # lays particles onto slanted mesh in phase space
emitthresh(n:integer,threshold:real,js:integer,iw:integer,
           ngridw:integer,ngridh:integer,tepsx:real,tepsy:real)
              subroutine # Calculates the emittance with thesholding.
        # --- Input:
        # ---   - n is number of thresholds
        # ---   - threshold are threshold values
        # ---   - js species of particles to include
        # ---   - iw z-window to select particles
        # ---   - ngridw number of grid points position is binned into
        # ---   - ngridh number of grid points velocity is binned into
        # --- Output:
        # ---   - tepsx,tepsy
emitellipse(xpshear:real,vxgam:real,npart:integer,xbar:real,xsqbar:real,
            xxpbar:real,xpsqbar:real,xpbar:real,npts:integer,vbeam:real,
            gamma:real,percent:real,emitt:real,
            xp:real,xprime:real,area:real,x2sum:real,xp2sum:real,xxpsum:real,
            xsum:real,xpsum:real,upercent:real,uemitt:real,nshell:integer)
              subroutine # Calculates the emittance versus percent of the
                         # beam current enclosed by nested ellipses in
                         # phase space.
unshear(xp:real,xprime:real,npart:integer,xsqbar:real,xpsqbar:real,
        xxpbar:real,xbar:real,xpbar:real,xun:real,xpun:real,vbeam:real,
        gamma:real)
              subroutine # Takes a phase ellipse which has non-zero <xx'> and
                         # rotates it to produce an ellipse with its axes
                         # along the coordinate axes and center at the same
                         # position as the original ellipse.
emitfrac(xp:real,uxp:real,uzp:real,np:integer,xbar:real,xpbar:real,
         xsqbar:real,xxpbar:real,xpsqbar:real,
         fracbin:real,emitbin:real,npts:integer,emitbinmax:real, 
         tx:real,txp:real,emitp:real,rwork:real,iwork:integer)
         subroutine # Calculates the emittance versus the fraction of the
                    # beam (live particles) particles in phase space 
                    # enclosed by nested emittance ellipses with the rms 
                    # equivalent beam.  Also returns single particle 
                    # emittances.
prin(xp:real,uxp:real,uzp:real,np:integer,xsqbar:real,xpsqbar:real,
     xxpbar:real,xbar:real,xpbar:real,tx:real,txp:real,
     txsqbar:real,txpsqbar:real,txxpbar:real)
         subroutine # Tranforms a phase ellipse which has non-zero <xx'>, 
                    # <x> and <x'> and translates and rotates it to 
                    # zero these moments in a mixed unit coordinate system.  
minidiag(it,tim:real,lmoments:logical)
              subroutine # Diagnostics that don't force a "special" step
stepid(it,time:real,zbeam:real)
              subroutine # Sets the pline3 string which is plotted at the bottom
                         # of every frame
thisstep(it,itcount,n) logical function
                         #
thiszbeam(zl:real,zr:real,zcount,n) logical function
                         #
savehist(time:real)
              subroutine # saves moments data to history arrays
prntpara(dxperp:real,dz:real)
              subroutine # Print various parameters to various files
psplots(freqflag:integer) 
              subroutine # Controls phase space plots
onedplts(freqflag:integer) subroutine # plots all 1d qtys w/ freqflag
tolabfrm(zcent:real,nn,x:real,z:real) subroutine
             # Converts data from WARP frame to lab frame.

*********** TopUtil:
# "Utility" subroutines at top level
copyarry(source:real,target:real,nwords)
             subroutine # Copies array from source to targer
dolabwn() logical function
                        # Checks if lab window is in beam frame
psumx(a:real,b:real,n)
             subroutine # b := partial sum of a
fnice(i,e10:real)
          real function   # makes nice contours values
alotpart()   subroutine # Allocate for particles and setup associated data
chckpart(is:integer,nlower:integer,nhigher:integer,lfullshft:logical)
             subroutine # Makes sure there is enough space for nn particles.
addpart(nn:integer,npid:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        gi:real,pid:real,
        is:integer,lallindomain:logical,zmmin:real,zmmax:real,
        lmomentum:logical)
             subroutine # Adds new particles to the simulation
clearpart(js:integer,fillmethod:integer)
             subroutine # Clears away lost particles.
processlostpart(is:integer,clearlostpart:integer,time:real,zbeam:real)
             subroutine # Processes lost particles (particles which have
                        # gaminv set to zero).
load2d(np,x:real,y:real,nx,ny,n:real,dx:real,dy:real)
             subroutine # Loads particles approximately into a 2-D distribution
shftpart(is:integer,ishft:integer) subroutine
                        # Moves particle data to end of species block.
copypart(it:integer,nn:integer,ii:integer;istart:integer) subroutine
                        # Copies particle data from one location to another
r2rev(lastrev:real)
          real function # bit reversed counter
                        # lastrev is `call by address', e.g. r2rev(&x)
rnrev(i:integer,nbase:integer)
          real function # i base `nbase' reversed.
rnrevarray(n:integer,x:real,i:integer,nbase:integer)
          subroutine    # Fills an array with uniform digit reversed rand numbs
rnorm()   real function # Gaussian random numbers.
rnormdig(i1,n,nbase1,nbase2,dx:real,x:real)
             subroutine # Gaussian random numbers via digit reversed.
rm()      real function # Pseudo-Gaussian random numbers (6 uniform nos.).
rma(a:real,n) subroutine # rma(&a,n) gives n Pseudo-Gaussian rand numbers.
sphere4(a:real,b:real,c:real,d:real,n:integer)
             subroutine # distrib pts on surf of 4d unit sphere
                        # a-d are arrays declared (1:n)
sphere4f(a:real,b:real,c:real,d:real,ig1,ig2,ig3)
             subroutine # distr pts on surf of 4d unit sphere
                        # a-d are arrays declared (1:n)
ssifa(a:real,lda,n,kpvt,info)
             subroutine # LINPACK matrix reduction routine # (in file UTIL.F)
ssidi(a:real,lda,n,kpvt,det:real,inert,work:real,job:real)
             subroutine # LINPACK matrix inverting routine # (in file UTIL.F)
ssisl(a:real,lda,n,kpvt,b:real)
             subroutine # LINPACK routine to solve a*x = b # (in file UTIL.F)
svdfit(x:real,y:real,ndata:integer,ndatap:integer,a:real,
       basis:real,ma:integer,map:integer,u:real,w:real,v:real,tmp:real,
       tol:real,chisq:real)    subroutine
       # matrix linear least squares fit using an input singular value 
       # decomposition of a basis matrix                   (in util.m) 
svdcmp(a:real,m:integer,n:integer,mp:integer,np:integer,
       w:real,v:real,tmp:real) subroutine
       # singular value decomposition of the matrix a      (in util.m)
svbksb(u:real,w:real,v:real,m:integer,n:integer,mp:integer,np:integer,
       b:real,x:real,tmp:real) subroutine
       # singular value back-substitution routine for solution of matrix 
       # problems                                          (in util.m)
vrpft2   (a:real, b:real, n, incp, lenv, lfd)  subroutine
         #  Vectorized Real Periodic Fourier Transform, 2 at a time
vrpfti2  (a:real, b:real, n, incp, lenv, lfd)  subroutine
         #  Vectorized Real Periodic Fourier Transform Inverse, 2 at a time
vcpft    (re:real, im:real, n, incp, signp:real, lenv, lfd) subroutine
         #  Vectorized Complex Periodic Fourier Transform
writarry(nn,arry:real,filename:string) subroutine # write array to file
wtime() real function # returns current absolute CPU time
wtimeon()  subroutine # turns timer on
wtimeoff() real function # returns time since last call to wtimeon or wtimeoff
wtremain() real function # returns the time remaining for the running job
                         # (T3E only, otherwise returns large number)
zbeamcom(zbeam:real) subroutine
       # Returns the center of mass in z of the beam calculated from the
       # particles.
zpartbnd(zmmax:real,zmmin:real,dz:real,zgrid:real) subroutine
       # Enforces axial particle boundary conditions
reorgparticles() subroutine
       # Reorganizes particles for the parallel version

******* Parallel:
nslaves       integer /0/         # Number of slaves
my_index      integer   +parallel # Processor index to array of task ids
grid_overlap  integer       +dump # Overlap of field grid in processors
slavenp       integer /0/   +dump # Value of npmax that slave is to use
maxslaves     integer /512/ +dump # Max numer of slaves
lautodecomp   logical /.true./    # When false, the domain decompostion for the
                                  # particles is supplied by the user.
lfsautodecomp logical /.true./    # When false, the domain decompostion for the
                                  # field solver is supplied by the user. Can
                                  # only be done for fstype == 3.
xynppgroup    integer /16/        # For slice field solver, number of process
                                  # in each group which cooperatively does
                                  # a field solve.
zslave(0:maxslaves-1)    _real +dump    # User supplied weighting for the domain
                                        # decomposition of the particles.
izslave(0:maxslaves-1)   _integer +dump # starting iz for each slave
izfsslave(0:maxslaves-1) _integer +dump # starting iz for which each slave does
                                        # a field solve calculation
nzslave(0:maxslaves-1)   _integer +dump # number of z grid cells for each slave
nzfsslave(0:maxslaves-1) _integer +dump # number of z grid cells for which each
                                        # slave does a field solve calculation
zmslmin(0:maxslaves-1)   _real    +dump # Mesh Z minimum for each slave
zmslmax(0:maxslaves-1)   _real    +dump # Mesh Z maximum for each slave
izpslave(0:maxslaves-1)  _integer +dump # Starting iz of particle extent
nzpslave(0:maxslaves-1)  _integer +dump # Number of Z cells of particle extent
zpslmin(0:maxslaves-1)   _real    +dump # Particle Z minimum for each slave
zpslmax(0:maxslaves-1)   _real    +dump # Particle Z maximum for each slave

******* Databuffers:
# Primarily used as data buffers for message passing
b1size integer /1/
buffer1(b1size) _real
b2size integer /1/
buffer2(b2size) _real
b3size integer /1/
buffer3(b3size) _real
b4size integer /1/
buffer4(b4size) _real
ib1size integer /1/
ibuffer1(ib1size) _integer
ib2size integer /1/
ibuffer2(ib2size) _integer
ib3size integer /1/
ibuffer3(ib3size) _integer
ib4size integer /1/
ibuffer4(ib4size) _integer
b2d1 integer /0/
b2d2 integer /0/
buffer2d(0:b2d1,0:b2d2,2) _real

******* SemiTransparentDisc dump:
# semitransparent disc data
n_STdiscs             integer /0/         # Number of semitransparent discs
z_STdiscs(n_STdiscs)  _real               # Positions of semitransparent discs
r_STdiscs(n_STdiscs)  _real               # Radii of semitransparent discs
t_STdiscs(n_STdiscs)  _real               # Transparency of semitransparent discs
semitransparent_disc(dz:real) subroutine # Randomly absorb particles passing through disc based on data in group SemiTransparentDisc.
					 # dz is the maximum distance traveled by particles in one time step.

*********** Temperatures dump:
nstmp integer /1/ # nb species
l_temp                 logical /.false./ # compute temperature in true
l_temp_collapseinz     logical /.false./ # collapse z-slices in temperature calculations
l_temp_rmcorrelations  logical /.true./  # remove x*v correlations in temperature calculations
l_accumulate_temperatures logical /.false./ # allows accumulation of data for temperatures
nxtslices    integer    /0/   # nb cells in x of temperature slices
nytslices    integer    /0/   # nb cells in y of temperature slices
nztslices    integer    /1/   # nb temperature slices
nxtslicesc   integer    /0/   # nb cells in x of temperature slices (correlation variables)
nytslicesc   integer    /0/   # nb cells in y of temperature slices (correlation variables)
nztslicesc   integer    /1/   # nb temperature slices               (correlation variables)
tslicexmin(nztslices) _real   # min in x of temperature slices
tslicexmax(nztslices) _real   # max in x of temperature slices
tsliceymin(nztslices) _real   # min in y of temperature slices
tsliceymax(nztslices) _real   # max in y of temperature slices
tslicezmin(nztslices) _real   # min in z of temperature slices
tslicezmax(nztslices) _real   # max in z of temperature slices
dxti(nztslices)       _real   # inverse mesh size in x for temperature slices
dyti(nztslices)       _real   # inverse mesh size in y for temperature slices
pnumt(0:nxtslices,0:nytslices,nztslices)      _real # nb particles    in temperature slices
pnumtw(0:nxtslices,0:nytslices,nztslices)     _real # weights         in temperature slices
vxbart(0:nxtslices,0:nytslices,nztslices)     _real # average of vx   in temperature slices
vybart(0:nxtslices,0:nytslices,nztslices)     _real # average of vy   in temperature slices
vzbart(0:nxtslices,0:nytslices,nztslices)     _real # average of vz   in temperature slices
vxsqbart(0:nxtslices,0:nytslices,nztslices)   _real # average of vx^2 in temperature slices
vysqbart(0:nxtslices,0:nytslices,nztslices)   _real # average of vy^2 in temperature slices
vzsqbart(0:nxtslices,0:nytslices,nztslices)   _real # average of vz^2 in temperature slices
xbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of x    in temperature slices
ybart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of y    in temperature slices
zbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of z    in temperature slices
xsqbart(0:nxtslicesc,0:nytslicesc,nztslicesc) _real # average of x^2  in temperature slices
ysqbart(0:nxtslicesc,0:nytslicesc,nztslicesc) _real # average of y^2  in temperature slices
zsqbart(0:nxtslicesc,0:nytslicesc,nztslicesc) _real # average of z^2  in temperature slices
xvxbart(0:nxtslicesc,0:nytslicesc,nztslicesc) _real # average of x*vx in temperature slices
yvybart(0:nxtslicesc,0:nytslicesc,nztslicesc) _real # average of y*vy in temperature slices
zvzbart(0:nxtslicesc,0:nytslicesc,nztslicesc) _real # average of z*vz in temperature slices
tempx(0:nxtslices,0:nytslices,nztslices,nstmp) _real # x-temperature  in temperature slices
tempy(0:nxtslices,0:nytslices,nztslices,nstmp) _real # y-temperature  in temperature slices
tempz(0:nxtslices,0:nytslices,nztslices,nstmp) _real # z-temperature  in temperature slices
tempxz(nztslices,nstmp) _real # x-temperature  in temperature slices
tempyz(nztslices,nstmp) _real # y-temperature  in temperature slices
tempzz(nztslices,nstmp) _real # z-temperature  in temperature slices
nztlocator integer /1/ # nb meshes locator temperature slices
ntlmax     integer /1/ # max nb of slices in locator temperature slices
ntl(nztlocator) _integer # nb of slices in locator temperature slices
tslice_locator(nztlocator,ntlmax) _integer # locator array for temperature array
tloc_zmin  real # min locator temperature array 
tloc_zmax  real # max locator temperature array 
tloc_dzi   real # inverse mesh size locator temperature array 
gett(is:integer,itask:integer) subroutine # get temperature for species is 
setregulartgrid(nx:integer,ny:integer,nz:integer,
                xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real,
                dz:real,nzloc:integer,lcollapse:logical,lcorrel:logical) subroutine 
                # set slices regularly in z between zmin and zmax.
                # nzloc refers to the number of nodes of a lookup table for fast 
                # localization of temperature slices. In general setting nzloc=w3d.nz is good.
