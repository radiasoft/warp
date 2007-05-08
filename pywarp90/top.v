top
#@(#) File TOP.V, version $Revision: 3.202 $, $Date: 2007/05/08 16:44:13 $
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
verstop character*19 /"$Revision: 3.202 $"/ # Global common version, set by CVS

*********** Machine_param:
wordsize integer /64/ # Wordsize on current machine--used in bas.wrp
largepos real    /LARGEPOS/ # Large positive number
smallpos real    /SMALLPOS/ # Small positive number

*********** GlobalVars:
nparpgrp  integer /NPARPGRP/ +dump # Number of particles per group. Effects
                                   # size of temporaries and cache use.
dirichlet integer /0/        # constant for boundary condition (constant potential)
neumann   integer /1/        # constant for boundary condition (derivative = 0)
periodic  integer /2/        # constant for boundary condition 
absorb    integer /0/        # constant for particle absorption at boundaries
reflect   integer /1/        # constant for particle reflection at boundaries

*********** DebugFlags dump:
debug logical /.false./ # When true, extensive debugging is done.

*********** Timers dump parallel:
starttime       real /0./ -dump # CPU start time (in seconds)
starttimedump   real /0./ # CPU start time (in seconds)
gentime         real /0./ # CPU for generate (in seconds)
steptime        real /0./ # Total CPU run time minus gentime (in seconds)
plottime        real /0./ # Time making automatic plots (in seconds)
momentstime     real /0./ # Time to calculate the moments (in seconds)
fstime          real /0./ # Time to do the field solve, including beforefs and
                          # afterfs (in seconds)
latticetime     real /0./ # Time to apply the fields from the lattice
dumptime        real /0./ # Time to do the data dumps, using the dump command
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
negrd     integer /NELEMENT/  # No. of elements with E field on a 3-D grid
nbgrd     integer /NELEMENT/  # No. of elements with B field on a 3-D grid
npgrd     integer /NELEMENT/  # No. of elements with potential on a 3-D grid
nbsqgrad  integer /NELEMENT/  # No. of elements with grad B^2 field on a 3-D grid
drftzs(0:ndrft)   _real [m]   # Z's of drift starts
drftze(0:ndrft)   _real [m]   # Z's of drift ends
drftap(0:ndrft)   _real [m]   # Aperture in drifts
drftax(0:ndrft)   _real [m]   # Aperture in drifts in x
drftay(0:ndrft)   _real [m]   # Aperture in drifts in y
drftox(0:ndrft)   _real [m]   # X-offsets of drifts
drftoy(0:ndrft)   _real [m]   # Y-offsets of drifts
drftol(0:ndrft)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
bendzs(0:nbend)   _real [m]   # Z's of bend starts
bendze(0:nbend)   _real [m]   # Z's of bend ends
bendrc(0:nbend)   _real [m]   # Radii of curvature of bends
bendap(0:nbend)   _real [m]   # Aperture in bends
bendax(0:nbend)   _real [m]   # Aperture in bends in x
benday(0:nbend)   _real [m]   # Aperture in bends in y
bendox(0:nbend)   _real [m]   # X-offsets of bends (for aperture only)
bendoy(0:nbend)   _real [m]   # Y-offsets of bends (for aperture only)
bendol(0:nbend)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
dipozs(0:ndipo)   _real [m]   # Z's of dipo starts (set from bendzs if =dipoze)
dipoze(0:ndipo)   _real [m]   # Z's of dipo ends   (set from bendze if =dipozs)
dipoby(0:ndipo)   _real [T]   # By's of dipos (set from bendrc if 0 & dipoex=0)
dipobx(0:ndipo)   _real [T]   # Bx's of dipos
dipota(0:ndipo)   _real [1]   # Tan of dipo entry face angle (auto-set if 0)
dipotb(0:ndipo)   _real [1]   # Tan of dipo exit  face angle (auto-set if 0)
dipoex(0:ndipo)   _real [V/m] # Ex's of dipos
dipoey(0:ndipo)   _real [V/m] # Ey's of dipos
dipoap(0:ndipo)   _real [m]   # Aperture in dipoles
dipoax(0:ndipo)   _real [m]   # Aperture in dipoles in x
dipoay(0:ndipo)   _real [m]   # Aperture in dipoles in y
dipox1(0:ndipo)   _real [m]   # X location of first dipole plates
dipox2(0:ndipo)   _real [m]   # X location of second dipole plates
dipov1(0:ndipo)   _real [V]   # Voltage of first dipole plates
dipov2(0:ndipo)   _real [V]   # Voltage of second dipole plates
dipol1(0:ndipo)   _real [m]   # Length of first dipole plates
dipol2(0:ndipo)   _real [m]   # Length of second dipole plates
dipow1(0:ndipo)   _real [m]   # Width of first dipole plates
dipow2(0:ndipo)   _real [m]   # Width of second dipole plates
dipool(0:ndipo)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
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
quadax(0:nquad)   _real [m]   # Aperture of quad in x
quaday(0:nquad)   _real [m]   # Aperture of quad in y
quadrr(0:nquad)   _real [m]   # Radius of electrostatic quadrupole rod
quadrl(0:nquad)   _real [m]   # Length of electrostatic quadrupole rod
quadgl(0:nquad)   _real [m]   # Length of electrostatic quadrupole gap
quadgp(0:nquad)   _real [ ]   # Gap position of ESQ, only sign is used
quadpw(0:nquad)   _real [m]   # End plate width of electrostatic quadrupole
quadpa(0:nquad)   _real [m]   # End plate aperture of electrostatic quadrupole
quadpr(0:nquad)   _real [m]   # End plate max radius
quadsl(0:nquad)   _real [m]   # Slant of rod which makes it a cone
quadol(0:nquad)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
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
sextol(0:nsext)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
helezs(0:nhele)   _real   [m] # Z's of hard-edge (h.e.) multipole starts
heleze(0:nhele)   _real   [m] # Z's of hard-edge (h.e.) multipole ends
heleap(0:nhele)   _real [m]   # Aperture in hard-edge elements
heleax(0:nhele)   _real [m]   # Aperture in hard-edge elements in x
heleay(0:nhele)   _real [m]   # Aperture in hard-edge elements in y
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
heleol(0:nhele)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
emltzs(0:nemlt)   _real [m]   # Z's of electric multipole element starts
emltze(0:nemlt)   _real [m]   # Z's of electric multipole element ends
emltap(0:nemlt)   _real [m]   # Aperture in electric multipole elements
emltax(0:nemlt)   _real [m]   # Aperture in electric multipole elements in x
emltay(0:nemlt)   _real [m]   # Aperture in electric multipole elements in y
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
emltol(0:nemlt)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
mmltzs(0:nmmlt)   _real [m]   # Z's of magnetic multipole element starts
mmltze(0:nmmlt)   _real [m]   # Z's of magnetic multipole element ends
mmltap(0:nmmlt)   _real [m]   # Aperture in magnetic multipole elements
mmltax(0:nmmlt)   _real [m]   # Aperture in magnetic multipole elements in x
mmltay(0:nmmlt)   _real [m]   # Aperture in magnetic multipole elements in y
mmltas(0:nmmlt)   _real [m]   # Z start's of aperture in mag multipole elements
mmltae(0:nmmlt)   _real [m]   # Z end's of aperture in mag multipole elements
mmltph(0:nmmlt)   _real [rad] # Phase angle of magnetic multipole element field
mmltsf(0:nmmlt)   _real [1] /0./ # Scale factor for magnetic multipole element
                                 # Field is scaled by (mmltsc+mmltsf)
mmltsc(0:nmmlt)   _real [1] /1./ # Scale factor for magnetic multipole element
                                 # Field is scaled by (mmltsc+mmltsf)
mmltid(0:nmmlt)   _integer    # Index of magnetic multipole dataset
mmltox(0:nmerr)   _real [m]   # Offset in x of magnetic multipole centers
mmltoy(0:nmerr)   _real [m]   # Offset in y of magnetic multipole centers
mmltol(0:nmmlt)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
acclzs(0:naccl)   _real [m]   # Z's of acceleration gap starts
acclze(0:naccl)   _real [m]   # Z's of acceleration gap ends
acclez(0:naccl)   _real [V/m] # Ez's of acceleration gaps, constant part
acclap(0:naccl)   _real [m]   # Aperture in acceleration gaps
acclax(0:naccl)   _real [m]   # Aperture in acceleration gaps in x
acclay(0:naccl)   _real [m]   # Aperture in acceleration gaps in y
acclox(0:naccl)   _real [m]   # Offset in x of accl elements
accloy(0:naccl)   _real [m]   # Offset in y of accl elements
acclxw(0:naccl) _real [V/m^2] # Weights for linear x components of gap Ez's
# Acceleration gap Ez's are of the form, E_z = acclez + acclxw*x.
acclsw(0:naccl)   _real [1]   # Switch for grid accel by gaps (0 = on,1 = off)
acclet(0:ntaccl,0:naccl) _real [V/m] # Ez's of acceleration gaps as a function
                                     # of time.
acclts(0:naccl)   _real [t]   # Time of start of gap field in acclet
accldt(0:naccl)   _real [t]   # Delta t for gap field data in acclet
acclol(0:naccl)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
egrdzs(0:negrd)   _real [m]   # Z starts of 3-D grid of E field data (EGRDdata)
egrdze(0:negrd)   _real [m]   # Z ends of 3-D grid of E field data (EGRDdata)
egrdxs(0:negrd)   _real [m]   # X starts of 3-D grid of E field data (EGRDdata)
egrdys(0:negrd)   _real [m]   # Y starts of 3-D grid of E field data (EGRDdata)
egrdap(0:negrd)   _real [m]   # Aperture in egrd elements
egrdax(0:negrd)   _real [m]   # Aperture in egrd elements in x
egrday(0:negrd)   _real [m]   # Aperture in egrd elements in y
egrdox(0:negrd)   _real [m]   # Offset in x of egrd elements
egrdoy(0:negrd)   _real [m]   # Offset in y of egrd elements
egrdph(0:negrd)   _real [rad] # Phase angle of egrd elements 
egrdsp(0:negrd)   _real [1]   # Sine   of egrdph (auto-set) 
egrdcp(0:negrd)   _real [1]   # Cosine of egrdph (auto-set) 
egrdid(0:negrd)   _integer    # Index of to 3-D E field data sets (EGRDdata)
egrdsf(0:negrd) _real [1] /0./ # Scale factor to multiply 3-D E field data set
                               # EGRDdata. Field is scaled by (egrdsc+egrdsf)
egrdsc(0:negrd) _real [1] /1./ # Scale factor to multiply 3-D E field data set
                               # EGRDdata. Field is scaled by (egrdsc+egrdsf)
egrdsy(0:negrd) _integer /0/   # Level of symmetry in the egrd data.
                               # (0, no symmetry; 2, quadrupole)
                               # Defaul is no symmetry.
egrdol(0:negrd)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
bgrdzs(0:nbgrd)   _real [m]   # Z starts of 3-D grid of B field data (BGRDdata)
bgrdze(0:nbgrd)   _real [m]   # Z ends of 3-D grid of B field data (BGRDdata)
bgrdxs(0:nbgrd)   _real [m]   # X starts of 3-D grid of B field data (BGRDdata)
bgrdys(0:nbgrd)   _real [m]   # Y starts of 3-D grid of B field data (BGRDdata)
bgrdap(0:nbgrd)   _real [m]   # Aperture in bgrd elements
bgrdax(0:nbgrd)   _real [m]   # Aperture in bgrd elements in x
bgrday(0:nbgrd)   _real [m]   # Aperture in bgrd elements in y
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
bgrdol(0:nbgrd)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
pgrdzs(0:npgrd)   _real [m]   # Z starts of 3-D grid of potential data(PGRDdata)
pgrdze(0:npgrd)   _real [m]   # Z ends of 3-D grid of potential data(PGRDdata)
pgrdxs(0:npgrd)   _real [m]   # X starts of 3-D grid of potential data(PGRDdata)
pgrdys(0:npgrd)   _real [m]   # Y starts of 3-D grid of potential data(PGRDdata)
pgrdap(0:npgrd)   _real [m]   # Aperture in pgrd elements
pgrdax(0:npgrd)   _real [m]   # Aperture in pgrd elements in x
pgrday(0:npgrd)   _real [m]   # Aperture in pgrd elements in y
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
pgrdol(0:npgrd)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
bsqgradzs(0:nbsqgrad)   _real [m]   # Z starts of 3-D grid of grad B^2 field data (BSQGRADdata)
bsqgradze(0:nbsqgrad)   _real [m]   # Z ends of 3-D grid of grad B^2 field data (BSQGRADdata)
bsqgradxs(0:nbsqgrad)   _real [m]   # X starts of 3-D grid of grad B^2 field data (BSQGRADdata)
bsqgradys(0:nbsqgrad)   _real [m]   # Y starts of 3-D grid of grad B^2 field data (BSQGRADdata)
bsqgradap(0:nbsqgrad)   _real [m]   # Aperture in bsqgrad elements
bsqgradax(0:nbsqgrad)   _real [m]   # Aperture in bsqgrad elements in x
bsqgraday(0:nbsqgrad)   _real [m]   # Aperture in bsqgrad elements in y
bsqgradox(0:nbsqgrad)   _real [m]   # Offset in x of bsqgrad elements
bsqgradoy(0:nbsqgrad)   _real [m]   # Offset in y of bsqgrad elements
bsqgradph(0:nbsqgrad)   _real [rad] # Phase angle of bsqgrad elements 
bsqgradsp(0:nbsqgrad)   _real [1]   # Sine   of bsqgradph (auto-set) 
bsqgradcp(0:nbsqgrad)   _real [1]   # Cosine of bsqgradph (auto-set) 
bsqgradid(0:nbsqgrad)   _integer    # Index of to 3-D grad B^2 field data sets (BSQGRADdata)
bsqgradsf(0:nbsqgrad) _real [1] /0./ # Scale factor to multiply 3-D grad B^2 field data set
                               # BSQGRADdata. Field is scaled by (bsqgradsc+bsqgradsf)
bsqgradsc(0:nbsqgrad) _real [1] /1./ # Scale factor to multiply 3-D grad B^2 field data set
                               # BSQGRADdata. Field is scaled by (bsqgradsc+bsqgradsf)
bsqgradsy(0:nbsqgrad) _integer /0/   # Level of symmetry in the bsqgrad data.
                               # (0, no symmetry; 2, quadrupole)
                               # Defaul is no symmetry.
bsqgradol(0:nbsqgrad)   _integer    # Overlap level of the element (autoset).
                              # Set to -1 to ignore overlaps.
drfts     logical             # Flag for existence of drfts (auto set)
bends     logical             # Flag for existence of bends (auto set)
dipos     logical             # Flag for existence of dipos (auto set)
quads     logical             # Flag for existence of quads (auto set)
sexts     logical             # Flag for existence of sexts (auto set)
heles     logical             # Flag for exist. of hard-edge mults (auto set)
accls     logical             # Flag for existence of accel (auto set)
emlts     logical             # Flag for existence of emlts (auto set)
mmlts     logical             # Flag for existence of mmlts (auto set)
egrds     logical             # Flag for existence of egrds (auto set)
bgrds     logical             # Flag for existence of bgrds (auto set)
pgrds     logical             # Flag for existence of pgrds (auto set)
bsqgrads  logical             # Flag for existence of bsqgrads (auto set)
diposet   logical  /.true./   # Auto-set dipoles from bend locations and radii 
ldipocond logical  /.false./  # Auto generate flat dipole plates from the
                              # dipole parameters

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

******** EGRDdata dump:
# Data for the 3-D E field lattice element
egrdnx integer /0/ # Number of X cells
egrdny integer /0/ # Number of Y cells
egrdnz integer /0/ # Number of Z cells
egrdns integer /0/ # Number of data sets
egrdnc integer /0/ # Number of components stored in egrd
egrddx(egrdns)                      _real [m]   # X cell size
egrddy(egrdns)                      _real [m]   # Y cell size
egrddz(egrdns)                      _real [m]   # Z cell size
egrddxi(egrdns)                     _real [1/m] # 1 over X cell size (autoset)
egrddyi(egrdns)                     _real [1/m] # 1 over Y cell size (autoset)
egrddzi(egrdns)                     _real [1/m] # 1 over Z cell size (autoset)
egrdex(0:egrdnx,0:egrdny,0:egrdnz,egrdns) _real [T] # Ex
egrdey(0:egrdnx,0:egrdny,0:egrdnz,egrdns) _real [T] # Ey
egrdez(0:egrdnx,0:egrdny,0:egrdnz,egrdns) _real [T] # Ez

******** BGRDdata dump:
# Data for the 3-D B field lattice element
bgrdnx integer /0/ # Number of X cells
bgrdny integer /0/ # Number of Y cells
bgrdnz integer /0/ # Number of Z cells
bgrdns integer /0/ # Number of data sets
bgrdnc integer /0/ # Number of components stored in bgrd
bgrddx(bgrdns)                      _real [m]   # X cell size
bgrddy(bgrdns)                      _real [m]   # Y cell size
bgrddz(bgrdns)                      _real [m]   # Z cell size
bgrddxi(bgrdns)                     _real [1/m] # 1 over X cell size (autoset)
bgrddyi(bgrdns)                     _real [1/m] # 1 over Y cell size (autoset)
bgrddzi(bgrdns)                     _real [1/m] # 1 over Z cell size (autoset)
bgrdbx(0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns) _real [T] # Bx
bgrdby(0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns) _real [T] # By
bgrdbz(0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns) _real [T] # Bz

******** BSQGRADdata dump:
# Data for the 3-D grad B^2 field 'lattice element'
bsqgradnx integer /0/ # Number of X cells
bsqgradny integer /0/ # Number of Y cells
bsqgradnz integer /0/ # Number of Z cells
bsqgradns integer /0/ # Number of data sets
bsqgradnc integer /0/ # Number of components stored in bsqgrad
bsqgradntemp integer /0/ # Number of temporary components stored for
                         # automatic calculation of grad B dot B
bsqgraddx(bsqgradns)                      _real [m]   # X cell size
bsqgraddy(bsqgradns)                      _real [m]   # Y cell size
bsqgraddz(bsqgradns)                      _real [m]   # Z cell size
bsqgraddxi(bsqgradns)                     _real [1/m] # 1 over X cell size (autoset)
bsqgraddyi(bsqgradns)                     _real [1/m] # 1 over Y cell size (autoset)
bsqgraddzi(bsqgradns)                     _real [1/m] # 1 over Z cell size (autoset)
bsqgrad(bsqgradnc,0:bsqgradnx,0:bsqgradny,0:bsqgradnz,bsqgradns) _real
  # B field stored in interleaved array. First dimension is changable,
  # allowing option of storing additional data.
bsqgradtemp(0:bsqgradnx,0:bsqgradny,0:bsqgradnz,0:bsqgradntemp-1,bsqgradns) _real
calculatebsqgrad() subroutine # Calculates bsqgrad from the current B fields

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
ndrftol            integer /0/ # Maximum level of overlapping drft elements
odrftoi(0:ndrft)   _integer     # Overlap indices for drft elements
odrftio(0:ndrft)   _integer     # Overlap indices for drft elements
odrftii(ndrftol)   _integer     # Overlap indices for drft elements
odrftnn(ndrftol)   _integer     # Number of drft elements in overlap levels
nbendol            integer /0/ # Maximum level of overlapping bend elements
obendoi(0:nbend)   _integer     # Overlap indices for bend elements
obendio(0:nbend)   _integer     # Overlap indices for bend elements
obendii(nbendol)   _integer     # Overlap indices for bend elements
obendnn(nbendol)   _integer     # Number of bend elements in overlap levels
ndipool            integer /0/ # Maximum level of overlapping dipo elements
odipooi(0:ndipo)   _integer     # Overlap indices for dipo elements
odipoio(0:ndipo)   _integer     # Overlap indices for dipo elements
odipoii(ndipool)   _integer     # Overlap indices for dipo elements
odiponn(ndipool)   _integer     # Number of dipo elements in overlap levels
nquadol            integer /0/ # Maximum level of overlapping quad elements
oquadoi(0:nquad)   _integer     # Overlap indices for quad elements
oquadio(0:nquad)   _integer     # Overlap indices for quad elements
oquadii(nquadol)   _integer     # Overlap indices for quad elements
oquadnn(nquadol)   _integer     # Number of quad elements in overlap levels
nsextol            integer /0/ # Maximum level of overlapping sext elements
osextoi(0:nsext)   _integer     # Overlap indices for sext elements
osextio(0:nsext)   _integer     # Overlap indices for sext elements
osextii(nsextol)   _integer     # Overlap indices for sext elements
osextnn(nsextol)   _integer     # Number of sext elements in overlap levels
nheleol            integer /0/ # Maximum level of overlapping hele elements
oheleoi(0:nhele)   _integer     # Overlap indices for hele elements
oheleio(0:nhele)   _integer     # Overlap indices for hele elements
oheleii(nheleol)   _integer     # Overlap indices for hele elements
ohelenn(nheleol)   _integer     # Number of hele elements in overlap levels
nemltol            integer /0/ # Maximum level of overlapping emlt elements
oemltoi(0:nemlt)   _integer     # Overlap indices for emlt elements
oemltio(0:nemlt)   _integer     # Overlap indices for emlt elements
oemltii(nemltol)   _integer     # Overlap indices for emlt elements
oemltnn(nemltol)   _integer     # Number of emlt elements in overlap levels
nmmltol            integer /0/ # Maximum level of overlapping mmlt elements
ommltoi(0:nmmlt)   _integer     # Overlap indices for mmlt elements
ommltio(0:nmmlt)   _integer     # Overlap indices for mmlt elements
ommltii(nmmltol)   _integer     # Overlap indices for mmlt elements
ommltnn(nmmltol)   _integer     # Number of mmlt elements in overlap levels
nacclol            integer /0/ # Maximum level of overlapping accl elements
oaccloi(0:naccl)   _integer     # Overlap indices for accl elements
oacclio(0:naccl)   _integer     # Overlap indices for accl elements
oacclii(nacclol)   _integer     # Overlap indices for accl elements
oacclnn(nacclol)   _integer     # Number of accl elements in overlap levels
negrdol            integer /0/ # Maximum level of overlapping egrd elements
oegrdoi(0:negrd)   _integer     # Overlap indices for egrd elements
oegrdio(0:negrd)   _integer     # Overlap indices for egrd elements
oegrdii(negrdol)   _integer     # Overlap indices for egrd elements
oegrdnn(negrdol)   _integer     # Number of egrd elements in overlap levels
nbgrdol            integer /0/ # Maximum level of overlapping bgrd elements
obgrdoi(0:nbgrd)   _integer     # Overlap indices for bgrd elements
obgrdio(0:nbgrd)   _integer     # Overlap indices for bgrd elements
obgrdii(nbgrdol)   _integer     # Overlap indices for bgrd elements
obgrdnn(nbgrdol)   _integer     # Number of bgrd elements in overlap levels
npgrdol            integer /0/ # Maximum level of overlapping pgrd elements
opgrdoi(0:npgrd)   _integer     # Overlap indices for pgrd elements
opgrdio(0:npgrd)   _integer     # Overlap indices for pgrd elements
opgrdii(npgrdol)   _integer     # Overlap indices for pgrd elements
opgrdnn(npgrdol)   _integer     # Number of pgrd elements in overlap levels
nbsqgradol            integer /0/ # Maximum level of overlapping bsqgrad elements
obsqgradoi(0:nbsqgrad)   _integer     # Overlap indices for bsqgrad elements
obsqgradio(0:nbsqgrad)   _integer     # Overlap indices for bsqgrad elements
obsqgradii(nbsqgradol)   _integer     # Overlap indices for bsqgrad elements
obsqgradnn(nbsqgradol)   _integer     # Number of bsqgrad elements in overlap levels
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
cegrdzs(0:nzlmax,negrdol)    _real [m]     # by z, Z's of 3-D E field start
cegrdze(0:nzlmax,negrdol)    _real [m]     # by z, Z's of 3-D E field end
cegrdid(0:nzlmax,negrdol) _integer [1]     # by z, Index of egrd arrays
cbgrdzs(0:nzlmax,nbgrdol)    _real [m]     # by z, Z's of 3-D B field start
cbgrdze(0:nzlmax,nbgrdol)    _real [m]     # by z, Z's of 3-D B field end
cbgrdid(0:nzlmax,nbgrdol) _integer [1]     # by z, Index of bgrd arrays
cpgrdzs(0:nzlmax,npgrdol)    _real [m]     # by z, Z's of 3-D potential start
cpgrdze(0:nzlmax,npgrdol)    _real [m]     # by z, Z's of 3-D potential end
cpgrdid(0:nzlmax,npgrdol) _integer [1]     # by z, Index of pgrd arrays
cbsqgradzs(0:nzlmax,nbsqgradol)    _real [m]     # by z, Z's of 3-D B field start
cbsqgradze(0:nzlmax,nbsqgradol)    _real [m]     # by z, Z's of 3-D B field end
cbsqgradid(0:nzlmax,nbsqgradol) _integer [1]     # by z, Index of bsqgrad arrays
lindrft(0:ndrftol)        _logical         # Flag for when drft element in mesh
linbend                    logical         # Flag for when bend element in mesh
lindipo(0:ndipool)        _logical         # Flag for when dipo element in mesh
linquad(0:nquadol)        _logical         # Flag for when quad element in mesh
linsext(0:nsextol)        _logical         # Flag for when sext element in mesh
linhele(0:nheleol)        _logical         # Flag for when hele element in mesh
linemlt(0:nemltol)        _logical         # Flag for when emlt element in mesh
linmmlt(0:nmmltol)        _logical         # Flag for when mmlt element in mesh
linaccl(0:nacclol)        _logical         # Flag for when accl element in mesh
linegrd(0:negrdol)        _logical         # Flag for when egrd element in mesh
linbgrd(0:nbgrdol)        _logical         # Flag for when bgrd element in mesh
linpgrd(0:npgrdol)        _logical         # Flag for when pgrd element in mesh
linbsqgrad(0:nbsqgradol)        _logical         # Flag for when bsqgrad element in mesh

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
lspeciesmoments           logical /.true./
   # When true, the various moments are calculated separately for each
   # species, as well as for all of the species combined. The species list
   # will be the last index in the moments arrays, with the last element
   # being all species combined.
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
nskipcol integer /5/       # number of particles to skip in color plots
nstepcol integer /10/      # number of sweeps to plots colot plots
ncolor   integer /5/       # number of colors in color plots
iscolor  integer /2/       # starting color number in color plots
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
bx0                       real [T]   /0./
   # Uniform B field in x
by0                       real [T]   /0./
   # Uniform B field in y
bz0                       real [T]   /0./
   # Uniform axial "guide" B field
ex0                       real /0./
   # Uniform E field in x
ey0                       real /0./
   # Uniform E field in y
ez0                       real /0./
   # Linear Ez field for bunching force
vbeamfrm                  real [M/S]
   # Velocity of the beam frame.  Normally the same as vbeam.
allspecl                  logical  /.false./
   # flag for making all time steps "special" ones
depos                     character*8 /"direct1"/
   # Specifies charge deposition algorithm, "scalar", "vector", "vector1",
   # "direct", "direct1", "dspline2".
laccumulate_rho           logical /.false./
   # When true, rho is accumulated over multiple time steps. (i.e. it is not
   # zeroed each time step before the deposition.)
gamadv                    character*8 /"stndrd"/
   # Specifies type of gamma advance, "stndrd", "fast 1", or "fast 2"
lvdadz                    logical /.false./
   # sets on/off calculation of e=dA/dt=vdA/dz
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
   # Specifies type of field solve
   # -1: none
   #  0:sine-sine-periodic FFT,
   #  1: 8-fold symmetric capacity matrix in kz space,
   #  2: capacity matrix for quadrupoles,
   #  3: SOR (obsolete, use 7, multigrid),
   #  4: 2d sine-sine FFT + tridiag in z,
   #  5: general capacity matrix in kz space,
   #  6: general capacity matrix,
   #  7: multigrid solver,
   #  8: parallel solver (in development, don't use),
   #  9: parallel solver (in development, don't use),
   # 10: RZ multigrid solver
   # 11: Chombo AMR/multigrid solver
   # 12: Use field solver registered in python
   # 13: 3d multigrid with Boltzmann electrons
bfstype                   integer /-1/
   # Specifies type of field solver to use to calculate the B fields
   # -1: none
   #  0:sine-sine-periodic FFT,
   #  4: 2d sine-sine FFT + tridiag in z,
   #  7: multigrid solver,
   # 12: Use field solver registered in python
ztransformtype            integer /0/
   # Specifies the type of z transform
   #  0: periodic FFT
   #  1: tridiagonal matrix solve
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
   # When true, zbeam follows the beam center of mass (minus zbeamcomoffset).
zbeamcomoffset            real /0./
   # Offset of zbeam relative to the beam center of mass when lbeamcom is true.
relativity                integer /0/
   # Level of relativitistic corrections.
   #  1: scale transverse self E-field by 1/gamma**2
lcallfetchb               logical /.false./
   # When true, the fetchb routine is called during the time step
clearlostpart             integer /1/
   # When 0, do not clear lost particles, when 1, swap lost particles with
   # ones at the end of the array (much faster), when 2, shift particles to
   # fill in the gaps (much slower)
courantmin                     real    /0.25/
   # minimum fraction of particle courant condition under which particle time step is multiplied by two
courantmax                     real    /0.5/
   # maximum fraction of particle courant condition above which particle time step is divided by two
defieldmax                     real    /0.5/
   # maximum fraction of electric field change in last time step allowed 
dbfieldmax                     real    /0.5/
   # maximum fraction of magnetic field change in last time step allowed 

*********** InPart dump:
# Particle input quantities (input qtys)
ns              integer         /1/
   # Number of species
np_s(ns)       _integer [1]     /0/
   # Number of particles by species: if zero, set to npmax; if npmax is zero,
   # it is set to sum(np_s).
idkinterp(ns)  _integer        /0/ 
   # drift-kinetic interpolation parameter.  0 for full particle
   # dynamics.  1 for interpolation.  2 for pure drift kinetics
a0_s(ns)       _real    [m]     /LARGEPOS/
   # Initial beam width in x by species: if LARGEPOS, set to a0
   # ;if a0 is zero, it is set to a0_s(1)
ap0_s(ns)      _real    [rad]   /LARGEPOS/
   # Initial beam envelope vx/vz by species: if LARGEPOS, set to ap0
   # ;if ap0 is zero, it is set to ap0_s(1)
b0_s(ns)       _real    [m]     /LARGEPOS/
   # Initial beam width in y by species: if LARGEPOS, set to b0
   # ;if b0 is zero, it is set to b0_s(1)
bp0_s(ns)      _real    [rad]   /LARGEPOS/
   # Initial beam envelope vy/vz by species: if LARGEPOS, set to bp0
   # ;if bp0 is zero, it is set to bp0_s(1)
x0_s(ns)       _real    [m]     /LARGEPOS/
   # Initial beam centroid in x by species: if LARGEPOS, set to x0
   # ;if x0 is zero, it is set to x0_s(1)
xp0_s(ns)      _real    [rad]   /LARGEPOS/
   # Initial beam centroid vx/vz by species: if LARGEPOS, set to xp0
   # ;if xp0 is zero, it is set to xp0_s(1)
y0_s(ns)       _real    [m]     /LARGEPOS/
   # Initial beam centroid in y by species: if LARGEPOS, set to y0
   # ;if y0 is zero, it is set to y0_s(1)
yp0_s(ns)      _real    [rad]   /LARGEPOS/
   # Initial beam centroid vy/vz by species: if LARGEPOS, set to yp0
   # ;if yp0 is zero, it is set to yp0_s(1)
aion_s(ns)     _real    [1]     /LARGEPOS/
   # Atomic number by species: if LARGEPOS, set to aion
   # ;if aion is zero, it is set to aion_s(1)
ekin_s(ns)     _real    [V]     /LARGEPOS/
   # Particle energy by species: if LARGEPOS, set to ekin
   # ;if ekin is zero, it is set to ekin_s(1)
emit_s(ns)     _real    [m-rad] /LARGEPOS/
   # Emittance by species: if LARGEPOS, set to emit
   # ;if emit is zero, it is set to emit_s(1)
emitx_s(ns)    _real    [m-rad] /LARGEPOS/
   # X-Emittance by species: if LARGEPOS, set to xemit
   # ;if emitx is zero, it is set to emitx_s(1)
emity_s(ns)    _real    [m-rad] /LARGEPOS/
   # Y-Emittance by species: if LARGEPOS, set to yemit
   # ;if emit is zero, it is set to emity_s(1)
emitn_s(ns)    _real    [m-rad] /LARGEPOS/
   # Normalized emittance by species: if LARGEPOS, set to emitn
   # ;if emitn is zero, it is set to emitn_s(1)
emitnx_s(ns)   _real    [m-rad] /LARGEPOS/
   # Normalized X-emittance by species: if LARGEPOS, set to emitnx
   # ; if emitnx is zero, it is set to emitnx_s(1)
emitny_s(ns)   _real    [m-rad] /LARGEPOS/
   # Normalized Y-emittance by species: if LARGEPOS, set to emitny
   # ; if emitny is zero, it is set to emitny_s(1)
ibeam_s(ns)    _real    [A]     /LARGEPOS/
   # Beam current by species: if LARGEPOS, set to ibeam
   # ;if ibeam is zero, it is set to ibeam_s(1)
zion_s(ns)     _real    [1]     /LARGEPOS/
   # Charge state by species: if LARGEPOS, set to zion
   # ;if zion is zero, it is set to zion_s(1)
straight_s(ns) _real    [1]     /LARGEPOS/
   # Percent of beam that isn't cigar: if LARGEPOS, set to straight
   # ;if straight is zero, it is set to straight_s
vbeam_s(ns)    _real    [m/s]   /LARGEPOS/
   # Particle axial velocity by species: if LARGEPOS, set to vbeam
   # ;if vbeam is zero, it is set to vbeam_s(1)
vtilt_s(ns)    _real    [1]     /LARGEPOS/
   # Axial velocity tilt by species: if LARGEPOS, set to vtilt
   # ;if vtilt is zero, it is set to vtilt_s(1)
vthperp_s(ns)  _real    [m/s]   /LARGEPOS/
   # Transverse thermal spread by species: if LARGEPOS, set to vthperp
   # ;if vthperp is zero, it is set to vthperp_s(1)
vthz_s(ns)     _real    [m/s]   /LARGEPOS/
   # Axial thermal spread by species: if LARGEPOS, set to vthz
   # ;if vthz is zero, it is set to vthz_s(1)
zimin_s(ns)    _real    [m]     /LARGEPOS/
   # Minimum initial z of beam by species: if LARGEPOS, set to zimin
   # ;if zimin is zero, it is set to zimin_s(1)
zimax_s(ns)    _real    [m]     /LARGEPOS/
   # Maximum initial z of beam by species: if LARGEPOS, set to zimax
   # ;if zimax is zero, it is set to zimax_s(1)
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

efetch(ns)     _integer /1/
   # Specifies type of fetch used to set internal E field
   # 1: direct calculation from phi
   # 2: direct calculation from phi, vectorization over 8 grid cell corners
   #    (not recommended)
   # 3: indirect calculation from pre-calculated E (finite differences of phi)
   # 4: energy conserving - uses nearest grid point interpolation
   # 5: same as 1, but checks for particles out of bounds
   # Method 3 is generally fastest but requires lots of extra storage space.
   # Next best is 1.

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
prwelip         real   [1]    /1./
   # Ellipticity of cylindrical wall (ay = prwelip*ax)
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
xpinject(ninject) _real    # Vx/Vz of injected particles
ypinject(ninject) _real    # Vy/Vz of injected particles
thetainject(ninject) _real # theta angle of injection source
phiinject(ninject)   _real # phi angle of injection source
vzinject(ninject,ns) _real /0./ [m/s] # Starting velocity.
                                      # For inject == 1, autoset to vbeam
finject(ninject,ns) _real  # Species fraction for each source
winject(ninject,ns) _real /1./ # Scale factor on the particle weight when
                               # weighted particles are used (when wpid > 0)
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
rnpinje_s(ns) _real        # Number of particles injected each step
                           # for each species (do not need to be integers)
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
linj_efromgrid logical /.false./ # When true, the special calculation of the
                                 # E field from the virtual surface is skipped
                                 # and the field from the full field grid is
                                 # used. It is not recommended that this flag
                                 # be set.
linj_rectangle logical /.false./ # Flags whether injection is over a
                                 # rectangular region and not elliptical  
linj_includebz0 logical /.false./ # Include the solenoidal field bz0, giving
                                  # the particles a vtheta corresponding to the
                                  # Larmor frequency

ntinj          integer /0/ # Number of transverse emitting surfaces.
ztinjmn(ntinj)    _real # Minimum z extent of transverse injection.
ztinjmx(ntinj)    _real # Maximum z extent of transverse injection.
dztinj(ntinj)     _real    # Step size in z of the transverse injection arrays
nztinj(ntinj) _integer     # Number of longitudinal points for transverse injection
nztmax         integer /0/ # Maximum length (in z grid cells) of
                           # transverse injection surfaces
atinject(0:nztmax,ntinj) _real /-1./
                           # X Radius as a function of z of transverse
                           # emitting surface.
btinject(0:nztmax,ntinj) _real /-1./
                           # Y Radius as a function of z of transverse
                           # emitting surface.
ltinj_outward(ntinj) _logical /1/ # When true, the direction of injection is
                                  # outward, otherwise inward.
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
wtinject(ntinj,ns) _real /1./ # Scale factor on the particle weight when
                              # weighted particles are used (when wpid > 0)
ntinject(ns)      _integer # Number of particles injected off of the
                           # transverse injection sources.
jmaxtinj(ntinj)  _real /LARGEPOS/ # Maximum injected current for transverse injection.
tinj_phi(0:nttinjmax,0:nztmax,ntinj) _real
   # Grid holding the potential drop for transverse emitting surfaces
tinjprev(0:nttinjmax,0:nztmax,ntinj) _real
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
dzz                     real  [m]      # Z_arrays mesh grid cell size
dzzi                    real  [m]      # Z_arrays mesh grid cell size inverse
zzmin                   real  [m]      # Z_arrays mesh maximum in z
zzmax                   real  [m]      # Z_arrays mesh minimum in z
nzzarr            integer /0/ [1]      # Length of arrays in group Z_arrays
nszarr            integer /0/          # Number of species z data is calculated
                                       # for. Defaults to zero, unless
                                       # lspeciesmoments is true, then it
                                       # defaults to top.ns.
                                       # Should always be same as nszmmnt.
zplmesh(0:nzzarr)      _real  [m]      # Z mesh used for qtys in group Z_arrays
egap(0:nzzarr)         _real  [V/m]    # Gap electric field (smeared in z)
linechg(0:nzzarr)      _real  [C/m]    # Line charge density
vzofz(0:nzzarr)        _real  [m/s]    # Mean axial speed vs z
rhoax(0:nzzarr)        _real  [C/m^3]  # charge density on axis
phiax(0:nzzarr)        _real  [V]      # potential on axis
ezax(0:nzzarr)         _real  [V/m]    # space charge E field on axis
eearsofz(0:nzzarr)     _real  [V/m]    # confining Eears, as a function of z
prwallz(0:nzzarr)      _real  [m]      # Radius at which particles are absorbed
prwallxz(0:nzzarr)     _real  [m]      # X of center of cylindrical wall
prwallyz(0:nzzarr)     _real  [m]      # Y of center of cylindrical wall
prwelipz(0:nzzarr)     _real  [1] /1./ # Ellipticity of cylindrical wall
                                       # (ay = prwelipz*ax)
lamkreal(0:nzzarr)     _real  [C/m]    # Real part of FFT of lambda
lamkimag(0:nzzarr)     _real  [C/m]    # Imaginary part of FFT of lambda
curr(0:nzzarr,0:nszarr)     _real [A]  # Beam current
lostpars(0:nzzarr,0:nszarr) _integer   # number of lost particles by zcells


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
nzwind       integer /0/ # Actual number of z windows used for emittance calcs
nswind       integer /0/ # Number of species window moments data is calculated
                         # for. Defaults to zero, unless lspeciesmoments is
                         # true, then it defaults to top.ns.
                         # Should always be same as nszmmnt.
pnum(0:nzwind,0:nswind)    _real [1]     # Total no. of (physical) ions in window
xbar(0:nzwind,0:nswind)    _real [m]     # Mean X coordinate in window
ybar(0:nzwind,0:nswind)    _real [m]     # Mean Y coordinate in window
zbar(0:nzwind,0:nswind)    _real [m]     # Mean axial location in window
xpbar(0:nzwind,0:nswind)   _real [1]     # Mean X' in window
ypbar(0:nzwind,0:nswind)   _real [1]     # Mean Y' in window
vxbar(0:nzwind,0:nswind)   _real [m/s]   # Mean Vx in window
vybar(0:nzwind,0:nswind)   _real [m/s]   # Mean Vy in window
vzbar(0:nzwind,0:nswind)   _real [m/s]   # Mean Vz in window
xybar(0:nzwind,0:nswind)   _real [m**2]  # Mean product of X and Y in window
xypbar(0:nzwind,0:nswind)  _real [m]     # Mean product of X  and Y'
yxpbar(0:nzwind,0:nswind)  _real [m]     # Mean product of Y  and X'
xpypbar(0:nzwind,0:nswind) _real [1]     # Mean product of X' and Y'
xsqbar(0:nzwind,0:nswind)  _real [m**2]  # Mean X-squared in window
ysqbar(0:nzwind,0:nswind)  _real [m**2]  # Mean Y-squared in window
zsqbar(0:nzwind,0:nswind)  _real [m**2]  # Mean Z-squared in window
xpsqbar(0:nzwind,0:nswind) _real [1]     # Mean X' squared in window
ypsqbar(0:nzwind,0:nswind) _real [1]     # Mean Y' squared in window
vxsqbar(0:nzwind,0:nswind) _real [m/s]   # Mean Vx squared in window
vysqbar(0:nzwind,0:nswind) _real [m/s]   # Mean Vy squared in window
vzsqbar(0:nzwind,0:nswind) _real [m/s]   # Mean Vz squared in window
xxpbar(0:nzwind,0:nswind)  _real [m]     # Mean product of X and X' in window
yypbar(0:nzwind,0:nswind)  _real [m]     # Mean product of Y and Y' in window
zvzbar(0:nzwind,0:nswind)  _real [m]     # Mean product of Z and Vz in window
xvzbar(0:nzwind,0:nswind)  _real [m]     # Mean product of X and Vz in window
yvzbar(0:nzwind,0:nswind)  _real [m]     # Mean product of X and Vz in window
vxvzbar(0:nzwind,0:nswind) _real [m]     # Mean product of Vx and Vz in window
vyvzbar(0:nzwind,0:nswind) _real [m]     # Mean product of Vy and Vz in window
xrms(0:nzwind,0:nswind)    _real [m]     # RMS X in window
yrms(0:nzwind,0:nswind)    _real [m]     # RMS Y in window
zrms(0:nzwind,0:nswind)    _real [m]     # RMS Z in window
rrms(0:nzwind,0:nswind)    _real [m]     # RMS R in window
xprms(0:nzwind,0:nswind)   _real [m]     # RMS X' in window
yprms(0:nzwind,0:nswind)   _real [m]     # RMS Y' in window
epsx(0:nzwind,0:nswind)    _real [m-rad] # X-X' emittance
epsy(0:nzwind,0:nswind)    _real [m-rad] # Y-Y' emittance
epsz(0:nzwind,0:nswind)    _real [m-rad] # Z-Z' emittance
epsnx(0:nzwind,0:nswind)   _real [mm-mrad] # X-X' normalized emittance
epsny(0:nzwind,0:nswind)   _real [mm-mrad] # Y-Y' normalized emittance
epsnz(0:nzwind,0:nswind)   _real [mm-mrad] # Z-Z' normalized emittance
epsg(0:nzwind,0:nswind)    _real [m-rad] # Generalized emittance
epsh(0:nzwind,0:nswind)    _real [m-rad] # Generalized emittance
epsng(0:nzwind,0:nswind)   _real [mm-mrad] # Generalized normalized emittance
epsnh(0:nzwind,0:nswind)   _real [mm-mrad] # Generalized normalized emittance
vxrms(0:nzwind,0:nswind)   _real [m/s]   # True RMS Vx in window
vyrms(0:nzwind,0:nswind)   _real [m/s]   # True RMS Vy in window
vzrms(0:nzwind,0:nswind)   _real [m/s]   # True RMS Vz in window
rhomid(0:nzwind)           _real [C/m^3] # Charge dens. on axis at ctr of window
rhomax(0:nzwind)           _real [C/m^3] # Charge dens. max-over-X,Y at ctr of win.

*********** Z_Moments dump:
# Particle and field moment data (including emittances) at current timestep
# as a function of Z 
zmmntmax         real         # Moments grid maximum in Z
zmmntmin         real         # Moments grid minimum in Z
nzmmnt           integer /0/  # Number of points in z moments grid
nszmmnt          integer /0/  # Number of species z moments data is
                              # calculated for. Defaults to zero, unless
                              # lspeciesmoments is true, then it defaults
                              # to top.ns.
dzm              real         # Moments grid cell size
dzmi             real         # Moments grid cell size inverse
numzmmnt         integer /NUMZMMNT/ # Number of moments calculated
zmmntdtextmax    real /LARGEPOS/ # Cutoff of time step for extrapolation of
                                 # particles, in units of top.dt.
zmntmesh(0:nzmmnt)  _real [m]    # Z mesh associated with Z moments
zmmntsq(0:nszmmnt)  _real [kg]   # Particle charge of species associated with
                                 # the moments
zmmntsm(0:nszmmnt)  _real [kg]   # Particle mass of species associated with
                                 # the moments
zmmntsw(0:nszmmnt)  _real [kg]   # Particle weight of species associated with
                                 # the moments
zmomentscalculated(0:nszmmnt) _logical     # Set to true if the moments for
                                           # that species has been calculated
pnumz(0:nzmmnt,0:nszmmnt)    _real [1]     # No. of (physical) ions at grid point
xbarz(0:nzmmnt,0:nszmmnt)    _real [m]     # Mean X coordinate at grid point
ybarz(0:nzmmnt,0:nszmmnt)    _real [m]     # Mean Y coordinate at grid point
zbarz(0:nzmmnt,0:nszmmnt)    _real [m]     # Mean axial location at grid point
xpbarz(0:nzmmnt,0:nszmmnt)   _real [1]     # Mean X' at grid point
ypbarz(0:nzmmnt,0:nszmmnt)   _real [1]     # Mean Y' at grid point
vxbarz(0:nzmmnt,0:nszmmnt)   _real [m/s]   # Mean Vx at grid point
vybarz(0:nzmmnt,0:nszmmnt)   _real [m/s]   # Mean Vy at grid point
vzbarz(0:nzmmnt,0:nszmmnt)   _real [m/s]   # Mean Vz at grid point
xybarz(0:nzmmnt,0:nszmmnt)   _real [m**2]  # Mean product of X  and Y  at grid point
xypbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of X  and Y' at grid point
yxpbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of Y  and X' at grid point
xpypbarz(0:nzmmnt,0:nszmmnt) _real [1]     # Mean product of X' and Y' at grid point
xsqbarz(0:nzmmnt,0:nszmmnt)  _real [m**2]  # Mean X-squared at grid point
ysqbarz(0:nzmmnt,0:nszmmnt)  _real [m**2]  # Mean Y-squared at grid point
zsqbarz(0:nzmmnt,0:nszmmnt)  _real [m**2]  # Mean Z-squared at grid point
xpsqbarz(0:nzmmnt,0:nszmmnt) _real [1]     # Mean X' squared at grid point
ypsqbarz(0:nzmmnt,0:nszmmnt) _real [1]     # Mean Y' squared at grid point
vxsqbarz(0:nzmmnt,0:nszmmnt) _real [m/s]   # Mean Vx squared at grid point
vysqbarz(0:nzmmnt,0:nszmmnt) _real [m/s]   # Mean Vy squared at grid point
vzsqbarz(0:nzmmnt,0:nszmmnt) _real [m/s]   # Mean Vz squared at grid point
xxpbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of X and X' at grid point
yypbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of Y and Y' at grid point
zvzbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of Z and Vz at grid point
xvzbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of X and Vz at grid point
yvzbarz(0:nzmmnt,0:nszmmnt)  _real [m]     # Mean product of X and Vz at grid point
vxvzbarz(0:nzmmnt,0:nszmmnt) _real [m]     # Mean product of Vx and Vz at grid point
vyvzbarz(0:nzmmnt,0:nszmmnt) _real [m]     # Mean product of Vy and Vz at grid point
xrmsz(0:nzmmnt,0:nszmmnt)    _real [m]     # RMS X at grid point
yrmsz(0:nzmmnt,0:nszmmnt)    _real [m]     # RMS Y at grid point
zrmsz(0:nzmmnt,0:nszmmnt)    _real [m]     # RMS Z at grid point
rrmsz(0:nzmmnt,0:nszmmnt)    _real [m]     # RMS R at grid point
xprmsz(0:nzmmnt,0:nszmmnt)   _real [m]     # RMS X' at grid point
yprmsz(0:nzmmnt,0:nszmmnt)   _real [m]     # RMS Y' at grid point
epsxz(0:nzmmnt,0:nszmmnt)    _real [m-rad] # X-X' emittance at grid point
epsyz(0:nzmmnt,0:nszmmnt)    _real [m-rad] # Y-Y' emittance at grid point
epszz(0:nzmmnt,0:nszmmnt)    _real [m-rad] # Z-Z' emittance at grid point
epsnxz(0:nzmmnt,0:nszmmnt)   _real [mm-mrad] # X-X' normalized emittance at grid point
epsnyz(0:nzmmnt,0:nszmmnt)   _real [mm-mrad] # Y-Y' normalized emittance at grid point
epsnzz(0:nzmmnt,0:nszmmnt)   _real [mm-mrad] # Z-Z' normalized emittance at grid point
epsgz(0:nzmmnt,0:nszmmnt)    _real [m-rad]   # Generalized emittance on grid
epshz(0:nzmmnt,0:nszmmnt)    _real [m-rad]   # Generalized emittance on grid
epsngz(0:nzmmnt,0:nszmmnt)   _real [mm-mrad] # Generalized normalized emittance on grid
epsnhz(0:nzmmnt,0:nszmmnt)   _real [mm-mrad] # Generalized normalized emittance on grid
vxrmsz(0:nzmmnt,0:nszmmnt)   _real [m/s]   # True RMS Vx at grid point
vyrmsz(0:nzmmnt,0:nszmmnt)   _real [m/s]   # True RMS Vy at grid point
vzrmsz(0:nzmmnt,0:nszmmnt)   _real [m/s]   # True RMS Vz at grid point
rhomidz(0:nzmmnt)            _real [C/m^3] # Charge dens. on axis at grid point
rhomaxz(0:nzmmnt)            _real [C/m^3] # Charge dens. max-over-X,Y at grid point
tempmaxp(6,0:nszmmnt)                   _real # Temporary work array
tempminp(6,0:nszmmnt)                   _real # Temporary work array
tempzmmnts0(NUMZMMNT,0:nszmmnt)         _real # Temporary work array
tempzmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt) _real # Temporary work array

********** Lab_Moments dump:
# Particle moment data as a function of time at locations in the lab frame
nlabwn  integer   /50/ # number of lab windows
nslabwn integer   /0/  # Number of species lab moments data is
                       # calculated for. Defaults to zero, unless
                       # lspeciesmoments is true, then it defaults to top.ns.
                       # Should always be same as nszmmnt.
zlw(nlabwn) _real [m] /LARGEPOS/ # z for lab windows
iflabwn integer /1/ # turns on lab window moments (0 off; 1 on)
itlabwn integer /0/ # Sets how often the lab moments are calculated
ntlabwn integer     # Maximum number of times lab frame moments are calculated
ilabwn(nlabwn,0:nslabwn) _integer # Number of times lab frame moments have been calculated
timelw(ntlabwn,nlabwn,0:nslabwn)     _real # Time in lab frame
pnumlw(ntlabwn,nlabwn,0:nslabwn)     _real # Number of particles in lab frame
xbarlw(ntlabwn,nlabwn,0:nslabwn)     _real # X bar in lab frame
ybarlw(ntlabwn,nlabwn,0:nslabwn)     _real # Y bar in lab frame
vzbarlw(ntlabwn,nlabwn,0:nslabwn)    _real # Vz bar in lab frame
epsxlw(ntlabwn,nlabwn,0:nslabwn)     _real # X emittance in lab frame
epsylw(ntlabwn,nlabwn,0:nslabwn)     _real # Y emittance in lab frame
epszlw(ntlabwn,nlabwn,0:nslabwn)     _real # Z emittance in lab frame
vxrmslw(ntlabwn,nlabwn,0:nslabwn)    _real # Vx RMS in lab frame
vyrmslw(ntlabwn,nlabwn,0:nslabwn)    _real # Vy RMS in lab frame
vzrmslw(ntlabwn,nlabwn,0:nslabwn)    _real # Vz RMS in lab frame
xrmslw(ntlabwn,nlabwn,0:nslabwn)     _real # X RMS in lab frame
yrmslw(ntlabwn,nlabwn,0:nslabwn)     _real # Y RMS in lab frame
rrmslw(ntlabwn,nlabwn,0:nslabwn)     _real # R RMS in lab frame
xxpbarlw(ntlabwn,nlabwn,0:nslabwn)   _real # XX' bar in lab frame
yypbarlw(ntlabwn,nlabwn,0:nslabwn)   _real # YY' bar in lab frame
currlw(ntlabwn,nlabwn,0:nslabwn)     _real # Current in lab frame
lostparslw(ntlabwn,nlabwn,0:nslabwn) _real # Number of lost particles in lab frame
linechglw(ntlabwn,nlabwn,0:nslabwn)  _real # Line-charge in lab frame

*********** Moments dump:
# Scalar moments of general interest
nsmmnt integer /0/  # Number of species moments data is
                    # calculated for. Defaults to zero, unless
                    # lspeciesmoments is true, then it defaults to top.ns.
                    # Should always be same as nszmmnt.
ese                real [J]      # Electrostatic energy
ek(0:nsmmnt)      _real [J]      # Total Kinetic energy
ekzmbe(0:nsmmnt)  _real [J]      # Total Z Kinetic energy minus beam energy
                      # 1/2m(vzsqbar - vbeam**2)
ekzbeam(0:nsmmnt) _real [J]      # Z Kinetic energy in the beam frame
                      # 1/2m ave[(vz - vbeam)**2]
ekperp(0:nsmmnt)  _real [J]      # Perp Kinetic energy
pz(0:nsmmnt)      _real [kg-m/s] # Total axial momentum (subtracting out Vbeam)
xmaxp(0:nsmmnt)   _real [m]    # Maximum X  over particles (set intermittently)
xminp(0:nsmmnt)   _real [m]    # Minimum X  over particles (set intermittently)
ymaxp(0:nsmmnt)   _real [m]    # Maximum Y  over particles (set intermittently)
yminp(0:nsmmnt)   _real [m]    # Minimum Y  over particles (set intermittently)
zmaxp(0:nsmmnt)   _real [m]    # Maximum Z  over particles (set intermittently)
zminp(0:nsmmnt)   _real [m]    # Minimum Z  over particles (set intermittently)
vxmaxp(0:nsmmnt)  _real [m/s]  # Maximum Vx over particles (set intermittently)
vxminp(0:nsmmnt)  _real [m/s]  # Minimum Vx over particles (set intermittently)
vymaxp(0:nsmmnt)  _real [m/s]  # Maximum Vy over particles (set intermittently)
vyminp(0:nsmmnt)  _real [m/s]  # Minimum Vy over particles (set intermittently)
vzmaxp(0:nsmmnt)  _real [m/s]  # Maximum Vz over particles (set intermittently)
vzminp(0:nsmmnt)  _real [m/s]  # Minimum Vz over particles (set intermittently)

*********** Hist dump history:
# History data
jhist                          integer  /-1/
   # pointer to current entry in history arrays
lenhist                        integer  /0/
   # length of the history arrays
nshist integer /0/  # Number of species history data is save.
                    # Defaults to zero, unless lspeciesmoments is true, then
                    # it defaults to top.ns.
                    # Should always be same as nszmmnt.
thist(0:lenhist)              _real [s]     limited (0:jhist)
   # Times at which data is saved
hzbeam(0:lenhist)             _real [m]     limited (0:jhist)
   # Beam frame location in lab frame
hvbeam(0:lenhist)             _real [m/s]   limited (0:jhist)
   # Beam frame velocity
hbmlen(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # RMS beam length
hefld(0:lenhist)              _real [J]     limited (0:jhist)
   # Field energy
hekzmbe(0:lenhist,0:nshist)   _real [J]     limited (0:jhist,0:nshist)
   # Total Z Kinetic energy minus beam energy
hekzbeam(0:lenhist,0:nshist)  _real [J]     limited (0:jhist,0:nshist)
   # Z Kinetic energy in the beam frame
hekperp(0:lenhist,0:nshist)   _real [J]     limited (0:jhist,0:nshist)
   # Perp Kinetic energy
hxmaxp(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # History of maximum X over particles
hxminp(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # History of minimum X over particles
hymaxp(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # History of maximum Y over particles
hyminp(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # History of minimum Y over particles
hzmaxp(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # History of maximum Z over particles
hzminp(0:lenhist,0:nshist)    _real [m]     limited (0:jhist,0:nshist)
   # History of minimum Z over particles
hvxmaxp(0:lenhist,0:nshist)   _real [m/s]   limited (0:jhist,0:nshist)
   # History of maximum Vx over particles
hvxminp(0:lenhist,0:nshist)   _real [m/s]   limited (0:jhist,0:nshist)
   # History of minimum Vx over particles
hvymaxp(0:lenhist,0:nshist)   _real [m/s]   limited (0:jhist,0:nshist)
   # History of maximum Vy over particles
hvyminp(0:lenhist,0:nshist)   _real [m/s]   limited (0:jhist,0:nshist)
   # History of minimum Vy over particles
hvzmaxp(0:lenhist,0:nshist)   _real [m/s]   limited (0:jhist,0:nshist)
   # History of maximum Vz over particles
hvzminp(0:lenhist,0:nshist)   _real [m/s]   limited (0:jhist,0:nshist)
   # History of minimum Vz over particles
hepsx(0:nzwind,0:lenhist,0:nshist)     _real [m-r]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # X-X' emittance by window as a function of time
hepsy(0:nzwind,0:lenhist,0:nshist)     _real [m-r]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Y-Y' emittance by window as a function of time
hepsz(0:nzwind,0:lenhist,0:nshist)     _real [m-r]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Z-Z' emittance by window as a function of time
hepsnx(0:nzwind,0:lenhist,0:nshist)    _real [mm-mr]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # X-X' normalized emittance by window as a function of time
hepsny(0:nzwind,0:lenhist,0:nshist)    _real [mm-mr]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Y-Y' normalized emittance by window as a function of time
hepsnz(0:nzwind,0:lenhist,0:nshist)    _real [mm-mr]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Z-Z' normalized emittance by window as a function of time
hepsg(0:nzwind,0:lenhist,0:nshist)     _real [m-r]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Generalized emittance by window as a function of time
hepsh(0:nzwind,0:lenhist,0:nshist)     _real [m-r]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Generalized emittance by window as a function of time
hepsng(0:nzwind,0:lenhist,0:nshist)    _real [mm-mr]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Generalized normalized emittance by window as a function of time
hepsnh(0:nzwind,0:lenhist,0:nshist)    _real [mm-mr]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Generalized normalized emittance by window as a function of time
hpnum(0:nzwind,0:lenhist,0:nshist)     _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Number of particles in each z window (species 1,0:nshist)
hrhomid(0:nzwind,0:lenhist)            _real [C/m^3]
   limited (0:nzwind,0:jhist)          +winhist
   # Charge density on axis at center of z window as a fcn of time
hrhomax(0:nzwind,0:lenhist)            _real [C/m^3]
   limited (0:nzwind,0:jhist)          +winhist
   # Charge dens. max-over-x,y at ctr of z window as a fcn of time
hxbar(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True mean x by window as a function of time
hybar(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True mean y by window as a function of time
hzbar(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True mean z by window as a function of time
hxybar(0:nzwind,0:lenhist,0:nshist)    _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True mean xy by window as a function of time
hxrms(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS x by window as a function of time
hyrms(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS y by window as a function of time
hrrms(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS r by window as a function of time
hzrms(0:nzwind,0:lenhist,0:nshist)     _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS z by window as a function of time
hxprms(0:nzwind,0:lenhist,0:nshist)    _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS x' by window as a function of time
hyprms(0:nzwind,0:lenhist,0:nshist)    _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS y' by window as a function of time
hxsqbar(0:nzwind,0:lenhist,0:nshist)   _real [m^2]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # x squared bar by window as a function of time
hysqbar(0:nzwind,0:lenhist,0:nshist)   _real [m^2]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # y squared bar by window as a function of time
hvxbar(0:nzwind,0:lenhist,0:nshist)    _real [m/s]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean vx by window as a function of time
hvybar(0:nzwind,0:lenhist,0:nshist)    _real [m/s]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean vy by window as a function of time
hvzbar(0:nzwind,0:lenhist,0:nshist)    _real [m/s]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean vz by window as a function of time
hxpbar(0:nzwind,0:lenhist,0:nshist)    _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean x' by window as a function of time
hypbar(0:nzwind,0:lenhist,0:nshist)    _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean y' by window as a function of time
hvxrms(0:nzwind,0:lenhist,0:nshist)    _real [m/s]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS vx by window as a function of time
hvyrms(0:nzwind,0:lenhist,0:nshist)    _real [m/s]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS vy by window as a function of time
hvzrms(0:nzwind,0:lenhist,0:nshist)    _real [m/s]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # True RMS vz by window as a function of time
hxpsqbar(0:nzwind,0:lenhist,0:nshist)  _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean x' squared by window as a function of time
hypsqbar(0:nzwind,0:lenhist,0:nshist)  _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean y' squared by window as a function of time
hxxpbar(0:nzwind,0:lenhist,0:nshist)   _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean x * x' by window as a function of time
hyypbar(0:nzwind,0:lenhist,0:nshist)   _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean y  * y' by window as a function of time
hxypbar(0:nzwind,0:lenhist,0:nshist)   _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean x  * y' by window as a function of time 
hyxpbar(0:nzwind,0:lenhist,0:nshist)   _real [m]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean y  * x' by window as a function of time 
hxpypbar(0:nzwind,0:lenhist,0:nshist)  _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean x' * y' by window as a function of time
hxvzbar(0:nzwind,0:lenhist,0:nshist)   _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean x * vz by window as a function of time
hyvzbar(0:nzwind,0:lenhist,0:nshist)   _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean y * vz by window as a function of time
hvxvzbar(0:nzwind,0:lenhist,0:nshist)  _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean vx * vz by window as a function of time
hvyvzbar(0:nzwind,0:lenhist,0:nshist)  _real [1]
   limited (0:nzwind,0:jhist,0:nshist) +winhist
   # Mean vy * vz by window as a function of time
lhlinechg logical /.true./   # Turns on history of line charge
ihlinechg integer /1/        # Multiplier for hlinechg memory size (autoset)
hlinechg(0:nzzarr*ihlinechg,0:lenhist) _real [C/m]
            limited (0:nzzarr,0:jhist)
            +zhist           # Line charge density vs. space and time
lhvzofz logical /.true./     # Turns on history of vz
ihvzofz integer /1/          # Multiplier for hvzofz memory size (autoset)
hvzofz(0:nzzarr*ihvzofz,0:lenhist)  _real [m/s]
            limited (0:nzzarr,0:jhist)
            +zhist           # Vz versus space and time
lhcurrz logical /.false./    # Turns on history of current
ihcurrz integer /1/          # Multiplier for hcurrz memory size (autoset)
hcurrz(0:nzzarr*ihcurrz,0:lenhist,0:nszarr)  _real [m/s]
            limited (0:nzzarr,0:jhist,0:nszarr)
            +zhist           # Current versus space and time
lhpnumz logical /.false./    # Turns on history of particle number
ihpnumz integer /0 /         # Multiplier for hpnumz memory size (autoset)
hpnumz(0:nzmmnt*ihpnumz,0:lenhist,0:nshist)  _real [m-r]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Simulation particle number
lhepsxz logical /.false./    # Turns on history of X emittance
ihepsxz integer /0 /         # Multiplier for hepsxz memory size (autoset)
hepsxz(0:nzmmnt*ihepsxz,0:lenhist,0:nshist)  _real [m-r]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X emittance versus space and time
lhepsyz logical /.false./    # Turns on history of Y emittance
ihepsyz integer /0 /         # Multiplier for hepsyz memory size (autoset)
hepsyz(0:nzmmnt*ihepsyz,0:lenhist,0:nshist)  _real [m-r]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y emittance versus space and time
lhepsnxz logical /.false./   # Turns on history of X normalized emittance
ihepsnxz integer /0 /        # Multiplier for hepsnxz memory size (autoset)
hepsnxz(0:nzmmnt*ihepsnxz,0:lenhist,0:nshist)  _real [mm-mrad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X normalized emittance versus space and time
lhepsnyz logical /.false./   # Turns on history of Y normalized emittance
ihepsnyz integer /0 /        # Multiplier for hepsnyz memory size (autoset)
hepsnyz(0:nzmmnt*ihepsnyz,0:lenhist,0:nshist)  _real [mm-mrad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y normalized emittance versus space and time
lhepsgz logical /.false./    # Turns on history of Generalized emittance
ihepsgz integer /0 /         # Multiplier for hepsgz memory size (autoset)
hepsgz(0:nzmmnt*ihepsgz,0:lenhist,0:nshist)  _real [m-rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Generalized emittance versus space and time
lhepshz logical /.false./    # Turns on history of Generalized emittance
ihepshz integer /0 /         # Multiplier for hepshz memory size (autoset)
hepshz(0:nzmmnt*ihepshz,0:lenhist,0:nshist)  _real [m-rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Generalized emittance versus space and time
lhepsngz logical /.false./   # Turns on history of Generalized nrmlzd emittance
ihepsngz integer /0 /        # Multiplier for hepsngz memory size (autoset)
hepsngz(0:nzmmnt*ihepsngz,0:lenhist,0:nshist)  _real [mm-mrad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Generalized nrmlzd emittance versus space andtime
lhepsnhz logical /.false./   # Turns on history of Generalized nrmlzd emittance
ihepsnhz integer /0 /        # Multiplier for hepsnhz memory size (autoset)
hepsnhz(0:nzmmnt*ihepsnhz,0:lenhist,0:nshist)  _real [mm-mrad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Generalized nrmlzd emittance versus space andtime
lhxbarz logical /.false./    # Turns on history of X bar
ihxbarz integer /0 /         # Multiplier for hxbarz memory size (autoset)
hxbarz(0:nzmmnt*ihxbarz,0:lenhist,0:nshist)  _real [m]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X bar versus space and time
lhybarz logical /.false./    # Turns on history of Y bar
ihybarz integer /0 /         # Multiplier for hybarz memory size (autoset)
hybarz(0:nzmmnt*ihybarz,0:lenhist,0:nshist)  _real [m]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y bar versus space and time
lhxybarz logical /.false./   # Turns on history of XY bar
ihxybarz integer /0 /        # Multiplier for hxybarz memory size (autoset)
hxybarz(0:nzmmnt*ihxybarz,0:lenhist,0:nshist)  _real [m**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # XY bar versus space and time
lhxrmsz logical /.false./    # Turns on history of X rms
ihxrmsz integer /0 /         # Multiplier for hxrmsz memory size (autoset)
hxrmsz(0:nzmmnt*ihxrmsz,0:lenhist,0:nshist)  _real [m]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X rms versus space and time
lhyrmsz logical /.false./    # Turns on history of Y rms
ihyrmsz integer /0 /         # Multiplier for hyrmsz memory size (autoset)
hyrmsz(0:nzmmnt*ihyrmsz,0:lenhist,0:nshist)  _real [m]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y rms versus space and time
lhrrmsz logical /.false./    # Turns on history of X rms
ihrrmsz integer /0 /         # Multiplier for hrrmsz memory size (autoset)
hrrmsz(0:nzmmnt*ihrrmsz,0:lenhist,0:nshist)  _real [m]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X rms versus space and time
lhxprmsz logical /.false./   # Turns on history of X' rms
ihxprmsz integer /0 /        # Multiplier for hxprmsz memory size (autoset)
hxprmsz(0:nzmmnt*ihxprmsz,0:lenhist,0:nshist)  _real [rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X' rms versus space and time
lhyprmsz logical /.false./   # Turns on history of Y' rms
ihyprmsz integer /0 /        # Multiplier for hyprmsz memory size (autoset)
hyprmsz(0:nzmmnt*ihyprmsz,0:lenhist,0:nshist)  _real [rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y' rms versus space and time
lhxsqbarz logical /.false./  # Turns on history of X**2 bar
ihxsqbarz integer /0 /       # Multiplier for hxsqbarz memory size (autoset)
hxsqbarz(0:nzmmnt*ihxsqbarz,0:lenhist,0:nshist)  _real [m**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X**2 bar versus space and time
lhysqbarz logical /.false./  # Turns on history of Y**2 bar
ihysqbarz integer /0 /       # Multiplier for hysqbarz memory size (autoset)
hysqbarz(0:nzmmnt*ihysqbarz,0:lenhist,0:nshist)  _real [m**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y**2 bar versus space and time
lhvxbarz logical /.false./   # Turns on history of Vx bar
ihvxbarz integer /0 /        # Multiplier for hvxbarz memory size (autoset)
hvxbarz(0:nzmmnt*ihvxbarz,0:lenhist,0:nshist)  _real [m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Vx bar versus space and time
lhvybarz logical /.false./   # Turns on history of Vy bar
ihvybarz integer /0 /        # Multiplier for hvybarz memory size (autoset)
hvybarz(0:nzmmnt*ihvybarz,0:lenhist,0:nshist)  _real [m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Vy bar versus space and time
lhvzbarz logical /.false./   # Turns on history of Vz bar
ihvzbarz integer /0 /        # Multiplier for hvzbarz memory size (autoset)
hvzbarz(0:nzmmnt*ihvzbarz,0:lenhist,0:nshist)  _real [m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Vz bar versus space and time
lhxpbarz logical /.false./   # Turns on history of X' bar
ihxpbarz integer /0 /        # Multiplier for hxpbarz memory size (autoset)
hxpbarz(0:nzmmnt*ihxpbarz,0:lenhist,0:nshist)  _real [rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X' bar versus space and time
lhypbarz logical /.false./   # Turns on history of Y' bar
ihypbarz integer /0 /        # Multiplier for hypbarz memory size (autoset)
hypbarz(0:nzmmnt*ihypbarz,0:lenhist,0:nshist)  _real [rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y' bar versus space and time
lhvxrmsz logical /.false./   # Turns on history of Vx rms
ihvxrmsz integer /0 /        # Multiplier for hvxrmsz memory size (autoset)
hvxrmsz(0:nzmmnt*ihvxrmsz,0:lenhist,0:nshist)  _real [m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Vx rms versus space and time
lhvyrmsz logical /.false./   # Turns on history of Vy rms
ihvyrmsz integer /0 /        # Multiplier for hvyrmsz memory size (autoset)
hvyrmsz(0:nzmmnt*ihvyrmsz,0:lenhist,0:nshist)  _real [m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Vy rms versus space and time
lhvzrmsz logical /.false./   # Turns on history of Vz rms
ihvzrmsz integer /0 /        # Multiplier for hvzrmsz memory size (autoset)
hvzrmsz(0:nzmmnt*ihvzrmsz,0:lenhist,0:nshist)  _real [m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Vz rms versus space and time
lhxpsqbarz logical /.false./ # Turns on history of X'**2 bar
ihxpsqbarz integer /0 /      # Multiplier for hxpsqbarz memory size (autoset)
hxpsqbarz(0:nzmmnt*ihxpsqbarz,0:lenhist,0:nshist)  _real [rad**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X'**2 bar versus space and time
lhypsqbarz logical /.false./ # Turns on history of Y'**2 bar
ihypsqbarz integer /0 /      # Multiplier for hypsqbarz memory size (autoset)
hypsqbarz(0:nzmmnt*ihypsqbarz,0:lenhist,0:nshist)  _real [rad**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # Y'**2 bar versus space and time
lhxxpbarz logical /.false./  # Turns on history of XX' bar
ihxxpbarz integer /0 /       # Multiplier for hxxpbarz memory size (autoset)
hxxpbarz(0:nzmmnt*ihxxpbarz,0:lenhist,0:nshist)  _real [m-rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # XX' bar versus space and time
lhyypbarz logical /.false./  # Turns on history of YY' bar
ihyypbarz integer /0 /       # Multiplier for hyypbarz memory size (autoset)
hyypbarz(0:nzmmnt*ihyypbarz,0:lenhist,0:nshist)  _real [m-rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # YY' bar versus space and time
lhxypbarz logical /.false./  # Turns on history of XY' bar
ihxypbarz integer /0 /       # Multiplier for hxypbarz memory size (autoset)
hxypbarz(0:nzmmnt*ihxypbarz,0:lenhist,0:nshist)  _real [m-rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # XY' bar versus space and time
lhyxpbarz logical /.false./  # Turns on history of YX' bar
ihyxpbarz integer /0 /       # Multiplier for hyxpbarz memory size (autoset)
hyxpbarz(0:nzmmnt*ihyxpbarz,0:lenhist,0:nshist)  _real [m-rad]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # YX' bar versus space and time
lhxpypbarz logical /.false./ # Turns on history of X'Y' bar
ihxpypbarz integer /0 /      # Multiplier for hxpypbarz memory size (autoset)
hxpypbarz(0:nzmmnt*ihxpypbarz,0:lenhist,0:nshist)  _real [rad**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # X'Y' bar versus space and time
lhxvzbarz logical /.false./  # Turns on history of XVz bar
ihxvzbarz integer /0 /       # Multiplier for hxvzbarz memory size (autoset)
hxvzbarz(0:nzmmnt*ihxvzbarz,0:lenhist,0:nshist)  _real [m*m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # XVz bar versus space and time
lhyvzbarz logical /.false./  # Turns on history of YVz bar
ihyvzbarz integer /0 /       # Multiplier for hyvzbarz memory size (autoset)
hyvzbarz(0:nzmmnt*ihyvzbarz,0:lenhist,0:nshist)  _real [m*m/s]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # YVz bar versus space and time
lhvxvzbarz logical /.false./ # Turns on history of VxVz bar
ihvxvzbarz integer /0 /      # Multiplier for hvxvzbarz memory size (autoset)
hvxvzbarz(0:nzmmnt*ihvxvzbarz,0:lenhist,0:nshist)  _real [(m/s)**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # VxVz bar versus space and time
lhvyvzbarz logical /.false./ # Turns on history of VyVz bar
ihvyvzbarz integer /0 /      # Multiplier for hvyvzbarz memory size (autoset)
hvyvzbarz(0:nzmmnt*ihvyvzbarz,0:lenhist,0:nshist)  _real [(m/s)**2]
            limited (0:nzmmnt,0:jhist,0:nshist)
            +zhist           # VyVz bar versus space and time

*********** Particles dump:
# Dynamic particle arrays, and related data
npmax  integer    /0/ +parallel  # Maximum no. of particles
                                 # (user input for some loadings)
pgroup ParticleGroup +parallel # Main group holding the particles
nplive integer    /0/  # No. of "live" particles
npid   integer    /0/  # number of columns for pid.
spid   integer    /0/  # position of particles SSNs in array pid
                       # (FORTRAN indexed: based 1)
wpid   integer    /0/  # position of particle weights in array pid
                       # (FORTRAN indexed: based 1)
tpid   integer    /0/  # position of particle creation time in array pid
                       # (FORTRAN indexed: based 1)
rpid   integer    /0/  # position of particle initial radius in array pid
                       # (FORTRAN indexed: based 1)
xbirthpid integer /0/  # position of particle birth x in array pid
                       # (FORTRAN indexed: based 1)
ybirthpid integer /0/  # position of particle birth y in array pid
                       # (FORTRAN indexed: based 1)
zbirthpid integer /0/  # position of particle birth z in array pid
                       # (FORTRAN indexed: based 1)
uxbirthpid integer /0/ # position of particle birth ux in array pid
                       # (FORTRAN indexed: based 1)
uybirthpid integer /0/ # position of particle birth uy in array pid
                       # (FORTRAN indexed: based 1)
uzbirthpid integer /0/ # position of particle birth uz in array pid
                       # (FORTRAN indexed: based 1)
xoldpid   integer /0/  # position of particles previous position in x in array pid
                       # (FORTRAN indexed: based 1)                  
yoldpid   integer /0/  # position of particles previous position in y in array pid
                       # (FORTRAN indexed: based 1)                  
zoldpid   integer /0/  # position of particles previous position in z in array pid
uxoldpid   integer /0/  # position of particles previous x velocity in array pid
                       # (FORTRAN indexed: based 1)                  
uyoldpid   integer /0/  # position of particles previous y velocity in array pid
                       # (FORTRAN indexed: based 1)                  
uzoldpid   integer /0/  # position of particles previous z velocity in array pid
                       # (FORTRAN indexed: based 1)                  
vdxoldpid   integer /0/  # position of particles previous x drift velocity in array pid
                       # (FORTRAN indexed: based 1)                  
vdyoldpid   integer /0/  # position of particles previous y drift velocity in array pid
                       # (FORTRAN indexed: based 1)                  
vdzoldpid   integer /0/  # position of particles previous z drift velocity in array pid
                       # (FORTRAN indexed: based 1)                  
dxpid   integer /0/    # position of particles dx grid cell size in array pid
                       # (FORTRAN indexed: based 1)                  
dypid   integer /0/    # position of particles dy grid cell size in array pid
                       # (FORTRAN indexed: based 1)                  
dzpid   integer /0/    # position of particles dz grid cell size in array pid
                       # (FORTRAN indexed: based 1)                  
bxpid   integer /0/  # position of particles Bx field in array pid
                       # (FORTRAN indexed: based 1)                  
bypid   integer /0/  # position of particles By field in array pid
                       # (FORTRAN indexed: based 1)                  
bzpid   integer /0/  # position of particles Bz field in array pid
                       # (FORTRAN indexed: based 1)                  
bxoldpid   integer /0/  # position of particles old Bx field in array pid
                       # (FORTRAN indexed: based 1)                  
byoldpid   integer /0/  # position of particles old By field in array pid
                       # (FORTRAN indexed: based 1)                  
bzoldpid   integer /0/  # position of particles old Bz field in array pid
                       # (FORTRAN indexed: based 1)                  
chdtspid integer /0/   # position of particle pid for dts change
                       # (FORTRAN indexed: based 1)
uparBopid  integer /0/ # position of particles uparallel * B in array pid

ssn    integer   /1/   # next particles 'social security number' available

lsaveoldpos logical /.false./ # Flag setting whether old particle positions are saved
particlesortbyindex(pgroup:ParticleGroup,pindex(np):integer,pindexmin:integer,
                    ipmin:integer,np:integer,
                    nn:integer,npblock(nn):integer) subroutine
                       # Sorts particles based on an inputted index

%%%%%%%%%%% ParticleGroup:
# Dynamic particle arrays, and related data
ns     integer    /1/  # Number of species
npmax  integer    /0/  # Size of data arrays
npid   integer    /0/  # number of columns for pid.
ipmax_s(0:ns)  _integer [1] /0/
   # Maximum index of particles of each species
   # Index of 0 is used as a guard and always has a value of zero.
sm(ns) _real [kg] /0./ # Species mass
sq(ns) _real [C]  /0./ # Species charge
sw(ns) _real [1]  /0./ # Species weight
                       # (real particles per simulation particles)
ins(ns)  _integer /1/  # Index of first particle in species
nps(ns)  _integer /0/  # Number of particles in species
sid(0:ns-1) _integer /-1/ # Global species index for each species
ndts(0:ns-1) _integer /1/  # Stride for time step advance for each species
ldts(0:ns-1) _logical /1/
lvdts(0:ns-1)  _logical /1/
iselfb(0:ns-1) _integer /0/ # Group number for particles that are affected by
                            # their own magnetic field, using the 1/gamma**2
                            # approximation. The correction is not applied to
                            # group number -1.
fselfb(0:ns-1) _real   /0./ # The scaling factor, vz.
dtscale(ns) _real /1./ # Scale factor applied to time step size for each
                       # species. Only makes sense in steaday and and
                       # transverse slice modes.
limplicit(0:ns-1) _logical /0/ # Flags implicit particle species
iimplicit(0:ns-1) _integer /-1/ # Group number for implicit particles
zshift(ns) _real /0./
lebcancel        logical   /.false./ # turns on/off cancellation of E+VxB before V push
gaminv(npmax)   _real [1]  /1./ # inverse relativistic gamma factor
xp(npmax)       _real [m]       # X-positions of particles
yp(npmax)       _real [m]       # Y-positions of particles
zp(npmax)       _real [m]       # Z-positions of particles
uxp(npmax)      _real [m/s]     # gamma * X-velocities of particles
uyp(npmax)      _real [m/s]     # gamma * Y-velocities of particles
uzp(npmax)      _real [m/s]     # gamma * Z-velocities of particles
ex(npmax)       _real [m]       # Ex of particles
ey(npmax)       _real [m]       # Ey of particles
ez(npmax)       _real [m]       # Ez of particles
bx(npmax)       _real [m]       # Bx of particles
by(npmax)       _real [m]       # By of particles
bz(npmax)       _real [m]       # Bz of particles
pid(npmax,npid) _real [1]       # Particle ID - used for various purposes

*********** Subcycling dump:
ndtsaveraging integer /0/ # Sets the type of averaging to do when using
                          # subcycling.
                          # When 0, no averaging is done, the rho from
                          # each ndts group closest in time is used.
                          # When 1, use linear weighting of the old and new
                          # rho of the slow particles.
                          # When 2, average rho of the fast particles over
                          # n-1/2 to n+1/2, and the two halves of the velocity
                          # advance use the same self field. This does not work
                          # with even ndts step factors. Not complete.
                          # When 3, the fields for the two velocity halves are
                          # different, the first averaging from n-1/2 to n,
                          # the second from n to n+1/2. Not complete.
ndtsmax integer /1/       # Maximum step size factor
nsndts integer /0/        # Number of different step size factors
nsndtsphi integer /0/     # The number of extra copies of the phi array.
                          # Either 1 or same as nsndts depending on the type
                          # of subcycling.
ndtstorho(ndtsmax) _integer /-1/ # Converts the step size factor to the
                                 # rho array index to use.
ndts(0:nsndts-1) _integer /1/    # Timestep stride for subcycling groups
ldts(0:nsndts-1) _logical /1/    # Set to true when the particle positions
                                 # have been advanced (and the rho updated)
lvdts(0:nsndts-1) _logical /1/   # Set to true when the particle velocities
                                 # have been advanced
itndts(0:nsndts-1) _integer /0/  # The time level of the position of the
                                 # ndts group
zgridndts(0:nsndts-1) _real      # z location of the grid for each ndts group
nrhopndtscopies integer /1/ # Number of copies of rho for each ndts group
                           # It defaults to 1 which is what is needed if
                           # there are only groups with ndts==1. Otherwise
                           # it will be 2.
setupSubcycling(pgroup:ParticleGroup):integer) subroutine
                        # Sets up data for particle subcycling
getnsndtsforsubcycling()
             integer function # Get number of ndts groups to loop over


*********** SelfB:
nsselfb integer /0/     # Number of groups with different gammas that
                        # require the 1/gamma**2 correction.
fselfb(0:nsselfb-1) _real /0./ # The scaling factor, vz.
iselfb(0:ns-1) _integer /0/ # Group number for particles that are affected by
                            # their own magnetic field, using the 1/gamma**2
                            # approximation. Note that by default, group 0
                            # has gamma == 1.
setupSelfB(pgroup:ParticleGroup) subroutine
                        # Sets up data for particle requiring self B correction

*********** ImplicitModule:
nsimplicit integer /0/  # Number of implicit groups with different q/m values
implicitfactor(0:ns-1) _real # Coefficients for chi for each implicit group
setupImplicit(pgroup:ParticleGroup) subroutine
                        # Sets up the bookkeepping for implicit particle species

*********** Scraped_Particles dump:
# Arrays for scraped particles
scr_np             integer  /0/   # Total no. of scraped particles. 
scr_npmax          integer  /0/   # Max. no. of scraped particles.
scr_npbunch        integer  /10000/ # Number of particles in a bunch
scr_ns             integer  /0/   # number of species
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
npidlost            integer /1/ # Number of columns in pidlist
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
pidlost(npmaxlost,npidlost) _real [1] # Particle ID of lost particles

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
lresetparticlee           logical /.true./
   # When true, the particle's E field is reset to zero at the beginning
   # of each step.

*********** ExtPart dump:
# Arrays that hold particle data that is extrapolated to grid cell centers
# in the GETZMMNT and GETEXTRAPOLATEDPARTICLES routines.
nepwin                  integer /0/ # Number of grid locations (windows)
izepwin(nepwin)        _integer /-1/ # List of grid locations (indx of zmntmesh)
zzepwin(nepwin)        _real        # List of lab frame locations
wzepwin(nepwin)        _real        # List of lab frame widths
nepmax                  integer /0/ # Maximum number of particles
npidepmax               integer /0/ # Max number of columns in pidep
lepsaveonce             logical /.false./ # When true, each particle is saved
                                          # at most only once in each window
                                          # This only works if the windows
                                          # do not overlap each other!
epclosenessfactor       real /0.75/ # When saving particles only once,
                                    # particle is saved if the extrapolation
                                    # factor is this much less than the
                                    # estimated extrapolation factor for the
                                    # next step.
epflagpid               integer /0/ # Position in the array pid where the
                                    # saveonce flag is stored.
                                    # (FORTRAN indexed: based 1)
nep(nepwin,ns)         _integer     # Number of particles in each grid location
tep(nepmax,nepwin,ns)  _real        # time of particles at grid cell centers
xep(nepmax,nepwin,ns)  _real        # X coordinates at grid cell centers
yep(nepmax,nepwin,ns)  _real        # Y coordinates at grid cell centers
uxep(nepmax,nepwin,ns) _real        # X velocities at grid cell centers
uyep(nepmax,nepwin,ns) _real        # Y velocities at grid cell centers
uzep(nepmax,nepwin,ns) _real        # Z velocities at grid cell centers
pidep(nepmax,npidepmax,nepwin,ns) _real # Particle ID 

*********** TopPhys:
# "Physics" subroutines at top level
derivqty()  subroutine # Calculates global derived qtys.
getzmmnt(np,xp(np):real,yp(np):real,zp(np):real,
         uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,
         q:real,m:real,w:real,dt:real,dtscale:real,itask,nplive,
         uxpo(np):real,uypo(np):real,uzpo(np):real,
         is:integer,isid:integer,ismax:integer,
         maxp:real,minp:real,zmmnts0:real,zmmnts:real)
            subroutine # Sets moments as a function of z for species 1
getzmmnt_weights(np,xp(np):real,yp(np):real,zp(np):real,
         uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,
         wp(np):real,q:real,m:real,w:real,dt:real,dtscale:real,itask,nplive,
         uxpo:real,uypo:real,uzpo:real,is:integer,isid:integer,ismax:integer,
         maxp:real,minp:real,zmmnts0:real,zmmnts:real)
            subroutine # Sets moments as a function of z for species 1 with variables weights
getlabwn()  subroutine # Gets the lab window moments from the z moments
getvzofz()  subroutine
setgamma(lrelativ)
            subroutine # Converts v to u, sets gammainv for all ptcls
gammaadv(np,gaminv(np):real,uxp(np):real,uyp(np):real,uzp(np):real,
         gamadv:string,lrelativ)
            subroutine # Advances gamma
resetlat()  subroutine # Resizes lattice arrays to their true lengths
resetlatdrft() subroutine # Resizes drft lattice arrays to their true lengths
resetlatbend() subroutine # Resizes bend lattice arrays to their true lengths
resetlatdipo() subroutine # Resizes dipo lattice arrays to their true lengths
resetlatquad() subroutine # Resizes quad lattice arrays to their true lengths
resetlatsext() subroutine # Resizes sext lattice arrays to their true lengths
resetlathele() subroutine # Resizes hele lattice arrays to their true lengths
resetlataccl() subroutine # Resizes accl lattice arrays to their true lengths
resetlatemlt() subroutine # Resizes emlt lattice arrays to their true lengths
resetlatmmlt() subroutine # Resizes mmlt lattice arrays to their true lengths
resetlategrd() subroutine # Resizes egrd lattice arrays to their true lengths
resetlatbgrd() subroutine # Resizes bgrd lattice arrays to their true lengths
resetlatpgrd() subroutine # Resizes pgrd lattice arrays to their true lengths
resetlatbsqgrad() subroutine # Resizes bsqgrad lattice arrays to their true lengths
setlatt()   subroutine # Sets lattice pointers for the current beam location
setlattzt(zbeam:real,time:real,fstype:integer)
            subroutine # Sets lattice pointers at zbeam and time
species()   subroutine # Sets species related arrays.
zgapcorr(np:integer,zp(np):real,xp(np):real,uzp(np):real,gaminv(np):real,
         dtl:real,dtr:real,dt:real,m:real,q:real,time:real)
            subroutine # Adds z correction term for accelerating gap residence
                       # correction.
acclbfrm(zcorrection)
            subroutine # Gets the acceleration residence correction for
                       # the beam frame

*********** TopDiag:
# Subroutines in package TOP
setgrid1d(np:integer,x(np):real,nx:integer,grid(0:nx):real,xmin:real,xmax:real)
        subroutine
        # Deposits data onto a 1-D grid.
deposgrid1d(itask:integer,np:integer,x(np):real,z(np):real,nx:integer,
            grid(0:nx):real,gridcount(0:nx):real,xmin:real,xmax:real)
        subroutine
        # Deposits data onto a 1-D grid.
getgrid1d(np:integer,x(np):real,z(np):real,nx:integer,grid(0:nx):real,
          xmin:real,xmax:real)
        subroutine
        # Gathers data from a 1-D grid.
setgrid2d(np:integer,x(np):real,y(np):real,nx:integer,ny:integer,
          grid(0:nx,0:ny):real,
          xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Deposits uniform data onto a 2-D grid.
deposgrid2d(itask:integer,np:integer,x(np):real,y(np):real,z(np):real,
            nx:integer,ny:integer,
            grid(0:nx,0:ny):real,gridcount(0:nx,0:ny):real,
            xmin:real,xmax:real,ymin:real,ymax:real)
        subroutine
        # Deposits data onto a 2-D grid.
setgrid2dw(np:integer,x(np):real,y(np):real,w(np):real,
           nx:integer,ny:integer,grid(0:nx,0:ny):real,
           xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Deposits uniform data onto a 2-D grid.
deposgrid2dw(itask:integer,np:integer,x(np):real,y(np):real,z(np):real,
             w(np):real,nx:integer,ny:integer,
            grid(0:nx,0:ny):real,gridcount(0:nx,0:ny):real,
            xmin:real,xmax:real,ymin:real,ymax:real)
        subroutine
        # Deposits data onto a 2-D grid.
deposgridrzvect(itask:integer,np:integer,x(np):real,y(np):real,z(np):real,
            vx(np):real,vy(np):real,vz(np):real,
            w(np):real,nx:integer,ny:integer,
            grid(0:nx,0:ny,3):real,gridcount(0:nx,0:ny):real,
            xmin:real,xmax:real,ymin:real,ymax:real)
        subroutine
        # Deposits velocities onto a 2-D radial grid.
getgrid2d(np:integer,x(np):real,y(np):real,z(np):real,
          nx:integer,ny:integer,grid(0:nx,0:ny):real,
          xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Gathers data from a 2-D grid.
getgridngp2d(np:integer,x(np):real,y(np):real,z(np):real,
             nx:integer,ny:integer,grid(0:nx,0:ny):real,
             xmin:real,xmax:real,ymin:real,ymax:real) subroutine
        # Gathers data from a 2-D grid using nearest grid point
setgrid3d(np:integer,x(np):real,y(np):real,z(np):real,
          nx:integer,ny:integer,nz:integer,
          grid(0:nx,0:ny,0:nz):real,
          xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real)
          subroutine
        # Deposits uniform data onto a 3-D grid.
deposgrid3d(itask:integer,np:integer,x(np):real,y(np):real,z(np):real,
            q(np):real,
            nx:integer,ny:integer,nz:integer,
            grid(0:nx,0:ny,0:nz):real,gridcount(0:nx,0:ny,0:nz):real,
            xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real)
        subroutine
        # Deposits data onto a 3-D grid.
deposgrid3dvect(itask:integer,np:integer,x(np):real,y(np):real,z(np):real,
            vx(np):real,vy(np):real,vz(np):real,w(np):real,
            nx:integer,ny:integer,nz:integer,
            grid(0:nx,0:ny,0:nz,3):real,gridcount(0:nx,0:ny,0:nz):real,
            xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real)
        subroutine
        # Deposits velocities onto a 3-D grid.
getgrid3d(np:integer,x(np):real,y(np):real,z(np):real,f(np):real,
          nx:integer,ny:integer,nz:integer,grid(0:nx,0:ny,0:nz):real,
          xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real,
          l2symtry:logical,l4symtry:logical) subroutine
        # Gathers data from a 3-D grid.
getgridngp3d(np:integer,x(np):real,y(np):real,z(np):real,f(np):real,
             nx:integer,ny:integer,nz:integer,grid(0:nx,0:ny,0:nz):real,
             xmin:real,xmax:real,ymin:real,ymax:real,zmin:real,zmax:real,
             zgrid:real,l2symtry:logical,l4symtry:logical) subroutine
        # Gathers data from a 3-D grid using nearest grid point.



grid2grid(unew(0:nxnew,0:nynew):real,nxnew:integer,nynew:integer,
          xminnew:real,xmaxnew:real,yminnew:real,ymaxnew:real,
          uold(0:nxold,0:nyold):real,nxold:integer,nyold:integer,
          xminold:real,xmaxold:real,yminold:real,ymaxold:real) subroutine
        # project field from one grid to another
gridtogrid3d(nxin:integer,nyin:integer,nzin:integer,
             xminin:real,xmaxin:real,yminin:real,ymaxin:real,
             zminin:real,zmaxin:real,
             gridin(0:nxin,0:nyin,0:nzin):real,
             nxout:integer,nyout:integer,nzout:integer,
             xminout:real,xmaxout:real,yminout:real,ymaxout:real,
             zminout:real,zmaxout:real,
             gridout(0:nxout,0:nyout,0:nzout):real) subroutine
        # Linearly interpolates from one grid to another. This will also work
        # for 2d and 1d arrays if the n's are set to zero.
sum_neighbors3d(fin(nx+1,ny+1,nz+1):integer,fout(nx+1,ny+1,nz+1):integer,nx:integer,ny:integer,nz:integer) subroutine # sum neighbouring cells
take2dint(a(0:n1-1,0:n2-1):integer,n1:integer,n2:integer,
          i(n):integer,j(n):integer,n:integer,b(n):integer) subroutine
getpsgrd(np,xp(np):real,uxp(np):real,nw,nh,psgrd(0:nw,0:nh):real,
         wmin:real,wmax:real,hmin:real,
         hmax:real,zl:real,zr:real,zp(np):real,uzp(np):real,slope:real)
              subroutine # lays particles onto slanted mesh in phase space
emitthresh(pgroup:ParticleGroup,n:integer,threshold(n):real,js:integer,
           iw:integer,
           ngridw:integer,ngridh:integer,tepsx(n):real,tepsy(n):real)
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
emitfrac(xp(np):real,uxp(np):real,uzp(np):real,np:integer,xbar:real,xpbar:real,
         xsqbar:real,xxpbar:real,xpsqbar:real,
         fracbin(0:npts):real,emitbin(0:npts):real,npts:integer,emitbinmax:real, 
         tx(np):real,txp(np):real,emitp(np):real,
         rwork(0:npts):real,iwork:integer)
         subroutine # Calculates the emittance versus the fraction of the
                    # beam (live particles) particles in phase space 
                    # enclosed by nested emittance ellipses with the rms 
                    # equivalent beam.  Also returns single particle 
                    # emittances.
prin(xp(np):real,uxp(np):real,uzp(np):real,np:integer,xsqbar:real,xpsqbar:real,
     xxpbar:real,xbar:real,xpbar:real,tx(np):real,txp(np):real,
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
thiszbeam(zl:real,zr:real,control(ncontrol):real,ncontrol:integer)
              logical function
                         #
savehist(time:real)
              subroutine # saves moments data to history arrays
prntpara(dxperp:real,dz:real)
              subroutine # Print various parameters to various files
psplots(freqflag:integer) 
              subroutine # Controls phase space plots
onedplts(freqflag:integer) subroutine # plots all 1d qtys w/ freqflag
tolabfrm(zcent:real,nn,x(nn):real,z(nn):real) subroutine
             # Converts data from WARP frame to lab frame.

*********** TopUtil:
# "Utility" subroutines at top level
getfortantruefalse(tval:logical,fval:logical) subroutine
                        # Returns the fortran values of true and false
nextpid() integer function
                        # Returns the next value of npid. Note that this
                        # function should be used rather that directly
                        # changing npid.
dolabwn() logical function
                        # Checks if lab window is in beam frame
setuppgroup(pgroup:ParticleGroup) subroutine
                        # Does basis setup of pgroup data, setting ns and npid
alotpart(pgroup:ParticleGroup)
             subroutine # Allocate space for particles and setup 
                        # associated data
chckpart(pgroup:ParticleGroup,is:integer,nlower:integer,nhigher:integer,
         lfullshft:logical)
             subroutine # Makes sure there is enough space for nn particles.
addpart(pgroup:ParticleGroup,nn:integer,npid:integer,
        x(nn):real,y(nn):real,z(nn):real,vx(nn):real,vy(nn):real,vz(nn):real,
        gi(nn):real,pid(nn,npid):real,
        is:integer,lallindomain:logical,zmmin:real,zmmax:real,
        lmomentum:logical)
             subroutine # Adds new particles to the simulation
clearpart(pgroup:ParticleGroup,js:integer,fillmethod:integer)
             subroutine # Clears away lost particles.
shrinkpart(pgroup:ParticleGroup)
             subroutine # Removes unused space in the particle arrays
particlesortyzwithcopy(pgroup:ParticleGroup,dy:real,dz:real,
                       ymmin:real,zmmin:real,ny:integer,nz:integer)
             subroutine # Sorts particles via a full copy
processlostpart(pgroup:ParticleGroup,is:integer,clearlostpart:integer,
                time:real,zbeam:real)
             subroutine # Processes lost particles (particles which have
                        # gaminv set to zero).
alotlostpart() subroutine # Allocate space for saving lost particles
chcklostpart(is:integer,nlower:integer,nhigher:integer) 
             subroutine # Make sure that there is enough space in the lost particle arrays for nlower
                        # new particles below and nhigher above the lost particles.  Returns if
                        # there is already enough space above and below.  If there is enough total
                        # space but not enough room above or below, the lost particles are shifted
                        # appropriately. If there is not enough space, add more to the arrays.
                        # Particle data is shifted appropriately.
checkparticlegroup(pgroup:ParticleGroup,is:integer,
                   nlower:integer,nhigher:integer)
             subroutine # Makes sure there is enough space for nn particles in
                        # the given particle group.
copygrouptogroup(pgroupin:ParticleGroup,nn:integer,ii:integer,istart:integer,
                pgroupout:ParticleGroup,it:integer)
             subroutine # Copies particle data between groups
load2d(np,x(np):real,y(np):real,nx,ny,n(0:nx,0:ny):real,dx:real,dy:real)
             subroutine # Loads particles approximately into a 2-D distribution
shftpart(pgroup:ParticleGroup,is:integer,ishft:integer) subroutine
                        # Moves particle data to end of species group.
copypart(pgroup:ParticleGroup,it:integer,nn:integer,ii:integer;istart:integer)
             subroutine
                        # Copies particle data from one location to another
wrandom(method:string,ii:integer,idig:integer,ifib1:integer,ifib2:integer)
          real function # Returns random number using the given method
wrandomgauss(method:string,ii:integer,idig1:integer,idig2:integer,
             ifib1:integer,ifib2:integer,lefficient:logical)
          real function # Returns Gaussian random number using the given method
r2rev(lastrev:real)
          real function # bit reversed counter
                        # lastrev is `call by address', e.g. r2rev(&x)
rnrev(i:integer,nbase:integer)
          real function # i base `nbase' reversed.
rnrevarray(n:integer,x(n):real,i:integer,nbase:integer)
          subroutine    # Fills an array with uniform digit reversed rand numbs
rnorm()   real function # Gaussian random numbers.
rnormdig(i1,n,nbase1,nbase2,dx:real,x(n):real)
             subroutine # Gaussian random numbers via digit reversed.
rm()      real function # Pseudo-Gaussian random numbers (6 uniform nos.).
rma(a(n):real,n) subroutine # rma(&a,n) gives n Pseudo-Gaussian rand numbers.
dsifa(a:real,lda,n,kpvt,info)
             subroutine # LINPACK matrix reduction routine # (in file UTIL.F)
dsidi(a:real,lda,n,kpvt,det:real,inert,work:real,job:real)
             subroutine # LINPACK matrix inverting routine # (in file UTIL.F)
dsisl(a:real,lda,n,kpvt,b:real)
             subroutine # LINPACK routine to solve a*x = b # (in file UTIL.F)
dvdfit(x:real,y:real,ndata:integer,ndatap:integer,a:real,
       basis:real,ma:integer,map:integer,u:real,w:real,v:real,tmp:real,
       tol:real,chisq:real)    subroutine
       # matrix linear least squares fit using an input singular value 
       # decomposition of a basis matrix                   (in util.m) 
dvdcmp(a:real,m:integer,n:integer,mp:integer,np:integer,
       w:real,v:real,tmp:real) subroutine
       # singular value decomposition of the matrix a      (in util.m)
dvbksb(u:real,w:real,v:real,m:integer,n:integer,mp:integer,np:integer,
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
getbeamcom(pgroup:ParticleGroup) real function
       # Returns the center of mass in z of the beam calculated from the
       # particles.
zpartbnd(pgroup:ParticleGroup,zmmax:real,zmmin:real,dz:real)
       subroutine
       # Enforces axial particle boundary conditions
zpartbndwithdata(n:integer,z(n):real,uz(n):real,gaminv(n):real,
                 zmmax:real,zmmin:real,dz:real,zgrid:real) subroutine
       # Enforces axial particle boundary conditions
reorgparticles(pgroup:ParticleGroup) subroutine
       # Reorganizes particles for the parallel version

******* Parallel:
nslaves       integer /0/         # Number of slaves
my_index      integer   +parallel # Processor index to array of task ids
grid_overlap  integer       +dump # Overlap of field grid in processors
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
izfsslave(0:maxslaves-1) _integer +dump # starting iz for which each slave does
                                        # a field solve calculation
nzfsslave(0:maxslaves-1) _integer +dump # number of z grid cells for which each
                                        # slave does a field solve calculation
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
n_STdiscs             integer /0/  # Number of semitransparent discs
z_STdiscs(n_STdiscs)  _real        # Positions of semitransparent discs
r_STdiscs(n_STdiscs)  _real        # Radii of semitransparent discs
t_STdiscs(n_STdiscs)  _real        # Transparency of semitransparent discs
semitransparent_disc(dz:real) subroutine
    # Randomly absorb particles passing through disc based on data in group
    # SemiTransparentDisc. dz is the maximum distance traveled by particles
    # in one time step.
lenvfix_STdisc(n_STdiscs)  _logical
    # When true, the perveance is scaled as the beam nears the disc to model
    # the shorting out of the self fields.
envfixzscale_STdisc(n_STdiscs) _real /LARGEPOS/ [m]
    # Length (in meters) over which the perveance correction is applied
    # to the envelope calculation.
envfixsscale_STdisc(n_STdiscs) _real /LARGEPOS/ [1]
    # Length (in beam radii) over which the perveance correction is applied
    # to the envelope calculation.

*********** Temperatures dump:
nstemp  integer /1/                                   # nb species
evolt   integer /1/                                   # temperature units in evolt
joule   integer /2/                                   # temperature units in joule
kelvin  integer /3/                                   # temperature units in kelvin
t_units integer /1/                                   # temperature units
l_temp_collapseinz     logical /.false./              # collapse z-slices in temperature calculations
l_temp_rmcorrelations  logical /.true./               # remove x*v correlations in temperature calculations
nxtslices    integer    /0/                           # nb cells in x of temperature slices
nytslices    integer    /0/                           # nb cells in y of temperature slices
nztslices    integer    /1/                           # nb temperature slices
nxtslicesc   integer    /0/                           # nb cells in x of temperature slices (correlation variables)
nytslicesc   integer    /0/                           # nb cells in y of temperature slices (correlation variables)
nztslicesc   integer    /1/                           # nb temperature slices               (correlation variables)
tslicexmin(nztslices) _real                           # min in x of temperature slices
tslicexmax(nztslices) _real                           # max in x of temperature slices
tsliceymin(nztslices) _real                           # min in y of temperature slices
tsliceymax(nztslices) _real                           # max in y of temperature slices
tslicezmin(nztslices) _real                           # min in z of temperature slices
tslicezmax(nztslices) _real                           # max in z of temperature slices
dxti(nztslices)       _real                           # inverse mesh size in x for temperature slices
dyti(nztslices)       _real                           # inverse mesh size in y for temperature slices
pnumt(0:nxtslices,0:nytslices,nztslices)        _real # nb particles    in temperature slices
pnumtw(0:nxtslices,0:nytslices,nztslices)       _real # weights         in temperature slices
vxbart(0:nxtslices,0:nytslices,nztslices)       _real # average of vx   in temperature slices
vybart(0:nxtslices,0:nytslices,nztslices)       _real # average of vy   in temperature slices
vzbart(0:nxtslices,0:nytslices,nztslices)       _real # average of vz   in temperature slices
vxsqbart(0:nxtslices,0:nytslices,nztslices)     _real # average of vx^2 in temperature slices
vysqbart(0:nxtslices,0:nytslices,nztslices)     _real # average of vy^2 in temperature slices
vzsqbart(0:nxtslices,0:nytslices,nztslices)     _real # average of vz^2 in temperature slices
xbart(0:nxtslicesc,0:nytslicesc,nztslicesc)     _real # average of x    in temperature slices
ybart(0:nxtslicesc,0:nytslicesc,nztslicesc)     _real # average of y    in temperature slices
zbart(0:nxtslicesc,0:nytslicesc,nztslicesc)     _real # average of z    in temperature slices
xsqbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of x^2  in temperature slices
ysqbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of y^2  in temperature slices
zsqbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of z^2  in temperature slices
xvxbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of x*vx in temperature slices
yvybart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of y*vy in temperature slices
zvzbart(0:nxtslicesc,0:nytslicesc,nztslicesc)   _real # average of z*vz in temperature slices
tempx(0:nxtslices,0:nytslices,nztslices,nstemp) _real # x-temperature  in temperature slices
tempy(0:nxtslices,0:nytslices,nztslices,nstemp) _real # y-temperature  in temperature slices
tempz(0:nxtslices,0:nytslices,nztslices,nstemp) _real # z-temperature  in temperature slices
tempxz(nztslices,nstemp) _real                        # x-temperature  in temperature slices
tempyz(nztslices,nstemp) _real                        # y-temperature  in temperature slices
tempzz(nztslices,nstemp) _real                        # z-temperature  in temperature slices
nztlocator integer /1/                                # nb meshes locator temperature slices
ntlmax     integer /1/                                # max nb of slices in locator temperature slices
ntl(nztlocator) _integer                              # nb of slices in locator temperature slices
tslice_locator(nztlocator,ntlmax) _integer            # locator array for temperature array
tloc_zmin  real                                       # min locator temperature array 
tloc_zmax  real                                       # max locator temperature array 
tloc_dzi   real                                       # inverse mesh size locator temperature array 
gett(is:integer,lrtheta:logical,
     l2symtry:logical,l4symtry:logical) subroutine    # get temperature for species is (r and theta temp. if lrtheta=true) 
setregulartgrid(nx:integer,ny:integer,nz:integer,
                xmin:real,xmax:real,ymin:real,
                ymax:real,zmin:real,zmax:real,
                dz:real,nzloc:integer,
                lcollapse:logical,
                lcorrel:logical) subroutine           # Set slices regularly in z in a box delimited by (xmin,xmax,ymin,ymax,zmin,zmax).
                                                      # nzloc refers to the number of nodes of a lookup table for fast 
                                                      # localization of temperature slices. In general, set nzloc=w3d.nz.
impact_ion(is1:integer,is2:integer,nbp:real,w:real,
           shiftx:real,shifty:real,shiftz:real,
           deltax:real,deltay:real,deltaz:real,condid:integer)
           subroutine # create particles of species is2 from impact of particles if species is1 
                      # nbp particles are created for each incident particle
                      # w: energy emitted particles (all in -Vz)
                      # We have
                      # xpemit = xpabsorb + shiftx + (wranf()-0.5)*deltax 
                      # ypemit = ypabsorb + shifty + (wranf()-0.5)*deltay 
                      # zpemit = zpabsorb + shiftz + (wranf()-0.5)*deltaz 
chgparticlesdts(pgroup:ParticleGroup) subroutine 

******** Subtimerstop:
ltoptimesubs logical /.false./
timezbendcor                   real /0./
timecalculatebsqgrad           real /0./
timealotpart                   real /0./
timechckpart                   real /0./
timeshftpartwork               real /0./
timecopypart                   real /0./
timeaddpart                    real /0./
timeclearpart                  real /0./
timeshrinkpart                 real /0./
timeprocesslostpart            real /0./
timealotlostpart               real /0./
timechcklostpart               real /0./
timeshftlostpart               real /0./
timecheckparticlegroup         real /0./
timeshiftparticlegroup         real /0./
timecopyparttogroup            real /0./
timecopygrouptopart            real /0./
timezpartbnd_slave             real /0./
timereorgparticles_parallel    real /0./
timecheckzpartbnd              real /0./
timeparallel_sum_mmnts         real /0./
timeparallel_sum_temperature   real /0./
timeparallelsumrealarray       real /0./
timeparallelsumintegerarray    real /0./
timeparallelmaxrealarray       real /0./
timeparallelminrealarray       real /0./
timeparallellor                real /0./
timeparallelbroadcastrealarray real /0./
timeparallelbarrier            real /0./
timeparallelnonzerorealarray   real /0./

