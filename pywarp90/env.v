env
#@(#) File ENV.V, version $Revision: 3.4 $, $Date: 2001/05/15 01:51:02 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package ENV of the WARP code.
# Envelope data generated from this package can be used in WARP3d.   
# ENV - envelope code
# Alex Friedman, LLNL, (510) 422-0827  friedman@hif.llnl.gov 
# Dave Grote,    LLNL, (510) 423-7194  dave@hif.llnl.gov 
# Steve Lund,    LLNL, (510) 423-4463  lund@hif.llnl.gov 
{
}

*********** ENVversion:
# Version control for env package
versenv character*19 /"$Revision: 3.4 $"/  # Current code version, set by CVS

*********** ENVvars dump:
# Variables needed by the package ENV
dzenv   real /0.001/ [m] # Step size in envelope calculation
zl      real /-1.2/  [m] # Starting z for envelope calc
zu      real /3./    [m] # Maximum  z for envelope calc
tunezs  real /0./    [m] # Start of region over which tune is calculated;
                         # defaults to zl
tuneze  real /0./    [m] # End of region over which tune is calculated;
                         # defaults to zl+tunelen
sig0x   real [deg] # Undepressed particle x-phase advance over tunezs -> tuneze
sigx    real [deg] # Depressed   particle x-phase advance over tunezs -> tuneze 
sig0y   real [deg] # Undepressed particle y-phase advance over tunezs -> tuneze
sigy    real [deg] # Depressed   particle y-phase advance over tunezs -> tuneze 
deltaa  real [m]   # Change in a  over tunezs -> tuneze 
deltaap real [1]   # Change in ap over tunezs -> tuneze 
deltab  real [m]   # Change in b  over tunezs -> tuneze 
deltabp real [1]   # Change in bp over tunezs -> tuneze 
genprv  real [?]   # Beam generalized perveance
nenv         integer  [1]     # Number of envelope points - 1
aenv(0:nenv)   _real  [m]     # Computed beam width in x
apenv(0:nenv)  _real  [1]     # Computed slope in x
benv(0:nenv)   _real  [m]     # Computed beam width in y
bpenv(0:nenv)  _real  [1]     # Computed slope in y
vzenv(0:nenv)  _real  [m/s]   # Computed axial velocity
fqxenv(0:nenv) _real  [1/m^2] # Quad    X-focusing force 
fqyenv(0:nenv) _real  [1/m^2] # Quad    Y-focusing force 
fuxenv(0:nenv) _real  [1/m^2] # Uniform X-focusing force
fuyenv(0:nenv) _real  [1/m^2] # Uniform Y-focusing force 
xenv(0:nenv)   _real  [m]     # X centroid of envelope
yenv(0:nenv)   _real  [m]     # Y centroid of envelope
xpenv(0:nenv)  _real  [m]     # X' centroid of envelope
ypenv(0:nenv)  _real  [m]     # Y' centroid of envelope
xorb(0:nenv)   _real  [m]     # X of main particle orbit
xporb(0:nenv)  _real  [m]     # X' of main particle orbit
yorb(0:nenv)   _real  [m]     # Y of main particle orbit
yporb(0:nenv)  _real  [m]     # Y' of main particle orbit
zenv(0:nenv)   _real  [m]     # Z coordinate of envelope point
lenvout logical /.true./ # Sets whether data is printed out after an env step
envtime        real   [s]     # CPU time for envelope calculation
iesemltq integer # Obsolete
iesemltu integer # Obsolete
imsmmltq integer # Obsolete

*********** ENVfofz dump:
# Variables which can vary as a function of z.  Data stored on a separate
# mesh, independent of mesh that envelope data is calculated on.
lefofz logical      /.false./   # Sets whether any variables are specified as a
                                # function of z.
dzefofz real              [m]   # Mesh cell size, defaults to dzenv
nzefofz integer     /0/         # Number of mesh points, defaults to nenv
zlefofz real        /0./  [m]   # Lower bound of mesh, defaults to zl.
zuefofz real        /0./  [m]   # Upper bound of mesh, defaults to zu.
libeame_z logical   /.false./   # Sets whether current varies as a function of z
ibeame_z(0:nzefofz) _real [A]   # Current as a function of z.
lemitne_z logical    /.false./  # Sets whether normalized emittance varies
                                # as function of z
emitnxe_z(0:nzefofz)  _real [A] # Normalized X-Emittance as a function of z
emitnye_z(0:nzefofz)  _real [A] # Normalized Y-Emittance as a function of z 

*********** ENVtune local:
# Position and angles of test particles at the start and end of the region
# over which the x- and y- tunes are calculated and associated variables 
# used in the tune calculations.   Block is included for possible debugging 
# use.  
astart    real [m] # Size  of X envelope at tunezs
apstart   real [m] # Angle of X envelope at tunezs
aend      real [m] # Size  of X envelope at tuneze
apend     real [m] # Angle of X envelope at tuneze
bstart    real [m] # Size  of Y envelope at tunezs
bpstart   real [m] # Angle of Y envelope at tunezs
bend      real [m] # Size  of Y envelope at tuneze
bpend     real [m] # Angle of Y envelope at tuneze
x1start   real [m] # Position of x-orbit 1 (with space charge) at tunezs
xp1start  real [m] # Angle    of x-orbit 1 (with space charge) at tunezs
x1end     real [m] # Position of x-orbit 1 (with space charge) at tuneze
xp1end    real [m] # Angle    of x-orbit 1 (with space charge) at tuneze
x2start   real [m] # Position of x-orbit 2 (with space charge) at tunezs
xp2start  real [m] # Angle    of x-orbit 2 (with space charge) at tunezs
x2end     real [m] # Position of x-orbit 2 (with space charge) at tuneze
xp2end    real [m] # Angle    of x-orbit 2 (with space charge) at tuneze
x01start  real [m] # Position of x-orbit 3 (without space charge) at tunezs
x0p1start real [m] # Angle    of x-orbit 3 (without space charge) at tunezs
x01end    real [m] # Position of x-orbit 3 (without space charge) at tuneze
x0p1end   real [m] # Angle    of x-orbit 3 (without space charge) at tuneze
x02start  real [m] # Position of x-orbit 4 (without space charge) at tunezs
x0p2start real [m] # Angle    of x-orbit 4 (without space charge) at tunezs
x02end    real [m] # Position of x-orbit 4 (without space charge) at tuneze
x0p2end   real [m] # Angle    of x-orbit 4 (without space charge) at tuneze
y1start   real [m] # Position of y-orbit 1 (with space charge) at tunezs
yp1start  real [m] # Angle    of y-orbit 1 (with space charge) at tunezs
y1end     real [m] # Position of y-orbit 1 (with space charge) at tuneze
yp1end    real [m] # Angle    of y-orbit 1 (with space charge) at tuneze
y2start   real [m] # Position of y-orbit 2 (with space charge) at tunezs
yp2start  real [m] # Angle    of y-orbit 2 (with space charge) at tunezs
y2end     real [m] # Position of y-orbit 2 (with space charge) at tuneze
yp2end    real [m] # Angle    of y-orbit 2 (with space charge) at tuneze
y01start  real [m] # Position of y-orbit 3 (without space charge) at tunezs
y0p1start real [m] # Angle    of y-orbit 3 (without space charge) at tunezs
y01end    real [m] # Position of y-orbit 3 (without space charge) at tuneze
y0p1end   real [m] # Angle    of y-orbit 3 (without space charge) at tuneze
y02start  real [m] # Position of y-orbit 4 (without space charge) at tunezs
y0p2start real [m] # Angle    of y-orbit 4 (without space charge) at tunezs
y02end    real [m] # Position of y-orbit 4 (without space charge) at tuneze
y0p2end   real [m] # Angle    of y-orbit 4 (without space charge) at tuneze

*********** Match:
mtch_tol    real    /1.e-4/   # Minimal allowed change in values in mtch
sig_desr    real    /20./     # Desired value of depressed tune shift, sigma
mtchdlta(6) real    /6*1.e-2/ # Initial deltas (fractional change) in the
                              # order a0, ap0, b0, bp0, ibeam, emit
wdeltaa     real    /1./      # Weight of deltaa in calculating error
wdeltab     real    /1./      # Weight of deltab in calculating error
wdeltaap    real    /1./      # Weight of deltaap in calculating error
wdeltabp    real    /1./      # Weight of deltabp in calculating error
wsig        real    /1.e-2/   # Weight of (sigma-sig_desr) in calculating error
mtch(7,6)   real              # Seven values for which the envelope is solved
mtcherrs(7) real              # Error in each of the seven cases
mtchiter    integer /100/     # Maximum number of allowed iterations
mtchinit()  subroutine        # Routine to initialize mtch and mtcherrs
                              # Must be called before call to envmatch
envmatch()  subroutine        # Routine to do the search

*********** ENVsubs:
#  Callable subroutines in the ENV package
envgen() subroutine
envexe() subroutine
envx()   subroutine
         # BASIS-level interface to ENVELOPE, using ENV database variables.
envxport(np,z:real,a:real,ap:real,b:real,bp:real,
         x:real,xp:real,y:real,yp:real) integer function
         # Export routine with envelope data at z(1:np)


