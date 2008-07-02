env
#@(#) File ENV.V, version $Revision: 3.35 $, $Date: 2008/07/02 19:50:41 $
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
versenv character*19 /"$Revision: 3.35 $"/  # Current code version, set by CVS

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
nenv            integer  [1]     # Number of envelope points - 1
aenv(0:nenv)      _real  [m]     # Computed beam width in x
apenv(0:nenv)     _real  [1]     # Computed slope in x
benv(0:nenv)      _real  [m]     # Computed beam width in y
bpenv(0:nenv)     _real  [1]     # Computed slope in y
vzenv(0:nenv)     _real  [m/s]   # Computed axial velocity
fqxenv(0:nenv)    _real  [1/m^2] # Quad    X-focusing force 
fqyenv(0:nenv)    _real  [1/m^2] # Quad    Y-focusing force 
fuxenv(0:nenv)    _real  [1/m^2] # Uniform X-focusing force
fuyenv(0:nenv)    _real  [1/m^2] # Uniform Y-focusing force 
xenv(0:nenv)      _real  [m]     # X centroid of envelope
yenv(0:nenv)      _real  [m]     # Y centroid of envelope
xpenv(0:nenv)     _real  [1]     # X' centroid of envelope
ypenv(0:nenv)     _real  [1]     # Y' centroid of envelope
xorb(0:nenv)      _real  [m]     # X of main particle orbit
xporb(0:nenv)     _real  [1]     # X' of main particle orbit
yorb(0:nenv)      _real  [m]     # Y of main particle orbit
yporb(0:nenv)     _real  [1]     # Y' of main particle orbit
xxenv(0:nenv)     _real  [m^2]   # Average value of XX   w.r.t. beam centroid
xxpenv(0:nenv)    _real  [m]     # Average value of XX'  w.r.t. beam centroid
xpxpenv(0:nenv)   _real  [1]     # Average value of X'X' w.r.t. beam centroid
yyenv(0:nenv)     _real  [m^2]   # Average value of YY   w.r.t. beam centroid
yypenv(0:nenv)    _real  [m]     # Average value of YY'  w.r.t. beam centroid
ypypenv(0:nenv)   _real  [1]     # Average value of Y'Y' w.r.t. beam centroid
xyenv(0:nenv)     _real  [m^2]   # Average value of XY   w.r.t. beam centroid
xpyenv(0:nenv)    _real  [m]     # Average value of X'Y  w.r.t. beam centroid
xypenv(0:nenv)    _real  [m]     # Average value of XY'  w.r.t. beam centroid
xpypenv(0:nenv)   _real  [1]     # Average value of X'Y' w.r.t. beam centroid
bphenv(0:nenv)    _real  [1]     # Skew-coupled beam phase angle
rxenv(0:nenv)     _real  [m]     # x-plane beam edge radius in rotating frame
ryenv(0:nenv)     _real  [m]     # y-plane beam edge radius in rotating frame
emitxenv(0:nenv)  _real  []      # x-plane emittance
emityenv(0:nenv)  _real  []      # y-plane emittance
emitnxenv(0:nenv) _real  []      # normalized x-plane emittance
emitnyenv(0:nenv) _real  []      # normalized y-plane emittance
emitngenv(0:nenv) _real  []      # normalized emittance-like invariant
emitnhenv(0:nenv) _real  []      # normalized emittance-like invariant
zenv(0:nenv)      _real  [m]     # Z coordinate of envelope point
lenvout logical /.true./ # Sets whether data is printed out after an env step
envtime           real   [s]     # CPU time for envelope calculation
llarmorframe   logical  /.false./  # Determines whether envelope radii, angles, 
			         # and phase advances are in Larmor frame
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
         # Python-level interface to ENVELOPE, using ENV database variables.
envxport(np,z(np):real,a(np):real,ap(np):real,b(np):real,bp(np):real,
         x(np):real,xp(np):real,y(np):real,yp(np):real,vz(np):real,
         emitx(np):real,emity(np):real,ibeam(np):real) integer function
         # Export routine with envelope data at z(1:np)
RK4HillSolve(karray(2*numsteps+1):real,si:real,sf:real,xi:real,xpi:real,
	     numsteps:integer,xxparray(numsteps+1,2):real) subroutine
	 # Solves Hill's Equation on interval [si,sf] subject to initial 
	 # conditions xi and xpi via the RK4 method
kappax(z:real) real function  # x-plane lattice focusing function
kappay(z:real) real function  # y-plane lattice focusing function
kappaxvec(zarray(n):real,karray(n):real,n:integer) subroutine
	 # Calculates kappax at each point in zarray and stores those values 
	 # in karray
kappayvec(zarray(n):real,karray(n):real,n:integer) subroutine
	 # Calculates kappay at each point in zarray and stores those values 
	 # in karray
ludcmp(a(np,np):real,n:integer,np:integer,indx(n):integer,d:real) subroutine
	 # Performs LU decomposition on matrix a
lubksb(a(np,np):real,n:integer,np:integer,indx(n):integer,b(n):real) subroutine
	 # Performs LU backsubstitution on the LU decomposition of a
matinv(a(np,np):real,n:integer,np:integer,indx(n):integer,d:real,
	 y(np,np):real) subroutine
	 # Computes the inverse of matrix a using ludcmp and lubksb
matprod(a(i,j):real,b(j,k):real,c(i,k):real,i:integer,j:integer,k:integer) subroutine
	 # Multiplies matrices a and b to yield matrix c
matvecprod(a(i,j):real,b(j):real,c(i):real,i:integer,j:integer) subroutine
	 # Multiplies matrix a and vector b to yield vector c