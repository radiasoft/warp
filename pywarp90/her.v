her
#@(#) File HER.V, version $Revision: 3.3 $, $Date: 2001/05/30 00:05:36 $
# Copyright (c) 1999, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for HERMES.
# Michiel de Hoon,    LBNL, (510) 486-5157  mdehoon@lbl.gov
{
}

*********** HERversion:
# Version control for her package
versher character*19 /"$Revision: 3.3 $"/  # Current code version, set by CVS

*********** HERvars dump:
# Variables needed by the package HER
dther real /0.0/  [m]    # Time step size
zlher real /0.0/  [m]    # Starting z for HERMES calc
zuher real /0.0/  [m]    # Maximum  z for HERMES calc
nther integer /0/        # Total number of time steps
niz       integer [1]    # Number of point along beam
aher(niz)   _real [m]    # Width in x (1)[0]
apher(niz)  _real [1]    # Slope in x (2)[1]
bher(niz)   _real [m]    # Width in y (3)[2]
bpher(niz)  _real [1]    # Slope in y (4)[3]
xher(niz)   _real [m]    # Centroid in x (5)[4]
xpher(niz)  _real [m]    # Centroid slope in x (6)[5]
yher(niz)   _real [m]    # Centroid in y (7)[6]
ypher(niz)  _real [m]    # Centroid slope in y (8)[7]
sher(niz)   _real [m]    # Position (9)[8]
vzher(niz)  _real [1]    # Axial velocity over clight (10)[9]
enxher(niz) _real [pi-m-rad] # Normalized X emittance (11)[10]
enyher(niz) _real [pi-m-rad] # Normalized Y emittance (12)[11]
cur(niz)    _real [Amps] # Current (13)[12]
dq(niz)     _real [?]    # Charge per slice (14)[13]
den(niz)    _real        # Line-charge density times clight (15)[14]
var(16,niz) _real        # Copy of envelope data, all in one place.
fviscous   real  /0.0/ [1]    # Typical shock width as a fraction of the beam length

*********** HERgeneral dump:
lherout    logical  /.true./ # Print diagnostic output
hertime    real              # Total runtime

*********** HERflags dump:
icharge    integer /1/ # The space-charge model...
		       # 0: uses simpel g-factor model with a fixed pipe/beam radius of 1.6
                       # 1: uses simple g-factor space-charge model
                       # 2: includes envelope variation in space-charge model
                       # 7: uses Warp's RZ solver to find the longitudinal field
lezbeam    logical /.true./ # Turns on use of axial self-field
lperveance logical /.true./ # Turns on use of perveance (transverse self-field)
lemittance logical /.true./ # Turns on use of emittance
lallez     logical /.true./ # Turns on use of all axial fields
lezcenter  logical /.false./ # If true, the Ez field at the beam center is used
			     # If false, Ez is averaged transversely
			     # For icharge == 7 only
lviscous   logical /.false./ # If true, artificial viscosity is added			     
lsavehist  logical /.true./ # Turns on saving of history of envelope
iprofile   integer /0/ # line-charge profile flag
                       # 0 uses hyperbolic-tangent profile
                       # 1 uses flat-top profile with linear fall-off
                       # 2 uses flat-top profile with quadratic fall-off
                       # 3 uses flat-top profile with cubic fall-off
islice     integer /0/ # reference slice number.
                       # the location of this slice is compared to zuher to determine if the run has finished.
nsteps     integer /0/ # if nonzero, nsteps steps are taken; ignored otherwise

*********** HERfield:
# Fields and forces at the current time step
nizfield             integer # Size of arrays
dedx(nizfield)       _real # From quadrupole
dbdx(nizfield)       _real # From quadrupole
ezbeam(nizfield)     _real [V/m] # Self axial electric field
rpipe(nizfield)      _real # Pipe radius at the z-location of the different slices

*********** HERvartmp:
# These are used as temporary space in the integration of the envelope
# equations.
nizvartmp         integer
dy(16,nizvartmp)  _real

*********** HERtmp:
niztmp           integer
gtemp(niztmp)    _real # Temp space for getselffield, G-factor
denmid(niztmp)   _real # Temp space for getselffield
dden(niztmp)     _real # Temp space for getselffield
denv(niztmp)     _real # Temp space for getselffield
rad(niztmp)      _real # Temp space for getselffield
deval(niztmp)    _real # Temp space for getselffield, Bessel model

*********** HERbessel:
nbessel           integer /20/ # Number of terms to use in the Bessel expansion
besselzero(1024)   _real # Zeros of the Bessel function J0
besselfactor(1024) _real # (J1(x))^(-2) in which x is a zero of the Bessel
			 # function J0

*********** HERhist dump:
lhher integer
nhher integer /1/
nizhist  integer # Number of slices whose history is saved
jhher integer /-1/
hther(0:lhher)        _real [s]   limited(0:jhher)     # Time
haher(nizhist,0:lhher)   _real [m]   limited(nizhist,0:jhher)
hapher(nizhist,0:lhher)  _real [m]   limited(nizhist,0:jhher)
hbher(nizhist,0:lhher)   _real [m]   limited(nizhist,0:jhher)
hbpher(nizhist,0:lhher)  _real [m]   limited(nizhist,0:jhher)
hxher(nizhist,0:lhher)   _real [m]   limited(nizhist,0:jhher) # X centroid in x
hxpher(nizhist,0:lhher)  _real [m]   limited(nizhist,0:jhher) # X centroid slope
hyher(nizhist,0:lhher)   _real [m]   limited(nizhist,0:jhher) # Y centroid in y
hypher(nizhist,0:lhher)  _real [m]   limited(nizhist,0:jhher) # Y centroid slope
hsher(nizhist,0:lhher)   _real [m]   limited(nizhist,0:jhher) # Position
hvzher(nizhist,0:lhher)  _real [m/s] limited(nizhist,0:jhher) # Axial velocity
henxher(nizhist,0:lhher) _real []    limited(nizhist,0:jhher) #
henyher(nizhist,0:lhher) _real []    limited(nizhist,0:jhher) #
hcur(nizhist,0:lhher)    _real []    limited(nizhist,0:jhher) #
hdq(nizhist,0:lhher)     _real []    limited(nizhist,0:jhher) #
hden(nizhist,0:lhher)    _real []    limited(nizhist,0:jhher) #

*********** HERsubs:
#  Callable subroutines in the HER package
herinit()  subroutine
hergen()   subroutine
herexe()   subroutine
herx ()    subroutine #  BASIS-level interface to HERMES
herrun(niz:integer,y:real,dsher:real,nther:integer,
       aion:real,zion:real,zlher:real,zuher:real,fviscous:real,islice:integer,
       nsteps:integer,icharge:integer,lezbeam:logical,lperveance:logical,
       lemittance:logical,lallez:logical,lezcenter:logical,lviscous:logical,
       lsavehist:logical,lfail:logical) subroutine
  # Runs HERMES kernel
getappliedfield(t:real,dt:real,y:real,niz:integer) subroutine
  # Gather applied fields
getderivs(t:real,dt:real,y:real,dy:real,niz:integer,
       aion:real,zion:real,fviscous:real,
       icharge:integer,lezbeam:logical,lperveance:logical,lemittance:logical,
       lallez:logical,lezcenter:logical,lviscous:logical,lfail:logical) subroutine
  # Calculates new envelope quantities based using applied and self fields
onestep(y:real,dy:real,niz:integer,t:real,dt:real,
        aion:real,zion:real,fviscous:real,
        icharge:integer,lezbeam:logical,lperveance:logical,
        lemittance:logical,lallez:logical,lezcenter:logical,
        lviscous:logical,lfail:logical) subroutine
  # One complete isochronous leapfrog integration step
getselffield(niz:integer,y:real,eval:real,rpipe:real,icharge:integer,
          lezbeam:logical,lezcenter:logical,lfail:logical) subroutine
  # Calculates self axial field
gethertmp(y:real,niz:integer,rpipe:real,icharge:integer,lfail:logical) subroutine
  # Calculates the variables in the HERtmp group
savehisther(t:real,y:real,niz:integer,nther:integer,lsavehist:logical)
    subroutine
  # Checks if history needs to be save and saves it
savehisther1(t:real,y:real,niz:integer) subroutine
  # Copies data into history arrays (experts only)
readbessel() subroutine
  # Initialized the arrays besselzero and besselfactor.
