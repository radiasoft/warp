her
#@(#) File HER.V, version $Revision: 3.6 $, $Date: 2001/07/13 00:37:34 $
# Copyright (c) 1999, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for HERMES.
# Michiel de Hoon,    LBNL, (510) 486-5157  mdehoon@lbl.gov
{
}

*********** HERversion:
# Version control for her package
versher       character*19 /"$Revision: 3.6 $"/  # Current code version, set by SCCS

*********** HERvars dump:
# Variables needed by the package HER
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
                       # 0: uses simpel g-factor model with a fixed g-factor (initialized to 2 ln(1.6) in the generate)
                       # 1: uses simple g-factor space-charge model
                       # 2: includes envelope variation in space-charge model
                       # 3: g-factor model with end effects (placeholder; not implemented)
                       # 4: g-factor model with end effects (placeholder; not implemented)
                       # 5: Bessel expansion with uniform line charge density in each slice
                       # 6: Bessel expansion with smoothly varying line charge density in each slice
                       # 7: uses Warp's RZ solver to find the longitudinal field
lezbeam    logical /.true./ # Turns on use of axial self-field
lperveance logical /.true./ # Turns on use of perveance (transverse self-field)
lemittance logical /.true./ # Turns on use of emittance
lallez     logical /.true./ # Turns on use of all axial fields
lezcenter  logical /.false./ # If true, the Ez field at the beam center is used
                             # If false, Ez is averaged transversely
                             # For icharge == 7 only
lcurgrid   logical /.false./ # If true, calculate the line charge density and
                             # the current using the grid
                             # For icharge == 7 only
lviscous   logical /.false./ # If true, artificial viscosity is added
lsavehist  logical /.true./ # Turns on saving of history of envelope
iprofile   integer /0/ # line-charge profile flag
                       # 0 uses hyperbolic-tangent profile
                       # 1 uses flat-top profile with linear fall-off
                       # 2 uses flat-top profile with quadratic fall-off
                       # 3 uses flat-top profile with cubic fall-off
lfail      logical /.false./ # Set to true if a time step fails

*********** HERfield:
# Fields and forces at the current time step, for the different slices
nizfield             integer # Size of arrays
dedx(nizfield)       _real # From quadrupole
dbdx(nizfield)       _real # From quadrupole
er(nizfield)         _real # Radial electric field (n = 0, v = 0)
ez(nizfield)         _real # Longitudinal electric field (n = 0, v = 0)
br(nizfield)         _real # Radial magnetic field (n = 0, v = 0)
bz(nizfield)         _real # Longitudinal magnetic field (n = 0, v = 0)
# Note: er,bz,br are ignored for now
ezbeam(nizfield)     _real [V/m] # Self axial electric field
fx(nizfield)         _real # From image
fy(nizfield)         _real # From image
gxx(nizfield)        _real # From image
gxy(nizfield)        _real # From image
gyx(nizfield)        _real # From image
gyy(nizfield)        _real # From image
ipipetype            integer /0/ # Defines the pipe shape used for image force calculations
                                 # 0 round pipe (default)
                                 # 1 electric quadrupole rods (hyperbolic sufaces)
                                 # 2 dipole plates (parallel plates)
rpipe(nizfield)      _real [m] # Pipe radius at the z-location of the different slices
pipeox(nizfield)   _real /0.0/ [m] # Horizontal offset of the pipe (needed for images)
pipeoy(nizfield)   _real /0.0/ [m] # Vertical offset of the pipe (needed for images)
pipeph(nizfield)   _real /0.0/ [rad] # Rotation angle of the pipe (needed for images)

*********** HERvartmp:
# These are used as temporary space in the integration of the envelope
# equations.
nizvartmp         integer
dy(16,nizvartmp)  _real

*********** HERtmp:
niztmp           integer
gtemp(niztmp)    _real # Temp space for fieldsolhr, g-factor
denmid(niztmp)   _real # Temp space for fieldsolhr
dden(niztmp)     _real # Temp space for fieldsolhr
denv(niztmp)     _real # Temp space for fieldsolhr
rad(niztmp)      _real # Temp space for fieldsolhr
deval(niztmp)    _real # Temp space for fieldsolhr, Bessel model

*********** HERbessel:
nbessel           integer /20/ # Number of terms to use in the Bessel expansion
besselzero(1024)   _real # Zeros of the Bessel function J0
besselfactor(1024) _real # (J1(x))^(-2) in which x is a zero of the Bessel
                         # function J0

*********** HERhist dump:
lhher integer
nhher integer /1/
jhher integer /-1/
nizhist integer
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
  # Runs HERMES kernel
extebher(t:real,dt:real,y:real,niz:integer) subroutine
  # Gather applied fields
getderivs(t:real,dt:real,y:real,dy:real,niz:integer,
       aion:real,zion:real,fviscous:real,
       icharge:integer,iimage:integer,lezbeam:logical,lperveance:logical,lemittance:logical,
       lallez:logical,lezcenter:logical,lviscous:logical,lcurgrid:logical,lfail:logical) subroutine
  # Calculates new envelope quantities using applied and self fields
stephr(y:real,dy:real,niz:integer,t:real,dt:real,
        aion:real,zion:real,fviscous:real,
        icharge:integer,iimage:integer,lezbeam:logical,lperveance:logical,
        lemittance:logical,lallez:logical,lezcenter:logical,
        lviscous:logical,lcurgrid:logical,lfail:logical) subroutine
  # One complete isochronous leapfrog integration step
setrhohr(niz:integer,y:real,rpipe:real,icharge:integer,lfail:logical) subroutine 
  # Loads the charge onto the RZ grid
fieldsolhr(niz:integer,y:real,eval:real,rpipe:real,icharge:integer,
          lezbeam:logical,lezcenter:logical,lfail:logical) subroutine
  # Calculates self axial field
getradius(y:real,niz:integer,rpipe:real,icharge:integer) subroutine
  # Calculates the beam radius
getgfactor(niz:integer,rpipe:real,icharge:integer) subroutine
  # Calculates the g-factor
getcurrent(y:real,niz:integer,rpipe:real,icharge:integer,lcurgrid:logical,
           lfail:logical) subroutine
  # Calculates the current at the slices boundaries
savehermesvars(t:real,y:real,niz:integer,it:integer,lsavehist:logical)
    subroutine
  # Checks if history needs to be saved and saves it
getimage(y:real,niz:integer,iimage:integer,lfail:logical) subroutine
  # Calculates the image effects on the transverse envelope
copyhermesvars() subroutine
   # Copies the her.var variables into the separated arrays
resizehermeshist() subroutine
   # Resizes the history arrays to only those entries that are used.
