cir
#@(#) File CIR.V, version $Revision: 3.5 $, $Date: 2001/05/29 21:20:07 $
# Copyright (c) 1999, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for CIRCE.
# Dave Grote,    LLNL, (510) 423-7194  dave@hif.llnl.gov 
# Bill Sharp,    LLNL 
{
}

*********** CIRversion:
# Version control for cir package
verscir  character*19 /"$Revision: 3.5 $"/  # Current code version, set by CVS

*********** CIRvars dump:
# Variables needed by the package CIR
dscir real /0.0/  [m]    # Step size in envelope calculation
zlcir real /0.0/  [m]    # Starting z for CIRCE calc
zucir real /0.0/  [m]    # Maximum  z for CIRCE calc
nscir integer /0/        # Total number of time steps
nit       integer [1]    # Number of point along beam
acir(nit)   _real [m]    # Width in x (1)[0]
apcir(nit)  _real [1]    # Slope in x (2)[1]
bcir(nit)   _real [m]    # Width in y (3)[2]
bpcir(nit)  _real [1]    # Slope in y (4)[3]
xcir(nit)   _real [m]    # Centroid in x (5)[4]
xpcir(nit)  _real [m]    # Centroid slope in x (6)[5]
ycir(nit)   _real [m]    # Centroid in y (7)[6]
ypcir(nit)  _real [m]    # Centroid slope in y (8)[7]
tcir(nit)   _real [s]    # Time (9)[8]
vzcir(nit)  _real [1]    # Axial velocity over clight (10)[9]
enxcir(nit) _real [pi-m-rad] # Normalized X emittance (11)[10]
enycir(nit) _real [pi-m-rad] # Normalized Y emittance (12)[11]
cur(nit)    _real [Amps] # Current (13)[12]
dq(nit)     _real [?]    # Charge per slice (14)[13]
den(nit)    _real        # Line-charge density times clight (15)[14]
var(16,nit) _real        # Copy of envelope data, all in one place.

*********** CIRgeneral dump:
lcirout    logical  /.true./ # Print diagnostic output
cirtime    real              # Total runtime
footrise    real
tfoot    real
curfoot    real
tiltmid    real /0.5/
tiltlen    real /10.0/

*********** CIRflags dump:
ltdependent logical /.false./ # When true, CIRCE uses follows a beam
                              # using time as the independent variable
                              # instead of position
icharge    integer /1/ # The space-charge model...
                       # 1: uses simple g-factor space-charge model
                       # 2: includes envelope variation in space-charge model
                       # 3: includes end effects in space-charge model
                       # 4: uses a more accurate model of end effects
                       # 5: uses a Bessel-series model of end effects
lezbeam    logical /.true./ # Turns on use of axial self-field
lperveance logical /.true./ # Turns on use of perveance (transverse self-field)
lemittance logical /.true./ # Turns on use of emittance
lallez     logical /.true./ # Turns on use of all axial fields
llinear    logical /.true./ # Only includes linear turns in bends
limage     logical /.true./ # Turns on use of image fields
lendzero   logical /.true./ # If true, a zero radius at the first and last
                            # slice is kept zero
lsavehist  logical /.true./ # Turns on saving of history of envelope
icurload   integer /0/ # current-load flag
                       # 0 sets up line-charge profile as function of time
                       # 1 sets up line-charge profile as function of s
iprofile   integer /0/ # line-charge profile flag
                       # 0 uses hyperbolic-tangent profile
                       # 1 uses flat-top profile with linear fall-off
                       # 2 uses flat-top profile with quadratic fall-off
                       # 3 uses flat-top profile with cubic fall-off
                       # 4 uses flat-top profile with Gaussian fall-off
                       # 5 uses flat-top profile with Fermi-function fall-off
ivelload   integer /0/ # velocity-tilt flag
                       # 0 sets up an initial tilt in beam velocity vs time
                       # 1 sets up an initial tilt in beam velocity vs position
igrid      integer /0/ # controls initial spacing of arrival-time values
                       # 0 uniform layout
                       # 1 cell size decreasing near ends
iemload    integer /0/ # emittance-load flag
                       # 0 gives uniform emittance
                       # 1 gives uniform transverse temperature
                       # 2 gives uniform charge density
                       # 3 gives uniform beam radius
lfixed     logical /.false./ # if true, the correct expression for d(lambda)/dz is used

*********** CIRfield:
# Applid fields at the current location
nitfield    integer # Size of arrays
ex(nitfield)     _real # From dipole
ey(nitfield)     _real # From dipole
ez(nitfield)     _real # From accelation gap
ears(nitfield)   _real # From axially confining ear field
bx(nitfield)     _real # From dipole
by(nitfield)     _real # From dipole
bz(nitfield)     _real # (Unused)
dedx(nitfield)   _real # From quadrupole
dbdx(nitfield)   _real # From quadrupole
fx(nitfield)     _real # From image
fy(nitfield)     _real # From image
gxx(nitfield)    _real # From image
gxy(nitfield)    _real # From image
gyx(nitfield)    _real # From image
gyy(nitfield)    _real # From image
ezbeam(nitfield) _real [V/m] # Self axial electric field
yflderr      real
bendcurv     real /0.0/
ipipetype    integer /0/
rpipe        real
xpipeerr     real /0.0/
ypipeerr     real /0.0/
sinrot       real /0.0/
cosrot       real /1.0/

*********** CIRvartmp:
# These are used as temporary space in the integration of the envelope
# equations.
nitvartmp         integer
yt(16,nitvartmp)  _real
dy1(16,nitvartmp) _real
dy2(16,nitvartmp) _real
dy3(16,nitvartmp) _real

*********** CIRtmp:
nittmp           integer
gtemp(nittmp)    _real # Temp space for getezbeam, G-factor
curmid(nittmp)   _real # Temp space for getezbeam
denmid(nittmp)   _real # Temp space for getezbeam
dden(nittmp)     _real # Temp space for getezbeam
denv(nittmp)     _real # Temp space for getezbeam
rad(nittmp)      _real # Temp space for getezbeam
zeta(nittmp)     _real # Temp space for getezbeam
sum0(nittmp)     _real # Temp space for getezbeam
sum1(nittmp)     _real # Temp space for getezbeam
dlna(nittmp)     _real # Temp space for getezbeam
dlnb(nittmp)     _real # Temp space for getezbeam
etemp(nittmp)    _real # Temp space for getezbeam
ftemp(nittmp)    _real # Temp space for getezbeam
ezval(nittmp)    _real # Temp space for setears
dummy(nittmp)    _real # Temp space for setears

*********** CIRmatch:
varfract            real /0.05/
varstepmin          real /1.0e-08/
nitmatch            integer
varold(8,nitmatch) _real
varnew(8,nitmatch) _real
del(8)              real
der(8,8,nitmatch)  _real
delvar(8,nitmatch) _real
ratio(nitmatch)    _real

*********** CIRbessdata dump:
nmax integer /150/
#besszero real /2.408255577/
besszeros(10) real /2.4048255577, 5.5201, 8.6537, 11.7915, 14.9309, 18.0711, 21.2116, 24.3525, 27.4935, 30.6346/
bessdenom(10) real /1.5587, 3.5280, 5.5182, 7.5134, 9.5106, 11.5088, 13.5075, 15.5065, 17.5058, 19.5052/

*********** CIRhist dump:
lhcir integer
nhcir integer
nithist  integer # Number of slices whose history is saved
jhcir integer /-1/
hscir(0:lhcir)        _real [m]   limited(0:jhcir)     # Beam location
hacir(nithist,0:lhcir)   _real [m]   limited(nithist,0:jhcir)
hapcir(nithist,0:lhcir)  _real [m]   limited(nithist,0:jhcir)
hbcir(nithist,0:lhcir)   _real [m]   limited(nithist,0:jhcir)
hbpcir(nithist,0:lhcir)  _real [m]   limited(nithist,0:jhcir)
hxcir(nithist,0:lhcir)   _real [m]   limited(nithist,0:jhcir) # X centroid in x
hxpcir(nithist,0:lhcir)  _real [m]   limited(nithist,0:jhcir) # X centroid slope
hycir(nithist,0:lhcir)   _real [m]   limited(nithist,0:jhcir) # Y centroid in y
hypcir(nithist,0:lhcir)  _real [m]   limited(nithist,0:jhcir) # Y centroid slope
htcir(nithist,0:lhcir)   _real [s]   limited(nithist,0:jhcir) # Time
hvzcir(nithist,0:lhcir)  _real [m/s] limited(nithist,0:jhcir) # Axial velocity
henxcir(nithist,0:lhcir) _real []    limited(nithist,0:jhcir) #
henycir(nithist,0:lhcir) _real []    limited(nithist,0:jhcir) #
hcur(nithist,0:lhcir)    _real []    limited(nithist,0:jhcir) #
hdq(nithist,0:lhcir)     _real []    limited(nithist,0:jhcir) #
hden(nithist,0:lhcir)    _real []    limited(nithist,0:jhcir) #

*********** CIRsetacc:
# Temporary arrays for setacc
nitacc integer
ct0(nitacc)      _real
ctold(nitacc)    _real
ctnext(nitacc)   _real
betinit(nitacc)  _real
betaold(nitacc)  _real
betanext(nitacc) _real
curinit(nitacc)  _real
sumslice(nitacc) _real
tout(nitacc)     _real
accout(nitacc)   _real
earout(nitacc)   _real
betamid_set       real
deltbeam_set      real
iaccl             integer

*********** CIRsubs:
#  Callable subroutines in the CIR package
cirinit()  subroutine
cirgen()   subroutine
cirexe()   subroutine
cirx ()    subroutine #  BASIS-level interface to CIRCE
initbeam() subroutine # Calls setbeam using WARP database
cirrun(nit:integer,y:real,zlcir:real,zucir:real,dscir:real,nscir:integer,
       aion:real,zion:real,
       icharge:integer,lezbeam:logical,lperveance:logical,lemittance:logical,
       lallez:logical,llinear:logical,limage:logical,lendzero:logical,
       lsavehist:logical,lfixed:logical) subroutine
  # Runs CIRCE kernel
setstep(s:real,ds:real,zlcir:real,zucir:real,dscir:real) subroutine
  # Sets step size, looking for beginning of the next element
extebcir(z:real,dz:real,y:real,nit:integer) subroutine
  # Gather applied fields
derivs(s:real,ds:real,dsfrac:real,y:real,dy:real,nit:integer,
       aion:real,zion:real,
       icharge:integer,lezbeam,lperveance,lemittance,lallez,llinear,
       limage:logical,lendzero:logical,lfixed:logical) subroutine
  # Calculates new envelope quantities based using applied and self fields
integrt(y:real,yt:real,dy1:real,dy2:real,dy3:real,nit:integer,s:real,ds:real,
        aion:real,zion:real,
        icharge:integer,lezbeam,lperveance,lemittance,lallez,llinear,
        limage:logical,lendzero:logical,lfixed:logical) subroutine
  # Complete Runga-Kutta integration step
getimage(y:real,nit:integer,llinear,limage:logical,lendzero:logical) subroutine
  # Calculate images fields from surrounding pipe
getezbeam(nit:integer,y:real,eval:real,rpipe:real,icharge:integer,
          lezbeam:logical,lendzero:logical,lfixed:logical,lfailed:logical) subroutine
  # Calculates self axial field
savehistcir(s:real,y:real,nit:integer,nscir:integer,lsavehist:logical)
    subroutine
  # Checks if history needs to be save and saves it
savehistcir1(s:real,y:real,nit:integer) subroutine
  # Copies data into history arrays (experts only)
setcur(nit:integer,y:real,
       curmax:real,currise:real,curfall:real,deltbeam:real,footrise:real,
       tfoot:real,curfoot:real,
       icurload:integer,iprofile:integer,
       beta0:real,tilt:real,tiltmid:real,tiltlen:real,
       ivelload:integer) subroutine
  # Loads a current (and velocity) based on the given profile and tilts.
  # Calls setbeta to load the velocity
setbeta(beta0:real,y:real,nit:integer,
        tilt:real,tiltmid:real,tiltlen:real,
        ivelload:integer) subroutine
  # Loads a velocity based on the given tilt
settval(nit:integer,y:real,deltbeam:real,time0:real,igrid:integer) subroutine
  # Loads the initial time values

getphase(sig0:real,sig1:real,hlp:real,betaval:real,curval:real,tval:real,
         sval:real,emitx:real,emity:real,aion:real,zion:real,
         lperveance:logical) subroutine
  # Calculates sig0 and sig at one slice
getphaset(nit:integer,sig0:real,sig1:real,hlp:real,sval:real,y:real,
          aion:real,zion:real,lperveance:logical) subroutine
  # Calculates sig0 and sig at each slice
setears(nit:integer,y:real,ntval:integer,earval:real,tval:real,z:real,
        iearset:integer,iearform:integer,icharge:integer,ieargap:integer,
        lendzero:logical,lfailed:logical) subroutine
  # Calculates the ear field
setemit(nit:integer,y:real,iemload:integer,lperveance:logical) subroutine
  # Sets emittance
getlen(nit:integer,y:real,beamlen:real) subroutine
  # Get beam length
gettemp(nit:integer,y:real,s:real,beamtemp:real,delv:real,aion:real,zion:real,
        icharge:integer,lezbeam:logical,lperveance:logical,lemittance:logical,
        lallez:logical,llinear:logical,limage:logical,lendzero:logical,lfixed:logical)
        subroutine
  # Get beam temperature
resetbeta(nit:integer,y:real,inout:string) subroutine
  # Adjust axial velocity in entrance and exit of a dipole
setbeam(nit:integer,y:real,a0:real,b0:real,ap0:real,bp0:real,
        x0:real,y0:real,xp0:real,yp0:real,
        ibeam:real,currise:real,curfall:real,deltbeam:real,
        footrise:real,tfoot:real,curfoot:real,
        beta:real,gammabar:real,tilt:real,tiltmid:real,tiltlen:real,
        emitx:real,emity:real,time:real,zbeam:real,aion:real,zion:real,
        icurload:integer,iprofile:integer,ivelload:integer,igrid:integer,
        iemload:integer,lperveance:logical,lendzero:logical) subroutine
  # Loads all of the data needed for a beam.
setacc(schedule:string,nit:integer,y:real,ekinstart:real,ekinend:real,
       zlcir:real,zucir:real,deltbeam:real,curmax:real,currise:real,
       curfall:real,footrise:real,tfoot:real,curfoot:real,icurload:integer,
       iprofile:integer,beta0:real,tilt:real,tiltmid:real,tiltlen:real,
       ivelload:integer,lperveance:logical,aion:real,zion:real) subroutine
  # Sets accelerating schedules for the gaps



