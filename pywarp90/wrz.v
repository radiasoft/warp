wrz
#@(#) File WRZ.V, version $Revision: 3.8 $, $Date: 2002/11/05 22:36:44 $ 
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package WRZ of code WARP
# RZ - PIC package of rz particle code
# Alex Friedman, LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194
# Debbie Callahan, LLNL, (510)423-5926
 
*********** WRZversion dump:
# Quantities associated with version control 
verswrz character*19 /"$Revision: 3.8 $"/ # Current code version, set by CVS
 
*********** InPltCtlrz dump:
# Controls for when the various plots are made
icrhorz     integer /SELDOM/ # rho contours 
icphirz     integer /SELDOM/ # phi contours 
ipphivz     integer /SELDOM/ # phi on axis and Vz vs Z for videos
ipphivzc    integer /SELDOM/ # phi on axis and Vz vs Z (color)

*********** InGenrz dump:
# General parameters which control the mechanics of the run (input qtys)
filt(5)                   real  [1]  /5*0./   # filtering coeffs for fieldsolve
rwallfac                  real  [1]  /1./     # factor for g factor in wrzgen
   # Filtering coefficients for fieldsolve
fftdiag                logical /.false./
vecrho                    logical /.false./# Controls use of vector depos.
frcetype                  character*8 /"regress"/
    # Method used to calculate external radial force--regress or analytic
eta                       real  [V/A/m] /0.0/ # Surface resistivity (Ohm/meter)
taurc                     real  [s]     /0.0/ # RC time
lbeforefs  logical    /.false./  # Turns on call to basis function "beforefs"
lafterfs   logical    /.false./  # Turns on call to basis function "afterfs"
 
*********** InPartrz dump:
# Particle input quantities (input qtys)
distrbtn                  character*8  /"K-V"/
   # particle distribution, either "semigaus" or "K-V"
distr_l                   character*8  /"neuffer"/
   # longitudinal velocity distribution for cigar load: either "neuffer"
   # for hard edged distribution (Vlasov equilibrium), or "gaussian" for
   # gaussian distribution with same rms variation in z.
cigarld                   logical      /.false./
   # specifies whether or not to do a cigar load (finite beam)
xrandom                   character*8  /"pseudo"/
   # random numbers used for x,y, and z. "pseudo","fibonacc","digitrev","grid"
vtrandom                  character*8  /"pseudo"/
   # random numbers used for vx and vy. "pseudo" or "digitrev"
vzrandom                  character*8  /"pseudo"/
   # random numbers used for vz. "pseudo" or "digitrev"
ldprfile                  character*8  /"streamls"/
   # load profile "streamls" or "stripes"
cylinder                  logical      /.false./
   # specifies whether or not to load a cylinder
hollow                    integer    /0/
   # specifies type of hollow beam (0:none, 1:linear in r^2, 2:n(r)~(1+(h-1)r^2))
hollow_h                  real       /.5/
   # Hollowness factor used for hollow=2, Note: cannot be 1.
#shell                     logical      /.false./
   # specifies whether or not to load a shell
   # variable is kept in Fieldsrz 
kzpert                    real /0.0/
   # kz perturbation
vzpert                    real /0.0/
   # uz perturbation
phipert                   real /0.0/
   # phase shift for perturbation
zjig                      real  [1]    /0./
   # Controls "jiggling" of initial positions in z, for grid loading
nxstripe                  integer  /20/
   # Number of x stripes on which particles are loaded, for grid loading
nystripe                  integer  /20/
   # Number of y stripes on which particles are loaded, for grid loading
nzstripe                  integer  /24/
   # Number of z stripes on which particles are loaded, for grid loading
nfibgrps                  integer  /4/
   # Number used fibonacci groups
fibg1                     integer  /1958/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
fibg2                     integer  /252/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
fibg3                     integer  /414/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
fibg4                     integer  /0/
   # Number used for fibonacci loading (See S.K. Zaremba ref. in setptcls)
dig1                      integer /2/
   # first base used for digit reversed loading
dig2                      integer /3/
   # second base used for digit reversed loading
dig3                      integer /5/
   # third base used for digit reversed loading
dig4                      integer /11/
   # fourth base used for digit reversed loading
dig5                      integer /17/
   # fifth base used for digit reversed loading
 
*********** InMeshrz dump:
# Mesh specifications (input qtys)
nr                        integer /16/    #  Mesh points are 0,...,nr
nz                        integer /32/    #  Mesh points are 0,...,nz
rmmin                     real   [m]  / .00/   #  Lower limit of mesh
rmmax                     real   [m]  / .05/   #  Upper limit of mesh
zmmin                     real   [m]  /-1.0/   #  Lower limit of mesh
zmmax                     real   [m]  / 1.0/   #  Upper limit of mesh
 
*********** Sortrz dump:
ignmax                    integer /0/     #  Number of groups-psort sets it
npsort                    integer /0/     # number of particles-set in wrzgen
igrp(npsort)              _integer
ipgrp(npsort)             _integer
npic(0:nz)                _integer
nig(ignmax)               _integer
ipnt(ignmax+1)            _integer
ipindex(ignmax)           _integer
index(0:3,nz/2+1)         _integer
s(0:3,nz/2+1)             _real

*********** Fieldsrz dump:
# Large arrays: potential and charge density, plot scratch, etc.
ibc                  character*8 /"ereq0"/ # Boundary condition for fieldsolve
                                           # ereq0--E_r=0 for r > rmmax
                                           # bessel--E_r~K_1(kr) for r > rmmax
eoffrz                    real   [V/m]     
vphase                    real   [m/s]     # Calculated phase velocity of wave
rmesh(0:nr)               _real  [m]       # R coordinates of mesh points
rhalf(0:nr)               _real  [m]       # R coordinate on half mesh point
zmesh(0:nz)               _real  [m]       # Z coordinates of mesh points
scrtch(0:nz,0:nr)         _real            # Scratch for plots
phi(0:nr,-1:nz+1)         _real  [V]       # Electrostatic potential
rho(0:nr,0:nz)            _real  [C/m**3]  # Charge density
    # Must come after scrtch, to protect lost particle phi lookup
rhov(0:nr,0:nz)           _real  [C/m**3]  # temporary rho
currv(0:nz)               _real            # temporary current
attz(0:nz/2)              _real            # Attenuation factor as fcn. of kz
kzsq(0:nz)                _real  [1/m**2]  # Discrete analog to kz^2/4Pi
force(0:nr,0:nz)          _real  [N]       # Radial external force
erfld(0:nr,0:nz)          _real  [V/m]     # R component of E field
ezfld(0:nr,0:nz)          _real  [V/m]     # Z component of E field
schrg(0:nz)               _real  [C]       # Surface charge for resistive wall
rfsmat(0:nr,3,0:nz)       _real            # Tridiagonal matrix for field solve
rwwork(0:nz)              _real            # Workspace for resistive wall
rwwork2(0:nz)             _real            # Workspace for resistive wall
phikold(0:nz)             _real            # FFT of phi at wall at last timestep
eearsrz(0:nr,0:nz)        _real [V/m]      # Eears over whole beam bunch
eearsave(0:nz)            _real [V/m]      # Ears saved for intermittent ears
ezext(0:nz)               _real [V/m]      # Ext. E field for resistive wall
curr0(0:nz)               _real [A]        # Current at t=0 for capacitance
phiresist(0:nz)           _real [V]        # Int I*eta for phi on axis plots

*********** LinearBeam dump:
shell                     logical      /.false./
   # specifies whether or not to load a shell
linbeam                   logical      /.false./
   # specifies whether or not to load a 'linear' beam

*********** InFluidrz:
nfrstrz                   integer  /0/     # starting time step
ihist                     integer  /0/
ikmode                    integer  /0/
nhistrz                   integer  /100/


*********** Fluidsrz dump:
lamfld(0:nz)              _complex [C/m]   # Line charge for fluid code
vzfld(0:nz)               _complex [m/s]   # Fluid velocity
e1fld(0:nz)               _complex [V/m]   # Electric Field
lamold(0:nz)              _complex [C/m]   # Line charge at last time step
ikdt(0:nz)                _complex         # i*kz*dt
arryr(0:nz)               _real            # scratch array
arryi(0:nz)               _real            # scratch array
thistrz(0:nhistrz)        _real 
lamrhist(0:nhistrz)       _real
lamihist(0:nhistrz)       _real

*********** Picglbrz dump:
# Globally useful quantities for PIC simulation
dr                        real   [m]  /0./  #  mesh spacing in r
dz                        real   [m]  /0./  #  mesh spacing in z
nmrz                      integer /0/
   # larger of nr, nz
nrz                      integer /0/
   # size of a field array, (nr+1)*(nz+1)

*********** Setpworkrz:
# Scratch arrays for subroutine setptcls
npgrp            integer /0/
indx(npgrp)     _integer
xt(npgrp)       _real
yt(npgrp)       _real
zt(npgrp)       _real
uxt(npgrp)      _real
uyt(npgrp)      _real
uzt(npgrp)      _real
perpscal(npgrp) _real
at(npgrp)       _real
apt(npgrp)      _real
bt(npgrp)       _real
bpt(npgrp)      _real
xct(npgrp)      _real
xpct(npgrp)     _real
yct(npgrp)      _real
ypct(npgrp)     _real

*********** FFTDiagrz dump:
# Arrays used for the Fourier transform diagnostic
indxk                  integer /0/
nzfft                  integer /0/
xreal(0:nzfft)         _real
ximag(0:nzfft)         _real

*********** WRZsubs:
# Subroutines in package RZ
extbrz (np,rp:real,tp:real,zp:real,uzp:real,gaminv:real,br:real,bt:real,
        bz:real)
             subroutine # Sets external B field
exterz (np,zp:real,uzp:real,zbeam:real,rp:real,ez0:real,er:real,et:real,
        ez:real,dedr:real)
             subroutine # Sets external E field
fixrhorz()   subroutine # "Fixes" charge density by 1-beta^2, curvature
geteserz()   subroutine # Computes electrostatic energy
getforce()   subroutine # Computes external force to balance space-charge force 
gtlchgrz()   subroutine # Computes line charge density
padvncrz (center:string)
             subroutine # Advances particles and rho
perphirz(phi:real,nr,nz)
             subroutine # Equates end slices of phi for periodicity
perrhorz(rho:real,nr,nz)
             subroutine # Sums end slices of rho for periodicity
prntparz(lprntpara:logical)
             subroutine # Prints out rz specific stuff (like prntpara())
pushrz (center:string,np,rp:real,tp:real,zp:real,urp:real,utp:real,
        uzp:real,gaminv:real,er:real,et:real,ez:real,br:real,bt:real,
        bz:real,frext:real,q:real,m:real,dt:real,gamadv:real)
             subroutine # Particle advance
seterz (itask,phi:real,force:real,np,rp:real,zp:real,zbeam:real,
        rmmin:real,zmmin:real,zmmax:real,dr:real,dz:real,nr,nz,
        efetch:integer,
        er:real,et:real,ez:real,frext:real,erfld:real,
        ezfld:real,eoffrz:real)
             subroutine # Sets internal E field
setrhorz (rho1d:real,np,rp:real,zp:real,zbeam:real,uzp:real,q:real,
          wght:real,depos)
             subroutine # Computes charge density
stckyrz(np,rp:real,rmmax:real,rmmin:real,dr:real,zp:real,zmmin:real,
        dz:real,urp:real,utp:real,uzp:real,zgrid:real,ipmin,it)
             subroutine # Enforces sticky r walls
stptclrz()   subroutine # Particle initializer
fieldsolrz(iwhich:integer)
             subroutine # Full Piosson solver, calls vprz.
vprz(iwhich) subroutine # The rz Poisson solver

pltfldrz(fld:string,freqflag:integer)
             subroutine # Controls field plots
setsc(schrg:real,a0:real,rmesh:real,rho:real,nr,nz,dz:real,eta:real)
             subroutine # Sets surface charge for resistive wall
sezaxrz()    subroutine # Sets E_z on axis for plotting by onedplts
sphiaxrz()   subroutine # Sets Phi on axis for plotting by onedplts
srhoaxrz()   subroutine # Sets Rho on axis for plotting by onedplts
fldrz(nstep) subroutine # 1-d fluid code
