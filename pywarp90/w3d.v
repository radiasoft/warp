w3d
#@(#) File W3D.V, version $Revision: 3.19 $, $Date: 2001/08/14 20:32:09 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package W3D of code WARP
# 3D - PIC package of 3d particle code
# Alex Friedman, LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194

*********** W3Dversion:
# Quantities associated with version control 
versw3d character*19 /"$Revision: 3.19 $"/ # Current code version, set by CVS

*********** InPltCtl3d dump:
# Controls for when the various plots are made
# Elements 0,1,2,3,... control plots for corresponding window
icrhoxy(0:NWINDOWS) integer  # rho contours zwindows
icrhozx(0:NWINDOWS) integer  # rho contours ywindows
icrhozy(0:NWINDOWS) integer  # rho contours xwindows
icphixy(0:NWINDOWS) integer  # phi contours zwindows
icphizx(0:NWINDOWS) integer  # phi contours ywindows
icphizy(0:NWINDOWS) integer  # phi contours xwindows
icrhoxy4            integer  # rho contours 4 zwindows / page 
icrhozx4            integer  # rho contours 4 ywindows / page
icrhozy4            integer  # rho contours 4 xwindows / page
icphixy4            integer  # phi contours 4 zwindows / page
icphizx4            integer  # phi contours 4 ywindows / page
icphizy4            integer  # phi contours 4 xwindows / page

*********** InGen3d dump:
# General parameters which control the mechanics of the run (input qtys)
filt(5,3)            real  [1]  /15*0./    # Filtering coeffs for fieldsolve
rwallfac             real  [1]  /1./       # factor for g factor in w3dgen
bndftol              real  [1]  /1.e-6/    # relative tol for b-s-f iteration
bndfflag             logical    /.true./   # flag, bent-self-field iteration
bnprflag             logical    /.true./   # flag, printing b-s-f iterations
bnezflag             logical    /.true./   # flag, bend Ez correction
bnjtflag             logical    /.true./   # flag, bend "jump term" correction
bndfitmx             integer    /15/       # maximum number of b-s-f iterations
pipeshpe             character*8 /" "/
      # Shape of pipe for Kz capacity matrix (default circle; "quad" hyperbola)
      # Shape of pipe for 3D capacity matrix (default circle; "hyp" hyperbola)
l2symtry             logical    /.false./  # Turns on 2-fold symmetry
l4symtry             logical    /.false./  # Turns on 4-fold symmetry
lbeforefs  logical    /.false./  # Turns on call to basis function "beforefs"
lafterfs   logical    /.false./  # Turns on call to basis function "afterfs"
izfsmin    integer    /0/ +parallel # Left boundary for partial field solve.
izfsmax    integer    /0/ +parallel # Right boundary for partial field solve.
solvergeom integer    /0/  # Geometry of field solver
XYZgeom    integer    /0/  # 3D-XYZ geometry will be used if 3Dsolver=XYZgeom
RZgeom     integer    /1/  # axisymmetric RZ geometry will be used if 3Dsolver=RZgeom

 
*********** InDiag3d dump:
lgetese3d logical /.true./ # Sets whether electrostatic-energy is calculated,
                           # the product of rho and phi.
                           # This calculation is expense so this flag was
                           # added to turn it off if desired.
lrhodia3d logical /.true./ # Sets whether rho diagnostics are done,
                           # calculating rhomid and rhomax along the axis.
                           # This calculation is expense so this flag was
                           # added to turn it off if desired.
lpltfld3d logical /.false./ # When true, the compiled pltfld3d calls the
                            # interpreter function 'pltfld3d' with the
                            # appropriate arguments and then returns.

*********** InPart3d dump:
# Particle input quantities (input qtys)
zjig                      real  [1]    /0./
   # Controls "jiggling" of initial positions in z, for grid loading
distrbtn                  character*8  /"none"/
   # particle distribution, either "semigaus" or "K-V"
distr_l                   character*8  /"neuffer"/
   # longitudinal velocity distribution for cigar load: either "neuffer"
   # for hard edged distribution (Vlasov equilibrium), or "gaussian" for
   # gaussian distribution with same rms variation in z.
distr_t                   character*8  /"gaussian"/
   # transverse velocity distribution: either "uniform"
   # for uniform distribution, or "gaussian" for gaussian distribution.
cigarld                   logical      /.false./
   # specifies whether or not to do a cigar load (finite beam)
xrandom                   character*8 /"pseudo"/
   # random numbers used for x,y, and z. "pseudo","fibonacc","digitrev","grid"
vtrandom                  character*8
   # random numbers used for vx and vy. "pseudo" or "digitrev"
vzrandom                  character*8
   # random numbers used for vz. "pseudo" or "digitrev"
ldprfile                  character*8
   # load profile "polar", "streamls" or "stripes"
cylinder                  logical      /.false./
   # specifies whether or not to load a cylinder
hollow                    integer      /0/
   # specifies type of hollow beam(0:none, 1:linear in r^2, 2:n(r)~(h+(1-h)r^2))
hollow_h                  real     /.5/
   # Hollowness factor used for hollow=2, Note: cannot be 1.
nxstripe                  integer  /0/
   # Number of x stripes on which particles are loaded, for grid loading
nystripe                  integer  /0/
   # Number of y stripes on which particles are loaded, for grid loading
nzstripe                  integer  /0/
   # Number of z stripes on which particles are loaded, for grid loading
nfibgrps                  integer  /0/
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
dig4                      integer /7/
   # fourth base used for digit reversed loading
dig5                      integer /11/
   # fifth base used for digit reversed loading
dig6                      integer /13/
   # sixth base used for digit reversed loading
dig7                      integer /17/
   # seventh base used for digit reversed loading
dig8                      integer /23/
   # eigth base used for digit reversed loading
nrdist                    integer
   # number of data points for arbitrary particle distribution in r
rdist(0:nrdist)           _real
   # User specified arbitrary particle distribution in r
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at rmin and 1 being at rmax.
nrmrdist(0:nrdist)          _real
   # Normalized rdist. sum(nrmrdist)=1. Calculated automatically.
intrdist(0:nrdist)         _real
   # Integral of rdist.  Calculated automatically.
nvrdist                   integer
   # number of data points for arbitrary particle velocity distribution in r
vthrofr(0:nvrdist)         _real
   # User specified radial thermal velocity as a function of radius.
vrbarofr(0:nvrdist)         _real
   # User specified average radial velocity as a function of radius.
nzdist                    integer
   # number of data points for arbitrary particle distribution in z
zdist(0:nzdist)           _real
   # User specified arbitrary particle distribution in z
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at zimin and 1 being at zimax.
nrmzdist(0:nzdist)          _real
   # Normalized zdist. sum(nrmzdist)=1. Calculated automatically.
intzdist(0:nzdist)         _real
   # Integral of zdist.  Calculated automatically.
nenvofz                    integer
   # Number of data points for the axially varying envelope
aofz(0:nenvofz)            _real
   # User specified axially varying X envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
bofz(0:nenvofz)            _real
   # User specified axially varying Y envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
apofz(0:nenvofz)           _real
   # User specified axially varying X' envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
bpofz(0:nenvofz)           _real
   # User specified axially varying Y' envelope.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
xofz(0:nenvofz)            _real
   # User specified axially varying X centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
yofz(0:nenvofz)            _real
   # User specified axially varying Y centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meters.
xpofz(0:nenvofz)           _real
   # User specified axially varying X' centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
ypofz(0:nenvofz)           _real
   # User specified axially varying Y' centroid
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of radians.
nvbeamofz                  integer
   # Number of data points for the axially varying axial velocity
vbeamofz(0:nvbeamofz)      _real
   # User specified axially varying axial velocity
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at zimin and 1 being at zimax.
nvthzofz                   integer
   # Number of data points for the axially varying axial velocity spread
vthzofz(0:nvthzofz)        _real
   # User specified axially varying axial velocity spread.
   # Data is assumed to be uniformly spaced from zimin to zimax.
   # Data is assumed to have units of meter/second.
nemitofz                   integer
   # Number of data points for axially varying transverse emittances 
   # emitxofz and emityofz. 
emitxofz(0:nemitofz)       _real
   # User specified axially varying transverse x-emittances
   # Data is assumed to be uniformly spaced in a normalized range from 0 to 1
   # with 0 being at zimin and 1 being at zimax.
emityofz(0:nemitofz)       _real
   # User specified axially varying transverse y-emittances
   # with same mesh as emitxofz.  


*********** InMesh3d dump:
# Mesh specifications (input qtys)
xmmin  real  [m]  /0./            #  Lower limit of mesh
xmmax  real  [m]  /0./            #  Upper limit of mesh
ymmin  real  [m]  /0./            #  Lower limit of mesh
ymmax  real  [m]  /0./            #  Upper limit of mesh
zmmin  real  [m]  /0./ +parallel  #  Lower limit of mesh
zmmax  real  [m]  /0./ +parallel  #  Upper limit of mesh
nx     integer    /2/             #  Mesh points are 0,...,nx
ny     integer    /2/             #  Mesh points are 0,...,ny
nz     integer    /2/  +parallel  #  Mesh points are 0,...,nz

*********** Damped_eom dump:
# All quantities associated with the damped mover algorithm
eomdamp                   real   [1]  /0./     # EOM Damping param and switch 
itdamp                    integer     /5/      # Timestep when damping starts
npdamp                    integer # Length of particle arrays for damping qtys
exold(npdamp)             _real   # Electric field at particle one step back
eyold(npdamp)             _real   # Electric field at particle one step back
ezold(npdamp)             _real   # Electric field at particle one step back
exlag(npdamp)             _real   # Electric field, lag-averaged for damping
eylag(npdamp)             _real   # Electric field, lag-averaged for damping
ezlag(npdamp)             _real   # Electric field, lag-averaged for damping

*********** Fields3d:
# Large arrays: potential and charge density, plot scratch, etc.
bndferr                real    [1]     # Error after bent-self-field solve
bndfit                 integer [1]     # Iters used, bent-self-field solve
xmesh(0:nx)            _real [m] +dump # X coordinates of mesh points
ymesh(0:ny)            _real [m] +dump # Y coordinates of mesh points
zmesh(0:nz)            _real [m] +dump +parallel # Z coordinates of mesh points
nmxy                   integer /0/ +dump # larger of nx, ny
nmxyz                  integer /0/ +dump +parallel # largest of nx, ny, nz
nzfull                 integer /0/ +dump # Full size of nz
izextra                integer /1/ +dump # Amount of extra space at end of phi
scrtch(0:nmxyz,0:nmxy) _real           # Scratch for fieldsolve, plots
phi(0:nx,0:ny,-1:nz+izextra) _real [V] +parallel # Electrostatic potential
rho(0:nx,0:ny,0:nz)    _real [C/m**3] +parallel # Charge density
attx(0:nx-1)           _real           # Attenuation factor as fcn. of kx
atty(0:ny-1)           _real           # Attenuation factor as fcn. of ky
attz(0:nzfull)         _real           # Attenuation factor as fcn. of kz
kxsq(0:nx-1)           _real [1/m**2]  # Discrete analog to kx^2/4Pi
kysq(0:ny-1)           _real [1/m**2]  # Discrete analog to ky^2/4Pi
kzsq(0:nzfull)         _real [1/m**2]  # Discrete analog to kz^2/4Pi
rstar(-1:nz+1)         _real [m] +dump # Radius of curv of refrnce orbit
phiprv(0:nx,0:nz)      _real [V]       # Prev phi at y_mid, for error test
phisav(0:nx,-1:nz)     _real [V]       # Phi at current y slice (scratch) 
xywork(0:nx,0:ny)      _real           # Work space for transverse FFTs
zwork(2,0:nx,0:nzfull) _real           # Work space used to optimize vsftz

*********** Efields3d:
nx_selfe integer /0/ +dump           # Same as nx
ny_selfe integer /0/ +dump           # Same as ny
nz_selfe integer /0/ +dump +parallel # Same as nz
selfe(3,0:nx_selfe,0:ny_selfe,0:nz_selfe) _real [V/m] # Self E field,
 # calculated from phi via finite difference. Only used when top.efetch = 3

*********** Picglb3d dump:
# Globally useful quantities for PIC simulation
dx                        real   [m]  /0./  #  mesh spacing in x
dy                        real   [m]  /0./  #  mesh spacing in y
dz                        real   [m]  /0./  #  mesh spacing in z
nxyz                      integer /0/
   # size of a field array, (nx+1)*(ny+1)*(nz+1)
ix_axis                   integer [1]       # x location of axis in mesh
iy_axis                   integer [1]       # y location of axis in mesh
iz_axis                   integer [1]       # z location of axis in mesh

*********** InjectVars3d dump:
inj_nx               integer  # size injection arrays in x
inj_ny               integer  # size injection arrays in y
inj_ninj             integer  # Auto set to either 1 or ninject. Set to 1
                              # when no emitting source are within two grid
                              # cells of each other.
inj_ns               integer  # Auto set to either 1 or ns. Set to 1
                              # when only one species is being injected from
                              # each source.
inj_d                real /1./
   # Distance from surface where phi is fetched. In units of dz.
inj_f                real /1./
   # Scaling factor on the number of particles injected.
inj_zstart           real /0./
   # Starting location relative to the emitting surface location.
inj_grid(0:inj_nx,0:inj_ny) _real
   # Grid giving axial field grid location of injection sources
inj_angl(0:inj_nx,0:inj_ny) _real
   # Grid giving angle of injection sources for each transverse location
inj_np(0:inj_nx,0:inj_ny,inj_ninj,inj_ns)   _real
   # Grid holding number of particles injected on the current time step.
inj_prev(0:inj_nx,0:inj_ny,inj_ninj,inj_ns) _real
   # Grid holding number of particles injected on previous time step.
inj_id(0:inj_nx,0:inj_ny)   _integer
   # Grid holding the number of the injection source for each transvere cell
inj_rho(0:inj_nx,0:inj_ny,inj_ninj)  _real
   # Surface charge density at the emitting surface.
inj_phi(0:inj_nx,0:inj_ny)  _real
   # Electrostatic potential at the emitting surface.
inj_area(0:inj_nx,0:inj_ny,inj_ninj)  _real
   # Inverse of the fraction of the grid cell's area within the emitting surface

*********** Setpwork3d:
# Scratch arrays for subroutine setptcls
npgrp             integer /0/
indx(npgrp)      _integer
xt(npgrp)        _real
yt(npgrp)        _real
zt(npgrp)        _real
rt(npgrp)        _real
tt(npgrp)        _real
uxt(npgrp)       _real
uyt(npgrp)       _real
uzt(npgrp)       _real
perpscal(npgrp)  _real
at(npgrp)        _real
apt(npgrp)       _real
bt(npgrp)        _real
bpt(npgrp)       _real
xct(npgrp)       _real
xpct(npgrp)      _real
yct(npgrp)       _real
ypct(npgrp)      _real

*********** Multipole dump:
# Electrostatic multipole moments of the electrostatic potential
nmom                  integer  /1/  # number of terms in multipole expansion  
lazi(nmom)            _integer      # azimuthal mode numbers of moments 
irpow(nmom)           _integer      # radial powers of moments 
rmomcex(nmom,0:nz)    _real    [1]  # multipoles, cos(theta) expansion terms 
rmomsex(nmom,0:nz)    _real    [1]  # multipoles, sin(theta) expansion terms

*********** Apertures dump:
# This group contains data describing the size and location of a circular
# aperture which is used to scrape particles. At the aperture, the transverse
# E field is calculated seperately, to explicitly include the conductor.
napertures integer /0/          # Number of special apertures
aper_zmax  integer /1/          # Maximum z length of the apertures
aper_zs(napertures)   _real [m] # Z start in lab frame of special apertures.
aper_ze(napertures)   _real [m] # Z end in lab frame of special apertures.
aper_rad(napertures)  _real [m] # Radius of special apertures.
aper_x(napertures)    _real [m] # X of special apertures center.
aper_y(napertures)    _real [m] # Y of special apertures center.
aper_volt(napertures) _real [V] # Voltage on special apertures.
aper_ex(0:nx,0:ny,-1:aper_zmax,napertures) _real [V/m]
                            # Calculated Ex on special apertures.
aper_ey(0:nx,0:ny,-1:aper_zmax,napertures) _real [V/m]
                            # Calculated Ey on special apertures.

*********** W3Dsubs:
# Subroutines in package 3D
w3dgen() subroutine
w3dexe() subroutine
divxy(iz,ndiv,divx:real,divy:real,divvx:real,divvx2:real,divvy:real,
      divvy2:real,wnpx:real,wnpy:real,itask)
             subroutine # calculates RMS vx and vy versus x and y
exteb3d(np,xp:real,yp:real,zp:real,uzp:real,gaminv:real,dtl:real,dtr:real,
        bx:real,by:real,bz:real,ex:real,ey:real,ez:real,
        m:real,q:real,bendres:real,bendradi:real,gammabar:real,dt:real)
             subroutine # Sets external E and B fields
othere3d(np,xp:real,yp:real,zp:real,zbeam:real,zimax:real,zimin:real,
         straight:real,ifeears,eears:real,eearsofz:real,dzzi:real,nzzarr,
         zzmin:real,dedr:real,dexdx:real,deydy:real,dbdr:real,
         ex:real,ey:real,ez:real,
         bx:real,by:real,bz:real)
             subroutine # Sets external E field
getese3d()   subroutine # Computes electrostatic energy
gtlchg3d()   subroutine # Computes line charge density
inject3d(itask:integer)
             subroutine # Injection routine
injctint()   subroutine # Initialization for injection
fill_inj(dx:real,dy:real,dz:real,ix_axis:integer,iy_axis:integer)
             subroutine # Initializes arrays describing the geometry of the
                        # emitting surface. Automatically called.
inj_sete3d(phi:real,np:integer,xp:real,yp:real,zp:real,zgrid:real,
           xmmin:real,ymmin:real,zmmin:real,dx:real,dy:real,dz:real,
           nx:integer,ny:integer,nz:integer,ex:real,ey:real,ez)
             subroutine # Calculate the E field for particles near the
                        # emitting surface.
loadrho3d(ins,nps,is,lzero:logical) 
             subroutine # Provides a simple interface to the charge density
                        # loading routine setrho3d
padvnc3d(center:string)
             subroutine # Advances particles and rho
perphi3d(phi:real,nx,ny,nz)
             subroutine # Equates end slices of phi for periodicity
perrho3d(rho:real,nx,ny,nz,periinz:logical)
             subroutine # Sums end slices of rho for periodicity
prntpa3d(lprntpara:logical)
             subroutine # Prints out 3d specific stuff (like prntpara())
bendez3d(np,xp:real,zp:real,ez:real,bendres:real,bendradi:real,
         bends:logical,bnezflag:logical,linbend:logical)
             subroutine #  Corrects axial electric field for warped geometry
zbendcor(np, xp:real, zp:real, uxp:real,uzp:real, gaminv:real, ddt:real,
         bendres:real,bendradi:real, bends:logical,linbend:logical)
             subroutine # Applies correction to z-advance for bends
epush3d(np,uxp:real,uyp:real,uzp:real,ex:real,ey:real,ez:real,q:real,m:real,
        dt:real)
             subroutine # Particle velocity advance from E field
bpush3d (np,uxp:real,uyp:real,uzp:real,gaminv:real,bx:real,by:real,bz:real,
         q:real,m:real,dt:real,bpush:real)
             subroutine # Particle velocity advance from B field
xpush3d (np,xp:real,yp:real,zp:real,uxp:real,uyp:real,uzp:real,gaminv:real,
         dt:real)
             subroutine # Particle position advance
seteears ()  subroutine # Sets eearsofz, the axial confining field
sete3d (phi1d:real,selfe:real,np,xp:real,yp:real,zp:real,zgrid:real,
        xmmin:real,ymmin:real,zmmin:real,dx:real,dy:real,dz:real,nx,ny,nz,
        efetch:integer,ex:real,ey:real,ez:real,
        l2symtry:logical,l4symtry:logical)
             subroutine # Sets internal E field
getselfe3d(phi:real,nx:integer,ny:integer,nz:integer,
           selfe:real,nx_selfe:integer,ny_selfe:integer,nz_selfe:integer,
           dx:real,dy:real,dz:real)
             subroutine # Calculates the self-E via finite difference of phi
setrho3d (rho1d:real,np,xp:real,yp:real,zp:real,zgrid:real,uzp:real,q:real,
          wght:real,depos:string)
             subroutine # Computes charge density
sezax3d()    subroutine # Sets EZAX, Ez on axix
sphiax3d()   subroutine # Sets PHIAX, E. S. potential on axis
srhoax3d()   subroutine # Sets RHOAX, charge density on axis
rhodia3d()   subroutine # Sets rhomid and rhomax diagnostics
stckxy3d(np,xp:real,xmmax:real,xmmin:real,dx:real,yp:real,ymmax:real,ymmin:real,
         dy:real,zp:real,zmmin:real,dz:real,uxp:real,uyp:real,uzp:real,
         zgrid:real,zbeam:real,nzzarr,prwallz:real,prwallxz:real,prwallyz:real,
         prwelips:real,
         l2symtry:logical,l4symtry:logical,lostpars:real,zzmin:real,dzzi:real)
             subroutine # Enforces sticky x and y walls
stptcl3d()   subroutine # Particle initializer
setrstar(rstar:real,nz:integer,dz:real,zmmin:real,zgrid:real)
             subroutine # Loads radius of reference orbit into rstar array 
fieldsol3d(iwhich) subroutine # Bent-self-field iterative solver
vp3d(iwhich) subroutine # The 3d Poisson solver

pltfld3d(fld:string,freqflag:integer)
             subroutine # Controls field plots

multpole(lmod:integer,nlmod:integer,irpowmx:integer,
         lcosex:logical,lsinex:logical,aper:real,xcen:real,ycen:real,
         nmult:integer,nres:integer,tol:real) 
            subroutine # calculate the multipole moments of the potential
inj_smoother(nx:integer,ny:integer,inj_phi:real,dx:real,dy:real,
             xmmin:real,ymmin:real,x0:real,y0:real,a0:real,b0:real,
             inj_nsmooth:integer) subroutine
getinj_phi(nx:integer,ny:integer,nz:integer,phi:real,dx:real,dy:real,dz:real,
           xmmin:real,ymmin:real) subroutine




******** Subtimers:
lw3dtimesubs logical /.false./
timew3dinit real /0./
timew3dvers real /0./
timew3dgen real /0./
timew3dexe real /0./
timew3dfin real /0./
timestep3d real /0./
timeexteb3d real /0./
timeothere3d real /0./
timegetese3d real /0./
timegtlchg3d real /0./
timeseteears real /0./
timepadvnc3d real /0./
timeperphi3d real /0./
timeperrho3d real /0./
timeepush3d real /0./
timeepusht3d real /0./
timebpush3d real /0./
timebpusht3d real /0./
timexpush3d real /0./
timexpusht3d real /0./
timesete3d_relativity real /0./
timeedamp real /0./
timegetbend real /0./
timebendez3d real /0./
timezbendcor real /0./
timesete3d real /0./
timegetselfe3d real /0./
timestptcl3d real /0./
timesetrho3d real /0./
timeloadrho3d real /0./
timestckxy3d real /0./
timesetrstar real /0./
timeinject3d real /0./
timeinjctint real /0./
timefill_inj real /0./
timeinj_sete3d real /0./
timeinj_setrho3d real /0./
timesete3d_aperture real /0./
timeset_aperture_e real /0./
timefieldsol3d real /0./
timevp3d real /0./
