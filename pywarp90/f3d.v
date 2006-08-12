f3d
#@(#) File F3D.V, version $Revision: 3.153 $, $Date: 2006/08/12 04:30:29 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package F3D of code WARP6
# fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
{
LARGEPOS = 1.0e+36 # This must be the same as in top.v
}

*********** F3Dversion:
versf3d character*19 /"$Revision: 3.153 $"/#  Code version version is set by CVS

*********** F3Dvars:
# Variables needed by the test driver of package F3D
lx                        real /0.7/   #  System length in x (arbitrary units)
ly                        real /1.3/   #  System length in y (arbitrary units)
lz                        real /1.6/   #  System length in z (arbitrary units)
filt(5,3)                 real /15*0./ #  Spatial filtering coefficients
ibc                       integer      #  Boundary conditions-future use
nx                        integer /2/  #  Mesh points are 0,...,nx
ny                        integer /2/  #  Mesh points are 0,...,ny
nz                        integer /2/  #  Mesh points are 0,...,nz
a(0:nx,0:ny)              _real        #  2d charge density, potential array
b(0:nx,0:ny,0:nz)         _real        #  Charge density, potential array
bsav(0:nx,0:ny,0:nz)      _real        #  "Save" array for b
attx(0:nx-1)              _real        #  Attenuation factor as fcn. of kx
atty(0:ny-1)              _real        #  Attenuation factor as fcn. of ky
attz(0:nz)                _real        #  Attenuation factor as fcn. of kz
kxsq(0:nx-1)              _real        #  Discrete analog to kx^2/4Pi
kysq(0:ny-1)              _real        #  Discrete analog to ky^2/4Pi
kzsq(0:nz)                _real        #  Discrete analog to kz^2/4Pi
work(nx+ny-4)             _real        #  Workspace for sine transforms
work2d(1:nx+2,1:ny+2)     _real        #  Workspace for 2d fieldsolver
    #  Note: should be ((nx+2)*(ny+2)-7), but MAC produces bad code !
zwork(2,0:nx,0:nz)        _real        #  Workspace to optimize vpftz

******* Transpose_work_space:
phi_trnspssize /0/       integer
phi_trnsps(phi_trnspssize)  _real

******** CapMat3d dump:
nc3dmax         integer /0/ # Maximum number of points in conductors
nc3d            integer /0/ # Number of points within conductors
nc3dz           integer /0/ # Number of z grid points (generally = w3d.nz)
nc3dz2          integer /0/ # Number of matrices minus 1 (for capacity
                            # matrices in kz space)
                            # For serial, nc3dz2 = nz/2 (takes advantage
                            # of equality of matrices for slices iz and nz-iz).
                            # For parallel, nc3dz2 = nz-1.
xcond3d(nc3dmax)  _real [m] # X coordinate of points in conductors
ycond3d(nc3dmax)  _real [m] # Y coordinate of points in conductors
zcond3d(nc3dmax)  _real [m] # Z coordinate of points in conductors
vcond3d(nc3dmax)  _real [V] # voltage of points in conductors
pcond3d(nc3dmax,0:nc3dz)  _real [V] # actual potential on points in conductors
                                    # (auto set)
qcond3d(nc3dmax,0:nc3dz)  _real [C] # induced charge on points in conductors
                                    # (auto set)
cmat3d(nc3d,nc3d,0:nc3dz2) _real     # Capacity matrix (auto set)
kpvt3d(nc3d)      _integer  # Pivot points for matrix solve (auto set)

******** Pipe3d dump:
piperadi real    /0./    # Radius of pipe
pipenz   integer /0/     # Length of pipe grid
pipenz2  integer /0/     # Length of pipe grid (see comments on nc3dz2)
pipen    integer /0/     # Number of points on conductor in each kz slice
pipe8th  integer /0/     # Number of points in each eighth plus one
pipex(pipen+1)                  _real # X coordinates of points
pipey(pipen+1)                  _real # Y coordinates of points
pipephi(pipen,0:pipenz-1)       _real # Potential at each point
pipeq(pipen,0:pipenz-1)         _real # induced charge at each point
cap3d(pipen,pipe8th,0:pipenz2)  _real # Capacity matrix for kz slices
       # Reduced in size by using eight fold symmetry of circle in square
kpvt(pipen)                  _integer # Pivot points for matrix solve

******** Capvars dump:
qdslclen              real  # length of quad z slice
quadcent              real  # distance of quad centers from pipe center
quadradi              real  # radius of quad conductors
loadquad              logical /.true./ # should quad points be loaded?
nzquad                integer /0/  # number of z slices in conductor
nzpts                 integer /0/  # number of points per z slice in conductor
nendquad              integer /0/  # number of points in ends of conductor
ncndpts               integer /1/  # number of points in one conductor
numends               integer /0/  # number of starts of quads in grid
nendmax               integer /20/ # max number of starts of quads in grid
quadx(ncndpts,4)     _real  # x locations of four conductors
quady(ncndpts,4)     _real  # y locations of four conductors
quadz(ncndpts)       _real  # z locations of conductors
quadv(ncndpts,4)     _real  # relative sign of voltage
quadrho(ncndpts,4,nendmax) _real limited (ncndpts,4,numends)
                            # induced charges on conductors
quadcap(4*ncndpts,4*ncndpts) _real # capacity matrix, in form of reduced inverse
quadend(nendmax)       _real limited (numends)   # locations of quad starts
quadvlt(nendmax)       _real limited (numends)   # voltage of quad
vshift(nendmax)        _real limited (numends)   # shift in focusing voltage
kkkk(4*ncndpts)          _integer # Pivot points for matrix solve

*********** PSOR3d dump:
sorparam  real    /1.88/ # Overrelaxation parameter, default value is a
                         # reasonable choice, but it should be optimized for
                         # the particular configuration.
soriter   integer /0/    # Number of iterations needed to reach tolerance.
sormaxit  integer /2000/ # Maximum number of iterations allowed.
sortol    real   /1.e-6/ # Tolerance for convergence
sorerror  real           # Maximum error after convergence
isorerr   integer /5/    # Frequency of error calculation
lchebshv  logical /.false./ # Turns on Chebyshev acceleration
pjacobi   real    /.99/  # spectral radius of Jacobi iteration, used only for
                         # Chebyshev acceleration
gridmode    integer /0/ # Mode of grid motion, if 0, then normal, if 1, then
                        # grid is assumed not to move.  slice arrays only set
                        # when iwhich is set to -2.
zparity   integer /0/ +parallel # iz parity, used in the parallel version so that the
                       # parity of the subgrid conductor points can be
                       # relative to the full grid.
dxfine    real         # Size of transverse grid cells are finest level
lplates  logical /.true./ # Sets whether or not quadruple endplates are included
rodfract        real /1./ # Fraction of quadrupole rod which is used

%%%%%%%%%% ConductorInteriorType:
nmax           integer /0/ # Maximum number of points in conductor
n              integer /0/ # Number of points within conductors
indx(0:2,nmax)  _integer # Coordinates of points in conductor
volt(nmax)    _real    # Voltage of points in conductor
numb(nmax)    _integer # Number of the conductor the points are in
ilevel(nmax)  _integer /-1/ # Coarseness level at which the point is on grid
istart(0:100)  integer /1/ # Start of the conductor points for each MG level

%%%%%%%%%% ConductorSubGridType:
nmax            integer /0/ # Maximum number of points for sub-grid boundaries
n               integer /0/ # Number of points for sub-grid boundaries
prevphi(nmax)  _real    # Saves phi for sub-grid boundaries
indx(0:2,nmax)   _integer # Location of points for sub-grid boundaries
dels(0:5,nmax)   _real    # Distances to the surface - stored in the order
                          # mx, px, my, py, mz, pz
volt(0:5,nmax) _real    # Voltage of points for sub-grid boundaries
numb(0:5,nmax) _integer # ID of the conductor the points are in
ilevel(nmax)   _integer /-1/ # Coarseness level at which the point is on grid
istart(0:100)   integer /1/ # Start of the conductor data for each MG level

%%%%%%%%%% ConductorType:
interior ConductorInteriorType   # Interior of the conductors
evensubgrid ConductorSubGridType # Even subgrid data for conductors
oddsubgrid ConductorSubGridType  # Odd subgrid data for conductors
levels          integer  # Number of coarsening levels
leveliz(0:100)  integer  # List of iz for the levels of coarsening
levelnz(0:100)  integer  # List of nz for the levels of coarsening
levellx(0:100)  real /1/ # List of coarsening factors in x
levelly(0:100)  real /1/ # List of coarsening factors in y
levellz(0:100)  real /1/ # List of coarsening factors in z
fuzzsign integer /-1/    # When -1, subgrid points with distances == 1 are
                         # skipped, when +1 not skipped.

*********** Conductor3d dump parallel:
conductors ConductorType # Default data structure for conductor data
lcndbndy logical /.true./ # Turns on sub-grid boundaries
icndbndy integer /1/      # Type of interpolant to use for sub-grid boundaries
                          # 1 egun style
                          # 2 EBC style (non-centered finite-difference)
laddconductor logical /.false./ # When true, the python function
                          # calladdconductor is called at the beginning of the 
                          # field solve.
checkconductors(nx:integer,ny:integer,nz:integer,nzfull:integer,
                dx:real,dy:real,dz:real,conductors:ConductorType,
                my_index:integer,nslaves:integer,
                izfsslave:integer,nzfsslave:integer) subroutine

*********** MGLevels3d:
mglevels              integer /0/  # Number of coarsening levels
mglevelsnx(0:100)     integer      # List of nx for the levels of coarsening
mglevelsny(0:100)     integer      # List of ny for the levels of coarsening
mglevelsnzfull(0:100) integer      # List of nzfull for the levels of coarsening
mglevelsiz(0:100)     integer      # List of iz for the levels of coarsening
mglevelsnz(0:100)     integer      # List of nz for the levels of coarsening
mglevelslx(0:100)     real /101*1/ # List of coarsening factors in x
mglevelsly(0:100)     real /101*1/ # List of coarsening factors in y
mglevelslz(0:100)     real /101*1/ # List of coarsening factors in z
mglevelspart(0:100)   logical      # List of flags for whether full or partial
                                   # coarsening is done: 0 is full, 1 is partial
*********** Multigrid3d dump:
mgparam     real    /1.2/ # Acceleration parameter for multigrid fieldsolver
mgmaxiters  integer /100/ # Maximum number of iterations
mgmaxlevels integer /101/ # Minimum grid size in x-y to coarsen to
mgiters     integer       # Actual number of iterations
mgtol       real  /1.e-6/ # Absolute tolerance in change in last iteration
mgerror     real          # Maximum error after convergence
mgform      integer /1/   # When 1, MG operates on phi (and rho),
                          # when 2, MG operates on error (and residual)
downpasses  integer /1/   # Number of downpasses
uppasses    integer /1/   # Number of uppasses
bounds(0:5) integer /6*0/ # Type of boundaries at edge of mesh, in order of
                          # lower, upper for x, y, z.
mggoodnumbers(56) integer /2,4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,
                           80,96,112,128,160,192,224,256,320,384,448,512,
                           640,768,896,1024,1280,1536,1792,2048,2560,3072,
                           3584,4096,5120,6144,7168,8192,10240,12288,14336,
                           16384,20480,24576,28672,32768,40960,49152,57344,
                           65536/
                         # A list of good numbers to use for the grid
                         # dimension. This is an ordered list of powers of two
                         # times 1, 3, 5, and 7.
lbuildquads logical /.true./ # When true, quad elements are constructed
                             # automatically if defined.
getmglevels(nx:integer,ny:integer,nz:integer,nzfull:integer,
            dx:real,dy:real,dz:real,conductors:ConductorType,
            my_index:integer,nslaves:integer,
            izfsslave:integer,nzfsslave:integer)
   subroutine
   # Calculates levels of coarsening. Note that mglevels
   # must be zero when calling this routine.
multigrid3df(iwhich:integer,nx:integer,ny:integer,nz:integer,nzfull:integer,
             dx:real,dy:real,dz:real,
             phi(0:nx,0:ny,-1:nz+1):real,rho(0:nx,0:ny,0:nz):real,
             rstar(-1:nz+1):real,linbend:logical,
             bound0:integer,boundnz:integer,boundxy:integer,
             l2symtry:logical,l4symtry:logical,
             xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real)
   subroutine
   # Solves Poisson's equation using the multigrid method. This uses variables
   # from the f3d package to control the iterations and conductors.
multigrid3dsolve(iwhich:integer,nx:integer,ny:integer,nz:integer,nzfull:integer,
                 dx:real,dy:real,dz:real,
                 phi(0:nx,0:ny,-1:nz+1):real,rho(0:nx,0:ny,0:nz):real,
                 rstar(-1:nz+1):real,linbend:logical,bounds:integer,
                 xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real,
                 mgparam:real,mgform:integer,mgiters:integer,
                 mgmaxiters:integer,mgmaxlevels:integer,mgerror:real,mgtol:real,
                 downpasses:integer,uppasses:integer,
                 lcndbndy:logical,laddconductor:logical,icndbndy:integer,
                 lbuildquads:logical,gridmode:integer,conductors:ConductorType,
                 my_index:integer,nslaves:integer,
                 izfsslave:integer,nzfsslave:integer)
   subroutine
   # Solves Poisson's equation using the multigrid method. All input is
   # through the argument list.
residual(nx:integer,ny:integer,nz:integer,nzfull:integer,
         dxsqi:real,dysqi:real,dzsqi:real,
         phi(0:nx,0:ny,-1:nz+1):real,rho(0:nx,0:ny,0:nz):real,
         res(-1:nx+1,-1:ny+1,-3:nz+3):real,
         mglevel:integer,bounds:integer,
         mgparam:real,mgform:integer,mgform2init:logical,
         lcndbndy:logical,icndbndy:integer,conductors:ConductorType)
   subroutine
   # Calculates the residual
restrict2d(nx:integer,ny:integer,nz:integer,res(-1:nx+1,-1:ny+1,-3:nz+3):real,
           nxcoarse:integer,nycoarse:integer,
           rhocoarse(0:nxcoarse,0:nycoarse,0:nz):real,
           ff:integer,bounds:integer)
   subroutine
   # Restricts phi in 2 transverse dimensions
expand2d(nx:integer,ny:integer,nz:integer,phi(0:nx,0:ny,-1:nz+1):real,
         nxcoarse:integer,nycoarse:integer,phicoarse(0:nx,0:ny,-1:nz+1):real,
         bounds(0:5):integer)
   subroutine
   # Expands phi in 2 transverse dimensiosn
restrict3d(nx:integer,ny:integer,nz:integer,nzfull:integer,
           res(-1:nx+1,-1:ny+1,-3:nz+3):real,
           nxcoarse:integer,nycoarse:integer,nzcoarse:integer,
           nzfullcoarse:integer,
           rhocoarse(0:nxcoarse,0:nycoarse,0:nzcoarse):real,
           ff:real,bounds(0:5):integer,boundscoarse(0:5):integer,
           lzoffset:integer)
   subroutine
   # Restricts phi in 3 dimensions
expand3d(nx:integer,ny:integer,nz:integer,nzfull:integer,
         phi(0:nx,0:ny,-1:nz+1):real,
         nxcoarse:integer,nycoarse:integer,nzcoarse:integer,
         nzfullcoarse:integer,
         phicoarse(0:nxcoarse,0:nycoarse,-1:nzcoarse+1):real,
         bounds(0:5):integer,lzoffset:integer)
   subroutine
   # Expands phi in 3 dimensiosn
sorhalfpass3d(parity:integer,mglevel:integer,
              nx:integer,ny:integer,nz:integer,nzfull:integer,
              phi(0:nx,0:ny,-1:nz+1):real,rho(0:nx,0:ny,0:nz):real,
              rstar(-1:nz+1):real,dxsqi:real,dysqi:real,dzsqi:real,
              linbend:logical,bendx((nx+1)*(ny+1)):real,bounds(0:5):integer,
              mgparam:real,mgform:integer,
              lcndbndy:logical,icndbndy:integer,conductors:ConductorType)
   subroutine
   # Performs one pass of SOR relaxation, either even or odd.
cond_potmg(interior:ConductorInteriorType,nx:integer,ny:integer,nz:integer,
           phi(0:nx,0:ny,-1:nz+1):real,
           mglevel:integer,mgform:integer,mgform2init:logical)
    subroutine
    # Sets voltage on interior of conductors
cond_zerorhointerior(interior:ConductorInteriorType,nx:integer,ny:integer,nz:integer,
                     rho(0:nx,0:ny,0:nz):real)
    subroutine
    # Sets rho to zero inside conductor points.
cond_sumrhointerior(interior:ConductorInteriorType,
                    nx:integer,ny:integer,nz:integer,rho(0:nx,0:ny,0:nz):real,
                    ixmin:integer,ixmax:integer,iymin:integer,iymax:integer,
                    izmin:integer,izmax:integer) real function
subcond_sumrhointerior(rhosum:real,interior:ConductorInteriorType,
                    nx:integer,ny:integer,nz:integer,rho(0:nx,0:ny,0:nz):real,
                    ixmin:integer,ixmax:integer,iymin:integer,iymax:integer,
                    izmin:integer,izmax:integer) subroutine

*********** MultigridBE3d dump:
multigridbe3df(iwhich:integer,nx:integer,ny:integer,nz:integer,nzfull:integer,
             dx:real,dy:real,dz:real,phi:real,rho:real,
             rstar:real,linbend:logical,
             bound0:integer,boundnz:integer,boundxy:integer,
             l2symtry:logical,l4symtry:logical,
             xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real,
             iondensity:real,electrontemperature:real,plasmapotential:real,
             electrondensitymaxscale:real)
   subroutine
   # Solves Poisson's equation using the multigrid method, including the
   # Boltzmann electron term. This uses variables
   # from the f3d package to control the iterations and conductors.
multigridbe3dsolve(iwhich:integer,
             nx:integer,ny:integer,nz:integer,nzfull:integer,
             dx:real,dy:real,dz:real,phi:real,rho:real,
             rstar:real,linbend:logical,bounds:integer,
             xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real,
             mgparam:real,mgiters:integer,mgmaxiters:integer,
             mgmaxlevels:integer,mgerror:real,mgtol:real,
             downpasses:integer,uppasses:integer,
             lcndbndy:logical,laddconductor:logical,icndbndy:integer,
             lbuildquads:logical,gridmode:integer,conductors:ConductorType,
             my_index:integer,nslaves:integer,
             izfsslave:integer,nzfsslave:integer,
             iondensity:real,electrontemperature:real,plasmapotential:real,
             electrondensitymaxscale:real)
   subroutine
   # Solves Poisson's equation using the multigrid method, including the
   # Boltzmann electron term. All input is through the argument list.
lphibe(nx:integer,ny:integer,nz:integer,nzfull:integer,dxsqi:real,dysqi:real,dzsqi:real,phi:real,res:real,
       mglevel:integer,bounds:integer,mgparam:real,
       lcndbndy:logical,icndbndy:integer,conductors:ConductorType,
       iondensity:real,electrontemperature:real,plasmapotential:real,
       electrondensitymaxscale:real) subroutine
restrictbe3d(nx:integer,ny:integer,nz:integer,nzfull:integer,u:real,delz:integer,
             nxcoarse:integer,nycoarse:integer,nzcoarse:integer,nzfullcoarse:integer,ucoarse:real,
             bounds:integer,boundscoarse:integer,lzoffset:integer,
             my_index:integer,nslaves:integer,izfsslavec:integer,nzfsslavec:integer,
             whosendingleftc:integer,izsendingleftc:integer,
             whosendingrightc:integer,izsendingrightc:integer) subroutine
expandbe3d(nx:integer,ny:integer,nz:integer,nzfull:integer,phi:real,
           nxcoarse:integer,nycoarse:integer,nzcoarse:integer,nzfullcoarse:integer,phicoarse:real,
           bounds:integer,lzoffset:integer) subroutine
relaxbe3d(mglevel:integer,nx:integer,ny:integer,nz:integer,nzfull:integer,phi:real,rho:real,rstar:real,
          dxsqi:real,dysqi:real,dzsqi:real,linbend:logical,bendx:real,bounds:integer,
          mgparam:real,lcndbndy:logical,icndbndy:integer,conductors:ConductorType,
          my_index:integer,nslaves:integer,izfsslave:integer,nzfsslave:integer,
          whosendingleft:integer,izsendingleft:integer,
          whosendingright:integer,izsendingright:integer,
          iondensity:real,electrontemperature:real,plasmapotential:real,
          electrondensitymaxscale:real) subroutine
cond_potmgbe(interior:ConductorInteriorType,nx:integer,ny:integer,nz:integer,phi:real,mglevel:integer)
          subroutine
cond_potmgbezero(interior:ConductorInteriorType,nx:integer,ny:integer,nz:integer,u:real,mglevel:integer,
                 delt:integer,delz:integer) subroutine
condbndymgbe(subgrid:ConductorSubGridType,nx:integer,ny:integer,nz:integer,phi:real,rho:real,dxsqi:real,dysqi:real,dzsqi:real,
             mgparam:real,bounds:integer,mglevel:integer,icndbndy:integer,
             iondensity:real,electrontemperature:real,plasmapotential:real,
             electrondensitymaxscale:real) subroutine
copyphiwithguardcells(nx:integer,ny:integer,nz:integer,nzfull:integer,phiin:real,phiout:real,bounds:integer)
subroutine
copyrhowithguardcells(nx:integer,ny:integer,nz:integer,nzfull:integer,rhoin:real,rhoout:real,bounds:integer)
subroutine
vcyclebe(mglevel:integer,nx:integer,ny:integer,nz:integer,nzfull:integer,dx:real,dy:real,dz:real,
         phi:real,rho:real,rstar:real,linbend:logical,bendx:real,bounds:integer,
         mgparam:real,
         mgmaxlevels:integer,downpasses:integer,uppasses:integer,
         lcndbndy:logical,icndbndy:integer,conductors:ConductorType,
         my_index:integer,nslaves:integer,izfsslave:integer,nzfsslave:integer,
         iondensity:real,electrontemperature:real,plasmapotential:real,
         electrondensitymaxscale:real) subroutine

*********** Multigrid3d_work:
# Temporary variables and array used by subgrid_sor_to_mg
wnx integer
wny integer
wnz integer
iii(0:wnx,0:wny,0:wnz) _integer


%%%%%%%%%%% BFieldGridType:
nx     integer  /0/  # Number of grid cells in x in B grid
ny     integer  /0/  # Number of grid cells in y in B grid
nz     integer  /0/  # Number of grid cells in z in B grid
nzfull integer  /0/  # Number of grid cells in z in B grid
dx     real [m] /0./ # x grid cell size in B grid
dy     real [m] /0./ # y grid cell size in B grid
dz     real [m] /0./ # z grid cell size in B grid
xmmin  real [m] /0./ # X lower limit of mesh
xmmax  real [m] /0./ # X upper limit of mesh
ymmin  real [m] /0./ # Y lower limit of mesh
ymmax  real [m] /0./ # Y upper limit of mesh
zmmin  real [m] /0./ # Z lower limit of mesh
zmmax  real [m] /0./ # Z upper limit of mesh
zmminglobal real [m] # Global value of zmmin
zmmaxglobal real [m] # Global value of zmmax

lcndbndy logical /.true./ # Turns on sub-grid boundaries
icndbndy integer /1/      # Type of interpolant to use for sub-grid boundaries
                          # 1 egun style
                          # 2 EBC style (non-centered finite-difference)
laddconductor logical /.false./ # When true, the python function
                          # calladdconductor is called at the beginning of the
                          # field solve.
mgparam(0:2)     real    /1.2/ # Acceleration parameter for multigrid solver
mgmaxiters(0:2)  integer /100/ # Maximum number of multigrid iterations
mgmaxlevels(0:2) integer /101/ # Minimum grid size in x-y to coarsen to
mgiters(0:2)     integer       # Actual number of multigrid iterations
mgtol(0:2)       real    /0./  # Absolute tolerance in change in last iteration
mgerror(0:2)     real          # Maximum error after convergence
mgform(0:2)      integer /1/   # When 1, MG operates on phi (and rho),
                          # when 2, MG operates on error (and residual)
downpasses(0:2)  integer /1/   # Number of downpasses
uppasses(0:2)    integer /1/   # Number of uppasses
bounds(0:5) integer  # Boundary conditions on grid surfaces



lcylindrical logical /.false./ # When true, signifies that cylindrical
                               # coordinates are being used, which means
                               # that 0 is r, 1 is theta, and 2 is z.
lusevectorpotential logical /.true./ # When true, the vector potential A is
                                      # solver for, solving del^2 A = -mu0 J.
                                      # Otherwise solve del^2 B = -mu0 curl J.
lanalyticbtheta logical /.false./ # When true, Btheta is calculated from Jz
                                  # using Btheta = mu0*Iz(r)/(2*pi*r) where
                                  # Iz(r) is the total current inside the
                                  # radius r.
                                  # Warning: this option seems to give
                                  # unstable results.

j(0:2,0:nx,0:ny,0:nz) _real # Current density
b(0:2,0:nx,0:ny,0:nz) _real # B field, calculated from B = del cross A
a(0:2,-1:nx+1,-1:ny+1,-1:nz+1) _real
  # Vector magnetic potential, calculated from del sq A = J

attx(0:nx-1)     _real           # Attenuation factor as fcn. of kx
atty(0:ny-1)     _real           # Attenuation factor as fcn. of ky
attz(0:nzfull)   _real           # Attenuation factor as fcn. of kz
kxsq(0:nx-1)     _real [1/m**2]  # Discrete analog to kx^2/4Pi
kysq(0:ny-1)     _real [1/m**2]  # Discrete analog to ky^2/4Pi
kzsq(0:nzfull)   _real [1/m**2]  # Discrete analog to kz^2/4Pi

rstar(-1:nz+1)         _real [m] # Radius of curv of reference orbit
scrtch(2*nx+2*ny)      _real     # Scratch for fieldsolve
xywork(2,0:nx,0:ny)    _real     # Work space for transverse FFTs
zwork(2,0:nx,0:nzfull) _real     # Work space used to optimize vsftz

nsjtmp integer /0/
jsjtmp(0:nsjtmp-1) _integer /-1/ #
nsndtsj integer /0/
jtmp(3,0:nx,0:ny,0:nz,0:nsndtsj-1)  _real
             # Temporary copy of the current density from the particles.

conductors ConductorType # Default data structure for conductor data

*********** BFieldGrid:
bfield BFieldGridType
bfieldp BFieldGridType
init_bfieldsolver(bfstype:integer) subroutine # Initializes the B-field solver
bfieldsol3d(iwhich) subroutine # Self B-field solver
loadj3d(pgroup:ParticleGroup,ins:integer,nps:integer,is:integer,lzero:logical)
             subroutine # Provides a simple interface to the current density
                        # loading routine setj3d
setaboundaries3d(bfield:BFieldGridType)
             subroutine #
perj3d(j:real,nx:integer,ny:integer,nz:integer,bound0:integer,boundxy:integer)
             subroutine #
setb3d(bfield:BFieldGridType,np:integer,xp:real,yp:real,zp:real,zgrid:real,
       bx:real,by:real,bz:real,l2symtry:logical,l4symtry:logical)
             subroutine #
fetchafrompositions3d(np:integer,xp:real,yp:real,zp:real,a:real,zgrid:real,
                      bfield:BFieldGridType,l2symtry:logical,l4symtry:logical)
             subroutine #
getbfroma3d(bfield:BFieldGridType)
             subroutine #
curl3d(a:real,b:real,nx:integer,ny:integer,nz:integer,dx:real,dy:real,dz:real,
       xmmin:real,lcylindrical:logical,
       adelx:integer,adelz:integer,bdelx:integer,bdelz:integer)
             subroutine # Calculates the curl of the input array a, putting
                        # the result into b.
setj3d(bfield:BFieldGridType,j1d:real,np:integer,xp:real,yp:real,zp:real,
       zgrid:real,uxp:real,uyp:real,uzp:real,gaminv:real,q:real,
       nw:integer,wght:real,depos:string,l2symtry:logical,l4symtry:logical)
             subroutine # Computes current density
getjforfieldsolve()
             subroutine #
getjforfieldsolve3d(bfield:BFieldGridType,bfieldp:BFieldGridType,
                    my_index:integer,nslaves:integer,izfsslave:integer,
                    nzfsslave:integer,izpslave:integer,nzpslave:integer)
             subroutine #
setupbfieldsforparticles3d(ns:integer,ndts:integer,it:integer,
                           bfield:BFieldGridType,bfieldp:BFieldGridType)
             subroutine #
fetchb3dfrompositions(is:integer,n:integer,x(n):real,y(n):real,z(n):real,
                      bx(n):real,by(n):real,bz(n):real)
             subroutine #
fetcha(n:integer,x(n):real,y(n):real,z(n):real,a(n):real)
             subroutine #
getbforparticles()
             subroutine #
bvp3d(iwhich:integer,bfstype:integer)
             subroutine #

*********** AMR3droutines:
gatherrhofromchild(rho:real,nn:integer,childrho:real,cnn:integer,
                   l:integer,u:integer,fulllower:integer,
                   childlower:integer,childupper:integer,
                   r:integer,weights:real,
                   dobounds:integer,bounds:integer,rootdims:integer)
      subroutine
gatherphifromparents(phi:real,nn:integer,l:integer,u:integer,fulllower:integer,
                     parentphi:real,pnn:integer,parentlower:integer,r:integer)
      subroutine
gatherafromparents(a:real,nn:integer,l:integer,u:integer,fulllower:integer,
                     parenta:real,pnn:integer,parentlower:integer,r:integer)
      subroutine
gatherjfromchild(j:real,nn:integer,childj:real,cnn:integer,
                   l:integer,u:integer,fulllower:integer,
                   childlower:integer,childupper:integer,
                   r:integer,weights:real,
                   lcylinderical:logical,radius:real,cradius:real,
                   dobounds:integer,bounds:integer,rootdims:integer)
      subroutine

*********** Surface_of_Rev dump:
srfrv_pernz              integer /0/ # Number of points per nz for tablized data
srfrv_z                  real # Value of z passed to srfrv_f
srfrv_r                  real # Value of r returned by srfrv_f
lsrlinr                  logical /.false./ # Use piecewise-linear curve
npnts_sr                 integer /0/ # Number points in piecewise-linear curve
z_sr(npnts_sr)           _real # Z of points in piecewise-linear curve
r_sr(npnts_sr)           _real # R of points in piecewise-linear curve
rad_sr(npnts_sr-1)       _real /LARGEPOS/ # Radius of curvature of curve arc
zc_sr(npnts_sr-1)        _real # Z center of circle
rc_sr(npnts_sr-1)        _real # R center of circle
lsrminlinr               logical /.false./ # Use piecewise-linear curve for rmin
npnts_srmin              integer /0/ # Number points in piecewise-linear curve
z_srmin(npnts_srmin)     _real # Z of points in piecewise-linear curve
r_srmin(npnts_srmin)     _real # R of points in piecewise-linear curve
rad_srmin(npnts_srmin-1) _real /LARGEPOS/ # Radius of curvature of curve arc
zc_srmin(npnts_srmin-1)  _real # Z center of circle
rc_srmin(npnts_srmin-1)  _real # R center of circle
lsrmaxlinr               logical /.false./ # Use piecewise-linear curve for rmax
npnts_srmax              integer /0/ # Number points in piecewise-linear curve
z_srmax(npnts_srmax)     _real # Z of points in piecewise-linear curve
r_srmax(npnts_srmax)     _real # R of points in piecewise-linear curve
rad_srmax(npnts_srmax-1) _real /LARGEPOS/ # Radius of curvature of curve arc
zc_srmax(npnts_srmax-1)  _real # Z center of circle
rc_srmax(npnts_srmax-1)  _real # R center of circle
lsrfindrextremum         logical /.false./ # When true, an extra search is
  # done to find an r extremum that lies between grid points. Set this to true
  # if there are structures taller than dr but shorter than dz.

*********** LantzSolverTemp:
nxlan integer
nylan integer
nzlan integer
nzfulllan integer
nxtranlan integer
nytranlan integer
nztranlan integer
alan(0:nxlan,0:nylan,0:nzlan) _real
blan(0:nxlan,0:nylan,0:nzlan) _real
clan(0:nxlan,0:nylan,0:nzlan) _real
atranlan(0:nxtranlan,0:nytranlan-1,0:nzfulllan) _real
btranlan(0:nxtranlan,0:nytranlan-1,0:nzfulllan) _real
ctranlan(0:nxtranlan,0:nytranlan-1,0:nzfulllan) _real
dtranlan(0:nxlan,0:nylan,0:2) _real

*********** PSOR3d_subs:
psor3df(iwhich,nx,ny,nz,phi:real,rho:real,phi1d:real,rho1d:real,rstar:real,
       dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,zbeam:real,
       zgrid:real,linbends:logical,bound0:integer,boundnz:integer,
       boundxy:integer,l2symtry:logical,l4symtry:logical,lzerophiedge:logical,
       scrtch:real,izfsmin:integer,izfsmax:integer)
     subroutine # PSOR field solver
psorinit(nx,ny,nz,dx:real,dy:real,dz:real,l2symtry:logical,l4symtry:logical)
     subroutine # Initialize arrays that hold conductor points
cond_pot(nx,ny,nz,phi:real)
     subroutine # Sets potential in phi to desired potential on conductors
setcndtr(xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real,nx,ny,nz,
         dx:real,dy:real,dz:real,bound0:integer,boundnz:integer,boundxy:integer,
         l2symtry:logical,l4symtry:logical)
     subroutine # Calculates conductor locations
srfrvoutf(rofzfunc:string,volt:real,zmin:real,zmax:real,
         xcent:real,ycent:real,rmax:real,lfill:logical,
         xmin:real,xmax:real,ymin:real,ymax:real,lshell:logical,
         zmmin:real,zmmax:real,zbeam:real,
         dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,
         ix_axis:integer,iy_axis:integer,xmesh:real,ymesh:real,
         l2symtry:logical,l4symtry:logical,condid:integer)
     subroutine # Set conductor points for a conductor that is the
                # outside of a surface of revolution.
srfrvinf(rofzfunc:string,volt:real,zmin:real,zmax:real,
         xcent:real,ycent:real,rmin:real,lfill:logical,
         xmin:real,xmax:real,ymin:real,ymax:real,lshell:logical,
         zmmin:real,zmmax:real,zbeam:real,
         dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,
         ix_axis:integer,iy_axis:integer,xmesh:real,ymesh:real,
         l2symtry:logical,l4symtry:logical,condid:integer)
     subroutine # Set conductor points for a conductor that is the
                # inside of a surface of revolution.
srfrvinoutf(rminofz:string,rmaxofz:string,volt:real,zmin:real,zmax:real,
         xcent:real,ycent:real,lzend:logical,
         xmin:real,xmax:real,ymin:real,ymax:real,lshell:logical,
         zmmin:real,zmmax:real,zbeam:real,
         dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,
         ix_axis:integer,iy_axis:integer,xmesh:real,ymesh:real,
         l2symtry:logical,l4symtry:logical,condid:integer)
     subroutine # Set conductor points for a conductor that is between two
                # surfaces of revolution.
srfrv_f(zz:real,rofzfunc:string,icase:integer,izflag:integer)
     real function # Function used to calculate r of z given some input.

*********** F3Dsubs:
#  Callable subroutines in the F3D package
vpois3d  (iwhich:integer,a:real,ak:real,kxsq:real,kysq:real,kzsq:real,
         attx:real,atty:real,attz:real,filt:real,
         lx:real,ly:real,lz:real,nx:integer,ny:integer,nz:integer,
         nzfull:integer,
         w:real,xywork:real,zwork:real,ibc:integer,
         l2symtry:logical,l4symtry:logical,
         bound0:integer,boundnz:integer,boundxy:integer)
     subroutine #  The 3d Poisson solver
vp3x(iwhich) subroutine
     # Python-level interface to VPOIS3d, using FS3 database variables
     #  The user program should declare a similar subroutine w/ its vars.
vpois2d  (iwhich, a:real, ak:real, kxsq:real, kysq:real, attx:real, atty:real,
          filt:real, lx:real, ly:real, nx, ny, work:real, xywork:real, ibc,
          l2symtry:logical,l4symtry:logical)
     subroutine #  The 2d Sine-Sine Poisson solver
vsftx    (a:real, work:real, cp:real, cm:real, nx, ny, isetup)
     subroutine #  Vectorized Sine Fourier Transform in X
vsfty    (a:real, work:real, cp:real, cm:real, nx, ny, isetup)
     subroutine #  Vectorized Sine Fourier Transform in Y
vpftx    (nx,ny,nz,a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in X
vpfty    (nx,ny,nz,a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Y
vpftz    (nx,ny,nz,a:real,norm:real,esx,esy,nxy,zwork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Z
vpftxi    (nx,ny,nz,a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in X
vpftyi    (nx,ny,nz,a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Y
vpftzi    (nx,ny,nz,a:real,norm:real,esx,esy,nxy,zwork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Z
cosqx(a:real,w:real,c:real,nx:integer,ny:integer,isign:integer) subroutine
cosqy(a:real,w:real,c:real,nx:integer,ny:integer,isign:integer) subroutine
vcpft(r:real,i:real,n:integer,incp:integer,signp:integer,lenv:integer,
      lfd:integer) subroutine
attenuate(nx:integer,ny:integer,nz:integer,a:real,
          attx:real,atty:real,attz:real,
          ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
unattenuate(nx:integer,ny:integer,nz:integer,a:real,
            attx:real,atty:real,attz:real,
            ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
rhotophi(nx:integer,ny:integer,nz:integer,a:real,kxsq:real,kysq:real,kzsq:real,
         ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
phitorho(nx:integer,ny:integer,nz:integer,a:real,kxsq:real,kysq:real,kzsq:real,
         ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
pipe3df(iwhich, pipeshpe:real, rho:real, phi:real, kxsq:real, kysq:real,
        kzsq:real, attx:real, atty:real, attz:real, filt:real,
        xlen:real, ylen:real, zlen:real, nx:integer, ny:integer, nz:integer,
        nzfull:integer, scrtch:real,
        l2symtry:logical,l4symtry:logical)
     subroutine #  External interface to capacity field solver
pipest3d(pipeshpe:real, cap:real, phi:real, kxsq:real, kysq:real, kzsq:real,
         attx:real, atty:real, attz:real, filt:real, xlen:real,
         ylen:real, zlen:real,
         nx:integer,ny:integer,nz:integer,nzfull:integer,
         scrtch:real, piperadi:real,
         pipen:real,pipe8th:real,pipex:real,pipey:real,cap3d:real,kpvt,
         l2symtry:logical,l4symtry:logical)
     subroutine #  Sets up capacity matrix
setcap3d(ncndpts,quadx:real,quady:real,quadz:real,quadcap:real,kpvt:integer,
         nx:integer,ny:integer,nz:integer,nzfull:integer,
         phi:real,work:real,kxsq:real,kysq:real,kzsq:real,attx:real,atty:real,
         attz:real,xlen:real,ylen:real,zlen:real,filt:real,
         l2symtry:logical,l4symtry:logical)
     subroutine # Sets up 3-d cap matrix
capfs(iwhich,nx:integer,ny:integer,nz:integer,nzfull:integer,
      phi:real,rho:real,xlen:real,ylen:real,zlen:real,
      kxsq:real,kysq:real,kzsq:real,attx:real,atty:real,attz:real,
      filt:real,work:real,workz:real,xmmax:real,
      pipeshpe:string,periinz:logical,l2symtry:logical,l4symtry:logical)
     subroutine # 3-d field solve using cap matrices
cndctr3d(nx,ny,xmmax:real,dx:real,dy:real,quadradi:real,quadcent:real,ncndpts,
         quadx:real,quady:real,quadz:real,quadv:real,nzpts,nendquad,nzquad,
         qdslclen:real,pipeshpe:string,loadquad:logical)
     subroutine # Finds points on quadrupole conductor surfaces
findqdnd(nquad:integer,quadzs:real,quadde:real,quadvx:real,quadvy:real,
         zmmin:real,zgrid:real,nz:integer,dz:real,
         numends:integer,nendmax:integer,quadend:real,quadvlt:real,
         vshift:real,quadcent:real,quadradi:real)
     subroutine # Finds starts of quads that are on main grid
vcap3d(iwhich,rho:real,phi:real,kxsq:real,kysq:real,kzsq:real,attx:real,
       atty:real,attz:real,filt:real,xlen:real,ylen:real,zlen:real,
       nx:integer,ny:integer,nz:integer,nzfull:integer,
       scrtch:real,xmmax:real,zmmin:real,zgrid:real,
       pipeshpe:string,periinz:logical,l2symtry:logical,l4symtry:logical)
     subroutine # External routine for capacity matrix field solve

******** ConductorGeometryGenerators:
solvequartic(a0:real,a1:real,a2:real,a3:real,x1:complex,x2:complex,x3:complex,x4:complex) subroutine
setconductorparity(nn:integer,ix:integer,iy:integer,iz:integer,
                   dels:real,parity:integer,fuzz:real,fuzzsign:integer,
                   dfill:real) subroutine
zplaneconductorf(zcent:real,zsign:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zplaneconductord(zcent:real,zsign:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zplaneintercept(zcent:real,zsign:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
planeconductorf(z0:real,zsign:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
planeconductord(z0:real,zsign:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
planeintercept(z0:real,zsign:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
boxconductorf(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
boxconductord(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
boxintercept(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
cylinderconductorf(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
cylinderconductord(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
cylinderintercept(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
cylindersconductorf(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
cylindersconductord(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
cylindersintercept(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zcylinderconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zcylinderconductord(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zcylinderintercept(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zroundedcylinderconductorf(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zroundedcylinderconductord(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zroundedcylinderintercept(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zcylinderoutconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zcylinderoutconductord(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zcylinderoutintercept(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zroundedcylinderoutconductorf(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zroundedcylinderoutconductord(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zroundedcylinderoutintercept(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
sphereconductorf(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
sphereconductord(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
sphereintercept(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
coneconductorf(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
coneconductord(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
coneintercept(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
conesconductorf(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
conesconductord(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
conesintercept(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
ztorusconductorf(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
ztorusconductord(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
ztorusintercept(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
beamletplateconductorf(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
beamletplateconductord(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
beamletplateintercept(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zsrfrvoutconductorf(rofzfunc:string,zmin:real,zmax:real,rmax:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zsrfrvoutconductord(rofzfunc:string,zmin:real,zmax:real,rmax:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zsrfrvoutintercept(rofzfunc:string,zmin:real,zmax:real,rmax:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zsrfrvinconductorf(rofzfunc:string,zmin:real,zmax:real,rmin:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zsrfrvinconductord(rofzfunc:string,zmin:real,zmax:real,rmin:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zsrfrvinintercept(rofzfunc:string,zmin:real,zmax:real,rmin:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zsrfrvinoutconductorf(rminofz:string,rmaxofz:string,zmin:real,zmax:real,
        griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zsrfrvinoutconductord(rminofz:string,rmaxofz:string,zmin:real,zmax:real,
        griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,distance:real) subroutine
zsrfrvinoutintercept(rminofz:string,rmaxofz:string,zmin:real,zmax:real,
        griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine


******** ConductorGeometryVisualization:
maxtriangles integer /0/
ntriangles integer /0/
triangles(0:2,0:2,maxtriangles) _real
normals(0:2,0:2,maxtriangles) _real
connections(0:2,maxtriangles) _integer
maxpoints integer /0/
npoints integer
points(0:2,maxpoints) _real
pnormals(0:2,maxpoints) _real
getconductorfacets(nc:integer,icnd:integer,dels:real,
                   gridn:integer,griddd:real,gridmin:real) subroutine
getconductorsnewfacet(ix:integer,iy:integer,iz:integer,oo:integer,
                      parity:integer,gridn:integer,iii:integer,nc:integer,
                      dels:real,gridmin:real,griddd:real,pp:real,npp:integer)
                      subroutine
conductorsmoothshading() subroutine

******** Subtimersf3d:
lf3dtimesubs logical /.false./
timemultigrid3dsolve real /0./
timegatherrhofromchild real /0./
timegatherphifromparents real /0./
timegatherafromparents real /0./
timegatherjfromchild real /0./

timeexchange_phi           real /0./
timetranspose              real /0./
timetransposei             real /0./
timelantzsolver            real /0./
timegeneraltridiag         real /0./
timeparalleltridiag        real /0./
timeparallelgatherall      real /0./
timemgdividenz             real /0./
timemggetexchangepes       real /0./
timemgexchange_phi         real /0./
timemgexchange_phiperiodic real /0./
timemgexchange_rho         real /0./
timeprintarray3d           real /0./

timepera3d                     real /0./
timeperj3d                     real /0./
timesetb3d                     real /0./
timefetchafrompositions3d      real /0./
timegetbfroma3d                real /0./
timesetj3d                     real /0./
timeloadj3d                    real /0./
timefetchb3dfrompositions      real /0./
timebfieldsol3d                real /0./
timebvp3d                      real /0./
timesumjondomainboundaries     real /0./
timeperj3d_slave               real /0./
timegetjforfieldsolve3d        real /0./
timepera3d_slave               real /0./
timegetbforparticles3d         real /0./
timegetaforfields3d            real /0./

