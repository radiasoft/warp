f3d
#@(#) File F3D.V, version $Revision: 3.222 $, $Date: 2011/10/05 21:14:29 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package F3D of code WARP6
# fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
{
LARGEPOS = 1.0e+36 # This must be the same as in top.v
}

*********** F3Dversion:
versf3d character*19 /"$Revision: 3.222 $"/#  Code version version is set by CVS

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
a(-1:nx+1,-1:ny+1)        _real        #  2d charge density, potential array
b(-1:nx+1,-1:ny+1,-1:nz+1) _real        #  Charge density, potential array
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
phi_trnspsnx /0/ integer
phi_trnspsny /0/ integer
phi_trnspsnz /0/ integer
phi_trnspsnxguardphi /1/ integer
phi_trnspsnyguardphi /1/ integer
phi_trnspsnzguardphi /1/ integer
phi_trnsps(-phi_trnspsnxguardphi:phi_trnspsnx+phi_trnspsnxguardphi,-phi_trnspsnyguardphi:phi_trnspsny+phi_trnspsnyguardphi,-phi_trnspsnzguardphi:phi_trnspsnz+phi_trnspsnzguardphi) _real

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
qdslclen              real  /0./ # length of quad z slice
quadcent              real  /0./ # distance of quad centers from pipe center
quadradi              real  /0./ # radius of quad conductors
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

%%%%%%%%%% ConductorInteriorType:
nmax           integer  /0/  # Maximum number of points in conductor
n              integer  /0/  # Number of points within conductors
indx(0:2,nmax) _integer      # Coordinates of points in conductor
volt(nmax)     _real         # Voltage of points in conductor
numb(nmax)     _integer      # Number of the conductor the points are in
ilevel(nmax)   _integer /-1/ # Coarseness level at which the point is on grid
istart(0:100)  integer  /1/  # Start of the conductor points for each MG level

%%%%%%%%%% ConductorSubGridType:
nmax             integer  /0/  # Maximum number of points for sub-grid boundaries
n                integer  /0/  # Number of points for sub-grid boundaries
prevphi(nmax)    _real         # Saves phi for sub-grid boundaries
indx(0:2,nmax)   _integer      # Location of points for sub-grid boundaries
dels(0:5,nmax)   _real   /2./  # Distances to the surface - stored in the order
                               # mx, px, my, py, mz, pz
volt(0:5,nmax)   _real         # Voltage of points for sub-grid boundaries
numb(0:5,nmax)   _integer      # ID of the conductor the points are in
efield0(0:2,nmax) _real        # E-field at the point
efieldd(0:5,nmax) _real        # E-field along each of the six axis directions
ilevel(nmax)     _integer /-1/ # Coarseness level at which the point is on grid
istart(0:100)    integer  /1/  # Start of the conductor data for each MG level

%%%%%%%%%% MGCoefficients:
nx integer /0/
ny integer /0/
nz integer /0/
data(0:8,0:nx,0:ny,0:nz) _real
coarser _MGCoefficients
finer   _MGCoefficients

%%%%%%%%%% ConductorType:
interior ConductorInteriorType   # Interior of the conductors
evensubgrid ConductorSubGridType # Even subgrid data for conductors
oddsubgrid ConductorSubGridType  # Odd subgrid data for conductors
levels          integer  # Number of coarsening levels
levelix(0:100)  integer  # List of ix for the levels of coarsening
leveliy(0:100)  integer  # List of iy for the levels of coarsening
leveliz(0:100)  integer  # List of iz for the levels of coarsening
levelnx(0:100)  integer  # List of nx for the levels of coarsening
levelny(0:100)  integer  # List of ny for the levels of coarsening
levelnz(0:100)  integer  # List of nz for the levels of coarsening
levellx(0:100)  real /1/ # List of coarsening factors in x
levelly(0:100)  real /1/ # List of coarsening factors in y
levellz(0:100)  real /1/ # List of coarsening factors in z
fuzzsign integer /-1/    # When -1, subgrid points with distances == 1 are
                         # skipped, when +1 not skipped.
lcorrectede logical /.false./ # When true, the E field near conductors is
                              # calculated using a reduced finite
                              # difference. This only works when using the
                              # MultiGrid3D field solver that is created and
                              # registered in python. It will not work with
                              # fstype=7 or any of the FFT solvers.
nxct integer /0/
nyct integer /0/
nzct integer /0/
icgrid(0:nxct,0:nyct,0:nzct) _integer # Used to determine which conductor
                         # point is at each grid point when lcorrectede is on.
coeffs         _MGCoefficients # Precalculated coefficients, used when
                               # lprecalccoeffs is true.

*********** Conductor3d dump parallel:
conductors ConductorType # Default data structure for conductor data
lcndbndy logical /.true./ # Turns on sub-grid boundaries
icndbndy integer /2/      # Type of interpolant to use for sub-grid boundaries
                          # 1 egun style
                          # 2 EBC style (non-centered finite-difference)
lcorrectede logical /.false./ # When true, the E field near conductors is
                              # calculated using a reduced finite
                              # difference. This only works when using the
                              # MultiGrid3D field solver that is created and
                              # registered in python. It will not work with
                              # fstype=7 or any of the FFT solvers.
laddconductor logical /.false./ # When true, the python function
                          # calladdconductor is called at the beginning of the 
                          # field solve.
lprecalccoeffs logical /.false./ # When true, the finite difference coefficients
                                 # are precalculated and saved on a mesh. This
                                 # is faster but uses more memory.
checkconductors(nx:integer,ny:integer,nz:integer,
                nxlocal:integer,nylocal:integer,nzlocal:integer,
                dx:real,dy:real,dz:real,
                conductors:ConductorType,fsdecomp:Decomposition) subroutine

*********** MGLevels3d:
mglevels              integer /0/  # Number of coarsening levels
mglevelsnx(0:100)     integer      # List of nx for the levels of coarsening
mglevelsny(0:100)     integer      # List of ny for the levels of coarsening
mglevelsnz(0:100) integer          # List of nz for the levels of coarsening
mglevelsix(0:100)     integer      # List of ix for the levels of coarsening
mglevelsiy(0:100)     integer      # List of iy for the levels of coarsening
mglevelsiz(0:100)     integer      # List of iz for the levels of coarsening
mglevelslx(0:100)     real /101*1/ # List of coarsening factors in x
mglevelsly(0:100)     real /101*1/ # List of coarsening factors in y
mglevelslz(0:100)     real /101*1/ # List of coarsening factors in z
mglevelspart(0:100)   logical      # List of flags for whether full or partial
                                   # coarsening is done: 0 is full, 1 is partial
*********** Multigrid3d_diagnostic dump:
lprintmgphimaxchange logical /.false./ # When true, prints the maximum change in phi
                                       # at each level of the multigrid field solver.
lprintmgarraysumdiagnostic logical /.false./

*********** Multigrid3d dump:
gridmode    integer /0/ # Mode of grid motion. In most cases this will be
                        # autoset. Use with care.
                        # if 0, then grid is assumed to be moving, and the
                        # conductor information is regenerated every time step.
                        # This is an expensive operation and should be turned
                        # off is not needed.
                        # if 1, then grid is assumed not to move. The
                        # conductor information is only generated once.
mgparam     real    /1.2/ # Acceleration parameter for multigrid fieldsolver.
                          # This can be automatically set by calling
                          # find_mgparam(). It must be greater than zero and
                          # less than two.
mgmaxiters  integer /100/ # Maximum number of V-cycle iterations
mgmaxlevels integer /101/ # Maximum numer of levels of coarsening
mgiters     integer       # Actual number of iterations
mgtol       real  /1.e-6/ # Absolute tolerance in change in last iteration
                          # Note that this has units of volts.
mgerror     real          # Maximum error after convergence
mgform      integer /1/   # Do not use - not fully supported
                          # When 1, MG operates on phi (and rho),
                          # when 2, MG operates on error (and residual)
downpasses  integer /1/   # Number of downpasses, autoset by find_mgparam()
uppasses    integer /1/   # Number of uppasses, autoset by find_mgparam()
mgverbose   integer /1/   # Level of verbosity of multigrid solver
mgntverbose integer /1/   # Time step period when convergence information
                          # is printed
bounds(0:5) integer /6*0/ # Type of boundaries at edge of mesh, in order of
                          # lower, upper for x, y, z. This is generally
                          # autoset, but can be changed for special purposes.
mggoodnumbers(56) integer /2,4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,
                           80,96,112,128,160,192,224,256,320,384,448,512,
                           640,768,896,1024,1280,1536,1792,2048,2560,3072,
                           3584,4096,5120,6144,7168,8192,10240,12288,14336,
                           16384,20480,24576,28672,32768,40960,49152,57344,
                           65536/
                         # This is no longer relevant and should not be used.
                         # This array is only keep for legacy purposes.
mgscaleserial real /1.e39/ # Do not use - still being tested
mggrid_overlap integer /0/ # When 1, extra overlap is taken advantage of to
                           # reduce the amount of parallel communication
                           # Do not use, except for special purposes.
mgcoarsening integer /2/ # Amount of coarsening to do for each multigrid level
                         # Not fully supported. 2 is almost always the best
                         # value anyway.
mgusempistate logical /.true./
getmglevels(nx:integer,ny:integer,nz:integer,
            nxlocal:integer,nylocal:integer,nzlocal:integer,
            dx:real,dy:real,dz:real,
            conductors:ConductorType,fsdecomp:Decomposition)
   subroutine
   # Calculates levels of coarsening. Note that mglevels
   # must be zero when calling this routine.
multigrid3df(iwhich:integer,nx:integer,ny:integer,nz:integer,
             nxlocal:integer,nylocal:integer,nzlocal:integer,
             nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
             nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
             dx:real,dy:real,dz:real,
             phi(-nxguardphi:nxlocal+nxguardphi,
                 -nyguardphi:nylocal+nyguardphi,
                 -nzguardphi:nzlocal+nzguardphi):real,
             rho(-nxguardrho:nxlocal+nxguardrho,
                 -nyguardrho:nylocal+nyguardrho,
                 -nzguardrho:nzlocal+nzguardrho):real,
             rstar(-1:nzlocal+1):real,linbend:logical,
             bound0:integer,boundnz:integer,boundxy:integer,
             l2symtry:logical,l4symtry:logical,
             xmmin:real,ymmin:real,zmmin:real)
   subroutine
   # Solves Poisson's equation using the multigrid method. This uses variables
   # from the f3d package to control the iterations and conductors.
multigrid3dsolve(iwhich:integer,nx:integer,ny:integer,nz:integer,
                 nxlocal:integer,nylocal:integer,nzlocal:integer,
                 nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
                 nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                 dx:real,dy:real,dz:real,
                 phi(-nxguardphi:nxlocal+nxguardphi,
                     -nyguardphi:nylocal+nyguardphi,
                     -nzguardphi:nzlocal+nzguardphi):real,
                 rho(-nxguardrho:nxlocal+nxguardrho,
                     -nyguardrho:nylocal+nyguardrho,
                     -nzguardrho:nzlocal+nzguardrho):real,
                 rstar(-1:nzlocal+1):real,linbend:logical,bounds(0:5):integer,
                 xmmin:real,ymmin:real,zmmin:real,
                 mgparam:real,mgform:integer,mgiters:integer,
                 mgmaxiters:integer,mgmaxlevels:integer,mgerror:real,mgtol:real,
                 mgverbose:integer,
                 downpasses:integer,uppasses:integer,
                 lcndbndy:logical,laddconductor:logical,icndbndy:integer,
                 gridmode:integer,conductors:ConductorType,
                 lprecalccoeffs:logical,fsdecomp:Decomposition)
   subroutine
   # Solves Poisson's equation using the multigrid method. All input is
   # through the argument list.
multigrid2dsolve(iwhich:integer,nx:integer,nz:integer,
                 nxlocal:integer,nzlocal:integer,
                 nxguardphi:integer,nzguardphi:integer,
                 nxguardrho:integer,nzguardrho:integer,
                 dx:real,dz:real,
                 phi(-nxguardphi:nxlocal+nxguardphi,
                     -nzguardphi:nzlocal+nzguardphi):real,
                 rho(-nxguardrho:nxlocal+nxguardrho,
                     -nzguardrho:nzlocal+nzguardrho):real,
                 bounds(0:5):integer,
                 xmminlocal:real,
                 mgparam:real,mgform:integer,mgiters:integer,mgmaxiters:integer,
                 mgmaxlevels:integer,mgerror:real,mgtol:real,mgverbose:integer,
                 downpasses:integer,uppasses:integer,
                 lcndbndy:logical,laddconductor:logical,icndbndy:integer,
                 gridmode:integer,conductors:ConductorType,lrz:logical,
                 lmagnetostatic:logical,
                 fsdecomp:Decomposition)
   subroutine
   # Solves Poisson's equation in 2D using the multigrid method. All input is
   # through the argument list.
multigrid2ddielectricsolve(iwhich:integer,nx:integer,nz:integer,
                 nxlocal:integer,nzlocal:integer,
                 nxguardphi:integer,nzguardphi:integer,
                 nxguardrho:integer,nzguardrho:integer,
                 dx:real,dz:real,
                 phi(-nxguardphi:nxlocal+nxguardphi,
                     -nzguardphi:nzlocal+nzguardphi):real,
                 rho(-nxguardrho:nxlocal+nxguardrho,
                     -nzguardrho:nzlocal+nzguardrho):real,
                 epsilon(0:nxlocal+1,0:nzlocal+1):real,
                 bounds(0:5):integer,xmminlocal:real,
                 mgparam:real,mgiters:integer,mgmaxiters:integer,
                 mgmaxlevels:integer,mgerror:real,mgtol:real,mgverbose:integer,
                 downpasses:integer,uppasses:integer,
                 lcndbndy:logical,laddconductor:logical,
                 gridmode:integer,conductors:ConductorType,lrz:logical,
                 fsdecomp:Decomposition)
   subroutine
   # Solves Poisson's equation in 2D with a variable dielectric using the
   # multigrid method. All input is through the argument list.
mgsolveimplicites3d(iwhich:integer,nx:integer,ny:integer,nz:integer,
                    nxlocal:integer,nylocal:integer,nzlocal:integer,
                    dx:real,dy:real,dz:real,phi:real,rho:real,
                    ns:integer,qomdt:real,chi0:real,rstar:real,linbend:logical,
                    bounds(0:5):integer,
                    xmminlocal:real,ymminlocal:real,zmminlocal:real,
                    zgrid:real,mgparam:real,
                    mgiters:integer,mgmaxiters:integer,mgmaxlevels:integer,
                    mgerror:real,mgtol:real,mgverbose:integer,
                    downpasses:integer,uppasses:integer,lcndbndy:logical,
                    laddconductor:logical,icndbndy:integer,
                    gridmode:integer,conductors:ConductorType,
                    fsdecomp:Decomposition)
   subroutine
   # Solves Poisson's equation in 3D using the multigrid method including
   # the implicit susceptability. All input is through the argument list.
residual3d(nxlocal:integer,nylocal:integer,nzlocal:integer,
           nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
           nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
           nxguardres:integer,nyguardres:integer,nzguardres:integer,
           dxsqi:real,dysqi:real,dzsqi:real,
           phi(-nxguardphi:nxlocal+nxguardphi,
               -nyguardphi:nylocal+nyguardphi,
               -nzguardphi:nzlocal+nzguardphi):real,
           rho(-nxguardrho:nxlocal+nxguardrho,
               -nyguardrho:nylocal+nyguardrho,
               -nzguardrho:nzlocal+nzguardrho):real,
           res(-nxguardres:nxlocal+nxguardres,
               -nyguardres:nylocal+nyguardres,
               -nzguardres:nzlocal+nzguardres):real,
           mglevel:integer,bounds(0:5):integer,
           mgparam:real,mgform:integer,mgform2init:logical,
           lcndbndy:logical,icndbndy:integer,conductors:ConductorType,
           lprecalccoeffs:logical)
   subroutine
   # Calculates the residual for 3d arrays
residual2d(nxlocal:integer,nzlocal:integer,
           nxguardphi:integer,nzguardphi:integer,
           nxguardrho:integer,nzguardrho:integer,
           nxguardres:integer,nzguardres:integer,
           dxsqi:real,dzsqi:real,
           xminodx:real,lrz:logical,phi:real,rho:real,res:real,
           mglevel:integer,bounds(0:5):integer,
           lcndbndy:logical,icndbndy:integer,conductors:ConductorType)
   subroutine
   # Calculates the residual for 2d arrays
residual2ddielectric(nxlocal:integer,nzlocal:integer,
           nxguardphi:integer,nzguardphi:integer,
           nxguardrho:integer,nzguardrho:integer,
           nxguardres:integer,nzguardres:integer,
           phi:real,rho:real,epsilon:real,res:real,dx:real,dz:real,
           mglevel:integer,localbounds(0:5):integer,lcndbndy:logical,
           conductors:ConductorType,xminodx:real,lrz:logical)
   subroutine
   # Calculates the residual for 2d arrays with a variable dielectric
restrict3d(nx:integer,ny:integer,nz:integer,
           nxlocal:integer,nylocal:integer,nzlocal:integer,
           nxguardres:integer,nyguardres:integer,nzguardres:integer,
           res(-nxguardres:nxlocal+nxguardres,
               -nyguardres:nylocal+nyguardres,
               -nzguardres:nzlocal+nzguardres):real,
           nxcoarse:integer,nycoarse:integer,nzcoarse:integer,
           nxlocalcoarse:integer,nylocalcoarse:integer,nzlocalcoarse:integer,
           rhocoarse(0:nxlocalcoarse,0:nylocalcoarse,0:nzlocalcoarse):real,
           ff:real,bounds(0:5):integer,boundscoarse(0:5):integer,
           lxoffset:integer,lyoffset:integer,lzoffset:integer)
   subroutine
   # Restricts phi in 3 dimensions
expand3d(nx:integer,ny:integer,nz:integer,
         nxlocal:integer,nylocal:integer,nzlocal:integer,
         nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
         phi(-nxguardphi:nxlocal+nxguardphi,
             -nyguardphi:nylocal+nyguardphi,
             -nzguardphi:nzlocal+nzguardphi):real,
         nxcoarse:integer,nycoarse:integer,nzcoarse:integer,
         nxlocalcoarse:integer,nylocalcoarse:integer,nzlocalcoarse:integer,
         phicoarse(-1:nxlocalcoarse+1,
                   -1:nylocalcoarse+1,
                   -1:nzlocalcoarse+1):real,
         bounds(0:5):integer,lxoffset:integer,lyoffset:integer,lzoffset:integer,
         conductors:ConductorType,lprecalccoeffs:logical)
   subroutine
   # Expands phi in 3 dimensiosn
sorhalfpass3d(parity:integer,mglevel:integer,
              nxlocal:integer,nylocal:integer,nzlocal:integer,
              nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
              nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
              phi(-nxguardphi:nxlocal+nxguardphi,
                  -nyguardphi:nylocal+nyguardphi,
                  -nzguardphi:nzlocal+nzguardphi):real,
              rho(-nxguardrho:nxlocal+nxguardrho,
                  -nyguardrho:nylocal+nyguardrho,
                  -nzguardrho:nzlocal+nzguardrho):real,
              rstar(-1:nzlocal+1):real,dxsqi:real,dysqi:real,dzsqi:real,
              linbend:logical,bendx(-1:nxlocal):real,bounds(0:5):integer,
              mgparam:real,mgform:integer,
              lcndbndy:logical,icndbndy:integer,conductors:ConductorType,
              lprecalccoeffs:logical)
   subroutine
   # Performs one pass of SOR relaxation, either even or odd.
cond_potmg(interior:ConductorInteriorType,
           nxlocal:integer,nylocal:integer,nzlocal:integer,
           nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
           phi(-nxguardphi:nxlocal+nxguardphi,
               -nyguardphi:nylocal+nyguardphi,
               -nzguardphi:nzlocal+nzguardphi):real,
           mglevel:integer,mgform:integer,mgform2init:logical)
    subroutine
    # Sets voltage on interior of conductors
cond_zerorhointerior(interior:ConductorInteriorType,
                     nxlocal:integer,nylocal:integer,nzlocal:integer,
                     nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                     rho(-nxguardrho:nxlocal+nxguardrho,
                         -nyguardrho:nylocal+nyguardrho,
                         -nzguardrho:nzlocal+nzguardrho):real)
    subroutine
    # Sets rho to zero inside conductor points.
cond_sumrhointerior(interior:ConductorInteriorType,
                    nxlocal:integer,nylocal:integer,nzlocal:integer,
                    nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                    rho(-nxguardrho:nxlocal+nxguardrho,
                        -nyguardrho:nylocal+nyguardrho,
                        -nzguardrho:nzlocal+nzguardrho):real,
                    ixmin:integer,ixmax:integer,iymin:integer,iymax:integer,
                    izmin:integer,izmax:integer) real function
subcond_sumrhointerior(rhosum:real,interior:ConductorInteriorType,
                       nxlocal:integer,nylocal:integer,nzlocal:integer,
                       nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
                       rho(-nxguardrho:nxlocal+nxguardrho,
                           -nyguardrho:nylocal+nyguardrho,
                           -nzguardrho:nzlocal+nzguardrho):real,
                       ixmin:integer,ixmax:integer,iymin:integer,iymax:integer,
                       izmin:integer,izmax:integer) subroutine

setupconductorfielddata(nx:integer,ny:integer,nz:integer,
                        nxlocal:integer,nylocal:integer,nzlocal:integer,
                        dx:real,dy:real,dz:real,
                        conductors:ConductorType,fsdecomp:Decomposition)
                    subroutine
getefieldatconductorsubgrid(conductors:ConductorType,
                     dx:real,dy:real,dz:real,
                     nxlocal:integer,nylocal:integer,nzlocal:integer,
                     nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
                     phi(-nxguardphi:nxlocal+nxguardphi,
                         -nyguardphi:nylocal+nyguardphi,
                         -nzguardphi:nzlocal+nzguardphi):real,
                     bounds(0:5):integer)
                    subroutine
fixefieldatconductorpoints(conductors:ConductorType,
                           dx:real,dy:real,dz:real,
                           nxlocal:integer,nylocal:integer,nzlocal:integer,
                           nxguarde:integer,nyguarde:integer,nzguarde:integer,
                           field(-nxguarde:nxlocal+nxguarde,
                                 -nyguarde:nylocal+nyguarde,
                                 -nzguarde:nzlocal+nzguarde):real)
                    subroutine
    # Fix the two point finite difference operator at conductor points.
sete3dwithconductor(conductors:ConductorType,
                    np:integer,xp(np):real,yp(np):real,zp(np):real,
                    efetch:integer,ex(np):real,ey(np):real,ez(np):real,
                    zgrid:real,xmmin:real,ymmin:real,zmmin:real,
                    dx:real,dy:real,dz:real,
                    nxlocal:integer,nylocal:integer,nzlocal:integer,
                    nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
                    nxguarde:integer,nyguarde:integer,nzguarde:integer,
                    phi(-nxguardphi:nxlocal+nxguardphi,
                        -nyguardphi:nylocal+nyguardphi,
                        -nzguardphi:nzlocal+nzguardphi):real,
                    selfe:real,
                    l2symtry:logical,l4symtry:logical,lcylindrical:logical)
                    subroutine
sete3dongridwithconductor(conductors:ConductorType,phi:real,
                      dx:real,dy:real,dz:real,
                      nx:integer,ny:integer,nz:integer,
                      ex:real,ey:real,ez:real,
                      nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
                      bounds(0:5):integer)
                      subroutine

*********** MultigridBE3d dump:
multigridbe3df(iwhich:integer,nx:integer,ny:integer,nzlocal:integer,nz:integer,
             dx:real,dy:real,dz:real,phi:real,rho:real,
             rstar:real,linbend:logical,
             bound0:integer,boundnz:integer,boundxy:integer,
             l2symtry:logical,l4symtry:logical,
             xmmin:real,ymmin:real,zmmin:real,
             iondensity:real,electrontemperature:real,plasmapotential:real,
             electrondensitymaxscale:real)
   subroutine
   # Solves Poisson's equation using the multigrid method, including the
   # Boltzmann electron term. This uses variables
   # from the f3d package to control the iterations and conductors.
multigridbe3dsolve(iwhich:integer,
             nx:integer,ny:integer,nz:integer,
             dx:real,dy:real,dz:real,phi:real,rho:real,
             rstar:real,linbend:logical,bounds(0:5):integer,
             xmmin:real,ymmin:real,zmmin:real,
             mgparam:real,mgiters:integer,mgmaxiters:integer,
             mgmaxlevels:integer,mgerror:real,mgtol:real,mgverbose:integer,
             downpasses:integer,uppasses:integer,
             lcndbndy:logical,laddconductor:logical,icndbndy:integer,
             gridmode:integer,conductors:ConductorType,
             iondensitygrid3d:Grid3dtype,fsdecomp:Decomposition)
   subroutine
   # Solves Poisson's equation using the multigrid method, including the
   # Boltzmann electron term. All input is through the argument list.
setupiondensitygrid3d(xmmin:real,ymmin:real,zmmin:real,dx:real,dy:real,dz:real,
                      nx:integer,ny:integer,nzlocal:integer,
                      rho:real,iondensitygrid3d:Grid3dtype)
   subroutine # Sets up the iondensity grid for the 3d Boltzman electron
              # field solver.


*********** Multigrid3d_work:
# Temporary variables and array used by subgrid_sor_to_mg
wnx integer
wny integer
wnz integer
iii(0:wnx,0:wny,0:wnz) _integer


%%%%%%%%%%% BFieldGridType:
nx      integer  /0/  # Number of grid cells in x in B grid
ny      integer  /0/  # Number of grid cells in y in B grid
nz      integer  /0/  # Number of grid cells in z in B grid
nxlocal integer  /0/  # Number of grid cells in x in B grid
nylocal integer  /0/  # Number of grid cells in y in B grid
nzlocal integer  /0/  # Number of grid cells in z in B grid
nxguarda integer /1/  # Number of guard cells in x for the potential
nyguarda integer /1/  # Number of guard cells in y for the potential
nzguarda integer /1/  # Number of guard cells in z for the potential
nxguardj integer /0/  # Number of guard cells in x for the current density
nyguardj integer /0/  # Number of guard cells in y for the current density
nzguardj integer /0/  # Number of guard cells in z for the current density
nxguardb integer /0/  # Number of guard cells in x for the B field
nyguardb integer /0/  # Number of guard cells in y for the B field
nzguardb integer /0/  # Number of guard cells in z for the B field
dx      real [m] /0./ # x grid cell size in B grid
dy      real [m] /0./ # y grid cell size in B grid
dz      real [m] /0./ # z grid cell size in B grid
xmmin   real [m] /0./ # X lower limit of mesh
xmmax   real [m] /0./ # X upper limit of mesh
ymmin   real [m] /0./ # Y lower limit of mesh
ymmax   real [m] /0./ # Y upper limit of mesh
zmmin   real [m] /0./ # Z lower limit of mesh
zmmax   real [m] /0./ # Z upper limit of mesh
xmminlocal real [m]  # Local value of xmmin
xmmaxlocal real [m]  # Local value of xmmax
ymminlocal real [m]  # Local value of ymmin
ymmaxlocal real [m]  # Local value of ymmax
zmminlocal real [m]  # Local value of zmmin
zmmaxlocal real [m]  # Local value of zmmax

lcndbndy logical /.true./ # Turns on sub-grid boundaries
icndbndy integer /2/      # Type of interpolant to use for sub-grid boundaries
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
mgverbose(0:2)   integer /1/   # Level of verbosity of bfield MG solver
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

j(0:2,-nxguardj:nxlocal+nxguardj,
      -nyguardj:nylocal+nyguardj,
      -nzguardj:nzlocal+nzguardj) _real # Current density
b(0:2,-nxguardb:nxlocal+nxguardb,
      -nyguardb:nylocal+nyguardb,
      -nzguardb:nzlocal+nzguardb) _real # B field, calculated from B = del cross A
a(0:2,-nxguarda:nxlocal+nxguarda,
      -nyguarda:nylocal+nyguarda,
      -nzguarda:nzlocal+nzguarda) _real # Magnetic vector potential
nxold     integer  /0/  # Number of grid cells in x in B grid
nyold     integer  /0/  # Number of grid cells in y in B grid
nzold     integer  /0/  # Number of grid cells in z in B grid
aold(0:2,-nxguarda:nxold+nxguarda,
         -nyguarda:nyold+nyguarda,
         -nzguarda:nzold+nzguarda) _real
  # Magnetic vector potential, calculated from del sq A = J

attx(0:nxlocal-1)     _real           # Attenuation factor as fcn. of kx
atty(0:nylocal-1)     _real           # Attenuation factor as fcn. of ky
attz(0:nz)   _real           # Attenuation factor as fcn. of kz
kxsq(0:nxlocal-1)     _real [1/m**2]  # Discrete analog to kx^2/4Pi
kysq(0:nylocal-1)     _real [1/m**2]  # Discrete analog to ky^2/4Pi
kzsq(0:nz)   _real [1/m**2]  # Discrete analog to kz^2/4Pi

rstar(-1:nzlocal+1)         _real [m] # Radius of curv of reference orbit
scrtch(2*nxlocal+2*nylocal)      _real     # Scratch for fieldsolve
xywork(2,0:nxlocal,0:nylocal)    _real     # Work space for transverse FFTs
zwork(2,0:nxlocal,0:nz) _real     # Work space used to optimize vsftz

nsjtmp integer /0/
jsjtmp(0:nsjtmp-1) _integer /-1/ #
nsndtsj integer /0/
jtmp(3,-nxguardj:nxlocal+nxguardj,
       -nyguardj:nylocal+nyguardj,
       -nzguardj:nzlocal+nzguardj,0:nsndtsj-1) _real
             # Temporary copy of the current density from the particles.

conductors ConductorType # Default data structure for conductor data

*********** BFieldGrid dump:
bfield BFieldGridType
bfieldp BFieldGridType
ljfinalized logical /.true./
                 # Flags whether the boundary conditions and other processing
                 # of bfield.j has been done, so that its ready to be used in a
                 # field solve. Anytime bfieldp.j is changed, this flag is set
                 # to false. Just before a field solve, the flag is checked,
                 # and if needed, the appropriate operations done.

init_bfieldsolver(bfstype:integer) subroutine # Initializes the B-field solver
bfieldsol3d(iwhich) subroutine # Self B-field solver
loadj3d(pgroup:ParticleGroup,ins:integer,nps:integer,is:integer,lzero:logical)
             subroutine # Provides a simple interface to the current density
                        # loading routine setj3d
applyjboundaryconditions(bfield:BFieldGridType,
                         boundxy:integer,bound0:integer,boundnz:integer)
             subroutine #
setb3d(b(0:2,0:nxlocal,0:nylocal,0:nzlocal):real,
       np:integer,xp:real,yp:real,zp:real,zgrid:real,
       bx:real,by:real,bz:real,nxlocal:integer,nylocal:integer,nzlocal:integer,
       nxguardb:integer,nyguardb:integer,nzguardb:integer,
       dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,
       l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine #
fetchafrompositions3d(a(0:2,-nxguarda:nxlocal+nxguarda,
                            -nyguarda:nylocal+nyguarda,
                            -nzguarda:nzlocal+nzguarda):real,
                      np:integer,
                      xp(np):real,yp(np):real,zp(np):real,ap(0:2,np):real,
                      zgrid:real,
                      nxlocal:integer,nylocal:integer,nzlocal:integer,
                      nxguarda:integer,nyguarda:integer,nzguarda:integer,
                      dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,
                      l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine #
getbfroma3d(a(0:2,-nxguarda:nxlocal+nxguarda,
                  -nyguarda:nylocal+nyguarda,
                  -nzguarda:nzlocal+nzguarda):real,
            b(0:2,-nxguardb:nxlocal+nxguardb,
                  -nyguardb:nylocal+nyguardb,
                  -nzguardb:nzlocal+nzguardb):real,
            nxlocal:integer,nylocal:integer,nzlocal:integer,
            nxguarda:integer,nyguarda:integer,nzguarda:integer,
            nxguardb:integer,nyguardb:integer,nzguardb:integer,
            dx:real,dy:real,dz:real,
            xmmin:real,lcylindrical:logical,lusevectorpotential:logical)
             subroutine #
curl3d(a:real,b:real,nx:integer,ny:integer,nzlocal:integer,
       nxguarda:integer,nyguarda:integer,nzguarda:integer,
       nxguardb:integer,nyguardb:integer,nzguardb:integer,
       dx:real,dy:real,dz:real,
       xmmin:real,lcylindrical:logical)
             subroutine # Calculates the curl of the input array a, putting
                        # the result into b.
setj3d(j(0:2,-nxguardj:nxlocal+nxguardj,
             -nyguardj:nylocal+nyguardj,
             -nzguardj:nzlocal+nzguardj):real,
       j1d:real,np:integer,xp:real,yp:real,zp:real,
       zgrid:real,uxp:real,uyp:real,uzp:real,gaminv:real,q:real,w:real,
       nw:integer,wghtp:real,depos:string,
       nxlocal:integer,nylocal:integer,nzlocal:integer,
       nxguardj:integer,nyguardj:integer,nzguardj:integer,
       dx:real,dy:real,dz:real,xmmin:real,ymmin:real,zmmin:real,
       l2symtry:logical,l4symtry:logical,lcylindrical:logical)
             subroutine # Computes current density
getjforfieldsolve()
             subroutine #
setjforfieldsolve3d(nxlocal:integer,nylocal:integer,nzlocal:integer,
                    j(0:2,-nxguardj:nxlocal+nxguardj,
                          -nyguardj:nylocal+nyguardj,
                          -nzguardj:nzlocal+nzguardj):real,
                    nxp:integer,nyp:integer,nzp:integer,
                    jp(0:2,-nxguardj:nxp+nxguardj,
                           -nyguardj:nyp+nyguardj,
                           -nzguardj:nzp+nzguardj):real,
                    nxguardj:integer,nyguardj:integer,nzguardj:integer,
                    fsdecomp:Decomposition,ppdecomp:Decomposition)
             subroutine #
setupbfieldsforparticles3d(ns:integer,ndts:integer,it:integer,
                           bfield:BFieldGridType,bfieldp:BFieldGridType)
             subroutine #
fetchb3dfrompositions(jsid:integer,ndts:integer,n:integer,x(n):real,y(n):real,z(n):real,
                      bx(n):real,by(n):real,bz(n):real)
             subroutine #
fetcha(n:integer,x(n):real,y(n):real,z(n):real,a(n):real)
             subroutine #
getbforparticles()
             subroutine #
bvp3d(iwhich:integer,bfstype:integer)
             subroutine #
getanalyticbtheta(b(0:2,-nxguardb:nxlocal+nxguardb,
                        -nyguardb:nylocal+nyguardb,
                        -nzguardb:nzlocal+nzguardb):real,
                  j(0:2,-nxguardj:nxlocal+nxguardj,
                        -nyguardj:nylocal+nyguardj,
                        -nzguardj:nzlocal+nzguardj):real,
                  nxlocal:integer,nylocal:integer,nzlocal:integer,
                  nxguardb:integer,nyguardb:integer,nzguardb:integer,
                  nxguardj:integer,nyguardj:integer,nzguardj:integer,
                  dx:real,xmmin:real)
             subroutine # Calculate and analytic Btheta

*********** AMR3droutines:
gathersourcefromchild(rho:real,nc:integer,ng:integer,nn:integer,
                      childrho:real,cnn:integer,
                      l:integer,u:integer,fulllower:integer,
                      childlower:integer,childupper:integer,
                      r:integer,weights:real,
                      radius:real,cradius:real,lcylinderical:logical)
      subroutine
gatherpotentialfromparents(potential:real,nc:integer,ng:integer,nn:integer,
                           l:integer,u:integer,fulllower:integer,
                           parentpotential:real,pnn:integer,
                           parentlower:integer,r:integer)
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
nzlocallan integer
nxtranlan integer
nytranlan integer
nztranlan integer
alan(0:nxlan,0:nylan,0:nzlocallan) _real
blan(0:nxlan,0:nylan,0:nzlocallan) _real
clan(0:nxlan,0:nylan,0:nzlocallan) _real
atranlan(0:nxtranlan,0:nytranlan-1,0:nzlan) _real
btranlan(0:nxtranlan,0:nytranlan-1,0:nzlan) _real
ctranlan(0:nxtranlan,0:nytranlan-1,0:nzlan) _real
dtranlan(0:nxlan,0:nylan,0:2) _real

*********** F3Dsubs:
#  Callable subroutines in the F3D package
vpois3d(iwhich:integer,a:real,ak:real,kxsq:real,kysq:real,kzsq:real,
        attx:real,atty:real,attz:real,filt:real,
        lx:real,ly:real,lz:real,nx:integer,ny:integer,nzlocal:integer,
        nz:integer,
        nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
        w:real,xywork:real,zwork:real,ibc:integer,
        l2symtry:logical,l4symtry:logical,
        bound0:integer,boundnz:integer,boundxy:integer)
     subroutine #  The 3d Poisson solver
vp3x(iwhich) subroutine
     # Python-level interface to VPOIS3d, using FS3 database variables
     #  The user program should declare a similar subroutine w/ its vars.
vpois2d(iwhich,a:real,ak:real,kxsq:real,kysq:real,attx:real,atty:real,
        filt:real,lx:real,ly:real,nx:integer,ny:integer,
        nxguardphi:integer,nyguardphi:integer,
        work:real,xywork:real,ibc,
        l2symtry:logical,l4symtry:logical)
     subroutine #  The 2d Sine-Sine Poisson solver
vsftx(a:real,work:real,cp:real,cm:real,nx:integer,ny:integer,
      nxguardphi:integer,nyguardphi:integer,isetup)
     subroutine #  Vectorized Sine Fourier Transform in X
vsfty(a:real,work:real,cp:real,cm:real,nx:integer,ny:integer,
      nxguardphi:integer,nyguardphi:integer,isetup)
     subroutine #  Vectorized Sine Fourier Transform in Y
vpftx(nx:integer,ny:integer,nz:integer,
      nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
      a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in X
vpfty(nx:integer,ny:integer,nz:integer,
      nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
      a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Y
vpftz(nx:integer,ny:integer,nz:integer,
      nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
      a:real,norm:real,esx,esy,nxy,zwork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Z
vpftxi(nx:integer,ny:integer,nz:integer,
       nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
       a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in X
vpftyi(nx:integer,ny:integer,nz:integer,
       nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
       a:real,norm:real,esx,esy,nxy,xywork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Y
vpftzi(nx:integer,ny:integer,nz:integer,
       nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
       a:real,norm:real,esx,esy,nxy,zwork:real)
     subroutine #  Vectorized Real Periodic Fourier Transform in Z
cosqx(a:real,w:real,c:real,nx:integer,ny:integer,
      nxguardphi:integer,nyguardphi:integer,
      isetup:integer,isign:integer) subroutine
cosqy(a:real,w:real,c:real,nx:integer,ny:integer,
      nxguardphi:integer,nyguardphi:integer,
      isetup:integer,isign:integer) subroutine
vcpft(r:real,i:real,n:integer,incp:integer,signp:integer,lenv:integer,
      lfd:integer) subroutine
attenuate(nx:integer,ny:integer,nz:integer,
          nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
          a:real,
          attx:real,atty:real,attz:real,
          ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
unattenuate(nx:integer,ny:integer,nz:integer,
            nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
            a:real,
            attx:real,atty:real,attz:real,
            ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
rhotophi(nx:integer,ny:integer,nz:integer,
         nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
         a:real,kxsq:real,kysq:real,kzsq:real,
         ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
phitorho(nx:integer,ny:integer,nz:integer,
         nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
         a:real,kxsq:real,kysq:real,kzsq:real,
         ikxmin:integer,ikymin:integer,esx:integer,esy:integer,esz:integer)
     subroutine
pipe3df(iwhich,pipeshpe:real,rho:real,phi:real,kxsq:real,kysq:real,
        kzsq:real,attx:real,atty:real,attz:real,filt:real,
        xlen:real,ylen:real,zlen:real,nx:integer,ny:integer,nzlocal:integer,
        nz:integer,
        nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
        scrtch:real,
        l2symtry:logical,l4symtry:logical)
     subroutine #  External interface to capacity field solver
pipest3d(pipeshpe:real,cap:real,phi:real,kxsq:real,kysq:real,kzsq:real,
         attx:real,atty:real,attz:real,filt:real,xlen:real,
         ylen:real,zlen:real,
         nx:integer,ny:integer,nzlocal:integer,nz:integer,
         nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
         scrtch:real,piperadi:real,
         pipen:real,pipe8th:real,pipex:real,pipey:real,cap3d:real,kpvt:integer,
         l2symtry:logical,l4symtry:logical)
     subroutine #  Sets up capacity matrix
setcap3d(ncndpts,quadx:real,quady:real,quadz:real,quadcap:real,kpvt:integer,
         nx:integer,ny:integer,nzlocal:integer,nz:integer,
         nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
         phi:real,work:real,kxsq:real,kysq:real,kzsq:real,attx:real,atty:real,
         attz:real,xlen:real,ylen:real,zlen:real,filt:real,
         l2symtry:logical,l4symtry:logical)
     subroutine # Sets up 3-d cap matrix
capfs(iwhich,nx:integer,ny:integer,nzlocal:integer,nz:integer,
      nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
      nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
      phi:real,rho:real,xlen:real,ylen:real,zlen:real,
      kxsq:real,kysq:real,kzsq:real,attx:real,atty:real,attz:real,
      filt:real,work:real,workz:real,xmmax:real,
      pipeshpe:string,periinz:logical,l2symtry:logical,l4symtry:logical)
     subroutine # 3-d field solve using cap matrices
cndctr3d(nx:integer,ny:integer,xmmax:real,dx:real,dy:real,
         quadradi:real,quadcent:real,ncndpts,
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
       nx:integer,ny:integer,nzlocal:integer,nz:integer,
       nxguardphi:integer,nyguardphi:integer,nzguardphi:integer,
       nxguardrho:integer,nyguardrho:integer,nzguardrho:integer,
       scrtch:real,xmmax:real,zmmin:real,zgrid:real,
       pipeshpe:string,periinz:logical,l2symtry:logical,l4symtry:logical)
     subroutine # External routine for capacity matrix field solve
     
%%%%%%%% ConductorInterceptType:
mglevel integer
xmmin real
ymmin real
zmmin real
dx real
dy real
dz real
nx integer
ny integer
nz integer
ix integer
iy integer
iz integer
nxicpt integer
nyicpt integer
nzicpt integer
xintercepts(nxicpt,0:ny,0:nz) _real
yintercepts(nyicpt,0:nx,0:nz) _real
zintercepts(nzicpt,0:nx,0:ny) _real
xvoltages((nxicpt+1)/2,0:ny,0:nz) _real
yvoltages((nyicpt+1)/2,0:nx,0:nz) _real
zvoltages((nzicpt+1)/2,0:nx,0:ny) _real
xcondids((nxicpt+1)/2,0:ny,0:nz) _integer
ycondids((nyicpt+1)/2,0:nx,0:nz) _integer
zcondids((nzicpt+1)/2,0:nx,0:ny) _integer

******** ConductorGeometryGenerators:
conductordelfromintercepts(intercepts:ConductorInterceptType,
                           conductors:ConductorType,
                           dfill:real,fuzz:real) subroutine
intercepts_or(ileft:ConductorInterceptType,iright:ConductorInterceptType,
              iresult:ConductorInterceptType) subroutine
intercepts_and(ileft:ConductorInterceptType,iright:ConductorInterceptType,
               iresult:ConductorInterceptType) subroutine
intercepts_not(ileft:ConductorInterceptType,
               iresult:ConductorInterceptType) subroutine
zplaneconductorfnew(z0:real,zsign:real,xcent:real,ycent:real,zcent:real,
                    intercepts:ConductorInterceptType,fuzz:real) subroutine
planeconductorfnew(z0:real,zsign:real,theta:real,phi:real,
                   xcent:real,ycent:real,zcent:real,
                   intercepts:ConductorInterceptType,fuzz:real) subroutine
zcylinderconductorfnew(rad:real,length:real,xcent:real,ycent:real,zcent:real,
                       intercepts:ConductorInterceptType,fuzz:real) subroutine
zroundedcylinderconductorfnew(rad:real,length:real,rad2:real,
                              xcent:real,ycent:real,zcent:real,
                       intercepts:ConductorInterceptType,fuzz:real) subroutine
zcylinderoutconductorfnew(rad:real,length:real,xcent:real,ycent:real,zcent:real,
                       intercepts:ConductorInterceptType,fuzz:real) subroutine
zroundedcylinderoutconductorfnew(rad:real,length:real,rad2:real,
                                 xcent:real,ycent:real,zcent:real,
                       intercepts:ConductorInterceptType,fuzz:real) subroutine
boxconductorfnew(xsize:real,ysize:real,zsize:real,
                 xcent:real,ycent:real,zcent:real,
                 intercepts:ConductorInterceptType,fuzz:real) subroutine
sphereconductorfnew(rad:real,xcent:real,ycent:real,zcent:real,
                    intercepts:ConductorInterceptType,fuzz:real) subroutine
zsrfrvconductorfnew(nn:integer,rsrf(nn):real,zsrf(nn):real,
                    rad(nn-1):real,rc(nn-1):real,zc(nn-1):real,
                    xcent:real,ycent:real,zcent:real,
                    intercepts:ConductorInterceptType,fuzz:real) subroutine
zsrfrvintercept(nn:integer,rsrf(nn):real,zsrf(nn):real,
                rad(nn-1):real,rc(nn-1):real,zc(nn-1):real,
                xcent:real,ycent:real,zcent:real,
                n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
                xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine

solvequartic(a0:real,a1:real,a2:real,a3:real,x1:complex,x2:complex,x3:complex,x4:complex) subroutine
setconductorparity(nn:integer,ix:integer,iy:integer,iz:integer,
                   dels:real,parity:integer,fuzz0:real,fuzz1:real,fuzzsign:integer,
                   dfill:real) subroutine
zplaneconductorf(z0:real,zsign:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zplaneconductord(z0:real,zsign:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zplaneintercept(z0:real,zsign:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
planeconductorf(z0:real,zsign:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
planeconductord(z0:real,zsign:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
planeintercept(z0:real,zsign:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
boxconductorf(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
boxconductord(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
boxintercept(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
cylinderconductorf(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
cylinderconductord(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
cylinderintercept(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
cylindersconductorf(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
cylindersconductord(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
cylindersintercept(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zcylinderconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zcylinderconductord(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zcylinderintercept(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zroundedcylinderconductorf(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zroundedcylinderconductord(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zroundedcylinderintercept(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zcylinderoutconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zcylinderoutconductord(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zcylinderoutintercept(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zroundedcylinderoutconductorf(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zroundedcylinderoutconductord(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zroundedcylinderoutintercept(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
sphereconductorf(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
sphereconductord(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
sphereintercept(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
coneconductorf(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
coneconductord(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
coneintercept(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
conesconductorf(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
conesconductord(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
conesintercept(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
annulusconductorf(rmin:real,rmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
annulusconductord(rmin:real,rmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
annulusintercept(rmin:real,rmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
ztorusconductorf(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
ztorusconductord(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
ztorusintercept(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zgridconductorf(xcellsize:real,ycellsize:real,length:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zgridconductord(xcellsize:real,ycellsize:real,length:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zgridintercept(xcellsize:real,ycellsize:real,length:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
beamletplateconductorf(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
beamletplateconductord(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
beamletplateintercept(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zsrfrvconductord(nn:integer,rsrf(nn):real,zsrf(nn):real,
        rad(nn-1):real,rc(nn-1):real,zc(nn-1):real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zsrfrvoutconductorf(lrofzfunc:logical,zmin:real,zmax:real,rmax:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zsrfrvoutconductord(lrofzfunc:logical,zmin:real,zmax:real,rmax:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zsrfrvoutintercept(lrofzfunc:logical,zmin:real,zmax:real,rmax:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zsrfrvinconductorf(lrofzfunc:logical,zmin:real,zmax:real,rmin:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zsrfrvinconductord(lrofzfunc:logical,zmin:real,zmax:real,rmin:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zsrfrvinintercept(lrofzfunc:logical,zmin:real,zmax:real,rmin:real,griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
zsrfrvinoutconductorf(lrminofz:logical,lrmaxofz:logical,zmin:real,zmax:real,
        griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,delmx(n):real,delpx(n):real,
        delmy(n):real,delpy(n):real,delmz(n):real,delpz(n):real,
        fuzz:real) subroutine
zsrfrvinoutconductord(lrminofz:logical,lrmaxofz:logical,zmin:real,zmax:real,
        griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,distance(n):real) subroutine
zsrfrvinoutintercept(lrminofz:logical,lrmaxofz:logical,zmin:real,zmax:real,
        griddz:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x(n):real,y(n):real,z(n):real,vx:real,vy:real,vz:real,
        xi:real,yi:real,zi:real,itheta:real,iphi:real) subroutine
srfrv_f(zz:real,lrofzfunc:logical,icase:integer,izflag:integer)
     real function # Function used to calculate r of z given some input.


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
timemultigrid2dsolve real /0./
timemultigridrzsolve real /0./
timegathersourcefromchild real /0./
timegatherpotentialfromparents real /0./
timegatherafromparents real /0./
timegatherjfromchild real /0./
timecheckconductors real /0./

timevcycle real /0./
timerestrict2d real /0./
timerestrict3d real /0./
timeexpand2d real /0./
timeexpand3d real /0./
timesorpass3d real /0./
timesorhalfpass3d real /0./
timesorhalfpass3d1 real /0./
timesorhalfpass3d2 real /0./
timesorhalfpass3d3 real /0./
timesorhalfpass3d4 real /0./
timesorhalfpass3d5 real /0./
timemgsor_loop_bend real /0./
timecond_potmg real /0./
timecond_potmgres real /0./
timecond_zerorhointerior real /0./
timecond_sumrhointerior real /0./
timesubcond_sumrhointerior real /0./
timecondbndymg real /0./
timecondbndymgint real /0./
timecondbndyres real /0./
timeresidual real /0./

timeexchange_phi            real /0./
timetranspose               real /0./
timetransposei              real /0./
timelantzsolver             real /0./
timegeneraltridiag          real /0./
timeparalleltridiag         real /0./
timeparallelgatherall       real /0./
timemgdividenz              real /0./
timemggetexchangepes        real /0./
timegetmessagedata          real /0./
timemgexchange_phi          real /0./
timemgexchange_phi_periodic real /0./
timemgexchange_res          real /0./
timemgexchange_phiupdate    real /0./
timeprintarray3d            real /0./

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
timegetbforparticles3d         real /0./
timegetaforfields3d            real /0./

timegetbfieldsongrid3d real /0./
timegetcoefficientsimplicites3d real /0./
timegetcoefficientssubgridimplicites3d real /0./
timemgsolveimplicites3d real /0./
timevcycleimplicites3d real /0./
timerelaximplicites3d real /0./
timeresidualimplicites3d real /0./
timeaverageperiodicphi3d real /0./
timeapplyboundaryconditions3d real /0./
timeapplyparallelboundaryconditions3d real /0./


****** RRRRRR:
rrr(:,:,:) _real
ppp(:,:,:) _real
qqq(:,:,:) _real
ncalculate integer /0/
