f3d
#@(#) File F3D.V, version $Revision: 3.69 $, $Date: 2003/03/25 22:55:02 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package F3D of code WARP6
# fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
{
}

*********** F3Dversion:
versf3d character*19 /"$Revision: 3.69 $"/#  Code version version is set by CVS

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
zwork(2,0:nx,0:nz)        _real        #  Workspace to optimize vsftz

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

*********** PSOR3d dump parallel:
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
lzerophiedge logical /.true./ # When true and when gridmode == 0, the edge of
                       # the phi array is zeroed out. This clears phi at any
                       # conductor points on the edge of the mesh.
bound0    integer /0/  # Type of boundary condition at plane z=0
                       # 0 is constant potential, 1 is zero normal derivative,
                       # and 2 is periodic
boundnz   integer /0/  # Type of boundary condition at plane z=nz
                       # 0 is constant potential, 1 is zero normal derivative,
                       # and 2 is periodic
boundxy   integer /0/  # Type of boundary condition at sides
                       # 0 is constant potential, 1 is zero normal derivative,
                       # and 2 is periodic
zparity   integer /0/  # iz parity, used in the parallel version so that the
                       # parity of the subgrid conductor points can be
                       # relative to the full grid.
dxfine    real         # Size of transverse grid cells are finest level
lplates          logical  # Sets whether or not quadruple endplates are included
rodfract        real /1./ # Fraction of quadrupole rod which is used

lcndbndy logical /.true./ # Turns on sub-grid boundaries
icndbndy integer /1/      # Type of interpolant to use for sub-grid boundaries
                          # 1 egun style
                          # 2 EBC style (non-centered finite-difference)

*********** Conductor3d dump parallel:
icstart(0:100)  integer -dump # Start of the conductor points for each MG level
ecstart(0:100)  integer -dump # Start of the even conductor points for each MG level
ocstart(0:100)  integer -dump # Start of the odd conductor points for each MG level
ncondmax          integer # Maximum number of points in conductor
ncond             integer # Number of points within conductors
ixcond(ncondmax) _integer # X coordinate of points in conductor
iycond(ncondmax) _integer # Y coordinate of points in conductor
izcond(ncondmax) _integer # Y coordinate of points in conductor
condvolt(ncondmax) _real  # voltage of points in conductor
condnumb(ncondmax) _integer # Number of the conductor the points are in
icondlevel(ncondmax) _integer /-1/ # Coarseness level at which the point is on grid
icondlxy(ncondmax) _integer # Obsolete array, only used to recover old datasets
icondlz(ncondmax) _integer # Obsolete array, only used to recover old datasets

fuzzsign      integer /-1/ # When -1, subgrid points with distances == 1 are
                           # skipped, when +1 not skipped.
ncndmax       integer /0/   # Maximum number of points for sub-grid boundaries
necndbdy      integer /0/   # Number of points for even sub-grid boundaries
ecndpvph(ncndmax)     _real -dump # Saves phi for even sub-grid boundaries
iecndx  (ncndmax)  _integer # location of points for even sub-grid boundaries
iecndy  (ncndmax)  _integer # location of points for even sub-grid boundaries
iecndz  (ncndmax)  _integer # location of points for even sub-grid boundaries
ecdelmx (ncndmax)     _real # distance in x of surface with lower x, even
ecdelmy (ncndmax)     _real # distance in y of surface with lower y, even
ecdelmz (ncndmax)     _real # distance in z of surface with lower z, even
ecdelpx (ncndmax)     _real # distance in x of surface with higher x, even
ecdelpy (ncndmax)     _real # distance in y of surface with higher y, even
ecdelpz (ncndmax)     _real # distance in z of surface with higher z, even
ecvolt  (ncndmax)     _real # voltage of points for even sub-grid boundaries
ecvoltmx(ncndmax)     _real # Voltage on conductor in minus x direction, even
ecvoltpx(ncndmax)     _real # Voltage on conductor in plus  x direction, even
ecvoltmy(ncndmax)     _real # Voltage on conductor in minus y direction, even
ecvoltpy(ncndmax)     _real # Voltage on conductor in plus  y direction, even
ecvoltmz(ncndmax)     _real # Voltage on conductor in minus z direction, even
ecvoltpz(ncndmax)     _real # Voltage on conductor in plus  z direction, even
ecnumb  (ncndmax)  _integer # Number of the conductor the even points are in
ecnumbmx(ncndmax)  _integer # Number of the conductor in minus x direction, even
ecnumbpx(ncndmax)  _integer # Number of the conductor in plus  x direction, even
ecnumbmy(ncndmax)  _integer # Number of the conductor in minus y direction, even
ecnumbpy(ncndmax)  _integer # Number of the conductor in plus  y direction, even
ecnumbmz(ncndmax)  _integer # Number of the conductor in minus z direction, even
ecnumbpz(ncndmax)  _integer # Number of the conductor in plus  z direction, even
iecndlevel(ncndmax)  _integer /-1/ # Coarseness level at which the point is on grid
iecndlxy(ncondmax) _integer # Obsolete array, only used to recover old datasets
iecndlz(ncondmax) _integer # Obsolete array, only used to recover old datasets

nocndbdy      integer /0/   # Number of points for odd sub-grid boundaries
ocndpvph(ncndmax)     _real -dump # Saves phi for odd sub-grid boundaries
iocndx  (ncndmax)  _integer # location of points for odd sub-grid boundaries
iocndy  (ncndmax)  _integer # location of points for odd sub-grid boundaries
iocndz  (ncndmax)  _integer # location of points for odd sub-grid boundaries
ocdelmx (ncndmax)     _real # distance in x of surface with lower x, odd
ocdelmy (ncndmax)     _real # distance in y of surface with lower y, odd
ocdelmz (ncndmax)     _real # distance in z of surface with lower z, odd
ocdelpx (ncndmax)     _real # distance in x of surface with higher x, odd
ocdelpy (ncndmax)     _real # distance in y of surface with higher y, odd
ocdelpz (ncndmax)     _real # distance in z of surface with higher z, odd
ocvolt  (ncndmax)     _real # voltage of points for odd sub-grid boundaries
ocvoltmx(ncndmax)     _real # Voltage on conductor in minus x direction, odd
ocvoltpx(ncndmax)     _real # Voltage on conductor in plus  x direction, odd
ocvoltmy(ncndmax)     _real # Voltage on conductor in minus y direction, odd
ocvoltpy(ncndmax)     _real # Voltage on conductor in plus  y direction, odd
ocvoltmz(ncndmax)     _real # Voltage on conductor in minus z direction, odd
ocvoltpz(ncndmax)     _real # Voltage on conductor in plus  z direction, odd
ocnumb  (ncndmax)  _integer # Number of the conductor the odd points are in
ocnumbmx(ncndmax)  _integer # Number of the conductor in minus x direction, odd
ocnumbpx(ncndmax)  _integer # Number of the conductor in plus  x direction, odd
ocnumbmy(ncndmax)  _integer # Number of the conductor in minus y direction, odd
ocnumbpy(ncndmax)  _integer # Number of the conductor in plus  y direction, odd
ocnumbmz(ncndmax)  _integer # Number of the conductor in minus z direction, odd
ocnumbpz(ncndmax)  _integer # Number of the conductor in plus  z direction, odd
iocndlevel(ncndmax)  _integer /-1/ # Coarseness level at which the point is on grid
iocndlxy(ncondmax) _integer # Obsolete array, only used to recover old datasets
iocndlz(ncondmax) _integer # Obsolete array, only used to recover old datasets
checkconductors(nx:integer,ny:integer,nz:integer,nzfull:integer,dx:real,dy:real,dz:real,l2symtry:logical,l4symtry:logical) subroutine

*********** Multigrid3d dump:
mgparam    real    /1.2/ # Acceleration parameter for multigrid fieldsolver
mgmaxiters integer /100/ # Maximum number of iterations
mgiters    integer       # Actual number of iterations
mgtol      real  /1.e-6/ # Absolute tolerance in change in last iteration
mgmaxlevels integer /101/   # Minimum grid size in x-y to coarsen to
mgform     integer /1/   # When 1, MG operates on phi (and rho),
                         # when 2, MG operates on error (and residual)
downpasses integer /1/   # Number of downpasses
uppasses   integer /1/   # Number of uppasses
mglevels   integer /0/      # Number of coarsening levels
mglevelsnx(0:100) integer     # List of nx for the levels of coarsening
mglevelsny(0:100) integer     # List of ny for the levels of coarsening
mglevelsnzfull(0:100) integer # List of nzfull for the levels of coarsening
mglevelsiz(0:100) integer     # List of iz for the levels of coarsening
mglevelsnz(0:100) integer     # List of nz for the levels of coarsening
mglevelslx(0:100) real /101*1/ # List of coarsening factors in x
mglevelsly(0:100) real /101*1/ # List of coarsening factors in y
mglevelslz(0:100) real /101*1/ # List of coarsening factors in z
mglevelspart(0:100) logical    # List of flags for whether full or partial
                               # coarsening is done: 0 is full, 1 is partial
mggoodnumbers(40) integer /4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,80,96,112,
                           128,160,192,224,256,320,384,448,512,640,768,896,
                           1024,1280,1536,1792,2048,2560,3072,3584,5120,7168/
                         # A list of good numbers to use for the grid
                         # dimension. This is an ordered list of powers of two
                         # times 1, 3, 5, and 7.
tempsize   integer       # Size of work space (autoset)
phi_temp(tempsize) _real # Work space holding phi on all grid levels
rho_temp(tempsize) _real # Work space holding source on all grid levels
subgrid_sor_to_mg(nx:integer,ny:integer,nz:integer,dx:real,dy:real,dz:real,
                  l2symtry:logical,l4symtry:logical) subroutine
  # Converts a set of points generated for the SOR fieldsolver into the set of
  # points needed for the multigrid fieldsolver.
setmglevels(nx:integer,ny:integer,nz:integer,nzfull:integer,
            dx:real,dy:real,dz:real)
            subroutine   # Calculates levels of coarsening. Note that mglevels
                         # must be zero when calling this routine.

*********** Multigrid3d_work:
# Temporary variables and array used by subgrid_sor_to_mg
wnx integer
wny integer
wnz integer
iii(0:wnx,0:wny,0:wnz) _integer

*********** Surface_of_Rev dump:
srfrv_pernz            integer  /0/ # Number of points per nz for tablized data
srfrv_z                real # Value of z passed to srfrv_f
srfrv_r                real # Value of r returned by srfrv_f
lsrlinr                logical /.false./ # Use piecewise-linear curve.
npnts_sr               integer /0/ # Number points in piecewise-linear curve.
z_sr(npnts_sr)         _real # Z of points in piecewise-linear curve
r_sr(npnts_sr)         _real # R of points in piecewise-linear curve
lsrminlinr             logical /.false./ # Use piecewise-linear curve for rmin.
npnts_srmin            integer /0/ # Number points in piecewise-linear curve.
z_srmin(npnts_srmin)   _real # Z of points in piecewise-linear curve
r_srmin(npnts_srmin)   _real # R of points in piecewise-linear curve
lsrmaxlinr             logical /.false./ # Use piecewise-linear curve for rmax.
npnts_srmax            integer /0/ # Number points in piecewise-linear curve.
z_srmax(npnts_srmax)   _real # Z of points in piecewise-linear curve
r_srmax(npnts_srmax)   _real # R of points in piecewise-linear curve

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
       zgrid:real,linbends:logical,l2symtry:logical,l4symtry:logical,
       scrtch:real,izfsmin:integer,izfsmax:integer)
     subroutine # PSOR field solver
psorinit(nx,ny,nz,dx:real,dy:real,dz:real,l2symtry:logical,l4symtry:logical)
     subroutine # Initialize arrays that hold conductor points
cond_pot(nx,ny,nz,phi:real)
     subroutine # Sets potential in phi to desired potential on conductors
setcndtr(xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real,nx,ny,nz,
         dx:real,dy:real,dz:real,l2symtry:logical,l4symtry:logical)
     subroutine # Calculates conductor locations
rodpoint(ixmin,ixmax,iymin,iymax,dx:real,dy:real,rodx:real,rody:real,
         lz_in_rod:logical,rminsq:real,rodrrsq:real,ix_axis,iy_axis,rrr:real,
         rodmxsq:real,delz_in:real,quadsl:real,rodap:real,fuzz:real,qzs:real,
         dz:real,iz,voltage:real)
     subroutine # Set conductor points for general rod.
platepnt(ixmin:integer,ixmax:integer,iymin:integer,iymax:integer,
         ix_axis:integer,iy_axis:integer,dx:real,dy:real,aper:real,rmax:real,
         vvv:real,xoff:real,yoff:real,delz_in:real,iz:integer,
         lz_in_plate:logical,fuzz:real)
     subroutine # Set conductor points for a transverse plate
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
         l2symtry:logical,l4symtry:logical)
     subroutine #  The 3d Poisson solver
vp3x(iwhich) subroutine
     # BASIS-level interface to VPOIS3d, using FS3 database variables
     #  The user program should declare a similar subroutine w/ its vars.
vpois2d  (iwhich, a:real, ak:real, kxsq:real, kysq:real, attx:real, atty:real,
          filt:real, lx:real, ly:real, nx, ny, work:real, xywork:real, ibc,
          l2symtry:logical,l4symtry:logical)
     subroutine #  The 2d Sine-Sine Poisson solver
vsftx    (a:real, work:real, cp:real, cm:real, nx, ny, isetup)
     subroutine #  Vectorized Sine Fourier Transform in X
vsfty    (a:real, work:real, cp:real, cm:real, nx, ny, isetup)
     subroutine #  Vectorized Sine Fourier Transform in Y
cosqx(a:real,w:real,c:real,nx:integer,ny:integer,isign:integer) subroutine
cosqy(a:real,w:real,c:real,nx:integer,ny:integer,isign:integer) subroutine
vcpft(r:real,i:real,n:integer,incp:integer,signp:integer,lenv:integer,
      lfd:integer) subroutine
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
setconductorparity(nn:integer,ix:integer,iy:integer,iz:integer,
                   dels:real,parity:integer,fuzz:real,fuzzsign:integer,
                   dfill:real) subroutine
zplaneconductorf(zcent:real,zsign,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
boxconductorf(xsize:real,ysize:real,zsize:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
cylinderconductorf(rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
cylindersconductorf(ncylinders:integer,rad:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
xcylinderconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
ycylinderconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zcylinderconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zroundedcylinderconductorf(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zcylinderoutconductorf(rad:real,length:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zroundedcylinderoutconductorf(rad:real,length:real,rad2:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
sphereconductorf(rad:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
coneconductorf(r_zmin:real,r_zmax:real,length:real,theta:real,phi:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
conesconductorf(ncones:integer,r_zmin:real,r_zmax:real,length:real,
        theta:real,phi:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zconeconductorf(r_zmin:real,r_zmax:real,length:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
zconeoutconductorf(r_zmin:real,r_zmax:real,length:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
ztorusconductorf(r1:real,r2:real,xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine
beamletplateconductorf(za:real,zb:real,z0:real,thickness:real,
        xcent:real,ycent:real,zcent:real,
        n:integer,x:real,y:real,z:real,delmx:real,delpx:real,
        delmy:real,delpy:real,delmz:real,delpz:real,fuzz:real) subroutine

******** ConductorGeometryVisualization:
maxtriangles integer/0/
ntriangles integer /0/
triangles(0:2,0:2,maxtriangles) _real
normals(0:2,0:2,maxtriangles) _real
getconductorfacets(nc:integer,icnd:integer,dels:real,
                   gridn:integer,griddd:real,gridmin:real) subroutine
getconductorsnewfacet(ix:integer,iy:integer,iz:integer,oo:integer,
                      parity:integer,gridn:integer,iii:integer,nc:integer,
                      dels:real,gridmin:real,griddd:real,pp:real,npp:integer)
                      subroutine
conductorsmoothshading() subroutine
