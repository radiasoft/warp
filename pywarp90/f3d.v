f3d
#@(#) File F3D.V, version $Revision: 3.11 $, $Date: 2001/08/31 22:36:40 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package F3D of code WARP6
# fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
{
}

*********** F3Dversion:
versf3d character*19 /"$Revision: 3.11 $"/#  Code version version is set by CVS

*********** F3Dvars:
# Variables needed by the test driver of package F3D
lx                        real /0.7/   #  System length in x (arbitrary units)
ly                        real /1.3/   #  System length in y (arbitrary units)
lz                        real /1.6/   #  System length in z (arbitrary units)
filt(5,3)                 real /15*0./ #  Spatial filtering coefficients
ibc                  character*8 /" "/ #  Boundary conditions-future use
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
nsorerr   integer -dump  # number of points used for error calculation
sorerrar(nsorerr) _real -dump # Error array
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
dxpsor    real         # Size of transverse grid cells.
nxpsor    integer /0/  # Maximum number of transverse grid points
nzpsor    integer      # Number of z grid points
boundarr(0:nxpsor,4,0:nzpsor) _real -dump # Array used to preserve boundaries
                                          # during iterations
ncondmax          integer # Maximum number of points in conductor
ncond             integer # Number of points within conductors
ixcond(ncondmax) _integer # X coordinate of points in conductor
iycond(ncondmax) _integer # Y coordinate of points in conductor
izcond(ncondmax) _integer # Y coordinate of points in conductor
condvolt(ncondmax) _real  # voltage of points in conductor
lplates          logical  # Sets whether or not quadruple endplates are included
rodfract        real /1./ # Fraction of quadrupole rod which is used
islctype(-1:nzpsor+1) _integer -dump # Type of structure at each z slice,
                          # 0: drift, 1:quads on x axis, 2: quads on both axis,
                          # 3: quads on y axis, 4: plate left of quad, 5: plate
                          # right of quad
slicvx(-1:nzpsor+1) _real -dump # Voltage of quads on x axis
slicvy(-1:nzpsor+1) _real -dump # Voltage of quads on y axis

lcndbndy logical /.true./ # Turns on sub-grid boundaries
icndbndy integer /1/      # Type of interpolant to use for sub-grid boundaries
                          # 1 egun style
                          # 2 EBC style (non-centered finite-difference)
ncndmax       integer /0/ # Maximum number of points for sub-grid boundaries
necndbdy          integer # Number of points for even sub-grid boundaries
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
nocndbdy        integer # Number of points for odd sub-grid boundaries
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

*********** MultigridConductor3d dump parallel:
ecvoltmx(ncndmax)     _real # Voltage on conductor in minus x direction
ecvoltpx(ncndmax)     _real # Voltage on conductor in plus  x direction
ecvoltmy(ncndmax)     _real # Voltage on conductor in minus y direction
ecvoltpy(ncndmax)     _real # Voltage on conductor in plus  y direction
ecvoltmz(ncndmax)     _real # Voltage on conductor in minus z direction
ecvoltpz(ncndmax)     _real # Voltage on conductor in plus  z direction
ocvoltmx(ncndmax)     _real # Voltage on conductor in minus x direction
ocvoltpx(ncndmax)     _real # Voltage on conductor in plus  x direction
ocvoltmy(ncndmax)     _real # Voltage on conductor in minus y direction
ocvoltpy(ncndmax)     _real # Voltage on conductor in plus  y direction
ocvoltmz(ncndmax)     _real # Voltage on conductor in minus z direction
ocvoltpz(ncndmax)     _real # Voltage on conductor in plus  z direction
icondlxy(ncondmax) _integer # Coarseness level at which the point is on grid
icondlz (ncondmax) _integer # Coarseness level at which the point is on grid
iecndlxy(ncndmax)  _integer # Coarseness level at which the point is on grid
iecndlz (ncndmax)  _integer # Coarseness level at which the point is on grid
iocndlxy(ncndmax)  _integer # Coarseness level at which the point is on grid
iocndlz (ncndmax)  _integer # Coarseness level at which the point is on grid

*********** Multigrid3d dump:
mgparam    real    /1.2/ # Acceleration parameter for multigrid fieldsolver
mgmaxiters integer /100/ # Maximum number of iterations
mgiters    integer       # Actual number of iterations
mgtol      real  /1.e-6/ # Absolute tolerance in change in last iteration
downpasses integer /1/   # Number of downpasses
uppasses   integer /1/   # Number of uppasses
tempsize   integer       # Size of work space (autoset)
mggoodnumbers(40) integer /4,6,8,10,12,14,16,20,24,28,32,40,48,56,64,80,96,112,
                           128,160,192,224,256,320,384,448,512,640,768,896,
                           1024,1280,1536,1792,2048,2560,3072,3584,5120,7168/
                         # A list of good numbers to use for the grid
                         # dimension. This is an ordered list of powers of two
                         # times 1, 3, 5, and 7.
phi_temp(tempsize) _real # Work space holding phi on all grid levels
rho_temp(tempsize) _real # Work space holding source on all grid levels
conductor_data_level(nx:integer,ny:integer,nz:integer,dx:real,dy:real,dz:real)
  subroutine
  # Calculates level of coarseness at which all of the conductor points
  # are on the grid.
subgrid_sor_to_mg(nx:integer,ny:integer,nz:integer,dx:real,dy:real,dz:real,
                  l2symtry:logical,l4symtry:logical) subroutine
  # Converts a set of points generated for the SOR fieldsolver into the set of
  # points needed for the multigrid fieldsolver.

*********** Multigrid3d_res:
# Work array which holds the residual.
res_size integer
res(res_size) _real

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
srfrvout(rofzfunc:string,volt:real,zmin:real,zmax:real,
         xcent:real,ycent:real,rmax:real,lfill:logical,
         xmin:real,xmax:real,ymin:real,ymax:real,lshell:logical,
         zmmin:real,zmmax:real,zbeam:real,
         dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,
         ix_axis:integer,iy_axis:integer,xmesh:real,ymesh:real,
         l2symtry:logical,l4symtry:logical)
     subroutine # Set conductor points for a conductor that is the
                # outside of a surface of revolution.
srfrvin (rofzfunc:string,volt:real,zmin:real,zmax:real,
         xcent:real,ycent:real,rmin:real,lfill:logical,
         xmin:real,xmax:real,ymin:real,ymax:real,lshell:logical,
         zmmin:real,zmmax:real,zbeam:real,
         dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,
         ix_axis:integer,iy_axis:integer,xmesh:real,ymesh:real,
         l2symtry:logical,l4symtry:logical)
     subroutine # Set conductor points for a conductor that is the
                # inside of a surface of revolution.
srfrvinout(rminofz:string,rmaxofz:string,volt:real,zmin:real,zmax:real,
         xcent:real,ycent:real,lzend:logical,
         xmin:real,xmax:real,ymin:real,ymax:real,lshell:logical,
         zmmin:real,zmmax:real,zbeam:real,
         dx:real,dy:real,dz:real,nx:integer,ny:integer,nz:integer,
         ix_axis:integer,iy_axis:integer,xmesh:real,ymesh:real,
         l2symtry:logical,l4symtry:logical)
     subroutine # Set conductor points for a conductor that is between two
                # surfaces of revolution.

*********** F3Dsubs:
#  Callable subroutines in the F3D package
vpois3d  (iwhich,a:real,ak:real,kxsq:real,kysq:real,kzsq:real,
         attx:real,atty:real,attz:real,filt:real,
         lx:real,ly:real,lz:real,nx,ny,nz,work:real,ibc,
         l2symtry:logical,l4symtry:logical)
     subroutine #  The 3d Poisson solver
vp3x(iwhich) subroutine
     # BASIS-level interface to VPOIS3d, using FS3 database variables
     #  The user program should declare a similar subroutine w/ its vars.
vpois2d  (iwhich, a:real, ak:real, kxsq:real, kysq:real, attx:real, atty:real,
          filt:real, lx:real, ly:real, nx, ny, work:real, ibc,
          l2symtry:logical,l4symtry:logical)
     subroutine #  The 2d Sine-Sine Poisson solver
vp2x     (iwhich) subroutine
     #  BASIS-level interface to VPOIS2d, using FS3 database variables
     #  The user program should declare a similar subroutine w/ its vars.
vsftx    (a:real, work:real, cp:real, cm:real, nx, ny, isetup)
     subroutine #  Vectorized Sine Fourier Transform in X
vsfty    (a:real, work:real, cp:real, cm:real, nx, ny, isetup)
     subroutine #  Vectorized Sine Fourier Transform in Y
pipe3df  (iwhich, pipeshpe:real, rho:real, phi:real, kxsq:real, kysq:real,
          kzsq:real, attx:real, atty:real, attz:real, filt:real,
          xlen:real, ylen:real, zlen:real, nx, ny, nz, scrtch:real,
          l2symtry:logical,l4symtry:logical)
     subroutine #  External interface to capacity field solver
pipest3d (pipeshpe:real, cap:real, phi:real, kxsq:real, kysq:real, kzsq:real,
          attx:real, atty:real, attz:real, filt:real, xlen:real,
          ylen:real, zlen:real, nx, ny, nz, scrtch:real, piperadi:real,
          pipen:real,pipe8th:real,pipex:real,pipey:real,cap3d:real,kpvt,
          l2symtry:logical,l4symtry:logical)
     subroutine #  Sets up capacity matrix
setcap3d (ncndpts,quadx:real,quady:real,quadz:real,quadcap:real,kpvt,nx,ny,nz,
          phi:real,work:real,kxsq:real,kysq:real,kzsq:real,attx:real,atty:real,
          attz:real,xlen:real,ylen:real,zlen:real,filt:real,
          l2symtry:logical,l4symtry:logical)
     subroutine # Sets up 3-d cap matrix
capfs    (iwhich,nx,ny,nz,phi:real,rho:real,xlen:real,ylen:real,zlen:real,
          kxsq:real,kysq:real,kzsq:real,attx:real,atty:real,attz:real,
          filt:real,work:real,workz:real,xmmax:real,
          pipeshpe:string,periinz:logical,l2symtry:logical,l4symtry:logical)
     subroutine # 3-d field solve using cap matrices
cndctr3d (nx,ny,xmmax:real,dx:real,dy:real,quadradi:real,quadcent:real,ncndpts,
          quadx:real,quady:real,quadz:real,quadv:real,nzpts,nendquad,nzquad,
          qdslclen:real,pipeshpe:string,loadquad:logical)
     subroutine # Finds points on quadrupole conductor surfaces
findqdnd(nquad:integer,quadzs:real,quadde:real,quadvx:real,quadvy:real,
         zmmin:real,zgrid:real,nz:integer,dz:real,
         numends:integer,nendmax:integer,quadend:real,quadvlt:real,
         vshift:real,quadcent:real,quadradi:real)
     subroutine # Finds starts of quads that are on main grid
vcap3d   (iwhich,rho:real,phi:real,kxsq:real,kysq:real,kzsq:real,attx:real,
          atty:real,attz:real,filt:real,xlen:real,ylen:real,zlen:real,
          nx,ny,nz,scrtch:real,xmmax:real,zmmin:real,zgrid:real,
          pipeshpe:string,periinz:logical,l2symtry:logical,l4symtry:logical)
     subroutine # External routine for capacity matrix field solve
