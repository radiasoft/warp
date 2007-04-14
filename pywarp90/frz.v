frz
#@(#) File FRZ.V, version $Revision: 3.57 $, $Date: 2007/04/14 00:08:04 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package FRZ of code WARP6
# FRZ - fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
# Debbie Callahan, LLNL, (510)423-5926
{
}

*********** FRZversion:
versfrz character*19 /"$Revision: 3.57 $"/#  Code version set by CVS

*********** FRZvars:
# Variables needed by the test driver of package FRZ
ibc                       integer /0/  #  Boundary conditions-future use
lr                        real /0.7/   #  System length in r (arbitrary units)
lz                        real /1.6/   #  System length in z (arbitrary units)
eta                       real /0.0/   #  Resistivity of wall (ohm/m)
taurc                     real /0.0/   #  RC time
filt(5)                   real /5*0./  #  Spatial filtering coefficients
vbeam                     real /0./    #  Beam velocity
nr                        integer /64/ #  Mesh points are 0,...,nr
nz                        integer /64/ #  Mesh points are 0,...,nz
b(0:nr,0:nz)              _real        #  Charge density, potential array
bsav(0:nr,0:nz)           _real        #  "Save" array for b
schrg(0:nz)               _real        #  Surface charge for resistive wall
attz(0:nz/2)              _real        #  Attenuation factor as fcn. of kz
kzsq(0:nz)                _real        #  Discrete analog to dr*dr*kz^2
r(0:nr)                   _real        #  radial distance
rfsmat(0:nr,3,0:nz)       _real        #  Tridi solve matrix
scrtch(0:nz)              _real        #  workspace for resistive wall
scrtch2(0:nz)             _real        #  workspace for resistive wall
phikold(0:nz)             _real        #  FT of phi at old time step for C
err1(0:nr,0:nz)           _real

*********** FRZmgrid dump:
mgridrz_accuracy          real /1.e-3/  # average accuracy of multigrid solver
mgridrz_ncmax             integer /100/ # maximum number of full multigrid 
                                        # cycles
mgridrz_npre              integer /4/   # number of relaxations steps before 
                                        # coarsening, in multigrid solver
mgridrz_npost             integer /4/   # number of relaxations steps after 
                                        # coarsening, in multigrid solver  
mgridrz_ncycles           integer /2/   # number of multigrid cycles per level
mgridrz_nlevels_max       integer /100/ # maximum number of multigrid levels
mgridrz_levels_min        integer /1/   # lowest level of coarsening
mgridrz_nmeshmin          integer /8/   # minimum number of meshes in each direction at coarsest level
mgridrz_mgparam           real /1.8/    # SOR parameter
mgridrz_workfact          integer /4/   # weight factor for grid merging/procs 
mgridrz_mgiters           real /0/      # actual number of iterations for a solve
mgridrz_sub_accuracy      real /1.e-1/  # average accuracy for a sublevel
mgridrz_deform            logical /.false./ # flag for use of elliptic deformation
mgridrz_nx                integer  /0/  # 
mgridrz_ny                integer  /0/  # 
mgridrz_nz                integer  /0/  # 
mgridrz_xfact(0:mgridrz_nz+2) _real     # array for deformation factor in X
mgridrz_yfact(0:mgridrz_nz+2) _real     # array for deformation factor in Y
mgridrz_phi3d(0:mgridrz_nx,0:mgridrz_ny,-1:mgridrz_nz) _real # array containing '3D' phi
mgridrz_ngrids            integer  /1/  # number of grids (useful when using mesh refinement)
mgridrz_grid_is(mgridrz_ngrids)   _integer # array id grid associated with grid species 
ngrids                    integer  /1/  # number of grids (includes base grid + patches)
nrg(ngrids)               _integer      # number of mesh in R for each grid
nzg(ngrids)               _integer      # number of mesh in Z for each grid
drg(ngrids)               _real         # number of mesh in R for each grid
dzg(ngrids)               _real         # number of mesh in Z for each grid
lverbose                  integer  /1/  # level of verbosity (0=no output; 1=low; 2=medium; 3=high)
l_change_grid             logical  /.false./ # change grid patches during step if true
ngrids_cg                 integer  /0/       # nb grid patches to change
id_cg(ngrids_cg,2)        _integer           # IDs grid patches to change
nr_cg(ngrids_cg)          _integer           # new nr grid patch change
nz_cg(ngrids_cg)          _integer           # new nz grid patch change
dr_cg(ngrids_cg)          _real              # new dr grid patch change
dz_cg(ngrids_cg)          _real              # new dz grid patch change
rmin_cg(ngrids_cg)        _real              # new rmin grid patch change
zmin_cg(ngrids_cg)        _real              # new zmin grid patch change
transit_min_r_cg(ngrids_cg) _integer  # new min guard in r grid patch change (helps reduce spurious force)
transit_max_r_cg(ngrids_cg) _integer  # new max guard in r grid patch change (helps reduce spurious force)
transit_min_z_cg(ngrids_cg) _integer  # new min guard in z grid patch change (helps reduce spurious force)
transit_max_z_cg(ngrids_cg) _integer  # new max guard in z grid patch change (helps reduce spurious force)
l_change_loc_part         logical  /.false./ # if true, loc_part is changed according to array rmc at next step
nz_rmc                    integer            # size array rmc
rmc(nz_rmc+1)             _integer           # minimum radius for field gathering  
l_get_field_from_base     logical  /.false./ # if true, gather field from base grid only
l_get_injphi_from_base    logical  /.false./ # if true, gather injphi from base grid only
l_dep_rho_on_base         logical  /.false./ # if true, deposit rho on base grid only
l_distribute              logical  /.true./  # if true, distribute rho between high level to low level patches
nguardx                   integer /1/        # number of guard cell in x/r
nguardz                   integer /1/        # number of guard cell in z
grids_nids                integer            # number of grid IDs
n_avail_ids               integer            # number of available IDs
avail_ids(100)            _integer           # array of grid IDs
level_del_grid            integer            # level of deleted grid

*********** multigrid_common_base dump:
l_mgridrz_debug  logical  /.false./
l_timing_rz      logical  /.true./
restrictwbnd     logical  /.true./
vlocs            logical  /.false./
l_print_timing   logical  /.false./
l_jump           logical  /.false./
t_relax          real
t_restrict       real
t_expand         real
t_before         real
t_apply_voltage  real
t_updateguard    real
t_allocate       real
inveps0          real
o_line           character*(240) # line for output
bndy_allocated   logical  /.false./ # flag set to true if bndy is allocated
ixlbndi          integer  /0/
ixrbndi          integer  /0/
izlbndi          integer  /0/
izrbndi          integer  /0/
ixlbnd           integer  /0/
ixrbnd           integer  /0/
izlbnd           integer  /0/
izrbnd           integer  /0/
bnd_method       integer  /0/
nlevels          integer      # number of multigrid levels
level            integer      # current multigrid level
nb_iters         integer      # actual number of iterations used for a solve

*********** FRZmgrid_ptrs dump:
basegrid _GRIDtype  # primary grid for RZ solver

*********** BWorkRZ dump:
bworkgrid _GRIDtype  # Working grid for the RZ B field solver

*********** InjectVars_eq dump:
# variables and functions needed for getting voltage risetime from assumption of 
# constant injection (works with inj_d=2).
inj_phi_eq           real # Electrostatic potential at the emitting surface at equilibrium.
v_max                real /0./
l_find_rise_time     logical /.false./
afact                real /1./
calc_a               integer  /1/  # determines way of calculating voltage factor for rise time
init_gridinit() subroutine #

*********** FRZsubs:
#  Callable subroutines in the FRZ package
vpoisrz  (iwhich, a:real, ak:real, kzsq:real,schrg:real, eta:real,
          phikold:real, taurc:real,
         attz:real, filt:real, dt:real, vbeam:real,
         lr:real, lz:real, nr, nz, rfsmat:real, 
         scrtch:real, scrtch2:real, ibc)   integer function
         #  The RZ Poisson solver
vprzx     (iwhich,dt:real)                                          subroutine
         #  Python-level interface to VPOISRZ, using FRZ database variables
         #  The user program should declare a similar subroutine w/ its vars.
advsc    (schrg:real,eta:real,a:real,kzsq:real,dt:real,dr:real,dz:real,
          nr,nz,phikold:real,taurc:real)       subroutine
         #  Routine that advances the surface charge from time t
         #  to time t+dt in Fourier space
advect   (schrg:real, vbeam:real, dt:real, nz, dz, tmp:real)        subroutine
         #  Advects surface charge with the moving window.
multigridrzf(iwhich:integer,phi:real,rho:real,nx:integer,nz:integer) subroutine
         # Multigrid Poisson solver (using "full-multigrid" method, i.e. the 
         # solution is calculatedd at each level using a multigrid procedure 
         # and used as an approximated solution to start the calculation at 
         # the next level)
solve_mgridrz(grid:GRIDtype,accuracy:real,fromup:logical) subroutine
         # Deeper interface to the RZ multigrid solver. This is called by
         # multigridrzf.
setmglevels_rz(grid:GRIDtype) subroutine
         # set mglevels in f3d arrays from RZ solver structure
get_cond_rz(grid:integer) subroutine
         # get internal conductors locations from RZ multigrid solver
get_cond_rz_level(grid:integer,level:integer) subroutine
         # get internal conductors locations from RZ multigrid solver
setconductorvoltagerz(volt:real,nz:integer,zmmin:real,dz:real,discrete:logical,
                      condid:integer)
         subroutine
         # set voltage on conductors from a z-grid
setconductorvoltagerz_id(id:integer,volt:real) subroutine
         # set voltage on conductor given its ID
cond_sumrhointerior2d(rhosum:real,grid:GRIDtype,nx:integer,nz:integer,
                      rho(0:nx,0:nz):real,ixmin:integer,ixmax:integer,
                      izmin:integer,izmax:integer,dr:real,rmmin:real) subroutine
         # Sums rho in the interior of the conductor
calcfact_deform(dz:real,zmin:real,
                xfact:real,yfact:real,nz:integer,ns:integer,is:integer,
                ins:integer,nps:integer,ws:real,zgrid:real) subroutine
         # computes factors for elliptical deformation in X and Y planes
init_base(nr:integer,nz:integer,dr:real,dz:real,rmin:real,zmin:real) subroutine
         # initializes the base grid for RZ solver
del_base() subroutine
         # removes the base grid
set_basegrid() subroutine
nullify_basegrid() subroutine
mk_grids_ptr() subroutine
add_subgrid(id:integer,nr:integer,nz:integer,dr:real,dz:real,
            rmin:real,zmin:real,
            transit_min_r:integer,transit_max_r:integer,
            transit_min_z:integer,transit_max_z:integer) subroutine
         # add a subgrid to the grid id
del_subgrid(id:integer) subroutine
         # delete a subgrid and all its 'children'
init_gridbnd(g:GRIDtype) subroutine
         # Initialize the boundary instance
del_conductors() subroutine
         # delete all conductors data on RZ grids
get_phi_subgrid(id:integer,phi:real,nr:integer,nz:integer) subroutine
         # get the potential of grid id
set_rho_rz(rho:real,nr:integer,nz:integer,id:integer) subroutine
         # set rho of grid id       
mix_rho_rz(rho:real,nr:integer,nz:integer,id:integer,fmix:real) subroutine
         # set rho of grid id       
get_rho_rz(rho:real,nr:integer,nz:integer,id:integer,rhop:integer) subroutine
         # get rho of grid id
reset_rzmgrid_rho() subroutine
         # sets rho to zero.
rhoweightrz(xp:real,yp:real,zp:real,np:integer,q:real,nr:integer,nz:integer,
            dr:real,dz:real,rgrid:real,zgrid:real) subroutine
rhoweightrzgrid(grid:GRIDtype,xp(np):real,yp(np):real,zp(np):real,np:integer,
                q:real,nr:integer,nz:integer,dr:real,dz:real,
                rgrid:real,zgrid:real) subroutine
         # deposits rho on the specified grid object
rhoweightrzgrid_weights(grid:GRIDtype,xp(np):real,yp(np):real,zp(np):real,
                        w(np):real,np:integer,
                        q:real,nr:integer,nz:integer,dr:real,dz:real,
                        rgrid:real,zgrid:real) subroutine
         # deposits rho from weighted particles on the specified grid object
fieldweightz(zp:real,ez:real,np:integer,zgrid:real) subroutine
fieldweightrz(xp:real,yp:real,zp:real,ex:real,ey:real,ez:real,np:integer,zgrid:real,efetch:integer) subroutine
fieldweightxz(xp:real,zp:real,ex:real,ez:real,np:integer,zgrid:real,efetch:integer) subroutine
dep_rho_rz(is:integer,rho:real,nr:integer,nz:integer,dr:real,dz:real,
           xmin:real,zmin:real) subroutine
         # makes rho deposition on RZ grid
distribute_rho_rz() subroutine
         # recursively distributes rho from fine patches to coarse patches
find_mgparam_rz(lsavephi:logical) subroutine
         # RZ version of find_mgparam. Does the search for each subgrid.
find_mgparam_rz_1g(grid:GRIDtype) subroutine
         # RZ version of find_mgparam for one grid only.
gchange_rhop_phip_rz() subroutine
         # reallocate rhop and phip arrays
install_conductors_rz(conductors:ConductorType,grid:GRIDtype) subroutine
         # install conductors data into RZ arrays
set_basegrid_phi() subroutine
	 # set phi on basegrid using w3d.phi 
setbnd_subgrid_to_inj_d() subroutine
         # set indices for force gathering (locpart) to coarser grid 
         # when grid point more than inj_d*grid%dz away from emitting surface
clean_conductor_interior() subroutine
         # removes conductor interior points from calculation 
build_vlocs() subroutine
set_patches_around_emitter(id:integer,np:integer,ij:integer,nz:integer,
            transit_min_r:integer,transit_max_r:integer,
            transit_min_z:integer,transit_max_z:integer) subroutine
         # add patches around emitter
test_subgrid_rz() subroutine
test_subgrid_xz() subroutine
sum_neighbors(fin:integer,fout:integer,nx:integer,ny:integer) subroutine
         # returns sum of neighboring cells (inculding itself)
adjust_lpfd(f:integer,nr:integer,nz:integer,rmin:real,rmax:real,zmin:real,zmax:real) subroutine
         # adjust loc_part_fd arrays according to array (f) used to generate grids.
setphirz(np:integer,xp:real,yp:real,zp:real,p:real,zgrid:real) subroutine
         # get phi in p from RZ grid at locations [x,y,z]

setphixz(np:integer,xp:real,yp:real,zp:real,p:real,zgrid:real) subroutine
         # get phi in p from XZ grid at locations [x,y,z]

init_bworkgrid(nr:integer,nz:integer,dr:real,dz:real,rmin:real,zmin:real,
               bounds(0:5):integer,
               u(0:nr+2,0:nz+2):real,rho(nr+1,nz+1):real) subroutine
         # Sets up the GRDItype object used during the RZ magnetic field solve.
multigridrzb(iwhich:integer,iaxis:integer,u0(0:nr0+2,0:nz0+2):real,
             rho0(nr0+1,nz0+1):real,nr0:integer,nz0:integer,accuracy:real)
         subroutine
         # Does one part of the RZ magnetic field solve
lphiberz(nx:integer,nz:integer,nzfull:integer,dxsqi:real,dzsqi:real,
         phi:real,res:real,mglevel:integer,bounds:integer,mgparam:integer,
         lcndbndy:logical,icndbndy:integer,conductors:ConductorType,
         iondensity:real,electrontemperature:real,plasmapotential:real,
         electrondensitymaxscale:real,epsilon:real,epszfacdzsqi:real,
         inormtoavboltzfac:logical) subroutine
restrictberz(nx:integer,nz:integer,nzfull:integer,u:real,delx:integer,delz:integer,
             nxcoarse:integer,nzcoarse:integer,nzfullcoarse:integer,ucoarse:real,
             delcx:integer,delcz:integer,
             bounds:integer,boundscoarse:integer,lzoffset:integer,
             withoutzeros:logical) subroutine
getanalyticbtheta(bfield:BFieldGridType) subroutine
         # Performs and optional ananlytic calculation of Btheta
updateguardcells2d() subroutine
         # update guard cells of 2D multigrid solver

%%%%%%%% CONDtype:
# structure for potential calculation close to conductors.
# The stencil for the iterative calculation of the potential f is given by
#     i = 1 -> nbbnd
#        j=jj(i); l=kk(i)
#        f(j,l) = dt*(cf0(i)   * f(j  ,l  )
#               +     cfxp(i)  * f(j+1,l  )
#               +     cfxm(i)  * f(j-1,l  )
#               +     cfzp(i)  * f(j  ,l+1)
#               +     cfzm(i)  * f(j  ,l-1)
#               +     phi0xp(i)
#               +     phi0xm(i)
#               +     phi0zp(i)
#               +     phi0zm(i)
#               +     rhs(j,l))
# for the ith grid point being close to the considered conductor and located
# at (j,l) on the grid. If there is a conductor boundary lying for example
# between j and j+1, cfxp(i) is set to zero while phi0xp(i) is set to
# the voltage of the nearest conductor.
ncond            integer  # number of nodes inside conductor
nbbnd            integer  # number of "red" nodes near conductor (for red-black gauss-seidel)
nbbndred         integer  # number of nodes inside conductor
voltage(ncond)  _real     # conductors voltage
condid(ncond)   _integer  # conductors ID
cf0(nbbnd)      _real     # stencil coefficients for relaxation iteration at node  
cfxm(nbbnd)     _real     # stencil coefficients for relaxation iteration at j-1
cfxp(nbbnd)     _real     # stencil coefficients for relaxation iteration at j+1
cfzm(nbbnd)     _real     # stencil coefficients for relaxation iteration at l-1
cfzp(nbbnd)     _real     # stencil coefficients for relaxation iteration at l+1
dt(nbbnd)       _real     # overall coefficient for relaxation
phi0(nbbnd)     _real     # phi in conductor at j,l
phi0xm(nbbnd)   _real     # phi in conductor at j-1 (phi0xm=cfxm*volt0xm)
phi0xp(nbbnd)   _real     # phi in conductor at j+1 (phi0xp=cfxp*volt0xp)
phi0zm(nbbnd)   _real     # phi in conductor at l-1 (phi0zm=cfzm*volt0zm)
phi0zp(nbbnd)   _real     # phi in conductor at l+1 (phi0zp=cfzp*volt0zp)
volt0xm(nbbnd)  _real     # voltage in conductor at j-1
volt0xp(nbbnd)  _real     # voltage in conductor at j+1
volt0zm(nbbnd)  _real     # voltage in conductor at l-1
volt0zp(nbbnd)  _real     # voltage in conductor at l+1
condidxm(nbbnd) _integer  # conductor ID at j-1
condidxp(nbbnd) _integer  # conductor ID at j+1
condidzm(nbbnd) _integer  # conductor ID at l-1
condidzp(nbbnd) _integer  # conductor ID at l+1
dxm(nbbnd)      _real     # distance from node to conductor at j-1
dxp(nbbnd)      _real     # distance from node to conductor at j+1
dzm(nbbnd)      _real     # distance from node to conductor at l-1
dzp(nbbnd)      _real     # distance from node to conductor at l+1
jj(nbbnd)       _integer  # location node in x points near conductor
kk(nbbnd)       _integer  # location node in z points near conductor
docalc(nbbnd)   _logical  # performs calculation if set to .true.
jcond(ncond)    _integer  # location node in x points in conductor
kcond(ncond)    _integer  # location node in z points in conductor
next            _CONDtype # next conductor element in linked list
prev            _CONDtype # previous conductor element in linked list

%%%%%%%% BNDtype:
nb_conductors     integer
nr                integer 
nz                integer
nvlocs            integer
nvlocsred         integer
izlbnd            integer
izrbnd            integer
nworkpproc        integer
l_powerof2        logical
l_merged          logical
dr                real
dz                real
zmin              real
zmax              real
v(1:nr+1,1:nz+1)  _integer  
vlocs_j(nvlocs)   _integer
vlocs_k(nvlocs)   _integer
cndfirst          _CONDtype
cndlast           _CONDtype
next              _BNDtype
prev              _BNDtype

%%%%%%%% OVERLAPtype:
gid(1)               _integer +fassign
rmin                  real
rmax                  real
zmin                  real
zmax                  real
index                 integer
nr                    integer
nz                    integer
rho(1:nr+1,1:nz+1,2) _real +fassign
next                 _OVERLAPtype

%%%%%%%% GRIDtype:
nguardx                   integer 
nguardz                   integer  
gid(1) _integer 
nlevels integer
nr integer
nz integer
nzp integer
nrpar integer
nzpar integer
jmin integer
jmax integer
lmin integer
lmax integer
ixlbnd integer
ixrbnd integer
izlbnd integer
izrbnd integer
rmin real
rmax real
xmin real
xmax real
zmin real
zmax real
dr real
dz real
invdr real
invdz real
zminp real
invvol(1:nr+1) _real
rho(1:nr+1,1:nz+1) _real
phi(1-nguardx:nr+nguardx+1,1-nguardz:nz+nguardz+1) _real        # potential
rhop(1:nrpar+1,1:nzpar+1) _real
phip(1-nguardx:nrpar+nguardx+1,1-nguardz:nzpar+nguardz+1) _real #
rhominr integer
rhomaxr integer
rhominz integer
rhomaxz integer
loc_part(1:nr+1,1:nzp+1) _integer
loc_part_fd(1:nr+1,1:nzp+1) _integer
npre integer
npost integer
ncycles integer
ncmax integer
npmin integer
transit_min_r integer
transit_max_r integer
transit_min_z integer
transit_max_z integer
mgparam real
lmagnetostatic logical /.false./ # When true, includes extra terms in the Poisson equations for Ar and Atheta
bndfirst _BNDtype
bndlast  _BNDtype
next _GRIDtype
prev _GRIDtype
down _GRIDtype
up   _GRIDtype
neighbors _OVERLAPtype
parents   _OVERLAPtype
children  _OVERLAPtype

%%%%%%%% GRDPTRtype:
grid _GRIDtype

$******** SuperLUInterface:
$superlu_dgssv(n:integer,nnz:integer,nrhs:integer,
$              values:real,rowind:integer,colptr:integer,
$              b:real,info:integer)
$             subroutine

******** ImplicitMG2D:
coeffs1(:,:,:,:) _real
chi01(:,:,:) _real
mgsolveimplicites2d(iwhich:integer,nx:integer,nz:integer,nzfull:integer,
                    dx:real,dz:real,phi:real,rho:real,ns:integer,qomdt:real,chi0:real,bounds:integer,
                    xmmin:real,zmmin:real,zmminglobal:real,zbeam:real,zgrid:real,
                    mgparam:real,mgiters:integer,mgmaxiters:integer,
                    mgmaxlevels:integer,mgerror:real,mgtol:real,mgverbose:integer,
                    downpasses:integer,uppasses:integer,
                    lcndbndy:logical,laddconductor:logical,icndbndy:integer,lbuildquads:logical,
                    gridmode:integer,conductors:ConductorType,lrz:logical,
                    my_index:integer,nslaves:integer,izfsslave:integer,nzfsslave:integer)
            subroutine
residuales2d(nx:integer,nz:integer,nzfull:integer,phi:real,rho:real,coeffs:real,res:real,
             mglevel:integer,bounds(0:5):integer,mgparam:real)
            subroutine

******** Subtimersfrz:
lfrztimesubs logical /.false./
timemultigridberzsolve   real /0./
timemgsolveimplicites2d real /0./

