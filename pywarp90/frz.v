frz
#@(#) File FRZ.V, version $Revision: 3.21 $, $Date: 2002/08/29 17:49:35 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package FRZ of code WARP6
# FRZ - fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
# Debbie Callahan, LLNL, (510)423-5926
{
}

*********** FRZversion:
versfrz character*19 /"$Revision: 3.21 $"/#  Code version set by CVS

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
mgridrz_nz                integer       # 
mgridrz_xfact(0:mgridrz_nz) _real       # array for deformation factor in X
mgridrz_yfact(0:mgridrz_nz) _real       # array for deformation factor in Y
mgridrz_ngrids            integer  /1/  # number of grids (useful when using mesh refinement)
mgridrz_grid_is(mgridrz_ngrids)   _integer # array id grid associated with grid species 

*********** InjectVars_eq dump:
# variables and functions needed for getting voltage risetime from assumption of 
# constant injection (works with inj_d=2).
inj_phi_eq  _real # Electrostatic potential at the emitting surface at equilibrium.
v_max                real /0./
l_find_rise_time     logical /.false./
afact                real /1./
calc_a               integer  /1/  # determines way of calculating voltage factor for rise time
l_verbose              logical /.true./
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
         #  BASIS-level interface to VPOISRZ, using FRZ database variables
         #  The user program should declare a similar subroutine w/ its vars.
advsc    (schrg:real,eta:real,a:real,kzsq:real,dt:real,dr:real,dz:real,
          nr,nz,phikold:real,taurc:real)       subroutine
         #  Routine that advances the surface charge from time t
         #  to time t+dt in Fourier space
advect   (schrg:real, vbeam:real, dt:real, nz, dz, tmp:real)        subroutine
         #  Advects surface charge with the moving window.
multigridrzf(phi:real,rho:real,nx:integer,nz:integer,dx:real,dz:real,
             mgridrz_accuracy:real) subroutine
         # Multigrid Poisson solver (using "full-multigrid" method, i.e. the 
         # solution is calculatedd at each level using a multigrid procedure 
         # and used as an approximated solution to start the calculation at 
         # the next level)
save_bndstructure_rz(filename:string) subroutine
         # save internal conductor boundary coefficients for each multigrid
         # level
read_bndstructure_rz(filename:string) subroutine
         # read internal conductor boundary coefficients for each multigrid
         # level
get_cond_rz(grid:integer,level:integer) subroutine
         # get internal conductors locations from RZ multigrid solver
setconductorvoltagerz(volt:real,nz:integer,zmmin:real,dz:real,discrete:logical)
         subroutine
         # set voltage on conductors from a z-grid
setconductorvoltagerz_id(id:integer,volt:real) subroutine
         # set voltage on conductor given its ID
calcfact_deform(dz:real,zmin:real,
                xfact:real,yfact:real,nz:integer,ns:integer,is:integer,
                ins:integer,nps:integer,ws:real) subroutine
         # computes factors for elliptical deformation in X and Y planes
init_base(nr:integer,nz:integer,dr:real,dz:real,rmin:real,zmin:real) subroutine
         # initializes the base grid
del_base() subroutine
         # removes the base grid
add_subgrid(id:integer,nr:integer,nz:integer,dr:real,dz:real,
            rmin:real,zmin:real,
            guard_min_r:integer,guard_max_r:integer,
            guard_min_z:integer,guard_max_z:integer) subroutine
         # add a subgrid to the grid id
get_phi_subgrid(id:integer,phi:real,nr:integer,nz:integer) subroutine
         # get the potential of grid id
get_array_subgrid(id:integer,phi:real,nr:integer,nz:integer,which:string) subroutine
         # get the potential of grid id
set_rho_rz(rho:real,nr:integer,nz:integer,id:integer) subroutine
         # set rho of grid id       
get_rho_rz(rho:real,nr:integer,nz:integer,id:integer,rhop:integer) subroutine
         # get rho of grid id
reset_rzmgrid_rho() subroutine
         # sets rho to zero.
find_mgparam_rz() subroutine
         # RZ version of find_mgparam. Does the search for each subgrid.
gchange_rhop_phip_rz() subroutine
         # reallocate rhop and phip arrays
install_conductors_rz() subroutine
         # install conductors data into RZ arrays
set_basegrid_phi() subroutine
	 # set phi on basegrid using w3d.phi 
setbnd_subgrid_to_inj_d() subroutine
         # set indices for force gathering (locpart) to coarser grid 
         # when grid point more than inj_d*grid%dz away from emitting surface
