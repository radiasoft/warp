frz
#@(#) File FRZ.V, version $Revision: 3.7 $, $Date: 2001/10/22 17:12:46 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package FRZ of code WARP6
# FRZ - fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
# Debbie Callahan, LLNL, (510)423-5926
{
}

*********** FRZversion:
versfrz character*19 /"$Revision: 3.7 $"/#  Code version set by CVS

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

*********** FRZmgrid:
mgridrz_accuracy          real /1.e-8/  # average accuracy of multigrid solver
mgridrz_ncmax             integer /100/ # maximum number of full multigrid 
                                        # cycles
mgridrz_npre              integer /2/   # number of relaxations steps before 
                                        # coarsening, in multigrid solver
mgridrz_npost             integer /2/   # number of relaxations steps after 
                                        # coarsening, in multigrid solver  
mgridrz_ncycles           integer /2/   # number of multigrid cycles per level
mgridrz_nlevels_max       integer /100/ # maximum number of multigrid levels
mgridrz_nrecurs_min       integer /1/   # minimum level for multigrid recursion
mgridrz_nmeshmin          integer /8/   # minimum number of meshes in ech direction at coarsest level
mgridrz_sub_accuracy      real /1.e-4/  # average accuracy for a sublevel
mgridrz_deform            logical /.false./ # flag for use of elliptic deformation
mgridrz_nz                integer       # 
mgridrz_xfact(0:mgridrz_nz) _real         # array for deformation factor in X
mgridrz_yfact(0:mgridrz_nz) _real         # array for deformation factor in Y

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
             boundxy:integer,bound0:integer,boundnz:integer,
             mgridrz_accuracy:real,mgridrz_ncmax:integer,mgridrz_npre:integer,
             mgridrz_npost:integer,mgridrz_ncycles:integer) subroutine
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
#calcfact_deform(xp:real,yp:real,zp:real,np:integer,dz:real,zmin:real,
calcfact_deform(dz:real,zmin:real,
                xfact:real,yfact:real,nz:integer,ns:integer,is:integer,
                ins:integer,nps:integer,ws:real) subroutine
         # computes factors for elliptical deformation in X and Y planes
