frz
#@(#) File FRZ.V, version $Revision: 3.3 $, $Date: 2001/04/30 20:04:40 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package FRZ of code WARP6
# FRZ - fieldsolver package and test driver
# Alex Friedman, LLNL, (510)422-0827
# Debbie Callahan, LLNL, (510)423-5926
{
}

*********** FRZversion:
versfrz character*19 /"$Revision: 3.3 $"/#  Code version set by CVS

*********** FRZvars:
# Variables needed by the test driver of package FRZ
ibc                  character*8 /" "/ #  Boundary conditions-future use
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
