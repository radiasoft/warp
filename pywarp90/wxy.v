wxy
#@(#) File WXY.V, version $Revision: 3.19 $, $Date: 2002/04/24 21:37:09 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package WXY of code WARP
# XY - PIC package of xy particle code
# Alex Friedman, LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194

*********** WXYversion:
# Quantities associated with version control 
verswxy character*19 /"$Revision: 3.19 $"/ # Current code version, set by CVS

*********** Particlesxy dump parallel:
npmaxxy    integer # Number of particles - same as npmax from TOP
dtp(npmaxxy) _real # Time step size for each particle.
                   # Note that ultimately, this should be included somehow in
                   # the top group Particles

*********** InGenxy dump:
ds       real              # Axial step size, defaults to vbeam*dt
lvzchang logical /.true./  # When true, fancy algorithm is used to estimate dt.
niter_dt integer /4/       # Number of iterations in the calculation of dt
lthick   logical /.false./ # Sets whether to do thick slice model.
                           # When false, does thin slice model.
lvp3d    logical /.false./ # Sets whether or not to use the 3-D field solver
ldiag    logical /.true./  # When false, no diagnostics are done
lexbend  logical /.true./  # When true, use exact transformation in bend, when
                           # false, use same transformation as in 3-D code.

*********** WXYsubs:
# Subroutines in package XY
wxygen() subroutine
wxyexe() subroutine
extebxy(np,xp:real,yp:real,zp:real,uzp:real,gaminv:real,dtl:real,dtr:real,
        bx:real,by:real,bz:real,ex:real,ey:real,ez:real,
        m:real,q:real,bendres:real,bendradi:real,lexbend:logical,gammabar:real,
        zbeam:real,vbeam:real,dt:real,time:real)
             subroutine # Sets external E and B fields
otherexy(np,xp:real,yp:real,dedr:real,dexdx:real,deydy:real,dbdr:real,
         ex:real,ey:real,ez:real,
         bx:real,by:real,bz:real)
             subroutine # Sets external E field
padvncxy(center:string)
             subroutine # Advances particles and rho
fixrhoxy(rho:real,nx,ny,nz,periinz:logical,lthick:logical)
             subroutine # Sums end slices of rho for periodicity
bendezxy(np,xp:real,zp:real,ez:real,bendres:real,bendradi:real,
          bends:logical,bnezflag:logical,linbend:logical)
             subroutine # Corrects axial electric field for warped geometry
epushxy(np,uxp:real,uyp:real,uzp:real,ex:real,ey:real,ez:real,q:real,m:real,
        dtp:real,fdt:real)
             subroutine # Particle velocity advance from E field
bpushxy(np,uxp:real,uyp:real,uzp:real,gaminv:real,bx:real,by:real,bz:real,
        q:real,m:real,dtp:real,fdt:real,bpush:real)
             subroutine # Particle velocity advance from B field
xpushxy(np,xp:real,yp:real,zp:real,uxp:real,uyp:real,uzp:real,gaminv:real,
        dtp:real)
             subroutine # Particle position advance
setcurrxy(curr:real,np,zp:real,uzp:real,gaminv:real,q:real,wght:real,
          zbeam:real,dzz:real,zzmin:real,dz:real,lthick:logical)
             subroutine # Computes current
setrhoxy(rho1d:real,np:integer,xp:real,yp:real,zp:real,zgrid:real,
         uzp:real,gaminv:real,q:real,wght:real)
             subroutine # Computes charge density
fieldsolxy(iwhich:integer)
             subroutine # Complete field solve
vpxy(iwhich:integer)
             subroutine # Call field solver
loadrhoxy(ins,nps,is,lzero:logical)
             subroutine # Simple interface to setrhoxy


