cho
#@(#) File CHO.V, version $Revision: 1.6 $, $Date: 2002/08/28 16:15:48 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package cho of code WARP6
# adaptive mesh refinement package (coupled with CHOMBO package)
# Jean-Luc Vay, LBNL, (510)486-4934
{
}

*********** CHOversion:
verscho character*19 /"$Revision: 1.6 $"/#  Code version version is set by CVS

*********** CHOHandle:
cho_handle integer # Handle to the ChomboPIC package
cho_status integer # Status from most recent call to ChomboPIC
cho_bunchid integer # Particle bunch id

*********** CHOInput dump:
cho_nlevels  integer /1/ # Number of levels in Chombo hierarchy
cho_refratio integer /2/ # Refinement ratio
cho_maxparticlespercell integer /0/ # Maximum particles per cell
cho_tol real /1.e-8/ # Tolerance on field solution
cho_bcflags(3,2) integer /6*0/ # Boundary conditions. 0 Dirichlet, 1 Nuemman

*********** CHOsubs:
#setamrgrids(nx:integer, ny:integer, nz:integer, dx:real, xmmin:real,
#           ymmin:real, zmmin:real,numLevels:integer, refratio:integer,
#           i:integer, j:integer, k:integer, level:integer, numtags:integer,
#           bcxlo:integer, bcxhi:integer, bcylo:integer, bcyhi:integer,
#           bczlo:integer, bczhi:integer) subroutine
#returnphi(phiout:real) subroutine
#returnrho(rhoout:real) subroutine
cho_solve3d(iwhich:integer,nx:integer,ny:integer,nz:integer,nzfull:integer,
            dx:real,dy:real,dz:real,l2symtry,l4symtry,
            xmmin:real,ymmin:real,zmmin:real,zbeam:real,zgrid:real) subroutine
cho_setrho3d(np:integer,xp:real,yp:real,zp:real,zgrid:real,sq:real,sw:real,
             js:integer,ip:integer) subroutine
cho_gete3d(np:integer,xp:real,yp:real,zp:real,zgrid:real,
           ex:real,ey:real,ez:real,js:integer,ip:integer) subroutine
cho_getphi3d(np:integer,xp:real,yp:real,zp:real,zgrid:real,p:real,
             js:integer,ip:integer) subroutine
cho_getrho3d(np:integer,xp:real,yp:real,zp:real,zgrid:real,r:real,
             js:integer,ip:integer) subroutine

reachablenodes(dx:real,mask:integer,xlo:integer,ylo:integer,zlo:integer,
               xhi:integer,yhi:integer,zhi:integer,ncomp:integer) subroutine
coverednodes(dx:real,mask:integer,xlo:integer,ylo:integer,zlo:integer,
             xhi:integer,yhi:integer,zhi:integer) subroutine
nodalcoefficients(dx:real,coeffs:real,xlo:integer,ylo:integer,zlo:integer,
                  xhi:integer,yhi:integer,zhi:integer,ncomp:integer) subroutine
