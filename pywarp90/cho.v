cho
#@(#) File CHO.V, version $Revision: 1.14 $, $Date: 2003/05/22 00:13:15 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package cho of code WARP6
# adaptive mesh refinement package (coupled with CHOMBO package)
# Jean-Luc Vay, LBNL, (510)486-4934
{
}

*********** CHOversion:
verscho character*19 /"$Revision: 1.14 $"/#  Code version version is set by CVS

*********** CHOHandle:
cho_handle integer /0/ # Handle to the ChomboPIC package
cho_status integer /0/ # Status from most recent call to ChomboPIC
cho_bunchid integer # Particle bunch id

*********** CHOInput dump:
cho_nlevels  integer /1/ # Number of levels in Chombo hierarchy
cho_refratio integer /2/ # Refinement ratio
cho_maxparticlespercell integer /100000/ # Maximum particles per cell
cho_tagsbuffercells integer /0/ # Depth of buffer around refined cells
cho_fillratio real /0.8/ # Sets blockiness of refined regions, between 0 and 1
cho_tol real /1.e-8/ # Tolerance on field solution
cho_bcflags(3,2) integer /6*0/ # Boundary conditions. 0 Dirichlet, 1 Nuemman
cho_bcvals(3,2) real /6*0./ # Boundary values for grid boundary faces
cho_chargeinterptype integer /1/ # Charge interpolation type, 0 for NGP, 1 for
                                 # linear
cho_debug integer /1/ # Sets amount of debugging output and checking
lcho_allnodesreachable logical /.true./ # When true, all nodes are considered reachable

*********** CHOsubs:
cho_setamrgrids(nx:integer, ny:integer, nz:integer, dx:real, xmmin:real,
           ymmin:real, zmmin:real,numLevels:integer, refratio:integer,
           i:integer, j:integer, k:integer, level:integer, numtags:integer,
           bcxlo:integer, bcxhi:integer, bcylo:integer, bcyhi:integer,
           bczlo:integer, bczhi:integer) subroutine
returnphi(phiout:real) subroutine
returnphic(phiout:real) subroutine
returnrho(rhoout:real) subroutine
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
