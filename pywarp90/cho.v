cho
#@(#) File CHO.V, version $Revision: 1.2 $, $Date: 2001/12/07 23:46:44 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package cho of code WARP6
# adaptive mesh refinement package (coupled with CHOMBO package)
# Jean-Luc Vay, LBNL, (510)486-4934
{
}

*********** CHOversion:
verscho character*19 /"$Revision: 1.2 $"/#  Code version version is set by CVS

*********** CHOsubs:
#amrnx integer
#amrny integer
#amrnz integer
#amrnx_ref integer
#amrny_ref integer
#amrnz_ref integer
#phi0(0:amrnx,0:amrny,0:amrnz) _real
#phi1(0:amrnx,0:amrny,0:amrnz) _real
#rho0(0:amrnx,0:amrny,0:amrnz) _real
#rho1(-1:amrnx+1,-1:amrny+1,-1:amrnz+1) _real
#phi0_ref(0:amrnx_ref,0:amrny_ref,0:amrnz_ref) _real
#rho0_ref(0:amrnx_ref,0:amrny_ref,0:amrnz_ref) _real
#phi0xy(0:amrnx,0:amrny) _real
#phi1xy(0:amrnx,0:amrny) _real
#rho0xy(0:amrnx,0:amrny) _real
#rho1xy(-1:amrnx+1,-1:amrny+1) _real
#rho1xyn(0:amrnx,0:amrny) _real
#phi0xy_ref(0:amrnx_ref,0:amrny_ref) _real
#rho0xy_ref(0:amrnx_ref,0:amrny_ref) _real
#amrfieldsolve(phi0:real, phi1:real, rhs0:real, rhs1:real,
#              mx0:integer, mx1:integer, mx2:integer,
#              len0:integer, len1:integer, len2:integer, 
#              off0:integer, off1:integer, off2:integer, 
#              l0:integer, l1:integer, l2:integer,
#              refrat:integer, dx:real)
#     subroutine # External CHOMBO routine for AMR field solve
#amrfieldsolve1grid(phi0:real, rhs0:real,
#              mx0:integer, mx1:integer, mx2:integer,
#              len0:integer, len1:integer, len2:integer, 
#              dx:real)
#     subroutine # External CHOMBO routine for AMR field solve
#getrhs(rhs0:real, rhs1:real,
#       len0:integer, len1:integer, len2:integer, 
#       off0:integer, off1:integer, off2:integer, 
#       l0:integer,   l1:integer,   l2:integer,
#       refrat:integer, dx:real)
#     subroutine # External CHOMBO routine for AMR field solve testing
setamrgrids(nx:integer, ny:integer, nz:integer, dx:real, xmmin:real,
           ymmin:real, zmmin:real,numLevels:integer, refratio:integer,
           i:integer, j:integer, k:integer, level:integer, numtags:integer) subroutine

