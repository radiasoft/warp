em3d
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for the 3-D EM solver of code WARP
# Jean-Luc Vay,   LBNL, (510)486-4934
# David P. Grote, LLNL, (510)423-7194
# Alex Friedman,  LLNL, (510)422-0827

*********** EM3D_APML:
pml              integer /1/
pml_sadjusted    integer /2/
apml_exponential integer /3/
apml_hybrid      integer /4/
apml_ssa         integer /5/
apml_lwa         integer /6/

*********** EM3D_bnd dump:
l_pml_cummer  logical    /.false./
s_max_init       real    /4./
s_max_x          real
s_max_y          real
s_delta          real    /5./
sb_coef          real    /0./
nn               real    /2./
bnd_cond         integer /6/

*********** EM3D_FIELDobjects dump:
l_onegrid                    logical /.true./
l_elaser_out_plane           logical /.false./
l_moving_window              logical /.false./
l_noinputfield               logical /.false./
l_copyfields                 logical /.false./
l_smoothdensity              logical /.false./
ntamp_apr                    integer /4/
rap                          integer /1/
ndelta_t                     integer /1/
nxpatch                      integer /1/
nypatch                      integer /1/
ixpatch                      integer /0/
iypatch                      integer /0/
ntamp_scatter                integer /2/
ntamp_gather                 integer /4/
transition_zone              real /0./ # length of zone for linear transition from coarse to fine force (in coarse cell units)
tmin_moving_main_window      real /0./

init_3dem_block(b:EM3D_BLOCKtype,nx:integer,ny:integer,nz:integer,
                nbndx:integer,nbndy:integer,nbndz:integer,
                nxguard:integer,nyguard:integer,nzguard:integer,
                dt:real,dx:real,dy:real,dz:real,clight:real,mu0:real,
                xmin:real,ymin:real,zmin:real,xlb:real,ylb:real,zlb:real,
                xrb:real,yrb:real,zrb:real) subroutine
push_em3d_e(f:EM3D_FIELDtype,dt:real) subroutine
push_em3d_b(f:EM3D_FIELDtype,dt:real) subroutine
push_em3d_ef(f:EM3D_FIELDtype,dt:real) subroutine
push_em3d_f(f:EM3D_FIELDtype,dt:real) subroutine
push_em3d_block(f:EM3D_BLOCKtype,dt:real,which:integer) subroutine
push_em3d_eef(f:EM3D_BLOCKtype,dt:real,which:integer,l_pushf:logical) subroutine
push_em3d_bf(f:EM3D_BLOCKtype,dt:real,which:integer,l_pushf:logical) subroutine
init_splitfield(sf:EM3D_SPLITYEEFIELDtype, 
                nx:integer,ny:integer,nz:integer, 
                nxguard:integer,nyguard:integer,nzguard:integer, 
                dt:real,dx:real,dy:real,dz:real,clight:real,
                lsx:integer,lsy:integer,lsz:integer, 
                nnx:integer, smaxx:real, sdeltax:real, 
                nny:integer, smaxy:real, sdeltay:real, 
                nnz:integer, smaxz:real, sdeltaz:real) subroutine
depose_jxjyjz_esirkepov_linear_serial(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w(n):real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           l_particles_weight:logical)
                           subroutine
depose_rho_linear_serial(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w(n):real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           l_particles_weight:logical)
                           subroutine
getf3d_linear(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               exg:real,eyg:real,ezg:real)
                           subroutine
yee2node3d(f:EM3D_YEEFIELDtype) subroutine
node2yee3d(f:EM3D_YEEFIELDtype) subroutine
em3d_exchange_j(b:EM3D_BLOCKtype) subroutine
add_current_slice_3d(f:EM3D_YEEFIELDtype,i:integer) subroutine

%%%%%%%% EM3D_SPLITYEEFIELDtype:
fieldtype integer /-2/
nx integer
ny integer
nz integer
nxguard integer /1/
nyguard integer /1/
nzguard integer /1/
dx real
dy real
dz real
dxi real
dyi real
dzi real
dt real
clight real
lsx integer
nnx integer
smaxx real
sdeltax real
lsy integer
nny integer
smaxy real
sdeltay real
lsz integer
nnz integer
smaxz real
sdeltaz real
xlbnd integer /0/
xrbnd integer /0/
ylbnd integer /0/
yrbnd integer /0/
zlbnd integer /0/
zrbnd integer /0/
exx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
exy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
exz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
afx(-nxguard:nx+nxguard) _real
bpfx(-nxguard:nx+nxguard) _real
bmfx(-nxguard:nx+nxguard) _real
agx(-nxguard:nx+nxguard) _real
bpgx(-nxguard:nx+nxguard) _real
bmgx(-nxguard:nx+nxguard) _real
afy(-nyguard:ny+nyguard) _real
bpfy(-nyguard:ny+nyguard) _real
bmfy(-nyguard:ny+nyguard) _real
agy(-nyguard:ny+nyguard) _real
bpgy(-nyguard:ny+nyguard) _real
bmgy(-nyguard:ny+nyguard) _real
afz(-nzguard:nz+nzguard) _real
bpfz(-nzguard:nz+nzguard) _real
bmfz(-nzguard:nz+nzguard) _real
agz(-nzguard:nz+nzguard) _real
bpgz(-nzguard:nz+nzguard) _real
bmgz(-nzguard:nz+nzguard) _real

%%%%%%%% EM3D_YEEFIELDtype:
fieldtype integer /-1/
nx integer
ny integer
nz integer
nxguard integer /1/
nyguard integer /1/
nzguard integer /1/
ntimes integer /1/
xmin real
ymin real
zmin real
xmax real
ymax real
zmax real
dx real
dy real
dz real
dxi real
dyi real
dzi real
clight real
mu0    real
xlbnd integer /0/
xrbnd integer /0/
ylbnd integer /0/
yrbnd integer /0/
zlbnd integer /0/
zrbnd integer /0/
Ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
By(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
F(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Rho(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Rhoarray(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,ntimes) _real
J(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) _real
Jarray(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3,ntimes) _real

%%%%%%%% EM3D_FIELDtype:
fieldtype integer /0/
yf _EM3D_YEEFIELDtype
syf _EM3D_SPLITYEEFIELDtype
proc integer /0/
xl _EM3D_FIELDtype
xr _EM3D_FIELDtype
yl _EM3D_FIELDtype
yr _EM3D_FIELDtype
zl _EM3D_FIELDtype
zr _EM3D_FIELDtype

%%%%%%%% EM3D_BLOCKtype:
nx integer
ny integer
nz integer
nbndx integer
nbndy integer
nbndz integer
xmin real
ymin real
zmin real
xmax real
ymax real
zmax real
dx real
dy real
dz real
dxi real
dyi real
dzi real
xlbnd                      integer /0/
xrbnd                      integer /0/
ylbnd                      integer /0/
yrbnd                      integer /0/
zlbnd                      integer /0/
zrbnd                      integer /0/
core _EM3D_FIELDtype
sidexl _EM3D_FIELDtype 
sidexr _EM3D_FIELDtype
sideyl _EM3D_FIELDtype
sideyr _EM3D_FIELDtype
sidezl _EM3D_FIELDtype
sidezr _EM3D_FIELDtype
edgexlyl _EM3D_FIELDtype
edgexryl _EM3D_FIELDtype
edgexlyr _EM3D_FIELDtype
edgexryr _EM3D_FIELDtype
edgexlzl _EM3D_FIELDtype
edgexlzr _EM3D_FIELDtype
edgexrzl _EM3D_FIELDtype
edgexrzr _EM3D_FIELDtype
edgeylzl _EM3D_FIELDtype
edgeyrzl _EM3D_FIELDtype
edgeylzr _EM3D_FIELDtype
edgeyrzr _EM3D_FIELDtype
cornerxlylzl  _EM3D_FIELDtype
cornerxrylzl  _EM3D_FIELDtype
cornerxlyrzl  _EM3D_FIELDtype
cornerxryrzl  _EM3D_FIELDtype
cornerxlylzr  _EM3D_FIELDtype
cornerxrylzr  _EM3D_FIELDtype
cornerxlyrzr  _EM3D_FIELDtype
cornerxryrzr  _EM3D_FIELDtype
