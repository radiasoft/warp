em2d
#@(#) File EM2D.V, version $Revision: 1.15 $, $Date: 2007/12/18 19:33:10 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package TOP of code WARP
# TOP - all variables and code which are needed by more than one package.
#       --> These common blocks are available to all packages <--
# Alex Friedman,  LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194
# Jean-Luc Vay,   LBNL, (510)486-4934

*********** EM2D_APML:
pml              integer /1/
pml_sadjusted    integer /2/
apml_exponential integer /3/
apml_hybrid      integer /4/
apml_ssa         integer /5/
apml_lwa         integer /6/

*********** EM2D_bnd dump:
l_pml_cummer  logical    /.false./
s_max_init       real    /4./
s_max_x          real
s_max_y          real
s_delta          real    /5./
sb_coef          real    /0./
nn               real    /2./
bnd_cond         integer /6/
bndexeybz        _type_bnd
bndexeybzc       _type_bnd
bndexeybzf       _type_bnd
bndbxbyez        _type_bnd
bndbxbyezc       _type_bnd
bndbxbyezf       _type_bnd
bndexeybz_c        _type_bnd_cummer
bndexeybzc_c       _type_bnd_cummer
bndexeybzf_c       _type_bnd_cummer
bndbxbyez_c        _type_bnd_cummer
bndbxbyezc_c       _type_bnd_cummer
bndbxbyezf_c       _type_bnd_cummer

*********** EM2D_FIELDobjects dump:
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

init_fields(f:EM2D_FIELDtype, nx:integer, ny:integer, 
           nbndx:integer, nbndy:integer, 
           dtm:real, dx:real, dy:real, clight:real, mu0:real, 
           xmin:real, ymin:real, rap:integer, 
           xlb:integer,ylb:integer,xrb:integer,yrb:integer) subroutine
push_em_e(f:EM2D_FIELDtype,dt:real) subroutine
push_em_b(f:EM2D_FIELDtype,dt:real) subroutine
depose_current_em2d(n:integer,x(n):real,y(n):real,
                    ux(n):real,uy(n):real,uz(n):real,gaminv(n):real,
                    w(n):real,q:real,dt:real,l_particles_weight:logical,
                    field:EM2D_FIELDtype,fpatchfine:EM2D_FIELDtype) subroutine
em2d_geteb2d_linear_serial(n:integer,x(n):real,y(n):real,
                           ex(n):real,ey(n):real,ez(n):real,
                           bx(n):real,by(n):real,bz(n):real,
                           xmin:real,ymin:real,dx:real,dy:real,
                           nx:integer,ny:integer,
                           exg(0:nx+3,0:ny+2):real,eyg(0:nx+3,0:ny+2):real,
                           ezg(0:nx+3,0:ny+2):real,
                           bxg(0:nx+3,0:ny+2):real,byg(0:nx+3,0:ny+2):real,
                           bzg(0:nx+3,0:ny+2):real) subroutine
em2d_getf2d_linear_serial(n:integer,x(n):real,y(n):real,
                          fx(n):real,fy(n):real,fz(n):real,
                          xmin:real,ymin:real,dx:real,dy:real,
                          nx:integer,ny:integer,
                          fxg(0:nx+3,0:ny+2):real,fyg(0:nx+3,0:ny+2):real,
                          fzg(0:nx+3,0:ny+2):real) subroutine
em2d_depose_jxjy_esirkepov_linear_serial(j:real,
                           n:integer,x(n):real,y(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w(n):real,q:real,
                           xmin:real,ymin:real,dt:real,dx:real,dy:real,
                           nx:integer,ny:integer,l_particles_weight:logical)
                           subroutine
getf_em2d(n:integer,x(n):real,y(n):real,
           fx(n):real,fy(n):real,fz(n):real,
           field:EM2D_FIELDtype,fpatchfine:EM2D_FIELDtype,
           which:integer) subroutine   
em2d_step() subroutine
griuni(f:EM2D_FIELDtype) subroutine
grimax(f:EM2D_FIELDtype) subroutine
smooth2d_lindman(q(0:nx+3,0:ny+2),nx,ny) subroutine
move_window_field(f:EM2D_FIELDtype) subroutine
project_j(f:EM2D_FIELDtype,fc:EM2D_FIELDtype,ff:EM2D_FIELDtype) subroutine
set_substitute_fields(field:EM2D_FIELDtype,fpatchcoarse:EM2D_FIELDtype,fpatchfine:EM2D_FIELDtype) subroutine
bndijk(f:EM2D_FIELDtype,j:integer,k:integer) integer function
add_current_slice(f:EM2D_FIELDtype,i:integer) subroutine

%%%%%%%% type_bnd:
n integer
nx integer
ny integer
nbndx integer
nbndy integer
n1x integer
nbot integer
nint integer
ntop integer
nbot1 integer
nbot2 integer
ntop1 integer
ntop2 integer
Ex(1:n) _real
Ey(1:n) _real
Bzx(1:n) _real
Bzy(1:n) _real
aEx(1:n) _real
bEx(1:n) _real
cEx(1:n) _real
aEy(1:n) _real
bEy(1:n) _real
cEy(1:n) _real
aBzx(1:n) _real
bBzx(1:n) _real
cBzx(1:n) _real
aBzy(1:n) _real
bBzy(1:n) _real
cBzy(1:n) _real

%%%%%%%% type_bnd_cummer:
n integer
nx integer
ny integer
nbndx integer
nbndy integer
n1x integer
nbot integer
nint integer
ntop integer
nbot1 integer
nbot2 integer
ntop1 integer
ntop2 integer
Ex(1:n) _real
Ey(1:n) _real
Bz(1:n) _real
Extild(1:n) _real
Eytild(1:n) _real
Bzxtild(1:n) _real
Bzytild(1:n) _real
aEx(1:n) _real
bEx(1:n) _real
cEx(1:n) _real
aEy(1:n) _real
bEy(1:n) _real
cEy(1:n) _real
aBz(1:n) _real
bBzx(1:n) _real
cBzx(1:n) _real
bBzy(1:n) _real
cBzy(1:n) _real
aExtild(1:n) _real
aEytild(1:n) _real
aBzxtild(1:n) _real
aBzytild(1:n) _real
bExtild(1:n) _real
bEytild(1:n) _real
bBzxtild(1:n) _real
bBzytild(1:n) _real

%%%%%%%% EM2D_FIELDtype:
nx integer
ny integer
nxi integer
nyi integer
nxcoeffs integer
nycoeffs integer
nxcopy integer 
nycopy integer 
ntimes integer /1/
xmin real
ymin real
xmax real
ymax real
xminpatch_scatter            real /0./
xmaxpatch_scatter            real /0./
yminpatch_scatter            real /0./
ymaxpatch_scatter            real /0./
xminpatch_gather             real /0./
xmaxpatch_gather             real /0./
yminpatch_gather             real /0./
ymaxpatch_gather             real /0./
rap integer
dx real
dy real
dxi real
dyi real
xlbound                      integer /0/
xrbound                      integer /0/
ylbound                      integer /0/
yrbound                      integer /0/
l_apply_pml        logical /.true./
l_add_source       logical /.true./
l_overcycle_ions   logical /.true./
l_addpatchresidual logical /.false./
l_usecoeffs        logical /.false./
ntemp integer
ipulse integer /1/
npulse integer
js integer
testc real
sinteta real
cst1 real
cst2 real
cj(2) _real
clight real
mu0    real
Ex(0:nx+3,0:ny+2) _real
Ey(0:nx+3,0:ny+2) _real
Ez(0:nx+3,0:ny+2) _real
Bx(0:nx+3,0:ny+2) _real
By(0:nx+3,0:ny+2) _real
Bz(0:nx+3,0:ny+2) _real
aEx(0:nxcoeffs+3,0:nycoeffs+2) _real
bEx(0:nxcoeffs+3,0:nycoeffs+2) _real
cEx(0:nxcoeffs+3,0:nycoeffs+2) _real
dEx(0:nxcoeffs+3,0:nycoeffs+2) _real
aEy(0:nxcoeffs+3,0:nycoeffs+2) _real
bEy(0:nxcoeffs+3,0:nycoeffs+2) _real
cEy(0:nxcoeffs+3,0:nycoeffs+2) _real
dEy(0:nxcoeffs+3,0:nycoeffs+2) _real
aEz(0:nxcoeffs+3,0:nycoeffs+2) _real
bEzx(0:nxcoeffs+3,0:nycoeffs+2) _real
cEzx(0:nxcoeffs+3,0:nycoeffs+2) _real
bEzy(0:nxcoeffs+3,0:nycoeffs+2) _real
cEzy(0:nxcoeffs+3,0:nycoeffs+2) _real
dEz(0:nxcoeffs+3,0:nycoeffs+2) _real
aBx(0:nxcoeffs+3,0:nycoeffs+2) _real
bBx(0:nxcoeffs+3,0:nycoeffs+2) _real
cBx(0:nxcoeffs+3,0:nycoeffs+2) _real
aBy(0:nxcoeffs+3,0:nycoeffs+2) _real
bBy(0:nxcoeffs+3,0:nycoeffs+2) _real
cBy(0:nxcoeffs+3,0:nycoeffs+2) _real
aBz(0:nxcoeffs+3,0:nycoeffs+2) _real
bBzx(0:nxcoeffs+3,0:nycoeffs+2) _real
cBzx(0:nxcoeffs+3,0:nycoeffs+2) _real
bBzy(0:nxcoeffs+3,0:nycoeffs+2) _real
cBzy(0:nxcoeffs+3,0:nycoeffs+2) _real
Excopy(0:nxcopy+3,0:nycopy+2) _real
Eycopy(0:nxcopy+3,0:nycopy+2) _real
Bzcopy(0:nxcopy+3,0:nycopy+2) _real
J(0:nx+3,0:ny+2,3) _real
Jarray(0:nx+3,0:ny+2,3,ntimes) _real
nxfsum integer /0/
nyfsum integer /0/
Exfsum(0:nxfsum+3,0:nyfsum+2) _real
Eyfsum(0:nxfsum+3,0:nyfsum+2) _real
Ezfsum(0:nxfsum+3,0:nyfsum+2) _real
Bxfsum(0:nxfsum+3,0:nyfsum+2) _real
Byfsum(0:nxfsum+3,0:nyfsum+2) _real
Bzfsum(0:nxfsum+3,0:nyfsum+2) _real
rm1(ny) _real
rm2(ny) _real
Bz_in(0:ny+2) _real
Ey_in(0:ny+2) _real
Ex_in(0:ny+2) _real
Ez_in(0:ny+2) _real
By_in(0:ny+2) _real
Bx_in(0:ny+2) _real
nxs integer
nys integer
dirprop integer /0/ # if 1 or -1, will suppress forward/backward emission
Ez_s(0:nxs+3,0:nys+2) _real
Ezx_s(0:nxs+3,0:nys+2) _real
Bz_s(0:nxs+3,0:nys+2) _real
Bzx_s(0:nxs+3,0:nys+2) _real
Ey_sbnd(0:nxs+3,0:nys+2) _real
By_sbnd(0:nxs+3,0:nys+2) _real
laser_profile(0:ny+3) _real
temp(0:ntemp) _real
tpulse(0:npulse+1) _real
pulse(0:npulse+1) _real
bndexeybz _type_bnd
bndbxbyez _type_bnd
bndexeybz_cummer _type_bnd_cummer
bndbxbyez_cummer _type_bnd_cummer
