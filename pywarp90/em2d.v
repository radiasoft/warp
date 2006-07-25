em2d
#@(#) File EM2D.V, version $Revision: 1.5 $, $Date: 2006/07/25 17:31:19 $
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

*********** EM2D_system dump:
mksa                         integer /0/
oneoverk0                    integer /1/
inputunits                   integer /0/
idebe                        integer /1/
idebi                        integer /1/
ifine                        integer /0/
ifini                        integer /0/
nbpart                       integer /0/
nsteps                       integer /1/
itime                        integer /0/
dtfact                         real /0.995/
delta_t                        real /0./
tmin_moving_main_window      real /0./
anormr	                     real /1./
anormj                       real /1./
xlong                          real /1./
ylong                          real /1./
lambda_laser                   real /1./ 
xgauss                         real /1./
ampfac                         real /1./
nmaxw                        integer /1/
ntamp_scatter                integer /2/
ntamp_gather                 integer /4/
xlbound                      integer /0/
xrbound                      integer /0/
ylbound                      integer /0/
yrbound                      integer /0/
pulse_type                   character*8 /"user_def"/
stretch                      real    /1./
misme                        real    /1836./
tiste                        real    /0.1/
vte                          real    /0./
nstepe                       integer /1/
xvide                        real    /3./
npcell_max                   integer /1/
tlance                       real    /1./
teta                         real    /0./
ksmth                        logical /.false./
rho1                         real    /0./ 
rho2                         real    /0./ 
rho3                         real    /0./ 
rho4                         real    /0./ 
rho5                         real    /0./ 
rhomax                       real    /0./ 
xrh1                         real    /0./ 
xrh2                         real    /0./ 
xrh3                         real    /0./ 
xrh4                         real    /0./ 
xrh5                         real    /0./ 
yrh1                         real    /0./ 
yrh2                         real    /0./ 
angle                        real    /0./
cyl                          logical /.false./
xcentre                      real    /0./
ycentre                      real    /0./
iboundary                    integer /1/
proba_cool                   real    /0./
vte_cool                     real    /0./
trapeze                      logical /.false./
esirkepov                    logical /.true./
cyl_cool                     logical /.false./
cyl_cool_centre              logical /.false./
del_collect                  real    /0./
icone                        integer /0/
l_residue                    logical /.false./
type_residue                 integer /1/
l_verbose                    logical /.false./
l_readinputnamelists         logical /.false./
l_getfieldfrombase           logical /.false./
l_absorbx                    logical /.true./
l_absorby                    logical /.true./
l_add_particles              logical /.false./
l_add_particles_random       logical /.false./
l_particles_weight           logical /.false./
l_conserve_global_neutrality logical /.false./ # electrons are reflected at bnd if no ion has left
l_linear_new                 logical /.true./
l_overcycle_ions             logical /.true./
l_push_ions                  logical /.true./

*********** EM2D_FIELDobjects dump:
field        _EM2D_FIELDtype
fpatchcoarse _EM2D_FIELDtype
fpatchfine   _EM2D_FIELDtype
l_onegrid                    logical /.true./
l_elaser_out_plane           logical /.false./
l_moving_window              logical /.false./
l_noinputfield               logical /.false./
l_copyfields                 logical /.false./
ixpatch integer
iypatch integer
ntamp_apr                    integer /4/
rap                          integer /1/
ndelta_t                     integer /1/
xminpatch_scatter            real /0./
xmaxpatch_scatter            real /0./
yminpatch_scatter            real /0./
ymaxpatch_scatter            real /0./
xminpatch_gather             real /0./
xmaxpatch_gather             real /0./
yminpatch_gather             real /0./
ymaxpatch_gather             real /0./
transition_zone              real /0./ # length of zone for linear transition from coarse to fine force (in coarse cell units)
nxfsum integer /0/
nyfsum integer /0/
Exfsum(0:nxfsum+3,0:nyfsum+2) _real
Eyfsum(0:nxfsum+3,0:nyfsum+2) _real
Ezfsum(0:nxfsum+3,0:nyfsum+2) _real
Bxfsum(0:nxfsum+3,0:nyfsum+2) _real
Byfsum(0:nxfsum+3,0:nyfsum+2) _real
Bzfsum(0:nxfsum+3,0:nyfsum+2) _real

*********** EM2D_main dump:
# Variables needed by the main routines of package EMI
nxmain integer /0/
nymain integer /0/
nplanmain integer /0/
rapm integer /0/
nxpload integer /0/
nypload integer /0/
npload(0:nxpload,0:nypload) _real 
xminpload real /0./
xmaxpload real /0./
yminpload real /0./
ymaxpload real /0./
isave_freq integer /100/
l_savefields logical /.false./
l_noacceleration logical /.false./
l_nodeposition logical /.false./
datapath character*(80) /" "/
rhoaux(0:nxmain+3,0:nymain+2,nplanmain) _real

init_fields(f:EM2D_FIELDtype, nx:integer, ny:integer, 
           nbndx:integer, nbndy:integer, 
           dtm:real, dx:real, dy:real, xmin:real, 
           ymin:real, rap:integer, 
           xlb:integer,ylb:integer,xrb:integer,yrb:integer) subroutine
push_em_e(f:EM2D_FIELDtype,dt:real) subroutine
push_em_b(f:EM2D_FIELDtype,dt:real) subroutine
depose_current_em2d(np:integer,xp(np):real,yp(np):real,
                    uxp(np):real,uyp(np):real,uzp(np):real,gaminv(np):real,
                    w(np):real,q:real,dt:real,l_particles_weight:logical) subroutine
geteb_em2d(np:integer,xp(np):real,yp(np):real,
           ex(np):real,ey(np):real,ez(np):real,
           bx(np):real,by(np):real,bz(np):real) subroutine   
em2d_step() subroutine
griuni(f:EM2D_FIELDtype) subroutine
grimax(f:EM2D_FIELDtype) subroutine

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

%%%%%%%% EM2D_FIELDtype:
nx integer
ny integer
nxi integer
nyi integer
nxcopy integer 
nycopy integer 
xmin real
ymin real
rap integer
dx real
dy real
dxi real
dyi real
xlbound                      integer /0/
xrbound                      integer /0/
ylbound                      integer /0/
yrbound                      integer /0/
l_apply_pml logical /.true./
l_add_source logical /.true./
l_overcycle_ions logical /.true./
l_addpatchresidual           logical /.false./
ntemp integer
ipulse integer /1/
npulse integer
js integer
testc real
sinteta real
cst1 real
cst2 real
cj(2) _real
Ex(0:nx+3,0:ny+2) _real
Ey(0:nx+3,0:ny+2) _real
Ez(0:nx+3,0:ny+2) _real
Bx(0:nx+3,0:ny+2) _real
By(0:nx+3,0:ny+2) _real
Bz(0:nx+3,0:ny+2) _real
Excopy(0:nxcopy+3,0:nycopy+2) _real
Eycopy(0:nxcopy+3,0:nycopy+2) _real
Bzcopy(0:nxcopy+3,0:nycopy+2) _real
J(0:nx+3,0:ny+2,3) _real
rm1(ny) _real
rm2(ny) _real
Bz_in(0:ny+2) _real
Ey_in(0:ny+2) _real
Ex_in(0:ny+2) _real
Ez_in(0:nx+3) _real
By_in(0:nx+3) _real
Bx_in(0:nx+3) _real
profx(0:nx+3) _real
profy(0:ny+2) _real
temp(0:ntemp) _real
tpulse(0:npulse+1) _real
pulse(0:npulse+1) _real
bndexeybz _type_bnd
bndbxbyez _type_bnd

