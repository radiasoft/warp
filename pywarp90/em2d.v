emi
#@(#) File EMI.V, version $Revision: 1.1 $, $Date: 2006/07/18 17:21:07 $
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for package TOP of code WARP
# TOP - all variables and code which are needed by more than one package.
#       --> These common blocks are available to all packages <--
# Alex Friedman,  LLNL, (510)422-0827
# David P. Grote, LLNL, (510)423-7194
# Jean-Luc Vay,   LBNL, (510)486-4934

*********** APML:
pml              integer /1/
pml_sadjusted    integer /2/
apml_exponential integer /3/
apml_hybrid      integer /4/
apml_ssa         integer /5/
apml_lwa         integer /6/

*********** EMIbnd dump:
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

*********** EMIsystem dump:
mksa                         integer /0/
oneoverk0                    integer /1/
inputunits                   integer /0/
absorb                       integer /0/
reflect                      integer /1/
periodic                     integer /2/
idebe                        integer /1/
idebi                        integer /1/
ifine                        integer /0/
ifini                        integer /0/
nbpart                       integer /0/
nsteps                       integer /1/
itime                        integer /0/
time                           real /0./
dt                             real /0.5/
dtfact                         real /0.995/
delta_t                        real /0./
ndelta_t                     integer /1/
tmin_moving_main_window      real /0./
anormr	                     real /1./
anormj                       real /1./
nx                           integer /3/
ny                           integer /2/
xlong                          real /1./
ylong                          real /1./
lambda_laser                   real /1./ 
xgauss                         real /1./
ampfac                         real /1./
nmaxw                        integer /1/
rap                          integer /1/
ntamp_scatter                integer /2/
ntamp_gather                 integer /4/
ntamp_apr                    integer /4/
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
l_onegrid                    logical /.true./
l_residue                    logical /.false./
type_residue                 integer /1/
l_verbose                    logical /.false./
l_readinputnamelists         logical /.false./
l_getfieldfrombase           logical /.false./
l_noinputfield               logical /.false./
l_copyfields                 logical /.false./
l_moving_window              logical /.false./
l_absorbx                    logical /.true./
l_absorby                    logical /.true./
l_add_particles              logical /.false./
l_add_particles_random       logical /.false./
l_particles_weight           logical /.false./
l_conserve_global_neutrality logical /.false./ # electrons are reflected at bnd if no ion has left
l_linear_new                 logical /.true./
l_overcycle_ions             logical /.true./
l_push_ions                  logical /.true./
l_elaser_out_plane           logical /.false./

*********** EMIFIELDobjects dump:
field _EMIFIELDtype

*********** EMIIONobjects dump:
ion _EMIIONtype

*********** EMImain dump:
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
transition_zone real /0/ # length of zone for linear transition from coarse to fine force (in coarse cell units)
rhoaux(0:nxmain+3,0:nymain+2,nplanmain) _real
nxfsum integer /0/
nyfsum integer /0/
nxfsumi integer /0/
nyfsumi integer /0/
Exfsum(0:nxfsum+3,0:nyfsum+2) _real
Eyfsum(0:nxfsum+3,0:nyfsum+2) _real
Ezfsum(0:nxfsum+3,0:nyfsum+2) _real
Bxfsum(0:nxfsum+3,0:nyfsum+2) _real
Byfsum(0:nxfsum+3,0:nyfsum+2) _real
Bzfsum(0:nxfsum+3,0:nyfsum+2) _real
Exfsumi(0:nxfsumi+3,0:nyfsumi+2) _real
Eyfsumi(0:nxfsumi+3,0:nyfsumi+2) _real
Ezfsumi(0:nxfsumi+3,0:nyfsumi+2) _real
Bxfsumi(0:nxfsumi+3,0:nyfsumi+2) _real
Byfsumi(0:nxfsumi+3,0:nyfsumi+2) _real
Bzfsumi(0:nxfsumi+3,0:nyfsumi+2) _real

emi2d_init()  subroutine # 
emi2d_init2()  subroutine # 
propalaser() subroutine #
emi2d_step()  subroutine # 
deprho() subroutine #
dep_electron() subroutine #
dep_ion() subroutine #
step_field() subroutine #
step_field_simple(nt) subroutine #
getex(it:integer) subroutine #
getey(it:integer) subroutine #
getbz(it:integer) subroutine #
getrhojxjy(it:integer) subroutine #
smooth2(q:real,nx:integer,ny:integer) subroutine
check_divemrho() subroutine
dep_rho() subroutine
grimaxall() subroutine
griuniall() subroutine
maxwell(which:integer,t:real) subroutine
gchange_electrons(nbpart:integer) subroutine
move_electrons() subroutine
initfields(which:integer, nx:integer, ny:integer, 
           nbndx:integer, nbndy:integer, 
           dtm:real, dx:real, dy:real, xmin:real, 
           ymin:real, rap:integer, nvect:integer) subroutine
push_em_e(which:integer,dtsdx:real,dtsdy:real) subroutine
push_em_b(which:integer,dtsdx:real,dtsdy:real) subroutine
create_bnd(b:type_bnd, nx:integer, ny:integer, nbndx:integer,
           nbndy:integer, dt:real, dx:real, dy:real) subrountine

%%%%%%%% EMIFIELDtype:
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
l_apply_pml logical /.true./
l_add_source logical /.true./
l_overcycle_ions logical /.true./
ntemp integer
ipulse integer /1/
npulse integer
js integer
testc real
sinteta real
cst1 real
cst2 real
cj(2) _real
potfax real
potfay real
ngmax integer 
ngmaxf2 integer
ngmaxf3 integer
n4x integer
n4xp1 integer
n4xp2 integer
n4xp3 integer
n4xf2 integer
n4xf2p1 integer
n4xf2p2 integer
n4xf2p3 integer
n4xf3 integer
n4xf3p1 integer
n4xf3p2 integer
n4xf3p3 integer
nvect integer
id1(nvect*41) _integer
id2(nvect*17) _integer
id3(nvect*49) _integer
id4(nvect*19) _integer
sd(nvect*41) _real
Ex(0:nx+3,0:ny+2) _real
Ey(0:nx+3,0:ny+2) _real
Ez(0:nx+3,0:ny+2) _real
Bx(0:nx+3,0:ny+2) _real
By(0:nx+3,0:ny+2) _real
Bz(0:nx+3,0:ny+2) _real
Excopy(0:nxcopy+3,0:nycopy+2) _real
Eycopy(0:nxcopy+3,0:nycopy+2) _real
Bzcopy(0:nxcopy+3,0:nycopy+2) _real
Rhojxjy(0:nx+3,0:ny+2,4) _real
Exi(0:nxi+3,0:nyi+2) _real 
Eyi(0:nxi+3,0:nyi+2) _real
Ezi(0:nxi+3,0:nyi+2) _real
Bxi(0:nxi+3,0:nyi+2) _real 
Byi(0:nxi+3,0:nyi+2) _real
Bzi(0:nxi+3,0:nyi+2) _real
Rhoi(0:nxi+3,0:nyi+2) _real
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

%%%%%%%% EMIIONtype:
nbpart integer
nbpartmax integer
nvect integer
rxi(1:nbpartmax) _real
ryi(1:nbpartmax) _real
vxi(1:nbpartmax) _real
vyi(1:nbpartmax) _real
vzi(1:nbpartmax) _real
weight(1:nbpartmax) _real
wexn(1:nvect) _real
weyn(1:nvect) _real
wbzn(1:nvect) _real
jd(nvect+1,16) _integer
w (nvect+1,16) _real
isort_tab1(nvect) _integer 
isort_tab2(nvect) _integer
isort_tab3(nvect) _integer
isort_tab4(nvect) _integer
isort_tab5(nvect) _integer
isort_tab6(nvect) _integer
npi(nvect) _integer
npie(nvect) _integer

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
