em3d
# Copyright (c) 1990-1998, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
# This is the parameter and variable database for the 3-D EM solver of code WARP
# Jean-Luc Vay,   LBNL, (510)486-4934
# David P. Grote, LLNL, (510)423-7194
# Alex Friedman,  LLNL, (510)422-0827

*********** EM3D_APML:
pml_user         integer /0/
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
bnd_cond         integer /2/

*********** EM3D_kyee dump:
alphax real /0.58333333333333337/  # 7./12.
betaxy real /0.083333333333333329/ # 1./12.
betaxz real /0.083333333333333329/ # 1./12.
gammax real /0.020833333333333332/ # 1./48.
alphay real /0.58333333333333337/  # 7./12.
betayx real /0.083333333333333329/ # 1./12.
betayz real /0.083333333333333329/ # 1./12.
gammay real /0.020833333333333332/ # 1./48.
alphaz real /0.58333333333333337/  # 7./12.
betazx real /0.083333333333333329/ # 1./12.
betazy real /0.083333333333333329/ # 1./12.
gammaz real /0.020833333333333332/ # 1./48.
deltaz real /0.000000000000000000/ # for the lehe solver

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
otherproc                    integer /10/
otherblock                   integer /11/
push_em3d_e(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_b(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_ef(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_f(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_phi(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_a(f:EM3D_YEEFIELDtype,dt:real) subroutine
push_em3d_block(f:EM3D_BLOCKtype,dt:real,which:integer) subroutine
push_em3d_eef(f:EM3D_BLOCKtype,dt:real,which:integer,l_pushf:logical,l_pushpot:logical) subroutine
push_em3d_bf(f:EM3D_BLOCKtype,dt:real,which:integer,l_pushf:logical,l_pushpot:logical) subroutine
push_em3d_blockbnde(b:EM3D_BLOCKtype,dt:real,which:integer) subroutine
push_em3d_blockbndb(b:EM3D_BLOCKtype,dt:real,which:integer) subroutine
push_em3d_blockbndef(b:EM3D_BLOCKtype,dt:real,which:integer) subroutine
push_em3d_blockbndf(b:EM3D_BLOCKtype,dt:real,which:integer) subroutine
scale_em3d_split_fields(sf:EM3D_SPLITYEEFIELDtype,dt:real,l_pushf:logical) subroutine
scale_em3d_bnd_fields(b:EM3D_BLOCKtype,dt:real,l_pushf:logical,l_pushg:logical) subroutine
init_splitfield(sf:EM3D_SPLITYEEFIELDtype, 
                nx:integer,ny:integer,nz:integer, 
                nxguard:integer,nyguard:integer,nzguard:integer, 
                dt:real,dx:real,dy:real,dz:real,
                xmin:real,ymin:real,zmin:real,clight:real,
                lsx:integer,lsy:integer,lsz:integer, 
                nnx:integer, smaxx:real, sdeltax:real, 
                nny:integer, smaxy:real, sdeltay:real, 
                nnz:integer, smaxz:real, sdeltaz:real, 
                l_2dxz:logical, l_1dz:logical, l_2drz:logical,
                norderx:integer,nordery:integer,norderz:integer,
                xcoefs(norderx/2):real,ycoefs(nordery/2):real,zcoefs(norderz/2):real,
                l_nodalgrid:logical,pml_method:integer) subroutine
depose_jxjyjz_esirkepov_linear_serial(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           l_particles_weight:logical)
                           subroutine
depose_jxjyjz_esirkepov_n(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           vx(n):real,vy(n):real,vz(n):real,gaminv(n):real,
                           w:real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           nox:integer,noy:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
depose_jxjyjz_pxpypz_esirkepov_linear_serial(cj:real,mp:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,m:real,
                           xmin:real,ymin:real,zmin:real,
                           dt:real,dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           l_particles_weight:logical,
                           l_relativ:logical)
                           subroutine
depose_rho_linear_serial(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w(n):real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           l_particles_weight:logical)
                           subroutine
depose_rho_n(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w:real,q:real,
                           xmin:real,ymin:real,zmin:real,
                           dx:real,dy:real,dz:real,
                           nx:integer,ny:integer,nz:integer,
                           nxguard:integer,nyguard:integer,nzguard:integer,
                           nox:integer,noy:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
depose_rho_n_2dxz(rho:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w:real,q:real,
                           xmin:real,zmin:real,
                           dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical,l_2drz:logical,
		           type_rz_depose:integer)
                           subroutine
depose_rho_n_2d_circ(rho:real,rho_circ:complex,circ_n:integer,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           w:real,q:real,
                           xmin:real,zmin:real,
                           dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
		           type_rz_depose:integer)
                           subroutine
depose_j_n_1dz(cj:real,
                           n:integer,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           zmin:real,
                           dt:real,dz:real,
                           nz:integer,
                           nzguard:integer,
                           noz:integer,
                           l_particles_weight:logical)
                           subroutine
depose_j_n_2dxz(cj:real,
                           n:integer,x(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,
                           dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
depose_j_n_2dxz_spectral(cj:real,
                           n:integer,x(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,
                           dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
depose_drhoodt_n_2dxz(drhoodt:real,
                           n:integer,x(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,
                           dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
getf3d_linear(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               exg:real,eyg:real,ezg:real)
                           subroutine
getf3d_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         xmin:real,ymin:real,zmin:real,
         dx:real,dy:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noy:integer,noz:integer,
         exg:real,eyg:real,ezg:real,l4symtry:logical)
                           subroutine
averagef3d_rz(nx:integer,ny:integer,nz:integer,
              nxguard:integer,nyguard:integer,nzguard:integer,
              fxg:real,fyg:real,fzg:real,ntheta:integer) subroutine
getf1dz_n(n:integer,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         zmin:real,
         dz:real,
         nz:integer,
         nzguard:integer,
         noz:integer,
         exg:real,eyg:real,ezg:real)
                           subroutine
getf2dxz_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         xmin:real,zmin:real,
         dx:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noz:integer,
         exg:real,eyg:real,ezg:real,l4symtry:logical,l_2drz:logical)
                           subroutine
getfs2dxz_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         fs(n):real,
         xmin:real,zmin:real,
         dx:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noz:integer,
         fsg:real,l4symtry:logical,l_2drz:logical)
                           subroutine
getf2drz_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         xmin:real,zmin:real,
         dx:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noz:integer,
         exg:real,eyg:real,ezg:real,l4symtry:logical,l_2drz:logical)
                           subroutine
getf2drz_circ_n(n:integer,xp(n):real,yp(n):real,zp(n):real,
         ex(n):real,ey(n):real,ez(n):real,
         xmin:real,zmin:real,
         dx:real,dz:real,
         nx:integer,ny:integer,nz:integer,
         nxguard:integer,nyguard:integer,nzguard:integer,
         nox:integer,noz:integer,
	 exg:real,eyg:real,ezg:real,
         exg_circ:complex,eyg_circ:complex,ezg_circ:complex,circ_m:integer)
                           subroutine
gete3d_linear_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               exg:real,eyg:real,ezg:real)
                           subroutine
getb3d_linear_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               bxg:real,byg:real,bzg:real)
                           subroutine
geteb3d_linear_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               exg:real,eyg:real,ezg:real,
               bxg:real,byg:real,bzg:real)
                           subroutine
gete3d_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               nox:integer,noy:integer,noz:integer,
               exg:real,eyg:real,ezg:real,
                           l4symtry:logical,
                           l_lower_order_in_v:logical)
                           subroutine
getb3d_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,ymin:real,zmin:real,
               dx:real,dy:real,dz:real,
               nx:integer,ny:integer,nz:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               nox:integer,noy:integer,noz:integer,
               bxg:real,byg:real,bzg:real,
                           l4symtry:logical,
                           l_lower_order_in_v:logical)
                           subroutine
gete2dxz_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               xmin:real,zmin:real,
               dx:real,dz:real,
               nx:integer,nz:integer,
               nxguard:integer,nzguard:integer,
               nox:integer,noz:integer,
               exg:real,eyg:real,ezg:real,
                           l4symtry:logical,l_2drz:logical,
                           l_lower_order_in_v:logical)
                           subroutine
getb2dxz_n_energy_conserving(n:integer,xp(n):real,yp(n):real,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               xmin:real,zmin:real,
               dx:real,dz:real,
               nx:integer,nz:integer,
               nxguard:integer,nzguard:integer,
               nox:integer,noz:integer,
               bxg:real,byg:real,bzg:real,
                           l4symtry:logical,l_2drz:logical,
                           l_lower_order_in_v:logical)
                           subroutine
gete1dz_n_energy_conserving(n:integer,zp(n):real,
               ex(n):real,ey(n):real,ez(n):real,
               zmin:real,
               dz:real,
               nz:integer,
               nzguard:integer,
               noz:integer,
               exg:real,eyg:real,ezg:real,
                           l_lower_order_in_v:logical)
                           subroutine
getb1dz_n_energy_conserving(n:integer,zp(n):real,
               bx(n):real,by(n):real,bz(n):real,
               zmin:real,
               dz:real,
               nz:integer,
               nzguard:integer,
               noz:integer,
               bxg:real,byg:real,bzg:real,
                           l_lower_order_in_v:logical)
                           subroutine
yee2node3d(f:EM3D_YEEFIELDtype) subroutine
node2yee3d(f:EM3D_YEEFIELDtype) subroutine
Jyee2node3d(f:EM3D_YEEFIELDtype) subroutine
em3d_exchange_e(b:EM3D_BLOCKtype) subroutine
em3d_exchange_b(b:EM3D_BLOCKtype) subroutine
em3d_exchange_f(b:EM3D_BLOCKtype) subroutine
em3d_exchange_j(b:EM3D_BLOCKtype) subroutine
em3d_exchange_rho(b:EM3D_BLOCKtype) subroutine
add_current_slice_3d(f:EM3D_YEEFIELDtype,i:integer) subroutine
add_rho_slice_3d(f:EM3D_YEEFIELDtype,i:integer) subroutine
set_incond(f:EM3D_YEEFIELDtype,n:integer,indx(3,n):integer) subroutine
set_macroscopic_coefs_on_yee(f:EM3D_YEEFIELDtype,n:integer,indx(3,n):integer,
                             sigma:real,epsi:real,mu:real) subroutine
em3d_applybc_rho(f:EM3D_YEEFIELDtype,xlbnd:integer,xrbnd:integer,
                                     ylbnd:integer,yrbnd:integer,
                                     zlbnd:integer,zrbnd:integer,
		                     type_rz_depose:integer) subroutine
em3d_applybc_j(f:EM3D_YEEFIELDtype,xlbnd:integer,xrbnd:integer,
                                   ylbnd:integer,yrbnd:integer,
                                   zlbnd:integer,zrbnd:integer,
		                   type_rz_depose:integer) subroutine
project_jxjyjz(jfine:real,jcoarse:real,jcoarse_mother:real,
               nxf:integer,nyf:integer,nzf:integer,
               nxc:integer,nyc,nzc:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               rapx:integer,rapy:integer,rapz:integer,
               ixc:integer,iyc:integer,izc:integer,l_2dxz:logical,
               icycle:integer,novercycle:integer) subroutine
project_rho(rhofine:real,rhocoarse:real,rhocoarse_mother:real,
               nxf:integer,nyf:integer,nzf:integer,
               nxc:integer,nyc,nzc:integer,
               nxguard:integer,nyguard:integer,nzguard:integer,
               rapx:integer,rapy:integer,rapz:integer,
               ixc:integer,iyc:integer,izc:integer,l_2dxz:logical) subroutine
apply_dmask(rho:real,jc:real,dmaskx:real,dmasky:real,dmaskz:real,
            bounds(6):integer,nguarddepos(3):integer,ntrans(3):integer,
            nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,
            l_pushf:logical,l_2dxz:logical) subroutine
addsubstractfields(child:EM3D_BLOCKtype,child_coarse:EM3D_BLOCKtype,
                   parent:EM3D_BLOCKtype,lc(3):integer,ref(3):integer,l_2dxz:logical) subroutine
addsubstractfields_nodal(child:EM3D_BLOCKtype,child_coarse:EM3D_BLOCKtype,
                   parent:EM3D_BLOCKtype,lc(3):integer,ref(3):integer,l_2dxz:logical) subroutine
shift_em3dblock_ncells_x(b:EM3D_BLOCKtype,n:integer) subroutine
shift_em3dblock_ncells_z(b:EM3D_BLOCKtype,n:integer) subroutine
depose_jxjy_esirkepov_linear_serial_2d(j:real,
                           n:integer,x(n):real,y(n):real,
                           xold(n):real,yold(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,ymin:real,dt:real,dx:real,dy:real,
                           nx:integer,ny:integer,l_particles_weight:logical)
                           subroutine
depose_jxjyjz_esirkepov_n_2d(j:real,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical,l_2drz:logical,
		           type_rz_depose:integer)
                           subroutine
depose_jxjyjz_esirkepov_n_2d_circ(j:real,j_circ:complex,circ_m:integer,
                           n:integer,x(n):real,y(n):real,z(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
		           type_rz_depose:integer)
                           subroutine
depose_jxjyjz_villasenor_n_2d(j:real,
                           n:integer,x(n):real,y(n):real,
                           ux(n):real,uy(n):real,uz(n):real,
                           gaminv(n):real,w:real,q:real,
                           xmin:real,zmin:real,dt:real,dx:real,dz:real,
                           nx:integer,nz:integer,
                           nxguard:integer,nzguard:integer,
                           nox:integer,noz:integer,
                           l_particles_weight:logical,
                           l4symtry:logical)
                           subroutine
setebp(emblock:EM3D_YEEFIELDtype,icycle:integer,novercycle:integer) subroutine
getdive(ex:real,ey:real,ez:real,dive:real,dx:real,dy:real,dz:real,
        nx:integer,ny:integer,nz:integer,nxguard:integer,nyguard:integer,nzguard:integer,
        xmin:real,
        l_2dxz:logical,l_2drz:logical) subroutine

%%%%%%%% EM3D_SPLITYEEFIELDtype:
fieldtype integer /-2/
stencil integer /0/ # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F
pml_method integer /1/
l_nodalgrid logical /.false./    # specifies whether FDTD calculation is on nodal (true) or staggerd (false) grid
nx integer
ny integer
nz integer
nxguard integer /1/
nyguard integer /1/
nzguard integer /1/
nxes integer /0/
nyes integer /0/
nzes integer /0/
nxbs integer /1/
nybs integer /1/
nzbs integer /1/
ixmin integer /0/ # position of first node of grid interior in the x direction (FORTRAN indexing)
iymin integer /0/ # position of first node of grid interior in the y direction (FORTRAN indexing)
izmin integer /0/ # position of first node of grid interior in the z direction (FORTRAN indexing)
ixmax integer /obj__%nx/ # position of last node of grid interior in the x direction (FORTRAN indexing)
iymax integer /obj__%ny/ # position of last node of grid interior in the y direction (FORTRAN indexing)
izmax integer /obj__%nz/ # position of last node of grid interior in the z direction (FORTRAN indexing)
ixming integer /-obj__%nxguard/ # position of first node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iyming integer /-obj__%nyguard/ # position of first node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izming integer /-obj__%nzguard/ # position of first node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
ixmaxg integer /obj__%ixmax+obj__%nxguard/ # position of last node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iymaxg integer /obj__%iymax+obj__%nyguard/ # position of last node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izmaxg integer /obj__%izmax+obj__%nzguard/ # position of last node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
jxmin integer /0/ # position of first node of grid interior in the x direction (Python indexing)
jymin integer /0/ # position of first node of grid interior in the y direction (Python indexing)
jzmin integer /0/ # position of first node of grid interior in the z direction (Python indexing)
jxmax integer /0/ # position of last node of grid interior in the x direction (Python indexing)
jymax integer /0/ # position of last node of grid interior in the y direction (Python indexing)
jzmax integer /0/ # position of last node of grid interior in the z direction (Python indexing)
jxming integer /0/ # position of first node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jyming integer /0/ # position of first node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzming integer /0/ # position of first node of entire grid (interior+guard nodes) in the z direction (Python indexing)
jxmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jymaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the z direction (Python indexing)
nxpo integer /0/
nypo integer /0/
nzpo integer /0/
nconds integer /0/
nxcond integer /0/
nycond integer /0/
nzcond integer /0/
dx real
dy real
dz real
dxi real
dyi real
dzi real
dt real
xmin real
xmax real
ymin real
ymax real
zmin real
zmax real
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
l_1dz logical /.False./
l_2dxz logical /.False./
l_2drz logical /.False./
exx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
exy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
exz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
eyz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ezz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bxz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
byz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
bzz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
fz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
gx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
gy(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
gz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
ax(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
ay(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
az(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
phi(1:3,-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
incond(-nxguard:nxcond+nxguard,-nyguard:nycond+nyguard,-nzguard:nzcond+nzguard) _logical
afx(-nxguard:nx+nxguard) _real
bpfx(-nxguard:nx+nxguard) _real
bmfx(-nxguard:nx+nxguard) _real
sfx(-nxguard:nx+nxguard) _real
agx(-nxguard:nx+nxguard) _real
bpgx(-nxguard:nx+nxguard) _real
bmgx(-nxguard:nx+nxguard) _real
sgx(-nxguard:nx+nxguard) _real
afy(-nyguard:ny+nyguard) _real
bpfy(-nyguard:ny+nyguard) _real
bmfy(-nyguard:ny+nyguard) _real
sfy(-nyguard:ny+nyguard) _real
agy(-nyguard:ny+nyguard) _real
bpgy(-nyguard:ny+nyguard) _real
bmgy(-nyguard:ny+nyguard) _real
sgy(-nyguard:ny+nyguard) _real
afz(-nzguard:nz+nzguard) _real
bpfz(-nzguard:nz+nzguard) _real
bmfz(-nzguard:nz+nzguard) _real
sfz(-nzguard:nz+nzguard) _real
agz(-nzguard:nz+nzguard) _real
bpgz(-nzguard:nz+nzguard) _real
bmgz(-nzguard:nz+nzguard) _real
sgz(-nzguard:nz+nzguard) _real
norderx integer /2/ # order of finite-difference approximation in x
nordery integer /2/ # order of finite-difference approximation in y
norderz integer /2/ # order of finite-difference approximation in z
xcoefs(norderx/2) _real # coefficients of finite-difference stencil in x
ycoefs(nordery/2) _real # coefficients of finite-difference stencil in x
zcoefs(norderz/2) _real # coefficients of finite-difference stencil in x

%%%%%%%% EM3D_YEEFIELDtype:
fieldtype integer /-1/
stencil integer /0/ # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F
spectral logical /.false./
circ_m integer /0/ # max number of azimuthal modes
l_nodalgrid logical /.false./    # specifies whether FDTD calculation is on nodal (true) or staggerd (false) grid
l_nodecentered logical /.false./ # specifies whether field data have been gathered at nodes (when computing with staggered "Yee" grid)
l_macroscopic logical /.false./
nx integer /0/ # nb of mesh cells of grid interior in the x direction
ny integer /0/ # nb of mesh cells of grid interior in the y direction
nz integer /0/ # nb of mesh cells of grid interior in the z direction
nxs integer /0/ # nb of mesh cells of grid interior in the x direction for conductivity array
nys integer /0/ # nb of mesh cells of grid interior in the y direction for conductivity array
nzs integer /0/ # nb of mesh cells of grid interior in the z direction for conductivity array
nxguard integer /1/ # nb of guard cells in the x direction
nyguard integer /1/ # nb of guard cells in the y direction
nzguard integer /1/ # nb of guard cells in the z direction
nxes integer /0/
nyes integer /0/
nzes integer /0/
nxbs integer /1/
nybs integer /1/
nzbs integer /1/
ixmin integer /0/ # position of first node of grid interior in the x direction (FORTRAN indexing)
iymin integer /0/ # position of first node of grid interior in the y direction (FORTRAN indexing)
izmin integer /0/ # position of first node of grid interior in the z direction (FORTRAN indexing)
ixmax integer /obj__%nx/ # position of last node of grid interior in the x direction (FORTRAN indexing)
iymax integer /obj__%ny/ # position of last node of grid interior in the y direction (FORTRAN indexing)
izmax integer /obj__%nz/ # position of last node of grid interior in the z direction (FORTRAN indexing)
ixming integer /-obj__%nxguard/ # position of first node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iyming integer /-obj__%nyguard/ # position of first node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izming integer /-obj__%nzguard/ # position of first node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
ixmaxg integer /obj__%ixmax+obj__%nxguard/ # position of last node of entire grid (interior+guard nodes) in the x direction (FORTRAN indexing)
iymaxg integer /obj__%iymax+obj__%nyguard/ # position of last node of entire grid (interior+guard nodes) in the y direction (FORTRAN indexing)
izmaxg integer /obj__%izmax+obj__%nzguard/ # position of last node of entire grid (interior+guard nodes) in the z direction (FORTRAN indexing)
jxmin integer /0/ # position of first node of grid interior in the x direction (Python indexing)
jymin integer /0/ # position of first node of grid interior in the y direction (Python indexing)
jzmin integer /0/ # position of first node of grid interior in the z direction (Python indexing)
jxmax integer /0/ # position of last node of grid interior in the x direction (Python indexing)
jymax integer /0/ # position of last node of grid interior in the y direction (Python indexing)
jzmax integer /0/ # position of last node of grid interior in the z direction (Python indexing)
jxming integer /0/ # position of first node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jyming integer /0/ # position of first node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzming integer /0/ # position of first node of entire grid (interior+guard nodes) in the z direction (Python indexing)
jxmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the x direction (Python indexing)
jymaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the y direction (Python indexing)
jzmaxg integer /0/ # position of last node of entire grid (interior+guard nodes) in the z direction (Python indexing)
nxp integer /0/
nyp integer /0/
nzp integer /0/
nxdamp integer /0/
nydamp integer /0/
nzdamp integer /0/
nxpnext integer /0/
nypnext integer /0/
nzpnext integer /0/
nxr integer /0/
nyr integer /0/
nzr integer /0/
nxf integer /0/
nyf integer /0/
nzf integer /0/
nxg integer /0/
nyg integer /0/
nzg integer /0/
nxpo integer /0/
nypo integer /0/
nzpo integer /0/
nxmp integer /0/
nymp integer /0/
nzmp integer /0/
nxdrho integer /0/
nydrho integer /0/
nzdrho integer /0/
nxdrhoguard integer /0/
nydrhoguard integer /0/
nzdrhoguard integer /0/
ntimes integer /1/
nconds integer /0/
nxcond integer /0/
nycond integer /0/
nzcond integer /0/
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
theta_damp real /0./
l_1dz logical /.False./
l_2dxz logical /.False./
l_2drz logical /.False./
sigmae real /0./ # coefficient for extended solver
sigmab real /0./ # coefficient for extended solver
Ex(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Ey(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Ez(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Bx(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
By(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Bz(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
Sigmax(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # conductivity at Ex grid location
Sigmay(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # conductivity at Ey grid location
Sigmaz(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # conductivity at Ez grid location
Epsix(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # permittivity at Ex grid location
Epsiy(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # permittivity at Ey grid location
Epsiz(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # permittivity at Ez grid location
Mux(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # permeability at Ex grid location
Muy(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # permeability at Ey grid location
Muz(-nxguard:nxs+nxguard,-nyguard:nys+nyguard,-nzguard:nzs+nzguard) _real # permeability at Ez grid location
Exp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Eyp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Ezp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Bxp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Byp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Bzp(-nxguard:nxp+nxguard,-nyguard:nyp+nyguard,-nzguard:nzp+nzguard) _real
Expnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Eypnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Ezpnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Bxpnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Bypnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
Bzpnext(-nxguard:nxpnext+nxguard,-nyguard:nypnext+nyguard,-nzguard:nzpnext+nzguard) _real
F(-nxguard:nxf+nxguard,-nyguard:nyf+nyguard,-nzguard:nzf+nzguard) _real
G(-nxguard:nxg+nxguard,-nyguard:nyg+nyguard,-nzguard:nzg+nzguard) _real
Rho(-nxguard:nxr+nxguard,-nyguard:nyr+nyguard,-nzguard:nzr+nzguard) _real
Rhoold(-nxguard:nxr+nxguard,-nyguard:nyr+nyguard,-nzguard:nzr+nzguard) _real
Rhoarray(-nxguard:nxr+nxguard,-nyguard:nyr+nyguard,-nzguard:nzr+nzguard,ntimes) _real
DRhoodt(-nxdrhoguard:nxdrho+nxdrhoguard,-nydrhoguard:nydrho+nydrhoguard,-nzdrhoguard:nzdrho+nzdrhoguard) _real
J(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) _real
Jarray(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3,ntimes) _real
incond(-nxguard:nxcond+nxguard,-nyguard:nycond+nyguard,-nzguard:nzcond+nzguard) _logical
Ex_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Ey_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Ez_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Bx_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
By_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Bz_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Exp_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Eyp_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Ezp_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Bxp_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Byp_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Bzp_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
F_circ(-nxguard:nxf+nxguard,-nzguard:nzf+nzguard,1:circ_m) _complex
Rho_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) _complex
Rhoarray_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,ntimes,1:circ_m) _complex
J_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,3,1:circ_m) _complex
Jarray_circ(-nxguard:nx+nxguard,-nzguard:nz+nzguard,3,ntimes,1:circ_m) _complex
Ax(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Ay(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Az(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Phi(-nxguard:nxpo+nxguard,-nyguard:nypo+nyguard,-nzguard:nzpo+nzguard) _real
Mp(-nxguard:nxmp+nxguard,-nyguard:nymp+nyguard,-nzguard:nzmp+nzguard,3) _real
Exold(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Eyold(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Ezold(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Exbar(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Eybar(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Ezbar(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Excp(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Eycp(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
Ezcp(-nxguard:nxdamp+nxguard,-nyguard:nydamp+nyguard,-nzguard:nzdamp+nzguard) _real
E_inx_pos integer /-1/
E_inx_angle real  /0./
E_inx(-nyguard:ny+nyguard,-nzguard:nz+nzguard) _real
E_iny_pos integer /-1/
E_iny_angle real  /0./
E_iny(-nxguard:nx+nxguard,-nzguard:nz+nzguard) _real
E_inz_pos real /0./
E_inz_vel real /0./
E_inz_angle real  /0./
Ex_inz(-nxguard:nx+nxguard,-nyguard:ny+nyguard) _real
Ey_inz(-nxguard:nx+nxguard,-nyguard:ny+nyguard) _real
Ez_inz(-nxguard:nx+nxguard,-nyguard:ny+nyguard) _real
dmaskx(-nxguard:nx+nxguard) _real
dmasky(-nyguard:ny+nyguard) _real
dmaskz(-nzguard:nz+nzguard) _real
ex_stencil(0:4) _real
by_stencil(0:4) _real
norderx integer /2/ # order of finite-difference approximation in x
nordery integer /2/ # order of finite-difference approximation in y
norderz integer /2/ # order of finite-difference approximation in z
xcoefs(norderx/2) _real # coefficients of finite-difference stencil in x
ycoefs(nordery/2) _real # coefficients of finite-difference stencil in x
zcoefs(norderz/2) _real # coefficients of finite-difference stencil in x

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
nxguard integer 
nyguard integer 
nzguard integer 
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
