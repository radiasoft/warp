#include "top.h"
!     Last change:  JLV   3 Jun 2004    0:17 am
!************* MODULE field  **********************************************

module mod_emfield3d
#ifdef MPIPARALLEL
use Parallel
use mpirz
#endif
use EM3D_BLOCKtypemodule
use EM3D_YEEFIELDtypemodule
use EM3D_SPLITYEEFIELDtypemodule
USE EM2D_FIELDobjects
use EM3D_bnd
!use GlobalVars
!use Picglb

implicit none

integer, parameter :: dirichlet=0, neumann=1, periodic=2, openbc=3, yeefield=-1 , splityeefield=-2

INTEGER(ISZ), parameter :: pml = 1, &
                      pml_sadjusted = 2, &
                      apml_exponential = 3, &
                      apml_hybrid = 4, &
                      apml_ssa = 5, &
                      apml_lwa = 6

#ifndef MPIPARALLEL
integer, parameter :: my_index=0
#endif

contains

  subroutine set_bndcoeffsem3d(sf,dt,which)
    TYPE(EM3D_SPLITYEEFIELDtype) :: sf
    real(kind=8) :: sigmax(-sf%nxguard:sf%nx+sf%nxguard),  sigmax_next(-sf%nxguard:sf%nx+sf%nxguard), &
                   lsigmax(-sf%nxguard:sf%nx+sf%nxguard), lsigmax_next(-sf%nxguard:sf%nx+sf%nxguard)
    real(kind=8) :: sigmay(-sf%nyguard:sf%ny+sf%nyguard),  sigmay_next(-sf%nyguard:sf%ny+sf%nyguard), &
                   lsigmay(-sf%nyguard:sf%ny+sf%nyguard), lsigmay_next(-sf%nyguard:sf%ny+sf%nyguard)
    real(kind=8) :: sigmaz(-sf%nzguard:sf%nz+sf%nzguard),  sigmaz_next(-sf%nzguard:sf%nz+sf%nzguard), &
                   lsigmaz(-sf%nzguard:sf%nz+sf%nzguard), lsigmaz_next(-sf%nzguard:sf%nz+sf%nzguard)
    integer(ISZ) :: j,which
    real(kind=8) :: dt
    
    sf%afx = 1.
    sf%agx = 1.
    sf%afy = 1.
    sf%agy = 1.
    sf%afz = 1.
    sf%agz = 1.
    if (which==0) then
      sf%bpfx = sf%clight*dt/sf%dx
      sf%bmfx = -sf%clight*dt/sf%dx
      sf%bpfy = sf%clight*dt/sf%dy
      sf%bmfy = -sf%clight*dt/sf%dy
      sf%bpfz = sf%clight*dt/sf%dz
      sf%bmfz = -sf%clight*dt/sf%dz
      sf%bpgx = sf%clight*dt/sf%dx
      sf%bmgx = -sf%clight*dt/sf%dx
      sf%bpgy = sf%clight*dt/sf%dy
      sf%bmgy = -sf%clight*dt/sf%dy
      sf%bpgz = sf%clight*dt/sf%dz
      sf%bmgz = -sf%clight*dt/sf%dz
      sf%dt = dt
    else
      sf%bpfx = 0.5*sf%clight*dt/sf%dx
      sf%bmfx = -0.5*sf%clight*dt/sf%dx
      sf%bpfy = 0.5*sf%clight*dt/sf%dy
      sf%bmfy = -0.5*sf%clight*dt/sf%dy
      sf%bpfz = 0.5*sf%clight*dt/sf%dz
      sf%bmfz = -0.5*sf%clight*dt/sf%dz
      sf%bpgx = 0.5*sf%clight*dt/sf%dx
      sf%bmgx = -0.5*sf%clight*dt/sf%dx
      sf%bpgy = 0.5*sf%clight*dt/sf%dy
      sf%bmgy = -0.5*sf%clight*dt/sf%dy
      sf%bpgz = 0.5*sf%clight*dt/sf%dz
      sf%bmgz = -0.5*sf%clight*dt/sf%dz
      sf%dt = 0.5*dt
    end if
    if (sf%lsx/=0) then
      sigmax=0.
      sigmax_next=0.
      do j = 0, sf%nx
        sigmax(j)      = (sf%smaxx/sf%dx)*(REAL(j,8)/sf%sdeltax)**sf%nnx
        sigmax_next(j) = (sf%smaxx/sf%dx)*((REAL(j,8)+0.5)/sf%sdeltax)**sf%nnx
        lsigmax(sf%nx-j) = sigmax(j)
        lsigmax_next(sf%nx-j-1) = sigmax_next(j)
      end do
      select case(sf%lsx)
        case(1)
          do j = 0, sf%nx-1
!            write(0,*) 'sigma,sigma_next',sigmax(j),sigmax_next(j),sigmax(j+1)
            call assign_coefs(bnd_cond,sf%afx(j),sf%bpfx(j),sf%bmfx(j),sf%clight*dt,sf%dx,sigmax(j),sigmax_next(j),sb_coef,which)
            call assign_coefs(bnd_cond,sf%agx(j),sf%bpgx(j),sf%bmgx(j),sf%clight*dt,sf%dx,sigmax_next(j),sigmax(j+1),sb_coef,which)
          end do
        case(-1)
          do j = sf%nx, 1, -1
           call assign_coefs(bnd_cond,sf%afx(j),sf%bmfx(j),sf%bpfx(j),sf%clight*dt,sf%dx,lsigmax(j),lsigmax_next(j-1),sb_coef,which)
           call assign_coefs(bnd_cond,sf%agx(j-1),sf%bmgx(j-1),sf%bpgx(j-1),sf%clight*dt,sf%dx,lsigmax_next(j-1),lsigmax(j-1), &
                             sb_coef,which)
          end do
          sf%bmfx(1:sf%nx)=-sf%bmfx(1:sf%nx)
          sf%bpfx(1:sf%nx)=-sf%bpfx(1:sf%nx)
          sf%bmgx(0:sf%nx-1)=-sf%bmgx(0:sf%nx-1)
          sf%bpgx(0:sf%nx-1)=-sf%bpgx(0:sf%nx-1)
       end select
    end if

    if (sf%lsy/=0) then
      sigmay=0.
      sigmay_next=0.
      do j = 0, sf%ny
        sigmay(j)      = (sf%smaxy/sf%dy)*(REAL(j,8)/sf%sdeltay)**sf%nny
        sigmay_next(j) = (sf%smaxy/sf%dy)*((REAL(j,8)+0.5)/sf%sdeltay)**sf%nny
        lsigmay(sf%ny-j) = sigmay(j)
        lsigmay_next(sf%ny-j-1) = sigmay_next(j)
      end do
      select case(sf%lsy)
        case(1)
          do j = 0, sf%ny-1
            call assign_coefs(bnd_cond,sf%afy(j),sf%bpfy(j),sf%bmfy(j),sf%clight*dt,sf%dy,sigmay(j),sigmay_next(j),sb_coef,which)
            call assign_coefs(bnd_cond,sf%agy(j),sf%bpgy(j),sf%bmgy(j),sf%clight*dt,sf%dy,sigmay_next(j),sigmay(j+1),sb_coef,which)
          end do
        case(-1)
          do j = sf%ny, 1, -1
           call assign_coefs(bnd_cond,sf%afy(j),sf%bmfy(j),sf%bpfy(j),sf%clight*dt,sf%dy,lsigmay(j),lsigmay_next(j-1),sb_coef,which)
           call assign_coefs(bnd_cond,sf%agy(j-1),sf%bmgy(j-1),sf%bpgy(j-1),sf%clight*dt,sf%dy,lsigmay_next(j-1),lsigmay(j-1), &
                             sb_coef,which)
          end do
          sf%bmfy(1:sf%ny)=-sf%bmfy(1:sf%ny)
          sf%bpfy(1:sf%ny)=-sf%bpfy(1:sf%ny)
          sf%bmgy(0:sf%ny-1)=-sf%bmgy(0:sf%ny-1)
          sf%bpgy(0:sf%ny-1)=-sf%bpgy(0:sf%ny-1)
       end select
    end if

    if (sf%lsz/=0) then
      sigmaz=0.
      sigmaz_next=0.
      do j = 0, sf%nz
        sigmaz(j)      = (sf%smaxz/sf%dz)*(REAL(j,8)/sf%sdeltaz)**sf%nnz
        sigmaz_next(j) = (sf%smaxz/sf%dz)*((REAL(j,8)+0.5)/sf%sdeltaz)**sf%nnz
        lsigmaz(sf%nz-j) = sigmaz(j)
        lsigmaz_next(sf%nz-j-1) = sigmaz_next(j)
      end do
      select case(sf%lsz)
        case(1)
          do j = 0, sf%nz-1
            call assign_coefs(bnd_cond,sf%afz(j),sf%bpfz(j),sf%bmfz(j),sf%clight*dt,sf%dz,sigmaz(j),sigmaz_next(j),sb_coef,which)
            call assign_coefs(bnd_cond,sf%agz(j),sf%bpgz(j),sf%bmgz(j),sf%clight*dt,sf%dz,sigmaz_next(j),sigmaz(j+1),sb_coef,which)
          end do
        case(-1)
          do j = sf%nz, 1, -1
           call assign_coefs(bnd_cond,sf%afz(j),sf%bmfz(j),sf%bpfz(j),sf%clight*dt,sf%dz,lsigmaz(j),lsigmaz_next(j-1),sb_coef,which)
           call assign_coefs(bnd_cond,sf%agz(j-1),sf%bmgz(j-1),sf%bpgz(j-1),sf%clight*dt,sf%dz,lsigmaz_next(j-1),lsigmaz(j-1), &
                             sb_coef,which)
          end do
          sf%bmfz(1:sf%nz)=-sf%bmfz(1:sf%nz)
          sf%bpfz(1:sf%nz)=-sf%bpfz(1:sf%nz)
          sf%bmgz(0:sf%nz-1)=-sf%bmgz(0:sf%nz-1)
          sf%bpgz(0:sf%nz-1)=-sf%bpgz(0:sf%nz-1)
       end select
    end if

    return 
  end subroutine set_bndcoeffsem3d
    

!************* SUBROUTINE assign_coefs  ********

subroutine assign_coefs(bnd_cond,a,bp,bm,dt,dx,sigma,sigma_next,coef_sigmab,which)
implicit none
REAL(kind=8), INTENT(OUT) :: a,bp,bm
REAL(kind=8), INTENT(IN) :: dt,dx,sigma,sigma_next,coef_sigmab
INTEGER(ISZ),INTENT(IN) :: bnd_cond, which

REAL(kind=8) :: sigma_local, sigmab, sigmab_next, tp, tpp, tm, tmm, g, gp, gm

  tp  = EXP(-sigma*0.5*dx)
  tpp  = EXP(-sigma_next*0.5*dx)
  select case (bnd_cond)
    case (pml,pml_sadjusted)
       IF(bnd_cond==pml) then
        sigma_local = sigma
      else
        sigma_local = MIN(1.e15,abs(tpp-1./tp)/dx)
      END if
      IF(sigma_local == 0.) then
      ! --- end of mesh
        a  =  1.
        bp =  dt / dx
      else
        a  =  EXP(-sigma_local*dt)
        bp =  (1.-a)/(sigma_local*dx)
      END if
      bm =  -bp
    case (apml_exponential)
      sigmab = coef_sigmab*sigma
      IF(sigma == 0.) then
      ! --- end of mesh
        a  =  1.
        bp =  dt / dx
        bm =  -bp
      else
        a  =  EXP(-sigma*dt)
        IF(sigmab==0.) then
          bp =  (1.-a)/(sigma*dx)
          bm =  -bp
        else
          bp =  (sigmab/sigma)*(1.-a)/(1.-EXP(-sigmab*dx))
          bm =  -EXP(-sigmab*dx)*bp
        END if
      END if
    case (apml_hybrid)
      bp = dt/dx
      bm = -dt/dx*(1.+((dx-dt)/(dx+dt))*(1.-tpp))
      a  = 1.+bp*tpp+bm
      bm = tp*bm
    case (apml_ssa,apml_lwa)
      sigmab      = coef_sigmab*sigma
      sigmab_next = coef_sigmab*sigma_next
      tp  = EXP(-(sigma+sigmab)*0.5*dx)
      tmm = EXP(-(sigma-sigmab)*0.5*dx)
      tpp = EXP(-(sigma_next+sigmab_next)*0.5*dx)
      tm  = EXP(-(sigma_next-sigmab_next)*0.5*dx)
      IF(bnd_cond==apml_ssa) then
        g = dx/dt
        a  = -(1-g*(tm+tp)-g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)/(1+g*(tm+tp)+g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)
        bp = 2*tm*(1.+tp*tmm)/(1.+g*(tm+tp)+g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)
        bm = -2*tp*(1.+tm*tpp)/(1.+g*(tm+tp)+g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)
      else
        gp = dx/dt
        gm = gp
        a = -(1-gm-(gp+gm)*tm*tpp-(gp+1)*tm*tmm*tpp*tp)/(1+gm+(gp+gm)*tm*tpp+(gp-1)*tm*tmm*tpp*tp)
        bp = 2*tm*(1.+tp*tmm)/(1+gm+(gp+gm)*tm*tpp+(gp-1)*tm*tmm*tpp*tp)
        bm = -2*tp*(1.+tm*tpp)/(1+gm+(gp+gm)*tm*tpp+(gp-1)*tm*tmm*tpp*tp)
      END if
    case default
      write(0,*) 'Error in assign_coefs: bnd_cond out fo bounds'
  end select

  select case (which)
    case (0)
      ! full time step, do nothing
    case (1)
      ! first half time step
      bp = bp*0.5
      bm = bm*0.5
      a  = (1.+a)*0.5
    case (2)
      ! second half time step
      bp = bp/(1.+a)
      bm = bm/(1.+a)
      a  = 2.*a/(1.+a)
    case default
      write(0,*) 'Error in assign_coefs: which out fo bounds'
      stop
  end select

END subroutine assign_coefs

end module mod_emfield3d

  subroutine init_splitfield(sf, nx, ny, nz, nxguard, nyguard, nzguard, dt, dx, dy, dz, clight, lsx, lsy, lsz, &
                             nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_2dxz)
    use mod_emfield3d
    TYPE(EM3D_SPLITYEEFIELDtype) :: sf
    INTEGER(ISZ), INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard, nnx, nny, nnz, lsx, lsy, lsz
    REAL(kind=8), INTENT(IN) :: dt, dx, dy, dz, clight, smaxx, smaxy, smaxz, sdeltax, sdeltay, sdeltaz
    integer(ISZ) :: j
    logical(ISZ) :: l_2dxz
    
    sf%nx = nx
    sf%ny = ny
    sf%nz = nz
    sf%nxguard = nxguard
    sf%nyguard = nyguard
    sf%nzguard = nzguard
    sf%dx = dx
    sf%dy = dy
    sf%dz = dz
    sf%dxi = 1./dx
    sf%dyi = 1./dy
    sf%dzi = 1./dz
    sf%lsx = lsx
    sf%lsy = lsy
    sf%lsz = lsz
    sf%nnx = nnx
    sf%nny = nny
    sf%nnz = nnz
    sf%smaxx = smaxx
    sf%smaxy = smaxy
    sf%smaxz = smaxz
    sf%sdeltax = sdeltax
    sf%sdeltay = sdeltay
    sf%sdeltaz = sdeltaz
    sf%clight=clight
  ! set min/max of cells positions with FORTRAN indexing
    sf%ixmin = 0
    sf%iymin = 0
    sf%izmin = 0
    sf%ixmax =  sf%nx
    sf%iymax =  sf%ny
    sf%izmax =  sf%nz
    sf%ixming = - sf%nxguard
    sf%iyming = - sf%nyguard
    sf%izming = - sf%nzguard
    sf%ixmaxg =  sf%ixmax+ sf%nxguard
    sf%iymaxg =  sf%iymax+ sf%nyguard
    sf%izmaxg =  sf%izmax+ sf%nzguard
  ! set min/max of cells positions with Python indexing
    sf%jxmin =  sf%ixmin- sf%ixming
    sf%jymin =  sf%iymin- sf%iyming
    sf%jzmin =  sf%izmin- sf%izming
    sf%jxmax =  sf%ixmax- sf%ixming
    sf%jymax =  sf%iymax- sf%iyming
    sf%jzmax =  sf%izmax- sf%izming
    sf%jxming = 0
    sf%jyming = 0
    sf%jzming = 0
    sf%jxmaxg =  sf%ixmaxg- sf%ixming
    sf%jymaxg =  sf%iymaxg- sf%iyming
    sf%jzmaxg =  sf%izmaxg- sf%izming
    sf%l_2dxz = sf%l_2dxz
    call EM3D_SPLITYEEFIELDtypeallot(sf)

    return 
  end subroutine init_splitfield
    
subroutine push_em3d_e(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

if(f%nconds>0) then 
  call push_em3d_conde(f,dt)
  return
end if

dtsdx = f%clight**2*dt/f%dx
dtsdy = f%clight**2*dt/f%dy
dtsdz = f%clight**2*dt/f%dz
mudt  = f%mu0*f%clight**2*dt

if (f%stencil==0 .or. f%stencil==1) then
  call push_em3d_evec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                      mudt,dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%E_inz_pos,f%Ex_inz,f%Ey_inz,f%l_2dxz)
else
  call push_em3d_kyeevec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                         mudt,dtsdx,dtsdy,dtsdz, &
                         f%nx,f%ny,f%nz, &
                         f%nxguard,f%nyguard,f%nzguard,f%E_inz_pos,f%Ex_inz,f%Ey_inz,f%l_2dxz)
end if

return
end subroutine push_em3d_e

subroutine push_em3d_evec(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz, &
                          nxguard,nyguard,nzguard,e_inz_pos,Ex_inz,Ey_inz,l_2dxz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard,E_inz_pos
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard) :: Ex_inz,Ey_inz
real(kind=8), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not. l_2dxz) then
  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l)   - Bz(j,k-1,l  )) &
                            - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
   end do
  end do

  ! advance Ey
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                            - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
   end do
  end do

  ! --- add laser field
  if (E_inz_pos>-1) then
    l=E_inz_pos
    do k = 0, ny-1
      do j = 0, nx
        Ex(j,k,l) = Ex(j,k,l) + Ex_inz(j,k)/2
        Ey(j,k,l) = Ey(j,k,l) + Ey_inz(j,k)/2
!        Ex(j,k,l) = Ey_inz(j,k)
!        Ey(j,k,l) = Ex_inz(j,k)
      end do
    end do
  end if

else

  k = 0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
  end do

  ! advance Ey
  do l = 0, nz
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
  end do

  ! --- add laser field
  if (E_inz_pos>-1) then
    l=E_inz_pos
      do j = 0, nx
        Ex(j,:,l) = Ex(j,:,l) + Ex_inz(j,:)/2
        Ey(j,:,l) = Ey(j,:,l) + Ey_inz(j,:)/2
!        Ex(j,k,l) = Ex_inz(j,k)
!        Ey(j,k,l) = Ey_inz(j,k)
      end do
  end if
end if


return
end subroutine push_em3d_evec

subroutine push_em3d_kyeevec(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz, &
                             nx,ny,nz,nxguard,nyguard,nzguard,e_inz_pos,Ex_inz,Ey_inz,l_2dxz)
use EM3D_kyee
implicit none
integer :: nx,ny,nz,nxguard,nyguard,nzguard,E_inz_pos
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard) :: Ex_inz,Ey_inz
logical(ISZ) :: l_2dxz

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt,E_inz_angle

if (.not.l_2dxz) then
  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + alphay*dtsdy * (Bz(j  ,k,l  ) - Bz(j  ,k-1,l  )) &
                            + betay *dtsdy * (Bz(j+1,k,l  ) - Bz(j+1,k-1,l  ) &
                                           +  Bz(j-1,k,l  ) - Bz(j-1,k-1,l  ) &
                                           +  Bz(j  ,k,l+1) - Bz(j  ,k-1,l+1) &
                                           +  Bz(j  ,k,l-1) - Bz(j  ,k-1,l-1))&
                            + gammay*dtsdy * (Bz(j+1,k,l+1) - Bz(j+1,k-1,l+1) &
                                           +  Bz(j-1,k,l+1) - Bz(j-1,k-1,l+1) &
                                           +  Bz(j+1,k,l-1) - Bz(j+1,k-1,l-1) &
                                           +  Bz(j-1,k,l-1) - Bz(j-1,k-1,l-1)) &
                            - alphaz*dtsdz * (By(j  ,k  ,l) - By(j  ,k  ,l-1)) &
                            - betaz *dtsdz * (By(j+1,k  ,l) - By(j+1,k  ,l-1)  &
                                           +  By(j-1,k  ,l) - By(j-1,k  ,l-1)  &
                                           +  By(j  ,k+1,l) - By(j  ,k+1,l-1)  &
                                           +  By(j  ,k-1,l) - By(j  ,k-1,l-1)) &
                            - gammaz*dtsdz * (By(j+1,k+1,l) - By(j+1,k+1,l-1)  &
                                           +  By(j-1,k+1,l) - By(j-1,k+1,l-1)  &
                                           +  By(j+1,k-1,l) - By(j+1,k-1,l-1)  &
                                           +  By(j-1,k-1,l) - By(j-1,k-1,l-1)) &
                            - 0.5*(alphay+alphaz)*mudt * CJ(j,k,l,1) &
                                     - 0.5*betay *mudt * (CJ(j+1,k,l  ,1) &
                                                     +CJ(j-1,k,l  ,1) &
                                                     +CJ(j  ,k,l+1,1) &
                                                     +CJ(j  ,k,l-1,1)) &
                                     - 0.5*gammay*mudt * (CJ(j+1,k,l+1,1) &
                                                     +CJ(j-1,k,l+1,1) &
                                                     +CJ(j+1,k,l-1,1) &
                                                     +CJ(j-1,k,l-1,1)) &
                                     - 0.5*betaz *mudt * (CJ(j+1,k  ,l,1) &
                                                     +CJ(j-1,k  ,l,1) &
                                                     +CJ(j  ,k+1,l,1) &
                                                     +CJ(j  ,k-1,l,1)) &
                                     - 0.5*gammaz*mudt * (CJ(j+1,k+1,l,1) &
                                                     +CJ(j-1,k+1,l,1) &
                                                     +CJ(j+1,k-1,l,1) &
                                                     +CJ(j-1,k-1,l,1))
    end do
   end do
  end do

  ! advance Ey
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) - alphax*dtsdx * (Bz(j,k  ,l  ) - Bz(j-1,k  ,l  )) &
                                - betax *dtsdx * (Bz(j,k+1,l  ) - Bz(j-1,k+1,l  ) &
                                               +  Bz(j,k-1,l  ) - Bz(j-1,k-1,l  ) &
                                               +  Bz(j,k  ,l+1) - Bz(j-1,k  ,l+1) &
                                               +  Bz(j,k  ,l-1) - Bz(j-1,k  ,l-1)) &
                                - gammax*dtsdx * (Bz(j,k+1,l+1) - Bz(j-1,k+1,l+1) &
                                               +  Bz(j,k-1,l+1) - Bz(j-1,k-1,l+1) &
                                               +  Bz(j,k+1,l-1) - Bz(j-1,k+1,l-1) &
                                               +  Bz(j,k-1,l-1) - Bz(j-1,k-1,l-1)) &                              
                                + alphaz*dtsdz * (Bx(j  ,k  ,l) - Bx(j  ,k  ,l-1)) &
                                + betaz *dtsdz * (Bx(j+1,k  ,l) - Bx(j+1,k  ,l-1) &
                                               +  Bx(j-1,k  ,l) - Bx(j-1,k  ,l-1) &
                                               +  Bx(j  ,k+1,l) - Bx(j  ,k+1,l-1) &
                                               +  Bx(j  ,k-1,l) - Bx(j  ,k-1,l-1)) &
                                + gammaz*dtsdz * (Bx(j+1,k+1,l) - Bx(j+1,k+1,l-1) &
                                               +  Bx(j-1,k+1,l) - Bx(j-1,k+1,l-1) &
                                               +  Bx(j+1,k-1,l) - Bx(j+1,k-1,l-1) &
                                               +  Bx(j-1,k-1,l) - Bx(j-1,k-1,l-1)) &
                                - 0.5*(alphax+alphaz)*mudt * CJ(j,k,l,2) &
                                        - 0.5*betax *mudt * (CJ(j,k+1,l  ,2) &
                                                      +  CJ(j,k-1,l  ,2) &
                                                      +  CJ(j,k  ,l+1,2) &
                                                      +  CJ(j,k  ,l-1,2)) &
                                        - 0.5*gammax*mudt * (CJ(j,k+1,l+1,2) &
                                                      +  CJ(j,k-1,l+1,2) &
                                                      +  CJ(j,k+1,l-1,2) &
                                                      +  CJ(j,k-1,l-1,2)) &
                                        - 0.5*betaz *mudt * (CJ(j  ,k+1,l,2) &
                                                      +  CJ(j  ,k-1,l,2) &
                                                      +  CJ(j+1,k  ,l,2) &
                                                      +  CJ(j-1,k  ,l,2)) &
                                        - 0.5*gammaz*mudt * (CJ(j+1,k+1,l,2) &
                                                      +  CJ(j+1,k-1,l,2) &
                                                      +  CJ(j-1,k+1,l,2) &
                                                      +  CJ(j-1,k-1,l,2))
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + alphax*dtsdx * (By(j,k  ,l  ) - By(j-1,k  ,l  )) &
                                + betax *dtsdx * (By(j,k+1,l  ) - By(j-1,k+1,l  ) &
                                                 +  By(j,k-1,l  ) - By(j-1,k-1,l  ) &
                                                 +  By(j,k  ,l+1) - By(j-1,k  ,l+1) &
                                                 +  By(j,k  ,l-1) - By(j-1,k  ,l-1)) &
                                + gammax*dtsdx * (By(j,k+1,l+1) - By(j-1,k+1,l+1) &
                                                 +  By(j,k-1,l+1) - By(j-1,k-1,l+1) &
                                                 +  By(j,k+1,l-1) - By(j-1,k+1,l-1) &
                                                 +  By(j,k-1,l-1) - By(j-1,k-1,l-1)) &
                                - alphay*dtsdy * (Bx(j  ,k,l  ) - Bx(j  ,k-1,l  )) &
                                - betay *dtsdy * (Bx(j+1,k,l  ) - Bx(j+1,k-1,l  ) &
                                                 +  Bx(j-1,k,l  ) - Bx(j-1,k-1,l  ) &
                                                 +  Bx(j  ,k,l+1) - Bx(j  ,k-1,l+1) &
                                                 +  Bx(j  ,k,l-1) - Bx(j  ,k-1,l-1)) &
                                - gammay*dtsdy * (Bx(j+1,k,l+1) - Bx(j+1,k-1,l+1) &
                                                 +  Bx(j-1,k,l+1) - Bx(j-1,k-1,l+1) &
                                                 +  Bx(j+1,k,l-1) - Bx(j+1,k-1,l-1) &
                                                 +  Bx(j-1,k,l-1) - Bx(j-1,k-1,l-1)) &
                                - 0.5*(alphax+alphay)*mudt * CJ(j,k,l,3) &
                                        - 0.5*betax *mudt * (CJ(j,k+1,l  ,3) &
                                                      +  CJ(j,k-1,l  ,3) &
                                                      +  CJ(j,k  ,l+1,3) &
                                                      +  CJ(j,k  ,l-1,3)) &
                                        - 0.5*gammax*mudt * (CJ(j,k+1,l+1,3) &
                                                      +  CJ(j,k-1,l+1,3) &
                                                      +  CJ(j,k+1,l-1,3) &
                                                      +  CJ(j,k-1,l-1,3)) &
                                        - 0.5*betay *mudt * (CJ(j+1,k,l  ,3) &
                                                      +  CJ(j-1,k,l  ,3) &
                                                      +  CJ(j  ,k,l+1,3) &
                                                      +  CJ(j  ,k,l-1,3)) &
                                        - 0.5*gammay*mudt * (CJ(j+1,k,l+1,3) &
                                                      +  CJ(j-1,k,l+1,3) &
                                                      +  CJ(j+1,k,l-1,3) &
                                                      +  CJ(j-1,k,l-1,3))
    end do
   end do
  end do

  ! --- add laser field
  if (E_inz_pos>-1) then
    l=E_inz_pos
    do k = 0, ny-1
      do j = 0, nx
        Ex(j,k,l) = Ex(j,k,l) + Ex_inz(j,k)/2
        Ey(j,k,l) = Ey(j,k,l) + Ey_inz(j,k)/2
!        Ex(j,k,l) = Ey_inz(j,k)
!        Ey(j,k,l) = Ex_inz(j,k)
      end do
    end do
  end if

else
  k = 0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) -     alphaz*dtsdz * (By(j  ,k  ,l) - By(j  ,k  ,l-1)) &
                            - 2.5*betaz *dtsdz * (By(j+1,k  ,l) - By(j+1,k  ,l-1)  &
                                               +  By(j-1,k  ,l) - By(j-1,k  ,l-1))  &
                            - alphaz*mudt       * CJ(j,k,l,1) &
                            - 2.5*betaz*mudt    * (CJ(j+1,k  ,l,1)+CJ(j-1,k  ,l,1) )
    end do
  end do

  ! advance Ey
  do l = 0, nz
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) - alphax*dtsdx * (Bz(j,k  ,l  ) - Bz(j-1,k  ,l  )) &
                                - betax *dtsdx * (Bz(j,k  ,l+1) - Bz(j-1,k  ,l+1) &
                                               +  Bz(j,k  ,l-1) - Bz(j-1,k  ,l-1)) &
                                + alphaz*dtsdz * (Bx(j  ,k  ,l) - Bx(j  ,k  ,l-1)) &
                                + betaz *dtsdz * (Bx(j+1,k  ,l) - Bx(j+1,k  ,l-1) &
                                               +  Bx(j-1,k  ,l) - Bx(j-1,k  ,l-1)) &
                                - 0.5*(alphax+alphaz)*mudt * CJ(j,k,l,2) &
                                        - 0.5*2.5*betax *mudt * (CJ(j,k  ,l+1,2) &
                                                      +  CJ(j,k  ,l-1,2)) &
                                        - 0.5*2.5*betaz *mudt * (CJ(j+1,k  ,l,2) &
                                                      +  CJ(j-1,k  ,l,2)) 
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + alphax*dtsdx * (By(j,k  ,l  ) - By(j-1,k  ,l  )) &
                                + betax *dtsdx * ( By(j,k  ,l+1) - By(j-1,k  ,l+1) &
                                                 +  By(j,k  ,l-1) - By(j-1,k  ,l-1)) &
                                - alphax*mudt * CJ(j,k,l,3) &
                                        - 2.5*betax *mudt * (CJ(j,k  ,l+1,3) &
                                                      +  CJ(j,k  ,l-1,3))  
    end do
  end do

  ! --- add laser field
  if (E_inz_pos>-1) then
    l=E_inz_pos
      do j = 0, nx
        Ex(j,k,l) = Ex(j,k,l) + Ex_inz(j,k)/2
        Ey(j,k,l) = Ey(j,k,l) + Ey_inz(j,k)/2
!        Ex(j,k,l) = Ey_inz(j,k)
!        Ey(j,k,l) = Ex_inz(j,k)
      end do
  end if

end if

return
end subroutine push_em3d_kyeevec

subroutine push_em3d_b(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz

dtsdx = dt/f%dx
dtsdy = dt/f%dy
dtsdz = dt/f%dz

if (f%stencil==0 .or. f%stencil==2) then
  call push_em3d_bvec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
else
  call push_em3d_kyeebvec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
end if

return
end subroutine push_em3d_b

subroutine push_em3d_bvec(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  ! advance Bx
  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j,k+1,l  ) - Ez(j,k,l)) &
                            + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
    end do
   end do
  end do

  ! advance By
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &  
                            - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
    end do
   end do
  end do

  ! advance Bz 
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j,k,l)) &
                            + dtsdy * (Ex(j,k+1,l) - Ex(j,k,l))
    end do
   end do
  end do

else

  k=0
  ! advance Bx
  do l = 0, nz-1
    do j = 0, nx
      Bx(j,k,l) = Bx(j,k,l) + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
    end do
  end do

  ! advance By
  do l = 0, nz-1
    do j = 0, nx-1
      By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &  
                            - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
    end do
  end do

  ! advance Bz 
  do l = 0, nz
    do j = 0, nx-1
      Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j,k,l)) 
    end do
  end do

end if

return
end subroutine push_em3d_bvec

subroutine push_em3d_kyeebvec(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
use EM3D_kyee
implicit none
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  ! advance Bx
  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      Bx(j,k,l) = Bx(j,k,l) - alphay*dtsdy * (Ez(j  ,k+1,l  ) - Ez(j  ,k  ,l  )) &
                            -  betay*dtsdy * (Ez(j+1,k+1,l  ) - Ez(j+1,k  ,l  ) &
                                           +  Ez(j-1,k+1,l  ) - Ez(j-1,k  ,l  ) &
                                           +  Ez(j  ,k+1,l+1) - Ez(j  ,k  ,l+1) &
                                           +  Ez(j  ,k+1,l-1) - Ez(j  ,k  ,l-1)) &
                            - gammay*dtsdy * (Ez(j+1,k+1,l+1) - Ez(j+1,k  ,l+1) &
                                           +  Ez(j-1,k+1,l+1) - Ez(j-1,k  ,l+1) &
                                           +  Ez(j+1,k+1,l-1) - Ez(j+1,k  ,l-1) &
                                           +  Ez(j-1,k+1,l-1) - Ez(j-1,k  ,l-1)) &
                            + alphaz*dtsdz * (Ey(j  ,k  ,l+1) - Ey(j  ,k  ,l  )) &
                            +  betaz*dtsdz * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  ) &
                                           +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  ) &
                                           +  Ey(j  ,k+1,l+1) - Ey(j  ,k+1,l  ) &
                                           +  Ey(j  ,k-1,l+1) - Ey(j  ,k-1,l  )) &
                            + gammaz*dtsdz * (Ey(j+1,k+1,l+1) - Ey(j+1,k+1,l  ) &
                                           +  Ey(j-1,k+1,l+1) - Ey(j-1,k+1,l  ) &
                                           +  Ey(j+1,k-1,l+1) - Ey(j+1,k-1,l  ) &
                                           +  Ey(j-1,k-1,l+1) - Ey(j-1,k-1,l  )) 
    end do
   end do
  end do

  ! advance By
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      By(j,k,l) = By(j,k,l) + alphax*dtsdx * (Ez(j+1,k  ,l  ) - Ez(j  ,k  ,l  )) &  
                            +  betax*dtsdx * (Ez(j+1,k+1,l  ) - Ez(j  ,k+1,l  ) &
                                           +  Ez(j+1,k-1,l  ) - Ez(j  ,k-1,l  ) &
                                           +  Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1) &
                                           +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                            + gammax*dtsdx * (Ez(j+1,k+1,l+1) - Ez(j  ,k+1,l+1) &
                                           +  Ez(j+1,k-1,l+1) - Ez(j  ,k-1,l+1) &
                                           +  Ez(j+1,k+1,l-1) - Ez(j  ,k+1,l-1) &
                                           +  Ez(j+1,k-1,l-1) - Ez(j  ,k-1,l-1)) &
                            - alphaz*dtsdz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                            -  betaz*dtsdz * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  ) &
                                           +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  ) &
                                           +  Ex(j  ,k+1,l+1) - Ex(j  ,k+1,l  ) &
                                           +  Ex(j  ,k-1,l+1) - Ex(j  ,k-1,l  )) &
                            - gammaz*dtsdz * (Ex(j+1,k+1,l+1) - Ex(j+1,k+1,l  ) &
                                           +  Ex(j-1,k+1,l+1) - Ex(j-1,k+1,l  ) &
                                           +  Ex(j+1,k-1,l+1) - Ex(j+1,k-1,l  ) &
                                           +  Ex(j-1,k-1,l+1) - Ex(j-1,k-1,l  )) 
    end do
   end do
  end do

  ! advance Bz 
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      Bz(j,k,l) = Bz(j,k,l) - alphax*dtsdx * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  )) &
                            -  betax*dtsdx * (Ey(j+1,k+1,l  ) - Ey(j  ,k+1,l  ) &
                                           +  Ey(j+1,k-1,l  ) - Ey(j  ,k-1,l  ) &
                                           +  Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1) &
                                           +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) &
                            - gammax*dtsdx * (Ey(j+1,k+1,l+1) - Ey(j  ,k+1,l+1) &
                                           +  Ey(j+1,k-1,l+1) - Ey(j  ,k-1,l+1) &
                                           +  Ey(j+1,k+1,l-1) - Ey(j  ,k+1,l-1) &
                                           +  Ey(j+1,k-1,l-1) - Ey(j  ,k-1,l-1)) &
                            + alphay*dtsdy * (Ex(j  ,k+1,l  ) - Ex(j  ,k  ,l  )) &
                            +  betay*dtsdy * (Ex(j+1,k+1,l  ) - Ex(j+1,k  ,l  ) &
                                           +  Ex(j-1,k+1,l  ) - Ex(j-1,k  ,l  ) &
                                           +  Ex(j  ,k+1,l+1) - Ex(j  ,k  ,l+1) &
                                           +  Ex(j  ,k+1,l-1) - Ex(j  ,k  ,l-1)) &
                            + gammay*dtsdy * (Ex(j+1,k+1,l+1) - Ex(j+1,k  ,l+1) &
                                           +  Ex(j-1,k+1,l+1) - Ex(j-1,k  ,l+1) &
                                           +  Ex(j+1,k+1,l-1) - Ex(j+1,k  ,l-1) &
                                           +  Ex(j-1,k+1,l-1) - Ex(j-1,k  ,l-1)) 
    end do
   end do
  end do

else

  k=0
  ! advance Bx
  do l = 0, nz-1
    do j = 0, nx
      Bx(j,k,l) = Bx(j,k,l) + alphaz*dtsdz * (Ey(j  ,k  ,l+1) - Ey(j  ,k  ,l  )) &
                            +  2.5*betaz*dtsdz * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  ) &
                                           +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  )) 
    end do
  end do

  ! advance By
  do l = 0, nz-1
    do j = 0, nx-1
      By(j,k,l) = By(j,k,l) + alphax*dtsdx * (Ez(j+1,k  ,l  ) - Ez(j  ,k  ,l  )) &  
                            +  2.5*betax*dtsdx * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1) &
                                           +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                            - alphaz*dtsdz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                            -  2.5*betaz*dtsdz * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  ) &
                                           +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  )) 
    end do
  end do

  ! advance Bz 
  do l = 0, nz
    do j = 0, nx-1
      Bz(j,k,l) = Bz(j,k,l) - alphax*dtsdx * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  )) &
                            -  2.5*betax*dtsdx * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1) &
                                           +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) 
    end do
  end do

end if

return
end subroutine push_em3d_kyeebvec

subroutine push_em3d_f(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,dtsepsi

if(f%nconds>0) then 
  call push_em3d_condf(f,dt)
  return
end if

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz
dtsepsi = f%mu0*f%clight**3*dt

if (f%stencil==0 .or. f%stencil==1) then
  call push_em3d_fvec(f%ex,f%ey,f%ez,f%f, f%rho, &
                      dtsepsi,dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
else
  call push_em3d_kyeefvec(f%ex,f%ey,f%ez,f%f, f%rho, &
                      dtsepsi,dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
endif

end subroutine push_em3d_f

subroutine push_em3d_fvec(ex,ey,ez,f,rho,dtsepsi,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f,rho
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz,dtsepsi
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
      F(j,k,l) = F(j,k,l) + dtsdx * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                          + dtsdy * (Ey(j,k,l) - Ey(j  ,k-1,l  )) &
                          + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                          - dtsepsi * Rho(j,k,l)
    end do
   end do
  end do

else

  k=0
  do l = 0, nz
    do j = 0, nx
      F(j,k,l) = F(j,k,l) + dtsdx * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                          + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                          - dtsepsi * Rho(j,k,l)
    end do
  end do

end if

return
end subroutine push_em3d_fvec

subroutine push_em3d_kyeefvec(ex,ey,ez,f,rho,dtsepsi,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
use EM3D_kyee
implicit none
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f,rho
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz,dtsepsi
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
      F(j,k,l) = F(j,k,l) + alphax*dtsdx * (Ex(j  ,k  ,l  ) - Ex(j-1,k  ,l  )) &
                          +  betax*dtsdx * (Ex(j  ,k+1,l  ) - Ex(j-1,k+1,l  ) &
                                         +  Ex(j  ,k-1,l  ) - Ex(j-1,k-1,l  ) &
                                         +  Ex(j  ,k  ,l+1) - Ex(j-1,k  ,l+1) &
                                         +  Ex(j  ,k  ,l-1) - Ex(j-1,k  ,l-1)) &
                          + gammax*dtsdx * (Ex(j  ,k+1,l+1) - Ex(j-1,k+1,l+1) &
                                         +  Ex(j  ,k-1,l+1) - Ex(j-1,k-1,l+1) &
                                         +  Ex(j  ,k+1,l-1) - Ex(j-1,k+1,l-1) &
                                         +  Ex(j  ,k-1,l-1) - Ex(j-1,k-1,l-1)) &
                          + alphay*dtsdy * (Ey(j  ,k  ,l  ) - Ey(j  ,k-1,l  )) &
                          +  betay*dtsdy * (Ey(j+1,k  ,l  ) - Ey(j+1,k-1,l  ) &
                                         +  Ey(j-1,k  ,l  ) - Ey(j-1,k-1,l  ) &
                                         +  Ey(j  ,k  ,l+1) - Ey(j  ,k-1,l+1) &
                                         +  Ey(j  ,k  ,l-1) - Ey(j  ,k-1,l-1)) &
                          + gammay*dtsdy * (Ey(j+1,k  ,l+1) - Ey(j+1,k-1,l+1) &
                                         +  Ey(j-1,k  ,l+1) - Ey(j-1,k-1,l+1) &
                                         +  Ey(j+1,k  ,l-1) - Ey(j+1,k-1,l-1) &
                                         +  Ey(j-1,k  ,l-1) - Ey(j-1,k-1,l-1)) &
                          + alphaz*dtsdz * (Ez(j  ,k  ,l  ) - Ez(j  ,k  ,l-1)) &
                          +  betaz*dtsdz * (Ez(j+1,k  ,l  ) - Ez(j+1,k  ,l-1) &
                                         +  Ez(j-1,k  ,l  ) - Ez(j-1,k  ,l-1) &
                                         +  Ez(j  ,k+1,l  ) - Ez(j  ,k+1,l-1) &
                                         +  Ez(j  ,k-1,l  ) - Ez(j  ,k-1,l-1)) &
                          + gammaz*dtsdz * (Ez(j+1,k+1,l  ) - Ez(j+1,k+1,l-1) &
                                         +  Ez(j-1,k+1,l  ) - Ez(j-1,k+1,l-1) &
                                         +  Ez(j+1,k-1,l  ) - Ez(j+1,k-1,l-1) &
                                         +  Ez(j-1,k-1,l  ) - Ez(j-1,k-1,l-1)) &

                          - dtsepsi/3. * ( (alphax+alphay+alphaz)* Rho(j,k,l) &
                                        +   betax*(Rho(j  ,k+1,l  )+Rho(j  ,k-1,l  )+Rho(j  ,k  ,l+1)+Rho(j  ,k  ,l-1)) &
                                        +  gammax*(Rho(j  ,k+1,l+1)+Rho(j  ,k-1,l+1)+Rho(j  ,k+1,l-1)+Rho(j  ,k-1,l-1)) &
                                        +   betay*(Rho(j+1,k  ,l  )+Rho(j-1,k  ,l  )+Rho(j  ,k  ,l+1)+Rho(j  ,k  ,l-1)) &
                                        +  gammay*(Rho(j+1,k  ,l+1)+Rho(j-1,k  ,l+1)+Rho(j+1,k  ,l-1)+Rho(j-1,k  ,l-1)) &
                                        +   betaz*(Rho(j+1,k  ,l  )+Rho(j-1,k  ,l  )+Rho(j  ,k+1,l  )+Rho(j  ,k-1,l  )) &
                                        +  gammaz*(Rho(j+1,k+1,l  )+Rho(j-1,k+1,l  )+Rho(j+1,k-1,l  )+Rho(j-1,k-1,l  )) )
    end do
   end do
  end do

else

  k=0
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
      F(j,k,l) = F(j,k,l) + alphax*dtsdx * (Ex(j  ,k  ,l  ) - Ex(j-1,k  ,l  )) &
                          +  2.5*betax*dtsdx * (Ex(j  ,k  ,l+1) - Ex(j-1,k  ,l+1) &
                                         +  Ex(j  ,k  ,l-1) - Ex(j-1,k  ,l-1)) &
                          + alphaz*dtsdz * (Ez(j  ,k  ,l  ) - Ez(j  ,k  ,l-1)) &
                          +  2.5*betaz*dtsdz * (Ez(j+1,k  ,l  ) - Ez(j+1,k  ,l-1) &
                                         +  Ez(j-1,k  ,l  ) - Ez(j-1,k  ,l-1)) &
                          - dtsepsi/2. * ( (alphax+alphaz)* Rho(j,k,l) &
                                        +   2.5*betax*(Rho(j  ,k  ,l+1)+Rho(j  ,k  ,l-1)) &
                                        +   2.5*betaz*(Rho(j+1,k  ,l  )+Rho(j-1,k  ,l  )) )
    end do
   end do
  end do

end if

return
end subroutine push_em3d_kyeefvec

subroutine push_em3d_ef(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

if(f%nconds>0) then 
  call push_em3d_condef(f,dt)
  return
end if

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz

if (f%stencil==0 .or. f%stencil==2) then
  call push_em3d_efvec(f%ex,f%ey,f%ez,f%f, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
else
  call push_em3d_kyeeefvec(f%ex,f%ey,f%ez,f%f, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
endif

return
end subroutine push_em3d_ef

subroutine push_em3d_efvec(ex,ey,ez,f,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + dtsdx * (F(j+1,k,l) - F(j,k,l)) 
    end do
   end do
  end do

  ! advance Ey
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) + dtsdy * (F(j,k+1,l) - F(j,k,l))
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + dtsdz * (F(j,k,l+1) - F(j,k,l)) 
    end do
   end do
  end do

else

  k=0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + dtsdx * (F(j+1,k,l) - F(j,k,l)) 
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + dtsdz * (F(j,k,l+1) - F(j,k,l)) 
    end do
  end do

end if

return
end subroutine push_em3d_efvec

subroutine push_em3d_kyeeefvec(ex,ey,ez,f,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
use EM3D_kyee
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + alphax*dtsdx * (F(j+1,k  ,l  ) - F(j  ,k  ,l  )) &
                            +  betax*dtsdx * (F(j+1,k+1,l  ) - F(j  ,k+1,l  ) &
                                           +  F(j+1,k-1,l  ) - F(j  ,k-1,l  ) &
                                           +  F(j+1,k  ,l+1) - F(j  ,k  ,l+1) &
                                           +  F(j+1,k  ,l-1) - F(j  ,k  ,l-1)) &
                            + gammax*dtsdx * (F(j+1,k+1,l+1) - F(j  ,k+1,l+1) &
                                           +  F(j+1,k-1,l+1) - F(j  ,k-1,l+1) &
                                           +  F(j+1,k+1,l-1) - F(j  ,k+1,l-1) &
                                           +  F(j+1,k-1,l-1) - F(j  ,k-1,l-1)) 
    end do
   end do
  end do

  ! advance Ey
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) + alphay*dtsdy * (F(j  ,k+1,l  ) - F(j  ,k  ,l  )) &
                            +  betay*dtsdy * (F(j+1,k+1,l  ) - F(j+1,k  ,l  ) &
                                           +  F(j-1,k+1,l  ) - F(j-1,k  ,l  ) &
                                           +  F(j  ,k+1,l+1) - F(j  ,k  ,l+1) &
                                           +  F(j  ,k+1,l-1) - F(j  ,k  ,l-1)) &
                            + gammay*dtsdy * (F(j+1,k+1,l+1) - F(j+1,k  ,l+1) &
                                           +  F(j-1,k+1,l+1) - F(j-1,k  ,l+1) &
                                           +  F(j+1,k+1,l-1) - F(j+1,k  ,l-1) &
                                           +  F(j-1,k+1,l-1) - F(j-1,k  ,l-1)) 
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + alphaz*dtsdz * (F(j  ,k  ,l+1) - F(j  ,k  ,l  )) &
                            +  betaz*dtsdz * (F(j+1,k  ,l+1) - F(j+1,k  ,l  ) &
                                           +  F(j-1,k  ,l+1) - F(j-1,k  ,l  ) &
                                           +  F(j  ,k+1,l+1) - F(j  ,k+1,l  ) &
                                           +  F(j  ,k-1,l+1) - F(j  ,k-1,l  )) &
                            + gammaz*dtsdz * (F(j+1,k+1,l+1) - F(j+1,k+1,l  ) &
                                           +  F(j-1,k+1,l+1) - F(j-1,k+1,l  ) &
                                           +  F(j+1,k-1,l+1) - F(j+1,k-1,l  ) &
                                           +  F(j-1,k-1,l+1) - F(j-1,k-1,l  )) 
    end do
   end do
  end do

else

  k=0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + alphax*dtsdx * (F(j+1,k  ,l  ) - F(j  ,k  ,l  )) &
                            +  betax*dtsdx * (F(j+1,k  ,l+1) - F(j  ,k  ,l+1) &
                                           +  F(j+1,k  ,l-1) - F(j  ,k  ,l-1))  
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + alphaz*dtsdz * (F(j  ,k  ,l+1) - F(j  ,k  ,l  )) &
                            +  betaz*dtsdz * (F(j+1,k  ,l+1) - F(j+1,k  ,l  ) &
                                           +  F(j-1,k  ,l+1) - F(j-1,k  ,l  ))
    end do
  end do

end if

return
end subroutine push_em3d_kyeeefvec

subroutine push_em3d_phi(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,dtsepsi

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz

  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx
      f%Phi(j,k,l) = f%Phi(j,k,l) - dtsdx * (f%Ax(j,k,l) - f%Ax(j-1,k  ,l  )) &
                                  - dtsdy * (f%Ay(j,k,l) - f%Ay(j  ,k-1,l  )) &
                                  - dtsdz * (f%Az(j,k,l) - f%Az(j  ,k  ,l-1)) 
    end do
   end do
  end do

end subroutine push_em3d_phi

subroutine push_em3d_a(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

  dtsdx = f%clight*dt/f%dx
  dtsdy = f%clight*dt/f%dy
  dtsdz = f%clight*dt/f%dz

  ! advance Ex
  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx-1
      f%Ax(j,k,l) = f%Ax(j,k,l) - dtsdx * (f%Phi(j+1,k,l) - f%Phi(j,k,l)) - dt*f%Ex(j,k,l)
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      f%Ay(j,k,l) = f%Ay(j,k,l) - dtsdy * (f%Phi(j,k+1,l) - f%Phi(j,k,l)) - dt*f%Ey(j,k,l)
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      f%Az(j,k,l) = f%Az(j,k,l) - dtsdz * (f%Phi(j,k,l+1) - f%Phi(j,k,l)) - dt*f%Ez(j,k,l)
    end do
   end do
  end do

return
end subroutine push_em3d_a

subroutine push_em3d_conde(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

dtsdx = f%clight**2*dt/f%dx
dtsdy = f%clight**2*dt/f%dy
dtsdz = f%clight**2*dt/f%dz
mudt  = f%mu0*f%clight**2*dt


  ! advance Ex
  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx-1
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j+1,k,l)) &
      f%Ex(j,k,l) = f%Ex(j,k,l) + dtsdy * (f%Bz(j,k,l)   - f%Bz(j,k-1,l  )) &
                                - dtsdz * (f%By(j,k,l)   - f%By(j,k  ,l-1)) &
                                - mudt  * f%J(j,k,l,1)
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k+1,l)) &
      f%Ey(j,k,l) = f%Ey(j,k,l) - dtsdx * (f%Bz(j,k,l)   - f%Bz(j-1,k,l)) &
                                + dtsdz * (f%Bx(j,k,l)   - f%Bx(j,k,l-1)) &
                                - mudt  * f%J(j,k,l,2)
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k,l+1)) &
      f%Ez(j,k,l) = f%Ez(j,k,l) + dtsdx * (f%By(j,k,l) - f%By(j-1,k  ,l)) &
                                - dtsdy * (f%Bx(j,k,l) - f%Bx(j  ,k-1,l)) &
                                - mudt  * f%J(j,k,l,3)
    end do
   end do
  end do

return
end subroutine push_em3d_conde

subroutine push_em3d_condef(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

  dtsdx = f%clight*dt/f%dx
  dtsdy = f%clight*dt/f%dy
  dtsdz = f%clight*dt/f%dz

  ! advance Ex
  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx-1
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j+1,k,l)) &
      f%Ex(j,k,l) = f%Ex(j,k,l) + dtsdx * (f%F(j+1,k,l) - f%F(j,k,l)) 
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k+1,l)) &
      f%Ey(j,k,l) = f%Ey(j,k,l) + dtsdy * (f%F(j,k+1,l) - f%F(j,k,l))
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k,l+1)) &
      f%Ez(j,k,l) = f%Ez(j,k,l) + dtsdz * (f%F(j,k,l+1) - f%F(j,k,l)) 
    end do
   end do
  end do

return
end subroutine push_em3d_condef

subroutine push_em3d_condf(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,dtsepsi

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz
dtsepsi = f%mu0*f%clight**3*dt

  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx
      if (.not.f%incond(j,k,l)) &
      f%F(j,k,l) = f%F(j,k,l) + dtsdx * (f%Ex(j,k,l) - f%Ex(j-1,k  ,l  )) &
                              + dtsdy * (f%Ey(j,k,l) - f%Ey(j  ,k-1,l  )) &
                              + dtsdz * (f%Ez(j,k,l) - f%Ez(j  ,k  ,l-1)) &
                              - dtsepsi * f%Rho(j,k,l)
    end do
   end do
  end do

end subroutine push_em3d_condf

subroutine push_em3d_splite(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  if (sf%stencil==0 .or. sf%stencil==1) then
  call push_em3d_splitevec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                           sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                           sf%afx,sf%afy,sf%afz, &
                           sf%bpfx,sf%bpfy,sf%bpfz, &
                           sf%bmfx,sf%bmfy,sf%bmfz,sf%l_2dxz)
  else
    write(0,*) 'splite extended pml not implemented'
    stop
  end if
  
  return
end subroutine push_em3d_splite

subroutine push_em3d_splitevec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz,l_2dxz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx,exy,exz, &
                                                                                                       eyx,eyy,eyz, &
                                                                                                       ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      exy(j,k,l) = afy(k)*exy(j,k,l) + bpfy(k)*(bzx(j,k,l)+bzy(j,k,l))  &
                                     + bmfy(k)*(bzx(j,k-1,l)+bzy(j,k-1,l)) !- 0.5_8*dt*j(j,k,l,1)
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      exz(j,k,l) = afz(l)*exz(j,k,l) - bpfz(l)*(byx(j,k,l)+byz(j,k,l))  &
                                     - bmfz(l)*(byx(j,k,l-1)+byz(j,k,l-1)) !- 0.5_8*dt*j(j,k,l,1)
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      eyx(j,k,l) = afx(j)*eyx(j,k,l) - bpfx(j)*(bzx(j,k,l)+bzy(j,k,l))  &
                                     - bmfx(j)*(bzx(j-1,k,l)+bzy(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,2)
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      eyz(j,k,l) = afz(l)*eyz(j,k,l) + bpfz(l)*(bxy(j,k,l)+bxz(j,k,l))  &
                                     + bmfz(l)*(bxy(j,k,l-1)+bxz(j,k,l-1)) !- 0.5_8*dt*j(j,k,l,2)
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      ezx(j,k,l) = afx(j)*ezx(j,k,l) + bpfx(j)*(byx(j,k,l)+byz(j,k,l))  &
                                     + bmfx(j)*(byx(j-1,k,l)+byz(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,3)
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      ezy(j,k,l) = afy(k)*ezy(j,k,l) - bpfy(k)*(bxy(j,k,l)+bxz(j,k,l))  &
                                     - bmfy(k)*(bxy(j,k-1,l)+bxz(j,k-1,l)) !- 0.5_8*dt*j(j,k,l,3)
    end do
   end do
  end do

else
  k = 0

  do l = 0, nz
    do j = 0, nx-1
      exz(j,k,l) = afz(l)*exz(j,k,l) - bpfz(l)*(byx(j,k,l)+byz(j,k,l))  &
                                     - bmfz(l)*(byx(j,k,l-1)+byz(j,k,l-1)) !- 0.5_8*dt*j(j,k,l,1)
    end do
  end do

  do l = 0, nz
    do j = 0, nx
      eyx(j,k,l) = afx(j)*eyx(j,k,l) - bpfx(j)*(bzx(j,k,l))  &
                                     - bmfx(j)*(bzx(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,2)
    end do
  end do

  do l = 0, nz
    do j = 0, nx
      eyz(j,k,l) = afz(l)*eyz(j,k,l) + bpfz(l)*(bxz(j,k,l))  &
                                     + bmfz(l)*(bxz(j,k,l-1)) !- 0.5_8*dt*j(j,k,l,2)
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx
      ezx(j,k,l) = afx(j)*ezx(j,k,l) + bpfx(j)*(byx(j,k,l)+byz(j,k,l))  &
                                     + bmfx(j)*(byx(j-1,k,l)+byz(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,3)
    end do
  end do

end if

  return
end subroutine push_em3d_splitevec

subroutine push_em3d_splitef(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  if (sf%stencil==0 .or. sf%stencil==2) then
    call push_em3d_splitefvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%exx,sf%eyy,sf%ezz, &
                             sf%fx,sf%fy,sf%fz, &
                             sf%agx,sf%agy,sf%agz, &
                             sf%bpgx,sf%bpgy,sf%bpgz, &
                             sf%bmgx,sf%bmgy,sf%bmgz,sf%l_2dxz)
  else
    call push_em3d_splitkyeeefvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%exx,sf%eyy,sf%ezz, &
                             sf%fx,sf%fy,sf%fz, &
                             sf%agx,sf%agy,sf%agz, &
                             sf%bpgx,sf%bpgy,sf%bpgz, &
                             sf%bmgx,sf%bmgy,sf%bmgz,sf%l_2dxz)
  end if

  return
end subroutine push_em3d_splitef

subroutine push_em3d_splitefvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,eyy,ezz,fx,fy,fz, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz,l_2dxz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx,eyy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      exx(j,k,l) = agx(j)*exx(j,k,l) + bpgx(j)*( fx(j+1,k,l) + fy(j+1,k,l) + fz(j+1,k,l) ) &
                                     + bmgx(j)*( fx(j  ,k,l) + fy(j  ,k,l) + fz(j  ,k,l) )
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      eyy(j,k,l) = agy(k)*eyy(j,k,l) + bpgy(k)*( fx(j,k+1,l) + fy(j,k+1,l) + fz(j,k+1,l) ) &
                                     + bmgy(k)*( fx(j,k  ,l) + fy(j,k  ,l) + fz(j,k  ,l) )
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      ezz(j,k,l) = agz(l)*ezz(j,k,l) + bpgz(l)*( fx(j,k,l+1) + fy(j,k,l+1) + fz(j,k,l+1) ) &
                                     + bmgz(l)*( fx(j,k,l)   + fy(j,k,l  ) + fz(j,k,l  ) )
    end do
   end do
  end do

else
  k = 0
  do l = 0, nz
    do j = 0, nx-1
      exx(j,k,l) = agx(j)*exx(j,k,l) + bpgx(j)*( fx(j+1,k,l) + fz(j+1,k,l) ) &
                                     + bmgx(j)*( fx(j  ,k,l) + fz(j  ,k,l) )
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx
      ezz(j,k,l) = agz(l)*ezz(j,k,l) + bpgz(l)*( fx(j,k,l+1) + fz(j,k,l+1) ) &
                                     + bmgz(l)*( fx(j,k,l)   + fz(j,k,l  ) )
    end do
  end do

end if

  return
end subroutine push_em3d_splitefvec

subroutine push_em3d_splitkyeeefvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,eyy,ezz,fx,fy,fz, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz,l_2dxz)
use EM3D_kyee
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx,eyy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      exx(j,k,l) = agx(j)*exx(j,k,l) + bpgx(j)*alphax * ( fx(j+1,k  ,l  ) + fy(j+1,k  ,l  ) + fz(j+1,k  ,l  ) ) &
                                     + bpgx(j)* betax * ( fx(j+1,k+1,l  ) + fy(j+1,k+1,l  ) + fz(j+1,k+1,l  )  &
                                                      +   fx(j+1,k-1,l  ) + fy(j+1,k-1,l  ) + fz(j+1,k-1,l  )  &
                                                      +   fx(j+1,k  ,l+1) + fy(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                      +   fx(j+1,k  ,l-1) + fy(j+1,k  ,l-1) + fz(j+1,k  ,l-1))  &
                                     + bpgx(j)*gammax * ( fx(j+1,k+1,l+1) + fy(j+1,k+1,l+1) + fz(j+1,k+1,l+1)  &
                                                      +   fx(j+1,k-1,l+1) + fy(j+1,k-1,l+1) + fz(j+1,k-1,l+1)  &
                                                      +   fx(j+1,k+1,l-1) + fy(j+1,k+1,l-1) + fz(j+1,k+1,l-1)  &
                                                      +   fx(j+1,k-1,l-1) + fy(j+1,k-1,l-1) + fz(j+1,k-1,l-1)) &
                                     + bmgx(j)*alphax * ( fx(j  ,k  ,l  ) + fy(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgx(j)* betax * ( fx(j  ,k+1,l  ) + fy(j  ,k+1,l  ) + fz(j  ,k+1,l  )  &
                                                      +   fx(j  ,k-1,l  ) + fy(j  ,k-1,l  ) + fz(j  ,k-1,l  )  &
                                                      +   fx(j  ,k  ,l+1) + fy(j  ,k  ,l+1) + fz(j  ,k  ,l+1)  &
                                                      +   fx(j  ,k  ,l-1) + fy(j  ,k  ,l-1) + fz(j  ,k  ,l-1))  &
                                     + bmgx(j)*gammax * ( fx(j  ,k+1,l+1) + fy(j  ,k+1,l+1) + fz(j  ,k+1,l+1)  &
                                                      +   fx(j  ,k-1,l+1) + fy(j  ,k-1,l+1) + fz(j  ,k-1,l+1)  &
                                                      +   fx(j  ,k+1,l-1) + fy(j  ,k+1,l-1) + fz(j  ,k+1,l-1)  &
                                                      +   fx(j  ,k-1,l-1) + fy(j  ,k-1,l-1) + fz(j  ,k-1,l-1))
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      eyy(j,k,l) = agy(k)*eyy(j,k,l) + bpgy(k)*alphay * ( fx(j  ,k+1,l  ) + fy(j  ,k+1,l  ) + fz(j  ,k+1,l  ) ) &
                                     + bpgy(k)* betay * ( fx(j+1,k+1,l  ) + fy(j+1,k+1,l  ) + fz(j+1,k+1,l  ) &
                                                      +   fx(j-1,k+1,l  ) + fy(j-1,k+1,l  ) + fz(j-1,k+1,l  ) &
                                                      +   fx(j  ,k+1,l+1) + fy(j  ,k+1,l+1) + fz(j  ,k+1,l+1) &
                                                      +   fx(j  ,k+1,l-1) + fy(j  ,k+1,l-1) + fz(j  ,k+1,l-1)) &
                                     + bpgy(k)*gammay * ( fx(j+1,k+1,l+1) + fy(j+1,k+1,l+1) + fz(j+1,k+1,l+1) &
                                                      +   fx(j-1,k+1,l+1) + fy(j-1,k+1,l+1) + fz(j-1,k+1,l+1) &
                                                      +   fx(j+1,k+1,l-1) + fy(j+1,k+1,l-1) + fz(j+1,k+1,l-1) &
                                                      +   fx(j-1,k+1,l-1) + fy(j-1,k+1,l-1) + fz(j-1,k+1,l-1)) &
                                     + bmgy(k)*alphay * ( fx(j  ,k  ,l  ) + fy(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgy(k)* betay * ( fx(j+1,k  ,l  ) + fy(j+1,k  ,l  ) + fz(j+1,k  ,l  ) &
                                                      +   fx(j-1,k  ,l  ) + fy(j-1,k  ,l  ) + fz(j-1,k  ,l  ) &
                                                      +   fx(j  ,k  ,l+1) + fy(j  ,k  ,l+1) + fz(j  ,k  ,l+1) &
                                                      +   fx(j  ,k  ,l-1) + fy(j  ,k  ,l-1) + fz(j  ,k  ,l-1)) &
                                     + bmgy(k)*gammay * ( fx(j+1,k  ,l+1) + fy(j+1,k  ,l+1) + fz(j+1,k  ,l+1) &
                                                      +   fx(j-1,k  ,l+1) + fy(j-1,k  ,l+1) + fz(j-1,k  ,l+1) &
                                                      +   fx(j+1,k  ,l-1) + fy(j+1,k  ,l-1) + fz(j+1,k  ,l-1) &
                                                      +   fx(j-1,k  ,l-1) + fy(j-1,k  ,l-1) + fz(j-1,k  ,l-1)) 
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      ezz(j,k,l) = agz(l)*ezz(j,k,l) + bpgz(l)*alphaz * ( fx(j  ,k  ,l+1) + fy(j  ,k  ,l+1) + fz(j  ,k  ,l+1) ) &
                                     + bpgz(l)* betaz * ( fx(j+1,k  ,l+1) + fy(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                      +   fx(j-1,k  ,l+1) + fy(j-1,k  ,l+1) + fz(j-1,k  ,l+1)  &
                                                      +   fx(j  ,k+1,l+1) + fy(j  ,k+1,l+1) + fz(j  ,k+1,l+1)  &
                                                      +   fx(j  ,k-1,l+1) + fy(j  ,k-1,l+1) + fz(j  ,k-1,l+1))  &
                                     + bpgz(l)*gammaz * ( fx(j+1,k+1,l+1) + fy(j+1,k+1,l+1) + fz(j+1,k+1,l+1)  &
                                                      +   fx(j-1,k+1,l+1) + fy(j-1,k+1,l+1) + fz(j-1,k+1,l+1)  &
                                                      +   fx(j+1,k-1,l+1) + fy(j+1,k-1,l+1) + fz(j+1,k-1,l+1)  &
                                                      +   fx(j-1,k-1,l+1) + fy(j-1,k-1,l+1) + fz(j-1,k-1,l+1))  &
                                     + bmgz(l)*alphaz * ( fx(j  ,k  ,l  ) + fy(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgz(l)* betaz * ( fx(j+1,k  ,l  ) + fy(j+1,k  ,l  ) + fz(j+1,k  ,l  )  &
                                                      +   fx(j-1,k  ,l  ) + fy(j-1,k  ,l  ) + fz(j-1,k  ,l  )  &
                                                      +   fx(j  ,k+1,l  ) + fy(j  ,k+1,l  ) + fz(j  ,k+1,l  )  &
                                                      +   fx(j  ,k-1,l  ) + fy(j  ,k-1,l  ) + fz(j  ,k-1,l  ))  &
                                     + bmgz(l)*gammaz * ( fx(j+1,k+1,l  ) + fy(j+1,k+1,l  ) + fz(j+1,k+1,l  )  &
                                                      +   fx(j-1,k+1,l  ) + fy(j-1,k+1,l  ) + fz(j-1,k+1,l  )  &
                                                      +   fx(j+1,k-1,l  ) + fy(j+1,k-1,l  ) + fz(j+1,k-1,l  )  &
                                                      +   fx(j-1,k-1,l  ) + fy(j-1,k-1,l  ) + fz(j-1,k-1,l  ))  
    end do
   end do
  end do

else
  k = 0

  do l = 0, nz
    do j = 0, nx-1
      exx(j,k,l) = agx(j)*exx(j,k,l) + bpgx(j)*    alphax * ( fx(j+1,k  ,l  ) + fz(j+1,k  ,l  ) ) &
                                     + bpgx(j)* 2.5*betax * ( fx(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                          +   fx(j+1,k  ,l-1) + fz(j+1,k  ,l-1))  &
                                     + bmgx(j)*    alphax * ( fx(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgx(j)* 2.5*betax * ( fx(j  ,k  ,l+1) + fz(j  ,k  ,l+1)  &
                                                          +   fx(j  ,k  ,l-1) + fz(j  ,k  ,l-1))  
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx
      ezz(j,k,l) = agz(l)*ezz(j,k,l) + bpgz(l)*    alphaz * ( fx(j  ,k  ,l+1) + fz(j  ,k  ,l+1) ) &
                                     + bpgz(l)* 2.5*betaz * ( fx(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                          +   fx(j-1,k  ,l+1) + fz(j-1,k  ,l+1))  &
                                     + bmgz(l)*    alphaz * ( fx(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgz(l)* 2.5*betaz * ( fx(j+1,k  ,l  ) + fz(j+1,k  ,l  )  &
                                                          +   fx(j-1,k  ,l  ) + fz(j-1,k  ,l  ))  
    end do
  end do

end if

  return
end subroutine push_em3d_splitkyeeefvec

subroutine push_em3d_splitb(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  if (sf%stencil==0 .or. sf%stencil==2) then
    call push_em3d_splitbvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                             sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                             sf%agx,sf%agy,sf%agz, &
                             sf%bpgx,sf%bpgy,sf%bpgz, &
                             sf%bmgx,sf%bmgy,sf%bmgz,sf%l_2dxz)
  else
    call push_em3d_splitkyeebvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                             sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                             sf%agx,sf%agy,sf%agz, &
                             sf%bpgx,sf%bpgy,sf%bpgz, &
                             sf%bmgx,sf%bmgy,sf%bmgz,sf%l_2dxz)
  end if

  return
end subroutine push_em3d_splitb

subroutine push_em3d_splitbvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz,l_2dxz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz
!real(kind=8), dimension(-nxguard:,-nyguard:,-nzguard:), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
!real(kind=8), dimension(-nxguard:,-nyguard:,-nzguard:), intent(in) :: exx,exy,exz, &
!                                                                                                    eyx,eyy,eyz, &
!                                                                                                    ezx,ezy,ezz
!real(kind=8), dimension(-nxguard:), intent(in) :: agx,bpgx,bmgx
!real(kind=8), dimension(-nyguard:), intent(in) :: agy,bpgy,bmgy
!real(kind=8), dimension(-nzguard:), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxy(j,k,l) = agy(k)*bxy(j,k,l) - bpgy(k)*(ezx(j,k+1,l  )+ezy(j,k+1,l  )+ezz(j,k+1,l  )) &
                                     - bmgy(k)*(ezx(j,k  ,l  )+ezy(j,k  ,l  )+ezz(j,k  ,l  ))
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*(eyx(j,k  ,l+1)+eyy(j,k  ,l+1)+eyz(j,k  ,l+1)) &
                                     + bmgz(l)*(eyx(j,k  ,l  )+eyy(j,k  ,l  )+eyz(j,k  ,l  ))
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*(ezx(j+1,k,l  )+ezy(j+1,k,l  )+ezz(j+1,k,l  )) &
                                     + bmgx(j)*(ezx(j  ,k,l  )+ezy(j  ,k,l  )+ezz(j  ,k,l  ))
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*(exx(j  ,k,l+1)+exy(j  ,k,l+1)+exz(j  ,k,l+1)) &
                                     - bmgz(l)*(exx(j  ,k,l  )+exy(j  ,k,l  )+exz(j  ,k,l  ))
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*(eyx(j+1,k  ,l)+eyy(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                     - bmgx(j)*(eyx(j  ,k  ,l)+eyy(j  ,k  ,l)+eyz(j  ,k  ,l))
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      bzy(j,k,l) = agy(k)*bzy(j,k,l) + bpgy(k)*(exx(j  ,k+1,l)+exy(j  ,k+1,l)+exz(j  ,k+1,l)) &
                                     + bmgy(k)*(exx(j  ,k  ,l)+exy(j  ,k  ,l)+exz(j  ,k  ,l))
    end do
   end do
  end do

else
  k = 0
  do l = 0, nz-1
    do j = 0, nx
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*(eyx(j,k  ,l+1)+eyz(j,k  ,l+1)) &
                                     + bmgz(l)*(eyx(j,k  ,l  )+eyz(j,k  ,l  ))
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*(ezx(j+1,k,l  )+ezz(j+1,k,l  )) &
                                     + bmgx(j)*(ezx(j  ,k,l  )+ezz(j  ,k,l  ))
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*(exx(j  ,k,l+1)+exz(j  ,k,l+1)) &
                                     - bmgz(l)*(exx(j  ,k,l  )+exz(j  ,k,l  ))
    end do
  end do

  do l = 0, nz
    do j = 0, nx-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*(eyx(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                     - bmgx(j)*(eyx(j  ,k  ,l)+eyz(j  ,k  ,l))
    end do
  end do

end if

  return
end subroutine push_em3d_splitbvec


subroutine push_em3d_splitkyeebvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz,l_2dxz)
use EM3D_kyee
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxy(j,k,l) = agy(k)*bxy(j,k,l) - bpgy(k)*alphay * (ezx(j  ,k+1,l  )+ezy(j  ,k+1,l  )+ezz(j  ,k+1,l  )) &
                                     - bpgy(k)* betay * (ezx(j+1,k+1,l  )+ezy(j+1,k+1,l  )+ezz(j+1,k+1,l  ) &
                                                      +  ezx(j-1,k+1,l  )+ezy(j-1,k+1,l  )+ezz(j-1,k+1,l  ) &
                                                      +  ezx(j  ,k+1,l+1)+ezy(j  ,k+1,l+1)+ezz(j  ,k+1,l+1) &
                                                      +  ezx(j  ,k+1,l-1)+ezy(j  ,k+1,l-1)+ezz(j  ,k+1,l-1)) &
                                     - bpgy(k)*gammay * (ezx(j+1,k+1,l+1)+ezy(j+1,k+1,l+1)+ezz(j+1,k+1,l+1) &
                                                      +  ezx(j-1,k+1,l+1)+ezy(j-1,k+1,l+1)+ezz(j-1,k+1,l+1) &
                                                      +  ezx(j+1,k+1,l-1)+ezy(j+1,k+1,l-1)+ezz(j+1,k+1,l-1) &
                                                      +  ezx(j-1,k+1,l-1)+ezy(j-1,k+1,l-1)+ezz(j-1,k+1,l-1)) &
                                     - bmgy(k)*alphay * (ezx(j  ,k  ,l  )+ezy(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     - bmgy(k)* betay * (ezx(j+1,k  ,l  )+ezy(j+1,k  ,l  )+ezz(j+1,k  ,l  ) &
                                                      +  ezx(j-1,k  ,l  )+ezy(j-1,k  ,l  )+ezz(j-1,k  ,l  ) &
                                                      +  ezx(j  ,k  ,l+1)+ezy(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
                                                      +  ezx(j  ,k  ,l-1)+ezy(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) &
                                     - bmgy(k)*gammay * (ezx(j+1,k  ,l+1)+ezy(j+1,k  ,l+1)+ezz(j+1,k  ,l+1) &
                                                      +  ezx(j-1,k  ,l+1)+ezy(j-1,k  ,l+1)+ezz(j-1,k  ,l+1) &
                                                      +  ezx(j+1,k  ,l-1)+ezy(j+1,k  ,l-1)+ezz(j+1,k  ,l-1) &
                                                      +  ezx(j-1,k  ,l-1)+ezy(j-1,k  ,l-1)+ezz(j-1,k  ,l-1)) 
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*alphaz * (eyx(j  ,k  ,l+1)+eyy(j  ,k  ,l+1)+eyz(j  ,k  ,l+1)) &
                                     + bpgz(l)* betaz * (eyx(j+1,k  ,l+1)+eyy(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                      +  eyx(j-1,k  ,l+1)+eyy(j-1,k  ,l+1)+eyz(j-1,k  ,l+1) &
                                                      +  eyx(j  ,k+1,l+1)+eyy(j  ,k+1,l+1)+eyz(j  ,k+1,l+1) &
                                                      +  eyx(j  ,k-1,l+1)+eyy(j  ,k-1,l+1)+eyz(j  ,k-1,l+1)) &
                                     + bpgz(l)*gammaz * (eyx(j+1,k+1,l+1)+eyy(j+1,k+1,l+1)+eyz(j+1,k+1,l+1) &
                                                      +  eyx(j-1,k+1,l+1)+eyy(j-1,k+1,l+1)+eyz(j-1,k+1,l+1) &
                                                      +  eyx(j+1,k-1,l+1)+eyy(j+1,k-1,l+1)+eyz(j+1,k-1,l+1) &
                                                      +  eyx(j-1,k-1,l+1)+eyy(j-1,k-1,l+1)+eyz(j-1,k-1,l+1)) &
                                     + bmgz(l)*alphaz * (eyx(j  ,k  ,l  )+eyy(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     + bmgz(l)* betaz * (eyx(j+1,k  ,l  )+eyy(j+1,k  ,l  )+eyz(j+1,k  ,l  ) &
                                                      +  eyx(j-1,k  ,l  )+eyy(j-1,k  ,l  )+eyz(j-1,k  ,l  ) &
                                                      +  eyx(j  ,k+1,l  )+eyy(j  ,k+1,l  )+eyz(j  ,k+1,l  ) &
                                                      +  eyx(j  ,k-1,l  )+eyy(j  ,k-1,l  )+eyz(j  ,k-1,l  )) &
                                     + bmgz(l)*gammaz * (eyx(j+1,k+1,l  )+eyy(j+1,k+1,l  )+eyz(j+1,k+1,l  ) &
                                                      +  eyx(j-1,k+1,l  )+eyy(j-1,k+1,l  )+eyz(j-1,k+1,l  ) &
                                                      +  eyx(j+1,k-1,l  )+eyy(j+1,k-1,l  )+eyz(j+1,k-1,l  ) &
                                                      +  eyx(j-1,k-1,l  )+eyy(j-1,k-1,l  )+eyz(j-1,k-1,l  )) 
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*alphax * (ezx(j+1,k  ,l  )+ezy(j+1,k  ,l  )+ezz(j+1,k  ,l  )) &
                                     + bpgx(j)* betax * (ezx(j+1,k+1,l  )+ezy(j+1,k+1,l  )+ezz(j+1,k+1,l  ) &
                                                      +  ezx(j+1,k-1,l  )+ezy(j+1,k-1,l  )+ezz(j+1,k-1,l  ) &
                                                      +  ezx(j+1,k  ,l+1)+ezy(j+1,k  ,l+1)+ezz(j+1,k  ,l+1) &
                                                      +  ezx(j+1,k  ,l-1)+ezy(j+1,k  ,l-1)+ezz(j+1,k  ,l-1)) &
                                     + bpgx(j)*gammax * (ezx(j+1,k+1,l+1)+ezy(j+1,k+1,l+1)+ezz(j+1,k+1,l+1) &
                                                      +  ezx(j+1,k-1,l+1)+ezy(j+1,k-1,l+1)+ezz(j+1,k-1,l+1) &
                                                      +  ezx(j+1,k+1,l-1)+ezy(j+1,k+1,l-1)+ezz(j+1,k+1,l-1) &
                                                      +  ezx(j+1,k-1,l-1)+ezy(j+1,k-1,l-1)+ezz(j+1,k-1,l-1)) &
                                     + bmgx(j)*alphax * (ezx(j  ,k  ,l  )+ezy(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     + bmgx(j)* betax * (ezx(j  ,k+1,l  )+ezy(j  ,k+1,l  )+ezz(j  ,k+1,l  ) &
                                                      +  ezx(j  ,k-1,l  )+ezy(j  ,k-1,l  )+ezz(j  ,k-1,l  ) &
                                                      +  ezx(j  ,k  ,l+1)+ezy(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
                                                      +  ezx(j  ,k  ,l-1)+ezy(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) &
                                     + bmgx(j)*gammax * (ezx(j  ,k+1,l+1)+ezy(j  ,k+1,l+1)+ezz(j  ,k+1,l+1) &
                                                      +  ezx(j  ,k-1,l+1)+ezy(j  ,k-1,l+1)+ezz(j  ,k-1,l+1) &
                                                      +  ezx(j  ,k+1,l-1)+ezy(j  ,k+1,l-1)+ezz(j  ,k+1,l-1) &
                                                      +  ezx(j  ,k-1,l-1)+ezy(j  ,k-1,l-1)+ezz(j  ,k-1,l-1)) 
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*alphaz * (exx(j  ,k  ,l+1)+exy(j  ,k  ,l+1)+exz(j  ,k  ,l+1)) &
                                     - bpgz(l)* betaz * (exx(j+1,k  ,l+1)+exy(j+1,k  ,l+1)+exz(j+1,k  ,l+1) &
                                                      +  exx(j-1,k  ,l+1)+exy(j-1,k  ,l+1)+exz(j-1,k  ,l+1) &
                                                      +  exx(j  ,k+1,l+1)+exy(j  ,k+1,l+1)+exz(j  ,k+1,l+1) &
                                                      +  exx(j  ,k-1,l+1)+exy(j  ,k-1,l+1)+exz(j  ,k-1,l+1)) &
                                     - bpgz(l)*gammaz * (exx(j+1,k+1,l+1)+exy(j+1,k+1,l+1)+exz(j+1,k+1,l+1) &
                                                      +  exx(j-1,k+1,l+1)+exy(j-1,k+1,l+1)+exz(j-1,k+1,l+1) &
                                                      +  exx(j+1,k-1,l+1)+exy(j+1,k-1,l+1)+exz(j+1,k-1,l+1) &
                                                      +  exx(j-1,k-1,l+1)+exy(j-1,k-1,l+1)+exz(j-1,k-1,l+1)) &
                                     - bmgz(l)*alphaz * (exx(j  ,k  ,l  )+exy(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     - bmgz(l)* betaz * (exx(j+1,k  ,l  )+exy(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                      +  exx(j-1,k  ,l  )+exy(j-1,k  ,l  )+exz(j-1,k  ,l  ) &
                                                      +  exx(j  ,k+1,l  )+exy(j  ,k+1,l  )+exz(j  ,k+1,l  ) &
                                                      +  exx(j  ,k-1,l  )+exy(j  ,k-1,l  )+exz(j  ,k-1,l  )) &
                                     - bmgz(l)*gammaz * (exx(j+1,k+1,l  )+exy(j+1,k+1,l  )+exz(j+1,k+1,l  ) &
                                                      +  exx(j-1,k+1,l  )+exy(j-1,k+1,l  )+exz(j-1,k+1,l  ) &
                                                      +  exx(j+1,k-1,l  )+exy(j+1,k-1,l  )+exz(j+1,k-1,l  ) &
                                                      +  exx(j-1,k-1,l  )+exy(j-1,k-1,l  )+exz(j-1,k-1,l  )) 

    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*alphax * (eyx(j+1,k  ,l  )+eyy(j+1,k  ,l  )+eyz(j+1,k  ,l  )) &
                                     - bpgx(j)* betax * (eyx(j+1,k+1,l  )+eyy(j+1,k+1,l  )+eyz(j+1,k+1,l  ) &
                                                      +  eyx(j+1,k-1,l  )+eyy(j+1,k-1,l  )+eyz(j+1,k-1,l  ) &
                                                      +  eyx(j+1,k  ,l+1)+eyy(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                      +  eyx(j+1,k  ,l-1)+eyy(j+1,k  ,l-1)+eyz(j+1,k  ,l-1)) &
                                     - bpgx(j)*gammax * (eyx(j+1,k+1,l+1)+eyy(j+1,k+1,l+1)+eyz(j+1,k+1,l+1) &
                                                      +  eyx(j+1,k-1,l+1)+eyy(j+1,k-1,l+1)+eyz(j+1,k-1,l+1) &
                                                      +  eyx(j+1,k+1,l-1)+eyy(j+1,k+1,l-1)+eyz(j+1,k+1,l-1) &
                                                      +  eyx(j+1,k-1,l-1)+eyy(j+1,k-1,l-1)+eyz(j+1,k-1,l-1)) &
                                     - bmgx(j)*alphax * (eyx(j  ,k  ,l  )+eyy(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     - bmgx(j)* betax * (eyx(j  ,k+1,l  )+eyy(j  ,k+1,l  )+eyz(j  ,k+1,l  ) &
                                                      +  eyx(j  ,k-1,l  )+eyy(j  ,k-1,l  )+eyz(j  ,k-1,l  ) &
                                                      +  eyx(j  ,k  ,l+1)+eyy(j  ,k  ,l+1)+eyz(j  ,k  ,l+1) &
                                                      +  eyx(j  ,k  ,l-1)+eyy(j  ,k  ,l-1)+eyz(j  ,k  ,l-1)) &
                                     - bmgx(j)*gammax * (eyx(j  ,k+1,l+1)+eyy(j  ,k+1,l+1)+eyz(j  ,k+1,l+1) &
                                                      +  eyx(j  ,k-1,l+1)+eyy(j  ,k-1,l+1)+eyz(j  ,k-1,l+1) &
                                                      +  eyx(j  ,k+1,l-1)+eyy(j  ,k+1,l-1)+eyz(j  ,k+1,l-1) &
                                                      +  eyx(j  ,k-1,l-1)+eyy(j  ,k-1,l-1)+eyz(j  ,k-1,l-1)) 
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      bzy(j,k,l) = agy(k)*bzy(j,k,l) + bpgy(k)*alphay * (exx(j  ,k+1,l  )+exy(j  ,k+1,l  )+exz(j  ,k+1,l  )) &
                                     + bpgy(k)* betay * (exx(j+1,k+1,l  )+exy(j+1,k+1,l  )+exz(j+1,k+1,l  ) &
                                                      +  exx(j-1,k+1,l  )+exy(j-1,k+1,l  )+exz(j-1,k+1,l  ) &
                                                      +  exx(j  ,k+1,l+1)+exy(j  ,k+1,l+1)+exz(j  ,k+1,l+1) &
                                                      +  exx(j  ,k+1,l-1)+exy(j  ,k+1,l-1)+exz(j  ,k+1,l-1)) &
                                     + bpgy(k)*gammay * (exx(j+1,k+1,l+1)+exy(j+1,k+1,l+1)+exz(j+1,k+1,l+1) &
                                                      +  exx(j-1,k+1,l+1)+exy(j-1,k+1,l+1)+exz(j-1,k+1,l+1) &
                                                      +  exx(j+1,k+1,l-1)+exy(j+1,k+1,l-1)+exz(j+1,k+1,l-1) &
                                                      +  exx(j-1,k+1,l-1)+exy(j-1,k+1,l-1)+exz(j-1,k+1,l-1)) &
                                     + bmgy(k)*alphay * (exx(j  ,k  ,l  )+exy(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     + bmgy(k)* betay * (exx(j+1,k  ,l  )+exy(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                      +  exx(j-1,k  ,l  )+exy(j-1,k  ,l  )+exz(j-1,k  ,l  ) &
                                                      +  exx(j  ,k  ,l+1)+exy(j  ,k  ,l+1)+exz(j  ,k  ,l+1) &
                                                      +  exx(j  ,k  ,l-1)+exy(j  ,k  ,l-1)+exz(j  ,k  ,l-1)) &
                                     + bmgy(k)*gammay * (exx(j+1,k  ,l+1)+exy(j+1,k  ,l+1)+exz(j+1,k  ,l+1) &
                                                      +  exx(j-1,k  ,l+1)+exy(j-1,k  ,l+1)+exz(j-1,k  ,l+1) &
                                                      +  exx(j+1,k  ,l-1)+exy(j+1,k  ,l-1)+exz(j+1,k  ,l-1) &
                                                      +  exx(j-1,k  ,l-1)+exy(j-1,k  ,l-1)+exz(j-1,k  ,l-1)) 
    end do
   end do
  end do

else
  k = 0

  do l = 0, nz-1
    do j = 0, nx
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*alphaz * (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1)) &
                                     + bpgz(l)* betaz * (eyx(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                      +  eyx(j-1,k  ,l+1)+eyz(j-1,k  ,l+1)) &
                                     + bmgz(l)*alphaz * (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     + bmgz(l)* betaz * (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  ) &
                                                      +  eyx(j-1,k  ,l  )+eyz(j-1,k  ,l  )) 
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*    alphax * (ezx(j+1,k  ,l  )+ezz(j+1,k  ,l  )) &
                                     + bpgx(j)* 2.5*betax * (ezx(j+1,k  ,l+1)+ezz(j+1,k  ,l+1) &
                                                          +  ezx(j+1,k  ,l-1)+ezz(j+1,k  ,l-1)) &
                                     + bmgx(j)*    alphax * (ezx(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     + bmgx(j)* 2.5*betax * (ezx(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
                                                          +  ezx(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) 
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*    alphaz * (exx(j  ,k  ,l+1)+exz(j  ,k  ,l+1)) &
                                     - bpgz(l)* 2.5*betaz * (exx(j+1,k  ,l+1)+exz(j+1,k  ,l+1) &
                                                          +  exx(j-1,k  ,l+1)+exz(j-1,k  ,l+1)) &
                                     - bmgz(l)*    alphaz * (exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     - bmgz(l)* 2.5*betaz * (exx(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                          +  exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) 

    end do
  end do

  do l = 0, nz
    do j = 0, nx-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*    alphax * (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  )) &
                                     - bpgx(j)* 2.5*betax * (eyx(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                          +  eyx(j+1,k  ,l-1)+eyz(j+1,k  ,l-1)) &
                                     - bmgx(j)*    alphax * (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     - bmgx(j)* 2.5*betax * (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1) &
                                                          +  eyx(j  ,k  ,l-1)+eyz(j  ,k  ,l-1)) 
    end do
  end do

end if

  return
end subroutine push_em3d_splitkyeebvec


subroutine push_em3d_splitf(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  call push_em3d_splitfvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                           sf%fx,sf%fy,sf%fz, &
                           sf%afx,sf%afy,sf%afz, &
                           sf%bpfx,sf%bpfy,sf%bpfz, &
                           sf%bmfx,sf%bmfy,sf%bmfz,sf%l_2dxz)

  return
end subroutine push_em3d_splitf

subroutine push_em3d_splitfvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,fx,fy,fz, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz,l_2dxz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
     fx(j,k,l) = afx(j)*fx(j,k,l) + bpfx(j)*(exx(j  ,k  ,l  )+exy(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                  + bmfx(j)*(exx(j-1,k  ,l  )+exy(j-1,k  ,l  )+exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
     fy(j,k,l) = afy(k)*fy(j,k,l) + bpfy(k)*(eyx(j  ,k  ,l  )+eyy(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                  + bmfy(k)*(eyx(j  ,k-1,l  )+eyy(j  ,k-1,l  )+eyz(j  ,k-1,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
     fz(j,k,l) = afz(l)*fz(j,k,l) + bpfz(l)*(ezx(j  ,k  ,l  )+ezy(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                  + bmfz(l)*(ezx(j  ,k  ,l-1)+ezy(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
   end do
  end do

else
  k = 0
  do l = 0, nz
    do j = 0, nx
     fx(j,k,l) = afx(j)*fx(j,k,l) + bpfx(j)*(exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                  + bmfx(j)*(exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
  end do

  do l = 0, nz
    do j = 0, nx
     fz(j,k,l) = afz(l)*fz(j,k,l) + bpfz(l)*(ezx(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                  + bmfz(l)*(ezx(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
  end do

end if

  return
end subroutine push_em3d_splitfvec

subroutine push_em3d_splita(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  do l = 0, sf%nz
   do k = 0, sf%ny
    do j = 0, sf%nx-1
      sf%ax(j,k,l) = sf%agx(j)*sf%ax(j,k,l) - sf%bpgx(j)*( sf%phi(1,j+1,k,l) + sf%phi(2,j+1,k,l) + sf%phi(3,j+1,k,l) ) &
                                            - sf%bmgx(j)*( sf%phi(1,j  ,k,l) + sf%phi(2,j  ,k,l) + sf%phi(3,j  ,k,l) )
    end do
   end do
  end do

  do l = 0, sf%nz
   do k = 0, sf%ny-1
    do j = 0, sf%nx
      sf%ay(j,k,l) = sf%agy(k)*sf%ay(j,k,l) - sf%bpgy(k)*( sf%phi(1,j,k+1,l) + sf%phi(2,j,k+1,l) + sf%phi(3,j,k+1,l) ) &
                                            - sf%bmgy(k)*( sf%phi(1,j,k  ,l) + sf%phi(2,j,k  ,l) + sf%phi(3,j,k  ,l) )
    end do
   end do
  end do

  do l = 0, sf%nz-1
   do k = 0, sf%ny
    do j = 0, sf%nx
      sf%az(j,k,l) = sf%agz(l)*sf%az(j,k,l) - sf%bpgz(l)*( sf%phi(1,j,k,l+1) + sf%phi(2,j,k,l+1) + sf%phi(3,j,k,l+1) ) &
                                            - sf%bmgz(l)*( sf%phi(1,j,k,l)   + sf%phi(2,j,k,l  ) + sf%phi(3,j,k,l  ) )
    end do
   end do
  end do

  return
end subroutine push_em3d_splita

subroutine push_em3d_splitphi(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  do l = 0, sf%nz
   do k = 0, sf%ny
    do j = 0, sf%nx
     sf%phi(1,j,k,l) = sf%afx(j)*sf%phi(1,j,k,l) - sf%bpfx(j)*sf%ax(j,k,l) - sf%bmfx(j)*sf%ax(j-1,k  ,l  ) 
     sf%phi(2,j,k,l) = sf%afy(k)*sf%phi(2,j,k,l) - sf%bpfy(k)*sf%ay(j,k,l) - sf%bmfy(k)*sf%ay(j  ,k-1,l  ) 
     sf%phi(3,j,k,l) = sf%afz(l)*sf%phi(3,j,k,l) - sf%bpfz(l)*sf%az(j,k,l) - sf%bmfz(l)*sf%az(j  ,k  ,l-1) 
    end do
   end do
  end do

  return
end subroutine push_em3d_splitphi

subroutine push_em3d_block(b,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l

  call push_em3d_f(b%core%yf,dt*0.5)
  call push_em3d_b(b%core%yf,dt*0.5)
  call push_em3d_blockbndf(b,dt,2)
  call push_em3d_blockbndb(b,dt,2)

  call em3d_exchange_f(b)
  call em3d_exchange_b(b)

  call push_em3d_ef(b%core%yf,dt)
  call push_em3d_e(b%core%yf,dt)
  call push_em3d_blockbndef(b,dt,0)
  call push_em3d_blockbnde(b,dt,0)
  
  call em3d_exchange_e(b)

  call push_em3d_f(b%core%yf,dt*0.5)
  call push_em3d_b(b%core%yf,dt*0.5)
  call push_em3d_blockbndf(b,dt,1)
  call push_em3d_blockbndb(b,dt,1)

  call em3d_exchange_f(b)
  call em3d_exchange_b(b)

  return
end subroutine push_em3d_block

subroutine push_em3d_eef(b,dt,which,l_pushf,l_pushpot)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
logical(ISZ) :: l_pushf,l_pushpot

INTEGER(ISZ) :: j, k, l,which
  
  if (which==0) then
    if(l_pushpot) call push_em3d_phi(b%core%yf,dt)
    if(l_pushf) call push_em3d_ef(b%core%yf,dt)
    call push_em3d_e(b%core%yf,dt)
  else
    if(l_pushpot) call push_em3d_phi(b%core%yf,dt*0.5)
    if(l_pushf) call push_em3d_ef(b%core%yf,dt*0.5)
    call push_em3d_e(b%core%yf,dt*0.5)
  end if

!  if(l_pushpot) call push_em3d_blockbndphi(b,dt,which)
  if(l_pushf) call push_em3d_blockbndef(b,dt,which)
  call push_em3d_blockbnde(b,dt,which)
  
  call em3d_applybc_e(b%core%yf, &
                      b%xlbnd, &
                      b%xrbnd, &
                      b%ylbnd, &
                      b%yrbnd, &
                      b%zlbnd, &
                      b%zrbnd)
  ! --- need to exchange e even if not pushing f, for calculation of e at nodes
  call em3d_exchange_e(b)

  return
end subroutine push_em3d_eef

subroutine push_em3d_bf(b,dt,which,l_pushf,l_pushpot)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
logical(ISZ) :: l_pushf,l_pushpot

INTEGER(ISZ) :: j, k, l,which

  if (which==0) then
    if(l_pushf) call push_em3d_f(b%core%yf,dt)
    call push_em3d_b(b%core%yf,dt)
  else
    if(l_pushf) call push_em3d_f(b%core%yf,dt*0.5)
    call push_em3d_b(b%core%yf,dt*0.5)
  endif

  if(l_pushf) call push_em3d_blockbndf(b,dt,which)
  call push_em3d_blockbndb(b,dt,which)

  if(l_pushf) call em3d_exchange_f(b)
  call em3d_applybc_b(b%core%yf, &
                      b%xlbnd, &
                      b%xrbnd, &
                      b%ylbnd, &
                      b%yrbnd, &
                      b%zlbnd, &
                      b%zrbnd)
  call em3d_exchange_b(b)

  return
end subroutine push_em3d_bf

subroutine push_em3d_blockbnde(b,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
integer(ISZ) :: which

  

  if(b%xlbnd==openbc) call push_em3d_splite(b%sidexl%syf,dt,which)
  if(b%xrbnd==openbc) call push_em3d_splite(b%sidexr%syf,dt,which)
  if(b%ylbnd==openbc) call push_em3d_splite(b%sideyl%syf,dt,which)
  if(b%yrbnd==openbc) call push_em3d_splite(b%sideyr%syf,dt,which)
  if(b%zlbnd==openbc) call push_em3d_splite(b%sidezl%syf,dt,which)
  if(b%zrbnd==openbc) call push_em3d_splite(b%sidezr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splite(b%edgexlyl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splite(b%edgexryl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splite(b%edgexlyr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splite(b%edgexryr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%edgexlzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%edgexrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%edgexlzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%edgexrzr%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%edgeylzl%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%edgeyrzl%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%edgeylzr%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%edgeyrzr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%cornerxlylzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%cornerxrylzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%cornerxlyrzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splite(b%cornerxryrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%cornerxlylzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%cornerxrylzr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%cornerxlyrzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splite(b%cornerxryrzr%syf,dt,which)


  if (associated(b%sidexl%syf) .and. b%sidexl%proc==my_index) &
  call em3d_applybc_splite(b%sidexl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sidexr%syf) .and. b%sidexr%proc==my_index) &
  call em3d_applybc_splite(b%sidexr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sideyl%syf) .and. b%sideyl%proc==my_index) &
  call em3d_applybc_splite(b%sideyl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sideyr%syf) .and. b%sideyr%proc==my_index) &
  call em3d_applybc_splite(b%sideyr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sidezl%syf) .and. b%sidezl%proc==my_index) &
  call em3d_applybc_splite(b%sidezl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sidezr%syf) .and. b%sidezr%proc==my_index) &
  call em3d_applybc_splite(b%sidezr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)

  if (associated(b%edgexlyl%syf) .and. b%edgexlyl%proc==my_index) &
  call em3d_applybc_splite(b%edgexlyl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexryl%syf) .and. b%edgexryl%proc==my_index) &
  call em3d_applybc_splite(b%edgexryl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexlyr%syf) .and. b%edgexlyr%proc==my_index) &
  call em3d_applybc_splite(b%edgexlyr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexryr%syf) .and. b%edgexryr%proc==my_index) &
  call em3d_applybc_splite(b%edgexryr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexlzl%syf) .and. b%edgexlzl%proc==my_index) &
  call em3d_applybc_splite(b%edgexlzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexrzl%syf) .and. b%edgexrzl%proc==my_index)  &
  call em3d_applybc_splite(b%edgexrzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexlzr%syf) .and. b%edgexlzr%proc==my_index)  &
  call em3d_applybc_splite(b%edgexlzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexrzr%syf) .and. b%edgexrzr%proc==my_index) &
  call em3d_applybc_splite(b%edgexrzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeylzl%syf) .and. b%edgeylzl%proc==my_index) &
  call em3d_applybc_splite(b%edgeylzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeyrzl%syf) .and. b%edgeyrzl%proc==my_index) &
  call em3d_applybc_splite(b%edgeyrzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeylzr%syf) .and. b%edgeylzr%proc==my_index) &
  call em3d_applybc_splite(b%edgeylzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeyrzr%syf) .and. b%edgeyrzr%proc==my_index) &
  call em3d_applybc_splite(b%edgeyrzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)

  return
end subroutine push_em3d_blockbnde

subroutine push_em3d_blockbndb(b,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
integer(ISZ) :: which

  if(b%xlbnd==openbc) call push_em3d_splitb(b%sidexl%syf,dt,which)
  if(b%xrbnd==openbc) call push_em3d_splitb(b%sidexr%syf,dt,which)
  if(b%ylbnd==openbc) call push_em3d_splitb(b%sideyl%syf,dt,which)
  if(b%yrbnd==openbc) call push_em3d_splitb(b%sideyr%syf,dt,which)
  if(b%zlbnd==openbc) call push_em3d_splitb(b%sidezl%syf,dt,which)
  if(b%zrbnd==openbc) call push_em3d_splitb(b%sidezr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splitb(b%edgexlyl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splitb(b%edgexryl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splitb(b%edgexlyr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splitb(b%edgexryr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%edgexlzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%edgexrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%edgexlzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%edgexrzr%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%edgeylzl%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%edgeyrzl%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%edgeylzr%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%edgeyrzr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%cornerxlylzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%cornerxrylzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%cornerxlyrzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitb(b%cornerxryrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%cornerxlylzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%cornerxrylzr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%cornerxlyrzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitb(b%cornerxryrzr%syf,dt,which)

  if (associated(b%sidexl%syf) .and. b%sidexl%proc==my_index) &
  call em3d_applybc_splitb(b%sidexl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sidexr%syf) .and. b%sidexr%proc==my_index) &
  call em3d_applybc_splitb(b%sidexr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sideyl%syf) .and. b%sideyl%proc==my_index) &
  call em3d_applybc_splitb(b%sideyl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sideyr%syf) .and. b%sideyr%proc==my_index) &
  call em3d_applybc_splitb(b%sideyr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sidezl%syf) .and. b%sidezl%proc==my_index) &
  call em3d_applybc_splitb(b%sidezl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%sidezr%syf) .and. b%sidezr%proc==my_index) &
  call em3d_applybc_splitb(b%sidezr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)

  if (associated(b%edgexlyl%syf) .and. b%edgexlyl%proc==my_index) &
  call em3d_applybc_splitb(b%edgexlyl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexryl%syf) .and. b%edgexryl%proc==my_index) &
  call em3d_applybc_splitb(b%edgexryl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexlyr%syf) .and. b%edgexlyr%proc==my_index) &
  call em3d_applybc_splitb(b%edgexlyr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexryr%syf) .and. b%edgexryr%proc==my_index) &
  call em3d_applybc_splitb(b%edgexryr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexlzl%syf) .and. b%edgexlzl%proc==my_index) &
  call em3d_applybc_splitb(b%edgexlzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexrzl%syf) .and. b%edgexrzl%proc==my_index)  &
  call em3d_applybc_splitb(b%edgexrzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexlzr%syf) .and. b%edgexlzr%proc==my_index)  &
  call em3d_applybc_splitb(b%edgexlzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgexrzr%syf) .and. b%edgexrzr%proc==my_index) &
  call em3d_applybc_splitb(b%edgexrzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeylzl%syf) .and. b%edgeylzl%proc==my_index) &
  call em3d_applybc_splitb(b%edgeylzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeyrzl%syf) .and. b%edgeyrzl%proc==my_index) &
  call em3d_applybc_splitb(b%edgeyrzl%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeylzr%syf) .and. b%edgeylzr%proc==my_index) &
  call em3d_applybc_splitb(b%edgeylzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)
  if (associated(b%edgeyrzr%syf) .and. b%edgeyrzr%proc==my_index) &
  call em3d_applybc_splitb(b%edgeyrzr%syf,b%xlbnd,b%xrbnd,b%ylbnd,b%yrbnd,b%zlbnd,b%zrbnd)

  return
end subroutine push_em3d_blockbndb

subroutine push_em3d_blockbndef(b,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
integer(ISZ) :: which

  if(b%xlbnd==openbc) call push_em3d_splitef(b%sidexl%syf,dt,which)
  if(b%xrbnd==openbc) call push_em3d_splitef(b%sidexr%syf,dt,which)
  if(b%ylbnd==openbc) call push_em3d_splitef(b%sideyl%syf,dt,which)
  if(b%yrbnd==openbc) call push_em3d_splitef(b%sideyr%syf,dt,which)
  if(b%zlbnd==openbc) call push_em3d_splitef(b%sidezl%syf,dt,which)
  if(b%zrbnd==openbc) call push_em3d_splitef(b%sidezr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splitef(b%edgexlyl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splitef(b%edgexryl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splitef(b%edgexlyr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splitef(b%edgexryr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%edgexlzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%edgexrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%edgexlzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%edgexrzr%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%edgeylzl%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%edgeyrzl%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%edgeylzr%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%edgeyrzr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%cornerxlylzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%cornerxrylzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%cornerxlyrzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitef(b%cornerxryrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%cornerxlylzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%cornerxrylzr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%cornerxlyrzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitef(b%cornerxryrzr%syf,dt,which)

  return
end subroutine push_em3d_blockbndef

subroutine push_em3d_blockbndf(b,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
integer(ISZ) :: which

  if(b%xlbnd==openbc) call push_em3d_splitf(b%sidexl%syf,dt,which)
  if(b%xrbnd==openbc) call push_em3d_splitf(b%sidexr%syf,dt,which)
  if(b%ylbnd==openbc) call push_em3d_splitf(b%sideyl%syf,dt,which)
  if(b%yrbnd==openbc) call push_em3d_splitf(b%sideyr%syf,dt,which)
  if(b%zlbnd==openbc) call push_em3d_splitf(b%sidezl%syf,dt,which)
  if(b%zrbnd==openbc) call push_em3d_splitf(b%sidezr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splitf(b%edgexlyl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call push_em3d_splitf(b%edgexryl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splitf(b%edgexlyr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call push_em3d_splitf(b%edgexryr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%edgexlzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%edgexrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%edgexlzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%edgexrzr%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%edgeylzl%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%edgeyrzl%syf,dt,which)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%edgeylzr%syf,dt,which)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%edgeyrzr%syf,dt,which)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%cornerxlylzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%cornerxrylzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%cornerxlyrzl%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call push_em3d_splitf(b%cornerxryrzl%syf,dt,which)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%cornerxlylzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%cornerxrylzr%syf,dt,which)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%cornerxlyrzr%syf,dt,which)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call push_em3d_splitf(b%cornerxryrzr%syf,dt,which)

  return
end subroutine push_em3d_blockbndf

subroutine shift_em3dblock_ncells_z(b,n)
use mod_emfield3d
implicit none
TYPE(EM3D_BLOCKtype) :: b
integer(ISZ):: n

  ! --- shift core 
  call shift_em3df_ncells_z(b%core%yf,b%zlbnd,b%zrbnd,n)
  ! --- shift sides
  if(b%xlbnd==openbc) call shift_em3dsplitf_ncells_z(b%sidexl%syf,b%zlbnd,b%zrbnd,n)
  if(b%xrbnd==openbc) call shift_em3dsplitf_ncells_z(b%sidexr%syf,b%zlbnd,b%zrbnd,n)
  if(b%ylbnd==openbc) call shift_em3dsplitf_ncells_z(b%sideyl%syf,b%zlbnd,b%zrbnd,n)
  if(b%yrbnd==openbc) call shift_em3dsplitf_ncells_z(b%sideyr%syf,b%zlbnd,b%zrbnd,n)
  if(b%zlbnd==openbc) call shift_em3dsplitf_ncells_z(b%sidezl%syf,openbc,openbc,n)
  if(b%zrbnd==openbc) call shift_em3dsplitf_ncells_z(b%sidezr%syf,openbc,openbc,n)
  ! --- shift edges
  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexlyl%syf,b%zlbnd,b%zrbnd,n)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexryl%syf,b%zlbnd,b%zrbnd,n)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexlyr%syf,b%zlbnd,b%zrbnd,n)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexryr%syf,b%zlbnd,b%zrbnd,n)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexlzl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexrzl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexlzr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgexrzr%syf,openbc,openbc,n)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgeylzl%syf,openbc,openbc,n)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgeyrzl%syf,openbc,openbc,n)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgeylzr%syf,openbc,openbc,n)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_z(b%edgeyrzr%syf,openbc,openbc,n)
  ! --- shift corners
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxlylzl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxrylzl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxlyrzl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxryrzl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxlylzr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxrylzr%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxlyrzr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_z(b%cornerxryrzr%syf,openbc,openbc,n)
  
  return
end subroutine shift_em3dblock_ncells_z

subroutine shift_em3df_ncells_z(f,zl,zr,n)
use mod_emfield3d
implicit none
TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ):: n,i,it,zl,zr
  call shift_3darray_ncells_z(f%Ex,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Ey,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Ez,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Bx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%By,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Bz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  do it=1,f%ntimes
    do i=1,3
      call shift_3darray_ncells_z(f%Jarray(-f%nxguard,-f%nyguard,-f%nzguard,i,it), &
                                    f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    end do
  end do
  if (f%nxf>0) then
    call shift_3darray_ncells_z(f%F,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    call shift_3darray_ncells_z(f%Rhoold,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Rho,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    do it=1,f%ntimes
!      call shift_3darray_ncells_z(f%Rhoarray(-f%nxguard,-f%nyguard,-f%nzguard,it), &
!                                    f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    end do
  end if
  if (f%nxpo>0) then
    call shift_3darray_ncells_z(f%Ax,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Ay,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Az,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Phi,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  end if
  
  return
end subroutine shift_em3df_ncells_z

subroutine shift_em3dsplitf_ncells_z(f,zl,zr,n)
use mod_emfield3d
implicit none
TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ):: n,zl,zr

  call shift_3darray_ncells_z(f%Exx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Exy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Exz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Eyx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Eyy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Eyz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Ezx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Ezy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Ezz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Bxy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Bxz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Byx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Byz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Bzx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Bzy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Fx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Fy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  call shift_3darray_ncells_z(f%Fz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  if (f%nxpo>0) then
    call shift_3darray_ncells_z(f%Ax,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Ay,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Az,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Phi,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  end if

  return
end subroutine shift_em3dsplitf_ncells_z

subroutine shift_3darray_ncells_z(f,nx,ny,nz,nxguard,nyguard,nzguard,zl,zr,n)
#ifdef MPIPARALLEL
use mpirz
#endif
implicit none
integer(ISZ) :: nx,ny,nz,nxguard,nyguard,nzguard,n,zl,zr
integer(ISZ), parameter:: otherproc=10, ibuf = 25
real(kind=8) :: f(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) 

f(:,:,-nzguard:nz+nzguard-n) = f(:,:,-nzguard+n:nz+nzguard)
if (zr/=otherproc) f(:,:,nz-n+1:nz+nzguard) = 0.

#ifdef MPIPARALLEL
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,:,-n+1:nzguard)),ibuf)
     call mympi_pack(f(:,:,-n+1:nzguard),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if    
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,:,nz-n+1:nz+nzguard)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,:,nz:nz+nzguard) = reshape(mpi_unpack_real_array( size(f(:,:,nz-n+1:nz+nzguard)),ibuf), &
                                                         shape(f(:,:,nz-n+1:nz+nzguard)))
  end if
#endif

  return
end subroutine shift_3darray_ncells_z

subroutine shift_3darray_ncells_zold(f,nx,ny,nz,nxguard,nyguard,nzguard,zl,zr,n)
#ifdef MPIPARALLEL
use mpirz
#endif
implicit none
integer(ISZ) :: nx,ny,nz,nxguard,nyguard,nzguard,n,zl,zr
integer(ISZ), parameter:: otherproc=10, ibuf = 25
real(kind=8) :: f(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) 

f(:,:,-nzguard:nz+nzguard-n) = f(:,:,-nzguard+n:nz+nzguard)
if (zr/=otherproc) f(:,:,nz+nzguard-n+1:nz+nzguard) = 0.

#ifdef MPIPARALLEL
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,:,0:nzguard)),ibuf)
     call mympi_pack(f(:,:,0:nzguard),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if    
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,:,nz:nz+nzguard)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,:,nz:nz+nzguard) = reshape(mpi_unpack_real_array( size(f(:,:,nz:nz+nzguard)),ibuf), &
                                                         shape(f(:,:,nz:nz+nzguard)))
  end if
#endif

  return
end subroutine shift_3darray_ncells_zold

subroutine em3d_applybc_e(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%ex(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = -f%ex(f%ixmin:f%ixmin+f%nxguard-1,:,:)
     f%ey(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = f%ey(f%ixmin+1:f%ixmin+f%nxguard,:,:)
     f%ez(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = f%ez(f%ixmin+1:f%ixmin+f%nxguard,:,:)
  end if
  if (xrbnd==neumann) then
     f%ex(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:) = -f%ex(f%ixmax-f%nxguard-1:f%ixmax,:,:)
     f%ey(f%ixmax+1:f%ixmax+f%nxguard:-1,:,:) = f%ey(f%ixmax-f%nxguard:f%ixmax-1,:,:)
     f%ez(f%ixmax+1:f%ixmax+f%nxguard:-1,:,:) = f%ez(f%ixmax-f%nxguard:f%ixmax-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%ey(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = -f%ey(:,f%iymin:f%iymin+f%nyguard-1,:)
     f%ex(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = f%ex(:,f%iymin+1:f%iymin+f%nyguard,:)
     f%ez(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = f%ez(:,f%iymin+1:f%iymin+f%nyguard,:)
  end if
  if (yrbnd==neumann) then
     f%ey(:,f%iymax:f%iymax+f%nyguard-1:-1,:) = -f%ey(:,f%iymax-f%nyguard-1:f%iymax,:)
     f%ex(:,f%iymax+1:f%iymax+f%nyguard:-1,:) = f%ex(:,f%iymax-f%nyguard:f%iymax-1,:)
     f%ez(:,f%iymax+1:f%iymax+f%nyguard:-1,:) = f%ez(:,f%iymax-f%nyguard:f%iymax-1,:)
  end if

  return
end subroutine em3d_applybc_e

subroutine em3d_applybc_b(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%bx(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = f%bx(f%ixmin+1:f%ixmin+f%nxguard,:,:)
     f%by(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = -f%by(f%ixmin:f%ixmin+f%nxguard-1,:,:)
     f%bz(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = -f%bz(f%ixmin:f%ixmin+f%nxguard-1,:,:)
  end if
  if (xrbnd==neumann) then
     f%bx(f%ixmax+1:f%ixmax+f%nxguard:-1,:,:) = f%bx(f%ixmax-f%nxguard:f%ixmax-1,:,:)
     f%by(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:) = -f%by(f%ixmax-f%nxguard-1:f%ixmax,:,:)
     f%bz(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:) = -f%bz(f%ixmax-f%nxguard-1:f%ixmax,:,:)
  end if

  if (ylbnd==neumann) then
     f%by(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = f%by(:,f%iymin+1:f%iymin+f%nyguard,:)
     f%bx(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = -f%bx(:,f%iymin:f%iymin+f%nyguard-1,:)
     f%bz(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = -f%bz(:,f%iymin:f%iymin+f%nyguard-1,:)
  end if
  if (yrbnd==neumann) then
     f%by(:,f%iymax+1:f%iymax+f%nyguard:-1,:) = f%by(:,f%iymax-f%nyguard:f%iymax-1,:)
     f%bx(:,f%iymax:f%iymax+f%nyguard-1:-1,:) = -f%bx(:,f%iymax-f%nyguard-1:f%iymax,:)
     f%bz(:,f%iymax:f%iymax+f%nyguard-1:-1,:) = -f%bz(:,f%iymax-f%nyguard-1:f%iymax,:)
  end if

  return
end subroutine em3d_applybc_b

subroutine em3d_applybc_j(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) = f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) + f%j(f%ixmin:f%ixmin-f%nxguard:-1,:,:,2:3)
     f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) = f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) - f%j(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:,1)
  end if
  if (xrbnd==neumann) then
     f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) = f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) + f%j(f%ixmax:f%ixmax+f%nxguard:-1,:,:,2:3)
     f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) = f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) - f%j(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:,1)
  end if

  if (ylbnd==neumann) then
     f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) = f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) &
                                              + f%j(:,f%iymin:f%iymin-f%nyguard:-1,:,1:3:2)
     f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) = f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) - f%j(:,f%iymin-1:f%iymin-f%nyguard:-1,:,2)
  end if
  if (yrbnd==neumann) then
     f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) = f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) &
                                              + f%j(:,f%iymax:f%iymax+f%nyguard:-1,:,1:3:2)
     f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) = f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) - f%j(:,f%iymax:f%iymax+f%nyguard-1:-1,:,2)
  end if

  return
end subroutine em3d_applybc_j

subroutine em3d_applybc_rho(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) = f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) + f%rho(f%ixmin:f%ixmin-f%nxguard:-1,:,:)
  end if
  if (xrbnd==neumann) then
     f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) = f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) + f%rho(f%ixmax:f%ixmax+f%nxguard:-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%rho(:,f%iymin:f%iymin+f%nyguard,:) = f%rho(:,f%iymin:f%iymin+f%nyguard,:) + f%rho(:,f%iymin:f%iymin-f%nyguard:-1,:)
  end if
  if (yrbnd==neumann) then
     f%rho(:,f%iymax-f%nyguard:f%iymax,:) = f%rho(:,f%iymax-f%nyguard:f%iymax,:) + f%rho(:,f%iymax:f%iymax+f%nyguard:-1,:)
  end if

  return
end subroutine em3d_applybc_rho

subroutine em3d_applybc_splite(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%exx(f%ixmin-1,:,:) = -f%exx(f%ixmin  ,:,:)
     f%exy(f%ixmin-1,:,:) = -f%exy(f%ixmin  ,:,:)
     f%exz(f%ixmin-1,:,:) = -f%exz(f%ixmin  ,:,:)
  end if
  if (xrbnd==neumann) then
     f%exx(f%ixmax  ,:,:) = -f%exx(f%ixmax-1,:,:)
     f%exy(f%ixmax  ,:,:) = -f%exy(f%ixmax-1,:,:)
     f%exz(f%ixmax  ,:,:) = -f%exz(f%ixmax-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%eyx(:,f%iymin-1,:) = -f%eyx(:,f%iymin  ,:)
     f%eyy(:,f%iymin-1,:) = -f%eyy(:,f%iymin  ,:)
     f%eyz(:,f%iymin-1,:) = -f%eyz(:,f%iymin  ,:)
  end if
  if (yrbnd==neumann) then
     f%eyx(:,f%iymax  ,:) = -f%eyx(:,f%iymax-1,:)
     f%eyy(:,f%iymax  ,:) = -f%eyy(:,f%iymax-1,:)
     f%eyz(:,f%iymax  ,:) = -f%eyz(:,f%iymax-1,:)
  end if

  return
end subroutine em3d_applybc_splite

subroutine em3d_applybc_splitb(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%byx(f%ixmin-1,:,:) = -f%byx(f%ixmin  ,:,:)
     f%byz(f%ixmin-1,:,:) = -f%byz(f%ixmin  ,:,:)
     f%bzx(f%ixmin-1,:,:) = -f%bzx(f%ixmin  ,:,:)
     f%bzy(f%ixmin-1,:,:) = -f%bzy(f%ixmin  ,:,:)
  end if
  if (xrbnd==neumann) then
     f%byx(f%ixmax  ,:,:) = -f%byx(f%ixmax-1,:,:)
     f%byz(f%ixmax  ,:,:) = -f%byz(f%ixmax-1,:,:)
     f%bzx(f%ixmax  ,:,:) = -f%bzx(f%ixmax-1,:,:)
     f%bzy(f%ixmax  ,:,:) = -f%bzy(f%ixmax-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%bxy(:,f%iymin-1,:) = -f%bxy(:,f%iymin  ,:)
     f%bxz(:,f%iymin-1,:) = -f%bxz(:,f%iymin  ,:)
     f%bzx(:,f%iymin-1,:) = -f%bzx(:,f%iymin  ,:)
     f%bzy(:,f%iymin-1,:) = -f%bzy(:,f%iymin  ,:)
  end if
  if (yrbnd==neumann) then
     f%bxy(:,f%iymax  ,:) = -f%bxy(:,f%iymax-1,:)
     f%bxz(:,f%iymax  ,:) = -f%bxz(:,f%iymax-1,:)
     f%bzx(:,f%iymax  ,:) = -f%bzx(:,f%iymax-1,:)
     f%bzy(:,f%iymax  ,:) = -f%bzy(:,f%iymax-1,:)
  end if

  return
end subroutine em3d_applybc_splitb

subroutine em3d_exchange_bnde_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
integer(ISZ)   ::ix,ibuf
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          ! --- case lower yee, upper yee
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 1
            call mpi_packbuffer_init((3*yfu%nxguard-2)*size(yfu%ez(0,:,:)),ibuf)
            do ix = yfu%ixmin,yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%ex(ix,:,:),ibuf)
            end do
            if (yfu%nxguard>1) then
              do ix = yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
                call mympi_pack(yfu%ey(ix,:,:),ibuf)
              end do
              do ix = yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
                call mympi_pack(yfu%ez(ix,:,:),ibuf)
              end do
            end if
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 2
            call mpi_packbuffer_init((3*yfl%nxguard-2)*size(yfl%ez(0,:,:)),ibuf)
            do ix =yfl%ixmax-yfl%nxguard,yfl%ixmax-1
              call mympi_pack(yfl%ex(ix,:,:),ibuf)
            end do
            if (yfl%nxguard>1) then
              do ix = yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
                call mympi_pack(yfl%ey(ix,:,:),ibuf)
              end do
              do ix = yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
                call mympi_pack(yfl%ez(ix,:,:),ibuf)
              end do
            end if
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
            yfl%ex(yfl%ixmax  :yfl%ixmaxg-1,:,:) = yfu%ex(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)
            yfl%ey(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%ey(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
            yfl%ez(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%ez(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)

            yfu%ex(yfu%ixming  :yfu%ixmin-1,:,:) = yfl%ex(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)
            yfu%ey(yfu%ixming+1:yfu%ixmin-1,:,:) = yfl%ey(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
            yfu%ez(yfu%ixming+1:yfu%ixmin-1,:,:) = yfl%ez(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          ! --- case lower yee, upper split yee
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%ex(yfl%ixmax  :yfl%ixmaxg-1,:,:) = syfu%exx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:) &
                                               + syfu%exy(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:) &
                                               + syfu%exz(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          yfl%ey(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = syfu%eyx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                               + syfu%eyy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                               + syfu%eyz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          yfl%ez(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = syfu%ezx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                               + syfu%ezy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                               + syfu%ezz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfu%exx(syfu%ixming  :syfu%ixmin-1,:,:) = 0.
          syfu%exy(syfu%ixming  :syfu%ixmin-1,:,:) = 0.
          syfu%exz(syfu%ixming  :syfu%ixmin-1,:,:) = yfl%ex(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)
          syfu%eyx(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%eyy(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%eyz(syfu%ixming+1:syfu%ixmin-1,:,:) = yfl%ey(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
          syfu%ezx(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%ezy(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%ezz(syfu%ixming+1:syfu%ixmin-1,:,:) = yfl%ez(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
          ! --- case lower split yee, upper yee
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%ex(yfu%ixming  :yfu%ixmin-1,:,:)      = syfl%exx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:) &
                                                    + syfl%exy(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:) &
                                                    + syfl%exz(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          yfu%ey(yfu%ixming+1:yfu%ixmin-1,:,:)      = syfl%eyx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                                    + syfl%eyy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                                    + syfl%eyz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          yfu%ez(yfu%ixming+1:yfu%ixmin-1,:,:)      = syfl%ezx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                                    + syfl%ezy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                                    + syfl%ezz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfl%exx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = 0.
          syfl%exy(syfl%ixmax  :syfl%ixmaxg-1,:,:) = 0.
          syfl%exz(syfl%ixmax  :syfl%ixmaxg-1,:,:) = yfu%ex(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)
          syfl%eyx(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%eyy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%eyz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = yfu%ey(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
          syfl%ezx(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%ezy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%ezz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = yfu%ez(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
        case(splityeefield)
          ! --- case lower split yee, upper split yee
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 3
            call mpi_packbuffer_init( 6*int(size(syfu%ezx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:))) &
                                    + 3*int(size(syfu%exx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:))) ,ibuf)

            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%exx(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%exy(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%exz(ix,:,:),ibuf)
            end do

            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%eyx(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%eyy(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%eyz(ix,:,:),ibuf)
            end do

            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%ezx(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%ezy(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%ezz(ix,:,:),ibuf)
            end do

            call mpi_isend_pack(fl%proc,3,ibuf)

            mpireqpnt=mpireqpnt+1
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 4
            call mpi_packbuffer_init(6*int(size(syfl%ezx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:))) &
                                    +3*int(size(syfl%exx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:))),ibuf)

            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%exx(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%exy(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%exz(ix,:,:),ibuf)
            end do

            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%eyx(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%eyy(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%eyz(ix,:,:),ibuf)
            end do

            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%ezx(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%ezy(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%ezz(ix,:,:),ibuf)
            end do

            call mpi_isend_pack(fu%proc,4,ibuf)

            mpireqpnt=mpireqpnt+1
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%exx(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%exx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfu%exy(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%exy(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfu%exz(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%exz(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfu%eyx(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%eyx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%eyy(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%eyy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%eyz(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%eyz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%ezx(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%ezx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%ezy(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%ezy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%ezz(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%ezz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)

          syfl%exx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%exx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%exy(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%exy(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%exz(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%exz(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%eyx(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%eyx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%eyy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%eyy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%eyz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%eyz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%ezx(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%ezx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%ezy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%ezy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%ezz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%ezz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bnde_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bnde_xrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 5
            call mpi_packbuffer_init((3*yfu%nxguard-2)*size(yfu%Ez(0,:,:)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do ix = yfu%ixming,yfu%ixmin-1
              yfu%ex(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%Ex(0,:,:)),ibuf),shape(yfu%Ex(0,:,:)))
            end do
            if (yfu%nxguard>1) then
              do ix = yfu%ixming+1,yfu%ixmin-1
                yfu%ey(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%Ez(0,:,:)),ibuf),shape(yfu%Ez(0,:,:)))
              end do
              do ix = yfu%ixming+1,yfu%ixmin-1
                yfu%ez(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%Ez(0,:,:)),ibuf),shape(yfu%Ez(0,:,:)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 6
            call mpi_packbuffer_init((3*yfl%nxguard-2)*size(yfl%ez(yfl%ixmin,:,:)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do ix = yfl%ixmax,yfl%ixmaxg-1
              yfl%ex(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(0,:,:)),ibuf),shape(yfl%Ez(0,:,:)))
            end do
            if (yfl%nxguard>1) then
              do ix = yfl%ixmax+1,yfl%ixmaxg-1
                yfl%ey(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(0,:,:)),ibuf),shape(yfl%Ez(0,:,:)))
              end do
              do ix = yfl%ixmax+1,yfl%ixmaxg-1
                yfl%ez(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(0,:,:)),ibuf),shape(yfl%Ez(0,:,:)))
              end do
            end if
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 7
            call mpi_packbuffer_init(6*size(syfu%ezx(syfu%ixming+1:syfu%ixmin-1,:,:)) &
                                    +3*size(syfu%ezx(syfu%ixming  :syfu%ixmin-1,:,:)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)

            do ix = syfu%ixming,syfu%ixmin-1
              syfu%exx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%exx(ix,:,:)),ibuf),shape(syfu%exx(ix,:,:)))
            end do
            do ix = syfu%ixming,syfu%ixmin-1
              syfu%exy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%exy(ix,:,:)),ibuf),shape(syfu%exy(ix,:,:)))
            end do
            do ix = syfu%ixming,syfu%ixmin-1
              syfu%exz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%exz(ix,:,:)),ibuf),shape(syfu%exz(ix,:,:)))
            end do

            do ix = syfu%ixming+1,syfu%ixmin-1
              syfu%eyx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%eyx(ix,:,:)),ibuf),shape(syfu%eyx(ix,:,:)))
            end do
            do ix = syfu%ixming+1,syfu%ixmin-1
              syfu%eyy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%eyy(ix,:,:)),ibuf),shape(syfu%eyy(ix,:,:)))
            end do
            do ix = syfu%ixming+1,syfu%ixmin-1
              syfu%eyz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%eyz(ix,:,:)),ibuf),shape(syfu%eyz(ix,:,:)))
            end do

            do ix = syfu%ixming+1,syfu%ixmin-1
              syfu%ezx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%ezx(ix,:,:)),ibuf),shape(syfu%ezx(ix,:,:)))
            end do
            do ix = syfu%ixming+1,syfu%ixmin-1
              syfu%ezy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%ezy(ix,:,:)),ibuf),shape(syfu%ezy(ix,:,:)))
            end do
            do ix = syfu%ixming+1,syfu%ixmin-1
              syfu%ezz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%ezz(ix,:,:)),ibuf),shape(syfu%ezz(ix,:,:)))
            end do

          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 8
            call mpi_packbuffer_init(6*size(syfl%ezx(syfl%ixmax+1:syfl%ixmaxg-1,:,:)) &
                                    +3*size(syfl%ezx(syfl%ixmax  :syfl%ixmaxg-1,:,:)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%exx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%exx(ix,:,:)),ibuf),shape(syfl%exx(ix,:,:)))
            end do
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%exy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%exy(ix,:,:)),ibuf),shape(syfl%exy(ix,:,:)))
            end do
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%exz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%exz(ix,:,:)),ibuf),shape(syfl%exz(ix,:,:)))
            end do
            
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%eyx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%eyx(ix,:,:)),ibuf),shape(syfl%eyx(ix,:,:)))
            end do
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%eyy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%eyy(ix,:,:)),ibuf),shape(syfl%eyy(ix,:,:)))
            end do
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%eyz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%eyz(ix,:,:)),ibuf),shape(syfl%eyz(ix,:,:)))
            end do
            
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%ezx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%ezx(ix,:,:)),ibuf),shape(syfl%ezx(ix,:,:)))
            end do
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%ezy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%ezy(ix,:,:)),ibuf),shape(syfl%ezy(ix,:,:)))
            end do
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%ezz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%ezz(ix,:,:)),ibuf),shape(syfl%ezz(ix,:,:)))
            end do

          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bnde_xrecv
#endif

subroutine em3d_exchange_bndb_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
integer(ISZ) :: ibuf,ix
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 9
            call mpi_packbuffer_init((3*yfu%nxguard-1)*size(yfu%by(0,:,:)),ibuf)
            do ix=yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%bx(ix,:,:),ibuf)
            end do
            do ix=yfu%ixmin, yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%by(ix,:,:),ibuf)
            end do
            do ix=yfu%ixmin,yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%bz(ix,:,:),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 10
            call mpi_packbuffer_init((3*yfl%nxguard-1)*size(yfl%by(0,:,:)),ibuf)
            do ix=yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
              call mympi_pack(yfl%bx(ix,:,:),ibuf)
            end do
            do ix=yfl%ixmax-yfl%nxguard,yfl%ixmax-1
              call mympi_pack(yfl%by(ix,:,:),ibuf)
            end do
            do ix=yfl%ixmax-yfl%nxguard,yfl%ixmax-1
              call mympi_pack(yfl%bz(ix,:,:),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
          yfl%bx(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%bx(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
          yfl%by(yfl%ixmax  :yfl%ixmaxg-1,:,:) = yfu%by(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)
          yfl%bz(yfl%ixmax  :yfl%ixmaxg-1,:,:) = yfu%bz(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)

          yfu%bx(yfu%ixming+1:yfu%ixmin-1,:,:) = yfl%bx(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
          yfu%by(yfu%ixming  :yfu%ixmin-1,:,:) = yfl%by(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)
          yfu%bz(yfu%ixming  :yfu%ixmin-1,:,:) = yfl%bz(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)
          
#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=fu%proc) return
          yfl%bx(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = (syfu%bxy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                               +  syfu%bxz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:))/syfu%clight
          yfl%by(yfl%ixmax  :yfl%ixmaxg-1,:,:) = (syfu%byx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:) &
                                               +  syfu%byz(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:))/syfu%clight
          yfl%bz(yfl%ixmax  :yfl%ixmaxg-1,:,:) = (syfu%bzx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:) &
                                               +  syfu%bzy(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:))/syfu%clight
          syfu%bxy(syfu%ixming+1:syfu%ixmin-1,:,:) = yfl%bx(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)*syfu%clight
          syfu%bxz(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%byx(syfu%ixming  :syfu%ixmin-1,:,:) = yfl%by(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)*syfu%clight
          syfu%byz(syfu%ixming  :syfu%ixmin-1,:,:) = 0.
          syfu%bzx(syfu%ixming  :syfu%ixmin-1,:,:) = yfl%bz(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)*syfu%clight
          syfu%bzy(syfu%ixming  :syfu%ixmin-1,:,:) = 0.
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=fu%proc) return
          yfu%bx(yfu%ixming+1:yfu%ixmin-1,:,:) = (syfl%bxy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                               +  syfl%bxz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:))/syfl%clight
          yfu%by(yfu%ixming  :yfu%ixmin-1,:,:) = (syfl%byx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:) &
                                               +  syfl%byz(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:))/syfl%clight
          yfu%bz(yfu%ixming  :yfu%ixmin-1,:,:) = (syfl%bzx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:) &
                                               +  syfl%bzy(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:))/syfl%clight
          syfl%bxy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = yfu%bx(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)*syfl%clight
          syfl%bxz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%byx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = yfu%by(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)*syfl%clight
          syfl%byz(syfl%ixmax  :syfl%ixmaxg-1,:,:) = 0.
          syfl%bzx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = yfu%bz(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)*syfl%clight
          syfl%bzy(syfl%ixmax  :syfl%ixmaxg-1,:,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 11
            call mpi_packbuffer_init(4*size(syfu%byx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)) &
                                    +2*size(syfu%byx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)),ibuf)

            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%bxy(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%bxz(ix,:,:),ibuf)
            end do
            
            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%byx(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%byz(ix,:,:),ibuf)
            end do
            
            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%bzx(ix,:,:),ibuf)
            end do
            do ix=syfu%ixmin,syfu%ixmin+syfu%nxguard-1
              call mympi_pack(syfu%bzy(ix,:,:),ibuf)
            end do

            call mpi_isend_pack(fl%proc,3,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 12
            call mpi_packbuffer_init(4*size(syfl%byx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)) &
                                    +2*size(syfl%byx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)),ibuf)

            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%bxy(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
              call mympi_pack(syfl%bxz(ix,:,:),ibuf)
            end do

            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%byx(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%byz(ix,:,:),ibuf)
            end do

            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%bzx(ix,:,:),ibuf)
            end do
            do ix=syfl%ixmax-syfl%nxguard,syfl%ixmax-1
              call mympi_pack(syfl%bzy(ix,:,:),ibuf)
            end do

            call mpi_isend_pack(fu%proc,4,ibuf)
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%bxy(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%bxy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%bxz(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%bxz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%byx(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%byx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfu%byz(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%byz(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfu%bzx(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%bzx(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfu%bzy(syfu%ixming  :syfu%ixmin-1,:,:) = syfl%bzy(syfl%ixmax-syfl%nxguard  :syfl%ixmax-1,:,:)
          syfl%bxy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%bxy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%bxz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%bxz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%byx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%byx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%byz(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%byz(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%bzx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%bzx(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%bzy(syfl%ixmax  :syfl%ixmaxg-1,:,:) = syfu%bzy(syfu%ixmin  :syfu%ixmin+syfu%nxguard-1,:,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndb_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndb_xrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ibuf,ix

          if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 13
            call mpi_packbuffer_init((3*yfu%nxguard-1)*size(yfu%bx(0,:,:)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do ix=yfu%ixming+1,yfu%ixmin-1
              yfu%bx(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%bx(ix,:,:)),ibuf),shape(yfu%bx(ix,:,:)))
            end do
            do ix=yfu%ixming,yfu%ixmin-1
              yfu%by(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%by(ix,:,:)),ibuf),shape(yfu%by(ix,:,:)))
            end do
            do ix=yfu%ixming,yfu%ixmin-1
              yfu%bz(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%bz(ix,:,:)),ibuf),shape(yfu%bz(ix,:,:)))
            end do
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 14
            call mpi_packbuffer_init((3*yfl%nxguard-1)*size(yfl%bx(0,:,:)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do ix=yfl%ixmax+1,yfl%ixmaxg-1
              yfl%bx(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%bx(ix,:,:)),ibuf),shape(yfl%bx(ix,:,:)))
            end do
            do ix=yfl%ixmax,yfl%ixmaxg-1
              yfl%by(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%by(ix,:,:)),ibuf),shape(yfl%by(ix,:,:)))
            end do
            do ix=yfl%ixmax,yfl%ixmaxg-1
              yfl%bz(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%bz(ix,:,:)),ibuf),shape(yfl%bz(ix,:,:)))
            end do
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ibuf = 15
            call mpi_packbuffer_init(4*size(syfu%bxy(syfu%ixming  :syfu%ixmin-1,:,:)) &
                                    +2*size(syfu%bxy(syfu%ixming+1:syfu%ixmin-1,:,:)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)
            ! --- recv data from down in z
            do ix=syfu%ixming+1,syfu%ixmin-1
              syfu%bxy(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfu%bxy(ix,:,:)),ibuf),shape(syfu%bxy(ix,:,:)))
            end do
            do ix=syfu%ixming+1,syfu%ixmin-1
              syfu%bxz(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfu%bxz(ix,:,:)),ibuf),shape(syfu%bxz(ix,:,:)))
            end do
            do ix=syfu%ixming,syfu%ixmin-1
              syfu%byx(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfu%byx(ix,:,:)),ibuf),shape(syfu%byx(ix,:,:)))
            end do
            do ix=syfu%ixming,syfu%ixmin-1
              syfu%byz(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfu%byz(ix,:,:)),ibuf),shape(syfu%byz(ix,:,:)))
            end do
            do ix=syfu%ixming,syfu%ixmin-1
              syfu%bzx(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfu%bzx(ix,:,:)),ibuf),shape(syfu%bzx(ix,:,:)))
            end do
            do ix=syfu%ixming,syfu%ixmin-1
              syfu%bzy(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfu%bzy(ix,:,:)),ibuf),shape(syfu%bzy(ix,:,:)))
            end do
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 16
            call mpi_packbuffer_init(4*size(syfl%bxy(syfl%ixmax  :syfl%ixmaxg-1,:,:)) &
                                    +2*size(syfl%bxy(syfl%ixmax+1:syfl%ixmaxg-1,:,:)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%bxy(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfl%bxy(ix,:,:)),ibuf),shape(syfl%bxy(ix,:,:)))
            end do
            do ix=syfl%ixmax+1,syfl%ixmaxg-1
              syfl%bxz(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfl%bxz(ix,:,:)),ibuf),shape(syfl%bxz(ix,:,:)))
            end do
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%byx(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfl%byx(ix,:,:)),ibuf),shape(syfl%byx(ix,:,:)))
            end do
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%byz(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfl%byz(ix,:,:)),ibuf),shape(syfl%byz(ix,:,:)))
            end do
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%bzx(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfl%bzx(ix,:,:)),ibuf),shape(syfl%bzx(ix,:,:)))
            end do
            do ix=syfl%ixmax,syfl%ixmaxg-1
              syfl%bzy(ix,:,:) =  reshape(mpi_unpack_real_array( size(syfl%bzy(ix,:,:)),ibuf),shape(syfl%bzy(ix,:,:)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndb_xrecv
#endif

subroutine em3d_exchange_bndf_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
integer(ISZ)   ::ix,ibuf
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          ! --- case lower yee, upper yee
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 1
            if (yfu%nxguard>1) then
              call mpi_packbuffer_init((yfu%nxguard-1)*size(yfu%f(0,:,:)),ibuf)
              do ix = yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
                call mympi_pack(yfu%f(ix,:,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,1,ibuf)
            end if
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 2
            if (yfl%nxguard>1) then
              call mpi_packbuffer_init((yfl%nxguard-1)*size(yfl%f(0,:,:)),ibuf)
              do ix = yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
                call mympi_pack(yfl%f(ix,:,:),ibuf)
              end do
              call mpi_isend_pack(fu%proc,2,ibuf)
            end if
          else
#endif
            yfl%f(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%f(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
            yfu%f(yfu%ixming+1:yfu%ixmin-1,:,:) = yfl%f(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          ! --- case lower yee, upper split yee
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%f(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = syfu%fx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                              + syfu%fy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:) &
                                              + syfu%fz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfu%fx(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%fy(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%fz(syfu%ixming+1:syfu%ixmin-1,:,:) = yfl%f(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
          ! --- case lower split yee, upper yee
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%f(yfu%ixming+1:yfu%ixmin-1,:,:)      = syfl%fx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                                   + syfl%fy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:) &
                                                   + syfl%fz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfl%fx(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%fy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%fz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = yfu%f(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
        case(splityeefield)
          ! --- case lower split yee, upper split yee
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            if (syfu%nxguard>1) then
              ibuf = 3
              call mpi_packbuffer_init( 3*int(size(syfu%fx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:))) ,ibuf)
              do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
                call mympi_pack(syfu%fx(ix,:,:),ibuf)
              end do
              do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
                call mympi_pack(syfu%fy(ix,:,:),ibuf)
              end do
              do ix=syfu%ixmin+1,syfu%ixmin+syfu%nxguard-1
                call mympi_pack(syfu%fz(ix,:,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,3,ibuf)
              mpireqpnt=mpireqpnt+1
            end if
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            if (syfl%nxguard>1) then
              ibuf = 4
              call mpi_packbuffer_init(3*int(size(syfl%fx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:))) ,ibuf)
              do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
                call mympi_pack(syfl%fx(ix,:,:),ibuf)
              end do
              do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
                call mympi_pack(syfl%fy(ix,:,:),ibuf)
              end do
              do ix=syfl%ixmax-syfl%nxguard+1,syfl%ixmax-1
                call mympi_pack(syfl%fz(ix,:,:),ibuf)
              end do
              call mpi_isend_pack(fu%proc,4,ibuf)
              mpireqpnt=mpireqpnt+1
            end if
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%fx(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%fx(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%fy(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%fy(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)
          syfu%fz(syfu%ixming+1:syfu%ixmin-1,:,:) = syfl%fz(syfl%ixmax-syfl%nxguard+1:syfl%ixmax-1,:,:)

          syfl%fx(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%fx(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%fy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%fy(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
          syfl%fz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = syfu%fz(syfu%ixmin+1:syfu%ixmin+syfu%nxguard-1,:,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndf_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndf_xrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            if (yfu%nxguard>1) then
              ibuf = 5
              call mpi_packbuffer_init((yfu%nxguard-1)*size(yfu%Ez(0,:,:)),ibuf)
              call mpi_recv_pack(fl%proc,2,ibuf)
              do ix = yfu%ixming+1,yfu%ixmin-1
                yfu%f(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%F(0,:,:)),ibuf),shape(yfu%F(0,:,:)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            if (yfl%nxguard>1) then
              ibuf = 6
              call mpi_packbuffer_init((yfl%nxguard-1)*size(yfl%ez(yfl%ixmin,:,:)),ibuf)
              call mpi_recv_pack(fu%proc,1,ibuf)
              do ix = yfl%ixmax+1,yfl%ixmaxg-1
                yfl%f(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%F(0,:,:)),ibuf),shape(yfl%F(0,:,:)))
              end do
            end if
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            if (syfu%nxguard>1) then
              ibuf = 7
              call mpi_packbuffer_init(3*size(syfu%ezx(syfu%ixming+1:syfu%ixmin-1,:,:)),ibuf)
              call mpi_recv_pack(fl%proc,4,ibuf)

              do ix = syfu%ixming+1,syfu%ixmin-1
                syfu%fx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%fx(ix,:,:)),ibuf),shape(syfu%fx(ix,:,:)))
              end do
              do ix = syfu%ixming+1,syfu%ixmin-1
                syfu%fy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%fy(ix,:,:)),ibuf),shape(syfu%fy(ix,:,:)))
              end do
              do ix = syfu%ixming+1,syfu%ixmin-1
                syfu%fz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfu%fz(ix,:,:)),ibuf),shape(syfu%fz(ix,:,:)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            if (syfl%nxguard>1) then
              ibuf = 8
              call mpi_packbuffer_init(3*size(syfl%ezx(syfl%ixmax+1:syfl%ixmaxg-1,:,:)),ibuf)
              call mpi_recv_pack(fu%proc,3,ibuf)
              do ix=syfl%ixmax+1,syfl%ixmaxg-1
                syfl%fx(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%fx(ix,:,:)),ibuf),shape(syfl%fx(ix,:,:)))
              end do
              do ix=syfl%ixmax+1,syfl%ixmaxg-1
                syfl%fy(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%fy(ix,:,:)),ibuf),shape(syfl%fy(ix,:,:)))
              end do
              do ix=syfl%ixmax+1,syfl%ixmaxg-1
                syfl%fz(ix,:,:) = reshape(mpi_unpack_real_array( size(syfl%fz(ix,:,:)),ibuf),shape(syfl%fz(ix,:,:)))
              end do
            end if
          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bndf_xrecv
#endif

subroutine em3d_exchange_bndj_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix,ibuf

#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then

            ! --- send data down in z
            ibuf = 17
            call mpi_packbuffer_init((3*yfu%nxguard+2)*size(yfu%J(0,:,:,1)),ibuf)
            do ix = yfu%ixming,yfu%ixmin-1
              call mympi_pack(yfu%J(ix,:,:,1),ibuf)
            end do
            do ix = yfu%ixming,yfu%ixmin
              call mympi_pack(yfu%J(ix,:,:,2),ibuf)
              call mympi_pack(yfu%J(ix,:,:,3),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            ibuf = 18
            call mpi_packbuffer_init((3*yfl%nxguard+2)*size(yfl%J(0,:,:,1)),ibuf)
            do ix = yfl%ixmax, yfl%ixmaxg-1
              call mympi_pack(yfl%J(ix,:,:,1),ibuf)
            end do
            do ix = yfl%ixmax, yfl%ixmaxg
              call mympi_pack(yfl%J(ix,:,:,2),ibuf)
              call mympi_pack(yfl%J(ix,:,:,3),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)

          else
#endif
          ix = yfu%nxguard
          yfu%J(0:ix-1,:,:,1)    = yfu%J(0:ix-1,:,:,1)    + yfl%J(yfl%nx:yfl%nx+ix-1,:,:,1) 
          yfu%J(0:ix,:,:,2:3)    = yfu%J(0:ix,:,:,2:3)    + yfl%J(yfl%nx:yfl%nx+ix,:,:,2:3)

          yfl%J(yfl%nx-ix:yfl%nx-1,:,:,:) = yfl%J(yfl%nx-ix:yfl%nx-1,:,:,:) + yfu%J(-ix:-1,:,:,:)
          yfl%J(yfl%nx,:,:,2:3) = yfu%J(0,:,:,2:3)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndj_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_xrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            ibuf = 19
            call mpi_packbuffer_init((3*yfu%nxguard+2)*size(yfu%J(0,:,:,1)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do ix = 0,yfu%nxguard-1
              yfu%J(ix,:,:,1) = yfu%J(ix,:,:,1) + reshape(mpi_unpack_real_array( size(yfu%J(0,:,:,1)),ibuf), &
                                                                                shape(yfu%J(0,:,:,1)))
            end do
            do ix = 0,yfu%nxguard
              yfu%J(ix,:,:,2) = yfu%J(ix,:,:,2) + reshape(mpi_unpack_real_array( size(yfu%J(0,:,:,1)),ibuf), &
                                                                                shape(yfu%J(0,:,:,1)))
              yfu%J(ix,:,:,3) = yfu%J(ix,:,:,3) + reshape(mpi_unpack_real_array( size(yfu%J(0,:,:,1)),ibuf), &
                                                                                shape(yfu%J(0,:,:,1)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            ibuf = 20
            call mpi_packbuffer_init((3*yfl%nxguard+2)*size(yfl%J(0,:,:,1)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do ix = -yfl%nxguard,-1
              yfl%J(yfl%nx+ix,:,:,1) = yfl%J(yfl%nx+ix,:,:,1) + reshape(mpi_unpack_real_array( size(yfl%J(yfl%nx-1,:,:,1)),ibuf),&
                                                                                              shape(yfl%J(yfl%nx-1,:,:,1)))
            end do
            do ix = -yfl%nxguard,0
              yfl%J(yfl%nx+ix,:,:,2) = yfl%J(yfl%nx+ix,:,:,2) + reshape(mpi_unpack_real_array( size(yfl%J(yfl%nx-1,:,:,2)),ibuf),&
                                                                                              shape(yfl%J(yfl%nx-1,:,:,2)))
              yfl%J(yfl%nx+ix,:,:,3) = yfl%J(yfl%nx+ix,:,:,3) + reshape(mpi_unpack_real_array( size(yfl%J(yfl%nx-1,:,:,3)),ibuf),&
                                                                                              shape(yfl%J(yfl%nx-1,:,:,3)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndj_xrecv
#endif

subroutine em3d_exchange_bndrho_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix,ibuf

#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ibuf=21
            ! --- send data down in z
            call mpi_packbuffer_init(size(yfu%rho(0:yfu%nxguard,:,:)),ibuf)
            do ix = -yfu%nxguard,0
              call mympi_pack(yfu%rho(ix,:,:),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            ibuf = 22
            call mpi_packbuffer_init(size(yfl%rho(0:yfl%nxguard,:,:)),ibuf)
            do ix = 0,yfl%nxguard
              call mympi_pack(yfl%rho(yfl%nx+ix,:,:),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)

          else
#endif
          ix = yfu%nxguard
          yfu%Rho(0:ix,:,:)      = yfu%Rho(0:ix,:,:)      + yfl%Rho(yfl%nx:yfl%nx+ix,:,:)
          yfl%Rho(yfl%nx-ix:yfl%nx-1,:,:) = yfl%Rho(yfl%nx-ix:yfl%nx-1,:,:) + yfu%Rho(-ix:-1,:,:)
          yfl%Rho(yfl%nx,:,:)   = yfu%Rho(0,:,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndrho_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_xrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix, ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            ibuf = 23
            call mpi_packbuffer_init(size(yfu%rho(0:yfu%nxguard,:,:)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do ix = 0,yfu%nxguard
              yfu%rho(ix,:,:) = yfu%rho(ix,:,:) + reshape(mpi_unpack_real_array( size(yfu%rho(0,:,:)),ibuf), &
                                                                                    shape(yfu%rho(0,:,:)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            ibuf = 24
            call mpi_packbuffer_init(size(yfl%rho(0:yfl%nxguard,:,:)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do ix = -yfl%nxguard,0
              yfl%rho(yfl%nx+ix,:,:) = yfl%rho(yfl%nx+ix,:,:) + reshape(mpi_unpack_real_array(size(yfl%rho(yfl%nx,:,:)),ibuf),&
                                                                                            shape(yfl%rho(yfl%nx,:,:)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndrho_xrecv
#endif

subroutine em3d_exchange_bnde_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
integer(ISZ)   ::iy,ibuf
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          ! --- case lower yee, upper yee
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 1
            call mpi_packbuffer_init((3*yfu%nyguard-2)*size(yfu%ez(:,0,:)),ibuf)
            if (yfu%nyguard>1) then
              do iy = yfu%iymin+1,yfu%iymin+yfu%nyguard-1
                call mympi_pack(yfu%ex(:,iy,:),ibuf)
              end do
            end if
            do iy = yfu%iymin,yfu%iymin+yfu%nyguard-1
              call mympi_pack(yfu%ey(:,iy,:),ibuf)
            end do
            if (yfu%nyguard>1) then
              do iy = yfu%iymin+1,yfu%iymin+yfu%nyguard-1
                call mympi_pack(yfu%ez(:,iy,:),ibuf)
              end do
            end if
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 2
            call mpi_packbuffer_init((3*yfl%nyguard-2)*size(yfl%ez(:,0,:)),ibuf)
            if (yfl%nyguard>1) then
              do iy = yfl%iymax-yfl%nyguard+1,yfl%iymax-1
                call mympi_pack(yfl%ex(:,iy,:),ibuf)
              end do
            end if
            do iy =yfl%iymax-yfl%nyguard,yfl%iymax-1
              call mympi_pack(yfl%ey(:,iy,:),ibuf)
            end do
            if (yfl%nyguard>1) then
              do iy = yfl%iymax-yfl%nyguard+1,yfl%iymax-1
                call mympi_pack(yfl%ez(:,iy,:),ibuf)
              end do
            end if
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
            yfl%ex(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%ex(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
            yfl%ey(:,yfl%iymax  :yfl%iymaxg-1,:) = yfu%ey(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)
            yfl%ez(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%ez(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)

            yfu%ex(:,yfu%iyming+1:yfu%iymin-1,:) = yfl%ex(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)
            yfu%ey(:,yfu%iyming  :yfu%iymin-1,:) = yfl%ey(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)
            yfu%ez(:,yfu%iyming+1:yfu%iymin-1,:) = yfl%ez(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          ! --- case lower yee, upper split yee
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%ex(:,yfl%iymax+1:yfl%iymaxg-1,:) = syfu%exx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                               + syfu%exy(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                               + syfu%exz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          yfl%ey(:,yfl%iymax  :yfl%iymaxg-1,:) = syfu%eyx(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:) &
                                               + syfu%eyy(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:) &
                                               + syfu%eyz(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          yfl%ez(:,yfl%iymax+1:yfl%iymaxg-1,:) = syfu%ezx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                               + syfu%ezy(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                               + syfu%ezz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfu%exx(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%exy(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%exz(:,syfu%iyming+1:syfu%iymin-1,:) = yfl%ex(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)
          syfu%eyx(:,syfu%iyming  :syfu%iymin-1,:) = 0.
          syfu%eyy(:,syfu%iyming  :syfu%iymin-1,:) = 0.
          syfu%eyz(:,syfu%iyming  :syfu%iymin-1,:) = yfl%ey(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)
          syfu%ezx(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%ezy(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%ezz(:,syfu%iyming+1:syfu%iymin-1,:) = yfl%ez(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)

      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
          ! --- case lower split yee, upper yee
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%ex(:,yfu%iyming+1:yfu%iymin-1,:)      = syfl%exx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                                    + syfl%exy(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                                    + syfl%exz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          yfu%ey(:,yfu%iyming  :yfu%iymin-1,:)      = syfl%eyx(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:) &
                                                    + syfl%eyy(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:) &
                                                    + syfl%eyz(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          yfu%ez(:,yfu%iyming+1:yfu%iymin-1,:)      = syfl%ezx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                                    + syfl%ezy(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                                    + syfl%ezz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfl%exx(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%exy(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%exz(:,syfl%iymax+1:syfl%iymaxg-1,:) = yfu%ex(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
          syfl%eyx(:,syfl%iymax  :syfl%iymaxg-1,:) = 0.
          syfl%eyy(:,syfl%iymax  :syfl%iymaxg-1,:) = 0.
          syfl%eyz(:,syfl%iymax  :syfl%iymaxg-1,:) = yfu%ey(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)
          syfl%ezx(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%ezy(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%ezz(:,syfl%iymax+1:syfl%iymaxg-1,:) = yfu%ez(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
        case(splityeefield)
          ! --- case lower split yee, upper split yee
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 3
            call mpi_packbuffer_init( 6*int(size(syfu%exx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:))) &
                                    + 3*int(size(syfu%ezx(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:))) ,ibuf)

            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%exx(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%exy(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%exz(:,iy,:),ibuf)
            end do

            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%eyx(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%eyy(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%eyz(:,iy,:),ibuf)
            end do

            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%ezx(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%ezy(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%ezz(:,iy,:),ibuf)
            end do

            call mpi_isend_pack(fl%proc,3,ibuf)

            mpireqpnt=mpireqpnt+1
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 4
            call mpi_packbuffer_init(6*int(size(syfl%ezx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:))) &
                                    +3*int(size(syfl%ezx(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:))),ibuf)

            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%exx(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%exy(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%exz(:,iy,:),ibuf)
            end do

            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%eyx(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%eyy(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%eyz(:,iy,:),ibuf)
            end do

            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%ezx(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%ezy(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%ezz(:,iy,:),ibuf)
            end do

            call mpi_isend_pack(fu%proc,4,ibuf)

            mpireqpnt=mpireqpnt+1
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%exx(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%exx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%exy(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%exy(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%exz(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%exz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%eyx(:,syfu%iyming  :syfu%iymin-1,:) = syfl%eyx(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfu%eyy(:,syfu%iyming  :syfu%iymin-1,:) = syfl%eyy(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfu%eyz(:,syfu%iyming  :syfu%iymin-1,:) = syfl%eyz(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfu%ezx(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%ezx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%ezy(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%ezy(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%ezz(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%ezz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)

          syfl%exx(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%exx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%exy(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%exy(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%exz(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%exz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%eyx(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%eyx(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          syfl%eyy(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%eyy(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          syfl%eyz(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%eyz(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          syfl%ezx(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%ezx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%ezy(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%ezy(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%ezz(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%ezz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bnde_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bnde_yrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 5
            call mpi_packbuffer_init((3*yfu%nyguard-2)*size(yfu%Ez(:,0,:)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            if (yfu%nyguard>1) then
              do iy = yfu%iyming+1,yfu%iymin-1
                yfu%ex(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,0,:)),ibuf),shape(yfu%Ez(:,0,:)))
              end do
            end if
            do iy = yfu%iyming,yfu%iymin-1
              yfu%ey(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,0,:)),ibuf),shape(yfu%Ez(:,0,:)))
            end do
            if (yfu%nyguard>1) then
              do iy = yfu%iyming+1,yfu%iymin-1
                yfu%ez(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,0,:)),ibuf),shape(yfu%Ez(:,0,:)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 6
            call mpi_packbuffer_init((3*yfl%nyguard-2)*size(yfl%ez(:,yfl%iymin,:)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            if (yfl%nyguard>1) then
              do iy = yfl%iymax+1,yfl%iymaxg-1
                yfl%ex(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,0,:)),ibuf),shape(yfl%Ez(:,0,:)))
              end do
            end if
            do iy = yfl%iymax,yfl%iymaxg-1
              yfl%ey(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,0,:)),ibuf),shape(yfl%Ez(:,0,:)))
            end do
            if (yfl%nyguard>1) then
              do iy = yfl%iymax+1,yfl%iymaxg-1
                yfl%ez(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,0,:)),ibuf),shape(yfl%Ez(:,0,:)))
              end do
            end if
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 7
            call mpi_packbuffer_init(6*size(syfu%ezx(:,syfu%iyming+1:syfu%iymin-1,:)) &
                                    +3*size(syfu%ezx(:,syfu%iyming  :syfu%iymin-1,:)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)

            do iy = syfu%iyming+1,syfu%iymin-1
              syfu%exx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%exx(:,iy,:)),ibuf),shape(syfu%exx(:,iy,:)))
            end do
            do iy = syfu%iyming+1,syfu%iymin-1
              syfu%exy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%exy(:,iy,:)),ibuf),shape(syfu%exy(:,iy,:)))
            end do
            do iy = syfu%iyming+1,syfu%iymin-1
              syfu%exz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%exz(:,iy,:)),ibuf),shape(syfu%exz(:,iy,:)))
            end do

            do iy = syfu%iyming,syfu%iymin-1
              syfu%eyx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%eyx(:,iy,:)),ibuf),shape(syfu%eyx(:,iy,:)))
            end do
            do iy = syfu%iyming,syfu%iymin-1
              syfu%eyy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%eyy(:,iy,:)),ibuf),shape(syfu%eyy(:,iy,:)))
            end do
            do iy = syfu%iyming,syfu%iymin-1
              syfu%eyz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%eyz(:,iy,:)),ibuf),shape(syfu%eyz(:,iy,:)))
            end do

            do iy = syfu%iyming+1,syfu%iymin-1
              syfu%ezx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%ezx(:,iy,:)),ibuf),shape(syfu%ezx(:,iy,:)))
            end do
            do iy = syfu%iyming+1,syfu%iymin-1
              syfu%ezy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%ezy(:,iy,:)),ibuf),shape(syfu%ezy(:,iy,:)))
            end do
            do iy = syfu%iyming+1,syfu%iymin-1
              syfu%ezz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%ezz(:,iy,:)),ibuf),shape(syfu%ezz(:,iy,:)))
            end do

          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 8
            call mpi_packbuffer_init(6*size(syfl%ezx(:,syfl%iymax+1:syfl%iymaxg-1,:)) &
                                    +3*size(syfl%ezx(:,syfl%iymax  :syfl%iymaxg-1,:)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%exx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%exx(:,iy,:)),ibuf),shape(syfl%exx(:,iy,:)))
            end do
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%exy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%exy(:,iy,:)),ibuf),shape(syfl%exy(:,iy,:)))
            end do
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%exz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%exz(:,iy,:)),ibuf),shape(syfl%exz(:,iy,:)))
            end do
            
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%eyx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%eyx(:,iy,:)),ibuf),shape(syfl%eyx(:,iy,:)))
            end do
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%eyy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%eyy(:,iy,:)),ibuf),shape(syfl%eyy(:,iy,:)))
            end do
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%eyz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%eyz(:,iy,:)),ibuf),shape(syfl%eyz(:,iy,:)))
            end do
            
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%ezx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%ezx(:,iy,:)),ibuf),shape(syfl%ezx(:,iy,:)))
            end do
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%ezy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%ezy(:,iy,:)),ibuf),shape(syfl%ezy(:,iy,:)))
            end do
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%ezz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%ezz(:,iy,:)),ibuf),shape(syfl%ezz(:,iy,:)))
            end do

          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bnde_yrecv
#endif

subroutine em3d_exchange_bndb_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
integer(ISZ) :: ibuf,iy
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 9
            call mpi_packbuffer_init((3*yfu%nyguard-1)*size(yfu%by(:,0,:)),ibuf)
            do iy=yfu%iymin,yfu%iymin+yfu%nyguard-1
              call mympi_pack(yfu%bx(:,iy,:),ibuf)
            end do
            do iy=yfu%iymin+1, yfu%iymin+yfu%nyguard-1
              call mympi_pack(yfu%by(:,iy,:),ibuf)
            end do
            do iy=yfu%iymin,yfu%iymin+yfu%nyguard-1
              call mympi_pack(yfu%bz(:,iy,:),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 10
            call mpi_packbuffer_init((3*yfl%nyguard-1)*size(yfl%by(:,0,:)),ibuf)
            do iy=yfl%iymax-yfl%nyguard,yfl%iymax-1
              call mympi_pack(yfl%bx(:,iy,:),ibuf)
            end do
            do iy=yfl%iymax-yfl%nyguard+1,yfl%iymax-1
              call mympi_pack(yfl%by(:,iy,:),ibuf)
            end do
            do iy=yfl%iymax-yfl%nyguard,yfl%iymax-1
              call mympi_pack(yfl%bz(:,iy,:),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
          yfl%bx(:,yfl%iymax  :yfl%iymaxg-1,:) = yfu%bx(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)
          yfl%by(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%by(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
          yfl%bz(:,yfl%iymax  :yfl%iymaxg-1,:) = yfu%bz(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)

          yfu%bx(:,yfu%iyming  :yfu%iymin-1,:) = yfl%bx(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)
          yfu%by(:,yfu%iyming+1:yfu%iymin-1,:) = yfl%by(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)
          yfu%bz(:,yfu%iyming  :yfu%iymin-1,:) = yfl%bz(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)
          
#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=fu%proc) return
          yfl%bx(:,yfl%iymax  :yfl%iymaxg-1,:) = (syfu%bxy(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:) &
                                               +  syfu%bxz(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:))/syfu%clight
          yfl%by(:,yfl%iymax+1:yfl%iymaxg-1,:) = (syfu%byx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                               +  syfu%byz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:))/syfu%clight
          yfl%bz(:,yfl%iymax  :yfl%iymaxg-1,:) = (syfu%bzx(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:) &
                                               +  syfu%bzy(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:))/syfu%clight
          syfu%bxy(:,syfu%iyming  :syfu%iymin-1,:) = yfl%bx(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)*syfu%clight
          syfu%bxz(:,syfu%iyming  :syfu%iymin-1,:) = 0.
          syfu%byx(:,syfu%iyming+1:syfu%iymin-1,:) = yfl%by(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)*syfu%clight
          syfu%byz(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%bzx(:,syfu%iyming  :syfu%iymin-1,:) = yfl%bz(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)*syfu%clight
          syfu%bzy(:,syfu%iyming  :syfu%iymin-1,:) = 0.
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=fu%proc) return
          yfu%bx(:,yfu%iyming  :yfu%iymin-1,:) = (syfl%bxy(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:) &
                                               +  syfl%bxz(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:))/syfl%clight
          yfu%by(:,yfu%iyming+1:yfu%iymin-1,:) = (syfl%byx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                               +  syfl%byz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:))/syfl%clight
          yfu%bz(:,yfu%iyming  :yfu%iymin-1,:) = (syfl%bzx(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:) &
                                               +  syfl%bzy(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:))/syfl%clight
          syfl%bxy(:,syfl%iymax  :syfl%iymaxg-1,:) = yfu%bx(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)*syfl%clight
          syfl%bxz(:,syfl%iymax  :syfl%iymaxg-1,:) = 0.
          syfl%byx(:,syfl%iymax+1:syfl%iymaxg-1,:) = yfu%by(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)*syfl%clight
          syfl%byz(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%bzx(:,syfl%iymax  :syfl%iymaxg-1,:) = yfu%bz(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)*syfl%clight
          syfl%bzy(:,syfl%iymax  :syfl%iymaxg-1,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 11
            call mpi_packbuffer_init(4*size(syfu%byx(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)) &
                                    +2*size(syfu%byx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)),ibuf)

            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%bxy(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%bxz(:,iy,:),ibuf)
            end do
            
            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%byx(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%byz(:,iy,:),ibuf)
            end do
            
            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%bzx(:,iy,:),ibuf)
            end do
            do iy=syfu%iymin,syfu%iymin+syfu%nyguard-1
              call mympi_pack(syfu%bzy(:,iy,:),ibuf)
            end do

            call mpi_isend_pack(fl%proc,3,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 12
            call mpi_packbuffer_init(4*size(syfl%byx(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)) &
                                    +2*size(syfl%byx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)),ibuf)

            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%bxy(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%bxz(:,iy,:),ibuf)
            end do

            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%byx(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
              call mympi_pack(syfl%byz(:,iy,:),ibuf)
            end do

            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%bzx(:,iy,:),ibuf)
            end do
            do iy=syfl%iymax-syfl%nyguard,syfl%iymax-1
              call mympi_pack(syfl%bzy(:,iy,:),ibuf)
            end do

            call mpi_isend_pack(fu%proc,4,ibuf)
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%bxy(:,syfu%iyming  :syfu%iymin-1,:) = syfl%bxy(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfu%bxz(:,syfu%iyming  :syfu%iymin-1,:) = syfl%bxz(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfu%byx(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%byx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%byz(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%byz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%bzx(:,syfu%iyming  :syfu%iymin-1,:) = syfl%bzx(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfu%bzy(:,syfu%iyming  :syfu%iymin-1,:) = syfl%bzy(:,syfl%iymax-syfl%nyguard  :syfl%iymax-1,:)
          syfl%bxy(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%bxy(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          syfl%bxz(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%bxz(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          syfl%byx(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%byx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%byz(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%byz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%bzx(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%bzx(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
          syfl%bzy(:,syfl%iymax  :syfl%iymaxg-1,:) = syfu%bzy(:,syfu%iymin  :syfu%iymin+syfu%nyguard-1,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndb_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndb_yrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ibuf,iy

          if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 13
            call mpi_packbuffer_init((3*yfu%nyguard-1)*size(yfu%bx(:,0,:)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iy=yfu%iyming,yfu%iymin-1
              yfu%bx(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%bx(:,iy,:)),ibuf),shape(yfu%bx(:,iy,:)))
            end do
            do iy=yfu%iyming+1,yfu%iymin-1
              yfu%by(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%by(:,iy,:)),ibuf),shape(yfu%by(:,iy,:)))
            end do
            do iy=yfu%iyming,yfu%iymin-1
              yfu%bz(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%bz(:,iy,:)),ibuf),shape(yfu%bz(:,iy,:)))
            end do
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 14
            call mpi_packbuffer_init((3*yfl%nyguard-1)*size(yfl%bx(:,0,:)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iy=yfl%iymax,yfl%iymaxg-1
              yfl%bx(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%bx(:,iy,:)),ibuf),shape(yfl%bx(:,iy,:)))
            end do
            do iy=yfl%iymax+1,yfl%iymaxg-1
              yfl%by(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%by(:,iy,:)),ibuf),shape(yfl%by(:,iy,:)))
            end do
            do iy=yfl%iymax,yfl%iymaxg-1
              yfl%bz(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%bz(:,iy,:)),ibuf),shape(yfl%bz(:,iy,:)))
            end do
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ibuf = 15
            call mpi_packbuffer_init(4*size(syfu%bxy(:,syfu%iyming  :syfu%iymin-1,:)) &
                                    +2*size(syfu%bxy(:,syfu%iyming+1:syfu%iymin-1,:)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)
            ! --- recv data from down in z
            do iy=syfu%iyming,syfu%iymin-1
              syfu%bxy(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfu%bxy(:,iy,:)),ibuf),shape(syfu%bxy(:,iy,:)))
            end do
            do iy=syfu%iyming,syfu%iymin-1
              syfu%bxz(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfu%bxz(:,iy,:)),ibuf),shape(syfu%bxz(:,iy,:)))
            end do
            do iy=syfu%iyming+1,syfu%iymin-1
              syfu%byx(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfu%byx(:,iy,:)),ibuf),shape(syfu%byx(:,iy,:)))
            end do
            do iy=syfu%iyming+1,syfu%iymin-1
              syfu%byz(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfu%byz(:,iy,:)),ibuf),shape(syfu%byz(:,iy,:)))
            end do
            do iy=syfu%iyming,syfu%iymin-1
              syfu%bzx(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfu%bzx(:,iy,:)),ibuf),shape(syfu%bzx(:,iy,:)))
            end do
            do iy=syfu%iyming,syfu%iymin-1
              syfu%bzy(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfu%bzy(:,iy,:)),ibuf),shape(syfu%bzy(:,iy,:)))
            end do
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 16
            call mpi_packbuffer_init(4*size(syfl%bxy(:,syfl%iymax  :syfl%iymaxg-1,:)) &
                                    +2*size(syfl%bxy(:,syfl%iymax+1:syfl%iymaxg-1,:)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%bxy(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfl%bxy(:,iy,:)),ibuf),shape(syfl%bxy(:,iy,:)))
            end do
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%bxz(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfl%bxz(:,iy,:)),ibuf),shape(syfl%bxz(:,iy,:)))
            end do
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%byx(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfl%byx(:,iy,:)),ibuf),shape(syfl%byx(:,iy,:)))
            end do
            do iy=syfl%iymax+1,syfl%iymaxg-1
              syfl%byz(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfl%byz(:,iy,:)),ibuf),shape(syfl%byz(:,iy,:)))
            end do
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%bzx(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfl%bzx(:,iy,:)),ibuf),shape(syfl%bzx(:,iy,:)))
            end do
            do iy=syfl%iymax,syfl%iymaxg-1
              syfl%bzy(:,iy,:) =  reshape(mpi_unpack_real_array( size(syfl%bzy(:,iy,:)),ibuf),shape(syfl%bzy(:,iy,:)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndb_yrecv
#endif

subroutine em3d_exchange_bndf_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
integer(ISZ)   ::iy,ibuf
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          ! --- case lower yee, upper yee
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 1
            if (yfu%nyguard>1) then
              call mpi_packbuffer_init((yfu%nyguard-1)*size(yfu%f(:,0,:)),ibuf)
              do iy = yfu%iymin+1,yfu%iymin+yfu%nyguard-1
                call mympi_pack(yfu%f(:,iy,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,1,ibuf)
            end if
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 2
            if (yfl%nyguard>1) then
              call mpi_packbuffer_init((yfl%nyguard-1)*size(yfl%f(:,0,:)),ibuf)
              do iy = yfl%iymax-yfl%nyguard+1,yfl%iymax-1
                call mympi_pack(yfl%f(:,iy,:),ibuf)
              end do
              call mpi_isend_pack(fu%proc,2,ibuf)
            end if
          else
#endif
            yfl%f(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%f(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
            yfu%f(:,yfu%iyming+1:yfu%iymin-1,:) = yfl%f(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          ! --- case lower yee, upper split yee
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%f(:,yfl%iymax+1:yfl%iymaxg-1,:) = syfu%fx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                              + syfu%fy(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:) &
                                              + syfu%fz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfu%fx(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%fy(:,syfu%iyming+1:syfu%iymin-1,:) = 0.
          syfu%fz(:,syfu%iyming+1:syfu%iymin-1,:) = yfl%f(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)

      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
          ! --- case lower split yee, upper yee
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%f(:,yfu%iyming+1:yfu%iymin-1,:)      = syfl%fx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                                   + syfl%fy(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:) &
                                                   + syfl%fz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfl%fx(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%fy(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%fz(:,syfl%iymax+1:syfl%iymaxg-1,:) = yfu%f(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
        case(splityeefield)
          ! --- case lower split yee, upper split yee
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            if (syfu%nyguard>1) then
              ibuf = 3
              call mpi_packbuffer_init( 3*int(size(syfu%fx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:))) ,ibuf)
              do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
                call mympi_pack(syfu%fx(:,iy,:),ibuf)
              end do
              do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
                call mympi_pack(syfu%fy(:,iy,:),ibuf)
              end do
              do iy=syfu%iymin+1,syfu%iymin+syfu%nyguard-1
                call mympi_pack(syfu%fz(:,iy,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,3,ibuf)
              mpireqpnt=mpireqpnt+1
            end if
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            if (syfl%nyguard>1) then
              ibuf = 4
              call mpi_packbuffer_init(3*int(size(syfl%fx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:))) ,ibuf)
              do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
                call mympi_pack(syfl%fx(:,iy,:),ibuf)
              end do
              do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
                call mympi_pack(syfl%fy(:,iy,:),ibuf)
              end do
              do iy=syfl%iymax-syfl%nyguard+1,syfl%iymax-1
                call mympi_pack(syfl%fz(:,iy,:),ibuf)
              end do
              call mpi_isend_pack(fu%proc,4,ibuf)
              mpireqpnt=mpireqpnt+1
            end if
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%fx(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%fx(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%fy(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%fy(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)
          syfu%fz(:,syfu%iyming+1:syfu%iymin-1,:) = syfl%fz(:,syfl%iymax-syfl%nyguard+1:syfl%iymax-1,:)

          syfl%fx(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%fx(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%fy(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%fy(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
          syfl%fz(:,syfl%iymax+1:syfl%iymaxg-1,:) = syfu%fz(:,syfu%iymin+1:syfu%iymin+syfu%nyguard-1,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndf_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndf_yrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            if (yfu%nyguard>1) then
              ibuf = 5
              call mpi_packbuffer_init((yfu%nyguard-1)*size(yfu%Ez(:,0,:)),ibuf)
              call mpi_recv_pack(fl%proc,2,ibuf)
              do iy = yfu%iyming+1,yfu%iymin-1
                yfu%f(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%F(:,0,:)),ibuf),shape(yfu%F(:,0,:)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            if (yfl%nyguard>1) then
              ibuf = 6
              call mpi_packbuffer_init((yfl%nyguard-1)*size(yfl%ez(:,yfl%iymin,:)),ibuf)
              call mpi_recv_pack(fu%proc,1,ibuf)
              do iy = yfl%iymax+1,yfl%iymaxg-1
                yfl%f(:,iy,:) = reshape(mpi_unpack_real_array( size(yfl%F(:,0,:)),ibuf),shape(yfl%F(:,0,:)))
              end do
            end if
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            if (syfu%nyguard>1) then
              ibuf = 7
              call mpi_packbuffer_init(3*size(syfu%ezx(:,syfu%iyming+1:syfu%iymin-1,:)),ibuf)
              call mpi_recv_pack(fl%proc,4,ibuf)

              do iy = syfu%iyming+1,syfu%iymin-1
                syfu%fx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%fx(:,iy,:)),ibuf),shape(syfu%fx(:,iy,:)))
              end do
              do iy = syfu%iyming+1,syfu%iymin-1
                syfu%fy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%fy(:,iy,:)),ibuf),shape(syfu%fy(:,iy,:)))
              end do
              do iy = syfu%iyming+1,syfu%iymin-1
                syfu%fz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfu%fz(:,iy,:)),ibuf),shape(syfu%fz(:,iy,:)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            if (syfl%nyguard>1) then
              ibuf = 8
              call mpi_packbuffer_init(3*size(syfl%ezx(:,syfl%iymax+1:syfl%iymaxg-1,:)),ibuf)
              call mpi_recv_pack(fu%proc,3,ibuf)
              do iy=syfl%iymax+1,syfl%iymaxg-1
                syfl%fx(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%fx(:,iy,:)),ibuf),shape(syfl%fx(:,iy,:)))
              end do
              do iy=syfl%iymax+1,syfl%iymaxg-1
                syfl%fy(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%fy(:,iy,:)),ibuf),shape(syfl%fy(:,iy,:)))
              end do
              do iy=syfl%iymax+1,syfl%iymaxg-1
                syfl%fz(:,iy,:) = reshape(mpi_unpack_real_array( size(syfl%fz(:,iy,:)),ibuf),shape(syfl%fz(:,iy,:)))
              end do
            end if
          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bndf_yrecv
#endif

subroutine em3d_exchange_bndj_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy,ibuf

#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then

            ! --- send data down in z
            ibuf = 17
            call mpi_packbuffer_init((3*yfu%nyguard+2)*size(yfu%J(:,-1,:,1)),ibuf)
            do iy = yfu%iyming,yfu%iymin
              call mympi_pack(yfu%J(:,iy,:,1),ibuf)
            end do
            do iy = yfu%iyming,yfu%iymin-1
              call mympi_pack(yfu%J(:,iy,:,2),ibuf)
            end do
            do iy = yfu%iyming,yfu%iymin
              call mympi_pack(yfu%J(:,iy,:,3),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            ibuf = 18
            call mpi_packbuffer_init((3*yfl%nyguard+2)*size(yfl%J(:,0,:,1)),ibuf)
            do iy = yfl%iymax, yfl%iymaxg
              call mympi_pack(yfl%J(:,iy,:,1),ibuf)
            end do
            do iy = yfl%iymax, yfl%iymaxg-1
              call mympi_pack(yfl%J(:,iy,:,2),ibuf)
            end do
            do iy = yfl%iymax, yfl%iymaxg
              call mympi_pack(yfl%J(:,iy,:,3),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)

          else
#endif
          iy = yfu%nyguard
          yfu%J(:,0:iy  ,:,1:3:2) = yfu%J(:,0:iy,  :,1:3:2) + yfl%J(:,yfl%ny:yfl%ny+iy  ,:,1:3:2)
          yfu%J(:,0:iy-1,:,2    ) = yfu%J(:,0:iy-1,:,2    ) + yfl%J(:,yfl%ny:yfl%ny+iy-1,:,2    ) 

          yfl%J(:,yfl%ny-iy:yfl%ny-1,:,:) = yfl%J(:,yfl%ny-iy:yfl%ny-1,:,:) + yfu%J(:,-iy:-1,:,:)
          yfl%J(:,yfl%ny,:,1:3:2) = yfu%J(:,0,:,1:3:2)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndj_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_yrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            ibuf = 19
            call mpi_packbuffer_init((3*yfu%nyguard+2)*size(yfu%J(:,0,:,1)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iy = 0,yfu%nyguard
              yfu%J(:,iy,:  ,1) = yfu%J(:,iy,:  ,1) + reshape(mpi_unpack_real_array( size(yfu%J(:,0,:,1)),ibuf), &
                                                                                    shape(yfu%J(:,0,:,1)))
            end do
            do iy = 0,yfu%nyguard-1
              yfu%J(:,iy,:  ,2) = yfu%J(:,iy,:  ,2) + reshape(mpi_unpack_real_array( size(yfu%J(:,0,:,1)),ibuf), &
                                                                                    shape(yfu%J(:,0,:,1)))
            end do
            do iy = 0,yfu%nyguard
              yfu%J(:,iy,:  ,3) = yfu%J(:,iy,:  ,3) + reshape(mpi_unpack_real_array( size(yfu%J(:,0,:,1)),ibuf), &
                                                                                    shape(yfu%J(:,0,:,1)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            ibuf = 20
            call mpi_packbuffer_init((3*yfl%nyguard+2)*size(yfl%J(:,0,:,1)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iy = -yfl%nyguard,0
              yfl%J(:,yfl%ny+iy,:,1) = yfl%J(:,yfl%ny+iy,:,1) + reshape(mpi_unpack_real_array( size(yfl%J(:,yfl%ny-1,:,1)),ibuf),&
                                                                                              shape(yfl%J(:,yfl%ny-1,:,1)))
            end do
            do iy = -yfl%nyguard,-1
              yfl%J(:,yfl%ny+iy,:,2) = yfl%J(:,yfl%ny+iy,:,2) + reshape(mpi_unpack_real_array( size(yfl%J(:,yfl%ny-1,:,2)),ibuf),&
                                                                                              shape(yfl%J(:,yfl%ny-1,:,2)))
            end do
            do iy = -yfl%nyguard,0
              yfl%J(:,yfl%ny+iy,:,3) = yfl%J(:,yfl%ny+iy,:,3) + reshape(mpi_unpack_real_array( size(yfl%J(:,yfl%ny-1,:,3)),ibuf),&
                                                                                              shape(yfl%J(:,yfl%ny-1,:,3)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndj_yrecv
#endif

subroutine em3d_exchange_bndrho_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy,ibuf

#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ibuf=21
            ! --- send data down in y
            call mpi_packbuffer_init(size(yfu%rho(:,0:yfu%nyguard,:)),ibuf)
            do iy = -yfu%nyguard,0
              call mympi_pack(yfu%rho(:,iy,:),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in y
            ibuf = 22
            call mpi_packbuffer_init(size(yfl%rho(:,0:yfl%nyguard,:)),ibuf)
            do iy = 0,yfl%nyguard
              call mympi_pack(yfl%rho(:,yfl%ny+iy,:),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)

          else
#endif
          iy = yfu%nyguard
          yfu%Rho(:,0:iy,:)      = yfu%Rho(:,0:iy,:)      + yfl%Rho(:,yfl%ny:yfl%ny+iy,:)
          yfl%Rho(:,yfl%ny-iy:yfl%ny-1,:) = yfl%Rho(:,yfl%ny-iy:yfl%ny-1,:) + yfu%Rho(:,-iy:-1,:)
          yfl%Rho(:,yfl%ny,:)   = yfu%Rho(:,0,:)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndrho_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_yrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy, ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            ibuf = 23
            call mpi_packbuffer_init(size(yfu%rho(:,0:yfu%nyguard,:)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iy = 0,yfu%nyguard
              yfu%rho(:,iy,:) = yfu%rho(:,iy,:) + reshape(mpi_unpack_real_array( size(yfu%rho(:,0,:)),ibuf), &
                                                                                shape(yfu%rho(:,0,:)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            ibuf = 24
            call mpi_packbuffer_init(size(yfl%rho(:,0:yfl%nyguard,:)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iy = -yfl%nyguard,0
              yfl%rho(:,yfl%ny+iy,:) = yfl%rho(:,yfl%ny+iy,:) + reshape(mpi_unpack_real_array(size(yfl%rho(:,yfl%ny,:)),ibuf),&
                                                                                             shape(yfl%rho(:,yfl%ny,:)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndrho_yrecv
#endif

subroutine em3d_exchange_bnde_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
integer(ISZ)   ::iz,ibuf
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          ! --- case lower yee, upper yee
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 1
            call mpi_packbuffer_init((3*yfu%nzguard-2)*size(yfu%ez(:,:,0)),ibuf)
            do iz = yfu%izmin,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%ez(:,:,iz),ibuf)
            end do
            if (yfu%nzguard>1) then
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                call mympi_pack(yfu%ex(:,:,iz),ibuf)
              end do
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                call mympi_pack(yfu%ey(:,:,iz),ibuf)
              end do
            end if
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 2
            call mpi_packbuffer_init((3*yfl%nzguard-2)*size(yfl%ez(:,:,0)),ibuf)
            do iz =yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%ez(:,:,iz),ibuf)
            end do
            if (yfl%nzguard>1) then
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                call mympi_pack(yfl%ex(:,:,iz),ibuf)
              end do
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                call mympi_pack(yfl%ey(:,:,iz),ibuf)
              end do
            end if
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
            yfl%ex(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%ex(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
            yfl%ey(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%ey(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
            yfl%ez(:,:,yfl%izmax  :yfl%izmaxg-1) = yfu%ez(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)

            yfu%ex(:,:,yfu%izming+1:yfu%izmin-1) = yfl%ex(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
            yfu%ey(:,:,yfu%izming+1:yfu%izmin-1) = yfl%ey(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
            yfu%ez(:,:,yfu%izming  :yfu%izmin-1) = yfl%ez(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          ! --- case lower yee, upper split yee
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%ex(:,:,yfl%izmax+1:yfl%izmaxg-1) = syfu%exx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                               + syfu%exy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                               + syfu%exz(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          yfl%ey(:,:,yfl%izmax+1:yfl%izmaxg-1) = syfu%eyx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                               + syfu%eyy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                               + syfu%eyz(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          yfl%ez(:,:,yfl%izmax  :yfl%izmaxg-1) = syfu%ezx(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1) &
                                               + syfu%ezy(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1) &
                                               + syfu%ezz(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfu%exx(:,:,syfu%izming+1:syfu%izmin-1) = 0.
          syfu%exy(:,:,syfu%izming+1:syfu%izmin-1) = 0.
          syfu%exz(:,:,syfu%izming+1:syfu%izmin-1) = yfl%ex(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
          syfu%eyx(:,:,syfu%izming+1:syfu%izmin-1) = 0.
          syfu%eyy(:,:,syfu%izming+1:syfu%izmin-1) = 0.
          syfu%eyz(:,:,syfu%izming+1:syfu%izmin-1) = yfl%ey(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
          syfu%ezx(:,:,syfu%izming  :syfu%izmin-1) = 0.
          syfu%ezy(:,:,syfu%izming  :syfu%izmin-1) = 0.
          syfu%ezz(:,:,syfu%izming  :syfu%izmin-1) = yfl%ez(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)

      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
          ! --- case lower split yee, upper yee
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%ex(:,:,yfu%izming+1:yfu%izmin-1)      = syfl%exx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                                    + syfl%exy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                                    + syfl%exz(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          yfu%ey(:,:,yfu%izming+1:yfu%izmin-1)      = syfl%eyx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                                    + syfl%eyy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                                    + syfl%eyz(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          yfu%ez(:,:,yfu%izming  :yfu%izmin-1)      = syfl%ezx(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1) &
                                                    + syfl%ezy(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1) &
                                                    + syfl%ezz(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfl%exx(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
          syfl%exy(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
          syfl%exz(:,:,syfl%izmax+1:syfl%izmaxg-1) = yfu%ex(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
          syfl%eyx(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
          syfl%eyy(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
          syfl%eyz(:,:,syfl%izmax+1:syfl%izmaxg-1) = yfu%ey(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
          syfl%ezx(:,:,syfl%izmax  :syfl%izmaxg-1) = 0.
          syfl%ezy(:,:,syfl%izmax  :syfl%izmaxg-1) = 0.
          syfl%ezz(:,:,syfl%izmax  :syfl%izmaxg-1) = yfu%ez(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)
        case(splityeefield)
          ! --- case lower split yee, upper split yee
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 3
            call mpi_packbuffer_init( 6*int(size(syfu%exx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1))) &
                                    + 3*int(size(syfu%ezx(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1))) ,ibuf)

            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%exx(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%exy(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%exz(:,:,iz),ibuf)
            end do

            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%eyx(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%eyy(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%eyz(:,:,iz),ibuf)
            end do

            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%ezx(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%ezy(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%ezz(:,:,iz),ibuf)
            end do

            call mpi_isend_pack(fl%proc,3,ibuf)

            mpireqpnt=mpireqpnt+1
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 4
            call mpi_packbuffer_init(6*int(size(syfl%ezx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1))) &
                                    +3*int(size(syfl%ezx(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1))),ibuf)

            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%exx(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%exy(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%exz(:,:,iz),ibuf)
            end do

            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%eyx(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%eyy(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%eyz(:,:,iz),ibuf)
            end do

            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%ezx(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%ezy(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%ezz(:,:,iz),ibuf)
            end do

            call mpi_isend_pack(fu%proc,4,ibuf)

            mpireqpnt=mpireqpnt+1
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%exx(:,:,syfu%izming+1:syfu%izmin-1) = syfl%exx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%exy(:,:,syfu%izming+1:syfu%izmin-1) = syfl%exy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%exz(:,:,syfu%izming+1:syfu%izmin-1) = syfl%exz(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%eyx(:,:,syfu%izming+1:syfu%izmin-1) = syfl%eyx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%eyy(:,:,syfu%izming+1:syfu%izmin-1) = syfl%eyy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%eyz(:,:,syfu%izming+1:syfu%izmin-1) = syfl%eyz(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%ezx(:,:,syfu%izming  :syfu%izmin-1) = syfl%ezx(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfu%ezy(:,:,syfu%izming  :syfu%izmin-1) = syfl%ezy(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfu%ezz(:,:,syfu%izming  :syfu%izmin-1) = syfl%ezz(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)

          syfl%exx(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%exx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%exy(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%exy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%exz(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%exz(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%eyx(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%eyx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%eyy(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%eyy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%eyz(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%eyz(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%ezx(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%ezx(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfl%ezy(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%ezy(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfl%ezz(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%ezz(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bnde_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bnde_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 5
            call mpi_packbuffer_init((3*yfu%nzguard-2)*size(yfu%Ez(:,:,0)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iz = yfu%izming,yfu%izmin-1
              yfu%ez(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),ibuf),shape(yfu%Ez(:,:,0)))
            end do
            if (yfu%nzguard>1) then
              do iz = yfu%izming+1,yfu%izmin-1
                yfu%ex(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),ibuf),shape(yfu%Ez(:,:,0)))
              end do
              do iz = yfu%izming+1,yfu%izmin-1
                yfu%ey(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),ibuf),shape(yfu%Ez(:,:,0)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 6
            call mpi_packbuffer_init((3*yfl%nzguard-2)*size(yfl%ez(:,:,yfl%izmin)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iz = yfl%izmax,yfl%izmaxg-1
              yfl%ez(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),ibuf),shape(yfl%Ez(:,:,0)))
            end do
            if (yfl%nzguard>1) then
              do iz = yfl%izmax+1,yfl%izmaxg-1
                yfl%ex(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),ibuf),shape(yfl%Ez(:,:,0)))
              end do
              do iz = yfl%izmax+1,yfl%izmaxg-1
                yfl%ey(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),ibuf),shape(yfl%Ez(:,:,0)))
              end do
            end if
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 7
            call mpi_packbuffer_init(6*size(syfu%ezx(:,:,syfu%izming+1:syfu%izmin-1)) &
                                    +3*size(syfu%ezx(:,:,syfu%izming  :syfu%izmin-1)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)

            do iz = syfu%izming+1,syfu%izmin-1
              syfu%exx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%exx(:,:,iz)),ibuf),shape(syfu%exx(:,:,iz)))
            end do
            do iz = syfu%izming+1,syfu%izmin-1
              syfu%exy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%exy(:,:,iz)),ibuf),shape(syfu%exy(:,:,iz)))
            end do
            do iz = syfu%izming+1,syfu%izmin-1
              syfu%exz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%exz(:,:,iz)),ibuf),shape(syfu%exz(:,:,iz)))
            end do

            do iz = syfu%izming+1,syfu%izmin-1
              syfu%eyx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%eyx(:,:,iz)),ibuf),shape(syfu%eyx(:,:,iz)))
            end do
            do iz = syfu%izming+1,syfu%izmin-1
              syfu%eyy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%eyy(:,:,iz)),ibuf),shape(syfu%eyy(:,:,iz)))
            end do
            do iz = syfu%izming+1,syfu%izmin-1
              syfu%eyz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%eyz(:,:,iz)),ibuf),shape(syfu%eyz(:,:,iz)))
            end do

            do iz = syfu%izming,syfu%izmin-1
              syfu%ezx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%ezx(:,:,iz)),ibuf),shape(syfu%ezx(:,:,iz)))
            end do
            do iz = syfu%izming,syfu%izmin-1
              syfu%ezy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%ezy(:,:,iz)),ibuf),shape(syfu%ezy(:,:,iz)))
            end do
            do iz = syfu%izming,syfu%izmin-1
              syfu%ezz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%ezz(:,:,iz)),ibuf),shape(syfu%ezz(:,:,iz)))
            end do

          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 8
            call mpi_packbuffer_init(6*size(syfl%ezx(:,:,syfl%izmax+1:syfl%izmaxg-1)) &
                                    +3*size(syfl%ezx(:,:,syfl%izmax  :syfl%izmaxg-1)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%exx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%exx(:,:,iz)),ibuf),shape(syfl%exx(:,:,iz)))
            end do
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%exy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%exy(:,:,iz)),ibuf),shape(syfl%exy(:,:,iz)))
            end do
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%exz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%exz(:,:,iz)),ibuf),shape(syfl%exz(:,:,iz)))
            end do
            
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%eyx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%eyx(:,:,iz)),ibuf),shape(syfl%eyx(:,:,iz)))
            end do
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%eyy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%eyy(:,:,iz)),ibuf),shape(syfl%eyy(:,:,iz)))
            end do
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%eyz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%eyz(:,:,iz)),ibuf),shape(syfl%eyz(:,:,iz)))
            end do
            
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%ezx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%ezx(:,:,iz)),ibuf),shape(syfl%ezx(:,:,iz)))
            end do
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%ezy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%ezy(:,:,iz)),ibuf),shape(syfl%ezy(:,:,iz)))
            end do
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%ezz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%ezz(:,:,iz)),ibuf),shape(syfl%ezz(:,:,iz)))
            end do

          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bnde_zrecv
#endif

subroutine em3d_exchange_bndb_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
integer(ISZ) :: ibuf,iz
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 9
            call mpi_packbuffer_init((3*yfu%nzguard-1)*size(yfu%by(:,:,0)),ibuf)
            do iz=yfu%izmin,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%bx(:,:,iz),ibuf)
            end do
            do iz=yfu%izmin, yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%by(:,:,iz),ibuf)
            end do
            do iz=yfu%izmin+1,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%bz(:,:,iz),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 10
            call mpi_packbuffer_init((3*yfl%nzguard-1)*size(yfl%by(:,:,0)),ibuf)
            do iz=yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%bx(:,:,iz),ibuf)
            end do
            do iz=yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%by(:,:,iz),ibuf)
            end do
            do iz=yfl%izmax-yfl%nzguard+1,yfl%izmax-1
              call mympi_pack(yfl%bz(:,:,iz),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
          yfl%bx(:,:,yfl%izmax  :yfl%izmaxg-1) = yfu%bx(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)
          yfl%by(:,:,yfl%izmax  :yfl%izmaxg-1) = yfu%by(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)
          yfl%bz(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%bz(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)

          yfu%bx(:,:,yfu%izming  :yfu%izmin-1) = yfl%bx(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)
          yfu%by(:,:,yfu%izming  :yfu%izmin-1) = yfl%by(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)
          yfu%bz(:,:,yfu%izming+1:yfu%izmin-1) = yfl%bz(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
          
#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=fu%proc) return
          yfl%bx(:,:,yfl%izmax  :yfl%izmaxg-1) = (syfu%bxy(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1) &
                                               +  syfu%bxz(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1))/syfu%clight
          yfl%by(:,:,yfl%izmax  :yfl%izmaxg-1) = (syfu%byx(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1) &
                                               +  syfu%byz(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1))/syfu%clight
          yfl%bz(:,:,yfl%izmax+1:yfl%izmaxg-1) = (syfu%bzx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                               +  syfu%bzy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1))/syfu%clight
          syfu%bxy(:,:,syfu%izming  :syfu%izmin-1) = yfl%bx(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)*syfu%clight
          syfu%bxz(:,:,syfu%izming  :syfu%izmin-1) = 0.
          syfu%byx(:,:,syfu%izming  :syfu%izmin-1) = 0.
          syfu%byz(:,:,syfu%izming  :syfu%izmin-1) = yfl%by(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)*syfu%clight
          syfu%bzx(:,:,syfu%izming+1:syfu%izmin-1) = yfl%bz(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)*syfu%clight
          syfu%bzy(:,:,syfu%izming+1:syfu%izmin-1) = 0.
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=fu%proc) return
          yfu%bx(:,:,yfu%izming  :yfu%izmin-1) = (syfl%bxy(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1) &
                                               +  syfl%bxz(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1))/syfl%clight
          yfu%by(:,:,yfu%izming  :yfu%izmin-1) = (syfl%byx(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1) &
                                               +  syfl%byz(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1))/syfl%clight
          yfu%bz(:,:,yfu%izming+1:yfu%izmin-1) = (syfl%bzx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                               +  syfl%bzy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1))/syfl%clight
          syfl%bxy(:,:,syfl%izmax  :syfl%izmaxg-1) = yfu%bx(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)*syfl%clight
          syfl%bxz(:,:,syfl%izmax  :syfl%izmaxg-1) = 0.
          syfl%byx(:,:,syfl%izmax  :syfl%izmaxg-1) = 0.
          syfl%byz(:,:,syfl%izmax  :syfl%izmaxg-1) = yfu%by(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)*syfl%clight
          syfl%bzx(:,:,syfl%izmax+1:syfl%izmaxg-1) = yfu%bz(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)*syfl%clight
          syfl%bzy(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 11
            call mpi_packbuffer_init(4*size(syfu%byx(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)) &
                                    +2*size(syfu%byx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)),ibuf)

            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%bxy(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%bxz(:,:,iz),ibuf)
            end do
            
            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%byx(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%byz(:,:,iz),ibuf)
            end do
            
            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%bzx(:,:,iz),ibuf)
            end do
            do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
              call mympi_pack(syfu%bzy(:,:,iz),ibuf)
            end do

            call mpi_isend_pack(fl%proc,3,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 12
            call mpi_packbuffer_init(4*size(syfl%byx(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)) &
                                    +2*size(syfl%byx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)),ibuf)

            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%bxy(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%bxz(:,:,iz),ibuf)
            end do

            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%byx(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard,syfl%izmax-1
              call mympi_pack(syfl%byz(:,:,iz),ibuf)
            end do

            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%bzx(:,:,iz),ibuf)
            end do
            do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
              call mympi_pack(syfl%bzy(:,:,iz),ibuf)
            end do

            call mpi_isend_pack(fu%proc,4,ibuf)
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%bxy(:,:,syfu%izming  :syfu%izmin-1) = syfl%bxy(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfu%bxz(:,:,syfu%izming  :syfu%izmin-1) = syfl%bxz(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfu%byx(:,:,syfu%izming  :syfu%izmin-1) = syfl%byx(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfu%byz(:,:,syfu%izming  :syfu%izmin-1) = syfl%byz(:,:,syfl%izmax-syfl%nzguard  :syfl%izmax-1)
          syfu%bzx(:,:,syfu%izming+1:syfu%izmin-1) = syfl%bzx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%bzy(:,:,syfu%izming+1:syfu%izmin-1) = syfl%bzy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfl%bxy(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%bxy(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfl%bxz(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%bxz(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfl%byx(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%byx(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfl%byz(:,:,syfl%izmax  :syfl%izmaxg-1) = syfu%byz(:,:,syfu%izmin  :syfu%izmin+syfu%nzguard-1)
          syfl%bzx(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%bzx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%bzy(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%bzy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndb_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndb_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ibuf,iz

          if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 13
            call mpi_packbuffer_init((3*yfu%nzguard-1)*size(yfu%bx(:,:,0)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iz=yfu%izming,yfu%izmin-1
              yfu%bx(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%bx(:,:,iz)),ibuf),shape(yfu%bx(:,:,iz)))
            end do
            do iz=yfu%izming,yfu%izmin-1
              yfu%by(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%by(:,:,iz)),ibuf),shape(yfu%by(:,:,iz)))
            end do
            do iz=yfu%izming+1,yfu%izmin-1
              yfu%bz(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%bz(:,:,iz)),ibuf),shape(yfu%bz(:,:,iz)))
            end do
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 14
            call mpi_packbuffer_init((3*yfl%nzguard-1)*size(yfl%bx(:,:,0)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iz=yfl%izmax,yfl%izmaxg-1
              yfl%bx(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%bx(:,:,iz)),ibuf),shape(yfl%bx(:,:,iz)))
            end do
            do iz=yfl%izmax,yfl%izmaxg-1
              yfl%by(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%by(:,:,iz)),ibuf),shape(yfl%by(:,:,iz)))
            end do
            do iz=yfl%izmax+1,yfl%izmaxg-1
              yfl%bz(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%bz(:,:,iz)),ibuf),shape(yfl%bz(:,:,iz)))
            end do
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ibuf = 15
            call mpi_packbuffer_init(4*size(syfu%bxy(:,:,syfu%izming  :syfu%izmin-1)) &
                                    +2*size(syfu%bxy(:,:,syfu%izming+1:syfu%izmin-1)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)
            ! --- recv data from down in z
            do iz=syfu%izming,syfu%izmin-1
              syfu%bxy(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfu%bxy(:,:,iz)),ibuf),shape(syfu%bxy(:,:,iz)))
            end do
            do iz=syfu%izming,syfu%izmin-1
              syfu%bxz(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfu%bxz(:,:,iz)),ibuf),shape(syfu%bxz(:,:,iz)))
            end do
            do iz=syfu%izming,syfu%izmin-1
              syfu%byx(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfu%byx(:,:,iz)),ibuf),shape(syfu%byx(:,:,iz)))
            end do
            do iz=syfu%izming,syfu%izmin-1
              syfu%byz(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfu%byz(:,:,iz)),ibuf),shape(syfu%byz(:,:,iz)))
            end do
            do iz=syfu%izming+1,syfu%izmin-1
              syfu%bzx(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfu%bzx(:,:,iz)),ibuf),shape(syfu%bzx(:,:,iz)))
            end do
            do iz=syfu%izming+1,syfu%izmin-1
              syfu%bzy(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfu%bzy(:,:,iz)),ibuf),shape(syfu%bzy(:,:,iz)))
            end do
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 16
            call mpi_packbuffer_init(4*size(syfl%bxy(:,:,syfl%izmax  :syfl%izmaxg-1)) &
                                    +2*size(syfl%bxy(:,:,syfl%izmax+1:syfl%izmaxg-1)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%bxy(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfl%bxy(:,:,iz)),ibuf),shape(syfl%bxy(:,:,iz)))
            end do
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%bxz(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfl%bxz(:,:,iz)),ibuf),shape(syfl%bxz(:,:,iz)))
            end do
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%byx(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfl%byx(:,:,iz)),ibuf),shape(syfl%byx(:,:,iz)))
            end do
            do iz=syfl%izmax,syfl%izmaxg-1
              syfl%byz(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfl%byz(:,:,iz)),ibuf),shape(syfl%byz(:,:,iz)))
            end do
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%bzx(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfl%bzx(:,:,iz)),ibuf),shape(syfl%bzx(:,:,iz)))
            end do
            do iz=syfl%izmax+1,syfl%izmaxg-1
              syfl%bzy(:,:,iz) =  reshape(mpi_unpack_real_array( size(syfl%bzy(:,:,iz)),ibuf),shape(syfl%bzy(:,:,iz)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndb_zrecv
#endif

subroutine em3d_exchange_bndf_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
integer(ISZ)   ::iz,ibuf
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          ! --- case lower yee, upper yee
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 1
            if (yfu%nzguard>1) then
              call mpi_packbuffer_init((yfu%nzguard-1)*size(yfu%f(:,:,0)),ibuf)
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                call mympi_pack(yfu%f(:,:,iz),ibuf)
              end do
              call mpi_isend_pack(fl%proc,1,ibuf)
            end if
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 2
            if (yfl%nzguard>1) then
              call mpi_packbuffer_init((yfl%nzguard-1)*size(yfl%f(:,:,0)),ibuf)
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                call mympi_pack(yfl%f(:,:,iz),ibuf)
              end do
              call mpi_isend_pack(fu%proc,2,ibuf)
            end if
          else
#endif
            yfl%f(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%f(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
            yfu%f(:,:,yfu%izming+1:yfu%izmin-1) = yfl%f(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          ! --- case lower yee, upper split yee
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%f(:,:,yfl%izmax+1:yfl%izmaxg-1) = syfu%fx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                              + syfu%fy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1) &
                                              + syfu%fz(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfu%fx(:,:,syfu%izming+1:syfu%izmin-1) = 0.
          syfu%fy(:,:,syfu%izming+1:syfu%izmin-1) = 0.
          syfu%fz(:,:,syfu%izming+1:syfu%izmin-1) = yfl%f(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)

      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
          ! --- case lower split yee, upper yee
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%f(:,:,yfu%izming+1:yfu%izmin-1)      = syfl%fx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                                   + syfl%fy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1) &
                                                   + syfl%fz(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfl%fx(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
          syfl%fy(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
          syfl%fz(:,:,syfl%izmax+1:syfl%izmaxg-1) = yfu%f(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
        case(splityeefield)
          ! --- case lower split yee, upper split yee
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            if (syfu%nzguard>1) then
              ibuf = 3
              call mpi_packbuffer_init( 3*int(size(syfu%fx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1))) ,ibuf)
              do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
                call mympi_pack(syfu%fx(:,:,iz),ibuf)
              end do
              do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
                call mympi_pack(syfu%fy(:,:,iz),ibuf)
              end do
              do iz=syfu%izmin+1,syfu%izmin+syfu%nzguard-1
                call mympi_pack(syfu%fz(:,:,iz),ibuf)
              end do
              call mpi_isend_pack(fl%proc,3,ibuf)
              mpireqpnt=mpireqpnt+1
            end if
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            if (syfl%nzguard>1) then
              ibuf = 4
              call mpi_packbuffer_init(3*int(size(syfl%fx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1))) ,ibuf)
              do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
                call mympi_pack(syfl%fx(:,:,iz),ibuf)
              end do
              do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
                call mympi_pack(syfl%fy(:,:,iz),ibuf)
              end do
              do iz=syfl%izmax-syfl%nzguard+1,syfl%izmax-1
                call mympi_pack(syfl%fz(:,:,iz),ibuf)
              end do
              call mpi_isend_pack(fu%proc,4,ibuf)
              mpireqpnt=mpireqpnt+1
            end if
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%fx(:,:,syfu%izming+1:syfu%izmin-1) = syfl%fx(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%fy(:,:,syfu%izming+1:syfu%izmin-1) = syfl%fy(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)
          syfu%fz(:,:,syfu%izming+1:syfu%izmin-1) = syfl%fz(:,:,syfl%izmax-syfl%nzguard+1:syfl%izmax-1)

          syfl%fx(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%fx(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%fy(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%fy(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
          syfl%fz(:,:,syfl%izmax+1:syfl%izmaxg-1) = syfu%fz(:,:,syfu%izmin+1:syfu%izmin+syfu%nzguard-1)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndf_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndf_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            if (yfu%nzguard>1) then
              ibuf = 5
              call mpi_packbuffer_init((yfu%nzguard-1)*size(yfu%Ez(:,:,0)),ibuf)
              call mpi_recv_pack(fl%proc,2,ibuf)
              do iz = yfu%izming+1,yfu%izmin-1
                yfu%f(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%F(:,:,0)),ibuf),shape(yfu%F(:,:,0)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            if (yfl%nzguard>1) then
              ibuf = 6
              call mpi_packbuffer_init((yfl%nzguard-1)*size(yfl%ez(:,:,yfl%izmin)),ibuf)
              call mpi_recv_pack(fu%proc,1,ibuf)
              do iz = yfl%izmax+1,yfl%izmaxg-1
                yfl%f(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%F(:,:,0)),ibuf),shape(yfl%F(:,:,0)))
              end do
            end if
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            if (syfu%nzguard>1) then
              ibuf = 7
              call mpi_packbuffer_init(3*size(syfu%ezx(:,:,syfu%izming+1:syfu%izmin-1)),ibuf)
              call mpi_recv_pack(fl%proc,4,ibuf)

              do iz = syfu%izming+1,syfu%izmin-1
                syfu%fx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%fx(:,:,iz)),ibuf),shape(syfu%fx(:,:,iz)))
              end do
              do iz = syfu%izming+1,syfu%izmin-1
                syfu%fy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%fy(:,:,iz)),ibuf),shape(syfu%fy(:,:,iz)))
              end do
              do iz = syfu%izming+1,syfu%izmin-1
                syfu%fz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfu%fz(:,:,iz)),ibuf),shape(syfu%fz(:,:,iz)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            if (syfl%nzguard>1) then
              ibuf = 8
              call mpi_packbuffer_init(3*size(syfl%ezx(:,:,syfl%izmax+1:syfl%izmaxg-1)),ibuf)
              call mpi_recv_pack(fu%proc,3,ibuf)
              do iz=syfl%izmax+1,syfl%izmaxg-1
                syfl%fx(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%fx(:,:,iz)),ibuf),shape(syfl%fx(:,:,iz)))
              end do
              do iz=syfl%izmax+1,syfl%izmaxg-1
                syfl%fy(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%fy(:,:,iz)),ibuf),shape(syfl%fy(:,:,iz)))
              end do
              do iz=syfl%izmax+1,syfl%izmaxg-1
                syfl%fz(:,:,iz) = reshape(mpi_unpack_real_array( size(syfl%fz(:,:,iz)),ibuf),shape(syfl%fz(:,:,iz)))
              end do
            end if
          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bndf_zrecv
#endif

subroutine em3d_exchange_bndj_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz,ibuf

#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then

            ! --- send data down in z
            ibuf = 17
            call mpi_packbuffer_init((3*yfu%nzguard+2)*size(yfu%J(:,:,-1,1)),ibuf)
            do iz = yfu%izming,yfu%izmin
              call mympi_pack(yfu%J(:,:,iz,1),ibuf)
              call mympi_pack(yfu%J(:,:,iz,2),ibuf)
            end do
            do iz = yfu%izming,yfu%izmin-1
              call mympi_pack(yfu%J(:,:,iz,3),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            ibuf = 18
            call mpi_packbuffer_init((3*yfl%nzguard+2)*size(yfl%J(:,:,0,1)),ibuf)
            do iz = yfl%izmax, yfl%izmaxg
              call mympi_pack(yfl%J(:,:,iz,1),ibuf)
              call mympi_pack(yfl%J(:,:,iz,2),ibuf)
            end do
            do iz = yfl%izmax, yfl%izmaxg-1
              call mympi_pack(yfl%J(:,:,iz,3),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)

          else
#endif
          iz = yfu%nzguard
          yfu%J(:,:,0:iz,1:2)    = yfu%J(:,:,0:iz,1:2)    + yfl%J(:,:,yfl%nz:yfl%nz+iz,1:2)
          yfu%J(:,:,0:iz-1,3)    = yfu%J(:,:,0:iz-1,3)    + yfl%J(:,:,yfl%nz:yfl%nz+iz-1,3) 

          yfl%J(:,:,yfl%nz-iz:yfl%nz-1,:) = yfl%J(:,:,yfl%nz-iz:yfl%nz-1,:) + yfu%J(:,:,-iz:-1,:)
          yfl%J(:,:,yfl%nz,1:2) = yfu%J(:,:,0,1:2)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndj_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz,ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            ibuf = 19
            call mpi_packbuffer_init((3*yfu%nzguard+2)*size(yfu%J(:,:,0,1)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iz = 0,yfu%nzguard
              yfu%J(:,:,iz  ,1) = yfu%J(:,:,iz  ,1) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),ibuf), &
                                                                                     shape(yfu%J(:,:,0,1)))
              yfu%J(:,:,iz  ,2) = yfu%J(:,:,iz  ,2) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),ibuf), &
                                                                                     shape(yfu%J(:,:,0,1)))
            end do
            do iz = 0,yfu%nzguard-1
              yfu%J(:,:,iz  ,3) = yfu%J(:,:,iz  ,3) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),ibuf), &
                                                                                     shape(yfu%J(:,:,0,1)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            ibuf = 20
            call mpi_packbuffer_init((3*yfl%nzguard+2)*size(yfl%J(:,:,0,1)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iz = -yfl%nzguard,0
              yfl%J(:,:,yfl%nz+iz,1) = yfl%J(:,:,yfl%nz+iz,1) + reshape(mpi_unpack_real_array( size(yfl%J(:,:,yfl%nz-1,1)),ibuf),&
                                                                                            shape(yfl%J(:,:,yfl%nz-1,1)))
              yfl%J(:,:,yfl%nz+iz,2) = yfl%J(:,:,yfl%nz+iz,2) + reshape(mpi_unpack_real_array( size(yfl%J(:,:,yfl%nz-1,2)),ibuf),&
                                                                                            shape(yfl%J(:,:,yfl%nz-1,2)))
            end do
            do iz = -yfl%nzguard,-1
              yfl%J(:,:,yfl%nz+iz,3) = yfl%J(:,:,yfl%nz+iz,3) + reshape(mpi_unpack_real_array( size(yfl%J(:,:,yfl%nz-1,3)),ibuf),&
                                                                                            shape(yfl%J(:,:,yfl%nz-1,3)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndj_zrecv
#endif

subroutine em3d_exchange_bndrho_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz,ibuf

#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ibuf=21
            ! --- send data down in z
            call mpi_packbuffer_init(size(yfu%rho(:,:,0:yfu%nzguard)),ibuf)
            do iz = -yfu%nzguard,0
              call mympi_pack(yfu%rho(:,:,iz  ),ibuf)
            end do
            call mpi_isend_pack(fl%proc,1,ibuf)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            ibuf = 22
            call mpi_packbuffer_init(size(yfl%rho(:,:,0:yfl%nzguard)),ibuf)
            do iz = 0,yfl%nzguard
              call mympi_pack(yfl%rho(:,:,yfl%nz+iz  ),ibuf)
            end do
            call mpi_isend_pack(fu%proc,2,ibuf)

          else
#endif
          iz = yfu%nzguard
          yfu%Rho(:,:,0:iz)      = yfu%Rho(:,:,0:iz)      + yfl%Rho(:,:,yfl%nz:yfl%nz+iz)
          yfl%Rho(:,:,yfl%nz-iz:yfl%nz-1) = yfl%Rho(:,:,yfl%nz-iz:yfl%nz-1) + yfu%Rho(:,:,-iz:-1)
          yfl%Rho(:,:,yfl%nz)   = yfu%Rho(:,:,0)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndrho_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz, ibuf

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            ibuf = 23
            call mpi_packbuffer_init(size(yfu%rho(:,:,0:yfu%nzguard)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            do iz = 0,yfu%nzguard
              yfu%rho(:,:,iz  ) = yfu%rho(:,:,iz  ) + reshape(mpi_unpack_real_array( size(yfu%rho(:,:,0)),ibuf), &
                                                                                    shape(yfu%rho(:,:,0)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            ibuf = 24
            call mpi_packbuffer_init(size(yfl%rho(:,:,0:yfl%nzguard)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            do iz = -yfl%nzguard,0
              yfl%rho(:,:,yfl%nz+iz  ) = yfl%rho(:,:,yfl%nz+iz  ) + reshape(mpi_unpack_real_array(size(yfl%rho(:,:,yfl%nz)),ibuf),&
                                                                                            shape(yfl%rho(:,:,yfl%nz  )))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndrho_zrecv
#endif

subroutine em3d_exchange_e(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
 
  ! --- X
  ! core<--->sides
  call em3d_exchange_bnde_x(b%core,   b%sidexr)
  call em3d_exchange_bnde_x(b%sidexl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%sidexl, b%core)
  call em3d_exchange_bnde_xrecv(b%core,   b%sidexr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bnde_x(b%sideyl,   b%edgexryl)
  call em3d_exchange_bnde_x(b%edgexlyl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%edgexlyl, b%sideyl)
  call em3d_exchange_bnde_xrecv(b%sideyl,   b%edgexryl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_x(b%sideyr,   b%edgexryr)
  call em3d_exchange_bnde_x(b%edgexlyr, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%edgexlyr, b%sideyr)
  call em3d_exchange_bnde_xrecv(b%sideyr,   b%edgexryr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_x(b%sidezl,   b%edgexrzl)
  call em3d_exchange_bnde_x(b%edgexlzl, b%sidezl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%edgexlzl, b%sidezl)
  call em3d_exchange_bnde_xrecv(b%sidezl,   b%edgexrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_x(b%sidezr,   b%edgexrzr)
  call em3d_exchange_bnde_x(b%edgexlzr, b%sidezr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%edgexlzr, b%sidezr)
  call em3d_exchange_bnde_xrecv(b%sidezr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bnde_x(b%edgeylzl,     b%cornerxrylzl)
  call em3d_exchange_bnde_x(b%cornerxlylzl, b%edgeylzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%cornerxlylzl, b%edgeylzl)
  call em3d_exchange_bnde_xrecv(b%edgeylzl,     b%cornerxrylzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_x(b%edgeyrzl,     b%cornerxryrzl)
  call em3d_exchange_bnde_x(b%cornerxlyrzl, b%edgeyrzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%cornerxlyrzl, b%edgeyrzl)
  call em3d_exchange_bnde_xrecv(b%edgeyrzl,     b%cornerxryrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_x(b%edgeylzr,     b%cornerxrylzr)
  call em3d_exchange_bnde_x(b%cornerxlylzr, b%edgeylzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%cornerxlylzr, b%edgeylzr)
  call em3d_exchange_bnde_xrecv(b%edgeylzr,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_x(b%edgeyrzr,     b%cornerxryrzr)
  call em3d_exchange_bnde_x(b%cornerxlyrzr, b%edgeyrzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_xrecv(b%cornerxlyrzr, b%edgeyrzr)
  call em3d_exchange_bnde_xrecv(b%edgeyrzr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bnde_y(b%core,   b%sideyr)
  call em3d_exchange_bnde_y(b%sideyl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%sideyl, b%core)
  call em3d_exchange_bnde_yrecv(b%core,   b%sideyr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bnde_y(b%sidexl,   b%edgexlyr)
  call em3d_exchange_bnde_y(b%edgexlyl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%edgexlyl, b%sidexl)
  call em3d_exchange_bnde_yrecv(b%sidexl,   b%edgexlyr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_y(b%sidexr,   b%edgexryr)
  call em3d_exchange_bnde_y(b%edgexryl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%edgexryl, b%sidexr)
  call em3d_exchange_bnde_yrecv(b%sidexr,   b%edgexryr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_y(b%sidezl,   b%edgeyrzl)
  call em3d_exchange_bnde_y(b%edgeylzl, b%sidezl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%edgeylzl, b%sidezl)
  call em3d_exchange_bnde_yrecv(b%sidezl,   b%edgeyrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_y(b%sidezr,   b%edgeyrzr)
  call em3d_exchange_bnde_y(b%edgeylzr, b%sidezr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%edgeylzr, b%sidezr)
  call em3d_exchange_bnde_yrecv(b%sidezr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bnde_y(b%edgexlzl,     b%cornerxlyrzl)
  call em3d_exchange_bnde_y(b%cornerxlylzl, b%edgexlzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%cornerxlylzl, b%edgexlzl)
  call em3d_exchange_bnde_yrecv(b%edgexlzl,     b%cornerxlyrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_y(b%edgexrzl,     b%cornerxryrzl)
  call em3d_exchange_bnde_y(b%cornerxrylzl, b%edgexrzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%cornerxrylzl, b%edgexrzl)
  call em3d_exchange_bnde_yrecv(b%edgexrzl,     b%cornerxryrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_y(b%edgexlzr,     b%cornerxlyrzr)
  call em3d_exchange_bnde_y(b%cornerxlylzr, b%edgexlzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%cornerxlylzr, b%edgexlzr)
  call em3d_exchange_bnde_yrecv(b%edgexlzr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_y(b%edgexrzr,     b%cornerxryrzr)
  call em3d_exchange_bnde_y(b%cornerxrylzr, b%edgexrzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_yrecv(b%cornerxrylzr, b%edgexrzr)
  call em3d_exchange_bnde_yrecv(b%edgexrzr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bnde_z(b%core,   b%sidezr)
  call em3d_exchange_bnde_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%sidezl, b%core)
  call em3d_exchange_bnde_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bnde_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bnde_z(b%edgexlzl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgexlzl, b%sidexl)
  call em3d_exchange_bnde_zrecv(b%sidexl,   b%edgexlzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bnde_z(b%edgexrzl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgexrzl, b%sidexr)
  call em3d_exchange_bnde_zrecv(b%sidexr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bnde_z(b%edgeylzl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgeylzl, b%sideyl)
  call em3d_exchange_bnde_zrecv(b%sideyl,   b%edgeylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bnde_z(b%edgeyrzl, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgeyrzl, b%sideyr)
  call em3d_exchange_bnde_zrecv(b%sideyr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bnde_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bnde_z(b%cornerxlylzl, b%edgexlyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bnde_zrecv(b%edgexlyl,     b%cornerxlylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bnde_z(b%cornerxrylzl, b%edgexryl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bnde_zrecv(b%edgexryl,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bnde_z(b%cornerxlyrzl, b%edgexlyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bnde_zrecv(b%edgexlyr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bnde_z(b%cornerxryrzl, b%edgexryr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxryrzl, b%edgexryr)
  call em3d_exchange_bnde_zrecv(b%edgexryr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_e

subroutine em3d_exchange_b(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
 
  ! --- X
  ! core<--->sides
  call em3d_exchange_bndb_x(b%core,   b%sidexr)
  call em3d_exchange_bndb_x(b%sidexl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%sidexl, b%core)
  call em3d_exchange_bndb_xrecv(b%core,   b%sidexr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndb_x(b%sideyl,   b%edgexryl)
  call em3d_exchange_bndb_x(b%edgexlyl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%edgexlyl, b%sideyl)
  call em3d_exchange_bndb_xrecv(b%sideyl,   b%edgexryl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_x(b%sideyr,   b%edgexryr)
  call em3d_exchange_bndb_x(b%edgexlyr, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%edgexlyr, b%sideyr)
  call em3d_exchange_bndb_xrecv(b%sideyr,   b%edgexryr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_x(b%sidezl,   b%edgexrzl)
  call em3d_exchange_bndb_x(b%edgexlzl, b%sidezl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%edgexlzl, b%sidezl)
  call em3d_exchange_bndb_xrecv(b%sidezl,   b%edgexrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_x(b%sidezr,   b%edgexrzr)
  call em3d_exchange_bndb_x(b%edgexlzr, b%sidezr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%edgexlzr, b%sidezr)
  call em3d_exchange_bndb_xrecv(b%sidezr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndb_x(b%edgeylzl,     b%cornerxrylzl)
  call em3d_exchange_bndb_x(b%cornerxlylzl, b%edgeylzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%cornerxlylzl, b%edgeylzl)
  call em3d_exchange_bndb_xrecv(b%edgeylzl,     b%cornerxrylzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_x(b%edgeyrzl,     b%cornerxryrzl)
  call em3d_exchange_bndb_x(b%cornerxlyrzl, b%edgeyrzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%cornerxlyrzl, b%edgeyrzl)
  call em3d_exchange_bndb_xrecv(b%edgeyrzl,     b%cornerxryrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_x(b%edgeylzr,     b%cornerxrylzr)
  call em3d_exchange_bndb_x(b%cornerxlylzr, b%edgeylzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%cornerxlylzr, b%edgeylzr)
  call em3d_exchange_bndb_xrecv(b%edgeylzr,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_x(b%edgeyrzr,     b%cornerxryrzr)
  call em3d_exchange_bndb_x(b%cornerxlyrzr, b%edgeyrzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_xrecv(b%cornerxlyrzr, b%edgeyrzr)
  call em3d_exchange_bndb_xrecv(b%edgeyrzr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndb_y(b%core,   b%sideyr)
  call em3d_exchange_bndb_y(b%sideyl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%sideyl, b%core)
  call em3d_exchange_bndb_yrecv(b%core,   b%sideyr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndb_y(b%sidexl,   b%edgexlyr)
  call em3d_exchange_bndb_y(b%edgexlyl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%edgexlyl, b%sidexl)
  call em3d_exchange_bndb_yrecv(b%sidexl,   b%edgexlyr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_y(b%sidexr,   b%edgexryr)
  call em3d_exchange_bndb_y(b%edgexryl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%edgexryl, b%sidexr)
  call em3d_exchange_bndb_yrecv(b%sidexr,   b%edgexryr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_y(b%sidezl,   b%edgeyrzl)
  call em3d_exchange_bndb_y(b%edgeylzl, b%sidezl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%edgeylzl, b%sidezl)
  call em3d_exchange_bndb_yrecv(b%sidezl,   b%edgeyrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_y(b%sidezr,   b%edgeyrzr)
  call em3d_exchange_bndb_y(b%edgeylzr, b%sidezr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%edgeylzr, b%sidezr)
  call em3d_exchange_bndb_yrecv(b%sidezr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndb_y(b%edgexlzl,     b%cornerxlyrzl)
  call em3d_exchange_bndb_y(b%cornerxlylzl, b%edgexlzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%cornerxlylzl, b%edgexlzl)
  call em3d_exchange_bndb_yrecv(b%edgexlzl,     b%cornerxlyrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_y(b%edgexrzl,     b%cornerxryrzl)
  call em3d_exchange_bndb_y(b%cornerxrylzl, b%edgexrzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%cornerxrylzl, b%edgexrzl)
  call em3d_exchange_bndb_yrecv(b%edgexrzl,     b%cornerxryrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_y(b%edgexlzr,     b%cornerxlyrzr)
  call em3d_exchange_bndb_y(b%cornerxlylzr, b%edgexlzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%cornerxlylzr, b%edgexlzr)
  call em3d_exchange_bndb_yrecv(b%edgexlzr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_y(b%edgexrzr,     b%cornerxryrzr)
  call em3d_exchange_bndb_y(b%cornerxrylzr, b%edgexrzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_yrecv(b%cornerxrylzr, b%edgexrzr)
  call em3d_exchange_bndb_yrecv(b%edgexrzr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndb_z(b%core,   b%sidezr)
  call em3d_exchange_bndb_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndb_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndb_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bndb_z(b%edgexlzl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgexlzl, b%sidexl)
  call em3d_exchange_bndb_zrecv(b%sidexl,   b%edgexlzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bndb_z(b%edgexrzl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgexrzl, b%sidexr)
  call em3d_exchange_bndb_zrecv(b%sidexr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bndb_z(b%edgeylzl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgeylzl, b%sideyl)
  call em3d_exchange_bndb_zrecv(b%sideyl,   b%edgeylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bndb_z(b%edgeyrzl, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgeyrzl, b%sideyr)
  call em3d_exchange_bndb_zrecv(b%sideyr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndb_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bndb_z(b%cornerxlylzl, b%edgexlyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bndb_zrecv(b%edgexlyl,     b%cornerxlylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bndb_z(b%cornerxrylzl, b%edgexryl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bndb_zrecv(b%edgexryl,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bndb_z(b%cornerxlyrzl, b%edgexlyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bndb_zrecv(b%edgexlyr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bndb_z(b%cornerxryrzl, b%edgexryr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxryrzl, b%edgexryr)
  call em3d_exchange_bndb_zrecv(b%edgexryr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_b

subroutine em3d_exchange_f(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
 
  ! --- X
  ! core<--->sides
  call em3d_exchange_bndf_x(b%core,   b%sidexr)
  call em3d_exchange_bndf_x(b%sidexl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%sidexl, b%core)
  call em3d_exchange_bndf_xrecv(b%core,   b%sidexr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndf_x(b%sideyl,   b%edgexryl)
  call em3d_exchange_bndf_x(b%edgexlyl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%edgexlyl, b%sideyl)
  call em3d_exchange_bndf_xrecv(b%sideyl,   b%edgexryl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_x(b%sideyr,   b%edgexryr)
  call em3d_exchange_bndf_x(b%edgexlyr, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%edgexlyr, b%sideyr)
  call em3d_exchange_bndf_xrecv(b%sideyr,   b%edgexryr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_x(b%sidezl,   b%edgexrzl)
  call em3d_exchange_bndf_x(b%edgexlzl, b%sidezl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%edgexlzl, b%sidezl)
  call em3d_exchange_bndf_xrecv(b%sidezl,   b%edgexrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_x(b%sidezr,   b%edgexrzr)
  call em3d_exchange_bndf_x(b%edgexlzr, b%sidezr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%edgexlzr, b%sidezr)
  call em3d_exchange_bndf_xrecv(b%sidezr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndf_x(b%edgeylzl,     b%cornerxrylzl)
  call em3d_exchange_bndf_x(b%cornerxlylzl, b%edgeylzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%cornerxlylzl, b%edgeylzl)
  call em3d_exchange_bndf_xrecv(b%edgeylzl,     b%cornerxrylzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_x(b%edgeyrzl,     b%cornerxryrzl)
  call em3d_exchange_bndf_x(b%cornerxlyrzl, b%edgeyrzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%cornerxlyrzl, b%edgeyrzl)
  call em3d_exchange_bndf_xrecv(b%edgeyrzl,     b%cornerxryrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_x(b%edgeylzr,     b%cornerxrylzr)
  call em3d_exchange_bndf_x(b%cornerxlylzr, b%edgeylzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%cornerxlylzr, b%edgeylzr)
  call em3d_exchange_bndf_xrecv(b%edgeylzr,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_x(b%edgeyrzr,     b%cornerxryrzr)
  call em3d_exchange_bndf_x(b%cornerxlyrzr, b%edgeyrzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_xrecv(b%cornerxlyrzr, b%edgeyrzr)
  call em3d_exchange_bndf_xrecv(b%edgeyrzr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndf_y(b%core,   b%sideyr)
  call em3d_exchange_bndf_y(b%sideyl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%sideyl, b%core)
  call em3d_exchange_bndf_yrecv(b%core,   b%sideyr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndf_y(b%sidexl,   b%edgexlyr)
  call em3d_exchange_bndf_y(b%edgexlyl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%edgexlyl, b%sidexl)
  call em3d_exchange_bndf_yrecv(b%sidexl,   b%edgexlyr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_y(b%sidexr,   b%edgexryr)
  call em3d_exchange_bndf_y(b%edgexryl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%edgexryl, b%sidexr)
  call em3d_exchange_bndf_yrecv(b%sidexr,   b%edgexryr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_y(b%sidezl,   b%edgeyrzl)
  call em3d_exchange_bndf_y(b%edgeylzl, b%sidezl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%edgeylzl, b%sidezl)
  call em3d_exchange_bndf_yrecv(b%sidezl,   b%edgeyrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_y(b%sidezr,   b%edgeyrzr)
  call em3d_exchange_bndf_y(b%edgeylzr, b%sidezr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%edgeylzr, b%sidezr)
  call em3d_exchange_bndf_yrecv(b%sidezr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndf_y(b%edgexlzl,     b%cornerxlyrzl)
  call em3d_exchange_bndf_y(b%cornerxlylzl, b%edgexlzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%cornerxlylzl, b%edgexlzl)
  call em3d_exchange_bndf_yrecv(b%edgexlzl,     b%cornerxlyrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_y(b%edgexrzl,     b%cornerxryrzl)
  call em3d_exchange_bndf_y(b%cornerxrylzl, b%edgexrzl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%cornerxrylzl, b%edgexrzl)
  call em3d_exchange_bndf_yrecv(b%edgexrzl,     b%cornerxryrzl)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_y(b%edgexlzr,     b%cornerxlyrzr)
  call em3d_exchange_bndf_y(b%cornerxlylzr, b%edgexlzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%cornerxlylzr, b%edgexlzr)
  call em3d_exchange_bndf_yrecv(b%edgexlzr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_y(b%edgexrzr,     b%cornerxryrzr)
  call em3d_exchange_bndf_y(b%cornerxrylzr, b%edgexrzr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_yrecv(b%cornerxrylzr, b%edgexrzr)
  call em3d_exchange_bndf_yrecv(b%edgexrzr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndf_z(b%core,   b%sidezr)
  call em3d_exchange_bndf_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndf_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndf_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bndf_z(b%edgexlzl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgexlzl, b%sidexl)
  call em3d_exchange_bndf_zrecv(b%sidexl,   b%edgexlzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bndf_z(b%edgexrzl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgexrzl, b%sidexr)
  call em3d_exchange_bndf_zrecv(b%sidexr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bndf_z(b%edgeylzl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgeylzl, b%sideyl)
  call em3d_exchange_bndf_zrecv(b%sideyl,   b%edgeylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bndf_z(b%edgeyrzl, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgeyrzl, b%sideyr)
  call em3d_exchange_bndf_zrecv(b%sideyr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndf_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bndf_z(b%cornerxlylzl, b%edgexlyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bndf_zrecv(b%edgexlyl,     b%cornerxlylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bndf_z(b%cornerxrylzl, b%edgexryl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bndf_zrecv(b%edgexryl,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bndf_z(b%cornerxlyrzl, b%edgexlyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bndf_zrecv(b%edgexlyr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bndf_z(b%cornerxryrzl, b%edgexryr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxryrzl, b%edgexryr)
  call em3d_exchange_bndf_zrecv(b%edgexryr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_f

subroutine em3d_exchange_j(b)
use mod_emfield3d
implicit none
TYPE(EM3D_BLOCKtype) :: b
#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2)
#endif

  ! --- X
  ! core<--->sides
  call em3d_exchange_bndj_x(b%core,   b%sidexr)
  if(b%xrbnd /= periodic) call em3d_exchange_bndj_x(b%sidexl, b%core)
#ifdef MPIPARALLEL
  if(b%xrbnd /= periodic) call em3d_exchange_bndj_xrecv(b%sidexl, b%core)
  call em3d_exchange_bndj_xrecv(b%core,   b%sidexr)
  call mpi_waitall_requests()
#endif

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndj_y(b%core,   b%sideyr)
  if(b%yrbnd /= periodic) call em3d_exchange_bndj_y(b%sideyl, b%core)
#ifdef MPIPARALLEL
  if(b%yrbnd /= periodic) call em3d_exchange_bndj_yrecv(b%sideyl, b%core)
  call em3d_exchange_bndj_yrecv(b%core,   b%sideyr)
  call mpi_waitall_requests()
#endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndj_z(b%core,   b%sidezr)
  if(b%zrbnd /= periodic) call em3d_exchange_bndj_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  if(b%zrbnd /= periodic) call em3d_exchange_bndj_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndj_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_j

subroutine em3d_exchange_rho(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b

  ! --- X
  ! core<--->sides
  call em3d_exchange_bndrho_x(b%core,   b%sidexr)
  if(b%xrbnd /= periodic) call em3d_exchange_bndrho_x(b%sidexl, b%core)
#ifdef MPIPARALLEL
  if(b%xrbnd /= periodic) call em3d_exchange_bndrho_xrecv(b%sidexl, b%core)
  call em3d_exchange_bndrho_xrecv(b%core,   b%sidexr)
  call mpi_waitall_requests()
#endif

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndrho_y(b%core,   b%sideyr)
  if(b%yrbnd /= periodic) call em3d_exchange_bndrho_y(b%sideyl, b%core)
#ifdef MPIPARALLEL
  if(b%yrbnd /= periodic) call em3d_exchange_bndrho_yrecv(b%sideyl, b%core)
  call em3d_exchange_bndrho_yrecv(b%core,   b%sideyr)
  call mpi_waitall_requests()
#endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndrho_z(b%core,   b%sidezr)
  if(b%zrbnd /= periodic) call em3d_exchange_bndrho_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  if(b%zrbnd /= periodic) call em3d_exchange_bndrho_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndrho_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif
  return
end subroutine em3d_exchange_rho

subroutine yee2node3d(f)
! puts EM value from Yee grid to nodes
use mod_emfield3d
implicit none
TYPE(EM3D_YEEFIELDtype) :: f

INTEGER :: j,k,l
!return
  do l=-f%nzguard,f%nz+f%nzguard
    do k=-f%nyguard,f%ny+f%nyguard
      do j=f%nx+f%nxguard-1,-f%nxguard+1,-1
        f%ex(j,k,l)=0.5*(f%ex(j,k,l)+f%ex(j-1,k,l))
        f%by(j,k,l)=0.5*(f%by(j,k,l)+f%by(j-1,k,l))
        f%bz(j,k,l)=0.5*(f%bz(j,k,l)+f%bz(j-1,k,l))
      enddo
    enddo
  enddo

  do l=-f%nzguard,f%nz+f%nzguard
    do k=f%ny+f%nyguard-1,-f%nyguard+1,-1
      do j=-f%nxguard,f%nx+f%nxguard
        f%ey(j,k,l)=0.5*(f%ey(j,k,l)+f%ey(j,k-1,l))
        f%bz(j,k,l)=0.5*(f%bz(j,k,l)+f%bz(j,k-1,l))
        f%bx(j,k,l)=0.5*(f%bx(j,k,l)+f%bx(j,k-1,l))
      enddo
    enddo
  enddo

  do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
    do k=-f%nyguard,f%ny+f%nyguard
      do j=-f%nxguard,f%nx+f%nxguard
        f%ez(j,k,l)=0.5*(f%ez(j,k,l)+f%ez(j,k,l-1))
        f%bx(j,k,l)=0.5*(f%bx(j,k,l)+f%bx(j,k,l-1))
        f%by(j,k,l)=0.5*(f%by(j,k,l)+f%by(j,k,l-1))
      enddo
    enddo
  enddo

  return
end subroutine yee2node3d

subroutine node2yee3d(f)
! puts EM field back from node to Yee grid
use mod_emfield3d
implicit none
TYPE(EM3D_YEEFIELDtype) :: f

INTEGER :: j,k,l
!return

  do l=-f%nzguard+1,f%nz+f%nzguard-1
    do k=-f%nyguard,f%ny+f%nyguard
      do j=-f%nxguard,f%nx+f%nxguard
        f%ez(j,k,l)=2.*f%ez(j,k,l)-f%ez(j,k,l-1)
        f%bx(j,k,l)=2.*f%bx(j,k,l)-f%bx(j,k,l-1)
        f%by(j,k,l)=2.*f%by(j,k,l)-f%by(j,k,l-1)
      enddo
    enddo
  enddo

  do l=-f%nzguard,f%nz+f%nzguard
    do k=-f%nyguard+1,f%ny+f%nyguard-1
      do j=-f%nxguard,f%nx+f%nxguard
        f%ey(j,k,l)=2.*f%ey(j,k,l)-f%ey(j,k-1,l)
        f%bz(j,k,l)=2.*f%bz(j,k,l)-f%bz(j,k-1,l)
        f%bx(j,k,l)=2.*f%bx(j,k,l)-f%bx(j,k-1,l)
      enddo
    enddo
  enddo

  do l=-f%nzguard,f%nz+f%nzguard
    do k=-f%nyguard,f%ny+f%nyguard
      do j=-f%nxguard+1,f%nx+f%nxguard-1
        f%ex(j,k,l)=2.*f%ex(j,k,l)-f%ex(j-1,k,l)
        f%by(j,k,l)=2.*f%by(j,k,l)-f%by(j-1,k,l)
        f%bz(j,k,l)=2.*f%bz(j,k,l)-f%bz(j-1,k,l)
      enddo
    enddo
  enddo

  return
end subroutine node2yee3d

subroutine add_current_slice_3d(f,i)
use mod_emfield3d
TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: i
  
  f%Jarray(:,:,:,:,i) = f%Jarray(:,:,:,:,i) + f%Jarray(:,:,:,:,i+1)

end subroutine add_current_slice_3d

subroutine add_rho_slice_3d(f,i)
use mod_emfield3d
TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: i
  
  f%Rhoarray(:,:,:,i) = f%Rhoarray(:,:,:,i) + f%Rhoarray(:,:,:,i+1)

end subroutine add_rho_slice_3d

subroutine set_incond(f,n,indx)
use mod_emfield3d
TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: i,n,indx(3,n)
  do i=1,n
    f%incond(indx(1,i),indx(2,i),indx(3,i)) = .true.
  end do  

end subroutine set_incond
