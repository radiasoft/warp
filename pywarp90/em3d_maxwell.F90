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

INTEGER(ISZ), parameter :: pml_user = 0, &
                      pml = 1, &
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
    integer(ISZ) :: j,which,jmin
    real(kind=8) :: dt
    
    if (bnd_cond==pml_user) return
    
    jmin = 1 ! number of cells after which to start absorbing 
    
    sf%afx = 1.
    sf%agx = 1.
    sf%afy = 1.
    sf%agy = 1.
    sf%afz = 1.
    sf%agz = 1.
    sf%sfx = 1.
    sf%sgx = 1.
    sf%sfy = 1.
    sf%sgy = 1.
    sf%sfz = 1.
    sf%sgz = 1.
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
    
    if (bnd_cond==0) return
    
    if (sf%lsx/=0) then
      sigmax=0.
      sigmax_next=0.
      lsigmax=0.
      lsigmax_next=0.
      do j = 0, sf%nx+sf%nxguard
        if (j>=jmin) then
          sigmax(j)      = (sf%smaxx/sf%dx)*(REAL(j-jmin,8)/sf%sdeltax)**sf%nnx
          if (sf%l_nodalgrid) then
            sigmax_next(j) = (sf%smaxx/sf%dx)*((REAL(j-jmin,8)+1.)/sf%sdeltax)**sf%nnx
          else
            sigmax_next(j) = (sf%smaxx/sf%dx)*((REAL(j-jmin,8)+0.5)/sf%sdeltax)**sf%nnx
          end if
        end if
        if (j<sf%nx+sf%nxguard) then
          lsigmax(sf%nx-j) = sigmax(j)
          lsigmax_next(sf%nx-j-1) = sigmax_next(j)
        end if
      end do

      select case(sf%lsx)
        case(1)
          do j = 0, sf%nx+sf%nxguard-1
           if (sf%l_nodalgrid) then
            call assign_coefs(bnd_cond,sf%afx(j),sf%bpfx(j),sf%bmfx(j),sf%sfx(j),sf%clight*dt,sf%dx, &
                              sigmax(j),sigmax_next(j+1), &
                              sb_coef,which,sf%pml_method)
           else
            call assign_coefs(bnd_cond,sf%afx(j),sf%bpfx(j),sf%bmfx(j),sf%sfx(j),sf%clight*dt,sf%dx, &
                              sigmax(j),sigmax_next(j), &
                              sb_coef,which,sf%pml_method)
           end if
          call assign_coefs(bnd_cond,sf%agx(j),sf%bpgx(j),sf%bmgx(j),sf%sgx(j),sf%clight*dt,sf%dx, &
                              sigmax_next(j),sigmax(j+1), &
                              sb_coef,which,sf%pml_method)
          end do
        case(-1)
          do j = sf%nx, -sf%nxguard+1, -1
            call assign_coefs(bnd_cond,sf%afx(j),sf%bmfx(j),sf%bpfx(j),sf%sfx(j),sf%clight*dt,sf%dx, &
                              lsigmax(j),lsigmax_next(j-1), &
                              sb_coef,which,sf%pml_method)
           if (sf%l_nodalgrid) then
            call assign_coefs(bnd_cond,sf%agx(j),sf%bmgx(j),sf%bpgx(j),sf%sgx(j),sf%clight*dt,sf%dx, &
                              lsigmax_next(j-1),lsigmax(j), &
                              sb_coef,which,sf%pml_method)
           else
            call assign_coefs(bnd_cond,sf%agx(j-1),sf%bmgx(j-1),sf%bpgx(j-1),sf%sgx(j-1),sf%clight*dt,sf%dx, &
                              lsigmax_next(j-1),lsigmax(j-1), &
                              sb_coef,which,sf%pml_method)
           end if
          end do
          if (sf%pml_method==1) then
            sf%bmfx(1:sf%nx+sf%nxguard)=-sf%bmfx(1:sf%nx+sf%nxguard)
            sf%bpfx(1:sf%nx+sf%nxguard)=-sf%bpfx(1:sf%nx+sf%nxguard)
            sf%bmgx(0:sf%nx+sf%nxguard-1)=-sf%bmgx(0:sf%nx+sf%nxguard-1)
            sf%bpgx(0:sf%nx+sf%nxguard-1)=-sf%bpgx(0:sf%nx+sf%nxguard-1)
          end if
       end select
    end if


    if (sf%lsy/=0) then
      sigmay=0.
      sigmay_next=0.
      lsigmay=0.
      lsigmay_next=0.
      do j = jmin, sf%ny+sf%nyguard
        if (j>=jmin) then
          sigmay(j)      = (sf%smaxy/sf%dy)*(REAL(j-jmin,8)/sf%sdeltay)**sf%nny
          if (sf%l_nodalgrid) then
            sigmay_next(j) = (sf%smaxy/sf%dy)*((REAL(j-jmin,8)+1.)/sf%sdeltay)**sf%nny
          else
            sigmay_next(j) = (sf%smaxy/sf%dy)*((REAL(j-jmin,8)+0.5)/sf%sdeltay)**sf%nny
          end if
        end if
        if (j<sf%ny+sf%nyguard) then
          lsigmay(sf%ny-j) = sigmay(j)
          lsigmay_next(sf%ny-j-1) = sigmay_next(j)
        end if
      end do
      select case(sf%lsy)
        case(1)
          do j = 0, sf%ny+sf%nyguard-1
           if (sf%l_nodalgrid) then
            call assign_coefs(bnd_cond,sf%afy(j),sf%bpfy(j),sf%bmfy(j),sf%sfy(j),sf%clight*dt,sf%dy, &
                              sigmay(j),sigmay_next(j+1), &
                              sb_coef,which,sf%pml_method)
           else
            call assign_coefs(bnd_cond,sf%afy(j),sf%bpfy(j),sf%bmfy(j),sf%sfy(j),sf%clight*dt,sf%dy, &
                              sigmay(j),sigmay_next(j), &
                              sb_coef,which,sf%pml_method)
           end if
            call assign_coefs(bnd_cond,sf%agy(j),sf%bpgy(j),sf%bmgy(j),sf%sgy(j),sf%clight*dt,sf%dy, &
                              sigmay_next(j),sigmay(j+1), &
                              sb_coef,which,sf%pml_method)
          end do
        case(-1)
          do j = sf%ny, -sf%nyguard+1, -1
            call assign_coefs(bnd_cond,sf%afy(j),sf%bmfy(j),sf%bpfy(j),sf%sfy(j),sf%clight*dt,sf%dy, &
                              lsigmay(j),lsigmay_next(j-1), &
                              sb_coef,which,sf%pml_method)
           if (sf%l_nodalgrid) then
            call assign_coefs(bnd_cond,sf%agy(j),sf%bmgy(j),sf%bpgy(j),sf%sgy(j),sf%clight*dt,sf%dy, &
                              lsigmay_next(j-1),lsigmay(j), &
                              sb_coef,which,sf%pml_method)
           else
            call assign_coefs(bnd_cond,sf%agy(j-1),sf%bmgy(j-1),sf%bpgy(j-1),sf%sgy(j-1),sf%clight*dt,sf%dy, &
                              lsigmay_next(j-1),lsigmay(j-1), &
                              sb_coef,which,sf%pml_method)
           end if
          end do
          if (sf%pml_method==1) then
            sf%bmfy(1:sf%ny+sf%nyguard)=-sf%bmfy(1:sf%ny+sf%nyguard)
            sf%bpfy(1:sf%ny+sf%nyguard)=-sf%bpfy(1:sf%ny+sf%nyguard)
            sf%bmgy(0:sf%ny+sf%nyguard-1)=-sf%bmgy(0:sf%ny+sf%nyguard-1)
            sf%bpgy(0:sf%ny+sf%nyguard-1)=-sf%bpgy(0:sf%ny+sf%nyguard-1)
          end if
       end select
    end if

    if (sf%lsz/=0) then
      sigmaz=0.
      sigmaz_next=0.
      lsigmaz=0.
      lsigmaz_next=0.
      do j = jmin, sf%nz+sf%nzguard
        if (j>=jmin) then
          sigmaz(j)      = (sf%smaxz/sf%dz)*(REAL(j-jmin,8)/sf%sdeltaz)**sf%nnz
          if (sf%l_nodalgrid) then
            sigmaz_next(j) = (sf%smaxz/sf%dz)*((REAL(j-jmin,8)+1.)/sf%sdeltaz)**sf%nnz
          else
            sigmaz_next(j) = (sf%smaxz/sf%dz)*((REAL(j-jmin,8)+0.5)/sf%sdeltaz)**sf%nnz
          end if
        end if
        if (j<sf%nz+sf%nzguard) then
          lsigmaz(sf%nz-j) = sigmaz(j)
          lsigmaz_next(sf%nz-j-1) = sigmaz_next(j)
        end if
      end do
      select case(sf%lsz)
        case(1)
          do j = 0, sf%nz+sf%nzguard-1
           if (sf%l_nodalgrid) then
            call assign_coefs(bnd_cond,sf%afz(j),sf%bpfz(j),sf%bmfz(j),sf%sfz(j),sf%clight*dt,sf%dz, &
                              sigmaz(j),sigmaz_next(j+1), &
                              sb_coef,which,sf%pml_method)
           else
            call assign_coefs(bnd_cond,sf%afz(j),sf%bpfz(j),sf%bmfz(j),sf%sfz(j),sf%clight*dt,sf%dz, &
                              sigmaz(j),sigmaz_next(j), &
                              sb_coef,which,sf%pml_method)
           end if
            call assign_coefs(bnd_cond,sf%agz(j),sf%bpgz(j),sf%bmgz(j),sf%sgz(j),sf%clight*dt,sf%dz,sigmaz_next(j),sigmaz(j+1), &
                              sb_coef,which,sf%pml_method)
          end do
        case(-1)
          do j = sf%nz, -sf%nzguard+1, -1
            call assign_coefs(bnd_cond,sf%afz(j),sf%bmfz(j),sf%bpfz(j),sf%sfz(j),sf%clight*dt,sf%dz, &
                              lsigmaz(j),lsigmaz_next(j-1), &
                              sb_coef,which,sf%pml_method)
           if (sf%l_nodalgrid) then
            call assign_coefs(bnd_cond,sf%agz(j),sf%bmgz(j),sf%bpgz(j),sf%sgz(j),sf%clight*dt,sf%dz, &
                              lsigmaz_next(j-1),lsigmaz(j), &
                              sb_coef,which,sf%pml_method)
           else
            call assign_coefs(bnd_cond,sf%agz(j-1),sf%bmgz(j-1),sf%bpgz(j-1),sf%sgz(j-1),sf%clight*dt,sf%dz, &
                              lsigmaz_next(j-1),lsigmaz(j-1), &
                              sb_coef,which,sf%pml_method)
           end if
          end do
          if (sf%pml_method==1) then
            sf%bmfz(1:sf%nz+sf%nzguard)=-sf%bmfz(1:sf%nz+sf%nzguard)
            sf%bpfz(1:sf%nz+sf%nzguard)=-sf%bpfz(1:sf%nz+sf%nzguard)
            sf%bmgz(0:sf%nz+sf%nzguard-1)=-sf%bmgz(0:sf%nz+sf%nzguard-1)
            sf%bpgz(0:sf%nz+sf%nzguard-1)=-sf%bpgz(0:sf%nz+sf%nzguard-1)
          end if
      end select
    end if

    return 
  end subroutine set_bndcoeffsem3d
    

!************* SUBROUTINE assign_coefs  ********

subroutine assign_coefs(bnd_cond,a,bp,bm,s,dt,dx,sigma,sigma_next,coef_sigmab,which,method)
implicit none
REAL(kind=8), INTENT(OUT) :: a,bp,bm,s
REAL(kind=8), INTENT(IN) :: dt,dx,sigma,sigma_next,coef_sigmab
INTEGER(ISZ),INTENT(IN) :: bnd_cond, which, method

REAL(kind=8) :: sigma_local, sigmab, sigmab_next, tp, tpp, tm, tmm, g, gp, gm

  if (method==2 .and. (bnd_cond/=pml .and. bnd_cond/=pml_sadjusted)) then
    write(0,*) 'Error: bnd_cond must be either pml or pml_sadjusted with method=2.'
    call abort()
  end if
  tp  = EXP(-sigma*0.5*dx)
  tpp  = EXP(-sigma_next*0.5*dx)
  select case (bnd_cond)
    case (pml_user)
      return
    case (pml,pml_sadjusted)
      IF(bnd_cond==pml) then
        sigma_local = sigma
      else
        sigma_local = MIN(1.e15,abs(tpp-1./tp)/dx)
      END if
      if (method==1) then
        IF(sigma_local == 0.) then
        ! --- end of mesh
          a  =  1.
          bp =  dt / dx
        else
          a  =  EXP(-sigma_local*dt) ! one can use the exponential intergration or
!                                       direct differenciation as below for about the same effect.
!          a  = (1.-sigma_local*dt/2)/(1.+sigma_local*dt/2)
          bp =  (1.-a)/(sigma_local*dx)
        END if
        bm =  -bp
      else
        s  =  EXP(-sigma_local*dt)
        
      end if
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

  if (method==1) then
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
  end if

END subroutine assign_coefs

end module mod_emfield3d

  subroutine init_splitfield(sf, nx, ny, nz, nxguard, nyguard, nzguard, dt, dx, dy, dz, xmin, ymin, zmin, clight, lsx, lsy, lsz, &
                             nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz, l_1dz, l_2dxz, l_2drz, &
                             norderx,nordery,norderz,xcoefs,ycoefs,zcoefs, &
                             l_nodalgrid, pml_method)
    use mod_emfield3d
    TYPE(EM3D_SPLITYEEFIELDtype) :: sf
    INTEGER(ISZ), INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard, nnx, nny, nnz, lsx, lsy, lsz, pml_method, &     
                                norderx,nordery,norderz
    REAL(kind=8), INTENT(IN) :: dt, dx, dy, dz, clight, smaxx, smaxy, smaxz, sdeltax, sdeltay, sdeltaz, xmin, ymin, zmin, &
                                xcoefs(norderx/2),ycoefs(nordery/2),zcoefs(norderz/2)
    integer(ISZ) :: j
    logical(ISZ) :: l_1dz, l_2dxz, l_2drz, l_nodalgrid
    
    sf%pml_method = pml_method
    sf%l_nodalgrid = l_nodalgrid
    sf%nx = nx
    sf%ny = ny
    sf%nz = nz
    sf%nxguard = nxguard
    sf%nyguard = nyguard
    sf%nzguard = nzguard
    sf%dx = dx
    sf%dy = dy
    sf%dz = dz
    sf%xmin = xmin
    sf%xmax = xmin+nx*dx
    sf%ymax = ymin+ny*dy
    sf%zmax = zmin+nz*dz
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
    sf%l_1dz = l_1dz
    sf%l_2dxz = l_2dxz
    sf%l_2drz = l_2drz
    sf%norderx = norderx
    sf%nordery = nordery
    sf%norderz = norderz
    call EM3D_SPLITYEEFIELDtypeallot(sf)
    sf%xcoefs = xcoefs
    sf%ycoefs = ycoefs
    sf%zcoefs = zcoefs

    return 
  end subroutine init_splitfield
    
subroutine push_em3d_e(f,dt)
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

!if (f%spectral) then
!  call callpythonfunc("push_spectral","em3dsolverFFT")
!  return
!end if

if (f%theta_damp/=0.) then
!  f%exbar = f%theta_damp*f%exbar + f%exold
!  f%eybar = f%theta_damp*f%eybar + f%eyold
!  f%ezbar = f%theta_damp*f%ezbar + f%ezold
  f%exold = f%ex
  f%eyold = f%ey
  f%ezold = f%ez
  f%exbar = (1.-0.5*f%theta_damp)*f%ex+0.5*f%theta_damp*f%exbar
  f%eybar = (1.-0.5*f%theta_damp)*f%ey+0.5*f%theta_damp*f%eybar
  f%ezbar = (1.-0.5*f%theta_damp)*f%ez+0.5*f%theta_damp*f%ezbar
end if

select case(f%stencil)
   ! Choose the kind of stencil that is used for the E push

case(0,1,3) ! Yee stencil on the E push
   ! (Note : Yee stencil on the B push in case 0
   ! Cole-Karkkainen stencil on the B push in case 1
   ! Lehe stencil on the B push in case 3)
 if (f%sigmae==0.) then
  if(f%nconds>0) then 
   if (f%l_macroscopic) then
      call push_em3d_evec_macroscopic(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                       mudt,dtsdx,dtsdy,dtsdz, &
                       f%nx,f%ny,f%nz, &
                       f%nxguard,f%nyguard,f%nzguard, &
                       f%l_2dxz,f%l_2drz,f%xmin,f%zmin,f%dx,f%dz, &
                       f%sigmax,f%sigmay,f%sigmaz, &
                       f%epsix,f%epsiy,f%epsiz, &
                       f%mux,f%muy,f%muz)
   else
      call push_em3d_evec_cond(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                       mudt,dtsdx,dtsdy,dtsdz, &
                       f%nx,f%ny,f%nz, &
                       f%nxguard,f%nyguard,f%nzguard, &
                       f%l_2dxz,f%l_2drz,f%xmin,f%zmin,f%dx,f%dz,f%incond)
   endif
  else
   if ((f%norderx==2) .and. (f%nordery==2) .and. (f%norderz==2) .and. .not. f%l_nodalgrid) then
     call push_em3d_evec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                          mudt,dtsdx,dtsdy,dtsdz, &
                          f%nx,f%ny,f%nz, &
                          f%nxguard,f%nyguard,f%nzguard, &
                          f%nxes,f%nyes,f%nzes, &
                          f%l_1dz,f%l_2dxz,f%l_2drz,f%xmin,f%zmin, &
                          f%dx,f%dy,f%dz,f%clight)
     if (f%circ_m>0) &
       call push_em3d_evec_circ(f%ex_circ,f%ey_circ,f%ez_circ, &
                                f%bx_circ,f%by_circ,f%bz_circ,f%J_circ, &
                                mudt,dtsdx,dtsdz,f%nx,f%nz,f%nxguard,f%nzguard, &
                                f%xmin,f%zmin,f%dx,f%dz,f%clight,f%circ_m)
    else
     call push_em3d_evec_norder(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                          mudt,dtsdx*f%xcoefs,dtsdy*f%ycoefs,dtsdz*f%zcoefs, &
                          f%nx,f%ny,f%nz, &
                          f%norderx,f%nordery,f%norderz, &
                          f%nxguard,f%nyguard,f%nzguard, &
                          f%nxes,f%nyes,f%nzes, &
                          f%l_1dz,f%l_2dxz,f%l_2drz,f%l_nodalgrid, &
                          f%xmin,f%zmin, &
                          f%dx,f%dy,f%dz,f%clight)
    end if
  end if
 endif

case(2)  ! Cole-Karkkainen stencil on the E push (Note : Yee stencil on the B push )
  call push_em3d_kyeevec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                         mudt,dtsdx,dtsdy,dtsdz, &
                         f%nx,f%ny,f%nz, &
                         f%nxguard,f%nyguard,f%nzguard,f%E_inz_pos,f%Ex_inz,f%Ey_inz,f%l_2dxz,f%zmin,f%dz)

end select

return
end subroutine push_em3d_e

subroutine push_em3d_evec(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz, &
                          nxguard,nyguard,nzguard,nxs,nys,nzs, &
                          l_1dz,l_2dxz,l_2drz,xmin,zmin,dx,dy,dz,clight)
integer :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz,xmin,zmin,dx,dy,dz,clight
integer(ISZ) :: j,k,l
logical(ISZ) :: l_1dz,l_2dxz,l_2drz
real(kind=8) :: w,zlaser,rd,ru
!real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard) :: Ez_inz

! --- NOTE: if l_2drz is TRUE, then l_2dxz is TRUE
if (.not. l_2dxz) then ! --- 3D XYZ
  ! advance Ex
  do l = -nzs, nz+nzs
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-1
      Ex(j,k,l) = Ex(j,k,l) + dtsdy * (Bz(j,k,l)   - Bz(j,k-1,l  )) &
                            - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
   end do
  end do

  ! advance Ey
  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
   end do
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                            - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
   end do
  end do

else ! --- now 1D Z, 2D XZ or RZ

 if (l_1dz) then ! 1D Z

  j = 0
  k = 0
  ! advance Ex
  do l = -nzs, nz+nzs
      Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
  end do

  ! advance Ey
  do l = -nzs, nz+nzs
      Ey(j,k,l) = Ey(j,k,l) + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-1
      Ez(j,k,l) = Ez(j,k,l) - mudt  * CJ(j,k,l,3)
  end do

 else if (.not. l_2drz) then ! 2D XZ

  k = 0
  ! advance Ex
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-1
      Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
  end do

  ! advance Ey
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
  end do

 else ! l_2drz=True

  k = 0
  ! advance Er
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-1
      Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
  end do

  ! advance Etheta
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
     if (j/=0) &
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
    j = 0
    if (xmin/=0.) then
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end if
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs
     if (j/=0) then
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
     end if
    end do
    j = 0
    if (xmin==0.) then
      Ez(j,k,l) = Ez(j,k,l) + 4.*dtsdx * By(j,k,l)  &
                            - mudt  * CJ(j,k,l,3)
    else
      ru = 1.+0.5/(xmin/dx)
      rd = 1.-0.5/(xmin/dx)
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    endif
  end do
 end if
end if


return
end subroutine push_em3d_evec

subroutine push_em3d_evec_norder(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz, &
                          norderx,nordery,norderz, &
                          nxguard,nyguard,nzguard,nxs,nys,nzs, &
                          l_1dz,l_2dxz,l_2drz,l_nodalgrid, &
                          xmin,zmin,dx,dy,dz,clight)
integer :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN) :: mudt,dtsdx(norderx/2),dtsdy(norderx/2),dtsdz(norderx/2),xmin,zmin,dx,dy,dz,clight
integer(ISZ) :: i,j,k,l,ist
logical(ISZ) :: l_1dz,l_2dxz,l_2drz,l_nodalgrid
real(kind=8) :: w,zlaser,rd,ru
!real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard) :: Ez_inz

if (l_nodalgrid) then
  ist = 0
else
  ist = 1
end if

! --- NOTE: if l_2drz is TRUE, then l_2dxz is TRUE
if (.not. l_2dxz) then ! --- 3D XYZ
  ! advance Ex
  do l = -nzs, nz+nzs
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-ist
      Ex(j,k,l) = Ex(j,k,l) - mudt  * CJ(j,k,l,1)
      do i = 1, nordery/2
        Ex(j,k,l) = Ex(j,k,l) + dtsdy(i) * (Bz(j,k+i-ist,l)   - Bz(j,k-i,l  ))
      end do
      do i = 1, norderz/2
        Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i-ist)   - By(j,k  ,l-i)) 
      end do
    end do
   end do
  end do

  ! advance Ey
  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs
      Ey(j,k,l) = Ey(j,k,l) - mudt  * CJ(j,k,l,2)
      do i = 1, norderx/2
        Ey(j,k,l) = Ey(j,k,l) - dtsdx(i) * (Bz(j+i-ist,k,l)   - Bz(j-i,k,l))
      end do
      do i = 1, norderz/2
        Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist)   - Bx(j,k,l-i))
      end do
    end do
   end do
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs
      Ez(j,k,l) = Ez(j,k,l) - mudt  * CJ(j,k,l,3)
      do i = 1, norderx/2
        Ez(j,k,l) = Ez(j,k,l) + dtsdx(i) * (By(j+i-ist,k,l) - By(j-i,k  ,l))
      end do
      do i = 1, nordery/2
        Ez(j,k,l) = Ez(j,k,l) - dtsdy(i) * (Bx(j,k+i-ist,l) - Bx(j  ,k-i,l))
      end do
    end do
   end do
  end do

else ! --- now 1D Z, 2D XZ or RZ

 if (l_1dz) then ! 1D Z

  j = 0
  k = 0
  ! advance Ex
  do l = -nzs, nz+nzs
      Ex(j,k,l) = Ex(j,k,l) - mudt  * CJ(j,k,l,1)
      do i = 1, norderz/2
        Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i-ist)   - By(j,k  ,l-i))
      end do
  end do

  ! advance Ey
  do l = -nzs, nz+nzs
      Ey(j,k,l) = Ey(j,k,l) - mudt  * CJ(j,k,l,2)
      do i = 1, norderz/2
        Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist)   - Bx(j,k,l-i)) 
      end do
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-1
      Ez(j,k,l) = Ez(j,k,l) - mudt  * CJ(j,k,l,3)
  end do

 else if (.not. l_2drz) then ! 2D XZ

  k = 0
  ! advance Ex
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-ist
      Ex(j,k,l) = Ex(j,k,l) - mudt  * CJ(j,k,l,1)
      do i = 1, norderz/2
        Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i-ist)   - By(j,k  ,l-i)) 
      end do
    end do
  end do

  ! advance Ey
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
      Ey(j,k,l) = Ey(j,k,l) - mudt  * CJ(j,k,l,2)
      do i = 1, norderx/2
        Ey(j,k,l) = Ey(j,k,l) - dtsdx(i) * (Bz(j+i-ist,k,l)   - Bz(j-i,k,l)) 
      end do
      do i = 1, norderz/2 
        Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist)   - Bx(j,k,l-i)) 
      end do
    end do
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs
      Ez(j,k,l) = Ez(j,k,l) - mudt  * CJ(j,k,l,3)
      do i = 1, norderx/2
        Ez(j,k,l) = Ez(j,k,l) + dtsdx(i) * (By(j+i-ist,k,l) - By(j-i,k  ,l))
      end do
    end do
  end do

 else ! l_2drz=True

  if (norderx.ne.2) then
    write(0,*) 'Error: norderx>2 not supported in RZ axisymmetric mode'
    call abort()
  end if

  k = 0
  ! advance Er
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-ist
      Ex(j,k,l) = Ex(j,k,l) - mudt  * CJ(j,k,l,1)
      do i = 1, norderz/2
        Ex(j,k,l) = Ex(j,k,l) - dtsdz(i) * (By(j,k,l+i-ist)   - By(j,k  ,l-i))
      end do
    end do
  end do

  ! advance Etheta
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs
     if (j/=0) then
        Ey(j,k,l) = Ey(j,k,l) - dtsdx(1) * (Bz(j,k,l) - Bz(j-1,k,l)) &
                              - mudt  * CJ(j,k,l,2)
        do i = 1, norderz/2
           Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist) - Bx(j,k,l-i)) 
        end do
      end if
    end do
    j = 0
    if (xmin/=0.) then
      Ey(j,k,l) = Ey(j,k,l) - dtsdx(1) * (Bz(j,k,l) - Bz(j-1,k,l)) &
                            - mudt  * CJ(j,k,l,2)
      do i = 1, norderz/2
         Ey(j,k,l) = Ey(j,k,l) + dtsdz(i) * (Bx(j,k,l+i-ist) - Bx(j,k,l-i))
      end do
    end if
  end do

  ! advance Ez 
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs
     if (j/=0) then
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      Ez(j,k,l) = Ez(j,k,l) + dtsdx(1) * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
     end if
    end do
    j = 0
    if (xmin==0.) then
      Ez(j,k,l) = Ez(j,k,l) + 4.*dtsdx(1) * By(j,k,l)  &
                            - mudt  * CJ(j,k,l,3)
    else
      ru = 1.+0.5/(xmin/dx)
      rd = 1.-0.5/(xmin/dx)
      Ez(j,k,l) = Ez(j,k,l) + dtsdx(1) * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    endif
  end do
 end if
end if


return
end subroutine push_em3d_evec_norder

subroutine push_em3d_evec_circ(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdz,nx,nz, &
                          nxguard,nzguard, &
                          xmin,zmin,dx,dz,clight,circ_m)
integer :: nx,nz,nxguard,nzguard,circ_m
complex(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m) :: ex,ey,ez,bx,by,bz
complex(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,3,1:circ_m) :: CJ
real(kind=8), intent(IN) :: mudt,dtsdx,dtsdz,xmin,zmin,dx,dz,clight
integer(ISZ) :: j,l,m
real(kind=8) :: w,r,rd,ru,dt
complex(kind=8) :: i=(0.,1.)

  ! ===============================
  !             2-D RZ multipole
  ! ===============================

  dt = dtsdx*dx
  do m = 1, circ_m

     ! advance Er
     do l = 0, nz
        do j = 0, nx-1
           r = xmin+j*dx+0.5*dx
           Ex(j,l,m) = Ex(j,l,m) - i*m*dt*Bz(j,l,m)/r &
                - dtsdz * (By(j,l,m)   - By(j  ,l-1,m)) &
                - mudt  * CJ(j,l,1,m)
        end do
     end do

     ! advance Etheta
     do l = 0, nz
        do j = 0, nx
           if ( j==0 .and. xmin==0 ) then
              if( .not. m == 1 ) then
                 ! Etheta should remain 0 on axis, for modes different than m=1
                 Ey(j,l,m) = 0
              else ! Mode m=1
                 ! The bulk equation could in principle be used here since it does not diverge
                 ! on axis. However, it typically gives poore results e.g. for the propagation
                 ! of a laser pulse (The field is spuriously reduced on axis.) For this reason
                 ! a modified on-axis condition is used here : we use the fact that
                 ! Etheta(r=0,m=1) should equal -iEr(r=0,m=1), for the fields Ex and Ey to be
                 ! independent of theta at r=0. Now with linear interpolation :
                 ! Er(r=0,m=1) = 0.5*[Er(r=dr/2,m=1)+Er(r=-dr/2,m=1)]
                 ! And using the rule applying for the guards cells (see em3d_applybc_e)
                 ! Er(r=-dr/2,m=1) = Er(r=dr/2,m=1). Thus :
                 Ey(j,l,m) = -i*Ex(j,l,m)
              endif
           else
              ! Equation used in the bulk of the grid
              Ey(j,l,m) = Ey(j,l,m) - dtsdx * (Bz(j,l,m) - Bz(j-1,l,m)) &
                   + dtsdz * (Bx(j,l,m) - Bx(j,l-1,m)) &
                   - mudt * CJ(j,l,2,m)
           endif
        end do
     end do
     
     ! advance Ez 
     do l = 0, nz-1
        do j = 0, nx
           if ( j==0 .and. xmin==0 ) then
              ! Ez should remain 0 on axis, for modes with m>0,
              ! but the bulk equation does not necessarily ensure this.
              Ez(j,l,m) = 0.
           else
              ! Equation used in the bulk of the grid
              ru = 1.+0.5/(xmin/dx+j)
              rd = 1.-0.5/(xmin/dx+j)
              r = xmin+j*dx
              Ez(j,l,m) = Ez(j,l,m) + dtsdx * (ru*By(j,l,m) - rd*By(j-1  ,l,m)) &
                   + i*m*dt*Bx(j,l,m)/r &
                   - mudt  * CJ(j,l,3,m)
           end if
        end do
     end do

  end do
return
end subroutine push_em3d_evec_circ

subroutine push_em3d_evec_cond(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz, &
                          nxguard,nyguard,nzguard,l_2dxz,l_2drz,xmin,zmin,dx,dz,incond)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
logical(ISZ), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: incond
real(kind=8), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz,xmin,zmin,dx,dz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz,l_2drz
real(kind=8) :: w,zlaser,rd,ru

! --- NOTE: if l_2drz is TRUE, then l_2dxz is TRUE
if (.not. l_2dxz) then ! --- 3D XYZ
  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      if (.not.incond(j,k,l) .or. .not.incond(j+1,k,l)) &
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
      if (.not.incond(j,k,l) .or. .not.incond(j,k+1,l)) &
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
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                            - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
   end do
  end do

else ! --- now 2D XZ or RZ

 if (.not. l_2drz) then ! 2D XZ

  k = 0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      if (.not.incond(j,k,l) .or. .not.incond(j+1,k,l)) &
      Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
  end do

  ! advance Ey
  do l = 0, nz
    do j = 0, nx
      if (.not.incond(j,k,l)) &
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
  end do

 else ! l_2drz=True

  k = 0
  ! advance Er
  do l = 0, nz
    do j = 0, nx-1
      if (.not.incond(j,k,l) .or. .not.incond(j+1,k,l)) &
      Ex(j,k,l) = Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                            - mudt  * CJ(j,k,l,1)
    end do
  end do

  ! advance Etheta
  do l = 0, nz
    do j = 1, nx
      if (.not.incond(j,k,l)) &
      Ey(j,k,l) = Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                            + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                            - mudt  * CJ(j,k,l,2)
    end do
    j = 0
    if (.not.incond(j,k,l)) &
    Ey(j,k,l) = Ey(j,k,l) - 2.*dtsdx * Bz(j,k,l) &
                          + dtsdz * (Bx(j,k,l)    - Bx(j,k,l-1)) &
                          - mudt  * CJ(j,k,l,2)
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 1, nx
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    end do
    j = 0
    if (xmin==0.) then
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + 4.*dtsdx * By(j,k,l)  &
                            - mudt  * CJ(j,k,l,3)
    else
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                            - mudt  * CJ(j,k,l,3)
    end if
  end do
 end if
end if


return
end subroutine push_em3d_evec_cond

subroutine push_em3d_evec_macroscopic(ex,ey,ez,bx,by,bz,CJ,mu0dt0,dt0sdx,dt0sdy,dt0sdz,nx,ny,nz, &
                          nxguard,nyguard,nzguard,l_2dxz,l_2drz,xmin,zmin,dx,dz, &
                          sigmax,sigmay,sigmaz,epsix,epsiy,epsiz,mux,muy,muz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: sigmax,sigmay,sigmaz, &
                                                                                                    epsix,epsiy,epsiz,mux,muy,muz
real(kind=8), intent(IN) :: mu0dt0,dt0sdx,dt0sdy,dt0sdz,xmin,zmin,dx,dz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz,l_2drz
real(kind=8) :: w,zlaser,rd,ru,a,b,mu0dt,dtsdx,dtsdy,dtsdz

! --- NOTE: if l_2drz is TRUE, then l_2dxz is TRUE
if (.not. l_2dxz) then ! --- 3D XYZ
  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
     if (sigmax(j,k,l)<0.) then
       Ex(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmax(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsix(j,k,l)*b)
      dtsdy = dt0sdy/(epsix(j,k,l)*mux(j,k,l)*b)
      dtsdz = dt0sdz/(epsix(j,k,l)*mux(j,k,l)*b)
      Ex(j,k,l) = a*Ex(j,k,l) + dtsdy * (Bz(j,k,l)   - Bz(j,k-1,l  )) &
                              - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                              - mu0dt  * CJ(j,k,l,1)
     end if
    end do
   end do
  end do

  ! advance Ey
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
     if (sigmay(j,k,l)<0.) then
       Ey(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmay(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiy(j,k,l)*b)
      dtsdx = dt0sdx/(epsiy(j,k,l)*muy(j,k,l)*b)
      dtsdz = dt0sdz/(epsiy(j,k,l)*muy(j,k,l)*b)
      Ey(j,k,l) = a*Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                              + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                              - mu0dt  * CJ(j,k,l,2)
     end if
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
     if (sigmaz(j,k,l)<0.) then
       Ez(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmaz(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiz(j,k,l)*b)
      dtsdx = dt0sdx/(epsiz(j,k,l)*muz(j,k,l)*b)
      dtsdy = dt0sdy/(epsiz(j,k,l)*muz(j,k,l)*b)
      Ez(j,k,l) = a*Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                              - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                              - mu0dt  * CJ(j,k,l,3)
     end if
    end do
   end do
  end do

else ! --- now 2D XZ or RZ

 if (.not. l_2drz) then ! 2D XZ

  k = 0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
     if (sigmax(j,k,l)<0.) then
       Ex(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmax(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsix(j,k,l)*b)
      dtsdz = dt0sdz/(epsix(j,k,l)*mux(j,k,l)*b)
      Ex(j,k,l) = a*Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                              - mu0dt  * CJ(j,k,l,1)
     end if
    end do
  end do

  ! advance Ey
  do l = 0, nz
    do j = 0, nx
      if (sigmay(j,k,l)<0.) then
       Ey(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmay(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiy(j,k,l)*b)
      dtsdx = dt0sdx/(epsiy(j,k,l)*muy(j,k,l)*b)
      dtsdz = dt0sdz/(epsiy(j,k,l)*muy(j,k,l)*b)
      Ey(j,k,l) = a*Ey(j,k,l) - dtsdx * (Bz(j,k,l)   - Bz(j-1,k,l)) &
                              + dtsdz * (Bx(j,k,l)   - Bx(j,k,l-1)) &
                              - mu0dt  * CJ(j,k,l,2)
     end if
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      a = 0.5*mu0dt0*sigmaz(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiz(j,k,l)*b)
      dtsdx = dt0sdx/(epsiz(j,k,l)*muz(j,k,l)*b)
     if (sigmaz(j,k,l)<0.) then
       Ez(j,k,l)=0.
     else
      Ez(j,k,l) = a*Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                              - mu0dt  * CJ(j,k,l,3)
     end if
    end do
  end do

 else ! l_2drz=True

  k = 0
  ! advance Er
  do l = 0, nz
    do j = 0, nx-1
     if (sigmax(j,k,l)<0.) then
       Ex(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmax(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsix(j,k,l)*b)
      dtsdy = dt0sdy/(epsix(j,k,l)*mux(j,k,l)*b)
      dtsdz = dt0sdz/(epsix(j,k,l)*mux(j,k,l)*b)
      Ex(j,k,l) = a*Ex(j,k,l) - dtsdz * (By(j,k,l)   - By(j,k  ,l-1)) &
                              - mu0dt  * CJ(j,k,l,1)
     end if
    end do
  end do

  ! advance Etheta
  do l = 0, nz
    do j = 1, nx
     if (sigmay(j,k,l)<0.) then
       Ey(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmay(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiy(j,k,l)*b)
      dtsdx = dt0sdx/(epsiy(j,k,l)*muy(j,k,l)*b)
      dtsdz = dt0sdz/(epsiy(j,k,l)*muy(j,k,l)*b)
      Ey(j,k,l) = a*Ey(j,k,l) - dtsdx * (Bz(j,k,l) - Bz(j-1,k,l)) &
                              + dtsdz * (Bx(j,k,l) - Bx(j,k,l-1)) &
                              - mu0dt  * CJ(j,k,l,2)
     end if
    end do
    j = 0
    a = 0.5*mu0dt0*sigmay(j,k,l)
    b = (1.+a)
    a = (1-a)/(1+a)
    mu0dt = mu0dt0/(epsiy(j,k,l)*b)
    dtsdx = dt0sdx/(epsiy(j,k,l)*muy(j,k,l)*b)
    dtsdz = dt0sdz/(epsiy(j,k,l)*muy(j,k,l)*b)
    if (sigmay(j,k,l)<0.) then
       Ey(j,k,l)=0.
    else
      Ey(j,k,l) = a*Ey(j,k,l) - 2.*dtsdx * Bz(j,k,l) &
                              + dtsdz * (Bx(j,k,l)    - Bx(j,k,l-1)) &
                              - mu0dt  * CJ(j,k,l,2)
    end if
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 1, nx
     if (sigmaz(j,k,l)<0.) then
       Ez(j,k,l)=0.
     else
      a = 0.5*mu0dt0*sigmaz(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiz(j,k,l)*b)
      dtsdx = dt0sdx/(epsiz(j,k,l)*muz(j,k,l)*b)
      dtsdy = dt0sdy/(epsiz(j,k,l)*muz(j,k,l)*b)
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      Ez(j,k,l) = a*Ez(j,k,l) + dtsdx * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                              - mu0dt  * CJ(j,k,l,3)
     end if
    end do
    j = 0
    if (sigmaz(j,k,l)<0.) then
      Ez(j,k,l)=0.
    else
     if (xmin==0.) then
      a = 0.5*mu0dt0*sigmaz(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiz(j,k,l)*b)
      dtsdx = dt0sdx/(epsiz(j,k,l)*muz(j,k,l)*b)
      dtsdy = dt0sdy/(epsiz(j,k,l)*muz(j,k,l)*b)
      Ez(j,k,l) = a*Ez(j,k,l) + 4.*dtsdx * By(j,k,l)  &
                              - mu0dt  * CJ(j,k,l,3)
     else
      a = 0.5*mu0dt0*sigmaz(j,k,l)
      b = (1.+a)
      a = (1-a)/(1+a)
      mu0dt = mu0dt0/(epsiz(j,k,l)*b)
      dtsdx = dt0sdx/(epsiz(j,k,l)*muz(j,k,l)*b)
      dtsdy = dt0sdy/(epsiz(j,k,l)*muz(j,k,l)*b)
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      Ez(j,k,l) = a*Ez(j,k,l) + dtsdx * (ru*By(j,k,l) - rd*By(j-1,k  ,l)) &
                              - mu0dt  * CJ(j,k,l,3)
     end if
    end if
  end do
 end if
end if

return
end subroutine push_em3d_evec_macroscopic

subroutine push_em3d_kyeevec(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz, &
                             nx,ny,nz,nxguard,nyguard,nzguard,e_inz_pos,Ex_inz,Ey_inz,l_2dxz,zmin,dz)
use EM3D_kyee
implicit none
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN) :: E_inz_pos,zmin,dz
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard) :: Ex_inz,Ey_inz
logical(ISZ) :: l_2dxz

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt,E_inz_angle,zlaser,w

if (.not.l_2dxz) then
  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) + alphay*dtsdy * (Bz(j  ,k,l  ) - Bz(j  ,k-1,l  )) &
                            + betayx*dtsdy * (Bz(j+1,k,l  ) - Bz(j+1,k-1,l  ) &
                                           +  Bz(j-1,k,l  ) - Bz(j-1,k-1,l  )) &
                            + betayz*dtsdy * (Bz(j  ,k,l+1) - Bz(j  ,k-1,l+1) &
                                           +  Bz(j  ,k,l-1) - Bz(j  ,k-1,l-1))&
                            + gammay*dtsdy * (Bz(j+1,k,l+1) - Bz(j+1,k-1,l+1) &
                                           +  Bz(j-1,k,l+1) - Bz(j-1,k-1,l+1) &
                                           +  Bz(j+1,k,l-1) - Bz(j+1,k-1,l-1) &
                                           +  Bz(j-1,k,l-1) - Bz(j-1,k-1,l-1)) &
                            - alphaz*dtsdz * (By(j  ,k  ,l) - By(j  ,k  ,l-1)) &
                            + betayx*dtsdy * (By(j+1,k  ,l) - By(j+1,k  ,l-1)  &
                                           +  By(j-1,k  ,l) - By(j-1,k  ,l-1))  &
                            + betayz*dtsdy * (By(j  ,k+1,l) - By(j  ,k+1,l-1)  &
                                           +  By(j  ,k-1,l) - By(j  ,k-1,l-1)) &
                            - gammaz*dtsdz * (By(j+1,k+1,l) - By(j+1,k+1,l-1)  &
                                           +  By(j-1,k+1,l) - By(j-1,k+1,l-1)  &
                                           +  By(j+1,k-1,l) - By(j+1,k-1,l-1)  &
                                           +  By(j-1,k-1,l) - By(j-1,k-1,l-1)) &
                            - 0.5*(alphay+alphaz)*mudt * CJ(j,k,l,1) &
                                     - 0.5*betayx *mudt * (CJ(j+1,k,l  ,1) &
                                                     +CJ(j-1,k,l  ,1)) &
                                     - 0.5*betayz *mudt * (CJ(j  ,k,l+1,1) &
                                                     +CJ(j  ,k,l-1,1)) &
                                     - 0.5*gammay*mudt * (CJ(j+1,k,l+1,1) &
                                                     +CJ(j-1,k,l+1,1) &
                                                     +CJ(j+1,k,l-1,1) &
                                                     +CJ(j-1,k,l-1,1)) &
                                     - 0.5*betazx *mudt * (CJ(j+1,k  ,l,1) &
                                                     +CJ(j-1,k  ,l,1)) &
                                     - 0.5*betazy *mudt * (CJ(j  ,k+1,l,1) &
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
                                - betaxy *dtsdx * (Bz(j,k+1,l  ) - Bz(j-1,k+1,l  ) &
                                               +  Bz(j,k-1,l  ) - Bz(j-1,k-1,l  )) &
                                - betaxz *dtsdx * (Bz(j,k  ,l+1) - Bz(j-1,k  ,l+1) &
                                               +  Bz(j,k  ,l-1) - Bz(j-1,k  ,l-1)) &
                                - gammax*dtsdx * (Bz(j,k+1,l+1) - Bz(j-1,k+1,l+1) &
                                               +  Bz(j,k-1,l+1) - Bz(j-1,k-1,l+1) &
                                               +  Bz(j,k+1,l-1) - Bz(j-1,k+1,l-1) &
                                               +  Bz(j,k-1,l-1) - Bz(j-1,k-1,l-1)) &                              
                                + alphaz*dtsdz * (Bx(j  ,k  ,l) - Bx(j  ,k  ,l-1)) &
                                + betazx *dtsdz * (Bx(j+1,k  ,l) - Bx(j+1,k  ,l-1) &
                                               +  Bx(j-1,k  ,l) - Bx(j-1,k  ,l-1)) &
                                + betazy *dtsdz * (Bx(j  ,k+1,l) - Bx(j  ,k+1,l-1) &
                                               +  Bx(j  ,k-1,l) - Bx(j  ,k-1,l-1)) &
                                + gammaz*dtsdz * (Bx(j+1,k+1,l) - Bx(j+1,k+1,l-1) &
                                               +  Bx(j-1,k+1,l) - Bx(j-1,k+1,l-1) &
                                               +  Bx(j+1,k-1,l) - Bx(j+1,k-1,l-1) &
                                               +  Bx(j-1,k-1,l) - Bx(j-1,k-1,l-1)) &
                                - 0.5*(alphax+alphaz)*mudt * CJ(j,k,l,2) &
                                        - 0.5*betaxy *mudt * (CJ(j,k+1,l  ,2) &
                                                      +  CJ(j,k-1,l  ,2)) &
                                        - 0.5*betaxz *mudt * (CJ(j,k  ,l+1,2) &
                                                      +  CJ(j,k  ,l-1,2)) &
                                        - 0.5*gammax*mudt * (CJ(j,k+1,l+1,2) &
                                                      +  CJ(j,k-1,l+1,2) &
                                                      +  CJ(j,k+1,l-1,2) &
                                                      +  CJ(j,k-1,l-1,2)) &
                                        - 0.5*betazx *mudt * (CJ(j  ,k+1,l,2) &
                                                      +  CJ(j  ,k-1,l,2)) &
                                        - 0.5*betazy *mudt * (CJ(j+1,k  ,l,2) &
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
                                + betaxy*dtsdx * (By(j,k+1,l  ) - By(j-1,k+1,l  ) &
                                                 +  By(j,k-1,l  ) - By(j-1,k-1,l  )) &
                                + betaxz*dtsdx * (By(j,k  ,l+1) - By(j-1,k  ,l+1) &
                                                 +  By(j,k  ,l-1) - By(j-1,k  ,l-1)) &
                                + gammax*dtsdx * (By(j,k+1,l+1) - By(j-1,k+1,l+1) &
                                                 +  By(j,k-1,l+1) - By(j-1,k-1,l+1) &
                                                 +  By(j,k+1,l-1) - By(j-1,k+1,l-1) &
                                                 +  By(j,k-1,l-1) - By(j-1,k-1,l-1)) &
                                - alphay*dtsdy * (Bx(j  ,k,l  ) - Bx(j  ,k-1,l  )) &
                                + betayx*dtsdy * (Bx(j+1,k,l  ) - Bx(j+1,k-1,l  ) &
                                                 +  Bx(j-1,k,l  ) - Bx(j-1,k-1,l  )) &
                                + betayz*dtsdy * (Bx(j  ,k,l+1) - Bx(j  ,k-1,l+1) &
                                                 +  Bx(j  ,k,l-1) - Bx(j  ,k-1,l-1)) &
                                - gammay*dtsdy * (Bx(j+1,k,l+1) - Bx(j+1,k-1,l+1) &
                                                 +  Bx(j-1,k,l+1) - Bx(j-1,k-1,l+1) &
                                                 +  Bx(j+1,k,l-1) - Bx(j+1,k-1,l-1) &
                                                 +  Bx(j-1,k,l-1) - Bx(j-1,k-1,l-1)) &
                                - 0.5*(alphax+alphay)*mudt * CJ(j,k,l,3) &
                                        - 0.5*betaxy *mudt * (CJ(j,k+1,l  ,3) &
                                                      +  CJ(j,k-1,l  ,3)) &
                                        - 0.5*betaxz *mudt * (CJ(j,k  ,l+1,3) &
                                                      +  CJ(j,k  ,l-1,3)) &
                                        - 0.5*gammax*mudt * (CJ(j,k+1,l+1,3) &
                                                      +  CJ(j,k-1,l+1,3) &
                                                      +  CJ(j,k+1,l-1,3) &
                                                      +  CJ(j,k-1,l-1,3)) &
                                        - 0.5*betayx *mudt * (CJ(j+1,k,l  ,3) &
                                                      +  CJ(j-1,k,l  ,3)) &
                                        - 0.5*betayz *mudt * (CJ(j  ,k,l+1,3) &
                                                      +  CJ(j  ,k,l-1,3)) &
                                        - 0.5*gammay*mudt * (CJ(j+1,k,l+1,3) &
                                                      +  CJ(j-1,k,l+1,3) &
                                                      +  CJ(j+1,k,l-1,3) &
                                                      +  CJ(j-1,k,l-1,3))
    end do
   end do
  end do

else
  k = 0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      Ex(j,k,l) = Ex(j,k,l) -     alphaz*dtsdz * (By(j  ,k  ,l) - By(j  ,k  ,l-1)) &
                            -     betazx*dtsdz * (By(j+1,k  ,l) - By(j+1,k  ,l-1)  &
                                               +  By(j-1,k  ,l) - By(j-1,k  ,l-1))  &
                            - alphaz*mudt       * CJ(j,k,l,1) &
                            -     betazx*mudt    * (CJ(j+1,k  ,l,1)+CJ(j-1,k  ,l,1) )
    end do
  end do

  ! advance Ey
  do l = 0, nz
    do j = 0, nx
      Ey(j,k,l) = Ey(j,k,l) - alphax*dtsdx * (Bz(j,k  ,l  ) - Bz(j-1,k  ,l  )) &
                                - betaxz*dtsdx * (Bz(j,k  ,l+1) - Bz(j-1,k  ,l+1) &
                                               +  Bz(j,k  ,l-1) - Bz(j-1,k  ,l-1)) &
                                + alphaz*dtsdz * (Bx(j  ,k  ,l) - Bx(j  ,k  ,l-1)) &
                                + betazx*dtsdz * (Bx(j+1,k  ,l) - Bx(j+1,k  ,l-1) &
                                               +  Bx(j-1,k  ,l) - Bx(j-1,k  ,l-1)) &
                                - 0.5*(alphax+alphaz)*mudt * CJ(j,k,l,2) &
                                        - 0.5*    betaxz *mudt * (CJ(j,k  ,l+1,2) &
                                                      +  CJ(j,k  ,l-1,2)) &
                                        - 0.5*    betazx *mudt * (CJ(j+1,k  ,l,2) &
                                                      +  CJ(j-1,k  ,l,2)) 
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + alphax*dtsdx * (By(j,k  ,l  ) - By(j-1,k  ,l  )) &
                                + betaxz *dtsdx * ( By(j,k  ,l+1) - By(j-1,k  ,l+1) &
                                                 +  By(j,k  ,l-1) - By(j-1,k  ,l-1)) &
                                - alphax*mudt * CJ(j,k,l,3) &
                                        -     betaxz *mudt * (CJ(j,k  ,l+1,3) &
                                                      +  CJ(j,k  ,l-1,3))  
    end do
  end do

end if

return
end subroutine push_em3d_kyeevec

subroutine push_em3d_b(f,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8) :: dt

INTEGER :: j, k, l, which
real(kind=8) :: dtsdx,dtsdy,dtsdz,alpha,beta,gamma

dtsdx = dt/f%dx
dtsdy = dt/f%dy
dtsdz = dt/f%dz

!if (f%spectral) return

if (f%theta_damp/=0.) then
  f%excp = f%ex
  f%eycp = f%ey
  f%ezcp = f%ez
  f%ex = (1.+0.25*f%theta_damp)*f%ex-0.5*f%exold+(0.5-0.25*f%theta_damp)*f%exbar
  f%ey = (1.+0.25*f%theta_damp)*f%ey-0.5*f%eyold+(0.5-0.25*f%theta_damp)*f%eybar
  f%ez = (1.+0.25*f%theta_damp)*f%ez-0.5*f%ezold+(0.5-0.25*f%theta_damp)*f%ezbar
  alpha = (1.+0.5*f%theta_damp)
  beta  = -(1.-0.5*f%theta_damp)*f%theta_damp
  gamma = 0.5*(1.-f%theta_damp)**2*f%theta_damp
!  f%ex = alpha*f%ex + beta*f%exold + gamma*f%exbar
!  f%ey = alpha*f%ey + beta*f%eyold + gamma*f%eybar
!  f%ez = alpha*f%ez + beta*f%ezold + gamma*f%ezbar
end if

select case (f%stencil)
   ! Choose the kind of stencil that is used for the B push
   
case( 0, 2 ) ! Standard Yee stencil on B push
  ! (Note: Cole-Karkkainen stencil on E push in case 2)
  if (f%sigmab==0.) then
   if ((f%norderx==2) .and. (f%nordery==2) .and. (f%norderz==2) .and. .not. f%l_nodalgrid) then
    call push_em3d_bvec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
                        dtsdx,dtsdy,dtsdz, &
                        f%dx,f%dy,f%dz, &
                        f%xmin,f%ymin,f%zmin, &
                        f%nx,f%ny,f%nz, &
                        f%nxguard,f%nyguard,f%nzguard, &
                        f%nxbs,f%nybs,f%nzbs, &
                        f%l_1dz,f%l_2dxz,f%l_2drz)
    if (f%circ_m>0) &
      call push_em3d_bvec_circ(f%ex_circ,f%ey_circ,f%ez_circ, &
                               f%bx_circ,f%by_circ,f%bz_circ, &
                               dtsdx,dtsdz,f%dx,f%dz,f%xmin,f%zmin, &
                               f%nx,f%nz, f%nxguard,f%nzguard,f%circ_m)
  else
    call push_em3d_bvec_norder(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
                        dtsdx*f%xcoefs,dtsdy*f%ycoefs,dtsdz*f%zcoefs, &
                        f%dx,f%dy,f%dz, &
                        f%xmin,f%ymin,f%zmin, &
                        f%nx,f%ny,f%nz, &
                        f%nxguard,f%nyguard,f%nzguard, &
                        f%norderx,f%nordery,f%norderz, &
                        f%nxbs,f%nybs,f%nzbs, &
                        f%l_1dz,f%l_2dxz,f%l_2drz,f%l_nodalgrid)
  end if
endif

case( 1 ) ! Cole-Karkkainen stencil on B push (Note: Yee stencil on E push)
  call push_em3d_kyeebvec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)

case( 3 ) ! Lehe stencil on B push (Note: Yee stencil on E push)
   call push_em3d_lehebvec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
        dtsdx,dtsdy,dtsdz, &
        f%nx,f%ny,f%nz, &
        f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
   
end select

if (f%theta_damp/=0.) then
  f%ex = f%excp
  f%ey = f%eycp
  f%ez = f%ezcp
end if

return
end subroutine push_em3d_b

subroutine push_em3d_bvec(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,dx,dy,dz, &
                          xmin,ymin,zmin,nx,ny,nz,nxguard,nyguard,nzguard, &
                          nxs,nys,nzs,l_1dz,l_2dxz,l_2drz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz,xmin,ymin,zmin,dx,dy,dz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_1dz,l_2dxz,l_2drz
real(kind=8) :: rd, ru

if (.not.l_2dxz) then

  ! advance Bx
  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs
      Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j,k+1,l  ) - Ez(j,k,l)) &
                            + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
    end do
   end do
  end do

  ! advance By
  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-1
      By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &  
                            - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
    end do
   end do
  end do

  ! advance Bz 
  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs-1
      Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j,k,l)) &
                            + dtsdy * (Ex(j,k+1,l) - Ex(j,k,l))
    end do
   end do
  end do

else
 if (l_1dz) then
  j=0
  k=0
  ! advance Bx
  do l = -nzs, nz+nzs-1
      Bx(j,k,l) = Bx(j,k,l) + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
  end do

  ! advance By
  do l = -nzs, nz+nzs-1
      By(j,k,l) = By(j,k,l) - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
  end do

 else if (.not. l_2drz) then
  k=0
  ! advance Bx
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs
      Bx(j,k,l) = Bx(j,k,l) + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
    end do
  end do

  ! advance By
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs-1
      By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &  
                            - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
    end do
  end do

  ! advance Bz 
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-1
      Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j,k,l)) 
    end do
  end do

 else ! l_2drz = True

  k=0
  ! advance Br
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs
      Bx(j,k,l) = Bx(j,k,l) + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
    end do
  end do

  ! advance Btheta
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs-1
      By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &  
                            - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
    end do
  end do

  ! advance Bz 
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-1
      ru = 1.+0.5/(xmin/dx+j+0.5)
      rd = 1.-0.5/(xmin/dx+j+0.5)
      Bz(j,k,l) = Bz(j,k,l) - dtsdx * (ru*Ey(j+1,k,l) - rd*Ey(j,k,l)) 
    end do
   end do

 end if
end if

return
end subroutine push_em3d_bvec

subroutine push_em3d_bvec_norder(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,dx,dy,dz, &
                          xmin,ymin,zmin,nx,ny,nz,nxguard,nyguard,nzguard, &
                          norderx,nordery,norderz, &
                          nxs,nys,nzs,l_1dz,l_2dxz,l_2drz,l_nodalgrid)
integer :: nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs,norderx,nordery,norderz
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN) :: dtsdx(norderx/2),dtsdy(norderx/2),dtsdz(norderx/2),xmin,ymin,zmin,dx,dy,dz
integer(ISZ) :: i,j,k,l,ist
logical(ISZ) :: l_1dz,l_2dxz,l_2drz,l_nodalgrid
real(kind=8) :: rd, ru

if (l_nodalgrid) then
  ist = 0
else
  ist = 1
end if

if (.not.l_2dxz) then

  ! advance Bx
  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs
      do i = 1, nordery/2
        Bx(j,k,l) = Bx(j,k,l) - dtsdy(i) * (Ez(j,k+i,l  ) - Ez(j,k-i+ist,l))
      end do
      do i = 1, norderz/2
        Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,  l+i) - Ey(j,k,l-i+ist))
      end do
    end do
   end do
  end do

  ! advance By
  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-ist
      do i = 1, norderx/2
        By(j,k,l) = By(j,k,l) + dtsdx(i) * (Ez(j+i,k,l  ) - Ez(j-i+ist,k,l)) 
      end do
      do i = 1, norderz/2
        By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j  ,k,l+i) - Ex(j,k,l-i+ist)) 
      end do
    end do
   end do
  end do

  ! advance Bz 
  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs-ist
      do i = 1, norderx/2
        Bz(j,k,l) = Bz(j,k,l) - dtsdx(i) * (Ey(j+i,k,l) - Ey(j-i+ist,k,l))
      end do
      do i = 1, nordery/2
        Bz(j,k,l) = Bz(j,k,l) + dtsdy(i) * (Ex(j,k+i,l) - Ex(j,k-i+ist,l))
      end do
    end do
   end do
  end do

else
 if (l_1dz) then
  j=0
  k=0
  ! advance Bx
  do l = -nzs, nz+nzs-ist
    do i = 1, norderz/2
      Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,  l+i) - Ey(j,k,l-i+ist))
    end do
  end do

  ! advance By
  do l = -nzs, nz+nzs-ist
    do i = 1, norderz/2
      By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j  ,k,l+i) - Ex(j,k,l-i+ist)) 
    end do
  end do

 else if (.not. l_2drz) then
  k=0
  ! advance Bx
  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs
      do i = 1, norderz/2
        Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,  l+i) - Ey(j,k,l-i+ist))
      end do
    end do
  end do

  ! advance By
  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs-ist
      do i = 1, norderx/2
        By(j,k,l) = By(j,k,l) + dtsdx(i) * (Ez(j+i,k,l  ) - Ez(j-i+ist,k,l))  
      end do
      do i = 1, norderz/2
        By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j  ,k,l+i) - Ex(j,k,l-i+ist)) 
      end do
    end do
  end do

  ! advance Bz 
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-ist
      do i = 1, norderx/2
        Bz(j,k,l) = Bz(j,k,l) - dtsdx(i) * (Ey(j+i,k,l) - Ey(j-i+ist,k,l)) 
      end do
    end do
  end do

 else ! l_2drz = True

  if (norderx.ne.2) then
    write(0,*) 'Error: norderx>2 not supported in RZ axisymmetric mode'
    call abort()
  end if

  k=0
  ! advance Br
  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs
      do i = 1, norderz/2
        Bx(j,k,l) = Bx(j,k,l) + dtsdz(i) * (Ey(j,k,  l+i) - Ey(j,k,l-i+ist))
      end do
    end do
  end do

  ! advance Btheta
  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs-1
      By(j,k,l) = By(j,k,l) + dtsdx(1) * (Ez(j+1,k,l  ) - Ez(j,k,l))   
      do i = 1, norderz/2
        By(j,k,l) = By(j,k,l) - dtsdz(i) * (Ex(j  ,k,l+i) - Ex(j,k,l-i+ist)) 
      end do
    end do
  end do

  ! advance Bz 
  do l = -nzs, nz+nzs
    do j = -nxs, nx+nxs-1
      ru = 1.+0.5/(xmin/dx+j+0.5)
      rd = 1.-0.5/(xmin/dx+j+0.5)
      Bz(j,k,l) = Bz(j,k,l) - dtsdx(1) * (ru*Ey(j+1,k,l) - rd*Ey(j,k,l)) 
    end do
   end do

 end if
end if

return
end subroutine push_em3d_bvec_norder

subroutine push_em3d_bvec_circ(ex,ey,ez,bx,by,bz,dtsdx,dtsdz,dx,dz, &
                          xmin,zmin,nx,nz,nxguard,nzguard,circ_m)
integer :: nx,nz,nxguard,nzguard,circ_m
complex(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,circ_m) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN) :: dtsdx,dtsdz,xmin,zmin,dx,dz
integer(ISZ) :: j,l,m
real(kind=8) :: rd, ru, r, dt
complex(kind=8) :: i=(0.,1.)

  dt = dtsdx*dx
  do m = 1, circ_m

     ! advance Bx
     do l = 0, nz-1
        do j = 0,nx
           if (j==0 .and. xmin==0) then
              ! On axis
              if (.not. m == 1) then
                 ! Br should remain 0 on axis, for modes different than m=1,
                 ! but the bulk equation does not necessarily ensure this.
                 Bx(j,l,m) = 0.
              else
                 ! For the mode m = 1, the bulk equation diverges on axis
                 ! (due to the 1/r terms). The following expressions regularize
                 ! these divergences by assuming, on axis :
                 ! Ez/r = 0/r + dEz/dr
                 Bx(j,l,m) = Bx(j,l,m) + i*m*dt*Ez(j+1,l,m)/dx &
                      + dtsdz * (Ey(j,  l+1,m) - Ey(j,l,m))
              endif
           else
              ! Equations in the bulk of the grid
              r = xmin+j*dx
              Bx(j,l,m) = Bx(j,l,m) + i*m*dt*Ez(j,l,m)/r &
                   + dtsdz * (Ey(j,  l+1,m) - Ey(j,l,m))
           endif
        end do
     end do

     ! advance Btheta
     do l = 0, nz-1
        do j = 0, nx-1
           By(j,l,m) = By(j,l,m) + dtsdx * (Ez(j+1,l  ,m) - Ez(j,l,m)) &  
                - dtsdz * (Ex(j,l+1,m) - Ex(j,l,m)) 
        end do
     end do

     ! advance Bz 
     do l = 0, nz
        do j = 0, nx-1
           r  = xmin+j*dx+0.5*dx
           ru = 1.+0.5/(xmin/dx+j+0.5)
           rd = 1.-0.5/(xmin/dx+j+0.5)
           Bz(j,l,m) = Bz(j,l,m) - dtsdx * (ru*Ey(j+1,l,m) - rd*Ey(j,l,m)) &
                - i*m*dt*Ex(j,l,m)/r
        end do
     end do
     
  end do
return
end subroutine push_em3d_bvec_circ

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
                            - betayx*dtsdy * (Ez(j+1,k+1,l  ) - Ez(j+1,k  ,l  ) &
                                           +  Ez(j-1,k+1,l  ) - Ez(j-1,k  ,l  )) &
                            - betayz*dtsdy * (Ez(j  ,k+1,l+1) - Ez(j  ,k  ,l+1) &
                                           +  Ez(j  ,k+1,l-1) - Ez(j  ,k  ,l-1)) &
                            - gammay*dtsdy * (Ez(j+1,k+1,l+1) - Ez(j+1,k  ,l+1) &
                                           +  Ez(j-1,k+1,l+1) - Ez(j-1,k  ,l+1) &
                                           +  Ez(j+1,k+1,l-1) - Ez(j+1,k  ,l-1) &
                                           +  Ez(j-1,k+1,l-1) - Ez(j-1,k  ,l-1)) &
                            + alphaz*dtsdz * (Ey(j  ,k  ,l+1) - Ey(j  ,k  ,l  )) &
                            + betazx*dtsdz * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  ) &
                                           +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  )) &
                            + betazy*dtsdz * (Ey(j  ,k+1,l+1) - Ey(j  ,k+1,l  ) &
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
                            + betaxy*dtsdx * (Ez(j+1,k+1,l  ) - Ez(j  ,k+1,l  ) &
                                           +  Ez(j+1,k-1,l  ) - Ez(j  ,k-1,l  )) &
                            + betaxz*dtsdx * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1) &
                                           +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                            + gammax*dtsdx * (Ez(j+1,k+1,l+1) - Ez(j  ,k+1,l+1) &
                                           +  Ez(j+1,k-1,l+1) - Ez(j  ,k-1,l+1) &
                                           +  Ez(j+1,k+1,l-1) - Ez(j  ,k+1,l-1) &
                                           +  Ez(j+1,k-1,l-1) - Ez(j  ,k-1,l-1)) &
                            - alphaz*dtsdz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                            - betazx*dtsdz * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  ) &
                                           +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  )) &
                            - betazy*dtsdz * (Ex(j  ,k+1,l+1) - Ex(j  ,k+1,l  ) &
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
                            - betaxy*dtsdx * (Ey(j+1,k+1,l  ) - Ey(j  ,k+1,l  ) &
                                           +  Ey(j+1,k-1,l  ) - Ey(j  ,k-1,l  )) &
                            - betaxz*dtsdx * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1) &
                                           +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) &
                            - gammax*dtsdx * (Ey(j+1,k+1,l+1) - Ey(j  ,k+1,l+1) &
                                           +  Ey(j+1,k-1,l+1) - Ey(j  ,k-1,l+1) &
                                           +  Ey(j+1,k+1,l-1) - Ey(j  ,k+1,l-1) &
                                           +  Ey(j+1,k-1,l-1) - Ey(j  ,k-1,l-1)) &
                            + alphay*dtsdy * (Ex(j  ,k+1,l  ) - Ex(j  ,k  ,l  )) &
                            + betayx*dtsdy * (Ex(j+1,k+1,l  ) - Ex(j+1,k  ,l  ) &
                                           +  Ex(j-1,k+1,l  ) - Ex(j-1,k  ,l  )) &
                            + betayz*dtsdy * (Ex(j  ,k+1,l+1) - Ex(j  ,k  ,l+1) &
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
      Bx(j,k,l) = Bx(j,k,l) +    alphaz*dtsdz * (Ey(j  ,k  ,l+1) - Ey(j  ,k  ,l  )) &
                            +    betazx*dtsdz * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  ) &
                                              +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  )) 
    end do
  end do

  ! advance By
  do l = 0, nz-1
    do j = 0, nx-1
      By(j,k,l) = By(j,k,l) +    alphax*dtsdx * (Ez(j+1,k  ,l  ) - Ez(j  ,k  ,l  )) &  
                            +    betaxz*dtsdx * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1) &
                                              +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                            -    alphaz*dtsdz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                            -    betazx*dtsdz * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  ) &
                                              +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  )) 
    end do
  end do

  ! advance Bz 
  do l = 0, nz
    do j = 0, nx-1
      Bz(j,k,l) = Bz(j,k,l) -    alphax*dtsdx * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  )) &
                            -    betaxz*dtsdx * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1) &
                                              +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) 
    end do
  end do

end if

return
end subroutine push_em3d_kyeebvec

subroutine push_em3d_lehebvec(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
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
              Bx(j,k,l) = Bx(j,k,l) &
                   - alphay*dtsdy * (Ez(j  ,k+1,l  ) - Ez(j  ,k  ,l  )) &
                   - betayx*dtsdy * (Ez(j+1,k+1,l  ) - Ez(j+1,k  ,l  ) &
                   +  Ez(j-1,k+1,l  ) - Ez(j-1,k  ,l  )) &
                   - betayz*dtsdy * (Ez(j  ,k+1,l+1) - Ez(j  ,k  ,l+1) &
                   +  Ez(j  ,k+1,l-1) - Ez(j  ,k  ,l-1)) &
                   + alphaz*dtsdz * (Ey(j  ,k  ,l+1) - Ey(j  ,k  ,l  )) &
                   + betazx*dtsdz * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  ) &
                   +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  )) &
                   + betazy*dtsdz * (Ey(j  ,k+1,l+1) - Ey(j  ,k+1,l  ) &
                   +  Ey(j  ,k-1,l+1) - Ey(j  ,k-1,l  )) &
                   + deltaz*dtsdz * (Ey(j  ,k  ,l+2) - Ey(j  ,k  ,l-1))
           end do
        end do
     end do

     ! advance By
     do l = 0, nz-1
        do k = 0, ny
           do j = 0, nx-1
              By(j,k,l) = By(j,k,l) &
                   + alphax*dtsdx * (Ez(j+1,k  ,l  ) - Ez(j  ,k  ,l  )) &  
                   + betaxy*dtsdx * (Ez(j+1,k+1,l  ) - Ez(j  ,k+1,l  ) &
                   +  Ez(j+1,k-1,l  ) - Ez(j  ,k-1,l  )) &
                   + betaxz*dtsdx * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1) &
                   +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                   - alphaz*dtsdz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                   - betazx*dtsdz * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  ) &
                   +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  )) &
                   - betazy*dtsdz * (Ex(j  ,k+1,l+1) - Ex(j  ,k+1,l  ) &
                   +  Ex(j  ,k-1,l+1) - Ex(j  ,k-1,l  )) &
                   - deltaz*dtsdz * (Ex(j  ,k  ,l+2) - Ex(j  ,k  ,l-1)) 
           end do
        end do
     end do

     ! advance Bz 
     do l = 0, nz
        do k = 0, ny-1
           do j = 0, nx-1
              Bz(j,k,l) = Bz(j,k,l) &
                   - alphax*dtsdx * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  )) &
                   - betaxy*dtsdx * (Ey(j+1,k+1,l  ) - Ey(j  ,k+1,l  ) &
                   +  Ey(j+1,k-1,l  ) - Ey(j  ,k-1,l  )) &
                   - betaxz*dtsdx * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1) &
                   +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) &
                   + alphay*dtsdy * (Ex(j  ,k+1,l  ) - Ex(j  ,k  ,l  )) &
                   + betayx*dtsdy * (Ex(j+1,k+1,l  ) - Ex(j+1,k  ,l  ) &
                   +  Ex(j-1,k+1,l  ) - Ex(j-1,k  ,l  )) &
                   + betayz*dtsdy * (Ex(j  ,k+1,l+1) - Ex(j  ,k  ,l+1) &
                   +  Ex(j  ,k+1,l-1) - Ex(j  ,k  ,l-1)) 
           end do
        end do
     end do

  else

     k=0
     ! advance Bx
     do l = 0, nz-1
        do j = 0, nx
           Bx(j,k,l) = Bx(j,k,l) &
                +    alphaz*dtsdz * (Ey(j  ,k  ,l+1) - Ey(j  ,k  ,l  )) &
                +    betazx*dtsdz * (Ey(j+1,k  ,l+1) - Ey(j+1,k  ,l  ) &
                +  Ey(j-1,k  ,l+1) - Ey(j-1,k  ,l  )) &
                +    deltaz*dtsdz * (Ey(j  ,k  ,l+2) - Ey(j  ,k  ,l-1))
        end do
     end do

     ! advance By
     do l = 0, nz-1
        do j = 0, nx-1
           By(j,k,l) = By(j,k,l) +    alphax*dtsdx * (Ez(j+1,k  ,l  ) - Ez(j  ,k  ,l  )) &  
                +    betaxz*dtsdx * (Ez(j+1,k  ,l+1) - Ez(j  ,k  ,l+1) &
                +  Ez(j+1,k  ,l-1) - Ez(j  ,k  ,l-1)) &
                -    alphaz*dtsdz * (Ex(j  ,k  ,l+1) - Ex(j  ,k  ,l  )) &
                -    betazx*dtsdz * (Ex(j+1,k  ,l+1) - Ex(j+1,k  ,l  ) &
                +  Ex(j-1,k  ,l+1) - Ex(j-1,k  ,l  )) &
                -    deltaz*dtsdz * (Ex(j  ,k  ,l+2) - Ex(j  ,k  ,l-1)) 
        end do
     end do

     ! advance Bz 
     do l = 0, nz
        do j = 0, nx-1
           Bz(j,k,l) = Bz(j,k,l) -    alphax*dtsdx * (Ey(j+1,k  ,l  ) - Ey(j  ,k  ,l  )) &
                -    betaxz*dtsdx * (Ey(j+1,k  ,l+1) - Ey(j  ,k  ,l+1) &
                +  Ey(j+1,k  ,l-1) - Ey(j  ,k  ,l-1)) 
        end do
     end do

  end if

  return
end subroutine push_em3d_lehebvec

subroutine push_em3d_f(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,dtsepsi

!if (f%spectral) return

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz
dtsepsi = f%mu0*f%clight**3*dt

select case (f%stencil)
     ! Choose the kind of stencil that is used for the E correction
     ! (Second step of the propagative Poisson correction)

case( 0,2 ) ! Yee stencil
 if(f%nconds>0) then 
  call push_em3d_fvec_cond(f%ex,f%ey,f%ez,f%f, f%rho, &
                      dtsepsi,dtsdx,dtsdy,dtsdz, &
                      f%dx,f%dy,f%dz, &
                      f%nx,f%ny,f%nz, &
                      f%xmin,f%ymin,f%zmin, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz,f%l_2drz,f%incond)
 else
   call push_em3d_fvec(f%ex,f%ey,f%ez,f%f, f%rho, &
                      dtsepsi,dtsdx,dtsdy,dtsdz, &
                      f%dx,f%dy,f%dz, &
                      f%nx,f%ny,f%nz, &
                      f%xmin,f%ymin,f%zmin, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz,f%l_2drz)
   if (f%circ_m>0) &
     call push_em3d_fvec_circ(f%ex_circ,f%ey_circ,f%ez_circ,f%f_circ, f%rho_circ, &
                      dtsepsi,dtsdx,dtsdz, &
                      f%dx,f%dz, &
                      f%nx,f%nz, &
                      f%xmin,f%zmin, &
                      f%nxguard,f%nzguard,f%circ_m)
 end if

case( 1 ) ! Cole-Karkkainen stencil (since the B push uses the Cole-Karkkainen stencil,
   ! and since this Poisson correction should not modify curl(E) in the B push)
  call push_em3d_kyeefvec(f%ex,f%ey,f%ez,f%f, f%rho, &
                      dtsepsi,dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)

case( 3 ) ! Lehe stencil (since the B push uses Lehe stencil,
     ! and since this Poisson correction should not modify curl(E) in the B push)
     call push_em3d_leheefvec(f%ex,f%ey,f%ez,f%f, &
          dtsdx,dtsdy,dtsdz, &
          f%nx,f%ny,f%nz, &
          f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
  
end select

end subroutine push_em3d_f

subroutine push_em3d_fvec(ex,ey,ez,f,rho,dtsepsi,dtsdx,dtsdy,dtsdz,dx,dy,dz,nx,ny,nz, &
                          xmin,ymin,zmin,nxguard,nyguard,nzguard,l_2dxz,l_2drz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f,rho
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz,dtsepsi,xmin,ymin,zmin,dx,dy,dz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz,l_2drz
real(kind=8) :: ru,rd

if (.not.l_2dxz) then
  ! --- 3D XYZ
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
 if (.not.l_2drz) then
  ! --- 2D XZ
  k=0
  do l = 0, nz
    do j = 0, nx
      F(j,k,l) = F(j,k,l) + dtsdx * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                          + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                          - dtsepsi * Rho(j,k,l)
    end do
  end do
 else
  ! --- 2D RZ (axisymmetric)
  k=0
  do l = 0, nz
    do j = 0,nx
       if (j==0 .and. xmin==0.) then
          ! the bulk equation diverges on axis
          ! (due to the 1/r terms). The following expressions regularize
          ! these divergences.
          F(j,k,l) = F(j,k,l) + 4.*dtsdx * Ex(j,k,l) &
               + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
               - dtsepsi * Rho(j,k,l)
       else
          ru = 1.+0.5/(xmin/dx+j)
          rd = 1.-0.5/(xmin/dx+j)
          F(j,k,l) = F(j,k,l) + dtsdx * (ru*Ex(j,k,l) - rd*Ex(j-1,k,l)) &
               + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
               - dtsepsi * Rho(j,k,l)
       end if
       
    end do
 end do
end if
end if

return
end subroutine push_em3d_fvec

subroutine push_em3d_fvec_circ(ex,ey,ez,f,rho,dtsepsi,dtsdx,dtsdz,dx,dz,nx,nz, &
                          xmin,zmin,nxguard,nzguard,circ_m)
integer :: nx,ny,nz,nxguard,nyguard,nzguard,circ_m
complex(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,circ_m) :: ex,ey,ez,f,rho
real(kind=8), intent(IN) :: dtsdx,dtsdz,dtsepsi,xmin,zmin,dx,dz
integer(ISZ) :: j,l,m
real(kind=8) :: ru,rd,r,dt
complex(kind=8) :: i=(0.,1.)

  dt=dtsdx*dx
  do m = 1, circ_m

     ! --- 2D RZ (axisymmetric)
     do l = 0, nz
        do j = 0, nx
           if (j==0 .and. xmin==0) then
              ! F should remain 0 on axis, for modes different than m=0,
              ! but the bulk equation does not necessarily ensure this.
              F(j,l,m) = 0.
           else
              ! Equations for the bulk of the grid
              ru = 1.+0.5/(xmin/dx+j)
              rd = 1.-0.5/(xmin/dx+j)
              r = xmin+j*dx
              F(j,l,m) = F(j,l,m) + dtsdx * (ru*Ex(j,l,m) - rd*Ex(j-1,l  , m)) &
                   - i*m*dt*Ey(j,l,m)/r &
                   + dtsdz * (Ez(j,l,m) - Ez(j  ,l-1, m)) &
                   - dtsepsi * Rho(j,l,m)
           endif
        end do

     end do
  end do

return
end subroutine push_em3d_fvec_circ

subroutine push_em3d_fvec_cond(ex,ey,ez,f,rho,dtsepsi,dtsdx,dtsdy,dtsdz,dx,dy,dz,nx,ny,nz, &
                          xmin,ymin,zmin,nxguard,nyguard,nzguard,l_2dxz,l_2drz,incond)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f,rho
logical(ISZ), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: incond
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz,dtsepsi,xmin,ymin,zmin,dx,dy,dz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz,l_2drz
real(kind=8) :: ru,rd

if (.not.l_2dxz) then
  ! --- 3D XYZ
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
      if (.not.incond(j,k,l)) &
        F(j,k,l) = F(j,k,l) + dtsdx * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                            + dtsdy * (Ey(j,k,l) - Ey(j  ,k-1,l  )) &
                            + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                            - dtsepsi * Rho(j,k,l)
    end do
   end do
  end do

else
 if (.not.l_2drz) then
  ! --- 2D XZ
  k=0
  do l = 0, nz
    do j = 0, nx
      if (.not.incond(j,k,l)) &
        F(j,k,l) = F(j,k,l) + dtsdx * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                            + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                            - dtsepsi * Rho(j,k,l)
    end do
  end do
 else
  ! --- 2D RZ (axisymmetric)
  k=0
  do l = 0, nz
    j = 0
    if (xmin==0.) then
      if (.not.incond(j,k,l)) & 
        F(j,k,l) = F(j,k,l) + 4.*dtsdx * Ex(j,k,l) &
                            + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                            - dtsepsi * Rho(j,k,l)
    else
      ru = 1.+0.5/(xmin/dx)
      rd = 1.-0.5/(xmin/dx)
      if (.not.incond(j,k,l)) &
        F(j,k,l) = F(j,k,l) + dtsdx * (ru*Ex(j,k,l) - rd*Ex(j-1,k  ,l  )) &
                            + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                            - dtsepsi * Rho(j,k,l)
    end if
    do j = 1, nx
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      if (.not.incond(j,k,l)) & 
        F(j,k,l) = F(j,k,l) + dtsdx * (ru*Ex(j,k,l) - rd*Ex(j-1,k  ,l  )) &
                            + dtsdz * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) &
                            - dtsepsi * Rho(j,k,l)
    end do
  end do
 end if
end if

return
end subroutine push_em3d_fvec_cond

subroutine getdive(ex,ey,ez,dive,dx,dy,dz,nx,ny,nz,nxguard,nyguard,nzguard,xmin,l_2dxz,l_2drz)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,dive
real(kind=8), intent(IN) :: dx,dy,dz,xmin
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz,l_2drz
real(kind=8) :: ru,rd,dxi,dyi,dzi

dxi = 1./dx
dyi = 1./dy
dzi = 1./dz
if (.not.l_2dxz) then
  ! --- 3D XYZ
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
      dive(j,k,l) = dive(j,k,l) + dxi * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                                + dyi * (Ey(j,k,l) - Ey(j  ,k-1,l  )) &
                                + dzi * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) 
    end do
   end do
  end do

else
 if (.not.l_2drz) then
  ! --- 2D XZ
  k=0
  do l = 0, nz
    do j = 0, nx
      dive(j,k,l) = dive(j,k,l) + dxi * (Ex(j,k,l) - Ex(j-1,k  ,l  )) &
                                + dzi * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) 
    end do
  end do
 else
  ! --- 2D RZ (axisymmetric)
  k=0
  do l = 0, nz
    j = 0
    if (xmin==0.) then
      dive(j,k,l) = dive(j,k,l) + 4.*dxi * Ex(j,k,l) &
                                + dzi * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) 
    else
      ru = 1.+0.5/(xmin/dx)
      rd = 1.-0.5/(xmin/dx)
      dive(j,k,l) = dive(j,k,l) + dxi * (ru*Ex(j,k,l) - rd*Ex(j-1,k  ,l  )) &
                                + dzi * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) 
    end if
    do j = 1, nx
      ru = 1.+0.5/(xmin/dx+j)
      rd = 1.-0.5/(xmin/dx+j)
      dive(j,k,l) = dive(j,k,l) + dxi * (ru*Ex(j,k,l) - rd*Ex(j-1,k  ,l  )) &
                                + dzi * (Ez(j,k,l) - Ez(j  ,k  ,l-1)) 
    end do
  end do
 end if
end if

return
end subroutine getdive

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
                          + betaxy*dtsdx * (Ex(j  ,k+1,l  ) - Ex(j-1,k+1,l  ) &
                                         +  Ex(j  ,k-1,l  ) - Ex(j-1,k-1,l  )) &
                          + betaxz*dtsdx * (Ex(j  ,k  ,l+1) - Ex(j-1,k  ,l+1) &
                                         +  Ex(j  ,k  ,l-1) - Ex(j-1,k  ,l-1)) &
                          + gammax*dtsdx * (Ex(j  ,k+1,l+1) - Ex(j-1,k+1,l+1) &
                                         +  Ex(j  ,k-1,l+1) - Ex(j-1,k-1,l+1) &
                                         +  Ex(j  ,k+1,l-1) - Ex(j-1,k+1,l-1) &
                                         +  Ex(j  ,k-1,l-1) - Ex(j-1,k-1,l-1)) &
                          + alphay*dtsdy * (Ey(j  ,k  ,l  ) - Ey(j  ,k-1,l  )) &
                          + betayx*dtsdy * (Ey(j+1,k  ,l  ) - Ey(j+1,k-1,l  ) &
                                         +  Ey(j-1,k  ,l  ) - Ey(j-1,k-1,l  )) &
                          + betayz*dtsdy * (Ey(j  ,k  ,l+1) - Ey(j  ,k-1,l+1) &
                                         +  Ey(j  ,k  ,l-1) - Ey(j  ,k-1,l-1)) &
                          + gammay*dtsdy * (Ey(j+1,k  ,l+1) - Ey(j+1,k-1,l+1) &
                                         +  Ey(j-1,k  ,l+1) - Ey(j-1,k-1,l+1) &
                                         +  Ey(j+1,k  ,l-1) - Ey(j+1,k-1,l-1) &
                                         +  Ey(j-1,k  ,l-1) - Ey(j-1,k-1,l-1)) &
                          + alphaz*dtsdz * (Ez(j  ,k  ,l  ) - Ez(j  ,k  ,l-1)) &
                          + betazx*dtsdz * (Ez(j+1,k  ,l  ) - Ez(j+1,k  ,l-1) &
                                         +  Ez(j-1,k  ,l  ) - Ez(j-1,k  ,l-1)) &
                          + betazy*dtsdz * (Ez(j  ,k+1,l  ) - Ez(j  ,k+1,l-1) &
                                         +  Ez(j  ,k-1,l  ) - Ez(j  ,k-1,l-1)) &
                          + gammaz*dtsdz * (Ez(j+1,k+1,l  ) - Ez(j+1,k+1,l-1) &
                                         +  Ez(j-1,k+1,l  ) - Ez(j-1,k+1,l-1) &
                                         +  Ez(j+1,k-1,l  ) - Ez(j+1,k-1,l-1) &
                                         +  Ez(j-1,k-1,l  ) - Ez(j-1,k-1,l-1)) &

                          - dtsepsi/3. * ( (alphax+alphay+alphaz)* Rho(j,k,l) &
                                        +  betaxy*(Rho(j  ,k+1,l  )+Rho(j  ,k-1,l  )) &
                                        +  betaxz*(Rho(j  ,k  ,l+1)+Rho(j  ,k  ,l-1)) &
                                        +  gammax*(Rho(j  ,k+1,l+1)+Rho(j  ,k-1,l+1)+Rho(j  ,k+1,l-1)+Rho(j  ,k-1,l-1)) &
                                        +  betayx*(Rho(j+1,k  ,l  )+Rho(j-1,k  ,l  )) &
                                        +  betayz*(Rho(j  ,k  ,l+1)+Rho(j  ,k  ,l-1)) &
                                        +  gammay*(Rho(j+1,k  ,l+1)+Rho(j-1,k  ,l+1)+Rho(j+1,k  ,l-1)+Rho(j-1,k  ,l-1)) &
                                        +  betazx*(Rho(j+1,k  ,l  )+Rho(j-1,k  ,l  )) &
                                        +  betazy*(Rho(j  ,k+1,l  )+Rho(j  ,k-1,l  )) &
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
                          +      betaxz*dtsdx * (Ex(j  ,k  ,l+1) - Ex(j-1,k  ,l+1) &
                                         +  Ex(j  ,k  ,l-1) - Ex(j-1,k  ,l-1)) &
                          + alphaz*dtsdz * (Ez(j  ,k  ,l  ) - Ez(j  ,k  ,l-1)) &
                          +      betazx*dtsdz * (Ez(j+1,k  ,l  ) - Ez(j+1,k  ,l-1) &
                                         +  Ez(j-1,k  ,l  ) - Ez(j-1,k  ,l-1)) &
                          - dtsepsi/2. * ( (alphax+alphaz)* Rho(j,k,l) &
                                        +       betaxz*(Rho(j  ,k  ,l+1)+Rho(j  ,k  ,l-1)) &
                                        +       betazx*(Rho(j+1,k  ,l  )+Rho(j-1,k  ,l  )) )
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

!if (f%spectral) return

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz

if (f%stencil==0 .or. f%stencil==2) then
  if(f%nconds>0) then 
    call push_em3d_efvec_cond(f%ex,f%ey,f%ez,f%f, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz,f%incond)
  else
    call push_em3d_efvec(f%ex,f%ey,f%ez,f%f, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard,f%l_2dxz)
    if (f%circ_m>0) &
      call push_em3d_efvec_circ(f%ex_circ,f%ey_circ,f%ez_circ,f%f_circ, &
                        dtsdx,dtsdz, &
                        f%nx,f%nz, &
                        f%nxguard,f%nzguard,f%xmin,f%dx,f%circ_m)

  end if
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

subroutine push_em3d_efvec_circ(ex,ey,ez,f,dtsdx,dtsdz,nx,nz,nxguard,nzguard,xmin,dx,circ_m)
integer :: nx,nz,nxguard,nzguard,circ_m
complex(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nzguard:nz+nzguard,circ_m) :: ex,ey,ez,f
real(kind=8), intent(IN) :: dtsdx,dtsdz,xmin,dx
integer(ISZ) :: j,l,m,jmin
complex(kind=8) :: i=(0.,1.)
real(kind=8) :: dt,r

  ! ===============================
  !             2-D RZ multipole
  ! ===============================

  dt = dtsdx*dx
  do m = 1, circ_m

     ! advance Ex
     do l = 0, nz
        do j = 0, nx-1
           Ex(j,l,m) = Ex(j,l,m) + dtsdx * (F(j+1,l,m) - F(j,l,m)) 
        end do
     end do

     ! advance Ey
     do l = 0, nz
        do j = 0, nx
           if (j==0 .and. xmin==0.) then
              ! On axis, the bulk equations diverge (due to
              ! the 1/r terms). The following expression
              ! regularizes this divergence by assuming
              ! F/r = 0/r + dF/dr on axis
              Ey(j,l,m) = Ey(j,l,m) - i*m*dt*( F(j+1,l,m) - F(j,l,m) )/dx
           else
              ! Equation for the bulk of the grid
              r = xmin+j*dx
              Ey(j,l,m) = Ey(j,l,m) - i*m*dt*F(j,l,m)/r
           end if
        end do
     end do

     ! advance Ez 
     do l = 0, nz-1
        do j = 0, nx
           Ez(j,l,m) = Ez(j,l,m) + dtsdz * (F(j,l+1,m) - F(j,l,m)) 
        end do
     end do

  end do

return
end subroutine push_em3d_efvec_circ

subroutine push_em3d_efvec_cond(ex,ey,ez,f,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz,incond)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f
logical(ISZ), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: incond
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      if (.not.incond(j,k,l) .or. .not.incond(j+1,k,l)) &
      Ex(j,k,l) = Ex(j,k,l) + dtsdx * (F(j+1,k,l) - F(j,k,l)) 
    end do
   end do
  end do

  ! advance Ey
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      if (.not.incond(j,k,l) .or. .not.incond(j,k+1,l)) &
      Ey(j,k,l) = Ey(j,k,l) + dtsdy * (F(j,k+1,l) - F(j,k,l))
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + dtsdz * (F(j,k,l+1) - F(j,k,l)) 
    end do
   end do
  end do

else

  k=0
  ! advance Ex
  do l = 0, nz
    do j = 0, nx-1
      if (.not.incond(j,k,l) .or. .not.incond(j+1,k,l)) &
      Ex(j,k,l) = Ex(j,k,l) + dtsdx * (F(j+1,k,l) - F(j,k,l)) 
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      if (.not.incond(j,k,l) .or. .not.incond(j,k,l+1)) &
      Ez(j,k,l) = Ez(j,k,l) + dtsdz * (F(j,k,l+1) - F(j,k,l)) 
    end do
  end do

end if

return
end subroutine push_em3d_efvec_cond

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
                            + betaxy*dtsdx * (F(j+1,k+1,l  ) - F(j  ,k+1,l  ) &
                                           +  F(j+1,k-1,l  ) - F(j  ,k-1,l  )) &
                            + betaxz*dtsdx * (F(j+1,k  ,l+1) - F(j  ,k  ,l+1) &
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
                            + betayx*dtsdy * (F(j+1,k+1,l  ) - F(j+1,k  ,l  ) &
                                           +  F(j-1,k+1,l  ) - F(j-1,k  ,l  )) &
                            + betayz*dtsdy * (F(j  ,k+1,l+1) - F(j  ,k  ,l+1) &
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
                            + betazx*dtsdz * (F(j+1,k  ,l+1) - F(j+1,k  ,l  ) &
                                           +  F(j-1,k  ,l+1) - F(j-1,k  ,l  )) &
                            + betazy*dtsdz * (F(j  ,k+1,l+1) - F(j  ,k+1,l  ) &
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
                            + betaxz*dtsdx * (F(j+1,k  ,l+1) - F(j  ,k  ,l+1) &
                                           +  F(j+1,k  ,l-1) - F(j  ,k  ,l-1))  
    end do
  end do

  ! advance Ez 
  do l = 0, nz-1
    do j = 0, nx
      Ez(j,k,l) = Ez(j,k,l) + alphaz*dtsdz * (F(j  ,k  ,l+1) - F(j  ,k  ,l  )) &
                            + betazx*dtsdz * (F(j+1,k  ,l+1) - F(j+1,k  ,l  ) &
                                           +  F(j-1,k  ,l+1) - F(j-1,k  ,l  ))
    end do
  end do

end if

return
end subroutine push_em3d_kyeeefvec

subroutine push_em3d_leheefvec(ex,ey,ez,f,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard,l_2dxz)
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
              Ex(j,k,l) = Ex(j,k,l) &
                   + alphax*dtsdx * (F(j+1,k  ,l  ) - F(j  ,k  ,l  )) &
                   + betaxy*dtsdx * (F(j+1,k+1,l  ) - F(j  ,k+1,l  ) &
                   +  F(j+1,k-1,l  ) - F(j  ,k-1,l  )) &
                   + betaxz*dtsdx * (F(j+1,k  ,l+1) - F(j  ,k  ,l+1) &
                   +  F(j+1,k  ,l-1) - F(j  ,k  ,l-1)) 
           end do
        end do
     end do

     ! advance Ey
     do l = 0, nz
        do k = 0, ny-1
           do j = 0, nx
              Ey(j,k,l) = Ey(j,k,l) &
                   + alphay*dtsdy * (F(j  ,k+1,l  ) - F(j  ,k  ,l  )) &
                   + betayx*dtsdy * (F(j+1,k+1,l  ) - F(j+1,k  ,l  ) &
                   +  F(j-1,k+1,l  ) - F(j-1,k  ,l  )) &
                   + betayz*dtsdy * (F(j  ,k+1,l+1) - F(j  ,k  ,l+1) &
                   +  F(j  ,k+1,l-1) - F(j  ,k  ,l-1)) 
           end do
        end do
     end do

     ! advance Ez 
     do l = 0, nz-1
        do k = 0, ny
           do j = 0, nx
              Ez(j,k,l) = Ez(j,k,l) &
                   + alphaz*dtsdz * (F(j  ,k  ,l+1) - F(j  ,k  ,l  )) &
                   + betazx*dtsdz * (F(j+1,k  ,l+1) - F(j+1,k  ,l  ) &
                   +  F(j-1,k  ,l+1) - F(j-1,k  ,l  )) &
                   + betazy*dtsdz * (F(j  ,k+1,l+1) - F(j  ,k+1,l  ) &
                   +  F(j  ,k-1,l+1) - F(j  ,k-1,l  )) &
                   + deltaz*dtsdz * (F(j  ,k  ,l+2) - F(j  ,k  ,l-1)) 
           end do
        end do
     end do

  else

     k=0
     ! advance Ex
     do l = 0, nz
        do j = 0, nx-1
           Ex(j,k,l) = Ex(j,k,l) &
                + alphax*dtsdx * (F(j+1,k  ,l  ) - F(j  ,k  ,l  )) &
                + betaxz*dtsdx * (F(j+1,k  ,l+1) - F(j  ,k  ,l+1) &
                +  F(j+1,k  ,l-1) - F(j  ,k  ,l-1))  
        end do
     end do

     ! advance Ez 
     do l = 0, nz-1
        do j = 0, nx
           Ez(j,k,l) = Ez(j,k,l) &
                + alphaz*dtsdz * (F(j  ,k  ,l+1) - F(j  ,k  ,l  )) &
                + betazx*dtsdz * (F(j+1,k  ,l+1) - F(j+1,k  ,l  ) &
                +  F(j-1,k  ,l+1) - F(j-1,k  ,l  )) &
                + deltaz*dtsdz * (F(j  ,k  ,l+2) - F(j  ,k  ,l-1)) 
        end do
     end do

  end if

  return
end subroutine push_em3d_leheefvec

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

subroutine push_em3d_splite(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  if (sf%stencil==0 .or. sf%stencil==1) then
   if ((sf%norderx==2) .and. (sf%nordery==2) .and. (sf%norderz==2) .and. .not. sf%l_nodalgrid) then
    call push_em3d_splitevec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%exy,sf%exz,sf%eyx,sf%eyz,sf%ezx,sf%ezy, &
                             sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                             sf%afx,sf%afy,sf%afz, &
                             sf%bpfx,sf%bpfy,sf%bpfz, &
                             sf%bmfx,sf%bmfy,sf%bmfz,sf%l_1dz,sf%l_2dxz,sf%l_2drz, &
                             sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
   else
    call push_em3d_splitevec_norder(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%norderx,sf%nordery,sf%norderz, &
                             sf%xcoefs,sf%ycoefs,sf%zcoefs, &
                             sf%exy,sf%exz,sf%eyx,sf%eyz,sf%ezx,sf%ezy, &
                             sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                             sf%afx,sf%afy,sf%afz, &
                             sf%bpfx,sf%bpfy,sf%bpfz, &
                             sf%bmfx,sf%bmfy,sf%bmfz,sf%l_1dz,sf%l_2dxz,sf%l_2drz,sf%l_nodalgrid, &
                             sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
   end if
   if (sf%nconds>0) then 
       call push_em3d_splite_setcond(sf%nx,sf%ny,sf%nz,sf%nxcond,sf%nycond,sf%nzcond,sf%nxguard,sf%nyguard,sf%nzguard, &
                                     sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz,sf%incond,sf%l_2dxz,sf%l_2drz, &
                                     sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
   end if
  else
    write(0,*) 'splite extended pml not implemented'
    stop
  end if
  
  return
end subroutine push_em3d_splite

subroutine push_em3d_splitevec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exy,exz,eyx,eyz,ezx,ezy,bxy,byx,bzx,bxz,byz,bzy, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz,l_1dz,l_2dxz,l_2drz, &
                               xmin,ymin,zmin,dx,dy,dz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exy,exz, &
                                                                                                       eyx,eyz, &
                                                                                                       ezx,ezy
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz
real(kind=8), intent(in) :: xmin,ymin,zmin,dx,dy,dz

INTEGER :: j, k, l
logical(ISZ) :: l_1dz,l_2dxz,l_2drz
real(8) :: ru,rd

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

 if (l_1dz) then
  j = 0

  do l = 0, nz
      exz(j,k,l) = afz(l)*exz(j,k,l) - bpfz(l)*(byx(j,k,l)+byz(j,k,l))  &
                                     - bmfz(l)*(byx(j,k,l-1)+byz(j,k,l-1)) !- 0.5_8*dt*j(j,k,l,1)
  end do

  do l = 0, nz
      eyz(j,k,l) = afz(l)*eyz(j,k,l) + bpfz(l)*(bxz(j,k,l))  &
                                     + bmfz(l)*(bxz(j,k,l-1)) !- 0.5_8*dt*j(j,k,l,2)
  end do

 else if (.not. l_2drz) then

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

 else ! l_2drz=True

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

  if (xmin==0.) then
    do l = 0, nz-1
      j = 0
      ezx(j,k,l) = afx(j)*ezx(j,k,l) + 4*bpfx(j)*(byx(j,k,l)+byz(j,k,l))  
      do j = 1, nx
        ru = 1.+0.5/j
        rd = 1.-0.5/j
        ezx(j,k,l) = afx(j)*ezx(j,k,l) + ru*bpfx(j)*(byx(j,k,l)+byz(j,k,l))  &
                                       + rd*bmfx(j)*(byx(j-1,k,l)+byz(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,3)
      end do
    end do
  else
    do l = 0, nz-1
      do j = 0, nx
        ru = 1.+0.5/(xmin/dx+j)
        rd = 1.-0.5/(xmin/dx+j)
        ezx(j,k,l) = afx(j)*ezx(j,k,l) + ru*bpfx(j)*(byx(j,k,l)+byz(j,k,l))  &
                                       + rd*bmfx(j)*(byx(j-1,k,l)+byz(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,3)
      end do
    end do
  end if 

 end if
 
end if

  return
end subroutine push_em3d_splitevec


subroutine push_em3d_splitevec_norder(nx,ny,nz,nxguard,nyguard,nzguard, &
                               norderx,nordery,norderz, &
                               xcoefs,ycoefs,zcoefs, &
                               exy,exz,eyx,eyz,ezx,ezy,bxy,byx,bzx,bxz,byz,bzy, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz,l_1dz,l_2dxz,l_2drz,l_nodalgrid, &
                               xmin,ymin,zmin,dx,dy,dz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard,norderx,nordery,norderz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exy,exz, &
                                                                                                       eyx,eyz, &
                                                                                                       ezx,ezy
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz
real(kind=8), intent(in) :: xmin,ymin,zmin,dx,dy,dz
real(kind=8), intent(in) :: xcoefs(norderx/2),ycoefs(nordery/2),zcoefs(norderz/2)

INTEGER :: i, j, k, l, ist
logical(ISZ) :: l_1dz,l_2dxz,l_2drz,l_nodalgrid
real(8) :: ru,rd

if (l_nodalgrid) then
  ist = 0
else
  ist = 1
end if

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-ist
     exy(j,k,l) = afy(k)*exy(j,k,l) 
     do i = 1, nordery/2
      exy(j,k,l) = exy(j,k,l) + ycoefs(i)*bpfy(k)*(bzx(j,k+i-ist,l)+bzy(j,k+i-ist,l))  &
                              + ycoefs(i)*bmfy(k)*(bzx(j,k-i,l)+bzy(j,k-i,l)) !- 0.5_8*dt*j(j,k,l,1)
     end do
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-ist
     exz(j,k,l) = afz(l)*exz(j,k,l)
     do i = 1, norderz/2
      exz(j,k,l) = exz(j,k,l) - zcoefs(i)*bpfz(l)*(byx(j,k,l+i-ist)+byz(j,k,l+i-ist))  &
                              - zcoefs(i)*bmfz(l)*(byx(j,k,l-i)+byz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,1)
     end do
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-ist
    do j = 0, nx
     eyx(j,k,l) = afx(j)*eyx(j,k,l)
     do i = 1, norderx/2
      eyx(j,k,l) = eyx(j,k,l) - xcoefs(i)*bpfx(j)*(bzx(j+i-ist,k,l)+bzy(j+i-ist,k,l))  &
                              - xcoefs(i)*bmfx(j)*(bzx(j-i,k,l)+bzy(j-i,k,l)) !- 0.5_8*dt*j(j,k,l,2)
     end do
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-ist
    do j = 0, nx
     eyz(j,k,l) = afz(l)*eyz(j,k,l)
     do i = 1, norderz/2
      eyz(j,k,l) = eyz(j,k,l) + zcoefs(i)*bpfz(l)*(bxy(j,k,l+i-ist)+bxz(j,k,l+i-ist))  &
                              + zcoefs(i)*bmfz(l)*(bxy(j,k,l-i)+bxz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,2)
     end do
    end do
   end do
  end do

  do l = 0, nz-ist
   do k = 0, ny
    do j = 0, nx
     ezx(j,k,l) = afx(j)*ezx(j,k,l)
     do i = 1, norderx/2
      ezx(j,k,l) = ezx(j,k,l) + xcoefs(i)*bpfx(j)*(byx(j+i-ist,k,l)+byz(j+i-ist,k,l))  &
                              + xcoefs(i)*bmfx(j)*(byx(j-i,k,l)+byz(j-i,k,l)) !- 0.5_8*dt*j(j,k,l,3)
     end do
    end do
   end do
  end do

  do l = 0, nz-ist
   do k = 0, ny
    do j = 0, nx
     ezy(j,k,l) = afy(k)*ezy(j,k,l)
     do i = 1, nordery/2
      ezy(j,k,l) = ezy(j,k,l) - ycoefs(i)*bpfy(k)*(bxy(j,k+i-ist,l)+bxz(j,k+i-ist,l))  &
                              - ycoefs(i)*bmfy(k)*(bxy(j,k-i,l)+bxz(j,k-i,l)) !- 0.5_8*dt*j(j,k,l,3)
     end do
    end do
   end do
  end do

else
 k = 0

 if (l_1dz) then
  j = 0

  do l = 0, nz
     exz(j,k,l) = afz(l)*exz(j,k,l)
     do i = 1, norderz/2
      exz(j,k,l) = exz(j,k,l) - zcoefs(i)*bpfz(l)*(byx(j,k,l+i-ist)+byz(j,k,l+i-ist))  &
                                     - zcoefs(i)*bmfz(l)*(byx(j,k,l-i)+byz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,1)
     end do
  end do

  do l = 0, nz
     eyz(j,k,l) = afz(l)*eyz(j,k,l)
     do i = 1, norderz/2
      eyz(j,k,l) = eyz(j,k,l) + zcoefs(i)*bpfz(l)*(bxz(j,k,l+i-ist))  &
                                     + zcoefs(i)*bmfz(l)*(bxz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,2)
     end do
  end do

 else if (.not. l_2drz) then

  do l = 0, nz
    do j = 0, nx-ist
     exz(j,k,l) = afz(l)*exz(j,k,l)
     do i = 1, norderz/2
      exz(j,k,l) = exz(j,k,l) - zcoefs(i)*bpfz(l)*(byx(j,k,l+i-ist)+byz(j,k,l+i-ist))  &
                              - zcoefs(i)*bmfz(l)*(byx(j,k,l-i)+byz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,1)
     end do
    end do
  end do

  do l = 0, nz
    do j = 0, nx
     eyx(j,k,l) = afx(j)*eyx(j,k,l)
     do i = 1, norderx/2
      eyx(j,k,l) = eyx(j,k,l) - xcoefs(i)*bpfx(j)*(bzx(j+i-ist,k,l))  &
                              - xcoefs(i)*bmfx(j)*(bzx(j-i,k,l)) !- 0.5_8*dt*j(j,k,l,2)
     end do
    end do
  end do

  do l = 0, nz
    do j = 0, nx
     eyz(j,k,l) = afz(l)*eyz(j,k,l)
     do i = 1, norderz/2
      eyz(j,k,l) = eyz(j,k,l) + zcoefs(i)*bpfz(l)*(bxz(j,k,l+i-ist))  &
                              + zcoefs(i)*bmfz(l)*(bxz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,2)
     end do
    end do
  end do

  do l = 0, nz-ist
    do j = 0, nx
     ezx(j,k,l) = afx(j)*ezx(j,k,l)
     do i = 1, norderx/2
      ezx(j,k,l) = ezx(j,k,l) + xcoefs(i)*bpfx(j)*(byx(j+i-ist,k,l)+byz(j+i-ist,k,l))  &
                              + xcoefs(i)*bmfx(j)*(byx(j-i,k,l)+byz(j-i,k,l)) !- 0.5_8*dt*j(j,k,l,3)
     end do
    end do
  end do

 else ! l_2drz=True

  do l = 0, nz
    do j = 0, nx-ist
     exz(j,k,l) = afz(l)*exz(j,k,l)
     do i = 1, norderz/2
      exz(j,k,l) = exz(j,k,l) - zcoefs(i)*bpfz(l)*(byx(j,k,l+i-ist)+byz(j,k,l+i-ist))  &
                                     - zcoefs(i)*bmfz(l)*(byx(j,k,l-i)+byz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,1)
     end do
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
     eyz(j,k,l) = afz(l)*eyz(j,k,l)
     do i = 1, norderz/2
      eyz(j,k,l) = eyz(j,k,l) + zcoefs(i)*bpfz(l)*(bxz(j,k,l+i-ist))  &
                                     + zcoefs(i)*bmfz(l)*(bxz(j,k,l-i)) !- 0.5_8*dt*j(j,k,l,2)
     end do
    end do
  end do

  if (xmin==0.) then
    do l = 0, nz-1
      j = 0
      ezx(j,k,l) = afx(j)*ezx(j,k,l) + 4*bpfx(j)*(byx(j,k,l)+byz(j,k,l))  
      do j = 1, nx
        ru = 1.+0.5/j
        rd = 1.-0.5/j
        ezx(j,k,l) = afx(j)*ezx(j,k,l) + ru*bpfx(j)*(byx(j,k,l)+byz(j,k,l))  &
                                       + rd*bmfx(j)*(byx(j-1,k,l)+byz(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,3)
      end do
    end do
  else
    do l = 0, nz-1
      do j = 0, nx
        ru = 1.+0.5/(xmin/dx+j)
        rd = 1.-0.5/(xmin/dx+j)
        ezx(j,k,l) = afx(j)*ezx(j,k,l) + ru*bpfx(j)*(byx(j,k,l)+byz(j,k,l))  &
                                       + rd*bmfx(j)*(byx(j-1,k,l)+byz(j-1,k,l)) !- 0.5_8*dt*j(j,k,l,3)
      end do
    end do
  end if 

 end if
 
end if

  return
end subroutine push_em3d_splitevec_norder

subroutine scale_em3d_splitevec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exy,exz,eyx,eyz,ezx,ezy,&
                               sfx,sfy,sfz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exy,exz, &
                                                                                                       eyx,eyz, &
                                                                                                       ezx,ezy
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: sfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: sfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: sfz

INTEGER :: j, k, l

  do l = -nzguard, nz+nzguard
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard-1
      exy(j,k,l) = sfy(k)*exy(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard-1
      exz(j,k,l) = sfz(l)*exz(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard
      eyx(j,k,l) = sfx(j)*eyx(j,k,l)
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard
      eyz(j,k,l) = sfz(l)*eyz(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard
      ezx(j,k,l) = sfx(j)*ezx(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard
      ezy(j,k,l) = sfy(k)*ezy(j,k,l) 
    end do
   end do
  end do

end subroutine scale_em3d_splitevec

subroutine push_em3d_splite_setcond(nx,ny,nz,nxcond,nycond,nzcond,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,incond,l_2dxz,l_2drz, &
                               xmin,ymin,zmin,dx,dy,dz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard,nxcond,nycond,nzcond
logical(ISZ), dimension(-nxguard:nxcond+nxguard,-nyguard:nycond+nyguard,-nzguard:nzcond+nzguard), intent(in) :: incond
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx,exy,exz, &
                                                                                                       eyx,eyy,eyz, &
                                                                                                       ezx,ezy,ezz
real(kind=8), intent(in) :: xmin,ymin,zmin,dx,dy,dz
logical(ISZ) :: l_2dxz,l_2drz

INTEGER :: j, k, l

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
      if (incond(j,k,l) .and. incond(j+1,k,l)) then
        exx(j,k,l) = 0.
        exy(j,k,l) = 0.
        exz(j,k,l) = 0.
      end if
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx
      if (incond(j,k,l) .and. incond(j,k+1,l)) then
        eyx(j,k,l) = 0.
        eyy(j,k,l) = 0.
        eyz(j,k,l) = 0.
      end if
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      if (incond(j,k,l) .and. incond(j,k,l+1)) then
        ezx(j,k,l) = 0.
        ezy(j,k,l) = 0.
        ezz(j,k,l) = 0.
      end if
    end do
   end do
  end do

else
  k = 0

  do l = 0, nz
    do j = 0, nx-1
      if (incond(j,k,l) .and. incond(j+1,k,l)) then
        exx(j,k,l) = 0.
        exy(j,k,l) = 0.
        exz(j,k,l) = 0.
      end if
    end do
  end do

  do l = 0, nz
    do j = 0, nx
      if (incond(j,k,l)) then
        eyx(j,k,l) = 0.
        eyy(j,k,l) = 0.
        eyz(j,k,l) = 0.
      end if
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx
      if (incond(j,k,l) .and. incond(j,k,l+1)) then
        ezx(j,k,l) = 0.
        ezy(j,k,l) = 0.
        ezz(j,k,l) = 0.
      end if
    end do
  end do
 
end if

  return
end subroutine push_em3d_splite_setcond

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
    if(sf%nconds>0) then 
       call push_em3d_splite_setcond(sf%nx,sf%ny,sf%nz,sf%nxcond,sf%nycond,sf%nzcond,sf%nxguard,sf%nyguard,sf%nzguard, &
                                     sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz,sf%incond,sf%l_2dxz,sf%l_2drz, &
                                     sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
    end if
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

subroutine scale_em3d_splitefvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,eyy,ezz,&
                               sgx,sgy,sgz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx, &
                                                                                                       eyy, &
                                                                                                       ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: sgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: sgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: sgz

INTEGER :: j, k, l
logical(ISZ) :: l_1dz,l_2dxz,l_2drz

  do l = -nzguard, nz+nzguard
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard-1
      exx(j,k,l) = sgx(j)*exx(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard
      eyy(j,k,l) = sgy(k)*eyy(j,k,l)
    end do
   end do
  end do


  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard
      ezz(j,k,l) = sgz(l)*ezz(j,k,l) 
    end do
   end do
  end do

end subroutine scale_em3d_splitefvec

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
      exx(j,k,l) = agx(j)*exx(j,k,l) + bpgx(j)*alphax * ( fx(j+1,k  ,l  ) + fy(j+1,k  ,l  ) + fz(j+1,k  ,l  )) &
                                     + bpgx(j)*betaxy * ( fx(j+1,k+1,l  ) + fy(j+1,k+1,l  ) + fz(j+1,k+1,l  )  &
                                                      +   fx(j+1,k-1,l  ) + fy(j+1,k-1,l  ) + fz(j+1,k-1,l  )) &
                                     + bpgx(j)*betaxz * ( fx(j+1,k  ,l+1) + fy(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                      +   fx(j+1,k  ,l-1) + fy(j+1,k  ,l-1) + fz(j+1,k  ,l-1)) &
                                     + bpgx(j)*gammax * ( fx(j+1,k+1,l+1) + fy(j+1,k+1,l+1) + fz(j+1,k+1,l+1)  &
                                                      +   fx(j+1,k-1,l+1) + fy(j+1,k-1,l+1) + fz(j+1,k-1,l+1)  &
                                                      +   fx(j+1,k+1,l-1) + fy(j+1,k+1,l-1) + fz(j+1,k+1,l-1)  &
                                                      +   fx(j+1,k-1,l-1) + fy(j+1,k-1,l-1) + fz(j+1,k-1,l-1)) &
                                     + bmgx(j)*alphax * ( fx(j  ,k  ,l  ) + fy(j  ,k  ,l  ) + fz(j  ,k  ,l  )) &
                                     + bmgx(j)*betaxy * ( fx(j  ,k+1,l  ) + fy(j  ,k+1,l  ) + fz(j  ,k+1,l  )  &
                                                      +   fx(j  ,k-1,l  ) + fy(j  ,k-1,l  ) + fz(j  ,k-1,l  )) &
                                     + bmgx(j)*betaxz * ( fx(j  ,k  ,l+1) + fy(j  ,k  ,l+1) + fz(j  ,k  ,l+1)  &
                                                      +   fx(j  ,k  ,l-1) + fy(j  ,k  ,l-1) + fz(j  ,k  ,l-1)) &
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
      eyy(j,k,l) = agy(k)*eyy(j,k,l) + bpgy(k)*alphay * ( fx(j  ,k+1,l  ) + fy(j  ,k+1,l  ) + fz(j  ,k+1,l  )) &
                                     + bpgy(k)*betayx * ( fx(j+1,k+1,l  ) + fy(j+1,k+1,l  ) + fz(j+1,k+1,l  )  &
                                                      +   fx(j-1,k+1,l  ) + fy(j-1,k+1,l  ) + fz(j-1,k+1,l  )) &
                                     + bpgy(k)*betayz * ( fx(j  ,k+1,l+1) + fy(j  ,k+1,l+1) + fz(j  ,k+1,l+1)  &
                                                      +   fx(j  ,k+1,l-1) + fy(j  ,k+1,l-1) + fz(j  ,k+1,l-1)) &
                                     + bpgy(k)*gammay * ( fx(j+1,k+1,l+1) + fy(j+1,k+1,l+1) + fz(j+1,k+1,l+1)  &
                                                      +   fx(j-1,k+1,l+1) + fy(j-1,k+1,l+1) + fz(j-1,k+1,l+1)  &
                                                      +   fx(j+1,k+1,l-1) + fy(j+1,k+1,l-1) + fz(j+1,k+1,l-1)  &
                                                      +   fx(j-1,k+1,l-1) + fy(j-1,k+1,l-1) + fz(j-1,k+1,l-1)) &
                                     + bmgy(k)*alphay * ( fx(j  ,k  ,l  ) + fy(j  ,k  ,l  ) + fz(j  ,k  ,l  )) &
                                     + bmgy(k)*betayx * ( fx(j+1,k  ,l  ) + fy(j+1,k  ,l  ) + fz(j+1,k  ,l  )  &
                                                      +   fx(j-1,k  ,l  ) + fy(j-1,k  ,l  ) + fz(j-1,k  ,l  )) &
                                     + bmgy(k)*betayz * ( fx(j  ,k  ,l+1) + fy(j  ,k  ,l+1) + fz(j  ,k  ,l+1)  &
                                                      +   fx(j  ,k  ,l-1) + fy(j  ,k  ,l-1) + fz(j  ,k  ,l-1)) &
                                     + bmgy(k)*gammay * ( fx(j+1,k  ,l+1) + fy(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                      +   fx(j-1,k  ,l+1) + fy(j-1,k  ,l+1) + fz(j-1,k  ,l+1)  &
                                                      +   fx(j+1,k  ,l-1) + fy(j+1,k  ,l-1) + fz(j+1,k  ,l-1)  &
                                                      +   fx(j-1,k  ,l-1) + fy(j-1,k  ,l-1) + fz(j-1,k  ,l-1)) 
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx
      ezz(j,k,l) = agz(l)*ezz(j,k,l) + bpgz(l)*alphaz * ( fx(j  ,k  ,l+1) + fy(j  ,k  ,l+1) + fz(j  ,k  ,l+1)) &
                                     + bpgz(l)*betazx * ( fx(j+1,k  ,l+1) + fy(j+1,k  ,l+1) + fz(j+1,k  ,l+1)  &
                                                      +   fx(j-1,k  ,l+1) + fy(j-1,k  ,l+1) + fz(j-1,k  ,l+1)) &
                                     + bpgz(l)*betazy * ( fx(j  ,k+1,l+1) + fy(j  ,k+1,l+1) + fz(j  ,k+1,l+1)  &
                                                      +   fx(j  ,k-1,l+1) + fy(j  ,k-1,l+1) + fz(j  ,k-1,l+1)) &
                                     + bpgz(l)*gammaz * ( fx(j+1,k+1,l+1) + fy(j+1,k+1,l+1) + fz(j+1,k+1,l+1)  &
                                                      +   fx(j-1,k+1,l+1) + fy(j-1,k+1,l+1) + fz(j-1,k+1,l+1)  &
                                                      +   fx(j+1,k-1,l+1) + fy(j+1,k-1,l+1) + fz(j+1,k-1,l+1)  &
                                                      +   fx(j-1,k-1,l+1) + fy(j-1,k-1,l+1) + fz(j-1,k-1,l+1)) &
                                     + bmgz(l)*alphaz * ( fx(j  ,k  ,l  ) + fy(j  ,k  ,l  ) + fz(j  ,k  ,l  )) &
                                     + bmgz(l)*betazx * ( fx(j+1,k  ,l  ) + fy(j+1,k  ,l  ) + fz(j+1,k  ,l  )  &
                                                      +   fx(j-1,k  ,l  ) + fy(j-1,k  ,l  ) + fz(j-1,k  ,l  )) &
                                     + bmgz(l)*betazy * ( fx(j  ,k+1,l  ) + fy(j  ,k+1,l  ) + fz(j  ,k+1,l  )  &
                                                      +   fx(j  ,k-1,l  ) + fy(j  ,k-1,l  ) + fz(j  ,k-1,l  )) &
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
                                     + bpgx(j)*    betaxz * ( fx(j+1,k  ,l+1) + fz(j+1,k  ,l+1)   &
                                                          +   fx(j+1,k  ,l-1) + fz(j+1,k  ,l-1))  &
                                     + bmgx(j)*    alphax * ( fx(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgx(j)*    betaxz * ( fx(j  ,k  ,l+1) + fz(j  ,k  ,l+1)   &
                                                          +   fx(j  ,k  ,l-1) + fz(j  ,k  ,l-1))  
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx
      ezz(j,k,l) = agz(l)*ezz(j,k,l) + bpgz(l)*    alphaz * ( fx(j  ,k  ,l+1) + fz(j  ,k  ,l+1) ) &
                                     + bpgz(l)*    betazx * ( fx(j+1,k  ,l+1) + fz(j+1,k  ,l+1)   &
                                                          +   fx(j-1,k  ,l+1) + fz(j-1,k  ,l+1))  &
                                     + bmgz(l)*    alphaz * ( fx(j  ,k  ,l  ) + fz(j  ,k  ,l  ) ) &
                                     + bmgz(l)*    betazx * ( fx(j+1,k  ,l  ) + fz(j+1,k  ,l  )   &
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
   if ((sf%norderx==2) .and. (sf%nordery==2) .and. (sf%norderz==2) .and. .not. sf%l_nodalgrid) then
    call push_em3d_splitbvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%nxbs,sf%nybs,sf%nzbs, &
                             sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                             sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                             sf%agx,sf%agy,sf%agz, &
                             sf%bpgx,sf%bpgy,sf%bpgz, &
                             sf%bmgx,sf%bmgy,sf%bmgz,sf%l_1dz,sf%l_2dxz,sf%l_2drz, &
                             sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
   else
    call push_em3d_splitbvec_norder(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%nxbs,sf%nybs,sf%nzbs, &
                             sf%norderx,sf%nordery,sf%norderz, &
                             sf%xcoefs,sf%ycoefs,sf%zcoefs, &
                             sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                             sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                             sf%agx,sf%agy,sf%agz, &
                             sf%bpgx,sf%bpgy,sf%bpgz, &
                             sf%bmgx,sf%bmgy,sf%bmgz,sf%l_1dz,sf%l_2dxz,sf%l_2drz,sf%l_nodalgrid, &
                             sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
   end if
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

subroutine push_em3d_splitbvec(nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz,l_1dz,l_2dxz,l_2drz, &
                               xmin,ymin,zmin,dx,dy,dz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxs,nys,nzs,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz
real(kind=8), intent(in) :: xmin,ymin,zmin,dx,dy,dz
!real(kind=8), dimension(-nxguard:,-nyguard:,-nzguard:), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
!real(kind=8), dimension(-nxguard:,-nyguard:,-nzguard:), intent(in) :: exx,exy,exz, &
!                                                                                                    eyx,eyy,eyz, &
!                                                                                                    ezx,ezy,ezz
!real(kind=8), dimension(-nxguard:), intent(in) :: agx,bpgx,bmgx
!real(kind=8), dimension(-nyguard:), intent(in) :: agy,bpgy,bmgy
!real(kind=8), dimension(-nzguard:), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l
logical(ISZ) :: l_1dz, l_2dxz,l_2drz
real(8) :: ru

if (.not.l_2dxz) then

  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs
      bxy(j,k,l) = agy(k)*bxy(j,k,l) - bpgy(k)*(ezx(j,k+1,l  )+ezy(j,k+1,l  )+ezz(j,k+1,l  )) &
                                     - bmgy(k)*(ezx(j,k  ,l  )+ezy(j,k  ,l  )+ezz(j,k  ,l  ))
    end do
   end do
  end do

  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*(eyx(j,k  ,l+1)+eyy(j,k  ,l+1)+eyz(j,k  ,l+1)) &
                                     + bmgz(l)*(eyx(j,k  ,l  )+eyy(j,k  ,l  )+eyz(j,k  ,l  ))
    end do
   end do
  end do

  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*(ezx(j+1,k,l  )+ezy(j+1,k,l  )+ezz(j+1,k,l  )) &
                                     + bmgx(j)*(ezx(j  ,k,l  )+ezy(j  ,k,l  )+ezz(j  ,k,l  ))
    end do
   end do
  end do

  do l = -nzs, nz+nzs-1
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*(exx(j  ,k,l+1)+exy(j  ,k,l+1)+exz(j  ,k,l+1)) &
                                     - bmgz(l)*(exx(j  ,k,l  )+exy(j  ,k,l  )+exz(j  ,k,l  ))
    end do
   end do
  end do

  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*(eyx(j+1,k  ,l)+eyy(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                     - bmgx(j)*(eyx(j  ,k  ,l)+eyy(j  ,k  ,l)+eyz(j  ,k  ,l))
    end do
   end do
  end do

  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-1
    do j = -nxs, nx+nxs-1
      bzy(j,k,l) = agy(k)*bzy(j,k,l) + bpgy(k)*(exx(j  ,k+1,l)+exy(j  ,k+1,l)+exz(j  ,k+1,l)) &
                                     + bmgy(k)*(exx(j  ,k  ,l)+exy(j  ,k  ,l)+exz(j  ,k  ,l))
    end do
   end do
  end do

else
  if (l_1dz) then
   j = 0
   k = 0
   do l = -nzs, nz+nzs-1
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*(eyx(j,k  ,l+1)+eyz(j,k  ,l+1)) &
                                     + bmgz(l)*(eyx(j,k  ,l  )+eyz(j,k  ,l  ))
   end do

   do l = -nzs, nz+nzs-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*(exx(j  ,k,l+1)+exz(j  ,k,l+1)) &
                                     - bmgz(l)*(exx(j  ,k,l  )+exz(j  ,k,l  ))
   end do

  else
   k = 0
  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*(eyx(j,k  ,l+1)+eyz(j,k  ,l+1)) &
                                     + bmgz(l)*(eyx(j,k  ,l  )+eyz(j,k  ,l  ))
    end do
   end do

  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*(ezx(j+1,k,l  )+ezz(j+1,k,l  )) &
                                     + bmgx(j)*(ezx(j  ,k,l  )+ezz(j  ,k,l  ))
    end do
   end do

  do l = -nzs, nz+nzs-1
    do j = -nxs, nx+nxs-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*(exx(j  ,k,l+1)+exz(j  ,k,l+1)) &
                                     - bmgz(l)*(exx(j  ,k,l  )+exz(j  ,k,l  ))
    end do
   end do

   if (l_2drz) then
    do l = -nzs, nz+nzs
      do j = -nxs, nx+nxs-1
        ru = (xmin+(j+1)*dx)/(xmin+j*dx+0.5*dx)
        bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*ru*(eyx(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                       - bmgx(j)*(eyx(j  ,k  ,l)+eyz(j  ,k  ,l))
      end do
    end do
   else
    do l = -nzs, nz+nzs
      do j = -nxs, nx+nxs-1
        bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*(eyx(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                       - bmgx(j)*(eyx(j  ,k  ,l)+eyz(j  ,k  ,l))
      end do
    end do
   end if
  end if

end if

  return
end subroutine push_em3d_splitbvec

subroutine push_em3d_splitbvec_norder(nx,ny,nz,nxguard,nyguard,nzguard,nxs,nys,nzs, &
                               norderx,nordery,norderz, &
                               xcoefs,ycoefs,zcoefs, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz,l_1dz,l_2dxz,l_2drz,l_nodalgrid, &
                               xmin,ymin,zmin,dx,dy,dz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxs,nys,nzs,nxguard,nyguard,nzguard,norderx,nordery,norderz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz
real(kind=8), intent(in) :: xmin,ymin,zmin,dx,dy,dz
real(kind=8), intent(in) :: xcoefs(norderx/2),ycoefs(nordery/2),zcoefs(norderz/2)

INTEGER :: i, j, k, l, ist
logical(ISZ) :: l_1dz, l_2dxz,l_2drz,l_nodalgrid
real(8) :: ru

if (l_nodalgrid) then
  ist = 0
else
  ist = 1
end if

if (.not.l_2dxz) then

  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs
     bxy(j,k,l) = agy(k)*bxy(j,k,l)
     do i = 1, nordery/2
      bxy(j,k,l) = bxy(j,k,l) - ycoefs(i)*bpgy(k)*(ezx(j,k+i,l  )+ezy(j,k+i,l  )+ezz(j,k+i,l  )) &
                              - ycoefs(i)*bmgy(k)*(ezx(j,k-i+ist,l  )+ezy(j,k-i+ist,l  )+ezz(j,k-i+ist,l  ))
     end do
    end do
   end do
  end do

  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs
     bxz(j,k,l) = agz(l)*bxz(j,k,l)
     do i = 1, norderx/2
      bxz(j,k,l) = bxz(j,k,l) + zcoefs(i)*bpgz(l)*(eyx(j,k  ,l+i)+eyy(j,k  ,l+i)+eyz(j,k  ,l+i)) &
                              + zcoefs(i)*bmgz(l)*(eyx(j,k  ,l-i+ist)+eyy(j,k  ,l-i+ist)+eyz(j,k  ,l-i+ist))
     end do
    end do
   end do
  end do

  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-ist
     byx(j,k,l) = agx(j)*byx(j,k,l)
     do i = 1, norderx/2
      byx(j,k,l) = byx(j,k,l) + xcoefs(i)*bpgx(j)*(ezx(j+i,k,l  )+ezy(j+i,k,l  )+ezz(j+i,k,l  )) &
                              + xcoefs(i)*bmgx(j)*(ezx(j-i+ist,k,l  )+ezy(j-i+ist,k,l  )+ezz(j-i+ist,k,l  ))
     end do
    end do
   end do
  end do

  do l = -nzs, nz+nzs-ist
   do k = -nys, ny+nys
    do j = -nxs, nx+nxs-ist
     byz(j,k,l) = agz(l)*byz(j,k,l)
     do i = 1, norderz/2
      byz(j,k,l) = byz(j,k,l) - zcoefs(i)*bpgz(l)*(exx(j  ,k,l+i)+exy(j  ,k,l+i)+exz(j  ,k,l+i)) &
                              - zcoefs(i)*bmgz(l)*(exx(j  ,k,l-i+ist)+exy(j  ,k,l-i+ist)+exz(j  ,k,l-i+ist))
     end do
    end do
   end do
  end do

  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs-ist
     bzx(j,k,l) = agx(j)*bzx(j,k,l)
     do i = 1, norderx/2
      bzx(j,k,l) = bzx(j,k,l) - xcoefs(i)*bpgx(j)*(eyx(j+i,k  ,l)+eyy(j+i,k  ,l)+eyz(j+i,k  ,l)) &
                              - xcoefs(i)*bmgx(j)*(eyx(j-i+ist,k  ,l)+eyy(j-i+ist,k  ,l)+eyz(j-i+ist,k  ,l))
     end do
    end do
   end do
  end do

  do l = -nzs, nz+nzs
   do k = -nys, ny+nys-ist
    do j = -nxs, nx+nxs-ist
     bzy(j,k,l) = agy(k)*bzy(j,k,l)
     do i = 1, nordery/2
      bzy(j,k,l) = bzy(j,k,l) + ycoefs(i)*bpgy(k)*(exx(j  ,k+i,l)+exy(j  ,k+i,l)+exz(j  ,k+i,l)) &
                              + ycoefs(i)*bmgy(k)*(exx(j  ,k-i+ist,l)+exy(j  ,k-i+ist,l)+exz(j  ,k-i+ist,l))
     end do
    end do
   end do
  end do

else
  if (l_1dz) then
   j = 0
   k = 0
   do l = -nzs, nz+nzs-ist
     bxz(j,k,l) = agz(l)*bxz(j,k,l)
     do i = 1, norderz/2
      bxz(j,k,l) = bxz(j,k,l) + zcoefs(i)*bpgz(l)*(eyx(j,k  ,l+i)+eyz(j,k  ,l+i)) &
                              + zcoefs(i)*bmgz(l)*(eyx(j,k  ,l-i+ist)+eyz(j,k  ,l-i+ist))
     end do
   end do

   do l = -nzs, nz+nzs-ist
     byz(j,k,l) = agz(l)*byz(j,k,l)
     do i = 1, norderz/2
      byz(j,k,l) = byz(j,k,l) - zcoefs(i)*bpgz(l)*(exx(j  ,k,l+i)+exz(j  ,k,l+i)) &
                              - zcoefs(i)*bmgz(l)*(exx(j  ,k,l-i+ist)+exz(j  ,k,l-i+ist))
     end do
   end do

  else
   k = 0
  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs
     bxz(j,k,l) = agz(l)*bxz(j,k,l)
     do i = 1, norderz/2
      bxz(j,k,l) = bxz(j,k,l) + zcoefs(i)*bpgz(l)*(eyx(j,k  ,l+i)+eyz(j,k  ,l+i)) &
                              + zcoefs(i)*bmgz(l)*(eyx(j,k  ,l-i+ist)+eyz(j,k  ,l-i+ist))
     end do
    end do
   end do

  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs-ist
     byx(j,k,l) = agx(j)*byx(j,k,l)
     do i = 1, norderx/2
      byx(j,k,l) = byx(j,k,l) + xcoefs(i)*bpgx(j)*(ezx(j+i,k,l  )+ezz(j+i,k,l  )) &
                              + xcoefs(i)*bmgx(j)*(ezx(j-i+ist,k,l  )+ezz(j-i+ist,k,l  ))
     end do
    end do
   end do

  do l = -nzs, nz+nzs-ist
    do j = -nxs, nx+nxs-ist
     byz(j,k,l) = agz(l)*byz(j,k,l)
     do i = 1, norderz/2
      byz(j,k,l) = byz(j,k,l) - zcoefs(i)*bpgz(l)*(exx(j  ,k,l+i)+exz(j  ,k,l+i)) &
                              - zcoefs(i)*bmgz(l)*(exx(j  ,k,l-i+ist)+exz(j  ,k,l-i+ist))
     end do
    end do
   end do

   if (l_2drz) then
    do l = -nzs, nz+nzs
      do j = -nxs, nx+nxs-1
        ru = (xmin+(j+1)*dx)/(xmin+j*dx+0.5*dx)
        bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*ru*(eyx(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                       - bmgx(j)*(eyx(j  ,k  ,l)+eyz(j  ,k  ,l))
      end do
    end do
   else
    do l = -nzs, nz+nzs
     do j = -nxs, nx+nxs-ist
       bzx(j,k,l) = agx(j)*bzx(j,k,l)
       do i = 1, norderx/2
        bzx(j,k,l) = bzx(j,k,l) - xcoefs(i)*bpgx(j)*(eyx(j+i,k  ,l)+eyz(j+i,k  ,l)) &
                                - xcoefs(i)*bmgx(j)*(eyx(j-i+ist,k  ,l)+eyz(j-i+ist,k  ,l))
       end do
     end do
    end do
   end if
  end if

end if

  return
end subroutine push_em3d_splitbvec_norder

subroutine scale_em3d_splitbvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               bxy,bxz,byx,byz,bzx,bzy,&
                               sgx,sgy,sgz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: bxy,bxz, &
                                                                                                       byx,byz, &
                                                                                                       bzx,bzy
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: sgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: sgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: sgz

INTEGER :: j, k, l

  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard
      bxy(j,k,l) = sgy(k)*bxy(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard
      bxz(j,k,l) = sgz(l)*bxz(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard-1
      byx(j,k,l) = sgx(j)*byx(j,k,l)
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard-1
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard-1
      byz(j,k,l) = sgz(l)*byz(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard-1
      bzx(j,k,l) = sgx(j)*bzx(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, max(0,ny+nyguard-1)
    do j = -nxguard, nx+nxguard-1
      bzy(j,k,l) = sgy(k)*bzy(j,k,l) 
    end do
   end do
  end do

end subroutine scale_em3d_splitbvec

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

real(8) :: b,c,dt

if (.not.l_2dxz) then

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxy(j,k,l) = agy(k)*bxy(j,k,l) - bpgy(k)*alphay * (ezx(j  ,k+1,l  )+ezy(j  ,k+1,l  )+ezz(j  ,k+1,l  )) &
                                     - bpgy(k)*betayx * (ezx(j+1,k+1,l  )+ezy(j+1,k+1,l  )+ezz(j+1,k+1,l  ) &
                                                      +  ezx(j-1,k+1,l  )+ezy(j-1,k+1,l  )+ezz(j-1,k+1,l  )) &
                                     - bpgy(k)*betayz * (ezx(j  ,k+1,l+1)+ezy(j  ,k+1,l+1)+ezz(j  ,k+1,l+1) &
                                                      +  ezx(j  ,k+1,l-1)+ezy(j  ,k+1,l-1)+ezz(j  ,k+1,l-1)) &
                                     - bpgy(k)*gammay * (ezx(j+1,k+1,l+1)+ezy(j+1,k+1,l+1)+ezz(j+1,k+1,l+1) &
                                                      +  ezx(j-1,k+1,l+1)+ezy(j-1,k+1,l+1)+ezz(j-1,k+1,l+1) &
                                                      +  ezx(j+1,k+1,l-1)+ezy(j+1,k+1,l-1)+ezz(j+1,k+1,l-1) &
                                                      +  ezx(j-1,k+1,l-1)+ezy(j-1,k+1,l-1)+ezz(j-1,k+1,l-1)) &
                                     - bmgy(k)*alphay * (ezx(j  ,k  ,l  )+ezy(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     - bmgy(k)*betayx * (ezx(j+1,k  ,l  )+ezy(j+1,k  ,l  )+ezz(j+1,k  ,l  ) &
                                                      +  ezx(j-1,k  ,l  )+ezy(j-1,k  ,l  )+ezz(j-1,k  ,l  )) &
                                     - bmgy(k)*betayz * (ezx(j  ,k  ,l+1)+ezy(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
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
                                     + bpgz(l)*betazx * (eyx(j+1,k  ,l+1)+eyy(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                      +  eyx(j-1,k  ,l+1)+eyy(j-1,k  ,l+1)+eyz(j-1,k  ,l+1)) &
                                     + bpgz(l)*betazy * (eyx(j  ,k+1,l+1)+eyy(j  ,k+1,l+1)+eyz(j  ,k+1,l+1) &
                                                      +  eyx(j  ,k-1,l+1)+eyy(j  ,k-1,l+1)+eyz(j  ,k-1,l+1)) &
                                     + bpgz(l)*gammaz * (eyx(j+1,k+1,l+1)+eyy(j+1,k+1,l+1)+eyz(j+1,k+1,l+1) &
                                                      +  eyx(j-1,k+1,l+1)+eyy(j-1,k+1,l+1)+eyz(j-1,k+1,l+1) &
                                                      +  eyx(j+1,k-1,l+1)+eyy(j+1,k-1,l+1)+eyz(j+1,k-1,l+1) &
                                                      +  eyx(j-1,k-1,l+1)+eyy(j-1,k-1,l+1)+eyz(j-1,k-1,l+1)) &
                                     + bmgz(l)*alphaz * (eyx(j  ,k  ,l  )+eyy(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     + bmgz(l)*betazx * (eyx(j+1,k  ,l  )+eyy(j+1,k  ,l  )+eyz(j+1,k  ,l  ) &
                                                      +  eyx(j-1,k  ,l  )+eyy(j-1,k  ,l  )+eyz(j-1,k  ,l  )) &
                                     + bmgz(l)*betazy * (eyx(j  ,k+1,l  )+eyy(j  ,k+1,l  )+eyz(j  ,k+1,l  ) &
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
                                     + bpgx(j)*betaxy * (ezx(j+1,k+1,l  )+ezy(j+1,k+1,l  )+ezz(j+1,k+1,l  ) &
                                                      +  ezx(j+1,k-1,l  )+ezy(j+1,k-1,l  )+ezz(j+1,k-1,l  )) &
                                     + bpgx(j)*betaxz * (ezx(j+1,k  ,l+1)+ezy(j+1,k  ,l+1)+ezz(j+1,k  ,l+1) &
                                                      +  ezx(j+1,k  ,l-1)+ezy(j+1,k  ,l-1)+ezz(j+1,k  ,l-1)) &
                                     + bpgx(j)*gammax * (ezx(j+1,k+1,l+1)+ezy(j+1,k+1,l+1)+ezz(j+1,k+1,l+1) &
                                                      +  ezx(j+1,k-1,l+1)+ezy(j+1,k-1,l+1)+ezz(j+1,k-1,l+1) &
                                                      +  ezx(j+1,k+1,l-1)+ezy(j+1,k+1,l-1)+ezz(j+1,k+1,l-1) &
                                                      +  ezx(j+1,k-1,l-1)+ezy(j+1,k-1,l-1)+ezz(j+1,k-1,l-1)) &
                                     + bmgx(j)*alphax * (ezx(j  ,k  ,l  )+ezy(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     + bmgx(j)*betaxy * (ezx(j  ,k+1,l  )+ezy(j  ,k+1,l  )+ezz(j  ,k+1,l  ) &
                                                      +  ezx(j  ,k-1,l  )+ezy(j  ,k-1,l  )+ezz(j  ,k-1,l  )) &
                                     + bmgx(j)*betaxz * (ezx(j  ,k  ,l+1)+ezy(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
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
                                     - bpgz(l)*betazx * (exx(j+1,k  ,l+1)+exy(j+1,k  ,l+1)+exz(j+1,k  ,l+1) &
                                                      +  exx(j-1,k  ,l+1)+exy(j-1,k  ,l+1)+exz(j-1,k  ,l+1)) &
                                     - bpgz(l)*betazy * (exx(j  ,k+1,l+1)+exy(j  ,k+1,l+1)+exz(j  ,k+1,l+1) &
                                                      +  exx(j  ,k-1,l+1)+exy(j  ,k-1,l+1)+exz(j  ,k-1,l+1)) &
                                     - bpgz(l)*gammaz * (exx(j+1,k+1,l+1)+exy(j+1,k+1,l+1)+exz(j+1,k+1,l+1) &
                                                      +  exx(j-1,k+1,l+1)+exy(j-1,k+1,l+1)+exz(j-1,k+1,l+1) &
                                                      +  exx(j+1,k-1,l+1)+exy(j+1,k-1,l+1)+exz(j+1,k-1,l+1) &
                                                      +  exx(j-1,k-1,l+1)+exy(j-1,k-1,l+1)+exz(j-1,k-1,l+1)) &
                                     - bmgz(l)*alphaz * (exx(j  ,k  ,l  )+exy(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     - bmgz(l)*betazx * (exx(j+1,k  ,l  )+exy(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                      +  exx(j-1,k  ,l  )+exy(j-1,k  ,l  )+exz(j-1,k  ,l  )) &
                                     - bmgz(l)*betazy * (exx(j  ,k+1,l  )+exy(j  ,k+1,l  )+exz(j  ,k+1,l  ) &
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
                                     - bpgx(j)*betaxy * (eyx(j+1,k+1,l  )+eyy(j+1,k+1,l  )+eyz(j+1,k+1,l  ) &
                                                      +  eyx(j+1,k-1,l  )+eyy(j+1,k-1,l  )+eyz(j+1,k-1,l  )) &
                                     - bpgx(j)*betaxz * (eyx(j+1,k  ,l+1)+eyy(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                      +  eyx(j+1,k  ,l-1)+eyy(j+1,k  ,l-1)+eyz(j+1,k  ,l-1)) &
                                     - bpgx(j)*gammax * (eyx(j+1,k+1,l+1)+eyy(j+1,k+1,l+1)+eyz(j+1,k+1,l+1) &
                                                      +  eyx(j+1,k-1,l+1)+eyy(j+1,k-1,l+1)+eyz(j+1,k-1,l+1) &
                                                      +  eyx(j+1,k+1,l-1)+eyy(j+1,k+1,l-1)+eyz(j+1,k+1,l-1) &
                                                      +  eyx(j+1,k-1,l-1)+eyy(j+1,k-1,l-1)+eyz(j+1,k-1,l-1)) &
                                     - bmgx(j)*alphax * (eyx(j  ,k  ,l  )+eyy(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     - bmgx(j)*betaxy * (eyx(j  ,k+1,l  )+eyy(j  ,k+1,l  )+eyz(j  ,k+1,l  ) &
                                                      +  eyx(j  ,k-1,l  )+eyy(j  ,k-1,l  )+eyz(j  ,k-1,l  )) &
                                     - bmgx(j)*betaxz * (eyx(j  ,k  ,l+1)+eyy(j  ,k  ,l+1)+eyz(j  ,k  ,l+1) &
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
                                     + bpgy(k)*betayx * (exx(j+1,k+1,l  )+exy(j+1,k+1,l  )+exz(j+1,k+1,l  ) &
                                                      +  exx(j-1,k+1,l  )+exy(j-1,k+1,l  )+exz(j-1,k+1,l  )) &
                                     + bpgy(k)*betayz * (exx(j  ,k+1,l+1)+exy(j  ,k+1,l+1)+exz(j  ,k+1,l+1) &
                                                      +  exx(j  ,k+1,l-1)+exy(j  ,k+1,l-1)+exz(j  ,k+1,l-1)) &
                                     + bpgy(k)*gammay * (exx(j+1,k+1,l+1)+exy(j+1,k+1,l+1)+exz(j+1,k+1,l+1) &
                                                      +  exx(j-1,k+1,l+1)+exy(j-1,k+1,l+1)+exz(j-1,k+1,l+1) &
                                                      +  exx(j+1,k+1,l-1)+exy(j+1,k+1,l-1)+exz(j+1,k+1,l-1) &
                                                      +  exx(j-1,k+1,l-1)+exy(j-1,k+1,l-1)+exz(j-1,k+1,l-1)) &
                                     + bmgy(k)*alphay * (exx(j  ,k  ,l  )+exy(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     + bmgy(k)*betayx * (exx(j+1,k  ,l  )+exy(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                      +  exx(j-1,k  ,l  )+exy(j-1,k  ,l  )+exz(j-1,k  ,l  )) &
                                     + bmgy(k)*betayz * (exx(j  ,k  ,l+1)+exy(j  ,k  ,l+1)+exz(j  ,k  ,l+1) &
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

  if (.false.) then
 do l = 0, nz-1
    do j = 0, nx
      b = 0.5*(bpgz(l)-bmgz(l))
      c = 0.5*(bpgz(l)+bmgz(l))
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + b*    alphaz * (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1)) &
                                     + b*    betazx * (eyx(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                    +  eyx(j-1,k  ,l+1)+eyz(j-1,k  ,l+1)) &
                                     - b*    alphaz * (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     - b*    betazx * (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  ) &
                                                    +  eyx(j-1,k  ,l  )+eyz(j-1,k  ,l  )) &
                                     + c* (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1)  &
                                     -    (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  ))) 
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      b = 0.5*(bpgx(j)-bmgx(j))
      c = 0.5*(bpgx(j)+bmgx(j))
      byx(j,k,l) = agx(j)*byx(j,k,l) + b*    alphax * (ezx(j+1,k  ,l  )+ezz(j+1,k  ,l  )) &
                                     + b*    betaxz * (ezx(j+1,k  ,l+1)+ezz(j+1,k  ,l+1) &
                                                          +  ezx(j+1,k  ,l-1)+ezz(j+1,k  ,l-1)) &
                                     - b*    alphax * (ezx(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     - b*    betaxz * (ezx(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
                                                          +  ezx(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) &
                                     + c*(ezx(j+1,k  ,l  )+ezz(j+1,k  ,l  ) &
                                     -   (ezx(j  ,k  ,l  )+ezz(j  ,k  ,l  )))
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      b = 0.5*(bpgz(l)-bmgz(l))
      c = 0.5*(bpgz(l)+bmgz(l))
      byz(j,k,l) = agz(l)*byz(j,k,l) - b*    alphaz * (exx(j  ,k  ,l+1)+exz(j  ,k  ,l+1)) &
                                     - b*    betazx * (exx(j+1,k  ,l+1)+exz(j+1,k  ,l+1) &
                                                          +  exx(j-1,k  ,l+1)+exz(j-1,k  ,l+1)) &
                                     + b*    alphaz * (exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     + b*    betazx * (exx(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                          +  exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) &
                                     - c*(exx(j  ,k  ,l+1)+exz(j  ,k  ,l+1) &
                                     -   (exx(j  ,k  ,l  )+exz(j  ,k  ,l  )))

    end do
  end do

  do l = 0, nz
    do j = 0, nx-1
      b = 0.5*(bpgx(j)-bmgx(j))
      c = 0.5*(bpgx(j)+bmgx(j))
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - b*    alphax * (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  )) &
                                     - b*    betaxz * (eyx(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                          +  eyx(j+1,k  ,l-1)+eyz(j+1,k  ,l-1)) &
                                     + b*    alphax * (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     + b*    betaxz * (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1) &
                                                          +  eyx(j  ,k  ,l-1)+eyz(j  ,k  ,l-1)) &
                                     - c* (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  ) &
                                     -    (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )))
    end do
  end do



  else


  do l = 0, nz-1
    do j = 0, nx
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*    alphaz * (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1)) &
                                     + bpgz(l)*    betazx * (eyx(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                          +  eyx(j-1,k  ,l+1)+eyz(j-1,k  ,l+1)) &
                                     + bmgz(l)*    alphaz * (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     + bmgz(l)*    betazx * (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  ) &
                                                          +  eyx(j-1,k  ,l  )+eyz(j-1,k  ,l  )) 
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*    alphax * (ezx(j+1,k  ,l  )+ezz(j+1,k  ,l  )) &
                                     + bpgx(j)*    betaxz * (ezx(j+1,k  ,l+1)+ezz(j+1,k  ,l+1) &
                                                          +  ezx(j+1,k  ,l-1)+ezz(j+1,k  ,l-1)) &
                                     + bmgx(j)*    alphax * (ezx(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                     + bmgx(j)*    betaxz * (ezx(j  ,k  ,l+1)+ezz(j  ,k  ,l+1) &
                                                          +  ezx(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) 
    end do
  end do

  do l = 0, nz-1
    do j = 0, nx-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*    alphaz * (exx(j  ,k  ,l+1)+exz(j  ,k  ,l+1)) &
                                     - bpgz(l)*    betazx * (exx(j+1,k  ,l+1)+exz(j+1,k  ,l+1) &
                                                          +  exx(j-1,k  ,l+1)+exz(j-1,k  ,l+1)) &
                                     - bmgz(l)*    alphaz * (exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                     - bmgz(l)*    betazx * (exx(j+1,k  ,l  )+exz(j+1,k  ,l  ) &
                                                          +  exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) 

    end do
  end do

  do l = 0, nz
    do j = 0, nx-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*    alphax * (eyx(j+1,k  ,l  )+eyz(j+1,k  ,l  )) &
                                     - bpgx(j)*    betaxz * (eyx(j+1,k  ,l+1)+eyz(j+1,k  ,l+1) &
                                                          +  eyx(j+1,k  ,l-1)+eyz(j+1,k  ,l-1)) &
                                     - bmgx(j)*    alphax * (eyx(j  ,k  ,l  )+eyz(j  ,k  ,l  )) &
                                     - bmgx(j)*    betaxz * (eyx(j  ,k  ,l+1)+eyz(j  ,k  ,l+1) &
                                                          +  eyx(j  ,k  ,l-1)+eyz(j  ,k  ,l-1)) 
    end do
  end do
  end if
  
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
                           sf%bmfx,sf%bmfy,sf%bmfz,sf%l_2dxz,sf%l_2drz, &
                           sf%xmin,sf%ymin,sf%zmin,sf%dx,sf%dy,sf%dz)
  if(sf%nconds>0) then 
     call push_em3d_splitf_setcond(sf%nx,sf%ny,sf%nz,sf%nxcond,sf%nycond,sf%nzcond,sf%nxguard,sf%nyguard,sf%nzguard, &
                                   sf%fx,sf%fy,sf%fz,sf%l_2dxz,sf%incond)
  end if

  return
end subroutine push_em3d_splitf

subroutine scale_em3d_split_fields(sf,dt,l_pushf)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt
logical:: l_pushf

  if (sf%pml_method==2) then
    call set_bndcoeffsem3d(sf,dt,0)
    call scale_em3d_splitbvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%bxy,sf%bxz,sf%byx,sf%byz,sf%bzx,sf%bzy, &
                           sf%sgx,sf%sgy,sf%sgz)
    call scale_em3d_splitevec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%exy,sf%exz,sf%eyx,sf%eyz,sf%ezx,sf%ezy, &
                           sf%sfx,sf%sfy,sf%sfz)
    if (l_pushf) then
      call scale_em3d_splitfvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%fx,sf%fy,sf%fz, &
                             sf%sfx,sf%sfy,sf%sfz)
      call scale_em3d_splitefvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                             sf%exx,sf%eyy,sf%ezz, &
                             sf%sgx,sf%sgy,sf%sgz)
    end if
  end if
  return
end subroutine scale_em3d_split_fields


subroutine push_em3d_splitfvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,fx,fy,fz, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz,l_2dxz,l_2drz, &
                               xmin,ymin,zmin,dx,dy,dz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz
real(kind=8), intent(in) :: xmin,ymin,zmin,dx,dy,dz

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz,l_2drz
real(8) :: ru,rd

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
  if (l_2drz) then
    if (xmin==0.) then
      do l = 0, nz
        j = 0
        fx(j,k,l) = afx(j)*fx(j,k,l) + 4*bpfx(j)*(exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) 
        do j = 1, nx
          ru = 1.+0.5/j
          rd = 1.-0.5/j
          fx(j,k,l) = afx(j)*fx(j,k,l) + ru*bpfx(j)*(exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                       + rd*bmfx(j)*(exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
        end do
      end do
    else
      do l = 0, nz
        do j = 0, nx
          ru = (xmin+j*dx+0.5*dx)/(xmin+j*dx)
          rd = (xmin+j*dx-0.5*dx)/(xmin+j*dx)
          fx(j,k,l) = afx(j)*fx(j,k,l) + ru*bpfx(j)*(exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                       + rd*bmfx(j)*(exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
        end do
      end do
    end if
  else
    do l = 0, nz
      do j = 0, nx
       fx(j,k,l) = afx(j)*fx(j,k,l) + bpfx(j)*(exx(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
                                    + bmfx(j)*(exx(j-1,k  ,l  )+exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
      end do
    end do
  end if

  do l = 0, nz
    do j = 0, nx
     fz(j,k,l) = afz(l)*fz(j,k,l) + bpfz(l)*(ezx(j  ,k  ,l  )+ezz(j  ,k  ,l  )) &
                                  + bmfz(l)*(ezx(j  ,k  ,l-1)+ezz(j  ,k  ,l-1)) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
  end do

end if

  return
end subroutine push_em3d_splitfvec

subroutine scale_em3d_splitfvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               fx,fy,fz,&
                               sfx,sfy,sfz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: fx, &
                                                                                                       fy, &
                                                                                                       fz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: sfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: sfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: sfz

INTEGER :: j, k, l

  do l = -nzguard, nz+nzguard
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard
      fx(j,k,l) = sfx(j)*fx(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard
      fy(j,k,l) = sfy(k)*fy(j,k,l) 
    end do
   end do
  end do

  do l = -nzguard, nz+nzguard
   do k = -nyguard, ny+nyguard
    do j = -nxguard, nx+nxguard
      fz(j,k,l) = sfz(l)*fz(j,k,l) 
    end do
   end do
  end do

end subroutine scale_em3d_splitfvec

subroutine push_em3d_splitf_setcond(nx,ny,nz,nxcond,nycond,nzcond,nxguard,nyguard,nzguard, &
                               fx,fy,fz, &
                               l_2dxz,incond)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxcond,nycond,nzcond,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: fx,fy,fz
logical(ISZ), dimension(-nxguard:nxcond+nxguard,-nyguard:nycond+nyguard,-nzguard:nzcond+nzguard), intent(inout) :: incond

INTEGER :: j, k, l
logical(ISZ) :: l_2dxz

if (.not.l_2dxz) then

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
     if (incond(j,k,l)) then 
       fx(j,k,l) = 0. 
       fy(j,k,l) = 0. 
       fz(j,k,l) = 0. 
      end if
    end do
   end do
  end do

else
  k = 0

    do l = 0, nz
      do j = 0, nx
        if (incond(j,k,l)) then 
          fx(j,k,l) = 0. 
          fz(j,k,l) = 0. 
        end if
      end do
    end do

end if

  return
end subroutine push_em3d_splitf_setcond

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
REAL(kind=8) :: dt

INTEGER :: j, k, l

  call push_em3d_f(b%core%yf,dt*0.5)
  call push_em3d_b(b%core%yf,dt*0.5,2)
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
  call push_em3d_b(b%core%yf,dt*0.5,1)
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
!  call em3d_exchange_e(b)

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
    call push_em3d_b(b%core%yf,dt,which)
  else
    if(l_pushf) call push_em3d_f(b%core%yf,dt*0.5)
    call push_em3d_b(b%core%yf,dt*0.5,which)
  endif

  if(l_pushf) call push_em3d_blockbndf(b,dt,which)
  call push_em3d_blockbndb(b,dt,which)

!  if(l_pushf) call em3d_exchange_f(b)
  call em3d_applybc_b(b%core%yf, &
                      b%xlbnd, &
                      b%xrbnd, &
                      b%ylbnd, &
                      b%yrbnd, &
                      b%zlbnd, &
                      b%zrbnd)
!  call em3d_exchange_b(b)

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

subroutine scale_em3d_bnd_fields(b,dt,l_pushf)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
logical(ISZ) :: l_pushf

!if (b%core%yf%spectral) return

  if(b%xlbnd==openbc) call scale_em3d_split_fields(b%sidexl%syf,dt,l_pushf)
  if(b%xrbnd==openbc) call scale_em3d_split_fields(b%sidexr%syf,dt,l_pushf)
  if(b%ylbnd==openbc) call scale_em3d_split_fields(b%sideyl%syf,dt,l_pushf)
  if(b%yrbnd==openbc) call scale_em3d_split_fields(b%sideyr%syf,dt,l_pushf)
  if(b%zlbnd==openbc) call scale_em3d_split_fields(b%sidezl%syf,dt,l_pushf)
  if(b%zrbnd==openbc) call scale_em3d_split_fields(b%sidezr%syf,dt,l_pushf)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call scale_em3d_split_fields(b%edgexlyl%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call scale_em3d_split_fields(b%edgexryl%syf,dt,l_pushf)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call scale_em3d_split_fields(b%edgexlyr%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call scale_em3d_split_fields(b%edgexryr%syf,dt,l_pushf)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%edgexlzl%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%edgexrzl%syf,dt,l_pushf)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%edgexlzr%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%edgexrzr%syf,dt,l_pushf)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%edgeylzl%syf,dt,l_pushf)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%edgeyrzl%syf,dt,l_pushf)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%edgeylzr%syf,dt,l_pushf)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%edgeyrzr%syf,dt,l_pushf)

  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%cornerxlylzl%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%cornerxrylzl%syf,dt,l_pushf)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%cornerxlyrzl%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) call scale_em3d_split_fields(b%cornerxryrzl%syf,dt,l_pushf)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%cornerxlylzr%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%cornerxrylzr%syf,dt,l_pushf)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%cornerxlyrzr%syf,dt,l_pushf)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) call scale_em3d_split_fields(b%cornerxryrzr%syf,dt,l_pushf)

end subroutine scale_em3d_bnd_fields

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
  if ( f%l_2drz .and. f%circ_m > 0 ) then
       call shift_circarray_ncells_z(f%Ex_circ,f%nx,f%nz,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
       call shift_circarray_ncells_z(f%Ey_circ,f%nx,f%nz,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
       call shift_circarray_ncells_z(f%Ez_circ,f%nx,f%nz,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
       call shift_circarray_ncells_z(f%Bx_circ,f%nx,f%nz,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
       call shift_circarray_ncells_z(f%By_circ,f%nx,f%nz,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
       call shift_circarray_ncells_z(f%Bz_circ,f%nx,f%nz,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
   endif
!  do it=1,f%ntimes
!    do i=1,3
!      call shift_3darray_ncells_z(f%Jarray(-f%nxguard,-f%nyguard,-f%nzguard,i,it), &
!                                    f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    end do
!  end do
   if (f%nxf>0) then
      call shift_3darray_ncells_z(f%F,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
      if ( f%l_2drz .and. f%circ_m > 0 ) then 
         call shift_circarray_ncells_z(f%F_circ,f%nxf,f%nzf,f%circ_m,f%nxguard,f%nzguard,zl,zr,n)
      endif
         !    call shift_3darray_ncells_z(f%Rhoold,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    call shift_3darray_ncells_z(f%Rho,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    do it=1,f%ntimes
!      call shift_3darray_ncells_z(f%Rhoarray(-f%nxguard,-f%nyguard,-f%nzguard,it), &
!                                    f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
!    end do
  end if
  if (f%nxcond>0) then
    call shift_3dlarray_ncells_z(f%incond,f%nxcond,f%nycond,f%nzcond,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
  end if
  if (f%nxpnext>0) then
    call shift_3darray_ncells_z(f%Expnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Eypnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Ezpnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Bxpnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Bypnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
    call shift_3darray_ncells_z(f%Bzpnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,zl,zr,n)
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
integer(ISZ), parameter:: otherproc=10, ibuf = 950
real(kind=8) :: f(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) 

if (n > 0) then
  f(:,:,-nzguard:nz+nzguard-n) = f(:,:,-nzguard+n:nz+nzguard)
  if (zr/=otherproc) f(:,:,nz+nzguard-n+1:) = 0.

#ifdef MPIPARALLEL
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,:,nzguard-n+1:nzguard)),ibuf)
     call mympi_pack(f(:,:,nzguard-n+1:nzguard),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,:,nz+nzguard-n+1:)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,:,nz+nzguard-n+1:) = reshape(mpi_unpack_real_array( size(f(:,:,nz+nzguard-n+1:)),ibuf), &
                                                           shape(f(:,:,nz+nzguard-n+1:)))
  end if
#endif

else if (n < 0) then
  f(:,:,-nzguard-n:nz+nzguard) = f(:,:,-nzguard:nz+nzguard+n)
  if (zr/=otherproc) f(:,:,:-nzguard-n-1) = 0.

#ifdef MPIPARALLEL
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,:,nz-nzguard:nz-nzguard-n-1)),ibuf)
     call mympi_pack(f(:,:,nz-nzguard:nz-nzguard-n-1),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,:,-nzguard:-nzguard-n-1)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,:,-nzguard:-nzguard-n-1) = reshape(mpi_unpack_real_array( size(f(:,:,-nzguard:-nzguard-n-1)),ibuf), &
                                                                 shape(f(:,:,-nzguard:-nzguard-n-1)))
  end if
#endif

endif

  return
end subroutine shift_3darray_ncells_z

subroutine shift_3dlarray_ncells_z(f,nx,ny,nz,nxguard,nyguard,nzguard,zl,zr,n)
#ifdef MPIPARALLEL
use mpirz
#endif
implicit none
integer(ISZ) :: nx,ny,nz,nxguard,nyguard,nzguard,n,zl,zr
integer(ISZ), parameter:: otherproc=10, ibuf = 950
logical(ISZ) :: f(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) 

if (n > 0) then
  f(:,:,-nzguard:nz+nzguard-n) = f(:,:,-nzguard+n:nz+nzguard)
  if (zr/=otherproc) f(:,:,nz+nzguard-n+1:) = .false.

#ifdef MPIPARALLEL
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,:,nzguard-n+1:nzguard)),ibuf)
     call mympi_pack(f(:,:,nzguard-n+1:nzguard),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,:,nz+nzguard-n+1:)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,:,nz+nzguard-n+1:) = reshape(mpi_unpack_logical_array( size(f(:,:,nz+nzguard-n+1:)),ibuf), &
                                                              shape(f(:,:,nz+nzguard-n+1:)))
  end if
#endif

else if (n < 0) then
  f(:,:,-nzguard-n:nz+nzguard) = f(:,:,-nzguard:nz+nzguard+n)
  if (zr/=otherproc) f(:,:,:-nzguard-n-1) = .false.

#ifdef MPIPARALLEL
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,:,nz-nzguard:nz-nzguard-n-1)),ibuf)
     call mympi_pack(f(:,:,nz-nzguard:nz-nzguard-n-1),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,:,-nzguard:-nzguard-n-1)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,:,-nzguard:-nzguard-n-1) = reshape(mpi_unpack_logical_array( size(f(:,:,-nzguard:-nzguard-n-1)),ibuf), &
                                                                    shape(f(:,:,-nzguard:-nzguard-n-1)))
  end if
#endif

endif

  return
end subroutine shift_3dlarray_ncells_z

subroutine shift_circarray_ncells_z(f,nx,nz,circ_m,nxguard,nzguard,zl,zr,n)
#ifdef MPIPARALLEL
use mpirz
#endif
implicit none
integer(ISZ) :: nx,nz,circ_m,nxguard,nzguard,n,zl,zr
integer(ISZ), parameter:: otherproc=10, ibuf = 950
complex(kind=8) :: f(-nxguard:nx+nxguard,-nzguard:nz+nzguard,1:circ_m)
complex(kind=8) :: i=(0.,1.)

if (n > 0) then
  f(:,-nzguard:nz+nzguard-n,:) = f(:,-nzguard+n:nz+nzguard,:)
  if (zr/=otherproc) f(:,nz+nzguard-n+1:,:) = 0.

#ifdef MPIPARALLEL
  ! Send and receive the real part of the array
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,nzguard-n+1:nzguard,:)),ibuf)
     call mympi_pack(dble(realpart(f(:,nzguard-n+1:nzguard,:))),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if    
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,nz+nzguard-n+1:,:)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,nz+nzguard-n+1:,:) = reshape(mpi_unpack_real_array( size(f(:,nz+nzguard-n+1:,:)),ibuf), &
                                                           shape(f(:,nz+nzguard-n+1:,:)))
 end if
 ! Send and receive the imaginary part of the array
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,nzguard-n+1:nzguard,:)),ibuf)
     call mympi_pack(dble(imagpart(f(:,nzguard-n+1:nzguard,:))),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if    
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,nz+nzguard-n+1:,:)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,nz+nzguard-n+1:,:) = f(:,nz+nzguard-n+1:,:) + i*reshape(mpi_unpack_real_array( &
         size(f(:,nz+nzguard-n+1:,:)),ibuf), shape(f(:,nz+nzguard-n+1:,:)))
 end if
#endif

else if (n < 0) then
  f(:,-nzguard-n:nz+nzguard,:) = f(:,-nzguard:nz+nzguard+n,:)
  if (zr/=otherproc) f(:,:-nzguard-n-1,:) = 0.

#ifdef MPIPARALLEL
  ! Send and receive the real part of the array
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,nz-nzguard:nz-nzguard-n-1,:)),ibuf)
     call mympi_pack(dble(realpart(f(:,nz-nzguard:nz-nzguard-n-1,:))),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,-nzguard:-nzguard-n-1,:)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,-nzguard:-nzguard-n-1,:) = reshape(mpi_unpack_real_array( size(f(:,-nzguard:-nzguard-n-1,:)),ibuf), &
                                                                    shape(f(:,-nzguard:-nzguard-n-1,:)))
  end if
 ! Send and receive the imaginary part of the array
  if (zl==otherproc) then
     call mpi_packbuffer_init(size(f(:,nz-nzguard:nz-nzguard-n-1,:)),ibuf)
     call mympi_pack(dble(imagpart(f(:,nz-nzguard:nz-nzguard-n-1,:))),ibuf)
     call mpi_send_pack(procneighbors(0,2),0,ibuf)
  end if
  if (zr==otherproc) then
    call mpi_packbuffer_init(size(f(:,-nzguard:-nzguard-n-1,:)),ibuf)
    call mpi_recv_pack(procneighbors(1,2),0,ibuf)
    f(:,-nzguard:-nzguard-n-1,:) = f(:,-nzguard:-nzguard-n-1,:) + i*reshape(mpi_unpack_real_array( &
         size(f(:,-nzguard:-nzguard-n-1,:)),ibuf), shape(f(:,-nzguard:-nzguard-n-1,:)))
  end if
#endif

endif

  return
end subroutine shift_circarray_ncells_z

subroutine shift_em3dblock_ncells_x(b,n)
use mod_emfield3d
implicit none
TYPE(EM3D_BLOCKtype) :: b
integer(ISZ):: n

  ! --- shift core 
  call shift_em3df_ncells_x(b%core%yf,b%xlbnd,b%xrbnd,n)
  ! --- shift sides
  if(b%xlbnd==openbc) call shift_em3dsplitf_ncells_x(b%sidexl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc) call shift_em3dsplitf_ncells_x(b%sidexr%syf,openbc,openbc,n)
  if(b%ylbnd==openbc) call shift_em3dsplitf_ncells_x(b%sideyl%syf,b%xlbnd,b%xrbnd,n)
  if(b%yrbnd==openbc) call shift_em3dsplitf_ncells_x(b%sideyr%syf,b%xlbnd,b%xrbnd,n)
  if(b%zlbnd==openbc) call shift_em3dsplitf_ncells_x(b%sidezl%syf,b%xlbnd,b%xrbnd,n)
  if(b%zrbnd==openbc) call shift_em3dsplitf_ncells_x(b%sidezr%syf,b%xlbnd,b%xrbnd,n)
  ! --- shift edges
  if(b%xlbnd==openbc .and. b%ylbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexlyl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexryl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexlyr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexryr%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexlzl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexrzl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexlzr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgexrzr%syf,openbc,openbc,n)
  if(b%ylbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgeylzl%syf,b%xlbnd,b%xrbnd,n)
  if(b%yrbnd==openbc .and. b%zlbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgeyrzl%syf,b%xlbnd,b%xrbnd,n)
  if(b%ylbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgeylzr%syf,b%xlbnd,b%xrbnd,n)
  if(b%yrbnd==openbc .and. b%zrbnd==openbc) call shift_em3dsplitf_ncells_x(b%edgeyrzr%syf,b%xlbnd,b%xrbnd,n)
  ! --- shift corners
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxlylzl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxrylzl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxlyrzl%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zlbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxryrzl%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxlylzr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%ylbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxrylzr%syf,openbc,openbc,n)
  if(b%xlbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxlyrzr%syf,openbc,openbc,n)
  if(b%xrbnd==openbc .and. b%yrbnd==openbc .and. b%zrbnd==openbc) &
  call shift_em3dsplitf_ncells_x(b%cornerxryrzr%syf,openbc,openbc,n)
  
  return
end subroutine shift_em3dblock_ncells_x

subroutine shift_em3df_ncells_x(f,xl,xr,n)
use mod_emfield3d
implicit none
TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ):: n,i,it,xl,xr
  call shift_3darray_ncells_x(f%Ex,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Ey,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Ez,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Bx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%By,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Bz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  do it=1,f%ntimes
    call shift_3darray_ncells_x(f%Jarray(:,:,:,1,it), &
                                    f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Jarray(:,:,:,2,it), &
                                    f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Jarray(:,:,:,3,it), &
                                    f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  end do
  if (f%nxf>0) then
    call shift_3darray_ncells_x(f%F,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Rhoold,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
!    call shift_3darray_ncells_x(f%Rho,f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    do it=1,f%ntimes
      call shift_3darray_ncells_x(f%Rhoarray(:,:,:,it), &
                                    f%nxf,f%nyf,f%nzf,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    end do
  end if
  if (f%nxpnext>0) then
    call shift_3darray_ncells_x(f%Expnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Eypnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Ezpnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Bxpnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Bypnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Bzpnext,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  end if
  if (f%nxpo>0) then
    call shift_3darray_ncells_x(f%Ax,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Ay,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Az,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Phi,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  end if
  
  return
end subroutine shift_em3df_ncells_x

subroutine shift_em3dsplitf_ncells_x(f,xl,xr,n)
use mod_emfield3d
implicit none
TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ):: n,xl,xr

  call shift_3darray_ncells_x(f%Exx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Exy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Exz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)

  call shift_3darray_ncells_x(f%Eyx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Eyy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Eyz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)

  call shift_3darray_ncells_x(f%Ezx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Ezy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Ezz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)

  call shift_3darray_ncells_x(f%Bxy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Bxz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)

  call shift_3darray_ncells_x(f%Byx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Byz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)

  call shift_3darray_ncells_x(f%Bzx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Bzy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)

  call shift_3darray_ncells_x(f%Fx,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Fy,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  call shift_3darray_ncells_x(f%Fz,f%nx,f%ny,f%nz,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  if (f%nxpo>0) then
    call shift_3darray_ncells_x(f%Ax,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Ay,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Az,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
    call shift_3darray_ncells_x(f%Phi,f%nxpo,f%nypo,f%nzpo,f%nxguard,f%nyguard,f%nzguard,xl,xr,n)
  end if

  return
end subroutine shift_em3dsplitf_ncells_x

subroutine shift_3darray_ncells_x(f,nx,ny,nz,nxguard,nyguard,nzguard,xl,xr,n)
#ifdef MPIPARALLEL
use mpirz
#endif
implicit none
integer(ISZ) :: nx,ny,nz,nxguard,nyguard,nzguard,n,xl,xr
integer(ISZ), parameter:: otherproc=10, ibuf = 950
real(kind=8) :: f(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) 

if (n>0) then
  f(-nxguard:nx+nxguard-n,:,:) = f(-nxguard+n:nx+nxguard,:,:)
  if (xr/=otherproc) f(nx+nxguard-n+1:,:,:) = 0.

#ifdef MPIPARALLEL
    if (xl==otherproc) then
       call mpi_packbuffer_init(size(f(nxguard-n+1:nxguard,:,:)),ibuf)
       call mympi_pack(f(nxguard-n+1:nxguard,:,:),ibuf)
       call mpi_send_pack(procneighbors(0,0),0,ibuf)
    end if    
    if (xr==otherproc) then
      call mpi_packbuffer_init(size(f(nx+nxguard-n+1:nx+nxguard,:,:)),ibuf)
      call mpi_recv_pack(procneighbors(1,0),0,ibuf)
      f(nx+nxguard-n+1:nx+nxguard,:,:) = reshape(mpi_unpack_real_array( size(f(nx+nxguard-n+1:nx+nxguard,:,:)),ibuf), &
                                                           shape(f(nx+nxguard-n+1:nx+nxguard,:,:)))
    end if
#endif
end if
if (n<0) then
  n=-n
  f(-nxguard+n:nx+nxguard,:,:) = f(-nxguard:nx+nxguard-n,:,:) 
  if (xl/=otherproc) f(-nxguard:-nxguard+n-1,:,:) = 0.

#ifdef MPIPARALLEL
    if (xr==otherproc) then
       call mpi_packbuffer_init(size(f(nx-nxguard:nx-nxguard+n,:,:)),ibuf)
       call mympi_pack(f(nx-nxguard:nx-nxguard+n,:,:),ibuf)
       call mpi_send_pack(procneighbors(1,0),0,ibuf)
    end if    
    if (xl==otherproc) then
      call mpi_packbuffer_init(size(f(:-nxguard+n-1,:,:)),ibuf)
      call mpi_recv_pack(procneighbors(0,0),0,ibuf)
      f(:-nxguard+n-1,:,:) = reshape(mpi_unpack_real_array( size(f(:-nxguard+n-1,:,:)),ibuf), &
                                                                   shape(f(:-nxguard+n-1,:,:)))
    end if
#endif
  n=-n
end if

  return
end subroutine shift_3darray_ncells_x

subroutine em3d_applybc_e(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd,ifact,m

  if (f%l_2drz .and. f%xmin==0.) then
     f%ex(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%ex(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     f%ey(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%ey(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     f%ez(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%ez(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     if (f%circ_m>0) then
        do m = 1 , f%circ_m
           if (mod(m,2)==0) then
              ifact=1   ! Modes with even index have even symmetry
           else
              ifact=-1  ! Modes with odd index have odd symmetry
           end if
           f%ez_circ(f%ixmin-f%nxguard:f%ixmin-1,:,m) = &
                ifact*f%ez_circ(f%ixmin+f%nxguard:f%ixmin+1:-1,:,m)
           ! Although the Er and Etheta fields (represented by ex_circ and ey_circ) have
           ! the same symmetry as ez_circ, one has to additionally take into account the
           ! factors cos(theta) and sin(theta) which are added when calculating the fields
           ! Ex and Ey on the macroparticles.
           ! These factors have an odd symmetry, which has to be taken into account **here**,
           ! since, when interpolating the fields onto the macroparticles, **theta is not
           ! changed into theta + pi when the guard cells (with r<0) are used**.
           f%ex_circ(f%ixmin-f%nxguard:f%ixmin-1,:,m) = &
                -ifact*f%ex_circ(f%ixmin+f%nxguard-1:f%ixmin:-1,:,m)
           f%ey_circ(f%ixmin-f%nxguard:f%ixmin-1,:,m) = &
                -ifact*f%ey_circ(f%ixmin+f%nxguard:f%ixmin+1:-1,:,m)
       end do 
    end if
  end if

  if (xlbnd==dirichlet) then
     f%ey(f%ixmin,:,:) = 0.
     f%ez(f%ixmin,:,:) = 0.
     f%ex(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%ex(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     f%ey(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%ey(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     f%ez(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%ez(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
  end if

  if (xrbnd==dirichlet) then
     f%ey(f%ixmax,:,:) = 0.
     f%ez(f%ixmax,:,:) = 0.
     f%ex(f%ixmax:f%ixmax+f%nxguard,:,:)   =  f%ex(f%ixmax-1:f%ixmax-f%nxguard-1:-1,:,:)
     f%ey(f%ixmax+1:f%ixmax+f%nxguard,:,:) = -f%ey(f%ixmax-1:f%ixmax-f%nxguard:-1,:,:)
     f%ez(f%ixmax+1:f%ixmax+f%nxguard,:,:) = -f%ez(f%ixmax-1:f%ixmax-f%nxguard:-1,:,:)
  end if

  if (ylbnd==dirichlet) then
     f%ex(:,f%iymin,:) = 0.
     f%ez(:,f%iymin,:) = 0.
     f%ex(:,f%iymin-f%nyguard:f%iymin-1,:) = -f%ex(:,f%iymin+f%nyguard:f%iymin+1:-1,:)
     f%ey(:,f%iymin-f%nyguard:f%iymin-1,:) =  f%ey(:,f%iymin+f%nyguard-1:f%iymin:-1,:)
     f%ez(:,f%iymin-f%nyguard:f%iymin-1,:) = -f%ez(:,f%iymin+f%nyguard:f%iymin+1:-1,:)
  end if

  if (yrbnd==dirichlet) then
     f%ex(:,f%iymax,:) = 0.
     f%ez(:,f%iymax,:) = 0.
     f%ex(:,f%iymax+1:f%iymax+f%nyguard,:) = -f%ex(:,f%iymax-1:f%iymax-f%nyguard:-1,:)
     f%ey(:,f%iymax:f%iymax+f%nyguard,:)   =  f%ey(:,f%iymax-1:f%iymax-f%nyguard-1:-1,:)
     f%ez(:,f%iymax+1:f%iymax+f%nyguard,:) = -f%ez(:,f%iymax-1:f%iymax-f%nyguard:-1,:)
  end if

  if (zlbnd==dirichlet) then
     f%ex(:,:,f%izmin) = 0.
     f%ey(:,:,f%izmin) = 0.
     f%ex(:,:,f%izmin-f%nzguard:f%izmin-1) = -f%ex(:,:,f%izmin+f%nzguard:f%izmin+1:-1)
     f%ey(:,:,f%izmin-f%nzguard:f%izmin-1) = -f%ey(:,:,f%izmin+f%nzguard:f%izmin+1:-1)
     f%ez(:,:,f%izmin-f%nzguard:f%izmin-1) =  f%ez(:,:,f%izmin+f%nzguard-1:f%izmin:-1)
  end if

  if (zrbnd==dirichlet) then
     f%ex(:,:,f%izmax) = 0.
     f%ey(:,:,f%izmax) = 0.
     f%ex(:,:,f%izmax+1:f%izmax+f%nzguard) = -f%ex(:,:,f%izmax-1:f%izmax-f%nzguard:-1)
     f%ey(:,:,f%izmax+1:f%izmax+f%nzguard) = -f%ey(:,:,f%izmax-1:f%izmax-f%nzguard:-1)
     f%ez(:,:,f%izmax:f%izmax+f%nzguard)   =  f%ez(:,:,f%izmax-1:f%izmax-f%nzguard-1:-1)
  end if

  if (xlbnd==neumann) then
     f%ex(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%ex(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     f%ey(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%ey(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     f%ez(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%ez(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
  end if

  if (xrbnd==neumann) then
     f%ex(f%ixmax:f%ixmax+f%nxguard,:,:)   = -f%ex(f%ixmax-1:f%ixmax-f%nxguard-1:-1,:,:)
     f%ey(f%ixmax+1:f%ixmax+f%nxguard,:,:) =  f%ey(f%ixmax-1:f%ixmax-f%nxguard:-1,:,:)
     f%ez(f%ixmax+1:f%ixmax+f%nxguard,:,:) =  f%ez(f%ixmax-1:f%ixmax-f%nxguard:-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%ex(:,f%iymin-f%nyguard:f%iymin-1,:) =  f%ex(:,f%iymin+f%nyguard:f%iymin+1:-1,:)
     f%ey(:,f%iymin-f%nyguard:f%iymin-1,:) = -f%ey(:,f%iymin+f%nyguard-1:f%iymin:-1,:)
     f%ez(:,f%iymin-f%nyguard:f%iymin-1,:) =  f%ez(:,f%iymin+f%nyguard:f%iymin+1:-1,:)
  end if

  if (yrbnd==neumann) then
     f%ex(:,f%iymax+1:f%iymax+f%nyguard,:) =  f%ex(:,f%iymax-1:f%iymax-f%nyguard:-1,:)
     f%ey(:,f%iymax:f%iymax+f%nyguard,:)   = -f%ey(:,f%iymax-1:f%iymax-f%nyguard-1:-1,:)
     f%ez(:,f%iymax+1:f%iymax+f%nyguard,:) =  f%ez(:,f%iymax-1:f%iymax-f%nyguard:-1,:)
  end if

  if (zlbnd==neumann) then
     f%ex(:,:,f%izmin-f%nzguard:f%izmin-1) =  f%ex(:,:,f%izmin+f%nzguard:f%izmin+1:-1)
     f%ey(:,:,f%izmin-f%nzguard:f%izmin-1) =  f%ey(:,:,f%izmin+f%nzguard:f%izmin+1:-1)
     f%ez(:,:,f%izmin-f%nzguard:f%izmin-1) = -f%ez(:,:,f%izmin+f%nzguard-1:f%izmin:-1)
  end if

  if (zrbnd==neumann) then
     f%ex(:,:,f%izmax+1:f%izmax+f%nzguard) =  f%ex(:,:,f%izmax-1:f%izmax-f%nzguard:-1)
     f%ey(:,:,f%izmax+1:f%izmax+f%nzguard) =  f%ey(:,:,f%izmax-1:f%izmax-f%nzguard:-1)
     f%ez(:,:,f%izmax:f%izmax+f%nzguard)   = -f%ez(:,:,f%izmax-1:f%izmax-f%nzguard-1:-1)
  end if

  return
end subroutine em3d_applybc_e

subroutine em3d_applybc_b(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd,ifact,m

  if (f%l_2drz .and. f%xmin==0.) then
     f%bx(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%bx(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     f%by(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%by(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     f%bz(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%bz(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     if (f%circ_m>0) then
        do m = 1 , f%circ_m
           if (mod(m,2)==0) then 
              ifact=1
           else
              ifact=-1
           end if
           ! Although the Br and Btheta fields (represented by bx_circ and by_circ) have
           ! the same symmetry as bz_circ, one has to additionally take into account the
           ! factors cos(theta) and sin(theta) which are added when calculating the fields
           ! Bx and By on the macroparticles.
           ! These factors have an odd symmetry, which has to be taken into account **here**,
           ! since, when interpolating the fields onto the macroparticles, **theta is not
           ! changed into theta + pi when the guard cells (with r<0) are used**.
           f%bx_circ(f%ixmin-f%nxguard:f%ixmin-1,:,m) = &
                -ifact*f%bx_circ(f%ixmin+f%nxguard:f%ixmin+1:-1,:,m)
           f%by_circ(f%ixmin-f%nxguard:f%ixmin-1,:,m) = &
                -ifact*f%by_circ(f%ixmin+f%nxguard-1:f%ixmin:-1,:,m)
           f%bz_circ(f%ixmin-f%nxguard:f%ixmin-1,:,m) = &
                ifact*f%bz_circ(f%ixmin+f%nxguard-1:f%ixmin:-1,:,m)
        end do
     end if
  end if

  if (xlbnd==dirichlet) then
     f%bx(f%ixmin,:,:) = 0.
     f%bx(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%bx(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     f%by(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%by(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     f%bz(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%bz(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
  end if

  if (xrbnd==dirichlet) then
     f%bx(f%ixmax,:,:) = 0.
     f%bx(f%ixmax+1:f%ixmax+f%nxguard,:,:) = -f%bx(f%ixmax-1:f%ixmax-f%nxguard:-1,:,:)
     f%by(f%ixmax:f%ixmax+f%nxguard,:,:)   =  f%by(f%ixmax-1:f%ixmax-f%nxguard-1:-1,:,:)
     f%bz(f%ixmax:f%ixmax+f%nxguard,:,:)   =  f%bz(f%ixmax-1:f%ixmax-f%nxguard-1:-1,:,:)
  end if

  if (ylbnd==dirichlet) then
     f%by(:,f%iymin,:) = 0.
     f%bx(:,f%iymin-f%nyguard:f%iymin-1,:) =  f%bx(:,f%iymin+f%nyguard-1:f%iymin:-1,:)
     f%by(:,f%iymin-f%nyguard:f%iymin-1,:) = -f%by(:,f%iymin+f%nyguard:f%iymin+1:-1,:)
     f%bz(:,f%iymin-f%nyguard:f%iymin-1,:) =  f%bz(:,f%iymin+f%nyguard-1:f%iymin:-1,:)
  end if

  if (yrbnd==dirichlet) then
     f%by(:,f%iymax,:) = 0.
     f%bx(:,f%iymax:f%iymax+f%nyguard,:)   =  f%bx(:,f%iymax-1:f%iymax-f%nyguard-1:-1,:)
     f%by(:,f%iymax+1:f%iymax+f%nyguard,:) = -f%by(:,f%iymax-1:f%iymax-f%nyguard:-1,:)
     f%bz(:,f%iymax:f%iymax+f%nyguard,:)   =  f%bz(:,f%iymax-1:f%iymax-f%nyguard-1:-1,:)
  end if

  if (zlbnd==dirichlet) then
     f%bz(:,:,f%izmin) = 0.
     f%bx(:,:,f%izmin-f%nzguard:f%izmin-1) =  f%bx(:,:,f%izmin+f%nzguard-1:f%izmin:-1)
     f%by(:,:,f%izmin-f%nzguard:f%izmin-1) =  f%by(:,:,f%izmin+f%nzguard-1:f%izmin:-1)
     f%bz(:,:,f%izmin-f%nzguard:f%izmin-1) = -f%bz(:,:,f%izmin+f%nzguard:f%izmin+1:-1)
  end if

  if (zrbnd==dirichlet) then
     f%bz(:,:,f%izmax) = 0.
     f%bx(:,:,f%izmax:f%izmax+f%nzguard)   =  f%bx(:,:,f%izmax-1:f%izmax-f%nzguard-1:-1)
     f%by(:,:,f%izmax:f%izmax+f%nzguard)   =  f%by(:,:,f%izmax-1:f%izmax-f%nzguard-1:-1)
     f%bz(:,:,f%izmax+1:f%izmax+f%nzguard) = -f%bz(:,:,f%izmax-1:f%izmax-f%nzguard:-1)
  end if

  if (xlbnd==neumann) then
     f%bx(f%ixmin-f%nxguard:f%ixmin-1,:,:) =  f%bx(f%ixmin+f%nxguard:f%ixmin+1:-1,:,:)
     f%by(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%by(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
     f%bz(f%ixmin-f%nxguard:f%ixmin-1,:,:) = -f%bz(f%ixmin+f%nxguard-1:f%ixmin:-1,:,:)
  end if

  if (xrbnd==neumann) then
     f%bx(f%ixmax+1:f%ixmax+f%nxguard,:,:) =  f%bx(f%ixmax-1:f%ixmax-f%nxguard:-1,:,:)
     f%by(f%ixmax:f%ixmax+f%nxguard,:,:)   = -f%by(f%ixmax-1:f%ixmax-f%nxguard-1:-1,:,:)
     f%bz(f%ixmax:f%ixmax+f%nxguard,:,:)   = -f%bz(f%ixmax-1:f%ixmax-f%nxguard-1:-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%bx(:,f%iymin-f%nyguard:f%iymin-1,:) = -f%bx(:,f%iymin+f%nyguard-1:f%iymin:-1,:)
     f%by(:,f%iymin-f%nyguard:f%iymin-1,:) =  f%by(:,f%iymin+f%nyguard:f%iymin+1:-1,:)
     f%bz(:,f%iymin-f%nyguard:f%iymin-1,:) = -f%bz(:,f%iymin+f%nyguard-1:f%iymin:-1,:)
  end if

  if (yrbnd==neumann) then
     f%bx(:,f%iymax:f%iymax+f%nyguard,:)   = -f%bx(:,f%iymax-1:f%iymax-f%nyguard-1:-1,:)
     f%by(:,f%iymax+1:f%iymax+f%nyguard,:) =  f%by(:,f%iymax-1:f%iymax-f%nyguard:-1,:)
     f%bz(:,f%iymax:f%iymax+f%nyguard,:)   = -f%bz(:,f%iymax-1:f%iymax-f%nyguard-1:-1,:)
  end if

  if (zlbnd==neumann) then
     f%bx(:,:,f%izmin-f%nzguard:f%izmin-1) = -f%bx(:,:,f%izmin+f%nzguard-1:f%izmin:-1)
     f%by(:,:,f%izmin-f%nzguard:f%izmin-1) = -f%by(:,:,f%izmin+f%nzguard-1:f%izmin:-1)
     f%bz(:,:,f%izmin-f%nzguard:f%izmin-1) =  f%bz(:,:,f%izmin+f%nzguard:f%izmin+1:-1)
  end if

  if (zrbnd==neumann) then
     f%bx(:,:,f%izmax:f%izmax+f%nzguard)   = -f%bx(:,:,f%izmax-1:f%izmax-f%nzguard-1:-1)
     f%by(:,:,f%izmax:f%izmax+f%nzguard)   = -f%by(:,:,f%izmax-1:f%izmax-f%nzguard-1:-1)
     f%bz(:,:,f%izmax+1:f%izmax+f%nzguard) =  f%bz(:,:,f%izmax-1:f%izmax-f%nzguard:-1)
  end if

  return
end subroutine em3d_applybc_b

subroutine em3d_applybc_j(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd,type_rz_depose)
  use mod_emfield3d
  use Constant
  implicit none
  
  TYPE(EM3D_YEEFIELDtype) :: f
  integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd,j,ifact, m,type_rz_depose
  real(8)::r,r1,r2
  complex(kind=8) :: I=(0.,1.)
  
  ! MODE 0 : Fetch the current deposited in the guards cells and add it to the grid (fold back)
  
  ! In rz geometry, for the guards cells below the axis
  if (f%l_2drz .and. f%xmin==0.) then
     ! Fields that are located on the boundary
     f%j(f%ixmin+1:f%ixmin+f%nxguard,:,:,2:3) = f%j(f%ixmin+1:f%ixmin+f%nxguard,:,:,2:3) &
          + f%j(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:,2:3)
     ! Fields that are located off the boundary 
     f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) = f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) &
          - f%j(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:,1)
  end if

  ! For guards at the lower and upper bound in x, in the Dirichlet case
  if (xlbnd==dirichlet) then
     f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) = f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) - f%j(f%ixmin:f%ixmin-f%nxguard:-1,:,:,2:3)
     f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) = f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) + f%j(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:,1)
  end if
  if (xrbnd==dirichlet) then
     f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) = f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) - f%j(f%ixmax+f%nxguard:f%ixmax:-1,:,:,2:3)
     f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) = f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) + f%j(f%ixmax+f%nxguard-1:f%ixmax:-1,:,:,1)
  end if

  ! For guards at the lower and upper bound in y, in the Dirichlet case
  if (ylbnd==dirichlet) then
     f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) = f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) &
                                              - f%j(:,f%iymin:f%iymin-f%nyguard:-1,:,1:3:2)
     f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) = f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) + f%j(:,f%iymin-1:f%iymin-f%nyguard:-1,:,2)
  end if
  if (yrbnd==dirichlet) then
     f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) = f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) &
                                              - f%j(:,f%iymax+f%nyguard:f%iymax:-1,:,1:3:2)
     f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) = f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) + f%j(:,f%iymax+f%nyguard-1:f%iymax:-1,:,2)
  end if

  ! For guards at the lower and upper bound in z, in the Dirichlet case
  if (zlbnd==dirichlet) then
     f%j(:,:,f%izmin:f%izmin+f%nzguard,1:2) = f%j(:,:,f%izmin:f%izmin+f%nzguard,1:2) - f%j(:,:,f%izmin:f%izmin-f%nzguard:-1,1:2)
     f%j(:,:,f%izmin:f%izmin+f%nzguard-1,3) = f%j(:,:,f%izmin:f%izmin+f%nzguard-1,3) + f%j(:,:,f%izmin-1:f%izmin-f%nzguard:-1,3)
  end if
  if (zrbnd==dirichlet) then
     f%j(:,:,f%izmax-f%nzguard:f%izmax,1:2) = f%j(:,:,f%izmax-f%nzguard:f%izmax,1:2) - f%j(:,:,f%izmax+f%nzguard:f%izmax:-1,1:2)
     f%j(:,:,f%izmax-f%nzguard:f%izmax-1,3) = f%j(:,:,f%izmax-f%nzguard:f%izmax-1,3) + f%j(:,:,f%izmax+f%nzguard-1:f%izmax:-1,3)
  end if

  ! For guards at the lower and upper bound in x, in the Neumann case
  if (xlbnd==neumann) then
     f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) = f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) + f%j(f%ixmin:f%ixmin-f%nxguard:-1,:,:,2:3)
     f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) = f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) - f%j(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:,1)
  end if
  if (xrbnd==neumann) then
     f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) = f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) + f%j(f%ixmax+f%nxguard:f%ixmax:-1,:,:,2:3)
     f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) = f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) - f%j(f%ixmax+f%nxguard-1:f%ixmax:-1,:,:,1)
  end if

  ! For guards at the lower and upper bound in y, in the Neumann case
  if (ylbnd==neumann) then
     f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) = f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) &
                                              + f%j(:,f%iymin:f%iymin-f%nyguard:-1,:,1:3:2)
     f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) = f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) - f%j(:,f%iymin-1:f%iymin-f%nyguard:-1,:,2)
  end if
  if (yrbnd==neumann) then
     f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) = f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) &
                                              + f%j(:,f%iymax+f%nyguard:f%iymax:-1,:,1:3:2)
     f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) = f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) - f%j(:,f%iymax+f%nyguard-1:f%iymax:-1,:,2)
  end if

  ! For guards at the lower and upper bound in z, in the Neumann case
  if (zlbnd==neumann) then
     f%j(:,:,f%izmin:f%izmin+f%nzguard,1:2) = f%j(:,:,f%izmin:f%izmin+f%nzguard,1:2) + f%j(:,:,f%izmin:f%izmin-f%nzguard:-1,1:2)
     f%j(:,:,f%izmin:f%izmin+f%nzguard-1,3) = f%j(:,:,f%izmin:f%izmin+f%nzguard-1,3) - f%j(:,:,f%izmin-1:f%izmin-f%nzguard:-1,3)
  end if
  if (zrbnd==neumann) then
     f%j(:,:,f%izmax-f%nzguard:f%izmax,1:2) = f%j(:,:,f%izmax-f%nzguard:f%izmax,1:2) + f%j(:,:,f%izmax+f%nzguard:f%izmax:-1,1:2)
     f%j(:,:,f%izmax-f%nzguard:f%izmax-1,3) = f%j(:,:,f%izmax-f%nzguard:f%izmax-1,3) - f%j(:,:,f%izmax+f%nzguard-1:f%izmax:-1,3)
  end if


  ! MODES > 0 : Fetch the current deposited in the guards cells and add it to the grid (fold back)
  if (f%circ_m>0) then

     if (f%l_2drz .and. f%xmin==0.) then
        do m = 1,f%circ_m
           if (mod(m,2)==0) then
              ifact=1
           else
              ifact=-1
           end if 
           f%j_circ(f%ixmin+1:f%ixmin+f%nxguard,:,2:3,m) = f%j_circ(f%ixmin+1:f%ixmin+f%nxguard,:,2:3,m) &
                + ifact*f%j_circ(f%ixmin-1:f%ixmin-f%nxguard:-1,:,2:3,m)
           f%j_circ(f%ixmin:f%ixmin+f%nxguard-1,:,1,m) = f%j_circ(f%ixmin:f%ixmin+f%nxguard-1,:,1,m) &
                - ifact*f%j_circ(f%ixmin-1:f%ixmin-f%nxguard:-1,:,1,m)
        end do
     end if

    if (xlbnd==dirichlet) then
     f%j_circ(f%ixmin:f%ixmin+f%nxguard,:,2:3,:) = f%j_circ(f%ixmin:f%ixmin+f%nxguard,:,2:3,:) &
                                                 - f%j_circ(f%ixmin:f%ixmin-f%nxguard:-1,:,2:3,:)
     f%j_circ(f%ixmin:f%ixmin+f%nxguard-1,:,1,:) = f%j_circ(f%ixmin:f%ixmin+f%nxguard-1,:,1,:) &
                                                 + f%j_circ(f%ixmin-1:f%ixmin-f%nxguard:-1,:,1,:)
    end if
    if (xrbnd==dirichlet) then
     f%j_circ(f%ixmax-f%nxguard:f%ixmax,:,2:3,:) = f%j_circ(f%ixmax-f%nxguard:f%ixmax,:,2:3,:) &
                                                 - f%j_circ(f%ixmax+f%nxguard:f%ixmax:-1,:,2:3,:)
     f%j_circ(f%ixmax-f%nxguard:f%ixmax-1,:,1,:) = f%j_circ(f%ixmax-f%nxguard:f%ixmax-1,:,1,:) &
                                                 + f%j_circ(f%ixmax+f%nxguard-1:f%ixmax:-1,:,1,:)
    end if

    if (zlbnd==dirichlet) then
     f%j_circ(:,f%izmin:f%izmin+f%nzguard,1:2,:) = f%j_circ(:,f%izmin:f%izmin+f%nzguard,1:2,:) &
                                                 - f%j_circ(:,f%izmin:f%izmin-f%nzguard:-1,1:2,:)
     f%j_circ(:,f%izmin:f%izmin+f%nzguard-1,3,:) = f%j_circ(:,f%izmin:f%izmin+f%nzguard-1,3,:) &
                                                 + f%j_circ(:,f%izmin-1:f%izmin-f%nzguard:-1,3,:)
    end if
    if (zrbnd==dirichlet) then
     f%j_circ(:,f%izmax-f%nzguard:f%izmax,1:2,:) = f%j_circ(:,f%izmax-f%nzguard:f%izmax,1:2,:) &
                                                 - f%j_circ(:,f%izmax+f%nzguard:f%izmax:-1,1:2,:)
     f%j_circ(:,f%izmax-f%nzguard:f%izmax-1,3,:) = f%j_circ(:,f%izmax-f%nzguard:f%izmax-1,3,:) &
                                                 + f%j_circ(:,f%izmax+f%nzguard-1:f%izmax:-1,3,:)
    end if

    if (xlbnd==neumann) then
     f%j_circ(f%ixmin:f%ixmin+f%nxguard,:,2:3,:) = f%j_circ(f%ixmin:f%ixmin+f%nxguard,:,2:3,:) &
                                                 + f%j_circ(f%ixmin:f%ixmin-f%nxguard:-1,:,2:3,:)
     f%j_circ(f%ixmin:f%ixmin+f%nxguard-1,:,1,:) = f%j_circ(f%ixmin:f%ixmin+f%nxguard-1,:,1,:) &
                                                 - f%j_circ(f%ixmin-1:f%ixmin-f%nxguard:-1,:,1,:)
    end if
    if (xrbnd==neumann) then
     f%j_circ(f%ixmax-f%nxguard:f%ixmax,:,2:3,:) = f%j_circ(f%ixmax-f%nxguard:f%ixmax,:,2:3,:) &
                                                 + f%j_circ(f%ixmax+f%nxguard:f%ixmax:-1,:,2:3,:)
     f%j_circ(f%ixmax-f%nxguard:f%ixmax-1,:,1,:) = f%j_circ(f%ixmax-f%nxguard:f%ixmax-1,:,1,:) &
                                                 - f%j_circ(f%ixmax+f%nxguard-1:f%ixmax:-1,:,1,:)
    end if

    if (zlbnd==neumann) then
     f%j_circ(:,f%izmin:f%izmin+f%nzguard,1:2,:) = f%j_circ(:,f%izmin:f%izmin+f%nzguard,1:2,:) &
                                                 + f%j_circ(:,f%izmin:f%izmin-f%nzguard:-1,1:2,:)
     f%j_circ(:,f%izmin:f%izmin+f%nzguard-1,3,:) = f%j_circ(:,f%izmin:f%izmin+f%nzguard-1,3,:) &
                                                 - f%j_circ(:,f%izmin-1:f%izmin-f%nzguard:-1,3,:)
    end if
    if (zrbnd==neumann) then
     f%j_circ(:,f%izmax-f%nzguard:f%izmax,1:2,:) = f%j_circ(:,f%izmax-f%nzguard:f%izmax,1:2,:) &
                                                 + f%j_circ(:,f%izmax+f%nzguard:f%izmax:-1,1:2,:)
     f%j_circ(:,f%izmax-f%nzguard:f%izmax-1,3,:) = f%j_circ(:,f%izmax-f%nzguard:f%izmax-1,3,:) &
                                                 - f%j_circ(:,f%izmax+f%nzguard-1:f%izmax:-1,3,:)
    end if
  end if

  ! In rz geometry, divide the current by the cell volume
  if (f%l_2drz) then

       ! -- Jr
      
       ! Since Jr is not cell centered in r, no need for distinction
       ! between on axis and off-axis factors
       do j=f%ixmin-f%nxguard,f%ixmax+f%nxguard
          r = abs(f%xmin+(float(j)+0.5)*f%dx)
          f%j(j,:,:,1) = f%j(j,:,:,1)/(2.*pi*r)
          if (f%circ_m>0) f%j_circ(j,:,1,:) = f%j_circ(j,:,1,:)/(2.*pi*r)
       end do

     ! -- Jtheta and Jz
     
     ! In the lower guard cells in x (exchanged with nearby processors)
     do j=f%ixmin-f%nxguard,f%ixmin-1
        r = abs(f%xmin+j*f%dx)
        f%j(j,:,:,2:3) = f%j(j,:,:,2:3)/(2.*pi*r)    ! Mode 0
        if (f%circ_m>0) f%j_circ(j,:,2:3,:) = f%j_circ(j,:,2:3,:)/(2.*pi*r)  ! Mode > 0
     end do
     
     ! On the lower boundary
     j = f%ixmin
     if (f%xmin==0.) then
        ! On axis
        ! Jz, mode 0
        if (type_rz_depose == 1) then ! Verboncoeur JCP 164, 421-427 (2001) : corrected volume
           f%j(j,:,:,3) = f%j(j,:,:,3)/(pi*f%dx/3.)
        else                          ! Standard volume
           f%j(j,:,:,3) = f%j(j,:,:,3)/(pi*f%dx/4.)
        endif
        ! Jz, modes > 0
        if (f%circ_m>0) f%j_circ(j,:,3,:) = 0. ! Mode > 0 : Jz is zero on axis.
        ! Jt, mode 0 and modes > 1 
        f%j(j,:,:,2)= 0. ! Mode 0 : Jt is zero on axis.
        if (f%circ_m>1) f%j_circ(j,:,2,2:) = 0. ! Modes > 1 : Jt = 0.
        ! Jt, mode 1
        if (f%circ_m>0) f%j_circ(j,:,2,1) = -I*f%j_circ(j,:,1,1)
        ! Because the previous line uses Jr, it is important that Jr be properly calculated first
     else
        ! Not the axis
        r = abs(f%xmin+j*f%dx)
        f%j(j,:,:,2:3) = f%j(j,:,:,2:3)/(2.*pi*r)    ! Mode 0
        if (f%circ_m>0) f%j_circ(j,:,2:3,:) = f%j_circ(j,:,2:3,:)/(2.*pi*r)  ! Mode > 0
     end if
     
     ! In the rest of the grid
     do j=f%ixmin+1,f%ixmax+f%nxguard
        r = abs(f%xmin+j*f%dx)
        f%j(j,:,:,2:3) = f%j(j,:,:,2:3)/(2.*pi*r)
        if (f%circ_m>0) f%j_circ(j,:,2:3,:) = f%j_circ(j,:,2:3,:)/(2.*pi*r)
     end do
  end if
  
  return
end subroutine em3d_applybc_j

subroutine em3d_applybc_rho(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd,type_rz_depose)
use mod_emfield3d
use Constant
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd,j,ifact,m,type_rz_depose
real(8)::r

  if (f%l_2drz .and. f%xmin==0.) then
     f%rho(f%ixmin+1:f%ixmin+f%nxguard,:,:) = f%rho(f%ixmin+1:f%ixmin+f%nxguard,:,:) + f%rho(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:)
  end if

  if (xlbnd==dirichlet) then
     f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) = f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) - f%rho(f%ixmin:f%ixmin-f%nxguard:-1,:,:)
  end if
  if (xrbnd==dirichlet) then
     f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) = f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) - f%rho(f%ixmax+f%nxguard:f%ixmax:-1,:,:)
  end if

  if (ylbnd==dirichlet) then
     f%rho(:,f%iymin:f%iymin+f%nyguard,:) = f%rho(:,f%iymin:f%iymin+f%nyguard,:) - f%rho(:,f%iymin:f%iymin-f%nyguard:-1,:)
  end if
  if (yrbnd==dirichlet) then
     f%rho(:,f%iymax-f%nyguard:f%iymax,:) = f%rho(:,f%iymax-f%nyguard:f%iymax,:) - f%rho(:,f%iymax+f%nyguard:f%iymax:-1,:)
  end if

  if (zlbnd==dirichlet) then
     f%rho(:,:,f%izmin:f%izmin+f%nzguard) = f%rho(:,:,f%izmin:f%izmin+f%nzguard) - f%rho(:,:,f%izmin:f%izmin-f%nzguard:-1)
  end if
  if (zrbnd==dirichlet) then
     f%rho(:,:,f%izmax-f%nzguard:f%izmax) = f%rho(:,:,f%izmax-f%nzguard:f%izmax) - f%rho(:,:,f%izmax+f%nzguard:f%izmax:-1)
  end if

  if (xlbnd==neumann) then
     f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) = f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) + f%rho(f%ixmin:f%ixmin-f%nxguard:-1,:,:)
  end if
  if (xrbnd==neumann) then
     f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) = f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) + f%rho(f%ixmax+f%nxguard:f%ixmax:-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%rho(:,f%iymin:f%iymin+f%nyguard,:) = f%rho(:,f%iymin:f%iymin+f%nyguard,:) + f%rho(:,f%iymin:f%iymin-f%nyguard:-1,:)
  end if
  if (yrbnd==neumann) then
     f%rho(:,f%iymax-f%nyguard:f%iymax,:) = f%rho(:,f%iymax-f%nyguard:f%iymax,:) + f%rho(:,f%iymax+f%nyguard:f%iymax:-1,:)
  end if

  if (zlbnd==neumann) then
     f%rho(:,:,f%izmin:f%izmin+f%nzguard) = f%rho(:,:,f%izmin:f%izmin+f%nzguard) + f%rho(:,:,f%izmin:f%izmin-f%nzguard:-1)
  end if
  if (zrbnd==neumann) then
     f%rho(:,:,f%izmax-f%nzguard:f%izmax) = f%rho(:,:,f%izmax-f%nzguard:f%izmax) + f%rho(:,:,f%izmax+f%nzguard:f%izmax:-1)
  end if

  if (f%circ_m>0) then
    if (f%l_2drz .and. f%xmin==0.) then
       do m= 1, f%circ_m
          if (mod(m,2)==0) then
             ifact=1
          else
             ifact=-1
          end if
          f%rho_circ(f%ixmin+1:f%ixmin+f%nxguard,:,m) = f%rho_circ(f%ixmin+1:f%ixmin+f%nxguard,:,m) &
                                                 + ifact*f%rho_circ(f%ixmin-1:f%ixmin-f%nxguard:-1,:,m)
       end do
    end if
    if (xlbnd==dirichlet) then
     f%rho_circ(f%ixmin:f%ixmin+f%nxguard,:,:) = f%rho_circ(f%ixmin:f%ixmin+f%nxguard,:,:) &
                                               - f%rho_circ(f%ixmin:f%ixmin-f%nxguard:-1,:,:)
    end if
    if (xrbnd==dirichlet) then
     f%rho_circ(f%ixmax-f%nxguard:f%ixmax,:,:) = f%rho_circ(f%ixmax-f%nxguard:f%ixmax,:,:) &
                                               - f%rho_circ(f%ixmax+f%nxguard:f%ixmax:-1,:,:)
    end if

    if (zlbnd==dirichlet) then
     f%rho_circ(:,f%izmin:f%izmin+f%nzguard,:) = f%rho_circ(:,f%izmin:f%izmin+f%nzguard,:) &
                                               - f%rho_circ(:,f%izmin:f%izmin-f%nzguard:-1,:)
    end if
    if (zrbnd==dirichlet) then
     f%rho_circ(:,f%izmax-f%nzguard:f%izmax,:) = f%rho_circ(:,f%izmax-f%nzguard:f%izmax,:) &
                                               - f%rho_circ(:,f%izmax+f%nzguard:f%izmax:-1,:)
    end if

    if (xlbnd==neumann) then
     f%rho_circ(f%ixmin:f%ixmin+f%nxguard,:,:) = f%rho_circ(f%ixmin:f%ixmin+f%nxguard,:,:) &
                                               + f%rho_circ(f%ixmin:f%ixmin-f%nxguard:-1,:,:)
    end if
    if (xrbnd==neumann) then
     f%rho_circ(f%ixmax-f%nxguard:f%ixmax,:,:) = f%rho_circ(f%ixmax-f%nxguard:f%ixmax,:,:) &
                                               + f%rho_circ(f%ixmax+f%nxguard:f%ixmax:-1,:,:)
    end if

    if (zlbnd==neumann) then
     f%rho_circ(:,f%izmin:f%izmin+f%nzguard,:) = f%rho_circ(:,f%izmin:f%izmin+f%nzguard,:) &
                                               + f%rho_circ(:,f%izmin:f%izmin-f%nzguard:-1,:)
    end if
    if (zrbnd==neumann) then
     f%rho_circ(:,f%izmax-f%nzguard:f%izmax,:) = f%rho_circ(:,f%izmax-f%nzguard:f%izmax,:) &
                                               + f%rho_circ(:,f%izmax+f%nzguard:f%izmax:-1,:)
    end if
  end if

  ! In rz geometry, divide the current by the cell volume
  if (f%l_2drz) then

     ! In the lower guard cells in x (exchanged with nearby processors)
     do j=f%ixmin-f%nxguard,f%ixmin-1
        r = abs(f%xmin+j*f%dx)
        f%rho(j,:,:) = f%rho(j,:,:)/(2.*pi*r)
        if (f%circ_m>0) f%rho_circ(j,:,:) = f%rho_circ(j,:,:)/(2.*pi*r)
     end do

     ! On the lower boundary
     j = f%ixmin
     if (f%xmin==0.) then
        ! On axis
        if (type_rz_depose == 1) then ! Verboncoeur JCP 164, 421-427 (2001) : corrected volumes
           f%rho(j,:,:) = f%rho(j,:,:)/(pi*f%dx/3.) 
        else                          ! Standard volume
           f%rho(j,:,:) = f%rho(j,:,:)/(pi*f%dx/4.)
        endif
        if (f%circ_m>0) f%rho_circ(j,:,:) = 0.   ! Modes with m>0 have 0 density on axis.
     else
        ! Not the axis
        r = abs(f%xmin+j*f%dx)
        f%rho(j,:,:) = f%rho(j,:,:)/(2.*pi*r)
        if (f%circ_m>0) f%rho_circ(j,:,:) = f%rho_circ(j,:,:)/(2.*pi*r)
     end if
     
     ! In the rest of the grid
     do j=f%ixmin+1,f%ixmax+f%nxguard
        r = abs(f%xmin+j*f%dx)
        f%rho(j,:,:) = f%rho(j,:,:)/(2.*pi*r)
        if (f%circ_m>0) f%rho_circ(j,:,:) = f%rho_circ(j,:,:)/(2.*pi*r)
     end do
     
  end if
  
  return
end subroutine em3d_applybc_rho

subroutine em3d_applybc_splite(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==dirichlet) then
     f%exx(f%ixmin-1,:,:) = f%exx(f%ixmin,:,:)
     f%eyx(f%ixmin,:,:) = 0.
     f%ezx(f%ixmin,:,:) = 0.
     f%exy(f%ixmin-1,:,:) = f%exy(f%ixmin,:,:)
     f%eyy(f%ixmin,:,:) = 0.
     f%ezy(f%ixmin,:,:) = 0.
     f%exz(f%ixmin-1,:,:) = f%exz(f%ixmin,:,:)
     f%eyz(f%ixmin,:,:) = 0.
     f%ezz(f%ixmin,:,:) = 0.
  end if
  if (xrbnd==dirichlet) then
     f%exx(f%ixmax,:,:) = f%exx(f%ixmax-1,:,:)
     f%eyx(f%ixmax,:,:) = 0.
     f%ezx(f%ixmax,:,:) = 0.
     f%exy(f%ixmax,:,:) = f%exy(f%ixmax-1,:,:)
     f%eyy(f%ixmax,:,:) = 0.
     f%ezy(f%ixmax,:,:) = 0.
     f%exz(f%ixmax,:,:) = f%exz(f%ixmax-1,:,:)
     f%eyz(f%ixmax,:,:) = 0.
     f%ezz(f%ixmax,:,:) = 0.
  end if

  if (ylbnd==dirichlet) then
     f%exx(:,f%iymin,:) = 0.
     f%eyx(:,f%iymin-1,:) = f%eyx(:,f%iymin,:)
     f%ezx(:,f%iymin,:) = 0.
     f%exy(:,f%iymin,:) = 0.
     f%eyy(:,f%iymin-1,:) = f%eyy(:,f%iymin,:)
     f%ezy(:,f%iymin,:) = 0.
     f%exz(:,f%iymin,:) = 0.
     f%eyz(:,f%iymin-1,:) = f%eyz(:,f%iymin,:)
     f%ezz(:,f%iymin,:) = 0.
  end if
  if (yrbnd==dirichlet) then
     f%exx(:,f%iymax,:) = 0.
     f%eyx(:,f%iymax,:) = f%eyx(:,f%iymax-1,:)
     f%ezx(:,f%iymax,:) = 0.
     f%exy(:,f%iymax,:) = 0.
     f%eyy(:,f%iymax,:) = f%eyy(:,f%iymax-1,:)
     f%ezy(:,f%iymax,:) = 0.
     f%exz(:,f%iymax,:) = 0.
     f%eyz(:,f%iymax,:) = f%eyz(:,f%iymax-1,:)
     f%ezz(:,f%iymax,:) = 0.
  end if

  if (zlbnd==dirichlet) then
     f%exx(:,:,f%izmin) = 0.
     f%eyx(:,:,f%izmin) = 0.
     f%ezx(:,:,f%izmin-1) = f%ezx(:,:,f%izmin)
     f%exy(:,:,f%izmin) = 0.
     f%eyy(:,:,f%izmin) = 0.
     f%ezy(:,:,f%izmin-1) = f%ezy(:,:,f%izmin)
     f%exz(:,:,f%izmin) = 0.
     f%eyz(:,:,f%izmin) = 0.
     f%ezz(:,:,f%izmin-1) = f%ezz(:,:,f%izmin)
  end if
  if (zrbnd==dirichlet) then
     f%exx(:,:,f%izmax) = 0.
     f%eyx(:,:,f%izmax) = 0.
     f%ezx(:,:,f%izmax) = f%ezx(:,:,f%izmax-1)
     f%exy(:,:,f%izmax) = 0.
     f%eyy(:,:,f%izmax) = 0.
     f%ezy(:,:,f%izmax) = f%ezx(:,:,f%izmax-1)
     f%exz(:,:,f%izmax) = 0.
     f%eyz(:,:,f%izmax) = 0.
     f%ezz(:,:,f%izmax) = f%ezx(:,:,f%izmax-1)
  end if

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

  if (zlbnd==neumann) then
     f%ezx(:,:,f%izmin-1) = -f%ezx(:,:,f%izmin  )
     f%ezy(:,:,f%izmin-1) = -f%ezy(:,:,f%izmin  )
     f%ezz(:,:,f%izmin-1) = -f%ezz(:,:,f%izmin  )
  end if
  if (zrbnd==neumann) then
     f%ezx(:,:,f%izmax  ) = -f%ezx(:,:,f%izmax-1)
     f%ezy(:,:,f%izmax  ) = -f%ezy(:,:,f%izmax-1)
     f%ezz(:,:,f%izmax  ) = -f%ezz(:,:,f%izmax-1)
  end if

  return
end subroutine em3d_applybc_splite

subroutine em3d_applybc_splitb(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==dirichlet) then
     f%bxy(f%ixmin,:,:) = 0.
     f%bxz(f%ixmin,:,:) = 0.
     f%byx(f%ixmin-1,:,:) = f%byx(f%ixmin,:,:)
     f%byz(f%ixmin-1,:,:) = f%byz(f%ixmin,:,:)
     f%bzx(f%ixmin-1,:,:) = f%bzx(f%ixmin,:,:)
     f%bzy(f%ixmin-1,:,:) = f%bzy(f%ixmin,:,:)
  end if
  if (xrbnd==dirichlet) then
     f%bxy(f%ixmax,:,:) = 0.
     f%bxz(f%ixmax,:,:) = 0.
     f%byx(f%ixmax,:,:) = f%byx(f%ixmax-1,:,:)
     f%byz(f%ixmax,:,:) = f%byz(f%ixmax-1,:,:)
     f%bzx(f%ixmax,:,:) = f%bzx(f%ixmax-1,:,:)
     f%bzy(f%ixmax,:,:) = f%bzy(f%ixmax-1,:,:)
  end if

  if (ylbnd==dirichlet) then
     f%bxy(:,f%iymin-1,:) = f%bxy(:,f%iymin,:)
     f%bxz(:,f%iymin-1,:) = f%bxz(:,f%iymin,:)
     f%byx(:,f%iymin,:) = 0.
     f%byz(:,f%iymin,:) = 0.
     f%bzx(:,f%iymin-1,:) = f%bzx(:,f%iymin,:)
     f%bzy(:,f%iymin-1,:) = f%bzy(:,f%iymin,:)
  end if
  if (yrbnd==dirichlet) then
     f%bxy(:,f%iymax,:) = f%bxy(:,f%iymax-1,:)
     f%bxz(:,f%iymax,:) = f%bxz(:,f%iymax-1,:)
     f%byx(:,f%iymax,:) = 0.
     f%byz(:,f%iymax,:) = 0.
     f%bzx(:,f%iymax,:) = f%bzx(:,f%iymax-1,:)
     f%bzy(:,f%iymax,:) = f%bzy(:,f%iymax-1,:)
  end if

  if (zlbnd==dirichlet) then
     f%bxy(:,:,f%izmin-1) = f%bxy(:,:,f%izmin)
     f%bxz(:,:,f%izmin-1) = f%bxz(:,:,f%izmin)
     f%byx(:,:,f%izmin-1) = f%byx(:,:,f%izmin)
     f%byz(:,:,f%izmin-1) = f%byz(:,:,f%izmin)
     f%bzx(:,:,f%izmin) = 0.
     f%bzy(:,:,f%izmin) = 0.
  end if
  if (zrbnd==dirichlet) then
     f%bxy(:,:,f%izmax) = f%bxy(:,:,f%izmax-1)
     f%bxz(:,:,f%izmax) = f%bxz(:,:,f%izmax-1)
     f%byx(:,:,f%izmax) = f%byx(:,:,f%izmax-1)
     f%byz(:,:,f%izmax) = f%byz(:,:,f%izmax-1)
     f%bzx(:,:,f%izmax) = 0.
     f%bzy(:,:,f%izmax) = 0.
  end if

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

  if (zlbnd==neumann) then
     f%bxy(:,:,f%izmin-1) = -f%bxy(:,:,f%izmin  )
     f%bxz(:,:,f%izmin-1) = -f%bxz(:,:,f%izmin  )
     f%byx(:,:,f%izmin-1) = -f%byx(:,:,f%izmin  )
     f%byz(:,:,f%izmin-1) = -f%byz(:,:,f%izmin  )
  end if
  if (zrbnd==neumann) then
     f%bxy(:,:,f%izmax  ) = -f%bxy(:,:,f%izmax-1)
     f%bxz(:,:,f%izmax  ) = -f%bxz(:,:,f%izmax-1)
     f%byx(:,:,f%izmax  ) = -f%byx(:,:,f%izmax-1)
     f%byz(:,:,f%izmax  ) = -f%byz(:,:,f%izmax-1)
  end if

  return
end subroutine em3d_applybc_splitb

subroutine em3d_exchange_bnde_x(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Send/exchange the electric field at the x boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  integer(ISZ)   ::ibuf
#ifdef MPIPARALLEL
  integer(ISZ)   ::ix, n_slices, bufsize
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
        if(l_mpiverbose) write(STDOUT,*) '-- sending e along x'

        if (fl%proc/=my_index) then   
           ! --- send data down in x

           ! Number of slices to communicate along x
           ! (yfu%nxguard for ex, yfu%nxguard-1 for ey and ez)
           n_slices = (3*yfu%nxguard-2)
           bufsize = n_slices*size(yfu%ez(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%ez_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the ex, ey and ez slices into that buffer
           ! ex
           do ix = yfu%ixmin,yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%ex(ix,:,:), ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%ex_circ(ix,:,:), ibuf)
           end do
           ! ey
           if (yfu%nxguard>1) then
              do ix = yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
                 call mympi_pack(yfu%ey(ix,:,:),ibuf)
                 if (yfu%circ_m > 0) call mympi_pack(yfu%ey_circ(ix,:,:),ibuf)
              end do
              ! ez
              do ix = yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
                 call mympi_pack(yfu%ez(ix,:,:),ibuf)
                 if (yfu%circ_m > 0) call mympi_pack(yfu%ez_circ(ix,:,:),ibuf)
              end do
           end if
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then
           ! --- send data up in x

           ! Number of slices to communicate along x
           ! (yfl%nxguard for ex, yfl%nxguard-1 for ey and ez)
           n_slices = (3*yfl%nxguard-2)
           bufsize = n_slices*size( yfl%ez(0,:,:))
           ! Check whether to also pack the circ arrays 
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%ez_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the ex, ey and ez slices into that buffer
           ! ex
           do ix =yfl%ixmax-yfl%nxguard,yfl%ixmax-1
              call mympi_pack(yfl%ex(ix,:,:),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%ex_circ(ix,:,:), ibuf)
           end do
           ! ey
           if (yfl%nxguard>1) then
              do ix = yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
                 call mympi_pack(yfl%ey(ix,:,:),ibuf)
                 if (yfl%circ_m > 0) call mympi_pack(yfl%ey_circ(ix,:,:), ibuf)
              end do
              ! ez
              do ix = yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
                 call mympi_pack(yfl%ez(ix,:,:),ibuf)
                 if (yfl%circ_m > 0) call mympi_pack(yfl%ez_circ(ix,:,:), ibuf)
              end do
           end if
           call mpi_isend_pack(fu%proc,2,ibuf)
        else
#endif
           ! Arrays are on the same processor, no need to send a buffer through mpi
           ! Instead exchange the data directly from array to array
           yfl%ex(yfl%ixmax  :yfl%ixmaxg-1,:,:) = &
                yfu%ex(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)
           yfl%ey(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = &
                yfu%ey(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
           yfl%ez(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = &
                yfu%ez(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)

           yfu%ex(yfu%ixming  :yfu%ixmin-1,:,:) = &
                yfl%ex(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)
           yfu%ey(yfu%ixming+1:yfu%ixmin-1,:,:) = &
                yfl%ey(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
           yfu%ez(yfu%ixming+1:yfu%ixmin-1,:,:) = &
                yfl%ez(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

           if (yfu%circ_m > 0) then
              yfl%ex_circ(yfl%ixmax  :yfl%ixmaxg-1,:,:) = &
                   yfu%ex_circ(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)
              yfl%ey_circ(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = &
                   yfu%ey_circ(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
              yfl%ez_circ(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = &
                   yfu%ez_circ(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)

              yfu%ex_circ(yfu%ixming  :yfu%ixmin-1,:,:) = &
                   yfl%ex_circ(yfl%ixmax-yfl%nxguard  :yfl%ixmax-1,:,:)
              yfu%ey_circ(yfu%ixming+1:yfu%ixmin-1,:,:) = &
                   yfl%ey_circ(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
              yfu%ez_circ(yfu%ixming+1:yfu%ixmin-1,:,:) = &
                   yfl%ez_circ(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
           endif

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
subroutine em3d_exchange_bnde_xrecv(fl,fu,ibuf)
  ! -------------------------------------------------
  ! Receive the electric field at the x boundary
  ! -------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ix, ibuf, bufsize, n_slices

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf
        
        if(l_mpiverbose) write(STDOUT,*) '-- receiving e along x'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in x

           ! Number of slices to communicate along x
           ! (yfu%nxguard for ex, yfu%nxguard-1 for ey and ez)
           n_slices = (3*yfu%nxguard-2)
           bufsize = n_slices*size(yfu%ez(yfu%ixmin,:,:))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%ez_circ(yfu%ixmin,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the ex, ey and ez slices into that buffer
           ! ex
           call mpi_recv_pack(fl%proc,2,ibuf)
           do ix = yfu%ixming,yfu%ixmin-1
              yfu%ex(ix,:,:) = reshape( mpi_unpack_real_array( size(yfu%Ex(0,:,:)), ibuf), shape(yfu%Ex(0,:,:)))
              if (yfu%circ_m > 0) &
                   yfu%ex_circ(ix,:,:) = &
                   reshape( mpi_unpack_complex_array( size(yfu%Ex_circ(0,:,:)), ibuf), shape(yfu%Ex_circ(0,:,:)))
           end do
           if (yfu%nxguard>1) then
              ! ey
              do ix = yfu%ixming+1,yfu%ixmin-1
                 yfu%ey(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%Ey(0,:,:)),ibuf),shape(yfu%Ey(0,:,:)))
                 if (yfu%circ_m > 0) &
                      yfu%ey_circ(ix,:,:) = &
                      reshape( mpi_unpack_complex_array( size(yfu%Ey_circ(0,:,:)), ibuf), shape(yfu%Ey_circ(0,:,:)))
              end do
              ! ez
              do ix = yfu%ixming+1,yfu%ixmin-1
                 yfu%ez(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%Ez(0,:,:)),ibuf),shape(yfu%Ez(0,:,:)))
                 if (yfu%circ_m > 0) &
                      yfu%ez_circ(ix,:,:) = &
                      reshape( mpi_unpack_complex_array( size(yfu%Ez_circ(0,:,:)), ibuf), shape(yfu%Ez_circ(0,:,:)))
              end do
           end if
           
        else if (fu%proc/=my_index) then
           ! --- recv data from up in x

           ! Number of slices to communicate along x
           ! (yfl%nxguard for ex, yfl%nxguard-1 for ey and ez)
           n_slices = (3*yfl%nxguard-2)
           bufsize = n_slices*size(yfl%ez(yfl%ixmin,:,:))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%ez_circ(yfl%ixmin,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the ex, ey and ez slices into that buffer
           ! ex
           call mpi_recv_pack(fu%proc,1,ibuf)
           do ix = yfl%ixmax,yfl%ixmaxg-1
              yfl%ex(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%Ex(0,:,:)),ibuf),shape(yfl%Ex(0,:,:)))
              if (yfl%circ_m > 0) &
                   yfl%ex_circ(ix,:,:) = &
                   reshape( mpi_unpack_complex_array( size(yfl%Ex_circ(0,:,:)), ibuf), shape(yfl%Ex_circ(0,:,:)))
           end do
           if (yfl%nxguard>1) then
              do ix = yfl%ixmax+1,yfl%ixmaxg-1
                 yfl%ey(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%Ey(0,:,:)),ibuf),shape(yfl%Ey(0,:,:)))
                 if (yfl%circ_m > 0) &
                      yfl%ey_circ(ix,:,:) = &
                   reshape( mpi_unpack_complex_array( size(yfl%Ey_circ(0,:,:)), ibuf), shape(yfl%Ey_circ(0,:,:)))
              end do
              do ix = yfl%ixmax+1,yfl%ixmaxg-1
                 yfl%ez(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%Ez(0,:,:)),ibuf),shape(yfl%Ez(0,:,:)))
                 if (yfl%circ_m > 0) &
                      yfl%ez_circ(ix,:,:) = &
                   reshape( mpi_unpack_complex_array( size(yfl%Ez_circ(0,:,:)), ibuf), shape(yfl%Ez_circ(0,:,:)))
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

subroutine em3d_exchange_bndb_x(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Send/exchange the magnetic field at the x boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ibuf, n_slices, bufsize
#ifdef MPIPARALLEL
  integer(MPIISZ)::mpirequest(2),mpierror
  integer(ISZ) :: ix
  if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf
#ifdef MPIPARALLEL
        if(l_mpiverbose) write(STDOUT,*) '-- sending b along x'

        if (fl%proc/=my_index) then
           ! --- send data down in x

           ! Number of slices to communicate along x
           ! (yfu%nxguard for by and bz, yfu%nxguard-1 for bx)
           n_slices = (3*yfu%nxguard-1)
           bufsize = n_slices*size(yfu%by(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%by_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the bx, by and bz slices into that buffer
           ! bx
           do ix=yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%bx(ix,:,:),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%bx_circ(ix,:,:),ibuf)
           end do
           ! by
           do ix=yfu%ixmin, yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%by(ix,:,:),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%by_circ(ix,:,:),ibuf)
           end do
           ! bz
           do ix=yfu%ixmin, yfu%ixmin+yfu%nxguard-1
              call mympi_pack(yfu%bz(ix,:,:),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%bz_circ(ix,:,:),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)
           
        else if (fu%proc/=my_index) then
           ! --- send data up in x

           ! Number of slices to communicate along z
           ! (yfl%nxguard for by and bz, yfl%nxguard-1 for bx)
           n_slices = (3*yfl%nxguard-1)
           bufsize = n_slices*size(yfl%by(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%by_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the bx, by and bz slices into that buffer
           ! bx
           do ix=yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
              call mympi_pack(yfl%bx(ix,:,:),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%bx_circ(ix,:,:),ibuf)
           end do
           ! by
           do ix=yfl%ixmax-yfl%nxguard,yfl%ixmax-1
              call mympi_pack(yfl%by(ix,:,:),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%by_circ(ix,:,:),ibuf)
           end do
           ! bz
           do ix=yfl%ixmax-yfl%nxguard,yfl%ixmax-1
              call mympi_pack(yfl%bz(ix,:,:),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%bz_circ(ix,:,:),ibuf)
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
          syfu%bxy(syfu%ixming+1:syfu%ixmin-1,:,:) = 0.
          syfu%bxz(syfu%ixming+1:syfu%ixmin-1,:,:) = yfl%bx(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)*syfu%clight
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
          syfl%bxy(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = 0.
          syfl%bxz(syfl%ixmax+1:syfl%ixmaxg-1,:,:) = yfu%bx(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)*syfl%clight
          syfl%byx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = yfu%by(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)*syfl%clight
          syfl%byz(syfl%ixmax  :syfl%ixmaxg-1,:,:) = 0.
          syfl%bzx(syfl%ixmax  :syfl%ixmaxg-1,:,:) = yfu%bz(yfu%ixmin  :yfu%ixmin+yfu%nxguard-1,:,:)*syfl%clight
          syfl%bzy(syfl%ixmax  :syfl%ixmaxg-1,:,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
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
subroutine em3d_exchange_bndb_xrecv(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Receive the magnetic field at the x boundary
  ! ---------------------------------------------------  
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ibuf,ix, n_slices, bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving b along x'

        if (fl%proc/=my_index) then
           ! --- recv data from down in x

           ! Number of slices to communicate along x
           ! (yfu%nxguard for by and bz, yfu%nxguard-1 for bx)
           n_slices = (3*yfu%nxguard-1)
           bufsize = n_slices*size(yfu%bx(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%bx_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the bx, by and bz slices from that buffer
           call mpi_recv_pack(fl%proc,2,ibuf)
           do ix=yfu%ixming+1,yfu%ixmin-1
              yfu%bx(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%bx(ix,:,:)),ibuf),shape(yfu%bx(ix,:,:)))
              if (yfu%circ_m > 0) &
                   yfu%bx_circ(ix,:,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%bx_circ(ix,:,:)),ibuf),shape(yfu%bx_circ(ix,:,:)))
           end do
           do ix=yfu%ixming,yfu%ixmin-1
              yfu%by(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%by(ix,:,:)),ibuf),shape(yfu%by(ix,:,:)))
              if (yfu%circ_m > 0) &
                   yfu%by_circ(ix,:,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%by_circ(ix,:,:)),ibuf),shape(yfu%by_circ(ix,:,:)))
           end do
           do ix=yfu%ixming,yfu%ixmin-1
              yfu%bz(ix,:,:) = reshape(mpi_unpack_real_array( size(yfu%bz(ix,:,:)),ibuf),shape(yfu%bz(ix,:,:)))
              if (yfu%circ_m > 0) &
                   yfu%bz_circ(ix,:,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%bz_circ(ix,:,:)),ibuf),shape(yfu%bz_circ(ix,:,:)))
           end do
           
        else if (fu%proc/=my_index) then
           ! --- recv data from up in z

           ! Number of slices to communicate along x
           ! (yfl%nxguard for by and bz, yfl%nxguard-1 for bx)
           n_slices = (3*yfl%nxguard-1)
           bufsize = n_slices*size(yfl%bx(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%bx_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the bx, by and bz slices from that buffer
           call mpi_recv_pack(fu%proc,1,ibuf)
           do ix=yfl%ixmax+1,yfl%ixmaxg-1
              yfl%bx(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%bx(ix,:,:)),ibuf),shape(yfl%bx(ix,:,:)))
              if (yfl%circ_m > 0) &
                   yfl%bx_circ(ix,:,:) = &
                   reshape(mpi_unpack_complex_array( size(yfl%bx_circ(ix,:,:)),ibuf),shape(yfl%bx_circ(ix,:,:)))
           end do
           do ix=yfl%ixmax,yfl%ixmaxg-1
              yfl%by(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%by(ix,:,:)),ibuf),shape(yfl%by(ix,:,:)))
              if (yfl%circ_m > 0) &
                   yfl%by_circ(ix,:,:) = &
                   reshape(mpi_unpack_complex_array( size(yfl%by_circ(ix,:,:)),ibuf),shape(yfl%by_circ(ix,:,:)))
           end do
           do ix=yfl%ixmax,yfl%ixmaxg-1
              yfl%bz(ix,:,:) = reshape(mpi_unpack_real_array( size(yfl%bz(ix,:,:)),ibuf),shape(yfl%bz(ix,:,:)))
              if (yfl%circ_m > 0) &
                   yfl%bz_circ(ix,:,:) = &
                   reshape(mpi_unpack_complex_array( size(yfl%bz_circ(ix,:,:)),ibuf),shape(yfl%bz_circ(ix,:,:)))
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

subroutine em3d_exchange_bndf_x(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Send/exchange the fields f at the x boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  integer(ISZ)   ::ibuf, n_slices, bufsize
#ifdef MPIPARALLEL
  integer(ISZ)   ::ix
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
        if(l_mpiverbose) write(STDOUT,*) '-- sending f along x'
        
        if (fl%proc/=my_index) then
           ! --- send data down in x

           if (yfu%nxguard>1) then
              ! Number of slices to communicate along x
              n_slices = (yfu%nxguard-1)
              bufsize = n_slices*size(yfu%f(0,:,:))
              ! Check whether to also pack the circ arrays
              if (yfu%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfu%f_circ(0,:,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Pack f
              do ix = yfu%ixmin+1,yfu%ixmin+yfu%nxguard-1
                 call mympi_pack(yfu%f(ix,:,:),ibuf)
                 if (yfu%circ_m > 0) call mympi_pack(yfu%f_circ(ix,:,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,1,ibuf)
           endif
           
        else if (fu%proc/=my_index) then
           ! --- send data up in x
           
           if (yfl%nxguard>1) then
              ! Number of slices to communicate along x
              n_slices = (yfl%nxguard-1)
              bufsize = n_slices*size(yfl%f(0,:,:))
              ! Check whether to also pack the circ arrays
              if (yfl%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfl%f_circ(0,:,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Pack f
              do ix = yfl%ixmax-yfl%nxguard+1,yfl%ixmax-1
                 call mympi_pack(yfl%f(ix,:,:),ibuf)
                 if (yfl%circ_m > 0) call mympi_pack(yfl%f_circ(ix,:,:),ibuf)
              end do
              call mpi_isend_pack(fu%proc,2,ibuf)
           end if
        else
#endif
           ! Arrays are on the same processor, no need to send a buffer through mpi
           ! Instead exchange the data directly from array to array
           yfl%f(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%f(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
           yfu%f(yfu%ixming+1:yfu%ixmin-1,:,:) = yfl%f(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

           if (yfu%circ_m > 0) then
              yfl%f_circ(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = &
                   yfu%f_circ(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
              yfu%f_circ(yfu%ixming+1:yfu%ixmin-1,:,:) = &
                   yfl%f_circ(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
           endif

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
subroutine em3d_exchange_bndf_xrecv(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ix,ibuf,n_slices,bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving f along x'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in x
           
           if (yfu%nxguard>1) then
              ! Number of slices to communicate along x
              n_slices = (yfu%nxguard-1)
              bufsize = n_slices*size(yfu%f(yfu%ixmin,:,:))
              ! Check whether to also pack the circ arrays
              if (yfu%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfu%f_circ(yfu%ixmin,:,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Receive f
              call mpi_recv_pack(fl%proc,2,ibuf)
              do ix = yfu%ixming+1,yfu%ixmin-1
                 yfu%f(ix,:,:) = reshape(mpi_unpack_real_array( &
                      size(yfu%F(0,:,:)),ibuf), shape(yfu%F(0,:,:)))
                 if (yfu%circ_m>0) &
                      yfu%f_circ(ix,:,:) = reshape(mpi_unpack_complex_array( &
                      size(yfu%F_circ(0,:,:)),ibuf), shape(yfu%F_circ(0,:,:)))
              end do
           end if
           
        else if (fu%proc/=my_index) then
           ! --- recv data from up in x

           if (yfl%nxguard>1) then

              ! Number of slices to communicate along x
              n_slices = (yfl%nxguard-1)
              bufsize = n_slices*size(yfl%f(yfl%ixmin,:,:))
              ! Check whether to also pack the circ arrays
              if (yfl%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfl%f_circ(yfl%ixmin,:,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Receive f
              call mpi_recv_pack(fu%proc,1,ibuf)
              do ix = yfl%ixmax+1,yfl%ixmaxg-1
                 yfl%f(ix,:,:) = reshape(mpi_unpack_real_array( &
                      size(yfl%F(0,:,:)),ibuf),shape(yfl%F(0,:,:)))
                 if (yfl%circ_m>0) &
                      yfl%f_circ(ix,:,:) = reshape(mpi_unpack_complex_array( &
                      size(yfl%F_circ(0,:,:)),ibuf), shape(yfl%F_circ(0,:,:)))
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

subroutine em3d_exchange_bndj_x(fl,fu,ibuf)
  ! --------------------------------------------
  ! Send/exchange the current at the x boundary
  ! --------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ix,ibuf,nguardinl,nguardinu, n_slices, bufsize

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

        if(l_mpiverbose) write(STDOUT,*) '-- sending j along x'
        
        if (fl%proc/=my_index) then
           ! --- send data down in x

           ! Number of slices to communicate along x
           nguardinu = yfu%nxguard
           n_slices = (3*(yfu%nzguard+nguardinu) + 2)
           bufsize = n_slices * size(yfu%J(0,:,:,1))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%J_circ(0,:,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the Jx, Jy and Jz slices into that buffer
           ! Jx
           do ix = -yfu%nxguard,-1+nguardinu
              call mympi_pack(yfu%J(ix,:,:,1),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%J_circ(ix,:,1,:),ibuf)
           end do
           ! Jy and Jz
           do ix = -yfu%nxguard,nguardinu
              call mympi_pack(yfu%J(ix,:,:,2),ibuf)
              call mympi_pack(yfu%J(ix,:,:,3),ibuf)
              if (yfu%circ_m > 0) then
                 call mympi_pack( yfu%J_circ(ix,:,2,:), ibuf)
                 call mympi_pack( yfu%J_circ(ix,:,3,:), ibuf)
              endif
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then
           ! --- send data up in x

           ! Number of slices to communicate along x
           nguardinl = yfl%nxguard
           n_slices = (3*(yfl%nxguard+nguardinl) + 2)
           bufsize = n_slices * size(yfl%J(0,:,:,1))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%J_circ(0,:,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the Jx, Jy and Jz slices into that buffer
           ! Jx
           do ix = yfl%nx-nguardinl, yfl%nx+yfl%nxguard-1
              call mympi_pack(yfl%J(ix,:,:,1),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%J_circ(ix,:,1,:),ibuf)
           end do
           ! Jy and Jz
           do ix = yfl%nx-nguardinl, yfl%nx+yfl%nxguard
              call mympi_pack(yfl%J(ix,:,:,2),ibuf)
              call mympi_pack(yfl%J(ix,:,:,3),ibuf)
              if (yfl%circ_m > 0) then
                 call mympi_pack( yfl%J_circ(ix,:,2,:),ibuf )
                 call mympi_pack( yfl%J_circ(ix,:,3,:),ibuf )
              endif
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)

        else
#endif
           ! Arrays are on the same processor, no need to send a buffer through mpi
           ! Instead exchange the data directly from array to array
           nguardinu = yfu%nxguard
           nguardinl = yfl%nxguard

           yfu%J(-nguardinu:yfu%nxguard-1,:,:,1)  = yfu%J(-nguardinu:yfu%nxguard-1,:,:,1)  &
                + yfl%J(yfl%nx-nguardinl:yfl%nx+yfl%nxguard-1,:,:,1) 
           yfu%J(-nguardinu:yfu%nxguard,:,:,2:3)  = yfu%J(-nguardinu:yfu%nxguard,:,:,2:3)  &
                + yfl%J(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,:,2:3)
           yfl%J(yfl%nx-nguardinl:yfl%nx+yfl%nxguard-1,:,:,1) = yfu%J(-nguardinu:yfu%nxguard-1,:,:,1)
           yfl%J(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,:,2:3) = yfu%J(-nguardinu:yfu%nxguard,:,:,2:3)

           if (yfu%circ_m > 0) then
              yfu%J_circ(-nguardinu:yfu%nxguard-1,:,1,:) = yfu%J_circ(-nguardinu:yfu%nxguard-1,:,1,:)  &
                   + yfl%J_circ(yfl%nx-nguardinl:yfl%nx+yfl%nxguard-1,:,1,:) 
              yfu%J_circ(-nguardinu:yfu%nxguard,:,2:3,:) = yfu%J_circ(-nguardinu:yfu%nxguard,:,2:3,:)  &
                   + yfl%J_circ(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,2:3,:)
              yfl%J_circ(yfl%nx-nguardinl:yfl%nx+yfl%nxguard-1,:,1,:) = yfu%J_circ(-nguardinu:yfu%nxguard-1,:,1,:)
              yfl%J_circ(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,2:3,:) = yfu%J_circ(-nguardinu:yfu%nxguard,:,2:3,:)
           endif
           
#ifdef MPIPARALLEL
        end if
#endif
     end select
  end select
  return
end subroutine em3d_exchange_bndj_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_xrecv(fl,fu,ibuf)
  ! -----------------------------------------
  ! Receive the current at the x boundary
  ! -----------------------------------------  
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ix,ibuf,nguardinl,nguardinu,n_slices,bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving j along x'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in x

           ! Number of slices to communicate along x
           nguardinu = yfu%nxguard
           n_slices = (3*(yfu%nxguard+nguardinu) + 2)
           bufsize = n_slices*size(yfu%J(0,:,:,1))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%J_circ(0,:,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the Jx, Jy and Jz slices from that buffer
           call mpi_recv_pack(fl%proc,2,ibuf)
           do ix = -nguardinu,yfu%nxguard-1
              yfu%J(ix,:,:,1) = yfu%J(ix,:,:,1) + reshape(mpi_unpack_real_array( size(yfu%J(0,:,:,1)),ibuf), &
                   shape(yfu%J(0,:,:,1)))
              if ( yfu%circ_m > 0 ) &
                   yfu%J_circ(ix,:,1,:) = yfu%J_circ(ix,:,1,:) + reshape( mpi_unpack_complex_array( &
                   size(yfu%J_circ(0,:,1,:)),ibuf), shape(yfu%J_circ(0,:,1,:)) )
           end do
           do ix = -nguardinu,yfu%nxguard
              yfu%J(ix,:,:,2) = yfu%J(ix,:,:,2) + reshape(mpi_unpack_real_array( size(yfu%J(0,:,:,1)),ibuf), &
                                                                                shape(yfu%J(0,:,:,1)))
              yfu%J(ix,:,:,3) = yfu%J(ix,:,:,3) + reshape(mpi_unpack_real_array( size(yfu%J(0,:,:,1)),ibuf), &
                   shape(yfu%J(0,:,:,1)))
              if ( yfu%circ_m > 0 ) then
                 yfu%J_circ(ix,:,2,:) = yfu%J_circ(ix,:,2,:) + reshape( mpi_unpack_complex_array( &
                      size(yfu%J_circ(0,:,2,:)),ibuf), shape(yfu%J_circ(0,:,2,:)) )
                 yfu%J_circ(ix,:,3,:) = yfu%J_circ(ix,:,3,:) + reshape( mpi_unpack_complex_array( &
                      size(yfu%J_circ(0,:,3,:)),ibuf), shape(yfu%J_circ(0,:,3,:)) )
              endif
           end do

        else if (fu%proc/=my_index) then
           ! --- recv data from up in x

           ! Number of slices to communicate along x
           nguardinl = yfl%nxguard
           n_slices = (3*(yfl%nxguard+nguardinl) + 2)
           bufsize = n_slices*size(yfl%J(0,:,:,1))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%J_circ(0,:,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the Jx, Jy and Jz slices from that buffer
           call mpi_recv_pack(fu%proc,1,ibuf)
           do ix = -yfl%nxguard,nguardinl-1
              yfl%J(yfl%nx+ix,:,:,1) = yfl%J(yfl%nx+ix,:,:,1) + reshape( &
                   mpi_unpack_real_array( size(yfl%J(yfl%nx-1,:,:,1)),ibuf), shape(yfl%J(yfl%nx-1,:,:,1)))
              if ( yfl%circ_m > 0 ) &
                   yfl%J_circ(yfl%nx+ix,:,1,:) = yfl%J_circ(yfl%nx+ix,:,1,:) + reshape( &
                   mpi_unpack_complex_array( size(yfl%J_circ(yfl%nx-1,:,1,:)),ibuf), shape(yfl%J_circ(yfl%nx-1,:,1,:)))
           end do
           do ix = -yfl%nxguard,nguardinl
              yfl%J(yfl%nx+ix,:,:,2) = yfl%J(yfl%nx+ix,:,:,2) + reshape( &
                   mpi_unpack_real_array( size(yfl%J(yfl%nx-1,:,:,2)),ibuf), shape(yfl%J(yfl%nx-1,:,:,2)))
              yfl%J(yfl%nx+ix,:,:,3) = yfl%J(yfl%nx+ix,:,:,3) + reshape( &
                   mpi_unpack_real_array( size(yfl%J(yfl%nx-1,:,:,3)),ibuf), shape(yfl%J(yfl%nx-1,:,:,3)))
              if ( yfl%circ_m > 0 ) then
                 yfl%J_circ(yfl%nx+ix,:,2,:) = yfl%J_circ(yfl%nx+ix,:,2,:) + reshape( &
                      mpi_unpack_complex_array( size(yfl%J_circ(yfl%nx-1,:,2,:)),ibuf), shape(yfl%J_circ(yfl%nx-1,:,2,:)))
                 yfl%J_circ(yfl%nx+ix,:,3,:) = yfl%J_circ(yfl%nx+ix,:,3,:) + reshape( &
                      mpi_unpack_complex_array( size(yfl%J_circ(yfl%nx-1,:,3,:)),ibuf), shape(yfl%J_circ(yfl%nx-1,:,3,:)))
              endif
           end do
        end if
     end select
  end select

  return
end subroutine em3d_exchange_bndj_xrecv
#endif

subroutine em3d_exchange_bndrho_x(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Send/exchange the field rho at the x boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ix,ibuf,nguardinl,nguardinu,n_slices,bufsize

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
        if(l_mpiverbose) write(STDOUT,*) '-- sending rho along x'
        
        if (fl%proc/=my_index) then
           ! --- send data down in x

           ! Number of slices to communicate along z
           nguardinu = yfu%nxguard
           n_slices =  yfu%nxguard + nguardinu + 1
           bufsize = n_slices*size(yfu%rho(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%rho_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Pack rho
           do ix = -yfu%nxguard,nguardinu
              call mympi_pack(yfu%rho(ix,:,:),ibuf)
              if (yfu%circ_m > 0) &
                   call mympi_pack(yfu%rho_circ(ix,:,:),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then
           ! --- send data up in z

           ! Number of slices to communicate along z
           nguardinl = yfl%nxguard
           n_slices =  yfl%nxguard + nguardinl + 1
           bufsize = n_slices*size(yfl%rho(yfl%nx,:,:))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%rho_circ(yfl%nx,:,:))
           ! Allocate a buffer array in mpibuffer           
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Pack rho
           do ix = yfl%nx-nguardinl,yfl%nx+yfl%nxguard
              call mympi_pack(yfl%rho(ix,:,:),ibuf)
              if (yfl%circ_m > 0) &
                   call mympi_pack(yfl%rho_circ(ix,:,:),ibuf)
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)

        else
#endif
           nguardinu = yfu%nxguard
           nguardinl = yfl%nxguard
           yfu%Rho(-nguardinu:yfu%nxguard,:,:) = yfu%Rho(-nguardinu:yfu%nxguard,:,:) &
                + yfl%Rho(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,:)
           yfl%Rho(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,:) = yfu%Rho(-nguardinu:yfu%nxguard,:,:)
           if ( yfu%circ_m>0 ) then
              yfu%Rho_circ(-nguardinu:yfu%nxguard,:,:) = yfu%Rho_circ(-nguardinu:yfu%nxguard,:,:) &
                   + yfl%Rho_circ(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,:)
              yfl%Rho_circ(yfl%nx-nguardinl:yfl%nx+yfl%nxguard,:,:) = &
                   yfu%Rho_circ(-nguardinu:yfu%nxguard,:,:)
           endif
#ifdef MPIPARALLEL
        end if
#endif
     end select
  end select
  return
end subroutine em3d_exchange_bndrho_x

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_xrecv(fl,fu,ibuf)
  ! --------------------------------------------
  ! Receive the field rho at the x boundary
  ! --------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ix, ibuf,nguardinl,nguardinu, n_slices, bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf
        
        if(l_mpiverbose) write(STDOUT,*) '-- receiving rho along x'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in x

           ! Number of slices to communicate along x
           nguardinu = yfu%nxguard
           n_slices =  yfu%nxguard + nguardinu + 1
           bufsize = n_slices*size(yfu%rho(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%rho_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Receive rho
           call mpi_recv_pack(fl%proc,2,ibuf)
           do ix = -nguardinu,yfu%nxguard
              yfu%rho(ix,:,:) = yfu%rho(ix,:,:) + reshape( &
                   mpi_unpack_real_array( size(yfu%rho(0,:,:)),ibuf), shape(yfu%rho(0,:,:)))
              if(yfu%circ_m>0) &
                 yfu%rho_circ(ix,:,:) = yfu%rho_circ(ix,:,:) + reshape( &
                 mpi_unpack_complex_array( size(yfu%rho_circ(0,:,:)),ibuf), &
                 shape(yfu%rho_circ(0,:,:)))
           end do

        else if (fu%proc/=my_index) then
           ! --- recv data from up in x

           ! Number of slices to communicate along x
           nguardinl = yfl%nxguard
           n_slices =  yfl%nxguard + nguardinl + 1
           bufsize = n_slices*size(yfl%rho(0,:,:))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%rho_circ(0,:,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Receive rho
           call mpi_recv_pack(fu%proc,1,ibuf)
           do ix = yfl%nx-yfl%nxguard,yfl%nx+nguardinl
              yfl%rho(ix,:,:) = yfl%rho(ix,:,:) + reshape( &
                   mpi_unpack_real_array(size(yfl%rho(yfl%nx,:,:)),ibuf),&
                   shape(yfl%rho(yfl%nx,:,:)))
              if(yfl%circ_m>0) &
                 yfl%rho_circ(ix,:,:) = yfl%rho_circ(ix,:,:) + reshape( &
                 mpi_unpack_complex_array( size(yfl%rho_circ(yfl%nx,:,:)),ibuf), &
                 shape(yfl%rho_circ(yfl%nx,:,:)))
           end do
        end if
     end select
  end select

  return
end subroutine em3d_exchange_bndrho_xrecv
#endif

subroutine em3d_exchange_bnde_y(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  integer(ISZ)   ::ibuf
#ifdef MPIPARALLEL
  integer(ISZ)   ::iy
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
subroutine em3d_exchange_bnde_yrecv(fl,fu,ibuf)
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

subroutine em3d_exchange_bndb_y(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ)   ::ibuf
#ifdef MPIPARALLEL
  integer(ISZ)   ::iy
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
          syfu%bxy(:,syfu%iyming  :syfu%iymin-1,:) = 0.
          syfu%bxz(:,syfu%iyming  :syfu%iymin-1,:) = yfl%bx(:,yfl%iymax-yfl%nyguard  :yfl%iymax-1,:)*syfu%clight
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
          syfl%bxy(:,syfl%iymax  :syfl%iymaxg-1,:) = 0.
          syfl%bxz(:,syfl%iymax  :syfl%iymaxg-1,:) = yfu%bx(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)*syfl%clight
          syfl%byx(:,syfl%iymax+1:syfl%iymaxg-1,:) = yfu%by(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)*syfl%clight
          syfl%byz(:,syfl%iymax+1:syfl%iymaxg-1,:) = 0.
          syfl%bzx(:,syfl%iymax  :syfl%iymaxg-1,:) = yfu%bz(:,yfu%iymin  :yfu%iymin+yfu%nyguard-1,:)*syfl%clight
          syfl%bzy(:,syfl%iymax  :syfl%iymaxg-1,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
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
subroutine em3d_exchange_bndb_yrecv(fl,fu,ibuf)
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

subroutine em3d_exchange_bndf_y(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  integer(ISZ)   ::ibuf
#ifdef MPIPARALLEL
  integer(ISZ)   ::iy
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
           if (yfu%nyguard>1) then
              call mpi_packbuffer_init((yfu%nyguard-1)*size(yfu%f(:,0,:)),ibuf)
              do iy = yfu%iymin+1,yfu%iymin+yfu%nyguard-1
                 call mympi_pack(yfu%f(:,iy,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,1,ibuf)
           end if
        else if (fu%proc/=my_index) then
           ! --- send data up in z
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
subroutine em3d_exchange_bndf_yrecv(fl,fu,ibuf)
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
              call mpi_packbuffer_init((yfu%nyguard-1)*size(yfu%Ez(:,0,:)),ibuf)
              call mpi_recv_pack(fl%proc,2,ibuf)
              do iy = yfu%iyming+1,yfu%iymin-1
                 yfu%f(:,iy,:) = reshape(mpi_unpack_real_array( size(yfu%F(:,0,:)),ibuf),shape(yfu%F(:,0,:)))
              end do
           end if
        else if (fu%proc/=my_index) then
           ! --- recv data from up in z
           if (yfl%nyguard>1) then
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

subroutine em3d_exchange_bndj_y(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iy,ibuf,nguardinl,nguardinu

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
           nguardinu = yfu%nyguard
           call mpi_packbuffer_init((3*(yfu%nyguard+nguardinu)+2)*size(yfu%J(:,-1,:,1)),ibuf)
           do iy = -yfu%nyguard,nguardinu
              call mympi_pack(yfu%J(:,iy,:,1),ibuf)
           end do
           do iy = -yfu%nyguard,nguardinu-1
              call mympi_pack(yfu%J(:,iy,:,2),ibuf)
           end do
           do iy = -yfu%nyguard,nguardinu
              call mympi_pack(yfu%J(:,iy,:,3),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then

           ! --- send data up in z
           nguardinl = yfl%nyguard
           call mpi_packbuffer_init((3*(yfl%nyguard+nguardinl)+2)*size(yfl%J(:,0,:,1)),ibuf)
           do iy = yfl%ny-nguardinl, yfl%ny+yfl%nyguard
              call mympi_pack(yfl%J(:,iy,:,1),ibuf)
           end do
           do iy = yfl%ny-nguardinl, yfl%ny+yfl%nyguard-1
              call mympi_pack(yfl%J(:,iy,:,2),ibuf)
           end do
           do iy = yfl%ny-nguardinl, yfl%ny+yfl%nyguard
              call mympi_pack(yfl%J(:,iy,:,3),ibuf)
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)

        else
#endif
          nguardinu = yfu%nyguard
          nguardinl = yfl%nyguard
          yfu%J(:,-nguardinu:yfu%nyguard  ,:,1:3:2) = yfu%J(:,-nguardinu:yfu%nyguard,  :,1:3:2) &
                                                    + yfl%J(:,yfl%ny-nguardinl:+yfl%ny+yfl%nyguard  ,:,1:3:2)
          yfu%J(:,-nguardinu:yfu%nyguard-1,:,2    ) = yfu%J(:,-nguardinu:yfu%nyguard-1,:,2    ) &
                                                    + yfl%J(:,yfl%ny-nguardinl:+yfl%ny+yfl%nyguard-1,:,2    ) 

           yfl%J(:,yfl%ny-nguardinl:+yfl%ny+yfl%nyguard  ,:,1:3:2) = yfu%J(:,-nguardinu:yfu%nyguard  ,:,1:3:2)
           yfl%J(:,yfl%ny-nguardinl:+yfl%ny+yfl%nyguard-1,:,2    ) = yfu%J(:,-nguardinu:yfu%nyguard-1,:,2    )
#ifdef MPIPARALLEL
        end if
#endif
     end select
  end select
  return
end subroutine em3d_exchange_bndj_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_yrecv(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iy,ibuf,nguardinl,nguardinu

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf
        if (fl%proc/=my_index) then

           ! --- recv data from down in z
           nguardinu = yfu%nyguard
           call mpi_packbuffer_init((3*(yfu%nyguard+nguardinu)+2)*size(yfu%J(:,0,:,1)),ibuf)
           call mpi_recv_pack(fl%proc,2,ibuf)
           do iy = -nguardinu,yfu%nyguard
              yfu%J(:,iy,:  ,1) = yfu%J(:,iy,:  ,1) + reshape(mpi_unpack_real_array( size(yfu%J(:,0,:,1)),ibuf), &
                                                                                    shape(yfu%J(:,0,:,1)))
            end do
            do iy = -nguardinu,yfu%nyguard-1
              yfu%J(:,iy,:  ,2) = yfu%J(:,iy,:  ,2) + reshape(mpi_unpack_real_array( size(yfu%J(:,0,:,1)),ibuf), &
                                                                                    shape(yfu%J(:,0,:,1)))
            end do
            do iy = -nguardinu,yfu%nyguard
              yfu%J(:,iy,:  ,3) = yfu%J(:,iy,:  ,3) + reshape(mpi_unpack_real_array( size(yfu%J(:,0,:,1)),ibuf), &
                                                                                    shape(yfu%J(:,0,:,1)))
            end do

        else if (fu%proc/=my_index) then

           ! --- recv data from up in z
           nguardinl = yfl%nyguard
           call mpi_packbuffer_init((3*(yfl%nyguard+nguardinl)+2)*size(yfl%J(:,0,:,1)),ibuf)
           call mpi_recv_pack(fu%proc,1,ibuf)
           do iy = -yfl%nyguard,nguardinl
              yfl%J(:,yfl%ny+iy,:,1) = yfl%J(:,yfl%ny+iy,:,1) + reshape(mpi_unpack_real_array( size(yfl%J(:,yfl%ny-1,:,1)),ibuf),&
                                                                                              shape(yfl%J(:,yfl%ny-1,:,1)))
            end do
            do iy = -yfl%nyguard,nguardinl-1
              yfl%J(:,yfl%ny+iy,:,2) = yfl%J(:,yfl%ny+iy,:,2) + reshape(mpi_unpack_real_array( size(yfl%J(:,yfl%ny-1,:,2)),ibuf),&
                                                                                              shape(yfl%J(:,yfl%ny-1,:,2)))
            end do
            do iy = -yfl%nyguard,nguardinl
              yfl%J(:,yfl%ny+iy,:,3) = yfl%J(:,yfl%ny+iy,:,3) + reshape(mpi_unpack_real_array( size(yfl%J(:,yfl%ny-1,:,3)),ibuf),&
                                                                                              shape(yfl%J(:,yfl%ny-1,:,3)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndj_yrecv
#endif

subroutine em3d_exchange_bndrho_y(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iy,ibuf,nguardinl,nguardinu

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
           nguardinu = yfu%nyguard
           ! --- send data down in y
           call mpi_packbuffer_init(size(yfu%rho(:,-yfu%nyguard:nguardinu,:)),ibuf)
           do iy = -yfu%nyguard,nguardinu
              call mympi_pack(yfu%rho(:,iy,:),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then

           ! --- send data up in y
           nguardinl = yfl%nyguard
           call mpi_packbuffer_init(size(yfl%rho(:,-nguardinl:yfl%nyguard,:)),ibuf)
           do iy = -nguardinl,yfl%nyguard
              call mympi_pack(yfl%rho(:,yfl%ny+iy,:),ibuf)
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)

        else
#endif
          nguardinu = yfu%nyguard
          nguardinl = yfl%nyguard
          yfu%Rho(:,-nguardinu:yfu%nyguard,:) = yfu%Rho(:,-nguardinu:yfu%nyguard,:)      &
                                              + yfl%Rho(:,-nguardinl:yfl%nyguard,:)
          yfl%Rho(:,-nguardinl:yfl%nyguard,:) = yfu%Rho(:,-nguardinu:yfu%nyguard,:)                                     
#ifdef MPIPARALLEL
        end if
#endif
     end select
  end select
  return
end subroutine em3d_exchange_bndrho_y

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_yrecv(fl,fu,ibuf)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iy, ibuf,nguardinl,nguardinu

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf
        if (fl%proc/=my_index) then

           ! --- recv data from down in z
           nguardinu = yfu%nyguard
           call mpi_packbuffer_init(size(yfu%rho(:,-nguardinu:yfu%nyguard,:)),ibuf)
           call mpi_recv_pack(fl%proc,2,ibuf)
           do iy = -nguardinu,yfu%nyguard
              yfu%rho(:,iy,:) = yfu%rho(:,iy,:) + reshape(mpi_unpack_real_array( size(yfu%rho(:,0,:)),ibuf), &
                   shape(yfu%rho(:,0,:)))
           end do

        else if (fu%proc/=my_index) then

           ! --- recv data from up in z
           nguardinl = yfl%nyguard
           call mpi_packbuffer_init(size(yfl%rho(:,yfl%ny-yfl%nyguard:yfl%ny+nguardinl,:)),ibuf)
           call mpi_recv_pack(fu%proc,1,ibuf)
           do iy = yfl%ny-yfl%nyguard,yfl%ny+nguardinl
              yfl%rho(:,iy,:) = yfl%rho(:,iy,:) + reshape(mpi_unpack_real_array(size(yfl%rho(:,yfl%ny,:)),ibuf),&
                                                                                             shape(yfl%rho(:,yfl%ny,:)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndrho_yrecv
#endif

subroutine em3d_exchange_bnde_z(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Send/exchange the electric field at the z boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  integer(ISZ)   ::ibuf, n_slices, bufsize
#ifdef MPIPARALLEL
  integer(ISZ)   ::iz
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
        
        if(l_mpiverbose) write(STDOUT,*) '-- sending e along z'
        
        if (fl%proc/=my_index) then
           ! --- send data down in z

           ! Number of slices to communicate along z
           ! (yfu%nzguard for ez, yfu%nzguard-1 for ex and ey)
           n_slices = (3*yfu%nzguard-2)
           bufsize = n_slices*size(yfu%ez(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%ez_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the ez, ex and ey slices into that buffer
           ! ez
           do iz = yfu%izmin,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%ez(:,:,iz),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%ez_circ(:,iz,:),ibuf)
           end do
           if (yfu%nzguard>1) then
              ! ex
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                 call mympi_pack(yfu%ex(:,:,iz),ibuf)
                 if (yfu%circ_m > 0) call mympi_pack(yfu%ex_circ(:,iz,:),ibuf)
              end do
              ! ey
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                 call mympi_pack(yfu%ey(:,:,iz),ibuf)
                 if (yfu%circ_m > 0) call mympi_pack(yfu%ey_circ(:,iz,:),ibuf)
              end do
           end if
           call mpi_isend_pack(fl%proc,1,ibuf)
           
        else if (fu%proc/=my_index) then
           ! --- send data up in z

           ! Number of slices to communicate along z
           ! (yfu%nzguard for ez, yfu%nzguard-1 for ex and ey)
           n_slices = (3*yfl%nzguard-2)
           bufsize = n_slices*size(yfl%ez(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%ez_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the ez, ex and ey slices into that buffer
           ! ez
           do iz =yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%ez(:,:,iz),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%ez_circ(:,iz,:),ibuf)
           end do
           if (yfl%nzguard>1) then
              ! ex
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                 call mympi_pack(yfl%ex(:,:,iz),ibuf)
                 if (yfl%circ_m > 0) call mympi_pack(yfl%ex_circ(:,iz,:),ibuf)
              end do
              ! ey
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                 call mympi_pack(yfl%ey(:,:,iz),ibuf)
                 if (yfl%circ_m > 0) call mympi_pack(yfl%ey_circ(:,iz,:),ibuf)
              end do
           end if
           call mpi_isend_pack(fu%proc,2,ibuf)
        else
#endif
           ! Arrays are on the same processor, no need to send a buffer through mpi
           ! Instead exchange the data directly from array to array
           yfl%ex(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%ex(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
           yfl%ey(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%ey(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
           yfl%ez(:,:,yfl%izmax  :yfl%izmaxg-1) = yfu%ez(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)

           yfu%ex(:,:,yfu%izming+1:yfu%izmin-1) = yfl%ex(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
           yfu%ey(:,:,yfu%izming+1:yfu%izmin-1) = yfl%ey(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
           yfu%ez(:,:,yfu%izming  :yfu%izmin-1) = yfl%ez(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)

           if (yfu%circ_m > 0) then
              yfl%ex_circ(:,yfl%izmax+1:yfl%izmaxg-1,:) = &
                   yfu%ex_circ(:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1,:)
              yfl%ey_circ(:,yfl%izmax+1:yfl%izmaxg-1,:) = &
                   yfu%ey_circ(:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1,:)
              yfl%ez_circ(:,yfl%izmax  :yfl%izmaxg-1,:) = &
                   yfu%ez_circ(:,yfu%izmin  :yfu%izmin+yfu%nzguard-1,:)
              
              yfu%ex_circ(:,yfu%izming+1:yfu%izmin-1,:) = &
                   yfl%ex_circ(:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1,:)
              yfu%ey_circ(:,yfu%izming+1:yfu%izmin-1,:) = &
                   yfl%ey_circ(:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1,:)
              yfu%ez_circ(:,yfu%izming  :yfu%izmin-1,:) = &
                   yfl%ez_circ(:,yfl%izmax-yfl%nzguard  :yfl%izmax-1,:)
          endif
           
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
subroutine em3d_exchange_bnde_zrecv(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Receive the electric field at the z boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iz,ibuf,n_slices,bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving e along z'

        if (fl%proc/=my_index) then
           ! --- recv data from down in z

           ! Number of slices to communicate along z
           ! (yfu%nzguard for ez, yfu%nzguard-1 for ex and ey)
           n_slices = (3*yfu%nzguard-2)
           bufsize = n_slices*size(yfu%Ez(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%Ez_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the ex, ey and ez slices into that buffer
           call mpi_recv_pack(fl%proc,2,ibuf)
           do iz = yfu%izming,yfu%izmin-1
              yfu%ez(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),ibuf),shape(yfu%Ez(:,:,0)))
              if (yfu%circ_m > 0) &
                   yfu%ez_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%Ez_circ(:,0,:)),ibuf),shape(yfu%Ez_circ(:,0,:)))
           end do
           if (yfu%nzguard>1) then
              do iz = yfu%izming+1,yfu%izmin-1
                 yfu%ex(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ex(:,:,0)),ibuf),shape(yfu%Ex(:,:,0)))
                 if (yfu%circ_m > 0) &
                      yfu%ex_circ(:,iz,:) = &
                      reshape(mpi_unpack_complex_array( size(yfu%Ex_circ(:,0,:)),ibuf),shape(yfu%Ex_circ(:,0,:)))
              end do
              do iz = yfu%izming+1,yfu%izmin-1
                 yfu%ey(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ey(:,:,0)),ibuf),shape(yfu%Ey(:,:,0)))
                 if (yfu%circ_m > 0) &
                      yfu%ey_circ(:,iz,:) = &
                      reshape(mpi_unpack_complex_array( size(yfu%Ey_circ(:,0,:)),ibuf),shape(yfu%Ey_circ(:,0,:)))
              end do
           end if
           
        else if (fu%proc/=my_index) then
           ! --- recv data from up in z
           
           ! Number of slices to communicate along z
           ! (yfu%nzguard for ez, yfu%nzguard-1 for ex and ey)
           n_slices = (3*yfl%nzguard-2)
           bufsize = n_slices*size(yfl%Ez(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%Ez_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the ex, ey and ez slices into that buffer
           call mpi_recv_pack(fu%proc,1,ibuf)
           do iz = yfl%izmax,yfl%izmaxg-1
              yfl%ez(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),ibuf),shape(yfl%Ez(:,:,0)))
              if (yfl%circ_m > 0) &
                   yfl%ez_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfl%Ez_circ(:,0,:)),ibuf),shape(yfl%Ez_circ(:,0,:)))
           end do
           if (yfl%nzguard>1) then
              do iz = yfl%izmax+1,yfl%izmaxg-1
                 yfl%ex(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ex(:,:,0)),ibuf),shape(yfl%Ex(:,:,0)))
                 if (yfl%circ_m > 0) &
                      yfl%ex_circ(:,iz,:) = &
                      reshape(mpi_unpack_complex_array( size(yfl%Ex_circ(:,0,:)),ibuf),shape(yfl%Ex_circ(:,0,:)))
              end do
              do iz = yfl%izmax+1,yfl%izmaxg-1
                 yfl%ey(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ey(:,:,0)),ibuf),shape(yfl%Ey(:,:,0)))
                 if (yfl%circ_m > 0) &
                      yfl%ey_circ(:,iz,:) = &
                      reshape(mpi_unpack_complex_array( size(yfl%Ey_circ(:,0,:)),ibuf),shape(yfl%Ey_circ(:,0,:)))
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

subroutine em3d_exchange_bndb_z(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Send/exchange the magnetic field at the z boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ)   ::ibuf, n_slices, bufsize
#ifdef MPIPARALLEL
  integer(ISZ)   ::iz
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
        if(l_mpiverbose) write(STDOUT,*) '-- sending b along z'

        if (fl%proc/=my_index) then
           ! --- send data down in z

           ! Number of slices to communicate along z
           ! (yfu%nzguard for bx and by, yfu%nzguard-1 for bz)
           n_slices = (3*yfu%nzguard-1)
           bufsize = n_slices*size(yfu%by(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%by_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the bx, by and bz slices into that buffer
           ! bx
           do iz=yfu%izmin,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%bx(:,:,iz),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%bx_circ(:,iz,:),ibuf)
           end do
           ! by
          do iz=yfu%izmin, yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%by(:,:,iz),ibuf)
              if (yfu%circ_m > 0) call mympi_pack(yfu%by_circ(:,iz,:),ibuf)           
           end do
           ! bz
           do iz=yfu%izmin+1,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%bz(:,:,iz),ibuf)
             if (yfu%circ_m > 0) call mympi_pack(yfu%bz_circ(:,iz,:),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)
           
        else if (fu%proc/=my_index) then
           ! --- send data up in z
           
           ! Number of slices to communicate along z
           ! (yfl%nzguard for bx and by, yfl%nzguard-1 for bz)
           n_slices = (3*yfl%nzguard-1)
           bufsize = n_slices*size(yfl%by(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%by_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the bx, by and bz slices into that buffer
           ! bx
           do iz=yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%bx(:,:,iz),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%bx_circ(:,iz,:),ibuf)
           end do
           ! by
           do iz=yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%by(:,:,iz),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%by_circ(:,iz,:),ibuf)
           end do
           ! bz
           do iz=yfl%izmax-yfl%nzguard+1,yfl%izmax-1
              call mympi_pack(yfl%bz(:,:,iz),ibuf)
              if (yfl%circ_m > 0) call mympi_pack(yfl%bz_circ(:,iz,:),ibuf)
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)
        else
#endif
           ! Arrays are on the same processor, no need to send a buffer through mpi
           ! Instead exchange the data directly from array to array
           yfl%bx(:,:,yfl%izmax  :yfl%izmaxg-1) = yfu%bx(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)
           yfl%by(:,:,yfl%izmax  :yfl%izmaxg-1) = yfu%by(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)
           yfl%bz(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%bz(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)

           yfu%bx(:,:,yfu%izming  :yfu%izmin-1) = yfl%bx(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)
           yfu%by(:,:,yfu%izming  :yfu%izmin-1) = yfl%by(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)
           yfu%bz(:,:,yfu%izming+1:yfu%izmin-1) = yfl%bz(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)

           if (yfu%circ_m > 0) then
              yfl%bx_circ(:,yfl%izmax+1:yfl%izmaxg-1,:) = &
                   yfu%bx_circ(:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1,:)
              yfl%by_circ(:,yfl%izmax+1:yfl%izmaxg-1,:) = &
                   yfu%by_circ(:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1,:)
              yfl%bz_circ(:,yfl%izmax  :yfl%izmaxg-1,:) = &
                   yfu%bz_circ(:,yfu%izmin  :yfu%izmin+yfu%nzguard-1,:)
              
              yfu%bx_circ(:,yfu%izming+1:yfu%izmin-1,:) = &
                   yfl%bx_circ(:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1,:)
              yfu%by_circ(:,yfu%izming+1:yfu%izmin-1,:) = &
                   yfl%by_circ(:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1,:)
              yfu%bz_circ(:,yfu%izming  :yfu%izmin-1,:) = &
                   yfl%bz_circ(:,yfl%izmax-yfl%nzguard  :yfl%izmax-1,:)
          endif
           
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
        syfu%bxy(:,:,syfu%izming  :syfu%izmin-1) = 0.
        syfu%bxz(:,:,syfu%izming  :syfu%izmin-1) = yfl%bx(:,:,yfl%izmax-yfl%nzguard  :yfl%izmax-1)*syfu%clight
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
        syfl%bxy(:,:,syfl%izmax  :syfl%izmaxg-1) = 0.
        syfl%bxz(:,:,syfl%izmax  :syfl%izmaxg-1) = yfu%bx(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)*syfl%clight
        syfl%byx(:,:,syfl%izmax  :syfl%izmaxg-1) = 0.
        syfl%byz(:,:,syfl%izmax  :syfl%izmaxg-1) = yfu%by(:,:,yfu%izmin  :yfu%izmin+yfu%nzguard-1)*syfl%clight
        syfl%bzx(:,:,syfl%izmax+1:syfl%izmaxg-1) = yfu%bz(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)*syfl%clight
        syfl%bzy(:,:,syfl%izmax+1:syfl%izmaxg-1) = 0.
     case(splityeefield)
        syfu=>fu%syf
#ifdef MPIPARALLEL
        if (fl%proc/=my_index) then
           ! --- send data down in z
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
subroutine em3d_exchange_bndb_zrecv(fl,fu,ibuf)
  ! ---------------------------------------------------
  ! Receive the magnetic field at the z boundary
  ! ---------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: ibuf,iz,n_slices,bufsize
  
  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving b along z'

        if (fl%proc/=my_index) then
           ! --- recv data from down in z

           ! Number of slices to communicate along z
           ! (yfu%nzguard for bx and by, yfu%nzguard-1 for bz)
           n_slices = (3*yfu%nzguard-1)
           bufsize = n_slices*size(yfu%by(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%by_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the bx, by and bz slices from that buffer
           call mpi_recv_pack(fl%proc,2,ibuf)
           do iz=yfu%izming,yfu%izmin-1
              yfu%bx(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%bx(:,:,iz)),ibuf),shape(yfu%bx(:,:,iz)))
              if (yfu%circ_m > 0) &
                   yfu%bx_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%bx_circ(:,iz,:)),ibuf),shape(yfu%bx_circ(:,iz,:)))
           end do
           do iz=yfu%izming,yfu%izmin-1
              yfu%by(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%by(:,:,iz)),ibuf),shape(yfu%by(:,:,iz)))
              if (yfu%circ_m > 0) &
                   yfu%by_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%by_circ(:,iz,:)),ibuf),shape(yfu%by_circ(:,iz,:)))
           end do
           do iz=yfu%izming+1,yfu%izmin-1
              yfu%bz(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%bz(:,:,iz)),ibuf),shape(yfu%bz(:,:,iz)))
              if (yfu%circ_m > 0) &
                   yfu%bz_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfu%bz_circ(:,iz,:)),ibuf),shape(yfu%bz_circ(:,iz,:)))
           end do
           
        else if (fu%proc/=my_index) then
           ! --- recv data from up in z

           ! Number of slices to communicate along z
           ! (yfl%nzguard for bx and by, yfl%nzguard-1 for bz)
           n_slices = (3*yfl%nzguard-1)
           bufsize = n_slices*size(yfl%by(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%by_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the bx, by and bz slices from that buffer
           call mpi_recv_pack(fu%proc,1,ibuf)
           do iz=yfl%izmax,yfl%izmaxg-1
              yfl%bx(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%bx(:,:,iz)),ibuf),shape(yfl%bx(:,:,iz)))
              if (yfl%circ_m > 0) &
                   yfl%bx_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfl%bx_circ(:,0,:)),ibuf), shape(yfl%bx_circ(:,0,:)))
           end do
           do iz=yfl%izmax,yfl%izmaxg-1
              yfl%by(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%by(:,:,iz)),ibuf),shape(yfl%by(:,:,iz)))
              if (yfl%circ_m > 0) &
                   yfl%by_circ(:,iz,:) = &
                   reshape(mpi_unpack_complex_array( size(yfl%by_circ(:,0,:)),ibuf),shape(yfl%by_circ(:,0,:)))
           end do
           do iz=yfl%izmax+1,yfl%izmaxg-1
              yfl%bz(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%bz(:,:,iz)),ibuf),shape(yfl%bz(:,:,iz)))
              if (yfl%circ_m > 0) &
                   yfl%bz_circ(:,iz,:) = & 
                   reshape(mpi_unpack_complex_array( size(yfl%bz_circ(:,0,:)),ibuf),shape(yfl%bz_circ(:,0,:)))
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

subroutine em3d_exchange_bndf_z(fl,fu,ibuf)
  ! -----------------------------------------------
  ! Send/exchange the field f at the z boundary
  ! -----------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  integer(ISZ)   ::ibuf, n_slices, bufsize
#ifdef MPIPARALLEL
  integer(ISZ)   ::iz
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

        if(l_mpiverbose) write(STDOUT,*) '-- sending f along z'
        
        if (fl%proc/=my_index) then
           ! --- send data down in z

           if (yfu%nzguard>1) then
              ! Number of slices to communicate along z
              n_slices = (yfu%nzguard-1)
              bufsize = n_slices*size(yfu%f(:,:,0))
              ! Check whether to also pack the circ arrays
              if (yfu%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfu%f_circ(:,0,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Pack f
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                 call mympi_pack(yfu%f(:,:,iz),ibuf)
                 if (yfu%circ_m > 0) &
                      call mympi_pack(yfu%f_circ(:,iz,:),ibuf)
              end do
              call mpi_isend_pack(fl%proc,1,ibuf)
           end if
           
        else if (fu%proc/=my_index) then
           ! --- send data up in z
           
           if (yfl%nzguard>1) then
              ! Number of slices to communicate along z
              n_slices = (yfl%nzguard-1)
              bufsize = n_slices*size(yfl%f(:,:,0))
              ! Check whether to also pack the circ arrays
              if (yfl%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfl%f_circ(:,0,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Pack f
              call mpi_packbuffer_init( bufsize, ibuf )
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                 call mympi_pack(yfl%f(:,:,iz),ibuf)
                 if (yfl%circ_m > 0) &
                      call mympi_pack(yfl%f_circ(:,iz,:),ibuf)
              end do
              call mpi_isend_pack(fu%proc,2,ibuf)
           end if
        else
#endif
           yfl%f(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%f(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
           yfu%f(:,:,yfu%izming+1:yfu%izmin-1) = yfl%f(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)

           if (yfu%circ_m > 0) then
              yfl%f_circ(:,yfl%izmax+1:yfl%izmaxg-1,:) = &
                   yfu%f(:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1,:)
              yfu%f_circ(:,yfu%izming+1:yfu%izmin-1,:) = &
                   yfl%f(:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1,:)
           endif

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
subroutine em3d_exchange_bndf_zrecv(fl,fu,ibuf)
  ! --------------------------------------------
  ! Receive the field f at the z boundary
  ! --------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iz,ibuf,n_slices,bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving f along z'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in z

           if (yfu%nzguard>1) then
              ! Number of slices to communicate along z
              n_slices = (yfu%nzguard-1)
              bufsize = n_slices*size(yfu%f(:,:,0))
              ! Check whether to also pack the circ arrays
              if (yfu%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfu%f_circ(:,0,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Receive f
              call mpi_recv_pack(fl%proc,2,ibuf)
              do iz = yfu%izming+1,yfu%izmin-1
                 yfu%f(:,:,iz) = reshape(mpi_unpack_real_array( &
                      size(yfu%F(:,:,0)),ibuf),shape(yfu%F(:,:,0)))
                 if (yfu%circ_m > 0) &
                      yfu%f_circ(:,iz,:) = reshape(mpi_unpack_complex_array( &
                      size(yfu%F_circ(:,0,:)),ibuf),shape(yfu%F_circ(:,0,:)))
              end do
           end if

        else if (fu%proc/=my_index) then
           ! --- recv data from up in z
           
           if (yfl%nzguard>1) then
              ! Number of slices to communicate along z
              n_slices = (yfl%nzguard-1)
              bufsize = n_slices*size(yfl%f(:,:,0))
              ! Check whether to also pack the circ arrays
              if (yfl%circ_m > 0) &
                   ! Factor 2 since a complex takes up twice more space
                   bufsize = bufsize + 2*n_slices*size(yfl%f_circ(:,0,:))
              ! Allocate a buffer array in mpibuffer
              call mpi_packbuffer_init( bufsize, ibuf )
              ! Receive f
              call mpi_recv_pack(fu%proc,1,ibuf)
              do iz = yfl%izmax+1,yfl%izmaxg-1
                 yfl%f(:,:,iz) = reshape(mpi_unpack_real_array( &
                      size(yfl%F(:,:,0)),ibuf),shape(yfl%F(:,:,0)))
                 if (yfl%circ_m > 0) &
                      yfl%f_circ(:,iz,:) = reshape(mpi_unpack_complex_array( &
                      size(yfl%F_circ(:,0,:)),ibuf),shape(yfl%F_circ(:,0,:)))
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

subroutine em3d_exchange_bndj_z(fl,fu,ibuf)
  ! ---------------------------------------------
  ! Send/exchange the current at the z boundary
  ! ---------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iz,ibuf,nguardinl,nguardinu,n_slices,bufsize

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

        if(l_mpiverbose) write(STDOUT,*) '-- sending j along z'
        
        if (fl%proc/=my_index) then
           ! --- send data down in z

           ! Number of slices to communicate along z
           nguardinu = yfu%nzguard
           n_slices = 3*(yfu%nzguard+nguardinu)+2
           bufsize = n_slices*size(yfu%J(:,:,0,1))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%J_circ(:,0,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the Jx, Jy and Jz slices into that buffer
           ! Jx
           do iz = -yfu%nzguard,nguardinu
              call mympi_pack(yfu%J(:,:,iz,1),ibuf)
              if (yfu%circ_m > 0) &
                   call mympi_pack(yfu%J_circ(:,iz,1,:),ibuf)
           end do
           ! Jy
           do iz = -yfu%nzguard,nguardinu
              call mympi_pack(yfu%J(:,:,iz,2),ibuf)
              if (yfu%circ_m > 0) &
                   call mympi_pack(yfu%J_circ(:,iz,2,:),ibuf)
           end do
           ! Jz
           do iz = -yfu%nzguard,nguardinu-1
              call mympi_pack(yfu%J(:,:,iz,3),ibuf)
              if (yfu%circ_m > 0) &
                   call mympi_pack(yfu%J_circ(:,iz,3,:),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then
           ! --- send data up in z

           ! Number of slices to communicate along z
           nguardinl = yfl%nzguard
           n_slices = 3*(yfl%nzguard+nguardinl)+2
           bufsize = n_slices*size(yfl%J(:,:,0,1))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%J_circ(:,0,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively pack the Jx, Jy and Jz slices into that buffer
           ! Jx
           do iz = yfl%nz-nguardinl, yfl%nz+yfl%nzguard
              call mympi_pack(yfl%J(:,:,iz,1),ibuf)
              if (yfl%circ_m > 0) &
                   call mympi_pack(yfl%J_circ(:,iz,1,:),ibuf)
           end do
           do iz = yfl%nz-nguardinl, yfl%nz+yfl%nzguard
              call mympi_pack(yfl%J(:,:,iz,2),ibuf)
              if (yfl%circ_m > 0) &
                   call mympi_pack(yfl%J_circ(:,iz,2,:),ibuf)
           end do
           do iz = yfl%nz-nguardinl, yfl%nz+yfl%nzguard-1
              call mympi_pack(yfl%J(:,:,iz,3),ibuf)
              if (yfl%circ_m > 0) &
                   call mympi_pack(yfl%J_circ(:,iz,3,:),ibuf)
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)

        else
#endif
           ! Arrays are on the same processor, no need to send a buffer through mpi
           ! Instead exchange the data directly from array to array
           nguardinl = yfl%nzguard
           nguardinu = yfu%nzguard

           yfu%J(:,:,-nguardinu:yfu%nzguard,1:2) = yfu%J(:,:,-nguardinu:yfu%nzguard,1:2) &
                + yfl%J(:,:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard,1:2)
           yfu%J(:,:,-nguardinu:yfu%nzguard-1,3) = yfu%J(:,:,-nguardinu:yfu%nzguard-1,3) &
                + yfl%J(:,:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard-1,3) 
           yfl%J(:,:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard,1:2) = yfu%J(:,:,-nguardinu:yfu%nzguard,1:2)
           yfl%J(:,:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard-1,3) = yfu%J(:,:,-nguardinu:yfu%nzguard-1,3)

           if (yfu%circ_m > 0) then
              yfu%J_circ(:,-nguardinu:yfu%nzguard,1:2,:) = yfu%J_circ(:,-nguardinu:yfu%nzguard,1:2,:) &
                   + yfl%J_circ(:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard,1:2,:)
              yfu%J_circ(:,-nguardinu:yfu%nzguard-1,3,:) = yfu%J_circ(:,-nguardinu:yfu%nzguard-1,3,:) &
                   + yfl%J_circ(:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard-1,3,:) 
              yfl%J_circ(:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard,1:2,:) = &
                   yfu%J_circ(:,-nguardinu:yfu%nzguard,1:2,:)
              yfl%J_circ(:,yfl%nz-nguardinl:yfl%nz+yfl%nzguard-1,3,:) = &
                   yfu%J_circ(:,-nguardinu:yfu%nzguard-1,3,:)
           endif
           
#ifdef MPIPARALLEL
        end if
#endif
     end select
  end select
  return
end subroutine em3d_exchange_bndj_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_zrecv(fl,fu,ibuf)
  ! -----------------------------------------
  ! Receive the current at the z boundary
  ! -----------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iz,ibuf,nguardinl,nguardinu, n_slices, bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf
        
        if(l_mpiverbose) write(STDOUT,*) '-- receiving j along z'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in z

           ! Number of slices to communicate along z
           nguardinu = yfu%nzguard
           n_slices = (3*(yfu%nzguard+nguardinu) + 2)
           bufsize = n_slices*size(yfu%J(:,:,0,1))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%J_circ(:,0,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the Jx, Jy and Jz slices from that buffer
           call mpi_recv_pack(fl%proc,2,ibuf)
           do iz = -nguardinu,yfu%nzguard
              yfu%J(:,:,iz,1) = yfu%J(:,:,iz,1) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),ibuf), &
                   shape(yfu%J(:,:,0,1)))
              if ( yfu%circ_m > 0 ) &
                   yfu%J_circ(:,iz,1,:) = yfu%J_circ(:,iz,1,:) + reshape( mpi_unpack_complex_array( &
                   size(yfu%J_circ(:,0,1,:)),ibuf), shape(yfu%J_circ(:,0,1,:)) )
           end do
           do iz = -nguardinu,yfu%nzguard
              yfu%J(:,:,iz,2) = yfu%J(:,:,iz,2) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),ibuf), &
                   shape(yfu%J(:,:,0,1)))
              if ( yfu%circ_m > 0 ) &
                   yfu%J_circ(:,iz,2,:) = yfu%J_circ(:,iz,2,:) + reshape( mpi_unpack_complex_array( &
                   size(yfu%J_circ(:,0,2,:)),ibuf), shape(yfu%J_circ(:,0,2,:)) )
           end do
           do iz = -nguardinu,yfu%nzguard-1
              yfu%J(:,:,iz,3) = yfu%J(:,:,iz,3) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),ibuf), &
                   shape(yfu%J(:,:,0,1)))
              if ( yfu%circ_m > 0 ) &
                   yfu%J_circ(:,iz,3,:) = yfu%J_circ(:,iz,3,:) + reshape( mpi_unpack_complex_array( &
                   size(yfu%J_circ(:,0,3,:)),ibuf), shape(yfu%J_circ(:,0,3,:)) )
           end do

        else if (fu%proc/=my_index) then
           ! --- recv data from up in z
           
           ! Number of slices to communicate along z
           nguardinl = yfl%nzguard
           n_slices = (3*(yfl%nzguard+nguardinl) + 2)
           bufsize = n_slices*size(yfl%J(:,:,0,1))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%J_circ(:,0,1,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Successively receive the Jx, Jy and Jz slices from that buffer
           call mpi_recv_pack(fu%proc,1,ibuf)
           do iz = -yfl%nzguard,nguardinl
              yfl%J(:,:,yfl%nz+iz,1) = yfl%J(:,:,yfl%nz+iz,1) + &
                   reshape(mpi_unpack_real_array( &
                   size(yfl%J(:,:,yfl%nz-1,1)),ibuf), shape(yfl%J(:,:,yfl%nz-1,1)))
              if ( yfl%circ_m > 0 ) &
                   yfl%J_circ(:,yfl%nz+iz,1,:) = yfl%J_circ(:,yfl%nz+iz,1,:) + &
                   reshape( mpi_unpack_complex_array( &
                   size(yfl%J_circ(:,0,1,:)),ibuf), shape(yfl%J_circ(:,0,1,:)) )
           end do
           do iz = -yfl%nzguard,nguardinl
              yfl%J(:,:,yfl%nz+iz,2) = yfl%J(:,:,yfl%nz+iz,2) + &
                   reshape(mpi_unpack_real_array( &
                   size(yfl%J(:,:,yfl%nz-1,2)),ibuf), shape(yfl%J(:,:,yfl%nz-1,2)))
              if ( yfl%circ_m > 0 ) &
                   yfl%J_circ(:,yfl%nz+iz,2,:) = yfl%J_circ(:,yfl%nz+iz,2,:) + &
                   reshape( mpi_unpack_complex_array( &
                   size(yfl%J_circ(:,0,2,:)),ibuf), shape(yfl%J_circ(:,0,2,:)) )
           end do
           do iz = -yfl%nzguard,nguardinl-1
              yfl%J(:,:,yfl%nz+iz,3) = yfl%J(:,:,yfl%nz+iz,3) + &
                   reshape(mpi_unpack_real_array( &
                   size(yfl%J(:,:,yfl%nz-1,3)),ibuf), shape(yfl%J(:,:,yfl%nz-1,3)))
              if ( yfl%circ_m > 0 ) &
                   yfl%J_circ(:,yfl%nz+iz,3,:) = yfl%J_circ(:,yfl%nz+iz,3,:) + &
                   reshape( mpi_unpack_complex_array( &
                   size(yfl%J_circ(:,0,3,:)),ibuf), shape(yfl%J_circ(:,0,3,:)) )
           end do
        end if
     end select
  end select

  return
end subroutine em3d_exchange_bndj_zrecv
#endif

subroutine em3d_exchange_bndrho_z(fl,fu,ibuf)
  ! ------------------------------------------------
  ! Send/exchange the field rho at the z boundary
  ! ------------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iz,ibuf,nguardinl,nguardinu, n_slices, bufsize

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
        if(l_mpiverbose) write(STDOUT,*) '-- sending rho along z'
        
        if (fl%proc/=my_index) then
           ! --- send data down in z

           ! Number of slices to communicate along z
           nguardinu = yfu%nzguard
           n_slices =  yfu%nzguard + nguardinu + 1
           bufsize = n_slices*size(yfu%rho(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%rho_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Pack rho
           do iz = -yfu%nzguard,nguardinu
              call mympi_pack(yfu%rho(:,:,iz  ),ibuf)
              if (yfu%circ_m > 0) &
                   call mympi_pack(yfu%rho_circ(:,iz,:),ibuf)
           end do
           call mpi_isend_pack(fl%proc,1,ibuf)

        else if (fu%proc/=my_index) then
           ! --- send data up in z
           
           ! Number of slices to communicate along z
           nguardinl = yfl%nzguard
           n_slices =  yfl%nzguard + nguardinl + 1
           bufsize = n_slices*size(yfl%rho(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%rho_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Pack rho
           do iz = -nguardinl,yfl%nzguard
              call mympi_pack(yfl%rho(:,:,yfl%nz+iz  ),ibuf)
              if (yfl%circ_m > 0) &
                   call mympi_pack(yfl%rho_circ(:,yfl%nz+iz,:),ibuf)
           end do
           call mpi_isend_pack(fu%proc,2,ibuf)

        else
#endif
           nguardinl = yfl%nzguard
           nguardinu = yfu%nzguard

           yfu%Rho(:,:,-nguardinu:yfu%nzguard) = yfu%Rho(:,:,-nguardinu:yfu%nzguard) &
                + yfl%Rho(:,:,yfl%nz-nguardinl:yfl%nz+yfu%nzguard)
           yfl%Rho(:,:,yfl%nz-nguardinl:yfl%nz+yfu%nzguard) = yfu%Rho(:,:,-nguardinu:yfu%nzguard)

           if (yfu%circ_m > 0) then
              yfu%Rho_circ(:,-nguardinu:yfu%nzguard,:) = yfu%Rho_circ(:,-nguardinu:yfu%nzguard,:) &
                   + yfl%Rho_circ(:,yfl%nz-nguardinl:yfl%nz+yfu%nzguard,:)
              yfl%Rho_circ(:,yfl%nz-nguardinl:yfl%nz+yfu%nzguard,:) = &
                   yfu%Rho_circ(:,-nguardinu:yfu%nzguard,:)
           endif
           
#ifdef MPIPARALLEL
        end if
#endif
     end select
  end select
  return
end subroutine em3d_exchange_bndrho_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_zrecv(fl,fu,ibuf)
  ! --------------------------------------------
  ! Receive the field rho at the z boundary
  ! --------------------------------------------
  use mod_emfield3d
  implicit none

  TYPE(EM3D_FIELDtype) :: fl, fu
  TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
  TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
  integer(ISZ) :: iz, ibuf,nguardinl,nguardinu, n_slices, bufsize

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
  case(yeefield)
     yfl=>fl%yf
     select case(fu%fieldtype)
     case(yeefield)
        yfu=>fu%yf

        if(l_mpiverbose) write(STDOUT,*) '-- receiving rho along z'
        
        if (fl%proc/=my_index) then
           ! --- recv data from down in z
           
           ! Number of slices to communicate along z
           nguardinu = yfu%nzguard
           n_slices =  yfu%nzguard + nguardinu + 1
           bufsize = n_slices*size(yfu%rho(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfu%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfu%rho_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Receive rho
           call mpi_recv_pack(fl%proc,2,ibuf)
           do iz = -nguardinu,yfu%nzguard
              yfu%rho(:,:,iz  ) = yfu%rho(:,:,iz) + reshape(mpi_unpack_real_array( &
                   size(yfu%rho(:,:,0)),ibuf), shape(yfu%rho(:,:,0)))
              if (yfu%circ_m > 0) &
                   yfu%rho_circ(:,iz,: ) = yfu%rho_circ(:,iz,:) + reshape(mpi_unpack_complex_array( &
                   size(yfu%rho_circ(:,0,:)),ibuf), shape(yfu%rho_circ(:,0,:)))
           end do

        else if (fu%proc/=my_index) then
           ! --- recv data from up in z

           ! Number of slices to communicate along z
           nguardinl = yfl%nzguard
           n_slices =  yfl%nzguard + nguardinl + 1
           bufsize = n_slices*size(yfl%rho(:,:,0))
           ! Check whether to also pack the circ arrays
           if (yfl%circ_m > 0) &
                ! Factor 2 since a complex takes up twice more space
                bufsize = bufsize + 2*n_slices*size(yfl%rho_circ(:,0,:))
           ! Allocate a buffer array in mpibuffer
           call mpi_packbuffer_init( bufsize, ibuf )
           ! Receive rho
           call mpi_recv_pack(fu%proc,1,ibuf)
           do iz = -yfl%nzguard,nguardinl
              yfl%rho(:,:,yfl%nz+iz  ) = yfl%rho(:,:,yfl%nz+iz  ) + reshape( &
                   mpi_unpack_real_array(size(yfl%rho(:,:,yfl%nz)),ibuf), shape(yfl%rho(:,:,yfl%nz  )))
              if (yfl%circ_m > 0) &
                   yfl%rho_circ(:,yfl%nz+iz,: ) = yfl%rho_circ(:,yfl%nz+iz,:) + reshape( &
                   mpi_unpack_complex_array(size(yfl%rho_circ(:,yfl%nz,:)),ibuf), &
                   shape(yfl%rho_circ(:,yfl%nz,:)))
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
  integer(ISZ) :: ibuf

  ibuf = 2

  ! --- X
  if (.not.b%core%yf%l_1dz) then
     ! core<--->sides
     call em3d_exchange_bnde_x(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%sidexl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%sidexl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! sides<--->edges
     call em3d_exchange_bnde_x(b%sideyl,   b%edgexryl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%edgexlyl, b%sideyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%edgexlyl, b%sideyl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%sideyl,   b%edgexryl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_x(b%sideyr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%edgexlyr, b%sideyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%edgexlyr, b%sideyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%sideyr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_x(b%sidezl,   b%edgexrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%edgexlzl, b%sidezl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%edgexlzl, b%sidezl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%sidezl,   b%edgexrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_x(b%sidezr,   b%edgexrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%edgexlzr, b%sidezr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%edgexlzr, b%sidezr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%sidezr,   b%edgexrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! edges<--->corners
     call em3d_exchange_bnde_x(b%edgeylzl,     b%cornerxrylzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%cornerxlylzl, b%edgeylzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%cornerxlylzl, b%edgeylzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%edgeylzl,     b%cornerxrylzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_x(b%edgeyrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%cornerxlyrzl, b%edgeyrzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%cornerxlyrzl, b%edgeyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%edgeyrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_x(b%edgeylzr,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%cornerxlylzr, b%edgeylzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%cornerxlylzr, b%edgeylzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%edgeylzr,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_x(b%edgeyrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_x(b%cornerxlyrzr, b%edgeyrzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_xrecv(b%cornerxlyrzr, b%edgeyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_xrecv(b%edgeyrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Y
  if (.not.b%core%yf%l_2dxz) then
     ! core<--->sides
     call em3d_exchange_bnde_y(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%sideyl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%sideyl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! sides<--->edges
     call em3d_exchange_bnde_y(b%sidexl,   b%edgexlyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%edgexlyl, b%sidexl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%edgexlyl, b%sidexl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%sidexl,   b%edgexlyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_y(b%sidexr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%edgexryl, b%sidexr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%edgexryl, b%sidexr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%sidexr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_y(b%sidezl,   b%edgeyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%edgeylzl, b%sidezl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%edgeylzl, b%sidezl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%sidezl,   b%edgeyrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_y(b%sidezr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%edgeylzr, b%sidezr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%edgeylzr, b%sidezr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%sidezr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! edges<--->corners
     call em3d_exchange_bnde_y(b%edgexlzl,     b%cornerxlyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%cornerxlylzl, b%edgexlzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%cornerxlylzl, b%edgexlzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%edgexlzl,     b%cornerxlyrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_y(b%edgexrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%cornerxrylzl, b%edgexrzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%cornerxrylzl, b%edgexrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%edgexrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_y(b%edgexlzr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%cornerxlylzr, b%edgexlzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%cornerxlylzr, b%edgexlzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%edgexlzr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bnde_y(b%edgexrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_y(b%cornerxrylzr, b%edgexrzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bnde_yrecv(b%cornerxrylzr, b%edgexrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bnde_yrecv(b%edgexrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bnde_z(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%sidezl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%sidezl, b%core, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bnde_z(b%sidexl,   b%edgexlzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%edgexlzl, b%sidexl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgexlzl, b%sidexl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%sidexl,   b%edgexlzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sidexr,   b%edgexrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%edgexrzl, b%sidexr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgexrzl, b%sidexr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%sidexr,   b%edgexrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sideyl,   b%edgeylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%edgeylzl, b%sideyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgeylzl, b%sideyl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%sideyl,   b%edgeylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sideyr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%edgeyrzl, b%sideyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgeyrzl, b%sideyr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%sideyr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bnde_z(b%edgexlyl,     b%cornerxlylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%cornerxlylzl, b%edgexlyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxlylzl, b%edgexlyl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%edgexlyl,     b%cornerxlylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexryl,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%cornerxrylzl, b%edgexryl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxrylzl, b%edgexryl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%edgexryl,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexlyr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%cornerxlyrzl, b%edgexlyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxlyrzl, b%edgexlyr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%edgexlyr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexryr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_z(b%cornerxryrzl, b%edgexryr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxryrzl, b%edgexryr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bnde_zrecv(b%edgexryr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_e

subroutine em3d_exchange_b(b)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_BLOCKtype) :: b
  integer(ISZ) :: ibuf

!  if (b%core%yf%spectral) return

  ibuf = 200

  ! --- X
  if (.not.b%core%yf%l_1dz) then
     ! core<--->sides
     call em3d_exchange_bndb_x(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%sidexl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%sidexl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! sides<--->edges
     call em3d_exchange_bndb_x(b%sideyl,   b%edgexryl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%edgexlyl, b%sideyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%edgexlyl, b%sideyl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%sideyl,   b%edgexryl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_x(b%sideyr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%edgexlyr, b%sideyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%edgexlyr, b%sideyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%sideyr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_x(b%sidezl,   b%edgexrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%edgexlzl, b%sidezl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%edgexlzl, b%sidezl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%sidezl,   b%edgexrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_x(b%sidezr,   b%edgexrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%edgexlzr, b%sidezr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%edgexlzr, b%sidezr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%sidezr,   b%edgexrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! edges<--->corners
     call em3d_exchange_bndb_x(b%edgeylzl,     b%cornerxrylzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%cornerxlylzl, b%edgeylzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%cornerxlylzl, b%edgeylzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%edgeylzl,     b%cornerxrylzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_x(b%edgeyrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%cornerxlyrzl, b%edgeyrzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%cornerxlyrzl, b%edgeyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%edgeyrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_x(b%edgeylzr,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%cornerxlylzr, b%edgeylzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%cornerxlylzr, b%edgeylzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%edgeylzr,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_x(b%edgeyrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_x(b%cornerxlyrzr, b%edgeyrzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_xrecv(b%cornerxlyrzr, b%edgeyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_xrecv(b%edgeyrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Y
  if (.not.b%core%yf%l_2dxz) then
     ! core<--->sides
     call em3d_exchange_bndb_y(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%sideyl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%sideyl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! sides<--->edges
     call em3d_exchange_bndb_y(b%sidexl,   b%edgexlyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%edgexlyl, b%sidexl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%edgexlyl, b%sidexl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%sidexl,   b%edgexlyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_y(b%sidexr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%edgexryl, b%sidexr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%edgexryl, b%sidexr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%sidexr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_y(b%sidezl,   b%edgeyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%edgeylzl, b%sidezl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%edgeylzl, b%sidezl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%sidezl,   b%edgeyrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_y(b%sidezr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%edgeylzr, b%sidezr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%edgeylzr, b%sidezr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%sidezr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! edges<--->corners
     call em3d_exchange_bndb_y(b%edgexlzl,     b%cornerxlyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%cornerxlylzl, b%edgexlzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%cornerxlylzl, b%edgexlzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%edgexlzl,     b%cornerxlyrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_y(b%edgexrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%cornerxrylzl, b%edgexrzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%cornerxrylzl, b%edgexrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%edgexrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_y(b%edgexlzr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%cornerxlylzr, b%edgexlzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%cornerxlylzr, b%edgexlzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%edgexlzr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndb_y(b%edgexrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_y(b%cornerxrylzr, b%edgexrzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndb_yrecv(b%cornerxrylzr, b%edgexrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndb_yrecv(b%edgexrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndb_z(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%sidezl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%sidezl, b%core, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndb_z(b%sidexl,   b%edgexlzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%edgexlzl, b%sidexl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgexlzl, b%sidexl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%sidexl,   b%edgexlzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sidexr,   b%edgexrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%edgexrzl, b%sidexr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgexrzl, b%sidexr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%sidexr,   b%edgexrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sideyl,   b%edgeylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%edgeylzl, b%sideyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgeylzl, b%sideyl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%sideyl,   b%edgeylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sideyr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%edgeyrzl, b%sideyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgeyrzl, b%sideyr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%sideyr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndb_z(b%edgexlyl,     b%cornerxlylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%cornerxlylzl, b%edgexlyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxlylzl, b%edgexlyl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%edgexlyl,     b%cornerxlylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexryl,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%cornerxrylzl, b%edgexryl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxrylzl, b%edgexryl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%edgexryl,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexlyr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%cornerxlyrzl, b%edgexlyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxlyrzl, b%edgexlyr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%edgexlyr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexryr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_z(b%cornerxryrzl, b%edgexryr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxryrzl, b%edgexryr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndb_zrecv(b%edgexryr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_b

subroutine em3d_exchange_f(b)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_BLOCKtype) :: b
  integer(ISZ) :: ibuf

  ibuf = 400

  ! --- X
  if (.not.b%core%yf%l_1dz) then
     ! core<--->sides
     call em3d_exchange_bndf_x(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%sidexl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%sidexl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! sides<--->edges
     call em3d_exchange_bndf_x(b%sideyl,   b%edgexryl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%edgexlyl, b%sideyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%edgexlyl, b%sideyl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%sideyl,   b%edgexryl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_x(b%sideyr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%edgexlyr, b%sideyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%edgexlyr, b%sideyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%sideyr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_x(b%sidezl,   b%edgexrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%edgexlzl, b%sidezl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%edgexlzl, b%sidezl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%sidezl,   b%edgexrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_x(b%sidezr,   b%edgexrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%edgexlzr, b%sidezr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%edgexlzr, b%sidezr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%sidezr,   b%edgexrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! edges<--->corners
     call em3d_exchange_bndf_x(b%edgeylzl,     b%cornerxrylzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%cornerxlylzl, b%edgeylzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%cornerxlylzl, b%edgeylzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%edgeylzl,     b%cornerxrylzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_x(b%edgeyrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%cornerxlyrzl, b%edgeyrzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%cornerxlyrzl, b%edgeyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%edgeyrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_x(b%edgeylzr,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%cornerxlylzr, b%edgeylzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%cornerxlylzr, b%edgeylzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%edgeylzr,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_x(b%edgeyrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_x(b%cornerxlyrzr, b%edgeyrzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_xrecv(b%cornerxlyrzr, b%edgeyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_xrecv(b%edgeyrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Y
  if (.not.b%core%yf%l_2dxz) then
     ! core<--->sides
     call em3d_exchange_bndf_y(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%sideyl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%sideyl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! sides<--->edges
     call em3d_exchange_bndf_y(b%sidexl,   b%edgexlyr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%edgexlyl, b%sidexl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%edgexlyl, b%sidexl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%sidexl,   b%edgexlyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_y(b%sidexr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%edgexryl, b%sidexr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%edgexryl, b%sidexr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%sidexr,   b%edgexryr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_y(b%sidezl,   b%edgeyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%edgeylzl, b%sidezl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%edgeylzl, b%sidezl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%sidezl,   b%edgeyrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_y(b%sidezr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%edgeylzr, b%sidezr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%edgeylzr, b%sidezr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%sidezr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     ! edges<--->corners
     call em3d_exchange_bndf_y(b%edgexlzl,     b%cornerxlyrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%cornerxlylzl, b%edgexlzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%cornerxlylzl, b%edgexlzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%edgexlzl,     b%cornerxlyrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_y(b%edgexrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%cornerxrylzl, b%edgexrzl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%cornerxrylzl, b%edgexrzl, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%edgexrzl,     b%cornerxryrzl, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_y(b%edgexlzr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%cornerxlylzr, b%edgexlzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%cornerxlylzr, b%edgexlzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%edgexlzr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
     call em3d_exchange_bndf_y(b%edgexrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_y(b%cornerxrylzr, b%edgexrzr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     call em3d_exchange_bndf_yrecv(b%cornerxrylzr, b%edgexrzr, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndf_yrecv(b%edgexrzr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndf_z(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%sidezl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%sidezl, b%core, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndf_z(b%sidexl,   b%edgexlzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%edgexlzl, b%sidexl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgexlzl, b%sidexl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%sidexl,   b%edgexlzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%sidexr,   b%edgexrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%edgexrzl, b%sidexr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgexrzl, b%sidexr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%sidexr,   b%edgexrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%sideyl,   b%edgeylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%edgeylzl, b%sideyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgeylzl, b%sideyl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%sideyl,   b%edgeylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%sideyr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%edgeyrzl, b%sideyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%edgeyrzl, b%sideyr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%sideyr,   b%edgeyrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndf_z(b%edgexlyl,     b%cornerxlylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%cornerxlylzl, b%edgexlyl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxlylzl, b%edgexlyl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%edgexlyl,     b%cornerxlylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%edgexryl,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%cornerxrylzl, b%edgexryl, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxrylzl, b%edgexryl, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%edgexryl,     b%cornerxrylzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%edgexlyr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%cornerxlyrzl, b%edgexlyr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxlyrzl, b%edgexlyr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%edgexlyr,     b%cornerxlyrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndf_z(b%edgexryr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_z(b%cornerxryrzl, b%edgexryr, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  call em3d_exchange_bndf_zrecv(b%cornerxryrzl, b%edgexryr, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndf_zrecv(b%edgexryr,     b%cornerxryrzr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_f

subroutine em3d_exchange_j(b)
  use mod_emfield3d
  implicit none
  TYPE(EM3D_BLOCKtype) :: b
  integer(ISZ) :: ibuf
#ifdef MPIPARALLEL
  integer(MPIISZ)::mpirequest(2)
#endif

  ibuf = 600

  ! --- X
  if (.not.b%core%yf%l_1dz) then
     ! core<--->sides
     call em3d_exchange_bndj_x(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     if(b%xrbnd /= periodic) call em3d_exchange_bndj_x(b%sidexl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     if(b%xrbnd /= periodic) call em3d_exchange_bndj_xrecv(b%sidexl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndj_xrecv(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Y
  if (.not.b%core%yf%l_2dxz) then
     ! core<--->sides
     call em3d_exchange_bndj_y(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     if(b%yrbnd /= periodic) call em3d_exchange_bndj_y(b%sideyl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     if(b%yrbnd /= periodic) call em3d_exchange_bndj_yrecv(b%sideyl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndj_yrecv(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndj_z(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  if(b%zrbnd /= periodic) call em3d_exchange_bndj_z(b%sidezl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  if(b%zrbnd /= periodic) call em3d_exchange_bndj_zrecv(b%sidezl, b%core, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndj_zrecv(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_j

subroutine em3d_exchange_rho(b)
  use mod_emfield3d
  implicit none

  TYPE(EM3D_BLOCKtype) :: b
  integer(ISZ) :: ibuf

  ibuf = 800

  ! --- X
  if (.not.b%core%yf%l_1dz) then
     ! core<--->sides
     call em3d_exchange_bndrho_x(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     if(b%xrbnd /= periodic) call em3d_exchange_bndrho_x(b%sidexl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     if(b%xrbnd /= periodic) call em3d_exchange_bndrho_xrecv(b%sidexl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndrho_xrecv(b%core,   b%sidexr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Y
  if (.not.b%core%yf%l_2dxz) then
     ! core<--->sides
     call em3d_exchange_bndrho_y(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     if(b%yrbnd /= periodic) call em3d_exchange_bndrho_y(b%sideyl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
     if(b%yrbnd /= periodic) call em3d_exchange_bndrho_yrecv(b%sideyl, b%core, ibuf); ibuf=ibuf+1
     call em3d_exchange_bndrho_yrecv(b%core,   b%sideyr, ibuf); ibuf=ibuf+1
     call mpi_waitall_requests()
#endif
  endif

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndrho_z(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
  if(b%zrbnd /= periodic) call em3d_exchange_bndrho_z(b%sidezl, b%core, ibuf); ibuf=ibuf+1
#ifdef MPIPARALLEL
  if(b%zrbnd /= periodic) call em3d_exchange_bndrho_zrecv(b%sidezl, b%core, ibuf); ibuf=ibuf+1
  call em3d_exchange_bndrho_zrecv(b%core,   b%sidezr, ibuf); ibuf=ibuf+1
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

  if (f%l_nodecentered) return

  if (.not.f%l_2dxz) then
     do l=-f%nzguard,f%nz+f%nzguard
        do k=-f%nyguard,f%ny+f%nyguard
           do j=f%nx+f%nxguard-1,-f%nxguard+1,-1
              f%exp(j,k,l)=0.5*(f%exp(j,k,l)+f%exp(j-1,k,l))
              f%byp(j,k,l)=0.5*(f%byp(j,k,l)+f%byp(j-1,k,l))
              f%bzp(j,k,l)=0.5*(f%bzp(j,k,l)+f%bzp(j-1,k,l))
           enddo
        enddo
     enddo

     do l=-f%nzguard,f%nz+f%nzguard
        do k=f%ny+f%nyguard-1,-f%nyguard+1,-1
           do j=-f%nxguard,f%nx+f%nxguard
              f%eyp(j,k,l)=0.5*(f%eyp(j,k,l)+f%eyp(j,k-1,l))
              f%bzp(j,k,l)=0.5*(f%bzp(j,k,l)+f%bzp(j,k-1,l))
              f%bxp(j,k,l)=0.5*(f%bxp(j,k,l)+f%bxp(j,k-1,l))
           enddo
        enddo
     enddo

     do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
        do k=-f%nyguard,f%ny+f%nyguard
           do j=-f%nxguard,f%nx+f%nxguard
              f%ezp(j,k,l)=0.5*(f%ezp(j,k,l)+f%ezp(j,k,l-1))
              f%bxp(j,k,l)=0.5*(f%bxp(j,k,l)+f%bxp(j,k,l-1))
              f%byp(j,k,l)=0.5*(f%byp(j,k,l)+f%byp(j,k,l-1))
           enddo
        enddo
     enddo

  else

     if (f%l_1dz) then
        j = 0
        k = 0
        do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
           f%ezp(j,k,l)=0.5*(f%ezp(j,k,l)+f%ezp(j,k,l-1))
           f%bxp(j,k,l)=0.5*(f%bxp(j,k,l)+f%bxp(j,k,l-1))
           f%byp(j,k,l)=0.5*(f%byp(j,k,l)+f%byp(j,k,l-1))
        enddo
     else
        k = 0
        do l=-f%nzguard,f%nz+f%nzguard
           do j=f%nx+f%nxguard-1,-f%nxguard+1,-1
              f%exp(j,k,l)=0.5*(f%exp(j,k,l)+f%exp(j-1,k,l))
              f%byp(j,k,l)=0.5*(f%byp(j,k,l)+f%byp(j-1,k,l))
              f%bzp(j,k,l)=0.5*(f%bzp(j,k,l)+f%bzp(j-1,k,l))
           enddo
        enddo

        do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
           do j=-f%nxguard,f%nx+f%nxguard
              f%ezp(j,k,l)=0.5*(f%ezp(j,k,l)+f%ezp(j,k,l-1))
              f%bxp(j,k,l)=0.5*(f%bxp(j,k,l)+f%bxp(j,k,l-1))
              f%byp(j,k,l)=0.5*(f%byp(j,k,l)+f%byp(j,k,l-1))
           enddo
        enddo

        if (f%circ_m>0) then
           do l=-f%nzguard,f%nz+f%nzguard
              do j=f%nx+f%nxguard-1,-f%nxguard+1,-1
                 f%exp_circ(j,l,:)=0.5*(f%exp_circ(j,l,:)+f%exp_circ(j-1,l,:))
                 f%byp_circ(j,l,:)=0.5*(f%byp_circ(j,l,:)+f%byp_circ(j-1,l,:))
                 f%bzp_circ(j,l,:)=0.5*(f%bzp_circ(j,l,:)+f%bzp_circ(j-1,l,:))
              enddo
           enddo

           do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
              do j=-f%nxguard,f%nx+f%nxguard
                 f%ezp_circ(j,l,:)=0.5*(f%ezp_circ(j,l,:)+f%ezp_circ(j,l-1,:))
                 f%bxp_circ(j,l,:)=0.5*(f%bxp_circ(j,l,:)+f%bxp_circ(j,l-1,:))
                 f%byp_circ(j,l,:)=0.5*(f%byp_circ(j,l,:)+f%byp_circ(j,l-1,:))
              enddo
           enddo
        end if

     endif
  endif

  f%l_nodecentered = .true.

  return
end subroutine yee2node3d

subroutine node2yee3d(f)
  ! puts EM field back from node to Yee grid
  use mod_emfield3d
  implicit none
  TYPE(EM3D_YEEFIELDtype) :: f

  INTEGER :: j,k,l
  !return

  if (.not.f%l_nodecentered) return

  if (.not.f%l_2dxz) then

     do l=-f%nzguard+1,f%nz+f%nzguard-1
        do k=-f%nyguard,f%ny+f%nyguard
           do j=-f%nxguard,f%nx+f%nxguard
              f%ezp(j,k,l)=2.*f%ezp(j,k,l)-f%ezp(j,k,l-1)
              f%bxp(j,k,l)=2.*f%bxp(j,k,l)-f%bxp(j,k,l-1)
              f%byp(j,k,l)=2.*f%byp(j,k,l)-f%byp(j,k,l-1)
           enddo
        enddo
     enddo

     do l=-f%nzguard,f%nz+f%nzguard
        do k=-f%nyguard+1,f%ny+f%nyguard-1
           do j=-f%nxguard,f%nx+f%nxguard
              f%eyp(j,k,l)=2.*f%eyp(j,k,l)-f%eyp(j,k-1,l)
              f%bzp(j,k,l)=2.*f%bzp(j,k,l)-f%bzp(j,k-1,l)
              f%bxp(j,k,l)=2.*f%bxp(j,k,l)-f%bxp(j,k-1,l)
           enddo
        enddo
     enddo

     do l=-f%nzguard,f%nz+f%nzguard
        do k=-f%nyguard,f%ny+f%nyguard
           do j=-f%nxguard+1,f%nx+f%nxguard-1
              f%exp(j,k,l)=2.*f%exp(j,k,l)-f%exp(j-1,k,l)
              f%byp(j,k,l)=2.*f%byp(j,k,l)-f%byp(j-1,k,l)
              f%bzp(j,k,l)=2.*f%bzp(j,k,l)-f%bzp(j-1,k,l)
           enddo
        enddo
     enddo

  else

     if (f%l_1dz) then
        j=0
        k=0
        do l=-f%nzguard+1,f%nz+f%nzguard-1
           f%ezp(j,k,l)=2.*f%ezp(j,k,l)-f%ezp(j,k,l-1)
           f%bxp(j,k,l)=2.*f%bxp(j,k,l)-f%bxp(j,k,l-1)
           f%byp(j,k,l)=2.*f%byp(j,k,l)-f%byp(j,k,l-1)
        enddo

     else
        k=0
        do l=-f%nzguard+1,f%nz+f%nzguard-1
           do j=-f%nxguard,f%nx+f%nxguard
              f%ezp(j,k,l)=2.*f%ezp(j,k,l)-f%ezp(j,k,l-1)
              f%bxp(j,k,l)=2.*f%bxp(j,k,l)-f%bxp(j,k,l-1)
              f%byp(j,k,l)=2.*f%byp(j,k,l)-f%byp(j,k,l-1)
           enddo
        enddo

        do l=-f%nzguard,f%nz+f%nzguard
           do j=-f%nxguard+1,f%nx+f%nxguard-1
              f%exp(j,k,l)=2.*f%exp(j,k,l)-f%exp(j-1,k,l)
              f%byp(j,k,l)=2.*f%byp(j,k,l)-f%byp(j-1,k,l)
              f%bzp(j,k,l)=2.*f%bzp(j,k,l)-f%bzp(j-1,k,l)
           enddo
        enddo

        if (f%circ_m>0) then
           do l=-f%nzguard+1,f%nz+f%nzguard-1
              do j=-f%nxguard,f%nx+f%nxguard
                 f%ezp_circ(j,l,:)=2.*f%ezp_circ(j,l,:)-f%ezp_circ(j,l-1,:)
                 f%bxp_circ(j,l,:)=2.*f%bxp_circ(j,l,:)-f%bxp_circ(j,l-1,:)
                 f%byp_circ(j,l,:)=2.*f%byp_circ(j,l,:)-f%byp_circ(j,l-1,:)
              enddo
           enddo

           do l=-f%nzguard,f%nz+f%nzguard
              do j=-f%nxguard+1,f%nx+f%nxguard-1
                 f%exp_circ(j,l,:)=2.*f%exp_circ(j,l,:)-f%exp_circ(j-1,l,:)
                 f%byp_circ(j,l,:)=2.*f%byp_circ(j,l,:)-f%byp_circ(j-1,l,:)
                 f%bzp_circ(j,l,:)=2.*f%bzp_circ(j,l,:)-f%bzp_circ(j-1,l,:)
              enddo
           enddo
        END IF

     endif

  endif

  f%l_nodecentered = .false.

  return
end subroutine node2yee3d

subroutine Jyee2node3d(f)
  ! puts EM value from Yee grid to nodes
  use mod_emfield3d
  implicit none
  TYPE(EM3D_YEEFIELDtype) :: f

  INTEGER :: j,k,l

  if (.not.f%l_2dxz) then
     do l=-f%nzguard,f%nz+f%nzguard
        do k=-f%nyguard,f%ny+f%nyguard
           do j=f%nx+f%nxguard-1,-f%nxguard+1,-1
              f%J(j,k,l,1)=0.5*(f%J(j,k,l,1)+f%J(j-1,k,l,1))
           enddo
        enddo
     enddo

     do l=-f%nzguard,f%nz+f%nzguard
        do k=f%ny+f%nyguard-1,-f%nyguard+1,-1
           do j=-f%nxguard,f%nx+f%nxguard
              f%J(j,k,l,2)=0.5*(f%J(j,k,l,2)+f%J(j,k-1,l,2))
           enddo
        enddo
     enddo

     do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
        do k=-f%nyguard,f%ny+f%nyguard
           do j=-f%nxguard,f%nx+f%nxguard
              f%J(j,k,l,3)=0.5*(f%J(j,k,l,3)+f%J(j,k,l-1,3))
           enddo
        enddo
     enddo

  else
     if (f%l_1dz) then
        j = 0
        k = 0
        do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
           f%J(j,k,l,3)=0.5*(f%J(j,k,l,3)+f%J(j,k,l-1,3))
        enddo
     else

        k = 0
        do l=-f%nzguard,f%nz+f%nzguard
           do j=f%nx+f%nxguard-1,-f%nxguard+1,-1
              f%J(j,k,l,1)=0.5*(f%J(j,k,l,1)+f%J(j-1,k,l,1))
           enddo
        enddo

        do l=f%nz+f%nzguard-1,-f%nzguard+1,-1
           do j=-f%nxguard,f%nx+f%nxguard
              f%J(j,k,l,3)=0.5*(f%J(j,k,l,3)+f%J(j,k,l-1,3))
           enddo
        enddo
     endif
  endif

  return
end subroutine Jyee2node3d

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

  if (f%l_2dxz) then
     do i=1,n
        f%incond(indx(1,i),0,indx(3,i)) = .true.
     end do
  else
     do i=1,n
        f%incond(indx(1,i),indx(2,i),indx(3,i)) = .true.
     end do
  end if

end subroutine set_incond

subroutine set_macroscopic_coefs_on_yee(f,n,indx,sigma,epsi,mu)
use mod_emfield3d
TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: i,n,indx(3,n)
real(kind=8), intent(in) :: sigma,epsi,mu
integer(ISZ) :: j,k,l

  f%incond(:,:,:) = .false.

  if (f%l_2dxz) then
    do i=1,n
      f%incond(indx(1,i),0,indx(3,i)) = .true.
    end do  
  else
    do i=1,n
      f%incond(indx(1,i),indx(2,i),indx(3,i)) = .true.
    end do  
  end if

  ! --- NOTE: if l_2drz is TRUE, then l_2dxz is TRUE
  if (.not. f%l_2dxz) then ! --- 3D XYZ
  ! advance Ex
  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx-1
      if (f%incond(j,k,l) .and. f%incond(j+1,k,l)) then
        f%sigmax(j,k,l) = sigma
        f%epsix(j,k,l) = epsi
        f%mux(j,k,l) = mu
      end if
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      if ( f%incond(j,k,l) .and. f%incond(j,k+1,l)) then
        f%sigmay(j,k,l) = sigma
        f%epsiy(j,k,l) = epsi
        f%muy(j,k,l) = mu
      end if
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      if ( f%incond(j,k,l) .and. f%incond(j,k,l+1)) then
        f%sigmaz(j,k,l) = sigma
        f%epsiz(j,k,l) = epsi
        f%muz(j,k,l) = mu
      end if
    end do
   end do
  end do

else ! --- now 2D XZ or RZ

  k = 0
  ! advance Ex
  do l = 0, f%nz
    do j = 0, f%nx-1
      if ( f%incond(j,k,l) .and. f%incond(j+1,k,l)) then
        f%sigmax(j,k,l) = sigma
        f%epsix(j,k,l) = epsi
        f%mux(j,k,l) = mu
      end if
    end do
  end do

  ! advance Ey
  do l = 0, f%nz
    do j = 0, f%nx
      if (f%incond(j,k,l)) then
        f%sigmay(j,k,l) = sigma
        f%epsiy(j,k,l) = epsi
        f%muy(j,k,l) = mu
      end if
    end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
    do j = 0, f%nx
      if ( f%incond(j,k,l) .and. f%incond(j,k,l+1)) then
        f%sigmaz(j,k,l) = sigma
        f%epsiz(j,k,l) = epsi
        f%muz(j,k,l) = mu
      end if
    end do
  end do

end if


return
end subroutine set_macroscopic_coefs_on_yee

