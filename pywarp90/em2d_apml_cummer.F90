#include "top.h"
!     Last change:  JLV  17 Oct 2003    4:33 pm


!*************  MODULE bnd  *****************
module mod_bnd_cummer
!use mod_EM2D_precision
!use mod_system
use EM2D_bnd
#ifdef MPIPARALLEL
use parallel
#endif
implicit none

INTEGER(ISZ), parameter :: pml = 1, &
                      pml_sadjusted = 2, &
                      apml_exponential = 3, &
                      apml_hybrid = 4, &
                      apml_ssa = 5, &
                      apml_lwa = 6

!INTEGER :: bnd_cond = pml

!REAL(kind=8) :: s_max_init = 4., &
!                               s_max_x, s_max_y, &
!                               s_delta = 5., &
!                               sb_coef = 0., &
!                               nn = 2.

!TYPE :: type_bnd
!  REAL(kind=8), POINTER, DIMENSION(:) :: Ex, Ey
!  REAL(kind=8), POINTER, DIMENSION(:) :: Bzx, Bzy
!  REAL(kind=8), POINTER, DIMENSION(:) :: aEx, bEx, cEx
!  REAL(kind=8), POINTER, DIMENSION(:) :: aEy, bEy, cEy
!  REAL(kind=8), POINTER, DIMENSION(:) :: aBzx, bBzx, cBzx
!  REAL(kind=8), POINTER, DIMENSION(:) :: aBzy, bBzy, cBzy
!  INTEGER :: nx, ny, nbndx, nbndy
!  integer :: n1x,nbot,nint,ntop,nbot1,nbot2,ntop1,ntop2
!END TYPE type_bnd

contains

!************* SUBROUTINE create_bnd  *******************************************************

subroutine create_bnd(b, nx, ny, nbndx, nbndy, dt, dx, dy, xbnd, ybnd)
  implicit none

  INTEGER(ISZ), INTENT(IN) :: nx, ny
  INTEGER(ISZ), INTENT(IN) :: nbndx, nbndy
  REAL(kind=8), INTENT(IN) :: dt, dx, dy
  integer(ISZ) :: xbnd, ybnd

  INTEGER :: n
  TYPE(type_bnd_cummer) :: b

  b%nbndx = nbndx
  b%nbndy = nbndy
!  if(l_moving_window .and. .not. l_elaser_out_plane) b%nbndx=0
!  if(l_moving_window .and.       l_elaser_out_plane) b%nbndy=0
  b%nx    = nx+2*b%nbndx
  b%ny    = ny+2*b%nbndy

 n = 2*(b%nx+1)*(b%nbndy+2)+2*(b%ny-2*b%nbndy-1)*(b%nbndx+2)

 b%n1x=b%nx+1
 b%nbot=b%n1x*(b%nbndy+1)
 b%nint=(b%nbndx+2)+(b%nbndx+1)
 b%ntop=b%nbot+(ny-2)*b%nint

 b%ntop1 = -b%n1x
 b%ntop2 = b%ntop - (b%ny-b%nbndy)*b%n1x
 b%nbot1 = b%nbot - (b%nbndy+2)*b%nint
 b%nbot2 = b%nbot - (b%nbndy+2)*b%nint - nx + 2

  b%n=n
  call type_bnd_cummerallot(b)

  b%Ex   = 0.
  b%Ey   = 0.
  b%Bz   = 0.
  b%Extild   = 0.
  b%Eytild   = 0.
  b%Bzxtild   = 0.
  b%Bzytild   = 0.
  b%aEx  = 0.
  b%aEy  = 0.
  b%aBz  = 0.
  b%bEx  = 0.
  b%bEy  = 0.
  b%bBzx  = 0.
  b%bBzy  = 0.
  b%cEx  = 0.
  b%cEy  = 0.
  b%cBzx  = 0.
  b%cBzy  = 0.
  b%aExtild = 0.
  b%aEytild = 0.
  b%aBzxtild = 0.
  b%aBzytild = 0.

call init_bnd(bnd=b, dt=dt, dx=dx, dy=dy, xbnd=xbnd, ybnd=ybnd)

return
end subroutine create_bnd


!************* SUBROUTINE init_bnd  ********

subroutine init_bnd(bnd, dt, dx, dy, xbnd, ybnd)
  use GlobalVars
  implicit none

  TYPE(type_bnd_cummer), INTENT(INOUT) :: bnd
  REAL(kind=8), INTENT(IN) :: dt, dx, dy
  integer(ISZ) :: xbnd, ybnd
  
  INTEGER(ISZ) :: i, j, k, jmin, kmin

  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: sigmax, sigmax_next, sigmay, sigmay_next

  bnd%aEx  =  1.
  bnd%aEy  =  1.
  bnd%aBz  =  1.

  bnd%bEx  =  dt/dy
  bnd%bEy  =  dt/dx
  bnd%bBzx  =  dt/dx
  bnd%bBzy  =  dt/dy

  bnd%cEx  = -dt/dy
  bnd%cEy  = -dt/dx
  bnd%cBzx  = -dt/dx
  bnd%cBzy  = -dt/dy

  bnd%aExtild  = 1.
  bnd%aEytild  = 1.
  bnd%aBzxtild = 1.
  bnd%aBzytild = 1.
  bnd%bExtild  = 1.
  bnd%bEytild  = 1.
  bnd%bBzxtild = 1.
  bnd%bBzytild = 1.

  s_max_x = s_max_init/dx
  s_max_y = s_max_init/dy

  WRITE(0,*) ' Boundary condition on fields:'
  select case (bnd_cond)
    case (pml)
      WRITE(0,*) '    pml'
    case (pml_sadjusted)
      WRITE(0,*) '    pml sigma-adjusted'
    case (apml_exponential)
      WRITE(0,*) '    apml-exponential'
    case (apml_hybrid)
      WRITE(0,*) '    apml-hybrid'
    case (apml_ssa)
      WRITE(0,*) '    apml-ssa'
    case (apml_lwa)
      WRITE(0,*) '    apml-lwa'
    case default
      WRITE(0,*) '    perfect metal'
  end select
  WRITE(0,'("      nn      = ",f5.2)') nn
  WRITE(0,'("      s_max   = ",f5.2)') s_max_x,s_max_y
  WRITE(0,'("      s_delta = ",f5.2)') s_delta
  WRITE(0,'("      sb_coef = ",f5.2)') sb_coef
  
  ALLOCATE(sigmax(1:bnd%nbndx+1),sigmax_next(1:bnd%nbndx+1),sigmay(1:bnd%nbndy+1),sigmay_next(1:bnd%nbndy+1))

       sigmax = 0.
       sigmax_next = 0.
       sigmay = 0.
       sigmay_next = 0.

       jmin = 2
       kmin = 2
       do j = jmin, bnd%nbndx+1
         sigmax(j)      = s_max_x*(REAL(j-jmin)/s_delta)**nn
         sigmax_next(j) = s_max_x*((REAL(j-jmin)+0.5)/s_delta)**nn
       end do
       do k = kmin, bnd%nbndy+1
         sigmay(k)      = s_max_y*(REAL(k-kmin)/s_delta)**nn
         sigmay_next(k) = s_max_y*((REAL(k-kmin)+0.5)/s_delta)**nn
       end do

     ! Bzy
     if(ybnd==absorb) then
       do k = 1, bnd%nbndy
         do j = 1, bnd%nx
           call assign_coefs(bnd_cond, &
                             bnd%aBzytild(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             bnd%bBzytild(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             bnd%bBzy(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             bnd%cBzy(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             dt,dy, &
                             sigmay_next(k), &
                             sigmay(k+1), &
                             sb_coef)
           bnd%aBzytild(ijk(bnd,j,bnd%nbndy+1-k)) =   bnd%aBzytild(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
           bnd%bBzytild(ijk(bnd,j,bnd%nbndy+1-k)) =   bnd%bBzytild(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
           bnd%bBzy(ijk(bnd,j,bnd%nbndy+1-k)) = - bnd%cBzy(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
           bnd%cBzy(ijk(bnd,j,bnd%nbndy+1-k)) = - bnd%bBzy(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
         end do
       end do
     end if

     ! Bzx

     if(xbnd==absorb) then
       do j = 1, bnd%nbndx
         do k = 1, bnd%ny
           call assign_coefs(bnd_cond, &
                             bnd%aBzxtild(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             bnd%bBzxtild(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             bnd%bBzx(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             bnd%cBzx(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             dt,dx, &
                             sigmax_next(j), &
                             sigmax(j+1), &
                             sb_coef)
           bnd%aBzxtild(ijk(bnd,bnd%nbndx+1-j,k)) =   bnd%aBzxtild(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
           bnd%bBzxtild(ijk(bnd,bnd%nbndx+1-j,k)) =   bnd%bBzxtild(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
           bnd%bBzx(ijk(bnd,bnd%nbndx+1-j,k)) = - bnd%cBzx(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
           bnd%cBzx(ijk(bnd,bnd%nbndx+1-j,k)) = - bnd%bBzx(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
         end do
       end do
     end if

       ! Exy
     if(ybnd==absorb) then
       do k = 1, bnd%nbndy
         do j = 1, bnd%nx
           call assign_coefs(bnd_cond, &
                             bnd%aExtild(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             bnd%bExtild(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             bnd%bEx(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             bnd%cEx(ijk(bnd,j,k+bnd%ny-bnd%nbndy)), &
                             dt,dy, &
                             sigmay(k), &
                             sigmay_next(k), &
                             sb_coef)
           bnd%aExtild(ijk(bnd,j,bnd%nbndy+2-k)) =  bnd%aExtild(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
           bnd%bExtild(ijk(bnd,j,bnd%nbndy+2-k)) =  bnd%bExtild(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
           bnd%bEx(ijk(bnd,j,bnd%nbndy+2-k)) = -bnd%cEx(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
           bnd%cEx(ijk(bnd,j,bnd%nbndy+2-k)) = -bnd%bEx(ijk(bnd,j,k+bnd%ny-bnd%nbndy))
         end do
       end do
     end if

       ! Eyx

     if(xbnd==absorb) then
       do j = 1, bnd%nbndx
         do k = 1, bnd%ny
           call assign_coefs(bnd_cond, &
                             bnd%aEytild(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             bnd%bEytild(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             bnd%bEy(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             bnd%cEy(ijk(bnd,j+bnd%nx-bnd%nbndx,k)), &
                             dt,dx, &
                             sigmax(j), &
                             sigmax_next(j), &
                             sb_coef)
           bnd%aEytild(ijk(bnd,bnd%nbndx+2-j,k)) =  bnd%aEytild(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
           bnd%bEytild(ijk(bnd,bnd%nbndx+2-j,k)) =  bnd%bEytild(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
           bnd%bEy(ijk(bnd,bnd%nbndx+2-j,k)) = -bnd%cEy(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
           bnd%cEy(ijk(bnd,bnd%nbndx+2-j,k)) = -bnd%bEy(ijk(bnd,j+bnd%nx-bnd%nbndx,k))
         end do
       end do
     end if
  DEALLOCATE(sigmax, sigmax_next, sigmay, sigmay_next)

return
end subroutine init_bnd

subroutine assign_coefs(bnd_cond,atild,btild,bp,bm,dt,dx,sigma,sigma_next,coef_sigmab)
implicit none
REAL(kind=8), INTENT(OUT) :: atild,btild,bp,bm
REAL(kind=8), INTENT(IN) :: dt,dx,sigma,sigma_next,coef_sigmab
INTEGER(ISZ),INTENT(IN) :: bnd_cond

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
        atild = 1.
        btild = 1.
      else
        atild =  EXP(-sigma_local*dt)
        btild = (1.-atild)/(sigma_local*dt)
!        atild =  (1.-sigma_local*dt/2)/(1.+sigma_local*dt/2)
!        btild = dt/(1.+sigma_local*dt/2)
      END if
      bp =  dt / dx
      bm =  -bp
    case (apml_exponential)
      sigmab = coef_sigmab*sigma
      IF(sigma == 0.) then
      ! --- end of mesh
        atild  =  1.
        bp =  dt / dx
        bm =  -bp
        btild = 1.
      else
        atild  =  EXP(-sigma*dt)
        btild = (1.-atild)/(sigma*dt)
        IF(sigmab==0.) then
          bp =  dt / dx
          bm =  -bp
        else
          bp =  dt*sigmab/(1.-EXP(-sigmab*dx))
          bm =  -EXP(-sigmab*dx)*bp
        END if
      END if
    case (apml_hybrid)
      write(0,*) 'apml cummer hybrid is not yet implemented'
      stop
      bp = dt/dx
      bm = -dt/dx*(1.+((dx-dt)/(dx+dt))*(1.-tpp))
      atild  = 1.+bp*tpp+bm
      bm = tp*bm
    case (apml_ssa,apml_lwa)
      write(0,*) 'apml cummer ssa/lwa are not yet implemented'
      stop
      sigmab      = coef_sigmab*sigma
      sigmab_next = coef_sigmab*sigma_next
      tp  = EXP(-(sigma+sigmab)*0.5*dx)
      tmm = EXP(-(sigma-sigmab)*0.5*dx)
      tpp = EXP(-(sigma_next+sigmab_next)*0.5*dx)
      tm  = EXP(-(sigma_next-sigmab_next)*0.5*dx)
      IF(bnd_cond==apml_ssa) then
        g = dx/dt
        atild  = -(1-g*(tm+tp)-g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)/(1+g*(tm+tp)+g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)
        bp = 2*tm*(1.+tp*tmm)/(1.+g*(tm+tp)+g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)
        bm = -2*tp*(1.+tm*tpp)/(1.+g*(tm+tp)+g*tm*tp*(tmm+tpp)-tm*tp*tmm*tpp)
      else
        gp = dx/dt
        gm = gp
        atild = -(1-gm-(gp+gm)*tm*tpp-(gp+1)*tm*tmm*tpp*tp)/(1+gm+(gp+gm)*tm*tpp+(gp-1)*tm*tmm*tpp*tp)
        bp = 2*tm*(1.+tp*tmm)/(1+gm+(gp+gm)*tm*tpp+(gp-1)*tm*tmm*tpp*tp)
        bm = -2*tp*(1.+tm*tpp)/(1+gm+(gp+gm)*tm*tpp+(gp-1)*tm*tmm*tpp*tp)
      END if
    case default
  end select
END subroutine assign_coefs

subroutine move_window_bnd(bnd,rap,l_elaser_out_plane)
  TYPE(type_bnd_cummer), INTENT(INOUT) :: bnd
  integer :: rap,j,k,ijk1,ijk2,ijk3,ijk4
  logical :: l_elaser_out_plane
  
       do k = 1, bnd%ny-rap
         do j = 1, bnd%nbndx
           if (l_elaser_out_plane) then
             ijk1 = ijk(bnd,bnd%nbndx+2-j,k)
             ijk2 = ijk(bnd,bnd%nbndx+2-j,k+rap)
             ijk3 = ijk(bnd,j+bnd%nx-bnd%nbndx,k)
             ijk4 = ijk(bnd,j+bnd%nx-bnd%nbndx,k+rap)
           else
             ijk1 = ijk(bnd,j,    bnd%nbndy+2-k)
             ijk2 = ijk(bnd,j+rap,bnd%nbndy+2-k)
             ijk3 = ijk(bnd,j,    k+bnd%ny-bnd%nbndy)
             ijk4 = ijk(bnd,j+rap,k+bnd%ny-bnd%nbndy)
           end if
           bnd%Ex(ijk1) = bnd%Ex(ijk2) 
           bnd%Ex(ijk3) = bnd%Ex(ijk4)
           bnd%Ey(ijk1) = bnd%Ey(ijk2) 
           bnd%Ey(ijk3) = bnd%Ey(ijk4)
           bnd%Bz(ijk1) = bnd%Bz(ijk2) 
           bnd%Bz(ijk3) = bnd%Bz(ijk4)
           bnd%Extild(ijk1) = bnd%Extild(ijk2) 
           bnd%Extild(ijk3) = bnd%Extild(ijk4)
           bnd%Eytild(ijk1) = bnd%Eytild(ijk2) 
           bnd%Eytild(ijk3) = bnd%Eytild(ijk4)
           bnd%Bzxtild(ijk1) = bnd%Bzxtild(ijk2) 
           bnd%Bzxtild(ijk3) = bnd%Bzxtild(ijk4)
           bnd%Bzytild(ijk1) = bnd%Bzytild(ijk2) 
           bnd%Bzytild(ijk3) = bnd%Bzytild(ijk4)
         end do
       end do
       do k = bnd%ny-rap+1,bnd%ny+1
         do j = 1, bnd%nbndx
           if (l_elaser_out_plane) then
             ijk1 = ijk(bnd,bnd%nbndx+2-j,k)
             ijk2 = ijk(bnd,j+bnd%nx-bnd%nbndx,k)
           else
             ijk1 = ijk(bnd,j,bnd%nbndy+2-k)
             ijk2 = ijk(bnd,j,k+bnd%ny-bnd%nbndy)
           end if
           bnd%Ex(ijk1) = 0.
           bnd%Ex(ijk2) = 0.
           bnd%Ey(ijk1) = 0.
           bnd%Ey(ijk2) = 0.
           bnd%Bz(ijk1) = 0.
           bnd%Bz(ijk2) = 0.
           bnd%Extild(ijk1) = 0.
           bnd%Extild(ijk2) = 0.
           bnd%Eytild(ijk1) = 0.
           bnd%Eytild(ijk2) = 0.
           bnd%Bzxtild(ijk1) = 0.
           bnd%Bzxtild(ijk2) = 0.
           bnd%Bzytild(ijk1) = 0.
           bnd%Bzytild(ijk2) = 0.
         end do
       end do

  return
end subroutine move_window_bnd

subroutine move_bnd(b)
  implicit none

  INTEGER :: j, k, jf, kf, jb, kb,jk,jk1,i,ntop,kmin,kmax
  TYPE(type_bnd_cummer), POINTER :: b
  real(kind=8) :: oldval
  real(kind=8):: bzrecv(b%nbndy,3),bztosend(b%nbndy,3)

#ifdef MPIPARALLEL
  include "mpif.h"
  integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpierror
  integer(MPIISZ):: mpirequest
  integer(MPIISZ):: w
  integer(MPIISZ):: messid 
#endif

  ! Bz
    do k = 1, b%nbndy
    jk1 = b%ntop1 + k * b%n1x

      do j = 1, b%nx
      jk = jk1 + j
        
        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(jk+b%n1x)  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(jk+1)		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)
        
      end do
    end do

    do k =   b%ny-b%nbndy+1, b%ny
    jk1 = b%ntop2 + k * b%n1x

      do j = 1, b%nx
      jk = jk1 + j

        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(jk+b%n1x)  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(jk+1)		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)
      end do
    end do

    do k = b%nbndy+2, b%ny-b%nbndy-2
      do j = 1, b%nbndx
        jk= b%nbot1 + j + k * b%nint

        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(jk+b%nint)  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(jk+1)		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)

      end do
    end do

    k=b%nbndy+1
    do j = 1, b%nbndx
    jk=ijk(b,j,k)

        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(ijk(b,j,k+1))  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(ijk(b,j+1,k))		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)

   enddo

    do k=b%ny-b%nbndy-1,b%ny-b%nbndy
    do j = 1, b%nbndx
    jk=ijk(b,j,k)

        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(ijk(b,j,k+1))  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(ijk(b,j+1,k))		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)

   enddo
   enddo


    do k = b%nbndy+2, b%ny-b%nbndy-2
     do j = b%nx-b%nbndx+1, b%nx
      jk = b%nbot2+ j + k * b%nint

        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(jk+b%nint)  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(jk+1)		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)
       end do
    end do

    k=b%nbndy+1
    do j = b%nx-b%nbndx+1, b%nx
    jk=ijk(b,j,k)
        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(ijk(b,j,k+1))  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(ijk(b,j+1,k))		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)
    enddo

    do k=b%ny-b%nbndy-1,b%ny-b%nbndy
    do j =  b%nx-b%nbndx+1, b%nx
    jk=ijk(b,j,k)
        oldval = b%Bz(jk)
        b%Bz(jk) = b%aBz(jk)  * b%Bz(jk) 	        &
                 + b%bBzy(jk) * b%Extild(ijk(b,j,k+1))  &
                 + b%cBzy(jk) * b%Extild(jk)        &
                 - b%bBzx(jk) * b%Eytild(ijk(b,j+1,k))		&
                 - b%cBzx(jk) * b%Eytild(jk)

        b%Bzxtild(jk) = b%aBzxtild(jk) * b%Bzxtild(jk) &
                      + b%bBzxtild(jk) * (b%Bz(jk)-oldval)
                      
        b%Bzytild(jk) = b%aBzytild(jk) * b%Bzytild(jk) &
                      + b%bBzytild(jk) * (b%Bz(jk)-oldval)
   enddo
   enddo
   
#ifdef MPIPARALLEL
  if(my_index>0) then
   do i = 1,2
    messid=100
!    write(0,*) my_index,' sends data to ',my_index-1
    j = b%nbndx+1
    if (i==1) then
      kmin = 1
      kmax = b%nbndy
      ntop = b%ntop1
    else
      kmin = b%ny-b%nbndy+1
      kmax = b%ny
      ntop = b%ntop2
    end if
    do k = kmin, kmax
      jk1 = ntop + k * b%n1x
      jk = jk1 + j
      Bztosend(k-kmin+1,1) = b%Bzxtild(jk)
      Bztosend(k-kmin+1,2) = b%Bzytild(jk)
      Bztosend(k-kmin+1,3) = b%Bz(jk)
    end do
    call MPI_ISEND(Bztosend,3*b%nbndy,MPI_DOUBLE_PRECISION, &
                   my_index-1,messid,MPI_COMM_WORLD,mpirequest,mpierror)
!    write(0,*) 'done'
    messid=101
!    write(0,*) my_index,' recv data from ',my_index-1
    call MPI_RECV(Bzrecv,3*b%nbndy,MPI_DOUBLE_PRECISION, &
                  my_index-1,messid,MPI_COMM_WORLD,mpistatus,mpierror)
    j = b%nbndx+0
    do k = kmin, kmax
      jk1 = ntop + k * b%n1x
      jk = jk1 + j
      b%Bzxtild(jk) = Bzrecv(k-kmin+1,1)
      b%Bzytild(jk) = Bzrecv(k-kmin+1,2)
      b%Bz(jk) = Bzrecv(k-kmin+1,3)
    end do
!    write(0,*) 'done'
   end do
  end if

  if(my_index<nslaves-1) then
   do i = 1,2
    messid=101
!    write(0,*) my_index,' sends data to ',my_index+1
    j = b%nx-b%nbndx+0
    if (i==1) then
      kmin = 1
      kmax = b%nbndy
      ntop = b%ntop1
    else
      kmin = b%ny-b%nbndy+1
      kmax = b%ny
      ntop = b%ntop2
    end if
    do k = kmin, kmax
!      jk1 = b%ntop1 + k * b%n1x
!      jk = jk1 + j
    jk=ijk(b,j,k)
      Bztosend(k-kmin+1,1) = b%Bzxtild(jk)
      Bztosend(k-kmin+1,2) = b%Bzytild(jk)
      Bztosend(k-kmin+1,3) = b%Bz(jk)
    end do
    call MPI_ISEND(Bztosend,3*b%nbndy,MPI_DOUBLE_PRECISION, &
                   my_index+1,messid,MPI_COMM_WORLD,mpirequest,mpierror)
!    write(0,*) 'done'
    messid=100
!    write(0,*) my_index,' recv data from ',my_index+1
    call MPI_RECV(Bzrecv,3*b%nbndy,MPI_DOUBLE_PRECISION, &
                  my_index+1,messid,MPI_COMM_WORLD,mpistatus,mpierror)
    j = b%nx-b%nbndx+1
    do k = kmin, kmax
!      jk1 = b%ntop1 + k * b%n1x
!      jk = jk1 + j
    jk=ijk(b,j,k)
      b%Bzxtild(jk) = Bzrecv(k-kmin+1,1)
      b%Bzytild(jk) = Bzrecv(k-kmin+1,2)
      b%Bz(jk) = Bzrecv(k-kmin+1,3)
    end do
!    write(0,*) 'done'
   end do
  end if
#endif

  ! Ex
    do k = 2, b%nbndy
    jk1 = b%ntop1 + k*b%n1x

      do j = 1, b%nx
      jk = jk1 + j

        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(jk-b%n1x)

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)

      end do
    end do

    do k = b%ny-b%nbndy+2, b%ny
    jk1 = b%ntop2 + k*b%n1x


      do j = 1, b%nx
      jk = jk1 + j

        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(jk-b%n1x)

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)
       end do
    end do

   do k = b%nbndy+3, b%ny-b%nbndy-1
    do j = 1, b%nbndx
      jk  = b%nbot1 + j + k*b%nint

        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(jk-b%nint)

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)
      end do
    end do

  do k=b%nbndy+1,b%nbndy+2
  do j = 1, b%nbndx
  jk=ijk(b,j,k)
        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(ijk(b,j,k-1))

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)
      end do
    end do

   do k=b%ny-b%nbndy,b%ny-b%nbndy+1
   do j = 1, b%nbndx
   jk=ijk(b,j,k)

        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(ijk(b,j,k-1))

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)
      end do
    end do

    do k = b%nbndy+3, b%ny-b%nbndy-1
     do j = b%nx-b%nbndx+1, b%nx
      jk  = b%nbot2 + j + k*b%nint

        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(jk-b%nint)

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)
        end do
    end do

  do k=b%nbndy+1,b%nbndy+2
  do j =b%nx-b%nbndx+1, b%nx
  jk=ijk(b,j,k)
        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(ijk(b,j,k-1))

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)
      end do
    end do

   do k=b%ny-b%nbndy,b%ny-b%nbndy+1
   do j = b%nx-b%nbndx+1, b%nx
   jk=ijk(b,j,k)
        oldval = b%Ex(jk)
        
        b%Ex(jk) = b%aEx(jk) * b%Ex(jk)  &
                 + b%bEx(jk) * b%Bzytild(jk) 	&
                 + b%cEx(jk) * b%Bzytild(ijk(b,j,k-1))

        b%Extild(jk) = b%aExtild(jk) * b%Extild(jk) &
                     + b%bExtild(jk) * (b%Ex(jk)-oldval)

      end do
    end do


  ! Ey

    do k = 1, b%nbndy
     do j = 2, b%nbndx
      jk=b%ntop1+j+k*b%n1x

        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)

      end do
    end do

    do k=b%nbndy+2,b%ny-b%nbndy-2
    do j = 2, b%nbndx
    jk=b%nbot1+j+k*b%nint
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)

      end do
     end do
     
     do k=b%ny-b%nbndy,b%ny
     do j = 2, b%nbndx
     jk=b%ntop2+j+k*b%n1x
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)

      end do
      end do

       k=b%nbndy+1
       do j = 2, b%nbndx
       jk=ijk(b,j,k)

        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(ijk(b,j-1,k))

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do

     k=b%ny-b%nbndy-1
      do j = 2, b%nbndx
      jk=ijk(b,j,k)

        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(ijk(b,j-1,k))

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do



    do k = 1, b%nbndy
     do j = b%nx-b%nbndx+2, b%nx
      jk=b%ntop1+j+k*b%n1x

        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do
     end do

    do k=b%nbndy+2,b%ny-b%nbndy-2
     do j = b%nx-b%nbndx+2, b%nx
      jk=b%nbot2+j+k*b%nint

        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do
     end do
     
     do k=b%ny-b%nbndy,b%ny
      do j = b%nx-b%nbndx+2, b%nx
       jk=b%ntop2+j+k*b%n1x

        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do
     end do

     k=b%nbndy+1
     do j = b%nx-b%nbndx+2, b%nx
       jk=ijk(b,j,k)
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(ijk(b,j-1,k))

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do

      k=b%ny-b%nbndy-1
      do j = b%nx-b%nbndx+2, b%nx
        jk=ijk(b,j,k)
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(ijk(b,j-1,k))

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
      end do
 
    if(.false.) then !l_moving_window) then
      do k = 1, b%nbndy
        jk1=b%ntop1+k*b%n1x

        do j = b%nbndx+2, b%nx-b%nbndx+1
          jk=jk1+j
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
        end do
      end do

      do k = b%ny-b%nbndy+1, b%ny
        jk1=b%ntop2+k*b%n1x

        do j = b%nbndx+2, b%nx-b%nbndx+1
          jk=jk1+j
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
        end do
      end do
    else
      do k = 1, b%nbndy
        jk1=b%ntop1+k*b%n1x

        do j = b%nbndx+1, b%nx-b%nbndx+1
          jk=jk1+j
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
        end do
      end do

      do k = b%ny-b%nbndy+1, b%ny
        jk1=b%ntop2+k*b%n1x

        do j = b%nbndx+1, b%nx-b%nbndx+1
          jk=jk1+j
        oldval = b%Ey(jk)
        
        b%Ey(jk) = b%aEy(jk) * b%Ey(jk)  &
                 - b%bEy(jk) * b%Bzxtild(jk) 	&
                 - b%cEy(jk) * b%Bzxtild(jk-1)

        b%Eytild(jk) = b%aEytild(jk) * b%Eytild(jk) &
                     + b%bEytild(jk) * (b%Ey(jk)-oldval)
        end do
      end do
    end if
return
end subroutine move_bnd

!************* FUNCTION ijk  *************

function ijk(bnd,j,k)
implicit none
TYPE(type_bnd_cummer), INTENT(INOUT) :: bnd
INTEGER(ISZ) :: ijk
INTEGER(ISZ), INTENT(IN) :: j, k

IF(k<1.or.k>bnd%ny+1.or.j<1.or.j>bnd%nx+1) then
  WRITE(0,*) 'Error ijk: ',j,k
  call kaboom("ijk: invalid input")
  return
END if

IF(k<=bnd%nbndy+1) then

  ijk = bnd%ntop1 + j  + k * bnd%n1x

ELSEIF(k>=bnd%ny-bnd%nbndy) then

  ijk = bnd%ntop2 +  j + k * bnd%n1x  
        		

ELSEIF(j<=bnd%nbndx+1) then

  ijk = bnd%nbot1 + j + k * bnd%nint 	
        

ELSEIF(j>=bnd%nx-bnd%nbndx) then

  ijk = bnd%nbot2 + j + k * bnd%nint
         
ELSE

  WRITE(0,'("Error in function ijk: j=",i8," k=",i8)') j, k
  call kaboom("ijk: invalid input")
  return

END if

return
end function ijk

end module mod_bnd_cummer

