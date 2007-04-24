#include "top.h"
!     Last change:  JLV   3 Jun 2004    0:17 am
!************* MODULE field  **********************************************

module mod_field

use EM2D_FIELDtypemodule
USE mod_bnd
USE EM2D_FIELDobjects
use GlobalVars
use Picglb
implicit none

TYPE bnd_pointer
  type(type_bnd), POINTER :: b
end type bnd_pointer
!type(bnd_pointer), dimension(3,2) :: bnds ! first dimension is for [main grid, coarse patch, fine patch]
!                                          ! second dimension is for [(Ex,Ey,Bz),(Bx,By,Ez)]
!INTEGER, parameter :: base=1, patchcoarse=2, patchfine=3

contains

subroutine champ_b(f,dt)
implicit none

INTEGER :: j, k
real(kind=8) :: dt
real(kind=8) :: dtsdx,dtsdy
real(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Exapr, Eyapr
TYPE(EM2D_FIELDtype) :: f

dtsdx = dt/f%dx
dtsdy = dt/f%dy

! advance Bx
do k = 1, f%ny+1
  do j = 0, f%nx+1
    f%Bx(j,k) = f%Bx(j,k) - dtsdy * (f%Ez(j,k)   - f%Ez(j,k-1)) 
  end do
end do

! advance By
do k = 0, f%ny+1
  do j = 1, f%nx+1
    f%By(j,k) = f%By(j,k) + dtsdx * (f%Ez(j,k)   - f%Ez(j-1,k)) 
  end do
end do

! advance Bz 
do k = 0, f%ny+1
  do j = 0, f%nx+1
    f%Bz(j,k) = f%Bz(j,k) - dtsdx * (f%Ey(j+1,k) - f%Ey(j,k)) &
                          + dtsdy * (f%Ex(j,k+1) - f%Ex(j,k))
  end do
end do

! add source term
!if (l_elaser_out_plane) then
!  f%cst1 = (1.-dt*f%clight/f%dy)/(1.+dt*f%clight/f%dy)
!  f%cst2 =  2.*dt/f%dy /(1.+dt*f%clight/f%dy)
!else
  f%cst1 = (1.-dt*f%clight/f%dx)/(1.+dt*f%clight/f%dx)
  f%cst2 =  2.*dt/f%dx /(1.+dt*f%clight/f%dx)
!end if
if (f%l_add_source) then
 if (.not.l_moving_window) then
  if(l_elaser_out_plane) then
  ! add contribution from Ez source (Ez_in)
!    k = f%js
!    do j = 0, f%nx+1
!      f%Bx(j,k) = f%Bx(j,k) - dtsdy * (-f%Ez_in(j))
!    end do

    ! Evaluate Bx and By source
!    do j = 1, f%nx+1
!      f%By_in(j) = f%By_in(j) + dtsdx * (f%Ez_in(j) - f%Ez_in(j-1))
!    end do

!    do j = 0, f%nx+1
!      f%Bx_in(j) = f%cst1 * f%Bx_in(j) + f%cst2 * f%Ez_in(j)
!    end do
 
    j = f%js
    do k = 0, f%ny+1
      f%By(j,k) = f%By(j,k) + dtsdx * (-f%Ez_in(k))
    end do

    ! Evaluate Bx and By source
    do k = 1, f%ny+1
      f%Bx_in(k) = f%Bx_in(k) - dtsdy * (f%Ez_in(k) - f%Ez_in(k-1))
    end do

    do k = 0, f%ny+1
      f%By_in(k) = f%cst1 * f%By_in(k) - f%cst2 * f%Ez_in(k)
    end do


  else
    j = f%js-1
    do k = 0, f%ny+1
      f%Bz(j,k) = f%Bz(j,k) - dtsdx * (-f%Ey_in(k)) 
    end do
  end if
 else
  if(.not. l_elaser_out_plane) then
   j = 0
   do k = 0, f%ny+1
     f%Bz(j,k) = f%Bz_in(k) 
   end do
  end if
 end if
end if

!IF(f%l_addpatchresidual) then
!  ALLOCATE(Exapr(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1))
!  ALLOCATE(Eyapr(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1))
!  Exapr = 0.
!  Eyapr = 0.
!  call project_ex(exfin=fpatchfine%ex(0:fpatchfine%nx+1,0:fpatchfine%ny+1), excoarse=Exapr,rap=rap)
!  call project_ey(eyfin=fpatchfine%ey(0:fpatchfine%nx+1,0:fpatchfine%ny+1), eycoarse=Eyapr,rap=rap)
!  Exapr = Exapr - fpatchcoarse%ex(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1)
!  Eyapr = Eyapr - fpatchcoarse%ey(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1)
!  j = ntamp_apr
!  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
!    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) - dtsdx * Eyapr(j,k)
!  END do
!  j = fpatchcoarse%nx+1-ntamp_apr
!  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
!    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) + dtsdx * Eyapr(j+1,k)
!  END do
!  k = ntamp_apr
!  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
!    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) + dtsdy * Exapr(j,k)
!  END do
!  k = fpatchcoarse%ny+1-ntamp_apr
!  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
!    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) - dtsdy * Exapr(j,k+1)
!  END do
!  DEALLOCATE(Exapr,Eyapr)
!END if

end subroutine champ_b


subroutine champ_e(f,dt)
implicit none

TYPE(EM2D_FIELDtype) :: f
INTEGER :: j, k
real(kind=8) :: dt,dtsdx,dtsdy,mudt
real(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Bzapr

dtsdx = f%clight**2*dt/f%dx
dtsdy = f%clight**2*dt/f%dy
mudt  = f%mu0*f%clight**2*dt

if (f%l_apply_pml) then
  call exchange_bnd_field2(f%bndbxbyez,f)
END IF

if (f%l_apply_pml) then 
  call move_bnd(f%bndexeybz)
  call move_bnd(f%bndbxbyez)
end if

! advance Ex
do k = 1, f%ny+1
  do j = 0, f%nx+1
    f%Ex(j,k) = f%Ex(j,k) + dtsdy * (f%Bz(j,k)   - f%Bz(j,k-1)) &
                          - mudt  * f%J(j,k,1)
  end do
end do

! advance Ey
do k = 0, f%ny+1
  do j = 1, f%nx+1
    f%Ey(j,k) = f%Ey(j,k) - dtsdx * (f%Bz(j,k)   - f%Bz(j-1,k)) &
                          - mudt  * f%J(j,k,2)
  end do
end do

! advance Ez 
do k = 0, f%ny+1
  do j = 0, f%nx+1
    f%Ez(j,k) = f%Ez(j,k) + dtsdx * (f%By(j+1,k) - f%By(j,k)) &
                          - dtsdy * (f%Bx(j,k+1) - f%Bx(j,k)) &
                          - mudt  * f%J(j,k,3)
  end do
end do

!if (l_elaser_out_plane) then
!  f%cst1 = (1.-dt*f%clight/f%dy)/(1.+dt*f%clight/f%dy)
!  f%cst2 =  f%clight**2*2.*dt/f%dy /(1.+dt*f%clight/f%dy)
!else
  f%cst1 = (1.-dt*f%clight/f%dx)/(1.+dt*f%clight/f%dx)
  f%cst2 =  f%clight**2*2.*dt/f%dx /(1.+dt*f%clight/f%dx)
!end if
if (f%l_add_source .and. .not.l_moving_window .and. .not. l_elaser_out_plane) then
    ! Evaluate Ex and Ey source
    do k = 1, f%ny+1
      f%Ex_in(k) = f%Ex_in(k) + dtsdy * (f%Bz_in(k) - f%Bz_in(k-1))
    end do

    do k = 0, f%ny+1
      f%Ey_in(k) = f%cst1 * f%Ey_in(k) + f%cst2 * f%Bz_in(k)
    end do
end if

if (f%l_add_source) then
 if (.not.l_moving_window) then
  if(l_elaser_out_plane) then
!    ! substract contribution from Bx source (Bx_in)
!    k = f%js-1
!    do j = 0, f%nx+1
!      f%Ez(j,k) = f%Ez(j,k) - dtsdy * (-f%Bx_in(j)) 
!    end do
    ! substract contribution from Bx source (Bx_in)
    j = f%js-1
    do k = 0, f%ny+1
      f%Ez(j,k) = f%Ez(j,k) + dtsdx * (-f%By_in(k)) 
    end do
  else
    ! add contribution from Bz source (Bz_in)
    j = f%js
    do k = 0, f%ny+1
      f%Ey(j,k) = f%Ey(j,k) - dtsdx * (-f%Bz_in(k))
    end do
  end if
 else
  if(l_elaser_out_plane) then
!    k = 0
!    do j = 0, f%nx+1
!      f%Ez(j,k) = f%Ez_in(j) 
!    end do
    j = 0
    do k = 0, f%ny+1
      f%Ez(j,k) = f%Ez_in(k) 
    end do
  end if
 end if
end if

!IF(f%l_addpatchresidual) then
!  ALLOCATE(Bzapr(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1))
!  Bzapr = 0. 
!  call project_bz(bzfin=fpatchfine%bz(0:fpatchfine%nx+1,0:fpatchfine%ny+1), bzcoarse=Bzapr,rap=rap)
!  Bzapr = Bzapr - fpatchcoarse%bz(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1)
!  
!  j = ntamp_apr
!  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
!    f%Ey(j+ixpatch,k+iypatch) = f%Ey(j+ixpatch,k+iypatch) - dtsdx * Bzapr(j,k)
!  END do
!  j = fpatchcoarse%nx+2-ntamp_apr
!  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
!    f%Ey(j+ixpatch,k+iypatch) = f%Ey(j+ixpatch,k+iypatch) + dtsdx * Bzapr(j-1,k)
!  END do
!  k = ntamp_apr
!  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
!    f%Ex(j+ixpatch,k+iypatch) = f%Ex(j+ixpatch,k+iypatch) + dtsdy * Bzapr(j,k)
!  END do
!  k = fpatchcoarse%ny+2-ntamp_apr
!  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
!    f%Ex(j+ixpatch,k+iypatch) = f%Ex(j+ixpatch,k+iypatch) - dtsdy * Bzapr(j,k-1)
!  END do
!  DEALLOCATE(Bzapr)
!END if

if (f%l_apply_pml) then
  call exchange_bnd_field(f%bndexeybz,f)
END IF

return
end subroutine champ_e

!************* SUBROUTINE exchange_bnd_field  ****************************************

subroutine exchange_bnd_field(b, f)
 implicit none

TYPE(type_bnd  ) :: b
TYPE(EM2D_FIELDtype) :: f

INTEGER :: jb, kb, jf, kf,jk1,jk

if(f%ylbound==periodic) then
  f%Ex(0:f%nx+1,0)      = f%Ex(0:f%nx+1,f%ny+1)
  f%Ex(0:f%nx+1,f%ny+2) = f%Ex(0:f%nx+1,1)

!else if(.not. (l_moving_window .and. l_elaser_out_plane)) then
else 
  kf = 0
  kb = kf + b%nbndy
  jk1=b%ntop1+kb*b%n1x
  do jf = 0, f%nx+1
    jb = jf+b%nbndx
    jk=jk1+jb
    f%Ex(jf,kf) = b%Ex(jk)
  end do

  kf = f%ny+2
  kb = kf + b%nbndy
  jk1=b%ntop2+kb*b%n1x
  do jf = 0, f%nx+1
    jb = jf+b%nbndx
    jk=jk1+jb
    f%Ex(jf,kf) = b%Ex(jk)
  end do

  kf = 1
  kb = kf + b%nbndy
  jk1=b%ntop1+kb*b%n1x
  do jf = 1, f%nx
    jb = jf+b%nbndx
    jk=jk1+jb
    b%Ex(jk) = f%Ex(jf,kf)
  end do

  kf = f%ny+1
  kb = kf + b%nbndy
  jk1=b%ntop2+kb*b%n1x
  do jf = 1, f%nx
    jb = jf+b%nbndx
    jk=jk1+jb
    b%Ex(jk) = f%Ex(jf,kf)
  end do

  kf = 0
  kb = kf + b%nbndy
  jk1=b%ntop1+kb*b%n1x
  do jf = 1, f%nx+1
    jb = jf+b%nbndx
    jk = jk1+jb
    b%Ey(jk) = f%Ey(jf,kf)
  end do

  kf = f%ny+1
  kb = kf + b%nbndy
  jk1=b%ntop2+kb*b%n1x
  do jf = 1, f%nx+1
    jb = jf+b%nbndx
    jk = jk1+jb
    b%Ey(jk) = f%Ey(jf,kf)
  end do
end if

if(f%xlbound==periodic) then
  f%Ey(0,0:f%ny+1)      = f%Ey(f%nx+1,0:f%ny+1)
  f%Ey(f%nx+2,0:f%ny+1) = f%Ey(1,0:f%ny+1)

else if(.not. (l_moving_window .and. (.not.l_elaser_out_plane))) then
  jf = 0
  jb = jf + b%nbndx
  jk1=b%ntop1+jb
  do kf=0,1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%Ey(jf,kf)=b%Ey(jk)
  enddo
  jk1=b%nbot1+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    f%Ey(jf,kf)=b%Ey(jk)
  end do
  jk1=b%ntop2+jb
  do kf=f%ny,f%ny+1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%Ey(jf,kf)=b%Ey(jk)
  enddo


  jf = f%nx+2
  jb = jf + b%nbndx
  jk1=b%ntop1+jb
  do kf=0,1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%Ey(jf,kf)=b%Ey(jk)
  enddo
  jk1=b%nbot2+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    f%Ey(jf,kf)=b%Ey(jk)
  end do
  jk1=b%ntop2+jb
  do kf=f%ny,f%ny+1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%Ey(jf,kf)=b%Ey(jk)
  enddo

  jf = 0
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
  if(f%js==1) then
    b%Ex(jk) = f%Ex(jf,kf) - f%Ex_in(kf)
  else
    b%Ex(jk) = f%Ex(jf,kf) 
  end if

  jk1=b%nbot1+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    if(f%js==1) then
      b%Ex(jk) = f%Ex(jf,kf) - f%Ex_in(kf)
    else
      b%Ex(jk) = f%Ex(jf,kf) 
    end if
  end do
  jk1=b%ntop2+jb
  do kf = f%ny, f%ny+1
    kb = kf+b%nbndy
    jk=jk1+kb*b%n1x
    if(f%js==1) then
      b%Ex(jk) = f%Ex(jf,kf) - f%Ex_in(kf)
    else
      b%Ex(jk) = f%Ex(jf,kf) 
    end if
  end do

  jf = f%nx+1
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
  b%Ex(jk) = f%Ex(jf,kf)

  jk1=b%nbot2+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    b%Ex(jk) = f%Ex(jf,kf)
  end do
  jk1=b%ntop2+jb
  do kf = f%ny, f%ny+1
    kb = kf+b%nbndy
    jk=jk1+kb*b%n1x
    b%Ex(jk) = f%Ex(jf,kf)
  end do

  jf = 1
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
    if(f%js==1) then
      b%Ey(jk) = f%Ey(jf,kf) - f%Ey_in(kf)
    else
      b%Ey(jk) = f%Ey(jf,kf) 
    end if
  jk1=b%nbot1+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    if(f%js==1) then
      b%Ey(jk) = f%Ey(jf,kf) - f%Ey_in(kf)
    else
      b%Ey(jk) = f%Ey(jf,kf) 
    end if
  enddo
  jk1=b%ntop2+jb
  kf = f%ny
  kb = kf+b%nbndy
  jk=jk1+kb*b%n1x
  if(f%js==1) then
    b%Ey(jk) = f%Ey(jf,kf) - f%Ey_in(kf)
  else
    b%Ey(jk) = f%Ey(jf,kf)
  end if

  jf = f%nx+1
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
  b%Ey(jk) = f%Ey(jf,kf) 
  jk1=b%nbot2+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    b%Ey(jk) = f%Ey(jf,kf)
  enddo
  jk1=b%ntop2+jb
  kf = f%ny
  kb = kf+b%nbndy
  jk=jk1+kb*b%n1x
  b%Ey(jk) = f%Ey(jf,kf)
end if

return
END subroutine exchange_bnd_field

subroutine exchange_bnd_field2(b, f)
 implicit none

TYPE(type_bnd  ) :: b
TYPE(EM2D_FIELDtype) :: f

INTEGER :: jb, kb, jf, kf,jk1,jk

if(f%ylbound==periodic) then
  f%Bx(0:f%nx+1,0)      = f%Bx(0:f%nx+1,f%ny+1)
  f%Bx(0:f%nx+1,f%ny+2) = f%Bx(0:f%nx+1,1)

!else if(.not. (l_moving_window .and. l_elaser_out_plane)) then
else 
kf = 0
kb = kf + b%nbndy
jk1=b%ntop1+kb*b%n1x
do jf = 0, f%nx+1
  jb = jf+b%nbndx
  jk=jk1+jb
  f%Bx(jf,kf) = b%Ex(jk)
end do

kf = f%ny+2
kb = kf + b%nbndy
jk1=b%ntop2+kb*b%n1x
do jf = 0, f%nx+1
  jb = jf+b%nbndx
  jk=jk1+jb
  f%Bx(jf,kf) = b%Ex(jk)
end do

kf = 1
kb = kf + b%nbndy
jk1=b%ntop1+kb*b%n1x
do jf = 1, f%nx
  jb = jf+b%nbndx
  jk=jk1+jb
  if (f%js==1) then
    b%Ex(jk) = f%Bx(jf,kf) - f%Bx_in(jf)
  else
    b%Ex(jk) = f%Bx(jf,kf)
  end if
end do

kf = f%ny+1
kb = kf + b%nbndy
jk1=b%ntop2+kb*b%n1x
do jf = 1, f%nx
  jb = jf+b%nbndx
  jk=jk1+jb
  b%Ex(jk) = f%Bx(jf,kf)
end do
end if

!if(.not. (l_moving_window .and. l_elaser_out_plane)) then
kf = 0
kb = kf + b%nbndy
jk1=b%ntop1+kb*b%n1x
do jf = 1, f%nx+1
  jb = jf+b%nbndx
  jk = jk1+jb
  if(f%js==1) then
    b%Ey(jk) = f%By(jf,kf) - f%By_in(jf)
  else
    b%Ey(jk) = f%By(jf,kf)
  end if
end do

kf = f%ny+1
kb = kf + b%nbndy
jk1=b%ntop2+kb*b%n1x
do jf = 1, f%nx+1
  jb = jf+b%nbndx
  jk = jk1+jb
  b%Ey(jk) = f%By(jf,kf)
end do
!end if

if(f%xlbound==periodic) then
  f%By(0,0:f%ny+1)    = f%By(f%nx+1,0:f%ny+1)
  f%By(f%nx+2,0:f%ny+1) = f%By(1,0:f%ny+1)

else if(.not. (l_moving_window .and. (.not. l_elaser_out_plane))) then
  jf = 0
  jb = jf + b%nbndx
  jk1=b%ntop1+jb
  do kf=0,1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%By(jf,kf)=b%Ey(jk)
  enddo
  jk1=b%nbot1+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    f%By(jf,kf)=b%Ey(jk)
  end do
  jk1=b%ntop2+jb
  do kf=f%ny,f%ny+1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%By(jf,kf)=b%Ey(jk)
  enddo

  jf = f%nx+2
  jb = jf + b%nbndx
  jk1=b%ntop1+jb
  do kf=0,1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%By(jf,kf)=b%Ey(jk)
  enddo
  jk1=b%nbot2+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    f%By(jf,kf)=b%Ey(jk)
  end do
  jk1=b%ntop2+jb
  do kf=f%ny,f%ny+1
    kb=kf+b%nbndy
    jk=jk1+kb*b%n1x
    f%By(jf,kf)=b%Ey(jk)
  enddo

   jf = 0
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
    if(f%js==1) then
      b%Ex(jk) = f%Bx(jf,kf) - f%Bx_in(jf)
    else
      b%Ex(jk) = f%Bx(jf,kf) 
    end if

  jk1=b%nbot1+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    b%Ex(jk) = f%Bx(jf,kf) 
  end do
  jk1=b%ntop2+jb
  do kf = f%ny, f%ny+1
    kb = kf+b%nbndy
    jk=jk1+kb*b%n1x
    b%Ex(jk) = f%Bx(jf,kf) 
  end do

  jf = f%nx+1
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
    if(f%js==1) then
      b%Ex(jk) = f%Bx(jf,kf) - f%Bx_in(jf)
    else
      b%Ex(jk) = f%Bx(jf,kf)
    end if

  jk1=b%nbot2+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    b%Ex(jk) = f%Bx(jf,kf)
  end do
  jk1=b%ntop2+jb
  do kf = f%ny, f%ny+1
    kb = kf+b%nbndy
    jk=jk1+kb*b%n1x
    b%Ex(jk) = f%Bx(jf,kf)
  end do

  jf = 1
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
  b%Ey(jk) = f%By(jf,kf) 
  jk1=b%nbot1+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    b%Ey(jk) = f%By(jf,kf) 
  enddo
  jk1=b%ntop2+jb
  kf = f%ny
  kb = kf+b%nbndy
  jk=jk1+kb*b%n1x
  b%Ey(jk) = f%By(jf,kf)

  jf = f%nx+1
  jb = jf + b%nbndx
  kf=1
  jk=b%ntop1+jb+b%n1x*(kf+b%nbndy)
  b%Ey(jk) = f%By(jf,kf) 
  jk1=b%nbot2+jb
  do kf = 2, f%ny-1
    kb = kf+b%nbndy
    jk=jk1+kb*b%nint
    b%Ey(jk) = f%By(jf,kf)
  enddo
  jk1=b%ntop2+jb
  kf = f%ny
  kb = kf+b%nbndy
  jk=jk1+kb*b%n1x
  b%Ey(jk) = f%By(jf,kf)
end if

return
END subroutine exchange_bnd_field2

!*************************  SUBROUTINE griuni****************************************

subroutine griuni(f)

implicit none
TYPE(EM2D_FIELDtype) :: f

INTEGER :: which,i,j

!     ce sous programme met tous les champs sur une 
!     seule grille
      do  i=f%nx+1,1,-1
      do  j=0,f%ny+1
      f%ex(i,j)=0.5*(f%ex(i,j)+f%ex(i-1,j))
      f%ez(i,j)=0.5*(f%ez(i,j)+f%ez(i-1,j))
      f%bx(i,j)=0.5*(f%bx(i,j)+f%bx(i-1,j))
      f%bz(i,j)=0.5*(f%bz(i,j)+f%bz(i-1,j))
      enddo
      enddo

      do j=f%ny+1,1,-1
      do i=1,f%nx+1
      f%ey(i,j)=0.5*(f%ey(i,j)+f%ey(i,j-1))
      f%ez(i,j)=0.5*(f%ez(i,j)+f%ez(i,j-1))
      f%by(i,j)=0.5*(f%by(i,j)+f%by(i,j-1))
      f%bz(i,j)=0.5*(f%bz(i,j)+f%bz(i,j-1))
      enddo
      enddo

      return
 end subroutine griuni


!*************************  SUBROUTINE griuni****************************************

      subroutine grimax(f)

implicit none
TYPE(EM2D_FIELDtype) :: f

INTEGER :: which,i,j

!     ce sous programme defait le travail de griuni

      do j=1,f%ny+1
      do i=1,f%nx+1
      f%ey(i,j)=2.*f%ey(i,j)-f%ey(i,j-1)
      f%ez(i,j)=2.*f%ez(i,j)-f%ez(i,j-1)
      f%by(i,j)=2.*f%by(i,j)-f%by(i,j-1)
      f%bz(i,j)=2.*f%bz(i,j)-f%bz(i,j-1)
      enddo
      enddo

      do  i=1,f%nx+1
      do  j=0,f%ny+1
      f%ex(i,j)=2.*f%ex(i,j)-f%ex(i-1,j)
      f%ez(i,j)=2.*f%ez(i,j)-f%ez(i-1,j)
      f%bx(i,j)=2.*f%bx(i,j)-f%bx(i-1,j)
      f%bz(i,j)=2.*f%bz(i,j)-f%bz(i-1,j)
      enddo
      enddo

     
      return
      end subroutine grimax

 subroutine smooth(f,q,nx,ny)
 implicit none

 integer :: nx,ny,ns,i1,i2,j1,j2,is,i,j

 real(kind=8), dimension(:,:) :: q
 real(kind=8), dimension(5) :: cs,ds,dc

 data cs /4*.25,-1.25/,ds/4*.5,3.5/,ns/5/
 data dc /4*2.,-2.8/

 TYPE(EM2D_FIELDtype) :: f


      i1=0
      i2=nx+2
      j1=0
      j2=ny+2

      f%temp=0.

!     x smoothing

      do 110 is=1,ns
      do  i=2,nx-1,2
!cdir nodep
      do  j=1,ny+1
      f%temp(j+j1)=q(i-1,j)+dc(is)*q(i,j)+q(i+1,j)
      q(i-1,j)=cs(is)*f%temp(j+j2)
      f%temp(j+j2)=q(i,j)+dc(is)*q(i+1,j)+q(i+2,j)
      q(i,j)=cs(is)*f%temp(j+j1)

      enddo
      enddo

      do  j=1,ny+1
      q(nx,j)=cs(is)*f%temp(j+j2)
      q(1,j)=0.
      q(nx+1,j)=0.
      enddo

 110  continue

!     y smoothing
!     -----------
      do 160 is=1,ns

      do j=2,ny-1,2

!cdir nodep
      do  i=1,nx+1
      f%temp(i+i1)=q(i,j-1)+dc(is)*q(i,j)+q(i,j+1)
      q(i,j-1)=cs(is)*f%temp(i+i2)
      f%temp(i+i2)=q(i,j)+dc(is)*q(i,j+1)+q(i,j+2)
      q(i,j)=cs(is)*f%temp(i+i1)
      enddo
      enddo

      do  i=1,nx+1
      q(i,ny)=cs(is)*f%temp(i+i2)
      q(i,1)=0.
      q(i,ny+1)=0.
      enddo

 160  continue

      return
      end subroutine smooth
   subroutine project_jxjy(jxjyfin,jxjycoarse,rap)
   ! Routine de projection des J d'une grille fine sur une grille coarsesiere.
   ! Soit nx*ny la taille de la grille coarsesiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive jxjycoarse(0:nx+1,0:ny+1,2) et jxjyfin(0:rap*nx+1,0:rap*ny+1,2).
   real(kind=8), DIMENSION(0:,0:,:) :: jxjyfin,jxjycoarse
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, j, k, jg, kg
   real(kind=8) :: w, invrapvol

!      invrapvol = 1./rap!**2
      invrapvol = 1./rap**2

      nxfin = SIZE(jxjyfin,1)-2
      nyfin = SIZE(jxjyfin,2)-2
      nxcoarse = SIZE(jxjycoarse,1)-2
      nycoarse = SIZE(jxjycoarse,2)-2

      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in project_jxjy: rap does not match grid sizes.")
        return
      END if

      do k = 0, nyfin
        kg = k/rap
        w = REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          jxjycoarse(jg,kg,  1) = jxjycoarse(jg,kg,  1) + (1.-w)*jxjyfin(j,k,1)*invrapvol
          jxjycoarse(jg,kg+1,1) = jxjycoarse(jg,kg+1,1) +     w *jxjyfin(j,k,1)*invrapvol
        end do
      end do

      do k = 0, nyfin-1
        kg = k/rap
        do j = 0, nxfin
          jg = j/rap
          w = REAL(MOD(j,rap))/rap
          jxjycoarse(jg,  kg,2) = jxjycoarse(jg  ,kg,2) + (1.-w)*jxjyfin(j,k,2)*invrapvol
          jxjycoarse(jg+1,kg,2) = jxjycoarse(jg+1,kg,2) +     w *jxjyfin(j,k,2)*invrapvol
        end do
      end do

   end subroutine project_jxjy

   subroutine interpol_jxjy(jxjyfin,jxjycoarse,rap)
   ! Routine d'interpolation-soustraction des J d'une grille coarsesiere sur une grille fine.
   ! Soit nx*ny la taille de la grille coarsesiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive jxjycoarse(0:nx+1,0:ny+1,2) et jxjyfin(0:rap*nx+1,0:rap*ny+1,2).
   real(kind=8), DIMENSION(0:,0:,:) :: jxjyfin,jxjycoarse
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, j, k, jg, kg
   real(kind=8) :: w
   logical :: l_onfinegridfirst = .true.
   real(kind=8), DIMENSION(:,:,:), allocatable :: jxjyfin_new

      nxfin = SIZE(jxjyfin,1)-2
      nyfin = SIZE(jxjyfin,2)-2
      nxcoarse = SIZE(jxjycoarse,1)-2
      nycoarse = SIZE(jxjycoarse,2)-2

      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in interpol_jxjy: rap does not match grid sizes.")
        return
      END if

      if(.not.l_onfinegridfirst) then
        do k = 0, nyfin
          kg = k/rap
          w = REAL(MOD(k,rap))/rap
          do j = 0, nxfin-1
            jg = j/rap
            jxjyfin(j,k,1) = jxjyfin(j,k,1)-((1.-w)*jxjycoarse(jg,kg,1)+w*jxjycoarse(jg,kg+1,1))/rap
          end do
        end do

        do k = 0, nyfin-1
          kg = k/rap
          do j = 0, nxfin
            jg = j/rap
            w = REAL(MOD(j,rap))/rap
            jxjyfin(j,k,2) = jxjyfin(j,k,2)-((1.-w)*jxjycoarse(jg,kg,2)+w*jxjycoarse(jg+1,kg,2))/rap
          end do
        end do
      else
        allocate(jxjyfin_new(0:nxfin+1,0:nyfin+1,2))
	jxjyfin_new=0.
        do k = 0, nyfin
          do j = 0, nxfin-1
!            jxjyfin_new(j,  k,  1) = jxjyfin_new(j,  k,  1) +      jxjyfin(j,k,1)
            jxjyfin_new(j+1,k,  1) = jxjyfin_new(j+1,k,  1) - 0.25*jxjyfin(j,k,1)
            jxjyfin_new(j-1,k,  1) = jxjyfin_new(j-1,k,  1) - 0.25*jxjyfin(j,k,1)
            jxjyfin_new(j,  k+1,1) = jxjyfin_new(j,  k+1,1) - 0.25*jxjyfin(j,k,1)
            jxjyfin_new(j,  k-1,1) = jxjyfin_new(j,  k-1,1) - 0.25*jxjyfin(j,k,1)
          end do
        end do

        do k = 0, nyfin-1
          do j = 0, nxfin
!            jxjyfin_new(j,  k,  2) = jxjyfin_new(j,  k,  2) +      jxjyfin(j,k,2)
            jxjyfin_new(j+1,k,  2) = jxjyfin_new(j+1,k,  2) - 0.25*jxjyfin(j,k,2)
            jxjyfin_new(j-1,k,  2) = jxjyfin_new(j-1,k,  2) - 0.25*jxjyfin(j,k,2)
            jxjyfin_new(j,  k+1,2) = jxjyfin_new(j,  k+1,2) - 0.25*jxjyfin(j,k,2)
            jxjyfin_new(j,  k-1,2) = jxjyfin_new(j,  k-1,2) - 0.25*jxjyfin(j,k,2)
          end do
        end do	
	
	jxjyfin = 0.5*jxjyfin_new
	deallocate(jxjyfin_new)
      end if
   end subroutine interpol_jxjy

   subroutine project_ex(exfin,excoarse,rap)
   ! Routine de projection des J d'une grille fine sur une grille coarsesiere.
   ! Soit nx*ny la taille de la grille coarsesiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive excoarse(0:nx+1,0:ny+1,2) et exfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: exfin,excoarse
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, j, k, jg, kg
   real(kind=8) :: w, invrapvol

      invrapvol = 1./rap**2

      nxfin = SIZE(exfin,1)-2
      nyfin = SIZE(exfin,2)-2
      nxcoarse = SIZE(excoarse,1)-2
      nycoarse = SIZE(excoarse,2)-2

      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in project_ex: rap does not match grid sizes.")
        return
      END if

      do k = 0, nyfin
        kg = k/rap
        w = REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          excoarse(jg,kg  ) = excoarse(jg,kg  ) + (1.-w)*exfin(j,k)*invrapvol
          excoarse(jg,kg+1) = excoarse(jg,kg+1) +     w *exfin(j,k)*invrapvol
        end do
      end do

   end subroutine project_ex

   subroutine project_ey(eyfin,eycoarse,rap)
   ! Routine de projection des J d'une grille fine sur une grille coarsesiere.
   ! Soit nx*ny la taille de la grille coarsesiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive eycoarse(0:nx+1,0:ny+1,2) et eyfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: eyfin,eycoarse
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, j, k, jg, kg
   real(kind=8) :: w, invrapvol

      invrapvol = 1./rap**2

      nxfin = SIZE(eyfin,1)-2
      nyfin = SIZE(eyfin,2)-2
      nxcoarse = SIZE(eycoarse,1)-2
      nycoarse = SIZE(eycoarse,2)-2

      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in project_ey: rap does not match grid sizes.")
        return
      END if

      do k = 0, nyfin-1
        kg = k/rap
        do j = 0, nxfin
          jg = j/rap
          w = REAL(MOD(j,rap))/rap
          eycoarse(jg,  kg) = eycoarse(jg  ,kg) + (1.-w)*eyfin(j,k)*invrapvol
          eycoarse(jg+1,kg) = eycoarse(jg+1,kg) +     w *eyfin(j,k)*invrapvol
        end do
      end do

   end subroutine project_ey

   subroutine project_rho(bzfin,bzcoarse,rap)
   ! Routine de projection de Bz d'une grille fine sur une grille coarsesiere.
   ! Soit nx*ny la taille de la grille coarsesiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive bzcoarse(0:nx+1,0:ny+1) et bzfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: bzfin,bzcoarse
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, j, k, jg, kg
   real(kind=8) :: wx, wy, invrapvol,q

      invrapvol = 1./rap**2

      nxfin = SIZE(bzfin,1)-2
      nyfin = SIZE(bzfin,2)-2
      nxcoarse = SIZE(bzcoarse,1)-2
      nycoarse = SIZE(bzcoarse,2)-2

      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in project_bz: rap does not match grid sizes.")
        return
      END if

      do k = 0, nyfin-1
        kg = k/rap
        wy = REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          wx = REAL(MOD(j,rap))/rap
          q = bzfin(j,k)*invrapvol
          bzcoarse(jg,  kg  ) = bzcoarse(jg,  kg  ) + (1.-wx) * (1.-wy) * q
          bzcoarse(jg+1,kg  ) = bzcoarse(jg+1,kg  ) +     wx  * (1.-wy) * q
          bzcoarse(jg  ,kg+1) = bzcoarse(jg,  kg+1) + (1.-wx) *     wy  * q
          bzcoarse(jg+1,kg+1) = bzcoarse(jg+1,kg+1) +     wx  *     wy  * q
        end do
      end do

   end subroutine project_rho

   subroutine project_jz(bzfin,bzcoarse,rap)
   ! Routine de projection de Bz d'une grille fine sur une grille coarsesiere.
   ! Soit nx*ny la taille de la grille coarsesiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive bzcoarse(0:nx+1,0:ny+1) et bzfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: bzfin,bzcoarse
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, j, k, jg, kg
   real(kind=8) :: wx, wy, invrapvol,q

      invrapvol = 1./rap**2

      nxfin = SIZE(bzfin,1)-2
      nyfin = SIZE(bzfin,2)-2
      nxcoarse = SIZE(bzcoarse,1)-2
      nycoarse = SIZE(bzcoarse,2)-2

      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in project_bz: rap does not match grid sizes.")
        return
      END if

      do k = 0, nyfin-1
        kg = k/rap
        wy = 0.25+REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          wx = 0.25+REAL(MOD(j,rap))/rap
          q = bzfin(j,k)*invrapvol
          bzcoarse(jg,  kg  ) = bzcoarse(jg,  kg  ) + (1.-wx) * (1.-wy) * q
          bzcoarse(jg+1,kg  ) = bzcoarse(jg+1,kg  ) +     wx  * (1.-wy) * q
          bzcoarse(jg  ,kg+1) = bzcoarse(jg,  kg+1) + (1.-wx) *     wy  * q
          bzcoarse(jg+1,kg+1) = bzcoarse(jg+1,kg+1) +     wx  *     wy  * q
        end do
      end do

   end subroutine project_jz

   subroutine interpol(ffin,fcoarse,rap)
   real(kind=8), DIMENSION(0:,0:) :: ffin,fcoarse
   INTEGER :: rap
   INTEGER :: nxfin, nyfin, nxcoarse, nycoarse, jf, kf, jg, kg
   real(kind=8) :: wj, wk

      nxfin = SIZE(ffin,1)-2
      nyfin = SIZE(ffin,2)-2
      nxcoarse = SIZE(fcoarse,1)-2
      nycoarse = SIZE(fcoarse,2)-2
      IF(nxcoarse*rap/=nxfin .OR. nycoarse*rap/=nyfin) then
        call kaboom("Error in interpol: rap does not match grid sizes.")
        return
      END if

      do kf = 0, nyfin
        kg = kf/rap
        wk = REAL(MOD(kf,rap))/rap
        do jf = 0, nxfin
          jg = jf/rap
          wj = REAL(MOD(jf,rap))/rap
          ffin(jf,kf) = ffin(jf,kf) + (1.-wj)*(1.-wk)*fcoarse(jg,  kg)   &
                                    +     wj *(1.-wk)*fcoarse(jg+1,kg)   &
                                    + (1.-wj)*    wk *fcoarse(jg,  kg+1) &
                                    +     wj *    wk *fcoarse(jg+1,kg+1)
        end do
      end do


   END subroutine interpol

subroutine shift_fields_1cell(f)
TYPE(EM2D_FIELDtype) :: f
    integer(4):: ix,iy

!   --- Note that the loops are written out since they seem to give
!   --- the intel compiler fits.
    if (l_elaser_out_plane .and. .false.) then
      do iy=0,f%ny+2-f%rap
        do ix=0,f%nx+3
          f%Ex(ix,iy) = f%Ex(ix,iy+f%rap)
          f%Ey(ix,iy) = f%Ey(ix,iy+f%rap)
          f%Ez(ix,iy) = f%Ez(ix,iy+f%rap)
          f%Bx(ix,iy) = f%Bx(ix,iy+f%rap)
          f%By(ix,iy) = f%By(ix,iy+f%rap)
          f%Bz(ix,iy) = f%Bz(ix,iy+f%rap)
          f%J(ix,iy,1) = f%J(ix,iy+f%rap,1)
          f%J(ix,iy,2) = f%J(ix,iy+f%rap,2)
          f%J(ix,iy,3) = f%J(ix,iy+f%rap,3)
        enddo
      enddo

      do iy=f%ny+2-f%rap+1,f%ny+2
        do ix=0,f%nx+3
          f%Ex(ix,iy) = 0.
          f%Ey(ix,iy) = 0.
          f%Ez(ix,iy) = 0.
          f%Bx(ix,iy) = 0.
          f%By(ix,iy) = 0.
          f%Bz(ix,iy) = 0.
          f%J(ix,iy,1) = 0.
          f%J(ix,iy,2) = 0.
          f%J(ix,iy,3) = 0.
        enddo
      enddo

    else
      do iy=0,f%ny+2
        do ix=0,f%nx+3-f%rap
          f%Ex(ix,iy) = f%Ex(ix+f%rap,iy)
          f%Ey(ix,iy) = f%Ey(ix+f%rap,iy)
          f%Ez(ix,iy) = f%Ez(ix+f%rap,iy)
          f%Bx(ix,iy) = f%Bx(ix+f%rap,iy)
          f%By(ix,iy) = f%By(ix+f%rap,iy)
          f%Bz(ix,iy) = f%Bz(ix+f%rap,iy)
          f%J(ix,iy,1) = f%J(ix+f%rap,iy,1)
          f%J(ix,iy,2) = f%J(ix+f%rap,iy,2)
          f%J(ix,iy,3) = f%J(ix+f%rap,iy,3)
        enddo
      enddo

      do iy=0,f%ny+2
        do ix=f%nx+3-f%rap+1,f%nx+3
          f%Ex(ix,iy) = 0.
          f%Ey(ix,iy) = 0.
          f%Ez(ix,iy) = 0.
          f%Bx(ix,iy) = 0.
          f%By(ix,iy) = 0.
          f%Bz(ix,iy) = 0.
          f%J(ix,iy,1) = 0.
          f%J(ix,iy,2) = 0.
          f%J(ix,iy,3) = 0.
        enddo
      enddo
    end if

end subroutine shift_fields_1cell
end module mod_field

subroutine move_window_field(f)
use EM2D_FIELDtypemodule
use mod_field,Only: l_elaser_out_plane,shift_fields_1cell
use mod_bnd
TYPE(EM2D_FIELDtype):: f
  call shift_fields_1cell(f)
  call move_window_bnd(f%bndexeybz,f%rap,.false.)!l_elaser_out_plane)
  call move_window_bnd(f%bndbxbyez,f%rap,.false.)!l_elaser_out_plane)
  if (l_elaser_out_plane .and. .false.) then
    f%ymin  =f%ymin  +f%dy*f%rap
    f%ymax  =f%ymax  +f%dy*f%rap
  else
    f%xmin  =f%xmin  +f%dx*f%rap
    f%xmax  =f%xmax  +f%dx*f%rap
  end if
end subroutine move_window_field

subroutine smooth2(q,nx,ny)
 implicit none

 integer :: nx,ny,ns,i1,i2,j1,j2,is,i,j,ntemp

 real(kind=8), dimension(0:nx+2,0:ny+2) :: q
 real(kind=8), dimension(5) :: cs,ds,dc
 real(kind=8), dimension(:), ALLOCATABLE :: temp

 data cs /4*.25,-1.25/,ds/4*.5,3.5/,ns/5/
 data dc /4*2.,-2.8/

     ntemp = 2*max(nx,ny)+4
     ALLOCATE(temp(0:ntemp))


      i1=0
      i2=nx+2
      j1=0
      j2=ny+2

      temp=0.

!     x smoothing

      do 110 is=1,ns
      do  i=2,nx-1,2
!cdir nodep
      do  j=1,ny+1
      temp(j+j1)=q(i-1,j)+dc(is)*q(i,j)+q(i+1,j)
      q(i-1,j)=cs(is)*temp(j+j2)
      temp(j+j2)=q(i,j)+dc(is)*q(i+1,j)+q(i+2,j)
      q(i,j)=cs(is)*temp(j+j1)

      enddo
      enddo

      do  j=1,ny+1
      q(nx,j)=cs(is)*temp(j+j2)
      q(1,j)=0.
      q(nx+1,j)=0.
      enddo

 110  continue

!     y smoothing
!     -----------
      do 160 is=1,ns

      do j=2,ny-1,2

!cdir nodep
      do  i=1,nx+1
      temp(i+i1)=q(i,j-1)+dc(is)*q(i,j)+q(i,j+1)
      q(i,j-1)=cs(is)*temp(i+i2)
      temp(i+i2)=q(i,j)+dc(is)*q(i,j+1)+q(i,j+2)
      q(i,j)=cs(is)*temp(i+i1)
      enddo
      enddo

      do  i=1,nx+1
      q(i,ny)=cs(is)*temp(i+i2)
      q(i,1)=0.
      q(i,ny+1)=0.
      enddo

 160  continue
      DEALLOCATE(temp)

      return
      end subroutine smooth2


!subroutine initfields(f,nx, ny, nbndx, nbndy, dtm, dx, dy, xmin, ymin, rap, xlb, ylb, xrb, yrb)
!use mod_field, only:init_fields, EM2D_FIELDtype
!implicit none

!TYPE(EM2D_FIELDtype), pointer :: f
!INTEGER(ISZ), INTENT(IN) :: nx, ny, rap
!INTEGER(ISZ), INTENT(IN) :: nbndx, nbndy, xlb, ylb, xrb, yrb
!REAL(kind=8), INTENT(IN) :: dtm, dx, dy, xmin, ymin

! call init_fields(f,nx, ny, nbndx, nbndy, dtm, dx, dy, xmin, ymin, rap, xlb, ylb, xrb, yrb)

! return
!end subroutine initfields


!************* SUBROUTINE init_fields  *************************************************
subroutine init_fields(f,nx, ny, nbndx, nbndy, dt, dx, dy, clight, mu0, xmin, ymin, rap, xlb, ylb, xrb, yrb)
use mod_bnd
use mod_field, only:EM2D_FIELDtype, l_copyfields, l_elaser_out_plane
implicit none

TYPE(EM2D_FIELDtype) :: f
INTEGER(ISZ), INTENT(IN) :: nx, ny, rap
INTEGER(ISZ), INTENT(IN) :: nbndx, nbndy, xlb, ylb, xrb, yrb
REAL(kind=8), INTENT(IN) :: dt, dx, dy, xmin, ymin, clight, mu0
INTEGER :: k,m

!f => NewEM2D_FIELDType()
f%bndexeybz => Newtype_bnd()
f%bndbxbyez => Newtype_bnd()

f%l_apply_pml=.true.
f%nx = nx
f%ny = ny
f%xmin = xmin
f%ymin = ymin
f%rap = rap
f%dx = dx
f%dy = dy
f%xmax = xmin+f%dx*f%nx
f%ymax = ymin+f%dy*f%ny
f%dxi = 1./dx
f%dyi = 1./dy
f%npulse=300
f%ntemp = 2*max(nx,ny)+4
f%clight = clight
f%mu0    = mu0
f%js = 1 ! position of the source
if(rap>1) then
  f%nxfsum = nx
  f%nyfsum = ny
end if

  IF(l_copyfields) then
    f%nxcopy = f%nx
    f%nycopy = f%ny
  ELSE
    f%nxcopy = 0
    f%nycopy = 0
  END if
!  call EM2D_FIELDtypeallot(f)

call create_bnd(f%bndexeybz, nx, ny, nbndx=nbndx, nbndy=nbndy, dt=dt*clight, dx=dx, dy=dy, xbnd=xlb, ybnd=ylb)
call create_bnd(f%bndbxbyez, nx, ny, nbndx=nbndx, nbndy=nbndy, dt=dt*clight, dx=dx, dy=dy, xbnd=xlb, ybnd=ylb)
  
  call EM2D_FIELDtypeallot(f)


	f%Ex = 0.
	f%Ey = 0.
	f%Ez = 0.
	f%Bx = 0.
	f%By = 0.
	f%Bz = 0.
	
	f%J = 0.

	f%Bz_in = 0.
	f%Ey_in = 0.
	f%Ex_in = 0.
	f%Ez_in = 0.
	f%By_in = 0.
	f%Bx_in = 0.
      f%pulse=0.
      f%tpulse=0.

f%xlbound = xlb
f%xrbound = xrb
f%ylbound = ylb
f%yrbound = yrb

return

END subroutine init_fields

subroutine push_em_e(f,dt)
use mod_field, only: champ_e, EM2D_FIELDtype
implicit none

TYPE(EM2D_FIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

call champ_e(f,dt)

return
end subroutine push_em_e

subroutine push_em_b(f,dt)
use mod_field, only: champ_b, EM2D_FIELDtype
implicit none

TYPE(EM2D_FIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

call champ_b(f,dt)

return
end subroutine push_em_b

!*************************  SUBROUTINE griuni****************************************

subroutine griuni(f)
! average fields at nodes locations

use mod_field, only: EM2D_FIELDtype
implicit none
INTEGER :: i,j
TYPE(EM2D_FIELDtype) :: f

  do  i=f%nx+1,1,-1
    do  j=0,f%ny+1
      f%ex(i,j)=0.5*(f%ex(i,j)+f%ex(i-1,j))
      f%ez(i,j)=0.5*(f%ez(i,j)+f%ez(i-1,j))
      f%bx(i,j)=0.5*(f%bx(i,j)+f%bx(i-1,j))
      f%bz(i,j)=0.5*(f%bz(i,j)+f%bz(i-1,j))
    enddo
  enddo

  do j=f%ny+1,1,-1
    do i=1,f%nx+1
      f%ey(i,j)=0.5*(f%ey(i,j)+f%ey(i,j-1))
      f%ez(i,j)=0.5*(f%ez(i,j)+f%ez(i,j-1))
      f%by(i,j)=0.5*(f%by(i,j)+f%by(i,j-1))
      f%bz(i,j)=0.5*(f%bz(i,j)+f%bz(i,j-1))
    enddo
  enddo

  return
end subroutine griuni


!*************************  SUBROUTINE grimax****************************************

subroutine grimax(f)
! undo of griuni (puts E and B staggered values back)

use mod_field, only: EM2D_FIELDtype
implicit none
INTEGER :: i,j
TYPE(EM2D_FIELDtype) :: f

  do j=1,f%ny+1
    do i=1,f%nx+1
      f%ey(i,j)=2.*f%ey(i,j)-f%ey(i,j-1)
      f%ez(i,j)=2.*f%ez(i,j)-f%ez(i,j-1)
      f%by(i,j)=2.*f%by(i,j)-f%by(i,j-1)
      f%bz(i,j)=2.*f%bz(i,j)-f%bz(i,j-1)
    enddo
  enddo

  do i=1,f%nx+1
    do  j=0,f%ny+1
      f%ex(i,j)=2.*f%ex(i,j)-f%ex(i-1,j)
      f%ez(i,j)=2.*f%ez(i,j)-f%ez(i-1,j)
      f%bx(i,j)=2.*f%bx(i,j)-f%bx(i-1,j)
      f%bz(i,j)=2.*f%bz(i,j)-f%bz(i-1,j)
    enddo
  enddo
    
  return
end subroutine grimax

subroutine project_j(f,fc,ff)
   use mod_field, only: EM2D_FIELDtype, project_jxjy, project_jz
   implicit none
   
   TYPE(EM2D_FIELDtype) :: f,fc, ff
   integer :: nxpatch, nypatch, ixpatch, iypatch

   nxpatch = fc%nx
   nypatch = fc%ny
   ixpatch = nint((fc%xmin-f%xmin)/f%dx)
   iypatch = nint((fc%ymin-f%ymin)/f%dy)
   call project_jxjy(        ff%J(1:ff%nx+2, 1:ff%ny+2, 1:2), fc%J(1:fc%nx+2, 1:fc%ny+2, 1:2), ff%rap)
   call project_jz(ff%J(0:ff%nx+1, 0:ff%ny+1, 3),   fc%J(0:fc%nx+1, 0:fc%ny+1, 3),   ff%rap)
   f%J(ixpatch:ixpatch+nxpatch+3,iypatch:iypatch+nypatch+2,:) = &
   f%J(ixpatch:ixpatch+nxpatch+3,iypatch:iypatch+nypatch+2,:) + fc%J

   return
end subroutine project_j

subroutine set_substitute_fields(field,fpatchcoarse,fpatchfine)
use mod_field, only: EM2D_FIELDtype, interpol
implicit none

TYPE(EM2D_FIELDtype) :: field,fpatchcoarse,fpatchfine

integer :: nxfsum, nyfsum, nxpatch, nypatch, ixpatch, iypatch, rap

  nxfsum = fpatchfine%nxfsum
  nyfsum = fpatchfine%nyfsum
  nxpatch = fpatchcoarse%nx
  nypatch = fpatchcoarse%ny
  ixpatch = nint((fpatchcoarse%xmin-field%xmin)/field%dx)+1
  iypatch = nint((fpatchcoarse%ymin-field%ymin)/field%dy)+1
  rap = fpatchfine%rap
  
    fpatchfine%exfsum = fpatchfine%ex
    fpatchfine%eyfsum = fpatchfine%ey 
    fpatchfine%ezfsum = fpatchfine%ez
    fpatchfine%bxfsum = fpatchfine%bx
    fpatchfine%byfsum = fpatchfine%by 
    fpatchfine%bzfsum = fpatchfine%bz

    call interpol(ffin=fpatchfine%exfsum(1:nxfsum+2,1:nyfsum+2),fcoarse=-fpatchcoarse%ex(1:nxpatch+2,1:nypatch+2),rap=rap)
    call interpol(ffin=fpatchfine%eyfsum(1:nxfsum+2,1:nyfsum+2),fcoarse=-fpatchcoarse%ey(1:nxpatch+2,1:nypatch+2),rap=rap)
    call interpol(ffin=fpatchfine%ezfsum(1:nxfsum+2,1:nyfsum+2),fcoarse=-fpatchcoarse%ez(1:nxpatch+2,1:nypatch+2),rap=rap)
    call interpol(ffin=fpatchfine%bxfsum(1:nxfsum+2,1:nyfsum+2),fcoarse=-fpatchcoarse%bx(1:nxpatch+2,1:nypatch+2),rap=rap)
    call interpol(ffin=fpatchfine%byfsum(1:nxfsum+2,1:nyfsum+2),fcoarse=-fpatchcoarse%by(1:nxpatch+2,1:nypatch+2),rap=rap)
    call interpol(ffin=fpatchfine%bzfsum(1:nxfsum+2,1:nyfsum+2),fcoarse=-fpatchcoarse%bz(1:nxpatch+2,1:nypatch+2),rap=rap)

    call interpol(ffin=fpatchfine%exfsum(1:nxfsum+2,1:nyfsum+2), &
                  fcoarse=field%ex(ixpatch:ixpatch+nxpatch+1,iypatch:iypatch+nypatch+1),rap=rap)
    call interpol(ffin=fpatchfine%eyfsum(1:nxfsum+2,1:nyfsum+2), &
                  fcoarse=field%ey(ixpatch:ixpatch+nxpatch+1,iypatch:iypatch+nypatch+1),rap=rap)
    call interpol(ffin=fpatchfine%ezfsum(1:nxfsum+2,1:nyfsum+2), &
                  fcoarse=field%ez(ixpatch:ixpatch+nxpatch+1,iypatch:iypatch+nypatch+1),rap=rap)
    call interpol(ffin=fpatchfine%bxfsum(1:nxfsum+2,1:nyfsum+2), &
                  fcoarse=field%bx(ixpatch:ixpatch+nxpatch+1,iypatch:iypatch+nypatch+1),rap=rap)
    call interpol(ffin=fpatchfine%byfsum(1:nxfsum+2,1:nyfsum+2), &
                  fcoarse=field%by(ixpatch:ixpatch+nxpatch+1,iypatch:iypatch+nypatch+1),rap=rap)
    call interpol(ffin=fpatchfine%bzfsum(1:nxfsum+2,1:nyfsum+2), &
                  fcoarse=field%bz(ixpatch:ixpatch+nxpatch+1,iypatch:iypatch+nypatch+1),rap=rap)

  return
end subroutine set_substitute_fields
