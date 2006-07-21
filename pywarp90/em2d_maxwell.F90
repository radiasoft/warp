#include "top.h"
!     Last change:  JLV   3 Jun 2004    0:17 am
!************* MODULE field  **********************************************

module mod_field

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
if (f%l_add_source) then
 if (.not.l_moving_window) then
  if(l_elaser_out_plane) then
  ! add contribution from Ez source (Ez_in)
    k = f%js
    do j = 0, f%nx+1
      f%Bx(j,k) = f%Bx(j,k) - dtsdy * (-f%Ez_in(j))
    end do

    ! Evaluate Bx and By source
    do j = 1, f%nx+1
      f%By_in(j) = f%By_in(j) + dtsdx * (f%Ez_in(j) - f%Ez_in(j-1))
    end do

    do j = 0, f%nx+1
      f%Bx_in(j) = f%cst1 * f%Bx_in(j) + f%cst2 * f%Ez_in(j)
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

IF(f%l_addpatchresidual) then
  ALLOCATE(Exapr(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1))
  ALLOCATE(Eyapr(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1))
  Exapr = 0.
  Eyapr = 0.
  call project_ex(exfin=fpatchfine%ex(0:fpatchfine%nx+1,0:fpatchfine%ny+1), exgros=Exapr,rap=rap)
  call project_ey(eyfin=fpatchfine%ey(0:fpatchfine%nx+1,0:fpatchfine%ny+1), eygros=Eyapr,rap=rap)
  Exapr = Exapr - fpatchcoarse%ex(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1)
  Eyapr = Eyapr - fpatchcoarse%ey(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1)
  j = ntamp_apr
  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) - dtsdx * Eyapr(j,k)
  END do
  j = fpatchcoarse%nx+1-ntamp_apr
  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) + dtsdx * Eyapr(j+1,k)
  END do
  k = ntamp_apr
  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) + dtsdy * Exapr(j,k)
  END do
  k = fpatchcoarse%ny+1-ntamp_apr
  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
    f%Bz(j+ixpatch,k+iypatch) = f%Bz(j+ixpatch,k+iypatch) - dtsdy * Exapr(j,k+1)
  END do
  DEALLOCATE(Exapr,Eyapr)
END if

end subroutine champ_b


subroutine champ_e(f,dt)
implicit none

TYPE(EM2D_FIELDtype) :: f
INTEGER :: j, k
real(kind=8) :: dt,dtsdx,dtsdy
real(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Bzapr

dtsdx = dt/f%dx
dtsdy = dt/f%dy

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
                          - f%J(j,k,1)
  end do
end do

! advance Ey
do k = 0, f%ny+1
  do j = 1, f%nx+1
    f%Ey(j,k) = f%Ey(j,k) - dtsdx * (f%Bz(j,k)   - f%Bz(j-1,k)) &
                          - f%J(j,k,2)
  end do
end do

! advance Ez 
do k = 0, f%ny+1
  do j = 0, f%nx+1
    f%Ez(j,k) = f%Ez(j,k) + dtsdx * (f%By(j+1,k) - f%By(j,k)) &
                          - dtsdy * (f%Bx(j,k+1) - f%Bx(j,k)) &
                          - f%J(j,k,3)
  end do
end do

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
    ! substract contribution from Bx source (Bx_in)
    k = f%js-1
    do j = 0, f%nx+1
      f%Ez(j,k) = f%Ez(j,k) - dtsdy * (-f%Bx_in(j)) 
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
    k = 0
    do j = 0, f%nx+1
      f%Ez(j,k) = f%Ez_in(j) 
    end do
  end if
 end if
end if

IF(f%l_addpatchresidual) then
  ALLOCATE(Bzapr(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1))
  Bzapr = 0. 
  call project_bz(bzfin=fpatchfine%bz(0:fpatchfine%nx+1,0:fpatchfine%ny+1), bzgros=Bzapr,rap=rap)
  Bzapr = Bzapr - fpatchcoarse%bz(0:fpatchcoarse%nx+1,0:fpatchcoarse%ny+1)
  
  j = ntamp_apr
  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
    f%Ey(j+ixpatch,k+iypatch) = f%Ey(j+ixpatch,k+iypatch) - dtsdx * Bzapr(j,k)
  END do
  j = fpatchcoarse%nx+2-ntamp_apr
  do k = ntamp_apr, fpatchcoarse%ny+1-ntamp_apr
    f%Ey(j+ixpatch,k+iypatch) = f%Ey(j+ixpatch,k+iypatch) + dtsdx * Bzapr(j-1,k)
  END do
  k = ntamp_apr
  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
    f%Ex(j+ixpatch,k+iypatch) = f%Ex(j+ixpatch,k+iypatch) + dtsdy * Bzapr(j,k)
  END do
  k = fpatchcoarse%ny+2-ntamp_apr
  do j = ntamp_apr, fpatchcoarse%nx+1-ntamp_apr
    f%Ex(j+ixpatch,k+iypatch) = f%Ex(j+ixpatch,k+iypatch) - dtsdy * Bzapr(j,k-1)
  END do
  DEALLOCATE(Bzapr)
END if

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

else if(.not. (l_moving_window .and. l_elaser_out_plane)) then
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

else if(.not. (l_moving_window .and. l_elaser_out_plane)) then
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

if(.not. (l_moving_window .and. l_elaser_out_plane)) then
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
end if

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
   subroutine project_jxjy(jxjyfin,jxjygros,rap)
   ! Routine de projection des J d'une grille fine sur une grille grossiere.
   ! Soit nx*ny la taille de la grille grossiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive jxjygros(0:nx+1,0:ny+1,2) et jxjyfin(0:rap*nx+1,0:rap*ny+1,2).
   real(kind=8), DIMENSION(0:,0:,:) :: jxjyfin,jxjygros
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxgros, nygros, j, k, jg, kg
   real(kind=8) :: w, invrapvol

      invrapvol = 1./rap!**2

      nxfin = SIZE(jxjyfin,1)-2
      nyfin = SIZE(jxjyfin,2)-2
      nxgros = SIZE(jxjygros,1)-2
      nygros = SIZE(jxjygros,2)-2

      IF(nxgros*rap/=nxfin .OR. nygros*rap/=nyfin) then
        WRITE(0,*) 'Error in project_jxjy: rap does not match grid sizes.'
        stop
      END if

      do k = 0, nyfin
        kg = k/rap
        w = REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          jxjygros(jg,kg,  1) = jxjygros(jg,kg,  1) + (1.-w)*jxjyfin(j,k,1)*invrapvol
          jxjygros(jg,kg+1,1) = jxjygros(jg,kg+1,1) +     w *jxjyfin(j,k,1)*invrapvol
        end do
      end do

      do k = 0, nyfin-1
        kg = k/rap
        do j = 0, nxfin
          jg = j/rap
          w = REAL(MOD(j,rap))/rap
          jxjygros(jg,  kg,2) = jxjygros(jg  ,kg,2) + (1.-w)*jxjyfin(j,k,2)*invrapvol
          jxjygros(jg+1,kg,2) = jxjygros(jg+1,kg,2) +     w *jxjyfin(j,k,2)*invrapvol
        end do
      end do

   end subroutine project_jxjy

   subroutine interpol_jxjy(jxjyfin,jxjygros,rap)
   ! Routine d'interpolation-soustraction des J d'une grille grossiere sur une grille fine.
   ! Soit nx*ny la taille de la grille grossiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive jxjygros(0:nx+1,0:ny+1,2) et jxjyfin(0:rap*nx+1,0:rap*ny+1,2).
   real(kind=8), DIMENSION(0:,0:,:) :: jxjyfin,jxjygros
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxgros, nygros, j, k, jg, kg
   real(kind=8) :: w
   logical :: l_onfinegridfirst = .true.
   real(kind=8), DIMENSION(:,:,:), allocatable :: jxjyfin_new

      nxfin = SIZE(jxjyfin,1)-2
      nyfin = SIZE(jxjyfin,2)-2
      nxgros = SIZE(jxjygros,1)-2
      nygros = SIZE(jxjygros,2)-2

      IF(nxgros*rap/=nxfin .OR. nygros*rap/=nyfin) then
        WRITE(0,*) 'Error in interpol_jxjy: rap does not match grid sizes.'
!        WRITE(0,*) nxgros,rap,nxfin,nygros,rap,nyfin
        stop
      END if

      if(.not.l_onfinegridfirst) then
        do k = 0, nyfin
          kg = k/rap
          w = REAL(MOD(k,rap))/rap
          do j = 0, nxfin-1
            jg = j/rap
            jxjyfin(j,k,1) = jxjyfin(j,k,1)-((1.-w)*jxjygros(jg,kg,1)+w*jxjygros(jg,kg+1,1))/rap
          end do
        end do

        do k = 0, nyfin-1
          kg = k/rap
          do j = 0, nxfin
            jg = j/rap
            w = REAL(MOD(j,rap))/rap
            jxjyfin(j,k,2) = jxjyfin(j,k,2)-((1.-w)*jxjygros(jg,kg,2)+w*jxjygros(jg+1,kg,2))/rap
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

   subroutine project_ex(exfin,exgros,rap)
   ! Routine de projection des J d'une grille fine sur une grille grossiere.
   ! Soit nx*ny la taille de la grille grossiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive exgros(0:nx+1,0:ny+1,2) et exfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: exfin,exgros
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxgros, nygros, j, k, jg, kg
   real(kind=8) :: w, invrapvol

      invrapvol = 1./rap**2

      nxfin = SIZE(exfin,1)-2
      nyfin = SIZE(exfin,2)-2
      nxgros = SIZE(exgros,1)-2
      nygros = SIZE(exgros,2)-2

      IF(nxgros*rap/=nxfin .OR. nygros*rap/=nyfin) then
        WRITE(0,*) 'Error in project_ex: rap does not match grid sizes.'
        stop
      END if

      do k = 0, nyfin
        kg = k/rap
        w = REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          exgros(jg,kg  ) = exgros(jg,kg  ) + (1.-w)*exfin(j,k)*invrapvol
          exgros(jg,kg+1) = exgros(jg,kg+1) +     w *exfin(j,k)*invrapvol
        end do
      end do

   end subroutine project_ex

   subroutine project_ey(eyfin,eygros,rap)
   ! Routine de projection des J d'une grille fine sur une grille grossiere.
   ! Soit nx*ny la taille de la grille grossiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive eygros(0:nx+1,0:ny+1,2) et eyfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: eyfin,eygros
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxgros, nygros, j, k, jg, kg
   real(kind=8) :: w, invrapvol

      invrapvol = 1./rap**2

      nxfin = SIZE(eyfin,1)-2
      nyfin = SIZE(eyfin,2)-2
      nxgros = SIZE(eygros,1)-2
      nygros = SIZE(eygros,2)-2

      IF(nxgros*rap/=nxfin .OR. nygros*rap/=nyfin) then
        WRITE(0,*) 'Error in project_ey: rap does not match grid sizes.'
        stop
      END if

      do k = 0, nyfin-1
        kg = k/rap
        do j = 0, nxfin
          jg = j/rap
          w = REAL(MOD(j,rap))/rap
          eygros(jg,  kg) = eygros(jg  ,kg) + (1.-w)*eyfin(j,k)*invrapvol
          eygros(jg+1,kg) = eygros(jg+1,kg) +     w *eyfin(j,k)*invrapvol
        end do
      end do

   end subroutine project_ey

   subroutine project_bz(bzfin,bzgros,rap)
   ! Routine de projection de Bz d'une grille fine sur une grille grossiere.
   ! Soit nx*ny la taille de la grille grossiere en nombre de mailles.
   ! On suppose que la maille fine a un nombre de mailles rap*nx*rap*ny.
   ! On passe une maille de plus pour chaque tableau a la limite superieure
   ! pour chaque dimension de facon a pouvoir utiliser la fonction modulo
   ! pour calculer les poids. Les tailles des tableaux passes sont
   ! en definitive bzgros(0:nx+1,0:ny+1) et bzfin(0:rap*nx+1,0:rap*ny+1).
   real(kind=8), DIMENSION(0:,0:) :: bzfin,bzgros
   INTEGER :: rap

   INTEGER :: nxfin, nyfin, nxgros, nygros, j, k, jg, kg
   real(kind=8) :: wx, wy, invrapvol,q

      invrapvol = 1./rap**2

      nxfin = SIZE(bzfin,1)-2
      nyfin = SIZE(bzfin,2)-2
      nxgros = SIZE(bzgros,1)-2
      nygros = SIZE(bzgros,2)-2

      IF(nxgros*rap/=nxfin .OR. nygros*rap/=nyfin) then
        WRITE(0,*) 'Error in project_bz: rap does not match grid sizes.'
        stop
      END if

      do k = 0, nyfin-1
        kg = k/rap
        wy = REAL(MOD(k,rap))/rap
        do j = 0, nxfin-1
          jg = j/rap
          wx = REAL(MOD(j,rap))/rap
          q = bzfin(j,k)*invrapvol
          bzgros(jg,  kg  ) = bzgros(jg,  kg  ) + (1.-wx) * (1.-wy) * q
          bzgros(jg+1,kg  ) = bzgros(jg+1,kg  ) +     wx  * (1.-wy) * q
          bzgros(jg  ,kg+1) = bzgros(jg,  kg+1) + (1.-wx) *     wy  * q
          bzgros(jg+1,kg+1) = bzgros(jg+1,kg+1) +     wx  *     wy  * q
        end do
      end do

   end subroutine project_bz

   subroutine interpol(ffin,fgros,rap)
   real(kind=8), DIMENSION(0:,0:) :: ffin,fgros
   INTEGER :: rap
   INTEGER :: nxfin, nyfin, nxgros, nygros, jf, kf, jg, kg
   real(kind=8) :: wj, wk

      nxfin = SIZE(ffin,1)-2
      nyfin = SIZE(ffin,2)-2
      nxgros = SIZE(fgros,1)-2
      nygros = SIZE(fgros,2)-2

      IF(nxgros*rap/=nxfin .OR. nygros*rap/=nyfin) then
        WRITE(0,*) 'Error in interpol: rap does not match grid sizes.'
        stop
      END if

      do kf = 0, nyfin
        kg = kf/rap
        wk = REAL(MOD(kf,rap))/rap
        do jf = 0, nxfin
          jg = jf/rap
          wj = REAL(MOD(jf,rap))/rap
          ffin(jf,kf) = ffin(jf,kf) + (1.-wj)*(1.-wk)*fgros(jg,  kg)   &
                                    +     wj *(1.-wk)*fgros(jg+1,kg)   &
                                    + (1.-wj)*    wk *fgros(jg,  kg+1) &
                                    +     wj *    wk *fgros(jg+1,kg+1)
        end do
      end do


   END subroutine interpol

end module mod_field

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
subroutine init_fields(f,nx, ny, nbndx, nbndy, dtm, dx, dy, xmin, ymin, rap, xlb, ylb, xrb, yrb)
use mod_bnd
use mod_field, only:EM2D_FIELDtype, l_copyfields, l_elaser_out_plane
implicit none

TYPE(EM2D_FIELDtype) :: f
INTEGER(ISZ), INTENT(IN) :: nx, ny, rap
INTEGER(ISZ), INTENT(IN) :: nbndx, nbndy, xlb, ylb, xrb, yrb
REAL(kind=8), INTENT(IN) :: dtm, dx, dy, xmin, ymin
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
f%dxi = 1./dx
f%dyi = 1./dy
f%npulse=300
f%ntemp = 2*max(nx,ny)+4
if (l_elaser_out_plane) then
  f%cst1 = (1.-dtm/dy)/(1.+dtm/dy)
  f%cst2 =  2.*dtm/dy /(1.+dtm/dy)
else
  f%cst1 = (1.-dtm/dx)/(1.+dtm/dx)
  f%cst2 =  2.*dtm/dx /(1.+dtm/dx)
end if

f%js = 1 ! position of the source

  IF(l_copyfields) then
    f%nxcopy = f%nx
    f%nycopy = f%ny
  ELSE
    f%nxcopy = 0
    f%nycopy = 0
  END if
!  call EM2D_FIELDtypeallot(f)

call create_bnd(f%bndexeybz, nx, ny, nbndx=10, nbndy=10, dt=dtm, dx=dx, dy=dy, xbnd=xlb, ybnd=ylb)
call create_bnd(f%bndbxbyez, nx, ny, nbndx=10, nbndy=10, dt=dtm, dx=dx, dy=dy, xbnd=xlb, ybnd=ylb)
  
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

