!     Last change:  JLV   6 Aug 2001    4:28 pm
#include "top.h"
module multigridrz
USE constant
TYPE conductor_type
  REAL(8) :: voltage
  REAL(8), POINTER, DIMENSION(:) :: phi0xm,phi0xp,phi0ym,phi0yp
  REAL(8), POINTER, DIMENSION(:) :: cf0, cfxp, cfxm, cfyp, cfym, cfrhs
  REAL(8), POINTER, DIMENSION(:) :: rphi0xm,rphi0xp,rphi0ym,rphi0yp
  REAL(8), POINTER, DIMENSION(:) :: rcf0, rcfxp, rcfxm, rcfyp, rcfym, rcfrhs
  REAL(8), POINTER, DIMENSION(:) :: dxm,dxp,dym,dyp
  INTEGER(ISZ), POINTER, DIMENSION(:) :: jj, kk, jcond, kcond
  LOGICAL(ISZ), POINTER :: docalc(:)
  INTEGER(ISZ) :: nbbnd,nbbndred,ncond
  TYPE(conductor_type), POINTER :: next,prev
end TYPE conductor_type

TYPE bndptr
  LOGICAL(ISZ), POINTER:: v(:,:)    ! vacuum points
  REAL(8) :: dx
  REAL(8) :: dy
  INTEGER(ISZ) :: nb_conductors
  TYPE(conductor_type), POINTER :: cnd,first
END TYPE bndptr

TYPE(bndptr), pointer :: bndy(:)
LOGICAL(ISZ) :: bndy_allocated=.false.

REAL(8), parameter :: coefrelax=1._8, coefrelaxbnd=1._8

INTEGER(ISZ), parameter :: dirichlet=0, neumann=1, periodic=2
INTEGER(ISZ) :: ixlbnd=dirichlet, ixrbnd=dirichlet, izlbnd=dirichlet, izrbnd=dirichlet

INTEGER(ISZ), parameter :: egun=0, ecb=1
INTEGER(ISZ) :: bnd_method=egun

contains

subroutine init_bnd(nlevels,nx,ny)
implicit none
INTEGER(ISZ), INTENT(IN) :: nlevels, nx, ny

INTEGER(ISZ) :: i, nx2, ny2

  allocate(bndy(nlevels))
  bndy_allocated=.true.

  nx2=nx; ny2=ny
  do i = nlevels, 1, -1
    NULLIFY(bndy(i)%first)
    bndy(i)%nb_conductors = 0
    bndy(i)%dx=1._8; bndy(i)%dy=1._8 
    ALLOCATE(bndy(i)%v(nx+1,ny+1))
    bndy(i)%v(:,:)=.true.
    nx2=nx2/2; ny2=ny2/2
  end do

end subroutine init_bnd

subroutine init_bnd_sublevel(bndl,nbnd,ncond)
implicit none
TYPE(bndptr), INTENT(IN OUT) :: bndl
INTEGER(ISZ), INTENT(IN) :: nbnd, ncond
TYPE(conductor_type), POINTER :: cnd

ALLOCATE(cnd)
IF(.not.associated(bndl%first)) then
  bndl%first => cnd
  bndl%cnd => cnd
  bndl%nb_conductors = 1
else
  cnd%prev => bndl%cnd
  bndl%cnd%next => cnd
  bndl%cnd => cnd
  bndl%nb_conductors = bndl%nb_conductors + 1
END if

bndl%cnd%nbbnd = nbnd
bndl%cnd%ncond = ncond
IF(nbnd>0) then
  ALLOCATE(bndl%cnd%phi0xm(nbnd),bndl%cnd%phi0xp(nbnd),bndl%cnd%phi0ym(nbnd),bndl%cnd%phi0yp(nbnd))
  ALLOCATE(bndl%cnd%cf0(nbnd),bndl%cnd%cfxp(nbnd),bndl%cnd%cfxm(nbnd),bndl%cnd%cfyp(nbnd),bndl%cnd%cfym(nbnd),bndl%cnd%cfrhs(nbnd))
  ALLOCATE(bndl%cnd%rphi0xm(nbnd),bndl%cnd%rphi0xp(nbnd),bndl%cnd%rphi0ym(nbnd),bndl%cnd%rphi0yp(nbnd))
  ALLOCATE(bndl%cnd%rcf0(nbnd),bndl%cnd%rcfxp(nbnd),bndl%cnd%rcfxm(nbnd), &
           bndl%cnd%rcfyp(nbnd),bndl%cnd%rcfym(nbnd),bndl%cnd%rcfrhs(nbnd))
  ALLOCATE(bndl%cnd%dxm(nbnd),bndl%cnd%dxp(nbnd),bndl%cnd%dym(nbnd),bndl%cnd%dyp(nbnd))
  ALLOCATE(bndl%cnd%jj(nbnd),bndl%cnd%kk(nbnd),bndl%cnd%docalc(nbnd))
  bndl%cnd%docalc=.false.
END if
IF(ncond>0) then
  ALLOCATE(bndl%cnd%jcond(ncond),bndl%cnd%kcond(ncond))
END if

end subroutine init_bnd_sublevel

function expandwguard(f)
implicit none

REAL(8), DIMENSION(0:,0:), INTENT(IN) :: f
REAL(8), DIMENSION(0:2*(SIZE(f,1)-2)-1+1,0:2*(SIZE(f,2)-2)-1+1) :: expandwguard

INTEGER(ISZ) :: j, l, nxe, nye

nxe = 2*(SIZE(f,1)-2)-2
nye = 2*(SIZE(f,2)-2)-2

do l = 1, nye+1, 2
  do j = 1, nxe+1, 2
    expandwguard(j,l) = f(j/2+1,l/2+1)
  end do
end do
do l = 1, nye+1, 2
  do j = 2, nxe, 2
    expandwguard(j,l) = 0.5*(f(j/2,l/2+1)+f(j/2+1,l/2+1))
  end do
end do
do l = 2, nye, 2
  do j = 1, nxe+1, 2
    expandwguard(j,l) = 0.5*(f(j/2+1,l/2)+f(j/2+1,l/2+1))
  end do
end do
do l = 2, nye, 2
  do j = 2, nxe, 2
    expandwguard(j,l) = 0.25*(f(j/2,l/2)+f(j/2,l/2+1)+f(j/2+1,l/2)+f(j/2+1,l/2+1))
  end do
end do

expandwguard(0,:)=0._8
expandwguard(nxe+2,:)=0._8
expandwguard(:,0)=0._8
expandwguard(:,nye+2)=0._8

return
end function expandwguard

subroutine interp(unew, uold, xminold, xmaxold, yminold, ymaxold, xminnew, xmaxnew, yminnew, ymaxnew)
implicit none
REAL(8), DIMENSION(:,:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(:,:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, yminold, ymaxold, xminnew, xmaxnew, yminnew, ymaxnew

INTEGER(ISZ) :: nxnew, nynew, nxold, nyold
INTEGER(ISZ) :: jold, kold, jnew, knew
REAL(8) :: x, y, dxold, dyold, dxnew, dynew, ddx, ddy, rap

nxnew = SIZE(unew,1) - 1
nynew = SIZE(unew,2) - 1
nxold = SIZE(uold,1) - 1
nyold = SIZE(uold,2) - 1

dxnew = (xmaxnew-xminnew) / nxnew
dynew = (ymaxnew-yminnew) / nynew
dxold = (xmaxold-xminold) / nxold
dyold = (ymaxold-yminold) / nyold

rap = dxold*dyold / (dxnew*dynew)

do kold = 1, nyold+1
  y = yminold + (kold-1)*dyold
  knew = MIN(1 + INT((y-yminnew) / dynew), nynew)
  ddy = (y-yminnew)/dynew-real(knew-1)
  do jold = 1, nxold+1
    x = xminold + (jold-1)*dxold
    jnew = MIN(1 + INT((x-xminnew) / dxnew), nxnew)
    ddx = (x-xminnew)/dxnew-real(jnew-1)
    unew(jnew,knew)     = unew(jnew,knew)     + uold(jold,kold) * rap * (1.-ddx) * (1.-ddy)
    unew(jnew+1,knew)   = unew(jnew+1,knew)   + uold(jold,kold) * rap * ddx * (1.-ddy)
    unew(jnew,knew+1)   = unew(jnew,knew+1)   + uold(jold,kold) * rap * (1.-ddx) * ddy
    unew(jnew+1,knew+1) = unew(jnew+1,knew+1) + uold(jold,kold) * rap * ddx * ddy
  end do
end do

return
END subroutine interp

subroutine interp_bndwguard(unew, uold, xminold, xmaxold, yminold, ymaxold, xminnew, xmaxnew, yminnew, ymaxnew)
implicit none
REAL(8), DIMENSION(0:,0:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, yminold, ymaxold, xminnew, xmaxnew, yminnew, ymaxnew

INTEGER(ISZ) :: nxnew, nynew, nxold, nyold
INTEGER(ISZ) :: jold, kold, jnew, knew
REAL(8) :: x, y, dxold, dyold, dxnew, dynew, ddx, ddy

nxnew = SIZE(unew,1)-2 - 1
nynew = SIZE(unew,2)-2 - 1
nxold = SIZE(uold,1)-2 - 1
nyold = SIZE(uold,2)-2 - 1

dxnew = (xmaxnew-xminnew) / nxnew
dynew = (ymaxnew-yminnew) / nynew
dxold = (xmaxold-xminold) / nxold
dyold = (ymaxold-yminold) / nyold

knew = 1
kold = 1
  do jnew = 1, nxnew+1
    x = xminnew+(jnew-1)*dxnew
    jold = MIN(1 + INT((x-xminold) / dxold), nxold)
    ddx = (x-xminold)/dxold-REAL(jold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddx) &
                    + uold(jold+1,kold)   * ddx
  end do
knew = nynew+1
kold = nyold+1
  do jnew = 1, nxnew+1
    x = xminnew+(jnew-1)*dxnew
    jold = MIN(1 + INT((x-xminold) / dxold), nxold)
    ddx = (x-xminold)/dxold-REAL(jold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddx) &
                    + uold(jold+1,kold)   * ddx
  end do
jnew = 1
jold = 1
  do knew = 2, nynew
    y = yminnew+(knew-1)*dynew
    kold = MIN(1 + INT((y-yminold) / dyold), nyold)
    ddy = (y-yminold)/dyold-REAL(kold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddy) &
                    + uold(jold,kold+1)   * ddy
  end do
jnew = nxnew+1
jold = nxold+1
  do knew = 2, nynew
    y = yminnew+(knew-1)*dynew
    kold = MIN(1 + INT((y-yminold) / dyold), nyold)
    ddy = (y-yminold)/dyold-REAL(kold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddy) &
                    + uold(jold,kold+1)   * ddy
  end do

return
END subroutine interp_bndwguard

function restrict(f)
implicit none

REAL(8), DIMENSION(1:,1:), INTENT(IN) :: f
REAL(8), DIMENSION(SIZE(f,1)/2+1, SIZE(f,2)/2+1) :: restrict

INTEGER(ISZ) :: nx, ny, nxi, nyi, j, l, jj, ll

nxi = SIZE(f,1)-1
nyi = SIZE(f,2)-1
nx = nxi/2
ny = nyi/2

do l = 2, ny
  do j = 2, nx
    jj = 2*j-1
    ll = 2*l-1
    restrict(j,l) = 0.25*f(jj,ll) &
                  + 0.125*(f(jj-1,ll)+f(jj,ll-1)+f(jj+1,ll)+f(jj,ll+1))  &
                  + 0.0625*(f(jj-1,ll-1)+f(jj+1,ll-1)+f(jj+1,ll+1)+f(jj-1,ll+1))
  end do
end do
restrict(1:nx+1, 1)    = f(1:nxi+1:2, 1)
restrict(1:nx+1, ny+1) = f(1:nxi+1:2, nyi+1)
restrict(1,      2:ny) = f(1,         3:nyi-1:2)
restrict(nx+1,   2:ny) = f(nxi+1,     3:nyi-1:2)

return
end function restrict

subroutine relaxbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,nc,voltfact)
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc
REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nx+2,0:ny+2)
REAL(8), INTENT(IN) :: rhs(1:,1:)!rhs(nx+1,ny+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic
REAL(8) :: dt, dt0, r
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs

! define CFL
dt = coefrelax/(2._8/dr**2+2._8/dz**2)
dt0 = coefrelax/(4._8/dr**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cf0 = 1._8-2._8*dt/dr**2-2._8*cfz
cfrhs = -dt
do j = 2, nr+1
  cfrp(j) = dt * (1._8+0.5_8/REAL(j-1,8)) / dr**2
  cfrm(j) = dt * (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do i = 1, bnd%cnd%ncond
    f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = voltfact*bnd%cnd%voltage
  end do
end do

do i = 1, nc

lsw = 1
do redblack = 1, 2
  jsw = lsw
  do ic = 1, bnd%nb_conductors
    IF(ic==1) then
      bnd%cnd => bnd%first
    else
      bnd%cnd => bnd%cnd%next
    END if

    IF(redblack==1) THEN !red
      iil=1
      iiu=bnd%cnd%nbbndred
    else !black
      iil=bnd%cnd%nbbndred+1
      iiu=bnd%cnd%nbbnd
    ENDif
    do ii = iil, iiu
      j = bnd%cnd%jj(ii)
      l = bnd%cnd%kk(ii)
      IF(j==1) then
!        f(j,l) =  &
        IF(bnd%cnd%docalc(ii)) &
        f(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l) &
               + bnd%cnd%cfyp(ii)*f(j,l+1)+bnd%cnd%cfym(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0ym(ii)+bnd%cnd%phi0yp(ii)) &
               + bnd%cnd%cfrhs(ii)*rhs(j,l)
      else
!        f(j,l) =  &
        IF(bnd%cnd%docalc(ii)) &
        f(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfyp(ii)*f(j,l+1)+bnd%cnd%cfym(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0ym(ii)+bnd%cnd%phi0yp(ii)) &
               + bnd%cnd%cfrhs(ii)*rhs(j,l)
      END if
    ENDDO
  END do
  do l = 1, nz+1
    IF(jsw==2) then! origin
      j = 1
!      f(j,l) =  &
      IF(bnd%v(j,l)) &
      f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                 + 4._8*dt0*f(j+1,l)/dr**2   &
                                 + (dt0/dz**2)*(f(j,l+1)+f(j,l-1)) &
                                 - dt0*rhs(j,l)
    END if
    do j = jsw+1, nr+1, 2
!        f(j,l) =  &
      IF(bnd%v(j,l)) &
        f(j,l) = cf0 * f(j,l) &
                                 + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 + cfrhs*rhs(j,l)
    end do
    jsw = 3-jsw
  end do
  lsw = 3-lsw
END do !redblack=1, 2

END do !i=1, nc

return
END subroutine relaxbndrzwguard

function residbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,voltfact)
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: f(0:,0:)!f(0:nx+2,0:ny+2)
REAL(8), INTENT(IN OUT) :: rhs(:,:)!rhs(nx+1,ny+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dr, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2) :: residbndrzwguard

INTEGER(ISZ) :: i, j, l, ii, jsw, ksw, redblack, ic
REAL(8) :: r
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz

cfz = -1._8 / dz**2
cf0 = 2._8/dr**2-2._8*cfz
do j = 2, nr+1
  cfrp(j) = - (1._8+0.5_8/REAL(j-1,8)) / dr**2
  cfrm(j) = - (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

residbndrzwguard = 0._8

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
    do i = 1, bnd%cnd%ncond
      rhs(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = 0._8
  end do
END do

do l = 1, nz+1
  j = 1
  IF(bnd%v(j,l)) &
  residbndrzwguard(j,l) = (cf0+2._8/dr**2) * f(j,l) - 4._8*f(j+1,l)/dr**2   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 + rhs(j,l)
  do j = 2, nr+1
     IF(bnd%v(j,l)) &
       residbndrzwguard(j,l) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                                      + cfz*(f(j,l+1)+f(j,l-1)) &
                                      + rhs(j,l)
  end do
end do

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do ii = 1, bnd%cnd%nbbnd
    j = bnd%cnd%jj(ii)
    l = bnd%cnd%kk(ii)
    IF(j==1) then
      IF(bnd%cnd%docalc(ii)) &
      residbndrzwguard(j,l) = bnd%cnd%rcf0(ii)*f(j,l) &
                        + bnd%cnd%rcfxp(ii)*f(j+1,l) &
                        + bnd%cnd%rcfyp(ii)*f(j,l+1)+bnd%cnd%rcfym(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%rphi0xp(ii) &
                        + bnd%cnd%rphi0yp(ii)+bnd%cnd%rphi0ym(ii)) &
                        + bnd%cnd%rcfrhs(ii)*rhs(j,l)
    else
      IF(bnd%cnd%docalc(ii)) &
      residbndrzwguard(j,l) = bnd%cnd%rcf0(ii)*f(j,l) &
                        + bnd%cnd%rcfxp(ii)*f(j+1,l)+bnd%cnd%rcfxm(ii)*f(j-1,l) &
                        + bnd%cnd%rcfyp(ii)*f(j,l+1)+bnd%cnd%rcfym(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%rphi0xp(ii)+bnd%cnd%rphi0xm(ii) &
                        + bnd%cnd%rphi0yp(ii)+bnd%cnd%rphi0ym(ii)) &
                        + bnd%cnd%rcfrhs(ii)*rhs(j,l)
    END if
  ENDDO
END do

return
end function residbndrzwguard

RECURSIVE subroutine mgbndrzwguard(j, u, rhs, bnd, nr, nz, dr, dz, npre, npost, ncycle, sub)
implicit none

INTEGER(ISZ), INTENT(IN) :: j, nr, nz, npre, npost, ncycle
REAL(8), DIMENSION(:,:), INTENT(IN OUT) :: u
REAL(8), DIMENSION(:,:), INTENT(IN OUT) :: rhs
REAL(8) :: dr, dz
TYPE(bndptr) :: bnd(:)
LOGICAL(ISZ), INTENT(IN) :: sub

REAL(8), DIMENSION((SIZE(u,1)-2+1)/2,(SIZE(u,2)-2+1)/2) :: res
REAL(8), DIMENSION(0:(SIZE(u,1)-2+1)/2+1,0:(SIZE(u,2)-2+1)/2+1) :: v
INTEGER(ISZ) :: i

IF(j<=3) then
  call updateguardcellsrz(f=u)
  IF(sub) then
    call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=0._8)
  else
    call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=1._8)
  END if
else
  call updateguardcellsrz(f=u)
  IF(sub) then
    call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=0._8)
  else
    call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=1._8)
  END if
  IF(sub) then
    res = restrict(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=0._8))
  else
    res = restrict(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=1._8))
  END if
  call apply_voltage(res,bnd(j-1),0._8)
  v = 0.0_8
  do i = 1, ncycle  !(1=V cycles, 2=W cycle)
    call mgbndrzwguard(j=j-1, u=v, rhs=res, bnd=bnd(1:j-1), nr=nr/2, nz=nz/2, dr=2._8*dr, dz=2._8*dz, npre=npre, npost=npost, &
                   ncycle=ncycle, sub=.true.)
  end do
  call apply_voltagewguard(v,bnd(j-1),0._8)
  u = u + expandwguard(v)
  IF(sub) then
    call apply_voltagewguard(u,bnd(j),0._8)
  else
    call apply_voltagewguard(u,bnd(j),1._8)
  END if
  call updateguardcellsrz(f=u)
  IF(sub) then
    call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npost,voltfact=0._8)
  else
    call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npost,voltfact=1._8)
  END if
END if

return
end subroutine mgbndrzwguard

subroutine updateguardcellsrz(f)
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)

INTEGER(ISZ) :: ixmax, izmax

ixmax=SIZE(f,1)
izmax=SIZE(f,2)

select case (ixrbnd)
    case (dirichlet)
      f(ixmax,:) = -f(ixmax-2,:)
    case (neumann)
      f(ixmax,:) = f(ixmax-2,:)
    case default
end select
select case (izlbnd)
    case (dirichlet)
      f(:,1) = -f(:,3)
    case (neumann)
      f(:,1) = f(:,3)
    case (periodic)
      f(:,1) = f(:,izmax-1)
    case default
end select
select case (izrbnd)
    case (dirichlet)
      f(:,izmax) = -f(:,izmax-2)
    case (neumann)
      f(:,izmax) = f(:,izmax-2)
    case (periodic)
      f(:,izmax) = f(:,2)
    case default
end select

end subroutine updateguardcellsrz

subroutine apply_voltage(f,bnd,coef_voltage)
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: coef_voltage

INTEGER(ISZ) :: ic, i

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
    do i = 1, bnd%cnd%ncond
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = coef_voltage*bnd%cnd%voltage
  end do
END do

return
end subroutine apply_voltage

subroutine apply_voltagewguard(f,bnd,coef_voltage)
implicit none
REAL(8),INTENT(IN OUT) :: f(0:,0:)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: coef_voltage

INTEGER(ISZ) :: ic, i

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
    do i = 1, bnd%cnd%ncond
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = coef_voltage*bnd%cnd%voltage
  end do
END do

return
end subroutine apply_voltagewguard

subroutine slvfldrz_bnd(u,rhoinit,bnd,nrp,nzp,length_r,length_z,accuracy,nc,npre,npost,ncycle)
implicit none

INTEGER(ISZ), INTENT(IN) :: nrp, nzp, nc, npre, npost, ncycle
REAL(8), INTENT(IN OUT) :: u(:,:)!u(0:(2**nrp)+2,0:(2**nzp)+2)
REAL(8), INTENT(IN) :: rhoinit(:,:)!rhoinit((2**nrp)+1,(2**nzp)+1)
REAL(8), INTENT(IN) :: length_r, length_z, accuracy
TYPE(bndptr), DIMENSION(:) :: bnd

REAL(8), parameter :: pi = 3.14159265358979_8

INTEGER(ISZ) :: i, j, l, np, nrc, nzc
REAL(8), allocatable, DIMENSION(:,:) :: f, fold
REAL(8) :: dr, dz
TYPE rhoptr
  REAL(8), POINTER:: a(:,:)
END TYPE rhoptr
TYPE(rhoptr), ALLOCATABLE :: rho(:)
REAL(8) :: r, theta, r1, r2, x1, x2, y1, y2, xc, yc, average_residue, average_residue_init

np = SIZE(bnd)

nrc = 2**(nrp-np+1)
nzc = 2**(nzp-np+1)
dr = length_r/nrc
dz = length_z/nzc

ALLOCATE(f(0:nrc+2,0:nzc+2), rho(np))
ALLOCATE(rho(np)%a(2**nrp+1,2**nzp+1))

rho(np)%a = rhoinit
f = 0._8
call interp(unew=f(1:nrc+1,1:nzc+1),uold=u, &
            xminold=0._8, xmaxold=length_r, yminold=0._8, ymaxold=length_z, &
            xminnew=0._8, xmaxnew=length_r, yminnew=0._8, ymaxnew=length_z)

do i = np-1, 1, -1
  nrc = 2**(nrp-np+i)
  nzc = 2**(nzp-np+i)
  dr = length_r/nrc
  dz = length_z/nzc
  ALLOCATE(rho(i)%a(nrc+1,nzc+1))
  rho(i)%a = 0.
  rho(i)%a = restrict(rho(i+1)%a)
end do

nrc = 2**(nrp-np+1)
nzc = 2**(nzp-np+1)

do i = 2, np
  ALLOCATE(fold(SIZE(f,1),SIZE(f,2)))
  fold = f
  DEALLOCATE(f)
  nrc = 2**(nrp-np+i)
  nzc = 2**(nzp-np+i)
  dr = length_r/nrc
  dz = length_z/nzc
  ALLOCATE(f(0:nrc+2,0:nzc+2))
  f = expandwguard(fold)
  call interp_bndwguard(unew=f,uold=u, &
                  xminold=0._8, xmaxold=length_r, yminold=0._8, ymaxold=length_r, &
                  xminnew=0._8, xmaxnew=length_z, yminnew=0._8, ymaxnew=length_z)
  DEALLOCATE(fold)
  average_residue_init=SUM(ABS(residbndrzwguard(f=RESHAPE((/(0._8*j,j=1,(nrc+3)*(nzc+3))/), (/nrc+3,nzc+3/)) &
                           ,rhs=rho(i)%a,bnd=bnd(i),nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8)))/(nrc*nzc)
  if(average_residue_init==0._8) cycle
  do  j = 1, nc
    call mgbndrzwguard(j=i,u=f,rhs=rho(i)%a,bnd=bnd,nr=nrc,nz=nzc,dr=dr,dz=dz,npre=npre,npost=npost,ncycle=ncycle,sub=.false.)
    average_residue=SUM(ABS(residbndrzwguard(f=f,rhs=rho(i)%a,bnd=bnd(i),nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8)))/(nrc*nzc)
     IF(i==np) then
       IF(average_residue/average_residue_init <= accuracy) exit
     else
       IF(average_residue/average_residue_init <= 1.e-2) exit
     END if
  end do
end do

WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') average_residue/average_residue_init,j
u = f
call updateguardcellsrz(f=u)

do i = np, 1, -1
  DEALLOCATE(rho(i)%a)
end do
DEALLOCATE(f, rho)

return
end subroutine slvfldrz_bnd

END module multigridrz

subroutine slvfld1gridrz_bnd(u0,rho0,nr0,nz0,dr0,dz0,rxbnd,lzbnd,rzbnd,accuracy,ncmax,npre,npost,ncycle)
USE multigridrz
implicit none
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN) :: rho0(nr0+1,nz0+1)
INTEGER(ISZ), INTENT(IN) :: nr0, nz0, rxbnd, lzbnd, rzbnd,ncmax,npre,npost,ncycle
REAL(8), INTENT(IN) :: dr0, dz0,accuracy

INTEGER(ISZ) :: nrp0, nzp0, j, l, nrc, nzc, np, i

REAL(8) :: wtime,timeinit
real(8) :: u(0:nr0+2,0:nz0+2)

nrp0 = INT(LOG(REAL(nr0+1))/LOG(2.))
nzp0 = INT(LOG(REAL(nz0+1))/LOG(2.))

IF(2**nrp0/=nr0.or.2**nzp0/=nz0) THEN
  WRITE(0,*) 'Error in subroutine solvefieldjlv1gridrz_bnd: all dimensions must be a power of two (number of cells!).'
  WRITE(0,*) 'Aborting...'
  call abort()
END if

nrc = nr0
nzc = nz0
np=1
do i = 1, 1000
  nrc = nrc/2
  nzc = nzc/2
  IF(nrc<=1.or.nzc<=1) exit
  np = np + 1
end do
if(.not.bndy_allocated) call init_bnd(np+1,nr0,nz0)

ixrbnd = rxbnd
izlbnd = lzbnd
izrbnd = rzbnd
u(1:nr0+1,:)=u0(1:nr0+1,:)
u(0,:) = 0._8
u(nr0+2,:)=0._8
  call slvfldrz_bnd(u=u,rhoinit=-rho0/eps0,bnd=bndy,nrp=nrp0,nzp=nzp0, &
                  length_r=nr0*dr0,length_z=nz0*dz0, &
                  accuracy=accuracy,nc=ncmax,npre=npre,npost=npost,ncycle=ncycle)
u0(1:nr0+1,:)=u(1:nr0+1,:)

return
end subroutine slvfld1gridrz_bnd

subroutine srfrvoutrz(rofzfunc,volt,zmin,zmax,xcent,rmax,lfill, &
                       xmin,xmax,lshell, &
                       zmmin,zmmax,zbeam,dx,dz,nx,nz, &
                       ix_axis,xmesh)
use PSOR3d
USE multigridrz
implicit none
character(*) rofzfunc
real(kind=8):: volt,zmin,zmax,xcent,rmax
LOGICAL(ISZ):: lfill,lshell
real(kind=8):: xmin,xmax
real(kind=8):: zmmin,zmmax,zbeam,dx,dz
INTEGER(ISZ):: nx,nz,ix_axis
real(kind=8):: xmesh(0:nx)

INTEGER(ISZ) :: i,ii,iii,iv,iiv,nxbnd,nzbnd,nrc,nzc,np,nrp0,nzp0, lfill_tmp,j,l
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz,voltxm,voltxp,voltzm,voltzp

TYPE(conductor_type), POINTER :: cndpnt

nrc = nx
nzc = nz
np=1
do i = 1, 1000
  nrc = nrc/2
  nzc = nzc/2
  IF(nrc<=1.or.nzc<=1) exit
  np = np + 1
end do
IF(.not.bndy_allocated) call init_bnd(np+1,nx,nz)

nrp0 = INT(LOG(REAL(nx+1))/LOG(2.))
nzp0 = INT(LOG(REAL(nz+1))/LOG(2.))

IF(2**nrp0/=nx.or.2**nzp0/=nz) THEN
  WRITE(0,*) 'Error in subroutine bnd_init: all dimensions must be a power of two (number of cells!).'
  WRITE(0,*) 'Aborting...'
  call abort()
END if

nrc = 2*nx
nzc = 2*nz
do i = np+1,1,-1
  nrc = nrc/2
  nzc = nzc/2

  nxbnd=nrc+1
  nzbnd=nzc+1
  bndy(i)%dx = (dx*nx)/nrc
  bndy(i)%dy = (dz*nz)/nzc

  necndbdy=0
  nocndbdy=0
  ncond = 0

  call srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,rmax,lfill, &
                       xmin,xmax,lshell, &
                       zmmin,zmmax,zbeam,bndy(i)%dx,bndy(i)%dy,nrc,nzc, &
                       ix_axis,xmesh,nz/(nzbnd-1))

  call init_bnd_sublevel(bndy(i),necndbdy+nocndbdy,ncond)

  bndy(i)%cnd%nbbndred = necndbdy
  bndy(i)%cnd%nbbnd = necndbdy + nocndbdy

  bndy(i)%cnd%ncond = ncond
  bndy(i)%cnd%voltage = volt
  do ii=1,ncond
    bndy(i)%cnd%jcond(ii) = ixcond(ii)+1
    bndy(i)%cnd%kcond(ii) = izcond(ii)+1
    bndy(i)%v(bndy(i)%cnd%jcond(ii),bndy(i)%cnd%kcond(ii)) = .false.
  end do

  do ii = 1, necndbdy
   bndy(i)%cnd%jj(ii)  = iecndx(ii)+1
   bndy(i)%cnd%kk(ii)  = iecndz(ii)+1
   IF(iecndx(ii)>=0 .and. iecndx(ii)<=nxbnd-1 .and. iecndz(ii)>=0 .and. iecndz(ii)<=nzbnd-1) bndy(i)%cnd%docalc(ii)=.true.
   bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii)) = .false.
   dxm = MIN(1._8,ecdelmx(ii))*bndy(i)%dx
   dxp = MIN(1._8,ecdelpx(ii))*bndy(i)%dx
   dzm = MIN(1._8,ecdelmz(ii))*bndy(i)%dy
   dzp = MIN(1._8,ecdelpz(ii))*bndy(i)%dy
   voltxm=volt
   voltxp=volt
   voltzm=volt
   voltzp=volt
   IF(.not.bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii)) ) then
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=1,cndpnt%nbbndred
         IF(bndy(i)%cnd%jj(ii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(ii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)/=bndy(i)%dx) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dx) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dym(iiv)/=bndy(i)%dy) then
             dzm=cndpnt%dym(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dyp(iiv)/=bndy(i)%dy) then
             dzp=cndpnt%dyp(iiv)
             voltzp=cndpnt%voltage
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(ii)=dxm
   bndy(i)%cnd%dxp(ii)=dxp
   bndy(i)%cnd%dym(ii)=dzm
   bndy(i)%cnd%dyp(ii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dx
       dzz=bndy(i)%dy
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(bndy(i)%cnd%jj(ii)==1) then
     rp = 0.5_8*bndy(i)%dx
     dt = coefrelaxbnd/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(ii) = 4._8*dt/(dxp*dxx)
     bndy(i)%cnd%cfym(ii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(ii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(ii) = -dt
     bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfym(ii)-bndy(i)%cnd%cfyp(ii)
     bndy(i)%cnd%rcfxp(ii) = -4._8/(dxp*dxx)
     bndy(i)%cnd%rcfym(ii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(ii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(ii) = 1._8
     bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxp(ii)-bndy(i)%cnd%rcfym(ii)-bndy(i)%cnd%rcfyp(ii)
   else
     r = (bndy(i)%cnd%jj(ii)-1)*bndy(i)%dx
     rm = r-0.5_8*dxx
     rp = r+0.5_8*dxx
     dt = coefrelaxbnd/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = dt*rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(ii) = dt*rp/(r*dxp*dxx)
     bndy(i)%cnd%cfym(ii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(ii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(ii) = -dt
     bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxm(ii)-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfym(ii)-bndy(i)%cnd%cfyp(ii)
     bndy(i)%cnd%rcfxm(ii) = -1._8*rm/(r*dxm*dxx)
     bndy(i)%cnd%rcfxp(ii) = -1._8*rp/(r*dxp*dxx)
     bndy(i)%cnd%rcfym(ii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(ii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(ii) = 1._8
     bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxm(ii)-bndy(i)%cnd%rcfxp(ii)-bndy(i)%cnd%rcfym(ii)-bndy(i)%cnd%rcfyp(ii)
   END if
     IF(ecdelmx(ii)>=1.) then
       bndy(i)%cnd%phi0xm(ii)=0._8
       bndy(i)%cnd%rphi0xm(ii)=0._8
     else
       bndy(i)%cnd%phi0xm(ii)=bndy(i)%cnd%cfxm(ii)*voltxm
       bndy(i)%cnd%cfxm(ii)=0._8
       bndy(i)%cnd%rphi0xm(ii)=bndy(i)%cnd%rcfxm(ii)*voltxm
       bndy(i)%cnd%rcfxm(ii)=0._8
     END if
     IF(ecdelpx(ii)>=1.) then
       bndy(i)%cnd%phi0xp(ii)=0._8
       bndy(i)%cnd%rphi0xp(ii)=0._8
     else
       bndy(i)%cnd%phi0xp(ii)=bndy(i)%cnd%cfxp(ii)*voltxp
       bndy(i)%cnd%cfxp(ii)=0._8
       bndy(i)%cnd%rphi0xp(ii)=bndy(i)%cnd%rcfxp(ii)*voltxp
       bndy(i)%cnd%rcfxp(ii)=0._8
     END if
     IF(ecdelmz(ii)>=1.) then
       bndy(i)%cnd%phi0ym(ii)=0._8
       bndy(i)%cnd%rphi0ym(ii)=0._8
     else
       bndy(i)%cnd%phi0ym(ii)=bndy(i)%cnd%cfym(ii)*voltzm
       bndy(i)%cnd%cfym(ii)=0._8
       bndy(i)%cnd%rphi0ym(ii)=bndy(i)%cnd%rcfym(ii)*voltzm
       bndy(i)%cnd%rcfym(ii)=0._8
     END if
     IF(ecdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0yp(ii)=0._8
       bndy(i)%cnd%rphi0yp(ii)=0._8
     else
       bndy(i)%cnd%phi0yp(ii)=bndy(i)%cnd%cfyp(ii)*voltzp
       bndy(i)%cnd%cfyp(ii)=0._8
       bndy(i)%cnd%rphi0yp(ii)=bndy(i)%cnd%rcfyp(ii)*voltzp
       bndy(i)%cnd%rcfyp(ii)=0._8
     END if
  end do
  do ii = 1, nocndbdy
   iii=necndbdy+ii
   bndy(i)%cnd%jj(iii)  = iocndx(ii)+1
   bndy(i)%cnd%kk(iii)  = iocndz(ii)+1
   IF(iocndx(ii)>=0 .and. iocndx(ii)<=nxbnd-1 .and. iocndz(ii)>=0 .and. iocndz(ii)<=nzbnd-1) bndy(i)%cnd%docalc(iii)=.true.

   bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) = .false.
   dxm = MIN(1._8,ocdelmx(ii))*bndy(i)%dx
   dxp = MIN(1._8,ocdelpx(ii))*bndy(i)%dx
   dzm = MIN(1._8,ocdelmz(ii))*bndy(i)%dy
   dzp = MIN(1._8,ocdelpz(ii))*bndy(i)%dy
   voltxm=volt
   voltxp=volt
   voltzm=volt
   voltzp=volt
   IF(.not.bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) ) then
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=1,cndpnt%nbbndred
         IF(bndy(i)%cnd%jj(iii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(iii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)/=bndy(i)%dx) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dx) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dym(iiv)/=bndy(i)%dy) then
             dzm=cndpnt%dym(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dyp(iiv)/=bndy(i)%dy) then
             dzp=cndpnt%dyp(iiv)
             voltzp=cndpnt%voltage
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(iii)=dxm
   bndy(i)%cnd%dxp(iii)=dxp
   bndy(i)%cnd%dym(iii)=dzm
   bndy(i)%cnd%dyp(iii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dx
       dzz=bndy(i)%dy
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(bndy(i)%cnd%jj(iii)==1) then
     rp = 0.5_8*bndy(i)%dx
     dt = coefrelaxbnd/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(iii) = 4._8*dt/(dxp*dxx)
     bndy(i)%cnd%cfym(iii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(iii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(iii) = -dt
     bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfym(iii)-bndy(i)%cnd%cfyp(iii)
     bndy(i)%cnd%rcfxp(iii) = -4._8/(dxp*dxx)
     bndy(i)%cnd%rcfym(iii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(iii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(iii) = 1._8
     bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxp(iii)-bndy(i)%cnd%rcfym(iii)-bndy(i)%cnd%rcfyp(iii)
   else
     r = (bndy(i)%cnd%jj(iii)-1)*bndy(i)%dx
     rm = r-0.5_8*dxm
     rp = r+0.5_8*dxp
     dt = coefrelaxbnd/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(iii) = dt*rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(iii) = dt*rp/(r*dxp*dxx)
     bndy(i)%cnd%cfym(iii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(iii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(iii) = -dt
     bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxm(iii)-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfym(iii)-bndy(i)%cnd%cfyp(iii)
     bndy(i)%cnd%rcfxm(iii) = -1._8*rm/(r*dxm*dxx)
     bndy(i)%cnd%rcfxp(iii) = -1._8*rp/(r*dxp*dxx)
     bndy(i)%cnd%rcfym(iii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(iii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(iii) = 1._8
     bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxm(iii)-bndy(i)%cnd%rcfxp(iii)-bndy(i)%cnd%rcfym(iii)-bndy(i)%cnd%rcfyp(iii)
   END if
     IF(ocdelmx(ii)>=1.) then
       bndy(i)%cnd%phi0xm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=0._8
     else
       bndy(i)%cnd%phi0xm(iii)=bndy(i)%cnd%cfxm(iii)*voltxm
       bndy(i)%cnd%cfxm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=bndy(i)%cnd%rcfxm(iii)*voltxm
       bndy(i)%cnd%rcfxm(iii)=0._8
     END if
     IF(ocdelpx(ii)>=1.) then
       bndy(i)%cnd%phi0xp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=0._8
     else
       bndy(i)%cnd%phi0xp(iii)=bndy(i)%cnd%cfxp(iii)*voltxp
       bndy(i)%cnd%cfxp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=bndy(i)%cnd%rcfxp(iii)*voltxp
       bndy(i)%cnd%rcfxp(iii)=0._8
     END if
     IF(ocdelmz(ii)>=1.) then
       bndy(i)%cnd%phi0ym(iii)=0._8
       bndy(i)%cnd%rphi0ym(iii)=0._8
     else
       bndy(i)%cnd%phi0ym(iii)=bndy(i)%cnd%cfym(iii)*voltzm
       bndy(i)%cnd%cfym(iii)=0._8
       bndy(i)%cnd%rphi0ym(iii)=bndy(i)%cnd%rcfym(iii)*voltzm
       bndy(i)%cnd%rcfym(iii)=0._8
     END if
     IF(ocdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0yp(iii)=0._8
       bndy(i)%cnd%rphi0yp(iii)=0._8
     else
       bndy(i)%cnd%phi0yp(iii)=bndy(i)%cnd%cfyp(iii)*voltzp
       bndy(i)%cnd%cfyp(iii)=0._8
       bndy(i)%cnd%rphi0yp(iii)=bndy(i)%cnd%rcfyp(iii)*voltzp
       bndy(i)%cnd%rcfyp(iii)=0._8
     END if
  end do
end do
end subroutine srfrvoutrz

subroutine srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,   &
                        lzend,xmin,xmax,lshell,             &
                        zmmin,zmmax,zbeam,dx,dz,nx,nz,      &
                        ix_axis,xmesh)
use PSOR3d
USE multigridrz
implicit none
character(*) rminofz,rmaxofz
real(kind=8):: volt,zmin,zmax,xcent
LOGICAL(ISZ):: lzend,lshell
real(kind=8):: xmin,xmax
real(kind=8):: zmmin,zmmax,zbeam,dx,dz
INTEGER(ISZ):: nx,nz,ix_axis
real(kind=8):: xmesh(0:nx)

INTEGER(ISZ) :: i,ii,iii,iv,iiv,nxbnd,nzbnd,nrc,nzc,np,nrp0,nzp0, lfill_tmp,j,l
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz,voltxm,voltxp,voltzm,voltzp

TYPE(conductor_type), POINTER :: cndpnt


nrc = nx
nzc = nz
np=1
do i = 1, 1000
  nrc = nrc/2
  nzc = nzc/2
  IF(nrc<=1.or.nzc<=1) exit
  np = np + 1
end do
IF(.not.bndy_allocated) call init_bnd(np+1,nx,nz)

nrp0 = INT(LOG(REAL(nx+1))/LOG(2.))
nzp0 = INT(LOG(REAL(nz+1))/LOG(2.))

IF(2**nrp0/=nx.or.2**nzp0/=nz) THEN
  WRITE(0,*) 'Error in subroutine bnd_init: all dimensions must be a power of two (number of cells!).'
  WRITE(0,*) 'Aborting...'
  call abort()
END if

nrc = 2*nx
nzc = 2*nz
do i = np+1,1,-1
  nrc = nrc/2
  nzc = nzc/2

  nxbnd=nrc+1
  nzbnd=nzc+1
  bndy(i)%dx = (dx*nx)/nrc
  bndy(i)%dy = (dz*nz)/nzc

  necndbdy=0
  nocndbdy=0
  ncond = 0

  call srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,   &
                        lzend,xmin,xmax,lshell,             &
                        zmmin,zmmax,zbeam,bndy(i)%dx,bndy(i)%dy,nrc,nzc, &
                        ix_axis,xmesh)

  call init_bnd_sublevel(bndy(i),necndbdy+nocndbdy,ncond)

  bndy(i)%cnd%nbbndred = necndbdy
  bndy(i)%cnd%nbbnd = necndbdy + nocndbdy

  bndy(i)%cnd%ncond = ncond
  bndy(i)%cnd%voltage = volt
  do ii=1,ncond
    bndy(i)%cnd%jcond(ii) = ixcond(ii)+1
    bndy(i)%cnd%kcond(ii) = izcond(ii)+1
    bndy(i)%v(bndy(i)%cnd%jcond(ii),bndy(i)%cnd%kcond(ii)) = .false.
  end do

  do ii = 1, necndbdy
   bndy(i)%cnd%jj(ii)  = iecndx(ii)+1
   bndy(i)%cnd%kk(ii)  = iecndz(ii)+1
   IF(iecndx(ii)>=0 .and. iecndx(ii)<=nxbnd-1 .and. iecndz(ii)>=0 .and. iecndz(ii)<=nzbnd-1) bndy(i)%cnd%docalc(ii)=.true.
   bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii)) = .false.
   dxm = MIN(1._8,ecdelmx(ii))*bndy(i)%dx
   dxp = MIN(1._8,ecdelpx(ii))*bndy(i)%dx
   dzm = MIN(1._8,ecdelmz(ii))*bndy(i)%dy
   dzp = MIN(1._8,ecdelpz(ii))*bndy(i)%dy
   voltxm=volt
   voltxp=volt
   voltzm=volt
   voltzp=volt
   IF(.not.bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii)) ) then
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=1,cndpnt%nbbndred
         IF(bndy(i)%cnd%jj(ii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(ii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)/=bndy(i)%dx) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dx) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dym(iiv)/=bndy(i)%dy) then
             dzm=cndpnt%dym(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dyp(iiv)/=bndy(i)%dy) then
             dzp=cndpnt%dyp(iiv)
             voltzp=cndpnt%voltage
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(ii)=dxm
   bndy(i)%cnd%dxp(ii)=dxp
   bndy(i)%cnd%dym(ii)=dzm
   bndy(i)%cnd%dyp(ii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dx
       dzz=bndy(i)%dy
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(bndy(i)%cnd%jj(ii)==1) then
     rp = 0.5_8*bndy(i)%dx
     dt = coefrelaxbnd/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(ii) = 4._8*dt/(dxp*dxx)
     bndy(i)%cnd%cfym(ii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(ii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(ii) = -dt
     bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfym(ii)-bndy(i)%cnd%cfyp(ii)
     bndy(i)%cnd%rcfxp(ii) = -4._8/(dxp*dxx)
     bndy(i)%cnd%rcfym(ii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(ii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(ii) = 1._8
     bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxp(ii)-bndy(i)%cnd%rcfym(ii)-bndy(i)%cnd%rcfyp(ii)
   else
     r = (bndy(i)%cnd%jj(ii)-1)*bndy(i)%dx
     rm = r-0.5_8*dxx
     rp = r+0.5_8*dxx
     dt = coefrelaxbnd/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = dt*rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(ii) = dt*rp/(r*dxp*dxx)
     bndy(i)%cnd%cfym(ii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(ii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(ii) = -dt
     bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxm(ii)-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfym(ii)-bndy(i)%cnd%cfyp(ii)
     bndy(i)%cnd%rcfxm(ii) = -1._8*rm/(r*dxm*dxx)
     bndy(i)%cnd%rcfxp(ii) = -1._8*rp/(r*dxp*dxx)
     bndy(i)%cnd%rcfym(ii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(ii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(ii) = 1._8
     bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxm(ii)-bndy(i)%cnd%rcfxp(ii)-bndy(i)%cnd%rcfym(ii)-bndy(i)%cnd%rcfyp(ii)
   END if
     IF(ecdelmx(ii)>=1.) then
       bndy(i)%cnd%phi0xm(ii)=0._8
       bndy(i)%cnd%rphi0xm(ii)=0._8
     else
       bndy(i)%cnd%phi0xm(ii)=bndy(i)%cnd%cfxm(ii)*voltxm
       bndy(i)%cnd%cfxm(ii)=0._8
       bndy(i)%cnd%rphi0xm(ii)=bndy(i)%cnd%rcfxm(ii)*voltxm
       bndy(i)%cnd%rcfxm(ii)=0._8
     END if
     IF(ecdelpx(ii)>=1.) then
       bndy(i)%cnd%phi0xp(ii)=0._8
       bndy(i)%cnd%rphi0xp(ii)=0._8
     else
       bndy(i)%cnd%phi0xp(ii)=bndy(i)%cnd%cfxp(ii)*voltxp
       bndy(i)%cnd%cfxp(ii)=0._8
       bndy(i)%cnd%rphi0xp(ii)=bndy(i)%cnd%rcfxp(ii)*voltxp
       bndy(i)%cnd%rcfxp(ii)=0._8
     END if
     IF(ecdelmz(ii)>=1.) then
       bndy(i)%cnd%phi0ym(ii)=0._8
       bndy(i)%cnd%rphi0ym(ii)=0._8
     else
       bndy(i)%cnd%phi0ym(ii)=bndy(i)%cnd%cfym(ii)*voltzm
       bndy(i)%cnd%cfym(ii)=0._8
       bndy(i)%cnd%rphi0ym(ii)=bndy(i)%cnd%rcfym(ii)*voltzm
       bndy(i)%cnd%rcfym(ii)=0._8
     END if
     IF(ecdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0yp(ii)=0._8
       bndy(i)%cnd%rphi0yp(ii)=0._8
     else
       bndy(i)%cnd%phi0yp(ii)=bndy(i)%cnd%cfyp(ii)*voltzp
       bndy(i)%cnd%cfyp(ii)=0._8
       bndy(i)%cnd%rphi0yp(ii)=bndy(i)%cnd%rcfyp(ii)*voltzp
       bndy(i)%cnd%rcfyp(ii)=0._8
     END if
  end do
  do ii = 1, nocndbdy
   iii=necndbdy+ii
   bndy(i)%cnd%jj(iii)  = iocndx(ii)+1
   bndy(i)%cnd%kk(iii)  = iocndz(ii)+1
   IF(iocndx(ii)>=0 .and. iocndx(ii)<=nxbnd-1 .and. iocndz(ii)>=0 .and. iocndz(ii)<=nzbnd-1) bndy(i)%cnd%docalc(iii)=.true.
   bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) = .false.
   dxm = MIN(1._8,ocdelmx(ii))*bndy(i)%dx
   dxp = MIN(1._8,ocdelpx(ii))*bndy(i)%dx
   dzm = MIN(1._8,ocdelmz(ii))*bndy(i)%dy
   dzp = MIN(1._8,ocdelpz(ii))*bndy(i)%dy
   voltxm=volt
   voltxp=volt
   voltzm=volt
   voltzp=volt
   IF(.not.bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) ) then
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=1,cndpnt%nbbndred
         IF(bndy(i)%cnd%jj(iii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(iii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)/=bndy(i)%dx) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dx) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dym(iiv)/=bndy(i)%dy) then
             dzm=cndpnt%dym(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dyp(iiv)/=bndy(i)%dy) then
             dzp=cndpnt%dyp(iiv)
             voltzp=cndpnt%voltage
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(iii)=dxm
   bndy(i)%cnd%dxp(iii)=dxp
   bndy(i)%cnd%dym(iii)=dzm
   bndy(i)%cnd%dyp(iii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dx
       dzz=bndy(i)%dy
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(bndy(i)%cnd%jj(iii)==1) then
     rp = 0.5_8*bndy(i)%dx
     dt = coefrelaxbnd/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(iii) = 4._8*dt/(dxp*dxx)
     bndy(i)%cnd%cfym(iii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(iii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(iii) = -dt
     bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfym(iii)-bndy(i)%cnd%cfyp(iii)
     bndy(i)%cnd%rcfxp(iii) = -4._8/(dxp*dxx)
     bndy(i)%cnd%rcfym(iii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(iii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(iii) = 1._8
     bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxp(iii)-bndy(i)%cnd%rcfym(iii)-bndy(i)%cnd%rcfyp(iii)
   else
     r = (bndy(i)%cnd%jj(iii)-1)*bndy(i)%dx
     rm = r-0.5_8*dxm
     rp = r+0.5_8*dxp
     dt = coefrelaxbnd/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(iii) = dt*rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(iii) = dt*rp/(r*dxp*dxx)
     bndy(i)%cnd%cfym(iii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfyp(iii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(iii) = -dt
     bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxm(iii)-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfym(iii)-bndy(i)%cnd%cfyp(iii)
     bndy(i)%cnd%rcfxm(iii) = -1._8*rm/(r*dxm*dxx)
     bndy(i)%cnd%rcfxp(iii) = -1._8*rp/(r*dxp*dxx)
     bndy(i)%cnd%rcfym(iii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfyp(iii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(iii) = 1._8
     bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxm(iii)-bndy(i)%cnd%rcfxp(iii)-bndy(i)%cnd%rcfym(iii)-bndy(i)%cnd%rcfyp(iii)
   END if
     IF(ocdelmx(ii)>=1.) then
       bndy(i)%cnd%phi0xm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=0._8
     else
       bndy(i)%cnd%phi0xm(iii)=bndy(i)%cnd%cfxm(iii)*voltxm
       bndy(i)%cnd%cfxm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=bndy(i)%cnd%rcfxm(iii)*voltxm
       bndy(i)%cnd%rcfxm(iii)=0._8
     END if
     IF(ocdelpx(ii)>=1.) then
       bndy(i)%cnd%phi0xp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=0._8
     else
       bndy(i)%cnd%phi0xp(iii)=bndy(i)%cnd%cfxp(iii)*voltxp
       bndy(i)%cnd%cfxp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=bndy(i)%cnd%rcfxp(iii)*voltxp
       bndy(i)%cnd%rcfxp(iii)=0._8
     END if
     IF(ocdelmz(ii)>=1.) then
       bndy(i)%cnd%phi0ym(iii)=0._8
       bndy(i)%cnd%rphi0ym(iii)=0._8
     else
       bndy(i)%cnd%phi0ym(iii)=bndy(i)%cnd%cfym(iii)*voltzm
       bndy(i)%cnd%cfym(iii)=0._8
       bndy(i)%cnd%rphi0ym(iii)=bndy(i)%cnd%rcfym(iii)*voltzm
       bndy(i)%cnd%rcfym(iii)=0._8
     END if
     IF(ocdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0yp(iii)=0._8
       bndy(i)%cnd%rphi0yp(iii)=0._8
     else
       bndy(i)%cnd%phi0yp(iii)=bndy(i)%cnd%cfyp(iii)*voltzp
       bndy(i)%cnd%cfyp(iii)=0._8
       bndy(i)%cnd%rphi0yp(iii)=bndy(i)%cnd%rcfyp(iii)*voltzp
       bndy(i)%cnd%rcfyp(iii)=0._8
     END if
  end do
end do
end subroutine srfrvinoutrz

!     ******************************************************************
!     *
!     *                        SUBROUTINE RHOWEIGHTRZ
!     *
!     ******************************************************************


subroutine rhoweightrz(xp,yp,zp,np,q,rho,nr,nz,dr,dz,zmin)
USE constant
implicit none

REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(0:nr,0:nz), INTENT(INOUT) :: rho
REAL(8), INTENT(IN) :: q, dr, dz, zmin
INTEGER(ISZ), INTENT(IN) :: np, nr, nz

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr)
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp

  invdr = 1._8/dr
  invdz = 1._8/dz

  ! computes divider by cell volumes to get density
  j = 0
  ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
  ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
  ! and Verboncoeur, J. of Comp. Phys.,
  invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
  do j = 1, nr
    invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr * dz)
  end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(zp(i)<zmin) cycle
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr
    zpos = (zp(i)-zmin)*invdz
    jn = INT(rpos)
    ln = INT(zpos)
    ddr = rpos-REAL(jn)
    ddz = zpos-REAL(ln)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    rho(jn, ln)  = rho(jn, ln)  + q * oddr * oddz * invvol(jn)
    rho(jnp,ln)  = rho(jnp,ln)  + q *  ddr * oddz * invvol(jnp)
    rho(jn, lnp) = rho(jn, lnp) + q * oddr *  ddz * invvol(jn)
    rho(jnp,lnp) = rho(jnp,lnp) + q *  ddr *  ddz * invvol(jnp)
  end do

  return
END SUBROUTINE RHOWEIGHTRZ

subroutine rhobndrz(rho,nr,nz)
USE multigridrz
implicit none

REAL(8), dimension(0:nr,0:nz), INTENT(IN OUT) :: rho
INTEGER(ISZ), INTENT(IN) :: nr, nz

INTEGER(ISZ) :: irmax, izmax

  select case (ixrbnd)
    case (dirichlet)
    case (neumann)
      rho(nr,:) = 2._8*rho(nr,:)
    case default
  end select
  select case (izlbnd)
    case (dirichlet)
    case (neumann)
      rho(:,0) = 2._8*rho(:,0)
    case default
  end select
  select case (izrbnd)
    case (dirichlet)
    case (neumann)
      rho(:,nz) = 2._8*rho(:,nz)
    case default
  end select
  
  return
end subroutine rhobndrz

subroutine fieldweightrz(xp,yp,zp,ex,ey,ez,np,phi,nr,nz,dr,dz,zmin)
USE constant
implicit none

REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), INTENT(IN) :: dr, dz, zmin
INTEGER(ISZ), INTENT(IN) :: np, nr, nz

REAL(8) :: invdr, invdz, rpos, zpos, invrpos, ddr, ddz, oddr, oddz, er
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp
REAL(8) :: e(2,0:nr,0:nz)

  invdr = 0.5_8/dr
  invdz = 0.5_8/dz

! compute electric field e from phi
 ! interior
  do l = 0, nz
    do j = 1, nr-1
      e(1,j,l) = invdr * (phi(j-1,l)-phi(j+1,l))
      e(2,j,l) = invdz * (phi(j,l-1)-phi(j,l+1))
    end do
  end do
 ! sides
  j = 0
  do l = 1, nz-1
    e(1,j,l)= 0.
    e(2,j,l)=      invdz * (phi(j,l-1)-phi(j,l+1))
  end do
  j = nr
  do l = 1, nz-1
    e(1,j,l)= 2._8*invdr * (phi(j-1,l)-phi(j,l)  )
    e(2,j,l)=      invdz * (phi(j,l-1)-phi(j,l+1))
  end do
 ! corners
  j=0;l=0
    e(1,j,l)= 0.
    e(2,j,l)= 2._8*invdz * (phi(j,l)  -phi(j,l+1))
  j=nr;l=0
    e(1,j,l)= 2._8*invdr * (phi(j-1,l)-phi(j,l)  )
    e(2,j,l)= 2._8*invdz * (phi(j,l)  -phi(j,l+1))
  j=0;l=nz
    e(1,j,l)= 0.
    e(2,j,l)= 2._8*invdz * (phi(j,l-1)-phi(j,l)  )
  j=nr;l=nz
    e(1,j,l)= 2._8*invdr * (phi(j-1,l)-phi(j,l)  )
    e(2,j,l)= 2._8*invdz * (phi(j,l-1)-phi(j,l)  )


  invdr = 1._8/dr
  invdz = 1._8/dz

  ! make field deposition using CIC weighting
  do i = 1, np
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr
    zpos = (zp(i)-zmin)*invdz
    jn = INT(rpos)
    ln = INT(zpos)
    ddr = rpos-REAL(jn)
    ddz = zpos-REAL(ln)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    er = oddr * oddz * e(1,jn ,ln)  &
       + ddr  * oddz * e(1,jnp,ln)  &
       + oddr * ddz  * e(1,jn ,lnp) &
       + ddr  * ddz  * e(1,jnp,lnp)
    IF(rpos>1.e-10) then
      invrpos=invdr/rpos
      ex(i) = er*xp(i)*invrpos
      ey(i) = er*yp(i)*invrpos
    else
      ex(i) = er
      ey(i) = 0._8
    END if
    ez(i) = oddr * oddz * e(2,jn ,ln)  &
          + ddr  * oddz * e(2,jnp,ln)  &
          + oddr * ddz  * e(2,jn ,lnp) &
          + ddr  * ddz  * e(2,jnp,lnp)
  END do

  return
end subroutine fieldweightrz

subroutine save_bndstructure_rz(filename)
use multigridrz
implicit none
CHARACTER(*) :: filename
INTEGER(ISZ) :: i,ic

  OPEN(10,FILE=filename,STATUS='unknown')
    WRITE(10,*) SIZE(bndy)
    do i = 1, SIZE(bndy)
      WRITE(10,*) SIZE(bndy(i)%v,1),SIZE(bndy(i)%v,2)
      WRITE(10,*) bndy(i)%v
      WRITE(10,*) bndy(i)%dx, bndy(i)%dy, bndy(i)%nb_conductors
      do ic = 1, bndy(i)%nb_conductors
        IF(ic==1) then
          bndy(i)%cnd => bndy(i)%first
        else
          bndy(i)%cnd => bndy(i)%cnd%next
        END if
          WRITE(0,*) i,ic,SIZE(bndy),bndy(i)%cnd%nbbnd,bndy(i)%cnd%ncond, bndy(i)%cnd%nbbndred
          WRITE(10,*) bndy(i)%cnd%nbbnd,bndy(i)%cnd%ncond
          WRITE(10,*) bndy(i)%cnd%nbbndred
          WRITE(10,*) bndy(i)%cnd%voltage
        IF(bndy(i)%cnd%nbbnd>0) then
          WRITE(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp,bndy(i)%cnd%phi0ym,bndy(i)%cnd%phi0yp
          WRITE(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, bndy(i)%cnd%cfyp, bndy(i)%cnd%cfym, bndy(i)%cnd%cfrhs
          WRITE(10,*) bndy(i)%cnd%rphi0xm,bndy(i)%cnd%rphi0xp,bndy(i)%cnd%rphi0ym,bndy(i)%cnd%rphi0yp
          WRITE(10,*) bndy(i)%cnd%rcf0, bndy(i)%cnd%rcfxp, bndy(i)%cnd%rcfxm, &
                      bndy(i)%cnd%rcfyp, bndy(i)%cnd%rcfym, bndy(i)%cnd%rcfrhs
          WRITE(10,*) bndy(i)%cnd%dxm,bndy(i)%cnd%dxp,bndy(i)%cnd%dym,bndy(i)%cnd%dyp
          WRITE(10,*) bndy(i)%cnd%jj, bndy(i)%cnd%kk   
          WRITE(10,*) bndy(i)%cnd%docalc
        END if
        IF(bndy(i)%cnd%ncond>0) then
          WRITE(10,*) bndy(i)%cnd%jcond, bndy(i)%cnd%kcond
        END if
      end do
    end do
  CLOSE(10)

return
end subroutine save_bndstructure_rz

subroutine read_bndstructure_rz(filename)
use multigridrz
implicit none
CHARACTER(*), INTENT(IN) :: filename
INTEGER(ISZ) :: nbndy,nx,ny,nbbnd,ncond,i,ic

  OPEN(10,FILE=filename,STATUS='unknown')
    read(10,*) nbndy
    ALLOCATE(bndy(nbndy))
    do i = 1, nbndy
      read(10,*) nx,ny
      NULLIFY(bndy(i)%first)
      ALLOCATE(bndy(i)%v(nx,ny))
      read(10,*) bndy(i)%v
      read(10,*) bndy(i)%dx, bndy(i)%dy, bndy(i)%nb_conductors
      do ic = 1, bndy(i)%nb_conductors
        read(10,*) nbbnd,ncond
        call init_bnd_sublevel(bndy(i),nbbnd,ncond)
        read(10,*) bndy(i)%cnd%nbbndred
        read(10,*) bndy(i)%cnd%voltage
        WRITE(0,*) i,ic,nbndy,nbbnd,ncond,bndy(i)%cnd%nbbndred
        IF(bndy(i)%cnd%nbbnd>0) then
          read(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp,bndy(i)%cnd%phi0ym,bndy(i)%cnd%phi0yp
          read(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, bndy(i)%cnd%cfyp, bndy(i)%cnd%cfym, bndy(i)%cnd%cfrhs
          read(10,*) bndy(i)%cnd%rphi0xm,bndy(i)%cnd%rphi0xp,bndy(i)%cnd%rphi0ym,bndy(i)%cnd%rphi0yp
          read(10,*) bndy(i)%cnd%rcf0, bndy(i)%cnd%rcfxp, bndy(i)%cnd%rcfxm, &
                     bndy(i)%cnd%rcfyp, bndy(i)%cnd%rcfym, bndy(i)%cnd%rcfrhs
          read(10,*) bndy(i)%cnd%dxm,bndy(i)%cnd%dxp,bndy(i)%cnd%dym,bndy(i)%cnd%dyp
          read(10,*) bndy(i)%cnd%jj, bndy(i)%cnd%kk
          read(10,*) bndy(i)%cnd%docalc
        END if
        IF(bndy(i)%cnd%ncond>0) then
          read(10,*) bndy(i)%cnd%jcond, bndy(i)%cnd%kcond
        END if
      end do
    end do
  CLOSE(10)

return
end subroutine read_bndstructure_rz


