!     Last change:  JLV  14 Jan 2002    2:07 pm
#include "top.h"

module multigrid_common

USE constant
USE PSOR3d, ONLY:boundxy,bound0,boundnz
USE FRZmgrid
#ifdef MPIPARALLEL
  use Parallel
  use mpirz
#endif

LOGICAL :: restrictwbnd = .true.

LOGICAL :: l_mgridrz_debug=.false.

TYPE conductor_type
! structure for potential calculation close to conductors.
! The stencil for the iterative calculation of the potential f is given by
!        f(j,l) = cf0(i)   * f(j  ,l  )
!               + cfxp(i)  * f(j+1,l  )
!               + cfxm(i)  * f(j-1,l  )
!               + cfzp(i)  * f(j  ,l+1)
!               + cfzm(i)  * f(j  ,l-1)
!               + phi0xp(i)
!               + phi0xm(i)
!               + phi0zp(i)
!               + phi0zm(i)
!               + cfrhs(i) * rhs(j,l)
! for the ith grid point being close to the considered conductor and located
! at (j,l) on the grid. If there is a conductor boundary lying for example
! between j and j+1, cfxp(i) is set to zero while phi0xp(i) is set to
! the voltage of the nearest conductor.
  ! conductor voltage
  REAL(8), POINTER, DIMENSION(:) :: voltage
  ! stencil coefficients for relaxation iteration
  REAL(8), POINTER, DIMENSION(:) :: cf0, cfxp, cfxm, cfyp, cfym, cfzp, cfzm, cfrhs
  REAL(8), POINTER, DIMENSION(:) :: phi0xm, phi0xp, phi0ym, phi0yp, phi0zm, phi0zp
  REAL(8), POINTER, DIMENSION(:) :: volt0xm, volt0xp, volt0ym, volt0yp, volt0zm, volt0zp
  ! stencil coefficients for residue calculation
  REAL(8), POINTER, DIMENSION(:) :: rcf0, rcfxp, rcfxm, rcfyp, rcfym, rcfzp, rcfzm, rcfrhs
  REAL(8), POINTER, DIMENSION(:) :: rphi0xm,rphi0xp,rphi0ym,rphi0yp,rphi0zm,rphi0zp
  ! distances from grid node to conductor
  REAL(8), POINTER, DIMENSION(:) :: dxm,dxp,dym,dyp,dzm,dzp
  ! locations of nodes near conductor
  INTEGER(ISZ), POINTER, DIMENSION(:) :: jj, kk, ll
  ! locations of nodes inside conductor
  INTEGER(ISZ), POINTER, DIMENSION(:) :: jcond, kcond, lcond
  ! logical array. Calculation will be performed only when docalc=.true.
  LOGICAL(ISZ), POINTER :: docalc(:)
  INTEGER(ISZ) :: nbbnd, &      ! number of nodes near conductor
                  nbbndred, &   ! number of "red" nodes near conductor (for red-black gauss-seidel)
                  ncond         ! number of nodes inside conductor
  ! next and previous linked-list elements
  TYPE(conductor_type), POINTER :: next,prev
end TYPE conductor_type

LOGICAL(ISZ) :: bndy_allocated=.FALSE. ! flag set to true if bndy is allocated

REAL(8), parameter :: coefrelax=1._8, coefrelaxbnd=1._8 ! coefficients for relaxation steps

INTEGER(ISZ), parameter :: dirichlet=0, neumann=1, periodic=2, othertask=-1  ! boundary condition types
INTEGER(ISZ) :: ixlbnd=dirichlet, ixrbnd=dirichlet, iylbnd=dirichlet, iyrbnd=dirichlet, izlbnd=dirichlet, izrbnd=dirichlet  ! boundary conditions for each side

INTEGER(ISZ) :: v_vacuum=0, v_cond=1, v_bnd=2
INTEGER(ISZ), parameter :: egun=0, ecb=1
INTEGER(ISZ) :: bnd_method=egun

INTEGER(ISZ) :: nlevels ! number of multigrid levels
INTEGER(ISZ) :: level ! current multigrid level

#ifdef MPIPARALLEL
  INTEGER(ISZ) :: nzfine, nworkpproc, workfact=3
#endif

end module multigrid_common

module multigridrz
! module containing RZ multigrid solver

USE multigrid_common

TYPE bndptr
! type for potential calculation
  ! mask array on the grid, = .true. on nodes in vacuum (i.e. not near or in conductors)
  INTEGER(ISZ), POINTER :: v(:,:)
  LOGICAL(ISZ) :: l_powerof2 ! set to true if all grid dimensions (in number of meshes) are a power of two
  LOGICAL(ISZ) :: l_merged ! set to true if grids between tasks have been merged
  INTEGER(ISZ) :: nr, nz
  ! mesh sizes
  REAL(8) :: dr, dz
  ! number of conductors
  INTEGER(ISZ) :: nb_conductors
  ! linked-list of conductors
  TYPE(conductor_type), POINTER :: cnd,first
  INTEGER(ISZ) :: izlbnd, izrbnd  ! boundary conditions for each side
#ifdef MPIPARALLEL
  INTEGER(ISZ) :: nworkpproc
#endif
END TYPE bndptr

TYPE(bndptr), pointer :: bndy(:)

contains

subroutine init_bnd(nr,nz,dr,dz)
! intializes grid quantities according to the number of multigrid levels and grid sizes nx and nz.
USE InGen3d, ONLY:l2symtry, l4symtry
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr, dz

INTEGER(ISZ) :: i, nrp0, nzp0, nrc, nzc
REAL(8) :: drc, dzc

#ifdef MPIPARALLEL
  nzfine = nz
#endif

  ixrbnd = boundxy
  izlbnd = bound0
  izrbnd = boundnz

  nrc = nr
  nzc = nz
  drc = dr
  dzc = dz
  nlevels    = 1
#ifdef MPIPARALLEL
  nworkpproc = 1
#endif
  do WHILE(nrc>mgridrz_nmeshmin.or.nzc>mgridrz_nmeshmin)
    call evalnewgrid(nrc,nzc,drc,dzc)
    nlevels = nlevels + 1
#ifdef MPIPARALLEL
    IF(nworkpproc*nrc*nzc<=workfact*nzfine) then
      nlevels = nlevels + 1
      nworkpproc = nworkpproc*2
    END if
#endif
  end do

  nlevels=MIN(nlevels,mgridrz_nlevels_max)

  allocate(bndy(nlevels))
  bndy_allocated=.true.

  nrc = nr
  nzc = nz
  drc = dr
  dzc = dz
#ifdef MPIPARALLEL
  nworkpproc = 1
#endif
  do i = nlevels, 1, -1
    bndy(i)%l_merged=.false.
    bndy(i)%izlbnd=izlbnd
    bndy(i)%izrbnd=izrbnd
    IF(i/=nlevels) then
      call evalnewgrid(nrc,nzc,drc,dzc)
#ifdef MPIPARALLEL
      IF(nworkpproc<nslaves.and.nworkpproc*nrc*nzc<=workfact*nzfine) then
        nworkpproc = nworkpproc*2
        bndy(i)%l_merged=.true.
      END if
#endif
    END if
#ifdef MPIPARALLEL
    bndy(i)%nworkpproc = nworkpproc
    IF(my_index-nworkpproc>=0) bndy(i)%izlbnd = -(my_index-nworkpproc+1)
    IF(my_index+nworkpproc<nslaves) bndy(i)%izrbnd = -(my_index+nworkpproc+1)
    bndy(i)%nr = nrc
    bndy(i)%nz = nworkpproc*nzc
#else
    bndy(i)%nr = nrc
    bndy(i)%nz = nzc
#endif
    bndy(i)%dr = drc
    bndy(i)%dz = dzc
    NULLIFY(bndy(i)%first)
    bndy(i)%nb_conductors = 0
    ALLOCATE(bndy(i)%v(bndy(i)%nr+1,bndy(i)%nz+1))
    bndy(i)%v(:,:)=v_vacuum
  end do
#ifdef MPIPARALLEL
  IF(nslaves>1) then
    do i = 1, nlevels
      bndy(i)%l_powerof2 = .false.
    END do
    do i = nlevels,1,-1
    END do
    return
  END if
#endif
  do i = nlevels, 2, -1
    nrp0 = INT(LOG(REAL(bndy(i)%nr+1))/LOG(2.))
    nzp0 = INT(LOG(REAL(bndy(i)%nz+1))/LOG(2.))
!    nrp0 = LOG(REAL(bndy(i)%nr))/LOG(2.)+0.5
!    nzp0 = LOG(REAL(bndy(i)%nz))/LOG(2.)+0.5
    IF(2**nrp0/=bndy(i)%nr.OR.2**nzp0/=bndy(i)%nz.OR.bndy(i-1)%nr==bndy(i)%nr.OR.bndy(i-1)%nz==bndy(i)%nz) THEN
      bndy(i)%l_powerof2=.false.
    else
      bndy(i)%l_powerof2=.true.
    END if
    IF(l_mgridrz_debug) WRITE(0,*) i,bndy(i)%l_powerof2,2**nrp0,bndy(i)%nr,2**nzp0,bndy(i)%nz
  end do
  bndy(1)%l_powerof2=.false.

  IF(mgridrz_deform) then
    mgridrz_nz = nz
    call gchange("FRZmgrid",0)
    mgridrz_xfact = 1._8
    mgridrz_yfact = 1._8
  END if

  do i = nlevels, 1, -1
    WRITE(0,*) 'level,nx,nz = ',i,bndy(i)%nr,bndy(i)%nz
  END do
  return
end subroutine init_bnd

subroutine evalnewgrid(nr,nz,dr,dz)
! evaluate nr and nz at coarser level
INTEGER(ISZ), INTENT(IN OUT) :: nr, nz
REAL(8), INTENT(IN OUT) :: dr,dz

REAL(8) :: rap
INTEGER :: nrnew, nznew

  rap = dr/dz
  IF(rap>4._8/3._8.or.nr<=mgridrz_nmeshmin) then
    nznew = MAX(mgridrz_nmeshmin,nz/2)
    dz = dz * REAL(nz,8)/REAL(nznew,8)
    nz=nznew
  ELSE IF(rap<2._8/3._8.or.nz<=mgridrz_nmeshmin) then
    nrnew = MAX(mgridrz_nmeshmin,nr/2)
    dr = dr * REAL(nr,8)/REAL(nrnew,8)
    nr=nrnew
  ELSE
    nrnew = MAX(mgridrz_nmeshmin,nr/2)
    dr = dr * REAL(nr,8)/REAL(nrnew,8)
    nznew = MAX(mgridrz_nmeshmin,nz/2)
    dz = dz * REAL(nz,8)/REAL(nznew,8)
    nr=nrnew
    nz=nznew
  END if

  return
end subroutine evalnewgrid

subroutine init_bnd_sublevel(bndl,nbnd,ncond)
! initializes quantities for one grid level
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
  ALLOCATE(bndl%cnd%phi0xm(nbnd),bndl%cnd%phi0xp(nbnd),bndl%cnd%phi0zm(nbnd),bndl%cnd%phi0zp(nbnd))
  ALLOCATE(bndl%cnd%volt0xm(nbnd),bndl%cnd%volt0xp(nbnd),bndl%cnd%volt0zm(nbnd),bndl%cnd%volt0zp(nbnd))
  ALLOCATE(bndl%cnd%cf0(nbnd),bndl%cnd%cfxp(nbnd),bndl%cnd%cfxm(nbnd),bndl%cnd%cfzp(nbnd),bndl%cnd%cfzm(nbnd),bndl%cnd%cfrhs(nbnd))
  ALLOCATE(bndl%cnd%rphi0xm(nbnd),bndl%cnd%rphi0xp(nbnd),bndl%cnd%rphi0zm(nbnd),bndl%cnd%rphi0zp(nbnd))
  ALLOCATE(bndl%cnd%rcf0(nbnd),bndl%cnd%rcfxp(nbnd),bndl%cnd%rcfxm(nbnd), &
           bndl%cnd%rcfzp(nbnd),bndl%cnd%rcfzm(nbnd),bndl%cnd%rcfrhs(nbnd))
  ALLOCATE(bndl%cnd%dxm(nbnd),bndl%cnd%dxp(nbnd),bndl%cnd%dzm(nbnd),bndl%cnd%dzp(nbnd))
  ALLOCATE(bndl%cnd%jj(nbnd),bndl%cnd%kk(nbnd),bndl%cnd%docalc(nbnd))
  bndl%cnd%docalc=.false.
END if
IF(ncond>0) then
  ALLOCATE(bndl%cnd%jcond(ncond),bndl%cnd%kcond(ncond),bndl%cnd%voltage(ncond))
END if

end subroutine init_bnd_sublevel

function expandwguard(f)
! expand field from a grid to finer one. Each dimension is assumed to be a power of two.
implicit none

REAL(8), DIMENSION(0:,0:), INTENT(IN) :: f
REAL(8), DIMENSION(0:2*(SIZE(f,1)-2),0:2*(SIZE(f,2)-2)) :: expandwguard

INTEGER(ISZ) :: j, l, nxe, nze


IF(l_mgridrz_debug) WRITE(0,*) 'enter expand, level = ',level

nxe = 2*(SIZE(f,1)-3)
nze = 2*(SIZE(f,2)-3)

do l = 1, nze+1, 2
  do j = 1, nxe+1, 2
    expandwguard(j,l) = f(j/2+1,l/2+1)
  end do
end do
do l = 1, nze+1, 2
  do j = 2, nxe, 2
    expandwguard(j,l) = 0.5*(f(j/2,l/2+1)+f(j/2+1,l/2+1))
  end do
end do
do l = 2, nze, 2
  do j = 1, nxe+1, 2
    expandwguard(j,l) = 0.5*(f(j/2+1,l/2)+f(j/2+1,l/2+1))
  end do
end do
do l = 2, nze, 2
  do j = 2, nxe, 2
    expandwguard(j,l) = 0.25*(f(j/2,l/2)+f(j/2,l/2+1)+f(j/2+1,l/2)+f(j/2+1,l/2+1))
  end do
end do

expandwguard(0,:)=0._8
expandwguard(nxe+2,:)=0._8
expandwguard(:,0)=0._8
expandwguard(:,nze+2)=0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
end function expandwguard

function expandwguard_any(uold, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, xminold,  xmaxold,  zminold,  zmaxold)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguard_any(0:nxnew+2,0:nznew+2)

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jnew, knew, j, k, jp, kp
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+1), koldp(nznew+1)

IF(l_mgridrz_debug) WRITE(0,*) 'enter expand, level = ',level

nxold = SIZE(uold,1) - 1 - 2
nzold = SIZE(uold,2) - 1 - 2

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

do knew = 1, nznew+1
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = 1+MIN(nzold,INT(zz))
  koldp(knew) = kold(knew) + 1
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
do jnew = 1, nxnew+1
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
do knew = 1, nznew+1
  k = kold(knew)
  kp = koldp(knew)
  delz = ddz(knew)
  odelz = oddz(knew)
  do jnew = 1, nxnew+1
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
    expandwguard_any(jnew,knew) = uold(j, k)  * odelx * odelz &
                                + uold(jp,k)  * delx  * odelz &
                                + uold(j, kp) * odelx * delz &
                                + uold(jp,kp) * delx  * delz
  end do
END do

expandwguard_any(0,:) = 0._8
expandwguard_any(nxnew+2,:) = 0._8
expandwguard_any(:,0) = 0._8
expandwguard_any(:,nznew+2) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END function expandwguard_any

function expandwguardandbnd_any(uold, bnd, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, xminold,  xmaxold,  zminold,  zmaxold)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguardandbnd_any(0:nxnew+2,0:nznew+2)
TYPE(bndptr) :: bnd

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jnew, knew, j, k, jp, kp
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+1), koldp(nznew+1)

IF(l_mgridrz_debug) WRITE(0,*) 'enter expand, level = ',level

expandwguardandbnd_any = 0._8

nxold = SIZE(uold,1) - 1 - 2
nzold = SIZE(uold,2) - 1 - 2

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

do knew = 1, nznew+1
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = 1+MIN(nzold,INT(zz))
  koldp(knew) = kold(knew) + 1
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
do jnew = 1, nxnew+1
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
do knew = 1, nznew+1
  k = kold(knew)
  kp = koldp(knew)
  delz = ddz(knew)
  odelz = oddz(knew)
  do jnew = 1, nxnew+1
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
!    IF(.NOT.bnd%v(jnew,knew)==v_vacuum) cycle
    expandwguardandbnd_any(jnew,knew) = uold(j, k)  * odelx * odelz &
                                + uold(jp,k)  * delx  * odelz &
                                + uold(j, kp) * odelx * delz &
                                + uold(jp,kp) * delx  * delz
  end do
END do

expandwguardandbnd_any(0,:) = 0._8
expandwguardandbnd_any(nxnew+2,:) = 0._8
expandwguardandbnd_any(:,0) = 0._8
expandwguardandbnd_any(:,nznew+2) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END function expandwguardandbnd_any

function restrict_pof2(f)
! restrict field from one grid to a coarser one. Each dimension is assumed to be a power of two.
implicit none

REAL(8), DIMENSION(1:,1:), INTENT(IN) :: f
REAL(8), DIMENSION(SIZE(f,1)/2+1, SIZE(f,2)/2+1) :: restrict_pof2

INTEGER(ISZ) :: nx, nz, nxi, nzi, j, l, jj, ll
REAL(8) :: flz, frz, aflz1, aflz2, afrz1, afrz2

IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict_pof2, level = ',level

nxi = SIZE(f,1)-1
nzi = SIZE(f,2)-1
nx = nxi/2
nz = nzi/2

flz = 1._8
frz = 1._8
aflz1 = 4._8/3._8
afrz1 = 4._8/3._8
aflz2 = 16._8/9._8
afrz2 = 16._8/9._8

#ifdef MPIPARALLEL
  IF(bndy(level-1)%izlbnd<0) then
    flz=0.5_8
    aflz1 = 1._8
    aflz2 = 4._8/3._8
  endif
  IF(bndy(level-1)%izrbnd<0) then
    frz=0.5_8
    afrz1 = 1._8
    afrz2 = 4._8/3._8
  endif
#endif

do l = 2, nz
  do j = 2, nx
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = 0.25*f(jj,ll) &
                       + 0.125*(f(jj-1,ll)+f(jj,ll-1)+f(jj+1,ll)+f(jj,ll+1))  &
                       + 0.0625*(f(jj-1,ll-1)+f(jj+1,ll-1)+f(jj+1,ll+1)+f(jj-1,ll+1))
  end do
end do
l = 1
  do j = 2, nx
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = aflz1*(0.25*flz*f(jj,ll) &
                       + 0.125*(flz*f(jj-1,ll)+flz*f(jj+1,ll)+f(jj,ll+1))  &
                       + 0.0625*(f(jj+1,ll+1)+f(jj-1,ll+1)))
  end do
l = nz+1
  do j = 2, nx
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = afrz1*(0.25*frz*f(jj,ll) &
                       + 0.125*(frz*f(jj-1,ll)+f(jj,ll-1)+frz*f(jj+1,ll))  &
                       + 0.0625*(f(jj-1,ll-1)+f(jj+1,ll-1)))
  end do
j = 1
  do l = 2, nz
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = (4._8/3._8)*(0.25*f(jj,ll) &
                       + 0.125*(f(jj,ll-1)+f(jj+1,ll)+f(jj,ll+1))  &
                       + 0.0625*(f(jj+1,ll-1)+f(jj+1,ll+1)))
  end do
j = nx+1
  do l = 2, nz
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = (4._8/3._8)*(0.25*f(jj,ll) &
                       + 0.125*(f(jj,ll-1)+f(jj-1,ll)+f(jj,ll+1))  &
                       + 0.0625*(f(jj-1,ll-1)+f(jj-1,ll+1)))
  end do
j = 1
l = 1
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = aflz2*(0.25*flz*f(jj,ll) &
                       + 0.125*(flz*f(jj+1,ll)+f(jj,ll+1))  &
                       + 0.0625*(f(jj+1,ll+1)))
j = 1
l = nz+1
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = afrz2*(0.25*frz*f(jj,ll) &
                       + 0.125*(frz*f(jj+1,ll)+f(jj,ll-1))  &
                       + 0.0625*(f(jj+1,ll-1)))
j = nx+1
l = 1
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = aflz2*(0.25*flz*f(jj,ll) &
                       + 0.125*(flz*f(jj-1,ll)+f(jj,ll+1))  &
                       + 0.0625*(f(jj-1,ll+1)))
j = nx+1
l = nz+1
    jj = 2*j-1
    ll = 2*l-1
    restrict_pof2(j,l) = afrz2*(0.25*frz*f(jj,ll) &
                       + 0.125*(frz*f(jj-1,ll)+f(jj,ll-1))  &
                       + 0.0625*(f(jj-1,ll-1)))

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict_pof2, level = ',level

return
end function restrict_pof2

function restrict(uold, nxnew, nznew, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: restrict(1:nxnew+1,1:nznew+1),rap(1:nxnew+1,1:nznew+1)

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jold, kold, j, k, jp, kp
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

IF(bndy(level)%l_powerof2) then
  restrict = restrict_pof2(uold)
  return
END if

IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict, level = ',level

nxold = SIZE(uold,1) - 1
nzold = SIZE(uold,2) - 1

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1))

invdxnew = nxnew / (xmaxnew-xminnew)
invdznew = nznew / (zmaxnew-zminnew)
dxold = (xmaxold-xminold) / nxold
dzold = (zmaxold-zminold) / nzold

restrict=0._8
rap = 0._8

do kold = 1, nzold+1
  z = zminold + (kold-1)*dzold
  knew(kold) = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  ddz(kold) = (z-zminnew) * invdznew-real(knew(kold)-1)
  knewp(kold) = knew(kold)+1
  oddz(kold) = 1._8-ddz(kold)
END do

#ifdef MPIPARALLEL
    IF(bndy(level)%izlbnd<0) then
      ddz(1) = 0.5_8*ddz(1)
      oddz(1) = 0.5_8*oddz(1)
    END if
    IF(bndy(level)%izrbnd<0) then
      ddz(nzold+1) = 0.5_8*ddz(nzold+1)
      oddz(nzold+1) = 0.5_8*oddz(nzold+1)
    END if
#endif

do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = 1._8-ddx(jold)
END do

do kold = 1, nzold+1
  k = knew(kold)
  kp = knewp(kold)
  delz  =  ddz(kold)
  odelz = oddz(kold)
  do jold = 1, nxold+1
    j = jnew(jold)
    jp = jnewp(jold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    restrict(j,k)   = restrict(j,k)   + uold(jold,kold) * odelx * odelz
    restrict(jp,k)  = restrict(jp,k)  + uold(jold,kold) * delx  * odelz
    restrict(j,kp)  = restrict(j,kp)  + uold(jold,kold) * odelx * delz
    restrict(jp,kp) = restrict(jp,kp) + uold(jold,kold) * delx  * delz
    rap(j,k)   = rap(j,k)   + odelx * odelz
    rap(jp,k)  = rap(jp,k)  + delx  * odelz
    rap(j,kp)  = rap(j,kp)  + odelx * delz
    rap(jp,kp) = rap(jp,kp) + delx  * delz
  end do
end do

do k = 1, nznew+1
  do j = 1, nxnew+1
    IF(rap(j,k)/=0._8) restrict(j,k)   = restrict(j,k)   / rap(j,k)
  end do
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END function restrict

function restrict_wbnd(uold, bnd, nxnew, nznew, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: restrict_wbnd(1:nxnew+1,1:nznew+1),rap(1:nxnew+1,1:nznew+1)
TYPE(bndptr) :: bnd

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jold, kold, j, k, jp, kp, ic, ii, itot
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz, dxnew, dznew, u1, u2, u3, u4, q
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp
LOGICAL(ISZ) :: l_dxm, l_dxp, l_dzm, l_dzp
INTEGER(ISZ),allocatable :: cnt(:,:)

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold)
!  return
!END if
itot = 0
IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict, level = ',level

nxold = SIZE(uold,1) - 1
nzold = SIZE(uold,2) - 1

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1))


ALLOCATE(cnt(nxold+1,nzold+1))
invdxnew = nxnew / (xmaxnew-xminnew)
invdznew = nznew / (zmaxnew-zminnew)
dxnew = 1./invdxnew
dznew = 1./invdznew
dxold = (xmaxold-xminold) / nxold
dzold = (zmaxold-zminold) / nzold

restrict_wbnd=0._8
rap = 0._8
cnt= 0

do kold = 1, nzold+1
  z = zminold + (kold-1)*dzold
  knew(kold) = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  ddz(kold) = (z-zminnew) * invdznew-real(knew(kold)-1)
  knewp(kold) = knew(kold)+1
  oddz(kold) = 1._8-ddz(kold)
END do

#ifdef MPIPARALLEL
    IF(bndy(level)%izlbnd<0) then
      ddz(1) = 0.5_8*ddz(1)
      oddz(1) = 0.5_8*oddz(1)
    END if
    IF(bndy(level)%izrbnd<0) then
      ddz(nzold+1) = 0.5_8*ddz(nzold+1)
      oddz(nzold+1) = 0.5_8*oddz(nzold+1)
    END if
#endif

do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = 1._8-ddx(jold)
END do

do kold = 1, nzold+1
  k = knew(kold)
  kp = knewp(kold)
  delz  =  ddz(kold)
  odelz = oddz(kold)
  do jold = 1, nxold+1
#ifdef MPIPARALLEL
    IF((bndy(level)%izlbnd<0.and.kold>1).or.(bndy(level)%izrbnd<0.and.kold<nzold+1)) then
      IF(.NOT.(bnd%v(jold,kold)==v_vacuum.or.jold==1.or.jold==nxold+1)) cycle
    END if
#else
!    IF(.NOT.(bnd%v(jold,kold).or.jold==1.or.jold==nxold+1.or.kold==1.or.kold==nzold+1)) cycle
    IF(.NOT.bnd%v(jold,kold)==v_vacuum) cycle
#endif
    j = jnew(jold)
    jp = jnewp(jold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    u1 = odelx * odelz
    u2 = delx  * odelz
    u3 = odelx * delz
    u4 = delx  * delz
    q = uold(jold,kold)
    itot = itot + 1
    cnt(jold,kold) = cnt(jold,kold) + 1
    restrict_wbnd(j,k)   = restrict_wbnd(j,k)   + q * u1
    restrict_wbnd(jp,k)  = restrict_wbnd(jp,k)  + q * u2
    restrict_wbnd(j,kp)  = restrict_wbnd(j,kp)  + q * u3
    restrict_wbnd(jp,kp) = restrict_wbnd(jp,kp) + q * u4
    rap(j,k)   = rap(j,k)   + u1
    rap(jp,k)  = rap(jp,k)  + u2
    rap(j,kp)  = rap(j,kp)  + u3
    rap(jp,kp) = rap(jp,kp) + u4 
  end do
end do

!GO TO 10
do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do ii = 1, bnd%cnd%nbbnd
    jold = bnd%cnd%jj(ii)
    kold = bnd%cnd%kk(ii)
    IF(.NOT.(bnd%v(jold,kold)==v_bnd.and.bnd%cnd%docalc(ii))) cycle
    j = jnew(jold)
    jp = jnewp(jold)
    k = knew(kold)
    kp = knewp(kold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    delz  =  ddz(kold)
    odelz = oddz(kold)
    u1 = odelx * odelz
    u2 = delx  * odelz
    u3 = odelx * delz
    u4 = delx  * delz
    l_dxm = .true.
    l_dxp = .true.
    l_dzm = .true.
    l_dzp = .true.
!GO TO 10
    IF(bnd%cnd%dxm(ii)<delx*dxnew)  l_dxm = .false.
    IF(bnd%cnd%dxp(ii)<odelx*dxnew) l_dxp = .false.
    IF(bnd%cnd%dzm(ii)<delz*dznew)  l_dzm = .false.
    IF(bnd%cnd%dzp(ii)<odelz*dznew) l_dzp = .false.
    IF((l_dxm.and.l_dxp).OR.(l_dzm.and.l_dzp)) then
      cycle
    else
      IF(l_dxm) then
!        u2=u2+u1
        u1=0
!        u4=u4+u3
        u3=0
      END if
      IF(l_dxp) then
!        u1=u1+u2
        u2=0
!        u3=u3+u4
        u4=0
      END if
      IF(l_dzm) then
!        u3=u3+u1
        u1=0
!        u4=u4+u2
        u2=0
      END if
      IF(l_dzp) then
!        u1=u1+u3
        u3=0
!        u2=u2+u4
        u4=0
      END if
    END if
!10 continue
    q = uold(jold,kold)
    itot = itot + 1
    cnt(jold,kold) = cnt(jold,kold) + 1
    restrict_wbnd(j,k)   = restrict_wbnd(j,k)   + q * u1
    restrict_wbnd(jp,k)  = restrict_wbnd(jp,k)  + q * u2
    restrict_wbnd(j,kp)  = restrict_wbnd(j,kp)  + q * u3
    restrict_wbnd(jp,kp) = restrict_wbnd(jp,kp) + q * u4
    rap(j,k)   = rap(j,k)   + u1
    rap(jp,k)  = rap(jp,k)  + u2
    rap(j,kp)  = rap(j,kp)  + u3
    rap(jp,kp) = rap(jp,kp) + u4 

  ENDDO
END do
!10 continue
do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do ii = 1, bnd%cnd%ncond
    jold = bnd%cnd%jcond(ii)
    kold = bnd%cnd%kcond(ii)
    j = jnew(jold)
    jp = jnewp(jold)
    k = knew(kold)
    kp = knewp(kold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    delz  =  ddz(kold)
    odelz = oddz(kold)
    u1 = odelx * odelz
    u2 = delx  * odelz
    u3 = odelx * delz
    u4 = delx  * delz
    q = uold(jold,kold)
    itot = itot + 1
    cnt(jold,kold) = cnt(jold,kold) + 1
    restrict_wbnd(j,k)   = restrict_wbnd(j,k)   + q * u1
    restrict_wbnd(jp,k)  = restrict_wbnd(jp,k)  + q * u2
    restrict_wbnd(j,kp)  = restrict_wbnd(j,kp)  + q * u3
    restrict_wbnd(jp,kp) = restrict_wbnd(jp,kp) + q * u4
    rap(j,k)   = rap(j,k)   + u1
    rap(jp,k)  = rap(jp,k)  + u2
    rap(j,kp)  = rap(j,kp)  + u3
    rap(jp,kp) = rap(jp,kp) + u4 
  end do
END do
!10 continue

do k = 1, nznew+1
  do j = 1, nxnew+1
    IF(rap(j,k)/=0._8) restrict_wbnd(j,k)   = restrict_wbnd(j,k)   / rap(j,k)
  end do
end do
do k = 1, nzold+1
  do j = 1, nxold+1
    IF(cnt(j,k)>1) then
      do kold = 1, nzold+1
        do jold = 1, nxold+1
          IF(.NOT.bnd%v(jold,kold)==v_vacuum) cycle
          IF(j==jold.and.k==kold) WRITE(0,*) level,': #1# ',j,k,cnt(j,k)
        END do
      END do
      do ic = 1, bnd%nb_conductors
        IF(ic==1) then
          bnd%cnd => bnd%first
        else
          bnd%cnd => bnd%cnd%next
        END if
        do ii = 1, bnd%cnd%nbbnd
          jold = bnd%cnd%jj(ii)
          kold = bnd%cnd%kk(ii)
          IF(.NOT.(bnd%v(jold,kold)==v_bnd.and.bnd%cnd%docalc(ii))) cycle
          IF(j==jold.and.k==kold) WRITE(0,*) level,': #2# ',j,k,cnt(j,k)
        END do
        do ii = 1, bnd%cnd%ncond
          jold = bnd%cnd%jcond(ii)
          kold = bnd%cnd%kcond(ii)
          IF(j==jold.and.k==kold) WRITE(0,*) level,': #3# ',j,k,cnt(j,k)
        END do
      END do
    END if
  END do
END do
DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp, cnt)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END function restrict_wbnd

subroutine interp_bndwguard(unew, uold, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! interpolate boundary values from two conformal grids. Grid is assumed to have guard cells.
implicit none
REAL(8), DIMENSION(0:,0:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew

INTEGER(ISZ) :: nxnew, nznew, nxold, nzold
INTEGER(ISZ) :: jold, kold, jnew, knew
REAL(8) :: x, z, dxold, dzold, dxnew, dznew, ddx, ddz

nxnew = SIZE(unew,1)-2 - 1
nznew = SIZE(unew,2)-2 - 1
nxold = SIZE(uold,1)-2 - 1
nzold = SIZE(uold,2)-2 - 1

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
dxold = (xmaxold-xminold) / nxold
dzold = (zmaxold-zminold) / nzold

knew = 1
kold = 1
  do jnew = 1, nxnew+1
    x = xminnew+(jnew-1)*dxnew
    jold = MIN(1 + INT((x-xminold) / dxold), nxold)
    ddx = (x-xminold)/dxold-REAL(jold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddx) &
                    + uold(jold+1,kold)   * ddx
  end do
knew = nznew+1
kold = nzold+1
  do jnew = 1, nxnew+1
    x = xminnew+(jnew-1)*dxnew
    jold = MIN(1 + INT((x-xminold) / dxold), nxold)
    ddx = (x-xminold)/dxold-REAL(jold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddx) &
                    + uold(jold+1,kold)   * ddx
  end do
jnew = 1
jold = 1
  do knew = 2, nznew
    z = zminnew+(knew-1)*dznew
    kold = MIN(1 + INT((z-zminold) / dzold), nzold)
    ddz = (z-zminold)/dzold-REAL(kold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddz) &
                    + uold(jold,kold+1)   * ddz
  end do
jnew = nxnew+1
jold = nxold+1
  do knew = 2, nznew
    z = zminnew+(knew-1)*dznew
    kold = MIN(1 + INT((z-zminold) / dzold), nzold)
    ddz = (z-zminold)/dzold-REAL(kold-1)
    unew(jnew,knew) = uold(jold,kold)     * (1.-ddz) &
                    + uold(jold,kold+1)   * ddz
  end do

return
END subroutine interp_bndwguard

subroutine relaxbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,nc,voltfact)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc
REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(1:,1:)!rhs(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic, nrf, nzi, nzf
REAL(8) :: dt, dt0
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs

IF(l_mgridrz_debug) WRITE(0,*) 'enter relax, level = ',level

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

IF(ixrbnd==dirichlet) then
  nrf=nr-1
else
  nrf=nr
END if
IF(bndy(level)%izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(bndy(level)%izrbnd==dirichlet) then
  nzf=nz-1
else
  nzf=nz
END if

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
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + bnd%cnd%cfrhs(ii)*rhs(j,l)
      else
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + bnd%cnd%cfrhs(ii)*rhs(j,l)
      END if
    ENDDO
  END do
  do l = nzi, nzf+1
    IF(jsw==2) then! origin
      j = 1
      IF(bnd%v(j,l)==v_vacuum) &
      f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                 + 4._8*dt0*f(j+1,l)/dr**2   &
                                 + (dt0/dz**2)*(f(j,l+1)+f(j,l-1)) &
                                 - dt0*rhs(j,l)
    END if
    do j = jsw+1, nrf+1, 2
      IF(bnd%v(j,l)==v_vacuum) &
        f(j,l) = cf0 * f(j,l) &
                                 + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 + cfrhs*rhs(j,l)
    end do
    jsw = 3-jsw
  end do
  lsw = 3-lsw
END do !redblack=1, 2

call updateguardcellsrz(f=f,level=level)

#ifdef MPIPARALLEL
  call exchange_fbndz(f,level)
#endif

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbndrzwguard

#ifdef MPIPARALLEL
  subroutine exchange_fbndz(f,level)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: level

    INTEGER(ISZ) :: nr,nz
    INTEGER(ISZ) :: p_up, p_down

    p_up   = -bndy(level)%izrbnd-1
    p_down = -bndy(level)%izlbnd-1

    nr = SIZE(f,1)-3
    nz = SIZE(f,2)-3

    ! send
    IF(bndy(level)%izlbnd<0) call mpi_send_real_array(f(:,2), p_down, 0)
    IF(bndy(level)%izrbnd<0) call mpi_send_real_array(f(:,nz), p_up, 0)

    ! receive
    IF(bndy(level)%izrbnd<0) f(:,nz+2) = mpi_recv_real_array(SIZE(f(:,nz)),p_up,0)
    IF(bndy(level)%izlbnd<0) f(:,0)    = mpi_recv_real_array(SIZE(f(:,0 )),p_down,0)

  end subroutine exchange_fbndz
  subroutine exchange_rhobndz(rho,level)
    REAL(8), INTENT(IN OUT) :: rho(1:,1:)!rho(1:nr+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: level

    INTEGER(ISZ) :: nr,nz
    INTEGER(ISZ) :: p_up, p_down

    p_up   = -bndy(level)%izrbnd-1
    p_down = -bndy(level)%izlbnd-1

    nr = SIZE(rho,1)-1
    nz = SIZE(rho,2)-1

    ! send
    IF(bndy(level)%izlbnd<0) call mpi_send_real_array(rho(:,1), p_down, 1)
    IF(bndy(level)%izrbnd<0) call mpi_send_real_array(rho(:,nz+1), p_up, 1)

    ! receive
    IF(bndy(level)%izrbnd<0) rho(:,nz+1) = 0.5_8*rho(:,nz+1) + 0.5_8*mpi_recv_real_array(SIZE(rho(:,nz+1)),p_up,1)
    IF(bndy(level)%izlbnd<0) rho(:,1)    = 0.5_8*rho(:,1)    + 0.5_8*mpi_recv_real_array(SIZE(rho(:,1 ))  ,p_down,1)

  end subroutine exchange_rhobndz
  subroutine merge_work(f,level)
    REAL(8), INTENT(IN OUT) :: f(1:,1:)!f(1:nr+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: level

    INTEGER(ISZ) :: nz, p_up, p_down

    nz     = bndy(level-1)%nz
    p_up   = -bndy(level)%izrbnd-1
    p_down = -bndy(level)%izlbnd-1

    IF(MOD(my_index/bndy(level)%nworkpproc,2)==0) then
    ! send up
      call mpi_send_real_array(PACK(f(:,1:nz/2),     .TRUE.), p_up,   2)
    ! receive up
      f(:,nz/2+2:nz+1) = RESHAPE(mpi_recv_real_array(SIZE(f(:,nz/2+2:nz+1)), p_up,2), &
                                                    SHAPE(f(:,nz/2+2:nz+1)))
    else
    ! send down
      call mpi_send_real_array(PACK(f(:,nz/2+2:nz+1),.TRUE.), p_down, 2)
    ! receive down
      f(:,1:nz/2)      = RESHAPE(mpi_recv_real_array(SIZE(f(:,1:nz/2)),       p_down,2), &
                                                    SHAPE(f(:,1:nz/2)))
    END if

  end subroutine merge_work
#endif

function residbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,voltfact,l_zerolastz)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: f(0:,0:)!f(0:nx+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(:,:)!rhs(nx+1,nz+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dr, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2) :: residbndrzwguard
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, l, ii, ic, nrf, nzi, nzf
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz

IF(l_mgridrz_debug) WRITE(0,*) 'enter resid, level = ',level

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
    residbndrzwguard(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = 0._8
  end do
END do

IF(ixrbnd==dirichlet) then
  nrf=nr-1
else
  nrf=nr
END if
IF(bndy(level)%izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(bndy(level)%izrbnd==dirichlet) then
  nzf=nz-1
else
  nzf=nz
END if

do l = nzi, nzf+1
  j = 1
  IF(bnd%v(j,l)==v_vacuum) &
  residbndrzwguard(j,l) = (cf0+2._8/dr**2) * f(j,l) - 4._8*f(j+1,l)/dr**2   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 + rhs(j,l)
  do j = 2, nrf+1
     IF(bnd%v(j,l)==v_vacuum) &
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
      IF(bnd%v(j,l)==v_bnd.and.bnd%cnd%docalc(ii)) &
      residbndrzwguard(j,l) = bnd%cnd%rcf0(ii)*f(j,l) &
                        + bnd%cnd%rcfxp(ii)*f(j+1,l) &
                        + bnd%cnd%rcfzp(ii)*f(j,l+1)+bnd%cnd%rcfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%rphi0xp(ii) &
                        + bnd%cnd%rphi0zp(ii)+bnd%cnd%rphi0zm(ii)) &
                        + bnd%cnd%rcfrhs(ii)*rhs(j,l)
    else
      IF(bnd%v(j,l)==v_bnd.and.bnd%cnd%docalc(ii)) &
      residbndrzwguard(j,l) = bnd%cnd%rcf0(ii)*f(j,l) &
                        + bnd%cnd%rcfxp(ii)*f(j+1,l)+bnd%cnd%rcfxm(ii)*f(j-1,l) &
                        + bnd%cnd%rcfzp(ii)*f(j,l+1)+bnd%cnd%rcfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%rphi0xp(ii)+bnd%cnd%rphi0xm(ii) &
                        + bnd%cnd%rphi0zp(ii)+bnd%cnd%rphi0zm(ii)) &
                        + bnd%cnd%rcfrhs(ii)*rhs(j,l)
    END if
  ENDDO
END do

IF(l_zerolastz) residbndrzwguard(:,nz+1) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit resid, level = ',level

return
end function residbndrzwguard

RECURSIVE subroutine mgbndrzwguard(j, u, rhs, bnd, nr, nz, dr, dz, npre, npost, ncycle, sub, relax_only, npmin)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: j, nr, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(:,:), INTENT(IN OUT) :: u
REAL(8), DIMENSION(:,:), INTENT(IN) :: rhs
REAL(8) :: dr, dz
TYPE(bndptr) :: bnd(:)
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(:,:), allocatable :: res, v
INTEGER(ISZ) :: i
INTEGER :: nrnext, nznext, nzresmin, nzresmax, nzres
REAL(8) :: drnext, dznext, voltf

REAL(8) :: error1, error2

level = j

IF(l_mgridrz_debug) WRITE(0,*) 'enter mg, level = ',level

IF(sub) then
  voltf = 0._8
else
  voltf = 1._8
END if

IF(j<=npmin .or. relax_only) then
!error1 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,level=j)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf)
!error2 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'relax only, level = ',level,' error1/2 = ',error1,error2
else
  nrnext = bnd(j-1)%nr
  nznext = bnd(j-1)%nz
  drnext = bnd(j-1)%dr
  dznext = bnd(j-1)%dz
  ALLOCATE(res(nrnext+1,nznext+1),v(0:nrnext+2,0:nznext+2))
!error1 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,level=j)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf)
!error2 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'relax 1, level = ',level,' error1/2 = ',error1,error2
!error1 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  IF(bnd(level-1)%l_merged) then
#ifdef MPIPARALLEL
    IF(MOD(my_index/bnd(level)%nworkpproc,2)==0) then
      nzresmin = 1
      nzresmax = nznext/2+1
    else
      nzresmin = nznext/2+1
      nzresmax = nznext+1
    END if
    nzres = nznext/2
#endif
  else
    nzresmin = 1
    nzresmax = nznext+1
    nzres = nznext
  END if
  IF(restrictwbnd) then
  res(:,nzresmin:nzresmax) = restrict_wbnd( &
                             residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.), &
                             bnd(j),nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  else
  res(:,nzresmin:nzresmax) = restrict( &
                             residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.), &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
#ifdef MPIPARALLEL
  IF(bnd(level-1)%l_merged) then
    call merge_work(res,level)
  else
    call exchange_rhobndz(res,level-1)
  END if
#endif
  call apply_voltage(res,bnd(j-1),0._8)
  v = 0.0_8
  do i = 1, ncycle  !(1=V cycles, 2=W cycle)
    call mgbndrzwguard(j=j-1, u=v, rhs=res, bnd=bnd(1:j-1), nr=nrnext, nz=nznext, dr=drnext, dz=dznext, npre=npre, npost=npost, &
                   ncycle=ncycle, sub=.TRUE., relax_only=.FALSE., npmin=npmin)
    level = j
  end do
  call apply_voltagewguard(v,bnd(j-1),0._8)
call updateguardcellsrz(f=v,level=j-1)
!  IF(bnd(j)%l_powerof2) then
!    u = u + expandwguard(v(:,nzresmin-1:nzresmax+1))
!  else
  IF(restrictwbnd) then
    u = u + expandwguardandbnd_any(v(:,nzresmin-1:nzresmax+1),bnd(j),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  else
    u = u + expandwguard_any(v(:,nzresmin-1:nzresmax+1),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
!  END if
#ifdef MPIPARALLEL
  call exchange_fbndz(u,level)
#endif
!error2 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'V cycle, level = ',level,' error1/2 = ',error1,error2
!error1 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,level=j)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npost,voltfact=voltf)
!error2 = SUM(ABS(residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'relax 2, level = ',level,' error1/2 = ',error1,error2
  DEALLOCATE(res,v)
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit mg, level = ',level

return
end subroutine mgbndrzwguard

subroutine updateguardcellsrz(f,level)
! update guard cells values according to boundary conditions.
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)
INTEGER(ISZ) :: level

INTEGER(ISZ) :: ixmax, izmax

ixmax=SIZE(f,1)
izmax=SIZE(f,2)

select case (ixrbnd)
    case (dirichlet)
      f(ixmax,:) = 2._8*f(ixmax-1,:) - f(ixmax-2,:)
    case (neumann)
      f(ixmax,:) = f(ixmax-2,:)
    case default
end select
select case (bndy(level)%izlbnd)
    case (dirichlet)
      f(:,1) = 2._8*f(:,2) - f(:,3)
    case (neumann)
      f(:,1) = f(:,3)
    case (periodic)
      f(:,1) = f(:,izmax-1)
    case default
end select
select case (bndy(level)%izrbnd)
    case (dirichlet)
      f(:,izmax) = 2._8*f(:,izmax-1) - f(:,izmax-2)
    case (neumann)
      f(:,izmax) = f(:,izmax-2)
    case (periodic)
      f(:,izmax) = f(:,2)
    case default
end select

end subroutine updateguardcellsrz

subroutine apply_voltage(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors
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
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = coef_voltage*bnd%cnd%voltage(i)
  end do
END do

return
end subroutine apply_voltage

subroutine apply_voltagewguard(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors. Grid is assumed to have guard cells.
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
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = coef_voltage*bnd%cnd%voltage(i)
  end do
END do

return
end subroutine apply_voltagewguard

subroutine solve_full_multigridrz(u,rhoinit,bnd,nr0,nz0,length_r,length_z,accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nr0, nz0, &       ! number of meshes in R and Z (excluing guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle            ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
INTEGER(ISZ), INTENT(IN OUT) :: nrecurs_min   ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:)             ! u(0:nr0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:)           ! rhoinit(nr0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_r, length_z, &  ! grid lengths in R and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

INTEGER(ISZ) :: i, j, nrc, nzc, npmin
REAL(8), allocatable, DIMENSION(:,:) :: f, fold
REAL(8) :: dr, dz
TYPE rhoptr
  REAL(8), POINTER:: a(:,:)
END TYPE rhoptr
TYPE(rhoptr), ALLOCATABLE :: rho(:)
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged
INTEGER(ISZ):: nzresmin,nzresmax,nzres

do_calc=.true.
has_diverged = .false.

nrc = bnd(nlevels)%nr
nzc = bnd(nlevels)%nz

ALLOCATE(f(0:bnd(1)%nr+2,0:bnd(1)%nz+2), rho(nlevels))
ALLOCATE(rho(nlevels)%a(nrc+1,nzc+1))
rho(nlevels)%a = rhoinit
f = 0._8
!f(1:bnd(1)%nr+1,1:bnd(1)%nz+1) = restrict(uold=u, nxnew=bnd(1)%nr, nznew=bnd(1)%nz, &
!                     xminold=0._8, xmaxold=length_r, zminold=0._8, zmaxold=length_z, &
!                     xminnew=0._8, xmaxnew=length_r, zminnew=0._8, zmaxnew=length_z)

do i = nlevels-1, 1, -1
  ALLOCATE(rho(i)%a(bnd(i)%nr+1,bnd(i)%nz+1))
  rho(i)%a = 0.
  level = i
  IF(bnd(level)%l_merged) then
#ifdef MPIPARALLEL
    IF(MOD(my_index/bnd(level+1)%nworkpproc,2)==0) then
      nzresmin = 1
      nzresmax = bnd(level)%nz/2+1
    else
      nzresmin = bnd(level)%nz/2+1
      nzresmax = bnd(level)%nz+1
    END if
    nzres = bnd(level)%nz/2
#endif
  else
    nzresmin = 1
    nzresmax = bnd(level)%nz+1
    nzres = bnd(level)%nz
  END if
  level = i+1
  IF(restrictwbnd) then
  rho(i)%a(:,nzresmin:nzresmax) = restrict_wbnd(rho(i+1)%a,bnd(level),bnd(i)%nr,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  else
  rho(i)%a(:,nzresmin:nzresmax) = restrict(rho(i+1)%a,bnd(i)%nr,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
#ifdef MPIPARALLEL
  IF(bnd(level-1)%l_merged) then
    call merge_work(rho(i)%a,level)
  else
    call exchange_rhobndz(rho(i)%a,level-1)
  END if
#endif
end do

npmin=MIN(nrecurs_min,MAX(nlevels,1))
main_loop: do i = 1, nlevels
  IF(l_mgridrz_debug) WRITE(0,*) 'begin main loop, level = ',i
  level = i
  nrc = bnd(i)%nr
  nzc = bnd(i)%nz
  dr = bnd(i)%dr
  dz = bnd(i)%dz
  IF(i>1) then
    IF(ALLOCATED(fold)) DEALLOCATE(fold)
    ALLOCATE(fold(0:bnd(i-1)%nr+2,0:bnd(i-1)%nz+2))
    fold = f
    DEALLOCATE(f)
    ALLOCATE(f(0:nrc+2,0:nzc+2))
    IF(bnd(i-1)%l_merged) then
#ifdef MPIPARALLEL
      IF(MOD(my_index/bnd(i)%nworkpproc,2)==0) then
        nzresmin = 1
        nzresmax = bnd(i-1)%nz/2+1
      else
        nzresmin = bnd(i-1)%nz/2+1
        nzresmax = bnd(i-1)%nz+1
      END if
      nzres = bnd(i-1)%nz/2
#endif
    else
      nzresmin = 1
      nzresmax = bnd(i-1)%nz+1
      nzres = bnd(i-1)%nz
    END if
!    IF(bnd(i)%l_powerof2) then
!      f = expandwguard(fold(:,nzresmin-1:nzresmax+1))
!   else
!      f = expandwguard_any(fold(:,nzresmin-1:nzresmax+1),nrc,nzc, &
!                      xminold=0._8, xmaxold=length_r, zminold=0._8, zmaxold=length_z, &
!                      xminnew=0._8, xmaxnew=length_r, zminnew=0._8, zmaxnew=length_z)
!    END if
      f = expandwguardandbnd_any(fold(:,nzresmin-1:nzresmax+1),bnd(i),nrc,nzc, &
                      xminold=0._8, xmaxold=length_r, zminold=0._8, zmaxold=length_z, &
                      xminnew=0._8, xmaxnew=length_r, zminnew=0._8, zmaxnew=length_z)
#ifdef MPIPARALLEL
    call exchange_fbndz(f,i)
#endif
    call apply_voltagewguard(f,bnd(i),1._8)
    DEALLOCATE(fold)
  END if
  ALLOCATE(fold(0:nrc+2,0:nzc+2))
  fold = f
!f = 0. ! for debugging purpose only
!  npmin=MIN(nrecurs_min,MAX(nlevels-1,1))
  do_calc=.true.
  do while(do_calc)
    average_residue_init=SUM(ABS(residbndrzwguard(f=RESHAPE((/(0._8*j,j=1,(nrc+3)*(nzc+3))/), (/nrc+3,nzc+3/)) &
                             ,rhs=rho(i)%a,bnd=bnd(i),nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nrc*nzc)
#ifdef MPIPARALLEL
    average_residue_init = mpi_global_compute_real(average_residue_init,MPI_SUM)/nslaves
#endif
    if(average_residue_init==0._8) then
      average_residue=0.
      average_residue_init=1.
      cycle main_loop
    END if
    average_residue_prev=average_residue_init
    do  j = 1, nc
      call mgbndrzwguard(j=i,u=f,rhs=rho(i)%a,bnd=bnd,nr=nrc,nz=nzc,dr=dr,dz=dz,npre=npre,npost=npost,ncycle=ncycle,sub=.FALSE., &
      relax_only=.false.,npmin=npmin)
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
      average_residue=SUM(ABS(residbndrzwguard(f=f,rhs=rho(i)%a,bnd=bnd(i), &
                              nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nrc*nzc)
#ifdef MPIPARALLEL
      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(average_residue/average_residue_prev>=1. .and. average_residue/average_residue_init>=1.) then
        IF(npmin<i) then
          has_diverged = .true.
          WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          f=fold
          do_calc=.true.
          npmin=npmin+1
          WRITE(0,*) '        trying npmin = ',npmin
          exit
        else
          WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          f=fold
          do_calc=.false.
          exit
        END if
        exit
      END if
      IF(i==nlevels) then
        IF(average_residue/average_residue_init <= accuracy) then
          do_calc=.false.
          exit
        END if
      else
        IF(average_residue/average_residue_init <= sub_accuracy) then
          do_calc=.false.
          exit
        END if
      END if
      average_residue_prev=average_residue
    end do
    IF(j>=nc) do_calc=.false.
  END do
end do main_loop

WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') average_residue/average_residue_init,j
u = f
call updateguardcellsrz(f=u,level=nlevels)

do i = nlevels, 1, -1
  DEALLOCATE(rho(i)%a)
end do
DEALLOCATE(f, rho)

IF(nrecurs_min/=npmin) then
  WRITE(0,'("WARNING multigridrz, nrecurs_min = ",i2,", npmin = ",i2,". Setting nrecurs_min = ",i2,".")') nrecurs_min, npmin, npmin
  nrecurs_min = npmin
END if

return
end subroutine solve_full_multigridrz

subroutine solve_multigridrz(u,rhoinit,bnd,nr0,nz0,length_r,length_z,accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nr0, nz0, &       ! number of meshes in R and Z (excluing guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle            ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
INTEGER(ISZ), INTENT(IN OUT) :: nrecurs_min   ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:)             ! u(0:nr0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:)           ! rhoinit(nr0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_r, length_z, &  ! grid lengths in R and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

INTEGER(ISZ) :: i, j, nrc, nzc, npmin
!REAL(8), allocatable, DIMENSION(:,:) :: f!, fold
REAL(8) :: uold(SIZE(u,1),SIZE(u,2))
REAL(8) :: dr, dz
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged
INTEGER(ISZ):: nzresmin,nzresmax,nzres

do_calc=.true.
has_diverged = .false.

nrc = bnd(nlevels)%nr
nzc = bnd(nlevels)%nz

!ALLOCATE(f(0:bnd(nlevels)%nr+2,0:bnd(nlevels)%nz+2))
!u = 0._8

npmin=MIN(nrecurs_min,MAX(nlevels,1))
  level = nlevels
  i = nlevels
  nrc = bnd(i)%nr
  nzc = bnd(i)%nz
  dr = bnd(i)%dr
  dz = bnd(i)%dz
!  ALLOCATE(fold(0:bnd(nlevels)%nr+2,0:bnd(nlevels)%nz+2))
  do_calc=.true.
  do while(do_calc)
    average_residue_init=SUM(ABS(residbndrzwguard(f=RESHAPE((/(0._8*j,j=1,(nrc+3)*(nzc+3))/), (/nrc+3,nzc+3/)) &
                             ,rhs=rhoinit,bnd=bnd(i),nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nrc*nzc)
#ifdef MPIPARALLEL
    average_residue_init = mpi_global_compute_real(average_residue_init,MPI_SUM)/nslaves
#endif
    if(average_residue_init==0._8) then
      average_residue=0.
      average_residue_init=1.
!      cycle main_loop
    END if
    average_residue_prev=average_residue_init
    do  j = 1, nc
      uold = u
      call mgbndrzwguard(j=i,u=u,rhs=rhoinit,bnd=bnd,nr=nrc,nz=nzc,dr=dr,dz=dz,npre=npre,npost=npost,ncycle=ncycle,sub=.FALSE., &
      relax_only=.false.,npmin=npmin)
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
      average_residue=SUM(ABS(residbndrzwguard(f=u,rhs=rhoinit,bnd=bnd(i), &
                              nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nrc*nzc)
#ifdef MPIPARALLEL
      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(average_residue/average_residue_prev>=1. .and. average_residue/average_residue_init>=1.) then
        IF(npmin<i) then
          has_diverged = .true.
          WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.true.
          npmin=npmin+1
          WRITE(0,*) '        trying npmin = ',npmin
          exit
        else
          WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.false.
          exit
        END if
        exit
      END if
      IF(i==nlevels) then
        IF(average_residue/average_residue_init <= accuracy) then
          do_calc=.false.
          exit
        END if
      else
        IF(average_residue/average_residue_init <= sub_accuracy) then
          do_calc=.false.
          exit
        END if
      END if
      average_residue_prev=average_residue
    end do
    IF(j>=nc) do_calc=.false.
  end do

WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') average_residue/average_residue_init,j
!u = f
call updateguardcellsrz(f=u,level=nlevels)

IF(nrecurs_min/=npmin) then
  WRITE(0,'("WARNING multigridrz, nrecurs_min = ",i2,", npmin = ",i2,". Setting nrecurs_min = ",i2,".")') nrecurs_min, npmin, npmin
  nrecurs_min = npmin
END if

return
end subroutine solve_multigridrz

subroutine solve_multigridrz_jlv2(u,rhoinit,bnd,nr0,nz0,length_r,length_z, &
                             accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nr0, nz0, &       ! number of meshes in R and Z (excluing guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle            ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
INTEGER(ISZ), INTENT(IN OUT) :: nrecurs_min   ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:)             ! u(0:nr0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:)           ! rhoinit(nr0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_r, length_z, &  ! grid lengths in X, Y and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

INTEGER(ISZ) :: i, j, nrc, nzc, npmin
!REAL(8), allocatable, DIMENSION(:,:) :: f!, fold
REAL(8) :: uold(SIZE(u,1),SIZE(u,2)), maxerr, maxerr_old
REAL(8) :: dr, dz
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged

do_calc=.true.
has_diverged = .false.

nrc = bnd(nlevels)%nr
nzc = bnd(nlevels)%nz

!ALLOCATE(f(0:bnd(nlevels)%nx+2,0:bnd(nlevels)%nz+2))
!u = 0._8

npmin=MIN(nrecurs_min,MAX(nlevels,1))
  level = nlevels
  i = nlevels
  nrc = bnd(i)%nr
  nzc = bnd(i)%nz
  dr = bnd(i)%dr
  dz = bnd(i)%dz
!  ALLOCATE(fold(0:bnd(nlevels)%nx+2,0:bnd(nlevels)%nz+2))
  do_calc=.true.
  do while(do_calc)
!    average_residue_init=SUM(ABS(residbnd3dwguard(f=RESHAPE((/(0._8*j,j=1,(nxc+3)*(nyc+3)*(nzc+3))/), (/nxc+3,nyc+3,nzc+3/)) &
!                   ,rhs=rhoinit,bnd=bnd(i),nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef PARALLEL
!    average_residue_init = mpi_global_compute_real(average_residue_init,MPI_SUM)/nslaves
#endif
!    if(average_residue_init==0._8) then
!      average_residue=0.
!      average_residue_init=1.
!      cycle main_loop
!    END if
!    average_residue_prev=average_residue_init
    maxerr = 1.
    do  j = 1, nc
  !    fold = f
      uold=u
      call mgbndrzwguard(j=i,u=u,rhs=rhoinit,bnd=bnd,nr=nrc,nz=nzc,dr=dr,dz=dz, &
                         npre=npre,npost=npost,ncycle=ncycle,sub=.FALSE., &
      relax_only=.false.,npmin=npmin)
      maxerr_old = maxerr
      maxerr = maxval(abs(u-uold))
!      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
!      average_residue=SUM(ABS(residbnd3dwguard(f=u,rhs=rhoinit,bnd=bnd(i), &
!                              nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef PARALLEL
!      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
!      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(maxerr/maxerr_old>=1..and.j>1) then
        IF(npmin<i) then
          has_diverged = .true.
          WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',maxerr_old
          WRITE(0,*) '        average current residue = ',maxerr
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.true.
          npmin=npmin+1
          WRITE(0,*) '        trying npmin = ',npmin
          exit
        else
          WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',maxerr_old
          WRITE(0,*) '        average current residue = ',maxerr    
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.false.
          exit
        END if
        exit
      END if
      IF(maxerr <= accuracy) then
        do_calc=.false.
        exit
      END if
    end do
    IF(j>=nc) do_calc=.false.
  end do

WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') maxerr,j
!u = f
call updateguardcellsrz(f=u,level=nlevels)

IF(nrecurs_min/=npmin) then
  WRITE(0,'("WARNING multigridrz, nrecurs_min = ",i2,", npmin = ",i2,". Setting nrecurs_min = ",i2,".")') nrecurs_min, npmin, npmin
  nrecurs_min = npmin
END if

return
end subroutine solve_multigridrz_jlv2

END module multigridrz

subroutine multigridrzf(iwhich,u0,rho0,nr0,nz0,dr0,dz0,accuracy,ncmax,npre,npost,ncycle)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich,nr0, nz0, ncmax,npre,npost,ncycle
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN) :: rho0(nr0+1,nz0+1)
REAL(8), INTENT(IN) :: dr0, dz0, accuracy

real(8) :: u(0:nr0+2,0:nz0+2)

WRITE(0,*) "JFDSKFJDSDLFJ"
IF(ncmax==0) return

  IF(iwhich==0.or.iwhich==1) then
    if(.not.bndy_allocated) call init_bnd(nr0,nz0,dr0,dz0)
  END if

  IF(iwhich==1) return

  u(1:nr0+1,:)=u0(1:nr0+1,:)
  u(0,:) = u(2,:)! 0._8
  u(nr0+2,:)=u(nr0+1,:)!0._8

! call solve_full_multigridrz(u=u,rhoinit=-rho0/eps0,bnd=bndy,nr0=nr0,nz0=nz0, &
!                         length_r=nr0*dr0,length_z=nz0*dz0, &
!                         accuracy=accuracy,sub_accuracy=mgridrz_sub_accuracy, &
!                        nc=ncmax,npre=npre,npost=npost,ncycle=ncycle,nrecurs_min=mgridrz_nrecurs_min)

! call solve_multigridrz(u=u,rhoinit=-rho0/eps0,bnd=bndy,nr0=nr0,nz0=nz0, &
!                        length_r=nr0*dr0,length_z=nz0*dz0, &
!                         accuracy=accuracy,sub_accuracy=mgridrz_sub_accuracy, &
!                        nc=ncmax,npre=npre,npost=npost,ncycle=ncycle,nrecurs_min=mgridrz_nrecurs_min)
 call solve_multigridrz_jlv2(u=u,rhoinit=-rho0/eps0,bnd=bndy,nr0=nr0,nz0=nz0, &
                        length_r=nr0*dr0,length_z=nz0*dz0, &
                         accuracy=accuracy,sub_accuracy=mgridrz_sub_accuracy, &
                        nc=ncmax,npre=npre,npost=npost,ncycle=ncycle,nrecurs_min=mgridrz_nrecurs_min)

  u0(1:nr0+1,:)=u(1:nr0+1,:)

return
end subroutine multigridrzf

subroutine srfrvoutrz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                      xmin,xmax,ymin,ymax,lshell,                      &
                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,             &
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
! call subroutine srfrvout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use PSOR3d
USE multigridrz
implicit none
character(*) rofzfunc
real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
LOGICAL(ISZ):: lfill,lshell,l2symtry,l4symtry
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)

INTEGER(ISZ) :: i,nrc,nzc
REAL(8) :: drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(.not.bndy_allocated) call init_bnd(nx,nz,dx,dz)

do i = nlevels,1,-1
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef MPIPARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
  zmax_in = zmin_in + bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
  zmax_in = zmmax
#endif

  call srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                      xmin,xmax,ymin,ymax,lshell,                      &
                      zmmin,zmmax,zbeam,drc,drc,dzc,nrc,ny,nzc,             &
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

  call addconductors_rz(i,nrc,nzc,drc,dzc)

end do
end subroutine srfrvoutrz

subroutine srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,   &
                           lzend,xmin,xmax,ymin,ymax,lshell,          &
                           zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,       &
                           ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use PSOR3d
USE multigridrz
implicit none
character(*) rminofz,rmaxofz
real(kind=8):: volt,zmin,zmax,xcent,ycent
LOGICAL(ISZ):: lzend,lshell,l2symtry,l4symtry
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)

INTEGER(ISZ) :: i,nrc,nzc
REAL(8) :: drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(.not.bndy_allocated) call init_bnd(nx,nz,dx,dz)

do i = nlevels,1,-1
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef MPIPARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
  zmax_in = zmin_in + bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
  zmax_in = zmmax
#endif

  call srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,  &
                     lzend,xmin,xmax,ymin,ymax,lshell,                &
                     zmin_in,zmax_in,zbeam,drc,drc,dzc,nrc,ny,nzc,             &
                     ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

  call addconductors_rz(i,nrc,nzc,drc,dzc)

end do
end subroutine srfrvinoutrz

subroutine setcndtrrz(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz, &
                      l2symtry,l4symtry)
use PSOR3d
USE multigridrz
integer(ISZ):: nx,ny,nz
real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
logical(ISZ):: l2symtry,l4symtry

INTEGER(ISZ) :: i,nrc,nzc
REAL(8) :: drc,dzc,zmin_in

IF(.not.bndy_allocated) call init_bnd(nx,nz,dx,dz)

do i = nlevels,1,-1
  level = i
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef MPIPARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
#endif

  call setcndtr_rz(xmmin,ymmin,zmin_in,zbeam,zgrid,nrc,ny,nzc,drc,drc,dzc, &
                   l2symtry,l4symtry)

  call addconductors_rz(i,nrc,nzc,drc,dzc)

end do

END subroutine setcndtrrz

!subroutine add_egun_data(volt,npts,aj,al,ax,az)
!implicit none

!INTEGER(ISZ), INTENT(IN) :: npts, aj(npts), al(npts)
!REAL(8), INTENT(IN) :: volt, ax(npts), az(npts)

!INTEGER(ISZ) :: i

!do i = 1, npts
!  IF(ax(i)>0.) then
!    delmx=ax(i)
!    delpx=2.0
!  else
!    delpx=ax(i)
!    delmx=2.0
!  END if
!  IF(az(i)>0.) then
!    delmz=az(i)
!    delpz=2.0
!  else
!    delpz=az(i)
!    delmz=2.0
!  END if
!  if (mod(aj(i)+al(i)+zparity,2) == 0) then
!!   --- even
!    necndbdy = necndbdy + 1
!    iecndx(necndbdy) = ix
!    iecndy(necndbdy) = iy
!    iecndz(necndbdy) = iz
!    ecdelmx(necndbdy) = delmx
!    ecdelmy(necndbdy) = delmy
!    ecdelmz(necndbdy) = delmz
!    ecdelpx(necndbdy) = delpx
!    ecdelpy(necndbdy) = delpy
!    ecdelpz(necndbdy) = delpz
!    ecvolt(necndbdy) = volt
!  else
!!   --- odd
!    nocndbdy = nocndbdy + 1
!    iocndx(nocndbdy) = ix
!    iocndy(nocndbdy) = iy
!    iocndz(nocndbdy) = iz
!    ocdelmx(nocndbdy) = delmx
!    ocdelmy(nocndbdy) = delmy
!    ocdelmz(nocndbdy) = delmz
!    ocdelpx(nocndbdy) = delpx
!    ocdelpy(nocndbdy) = delpy
!    ocdelpz(nocndbdy) = delpz
!    ocvolt(nocndbdy) = volt
!  endif
!end do
!
!
!return
!end subroutine add_egun_data

subroutine addconductors_rz(i,nrc,nzc,drc,dzc)
use PSOR3d
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: nrc,nzc,i
REAL(8), INTENT(IN) :: drc,dzc

INTEGER(ISZ) :: ii,iii,iv,iiv,nxbnd,nzbnd
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz

TYPE(conductor_type), POINTER :: cndpnt

  nxbnd=nrc+1
  nzbnd=nzc+1

  call init_bnd_sublevel(bndy(i),necndbdy+nocndbdy,ncond)

  bndy(i)%cnd%nbbndred = necndbdy
  bndy(i)%cnd%nbbnd = necndbdy + nocndbdy

  bndy(i)%cnd%ncond = ncond
  iii=0
  do ii=1,ncond
    iii = iii + 1
    bndy(i)%cnd%jcond(iii) = ixcond(ii)+1
    bndy(i)%cnd%kcond(iii) = izcond(ii)+1
    IF(bndy(i)%v(bndy(i)%cnd%jcond(iii),bndy(i)%cnd%kcond(iii))==v_cond) then
      iii = iii-1
      bndy(i)%cnd%ncond = bndy(i)%cnd%ncond-1
      cycle
    END if
    bndy(i)%v(bndy(i)%cnd%jcond(iii),bndy(i)%cnd%kcond(iii)) = v_cond
    bndy(i)%cnd%voltage(iii) = condvolt(ii)
  end do

  do ii = 1, necndbdy
   bndy(i)%cnd%jj(ii)  = iecndx(ii)+1
   bndy(i)%cnd%kk(ii)  = iecndz(ii)+1
   IF(bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii))==v_cond) cycle
   IF(iecndx(ii)>=0 .and. iecndx(ii)<=nxbnd-1 .and. iecndz(ii)>=0 .and. iecndz(ii)<=nzbnd-1) then
     bndy(i)%cnd%docalc(ii)=.true.
   else
     cycle
   END if
   dxm = MIN(1._8,ecdelmx(ii))*bndy(i)%dr
   dxp = MIN(1._8,ecdelpx(ii))*bndy(i)%dr
   dzm = MIN(1._8,ecdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ecdelpz(ii))*bndy(i)%dz
   bndy(i)%cnd%volt0xm(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0xp(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0zm(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0zp(ii)=ecvolt(ii)
   IF(bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii))/=v_bnd ) then
     bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii)) = v_bnd
   else
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=1,cndpnt%nbbndred
         IF(bndy(i)%cnd%jj(ii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(ii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)<dxm) then
             dxm = cndpnt%dxm(iiv)
             bndy(i)%cnd%volt0xm(ii) = cndpnt%volt0xm(iiv)
           END if
           IF(cndpnt%dxp(iiv)<dxp) then
             dxp = cndpnt%dxp(iiv)
             bndy(i)%cnd%volt0xp(ii) = cndpnt%volt0xp(iiv)
           END if
           IF(cndpnt%dzm(iiv)<dzm) then
             dzm = cndpnt%dzm(iiv)
             bndy(i)%cnd%volt0zm(ii) = cndpnt%volt0zm(iiv)
           END if
           IF(cndpnt%dzp(iiv)<dzp) then
             dzp = cndpnt%dzp(iiv)
             bndy(i)%cnd%volt0zp(ii) = cndpnt%volt0zp(iiv)
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(ii)=dxm
   bndy(i)%cnd%dxp(ii)=dxp
   bndy(i)%cnd%dzm(ii)=dzm
   bndy(i)%cnd%dzp(ii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dr
       dzz=bndy(i)%dz
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(bndy(i)%cnd%jj(ii)==1) then
     rp = 0.5_8*bndy(i)%dr
     dt = coefrelaxbnd/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(ii) = 4._8*dt/(dxp*dxx)
     bndy(i)%cnd%cfzm(ii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfzp(ii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(ii) = -dt
     bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfzm(ii)-bndy(i)%cnd%cfzp(ii)
     bndy(i)%cnd%rcfxp(ii) = -4._8/(dxp*dxx)
     bndy(i)%cnd%rcfzm(ii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfzp(ii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(ii) = 1._8
     bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxp(ii)-bndy(i)%cnd%rcfzm(ii)-bndy(i)%cnd%rcfzp(ii)
   else
     r = (bndy(i)%cnd%jj(ii)-1)*bndy(i)%dr
     rm = r-0.5_8*dxx
     rp = r+0.5_8*dxx
     dt = coefrelaxbnd/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = dt*rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(ii) = dt*rp/(r*dxp*dxx)
     bndy(i)%cnd%cfzm(ii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfzp(ii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(ii) = -dt
     bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxm(ii)-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfzm(ii)-bndy(i)%cnd%cfzp(ii)
     bndy(i)%cnd%rcfxm(ii) = -1._8*rm/(r*dxm*dxx)
     bndy(i)%cnd%rcfxp(ii) = -1._8*rp/(r*dxp*dxx)
     bndy(i)%cnd%rcfzm(ii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfzp(ii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(ii) = 1._8
     bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxm(ii)-bndy(i)%cnd%rcfxp(ii)-bndy(i)%cnd%rcfzm(ii)-bndy(i)%cnd%rcfzp(ii)
   END if
     IF(dxm>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xm(ii)=0._8
       bndy(i)%cnd%rphi0xm(ii)=0._8
     else
       bndy(i)%cnd%phi0xm(ii)=bndy(i)%cnd%cfxm(ii)*bndy(i)%cnd%volt0xm(ii)
       bndy(i)%cnd%cfxm(ii)=0._8
       bndy(i)%cnd%rphi0xm(ii)=bndy(i)%cnd%rcfxm(ii)*bndy(i)%cnd%volt0xm(ii)
       bndy(i)%cnd%rcfxm(ii)=0._8
     END if
     IF(dxp>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xp(ii)=0._8
       bndy(i)%cnd%rphi0xp(ii)=0._8
     else
       bndy(i)%cnd%phi0xp(ii)=bndy(i)%cnd%cfxp(ii)*bndy(i)%cnd%volt0xp(ii)
       bndy(i)%cnd%cfxp(ii)=0._8
       bndy(i)%cnd%rphi0xp(ii)=bndy(i)%cnd%rcfxp(ii)*bndy(i)%cnd%volt0xp(ii)
       bndy(i)%cnd%rcfxp(ii)=0._8
     END if
     IF(dzm>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zm(ii)=0._8
       bndy(i)%cnd%rphi0zm(ii)=0._8
     else
       bndy(i)%cnd%phi0zm(ii)=bndy(i)%cnd%cfzm(ii)*bndy(i)%cnd%volt0zm(ii)
       bndy(i)%cnd%cfzm(ii)=0._8
       bndy(i)%cnd%rphi0zm(ii)=bndy(i)%cnd%rcfzm(ii)*bndy(i)%cnd%volt0zm(ii)
       bndy(i)%cnd%rcfzm(ii)=0._8
     END if
     IF(dzp>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zp(ii)=0._8
       bndy(i)%cnd%rphi0zp(ii)=0._8
     else
       bndy(i)%cnd%phi0zp(ii)=bndy(i)%cnd%cfzp(ii)*bndy(i)%cnd%volt0zp(ii)
       bndy(i)%cnd%cfzp(ii)=0._8
       bndy(i)%cnd%rphi0zp(ii)=bndy(i)%cnd%rcfzp(ii)*bndy(i)%cnd%volt0zp(ii)
       bndy(i)%cnd%rcfzp(ii)=0._8
     END if
  end do


  do ii = 1, nocndbdy
   iii=necndbdy+ii
   bndy(i)%cnd%jj(iii)  = iocndx(ii)+1
   bndy(i)%cnd%kk(iii)  = iocndz(ii)+1
   IF(bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii))==v_cond) cycle
   IF(iocndx(ii)>=0 .and. iocndx(ii)<=nxbnd-1 .and. iocndz(ii)>=0 .and. iocndz(ii)<=nzbnd-1) then
     bndy(i)%cnd%docalc(iii)=.true.
   else
     cycle
   endif
   dxm = MIN(1._8,ocdelmx(ii))*bndy(i)%dr
   dxp = MIN(1._8,ocdelpx(ii))*bndy(i)%dr
   dzm = MIN(1._8,ocdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ocdelpz(ii))*bndy(i)%dz
   bndy(i)%cnd%volt0xm(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0xp(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0zm(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0zp(iii)=ocvolt(ii)
   IF(bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii))/=v_bnd ) then
     bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) = v_bnd
   else
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=cndpnt%nbbndred+1,cndpnt%nbbnd
         IF(bndy(i)%cnd%jj(iii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(iii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)<dxm) then
             dxm = cndpnt%dxm(iiv)
             bndy(i)%cnd%volt0xm(iii) = cndpnt%volt0xm(iiv)
           END if
           IF(cndpnt%dxp(iiv)<dxp) then
             dxp = cndpnt%dxp(iiv)
             bndy(i)%cnd%volt0xp(iii) = cndpnt%volt0xp(iiv)
           END if
           IF(cndpnt%dzm(iiv)<dzm) then
             dzm = cndpnt%dzm(iiv)
             bndy(i)%cnd%volt0zm(iii) = cndpnt%volt0zm(iiv)
           END if
           IF(cndpnt%dzp(iiv)<dzp) then
             dzp = cndpnt%dzp(iiv)
             bndy(i)%cnd%volt0zp(iii) = cndpnt%volt0zp(iiv)
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(iii)=dxm
   bndy(i)%cnd%dxp(iii)=dxp
   bndy(i)%cnd%dzm(iii)=dzm
   bndy(i)%cnd%dzp(iii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dr
       dzz=bndy(i)%dz
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(bndy(i)%cnd%jj(iii)==1) then
     rp = 0.5_8*bndy(i)%dr
     dt = coefrelaxbnd/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(iii) = 4._8*dt/(dxp*dxx)
     bndy(i)%cnd%cfzm(iii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfzp(iii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(iii) = -dt
     bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfzm(iii)-bndy(i)%cnd%cfzp(iii)
     bndy(i)%cnd%rcfxp(iii) = -4._8/(dxp*dxx)
     bndy(i)%cnd%rcfzm(iii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfzp(iii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(iii) = 1._8
     bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxp(iii)-bndy(i)%cnd%rcfzm(iii)-bndy(i)%cnd%rcfzp(iii)
   else
     r = (bndy(i)%cnd%jj(iii)-1)*bndy(i)%dr
     rm = r-0.5_8*dxm
     rp = r+0.5_8*dxp
     dt = coefrelaxbnd/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(iii) = dt*rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(iii) = dt*rp/(r*dxp*dxx)
     bndy(i)%cnd%cfzm(iii) = dt/(dzm*dzz)
     bndy(i)%cnd%cfzp(iii) = dt/(dzp*dzz)
     bndy(i)%cnd%cfrhs(iii) = -dt
     bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxm(iii)-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfzm(iii)-bndy(i)%cnd%cfzp(iii)
     bndy(i)%cnd%rcfxm(iii) = -1._8*rm/(r*dxm*dxx)
     bndy(i)%cnd%rcfxp(iii) = -1._8*rp/(r*dxp*dxx)
     bndy(i)%cnd%rcfzm(iii) = -1._8/(dzm*dzz)
     bndy(i)%cnd%rcfzp(iii) = -1._8/(dzp*dzz)
     bndy(i)%cnd%rcfrhs(iii) = 1._8
     bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxm(iii)-bndy(i)%cnd%rcfxp(iii)-bndy(i)%cnd%rcfzm(iii)-bndy(i)%cnd%rcfzp(iii)
   END if
     IF(dxm>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=0._8
     else
       bndy(i)%cnd%phi0xm(iii)=bndy(i)%cnd%cfxm(iii)*bndy(i)%cnd%volt0xm(iii)
       bndy(i)%cnd%cfxm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=bndy(i)%cnd%rcfxm(iii)*bndy(i)%cnd%volt0xm(iii)
       bndy(i)%cnd%rcfxm(iii)=0._8
     END if
     IF(dxp>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=0._8
     else
       bndy(i)%cnd%phi0xp(iii)=bndy(i)%cnd%cfxp(iii)*bndy(i)%cnd%volt0xp(iii)
       bndy(i)%cnd%cfxp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=bndy(i)%cnd%rcfxp(iii)*bndy(i)%cnd%volt0xp(iii)
       bndy(i)%cnd%rcfxp(iii)=0._8
     END if
     IF(dzm>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=0._8
     else
       bndy(i)%cnd%phi0zm(iii)=bndy(i)%cnd%cfzm(iii)*bndy(i)%cnd%volt0zm(iii)
       bndy(i)%cnd%cfzm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=bndy(i)%cnd%rcfzm(iii)*bndy(i)%cnd%volt0zm(iii)
       bndy(i)%cnd%rcfzm(iii)=0._8
     END if
     IF(dzp>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=0._8
     else
       bndy(i)%cnd%phi0zp(iii)=bndy(i)%cnd%cfzp(iii)*bndy(i)%cnd%volt0zp(iii)
       bndy(i)%cnd%cfzp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=bndy(i)%cnd%rcfzp(iii)*bndy(i)%cnd%volt0zp(iii)
       bndy(i)%cnd%rcfzp(iii)=0._8
     END if
  end do

end subroutine addconductors_rz

!     ******************************************************************
!     *
!     *                        SUBROUTINE RHOWEIGHTRZ
!     *
!     ******************************************************************


subroutine rhoweightrz(xp,yp,zp,np,q,rho,nr,nz,dr,dz,zmin)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(0:nr,0:nz), INTENT(INOUT) :: rho
REAL(8), INTENT(IN) :: q, dr, dz, zmin

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr)
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp

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

!     ******************************************************************
!     *
!     *                        SUBROUTINE RHOWEIGHTRZ_DEFORM
!     *
!     ******************************************************************


subroutine rhoweightrz_deform(xp,yp,zp,np,q,rho,nr,nz,dr,dz,zmin,xfact,yfact)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(0:nr,0:nz), INTENT(INOUT) :: rho
REAL(8), INTENT(IN) :: q, dr, dz, zmin
REAL(8), DIMENSION(0:nz), INTENT(INOUT) :: xfact,yfact

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr)
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp

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
    zpos = (zp(i)-zmin)*invdz
    ln = INT(zpos)
    ddz = zpos-REAL(ln)
    oddz = 1._8-ddz
    rpos = SQRT(xfact(ln)*xp(i)*xfact(ln)*xp(i)+yfact(ln)*yp(i)*yfact(ln)*yp(i))*invdr
    jn = INT(rpos)
    IF(jn>=nr.or.ln<0.or.ln>=nz) cycle
    ddr = rpos-REAL(jn)
    oddr = 1._8-ddr
    jnp=jn+1
    lnp=ln+1
    rho(jn, ln)  = rho(jn, ln)  + q * oddr * oddz * invvol(jn)
    rho(jnp,ln)  = rho(jnp,ln)  + q *  ddr * oddz * invvol(jnp)
    rho(jn, lnp) = rho(jn, lnp) + q * oddr *  ddz * invvol(jn)
    rho(jnp,lnp) = rho(jnp,lnp) + q *  ddr *  ddz * invvol(jnp)
  end do

  return
END SUBROUTINE RHOWEIGHTRZ_DEFORM

subroutine rhobndrz(rho,nr,nz)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), dimension(0:nr,0:nz), INTENT(IN OUT) :: rho

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

subroutine fieldweightrz(xp,yp,zp,ex,ey,ez,np,phi,e,nr,nz,dr,dz,zmin,calcselfe)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), DIMENSION(1:2,0:nr,0:nz), INTENT(IN OUT) :: e
REAL(8), INTENT(IN) :: dr, dz, zmin
LOGICAL(ISZ), INTENT(IN) :: calcselfe

REAL(8) :: invdr, invdz, rpos, zpos, invrpos, ddr, ddz, oddr, oddz, er
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp

 IF(calcselfe) then
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
  END if

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

subroutine fieldweightrz_deform(xp,yp,zp,ex,ey,ez,np,phi,nr,nz,dr,dz,zmin,xfact,yfact,calcphi,phi3d,selfe)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), INTENT(IN) :: dr, dz, zmin
REAL(8), DIMENSION(0:nz), INTENT(INOUT) :: xfact,yfact
LOGICAL(ISZ) :: calcphi
REAL(8), INTENT(INOUT) :: phi3d(0:nr,0:nr,-1:nz+1), selfe(3,0:nr,0:nr,0:nz)

REAL(8) :: invdx, invdy, invdr, invdz, xpos, ypos, rpos, zpos, ddx, ddy, ddr, ddz, oddx, oddy, oddr, oddz, &
           w1, w2, w3, w4, w5, w6, w7, w8
INTEGER(ISZ) :: i, j, k, l, jn, kn, ln, jnp, knp, lnp


IF(calcphi) then
  phi3d = 0.
  selfe = 0.
  do l = -1, nz+1
    do k = 0, nr
      do j = 0, nr
        IF(l<=1.or.l==nz+1) then
          rpos = SQRT(REAL(j)**2+REAL(k)**2)
        else
          rpos = SQRT(j**2*xfact(l)+k**2*yfact(l))
        END if
        jn = INT(rpos)
        IF(jn>=nr) cycle
        ddr=rpos-REAL(jn)
        oddr=1._8-ddr
        phi3d(j,k,l) = oddr*phi(jn,l)+ddr*phi(jn+1,l)
      end do
    end do
  end do

  call getselfe3d(phi3d,nr,nr,nz,selfe,nr,nr,nz,dr,dr,dz)
!WRITE(0,*) 'sum(phi)',SUM(ABS(phi3d)),SUM(ABS(phi))
endif

  invdx = 1._8/dr
  invdy = 1._8/dr
  invdz = 1._8/dz

  ! make field deposition using CIC weighting
  do i = 1, np
    xpos = abs(xp(i))*invdx
    jn = INT(xpos)
    jnp = jn+1
    ddx = xpos-REAL(jn)
    oddx = 1._8-ddx

    ypos = abs(yp(i))*invdy
    kn = INT(ypos)
    knp = kn+1
    ddy = ypos-REAL(kn)
    oddy = 1._8-ddy

    zpos = (zp(i)-zmin)*invdz
    ln = INT(zpos)
    lnp = ln+1
    ddz = zpos-REAL(ln)
    oddz = 1._8-ddz

    w1 = oddx * oddy * oddz
    w2 = ddx  * oddy * oddz
    w3 = oddx * ddy  * oddz
    w4 = ddx  * ddy  * oddz
    w5 = oddx * oddy * ddz
    w6 = ddx  * oddy * ddz
    w7 = oddx * ddy  * ddz
    w8 = ddx  * ddy  * ddz 

    ex(i) = w1 * selfe(1,jn ,kn, ln)  &
          + w2 * selfe(1,jnp,kn, ln)  &
          + w3 * selfe(1,jn ,knp,ln)  &
          + w4 * selfe(1,jnp,knp,ln)  &
          + w5 * selfe(1,jn ,kn, lnp) &
          + w6 * selfe(1,jnp,kn, lnp) &
          + w7 * selfe(1,jn ,knp,lnp) &
          + w8 * selfe(1,jnp,knp,lnp)

    ey(i) = w1 * selfe(2,jn ,kn, ln)  &
          + w2 * selfe(2,jnp,kn, ln)  &
          + w3 * selfe(2,jn ,knp,ln)  &
          + w4 * selfe(2,jnp,knp,ln)  &
          + w5 * selfe(2,jn ,kn, lnp) &
          + w6 * selfe(2,jnp,kn, lnp) &
          + w7 * selfe(2,jn ,knp,lnp) &
          + w8 * selfe(2,jnp,knp,lnp)

    ez(i) = w1 * selfe(3,jn ,kn, ln)  &
          + w2 * selfe(3,jnp,kn, ln)  &
          + w3 * selfe(3,jn ,knp,ln)  &
          + w4 * selfe(3,jnp,knp,ln)  &
          + w5 * selfe(3,jn ,kn, lnp) &
          + w6 * selfe(3,jnp,kn, lnp) &
          + w7 * selfe(3,jn ,knp,lnp) &
          + w8 * selfe(3,jnp,knp,lnp)

    IF(xp(i)<0.) ex(i)=-ex(i)
    IF(yp(i)<0.) ey(i)=-ey(i)

  END do

  return
end subroutine fieldweightrz_deform

subroutine fieldweightrz_deform_old(xp,yp,zp,ex,ey,ez,np,phi,e,nr,nz,dr,dz,zmin,xfact,yfact,calcselfe)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), DIMENSION(1:2,0:nr,0:nz), INTENT(IN OUT) :: e
REAL(8), INTENT(IN) :: dr, dz, zmin
REAL(8), DIMENSION(0:nz), INTENT(INOUT) :: xfact,yfact
LOGICAL(ISZ), INTENT(IN) :: calcselfe

REAL(8) :: invdr, invdz, rpos, zpos, invrpos, ddr, ddz, oddr, oddz, er
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp

 IF(calcselfe) then
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
 END if

  invdr = 1._8/dr
  invdz = 1._8/dz

  ! make field deposition using CIC weighting
  do i = 1, np
    zpos = (zp(i)-zmin)*invdz
    ln = INT(zpos)
    ddz = zpos-REAL(ln)
    oddz = 1._8-ddz
    rpos = SQRT(xfact(ln)*xp(i)*xfact(ln)*xp(i)+yfact(ln)*yp(i)*yfact(ln)*yp(i))*invdr
    jn = INT(rpos)
    ddr = rpos-REAL(jn)
    oddr = 1._8-ddr
    jnp=jn+1
    lnp=ln+1
    er = oddr * oddz * e(1,jn ,ln)  &
       + ddr  * oddz * e(1,jnp,ln)  &
       + oddr * ddz  * e(1,jn ,lnp) &
       + ddr  * ddz  * e(1,jnp,lnp)
    IF(rpos>1.e-10) then
      invrpos=invdr/rpos
      ex(i) = er*xp(i)*xfact(ln)*invrpos
      ey(i) = er*yp(i)*yfact(ln)*invrpos
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
end subroutine fieldweightrz_deform_old

!subroutine calcfact_deform(xp,yp,zp,np,dz,zmin,xfact,yfact,nz,ns,is,ins,nps,ws)
subroutine calcfact_deform(dz,zmin,xfact,yfact,nz,ns,is,ins,nps,ws)
USE constant
USE Particles, ONLY: xp, yp, zp
implicit none

!INTEGER(ISZ), INTENT(IN) :: np, nz, ns
INTEGER(ISZ), INTENT(IN) :: nz, ns
!REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), INTENT(IN) :: dz, zmin, ws(ns)
REAL(8), DIMENSION(0:nz), INTENT(INOUT) :: xfact, yfact
INTEGER(ISZ), DIMENSION(ns), INTENT(IN) :: is, nps, ins

REAL(8), DIMENSION(0:nz) :: fweight, fnb
INTEGER(ISZ) :: i, ln, lnp, isp
REAL(8) :: ddz, oddz, wddz, woddz, xrms, yrms, invdz, zpos, xp2, yp2

  xfact = 0._8
  yfact = 0._8
  fweight = 0._8
  fnb = 0._8
  invdz = 1._8/dz

!WRITE(0,*) 'calcfact',dz,zmin,xfact,yfact,nz,ns,is,ins,nps,ws!isp,ns,ins(is(1)),nps(is(1)),ws(is(1))
!WRITE(0,*) 'calcfact',dz,zmin,nz,ns,is,ins,nps,ws!isp,ns,ins(is(1)),nps(is(1)),ws(is(1))
!return
  do isp = 1, ns
    do i = ins(is(isp)), ins(is(isp))+nps(is(isp))-1
      IF(zp(i)<zmin) cycle
      zpos = (zp(i)-zmin)*invdz
      ln = INT(zpos)
      lnp = ln+1
      ddz = zpos-REAL(ln)
      oddz = 1._8-ddz
      wddz = ws(is(isp))*ddz
      woddz = ws(is(isp))*oddz
      xp2 = xp(i)**2
      yp2 = yp(i)**2
      xfact(ln) = xfact(ln) + xp2 * woddz
      yfact(ln) = yfact(ln) + yp2 * woddz
      xfact(lnp) = xfact(lnp) + xp2 * wddz
      yfact(lnp) = yfact(lnp) + yp2 * wddz
      fweight(ln) = fweight(ln) + woddz
      fweight(lnp) = fweight(lnp) + wddz
      fnb(ln) = fnb(ln) + oddz
      fnb(lnp) = fnb(lnp) + ddz
    end do
  end do
 ! WRITE(0,*) 'sum ',SUM(fweight)/ws(is(1))
  do ln = 0, nz
    IF(fnb(ln)>25._8) then
      xrms = SQRT(xfact(ln)/fweight(ln))
      yrms = SQRT(yfact(ln)/fweight(ln))
 !     WRITE(0,*) 'rms',xrms,yrms
      xfact(ln) = 0.5_8*(xrms+yrms)/xrms
      yfact(ln) = 0.5_8*(xrms+yrms)/yrms
    else
      xfact(ln) = 1._8
      yfact(ln) = 1._8
    END if
!    WRITE(0,*) ln,xfact(ln),yfact(ln),xrms,yrms,fweight(ln),fnb(ln)
  end do

  return
end subroutine calcfact_deform

subroutine save_bndstructure_rz(filename)
use multigridrz
implicit none
CHARACTER(*) :: filename
INTEGER(ISZ) :: i,ic

  OPEN(10,FILE=filename,STATUS='unknown')
    WRITE(10,*) nlevels
    do i = 1, nlevels
      WRITE(10,*) bndy(i)%dr, bndy(i)%dz, bndy(i)%nb_conductors, bndy(i)%nr, bndy(i)%nz, bndy(i)%l_powerof2
      WRITE(10,*) bndy(i)%v
      do ic = 1, bndy(i)%nb_conductors
        IF(ic==1) then
          bndy(i)%cnd => bndy(i)%first
        else
          bndy(i)%cnd => bndy(i)%cnd%next
        END if
          WRITE(0,*) i,ic,nlevels,bndy(i)%cnd%nbbnd,bndy(i)%cnd%ncond, bndy(i)%cnd%nbbndred
          WRITE(10,*) bndy(i)%cnd%nbbnd,bndy(i)%cnd%ncond
          WRITE(10,*) bndy(i)%cnd%nbbndred
          WRITE(10,*) bndy(i)%cnd%voltage
        IF(bndy(i)%cnd%nbbnd>0) then
          WRITE(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp,bndy(i)%cnd%phi0zm,bndy(i)%cnd%phi0zp
          WRITE(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, bndy(i)%cnd%cfzp, bndy(i)%cnd%cfzm, bndy(i)%cnd%cfrhs
          WRITE(10,*) bndy(i)%cnd%rphi0xm,bndy(i)%cnd%rphi0xp,bndy(i)%cnd%rphi0zm,bndy(i)%cnd%rphi0zp
          WRITE(10,*) bndy(i)%cnd%rcf0, bndy(i)%cnd%rcfxp, bndy(i)%cnd%rcfxm, &
                      bndy(i)%cnd%rcfzp, bndy(i)%cnd%rcfzm, bndy(i)%cnd%rcfrhs
          WRITE(10,*) bndy(i)%cnd%dxm,bndy(i)%cnd%dxp,bndy(i)%cnd%dzm,bndy(i)%cnd%dzp
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
INTEGER(ISZ) :: nbbnd,ncond,i,ic,nbc

  OPEN(10,FILE=filename,STATUS='unknown')
    read(10,*) nlevels
    ALLOCATE(bndy(nlevels))
    bndy_allocated=.true.
    do i = 1, nlevels
      NULLIFY(bndy(i)%first)
      read(10,*) bndy(i)%dr, bndy(i)%dz, bndy(i)%nb_conductors, bndy(i)%nr, bndy(i)%nz, bndy(i)%l_powerof2
      ALLOCATE(bndy(i)%v(bndy(i)%nr+1,bndy(i)%nz+1))
      read(10,*) bndy(i)%v
      nbc = bndy(i)%nb_conductors
      do ic = 1, nbc
        read(10,*) nbbnd,ncond
        call init_bnd_sublevel(bndy(i),nbbnd,ncond)
        read(10,*) bndy(i)%cnd%nbbndred
        read(10,*) bndy(i)%cnd%voltage
        WRITE(0,*) i,ic,nlevels,nbbnd,ncond,bndy(i)%cnd%nbbndred
        IF(bndy(i)%cnd%nbbnd>0) then
          read(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp,bndy(i)%cnd%phi0zm,bndy(i)%cnd%phi0zp
          read(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, bndy(i)%cnd%cfzp, bndy(i)%cnd%cfzm, bndy(i)%cnd%cfrhs
          read(10,*) bndy(i)%cnd%rphi0xm,bndy(i)%cnd%rphi0xp,bndy(i)%cnd%rphi0zm,bndy(i)%cnd%rphi0zp
          read(10,*) bndy(i)%cnd%rcf0, bndy(i)%cnd%rcfxp, bndy(i)%cnd%rcfxm, &
                     bndy(i)%cnd%rcfzp, bndy(i)%cnd%rcfzm, bndy(i)%cnd%rcfrhs
          read(10,*) bndy(i)%cnd%dxm,bndy(i)%cnd%dxp,bndy(i)%cnd%dzm,bndy(i)%cnd%dzp
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

module multigrid3d_jlv

USE multigrid_common

TYPE bndptr
! type for potential calculation
  ! mask array on the grid, = .true. on nodes in vacuum (i.e. not near or in conductors)
  INTEGER(ISZ), POINTER :: v(:,:,:)
  LOGICAL(ISZ) :: l_powerof2 ! set to true if all grid dimensions (in number of meshes) are a power of two
  LOGICAL(ISZ) :: l_merged ! set to true if grids between tasks have been merged
  INTEGER(ISZ) :: nx, ny, nz
  ! mesh sizes
  REAL(8) :: dx, dy, dz
  ! number of conductors
  INTEGER(ISZ) :: nb_conductors
  ! linked-list of conductors
  TYPE(conductor_type), POINTER :: cnd,first
  INTEGER(ISZ) :: izlbnd, izrbnd  ! boundary conditions for each side
#ifdef MPIPARALLEL
  INTEGER(ISZ) :: nworkpproc
#endif
END TYPE bndptr

TYPE(bndptr), pointer :: bndy(:)

contains

subroutine init_bnd(nx,ny,nz,dx,dy,dz)
! intializes grid quantities according to the number of multigrid levels and grid sizes nx and nz.
USE InGen3d, ONLY:l2symtry, l4symtry
implicit none
INTEGER(ISZ), INTENT(IN) :: nx, ny, nz
REAL(8), INTENT(IN) :: dx, dy, dz

INTEGER(ISZ) :: i, nxp0, nyp0, nzp0, nxc, nyc, nzc
REAL(8) :: dxc, dyc, dzc

  IF(l_mgridrz_debug) WRITE(0,*) 'Enter init_bnd'

#ifdef MPIPARALLEL
  nzfine = nz
#endif

  ixlbnd = boundxy
  ixrbnd = boundxy
  iylbnd = boundxy
  iyrbnd = boundxy
  izlbnd = bound0
  izrbnd = boundnz

  IF(l2symtry) then
    ixlbnd = neumann
  END if
  IF(l4symtry) then
    ixlbnd = neumann
    iylbnd = neumann
  END if

  nxc = nx
  nyc = ny
  nzc = nz
  dxc = dx
  dyc = dy
  dzc = dz
  nlevels    = 1
#ifdef MPIPARALLEL
  nworkpproc = 1
#endif
  do WHILE(nxc>mgridrz_nmeshmin.or.nyc>mgridrz_nmeshmin.or.nzc>mgridrz_nmeshmin)
    call evalnewgrid(nxc,nyc,nzc,dxc,dyc,dzc)
    nlevels = nlevels + 1
    WRITE(0,*) nlevels,nxc,nyc,nzc,dxc,dyc,dzc
#ifdef MPIPARALLEL
    IF(nworkpproc*nxc*nyc*nzc<=workfact*nzfine) then
      nlevels = nlevels + 1
      nworkpproc = nworkpproc*2
    END if
#endif
  end do

  nlevels=MIN(nlevels,mgridrz_nlevels_max)

  allocate(bndy(nlevels))
  bndy_allocated=.true.

  nxc = nx
  nyc = ny
  nzc = nz
  dxc = dx
  dyc = dy
  dzc = dz
#ifdef MPIPARALLEL
  nworkpproc = 1
#endif
  do i = nlevels, 1, -1
    bndy(i)%l_merged=.false.
    bndy(i)%izlbnd=izlbnd
    bndy(i)%izrbnd=izrbnd
    IF(i/=nlevels) then
      call evalnewgrid(nxc,nyc,nzc,dxc,dyc,dzc)
#ifdef MPIPARALLEL
      IF(nworkpproc<nslaves.and.nworkpproc*nxc*nyc*nzc<=workfact*nzfine) then
        nworkpproc = nworkpproc*2
        bndy(i)%l_merged=.true.
      END if
#endif
    END if
#ifdef MPIPARALLEL
    bndy(i)%nworkpproc = nworkpproc
    IF(my_index-nworkpproc>=0) bndy(i)%izlbnd = -(my_index-nworkpproc+1)
    IF(my_index+nworkpproc<nslaves) bndy(i)%izrbnd = -(my_index+nworkpproc+1)
    bndy(i)%nx = nxc
    bndy(i)%ny = nyc
    bndy(i)%nz = nworkpproc*nzc
#else
    bndy(i)%nx = nxc
    bndy(i)%ny = nyc
    bndy(i)%nz = nzc
#endif
    bndy(i)%dx = dxc
    bndy(i)%dy = dyc
    bndy(i)%dz = dzc
    NULLIFY(bndy(i)%first)
    bndy(i)%nb_conductors = 0
    ALLOCATE(bndy(i)%v(bndy(i)%nx+1,bndy(i)%ny+1,bndy(i)%nz+1))
    bndy(i)%v(:,:,:)=v_vacuum
  end do
#ifdef MPIPARALLEL
  IF(nslaves>1) then
    do i = 1, nlevels
      bndy(i)%l_powerof2 = .false.
    END do
    do i = nlevels,1,-1
    END do
    return
  END if
#endif
  do i = nlevels, 2, -1
    nxp0 = INT(LOG(REAL(bndy(i)%nx+1))/LOG(2.))
    nyp0 = INT(LOG(REAL(bndy(i)%ny+1))/LOG(2.))
    nzp0 = INT(LOG(REAL(bndy(i)%nz+1))/LOG(2.))
!    nxp0 = LOG(REAL(bndy(i)%nx))/LOG(2.)+0.5
!    nzp0 = LOG(REAL(bndy(i)%nz))/LOG(2.)+0.5
    IF(2**nxp0/=bndy(i)%nx.OR.2**nyp0/=bndy(i)%ny.OR.2**nzp0/=bndy(i)%nz.OR. &
       bndy(i-1)%nx==bndy(i)%nx.OR.bndy(i-1)%ny==bndy(i)%ny.OR.bndy(i-1)%nz==bndy(i)%nz) THEN
      bndy(i)%l_powerof2=.false.
    else
      bndy(i)%l_powerof2=.true.
    END if
    IF(l_mgridrz_debug) WRITE(0,*) i,bndy(i)%l_powerof2,2**nxp0,bndy(i)%nx,2**nyp0,bndy(i)%ny,2**nzp0,bndy(i)%nz
  end do
  bndy(1)%l_powerof2=.false.

  IF(mgridrz_deform) then
    mgridrz_nz = nz
    call gchange("FRZmgrid",0)
    mgridrz_xfact = 1._8
    mgridrz_yfact = 1._8
  END if

  do i = nlevels, 1, -1
    WRITE(0,*) 'level,nx,ny,nz = ',i,bndy(i)%nx,bndy(i)%ny,bndy(i)%nz
  END do

  IF(l_mgridrz_debug) WRITE(0,*) 'Exit init_bnd'

  return
end subroutine init_bnd

subroutine evalnewgrid(nx,ny,nz,dx,dy,dz)
! evaluate nx, ny and nz at coarser level
INTEGER(ISZ), INTENT(IN OUT) :: nx, ny, nz
REAL(8), INTENT(IN OUT) :: dx,dy,dz

REAL(8) :: rapxy, rapxz, rapyz
INTEGER(ISZ) :: nxnew, nynew,nznew
LOGICAL(ISZ) :: l_divx, l_divy, l_divz

  rapxy = dx/dy
  rapxz = dx/dz
  rapyz = dy/dz

  l_divx = .true.
  l_divy = .true.
  l_divz = .true.

  IF(nx<=mgridrz_nmeshmin) l_divx = .false.
  IF(ny<=mgridrz_nmeshmin) l_divy = .false.
  IF(nz<=mgridrz_nmeshmin) l_divz = .false.

  IF(rapxy>4._8/3._8) then
    l_divx = .false.
  ELSE IF(rapxy<2._8/3._8) then
    l_divy = .false.
  END if

  IF(rapxz>4._8/3._8) then
    l_divx = .false.
  ELSE IF(rapxz<2._8/3._8) then
    l_divz = .false.
  END if

  IF(rapyz>4._8/3._8) then
    l_divy = .false.
  ELSE IF(rapyz<2._8/3._8) then
    l_divz = .false.
  END if

  IF(.not.l_divx.and..not.l_divy.AND..not.l_divz) then
    IF(nx>mgridrz_nmeshmin) l_divx = .true.
    IF(ny>mgridrz_nmeshmin) l_divy = .true.
    IF(nz>mgridrz_nmeshmin) l_divz = .true.
  END if

  IF(l_divx) then
    nxnew = MAX(mgridrz_nmeshmin,nx/2)
    dx = dx * REAL(nx,8)/REAL(nxnew,8)
    nx=nxnew
  END if
  IF(l_divy) then
    nynew = MAX(mgridrz_nmeshmin,ny/2)
    dy = dy * REAL(ny,8)/REAL(nynew,8)
    ny=nynew
  END if
  IF(l_divz) then
    nznew = MAX(mgridrz_nmeshmin,nz/2)
    dz = dz * REAL(nz,8)/REAL(nznew,8)
    nz=nznew
  END if

  return
end subroutine evalnewgrid

subroutine init_bnd_sublevel(bndl,nbnd,ncond)
! initializes quantities for one grid level
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
  ALLOCATE(bndl%cnd%phi0xm(nbnd),bndl%cnd%phi0xp(nbnd))
  ALLOCATE(bndl%cnd%phi0ym(nbnd),bndl%cnd%phi0yp(nbnd))
  ALLOCATE(bndl%cnd%phi0zm(nbnd),bndl%cnd%phi0zp(nbnd))
  ALLOCATE(bndl%cnd%volt0xm(nbnd),bndl%cnd%volt0xp(nbnd))
  ALLOCATE(bndl%cnd%volt0ym(nbnd),bndl%cnd%volt0yp(nbnd))
  ALLOCATE(bndl%cnd%volt0zm(nbnd),bndl%cnd%volt0zp(nbnd))
  ALLOCATE(bndl%cnd%cf0(nbnd),bndl%cnd%cfxp(nbnd),bndl%cnd%cfxm(nbnd), &
           bndl%cnd%cfyp(nbnd),bndl%cnd%cfym(nbnd), &
           bndl%cnd%cfzp(nbnd),bndl%cnd%cfzm(nbnd),bndl%cnd%cfrhs(nbnd))
  ALLOCATE(bndl%cnd%rphi0xm(nbnd),bndl%cnd%rphi0xp(nbnd), &
           bndl%cnd%rphi0ym(nbnd),bndl%cnd%rphi0yp(nbnd), &
           bndl%cnd%rphi0zm(nbnd),bndl%cnd%rphi0zp(nbnd))
  ALLOCATE(bndl%cnd%rcf0(nbnd),bndl%cnd%rcfxp(nbnd),bndl%cnd%rcfxm(nbnd), &
           bndl%cnd%rcfyp(nbnd),bndl%cnd%rcfym(nbnd), &
           bndl%cnd%rcfzp(nbnd),bndl%cnd%rcfzm(nbnd),bndl%cnd%rcfrhs(nbnd))
  ALLOCATE(bndl%cnd%dxm(nbnd),bndl%cnd%dxp(nbnd), &
           bndl%cnd%dym(nbnd),bndl%cnd%dyp(nbnd), &
           bndl%cnd%dzm(nbnd),bndl%cnd%dzp(nbnd))
  ALLOCATE(bndl%cnd%jj(nbnd),bndl%cnd%kk(nbnd),bndl%cnd%ll(nbnd),bndl%cnd%docalc(nbnd))
  bndl%cnd%docalc=.false.
END if
IF(ncond>0) then
  ALLOCATE(bndl%cnd%jcond(ncond),bndl%cnd%kcond(ncond),bndl%cnd%lcond(ncond),bndl%cnd%voltage(ncond))
END if

end subroutine init_bnd_sublevel

function expandwguard_any(uold, nxnew, nynew, nznew, xminnew, xmaxnew, yminnew, ymaxnew, zminnew, zmaxnew, &
                          xminold, xmaxold, yminold, ymaxold, zminold, zmaxold)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nynew, nznew
REAL(8), DIMENSION(0:,0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, yminold, ymaxold, zminold, zmaxold, xminnew, xmaxnew, yminnew, ymaxnew, zminnew, zmaxnew
REAL(8) :: expandwguard_any(0:nxnew+2,0:nynew+2,0:nznew+2)

INTEGER(ISZ) :: nxold, nyold, nzold
INTEGER(ISZ) :: jnew, knew, lnew, j, k, l, jp, kp, lp
REAL(8) :: x, y, z, xx, yy, zz, invdxold, invdyold, invdzold, dxnew, dynew, dznew, delx, dely, delz, odelx, odely, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddy(nynew+1), oddy(nynew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nynew+1), lold(nznew+1), joldp(nxnew+1), koldp(nynew+1), loldp(nznew+1)

IF(l_mgridrz_debug) WRITE(0,*) 'enter expand, level = ',level

nxold = SIZE(uold,1) - 1 - 2
nyold = SIZE(uold,2) - 1 - 2
nzold = SIZE(uold,3) - 1 - 2

dxnew = (xmaxnew-xminnew) / nxnew
dynew = (ymaxnew-yminnew) / nynew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdyold = REAL(nyold,8)/(ymaxold-yminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

do lnew = 1, nznew+1
  z = zminnew+(lnew-1)*dznew
  zz = (z-zminold) * invdzold
  lold(lnew) = 1+MIN(nzold,INT(zz))
  loldp(lnew) = lold(lnew) + 1
  ddz(lnew) = zz-(lold(lnew)-1)
  oddz(lnew) = 1.-ddz(lnew)
END do
do knew = 1, nynew+1
  y = yminnew+(knew-1)*dynew
  yy = (y-yminold) * invdyold
  kold(knew) = 1+MIN(nyold,INT(yy))
  koldp(knew) = kold(knew) + 1
  ddy(knew) = yy-(kold(knew)-1)
  oddy(knew) = 1.-ddy(knew)
END do
do jnew = 1, nxnew+1
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
do lnew = 1, nznew+1
  l = lold(lnew)
  lp = loldp(lnew)
  delz = ddz(lnew)
  odelz = oddz(lnew)
  do knew = 1, nynew+1
    k = kold(knew)
    kp = koldp(knew)
    dely = ddy(knew)
    odely = oddy(knew)
    do jnew = 1, nxnew+1
      j = jold(jnew)
      jp = joldp(jnew)
      delx = ddx(jnew)
      odelx = oddx(jnew)
      expandwguard_any(jnew,knew,lnew) = uold(j, k, l)  * odelx * odely * odelz &
                                       + uold(jp,k, l)  * delx  * odely * odelz &
                                       + uold(j, kp,l)  * odelx * dely  * odelz &
                                       + uold(jp,kp,l)  * delx  * dely  * odelz &
                                       + uold(j, k, lp) * odelx * odely * delz &
                                       + uold(jp,k, lp) * delx  * odely * delz &
                                       + uold(j, kp,lp) * odelx * dely  * delz &
                                       + uold(jp,kp,lp) * delx  * dely  * delz
    end do
  end do
END do

expandwguard_any(0,:,:) = 0._8
expandwguard_any(:,0,:) = 0._8
expandwguard_any(:,:,0) = 0._8
expandwguard_any(nxnew+2,:,:) = 0._8
expandwguard_any(:,nynew+2,:) = 0._8
expandwguard_any(:,:,nznew+2) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END function expandwguard_any

function restrict(uold, nxnew, nynew, nznew, xminold, xmaxold, yminold, ymaxold, zminold, zmaxold,  &
                        xminnew, xmaxnew, yminnew, ymaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nynew, nznew
REAL(8), DIMENSION(1:,1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, yminold, ymaxold, zminold, zmaxold, xminnew, xmaxnew, yminnew, ymaxnew, zminnew, zmaxnew
REAL(8) :: restrict(1:nxnew+1,1:nynew+1,1:nznew+1),rap(1:nxnew+1,1:nynew+1,1:nznew+1)

INTEGER(ISZ) :: nxold, nyold, nzold
INTEGER(ISZ) :: jold, kold, lold, j, k, l, jp, kp, lp
REAL(8) :: dxold, dyold, dzold, invdxnew, invdynew, invdznew, x, y, z, delx, dely, delz, odelx, odely, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddy, ddz, oddx, oddy, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, lnew, jnewp, knewp, lnewp

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold)
!  return
!END if

IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict, level = ',level

nxold = SIZE(uold,1) - 1
nyold = SIZE(uold,2) - 1
nzold = SIZE(uold,3) - 1

ALLOCATE(ddx(nxold+1), ddy(nyold+1), ddz(nzold+1), oddx(nxold+1), oddy(nyold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nyold+1), lnew(nzold+1), jnewp(nxold+1), knewp(nyold+1), lnewp(nzold+1))

invdxnew = nxnew / (xmaxnew-xminnew)
invdynew = nynew / (ymaxnew-yminnew)
invdznew = nznew / (zmaxnew-zminnew)
dxold = (xmaxold-xminold) / nxold
dyold = (ymaxold-yminold) / nyold
dzold = (zmaxold-zminold) / nzold

restrict=0._8
rap = 0._8

do lold = 1, nzold+1
  z = zminold + (lold-1)*dzold
  lnew(lold) = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  ddz(lold) = (z-zminnew) * invdznew-real(lnew(lold)-1)
  lnewp(lold) = lnew(lold)+1
  oddz(lold) = 1._8-ddz(lold)
END do

#ifdef MPIPARALLEL
    IF(bndy(level)%izlbnd<0) then
      ddz(1) = 0.5_8*ddz(1)
      oddz(1) = 0.5_8*oddz(1)
    END if
    IF(bndy(level)%izrbnd<0) then
      ddz(nzold+1) = 0.5_8*ddz(nzold+1)
      oddz(nzold+1) = 0.5_8*oddz(nzold+1)
    END if
#endif

do kold = 1, nyold+1
  y = yminold + (kold-1)*dyold
  knew(kold) = MIN(1 + INT((y-yminnew) * invdynew), nynew)
  ddy(kold) = (y-yminnew) * invdynew-real(knew(kold)-1)
  knewp(kold) = knew(kold)+1
  oddy(kold) = 1._8-ddy(kold)
END do
do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = 1._8-ddx(jold)
END do

do lold = 1, nzold+1
  l = lnew(lold)
  lp = lnewp(lold)
  delz  =  ddz(lold)
  odelz = oddz(lold)
  do kold = 1, nyold+1
    k = knew(kold)
    kp = knewp(kold)
    dely  =  ddy(kold)
    odely = oddy(kold)
    do jold = 1, nxold+1
      j = jnew(jold)
      jp = jnewp(jold)
      delx  =  ddx(jold)
      odelx = oddx(jold)
      restrict(j, k, l)  = restrict(j, k, l)  + uold(jold,kold,lold) * odelx * odely * odelz
      restrict(jp,k, l)  = restrict(jp,k, l)  + uold(jold,kold,lold) * delx  * odely * odelz
      restrict(j, kp,l)  = restrict(j, kp,l)  + uold(jold,kold,lold) * odelx * dely  * odelz
      restrict(jp,kp,l)  = restrict(jp,kp,l)  + uold(jold,kold,lold) * delx  * dely  * odelz
      restrict(j, k, lp) = restrict(j, k, lp) + uold(jold,kold,lold) * odelx * odely * delz
      restrict(jp,k, lp) = restrict(jp,k, lp) + uold(jold,kold,lold) * delx  * odely * delz
      restrict(j, kp,lp) = restrict(j, kp,lp) + uold(jold,kold,lold) * odelx * dely  * delz
      restrict(jp,kp,lp) = restrict(jp,kp,lp) + uold(jold,kold,lold) * delx  * dely  * delz
      rap(j, k, l)  = rap(j, k, l)  + odelx * odely * odelz
      rap(jp,k, l)  = rap(jp,k, l)  + delx  * odely * odelz
      rap(j, kp,l)  = rap(j, kp,l)  + odelx * dely  * odelz
      rap(jp,kp,l)  = rap(jp,kp,l)  + delx  * dely  * odelz
      rap(j, k, lp) = rap(j, k, lp) + odelx * odely * delz
      rap(jp,k, lp) = rap(jp,k, lp) + delx  * odely * delz
      rap(j, kp,lp) = rap(j, kp,lp) + odelx * dely  * delz
      rap(jp,kp,lp) = rap(jp,kp,lp) + delx  * dely  * delz  
    end do
  end do
end do

do l = 1, nznew+1
 do k = 1, nynew+1
  do j = 1, nxnew+1
    IF(rap(j,k,l)/=0._8) restrict(j,k,l) = restrict(j,k,l) / rap(j,k,l)
  end do
 end do
end do

DEALLOCATE(ddx, ddy, ddz, oddx, oddy, oddz, jnew, knew, lnew, jnewp, knewp, lnewp)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END function restrict

function restrict_wbnd(uold, bnd, nxnew, nynew, nznew, xminold, xmaxold, yminold, ymaxold, zminold, zmaxold, &
                       xminnew, xmaxnew, yminnew, ymaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nynew, nznew
REAL(8), DIMENSION(1:,1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, yminold, ymaxold, zminold, zmaxold, xminnew, xmaxnew, yminnew, ymaxnew, zminnew, zmaxnew
REAL(8) :: restrict_wbnd(1:nxnew+1,1:nynew+1,1:nznew+1),rap(1:nxnew+1,1:nynew+1,1:nznew+1)
TYPE(bndptr) :: bnd

INTEGER(ISZ) :: nxold, nyold, nzold
INTEGER(ISZ) :: jold, kold, lold, j, k, l, jp, kp, lp, ic, ii, itot
REAL(8) :: dxold, dyold, dzold, invdxnew, invdynew, invdznew, x, y, z, &
           delx, dely, delz, odelx, odely, odelz, dxnew, dynew, dznew, u1, u2, u3, u4, u5, u6, u7, u8, q
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddy, ddz, oddx, oddy, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, lnew, jnewp, knewp, lnewp
LOGICAL(ISZ) :: l_dxm, l_dxp, l_dym, l_dyp, l_dzm, l_dzp
INTEGER(ISZ),allocatable :: cnt(:,:,:)

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold)
!  return
!END if
itot = 0
IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict, level = ',level

nxold = SIZE(uold,1) - 1
nyold = SIZE(uold,2) - 1
nzold = SIZE(uold,3) - 1

ALLOCATE(ddx(nxold+1), ddy(nyold+1), ddz(nzold+1), oddx(nxold+1), oddy(nyold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nyold+1), lnew(nzold+1), jnewp(nxold+1), knewp(nyold+1), lnewp(nzold+1))

ALLOCATE(cnt(nxold+1,nyold+1,nzold+1))
invdxnew = nxnew / (xmaxnew-xminnew)
invdynew = nynew / (ymaxnew-yminnew)
invdznew = nznew / (zmaxnew-zminnew)
dxnew = 1./invdxnew
dynew = 1./invdynew
dznew = 1./invdznew
dxold = (xmaxold-xminold) / nxold
dyold = (ymaxold-yminold) / nyold
dzold = (zmaxold-zminold) / nzold

restrict_wbnd=0._8
rap = 0._8
cnt= 0

do lold = 1, nzold+1
  z = zminold + (lold-1)*dzold
  lnew(lold) = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  ddz(lold) = (z-zminnew) * invdznew-real(lnew(lold)-1)
  lnewp(lold) = lnew(lold)+1
  oddz(lold) = 1._8-ddz(lold)
END do

#ifdef MPIPARALLEL
    IF(bndy(level)%izlbnd<0) then
      ddz(1) = 0.5_8*ddz(1)
      oddz(1) = 0.5_8*oddz(1)
    END if
    IF(bndy(level)%izrbnd<0) then
      ddz(nzold+1) = 0.5_8*ddz(nzold+1)
      oddz(nzold+1) = 0.5_8*oddz(nzold+1)
    END if
#endif

do kold = 1, nyold+1
  y = yminold + (kold-1)*dyold
  knew(kold) = MIN(1 + INT((y-yminnew) * invdynew), nynew)
  ddy(kold) = (y-yminnew) * invdynew-real(knew(kold)-1)
  knewp(kold) = knew(kold)+1
  oddy(kold) = 1._8-ddy(kold)
END do
do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = 1._8-ddx(jold)
END do

do lold = 1, nzold+1
  l = lnew(lold)
  lp = lnewp(lold)
  delz  =  ddz(lold)
  odelz = oddz(lold)
  do kold = 1, nyold+1
    k = knew(kold)
    kp = knewp(kold)
    dely  =  ddy(kold)
    odely = oddy(kold)
    do jold = 1, nxold+1
#ifdef MPIPARALLEL
      IF((bndy(level)%izlbnd<0.and.lold>1).or.(bndy(level)%izrbnd<0.and.lold<nzold+1)) then
        IF(.NOT.(bnd%v(jold,kold,lold)==v_vacuum.or.jold==1.or.jold==nxold+1.or.kold==1.or.kold==nyold+1)) cycle
      END if
#else
!      IF(.NOT.(bnd%v(jold,kold).or.jold==1.or.jold==nxold+1.or.kold==1.or.kold==nzold+1)) cycle
      IF(.NOT.bnd%v(jold,kold,lold)==v_vacuum) cycle
#endif
    j = jnew(jold)
    jp = jnewp(jold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    u1 = odelx * odely * odelz
    u2 = delx  * odely * odelz
    u3 = odelx * dely  * odelz
    u4 = delx  * dely  * odelz
    u5 = odelx * odely * delz
    u6 = delx  * odely * delz
    u7 = odelx * dely  * delz
    u8 = delx  * dely  * delz
    q = uold(jold,kold,lold)
    itot = itot + 1
    cnt(jold,kold,lold) = cnt(jold,kold,lold) + 1
    restrict_wbnd(j, k, l)  = restrict_wbnd(j, k, l)  + q * u1
    restrict_wbnd(jp,k, l)  = restrict_wbnd(jp,k, l)  + q * u2
    restrict_wbnd(j, kp,l)  = restrict_wbnd(j, kp,l)  + q * u3
    restrict_wbnd(jp,kp,l)  = restrict_wbnd(jp,kp,l)  + q * u4
    restrict_wbnd(j, k, lp) = restrict_wbnd(j, k, lp) + q * u5
    restrict_wbnd(jp,k, lp) = restrict_wbnd(jp,k, lp) + q * u6
    restrict_wbnd(j, kp,lp) = restrict_wbnd(j, kp,lp) + q * u7
    restrict_wbnd(jp,kp,lp) = restrict_wbnd(jp,kp,lp) + q * u8
    rap(j, k, l)  = rap(j, k, l)  + u1
    rap(jp,k, l)  = rap(jp,k, l)  + u2
    rap(j, kp,l)  = rap(j, kp,l)  + u3
    rap(jp,kp,l)  = rap(jp,kp,l)  + u4
    rap(j, k, lp) = rap(j, k, lp) + u5
    rap(jp,k, lp) = rap(jp,k, lp) + u6
    rap(j, kp,lp) = rap(j, kp,lp) + u7
    rap(jp,kp,lp) = rap(jp,kp,lp) + u8
  end do
 end do
end do

!GO TO 10
do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do ii = 1, bnd%cnd%nbbnd
    jold = bnd%cnd%jj(ii)
    kold = bnd%cnd%kk(ii)
    lold = bnd%cnd%ll(ii)
    IF(.NOT.(bnd%v(jold,kold,lold)==v_bnd.and.bnd%cnd%docalc(ii))) cycle
    j = jnew(jold)
    jp = jnewp(jold)
    k = knew(kold)
    kp = knewp(kold)
    l = lnew(lold)
    lp = lnewp(lold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    dely  =  ddy(kold)
    odely = oddy(kold)
    delz  =  ddz(kold)
    odelz = oddz(kold)
    u1 = odelx * odely * odelz
    u2 = delx  * odely * odelz
    u3 = odelx * dely  * odelz
    u4 = delx  * dely  * odelz
    u5 = odelx * odely * delz
    u6 = delx  * odely * delz
    u7 = odelx * dely  * delz
    u8 = delx  * dely  * delz
    l_dxm = .true.
    l_dxp = .true.
    l_dym = .true.
    l_dyp = .true.
    l_dzm = .true.
    l_dzp = .true.
!GO TO 10
    IF(bnd%cnd%dxm(ii)<delx*dxnew)  l_dxm = .false.
    IF(bnd%cnd%dxp(ii)<odelx*dxnew) l_dxp = .false.
    IF(bnd%cnd%dym(ii)<dely*dynew)  l_dym = .false.
    IF(bnd%cnd%dyp(ii)<odely*dynew) l_dyp = .false.
    IF(bnd%cnd%dzm(ii)<delz*dznew)  l_dzm = .false.
    IF(bnd%cnd%dzp(ii)<odelz*dznew) l_dzp = .false.
    IF((l_dxm.and.l_dxp).OR.(l_dym.and.l_dyp).OR.(l_dzm.and.l_dzp)) then
      cycle
    else
      IF(l_dxm) then
!        u2=u2+u1
        u1=0
!        u4=u4+u3
        u3=0
!        u6=u6+u5
        u5=0
!        u8=u8+u7
        u7=0
      END if
      IF(l_dxp) then
!        u1=u1+u2
        u2=0
!        u3=u3+u4
        u4=0
!        u5=u5+u6
        u6=0
!        u7=u7+u8
        u8=0
      END if
      IF(l_dym) then
!        u3=u3+u1
        u1=0
!        u4=u4+u2
        u2=0
!        u7=u7+u5
        u5=0
!        u8=u8+u6
        u6=0
      END if
      IF(l_dyp) then
!        u1=u1+u3
        u3=0
!        u2=u2+u4
        u4=0
!        u5=u5+u7
        u7=0
!        u6=u6+u8
        u8=0
      END if
      IF(l_dzm) then
!        u5=u6+u1
        u1=0
!        u6=u6+u2
        u2=0
!        u7=u7+u3
        u3=0
!        u8=u8+u4
        u4=0
      END if
      IF(l_dzp) then
!        u1=u1+u5
        u5=0
!        u2=u2+u6
        u6=0
!        u3=u3+u7
        u7=0
!        u4=u4+u8
        u8=0
      END if
    END if
!10 continue
    q = uold(jold,kold,lold)
    itot = itot + 1
    cnt(jold,kold,lold) = cnt(jold,kold,lold) + 1
    restrict_wbnd(j, k, l)  = restrict_wbnd(j, k, l)  + q * u1
    restrict_wbnd(jp,k, l)  = restrict_wbnd(jp,k, l)  + q * u2
    restrict_wbnd(j, kp,l)  = restrict_wbnd(j, kp,l)  + q * u3
    restrict_wbnd(jp,kp,l)  = restrict_wbnd(jp,kp,l)  + q * u4
    restrict_wbnd(j, k, lp) = restrict_wbnd(j, k, lp) + q * u5
    restrict_wbnd(jp,k, lp) = restrict_wbnd(jp,k, lp) + q * u6
    restrict_wbnd(j, kp,lp) = restrict_wbnd(j, kp,lp) + q * u7
    restrict_wbnd(jp,kp,lp) = restrict_wbnd(jp,kp,lp) + q * u8
    rap(j, k, l)  = rap(j, k, l)  + u1
    rap(jp,k, l)  = rap(jp,k, l)  + u2
    rap(j, kp,l)  = rap(j, kp,l)  + u3
    rap(jp,kp,l)  = rap(jp,kp,l)  + u4
    rap(j, k, lp) = rap(j, k, lp) + u5
    rap(jp,k, lp) = rap(jp,k, lp) + u6
    rap(j, kp,lp) = rap(j, kp,lp) + u7
    rap(jp,kp,lp) = rap(jp,kp,lp) + u8

  ENDDO
END do
!10 continue
do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do ii = 1, bnd%cnd%ncond
    jold = bnd%cnd%jcond(ii)
    kold = bnd%cnd%kcond(ii)
    lold = bnd%cnd%lcond(ii)
    j = jnew(jold)
    jp = jnewp(jold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    k = knew(kold)
    kp = knewp(kold)
    dely  =  ddy(kold)
    odely = oddy(kold)
    l = lnew(lold)
    lp = lnewp(lold)
    delz  =  ddz(lold)
    odelz = oddz(lold)
    u1 = odelx * odely * odelz
    u2 = delx  * odely * odelz
    u3 = odelx * dely  * odelz
    u4 = delx  * dely  * odelz
    u5 = odelx * odely * delz
    u6 = delx  * odely * delz
    u7 = odelx * dely  * delz
    u8 = delx  * dely  * delz
    q = uold(jold,kold,lold)
    itot = itot + 1
    cnt(jold,kold,lold) = cnt(jold,kold,lold) + 1
    restrict_wbnd(j, k, l)  = restrict_wbnd(j, k, l)  + q * u1
    restrict_wbnd(jp,k, l)  = restrict_wbnd(jp,k, l)  + q * u2
    restrict_wbnd(j, kp,l)  = restrict_wbnd(j, kp,l)  + q * u3
    restrict_wbnd(jp,kp,l)  = restrict_wbnd(jp,kp,l)  + q * u4
    restrict_wbnd(j, k, lp) = restrict_wbnd(j, k, lp) + q * u5
    restrict_wbnd(jp,k, lp) = restrict_wbnd(jp,k, lp) + q * u6
    restrict_wbnd(j, kp,lp) = restrict_wbnd(j, kp,lp) + q * u7
    restrict_wbnd(jp,kp,lp) = restrict_wbnd(jp,kp,lp) + q * u8
    rap(j, k, l)  = rap(j, k, l)  + u1
    rap(jp,k, l)  = rap(jp,k, l)  + u2
    rap(j, kp,l)  = rap(j, kp,l)  + u3
    rap(jp,kp,l)  = rap(jp,kp,l)  + u4
    rap(j, k, lp) = rap(j, k, lp) + u5
    rap(jp,k, lp) = rap(jp,k, lp) + u6
    rap(j, kp,lp) = rap(j, kp,lp) + u7
    rap(jp,kp,lp) = rap(jp,kp,lp) + u8
  end do
END do
!10 continue

do l = 1, nznew+1
 do k = 1, nynew+1
  do j = 1, nxnew+1
    IF(rap(j,k,l)/=0._8) restrict_wbnd(j,k,l)   = restrict_wbnd(j,k,l)   / rap(j,k,l)
  end do
 end do
end do
do l = 1, nzold+1
 do k = 1, nyold+1
  do j = 1, nxold+1
    IF(cnt(j,k,l)>1) then
      do lold = 1, nzold+1
       do kold = 1, nyold+1
        do jold = 1, nxold+1
          IF(.NOT.bnd%v(jold,kold,lold)==v_vacuum) cycle
          IF(j==jold.and.k==kold.and.l==lold) WRITE(0,*) level,': #1# ',j,k,l,cnt(j,k,l)
        END do
       END do
      END do
      do ic = 1, bnd%nb_conductors
        IF(ic==1) then
          bnd%cnd => bnd%first
        else
          bnd%cnd => bnd%cnd%next
        END if
        do ii = 1, bnd%cnd%nbbnd
          jold = bnd%cnd%jj(ii)
          kold = bnd%cnd%kk(ii)
          lold = bnd%cnd%ll(ii)
          IF(.NOT.(bnd%v(jold,kold,lold)==v_bnd.and.bnd%cnd%docalc(ii))) cycle
          IF(j==jold.and.k==kold.and.l==lold) WRITE(0,*) level,': #2# ',j,k,l,cnt(j,k,l)
        END do
        do ii = 1, bnd%cnd%ncond
          jold = bnd%cnd%jcond(ii)
          kold = bnd%cnd%kcond(ii)
          lold = bnd%cnd%lcond(ii)
          IF(j==jold.and.k==kold.and.l==lold) WRITE(0,*) level,': #3# ',j,k,l,cnt(j,k,l)
        END do
      END do
    END if
  END do
 END do
END do
DEALLOCATE(ddx, ddy, ddz, oddx, oddy, oddz, jnew, knew, lnew, jnewp, knewp, lnewp, cnt)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END function restrict_wbnd

subroutine relaxbnd3dwguard(f,rhs,bnd,nx,ny,nz,dx,dy,dz,nc,voltfact)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nx, ny, nz, nc
REAL(8), INTENT(IN OUT) :: f(0:,0:,0:)!f(0:nx+2,0:ny+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(1:,1:,1:)!rhs(nx+1,ny+1,nz+1)
REAL(8), INTENT(IN) :: dx, dy, dz, voltfact
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, k, l, ii, jsw, ksw, lsw, redblack, iil, iiu, ic, nxi, nxf, nyi, nyf, nzi, nzf
REAL(8) :: dt
REAL(8) :: cf0, cfx, cfy, cfz, cfrhs

IF(l_mgridrz_debug) WRITE(0,*) 'enter relax, level = ',level

! define CFL
dt = coefrelax/(2._8/dx**2+2._8/dy**2+2._8/dz**2)
! define coefficients
cfx = dt / dx**2
cfy = dt / dy**2
cfz = dt / dz**2
cf0 = 1._8-2._8*(cfx+cfy+cfz)
cfrhs = -dt

IF(ixlbnd==dirichlet) then
  nxi=1
else
  nxi=0
END if
IF(ixrbnd==dirichlet) then
  nxf=nx-1
else
  nxf=nx
END if
IF(iylbnd==dirichlet) then
  nyi=2
else
  nyi=1
END if
IF(iyrbnd==dirichlet) then
  nyf=ny-1
else
  nyf=ny
END if
IF(bndy(level)%izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(bndy(level)%izrbnd==dirichlet) then
  nzf=nz-1
else
  nzf=nz
END if

do i = 1, nc

lsw = 1
do redblack = 1, 2
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
      k = bnd%cnd%kk(ii)
      l = bnd%cnd%ll(ii)
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,k,l)==v_bnd) &
        f(j,k,l) = bnd%cnd%cf0(ii)*f(j,k,l) &
                 + bnd%cnd%cfxp(ii)*f(j+1,k,l)+bnd%cnd%cfxm(ii)*f(j-1,k,l) &
                 + bnd%cnd%cfyp(ii)*f(j,k+1,l)+bnd%cnd%cfym(ii)*f(j,k-1,l) &
                 + bnd%cnd%cfzp(ii)*f(j,k,l+1)+bnd%cnd%cfzm(ii)*f(j,k,l-1) &
                 + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
                 +           bnd%cnd%phi0ym(ii)+bnd%cnd%phi0yp(ii) &
                 +           bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
                 + bnd%cnd%cfrhs(ii)*rhs(j,k,l)
    ENDDO
  END do
END do !redblack=1, 2
do redblack = 1, 2
  ksw = lsw
  do l = nzi, nzf+1
   jsw = ksw
   do k = nyi, nyf+1
    do j = jsw+nxi, nxf+1, 2
      IF(bnd%v(j,k,l)==v_vacuum) &
        f(j,k,l) = cf0 * f(j,k,l) &
                 + cfx*(f(j+1,k,l)+f(j-1,k,l))   &
                 + cfy*(f(j,k+1,l)+f(j,k-1,l)) &
                 + cfz*(f(j,k,l+1)+f(j,k,l-1)) &
                 + cfrhs*rhs(j,k,l)
    end do
    jsw = 3-jsw
   end do
   ksw = 3-ksw
  end do
  lsw = 3-lsw
END do !redblack=1, 2

call updateguardcells3d(f=f,level=level)

#ifdef MPIPARALLEL
  call exchange_fbndz(f,level)
#endif

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbnd3dwguard

#ifdef MPIPARALLEL
  subroutine exchange_fbndz(f,level)
    REAL(8), INTENT(IN OUT) :: f(0:,0:,0:)!f(0:nx+2,0:ny+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: level

    INTEGER(ISZ) :: nx,ny,nz
    INTEGER(ISZ) :: p_up, p_down

    p_up   = -bndy(level)%izrbnd-1
    p_down = -bndy(level)%izlbnd-1

    nx = SIZE(f,1)-3
    ny = SIZE(f,2)-3
    nz = SIZE(f,3)-3

    ! send
    IF(bndy(level)%izlbnd<0) call mpi_send_real_array(PACK(f(:,:,2),.TRUE.), p_down, 0)
    IF(bndy(level)%izrbnd<0) call mpi_send_real_array(PACK(f(:,:,nz),.TRUE.), p_up, 0)

    ! receive
    IF(bndy(level)%izrbnd<0) f(:,:,nz+2) = RESHAPE(mpi_recv_real_array(SIZE(f(:,:,nz)),p_up,0), &
                                                                      SHAPE(f(:,:,nz)))
    IF(bndy(level)%izlbnd<0) f(:,:,0)    = RESHAPE(mpi_recv_real_array(SIZE(f(:,:,0 )),p_down,0), &
                                                                      SHAPE(f(:,:,0)))

  end subroutine exchange_fbndz
  subroutine exchange_rhobndz(rho,level)
    REAL(8), INTENT(IN OUT) :: rho(1:,1:,1:)!rho(1:nx+1,1:ny+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: level

    INTEGER(ISZ) :: nx,ny,nz
    INTEGER(ISZ) :: p_up, p_down

    p_up   = -bndy(level)%izrbnd-1
    p_down = -bndy(level)%izlbnd-1

    nx = SIZE(rho,1)-1
    ny = SIZE(rho,2)-1
    nz = SIZE(rho,3)-1

    ! send
    IF(bndy(level)%izlbnd<0) call mpi_send_real_array(PACK(rho(:,:,1),.TRUE.), p_down, 1)
    IF(bndy(level)%izrbnd<0) call mpi_send_real_array(PACK(rho(:,:,nz+1),.TRUE.), p_up, 1)

    ! receive
    IF(bndy(level)%izrbnd<0) rho(:,:,nz+1) = 0.5_8*rho(:,:,nz+1) + 0.5_8*RESHAPE(mpi_recv_real_array(SIZE(rho(:,:,nz+1)),p_up,1), &
                                                                                                    SHAPE(rho(:,:,nz+1)))
    IF(bndy(level)%izlbnd<0) rho(:,:,1)    = 0.5_8*rho(:,:,1)    + 0.5_8*RESHAPE(mpi_recv_real_array(SIZE(rho(:,:,1 ))  ,p_down,1),&
                                                                                                    SHAPE(rho(:,:,1)))
  end subroutine exchange_rhobndz
  subroutine merge_work(f,level)
    REAL(8), INTENT(IN OUT) :: f(1:,1:,1:)!f(1:nx+1,1:ny+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: level

    INTEGER(ISZ) :: nz, p_up, p_down

    nz     = bndy(level-1)%nz
    p_up   = -bndy(level)%izrbnd-1
    p_down = -bndy(level)%izlbnd-1

    IF(MOD(my_index/bndy(level)%nworkpproc,2)==0) then
    ! send up
      call mpi_send_real_array(PACK(f(:,:,1:nz/2),     .TRUE.), p_up,   2)
    ! receive up
      f(:,:,nz/2+2:nz+1) = RESHAPE(mpi_recv_real_array(SIZE(f(:,:,nz/2+2:nz+1)), p_up,2), &
                                                      SHAPE(f(:,:,nz/2+2:nz+1)))
    else
    ! send down
      call mpi_send_real_array(PACK(f(:,:,nz/2+2:nz+1),.TRUE.), p_down, 2)
    ! receive down
      f(:,:,1:nz/2)      = RESHAPE(mpi_recv_real_array(SIZE(f(:,:,1:nz/2)),       p_down,2), &
                                                      SHAPE(f(:,:,1:nz/2)))
    END if

  end subroutine merge_work
#endif

function residbnd3dwguard(f,rhs,bnd,nx,ny,nz,dx,dy,dz,voltfact,l_zerolastz)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nx, ny, nz
REAL(8), INTENT(IN) :: f(0:,0:,0:)!f(0:nx+2,0:ny+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(:,:,:)!rhs(nx+1,ny+1,nz+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dx, dy, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2,SIZE(f,3)-2) :: residbnd3dwguard
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, k, l, ii, redblack, ic, nxi, nxf, nyi, nyf, nzi, nzf
REAL(8) :: cf0, cfx, cfy, cfz

IF(l_mgridrz_debug) WRITE(0,*) 'enter resid, level = ',level

cfx = -1._8 / dx**2
cfy = -1._8 / dy**2
cfz = -1._8 / dz**2
cf0 = -2._8*(cfx+cfy+cfz)

residbnd3dwguard = 0._8

!do ic = 1, bnd%nb_conductors
!  IF(ic==1) then
!    bnd%cnd => bnd%first
!  else
!    bnd%cnd => bnd%cnd%next
!  END if
!  do i = 1, bnd%cnd%ncond
!    residbnd3dwguard(bnd%cnd%jcond(i),bnd%cnd%kcond(i),bnd%cnd%lcond(i)) = 0._8
!  end do
!END do

IF(ixlbnd==dirichlet) then
  nxi=2
else
  nxi=1
END if
IF(ixrbnd==dirichlet) then
  nxf=nx-1
else
  nxf=nx
END if
IF(iylbnd==dirichlet) then
  nyi=2
else
  nyi=1
END if
IF(iyrbnd==dirichlet) then
  nyf=ny-1
else
  nyf=ny
END if
IF(bndy(level)%izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(bndy(level)%izrbnd==dirichlet) then
  nzf=nz-1
else
  nzf=nz
END if

do l = nzi, nzf+1
 do k = nyi, nyf+1
  do j = nxi, nxf+1
     IF(bnd%v(j,k,l)==v_vacuum) &
       residbnd3dwguard(j,k,l) = cf0 * f(j,k,l)    &
                               + cfx*(f(j+1,k,l)+f(j-1,k,l)) &
                               + cfy*(f(j,k+1,l)+f(j,k-1,l)) &
                               + cfz*(f(j,k,l+1)+f(j,k,l-1)) &
                               + rhs(j,k,l)
  end do
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
    k = bnd%cnd%kk(ii)
    l = bnd%cnd%ll(ii)
    IF(bnd%v(j,k,l)==v_bnd.and.bnd%cnd%docalc(ii)) &
      residbnd3dwguard(j,k,l) = bnd%cnd%rcf0(ii)*f(j,k,l) &
                              + bnd%cnd%rcfxp(ii)*f(j+1,k,l)+bnd%cnd%rcfxm(ii)*f(j-1,k,l) &
                              + bnd%cnd%rcfyp(ii)*f(j,k+1,l)+bnd%cnd%rcfym(ii)*f(j,k-1,l) &
                              + bnd%cnd%rcfzp(ii)*f(j,k,l+1)+bnd%cnd%rcfzm(ii)*f(j,k,l-1) &
                              + voltfact*(bnd%cnd%rphi0xp(ii)+bnd%cnd%rphi0xm(ii) &
                              +           bnd%cnd%rphi0yp(ii)+bnd%cnd%rphi0ym(ii) &
                              +           bnd%cnd%rphi0zp(ii)+bnd%cnd%rphi0zm(ii)) &
                              + bnd%cnd%rcfrhs(ii)*rhs(j,k,l)
  ENDDO
END do

IF(l_zerolastz) residbnd3dwguard(:,:,nz+1) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit resid, level = ',level

return
end function residbnd3dwguard

RECURSIVE subroutine mgbnd3dwguard(j, u, rhs, bnd, nx, ny, nz, dx, dy, dz, npre, npost, ncycle, sub, relax_only, npmin)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: j, nx, ny, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(:,:,:), INTENT(IN OUT) :: u
REAL(8), DIMENSION(:,:,:), INTENT(IN) :: rhs
REAL(8) :: dx, dy, dz
TYPE(bndptr) :: bnd(:)
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(:,:,:), allocatable :: res, v
INTEGER(ISZ) :: i
INTEGER :: nxnext, nynext, nznext, nzresmin, nyresmin, nzresmax, nzres
REAL(8) :: dxnext, dynext, dznext, voltf

REAL(8) :: error1, error2

level = j

IF(l_mgridrz_debug) WRITE(0,*) 'enter mg, level = ',level

IF(sub) then
  voltf = 0._8
else
  voltf = 1._8
END if

IF(j<=npmin .or. relax_only) then
!error1 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcells3d(f=u,level=j)
  call relaxbnd3dwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,ny=ny,nz=nz,dx=dx,dy=dy,dz=dz,nc=npre,voltfact=voltf)
!error2 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'relax only, level = ',level,' error1/2 = ',error1,error2
else
  nxnext = bnd(j-1)%nx
  nynext = bnd(j-1)%ny
  nznext = bnd(j-1)%nz
  dxnext = bnd(j-1)%dx
  dynext = bnd(j-1)%dy
  dznext = bnd(j-1)%dz
  ALLOCATE(res(nxnext+1,nynext+1,nznext+1),v(0:nxnext+2,0:nynext+2,0:nznext+2))
!error1 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcells3d(f=u,level=j)
  call relaxbnd3dwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,ny=ny,nz=nz,dx=dx,dy=dy,dz=dz,nc=npre,voltfact=voltf)
!error2 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'relax 1, level = ',level,' error1/2 = ',error1,error2
!error1 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  IF(bnd(level-1)%l_merged) then
#ifdef MPIPARALLEL
    IF(MOD(my_index/bnd(level)%nworkpproc,2)==0) then
      nzresmin = 1
      nzresmax = nznext/2+1
    else
      nzresmin = nznext/2+1
      nzresmax = nznext+1
    END if
    nzres = nznext/2
#endif
  else
    nzresmin = 1
    nzresmax = nznext+1
    nzres = nznext
  END if
  IF(restrictwbnd) then
  res(:,:,nzresmin:nzresmax) = restrict_wbnd( &
                             residbnd3dwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,ny=ny,nz=nz,dx=dx,dy=dy,dz=dz, &
                                              voltfact=voltf,l_zerolastz=.false.), &
                             bnd(j),nxnext,nynext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  else
  res(:,:,nzresmin:nzresmax) = restrict( &
                             residbnd3dwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,ny=ny,nz=nz,dx=dx,dy=dy,dz=dz, &
                                              voltfact=voltf,l_zerolastz=.false.), &
                             nxnext,nynext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
#ifdef MPIPARALLEL
  IF(bnd(level-1)%l_merged) then
    call merge_work(res,level)
  else
    call exchange_rhobndz(res,level-1)
  END if
#endif
  call apply_voltage(res,bnd(j-1),0._8)
  v = 0.0_8
  do i = 1, ncycle  !(1=V cycles, 2=W cycle)
    call mgbnd3dwguard(j=j-1, u=v, rhs=res, bnd=bnd(1:j-1), nx=nxnext, ny=nynext, nz=nznext, dx=dxnext, dy=dynext, dz=dznext, &
                       npre=npre, npost=npost, ncycle=ncycle, sub=.TRUE., relax_only=.FALSE., npmin=npmin)
    level = j
  end do
  call apply_voltagewguard(v,bnd(j-1),0._8)
call updateguardcells3d(f=v,level=j-1)
!  IF(bnd(j)%l_powerof2) then
!    u = u + expandwguard(v(:,nzresmin-1:nzresmax+1))
!  else
!  IF(restrictwbnd) then
!    u = u + expandwguardandbnd_any(v(:,:,nzresmin-1:nzresmax+1),bnd(j),nx,ny,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
!  else
    u = u + expandwguard_any(v(:,:,nzresmin-1:nzresmax+1),nx,ny,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
!  END if
!  END if
#ifdef MPIPARALLEL
  call exchange_fbndz(u,level)
#endif
!error2 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'V cycle, level = ',level,' error1/2 = ',error1,error2
!error1 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcells3d(f=u,level=j)
  call relaxbnd3dwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,ny=ny,nz=nz,dx=dx,dy=dy,dz=dz,nc=npost,voltfact=voltf)
!error2 = SUM(ABS(residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nx,nz=nz,dx=dx,dz=dz,voltfact=voltf,l_zerolastz=.false.)))
!IF(error2>error1) WRITE(0,*) 'relax 2, level = ',level,' error1/2 = ',error1,error2
  DEALLOCATE(res,v)
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit mg, level = ',level

return
end subroutine mgbnd3dwguard

subroutine updateguardcells3d(f,level)
! update guard cells values according to boundary conditions.
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:,:)
INTEGER(ISZ) :: level

INTEGER(ISZ) :: ixmax, iymax, izmax

ixmax=SIZE(f,1)
iymax=SIZE(f,2)
izmax=SIZE(f,3)

select case (ixlbnd)
    case (dirichlet)
      f(1,:,:) = 2._8*f(2,:,:)-f(3,:,:)
    case (neumann)
      f(1,:,:) = f(3,:,:)
    case default
end select
select case (ixrbnd)
    case (dirichlet)
      f(ixmax,:,:) = 2._8*f(ixmax-1,:,:) - f(ixmax-2,:,:)
    case (neumann)
      f(ixmax,:,:) = f(ixmax-2,:,:)
    case default
end select
select case (iylbnd)
    case (dirichlet)
      f(:,1,:) = 2._8*f(:,2,:) - f(:,3,:)
    case (neumann)
      f(:,1,:) = f(:,3,:)
    case default
end select
select case (iyrbnd)
    case (dirichlet)
      f(:,iymax,:) = 2._8*f(:,iymax-1,:) - f(:,iymax-2,:)
    case (neumann)
      f(:,iymax,:) = f(:,iymax-2,:)
    case default
end select
select case (bndy(level)%izlbnd)
    case (dirichlet)
      f(:,:,1) = 2._8*f(:,:,2) - f(:,:,3)
    case (neumann)
      f(:,:,1) = f(:,:,3)
    case (periodic)
      f(:,:,1) = f(:,:,izmax-1)
    case default
end select
select case (bndy(level)%izrbnd)
    case (dirichlet)
      f(:,:,izmax) = 2._8*f(:,:,izmax-1) - f(:,:,izmax-2)
    case (neumann)
      f(:,:,izmax) = f(:,:,izmax-2)
    case (periodic)
      f(:,:,izmax) = f(:,:,2)
    case default
end select

end subroutine updateguardcells3d

subroutine apply_voltage(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:,:)
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
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i),bnd%cnd%lcond(i)) = coef_voltage*bnd%cnd%voltage(i)
  end do
END do

return
end subroutine apply_voltage

subroutine apply_voltagewguard(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors. Grid is assumed to have guard cells.
implicit none
REAL(8),INTENT(IN OUT) :: f(0:,0:,0:)
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
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i),bnd%cnd%lcond(i)) = coef_voltage*bnd%cnd%voltage(i)
  end do
END do

return
end subroutine apply_voltagewguard

subroutine solve_full_multigrid3d_jlv(u,rhoinit,bnd,nx0,ny0,nz0,length_x,length_y,length_z, &
                                  accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nx0, ny0, nz0, &  ! number of meshes in X, Y and Z (excluding guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle            ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
INTEGER(ISZ), INTENT(IN OUT) :: nrecurs_min   ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:,:)           ! u(0:nx0+2,0:ny0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:,:)         ! rhoinit(nx0+1,ny0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_x, length_y, length_z, &  ! grid lengths in X, Y and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

INTEGER(ISZ) :: i, j, nxc, nyc, nzc, npmin
REAL(8), allocatable, DIMENSION(:,:,:) :: f, fold
REAL(8) :: dx, dy, dz
TYPE rhoptr
  REAL(8), POINTER:: a(:,:,:)
END TYPE rhoptr
TYPE(rhoptr), ALLOCATABLE :: rho(:)
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged
INTEGER(ISZ):: nzresmin,nzresmax,nzres

do_calc=.true.
has_diverged = .false.

nxc = bnd(nlevels)%nx
nyc = bnd(nlevels)%ny
nzc = bnd(nlevels)%nz

ALLOCATE(f(0:bnd(1)%nx+2,0:bnd(1)%ny+2,0:bnd(1)%nz+2), rho(nlevels))
ALLOCATE(rho(nlevels)%a(nxc+1,nyc+1,nzc+1))
rho(nlevels)%a = rhoinit
f = 0._8
!f(1:bnd(1)%nx+1,1:bnd(1)%nz+1) = restrict(uold=u, nxnew=bnd(1)%nx, nznew=bnd(1)%nz, &
!                     xminold=0._8, xmaxold=length_r, zminold=0._8, zmaxold=length_z, &
!                     xminnew=0._8, xmaxnew=length_r, zminnew=0._8, zmaxnew=length_z)

do i = nlevels-1, 1, -1
  ALLOCATE(rho(i)%a(bnd(i)%nx+1,bnd(i)%ny+1,bnd(i)%nz+1))
  rho(i)%a = 0.
  level = i
  IF(bnd(level)%l_merged) then
#ifdef MPIPARALLEL
    IF(MOD(my_index/bnd(level+1)%nworkpproc,2)==0) then
      nzresmin = 1
      nzresmax = bnd(level)%nz/2+1
    else
      nzresmin = bnd(level)%nz/2+1
      nzresmax = bnd(level)%nz+1
    END if
    nzres = bnd(level)%nz/2
#endif
  else
    nzresmin = 1
    nzresmax = bnd(level)%nz+1
    nzres = bnd(level)%nz
  END if
  level = i+1
  IF(restrictwbnd) then
  rho(i)%a(:,:,nzresmin:nzresmax) = restrict_wbnd(rho(i+1)%a,bnd(level),bnd(i)%nx,bnd(i)%ny,nzres, &
                                                  0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  else
  rho(i)%a(:,:,nzresmin:nzresmax) = restrict(rho(i+1)%a,bnd(i)%nx,bnd(i)%ny,nzres, &
                                             0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
#ifdef MPIPARALLEL
  IF(bnd(level-1)%l_merged) then
    call merge_work(rho(i)%a,level)
  else
    call exchange_rhobndz(rho(i)%a,level-1)
  END if
#endif
end do

npmin=MIN(nrecurs_min,MAX(nlevels,1))
main_loop: do i = 1, nlevels
  IF(l_mgridrz_debug) WRITE(0,*) 'begin main loop, level = ',i
  level = i
  nxc = bnd(i)%nx
  nyc = bnd(i)%ny
  nzc = bnd(i)%nz
  dx = bnd(i)%dx
  dy = bnd(i)%dy
  dz = bnd(i)%dz
  IF(i>1) then
    IF(ALLOCATED(fold)) DEALLOCATE(fold)
    ALLOCATE(fold(0:bnd(i-1)%nx+2,0:bnd(i-1)%ny+2,0:bnd(i-1)%nz+2))
    fold = f
    DEALLOCATE(f)
    ALLOCATE(f(0:nxc+2,0:nyc+2,0:nzc+2))
    IF(bnd(i-1)%l_merged) then
#ifdef MPIPARALLEL
      IF(MOD(my_index/bnd(i)%nworkpproc,2)==0) then
        nzresmin = 1
        nzresmax = bnd(i-1)%nz/2+1
      else
        nzresmin = bnd(i-1)%nz/2+1
        nzresmax = bnd(i-1)%nz+1
      END if
      nzres = bnd(i-1)%nz/2
#endif
    else
      nzresmin = 1
      nzresmax = bnd(i-1)%nz+1
      nzres = bnd(i-1)%nz
    END if
!    IF(bnd(i)%l_powerof2) then
!      f = expandwguard(fold(:,nzresmin-1:nzresmax+1))
!   else
      f = expandwguard_any(fold(:,:,nzresmin-1:nzresmax+1),nxc,nyc,nzc, &
                      xminold=0._8, xmaxold=length_x, yminold=0._8, ymaxold=length_y, zminold=0._8, zmaxold=length_z, &
                      xminnew=0._8, xmaxnew=length_x, yminnew=0._8, ymaxnew=length_y, zminnew=0._8, zmaxnew=length_z)
!    END if
!      f = expandwguardandbnd_any(fold(:,:,nzresmin-1:nzresmax+1),bnd(i),nxc,nyc,nzc, &
!                      xminold=0._8, xmaxold=length_x, yminold=0._8, ymaxold=length_y, zminold=0._8, zmaxold=length_z, &
!                      xminnew=0._8, xmaxnew=length_x, yminnew=0._8, ymaxnew=length_y, zminnew=0._8, zmaxnew=length_z)
#ifdef MPIPARALLEL
    call exchange_fbndz(f,i)
#endif
    call apply_voltagewguard(f,bnd(i),1._8)
    DEALLOCATE(fold)
  END if
  ALLOCATE(fold(0:nxc+2,0:nyc+2,0:nzc+2))
  fold = f
!f = 0. ! for debugging purpose only
!  npmin=MIN(nrecurs_min,MAX(nlevels-1,1))
  do_calc=.true.
  do while(do_calc)
    average_residue_init=SUM(ABS(residbnd3dwguard(f=RESHAPE((/(0._8*j,j=1,(nxc+3)*(nyc+3)*(nzc+3))/), (/nxc+3,nyc+3,nzc+3/)) &
                   ,rhs=rho(i)%a,bnd=bnd(i),nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef MPIPARALLEL
    average_residue_init = mpi_global_compute_real(average_residue_init,MPI_SUM)/nslaves
#endif
    if(average_residue_init==0._8) then
      average_residue=0.
      average_residue_init=1.
      cycle main_loop
    END if
    average_residue_prev=average_residue_init
    do  j = 1, nc
      call mgbnd3dwguard(j=i,u=f,rhs=rho(i)%a,bnd=bnd,nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz, &
                         npre=npre,npost=npost,ncycle=ncycle,sub=.FALSE., &
      relax_only=.false.,npmin=npmin)
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
      average_residue=SUM(ABS(residbnd3dwguard(f=f,rhs=rho(i)%a,bnd=bnd(i), &
                              nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef MPIPARALLEL
      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(average_residue/average_residue_prev>=1. .and. average_residue/average_residue_init>=1.) then
        IF(npmin<i) then
          has_diverged = .true.
          WRITE(0,*) 'WARNING multigrid3d, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          f=fold
          do_calc=.true.
          npmin=npmin+1
          WRITE(0,*) '        trying npmin = ',npmin
          exit
        else
          WRITE(0,*) 'WARNING multigrid3d, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          f=fold
          do_calc=.false.
          exit
        END if
        exit
      END if
      IF(i==nlevels) then
        IF(average_residue/average_residue_init <= accuracy) then
          do_calc=.false.
          exit
        END if
      else
        IF(average_residue/average_residue_init <= sub_accuracy) then
          do_calc=.false.
          exit
        END if
      END if
      average_residue_prev=average_residue
    end do
    IF(j>=nc) do_calc=.false.
  END do
end do main_loop

WRITE(0,'("multigrid3d: precision = ",e12.5, " after ",i5," iterations...")') average_residue/average_residue_init,j
u = f
call updateguardcells3d(f=u,level=nlevels)

do i = nlevels, 1, -1
  DEALLOCATE(rho(i)%a)
end do
DEALLOCATE(f, rho)

IF(nrecurs_min/=npmin) then
  WRITE(0,'("WARNING multigrid3d, nrecurs_min = ",i2,", npmin = ",i2,". Setting nrecurs_min = ",i2,".")') nrecurs_min, npmin, npmin
  nrecurs_min = npmin
END if

return
end subroutine solve_full_multigrid3d_jlv

subroutine solve_multigrid3d_jlv(u,rhoinit,bnd,nx0,ny0,nz0,length_x,length_y,length_z, &
                             accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nx0, ny0, nz0, &  ! number of meshes in R and Z (excluing guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle            ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
INTEGER(ISZ), INTENT(IN OUT) :: nrecurs_min   ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:,:)           ! u(0:nx0+2,0:ny0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:,:)         ! rhoinit(nx0+1,ny0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_x, length_y, length_z, &  ! grid lengths in X, Y and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

INTEGER(ISZ) :: i, j, nxc, nyc, nzc, npmin
!REAL(8), allocatable, DIMENSION(:,:) :: f!, fold
REAL(8) :: uold(SIZE(u,1),SIZE(u,2),SIZE(u,3))
REAL(8) :: dx, dy, dz
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged

do_calc=.true.
has_diverged = .false.

nxc = bnd(nlevels)%nx
nyc = bnd(nlevels)%ny
nzc = bnd(nlevels)%nz

!ALLOCATE(f(0:bnd(nlevels)%nx+2,0:bnd(nlevels)%nz+2))
!u = 0._8

npmin=MIN(nrecurs_min,MAX(nlevels,1))
  level = nlevels
  i = nlevels
  nxc = bnd(i)%nx
  nyc = bnd(i)%ny
  nzc = bnd(i)%nz
  dx = bnd(i)%dx
  dy = bnd(i)%dy
  dz = bnd(i)%dz
!  ALLOCATE(fold(0:bnd(nlevels)%nx+2,0:bnd(nlevels)%nz+2))
  do_calc=.true.
  do while(do_calc)
    average_residue_init=SUM(ABS(residbnd3dwguard(f=RESHAPE((/(0._8*j,j=1,(nxc+3)*(nyc+3)*(nzc+3))/), (/nxc+3,nyc+3,nzc+3/)) &
                   ,rhs=rhoinit,bnd=bnd(i),nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef MPIPARALLEL
    average_residue_init = mpi_global_compute_real(average_residue_init,MPI_SUM)/nslaves
#endif
    if(average_residue_init==0._8) then
      average_residue=0.
      average_residue_init=1.
!      cycle main_loop
    END if
    average_residue_prev=average_residue_init
    do  j = 1, nc
      uold = u
      call mgbnd3dwguard(j=i,u=u,rhs=rhoinit,bnd=bnd,nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz, &
                         npre=npre,npost=npost,ncycle=ncycle,sub=.FALSE., &
      relax_only=.false.,npmin=npmin)
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
      average_residue=SUM(ABS(residbnd3dwguard(f=u,rhs=rhoinit,bnd=bnd(i), &
                              nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef MPIPARALLEL
      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(average_residue/average_residue_prev>=1. .and. average_residue/average_residue_init>=1.) then
        IF(npmin<i) then
          has_diverged = .true.
          WRITE(0,*) 'WARNING multigrid3d, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.true.
          npmin=npmin+1
          WRITE(0,*) '        trying npmin = ',npmin
          exit
        else
          WRITE(0,*) 'WARNING multigrid3d, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',average_residue_init
          WRITE(0,*) '        average current residue = ',average_residue
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.false.
          exit
        END if
        exit
      END if
      IF(average_residue/average_residue_init <= accuracy) then
        do_calc=.false.
        exit
      END if
      average_residue_prev=average_residue
    end do
    IF(j>=nc) do_calc=.false.
  end do

WRITE(0,'("multigrid3d: precision = ",e12.5, " after ",i5," iterations...")') average_residue/average_residue_init,j
!u = f
call updateguardcells3d(f=u,level=nlevels)

IF(nrecurs_min/=npmin) then
  WRITE(0,'("WARNING multigrid3d, nrecurs_min = ",i2,", npmin = ",i2,". Setting nrecurs_min = ",i2,".")') nrecurs_min, npmin, npmin
  nrecurs_min = npmin
END if

return
end subroutine solve_multigrid3d_jlv

subroutine solve_multigrid3d_jlv2(u,rhoinit,bnd,nx0,ny0,nz0,length_x,length_y,length_z, &
                             accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nx0, ny0, nz0, &  ! number of meshes in R and Z (excluing guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle            ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
INTEGER(ISZ), INTENT(IN OUT) :: nrecurs_min   ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:,:)           ! u(0:nx0+2,0:ny0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:,:)         ! rhoinit(nx0+1,ny0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_x, length_y, length_z, &  ! grid lengths in X, Y and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

INTEGER(ISZ) :: i, j, nxc, nyc, nzc, npmin
!REAL(8), allocatable, DIMENSION(:,:) :: f!, fold
REAL(8) :: uold(SIZE(u,1),SIZE(u,2),SIZE(u,3)), maxerr, maxerr_old
REAL(8) :: dx, dy, dz
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged

do_calc=.true.
has_diverged = .false.

nxc = bnd(nlevels)%nx
nyc = bnd(nlevels)%ny
nzc = bnd(nlevels)%nz

!ALLOCATE(f(0:bnd(nlevels)%nx+2,0:bnd(nlevels)%nz+2))
!u = 0._8

npmin=MIN(nrecurs_min,MAX(nlevels,1))
  level = nlevels
  i = nlevels
  nxc = bnd(i)%nx
  nyc = bnd(i)%ny
  nzc = bnd(i)%nz
  dx = bnd(i)%dx
  dy = bnd(i)%dy
  dz = bnd(i)%dz
!  ALLOCATE(fold(0:bnd(nlevels)%nx+2,0:bnd(nlevels)%nz+2))
  do_calc=.true.
  do while(do_calc)
!    average_residue_init=SUM(ABS(residbnd3dwguard(f=RESHAPE((/(0._8*j,j=1,(nxc+3)*(nyc+3)*(nzc+3))/), (/nxc+3,nyc+3,nzc+3/)) &
!                   ,rhs=rhoinit,bnd=bnd(i),nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef MPIPARALLEL
!    average_residue_init = mpi_global_compute_real(average_residue_init,MPI_SUM)/nslaves
#endif
!    if(average_residue_init==0._8) then
!      average_residue=0.
!      average_residue_init=1.
!      cycle main_loop
!    END if
!    average_residue_prev=average_residue_init
    maxerr = 1.
    do  j = 1, nc
      uold = u
      uold=u
      call mgbnd3dwguard(j=i,u=u,rhs=rhoinit,bnd=bnd,nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz, &
                         npre=npre,npost=npost,ncycle=ncycle,sub=.FALSE., &
      relax_only=.false.,npmin=npmin)
      maxerr_old = maxerr
      maxerr = maxval(abs(u-uold))
!      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
!      average_residue=SUM(ABS(residbnd3dwguard(f=u,rhs=rhoinit,bnd=bnd(i), &
!                              nx=nxc,ny=nyc,nz=nzc,dx=dx,dy=dy,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nxc*nyc*nzc)
#ifdef MPIPARALLEL
!      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
!      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(maxerr/maxerr_old>=1..and.j>1) then
        IF(npmin<i) then
          has_diverged = .true.
          WRITE(0,*) 'WARNING multigrid3d, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',maxerr_old
          WRITE(0,*) '        average current residue = ',maxerr
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.true.
          npmin=npmin+1
          WRITE(0,*) '        trying npmin = ',npmin
          exit
        else
          WRITE(0,*) 'WARNING multigrid3d, calculation is diverging:'
          WRITE(0,*) '        average initial residue = ',maxerr_old
          WRITE(0,*) '        average current residue = ',maxerr    
          WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
          u=uold
          do_calc=.false.
          exit
        END if
        exit
      END if
      IF(maxerr <= accuracy) then
        do_calc=.false.
        exit
      END if
    end do
    IF(j>=nc) do_calc=.false.
  end do

WRITE(0,'("multigrid3d: precision = ",e12.5, " after ",i5," iterations...")') maxerr,j
!u = f
call updateguardcells3d(f=u,level=nlevels)

IF(nrecurs_min/=npmin) then
  WRITE(0,'("WARNING multigrid3d, nrecurs_min = ",i2,", npmin = ",i2,". Setting nrecurs_min = ",i2,".")') nrecurs_min, npmin, npmin
  nrecurs_min = npmin
END if

return
end subroutine solve_multigrid3d_jlv2

END module multigrid3d_jlv

subroutine multigrid3d_jlvf(iwhich,u0,rho0,nx0,ny0,nz0,dx0,dy0,dz0, &
                            accuracy,ncmax,npre,npost,ncycle)
USE multigrid3d_jlv
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich,nx0, ny0, nz0, ncmax,npre,npost,ncycle
REAL(8), INTENT(IN OUT) :: u0(1:nx0+1,1:ny0+1,0:nz0+2)
REAL(8), INTENT(IN) :: rho0(nx0+1,ny0+1,nz0+1)
REAL(8), INTENT(IN) :: dx0, dy0, dz0, accuracy

real(8) :: u(0:nx0+2,0:ny0+2,0:nz0+2)

INTEGER(ISZ) :: j,k,l,ic,i,ii

IF(ncmax==0) return

  IF(l_mgridrz_debug) WRITE(0,*) 'Enter multigrid3d_jlvf'

  IF(iwhich==0.or.iwhich==1) then
    if(.not.bndy_allocated) call init_bnd(nx0,ny0,nz0,dx0,dy0,dz0)
  END if

  IF(iwhich==1) return



GOTO 10
u0=0.
do ic = 1, bndy(nlevels)%nb_conductors
  IF(ic==1) then
    bndy(nlevels)%cnd => bndy(nlevels)%first
  else
    bndy(nlevels)%cnd => bndy(nlevels)%cnd%next
  END if
  do i = 1, bndy(nlevels)%cnd%ncond
    u0(bndy(nlevels)%cnd%jcond(i),bndy(nlevels)%cnd%kcond(i),bndy(nlevels)%cnd%lcond(i)) = 1._8
  end do
END do

do ic = 1, bndy(nlevels)%nb_conductors
  IF(ic==1) then
    bndy(nlevels)%cnd => bndy(nlevels)%first
  else
    bndy(nlevels)%cnd => bndy(nlevels)%cnd%next
  END if
  do ii = 1, bndy(nlevels)%cnd%nbbnd
    j = bndy(nlevels)%cnd%jj(ii)
    k = bndy(nlevels)%cnd%kk(ii)
    l = bndy(nlevels)%cnd%ll(ii)
    IF(bndy(nlevels)%v(j,k,l)==v_bnd.and.bndy(nlevels)%cnd%docalc(ii)) &
      u0(j,k,l) = 1.
  ENDDO
END do
return
10 continue

  u(1:nx0+1,1:ny0+1,:)=u0(1:nx0+1,1:ny0+1,:)
  u(0,:,:) = 0._8
  u(:,0,:)=0._8
  u(nx0+2,:,:)=0._8
  u(:,ny0+2,:)=0._8

! call solve_full_multigrid3d_jlv(u=u,rhoinit=-rho0/eps0,bnd=bndy,nx0=nx0,ny0=ny0,nz0=nz0, &
!                         length_x=nx0*dx0,length_y=ny0*dy0,length_z=nz0*dz0, &
!                         accuracy=accuracy,sub_accuracy=mgridrz_sub_accuracy, &
!                        nc=ncmax,npre=npre,npost=npost,ncycle=ncycle,nrecurs_min=mgridrz_nrecurs_min)
 call solve_multigrid3d_jlv2(u=u,rhoinit=-rho0/eps0,bnd=bndy,nx0=nx0,ny0=ny0,nz0=nz0, &
                         length_x=nx0*dx0,length_y=ny0*dy0,length_z=nz0*dz0, &
                         accuracy=accuracy,sub_accuracy=mgridrz_sub_accuracy, &
                        nc=ncmax,npre=npre,npost=npost,ncycle=ncycle,nrecurs_min=mgridrz_nrecurs_min)


  u0(1:nx0+1,1:ny0+1,:)=u(1:nx0+1,1:ny0+1,:)

  IF(l_mgridrz_debug) WRITE(0,*) 'Exit multigrid3d_jlvf'

return
end subroutine multigrid3d_jlvf

subroutine srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                      xmin,xmax,ymin,ymax,lshell,                      &
                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,             &
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
! call subroutine srfrvout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use PSOR3d
USE multigrid3d_jlv
implicit none
character(*) rofzfunc
real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
LOGICAL(ISZ):: lfill,lshell,l2symtry,l4symtry
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)

INTEGER(ISZ) :: i,nxc,nyc,nzc
REAL(8) :: dxc,dyc,dzc,zmin_in,zmax_in

IF(.not.bndy_allocated) call init_bnd(nx,ny,nz,dx,dy,dz)

do i = nlevels,1,-1
  nxc = bndy(i)%nx
  nyc = bndy(i)%ny
  nzc = bndy(i)%nz
  dxc = bndy(i)%dx
  dyc = bndy(i)%dy
  dzc = bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef MPIPARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
  zmax_in = zmin_in + bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
  zmax_in = zmmax
#endif

  call srfrvout_3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                      xmin,xmax,ymin,ymax,lshell,                      &
                      zmin_in,zmax_in,zbeam,dxc,dyc,dzc,nxc,nyc,nzc,             &
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

  call addconductors_3d(i,nxc,nyc,nzc)

end do
end subroutine srfrvout3d

subroutine srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,   &
                           lzend,xmin,xmax,ymin,ymax,lshell,          &
                           zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,       &
                           ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use PSOR3d
USE multigrid3d_jlv
implicit none
character(*) rminofz,rmaxofz
real(kind=8):: volt,zmin,zmax,xcent,ycent
LOGICAL(ISZ):: lzend,lshell,l2symtry,l4symtry
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)

INTEGER(ISZ) :: i,nxc,nyc,nzc
REAL(8) :: dxc,dyc,dzc,zmin_in,zmax_in

IF(.not.bndy_allocated) call init_bnd(nx,ny,nz,dx,dy,dz)

do i = nlevels,1,-1
  nxc = bndy(i)%nx
  nyc = bndy(i)%ny
  nzc = bndy(i)%nz
  dxc = bndy(i)%dx
  dyc = bndy(i)%dy
  dzc = bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef MPIPARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
  zmax_in = zmin_in + bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
  zmax_in = zmmax
#endif

  call srfrvinout_3d(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,  &
                     lzend,xmin,xmax,ymin,ymax,lshell,                &
                     zmin_in,zmax_in,zbeam,dxc,dyc,dzc,nxc,nyc,nzc,             &
                     ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry)

  call addconductors_3d(i,nxc,nyc,nzc)

end do
end subroutine srfrvinout3d

subroutine setcndtr3d(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz, &
                      l2symtry,l4symtry)
use PSOR3d
USE multigrid3d_jlv
integer(ISZ):: nx,ny,nz
real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
logical(ISZ):: l2symtry,l4symtry

INTEGER(ISZ) :: i,nxc,nyc,nzc
REAL(8) :: dxc,dyc,dzc,zmin_in

IF(.not.bndy_allocated) call init_bnd(nx,ny,nz,dx,dy,dz)

do i = nlevels,1,-1
  level = i
  nxc = bndy(i)%nx
  nyc = bndy(i)%ny
  nzc = bndy(i)%nz
  dxc = bndy(i)%dx
  dyc = bndy(i)%dy
  dzc = bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef MPIPARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
#endif

  call setcndtr_3d(xmmin,ymmin,zmin_in,zbeam,zgrid,nxc,nyc,nzc,dxc,dyc,dzc, &
                   l2symtry,l4symtry)

  call addconductors_3d(i,nxc,nyc,nzc)

end do

END subroutine setcndtr3d

subroutine addconductors_3d(i,nxc,nyc,nzc)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use PSOR3d
USE multigrid3d_jlv
implicit none

INTEGER(ISZ), INTENT(IN) :: nxc,nyc,nzc,i

INTEGER(ISZ) :: ii,iii,iv,iiv,nxbnd,nybnd,nzbnd,j,k,l
REAL(8) :: dt,dxm,dxp,dym,dyp,dzm,dzp,dxx,dyy,dzz

TYPE(conductor_type), POINTER :: cndpnt

  nxbnd=nxc+1
  nybnd=nyc+1
  nzbnd=nzc+1

  call init_bnd_sublevel(bndy(i),necndbdy+nocndbdy,ncond)

  bndy(i)%cnd%nbbndred = necndbdy
  bndy(i)%cnd%nbbnd = necndbdy + nocndbdy

  bndy(i)%cnd%ncond = ncond
  iii=0
  do ii=1,ncond
    iii = iii + 1
    j = ixcond(ii)+1
    k = iycond(ii)+1
    l = izcond(ii)+1
    bndy(i)%cnd%jcond(iii) = j
    bndy(i)%cnd%kcond(iii) = k
    bndy(i)%cnd%lcond(iii) = l
    IF(bndy(i)%v(j,k,l)==v_cond) then
      iii = iii-1
      bndy(i)%cnd%ncond = bndy(i)%cnd%ncond-1
      cycle
    END if
    bndy(i)%v(j,k,l) = v_cond
    bndy(i)%cnd%voltage(iii) = condvolt(ii)
  end do

  do ii = 1, necndbdy
   j = iecndx(ii)+1
   k = iecndy(ii)+1
   l = iecndz(ii)+1
   bndy(i)%cnd%jj(ii) = j
   bndy(i)%cnd%kk(ii) = k
   bndy(i)%cnd%ll(ii) = l
   IF(bndy(i)%v(j,k,l)==v_cond) cycle
   IF(j>=1 .and. j<=nxbnd .and. k>=1 .and. k<=nybnd .and. l>=1 .and. l<=nzbnd) then
     bndy(i)%cnd%docalc(ii)=.true.
   else
     cycle
   END if
   dxm = MIN(1._8,ecdelmx(ii))*bndy(i)%dx
   dxp = MIN(1._8,ecdelpx(ii))*bndy(i)%dx
   dym = MIN(1._8,ecdelmy(ii))*bndy(i)%dy
   dyp = MIN(1._8,ecdelpy(ii))*bndy(i)%dy
   dzm = MIN(1._8,ecdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ecdelpz(ii))*bndy(i)%dz
   bndy(i)%cnd%volt0xm(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0xp(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0ym(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0yp(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0zm(ii)=ecvolt(ii)
   bndy(i)%cnd%volt0zp(ii)=ecvolt(ii)
   IF(bndy(i)%v(j,k,l)/=v_bnd ) then
     bndy(i)%v(j,k,l) = v_bnd
   else
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=1,cndpnt%nbbndred
         IF(j==cndpnt%jj(iiv) .AND. k==cndpnt%kk(iiv) .AND. l==cndpnt%ll(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)<dxm) then
             dxm = cndpnt%dxm(iiv)
             bndy(i)%cnd%volt0xm(ii) = cndpnt%volt0xm(iiv)
           END if
           IF(cndpnt%dxp(iiv)<dxp) then
             dxp = cndpnt%dxp(iiv)
             bndy(i)%cnd%volt0xp(ii) = cndpnt%volt0xp(iiv)
           END if
           IF(cndpnt%dym(iiv)<dym) then
             dym = cndpnt%dym(iiv)
             bndy(i)%cnd%volt0ym(ii) = cndpnt%volt0ym(iiv)
           END if
           IF(cndpnt%dyp(iiv)<dyp) then
             dyp = cndpnt%dyp(iiv)
             bndy(i)%cnd%volt0yp(ii) = cndpnt%volt0yp(iiv)
           END if
           IF(cndpnt%dzm(iiv)<dzm) then
             dzm = cndpnt%dzm(iiv)
             bndy(i)%cnd%volt0zm(ii) = cndpnt%volt0zm(iiv)
           END if
           IF(cndpnt%dzp(iiv)<dzp) then
             dzp = cndpnt%dzp(iiv)
             bndy(i)%cnd%volt0zp(ii) = cndpnt%volt0zp(iiv)
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(ii)=dxm
   bndy(i)%cnd%dxp(ii)=dxp
   bndy(i)%cnd%dym(ii)=dym
   bndy(i)%cnd%dyp(ii)=dyp
   bndy(i)%cnd%dzm(ii)=dzm
   bndy(i)%cnd%dzp(ii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dx
       dyy=bndy(i)%dy
       dzz=bndy(i)%dz
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dyy=0.5_8*(dyp+dym)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   dt = coefrelaxbnd/((1._8/dxm+1._8/dxp)/dxx+(1._8/dym+1._8/dyp)/dyy+(1._8/dzm+1._8/dzp)/dzz)
   bndy(i)%cnd%cfxm(ii) = dt/(dxm*dxx)
   bndy(i)%cnd%cfxp(ii) = dt/(dxp*dxx)
   bndy(i)%cnd%cfym(ii) = dt/(dym*dyy)
   bndy(i)%cnd%cfyp(ii) = dt/(dyp*dyy)
   bndy(i)%cnd%cfzm(ii) = dt/(dzm*dzz)
   bndy(i)%cnd%cfzp(ii) = dt/(dzp*dzz)
   bndy(i)%cnd%cfrhs(ii) = -dt
   bndy(i)%cnd%cf0(ii)  = 1._8-bndy(i)%cnd%cfxm(ii)-bndy(i)%cnd%cfxp(ii) &
                              -bndy(i)%cnd%cfym(ii)-bndy(i)%cnd%cfyp(ii) &
                              -bndy(i)%cnd%cfzm(ii)-bndy(i)%cnd%cfzp(ii)
   bndy(i)%cnd%rcfxm(ii) = -1._8/(dxm*dxx)
   bndy(i)%cnd%rcfxp(ii) = -1._8/(dxp*dxx)
   bndy(i)%cnd%rcfym(ii) = -1._8/(dym*dyy)
   bndy(i)%cnd%rcfyp(ii) = -1._8/(dyp*dyy)
   bndy(i)%cnd%rcfzm(ii) = -1._8/(dzm*dzz)
   bndy(i)%cnd%rcfzp(ii) = -1._8/(dzp*dzz)
   bndy(i)%cnd%rcfrhs(ii) = 1._8
   bndy(i)%cnd%rcf0(ii)  = -bndy(i)%cnd%rcfxm(ii)-bndy(i)%cnd%rcfxp(ii) &
                           -bndy(i)%cnd%rcfym(ii)-bndy(i)%cnd%rcfyp(ii) &
                           -bndy(i)%cnd%rcfzm(ii)-bndy(i)%cnd%rcfzp(ii)
   IF(dxm>=bndy(i)%dx) then
     bndy(i)%cnd%phi0xm(ii)=0._8
     bndy(i)%cnd%rphi0xm(ii)=0._8
   else
     bndy(i)%cnd%phi0xm(ii)=bndy(i)%cnd%cfxm(ii)*bndy(i)%cnd%volt0xm(ii)
     bndy(i)%cnd%cfxm(ii)=0._8
     bndy(i)%cnd%rphi0xm(ii)=bndy(i)%cnd%rcfxm(ii)*bndy(i)%cnd%volt0xm(ii)
     bndy(i)%cnd%rcfxm(ii)=0._8
   END if
   IF(dxp>=bndy(i)%dx) then
     bndy(i)%cnd%phi0xp(ii)=0._8
     bndy(i)%cnd%rphi0xp(ii)=0._8
   else
     bndy(i)%cnd%phi0xp(ii)=bndy(i)%cnd%cfxp(ii)*bndy(i)%cnd%volt0xp(ii)
     bndy(i)%cnd%cfxp(ii)=0._8
     bndy(i)%cnd%rphi0xp(ii)=bndy(i)%cnd%rcfxp(ii)*bndy(i)%cnd%volt0xp(ii)
     bndy(i)%cnd%rcfxp(ii)=0._8
   END if
   IF(dym>=bndy(i)%dy) then
     bndy(i)%cnd%phi0ym(ii)=0._8
     bndy(i)%cnd%rphi0ym(ii)=0._8
   else
     bndy(i)%cnd%phi0ym(ii)=bndy(i)%cnd%cfym(ii)*bndy(i)%cnd%volt0ym(ii)
     bndy(i)%cnd%cfym(ii)=0._8
     bndy(i)%cnd%rphi0ym(ii)=bndy(i)%cnd%rcfym(ii)*bndy(i)%cnd%volt0ym(ii)
     bndy(i)%cnd%rcfym(ii)=0._8
   END if
   IF(dyp>=bndy(i)%dy) then
     bndy(i)%cnd%phi0yp(ii)=0._8
     bndy(i)%cnd%rphi0yp(ii)=0._8
   else
     bndy(i)%cnd%phi0yp(ii)=bndy(i)%cnd%cfyp(ii)*bndy(i)%cnd%volt0yp(ii)
     bndy(i)%cnd%cfyp(ii)=0._8
     bndy(i)%cnd%rphi0yp(ii)=bndy(i)%cnd%rcfyp(ii)*bndy(i)%cnd%volt0yp(ii)
     bndy(i)%cnd%rcfyp(ii)=0._8
   END if
   IF(dzm>=bndy(i)%dz) then
     bndy(i)%cnd%phi0zm(ii)=0._8
     bndy(i)%cnd%rphi0zm(ii)=0._8
   else
     bndy(i)%cnd%phi0zm(ii)=bndy(i)%cnd%cfzm(ii)*bndy(i)%cnd%volt0zm(ii)
     bndy(i)%cnd%cfzm(ii)=0._8
     bndy(i)%cnd%rphi0zm(ii)=bndy(i)%cnd%rcfzm(ii)*bndy(i)%cnd%volt0zm(ii)
     bndy(i)%cnd%rcfzm(ii)=0._8
   END if
   IF(dzp>=bndy(i)%dz) then
     bndy(i)%cnd%phi0zp(ii)=0._8
     bndy(i)%cnd%rphi0zp(ii)=0._8
   else
     bndy(i)%cnd%phi0zp(ii)=bndy(i)%cnd%cfzp(ii)*bndy(i)%cnd%volt0zp(ii)
     bndy(i)%cnd%cfzp(ii)=0._8
     bndy(i)%cnd%rphi0zp(ii)=bndy(i)%cnd%rcfzp(ii)*bndy(i)%cnd%volt0zp(ii)
     bndy(i)%cnd%rcfzp(ii)=0._8
   END if
  end do


  do ii = 1, nocndbdy
   iii=necndbdy+ii
   j = iocndx(ii)+1
   k = iocndy(ii)+1
   l = iocndz(ii)+1
   bndy(i)%cnd%jj(iii)  = j
   bndy(i)%cnd%kk(iii)  = k
   bndy(i)%cnd%ll(iii)  = l
   IF(bndy(i)%v(j,k,l)==v_cond) cycle
   IF(j>=1 .and. j<=nxbnd .and. k>=1 .and. k<=nybnd .and. l>=1 .and. l<=nzbnd) then
     bndy(i)%cnd%docalc(iii)=.true.
   else
     cycle
   endif
   dxm = MIN(1._8,ocdelmx(ii))*bndy(i)%dx
   dxp = MIN(1._8,ocdelpx(ii))*bndy(i)%dx
   dym = MIN(1._8,ocdelmy(ii))*bndy(i)%dy
   dyp = MIN(1._8,ocdelpy(ii))*bndy(i)%dy
   dzm = MIN(1._8,ocdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ocdelpz(ii))*bndy(i)%dz
   bndy(i)%cnd%volt0xm(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0xp(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0ym(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0yp(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0zm(iii)=ocvolt(ii)
   bndy(i)%cnd%volt0zp(iii)=ocvolt(ii)
   IF(bndy(i)%v(j,k,l)/=v_bnd ) then
     bndy(i)%v(j,k,l) = v_bnd
   else
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       do iiv=cndpnt%nbbndred+1,cndpnt%nbbnd
         IF(j==cndpnt%jj(iiv) .AND. k==cndpnt%kk(iiv) .AND. l==cndpnt%ll(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)<dxm) then
             dxm = cndpnt%dxm(iiv)
             bndy(i)%cnd%volt0xm(iii) = cndpnt%volt0xm(iiv)
           END if
           IF(cndpnt%dxp(iiv)<dxp) then
             dxp = cndpnt%dxp(iiv)
             bndy(i)%cnd%volt0xp(iii) = cndpnt%volt0xp(iiv)
           END if
           IF(cndpnt%dym(iiv)<dym) then
             dym = cndpnt%dym(iiv)
             bndy(i)%cnd%volt0ym(iii) = cndpnt%volt0ym(iiv)
           END if
           IF(cndpnt%dyp(iiv)<dyp) then
             dyp = cndpnt%dyp(iiv)
             bndy(i)%cnd%volt0yp(iii) = cndpnt%volt0yp(iiv)
           END if
           IF(cndpnt%dzm(iiv)<dzm) then
             dzm = cndpnt%dzm(iiv)
             bndy(i)%cnd%volt0zm(iii) = cndpnt%volt0zm(iiv)
           END if
           IF(cndpnt%dzp(iiv)<dzp) then
             dzp = cndpnt%dzp(iiv)
             bndy(i)%cnd%volt0zp(iii) = cndpnt%volt0zp(iiv)
           END if
         END if
       end do
     end do
   endif
   bndy(i)%cnd%dxm(iii)=dxm
   bndy(i)%cnd%dxp(iii)=dxp
   bndy(i)%cnd%dym(iii)=dym
   bndy(i)%cnd%dyp(iii)=dyp
   bndy(i)%cnd%dzm(iii)=dzm
   bndy(i)%cnd%dzp(iii)=dzp
   select case (bnd_method)
     case (egun)
       dxx=bndy(i)%dx
       dyy=bndy(i)%dy
       dzz=bndy(i)%dz
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dyy=0.5_8*(dyp+dym)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   dt = coefrelaxbnd/((1._8/dxm+1._8/dxp)/dxx+(1._8/dym+1._8/dyp)/dyy+(1._8/dzm+1._8/dzp)/dzz)
   bndy(i)%cnd%cfxm(iii) = dt/(dxm*dxx)
   bndy(i)%cnd%cfxp(iii) = dt/(dxp*dxx)
   bndy(i)%cnd%cfym(iii) = dt/(dym*dyy)
   bndy(i)%cnd%cfyp(iii) = dt/(dyp*dyy)
   bndy(i)%cnd%cfzm(iii) = dt/(dzm*dzz)
   bndy(i)%cnd%cfzp(iii) = dt/(dzp*dzz)
   bndy(i)%cnd%cfrhs(iii) = -dt
   bndy(i)%cnd%cf0(iii)  = 1._8-bndy(i)%cnd%cfxm(iii)-bndy(i)%cnd%cfxp(iii) &
                               -bndy(i)%cnd%cfym(iii)-bndy(i)%cnd%cfyp(iii) &
                               -bndy(i)%cnd%cfzm(iii)-bndy(i)%cnd%cfzp(iii)
   bndy(i)%cnd%rcfxm(iii) = -1._8/(dxm*dxx)
   bndy(i)%cnd%rcfxp(iii) = -1._8/(dxp*dxx)
   bndy(i)%cnd%rcfym(iii) = -1._8/(dym*dyy)
   bndy(i)%cnd%rcfyp(iii) = -1._8/(dyp*dyy)
   bndy(i)%cnd%rcfzm(iii) = -1._8/(dzm*dzz)
   bndy(i)%cnd%rcfzp(iii) = -1._8/(dzp*dzz)
   bndy(i)%cnd%rcfrhs(iii) = 1._8
   bndy(i)%cnd%rcf0(iii)  = -bndy(i)%cnd%rcfxm(iii)-bndy(i)%cnd%rcfxp(iii) &
                            -bndy(i)%cnd%rcfym(iii)-bndy(i)%cnd%rcfyp(iii) &
                            -bndy(i)%cnd%rcfzm(iii)-bndy(i)%cnd%rcfzp(iii)
     IF(dxm>=bndy(i)%dx) then
       bndy(i)%cnd%phi0xm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=0._8
     else
       bndy(i)%cnd%phi0xm(iii)=bndy(i)%cnd%cfxm(iii)*bndy(i)%cnd%volt0xm(iii)
       bndy(i)%cnd%cfxm(iii)=0._8
       bndy(i)%cnd%rphi0xm(iii)=bndy(i)%cnd%rcfxm(iii)*bndy(i)%cnd%volt0xm(iii)
       bndy(i)%cnd%rcfxm(iii)=0._8
     END if
     IF(dxp>=bndy(i)%dx) then
       bndy(i)%cnd%phi0xp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=0._8
     else
       bndy(i)%cnd%phi0xp(iii)=bndy(i)%cnd%cfxp(iii)*bndy(i)%cnd%volt0xp(iii)
       bndy(i)%cnd%cfxp(iii)=0._8
       bndy(i)%cnd%rphi0xp(iii)=bndy(i)%cnd%rcfxp(iii)*bndy(i)%cnd%volt0xp(iii)
       bndy(i)%cnd%rcfxp(iii)=0._8
     END if
     IF(dym>=bndy(i)%dy) then
       bndy(i)%cnd%phi0ym(iii)=0._8
       bndy(i)%cnd%rphi0ym(iii)=0._8
     else
       bndy(i)%cnd%phi0ym(iii)=bndy(i)%cnd%cfym(iii)*bndy(i)%cnd%volt0ym(iii)
       bndy(i)%cnd%cfym(iii)=0._8
       bndy(i)%cnd%rphi0ym(iii)=bndy(i)%cnd%rcfym(iii)*bndy(i)%cnd%volt0ym(iii)
       bndy(i)%cnd%rcfym(iii)=0._8
     END if
     IF(dyp>=bndy(i)%dy) then
       bndy(i)%cnd%phi0yp(iii)=0._8
       bndy(i)%cnd%rphi0yp(iii)=0._8
     else
       bndy(i)%cnd%phi0yp(iii)=bndy(i)%cnd%cfyp(iii)*bndy(i)%cnd%volt0yp(iii)
       bndy(i)%cnd%cfyp(iii)=0._8
       bndy(i)%cnd%rphi0yp(iii)=bndy(i)%cnd%rcfyp(iii)*bndy(i)%cnd%volt0yp(iii)
       bndy(i)%cnd%rcfyp(iii)=0._8
     END if
     IF(dzm>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=0._8
     else
       bndy(i)%cnd%phi0zm(iii)=bndy(i)%cnd%cfzm(iii)*bndy(i)%cnd%volt0zm(iii)
       bndy(i)%cnd%cfzm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=bndy(i)%cnd%rcfzm(iii)*bndy(i)%cnd%volt0zm(iii)
       bndy(i)%cnd%rcfzm(iii)=0._8
     END if
     IF(dzp>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=0._8
     else
       bndy(i)%cnd%phi0zp(iii)=bndy(i)%cnd%cfzp(iii)*bndy(i)%cnd%volt0zp(iii)
       bndy(i)%cnd%cfzp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=bndy(i)%cnd%rcfzp(iii)*bndy(i)%cnd%volt0zp(iii)
       bndy(i)%cnd%rcfzp(iii)=0._8
     END if
  end do

end subroutine addconductors_3d

subroutine save_bndstructure_3d(filename)
use multigrid3d_jlv
implicit none
CHARACTER(*) :: filename
INTEGER(ISZ) :: i,ic

  OPEN(10,FILE=filename,STATUS='unknown')
    WRITE(10,*) nlevels
    do i = 1, nlevels
      WRITE(10,*) bndy(i)%dx, bndy(i)%dy, bndy(i)%dz, bndy(i)%nb_conductors, bndy(i)%nx, bndy(i)%ny, bndy(i)%nz, bndy(i)%l_powerof2
      WRITE(10,*) bndy(i)%v
      do ic = 1, bndy(i)%nb_conductors
        IF(ic==1) then
          bndy(i)%cnd => bndy(i)%first
        else
          bndy(i)%cnd => bndy(i)%cnd%next
        END if
          WRITE(0,*) i,ic,nlevels,bndy(i)%cnd%nbbnd,bndy(i)%cnd%ncond, bndy(i)%cnd%nbbndred
          WRITE(10,*) bndy(i)%cnd%nbbnd,bndy(i)%cnd%ncond
          WRITE(10,*) bndy(i)%cnd%nbbndred
          WRITE(10,*) bndy(i)%cnd%voltage
        IF(bndy(i)%cnd%nbbnd>0) then
          WRITE(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp, &
                      bndy(i)%cnd%phi0ym,bndy(i)%cnd%phi0yp, &
                      bndy(i)%cnd%phi0zm,bndy(i)%cnd%phi0zp
          WRITE(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, &
                      bndy(i)%cnd%cfyp, bndy(i)%cnd%cfym, &
                      bndy(i)%cnd%cfzp, bndy(i)%cnd%cfzm, bndy(i)%cnd%cfrhs
          WRITE(10,*) bndy(i)%cnd%rphi0xm,bndy(i)%cnd%rphi0xp, &
                      bndy(i)%cnd%rphi0ym,bndy(i)%cnd%rphi0yp, &
                      bndy(i)%cnd%rphi0zm,bndy(i)%cnd%rphi0zp
          WRITE(10,*) bndy(i)%cnd%rcf0, bndy(i)%cnd%rcfxp, bndy(i)%cnd%rcfxm, &
                      bndy(i)%cnd%rcfyp, bndy(i)%cnd%rcfym, &
                      bndy(i)%cnd%rcfzp, bndy(i)%cnd%rcfzm, bndy(i)%cnd%rcfrhs
          WRITE(10,*) bndy(i)%cnd%dxm,bndy(i)%cnd%dxp, &
                      bndy(i)%cnd%dym,bndy(i)%cnd%dyp, &
                      bndy(i)%cnd%dzm,bndy(i)%cnd%dzp
          WRITE(10,*) bndy(i)%cnd%jj, bndy(i)%cnd%kk, bndy(i)%cnd%ll
          WRITE(10,*) bndy(i)%cnd%docalc
        END if
        IF(bndy(i)%cnd%ncond>0) then
          WRITE(10,*) bndy(i)%cnd%jcond, bndy(i)%cnd%kcond , bndy(i)%cnd%lcond
        END if
      end do
    end do
  CLOSE(10)

return
end subroutine save_bndstructure_3d

subroutine read_bndstructure_3d(filename)
use multigrid3d_jlv
implicit none
CHARACTER(*), INTENT(IN) :: filename
INTEGER(ISZ) :: nbbnd,ncond,i,ic,nbc

  OPEN(10,FILE=filename,STATUS='unknown')
    read(10,*) nlevels
    ALLOCATE(bndy(nlevels))
    bndy_allocated=.true.
    do i = 1, nlevels
      NULLIFY(bndy(i)%first)
      read(10,*) bndy(i)%dx, bndy(i)%dy, bndy(i)%dz, bndy(i)%nb_conductors, bndy(i)%nx, bndy(i)%ny, bndy(i)%nz, bndy(i)%l_powerof2
      ALLOCATE(bndy(i)%v(bndy(i)%nx+1,bndy(i)%ny+1,bndy(i)%nz+1))
      read(10,*) bndy(i)%v
      nbc = bndy(i)%nb_conductors
      do ic = 1, nbc
        read(10,*) nbbnd,ncond
        call init_bnd_sublevel(bndy(i),nbbnd,ncond)
        read(10,*) bndy(i)%cnd%nbbndred
        read(10,*) bndy(i)%cnd%voltage
        WRITE(0,*) i,ic,nlevels,nbbnd,ncond,bndy(i)%cnd%nbbndred
        IF(bndy(i)%cnd%nbbnd>0) then
          read(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp, &
                     bndy(i)%cnd%phi0ym,bndy(i)%cnd%phi0yp, &
                     bndy(i)%cnd%phi0zm,bndy(i)%cnd%phi0zp
          read(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, &
                     bndy(i)%cnd%cfyp, bndy(i)%cnd%cfym, &
                     bndy(i)%cnd%cfzp, bndy(i)%cnd%cfzm, bndy(i)%cnd%cfrhs
          read(10,*) bndy(i)%cnd%rphi0xm,bndy(i)%cnd%rphi0xp, &
                     bndy(i)%cnd%rphi0ym,bndy(i)%cnd%rphi0yp, &
                     bndy(i)%cnd%rphi0zm,bndy(i)%cnd%rphi0zp
          read(10,*) bndy(i)%cnd%rcf0, bndy(i)%cnd%rcfxp, bndy(i)%cnd%rcfxm, &
                     bndy(i)%cnd%rcfyp, bndy(i)%cnd%rcfym, &
                     bndy(i)%cnd%rcfzp, bndy(i)%cnd%rcfzm, bndy(i)%cnd%rcfrhs
          read(10,*) bndy(i)%cnd%dxm,bndy(i)%cnd%dxp, &
                     bndy(i)%cnd%dym,bndy(i)%cnd%dyp, &
                     bndy(i)%cnd%dzm,bndy(i)%cnd%dzp
          read(10,*) bndy(i)%cnd%jj, bndy(i)%cnd%kk, bndy(i)%cnd%ll
          read(10,*) bndy(i)%cnd%docalc
        END if
        IF(bndy(i)%cnd%ncond>0) then
          read(10,*) bndy(i)%cnd%jcond, bndy(i)%cnd%kcond, bndy(i)%cnd%lcond
        END if
      end do
    end do
  CLOSE(10)

return
end subroutine read_bndstructure_3d


