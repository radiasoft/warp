!     Last change:  JLV   5 Sep 2001    3:11 pm
#include "top.h"

module multigridrz
! module containing RZ multigrid solver
USE constant
USE PSOR3d, ONLY:boundxy,bound0,boundnz
USE FRZmgrid
#ifdef PARALLEL
  use Parallel
  use mpirz
#endif

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
  REAL(8) :: voltage
  ! stencil coefficients for relaxation iteration
  REAL(8), POINTER, DIMENSION(:) :: cf0, cfxp, cfxm, cfzp, cfzm, cfrhs
  REAL(8), POINTER, DIMENSION(:) :: phi0xm, phi0xp, phi0zm, phi0zp
  ! stencil coefficients for residue calculation
  REAL(8), POINTER, DIMENSION(:) :: rcf0, rcfxp, rcfxm, rcfzp, rcfzm, rcfrhs
  REAL(8), POINTER, DIMENSION(:) :: rphi0xm,rphi0xp,rphi0zm,rphi0zp
  ! distances from grid node to conductor
  REAL(8), POINTER, DIMENSION(:) :: dxm,dxp,dzm,dzp
  ! locations of nodes near conductor
  INTEGER(ISZ), POINTER, DIMENSION(:) :: jj, kk
  ! locations of nodes inside conductor
  INTEGER(ISZ), POINTER, DIMENSION(:) :: jcond, kcond
  ! logical array. Calculation will be performed only when docalc=.true.
  LOGICAL(ISZ), POINTER :: docalc(:)
  INTEGER(ISZ) :: nbbnd, &      ! number of nodes near conductor
                  nbbndred, &   ! number of "red" nodes near conductor (for red-black gauss-seidel)
                  ncond         ! number of nodes inside conductor
  ! next and previous linked-list elements
  TYPE(conductor_type), POINTER :: next,prev
end TYPE conductor_type

TYPE bndptr
! type for potential calculation
  ! mask array on the grid, = .true. on nodes in vacuum (i.e. not near or in conductors)
  LOGICAL(ISZ), POINTER :: v(:,:)
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
#ifdef PARALLEL
  INTEGER(ISZ) :: nworkpproc
#endif
END TYPE bndptr

TYPE(bndptr), pointer :: bndy(:)

LOGICAL(ISZ) :: bndy_allocated=.FALSE. ! flag set to true if bndy is allocated

REAL(8), parameter :: coefrelax=1._8, coefrelaxbnd=1._8 ! coefficients for relaxation steps
INTEGER, parameter :: nmeshmin=2

INTEGER(ISZ), parameter :: dirichlet=0, neumann=1, periodic=2, othertask=-1  ! boundary condition types
INTEGER(ISZ) :: ixlbnd=dirichlet, ixrbnd=dirichlet, izlbnd=dirichlet, izrbnd=dirichlet  ! boundary conditions for each side

INTEGER(ISZ), parameter :: egun=0, ecb=1
INTEGER(ISZ) :: bnd_method=egun

INTEGER(ISZ) :: nlevels ! number of multigrid levels
INTEGER(ISZ) :: level ! current multigrid level

#ifdef PARALLEL
  INTEGER(ISZ) :: nzfine, nworkpproc, workfact=3
#endif

contains

subroutine init_bnd(nr,nz,dr,dz)
! intializes grid quantities according to the number of multigrid levels and grid sizes nx and nz.
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr, dz

INTEGER(ISZ) :: i, nrp0, nzp0, nrc, nzc
REAL(8) :: drc, dzc

#ifdef PARALLEL
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
#ifdef PARALLEL
  nworkpproc = 1
#endif
  do WHILE(nrc>nmeshmin.or.nzc>nmeshmin)
    call evalnewgrid(nrc,nzc,drc,dzc)
    nlevels = nlevels + 1
#ifdef PARALLEL
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
#ifdef PARALLEL
  nworkpproc = 1
#endif
  do i = nlevels, 1, -1
    bndy(i)%l_merged=.false.
    bndy(i)%izlbnd=izlbnd
    bndy(i)%izrbnd=izrbnd
    IF(i/=nlevels) then
      call evalnewgrid(nrc,nzc,drc,dzc)
#ifdef PARALLEL
      IF(nworkpproc<nslaves.and.nworkpproc*nrc*nzc<=workfact*nzfine) then
        nworkpproc = nworkpproc*2
        bndy(i)%l_merged=.true.
      END if
#endif
    END if
#ifdef PARALLEL
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
    bndy(i)%v(:,:)=.true.
  end do
#ifdef PARALLEL
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


  return
end subroutine init_bnd

subroutine evalnewgrid(nr,nz,dr,dz)
! evaluate nr and nz at coarser level
INTEGER(ISZ), INTENT(IN OUT) :: nr, nz
REAL(8), INTENT(IN OUT) :: dr,dz

REAL(8) :: rap
INTEGER :: nrnew, nznew

  rap = dr/dz
  IF(rap>4._8/3._8.or.nr<=2) then
    nznew = MAX(2,nz/2)
    dz = dz * REAL(nz,8)/REAL(nznew,8)
    nz=nznew
  ELSE IF(rap<2._8/3._8.or.nz<=2) then
    nrnew = MAX(2,nr/2)
    dr = dr * REAL(nr,8)/REAL(nrnew,8)
    nr=nrnew
  ELSE
    nrnew = MAX(2,nr/2)
    dr = dr * REAL(nr,8)/REAL(nrnew,8)
    nznew = MAX(2,nz/2)
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
  ALLOCATE(bndl%cnd%cf0(nbnd),bndl%cnd%cfxp(nbnd),bndl%cnd%cfxm(nbnd),bndl%cnd%cfzp(nbnd),bndl%cnd%cfzm(nbnd),bndl%cnd%cfrhs(nbnd))
  ALLOCATE(bndl%cnd%rphi0xm(nbnd),bndl%cnd%rphi0xp(nbnd),bndl%cnd%rphi0zm(nbnd),bndl%cnd%rphi0zp(nbnd))
  ALLOCATE(bndl%cnd%rcf0(nbnd),bndl%cnd%rcfxp(nbnd),bndl%cnd%rcfxm(nbnd), &
           bndl%cnd%rcfzp(nbnd),bndl%cnd%rcfzm(nbnd),bndl%cnd%rcfrhs(nbnd))
  ALLOCATE(bndl%cnd%dxm(nbnd),bndl%cnd%dxp(nbnd),bndl%cnd%dzm(nbnd),bndl%cnd%dzp(nbnd))
  ALLOCATE(bndl%cnd%jj(nbnd),bndl%cnd%kk(nbnd),bndl%cnd%docalc(nbnd))
  bndl%cnd%docalc=.false.
END if
IF(ncond>0) then
  ALLOCATE(bndl%cnd%jcond(ncond),bndl%cnd%kcond(ncond))
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
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguard_any(0:nxnew+2,0:nznew+2)

INTEGER(ISZ) :: nxnew, nznew, nxold, nzold
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

#ifdef PARALLEL
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
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: restrict(1:nxnew+1,1:nznew+1),rap(1:nxnew+1,1:nznew+1)

INTEGER(ISZ) :: nxnew, nznew, nxold, nzold
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

#ifdef PARALLEL
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

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic
REAL(8) :: dt, dt0, r
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
        IF(bnd%cnd%docalc(ii)) &
        f(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + bnd%cnd%cfrhs(ii)*rhs(j,l)
      else
        IF(bnd%cnd%docalc(ii)) &
        f(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + bnd%cnd%cfrhs(ii)*rhs(j,l)
      END if
    ENDDO
  END do
  do l = 1, nz+1
    IF(jsw==2) then! origin
      j = 1
      IF(bnd%v(j,l)) &
      f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                 + 4._8*dt0*f(j+1,l)/dr**2   &
                                 + (dt0/dz**2)*(f(j,l+1)+f(j,l-1)) &
                                 - dt0*rhs(j,l)
    END if
    do j = jsw+1, nr+1, 2
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

#ifdef PARALLEL
  call exchange_fbndz(f,level)
#endif

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbndrzwguard

#ifdef PARALLEL
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
REAL(8), INTENT(IN OUT) :: rhs(:,:)!rhs(nx+1,nz+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dr, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2) :: residbndrzwguard
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, l, ii, jsw, ksw, redblack, ic
REAL(8) :: r
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
                        + bnd%cnd%rcfzp(ii)*f(j,l+1)+bnd%cnd%rcfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%rphi0xp(ii) &
                        + bnd%cnd%rphi0zp(ii)+bnd%cnd%rphi0zm(ii)) &
                        + bnd%cnd%rcfrhs(ii)*rhs(j,l)
    else
      IF(bnd%cnd%docalc(ii)) &
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
REAL(8), DIMENSION(:,:), INTENT(IN OUT) :: rhs
REAL(8) :: dr, dz
TYPE(bndptr) :: bnd(:)
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(:,:), allocatable :: res, v
INTEGER(ISZ) :: i
INTEGER :: nrnext, nznext, nzresmin, nzresmax, nzres
REAL(8) :: drnext, dznext, voltf

level = j

IF(l_mgridrz_debug) WRITE(0,*) 'enter mg, level = ',level

IF(sub) then
  voltf = 0._8
else
  voltf = 1._8
END if

IF(j<=npmin .or. relax_only) then
  call updateguardcellsrz(f=u,level=j)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf)
else
  nrnext = bnd(j-1)%nr
  nznext = bnd(j-1)%nz
  drnext = bnd(j-1)%dr
  dznext = bnd(j-1)%dz
  ALLOCATE(res(nrnext+1,nznext+1),v(0:nrnext+2,0:nznext+2))
  call updateguardcellsrz(f=u,level=j)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf)
  IF(bnd(level-1)%l_merged) then
#ifdef PARALLEL
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
  res(:,nzresmin:nzresmax) = restrict( &
                             residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false.), &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
#ifdef PARALLEL
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
  IF(bnd(j)%l_powerof2) then
    u = u + expandwguard(v(:,nzresmin-1:nzresmax+1))
  else
    u = u + expandwguard_any(v(:,nzresmin-1:nzresmax+1),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
#ifdef PARALLEL
  call exchange_fbndz(u,level)
#endif
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,level=j)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npost,voltfact=voltf)
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
      f(ixmax,:) = -f(ixmax-2,:)
    case (neumann)
      f(ixmax,:) = f(ixmax-2,:)
    case default
end select
select case (bndy(level)%izlbnd)
    case (dirichlet)
      f(:,1) = -f(:,3)
    case (neumann)
      f(:,1) = f(:,3)
    case (periodic)
      f(:,1) = f(:,izmax-1)
    case default
end select
select case (bndy(level)%izrbnd)
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
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = coef_voltage*bnd%cnd%voltage
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
      f(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = coef_voltage*bnd%cnd%voltage
  end do
END do

return
end subroutine apply_voltagewguard

subroutine solve_multigridrz(u,rhoinit,bnd,nr0,nz0,length_r,length_z,accuracy,sub_accuracy,nc,npre,npost,ncycle,nrecurs_min)
! solve field for u with density rhoinit.
implicit none

! input/output variables
INTEGER(ISZ), INTENT(IN) :: nr0, nz0, &       ! number of meshes in R and Z (excluing guard cells)
                            nc, &             ! maximum number of full-multigrid iterations
                            npre, npost, &    ! number of relaxations before and after multigrid level coarsening
                            ncycle, &         ! number of multigrid iterations (1: V-cycle, 2: VV cycle, etc.)
                            nrecurs_min       ! minimum level for recursion
REAL(8), INTENT(IN OUT) :: u(:,:)             ! u(0:nr0+2,0:nz0+2), potential
REAL(8), INTENT(IN) :: rhoinit(:,:)           ! rhoinit(nr0+1,nz0+1), charge density
REAL(8), INTENT(IN) :: length_r, length_z, &  ! grid lengths in R and Z
                       accuracy,           &  ! required average accuracy
                       sub_accuracy           ! required average accuracy at sublevel

! internal variables
TYPE(bndptr), DIMENSION(:) :: bnd

REAL(8), parameter :: pi = 3.14159265358979_8

INTEGER(ISZ) :: i, j, l, nrc, nzc, npmin
REAL(8), allocatable, DIMENSION(:,:) :: f, fold
REAL(8) :: dr, dz, drc, dzc, dr0, dz0
TYPE rhoptr
  REAL(8), POINTER:: a(:,:)
END TYPE rhoptr
TYPE(rhoptr), ALLOCATABLE :: rho(:)
REAL(8) :: average_residue, average_residue_init, average_residue_prev
LOGICAL :: do_calc, has_diverged
INTEGER(ISZ):: nzresmin,nzresmax,nzres,ia,ib

do_calc=.true.
has_diverged = .false.

dr0 = length_r/nr0
dz0 = length_z/nz0

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
#ifdef PARALLEL
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
  rho(i)%a(:,nzresmin:nzresmax) = restrict(rho(i+1)%a,bnd(i)%nr,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
#ifdef PARALLEL
  IF(bnd(level-1)%l_merged) then
    call merge_work(rho(i)%a,level)
  else
    call exchange_rhobndz(rho(i)%a,level-1)
  END if
#endif
end do

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
#ifdef PARALLEL
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
    IF(bnd(i)%l_powerof2) then
      f = expandwguard(fold(:,nzresmin-1:nzresmax+1))
    else
      f = expandwguard_any(fold(:,nzresmin-1:nzresmax+1),nrc,nzc, &
                      xminold=0._8, xmaxold=length_r, zminold=0._8, zmaxold=length_z, &
                      xminnew=0._8, xmaxnew=length_r, zminnew=0._8, zmaxnew=length_z)
    END if
#ifdef PARALLEL
    call exchange_fbndz(f,i)
#endif
    call apply_voltagewguard(f,bnd(i),1._8)
    DEALLOCATE(fold)
  END if
  ALLOCATE(fold(0:nrc+2,0:nzc+2))
  fold = f
!f = 0. ! for debugging purpose only
  npmin=MIN(nrecurs_min,MAX(nlevels-1,1))
  do_calc=.true.
  do while(do_calc)
    average_residue_init=SUM(ABS(residbndrzwguard(f=RESHAPE((/(0._8*j,j=1,(nrc+3)*(nzc+3))/), (/nrc+3,nzc+3/)) &
                             ,rhs=rho(i)%a,bnd=bnd(i),nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nrc*nzc)
#ifdef PARALLEL
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
      relax_only=.FALSE.,npmin=npmin)
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, level = ',level
      average_residue=SUM(ABS(residbndrzwguard(f=f,rhs=rho(i)%a,bnd=bnd(i), &
                              nr=nrc,nz=nzc,dr=dr,dz=dz,voltfact=1._8,l_zerolastz=.true.)))/(nrc*nzc)
#ifdef PARALLEL
      average_residue = mpi_global_compute_real(average_residue,MPI_SUM)/nslaves
#endif
      IF(l_mgridrz_debug) WRITE(0,*) 'evaluate residue, done.'
      IF(average_residue/average_residue_prev>=1. .and. average_residue/average_residue_init>=1.) then
        has_diverged = .true.
        WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
        WRITE(0,*) '        average initial residue = ',average_residue_init
        WRITE(0,*) '        average current residue = ',average_residue
        WRITE(0,*) '        level = ',i,' on ',nlevels,' levels; npmin = ',npmin
        f=fold
        do_calc=.true.
        IF(npmin>=i) do_calc=.false.
        npmin=npmin+1
        IF(do_calc) WRITE(0,*) '        trying npmin = ',npmin
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
    IF(.not.has_diverged) do_calc=.false.
  END do
end do main_loop

WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') average_residue/average_residue_init,j
u = f
call updateguardcellsrz(f=u,level=nlevels)

do i = nlevels, 1, -1
  DEALLOCATE(rho(i)%a)
end do
DEALLOCATE(f, rho)

return
end subroutine solve_multigridrz

END module multigridrz

subroutine multigridrzf(iwhich,u0,rho0,nr0,nz0,dr0,dz0,rxbnd,lzbnd,rzbnd,accuracy,ncmax,npre,npost,ncycle)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich,nr0, nz0, rxbnd, lzbnd, rzbnd,ncmax,npre,npost,ncycle
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN) :: rho0(nr0+1,nz0+1)
REAL(8), INTENT(IN) :: dr0, dz0, accuracy

INTEGER(ISZ) :: j, l, nrc, nzc, np, i
REAL(8) :: drc, dzc
REAL(8) :: wtime,timeinit
real(8) :: u(0:nr0+2,0:nz0+2)

  IF(iwhich==0.or.iwhich==1) then
    if(.not.bndy_allocated) call init_bnd(nr0,nz0,dr0,dz0)
  END if

  IF(iwhich==1) return

!  ixrbnd = rxbnd
!  izlbnd = lzbnd
!  izrbnd = rzbnd
  u(1:nr0+1,:)=u0(1:nr0+1,:)
  u(0,:) = 0._8
  u(nr0+2,:)=0._8

  call solve_multigridrz(u=u,rhoinit=-rho0/eps0,bnd=bndy,nr0=nr0,nz0=nz0, &
                         length_r=nr0*dr0,length_z=nz0*dz0, &
                         accuracy=accuracy,sub_accuracy=mgridrz_sub_accuracy, &
                         nc=ncmax,npre=npre,npost=npost,ncycle=ncycle,nrecurs_min=mgridrz_nrecurs_min)

  u0(1:nr0+1,:)=u(1:nr0+1,:)

return
end subroutine multigridrzf

subroutine srfrvoutrz(rofzfunc,volt,zmin,zmax,xcent,rmax,lfill, &
                       xmin,xmax,lshell, &
                       zmmin,zmmax,zbeam,dx,dz,nx,nz, &
                       ix_axis,xmesh)
! call subroutine srfrvout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
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

INTEGER(ISZ) :: i,ii,iii,iv,iiv,nxbnd,nzbnd,nrc,nzc,nrp0,nzp0, lfill_tmp,j,l
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz,voltxm,voltxp,voltzm,voltzp,drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(.not.bndy_allocated) call init_bnd(nx,nz,dx,dz)

do i = nlevels,1,-1
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz
  nxbnd=nrc+1
  nzbnd=nzc+1

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef PARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
  zmax_in = zmin_in + bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
  zmax_in = zmmax
#endif

  call srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,rmax,lfill, &
                   xmin,xmax,lshell, &
                   zmin_in,zmax_in,zbeam,bndy(i)%dr,bndy(i)%dz,nrc,nzc, &
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
   dxm = MIN(1._8,ecdelmx(ii))*bndy(i)%dr
   dxp = MIN(1._8,ecdelpx(ii))*bndy(i)%dr
   dzm = MIN(1._8,ecdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ecdelpz(ii))*bndy(i)%dz
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
           IF(cndpnt%dxm(iiv)/=bndy(i)%dr) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dr) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dzm(iiv)/=bndy(i)%dz) then
             dzm=cndpnt%dzm(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dzp(iiv)/=bndy(i)%dz) then
             dzp=cndpnt%dzp(iiv)
             voltzp=cndpnt%voltage
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
       bndy(i)%cnd%phi0zm(ii)=0._8
       bndy(i)%cnd%rphi0zm(ii)=0._8
     else
       bndy(i)%cnd%phi0zm(ii)=bndy(i)%cnd%cfzm(ii)*voltzm
       bndy(i)%cnd%cfzm(ii)=0._8
       bndy(i)%cnd%rphi0zm(ii)=bndy(i)%cnd%rcfzm(ii)*voltzm
       bndy(i)%cnd%rcfzm(ii)=0._8
     END if
     IF(ecdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0zp(ii)=0._8
       bndy(i)%cnd%rphi0zp(ii)=0._8
     else
       bndy(i)%cnd%phi0zp(ii)=bndy(i)%cnd%cfzp(ii)*voltzp
       bndy(i)%cnd%cfzp(ii)=0._8
       bndy(i)%cnd%rphi0zp(ii)=bndy(i)%cnd%rcfzp(ii)*voltzp
       bndy(i)%cnd%rcfzp(ii)=0._8
     END if
  end do
  do ii = 1, nocndbdy
   iii=necndbdy+ii
   bndy(i)%cnd%jj(iii)  = iocndx(ii)+1
   bndy(i)%cnd%kk(iii)  = iocndz(ii)+1
   IF(iocndx(ii)>=0 .and. iocndx(ii)<=nxbnd-1 .and. iocndz(ii)>=0 .and. iocndz(ii)<=nzbnd-1) bndy(i)%cnd%docalc(iii)=.true.

   bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) = .false.
   dxm = MIN(1._8,ocdelmx(ii))*bndy(i)%dr
   dxp = MIN(1._8,ocdelpx(ii))*bndy(i)%dr
   dzm = MIN(1._8,ocdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ocdelpz(ii))*bndy(i)%dz
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
           IF(cndpnt%dxm(iiv)/=bndy(i)%dr) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dr) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dzm(iiv)/=bndy(i)%dz) then
             dzm=cndpnt%dzm(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dzp(iiv)/=bndy(i)%dz) then
             dzp=cndpnt%dzp(iiv)
             voltzp=cndpnt%voltage
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
       bndy(i)%cnd%phi0zm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=0._8
     else
       bndy(i)%cnd%phi0zm(iii)=bndy(i)%cnd%cfzm(iii)*voltzm
       bndy(i)%cnd%cfzm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=bndy(i)%cnd%rcfzm(iii)*voltzm
       bndy(i)%cnd%rcfzm(iii)=0._8
     END if
     IF(ocdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0zp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=0._8
     else
       bndy(i)%cnd%phi0zp(iii)=bndy(i)%cnd%cfzp(iii)*voltzp
       bndy(i)%cnd%cfzp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=bndy(i)%cnd%rcfzp(iii)*voltzp
       bndy(i)%cnd%rcfzp(iii)=0._8
     END if
  end do
end do
end subroutine srfrvoutrz

subroutine srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,   &
                        lzend,xmin,xmax,lshell,             &
                        zmmin,zmmax,zbeam,dx,dz,nx,nz,      &
                        ix_axis,xmesh)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
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

INTEGER(ISZ) :: i,ii,iii,iv,iiv,nxbnd,nzbnd,nrc,nzc,nrp0,nzp0, lfill_tmp,j,l
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz,voltxm,voltxp,voltzm,voltzp,drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(.not.bndy_allocated) call init_bnd(nx,nz,dx,dz)

do i = nlevels,1,-1
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz
  nxbnd=nrc+1
  nzbnd=nzc+1

  necndbdy=0
  nocndbdy=0
  ncond = 0

#ifdef PARALLEL
  zmin_in = my_index / bndy(i)%nworkpproc * nzc * dzc
  zmax_in = zmin_in + bndy(i)%nworkpproc * nzc * dzc
#else
  zmin_in = zmmin
  zmax_in = zmmax
#endif

  call srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,   &
                        lzend,xmin,xmax,lshell,             &
                        zmin_in,zmax_in,zbeam,bndy(i)%dr,bndy(i)%dz,nrc,nzc, &
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
   dxm = MIN(1._8,ecdelmx(ii))*bndy(i)%dr
   dxp = MIN(1._8,ecdelpx(ii))*bndy(i)%dr
   dzm = MIN(1._8,ecdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ecdelpz(ii))*bndy(i)%dz
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
           IF(cndpnt%dxm(iiv)/=bndy(i)%dr) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dr) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dzm(iiv)/=bndy(i)%dz) then
             dzm=cndpnt%dzm(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dzp(iiv)/=bndy(i)%dz) then
             dzp=cndpnt%dzp(iiv)
             voltzp=cndpnt%voltage
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
       bndy(i)%cnd%phi0zm(ii)=0._8
       bndy(i)%cnd%rphi0zm(ii)=0._8
     else
       bndy(i)%cnd%phi0zm(ii)=bndy(i)%cnd%cfzm(ii)*voltzm
       bndy(i)%cnd%cfzm(ii)=0._8
       bndy(i)%cnd%rphi0zm(ii)=bndy(i)%cnd%rcfzm(ii)*voltzm
       bndy(i)%cnd%rcfzm(ii)=0._8
     END if
     IF(ecdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0zp(ii)=0._8
       bndy(i)%cnd%rphi0zp(ii)=0._8
     else
       bndy(i)%cnd%phi0zp(ii)=bndy(i)%cnd%cfzp(ii)*voltzp
       bndy(i)%cnd%cfzp(ii)=0._8
       bndy(i)%cnd%rphi0zp(ii)=bndy(i)%cnd%rcfzp(ii)*voltzp
       bndy(i)%cnd%rcfzp(ii)=0._8
     END if
  end do
  do ii = 1, nocndbdy
   iii=necndbdy+ii
   bndy(i)%cnd%jj(iii)  = iocndx(ii)+1
   bndy(i)%cnd%kk(iii)  = iocndz(ii)+1
   IF(iocndx(ii)>=0 .and. iocndx(ii)<=nxbnd-1 .and. iocndz(ii)>=0 .and. iocndz(ii)<=nzbnd-1) bndy(i)%cnd%docalc(iii)=.true.
   bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii)) = .false.
   dxm = MIN(1._8,ocdelmx(ii))*bndy(i)%dr
   dxp = MIN(1._8,ocdelpx(ii))*bndy(i)%dr
   dzm = MIN(1._8,ocdelmz(ii))*bndy(i)%dz
   dzp = MIN(1._8,ocdelpz(ii))*bndy(i)%dz
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
           IF(cndpnt%dxm(iiv)/=bndy(i)%dr) then
             dxm=cndpnt%dxm(iiv)
             voltxm=cndpnt%voltage
           END if
           IF(cndpnt%dxp(iiv)/=bndy(i)%dr) then
             dxp=cndpnt%dxp(iiv)
             voltxp=cndpnt%voltage
           END if
           IF(cndpnt%dzm(iiv)/=bndy(i)%dz) then
             dzm=cndpnt%dzm(iiv)
             voltzm=cndpnt%voltage
           END if
           IF(cndpnt%dzp(iiv)/=bndy(i)%dz) then
             dzp=cndpnt%dzp(iiv)
             voltzp=cndpnt%voltage
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
       bndy(i)%cnd%phi0zm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=0._8
     else
       bndy(i)%cnd%phi0zm(iii)=bndy(i)%cnd%cfzm(iii)*voltzm
       bndy(i)%cnd%cfzm(iii)=0._8
       bndy(i)%cnd%rphi0zm(iii)=bndy(i)%cnd%rcfzm(iii)*voltzm
       bndy(i)%cnd%rcfzm(iii)=0._8
     END if
     IF(ocdelpz(ii)>=1.) then
       bndy(i)%cnd%phi0zp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=0._8
     else
       bndy(i)%cnd%phi0zp(iii)=bndy(i)%cnd%cfzp(iii)*voltzp
       bndy(i)%cnd%cfzp(iii)=0._8
       bndy(i)%cnd%rphi0zp(iii)=bndy(i)%cnd%rcfzp(iii)*voltzp
       bndy(i)%cnd%rcfzp(iii)=0._8
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

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(0:nr,0:nz), INTENT(INOUT) :: rho
REAL(8), INTENT(IN) :: q, dr, dz, zmin

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

INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), dimension(0:nr,0:nz), INTENT(IN OUT) :: rho

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

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), INTENT(IN) :: dr, dz, zmin

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


