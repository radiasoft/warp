!     Last change:  JLV  26 Jun 2002   11:30 am
#include "top.h"

module multigrid_common

USE constant
USE PSOR3d, ONLY:boundxy,bound0,boundnz
USE FRZmgrid
#ifdef MPIPARALLEL
  use Parallel
  use mpirz
#endif

LOGICAL :: l_mgridrz_debug=.FALSE., l_timing_rz=.TRUE., restrictwbnd=.true.

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
  REAL(8), POINTER, DIMENSION(:) :: cf0, cfxp, cfxm, cfyp, cfym, cfzp, cfzm, dt
  REAL(8), POINTER, DIMENSION(:) :: phi0xm, phi0xp, phi0ym, phi0yp, phi0zm, phi0zp
  REAL(8), POINTER, DIMENSION(:) :: volt0xm, volt0xp, volt0ym, volt0yp, volt0zm, volt0zp
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

INTEGER(ISZ), parameter :: dirichlet=0, neumann=1, periodic=2, othertask=-1  ! boundary condition types
INTEGER(ISZ) :: ixlbnd=dirichlet, ixrbnd=dirichlet, iylbnd=dirichlet, iyrbnd=dirichlet, izlbnd=dirichlet, izrbnd=dirichlet  ! boundary conditions for each side

INTEGER(ISZ) :: v_vacuum=0, v_cond=1, v_bnd=2
INTEGER(ISZ), parameter :: egun=0, ecb=1
INTEGER(ISZ) :: bnd_method=egun

INTEGER(ISZ) :: nlevels ! number of multigrid levels
INTEGER(ISZ) :: level ! current multigrid level
INTEGER(ISZ) :: nb_iters ! actual number of iterations used for a solve

#ifdef MPIPARALLEL
  INTEGER(ISZ) :: nworkpproc, workfact=8
#else
  INTEGER(ISZ) :: my_index = 0
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
  REAL(8) :: dr, dz, zmin, zmax
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

TYPE grdptr
  INTEGER(ISZ) :: id,nlevels,nr,nz,jmin,jmax,lmin,lmax
  REAL(8) :: rmin,rmax,zmin,zmax,dr,dz,invdr,invdz,zminp
  REAL(8), POINTER, DIMENSION(:) :: invvol
  REAL(8), POINTER, DIMENSION(:,:) :: rho, phi
#ifdef MPIPARALLEL
  INTEGER(ISZ) :: nzp, rhominr, rhomaxr, rhominz, rhomaxz
  REAL(8), POINTER, DIMENSION(:,:) :: rhop, phip
#endif
  INTEGER(ISZ), POINTER, DIMENSION(:,:) :: loc_part, loc_part_field_dep
  INTEGER(ISZ) :: ixlbnd, ixrbnd, izlbnd, izrbnd  ! boundary conditions for each side
  TYPE(bndptr), pointer :: bnd(:)
  INTEGER(ISZ) :: npre, npost, ncycles, ncmax, npmin
  INTEGER(ISZ) :: guard_min_r,guard_max_r,guard_min_z,guard_max_z
  REAL(8) :: mgparam
  TYPE(grdptr), pointer :: next, prev, down, up
END TYPE grdptr
TYPE(grdptr), pointer :: grids, basegrid
TYPE grd_ptr
  TYPE(grdptr), POINTER :: grid
END TYPE grd_ptr
TYPE(grd_ptr), DIMENSION(:), ALLOCATABLE :: grids_ptr
INTEGER(ISZ) :: ngrids,grids_nids,n_avail_ids,avail_ids(100),level_del_grid

contains

subroutine init_basegrid(nr,nz,dr,dz,rmin,zmin)
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin
INTEGER(ISZ) :: i,j, nzp

  IF(associated(basegrid)) return

  ALLOCATE(grids)
  NULLIFY(grids%next,grids%prev,grids%down,grids%up,grids%bnd)
#ifdef MPIPARALLEL
  workfact = mgridrz_workfact
  nzp=nzpslave(my_index)
  ALLOCATE(grids%rho(nr+1,nz+1),grids%phi(0:nr+2,0:nz+2), &
           grids%rhop(nr+1,nzp+1),grids%phip(0:nr+2,0:nzp+2), &
           grids%loc_part(nr+1,nzp+1),grids%invvol(nr+1), &
           grids%loc_part_field_dep(nr+1,nzp+1),grids%invvol(nr+1))
#else
  ALLOCATE(grids%rho(nr+1,nz+1),grids%phi(0:nr+2,0:nz+2), &
           grids%loc_part(nr+1,nz+1),grids%invvol(nr+1), &
           grids%loc_part_field_dep(nr+1,nz+1),grids%invvol(nr+1))
#endif
  grids%id=1
  grids_nids=1
  grids%nr=nr
  grids%dr=dr
  grids%rmin=rmin
  grids%rmax=rmin+nr*dr
  grids%nz=nz
#ifdef MPIPARALLEL
  grids%nzp=nzp
!  grids%zminp=zpslmin(my_index)
  grids%zminp=zpslmin(0)+izpslave(my_index)*dz
#else
  grids%zminp=zmin
#endif
  grids%dz=dz
  grids%zmin=zmin
  grids%zmax=zmin+nz*dz
  grids%jmin=1
  grids%jmax=nr+1
  grids%lmin=1
  grids%lmax=nr+1
  grids%loc_part=1
  grids%loc_part_field_dep=1
  grids%mgparam = mgridrz_mgparam
  grids%npre = mgridrz_npre
  grids%npost = mgridrz_npost
  grids%ncycles = mgridrz_ncycles
  grids%ncmax = mgridrz_ncmax
  grids%npmin = mgridrz_levels_min
  grids%phi=0.
  grids%rho=0.
#ifdef MPIPARALLEL
  grids%phip=0.
  grids%rhop=0.
#endif
  grids%guard_min_r = 0
  grids%guard_max_r = 0
  grids%guard_min_z = 0
  grids%guard_max_z = 0
  basegrid => grids
  bndy => grids%bnd
  ngrids=1
  level_del_grid=0
  n_avail_ids=0
  avail_ids=-1
  grids%invdr = 1._8/dr
  grids%invdz = 1._8/dz
  ! computes divider by cell volumes to get density
  j = 1
  ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
  ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
  ! and Verboncoeur, J. of Comp. Phys.,
  grids%invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
  do j = 2, nr+1
    grids%invvol(j) = 1._8 / (2._8 * pi * real(j-1,8) * dr * dr * dz)
  end do

  call mk_grids_ptr()
  call init_bnd(grids%bnd,nr,nz,dr,dz,grids%zmin,grids%zmax)
  grids%nlevels=nlevels

  grids%ixrbnd = boundxy
  grids%izlbnd = bound0
  grids%izrbnd = boundnz

!  do i = 1,grids%nlevels, 1
!    grids%bnd(i)%izlbnd=grids%izlbnd
!    grids%bnd(i)%izrbnd=grids%izrbnd
!  END do

return
end subroutine init_basegrid

subroutine add_grid(grid,nr,nz,dri,dzi,rmini,zmini,guard_min_r,guard_max_r,guard_min_z,guard_max_z)
implicit none
TYPE(grdptr), pointer :: grid
INTEGER(ISZ), INTENT(IN) :: nr, nz, guard_min_r, guard_max_r, guard_min_z, guard_max_z
REAL(8), INTENT(IN) :: dri,dzi,rmini,zmini
TYPE(grdptr), pointer :: newgrid
INTEGER(ISZ) :: i,j,l,jmin,jmax,lmin,lmax
REAL(8) :: dr,dz,rmin,zmin

! adjust new grid boundaries to fall onto mother grid lines
! and recalculate mesh spacing for new grid

  jmin = 1+NINT((MAX(rmini,grid%rmin)-grid%rmin)/grid%dr)
  jmax = 1+NINT((MIN(rmini+nr*dri,grid%rmax)-grid%rmin)/grid%dr)
  lmin = 1+NINT((MAX(zmini,grid%zmin)-grid%zmin)/grid%dz)
  lmax = 1+NINT((MIN(zmini+nz*dzi,grid%zmax)-grid%zmin)/grid%dz)

  rmin = (jmin-1)*grid%dr
  zmin = (lmin-1)*grid%dz

  dr = ((jmax-1)*grid%dr-rmin)/nr
  dz = ((lmax-1)*grid%dz-zmin)/nz

! allocate new grid

  ALLOCATE(newgrid)
  NULLIFY(newgrid%next,newgrid%prev,newgrid%down,newgrid%up,newgrid%bnd)
  ALLOCATE(newgrid%rho(nr+1,nz+1), &
           newgrid%phi(0:nr+2,0:nz+2), &
           newgrid%loc_part(nr+1,nz+1), &
           newgrid%loc_part_field_dep(nr+1,nz+1), &
           newgrid%invvol(nr+1))

  IF(n_avail_ids==0) then
    newgrid%id=grids_nids+1
    grids_nids=grids_nids+1
  else
    newgrid%id=avail_ids(n_avail_ids)
    n_avail_ids=n_avail_ids-1
  END if
  newgrid%jmin=jmin
  newgrid%jmax=jmax
  newgrid%lmin=lmin
  newgrid%lmax=lmax
  newgrid%loc_part=newgrid%id
  newgrid%loc_part_field_dep=newgrid%id
  newgrid%phi=0.
  newgrid%rho=0.
  newgrid%nr=nr
  newgrid%dr=dr
  newgrid%rmin=rmin
  newgrid%rmax=rmin+nr*dr
  newgrid%nz=nz
  newgrid%dz=dz
  newgrid%zmin=zmin
  newgrid%zminp=zmin
  newgrid%zmax=zmin+nz*dz
  newgrid%mgparam = basegrid%mgparam
  newgrid%npre = basegrid%npre
  newgrid%npost = basegrid%npost
  newgrid%ncycles = basegrid%ncycles
  newgrid%ncmax = basegrid%ncmax
  newgrid%npmin = basegrid%npmin
  newgrid%guard_min_r = guard_min_r
  newgrid%guard_max_r = guard_max_r
  newgrid%guard_min_z = guard_min_z
  newgrid%guard_max_z = guard_max_z
  IF(associated(grid%down)) then
    IF(associated(grid%down%next)) then
      grid%down%next%prev => newgrid
      newgrid%next => grid%down%next
    END if
    newgrid%prev => grid%down
    grid%down%next => newgrid
  else
    grid%down=>newgrid
  END if
  newgrid%up => grid
  ngrids=ngrids+1
  mgridrz_ngrids = ngrids

! update grid pointers list array

  call mk_grids_ptr()

  newgrid%up%loc_part(jmin:jmax-1,lmin:lmax-1)=newgrid%id

  newgrid%up%loc_part_field_dep(jmin+guard_min_r:jmax-1-guard_max_r,lmin+guard_min_z:lmax-1-guard_max_z)=newgrid%id

! assign boundary types

  IF(ABS(newgrid%rmax-grid%rmax)<0.1*grid%dr) then
    newgrid%ixrbnd = grid%ixrbnd
  else
    newgrid%ixrbnd = dirichlet
  END if
  IF(ABS(newgrid%zmin-grid%zmin)<0.1*grid%dz) then
    newgrid%izlbnd = grid%izlbnd
  else
    newgrid%izlbnd = dirichlet
  END if
  IF(ABS(newgrid%zmax-grid%zmax)<0.1*grid%dz) then
    newgrid%izrbnd = grid%izrbnd
  else
    newgrid%izrbnd = dirichlet
  END if

! computes commodity quantities for charge deposition

  newgrid%invdr = 1._8/dr
  newgrid%invdz = 1._8/dz
  ! computes divider by cell volumes to get density
  j = 1
  ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
  ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
  ! and Verboncoeur, J. of Comp. Phys.,
  newgrid%invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
  do j = 2, nr+1
    newgrid%invvol(j) = 1._8 / (2._8 * pi * real(j-1,8) * dr * dr * dz)
  end do
  WRITE(0,'(" Add grid ID: ",i5)') newgrid%id
  call init_bnd(newgrid%bnd,nr,nz,dr,dz,newgrid%zmin,newgrid%zmax)
  newgrid%nlevels=nlevels
  do i = 1,newgrid%nlevels, 1
    newgrid%bnd(i)%izlbnd=newgrid%izlbnd
    newgrid%bnd(i)%izrbnd=newgrid%izrbnd
  END do

!  call print_structure(grid)

  return
end subroutine add_grid

RECURSIVE subroutine print_structure(grid)
implicit none
TYPE(grdptr), pointer :: grid

IF(level_del_grid==0) WRITE(0,*) 'ngrids = ',ngrids
level_del_grid=level_del_grid+1
WRITE(0,'("{",i5,":",i5,"}")') level_del_grid,grid%id
IF(associated(grid%next)) then
  WRITE(0,*) 'next'
  call print_structure(grid%next)
END if
IF(associated(grid%down)) then
  WRITE(0,*) 'down'
  call print_structure(grid%down)
END if
level_del_grid=level_del_grid-1

end subroutine print_structure

recursive subroutine del_grid(grid)
implicit none
TYPE(grdptr), pointer :: grid
INTEGER(ISZ) :: j,l,jmin,jmax,lmin,lmax

  grid%up%loc_part(grid%jmin:grid%jmax-1,grid%lmin:grid%lmax-1)=grid%up%id
  grid%up%loc_part_field_dep(grid%jmin+grid%guard_min_r:grid%jmax-1-grid%guard_max_r, &
                             grid%lmin+grid%guard_min_z:grid%lmax-1-grid%guard_max_z)=grid%up%id

  level_del_grid=level_del_grid+1
  IF(associated(grid%next)) then
    IF(associated(grid%prev)) then
      grid%next%prev => grid%prev
      grid%prev%next => grid%next
    else
      NULLIFY(grid%next%prev)
    END if
  END if
  IF(associated(grid%down)) then
    call del_grid(grid%down)
  else
    n_avail_ids=n_avail_ids+1
    avail_ids(n_avail_ids+1)=grid%id
    DEALLOCATE(grid%rho,grid%phi,grid%loc_part,grid%loc_part_field_dep,grid%invvol)
    deallocate(grid%next,grid%prev,grid%down,grid%up)
    DEALLOCATE(grid)
  END if
  ngrids=ngrids-1
  mgridrz_ngrids = ngrids
  IF(level_del_grid==1) call mk_grids_ptr()

  level_del_grid=level_del_grid-1

  return
end subroutine del_grid

subroutine mk_grids_ptr()
implicit none
INTEGER :: i

  IF(ALLOCATED(grids_ptr)) DEALLOCATE(grids_ptr)
  ALLOCATE(grids_ptr(ngrids))
  do i = 1, ngrids
    NULLIFY(grids_ptr(i)%grid)
  end do
  call assign_grids_ptr(basegrid)

  return
end subroutine mk_grids_ptr

RECURSIVE subroutine assign_grids_ptr(grid)
implicit none
TYPE(grdptr), pointer :: grid

  grids_ptr(grid%id)%grid => grid
  IF(associated(grid%down)) call assign_grids_ptr(grid%down)
  IF(associated(grid%next)) call assign_grids_ptr(grid%next)

return
end subroutine assign_grids_ptr

subroutine init_bnd(bndy,nr,nz,dr,dz,zmin,zmax)
! intializes grid quantities according to the number of multigrid levels and grid sizes nx and nz.
USE InGen3d, ONLY:l2symtry, l4symtry
USE Picglb3d, ONLY:dy
USE InMesh3d, ONLY:ny
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr, dz, zmin, zmax
TYPE(bndptr), pointer :: bndy(:)

INTEGER(ISZ) :: i, nrp0, nzp0, nrc, nzc, nrc_old, nzc_old
REAL(8) :: drc, dzc

!  dy = dr
  ny = 0

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
    nrc_old=nrc
    nzc_old=nzc
#ifdef MPIPARALLEL
    nzc = nzc * nslaves / nworkpproc
    call evalnewgrid(nrc,nzc,drc,dzc)
    nzc = nzc * nworkpproc / nslaves 
#else
    call evalnewgrid(nrc,nzc,drc,dzc)
#endif
    IF(nrc==nrc_old .AND. nzc==nzc_old) exit
    nlevels = nlevels + 1
#ifdef MPIPARALLEL
    IF(nslaves>1.and.nworkpproc<nslaves.and.nrc*nzc<=workfact*nr) then
      nworkpproc = nworkpproc*2
      nzc = nzc*2
      nzc_old = nworkpproc*2
    END if
! make sure that nz is even for parallel red/black Gauss-Seidel
    IF(nzc/2/=NINT(0.5*REAL(nzc))) then
      dzc = dzc * REAL(nzc,8)/REAL(nzc+1,8)
      nzc=nzc+1
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
!      call evalnewgrid(nrc,nzc,drc,dzc)
#ifdef MPIPARALLEL
    nzc = nzc * nslaves / nworkpproc
    call evalnewgrid(nrc,nzc,drc,dzc)
    nzc = nzc * nworkpproc / nslaves 
#else
    call evalnewgrid(nrc,nzc,drc,dzc)
#endif
#ifdef MPIPARALLEL
      IF(nslaves>1.and.nworkpproc<nslaves.and.nrc*nzc<=workfact*nr) then
        nworkpproc = nworkpproc*2
        bndy(i)%l_merged=.true.
        nzc = nzc*2
      END if
! make sure that nz is even for parallel red/black Gauss-Seidel
      IF(nzc/2/=NINT(0.5*REAL(nzc))) then
        dzc = dzc * REAL(nzc,8)/REAL(nzc+1,8)
        nzc=nzc+1
      END if
#endif
    END if
#ifdef MPIPARALLEL
    bndy(i)%zmin = zmslmin(int(my_index/nworkpproc)*nworkpproc)
    bndy(i)%zmax = zmslmax((1+int(my_index/nworkpproc))*nworkpproc-1)
#else
    bndy(i)%zmin = zmin
    bndy(i)%zmax = zmax
#endif

#ifdef MPIPARALLEL
    bndy(i)%nworkpproc = nworkpproc
    IF(my_index-nworkpproc>=0)      bndy(i)%izlbnd = -(my_index-nworkpproc+1)
    IF(my_index+nworkpproc<nslaves) bndy(i)%izrbnd = -(my_index+nworkpproc+1)
#endif
    bndy(i)%nr = nrc
    bndy(i)%nz = nzc
    bndy(i)%dr = drc
    bndy(i)%dz = dzc
    NULLIFY(bndy(i)%first)
    bndy(i)%nb_conductors = 0
    ALLOCATE(bndy(i)%v(bndy(i)%nr+1,bndy(i)%nz+1))
    bndy(i)%v(:,:)=v_vacuum
    if(my_index==0) write(0,*) i,nrc,nzc,drc,dzc,bndy(i)%l_merged
  end do
#ifdef MPIPARALLEL
  IF(nslaves>1) then
    do i = 1, nlevels
      bndy(i)%l_powerof2 = .false.
    END do
    do i = nlevels,1,-1
    END do
  END if
#else
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
#endif

  IF(mgridrz_deform) then
    mgridrz_nz = nz
    call gchange("FRZmgrid",0)
    mgridrz_xfact = 1._8
    mgridrz_yfact = 1._8
  END if

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
  ALLOCATE(bndl%cnd%phi0xm(nbnd),bndl%cnd%phi0xp(nbnd),bndl%cnd%phi0zm(nbnd),bndl%cnd%phi0zp(nbnd),bndl%cnd%dt(nbnd))
  ALLOCATE(bndl%cnd%volt0xm(nbnd),bndl%cnd%volt0xp(nbnd),bndl%cnd%volt0zm(nbnd),bndl%cnd%volt0zp(nbnd))
  ALLOCATE(bndl%cnd%cf0(nbnd),bndl%cnd%cfxp(nbnd),bndl%cnd%cfxm(nbnd),bndl%cnd%cfzp(nbnd),bndl%cnd%cfzm(nbnd))
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

function expandwguard_any(uold, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, xminold,  xmaxold,  zminold,  zmaxold, &
                          ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguard_any(0:nxnew+2,0:nznew+2)

INTEGER(ISZ) :: nxold, nzold, nrf, nzi, nzf
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

IF(ixrbnd==dirichlet) then
  nrf=nxnew
else
  nrf=nxnew+1
END if
!IF(bndy(level)%izlbnd==dirichlet) then
IF(izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
!IF(bndy(level)%izrbnd==dirichlet) then
IF(izrbnd==dirichlet) then
  nzf=nznew
else
  nzf=nznew+1
END if

do knew = nzi, nzf
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = 1+MIN(nzold,INT(zz))
  koldp(knew) = kold(knew) + 1
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
do jnew = 1, nrf
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
do knew = nzi, nzf
  k = kold(knew)
  kp = koldp(knew)
  delz = ddz(knew)
  odelz = oddz(knew)
  do jnew = 1, nrf
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

subroutine interpolate_any(unew, uold, nxnew, nznew, nxold, nzold, &
                           xminnew, xmaxnew, zminnew, zmaxnew, &
                           xminold,  xmaxold,  zminold,  zmaxold, &
                           ixrbnd, izlbnd, izrbnd, &
                           bnd_only, quad)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), DIMENSION(0:,0:), INTENT(IN OUT) :: unew
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
LOGICAL(ISZ), INTENT(IN) :: bnd_only, quad

INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, jm, km, jpp, kpp
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz, &
           s1x, s2x, s3x, s4x, s1z, s2z, s3z, s4z
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1)

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

do knew = 1, nznew+1
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = MAX(1,MIN(nzold,1+INT(zz)))
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
do jnew = 1, nxnew+1
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = MAX(1,MIN(nxold,1+INT(xx)))
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
IF(.not.quad) then
  IF(bnd_only) then
    do knew = 1, nznew+1, nznew
     IF((knew==1.and.izlbnd==dirichlet).OR.(knew==nznew+1.and.izrbnd==dirichlet)) then
     k = kold(knew)
     kp = k+1
     delz = ddz(knew)
     odelz = oddz(knew)
     do jnew = 1, nxnew+1
       j = jold(jnew)
       jp = j+1
       delx = ddx(jnew)
       odelx = oddx(jnew)
       unew(jnew,knew) = uold(j, k)  * odelx * odelz &
                       + uold(jp,k)  * delx  * odelz &
                       + uold(j, kp) * odelx * delz &
                       + uold(jp,kp) * delx  * delz
     end do
     END if
    end do
    IF(ixrbnd==dirichlet) then
     jnew = nxnew+1
     j = jold(jnew)
     jp = j+1
     delx = ddx(jnew)
     odelx = oddx(jnew)
     do knew = 2, nznew
       k = kold(knew)
       kp = k+1
       delz = ddz(knew)
       odelz = oddz(knew)
       unew(jnew,knew) = uold(j, k)  * odelx * odelz &
                       + uold(jp,k)  * delx  * odelz &
                       + uold(j, kp) * odelx * delz &
                       + uold(jp,kp) * delx  * delz
     END do
    END if
  else
    do knew = 1, nznew+1
     WRITE(0,*) 'knew',knew
     k = kold(knew)
     kp = k+1
     delz = ddz(knew)
     odelz = oddz(knew)
     do jnew = 1, nxnew+1
       j = jold(jnew)
       jp = j+1
       delx = ddx(jnew)
       odelx = oddx(jnew)
       unew(jnew,knew) = uold(j, k)  * odelx * odelz &
                       + uold(jp,k)  * delx  * odelz &
                       + uold(j, kp) * odelx * delz &
                       + uold(jp,kp) * delx  * delz
     end do
    END do
  END if
else
  IF(bnd_only) then
    do knew = 1, nznew+1, nznew
     IF((knew==1.and.izlbnd==dirichlet).OR.(knew==nznew+1.and.izrbnd==dirichlet)) then
      k   = kold(knew)
      kp  = k+1
      km  = k-1
      kpp = k+2
      delz  = ddz(knew)
      odelz = oddz(knew)
      s1z=0.166667*odelz**3
      s2z=0.666667-delz**2+0.5*delz**3
      s3z=0.666667-odelz**2+0.5*odelz**3
      s4z=0.166667*delz**3
      do jnew = 1, nxnew+1
        j   = jold(jnew)
        jp  = j+1
        jm  = j-1
        jpp = j+2
        delx = ddx(jnew)
        odelx = oddx(jnew)
        s1x=0.166667*odelx**3
        s2x=0.666667-delx**2+0.5*delx**3
        s3x=0.666667-odelx**2+0.5*odelx**3
        s4x=0.166667*delx**3
        unew(jnew,knew) = uold(jm, km ) * s1x*s1z &
                        + uold(j  ,km ) * s2x*s1z &
                        + uold(jp, km ) * s3x*s1z &
                        + uold(jpp,km ) * s4x*s1z &
                        + uold(jm, k  ) * s1x*s2z &
                        + uold(j  ,k  ) * s2x*s2z &
                        + uold(jp, k  ) * s3x*s2z &
                        + uold(jpp,k  ) * s4x*s2z &
                        + uold(jm, kp ) * s1x*s3z &
                        + uold(j  ,kp ) * s2x*s3z &
                        + uold(jp, kp ) * s3x*s3z &
                        + uold(jpp,kp ) * s4x*s3z &
                        + uold(jm, kpp) * s1x*s4z &
                        + uold(j  ,kpp) * s2x*s4z &
                        + uold(jp, kpp) * s3x*s4z &
                        + uold(jpp,kpp) * s4x*s4z
      end do
     END if
    end do
    IF(ixrbnd==dirichlet) then
      jnew = nxnew+1
      j   = jold(jnew)
      jp  = j+1
      jm  = j-1
      jpp = j+2
      delx = ddx(jnew)
      odelx = oddx(jnew)
      s1x=0.166667*odelx**3
      s2x=0.666667-delx**2+0.5*delx**3
      s3x=0.666667-odelx**2+0.5*odelx**3
      s4x=0.166667*delx**3
      do knew = 2, nznew
        k   = kold(knew)
        kp  = k+1
        km  = k-1
        kpp = k+2
        delz  = ddz(knew)
        odelz = oddz(knew)
        s1z=0.166667*odelz**3
        s2z=0.666667-delz**2+0.5*delz**3
        s3z=0.666667-odelz**2+0.5*odelz**3
        s4z=0.166667*delz**3
        unew(jnew,knew) = uold(jm, km ) * s1x*s1z &
                        + uold(j  ,km ) * s2x*s1z &
                        + uold(jp, km ) * s3x*s1z &
                        + uold(jpp,km ) * s4x*s1z &
                        + uold(jm, k  ) * s1x*s2z &
                        + uold(j  ,k  ) * s2x*s2z &
                        + uold(jp, k  ) * s3x*s2z &
                        + uold(jpp,k  ) * s4x*s2z &
                        + uold(jm, kp ) * s1x*s3z &
                        + uold(j  ,kp ) * s2x*s3z &
                        + uold(jp, kp ) * s3x*s3z &
                        + uold(jpp,kp ) * s4x*s3z &
                        + uold(jm, kpp) * s1x*s4z &
                        + uold(j  ,kpp) * s2x*s4z &
                        + uold(jp, kpp) * s3x*s4z &
                        + uold(jpp,kpp) * s4x*s4z
      END do
    END if
  else
    do knew = 1, nznew+1
      k   = kold(knew)
      kp  = k+1
      km  = k-1
      kpp = k+2
      delz  = ddz(knew)
      odelz = oddz(knew)
      s1z=0.166667*odelz**3
      s2z=0.666667-delz**2+0.5*delz**3
      s3z=0.666667-odelz**2+0.5*odelz**3
      s4z=0.166667*delz**3
      do jnew = 1, nxnew+1
        j   = jold(jnew)
        jp  = j+1
        jm  = j-1
        jpp = j+2
        delx = ddx(jnew)
        odelx = oddx(jnew)
        s1x=0.166667*odelx**3
        s2x=0.666667-delx**2+0.5*delx**3
        s3x=0.666667-odelx**2+0.5*odelx**3
        s4x=0.166667*delx**3
        unew(jnew,knew) = uold(jm, km ) * s1x*s1z &
                        + uold(j  ,km ) * s2x*s1z &
                        + uold(jp, km ) * s3x*s1z &
                        + uold(jpp,km ) * s4x*s1z &
                        + uold(jm, k  ) * s1x*s2z &
                        + uold(j  ,k  ) * s2x*s2z &
                        + uold(jp, k  ) * s3x*s2z &
                        + uold(jpp,k  ) * s4x*s2z &
                        + uold(jm, kp ) * s1x*s3z &
                        + uold(j  ,kp ) * s2x*s3z &
                        + uold(jp, kp ) * s3x*s3z &
                        + uold(jpp,kp ) * s4x*s3z &
                        + uold(jm, kpp) * s1x*s4z &
                        + uold(j  ,kpp) * s2x*s4z &
                        + uold(jp, kpp) * s3x*s4z &
                        + uold(jpp,kpp) * s4x*s4z
      end do
    END do
  END if
END if

return
END subroutine interpolate_any

subroutine combinewguard_any(unew, uold, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, xminold,  xmaxold,  zminold,  zmaxold)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), DIMENSION(0:,0:), INTENT(IN OUT) :: unew
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, jmin,jmax,kmin,kmax
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, oddx, ddz, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jold, kold, joldp, koldp

IF(l_mgridrz_debug) WRITE(0,*) 'enter combine, level = ',level

nxold = SIZE(uold,1) - 1 - 2
nzold = SIZE(uold,2) - 1 - 2

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

jmin=1+NINT((xminold-xminnew)/dxnew)
jmax=1+NINT((xmaxold-xminnew)/dxnew)
kmin=1+NINT((zminold-zminnew)/dznew)
kmax=1+NINT((zmaxold-zminnew)/dznew)
ALLOCATE(ddx(jmin:jmax), oddx(jmin:jmax), ddz(kmin:kmax), oddz(kmin:kmax), &
         jold(jmin:jmax), kold(kmin:kmax), joldp(jmin:jmax), koldp(kmin:kmax))
do knew = kmin, kmax
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = 1+MIN(nzold,INT(zz))
  koldp(knew) = kold(knew) + 1
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
do jnew = jmin, jmax
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
do knew = kmin,kmax
  k = kold(knew)
  kp = koldp(knew)
  delz = ddz(knew)
  odelz = oddz(knew)
  do jnew = jmin,jmax
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
    unew(jnew,knew) = uold(j, k)  * odelx * odelz &
                    + uold(jp,k)  * delx  * odelz &
                    + uold(j, kp) * odelx * delz &
                    + uold(jp,kp) * delx  * delz
  end do
END do
DEALLOCATE(ddx, oddx, ddz, oddz, jold, kold, joldp, koldp)

IF(l_mgridrz_debug) WRITE(0,*) 'exit combine, level = ',level

return
END subroutine combinewguard_any

function expandwguardandbnd_any(uold, bnd, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, &
                                                         xminold, xmaxold, zminold, zmaxold, &
                                                         ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguardandbnd_any(0:nxnew+2,0:nznew+2)
TYPE(bndptr) :: bnd

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, nrf, nzi, nzf
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

IF(ixrbnd==dirichlet) then
  nrf=nxnew
else
  nrf=nxnew+1
END if
!IF(bndy(level)%izlbnd==dirichlet) then
IF(izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
!IF(bndy(level)%izrbnd==dirichlet) then
IF(izrbnd==dirichlet) then
  nzf=nznew
else
  nzf=nznew+1
END if

do knew = nzi, nzf
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = 1+MIN(nzold,INT(zz))
  koldp(knew) = kold(knew) + 1
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
do jnew = 1, nrf
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
do knew = nzi, nzf
  k = kold(knew)
  kp = koldp(knew)
  delz = ddz(knew)
  odelz = oddz(knew)
  do jnew = 1, nrf
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
    IF(.NOT.bnd%v(jnew,knew)==v_vacuum) cycle
    expandwguardandbnd_any(jnew,knew) = uold(j, k)  * odelx * odelz &
                                      + uold(jp,k)  * delx  * odelz &
                                      + uold(j, kp) * odelx * delz &
                                      + uold(jp,kp) * delx  * delz
  end do
END do

!expandwguardandbnd_any(0,:) = 0._8
!expandwguardandbnd_any(nxnew+2,:) = 0._8
!expandwguardandbnd_any(:,0) = 0._8
!expandwguardandbnd_any(:,nznew+2) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END function expandwguardandbnd_any

function restrict_pof2(f, ixrbnd, izlbnd, izrbnd)
! restrict field from one grid to a coarser one. Each dimension is assumed to be a power of two.
implicit none

REAL(8), DIMENSION(1:,1:), INTENT(IN) :: f
INTEGER(ISZ), INTENT(IN) :: ixrbnd, izlbnd, izrbnd
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

function restrict(uold, nxnew, nznew, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew, &
                  ixrbnd, izlbnd, izrbnd)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: restrict(1:nxnew+1,1:nznew+1),rap(1:nxnew+1,1:nznew+1)

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jold, kold, j, k, jp, kp
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold, ixrbnd, izlbnd, izrbnd)
!  return
!END if

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
    IF(izlbnd<0) then
      ddz(1) = 0.5_8*ddz(1)
      oddz(1) = 0.5_8*oddz(1)
    END if
    IF(izrbnd<0) then
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

#ifdef MPIPARALLEL
    IF(izlbnd<0) then
      rap(:,1) = 2.*rap(:,1)
    END if
    IF(izrbnd<0) then
      rap(:,nznew+1) = 2.*rap(:,nznew+1)
    END if
#endif
do k = 1, nznew+1
  do j = 1, nxnew+1
    IF(rap(j,k)/=0._8) restrict(j,k)   = restrict(j,k)   / rap(j,k)
  end do
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END function restrict

subroutine deposit(unew, uold, invvolnew, invvolold, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! deposit rho from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
REAL(8), DIMENSION(1:,1:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), DIMENSION(1:), INTENT(IN) :: invvolnew, invvolold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew

INTEGER(ISZ) :: nxold, nzold, nxnew, nznew
INTEGER(ISZ) :: jold, kold, j, k, jp, kp
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz, volold
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

IF(l_mgridrz_debug) WRITE(0,*) 'enter deposit'

nxold = SIZE(uold,1) - 1
nzold = SIZE(uold,2) - 1
nxnew = SIZE(unew,1) - 1
nznew = SIZE(unew,2) - 1

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1), volold(nxold+1))

invdxnew = nxnew / (xmaxnew-xminnew)
invdznew = nznew / (zmaxnew-zminnew)
dxold = (xmaxold-xminold) / nxold
dzold = (zmaxold-zminold) / nzold

do kold = 1, nzold+1
  z = zminold + (kold-1)*dzold
  knew(kold) = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  ddz(kold) = (z-zminnew) * invdznew-real(knew(kold)-1)
  knewp(kold) = knew(kold)+1
  oddz(kold) = 1._8-ddz(kold)
END do

do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = 1._8-ddx(jold)
  volold(jold) = 1._8/invvolold(jold)
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
    unew(j,k)   = unew(j,k)   + uold(jold,kold) * odelx * odelz * invvolnew(j) *volold(jold)
    unew(jp,k)  = unew(jp,k)  + uold(jold,kold) * delx  * odelz * invvolnew(jp)*volold(jold)
    unew(j,kp)  = unew(j,kp)  + uold(jold,kold) * odelx * delz  * invvolnew(j) *volold(jold)
    unew(jp,kp) = unew(jp,kp) + uold(jold,kold) * delx  * delz  * invvolnew(jp)*volold(jold)
  end do
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp, volold)

IF(l_mgridrz_debug) WRITE(0,*) 'exit deposit'

return
END subroutine deposit

function restrict_wbnd(uold, bnd, nxnew, nznew, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew, &
                       ixrbnd, izlbnd, izrbnd)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixrbnd, izlbnd, izrbnd
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
    IF(izlbnd<0) then
      ddz(1) = 0.5_8*ddz(1)
      oddz(1) = 0.5_8*oddz(1)
    END if
    IF(izrbnd<0) then
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

GO TO 10
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
        u1=0._8
!        u4=u4+u3
        u3=0._8
      END if
      IF(l_dxp) then
!        u1=u1+u2
        u2=0._8
!        u3=u3+u4
        u4=0._8
      END if
      IF(l_dzm) then
!        u3=u3+u1
        u1=0._8
!        u4=u4+u2
        u2=0._8
      END if
      IF(l_dzp) then
!        u1=u1+u3
        u3=0._8
!        u2=u2+u4
        u4=0._8
      END if
    END if
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
10 continue

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

subroutine relaxbndrzwguard2(f,rhs,bnd,nr,nz,dr,dz,nc,voltfact,mgparam, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(1:,1:)!rhs(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact, mgparam
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic, nrf, nzi, nzf
REAL(8) :: dt, dt0
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs
REAL(8) :: bndcf0, bndcfxp, bndcfxm, bndcfzp, bndcfzm, bnddt, r, rm, rp, &
           dxx, dzz, dxm, dxp, dzm, dzp, bndphi0xm, bndphi0xp, bndphi0zm, bndphi0zp
IF(l_mgridrz_debug) WRITE(0,*) 'enter relax, level = ',level

! define CFL
dt = mgparam/(2._8/dr**2+2._8/dz**2)
dt0 = mgparam/(4._8/dr**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cf0 = 1._8-2._8*dt/dr**2-2._8*cfz
cfrhs = dt
do j = 2, nr+1
  cfrp(j) = dt * (1._8+0.5_8/REAL(j-1,8)) / dr**2
  cfrm(j) = dt * (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

IF(ixrbnd==dirichlet) then
  nrf=nr-1
else
  nrf=nr
END if
IF(izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(izrbnd==dirichlet) then
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
        f(j,l) = f(j,l) + mgparam*bnd%cnd%dt(ii)*( &
                 bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               - rhs(j,l))
      else
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) then
          dxm = bnd%cnd%dxm(ii)
          dxp = bnd%cnd%dxp(ii)
          dzm = bnd%cnd%dzm(ii)
          dzp = bnd%cnd%dzp(ii)
          select case (bnd_method)
            case (egun)
              dxx=dr
              dzz=dz
            case (ecb)
              dxx=0.5_8*(dxp+dxm)  !ecb
              dzz=0.5_8*(dzp+dzm)  !ecb
            case default
          end select
          r = (j-1)*bndy(i)%dr
          rm = r-0.5_8*dxx
          rp = r+0.5_8*dxx
          bnddt = 1._8/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
          bndcfxm = rm/(r*dxm*dxx)
          bndcfxp = rp/(r*dxp*dxx)
          bndcfzm = 1._8/(dzm*dzz)
          bndcfzp = 1._8/(dzp*dzz)
          bndcf0  = -bndcfxm-bndcfxp-bndcfzm-bndcfzp
          IF(dxm>=dr) then
            bndphi0xm=0._8
          else
            bndphi0xm=bndcfxm*bnd%cnd%volt0xm(ii)
            bndcfxm=0._8
          END if
          IF(dxp>=dr) then
            bndphi0xp=0._8
          else
            bndphi0xp=bndcfxp*bnd%cnd%volt0xp(ii)
            bndcfxp=0._8
          END if
          IF(dzm>=bndy(i)%dz) then
            bndphi0zm=0._8
          else
            bndphi0zm=bndcfzm*bnd%cnd%volt0zm(ii)
            bndcfzm=0._8
          END if
          IF(dzp>=bndy(i)%dz) then
            bndphi0zp=0._8
          else
            bndphi0zp=bndcfzp*bnd%cnd%volt0zp(ii)
            bndcfzp=0._8
          END if
          f(j,l) = f(j,l) + mgparam*bnddt*(bndcf0*f(j,l) &
                 + bndcfxp*f(j+1,l)+bndcfxm*f(j-1,l) &
                 + bndcfzp*f(j,l+1)+bndcfzm*f(j,l-1) &
                 + voltfact*(bndphi0xm+bndphi0xp &
                 + bndphi0zm+bndphi0zp) &
                 - rhs(j,l))
        END if
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
                                 - cfrhs*rhs(j,l)
    end do
    jsw = 3-jsw
  end do
  lsw = 3-lsw
END do !redblack=1, 2

call updateguardcellsrz(f=f, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

#ifdef MPIPARALLEL
  call exchange_fbndz(f,izlbnd,izrbnd)
#endif

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbndrzwguard2

subroutine relaxbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,nc,voltfact,mgparam, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(1:,1:)!rhs(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact, mgparam
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic, nrf, nzi, nzf
REAL(8) :: dt, dt0
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs

IF(l_mgridrz_debug) WRITE(0,*) 'enter relax, level = ',level

! define CFL
dt = mgparam/(2._8/dr**2+2._8/dz**2)
dt0 = mgparam/(4._8/dr**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cf0 = 1._8-2._8*dt/dr**2-2._8*cfz
cfrhs = dt
do j = 2, nr+1
  cfrp(j) = dt * (1._8+0.5_8/REAL(j-1,8)) / dr**2
  cfrm(j) = dt * (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

IF(ixrbnd==dirichlet) then
  nrf=nr-1
else
  nrf=nr
END if
IF(izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(izrbnd==dirichlet) then
  nzf=nz-1
else
  nzf=nz
END if

!  IF(my_index==0) WRITE(0,*) my_index,':e:',f(24,nzf+1)
!  IF(my_index==1) WRITE(0,*) my_index,':e:',f(24,1)
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

    IF(redblack==2) THEN !red
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
        f(j,l) = f(j,l) + mgparam*bnd%cnd%dt(ii)*( &
                 bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               - rhs(j,l))
      else
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*bnd%cnd%dt(ii)*( &
                 bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               - rhs(j,l))
      END if
    ENDDO
  END do
!  IF(my_index==0) WRITE(0,*) my_index,':c:',f(24,nzf+1)
!  IF(my_index==1) WRITE(0,*) my_index,':c:',f(24,1)
!  call exchange_fbndz(f,izlbnd,izrbnd)
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
                                 - cfrhs*rhs(j,l)
!      IF(j==24.and.my_index==0.and.l==nzf+1) &
!         WRITE(0,*)  my_index,':'&
!                                 , f(j+1,l),f(j-1,l)   &
!                                 , f(j,l+1),f(j,l-1) &
!                                 , rhs(j,l),bnd%v(j,l)
!      IF(j==24.and.my_index==1.and.l==1) &
!         WRITE(0,*)  my_index,':'&
!                                 , f(j+1,l),f(j-1,l)   &
!                                 , f(j,l+1),f(j,l-1) &
!                                 , rhs(j,l),bnd%v(j,l)

    end do
    jsw = 3-jsw
  end do
  lsw = 3-lsw
!  IF(my_index==0) WRITE(0,*) my_index,':v:',f(24,nzf+1)
!  IF(my_index==1) WRITE(0,*) my_index,':v:',f(24,1)

#ifdef MPIPARALLEL
!  call check_fbndz(f,izlbnd,izrbnd)
  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
!  call exchange_fbndz(f,izlbnd,izrbnd)
#endif

END do !redblack=1, 2
!  call exchange_fbndz(f,izlbnd,izrbnd)

call updateguardcellsrz(f=f, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

!IF(my_index==0 .AND. f(24,nzf+1)/=0.) stop
!IF(my_index==1 .AND. ANY(f/=0.)) stop

return
END subroutine relaxbndrzwguard

#ifdef MPIPARALLEL
  subroutine exchange_fbndz(f, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd

    INTEGER(ISZ) :: nr,nz
    INTEGER(ISZ) :: p_up, p_down
    integer(ISZ) :: mpi_req(2*nslaves+2),mpistatus(MPI_STATUS_SIZE,2*nslaves+2),mpierror,ir

    
!    write(0,*) my_index,':enter exchangefbnd'
    ir = 0

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nr = SIZE(f,1)-3
    nz = SIZE(f,2)-3

    ! send
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_isend(f(0,2),size(f(:,2)),mpi_double_precision,p_down,0,mpi_comm_world,mpi_req(ir),mpierror)
     end if
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' send to ',p_down,ir
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_isend(f(0,nz),size(f(:,nz)),mpi_double_precision,p_up,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' send to ',p_up,ir

    ! receive
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_irecv(f(0,nz+2),size(f(:,nz+2)),mpi_double_precision,p_up,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if     
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' recv from ',p_up,ir
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_irecv(f(0,0),size(f(:,0)),mpi_double_precision,p_down,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' recv from ',p_down,ir

    if(ir>0) call MPI_WAITALL(ir,mpi_req(1:ir),mpistatus(:,1:ir),mpierror)

  end subroutine exchange_fbndz


  subroutine exchange_fbndz_rb(f, izlbnd, izrbnd, izf)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd, izf

    INTEGER(ISZ) :: nr,nz
    INTEGER(ISZ) :: p_up, p_down
    integer(ISZ) :: mpi_req(4),mpistatus(MPI_STATUS_SIZE,4),mpierror,ir

    real(8), allocatable, dimension(:) :: ftmpd, ftmpu, ftmpds, ftmpus
    
!    write(0,*) my_index,':enter exchangefbnd',izf
    ir = 0

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nz = SIZE(f,2)-3
    if(izf==0) then
      nr = (SIZE(f,1)-1)/2+1
    else
      nr = (SIZE(f,1)-1)/2
    end if

    nr = SIZE(f(izf::2,2))

    allocate(ftmpd(izf:nr+izf-1),ftmpu(izf:nr+izf-1))
    allocate(ftmpds(izf:nr+izf-1),ftmpus(izf:nr+izf-1))

    ! send
    IF(izlbnd<0) then
      ir = ir + 1
!      WRITE(0,*) my_index,SIZE(ftmpds),SIZE(f(izf::2,2))
      ftmpds = f(izf::2,2) 
      call mpi_isend(ftmpds(izf),nr,mpi_double_precision,p_down,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' send to ',p_down,ir
    IF(izrbnd<0) then
      ir = ir + 1
!      WRITE(0,*) my_index,SIZE(ftmpus),SIZE(f(izf::2,nz))
      ftmpus = f(izf::2,nz)
      call mpi_isend(ftmpus(izf),nr,mpi_double_precision,p_up,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' send to ',p_up,ir

    ! receive
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_irecv(ftmpu(izf),nr,mpi_double_precision,p_up,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if     
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' recv from ',p_up,ir
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_irecv(ftmpd(izf),nr,mpi_double_precision,p_down,0,mpi_comm_world,mpi_req(ir),mpierror)
    end if
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' recv from ',p_down,ir

    if(ir>0) call MPI_WAITALL(ir,mpi_req(1),mpistatus(1,1),mpierror)

    IF(izrbnd<0) then
      f(izf::2,nz+2) = ftmpu
    end if     
    IF(izlbnd<0) then
      f(izf::2,0)    = ftmpd
    end if
    deallocate(ftmpd,ftmpu)
    deallocate(ftmpds,ftmpus)
!    write(0,*) my_index,':exit exchangefbnd'

  end subroutine exchange_fbndz_rb
  subroutine check_fbndz(f, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd

    INTEGER(ISZ) :: nr,nz,i
    INTEGER(ISZ) :: p_up, p_down, mpi_req

    REAL(8), DIMENSION(0:SIZE(f,1)-1) :: fr,fl

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nr = SIZE(f,1)-3
    nz = SIZE(f,2)-3

    ! send
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' send to ',p_down
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' send to ',p_up
    IF(izlbnd<0) call mpi_send_real_array(f(:,1), p_down, 0)
    IF(izrbnd<0) call mpi_send_real_array(f(:,nz+1), p_up, 0)

    ! receive
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' recv from ',p_up
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' recv from ',p_down
    IF(izrbnd<0) fr = mpi_recv_real_array(SIZE(f(:,nz)),p_up,0)
    IF(izlbnd<0) fl = mpi_recv_real_array(SIZE(f(:,0 )),p_down,0)

    IF(izrbnd<0) then
     do i = 1, nr+1
      IF(fr(i)/=f(i,nz+1)) WRITE(0,*) 'Error fr mismatch: level ',level,' procs ',my_index,p_up, ' i ',i,fr(i),f(i,nz+1)
     end do
    END if
    IF(izlbnd<0) then
     do i = 1, nr+1
      IF(fl(i)/=f(i,1)) WRITE(0,*) 'Error fl mismatch: level ',level,' procs ',p_down,my_index, ' i ',i,fl(i),f(i,1)
     end do
    END if

    call parallelbarrier()

  end subroutine check_fbndz
  subroutine exchange_resbndz(rho, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: rho(1:,1:)!rho(1:nr+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd

    INTEGER(ISZ) :: nr,nz
    INTEGER(ISZ) :: p_up, p_down
    integer(ISZ) :: mpi_req(4),mpistatus(MPI_STATUS_SIZE,4),mpierror,ir,j
    real(8), allocatable, dimension(:) :: fd, fu

!    write(0,*) my_index,':enter exchangeres'
    ir = 0

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nr = SIZE(rho,1)-1
    nz = SIZE(rho,2)-1

    ! send
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' send to ',p_down
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' send to ',p_up
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_isend_real_array(rho(:,1), p_down, 1,mpi_req(ir))
    end if
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_isend_real_array(rho(:,nz+1), p_up, 1,mpi_req(ir))
    end if

    ! receive
!    IF(bndy(level)%izrbnd<0) WRITE(0,*) my_index, ' recv from ',p_up
!    IF(bndy(level)%izlbnd<0) WRITE(0,*) my_index, ' recv from ',p_down
    IF(izrbnd<0) then
      ir = ir + 1
      allocate(fu(nr+1))
      call mpi_irecv(fu(1),nr+1,mpi_double_precision,p_up,1,mpi_comm_world,mpi_req(ir),mpierror)
    end if
    IF(izlbnd<0) then
      ir = ir + 1
      allocate(fd(nr+1))
      call mpi_irecv(fd(1),nr+1,mpi_double_precision,p_down,1,mpi_comm_world,mpi_req(ir),mpierror)
   end if

!    call parallelbarrier()
    if(ir>0) call MPI_WAITALL(ir,mpi_req(1),mpistatus(1,1),mpierror)
    IF(izrbnd<0) then
      rho(:,nz+1) = rho(:,nz+1) + fu(:)
      deallocate(fu)
    end if
    IF(izlbnd<0) then
      rho(:,1) = rho(:,1) + fd(:)
      deallocate(fd)
    end if

  end subroutine exchange_resbndz
   subroutine merge_work(f, level, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: f(1:,1:)!f(1:nr+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: level, izlbnd, izrbnd

    INTEGER(ISZ) :: nz, p_up, p_down, j, nr
    integer(ISZ) :: mpi_req(2),mpistatus(MPI_STATUS_SIZE,2),mpierror,ir
    real(8), allocatable, dimension(:,:) :: fd, fu

!    write(0,*) my_index,':enter merge'

    ir = 0

!    if(my_index==0) WRITE(0,*) my_index, 'enter merge, level = ',level
    nr     = size(f,1)-1
    nz     = bndy(level-1)%nz
    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    IF(MOD(my_index/bndy(level)%nworkpproc,2)==0) then
    ! send up
      ir = ir +1
      call mpi_isend(f(1,1),(nr+1)*(nz/2+1),mpi_real8,p_up,2,MPI_COMM_WORLD,mpi_req(ir),mpierror)
    ! receive up
      ir = ir +1
      allocate(fu(nr+1,nz/2+1:nz+1))
      call mpi_irecv(fu(1,nz/2+1),(nr+1)*(nz/2+1),mpi_real8,p_up,2,mpi_comm_world,mpi_req(ir),mpierror)
    else
    ! send down
      ir = ir +1
      call mpi_isend(f(1,nz/2+1),int(SIZE(f,1)*(nz/2+1)),mpi_real8,p_down,2,MPI_COMM_WORLD,mpi_req(ir),mpierror)
    ! receive down
      ir = ir +1
      allocate(fd(nr+1,1:nz/2+1))
      call mpi_irecv(fd(1,1),(nr+1)*(nz/2+1),mpi_real8,p_down,2,mpi_comm_world,mpi_req(ir),mpierror)
    END if

    if(ir>0) call MPI_WAITALL(ir,mpi_req(1),mpistatus(1,1),mpierror)

     IF(MOD(my_index/bndy(level)%nworkpproc,2)==0) then
       f(:,nz/2+2:nz+1) = fu(:,nz/2+2:nz+1)
       f(:,nz/2+1) = f(:,nz/2+1)+fu(:,nz/2+1)
       deallocate(fu)
     else
       f(:,1:nz/2) = fd(:,1:nz/2)
       f(:,nz/2+1) = f(:,nz/2+1)+fd(:,nz/2+1)
       deallocate(fd)
     end if

  end subroutine merge_work
#endif

function residbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,voltfact,l_zerolastz, ixrbnd, izlbnd, izrbnd)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN) :: f(0:,0:)!f(0:nx+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(:,:)!rhs(nx+1,nz+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dr, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2) :: residbndrzwguard
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, l, ii, ic, nrf, nzi, nzf
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz

IF(l_mgridrz_debug) WRITE(0,*) 'enter resid, level = ',level

cfz = 1._8 / dz**2
cf0 = -2._8/dr**2-2._8*cfz
do j = 2, nr+1
  cfrp(j) = (1._8+0.5_8/REAL(j-1,8)) / dr**2
  cfrm(j) = (1._8-0.5_8/REAL(j-1,8)) / dr**2
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
IF(izlbnd==dirichlet) then
  nzi=2
else
  nzi=1
END if
IF(izrbnd==dirichlet) then
  nzf=nz-1
else
  nzf=nz
END if

do l = nzi, nzf+1
  j = 1
  IF(bnd%v(j,l)==v_vacuum) &
  residbndrzwguard(j,l) = (cf0-2._8/dr**2) * f(j,l) + 4._8*f(j+1,l)/dr**2   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 - rhs(j,l)
  do j = 2, nrf+1
     IF(bnd%v(j,l)==v_vacuum) &
       residbndrzwguard(j,l) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                                      + cfz*(f(j,l+1)+f(j,l-1)) &
                                      - rhs(j,l)
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
      residbndrzwguard(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
                        + bnd%cnd%cfxp(ii)*f(j+1,l) &
                        + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%phi0xp(ii) &
                        + bnd%cnd%phi0zp(ii)+bnd%cnd%phi0zm(ii)) &
                        - rhs(j,l)
    else
      IF(bnd%v(j,l)==v_bnd.and.bnd%cnd%docalc(ii)) &
      residbndrzwguard(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
                        + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
                        + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%phi0xp(ii)+bnd%cnd%phi0xm(ii) &
                        + bnd%cnd%phi0zp(ii)+bnd%cnd%phi0zm(ii)) &
                        - rhs(j,l)
    END if
  ENDDO
END do

IF(l_zerolastz) residbndrzwguard(:,nz+1) = 0._8

IF(ixrbnd==dirichlet) then
  residbndrzwguard(nr+1,:) = 0.
END if
IF(izlbnd==dirichlet) then
  residbndrzwguard(:,1) = 0.
END if
IF(izrbnd==dirichlet) then
  residbndrzwguard(:,nz+1) = 0.
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit resid, level = ',level

return
end function residbndrzwguard

RECURSIVE subroutine mgbndrzwguard(j, u, rhs, bnd, nr, nz, dr, dz, npre, npost, ncycle, sub, relax_only, npmin, mgparam)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: j, nr, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(:,:), INTENT(IN OUT) :: u
REAL(8), DIMENSION(:,:), INTENT(IN) :: rhs
REAL(8) :: dr, dz, mgparam
TYPE(bndptr) :: bnd(:)
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(:,:), allocatable :: res, v
INTEGER(ISZ) :: i,jj,ll
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
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                        ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
else
  nrnext = bnd(j-1)%nr
  nznext = bnd(j-1)%nz
  drnext = bnd(j-1)%dr
  dznext = bnd(j-1)%dz
  ALLOCATE(res(nrnext+1,nznext+1),v(0:nrnext+2,0:nznext+2))
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                        ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
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
!  IF(restrictwbnd) then
!  res(:,nzresmin:nzresmax) = restrict_wbnd( &
!                             residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
!                             ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd), &
!                             bnd(j),nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
!                             ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
!  else
  IF(bnd(level-1)%l_merged) res=0.
  res(:,nzresmin:nzresmax) = restrict( &
                             residbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
                             ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd), &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
                             ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
!  END if
  call apply_voltage(res,bnd(j-1),0._8)
#ifdef MPIPARALLEL
  IF(bnd(level-1)%l_merged) call merge_work(res,level,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  call exchange_resbndz(rho=res,izlbnd=bnd(j-1)%izlbnd,izrbnd=bnd(j-1)%izrbnd)
#endif
  v = 0.0_8
  do i = 1, ncycle  !(1=V cycles, 2=W cycle)
    call mgbndrzwguard(j=j-1, u=v, rhs=-res, bnd=bnd(1:j-1), nr=nrnext, nz=nznext, dr=drnext, dz=dznext, npre=npre, npost=npost, &
                   ncycle=ncycle, sub=.TRUE., relax_only=.FALSE., npmin=npmin, mgparam=mgparam)
    level = j
  end do
  call apply_voltagewguard(v,bnd(j-1),0._8)
  call updateguardcellsrz(f=v,ixrbnd=ixrbnd,izlbnd=bnd(j-1)%izlbnd,izrbnd=bnd(j-1)%izrbnd)
!  IF(bnd(j)%l_powerof2) then
!    u = u + expandwguard(v(:,nzresmin-1:nzresmax+1))
!  else
 IF(restrictwbnd) then
    u = u + expandwguardandbnd_any(v(:,nzresmin-1:nzresmax+1),bnd(j),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
             ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  else
    u = u + expandwguard_any(v(:,nzresmin-1:nzresmax+1),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
                                ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  END if
!  END if
  call apply_voltagewguard(u,bnd(j),voltf)
  call updateguardcellsrz(f=u,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
#ifdef MPIPARALLEL
  call exchange_fbndz(u,bnd(j)%izlbnd,bnd(j)%izrbnd)
#endif
  call relaxbndrzwguard(f=u,rhs=rhs,bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npost,voltfact=voltf,mgparam=mgparam, &
                        ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  DEALLOCATE(res,v)
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit mg, level = ',level

return
end subroutine mgbndrzwguard

subroutine updateguardcellsrz(f, ixrbnd, izlbnd, izrbnd)
! update guard cells values according to boundary conditions.
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)
INTEGER(ISZ) :: ixrbnd, izlbnd, izrbnd

INTEGER(ISZ) :: ixmax, izmax

ixmax=SIZE(f,1)
izmax=SIZE(f,2)

f(1,:) = f(3,:)
select case (ixrbnd)
    case (dirichlet)
      f(ixmax,:) = f(ixmax-1,:)
    case (neumann)
      f(ixmax,:) = f(ixmax-2,:)
    case default
end select
select case (izlbnd)
    case (dirichlet)
      f(:,1) = f(:,2)
    case (neumann)
      f(:,1) = f(:,3)
    case (periodic)
      f(:,1) = f(:,izmax-1)
    case default
end select
select case (izrbnd)
    case (dirichlet)
      f(:,izmax) = f(:,izmax-1)
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

#ifdef MPIPARALLEL
!  call exchange_fbndz(f,izlbnd,izrbnd)
#endif

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

subroutine solve_multigridrz(grid,accuracy,l_for_timing)
! solve field for u with density rhoinit.
implicit none

! input/output variables
TYPE(grdptr) :: grid
REAL(8), INTENT(IN) :: accuracy  ! required average accuracy
LOGICAL(ISZ) :: l_for_timing

INTEGER(ISZ) :: i, j, nr, nz
REAL(8), allocatable, DIMENSION(:,:) :: uold, uinit
REAL(8) :: maxerr, maxerr_old
LOGICAL :: do_calc, has_diverged

nr = SIZE(grid%phi,1)
nz = SIZE(grid%phi,2)

ALLOCATE(uold(nr,nz))
IF(l_for_timing) then
  ALLOCATE(uinit(nr,nz))
  uinit = grid%phi
END if
do_calc=.true.
has_diverged = .false.

  level = nlevels
  do_calc=.true.
  do while(do_calc)
    maxerr = 1.
    do  j = 1, grid%ncmax
      uold=grid%phi
      call mgbndrzwguard(j=nlevels,u=grid%phi,rhs=-grid%rho/eps0,bnd=grid%bnd,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
                         npre=grid%npre,npost=grid%npost,ncycle=grid%ncycles,sub=.FALSE., relax_only=.false.,npmin=grid%npmin, &
!                         npre=grid%npre,npost=grid%npost,ncycle=grid%ncycles,sub=.FALSE., relax_only=.true.,npmin=grid%npmin, &
!                         npre=grid%npre,npost=grid%npost,ncycle=grid%ncycles,sub=.FALSE., relax_only=.false.,npmin=nlevels-1, &
                         mgparam=grid%mgparam)
      maxerr_old = maxerr
      maxerr = maxval(abs(grid%phi-uold))
#ifdef MPIPARALLEL
      maxerr = mpi_global_compute_real(maxerr,MPI_MAX)
 !     IF(my_index==0) WRITE(0,*) j,maxerr,'*****************************'
 !     IF(my_index==0) then
 !       OPEN(10,FILE='fout.dat',POSITION='append')
 !       WRITE(10,*) j,maxerr
 !       CLOSE(10)
 !     END if
#else
 !     WRITE(0,*) j,maxerr
#endif
      IF(maxerr <= accuracy) then
        do_calc=.false.
        exit
      END if
      IF(maxerr/maxerr_old>=1..and.j>1) then
#ifdef MPIPARALLEL
       IF(my_index==0) then
#endif
        WRITE(0,*) 'WARNING multigridrz, calculation is diverging:'
        WRITE(0,*) '        average initial residue = ',maxerr_old
        WRITE(0,*) '        average current residue = ',maxerr
        WRITE(0,*) '        trying npre and npost = ',grid%npre+1,' (also reset mgparam to 1.8)'
#ifdef MPIPARALLEL
      END if
#endif
        grid%npre  = grid%npre+1
        grid%npost = grid%npost+1
        grid%mgparam = 1.8
        IF(l_for_timing) then
          grid%phi=uinit
        else
          grid%phi=uold
        END if
        GOTO 10
        IF(grid%npmin<i) then
          has_diverged = .true.
          do_calc=.true.
          grid%npmin=grid%npmin+1
          WRITE(0,*) '        trying npmin = ',grid%npmin
        else
          do_calc=.false.
        END if
        10 continue
        IF(l_for_timing) exit
      END if
    end do
    IF(j>=grid%ncmax) do_calc=.false.
  end do

#ifdef MPIPARALLEL
 IF(my_index==0) &
  WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') maxerr,j
#else
  WRITE(0,'("multigridrz: precision = ",e12.5, " after ",i5," iterations...")') maxerr,j
#endif
nb_iters=j
IF(j<grid%ncmax.and..not.do_calc.and.maxerr >= accuracy) nb_iters=grid%ncmax
!IF(levels_min/=npmin) then
!  WRITE(0,'("WARNING multigridrz, levels_min = ",i2,", npmin = ",i2,". Setting levels_min = ",i2,".")') levels_min, npmin, npmin
!  levels_min = npmin
!END if

DEALLOCATE(uold)
IF(l_for_timing) DEALLOCATE(uinit)

return
end subroutine solve_multigridrz

recursive subroutine find_mgparam_rz_1grid(grid)
implicit none
TYPE(grdptr):: grid

REAL(8) :: nexttime, prevtime, prevparam
INTEGER(ISZ) :: npreinit, npostinit

  nlevels = grid%nlevels
  level = nlevels
  bndy => grid%bnd
  ixrbnd = grid%ixrbnd
  izlbnd = grid%izlbnd
  izrbnd = grid%izrbnd

  IF(associated(grid%up)) then
     CALL interpolate_any(unew=grid%phi,uold=grid%up%phi, &
                          nxnew=grid%nr, nznew=grid%nz, &
                          nxold=grid%up%nr, nzold=grid%up%nz, &
                          xminold=grid%up%rmin, xmaxold=grid%up%rmax, &
                          zminold=grid%up%zmin, zmaxold=grid%up%zmax, &
                          xminnew=grid%rmin, xmaxnew=grid%rmax, &
                          zminnew=grid%zmin, zmaxnew=grid%zmax, &
                          ixrbnd=grid%ixrbnd, &
                          izlbnd=grid%izlbnd, &
                          izrbnd=grid%izrbnd, &
                          bnd_only=.true., quad=.true.)
  END if

  npreinit = grid%npre
  npostinit = grid%npost
  ! --- Get initial field solve time
  nexttime = time_field_solve(grid)
  prevtime = 2*nexttime
  ! --- Loop, increasing the number of passes until the time is minimized.
  DO WHILE(nexttime < prevtime)
    prevparam = grid%mgparam
    prevtime = nexttime
    nexttime = ffind_mgparam(grid)
    IF(my_index==0) then
      WRITE(0,*) "Field solve time = ",nexttime
      write(0,*) "mgparam = ",grid%mgparam
      write(0,*) "npre    = ",grid%npre
      write(0,*) "npost   = ",grid%npost
    END if
    IF(nb_iters == grid%ncmax) prevtime=2*nexttime
    IF(nexttime < prevtime) then
      grid%npre  = grid%npre  + 1
      grid%npost = grid%npost + 1
    else
      ! --- Reset the values to the previous ones (which were the best)
      grid%mgparam = prevparam
      grid%npre  = MIN(npreinit,grid%npre-1)
      grid%npost = MIN(npostinit,grid%npost-1)
      ! --- Do some error checking first
      IF(grid%npre  == 0) grid%npre  = 1
      IF(grid%npost == 0) grid%npost = 1
    END if
  END do
  ! --- print error message if maximum iterations is reached.
  IF(nb_iters == grid%ncmax) then
    IF(my_index==0) then
      write(0,*) 'Notice: the maximum number of iterations has been reached, so '
      write(0,*) 'the values above are unlikely to be optimal. Try increasing the '
      write(0,*) 'tolerance, increasing the maximum number of iterations, or making a '
      write(0,*) 'better initial guess of mgparam.'
    END if
  else
    prevtime=findnrecursmin(grid,prevtime)
    IF(my_index==0) then
      write(0,*) "-----------------------------------------"
      write(0,*) "The optimized values:"
      write(0,*) "Field solve time = ",prevtime
      write(0,*) "frz.mgridrz_mgparam     = ",grid%mgparam
      write(0,*) "frz.mgridrz_npre        = ",grid%npre
      write(0,*) "frz.mgridrz_npost       = ",grid%npost
      write(0,*) "frz.mgridrz_levels_min  = ",grid%npmin
    END if
  END if

  IF(associated(grid%down)) call find_mgparam_rz_1grid(grid%down)
  IF(associated(grid%next)) call find_mgparam_rz_1grid(grid%next)

return
END subroutine find_mgparam_rz_1grid

function time_field_solve(grid)
implicit none
REAL(8) :: time_field_solve
TYPE(grdptr):: grid

INTEGER(ISZ) :: ixmax, izmin, izmax
REAL(8) :: beforetime, aftertime
REAL(8), EXTERNAL :: wtime

  ixmax = grid%nr+1
  izmin = 1
  izmax = grid%nz+1
  IF(grid%ixrbnd==dirichlet) ixmax=ixmax-1
  IF(grid%izlbnd==dirichlet) izmin=izmin+1
  IF(grid%izrbnd==dirichlet) izmax=izmax-1

  grid%phi(1:ixmax,izmin:izmax)=0.

  beforetime = wtime()
  call solve_multigridrz(grid=grid, accuracy=mgridrz_accuracy, l_for_timing=.true.)
  aftertime = wtime()
  time_field_solve = aftertime - beforetime
#ifdef MPIPARALLEL
  time_field_solve = mpi_global_compute_real(time_field_solve,MPI_MAX)
#endif

  return
END function time_field_solve

function ffind_mgparam(grid)
implicit none
REAL(8) :: ffind_mgparam
TYPE(grdptr):: grid

INTEGER(ISZ) :: icount, mgiters_prev, up_old, down_old, s
REAL(8) :: mgparam_prev, sincr, mgparam_init
REAL(8), EXTERNAL :: wranf

  icount = 0  ! iteration count

! --- Make sure that mgparam is between 0 and 2.
! --- If mgparam is less then zero, put mgparam closer to 2 since the
! --- optimal value is always closer to 2 than to 0.
  if (grid%mgparam <= 0.) grid%mgparam = max(1., 2. + grid%mgparam)

! --- If mgparam is greater than two, put it on the other side of two
! --- and reduce the increment.  This keeps mgparam near two.
  if (grid%mgparam > 2.) grid%mgparam = max(1., 4. - grid%mgparam)

  mgparam_init=grid%mgparam
 
! --- do initial field solve
  ffind_mgparam = time_field_solve(grid)

! --- set initail values for 'previous' quantities
  mgparam_prev = grid%mgparam
  mgiters_prev = nb_iters

! --- set initial increment for mgparam
  sincr = .05

! --- set mgiters to 0 so that while loop executes at least once
  nb_iters = 0

! --- increment mgparam so next field solve uses new mgparam
  grid%mgparam = grid%mgparam + sincr

! --- Execute while loop until two iterations give the same number of field
! --- solve iterations or until a maximum number of iterations has been
! --- reached.
  do while (mgiters_prev /= nb_iters .and. icount < 200)

!   --- print out current value of mgparam
    IF(my_index==0) write(0,*) "Best parameter so far = ", grid%mgparam

!   --- do field solve (which prints out number of field solve iterations)
    up_old = grid%npre
    down_old = grid%npost
    ffind_mgparam = time_field_solve(grid)

!   --- If field solve took more iterations than previous field solve, change
!   --- direction of the increment and reduce its size.  Reducing its size
!   --- removes the possibility of an infinite loop.
!   --- If a smaller number of field solve iterations was returned, then
!   --- reset the previous values and keep changing mgparam in the same
!   --- direction.  The previous number of iterations is saved in mgiters
!   --- temporarily to check if the two iterations had the same number
!   --- of field solver iterations.
    if (nb_iters > mgiters_prev) then
      sincr = - sincr/2.
      grid%mgparam = mgparam_prev + sincr
    else if (nb_iters < mgiters_prev) then
      s = mgiters_prev
      mgiters_prev = nb_iters
      nb_iters = s
      mgparam_prev = grid%mgparam
      grid%mgparam = mgparam_prev + sincr
    END if

!   --- Make sure that mgparam stays between 0.01 and 2.  .01 is used instead
!   --- of zero since when mgparam is too close to zero, misleading things
!   --- happen.
!   --- If mgparam is outside the range, start the iterations over at a
!   --- random place near the typical optimum value, 1.9. (1.8 for RZ solver)
    if (grid%mgparam <= 0.01 .or. 2. < grid%mgparam) then
      grid%mgparam = 1.8 + wranf()*.05
      sincr = .01
    END if

!   --- increment iteration counter
    icount = icount + 1

    if(grid%npre /= up_old .or. grid%npost/=down_old) then
      IF(my_index==0) write(0,*) "resetting ffind_mgparam"
      icount=0
      grid%npre = up_old
      grid%npost = down_old
      ffind_mgparam = time_field_solve(grid)
      mgparam_prev = mgparam_init
      mgiters_prev = nb_iters
      sincr = .05
!      nb_iters = 0
      grid%mgparam = grid%mgparam + sincr
    END if

  END do

! --- write(0,*) message if an optimal value wasn't found
  if (icount == 200) then
    IF(my_index==0) then
      write(0,*) "Warning: maximum number of iterations reached."
      write(0,*) "         The value of mgparam may not be optimal."
      write(0,*) "         Try increasing mgmaxit."
    END if
  END if

  return
END function ffind_mgparam

function findnrecursmin(grid,prevtime)
!Optimize levels_min, minimizing the fieldsolve time.
implicit none
REAL(8) :: findnrecursmin
TYPE(grdptr) :: grid
REAL(8), INTENT(IN) :: prevtime

REAL(8) :: nexttime, prvtime

  ! --- Get initial field solve time
  nexttime = prevtime
  prvtime = 2*nexttime
  ! --- Loop, increasing the number of passes until the time is minimized.
  do WHILE(nexttime < prvtime .and. grid%npmin < mgridrz_nlevels_max)
    prvtime = nexttime
    grid%npmin = grid%npmin + 1
    nexttime = time_field_solve(grid)
    IF(my_index==0) then
      write(0,*) "Field solve time = ",nexttime
      write(0,*) "frz.mgridrz_levels_min = ",grid%npmin
    END if
    IF(nb_iters == grid%ncmax) prvtime=2*nexttime
    IF(nexttime > prvtime) then
      ! --- Reset the values to the previous ones (which were the best)
      grid%npmin = grid%npmin - 1
      ! --- Do some error checking first
      IF(grid%npmin == 0) grid%npmin = 1
    END if
  END do

  findnrecursmin = prvtime

  return
END function findnrecursmin
END module multigridrz

subroutine multigridrzf(iwhich,u0,rho0,nr0,nz0,dr0,dz0,accuracy)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich, nr0, nz0
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN) :: rho0(nr0+1,nz0+1)
REAL(8), INTENT(IN) :: dr0, dz0, accuracy

  IF(mgridrz_ncmax==0) return

  IF(iwhich==1) return

!  call distribute_rho(basegrid)

  call solve_mgridrz(basegrid,accuracy)

  u0(1:nr0+1,:)=basegrid%phi(1:nr0+1,:)

return
end subroutine multigridrzf

subroutine distribute_rho_rz()
USE multigridrz
implicit none

call distribute_rho(basegrid)

return
END subroutine distribute_rho_rz

RECURSIVE subroutine distribute_rho(grid)
USE multigridrz
implicit none
TYPE(grdptr) :: grid

  IF(associated(grid%down)) call distribute_rho(grid%down)
  IF(associated(grid%next)) call distribute_rho(grid%next)
  IF(associated(grid%up)) then
    call deposit(unew=grid%up%rho, uold=grid%rho, invvolnew=grid%up%invvol, invvolold=grid%invvol, &
                 xminold=grid%rmin, xmaxold=grid%rmax, zminold=grid%zmin, zmaxold=grid%zmax, &
                 xminnew=grid%up%rmin, xmaxnew=grid%up%rmax, zminnew=grid%up%zmin, zmaxnew=grid%up%zmax)
  END if

return
END subroutine distribute_rho

RECURSIVE subroutine solve_mgridrz(grid,accuracy)
USE multigridrz
implicit none
TYPE(grdptr) :: grid
REAL(8), INTENT(IN) :: accuracy

INTEGER(ISZ) :: i

!    grid%mgparam=grid%mgparam+0.05

    IF(associated(grid%up)) then
      CALL interpolate_any(unew=grid%phi,uold=grid%up%phi, &
                           nxnew=grid%nr, nznew=grid%nz, &
                           nxold=grid%up%nr, nzold=grid%up%nz, &
                           xminold=grid%up%rmin, xmaxold=grid%up%rmax, &
                           zminold=grid%up%zmin, zmaxold=grid%up%zmax, &
                           xminnew=grid%rmin, xmaxnew=grid%rmax, &
                           zminnew=grid%zmin, zmaxnew=grid%zmax, &
                           ixrbnd=grid%ixrbnd, &
                           izlbnd=grid%izlbnd, &
                           izrbnd=grid%izrbnd, &
                           bnd_only=.true., quad=.false.)
!                           bnd_only=.true., quad=.true.)
    END if
    nlevels = grid%nlevels
    level = nlevels
    bndy => grid%bnd
    ixrbnd = grid%ixrbnd
    izlbnd = grid%izlbnd
    izrbnd = grid%izrbnd
    call solve_multigridrz(grid=grid, accuracy=accuracy, l_for_timing=.false.)
    IF(associated(grid%down)) call solve_mgridrz(grid%down,accuracy)
    IF(associated(grid%next)) call solve_mgridrz(grid%next,accuracy)

return
END subroutine solve_mgridrz

subroutine find_mgparam_rz()
USE multigridrz
implicit none

    call find_mgparam_rz_1grid(grid=basegrid)

return
END subroutine find_mgparam_rz

subroutine srfrvoutrz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                      xmin,xmax,ymin,ymax,lshell,                      &
                      zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,             &
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)
! call subroutine srfrvout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use Conductor3d
USE multigridrz
implicit none
character(*) rofzfunc
real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
LOGICAL(ISZ):: lfill,lshell,l2symtry,l4symtry
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)
integer(ISZ):: condid

INTEGER(ISZ) :: i,nrc,nzc,igrid
REAL(8) :: drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 bndy => grids_ptr(igrid)%grid%bnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = nlevels,1,-1
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz
  izlbnd = bndy(i)%izlbnd
  izrbnd = bndy(i)%izrbnd
  zmin_in = bndy(i)%zmin!-bndy(i)%dz
  zmax_in = bndy(i)%zmax!+bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

  call srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                      xmin,xmax,ymin,ymax,lshell,                      &
                      zmin_in,zmax_in,zbeam,drc,drc,dzc,nrc,ny,nzc,             &
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

  call addconductors_rz(i,nrc,nzc,drc,dzc,ixrbnd,izlbnd,izrbnd)

 end do
end do
return
end subroutine srfrvoutrz

subroutine srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,   &
                           lzend,xmin,xmax,ymin,ymax,lshell,          &
                           zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,       &
                           ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use Conductor3d
USE multigridrz
implicit none
character(*) rminofz,rmaxofz
real(kind=8):: volt,zmin,zmax,xcent,ycent
LOGICAL(ISZ):: lzend,lshell,l2symtry,l4symtry
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)
integer(ISZ):: condid

INTEGER(ISZ) :: i,nrc,nzc,igrid
REAL(8) :: drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 bndy => grids_ptr(igrid)%grid%bnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = nlevels,1,-1
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz
  izlbnd = bndy(i)%izlbnd
  izrbnd = bndy(i)%izrbnd
  zmin_in = bndy(i)%zmin!-bndy(i)%dz
  zmax_in = bndy(i)%zmax!+bndy(i)%dz

  necndbdy=0
  nocndbdy=0
  ncond = 0

  call srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,  &
                     lzend,xmin,xmax,ymin,ymax,lshell,                &
                     zmin_in,zmax_in,zbeam,drc,drc,dzc,nrc,ny,nzc,             &
                     ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,condid)

  call addconductors_rz(i,nrc,nzc,drc,dzc,ixrbnd,izlbnd,izrbnd)

 end do
end do
return
end subroutine srfrvinoutrz

subroutine setcndtrrz(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz, &
                      l2symtry,l4symtry)
use Conductor3d
USE multigridrz
integer(ISZ):: nx,ny,nz
real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
logical(ISZ):: l2symtry,l4symtry

INTEGER(ISZ) :: i,nrc,nzc,igrid
REAL(8) :: drc,dzc,rmin_in,zmin_in

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 bndy => grids_ptr(igrid)%grid%bnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = nlevels,1,-1
  level = i
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz
  izlbnd = bndy(i)%izlbnd
  izrbnd = bndy(i)%izrbnd
  rmin_in = grids_ptr(igrid)%grid%rmin
  zmin_in = bndy(i)%zmin

  necndbdy=0
  nocndbdy=0
  ncond = 0

  call setcndtr_rz(rmin_in,rmin_in,zmin_in,zbeam,zgrid,nrc,ny,nzc,drc,drc,dzc, &
                   l2symtry,l4symtry)

  call addconductors_rz(i,nrc,nzc,drc,dzc,ixrbnd,izlbnd,izrbnd)

 end do
end do

END subroutine setcndtrrz

subroutine addconductors_rz(i,nrc,nzc,drc,dzc,ixrbnd,izlbnd,izrbnd)
use Conductor3d
USE multigridrz, ONLY: conductor_type, bndy, dirichlet, v_cond, v_bnd, bnd_method, egun, ecb, init_bnd_sublevel
implicit none

INTEGER(ISZ), INTENT(IN) :: nrc,nzc,i,ixrbnd,izlbnd,izrbnd
REAL(8), INTENT(IN) :: drc,dzc

INTEGER(ISZ) :: ii,iii,iv,iiv,nxbnd,nzbndmin,nzbndmax
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz

TYPE(conductor_type), POINTER :: cndpnt

  nxbnd=nrc
  nzbndmin=0
  nzbndmax=nzc
  IF(ixrbnd==dirichlet) nxbnd=nxbnd-1
  IF(izlbnd==dirichlet) nzbndmin=nzbndmin+1
  IF(izrbnd==dirichlet) nzbndmax=nzbndmax-1

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
   IF(iecndx(ii)>=0 .and. iecndx(ii)<=nxbnd .and. iecndz(ii)>=nzbndmin .and. iecndz(ii)<=nzbndmax) then
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
     bndy(i)%cnd%dt(ii) = 1._8/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(ii) = 4._8/(dxp*dxx)
     bndy(i)%cnd%cfzm(ii) = 1._8/(dzm*dzz)
     bndy(i)%cnd%cfzp(ii) = 1._8/(dzp*dzz)
     bndy(i)%cnd%cf0(ii)  = -bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfzm(ii)-bndy(i)%cnd%cfzp(ii)
   else
     r = (bndy(i)%cnd%jj(ii)-1)*bndy(i)%dr
     rm = r-0.5_8*dxx
     rp = r+0.5_8*dxx
     bndy(i)%cnd%dt(ii) = 1._8/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(ii) = rp/(r*dxp*dxx)
     bndy(i)%cnd%cfzm(ii) = 1._8/(dzm*dzz)
     bndy(i)%cnd%cfzp(ii) = 1._8/(dzp*dzz)
     bndy(i)%cnd%cf0(ii)  = -bndy(i)%cnd%cfxm(ii)-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfzm(ii)-bndy(i)%cnd%cfzp(ii)
   END if
     IF(dxm>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xm(ii)=0._8
     else
       bndy(i)%cnd%phi0xm(ii)=bndy(i)%cnd%cfxm(ii)*bndy(i)%cnd%volt0xm(ii)
       bndy(i)%cnd%cfxm(ii)=0._8
     END if
     IF(dxp>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xp(ii)=0._8
     else
       bndy(i)%cnd%phi0xp(ii)=bndy(i)%cnd%cfxp(ii)*bndy(i)%cnd%volt0xp(ii)
       bndy(i)%cnd%cfxp(ii)=0._8
     END if
     IF(dzm>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zm(ii)=0._8
     else
       bndy(i)%cnd%phi0zm(ii)=bndy(i)%cnd%cfzm(ii)*bndy(i)%cnd%volt0zm(ii)
       bndy(i)%cnd%cfzm(ii)=0._8
     END if
     IF(dzp>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zp(ii)=0._8
     else
       bndy(i)%cnd%phi0zp(ii)=bndy(i)%cnd%cfzp(ii)*bndy(i)%cnd%volt0zp(ii)
       bndy(i)%cnd%cfzp(ii)=0._8
     END if
  end do


  do ii = 1, nocndbdy
   iii=necndbdy+ii
   bndy(i)%cnd%jj(iii)  = iocndx(ii)+1
   bndy(i)%cnd%kk(iii)  = iocndz(ii)+1
   IF(bndy(i)%v(bndy(i)%cnd%jj(iii),bndy(i)%cnd%kk(iii))==v_cond) cycle
   IF(iocndx(ii)>=0 .and. iocndx(ii)<=nxbnd .and. iocndz(ii)>=nzbndmin .and. iocndz(ii)<=nzbndmax) then
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
     bndy(i)%cnd%dt(iii) = 1._8/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxp(iii) = 4._8/(dxp*dxx)
     bndy(i)%cnd%cfzm(iii) = 1._8/(dzm*dzz)
     bndy(i)%cnd%cfzp(iii) = 1._8/(dzp*dzz)
     bndy(i)%cnd%cf0(iii)  = -bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfzm(iii)-bndy(i)%cnd%cfzp(iii)
   else
     r = (bndy(i)%cnd%jj(iii)-1)*bndy(i)%dr
     rm = r-0.5_8*dxm
     rp = r+0.5_8*dxp
     bndy(i)%cnd%dt(iii) = 1._8/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(iii) = rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(iii) = rp/(r*dxp*dxx)
     bndy(i)%cnd%cfzm(iii) = 1._8/(dzm*dzz)
     bndy(i)%cnd%cfzp(iii) = 1._8/(dzp*dzz)
     bndy(i)%cnd%cf0(iii)  = -bndy(i)%cnd%cfxm(iii)-bndy(i)%cnd%cfxp(iii)-bndy(i)%cnd%cfzm(iii)-bndy(i)%cnd%cfzp(iii)
   END if
     IF(dxm>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xm(iii)=0._8
     else
       bndy(i)%cnd%phi0xm(iii)=bndy(i)%cnd%cfxm(iii)*bndy(i)%cnd%volt0xm(iii)
       bndy(i)%cnd%cfxm(iii)=0._8
     END if
     IF(dxp>=bndy(i)%dr) then
       bndy(i)%cnd%phi0xp(iii)=0._8
     else
       bndy(i)%cnd%phi0xp(iii)=bndy(i)%cnd%cfxp(iii)*bndy(i)%cnd%volt0xp(iii)
       bndy(i)%cnd%cfxp(iii)=0._8
     END if
     IF(dzm>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zm(iii)=0._8
     else
       bndy(i)%cnd%phi0zm(iii)=bndy(i)%cnd%cfzm(iii)*bndy(i)%cnd%volt0zm(iii)
       bndy(i)%cnd%cfzm(iii)=0._8
     END if
     IF(dzp>=bndy(i)%dz) then
       bndy(i)%cnd%phi0zp(iii)=0._8
     else
       bndy(i)%cnd%phi0zp(iii)=bndy(i)%cnd%cfzp(iii)*bndy(i)%cnd%volt0zp(iii)
       bndy(i)%cnd%cfzp(iii)=0._8
     END if
  end do

end subroutine addconductors_rz

!     ******************************************************************
!     *
!     *                        SUBROUTINE RHOWEIGHTRZ
!     *
!     ******************************************************************


subroutine rhoweightrz(xp,yp,zp,uzp,np,q,nr,nz,dr,dz,zmin)
USE multigridrz
USE constant
USE Subtimers
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
REAL(8), INTENT(IN) :: q, dr, dz, zmin

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr)
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp
REAL(8):: substarttime,wtime

if (lw3dtimesubs) substarttime = wtime()

IF(np==0) return

IF(ngrids>1) then
  call rhoweightrznew(xp,yp,zp,uzp,np,q)
else
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
    IF(uzp(i)==0.) cycle
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr
    zpos = (zp(i)-basegrid%zminp)*invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
#ifdef MPIPARALLEL
    basegrid%rhop(jn, ln)  = basegrid%rhop(jn, ln)  + q * oddr * oddz * invvol(jn-1)
    basegrid%rhop(jnp,ln)  = basegrid%rhop(jnp,ln)  + q *  ddr * oddz * invvol(jnp-1)
    basegrid%rhop(jn, lnp) = basegrid%rhop(jn, lnp) + q * oddr *  ddz * invvol(jn-1)
    basegrid%rhop(jnp,lnp) = basegrid%rhop(jnp,lnp) + q *  ddr *  ddz * invvol(jnp-1)
    basegrid%rhominr = MIN(basegrid%rhominr,jn)
    basegrid%rhomaxr = MAX(basegrid%rhomaxr,jnp)
    basegrid%rhominz = MIN(basegrid%rhominz,ln)
    basegrid%rhomaxz = MAX(basegrid%rhomaxz,lnp)
#else
    basegrid%rho(jn, ln)  = basegrid%rho(jn, ln)  + q * oddr * oddz * invvol(jn-1)
    basegrid%rho(jnp,ln)  = basegrid%rho(jnp,ln)  + q *  ddr * oddz * invvol(jnp-1)
    basegrid%rho(jn, lnp) = basegrid%rho(jn, lnp) + q * oddr *  ddz * invvol(jn-1)
    basegrid%rho(jnp,lnp) = basegrid%rho(jnp,lnp) + q *  ddr *  ddz * invvol(jnp-1)
#endif
  end do
END if

if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
return
END SUBROUTINE RHOWEIGHTRZ

subroutine set_rho_rz(rho,nr,nz,id)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
REAL(8), DIMENSION(nr+1,nz+1), INTENT(IN) :: rho

  IF(id<1.or.id>ngrids) then
    WRITE(0,*) 'Error, id out of bounds'
    WRITE(0,*) 'given id = ',id,' while 1 < id < ',ngrids
    return
  END if
  IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
    WRITE(0,*) 'Error, dimensions should be the same: '
    WRITE(0,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2)
    WRITE(0,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
    return
  END if
  grids_ptr(id)%grid%rho=rho

return
end subroutine set_rho_rz

subroutine get_rho_rz(rho,nr,nz,id,rhop)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz,rhop
REAL(8), DIMENSION(nr+1,nz+1), INTENT(IN OUT) :: rho

  IF(id<1.or.id>ngrids) then
    WRITE(0,*) 'Error, id out of bounds'
    WRITE(0,*) 'given id = ',id,' while 1 < id < ',ngrids
    return
  END if
#ifdef MPIPARALLEL
  IF(rhop==1) then
    IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rhop,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rhop,2)) then
      WRITE(0,*) 'Error, dimensions should be the same: '
      WRITE(0,*) 'Nr, Nz for rhop    : ',SIZE(rho,1),SIZE(rho,2)
      WRITE(0,*) 'Nr, Nz for rhop(id): ',SIZE(grids_ptr(id)%grid%rhop,1),SIZE(grids_ptr(id)%grid%rhop,2)
      return
    END if
    rho=grids_ptr(id)%grid%rhop
  else
#endif
    IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
      WRITE(0,*) 'Error, dimensions should be the same: '
      WRITE(0,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2)
      WRITE(0,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
      return
    END if
    rho=grids_ptr(id)%grid%rho
#ifdef MPIPARALLEL
  END if
#endif
return
end subroutine get_rho_rz

subroutine rhoweightrznew(xp,yp,zp,uzp,np,q)
USE constant
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
REAL(8), INTENT(IN) :: q

REAL(8) :: rpos, zpos, ddr, ddz, oddr, oddz
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid

REAL(8), DIMENSION(:), ALLOCATABLE :: invdr, invdz, zmin

ALLOCATE(invdr(ngrids),invdz(ngrids),zmin(ngrids))

do igrid = 1, ngrids
  invdr(igrid) = grids_ptr(igrid)%grid%invdr
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  zmin (igrid) = grids_ptr(igrid)%grid%zminp
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    igrid = 1
    grids=>basegrid
    ingrid=.false.
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr(igrid)
    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part(jn,ln)
        grids=>grids_ptr(igrid)%grid
        rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr(igrid)
        zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
#ifdef MPIPARALLEL
    grids%rhop(jn, ln)  = grids%rhop(jn, ln)  + q * oddr * oddz * grids%invvol(jn)
    grids%rhop(jnp,ln)  = grids%rhop(jnp,ln)  + q *  ddr * oddz * grids%invvol(jnp)
    grids%rhop(jn, lnp) = grids%rhop(jn, lnp) + q * oddr *  ddz * grids%invvol(jn)
    grids%rhop(jnp,lnp) = grids%rhop(jnp,lnp) + q *  ddr *  ddz * grids%invvol(jnp)
#else
    grids%rho(jn, ln)  = grids%rho(jn, ln)  + q * oddr * oddz * grids%invvol(jn)
    grids%rho(jnp,ln)  = grids%rho(jnp,ln)  + q *  ddr * oddz * grids%invvol(jnp)
    grids%rho(jn, lnp) = grids%rho(jn, lnp) + q * oddr *  ddz * grids%invvol(jn)
    grids%rho(jnp,lnp) = grids%rho(jnp,lnp) + q *  ddr *  ddz * grids%invvol(jnp)
#endif
  end do

  DEALLOCATE(invdr,invdz,zmin)

  return
END SUBROUTINE RHOWEIGHTRZNEW

subroutine reset_rzmgrid_rho()
USE multigridrz
implicit none
INTEGER(ISZ) :: igrid

  do igrid = 1, ngrids
#ifdef MPIPARALLEL
    grids_ptr(igrid)%grid%rhop=0.
    grids_ptr(igrid)%grid%rhominr = grids_ptr(igrid)%grid%nr+2
    grids_ptr(igrid)%grid%rhomaxr = -1
    grids_ptr(igrid)%grid%rhominz = grids_ptr(igrid)%grid%nz+2
    grids_ptr(igrid)%grid%rhomaxz = -1
#else
    grids_ptr(igrid)%grid%rho=0.
#endif
  end do

return
end subroutine reset_rzmgrid_rho

!     ******************************************************************
!     *
!     *                        SUBROUTINE RHOWEIGHTRZ_DEFORM
!     *
!     ******************************************************************


subroutine rhoweightrz_deform(xp,yp,zp,uzp,np,q,rho,nr,nz,dr,dz,zmin,xfact,yfact)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
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
    IF(uzp(i)==0.) cycle
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

  IF(ixrbnd==neumann) rho(nr,:) = 2._8*rho(nr,:)
  IF(izlbnd==neumann) rho(:,0)  = 2._8*rho(:,0)
  IF(izrbnd==neumann) rho(:,nz) = 2._8*rho(:,nz)

  return
end subroutine rhobndrz


subroutine gchange_rhop_phip_rz()
USE multigridrz
implicit none
#ifdef MPIPARALLEL
INTEGER(ISZ) :: nzp

 IF(nzpslave(my_index)/=basegrid%nzp) then
  nzp=nzpslave(my_index)
  DEALLOCATE(basegrid%rhop,basegrid%phip)
  ALLOCATE(basegrid%rhop(basegrid%nr+1,nzp+1),basegrid%phip(0:basegrid%nr+2,0:nzp+2))
  basegrid%nzp=nzp
  basegrid%phip=0.
  basegrid%rhop=0.
 END if

  basegrid%zminp=zpslmin(0)+izpslave(my_index)*basegrid%dz
  call get_phip_from_phi(basegrid)
#endif
return
end subroutine gchange_rhop_phip_rz

#ifdef MPIPARALLEL

subroutine getrhoforfieldsolverz(nr,nz,rho)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: nr,nz
REAL(8), DIMENSION(0:nr,0:nz), INTENT(IN OUT) :: rho

  call get_rho_from_rhop(basegrid)

  return
end subroutine getrhoforfieldsolverz

subroutine get_rho_from_rhop(grid)
USE multigridrz
implicit none

TYPE(grdptr) :: grid

INTEGER(ISZ) :: i, ilg, iug, ilp, iup, il, iu, j, l, ii, ll, jex, lex

integer(ISZ):: mpistatus(MPI_STATUS_SIZE)
INTEGER(ISZ):: mpirequest,mpierror

logical(ISZ) :: testthis=.false.
INTEGER(ISZ), ALLOCATABLE :: wz(:)

grid%rho = 0.

IF(testthis) then
i=my_index
ilp = 1+izpslave(i)
iup = 1+izpslave(i)+nzpslave(i)
do l = ilp, iup
 do j = 1, basegrid%nr+1
   ll = l-ilp+1
   grid%rhop(j,ll) = l+j*10000
 end do
end do
IF(ANY(grid%rhop(:,:)==0.)) then
  WRITE(0,*) 'rhop = 0.'
  call abort()
END if
END if
! send slices of rhop to processors that need it
i=my_index
ilp = 1+izpslave(i)
iup = 1+izpslave(i)+nzpslave(i)
do i = 0, nslaves-1
  if (i/=my_index) then
    ilg = 1+izslave(i)
    iug = 1+izslave(i) +nzslave(i)
    il = MAX(ilg,ilp)-ilp+1
    iu = MIN(iug,iup)-ilp+1
    IF(il>iu) cycle
    call mpi_isend(grid%rhop(1,il),SIZE(grid%rhop(:,il:iu)),mpi_double_precision,i,0,MPI_COMM_WORLD,mpirequest,ierr)
  end if
end do

! recv slices of rhop from required processors
i=my_index
ilg = 1+izslave(i)
iug = 1+izslave(i) +nzslave(i)
do i = 0, nslaves-1
  ilp = 1+izpslave(i)
  iup = 1+izpslave(i)+nzpslave(i)
  il = MAX(ilg,ilp)-ilg+1
  iu = MIN(iug,iup)-ilg+1
  IF(il>iu) cycle
  if (i==my_index) then
    grid%rho(:,il:iu) = grid%rho(:,il:iu) + grid%rhop(:,il+ilg-ilp:iu+ilg-ilp)
  else
    grid%rho(:,il:iu) = grid%rho(:,il:iu) &
                        + RESHAPE(mpi_recv_real_array(SIZE(grid%rho(:,il:iu)), i, 0) &
                        , SHAPE(grid%rho(:,il:iu)))
  end if
end do

!call parallelbarrier()

IF(testthis) then
ALLOCATE(wz(izpslave(nslaves-1)+nzpslave(nslaves-1)+1))
wz = 0
do i = 0, nslaves-1
  ilp = 1+izpslave(i)
  iup = 1+izpslave(i)+nzpslave(i)
  wz(ilp:iup) = wz(ilp:iup)+1
END do
i=my_index
ilg = 1+izslave(i)
iug = 1+izslave(i)+nzslave(i)
ilp = 1+izpslave(i)
iup = 1+izpslave(i)+nzpslave(i)
do l = ilg, iug
 do j = 1, basegrid%nr+1
   ll = l-ilg+1
      jex = INT(grid%rho(j,ll)/10000)
      lex = INT(grid%rho(j,ll))-jex*10000
!      IF(jex/wz(l-ilg+ilp)/=j.or.lex/wz(l-ilg+ilp)/=l) WRITE(0,*) my_index,':',j,l,l-ilg+ilp,grid%rho(j,ll),wz(l-ilg+ilp)
      IF(jex/wz(l)/=j.or.lex/wz(l)/=l) WRITE(0,*) my_index,':',j,l,grid%rho(j,ll),wz(l)
 end do
end do
  DEALLOCATE(wz)
 call abort()
endif

end subroutine get_rho_from_rhop

subroutine getphiforparticlesrz()
USE multigridrz
implicit none

  call get_phip_from_phi(basegrid)

end subroutine getphiforparticlesrz

subroutine get_phip_from_phi(grid)
USE multigridrz
implicit none

TYPE(grdptr) :: grid

INTEGER(ISZ) :: i, ilg, iug, ilp, iup, il, iu, j,l, ll, jex, lex, testeq
integer(ISZ):: mpistatus(MPI_STATUS_SIZE)
INTEGER(ISZ):: mpirequest,mpierror

logical(ISZ) :: testthis=.false.

grid%phip = 0.

IF(testthis) then
i=my_index
ilg = 1+izslave(i)-1
iug = 1+izslave(i)+nzslave(i)+1
do l = ilg, iug
 do j = 0, basegrid%nr+2
   ll = l-ilg
   grid%phi(j,ll) = l+j*10000
 end do
end do
END if

! send slices of phi to processors that need it
i=my_index
ilg = 1+izslave(i)             - 1
iug = 1+izslave(i) +nzslave(i) + 1
do i = 0, nslaves-1
  if (i/=my_index) then
    ilp = 1+izpslave(i)             - 1
    iup = 1+izpslave(i)+nzpslave(i) + 1
    il = MAX(ilg,ilp)-ilg
    iu = MIN(iug,iup)-ilg
    IF(il>iu) cycle
    call mpi_isend(grid%phi(0,il),SIZE(grid%phi(:,il:iu)),mpi_double_precision,i,0,MPI_COMM_WORLD,mpirequest,ierr)
  end if
end do

! recv slices of phi from required processors
i=my_index
ilp = 1+izpslave(i)             - 1
iup = 1+izpslave(i)+nzpslave(i) + 1
do i = 0, nslaves-1
  ilg = 1+izslave(i)             - 1
  iug = 1+izslave(i) +nzslave(i) + 1
  il = MAX(ilg,ilp)-ilp
  iu = MIN(iug,iup)-ilp
  IF(il>iu) cycle
  if (i==my_index) then
    grid%phip(:,il:iu) = grid%phi(:,il+ilp-ilg:iu+ilp-ilg)
  else
    grid%phip(:,il:iu) = RESHAPE(mpi_recv_real_array(SIZE(grid%phip(:,il:iu)), i, 0) &
                        ,SHAPE(grid%phip(:,il:iu)))
  end if
end do

!call parallelbarrier()

IF(testthis) then
i=my_index
ilp = izpslave(i)
iup = izpslave(i)+nzpslave(i)
do l = ilp, iup
 do j = 0, basegrid%nr+2
   ll = l-ilp
      jex = INT(grid%phip(j,ll)/10000)
      lex = INT(grid%phip(j,ll))-jex*10000
      testeq = 0
      IF(jex/=j.or.lex/=l) WRITE(0,*) my_index,':',j,l,grid%phip(j,ll)
 end do
end do
 call abort()
endif

end subroutine get_phip_from_phi

#endif

subroutine fieldweightrzold(xp,yp,zp,uzp,ex,ey,ez,np,phi,e,nr,nz,dr,dz,zmin,calcselfe)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
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
    IF(uzp(i)==0.) cycle
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
end subroutine fieldweightrzold

subroutine fieldweightrz(xp,yp,zp,uzp,ex,ey,ez,np)
USE multigridrz
USE constant
USE subtimers
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez

REAL(8) :: rpos, zpos, invrpos, ddr, ddz, oddr, oddz, er
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid

REAL(8):: substarttime,wtime
if (lw3dtimesubs) substarttime = wtime()

IF(ngrids>1) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    igrid = 1
    grids => basegrid
    ingrid=.false.
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*grids%invdr
    zpos = (zp(i)-grids%zminp)*grids%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part_field_dep(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part_field_dep(jn,ln)
        grids=>grids_ptr(igrid)%grid
        rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*grids%invdr
        zpos = (zp(i)-grids%zminp)*grids%invdz
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    er = 0.5*(oddr * oddz * (grids%phip(jn-1,ln  )-grids%phip(jn+1,ln  ))  &
            + ddr  * oddz * (grids%phip(jn  ,ln  )-grids%phip(jn+2,ln  ))  &
            + oddr * ddz  * (grids%phip(jn-1,ln+1)-grids%phip(jn+1,ln+1))  &
            + ddr  * ddz  * (grids%phip(jn  ,ln+1)-grids%phip(jn+2,ln+1)))*grids%invdr
#else
    er = 0.5*(oddr * oddz * (grids%phi(jn-1,ln  )-grids%phi(jn+1,ln  ))  &
            + ddr  * oddz * (grids%phi(jn  ,ln  )-grids%phi(jn+2,ln  ))  &
            + oddr * ddz  * (grids%phi(jn-1,ln+1)-grids%phi(jn+1,ln+1))  &
            + ddr  * ddz  * (grids%phi(jn  ,ln+1)-grids%phi(jn+2,ln+1)))*grids%invdr
#endif
    IF(rpos>1.e-10) then
      invrpos=grids%invdr/rpos
      ex(i) = er*xp(i)*invrpos
      ey(i) = er*yp(i)*invrpos
    else
      ex(i) = er
      ey(i) = 0._8
    END if
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddr * oddz * (grids%phip(jn  ,ln-1)-grids%phip(jn  ,ln+1))  &
               + ddr  * oddz * (grids%phip(jn+1,ln-1)-grids%phip(jn+1,ln+1))  &
               + oddr * ddz  * (grids%phip(jn  ,ln  )-grids%phip(jn  ,ln+2))  &
               + ddr  * ddz  * (grids%phip(jn+1,ln  )-grids%phip(jn+1,ln+2)))*grids%invdz
#else
    ez(i) = 0.5*(oddr * oddz * (grids%phi(jn  ,ln-1)-grids%phi(jn  ,ln+1))  &
               + ddr  * oddz * (grids%phi(jn+1,ln-1)-grids%phi(jn+1,ln+1))  &
               + oddr * ddz  * (grids%phi(jn  ,ln  )-grids%phi(jn  ,ln+2))  &
               + ddr  * ddz  * (grids%phi(jn+1,ln  )-grids%phi(jn+1,ln+2)))*grids%invdz
#endif
  END do
else
  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    ingrid=.false.
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*basegrid%invdr
    zpos = (zp(i)-basegrid%zminp)*basegrid%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    er = 0.5*(oddr * oddz * (basegrid%phip(jn-1,ln  )-basegrid%phip(jn+1,ln  ))  &
            + ddr  * oddz * (basegrid%phip(jn  ,ln  )-basegrid%phip(jn+2,ln  ))  &
            + oddr * ddz  * (basegrid%phip(jn-1,ln+1)-basegrid%phip(jn+1,ln+1))  &
            + ddr  * ddz  * (basegrid%phip(jn  ,ln+1)-basegrid%phip(jn+2,ln+1)))*basegrid%invdr
#else
    er = 0.5*(oddr * oddz * (basegrid%phi(jn-1,ln  )-basegrid%phi(jn+1,ln  ))  &
            + ddr  * oddz * (basegrid%phi(jn  ,ln  )-basegrid%phi(jn+2,ln  ))  &
            + oddr * ddz  * (basegrid%phi(jn-1,ln+1)-basegrid%phi(jn+1,ln+1))  &
            + ddr  * ddz  * (basegrid%phi(jn  ,ln+1)-basegrid%phi(jn+2,ln+1)))*basegrid%invdr
#endif
    IF(rpos>1.e-10) then
      invrpos=basegrid%invdr/rpos
      ex(i) = er*xp(i)*invrpos
      ey(i) = er*yp(i)*invrpos
    else
      ex(i) = er
      ey(i) = 0._8
    END if
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddr * oddz * (basegrid%phip(jn  ,ln-1)-basegrid%phip(jn  ,ln+1))  &
               + ddr  * oddz * (basegrid%phip(jn+1,ln-1)-basegrid%phip(jn+1,ln+1))  &
               + oddr * ddz  * (basegrid%phip(jn  ,ln  )-basegrid%phip(jn  ,ln+2))  &
               + ddr  * ddz  * (basegrid%phip(jn+1,ln  )-basegrid%phip(jn+1,ln+2)))*basegrid%invdz
#else
    ez(i) = 0.5*(oddr * oddz * (basegrid%phi(jn  ,ln-1)-basegrid%phi(jn  ,ln+1))  &
               + ddr  * oddz * (basegrid%phi(jn+1,ln-1)-basegrid%phi(jn+1,ln+1))  &
               + oddr * ddz  * (basegrid%phi(jn  ,ln  )-basegrid%phi(jn  ,ln+2))  &
               + ddr  * ddz  * (basegrid%phi(jn+1,ln  )-basegrid%phi(jn+1,ln+2)))*basegrid%invdz
#endif
  END do
END if

  if (lw3dtimesubs) timesete3d = timesete3d + wtime() - substarttime
  return
end subroutine fieldweightrz

subroutine fieldweightrz_deform(xp,yp,zp,uzp,ex,ey,ez,np,phi,nr,nz,dr,dz,zmin,xfact,yfact,calcphi,phi3d,selfe)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
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
    IF(uzp(i)==0.) cycle
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

subroutine fieldweightrz_deform_old(xp,yp,zp,uzp,ex,ey,ez,np,phi,e,nr,nz,dr,dz,zmin,xfact,yfact,calcselfe)
USE constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
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
    IF(uzp(i)==0.) cycle
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
USE Particles, ONLY: xp, yp, zp, uzp
use FRZmgrid
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

  if(.not.mgridrz_deform) return

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
      IF(uzp(i)==0.) cycle
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

subroutine setphirz(np,xp,yp,zp,p)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: p

REAL(8) :: rpos, zpos, ddr, ddz, oddr, oddz
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid

! Collect phi using linear interpolation

IF(ngrids>1) then

  do i = 1, np
    igrid = 1
    grids => basegrid
    ingrid=.false.
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*grids%invdr
    zpos = (zp(i)-grids%zmin)*grids%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part(jn,ln)
        grids=>grids_ptr(igrid)%grid
        rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*grids%invdr
        zpos = (zp(i)-grids%zminp)*grids%invdz
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    p(i) = oddr * oddz * grids%phi(jn,  ln  )  &
         + ddr  * oddz * grids%phi(jn+1,ln  )  &
         + oddr * ddz  * grids%phi(jn,  ln+1)  &
         + ddr  * ddz  * grids%phi(jn+1,ln+1)
  END do
else
  do i = 1, np
    ingrid=.false.
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*basegrid%invdr
    zpos = (zp(i)-basegrid%zmin)*basegrid%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    p(i) = oddr * oddz * basegrid%phi(jn,  ln  )  &
         + ddr  * oddz * basegrid%phi(jn+1,ln  )  &
         + oddr * ddz  * basegrid%phi(jn,  ln+1)  &
         + ddr  * ddz  * basegrid%phi(jn+1,ln+1)
  END do
END if

  return
end subroutine setphirz

subroutine init_base(nr,nz,dr,dz,rmin,zmin)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: nr,nz
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin

 call init_basegrid(nr,nz,dr,dz,rmin,zmin)

 return
END subroutine init_base

subroutine add_subgrid(id,nr,nz,dr,dz,rmin,zmin,guard_min_r,guard_max_r,guard_min_z,guard_max_z)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz,guard_min_r,guard_max_r,guard_min_z,guard_max_z
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin

  IF(id<1 .or. id>ngrids) then
    WRITE(0,*) 'Fatal error in add_subgrid: id = ', id ,' WHILE id = (1,...,',ngrids,')'
    stop
  END if

  call add_grid(grids_ptr(id)%grid,nr,nz,dr,dz,rmin,zmin,guard_min_r,guard_max_r,guard_min_z,guard_max_z)

return
END subroutine add_subgrid

subroutine get_phi_subgrid(id,phi,nr,nz)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
REAL(8), DIMENSION(1:nr+1,1:nz+1) :: phi

  IF(id<1 .or. id>ngrids) then
    WRITE(0,*) 'Error in get_phi_subgrid: id = ', id ,' WHILE id = (1,...,',ngrids,')'
    WRITE(0,*) 'Returning Phi=0'
    phi(1:nr+1,1:nz+1) = 0.
  else
    phi(1:nr+1,1:nz+1) = grids_ptr(id)%grid%phi(1:nr+1,1:nz+1)
  END if

return
END subroutine get_phi_subgrid

subroutine get_array_subgrid(id,array,nr,nz,which)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
CHARACTER(*) :: which
REAL(8), DIMENSION(1:nr+1,1:nz+1) :: array

  IF(id<1 .or. id>ngrids) then
    WRITE(0,*) 'Error in get_phi_subgrid: id = ', id ,' WHILE id = (1,...,',ngrids,')'
    WRITE(0,*) 'Returning Array=0'
    array(1:nr+1,1:nz+1) = 0.
  else
    select case (which)
      case ("rho","r")
        array = grids_ptr(id)%grid%rho(1:nr+1,1:nz+1)
      case ("phi","p")
        array = grids_ptr(id)%grid%phi(1:nr+1,1:nz+1)
      case ("loc_part","lp")
        array = grids_ptr(id)%grid%loc_part(1:nr+1,1:nz+1)
      case ("loc_part_field_dep","lpfd")
        array = grids_ptr(id)%grid%loc_part_field_dep(1:nr+1,1:nz+1)
      case("res","RES")
        array = residbndrzwguard(f=grids_ptr(id)%grid%phi, &
                                 rhs=grids_ptr(id)%grid%rho, &
                                 bnd=grids_ptr(id)%grid%bnd(grids_ptr(id)%grid%nlevels), &
                                 nr=grids_ptr(id)%grid%nr, &
                                 nz=grids_ptr(id)%grid%nz, &
                                 dr=grids_ptr(id)%grid%dr, &
                                 dz=grids_ptr(id)%grid%dz, &
                                 voltfact=1., &
                                 l_zerolastz=.FALSE., &
                                 ixrbnd=ixrbnd, &
                                 izlbnd=grids_ptr(id)%grid%bnd(grids_ptr(id)%grid%nlevels)%izlbnd, &
                                 izrbnd=grids_ptr(id)%grid%bnd(grids_ptr(id)%grid%nlevels)%izrbnd)
      case default
        WRITE(0,*) which,' is not a valid option for get_array_subgrid.'
        WRITE(0,*) 'Valid options are: '&
                 //'  rho [abrv:r]' &
                 //'  phi [abrv:p]' &
                 //'  loc_part [abrv:lp]' &
                 //'  loc_part_field_dep [abrv:lpfd]'
    end select
  END if

return
END subroutine get_array_subgrid


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
          WRITE(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp,bndy(i)%cnd%phi0zm,bndy(i)%cnd%phi0zp,bndy(i)%cnd%dt
          WRITE(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, bndy(i)%cnd%cfzp, bndy(i)%cnd%cfzm
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
          read(10,*) bndy(i)%cnd%phi0xm,bndy(i)%cnd%phi0xp,bndy(i)%cnd%phi0zm,bndy(i)%cnd%phi0zp,bndy(i)%cnd%dt
          read(10,*) bndy(i)%cnd%cf0, bndy(i)%cnd%cfxp, bndy(i)%cnd%cfxm, bndy(i)%cnd%cfzp, bndy(i)%cnd%cfzm
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

subroutine get_cond_rz(igrid,ilevel)
USE multigridrz
USE Conductor3d
implicit none
INTEGER :: igrid,ilevel

INTEGER :: i,ic,icc,ice,ico
TYPE(bndptr), pointer :: bnd

 bnd => grids_ptr(igrid)%grid%bnd(grids_ptr(igrid)%grid%nlevels-ilevel+1)

 ncond    = 0
 necndbdy = 0
 nocndbdy = 0
 do ic = 1, bnd%nb_conductors
   IF(ic==1) then
     bnd%cnd => bnd%first
   else
     bnd%cnd => bnd%cnd%next
   END if
   ncond    = ncond    + bnd%cnd%ncond
   necndbdy = necndbdy + bnd%cnd%nbbndred
   nocndbdy = nocndbdy + bnd%cnd%nbbnd-bnd%cnd%nbbndred
 END do
 ncondmax = ncond
 ncndmax = max(necndbdy,nocndbdy)
 call gchange("Conductor3d",0)

 icc=0
 ice=0
 ico=0
 do ic = 1, bnd%nb_conductors
   IF(ic==1) then
     bnd%cnd => bnd%first
   else
     bnd%cnd => bnd%cnd%next
   END if
   do i = 1, bnd%cnd%ncond
     icc=icc+1
     ixcond(icc) = bnd%cnd%jcond(i)-1
     izcond(icc) = bnd%cnd%kcond(i)-1
     icondlevel(icc) = ic
   end do
   do i = 1, bnd%cnd%nbbndred
    IF(bnd%v(bnd%cnd%jj(i),bnd%cnd%kk(i))==v_bnd) then
     ice=ice+1
     iecndx(ice) = bnd%cnd%jj(i)-1
     iecndz(ice) = bnd%cnd%kk(i)-1
     ecdelmx(ice) = bnd%cnd%dxm(i)/bnd%dr
     ecdelpx(ice) = bnd%cnd%dxp(i)/bnd%dr
     ecdelmz(ice) = bnd%cnd%dzm(i)/bnd%dz
     ecdelpz(ice) = bnd%cnd%dzp(i)/bnd%dz
     iecndlevel(ice) = ic
    END if
   end do
   do i = bnd%cnd%nbbndred+1, bnd%cnd%nbbnd
    IF(bnd%v(bnd%cnd%jj(i),bnd%cnd%kk(i))==v_bnd) then
     ico=ico+1
     iocndx(ico) = bnd%cnd%jj(i)-1
     iocndz(ico) = bnd%cnd%kk(i)-1
     ocdelmx(ico) = bnd%cnd%dxm(i)/bnd%dr
     ocdelpx(ico) = bnd%cnd%dxp(i)/bnd%dr
     ocdelmz(ico) = bnd%cnd%dzm(i)/bnd%dz
     ocdelpz(ico) = bnd%cnd%dzp(i)/bnd%dz
     iocndlevel(ico) = ic
    END if
   end do
 END do
 necndbdy = ice
 nocndbdy = ico

return
end subroutine get_cond_rz

