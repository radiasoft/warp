!     Last change:  JLV   5 Nov 2002    2:17 pm
#include "top.h"

module multigrid_common

USE constant
USE PSOR3d, ONLY:boundxy,bound0,boundnz
USE InGen3d, ONLY:solvergeom,RZgeom,XYZgeom,XZgeom,Zgeom,l2symtry,l4symtry
USE FRZmgrid
#ifdef MPIPARALLEL
  use Parallel
  use mpirz
#endif

LOGICAL :: l_mgridrz_debug = .false., &
           l_timing_rz     = .TRUE.,  &
           restrictwbnd    = .true., &
           vlocs           = .FALSE., &
           l_print_timing  = .false., &
           l_jump          = .false.

REAL(8) :: t_relax, t_restrict, t_expand, t_before, t_apply_voltage, t_updateguard, t_allocate
REAL(8), EXTERNAL :: wtime

REAL(8) :: inveps0

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
  ! conductor ID
  INTEGER(ISZ), POINTER, DIMENSION(:) :: condid
  ! stencil coefficients for relaxation iteration
  REAL(8), POINTER, DIMENSION(:) :: cf0, cfxp, cfxm, cfyp, cfym, cfzp, cfzm, dt
  REAL(8), POINTER, DIMENSION(:) :: phi0xm, phi0xp, phi0ym, phi0yp, phi0zm, phi0zp
  ! voltages on conductors near grid
  REAL(8), POINTER, DIMENSION(:) :: volt0xm, volt0xp, volt0ym, volt0yp, volt0zm, volt0zp
  ! IDs of conductors
  INTEGER(ISZ), POINTER, DIMENSION(:) :: condidxm, condidxp, condidym, condidyp, condidzm, condidzp
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
INTEGER(ISZ) :: ixlbndi=dirichlet, ixrbndi=dirichlet, &
                iylbndi=dirichlet, iyrbndi=dirichlet, &
                izlbndi=dirichlet, izrbndi=dirichlet  ! boundary conditions for each side
INTEGER(ISZ) :: ixlbnd=dirichlet, ixrbnd=dirichlet, &
                iylbnd=dirichlet, iyrbnd=dirichlet, &
                izlbnd=dirichlet, izrbnd=dirichlet  ! boundary conditions for each side

INTEGER(ISZ), parameter :: v_vacuum=0, v_cond=1, v_bnd=2, v_dirichlet=3
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
  INTEGER(ISZ), POINTER, DIMENSION(:) :: vlocs_j, vlocs_k, lshift
  INTEGER(ISZ) :: nvlocs, nvlocsred
  LOGICAL(ISZ) :: l_lshift
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
  REAL(8) :: rmin,rmax,xmin,xmax,zmin,zmax,dr,dz,invdr,invdz,zminp
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
TYPE(grd_ptr), DIMENSION(:), ALLOCATABLE :: grids_ptr, gridinit

INTEGER(ISZ) :: grids_nids,     & ! number of grid IDs
                n_avail_ids,    & ! number of available IDs
                avail_ids(100), & ! array of grid IDs
                level_del_grid

INTEGER(ISZ) :: nguardx=1, nguardy=1, nguardz=1 ! number of guard cells in each direction


contains

subroutine init_basegrid(nr,nz,dr,dz,rmin,zmin)
USE Multigrid3d
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin
INTEGER(ISZ) :: i,j, nzp, mglevel

  inveps0 = 1./eps0

  IF(associated(basegrid)) return

  IF(solvergeom==RZgeom .OR. solvergeom==XZgeom) nguardy = 0
  IF(solvergeom==Zgeom) then
    nguardx = 0
    nguardy = 0
  END if

  ALLOCATE(grids)
  NULLIFY(grids%next,grids%prev,grids%down,grids%up,grids%bnd)
#ifdef MPIPARALLEL
  workfact = mgridrz_workfact
  nzp=nzpslave(my_index)
  ALLOCATE(grids%rho(nr+1,nz+1),grids%phi(1-nguardx:nr+1+nguardx,1-nguardz:nz+1+nguardz), &
           grids%rhop(nr+1,nzp+1),grids%phip(1-nguardx:nr+1+nguardx,1-nguardz:nzp+1+nguardz), &
           grids%loc_part(nr+1,nzp+1),grids%invvol(nr+1), &
           grids%loc_part_field_dep(nr+1,nzp+1),grids%invvol(nr+1))
#else
  ALLOCATE(grids%rho(nr+1,nz+1),grids%phi(1-nguardx:nr+1+nguardx,1-nguardz:nz+1+nguardz), &
           grids%loc_part(nr+1,nz+1),grids%invvol(nr+1), &
           grids%loc_part_field_dep(nr+1,nz+1),grids%invvol(nr+1))
#endif
WRITE(0,*) '******** Basegrid allocated.',nr,nz
  grids%id=1
  grids_nids=1
  grids%nr=nr
  grids%dr=dr
  grids%rmin=rmin
  grids%rmax=rmin+nr*dr
  grids%xmin=rmin
  grids%xmax=rmin+nr*dr
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
  IF(solvergeom==RZgeom) then
    ! computes divider by cell volumes to get density
    j = 1
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    grids%invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
    do j = 2, nr+1
      grids%invvol(j) = 1._8 / (2._8 * pi * real(j-1,8) * dr * dr * dz)
    end do
  else IF(solvergeom==XZgeom) then
    grids%invvol(:) = 1._8 / (dr * dz)
  END if

  ixlbndi = boundxy
  ixrbndi = boundxy
  izlbndi = bound0
  izrbndi = boundnz
  ixlbnd = boundxy
  ixrbnd = boundxy
  izlbnd = bound0
  izrbnd = boundnz
  IF(l2symtry .or. l4symtry .or. solvergeom==RZgeom) then
    ixlbndi=neumann
    ixlbnd =neumann
  endif
  grids%ixlbnd = ixlbndi
  grids%ixrbnd = ixrbndi
  grids%izlbnd = izlbndi
  grids%izrbnd = izrbndi

  IF(solvergeom==Zgeom) then
    grids%nlevels=nlevels
  else
    call init_bnd(grids%bnd,nr,nz,dr,dz,grids%zmin,grids%zmax)
    grids%nlevels=nlevels
  END if

!  do i = 1,grids%nlevels, 1
!    grids%bnd(i)%izlbnd=grids%izlbnd
!    grids%bnd(i)%izrbnd=grids%izrbnd
!  END do

  mglevels = basegrid%nlevels
  do i = 1, basegrid%nlevels
    mglevel = basegrid%nlevels - i
    mglevelsnx(mglevel) = basegrid%bnd(i)%nr
    mglevelsny(mglevel) = 0
#ifdef MPIPARALLEL
    mglevelsnzfull(mglevel) = basegrid%bnd(i)%nz*nslaves/basegrid%bnd(i)%nworkpproc
    mglevelsiz(mglevel) = INT(my_index/basegrid%bnd(i)%nworkpproc)*basegrid%bnd(i)%nz-1
#else
    mglevelsnzfull(mglevel) = basegrid%bnd(i)%nz
    mglevelsiz(mglevel) = 0
#endif
    mglevelsnz(mglevel) = basegrid%bnd(i)%nz
    IF(mglevel==0) then
      mglevelslx(mglevel) = 1.
      mglevelsly(mglevel) = 1.
      mglevelslz(mglevel) = 1.
    else
      mglevelslx(mglevel) = basegrid%bnd(i)%dr/basegrid%bnd(basegrid%nlevels)%dr
      mglevelsly(mglevel) = 1.
      mglevelslz(mglevel) = basegrid%bnd(i)%dz/basegrid%bnd(basegrid%nlevels)%dz
    END if
  end do
  do i = 1,basegrid%nlevels
    IF(basegrid%bnd(i)%izlbnd==dirichlet) basegrid%bnd(i)%v(:,1) = v_dirichlet
    IF(basegrid%bnd(i)%izrbnd==dirichlet) basegrid%bnd(i)%v(:,basegrid%bnd(i)%nz+1) = v_dirichlet
    IF(basegrid%ixlbnd==dirichlet) basegrid%bnd(i)%v(1,:) = v_dirichlet
    IF(basegrid%ixrbnd==dirichlet) basegrid%bnd(i)%v(basegrid%bnd(i)%nr+1,:) = v_dirichlet
  END do

  call mk_grids_ptr()

return
end subroutine init_basegrid

subroutine del_basegrid()
USE Multigrid3d
implicit none
INTEGER(ISZ) :: i

WRITE(0,*) 'enter del_basegrid'
  IF(.NOT.associated(basegrid)) then
    WRITE(0,*) 'Warning: call to del_basegrid while grid not associated'
    return
  end if

  do i = 1, basegrid%nlevels
    call del_cnds(basegrid%bnd(i))
  end do
  DEALLOCATE(basegrid%bnd)
#ifdef MPIPARALLEL
  DEALLOCATE(basegrid%rho,basegrid%phi, &
           basegrid%rhop,basegrid%phip, &
           basegrid%loc_part,basegrid%invvol, &
           basegrid%loc_part_field_dep,basegrid%invvol)
#else
  DEALLOCATE(basegrid%rho,basegrid%phi, &
           basegrid%loc_part,basegrid%invvol, &
           basegrid%loc_part_field_dep,basegrid%invvol)
#endif
  DEALLOCATE(basegrid)
WRITE(0,*) 'leave del_basegrid'

return
end subroutine del_basegrid

subroutine add_grid(mothergrid,nr,nz,dri,dzi,rmini,zmini,guard_min_r,guard_max_r,guard_min_z,guard_max_z,lshifts)
implicit none
TYPE(grdptr), pointer :: mothergrid
INTEGER(ISZ), INTENT(IN) :: nr, nz, guard_min_r, guard_max_r, guard_min_z, guard_max_z
REAL(8), INTENT(IN) :: dri,dzi,rmini,zmini
INTEGER(ISZ), INTENT(IN), OPTIONAL :: lshifts(nr+1)

TYPE(grdptr), pointer :: g  ! new grid
INTEGER(ISZ) :: i,j,l,ls,jmin,jmax,lmin,lmax,ratio_r,ratio_z,nzs
REAL(8) :: dr,dz,rmin,zmin

! adjust new grid boundaries to fall onto mother grid lines
! and recalculate mesh spacing for new grid

  jmin = 1 + NINT( (MAX(rmini,       mothergrid%rmin)-mothergrid%rmin) / mothergrid%dr)
  jmax = 1 + NINT( (MIN(rmini+nr*dri,mothergrid%rmax)-mothergrid%rmin) / mothergrid%dr)
  lmin = 1 + NINT( (MAX(zmini,       mothergrid%zmin)-mothergrid%zmin) / mothergrid%dz)
  lmax = 1 + NINT( (MIN(zmini+nz*dzi,mothergrid%zmax)-mothergrid%zmin) / mothergrid%dz)

  rmin = mothergrid%rmin + (jmin-1) * mothergrid%dr
  zmin = mothergrid%zmin + (lmin-1) * mothergrid%dz

  dr = (jmax-jmin) * mothergrid%dr / nr
  dz = (lmax-lmin) * mothergrid%dz / nz

! allocate new grid

  ALLOCATE(g)
  NULLIFY(g%next,g%prev,g%down,g%up,g%bnd)
  IF(.not. PRESENT(lshifts)) then
    nzs = 0
    ALLOCATE(g%rho(nr+1,nz+1), &
             g%phi(1-nguardx:nr+1+nguardx,1-nguardz:nz+1+nguardz), &
             g%loc_part(nr+1,nz+1), &
             g%loc_part_field_dep(nr+1,nz+1), &
             g%invvol(nr+1))
  else
    nzs = MAXVAL(lshifts)
    ALLOCATE(g%rho(nr+1,nz+nzs+1), &
             g%phi(1-nguardx:nr+1+nguardx,1-nguardz:nz+nzs+1+nguardz), &
             g%loc_part(nr+1,nz+nzs+1), &
             g%loc_part_field_dep(nr+1,nz+nzs+1), &
             g%invvol(nr+1))
  END if
  IF(n_avail_ids==0) then
    g%id=grids_nids+1
    grids_nids=grids_nids+1
  else
    g%id=avail_ids(n_avail_ids)
    n_avail_ids=n_avail_ids-1
  END if
  g%jmin=jmin
  g%jmax=jmax
  g%lmin=lmin
  g%lmax=lmax
  g%loc_part=g%id
  g%loc_part_field_dep=g%id
  g%phi=0.
  g%rho=0.
  g%nr=nr
  g%dr=dr
  g%rmin=rmin
  g%rmax=rmin+nr*dr
  g%xmin=rmin
  g%xmax=rmin+nr*dr
  g%nz=nz
  g%dz=dz
  g%zmin=zmin
  g%zminp=zmin
  g%zmax=zmin+nz*dz
  g%mgparam = basegrid%mgparam
  g%npre = basegrid%npre
  g%npost = basegrid%npost
  g%ncycles = basegrid%ncycles
  g%ncmax = basegrid%ncmax
  g%npmin = basegrid%npmin
  g%guard_min_r = guard_min_r
  g%guard_max_r = guard_max_r
  g%guard_min_z = guard_min_z
  g%guard_max_z = guard_max_z
  IF(associated(mothergrid%down)) then
    IF(associated(mothergrid%down%next)) then
      mothergrid%down%next%prev => g
      g%next => mothergrid%down%next
    END if
    g%prev => mothergrid%down
    mothergrid%down%next => g
  else
    mothergrid%down=>g
  END if
  g%up => mothergrid
  ngrids=ngrids+1
  mgridrz_ngrids = ngrids

! assign boundary types

  IF(ABS(g%xmin-mothergrid%xmin)<0.1*mothergrid%dr) then
    g%ixlbnd = mothergrid%ixlbnd
  else
    g%ixlbnd = dirichlet
  END if
  IF(ABS(g%rmax-mothergrid%rmax)<0.1*mothergrid%dr) then
    g%ixrbnd = mothergrid%ixrbnd
  else
    g%ixrbnd = dirichlet
  END if
  IF(ABS(g%zmin-mothergrid%zmin)<0.1*mothergrid%dz) then
    g%izlbnd = mothergrid%izlbnd
  else
    g%izlbnd = dirichlet
  END if
  IF(ABS(g%zmax-mothergrid%zmax)<0.1*mothergrid%dz) then
    g%izrbnd = mothergrid%izrbnd
  else
    g%izrbnd = dirichlet
  END if
  IF(PRESENT(lshifts)) then
    g%izlbnd = dirichlet
    g%izrbnd = dirichlet
  END if

! computes commodity quantities for charge deposition

  g%invdr = 1._8/dr
  g%invdz = 1._8/dz
  IF(solvergeom==RZgeom) then
    ! computes divider by cell volumes to get density
    j = 1
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    g%invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
    do j = 2, nr+1
      g%invvol(j) = 1._8 / (2._8 * pi * real(j-1,8) * dr * dr * dz)
    end do
  else IF(solvergeom==XZgeom) then
    g%invvol(:) = 1._8 / (dr * dz)
  END if


  WRITE(0,'(" Add grid ID: ",i5)') g%id
  IF(solvergeom==Zgeom) then
    g%nlevels=nlevels
  else
    IF(PRESENT(lshifts)) then
      call init_bnd(g%bnd,nr,nz,dr,dz,g%zmin,g%zmax,lshifts)
    else
      call init_bnd(g%bnd,nr,nz,dr,dz,g%zmin,g%zmax)
    END if
    g%nlevels=nlevels
  END if
  do i = 1,g%nlevels, 1
    g%bnd(i)%izlbnd=g%izlbnd
    g%bnd(i)%izrbnd=g%izrbnd
  END do
  do i = 1,g%nlevels
    IF(g%bnd(i)%izlbnd==dirichlet) g%bnd(i)%v(:,1) = v_dirichlet
    IF(g%bnd(i)%izrbnd==dirichlet) g%bnd(i)%v(:,g%bnd(i)%nz+1:) = v_dirichlet
    IF(g%ixlbnd==dirichlet)        g%bnd(i)%v(1,:) = v_dirichlet
    IF(g%ixrbnd==dirichlet)        g%bnd(i)%v(g%bnd(i)%nr+1,:) = v_dirichlet
  END do

! update grid pointers list array

  call mk_grids_ptr()

! initializes loc_part* arrays

  IF(solvergeom==Zgeom) then
    g%up%loc_part(1,lmin:lmax-1)=g%id
    g%up%loc_part_field_dep(1,lmin+guard_min_z:lmax-1-guard_max_z)=g%id
  else
    IF(.not. PRESENT(lshifts)) then
      g%up%loc_part(jmin:jmax-1,lmin:lmax-1)=g%id
      g%up%loc_part_field_dep(jmin+guard_min_r:jmax-1-guard_max_r,lmin+guard_min_z:lmax-1-guard_max_z)=g%id
    else
      ratio_r = NINT(mothergrid%dr/g%dr)
      ratio_z = NINT(mothergrid%dz/g%dz)
      do j = jmin, jmax-1
        ls = INT(lshifts(1+(j-jmin)*ratio_r)/ratio_z)-g%up%bnd(g%up%nlevels)%lshift(j)
        g%up%loc_part(j,lmin+ls:lmax-1+ls) = g%id
      end do
      do j = jmin+guard_min_r, jmax-1-guard_max_r
        ls = INT(lshifts(1+(j-jmin)*ratio_r)/ratio_z)-g%up%bnd(g%up%nlevels)%lshift(j)
        g%up%loc_part_field_dep(j,lmin+ls+guard_min_z:lmax-1+ls-guard_max_z) = g%id
      end do
    END if
  END if

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

recursive subroutine del_grid(g,next_too)
implicit none
TYPE(grdptr), pointer :: g
LOGICAL(ISZ),OPTIONAL :: next_too
INTEGER(ISZ):: i

  IF(associated(g%down)) call del_grid(g%down,.true.)
  IF(PRESENT(next_too)) then
    IF(associated(g%next)) call del_grid(g%next,.true.)
  END if

  WHERE(g%up%loc_part==g%id)
    g%up%loc_part = g%up%id
  END where
  WHERE(g%up%loc_part_field_dep==g%id)
    g%up%loc_part_field_dep = g%up%id
  END where

  level_del_grid=level_del_grid+1
  IF(associated(g%next)) then
    IF(associated(g%prev)) then
      g%next%prev => g%prev
      g%prev%next => g%next
    else
      NULLIFY(g%next%prev)
    END if
  END if

  n_avail_ids=n_avail_ids+1
  avail_ids(n_avail_ids)=g%id
  DEALLOCATE(g%rho,g%phi,g%loc_part,g%loc_part_field_dep,g%invvol)
!   nullify(g%next,g%prev,g%down,g%up)
  IF(solvergeom/=Zgeom) then
    do i = 1, g%nlevels
      call del_cnds(g%bnd(i))
      DEALLOCATE(g%bnd(i)%lshift)
    end do
    DEALLOCATE(g%bnd)
  END if
  NULLIFY(g%up%down)
  DEALLOCATE(g)

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
!  WRITE(0,*) 'call gchange ', ngrids
  call gchange("FRZmgrid",0)
!  WRITE(0,*) 'done.'
  do i = 1, ngrids
    nrg(i) = grids_ptr(i)%grid%nr
    nzg(i) = SIZE(grids_ptr(i)%grid%loc_part,2)-1
    drg(i) = grids_ptr(i)%grid%dr
    dzg(i) = grids_ptr(i)%grid%dz
  END do

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

subroutine init_bnd(bndy,nr,nz,dr,dz,zmin,zmax,lshifts)
! intializes grid quantities according to the number of multigrid levels and grid sizes nx and nz.
USE InGen3d, ONLY:l2symtry, l4symtry
USE Picglb3d, ONLY:dy
USE InMesh3d, ONLY:ny
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr, dz, zmin, zmax
TYPE(bndptr), pointer :: bndy(:)
INTEGER(ISZ),INTENT(IN), OPTIONAL :: lshifts(nr+1)

INTEGER(ISZ) :: i, nrp0, nzp0, nrc, nzc, nrc_old, nzc_old, j, jcoarse, lcoarse
REAL(8) :: drc, dzc

!  dy = dr
  ny = 0

  nrc = nr
  nzc = nz
  drc = dr
  dzc = dz
  nlevels    = 1
#ifdef MPIPARALLEL
  nworkpproc = 1
#endif
! first loop to compute the number of levels
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
    bndy(i)%izlbnd=izlbndi
    bndy(i)%izrbnd=izrbndi
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
    ALLOCATE(bndy(i)%lshift(0:bndy(i)%nr+2))
    IF(PRESENT(lshifts)) then
      bndy(i)%l_lshift = .true.
      IF(i==nlevels) then
        bndy(i)%lshift(1:bndy(i)%nr+1) = lshifts
      else
        bndy(i)%lshift(1:bndy(i)%nr+1) = 1000000
        do j = 1,bndy(i+1)%nr
          jcoarse = 1+INT(j*bndy(i+1)%dr/bndy(i)%dr)
          lcoarse = 1+INT(bndy(i+1)%lshift(j)*bndy(i+1)%dz/bndy(i)%dz)
          bndy(i)%lshift(jcoarse) = MIN(bndy(i)%lshift(jcoarse),lcoarse-1)
          bndy(i)%lshift(jcoarse+1) = MIN(bndy(i)%lshift(jcoarse+1),lcoarse-1)
          ! check if coarse grid long enough and adjust if necessary
          lcoarse = 1+INT((bndy(i+1)%nz+bndy(i+1)%lshift(jcoarse))*bndy(i+1)%dz/bndy(i)%dz)-bndy(i)%lshift(jcoarse)
          bndy(i)%nz = MAX(bndy(i)%nz,lcoarse)
          WRITE(0,*) 'j=',j,jcoarse,lcoarse,bndy(i)%lshift(jcoarse),bndy(i)%nz
        end do
        WRITE(0,*) 'level ',i,' nz = ',bndy(i)%nz
        bndy(i)%lshift(bndy(i)%nr+1) = bndy(i)%lshift(bndy(i)%nr)
      END if
      bndy(i)%lshift(0) = bndy(i)%lshift(1)
      bndy(i)%lshift(bndy(i)%nr+2) = bndy(i)%lshift(bndy(i)%nr+1)
    else
      bndy(i)%l_lshift = .false.
      bndy(i)%lshift = 0
    END if
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
bndl%cnd%nbbndred = nbnd
bndl%cnd%ncond = ncond
IF(nbnd>0) then
  ALLOCATE(bndl%cnd%phi0xm(nbnd),bndl%cnd%phi0xp(nbnd),bndl%cnd%phi0zm(nbnd),bndl%cnd%phi0zp(nbnd),bndl%cnd%dt(nbnd))
  ALLOCATE(bndl%cnd%volt0xm(nbnd),bndl%cnd%volt0xp(nbnd),bndl%cnd%volt0zm(nbnd),bndl%cnd%volt0zp(nbnd))
  ALLOCATE(bndl%cnd%condidxm(nbnd),bndl%cnd%condidxp(nbnd),bndl%cnd%condidzm(nbnd),bndl%cnd%condidzp(nbnd))
  ALLOCATE(bndl%cnd%cf0(nbnd),bndl%cnd%cfxp(nbnd),bndl%cnd%cfxm(nbnd),bndl%cnd%cfzp(nbnd),bndl%cnd%cfzm(nbnd))
  ALLOCATE(bndl%cnd%dxm(nbnd),bndl%cnd%dxp(nbnd),bndl%cnd%dzm(nbnd),bndl%cnd%dzp(nbnd))
  ALLOCATE(bndl%cnd%jj(nbnd),bndl%cnd%kk(nbnd),bndl%cnd%docalc(nbnd))
  bndl%cnd%docalc=.false.
END if
IF(ncond>0) then
  ALLOCATE(bndl%cnd%jcond(ncond),bndl%cnd%kcond(ncond),bndl%cnd%voltage(ncond),bndl%cnd%condid(ncond))
END if

end subroutine init_bnd_sublevel

subroutine del_cnds(bnd)
! initializes quantities for one grid level
implicit none
TYPE(bndptr), INTENT(IN OUT) :: bnd
TYPE(conductor_type), POINTER :: cnd,cndnext

WRITE(0,*) 'enter del_cnds'

IF(bnd%nb_conductors==0) return
cnd => bnd%first
do WHILE(associated(cnd%next))
  WRITE(0,*) '#1'
  cndnext => cnd%next
  WRITE(0,*) '#2'
  call del_cnd(cnd)
  WRITE(0,*) '#3'
  bnd%nb_conductors = bnd%nb_conductors - 1
  WRITE(0,*) '#4'
  cnd=>cndnext
  WRITE(0,*) '#5'
end do
WRITE(0,*) 'call cnd'
WRITE(0,*) associated(cnd)

call del_cnd(cnd)
bnd%nb_conductors = bnd%nb_conductors - 1
WRITE(0,*) 'leave del_cnds'

end subroutine del_cnds

subroutine del_cnd(cnd)
implicit none
TYPE(conductor_type), POINTER :: cnd

WRITE(0,*) 'enter del_cnd'
WRITE(0,*) '##1'
  IF(cnd%nbbnd>0) then
    DEALLOCATE(cnd%phi0xm,cnd%phi0xp,cnd%phi0zm,cnd%phi0zp,cnd%dt)
    DEALLOCATE(cnd%volt0xm,cnd%volt0xp,cnd%volt0zm,cnd%volt0zp)
    DEALLOCATE(cnd%condidxm,cnd%condidxp,cnd%condidzm,cnd%condidzp)
    DEALLOCATE(cnd%cf0,cnd%cfxp,cnd%cfxm,cnd%cfzp,cnd%cfzm)
    DEALLOCATE(cnd%dxm,cnd%dxp,cnd%dzm,cnd%dzp)
    DEALLOCATE(cnd%jj,cnd%kk,cnd%docalc)
  END if
WRITE(0,*) '##2'
  IF(cnd%ncond>0) then
    DEALLOCATE(cnd%jcond,cnd%kcond,cnd%voltage,cnd%condid)
  END if
WRITE(0,*) '##3'
  DEALLOCATE(cnd)
WRITE(0,*) 'leave del_cnd'

end subroutine del_cnd

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
                          ixlbnd, ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguard_any(0:nxnew+2,0:nznew+2)

INTEGER(ISZ) :: nxold, nzold, nri, nrf, nzi, nzf
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

IF(ixlbnd==dirichlet) then
  nri=2
else
  nri=1
END if
IF(solvergeom==RZgeom) nri=1
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
do jnew = nri, nrf
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

subroutine add_and_expand(unew, uold, bnd, nxnew, nznew, nxold, nzold, &
                          xminnew, xmaxnew, zminnew, zmaxnew, xminold,  xmaxold,  zminold,  zmaxold, &
                          ixlbnd, ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
TYPE(bndptr) :: bnd
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:nxnew+2,0:nznew+2), INTENT(INOUT) :: unew
REAL(8), DIMENSION(0:nxold+2,0:nzold+2), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew

INTEGER(ISZ) :: nri, nrf, nzi, nzf
INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, ii, ic
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+1), koldp(nznew+1)

IF(l_mgridrz_debug) WRITE(0,*) 'enter expand, level = ',level

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

IF(ixlbnd==dirichlet) then
  nri=2
else
  nri=1
END if
IF(solvergeom==RZgeom) nri=1
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
do jnew = nri, nrf
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do
IF(vlocs) then
  do ii = 1, bnd%nvlocs
    jnew = bnd%vlocs_j(ii)
    knew = bnd%vlocs_k(ii)
    k = kold(knew)
    kp = koldp(knew)
    delz = ddz(knew)
    odelz = oddz(knew)
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
    unew(jnew,knew) = unew(jnew,knew) + uold(j, k)  * odelx * odelz &
                                      + uold(jp,k)  * delx  * odelz &
                                      + uold(j, kp) * odelx * delz &
                                      + uold(jp,kp) * delx  * delz
  END do
else
  do knew = nzi, nzf
    k = kold(knew)
    kp = koldp(knew)
    delz = ddz(knew)
    odelz = oddz(knew)
    do jnew = 1, nrf
      IF(.NOT.bnd%v(jnew,knew)==v_vacuum) cycle
      j = jold(jnew)
      jp = joldp(jnew)
      delx = ddx(jnew)
      odelx = oddx(jnew)
      unew(jnew,knew) = unew(jnew,knew) + uold(j, k)  * odelx * odelz &
                                        + uold(jp,k)  * delx  * odelz &
                                        + uold(j, kp) * odelx * delz &
                                        + uold(jp,kp) * delx  * delz
    end do
  END do
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END subroutine add_and_expand

subroutine add_and_expand_wshift(unew, uold, snew, sold, bnd, nxnew, nznew, nzsnew, nxold, nzold, nzsold, &
                                 xminnew, zminnew, dxnew, dznew, xminold, zminold, dxold, dzold)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
TYPE(bndptr) :: bnd
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold, nzsnew, nzsold
REAL(8), DIMENSION(0:nxnew+2,0:nznew+nzsnew+2), INTENT(INOUT) :: unew
REAL(8), DIMENSION(0:nxold+2,0:nzold+nzsold+2), INTENT(IN) :: uold
INTEGER(ISZ), DIMENSION(0:nxnew+2), INTENT(IN) :: snew
INTEGER(ISZ), DIMENSION(0:nxold+2), INTENT(IN) :: sold
REAL(8), INTENT(IN) :: xminold, zminold, dxold, dzold, xminnew, zminnew, dxnew, dznew

INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, ii, ic, knews
REAL(8) :: x, z, xx, zz, invdxold, invdzold, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+nzsnew+1), oddz(nznew+nzsnew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+nzsnew+1), koldp(nznew+nzsnew+1)

IF(l_mgridrz_debug) WRITE(0,*) 'enter expand, level = ',level

invdxold = 1./dxold
invdzold = 1./dzold

do knew = 1, nznew+nzsnew+1
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
IF(vlocs) then
  do ii = 1, bnd%nvlocs
    jnew = bnd%vlocs_j(ii)
    knew = bnd%vlocs_k(ii)
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
    knews = knew+snew(jnew)
    k = kold(knews)
    kp = koldp(knews)
    delz = ddz(knews)
    odelz = oddz(knews)
    unew(jnew,knew) = unew(jnew,knew) + uold(j, k-sold(j))   * odelx * odelz &
                                      + uold(jp,k-sold(jp))  * delx  * odelz &
                                      + uold(j, kp-sold(j))  * odelx * delz &
                                      + uold(jp,kp-sold(jp)) * delx  * delz
  END do
else
  do knew = 1, nznew+1
    do jnew = 1, nxnew+1
      IF(.NOT.bnd%v(jnew,knew)==v_vacuum) cycle
      j = jold(jnew)
      jp = joldp(jnew)
      delx = ddx(jnew)
      odelx = oddx(jnew)
      knews = knew+snew(jnew)
      k = kold(knews)
      kp = koldp(knews)
      delz = ddz(knews)
      odelz = oddz(knews)
      unew(jnew,knew) = unew(jnew,knew) + uold(j, k-sold(j))   * odelx * odelz &
                                        + uold(jp,k-sold(jp))  * delx  * odelz &
                                        + uold(j, kp-sold(j))  * odelx * delz &
                                        + uold(jp,kp-sold(jp)) * delx  * delz
    end do
  END do
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END subroutine add_and_expand_wshift

function expandwguardandbnd_any(uold, bnd, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, &
                                                         xminold, xmaxold, zminold, zmaxold, &
                                                         ixlbnd, ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguardandbnd_any(0:nxnew+2,0:nznew+2)
TYPE(bndptr) :: bnd

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jnew, knew, i, j, k, jp, kp, nri, nrf, nzi, nzf
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

IF(ixlbnd==dirichlet) then
  nri=2
else
  nri=1
END if
IF(solvergeom==RZgeom) nri=1
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
do jnew = nri, nrf
  x = xminnew+(jnew-1)*dxnew
  xx = (x-xminold) * invdxold
  jold(jnew) = 1+MIN(nxold,INT(xx))
  joldp(jnew) = jold(jnew) + 1
  ddx(jnew) = xx-(jold(jnew)-1)
  oddx(jnew) = 1.-ddx(jnew)
END do

IF(vlocs) then

  do i = 1, bnd%nvlocs
    jnew = bnd%vlocs_j(i)
    knew = bnd%vlocs_k(i)
    k = kold(knew)
    kp = koldp(knew)
    delz = ddz(knew)
    odelz = oddz(knew)
    j = jold(jnew)
    jp = joldp(jnew)
    delx = ddx(jnew)
    odelx = oddx(jnew)
    expandwguardandbnd_any(jnew,knew) = uold(j, k)  * odelx * odelz &
                                      + uold(jp,k)  * delx  * odelz &
                                      + uold(j, kp) * odelx * delz &
                                      + uold(jp,kp) * delx  * delz
  END do

else

do knew = nzi, nzf
  k = kold(knew)
  kp = koldp(knew)
  delz = ddz(knew)
  odelz = oddz(knew)
  do jnew = nri, nrf
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

END if
!expandwguardandbnd_any(0,:) = 0._8
!expandwguardandbnd_any(nxnew+2,:) = 0._8
!expandwguardandbnd_any(:,0) = 0._8
!expandwguardandbnd_any(:,nznew+2) = 0._8

IF(l_mgridrz_debug) WRITE(0,*) 'exit expand, level = ',level

return
END function expandwguardandbnd_any

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
!     WRITE(0,*) 'knew',knew
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

subroutine interpolate_any_1d(unew, uold, nznew, nzold, &
                              zminnew, zmaxnew, &
                              zminold,  zmaxold, &
                              izlbnd, izrbnd, &
                              bnd_only, quad)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nznew, nzold, izlbnd, izrbnd
REAL(8), DIMENSION(0:), INTENT(IN) :: uold
REAL(8), DIMENSION(0:), INTENT(IN OUT) :: unew
REAL(8), INTENT(IN) :: zminold, zmaxold, zminnew, zmaxnew
LOGICAL(ISZ), INTENT(IN) :: bnd_only, quad

INTEGER(ISZ) :: knew, k, kp, km, kpp
REAL(8) :: z, zz, invdzold, dznew, delz, odelz, &
           s1z, s2z, s3z, s4z
REAL(8) :: ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: kold(nznew+1)

dznew = (zmaxnew-zminnew) / nznew
invdzold = REAL(nzold,8)/(zmaxold-zminold)

do knew = 1, nznew+1
  z = zminnew+(knew-1)*dznew
  zz = (z-zminold) * invdzold
  kold(knew) = MAX(1,MIN(nzold,1+INT(zz)))
  ddz(knew) = zz-(kold(knew)-1)
  oddz(knew) = 1.-ddz(knew)
END do
IF(.not.quad) then
  IF(bnd_only) then
    do knew = 1, nznew+1, nznew
     IF((knew==1.and.izlbnd==dirichlet).OR.(knew==nznew+1.and.izrbnd==dirichlet)) then
      k = kold(knew)
      kp = k+1
      delz = ddz(knew)
      odelz = oddz(knew)
      unew(knew) = uold(k)  * odelz &
                 + uold(kp) *  delz
     END if
    end do
  else
    do knew = 1, nznew+1
     k = kold(knew)
     kp = k+1
     delz = ddz(knew)
     odelz = oddz(knew)
     unew(knew) = uold(k)  * odelz &
                + uold(kp) * delz
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
      unew(knew) = uold(km ) * s1z &
                 + uold(k  ) * s2z &
                 + uold(kp ) * s3z &
                 + uold(kpp) * s4z
     END if
    end do
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
      unew(knew) = uold(km ) * s1z &
                 + uold(k  ) * s2z &
                 + uold(kp ) * s3z &
                 + uold(kpp) * s4z
    END do
  END if
END if

return
END subroutine interpolate_any_1d

subroutine interpolate_any_wshift(unew, uold, snew, sold, bnd, nxnew, nznew, nzsnew, nxold, nzold, nzsold, &
                                 xminnew, zminnew, dxnew, dznew, xminold, zminold, dxold, dzold)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
TYPE(bndptr) :: bnd
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold, nzsnew, nzsold
REAL(8), DIMENSION(0:nxnew+2,0:nznew+nzsnew+2), INTENT(INOUT) :: unew
REAL(8), DIMENSION(0:nxold+2,0:nzold+nzsold+2), INTENT(IN) :: uold
INTEGER(ISZ), DIMENSION(0:nxnew+2), INTENT(IN) :: snew
INTEGER(ISZ), DIMENSION(0:nxold+2), INTENT(IN) :: sold
REAL(8), INTENT(IN) :: xminold, zminold, dxold, dzold, xminnew, zminnew, dxnew, dznew

INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, ii, ic, knews
REAL(8) :: x, z, xx, zz, invdxold, invdzold, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+nzsnew+1), oddz(nznew+nzsnew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+nzsnew+1), koldp(nznew+nzsnew+1)

IF(l_mgridrz_debug) WRITE(0,*) 'enter interpolate_any_wshift, level = ',level

  invdxold = 1./dxold
  invdzold = 1./dzold

  do knew = 1, nznew+nzsnew+1
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
    do jnew = 1, nxnew+1
      j = jold(jnew)
      jp = joldp(jnew)
      delx = ddx(jnew)
      odelx = oddx(jnew)
      knews = knew+snew(jnew)
      k = kold(knews)
      kp = koldp(knews)
      delz = ddz(knews)
      odelz = oddz(knews)
      unew(jnew,knew) = unew(jnew,knew) + uold(j, k-sold(j))   * odelx * odelz &
                                        + uold(jp,k-sold(jp))  * delx  * odelz &
                                        + uold(j, kp-sold(j))  * odelx * delz &
                                        + uold(jp,kp-sold(jp)) * delx  * delz
    end do
  END do

IF(l_mgridrz_debug) WRITE(0,*) 'exit interpolate_any_wshift, level = ',level

return
END subroutine interpolate_any_wshift

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

!pgi$r nobounds
subroutine subrestrict(unew, uold, nxnew, nznew, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), DIMENSION(1:nxnew+1,1:nznew+1), INTENT(OUT) :: unew
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: rap(1:nxnew+1,1:nznew+1)

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

unew = 0._8
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
    unew(j,k)   = unew(j,k)   + uold(jold,kold) * odelx * odelz
    unew(jp,k)  = unew(jp,k)  + uold(jold,kold) * delx  * odelz
    unew(j,kp)  = unew(j,kp)  + uold(jold,kold) * odelx * delz
    unew(jp,kp) = unew(jp,kp) + uold(jold,kold) * delx  * delz
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
    IF(rap(j,k)/=0._8) unew(j,k)   = unew(j,k)   / rap(j,k)
  end do
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END subroutine subrestrict

subroutine subrestrict_wshift(unew, uold, snew, sold, nxnew, nznew, nzsnew, nxold, nzold, nzsold, &
                              xminold, zminold, dxold, dzold, xminnew, zminnew, dxnew, dznew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold, nzsnew, nzsold
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), DIMENSION(1:nxnew+1,1:nznew+nzsnew+1), INTENT(OUT) :: unew
REAL(8), DIMENSION(0:nxnew+2), INTENT(IN) :: snew
REAL(8), DIMENSION(0:nxold+2), INTENT(IN) :: sold
REAL(8), INTENT(IN) :: xminold, zminold, dxold, dzold, xminnew, zminnew, dxnew, dznew
REAL(8) :: rap(1:nxnew+1,1:nznew+1)

INTEGER(ISZ) :: jold, kold, j, k, jp, kp, kolds, ksj, ksjp, kpsj, kpsjp
REAL(8) :: invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold, ixrbnd, izlbnd, izrbnd)
!  return
!END if

IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict_wshift, level = ',level

ALLOCATE(ddx(nxold+1), ddz(nzold+nzsold+1), oddx(nxold+1), oddz(nzold+nzsold+1), &
         jnew(nxold+1), knew(nzold+nzsold+1), jnewp(nxold+1), knewp(nzold+nzsold+1))

invdxnew = 1./dxnew
invdznew = 1./dznew

unew = 0._8
rap = 0._8

do kold = 1, nzold+nzsold+1
  z = zminold + (kold-1)*dzold
  knew(kold) = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  ddz(kold) = (z-zminnew) * invdznew-real(knew(kold)-1)
  knewp(kold) = knew(kold)+1
  oddz(kold) = 1._8-ddz(kold)
END do

#ifdef MPIPARALLEL
!    IF(izlbnd<0) then
!      ddz(1) = 0.5_8*ddz(1)
!      oddz(1) = 0.5_8*oddz(1)
!    END if
!   IF(izrbnd<0) then
!     ddz(nzold+1) = 0.5_8*ddz(nzold+1)
!      oddz(nzold+1) = 0.5_8*oddz(nzold+1)
!    END if
#endif

do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = 1._8-ddx(jold)
END do

do kold = 1, nzold+1
  do jold = 1, nxold+1
    j = jnew(jold)
    jp = jnewp(jold)
    delx  =  ddx(jold)
    odelx = oddx(jold)
    kolds = kold+sold(jold)
    k = knew(kolds)
    kp = knewp(kolds)
    delz  =  ddz(kolds)
    odelz = oddz(kolds)
    ksj   = k  - snew(j)
    ksjp  = k  - snew(jp)
    kpsj  = kp - snew(j)
    kpsjp = kp - snew(jp)
    unew(j, ksj)   = unew(j, ksj)   + uold(jold,kold) * odelx * odelz
    unew(jp,ksjp)  = unew(jp,ksjp)  + uold(jold,kold) * delx  * odelz
    unew(j, kpsj)  = unew(j, kpsj)  + uold(jold,kold) * odelx * delz
    unew(jp,kpsjp) = unew(jp,kpsjp) + uold(jold,kold) * delx  * delz
    rap(j, ksj)   = rap(j, ksj)   + odelx * odelz
    rap(jp,ksjp)  = rap(jp,ksjp)  + delx  * odelz
    rap(j, kpsj)  = rap(j, kpsj)  + odelx * delz
    rap(jp,kpsjp) = rap(jp,kpsjp) + delx  * delz
  end do
end do

#ifdef MPIPARALLEL
!    IF(izlbnd<0) then
!      rap(:,1) = 2.*rap(:,1)
!    END if
!    IF(izrbnd<0) then
!      rap(:,nznew+1) = 2.*rap(:,nznew+1)
!    END if
#endif
do k = 1, nznew+1
  do j = 1, nxnew+1
    IF(rap(j,k)/=0._8) unew(j,k)   = unew(j,k)   / rap(j,k)
  end do
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict_wshift, level = ',level

return
END subroutine subrestrict_wshift

subroutine restrictlist(unew, uold, rhs, bnd, nxnew, nznew, nxold, nzold, voltfact, dr, dz, &
                        xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold
TYPE(bndptr) :: bnd
REAL(8), DIMENSION(0:nxold+2,0:nzold+2), INTENT(IN) :: uold
REAL(8), DIMENSION(1:nxold,1:nzold), INTENT(IN) :: rhs
REAL(8), DIMENSION(1:nxnew+1,1:nznew+1), INTENT(OUT) :: unew
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew, voltfact, dr, dz
REAL(8) :: rapp

INTEGER(ISZ) :: nlocs
INTEGER(ISZ) :: i, ic, jold, kold, j, k, jp, kp
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz, res
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp, jlocs, klocs

IF(l_mgridrz_debug) WRITE(0,*) 'enter restrict, level = ',level

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1))

invdxnew = nxnew / (xmaxnew-xminnew)
invdznew = nznew / (zmaxnew-zminnew)
dxold = (xmaxold-xminold) / nxold
dzold = (zmaxold-zminold) / nzold

unew = 0._8

nlocs = bnd%nvlocs
do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  nlocs = nlocs + bnd%cnd%nbbnd
END do

ALLOCATE(res(nlocs),jlocs(nlocs),klocs(nlocs))

call residbndrzwguard_list(res(1),jlocs(1),klocs(1),nlocs,f=uold(0,0),rhs=rhs(1,1), &
                           bnd=bnd,nr=nxold,nz=nzold,dr=dr,dz=dz,voltfact=voltfact)

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

rapp = REAL(nxnew*nznew)/REAL(nxold*nzold)
do i = 1, nlocs
  jold = jlocs(i)
  kold = klocs(i)
  k = knew(kold)
  kp = knewp(kold)
  delz  =  ddz(kold)
  odelz = oddz(kold)
  j = jnew(jold)
  jp = jnewp(jold)
  delx  =  ddx(jold)
  odelx = oddx(jold)
  unew(j,k)   = unew(j,k)   + rapp * res(i) * odelx * odelz
  unew(jp,k)  = unew(jp,k)  + rapp * res(i) * delx  * odelz
  unew(j,kp)  = unew(j,kp)  + rapp * res(i) * odelx * delz
  unew(jp,kp) = unew(jp,kp) + rapp * res(i) * delx  * delz
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp, res, jlocs, klocs)

IF(l_mgridrz_debug) WRITE(0,*) 'exit restrict, level = ',level

return
END subroutine restrictlist

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

subroutine deposit_z(unew, uold, invvolnew, invvolold, zminold, zmaxold, zminnew, zmaxnew)
! deposit rho from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
REAL(8), DIMENSION(1:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: invvolnew, invvolold
REAL(8), INTENT(IN) :: zminold, zmaxold, zminnew, zmaxnew

INTEGER(ISZ) :: nzold, nznew
INTEGER(ISZ) :: kold, k, kp
REAL(8) :: dzold, invdznew, z, delz, odelz, volold

IF(l_mgridrz_debug) WRITE(0,*) 'enter deposit'

nzold = SIZE(uold,1) - 1
nznew = SIZE(unew,1) - 1

invdznew = nznew / (zmaxnew-zminnew)
dzold = (zmaxold-zminold) / nzold

volold = 1._8/invvolold

do kold = 1, nzold+1
  z = zminold + (kold-1)*dzold
  k = MIN(1 + INT((z-zminnew) * invdznew), nznew)
  delz = (z-zminnew) * invdznew-real(k-1)
  kp = k+1
  odelz = 1._8-delz
  unew(k)   = unew(k)   + uold(kold) * odelz * invvolnew * volold
  unew(kp)  = unew(kp)  + uold(kold) * delz  * invvolnew * volold
end do

IF(l_mgridrz_debug) WRITE(0,*) 'exit deposit'

return
END subroutine deposit_z

subroutine deposit_wshift(unew, uold, snew, sold, invvolnew, invvolold, xminold, zminold, &
                          dxold, dzold, xminnew, zminnew, dxnew, dznew)
! deposit rho from one grid to a coarser one. Each dimension may have any number of cells.
! to be finished
implicit none
REAL(8), DIMENSION(1:,1:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), DIMENSION(1:), INTENT(IN) :: invvolnew, invvolold
INTEGER(ISZ), DIMENSION(0:), INTENT(IN) :: sold, snew
REAL(8), INTENT(IN) :: xminold, zminold, dxold, dzold, xminnew, zminnew, dxnew, dznew

INTEGER(ISZ) :: nxold, nzold, nxnew, nznew, nzsold, nzsnew
INTEGER(ISZ) :: jold, kold, j, k, jp, kp
REAL(8) :: invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz, volold
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

IF(l_mgridrz_debug) WRITE(0,*) 'enter deposit'

nxold = SIZE(uold,1) - 1
nzold = SIZE(uold,2) - 1
nxnew = SIZE(unew,1) - 1
nznew = SIZE(unew,2) - 1
nzsnew = MAXVAL(snew)
nzsold = MAXVAL(sold)

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1), volold(nxold+1))

invdxnew = 1./dxnew
invdznew = 1./dznew

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
END subroutine deposit_wshift

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

subroutine relaxbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,nc,voltfact,mgparam, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN OUT) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
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
cfrhs = dt*inveps0
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
               + rhs(j,l)*inveps0)
      else
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*bnd%cnd%dt(ii)*( &
                 bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
      END if
    ENDDO
  END do
  IF(vlocs) then
    IF(redblack==2) THEN !red
      iil=1
      iiu=bnd%nvlocsred
    else !black
      iil=bnd%nvlocsred+1
      iiu=bnd%nvlocs
    ENDif
    do ii = iil, iiu
      j = bnd%vlocs_j(ii)
      l = bnd%vlocs_k(ii)
      IF(j==1) then! origin
        f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
               + 4._8*dt0*f(j+1,l)/dr**2   &
               + (dt0/dz**2)*(f(j,l+1)+f(j,l-1)) &
               + dt0*rhs(j,l)*inveps0
      else
        f(j,l) = cf0 * f(j,l) &
               + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
               + cfz*(f(j,l+1)+f(j,l-1)) &
               + cfrhs*rhs(j,l)

      end if
    end do
  else
    do l = nzi, nzf+1
      IF(jsw==2) then! origin
        j = 1
        IF(bnd%v(j,l)==v_vacuum) &
        f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                   + 4._8*dt0*f(j+1,l)/dr**2   &
                                   + (dt0/dz**2)*(f(j,l+1)+f(j,l-1)) &
                                   + dt0*rhs(j,l)*inveps0
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
  END if

#ifdef MPIPARALLEL
  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
#endif

END do !redblack=1, 2

call updateguardcellsrz(f=f, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbndrzwguard

subroutine relaxbndrzwguard_jump(f,rhs,maxjump,curjump,bnd,nr,nz,dr,dz,nc,voltfact,mgparam, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc, ixrbnd, izlbnd, izrbnd, curjump
REAL(8), INTENT(IN OUT) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
INTEGER(ISZ), INTENT(IN) :: maxjump(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact, mgparam
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic, nrf, nzi, nzf, jump
REAL(8) :: dt, dt0, dtj,drj,dzj
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs

IF(l_mgridrz_debug) WRITE(0,*) 'enter relax, level = ',level

WRITE(0,*) curjump,nc

! define CFL
dt = mgparam/(2._8/dr**2+2._8/dz**2)
dt0 = mgparam/(4._8/dr**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cf0 = 1._8-2._8*dt/dr**2-2._8*cfz
cfrhs = dt*inveps0
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
               + rhs(j,l)*inveps0)
      else
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*bnd%cnd%dt(ii)*( &
                 bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               + bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
      END if
    ENDDO
  END do
  IF(vlocs) then
    IF(redblack==2) THEN !red
      iil=1
      iiu=bnd%nvlocsred
    else !black
      iil=bnd%nvlocsred+1
      iiu=bnd%nvlocs
    ENDif
    do ii = iil, iiu
      j = bnd%vlocs_j(ii)
      l = bnd%vlocs_k(ii)
      jump = MIN(curjump,maxjump(j,l))
      IF(j==1) then! origin
        f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
               + 4._8*dt0*f(j+jump,l)/dr**2   &
               + (dt0/dz**2)*(f(j,l+jump)+f(j,l-jump)) &
               + dt0*rhs(j,l)*inveps0
      else
        f(j,l) = cf0 * f(j,l) &
               + cfrp(j)*f(j+jump,l)+cfrm(j)*f(j-jump,l)   &
               + cfz*(f(j,l+jump)+f(j,l-jump)) &
               + cfrhs*rhs(j,l)

      end if
    end do
  else
    do l = nzi, nzf+1
      IF(jsw==2) then! origin
        j = 1
        jump = MIN(curjump,maxjump(j,l))
        IF(bnd%v(j,l)==v_vacuum) &
        f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                   + 4._8*dt0*f(j+jump,l)/dr**2   &
                                   + (dt0/dz**2)*(f(j,l+jump)+f(j,l-jump)) &
                                   + dt0*rhs(j,l)*inveps0
!                                   + (1./real(jump**2,8))*dt0*rhs(j,l)*inveps0
      END if
      do j = jsw+1, nrf+1, 2
        jump = MIN(curjump,maxjump(j,l))
        dt = mgparam/(2._8/dr**2+2._8/dz**2)
        IF(bnd%v(j,l)==v_vacuum) &
          f(j,l) = cf0 * f(j,l) &
                 + dt * (1._8+(REAL(jump,8)-0.5_8)/REAL(j-1,8)) / dr**2 * f(j+jump,l) &
                 + dt * (1._8-(REAL(jump,8)-0.5_8)/REAL(j-1,8)) / dr**2 * f(j+jump,l) &
                                   + cfz*(f(j,l+jump)+f(j,l-jump)) &
                                   + cfrhs*rhs(j,l)
!                                   + (1./real(jump**2,8))*cfrhs*rhs(j,l)

      end do
      jsw = 3-jsw
    end do
    lsw = 3-lsw
  END if

#ifdef MPIPARALLEL
  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
#endif

END do !redblack=1, 2

call updateguardcellsrz(f=f, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbndrzwguard_jump

!pgi$r nobounds
subroutine relaxbndxzwguard(f,rhs,bnd,nx,nz,dx,dz,nc,voltfact,mgparam, ixlbnd, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nx, nz, nc, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(1:,1:)!rhs(nr+1,nz+1)
REAL(8), INTENT(IN) :: dx, dz, voltfact, mgparam
TYPE(bndptr), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, redblack, iil, iiu, ic, nxi, nxf, nzi, nzf
REAL(8) :: dt, cf0, cfx, cfz, cfrhs

IF(l_mgridrz_debug) WRITE(0,*) 'enter relax, level = ',level

! define CFL
dt = mgparam/(2._8/dx**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cfx = dt / dx**2
cf0 = 1._8-2._8*(cfx+cfz)
cfrhs = dt*inveps0

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
        IF(bnd%cnd%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*bnd%cnd%dt(ii)*( &
                 bnd%cnd%cf0(ii)*f(j,l) &
               + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
               + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
               + voltfact*(bnd%cnd%phi0xm(ii)+bnd%cnd%phi0xp(ii) &
               +           bnd%cnd%phi0zm(ii)+bnd%cnd%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
    ENDDO
  END do
  do l = nzi, nzf+1
    do j = nxi+jsw-1, nxf+1, 2
      IF(bnd%v(j,l)==v_vacuum) &
        f(j,l) = cf0 * f(j,l) &
               + cfx*(f(j+1,l)+f(j-1,l))   &
               + cfz*(f(j,l+1)+f(j,l-1)) &
               + cfrhs*rhs(j,l)

    end do
    jsw = 3-jsw
  end do
  lsw = 3-lsw

#ifdef MPIPARALLEL
  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
#endif

END do !redblack=1, 2

call updateguardcellsrz(f=f, ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !i=1, nc

IF(l_mgridrz_debug) WRITE(0,*) 'exit relax, level = ',level

return
END subroutine relaxbndxzwguard

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
      call mpi_isend(rho(1,1),nr+1,mpi_double_precision,p_down,1,mpi_comm_world,mpi_req(ir),mpierror)
    end if
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_isend(rho(1,nz+1),nr+1,mpi_double_precision,p_up,1,mpi_comm_world,mpi_req(ir),mpierror)
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
REAL(8), INTENT(IN) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dr, dz,voltfact
REAL(8), DIMENSION(nr+1,nz+1) :: residbndrzwguard
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

IF(vlocs) then
  do ii = 1, bnd%nvlocs
    j = bnd%vlocs_j(ii)
    l = bnd%vlocs_k(ii)
    IF(j==1) then! origin
      residbndrzwguard(j,l) = (cf0-2._8/dr**2) * f(j,l) + 4._8*f(j+1,l)/dr**2   &
                            + cfz*(f(j,l+1)+f(j,l-1)) &
                            + rhs(j,l)*inveps0
    else
      residbndrzwguard(j,l) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                            + cfz*(f(j,l+1)+f(j,l-1)) &
                            + rhs(j,l)*inveps0
    END if
  enddo
else
 do l = nzi, nzf+1
  j = 1
  IF(bnd%v(j,l)==v_vacuum) &
  residbndrzwguard(j,l) = (cf0-2._8/dr**2) * f(j,l) + 4._8*f(j+1,l)/dr**2   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 + rhs(j,l)*inveps0

  do j = 2, nrf+1
     IF(bnd%v(j,l)==v_vacuum) &
       residbndrzwguard(j,l) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                                      + cfz*(f(j,l+1)+f(j,l-1)) &
                                      + rhs(j,l)*inveps0
  end do
 end do
END if

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
                        + rhs(j,l)*inveps0
    else
      IF(bnd%v(j,l)==v_bnd.and.bnd%cnd%docalc(ii)) &
      residbndrzwguard(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
                        + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
                        + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%phi0xp(ii)+bnd%cnd%phi0xm(ii) &
                        + bnd%cnd%phi0zp(ii)+bnd%cnd%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
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

subroutine residbndrzwguard_list(res,jlocs,klocs,nvlocs,f,rhs,bnd,nr,nz,dr,dz,voltfact)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nvlocs
REAL(8), INTENT(IN) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
REAL(8), INTENT(IN OUT) :: res(nvlocs)
INTEGER(ISZ), INTENT(IN OUT), dimension(nvlocs) :: jlocs, klocs
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dr, dz, voltfact

INTEGER(ISZ) :: i, j, l, ii, ic, nrf, nzi, nzf
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz

IF(l_mgridrz_debug) WRITE(0,*) 'enter resid, level = ',level

cfz = 1._8 / dz**2
cf0 = -2._8/dr**2-2._8*cfz
do j = 2, nr+1
  cfrp(j) = (1._8+0.5_8/REAL(j-1,8)) / dr**2
  cfrm(j) = (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

  do ii = 1, bnd%nvlocs
    j = bnd%vlocs_j(ii)
    l = bnd%vlocs_k(ii)
    jlocs(ii) = j
    klocs(ii) = l
    IF(j==1) then! origin
      res(ii) = (cf0-2._8/dr**2) * f(j,l) + 4._8*f(j+1,l)/dr**2   &
              + cfz*(f(j,l+1)+f(j,l-1)) &
              + rhs(j,l)*inveps0
    else
      res(ii) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
              + cfz*(f(j,l+1)+f(j,l-1)) &
              + rhs(j,l)*inveps0
    END if
  enddo
  i = bnd%nvlocs

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    bnd%cnd => bnd%first
  else
    bnd%cnd => bnd%cnd%next
  END if
  do ii = 1, bnd%cnd%nbbnd
    i = i+1
    j = bnd%cnd%jj(ii)
    l = bnd%cnd%kk(ii)
    jlocs(i) = j
    klocs(i) = l
    IF(bnd%cnd%docalc(ii)) then
    IF(j==1) then
      res(i) = bnd%cnd%cf0(ii)*f(j,l) &
                        + bnd%cnd%cfxp(ii)*f(j+1,l) &
                        + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%phi0xp(ii) &
                        + bnd%cnd%phi0zp(ii)+bnd%cnd%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
    else
      res(i) = bnd%cnd%cf0(ii)*f(j,l) &
                        + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
                        + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
                        + voltfact*(bnd%cnd%phi0xp(ii)+bnd%cnd%phi0xm(ii) &
                        + bnd%cnd%phi0zp(ii)+bnd%cnd%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
    END if
    else
      res(i) = 0.
    END if
  ENDDO
END do

IF(l_mgridrz_debug) WRITE(0,*) 'exit resid, level = ',level

return
end subroutine residbndrzwguard_list

function residbndxzwguard(f,rhs,bnd,nx,nz,dx,dz,voltfact,l_zerolastz, ixlbnd, ixrbnd, izlbnd, izrbnd)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nx, nz, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN) :: f(0:,0:)!f(0:nx+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(:,:)!rhs(nx+1,nz+1)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: dx, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2) :: residbndxzwguard
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, l, ii, ic, nxi, nxf, nzi, nzf
REAL(8) :: cf0, cfx, cfz

IF(l_mgridrz_debug) WRITE(0,*) 'enter resid, level = ',level

cfx = 1._8 / dx**2
cfz = 1._8 / dz**2
cf0 = -2._8*(cfx+cfz)

residbndxzwguard = 0._8

!do ic = 1, bnd%nb_conductors
!  IF(ic==1) then
!    bnd%cnd => bnd%first
!  else
!    bnd%cnd => bnd%cnd%next
!  END if
!  do i = 1, bnd%cnd%ncond
!    residbndxzwguard(bnd%cnd%jcond(i),bnd%cnd%kcond(i)) = 0._8
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
  do j = nxi, nxf+1
     IF(bnd%v(j,l)==v_vacuum) &
       residbndxzwguard(j,l) = cf0 * f(j,l) + cfx*(f(j+1,l)+f(j-1,l))   &
                                            + cfz*(f(j,l+1)+f(j,l-1)) &
                                            + rhs(j,l)*inveps0
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
    IF(bnd%v(j,l)==v_bnd.and.bnd%cnd%docalc(ii)) &
    residbndxzwguard(j,l) = bnd%cnd%cf0(ii)*f(j,l) &
                          + bnd%cnd%cfxp(ii)*f(j+1,l)+bnd%cnd%cfxm(ii)*f(j-1,l) &
                          + bnd%cnd%cfzp(ii)*f(j,l+1)+bnd%cnd%cfzm(ii)*f(j,l-1) &
                          + voltfact*(bnd%cnd%phi0xp(ii)+bnd%cnd%phi0xm(ii) &
                          + bnd%cnd%phi0zp(ii)+bnd%cnd%phi0zm(ii)) &
                          + rhs(j,l)*inveps0
  ENDDO
END do

IF(l_zerolastz) residbndxzwguard(:,nz+1) = 0._8

IF(ixlbnd==dirichlet) then
  residbndxzwguard(1,:) = 0.
END if
IF(ixrbnd==dirichlet) then
  residbndxzwguard(nx+1,:) = 0.
END if
IF(izlbnd==dirichlet) then
  residbndxzwguard(:,1) = 0.
END if
IF(izrbnd==dirichlet) then
  residbndxzwguard(:,nz+1) = 0.
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit resid, level = ',level

return
end function residbndxzwguard

!pgi$r nobounds
RECURSIVE subroutine mgbndrzwguard(j, u, rhs, bnd, nr, nz, dr, dz, npre, npost, ncycle, sub, relax_only, npmin, mgparam)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: j, nr, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(0:nr+2,0:nz+2), INTENT(IN OUT) :: u
REAL(8), DIMENSION(1:nr+1,1:nz+1), INTENT(IN) :: rhs
REAL(8) :: dr, dz, mgparam
TYPE(bndptr) :: bnd(:)
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(:,:), allocatable :: res, v
INTEGER(ISZ) :: i,jj,ll
INTEGER :: nrnext, nznext, nzresmin, nzresmax, nzres
REAL(8) :: drnext, dznext, voltf

REAL(8) :: error1, error2

IF(.not.sub) inveps0 = 1./eps0

level = j
IF(l_mgridrz_debug) WRITE(0,*) 'enter mg, level = ',level

IF(sub) then
  voltf = 0._8
else
  voltf = 1._8
END if

IF(j<=npmin .or. relax_only) then
    t_before = wtime()
  call apply_voltagewguard(u,bnd(j),voltf)
    t_apply_voltage = t_apply_voltage + wtime()-t_before
    t_before = wtime()
  call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
    t_updateguard = t_updateguard + wtime()-t_before
    t_before = wtime()
  IF(solvergeom==RZgeom) then
    call relaxbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
  else IF(solvergeom==XZgeom) then
    call relaxbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nr,nz=nz,dx=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
  END if
    t_relax = t_relax + wtime()-t_before
else
  nrnext = bnd(j-1)%nr
  nznext = bnd(j-1)%nz
  drnext = bnd(j-1)%dr
  dznext = bnd(j-1)%dz
    t_before = wtime()
  ALLOCATE(res(nrnext+1,nznext+1),v(0:nrnext+2,0:nznext+2))
    t_allocate = t_allocate + wtime()-t_before
    t_before = wtime()
  call apply_voltagewguard(u,bnd(j),voltf)
    t_apply_voltage = t_apply_voltage + wtime()-t_before
    t_before = wtime()
  call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
    t_updateguard = t_updateguard + wtime()-t_before
    t_before = wtime()
  IF(solvergeom==RZgeom) then
    call relaxbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
  else IF(solvergeom==XZgeom) then
    call relaxbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nr,nz=nz,dx=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
  END if
    t_relax = t_relax + wtime()-t_before
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
    t_before = wtime()
  IF(solvergeom==RZgeom) then
    IF(vlocs) then
      call restrictlist(res(1,nzresmin), u, rhs, bnd(j), nrnext, nzres, nr, nz, voltf,dr,dz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
    else
!      res(:,nzresmin:nzresmax) = restrict( &
!                              residbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
!                              ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd), &
!                              nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
    call subrestrict(res(1,nzresmin), &
                             residbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd(j),nr=nr,nz=nz, &
                             dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
                             ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd), &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
    END if
  else IF(solvergeom==XZgeom) then
    res(:,nzresmin:nzresmax) = restrict( &
                             residbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nr,nz=nz,dx=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
                             ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd), &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
  t_restrict = t_restrict + wtime()-t_before
!  END if
  t_before = wtime()
  call apply_voltage(res,bnd(j-1),0._8)
  t_apply_voltage = t_apply_voltage + wtime()-t_before
#ifdef MPIPARALLEL
  IF(bnd(level-1)%l_merged) call merge_work(res,level,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  call exchange_resbndz(rho=res,izlbnd=bnd(j-1)%izlbnd,izrbnd=bnd(j-1)%izrbnd)
#endif
  v = 0.0_8
  IF(.not.sub) inveps0 = 1.
  do i = 1, ncycle  !(1=V cycles, 2=W cycle)
!    call mgbndrzwguard(j=j-1, u=v(0,0), rhs=res(1,1), bnd=bnd(1:j-1), nr=nrnext, nz=nznext, dr=drnext, dz=dznext, npre=npre, npost=npost, &
    call mgbndrzwguard(j=j-1, u=v(0,0), rhs=res(1,1), bnd=bnd, nr=nrnext, nz=nznext, dr=drnext, dz=dznext, npre=npre, npost=npost, &
                   ncycle=ncycle, sub=.TRUE., relax_only=.FALSE., npmin=npmin, mgparam=mgparam)
    level = j
  end do
  IF(.not.sub) inveps0 = 1./eps0
  t_before = wtime()
  call apply_voltagewguard(v,bnd(j-1),0._8)
  t_apply_voltage = t_apply_voltage + wtime()-t_before
  t_before = wtime()
  call updateguardcellsrz(f=v,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j-1)%izlbnd,izrbnd=bnd(j-1)%izrbnd)
  t_updateguard = t_updateguard + wtime()-t_before
!  IF(bnd(j)%l_powerof2) then
!    u = u + expandwguard(v(:,nzresmin-1:nzresmax+1))
!  else
  t_before = wtime()
! IF(restrictwbnd) then
!    u = u + expandwguardandbnd_any(v(:,nzresmin-1:nzresmax+1),bnd(j),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
!             ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
!  else
!    u = u + expandwguard_any(v(:,nzresmin-1:nzresmax+1),nr,nz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
!                                ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
!  END if
!  END if
  call add_and_expand(u(0,0),v(0,nzresmin-1),bnd(j),nr,nz,nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
                                ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
  t_expand = t_expand + wtime()-t_before
  t_before = wtime()
  call apply_voltagewguard(u,bnd(j),voltf)
  t_apply_voltage = t_apply_voltage + wtime()-t_before
  t_before = wtime()
  call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd(j)%izlbnd,izrbnd=bnd(j)%izrbnd)
  t_updateguard = t_updateguard + wtime()-t_before
#ifdef MPIPARALLEL
  call exchange_fbndz(u,bnd(j)%izlbnd,bnd(j)%izrbnd)
#endif
  t_before = wtime()
  IF(solvergeom==RZgeom) then
    call relaxbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd(j),nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
  else IF(solvergeom==XZgeom) then
    call relaxbndxzwguard(f=u,rhs=rhs,bnd=bnd(j),nx=nr,nz=nz,dx=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd(j)%izlbnd, izrbnd=bnd(j)%izrbnd)
  END if
  t_relax = t_relax + wtime()-t_before
  t_before = wtime()
  DEALLOCATE(res,v)
  t_allocate = t_allocate + wtime()-t_before
END if

IF(l_mgridrz_debug) WRITE(0,*) 'exit mg, level = ',level

return
end subroutine mgbndrzwguard

subroutine mgbndrzwguard_jump(jmax, u, rhs, maxjump, bnd, nr, nz, dr, dz, npre, npost, ncycle, sub, relax_only, npmin, mgparam)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: jmax, nr, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(0:nr+2,0:nz+2), INTENT(IN OUT) :: u
REAL(8), DIMENSION(1:nr+1,1:nz+1), INTENT(IN) :: rhs
INTEGER(ISZ), DIMENSION(1:nr+1,1:nz+1), INTENT(IN) :: maxjump
REAL(8) :: dr, dz, mgparam
TYPE(bndptr) :: bnd
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(1:nr+1,1:nz+1) :: res
REAL(8), DIMENSION(0:nr+2,0:nz+2) :: v
INTEGER(ISZ) :: i,jj,ll,j
REAL(8) :: voltf

level = nlevels

IF(.not.sub) inveps0 = 1./eps0

IF(l_mgridrz_debug) WRITE(0,*) 'enter mg_jump'

voltf = 1._8
inveps0 = 1./eps0
res = residbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd,nr=nr,nz=nz, &
                       dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
                       ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
v = 0.
voltf = 0._8
inveps0 = 1.
do j = jmax, 1, -1
  t_before = wtime()
  call apply_voltagewguard(v,bnd,voltf)
  t_apply_voltage = t_apply_voltage + wtime()-t_before
  t_before = wtime()
  call updateguardcellsrz(f=v,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
  t_updateguard = t_updateguard + wtime()-t_before
  t_before = wtime()
  IF(j==1) then
    call relaxbndrzwguard(f=v(0,0),rhs=res(1,1), &
                          bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
!    call relaxbndrzwguard_jump(f=v(0,0),rhs=res(1,1), maxjump=maxjump(1,1), curjump=j, &
!                               bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
!                               ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  else
    call relaxbndrzwguard_jump(f=v(0,0),rhs=res(1,1), maxjump=maxjump(1,1), curjump=j, &
                               bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                               ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  END if
  t_relax = t_relax + wtime()-t_before
end do
u = u + v
voltf = 1._8
inveps0 = 1./eps0
t_before = wtime()
call apply_voltagewguard(u,bnd,voltf)
t_apply_voltage = t_apply_voltage + wtime()-t_before
t_before = wtime()
call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
t_updateguard = t_updateguard + wtime()-t_before

IF(l_mgridrz_debug) WRITE(0,*) 'exit mg_jump'

return
end subroutine mgbndrzwguard_jump

!pgi$r nobounds
subroutine updateguardcellsrz(f, ixlbnd, ixrbnd, izlbnd, izrbnd)
! update guard cells values according to boundary conditions.
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)
INTEGER(ISZ), optional :: ixlbnd
INTEGER(ISZ) :: ixrbnd, izlbnd, izrbnd

INTEGER(ISZ) :: ixmax, izmax

ixmax=SIZE(f,1)
izmax=SIZE(f,2)

IF(PRESENT(ixlbnd)) then
  select case (ixlbnd)
    case (dirichlet)
      f(1,:) = f(2,:)
    case (neumann)
      f(1,:) = f(3,:)
    case (periodic)
      f(1,:) = f(ixmax-1,:)
    case default
  end select
END if
select case (ixrbnd)
    case (dirichlet)
      f(ixmax,:) = f(ixmax-1,:)
    case (neumann)
      f(ixmax,:) = f(ixmax-2,:)
    case (periodic)
      f(ixmax,:) = f(2,2)
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

IF(solvergeom==RZgeom) f(1,:) = f(3,:)

end subroutine updateguardcellsrz

subroutine apply_voltage(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)
TYPE(bndptr) :: bnd
REAL(8), INTENT(IN) :: coef_voltage

INTEGER(ISZ) :: ic, i

return
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

!pgi$r nobounds
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

INTEGER(ISZ) :: i, ii, ic, j, l, nr, nz, nlocs
REAL(8), allocatable, DIMENSION(:,:) :: uold, uinit
REAL(8), allocatable, DIMENSION(:) :: uold_vlocs, uinit_vlocs
INTEGER(ISZ), allocatable, DIMENSION(:) :: jlocs, klocs
REAL(8) :: maxerr, maxerr_old, t_solve
LOGICAL :: do_calc, has_diverged

IF(l_jump) then
 call solve_multigridrz_jump(grid,accuracy,l_for_timing)
 return
END if

t_relax = 0.
t_restrict = 0.
t_expand = 0.
t_apply_voltage = 0.
t_updateguard = 0.
t_allocate = 0.
t_solve = wtime()

IF(solvergeom==Zgeom) then
  call solve_multigridz(grid)
  return
END if

nr = SIZE(grid%phi,1)
nz = SIZE(grid%phi,2)

IF(vlocs) then
  nlocs = grid%bnd(nlevels)%nvlocs
  do ic = 1, grid%bnd(nlevels)%nb_conductors
    IF(ic==1) then
      grid%bnd(nlevels)%cnd => grid%bnd(nlevels)%first
    else
      grid%bnd(nlevels)%cnd => grid%bnd(nlevels)%cnd%next
    END if
    nlocs = nlocs + grid%bnd(nlevels)%cnd%nbbnd
  END do
  ALLOCATE(uold_vlocs(nlocs),jlocs(nlocs),klocs(nlocs))
  do i = 1, grid%bnd(nlevels)%nvlocs
    j = grid%bnd(nlevels)%vlocs_j(i)
    l = grid%bnd(nlevels)%vlocs_k(i)
    jlocs(i) = j
    klocs(i) = l
  enddo
  i = grid%bnd(nlevels)%nvlocs
  do ic = 1, grid%bnd(nlevels)%nb_conductors
    IF(ic==1) then
      grid%bnd(nlevels)%cnd => grid%bnd(nlevels)%first
    else
      grid%bnd(nlevels)%cnd => grid%bnd(nlevels)%cnd%next
    END if
    do ii = 1, grid%bnd(nlevels)%cnd%nbbnd
      i = i+1
      j = grid%bnd(nlevels)%cnd%jj(ii)
      l = grid%bnd(nlevels)%cnd%kk(ii)
      jlocs(i) = j
      klocs(i) = l
    END do
  END do
  IF(l_for_timing) then
    ALLOCATE(uinit_vlocs(nlocs))
    do i = 1, nlocs
      uinit_vlocs(i) = grid%phi(jlocs(i),klocs(i))
    end do
  END if
else
  ALLOCATE(uold(nr,nz))
  IF(l_for_timing) then
    ALLOCATE(uinit(nr,nz))
    uinit = grid%phi
  END if
END if

do_calc=.true.
has_diverged = .false.

  level = nlevels
  do_calc=.true.
  do while(do_calc)
    maxerr = 1.
    do  j = 1, grid%ncmax
      IF(vlocs) then
        do i = 1, nlocs
          uold_vlocs(i) = grid%phi(jlocs(i),klocs(i))
        end do
      else
        uold=grid%phi
      END if
!      call mgbndrzwguard(j=nlevels,u=grid%phi(0,0),rhs=grid%rho(1,1),bnd=grid%bnd,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
      call mgbndrzwguard(j=nlevels,u=grid%phi,rhs=grid%rho,bnd=grid%bnd,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
                         npre=grid%npre,npost=grid%npost,ncycle=grid%ncycles,sub=.FALSE., relax_only=.false.,npmin=grid%npmin, &
                         mgparam=grid%mgparam)
      maxerr_old = maxerr
      IF(vlocs) then
        maxerr = 0.
        do i = 1, nlocs
          maxerr = max(maxerr,abs(grid%phi(jlocs(i),klocs(i))-uold_vlocs(i)))
        end do
      else
        maxerr = maxval(abs(grid%phi-uold))
      END if
#ifdef MPIPARALLEL
      maxerr = mpi_global_compute_real(maxerr,MPI_MAX)
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
        WRITE(0,*) '        initial maximum error = ',maxerr_old
        WRITE(0,*) '        current maximum error = ',maxerr
        WRITE(0,*) '        trying npre and npost = ',grid%npre+1!,' (also reset mgparam to 1.8)'
#ifdef MPIPARALLEL
      END if
#endif
        grid%npre  = grid%npre+1
        grid%npost = grid%npost+1
!        grid%mgparam = 1.8
        IF(vlocs) then
          IF(l_for_timing) then
            do i = 1, nlocs
              grid%phi(jlocs(i),klocs(i)) = uinit_vlocs(i)
            end do
          else
            do i = 1, nlocs
              grid%phi(jlocs(i),klocs(i)) = uold_vlocs(i)
            end do
          END if
        else
          IF(l_for_timing) then
            grid%phi=uinit
          else
            grid%phi=uold
          END if
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

t_solve = wtime() - t_solve
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

IF(vlocs) then
  DEALLOCATE(uold_vlocs,jlocs,klocs)
  IF(l_for_timing) DEALLOCATE(uinit_vlocs)
else
  DEALLOCATE(uold)
  IF(l_for_timing) DEALLOCATE(uinit)
END if
IF(.not.l_for_timing) then
  IF(l_print_timing) then
    WRITE(0,*) 'Time relax    = ',t_relax
    WRITE(0,*) 'Time restrict = ',t_restrict
    WRITE(0,*) 'Time expand   = ',t_expand
    WRITE(0,*) 'Time apply_voltage = ',t_apply_voltage
    WRITE(0,*) 'Time updateguard   = ',t_updateguard
    WRITE(0,*) 'Time allocate      = ',t_allocate      
    WRITE(0,*) 'Time total    = ',t_relax+t_restrict+t_expand+t_apply_voltage+t_updateguard+t_allocate
    WRITE(0,*) 'Time solve    = ',t_solve
  END if
END if

return
end subroutine solve_multigridrz

subroutine solve_multigridrz_jump(grid,accuracy,l_for_timing)
! solve field for u with density rhoinit.
implicit none

! input/output variables
TYPE(grdptr) :: grid
REAL(8), INTENT(IN) :: accuracy  ! required average accuracy
LOGICAL(ISZ) :: l_for_timing

INTEGER(ISZ) :: i, ii, ic, j, l, nr, nz, nlocs, jm, jp, lm, lp, jmax
REAL(8), allocatable, DIMENSION(:,:) :: uold, uinit
INTEGER(ISZ), DIMENSION(:,:), ALLOCATABLE :: maxjump
REAL(8) :: maxerr, maxerr_old, t_solve
LOGICAL :: do_calc, has_diverged

t_relax = 0.
t_restrict = 0.
t_expand = 0.
t_apply_voltage = 0.
t_updateguard = 0.
t_allocate = 0.
t_solve = wtime()

nr = grid%nr
nz = grid%nz

  ALLOCATE(uold(0:nr+2,0:nz+2),maxjump(nr+1,nz+1))
  IF(l_for_timing) then
    ALLOCATE(uinit(0:nr+2,0:nz+2))
    uinit = grid%phi
  END if

do_calc=.true.
has_diverged = .false.

  maxjump = 1
  do l = 3, nz-1
    do j = 3, nr-1
      i = 1
      jp = j + i
      jm = j - i
      lp = l + i
      lm = l - i
      do WHILE(jm>1 .AND. jp<nr+1 .and. lm>1 .AND. lp<nz+1 .AND. &
               grid%bnd(nlevels)%v(jp,l)==v_vacuum .AND. grid%bnd(nlevels)%v(jm,l)==v_vacuum .AND. &
               grid%bnd(nlevels)%v(j,lp)==v_vacuum .AND. grid%bnd(nlevels)%v(j,lm)==v_vacuum)
        i = i + 1
        jp = j + i
        jm = j - i
        lp = l + i
        lm = l - i
      end do
      maxjump(j,l) = i
    end do
  end do

!  grid%phi = 0.
!  grid%phi(1:nr+1,1:nz+1) = maxjump
!  return

  level = nlevels
  do_calc=.true.
  jmax = MIN(mgridrz_nlevels_max,MAXVAL(maxjump))
!  jmax = 1
  WRITE(0,*) 'jmax = ',jmax
  do while(do_calc)
    maxerr = 1.
    do  j = 1, mgridrz_ncmax
      uold=grid%phi
!      call mgbndrzwguard(j=nlevels,u=grid%phi(0,0),rhs=grid%rho(1,1),bnd=grid%bnd,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
      call mgbndrzwguard_jump(jmax=jmax,u=grid%phi,rhs=grid%rho,maxjump=maxjump, &
                         bnd=grid%bnd(nlevels),nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
                         npre=grid%npre,npost=grid%npost,ncycle=grid%ncycles,sub=.FALSE., relax_only=.false.,npmin=grid%npmin, &
                         mgparam=grid%mgparam)
      maxerr_old = maxerr
      maxerr = maxval(abs(grid%phi-uold))
#ifdef MPIPARALLEL
      maxerr = mpi_global_compute_real(maxerr,MPI_MAX)
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
        WRITE(0,*) '        initial maximum error = ',maxerr_old
        WRITE(0,*) '        current maximum error = ',maxerr
        WRITE(0,*) '        trying npre and npost = ',grid%npre+1!,' (also reset mgparam to 1.8)'
#ifdef MPIPARALLEL
      END if
#endif
        grid%npre  = grid%npre+1
        grid%npost = grid%npost+1
!        grid%mgparam = 1.8
        IF(l_for_timing) then
          grid%phi=uinit
        else
          grid%phi=uold
        END if
        IF(l_for_timing) exit
      END if
    end do
    IF(j>=grid%ncmax) do_calc=.false.
  end do

t_solve = wtime() - t_solve
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

DEALLOCATE(uold,maxjump)
IF(l_for_timing) DEALLOCATE(uinit)

IF(.not.l_for_timing) then
  IF(l_print_timing) then
    WRITE(0,*) 'Time relax    = ',t_relax
    WRITE(0,*) 'Time restrict = ',t_restrict
    WRITE(0,*) 'Time expand   = ',t_expand
    WRITE(0,*) 'Time apply_voltage = ',t_apply_voltage
    WRITE(0,*) 'Time updateguard   = ',t_updateguard
    WRITE(0,*) 'Time allocate      = ',t_allocate      
    WRITE(0,*) 'Time total    = ',t_relax+t_restrict+t_expand+t_apply_voltage+t_updateguard+t_allocate
    WRITE(0,*) 'Time solve    = ',t_solve
  END if
END if

return
end subroutine solve_multigridrz_jump

subroutine solve_multigridz(grid)
! solve field for u with density rhoinit.
implicit none

! input/output variables
TYPE(grdptr) :: grid

INTEGER(ISZ) :: j
REAL(8) :: V0, VL

V0 = grid%phi(1,1)
VL = grid%phi(1,grid%nz+1)

grid%phi(1,2) = 0
do j = 2, grid%nz
  grid%phi(1,2) = grid%phi(1,2) - (j-1)*grid%dz**2*grid%rho(1,grid%nz+2-j)*inveps0
end do
grid%phi(1,2) = ((V0*(grid%nz-1)+VL)-grid%phi(1,2))/grid%nz
do j = 2, grid%nz
  grid%phi(1,j+1) = 2.*grid%phi(1,j) - grid%phi(1,j-1) - grid%dz**2*grid%rho(1,j)*inveps0
end do
grid%phi(1,0) = 2.*grid%phi(1,1)-grid%phi(1,2)
grid%phi(1,grid%nz+2) = 2.*grid%phi(1,grid%nz+1)-grid%phi(1,grid%nz)

return
end subroutine solve_multigridz

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
USE InjectVars_eq
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich, nr0, nz0
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN OUT) :: rho0(nr0+1,nz0+1)
REAL(8), INTENT(IN) :: dr0, dz0, accuracy

  IF(mgridrz_ncmax==0) return

  IF(iwhich==1) return

  IF(l_find_rise_time) then
    call multigridrzf_risetime(iwhich,u0,rho0,nr0,nz0,accuracy)
    return
  END if

!  call distribute_rho(basegrid)
  IF(solvergeom/=Zgeom) then
    IF(solvergeom==XZgeom .AND. basegrid%ixlbnd==dirichlet) &
                                    basegrid%phi(1,:)           = u0(1,:)
    if (basegrid%ixrbnd==dirichlet) basegrid%phi(nr0+1,:)       = u0(nr0+1,:)
  END if
  if (basegrid%izlbnd==dirichlet) basegrid%phi(1:nr0+1,1)     = u0(1:nr0+1,1)
  if (basegrid%izrbnd==dirichlet) basegrid%phi(1:nr0+1,nz0+1) = u0(1:nr0+1,nz0+1)

  call solve_mgridrz(basegrid,accuracy)

  u0(1:nr0+1,:)=basegrid%phi(1:nr0+1,:)

return
end subroutine multigridrzf

subroutine multigridrzf_risetime(iwhich,u0,rho0,nr0,nz0,accuracy)
USE InGen
USE InPart
USE InMesh3d
USE Particles
USE InjectVars
USE InjectVars3D
USE InjectVars_eq, ONLY: inj_phi_eq,v_max,afact,calc_a,l_verbose
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich, nr0, nz0
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN) :: rho0(nr0+1,nz0+1)
REAL(8), INTENT(IN) :: accuracy

REAL(8) :: phi0, phiv, phiref, phirho, wtot
REAL(8), ALLOCATABLE, DIMENSION(:) :: weights
INTEGER(ISZ) :: i, j, max_j, nps_tmp(ns)
INTEGER(ISZ), parameter :: center=1,average_source=2,weighted_average_source=3,border=4

  IF(mgridrz_ncmax==0) return

  IF(iwhich==1) return

  IF(ninject>1) then
    WRITE(0,*) 'ninject>1 not supported by multigridrzf_risetime, stopping.'
    stop
  END if

! --- Calculate the charge density on the surface of the emitter.
!  if (inject == 3) then
!    if(solvergeom==XYZgeom) then
!      call inj_setrho3d(xmmin,ymmin,inj_dx,inj_dy,inj_dz,inj_nx,inj_ny,l2symtry,l4symtry)
!    elseif(solvergeom==RZgeom) then
! --- When using the RZ solver, inj_rho is forced to be
! --- four-fold symmetric.
!      call inj_setrho3d(xmmin,ymmin,inj_dx,inj_dy,inj_dz,inj_nx,inj_ny,.false.,.true.)
!    elseif(solvergeom==Zgeom) then
!      call inj_setrho3d_z(inj_dz,nz)
!    endif
!  endif

  do i = 1, ngrids
    grids_ptr(i)%grid%phi = gridinit(i)%grid%phi
  end do
  phi0 = v_max
  phiref = inj_phi_eq

  IF(inj_nz==1) then
    vinject = 0.
  else
    vinject = v_max
    nps_tmp = nps
    nps = 0
  END if
  call getinj_phi()
  IF(inj_nz/=1) then
    nps = nps_tmp
  END if
  select case (calc_a)
    case (center)
      phiv = vinject(1)-inj_phi(0,0,1)
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy)
      vinject=0.
      call getinj_phi()
      phirho = -inj_phi(0,0,1)
!      IF(inject==3) then
!        phirho = phirho + inj_dz*abs(inj_d(1))*0.5*inj_rho(0,0,1)*inj_dz*inveps0
!        phiref = phiref - inj_dz*abs(inj_d(1))*0.5*inj_rho_eq(0,0,1)*inj_dz*inveps0
!      END if
    case (average_source)
      max_j = 1+INT(ainject(1)/inj_dx)
      phiv = vinject(1)-SUM(inj_phi(0:max_j-1,0,1))/max_j
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy)
      vinject=0.
      call getinj_phi()
      phirho = -SUM(inj_phi(0:max_j-1,0,1))/max_j
    case (weighted_average_source)
      max_j = 1+INT(ainject(1)/inj_dx)
      ALLOCATE(weights(max_j))
      weights(1) = 0.25*pi*inj_dx**2
      do j = 2, max_j
        weights(j) = 2.*pi*(j-1)*inj_dx**2
      end do
      wtot = SUM(weights(1:max_j))
      phiv = vinject(1)-SUM(weights(1:max_j)*inj_phi(0:max_j-1,0,1))/wtot
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy)
      vinject=0.
      call getinj_phi()
      phirho = -SUM(weights(1:max_j)*inj_phi(0:max_j-1,0,1))/wtot
      DEALLOCATE(weights)
    case (border)
      max_j = 1+INT(ainject(1)/inj_dx)
      phiv = vinject(1)-inj_phi(max_j-2,0,1)
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy)
      vinject=0.
      call getinj_phi()
      phirho = -inj_phi(max_j-2,0,1)
    case default
  end select

  afact = (phiref+phirho)/(phi0-phiv)

!  basegrid%phi(1:nr0+1,:) = basegrid%phi(1:nr0+1,:) + a*phi_init(:,1,:)
  do i = 1, ngrids
    grids_ptr(i)%grid%phi = grids_ptr(i)%grid%phi + afact*gridinit(i)%grid%phi
  end do
  if(solvergeom==RZgeom) &
    call updateguardcellsrz(basegrid%phi, basegrid%ixlbnd, basegrid%ixrbnd, basegrid%izlbnd, basegrid%izrbnd)

  vinject=afact*phi0
  call getinj_phi()
  IF(l_verbose) then
    WRITE(0,*) 'a,inj_phi,inj_phi_eq',afact,inj_phi(0,0,1),inj_phi_eq
    WRITE(0,*) 'phi0, phiv, phiref, phirho',phi0, phiv, phiref, phirho
  END if

  u0(1:nr0+1,:)=basegrid%phi(1:nr0+1,:)

return
end subroutine multigridrzf_risetime

subroutine distribute_rho_rz()
USE multigridrz
implicit none

IF(.not.l_distribute) return
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
    IF(solvergeom==zgeom) then
      call deposit_z(unew=grid%up%rho(1,:), uold=grid%rho(1,:), invvolnew=1./grid%up%dz, invvolold=1./grid%dz, &
                     zminold=grid%zmin, zmaxold=grid%zmax, zminnew=grid%up%zmin, zmaxnew=grid%up%zmax)
    else
      call deposit(unew=grid%up%rho, uold=grid%rho, invvolnew=grid%up%invvol, invvolold=grid%invvol, &
                   xminold=grid%rmin, xmaxold=grid%rmax, zminold=grid%zmin, zmaxold=grid%zmax, &
                   xminnew=grid%up%rmin, xmaxnew=grid%up%rmax, zminnew=grid%up%zmin, zmaxnew=grid%up%zmax)
    END if
  END if

return
END subroutine distribute_rho

RECURSIVE subroutine solve_mgridrz(grid,accuracy)
USE multigridrz
implicit none
TYPE(grdptr) :: grid

TYPE(bndptr), POINTER :: bndlocal
REAL(8), INTENT(IN) :: accuracy

INTEGER(ISZ) :: i, ic

!    grid%mgparam=grid%mgparam+0.05

    IF(associated(grid%up)) then
      IF(solvergeom==Zgeom) then
        CALL interpolate_any_1d(unew=grid%phi(1,:),uold=grid%up%phi(1,:), &
                                nznew=grid%nz, nzold=grid%up%nz, &
                                zminold=grid%up%zmin, zmaxold=grid%up%zmax, &
                                zminnew=grid%zmin, zmaxnew=grid%zmax, &
                                izlbnd=grid%izlbnd, &
                                izrbnd=grid%izrbnd, &
                                bnd_only=.false., quad=.false.)
      else
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
                             bnd_only=.false., quad=.false.)
!                             bnd_only=.true., quad=.true.)
        bndlocal => grid%bnd(grid%nlevels)
        do ic = 1, bndlocal%nb_conductors
          IF(ic==1) then
            bndlocal%cnd => bndlocal%first
          else
            bndlocal%cnd => bndlocal%cnd%next
          END if
          do i = 1, bndlocal%cnd%ncond
            IF(bndlocal%cnd%condid(i)<0) bndlocal%cnd%voltage(i) = grid%phi(bndlocal%cnd%jcond(i),bndlocal%cnd%kcond(i))
          end do
        end do
      END if
    END if
    nlevels = grid%nlevels
    level = nlevels
    IF(solvergeom/=Zgeom) bndy => grid%bnd
    ixlbnd = grid%ixlbnd
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
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry_in,l4symtry_in,condid)
! call subroutine srfrvout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use Conductor3d
USE multigridrz
implicit none
character(*) rofzfunc
real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
LOGICAL(ISZ):: lfill,lshell,l2symtry_in,l4symtry_in
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis
real(kind=8):: xmesh(0:nx),ymesh(0:ny)
integer(ISZ):: condid

INTEGER(ISZ) :: i,nrc,nzc,igrid,jmin,jmax
REAL(8) :: drc,dzc,zmin_in,zmax_in
LOGICAL(ISZ) :: doloop

TYPE(conductor_type), POINTER :: cndpnt

IF(solvergeom==Zgeom) return

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 bndy => grids_ptr(igrid)%grid%bnd
 ixlbnd = grids_ptr(igrid)%grid%ixlbnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = nlevels,1,-1
  IF(.not.bndy(i)%l_lshift) then
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
                       ix_axis,iy_axis,xmesh,ymesh,l2symtry_in,l4symtry_in,condid)

   call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                         ncond, ixcond, izcond, condvolt, condnumb, &
                         necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                         ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                         ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                         nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz, &
                         ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                         ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)
  
  else
   jmin = 1
   jmax = 1
   doloop = .true.
   do WHILE(doloop)
    nrc = 1
    do WHILE(bndy(i)%lshift(jmax)==bndy(i)%lshift(jmax+1) .and. jmax<bndy(i)%nr)
     jmax = jmax + 1
     nrc = nrc+1
    end do
    nzc = bndy(i)%nz
    drc = bndy(i)%dr
    dzc = bndy(i)%dz
    izlbnd = bndy(i)%izlbnd
    izrbnd = bndy(i)%izrbnd
    zmin_in = bndy(i)%zmin+bndy(i)%lshift(jmin)*bndy(i)%dz
    zmax_in = bndy(i)%zmax+bndy(i)%lshift(jmin)*bndy(i)%dz

    necndbdy=0
    nocndbdy=0
    ncond = 0

    call srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,  &
                        xmin,xmax,ymin,ymax,lshell,                      &
                        zmin_in,zmax_in,zbeam,drc,drc,dzc,nrc,ny,nzc,             &
                        ix_axis,iy_axis,xmesh,ymesh,l2symtry_in,l4symtry_in,condid)

    call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                          ncond, ixcond, izcond, condvolt, condnumb, &
                          necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                          ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                          ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                          nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz, &
                          ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                          ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)

    IF(jmax==bndy(i)%nr) then
     doloop = .false.
    else
     jmin = jmax + 1
     jmax = jmax + 1
    END if
   end do
  END if
 end do
end do
necndbdy=0
nocndbdy=0
ncond = 0

return
end subroutine srfrvoutrz

subroutine srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,   &
                           lzend,xmin,xmax,ymin,ymax,lshell,          &
                           zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,       &
                           ix_axis,iy_axis,xmesh,ymesh,l2symtry_in,l4symtry_in,condid)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use Conductor3d
USE multigridrz
implicit none
character(*) rminofz,rmaxofz
real(kind=8):: volt,zmin,zmax,xcent,ycent
LOGICAL(ISZ):: lzend,lshell,l2symtry_in,l4symtry_in
real(kind=8):: xmin,xmax,ymin,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,iy_axis,jmin,jmax
real(kind=8):: xmesh(0:nx),ymesh(0:ny)
integer(ISZ):: condid
LOGICAL(ISZ) :: doloop

INTEGER(ISZ) :: i,nrc,nzc,igrid
REAL(8) :: drc,dzc,zmin_in,zmax_in

TYPE(conductor_type), POINTER :: cndpnt

IF(solvergeom==Zgeom) return

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 bndy => grids_ptr(igrid)%grid%bnd
 ixlbnd = grids_ptr(igrid)%grid%ixlbnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = nlevels,1,-1
  IF(.not.bndy(i)%l_lshift) then
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
                      ix_axis,iy_axis,xmesh,ymesh,l2symtry_in,l4symtry_in,condid)

   call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                         ncond, ixcond, izcond, condvolt, condnumb, &
                         necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                         ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                         ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                         nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz, &
                         ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                         ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)
  
  else
   jmin = 1
   jmax = 1
   doloop = .true.
   do WHILE(doloop)
    nrc = 1
    do WHILE(bndy(i)%lshift(jmax)==bndy(i)%lshift(jmax+1) .and. jmax<bndy(i)%nr)
     jmax = jmax + 1
     nrc = nrc+1
    end do
    nzc = bndy(i)%nz
    drc = bndy(i)%dr
    dzc = bndy(i)%dz
    izlbnd = bndy(i)%izlbnd
    izrbnd = bndy(i)%izrbnd
    zmin_in = bndy(i)%zmin+bndy(i)%lshift(jmin)*bndy(i)%dz
    zmax_in = bndy(i)%zmax+bndy(i)%lshift(jmin)*bndy(i)%dz

    necndbdy=0
    nocndbdy=0
    ncond = 0

    call srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,ycent,  &
                       lzend,xmin,xmax,ymin,ymax,lshell,                &
                       zmin_in,zmax_in,zbeam,drc,drc,dzc,nrc,ny,nzc,             &
                       ix_axis,iy_axis,xmesh,ymesh,l2symtry_in,l4symtry_in,condid)

    call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                          ncond, ixcond, izcond, condvolt, condnumb, &
                          necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                          ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                          ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                          nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz, &
                          ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                          ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)

    IF(jmax==bndy(i)%nr) then
     doloop = .false.
    else
     jmin = jmax + 1
     jmax = jmax + 1
    END if
   end do
  END if
 end do
end do
necndbdy=0
nocndbdy=0
ncond = 0
return
end subroutine srfrvinoutrz

subroutine setcndtrrz(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz, &
                      l2symtry_in,l4symtry_in)
use Conductor3d
USE multigridrz
integer(ISZ):: nx,ny,nz
real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
logical(ISZ):: l2symtry_in,l4symtry_in

INTEGER(ISZ) :: i,nrc,nzc,igrid,jmin,jmax
REAL(8) :: drc,dzc,rmin_in,zmin_in
LOGICAL(ISZ) :: doloop

IF(solvergeom==Zgeom) return

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 bndy => grids_ptr(igrid)%grid%bnd
 ixlbnd = grids_ptr(igrid)%grid%ixlbnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = nlevels,1,-1
  IF(.not.bndy(i)%l_lshift) then
   nrc = bndy(i)%nr
   nzc = bndy(i)%nz
   drc = bndy(i)%dr
   dzc = bndy(i)%dz
   izlbnd = bndy(i)%izlbnd
   izrbnd = bndy(i)%izrbnd
   zmin_in = bndy(i)%zmin!-bndy(i)%dz

   necndbdy=0
   nocndbdy=0
   ncond = 0

   call setcndtr_rz(rmin_in,rmin_in,zmin_in,zbeam,zgrid,nrc,ny,nzc,drc,drc,dzc, &
                    l2symtry_in,l4symtry_in)

   call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                         ncond, ixcond, izcond, condvolt, condnumb, &
                         necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                         ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                         ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                         nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz, &
                         ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                         ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)
  
  else
   jmin = 1
   jmax = 1
   doloop = .true.
   do WHILE(doloop)
    nrc = 1
    do WHILE(bndy(i)%lshift(jmax)==bndy(i)%lshift(jmax+1) .and. jmax<bndy(i)%nr)
     jmax = jmax + 1
     nrc = nrc+1
    end do
    nzc = bndy(i)%nz
    drc = bndy(i)%dr
    dzc = bndy(i)%dz
    izlbnd = bndy(i)%izlbnd
    izrbnd = bndy(i)%izrbnd
    zmin_in = bndy(i)%zmin+bndy(i)%lshift(jmin)*bndy(i)%dz

    necndbdy=0
    nocndbdy=0
    ncond = 0

    call setcndtr_rz(rmin_in,rmin_in,zmin_in,zbeam,zgrid,nrc,ny,nzc,drc,drc,dzc, &
                     l2symtry_in,l4symtry_in)

    call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                          ncond, ixcond, izcond, condvolt, condnumb, &
                          necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                          ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                          ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                          nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz, &
                          ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                          ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)

    IF(jmax==bndy(i)%nr) then
     doloop = .false.
    else
     jmin = jmax + 1
     jmax = jmax + 1
    END if
   end do
  END if
 end do
end do
necndbdy=0
nocndbdy=0
ncond = 0

END subroutine setcndtrrz

subroutine install_conductors_rz()
USE Conductor3d
USE Multigrid3d
USE multigridrz
implicit none

INTEGER(ISZ), DIMENSION(:), allocatable :: mg_ncond,mg_necndbdy, mg_nocndbdy
INTEGER(ISZ) :: i, ii, nrc, nzc, itmp
INTEGER(ISZ) :: ncondtmp, necndbdytmp, nocndbdytmp
REAL(8) :: drc, dzc

INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: ixcondtmp, izcondtmp, &
                                           iecndxtmp, iecndztmp, &
                                           iocndxtmp, iocndztmp, &
                                           condnumbtmp, &
                                           ecnumbmxtmp, ecnumbpxtmp, &
                                           ecnumbmztmp, ecnumbpztmp, &
                                           ocnumbmxtmp, ocnumbpxtmp, &
                                           ocnumbmztmp, ocnumbpztmp
REAL(8), ALLOCATABLE, DIMENSION(:) :: ecdelmxtmp, ecdelpxtmp, ecdelmztmp, ecdelpztmp, &
                                      ocdelmxtmp, ocdelpxtmp, ocdelmztmp, ocdelpztmp, &
                                      ecvoltmxtmp, ecvoltpxtmp, ecvoltmztmp, ecvoltpztmp, &
                                      ocvoltmxtmp, ocvoltpxtmp, ocvoltmztmp, ocvoltpztmp, &
                                      condvolttmp

IF(solvergeom==Zgeom) return

ALLOCATE(mg_ncond(basegrid%nlevels),mg_necndbdy(basegrid%nlevels), mg_nocndbdy(basegrid%nlevels))
mg_ncond = 0
mg_necndbdy = 0
mg_nocndbdy = 0

ixlbnd = basegrid%ixlbnd
ixrbnd = basegrid%ixrbnd
do i = 1, ncond
  ii = basegrid%nlevels-icondlevel(i)
  mg_ncond(ii) = mg_ncond(ii) + 1
end do
do i = 1, necndbdy
  ii = basegrid%nlevels-iecndlevel(i)
  mg_necndbdy(ii) = mg_necndbdy(ii) + 1
end do
do i = 1, nocndbdy
  ii = basegrid%nlevels-iocndlevel(i)
  mg_nocndbdy(ii) = mg_nocndbdy(ii) + 1
end do

 do i = nlevels,1,-1
  level = i
  nrc = bndy(i)%nr
  nzc = bndy(i)%nz
  drc = bndy(i)%dr
  dzc = bndy(i)%dz
  izlbnd = bndy(i)%izlbnd
  izrbnd = bndy(i)%izrbnd

  necndbdytmp = mg_necndbdy(i)
  nocndbdytmp = mg_nocndbdy(i)
  ncondtmp    = mg_ncond(i)

  ALLOCATE(ixcondtmp(ncondtmp),izcondtmp(ncondtmp),condvolttmp(ncondtmp),condnumbtmp(ncondtmp), &
           iecndxtmp(necndbdytmp),iecndztmp(necndbdytmp), &
           iocndxtmp(nocndbdytmp),iocndztmp(nocndbdytmp), &
           ecnumbmxtmp(necndbdytmp),ecnumbpxtmp(necndbdytmp), ecnumbmztmp(necndbdytmp),ecnumbpztmp(necndbdytmp), &
           ocnumbmxtmp(nocndbdytmp),ocnumbpxtmp(nocndbdytmp), ocnumbmztmp(nocndbdytmp),ocnumbpztmp(nocndbdytmp), &
           ecdelmxtmp(necndbdytmp),ecdelpxtmp(necndbdytmp), ecdelmztmp(necndbdytmp),ecdelpztmp(necndbdytmp), &
           ocdelmxtmp(nocndbdytmp),ocdelpxtmp(nocndbdytmp), ocdelmztmp(nocndbdytmp),ocdelpztmp(nocndbdytmp), &
           ecvoltmxtmp(necndbdytmp),ecvoltpxtmp(necndbdytmp), ecvoltmztmp(necndbdytmp),ecvoltpztmp(necndbdytmp), &
           ocvoltmxtmp(nocndbdytmp),ocvoltpxtmp(nocndbdytmp), ocvoltmztmp(nocndbdytmp),ocvoltpztmp(nocndbdytmp))
  itmp = 0
  do ii = 1, ncond
    IF(basegrid%nlevels-icondlevel(ii)==i) then
      itmp = itmp + 1
      ixcondtmp(itmp) = ixcond(ii)
      izcondtmp(itmp) = izcond(ii)
      condvolttmp(itmp) = condvolt(ii)
      condnumbtmp(itmp) = condnumb(ii)
    END if
  end do
  itmp = 0
  do ii = 1, necndbdy
    IF(basegrid%nlevels-iecndlevel(ii)==i) then
      itmp = itmp + 1
      iecndxtmp(itmp) = iecndx(ii)
      iecndztmp(itmp) = iecndz(ii)
      ecdelmxtmp(itmp) = ecdelmx(ii)
      ecdelpxtmp(itmp) = ecdelpx(ii)
      ecdelmztmp(itmp) = ecdelmz(ii)
      ecdelpztmp(itmp) = ecdelpz(ii)
      ecvoltmxtmp(itmp) = ecvoltmx(ii)
      ecvoltpxtmp(itmp) = ecvoltpx(ii)
      ecvoltmztmp(itmp) = ecvoltmz(ii)
      ecvoltpztmp(itmp) = ecvoltpz(ii)
      ecnumbmxtmp(itmp) = ecnumbmx(ii)
      ecnumbpxtmp(itmp) = ecnumbpx(ii)
      ecnumbmztmp(itmp) = ecnumbmz(ii)
      ecnumbpztmp(itmp) = ecnumbpz(ii)
    END if
  end do
  itmp = 0
  do ii = 1, nocndbdy
    IF(basegrid%nlevels-iocndlevel(ii)==i) then
      itmp = itmp + 1
      iocndxtmp(itmp) = iocndx(ii)
      iocndztmp(itmp) = iocndz(ii)
      ocdelmxtmp(itmp) = ocdelmx(ii)
      ocdelpxtmp(itmp) = ocdelpx(ii)
      ocdelmztmp(itmp) = ocdelmz(ii)
      ocdelpztmp(itmp) = ocdelpz(ii)
      ocvoltmxtmp(itmp) = ocvoltmx(ii)
      ocvoltpxtmp(itmp) = ocvoltpx(ii)
      ocvoltmztmp(itmp) = ocvoltmz(ii)
      ocvoltpztmp(itmp) = ocvoltpz(ii)
      ocnumbmxtmp(itmp) = ocnumbmx(ii)
      ocnumbpxtmp(itmp) = ocnumbpx(ii)
      ocnumbmztmp(itmp) = ocnumbmz(ii)
      ocnumbpztmp(itmp) = ocnumbpz(ii)
    END if
  end do
  call addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                        ncondtmp, ixcondtmp, izcondtmp, condvolttmp, condnumbtmp, &
                        necndbdytmp, iecndxtmp, iecndztmp, ecdelmxtmp, ecdelpxtmp, ecdelmztmp, ecdelpztmp, &
                        ecvoltmxtmp,ecvoltpxtmp, ecvoltmztmp,ecvoltpztmp, &
                        ecnumbmxtmp,ecnumbpxtmp, ecnumbmztmp,ecnumbpztmp, &
                        nocndbdytmp, iocndxtmp, iocndztmp, ocdelmxtmp, ocdelpxtmp, ocdelmztmp, ocdelpztmp, &
                        ocvoltmxtmp,ocvoltpxtmp, ocvoltmztmp,ocvoltpztmp, &
                        ocnumbmxtmp,ocnumbpxtmp, ocnumbmztmp,ocnumbpztmp)

  DEALLOCATE(ixcondtmp,izcondtmp,condvolttmp,condnumbtmp, &
             iecndxtmp,iecndztmp, &
             iocndxtmp,iocndztmp, &
             ecnumbmxtmp,ecnumbpxtmp, ecnumbmztmp,ecnumbpztmp, &
             ocnumbmxtmp,ocnumbpxtmp, ocnumbmztmp,ocnumbpztmp, &
             ecdelmxtmp,ecdelpxtmp, ecdelmztmp,ecdelpztmp, &
             ocdelmxtmp,ocdelpxtmp, ocdelmztmp,ocdelpztmp, &
             ecvoltmxtmp,ecvoltpxtmp, ecvoltmztmp,ecvoltpztmp, &
             ocvoltmxtmp,ocvoltpxtmp, ocvoltmztmp,ocvoltpztmp)

 end do

DEALLOCATE(mg_ncond,mg_necndbdy, mg_nocndbdy)

necndbdy=0
nocndbdy=0
ncond = 0

call get_cond_rz(1,1)

return
end subroutine install_conductors_rz


subroutine addconductors_rz(i,nrc,nzc,drc,dzc,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                            ncond, ixcond, izcond, condvolt, condnumb, &
                            necndbdy, iecndx, iecndz, ecdelmx, ecdelpx, ecdelmz, ecdelpz, &
                            ecvoltmx, ecvoltpx, ecvoltmz, ecvoltpz, &
                            ecnumbmx, ecnumbpx, ecnumbmz, ecnumbpz, &
                            nocndbdy, iocndx, iocndz, ocdelmx, ocdelpx, ocdelmz, ocdelpz,  &
                            ocvoltmx, ocvoltpx, ocvoltmz, ocvoltpz, &
                            ocnumbmx, ocnumbpx, ocnumbmz, ocnumbpz)
!use Conductor3d
USE InGen3d, ONLY:solvergeom,RZgeom,XYZgeom,XZgeom,Zgeom
USE multigridrz, ONLY: conductor_type, bndy, dirichlet, v_cond, v_bnd, v_dirichlet, bnd_method, egun, ecb, init_bnd_sublevel
implicit none

INTEGER(ISZ), INTENT(IN) :: nrc,nzc,i,ixlbnd,ixrbnd,izlbnd,izrbnd,ncond,necndbdy,nocndbdy
REAL(8), INTENT(IN) :: drc,dzc
INTEGER(ISZ), INTENT(IN) :: ixcond(ncond), izcond(ncond), &
                            iecndx(necndbdy), iecndz(necndbdy), &
                            iocndx(nocndbdy), iocndz(nocndbdy), &
                            condnumb(ncond), &
                            ecnumbmx(necndbdy), ecnumbpx(necndbdy), ecnumbmz(necndbdy), ecnumbpz(necndbdy), &
                            ocnumbmx(nocndbdy), ocnumbpx(nocndbdy), ocnumbmz(nocndbdy), ocnumbpz(nocndbdy)
REAL(8), INTENT(IN) :: condvolt(ncond), &
                       ecdelmx(necndbdy), ecdelpx(necndbdy), ecdelmz(necndbdy), ecdelpz(necndbdy), &
                       ocdelmx(nocndbdy), ocdelpx(nocndbdy), ocdelmz(nocndbdy), ocdelpz(nocndbdy), &
                       ecvoltmx(necndbdy), ecvoltpx(necndbdy), ecvoltmz(necndbdy), ecvoltpz(necndbdy), &
                       ocvoltmx(nocndbdy), ocvoltpx(nocndbdy), ocvoltmz(nocndbdy), ocvoltpz(nocndbdy)

INTEGER(ISZ) :: ii,iii,iv,iiv,nxbndmin,nxbndmax,nzbndmin,nzbndmax,iivmin,iivmax,ibnd,ne,no
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz

TYPE(conductor_type), POINTER :: cndpnt

  nxbndmin=0
  nxbndmax=nrc
  nzbndmin=0
  nzbndmax=nzc
  IF(ixlbnd==dirichlet) nxbndmin=nxbndmin+1
  IF(ixrbnd==dirichlet) nxbndmax=nxbndmax-1
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
    bndy(i)%cnd%condid(iii) = condnumb(ii)
  end do

  ii = 0
  ne = 0
  no = 0
  do ibnd = 1, necndbdy+nocndbdy
   ii = ii + 1
   IF(ibnd<=necndbdy) then
     iii = ibnd
     bndy(i)%cnd%jj(ii)  = iecndx(iii)+1
     bndy(i)%cnd%kk(ii)  = iecndz(iii)+1
   else
     iii = ibnd - necndbdy
     bndy(i)%cnd%jj(ii)  = iocndx(iii)+1
     bndy(i)%cnd%kk(ii)  = iocndz(iii)+1
   END if
   IF( (bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii))==v_cond) .or.                 &
     (bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii))==v_dirichlet) .or.              &
     (.not. (bndy(i)%cnd%jj(ii)>=nxbndmin+1 .and. bndy(i)%cnd%jj(ii)<=nxbndmax+1 .and. &
             bndy(i)%cnd%kk(ii)>=nzbndmin+1 .and. bndy(i)%cnd%kk(ii)<=nzbndmax+1))) then
      ii = ii - 1
      cycle
   END if
   IF(ibnd<=necndbdy) then
     ne = ne + 1
     dxm = MIN(1._8,ecdelmx(iii))*bndy(i)%dr
     dxp = MIN(1._8,ecdelpx(iii))*bndy(i)%dr
     dzm = MIN(1._8,ecdelmz(iii))*bndy(i)%dz
     dzp = MIN(1._8,ecdelpz(iii))*bndy(i)%dz
     bndy(i)%cnd%volt0xm(ii)=ecvoltmx(iii)
     bndy(i)%cnd%volt0xp(ii)=ecvoltpx(iii)
     bndy(i)%cnd%volt0zm(ii)=ecvoltmz(iii)
     bndy(i)%cnd%volt0zp(ii)=ecvoltpz(iii)
     bndy(i)%cnd%condidxm(ii)=ecnumbmx(iii)
     bndy(i)%cnd%condidxp(ii)=ecnumbpx(iii)
     bndy(i)%cnd%condidzm(ii)=ecnumbmz(iii)
     bndy(i)%cnd%condidzp(ii)=ecnumbpz(iii)
   else
     no = no + 1
     dxm = MIN(1._8,ocdelmx(iii))*bndy(i)%dr
     dxp = MIN(1._8,ocdelpx(iii))*bndy(i)%dr
     dzm = MIN(1._8,ocdelmz(iii))*bndy(i)%dz
     dzp = MIN(1._8,ocdelpz(iii))*bndy(i)%dz
     bndy(i)%cnd%volt0xm(ii)=ocvoltmx(iii)
     bndy(i)%cnd%volt0xp(ii)=ocvoltpx(iii)
     bndy(i)%cnd%volt0zm(ii)=ocvoltmz(iii)
     bndy(i)%cnd%volt0zp(ii)=ocvoltpz(iii)
     bndy(i)%cnd%condidxm(ii)=ocnumbmx(iii)
     bndy(i)%cnd%condidxp(ii)=ocnumbpx(iii)
     bndy(i)%cnd%condidzm(ii)=ocnumbmz(iii)
     bndy(i)%cnd%condidzp(ii)=ocnumbpz(iii)
   END if
   bndy(i)%cnd%docalc(ii)=.true.
   IF(bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii))/=v_bnd ) then
     bndy(i)%v(bndy(i)%cnd%jj(ii),bndy(i)%cnd%kk(ii)) = v_bnd
   else
     do iv=1, bndy(i)%nb_conductors-1
       IF(iv==1) then
         cndpnt => bndy(i)%first
       else
         cndpnt => cndpnt%next
       END if
       IF(ibnd<=necndbdy) then
         iivmin = 1
         iivmax = cndpnt%nbbndred
       else
         iivmin = cndpnt%nbbndred+1
         iivmax = cndpnt%nbbnd
       END if
       do iiv=iivmin,iivmax
         IF(bndy(i)%cnd%jj(ii)==cndpnt%jj(iiv) .AND. bndy(i)%cnd%kk(ii)==cndpnt%kk(iiv)) then
           cndpnt%docalc(iiv)=.false.
           IF(cndpnt%dxm(iiv)<dxm) then
             dxm = cndpnt%dxm(iiv)
             bndy(i)%cnd%volt0xm(ii) = cndpnt%volt0xm(iiv)
             bndy(i)%cnd%condidxm(ii) = cndpnt%condidxm(iiv)
           END if
           IF(cndpnt%dxp(iiv)<dxp) then
             dxp = cndpnt%dxp(iiv)
             bndy(i)%cnd%volt0xp(ii) = cndpnt%volt0xp(iiv)
             bndy(i)%cnd%condidxp(ii) = cndpnt%condidxp(iiv)
           END if
           IF(cndpnt%dzm(iiv)<dzm) then
             dzm = cndpnt%dzm(iiv)
             bndy(i)%cnd%volt0zm(ii) = cndpnt%volt0zm(iiv)
             bndy(i)%cnd%condidzm(ii) = cndpnt%condidzm(iiv)
           END if
           IF(cndpnt%dzp(iiv)<dzp) then
             dzp = cndpnt%dzp(iiv)
             bndy(i)%cnd%volt0zp(ii) = cndpnt%volt0zp(iiv)
             bndy(i)%cnd%condidzp(ii) = cndpnt%condidzp(iiv)
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
   IF(solvergeom==RZgeom) then
    IF(bndy(i)%cnd%jj(ii)==1) then
     bndy(i)%cnd%dt(ii) = 1._8/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = 0.
     bndy(i)%cnd%cfxp(ii) = 4._8/(dxp*dxx)
    else
     r = (bndy(i)%cnd%jj(ii)-1)*bndy(i)%dr
     select case (bnd_method)
       case (egun)
         rm = r-0.5_8*bndy(i)%dr
         rp = r+0.5_8*bndy(i)%dr
       case (ecb)
         rm = r-0.5_8*dxm
         rp = r+0.5_8*dxp
       case default
     end select
     bndy(i)%cnd%dt(ii) = 1._8/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = rm/(r*dxm*dxx)
     bndy(i)%cnd%cfxp(ii) = rp/(r*dxp*dxx)
    END if
   else ! (solvergeom==XZgeom)
     bndy(i)%cnd%dt(ii) = 1._8/((1._8/dxm+1._8/dxp)/dxx+(1._8/dzm+1._8/dzp)/dzz)
     bndy(i)%cnd%cfxm(ii) = 1._8/(dxm*dxx)
     bndy(i)%cnd%cfxp(ii) = 1._8/(dxp*dxx)
   END if
   bndy(i)%cnd%cfzm(ii) = 1._8/(dzm*dzz)
   bndy(i)%cnd%cfzp(ii) = 1._8/(dzp*dzz)
   bndy(i)%cnd%cf0(ii)  = -bndy(i)%cnd%cfxm(ii)-bndy(i)%cnd%cfxp(ii)-bndy(i)%cnd%cfzm(ii)-bndy(i)%cnd%cfzp(ii)
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
  bndy(i)%cnd%nbbndred = ne
  bndy(i)%cnd%nbbnd = ne + no

end subroutine addconductors_rz

!=============================================================================
subroutine gtlchgrz
USE multigridrz
use Subtimers3d
use Picglb
use InGen3d
use Picglb3d
use InMesh3d
use Fields3d
use Z_arrays
use InDiag3d

!  Calculates the line charge density from rho.

real(kind=8):: zz,wzg,dzi
integer(ISZ):: j,iz,izg
real(kind=8):: substarttime
if (lw3dtimesubs) substarttime = wtime()
if (.not. lgtlchg3d) return

dzi = 1./dz

!conversion factor to go from grid frame to beam frame
zz = zgrid + zmmin - zzmin - zbeam

do iz=0,nzzarr

  izg = (iz*dzz - zz)*dzi
  wzg = (iz*dzz - zz)*dzi - izg
  izg = izg+1

  if (1 <= izg .and. izg <= nz+1) then
   do j = 1, basegrid%nr+1
    linechg(iz) = basegrid%rho(j,izg)/(basegrid%invvol(j)*basegrid%dz)*(1. - wzg)
   end do
  else
    linechg(iz) = 0.
  endif
  if (1 <= izg+1 .and. izg+1 <= nz+1) then
   do j = 1, basegrid%nr+1
    linechg(iz) = linechg(iz) +  basegrid%rho(j,izg+1)/(basegrid%invvol(j)*basegrid%dz)*wzg
   end do
  endif

enddo

if (lw3dtimesubs) timegtlchg3d = timegtlchg3d + wtime() - substarttime
return
end subroutine gtlchgrz

subroutine dep_rho_rz(is,rho,nr,nz,dr,dz,xmin,zmin)
USE constant
use Particles
implicit none

INTEGER(ISZ), INTENT(IN) :: is, nr, nz
REAL(8), INTENT(IN OUT) :: rho(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, zmin, xmin

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr), invvolxz
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp
REAL(8):: q, qw
LOGICAL(ISZ) :: l_sym

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

  q = sq(is)*sw(is)

  ! make charge deposition using CIC weighting
  IF(wpid==0) then
   do i = ins(is), ins(is) + nps(is) - 1
    IF(uzp(i)==0.) cycle
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr
    zpos = (zp(i)-zmin)*invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    rho(jn, ln)  = rho(jn, ln)  + q * oddr * oddz * invvol(jn-1)
    rho(jnp,ln)  = rho(jnp,ln)  + q *  ddr * oddz * invvol(jnp-1)
    rho(jn, lnp) = rho(jn, lnp) + q * oddr *  ddz * invvol(jn-1)
    rho(jnp,lnp) = rho(jnp,lnp) + q *  ddr *  ddz * invvol(jnp-1)
   end do
  else
   do i = ins(is), ins(is) + nps(is) - 1
    IF(uzp(i)==0.) cycle
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))*invdr
    zpos = (zp(i)-zmin)*invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    qw = q*pid(i,wpid)
    rho(jn, ln)  = rho(jn, ln)  + qw * oddr * oddz * invvol(jn-1)
    rho(jnp,ln)  = rho(jnp,ln)  + qw *  ddr * oddz * invvol(jnp-1)
    rho(jn, lnp) = rho(jn, lnp) + qw * oddr *  ddz * invvol(jn-1)
    rho(jnp,lnp) = rho(jnp,lnp) + qw *  ddr *  ddz * invvol(jnp-1)
   end do
  END if
return
END SUBROUTINE dep_rho_rz

!     ******************************************************************
!     *
!     *                        SUBROUTINE RHOWEIGHTRZ
!     *
!     ******************************************************************


subroutine rhoweightrz(xp,yp,zp,uzp,np,q,nr,nz,dr,dz,xmin,zmin)
USE multigridrz
USE Subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
REAL(8), INTENT(IN) :: q, dr, dz, zmin, xmin

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr), invvolxz
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp
REAL(8):: substarttime
LOGICAL(ISZ) :: l_sym

if (lw3dtimesubs) substarttime = wtime()

IF(np==0) return

IF(solvergeom==RZgeom) then
 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightrz_meshref(xp,yp,zp,uzp,np,q)
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
else ! IF(solvergeom==XZgeom) then
 IF(l4symtry .OR. l2symtry) then
   l_sym = .true.
 else
   l_sym = .false.
 END if
 IF(ngrids>1) then
!  call rhoweightrznew(xp,yp,zp,uzp,np,q)
  WRITE(0,*) 'mesh refinement not yet supported in XZ.'
  stop
 else
  invdr = 1._8/dr
  invdz = 1._8/dz

  invvolxz = 1._8 / (dr*dz)

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    IF(l_sym) then
      rpos = abs(xp(i))*invdr
    else
      rpos = (xp(i)-xmin)*invdr
    END if
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
    basegrid%rhop(jn, ln)  = basegrid%rhop(jn, ln)  + q * oddr * oddz * invvolxz
    basegrid%rhop(jnp,ln)  = basegrid%rhop(jnp,ln)  + q *  ddr * oddz * invvolxz
    basegrid%rhop(jn, lnp) = basegrid%rhop(jn, lnp) + q * oddr *  ddz * invvolxz
    basegrid%rhop(jnp,lnp) = basegrid%rhop(jnp,lnp) + q *  ddr *  ddz * invvolxz
    basegrid%rhominr = MIN(basegrid%rhominr,jn)
    basegrid%rhomaxr = MAX(basegrid%rhomaxr,jnp)
    basegrid%rhominz = MIN(basegrid%rhominz,ln)
    basegrid%rhomaxz = MAX(basegrid%rhomaxz,lnp)
#else
    basegrid%rho(jn, ln)  = basegrid%rho(jn, ln)  + q * oddr * oddz * invvolxz
    basegrid%rho(jnp,ln)  = basegrid%rho(jnp,ln)  + q *  ddr * oddz * invvolxz
    basegrid%rho(jn, lnp) = basegrid%rho(jn, lnp) + q * oddr *  ddz * invvolxz
    basegrid%rho(jnp,lnp) = basegrid%rho(jnp,lnp) + q *  ddr *  ddz * invvolxz
#endif
  end do
 END if
END if

if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
return
END SUBROUTINE RHOWEIGHTRZ

subroutine rhoweightrz_weights(xp,yp,zp,uzp,w,np,q,nr,nz,dr,dz,xmin,zmin)
USE multigridrz
USE Subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp, w
REAL(8), INTENT(IN) :: q, dr, dz, zmin, xmin

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr), invvolxz, qw
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp
REAL(8):: substarttime
LOGICAL(ISZ) :: l_sym

if (lw3dtimesubs) substarttime = wtime()

IF(np==0) return

IF(solvergeom==RZgeom) then
 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightrz_meshref_weights(xp,yp,zp,uzp,w,np,q)
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
    qw = q*w(i)
#ifdef MPIPARALLEL
    basegrid%rhop(jn, ln)  = basegrid%rhop(jn, ln)  + qw * oddr * oddz * invvol(jn-1)
    basegrid%rhop(jnp,ln)  = basegrid%rhop(jnp,ln)  + qw *  ddr * oddz * invvol(jnp-1)
    basegrid%rhop(jn, lnp) = basegrid%rhop(jn, lnp) + qw * oddr *  ddz * invvol(jn-1)
    basegrid%rhop(jnp,lnp) = basegrid%rhop(jnp,lnp) + qw *  ddr *  ddz * invvol(jnp-1)
    basegrid%rhominr = MIN(basegrid%rhominr,jn)
    basegrid%rhomaxr = MAX(basegrid%rhomaxr,jnp)
    basegrid%rhominz = MIN(basegrid%rhominz,ln)
    basegrid%rhomaxz = MAX(basegrid%rhomaxz,lnp)
#else
    basegrid%rho(jn, ln)  = basegrid%rho(jn, ln)  + qw * oddr * oddz * invvol(jn-1)
    basegrid%rho(jnp,ln)  = basegrid%rho(jnp,ln)  + qw *  ddr * oddz * invvol(jnp-1)
    basegrid%rho(jn, lnp) = basegrid%rho(jn, lnp) + qw * oddr *  ddz * invvol(jn-1)
    basegrid%rho(jnp,lnp) = basegrid%rho(jnp,lnp) + qw *  ddr *  ddz * invvol(jnp-1)
#endif
  end do
 END if
else ! IF(solvergeom==XZgeom) then
 IF(l4symtry .OR. l2symtry) then
   l_sym = .true.
 else
   l_sym = .false.
 END if
 IF(ngrids>1) then
!  call rhoweightrznew(xp,yp,zp,uzp,np,q)
  WRITE(0,*) 'mesh refinement not yet supported in XZ.'
  stop
 else
  invdr = 1._8/dr
  invdz = 1._8/dz

  invvolxz = 1._8 / (dr*dz)

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    IF(l_sym) then
      rpos = abs(xp(i))*invdr
    else
      rpos = (xp(i)-xmin)*invdr
    END if
    zpos = (zp(i)-basegrid%zminp)*invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    qw = q*w(i)
#ifdef MPIPARALLEL
    basegrid%rhop(jn, ln)  = basegrid%rhop(jn, ln)  + qw * oddr * oddz * invvolxz
    basegrid%rhop(jnp,ln)  = basegrid%rhop(jnp,ln)  + qw *  ddr * oddz * invvolxz
    basegrid%rhop(jn, lnp) = basegrid%rhop(jn, lnp) + qw * oddr *  ddz * invvolxz
    basegrid%rhop(jnp,lnp) = basegrid%rhop(jnp,lnp) + qw *  ddr *  ddz * invvolxz
    basegrid%rhominr = MIN(basegrid%rhominr,jn)
    basegrid%rhomaxr = MAX(basegrid%rhomaxr,jnp)
    basegrid%rhominz = MIN(basegrid%rhominz,ln)
    basegrid%rhomaxz = MAX(basegrid%rhomaxz,lnp)
#else
    basegrid%rho(jn, ln)  = basegrid%rho(jn, ln)  + qw * oddr * oddz * invvolxz
    basegrid%rho(jnp,ln)  = basegrid%rho(jnp,ln)  + qw *  ddr * oddz * invvolxz
    basegrid%rho(jn, lnp) = basegrid%rho(jn, lnp) + qw * oddr *  ddz * invvolxz
    basegrid%rho(jnp,lnp) = basegrid%rho(jnp,lnp) + qw *  ddr *  ddz * invvolxz
#endif
  end do
 END if
END if

if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
return
END SUBROUTINE RHOWEIGHTRZ_weights

subroutine rhoweightz(zp,uzp,np,q,nz,dz)
USE multigridrz
USE Subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nz
REAL(8), DIMENSION(np), INTENT(IN) :: zp, uzp
REAL(8), INTENT(IN) :: q, dz

REAL(8) :: zpos, ddz, oddz
INTEGER(ISZ) :: i, ln, lnp, igrid
REAL(8):: substarttime
LOGICAL(ISZ) :: ingrid

REAL(8), DIMENSION(:), ALLOCATABLE :: invdz, zmin

if (lw3dtimesubs) substarttime = wtime()

IF(np==0) return

ALLOCATE(invdz(ngrids),zmin(ngrids))

do igrid = 1, ngrids
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  zmin (igrid) = grids_ptr(igrid)%grid%zminp
end do

IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    igrid = 1
    grids=>basegrid
    ingrid=.false.
    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part(1,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part(1,ln)
        grids=>grids_ptr(igrid)%grid
        zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        ln = 1+INT(zpos)
      END if
    end do
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    lnp=ln+1
#ifdef MPIPARALLEL
    grids%rhop(1,ln)  = grids%rhop(1,ln)  + q * oddz * invdz(igrid)
    grids%rhop(1,lnp) = grids%rhop(1,lnp) + q * ddz  * invdz(igrid)
    grids%rhominz = MIN(grids%rhominz,ln)
    grids%rhomaxz = MAX(grids%rhomaxz,lnp)
#else
    grids%rho(1,ln)  = grids%rho(1,ln)  + q * oddz * invdz(igrid)
    grids%rho(1,lnp) = grids%rho(1,lnp) + q * ddz  * invdz(igrid)
#endif
  end do
else
  ! make charge deposition using CIC weighting
  igrid = 1
  do i = 1, np
    IF(uzp(i)==0.) cycle
    zpos = (zp(i)-zmin(igrid))*invdz(1)
    ln = 1+INT(zpos)
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    lnp=ln+1
#ifdef MPIPARALLEL
    basegrid%rhop(1,ln)  = basegrid%rhop(1,ln)  + q * oddz * invdz(1)
    basegrid%rhop(1,lnp) = basegrid%rhop(1,lnp) + q * ddz  * invdz(1)
    basegrid%rhominz = MIN(basegrid%rhominz,ln)
    basegrid%rhomaxz = MAX(basegrid%rhomaxz,lnp)
#else
    basegrid%rho(1,ln)  = basegrid%rho(1,ln)  + q * oddz * invdz(1)
    basegrid%rho(1,lnp) = basegrid%rho(1,lnp) + q * ddz  * invdz(1)
#endif
  end do
END if

DEALLOCATE(invdz,zmin)

if (lw3dtimesubs) timesetrho3d = timesetrho3d + wtime() - substarttime
return
END SUBROUTINE RHOWEIGHTZ

subroutine set_rho_rz(rho,nr,nz,id)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
REAL(8), DIMENSION(nr+1,nz+1), INTENT(IN) :: rho

  IF(id<1.or.id>ngrids) then
    WRITE(0,*) 'Error, id out of bounds'
    WRITE(0,*) 'given id = ',id,' while 1 < id < ',ngrids
    stop
  END if
  IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
    WRITE(0,*) 'Error, dimensions should be the same: '
    WRITE(0,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2)
    WRITE(0,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
    stop
  END if
  grids_ptr(id)%grid%rho=rho

return
end subroutine set_rho_rz

subroutine mix_rho_rz(rho,nr,nz,id,fmix)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
REAL(8), DIMENSION(nr+1,nz+1), INTENT(IN) :: rho
REAL(8), INTENT(IN) :: fmix

  IF(id<1.or.id>ngrids) then
    WRITE(0,*) 'Error, id out of bounds'
    WRITE(0,*) 'given id = ',id,' while 1 < id < ',ngrids
    stop
  END if
  IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
    WRITE(0,*) 'Error, dimensions should be the same: '
    WRITE(0,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2)
    WRITE(0,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
    stop
  END if
  grids_ptr(id)%grid%rho=(1.-fmix)*grids_ptr(id)%grid%rho + fmix*rho

return
end subroutine mix_rho_rz

subroutine get_rho_rz(rho,nr,nz,id,rhop)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz,rhop
REAL(8), DIMENSION(nr+1,nz+1), INTENT(IN OUT) :: rho

IF(solvergeom==XZgeom) call rhobndrz()

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

subroutine get_rho_z(rho,nz,id,rhop)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nz,rhop
REAL(8), DIMENSION(nz+1), INTENT(IN OUT) :: rho

  IF(id<1.or.id>ngrids) then
    WRITE(0,*) 'Error, id out of bounds'
    WRITE(0,*) 'given id = ',id,' while 1 < id < ',ngrids
    return
  END if
#ifdef MPIPARALLEL
  IF(rhop==1) then
    IF(SIZE(rho)/=SIZE(grids_ptr(id)%grid%rhop,2)) then
      WRITE(0,*) 'Error, dimensions should be the same: '
      WRITE(0,*) 'Nz for rhop    : ',SIZE(rho,1)
      WRITE(0,*) 'Nz for rhop(id): ',SIZE(grids_ptr(id)%grid%rhop,2)
      return
    END if
    rho=grids_ptr(id)%grid%rhop(1,:)
  else
#endif
    IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,2)) then
      WRITE(0,*) 'Error, dimensions should be the same: '
      WRITE(0,*) 'Nz for rho    : ',SIZE(rho,1)
      WRITE(0,*) 'Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,2)
      return
    END if
    rho=grids_ptr(id)%grid%rho(1,:)
#ifdef MPIPARALLEL
  END if
#endif
return
end subroutine get_rho_z

subroutine rhoweightrz_meshref(xp,yp,zp,uzp,np,q)
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
END subroutine rhoweightrz_meshref

subroutine rhoweightrz_meshref_weights(xp,yp,zp,uzp,wp,np,q)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp, wp
REAL(8), INTENT(IN) :: q

REAL(8) :: rpos, zpos, ddr, ddz, oddr, oddz, qw
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
    qw = q*wp(i)
#ifdef MPIPARALLEL
    grids%rhop(jn, ln)  = grids%rhop(jn, ln)  + qw * oddr * oddz * grids%invvol(jn)
    grids%rhop(jnp,ln)  = grids%rhop(jnp,ln)  + qw *  ddr * oddz * grids%invvol(jnp)
    grids%rhop(jn, lnp) = grids%rhop(jn, lnp) + qw * oddr *  ddz * grids%invvol(jn)
    grids%rhop(jnp,lnp) = grids%rhop(jnp,lnp) + qw *  ddr *  ddz * grids%invvol(jnp)
#else
    grids%rho(jn, ln)  = grids%rho(jn, ln)  + qw * oddr * oddz * grids%invvol(jn)
    grids%rho(jnp,ln)  = grids%rho(jnp,ln)  + qw *  ddr * oddz * grids%invvol(jnp)
    grids%rho(jn, lnp) = grids%rho(jn, lnp) + qw * oddr *  ddz * grids%invvol(jn)
    grids%rho(jnp,lnp) = grids%rho(jnp,lnp) + qw *  ddr *  ddz * grids%invvol(jnp)
#endif
  end do

  DEALLOCATE(invdr,invdz,zmin)

  return
END subroutine rhoweightrz_meshref_weights

subroutine reset_rzmgrid_rho()
USE multigridrz
implicit none
INTEGER(ISZ) :: ig

  IF(l_change_grid) then
    do ig = 1, ngrids_cg
      call del_subgrid(id_cg(ig,1))
      call add_grid(grids_ptr(id_cg(ig,2))%grid, &
                    nr_cg(ig), &
                    nz_cg(ig), &
                    dr_cg(ig), &
                    dz_cg(ig), &
                    rmin_cg(ig), &
                    zmin_cg(ig), &
                    guard_min_r_cg(ig), &
                    guard_max_r_cg(ig), &
                    guard_min_z_cg(ig), &
                    guard_max_z_cg(ig))
    END do
    ngrids_cg = 0
    l_change_grid = .false.
  end if

  do ig = 1, ngrids
#ifdef MPIPARALLEL
    grids_ptr(ig)%grid%rhop=0.
    grids_ptr(ig)%grid%rhominr = grids_ptr(ig)%grid%nr+2
    grids_ptr(ig)%grid%rhomaxr = -1
    grids_ptr(ig)%grid%rhominz = grids_ptr(ig)%grid%nz+2
    grids_ptr(ig)%grid%rhomaxz = -1
#else
    grids_ptr(ig)%grid%rho=0.
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

subroutine rhobndrz()
USE multigridrz
implicit none

INTEGER(ISZ) :: igrid, nr, nz

  do igrid = 1, ngrids
    nr = grids_ptr(igrid)%grid%nr
    nz = grids_ptr(igrid)%grid%nz
#ifdef MPIPARALLEL
    if(solvergeom/=RZgeom) then
      IF(grids_ptr(igrid)%grid%ixlbnd==neumann) grids_ptr(igrid)%grid%rhop(1,:)    = 2._8*grids_ptr(igrid)%grid%rhop(1,:)
    END if
    IF(grids_ptr(igrid)%grid%ixrbnd==neumann) grids_ptr(igrid)%grid%rhop(nr+1,:) = 2._8*grids_ptr(igrid)%grid%rhop(nr+1,:)
    IF(grids_ptr(igrid)%grid%izlbnd==neumann) grids_ptr(igrid)%grid%rhop(:,1)    = 2._8*grids_ptr(igrid)%grid%rhop(:,1)
    IF(grids_ptr(igrid)%grid%izrbnd==neumann) grids_ptr(igrid)%grid%rhop(:,nz+1) = 2._8*grids_ptr(igrid)%grid%rhop(:,nz+1)
#else
    if(solvergeom/=RZgeom) then
      IF(grids_ptr(igrid)%grid%ixlbnd==neumann) grids_ptr(igrid)%grid%rho(1,:)    = 2._8*grids_ptr(igrid)%grid%rho(1,:)
    END if
    IF(grids_ptr(igrid)%grid%ixrbnd==neumann) grids_ptr(igrid)%grid%rho(nr+1,:) = 2._8*grids_ptr(igrid)%grid%rho(nr+1,:)
    IF(grids_ptr(igrid)%grid%izlbnd==neumann) grids_ptr(igrid)%grid%rho(:,1)    = 2._8*grids_ptr(igrid)%grid%rho(:,1)
    IF(grids_ptr(igrid)%grid%izrbnd==neumann) grids_ptr(igrid)%grid%rho(:,nz+1) = 2._8*grids_ptr(igrid)%grid%rho(:,nz+1)
#endif
  end do

  IF(bound0==periodic) then
#ifdef MPIPARALLEL
      WRITE(0,*) 'ERROR:periodicity in RZ not yet supported on parallel platform, aborting.'
      stop
#endif
    IF(ngrids>1) then
      WRITE(0,*) 'ERROR:periodicity in RZ not yet supported with mesh refinement, aborting.'
      stop
    END if
    basegrid%rho(:,1) = basegrid%rho(:,1) + basegrid%rho(:,basegrid%nz+1)
    basegrid%rho(:,basegrid%nz+1) = basegrid%rho(:,1)
  END if

  return
end subroutine rhobndrz

 subroutine perphirz()
 USE multigridrz
 use Subtimers3d
 real(kind=8):: substarttime
 if (lw3dtimesubs) substarttime = wtime()

!  Sets the slices on the exterior of phi for periodicity
!  sets slice at -1 equal to the slice at nz-1
!  sets slice at nz+1 equal to the slice at 1

#ifdef MPIPARALLEL
   call perphi3d_slave(basegrid%phi(1:basegrid%nr+1,0:basegrid%nz+2),basegrid%nr,0,basegrid%nz)
#endif

  if (lw3dtimesubs) timeperphi3d = timeperphi3d + wtime() - substarttime
  return
end subroutine perphirz

subroutine perrhorz()
USE multigridrz
implicit none

#ifdef MPIPARALLEL
   call perrho3d_slave(basegrid%rho(0,0),basegrid%nr,0,basegrid%nz,periodic)
#endif

end subroutine perrhorz

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

subroutine getrhoforfieldsolvez(nz,rho)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: nz
REAL(8), DIMENSION(0:nz), INTENT(IN OUT) :: rho

  call get_rho_from_rhop(basegrid)

  return
end subroutine getrhoforfieldsolvez

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
USE subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, uzp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez

REAL(8) :: rpos, zpos, invrpos, ddr, ddz, oddr, oddz, er
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid

REAL(8):: substarttime
if (lw3dtimesubs) substarttime = wtime()

IF(ngrids>1 .and. .not.l_get_field_from_base) then

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

subroutine fieldweightxz(xp,zp,uzp,ex,ez,np)
USE multigridrz
USE subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, zp, uzp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ez

REAL(8) :: rpos, zpos, invrpos, ddr, ddz, oddr, oddz
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid, l_sym

REAL(8):: substarttime
if (lw3dtimesubs) substarttime = wtime()

IF(l2symtry .OR. l4symtry) then
  l_sym = .true.
else
  l_sym = .false.
END if

IF(ngrids>1 .and. .not.l_get_field_from_base) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    igrid = 1
    grids => basegrid
    ingrid=.false.
    if(l_sym) then
      rpos = (ABS(xp(i))-grids%xmin)*grids%invdr
    else
      rpos = (xp(i)-grids%xmin)*grids%invdr
    end if
    zpos = (zp(i)-grids%zminp)*grids%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part_field_dep(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part_field_dep(jn,ln)
        grids=>grids_ptr(igrid)%grid
        if(l_sym) then
          rpos = (ABS(xp(i))-grids%xmin)*grids%invdr
        else
          rpos = (xp(i)-grids%xmin)*grids%invdr
        end if
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
    ex(i) = 0.5*(oddr * oddz * (grids%phip(jn-1,ln  )-grids%phip(jn+1,ln  ))  &
            + ddr  * oddz * (grids%phip(jn  ,ln  )-grids%phip(jn+2,ln  ))  &
            + oddr * ddz  * (grids%phip(jn-1,ln+1)-grids%phip(jn+1,ln+1))  &
            + ddr  * ddz  * (grids%phip(jn  ,ln+1)-grids%phip(jn+2,ln+1)))*grids%invdr
#else
    ex(i) = 0.5*(oddr * oddz * (grids%phi(jn-1,ln  )-grids%phi(jn+1,ln  ))  &
            + ddr  * oddz * (grids%phi(jn  ,ln  )-grids%phi(jn+2,ln  ))  &
            + oddr * ddz  * (grids%phi(jn-1,ln+1)-grids%phi(jn+1,ln+1))  &
            + ddr  * ddz  * (grids%phi(jn  ,ln+1)-grids%phi(jn+2,ln+1)))*grids%invdr
#endif
    IF(l_sym) then
      IF(xp(i)<0.) ex(i) = -ex(i)
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
    if(l_sym) then
      rpos = ABS(xp(i))*basegrid%invdr
    else
      rpos = (xp(i)-basegrid%xmin)*basegrid%invdr
    end if
    zpos = (zp(i)-basegrid%zminp)*basegrid%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    ex(i) = 0.5*(oddr * oddz * (basegrid%phip(jn-1,ln  )-basegrid%phip(jn+1,ln  ))  &
            + ddr  * oddz * (basegrid%phip(jn  ,ln  )-basegrid%phip(jn+2,ln  ))  &
            + oddr * ddz  * (basegrid%phip(jn-1,ln+1)-basegrid%phip(jn+1,ln+1))  &
            + ddr  * ddz  * (basegrid%phip(jn  ,ln+1)-basegrid%phip(jn+2,ln+1)))*basegrid%invdr
#else
    ex(i) = 0.5*(oddr * oddz * (basegrid%phi(jn-1,ln  )-basegrid%phi(jn+1,ln  ))  &
            + ddr  * oddz * (basegrid%phi(jn  ,ln  )-basegrid%phi(jn+2,ln  ))  &
            + oddr * ddz  * (basegrid%phi(jn-1,ln+1)-basegrid%phi(jn+1,ln+1))  &
            + ddr  * ddz  * (basegrid%phi(jn  ,ln+1)-basegrid%phi(jn+2,ln+1)))*basegrid%invdr
#endif
    IF(l_sym) then
      IF(xp(i)<0.) ex(i) = -ex(i)
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
end subroutine fieldweightxz

subroutine fieldweightz(zp,uzp,ez,np)
USE multigridrz
USE subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: zp, uzp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ez

REAL(8) :: zpos, ddz, oddz
INTEGER(ISZ) :: i, l, ln, lnp, igrid
LOGICAL(ISZ) :: ingrid

REAL(8):: substarttime
if (lw3dtimesubs) substarttime = wtime()

IF(ngrids>1 .and. .not.l_get_field_from_base) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    igrid = 1
    grids => basegrid
    ingrid=.false.
    zpos = (zp(i)-grids%zminp)*grids%invdz
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part_field_dep(1,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part_field_dep(1,ln)
        grids=>grids_ptr(igrid)%grid
        zpos = (zp(i)-grids%zminp)*grids%invdz
        ln = 1+INT(zpos)
      END if
    end do
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddz * (grids%phip(1,ln-1)-grids%phip(1,ln+1))  &
               + ddz  * (grids%phip(1,ln  )-grids%phip(1,ln+2)))*grids%invdz
#else
    ez(i) = 0.5*(oddz * (grids%phi(1,ln-1)-grids%phi(1,ln+1))  &
               + ddz  * (grids%phi(1,ln  )-grids%phi(1,ln+2)))*grids%invdz
#endif
  END do
else
  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(uzp(i)==0.) cycle
    ingrid=.false.
    zpos = (zp(i)-basegrid%zminp)*basegrid%invdz
    ln = 1+INT(zpos)
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddz * (basegrid%phip(1,ln-1)-basegrid%phip(1,ln+1))  &
               + ddz  * (basegrid%phip(1,ln  )-basegrid%phip(1,ln+2)))*basegrid%invdz
#else
    ez(i) = 0.5*(oddz * (basegrid%phi(1,ln-1)-basegrid%phi(1,ln+1))  &
               + ddz  * (basegrid%phi(1,ln  )-basegrid%phi(1,ln+2)))*basegrid%invdz
#endif
  END do
END if

  if (lw3dtimesubs) timesete3d = timesete3d + wtime() - substarttime
  return
end subroutine fieldweightz

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

IF(ngrids>1 .and. .not.l_get_injphi_from_base) then

  do i = 1, np
    igrid = 1
    grids => basegrid
    IF(zp(i)<grids%zmin.or.zp(i)>=grids%zmax) cycle
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
        IF(zp(i)<grids%zmin.or.zp(i)>=grids%zmax) cycle
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
    IF(zp(i)<basegrid%zmin.or.zp(i)>=basegrid%zmax) cycle
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

subroutine setphiz(np,zp,p)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: p

REAL(8) :: zpos, ddz, oddz
INTEGER(ISZ) :: i, l, ln, lnp, igrid
LOGICAL(ISZ) :: ingrid

! Collect phi using linear interpolation

IF(ngrids>1 .and. .not.l_get_injphi_from_base) then
  do i = 1, np
    igrid = 1
    grids => basegrid
    ingrid=.false.
    zpos = (zp(i)-grids%zmin)*grids%invdz
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(grids%loc_part(1,ln)==igrid) then
        ingrid=.true.
      else
        igrid = grids%loc_part(1,ln)
        grids=>grids_ptr(igrid)%grid
        zpos = (zp(i)-grids%zminp)*grids%invdz
        ln = 1+INT(zpos)
      END if
    end do
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    p(i) = oddz * grids%phi(1,ln  )  &
         + ddz  * grids%phi(1,ln+1)
  END do
else
  do i = 1, np
    zpos = (zp(i)-basegrid%zmin)
    IF(zpos<basegrid%zmin.or.zpos>=basegrid%zmax) cycle
    zpos = zpos*basegrid%invdz
    ln = 1+INT(zpos)
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    p(i) = oddz * basegrid%phi(1,ln  )  &
         + ddz  * basegrid%phi(1,ln+1)
  END do
END if

  return
end subroutine setphiz

!=============================================================================
subroutine setbnd_subgrid_to_inj_d()
use multigridrz
use InjectVars
use InjectVars3d

INTEGER(ISZ) :: ij, j, l, igrid, il
REAL(8) :: rs, rc, r, z
TYPE(grdptr), POINTER :: g

INTEGER(ISZ) :: nconds, n
INTEGER(ISZ), DIMENSION(:), ALLOCATABLE :: ixcond, izcond

 do ij = 1, ninject
   rs = ainject(ij)
   rc = rinject(ij)
   do igrid=1,ngrids
     g => grids_ptr(igrid)%grid
     do j = 1, g%nr+1
       r = (j-1)*g%dr
       IF(r>rs) exit
       ! we assume source emits forward
       z = zinject(ij) + rc - SQRT(rc**2-r**2) + inj_d(ij)*g%dz
       l = MAX(1,2+INT((z-g%zmin)/g%dz))
       g%loc_part_field_dep(j,l:) = g%id
     end do
   end do
 end do

 do ij = 1, ninject
   rs = ainject(ij)
   rc = rinject(ij)
! we begin at 2 because 1 is supposed to be associated with the basegrid
   do igrid=2,ngrids
     g => grids_ptr(igrid)%grid
     do il = 1, g%nlevels
       nconds = 0
       n = (g%bnd(il)%nr+1)*(g%bnd(il)%nz+1)
       ALLOCATE(ixcond(n),izcond(n))
       do j = 1, g%bnd(il)%nr+1
         r = MIN(rs,(j-1)*g%bnd(il)%dr)
         ! we assume source emits forward
         z = zinject(ij) + rc - SQRT(rc**2-r**2) + inj_d(ij)*g%bnd(il)%dz
         l = 4+MAX(1,2+INT((z-g%zmin)/g%bnd(il)%dz))
         do l = 4+MAX(1,2+INT((z-g%zmin)/g%bnd(il)%dz)), g%bnd(il)%nz+1
           IF(g%bnd(il)%v(j,l)==v_vacuum) then
             g%bnd(il)%v(j,l)=v_cond
             nconds=nconds+1
             ixcond(nconds)=j
             izcond(nconds)=l
           END if
         END do
       end do
       IF(nconds>0) then
         call init_bnd_sublevel(g%bnd(il),0,nconds)
         g%bnd(il)%cnd%jcond(1:nconds) = ixcond(1:nconds)
         g%bnd(il)%cnd%kcond(1:nconds) = izcond(1:nconds)
         g%bnd(il)%cnd%voltage(1:nconds) = 0.
         g%bnd(il)%cnd%condid(1:nconds) = -1
        END if
       DEALLOCATE(ixcond,izcond)
     end do
   end do
 end do

return
end subroutine setbnd_subgrid_to_inj_d

!=============================================================================
subroutine set_patches_around_emitter(id,np,ij,nzi,guard_min_r,guard_max_r,guard_min_z,guard_max_z)
use multigridrz
use InjectVars
use InjectVars3d

INTEGER(ISZ), INTENT(IN) :: id, & ! id of grid on which to add the patches
                            np, & ! number of patches
                            ij, & ! id of injection source
                            nzi, & ! size of patches in z (number of meshes)
                            guard_min_r, & ! number of guard cells at lower end in r for field gathering
                            guard_max_r, & ! number of guard cells at upper end in r for field gathering
                            guard_min_z, & ! number of guard cells at lower end in z for field gathering
                            guard_max_z    ! number of guard cells at upper end in z for field gathering

INTEGER(ISZ) :: i, j, l, igrid, il, nr, nz, l0
REAL(8) :: rs, rc, r, z, dr, dz, rmin, zmin
TYPE(grdptr), POINTER :: g
INTEGER(ISZ), DIMENSION(:), ALLOCATABLE :: lshifts

  nz = nzi

  IF(nz<mgridrz_nmeshmin) then
    nz = mgridrz_nmeshmin
  END if

  rs = ainject(ij)
  rc = rinject(ij)
  g => grids_ptr(id)%grid
  dr = g%dr
  dz = g%dz

  do i = 1, np
    il = g%nlevels
    ALLOCATE(lshifts(2*g%bnd(il)%nr+1))
    l0 = 1+INT((zinject(ij)-g%bnd(il)%zmin)/g%bnd(il)%dz)!-g%bnd(il)%lshift(1)
    do j = 1, g%bnd(il)%nr+1
      r = (j-1)*g%bnd(il)%dr
      IF(r>rs+(2+guard_max_r)*g%bnd(il)%dr) exit
      ! we assume source emits forward
      z = zinject(ij) + rc - SQRT(rc**2-r**2)
      l = 1+INT((z-g%bnd(il)%zmin)/g%bnd(il)%dz)!-g%bnd(il)%lshift(j)
      lshifts(1+(j-1)*2:2+(j-1)*2) = 2*(l-l0)
      IF(lshifts(1+(j-1)*2)<0) WRITE(0,*) 'error in set_patches_around_emitter: shifts < 0, shifts = ',lshifts(1+(j-1)*2)
    END do
    nr = 2*(j-1)
    lshifts(nr+1) = lshifts(nr)
    dr = 0.5*dr
    dz = 0.5*dz
    rmin = 0.
    zmin = g%bnd(il)%zmin+INT((zinject(ij)-g%bnd(il)%zmin)/g%bnd(il)%dz)*g%bnd(il)%dz
    WRITE(0,*) 'call add_grid'
    WRITE(0,*) 'nr = ',nr
    WRITE(0,*) 'nz = ',nz
    WRITE(0,*) 'dr = ',dr
    WRITE(0,*) 'dz = ',dz
    WRITE(0,*) 'zmin = ',zmin
    WRITE(0,*) 'lshifts = ',lshifts(1:nr+1)
    call add_grid(g,nr,nz,dr,dz,rmin,zmin,guard_min_r,guard_max_r,guard_min_z,guard_max_z,lshifts(1:nr+1))
    DEALLOCATE(lshifts)
    g => grids_ptr(ngrids)%grid
  end do

return
end subroutine set_patches_around_emitter

subroutine clean_conductor_interior()
USE multigridrz
implicit none

INTEGER(ISZ) :: i, ic, il, igrid, ix, iz, ncond
TYPE(grdptr), POINTER :: g

  do igrid=1,ngrids
    g => grids_ptr(igrid)%grid
    do il = 1, g%nlevels
      do ic = 1, g%bnd(il)%nb_conductors
        IF(ic==1) then
          g%bnd(il)%cnd => g%bnd(il)%first
        else
          g%bnd(il)%cnd => g%bnd(il)%cnd%next
        END if
        ncond = 0
        do i = 1, g%bnd(il)%cnd%ncond
          ix = g%bnd(il)%cnd%jcond(i)
          iz = g%bnd(il)%cnd%kcond(i)
          IF(.NOT.(ABS(g%bnd(il)%v(ix+1,iz+1))==v_cond .AND. &
                   ABS(g%bnd(il)%v(ix-1,iz+1))==v_cond .AND. &
                   ABS(g%bnd(il)%v(ix+1,iz-1))==v_cond .AND. &
                   ABS(g%bnd(il)%v(ix-1,iz-1))==v_cond)) then
            ncond = ncond + 1
            IF(i /= ncond) then
              g%bnd(il)%cnd%jcond(ncond)   = g%bnd(il)%cnd%jcond(i)
              g%bnd(il)%cnd%kcond(ncond)   = g%bnd(il)%cnd%kcond(i)
              g%bnd(il)%cnd%voltage(ncond) = g%bnd(il)%cnd%voltage(i)
              g%bnd(il)%cnd%condid(ncond)  = g%bnd(il)%cnd%condid(i)
            END if
          else
            g%bnd(il)%v(ix,iz) = -v_cond
          END if
        END do
        g%bnd(il)%cnd%ncond = ncond
      END do
    END do
  END do

  return
end subroutine clean_conductor_interior

subroutine build_vlocs()
USE multigridrz
implicit none

INTEGER(ISZ) :: il, igrid, j, l, ilred, ilblack, jmin, jmax, lmin, lmax
TYPE(grdptr), POINTER :: g

  vlocs = .true.

  do igrid=1,ngrids
    g => grids_ptr(igrid)%grid
    g%npmin = 1
    do il = 1, g%nlevels
      g%bnd(il)%nvlocs = 0
      g%bnd(il)%nvlocsred = 0
      do l = 1, g%bnd(il)%nz+1
        do j = 1, g%bnd(il)%nr+1
          IF(g%bnd(il)%v(j,l)==v_vacuum) then
            IF(MOD(j+l,2)==0) g%bnd(il)%nvlocsred = g%bnd(il)%nvlocsred + 1
            g%bnd(il)%nvlocs = g%bnd(il)%nvlocs + 1
          END if
        end do
      end do
      ALLOCATE(g%bnd(il)%vlocs_j(g%bnd(il)%nvlocs),g%bnd(il)%vlocs_k(g%bnd(il)%nvlocs))
      ilred = 0
      ilblack = g%bnd(il)%nvlocsred
      jmax = 1
      jmin = g%bnd(il)%nr+1
      lmax = 1
      lmin = g%bnd(il)%nz+1
      do l = 1, g%bnd(il)%nz+1
        do j = 1, g%bnd(il)%nr+1
          IF(g%bnd(il)%v(j,l)==v_vacuum) then
            jmin = MIN(jmin,j)
            jmax = MAX(jmax,j)
            lmin = MIN(lmin,l)
            lmax = MAX(lmax,l)
            IF(MOD(j+l,2)==0) then
              ilred = ilred + 1
              g%bnd(il)%vlocs_j(ilred) = j
              g%bnd(il)%vlocs_k(ilred) = l
            else
              ilblack = ilblack + 1
              g%bnd(il)%vlocs_j(ilblack) = j
              g%bnd(il)%vlocs_k(ilblack) = l
            END if
          END if
        end do
      end do
      IF(il<g%nlevels-2) then
        IF((jmax-jmin)<mgridrz_nmeshmin .OR. (lmax-lmin)<mgridrz_nmeshmin) g%npmin = g%npmin+1
      END if
    END do
  END do

  return
end subroutine build_vlocs

subroutine init_base(nr,nz,dr,dz,rmin,zmin)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: nr,nz
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin

 call init_basegrid(nr,nz,dr,dz,rmin,zmin)

 return
END subroutine init_base

subroutine del_base()
USE multigridrz

  call del_basegrid()

 return
END subroutine del_base

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

subroutine del_subgrid(id)
USE multigridrz
implicit none
INTEGER, INTENT(IN) :: id

  IF(id<1 .or. id>ngrids) then
    WRITE(0,*) 'Fatal error in add_subgrid: id = ', id ,' WHILE id = (1,...,',ngrids,')'
    stop
  END if

  call del_grid(grids_ptr(id)%grid)

  return
END subroutine del_subgrid

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

subroutine set_basegrid_phi()
USE multigridrz
USE Fields3D
implicit none

  basegrid%phi(1:basegrid%nr+1,:) = phi(:,0,:)

return
END subroutine set_basegrid_phi

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
      case ("v")
        array = grids_ptr(id)%grid%bnd(grids_ptr(id)%grid%nlevels)%v(1:nr+1,1:nz+1)
      case ("loc_part","lp")
        array = grids_ptr(id)%grid%loc_part(1:nr+1,1:nz+1)
      case ("loc_part_field_dep","lpfd")
        array = grids_ptr(id)%grid%loc_part_field_dep(1:nr+1,1:nz+1)
      case ("lshifts","ls")
        array(:,1) = grids_ptr(id)%grid%zmin &
                   + grids_ptr(id)%grid%bnd(grids_ptr(id)%grid%nlevels)%lshift(1:nr+1) &
                   * grids_ptr(id)%grid%dz
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
                 //'  res [no abrv]' &
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

 IF(solvergeom==Zgeom) return

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
     icondlevel(icc) = ilevel - 1
     condvolt(icc) = bnd%cnd%voltage(i)
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
     iecndlevel(ice) = ilevel - 1
     ecvolt(ice) = bnd%cnd%volt0xm(i)
     ecvoltmx(ice) = bnd%cnd%volt0xm(i)
     ecvoltpx(ice) = bnd%cnd%volt0xp(i)
     ecvoltmz(ice) = bnd%cnd%volt0zm(i)
     ecvoltpz(ice) = bnd%cnd%volt0zp(i)
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
     iocndlevel(ico) = ilevel - 1
     ocvolt(ico) = bnd%cnd%volt0xm(i)
     ocvoltmx(ico) = bnd%cnd%volt0xm(i)
     ocvoltpx(ico) = bnd%cnd%volt0xp(i)
     ocvoltmz(ico) = bnd%cnd%volt0zm(i)
     ocvoltpz(ico) = bnd%cnd%volt0zp(i)
    END if
   end do
 END do
 necndbdy = ice
 nocndbdy = ico

return
end subroutine get_cond_rz

subroutine setconductorvoltagerz(volt,nz,zmmin,dz,discrete)
USE multigridrz
implicit none
integer(ISZ):: nz
real(kind=8):: volt(0:nz)
real(kind=8):: zmmin,dz
logical(ISZ):: discrete

INTEGER :: igrid,i,iv,ic,icc,ice,ico
integer(ISZ):: iz
real(kind=8):: zz,wz,vv
TYPE(conductor_type), POINTER :: cndpnt
real(kind=8):: dxm,dxp,dzm,dzp,dxx,dzz,r,rm,rp

do igrid=1,ngrids
  nlevels=grids_ptr(igrid)%grid%nlevels
  bndy => grids_ptr(igrid)%grid%bnd
  do i = nlevels,1,-1

   do iv=1, bndy(i)%nb_conductors
     IF(iv==1) then
       cndpnt => bndy(i)%first
     else
       cndpnt => cndpnt%next
     END if

    do ic=1,cndpnt%ncond
      zz = grids_ptr(igrid)%grid%zmin + bndy(i)%dz*(cndpnt%kcond(ic)-1)
      if (zmmin <= zz .and. zz < zmmin + nz*dz) then
        iz = int(zz/dz)
        wz =     zz/dz - iz
        cndpnt%voltage(ic) = volt(iz)*(1.-wz) + volt(iz+1)*wz
      else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + bndy(i)%dz) then
        cndpnt%voltage(ic) = volt(nz)
      endif
    enddo

    do ic = 1,cndpnt%nbbnd
      zz = grids_ptr(igrid)%grid%zmin + bndy(i)%dz*(cndpnt%kk(ic)-1)
      if (zmmin <= zz .and. zz < zmmin + nz*dz) then
        iz = int(zz/dz)
        wz =     zz/dz - iz
        vv = volt(iz)*(1.-wz) + volt(iz+1)*wz
        if (cndpnt%dxm(ic) < bndy(i)%dr) cndpnt%volt0xm(ic) = vv
        if (cndpnt%dxp(ic) < bndy(i)%dr) cndpnt%volt0xp(ic) = vv
      else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + bndy(i)%dz) then
        vv = volt(nz)
        if (cndpnt%dxm(ic) < bndy(i)%dr) cndpnt%volt0xm(ic) = vv
        if (cndpnt%dxp(ic) < bndy(i)%dr) cndpnt%volt0xp(ic) = vv
      endif
      if (cndpnt%dzm(ic) < bndy(i)%dz) then
        zz = grids_ptr(igrid)%grid%zmin + bndy(i)%dz*(cndpnt%kk(ic)-1) &
             - cndpnt%dzm(ic)
        if (zmmin <= zz .and. zz < zmmin + nz*dz) then
          iz = int(zz/dz)
          wz =     zz/dz - iz
          if (discrete) wz = 0.
          cndpnt%volt0zm(ic) = volt(iz)*(1.-wz) + volt(iz+1)*wz
        else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + bndy(i)%dz) then
          cndpnt%volt0zm(ic) = volt(nz)
        endif
      endif
      if (cndpnt%dzp(ic) < bndy(i)%dz) then
        zz = grids_ptr(igrid)%grid%zmin + bndy(i)%dz*(cndpnt%kk(ic)-1) &
             + cndpnt%dzp(ic)
        if (zmmin <= zz .and. zz < zmmin + nz*dz) then
          iz = int(zz/dz)
          wz =     zz/dz - iz
          if (discrete) wz = 1.
          cndpnt%volt0zp(ic) = volt(iz)*(1.-wz) + volt(iz+1)*wz
        else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + bndy(i)%dz) then
          cndpnt%volt0zp(ic) = volt(nz)
        endif
      endif

      dxm = cndpnt%dxm(ic)
      dxp = cndpnt%dxp(ic)
      dzm = cndpnt%dzm(ic)
      dzp = cndpnt%dzp(ic)
      select case (bnd_method)
        case (egun)
          dxx=bndy(i)%dr
          dzz=bndy(i)%dz
        case (ecb)
          dxx=0.5_8*(dxp+dxm)  !ecb
          dzz=0.5_8*(dzp+dzm)  !ecb
        case default
      end select
      IF(solvergeom==RZgeom) then
       IF(cndpnt%jj(ic)==1) then
        cndpnt%cfxp(ic) = 4._8/(dxp*dxx)
       else
        r = (cndpnt%jj(ic)-1)*bndy(i)%dr
        select case (bnd_method)
          case (egun)
            rm = r-0.5_8*bndy(i)%dr
            rp = r+0.5_8*bndy(i)%dr
          case (ecb)
            rm = r-0.5_8*dxm
            rp = r+0.5_8*dxp
          case default
        end select
        cndpnt%cfxm(ic) = rm/(r*dxm*dxx)
        cndpnt%cfxp(ic) = rp/(r*dxp*dxx)
       END if
      else !IF(solvergeom==XZgeom) then
        cndpnt%cfxm(ic) = 1._8/(dxm*dxx)
        cndpnt%cfxp(ic) = 1._8/(dxp*dxx)
      END if
      cndpnt%cfzm(ic) = 1._8/(dzm*dzz)
      cndpnt%cfzp(ic) = 1._8/(dzp*dzz)
      IF(dxm>=bndy(i)%dr) then
        cndpnt%phi0xm(ic)=0._8
      else
        cndpnt%phi0xm(ic)=cndpnt%cfxm(ic)*cndpnt%volt0xm(ic)
        cndpnt%cfxm(ic)=0._8
      END if
      IF(dxp>=bndy(i)%dr) then
        cndpnt%phi0xp(ic)=0._8
      else
        cndpnt%phi0xp(ic)=cndpnt%cfxp(ic)*cndpnt%volt0xp(ic)
        cndpnt%cfxp(ic)=0._8
      END if
      IF(dzm>=bndy(i)%dz) then
        cndpnt%phi0zm(ic)=0._8
      else
        cndpnt%phi0zm(ic)=cndpnt%cfzm(ic)*cndpnt%volt0zm(ic)
        cndpnt%cfzm(ic)=0._8
      END if
      IF(dzp>=bndy(i)%dz) then
        cndpnt%phi0zp(ic)=0._8
      else
        cndpnt%phi0zp(ic)=cndpnt%cfzp(ic)*cndpnt%volt0zp(ic)
        cndpnt%cfzp(ic)=0._8
      END if

    enddo
   enddo
  enddo
enddo
return
end subroutine setconductorvoltagerz

subroutine setconductorvoltagerz_id(id,volt)
USE multigridrz
implicit none
integer(ISZ):: id
real(kind=8):: volt

INTEGER :: igrid,i,iv,ic,icc,ice,ico
integer(ISZ):: iz
real(kind=8):: zz,wz,vv
TYPE(conductor_type), POINTER :: cndpnt
real(kind=8):: dxm,dxp,dzm,dzp,dxx,dzz,r,rm,rp
LOGICAL(ISZ) :: l_change

do igrid=1,ngrids
  nlevels=grids_ptr(igrid)%grid%nlevels
  bndy => grids_ptr(igrid)%grid%bnd
  do i = nlevels,1,-1

   do iv=1, bndy(i)%nb_conductors
     IF(iv==1) then
       cndpnt => bndy(i)%first
     else
       cndpnt => cndpnt%next
     END if

    do ic=1,cndpnt%ncond
      IF(cndpnt%condid(ic)==id) cndpnt%voltage(ic) = volt
    enddo

    do ic = 1,cndpnt%nbbnd
     l_change = .false.
     if (cndpnt%dxm(ic) < bndy(i)%dr .and. cndpnt%condidxm(ic)==id) then
      l_change = .true.
      cndpnt%volt0xm(ic) = volt
     END if
     if (cndpnt%dxp(ic) < bndy(i)%dr .and. cndpnt%condidxp(ic)==id) then
      l_change = .true.
      cndpnt%volt0xp(ic) = volt
     END if
     if (cndpnt%dzm(ic) < bndy(i)%dz .and. cndpnt%condidzm(ic)==id) then
      l_change = .true.
      cndpnt%volt0zm(ic) = volt
     END if
     if (cndpnt%dzp(ic) < bndy(i)%dz .and. cndpnt%condidzp(ic)==id) then
      l_change = .true.
      cndpnt%volt0zp(ic) = volt
     END if

     IF(l_change) then
      dxm = cndpnt%dxm(ic)
      dxp = cndpnt%dxp(ic)
      dzm = cndpnt%dzm(ic)
      dzp = cndpnt%dzp(ic)
      select case (bnd_method)
        case (egun)
          dxx=bndy(i)%dr
          dzz=bndy(i)%dz
        case (ecb)
          dxx=0.5_8*(dxp+dxm)  !ecb
          dzz=0.5_8*(dzp+dzm)  !ecb
        case default
      end select
      IF(solvergeom==RZgeom) then
       IF(cndpnt%jj(ic)==1) then
        cndpnt%cfxp(ic) = 4._8/(dxp*dxx)
       else
        r = (cndpnt%jj(ic)-1)*bndy(i)%dr
        select case (bnd_method)
          case (egun)
            rm = r-0.5_8*bndy(i)%dr
            rp = r+0.5_8*bndy(i)%dr
          case (ecb)
            rm = r-0.5_8*dxm
            rp = r+0.5_8*dxp
          case default
        end select
        cndpnt%cfxm(ic) = rm/(r*dxm*dxx)
        cndpnt%cfxp(ic) = rp/(r*dxp*dxx)
       END if
      else !IF(solvergeom==XZgeom) then
        cndpnt%cfxm(ic) = 1._8/(dxm*dxx)
        cndpnt%cfxp(ic) = 1._8/(dxp*dxx)
      END if
      cndpnt%cfzm(ic) = 1._8/(dzm*dzz)
      cndpnt%cfzp(ic) = 1._8/(dzp*dzz)
      IF(dxm>=bndy(i)%dr) then
        cndpnt%phi0xm(ic)=0._8
      else
        cndpnt%phi0xm(ic)=cndpnt%cfxm(ic)*cndpnt%volt0xm(ic)
        cndpnt%cfxm(ic)=0._8
      END if
      IF(dxp>=bndy(i)%dr) then
        cndpnt%phi0xp(ic)=0._8
      else
        cndpnt%phi0xp(ic)=cndpnt%cfxp(ic)*cndpnt%volt0xp(ic)
        cndpnt%cfxp(ic)=0._8
      END if
      IF(dzm>=bndy(i)%dz) then
        cndpnt%phi0zm(ic)=0._8
      else
        cndpnt%phi0zm(ic)=cndpnt%cfzm(ic)*cndpnt%volt0zm(ic)
        cndpnt%cfzm(ic)=0._8
      END if
      IF(dzp>=bndy(i)%dz) then
        cndpnt%phi0zp(ic)=0._8
      else
        cndpnt%phi0zp(ic)=cndpnt%cfzp(ic)*cndpnt%volt0zp(ic)
        cndpnt%cfzp(ic)=0._8
      END if
     END if

    enddo
   enddo
  enddo
enddo
return
end subroutine setconductorvoltagerz_id

subroutine init_gridinit()
USE multigridrz
implicit none
INTEGER(ISZ) :: i

  ALLOCATE(gridinit(ngrids))
  do i = 1, ngrids
    ALLOCATE(gridinit(i)%grid)
    ALLOCATE(gridinit(i)%grid%phi(LBOUND(grids_ptr(i)%grid%phi,1):UBOUND(grids_ptr(i)%grid%phi,1), &
                                  LBOUND(grids_ptr(i)%grid%phi,2):UBOUND(grids_ptr(i)%grid%phi,2)))
    gridinit(i)%grid%phi = grids_ptr(i)%grid%phi
  end do

return
end subroutine init_gridinit

subroutine change_loc_part()
USE multigridrz
implicit none
!TYPE(grdptr), pointer :: g

INTEGER(ISZ):: i

  IF(.not. l_change_loc_part) return

  do i = 1, nz_rmc+1
    basegrid%loc_part_field_dep(1:rmc(i)-1,i) = 1
    basegrid%loc_part_field_dep(rmc(i):,i) = 2
  end do
  l_change_loc_part = .false.

END subroutine change_loc_part

