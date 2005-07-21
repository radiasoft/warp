!     Last change:  JLV   7 Jul 2004    3:28 pm
#include "top.h"

module multigrid_common

USE Constant
USE multigrid_common_base
USE GridBoundary3d
USE InGen3d, ONLY:solvergeom,RZgeom,XYZgeom,XZgeom,XYgeom,Zgeom,Rgeom,l2symtry,l4symtry
#ifdef MPIPARALLEL
  use Parallel
  use mpirz
#endif

REAL(8), EXTERNAL :: wtime

INTEGER(ISZ), parameter :: dirichlet=0,neumann=1,periodic=2,patchbnd=3,othertask=-1  ! boundary condition types
INTEGER(ISZ), parameter :: v_vacuum=0, v_cond=1, v_bnd=2, v_dirichlet=3
INTEGER(ISZ), parameter :: egun=0, ecb=1

#ifdef MPIPARALLEL
  INTEGER(ISZ) :: nworkpproc, workfact=8
  logical(ISZ) :: l_mpi_barriers=.true. 
#else
  INTEGER(ISZ) :: my_index = 0
#endif

end module multigrid_common

module multigridrz
! module containing RZ multigrid solver

USE multigrid_common
USE FRZmgrid
USE CONDtypemodule
USE BNDtypemodule
USE GRIDtypemodule
USE GRDPTRtypemodule
USE FRZmgrid_ptrs

TYPE(BNDtype), pointer :: bndcurrent

TYPE(GRDPTRtype), DIMENSION(:), ALLOCATABLE :: grids_ptr, gridinit

contains

subroutine init_basegrid(nr,nz,dr,dz,rmin,zmin)
!USE Multigrid3d
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin
INTEGER(ISZ) :: i,j, nzp
TYPE(GRIDtype), POINTER :: bg
TYPE(BNDtype), POINTER :: b

  if (lverbose>=1) then
    write(o_line,'("Init basegrid")')
    call remark(trim(o_line))
  endif

  IF(.not. associated(basegrid)) call set_basegrid()
  bg => basegrid

  inveps0 = 1./eps0

  IF(solvergeom==Zgeom) then
    nguardx = 0
  END if
  IF(solvergeom==Rgeom) then
    nguardz = 0
  END if

!  call gallot("FRZmgrid_ptrs",0)

#ifdef MPIPARALLEL
  workfact = mgridrz_workfact
  bg%nzp   = nzpslave(my_index)
  bg%nrpar = nr
  bg%nzpar = bg%nzp
#else
  bg%nzp   = nz
  bg%nrpar = 0
  bg%nzpar = 0
#endif
  grids_nids=1
  bg%nr=nr
  bg%dr=dr
  bg%rmin=rmin
  bg%rmax=rmin+nr*dr
  bg%xmin=rmin
  bg%xmax=rmin+nr*dr
  bg%nz=nz
#ifdef MPIPARALLEL
!  bg%zminp=zpslmin(my_index)
  bg%zminp=zpslmin(0)+izpslave(my_index)*dz
#else
  bg%zminp=zmin
#endif
  bg%dz=dz
  bg%zmin=zmin
  bg%zmax=zmin+nz*dz
  bg%jmin=1
  bg%jmax=nr+1
  bg%lmin=1
  bg%lmax=nr+1
  bg%nguardx = nguardx
  bg%nguardz = nguardz
  call GRIDtypeallot(bg)
  bg%gid=1
  bg%loc_part=1
  bg%loc_part_fd=1
  bg%mgparam = mgridrz_mgparam
  bg%npre = mgridrz_npre
  bg%npost = mgridrz_npost
  bg%ncycles = mgridrz_ncycles
  bg%ncmax = mgridrz_ncmax
  bg%npmin = mgridrz_levels_min
  bg%phi=0.
  bg%rho=0.
#ifdef MPIPARALLEL
  bg%phip=0.
  bg%rhop=0.
#endif
  bg%transit_min_r = 0
  bg%transit_max_r = 0
  bg%transit_min_z = 0
  bg%transit_max_z = 0
  ngrids=1
  level_del_grid=0
  n_avail_ids=0
  avail_ids=-1
  bg%invdr = 1._8/dr
  bg%invdz = 1._8/dz
  IF(solvergeom==RZgeom .or. solvergeom==Rgeom) then
    ! computes divider by cell volumes to get density
    IF(bg%rmin==0.) then
      j = 1
      ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
      ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
      ! and Verboncoeur, J. of Comp. Phys.,
      bg%invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr * dz))
      do j = 2, nr+1
        bg%invvol(j) = 1._8 / (2._8 * pi * real(j-1,8) * dr * dr * dz)
      end do
    else
      do j = 1, nr+1
        bg%invvol(j) = 1._8 / (2._8 * pi * (bg%rmin+real(j-1,8)*dr) * dr * dz)
      end do
    END if
    IF(solvergeom==Rgeom) bg%invvol = bg%invvol * dz
  else ! solvergeom==XZgeom or solvergeom==XYgeom
    bg%invvol(:) = 1._8 / (dr * dz)
  END if

  IF(solvergeom==RZgeom .or. solvergeom==XZgeom) then
    ixlbndi = boundxy
    ixrbndi = boundxy
    izlbndi = bound0
    izrbndi = boundnz
    ixlbnd = boundxy
    ixrbnd = boundxy
    izlbnd = bound0
    izrbnd = boundnz
    IF(l4symtry .or. (solvergeom==RZgeom .and. bg%rmin==0.)) then
      ixlbndi = neumann
      ixlbnd  = neumann
    END if
  else ! solvergeom==XYgeom
    ixlbndi = boundxy
    ixrbndi = boundxy
    izlbndi = boundxy
    izrbndi = boundxy
    ixlbnd = boundxy
    ixrbnd = boundxy
    izlbnd = boundxy
    izrbnd = boundxy
    IF(l4symtry) then
      ixlbndi = neumann
      ixlbnd  = neumann
    END if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      izlbndi = neumann
      izlbnd  = neumann
    END if
  END if
  bg%ixlbnd = ixlbndi
  bg%ixrbnd = ixrbndi
  bg%izlbnd = izlbndi
  bg%izrbnd = izrbndi

  IF(solvergeom==Zgeom .or. solvergeom==Rgeom) then
    bg%nlevels=nlevels
  else
    call init_bnd(bg,nr,nz,dr,dz,bg%zmin,bg%zmax)
    bg%nlevels=nlevels
  END if

!  do i = 1,bg%nlevels, 1
!    bg%bnd(i)%izlbnd=bg%izlbnd
!    bg%bnd(i)%izrbnd=bg%izrbnd
!  END do

  do i = 1, bg%nlevels
    IF(i==1) then
      b => bg%bndfirst
    else
      b => b%next
    END if
    IF(b%izlbnd==dirichlet)  b%v(:,1)      = v_dirichlet
    IF(b%izrbnd==dirichlet)  b%v(:,b%nz+1) = v_dirichlet
    IF(bg%ixlbnd==dirichlet) b%v(1,:)      = v_dirichlet
    IF(bg%ixrbnd==dirichlet) b%v(b%nr+1,:) = v_dirichlet
  END do
  call setmglevels_rz(bg)
  call mk_grids_ptr()

  if (lverbose>=1) then
    write(o_line,'("Exit init_basegrid")')
    call remark(trim(o_line))
  endif

return
end subroutine init_basegrid

subroutine add_grid(mothergrid,nr,nz,dri,dzi,rmini,zmini,transit_min_r,transit_max_r,transit_min_z,transit_max_z)
implicit none
TYPE(GRIDtype), pointer :: mothergrid
INTEGER(ISZ), INTENT(IN) :: nr, nz, transit_min_r, transit_max_r, transit_min_z, transit_max_z
REAL(8), INTENT(IN) :: dri,dzi,rmini,zmini

TYPE(GRIDtype), pointer :: g  ! new grid
TYPE(GRIDtype), pointer :: gtmp  !temporary grid pointer
INTEGER(ISZ) :: i,j,l,ls,jmin,jmax,lmin,lmax,ratio_r,ratio_z,nzs
TYPE(BNDtype), POINTER :: b
REAL(8) :: dr,dz,rmin,zmin

! adjust new grid boundaries to fall onto mother grid lines
! and recalculate mesh spacing for new grid

  jmin = 1 + NINT( (rmini       -mothergrid%rmin) / mothergrid%dr)
  jmax = 1 + NINT( (rmini+nr*dri-mothergrid%rmin) / mothergrid%dr)
  lmin = 1 + NINT( (zmini       -mothergrid%zmin) / mothergrid%dz)
  lmax = 1 + NINT( (zmini+nz*dzi-mothergrid%zmin) / mothergrid%dz)

  rmin = mothergrid%rmin + (jmin-1) * mothergrid%dr
  zmin = mothergrid%zmin + (lmin-1) * mothergrid%dz

  dr = (jmax-jmin) * mothergrid%dr / nr
  dz = (lmax-lmin) * mothergrid%dz / nz

  if (lverbose>0) then
    write(o_line,'("Adding refinement patch to grid ",i5)') mothergrid%gid(1)
    call remark(trim(o_line))
    write(o_line,'("Initial rmin, zmin, dr, dz ",4(" ",e12.5))') rmini,zmini,dri,dzi
    call remark(trim(o_line))
    write(o_line,'("New     rmin, zmin, dr, dz ",4(" ",e12.5))') rmin,zmin,dr,dz
    call remark(trim(o_line))
  end if

! Allocate new grid and initialize variables.

  g => NewGRIDtype()
  nzs = 0
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
  IF(rmin==0.) then
    g%transit_min_r = 0
  else
    g%transit_min_r = transit_min_r
  END if
  g%transit_max_r = transit_max_r
  g%transit_min_z = transit_min_z
  g%transit_max_z = transit_max_z
  g%nguardx = nguardx
  g%nguardz = nguardz
#ifdef MPIPARALLEL
  g%nzp   = nzpslave(my_index)
  g%nrpar = nr
  g%nzpar = g%nzp
#else
  g%nzp   = nz
  g%nrpar = 0
  g%nzpar = 0
#endif
  call GRIDtypeallot(g)
  g%gid=ngrids+1
!  IF(n_avail_ids==0) then
!    g%gid=grids_nids+1
!    grids_nids=grids_nids+1
!  else
!    g%gid=avail_ids(n_avail_ids)
!    n_avail_ids=n_avail_ids-1
!  END if
  IF(associated(mothergrid%down)) then
    IF(associated(mothergrid%down%next)) then
      mothergrid%down%next%prev => g
      g%next => mothergrid%down%next
    END if
    g%prev => mothergrid%down
    mothergrid%down%next => g
    IF(associated(g%prev%down)) g%down => g%prev%down
  else
    mothergrid%down=>g
    gtmp => mothergrid
    do WHILE(associated(gtmp%next))
      gtmp%next%down => g
      gtmp => gtmp%next
    end do
  END if
  g%up => mothergrid
  ngrids=ngrids+1
  mgridrz_ngrids = ngrids
  g%phi=0.
  g%rho=0.
  g%loc_part=g%gid(1)
  g%loc_part_fd=g%gid(1)

! Assign boundary types:
! boundary condition is the same as basegrid if boundary of
! new grid aligns with boundary of basegrid, and is patchbnd
! otherwise.

  IF(ABS(g%xmin-basegrid%xmin)<0.1*g%dr) then
    g%ixlbnd = basegrid%ixlbnd
    if (g%ixlbnd == dirichlet) g%ixlbnd = patchbnd
  else
    g%ixlbnd = patchbnd
  END if
  IF(ABS(g%rmax-basegrid%rmax)<0.1*g%dr) then
    g%ixrbnd = basegrid%ixrbnd
    if (g%ixrbnd == dirichlet) g%ixrbnd = patchbnd
  else
    g%ixrbnd = patchbnd
  END if
  IF(ABS(g%zmin-basegrid%zmin)<0.1*g%dz) then
    g%izlbnd = basegrid%izlbnd
    if (g%izlbnd == dirichlet) g%izlbnd = patchbnd
  else
    g%izlbnd = patchbnd
  END if
  IF(ABS(g%zmax-basegrid%zmax)<0.1*g%dz) then
    g%izrbnd = basegrid%izrbnd
    if (g%izrbnd == dirichlet) g%izrbnd = patchbnd
  else
    g%izrbnd = patchbnd
  END if

! computes commodity quantities for charge deposition

  g%invdr = 1._8/dr
  g%invdz = 1._8/dz
  IF(solvergeom==RZgeom .or. solvergeom==Rgeom) then
    ! computes divider by cell volumes to get density
    IF(g%rmin==0.) then
      j = 1
      ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
      ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
      ! and Verboncoeur, J. of Comp. Phys.,
      g%invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
      do j = 2, nr+1
        g%invvol(j) = 1._8 / (2._8 * pi * real(j-1,8) * dr * dr * dz)
      end do
    else
      do j = 1, nr+1
        g%invvol(j) = 1._8 / (2._8 * pi * (g%rmin+real(j-1,8)*dr) * dr * dz)
      end do
    END if
    IF(solvergeom==Rgeom) g%invvol = g%invvol * dz
  else ! solvergeom==XZgeom or solvergeom==XYgeom
    g%invvol(:) = 1._8 / (dr * dz)
  END if


  IF(lverbose>=3) then
    write(o_line,'(" Add grid ID: ",i5)') g%gid
    call remark(trim(o_line))
  END if
  IF(solvergeom==Zgeom .or.solvergeom==Rgeom) then
    g%nlevels=nlevels
  else
    call init_bnd(g,nr,nz,dr,dz,g%zmin,g%zmax)
    g%nlevels=nlevels
  END if
  do i = 1,g%nlevels
    IF(i==1) then
      b => g%bndfirst
    else
      b => b%next
    END if
    b%izlbnd=g%izlbnd
    b%izrbnd=g%izrbnd
    IF(b%izlbnd==dirichlet .or. b%izlbnd==patchbnd) b%v(:,1)       = v_dirichlet
    IF(b%izrbnd==dirichlet .or. b%izrbnd==patchbnd) b%v(:,b%nz+1:) = v_dirichlet
    IF(g%ixlbnd==dirichlet .or. g%ixlbnd==patchbnd) b%v(1,:)       = v_dirichlet
    IF(g%ixrbnd==dirichlet .or. g%ixrbnd==patchbnd) b%v(b%nr+1,:)  = v_dirichlet
  END do

! update grid pointers list array

  call mk_grids_ptr()

! initializes loc_part* arrays

  IF(solvergeom==Zgeom .or. solvergeom==Rgeom) then
    lmin=1+nint((g%zmin-g%up%zmin)/g%up%dz)
    lmax=1+nint((g%zmax-g%up%zmin)/g%up%dz)
    g%up%loc_part(1,lmin:lmax-1)=g%gid(1)
    g%up%loc_part_fd(1,lmin+transit_min_z:lmax-1-transit_max_z)=g%gid(1)
  else
!    g%up%loc_part(jmin:jmax-1,lmin:lmax-1)=g%gid(1)
!    g%up%loc_part_fd(jmin+transit_min_r:jmax-1-transit_max_r,lmin+transit_min_z:lmax-1-transit_max_z)=g%gid(1)
  END if

! setup neighbors
  if (lverbose>=3) call remark('set neighbors')
  call set_overlaps(g,'n')
! setup parents
  if (lverbose>=3) call remark('set parents')
  call set_overlaps(g,'p')
! setup children
  if (lverbose>=3) call remark('set children')
  call set_overlaps(g,'c')

!  call print_structure(grid)

  return
end subroutine add_grid

subroutine set_overlaps(g,which)
TYPE(GRIDtype), pointer :: g,g2
CHARACTER(1) :: which
TYPE(OVERLAPtype), POINTER :: n1,n2
LOGICAL :: doit,cond,test
REAL(8) :: rmin,rmax,zmin,zmax
INTEGER(ISZ) :: jmin, jmax, lmin, lmax

  select case (which)
      case ('p')
        g2 => g%up
      case ('n')
        g2 => g%up%down
      case ('c')
        IF(associated(g%down)) then
          g2 => g%down
        else
          return
        END if
      case default
        call remark("Argument WHICH has wrong value in subroutine set_overlaps.")
        call remark("Valid values are p, n or c (for parent, neighbor or child).")
        stop
  end select
  doit = .true.
  do while(doit)
     IF(g2%gid(1)/=g%gid(1)) then
       rmin = MAX(g%rmin,g2%rmin)
       rmax = MIN(g%rmax,g2%rmax)
       zmin = MAX(g%zmin,g2%zmin)
       zmax = MIN(g%zmax,g2%zmax)
!       cond = (rmax>rmin .or. ABS(g%rmin-g2%rmax)<0.5*MIN(g%dr,g2%dr) .or. ABS(g2%rmin-g%rmax)<0.5*MIN(g%dr,g2%dr)) .AND. &
!              (zmax>zmin .or. ABS(g%zmin-g2%zmax)<0.5*MIN(g%dz,g2%dz) .or. ABS(g2%zmin-g%zmax)<0.5*MIN(g%dz,g2%dz))
       jmin = 1+NINT((rmin-g%rmin)/g%dr)
       lmin = 1+NINT((zmin-g%zmin)/g%dz)
       jmax = 1+NINT((rmax-g%rmin)/g%dr)
       lmax = 1+NINT((zmax-g%zmin)/g%dz)
!       cond = rmax>rmin .AND. zmax>zmin
       cond = jmax>jmin .AND. lmax>lmin
       IF(cond) then
         n1 => NewOVERLAPtype();    n2 => NewOVERLAPtype()
         n1%rmin = rmin;            n2%rmin = rmin
         n1%rmax = rmax;            n2%rmax = rmax
         n1%zmin = zmin;            n2%zmin = zmin
         n1%zmax = zmax;            n2%zmax = zmax
         IF(which=='n') then
           n1%nr = NINT((rmax-rmin)/g%dr)
           n1%nz = NINT((zmax-zmin)/g%dz)
           n2%nr = n1%nr
           n2%nz = n1%nz
         else
           n1%nr=0
           n1%nz=0
           n2%nr=0
           n2%nz=0
         END if
         call OVERLAPtypeallot(n1)
         call OVERLAPtypeallot(n2)
         n1%gid => g2%gid
         n2%gid => g%gid
         select case (which)
           case ('p')
             call insert_parent(g,n1)
             call insert_child(g2,n2)
             call set_loc_part(g2,g,rmin,rmax,zmin,zmax)
             IF(lverbose>=3) then
               write(o_line,'(i5," is parent of ",i5)') g2%gid,g%gid
               call remark(trim(o_line))
             END if
           case ('n')
             n1%index = 1
             n2%index = 2
             n2%rho => n1%rho
             call insert_neighbor(g,n1)
             call insert_neighbor(g2,n2)
             IF(lverbose>=3) then
               write(o_line,'(i5," is neighbor of ",i5)') g2%gid,g%gid
               call remark(trim(o_line))
             END if
           case ('c')
             call insert_child(g,n1)
             call insert_parent(g2,n2)
             call set_loc_part(g,g2,rmin,rmax,zmin,zmax)
             IF(lverbose>=3) then
               write(o_line,'(i5," is child of ",i5)') g2%gid,g%gid
               call remark(trim(o_line))
             END if
         end select
       END if
     END if
     IF(associated(g2%next)) then
       g2 => g2%next
     else
       doit = .false.
     END if
  END do

end subroutine set_overlaps

subroutine set_loc_part(g,g2,rmin,rmax,zmin,zmax)
TYPE(GRIDtype), INTENT(IN OUT) :: g,g2
REAL(8), INTENT(IN) :: rmin,rmax,zmin,zmax

INTEGER(ISZ) :: jmin, jmax, lmin, lmax
INTEGER(ISZ) :: jminfd, jmaxfd, lminfd, lmaxfd

  jmin = 1+NINT((rmin-g%rmin)/g%dr)
  lmin = 1+NINT((zmin-g%zmin)/g%dz)
  jmax = 1+NINT((rmax-g%rmin)/g%dr)
  lmax = 1+NINT((zmax-g%zmin)/g%dz)
  WHERE(g%loc_part(jmin:jmax-1,lmin:lmax-1)==g%gid(1))
        g%loc_part(jmin:jmax-1,lmin:lmax-1)=g2%gid(1)
  END where
  jmin = 1+NINT((rmin+g2%transit_min_r*g2%dr-g%rmin)/g%dr)
  lmin = 1+NINT((zmin+g2%transit_min_z*g2%dz-g%zmin)/g%dz)
  jmax = 1+NINT((rmax-g2%transit_max_r*g2%dr-g%rmin)/g%dr)
  lmax = 1+NINT((zmax-g2%transit_max_z*g2%dz-g%zmin)/g%dz)
  jminfd = MAX(jmin,     1+g%transit_min_r)
  jmaxfd = MIN(jmax,g%nr+1-g%transit_max_r)
  lminfd = MAX(lmin,     1+g%transit_min_z)
  lmaxfd = MIN(lmax,g%nz+1-g%transit_max_z)
  WHERE(g%loc_part_fd(jminfd:jmaxfd-1,lminfd:lmaxfd-1)==g%gid(1))
        g%loc_part_fd(jminfd:jmaxfd-1,lminfd:lmaxfd-1)=g2%gid(1)
  END where

end subroutine set_loc_part

subroutine insert_neighbor(g,n)
TYPE(GRIDtype), pointer :: g
TYPE(OVERLAPtype), POINTER :: n

  IF(associated(g%neighbors)) then
    IF(associated(g%neighbors%next)) then
      n%next => g%neighbors%next
    END if
    g%neighbors%next => n
  else
    g%neighbors => n
  END if
end subroutine insert_neighbor

subroutine insert_parent(g,n)
TYPE(GRIDtype), pointer :: g
TYPE(OVERLAPtype), POINTER :: n

  IF(associated(g%parents)) then
    IF(associated(g%parents%next)) then
      n%next => g%parents%next
    END if
    g%parents%next => n
  else
    g%parents => n
  END if
end subroutine insert_parent

subroutine insert_child(g,n)
TYPE(GRIDtype), pointer :: g
TYPE(OVERLAPtype), POINTER :: n

  IF(associated(g%children)) then
    IF(associated(g%children%next)) then
      n%next => g%children%next
    END if
    g%children%next => n
  else
    g%children => n
  END if
end subroutine insert_child

subroutine exchange_rho_between_neighbors()
implicit none
TYPE(GRIDtype), pointer :: g
TYPE(OVERLAPtype), POINTER :: n
LOGICAL :: doitl, doitg, doitn
INTEGER(ISZ) :: i, phase, jmin, jmax, lmin, lmax

do phase = 1, 2
  doitg = .true.
  IF(associated(basegrid%down)) then
    doitl = .true.
    g => basegrid%down
  else
    doitl = .false.
  END if

  do WHILE(doitl)
    doitg = .true.
    do WHILE(doitg)
      IF(associated(g%neighbors)) then
        doitn = .true.
        n => g%neighbors
      else
        doitn = .false.
      END if
      do WHILE(doitn)
        i = n%index
        jmin = 1+NINT((n%rmin-g%rmin)/g%dr)
        lmin = 1+NINT((n%zmin-g%zmin)/g%dz)
        jmax = 1+NINT((n%rmax-g%rmin)/g%dr)
        lmax = 1+NINT((n%zmax-g%zmin)/g%dz)
        IF(phase==1) then
          n%rho(:,:,i) = g%rho(jmin:jmax,lmin:lmax)
        else
          i = 3-i
          g%rho(jmin:jmax,lmin:lmax) = g%rho(jmin:jmax,lmin:lmax) + n%rho(:,:,i)
        END if
        IF(associated(n%next)) then
          n => n%next
        else
          doitn = .false.
        END if
      end do
      IF(associated(g%next)) then
        g => g%next
      else
        doitg = .false.
      END if
    end do
    IF(associated(g%down)) then
      g => g%down
    else
      doitl = .false.
    END if
  end do
end do

end subroutine exchange_rho_between_neighbors

subroutine children_send_rho_to_parents()
implicit none
TYPE(GRIDtype), pointer :: g, gp
TYPE(OVERLAPtype), POINTER :: p
LOGICAL :: doitl, doitg, doitp
INTEGER(ISZ) :: i, jmin, jmax, lmin, lmax
REAL(8), DIMENSION(:,:), ALLOCATABLE :: rhotmp
LOGICAL(ISZ), DIMENSION(:,:), ALLOCATABLE :: rhodep

  IF(.not. associated(basegrid%down)) return

  ! point g to last level
  g => basegrid%down
  do WHILE(associated(g%down))
    g => g%down
  END do

  doitl = .true.
  do WHILE(doitl)
    doitg = .true.
    do WHILE(doitg)
      IF(associated(g%parents)) then
        doitp = .true.
        p => g%parents
!        gp => grids_ptr(p%gid(1))%grid
        ALLOCATE(rhodep(1:g%nr+1,1:g%nz+1))
        rhodep = .true.
      else
        doitp = .false.
      END if
      do WHILE(doitp)
        gp => grids_ptr(p%gid(1))%grid
        i = p%index
        jmin = 1+NINT((p%rmin-g%rmin)/g%dr)
        lmin = 1+NINT((p%zmin-g%zmin)/g%dz)
        jmax = 1+NINT((p%rmax-g%rmin)/g%dr)
        lmax = 1+NINT((p%zmax-g%zmin)/g%dz)
        ALLOCATE(rhotmp(1+jmax-jmin,1+lmax-lmin))
        WHERE(rhodep(jmin:jmax,lmin:lmax))
           rhotmp = g%rho(jmin:jmax,lmin:lmax)
        ELSEWHERE
           rhotmp = 0.
        END WHERE
        IF(solvergeom==RZgeom) then
          call deposit_rz(unew=gp%rho, uold=rhotmp, &
                          invvolnew=gp%invvol, invvolold=g%invvol(jmin:jmax), &
!                          xminold=g%rmin, xmaxold=g%rmax, zminold=g%zmin, zmaxold=g%zmax, &
                          xminold=g%rmin+(jmin-1)*g%dr, xmaxold=g%rmin+(jmax-1)*g%dr, &
                          zminold=g%zmin+(lmin-1)*g%dz, zmaxold=g%zmin+(lmax-1)*g%dz, &
                          xminnew=gp%rmin, xmaxnew=gp%rmax, zminnew=gp%zmin, zmaxnew=gp%zmax)
        else
          call deposit(unew=gp%rho, uold=rhotmp, &
!                     xminold=g%rmin, xmaxold=g%rmax, zminold=g%zmin, zmaxold=g%zmax, &
                     xminold=g%rmin+(jmin-1)*g%dr, xmaxold=g%rmin+(jmax-1)*g%dr, &
                     zminold=g%zmin+(lmin-1)*g%dz, zmaxold=g%zmin+(lmax-1)*g%dz, &
                     xminnew=gp%rmin, xmaxnew=gp%rmax, zminnew=gp%zmin, zmaxnew=gp%zmax)
        END if
        rhodep(jmin:jmax,lmin:lmax) = .false.
        DEALLOCATE(rhotmp)
        IF(associated(p%next)) then
          p => p%next
        else
          doitp = .false.
        END if
      end do
      IF(ALLOCATED(rhodep)) DEALLOCATE(rhodep)
      IF(associated(g%next)) then
        g => g%next
      else
        doitg = .false.
      END if
    end do
    IF(associated(g%up)) then
      g => g%up
    else
      doitl = .false.
    END if
  end do

end subroutine children_send_rho_to_parents

RECURSIVE subroutine print_structure(grid)
implicit none
TYPE(GRIDtype), pointer :: grid

IF(level_del_grid==0) then
  write(o_line,*) 'ngrids = ',ngrids
  call remark(trim(o_line))
END if
level_del_grid=level_del_grid+1
write(o_line,'("{",i5,":",i5,"}")') level_del_grid,grid%gid
call remark(trim(o_line))
IF(associated(grid%next)) then
  write(o_line,*) 'next'
  call remark(trim(o_line))
  call print_structure(grid%next)
END if
IF(associated(grid%down)) then
  write(o_line,*) 'down'
  call remark(trim(o_line))
  call print_structure(grid%down)
END if
level_del_grid=level_del_grid-1

end subroutine print_structure

recursive subroutine del_grid_old(g,next_too)
implicit none
TYPE(GRIDtype), pointer :: g
LOGICAL(ISZ),OPTIONAL :: next_too
INTEGER(ISZ):: i
TYPE(BNDtype), pointer :: b
TYPE(GRIDtype), pointer :: gup

  IF(associated(g%down)) call del_grid_old(g%down,.true.)
  IF(PRESENT(next_too)) then
    IF(associated(g%next)) call del_grid_old(g%next,.true.)
  END if

  WHERE(g%up%loc_part==g%gid(1))
    g%up%loc_part = g%up%gid(1)
  END where
  WHERE(g%up%loc_part_fd==g%gid(1))
    g%up%loc_part_fd = g%up%gid(1)
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
  avail_ids(n_avail_ids)=g%gid(1)
!   nullify(g%next,g%prev,g%down,g%up)
  IF(solvergeom/=Zgeom .and. solvergeom/=Rgeom) then
    do i = 1, g%nlevels
      IF(i==1) then
        b => g%bndlast
      else
        b => b%prev
        call ReleaseBNDtype(b%next)
      END if
      call del_cnds(b)
    end do
    call ReleaseBNDtype(b)
  END if
  NULLIFY(g%up%down)
  call ReleaseGRIDtype(g)

  ngrids=ngrids-1
  mgridrz_ngrids = ngrids
  IF(level_del_grid==1) call mk_grids_ptr()

  level_del_grid=level_del_grid-1

  return
end subroutine del_grid_old

RECURSIVE subroutine assign_grids_ptr(grid,godown)
implicit none
TYPE(GRIDtype), pointer :: grid
LOGICAL :: godown

  grids_ptr(grid%gid(1))%grid => grid
  IF(associated(grid%next)) call assign_grids_ptr(grid%next,.false.)
  IF(godown .and. associated(grid%down)) call assign_grids_ptr(grid%down,.true.)

return
end subroutine assign_grids_ptr

subroutine init_bnd(g,nr,nz,dr,dz,zmin,zmax)
! intializes grid quantities according to the number of multigrid levels and grid sizes nx and nz.
USE InGen3d, ONLY:l2symtry, l4symtry
implicit none
TYPE(GRIDtype), pointer :: g
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: dr, dz, zmin, zmax

INTEGER(ISZ) :: i, nrp0, nzp0, nrc, nzc, nrc_old, nzc_old, j, jcoarse, lcoarse
REAL(8) :: drc, dzc
TYPE(BNDtype), pointer :: b

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

  g%nlevels = nlevels

  g%bndfirst => NewBNDType()
  g%bndlast => g%bndfirst
  b => g%bndfirst
  do i = 2, nlevels
    b%next => NewBNDType()
    b%next%prev => b
    b=>b%next
  end do
  g%bndlast => b
  NULLIFY(g%bndfirst%prev)
  NULLIFY(b%next)

!  call GRIDtypechange(g)
!  bndy => g%bnd
!  do i = 1, nlevels
!    call initBNDtypepy(-1,g%bnd(i),g%bnd(i)%id)
!  end do
  bndy_allocated=.true.

  nrc = nr
  nzc = nz
  drc = dr
  dzc = dz
#ifdef MPIPARALLEL
  nworkpproc = 1
#endif
  do i = 1, nlevels
    IF(i==1) then
      b => g%bndfirst
    else
      b => b%next
    END if
    b%l_merged=.false.
    b%izlbnd=izlbndi
    b%izrbnd=izrbndi
    IF(i/=1) then
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
        b%l_merged=.true.
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
    b%zmin = zmslmin(int(my_index/nworkpproc)*nworkpproc)
    b%zmax = zmslmax((1+int(my_index/nworkpproc))*nworkpproc-1)
#else
    b%zmin = zmin
    b%zmax = zmax
#endif

#ifdef MPIPARALLEL
    b%nworkpproc = nworkpproc
    IF(my_index-nworkpproc>=0)      b%izlbnd = -(my_index-nworkpproc+1)
    IF(my_index+nworkpproc<nslaves) b%izrbnd = -(my_index+nworkpproc+1)
#endif
    b%nr = nrc
    b%nz = nzc
    b%dr = drc
    b%dz = dzc
    b%nvlocs = 0
    b%nb_conductors = 0
!    ALLOCATE(b%v(b%nr+1,b%nz+1))
    call BNDtypeallot(b)
    b%v(:,:)=v_vacuum
    IF(lverbose>=3) then
     if(my_index==0) then
      write(o_line,*) i,nrc,nzc,drc,dzc,b%l_merged
      call remark(trim(o_line))
     END if
    endif
  end do
#ifdef MPIPARALLEL
  IF(nslaves>1) then
    do i = 1, nlevels
      IF(i==1) then
        b => g%bndfirst
      else
        b => b%next
      END if
      b%l_powerof2 = .false.
    END do
  END if
#else
  do i = 1, nlevels-1
    IF(i==1) then
      b => g%bndfirst
    else
      b => b%next
    END if
    nrp0 = INT(LOG(REAL(b%nr+1))/LOG(2.))
    nzp0 = INT(LOG(REAL(b%nz+1))/LOG(2.))
!    nrp0 = LOG(REAL(b%nr))/LOG(2.)+0.5
!    nzp0 = LOG(REAL(b%nz))/LOG(2.)+0.5
    IF(2**nrp0/=b%nr.OR.2**nzp0/=b%nz.OR.b%next%nr==b%nr.OR.b%next%nz==b%nz) THEN
      b%l_powerof2=.false.
    else
      b%l_powerof2=.true.
    END if
    IF(l_mgridrz_debug) then
      write(o_line,*) i,b%l_powerof2,2**nrp0,b%nr,2**nzp0,b%nz
      call remark(trim(o_line))
    END if
  end do
  g%bndlast%l_powerof2=.false.
#endif

  IF(mgridrz_deform) then
    mgridrz_nx = nr
    mgridrz_ny = nr
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
TYPE(BNDtype), INTENT(IN OUT) :: bndl
INTEGER(ISZ), INTENT(IN) :: nbnd, ncond
TYPE(CONDtype), POINTER :: c

IF(.not.associated(bndl%cndfirst)) then
  bndl%cndfirst => NewCONDtype()
  bndl%cndlast => bndl%cndfirst
  bndl%nb_conductors = 1
else
  c => NewCONDtype()
  bndl%cndlast%next => c
  c%prev => bndl%cndlast
  bndl%cndlast => c
  bndl%nb_conductors = bndl%nb_conductors + 1
END if

bndl%cndlast%nbbnd = nbnd
bndl%cndlast%nbbndred = nbnd
bndl%cndlast%ncond = ncond
call CONDtypeallot(bndl%cndlast)
IF(nbnd>0) bndl%cndlast%docalc=.false.

end subroutine init_bnd_sublevel

subroutine del_grid(g)
implicit none
TYPE(GRIDtype), pointer :: g,mother

  if (associated(g%up)) then
    mother => g%up
    IF(associated(mother%down,g)) NULLIFY(mother%down)
    do while (associated(mother%next))
      mother => mother%next
      IF(associated(mother%down,g)) NULLIFY(mother%down)
    enddo
  end if
  if (associated(g%prev) .and. associated(g%next)) then
    g%prev%next => g%next
    g%next%prev => g%prev
  else if (associated(g%prev)) then
    NULLIFY(g%prev%next)
  else if (associated(g%next)) then
    NULLIFY(g%next%prev)
  endif

  IF(solvergeom/=Zgeom .and. solvergeom /=Rgeom) call del_grid_bnds(g)
  call del_overlaps(g)

  call ReleaseGRIDtype(g)

  return
end subroutine del_grid

subroutine del_grid_bnds(g)
TYPE(GRIDtype), pointer :: g
TYPE(BNDtype), POINTER :: b,bnext
INTEGER :: i

  b=>g%bndfirst
  NULLIFY(g%bndfirst)
  NULLIFY(g%bndlast)
  do WHILE(associated(b%next))
    bnext => b%next
    NULLIFY(b%next%prev)
    NULLIFY(b%next)
    call del_cnds(b)
    call ReleaseBNDtype(b)
    b => bnext
  end do
  call del_cnds(b)
  call ReleaseBNDtype(b)

  return
end subroutine del_grid_bnds

subroutine del_cnds(bnd)
! initializes quantities for one grid level
implicit none
TYPE(BNDtype), INTENT(IN OUT) :: bnd
TYPE(CONDtype), POINTER :: c,cprev

IF(bnd%nb_conductors==0) return
c => bnd%cndlast
NULLIFY(bnd%cndlast)
NULLIFY(bnd%cndfirst)
do WHILE(associated(c%prev))
  cprev => c%prev
  NULLIFY(c%prev%next)
  NULLIFY(c%prev)
  call ReleaseCONDtype(c)
  bnd%nb_conductors = bnd%nb_conductors - 1
  c=>cprev
end do
call del_cnd(c)
bnd%nb_conductors = bnd%nb_conductors - 1

IF(bnd%nb_conductors>0) then
  write(o_line,*) 'Error in del_cnds: nb_conductors>0.'
  call remark(trim(o_line))
  stop
END if


WHERE(bnd%v == v_cond .OR. bnd%v == v_bnd)
  bnd%v = v_vacuum
END WHERE

end subroutine del_cnds

subroutine del_cnd(c)
implicit none
TYPE(CONDtype), POINTER :: c

  call ReleaseCONDtype(c)

end subroutine del_cnd

subroutine del_overlaps(g)
implicit none
TYPE(GRIDtype) :: g

IF(associated(g%neighbors)) then
  call del_overlap(g%neighbors)
  NULLIFY(g%neighbors)
endif
IF(associated(g%parents)) then
  call del_overlap(g%parents)
  NULLIFY(g%parents)
endif
IF(associated(g%children)) then
  call del_overlap(g%children)
  NULLIFY(g%children)
endif

end subroutine del_overlaps

recursive subroutine del_overlap(o)
implicit none
TYPE(OVERLAPtype), pointer :: o

IF(associated(o%next)) then
  call del_overlap(o%next)
  NULLIFY(o%next)
endif
call ReleaseOVERLAPtype(o)

end subroutine del_overlap

function expandwguard(f)
! expand field from a grid to finer one. Each dimension is assumed to be a power of two.
implicit none

REAL(8), DIMENSION(0:,0:), INTENT(IN) :: f
REAL(8), DIMENSION(0:2*(SIZE(f,1)-2),0:2*(SIZE(f,2)-2)) :: expandwguard

INTEGER(ISZ) :: j, l, nxe, nze


IF(l_mgridrz_debug) then
  write(o_line,*) 'enter expand, level = ',level
  call remark(trim(o_line))
END if

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit expand, level = ',level
  call remark(trim(o_line))
END if

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter expand, level = ',level
  call remark(trim(o_line))
END if

nxold = SIZE(uold,1) - 1 - 2
nzold = SIZE(uold,2) - 1 - 2

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nri=2
else
  nri=1
END if
IF(solvergeom==RZgeom) nri=1
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nrf=nxnew
else
  nrf=nxnew+1
END if
!IF(bndy(level)%izlbnd==dirichlet) then
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
END if
!IF(bndy(level)%izrbnd==dirichlet) then
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit expand, level = ',level
  call remark(trim(o_line))
END if

return
END function expandwguard_any

subroutine add_and_expand(unew, uold, bnd, nxnew, nznew, nxold, nzold, &
                          xminnew, xmaxnew, zminnew, zmaxnew, xminold,  xmaxold,  zminold,  zmaxold, &
                          ixlbnd, ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
TYPE(BNDtype) :: bnd
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:nxnew+2,0:nznew+2), INTENT(INOUT) :: unew
REAL(8), DIMENSION(0:nxold+2,0:nzold+2), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew

INTEGER(ISZ) :: nri, nrf, nzi, nzf
INTEGER(ISZ) :: jnew, knew, j, k, jp, kp, ii, ic
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+1), koldp(nznew+1)

t_before = wtime()

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter expand, level = ',level
  call remark(trim(o_line))
END if

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nri=2
else
  nri=1
END if
IF(solvergeom==RZgeom) nri=1
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nrf=nxnew
else
  nrf=nxnew+1
END if
!IF(bndy(level)%izlbnd==dirichlet) then
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
END if
!IF(bndy(level)%izrbnd==dirichlet) then
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit expand, level = ',level
  call remark(trim(o_line))
END if

t_expand = t_expand + wtime()-t_before
return
END subroutine add_and_expand

function expandwguardandbnd_any(uold, bnd, nxnew, nznew, xminnew, xmaxnew, zminnew, zmaxnew, &
                                                         xminold, xmaxold, zminold, zmaxold, &
                                                         ixlbnd, ixrbnd, izlbnd, izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(0:,0:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8) :: expandwguardandbnd_any(0:nxnew+2,0:nznew+2)
TYPE(BNDtype) :: bnd

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jnew, knew, i, j, k, jp, kp, nri, nrf, nzi, nzf
REAL(8) :: x, z, xx, zz, invdxold, invdzold, dxnew, dznew, delx, delz, odelx, odelz
REAL(8) :: ddx(nxnew+1), oddx(nxnew+1), ddz(nznew+1), oddz(nznew+1)
INTEGER(ISZ) :: jold(nxnew+1), kold(nznew+1), joldp(nxnew+1), koldp(nznew+1)

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter expand, level = ',level
  call remark(trim(o_line))
END if
expandwguardandbnd_any = 0._8

nxold = SIZE(uold,1) - 1 - 2
nzold = SIZE(uold,2) - 1 - 2

dxnew = (xmaxnew-xminnew) / nxnew
dznew = (zmaxnew-zminnew) / nznew
invdxold = REAL(nxold,8)/(xmaxold-xminold)
invdzold = REAL(nzold,8)/(zmaxold-zminold)

IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nri=2
else
  nri=1
END if
IF(solvergeom==RZgeom) nri=1
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nrf=nxnew
else
  nrf=nxnew+1
END if
!IF(bndy(level)%izlbnd==dirichlet) then
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
END if
!IF(bndy(level)%izrbnd==dirichlet) then
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit expand, level = ',level
  call remark(trim(o_line))
END if

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
     IF((knew==1.AND.(izlbnd==dirichlet.or.izlbnd==patchbnd)).OR.(knew==nznew+1.AND.(izrbnd==dirichlet.or.izrbnd==patchbnd))) then
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
    IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
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
!     write(o_line,*) 'knew',knew
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
     IF((knew==1.AND.(izlbnd==dirichlet.or.izlbnd==patchbnd)).OR.(knew==nznew+1.AND.(izrbnd==dirichlet.or.izrbnd==patchbnd))) then
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
    IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
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
     IF((knew==1.AND.(izlbnd==dirichlet.or.izlbnd==patchbnd)).OR.(knew==nznew+1.AND.(izrbnd==dirichlet.or.izrbnd==patchbnd))) then
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
     IF((knew==1.AND.(izlbnd==dirichlet.or.izlbnd==patchbnd)).OR.(knew==nznew+1.AND.(izrbnd==dirichlet.or.izrbnd==patchbnd))) then
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

function restrict_pof2(f, ixrbnd, izlbnd, izrbnd)
! restrict field from one grid to a coarser one. Each dimension is assumed to be a power of two.
implicit none

REAL(8), DIMENSION(1:,1:), INTENT(IN) :: f
INTEGER(ISZ), INTENT(IN) :: ixrbnd, izlbnd, izrbnd
REAL(8), DIMENSION(SIZE(f,1)/2+1, SIZE(f,2)/2+1) :: restrict_pof2

INTEGER(ISZ) :: nx, nz, nxi, nzi, j, l, jj, ll
REAL(8) :: flz, frz, aflz1, aflz2, afrz1, afrz2

t_before = wtime()

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter restrict_pof2, level = ',level
  call remark(trim(o_line))
END if

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
  IF(bndcurrent%next%izlbnd<0) then
    flz=0.5_8
    aflz1 = 1._8
    aflz2 = 4._8/3._8
  endif
  IF(bndcurrent%next%izrbnd<0) then
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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit restrict_pof2, level = ',level
  call remark(trim(o_line))
END if

t_restrict = t_restrict + wtime()-t_before
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

t_before = wtime()

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold, ixrbnd, izlbnd, izrbnd)
!  return
!END if

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter restrict, level = ',level
  call remark(trim(o_line))
END if

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit restrict, level = ',level
  call remark(trim(o_line))
END if

t_restrict = t_restrict + wtime()-t_before
return
END function restrict

!pgi$r nobounds
subroutine subrestrict(unew, uold, nxnew, nznew, &
                       xminold, xmaxold, zminold, zmaxold,  &
                       xminnew, xmaxnew, zminnew, zmaxnew, &
                       izlbnd, izrbnd)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, izlbnd, izrbnd
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), DIMENSION(1:nxnew+1,1:nznew+1), INTENT(OUT) :: unew
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew
REAL(8),allocatable :: rap(:,:)

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jold, kold, j, k, jp, kp
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

t_before = wtime()

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold, ixrbnd, izlbnd, izrbnd)
!  return
!END if

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter restrict, level = ',level
  call remark(trim(o_line))
END if

nxold = SIZE(uold,1) - 1
nzold = SIZE(uold,2) - 1

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1))
ALLOCATE(rap(1:nxnew+1,1:nznew+1))

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
DEALLOCATE(rap)

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit restrict, level = ',level
  call remark(trim(o_line))
END if

t_restrict = t_restrict + wtime()-t_before
return
END subroutine subrestrict

subroutine restrictlist(unew, uold, rhs, bnd, nxnew, nznew, nxold, nzold, voltfact, dr, dz, &
                        xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! restrict field from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
INTEGER(ISZ), INTENT(IN) :: nxnew, nznew, nxold, nzold
TYPE(BNDtype) :: bnd
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
TYPE(CONDtype), POINTER :: c

t_before = wtime()

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter restrict, level = ',level
  call remark(trim(o_line))
END if

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
    c => bnd%cndfirst
  else
    c => c%next
  END if
  nlocs = nlocs + c%nbbnd
END do

ALLOCATE(res(nlocs),jlocs(nlocs),klocs(nlocs))

call residbndrzwguard_list(res(1),jlocs(1),klocs(1),nlocs,f=uold(0,0),rhs=rhs(1,1), &
                           bnd=bnd,nr=nxold,nz=nzold,dr=dr,dz=dz,rmin=0._8,voltfact=voltfact)

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit restrict, level = ',level
  call remark(trim(o_line))
END if

t_restrict = t_restrict + wtime()-t_before
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
TYPE(BNDtype) :: bnd

INTEGER(ISZ) :: nxold, nzold
INTEGER(ISZ) :: jold, kold, j, k, jp, kp, ic, ii, itot
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz, dxnew, dznew, u1, u2, u3, u4, q
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp
LOGICAL(ISZ) :: l_dxm, l_dxp, l_dzm, l_dzp
INTEGER(ISZ),allocatable :: cnt(:,:)
TYPE(CONDtype), POINTER :: c

t_before = wtime()

!IF(bndy(level)%l_powerof2) then
!  restrict = restrict_pof2(uold)
!  return
!END if
itot = 0
IF(l_mgridrz_debug) then
  write(o_line,*) 'enter restrict, level = ',level
  call remark(trim(o_line))
END if

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
    IF((bndcurrent%izlbnd<0.and.kold>1).or.(bndcurrent%izrbnd<0.and.kold<nzold+1)) then
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
    c => bnd%cndfirst
  else
    c => c%next
  END if
  do ii = 1, c%nbbnd
    jold = c%jj(ii)
    kold = c%kk(ii)
    IF(.NOT.(bnd%v(jold,kold)==v_bnd.and.c%docalc(ii))) cycle
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
    IF(c%dxm(ii)<delx*dxnew)  l_dxm = .false.
    IF(c%dxp(ii)<odelx*dxnew) l_dxp = .false.
    IF(c%dzm(ii)<delz*dznew)  l_dzm = .false.
    IF(c%dzp(ii)<odelz*dznew) l_dzp = .false.
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
    c => bnd%cndfirst
  else
    c => c%next
  END if
  do ii = 1, c%ncond
    jold = c%jcond(ii)
    kold = c%kcond(ii)
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
          IF(j==jold.and.k==kold) then
            write(o_line,*) level,': #1# ',j,k,cnt(j,k)
            call remark(trim(o_line))
          END if
        END do
      END do
      do ic = 1, bnd%nb_conductors
        IF(ic==1) then
          c => bnd%cndfirst
        else
          c => c%next
        END if
        do ii = 1, c%nbbnd
          jold = c%jj(ii)
          kold = c%kk(ii)
          IF(.NOT.(bnd%v(jold,kold)==v_bnd.and.c%docalc(ii))) cycle
          IF(j==jold.and.k==kold) then
            write(o_line,*) level,': #2# ',j,k,cnt(j,k)
            call remark(trim(o_line))
          END if
        END do
        do ii = 1, c%ncond
          jold = c%jcond(ii)
          kold = c%kcond(ii)
          IF(j==jold.and.k==kold) then
            write(o_line,*) level,': #3# ',j,k,cnt(j,k)
            call remark(trim(o_line))
          END if
        END do
      END do
    END if
  END do
END do
DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp, cnt)

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit restrict, level = ',level
  call remark(trim(o_line))
END if

t_restrict = t_restrict + wtime()-t_before
return
END function restrict_wbnd

subroutine deposit(unew, uold, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
! deposit rho from one grid to a coarser one. Each dimension may have any number of cells.
implicit none
REAL(8), DIMENSION(1:,1:), INTENT(IN OUT) :: unew
REAL(8), DIMENSION(1:,1:), INTENT(IN) :: uold
REAL(8), INTENT(IN) :: xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew

INTEGER(ISZ) :: nxold, nzold, nxnew, nznew
INTEGER(ISZ) :: jold, kold, j, k, jp, kp
REAL(8) :: dxold, dzold, invdxnew, invdznew, x, z, delx, delz, odelx, odelz
REAL(8), ALLOCATABLE, DIMENSION(:) :: ddx, ddz, oddx, oddz
INTEGER(ISZ), ALLOCATABLE, DIMENSION(:) :: jnew, knew, jnewp, knewp

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter deposit'
  call remark(trim(o_line))
END if

nxold = SIZE(uold,1) - 1
nzold = SIZE(uold,2) - 1
nxnew = SIZE(unew,1) - 1
nznew = SIZE(unew,2) - 1

ALLOCATE(ddx(nxold+1), ddz(nzold+1), oddx(nxold+1), oddz(nzold+1), &
         jnew(nxold+1), knew(nzold+1), jnewp(nxold+1), knewp(nzold+1))

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
  ddz(kold) = ddz(kold) * invdznew*dzold
  oddz(kold) = oddz(kold) * invdznew*dzold
END do

do jold = 1, nxold+1
  x = xminold + (jold-1)*dxold
  jnew(jold) = MIN(1 + INT((x-xminnew) * invdxnew), nxnew)
  ddx(jold) = (x-xminnew) * invdxnew-real(jnew(jold)-1)
  jnewp(jold) = jnew(jold)+1
  oddx(jold) = (1._8-ddx(jold))
  ddx(jold) = ddx(jold) * invdxnew*dxold
  oddx(jold) = oddx(jold) * invdxnew*dxold
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
  end do
end do

DEALLOCATE(ddx, ddz, oddx, oddz, jnew, knew, jnewp, knewp)

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit deposit'
  call remark(trim(o_line))
END if

return
END subroutine deposit

subroutine deposit_rz(unew, uold, invvolnew, invvolold, xminold, xmaxold, zminold, zmaxold, xminnew, xmaxnew, zminnew, zmaxnew)
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

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter deposit'
  call remark(trim(o_line))
END if

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit deposit'
  call remark(trim(o_line))
END if

return
END subroutine deposit_rz

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter deposit'
  call remark(trim(o_line))
END if

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

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit deposit'
  call remark(trim(o_line))
END if

return
END subroutine deposit_z

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

subroutine relaxbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,rmin,nc,voltfact,mgparam, ixlbnd, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN OUT) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact, mgparam, rmin
TYPE(BNDtype), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, lswinit, redblack, iil, iiu, ic, nri, nrf, nzi, nzf
REAL(8) :: dt, dt0
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs, r
TYPE(CONDtype), POINTER :: c

t_before = wtime()

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter relax, level = ',level
  call remark(trim(o_line))
END if

! define CFL
dt = mgparam/(2._8/dr**2+2._8/dz**2)
dt0 = mgparam/(4._8/dr**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cf0 = 1._8-2._8*dt/dr**2-2._8*cfz
cfrhs = dt*inveps0
do j = 2, nr+1
  r = rmin+REAL(j-1,8)*dr
  cfrp(j) = dt * (1._8+0.5_8*dr/r) / dr**2
  cfrm(j) = dt * (1._8-0.5_8*dr/r) / dr**2
!  cfrp(j) = dt * (1._8+0.5_8/REAL(j-1,8)) / dr**2
!  cfrm(j) = dt * (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

lswinit = 2
IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nri=2
else
  nri=1
  lswinit = 3-lswinit
END if
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nrf=nr-1
else
  nrf=nr
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
  lswinit = 3-lswinit
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
  nzf=nz-1
else
  nzf=nz
END if

do i = 1, nc

lsw = lswinit
do redblack = 1, 2
  jsw = lsw
  do ic = 1, bnd%nb_conductors
    IF(ic==1) then
      c => bnd%cndfirst
    else
      c => c%next
    END if

    IF(redblack==1) THEN !red
      iil=1
      iiu=c%nbbndred
    else !black
      iil=c%nbbndred+1
      iiu=c%nbbnd
    ENDif
    do ii = iil, iiu
      j = c%jj(ii)
      l = c%kk(ii)
      IF(j==1) then
        IF(c%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*c%dt(ii)*( &
                 c%cf0(ii)*f(j,l) &
               + c%cfxp(ii)*f(j+1,l) &
               + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
               + voltfact*(c%phi0xp(ii) &
               + c%phi0zm(ii)+c%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
      else
        IF(c%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*c%dt(ii)*( &
                 c%cf0(ii)*f(j,l) &
               + c%cfxp(ii)*f(j+1,l)+c%cfxm(ii)*f(j-1,l) &
               + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
               + voltfact*(c%phi0xm(ii)+c%phi0xp(ii) &
               + c%phi0zm(ii)+c%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
      END if
    ENDDO
  END do
  IF(vlocs) then
    IF(redblack==1) THEN !red
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
      IF(nri==1 .and. jsw==2) then! origin
        j = 1
        IF(bnd%v(j,l)==v_vacuum) &
        f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                   + 4._8*dt0*f(j+1,l)/dr**2   &
                                   + (dt0/dz**2)*(f(j,l+1)+f(j,l-1)) &
                                   + dt0*rhs(j,l)*inveps0
      END if
      do j = nri+jsw, nrf+1, 2
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
!  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
!  call exchange_fbndz_rb(f,izlbnd,izrbnd,(redblack-1))
  call exchange_fbndz_rb(f,izlbnd,izrbnd,lsw-1)
!  call exchange_fbndz(f,izlbnd,izrbnd)
#endif

call updateguardcellsrz(f=f, ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !redblack=1, 2

END do !i=1, nc

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit relax, level = ',level
  call remark(trim(o_line))
END if

t_relax = t_relax + wtime()-t_before
return
END subroutine relaxbndrzwguard

subroutine relaxbndrzwguard_jump(f,rhs,maxjump,curjump,bnd,nr,nz,dr,dz,rmin,nc,voltfact,mgparam, ixlbnd, ixrbnd, izlbnd, izrbnd)
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nc, ixlbnd, ixrbnd, izlbnd, izrbnd, curjump
REAL(8), INTENT(IN OUT) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
INTEGER(ISZ), INTENT(IN) :: maxjump(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, voltfact, mgparam, rmin
TYPE(BNDtype), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, lswinit, redblack, iil, iiu, ic, nri, nrf, nzi, nzf, jump
REAL(8) :: dt, dt0, dtj,drj,dzj, r
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, cfrhs
TYPE(CONDtype), pointer :: c

t_before = wtime()

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter relax, level = ',level
  call remark(trim(o_line))
END if

write(o_line,*) curjump,nc
call remark(trim(o_line))

! define CFL
dt = mgparam/(2._8/dr**2+2._8/dz**2)
dt0 = mgparam/(4._8/dr**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cf0 = 1._8-2._8*dt/dr**2-2._8*cfz
cfrhs = dt*inveps0
do j = 2, nr+1
  r = rmin+REAL(j-1,8)*dr
  cfrp(j) = dt * (1._8+0.5_8*dr/r) / dr**2
  cfrm(j) = dt * (1._8-0.5_8*dr/r) / dr**2
!  cfrp(j) = dt * (1._8+0.5_8/REAL(j-1,8)) / dr**2
!  cfrm(j) = dt * (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

lswinit = 2
IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nri=2
else
  nri=1
  lswinit = 3-lswinit
END if
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nrf=nr-1
else
  nrf=nr
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
  lswinit = 3-lswinit
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
  nzf=nz-1
else
  nzf=nz
END if

do i = 1, nc

lsw = lswinit
do redblack = 1, 2
  jsw = lsw
  do ic = 1, bnd%nb_conductors
    IF(ic==1) then
      c => bnd%cndfirst
    else
      c => c%next
    END if

    IF(redblack==1) THEN !red
      iil=1
      iiu=c%nbbndred
    else !black
      iil=c%nbbndred+1
      iiu=c%nbbnd
    ENDif
    do ii = iil, iiu
      j = c%jj(ii)
      l = c%kk(ii)
      IF(j==1) then
        IF(c%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*c%dt(ii)*( &
                 c%cf0(ii)*f(j,l) &
               + c%cfxp(ii)*f(j+1,l) &
               + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
               + voltfact*(c%phi0xp(ii) &
               + c%phi0zm(ii)+c%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
      else
        IF(c%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*c%dt(ii)*( &
                 c%cf0(ii)*f(j,l) &
               + c%cfxp(ii)*f(j+1,l)+c%cfxm(ii)*f(j-1,l) &
               + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
               + voltfact*(c%phi0xm(ii)+c%phi0xp(ii) &
               + c%phi0zm(ii)+c%phi0zp(ii)) &
               + rhs(j,l)*inveps0)
      END if
    ENDDO
  END do
  IF(vlocs) then
    IF(redblack==1) THEN !red
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
      IF(nri==1 .and. jsw==2) then! origin
        j = 1
        jump = MIN(curjump,maxjump(j,l))
        IF(bnd%v(j,l)==v_vacuum) &
        f(j,l) = (1._8-4._8*dt0/dr**2-2._8*dt0/dz**2) * f(j,l) &
                                   + 4._8*dt0*f(j+jump,l)/dr**2   &
                                   + (dt0/dz**2)*(f(j,l+jump)+f(j,l-jump)) &
                                   + dt0*rhs(j,l)*inveps0
!                                   + (1./real(jump**2,8))*dt0*rhs(j,l)*inveps0
      END if
      do j = nri+jsw, nrf+1, 2
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
!  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
!  call exchange_fbndz_rb(f,izlbnd,izrbnd,lsw-1)
  call exchange_fbndz_rb(f,izlbnd,izrbnd,(redblack-1))
#endif

END do !redblack=1, 2

call updateguardcellsrz(f=f, ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !i=1, nc

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit relax, level = ',level
  call remark(trim(o_line))
END if

t_relax = t_relax + wtime()-t_before
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
TYPE(BNDtype), INTENT(IN OUT) :: bnd

INTEGER(ISZ) :: i, j, l, ii, jsw, lsw, lswinit, redblack, iil, iiu, ic, nxi, nxf, nzi, nzf
REAL(8) :: dt, cf0, cfx, cfz, cfrhs
TYPE(CONDtype), pointer :: c

t_before = wtime()

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter relax, level = ',level
  call remark(trim(o_line))
END if

! define CFL
dt = mgparam/(2._8/dx**2+2._8/dz**2)
! define coefficients
cfz = dt / dz**2
cfx = dt / dx**2
cf0 = 1._8-2._8*(cfx+cfz)
cfrhs = dt*inveps0

lswinit = 1
IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nxi=2
else
  nxi=1
  lswinit = 3-lswinit
END if
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nxf=nx-1
else
  nxf=nx
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
  lswinit = 3-lswinit
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
  nzf=nz-1
else
  nzf=nz
END if

do i = 1, nc

lsw = lswinit
do redblack = 1, 2
  jsw = lsw
  do ic = 1, bnd%nb_conductors
    IF(ic==1) then
      c => bnd%cndfirst
    else
      c => c%next
    END if

    IF(redblack==1) THEN !red
      iil=1
      iiu=c%nbbndred
    else !black
      iil=c%nbbndred+1
      iiu=c%nbbnd
    ENDif
    do ii = iil, iiu
      j = c%jj(ii)
      l = c%kk(ii)
        IF(c%docalc(ii).and.bnd%v(j,l)==v_bnd) &
        f(j,l) = f(j,l) + mgparam*c%dt(ii)*( &
                 c%cf0(ii)*f(j,l) &
               + c%cfxp(ii)*f(j+1,l)+c%cfxm(ii)*f(j-1,l) &
               + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
               + voltfact*(c%phi0xm(ii)+c%phi0xp(ii) &
               +           c%phi0zm(ii)+c%phi0zp(ii)) &
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
!  call exchange_fbndz_rb(f,izlbnd,izrbnd,1-(redblack-1))
!  call exchange_fbndz_rb(f,izlbnd,izrbnd,(redblack-1))
  call exchange_fbndz_rb(f,izlbnd,izrbnd,lsw-1)
#endif

call updateguardcellsrz(f=f, ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=izlbnd, izrbnd=izrbnd)

END do !redblack=1, 2

END do !i=1, nc

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit relax, level = ',level
  call remark(trim(o_line))
END if

t_relax = t_relax + wtime()-t_before

return
END subroutine relaxbndxzwguard

#ifdef MPIPARALLEL
  subroutine exchange_fbndz(f, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd

    INTEGER(4) :: nr,nz
    INTEGER(4) :: p_up, p_down, nsize
    integer(4) :: mpi_req(2*nslaves+2),mpistatus(MPI_STATUS_SIZE,2*nslaves+2),mpierror,ir

    
!    write(o_line,*) my_index,':enter exchangefbnd'
    ir = 0

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nr = SIZE(f,1)
    nz = SIZE(f,2)-3

    ! send
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_isend(f(0,2),nr,mpi_double_precision,p_down,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_down,ir
!      call remark(trim(o_line))
     end if
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_isend(f(0,nz),nr,mpi_double_precision,p_up,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_up,ir
!      call remark(trim(o_line))
    end if

    ! receive
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_irecv(f(0,nz+2),nr,mpi_double_precision,p_up,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_up,ir
!      call remark(trim(o_line))
    end if     
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_irecv(f(0,0),nr,mpi_double_precision,p_down,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_down,ir
!      call remark(trim(o_line))
    end if

    if(ir>0 .and. l_mpi_barriers) call MPI_WAITALL(ir,mpi_req(1:ir),mpistatus(:,1:ir),mpierror)

  end subroutine exchange_fbndz


  subroutine exchange_fbndz_rb(f, izlbnd, izrbnd, izf)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd, izf

    INTEGER(4) :: nr,nz
    INTEGER(4) :: p_up, p_down
    integer(4) :: mpi_req(4),mpistatus(MPI_STATUS_SIZE,4),mpierror,ir

    real(8), allocatable, dimension(:) :: ftmpd, ftmpu, ftmpds, ftmpus
    
!    write(o_line,*) my_index,':enter exchangefbnd',izf
!    call remark(trim(o_line))
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
      ftmpds = f(izf::2,2)
      call mpi_isend(ftmpds(izf),nr,mpi_double_precision,p_down,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_down,ir
!      call remark(trim(o_line))
    end if
    IF(izrbnd<0) then
      ir = ir + 1
      ftmpus = f(izf::2,nz)
      call mpi_isend(ftmpus(izf),nr,mpi_double_precision,p_up,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_up,ir
!      call remark(trim(o_line))
    end if

    ! receive
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_irecv(ftmpu(izf),nr,mpi_double_precision,p_up,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_up,ir
!      call remark(trim(o_line))
    end if     
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_irecv(ftmpd(izf),nr,mpi_double_precision,p_down,0_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_down,ir
!      call remark(trim(o_line))
    end if

    if(ir>0) call MPI_WAITALL(ir,mpi_req(1),mpistatus(1,1),mpierror)

    IF(izrbnd<0) then
      f(izf::2,nz+2) = ftmpu
    end if     
    IF(izlbnd<0) then
      f(izf::2,0)    = ftmpd
    end if
    deallocate(ftmpd,ftmpu)
    deallocate(ftmpds,ftmpus)
!    write(o_line,*) my_index,':exit exchangefbnd'
!    call remark(trim(o_line))

  end subroutine exchange_fbndz_rb
  subroutine check_fbndz(f, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: f(0:,0:)!f(0:nr+2,0:nz+2)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd

    INTEGER(4) :: nr,nz,i
    INTEGER(ISZ) :: p_up, p_down

    REAL(8), DIMENSION(0:SIZE(f,1)-1) :: fr,fl

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nr = SIZE(f,1)-3
    nz = SIZE(f,2)-3

    ! send
    IF(izlbnd<0) call mpi_send_real_array(f(:,1), p_down, 0)
    IF(izrbnd<0) call mpi_send_real_array(f(:,nz+1), p_up, 0)
!    IF(izlbnd<0) write(o_line,*) my_index, ' send to ',p_down
!    if(izlbnd<0) call remark(trim(o_line))
!    IF(izrbnd<0) write(o_line,*) my_index, ' send to ',p_up
!    if(izrbnd<0) call remark(trim(o_line))

    ! receive
!    IF(bndy(level)%izrbnd<0) write(o_line,*) my_index, ' recv from ',p_up
!    IF(bndy(level)%izlbnd<0) write(o_line,*) my_index, ' recv from ',p_down
    IF(izrbnd<0) fr = mpi_recv_real_array(SIZE(f(:,nz)),p_up,0)
    IF(izlbnd<0) fl = mpi_recv_real_array(SIZE(f(:,0 )),p_down,0)
!    IF(izlbnd<0) write(o_line,*) my_index, ' recv from ',p_down
!    if(izlbnd<0) call remark(trim(o_line))
!    IF(izrbnd<0) write(o_line,*) my_index, ' recv from ',p_up
!    if(izrbnd<0) call remark(trim(o_line))

    IF(izrbnd<0) then
     do i = 1, nr+1
      IF(fr(i)/=f(i,nz+1)) then
        write(o_line,*) 'Error fr mismatch: level ',level,' procs ',my_index,p_up, ' i ',i,fr(i),f(i,nz+1)
        call remark(trim(o_line))
      END if
     end do
    END if
    IF(izlbnd<0) then
     do i = 1, nr+1
      IF(fl(i)/=f(i,1)) then
        write(o_line,*) 'Error fl mismatch: level ',level,' procs ',p_down,my_index, ' i ',i,fl(i),f(i,1)
        call remark(trim(o_line))
      END if
     end do
    END if

    if (l_mpi_barriers) call parallelbarrier()

  end subroutine check_fbndz
  subroutine exchange_resbndz(rho, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: rho(1:,1:)!rho(1:nr+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: izlbnd, izrbnd

    INTEGER(4) :: nr,nz
    INTEGER(4) :: p_up, p_down
    integer(4) :: mpi_req(4),mpistatus(MPI_STATUS_SIZE,4),mpierror,ir,j
    real(8), allocatable, dimension(:) :: fd, fu

!    write(o_line,*) my_index,':enter exchangeres'
!    call remark(trim(o_line))
    ir = 0

    p_up   = -izrbnd-1
    p_down = -izlbnd-1

    nr = SIZE(rho,1)-1
    nz = SIZE(rho,2)-1

    ! send
    IF(izlbnd<0) then
      ir = ir + 1
      call mpi_isend(rho(1,1),nr+1_4,mpi_double_precision,p_down,1_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_down
!      call remark(trim(o_line))
    end if
    IF(izrbnd<0) then
      ir = ir + 1
      call mpi_isend(rho(1,nz+1),nr+1_4,mpi_double_precision,p_up,1_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_up
!      call remark(trim(o_line))
    end if

    ! receive
    IF(izrbnd<0) then
      ir = ir + 1
      allocate(fu(nr+1))
      call mpi_irecv(fu(1),nr+1_4,mpi_double_precision,p_up,1_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_up
!      call remark(trim(o_line))
    end if
    IF(izlbnd<0) then
      ir = ir + 1
      allocate(fd(nr+1))
      call mpi_irecv(fd(1),nr+1_4,mpi_double_precision,p_down,1_4,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_down
!      call remark(trim(o_line))
   end if

!    call parallelbarrier()
    if(ir>0 .and. l_mpi_barriers) call MPI_WAITALL(ir,mpi_req(1),mpistatus(1,1),mpierror)
    IF(izrbnd<0) then
      rho(:,nz+1) = rho(:,nz+1) + fu(:)
      deallocate(fu)
    end if
    IF(izlbnd<0) then
      rho(:,1) = rho(:,1) + fd(:)
      deallocate(fd)
    end if
!    write(o_line,*) my_index,':exit exchangeres'
!    call remark(trim(o_line))

  end subroutine exchange_resbndz

   subroutine merge_work(f, level, izlbnd, izrbnd)
    REAL(8), INTENT(IN OUT) :: f(1:,1:)!f(1:nr+1,1:nz+1)
    INTEGER(ISZ), INTENT(IN) :: level, izlbnd, izrbnd

    INTEGER(4) :: nz, p_up, p_down, j, nr
    integer(4) :: mpi_req(2),mpistatus(MPI_STATUS_SIZE,2),mpierror,ir
    real(8), allocatable, dimension(:,:) :: fd, fu

!    write(o_line,*) my_index,':enter merge'
!    call remark(trim(o_line))
    ir = 0

    nr     = size(f,1)-1
    nz     = bndcurrent%next%nz
    p_up   = -izrbnd-1
    p_down = -izlbnd-1
    IF(MOD(my_index/bndcurrent%nworkpproc,2)==0) then
    ! send up
      ir = ir +1
!      call mpi_isend(f(1:nr+1,1:nz/2+1),(nr+1)*(nz/2+1),mpi_real8,p_up,2,MPI_COMM_WORLD,mpi_req(ir),mpierror)
      call mpi_isend(f(1,1),int((nr+1)*(nz/2+1),4),mpi_double_precision,p_up,2_4,MPI_COMM_WORLD,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_up
!      call remark(trim(o_line))
    ! receive up
      ir = ir +1
      allocate(fu(nr+1,nz/2+1:nz+1))
      fu=0.
      call mpi_irecv(fu(1,nz/2+1),int((nr+1)*(nz/2+1),4),mpi_double_precision,p_up,2_4,mpi_comm_world,mpi_req(ir),mpierror)
!      call mpi_irecv(fu,(nr+1)*(nz/2+1),mpi_real8,p_up,2,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_up
!      call remark(trim(o_line))
    else
    ! send down
      ir = ir +1
!      call mpi_isend(f(1:nr+1,nz/2+1:nz+1),(nr+1)*(nz/2+1),mpi_real8,p_down,2,MPI_COMM_WORLD,mpi_req(ir),mpierror)
      call mpi_isend(f(1,nz/2+1),int(SIZE(f,1)*(nz/2+1),4),mpi_double_precision,p_down,2_4,MPI_COMM_WORLD,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' send to ',p_down
!      call remark(trim(o_line))
    ! receive down
      ir = ir +1
      allocate(fd(nr+1,1:nz/2+1))
      fd=0.
      call mpi_irecv(fd(1,1),int((nr+1)*(nz/2+1),4),mpi_double_precision,p_down,2_4,mpi_comm_world,mpi_req(ir),mpierror)
!      call mpi_irecv(fd,(nr+1)*(nz/2+1),mpi_real8,p_down,2,mpi_comm_world,mpi_req(ir),mpierror)
!      write(o_line,*) my_index, ' recv from ',p_down
!      call remark(trim(o_line))
    END if

    if(ir>0 .and. l_mpi_barriers) call MPI_WAITALL(ir,mpi_req(1),mpistatus(1,1),mpierror)

     IF(MOD(my_index/bndcurrent%nworkpproc,2)==0) then
       f(:,nz/2+2:nz+1) = fu(:,nz/2+2:nz+1)
       f(:,nz/2+1) = f(:,nz/2+1)+fu(:,nz/2+1)
       deallocate(fu)
     else
       f(:,1:nz/2) = fd(:,1:nz/2)
       f(:,nz/2+1) = f(:,nz/2+1)+fd(:,nz/2+1)
       deallocate(fd)
     end if
!    write(o_line,*) my_index,':exit merge'
!    call remark(trim(o_line))

  end subroutine merge_work
#endif

subroutine residbndrzwguard(f,rhs,bnd,nr,nz,dr,dz,rmin,voltfact,l_zerolastz, ixrbnd, izlbnd, izrbnd,res)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
TYPE(BNDtype) :: bnd
REAL(8), INTENT(IN) :: dr, dz,voltfact, rmin
REAL(8), DIMENSION(nr+1,nz+1) :: res
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, l, ii, ic, nrf, nzi, nzf
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, r
TYPE(CONDtype), pointer :: c

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter resid, level = ',level
  call remark(trim(o_line))
END if

cfz = 1._8 / dz**2
cf0 = -2._8/dr**2-2._8*cfz
do j = 2, nr+1
  r = rmin+REAL(j-1,8)*dr
  cfrp(j) = (1._8+0.5_8*dr/r) / dr**2
  cfrm(j) = (1._8-0.5_8*dr/r) / dr**2
!  cfrp(j) = (1._8+0.5_8/REAL(j-1,8)) / dr**2
!  cfrm(j) = (1._8-0.5_8/REAL(j-1,8)) / dr**2
end do

res = 0._8

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    c => bnd%cndfirst
  else
    c => c%next
  END if
  do i = 1, c%ncond
    res(c%jcond(i),c%kcond(i)) = 0._8
  end do
END do

IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nrf=nr-1
else
  nrf=nr
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
  nzf=nz-1
else
  nzf=nz
END if

IF(vlocs) then
  do ii = 1, bnd%nvlocs
    j = bnd%vlocs_j(ii)
    l = bnd%vlocs_k(ii)
    IF(j==1) then! origin
      res(j,l) = (cf0-2._8/dr**2) * f(j,l) + 4._8*f(j+1,l)/dr**2   &
                            + cfz*(f(j,l+1)+f(j,l-1)) &
                            + rhs(j,l)*inveps0
    else
      res(j,l) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                            + cfz*(f(j,l+1)+f(j,l-1)) &
                            + rhs(j,l)*inveps0
    END if
  enddo
else
 do l = nzi, nzf+1
  j = 1
  IF(bnd%v(j,l)==v_vacuum) &
  res(j,l) = (cf0-2._8/dr**2) * f(j,l) + 4._8*f(j+1,l)/dr**2   &
                                 + cfz*(f(j,l+1)+f(j,l-1)) &
                                 + rhs(j,l)*inveps0

  do j = 2, nrf+1
     IF(bnd%v(j,l)==v_vacuum) &
       res(j,l) = cf0 * f(j,l) + cfrp(j)*f(j+1,l)+cfrm(j)*f(j-1,l)   &
                                      + cfz*(f(j,l+1)+f(j,l-1)) &
                                      + rhs(j,l)*inveps0
  end do
 end do
END if

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    c => bnd%cndfirst
  else
    c => c%next
  END if
  do ii = 1, c%nbbnd
    j = c%jj(ii)
    l = c%kk(ii)
    IF(j==1) then
      IF(bnd%v(j,l)==v_bnd.and.c%docalc(ii)) &
      res(j,l) = c%cf0(ii)*f(j,l) &
                        + c%cfxp(ii)*f(j+1,l) &
                        + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
                        + voltfact*(c%phi0xp(ii) &
                        + c%phi0zp(ii)+c%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
    else
      IF(bnd%v(j,l)==v_bnd.and.c%docalc(ii)) &
      res(j,l) = c%cf0(ii)*f(j,l) &
                        + c%cfxp(ii)*f(j+1,l)+c%cfxm(ii)*f(j-1,l) &
                        + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
                        + voltfact*(c%phi0xp(ii)+c%phi0xm(ii) &
                        + c%phi0zp(ii)+c%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
    END if
  ENDDO
END do

IF(l_zerolastz) res(:,nz+1) = 0._8

IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  res(nr+1,:) = 0.
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  res(:,1) = 0.
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
  res(:,nz+1) = 0.
END if

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit resid, level = ',level
  call remark(trim(o_line))
END if

return
end subroutine residbndrzwguard

subroutine residbndrzwguard_list(res,jlocs,klocs,nvlocs,f,rhs,bnd,nr,nz,dr,dz,rmin,voltfact)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nr, nz, nvlocs
REAL(8), INTENT(IN) :: f(0:nr+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(nr+1,nz+1)
REAL(8), INTENT(IN OUT) :: res(nvlocs)
INTEGER(ISZ), INTENT(IN OUT), dimension(nvlocs) :: jlocs, klocs
TYPE(BNDtype) :: bnd
REAL(8), INTENT(IN) :: dr, dz, voltfact, rmin

INTEGER(ISZ) :: i, j, l, ii, ic, nrf, nzi, nzf
REAL(8) :: cf0, cfrp(nr+1), cfrm(nr+1), cfz, r
TYPE(CONDtype), pointer :: c

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter resid, level = ',level
  call remark(trim(o_line))
END if

cfz = 1._8 / dz**2
cf0 = -2._8/dr**2-2._8*cfz
do j = 2, nr+1
  r = rmin+REAL(j-1,8)*dr
  cfrp(j) = (1._8+0.5_8*dr/r) / dr**2
  cfrm(j) = (1._8-0.5_8*dr/r) / dr**2
!  cfrp(j) = (1._8+0.5_8/REAL(j-1,8)) / dr**2
!  cfrm(j) = (1._8-0.5_8/REAL(j-1,8)) / dr**2
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
    c => bnd%cndfirst
  else
    c => c%next
  END if
  do ii = 1, c%nbbnd
    i = i+1
    j = c%jj(ii)
    l = c%kk(ii)
    jlocs(i) = j
    klocs(i) = l
    IF(c%docalc(ii)) then
    IF(j==1) then
      res(i) = c%cf0(ii)*f(j,l) &
                        + c%cfxp(ii)*f(j+1,l) &
                        + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
                        + voltfact*(c%phi0xp(ii) &
                        + c%phi0zp(ii)+c%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
    else
      res(i) = c%cf0(ii)*f(j,l) &
                        + c%cfxp(ii)*f(j+1,l)+c%cfxm(ii)*f(j-1,l) &
                        + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
                        + voltfact*(c%phi0xp(ii)+c%phi0xm(ii) &
                        + c%phi0zp(ii)+c%phi0zm(ii)) &
                        + rhs(j,l)*inveps0
    END if
    else
      res(i) = 0.
    END if
  ENDDO
END do

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit resid, level = ',level
  call remark(trim(o_line))
END if

return
end subroutine residbndrzwguard_list

function residbndxzwguard(f,rhs,bnd,nx,nz,dx,dz,voltfact,l_zerolastz, ixlbnd, ixrbnd, izlbnd, izrbnd)
! evaluate residue. Grid is assumed to have guard cells, but residue does not.
implicit none

INTEGER(ISZ), INTENT(IN) :: nx, nz, ixlbnd, ixrbnd, izlbnd, izrbnd
REAL(8), INTENT(IN) :: f(0:,0:)!f(0:nx+2,0:nz+2)
REAL(8), INTENT(IN) :: rhs(:,:)!rhs(nx+1,nz+1)
TYPE(BNDtype) :: bnd
REAL(8), INTENT(IN) :: dx, dz,voltfact
REAL(8), DIMENSION(SIZE(f,1)-2,SIZE(f,2)-2) :: residbndxzwguard
LOGICAL(ISZ) :: l_zerolastz

INTEGER(ISZ) :: i, j, l, ii, ic, nxi, nxf, nzi, nzf
REAL(8) :: cf0, cfx, cfz
TYPE(CONDtype), POINTER :: c

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter resid, level = ',level
  call remark(trim(o_line))
END if

cfx = 1._8 / dx**2
cfz = 1._8 / dz**2
cf0 = -2._8*(cfx+cfz)

residbndxzwguard = 0._8

!do ic = 1, bnd%nb_conductors
!  IF(ic==1) then
!    bnd%cnd => bnd%cndfirst
!  else
!    bnd%cnd => c%next
!  END if
!  do i = 1, c%ncond
!    residbndxzwguard(c%jcond(i),c%kcond(i)) = 0._8
!  end do
!END do

IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  nxi=2
else
  nxi=1
END if
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  nxf=nx-1
else
  nxf=nx
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
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
    c => bnd%cndfirst
  else
    c => c%next
  END if
  do ii = 1, c%nbbnd
    j = c%jj(ii)
    l = c%kk(ii)
    IF(bnd%v(j,l)==v_bnd.and.c%docalc(ii)) &
    residbndxzwguard(j,l) = c%cf0(ii)*f(j,l) &
                          + c%cfxp(ii)*f(j+1,l)+c%cfxm(ii)*f(j-1,l) &
                          + c%cfzp(ii)*f(j,l+1)+c%cfzm(ii)*f(j,l-1) &
                          + voltfact*(c%phi0xp(ii)+c%phi0xm(ii) &
                          + c%phi0zp(ii)+c%phi0zm(ii)) &
                          + rhs(j,l)*inveps0
  ENDDO
END do

IF(l_zerolastz) residbndxzwguard(:,nz+1) = 0._8

IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) then
  residbndxzwguard(1,:) = 0.
END if
IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) then
  residbndxzwguard(nx+1,:) = 0.
END if
IF(izlbnd==dirichlet .or. izlbnd==patchbnd) then
  residbndxzwguard(:,1) = 0.
END if
IF(izrbnd==dirichlet .or. izrbnd==patchbnd) then
  residbndxzwguard(:,nz+1) = 0.
END if

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit resid, level = ',level
  call remark(trim(o_line))
END if

return
end function residbndxzwguard

!pgi$r nobounds
RECURSIVE subroutine mgbndrzwguard(j, u, rhs, bnd, nr, nz, dr, dz, rmin, npre, npost, ncycle, sub, relax_only, npmin, mgparam)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: j, nr, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(0:nr+2,0:nz+2), INTENT(IN OUT) :: u
REAL(8), DIMENSION(1:nr+1,1:nz+1), INTENT(IN) :: rhs
REAL(8) :: dr, dz, mgparam, rmin
TYPE(BNDtype), pointer :: bnd
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(:,:), allocatable :: res, v, ressub
INTEGER(ISZ) :: i,jj,ll
INTEGER :: nrnext, nznext, nzresmin, nzresmax, nzres
REAL(8) :: drnext, dznext, voltf

IF(.not.sub) inveps0 = 1./eps0

level = j
bndcurrent => bnd
IF(l_mgridrz_debug) then
  write(o_line,*) 'enter mg, level = ',level
  call remark(trim(o_line))
END if

IF(sub) then
  voltf = 0._8
else
  voltf = 1._8
END if

IF(j<=npmin .or. relax_only) then
  call apply_voltagewguard(u,bnd,voltf)
  call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
  IF(solvergeom==RZgeom) then
    call relaxbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,rmin=rmin,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  else ! solvergeom==XZgeom or solvergeom==XYgeom
    call relaxbndxzwguard(f=u,rhs=rhs,bnd=bnd,nx=nr,nz=nz,dx=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  END if
else
  nrnext = bnd%next%nr
  nznext = bnd%next%nz
  drnext = bnd%next%dr
  dznext = bnd%next%dz
    t_before = wtime()
  ALLOCATE(res(nrnext+1,nznext+1),v(0:nrnext+2,0:nznext+2))
    t_allocate = t_allocate + wtime()-t_before
  call apply_voltagewguard(u,bnd,voltf)
  call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
  IF(solvergeom==RZgeom) then
    call relaxbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,rmin=rmin,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  else ! solvergeom==XZgeom or solvergeom==XYgeom
    call relaxbndxzwguard(f=u,rhs=rhs,bnd=bnd,nx=nr,nz=nz,dx=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  END if
  IF(bnd%next%l_merged) then
#ifdef MPIPARALLEL
    IF(MOD(my_index/bnd%nworkpproc,2)==0) then
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
!                             residbndrzwguard(f=u,rhs=rhs,bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
!                             ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd), &
!                             bnd,nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
!                             ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
!  else
  IF(bnd%next%l_merged) res=0.
  res = 0.
  IF(solvergeom==RZgeom) then
    IF(vlocs) then
      call restrictlist(res(1,nzresmin), u, rhs, bnd, nrnext, nzres, nr, nz, voltf,dr,dz,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
    else
      allocate(ressub(nr+1,nz+1))
      call residbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd,nr=nr,nz=nz, &
              dr=dr,dz=dz,rmin=rmin,voltfact=voltf,l_zerolastz=.false., &
              ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd,res=ressub)
      call subrestrict(res(1,nzresmin), &
                             ressub, &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8,bnd%izlbnd,bnd%izrbnd)
      deallocate(ressub)
    END if
  else ! solvergeom==XZgeom or solvergeom==XYgeom
    res(:,nzresmin:nzresmax) = restrict( &
                             residbndxzwguard(f=u,rhs=rhs,bnd=bnd,nx=nr,nz=nz,dx=dr,dz=dz,voltfact=voltf,l_zerolastz=.false., &
                             ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd), &
                             nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8)
  END if
!  END if
  call apply_voltage(res,bnd%next,0._8)
#ifdef MPIPARALLEL
  IF(bnd%next%l_merged) call merge_work(res,level,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
  call exchange_resbndz(rho=res,izlbnd=bnd%next%izlbnd,izrbnd=bnd%next%izrbnd)
#endif
  v = 0.0_8
  IF(.not.sub) inveps0 = 1.
  do i = 1, ncycle  !(1=V cycles, 2=W cycle)
    call mgbndrzwguard(j=j-1, u=v(0,0), rhs=res(1,1), bnd=bnd%next,  &
                       nr=nrnext, nz=nznext, dr=drnext, dz=dznext, rmin=rmin, npre=npre, npost=npost, &
                       ncycle=ncycle, sub=.TRUE., relax_only=.false., npmin=npmin, mgparam=mgparam)
    level = j
  end do
  IF(.not.sub) inveps0 = 1./eps0
  call apply_voltagewguard(v,bnd%next,0._8)
  call updateguardcellsrz(f=v,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%next%izlbnd,izrbnd=bnd%next%izrbnd)
  call add_and_expand(u(0,0),v(0,nzresmin-1),bnd,nr,nz,nrnext,nzres,0._8,1._8,0._8,1._8,0._8,1._8,0._8,1._8, &
                                ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
! expand field from grid to a finer one. Each dimension may have any number of cells.
  call apply_voltagewguard(u,bnd,voltf)
  call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
#ifdef MPIPARALLEL
  if (.not. bnd%l_merged) then
    call exchange_fbndz(u,bnd%izlbnd,bnd%izrbnd)
  end if
#endif
  IF(solvergeom==RZgeom) then
    call relaxbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,rmin=rmin,nc=npost,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  else ! solvergeom==XZgeom or solvergeom==XYgeom
    call relaxbndxzwguard(f=u,rhs=rhs,bnd=bnd,nx=nr,nz=nz,dx=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  END if
  t_before = wtime()
  DEALLOCATE(res,v)
  t_allocate = t_allocate + wtime()-t_before
END if

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit mg, level = ',level
  call remark(trim(o_line))
END if

return
end subroutine mgbndrzwguard

subroutine mgbndrzwguard_jump(jmax, u, rhs, maxjump, bnd, nr, nz, dr, dz, &
    rmin, npre, npost, ncycle, sub, relax_only, npmin, mgparam)
! performs a multigrid cycle. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ), INTENT(IN) :: jmax, nr, nz, npre, npost, ncycle, npmin
REAL(8), DIMENSION(0:nr+2,0:nz+2), INTENT(IN OUT) :: u
REAL(8), DIMENSION(1:nr+1,1:nz+1), INTENT(IN) :: rhs
INTEGER(ISZ), DIMENSION(1:nr+1,1:nz+1), INTENT(IN) :: maxjump
REAL(8) :: dr, dz, mgparam, rmin
TYPE(BNDtype), pointer :: bnd
LOGICAL(ISZ), INTENT(IN) :: sub, relax_only

REAL(8), DIMENSION(1:nr+1,1:nz+1) :: res
REAL(8), DIMENSION(0:nr+2,0:nz+2) :: v
INTEGER(ISZ) :: i,jj,ll,j
REAL(8) :: voltf

level = nlevels
bndcurrent => bnd

IF(.not.sub) inveps0 = 1./eps0

IF(l_mgridrz_debug) then
  write(o_line,*) 'enter mg_jump'
  call remark(trim(o_line))
END if

voltf = 1._8
inveps0 = 1./eps0
call residbndrzwguard(f=u(0,0),rhs=rhs(1,1),bnd=bnd,nr=nr,nz=nz, &
                 dr=dr,dz=dz,rmin=rmin,voltfact=voltf,l_zerolastz=.false., &
                 ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd,res=res)
v = 0.
voltf = 0._8
inveps0 = 1.
do j = jmax, 1, -1
  call apply_voltagewguard(v,bnd,voltf)
  call updateguardcellsrz(f=v,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)
  IF(j==1) then
    call relaxbndrzwguard(f=v(0,0),rhs=res(1,1), &
                          bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,rmin=rmin,nc=npre,voltfact=voltf,mgparam=mgparam, &
                          ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
!    call relaxbndrzwguard_jump(f=v(0,0),rhs=res(1,1), maxjump=maxjump(1,1), curjump=j, &
!                               bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,nc=npre,voltfact=voltf,mgparam=mgparam, &
!                               ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  else
    call relaxbndrzwguard_jump(f=v(0,0),rhs=res(1,1), maxjump=maxjump(1,1), curjump=j, &
                               bnd=bnd,nr=nr,nz=nz,dr=dr,dz=dz,rmin=rmin,nc=npre,voltfact=voltf,mgparam=mgparam, &
                               ixlbnd=ixlbnd, ixrbnd=ixrbnd, izlbnd=bnd%izlbnd, izrbnd=bnd%izrbnd)
  END if
end do
u = u + v
voltf = 1._8
inveps0 = 1./eps0
call apply_voltagewguard(u,bnd,voltf)
call updateguardcellsrz(f=u,ixlbnd=ixlbnd,ixrbnd=ixrbnd,izlbnd=bnd%izlbnd,izrbnd=bnd%izrbnd)

IF(l_mgridrz_debug) then
  write(o_line,*) 'exit mg_jump'
  call remark(trim(o_line))
END if

return
end subroutine mgbndrzwguard_jump

!pgi$r nobounds
subroutine updateguardcellsrz(f, ixlbnd, ixrbnd, izlbnd, izrbnd)
! update guard cells values according to boundary conditions.
implicit none
REAL(8),INTENT(IN OUT) :: f(0:,0:)
INTEGER(ISZ), optional :: ixlbnd
INTEGER(ISZ) :: ixrbnd, izlbnd, izrbnd

INTEGER(ISZ) :: ixmax, izmax

t_before = wtime()

ixmax=SIZE(f,1)
izmax=SIZE(f,2)

IF(PRESENT(ixlbnd)) then
  select case (ixlbnd)
    case (dirichlet)
      f(0,:) = 2.*f(1,:)-f(2,:)
    case (neumann)
      f(0,:) = f(2,:)
    case (periodic)
      f(0,:) = f(ixmax-2,:)
    case default
  end select
END if
select case (ixrbnd)
    case (dirichlet)
      f(ixmax-1,:) = 2.*f(ixmax-2,:)-f(ixmax-3,:)
    case (neumann)
      f(ixmax-1,:) = f(ixmax-3,:)
    case (periodic)
      f(ixmax-1,:) = f(1,:) 
    case default
end select
select case (izlbnd)
    case (dirichlet)
      f(:,0) = 2.*f(:,1)-f(:,2)
    case (neumann)
      f(:,0) = f(:,2)
    case (periodic)
      f(:,0) = f(:,izmax-2)
    case default
end select
select case (izrbnd)
    case (dirichlet)
      f(:,izmax-1) = 2.*f(:,izmax-2)-f(:,izmax-3)
    case (neumann)
      f(:,izmax-1) = f(:,izmax-3)
    case (periodic)
      f(:,izmax-1) = f(:,1)
    case default
end select

t_updateguard = t_updateguard + wtime()-t_before
end subroutine updateguardcellsrz

subroutine apply_voltage(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors
implicit none
REAL(8),INTENT(IN OUT) :: f(:,:)
TYPE(BNDtype) :: bnd
REAL(8), INTENT(IN) :: coef_voltage

INTEGER(ISZ) :: ic, i
TYPE(CONDtype), pointer :: c

return
t_before = wtime()
do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    c => bnd%cndfirst
  else
    c => c%next
  END if
    do i = 1, c%ncond
      f(c%jcond(i),c%kcond(i)) = coef_voltage*c%voltage(i)
  end do
END do

#ifdef MPIPARALLEL
!  call exchange_fbndz(f,izlbnd,izrbnd)
#endif

t_apply_voltage = t_apply_voltage + wtime()-t_before
return
end subroutine apply_voltage

!pgi$r nobounds
subroutine apply_voltagewguard(f,bnd,coef_voltage)
! assign voltage value at grid nodes located inside conductors. Grid is assumed to have guard cells.
implicit none
REAL(8),INTENT(IN OUT) :: f(0:,0:)
TYPE(BNDtype) :: bnd
REAL(8), INTENT(IN) :: coef_voltage

INTEGER(ISZ) :: ic, i
TYPE(CONDtype), pointer :: c
t_before = wtime()

do ic = 1, bnd%nb_conductors
  IF(ic==1) then
    c => bnd%cndfirst
  else
    c => c%next
  END if
    do i = 1, c%ncond
      f(c%jcond(i),c%kcond(i)) = coef_voltage*c%voltage(i)
  end do
END do

t_apply_voltage = t_apply_voltage + wtime()-t_before
return
end subroutine apply_voltagewguard

subroutine solve_multigridrz(grid,accuracy,l_for_timing)
! solve field for u with density rhoinit.
implicit none

! input/output variables
TYPE(GRIDtype) :: grid
REAL(8), INTENT(IN) :: accuracy  ! required average accuracy
LOGICAL(ISZ) :: l_for_timing

INTEGER(ISZ) :: i, ii, ic, j, l, nr, nz, nlocs
REAL(8), allocatable, DIMENSION(:,:) :: uold, uinit
REAL(8), allocatable, DIMENSION(:) :: uold_vlocs, uinit_vlocs
INTEGER(ISZ), allocatable, DIMENSION(:) :: jlocs, klocs
REAL(8) :: maxerr, maxerr_old, t_solve
LOGICAL :: do_calc, has_diverged
TYPE(CONDtype), POINTER :: c

nlevels = grid%nlevels
level = nlevels
ixlbnd = grid%ixlbnd
ixrbnd = grid%ixrbnd
izlbnd = grid%izlbnd
izrbnd = grid%izrbnd

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

nr = SIZE(grid%phi,1)
nz = SIZE(grid%phi,2)

IF(vlocs) then
  nlocs = grid%bndfirst%nvlocs
  do ic = 1, grid%bndfirst%nb_conductors
    IF(ic==1) then
      c => grid%bndfirst%cndfirst
    else
      c => c%next
    END if
    nlocs = nlocs + c%nbbnd
  END do
  ALLOCATE(uold_vlocs(nlocs),jlocs(nlocs),klocs(nlocs))
  do i = 1, grid%bndfirst%nvlocs
    j = grid%bndfirst%vlocs_j(i)
    l = grid%bndfirst%vlocs_k(i)
    jlocs(i) = j
    klocs(i) = l
  enddo
  i = grid%bndfirst%nvlocs
  do ic = 1, grid%bndfirst%nb_conductors
    IF(ic==1) then
      c => grid%bndfirst%cndfirst
    else
      c => c%next
    END if
    do ii = 1, c%nbbnd
      i = i+1
      j = c%jj(ii)
      l = c%kk(ii)
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
!  IF(l_for_timing) then
    ALLOCATE(uinit(nr,nz))
    uinit = grid%phi
!  END if
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
      call mgbndrzwguard(j=nlevels,u=grid%phi,rhs=grid%rho,bnd=grid%bndfirst,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
                         rmin=grid%rmin, &
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
      IF(maxerr/maxerr_old>1..and.j>1) then
#ifdef MPIPARALLEL
       IF(my_index==0) then
#endif
        write(o_line,*) 'WARNING multigridrz, calculation is diverging:';   call remark(trim(o_line))
        write(o_line,*) '        initial maximum error = ',maxerr_old;      call remark(trim(o_line))
        write(o_line,*) '        current maximum error = ',maxerr;          call remark(trim(o_line))
        IF(grid%npre<10) then
          write(o_line,*) '        trying npre and npost = ',grid%npre+1;     call remark(trim(o_line))!,' (also reset mgparam to 1.8)'
        else if(grid%npmin<grid%nlevels) then
          write(o_line,*) '        trying npmin+1/nlevels = ',grid%npmin+1,grid%nlevels;     call remark(trim(o_line))!,' (also reset mgparam to 1.8)'
        END if
#ifdef MPIPARALLEL
      END if
#endif
        IF(grid%npre<10) then
          grid%npre  = grid%npre+1
          grid%npost = grid%npost+1
        else if(grid%npmin<grid%nlevels) then
          grid%npmin=grid%npmin+1
        else
          WRITE(0,*) 'Convergence cannot be achieved on grid ',grid%gid(1)
          stop
        END if
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
            grid%phi=uinit!uold
          END if
        END if
!        IF(l_for_timing) exit
        exit
      END if
    end do
    IF(j>=grid%ncmax) do_calc=.false.
  end do

t_solve = wtime() - t_solve
#ifdef MPIPARALLEL
!IF(my_index==0) then
#endif
  IF(lverbose>=1 .and. grid%gid(1)==basegrid%gid(1)) then
    write(o_line,'("multigridrz: precision = ",e12.5, " after ",i5," iterations.")') maxerr,j; call remark(trim(o_line))
  END if
  IF(lverbose>=2 .and. grid%gid(1)/=basegrid%gid(1)) then
    write(o_line,'("multigridrz: precision = ",e12.5, " after ",i5," iterations.")') maxerr,j; call remark(trim(o_line))
    WRITE(o_line,'("grid Nb: ",i5,", Nr = ",i5,", Nz = ",i5)') grid%gid(1),grid%nr,grid%nz; call remark(trim(o_line))
  END if
nb_iters=j
IF(j<grid%ncmax.and..not.do_calc.and.maxerr >= accuracy) nb_iters=grid%ncmax

IF(vlocs) then
  DEALLOCATE(uold_vlocs,jlocs,klocs)
  IF(l_for_timing) DEALLOCATE(uinit_vlocs)
else
  DEALLOCATE(uold)
!  IF(l_for_timing) DEALLOCATE(uinit)
  DEALLOCATE(uinit)
END if
IF(.not.l_for_timing) then
  IF(l_print_timing) then
    write(o_line,*) 'Time relax    = ',t_relax;              call remark(trim(o_line))
    write(o_line,*) 'Time restrict = ',t_restrict;           call remark(trim(o_line))
    write(o_line,*) 'Time expand   = ',t_expand;             call remark(trim(o_line))
    write(o_line,*) 'Time apply_voltage = ',t_apply_voltage; call remark(trim(o_line))
    write(o_line,*) 'Time updateguard   = ',t_updateguard;   call remark(trim(o_line))
    write(o_line,*) 'Time allocate      = ',t_allocate;      call remark(trim(o_line))
    write(o_line,*) 'Time total    = ',t_relax+t_restrict+t_expand+t_apply_voltage+t_updateguard+t_allocate
    call remark(trim(o_line))
    write(o_line,*) 'Time solve    = ',t_solve;              call remark(trim(o_line))
  END if
END if

return
end subroutine solve_multigridrz

subroutine solve_multigridrz_jump(grid,accuracy,l_for_timing)
! solve field for u with density rhoinit.
implicit none

! input/output variables
TYPE(GRIDtype) :: grid
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
               grid%bndfirst%v(jp,l)==v_vacuum .AND. grid%bndfirst%v(jm,l)==v_vacuum .AND. &
               grid%bndfirst%v(j,lp)==v_vacuum .AND. grid%bndfirst%v(j,lm)==v_vacuum)
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
  write(o_line,*) 'jmax = ',jmax
  call remark(trim(o_line))
  do while(do_calc)
    maxerr = 1.
    do  j = 1, mgridrz_ncmax
      uold=grid%phi
!      call mgbndrzwguard(j=nlevels,u=grid%phi(0,0),rhs=grid%rho(1,1),bnd=grid%bnd,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
      call mgbndrzwguard_jump(jmax=jmax,u=grid%phi,rhs=grid%rho,maxjump=maxjump, &
                         rmin=grid%rmin, &
                         bnd=grid%bndfirst,nr=grid%nr,nz=grid%nz,dr=grid%dr,dz=grid%dz, &
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
        write(o_line,*) 'WARNING multigridrz, calculation is diverging:'; call remark(trim(o_line))
        write(o_line,*) '        initial maximum error = ',maxerr_old   ; call remark(trim(o_line))
        write(o_line,*) '        current maximum error = ',maxerr       ; call remark(trim(o_line))
        write(o_line,*) '        trying npre and npost = ',grid%npre+1  ; call remark(trim(o_line)) !,' (also reset mgparam to 1.8)'

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
  write(o_line,'("multigridrz: precision = ",e12.5, " after ",i5," iterations.")') maxerr,j
  call remark(trim(o_line))
#else
  write(o_line,'("multigridrz: precision = ",e12.5, " after ",i5," iterations.")') maxerr,j
  call remark(trim(o_line))
#endif
nb_iters=j
IF(j<grid%ncmax.and..not.do_calc.and.maxerr >= accuracy) nb_iters=grid%ncmax
!IF(levels_min/=npmin) then
!  write(o_line,'("WARNING multigridrz, levels_min = ",i2,", npmin = ",i2,". Setting levels_min = ",i2,".")') levels_min, npmin, npmin
!  levels_min = npmin
!END if

DEALLOCATE(uold,maxjump)
IF(l_for_timing) DEALLOCATE(uinit)

IF(.not.l_for_timing) then
  IF(l_print_timing) then
    write(o_line,*) 'Time relax    = ',t_relax
    call remark(trim(o_line))
    write(o_line,*) 'Time restrict = ',t_restrict
    call remark(trim(o_line))
    write(o_line,*) 'Time expand   = ',t_expand
    call remark(trim(o_line))
    write(o_line,*) 'Time apply_voltage = ',t_apply_voltage
    call remark(trim(o_line))
    write(o_line,*) 'Time updateguard   = ',t_updateguard
    call remark(trim(o_line))
    write(o_line,*) 'Time allocate      = ',t_allocate
    call remark(trim(o_line))
    write(o_line,*) 'Time total    = ',t_relax+t_restrict+t_expand+t_apply_voltage+t_updateguard+t_allocate
    call remark(trim(o_line))
    write(o_line,*) 'Time solve    = ',t_solve
    call remark(trim(o_line))
  END if
END if

return
end subroutine solve_multigridrz_jump

subroutine solve_multigridz(grid)
! solve potential phi from density rho on a 1-d regular grid.
implicit none

! input/output variables
TYPE(GRIDtype) :: grid

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

subroutine solve_multigridr(grid)
! solve potential phi from density rho on a 1-d regular radial grid.
! solution taken from Langdon and Birdsall, chap. 14-10
implicit none

! input/output variables
TYPE(GRIDtype) :: grid

INTEGER(ISZ) :: j
REAL(8) :: rj, er(1:grid%nr)

! first, compute electric field
er(1) = 0.25*grid%dr*grid%rho(1,1)*inveps0
do j = 2, grid%nr
  rj = REAL(j,8)
  er(j) = ((rj-0.5)*er(j-1)+j*grid%dr*grid%rho(j,1)*inveps0 )/(rj+0.5)
end do

! then, compute potential
do j = grid%nr, 1, -1
  grid%phi(j,1) = grid%phi(j+1,1)+grid%dr*er(j)
end do
grid%phi(0,1) = grid%phi(2,1)
grid%phi(grid%nr+2,1) = 2.*grid%phi(grid%nr+1,1)-grid%phi(grid%nr,1)

end subroutine solve_multigridr

recursive subroutine find_mgparam_rz_1grid(grid,lsavephi,l_gonext,l_godown)
implicit none
TYPE(GRIDtype):: grid
logical(ISZ):: lsavephi, l_godown, l_gonext
REAL(8) :: nexttime, prevtime, prevparam
INTEGER(ISZ) :: npreinit, npostinit
real(8),allocatable :: phisave(:,:)

  nlevels = grid%nlevels
  level = nlevels
  ixrbnd = grid%ixrbnd
  izlbnd = grid%izlbnd
  izrbnd = grid%izrbnd

  IF(grid%gid(1)/=basegrid%gid(1)) then
        call setphigridrz(grid%phi,                               &
                          grid%rmin-grid%nguardx*grid%dr, &
                          grid%zmin-grid%nguardz*grid%dz,         &
                          grid%dr,grid%dz,                        &
                          grid%nr+2*grid%nguardx,              &
                          grid%nz+2*grid%nguardz,.true.)
  END if

  if (lsavephi) then
    allocate(phisave(size(grid%phi,1),size(grid%phi,2)))
    phisave = grid%phi
  else
    allocate(phisave(1,1))
  endif

  npreinit = grid%npre
  npostinit = grid%npost
  ! --- Get initial field solve time
  nexttime = time_field_solve(grid,lsavephi,phisave)
  prevtime = 2*nexttime
  ! --- Loop, increasing the number of passes until the time is minimized.
  DO WHILE(nexttime < prevtime)
    prevparam = grid%mgparam
    prevtime = nexttime
    nexttime = ffind_mgparam(grid,lsavephi,phisave)
    IF(my_index==0) then
      write(o_line,*) "Field solve time = ",nexttime
      call remark(trim(o_line))
      write(o_line,*) "mgparam = ",grid%mgparam
      call remark(trim(o_line))
      write(o_line,*) "npre    = ",grid%npre
      call remark(trim(o_line))
      write(o_line,*) "npost   = ",grid%npost
      call remark(trim(o_line))
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
      write(o_line,*) 'Notice: the maximum number of iterations has been reached, so '
      call remark(trim(o_line))
      write(o_line,*) 'the values above are unlikely to be optimal. Try increasing the '
      call remark(trim(o_line))
      write(o_line,*) 'tolerance, increasing the maximum number of iterations, or making a '
      call remark(trim(o_line))
      write(o_line,*) 'better initial guess of mgparam.'
      call remark(trim(o_line))
    END if
  else
    prevtime=findnrecursmin(grid,prevtime,lsavephi,phisave)
    IF(my_index==0) then
      write(o_line,*) "-----------------------------------------"
      call remark(trim(o_line))
      write(o_line,*) "The optimized values:"
      call remark(trim(o_line))
      write(o_line,*) "Field solve time = ",prevtime
      call remark(trim(o_line))
      write(o_line,*) "frz.mgridrz_mgparam     = ",grid%mgparam
      call remark(trim(o_line))
      write(o_line,*) "frz.mgridrz_npre        = ",grid%npre
      call remark(trim(o_line))
      write(o_line,*) "frz.mgridrz_npost       = ",grid%npost
      call remark(trim(o_line))
      write(o_line,*) "frz.mgridrz_levels_min  = ",grid%npmin
      call remark(trim(o_line))
    END if
  END if

  deallocate(phisave)

  IF(associated(grid%next) .and. l_gonext) then
    call find_mgparam_rz_1grid(grid%next,lsavephi,.true.,.false.)
  END if
  IF(associated(grid%down) .and. l_godown) then
    call find_mgparam_rz_1grid(grid%down,lsavephi,.true.,.true.)
  END if

return
END subroutine find_mgparam_rz_1grid

function time_field_solve(grid,lsavephi,phisave)
implicit none
REAL(8) :: time_field_solve
TYPE(GRIDtype):: grid
logical(ISZ):: lsavephi
real(8):: phisave(:,:)

INTEGER(ISZ) :: ixmax, izmin, izmax
REAL(8) :: beforetime, aftertime
REAL(8), EXTERNAL :: wtime

  ixmax = grid%nr+1
  izmin = 1
  izmax = grid%nz+1
  IF(grid%ixrbnd==dirichlet .or. grid%ixrbnd==patchbnd) ixmax=ixmax-1
  IF(grid%izlbnd==dirichlet .or. grid%izlbnd==patchbnd) izmin=izmin+1
  IF(grid%izrbnd==dirichlet .or. grid%izrbnd==patchbnd) izmax=izmax-1

  if (lsavephi) then
    grid%phi = phisave
  else
    grid%phi(1:ixmax,izmin:izmax)=0.
  endif

  beforetime = wtime()
  call solve_multigridrz(grid=grid, accuracy=mgridrz_accuracy, l_for_timing=.true.)
  aftertime = wtime()
  time_field_solve = aftertime - beforetime
#ifdef MPIPARALLEL
  time_field_solve = mpi_global_compute_real(time_field_solve,MPI_MAX)
#endif

  return
END function time_field_solve

function ffind_mgparam(grid,lsavephi,phisave)
implicit none
REAL(8) :: ffind_mgparam
TYPE(GRIDtype):: grid
logical(ISZ):: lsavephi
real(8):: phisave(:,:)

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
  ffind_mgparam = time_field_solve(grid,lsavephi,phisave)

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
    IF(my_index==0) then
      write(o_line,*) "Best parameter so far = ", grid%mgparam
      call remark(trim(o_line))
    END if

!   --- do field solve (which prints out number of field solve iterations)
    up_old = grid%npre
    down_old = grid%npost
    ffind_mgparam = time_field_solve(grid,lsavephi,phisave)

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
      IF(my_index==0) then
        write(o_line,*) "resetting ffind_mgparam"
        call remark(trim(o_line))
      END if
      icount=0
      grid%npre = up_old
      grid%npost = down_old
      ffind_mgparam = time_field_solve(grid,lsavephi,phisave)
      mgparam_prev = mgparam_init
      mgiters_prev = nb_iters
      sincr = .05
!      nb_iters = 0
      grid%mgparam = grid%mgparam + sincr
    END if

  END do

! --- write(o_line,*) message if an optimal value wasn't found
  if (icount == 200) then
    IF(my_index==0) then
      write(o_line,*) "Warning: maximum number of iterations reached."
      call remark(trim(o_line))
      write(o_line,*) "         The value of mgparam may not be optimal."
      call remark(trim(o_line))
      write(o_line,*) "         Try increasing mgmaxit."
      call remark(trim(o_line))
    END if
  END if

  return
END function ffind_mgparam

function findnrecursmin(grid,prevtime,lsavephi,phisave)
!Optimize levels_min, minimizing the fieldsolve time.
implicit none
REAL(8) :: findnrecursmin
TYPE(GRIDtype) :: grid
REAL(8), INTENT(IN) :: prevtime
logical(ISZ):: lsavephi
real(8):: phisave(:,:)

REAL(8) :: nexttime, prvtime

  ! --- Get initial field solve time
  nexttime = prevtime
  prvtime = 2*nexttime
  ! --- Loop, increasing the number of passes until the time is minimized.
  do WHILE(nexttime < prvtime .and. grid%npmin < mgridrz_nlevels_max)
    prvtime = nexttime
    grid%npmin = grid%npmin + 1
    nexttime = time_field_solve(grid,lsavephi,phisave)
    IF(my_index==0) then
      write(o_line,*) "Field solve time = ",nexttime
      call remark(trim(o_line))
      write(o_line,*) "frz.mgridrz_levels_min = ",grid%npmin
      call remark(trim(o_line))
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

subroutine multigridrzf(iwhich,u0,rho0,nr0,nz0)
USE InjectVars_eq
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich, nr0, nz0
REAL(8), INTENT(IN OUT) :: u0(1:nr0+1,0:nz0+2)
REAL(8), INTENT(IN OUT) :: rho0(nr0+1,nz0+1)

  IF(mgridrz_ncmax==0) return

  IF(iwhich==1) return

  IF(l_find_rise_time) then
    call multigridrzf_risetime(iwhich,u0,rho0,nr0,nz0,mgridrz_accuracy)
    return
  END if

!  call distribute_rho(basegrid)
  IF(solvergeom==XZgeom) then
    if (basegrid%ixlbnd==dirichlet .or. basegrid%ixlbnd==patchbnd) basegrid%phi(1,:)           = u0(1,:)
  END if
  IF(solvergeom==RZgeom .OR. solvergeom==XZgeom) then
    if (basegrid%ixrbnd==dirichlet .or. basegrid%ixrbnd==patchbnd) basegrid%phi(nr0+1,:)       = u0(nr0+1,:)
    if (basegrid%izlbnd==dirichlet .or. basegrid%izlbnd==patchbnd) basegrid%phi(1:nr0+1,1)     = u0(1:nr0+1,1)
    if (basegrid%izrbnd==dirichlet .or. basegrid%izrbnd==patchbnd) basegrid%phi(1:nr0+1,nz0+1) = u0(1:nr0+1,nz0+1)
  END if
  IF(solvergeom==Zgeom) then
    basegrid%phi(1,1)     = u0(1,1)
    basegrid%phi(1,nz0+1) = u0(1,nz0+1)
  END if
  IF(solvergeom==Rgeom) then
    basegrid%phi(nr0+1,1) = u0(nr0+1,1)
  END if

  call solve_mgridrz(basegrid,mgridrz_accuracy,.true.)

  u0(1:nr0+1,1:nz0+1)=basegrid%phi(1:nr0+1,1:nz0+1)

return
end subroutine multigridrzf

subroutine multigridxyf2(iwhich,u0,rho0,nx0,ny0)
USE InjectVars_eq
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: iwhich, nx0, ny0
REAL(8), INTENT(IN OUT) :: u0(nx0+1,ny0+1)
REAL(8), INTENT(IN OUT) :: rho0(nx0+1,ny0+1)

  IF(mgridrz_ncmax==0) return

  IF(iwhich==1) return
!  call distribute_rho(basegrid)
  if (basegrid%ixlbnd==dirichlet .or. basegrid%ixlbnd==patchbnd) basegrid%phi(1,1:ny0+1)     = u0(1,:)
  if (basegrid%ixrbnd==dirichlet .or. basegrid%ixrbnd==patchbnd) basegrid%phi(nx0+1,1:ny0+1) = u0(nx0+1,:)
  if (basegrid%izlbnd==dirichlet .or. basegrid%izlbnd==patchbnd) basegrid%phi(1:nx0+1,1)     = u0(:,1)
  if (basegrid%izrbnd==dirichlet .or. basegrid%izrbnd==patchbnd) basegrid%phi(1:nx0+1,ny0+1) = u0(:,ny0+1)

  call solve_mgridrz(basegrid,mgridrz_accuracy,.true.)

  u0(1:nx0+1,1:ny0+1)=basegrid%phi(1:nx0+1,1:ny0+1)

return
end subroutine multigridxyf2

subroutine multigridrzf_risetime(iwhich,u0,rho0,nr0,nz0,accuracy)
USE InGen
USE InPart
USE InMesh3d
USE Particles
USE InjectVars
USE InjectVars3d
USE InjectVars_eq, ONLY: inj_phi_eq,v_max,afact,calc_a
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
    write(o_line,*) 'ninject>1 not supported by multigridrzf_risetime, stopping.'
    call remark(trim(o_line))
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

  vinject = v_max
  nps_tmp = nps
  nps = 0
  call getinj_phi()
  nps = nps_tmp
  select case (calc_a)
    case (center)
      phiv = inj_phi(0,0,1)
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy,.true.)
      vinject=0.
      call getinj_phi()
      phirho = inj_phi(0,0,1)
!      IF(inject==3) then
!        phirho = phirho + inj_dz*abs(inj_d(1))*0.5*inj_rho(0,0,1)*inj_dz*inveps0
!        phiref = phiref - inj_dz*abs(inj_d(1))*0.5*inj_rho_eq(0,0,1)*inj_dz*inveps0
!      END if
    case (average_source)
      max_j = 1+INT(ainject(1)/inj_dx)
      phiv = SUM(inj_phi(0:max_j-1,0,1))/max_j
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy,.true.)
      vinject=0.
      call getinj_phi()
      phirho = SUM(inj_phi(0:max_j-1,0,1))/max_j
    case (weighted_average_source)
      max_j = 1+INT(ainject(1)/inj_dx)
      ALLOCATE(weights(max_j))
      weights(1) = 0.25*pi*inj_dx**2
      do j = 2, max_j
        weights(j) = 2.*pi*(j-1)*inj_dx**2
      end do
      wtot = SUM(weights(1:max_j))
      phiv = SUM(weights(1:max_j)*inj_phi(0:max_j-1,0,1))/wtot
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy,.true.)
      vinject=0.
      call getinj_phi()
      phirho = SUM(weights(1:max_j)*inj_phi(0:max_j-1,0,1))/wtot
      DEALLOCATE(weights)
    case (border)
      max_j = 1+INT(ainject(1)/inj_dx)
      phiv = inj_phi(max_j-2,0,1)
      do i = 1, ngrids
        grids_ptr(i)%grid%phi=0.
      END do
      call solve_mgridrz(basegrid,accuracy,.true.)
      vinject=0.
      call getinj_phi()
      phirho = inj_phi(max_j-2,0,1)
    case default
  end select

  afact = (phiref-phirho)/phiv

!  basegrid%phi(1:nr0+1,:) = basegrid%phi(1:nr0+1,:) + a*phi_init(:,1,:)
  do i = 1, ngrids
    grids_ptr(i)%grid%phi = grids_ptr(i)%grid%phi + afact*gridinit(i)%grid%phi
  end do
  if(solvergeom==RZgeom) &
    call updateguardcellsrz(basegrid%phi, basegrid%ixlbnd, basegrid%ixrbnd, basegrid%izlbnd, basegrid%izrbnd)

  vinject=afact*phi0
  call getinj_phi()
  IF(lverbose>=3) then
    write(o_line,*) 'a,inj_phi,inj_phi_eq',afact,inj_phi(0,0,1),inj_phi_eq
    call remark(trim(o_line))
    write(o_line,*) 'phi0, phiv, phiref, phirho',phi0, phiv, phiref, phirho
    call remark(trim(o_line))
  END if

  u0(1:nr0+1,:)=basegrid%phi(1:nr0+1,:)

return
end subroutine multigridrzf_risetime

subroutine distribute_rho_rz()
USE multigridrz
implicit none

IF(.not.l_distribute) return
IF(solvergeom==Zgeom .or. solvergeom==Rgeom) then
  call distribute_rho(basegrid)
else
  call children_send_rho_to_parents()
END if
call exchange_rho_between_neighbors()

return
END subroutine distribute_rho_rz

RECURSIVE subroutine distribute_rho(grid)
USE multigridrz
implicit none
TYPE(GRIDtype) :: grid

  IF(associated(grid%down)) call distribute_rho(grid%down)
  IF(associated(grid%next)) call distribute_rho(grid%next)
  IF(associated(grid%up)) then
    IF(solvergeom==Zgeom) then
      call deposit_z(unew=grid%up%rho(1,:), uold=grid%rho(1,:), invvolnew=1./grid%up%dz, invvolold=1./grid%dz, &
                     zminold=grid%zmin, zmaxold=grid%zmax, zminnew=grid%up%zmin, zmaxnew=grid%up%zmax)
    else IF(solvergeom==RZgeom) then
      call deposit_rz(unew=grid%up%rho, uold=grid%rho, invvolnew=grid%up%invvol, invvolold=grid%invvol, &
                   xminold=grid%rmin, xmaxold=grid%rmax, zminold=grid%zmin, zmaxold=grid%zmax, &
                   xminnew=grid%up%rmin, xmaxnew=grid%up%rmax, zminnew=grid%up%zmin, zmaxnew=grid%up%zmax)
    else
      call deposit(unew=grid%up%rho, uold=grid%rho, &
                   xminold=grid%rmin, xmaxold=grid%rmax, zminold=grid%zmin, zmaxold=grid%zmax, &
                   xminnew=grid%up%rmin, xmaxnew=grid%up%rmax, zminnew=grid%up%zmin, zmaxnew=grid%up%zmax)
    END if
  END if

return
END subroutine distribute_rho

RECURSIVE subroutine solve_mgridrz(grid,accuracy,fromup)
USE multigridrz
implicit none
TYPE(GRIDtype) :: grid
REAL(8), INTENT(IN) :: accuracy
LOGICAL :: fromup

TYPE(BNDtype), POINTER :: bnd
TYPE(CONDtype), POINTER :: c

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
        call setphigridrz(grid%phi,                               &
                          grid%rmin-grid%nguardx*grid%dr, &
                          grid%zmin-grid%nguardz*grid%dz,         &
                          grid%dr,grid%dz,                        &
                          grid%nr+2*grid%nguardx,              &
                          grid%nz+2*grid%nguardz,.false.)
!        CALL interpolate_any(unew=grid%phi,uold=grid%up%phi, &
!                             nxnew=grid%nr, nznew=grid%nz, &
!                             nxold=grid%up%nr, nzold=grid%up%nz, &
!                             xminold=grid%up%rmin, xmaxold=grid%up%rmax, &
!                             zminold=grid%up%zmin, zmaxold=grid%up%zmax, &
!                             xminnew=grid%rmin, xmaxnew=grid%rmax, &
!                             zminnew=grid%zmin, zmaxnew=grid%zmax, &
!                             ixrbnd=grid%ixrbnd, &
!                             izlbnd=grid%izlbnd, &
!                             izrbnd=grid%izrbnd, &
!!                             bnd_only=.false., quad=.false.)
!                             bnd_only=.true., quad=.false.)
        bnd => grid%bndfirst
        do ic = 1, bnd%nb_conductors
          IF(ic==1) then
            c => bnd%cndfirst
          else
            c => c%next
          END if
          do i = 1, c%ncond
            IF(c%condid(i)<0) c%voltage(i) = grid%phi(c%jcond(i),c%kcond(i))
          end do
        end do
      END if
    END if
!    IF(.not. associated(grid%up)) call solve_multigridrz(grid=grid, accuracy=accuracy, l_for_timing=.false.)
    IF(solvergeom==Rgeom) then
      call solve_multigridr(grid=grid)
    ELSE IF(solvergeom==Zgeom) then
      call solve_multigridz(grid=grid)
    ELSE IF(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==XYgeom) then
      call solve_multigridrz(grid=grid, accuracy=accuracy, l_for_timing=.false.)
    END IF

    IF(associated(grid%next)) call solve_mgridrz(grid%next,accuracy,.false.)
    IF(fromup .and. associated(grid%down)) call solve_mgridrz(grid%down,accuracy,.true.)

return
END subroutine solve_mgridrz

subroutine find_mgparam_rz(lsavephi)
USE multigridrz
implicit none
logical(ISZ):: lsavephi

    call find_mgparam_rz_1grid(grid=basegrid,lsavephi=lsavephi,l_gonext=.TRUE.,l_godown=.TRUE.)

return
END subroutine find_mgparam_rz

subroutine find_mgparam_rz_1g(grid)
USE multigridrz
implicit none
TYPE(GRIDtype)::grid

 call find_mgparam_rz_1grid(grid=grid,lsavephi=.false.,l_gonext=.FALSE.,l_godown=.FALSE.)

return
END subroutine find_mgparam_rz_1g

subroutine srfrvoutrz(rofzfunc,volt,zmin,zmax,xcent,rmax,lfill,  &
                      xmin,xmax,lshell,                      &
                      zmmin,zmmax,zbeam,dx,dz,nx,nz,             &
                      ix_axis,xmesh,l2symtry_in,l4symtry_in,condid)
! call subroutine srfrvout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use Conductor3d
USE multigridrz
implicit none
character(*) rofzfunc
real(kind=8):: volt,zmin,zmax,xcent,rmax
LOGICAL(ISZ):: lfill,lshell,l2symtry_in,l4symtry_in
real(kind=8):: xmin,xmax,ymax
real(kind=8):: zmmin,zmmax,zbeam,dx,dz
INTEGER(ISZ):: nx,nz,ix_axis
real(kind=8):: xmesh(0:nx)
integer(ISZ):: condid

INTEGER(ISZ) :: i,nrc,nzc,igrid,jmin,jmax
REAL(8) :: drc,dzc,zmin_in,zmax_in
LOGICAL(ISZ) :: doloop

TYPE(CONDtype), POINTER :: c
TYPE(BNDtype), POINTER :: b
IF(solvergeom==Zgeom .or. solvergeom==Rgeom) return

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmin,zmmin)

do igrid=1,ngrids
 IF(zmin>grids_ptr(igrid)%grid%zmax .OR. zmax<grids_ptr(igrid)%grid%zmin .OR. &
    rmax<grids_ptr(igrid)%grid%rmin) cycle
 nlevels=grids_ptr(igrid)%grid%nlevels
 ixlbnd = grids_ptr(igrid)%grid%ixlbnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = 1, nlevels
  IF(i == 1) then
    b => grids_ptr(igrid)%grid%bndfirst
  else
    b => b%next
  END if
  nrc = b%nr
  nzc = b%nz
  drc = b%dr
  dzc = b%dz
  izlbnd = b%izlbnd
  izrbnd = b%izrbnd
  zmin_in = b%zmin!-b%dz
  zmax_in = b%zmax!+b%dz

  conductors%interior%n = 0
  conductors%evensubgrid%n = 0
  conductors%oddsubgrid%n = 0

  call srfrvout_rz(rofzfunc,volt,zmin,zmax,xcent,rmax,lfill,  &
                   grids_ptr(igrid)%grid%rmin,grids_ptr(igrid)%grid%rmax,lshell,                      &
                   zmin_in,zmax_in,zbeam,drc,dzc,nrc,nzc,             &
                   -NINT(grids_ptr(igrid)%grid%rmin/b%dr), &
                   xmesh,l2symtry_in,l4symtry_in,condid)

  call addconductors_rz(b,nrc,nzc,drc,dzc,grids_ptr(igrid)%grid%rmin, &
                        ixlbnd,ixrbnd,izlbnd,izrbnd,conductors)
  
 end do
end do
conductors%interior%n = 0
conductors%evensubgrid%n = 0
conductors%oddsubgrid%n = 0

return
end subroutine srfrvoutrz

subroutine srfrvinoutrz(rminofz,rmaxofz,volt,zmin,zmax,xcent,   &
                           lzend,xmin,xmax,lshell,          &
                           zmmin,zmmax,zbeam,dx,dz,nx,nz,       &
                           ix_axis,xmesh,l2symtry_in,l4symtry_in,condid)
! call subroutine srfrvinout_rz (which determines grid nodes near conductors and give
! directions and distances to conductors), initialize and assign coefficients
! for multigrid Poisson solver.
use Conductor3d
USE multigridrz
implicit none
character(*) rminofz,rmaxofz
real(kind=8):: volt,zmin,zmax,xcent
LOGICAL(ISZ):: lzend,lshell,l2symtry_in,l4symtry_in
real(kind=8):: xmin,xmax
real(kind=8):: zmmin,zmmax,zbeam,dx,dz
INTEGER(ISZ):: nx,ny,nz,ix_axis,jmin,jmax
real(kind=8):: xmesh(0:nx)
integer(ISZ):: condid
LOGICAL(ISZ) :: doloop

INTEGER(ISZ) :: i,nrc,nzc,igrid
REAL(8) :: drc,dzc,zmin_in,zmax_in

TYPE(CONDtype), POINTER :: c
TYPE(BNDtype), POINTER :: b

IF(solvergeom==Zgeom .or. solvergeom==Rgeom) return

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmin,zmmin)

do igrid=1,ngrids
 IF(zmin>grids_ptr(igrid)%grid%zmax .OR. zmax<grids_ptr(igrid)%grid%zmin .OR. &
    xmin>grids_ptr(igrid)%grid%rmax .OR. xmax<grids_ptr(igrid)%grid%rmin) cycle
 nlevels=grids_ptr(igrid)%grid%nlevels
 ixlbnd = grids_ptr(igrid)%grid%ixlbnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 do i = 1, nlevels
  IF(i == 1) then
    b => grids_ptr(igrid)%grid%bndfirst
  else
    b => b%next
  END if
  nrc = b%nr
  nzc = b%nz
  drc = b%dr
  dzc = b%dz
  izlbnd = b%izlbnd
  izrbnd = b%izrbnd
  zmin_in = b%zmin!-b%dz
  zmax_in = b%zmax!+b%dz

  conductors%interior%n = 0
  conductors%evensubgrid%n = 0
  conductors%oddsubgrid%n = 0

  call srfrvinout_rz(rminofz,rmaxofz,volt,zmin,zmax,xcent,  &
                     lzend,grids_ptr(igrid)%grid%rmin,grids_ptr(igrid)%grid%rmax,lshell,                &
                     zmin_in,zmax_in,zbeam,drc,dzc,nrc,nzc,             &
                     -NINT(grids_ptr(igrid)%grid%rmin/b%dr), &
                     xmesh,l2symtry_in,l4symtry_in,condid)

  call addconductors_rz(b,nrc,nzc,drc,dzc,grids_ptr(igrid)%grid%rmin, &
                        ixlbnd,ixrbnd,izlbnd,izrbnd,conductors)
  
 end do
end do
conductors%interior%n = 0
conductors%evensubgrid%n = 0
conductors%oddsubgrid%n = 0
return
end subroutine srfrvinoutrz

subroutine setcndtrrz(xmmin,zmmin,zbeam,zgrid,nx,nz,dx,dz, &
                      bound0_in,boundnz_in,boundxy_in,l2symtry_in,l4symtry_in)
USE multigridrz
use Conductor3d
integer(ISZ):: nx,nz
real(kind=8):: xmmin,zmmin,zbeam,zgrid,dx,dz
integer(ISZ):: bound0_in,boundnz_in,boundxy_in
logical(ISZ):: l2symtry_in,l4symtry_in

INTEGER(ISZ) :: i,nrc,nzc,igrid,jmin,jmax
REAL(8) :: drc,dzc,rmin_in,zmin_in
LOGICAL(ISZ) :: doloop
TYPE(BNDtype), POINTER :: b

IF(solvergeom==Zgeom .or. solvergeom==Rgeom) return

IF(.not.associated(basegrid)) call init_basegrid(nx,nz,dx,dz,xmmin,zmmin)

do igrid=1,ngrids
 nlevels=grids_ptr(igrid)%grid%nlevels
 ixlbnd = grids_ptr(igrid)%grid%ixlbnd
 ixrbnd = grids_ptr(igrid)%grid%ixrbnd
 rmin_in = grids_ptr(igrid)%grid%rmin
 do i = 1, nlevels
  IF(i == 1) then
    b => grids_ptr(igrid)%grid%bndfirst
  else
    b => b%next
  END if
  nrc = b%nr
  nzc = b%nz
  drc = b%dr
  dzc = b%dz
  izlbnd = b%izlbnd
  izrbnd = b%izrbnd
  zmin_in = b%zmin!-b%dz

  conductors%interior%n = 0
  conductors%evensubgrid%n = 0
  conductors%oddsubgrid%n = 0

  call setcndtr_rz(rmin_in,zmin_in,zbeam,zgrid,nrc,nzc,drc,dzc, &
                   bound0_in,boundnz_in,boundxy_in,l2symtry_in,l4symtry_in)

  call addconductors_rz(b,nrc,nzc,drc,dzc,grids_ptr(igrid)%grid%rmin, &
                        ixlbnd,ixrbnd,izlbnd,izrbnd,conductors)
  
 end do
end do
conductors%interior%n = 0
conductors%evensubgrid%n = 0
conductors%oddsubgrid%n = 0

END subroutine setcndtrrz

subroutine install_conductors_rz(conductors,grid)
USE Multigrid3d
USE multigridrz
use ConductorTypemodule
implicit none
TYPE(GRIDtype) :: grid
type(ConductorType):: conductors

INTEGER(ISZ), DIMENSION(:), allocatable :: mg_ncond,mg_necndbdy, mg_nocndbdy
INTEGER(ISZ) :: i, ii, nrc, nzc, itmp
INTEGER(ISZ) :: ncondtmp, necndbdytmp, nocndbdytmp
REAL(8) :: drc, dzc

type(ConductorType),pointer:: conductorstmp

TYPE(BNDtype), POINTER :: b

IF(solvergeom==Zgeom .or. solvergeom==Rgeom) return

ALLOCATE(mg_ncond(grid%nlevels),mg_necndbdy(grid%nlevels), mg_nocndbdy(grid%nlevels))
mg_ncond = 0
mg_necndbdy = 0
mg_nocndbdy = 0

ixlbnd = grid%ixlbnd
ixrbnd = grid%ixrbnd
do i = 1, conductors%interior%n
  ii = conductors%interior%ilevel(i) + 1
  mg_ncond(ii) = mg_ncond(ii) + 1
end do
do i = 1, conductors%evensubgrid%n
  ii = conductors%evensubgrid%ilevel(i) + 1
  mg_necndbdy(ii) = mg_necndbdy(ii) + 1
end do
do i = 1, conductors%oddsubgrid%n
  ii = conductors%oddsubgrid%ilevel(i) + 1
  mg_nocndbdy(ii) = mg_nocndbdy(ii) + 1
end do

! --- conductorstmp must be created this way since the gallot routine is called
! --- with it. The gallot needs to have the object accessible from python,
! --- which happens when the object is created this way. The alternative
! --- is to explicitly call the InitPyRef and DecRef routines for
! --- conductorstmp and all of the objects it refers to.
conductorstmp => NewConductorType()

 do i = 1, grid%nlevels
  IF(i == 1) then
    b => grid%bndfirst
  else
    b => b%next
  END if
  level = i
  nrc = b%nr
  nzc = b%nz
  drc = b%dr
  dzc = b%dz
  izlbnd = b%izlbnd
  izrbnd = b%izrbnd

  conductorstmp%interior%nmax    = mg_ncond(i)
  conductorstmp%evensubgrid%nmax = mg_necndbdy(i)
  conductorstmp%oddsubgrid%nmax  = mg_nocndbdy(i)
  call ConductorTypeallot(conductorstmp)

  conductorstmp%interior%n    = mg_ncond(i)
  conductorstmp%evensubgrid%n = mg_necndbdy(i)
  conductorstmp%oddsubgrid%n  = mg_nocndbdy(i)

  itmp = 0
  do ii = 1, conductors%interior%n
    IF(conductors%interior%ilevel(ii)+1==i) then
      itmp = itmp + 1
      conductorstmp%interior%indx(:,itmp) = conductors%interior%indx(:,ii)
      conductorstmp%interior%volt(itmp)   = conductors%interior%volt(ii)
      conductorstmp%interior%numb(itmp)   = conductors%interior%numb(ii)
    END if
  end do
  itmp = 0
  do ii = 1, conductors%evensubgrid%n
    IF(conductors%evensubgrid%ilevel(ii)+1==i) then
      itmp = itmp + 1
      conductorstmp%evensubgrid%indx(:,itmp) = conductors%evensubgrid%indx(:,ii)
      conductorstmp%evensubgrid%dels(:,itmp) = conductors%evensubgrid%dels(:,ii)
      conductorstmp%evensubgrid%volt(:,itmp) = conductors%evensubgrid%volt(:,ii)
      conductorstmp%evensubgrid%numb(:,itmp) = conductors%evensubgrid%numb(:,ii)
    END if
  end do
  itmp = 0
  do ii = 1, conductors%oddsubgrid%n
    IF(conductors%oddsubgrid%ilevel(ii)+1==i) then
      itmp = itmp + 1
      conductorstmp%oddsubgrid%indx(:,itmp) = conductors%oddsubgrid%indx(:,ii)
      conductorstmp%oddsubgrid%dels(:,itmp) = conductors%oddsubgrid%dels(:,ii)
      conductorstmp%oddsubgrid%volt(:,itmp) = conductors%oddsubgrid%volt(:,ii)
      conductorstmp%oddsubgrid%numb(:,itmp) = conductors%oddsubgrid%numb(:,ii)
    END if
  end do
  call addconductors_rz(b,nrc,nzc,drc,dzc,grid%rmin,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                        conductorstmp)
 end do

call ReleaseConductorType(conductorstmp)
DEALLOCATE(mg_ncond,mg_necndbdy, mg_nocndbdy)

conductors%interior%n = 0
conductors%evensubgrid%n = 0
conductors%oddsubgrid%n = 0

IF(grid%gid(1)==basegrid%gid(1)) call get_cond_rz(1)

return
end subroutine install_conductors_rz


subroutine addconductors_rz(b,nrc,nzc,drc,dzc,rmin,ixlbnd,ixrbnd,izlbnd,izrbnd, &
                            conductors)
USE InGen3d, ONLY:solvergeom,RZgeom,XYZgeom,XZgeom,XYgeom
USE multigridrz, ONLY: CONDtype, dirichlet, patchbnd, v_cond, v_bnd, v_dirichlet, bnd_method, egun, ecb, init_bnd_sublevel
USE BNDtypemodule
use ConductorTypemodule
implicit none

TYPE(BNDtype) :: b
INTEGER(ISZ), INTENT(IN) :: nrc,nzc,ixlbnd,ixrbnd,izlbnd,izrbnd
REAL(8), INTENT(IN) :: drc,dzc,rmin
type(ConductorType),intent(in):: conductors

INTEGER(ISZ) :: ii,iii,iv,iiv,nxbndmin,nxbndmax,nzbndmin,nzbndmax,iivmin,iivmax,ibnd,ne,no,kl
REAL(8) :: dt,dxm,dxp,dzm,dzp,r,rp,rm,dxx,dzz

TYPE(CONDtype), POINTER :: c

  IF(solvergeom==XYgeom) then
    kl = 1
  else
    kl = 2
  END if

  nxbndmin=0
  nxbndmax=nrc
  nzbndmin=0
  nzbndmax=nzc
  IF(ixlbnd==dirichlet .or. ixlbnd==patchbnd) nxbndmin=nxbndmin+1
  IF(ixrbnd==dirichlet .or. ixrbnd==patchbnd) nxbndmax=nxbndmax-1
  IF(izlbnd==dirichlet .or. izlbnd==patchbnd) nzbndmin=nzbndmin+1
  IF(izrbnd==dirichlet .or. izrbnd==patchbnd) nzbndmax=nzbndmax-1

  call init_bnd_sublevel(b,conductors%evensubgrid%n+conductors%oddsubgrid%n, &
                         conductors%interior%n)

  b%cndlast%nbbndred = conductors%evensubgrid%n
  b%cndlast%nbbnd = conductors%evensubgrid%n + conductors%oddsubgrid%n
  b%cndlast%ncond = conductors%interior%n
  iii=0
  do ii=1,conductors%interior%n
    iii = iii + 1
    b%cndlast%jcond(iii) = conductors%interior%indx(0 ,ii)+1
    b%cndlast%kcond(iii) = conductors%interior%indx(kl,ii)+1

!    IF(b%v(b%cndlast%jcond(iii),b%cndlast%kcond(iii))==v_cond) then
!      iii = iii-1
!      b%cndlast%ncond = b%cndlast%ncond-1
!      cycle
!    END if
    b%v(b%cndlast%jcond(iii),b%cndlast%kcond(iii)) = v_cond
    b%cndlast%voltage(iii) = conductors%interior%volt(ii)
    b%cndlast%condid(iii) = conductors%interior%numb(ii)
  end do

  ii = 0
  ne = 0
  no = 0
  b%cndlast%nbbndred = 0
  b%cndlast%nbbnd    = 0
  do ibnd = 1, conductors%evensubgrid%n+conductors%oddsubgrid%n
   ii = ii + 1
   IF(ibnd<=conductors%evensubgrid%n) then
     iii = ibnd
     b%cndlast%jj(ii)  = conductors%evensubgrid%indx(0 ,iii)+1
     b%cndlast%kk(ii)  = conductors%evensubgrid%indx(kl,iii)+1
   else
     iii = ibnd - conductors%evensubgrid%n
     b%cndlast%jj(ii)  = conductors%oddsubgrid%indx(0 ,iii)+1
     b%cndlast%kk(ii)  = conductors%oddsubgrid%indx(kl,iii)+1
   END if
   IF( (b%v(b%cndlast%jj(ii),b%cndlast%kk(ii))==v_cond) .or.                 &
       (b%v(b%cndlast%jj(ii),b%cndlast%kk(ii))==v_dirichlet) .or.              &
     (.not. (b%cndlast%jj(ii)>=nxbndmin+1 .and. b%cndlast%jj(ii)<=nxbndmax+1 .and. &
             b%cndlast%kk(ii)>=nzbndmin+1 .and. b%cndlast%kk(ii)<=nzbndmax+1))) then
      ii = ii - 1
      cycle
   END if
   IF(ibnd<=conductors%evensubgrid%n) then
     ne = ne + 1
     dxm = conductors%evensubgrid%dels(0     ,iii)*b%dr
     dxp = conductors%evensubgrid%dels(1     ,iii)*b%dr
     dzm = conductors%evensubgrid%dels(2*kl  ,iii)*b%dz
     dzp = conductors%evensubgrid%dels(2*kl+1,iii)*b%dz
     b%cndlast%volt0xm(ii)=conductors%evensubgrid%volt(0     ,iii)
     b%cndlast%volt0xp(ii)=conductors%evensubgrid%volt(1     ,iii)
     b%cndlast%volt0zm(ii)=conductors%evensubgrid%volt(2*kl  ,iii)
     b%cndlast%volt0zp(ii)=conductors%evensubgrid%volt(2*kl+1,iii)
     b%cndlast%condidxm(ii)=conductors%evensubgrid%numb(0     ,iii)
     b%cndlast%condidxp(ii)=conductors%evensubgrid%numb(1     ,iii)
     b%cndlast%condidzm(ii)=conductors%evensubgrid%numb(2*kl  ,iii)
     b%cndlast%condidzp(ii)=conductors%evensubgrid%numb(2*kl+1,iii)
   else
     no = no + 1
     dxm = conductors%oddsubgrid%dels(0     ,iii)*b%dr
     dxp = conductors%oddsubgrid%dels(1     ,iii)*b%dr
     dzm = conductors%oddsubgrid%dels(2*kl  ,iii)*b%dz
     dzp = conductors%oddsubgrid%dels(2*kl+1,iii)*b%dz
     b%cndlast%volt0xm(ii)=conductors%oddsubgrid%volt(0     ,iii)
     b%cndlast%volt0xp(ii)=conductors%oddsubgrid%volt(1     ,iii)
     b%cndlast%volt0zm(ii)=conductors%oddsubgrid%volt(2*kl  ,iii)
     b%cndlast%volt0zp(ii)=conductors%oddsubgrid%volt(2*kl+1,iii)
     b%cndlast%condidxm(ii)=conductors%oddsubgrid%numb(0     ,iii)
     b%cndlast%condidxp(ii)=conductors%oddsubgrid%numb(1     ,iii)
     b%cndlast%condidzm(ii)=conductors%oddsubgrid%numb(2*kl  ,iii)
     b%cndlast%condidzp(ii)=conductors%oddsubgrid%numb(2*kl+1,iii)
   END if
   b%cndlast%docalc(ii)=.true.
   IF(b%v(b%cndlast%jj(ii),b%cndlast%kk(ii))/=v_bnd ) then
      b%v(b%cndlast%jj(ii),b%cndlast%kk(ii)) = v_bnd
   else
     do iv=1, b%nb_conductors
       IF(iv==1) then
         c => b%cndfirst
       else
         c => c%next
       END if
       IF(ibnd<=conductors%evensubgrid%n) then
         iivmin = 1
         iivmax = c%nbbndred
       else
         iivmin = c%nbbndred+1
         iivmax = c%nbbnd
       END if
       do iiv=iivmin,iivmax
         IF(b%cndlast%jj(ii)==c%jj(iiv) .AND. b%cndlast%kk(ii)==c%kk(iiv)) then
           c%docalc(iiv)=.false.
           IF(c%dxm(iiv)<dxm) then
             dxm = c%dxm(iiv)
             b%cndlast%volt0xm(ii) = c%volt0xm(iiv)
             b%cndlast%condidxm(ii) = c%condidxm(iiv)
           END if
           IF(c%dxp(iiv)<dxp) then
             dxp = c%dxp(iiv)
             b%cndlast%volt0xp(ii) = c%volt0xp(iiv)
             b%cndlast%condidxp(ii) = c%condidxp(iiv)
           END if
           IF(c%dzm(iiv)<dzm) then
             dzm = c%dzm(iiv)
             b%cndlast%volt0zm(ii) = c%volt0zm(iiv)
             b%cndlast%condidzm(ii) = c%condidzm(iiv)
           END if
           IF(c%dzp(iiv)<dzp) then
             dzp = c%dzp(iiv)
             b%cndlast%volt0zp(ii) = c%volt0zp(iiv)
             b%cndlast%condidzp(ii) = c%condidzp(iiv)
           END if
         END if
       end do
     end do
   endif
   b%cndlast%dxm(ii)=dxm
   b%cndlast%dxp(ii)=dxp
   b%cndlast%dzm(ii)=dzm
   b%cndlast%dzp(ii)=dzp
   dxm = MIN(b%dr,dxm)
   dxp = MIN(b%dr,dxp)
   dzm = MIN(b%dz,dzm)
   dzp = MIN(b%dz,dzp)
   select case (bnd_method)
     case (egun)
       dxx=b%dr
       dzz=b%dz
     case (ecb)
       dxx=0.5_8*(dxp+dxm)  !ecb
       dzz=0.5_8*(dzp+dzm)  !ecb
     case default
   end select
   IF(solvergeom==RZgeom) then
    IF(b%cndlast%jj(ii)==1 .and. rmin==0.) then
     b%cndlast%dt(ii) = 1._8/(4._8/(dxp*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     b%cndlast%cfxm(ii) = 0.
     b%cndlast%cfxp(ii) = 4._8/(dxp*dxx)
    else
     r = rmin+(b%cndlast%jj(ii)-1)*b%dr
     select case (bnd_method)
       case (egun)
         rm = r-0.5_8*b%dr
         rp = r+0.5_8*b%dr
       case (ecb)
         rm = r-0.5_8*dxm
         rp = r+0.5_8*dxp
       case default
     end select
     b%cndlast%dt(ii) = 1._8/((rm/dxm+rp/dxp)/(r*dxx)+(1._8/dzm+1._8/dzp)/dzz)
     b%cndlast%cfxm(ii) = rm/(r*dxm*dxx)
     b%cndlast%cfxp(ii) = rp/(r*dxp*dxx)
    END if
   else ! (solvergeom==XZgeom)
     b%cndlast%dt(ii) = 1._8/((1._8/dxm+1._8/dxp)/dxx+(1._8/dzm+1._8/dzp)/dzz)
     b%cndlast%cfxm(ii) = 1._8/(dxm*dxx)
     b%cndlast%cfxp(ii) = 1._8/(dxp*dxx)
   END if
   b%cndlast%cfzm(ii) = 1._8/(dzm*dzz)
   b%cndlast%cfzp(ii) = 1._8/(dzp*dzz)
   b%cndlast%cf0(ii)  = -b%cndlast%cfxm(ii)-b%cndlast%cfxp(ii)-b%cndlast%cfzm(ii)-b%cndlast%cfzp(ii)
   IF(b%cndlast%dxm(ii)>=b%dr) then
     b%cndlast%phi0xm(ii)=0._8
   else
     b%cndlast%phi0xm(ii)=b%cndlast%cfxm(ii)*b%cndlast%volt0xm(ii)
     b%cndlast%cfxm(ii)=0._8
   END if
   IF(b%cndlast%dxp(ii)>=b%dr) then
     b%cndlast%phi0xp(ii)=0._8
   else
     b%cndlast%phi0xp(ii)=b%cndlast%cfxp(ii)*b%cndlast%volt0xp(ii)
     b%cndlast%cfxp(ii)=0._8
   END if
   IF(b%cndlast%dzm(ii)>=b%dz) then
     b%cndlast%phi0zm(ii)=0._8
   else
     b%cndlast%phi0zm(ii)=b%cndlast%cfzm(ii)*b%cndlast%volt0zm(ii)
     b%cndlast%cfzm(ii)=0._8
   END if
   IF(b%cndlast%dzp(ii)>=b%dz) then
     b%cndlast%phi0zp(ii)=0._8
   else
     b%cndlast%phi0zp(ii)=b%cndlast%cfzp(ii)*b%cndlast%volt0zp(ii)
     b%cndlast%cfzp(ii)=0._8
   END if
   IF(ibnd<=conductors%evensubgrid%n) b%cndlast%nbbndred = b%cndlast%nbbndred+1
   b%cndlast%nbbnd = b%cndlast%nbbnd + 1
  end do

  ! --- nbbnd has been changed above and since some points may be rejected,
  ! --- it may be less then its value when the allocate was done. To make
  ! --- the array lengths consistent, a change is done on the arrays.
  call CONDtypechange(b%cndlast)

end subroutine addconductors_rz

subroutine test_subgrid_rz()
USE multigridrz
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ) :: ic, i, j, l, ii, jsw, lsw, redblack,iil, iiu, nri, nrf, nzi, nzf, lswinit
TYPE(GRIDtype), POINTER :: g
TYPE(BNDtype), POINTER :: b
TYPE(CONDtype), POINTER :: c
LOGICAL :: cond
REAL(8) :: f


do i = 1, ngrids
IF(i==1) then
  g => basegrid
else
  IF(associated(g%next)) then
    g => g%next
  else
    g => g%down
  END if
END if
b => g%bndfirst
g%phi = 0.

lswinit = 2
IF(g%ixlbnd==dirichlet .or. g%ixlbnd==patchbnd) then
  nri=2
else
  nri=1
  lswinit = 3-lswinit
END if
IF(g%ixrbnd==dirichlet .or. g%ixrbnd==patchbnd) then
  nrf=g%nr-1
else
  nrf=g%nr
END if
IF(g%izlbnd==dirichlet .or. g%izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
  lswinit = 3-lswinit
END if
IF(g%izrbnd==dirichlet .or. g%izrbnd==patchbnd) then
  nzf=g%nz-1
else
  nzf=g%nz
END if

lsw = lswinit
do redblack = 1, 2
  jsw = lsw
  IF(lsw==1) then
    f = 1.
  else
    f = -1.
  END if
  do ic = 1, b%nb_conductors
    IF(ic==1) then
      c => b%cndfirst
    else
      c => c%next
    END if

    IF(redblack==1) THEN !red
      iil=1
      iiu=c%nbbndred
    else !black
      iil=c%nbbndred+1
      iiu=c%nbbnd
    ENDif
    do ii = iil, iiu
      j = c%jj(ii)
      l = c%kk(ii)
      cond = c%docalc(ii).and.b%v(j,l)==v_bnd
      if (cond) g%phi(j,l) = g%phi(j,l) + f
    ENDDO
  END do
  IF(vlocs) then
    IF(redblack==1) THEN !red
      iil=1
      iiu=b%nvlocsred
    else !black
      iil=b%nvlocsred+1
      iiu=b%nvlocs
    ENDif
    do ii = iil, iiu
      j = b%vlocs_j(ii)
      l = b%vlocs_k(ii)
      g%phi(j,l) = g%phi(j,l) + f
    end do
  else
    do l = nzi, nzf+1
      IF(nri==1 .and. jsw==2) then! origin
        j = 1
        g%phi(j,l) = g%phi(j,l) + f
      END if
      do j = nri+jsw, nrf+1, 2
        IF(b%v(j,l)==v_vacuum) &
        g%phi(j,l) = g%phi(j,l) + f
      end do
      jsw = 3-jsw
    end do
    lsw = 3-lsw
  END if

END do !redblack=1, 2
end do

END subroutine test_subgrid_rz

subroutine test_subgrid_xz()
USE multigridrz
! make a relaxation step. Grid is assumed to have guard cells.
implicit none

INTEGER(ISZ) :: ic, j, l, ii, jsw, lsw, redblack,iil, iiu, nri, nrf, nzi, nzf
TYPE(GRIDtype), POINTER :: g
TYPE(BNDtype), POINTER :: b
TYPE(CONDtype), POINTER :: c
LOGICAL :: cond
REAL(8) :: f

g => basegrid
b => g%bndfirst

  g%phi = 0.

lsw = 1
IF(g%ixlbnd==dirichlet .or. g%ixlbnd==patchbnd) then
  nri=2
else
  nri=1
  lsw = 3-lsw
END if
IF(g%ixrbnd==dirichlet .or. g%ixrbnd==patchbnd) then
  nrf=g%nr-1
else
  nrf=g%nr
END if
IF(g%izlbnd==dirichlet .or. izlbnd==patchbnd) then
  nzi=2
else
  nzi=1
  lsw = 3-lsw
END if
IF(g%izrbnd==dirichlet .or. izrbnd==patchbnd) then
  nzf=g%nz-1
else
  nzf=g%nz
END if

do redblack = 1, 2
  jsw = lsw
  IF(lsw==1) then
    f = 1.
  else
    f = -1.
  END if
  do ic = 1, b%nb_conductors
    IF(ic==1) then
      c => b%cndfirst
    else
      c => c%next
    END if

    IF(redblack==1) THEN !red
      iil=1
      iiu=c%nbbndred
    else !black
      iil=c%nbbndred+1
      iiu=c%nbbnd
    ENDif
    do ii = iil, iiu
      j = c%jj(ii)
      l = c%kk(ii)
      cond = c%docalc(ii).and.b%v(j,l)==v_bnd
      if (cond) g%phi(j,l) = g%phi(j,l) + f
    ENDDO
  END do
  IF(vlocs) then
    IF(redblack==1) THEN !red
      iil=1
      iiu=b%nvlocsred
    else !black
      iil=b%nvlocsred+1
      iiu=b%nvlocs
    ENDif
    do ii = iil, iiu
      j = b%vlocs_j(ii)
      l = b%vlocs_k(ii)
      g%phi(j,l) = g%phi(j,l) + f
    end do
  else
    do l = nzi, nzf+1
      do j = nri+jsw-1, nrf+1, 2
        IF(b%v(j,l)==v_vacuum) &
        g%phi(j,l) = g%phi(j,l) + f
      end do
      jsw = 3-jsw
    end do
    lsw = 3-lsw
  END if

END do !redblack=1, 2

END subroutine test_subgrid_xz

!=============================================================================
subroutine gtlchgrz
USE multigridrz
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
if (.not. lgtlchg3d) return

dzi = 1./dz

!conversion factor to go from grid frame to beam frame
zz = zgrid + zmmin - zzmin - zbeam

do iz=0,nzzarr

  izg = (iz*dzz - zz)*dzi
  wzg = (iz*dzz - zz)*dzi - izg
  izg = izg+1

  if (1 <= izg .and. izg <= nz+1) then
   linechg(iz) = 0.
   do j = 1, basegrid%nr+1
    linechg(iz) = linechg(iz) + basegrid%rho(j,izg)/(basegrid%invvol(j)*basegrid%dz)*(1. - wzg)
   end do
  else
    linechg(iz) = 0.
  endif
  if (1 <= izg+1 .and. izg+1 <= nz+1) then
   do j = 1, basegrid%nr+1
    linechg(iz) = linechg(iz) + basegrid%rho(j,izg+1)/(basegrid%invvol(j)*basegrid%dz)*wzg
   end do
  endif

enddo

return
end subroutine gtlchgrz

subroutine dep_rho_rz(is,rho,nr,nz,dr,dz,xmin,zmin)
USE Constant
use Particles
implicit none

INTEGER(ISZ), INTENT(IN) :: is, nr, nz
REAL(8), INTENT(IN OUT) :: rho(nr+1,nz+1)
REAL(8), INTENT(IN) :: dr, dz, zmin, xmin

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr), invvolxz
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp
REAL(8):: q, qw

  invdr = 1._8/dr
  invdz = 1._8/dz

  ! computes divider by cell volumes to get density
  IF(xmin==0.) then
    j = 0
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
    do j = 1, nr
      invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr * dz)
    end do
  else
    do j = 0, nr
      invvol(j) = 1._8 / (2._8 * pi * (xmin+real(j,8)*dr) * dr * dz)
    end do
  END if

  q = sq(is)*sw(is)

  ! make charge deposition using CIC weighting
  IF(wpid==0) then
   do i = ins(is), ins(is) + nps(is) - 1
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-xmin)*invdr
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
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-xmin)*invdr
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


subroutine rhoweightrz(xp,yp,zp,np,q,nr,nz,dr,dz,rgrid,zgrid)
USE multigridrz
USE FRZmgrid, ONLY: mgridrz_xfact, mgridrz_yfact
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), INTENT(IN) :: q, dr, dz, rgrid, zgrid

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr), invvolxz, zmin0, qw
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp
REAL(8):: substarttime

IF(np==0) return

if(mgridrz_deform) then
  call rhoweightrz_deform(xp(1),yp(1),zp(1),np, &
                          q,nr,nz,dr,dz, &
                          mgridrz_xfact, mgridrz_yfact, rgrid, zgrid)
  return
END if

IF(solvergeom==RZgeom) then
 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightrz_amr(xp,yp,zp,np,q,zgrid)
 else
  invdr = 1._8/dr
  invdz = 1._8/dz

  ! computes divider by cell volumes to get density
  IF(rgrid==0.) then
    j = 0
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
    do j = 1, nr
      invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr * dz)
    end do
  else
    do j = 0, nr
      invvol(j) = 1._8 / (2._8 * pi * (rgrid+real(j,8)*dr) * dr * dz)
    end do
  END if

  ! make charge deposition using CIC weighting
  do i = 1, np
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-rgrid)*invdr
    zpos = (zp(i)-basegrid%zminp-zgrid)*invdz
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
 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightxz_amr(xp,zp,np,q,zgrid)
 else
  invdr = 1._8/dr
  invdz = 1._8/dz

  invvolxz = 1._8 / (dr*dz)

  IF( solvergeom==XYgeom) then
    zmin0 = basegrid%zminp
  else ! solvergeom=XZgeom
    zmin0 = basegrid%zminp-zgrid
  END if

  qw = q
  if (solvergeom==XZgeom .and. l4symtry) qw=qw*0.5
  if (solvergeom==XYgeom) then
    if(l4symtry) then
      qw=qw*0.25
    elseif(l2symtry) then
      qw=qw*0.5
    end if
  end if

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(l4symtry) then
      rpos = (ABS(xp(i))-basegrid%xmin)*invdr
    else
      rpos = (xp(i)-basegrid%xmin)*invdr
    END if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      zpos = (ABS(zp(i))-zmin0)*invdz
    else
      zpos = (zp(i)-zmin0)*invdz
    END if
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
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

return
END SUBROUTINE RHOWEIGHTRZ

subroutine rhoweightrz_weights(xp,yp,zp,w,np,q,nr,nz,dr,dz,rgrid,zgrid)
USE multigridrz
USE Subtimers3d
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, w
REAL(8), INTENT(IN) :: q, dr, dz, rgrid, zgrid

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr), invvolxz, qw, q0, zmin0
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp

IF(np==0) return

IF(solvergeom==RZgeom) then
 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightrz_amr_weights(xp,yp,zp,w,np,q,zgrid)
 else
  invdr = 1._8/dr
  invdz = 1._8/dz

  ! computes divider by cell volumes to get density
  IF(rgrid==0.) then
    j = 0
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
    do j = 1, nr
      invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr * dz)
    end do
  else
    do j = 0, nr
      invvol(j) = 1._8 / (2._8 * pi * (rgrid+real(j,8)*dr) * dr * dz)
    end do
  END if

  ! make charge deposition using CIC weighting
  do i = 1, np
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-rgrid)*invdr
    zpos = (zp(i)-basegrid%zminp-zgrid)*invdz
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
 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightxz_amr_weights(xp,zp,w,np,q,zgrid)
 else
  invdr = 1._8/dr
  invdz = 1._8/dz

  invvolxz = 1._8 / (dr*dz)

  q0 = q
  if (solvergeom==XZgeom .and. l4symtry) q0=q0*0.5
  if (solvergeom==XYgeom) then
    if(l4symtry) then
      q0=q0*0.25
    elseif(l2symtry) then
      q0=q0*0.5
    end if
  end if

  IF( solvergeom==XYgeom) then
    zmin0 = basegrid%zminp
  else ! solvergeom=XZgeom
    zmin0 = basegrid%zminp-zgrid
  END if

  ! make charge deposition using CIC weighting
  do i = 1, np
    IF(l4symtry) then
      rpos = abs(xp(i))*invdr
    else
      rpos = (xp(i)-basegrid%xmin)*invdr
    END if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      zpos = (ABS(zp(i))-zmin0)*invdz
    else
      zpos = (zp(i)-zmin0)*invdz
    END if
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    qw = q0*w(i)
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

return
END SUBROUTINE RHOWEIGHTRZ_weights

subroutine rhoweightr(xp,yp,np,q,nr,dr,rgrid)
USE multigridrz
USE FRZmgrid, ONLY: mgridrz_xfact, mgridrz_yfact
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp
REAL(8), INTENT(IN) :: q, dr, rgrid

REAL(8) :: invdr, rpos, ddr, oddr, invvol(0:nr)
INTEGER(ISZ) :: i, j, jn, jnp
REAL(8):: substarttime

IF(np==0) return

!if(mgridrz_deform) then
!  call rhoweightrz_deform(xp(1),yp(1),zp(1),np, &
!                          q,nr,nz,dr,dz, &
!                          mgridrz_xfact, mgridrz_yfact, rgrid, zgrid)
!  return
!END if

 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightr_amr(xp,yp,np,q)
 else
  invdr = 1._8/dr

  ! computes divider by cell volumes to get density
  IF(rgrid==0.) then
    j = 0
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr))
    do j = 1, nr
      invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr)
    end do
  else
    do j = 0, nr
      invvol(j) = 1._8 / (2._8 * pi * (rgrid+real(j,8)*dr) * dr)
    end do
  END if

  ! make charge deposition using CIC weighting
  do i = 1, np
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-rgrid)*invdr
    jn = 1+INT(rpos)
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
    jnp=jn+1
#ifdef MPIPARALLEL
    basegrid%rhop(jn, 1)  = basegrid%rhop(jn, 1)  + q * oddr * invvol(jn-1)
    basegrid%rhop(jnp,1)  = basegrid%rhop(jnp,1)  + q *  ddr * invvol(jnp-1)
    basegrid%rhominr = MIN(basegrid%rhominr,jn)
    basegrid%rhomaxr = MAX(basegrid%rhomaxr,jnp)
#else
    basegrid%rho(jn, 1)  = basegrid%rho(jn, 1)  + q * oddr * invvol(jn-1)
    basegrid%rho(jnp,1)  = basegrid%rho(jnp,1)  + q *  ddr * invvol(jnp-1)
#endif
  end do
 END if

return
END SUBROUTINE RHOWEIGHTR

subroutine rhoweightr_weights(xp,yp,w,np,q,nr,dr,rgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, w
REAL(8), INTENT(IN) :: q, dr, rgrid

REAL(8) :: invdr, rpos, ddr, oddr, invvol(0:nr), qw
INTEGER(ISZ) :: i, j, jn, jnp

IF(np==0) return

 IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  call rhoweightr_amr_weights(xp,yp,w,np,q)
 else
  invdr = 1._8/dr

  ! computes divider by cell volumes to get density
  IF(rgrid==0.) then
    j = 0
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr))
    do j = 1, nr
      invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr)
    end do
  else
    do j = 0, nr
      invvol(j) = 1._8 / (2._8 * pi * (rgrid+real(j,8)*dr) * dr)
    end do
  END if

  ! make charge deposition using CIC weighting
  do i = 1, np
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-rgrid)*invdr
    jn = 1+INT(rpos)
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
    jnp=jn+1
    qw = q*w(i)
#ifdef MPIPARALLEL
    basegrid%rhop(jn, 1)  = basegrid%rhop(jn, 1)  + qw * oddr * invvol(jn-1)
    basegrid%rhop(jnp,1)  = basegrid%rhop(jnp,1)  + qw *  ddr * invvol(jnp-1)
    basegrid%rhominr = MIN(basegrid%rhominr,jn)
    basegrid%rhomaxr = MAX(basegrid%rhomaxr,jnp)
#else
    basegrid%rho(jn, 1)  = basegrid%rho(jn, 1)  + qw * oddr  * invvol(jn-1)
    basegrid%rho(jnp,1)  = basegrid%rho(jnp,1)  + qw *  ddr  * invvol(jnp-1)
#endif
  end do
 END if

return
END SUBROUTINE RHOWEIGHTR_weights

subroutine rhoweightz(zp,np,q,nz,dz,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nz
REAL(8), DIMENSION(np), INTENT(IN) :: zp
REAL(8), INTENT(IN) :: q, dz

REAL(8) :: zpos, ddz, oddz, zgrid
INTEGER(ISZ) :: i, ln, lnp, igrid
REAL(8):: substarttime
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdz, zmin

IF(np==0) return

ALLOCATE(invdz(ngrids),zmin(ngrids))

do igrid = 1, ngrids
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  zmin (igrid) = grids_ptr(igrid)%grid%zminp+zgrid
end do

IF(ngrids>1 .and. .not. l_dep_rho_on_base) then
  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(1,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(1,ln)
        g=>grids_ptr(igrid)%grid
        zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        ln = 1+INT(zpos)
      END if
    end do
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    lnp=ln+1
#ifdef MPIPARALLEL
    g%rhop(1,ln)  = g%rhop(1,ln)  + q * oddz * invdz(igrid)
    g%rhop(1,lnp) = g%rhop(1,lnp) + q * ddz  * invdz(igrid)
    g%rhominz = MIN(g%rhominz,ln)
    g%rhomaxz = MAX(g%rhomaxz,lnp)
#else
    g%rho(1,ln)  = g%rho(1,ln)  + q * oddz * invdz(igrid)
    g%rho(1,lnp) = g%rho(1,lnp) + q * ddz  * invdz(igrid)
#endif
  end do
else
  ! make charge deposition using CIC weighting
  igrid = 1
  do i = 1, np
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

return
END SUBROUTINE RHOWEIGHTZ

subroutine set_rho_rz(rho,nr,nz,id)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
REAL(8), DIMENSION(nr+1,nz+1), INTENT(IN) :: rho

  IF(id<1.or.id>ngrids) then
    write(o_line,*) 'Error, id out of bounds'
    call remark(trim(o_line))
    write(o_line,*) 'given id = ',id,' while 1 < id < ',ngrids
    call remark(trim(o_line))
    stop
  END if
  IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
    write(o_line,*) 'Error, dimensions should be the same: '
    call remark(trim(o_line))
    write(o_line,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2)
    call remark(trim(o_line))
    write(o_line,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
    call remark(trim(o_line))
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
    write(o_line,*) 'Error, id out of bounds'
    call remark(trim(o_line))
    write(o_line,*) 'given id = ',id,' while 1 < id < ',ngrids
    call remark(trim(o_line))
    stop
  END if
  IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
    write(o_line,*) 'Error, dimensions should be the same: '
    call remark(trim(o_line))
    write(o_line,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2)
    call remark(trim(o_line))
    write(o_line,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
    call remark(trim(o_line))
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

  IF(id<1.or.id>ngrids) then
    write(o_line,*) 'Error, id out of bounds';                  call remark(trim(o_line))
    write(o_line,*) 'given id = ',id,' while 1 < id < ',ngrids; call remark(trim(o_line))
    return
  END if
#ifdef MPIPARALLEL
  IF(rhop==1) then
    IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rhop,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rhop,2)) then
      write(o_line,*) 'Error, dimensions should be the same: ';        call remark(trim(o_line))
      write(o_line,*) 'Nr, Nz for rhop    : ',SIZE(rho,1),SIZE(rho,2); call remark(trim(o_line))
      write(o_line,*) 'Nr, Nz for rhop(id): ',SIZE(grids_ptr(id)%grid%rhop,1),SIZE(grids_ptr(id)%grid%rhop,2)
      call remark(trim(o_line))
      return
    END if
    rho=grids_ptr(id)%grid%rhop
  else
#endif
    IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1).or.SIZE(rho,2)/=SIZE(grids_ptr(id)%grid%rho,2)) then
      write(o_line,*) 'Error, dimensions should be the same: ';       call remark(trim(o_line))
      write(o_line,*) 'Nr, Nz for rho    : ',SIZE(rho,1),SIZE(rho,2); call remark(trim(o_line))
      write(o_line,*) 'Nr, Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1),SIZE(grids_ptr(id)%grid%rho,2)
      call remark(trim(o_line))
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
    write(o_line,*) 'Error, id out of bounds';                  call remark(trim(o_line))
    write(o_line,*) 'given id = ',id,' while 1 < id < ',ngrids; call remark(trim(o_line))
    return
  END if
#ifdef MPIPARALLEL
  IF(rhop==1) then
    IF(SIZE(rho)/=SIZE(grids_ptr(id)%grid%rhop,2)) then
      write(o_line,*) 'Error, dimensions should be the same: ';            call remark(trim(o_line))
      write(o_line,*) 'Nz for rhop    : ',SIZE(rho,1);                     call remark(trim(o_line))
      write(o_line,*) 'Nz for rhop(id): ',SIZE(grids_ptr(id)%grid%rhop,2); call remark(trim(o_line))
      return
    END if
    rho=grids_ptr(id)%grid%rhop(1,:)
  else
#endif
    IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,2)) then
      write(o_line,*) 'Error, dimensions should be the same: ';          call remark(trim(o_line))
      write(o_line,*) 'Nz for rho    : ',SIZE(rho,1);                    call remark(trim(o_line))
      write(o_line,*) 'Nz for rho(id): ',SIZE(grids_ptr(id)%grid%rho,2); call remark(trim(o_line))
      return
    END if
    rho=grids_ptr(id)%grid%rho(1,:)
#ifdef MPIPARALLEL
  END if
#endif
return
end subroutine get_rho_z

subroutine get_rho_r(rho,nr,id,rhop)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,rhop
REAL(8), DIMENSION(nr+1), INTENT(IN OUT) :: rho

  IF(id<1.or.id>ngrids) then
    write(o_line,*) 'Error, id out of bounds';                  call remark(trim(o_line))
    write(o_line,*) 'given id = ',id,' while 1 < id < ',ngrids; call remark(trim(o_line))
    return
  END if
  IF(SIZE(rho,1)/=SIZE(grids_ptr(id)%grid%rho,1)) then
    write(o_line,*) 'Error, dimensions should be the same: ';          call remark(trim(o_line))
    write(o_line,*) 'Nr for rho    : ',SIZE(rho,1);                    call remark(trim(o_line))
    write(o_line,*) 'Nr for rho(id): ',SIZE(grids_ptr(id)%grid%rho,1); call remark(trim(o_line))
    return
  END if
  rho=grids_ptr(id)%grid%rho(:,1)

 return
end subroutine get_rho_r

subroutine rhoweightrz_amr(xp,yp,zp,np,q,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), INTENT(IN) :: q

REAL(8) :: r, rpos, zpos, ddr, ddz, oddr, oddz, zgrid
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdr, invdz, zmin
#ifdef MPIPARALLEL
#define RHO g%rhop
#else
#define RHO g%rho
#endif

ALLOCATE(invdr(ngrids),invdz(ngrids),zmin(ngrids))

do igrid = 1, ngrids
  invdr(igrid) = grids_ptr(igrid)%grid%invdr
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  zmin (igrid) = grids_ptr(igrid)%grid%zminp+zgrid
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
    r = SQRT(xp(i)*xp(i)+yp(i)*yp(i))    
    rpos = (r-g%rmin)*invdr(igrid)
    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,ln)
        g=>grids_ptr(igrid)%grid
        rpos = (r-g%rmin)*invdr(igrid)
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
    RHO(jn, ln)  = RHO(jn, ln)  + q * oddr * oddz * g%invvol(jn)
    RHO(jnp,ln)  = RHO(jnp,ln)  + q *  ddr * oddz * g%invvol(jnp)
    RHO(jn, lnp) = RHO(jn, lnp) + q * oddr *  ddz * g%invvol(jn)
    RHO(jnp,lnp) = RHO(jnp,lnp) + q *  ddr *  ddz * g%invvol(jnp)
  end do

  DEALLOCATE(invdr,invdz,zmin)

  return
END subroutine rhoweightrz_amr

subroutine rhoweightrz_amr_weights(xp,yp,zp,wp,np,q,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp, wp
REAL(8), INTENT(IN) :: q, zgrid

REAL(8) :: rpos, zpos, ddr, ddz, oddr, oddz, qw
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdr, invdz, zmin
#ifdef MPIPARALLEL
#define RHO g%rhop
#else
#define RHO g%rho
#endif

ALLOCATE(invdr(ngrids),invdz(ngrids),zmin(ngrids))

do igrid = 1, ngrids
  invdr(igrid) = grids_ptr(igrid)%grid%invdr
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  zmin (igrid) = grids_ptr(igrid)%grid%zminp+zgrid
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-g%rmin)*invdr(igrid)
    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,ln)
        g=>grids_ptr(igrid)%grid
        rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-g%rmin)*invdr(igrid)
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
    RHO(jn, ln)  = RHO(jn, ln)  + qw * oddr * oddz * g%invvol(jn)
    RHO(jnp,ln)  = RHO(jnp,ln)  + qw *  ddr * oddz * g%invvol(jnp)
    RHO(jn, lnp) = RHO(jn, lnp) + qw * oddr *  ddz * g%invvol(jn)
    RHO(jnp,lnp) = RHO(jnp,lnp) + qw *  ddr *  ddz * g%invvol(jnp)
  end do

  DEALLOCATE(invdr,invdz,zmin)

  return
END subroutine rhoweightrz_amr_weights

subroutine rhoweightxz_amr(xp,zp,np,q,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, zp
REAL(8), INTENT(IN) :: q

REAL(8) :: xpos, zpos, ddx, ddz, oddx, oddz, zgrid, zmin0, qw
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdx, invdz, zmin, invvol
#ifdef MPIPARALLEL
#define RHO g%rhop
#else
#define RHO g%rho
#endif
ALLOCATE(invdx(ngrids),invdz(ngrids),zmin(ngrids),invvol(ngrids))

IF( solvergeom==XYgeom) then
  zmin0 = 0.
else ! solvergeom=XZgeom
  zmin0 = zgrid
END if

qw = q
if (solvergeom==XZgeom .and. l4symtry) qw=qw*0.5
if (solvergeom==XYgeom) then
  if(l4symtry) then
    qw=qw*0.25
  elseif(l2symtry) then
    qw=qw*0.5
  end if
end if

do igrid = 1, ngrids
  invdx(igrid) = grids_ptr(igrid)%grid%invdr
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  zmin (igrid) = grids_ptr(igrid)%grid%zminp+zmin0
  invvol(igrid) = invdx(igrid)*invdz(igrid)
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
!    xpos = (xp(i)-g%rmin)*invdx(igrid)
!    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    IF(l4symtry) then
      xpos = (ABS(xp(i))-g%xmin)*invdx(igrid)
    else
      xpos = (xp(i)-g%xmin)*invdx(igrid)
    END if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      zpos = (ABS(zp(i))-zmin(igrid))*invdz(igrid)
    else
      zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    END if
    jn = 1+INT(xpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,ln)
        g=>grids_ptr(igrid)%grid
!        xpos = (xp(i)-g%rmin)*invdx(igrid)
!        zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        IF(l4symtry) then
          xpos = (ABS(xp(i))-g%xmin)*invdx(igrid)
        else
          xpos = (xp(i)-g%xmin)*invdx(igrid)
        END if
        IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
          zpos = (ABS(zp(i))-zmin(igrid))*invdz(igrid)
        else
          zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        END if
        jn = 1+INT(xpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddx = xpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddx = 1._8-ddx
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    RHO(jn, ln)  = RHO(jn, ln)  + qw * oddx * oddz * invvol(igrid)
    RHO(jnp,ln)  = RHO(jnp,ln)  + qw *  ddx * oddz * invvol(igrid)
    RHO(jn, lnp) = RHO(jn, lnp) + qw * oddx *  ddz * invvol(igrid)
    RHO(jnp,lnp) = RHO(jnp,lnp) + qw *  ddx *  ddz * invvol(igrid)
  end do

  DEALLOCATE(invdx,invdz,zmin)

  return
END subroutine rhoweightxz_amr

subroutine rhoweightxz_amr_weights(xp,zp,wp,np,q,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, zp, wp
REAL(8), INTENT(IN) :: q, zgrid

REAL(8) :: xpos, zpos, ddx, ddz, oddx, oddz, qw, q0, zmin0
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdx, invdz, zmin, invvol
#ifdef MPIPARALLEL
#define RHO g%rhop
#else
#define RHO g%rho
#endif

ALLOCATE(invdx(ngrids),invdz(ngrids),zmin(ngrids),invvol(ngrids))

IF( solvergeom==XYgeom) then
  zmin0 = 0.
else ! solvergeom=XZgeom
  zmin0 = zgrid
END if

q0 = q
if (solvergeom==XZgeom .and. l4symtry) q0=q0*0.5
if (solvergeom==XYgeom) then
  if(l4symtry) then
    q0=q0*0.25
  elseif(l2symtry) then
    q0=q0*0.5
  end if
end if

do igrid = 1, ngrids
  invdx(igrid) = grids_ptr(igrid)%grid%invdr
  invdz(igrid) = grids_ptr(igrid)%grid%invdz
  invvol(igrid) = invdx(igrid)*invdz(igrid)
  zmin (igrid) = grids_ptr(igrid)%grid%zminp+zmin0
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
!    xpos = (xp(i)-g%rmin)*invdx(igrid)
!    zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    IF(l4symtry) then
      xpos = (ABS(xp(i))-g%xmin)*invdx(igrid)
    else
      xpos = (xp(i)-g%xmin)*invdx(igrid)
    END if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      zpos = (ABS(zp(i))-zmin(igrid))*invdz(igrid)
    else
      zpos = (zp(i)-zmin(igrid))*invdz(igrid)
    END if
    jn = 1+INT(xpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,ln)
        g=>grids_ptr(igrid)%grid
!        xpos = (xp(i)-g%rmin)*invdx(igrid)
!        zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        IF(l4symtry) then
          xpos = (ABS(xp(i))-g%xmin)*invdx(igrid)
        else
          xpos = (xp(i)-g%xmin)*invdx(igrid)
        END if
        IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
          zpos = (ABS(zp(i))-zmin(igrid))*invdz(igrid)
        else
          zpos = (zp(i)-zmin(igrid))*invdz(igrid)
        END if
        jn = 1+INT(xpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddx = xpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddx = 1._8-ddx
    oddz = 1._8-ddz
    jnp=jn+1
    lnp=ln+1
    qw = q0*wp(i)
    RHO(jn, ln)  = RHO(jn, ln)  + qw * oddx * oddz * invvol(igrid)
    RHO(jnp,ln)  = RHO(jnp,ln)  + qw *  ddx * oddz * invvol(igrid)
    RHO(jn, lnp) = RHO(jn, lnp) + qw * oddx *  ddz * invvol(igrid)
    RHO(jnp,lnp) = RHO(jnp,lnp) + qw *  ddx *  ddz * invvol(igrid)
  end do

  DEALLOCATE(invdx,invdz,zmin)

  return
END subroutine rhoweightxz_amr_weights

subroutine rhoweightr_amr(xp,yp,np,q)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp
REAL(8), INTENT(IN) :: q

REAL(8) :: rpos, ddr, oddr
INTEGER(ISZ) :: i, j, jn, jnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdr

ALLOCATE(invdr(ngrids))

do igrid = 1, ngrids
  invdr(igrid) = grids_ptr(igrid)%grid%invdr
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-g%rmin)*invdr(igrid)
    jn = 1+INT(rpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,1)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,1)
        g=>grids_ptr(igrid)%grid
        rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-g%rmin)*invdr(igrid)
        jn = 1+INT(rpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
    jnp=jn+1
#ifdef MPIPARALLEL
    g%rhop(jn, 1) = g%rhop(jn, 1) + q * oddr * g%invvol(jn)
    g%rhop(jnp,1) = g%rhop(jnp,1) + q *  ddr * g%invvol(jnp)
#else
    g%rho(jn, 1)  = g%rho(jn, 1)  + q * oddr * g%invvol(jn)
    g%rho(jnp,1)  = g%rho(jnp,1)  + q *  ddr * g%invvol(jnp)
#endif
  end do

  DEALLOCATE(invdr)

  return
END subroutine rhoweightr_amr

subroutine rhoweightr_amr_weights(xp,yp,wp,np,q)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, wp
REAL(8), INTENT(IN) :: q

REAL(8) :: rpos, ddr, oddr, qw
INTEGER(ISZ) :: i, j, jn, jnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8), DIMENSION(:), ALLOCATABLE :: invdr

ALLOCATE(invdr(ngrids))

do igrid = 1, ngrids
  invdr(igrid) = grids_ptr(igrid)%grid%invdr
end do

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g=>basegrid
    ingrid=.false.
    rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-g%rmin)*invdr(igrid)
    jn = 1+INT(rpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,1)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,1)
        g=>grids_ptr(igrid)%grid
        rpos = (SQRT(xp(i)*xp(i)+yp(i)*yp(i))-g%rmin)*invdr(igrid)
        jn = 1+INT(rpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
    jnp=jn+1
    qw = q*wp(i)
#ifdef MPIPARALLEL
    g%rhop(jn, 1) = g%rhop(jn, 1) + qw * oddr * g%invvol(jn)
    g%rhop(jnp,1) = g%rhop(jnp,1) + qw *  ddr * g%invvol(jnp)
#else
    g%rho(jn, 1)  = g%rho(jn, 1)  + qw * oddr * g%invvol(jn)
    g%rho(jnp,1)  = g%rho(jnp,1)  + qw *  ddr * g%invvol(jnp)
#endif
  end do

  DEALLOCATE(invdr)

  return
END subroutine rhoweightr_amr_weights

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
                    transit_min_r_cg(ig), &
                    transit_max_r_cg(ig), &
                    transit_min_z_cg(ig), &
                    transit_max_z_cg(ig))
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

subroutine rhoweightrz_deform(xp,yp,zp,np,q,nr,nz,dr,dz,xfact,yfact,rgrid,zgrid)
USE Constant
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), INTENT(IN) :: q, dr, dz, rgrid, zgrid
REAL(8), DIMENSION(0:nz+2), INTENT(INOUT) :: xfact,yfact

REAL(8) :: invdr, invdz, rpos, zpos, ddr, ddz, oddr, oddz, invvol(0:nr)
INTEGER(ISZ) :: i, j, jn, ln, jnp, lnp

  invdr = 1._8/dr
  invdz = 1._8/dz

  ! computes divider by cell volumes to get density
  IF(rgrid==0.) then
    j = 0
    ! the factor 0.75 corrects for overdeposition due to linear weighting (for uniform distribution)
    ! see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
    ! and Verboncoeur, J. of Comp. Phys.,
    invvol(j) = 0.75_8 / (pi * (0.5_8*0.5_8*dr*dr)*dz)
    do j = 1, nr
      invvol(j) = 1._8 / (2._8 * pi * real(j,8) * dr * dr * dz)
    end do
  else
    do j = 0, nr
      invvol(j) = 1._8 / (2._8 * pi * (rgrid+real(j,8)*dr) * dr * dz)
    end do
  END if

  ! make charge deposition using CIC weighting
  do i = 1, np
    zpos = (zp(i)-basegrid%zminp-zgrid)*invdz
    ln = 1+INT(zpos)
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    rpos = (SQRT(xfact(ln)*xp(i)*xfact(ln)*xp(i)+yfact(ln)*yp(i)*yfact(ln)*yp(i))-rgrid)*invdr
    jn = 1+INT(rpos)
    IF(jn>=nr+1.or.ln<1.or.ln>=nz+1) cycle
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
    jnp=jn+1
    lnp=ln+1
#ifdef MPIPARALLEL
    basegrid%rhop(jn, ln)  = basegrid%rhop(jn, ln)  + q * oddr * oddz * invvol(jn-1)
    basegrid%rhop(jnp,ln)  = basegrid%rhop(jnp,ln)  + q *  ddr * oddz * invvol(jnp-1)
    basegrid%rhop(jn, lnp) = basegrid%rhop(jn, lnp) + q * oddr *  ddz * invvol(jn-1)
    basegrid%rhop(jnp,lnp) = basegrid%rhop(jnp,lnp) + q *  ddr *  ddz * invvol(jnp-1)
#else
    basegrid%rho(jn, ln)  = basegrid%rho(jn, ln)  + q * oddr * oddz * invvol(jn-1)
    basegrid%rho(jnp,ln)  = basegrid%rho(jnp,ln)  + q *  ddr * oddz * invvol(jnp-1)
    basegrid%rho(jn, lnp) = basegrid%rho(jn, lnp) + q * oddr *  ddz * invvol(jn-1)
    basegrid%rho(jnp,lnp) = basegrid%rho(jnp,lnp) + q *  ddr *  ddz * invvol(jnp-1)
#endif
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
      write(o_line,*) 'ERROR:periodicity in RZ not yet supported on parallel platform, aborting.'
      call remark(trim(o_line))
      stop
#endif
    IF(ngrids>1) then
      write(o_line,*) 'ERROR:periodicity in RZ not yet supported with mesh refinement, aborting.'
      call remark(trim(o_line))
      stop
    END if
    basegrid%rho(:,1) = basegrid%rho(:,1) + basegrid%rho(:,basegrid%nz+1)
    basegrid%rho(:,basegrid%nz+1) = basegrid%rho(:,1)
  END if

  return
end subroutine rhobndrz

 subroutine perphirz()
 USE multigridrz
 real(kind=8):: substarttime

!  Sets the slices on the exterior of phi for periodicity
!  sets slice at -1 equal to the slice at nz-1
!  sets slice at nz+1 equal to the slice at 1

#ifdef MPIPARALLEL
   call perphi3d_slave(basegrid%phi(1:basegrid%nr+1,0:basegrid%nz+2),basegrid%nr,0,basegrid%nz)
#endif

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

TYPE(GRIDtype) :: grid

INTEGER(ISZ) :: i, ilg, iug, ilp, iup, il, iu, j, l, ii, ll, jex, lex

integer(4):: mpistatus(MPI_STATUS_SIZE)
INTEGER(4):: mpirequest,mpierror

logical(ISZ) :: testthis=.false.
INTEGER(ISZ), ALLOCATABLE :: wz(:)


!write(o_line,*) my_index,'Enter get_rho_from_rhop'
!call remark(trim(o_line))

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
  write(o_line,*) 'rhop = 0.'
  call remark(trim(o_line))
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
    call mpi_isend(grid%rhop(1,il),int(SIZE(grid%rhop(:,il:iu)),4),mpi_double_precision,int(i,4),0_4,MPI_COMM_WORLD,mpirequest,ierr)
!    write(o_line,*) my_index,'send to',i
!    call remark(trim(o_line))
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
!    write(o_line,*) my_index,'recv from',i
!    call remark(trim(o_line))
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
!      IF(jex/wz(l-ilg+ilp)/=j.or.lex/wz(l-ilg+ilp)/=l) write(o_line,*) my_index,':',j,l,l-ilg+ilp,grid%rho(j,ll),wz(l-ilg+ilp)
      IF(jex/wz(l)/=j.or.lex/wz(l)/=l) then
        write(o_line,*) my_index,':',j,l,grid%rho(j,ll),wz(l)
        call remark(trim(o_line))
      END if
 end do
end do
  DEALLOCATE(wz)
 call abort()
endif

if(l_mpi_barriers) call MPI_WAITALL(0_4,mpirequest,mpistatus(1),mpierror)

!write(o_line,*) my_index,'Exit get_rho_from_rhop'
!call remark(trim(o_line))

end subroutine get_rho_from_rhop

subroutine getphiforparticlesrz()
USE multigridrz
implicit none

  call get_phip_from_phi(basegrid)

end subroutine getphiforparticlesrz

subroutine get_phip_from_phi(grid)
USE multigridrz
implicit none

TYPE(GRIDtype) :: grid

INTEGER(ISZ) :: i, ilg, iug, ilp, iup, il, iu, j,l, ll, jex, lex, testeq
integer(4):: mpistatus(MPI_STATUS_SIZE)
INTEGER(4):: mpirequest,mpierror

logical(ISZ) :: testthis=.false.

!write(o_line,*) my_index,'Enter get_phip_from_phi'
!call remark(trim(o_line))

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
    call mpi_isend(grid%phi(0,il),int(SIZE(grid%phi(:,il:iu)),4),mpi_double_precision,int(i,4),0_4,MPI_COMM_WORLD,mpirequest,ierr)
!    write(o_line,*) my_index,'send to',i
!    call remark(trim(o_line))
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
!    write(o_line,*) my_index,'recv from',i
!    call remark(trim(o_line))
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
      IF(jex/=j.or.lex/=l) then
        write(o_line,*) my_index,':',j,l,grid%phip(j,ll)
        call remark(trim(o_line))
      END if
 end do
end do
 call abort()
endif

if(l_mpi_barriers) call MPI_WAITALL(0_4,mpirequest,mpistatus(1),mpierror)

!write(o_line,*) my_index,'Exit get_phip_from_phi'
!call remark(trim(o_line))

end subroutine get_phip_from_phi

#endif

subroutine fieldweightrzold(xp,yp,zp,ex,ey,ez,np,phi,e,nr,nz,dr,dz,zmin,calcselfe,zgrid)
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), DIMENSION(1:2,0:nr,0:nz), INTENT(IN OUT) :: e
REAL(8), INTENT(IN) :: dr, dz, zmin, zgrid
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
    zpos = (zp(i)-zmin-zgrid)*invdz
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

subroutine fieldweightrz(xp,yp,zp,ex,ey,ez,np,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8) :: zgrid

REAL(8) :: r, rpos, zpos, ddr, ddz, oddr, oddz, er
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8):: substarttime

IF(ngrids>1 .and. .not.l_get_field_from_base) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g => basegrid
    ingrid=.false.
    r = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
    rpos = (r-g%rmin)*g%invdr
    zpos = (zp(i)-g%zminp-zgrid)*g%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part_fd(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part_fd(jn,ln)
        g=>grids_ptr(igrid)%grid
        rpos = (r-g%rmin)*g%invdr
        zpos = (zp(i)-g%zminp-zgrid)*g%invdz
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    er = 0.5*(oddr * oddz * (g%phip(jn-1,ln  )-g%phip(jn+1,ln  ))  &
            + ddr  * oddz * (g%phip(jn  ,ln  )-g%phip(jn+2,ln  ))  &
            + oddr * ddz  * (g%phip(jn-1,ln+1)-g%phip(jn+1,ln+1))  &
            + ddr  * ddz  * (g%phip(jn  ,ln+1)-g%phip(jn+2,ln+1)))*g%invdr
#else
    er = 0.5*(oddr * oddz * (g%phi(jn-1,ln  )-g%phi(jn+1,ln  ))  &
            + ddr  * oddz * (g%phi(jn  ,ln  )-g%phi(jn+2,ln  ))  &
            + oddr * ddz  * (g%phi(jn-1,ln+1)-g%phi(jn+1,ln+1))  &
            + ddr  * ddz  * (g%phi(jn  ,ln+1)-g%phi(jn+2,ln+1)))*g%invdr
#endif
    IF(r*g%invdr>1.e-10) then
      ex(i) = er*xp(i)/r
      ey(i) = er*yp(i)/r
    else
      ex(i) = er
      ey(i) = 0._8
    END if
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddr * oddz * (g%phip(jn  ,ln-1)-g%phip(jn  ,ln+1))  &
               + ddr  * oddz * (g%phip(jn+1,ln-1)-g%phip(jn+1,ln+1))  &
               + oddr * ddz  * (g%phip(jn  ,ln  )-g%phip(jn  ,ln+2))  &
               + ddr  * ddz  * (g%phip(jn+1,ln  )-g%phip(jn+1,ln+2)))*g%invdz
#else
    ez(i) = 0.5*(oddr * oddz * (g%phi(jn  ,ln-1)-g%phi(jn  ,ln+1))  &
               + ddr  * oddz * (g%phi(jn+1,ln-1)-g%phi(jn+1,ln+1))  &
               + oddr * ddz  * (g%phi(jn  ,ln  )-g%phi(jn  ,ln+2))  &
               + ddr  * ddz  * (g%phi(jn+1,ln  )-g%phi(jn+1,ln+2)))*g%invdz
#endif
  END do
else
  ! make charge deposition using CIC weighting
  do i = 1, np
    ingrid=.false.
    r = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
    rpos = (r-basegrid%rmin)*basegrid%invdr
    zpos = (zp(i)-basegrid%zminp-zgrid)*basegrid%invdz
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
    IF(r*basegrid%invdr>1.e-10) then
      ex(i) = er*xp(i)/r
      ey(i) = er*yp(i)/r
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

  return
end subroutine fieldweightrz

subroutine fieldweightxz(xp,zp,ex,ez,np,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ez
REAL(8), INTENT(IN) :: zgrid

REAL(8) :: rpos, zpos, invrpos, ddr, ddz, oddr, oddz, zmin0
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

IF( solvergeom==XYgeom) then
  zmin0 = 0.
else ! solvergeom=XZgeom
  zmin0 = zgrid
END if

IF(ngrids>1 .and. .not.l_get_field_from_base) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g => basegrid
    ingrid=.false.
    if(l4symtry) then
      rpos = (ABS(xp(i))-g%xmin)*g%invdr
    else
      rpos = (xp(i)-g%xmin)*g%invdr
    end if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      zpos = (ABS(zp(i))-g%zminp-zmin0)*g%invdz
    else
      zpos = (zp(i)-g%zminp-zmin0)*g%invdz
    END if
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part_fd(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part_fd(jn,ln)
        g=>grids_ptr(igrid)%grid
        if(l4symtry) then
          rpos = (ABS(xp(i))-g%xmin)*g%invdr
        else
          rpos = (xp(i)-g%xmin)*g%invdr
        end if
        IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
          zpos = (ABS(zp(i))-g%zminp-zmin0)*g%invdz
        else
          zpos = (zp(i)-g%zminp-zmin0)*g%invdz
        END if
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    ex(i) = 0.5*(oddr * oddz * (g%phip(jn-1,ln  )-g%phip(jn+1,ln  ))  &
            + ddr  * oddz * (g%phip(jn  ,ln  )-g%phip(jn+2,ln  ))  &
            + oddr * ddz  * (g%phip(jn-1,ln+1)-g%phip(jn+1,ln+1))  &
            + ddr  * ddz  * (g%phip(jn  ,ln+1)-g%phip(jn+2,ln+1)))*g%invdr
#else
    ex(i) = 0.5*(oddr * oddz * (g%phi(jn-1,ln  )-g%phi(jn+1,ln  ))  &
            + ddr  * oddz * (g%phi(jn  ,ln  )-g%phi(jn+2,ln  ))  &
            + oddr * ddz  * (g%phi(jn-1,ln+1)-g%phi(jn+1,ln+1))  &
            + ddr  * ddz  * (g%phi(jn  ,ln+1)-g%phi(jn+2,ln+1)))*g%invdr
#endif
    IF(l4symtry) then
      IF(xp(i)<0.) ex(i) = -ex(i)
    END if
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddr * oddz * (g%phip(jn  ,ln-1)-g%phip(jn  ,ln+1))  &
               + ddr  * oddz * (g%phip(jn+1,ln-1)-g%phip(jn+1,ln+1))  &
               + oddr * ddz  * (g%phip(jn  ,ln  )-g%phip(jn  ,ln+2))  &
               + ddr  * ddz  * (g%phip(jn+1,ln  )-g%phip(jn+1,ln+2)))*g%invdz
#else
    ez(i) = 0.5*(oddr * oddz * (g%phi(jn  ,ln-1)-g%phi(jn  ,ln+1))  &
               + ddr  * oddz * (g%phi(jn+1,ln-1)-g%phi(jn+1,ln+1))  &
               + oddr * ddz  * (g%phi(jn  ,ln  )-g%phi(jn  ,ln+2))  &
               + ddr  * ddz  * (g%phi(jn+1,ln  )-g%phi(jn+1,ln+2)))*g%invdz
#endif
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      IF(zp(i)<0.) ez(i) = -ez(i)
    END if
  END do
else
  ! make charge deposition using CIC weighting
  do i = 1, np
    ingrid=.false.
    if(l4symtry) then
      rpos = (ABS(xp(i))-basegrid%xmin)*basegrid%invdr
    else
      rpos = (xp(i)-basegrid%xmin)*basegrid%invdr
    end if
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      zpos = (ABS(zp(i))-basegrid%zminp-zmin0)*basegrid%invdz
    else
      zpos = (zp(i)-basegrid%zminp-zmin0)*basegrid%invdz
    END if
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
    IF(l4symtry) then
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
    IF((l2symtry .or. l4symtry) .and. solvergeom==XYgeom) then
      IF(zp(i)<0.) ez(i) = -ez(i)
    END if
  END do
END if

  return
end subroutine fieldweightxz

subroutine fieldweightr(xp,yp,ex,ey,np)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey

REAL(8) :: r, rpos, ddr, oddr, er
INTEGER(ISZ) :: i, j, jn, jnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

REAL(8):: substarttime

IF(ngrids>1 .and. .not.l_get_field_from_base) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g => basegrid
    ingrid=.false.
    r = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
    rpos = (r-g%rmin)*g%invdr
    jn = 1+INT(rpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part_fd(jn,1)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part_fd(jn,1)
        g=>grids_ptr(igrid)%grid
        rpos = (r-g%rmin)*g%invdr
        jn = 1+INT(rpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
#ifdef MPIPARALLEL
    er = 0.5*(oddr * (g%phip(jn-1,1 )-g%phip(jn+1,1 ))  &
            + ddr  * (g%phip(jn  ,1 )-g%phip(jn+2,1 ))  )*g%invdr
#else
    er = 0.5*(oddr * (g%phi(jn-1,1 )-g%phi(jn+1,1 ))  &
            + ddr  * (g%phi(jn  ,1 )-g%phi(jn+2,1 ))  )*g%invdr
#endif
    IF(r*g%invdr>1.e-10) then
      ex(i) = er*xp(i)/r
      ey(i) = er*yp(i)/r
    else
      ex(i) = er
      ey(i) = 0._8
    END if
  END do
else
  ! make charge deposition using CIC weighting
  do i = 1, np
    r = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
    rpos = (r-basegrid%rmin)*basegrid%invdr
    jn = 1+INT(rpos)
    ddr = rpos-REAL(jn-1)
    oddr = 1._8-ddr
#ifdef MPIPARALLEL
    er = 0.5*(oddr * (basegrid%phip(jn-1,1  )-basegrid%phip(jn+1,1  ))  &
            + ddr  * (basegrid%phip(jn  ,1  )-basegrid%phip(jn+2,1  ))  )*basegrid%invdr
#else
    er = 0.5*(oddr * (basegrid%phi(jn-1,1  )-basegrid%phi(jn+1,1  ))  &
            + ddr  * (basegrid%phi(jn  ,1  )-basegrid%phi(jn+2,1  ))  )*basegrid%invdr
#endif
    IF(r*basegrid%invdr>1.e-10) then
      ex(i) = er*xp(i)/r
      ey(i) = er*yp(i)/r
    else
      ex(i) = er
      ey(i) = 0._8
    END if
  END do
END if

  return
end subroutine fieldweightr

subroutine fieldweightz(zp,ez,np,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ez
REAL(8), INTENT(IN) :: zgrid

REAL(8) :: zpos, ddz, oddz
INTEGER(ISZ) :: i, l, ln, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g


IF(ngrids>1 .and. .not.l_get_field_from_base) then

  ! make charge deposition using CIC weighting
  do i = 1, np
    igrid = 1
    g => basegrid
    ingrid=.false.
    zpos = (zp(i)-g%zminp-zgrid)*g%invdz
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part_fd(1,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part_fd(1,ln)
        g=>grids_ptr(igrid)%grid
        zpos = (zp(i)-g%zminp-zgrid)*g%invdz
        ln = 1+INT(zpos)
      END if
    end do
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
#ifdef MPIPARALLEL
    ez(i) = 0.5*(oddz * (g%phip(1,ln-1)-g%phip(1,ln+1))  &
               + ddz  * (g%phip(1,ln  )-g%phip(1,ln+2)))*g%invdz
#else
    ez(i) = 0.5*(oddz * (g%phi(1,ln-1)-g%phi(1,ln+1))  &
               + ddz  * (g%phi(1,ln  )-g%phi(1,ln+2)))*g%invdz
#endif
  END do
else
  ! make charge deposition using CIC weighting
  do i = 1, np
    ingrid=.false.
    zpos = (zp(i)-basegrid%zminp-zgrid)*basegrid%invdz
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

  return
end subroutine fieldweightz

subroutine setemgridrz(ipmin,ip,is,ex,ey,ez)
use InGen
use multigridrz
use FRZmgrid
use Particles
use Efields3d
use Picglb

integer(ISZ):: ipmin,ip,is
real(kind=8):: ex(ip),ey(ip),ez(ip)

  if(.not.mgridrz_deform) then
    call fieldweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ex,ey,ez,ip,zgridprv)
  else
    if(is==1 .and. ipmin==ins(is)) call calc_phi3d_from_phirz()
    call sete3d(mgridrz_phi3d(0,0,-1),selfe(1,0,0,0),ip, &
                xp(ipmin),yp(ipmin),zp(ipmin), &
                zgridprv,0.,0.,basegrid%zmin, &
                basegrid%dr,basegrid%dr,basegrid%dz, &
                mgridrz_nx,mgridrz_ny,mgridrz_nz, &
                efetch,ex,ey,ez,l2symtry,l4symtry)
  end if

end subroutine setemgridrz

subroutine calc_phi3d_from_phirz()
USE multigridrz
use FRZmgrid

INTEGER(ISZ) :: j, k, l, jn
REAL(8) :: ddr, oddr, rpos

  do l = -1, mgridrz_nz+1
    do k = 0, mgridrz_ny
      do j = 0, mgridrz_nx
        rpos = SQRT(j**2*mgridrz_xfact(l+1)**2+k**2*mgridrz_yfact(l+1)**2)
        jn = MIN(basegrid%nr,1+INT(rpos))
        ddr=rpos-REAL(jn-1)
        oddr=1._8-ddr
        mgridrz_phi3d(j,k,l) = oddr*basegrid%phi(jn,l) + ddr*basegrid%phi(jn+1,l)
      end do
    end do
  end do

  return
end subroutine calc_phi3d_from_phirz

subroutine fieldweightrz_deform_old2(xp,yp,zp,ex,ey,ez,np,phi,nr,nz,dr,dz,zmin,xfact,yfact,calcphi,phi3d,selfe,zgrid)
USE Constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), INTENT(IN) :: dr, dz, zmin, zgrid
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
!write(o_line,*) 'sum(phi)',SUM(ABS(phi3d)),SUM(ABS(phi))
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

    zpos = (zp(i)-zmin-zgrid)*invdz
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
end subroutine fieldweightrz_deform_old2

subroutine fieldweightrz_deform_old(xp,yp,zp,ex,ey,ez,np,phi,e,nr,nz,dr,dz,zmin,xfact,yfact,calcselfe,zgrid)
USE Constant
implicit none

INTEGER(ISZ), INTENT(IN) :: np, nr, nz
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: ex, ey, ez
REAL(8), DIMENSION(0:nr,-1:nz+1), INTENT(IN) :: phi
REAL(8), DIMENSION(1:2,0:nr,0:nz), INTENT(IN OUT) :: e
REAL(8), INTENT(IN) :: dr, dz, zmin, zgrid
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
    zpos = (zp(i)-zmin-zgrid)*invdz
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
subroutine calcfact_deform(dz,zmin,xfact,yfact,nz,ns,is,ins,nps,ws,zgrid)
USE Constant
USE Particles, ONLY: xp, yp, zp
use FRZmgrid
implicit none

!INTEGER(ISZ), INTENT(IN) :: np, nz, ns
INTEGER(ISZ), INTENT(IN) :: nz, ns
!REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), INTENT(IN) :: dz, zmin, ws(ns), zgrid
REAL(8), DIMENSION(0:nz+2), INTENT(INOUT) :: xfact, yfact
INTEGER(ISZ), DIMENSION(ns), INTENT(IN) :: is, nps, ins

REAL(8), DIMENSION(0:nz+2) :: fweight, fnb
INTEGER(ISZ) :: i, ln, lnp, isp
REAL(8) :: ddz, oddz, wddz, woddz, xrms, yrms, invdz, zpos, xp2, yp2

  if(.not.mgridrz_deform) return

  xfact = 0._8
  yfact = 0._8
  fweight = 0._8
  fnb = 0._8
  invdz = 1._8/dz

!write(o_line,*) 'calcfact',dz,zmin,xfact,yfact,nz,ns,is,ins,nps,ws!isp,ns,ins(is(1)),nps(is(1)),ws(is(1))
!write(o_line,*) 'calcfact',dz,zmin,nz,ns,is,ins,nps,ws!isp,ns,ins(is(1)),nps(is(1)),ws(is(1))
!return
  do isp = 1, ns
    do i = ins(is(isp)), ins(is(isp))+nps(is(isp))-1
      zpos = (zp(i)-zmin-zgrid)*invdz
      ln = 1+INT(zpos)
      lnp = ln+1
      ddz = zpos-REAL(ln-1)
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
 ! write(o_line,*) 'sum ',SUM(fweight)/ws(is(1))
  do ln = 0, nz+2
    IF(fnb(ln)>25._8) then
      xrms = SQRT(xfact(ln)/fweight(ln))
      yrms = SQRT(yfact(ln)/fweight(ln))
 !     write(o_line,*) 'rms',xrms,yrms
      xfact(ln) = 0.5_8*(xrms+yrms)/xrms
      yfact(ln) = 0.5_8*(xrms+yrms)/yrms
    else
      xfact(ln) = 1._8
      yfact(ln) = 1._8
    END if
!    write(o_line,*) ln,xfact(ln),yfact(ln),xrms,yrms,fweight(ln),fnb(ln)
  end do

  return
end subroutine calcfact_deform

subroutine setphirz(np,xp,yp,zp,p,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: xp, yp, zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: p
REAL(8), INTENT(IN) :: zgrid

REAL(8) :: r, rpos, zpos, ddr, ddz, oddr, oddz
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

! Collect phi using linear interpolation

IF(ngrids>1 .and. .not.l_get_injphi_from_base) then

  do i = 1, np
    igrid = 1
    g => basegrid
    r = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
    IF(r<g%rmin.or.r>=g%rmax.or.zp(i)<g%zmin+zgrid.or.zp(i)>=g%zmax+zgrid) cycle
    ingrid=.false.
    rpos = (r-g%rmin)*g%invdr
    zpos = (zp(i)-g%zmin-zgrid)*g%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,ln)
        g=>grids_ptr(igrid)%grid
        IF(r<g%rmin.or.r>=g%rmax.or.zp(i)<g%zmin+zgrid.or.zp(i)>=g%zmax+zgrid) cycle
        rpos = (r-g%rmin)*g%invdr
        zpos = (zp(i)-g%zminp-zgrid)*g%invdz
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    p(i) = oddr * oddz * g%phi(jn,  ln  )  &
         + ddr  * oddz * g%phi(jn+1,ln  )  &
         + oddr * ddz  * g%phi(jn,  ln+1)  &
         + ddr  * ddz  * g%phi(jn+1,ln+1)
  END do
else
  do i = 1, np
    rpos = SQRT(xp(i)*xp(i)+yp(i)*yp(i))
    IF(rpos<basegrid%rmin.or.rpos>=basegrid%rmax.or.zp(i)<basegrid%zmin+zgrid.or.zp(i)>=basegrid%zmax+zgrid) cycle
    ingrid=.false.
    rpos = (rpos-basegrid%rmin)*basegrid%invdr
    zpos = (zp(i)-basegrid%zmin-zgrid)*basegrid%invdz
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

subroutine setphiz(np,zp,p,zgrid)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: np
REAL(8), DIMENSION(np), INTENT(IN) :: zp
REAL(8), DIMENSION(np), INTENT(IN OUT) :: p
REAL(8), INTENT(IN) :: zgrid

REAL(8) :: zpos, ddz, oddz
INTEGER(ISZ) :: i, l, ln, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

! Collect phi using linear interpolation

IF(ngrids>1 .and. .not.l_get_injphi_from_base) then
  do i = 1, np
    igrid = 1
    g => basegrid
    ingrid=.false.
    zpos = (zp(i)-g%zmin-zgrid)*g%invdz
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(1,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(1,ln)
        g=>grids_ptr(igrid)%grid
        zpos = (zp(i)-g%zminp-zgrid)*g%invdz
        ln = 1+INT(zpos)
      END if
    end do
    ddz = zpos-REAL(ln-1)
    oddz = 1._8-ddz
    p(i) = oddz * g%phi(1,ln  )  &
         + ddz  * g%phi(1,ln+1)
  END do
else
  do i = 1, np
    zpos = (zp(i)-basegrid%zmin-zgrid)
    IF(zpos<basegrid%zmin+zgrid.or.zpos>=basegrid%zmax+zgrid) cycle
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

subroutine setphigridrz(p,rmin,zmin,dr,dz,nr,nz,bnd_only)
USE multigridrz
implicit none

INTEGER(ISZ), INTENT(IN) :: nr,nz
REAL(8), DIMENSION(0:nr,0:nz), INTENT(IN OUT) :: p
REAL(8), INTENT(IN) :: rmin, zmin, dr, dz
LOGICAL, INTENT(IN) :: bnd_only

REAL(8) :: r, z, rpos, zpos, ddr, ddz, oddr, oddz
INTEGER(ISZ) :: i, j, l, jn, ln, jnp, lnp, igrid
LOGICAL(ISZ) :: ingrid
TYPE(GRIDtype), pointer :: g

! Collect phi using linear interpolation

IF(ngrids>1 .and. .not.l_get_injphi_from_base) then
  do l = 0, nz!1, nz-1
   z = zmin+l*dz
   do j = 0, nr!1, nr-1
    IF(bnd_only) then
      IF(j>1 .and. j<nr-1 .and. l>1 .and. l<nz-1) cycle
    END if
    r = rmin+j*dr
    igrid = 1
    g => basegrid
!    IF(r<g%rmin.or.r>=g%rmax.or.z<g%zmin.or.z>=g%zmax) cycle
    ingrid=.false.
    ingrid=.true.
    rpos = (r-g%rmin)*g%invdr
    zpos = (z-g%zmin)*g%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    do WHILE(.not.ingrid)
      IF(g%loc_part(jn,ln)==igrid) then
        ingrid=.true.
      else
        igrid = g%loc_part(jn,ln)
        g=>grids_ptr(igrid)%grid
!        IF(r<g%rmin.or.r>=g%rmax.or.z<g%zmin.or.z>=g%zmax) cycle
        rpos = (r-g%rmin)*g%invdr
        zpos = (z-g%zminp)*g%invdz
        jn = 1+INT(rpos)
        ln = 1+INT(zpos)
      END if
    end do
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    p(j,l) = oddr * oddz * g%phi(jn,  ln  )  &
           + ddr  * oddz * g%phi(jn+1,ln  )  &
           + oddr * ddz  * g%phi(jn,  ln+1)  &
           + ddr  * ddz  * g%phi(jn+1,ln+1)
   END do
  END do
else
  do l = 0, nz
   z = zmin+l*dz
   do j = 0, nr
    r = rmin+j*dr
    IF(r<basegrid%rmin.or.r>=basegrid%rmax.or.z<basegrid%zmin.or.z>=basegrid%zmax) cycle
    ingrid=.false.
    rpos = (r-basegrid%rmin)*basegrid%invdr
    zpos = (z-basegrid%zmin)*basegrid%invdz
    jn = 1+INT(rpos)
    ln = 1+INT(zpos)
    ddr = rpos-REAL(jn-1)
    ddz = zpos-REAL(ln-1)
    oddr = 1._8-ddr
    oddz = 1._8-ddz
    p(j,l) = oddr * oddz * basegrid%phi(jn,  ln  )  &
           + ddr  * oddz * basegrid%phi(jn+1,ln  )  &
           + oddr * ddz  * basegrid%phi(jn,  ln+1)  &
           + ddr  * ddz  * basegrid%phi(jn+1,ln+1)
   END do
  END do
END if

  return
end subroutine setphigridrz

!=============================================================================
subroutine setbnd_subgrid_to_inj_d()
use multigridrz
use InjectVars
use InjectVars3d
USE Picglb

INTEGER(ISZ) :: ij, j, l, igrid, il
REAL(8) :: rs, rc, r, z
TYPE(GRIDtype), POINTER :: g

INTEGER(ISZ) :: nconds, n
INTEGER(ISZ), DIMENSION(:), ALLOCATABLE :: ixcond, izcond
TYPE(BNDtype), POINTER :: b

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
       l = MAX(1,2+INT((z-g%zmin-zgrid)/g%dz))
       g%loc_part_fd(j,l:) = g%gid(1)
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
       IF(il == 1) then
         b => g%bndfirst
       else
         b => b%next
       END if
       nconds = 0
       n = (b%nr+1)*(b%nz+1)
       ALLOCATE(ixcond(n),izcond(n))
       do j = 1, b%nr+1
         r = MIN(rs,(j-1)*b%dr)
         ! we assume source emits forward
         z = zinject(ij) + rc - SQRT(rc**2-r**2) + inj_d(ij)*b%dz
         l = 4+MAX(1,2+INT((z-g%zmin-zgrid)/b%dz))
         do l = 4+MAX(1,2+INT((z-g%zmin-zgrid)/b%dz)), b%nz+1
           IF(b%v(j,l)==v_vacuum) then
             b%v(j,l)=v_cond
             nconds=nconds+1
             ixcond(nconds)=j
             izcond(nconds)=l
           END if
         END do
       end do
       IF(nconds>0) then
         call init_bnd_sublevel(b,0,nconds)
         b%cndlast%jcond(1:nconds) = ixcond(1:nconds)
         b%cndlast%kcond(1:nconds) = izcond(1:nconds)
         b%cndlast%voltage(1:nconds) = 0.
         b%cndlast%condid(1:nconds) = -1
        END if
       DEALLOCATE(ixcond,izcond)
     end do
   end do
 end do

return
end subroutine setbnd_subgrid_to_inj_d

!=============================================================================
subroutine set_patches_around_emitter(id,np,ij,nzi,transit_min_r,transit_max_r,transit_min_z,transit_max_z)
use multigridrz
use InjectVars
use InjectVars3d

INTEGER(ISZ), INTENT(IN) :: id, & ! id of grid on which to add the patches
                            np, & ! number of patches
                            ij, & ! id of injection source
                            nzi, & ! size of patches in z (number of meshes)
                            transit_min_r, & ! number of guard cells at lower end in r for field gathering
                            transit_max_r, & ! number of guard cells at upper end in r for field gathering
                            transit_min_z, & ! number of guard cells at lower end in z for field gathering
                            transit_max_z    ! number of guard cells at upper end in z for field gathering

INTEGER(ISZ) :: i, j, l, igrid, nr, nz, l0
REAL(8) :: rs, rc, r, z, dr, dz, rmin, zmin
TYPE(GRIDtype), POINTER :: g
TYPE(BNDtype), POINTER :: b

  nz = nzi

  IF(nz<mgridrz_nmeshmin) then
    nz = mgridrz_nmeshmin
  END if

  rs = ainject(ij)
  rc = rinject(ij)
  g => grids_ptr(id)%grid
  b => g%bndfirst
  dr = g%dr
  dz = g%dz

  do i = 1, np
    l0 = 1+INT((zinject(ij)-b%zmin)/b%dz)
    do j = 1, b%nr+1
      r = (j-1)*b%dr
      IF(r>rs+(2+transit_max_r)*b%dr) exit
      ! we assume source emits forward
      z = zinject(ij) + rc - SQRT(rc**2-r**2)
      l = 1+INT((z-b%zmin)/b%dz)
    END do
    nr = 2*(j-1)
    dr = 0.5*dr
    dz = 0.5*dz
    rmin = 0.
    zmin = b%zmin+INT((zinject(ij)-b%zmin)/b%dz)*b%dz
    write(o_line,*) 'call add_grid'; call remark(trim(o_line))
    write(o_line,*) 'nr = ',nr;      call remark(trim(o_line))
    write(o_line,*) 'nz = ',nz;      call remark(trim(o_line))
    write(o_line,*) 'dr = ',dr;      call remark(trim(o_line))
    write(o_line,*) 'dz = ',dz;      call remark(trim(o_line))
    write(o_line,*) 'zmin = ',zmin;  call remark(trim(o_line))
    call add_grid(g,nr,nz,dr,dz,rmin,zmin,transit_min_r,transit_max_r,transit_min_z,transit_max_z)
    g => grids_ptr(ngrids)%grid
  end do

return
end subroutine set_patches_around_emitter

subroutine clean_conductor_interior()
USE multigridrz
implicit none

INTEGER(ISZ) :: i, ic, il, igrid, ix, iz, ncond
TYPE(GRIDtype), POINTER :: g
TYPE(BNDtype), POINTER :: b
TYPE(CONDtype), POINTER :: c

  do igrid=1,ngrids
    g => grids_ptr(igrid)%grid
    do il = 1, g%nlevels
      IF(il == 1) then
        b => g%bndfirst
      else
        b => b%next
      END if
      do ic = 1, b%nb_conductors
        IF(ic==1) then
          c => b%cndfirst
        else
          c => c%next
        END if
        ncond = 0
        do i = 1, c%ncond
          ix = c%jcond(i)
          iz = c%kcond(i)
          IF(.NOT.(ABS(b%v(ix+1,iz+1))==v_cond .AND. &
                   ABS(b%v(ix-1,iz+1))==v_cond .AND. &
                   ABS(b%v(ix+1,iz-1))==v_cond .AND. &
                   ABS(b%v(ix-1,iz-1))==v_cond)) then
            ncond = ncond + 1
            IF(i /= ncond) then
              c%jcond(ncond)   = c%jcond(i)
              c%kcond(ncond)   = c%kcond(i)
              c%voltage(ncond) = c%voltage(i)
              c%condid(ncond)  = c%condid(i)
            END if
          else
            b%v(ix,iz) = -v_cond
          END if
        END do
        c%ncond = ncond
      END do
    END do
  END do

  return
end subroutine clean_conductor_interior

subroutine build_vlocs()
USE multigridrz
implicit none

INTEGER(ISZ) :: il, igrid, j, l, ilred, ilblack, jmin, jmax, lmin, lmax
TYPE(GRIDtype), POINTER :: g
TYPE(BNDtype), POINTER :: b

  vlocs = .true.

  do igrid=1,ngrids
    g => grids_ptr(igrid)%grid
    g%npmin = 1
    do il = 1, g%nlevels
      IF(il == 1) then
        b => g%bndfirst
      else
        b => b%next
      END if
      b%nvlocs = 0
      b%nvlocsred = 0
      do l = 1, b%nz+1
        do j = 1, b%nr+1
          IF(b%v(j,l)==v_vacuum) then
            IF(MOD(j+l,2)==0) b%nvlocsred = b%nvlocsred + 1
            b%nvlocs = b%nvlocs + 1
          END if
        end do
      end do
!      ALLOCATE(b%vlocs_j(b%nvlocs),b%vlocs_k(b%nvlocs))
      call BNDtypechange(b)
      ilred = 0
      ilblack = b%nvlocsred
      jmax = 1
      jmin = b%nr+1
      lmax = 1
      lmin = b%nz+1
      do l = 1, b%nz+1
        do j = 1, b%nr+1
          IF(b%v(j,l)==v_vacuum) then
            jmin = MIN(jmin,j)
            jmax = MAX(jmax,j)
            lmin = MIN(lmin,l)
            lmax = MAX(lmax,l)
            IF(MOD(j+l,2)==0) then
              ilred = ilred + 1
              b%vlocs_j(ilred) = j
              b%vlocs_k(ilred) = l
            else
              ilblack = ilblack + 1
              b%vlocs_j(ilblack) = j
              b%vlocs_k(ilblack) = l
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

subroutine set_basegrid()
USE multigridrz
  IF(associated(basegrid)) return

  basegrid => NewGRIDType()

  return
end subroutine set_basegrid

subroutine nullify_basegrid()
USE multigridrz
  NULLIFY(basegrid)
  return
end subroutine nullify_basegrid

subroutine del_base()
USE multigridrz
implicit none

  IF(.NOT.associated(basegrid)) return
  IF(associated(basegrid)) then
    call del_grid(basegrid)
    NULLIFY(basegrid)
  endif

return
end subroutine del_base

subroutine mk_grids_ptr()
USE multigridrz
implicit none
INTEGER :: i

  IF(ALLOCATED(grids_ptr)) DEALLOCATE(grids_ptr)
  ALLOCATE(grids_ptr(ngrids))
  do i = 1, ngrids
    NULLIFY(grids_ptr(i)%grid)
  end do
  call assign_grids_ptr(basegrid,.true.)
!  write(o_line,*) 'call gchange ', ngrids
!  call remark(trim(o_line))
  call gchange("FRZmgrid",0)
!  write(o_line,*) 'done.'
!  call remark(trim(o_line))
  do i = 1, ngrids
    nrg(i) = grids_ptr(i)%grid%nr
    nzg(i) = SIZE(grids_ptr(i)%grid%loc_part,2)-1
    drg(i) = grids_ptr(i)%grid%dr
    dzg(i) = grids_ptr(i)%grid%dz
  END do

  return
end subroutine mk_grids_ptr

subroutine add_subgrid(id,nr,nz,dr,dz,rmin,zmin,transit_min_r,transit_max_r,transit_min_z,transit_max_z)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz,transit_min_r,transit_max_r,transit_min_z,transit_max_z
REAL(8), INTENT(IN) :: dr,dz,rmin,zmin

  IF(id<1 .or. id>ngrids) then
    write(o_line,*) 'Fatal error in add_subgrid: id = ', id ,' WHILE id = (1,..,',ngrids,')'
    call remark(trim(o_line))
    stop
  END if

  call add_grid(grids_ptr(id)%grid,nr,nz,dr,dz,rmin,zmin,transit_min_r,transit_max_r,transit_min_z,transit_max_z)

return
END subroutine add_subgrid

subroutine del_subgrid(id)
USE multigridrz
implicit none
INTEGER, INTENT(IN) :: id

  IF(id<1 .or. id>ngrids) then
    write(o_line,*) 'Fatal error in del_subgrid: id = ', id ,' WHILE id = (1,..,',ngrids,')'
    call remark(trim(o_line))
    return
  END if

  call del_grid(grids_ptr(id)%grid)
  NULLIFY(grids_ptr(id)%grid)

  return
END subroutine del_subgrid

subroutine del_conductors()
USE multigridrz
implicit none
TYPE(GRIDtype), pointer :: g
TYPE(BNDtype), POINTER :: b,bnext
INTEGER :: id

  do id = 1, ngrids
    g => grids_ptr(id)%grid
    b=>g%bndfirst
    call del_cnds(b)
    do WHILE(associated(b%next))
      b => b%next
      call del_cnds(b)
    end do
  end do

  return
end subroutine del_conductors

subroutine get_phi_subgrid(id,phi,nr,nz)
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: id,nr,nz
REAL(8), DIMENSION(1:nr+1,1:nz+1) :: phi

  IF(id<1 .or. id>ngrids) then
    write(o_line,*) 'Error in get_phi_subgrid: id = ', id ,' WHILE id = (1,..,',ngrids,')'
    call remark(trim(o_line))
    write(o_line,*) 'Returning Phi=0'
    call remark(trim(o_line))
    phi(1:nr+1,1:nz+1) = 0.
  else
    phi(1:nr+1,1:nz+1) = grids_ptr(id)%grid%phi(1:nr+1,1:nz+1)
  END if

return
END subroutine get_phi_subgrid

subroutine set_basegrid_phi()
USE multigridrz
USE Fields3d
implicit none

  basegrid%phi(1:basegrid%nr+1,:) = phi(:,0,:)

return
END subroutine set_basegrid_phi

subroutine setmglevels_rz(grid)
USE multigridrz
USE MGLevels3d
implicit none
INTEGER :: i,mglevel
TYPE(GRIDtype) :: grid
TYPE(BNDtype), pointer :: b

  mglevels = grid%nlevels
  do i = 1, grid%nlevels
    IF(i==1) then
      b => grid%bndfirst
    else
      b => b%next
    END if
    mglevel = i-1
    mglevelsnx(mglevel) = b%nr
    mglevelslx(mglevel) = b%dr/grid%bndfirst%dr
    if (solvergeom==XYgeom) then
      mglevelsny(mglevel) = b%nz
      mglevelsnzfull(mglevel) = 0
      mglevelsnz(mglevel) = 0
      mglevelsly(mglevel) = b%dz/grid%bndfirst%dz
      mglevelslz(mglevel) = 1.
    else
      mglevelsny(mglevel) = 0
#ifdef MPIPARALLEL
      mglevelsnzfull(mglevel) = b%nz*nslaves/b%nworkpproc
      mglevelsiz(mglevel) = INT(my_index/b%nworkpproc)*b%nz!-1
#else
      mglevelsnzfull(mglevel) = b%nz
      mglevelsiz(mglevel) = 0
#endif
      mglevelsnz(mglevel) = b%nz
      mglevelsly(mglevel) = 1.
      mglevelslz(mglevel) = b%dz/grid%bndfirst%dz
    end if
  END do

  return
end subroutine setmglevels_rz

subroutine get_cond_rz_level(igrid,ilevel)
USE multigridrz
USE Conductor3d
implicit none
INTEGER :: igrid,ilevel

INTEGER :: i,ic,icc,ice,ico
TYPE(BNDtype), pointer :: bnd
TYPE(CONDtype), pointer :: c

 IF(solvergeom==Zgeom .or. solvergeom==Rgeom) return

 IF(ilevel>grids_ptr(igrid)%grid%nlevels) then
   WRITE(0,*) 'Error in get_condrz: ilevel>nlevels.'
   return
 END if

 bnd => grids_ptr(igrid)%grid%bndfirst
 do i = 1, ilevel-1
   bnd => bnd%next
 end do
 
 conductors%interior%n = 0
 conductors%evensubgrid%n = 0
 conductors%oddsubgrid%n = 0
 do ic = 1, bnd%nb_conductors
   IF(ic==1) then
     c => bnd%cndfirst
   else
     c => c%next
   END if
   conductors%interior%n    = conductors%interior%n    + c%ncond
   conductors%evensubgrid%n = conductors%evensubgrid%n + c%nbbndred
   conductors%oddsubgrid%n  = conductors%oddsubgrid%n  + c%nbbnd-c%nbbndred
 END do
 if (conductors%interior%nmax < conductors%interior%n) &
     conductors%interior%nmax = conductors%interior%n
 if (conductors%evensubgrid%nmax < conductors%evensubgrid%n) &
     conductors%evensubgrid%nmax = conductors%evensubgrid%n
 if (conductors%oddsubgrid%nmax < conductors%oddsubgrid%n) &
     conductors%oddsubgrid%nmax = conductors%oddsubgrid%n
 call gchange("Conductor3d",0)

 icc=0
 ice=0
 ico=0
 do ic = 1, bnd%nb_conductors
   IF(ic==1) then
     c => bnd%cndfirst
   else
     c => c%next
   END if
   do i = 1, c%ncond
     icc=icc+1
     conductors%interior%indx(0,icc) = c%jcond(i)-1
     conductors%interior%indx(2,icc) = c%kcond(i)-1
     conductors%interior%ilevel(icc) = ilevel - 1
     conductors%interior%volt(icc) = c%voltage(i)
   end do
   do i = 1, c%nbbndred
    IF(bnd%v(c%jj(i),c%kk(i))==v_bnd) then
     ice=ice+1
     conductors%evensubgrid%indx(0,ice) = c%jj(i)-1
     conductors%evensubgrid%indx(2,ice) = c%kk(i)-1
     conductors%evensubgrid%dels(0,ice) = c%dxm(i)/bnd%dr
     conductors%evensubgrid%dels(1,ice) = c%dxp(i)/bnd%dr
     conductors%evensubgrid%dels(4,ice) = c%dzm(i)/bnd%dz
     conductors%evensubgrid%dels(5,ice) = c%dzp(i)/bnd%dz
     conductors%evensubgrid%ilevel(ice) = ilevel - 1
     conductors%evensubgrid%volt(0,ice) = c%volt0xm(i)
     conductors%evensubgrid%volt(1,ice) = c%volt0xp(i)
     conductors%evensubgrid%volt(4,ice) = c%volt0zm(i)
     conductors%evensubgrid%volt(5,ice) = c%volt0zp(i)
    END if
   end do
   do i = c%nbbndred+1, c%nbbnd
    IF(bnd%v(c%jj(i),c%kk(i))==v_bnd) then
     ico=ico+1
     conductors%oddsubgrid%indx(0,ico) = c%jj(i)-1
     conductors%oddsubgrid%indx(2,ico) = c%kk(i)-1
     conductors%oddsubgrid%dels(0,ico) = c%dxm(i)/bnd%dr
     conductors%oddsubgrid%dels(1,ico) = c%dxp(i)/bnd%dr
     conductors%oddsubgrid%dels(4,ico) = c%dzm(i)/bnd%dz
     conductors%oddsubgrid%dels(5,ico) = c%dzp(i)/bnd%dz
     conductors%oddsubgrid%ilevel(ico) = ilevel - 1
     conductors%oddsubgrid%volt(0,ico) = c%volt0xm(i)
     conductors%oddsubgrid%volt(1,ico) = c%volt0xp(i)
     conductors%oddsubgrid%volt(4,ico) = c%volt0zm(i)
     conductors%oddsubgrid%volt(5,ico) = c%volt0zp(i)
    END if
   end do
 END do
 conductors%evensubgrid%n = ice
 conductors%oddsubgrid%n = ico

return
end subroutine get_cond_rz_level

subroutine get_cond_rz(igrid)
USE multigridrz
USE Conductor3d
implicit none
INTEGER :: igrid

INTEGER :: i,il,ic,icc,ice,ico,kl
TYPE(BNDtype), pointer :: bnd
TYPE(CONDtype), pointer :: c

 IF(solvergeom==Zgeom .or. solvergeom==Rgeom) return

 IF(solvergeom==XYgeom) then
   kl = 1
 else
   kl = 2
 END if

 conductors%interior%n = 0
 conductors%evensubgrid%n = 0
 conductors%oddsubgrid%n = 0
 do il = 1, grids_ptr(igrid)%grid%nlevels
   IF(il == 1) then
     bnd => grids_ptr(igrid)%grid%bndfirst
   else
     bnd => bnd%next
   END if
   do ic = 1, bnd%nb_conductors
     IF(ic==1) then
       c => bnd%cndfirst
     else
       c => c%next
     END if
     conductors%interior%n    = conductors%interior%n    + c%ncond
     conductors%evensubgrid%n = conductors%evensubgrid%n + c%nbbndred
     conductors%oddsubgrid%n  = conductors%oddsubgrid%n  + c%nbbnd-c%nbbndred
   END do
 END do
 if (conductors%interior%nmax < conductors%interior%n) &
     conductors%interior%nmax = conductors%interior%n
 if (conductors%evensubgrid%nmax < conductors%evensubgrid%n) &
     conductors%evensubgrid%nmax = conductors%evensubgrid%n
 if (conductors%oddsubgrid%nmax < conductors%oddsubgrid%n) &
     conductors%oddsubgrid%nmax = conductors%oddsubgrid%n
 call gchange("Conductor3d",0)

 icc=0
 ice=0
 ico=0
 do il = 1, grids_ptr(igrid)%grid%nlevels
   IF(il == 1) then
     bnd => grids_ptr(igrid)%grid%bndfirst
   else
     bnd => bnd%next
   END if
   do ic = 1, bnd%nb_conductors
     IF(ic==1) then
       c => bnd%cndfirst
     else
       c => c%next
     END if
     do i = 1, c%ncond
       icc=icc+1
       conductors%interior%indx(0 ,icc) = c%jcond(i)-1
       conductors%interior%indx(kl,icc) = c%kcond(i)-1
       conductors%interior%ilevel(icc)  = il - 1
       conductors%interior%volt(icc)    = c%voltage(i)
     end do
     do i = 1, c%nbbndred
      IF(bnd%v(c%jj(i),c%kk(i))==v_bnd) then
       ice=ice+1
       conductors%evensubgrid%indx(0     ,ice) = c%jj(i)-1
       conductors%evensubgrid%indx(kl    ,ice) = c%kk(i)-1
       conductors%evensubgrid%dels(0     ,ice) = c%dxm(i)/bnd%dr
       conductors%evensubgrid%dels(1     ,ice) = c%dxp(i)/bnd%dr
       conductors%evensubgrid%dels(2*kl  ,ice) = c%dzm(i)/bnd%dz
       conductors%evensubgrid%dels(2*kl+1,ice) = c%dzp(i)/bnd%dz
       conductors%evensubgrid%ilevel(ice)      = il - 1
       conductors%evensubgrid%volt(0    ,ice) = c%volt0xm(i)
       conductors%evensubgrid%volt(1    ,ice) = c%volt0xp(i)
       conductors%evensubgrid%volt(2*kl ,ice) = c%volt0zm(i)
       conductors%evensubgrid%volt(2*kl+1,ice) = c%volt0zp(i)
      END if
     end do
     do i = c%nbbndred+1, c%nbbnd
      IF(bnd%v(c%jj(i),c%kk(i))==v_bnd) then
       ico=ico+1
       conductors%oddsubgrid%indx(0     ,ico) = c%jj(i)-1
       conductors%oddsubgrid%indx(kl    ,ico) = c%kk(i)-1
       conductors%oddsubgrid%dels(0     ,ico) = c%dxm(i)/bnd%dr
       conductors%oddsubgrid%dels(1     ,ico) = c%dxp(i)/bnd%dr
       conductors%oddsubgrid%dels(2*kl  ,ico) = c%dzm(i)/bnd%dz
       conductors%oddsubgrid%dels(2*kl+1,ico) = c%dzp(i)/bnd%dz
       conductors%oddsubgrid%ilevel(ico)      = il - 1
       conductors%oddsubgrid%volt(0     ,ico) = c%volt0xm(i)
       conductors%oddsubgrid%volt(1     ,ico) = c%volt0xp(i)
       conductors%oddsubgrid%volt(2*kl  ,ico) = c%volt0zm(i)
       conductors%oddsubgrid%volt(2*kl+1,ico) = c%volt0zp(i)
      END if
     end do
   END do
 END do
 conductors%evensubgrid%n = ice
 conductors%oddsubgrid%n = ico

return
end subroutine get_cond_rz

subroutine setconductorvoltagerz(volt,nz,zmmin,dz,discrete,id)
USE multigridrz
implicit none
integer(ISZ):: nz
real(kind=8):: volt(0:nz)
real(kind=8):: rmmin,zmmin,dz
logical(ISZ):: discrete
integer(ISZ):: id

! If the id is zero, then the voltage is applied to all conductors.
! Otherwise, it is applied only to the one specified.

INTEGER :: igrid,i,iv,ic,icc,ice,ico
integer(ISZ):: iz
real(kind=8):: zz,wz,vv
TYPE(CONDtype), POINTER :: c
real(kind=8):: dxm,dxp,dzm,dzp,dxx,dzz,r,rm,rp
TYPE(BNDtype), POINTER :: b
LOGICAL(ISZ) :: l_change

do igrid=1,ngrids
  nlevels=grids_ptr(igrid)%grid%nlevels
  rmmin = grids_ptr(igrid)%grid%rmin
  do i = 1, nlevels
   IF(i == 1) then
     b => grids_ptr(igrid)%grid%bndfirst
   else
     b => b%next
   END if
   do iv=1, b%nb_conductors
     IF(iv==1) then
       c => b%cndfirst
     else
       c => c%next
     END if

    do ic=1,c%ncond
      IF(c%condid(ic) /= id .and. id /= 0) cycle
      zz = grids_ptr(igrid)%grid%zmin + b%dz*(c%kcond(ic)-1)
      if (zmmin <= zz .and. zz < zmmin + nz*dz) then
        iz = int(zz/dz)
        wz =     zz/dz - iz
        c%voltage(ic) = volt(iz)*(1.-wz) + volt(iz+1)*wz
      else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + b%dz) then
        c%voltage(ic) = volt(nz)
      endif
    enddo

    do ic = 1,c%nbbnd

      l_change = .false.
      zz = grids_ptr(igrid)%grid%zmin + b%dz*(c%kk(ic)-1)
      if (zmmin <= zz .and. zz < zmmin + nz*dz) then
        iz = int(zz/dz)
        wz =     zz/dz - iz
        vv = volt(iz)*(1.-wz) + volt(iz+1)*wz
        if (c%dxm(ic) < b%dr .and. (c%condidxm(ic)==id .or. id == 0)) then
          l_change = .true.
          c%volt0xm(ic) = vv
        endif
        if (c%dxp(ic) < b%dr .and. (c%condidxp(ic)==id .or. id == 0)) then
          l_change = .true.
          c%volt0xp(ic) = vv
        endif
      else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + b%dz) then
        vv = volt(nz)
        if (c%dxm(ic) < b%dr .and. (c%condidxm(ic)==id .or. id == 0)) then
          l_change = .true.
          c%volt0xm(ic) = vv
        endif
        if (c%dxp(ic) < b%dr .and. (c%condidxp(ic)==id .or. id == 0)) then
          l_change = .true.
          c%volt0xp(ic) = vv
        endif
      endif
      if (c%dzm(ic) < b%dz .and. (c%condidzm(ic)==id .or. id == 0)) then
        zz = grids_ptr(igrid)%grid%zmin + b%dz*(c%kk(ic)-1) &
             - c%dzm(ic)
        if (zmmin <= zz .and. zz < zmmin + nz*dz) then
          iz = int(zz/dz)
          wz =     zz/dz - iz
          if (discrete) wz = 0.
          l_change = .true.
          c%volt0zm(ic) = volt(iz)*(1.-wz) + volt(iz+1)*wz
        else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + b%dz) then
          l_change = .true.
          c%volt0zm(ic) = volt(nz)
        endif
      endif
      if (c%dzp(ic) < b%dz .and. (c%condidzp(ic)==id .or. id == 0)) then
        zz = grids_ptr(igrid)%grid%zmin + b%dz*(c%kk(ic)-1) &
             + c%dzp(ic)
        if (zmmin <= zz .and. zz < zmmin + nz*dz) then
          iz = int(zz/dz)
          wz =     zz/dz - iz
          if (discrete) wz = 1.
          l_change = .true.
          c%volt0zp(ic) = volt(iz)*(1.-wz) + volt(iz+1)*wz
        else if (zmmin + nz*dz <= zz .and. zz < zmmin + nz*dz + b%dz) then
          l_change = .true.
          c%volt0zp(ic) = volt(nz)
        endif
      endif

      if (.not. l_change) cycle

      dxm = MIN(b%dr,c%dxm(ic))
      dxp = MIN(b%dr,c%dxp(ic))
      dzm = MIN(b%dz,c%dzm(ic))
      dzp = MIN(b%dz,c%dzp(ic))
      select case (bnd_method)
        case (egun)
          dxx=b%dr
          dzz=b%dz
        case (ecb)
          dxx=0.5_8*(dxp+dxm)  !ecb
          dzz=0.5_8*(dzp+dzm)  !ecb
        case default
      end select
      IF(solvergeom==RZgeom) then
       IF(c%jj(ic)==1 .and. rmmin==0.) then
        c%cfxp(ic) = 4._8/(dxp*dxx)
       else
        r = rmmin+(c%jj(ic)-1)*b%dr
        select case (bnd_method)
          case (egun)
            rm = r-0.5_8*b%dr
            rp = r+0.5_8*b%dr
          case (ecb)
            rm = r-0.5_8*dxm
            rp = r+0.5_8*dxp
          case default
        end select
        c%cfxm(ic) = rm/(r*dxm*dxx)
        c%cfxp(ic) = rp/(r*dxp*dxx)
       END if
      else !IF(solvergeom==XZgeom) then
        c%cfxm(ic) = 1._8/(dxm*dxx)
        c%cfxp(ic) = 1._8/(dxp*dxx)
      END if
      c%cfzm(ic) = 1._8/(dzm*dzz)
      c%cfzp(ic) = 1._8/(dzp*dzz)
      IF(dxm>=b%dr) then
        c%phi0xm(ic)=0._8
      else
        c%phi0xm(ic)=c%cfxm(ic)*c%volt0xm(ic)
        c%cfxm(ic)=0._8
      END if
      IF(dxp>=b%dr) then
        c%phi0xp(ic)=0._8
      else
        c%phi0xp(ic)=c%cfxp(ic)*c%volt0xp(ic)
        c%cfxp(ic)=0._8
      END if
      IF(dzm>=b%dz) then
        c%phi0zm(ic)=0._8
      else
        c%phi0zm(ic)=c%cfzm(ic)*c%volt0zm(ic)
        c%cfzm(ic)=0._8
      END if
      IF(dzp>=b%dz) then
        c%phi0zp(ic)=0._8
      else
        c%phi0zp(ic)=c%cfzp(ic)*c%volt0zp(ic)
        c%cfzp(ic)=0._8
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
real(kind=8):: dxm,dxp,dzm,dzp,dxx,dzz,r,rm,rp,rmmin
LOGICAL(ISZ) :: l_change
TYPE(BNDtype), POINTER :: b
TYPE(CONDtype), POINTER :: c

do igrid=1,ngrids
  nlevels=grids_ptr(igrid)%grid%nlevels
  rmmin = grids_ptr(igrid)%grid%rmin
  do i = 1, nlevels
   IF(i == 1) then
     b => grids_ptr(igrid)%grid%bndfirst
   else
     b => b%next
   END if
   do iv=1, b%nb_conductors
     IF(iv==1) then
       c => b%cndfirst
     else
       c => c%next
     END if

    do ic=1,c%ncond
      IF(c%condid(ic)==id) c%voltage(ic) = volt
    enddo

    do ic = 1,c%nbbnd
     l_change = .false.
     if (c%dxm(ic) < b%dr .and. c%condidxm(ic)==id) then
      l_change = .true.
      c%volt0xm(ic) = volt
     END if
     if (c%dxp(ic) < b%dr .and. c%condidxp(ic)==id) then
      l_change = .true.
      c%volt0xp(ic) = volt
     END if
     if (c%dzm(ic) < b%dz .and. c%condidzm(ic)==id) then
      l_change = .true.
      c%volt0zm(ic) = volt
     END if
     if (c%dzp(ic) < b%dz .and. c%condidzp(ic)==id) then
      l_change = .true.
      c%volt0zp(ic) = volt
     END if

     IF(l_change) then
      dxm = MIN(b%dr,c%dxm(ic))
      dxp = MIN(b%dr,c%dxp(ic))
      dzm = MIN(b%dz,c%dzm(ic))
      dzp = MIN(b%dz,c%dzp(ic))
      select case (bnd_method)
        case (egun)
          dxx=b%dr
          dzz=b%dz
        case (ecb)
          dxx=0.5_8*(dxp+dxm)  !ecb
          dzz=0.5_8*(dzp+dzm)  !ecb
        case default
      end select
      IF(solvergeom==RZgeom) then
       IF(c%jj(ic)==1 .and. rmmin==0.) then
        c%cfxp(ic) = 4._8/(dxp*dxx)
       else
        r = rmmin+(c%jj(ic)-1)*b%dr
        select case (bnd_method)
          case (egun)
            rm = r-0.5_8*b%dr
            rp = r+0.5_8*b%dr
          case (ecb)
            rm = r-0.5_8*dxm
            rp = r+0.5_8*dxp
          case default
        end select
        c%cfxm(ic) = rm/(r*dxm*dxx)
        c%cfxp(ic) = rp/(r*dxp*dxx)
       END if
      else !IF(solvergeom==XZgeom) then
        c%cfxm(ic) = 1._8/(dxm*dxx)
        c%cfxp(ic) = 1._8/(dxp*dxx)
      END if
      c%cfzm(ic) = 1._8/(dzm*dzz)
      c%cfzp(ic) = 1._8/(dzp*dzz)
      IF(dxm>=b%dr) then
        c%phi0xm(ic)=0._8
      else
        c%phi0xm(ic)=c%cfxm(ic)*c%volt0xm(ic)
        c%cfxm(ic)=0._8
      END if
      IF(dxp>=b%dr) then
        c%phi0xp(ic)=0._8
      else
        c%phi0xp(ic)=c%cfxp(ic)*c%volt0xp(ic)
        c%cfxp(ic)=0._8
      END if
      IF(dzm>=b%dz) then
        c%phi0zm(ic)=0._8
      else
        c%phi0zm(ic)=c%cfzm(ic)*c%volt0zm(ic)
        c%cfzm(ic)=0._8
      END if
      IF(dzp>=b%dz) then
        c%phi0zp(ic)=0._8
      else
        c%phi0zp(ic)=c%cfzp(ic)*c%volt0zp(ic)
        c%cfzp(ic)=0._8
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
!TYPE(GRIDtype), pointer :: g

INTEGER(ISZ):: i

  IF(.not. l_change_loc_part) return

  do i = 1, nz_rmc+1
    basegrid%loc_part_fd(1:rmc(i)-1,i) = 1
    basegrid%loc_part_fd(rmc(i):,i) = 2
  end do
  l_change_loc_part = .false.

END subroutine change_loc_part

subroutine adjust_lpfd(f,nr,nz,rmin,rmax,zmin,zmax)
! adjusts loc_part_fd according to data used to setup grids.
USE multigridrz
implicit none
INTEGER(ISZ), INTENT(IN) :: nr, nz
REAL(8), INTENT(IN) :: f(nr+1,nz+1), rmin, rmax, zmin, zmax

TYPE(GRIDtype), pointer :: g
INTEGER(ISZ) :: i, j, l, r, jf, lf
REAL(8) :: dr, dz

  dr = (rmax-rmin)/nr
  dz = (zmax-zmin)/nz

  do i = 1, ngrids
    if (i==1) then
      g => basegrid
    else
      IF(associated(g%next)) then
        g => g%next
      else
        g => g%down
      END if
    END if
    r = INT(dr/g%dr)
    do l = 1, g%nz+1
      lf = 1+INT(((g%zmin+(l-1)*g%dz+0.5*g%dz)-zmin)/dz)
      do j = 1, g%nr+1
        jf = 1+INT(((g%rmin+(j-1)*g%dr+0.5*g%dr)-rmin)/dr)
        IF(g%loc_part_fd(j,l)/=g%gid(1) .AND. f(jf,lf)==r) g%loc_part_fd(j,l)=g%gid(1)
      end do
    end do
  end do
end subroutine adjust_lpfd

subroutine sum_neighbors(fin,fout,nx,ny)
INTEGER(ISZ), INTENT(IN) :: nx, ny
INTEGER(ISZ) :: fin(nx+1,ny+1)
INTEGER(ISZ) :: fout(nx+1,ny+1)

INTEGER(ISZ) :: j,k

  do k = 2, ny
    do j = 2, nx
      fout(j,k) = SUM(fin(j-1:j+1,k-1:k+1))
    end do
  end do
  j = 1
  do k = 1, ny+1
    fout(j,k) = SUM(fin(MAX(1,j-1):MIN(nx+1,j+1),MAX(1,k-1):MIN(ny+1,k+1)))
  end do
  j = nx+1
  do k = 1, ny+1
    fout(j,k) = SUM(fin(MAX(1,j-1):MIN(nx+1,j+1),MAX(1,k-1):MIN(ny+1,k+1)))
  end do
  k = 1
  do j = 2, nx
    fout(j,k) = SUM(fin(MAX(1,j-1):MIN(nx+1,j+1),MAX(1,k-1):MIN(ny+1,k+1)))
  end do
  k = ny+1
  do j = 2, nx
    fout(j,k) = SUM(fin(MAX(1,j-1):MIN(nx+1,j+1),MAX(1,k-1):MIN(ny+1,k+1)))
  end do

end subroutine sum_neighbors

