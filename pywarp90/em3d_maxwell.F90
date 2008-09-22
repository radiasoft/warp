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
use EM3D_KYEEFIELDtypemodule
use EM3D_SPLITYEEFIELDtypemodule
USE mod_bnd
USE EM2D_FIELDobjects
!use GlobalVars
!use Picglb

implicit none

integer, parameter :: dirichlet=0, neumann=1, periodic=2, openbc=3, yeefield=-1 , splityeefield=-2
#ifndef MPIPARALLEL
integer, parameter :: my_index=0
#endif

!#ifdef MPIPARALLEL
!include "mpif.h"
!integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpierror
!integer(MPIISZ):: mpirequest
!integer(MPIISZ):: w
!integer(MPIISZ):: messid 
!#endif

TYPE bnd_pointer
  type(type_bnd), POINTER :: b
end type bnd_pointer
TYPE bnd_cummer_pointer
  type(type_bnd_cummer), POINTER :: b
end type bnd_cummer_pointer
!type(bnd_pointer), dimension(3,2) :: bnds ! first dimension is for [main grid, coarse patch, fine patch]
!                                          ! second dimension is for [(Ex,Ey,Bz),(Bx,By,Ez)]
!INTEGER, parameter :: base=1, patchcoarse=2, patchfine=3

contains

  subroutine set_bndcoeffsem3d(sf,dt,which)
    use mod_bnd
    TYPE(EM3D_SPLITYEEFIELDtype) :: sf
    real(kind=8) :: sigmax(-sf%nxguard:sf%nx+sf%nxguard), sigmax_next(-sf%nxguard:sf%nx+sf%nxguard), &
                  lsigmax(-sf%nxguard:sf%nx+sf%nxguard), lsigmax_next(-sf%nxguard:sf%nx+sf%nxguard)
    real(kind=8) :: sigmay(-sf%nyguard:sf%ny+sf%nyguard), sigmay_next(-sf%nyguard:sf%ny+sf%nyguard), &
                  lsigmay(-sf%nyguard:sf%ny+sf%nyguard), lsigmay_next(-sf%nyguard:sf%ny+sf%nyguard)
    real(kind=8) :: sigmaz(-sf%nzguard:sf%nz+sf%nzguard), sigmaz_next(-sf%nzguard:sf%nz+sf%nzguard), &
                  lsigmaz(-sf%nzguard:sf%nz+sf%nzguard), lsigmaz_next(-sf%nzguard:sf%nz+sf%nzguard)
    integer(ISZ) :: j,which
    real(kind=8) :: dt
    
    sf%afx = 1.
    sf%agx = 1.
    sf%afy = 1.
    sf%agy = 1.
    sf%afz = 1.
    sf%agz = 1.
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
!    return 
    if (sf%lsx/=0) then
      sigmax=0.
      sigmax_next=0.
      do j = 0, sf%nx
        sigmax(j)      = (sf%smaxx/sf%dx)*(REAL(j,8)/sf%sdeltax)**sf%nnx
        sigmax_next(j) = (sf%smaxx/sf%dx)*((REAL(j,8)+0.5)/sf%sdeltax)**sf%nnx
        lsigmax(sf%nx-j) = sigmax(j)
        lsigmax_next(sf%nx-j-1) = sigmax_next(j)
      end do
      select case(sf%lsx)
        case(1)
          do j = 0, sf%nx-1
            call assign_coefs(bnd_cond,sf%afx(j),sf%bpfx(j),sf%bmfx(j),sf%clight*dt,sf%dx,sigmax(j),sigmax_next(j),sb_coef,which)
            call assign_coefs(bnd_cond,sf%agx(j),sf%bpgx(j),sf%bmgx(j),sf%clight*dt,sf%dx,sigmax_next(j),sigmax(j+1),sb_coef,which)
          end do
        case(-1)
          do j = sf%nx, 1, -1
           call assign_coefs(bnd_cond,sf%afx(j),sf%bmfx(j),sf%bpfx(j),sf%clight*dt,sf%dx,lsigmax(j),lsigmax_next(j-1),sb_coef,which)
           call assign_coefs(bnd_cond,sf%agx(j-1),sf%bmgx(j-1),sf%bpgx(j-1),sf%clight*dt,sf%dx,lsigmax_next(j-1),lsigmax(j-1), &
                             sb_coef,which)
          end do
          sf%bmfx=-sf%bmfx
          sf%bpfx=-sf%bpfx
          sf%bmgx=-sf%bmgx
          sf%bpgx=-sf%bpgx
       end select
    end if

    if (sf%lsy/=0) then
      sigmay=0.
      sigmay_next=0.
      do j = 0, sf%ny
        sigmay(j)      = (sf%smaxy/sf%dy)*(REAL(j,8)/sf%sdeltay)**sf%nny
        sigmay_next(j) = (sf%smaxy/sf%dy)*((REAL(j,8)+0.5)/sf%sdeltay)**sf%nny
        lsigmay(sf%ny-j) = sigmay(j)
        lsigmay_next(sf%ny-j-1) = sigmay_next(j)
      end do
      select case(sf%lsy)
        case(1)
          do j = 0, sf%ny-1
            call assign_coefs(bnd_cond,sf%afy(j),sf%bpfy(j),sf%bmfy(j),sf%clight*dt,sf%dy,sigmay(j),sigmay_next(j),sb_coef,which)
            call assign_coefs(bnd_cond,sf%agy(j),sf%bpgy(j),sf%bmgy(j),sf%clight*dt,sf%dy,sigmay_next(j),sigmay(j+1),sb_coef,which)
          end do
        case(-1)
          do j = sf%ny, 1, -1
           call assign_coefs(bnd_cond,sf%afy(j),sf%bmfy(j),sf%bpfy(j),sf%clight*dt,sf%dy,lsigmay(j),lsigmay_next(j-1),sb_coef,which)
           call assign_coefs(bnd_cond,sf%agy(j-1),sf%bmgy(j-1),sf%bpgy(j-1),sf%clight*dt,sf%dy,lsigmay_next(j-1),lsigmay(j-1), &
                             sb_coef,which)
          end do
          sf%bmfy=-sf%bmfy
          sf%bpfy=-sf%bpfy
          sf%bmgy=-sf%bmgy
          sf%bpgy=-sf%bpgy
       end select
    end if

    if (sf%lsz/=0) then
      sigmaz=0.
      sigmaz_next=0.
      do j = 0, sf%nz
        sigmaz(j)      = (sf%smaxz/sf%dz)*(REAL(j,8)/sf%sdeltaz)**sf%nnz
        sigmaz_next(j) = (sf%smaxz/sf%dz)*((REAL(j,8)+0.5)/sf%sdeltaz)**sf%nnz
        lsigmaz(sf%nz-j) = sigmaz(j)
        lsigmaz_next(sf%nz-j-1) = sigmaz_next(j)
      end do
      select case(sf%lsz)
        case(1)
          do j = 0, sf%nz-1
            call assign_coefs(bnd_cond,sf%afz(j),sf%bpfz(j),sf%bmfz(j),sf%clight*dt,sf%dz,sigmaz(j),sigmaz_next(j),sb_coef,which)
            call assign_coefs(bnd_cond,sf%agz(j),sf%bpgz(j),sf%bmgz(j),sf%clight*dt,sf%dz,sigmaz_next(j),sigmaz(j+1),sb_coef,which)
          end do
        case(-1)
          do j = sf%nz, 1, -1
           call assign_coefs(bnd_cond,sf%afz(j),sf%bmfz(j),sf%bpfz(j),sf%clight*dt,sf%dz,lsigmaz(j),lsigmaz_next(j-1),sb_coef,which)
           call assign_coefs(bnd_cond,sf%agz(j-1),sf%bmgz(j-1),sf%bpgz(j-1),sf%clight*dt,sf%dz,lsigmaz_next(j-1),lsigmaz(j-1), &
                             sb_coef,which)
          end do
          sf%bmfz=-sf%bmfz
          sf%bpfz=-sf%bpfz
          sf%bmgz=-sf%bmgz
          sf%bpgz=-sf%bpgz
       end select
    end if

    return 
  end subroutine set_bndcoeffsem3d
    
end module mod_emfield3d

  subroutine init_splitfield(sf, nx, ny, nz, nxguard, nyguard, nzguard, dt, dx, dy, dz, clight, lsx, lsy, lsz, &
                             nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
    use mod_emfield3d
    TYPE(EM3D_SPLITYEEFIELDtype) :: sf
    INTEGER(ISZ), INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard, nnx, nny, nnz, lsx, lsy, lsz
    REAL(kind=8), INTENT(IN) :: dt, dx, dy, dz, clight, smaxx, smaxy, smaxz, sdeltax, sdeltay, sdeltaz
    integer(ISZ) :: j
    
    sf%nx = nx
    sf%ny = ny
    sf%nz = nz
    sf%nxguard = nxguard
    sf%nyguard = nyguard
    sf%nzguard = nzguard
    sf%dx = dx
    sf%dy = dy
    sf%dz = dz
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
    call EM3D_SPLITYEEFIELDtypeallot(sf)

    return 
  end subroutine init_splitfield
    

!************* SUBROUTINE init_fields  *************************************************
subroutine init_3dem_block(b, nx, ny, nz, nbndx, nbndy, nbndz, nxguard, nyguard, nzguard, dt, dx, dy, dz, clight, &
                           mu0, xmin, ymin, zmin, xlb, ylb, zlb, xrb, yrb, zrb)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
INTEGER(ISZ), INTENT(IN) :: nx, ny, nz, nxguard, nyguard, nzguard
INTEGER(ISZ), INTENT(IN) :: nbndx, nbndy, nbndz, xlb, ylb, zlb, xrb, yrb, zrb
REAL(kind=8), INTENT(IN) :: dt, dx, dy, dz, xmin, ymin, zmin, clight, mu0
INTEGER :: k,m, nnx, nny, nnz
real(kind=8) :: dtsdx, dtsdy, dtsdz, mudt, smaxx, sdeltax, smaxy, sdeltay, smaxz, sdeltaz

  b%nx = nx
  b%ny = ny
  b%nz = nz
  b%nbndx = nbndx
  b%nbndy = nbndy
  b%nbndz = nbndz
  b%xmin = xmin
  b%ymin = ymin
  b%zmin = zmin
  b%dx = dx
  b%dy = dy
  b%dz = dz
  b%xmax = xmin+dx*nx
  b%ymax = ymin+dy*ny
  b%zmax = zmin+dz*nz
  b%dxi = 1./dx
  b%dyi = 1./dy
  b%dzi = 1./dz
  b%xlbnd = xlb
  b%xrbnd = xrb
  b%ylbnd = ylb
  b%yrbnd = yrb
  b%zlbnd = zlb
  b%zrbnd = zrb

  allocate(b%core)
  allocate(b%core%yf)
  
  b%core%yf%nx = nx
  b%core%yf%ny = ny
  b%core%yf%nz = nz
  b%core%yf%nxguard = nxguard
  b%core%yf%nyguard = nyguard
  b%core%yf%nzguard = nzguard
  b%core%yf%xmin = xmin
  b%core%yf%ymin = ymin
  b%core%yf%zmin = zmin
  b%core%yf%dx = dx
  b%core%yf%dy = dy
  b%core%yf%dz = dz
  b%core%yf%xmax = xmin+dx*nx
  b%core%yf%ymax = ymin+dy*ny
  b%core%yf%zmax = zmin+dz*nz
  b%core%yf%dxi = 1./dx
  b%core%yf%dyi = 1./dy
  b%core%yf%dzi = 1./dz
  b%core%yf%clight = clight
  b%core%yf%mu0    = mu0

  call EM3D_YEEFIELDtypeallot(b%core%yf)
  
  nnx=nn
  nny=nn
  nnz=nn
  smaxx=s_max_init
  smaxy=s_max_init
  smaxz=s_max_init
  sdeltax=s_delta
  sdeltay=s_delta
  sdeltaz=s_delta

  allocate(b%sidexl,b%sidexr,b%sideyl,b%sideyr,b%sidezl,b%sidezr)
  allocate(b%edgexlyl,b%edgexryl,b%edgexlyr,b%edgexryr)
  allocate(b%edgeylzl,b%edgeyrzl,b%edgeylzr,b%edgeyrzr)
  allocate(b%edgexlzl,b%edgexrzl,b%edgexlzr,b%edgexrzr)
  allocate(b%cornerxlylzl,b%cornerxrylzl,b%cornerxlyrzl,b%cornerxlylzr,b%cornerxryrzl,b%cornerxrylzr,b%cornerxlyrzr,b%cornerxryrzr)

  allocate(b%sidexl%syf,b%sidexr%syf,b%sideyl%syf,b%sideyr%syf,b%sidezl%syf,b%sidezr%syf)
  allocate(b%edgexlyl%syf,b%edgexryl%syf,b%edgexlyr%syf,b%edgexryr%syf)
  allocate(b%edgeylzl%syf,b%edgeyrzl%syf,b%edgeylzr%syf,b%edgeyrzr%syf)
  allocate(b%edgexlzl%syf,b%edgexrzl%syf,b%edgexlzr%syf,b%edgexrzr%syf)
  allocate(b%cornerxlylzl%syf,b%cornerxrylzl%syf,b%cornerxlyrzl%syf,b%cornerxlylzr%syf, &
           b%cornerxryrzl%syf,b%cornerxrylzr%syf,b%cornerxlyrzr%syf,b%cornerxryrzr%syf)
  
! *** sides
! x
  call init_splitfield(b%sidexl%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%sidexr%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! y
  call init_splitfield(b%sideyl%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%sideyr%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! z
  call init_splitfield(b%sidezl%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%sidezr%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)

! *** edges
! xy
  call init_splitfield(b%edgexlyl%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexryl%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexlyr%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexryr%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 0, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! xz
  call init_splitfield(b%edgexlzl%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexrzl%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexlzr%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexrzr%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! yz
  call init_splitfield(b%edgeylzl%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgeyrzl%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgeylzr%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgeyrzr%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)

! *** corners
  call init_splitfield(b%cornerxlylzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxrylzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxlyrzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxryrzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1,-1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxlylzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxrylzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxlyrzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxryrzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 1, &
                       nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)


return

END subroutine init_3dem_block

subroutine push_em3d_e(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

if(f%nconds>0) then 
  call push_em3d_conde(f,dt)
  return
end if

dtsdx = f%clight**2*dt/f%dx
dtsdy = f%clight**2*dt/f%dy
dtsdz = f%clight**2*dt/f%dz
mudt  = f%mu0*f%clight**2*dt

  call push_em3d_evec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz,f%J, &
                      mudt,dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard)

return
end subroutine push_em3d_e

subroutine push_em3d_evec(ex,ey,ez,bx,by,bz,CJ,mudt,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard,3) :: CJ
real(kind=8), intent(IN) :: mudt,dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l

  ! advance Ex
  do l = 0, nz
   do k = 0, ny
    do j = 0, nx-1
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
      Ez(j,k,l) = Ez(j,k,l) + dtsdx * (By(j,k,l) - By(j-1,k  ,l)) &
                                - dtsdy * (Bx(j,k,l) - Bx(j  ,k-1,l)) &
                                - mudt  * CJ(j,k,l,3)
    end do
   end do
  end do

return
end subroutine push_em3d_evec

subroutine push_em3d_b(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz

dtsdx = dt/f%dx
dtsdy = dt/f%dy
dtsdz = dt/f%dz

  call push_em3d_bvec(f%ex,f%ey,f%ez,f%bx,f%by,f%bz, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard)

return
end subroutine push_em3d_b

subroutine push_em3d_bvec(ex,ey,ez,bx,by,bz,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,bx,by,bz
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l

  ! advance Bx
  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      Bx(j,k,l) = Bx(j,k,l) - dtsdy * (Ez(j,k+1,l  ) - Ez(j,k,l)) &
                                + dtsdz * (Ey(j,k,  l+1) - Ey(j,k,l))
    end do
   end do
  end do

  ! advance By
  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      By(j,k,l) = By(j,k,l) + dtsdx * (Ez(j+1,k,l  ) - Ez(j,k,l)) &  
                                - dtsdz * (Ex(j  ,k,l+1) - Ex(j,k,l)) 
    end do
   end do
  end do

  ! advance Bz 
  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      Bz(j,k,l) = Bz(j,k,l) - dtsdx * (Ey(j+1,k,l) - Ey(j,k,l)) &
                                + dtsdy * (Ex(j,k+1,l) - Ex(j,k,l))
    end do
   end do
  end do

return
end subroutine push_em3d_bvec

subroutine push_em3d_f(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,dtsepsi

if(f%nconds>0) then 
  call push_em3d_condf(f,dt)
  return
end if

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz
dtsepsi = f%mu0*f%clight**3*dt

  call push_em3d_fvec(f%ex,f%ey,f%ez,f%f, f%rho, &
                      dtsepsi,dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard)
  
end subroutine push_em3d_f

subroutine push_em3d_fvec(ex,ey,ez,f,rho,dtsepsi,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f,rho
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz,dtsepsi
integer(ISZ) :: j,k,l

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

return
end subroutine push_em3d_fvec

subroutine push_em3d_ef(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt

if(f%nconds>0) then 
  call push_em3d_condef(f,dt)
  return
end if

  dtsdx = f%clight*dt/f%dx
  dtsdy = f%clight*dt/f%dy
  dtsdz = f%clight*dt/f%dz

  call push_em3d_efvec(f%ex,f%ey,f%ez,f%f, &
                      dtsdx,dtsdy,dtsdz, &
                      f%nx,f%ny,f%nz, &
                      f%nxguard,f%nyguard,f%nzguard)

return
end subroutine push_em3d_ef

subroutine push_em3d_efvec(ex,ey,ez,f,dtsdx,dtsdy,dtsdz,nx,ny,nz,nxguard,nyguard,nzguard)
integer :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), intent(IN OUT), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard) :: ex,ey,ez,f
real(kind=8), intent(IN) :: dtsdx,dtsdy,dtsdz
integer(ISZ) :: j,k,l

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

return
end subroutine push_em3d_efvec

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

subroutine push_em3d_conde(f,dt)
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


  ! advance Ex
  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx-1
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j+1,k,l)) &
      f%Ex(j,k,l) = f%Ex(j,k,l) + dtsdy * (f%Bz(j,k,l)   - f%Bz(j,k-1,l  )) &
                                - dtsdz * (f%By(j,k,l)   - f%By(j,k  ,l-1)) &
                                - mudt  * f%J(j,k,l,1)
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k+1,l)) &
      f%Ey(j,k,l) = f%Ey(j,k,l) - dtsdx * (f%Bz(j,k,l)   - f%Bz(j-1,k,l)) &
                                + dtsdz * (f%Bx(j,k,l)   - f%Bx(j,k,l-1)) &
                                - mudt  * f%J(j,k,l,2)
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k,l+1)) &
      f%Ez(j,k,l) = f%Ez(j,k,l) + dtsdx * (f%By(j,k,l) - f%By(j-1,k  ,l)) &
                                - dtsdy * (f%Bx(j,k,l) - f%Bx(j  ,k-1,l)) &
                                - mudt  * f%J(j,k,l,3)
    end do
   end do
  end do

return
end subroutine push_em3d_conde

subroutine push_em3d_condef(f,dt)
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
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j+1,k,l)) &
      f%Ex(j,k,l) = f%Ex(j,k,l) + dtsdx * (f%F(j+1,k,l) - f%F(j,k,l)) 
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k+1,l)) &
      f%Ey(j,k,l) = f%Ey(j,k,l) + dtsdy * (f%F(j,k+1,l) - f%F(j,k,l))
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      if (.not.f%incond(j,k,l) .and. .not.f%incond(j,k,l+1)) &
      f%Ez(j,k,l) = f%Ez(j,k,l) + dtsdz * (f%F(j,k,l+1) - f%F(j,k,l)) 
    end do
   end do
  end do

return
end subroutine push_em3d_condef

subroutine push_em3d_condf(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,dtsepsi

dtsdx = f%clight*dt/f%dx
dtsdy = f%clight*dt/f%dy
dtsdz = f%clight*dt/f%dz
dtsepsi = f%mu0*f%clight**3*dt

  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx
      if (.not.f%incond(j,k,l)) &
      f%F(j,k,l) = f%F(j,k,l) + dtsdx * (f%Ex(j,k,l) - f%Ex(j-1,k  ,l  )) &
                              + dtsdy * (f%Ey(j,k,l) - f%Ey(j  ,k-1,l  )) &
                              + dtsdz * (f%Ez(j,k,l) - f%Ez(j  ,k  ,l-1)) &
                              - dtsepsi * f%Rho(j,k,l)
    end do
   end do
  end do

end subroutine push_em3d_condf

subroutine push_em3d_kyee_e(f,dt)
use mod_emfield3d
implicit none

TYPE(EM3D_KYEEFIELDtype) :: f
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l
real(kind=8) :: dtsdx,dtsdy,dtsdz,mudt, newj

dtsdx = f%clight**2*dt/f%dx
dtsdy = f%clight**2*dt/f%dy
dtsdz = f%clight**2*dt/f%dz
mudt  = 0.5*f%mu0*f%clight**2*dt


  ! advance Ex
  do l = 0, f%nz
   do k = 0, f%ny
    do j = 0, f%nx-1
      newj =  - (f%alphay+f%alphay)*mudt * f%J(j,k,l,1) &
                                - f%betay *mudt * (f%J(j+1,k,l  ,1) &
                                                +f%J(j-1,k,l  ,1) &
                                                +f%J(j  ,k,l+1,1) &
                                                +f%J(j  ,k,l-1,1)) &
                                - f%gammay*mudt * (f%J(j+1,k,l+1,1) &
                                                +f%J(j-1,k,l+1,1) &
                                                +f%J(j+1,k,l-1,1) &
                                                +f%J(j-1,k,l-1,1)) &
                                - f%betaz *mudt * (f%J(j+1,k  ,l,1) &
                                                +f%J(j-1,k  ,l,1) &
                                                +f%J(j  ,k+1,l,1) &
                                                +f%J(j  ,k-1,l,1)) &
                                - f%gammaz*mudt * (f%J(j+1,k+1,l,1) &
                                                +f%J(j-1,k+1,l,1) &
                                                +f%J(j+1,k-1,l,1) &
                                                +f%J(j-1,k-1,l,1))
      f%Ex(j,k,l) = f%Ex(j,k,l) + f%alphay*dtsdy * (f%Bz(j  ,k,l  ) - f%Bz(j  ,k-1,l  )) &
                                + f%betay *dtsdy * (f%Bz(j+1,k,l  ) - f%Bz(j+1,k-1,l  ) &
                                               +  f%Bz(j-1,k,l  ) - f%Bz(j-1,k-1,l  ) &
                                               +  f%Bz(j  ,k,l+1) - f%Bz(j  ,k-1,l+1) &
                                               +  f%Bz(j  ,k,l-1) - f%Bz(j  ,k-1,l-1))&
                                + f%gammay*dtsdy * (f%Bz(j+1,k,l+1) - f%Bz(j+1,k-1,l+1) &
                                               +  f%Bz(j-1,k,l+1) - f%Bz(j-1,k-1,l+1) &
                                               +  f%Bz(j+1,k,l-1) - f%Bz(j+1,k-1,l-1) &
                                               +  f%Bz(j-1,k,l-1) - f%Bz(j-1,k-1,l-1)) &
                                - f%alphaz*dtsdz * (f%By(j  ,k  ,l) - f%By(j  ,k  ,l-1)) &
                                - f%betaz *dtsdz * (f%By(j+1,k  ,l) - f%By(j+1,k  ,l-1)  &
                                               +  f%By(j-1,k  ,l) - f%By(j-1,k  ,l-1)  &
                                               +  f%By(j  ,k+1,l) - f%By(j  ,k+1,l-1)  &
                                               +  f%By(j  ,k-1,l) - f%By(j  ,k-1,l-1)) &
                                - f%gammaz*dtsdz * (f%By(j+1,k+1,l) - f%By(j+1,k+1,l-1)  &
                                               +  f%By(j-1,k+1,l) - f%By(j-1,k+1,l-1)  &
                                               +  f%By(j+1,k-1,l) - f%By(j+1,k-1,l-1)  &
                                               +  f%By(j-1,k-1,l) - f%By(j-1,k-1,l-1)) &
                                + 0.5*(newj+f%Jold(j,k,l,1))
      f%Jold(j,k,l,1) = newj
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      newj = - (f%alphax+f%alphaz)*mudt  * f%J(j,k,l,2) &
                                - f%betax *mudt * (f%J(j,k+1,l  ,2) &
                                              +  f%J(j,k-1,l  ,2) &
                                              +  f%J(j,k  ,l+1,2) &
                                              +  f%J(j,k  ,l-1,2)) &
                                - f%gammax*mudt * (f%J(j,k+1,l+1,2) &
                                              +  f%J(j,k-1,l+1,2) &
                                              +  f%J(j,k+1,l-1,2) &
                                              +  f%J(j,k-1,l-1,2)) &
                                - f%betaz *mudt * (f%J(j  ,k+1,l,2) &
                                              +  f%J(j  ,k-1,l,2) &
                                              +  f%J(j+1,k  ,l,2) &
                                              +  f%J(j-1,k  ,l,2)) &
                                - f%gammaz*mudt * (f%J(j+1,k+1,l,2) &
                                              +  f%J(j+1,k-1,l,2) &
                                              +  f%J(j-1,k+1,l,2) &
                                              +  f%J(j-1,k-1,l,2))
      f%Ey(j,k,l) = f%Ey(j,k,l) - f%alphax*dtsdx * (f%Bz(j,k  ,l  ) - f%Bz(j-1,k  ,l  )) &
                                - f%betax *dtsdx * (f%Bz(j,k+1,l  ) - f%Bz(j-1,k+1,l  ) &
                                               +  f%Bz(j,k-1,l  ) - f%Bz(j-1,k-1,l  ) &
                                               +  f%Bz(j,k  ,l+1) - f%Bz(j-1,k  ,l+1) &
                                               +  f%Bz(j,k  ,l-1) - f%Bz(j-1,k  ,l-1)) &
                                - f%gammax*dtsdx * (f%Bz(j,k+1,l+1) - f%Bz(j-1,k+1,l+1) &
                                               +  f%Bz(j,k-1,l+1) - f%Bz(j-1,k-1,l+1) &
                                               +  f%Bz(j,k+1,l-1) - f%Bz(j-1,k+1,l-1) &
                                               +  f%Bz(j,k-1,l-1) - f%Bz(j-1,k-1,l-1)) &                              
                                + f%alphaz*dtsdz * (f%Bx(j+1,k  ,l) - f%Bx(j+1,k  ,l-1) &
                                               +  f%Bx(j-1,k  ,l) - f%Bx(j-1,k  ,l-1) &
                                               +  f%Bx(j  ,k+1,l) - f%Bx(j  ,k+1,l-1) &
                                               +  f%Bx(j  ,k-1,l) - f%Bx(j  ,k-1,l-1)) &
                                + f%betaz *dtsdz * (f%Bx(j+1,k  ,l) - f%Bx(j+1,k  ,l-1) &
                                               +  f%Bx(j-1,k  ,l) - f%Bx(j-1,k  ,l-1) &
                                               +  f%Bx(j  ,k+1,l) - f%Bx(j  ,k+1,l-1) &
                                               +  f%Bx(j  ,k-1,l) - f%Bx(j  ,k-1,l-1)) &
                                + f%gammaz*dtsdz * (f%Bx(j+1,k+1,l) - f%Bx(j+1,k+1,l-1) &
                                               +  f%Bx(j-1,k+1,l) - f%Bx(j-1,k+1,l-1) &
                                               +  f%Bx(j+1,k-1,l) - f%Bx(j+1,k-1,l-1) &
                                               +  f%Bx(j-1,k-1,l) - f%Bx(j-1,k-1,l-1)) &
                                + 0.5*(newj+f%Jold(j,k,l,2))
      f%Jold(j,k,l,2) = newj
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      newj = - (f%alphax+f%alphay)*mudt * f%J(j,k,l,3) &
                                - f%betax *mudt * (f%J(j,k+1,l  ,3) &
                                                +  f%J(j,k-1,l  ,3) &
                                                +  f%J(j,k  ,l+1,3) &
                                                +  f%J(j,k  ,l-1,3)) &
                                - f%gammax*mudt * (f%J(j,k+1,l+1,3) &
                                                +  f%J(j,k-1,l+1,3) &
                                                +  f%J(j,k+1,l-1,3) &
                                                +  f%J(j,k-1,l-1,3)) &
                                - f%betay *mudt * (f%J(j+1,k,l  ,3) &
                                                +  f%J(j-1,k,l  ,3) &
                                                +  f%J(j  ,k,l+1,3) &
                                                +  f%J(j  ,k,l-1,3)) &
                                - f%gammay*mudt * (f%J(j+1,k,l+1,3) &
                                                +  f%J(j-1,k,l+1,3) &
                                                +  f%J(j+1,k,l-1,3) &
                                                +  f%J(j-1,k,l-1,3))
      f%Ez(j,k,l) = f%Ez(j,k,l) + f%alphax*dtsdx * (f%By(j,k  ,l  ) - f%By(j-1,k  ,l  )) &
                                + f%betax *dtsdx * (f%By(j,k+1,l  ) - f%By(j-1,k+1,l  ) &
                                                 +  f%By(j,k-1,l  ) - f%By(j-1,k-1,l  ) &
                                                 +  f%By(j,k  ,l+1) - f%By(j-1,k  ,l+1) &
                                                 +  f%By(j,k  ,l-1) - f%By(j-1,k  ,l-1)) &
                                + f%gammax*dtsdx * (f%By(j,k+1,l+1) - f%By(j-1,k+1,l+1) &
                                                 +  f%By(j,k-1,l+1) - f%By(j-1,k-1,l+1) &
                                                 +  f%By(j,k+1,l-1) - f%By(j-1,k+1,l-1) &
                                                 +  f%By(j,k-1,l-1) - f%By(j-1,k-1,l-1)) &
                                - f%alphay*dtsdy * (f%Bx(j  ,k,l  ) - f%Bx(j  ,k-1,l  )) &
                                - f%betay *dtsdy * (f%Bx(j+1,k,l  ) - f%Bx(j+1,k-1,l  ) &
                                                 +  f%Bx(j-1,k,l  ) - f%Bx(j-1,k-1,l  ) &
                                                 +  f%Bx(j  ,k,l+1) - f%Bx(j  ,k-1,l+1) &
                                                 +  f%Bx(j  ,k,l-1) - f%Bx(j  ,k-1,l-1)) &
                                - f%gammay*dtsdy * (f%Bx(j+1,k,l+1) - f%Bx(j+1,k-1,l+1) &
                                                 +  f%Bx(j-1,k,l+1) - f%Bx(j-1,k-1,l+1) &
                                                 +  f%Bx(j+1,k,l-1) - f%Bx(j+1,k-1,l-1) &
                                                 +  f%Bx(j-1,k,l-1) - f%Bx(j-1,k-1,l-1)) &
                                + 0.5*(newj+f%Jold(j,k,l,3))
      f%Jold(j,k,l,3) = newj
    end do
   end do
  end do

return
end subroutine push_em3d_kyee_e

subroutine push_em3d_splite(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  call push_em3d_splitevec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                           sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                           sf%afx,sf%afy,sf%afz, &
                           sf%bpfx,sf%bpfy,sf%bpfz, &
                           sf%bmfx,sf%bmfy,sf%bmfz)

  return
end subroutine push_em3d_splite

subroutine push_em3d_splitevec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx,exy,exz, &
                                                                                                       eyx,eyy,eyz, &
                                                                                                       ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz

INTEGER :: j, k, l

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

  return
end subroutine push_em3d_splitevec

subroutine push_em3d_splitef(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)
  call push_em3d_splitefvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%exx,sf%eyy,sf%ezz, &
                           sf%fx,sf%fy,sf%fz, &
                           sf%agx,sf%agy,sf%agz, &
                           sf%bpgx,sf%bpgy,sf%bpgz, &
                           sf%bmgx,sf%bmgy,sf%bmgz)


  return
end subroutine push_em3d_splitef

subroutine push_em3d_splitefvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,eyy,ezz,fx,fy,fz, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: exx,eyy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l

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

  return
end subroutine push_em3d_splitefvec

subroutine push_em3d_splitb(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  call push_em3d_splitbvec(sf%nx,sf%ny,sf%nz,sf%nxguard,sf%nyguard,sf%nzguard, &
                           sf%exx,sf%exy,sf%exz,sf%eyx,sf%eyy,sf%eyz,sf%ezx,sf%ezy,sf%ezz, &
                           sf%bxy,sf%byx,sf%bzx,sf%bxz,sf%byz,sf%bzy, &
                           sf%agx,sf%agy,sf%agz, &
                           sf%bpgx,sf%bpgy,sf%bpgz, &
                           sf%bmgx,sf%bmgy,sf%bmgz)

  return
end subroutine push_em3d_splitb

subroutine push_em3d_splitbvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,bxy,byx,bzx,bxz,byz,bzy, &
                               agx,agy,agz,bpgx,bpgy,bpgz,bmgx,bmgy,bmgz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: agx,bpgx,bmgx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: agy,bpgy,bmgy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: agz,bpgz,bmgz
!real(kind=8), dimension(-nxguard:,-nyguard:,-nzguard:), intent(inout) :: bxy,byx,bzx,bxz,byz,bzy
!real(kind=8), dimension(-nxguard:,-nyguard:,-nzguard:), intent(in) :: exx,exy,exz, &
!                                                                                                    eyx,eyy,eyz, &
!                                                                                                    ezx,ezy,ezz
!real(kind=8), dimension(-nxguard:), intent(in) :: agx,bpgx,bmgx
!real(kind=8), dimension(-nyguard:), intent(in) :: agy,bpgy,bmgy
!real(kind=8), dimension(-nzguard:), intent(in) :: agz,bpgz,bmgz

INTEGER :: j, k, l

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxy(j,k,l) = agy(k)*bxy(j,k,l) - bpgy(k)*(ezx(j,k+1,l  )+ezy(j,k+1,l  )+ezz(j,k+1,l  )) &
                                     - bmgy(k)*(ezx(j,k  ,l  )+ezy(j,k  ,l  )+ezz(j,k  ,l  ))
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny-1
    do j = 0, nx
      bxz(j,k,l) = agz(l)*bxz(j,k,l) + bpgz(l)*(eyx(j,k  ,l+1)+eyy(j,k  ,l+1)+eyz(j,k  ,l+1)) &
                                     + bmgz(l)*(eyx(j,k  ,l  )+eyy(j,k  ,l  )+eyz(j,k  ,l  ))
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      byx(j,k,l) = agx(j)*byx(j,k,l) + bpgx(j)*(ezx(j+1,k,l  )+ezy(j+1,k,l  )+ezz(j+1,k,l  )) &
                                     + bmgx(j)*(ezx(j  ,k,l  )+ezy(j  ,k,l  )+ezz(j  ,k,l  ))
    end do
   end do
  end do

  do l = 0, nz-1
   do k = 0, ny
    do j = 0, nx-1
      byz(j,k,l) = agz(l)*byz(j,k,l) - bpgz(l)*(exx(j  ,k,l+1)+exy(j  ,k,l+1)+exz(j  ,k,l+1)) &
                                     - bmgz(l)*(exx(j  ,k,l  )+exy(j  ,k,l  )+exz(j  ,k,l  ))
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      bzx(j,k,l) = agx(j)*bzx(j,k,l) - bpgx(j)*(eyx(j+1,k  ,l)+eyy(j+1,k  ,l)+eyz(j+1,k  ,l)) &
                                     - bmgx(j)*(eyx(j  ,k  ,l)+eyy(j  ,k  ,l)+eyz(j  ,k  ,l))
    end do
   end do
  end do

  do l = 0, nz
   do k = 0, ny-1
    do j = 0, nx-1
      bzy(j,k,l) = agy(k)*bzy(j,k,l) + bpgy(k)*(exx(j  ,k+1,l)+exy(j  ,k+1,l)+exz(j  ,k+1,l)) &
                                     + bmgy(k)*(exx(j  ,k  ,l)+exy(j  ,k  ,l)+exz(j  ,k  ,l))
    end do
   end do
  end do

  return
end subroutine push_em3d_splitbvec

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
                           sf%bmfx,sf%bmfy,sf%bmfz)

  return
end subroutine push_em3d_splitf

subroutine push_em3d_splitfvectest(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,fx,fy,fz, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz

!real(kind=8), dimension(-nxguard:nx+nxguard) :: exm,exp

INTEGER :: j, k, l

!  ex = exx+exy+exz

  do l = 0, nz
   do k = 0, ny
    do j = 0, nx
!     fx(j,k,l) = afx(j)*fx(j,k,l) + bpfx(j)*(exx(j  ,k  ,l  )+exy(j  ,k  ,l  )+exz(j  ,k  ,l  )) &
!                                  + bmfx(j)*(exx(j-1,k  ,l  )+exy(j-1,k  ,l  )+exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
     fy(j,k,l) = fy(j,k,l) + bpfy(1)*(eyx(j  ,k  ,l  )) &
                                  + bmfy(1)*(eyx(j  ,k-1,l  )) !- (1._8/3._8)*dt*rho(j,k,l)
    end do
   end do
  end do

  return
end subroutine push_em3d_splitfvectest

subroutine push_em3d_splitfvec(nx,ny,nz,nxguard,nyguard,nzguard, &
                               exx,exy,exz,eyx,eyy,eyz,ezx,ezy,ezz,fx,fy,fz, &
                               afx,afy,afz,bpfx,bpfy,bpfz,bmfx,bmfy,bmfz)
implicit none

integer(ISZ), INTENT(IN) :: nx,ny,nz,nxguard,nyguard,nzguard
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(inout) :: fx,fy,fz
real(kind=8), dimension(-nxguard:nx+nxguard,-nyguard:ny+nyguard,-nzguard:nz+nzguard), intent(in) :: exx,exy,exz, &
                                                                                                    eyx,eyy,eyz, &
                                                                                                    ezx,ezy,ezz
real(kind=8), dimension(-nxguard:nx+nxguard), intent(in) :: afx,bpfx,bmfx
real(kind=8), dimension(-nyguard:ny+nyguard), intent(in) :: afy,bpfy,bmfy
real(kind=8), dimension(-nzguard:nz+nzguard), intent(in) :: afz,bpfz,bmfz

INTEGER :: j, k, l

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

  return
end subroutine push_em3d_splitfvec

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
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l

  call push_em3d_f(b%core%yf,dt*0.5)
  call push_em3d_b(b%core%yf,dt*0.5)
  call push_em3d_blockbndf(b,dt,2)
  call push_em3d_blockbndb(b,dt,2)

  call em3d_exchange_b(b)

  call push_em3d_ef(b%core%yf,dt)
  call push_em3d_e(b%core%yf,dt)
  call push_em3d_blockbndef(b,dt,0)
  call push_em3d_blockbnde(b,dt,0)
  
  call em3d_exchange_e(b)

  call push_em3d_f(b%core%yf,dt*0.5)
  call push_em3d_b(b%core%yf,dt*0.5)
  call push_em3d_blockbndf(b,dt,1)
  call push_em3d_blockbndb(b,dt,1)

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
  call em3d_exchange_e(b)

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
    call push_em3d_b(b%core%yf,dt)
  else
    if(l_pushf) call push_em3d_f(b%core%yf,dt*0.5)
    call push_em3d_b(b%core%yf,dt*0.5)
  endif

  if(l_pushf) call push_em3d_blockbndf(b,dt,which)
  call push_em3d_blockbndb(b,dt,which)

  call em3d_applybc_b(b%core%yf, &
                      b%xlbnd, &
                      b%xrbnd, &
                      b%ylbnd, &
                      b%yrbnd, &
                      b%zlbnd, &
                      b%zrbnd)
  call em3d_exchange_b(b)

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

subroutine em3d_applybc_e(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%ex(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = -f%ex(f%ixmin:f%ixmin+f%nxguard-1,:,:)
     f%ey(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = f%ey(f%ixmin+1:f%ixmin+f%nxguard,:,:)
     f%ez(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = f%ez(f%ixmin+1:f%ixmin+f%nxguard,:,:)
  end if
  if (xrbnd==neumann) then
     f%ex(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:) = -f%ex(f%ixmax-f%nxguard-1:f%ixmax,:,:)
     f%ey(f%ixmax+1:f%ixmax+f%nxguard:-1,:,:) = f%ey(f%ixmax-f%nxguard:f%ixmax-1,:,:)
     f%ez(f%ixmax+1:f%ixmax+f%nxguard:-1,:,:) = f%ez(f%ixmax-f%nxguard:f%ixmax-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%ey(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = -f%ey(:,f%iymin:f%iymin+f%nyguard-1,:)
     f%ex(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = f%ex(:,f%iymin+1:f%iymin+f%nyguard,:)
     f%ez(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = f%ez(:,f%iymin+1:f%iymin+f%nyguard,:)
  end if
  if (yrbnd==neumann) then
     f%ey(:,f%iymax:f%iymax+f%nyguard-1:-1,:) = -f%ey(:,f%iymax-f%nyguard-1:f%iymax,:)
     f%ex(:,f%iymax+1:f%iymax+f%nyguard:-1,:) = f%ex(:,f%iymax-f%nyguard:f%iymax-1,:)
     f%ez(:,f%iymax+1:f%iymax+f%nyguard:-1,:) = f%ez(:,f%iymax-f%nyguard:f%iymax-1,:)
  end if

  return
end subroutine em3d_applybc_e

subroutine em3d_applybc_b(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%bx(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = f%bx(f%ixmin+1:f%ixmin+f%nxguard,:,:)
     f%by(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = -f%by(f%ixmin:f%ixmin+f%nxguard-1,:,:)
     f%bz(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:) = -f%bz(f%ixmin:f%ixmin+f%nxguard-1,:,:)
  end if
  if (xrbnd==neumann) then
     f%bx(f%ixmax+1:f%ixmax+f%nxguard:-1,:,:) = f%bx(f%ixmax-f%nxguard:f%ixmax-1,:,:)
     f%by(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:) = -f%by(f%ixmax-f%nxguard-1:f%ixmax,:,:)
     f%bz(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:) = -f%bz(f%ixmax-f%nxguard-1:f%ixmax,:,:)
  end if

  if (ylbnd==neumann) then
     f%by(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = f%by(:,f%iymin+1:f%iymin+f%nyguard,:)
     f%bx(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = -f%bx(:,f%iymin:f%iymin+f%nyguard-1,:)
     f%bz(:,f%iymin-1:f%iymin-f%nyguard:-1,:) = -f%bz(:,f%iymin:f%iymin+f%nyguard-1,:)
  end if
  if (yrbnd==neumann) then
     f%by(:,f%iymax+1:f%iymax+f%nyguard:-1,:) = f%by(:,f%iymax-f%nyguard:f%iymax-1,:)
     f%bx(:,f%iymax:f%iymax+f%nyguard-1:-1,:) = -f%bx(:,f%iymax-f%nyguard-1:f%iymax,:)
     f%bz(:,f%iymax:f%iymax+f%nyguard-1:-1,:) = -f%bz(:,f%iymax-f%nyguard-1:f%iymax,:)
  end if

  return
end subroutine em3d_applybc_b

subroutine em3d_applybc_j(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) = f%j(f%ixmin:f%ixmin+f%nxguard,:,:,2:3) + f%j(f%ixmin:f%ixmin-f%nxguard:-1,:,:,2:3)
     f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) = f%j(f%ixmin:f%ixmin+f%nxguard-1,:,:,1) - f%j(f%ixmin-1:f%ixmin-f%nxguard:-1,:,:,1)
  end if
  if (xrbnd==neumann) then
     f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) = f%j(f%ixmax-f%nxguard:f%ixmax,:,:,2:3) + f%j(f%ixmax:f%ixmax+f%nxguard:-1,:,:,2:3)
     f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) = f%j(f%ixmax-f%nxguard:f%ixmax-1,:,:,1) - f%j(f%ixmax:f%ixmax+f%nxguard-1:-1,:,:,1)
  end if

  if (ylbnd==neumann) then
     f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) = f%j(:,f%iymin:f%iymin+f%nyguard,:,1:3:2) &
                                              + f%j(:,f%iymin:f%iymin-f%nyguard:-1,:,1:3:2)
     f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) = f%j(:,f%iymin:f%iymin+f%nyguard-1,:,2) - f%j(:,f%iymin-1:f%iymin-f%nyguard:-1,:,2)
  end if
  if (yrbnd==neumann) then
     f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) = f%j(:,f%iymax-f%nyguard:f%iymax,:,1:3:2) &
                                              + f%j(:,f%iymax:f%iymax+f%nyguard:-1,:,1:3:2)
     f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) = f%j(:,f%iymax-f%nyguard:f%iymax-1,:,2) - f%j(:,f%iymax:f%iymax+f%nyguard-1:-1,:,2)
  end if

  return
end subroutine em3d_applybc_j

subroutine em3d_applybc_rho(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_YEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

  if (xlbnd==neumann) then
     f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) = f%rho(f%ixmin:f%ixmin+f%nxguard,:,:) + f%rho(f%ixmin:f%ixmin-f%nxguard:-1,:,:)
  end if
  if (xrbnd==neumann) then
     f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) = f%rho(f%ixmax-f%nxguard:f%ixmax,:,:) + f%rho(f%ixmax:f%ixmax+f%nxguard:-1,:,:)
  end if

  if (ylbnd==neumann) then
     f%rho(:,f%iymin:f%iymin+f%nyguard,:) = f%rho(:,f%iymin:f%iymin+f%nyguard,:) + f%rho(:,f%iymin:f%iymin-f%nyguard:-1,:)
  end if
  if (yrbnd==neumann) then
     f%rho(:,f%iymax-f%nyguard:f%iymax,:) = f%rho(:,f%iymax-f%nyguard:f%iymax,:) + f%rho(:,f%iymax:f%iymax+f%nyguard:-1,:)
  end if

  return
end subroutine em3d_applybc_rho

subroutine em3d_applybc_splite(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

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

  return
end subroutine em3d_applybc_splite

subroutine em3d_applybc_splitb(f,xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: f
integer(ISZ) :: xlbnd,xrbnd,ylbnd,yrbnd,zlbnd,zrbnd

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

  return
end subroutine em3d_applybc_splitb

subroutine em3d_exchange_bnde_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif
  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
           yfl%ex(yfl%ixmax:yfl%ixmaxg-1,:,:) = yfu%ex(yfu%ixmin:yfu%ixmin+yfu%nxguard-1,:,:)
           yfu%ex(yfu%ixming:yfu%ixmin-1,:,:) = yfl%ex(yfl%ixmax-yfl%nxguard:yfl%ixmax-1,:,:)

           yfl%ey(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%ey(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
           yfu%ey(yfu%ixming+1:yfu%ixmin-1  ,:,:) = yfl%ey(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

           yfl%ez(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%ez(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
           yfu%ez(yfu%ixming+1:yfu%ixmin-1  ,:,:) = yfl%ez(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)

        case(splityeefield)
          syfu=>fu%syf
          yfl%ex(yfl%nx,:,:) = syfu%exx(0,:,:)+syfu%exy(0,:,:)+syfu%exz(0,:,:)
          syfu%exx(-1,:,:) = yfl%ex(yfl%nx-1,:,:)
          syfu%exy(-1,:,:) = 0.
          syfu%exz(-1,:,:) = 0.
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%ex(-1,:,:) = syfl%exx(syfl%nx-1,:,:)+syfl%exy(syfl%nx-1,:,:)+syfl%exz(syfl%nx-1,:,:)
          syfl%exx(syfl%nx,:,:) = yfu%ex(0,:,:)
          syfl%exy(syfl%nx,:,:) = 0.
          syfl%exz(syfl%nx,:,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
          syfu%exx(-1,:,:)      = syfl%exx(syfl%nx-1,:,:)
          syfu%exy(-1,:,:)      = syfl%exy(syfl%nx-1,:,:)
          syfu%exz(-1,:,:)      = syfl%exz(syfl%nx-1,:,:)
          syfl%exx(syfl%nx,:,:) = syfu%exx(0,:,:)
          syfl%exy(syfl%nx,:,:) = syfu%exy(0,:,:)
          syfl%exz(syfl%nx,:,:) = syfu%exz(0,:,:)
      end select
  end select

  return
end subroutine em3d_exchange_bnde_x

subroutine em3d_exchange_bnde_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
           yfl%ey(:,yfl%iymax:yfl%iymaxg-1,:) = yfu%ey(:,yfu%iymin:yfu%iymin+yfu%nyguard-1,:)
           yfu%ey(:,yfu%iyming:yfu%iymin-1,:) = yfl%ey(:,yfl%iymax-yfl%nyguard:yfl%iymax-1,:)

           yfl%ex(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%ex(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
           yfu%ex(:,yfu%iyming+1:yfu%iymin-1  ,:) = yfl%ex(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)

           yfl%ez(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%ez(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
           yfu%ez(:,yfu%iyming+1:yfu%iymin-1  ,:) = yfl%ez(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)

        case(splityeefield)
          syfu=>fu%syf
          yfl%ey(:,yfl%ny,:) = syfu%eyx(:,0,:)+syfu%eyy(:,0,:)+syfu%eyz(:,0,:)
          syfu%eyx(:,-1,:) = 0.
          syfu%eyy(:,-1,:) = yfl%ey(:,yfl%ny-1,:)
          syfu%eyz(:,-1,:) = 0.
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%ey(:,-1,:) = syfl%eyx(:,syfl%ny-1,:)+syfl%eyy(:,syfl%ny-1,:)+syfl%eyz(:,syfl%ny-1,:)
          syfl%eyx(:,syfl%ny,:) = 0.
          syfl%eyy(:,syfl%ny,:) = yfu%ey(:,0,:)
          syfl%eyz(:,syfl%ny,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
          syfu%eyx(:,-1,:)      = syfl%eyx(:,syfl%ny-1,:)
          syfu%eyy(:,-1,:)      = syfl%eyy(:,syfl%ny-1,:)
          syfu%eyz(:,-1,:)      = syfl%eyz(:,syfl%ny-1,:)
          syfl%eyx(:,syfl%ny,:) = syfu%eyx(:,0,:)
          syfl%eyy(:,syfl%ny,:) = syfu%eyy(:,0,:)
          syfl%eyz(:,syfl%ny,:) = syfu%eyz(:,0,:)
      end select
  end select

  return
end subroutine em3d_exchange_bnde_y

subroutine em3d_exchange_bnde_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

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
          if (fl%proc/=my_index) then
            ! --- send data down in z
            call mpi_packbuffer_init((3*yfu%nzguard-2)*size(yfu%ez(:,:,0)),1)
            do iz = yfu%izmin,yfu%izmin+yfu%nzguard-1
              call mympi_pack(yfu%ez(:,:,iz),1)
            end do
            if (yfu%nzguard>1) then
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                call mympi_pack(yfu%ex(:,:,iz),1)
              end do
              do iz = yfu%izmin+1,yfu%izmin+yfu%nzguard-1
                call mympi_pack(yfu%ey(:,:,iz),1)
              end do
            end if
            call mpi_isend_pack(fl%proc,1,1)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            call mpi_packbuffer_init((3*yfl%nzguard-2)*size(yfl%ez(:,:,0)),2)
            do iz =yfl%izmax-yfl%nzguard,yfl%izmax-1
              call mympi_pack(yfl%ez(:,:,iz),2)
            end do
            if (yfl%nzguard>1) then
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                call mympi_pack(yfl%ex(:,:,iz),2)
              end do
              do iz = yfl%izmax-yfl%nzguard+1,yfl%izmax-1
                call mympi_pack(yfl%ey(:,:,iz),2)
              end do
            end if
            call mpi_isend_pack(fu%proc,2,2)
          else
#endif
            yfl%ez(:,:,yfl%izmax:yfl%izmaxg-1) = yfu%ez(:,:,yfu%izmin:yfu%izmin+yfu%nzguard-1)
            yfu%ez(:,:,yfu%izming:yfu%izmin-1) = yfl%ez(:,:,yfl%izmax-yfl%nzguard:yfl%izmax-1)

            yfl%ex(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%ex(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
            yfu%ex(:,:,yfu%izming+1:yfu%izmin-1  ) = yfl%ex(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)

            yfl%ey(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%ey(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
            yfu%ey(:,:,yfu%izming+1:yfu%izmin-1  ) = yfl%ey(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)

#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          if (fl%proc/=fu%proc) return
          syfu=>fu%syf
          yfl%ez(:,:,yfl%nz) = syfu%ezx(:,:,0)+syfu%ezy(:,:,0)+syfu%ezz(:,:,0)
          syfu%ezx(:,:,-1)   = 0.
          syfu%ezy(:,:,-1)   = 0.
          syfu%ezz(:,:,-1)   = yfl%ez(:,:,yfl%nz-1)
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          if (fl%proc/=fu%proc) return
          yfu=>fu%yf
          yfu%ez(:,:,-1)        = syfl%ezx(:,:,syfl%nz-1)+syfl%ezy(:,:,syfl%nz-1)+syfl%ezz(:,:,syfl%nz-1)
          syfl%ezx(:,:,syfl%nz) = 0.
          syfl%ezy(:,:,syfl%nz) = 0.
          syfl%ezz(:,:,syfl%nz) = yfu%ez(:,:,0)
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            call MPI_ISEND(syfu%ezx(-syfu%nxguard,-syfu%nyguard,0),int(size(syfu%ez(:,:,:,0)),MPIISZ),MPI_DOUBLE_PRECISION,&
                                   int(fl%proc,MPIISZ),3,MPI_COMM_WORLD,mpirequests(mpireqpnt+1),mpierror)
            mpireqpnt=mpireqpnt+1
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            call MPI_ISEND(syfl%ezx(-syfl%nxguard,-syfl%nyguard,syfl%nz-1),int(size(syfl%ez(:,:,:,syfl%nz-1)),MPIISZ),&
                           MPI_DOUBLE_PRECISION,int(fu%proc,MPIISZ),4,MPI_COMM_WORLD,mpirequests(mpireqpnt+1),mpierror)
            mpireqpnt=mpireqpnt+1
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%ezx(:,:,-1)      = syfl%ezx(:,:,syfl%nz-1)
          syfu%ezy(:,:,-1)      = syfl%ezy(:,:,syfl%nz-1)
          syfu%ezz(:,:,-1)      = syfl%ezz(:,:,syfl%nz-1)
          syfl%ezx(:,:,syfl%nz) = syfu%ezx(:,:,0)
          syfl%ezy(:,:,syfl%nz) = syfu%ezy(:,:,0)
          syfl%ezz(:,:,syfl%nz) = syfu%ezz(:,:,0)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bnde_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bnde_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            call mpi_packbuffer_init((3*yfu%nzguard-2)*size(yfu%Ez(:,:,0)),4)
            call mpi_recv_pack(fl%proc,2,4)
            do iz = yfu%izming,yfu%izmin-1
              yfu%ez(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),4),shape(yfu%Ez(:,:,0)))
            end do
            if (yfu%nzguard>1) then
              do iz = yfu%izming+1,yfu%izmin-1
                yfu%ex(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),4),shape(yfu%Ez(:,:,0)))
              end do
              do iz = yfu%izming+1,yfu%izmin-1
                yfu%ey(:,:,iz) = reshape(mpi_unpack_real_array( size(yfu%Ez(:,:,0)),4),shape(yfu%Ez(:,:,0)))
              end do
            end if
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            call mpi_packbuffer_init((3*yfl%nzguard-2)*size(yfl%ez(:,:,yfl%izmin)),3)
            call mpi_recv_pack(fu%proc,1,3)
            do iz = yfl%izmax,yfl%izmaxg-1
              yfl%ez(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),3),shape(yfl%Ez(:,:,0)))
            end do
            if (yfl%nzguard>1) then
              do iz = yfl%izmax+1,yfl%izmaxg-1
                yfl%ex(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),3),shape(yfl%Ez(:,:,0)))
              end do
              do iz = yfl%izmax+1,yfl%izmaxg-1
                yfl%ey(:,:,iz) = reshape(mpi_unpack_real_array( size(yfl%Ez(:,:,0)),3),shape(yfl%Ez(:,:,0)))
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
            syfu%ez(:,:,:,-1) = reshape( mpi_recv_real(size(syfu%ez(:,:,:,-1)),fl%proc,4) ,shape(syfu%ez(:,:,:,-1)))
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            syfl%ez(:,:,:,syfl%nz) = reshape(mpi_recv_real(size(syfl%ez(:,:,:,syfl%nz)),fu%proc,3),shape(syfl%ez(:,:,:,syfl%nz)))
          end if
      end select
  end select
!  call parallelbarrier()
  return
end subroutine em3d_exchange_bnde_zrecv
#endif

subroutine em3d_exchange_bndb_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%by(yfl%ixmax:yfl%ixmaxg-1,:,:) = yfu%by(yfu%ixmin:yfu%ixmin+yfu%nxguard-1,:,:)
          yfl%bz(yfl%ixmax:yfl%ixmaxg-1,:,:) = yfu%bz(yfu%ixmin:yfu%ixmin+yfu%nxguard-1,:,:)
          yfu%by(yfu%ixming:yfu%ixmin-1,:,:) = yfl%by(yfl%ixmax-yfl%nxguard:yfl%ixmax-1,:,:)
          yfu%bz(yfu%ixming:yfu%ixmin-1,:,:) = yfl%bz(yfl%ixmax-yfl%nxguard:yfl%ixmax-1,:,:)
          yfl%bx(yfl%ixmax+1:yfl%ixmaxg-1,:,:) = yfu%by(yfu%ixmin+1:yfu%ixmin+yfu%nxguard-1,:,:)
          yfu%bx(yfu%ixming+1:yfu%ixmin-1,:,:) = yfl%by(yfl%ixmax-yfl%nxguard+1:yfl%ixmax-1,:,:)
        case(splityeefield)
          syfu=>fu%syf
          yfl%by(yfl%nx,:,:) = (syfu%byx(0,:,:)+syfu%byz(0,:,:))/syfu%clight
          yfl%bz(yfl%nx,:,:) = (syfu%bzx(0,:,:)+syfu%bzy(0,:,:))/syfu%clight
          syfu%byx(-1,:,:) = yfl%by(yfl%nx-1,:,:)*syfu%clight
          syfu%byz(-1,:,:) = 0.
          syfu%bzx(-1,:,:) = yfl%bz(yfl%nx-1,:,:)*syfu%clight
          syfu%bzy(-1,:,:) = 0.
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%by(-1,:,:) = (syfl%byx(syfl%nx-1,:,:)+syfl%byz(syfl%nx-1,:,:))/syfl%clight
          yfu%bz(-1,:,:) = (syfl%bzx(syfl%nx-1,:,:)+syfl%bzy(syfl%nx-1,:,:))/syfl%clight
          syfl%byx(syfl%nx,:,:) = yfu%by(0,:,:)*syfl%clight
          syfl%byz(syfl%nx,:,:) = 0.
          syfl%bzx(syfl%nx,:,:) = yfu%bz(0,:,:)*syfl%clight
          syfl%bzy(syfl%nx,:,:) = 0.
        case(splityeefield)
          syfu=>fu%syf
          syfu%byx(-1,:,:) = syfl%byx(syfl%nx-1,:,:)
          syfu%byz(-1,:,:) = syfl%byz(syfl%nx-1,:,:)
          syfu%bzx(-1,:,:) = syfl%bzx(syfl%nx-1,:,:)
          syfu%bzy(-1,:,:) = syfl%bzy(syfl%nx-1,:,:)
          syfl%byx(syfl%nx,:,:) = syfu%byx(0,:,:)
          syfl%byz(syfl%nx,:,:) = syfu%byz(0,:,:)
          syfl%bzx(syfl%nx,:,:) = syfu%bzx(0,:,:)
          syfl%bzy(syfl%nx,:,:) = syfu%bzy(0,:,:)
      end select
  end select

  return
end subroutine em3d_exchange_bndb_x

subroutine em3d_exchange_bndb_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%bx(:,yfl%iymax:yfl%iymaxg-1,:) = yfu%bx(:,yfu%iymin:yfu%iymin+yfu%nyguard-1,:)
          yfl%bz(:,yfl%iymax:yfl%iymaxg-1,:) = yfu%bz(:,yfu%iymin:yfu%iymin+yfu%nyguard-1,:)
          yfu%bx(:,yfu%iyming:yfu%iymin-1,:) = yfl%bx(:,yfl%iymax-yfl%nyguard:yfl%iymax-1,:)
          yfu%bz(:,yfu%iyming:yfu%iymin-1,:) = yfl%bz(:,yfl%iymax-yfl%nyguard:yfl%iymax-1,:)
          yfl%by(:,yfl%iymax+1:yfl%iymaxg-1,:) = yfu%bz(:,yfu%iymin+1:yfu%iymin+yfu%nyguard-1,:)
          yfu%by(:,yfu%iyming+1:yfu%iymin-1,:) = yfl%bz(:,yfl%iymax-yfl%nyguard+1:yfl%iymax-1,:)
        case(splityeefield)
          syfu=>fu%syf
          yfl%bx(:,yfl%ny,:) = (syfu%bxy(:,0,:)+syfu%bxz(:,0,:))/syfu%clight
          yfl%bz(:,yfl%ny,:) = (syfu%bzx(:,0,:)+syfu%bzy(:,0,:))/syfu%clight
          syfu%bxy(:,-1,:) = yfl%bx(:,yfl%ny-1,:)*syfu%clight
          syfu%bxz(:,-1,:) = 0.
          syfu%bzx(:,-1,:) = 0.
          syfu%bzy(:,-1,:) = yfl%bz(:,yfl%ny-1,:)*syfu%clight
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%bx(:,-1,:) = (syfl%bxy(:,syfl%ny-1,:)+syfl%bxz(:,syfl%ny-1,:))/syfl%clight
          yfu%bz(:,-1,:) = (syfl%bzx(:,syfl%ny-1,:)+syfl%bzy(:,syfl%ny-1,:))/syfl%clight
          syfl%bxy(:,syfl%ny,:) = yfu%bx(:,0,:)*syfl%clight
          syfl%bxz(:,syfl%ny,:) = 0.
          syfl%bzx(:,syfl%ny,:) = 0.
          syfl%bzy(:,syfl%ny,:) = yfu%bz(:,0,:)*syfl%clight
        case(splityeefield)
          syfu=>fu%syf
          syfu%bxy(:,-1,:) = syfl%bxy(:,syfl%ny-1,:)
          syfu%bxz(:,-1,:) = syfl%bxz(:,syfl%ny-1,:)
          syfu%bzx(:,-1,:) = syfl%bzx(:,syfl%ny-1,:)
          syfu%bzy(:,-1,:) = syfl%bzy(:,syfl%ny-1,:)
          syfl%bxy(:,syfl%ny,:) = syfu%bxy(:,0,:)
          syfl%bxz(:,syfl%ny,:) = syfu%bxz(:,0,:)
          syfl%bzx(:,syfl%ny,:) = syfu%bzx(:,0,:)
          syfl%bzy(:,syfl%ny,:) = syfu%bzy(:,0,:)
      end select
  end select

  return
end subroutine em3d_exchange_bndb_y

subroutine em3d_exchange_bndb_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2),mpierror
integer(ISZ) :: ibuf
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
            ibuf = 5
            call mpi_packbuffer_init((3*yfu%nzguard-1)*size(yfu%by(:,:,0)),ibuf)
            call mympi_pack(yfu%bx(:,:,yfu%izmin:yfu%izmin+yfu%nzguard-1),ibuf)
            call mympi_pack(yfu%by(:,:,yfu%izmin:yfu%izmin+yfu%nzguard-1),ibuf)
            if (yfu%nzguard>1) &
            call mympi_pack(yfu%bz(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1),ibuf)
            call mpi_isend_pack(fl%proc,1,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 6
            call mpi_packbuffer_init((3*yfl%nzguard-1)*size(yfl%by(:,:,0)),ibuf)
            call mympi_pack(yfl%bx(:,:,yfl%izmax-yfl%nzguard:yfl%izmax-1),ibuf)
            call mympi_pack(yfl%by(:,:,yfl%izmax-yfl%nzguard:yfl%izmax-1),ibuf)
            if (yfl%nzguard>1) &
            call mympi_pack(yfl%bz(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1),ibuf)
            call mpi_isend_pack(fu%proc,2,ibuf)
          else
#endif
          yfl%bx(:,:,yfl%izmax:yfl%izmaxg-1) = yfu%bx(:,:,yfu%izmin:yfu%izmin+yfu%nzguard-1)
          yfl%by(:,:,yfl%izmax:yfl%izmaxg-1) = yfu%by(:,:,yfu%izmin:yfu%izmin+yfu%nzguard-1)
          yfu%bx(:,:,yfu%izming:yfu%izmin-1) = yfl%bx(:,:,yfl%izmax-yfl%nzguard:yfl%izmax-1)
          yfu%by(:,:,yfu%izming:yfu%izmin-1) = yfl%by(:,:,yfl%izmax-yfl%nzguard:yfl%izmax-1)
          yfl%bz(:,:,yfl%izmax+1:yfl%izmaxg-1) = yfu%by(:,:,yfu%izmin+1:yfu%izmin+yfu%nzguard-1)
          yfu%bz(:,:,yfu%izming+1:yfu%izmin-1) = yfl%by(:,:,yfl%izmax-yfl%nzguard+1:yfl%izmax-1)
#ifdef MPIPARALLEL
          end if
#endif
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=fu%proc) return
          yfl%bx(:,:,yfl%nz) = (syfu%bxy(:,:,0)+syfu%bxz(:,:,0))/syfu%clight
          yfl%by(:,:,yfl%nz) = (syfu%byx(:,:,0)+syfu%byz(:,:,0))/syfu%clight
          syfu%bxy(:,:,-1) = yfl%bx(:,:,yfl%nz-1)*syfu%clight
          syfu%bxz(:,:,-1) = 0.
          syfu%byx(:,:,-1) = 0.
          syfu%byz(:,:,-1) = yfl%by(:,:,yfl%nz-1)*syfu%clight
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=fu%proc) return
          yfu%bx(:,:,-1) = (syfl%bxy(:,:,syfl%nz-1)+syfl%bxz(:,:,syfl%nz-1))/syfl%clight
          yfu%by(:,:,-1) = (syfl%byx(:,:,syfl%nz-1)+syfl%byz(:,:,syfl%nz-1))/syfl%clight
          syfl%bxy(:,:,syfl%nz) = yfu%bx(:,:,0)*syfl%clight
          syfl%bxz(:,:,syfl%nz) = 0.
          syfl%byx(:,:,syfl%nz) = 0.
          syfl%byz(:,:,syfl%nz) = yfu%by(:,:,0)*syfl%clight
        case(splityeefield)
          syfu=>fu%syf
#ifdef MPIPARALLEL
          if (fl%proc/=my_index) then
            ! --- send data down in z
            ibuf = 7
            call mpi_packbuffer_init(2*size(syfu%by(:,:,:,0)),ibuf)
            call mympi_pack(syfu%bx(:,:,:,0),ibuf)
            call mympi_pack(syfu%by(:,:,:,0),ibuf)
            call mpi_isend_pack(fl%proc,3,ibuf)
          else if (fu%proc/=my_index) then
            ! --- send data up in z
            ibuf = 8
            call mpi_packbuffer_init(2*size(syfl%by(:,:,:,0)),ibuf)
            call mympi_pack(syfl%bx(:,:,:,syfl%nz-1),ibuf)
            call mympi_pack(syfl%by(:,:,:,syfl%nz-1),ibuf)
            call mpi_isend_pack(fu%proc,4,ibuf)
          else
#endif
          if (fl%proc/=fu%proc) return
          syfu%bxy(:,:,-1) = syfl%bxy(:,:,syfl%nz-1)
          syfu%bxz(:,:,-1) = syfl%bxz(:,:,syfl%nz-1)
          syfu%byx(:,:,-1) = syfl%byx(:,:,syfl%nz-1)
          syfu%byz(:,:,-1) = syfl%byz(:,:,syfl%nz-1)
          syfl%bxy(:,:,syfl%nz) = syfu%bxy(:,:,0)
          syfl%bxz(:,:,syfl%nz) = syfu%bxz(:,:,0)
          syfl%byx(:,:,syfl%nz) = syfu%byx(:,:,0)
          syfl%byz(:,:,syfl%nz) = syfu%byz(:,:,0)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select

  return
end subroutine em3d_exchange_bndb_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndb_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ibuf

          if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
            ! --- recv data from down in z
            ibuf = 9
            call mpi_packbuffer_init((3*yfu%nzguard-1)*size(yfu%bx(:,:,0)),ibuf)
            call mpi_recv_pack(fl%proc,2,ibuf)
            call submpi_unpack_real_array(yfu%bx(-yfu%nxguard,-yfu%nyguard,yfu%izming), &
                                     size(yfu%bx(:,:,yfu%izming:yfu%izmin-1)),ibuf)
            call submpi_unpack_real_array(yfu%by(-yfu%nxguard,-yfu%nyguard,yfu%izming), &
                                     size(yfu%by(:,:,yfu%izming:yfu%izmin-1)),ibuf)
            if (yfu%nzguard>1) &
            call submpi_unpack_real_array(yfu%bz(-yfu%nxguard,-yfu%nyguard,yfu%izming+1), &
                                     size(yfu%bz(:,:,yfu%izming+1:yfu%izmin-1)),ibuf)
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 10
            call mpi_packbuffer_init((3*yfl%nzguard-1)*size(yfl%bx(:,:,0)),ibuf)
            call mpi_recv_pack(fu%proc,1,ibuf)
            call submpi_unpack_real_array(yfl%bx(-yfl%nxguard,-yfl%nyguard,yfl%izmax), &
                                     size(yfl%bx(:,:,yfl%izmax:yfl%izmaxg-1)),ibuf)
            call submpi_unpack_real_array(yfl%by(-yfl%nxguard,-yfl%nyguard,yfl%izmax), &
                                     size(yfl%by(:,:,yfl%izmax:yfl%izmaxg-1)),ibuf)
            if (yfl%nzguard>1) &
            call submpi_unpack_real_array(yfl%bz(-yfl%nxguard,-yfl%nyguard,yfl%izmax+1), &
                                     size(yfl%bz(:,:,yfl%izmax+1:yfl%izmaxg-1)),ibuf)
          end if
      end select
    case(splityeefield)
      syfl=>fl%syf
      select case(fu%fieldtype)
        case(yeefield)
        case(splityeefield)
          syfu=>fu%syf
          if (fl%proc/=my_index) then
            ibuf = 11
            call mpi_packbuffer_init(2*size(syfu%bx(:,:,:,0)),ibuf)
            call mpi_recv_pack(fl%proc,4,ibuf)
            ! --- recv data from down in z
            call submpi_unpack_real_array(syfu%bxy(-syfu%nxguard,-syfu%nyguard,-1),size(syfu%bx(:,:,:,-1)),ibuf)
            call submpi_unpack_real_array(syfu%byx(-syfu%nxguard,-syfu%nyguard,-1),size(syfu%bx(:,:,:,-1)),ibuf)
          else if (fu%proc/=my_index) then
            ! --- recv data from up in z
            ibuf = 12
            call mpi_packbuffer_init(2*size(syfl%bx(:,:,:,0)),ibuf)
            call mpi_recv_pack(fu%proc,3,ibuf)
            call submpi_unpack_real_array(syfl%bxy(-syfl%nxguard,-syfl%nyguard,syfl%nz),size(syfl%bx(:,:,:,-1)),ibuf)
            call submpi_unpack_real_array(syfl%byx(-syfl%nxguard,-syfl%nyguard,syfl%nz),size(syfl%bx(:,:,:,-1)),ibuf)
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndb_zrecv
#endif

subroutine em3d_exchange_bndj_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix
#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf

!          yfu%J(0:1,:,:,2:3)    = yfu%J(0:1,:,:,2:3)    + yfl%J(yfl%nx:yfl%nx+1,:,:,2:3)
!          yfu%J(0,:,:,1)        = yfu%J(0,:,:,1)        + yfl%J(yfl%nx,:,:,1) 

!          yfl%J(yfl%nx-1,:,:,:) = yfl%J(yfl%nx-1,:,:,:) + yfu%J(-1,:,:,:)
!          yfl%J(yfl%nx,:,:,2:3) = yfu%J(0,:,:,2:3)

          ix = yfu%nxguard
          yfu%J(0:ix,:,:,2:3)    = yfu%J(0:ix,:,:,2:3)    + yfl%J(yfl%nx:yfl%nx+ix,:,:,2:3)
          yfu%J(0:ix-1,:,:,1)        = yfu%J(0:ix-1,:,:,1)        + yfl%J(yfl%nx:yfl%nx+ix-1,:,:,1) 

          yfl%J(yfl%nx-ix:yfl%nx-1,:,:,:) = yfl%J(yfl%nx-ix:yfl%nx-1,:,:,:) + yfu%J(-ix:-1,:,:,:)
          yfl%J(yfl%nx,:,:,2:3) = yfu%J(0,:,:,2:3)

      end select
  end select

  return
end subroutine em3d_exchange_bndj_x

subroutine em3d_exchange_bndj_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
!          yfu%J(:,0:1,:,1:3:2)    = yfu%J(:,0:1,:,1:3:2)  + yfl%J(:,yfl%ny:yfl%ny+1,:,1:3:2)
!          yfu%J(:,0,:,2)          = yfu%J(:,0,:,2)        + yfl%J(:,yfl%ny,:,2) 

!          yfl%J(:,yfl%ny-1,:,:)   = yfl%J(:,yfl%ny-1,:,:) + yfu%J(:,-1,:,:)
!          yfl%J(:,yfl%ny,:,1:3:2) = yfu%J(:,0,:,1:3:2)
          iy = yfu%nyguard
          yfu%J(:,0:iy,:,1:3:2)    = yfu%J(:,0:iy,:,1:3:2)    + yfl%J(:,yfl%ny:yfl%ny+iy,:,1:3:2)
          yfu%J(:,0:iy-1,:,2)      = yfu%J(:,0:iy-1,:,2)      + yfl%J(:,yfl%ny:yfl%ny+iy-1,:,2) 

          yfl%J(:,yfl%ny-iy:yfl%ny-1,:,:) = yfl%J(:,yfl%ny-iy:yfl%ny-1,:,:) + yfu%J(:,-iy:-1,:,:)
          yfl%J(:,yfl%ny,:,1:3:2) = yfu%J(:,0,:,1:3:2)
      end select
  end select

  return
end subroutine em3d_exchange_bndj_y

subroutine em3d_exchange_bndj_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz

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
            call mpi_packbuffer_init((3*yfu%nzguard+2)*size(yfu%J(:,:,-1,1)),1)
            do iz = yfu%izming,yfu%izmin
              call mympi_pack(yfu%J(:,:,iz,1),1)
              call mympi_pack(yfu%J(:,:,iz,2),1)
            end do
            do iz = yfu%izming,yfu%izmin-1
              call mympi_pack(yfu%J(:,:,iz,3),1)
            end do
            call mpi_isend_pack(fl%proc,1,1)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            call mpi_packbuffer_init((3*yfl%nzguard+2)*size(yfl%J(:,:,0,0)),2)
            do iz = yfl%izmax, yfl%izmaxg
              call mympi_pack(yfl%J(:,:,iz,1),2)
              call mympi_pack(yfl%J(:,:,iz,2),2)
            end do
            do iz = yfl%izmax, yfl%izmaxg-1
              call mympi_pack(yfl%J(:,:,iz,3),2)
            end do
            call mpi_isend_pack(fu%proc,2,2)

          else
#endif
!          yfu%J(:,:,0:1,1:2)    = yfu%J(:,:,0:1,1:2)    + yfl%J(:,:,yfl%nz:yfl%nz+1,1:2)
!          yfu%J(:,:,0,3)        = yfu%J(:,:,0,3)        + yfl%J(:,:,yfl%nz,3) 

!          yfl%J(:,:,yfl%nz-1,:) = yfl%J(:,:,yfl%nz-1,:) + yfu%J(:,:,-1,:)
!          yfl%J(:,:,yfl%nz,1:2) = yfu%J(:,:,0,1:2)
          iz = yfu%nzguard
          yfu%J(:,:,0:iz,1:2)    = yfu%J(:,:,0:iz,1:2)    + yfl%J(:,:,yfl%nz:yfl%nz+iz,1:2)
          yfu%J(:,:,0:iz-1,3)    = yfu%J(:,:,0:iz-1,3)    + yfl%J(:,:,yfl%nz:yfl%nz+iz-1,3) 

          yfl%J(:,:,yfl%nz-iz:yfl%nz-1,:) = yfl%J(:,:,yfl%nz-iz:yfl%nz-1,:) + yfu%J(:,:,-iz:-1,:)
          yfl%J(:,:,yfl%nz,1:2) = yfu%J(:,:,0,1:2)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndj_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndj_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            call mpi_packbuffer_init((3*yfu%nzguard+2)*size(yfu%J(:,:,0,1)),4)
            call mpi_recv_pack(fl%proc,2,4)
            do iz = 0,yfu%nzguard
              yfu%J(:,:,iz  ,1) = yfu%J(:,:,iz  ,1) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),4),shape(yfu%J(:,:,0,1)))
              yfu%J(:,:,iz  ,2) = yfu%J(:,:,iz  ,2) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),4),shape(yfu%J(:,:,0,1)))
            end do
            do iz = 0,yfu%nzguard-1
              yfu%J(:,:,iz  ,3) = yfu%J(:,:,iz  ,3) + reshape(mpi_unpack_real_array( size(yfu%J(:,:,0,1)),4),shape(yfu%J(:,:,0,1)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            call mpi_packbuffer_init((3*yfl%nzguard+2)*size(yfl%J(:,:,0,0)),3)
            call mpi_recv_pack(fu%proc,1,3)
            do iz = -yfl%nzguard,0
              yfl%J(:,:,yfl%nz+iz,1) = yfl%J(:,:,yfl%nz+iz,1) + reshape(mpi_unpack_real_array( size(yfl%J(:,:,yfl%nz-1,1)),3),&
                                                                                            shape(yfl%J(:,:,yfl%nz-1,1)))
              yfl%J(:,:,yfl%nz+iz,2) = yfl%J(:,:,yfl%nz+iz,2) + reshape(mpi_unpack_real_array( size(yfl%J(:,:,yfl%nz-1,2)),3),&
                                                                                            shape(yfl%J(:,:,yfl%nz-1,2)))
            end do
            do iz = -yfl%nzguard,-1
              yfl%J(:,:,yfl%nz+iz,3) = yfl%J(:,:,yfl%nz+iz,3) + reshape(mpi_unpack_real_array( size(yfl%J(:,:,yfl%nz-1,3)),3),&
                                                                                            shape(yfl%J(:,:,yfl%nz-1,3)))
            end do
          end if
      end select
  end select

  return
end subroutine em3d_exchange_bndj_zrecv
#endif

subroutine em3d_exchange_bndrho_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: ix

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          ix = yfu%nxguard
          yfu%Rho(0:ix,:,:)      = yfu%Rho(0:ix,:,:)      + yfl%Rho(yfl%nx:yfl%nx+ix,:,:)
          yfl%Rho(yfl%nx-ix:yfl%nx-1,:,:) = yfl%Rho(yfl%nx-ix:yfl%nx-1,:,:) + yfu%Rho(-ix:-1,:,:)
          yfl%Rho(yfl%nx,:,:)   = yfu%Rho(0,:,:)
      end select
  end select

  return
end subroutine em3d_exchange_bndrho_x

subroutine em3d_exchange_bndrho_y(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iy

#ifdef MPIPARALLEL
          if (fl%proc/=my_index .and. fu%proc/=my_index) return
#endif

 select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          iy = yfu%nyguard
          yfu%Rho(:,0:iy,:)      = yfu%Rho(:,0:iy,:)      + yfl%Rho(:,yfl%ny:yfl%ny+iy,:)
          yfl%Rho(:,yfl%ny-iy:yfl%ny-1,:) = yfl%Rho(:,yfl%ny-iy:yfl%ny-1,:) + yfu%Rho(:,-iy:-1,:)
          yfl%Rho(:,yfl%ny,:)   = yfu%Rho(:,0,:)
      end select
  end select

  return
end subroutine em3d_exchange_bndrho_y

subroutine em3d_exchange_bndrho_z(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz

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
            call mpi_packbuffer_init(size(yfu%rho(:,:,0:yfu%nzguard)),1)
            do iz = -yfu%nzguard,0
              call mympi_pack(yfu%rho(:,:,iz  ),1)
            end do
            call mpi_isend_pack(fl%proc,1,1)
            
          else if (fu%proc/=my_index) then

            ! --- send data up in z
            call mpi_packbuffer_init(size(yfl%rho(:,:,0:yfl%nzguard)),2)
            do iz = 0,yfl%nzguard
              call mympi_pack(yfl%rho(:,:,yfl%nz+iz  ),2)
            end do
            call mpi_isend_pack(fu%proc,2,2)

          else
#endif
          iz = yfu%nzguard
          yfu%Rho(:,:,0:iz)      = yfu%Rho(:,:,0:iz)      + yfl%Rho(:,:,yfl%nz:yfl%nz+iz)
          yfl%Rho(:,:,yfl%nz-iz:yfl%nz-1) = yfl%Rho(:,:,yfl%nz-iz:yfl%nz-1) + yfu%Rho(:,:,-iz:-1)
          yfl%Rho(:,:,yfl%nz)   = yfu%Rho(:,:,0)
#ifdef MPIPARALLEL
          end if
#endif
      end select
  end select
  return
end subroutine em3d_exchange_bndrho_z

#ifdef MPIPARALLEL
subroutine em3d_exchange_bndrho_zrecv(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu
integer(ISZ) :: iz

  if (fl%proc/=my_index .and. fu%proc/=my_index) return

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          if (fl%proc/=my_index) then
          
            ! --- recv data from down in z
            call mpi_packbuffer_init(size(yfu%rho(:,:,0:yfu%nzguard)),4)
            call mpi_recv_pack(fl%proc,2,4)
            do iz = 0,yfu%nzguard
              yfu%rho(:,:,iz  ) = yfu%rho(:,:,iz  ) + reshape(mpi_unpack_real_array( size(yfu%rho(:,:,0)),4),shape(yfu%rho(:,:,0)))
            end do

          else if (fu%proc/=my_index) then

            ! --- recv data from up in z
            call mpi_packbuffer_init(size(yfl%rho(:,:,0:yfl%nzguard)),3)
            call mpi_recv_pack(fu%proc,1,3)
            do iz = -yfl%nzguard,0
              yfl%rho(:,:,yfl%nz+iz  ) = yfl%rho(:,:,yfl%nz+iz  ) + reshape(mpi_unpack_real_array( size(yfl%rho(:,:,yfl%nz  )),3),&
                                                                                            shape(yfl%rho(:,:,yfl%nz  )))
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
 
  ! --- X
  ! core<--->sides
  call em3d_exchange_bnde_x(b%core,   b%sidexr)
  call em3d_exchange_bnde_x(b%sidexl, b%core)
  ! sides<--->edges
  call em3d_exchange_bnde_x(b%sideyl,   b%edgexryl)
  call em3d_exchange_bnde_x(b%edgexlyl, b%sideyl)
  call em3d_exchange_bnde_x(b%sideyr,   b%edgexryr)
  call em3d_exchange_bnde_x(b%edgexlyr, b%sideyr)
  call em3d_exchange_bnde_x(b%sidezl,   b%edgexrzl)
  call em3d_exchange_bnde_x(b%edgexlzl, b%sidezl)
  call em3d_exchange_bnde_x(b%sidezr,   b%edgexrzr)
  call em3d_exchange_bnde_x(b%edgexlzr, b%sidezr)
  ! edges<--->corners
  call em3d_exchange_bnde_x(b%edgeylzl,     b%cornerxrylzl)
  call em3d_exchange_bnde_x(b%cornerxlylzl, b%edgeylzl)
  call em3d_exchange_bnde_x(b%edgeyrzl,     b%cornerxryrzl)
  call em3d_exchange_bnde_x(b%cornerxlyrzl, b%edgeyrzl)
  call em3d_exchange_bnde_x(b%edgeylzr,     b%cornerxrylzr)
  call em3d_exchange_bnde_x(b%cornerxlylzr, b%edgeylzr)
  call em3d_exchange_bnde_x(b%edgeyrzr,     b%cornerxryrzr)
  call em3d_exchange_bnde_x(b%cornerxlyrzr, b%edgeyrzr)

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bnde_y(b%core,   b%sideyr)
  call em3d_exchange_bnde_y(b%sideyl, b%core)
  ! sides<--->edges
  call em3d_exchange_bnde_y(b%sidexl,   b%edgexlyr)
  call em3d_exchange_bnde_y(b%edgexlyl, b%sidexl)
  call em3d_exchange_bnde_y(b%sidexr,   b%edgexryr)
  call em3d_exchange_bnde_y(b%edgexryl, b%sidexr)
  call em3d_exchange_bnde_y(b%sidezl,   b%edgeyrzl)
  call em3d_exchange_bnde_y(b%edgeylzl, b%sidezl)
  call em3d_exchange_bnde_y(b%sidezr,   b%edgeyrzr)
  call em3d_exchange_bnde_y(b%edgeylzr, b%sidezr)
  ! edges<--->corners
  call em3d_exchange_bnde_y(b%edgexlzl,     b%cornerxlyrzl)
  call em3d_exchange_bnde_y(b%cornerxlylzl, b%edgexlzl)
  call em3d_exchange_bnde_y(b%edgexrzl,     b%cornerxryrzl)
  call em3d_exchange_bnde_y(b%cornerxrylzl, b%edgexrzl)
  call em3d_exchange_bnde_y(b%edgexlzr,     b%cornerxlyrzr)
  call em3d_exchange_bnde_y(b%cornerxlylzr, b%edgexlzr)
  call em3d_exchange_bnde_y(b%edgexrzr,     b%cornerxryrzr)
  call em3d_exchange_bnde_y(b%cornerxrylzr, b%edgexrzr)
  ! --- Z
  ! core<--->sides
  call em3d_exchange_bnde_z(b%core,   b%sidezr)
  call em3d_exchange_bnde_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%sidezl, b%core)
  call em3d_exchange_bnde_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bnde_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bnde_z(b%edgexlzl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgexlzl, b%sidexl)
  call em3d_exchange_bnde_zrecv(b%sidexl,   b%edgexlzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bnde_z(b%edgexrzl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgexrzl, b%sidexr)
  call em3d_exchange_bnde_zrecv(b%sidexr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bnde_z(b%edgeylzl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgeylzl, b%sideyl)
  call em3d_exchange_bnde_zrecv(b%sideyl,   b%edgeylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bnde_z(b%edgeyrzl, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%edgeyrzl, b%sideyr)
  call em3d_exchange_bnde_zrecv(b%sideyr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bnde_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bnde_z(b%cornerxlylzl, b%edgexlyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bnde_zrecv(b%edgexlyl,     b%cornerxlylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bnde_z(b%cornerxrylzl, b%edgexryl)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bnde_zrecv(b%edgexryl,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bnde_z(b%cornerxlyrzl, b%edgexlyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bnde_zrecv(b%edgexlyr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bnde_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bnde_z(b%cornerxryrzl, b%edgexryr)
#ifdef MPIPARALLEL
  call em3d_exchange_bnde_zrecv(b%cornerxryrzl, b%edgexryr)
  call em3d_exchange_bnde_zrecv(b%edgexryr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_e

subroutine em3d_exchange_b(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
 
  ! --- X
  ! core<--->sides
  call em3d_exchange_bndb_x(b%core,   b%sidexr)
  call em3d_exchange_bndb_x(b%sidexl, b%core)
  ! sides<--->edges
  call em3d_exchange_bndb_x(b%sideyl,   b%edgexryl)
  call em3d_exchange_bndb_x(b%edgexlyl, b%sideyl)
  call em3d_exchange_bndb_x(b%sideyr,   b%edgexryr)
  call em3d_exchange_bndb_x(b%edgexlyr, b%sideyr)
  call em3d_exchange_bndb_x(b%sidezl,   b%edgexrzl)
  call em3d_exchange_bndb_x(b%edgexlzl, b%sidezl)
  call em3d_exchange_bndb_x(b%sidezr,   b%edgexrzr)
  call em3d_exchange_bndb_x(b%edgexlzr, b%sidezr)
  ! edges<--->corners
  call em3d_exchange_bndb_x(b%edgeylzl,     b%cornerxrylzl)
  call em3d_exchange_bndb_x(b%cornerxlylzl, b%edgeylzl)
  call em3d_exchange_bndb_x(b%edgeyrzl,     b%cornerxryrzl)
  call em3d_exchange_bndb_x(b%cornerxlyrzl, b%edgeyrzl)
  call em3d_exchange_bndb_x(b%edgeylzr,     b%cornerxrylzr)
  call em3d_exchange_bndb_x(b%cornerxlylzr, b%edgeylzr)
  call em3d_exchange_bndb_x(b%edgeyrzr,     b%cornerxryrzr)
  call em3d_exchange_bndb_x(b%cornerxlyrzr, b%edgeyrzr)

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndb_y(b%core,   b%sideyr)
  call em3d_exchange_bndb_y(b%sideyl, b%core)
  ! sides<--->edges
  call em3d_exchange_bndb_y(b%sidexl,   b%edgexlyr)
  call em3d_exchange_bndb_y(b%edgexlyl, b%sidexl)
  call em3d_exchange_bndb_y(b%sidexr,   b%edgexryr)
  call em3d_exchange_bndb_y(b%edgexryl, b%sidexr)
  call em3d_exchange_bndb_y(b%sidezl,   b%edgeyrzl)
  call em3d_exchange_bndb_y(b%edgeylzl, b%sidezl)
  call em3d_exchange_bndb_y(b%sidezr,   b%edgeyrzr)
  call em3d_exchange_bndb_y(b%edgeylzr, b%sidezr)
  ! edges<--->corners
  call em3d_exchange_bndb_y(b%edgexlzl,     b%cornerxlyrzl)
  call em3d_exchange_bndb_y(b%cornerxlylzl, b%edgexlzl)
  call em3d_exchange_bndb_y(b%edgexrzl,     b%cornerxryrzl)
  call em3d_exchange_bndb_y(b%cornerxrylzl, b%edgexrzl)
  call em3d_exchange_bndb_y(b%edgexlzr,     b%cornerxlyrzr)
  call em3d_exchange_bndb_y(b%cornerxlylzr, b%edgexlzr)
  call em3d_exchange_bndb_y(b%edgexrzr,     b%cornerxryrzr)
  call em3d_exchange_bndb_y(b%cornerxrylzr, b%edgexrzr)

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndb_z(b%core,   b%sidezr)
  call em3d_exchange_bndb_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndb_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif
  ! sides<--->edges
  call em3d_exchange_bndb_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bndb_z(b%edgexlzl, b%sidexl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgexlzl, b%sidexl)
  call em3d_exchange_bndb_zrecv(b%sidexl,   b%edgexlzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bndb_z(b%edgexrzl, b%sidexr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgexrzl, b%sidexr)
  call em3d_exchange_bndb_zrecv(b%sidexr,   b%edgexrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bndb_z(b%edgeylzl, b%sideyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgeylzl, b%sideyl)
  call em3d_exchange_bndb_zrecv(b%sideyl,   b%edgeylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bndb_z(b%edgeyrzl, b%sideyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%edgeyrzl, b%sideyr)
  call em3d_exchange_bndb_zrecv(b%sideyr,   b%edgeyrzr)
  call mpi_waitall_requests()
#endif
  ! edges<--->corners
  call em3d_exchange_bndb_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bndb_z(b%cornerxlylzl, b%edgexlyl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bndb_zrecv(b%edgexlyl,     b%cornerxlylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bndb_z(b%cornerxrylzl, b%edgexryl)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bndb_zrecv(b%edgexryl,     b%cornerxrylzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bndb_z(b%cornerxlyrzl, b%edgexlyr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bndb_zrecv(b%edgexlyr,     b%cornerxlyrzr)
  call mpi_waitall_requests()
#endif
  call em3d_exchange_bndb_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bndb_z(b%cornerxryrzl, b%edgexryr)
#ifdef MPIPARALLEL
  call em3d_exchange_bndb_zrecv(b%cornerxryrzl, b%edgexryr)
  call em3d_exchange_bndb_zrecv(b%edgexryr,     b%cornerxryrzr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_b

subroutine em3d_exchange_j(b)
use mod_emfield3d
implicit none
TYPE(EM3D_BLOCKtype) :: b
#ifdef MPIPARALLEL
integer(MPIISZ)::mpirequest(2)
#endif
  ! --- X
  ! core<--->sides
  call em3d_exchange_bndj_x(b%core,   b%sidexr)
  if(b%xrbnd /= periodic) call em3d_exchange_bndj_x(b%sidexl, b%core)

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndj_y(b%core,   b%sideyr)
  if(b%yrbnd /= periodic) call em3d_exchange_bndj_y(b%sideyl, b%core)

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndj_z(b%core,   b%sidezr)
  if(b%zrbnd /= periodic) call em3d_exchange_bndj_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  if(b%zrbnd /= periodic) call em3d_exchange_bndj_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndj_zrecv(b%core,   b%sidezr)
  call mpi_waitall_requests()
#endif

  return
end subroutine em3d_exchange_j

subroutine em3d_exchange_rho(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b

  ! --- X
  ! core<--->sides
  call em3d_exchange_bndrho_x(b%core,   b%sidexr)
  if(b%xrbnd /= periodic) call em3d_exchange_bndrho_x(b%sidexl, b%core)

  ! --- Y
  ! core<--->sides
  call em3d_exchange_bndrho_y(b%core,   b%sideyr)
  if(b%yrbnd /= periodic) call em3d_exchange_bndrho_y(b%sideyl, b%core)

  ! --- Z
  ! core<--->sides
  call em3d_exchange_bndrho_z(b%core,   b%sidezr)
  if(b%zrbnd /= periodic) call em3d_exchange_bndrho_z(b%sidezl, b%core)
#ifdef MPIPARALLEL
  if(b%zrbnd /= periodic) call em3d_exchange_bndrho_zrecv(b%sidezl, b%core)
  call em3d_exchange_bndrho_zrecv(b%core,   b%sidezr)
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
!return
  do l=-1,f%nz
    do k=-1,f%ny
      do j=f%nx,0,-1
        f%ex(j,k,l)=0.5*(f%ex(j,k,l)+f%ex(j-1,k,l))
        f%by(j,k,l)=0.5*(f%by(j,k,l)+f%by(j-1,k,l))
        f%bz(j,k,l)=0.5*(f%bz(j,k,l)+f%bz(j-1,k,l))
      enddo
    enddo
  enddo

  do l=-1,f%nz
    do k=f%ny,0,-1
      do j=-1,f%nx
        f%ey(j,k,l)=0.5*(f%ey(j,k,l)+f%ey(j,k-1,l))
        f%bz(j,k,l)=0.5*(f%bz(j,k,l)+f%bz(j,k-1,l))
        f%bx(j,k,l)=0.5*(f%bx(j,k,l)+f%bx(j,k-1,l))
      enddo
    enddo
  enddo

  do l=f%nz,0,-1
    do k=-1,f%ny
      do j=-1,f%nx
        f%ez(j,k,l)=0.5*(f%ez(j,k,l)+f%ez(j,k,l-1))
        f%bx(j,k,l)=0.5*(f%bx(j,k,l)+f%bx(j,k,l-1))
        f%by(j,k,l)=0.5*(f%by(j,k,l)+f%by(j,k,l-1))
      enddo
    enddo
  enddo

  return
end subroutine yee2node3d

subroutine node2yee3d(f)
! puts EM field back from node to Yee grid
use mod_emfield3d
implicit none
TYPE(EM3D_YEEFIELDtype) :: f

INTEGER :: j,k,l
!return
  do l=0,f%nz
    do k=-1,f%ny
      do j=-1,f%nx
        f%ez(j,k,l)=2.*f%ez(j,k,l)-f%ez(j,k,l-1)
        f%bx(j,k,l)=2.*f%bx(j,k,l)-f%bx(j,k,l-1)
        f%by(j,k,l)=2.*f%by(j,k,l)-f%by(j,k,l-1)
      enddo
    enddo
  enddo

  do l=-1,f%nz
    do k=0,f%ny
      do j=-1,f%nx
        f%ey(j,k,l)=2.*f%ey(j,k,l)-f%ey(j,k-1,l)
        f%bz(j,k,l)=2.*f%bz(j,k,l)-f%bz(j,k-1,l)
        f%bx(j,k,l)=2.*f%bx(j,k,l)-f%bx(j,k-1,l)
      enddo
    enddo
  enddo

  do l=-1,f%nz
    do k=-1,f%ny
      do j=0,f%nx
        f%ex(j,k,l)=2.*f%ex(j,k,l)-f%ex(j-1,k,l)
        f%by(j,k,l)=2.*f%by(j,k,l)-f%by(j-1,k,l)
        f%bz(j,k,l)=2.*f%bz(j,k,l)-f%bz(j-1,k,l)
      enddo
    enddo
  enddo

  return
end subroutine node2yee3d

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
  do i=1,n
    f%incond(indx(1,i),indx(2,i),indx(3,i)) = .true.
  end do  

end subroutine set_incond
