#include "top.h"
!     Last change:  JLV   3 Jun 2004    0:17 am
!************* MODULE field  **********************************************

module mod_emfield3d
#ifdef MPIPARALLEL
use Parallel
#endif
use EM3D_BLOCKtypemodule
use EM3D_YEEFIELDtypemodule
use EM3D_SPLITYEEFIELDtypemodule
USE mod_bnd
USE mod_bnd_cummer, create_bnd_cummer => create_bnd, &
                    move_bnd_cummer => move_bnd, &
                    move_window_bnd_cummer => move_window_bnd , &
                    ijk_cummer => ijk
USE EM2D_FIELDobjects
!use GlobalVars
!use Picglb

implicit none

integer, parameter :: dirichlet=0, neumann=1, periodic=2, openbc=3, yeefield=-1 , splityeefield=-2

#ifdef MPIPARALLEL
include "mpif.h"
integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpierror
integer(MPIISZ):: mpirequest
integer(MPIISZ):: w
integer(MPIISZ):: messid 
#endif

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
    real(kind=8) :: sigmax(-sf%nxguard:sf%nx+sf%nxguard), sigmax_next(-sf%nxguard:sf%nx+sf%nxguard), lsigmax(-sf%nxguard:sf%nx+sf%nxguard), lsigmax_next(-sf%nxguard:sf%nx+sf%nxguard)
    real(kind=8) :: sigmay(-sf%nyguard:sf%ny+sf%nyguard), sigmay_next(-sf%nyguard:sf%ny+sf%nyguard), lsigmay(-sf%nyguard:sf%ny+sf%nyguard), lsigmay_next(-sf%nyguard:sf%ny+sf%nyguard)
    real(kind=8) :: sigmaz(-sf%nzguard:sf%nz+sf%nzguard), sigmaz_next(-sf%nzguard:sf%nz+sf%nzguard), lsigmaz(-sf%nzguard:sf%nz+sf%nzguard), lsigmaz_next(-sf%nzguard:sf%nz+sf%nzguard)
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
           call assign_coefs(bnd_cond,sf%agx(j-1),sf%bmgx(j-1),sf%bpgx(j-1),sf%clight*dt,sf%dx,lsigmax_next(j-1),lsigmax(j-1),sb_coef,which)
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
           call assign_coefs(bnd_cond,sf%agy(j-1),sf%bmgy(j-1),sf%bpgy(j-1),sf%clight*dt,sf%dy,lsigmay_next(j-1),lsigmay(j-1),sb_coef,which)
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
           call assign_coefs(bnd_cond,sf%agz(j-1),sf%bmgz(j-1),sf%bpgz(j-1),sf%clight*dt,sf%dz,lsigmaz_next(j-1),lsigmaz(j-1),sb_coef,which)
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
    call EM3D_SPLITYEEFIELDtypeallot(sf)

    return 
  end subroutine init_splitfield
    

!************* SUBROUTINE init_fields  *************************************************
subroutine init_3dem_block(b, nx, ny, nz, nbndx, nbndy, nbndz, nxguard, nyguard, nzguard, dt, dx, dy, dz, clight, mu0, xmin, ymin, zmin, xlb, ylb, zlb, xrb, yrb, zrb)
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
  call init_splitfield(b%sidexl%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%sidexr%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! y
  call init_splitfield(b%sideyl%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%sideyr%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! z
  call init_splitfield(b%sidezl%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%sidezr%syf,nbndx,ny,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)

! *** edges
! xy
  call init_splitfield(b%edgexlyl%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexryl%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexlyr%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexryr%syf,nbndx,nbndy,nz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 0, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! xz
  call init_splitfield(b%edgexlzl%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexrzl%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexlzr%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgexrzr%syf,nbndx,ny,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 0, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
! yz
  call init_splitfield(b%edgeylzl%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgeyrzl%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgeylzr%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%edgeyrzr%syf,nx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 0, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)

! *** corners
  call init_splitfield(b%cornerxlylzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxrylzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxlyrzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxryrzl%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1,-1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxlylzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxrylzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1,-1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxlyrzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight,-1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)
  call init_splitfield(b%cornerxryrzr%syf,nbndx,nbndy,nbndz,nxguard,nyguard,nzguard, dt, dx, dy, dz, clight, 1, 1, 1, nnx, smaxx, sdeltax, nny, smaxy, sdeltay, nnz, smaxz, sdeltaz)


return

END subroutine init_3dem_block

subroutine push_em3d_e(f,dt)
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
      f%Ez(j,k,l) = f%Ez(j,k,l) + dtsdx * (f%By(j,k,l) - f%By(j-1,k  ,l)) &
                                - dtsdy * (f%Bx(j,k,l) - f%Bx(j  ,k-1,l)) &
                                - mudt  * f%J(j,k,l,3)
    end do
   end do
  end do

return
end subroutine push_em3d_e

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

  ! advance Bx
  do l = 0, f%nz-1
   do k = 0, f%ny-1
    do j = 0, f%nx
      f%Bx(j,k,l) = f%Bx(j,k,l) - dtsdy * (f%Ez(j,k+1,l  ) - f%Ez(j,k,l)) &
                                + dtsdz * (f%Ey(j,k,  l+1) - f%Ey(j,k,l))
    end do
   end do
  end do

  ! advance By
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx-1
      f%By(j,k,l) = f%By(j,k,l) + dtsdx * (f%Ez(j+1,k,l  ) - f%Ez(j,k,l)) &  
                                - dtsdz * (f%Ex(j  ,k,l+1) - f%Ex(j,k,l)) 
    end do
   end do
  end do

  ! advance Bz 
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx-1
      f%Bz(j,k,l) = f%Bz(j,k,l) - dtsdx * (f%Ey(j+1,k,l) - f%Ey(j,k,l)) &
                                + dtsdy * (f%Ex(j,k+1,l) - f%Ex(j,k,l))
    end do
   end do
  end do

return
end subroutine push_em3d_b

subroutine push_em3d_f(f,dt)
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
      f%F(j,k,l) = f%F(j,k,l) + dtsdx * (f%Ex(j,k,l) - f%Ex(j-1,k  ,l  )) &
                              + dtsdy * (f%Ey(j,k,l) - f%Ey(j  ,k-1,l  )) &
                              + dtsdz * (f%Ez(j,k,l) - f%Ez(j  ,k  ,l-1)) &
                              - dtsepsi * f%Rho(j,k,l)
    end do
   end do
  end do

end subroutine push_em3d_f

subroutine push_em3d_ef(f,dt)
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
      f%Ex(j,k,l) = f%Ex(j,k,l) + dtsdx * (f%F(j+1,k,l) - f%F(j,k,l)) 
    end do
   end do
  end do

  ! advance Ey
  do l = 0, f%nz
   do k = 0, f%ny-1
    do j = 0, f%nx
      f%Ey(j,k,l) = f%Ey(j,k,l) + dtsdy * (f%F(j,k+1,l) - f%F(j,k,l))
    end do
   end do
  end do

  ! advance Ez 
  do l = 0, f%nz-1
   do k = 0, f%ny
    do j = 0, f%nx
      f%Ez(j,k,l) = f%Ez(j,k,l) + dtsdz * (f%F(j,k,l+1) - f%F(j,k,l)) 
    end do
   end do
  end do

return
end subroutine push_em3d_ef

subroutine push_em3d_splite(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  do l = 0, sf%nz
   do k = 0, sf%ny
    do j = 0, sf%nx-1
      sf%exy(j,k,l) = sf%afy(k)*sf%exy(j,k,l) + sf%bpfy(k)*(sf%bzx(j,k,l)+sf%bzy(j,k,l))  + sf%bmfy(k)*(sf%bzx(j,k-1,l)+sf%bzy(j,k-1,l)) !- 0.5_8*dt*sf%j(j,k,l,1)
      sf%exz(j,k,l) = sf%afz(l)*sf%exz(j,k,l) - sf%bpfz(l)*(sf%byx(j,k,l)+sf%byz(j,k,l))  - sf%bmfz(l)*(sf%byx(j,k,l-1)+sf%byz(j,k,l-1)) !- 0.5_8*dt*sf%j(j,k,l,1)
    end do
   end do
  end do

  do l = 0, sf%nz
   do k = 0, sf%ny-1
    do j = 0, sf%nx
      sf%eyx(j,k,l) = sf%afx(j)*sf%eyx(j,k,l) - sf%bpfx(j)*(sf%bzx(j,k,l)+sf%bzy(j,k,l))  - sf%bmfx(j)*(sf%bzx(j-1,k,l)+sf%bzy(j-1,k,l)) !- 0.5_8*dt*sf%j(j,k,l,2)
      sf%eyz(j,k,l) = sf%afz(l)*sf%eyz(j,k,l) + sf%bpfz(l)*(sf%bxy(j,k,l)+sf%bxz(j,k,l))  + sf%bmfz(l)*(sf%bxy(j,k,l-1)+sf%bxz(j,k,l-1)) !- 0.5_8*dt*sf%j(j,k,l,2)
    end do
   end do
  end do

  do l = 0, sf%nz-1
   do k = 0, sf%ny
    do j = 0, sf%nx
      sf%ezx(j,k,l) = sf%afx(j)*sf%ezx(j,k,l) + sf%bpfx(j)*(sf%byx(j,k,l)+sf%byz(j,k,l))  + sf%bmfx(j)*(sf%byx(j-1,k,l)+sf%byz(j-1,k,l)) !- 0.5_8*dt*sf%j(j,k,l,3)
      sf%ezy(j,k,l) = sf%afy(k)*sf%ezy(j,k,l) - sf%bpfy(k)*(sf%bxy(j,k,l)+sf%bxz(j,k,l))  - sf%bmfy(k)*(sf%bxy(j,k-1,l)+sf%bxz(j,k-1,l)) !- 0.5_8*dt*sf%j(j,k,l,3)
    end do
   end do
  end do

  return
end subroutine push_em3d_splite

subroutine push_em3d_splitef(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  do l = 0, sf%nz
   do k = 0, sf%ny
    do j = 0, sf%nx-1
      sf%exx(j,k,l) = sf%agx(j)*sf%exx(j,k,l) + sf%bpgx(j)*( sf%fx(j+1,k,l) + sf%fy(j+1,k,l) + sf%fz(j+1,k,l) ) &
                                              + sf%bmgx(j)*( sf%fx(j  ,k,l) + sf%fy(j  ,k,l) + sf%fz(j  ,k,l) )
    end do
   end do
  end do

  do l = 0, sf%nz
   do k = 0, sf%ny-1
    do j = 0, sf%nx
      sf%eyy(j,k,l) = sf%agy(k)*sf%eyy(j,k,l) + sf%bpgy(k)*( sf%fx(j,k+1,l) + sf%fy(j,k+1,l) + sf%fz(j,k+1,l) ) &
                                              + sf%bmgy(k)*( sf%fx(j,k  ,l) + sf%fy(j,k  ,l) + sf%fz(j,k  ,l) )
    end do
   end do
  end do

  do l = 0, sf%nz-1
   do k = 0, sf%ny
    do j = 0, sf%nx
      sf%ezz(j,k,l) = sf%agz(l)*sf%ezz(j,k,l) + sf%bpgz(l)*( sf%fx(j,k,l+1) + sf%fy(j,k,l+1) + sf%fz(j,k,l+1) ) &
                                              + sf%bmgz(l)*( sf%fx(j,k,l)   + sf%fy(j,k,l  ) + sf%fz(j,k,l  ) )
    end do
   end do
  end do

  return
end subroutine push_em3d_splitef

subroutine push_em3d_splitb(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  do l = 0, sf%nz-1
   do k = 0, sf%ny-1
    do j = 0, sf%nx
      sf%bxy(j,k,l) = sf%agy(k)*sf%bxy(j,k,l) - sf%bpgy(k)*(sf%ezx(j,k+1,l  )+sf%ezy(j,k+1,l  )+sf%ezz(j,k+1,l  )) &
                                              - sf%bmgy(k)*(sf%ezx(j,k  ,l  )+sf%ezy(j,k  ,l  )+sf%ezz(j,k  ,l  ))
      sf%bxz(j,k,l) = sf%agz(l)*sf%bxz(j,k,l) + sf%bpgz(l)*(sf%eyx(j,k  ,l+1)+sf%eyy(j,k  ,l+1)+sf%eyz(j,k  ,l+1)) &
                                              + sf%bmgz(l)*(sf%eyx(j,k  ,l  )+sf%eyy(j,k  ,l  )+sf%eyz(j,k  ,l  ))
    end do
   end do
  end do

  do l = 0, sf%nz-1
   do k = 0, sf%ny
    do j = 0, sf%nx-1
      sf%byx(j,k,l) = sf%agx(j)*sf%byx(j,k,l) + sf%bpgx(j)*(sf%ezx(j+1,k,l  )+sf%ezy(j+1,k,l  )+sf%ezz(j+1,k,l  )) &
                                              + sf%bmgx(j)*(sf%ezx(j  ,k,l  )+sf%ezy(j  ,k,l  )+sf%ezz(j  ,k,l  ))
      sf%byz(j,k,l) = sf%agz(l)*sf%byz(j,k,l) - sf%bpgz(l)*(sf%exx(j  ,k,l+1)+sf%exy(j  ,k,l+1)+sf%exz(j  ,k,l+1)) &
                                              - sf%bmgz(l)*(sf%exx(j  ,k,l  )+sf%exy(j  ,k,l  )+sf%exz(j  ,k,l  ))
    end do
   end do
  end do

  do l = 0, sf%nz
   do k = 0, sf%ny-1
    do j = 0, sf%nx-1
      sf%bzx(j,k,l) = sf%agx(j)*sf%bzx(j,k,l) - sf%bpgx(j)*(sf%eyx(j+1,k  ,l)+sf%eyy(j+1,k  ,l)+sf%eyz(j+1,k  ,l)) &
                                              - sf%bmgx(j)*(sf%eyx(j  ,k  ,l)+sf%eyy(j  ,k  ,l)+sf%eyz(j  ,k  ,l))
      sf%bzy(j,k,l) = sf%agy(k)*sf%bzy(j,k,l) + sf%bpgy(k)*(sf%exx(j  ,k+1,l)+sf%exy(j  ,k+1,l)+sf%exz(j  ,k+1,l)) &
                                              + sf%bmgy(k)*(sf%exx(j  ,k  ,l)+sf%exy(j  ,k  ,l)+sf%exz(j  ,k  ,l))
    end do
   end do
  end do

  return
end subroutine push_em3d_splitb

subroutine push_em3d_splitf(sf,dt,which)
use mod_emfield3d
implicit none

TYPE(EM3D_SPLITYEEFIELDtype) :: sf
REAL(kind=8), INTENT(IN) :: dt

INTEGER :: j, k, l,which

  call set_bndcoeffsem3d(sf,dt,which)

  do l = 0, sf%nz
   do k = 0, sf%ny
    do j = 0, sf%nx
     sf%fx(j,k,l) = sf%afx(j)*sf%fx(j,k,l) + sf%bpfx(j)*(sf%exx(j  ,k  ,l  )+sf%exy(j  ,k  ,l  )+sf%exz(j  ,k  ,l  )) &
                                           + sf%bmfx(j)*(sf%exx(j-1,k  ,l  )+sf%exy(j-1,k  ,l  )+sf%exz(j-1,k  ,l  )) !- (1._8/3._8)*dt*sf%rho(j,k,l)

     sf%fy(j,k,l) = sf%afy(k)*sf%fy(j,k,l) + sf%bpfy(k)*(sf%eyx(j  ,k  ,l  )+sf%eyy(j  ,k  ,l  )+sf%eyz(j  ,k  ,l  )) &
                                           + sf%bmfy(k)*(sf%eyx(j  ,k-1,l  )+sf%eyy(j  ,k-1,l  )+sf%eyz(j  ,k-1,l  )) !- (1._8/3._8)*dt*sf%rho(j,k,l)

     sf%fz(j,k,l) = sf%afz(l)*sf%fz(j,k,l) + sf%bpfz(l)*(sf%ezx(j  ,k  ,l  )+sf%ezy(j  ,k  ,l  )+sf%ezz(j  ,k  ,l  )) &
                                           + sf%bmfz(l)*(sf%ezx(j  ,k  ,l-1)+sf%ezy(j  ,k  ,l-1)+sf%ezz(j  ,k  ,l-1)) !- (1._8/3._8)*dt*sf%rho(j,k,l)
    end do
   end do
  end do

  return
end subroutine push_em3d_splitf

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

  return
end subroutine push_em3d_block

subroutine push_em3d_eef(b,dt,which,l_pushf)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
logical(ISZ) :: l_pushf

INTEGER(ISZ) :: j, k, l,which
  
  if (which==0) then
    if(l_pushf) call push_em3d_ef(b%core%yf,dt)
    call push_em3d_e(b%core%yf,dt)
  else
    if(l_pushf) call push_em3d_ef(b%core%yf,dt*0.5)
    call push_em3d_e(b%core%yf,dt*0.5)
  end if
  if(l_pushf) call push_em3d_blockbndef(b,dt,which)
  call push_em3d_blockbnde(b,dt,which)
  
  if(l_pushf) call em3d_exchange_e(b)

  return
end subroutine push_em3d_eef

subroutine push_em3d_bf(b,dt,which,l_pushf)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
REAL(kind=8), INTENT(IN) :: dt
logical(ISZ) :: l_pushf

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

subroutine em3d_exchange_bnde_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%ex(yfl%nx,:,:) = yfu%ex(0,:,:)
          yfu%ex(-1,:,:)     = yfl%ex(yfl%nx-1,:,:)
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

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%ey(:,yfl%ny,:) = yfu%ey(:,0,:)
          yfu%ey(:,-1,:)     = yfl%ey(:,yfl%ny-1,:)
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

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%ez(:,:,yfl%nz) = yfu%ez(:,:,0)
          yfu%ez(:,:,-1)     = yfl%ez(:,:,yfl%nz-1)
        case(splityeefield)
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
          yfu=>fu%yf
          yfu%ez(:,:,-1)        = syfl%ezx(:,:,syfl%nz-1)+syfl%ezy(:,:,syfl%nz-1)+syfl%ezz(:,:,syfl%nz-1)
          syfl%ezx(:,:,syfl%nz) = 0.
          syfl%ezy(:,:,syfl%nz) = 0.
          syfl%ezz(:,:,syfl%nz) = yfu%ez(:,:,0)
        case(splityeefield)
          syfu=>fu%syf
          syfu%ezx(:,:,-1)      = syfl%ezx(:,:,syfl%nz-1)
          syfu%ezy(:,:,-1)      = syfl%ezy(:,:,syfl%nz-1)
          syfu%ezz(:,:,-1)      = syfl%ezz(:,:,syfl%nz-1)
          syfl%ezx(:,:,syfl%nz) = syfu%ezx(:,:,0)
          syfl%ezy(:,:,syfl%nz) = syfu%ezy(:,:,0)
          syfl%ezz(:,:,syfl%nz) = syfu%ezz(:,:,0)
      end select
  end select

  return
end subroutine em3d_exchange_bnde_z

subroutine em3d_exchange_bndb_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%by(yfl%nx,:,:) = yfu%by(0,:,:)
          yfl%bz(yfl%nx,:,:) = yfu%bz(0,:,:)
          yfu%by(-1,:,:)     = yfl%by(yfl%nx-1,:,:)
          yfu%bz(-1,:,:)     = yfl%bz(yfl%nx-1,:,:)
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

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%bx(:,yfl%ny,:) = yfu%bx(:,0,:)
          yfl%bz(:,yfl%ny,:) = yfu%bz(:,0,:)
          yfu%bx(:,-1,:)     = yfl%bx(:,yfl%ny-1,:)
          yfu%bz(:,-1,:)     = yfl%bz(:,yfl%ny-1,:)
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

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfl%bx(:,:,yfl%nz) = yfu%bx(:,:,0)
          yfl%by(:,:,yfl%nz) = yfu%by(:,:,0)
          yfu%bx(:,:,-1)     = yfl%bx(:,:,yfl%nz-1)
          yfu%by(:,:,-1)     = yfl%by(:,:,yfl%nz-1)
        case(splityeefield)
          syfu=>fu%syf
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
          yfu%bx(:,:,-1) = (syfl%bxy(:,:,syfl%nz-1)+syfl%bxz(:,:,syfl%nz-1))/syfl%clight
          yfu%by(:,:,-1) = (syfl%byx(:,:,syfl%nz-1)+syfl%byz(:,:,syfl%nz-1))/syfl%clight
          syfl%bxy(:,:,syfl%nz) = yfu%bx(:,:,0)*syfl%clight
          syfl%bxz(:,:,syfl%nz) = 0.
          syfl%byx(:,:,syfl%nz) = 0.
          syfl%byz(:,:,syfl%nz) = yfu%by(:,:,0)*syfl%clight
        case(splityeefield)
          syfu=>fu%syf
          syfu%bxy(:,:,-1) = syfl%bxy(:,:,syfl%nz-1)
          syfu%bxz(:,:,-1) = syfl%bxz(:,:,syfl%nz-1)
          syfu%byx(:,:,-1) = syfl%byx(:,:,syfl%nz-1)
          syfu%byz(:,:,-1) = syfl%byz(:,:,syfl%nz-1)
          syfl%bxy(:,:,syfl%nz) = syfu%bxy(:,:,0)
          syfl%bxz(:,:,syfl%nz) = syfu%bxz(:,:,0)
          syfl%byx(:,:,syfl%nz) = syfu%byx(:,:,0)
          syfl%byz(:,:,syfl%nz) = syfu%byz(:,:,0)
      end select
  end select

  return
end subroutine em3d_exchange_bndb_z

subroutine em3d_exchange_bndj_x(fl,fu)
use mod_emfield3d
implicit none

TYPE(EM3D_FIELDtype) :: fl, fu
TYPE(EM3D_YEEFIELDtype), pointer :: yfl, yfu
TYPE(EM3D_SPLITYEEFIELDtype), pointer :: syfl, syfu

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%J(0,:,:,1)          = yfu%J(0,:,:,1)          + yfl%J(yfl%nx,:,:,1) 
          yfl%J(yfl%nx-1,:,:,1)   = yfl%J(yfl%nx-1,:,:,1)   + yfu%J(-1,:,:,1)
          yfu%J(0:1,:,:,2:3)      = yfu%J(0:1,:,:,2:3)      + yfl%J(yfl%nx:yfl%nx+1,:,:,2:3)
          yfl%J(yfl%nx-1,:,:,2:3) = yfl%J(yfl%nx-1,:,:,2:3) + yfu%J(-1,:,:,2:3)
          yfl%J(yfl%nx,:,:,2:3)   = yfu%J(0,:,:,2:3)
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

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%J(:,0,:,2)          = yfu%J(:,0,:,2)          + yfl%J(:,yfl%ny,:,2) 
          yfl%J(:,yfl%ny-1,:,2)   = yfl%J(:,yfl%ny-1,:,2)   + yfu%J(:,-1,:,2)
          yfu%J(:,0:1,:,1:3:2)      = yfu%J(:,0:1,:,1:3:2)      + yfl%J(:,yfl%ny:yfl%ny+1,:,1:3:2)
          yfl%J(:,yfl%ny-1,:,1:3:2) = yfl%J(:,yfl%ny-1,:,1:3:2) + yfu%J(:,-1,:,1:3:2)
          yfl%J(:,yfl%ny,:,1:3:2)   = yfu%J(:,0,:,1:3:2)
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

  select case(fl%fieldtype)
    case(yeefield)
      yfl=>fl%yf
      select case(fu%fieldtype)
        case(yeefield)
          yfu=>fu%yf
          yfu%J(:,:,0,3)          = yfu%J(:,:,0,3)          + yfl%J(:,:,yfl%nz,3) 
          yfl%J(:,:,yfl%nz-1,3)   = yfl%J(:,:,yfl%nz-1,3)   + yfu%J(:,:,-1,3)
          yfu%J(:,:,0:1,1:2)      = yfu%J(:,:,0:1,1:2)      + yfl%J(:,:,yfl%nz:yfl%nz+1,1:2)
          yfl%J(:,:,yfl%nz-1,1:2) = yfl%J(:,:,yfl%nz-1,1:2) + yfu%J(:,:,-1,1:2)
          yfl%J(:,:,yfl%nz,1:2)   = yfu%J(:,:,0,1:2)
      end select
  end select

  return
end subroutine em3d_exchange_bndj_z

subroutine em3d_exchange_e(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
!  return
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
  ! sides<--->edges
  call em3d_exchange_bnde_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bnde_z(b%edgexlzl, b%sidexl)
  call em3d_exchange_bnde_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bnde_z(b%edgexrzl, b%sidexr)
  call em3d_exchange_bnde_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bnde_z(b%edgeylzl, b%sideyl)
  call em3d_exchange_bnde_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bnde_z(b%edgeyrzl, b%sideyr)
!  end if
!  if (.false.) then
  ! edges<--->corners
  call em3d_exchange_bnde_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bnde_z(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bnde_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bnde_z(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bnde_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bnde_z(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bnde_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bnde_z(b%cornerxryrzl, b%edgexryr)

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
  ! sides<--->edges
  call em3d_exchange_bndb_z(b%sidexl,   b%edgexlzr)
  call em3d_exchange_bndb_z(b%edgexlzl, b%sidexl)
  call em3d_exchange_bndb_z(b%sidexr,   b%edgexrzr)
  call em3d_exchange_bndb_z(b%edgexrzl, b%sidexr)
  call em3d_exchange_bndb_z(b%sideyl,   b%edgeylzr)
  call em3d_exchange_bndb_z(b%edgeylzl, b%sideyl)
  call em3d_exchange_bndb_z(b%sideyr,   b%edgeyrzr)
  call em3d_exchange_bndb_z(b%edgeyrzl, b%sideyr)
  ! edges<--->corners
  call em3d_exchange_bndb_z(b%edgexlyl,     b%cornerxlylzr)
  call em3d_exchange_bndb_z(b%cornerxlylzl, b%edgexlyl)
  call em3d_exchange_bndb_z(b%edgexryl,     b%cornerxrylzr)
  call em3d_exchange_bndb_z(b%cornerxrylzl, b%edgexryl)
  call em3d_exchange_bndb_z(b%edgexlyr,     b%cornerxlyrzr)
  call em3d_exchange_bndb_z(b%cornerxlyrzl, b%edgexlyr)
  call em3d_exchange_bndb_z(b%edgexryr,     b%cornerxryrzr)
  call em3d_exchange_bndb_z(b%cornerxryrzl, b%edgexryr)

  return
end subroutine em3d_exchange_b

subroutine em3d_exchange_j(b)
use mod_emfield3d
implicit none

TYPE(EM3D_BLOCKtype) :: b
!  return
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

  return
end subroutine em3d_exchange_j

subroutine yee2node3d(f)
! puts EM value from Yee grid to nodes
use mod_emfield3d
implicit none
TYPE(EM3D_YEEFIELDtype) :: f

INTEGER :: j,k,l
!return
  do l=0,f%nz
    do k=0,f%ny
      do j=f%nx,0,-1
        f%ex(j,k,l)=0.5*(f%ex(j,k,l)+f%ex(j-1,k,l))
        f%by(j,k,l)=0.5*(f%by(j,k,l)+f%by(j-1,k,l))
        f%bz(j,k,l)=0.5*(f%bz(j,k,l)+f%bz(j-1,k,l))
      enddo
    enddo
  enddo

  do l=0,f%nz
    do k=f%ny,0,-1
      do j=0,f%nx
        f%ey(j,k,l)=0.5*(f%ey(j,k,l)+f%ey(j,k-1,l))
        f%bz(j,k,l)=0.5*(f%bz(j,k,l)+f%bz(j,k-1,l))
        f%bx(j,k,l)=0.5*(f%bx(j,k,l)+f%bx(j,k-1,l))
      enddo
    enddo
  enddo

  do l=f%nz,0,-1
    do k=0,f%ny
      do j=0,f%nx
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
    do k=0,f%ny
      do j=0,f%nx
        f%ez(j,k,l)=2.*f%ez(j,k,l)-f%ez(j,k,l-1)
        f%bx(j,k,l)=2.*f%bx(j,k,l)-f%bx(j,k,l-1)
        f%by(j,k,l)=2.*f%by(j,k,l)-f%by(j,k,l-1)
      enddo
    enddo
  enddo

  do l=0,f%nz
    do k=0,f%ny
      do j=0,f%nx
        f%ey(j,k,l)=2.*f%ey(j,k,l)-f%ey(j,k-1,l)
        f%bz(j,k,l)=2.*f%bz(j,k,l)-f%bz(j,k-1,l)
        f%bx(j,k,l)=2.*f%bx(j,k,l)-f%bx(j,k-1,l)
      enddo
    enddo
  enddo

  do l=0,f%nz
    do k=0,f%ny
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
