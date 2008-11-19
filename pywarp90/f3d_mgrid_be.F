#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID_BE.M,
c version $Revision: 1.43 $, $Date: 2008/11/19 18:29:49 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid field sovler which is part of the F3D
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigridbe3df(iwhich,nx,ny,nz,dx,dy,dz,phi,rho,
     &                          rstar,linbend,
     &                          bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                          xmmin,ymmin,zmmin,zbeam,zgrid)
      use GlobalVars
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use BoltzmannElectrons,Only: useriondensitygrid3d
      use Parallel
c     use ifcore
c     integer:: fsf
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      integer(ISZ):: bound0,boundnz,boundxy
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid

c     fsf = for_set_fpe(FPE_M_TRAP_INV+FPE_M_TRAP_DIV0+FPE_M_TRAP_OVF)

c     --- copy boundary positions from bound0, boundnz, and boundxy
      bounds(0) = boundxy
      bounds(1) = boundxy
      bounds(2) = boundxy
      bounds(3) = boundxy
      bounds(4) = bound0
      bounds(5) = boundnz
      if (l2symtry) then
        bounds(2) = neumann
        if (boundxy == 2) bounds(3) = neumann
      else if (l4symtry) then
        bounds(0) = neumann
        bounds(2) = neumann
        if (boundxy == 2) bounds(1) = neumann
        if (boundxy == 2) bounds(3) = neumann
      endif

      call setupiondensitygrid3d(xmmin,ymmin,zmmin,dx,dy,dz,nx,ny,nz,
     &                           rho,useriondensitygrid3d)

      call multigridbe3dsolve(iwhich,nx,ny,nz,dx,dy,dz,phi,rho,
     &                        rstar,linbend,bounds,
     &                        xmmin,ymmin,zmmin,zbeam,zgrid,
     &                        mgparam,mgiters,mgmaxiters,
     &                        mgmaxlevels,mgerror,mgtol,mgverbose,
     &                        downpasses,uppasses,
     &                        lcndbndy,laddconductor,icndbndy,
     &                        gridmode,conductors,useriondensitygrid3d,
     &                        fsdecomp)

      return
      end
c=============================================================================
      subroutine setupiondensitygrid3d(xmmin,ymmin,zmmin,dx,dy,dz,nx,ny,nz,
     &                                 rho,iondensitygrid3d)
      use BoltzmannElectrons
      use Grid3dtypemodule
      real(kind=8):: xmmin,ymmin,zmmin
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      type(Grid3dtype):: iondensitygrid3d

      integer(ISZ):: ixbemin,ixbemax,iybemin,iybemax,izbemin,izbemax
      integer(ISZ):: id,ix,iy,iz
      logical(ISZ):: lanyuseparticledensity,check

c     --- Check for luseparticledensity
      lanyuseparticledensity = .false.
      do id=1,nberegions
        if (luseparticledensity(id) > 0) then
          lanyuseparticledensity = .true.
          liondensitygrid3d(id) = .true.
c         --- Find the extrema of all of the regions
          iondensitygrid3d%xmin = min(xbemin(id),iondensitygrid3d%xmin)
          iondensitygrid3d%xmax = max(xbemax(id),iondensitygrid3d%xmax)
          iondensitygrid3d%ymin = min(ybemin(id),iondensitygrid3d%ymin)
          iondensitygrid3d%ymax = max(ybemax(id),iondensitygrid3d%ymax)
          iondensitygrid3d%zmin = min(zbemin(id),iondensitygrid3d%zmin)
          iondensitygrid3d%zmax = max(zbemax(id),iondensitygrid3d%zmax)
c         --- Make sure that are within the grid bounds
          iondensitygrid3d%xmin = max(xmmin      ,iondensitygrid3d%xmin)
          iondensitygrid3d%xmax = min(xmmin+nx*dx,iondensitygrid3d%xmax)
          iondensitygrid3d%ymin = max(ymmin      ,iondensitygrid3d%ymin)
          iondensitygrid3d%ymax = min(ymmin+ny*dy,iondensitygrid3d%ymax)
          iondensitygrid3d%zmin = max(zmmin      ,iondensitygrid3d%zmin)
          iondensitygrid3d%zmax = min(zmmin+nz*dz,iondensitygrid3d%zmax)
        endif
        logelecdenmaxscale(id) = log(electrondensitymaxscale(id))
      enddo

      if (lanyuseparticledensity) then
c       --- Force the min's and max's to coincide with field grid values
c       --- Note that the min's and max's are relative to the field grid,
c       --- so for example, the location if xmin in the lab frame is
c       --- iondensitygrid3d%xmin + grid%rmin. This is done so that
c       --- grid%rmin etc are not needed when iondensitygrid3d is used.
        ix = (floor((iondensitygrid3d%xmin - xmmin)/dx))
        iondensitygrid3d%xmin = dx*ix
        iy = (floor((iondensitygrid3d%ymin - ymmin)/dy))
        iondensitygrid3d%ymin = dy*iy
        iz = (floor((iondensitygrid3d%zmin - zmmin)/dz))
        iondensitygrid3d%zmin = dz*iz
        ix = (ceiling((iondensitygrid3d%xmax - xmmin)/dx))
        iondensitygrid3d%xmax = dx*ix
        iy = (ceiling((iondensitygrid3d%ymax - ymmin)/dy))
        iondensitygrid3d%ymax = dy*iy
        iz = (ceiling((iondensitygrid3d%zmax - zmmin)/dz))
        iondensitygrid3d%zmax = dz*iz
c       --- Now dx and dz can be set.
        iondensitygrid3d%dx = dx
        iondensitygrid3d%dy = dy
        iondensitygrid3d%dz = dz
c       --- This will calculate nx and nz and make sure the input is OK.
        call setupgrid3dtype(iondensitygrid3d,check)
        if (.not. check) then
           call kaboom("ERROR: luseparticledensity is being used but the input data is inconsistent")
           return
        endif
        
c       --- Now, copy the rho data into the iondensitygrid3d grid
        ixbemin = iondensitygrid3d%xmin/dx
        ixbemax = iondensitygrid3d%xmax/dx
        iybemin = iondensitygrid3d%ymin/dy
        iybemax = iondensitygrid3d%ymax/dy
        izbemin = iondensitygrid3d%zmin/dz
        izbemax = iondensitygrid3d%zmax/dz
        iondensitygrid3d%grid = rho(ixbemin:ixbemax,iybemin:iybemax,izbemin:izbemax)
c       iondensitygrid3d%grid = iondensitygrid3d%grid/eps0
      endif

      return
      end
c=============================================================================
      subroutine multigridbe3dsolve(iwhich,nx,ny,nz,dx,dy,dz,phi,rho,
     &                              rstar,linbend,bounds,
     &                              xmmin,ymmin,zmmin,zbeam,zgrid,
     &                              mgparam,mgiters,mgmaxiters,
     &                              mgmaxlevels,mgerror,mgtol,mgverbose,
     &                              downpasses,uppasses,
     &                              lcndbndy,laddconductor,icndbndy,
     &                              gridmode,conductors,iondensitygrid3d,
     &                              fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      use Grid3dtypemodule
      use Decompositionmodule
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels,mgverbose
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      type(Grid3dtype):: iondensitygrid3d
      type(Decomposition):: fsdecomp

c Use the multigrid method for solving Poisson's equation on a 3-D Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution and includes the Boltzmann electron term.
c
c When the grid cells are rectangular, semi-coarsening is done until the
c grid cell dimensions are roughly equal. Roughly equal means that
c   2/3 dx < dz < 4/3 dx
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3. Currently, it is still
c assumed that dx ~ dy and that semi-coarsening is not needed transversely.

      integer(ISZ):: i,ii,k,ix,iy,iz
      real(kind=8):: rs,x,r
      integer(ISZ),allocatable:: regionid(:,:,:)
      real(kind=8),allocatable:: phisave(:,:,:)
      real(kind=8):: bendx((nx+1)*(ny+1))
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- The parallel version does not yet work, it is just too complicated
c     --- and not needed yet.
      if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
        print*,"multigridbe3dsolve: does not yet work in parallel"
        call kaboom("multigridbe3dsolve: does not yet work in parallel")
        return
      endif

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

c     --- Determine the points that make up the conductor.  This takes extra
c     --- time and so should not be done if the grid is not moving in the lab
c     --- frame.  Set gridmode to 1 to avoid this call. The data is then
c     --- converted and expanded for the multigrid solver.
      if (gridmode == 0 .or. iwhich == -2) then
        conductors%interior%n = 0
        conductors%evensubgrid%n = 0
        conductors%oddsubgrid%n = 0
        if (laddconductor) call callpythonfunc("calladdconductor","controllers")
      endif
      call checkconductors(nx,ny,nz,nx,ny,nz,dx,dy,dz,conductors,fsdecomp)

!$OMP PARALLEL
!$OMP&PRIVATE(ii,i,k,rs,x,r,ix,iy,iz)

c     --- Prepare rho by dividing it by -eps0
      rho = -rho/eps0

      allocate(phisave(0:nx,0:ny,0:nz),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,ny,nz
        call kaboom("multigrid3dsolve: allocation error")
        return
      endif

c     --- Create and setup array which holds the BE region ids
      allocate(regionid(0:nx,0:ny,0:nz))
      call setupregionidsbe3d(regionid,nx,ny,nz,dx,dy,dz,
     &                        xmmin,ymmin,zmmin)

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- very few iterations are done.
      mgiters = 0
      mgerror = 2.*mgtol + 1.
      do while (mgerror > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phi(0:nx,0:ny,0:nz)

c       --- Do one vcycle.
        call vcyclebe(0,nx,ny,nz,dx,dy,dz,phi,rho,regionid,
     &                xmmin,ymmin,zmmin,
     &                rstar,linbend,bendx,bounds,mgparam,mgmaxlevels,
     &                downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &                iondensitygrid3d)

c       --- Calculate the change in phi.

c       --- This line seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop below.
c       mgerror = maxval(abs(phisave - phi))

        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nz
          do iy=0,ny
            do ix=0,nx
              mgerror = max(mgerror,abs(phisave(ix,iy,iz) - phi(ix,iy,iz)))
            enddo
          enddo
        enddo
!$OMP END DO

      enddo

c     --- Make a print out.
      if (mgverbose>=1 .or. mgerror > mgtol) then
        if (mgerror > mgtol) then
          call remark("MultigridBE: Maximum number of iterations reached")
        endif
        write(errline,20) mgerror,mgiters
  20    format("MultigridBE: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif
      
      deallocate(phisave)
      deallocate(regionid)

c     --- Undo the change of rho
      rho = -rho*eps0

!$OMP END PARALLEL

      if (lf3dtimesubs) timemultigrid3dsolve = timemultigrid3dsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcyclebe(mglevel,nx,ny,nz,dx,dy,dz,
     &                              phi,rho,regionid,xmmin,ymmin,zmmin,
     &                              rstar,linbend,bendx,
     $                              bounds,mgparam,
     &                              mgmaxlevels,downpasses,uppasses,
     &                              lcndbndy,icndbndy,conductors,
     &                              iondensitygrid3d)
      use ConductorTypemodule
c     use BoltzmannElectrons
      use Multigrid3d_diagnostic
      use Grid3dtypemodule
      integer(ISZ):: mglevel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1),rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      type(Grid3dtype):: iondensitygrid3d
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8),allocatable:: phicoarse(:,:,:),rhocoarse(:,:,:)
      real(kind=8),allocatable:: Lphi(:,:,:),Lphicoarse(:,:,:)
      real(kind=8),allocatable:: phicoarsesave(:,:,:)
      integer(ISZ),allocatable:: regionidcoarse(:,:,:)
      integer(ISZ):: i,iszone=1
      real(kind=8):: ff
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dycoarsesqi,dzcoarsesqi
      integer(ISZ):: allocerror

      dxsqi = 1./dx**2
      dysqi = 1./dy**2
      dzsqi = 1./dz**2

c     --- Do initial relaxation iterations
      do i=1,downpasses
        call relaxbe3d(mglevel,nx,ny,nz,phi,rho,regionid,rstar,
     &                 dxsqi,dysqi,dzsqi,dx,dy,dz,linbend,bendx,
     &                 bounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,iondensitygrid3d)
      enddo

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. ny >= 4 .and. nz >= 4 .and.
     &    mglevel < mgmaxlevels) then

        call getnextcoarselevel3d(nx,ny,nz,nx,ny,nz,dx,dy,dz,
     &                            nxcoarse,nycoarse,nzcoarse,
     &                            dxcoarse,dycoarse,dzcoarse)

        dxcoarsesqi = 1./dxcoarse**2
        dycoarsesqi = 1./dycoarse**2
        dzcoarsesqi = 1./dzcoarse**2

c       --- Alloate new work space
        allocate(phicoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          call kaboom("vcyclebe: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxcoarse,0:nycoarse,0:nzcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          call kaboom("vcyclebe: allocation error")
          return
        endif
        allocate(Lphi(-1:nx+1,-1:ny+1,-1:nz+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate Lphi to shape ",nx,ny,nz
          call kaboom("vcyclebe: allocation error")
          return
        endif
        allocate(Lphicoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate Lphicoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          call kaboom("vcyclebe: allocation error")
          return
        endif

c     --- Create and setup array which holds the BE region ids
        allocate(regionidcoarse(0:nxcoarse,0:nycoarse,0:nzcoarse),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate regionidcoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          call kaboom("vcyclebe: allocation error")
          return
        endif
        call setupregionidsbe3d(regionidcoarse,nxcoarse,nycoarse,nzcoarse,
     &                          dxcoarse,dycoarse,dzcoarse,
     &                          xmmin,ymmin,zmmin)

c       --- Calculate the coarsened phi
        phicoarse = 0.
        call restrictbe3d(nx,ny,nz,phi,1,1,1,
     &                    nxcoarse,nycoarse,nzcoarse,phicoarse,1,1,1,
     &                    bounds)
        call cond_potmgbe(conductors%interior,nxcoarse,nycoarse,nzcoarse,
     &                    phicoarse,mglevel+1)
        call applytransversebcbe3d(nxcoarse,nycoarse,nzcoarse,phicoarse,bounds)
        call applylongitudinalbcbe3d(nxcoarse,nycoarse,nzcoarse,
     &                               phicoarse,bounds)


c       --- Calculate the coarsened Lphi, putting it into rhocoarse
        call Lphibe(nx,ny,nz,dxsqi,dysqi,dzsqi,dx,dy,dz,phi,Lphi,regionid,
     &              mglevel,bounds,
     &              lcndbndy,icndbndy,conductors,iondensitygrid3d)
        Lphi(0:nx,0:ny,0:nz) = rho - Lphi(0:nx,0:ny,0:nz)
        call applytransversebcbe3d(nx,ny,nz,Lphi,bounds)
        call applylongitudinalbcbe3d(nx,ny,nz,Lphi,bounds)
        call restrictbe3d(nx,ny,nz,Lphi,1,1,1,
     &                    nxcoarse,nycoarse,nzcoarse,rhocoarse,0,0,0,
     &                    bounds)

        deallocate(Lphi)

c       --- Calculate L(R phi), adding it into rhocoarse
        call Lphibe(nxcoarse,nycoarse,nzcoarse,
     &              dxcoarsesqi,dycoarsesqi,dzcoarsesqi,
     &              dxcoarse,dycoarse,dzcoarse,
     &              phicoarse,Lphicoarse,regionidcoarse,
     &              mglevel+1,bounds,
     &              lcndbndy,icndbndy,conductors,iondensitygrid3d)
        rhocoarse = rhocoarse + Lphicoarse(0:nxcoarse,0:nycoarse,0:nzcoarse)
        deallocate(Lphicoarse)
        call cond_potmgbezero(conductors%interior,
     &                        nxcoarse,nycoarse,nzcoarse,
     &                        rhocoarse,mglevel+1,0,0)

c       --- Save the current coarsened phi since it is needed to apply
c       --- the corrections after relaxations.
        allocate(phicoarsesave(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate phicoarsesave to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          call kaboom("vcyclebe: allocation error")
          return
        endif
        phicoarsesave = phicoarse

c       --- Continue at the next coarsest level.
        call vcyclebe(mglevel+iszone,nxcoarse,nycoarse,nzcoarse,
     &                dxcoarse,dycoarse,dzcoarse,phicoarse,rhocoarse,
     &                regionidcoarse,xmmin,ymmin,zmmin,
     &                rstar,linbend,bendx,bounds,mgparam,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors,iondensitygrid3d)

c       --- Add in the correction term.
        phicoarse = phicoarse - phicoarsesave
        if (lprintmgphimaxchange) then
          print*,"Max change in phi = ",maxval(abs(phicoarse))," at MG level ",mglevel
        endif
        deallocate(phicoarsesave)
        call expandbe3d(nx,ny,nz,phi,
     &                  nxcoarse,nycoarse,nzcoarse,phicoarse,
     &                  bounds)
        call clampphitophimaxbe3d(nx,ny,nz,dx,dy,dz,phi,rho,regionid,
     &                            iondensitygrid3d)

        deallocate(phicoarse,rhocoarse)
        deallocate(regionidcoarse)

      endif

c     --- Do final relaxation passes.
      do i=1,uppasses
        call relaxbe3d(mglevel,nx,ny,nz,phi,rho,regionid,rstar,
     &                 dxsqi,dysqi,dzsqi,dx,dy,dz,linbend,bendx,
     &                 bounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,iondensitygrid3d)
        call clampphitophimaxbe3d(nx,ny,nz,dx,dy,dz,phi,rho,regionid,
     &                            iondensitygrid3d)
      enddo

      return
      end
c=============================================================================
      subroutine restrictbe3d(nx,ny,nz,u,delx,dely,delz,
     &                        nxcoarse,nycoarse,nzcoarse,ucoarse,
     &                        delcx,delcy,delcz,
     &                        bounds)

      integer(ISZ):: nx,ny,nz,delx,dely,delz,delcx,delcy,delcz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz)
      real(kind=8):: ucoarse(-delcx:nxcoarse+delcx,-delcy:nycoarse+delcy,-delcz:nzcoarse+delcz)
      integer(ISZ):: bounds(0:5)
      
c Restrict to a coarser grid.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixcoarse,iycoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      real(kind=8):: r,w,dx,dy,dz,dxi,dyi,dzi,wx(0:3),wy(0:3),wz(0:3)

c     --- Set the loop limits, always including edges.
      dx = 1.*nx/nxcoarse
      dxi = 1.*nxcoarse/nx
      if (ny > 0) then
        dy = 1.*ny/nycoarse
        dyi = 1.*nycoarse/ny
      else
        dy = 1.
        dyi = 1.
      endif
      dz = 1.*nz/nzcoarse
      dzi = 1.*nzcoarse/nz

c     --- Do the loops.
!$OMP DO
      do izcoarse=0,nzcoarse
        izmin = ((izcoarse-1)*nz + 4*nzcoarse)/nzcoarse-3
        izmax = ((izcoarse+1)*nz - 1)/nzcoarse
        if (izmin < -delz) izmin = -delz
        if (izmax > nz+delz) izmax = nz+delz

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - iz*dzi)
        enddo

        if (izcoarse == 0 .and. bounds(4) == 0) then
          izmin = 0
          izmax = 0
          wz(0) = 2.
        else if (izcoarse == nzcoarse .and. bounds(5) == 0) then
          izmin = nz
          izmax = nz
          wz(0) = 2.
        endif

        do iycoarse=0,nycoarse
          iymin = int(ceiling((iycoarse-1)*dy + 1.e-10))
          iymax = int(floor((iycoarse+1)*dy - 1.e-10))
          if (iymin < -dely) iymin = -dely
          if (iymax > ny+dely) iymax = ny+dely

          do iy=iymin,iymax
            wy(iy-iymin) = 1. - abs(iycoarse - iy*dyi)
          enddo

          if (iycoarse == 0 .and. bounds(2) == 0) then
            iymin = 0
            iymax = 0
            wy(0) = 2.
          else if (iycoarse == nycoarse .and. bounds(3) == 0) then
            iymin = ny
            iymax = ny
            wy(0) = 2.
          endif

          do ixcoarse=0,nxcoarse
            ixmin = int(ceiling((ixcoarse-1)*dx + 1.e-10))
            ixmax = int(floor((ixcoarse+1)*dx - 1.e-10))
            if (ixmin < -delx) ixmin = -delx
            if (ixmax > nx+delx) ixmax = nx+delx

            do ix=ixmin,ixmax
              wx(ix-ixmin) = 1. - abs(ixcoarse - ix*dxi)
            enddo

            if (ixcoarse == 0 .and. bounds(0) == 0) then
              ixmin = 0
              ixmax = 0
              wx(0) = 2.
            else if (ixcoarse == nxcoarse .and. bounds(1) == 0) then
              ixmin = nx
              ixmax = nx
              wx(0) = 2.
            endif

            r = 0.
            w = 0.
            do iz=izmin,izmax
              do iy=iymin,iymax
                do ix=ixmin,ixmax
                  r = r + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)*u(ix,iy,iz)
                  w = w + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)
                enddo
              enddo
            enddo
            if (w > 0.) then
              ucoarse(ixcoarse,iycoarse,izcoarse) = r/w
            else
              ucoarse(ixcoarse,iycoarse,izcoarse) = 0.
            endif

          enddo
        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine expandbe3d(nx,ny,nz,phi,
     &                      nxcoarse,nycoarse,nzcoarse,phicoarse,bounds)
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: phicoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1)
      integer(ISZ):: bounds(0:5)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: jx,jy,jz
      real(kind=8):: dx,dy,dz
      real(kind=8):: wx,wy,wz

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      iymin = 0
      iymax = ny
      izmin = 0
      izmax = nz
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx - 1
      if (bounds(2) == 0) iymin = 1
      if (bounds(3) == 0) iymax = ny - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz - 1

      dx = 1.*nxcoarse/nx
      dy = 1.*nycoarse/ny
      dz = 1.*nzcoarse/nz

!$OMP DO
      do iz=izmin,izmax
        jz = int((iz*nzcoarse)/nz)
        wz =  1.*(iz*nzcoarse)/nz - jz
        do iy=iymin,iymax
          jy = int(iy*dy)
          wy =     iy*dy - jy
          do ix=ixmin,ixmax
            jx = int(ix*dx)
            wx =     ix*dx - jx

            phi(ix,iy,iz) = phi(ix,iy,iz) +
     &             (1.-wx)*(1.-wy)*(1.-wz)*phicoarse(jx  ,jy  ,jz  ) +
     &                 wx *(1.-wy)*(1.-wz)*phicoarse(jx+1,jy  ,jz  ) +
     &             (1.-wx)*    wy *(1.-wz)*phicoarse(jx  ,jy+1,jz  ) +
     &                 wx *    wy *(1.-wz)*phicoarse(jx+1,jy+1,jz  ) +
     &             (1.-wx)*(1.-wy)*    wz *phicoarse(jx  ,jy  ,jz+1) +
     &                 wx *(1.-wy)*    wz *phicoarse(jx+1,jy  ,jz+1) +
     &             (1.-wx)*    wy *    wz *phicoarse(jx  ,jy+1,jz+1) +
     &                 wx *    wy *    wz *phicoarse(jx+1,jy+1,jz+1)
          enddo
        enddo
      enddo

      call applytransversebcbe3d(nx,ny,nz,phi,bounds)
      call applylongitudinalbcbe3d(nx,ny,nz,phi,bounds)

      return
      end
c=============================================================================
      subroutine relaxbe3d(mglevel,nx,ny,nz,phi,rho,regionid,rstar,
     &                     dxsqi,dysqi,dzsqi,dx,dy,dz,linbend,bendx,bounds,
     &                     mgparam,lcndbndy,icndbndy,conductors,
     &                     iondensitygrid3d)
      use Constant
      use ConductorTypemodule
      use BoltzmannElectrons
      use Grid3dtypemodule
      integer(ISZ):: mglevel,nx,ny,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1),rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      real(kind=8):: dxsqi,dysqi,dzsqi,dx,dy,dz
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      type(Grid3dtype):: iondensitygrid3d

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.  It makes calls to the routines which specify internal
c conductors. The routine also allows for a bent beam-pipe.
c
c The tranverse boundaries can either be held constant, have zero normal
c derivative, or be periodic.  When BOUNDXY is zero, the boundaries are held
c constant, when 1, they have zero normal derivative, and when 2, the
c boundaries are periodic.
c
c The longitudinal boundaries can either be held constant, have zero normal
c derivative, or be periodic.  When BOUND0 or BOUNDNZ is zero, the boundaries
c are held constant, when 1, they have zero normal derivative, and when 2, the
c boundaries are periodic.
c
c Note that loops over all directions assume that nx and ny are even.
c
c The arrangement of the loops was done to increase performance.  The entire
c grid is looped over as if it were a 1D array, ignoring boundaries.
c The boundaries are then reset, the previous value was destroyed.
c
c rstar(-1) and rstar(nz+1) are set based on the axial boundary conditions.

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ix,iy,iz,id
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ix1
      integer(ISZ):: i1,i2,ic
      real(kind=8):: iondensitytmp,rhoe,denom,Lphi,expo

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgbe(conductors%interior,nx,ny,nz,phi,mglevel)
      call condbndymgintbe3d(conductors%evensubgrid,nx,ny,nz,phi,bounds,mglevel)
      call condbndymgintbe3d(conductors%oddsubgrid ,nx,ny,nz,phi,bounds,mglevel)


c     --- Set starting and ending parity.
      s_parity = 0
      e_parity = 1

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      iymin = 0
      iymax = ny
      izmin = 0
      izmax = nz
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx - 1
      if (bounds(2) == 0) iymin = 1
      if (bounds(3) == 0) iymax = ny - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz - 1

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

c       --- Save values just outside conductor surfaces. Only save phi at the
c       --- subgrid points which are to be used at the current level of
c       --- grid refinement.
        if (lcndbndy) then
          if (parity == 0) then
            i1 = conductors%evensubgrid%istart(mglevel)
            i2 = conductors%evensubgrid%istart(mglevel+1)-1
            do ic = i1,i2
              ix = conductors%evensubgrid%indx(0,ic)
              iy = conductors%evensubgrid%indx(1,ic)
              iz = conductors%evensubgrid%indx(2,ic)
              conductors%evensubgrid%prevphi(ic) = phi(ix,iy,iz)
            enddo
          else
            i1 = conductors%oddsubgrid%istart(mglevel)
            i2 = conductors%oddsubgrid%istart(mglevel+1)-1
            do ic = i1,i2
              ix = conductors%oddsubgrid%indx(0,ic)
              iy = conductors%oddsubgrid%indx(1,ic)
              iz = conductors%oddsubgrid%indx(2,ic)
              conductors%oddsubgrid%prevphi(ic) = phi(ix,iy,iz)
            enddo
          endif
        endif

!$OMP DO
        do iz=izmin,izmax

          do iy=iymin,iymax

            ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
            do ix=ix1,ixmax,2

              rhoe = 0.
              denom = -2.*(dxsqi+dysqi+dzsqi)
              if (regionid(ix,iy,iz) > 0) then
                id = regionid(ix,iy,iz)
                call getlocaldens3d(ix,iy,iz,dx,dy,dz,id,iondensitygrid3d,
     &                              iondensitytmp)
                if (iondensitytmp .ne. 0. .and.
     &              electrontemperature(id) .ne. 0) then
                  expo = (phi(ix,iy,iz) - plasmapotential(id))/electrontemperature(id)
                  expo = min(expo,log(electrondensitymaxscale(id)))
                  rhoe = iondensitytmp*exp(expo)
                  denom = denom - rhoe/eps0/electrontemperature(id)
                endif
              endif
              Lphi =  (phi(ix-1,iy  ,iz  )+phi(ix+1,iy  ,iz  ))*dxsqi
     &             +  (phi(ix  ,iy-1,iz  )+phi(ix  ,iy+1,iz  ))*dysqi
     &             +  (phi(ix  ,iy  ,iz-1)+phi(ix  ,iy  ,iz+1))*dzsqi
     &              -  phi(ix,iy,iz)*2.*(dxsqi+dysqi+dzsqi) - rhoe/eps0
              phi(ix,iy,iz) = phi(ix,iy,iz) - mgparam*(Lphi - rho(ix,iy,iz))/denom
            enddo
          enddo
        enddo
!$OMP END DO

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
          if (parity == 0) then
           call condbndymgbe(conductors%evensubgrid,nx,ny,nz,phi,rho,regionid,
     &                       dxsqi,dysqi,dzsqi,dx,dy,dz,mgparam,bounds,
     &                       mglevel,icndbndy,iondensitygrid3d)
          endif
          if (parity == 1) then
           call condbndymgbe(conductors%oddsubgrid,nx,ny,nz,phi,rho,regionid,
     &                       dxsqi,dysqi,dzsqi,dx,dy,dz,mgparam,bounds,
     &                       mglevel,icndbndy,iondensitygrid3d)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgbe(conductors%interior,nx,ny,nz,phi,mglevel)
        call condbndymgintbe3d(conductors%evensubgrid,nx,ny,nz,phi,bounds,mglevel)
        call condbndymgintbe3d(conductors%oddsubgrid ,nx,ny,nz,phi,bounds,mglevel)

c       --- set phi in the guard planes
c       --- This must be done inside the loop over parities so that the
c       --- guard planes are updated with the most recent values.
        call applytransversebcbe3d(nx,ny,nz,phi,bounds)
        call applylongitudinalbcbe3d(nx,ny,nz,phi,bounds)

c     --- end of loop over even and odd points
      enddo

      return
      end
c=============================================================================
      subroutine cond_potmgbe(interior,nx,ny,nz,phi,mglevel)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)

c Set conductor points to the desired potential. The potential is used since
c at all levels, phi is being operated on directly.

      integer(ISZ):: ic,ix,iy,iz

!$OMP DO
      do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        phi(ix,iy,iz) = interior%volt(ic)
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine cond_potmgbezero(interior,nx,ny,nz,u,mglevel,delt,delz)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,ny,nz,mglevel,delt,delz
      real(kind=8):: u(-delt:nx+delt,-delt:ny+delt,-delz:nz+delz)

c Set data at conductor points to zero.

      integer(ISZ):: ic,ix,iy,iz

!$OMP DO
      do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        u(ix,iy,iz) = 0.
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine condbndymgbe(subgrid,nx,ny,nz,phi,rho,regionid,
     &                        dxsqi,dysqi,dzsqi,dx,dy,dz,
     &                        mgparam,bounds,mglevel,icndbndy,iondensitygrid3d)
      use Constant
      use ConductorSubGridTypemodule
      use BoltzmannElectrons
      use Grid3dtypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1), rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: dxsqi,dysqi,dzsqi,dx,dy,dz,mgparam
      integer(ISZ):: bounds(0:5),icndbndy
      type(Grid3dtype):: iondensitygrid3d

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold
c phi(i-+1)-phi(i) at minus and plus one in each direction.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c The Cx, Cy, and Cz hold the numerator of the coefficients of phi(i,j,k).
c The delx, dely, and delz hold the denominator of the coefficients of the
c full finite difference of phi.
c For icndbndy==1, these coefficients are just 1. For icndbndy==2, they
c include the dels.
 
      real(kind=8):: pijk,pxm,pym,pzm,pxp,pyp,pzp,denom,rhoe,Lphi,expo
      real(kind=8):: delx,dely,delz,Cx,Cy,Cz,ppp
      real(kind=8):: delxsqi,delysqi,delzsqi
      real(kind=8):: iondensitytmp
      integer(ISZ):: ic,ix,iy,iz,id
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

c       --- Set temporaries with initial values.
        pijk = subgrid%prevphi(ic)
        pxm = phi(ix-1,iy  ,iz  ) - pijk
        pxp = phi(ix+1,iy  ,iz  ) - pijk
        pym = phi(ix  ,iy-1,iz  ) - pijk
        pyp = phi(ix  ,iy+1,iz  ) - pijk
        pzm = phi(ix  ,iy  ,iz-1) - pijk
        pzp = phi(ix  ,iy  ,iz+1) - pijk
        delx = 1.
        dely = 1.
        delz = 1.
        Cx = 2.
        Cy = 2.
        Cz = 2.
        ppp = 1.

c       --- the point lower in x is inside the conductor
        if (0 < dels(0,ic) .and. dels(0,ic) < 1.) then
          pxm = (volt(0,ic) - pijk)/dels(0,ic)
          Cx = Cx - 1. + 1./dels(0,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(0,ic)
          ppp = min(ppp,dels(0,ic))
        else if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          pxm = 0.
          Cx = Cx - 1.
          delx = delx - 0.5 + (-dels(0,ic))
          if (-dels(0,ic) > 0.) then
            ppp = min(ppp,-dels(0,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point higher in x is inside the conductor
        if (0 < dels(1,ic) .and. dels(1,ic) < 1.) then
          pxp = (volt(1,ic) - pijk)/dels(1,ic)
          Cx = Cx - 1. + 1./dels(1,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(1,ic)
          ppp = min(ppp,dels(1,ic))
        else if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          pxp = 0.
          Cx = Cx - 1.
          delx = delx - 0.5 + (-dels(1,ic))
          if (-dels(1,ic) > 0.) then
            ppp = min(ppp,-dels(1,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point lower in y is inside the conductor
        if (0 < dels(2,ic) .and. dels(2,ic) < 1.) then
          pym = (volt(2,ic) - pijk)/dels(2,ic)
          Cy = Cy - 1. + 1./dels(2,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(2,ic)
          ppp = min(ppp,dels(2,ic))
        else if (-1. < dels(2,ic) .and. dels(2,ic) <= 0.) then
          pym = 0.
          Cy = Cy - 1.
          dely = dely - 0.5 + (-dels(2,ic))
          if (-dels(2,ic) > 0.) then
            ppp = min(ppp,-dels(2,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point higher in y is inside the conductor
        if (0 < dels(3,ic) .and. dels(3,ic) < 1.) then
          pyp = (volt(3,ic) - pijk)/dels(3,ic)
          Cy = Cy - 1. + 1./dels(3,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(3,ic)
          ppp = min(ppp,dels(3,ic))
        else if (-1. < dels(3,ic) .and. dels(3,ic) <= 0.) then
          pyp = 0.
          Cy = Cy - 1.
          dely = dely - 0.5 + (-dels(3,ic))
          if (-dels(3,ic) > 0.) then
            ppp = min(ppp,-dels(3,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point lower in z is inside the conductor
        if (0 < dels(4,ic) .and. dels(4,ic) < 1.) then
          pzm = (volt(4,ic) - pijk)/dels(4,ic)
          Cz = Cz - 1. + 1./dels(4,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(4,ic)
          ppp = min(ppp,dels(4,ic))
        else if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          pzm = 0.
          Cz = Cz - 1.
          delz = delz - 0.5 + (-dels(4,ic))
          if (-dels(4,ic) > 0.) then
            ppp = min(ppp,-dels(4,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point higher in z is inside the conductor
        if (0 < dels(5,ic) .and. dels(5,ic) < 1.) then
          pzp = (volt(5,ic) - pijk)/dels(5,ic)
          Cz = Cz - 1. + 1./dels(5,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(5,ic)
          ppp = min(ppp,dels(5,ic))
        else if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          pzp = 0.
          Cz = Cz - 1.
          delz = delz - 0.5 + (-dels(5,ic))
          if (-dels(5,ic) > 0.) then
            ppp = min(ppp,-dels(5,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
        if (mglevel == 0) ppp = 1.
c       --- calculate the new phi based on the boundary conditions
        delxsqi = dxsqi/dvnz(delx)
        delysqi = dysqi/dvnz(dely)
        delzsqi = dzsqi/dvnz(delz)
        rhoe = 0.
        denom = - Cx*delxsqi - Cy*delysqi - Cz*delzsqi
        if (regionid(ix,iy,iz) > 0) then
          id = regionid(ix,iy,iz)
          call getlocaldens3d(ix,iy,iz,dx,dy,dz,id,iondensitygrid3d,
     &                        iondensitytmp)
          if (iondensitytmp .ne. 0. .and. electrontemperature(id) .ne. 0) then
            expo = (pijk - plasmapotential(id))/electrontemperature(id)
            expo = min(expo,log(electrondensitymaxscale(id)))
            rhoe = iondensitytmp*exp(expo)
            denom = denom - rhoe/eps0/electrontemperature(id)
          endif
        endif
        Lphi = (pxm+pxp)*delxsqi + (pym+pyp)*delysqi + (pzm+pzp)*delzsqi
     &         -rhoe/eps0
        phi(ix,iy,iz) = pijk - mgparam*(Lphi*ppp - rho(ix,iy,iz))/denom
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine condbndymgintbe3d(subgrid,nx,ny,nz,phi,bounds,mglevel)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: bounds(0:5)

c Sets the potential on points just beyond a Neumann boundary condition.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

        if (0. > dels(0,ic) .and. dels(0,ic) > -1.) phi(ix-1,iy,iz) = phi(ix,iy,iz)
        if (0. > dels(1,ic) .and. dels(1,ic) > -1.) phi(ix+1,iy,iz) = phi(ix,iy,iz)
        if (0. > dels(2,ic) .and. dels(2,ic) > -1.) phi(ix,iy-1,iz) = phi(ix,iy,iz)
        if (0. > dels(3,ic) .and. dels(3,ic) > -1.) phi(ix,iy+1,iz) = phi(ix,iy,iz)
        if (0. > dels(4,ic) .and. dels(4,ic) > -1.) phi(ix,iy,iz-1) = phi(ix,iy,iz)
        if (0. > dels(5,ic) .and. dels(5,ic) > -1.) phi(ix,iy,iz+1) = phi(ix,iy,iz)

        if (0. == dels(0,ic)) phi(ix+1,iy,iz) = phi(ix,iy,iz)
        if (0. == dels(1,ic)) phi(ix-1,iy,iz) = phi(ix,iy,iz)
        if (0. == dels(2,ic)) phi(ix,iy+1,iz) = phi(ix,iy,iz)
        if (0. == dels(3,ic)) phi(ix,iy-1,iz) = phi(ix,iy,iz)
        if (0. == dels(4,ic)) phi(ix,iy,iz+1) = phi(ix,iy,iz)
        if (0. == dels(5,ic)) phi(ix,iy,iz-1) = phi(ix,iy,iz)

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine Lphibe(nx,ny,nz,dxsqi,dysqi,dzsqi,dx,dy,dz,phi,Lphi,regionid,
     &                  mglevel,bounds,
     &                  lcndbndy,icndbndy,conductors,iondensitygrid3d)
      use Constant
      use ConductorTypemodule
      use BoltzmannElectrons
      use Grid3dtypemodule
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dxsqi,dysqi,dzsqi,dx,dy,dz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: Lphi(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      integer(ISZ):: mglevel,bounds(0:5)
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      type(Grid3dtype):: iondensitygrid3d

c Calculate the Lphi on the grid including the BE term.
c Note that what is calculated here is not the residual since it does not
c include the source term. This will not go to zero when convergence
c is reached.

      real(kind=8):: rhoe,expo,iondensitytmp
      integer(ISZ):: ix,iy,iz,ic,id
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      iymin = 0
      iymax = ny
      izmin = 0
      izmax = nz
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx-1
      if (bounds(2) == 0) iymin = 1
      if (bounds(3) == 0) iymax = ny-1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

      Lphi = 0.
c     --- Calculate the Lphi.
!$OMP DO
      do iz=izmin,izmax

        do iy=iymin,iymax

          do ix=ixmin,ixmax

            rhoe = 0.
            if (regionid(ix,iy,iz) > 0) then
              id = regionid(ix,iy,iz)
              call getlocaldens3d(ix,iy,iz,dx,dy,dz,id,iondensitygrid3d,
     &                            iondensitytmp)
              if (iondensitytmp .ne. 0. .and.
     &            electrontemperature(id) .ne. 0) then
                expo = (phi(ix,iy,iz) - plasmapotential(id))/electrontemperature(id)
                expo = min(expo,log(electrondensitymaxscale(id)))
                rhoe = iondensitytmp*exp(expo)
              endif
            endif
            Lphi(ix,iy,iz) = (phi(ix-1,iy  ,iz  )+phi(ix+1,iy  ,iz  ))*dxsqi
     &                    +  (phi(ix  ,iy-1,iz  )+phi(ix  ,iy+1,iz  ))*dysqi
     &                    +  (phi(ix  ,iy  ,iz-1)+phi(ix  ,iy  ,iz+1))*dzsqi
     &                    -  phi(ix,iy,iz)*2.*(dxsqi+dysqi+dzsqi)
     &                    - rhoe/eps0
          enddo
        enddo
      enddo
!$OMP END DO

c     --- Zero the Lphi inside conductors.
      call cond_potmgbezero(conductors%interior,nx,ny,nz,Lphi,mglevel,1,1)

      if (lcndbndy) then
c       --- Calculate the Lphi near the conductor.
        call condbndyLphibe(conductors%evensubgrid,nx,ny,nz,phi,Lphi,regionid,
     &                     dxsqi,dysqi,dzsqi,dx,dy,dz,bounds,mglevel,icndbndy,
     &                     iondensitygrid3d)
        call condbndyLphibe(conductors%oddsubgrid,nx,ny,nz,phi,Lphi,regionid,
     &                     dxsqi,dysqi,dzsqi,dx,dy,dz,bounds,mglevel,icndbndy,
     &                     iondensitygrid3d)
      endif

c     --- Transverse boundaries
      if (bounds(0) == 0) Lphi(-1:0,:,:) = 0.
      if (bounds(1) == 0) Lphi(nx:nx+1,:,:) = 0.
      if (bounds(0) == 1)
     &  Lphi(-1,iymin:iymax,izmin:izmax) = Lphi(1,iymin:iymax,izmin:izmax)
      if (bounds(1) == 1)
     &  Lphi(nx+1,iymin:iymax,izmin:izmax) = Lphi(nx-1,iymin:iymax,izmin:izmax)
      if (bounds(0) == 2)
     &  Lphi(-1,iymin:iymax,izmin:izmax) = Lphi(nx-1,iymin:iymax,izmin:izmax)
      if (bounds(1) == 2)
     &  Lphi(nx:nx+1,iymin:iymax,izmin:izmax) = Lphi(0:1,iymin:iymax,izmin:izmax)

      if (bounds(2) == 0) Lphi(:,-1:0,:) = 0.
      if (bounds(3) == 0) Lphi(:,ny:ny+1,:) = 0.
      if (bounds(2) == 1) Lphi(:,-1,izmin:izmax) = Lphi(:,1,izmin:izmax)
      if (bounds(3) == 1) Lphi(:,ny+1,izmin:izmax) = Lphi(:,ny-1,izmin:izmax)
      if (bounds(2) == 2) Lphi(:,-1,izmin:izmax) = Lphi(:,ny-1,izmin:izmax)
      if (bounds(3) == 2) Lphi(:,ny:ny+1,izmin:izmax) = Lphi(:,0:1,izmin:izmax)

c     --- Longitudinal boundaries
      if (bounds(4) == 0) Lphi(:,:,-1:0) = 0.
      if (bounds(5) == 0) Lphi(:,:,nz:nz+1) = 0.
      if (bounds(4) == 1) Lphi(:,:,-1) = Lphi(:,:,1)
      if (bounds(5) == 1) Lphi(:,:,nz+1) = Lphi(:,:,nz-1)
      if (bounds(4) == 2 .and. nz == nz) Lphi(:,:,-1) = Lphi(:,:,nz-1)
      if (bounds(5) == 2 .and. nz == nz) Lphi(:,:,nz+1) = Lphi(:,:,1)

      return
      end
c=============================================================================
      subroutine condbndyLphibe(subgrid,nx,ny,nz,phi,Lphi,regionid,
     &                          dxsqi,dysqi,dzsqi,dx,dy,dz,
     &                          bounds,mglevel,icndbndy,iondensitygrid3d)
      use Constant
      use ConductorSubGridTypemodule
      use BoltzmannElectrons
      use Grid3dtypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1), Lphi(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: dxsqi,dysqi,dzsqi,dx,dy,dz
      integer(ISZ):: bounds(0:5),icndbndy
      type(Grid3dtype):: iondensitygrid3d

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold
c phi(i-+1)-phi(i) at minus and plus one in each direction.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c The Cx, Cy, and Cz hold the numerator of the coefficients of phi(i,j,k).
c The delx, dely, and delz hold the denominator of the coefficients of the
c full finite difference of phi.
c For icndbndy==1, these coefficients are just 1. For icndbndy==2, they
c include the dels.
 
      real(kind=8):: pijk,pxm,pym,pzm,pxp,pyp,pzp,denom,rhoe,ppp,expo
      real(kind=8):: iondensitytmp
      real(kind=8):: delx,dely,delz
      real(kind=8):: delxsqi,delysqi,delzsqi
      integer(ISZ):: ic,ix,iy,iz,id
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

c       --- Set temporaries with initial values.
        pijk = phi(ix,iy,iz)
        pxm = phi(ix-1,iy  ,iz  ) - pijk
        pxp = phi(ix+1,iy  ,iz  ) - pijk
        pym = phi(ix  ,iy-1,iz  ) - pijk
        pyp = phi(ix  ,iy+1,iz  ) - pijk
        pzm = phi(ix  ,iy  ,iz-1) - pijk
        pzp = phi(ix  ,iy  ,iz+1) - pijk
        delx = 1.
        dely = 1.
        delz = 1.
        ppp = 1.

c       --- the point lower in x is inside the conductor
        if (0 < dels(0,ic) .and. dels(0,ic) < 1.) then
          pxm = (volt(0,ic) - pijk)/dels(0,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(0,ic)
          ppp = min(ppp,dels(0,ic))
        else if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          pxm = 0.
          delx = delx - 0.5 + (-dels(0,ic))
          if (-dels(0,ic) > 0.) then
            ppp = min(ppp,-dels(0,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point higher in x is inside the conductor
        if (0 < dels(1,ic) .and. dels(1,ic) < 1.) then
          pxp = (volt(1,ic) - pijk)/dels(1,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(1,ic)
          ppp = min(ppp,dels(1,ic))
        else if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          pxp = 0.
          delx = delx - 0.5 + (-dels(1,ic))
          if (-dels(1,ic) > 0.) then
            ppp = min(ppp,-dels(1,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point lower in y is inside the conductor
        if (0 < dels(2,ic) .and. dels(2,ic) < 1.) then
          pym = (volt(2,ic) - pijk)/dels(2,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(2,ic)
          ppp = min(ppp,dels(2,ic))
        else if (-1. < dels(2,ic) .and. dels(2,ic) <= 0.) then
          pym = 0.
          dely = dely - 0.5 + (-dels(2,ic))
          if (-dels(2,ic) > 0.) then
            ppp = min(ppp,-dels(2,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point higher in y is inside the conductor
        if (0 < dels(3,ic) .and. dels(3,ic) < 1.) then
          pyp = (volt(3,ic) - pijk)/dels(3,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(3,ic)
          ppp = min(ppp,dels(3,ic))
        else if (-1. < dels(3,ic) .and. dels(3,ic) <= 0.) then
          pyp = 0.
          dely = dely - 0.5 + (-dels(3,ic))
          if (-dels(3,ic) > 0.) then
            ppp = min(ppp,-dels(3,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point lower in z is inside the conductor
        if (0 < dels(4,ic) .and. dels(4,ic) < 1.) then
          pzm = (volt(4,ic) - pijk)/dels(4,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(4,ic)
          ppp = min(ppp,dels(4,ic))
        else if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          pzm = 0.
          delz = delz - 0.5 + (-dels(4,ic))
          if (-dels(4,ic) > 0.) then
            ppp = min(ppp,-dels(4,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- the point higher in z is inside the conductor
        if (0 < dels(5,ic) .and. dels(5,ic) < 1.) then
          pzp = (volt(5,ic) - pijk)/dels(5,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(5,ic)
          ppp = min(ppp,dels(5,ic))
        else if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          pzp = 0.
          delz = delz - 0.5 + (-dels(5,ic))
          if (-dels(5,ic) > 0.) then
            ppp = min(ppp,-dels(5,ic))
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
c       --- calculate the new phi based on the boundary conditions
        delxsqi = dxsqi/dvnz(delx)
        delysqi = dysqi/dvnz(dely)
        delzsqi = dzsqi/dvnz(delz)
        rhoe = 0.
        if (regionid(ix,iy,iz) > 0) then
          id = regionid(ix,iy,iz)
          call getlocaldens3d(ix,iy,iz,dx,dy,dz,id,iondensitygrid3d,
     &                        iondensitytmp)
          if (iondensitytmp .ne. 0. .and. electrontemperature(id) .ne. 0) then
            expo = (pijk - plasmapotential(id))/electrontemperature(id)
            expo = min(expo,log(electrondensitymaxscale(id)))
            rhoe = iondensitytmp*exp(expo)
          endif
        endif
        Lphi(ix,iy,iz) = (+ (pxm+pxp)*delxsqi
     &                    + (pym+pyp)*delysqi
     &                    + (pzm+pzp)*delzsqi
     &                    - rhoe/eps0)*ppp

c       --- Zero out the Lphi in points just beyond a Neumann boundary
        if (0. > dels(0,ic) .and. dels(0,ic) >= -1.) Lphi(ix-1,iy,iz) = 0.
        if (0. > dels(1,ic) .and. dels(1,ic) >= -1.) Lphi(ix+1,iy,iz) = 0.
        if (0. > dels(2,ic) .and. dels(2,ic) >= -1.) Lphi(ix,iy-1,iz) = 0.
        if (0. > dels(3,ic) .and. dels(3,ic) >= -1.) Lphi(ix,iy+1,iz) = 0.
        if (0. > dels(4,ic) .and. dels(4,ic) >= -1.) Lphi(ix,iy,iz-1) = 0.
        if (0. > dels(5,ic) .and. dels(5,ic) >= -1.) Lphi(ix,iy,iz+1) = 0.

        if (0. == dels(0,ic)) Lphi(ix+1,iy,iz) = 0.
        if (0. == dels(1,ic)) Lphi(ix-1,iy,iz) = 0.
        if (0. == dels(2,ic)) Lphi(ix,iy+1,iz) = 0.
        if (0. == dels(3,ic)) Lphi(ix,iy-1,iz) = 0.
        if (0. == dels(4,ic)) Lphi(ix,iy,iz+1) = 0.
        if (0. == dels(5,ic)) Lphi(ix,iy,iz-1) = 0.

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine applytransversebcbe3d(nx,ny,nz,u,bounds)
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:ny+1,-1:nz+1)

      if (bounds(0) == 0) u(-1,:,:)   = 2.*u(0,:,:) - u(1,:,:)
      if (bounds(1) == 0) u(nx+1,:,:) = 2.*u(nx,:,:) - u(nx-1,:,:)
      if (bounds(0) == 1) u(-1,:,:)   = u(1,:,:)
      if (bounds(1) == 1) u(nx+1,:,:) = u(nx-1,:,:)
      if (bounds(0) == 2) u(-1,:,:)   = u(nx-1,:,:)
      if (bounds(1) == 2) u(nx+1,:,:) = u(1,:,:)

      if (bounds(2) == 0) u(:,-1,:)   = 2.*u(:,0,:) - u(:,1,:)
      if (bounds(3) == 0) u(:,ny+1,:) = 2.*u(:,ny,:) - u(:,ny-1,:)
      if (bounds(2) == 1) u(:,-1,:)   = u(:,1,:)
      if (bounds(3) == 1) u(:,ny+1,:) = u(:,ny-1,:)
      if (bounds(2) == 2) u(:,-1,:)   = u(:,ny-1,:)
      if (bounds(3) == 2) u(:,ny+1,:) = u(:,1,:)

      return
      end
c=============================================================================
      subroutine applylongitudinalbcbe3d(nx,ny,nz,u,bounds)
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:ny+1,-1:nz+1)
      if (bounds(4) == 0) u(:,:,-1)   = 2.*u(:,:,0) - u(:,:,1)
      if (bounds(5) == 0) u(:,:,nz+1) = 2.*u(:,:,nz) - u(:,:,nz-1)
      if (bounds(4) == 1) u(:,:,-1)   = u(:,:,1)
      if (bounds(5) == 1) u(:,:,nz+1) = u(:,:,nz-1)
      if (bounds(4) == 2) u(:,:,-1)   = u(:,:,nz-1)
      if (bounds(5) == 2) u(:,:,nz+1) = u(:,:,1)

      return
      end
c=============================================================================
      subroutine clampphitophimaxbe3d(nx,ny,nz,dx,dy,dz,phi,rho,regionid,
     &                                iondensitygrid3d)
      use Constant
      use BoltzmannElectrons
      use Grid3dtypemodule
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1),rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      type(Grid3dtype):: iondensitygrid3d

      integer(ISZ):: ix,iy,iz,id
      real(kind=8):: expo,phimax,rhomax,iondensitytmp

c     --- Clamp phi to phimax. The phimax is calculated assuming that
c     --- the max electron density will be the max of the iondensity and
c     --- and the ion particle density (rho).
      do iz=0,nz
        do iy=0,ny
          do ix=0,nx
            if (regionid(ix,iy,iz) > 0) then
              id = regionid(ix,iy,iz)
              call getlocaldens3d(ix,iy,iz,dx,dy,dz,id,iondensitygrid3d,
     &                            iondensitytmp)
              if (iondensitytmp == 0.) cycle
              if (luseparticledensity(id) > 0) then
                rhomax = max(electrondensitymaxscale(id),
     &                       -rho(ix,iy,iz)*eps0/iondensitytmp)
              else
                rhomax = electrondensitymaxscale(id)
              endif
              phimax = plasmapotential(id)+electrontemperature(id)*log(rhomax)
              if (phi(ix,iy,iz) > phimax) phi(ix,iy,iz) = phimax
            endif
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine setupregionidsbe3d(regionid,nx,ny,nz,dx,dy,dz,
     &                              xmmin,ymmin,zmmin)
      use BoltzmannElectrons
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz,xmmin,ymmin,zmmin

c Sets id of regions where the Boltzmann-electron distribution is used.

      integer(ISZ):: id,ixmin,ixmax,iymin,iymax,izmin,izmax

      regionid = 0
      do id=1,nberegions
        ixmin = (max(0.,xbemin(id) - xmmin))/dx
        ixmax = (min(nx*dx,xbemax(id) - xmmin))/dx
        iymin = (max(0.,ybemin(id) - ymmin))/dy
        iymax = (min(ny*dy,ybemax(id) - ymmin))/dy
        izmin = (max(0.,zbemin(id) - zmmin))/dz
        izmax = (min(nz*dz,zbemax(id) - zmmin))/dz
        regionid(ixmin:ixmax,iymin:iymax,izmin:izmax) = id
      enddo

      return
      end
c=============================================================================
      subroutine getlocaldens3d(ix,iy,iz,dx,dy,dz,id,iondensitygrid3d,
     &                          iondensitytmp)
      use BoltzmannElectrons,Only: liondensitygrid3d,iondensity
      use Grid3dtypemodule
      integer(ISZ):: ix,iy,iz,id
      real(kind=8):: dx,dy,dz
      type(Grid3dtype):: iondensitygrid3d
      real(kind=8):: iondensitytmp

c     generate interpolated local iondensity

      real(kind=8):: reixx,reiyy,reizz
      real(kind=8):: xx,yy,zz,wxx,wyy,wzz
      integer(ISZ):: ixx,iyy,izz

      iondensitytmp = 0.
      if (liondensitygrid3d(id)) then
        xx = ix*dx
        yy = iy*dy
        zz = iz*dz
        reixx = (xx - iondensitygrid3d%xmin)/iondensitygrid3d%dx
        reiyy = (yy - iondensitygrid3d%ymin)/iondensitygrid3d%dy
        reizz = (zz - iondensitygrid3d%zmin)/iondensitygrid3d%dz
        ixx = reixx
        iyy = reiyy
        izz = reizz
        wxx = reixx-ixx
        wyy = reiyy-iyy
        wzz = reizz-izz
        if (ixx+1 > iondensitygrid3d%nx) then
          ixx = iondensitygrid3d%nx-1
          wxx = 1.
        endif
        if (iyy+1 > iondensitygrid3d%ny) then
          iyy = iondensitygrid3d%ny-1
          wyy = 1.
        endif
        if (izz+1 > iondensitygrid3d%nz) then
          izz = iondensitygrid3d%nz-1
          wzz = 1.
        endif
        iondensitytmp =
     &   iondensitygrid3d%grid(ixx  ,iyy  ,izz  )*(1.-wxx)*(1.-wyy)*(1.-wzz) +
     &   iondensitygrid3d%grid(ixx+1,iyy  ,izz  )*    wxx *(1.-wyy)*(1.-wzz) +
     &   iondensitygrid3d%grid(ixx  ,iyy+1,izz  )*(1.-wxx)*    wyy *(1.-wzz) +
     &   iondensitygrid3d%grid(ixx+1,iyy+1,izz  )*    wxx *    wyy *(1.-wzz) +
     &   iondensitygrid3d%grid(ixx  ,iyy  ,izz+1)*(1.-wxx)*(1.-wyy)*    wzz  +
     &   iondensitygrid3d%grid(ixx+1,iyy  ,izz+1)*    wxx *(1.-wyy)*    wzz  +
     &   iondensitygrid3d%grid(ixx  ,iyy+1,izz+1)*(1.-wxx)*    wyy *    wzz  +
     &   iondensitygrid3d%grid(ixx+1,iyy+1,izz+1)*    wxx *    wyy *    wzz

      else
        iondensitytmp = iondensity(id)
      endif

      return
      end
c=============================================================================
