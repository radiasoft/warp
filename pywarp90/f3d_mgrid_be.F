#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID_BE.M,
c version $Revision: 1.30 $, $Date: 2006/11/14 18:16:58 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid field sovler which is part of the F3D
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigridbe3df(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                          rstar,linbend,
     &                          bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                          xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid)
      use GlobalVars
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use BoltzmannElectrons
      use Parallel
c     use ifcore
c     integer:: fsf
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      integer(ISZ):: bound0,boundnz,boundxy
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid

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

      call multigridbe3dsolve(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                        rstar,linbend,bounds,
     &                        xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid,
     &                        mgparam,mgiters,mgmaxiters,
     &                        mgmaxlevels,mgerror,mgtol,
     &                        downpasses,uppasses,
     &                        lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                        gridmode,conductors,
     &                        my_index,nslaves,izfsslave,nzfsslave)

      return
      end
c=============================================================================
      subroutine multigridbe3dsolve(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                              rstar,linbend,bounds,
     &                              xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid,
     &                              mgparam,mgiters,mgmaxiters,
     &                              mgmaxlevels,mgerror,mgtol,
     &                              downpasses,uppasses,
     &                              lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                              gridmode,conductors,
     &                              my_index,nslaves,izfsslave,nzfsslave)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      use BoltzmannElectrons
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor,lbuildquads
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

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
      real(kind=8),allocatable:: phiwithbounds(:,:,:)
      integer(ISZ),allocatable:: regionid(:,:,:)
      real(kind=8),allocatable:: phisave(:,:,:)
      real(kind=8):: bendx((nx+1)*(ny+1))
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

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
        if (lbuildquads) then
          call setcndtr(xmmin,ymmin,zmmin,zmminglobal,zbeam,zgrid,
     &                  nx,ny,nz,dx,dy,dz,
     &                  bounds(4),bounds(5),bounds(1),
     &                  (bounds(0)/=1.and.bounds(2)==1),
     &                  (bounds(0)==1.and.bounds(2)==1))
        endif
        if (laddconductor) call execuser("calladdconductor")
      endif
      call checkconductors(nx,ny,nz,nzfull,dx,dy,dz,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave)

!$OMP PARALLEL
!$OMP&PRIVATE(ii,i,k,rs,x,r,ix,iy,iz)

c     --- Prepare rho by dividing it by -eps0
      rho = -rho/eps0

      allocate(phiwithbounds(-1:nx+1,-1:ny+1,-1:nz+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phiwithbounds to shape ",nx,ny,nz
        stop
      endif
      allocate(phisave(0:nx,0:ny,0:nz),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,ny,nz
        stop
      endif

c     --- Create and setup array which holds the BE region ids
      allocate(regionid(0:nx,0:ny,0:nz))
      call setupregionidsbe3d(regionid,nx,ny,nz,dx,dy,dz,
     &                        xmmin,ymmin,zmmin)

c     --- Copy phi into phiwithbounds
      phiwithbounds = 0.
      phiwithbounds(0:nx,0:ny,-1:nz+1) = phi

c     --- Fill the transverse boundaries guard cells
      if (bounds(0) == 0) phiwithbounds(-1,:,:) = 2.*phiwithbounds(0,:,:)-phiwithbounds(1,:,:)
      if (bounds(1) == 0) phiwithbounds(nx+1,:,:) = 2.*phiwithbounds(nx,:,:)-phiwithbounds(nx-1,:,:)
      if (bounds(0) == 1) phiwithbounds(-1,:,:) = phiwithbounds(1,:,:)
      if (bounds(1) == 1) phiwithbounds(nx+1,:,:) = phiwithbounds(nx-1,:,:)
      if (bounds(0) == 2) phiwithbounds(-1,:,:) = phiwithbounds(nx-1,:,:)
      if (bounds(1) == 2) phiwithbounds(nx:nx+1,:,:) = phiwithbounds(0:1,:,:)

      if (bounds(2) == 0) phiwithbounds(:,-1,:) = 2.*phiwithbounds(:,0,:)-phiwithbounds(:,1,:)
      if (bounds(3) == 0) phiwithbounds(:,ny+1,:) = 2.*phiwithbounds(:,ny,:)-phiwithbounds(:,ny-1,:)
      if (bounds(2) == 1) phiwithbounds(:,-1,:) = phiwithbounds(:,1,:)
      if (bounds(3) == 1) phiwithbounds(:,ny+1,:) = phiwithbounds(:,ny-1,:)
      if (bounds(2) == 2) phiwithbounds(:,-1,:) = phiwithbounds(:,ny-1,:)
      if (bounds(3) == 2) phiwithbounds(:,ny:ny+1,:) = phiwithbounds(:,0:1,:)

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- very few iterations are done.
      mgiters = 0
      mgerror = 2.*mgtol + 1.
      do while (mgerror > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phiwithbounds(0:nx,0:ny,0:nz)

c       --- Do one vcycle.
        call vcyclebe(0,nx,ny,nz,nzfull,dx,dy,dz,phiwithbounds,rho,regionid,
     &                xmmin,ymmin,zmmin,
     &                rstar,linbend,bendx,bounds,mgparam,mgmaxlevels,
     &                downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &                my_index,nslaves,izfsslave,nzfsslave)

c       --- Calculate the change in phi.

c       --- This line seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop below.
c       mgerror = maxval(abs(phisave - phiwithbounds))

        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nz
          do iy=0,ny
            do ix=0,nx
              mgerror = max(mgerror,abs(phisave(ix,iy,iz) - phiwithbounds(ix,iy,iz)))
            enddo
          enddo
        enddo
!$OMP END DO

#ifdef MPIPARALLEL
        if (nslaves > 1) then
c         --- calculate global sorerror
          call parallelmaxrealarray(mgerror,1)
        endif
#endif

      enddo

c     --- Copy phi back
      phi = phiwithbounds(0:nx,0:ny,:)

c     --- Make a print out.
      if (mgerror > mgtol) then
        call remark("MultigridBE: Maximum number of iterations reached")
      endif
      write(errline,20) mgerror,mgiters
  20  format("MultigridBE: Error converged to ",1pe11.3," in ",i5," v-cycles")
      call remark(errline)

      deallocate(phiwithbounds)
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
      RECURSIVE subroutine vcyclebe(mglevel,nx,ny,nz,nzfull,dx,dy,dz,
     &                              phi,rho,regionid,xmmin,ymmin,zmmin,
     &                              rstar,linbend,bendx,
     $                              globalbounds,mgparam,
     &                              mgmaxlevels,downpasses,uppasses,
     &                              lcndbndy,icndbndy,conductors,
     &                              my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
      use BoltzmannElectrons
      use Multigrid3d_diagnostic
      integer(ISZ):: mglevel
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-1:nx,-1:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      logical(ISZ):: linbend
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dysqi,dzsqi
      logical(ISZ):: partialcoarsening
      real(kind=8),allocatable:: phicoarse(:,:,:),rhocoarse(:,:,:)
      real(kind=8),allocatable:: Lphi(:,:,:),Lphicoarse(:,:,:)
      real(kind=8),allocatable:: phicoarsesave(:,:,:)
      integer(ISZ),allocatable:: regionidcoarse(:,:,:)
      integer(ISZ):: i,iszone=1
      real(kind=8):: ff
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dycoarsesqi,dzcoarsesqi
      integer(ISZ):: localbounds(0:5),localboundsc(0:5)
      integer(ISZ):: lzoffsetall(0:nslaves-1),rzoffsetall(0:nslaves-1)
      integer(ISZ):: lzoffset,rzoffset
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleftc(0:nslaves-1), izsendingleftc(0:nslaves-1)
      integer(ISZ):: whosendingrightc(0:nslaves-1),izsendingrightc(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      integer(ISZ):: allocerror

      dxsqi = 1./dx**2
      dysqi = 1./dy**2
      dzsqi = 1./dz**2

      localbounds = globalbounds

#ifdef MPIPARALLEL
      lzoffsetall = 0
      rzoffsetall = 0
      call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                      globalbounds,nzfull,
     &                      lzoffsetall,rzoffsetall,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
      if (izfsslave(my_index) > 0) localbounds(4) = -1
      if (izfsslave(my_index)+nz < nzfull) localbounds(5) = -1
#endif

c     --- Do initial relaxation iterations
      do i=1,downpasses
        call relaxbe3d(mglevel,nx,ny,nz,nzfull,phi,rho,regionid,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 localbounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. ny >= 4 .and. nzfull >= 4 .and.
     &    mglevel < mgmaxlevels) then

c       --- Calculate the size of the next coarsest grid. If the current
c       --- size is twice an odd integer, the next coarsest size will be
c       --- n/2+1, gauranteeing that it is even.
c       --- Note that nzcoarse is only used if full-coarsening is done.
        nxcoarse = nx/2
        if (mod(nxcoarse,2) == 1) nxcoarse = nxcoarse + 1
        nycoarse = ny/2
        if (mod(nycoarse,2) == 1) nycoarse = nycoarse + 1
        nzcoarse = nz/2
        if (mod(nzcoarse,2) == 1) nzcoarse = nzcoarse + 1
        nzfullcoarse = nzfull/2
        if (mod(nzfullcoarse,2) == 1) nzfullcoarse = nzfullcoarse + 1

        partialcoarsening = (dz > 4./3.*dx)
        if (partialcoarsening) then
          nzcoarse = nz
          nzfullcoarse = nzfull
        endif

        dxcoarse = dx*nx/nxcoarse
        dycoarse = dy*ny/nycoarse
        dzcoarse = dz*nzfull/nzfullcoarse
        dxcoarsesqi = 1./dxcoarse**2
        dycoarsesqi = 1./dycoarse**2
        dzcoarsesqi = 1./dzcoarse**2

        localboundsc = globalbounds

#ifdef MPIPARALLEL
c       --- Find domains in coarser grid
        call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                  nzfull,nzfullcoarse)
c       --- Reset value to corrected one
        nzcoarse = nzfsslavec(my_index)
c       --- Difference between starts and ends of coarse and fine grids.
c       --- Should only be in the range 0-2.
        lzoffsetall = (nzfullcoarse*izfsslave-nzfull*izfsslavec)
        rzoffsetall = (nzfull*(izfsslavec + nzfsslavec) -
     &                 nzfullcoarse*(izfsslave + nzfsslave))
c       --- Note that the lzoffsetall and rzoffsetall can only be used in
c       --- MPIPARALLEL sections since they will be unallocated in the
c       --- serial code. So, separate scalars are used in code which is
c       --- used in the serial version.
        lzoffset = lzoffsetall(my_index)
        rzoffset = rzoffsetall(my_index)
c       --- Get processor with which to exchange data on coarse grid
        call mggetexchangepes(nslaves,izfsslavec,nzfsslavec,my_index,
     &                        globalbounds,nzfullcoarse,
     &                        lzoffsetall,rzoffsetall,
     &                        whosendingleftc,izsendingleftc,
     &                        whosendingrightc,izsendingrightc)
        if (izfsslavec(my_index) > 0) localboundsc(4) = -1
        if (izfsslavec(my_index) + nzcoarse < nzfullcoarse) localboundsc(5) = -1
#else
        lzoffset = 0
        rzoffset = 0
#endif

c       --- Alloate new work space
        allocate(phicoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          stop
        endif
        allocate(rhocoarse(0:nxcoarse,0:nycoarse,0:nzcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          stop
        endif
        allocate(Lphi(-1:nx+1,-1:ny+1,-1:nz+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate Lphi to shape ",nx,ny,nz
          stop
        endif
        allocate(Lphicoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate Lphicoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          stop
        endif

c     --- Create and setup array which holds the BE region ids
        allocate(regionidcoarse(0:nxcoarse,0:nycoarse,0:nzcoarse),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate regionidcoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          stop
        endif
        call setupregionidsbe3d(regionidcoarse,nxcoarse,nycoarse,nzcoarse,
     &                          dxcoarse,dycoarse,dzcoarse,
     &                          xmmin,ymmin,zmmin)

c       --- Calculate the coarsened phi
c       call copyphiwithguardcells(nx,ny,nz,nzfull,phi,Lphi,localbounds,
c    &                             my_index,nslaves,izfsslave,nzfsslave,
c    &                             whosendingleft,izsendingleft,
c    &                             whosendingright,izsendingright)
        phicoarse = 0.
        call restrictbe3d(nx,ny,nz,nzfull,phi,1,1,
     &                    nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,1,1,
     &                    localbounds,localboundsc,lzoffset,
     &                    my_index,nslaves,izfsslavec,nzfsslavec,
     &                    whosendingleftc,izsendingleftc,
     &                    whosendingrightc,izsendingrightc)
        call cond_potmgbe(conductors%interior,nxcoarse,nycoarse,nzcoarse,
     &                    phicoarse,mglevel+1)
        call applytransversebcbe3d(nxcoarse,nycoarse,nzcoarse,phicoarse,globalbounds)
        call applylongitudinalbcbe3d(nxcoarse,nycoarse,nzcoarse,nzfullcoarse,
     &                               phicoarse,globalbounds)


c       --- Calculate the coarsened Lphi, putting it into rhocoarse
        call Lphibe(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,Lphi,regionid,
     &              mglevel,localbounds,
     &              lcndbndy,icndbndy,conductors,
     &              my_index,nslaves,izfsslave,nzfsslave,
     &              whosendingleft,izsendingleft,
     &              whosendingright,izsendingright)
        Lphi(0:nx,0:ny,0:nz) = rho - Lphi(0:nx,0:ny,0:nz)
        call applytransversebcbe3d(nx,ny,nz,Lphi,globalbounds)
        call applylongitudinalbcbe3d(nx,ny,nz,nzfull,Lphi,globalbounds)
#ifdef MPIPARALLEL
        call copyrhowithguardcells(nx,ny,nz,nzfull,rho,Lphi,localbounds,
     &                             my_index,nslaves,izfsslave,nzfsslave,
     &                             whosendingleft,izsendingleft,
     &                             whosendingright,izsendingright)
#endif
        call restrictbe3d(nx,ny,nz,nzfull,Lphi,1,1,
     &                    nxcoarse,nycoarse,nzcoarse,nzfullcoarse,rhocoarse,0,0,
     &                    localbounds,localbounds,lzoffset,
     &                    my_index,nslaves,izfsslavec,nzfsslavec,
     &                    whosendingleftc,izsendingleftc,
     &                    whosendingrightc,izsendingrightc)

        deallocate(Lphi)

c       --- Calculate L(R phi), adding it into rhocoarse
        call Lphibe(nxcoarse,nycoarse,nzcoarse,nzfullcoarse,
     &              dxcoarsesqi,dycoarsesqi,dzcoarsesqi,
     &              phicoarse,Lphicoarse,regionidcoarse,
     &              mglevel+1,localbounds,
     &              lcndbndy,icndbndy,conductors,
     &              my_index,nslaves,izfsslavec,nzfsslavec,
     &              whosendingleftc,izsendingleftc,
     &              whosendingrightc,izsendingrightc)
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
          stop
        endif
        phicoarsesave = phicoarse

c       --- Continue at the next coarsest level.
        call vcyclebe(mglevel+iszone,nxcoarse,nycoarse,nzcoarse,nzfullcoarse,
     &                dxcoarse,dycoarse,dzcoarse,phicoarse,rhocoarse,
     &                regionidcoarse,xmmin,ymmin,zmmin,
     &                rstar,linbend,bendx,globalbounds,mgparam,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors,
     &                my_index,nslaves,izfsslavec,nzfsslavec)

c       --- Add in the correction term.
        phicoarse = phicoarse - phicoarsesave
        if (lprintmgphimaxchange) then
          print*,"Max change in phi = ",maxval(abs(phicoarse))," at MG level ",mglevel
        endif
        deallocate(phicoarsesave)
        call expandbe3d(nx,ny,nz,nzfull,phi,
     &                  nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                  localbounds,lzoffset,
     &                  my_index,nslaves,izfsslave,nzfsslave,
     &                  whosendingleft,izsendingleft,
     &                  whosendingright,izsendingright)
        call clampphitophimaxbe3d(nx,ny,nz,phi,rho,regionid)

        deallocate(phicoarse,rhocoarse)
        deallocate(regionidcoarse)

      endif

c     --- Do final relaxation passes.
      do i=1,uppasses
        call relaxbe3d(mglevel,nx,ny,nz,nzfull,phi,rho,regionid,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 localbounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
        call clampphitophimaxbe3d(nx,ny,nz,phi,rho,regionid)
      enddo

      return
      end
c=============================================================================
      subroutine restrictbe3d(nx,ny,nz,nzfull,u,delt,delz,
     &                        nxcoarse,nycoarse,nzcoarse,nzfullcoarse,ucoarse,
     &                        delcx,delcz,
     &                        bounds,boundscoarse,lzoffset,
     &                        my_index,nslaves,izfsslavec,nzfsslavec,
     &                        whosendingleftc,izsendingleftc,
     &                        whosendingrightc,izsendingrightc)

      integer(ISZ):: nx,ny,nz,nzfull,delt,delz,delcx,delcz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: u(-delt:nx+delt,-delt:ny+delt,-delz:nz+delz)
      real(kind=8):: ucoarse(-delcx:nxcoarse+delcx,-delcx:nycoarse+delcx,-delcz:nzcoarse+delcz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: boundscoarse(0:5)
      integer(ISZ):: lzoffset
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      integer(ISZ):: whosendingleftc(0:nslaves-1), izsendingleftc(0:nslaves-1)
      integer(ISZ):: whosendingrightc(0:nslaves-1),izsendingrightc(0:nslaves-1)
      
c Restrict to a coarser grid.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixcoarse,iycoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      real(kind=8):: r,w,dx,dy,dz,dxi,dyi,dzi,wx(0:3),wy(0:3),wz(0:3)

c     --- Set the loop limits, always including edges.
      dx = 1.*nx/nxcoarse
      dy = 1.*ny/nycoarse
      dz = 1.*nzfull/nzfullcoarse
      dxi = 1.*nxcoarse/nx
      dyi = 1.*nycoarse/ny
      dzi = 1.*nzfullcoarse/nzfull

c     --- Do the loops.
!$OMP DO
      do izcoarse=0,nzcoarse
        izmin = ((izcoarse-1)*nzfull - lzoffset + 4*nzfullcoarse)/nzfullcoarse-3
        izmax = ((izcoarse+1)*nzfull - lzoffset - 1)/nzfullcoarse
        if (izmin < -delz) izmin = -delz
        if (izmax > nz+delz) izmax = nz+delz

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzfullcoarse)*dzi)
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
          if (iymin < -delt) iymin = -delt
          if (iymax > ny+delt) iymax = ny+delt

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
            if (ixmin < -delt) ixmin = -delt
            if (ixmax > nx+delt) ixmax = nx+delt

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
      subroutine expandbe3d(nx,ny,nz,nzfull,phi,
     &                      nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                      bounds,lzoffset,
     &                      my_index,nslaves,izfsslave,nzfsslave,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
      integer(ISZ):: nx,ny,nz,nzfull
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: phicoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1)
      integer(ISZ):: lzoffset
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

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
      dz = 1.*nzfullcoarse/nzfull

!$OMP DO
      do iz=izmin,izmax
        jz = int((iz*nzfullcoarse + lzoffset)/nzfull)
        wz =  1.*(iz*nzfullcoarse + lzoffset)/nzfull - jz
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
      call applylongitudinalbcbe3d(nx,ny,nz,nzfull,phi,bounds)

      return
      end
c=============================================================================
      subroutine relaxbe3d(mglevel,nx,ny,nz,nzfull,phi,rho,regionid,rstar,
     &                     dxsqi,dysqi,dzsqi,linbend,bendx,bounds,
     &                     mgparam,lcndbndy,icndbndy,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright)
      use Constant
      use ConductorTypemodule
      use BoltzmannElectrons
      integer(ISZ):: mglevel,nx,ny,nz,nzfull
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1),rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      real(kind=8):: dxsqi,dysqi,dzsqi
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

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
      real(kind=8):: rhoe,denom,Lphi,expo

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgbe(conductors%interior,nx,ny,nz,phi,mglevel)
      call condbndymgintbe3d(conductors%evensubgrid,nx,ny,nz,phi,bounds,mglevel)
      call condbndymgintbe3d(conductors%oddsubgrid ,nx,ny,nz,phi,bounds,mglevel)


c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      s_parity = mod(izfsslave(my_index),2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

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
                if (iondensity(id) .ne. 0. .and.
     &              electrontemperature(id) .ne. 0) then
                  expo = (phi(ix,iy,iz) - plasmapotential(id))/electrontemperature(id)
                  expo = min(expo,log(electrondensitymaxscale(id)))
                  rhoe = iondensity(id)*exp(expo)
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
     &                       dxsqi,dysqi,dzsqi,mgparam,bounds,
     &                       mglevel,icndbndy)
          endif
          if (parity == 1) then
           call condbndymgbe(conductors%oddsubgrid,nx,ny,nz,phi,rho,regionid,
     &                       dxsqi,dysqi,dzsqi,mgparam,bounds,
     &                       mglevel,icndbndy)
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
        call applylongitudinalbcbe3d(nx,ny,nz,nzfull,phi,bounds)

#ifdef MPIPARALLEL
        call mgexchange_phi(nx,ny,nz,nzfull,phi,bounds,0,0,1,
     &                      my_index,nslaves,izfsslave,nzfsslave,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
        call mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,bounds,1,0,1,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c       --- Exchange phi in the z guard planes
        call mgexchange_phi(nx,ny,nz,nzfull,phi,bounds,-1,0,1,
     &                      my_index,nslaves,izfsslave,nzfsslave,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
#endif

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
     &                        dxsqi,dysqi,dzsqi,
     &                        mgparam,bounds,mglevel,icndbndy)
      use Constant
      use ConductorSubGridTypemodule
      use BoltzmannElectrons
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1), rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: dxsqi,dysqi,dzsqi,mgparam
      integer(ISZ):: bounds(0:5),icndbndy

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
          if (iondensity(id) .ne. 0. .and. electrontemperature(id) .ne. 0) then
            expo = (pijk - plasmapotential(id))/electrontemperature(id)
            expo = min(expo,log(electrondensitymaxscale(id)))
            rhoe = iondensity(id)*exp(expo)
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
      subroutine Lphibe(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,Lphi,regionid,
     &                  mglevel,bounds,
     &                  lcndbndy,icndbndy,conductors,
     &                  my_index,nslaves,izfsslave,nzfsslave,
     &                  whosendingleft,izsendingleft,
     &                  whosendingright,izsendingright)
      use Constant
      use ConductorTypemodule
      use BoltzmannElectrons
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: Lphi(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      integer(ISZ):: mglevel,bounds(0:5)
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c Calculate the Lphi on the grid including the BE term.
c Note that what is calculated here is not the residual since it does not
c include the source term. This will not go to zero when convergence
c is reached.

      real(kind=8):: rhoe,expo
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
              if (iondensity(id) .ne. 0. .and.
     &            electrontemperature(id) .ne. 0) then
                expo = (phi(ix,iy,iz) - plasmapotential(id))/electrontemperature(id)
                expo = min(expo,log(electrondensitymaxscale(id)))
                rhoe = iondensity(id)*exp(expo)
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
     &                     dxsqi,dysqi,dzsqi,bounds,mglevel,icndbndy)
        call condbndyLphibe(conductors%oddsubgrid,nx,ny,nz,phi,Lphi,regionid,
     &                     dxsqi,dysqi,dzsqi,bounds,mglevel,icndbndy)
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
      if (bounds(4) == 2 .and. nz == nzfull) Lphi(:,:,-1) = Lphi(:,:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) Lphi(:,:,nz+1) = Lphi(:,:,1)
#ifdef MPIPARALLEL
      call mgexchange_phi(nx,ny,nz,nzfull,Lphi,bounds,-1,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phi(nx,ny,nz,nzfull,Lphi,bounds,-2,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phi(nx,ny,nz,nzfull,Lphi,bounds,-3,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phiperiodic(nx,ny,nz,nzfull,Lphi,bounds,0,1,3,
     &                            my_index,nslaves,izfsslave,
     &                            whosendingleft,whosendingright)
#endif

      return
      end
c=============================================================================
      subroutine condbndyLphibe(subgrid,nx,ny,nz,phi,Lphi,regionid,
     &                          dxsqi,dysqi,dzsqi,bounds,mglevel,icndbndy)
      use Constant
      use ConductorSubGridTypemodule
      use BoltzmannElectrons
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1), Lphi(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)
      real(kind=8):: dxsqi,dysqi,dzsqi
      integer(ISZ):: bounds(0:5),icndbndy

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
          if (iondensity(id) .ne. 0. .and. electrontemperature(id) .ne. 0) then
            expo = (pijk - plasmapotential(id))/electrontemperature(id)
            expo = min(expo,log(electrondensitymaxscale(id)))
            rhoe = iondensity(id)*exp(expo)
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
      subroutine copyphiwithguardcells(nx,ny,nz,nzfull,phiin,phiout,bounds,
     &                                 my_index,nslaves,izfsslave,nzfsslave,
     &                                 whosendingleft,izsendingleft,
     &                                 whosendingright,izsendingright)
      integer(ISZ):: nx,ny,nz,nzfull,bounds(0:5)
      real(kind=8):: phiin(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: phiout(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

      phiout(:,:,-1:nz+1) = phiin

#ifdef MPIPARALLEL
      call mgexchange_phi(nx,ny,nz,nzfull,phiout,bounds,-2,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phi(nx,ny,nz,nzfull,phiout,bounds,-3,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
#endif

      return
      end
c=============================================================================
      subroutine copyrhowithguardcells(nx,ny,nz,nzfull,rhoin,rhoout,bounds,
     &                                 my_index,nslaves,izfsslave,nzfsslave,
     &                                 whosendingleft,izsendingleft,
     &                                 whosendingright,izsendingright)
      integer(ISZ):: nx,ny,nz,nzfull,bounds(0:5)
      real(kind=8):: rhoin(0:nx,0:ny,0:nz)
      real(kind=8):: rhoout(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

#ifdef MPIPARALLEL
      call mgexchange_phi(nx,ny,nz,nzfull,rhoout,bounds,-1,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phi(nx,ny,nz,nzfull,rhoout,bounds,-2,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phi(nx,ny,nz,nzfull,rhoout,bounds,-3,1,3,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phiperiodic(nx,ny,nz,nzfull,rhoout,bounds,1,1,3,
     &                            my_index,nslaves,izfsslave,
     &                            whosendingleft,whosendingright)
#endif

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
      subroutine applylongitudinalbcbe3d(nx,ny,nz,nzfull,u,bounds)
      integer(ISZ):: nx,ny,nz,nzfull
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:ny+1,-1:nz+1)
      if (bounds(4) == 0) u(:,:,-1)   = 2.*u(:,:,0) - u(:,:,1)
      if (bounds(5) == 0) u(:,:,nz+1) = 2.*u(:,:,nz) - u(:,:,nz-1)
      if (bounds(4) == 1) u(:,:,-1)   = u(:,:,1)
      if (bounds(5) == 1) u(:,:,nz+1) = u(:,:,nz-1)
      if (nz == nzfull) then
        if (bounds(4) == 2) u(:,:,-1)   = u(:,:,nz-1)
        if (bounds(5) == 2) u(:,:,nz+1) = u(:,:,1)
      endif
#ifdef MPIPARALLEL
c     call mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,bounds,1,0,1,
c    &                            my_index,nslaves,izfsslave,
c    &                            whosendingleft,whosendingright)
#endif

      return
      end
c=============================================================================
      subroutine clampphitophimaxbe3d(nx,ny,nz,phi,rho,regionid)
      use Constant
      use BoltzmannElectrons
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1),rho(0:nx,0:ny,0:nz)
      integer(ISZ):: regionid(0:nx,0:ny,0:nz)

      integer(ISZ):: ix,iy,iz,id
      real(kind=8):: expo,phimax,rhomax

c     --- Clamp phi to phimax. The phimax is calculated assuming that
c     --- the max electron density will be the max of the iondensity and
c     --- and the ion particle density (rho).
      do iz=0,nz
        do iy=0,ny
          do ix=0,nx
            if (regionid(ix,iy,iz) > 0) then
              id = regionid(ix,iy,iz)
              if (iondensity(id) == 0.) cycle
              rhomax = max(1.,-rho(ix,iy,iz)*eps0/iondensity(id))
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
