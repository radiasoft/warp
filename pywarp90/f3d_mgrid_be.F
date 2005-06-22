#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID_BE.M,
c version $Revision: 1.9 $, $Date: 2005/06/22 17:53:40 $
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
     &                          xmmin,ymmin,zmmin,zbeam,zgrid,
     &                          iondensity,electrontemperature,plasmapotential,
     &                          electrondensitymaxscale)
      use GlobalVars
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use Parallel
c     use ifcore
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      integer(ISZ):: bound0,boundnz,boundxy
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c     call for_set_fpe(FPE_M_TRAP_INV)

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
     &                        xmmin,ymmin,zmmin,zbeam,zgrid,
     &                        mgparam,mgiters,mgmaxiters,
     &                        mgmaxlevels,mgerror,mgtol,
     &                        downpasses,uppasses,
     &                        lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                        gridmode,conductors,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)

      return
      end
c=============================================================================
      subroutine multigridbe3dsolve(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                              rstar,linbend,bounds,
     &                              xmmin,ymmin,zmmin,zbeam,zgrid,
     &                              mgparam,mgiters,mgmaxiters,
     &                              mgmaxlevels,mgerror,mgtol,
     &                              downpasses,uppasses,
     &                              lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                              gridmode,conductors,
     &                              my_index,nslaves,izfsslave,nzfsslave,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor,lbuildquads
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

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
          call setcndtr(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
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

      allocate(phisave(0:nx,0:ny,-1:nz+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,ny,nz
        stop
      endif

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- very few iterations are done.
      mgiters = 0
      mgerror = 2.*mgtol + 1.
      do while (mgerror > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phi

c       --- Do one vcycle.
        call vcyclebe(0,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                rstar,linbend,bendx,bounds,mgparam,mgmaxlevels,
     &                downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &                my_index,nslaves,izfsslave,nzfsslave,
     &                iondensity,electrontemperature,plasmapotential,
     &                electrondensitymaxscale)

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

#ifdef MPIPARALLEL
        if (nslaves > 1) then
c         --- calculate global sorerror
          call parallelmaxrealarray(mgerror,1)
        endif
#endif

      enddo

c     --- For Dirichlet boundary conditions, copy data into guard planes
c     --- For other boundary conditions, the guard planes are used during
c     --- the solve are so are already set.
      if (bounds(4) == 0) phi(:,:,-1) = phi(:,:,0)
      if (bounds(5) == 0) phi(:,:,nz+1) = phi(:,:,nz)

c     --- Make a print out.
      if (mgerror > mgtol) then
        call remark("MultigridBE: Maximum number of iterations reached")
      endif
      write(errline,20) mgerror,mgiters
  20  format("MultigridBE: Error converged to ",1pe11.3," in ",i5," v-cycles")
      call remark(errline)

      deallocate(phisave)

c     --- Undo the change of rho
      rho = -rho*eps0

!$OMP END PARALLEL

      if (lf3dtimesubs) timemultigrid3dsolve = timemultigrid3dsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcyclebe(mglevel,nx,ny,nz,nzfull,dx,dy,dz,
     &                              phi,rho,rstar,linbend,bendx,bounds,
     &                              mgparam,
     &                              mgmaxlevels,downpasses,uppasses,
     &                              lcndbndy,icndbndy,conductors,
     &                              my_index,nslaves,izfsslave,nzfsslave,
     &                          iondensity,electrontemperature,plasmapotential,
     &                          electrondensitymaxscale)
      use ConductorTypemodule
      integer(ISZ):: mglevel
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dysqi,dzsqi
      logical(ISZ):: partialcoarsening
      real(kind=8),allocatable:: phicoarse(:,:,:),rhocoarse(:,:,:)
      real(kind=8),allocatable:: res(:,:,:),rescoarse(:,:,:)
      real(kind=8),allocatable:: phicoarsesave(:,:,:)
      integer(ISZ):: i,iszone=1
      real(kind=8):: ff
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dycoarsesqi,dzcoarsesqi
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

c     --- Do initial relaxation iterations
      do i=1,downpasses
        call relaxbe3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 bounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright,
     &                 iondensity,electrontemperature,plasmapotential,
     &                 electrondensitymaxscale)
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

        lzoffset = 0
        rzoffset = 0

c       --- Alloate new work space
        allocate(phicoarse(0:nxcoarse,0:nycoarse,-1:nzcoarse+1),
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
        allocate(res(-1:nx+1,-1:ny+1,-1:nz+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nx,ny,nz
          stop
        endif

c       --- Calculate the coarsened phi
        call copyphiwithguardcells(nx,ny,nz,nzfull,phi,res,bounds)
        call restrictbe3d(nx,ny,nz,nzfull,res,1,
     &                  nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                  bounds,bounds,lzoffset)
        call cond_potmgbe(conductors%interior,nxcoarse,nycoarse,nzcoarse,
     &                    phicoarse,mglevel+1)

c       --- Calculate the coarsened Lphi, putting it into rhocoarse
        call residualbe(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,res,
     &                  mglevel,bounds,mgparam,
     &                  lcndbndy,icndbndy,conductors,
     &                  iondensity,electrontemperature,plasmapotential,
     &                  electrondensitymaxscale)
        call copyrhowithguardcells(nx,ny,nz,nzfull,rho,res,bounds)
        call restrictbe3d(nx,ny,nz,nzfull,res,0,
     &                  nxcoarse,nycoarse,nzcoarse,nzfullcoarse,rhocoarse,
     &                  bounds,bounds,lzoffset)

        deallocate(res)

        allocate(rescoarse(-1:nxcoarse+1,-1:nycoarse+1,-1:nzcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclebe: allocation error ",allocerror,
     &           ": could not allocate rescoarse to shape ",
     &           nxcoarse,nycoarse,nzcoarse
          stop
        endif

c       --- Calculate L(R phi), adding it into rhocoarse
        call residualbe(nxcoarse,nycoarse,nzcoarse,nzfullcoarse,
     &                  dxcoarsesqi,dycoarsesqi,dzcoarsesqi,
     &                  phicoarse,rescoarse,
     &                  mglevel+1,bounds,mgparam,
     &                  lcndbndy,icndbndy,conductors,
     &                  iondensity,electrontemperature,plasmapotential,
     &                  electrondensitymaxscale)
        rhocoarse = rhocoarse + rescoarse(0:nx,0:ny,0:nz)
        deallocate(rescoarse)
        call cond_potmgbezero(conductors%interior,
     &                        nxcoarse,nycoarse,nzcoarse,
     &                        rhocoarse,mglevel+1,0,0)

c       --- Save the current coarsened phi since it is needed to apply
c       --- the corrections after relaxations.
        allocate(phicoarsesave(0:nxcoarse,0:nycoarse,-1:nzcoarse+1),
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
     &                rstar,linbend,bendx,bounds,mgparam,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors,
     &                my_index,nslaves,izfsslavec,nzfsslavec,
     &                iondensity,electrontemperature,plasmapotential,
     &                electrondensitymaxscale)

c       --- Add in the correction term.
        phicoarse = phicoarse - phicoarsesave
        call expandbe3d(nx,ny,nz,nzfull,phi,
     &                  nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                  bounds,lzoffset)

        deallocate(phicoarse,rhocoarse,phicoarsesave)

      endif

c     --- Do final relaxation passes.
      do i=1,uppasses
        call relaxbe3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 bounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright,
     &                 iondensity,electrontemperature,plasmapotential,
     &                 electrondensitymaxscale)
      enddo

      return
      end
c=============================================================================
      subroutine restrictbe3d(nx,ny,nz,nzfull,u,delz,
     &                        nxcoarse,nycoarse,nzcoarse,nzfullcoarse,ucoarse,
     &                        bounds,boundscoarse,lzoffset)
      integer(ISZ):: nx,ny,nz,nzfull,delz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: u(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: ucoarse(0:nxcoarse,0:nycoarse,-delz:nzcoarse+delz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: boundscoarse(0:5)
      integer(ISZ):: lzoffset
      
c Restrict to a coarser grid.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixcoarse,iycoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixcoarsemin,ixcoarsemax,iycoarsemin,iycoarsemax
      integer(ISZ):: izcoarsemin,izcoarsemax
      real(kind=8):: r,w,dx,dy,dz,dxi,dyi,dzi,wx(0:3),wy(0:3),wz(0:3)

c     --- Set the loop limits, always including edges.
      ixcoarsemin = 0
      ixcoarsemax = nxcoarse
      iycoarsemin = 0
      iycoarsemax = nycoarse
      izcoarsemin = 0
      izcoarsemax = nzcoarse

      dx = 1.*nx/nxcoarse
      dy = 1.*ny/nycoarse
      dz = 1.*nzfull/nzfullcoarse
      dxi = 1.*nxcoarse/nx
      dyi = 1.*nycoarse/ny
      dzi = 1.*nzfullcoarse/nzfull
      w = dxi*dyi*dzi

c     --- Do the loops.
!$OMP DO
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = ((izcoarse-1)*nzfull - lzoffset + 4*nzfullcoarse)/nzfullcoarse-3
        izmax = ((izcoarse+1)*nzfull - lzoffset - 1)/nzfullcoarse

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

        do iycoarse=iycoarsemin,iycoarsemax
          iymin = int(ceiling((iycoarse-1)*dy + 1.e-10))
          iymax = int(floor((iycoarse+1)*dy - 1.e-10))

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

          do ixcoarse=ixcoarsemin,ixcoarsemax
            ixmin = int(ceiling((ixcoarse-1)*dx + 1.e-10))
            ixmax = int(floor((ixcoarse+1)*dx - 1.e-10))

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
            do iz=izmin,izmax
              do iy=iymin,iymax
                do ix=ixmin,ixmax
                  r = r + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)*u(ix,iy,iz)
                enddo
              enddo
            enddo
            ucoarse(ixcoarse,iycoarse,izcoarse) = w*r

          enddo
        enddo
      enddo
!$OMP END DO

c     --- Set appropriate boundary values
      if (bounds(1) == 2) ucoarse(nxcoarse,:,:) = ucoarse(0,:,:)
      if (bounds(3) == 2) ucoarse(:,nycoarse,:) = ucoarse(:,0,:)

      if (delz == 1) then
        if (boundscoarse(4) == 0) ucoarse(:,:,-1) = ucoarse(:,:,0)
        if (boundscoarse(5) == 0) ucoarse(:,:,nzcoarse+1) = ucoarse(:,:,nzcoarse)
        if (boundscoarse(4) == 1) ucoarse(:,:,-1) = ucoarse(:,:,1)
        if (boundscoarse(5) == 1) ucoarse(:,:,nzcoarse+1) = ucoarse(:,:,nzcoarse-1)
        if (boundscoarse(4) == 2) ucoarse(:,:,-1) = ucoarse(:,:,nzcoarse-1)
        if (boundscoarse(5) == 2) ucoarse(:,:,nzcoarse+1) = ucoarse(:,:,1)
      endif

      return
      end
c=============================================================================
      subroutine expandbe3d(nx,ny,nz,nzfull,phi,
     &                      nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                      bounds,lzoffset)
      integer(ISZ):: nx,ny,nz,nzfull
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: phicoarse(0:nxcoarse,0:nycoarse,-1:nzcoarse+1)
      integer(ISZ):: lzoffset
      integer(ISZ):: bounds(0:5)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: jx,jy,jz
      real(kind=8):: dx,dy,dz
      real(kind=8):: wx,wy,wz

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 1
      ixmax = nx-1
      iymin = 1
      iymax = ny-1
      izmin = 0
      izmax = nz
      if (bounds(0) >  0) ixmin = 0
      if (bounds(2) >  0) iymin = 0
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

      if (bounds(1) > 0 .or. bounds(3) > 0) then
        do iz=izmin,izmax
          jz = int((iz*nzfullcoarse + lzoffset)/nzfull)
          wz =  1.*(iz*nzfullcoarse + lzoffset)/nzfull - jz
          phi(nx,ny,iz) = phi(nx,ny,iz) +
     &          (1.-wz)*phicoarse(nxcoarse,nycoarse,jz  ) +
     &              wz *phicoarse(nxcoarse,nycoarse,jz+1)
          if (bounds(1) > 0) then
            do iy=0,ny-1
              jy = int(iy*dy)
              wy =     iy*dy - jy
              phi(nx,iy,iz) = phi(nx,iy,iz) +
     &                       (1.-wy)*(1.-wz)*phicoarse(nxcoarse,jy  ,jz  ) +
     &                           wy *(1.-wz)*phicoarse(nxcoarse,jy+1,jz  ) +
     &                       (1.-wy)*    wz *phicoarse(nxcoarse,jy  ,jz+1) +
     &                           wy *    wz *phicoarse(nxcoarse,jy+1,jz+1)
            enddo
          endif
          if (bounds(1) > 0) then
            do ix=0,nx-1
              jx = int(ix*dx)
              wx =     ix*dx - jx
              phi(ix,ny,iz) = phi(ix,ny,iz) +
     &                       (1.-wx)*(1.-wz)*phicoarse(jx  ,nycoarse,jz  ) +
     &                           wx *(1.-wz)*phicoarse(jx+1,nycoarse,jz  ) +
     &                       (1.-wx)*    wz *phicoarse(jx  ,nycoarse,jz+1) +
     &                           wx *    wz *phicoarse(jx+1,nycoarse,jz+1)
            enddo
          endif
        enddo
      endif

      if (bounds(4) == 1) phi(:,:,-1) = phi(:,:,1)
      if (bounds(5) == 1) phi(:,:,nz+1) = phi(:,:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) phi(:,:,-1) = phi(:,:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) phi(:,:,nz+1) = phi(:,:,1)

      return
      end
c=============================================================================
      subroutine relaxbe3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                     dxsqi,dysqi,dzsqi,linbend,bendx,bounds,
     &                     mgparam,lcndbndy,icndbndy,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright,
     &                     iondensity,electrontemperature,plasmapotential,
     &                     electrondensitymaxscale)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: mglevel,nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz),rstar(-1:nz+1)
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
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

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
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1
      integer(ISZ):: ix1
      integer(ISZ):: i1,i2,ic
      real(kind=8):: rhoe,denom,Lphi

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgbe(conductors%interior,nx,ny,nz,phi,mglevel)

c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      s_parity = mod(izfsslave(my_index),2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 1
      ixmax = nx-1
      iymin = 1
      iymax = ny-1
      izmin = 0
      izmax = nz
      if (bounds(0) >  0) ixmin = 0
      if (bounds(1) == 1) ixmax = nx
      if (bounds(2) >  0) iymin = 0
      if (bounds(3) == 1) iymax = ny
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

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
          izm1 = iz - 1
          izp1 = iz + 1

          do iy=iymin,iymax
            iym1 = iy - 1
            if (iy==0 .and. bounds(2) == 1) iym1 = 1
            if (iy==0 .and. bounds(2) == 2) iym1 = ny-1
            iyp1 = iy + 1
            if (iy == ny .and. bounds(3) == 1)  iyp1 = ny-1
            if (iy == ny .and. bounds(3) == 2)  iyp1 = 1

            ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
            do ix=ix1,ixmax,2
              ixm1 = ix - 1
              if (ix==0 .and. bounds(0) == 1) ixm1 = 1
              if (ix==0 .and. bounds(0) == 2) ixm1 = nx-1
              ixp1 = ix + 1
              if (ix == nx .and. bounds(1) == 1)  ixp1 = nx-1
              if (ix == nx .and. bounds(1) == 2)  ixp1 = 1

              if (iondensity .ne. 0.) then
                rhoe = iondensity*exp((phi(ix,iy,iz)-plasmapotential)/
     &                                electrontemperature)
                rhoe = min(rhoe,electrondensitymaxscale*iondensity)
              else
                rhoe = 0.
              endif
              denom = -2.*(dxsqi+dysqi+dzsqi) - rhoe/eps0/electrontemperature
              Lphi =  (phi(ixm1,iy  ,iz  )+phi(ixp1,iy  ,iz  ))*dxsqi
     &             +  (phi(ix  ,iym1,iz  )+phi(ix  ,iyp1,iz  ))*dysqi
     &             +  (phi(ix  ,iy  ,izm1)+phi(ix  ,iy  ,izp1))*dzsqi
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
           call condbndymgbe(conductors%evensubgrid,nx,ny,nz,phi,rho,
     &                       dxsqi,dysqi,dzsqi,mgparam,bounds,
     &                       mglevel,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
          endif
          if (parity == 1) then
           call condbndymgbe(conductors%oddsubgrid,nx,ny,nz,phi,rho,
     &                       dxsqi,dysqi,dzsqi,mgparam,bounds,
     &                       mglevel,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgbe(conductors%interior,nx,ny,nz,phi,mglevel)

c     --- end of loop over even and odd points
      enddo

      if (bounds(1) == 2) phi(nx,:,:) = phi(0,:,:)
      if (bounds(3) == 2) phi(:,ny,:) = phi(:,0,:)

c     --- set phi in the z guard planes
      if (bounds(4) == 1) phi(:,:,-1) = phi(:,:,0)
      if (bounds(5) == 1) phi(:,:,nz+1)   = phi(:,:,nz)
      if (bounds(4) == 1) phi(:,:,-1) = phi(:,:,1)
      if (bounds(5) == 1) phi(:,:,nz+1)   = phi(:,:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) phi(:,:,-1) = phi(:,:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) phi(:,:,nz:nz+1) = phi(:,:,0:1)

      return
      end
c=============================================================================
      subroutine cond_potmgbe(interior,nx,ny,nz,phi,mglevel)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)

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
      subroutine condbndymgbe(subgrid,nx,ny,nz,phi,rho,dxsqi,dysqi,dzsqi,
     &                        mgparam,bounds,mglevel,icndbndy,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)
      use Constant
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), rho(0:nx,0:ny,0:nz)
      real(kind=8):: dxsqi,dysqi,dzsqi,mgparam
      integer(ISZ):: bounds(0:5),icndbndy
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

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
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: pijk,pxm,pym,pzm,pxp,pyp,pzp,denom,rhoe,Lphi
      real(kind=8):: delx,dely,delz,Cx,Cy,Cz
      real(kind=8):: delxsqi,delysqi,delzsqi
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1
        izp1 = iz + 1
        izm1 = iz - 1

        if (ixm1 == -1 .and. bounds(0) == 1) ixm1 = 1
        if (ixm1 == -1 .and. bounds(0) == 2) ixm1 = nx-1
        if (ixp1 == nx+1 .and. bounds(1) == 1) ixp1 = nx-1
        if (ixp1 == nx+1 .and. bounds(1) == 2) ixp1 = 1
        if (ixm1 == -1 .or. ixp1 == nx+1) cycle

        if (iym1 == -1 .and. bounds(2) == 1) iym1 = 1
        if (iym1 == -1 .and. bounds(2) == 2) iym1 = ny-1
        if (iyp1 == ny+1 .and. bounds(3) == 1) iyp1 = ny-1
        if (iyp1 == ny+1 .and. bounds(3) == 2) iyp1 = 1
        if (iym1 == -1 .or. iyp1 == ny+1) cycle

        if (iz == 0  .and. bounds(4) < 1) cycle
        if (iz == nz .and. bounds(5) < 1) cycle

c       --- Set temporaries with initial values.
        pijk = subgrid%prevphi(ic)
        pxm = phi(ixm1,iy   ,iz   ) - pijk
        pxp = phi(ixp1,iy   ,iz   ) - pijk
        pym = phi(ix   ,iym1,iz   ) - pijk
        pyp = phi(ix   ,iyp1,iz   ) - pijk
        pzm = phi(ix   ,iy   ,izm1) - pijk
        pzp = phi(ix   ,iy   ,izp1) - pijk
        delx = 1.
        dely = 1.
        delz = 1.
        Cx = 2.
        Cy = 2.
        Cz = 2.
        dosubgrid = .false.

c       --- the point lower in x is inside the conductor
        if (dels(0,ic) < 1.) then
          pxm = (volt(0,ic) - pijk)/dels(0,ic)
          Cx = Cx - 1. + 1./dels(0,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(0,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          pxp = (volt(1,ic) - pijk)/dels(1,ic)
          Cx = Cx - 1. + 1./dels(1,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(1,ic)
          dosubgrid = .true.
        endif
c       --- the point lower in y is inside the conductor
        if (dels(2,ic) < 1.) then
          pym = (volt(2,ic) - pijk)/dels(2,ic)
          Cy = Cy - 1. + 1./dels(2,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(2,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in y is inside the conductor
        if (dels(3,ic) < 1.) then
          pyp = (volt(3,ic) - pijk)/dels(3,ic)
          Cy = Cy - 1. + 1./dels(3,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(3,ic)
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (dels(4,ic) < 1.) then
          pzm = (volt(4,ic) - pijk)/dels(4,ic)
          Cz = Cz - 1. + 1./dels(4,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(4,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in z is inside the conductor
        if (dels(5,ic) < 1.) then
          pzp = (volt(5,ic) - pijk)/dels(5,ic)
          Cz = Cz - 1. + 1./dels(5,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(5,ic)
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          delxsqi = dxsqi/delx
          delysqi = dysqi/dely
          delzsqi = dzsqi/delz
          if (iondensity .ne. 0.) then
            rhoe = iondensity*exp((pijk-plasmapotential)/
     &                            electrontemperature)
            rhoe = min(rhoe,electrondensitymaxscale*iondensity)
          else
            rhoe = 0.
          endif
          denom = -(Cx*delxsqi+Cy*delysqi+Cz*delzsqi)
     &            -rhoe/eps0/electrontemperature
          Lphi = (pxm+pxp)*delxsqi + (pym+pyp)*delysqi + (pzm+pzp)*delzsqi
     &           -rhoe/eps0
          phi(ix,iy,iz) = pijk - mgparam*(Lphi - rho(ix,iy,iz))/denom
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine residualbe(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,res,
     &                      mglevel,bounds,mgparam,
     &                      lcndbndy,icndbndy,conductors,
     &                      iondensity,electrontemperature,plasmapotential,
     &                      electrondensitymaxscale)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: res(-1:nx+1,-1:ny+1,-1:nz+1)
      integer(ISZ):: mglevel,bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Calculate the residual on the grid including the BE term.

      real(kind=8):: rhoe
      integer(ISZ):: ix,iy,iz,ic
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 1
      ixmax = nx-1
      iymin = 1
      iymax = ny-1
      izmin = 0
      izmax = nz
      if (bounds(0) >  0) ixmin = 0
      if (bounds(1) == 1) ixmax = nx
      if (bounds(2) >  0) iymin = 0
      if (bounds(3) == 1) iymax = ny
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

      res = 0.
c     --- Calculate the residual.
!$OMP DO
      do iz=izmin,izmax
        izm1 = iz - 1
        izp1 = iz + 1

        do iy=iymin,iymax
          iym1 = iy - 1
          if (iy==0 .and. bounds(2) == 1) iym1 = 1
          if (iy==0 .and. bounds(2) == 2) iym1 = ny-1
          iyp1 = iy + 1
          if (iy == ny .and. bounds(3) == 1)  iyp1 = ny-1
          if (iy == ny .and. bounds(3) == 2)  iyp1 = 1

          do ix=ixmin,ixmax
            ixm1 = ix - 1
            if (ix==0 .and. bounds(0) == 1) ixm1 = 1
            if (ix==0 .and. bounds(0) == 2) ixm1 = nx-1
            ixp1 = ix + 1
            if (ix == nx .and. bounds(1) == 1)  ixp1 = nx-1
            if (ix == nx .and. bounds(1) == 2)  ixp1 = 1

            if (iondensity .ne. 0.) then
              rhoe = iondensity*exp((phi(ix,iy,iz)-plasmapotential)/
     &                              electrontemperature)
              rhoe = min(rhoe,electrondensitymaxscale*iondensity)
            else
              rhoe = 0.
            endif
            res(ix,iy,iz) = (phi(ixm1,iy  ,iz  )+phi(ixp1,iy  ,iz  ))*dxsqi
     &                   +  (phi(ix  ,iym1,iz  )+phi(ix  ,iyp1,iz  ))*dysqi
     &                   +  (phi(ix  ,iy  ,izm1)+phi(ix  ,iy  ,izp1))*dzsqi
     &                   -  phi(ix,iy,iz)*2.*(dxsqi+dysqi+dzsqi)
     &                   - rhoe/eps0

          enddo
        enddo
      enddo
!$OMP END DO

c     --- Zero the residual inside conductors.
      call cond_potmgbezero(conductors%interior,nx,ny,nz,res,mglevel,1,1)

      if (lcndbndy) then
c       --- Calculate the residual near the conductor.
        call condbndyresbe(conductors%evensubgrid,nx,ny,nz,phi,res,
     &                     dxsqi,dysqi,dzsqi,bounds,mglevel,icndbndy,
     &                     iondensity,electrontemperature,plasmapotential,
     &                     electrondensitymaxscale)
        call condbndyresbe(conductors%oddsubgrid,nx,ny,nz,phi,res,
     &                     dxsqi,dysqi,dzsqi,bounds,mglevel,icndbndy,
     &                     iondensity,electrontemperature,plasmapotential,
     &                     electrondensitymaxscale)
      endif

c     --- Transverse boundaries
      if (bounds(0) == 0) res(-1:0,:,:) = 0.
      if (bounds(1) == 0) res(nx:nx+1,:,:) = 0.
      if (bounds(0) == 1)
     &  res(-1,iymin:iymax,izmin:izmax) = res(1,iymin:iymax,izmin:izmax)
      if (bounds(1) == 1)
     &  res(nx+1,iymin:iymax,izmin:izmax) = res(nx-1,iymin:iymax,izmin:izmax)
      if (bounds(0) == 2)
     &  res(-1,iymin:iymax,izmin:izmax) = res(nx-1,iymin:iymax,izmin:izmax)
      if (bounds(1) == 2)
     &  res(nx:nx+1,iymin:iymax,izmin:izmax) = res(0:1,iymin:iymax,izmin:izmax)

      if (bounds(2) == 0) res(:,-1:0,:) = 0.
      if (bounds(3) == 0) res(:,ny:ny+1,:) = 0.
      if (bounds(2) == 1) res(:,-1,izmin:izmax) = res(:,1,izmin:izmax)
      if (bounds(3) == 1) res(:,ny+1,izmin:izmax) = res(:,ny-1,izmin:izmax)
      if (bounds(2) == 2) res(:,-1,izmin:izmax) = res(:,ny-1,izmin:izmax)
      if (bounds(3) == 2) res(:,ny:ny+1,izmin:izmax) = res(:,0:1,izmin:izmax)

c     --- Longitudinal boundaries
      if (bounds(4) == 0) res(:,:,-1:0) = 0.
      if (bounds(5) == 0) res(:,:,nz:nz+1) = 0.
      if (bounds(4) == 1) res(:,:,-1) = res(:,:,1)
      if (bounds(5) == 1) res(:,:,nz+1) = res(:,:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) res(:,:,-1) = res(:,:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) res(:,:,nz+1) = res(:,:,1)

      return
      end
c=============================================================================
      subroutine condbndyresbe(subgrid,nx,ny,nz,phi,res,dxsqi,dysqi,dzsqi,
     &                         bounds,mglevel,icndbndy,
     &                         iondensity,electrontemperature,plasmapotential,
     &                         electrondensitymaxscale)
      use Constant
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), res(-1:nx+1,-1:ny+1,-1:nz+1)
      real(kind=8):: dxsqi,dysqi,dzsqi
      integer(ISZ):: bounds(0:5),icndbndy
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

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
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: pijk,pxm,pym,pzm,pxp,pyp,pzp,denom,rhoe,Lphi
      real(kind=8):: delx,dely,delz
      real(kind=8):: delxsqi,delysqi,delzsqi
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1
        izp1 = iz + 1
        izm1 = iz - 1

        if (ixm1 == -1 .and. bounds(0) == 1) ixm1 = 1
        if (ixm1 == -1 .and. bounds(0) == 2) ixm1 = nx-1
        if (ixp1 == nx+1 .and. bounds(1) == 1) ixp1 = nx-1
        if (ixp1 == nx+1 .and. bounds(1) == 2) ixp1 = 1
        if (ixm1 == -1 .or. ixp1 == nx+1) cycle

        if (iym1 == -1 .and. bounds(2) == 1) iym1 = 1
        if (iym1 == -1 .and. bounds(2) == 2) iym1 = ny-1
        if (iyp1 == ny+1 .and. bounds(3) == 1) iyp1 = ny-1
        if (iyp1 == ny+1 .and. bounds(3) == 2) iyp1 = 1
        if (iym1 == -1 .or. iyp1 == ny+1) cycle

        if (iz == 0  .and. bounds(4) < 1) cycle
        if (iz == nz .and. bounds(5) < 1) cycle

c       --- Set temporaries with initial values.
        pijk = phi(ix,iy,iz)
        pxm = phi(ixm1,iy   ,iz   ) - pijk
        pxp = phi(ixp1,iy   ,iz   ) - pijk
        pym = phi(ix   ,iym1,iz   ) - pijk
        pyp = phi(ix   ,iyp1,iz   ) - pijk
        pzm = phi(ix   ,iy   ,izm1) - pijk
        pzp = phi(ix   ,iy   ,izp1) - pijk
        delx = 1.
        dely = 1.
        delz = 1.
        dosubgrid = .false.

c       --- the point lower in x is inside the conductor
        if (dels(0,ic) < 1.) then
          pxm = (volt(0,ic) - pijk)/dels(0,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(0,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          pxp = (volt(1,ic) - pijk)/dels(1,ic)
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(1,ic)
          dosubgrid = .true.
        endif
c       --- the point lower in y is inside the conductor
        if (dels(2,ic) < 1.) then
          pym = (volt(2,ic) - pijk)/dels(2,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(2,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in y is inside the conductor
        if (dels(3,ic) < 1.) then
          pyp = (volt(3,ic) - pijk)/dels(3,ic)
          if (icndbndy == 2) dely = dely - 0.5 + 0.5*dels(3,ic)
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (dels(4,ic) < 1.) then
          pzm = (volt(4,ic) - pijk)/dels(4,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(4,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in z is inside the conductor
        if (dels(5,ic) < 1.) then
          pzp = (volt(5,ic) - pijk)/dels(5,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(5,ic)
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          delxsqi = dxsqi/delx
          delysqi = dysqi/dely
          delzsqi = dzsqi/delz
          if (iondensity .ne. 0.) then
            rhoe = iondensity*exp((pijk-plasmapotential)/
     &                            electrontemperature)
            rhoe = min(rhoe,electrondensitymaxscale*iondensity)
          else
            rhoe = 0.
          endif
          res(ix,iy,iz) = + (pxm+pxp)*delxsqi
     &                    + (pym+pyp)*delysqi
     &                    + (pzm+pzp)*delzsqi
     &                    - rhoe/eps0
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
c=============================================================================
      subroutine copyphiwithguardcells(nx,ny,nz,nzfull,phiin,phiout,bounds)
      integer(ISZ):: nx,ny,nz,nzfull,bounds(0:5)
      real(kind=8):: phiin(0:nx,0:ny,-1:nz+1)
      real(kind=8):: phiout(-1:nx+1,-1:ny+1,-1:nz+1)

      phiout(0:nx,0:ny,-1:nz+1) = phiin

c     --- Transverse boundaries
      if (bounds(0) == 0) phiout(-1,0:ny,:) = phiout(0,0:ny,:)
      if (bounds(1) == 0) phiout(nx+1,0:ny,:) = phiout(nx,0:ny,:)
      if (bounds(0) == 1) phiout(-1,0:ny,:) = phiout(1,0:ny,:)
      if (bounds(1) == 1) phiout(nx+1,0:ny,:) = phiout(nx-1,0:ny,:)
      if (bounds(0) == 2) phiout(-1,0:ny,:) = phiout(nx-1,0:ny,:)
      if (bounds(1) == 2) phiout(nx:nx+1,0:ny,:) = phiout(0:1,0:ny,:)

      if (bounds(2) == 0) phiout(:,-1,:) = phiout(:,0,:)
      if (bounds(3) == 0) phiout(:,ny+1,:) = phiout(:,ny,:)
      if (bounds(2) == 1) phiout(:,-1,:) = phiout(:,1,:)
      if (bounds(3) == 1) phiout(:,ny+1,:) = phiout(:,ny-1,:)
      if (bounds(2) == 2) phiout(:,-1,:) = phiout(:,ny-1,:)
      if (bounds(3) == 2) phiout(:,ny:ny+1,:) = phiout(:,0:1,:)

c     --- Longitudinal boundaries
c     if (bounds(4) == 0) phiout(:,:,-1) = phiout(:,:,0)
c     if (bounds(5) == 0) phiout(:,:,nz+1) = phiout(:,:,nz)
c     if (bounds(4) == 1) phiout(:,:,-1) = phiout(:,:,1)
c     if (bounds(5) == 1) phiout(:,:,nz+1) = phiout(:,:,nz-1)
c     if (bounds(4) == 2 .and. nz == nzfull) phiout(:,:,-1) = phiout(:,:,nz-1)
c     if (bounds(5) == 2 .and. nz == nzfull) phiout(:,:,nz+1) = phiout(:,:,1)

      return
      end
c=============================================================================
      subroutine copyrhowithguardcells(nx,ny,nz,nzfull,rhoin,rhoout,bounds)
      integer(ISZ):: nx,ny,nz,nzfull,bounds(0:5)
      real(kind=8):: rhoin(0:nx,0:ny,0:nz)
      real(kind=8):: rhoout(-1:nx+1,-1:ny+1,-1:nz+1)

      rhoout(0:nx,0:ny,0:nz) = rhoin - rhoout(0:nx,0:ny,0:nz)

c     --- Transverse boundaries
      if (bounds(0) == 0) rhoout(-1,0:ny,0:nz) = rhoout(0,0:ny,0:nz)
      if (bounds(1) == 0) rhoout(nx+1,0:ny,0:nz) = rhoout(nx,0:ny,0:nz)
      if (bounds(0) == 1) rhoout(-1,0:ny,0:nz) = rhoout(1,0:ny,0:nz)
      if (bounds(1) == 1) rhoout(nx+1,0:ny,0:nz) = rhoout(nx-1,0:ny,0:nz)
      if (bounds(0) == 2) rhoout(-1,0:ny,0:nz) = rhoout(nx-1,0:ny,0:nz)
      if (bounds(1) == 2) rhoout(nx:nx+1,0:ny,0:nz) = rhoout(0:1,0:ny,0:nz)

      if (bounds(2) == 0) rhoout(:,-1,0:nz) = rhoout(:,0,0:nz)
      if (bounds(3) == 0) rhoout(:,ny+1,0:nz) = rhoout(:,ny,0:nz)
      if (bounds(2) == 1) rhoout(:,-1,0:nz) = rhoout(:,1,0:nz)
      if (bounds(3) == 1) rhoout(:,ny+1,0:nz) = rhoout(:,ny-1,0:nz)
      if (bounds(2) == 2) rhoout(:,-1,0:nz) = rhoout(:,ny-1,0:nz)
      if (bounds(3) == 2) rhoout(:,ny:ny+1,0:nz) = rhoout(:,0:1,0:nz)

c     --- Longitudinal boundaries
      if (bounds(4) == 0) rhoout(:,:,-1) = rhoout(:,:,0)
      if (bounds(5) == 0) rhoout(:,:,nz+1) = rhoout(:,:,nz)
      if (bounds(4) == 1) rhoout(:,:,-1) = rhoout(:,:,1)
      if (bounds(5) == 1) rhoout(:,:,nz+1) = rhoout(:,:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) rhoout(:,:,-1) = rhoout(:,:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) rhoout(:,:,nz+1) = rhoout(:,:,1)

      return
      end
c=============================================================================
c=============================================================================
