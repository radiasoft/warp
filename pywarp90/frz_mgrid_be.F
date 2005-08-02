#include "top.h"
c=============================================================================
c@(#) File FRZ_MGRID_BE.F,
c version $Revision: 1.1 $, $Date: 2005/08/02 21:19:12 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the RZ multigrid field sovler which is part of the FRZ
c  package of WARP - it includes Boltzmann electrons.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigridberzf(grid,accuracy)
      use GRIDtypemodule
      use BoltzmannElectrons
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use Parallel
c     use ifcore
      type(GRIDtype):: grid
      real(kind=8):: accuracy

c     call for_set_fpe(FPE_M_TRAP_INV)

c     --- copy boundary positions from bound0, boundnz, and boundxy
      bounds(0) = grid%ixlbnd
      bounds(1) = grid%ixrbnd
      bounds(4) = grid%izlbnd
      bounds(5) = grid%izrbnd

c     --- NOTICE: This solver only really works well in SOR mode, without
c     --- any coarsening. The problem seems to the restriction and expansion
c     --- of phi near the conducting boundaries. Attempts so far to fix the
c     --- problem have not yet worked completely.
c     --- So, for now, the nlevels parameter, grid%nlevels, is fixed to 0
c     --- for the production version of the code.
      call multigridberzsolve(grid%nr,grid%nz,grid%nz,grid%dr,grid%dz,
     &                        grid%phi,grid%rho,
     &                        bounds,
     &                        grid%rmin,grid%zmin,
     &                        grid%mgparam,grid%ncycles,grid%ncmax,
     &                        0,mgerror,accuracy,
     &                        grid%npre,grid%npost,
     &                        lcndbndy,laddconductor,icndbndy,
     &                        gridmode,conductors,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)

      return
      end
c=============================================================================
      subroutine multigridberzsolve(nx,nz,nzfull,dx,dz,phi,rho,
     &                              bounds,
     &                              xmmin,zmmin,
     &                              mgparam,mgiters,mgmaxiters,
     &                              mgmaxlevels,mgerror,mgtol,
     &                              downpasses,uppasses,
     &                              lcndbndy,laddconductor,icndbndy,
     &                              gridmode,conductors,
     &                              my_index,nslaves,izfsslave,nzfsslave,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)
      use Subtimersfrz
      use ConductorTypemodule
      use GRIDtypemodule
      use Constant
      integer(ISZ):: nx,nz,nzfull
      real(kind=8):: phi(-1:nx+1,-1:nz+1)
      real(kind=8):: rho(0:nx,0:nz)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,zmmin
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor
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
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3.

      integer(ISZ):: i,ii,k,ix,iz
      real(kind=8),allocatable:: phisave(:,:)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lfrztimesubs) substarttime = wtime()

c     --- Make sure that the conductor data is consistent.
      if (gridmode == 0) then
        conductors%interior%n = 0
        conductors%evensubgrid%n = 0
        conductors%oddsubgrid%n = 0
        if (laddconductor) call execuser("calladdconductor")
      endif
      call checkconductors(nx,0,nz,nzfull,dx,dx,dz,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave)

!$OMP PARALLEL
!$OMP&PRIVATE(ii,i,k,ix,iz)

c     --- Prepare rho by dividing it by -eps0
      rho = -rho/eps0

      allocate(phisave(-1:nx+1,-1:nz+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigridrzsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,nz
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
        call vcycleberz(0,nx,nz,nzfull,dx,dz,phi,rho,
     &                  bounds,mgparam,mgmaxlevels,
     &                  downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &                  my_index,nslaves,izfsslave,nzfsslave,
     &                  iondensity,electrontemperature,plasmapotential,
     &                  electrondensitymaxscale)

c       --- Calculate the change in phi.

c       --- This line seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop below.
c       mgerror = maxval(abs(phisave - phi))

        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nz
          do ix=0,nx
            mgerror = max(mgerror,abs(phisave(ix,iz) - phi(ix,iz)))
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
      if (bounds(0) == 0) phi(-1,:) = 2.*phi(0,:) - phi(1,:)
      if (bounds(1) == 0) phi(nx+1,:) = 2.*phi(nx,:) - phi(nx-1,:)
      if (bounds(4) == 0) phi(:,-1) = 2.*phi(:,0) - phi(:,1)
      if (bounds(5) == 0) phi(:,nz+1) = 2.*phi(:,nz) - phi(:,nz-1)

c     --- Make a print out.
      if (mgerror > mgtol) then
        call remark("MultigridBE: Maximum number of iterations reached")
      endif
      write(errline,20) mgerror,mgiters
  20  format("MultigridBE-RZ: Error converged to ",1pe11.3," in ",i5," v-cycles")
      call remark(errline)

      deallocate(phisave)

c     --- Undo the change of rho
      rho = -rho*eps0

!$OMP END PARALLEL

      if (lfrztimesubs) timemultigridrzsolve = timemultigridrzsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycleberz(mglevel,nx,nz,nzfull,dx,dz,
     &                                phi,rho,bounds,
     &                                mgparam,
     &                                mgmaxlevels,downpasses,uppasses,
     &                                lcndbndy,icndbndy,conductors,
     &                                my_index,nslaves,izfsslave,nzfsslave,
     &                          iondensity,electrontemperature,plasmapotential,
     &                          electrondensitymaxscale)
      use ConductorTypemodule
      use GRIDtypemodule
      integer(ISZ):: mglevel
      integer(ISZ):: nx,nz,nzfull
      real(kind=8):: dx,dz
      real(kind=8):: phi(-1:nx+1,-1:nz+1),rho(0:nx,0:nz)
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

      real(kind=8):: dxsqi,dzsqi
      logical(ISZ):: partialcoarsening
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: res(:,:),rescoarse(:,:)
      real(kind=8),allocatable:: phicoarsesave(:,:)
      integer(ISZ):: i,iszone=1
      real(kind=8):: ff
      integer(ISZ):: nxcoarse,nzcoarse,nzfullcoarse
      real(kind=8):: dxcoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dzcoarsesqi
      integer(ISZ):: lzoffsetall(0:nslaves-1),rzoffsetall(0:nslaves-1)
      integer(ISZ):: lzoffset,rzoffset
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleftc(0:nslaves-1), izsendingleftc(0:nslaves-1)
      integer(ISZ):: whosendingrightc(0:nslaves-1),izsendingrightc(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      integer(ISZ):: allocerror
      integer(ISZ):: i1,i2,ic,ix,iz
      logical(ISZ):: firstiter

      dxsqi = 1./dx**2
      dzsqi = 1./dz**2

c     --- Do initial relaxation iterations
      do i=1,downpasses
        call relaxberz(mglevel,nx,nz,nzfull,phi,rho,dxsqi,dzsqi,bounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright,
     &                 iondensity,electrontemperature,plasmapotential,
     &                 electrondensitymaxscale)
      enddo

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. nzfull >= 4 .and. mglevel < mgmaxlevels) then

c       --- Calculate the size of the next coarsest grid. If the current
c       --- size is twice an odd integer, the next coarsest size will be
c       --- n/2+1, gauranteeing that it is even.
c       --- Note that nzcoarse is only used if full-coarsening is done.
        nxcoarse = nx/2
        if (mod(nxcoarse,2) == 1) nxcoarse = nxcoarse + 1
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
        dzcoarse = dz*nzfull/nzfullcoarse
        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2

        lzoffset = 0
        rzoffset = 0

c       --- Create a new GRIDType instance to store the work arrays

c       --- Alloate new work space
        allocate(phicoarse(-1:nxcoarse+1,-1:nzcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxcoarse,nzcoarse
          stop
        endif
        allocate(rhocoarse(0:nxcoarse,0:nzcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxcoarse,nzcoarse
          stop
        endif

        allocate(res(-1:nx+1,-1:nz+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nx,nz
          stop
        endif
        allocate(rescoarse(-1:nxcoarse+1,-1:nzcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate rescoarse to shape ",
     &           nxcoarse,nzcoarse
          stop
        endif

c       --- Calculate the coarsened phi
        call restrictberz(nx,nz,nzfull,phi,1,1,
     &                    nxcoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                    bounds,bounds,lzoffset)
        call cond_potmgberz(conductors%interior,nxcoarse,nzcoarse,
     &                      phicoarse,mglevel+1)
        call applytransversebcberz(nxcoarse,nzcoarse,phicoarse,bounds)
        call applylongitudinalbcberz(nxcoarse,nzcoarse,nzfullcoarse,
     &                               phicoarse,bounds)

c       --- Calculate the coarsened Lphi, putting it into rhocoarse
        call residualberz(nx,nz,nzfull,dxsqi,dzsqi,phi,res,
     &                    mglevel,bounds,mgparam,
     &                    lcndbndy,icndbndy,conductors,
     &                    iondensity,electrontemperature,plasmapotential,
     &                    electrondensitymaxscale)
        res(0:nx,0:nz) = rho - res(0:nx,0:nz)
        call applytransversebcberz(nx,nz,res,bounds)
        call applylongitudinalbcberz(nx,nz,nzfull,res,bounds)
        call restrictberz(nx,nz,nzfull,res,0,0,
     &                    nxcoarse,nzcoarse,nzfullcoarse,rhocoarse,
     &                    bounds,bounds,lzoffset)
c       print*,"max r1 ",mglevel,maxval(abs(rhocoarse))

        deallocate(res)

c       --- Calculate L(R phi), adding it into rhocoarse
        call residualberz(nxcoarse,nzcoarse,nzfullcoarse,
     &                    dxcoarsesqi,dzcoarsesqi,
     &                    phicoarse,rescoarse,
     &                    mglevel+1,bounds,mgparam,
     &                    lcndbndy,icndbndy,conductors,
     &                    iondensity,electrontemperature,plasmapotential,
     &                    electrondensitymaxscale)
        rhocoarse = rhocoarse + rescoarse(0:nx,0:nz)
c       print*,"max r2 ",mglevel,maxval(abs(rhocoarse))
c       call cond_potmgbezerorz(conductors%interior,
c    &                          nxcoarse,nzcoarse,
c    &                          rhocoarse,mglevel+1,0,0)

c       --- Save the current coarsened phi since it is needed to apply
c       --- the corrections after relaxations.
        allocate(phicoarsesave(-1:nxcoarse+1,-1:nzcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate phicoarsesave to shape ",
     &           nxcoarse,nzcoarse
          stop
        endif
        phicoarsesave = phicoarse

c       --- Continue at the next coarsest level.
        call vcycleberz(mglevel+iszone,nxcoarse,nzcoarse,nzfullcoarse,
     &                  dxcoarse,dzcoarse,phicoarse,rhocoarse,
     &                  bounds,mgparam,
     &                  mgmaxlevels,downpasses,uppasses,
     &                  lcndbndy,icndbndy,conductors,
     &                  my_index,nslaves,izfsslavec,nzfsslavec,
     &                  iondensity,electrontemperature,plasmapotential,
     &                  electrondensitymaxscale)

c       --- Add in the correction term.
        rescoarse = phicoarse - phicoarsesave
        call expandberz(nx,nz,nzfull,phi,
     &                  nxcoarse,nzcoarse,nzfullcoarse,rescoarse,
     &                  bounds,lzoffset)


c       deallocate(phicoarse,rhocoarse,phicoarsesave)
        deallocate(rescoarse)
        deallocate(phicoarsesave)
        deallocate(phicoarse)
        deallocate(rhocoarse)

      endif

c     --- Do final relaxation passes.
      do i=1,uppasses
        call relaxberz(mglevel,nx,nz,nzfull,phi,rho,dxsqi,dzsqi,bounds,mgparam,
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
      subroutine restrictberz(nx,nz,nzfull,u,delx,delz,
     &                        nxcoarse,nzcoarse,nzfullcoarse,ucoarse,
     &                        bounds,boundscoarse,lzoffset)
      use ConductorTypemodule
      use ConductorSubGridTypemodule
      integer(ISZ):: nx,nz,nzfull,delx,delz
      integer(ISZ):: nxcoarse,nzcoarse,nzfullcoarse
      real(kind=8):: u(-1:nx+1,-1:nz+1)
      real(kind=8):: ucoarse(-delx:nxcoarse+delx,-delz:nzcoarse+delz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: boundscoarse(0:5)
      integer(ISZ):: lzoffset
      
c Restrict to a coarser grid.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixcoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: r,w,dx,dz,dxi,dzi,wx(0:3),wz(0:3)

      dx = 1.*nx/nxcoarse
      dz = 1.*nzfull/nzfullcoarse
      dxi = 1.*nxcoarse/nx
      dzi = 1.*nzfullcoarse/nzfull

c     --- Do the loops.
!$OMP DO
      do izcoarse=0,nzcoarse
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

          do ixcoarse=0,nxcoarse
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
            w = 0.
            do iz=izmin,izmax
                do ix=ixmin,ixmax
                  r = r + wx(ix-ixmin)*wz(iz-izmin)*u(ix,iz)
                  w = w + wx(ix-ixmin)*wz(iz-izmin)
                enddo
            enddo
            ucoarse(ixcoarse,izcoarse) = r/w

          enddo
      enddo
!$OMP END DO

c     --- Set appropriate boundary values
      if (delx == 1) then
        call applytransversebcberz(nxcoarse,nzcoarse,ucoarse,boundscoarse)
      endif

      if (delz == 1) then
        call applylongitudinalbcberz(nxcoarse,nzcoarse,nzfullcoarse,
     &                               ucoarse,boundscoarse)
      endif

      return
      end
c=============================================================================
      subroutine expandberz(nx,nz,nzfull,phi,
     &                      nxcoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                      bounds,lzoffset)
      integer(ISZ):: nx,nz,nzfull
      integer(ISZ):: nxcoarse,nzcoarse,nzfullcoarse
      real(kind=8):: phi(-1:nx+1,-1:nz+1)
      real(kind=8):: phicoarse(-1:nxcoarse+1,-1:nzcoarse+1)
      integer(ISZ):: lzoffset
      integer(ISZ):: bounds(0:5)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ix,iz
      integer(ISZ):: jx,jz
      real(kind=8):: dx,dz
      real(kind=8):: wx,wz

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      izmin = 0
      izmax = nz
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz - 1

      dx = 1.*nxcoarse/nx
      dz = 1.*nzfullcoarse/nzfull

!$OMP DO
      do iz=izmin,izmax
        jz = int((iz*nzfullcoarse + lzoffset)/nzfull)
        wz =  1.*(iz*nzfullcoarse + lzoffset)/nzfull - jz
          do ix=ixmin,ixmax
            jx = int(ix*dx)
            wx =     ix*dx - jx

            phi(ix,iz) = phi(ix,iz) +
     &             (1.-wx)*(1.-wz)*phicoarse(jx  ,jz  ) +
     &                 wx *(1.-wz)*phicoarse(jx+1,jz  ) +
     &             (1.-wx)*    wz *phicoarse(jx  ,jz+1) +
     &                 wx *    wz *phicoarse(jx+1,jz+1)
          enddo
      enddo

      call applytransversebcberz(nx,nz,phi,bounds)
      call applylongitudinalbcberz(nx,nz,nzfull,phi,bounds)

      return
      end
c=============================================================================
      subroutine relaxberz(mglevel,nx,nz,nzfull,phi,rho,
     &                     dxsqi,dzsqi,bounds,
     &                     mgparam,lcndbndy,icndbndy,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright,
     &                     iondensity,electrontemperature,plasmapotential,
     &                     electrondensitymaxscale)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: mglevel,nx,nz,nzfull
      real(kind=8):: phi(-1:nx+1,-1:nz+1),rho(0:nx,0:nz)
      real(kind=8):: dxsqi,dzsqi
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
c Note that loops over all directions assume that nx is even.
c
c The arrangement of the loops was done to increase performance.  The entire
c grid is looped over as if it were a 1D array, ignoring boundaries.
c The boundaries are then reset, the previous value was destroyed.

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ix,iz
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ix1
      integer(ISZ):: i1,i2,ic
      real(kind=8):: rhoe,denom,Lphi

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgberz(conductors%interior,nx,nz,phi,mglevel)

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
      izmin = 0
      izmax = nz
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx-1
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
              iz = conductors%evensubgrid%indx(2,ic)
              conductors%evensubgrid%prevphi(ic) = phi(ix,iz)
            enddo
          else
            i1 = conductors%oddsubgrid%istart(mglevel)
            i2 = conductors%oddsubgrid%istart(mglevel+1)-1
            do ic = i1,i2
              ix = conductors%oddsubgrid%indx(0,ic)
              iz = conductors%oddsubgrid%indx(2,ic)
              conductors%oddsubgrid%prevphi(ic) = phi(ix,iz)
            enddo
          endif
        endif

!$OMP DO
        do iz=izmin,izmax

            ix1 = ixmin + mod(ixmin + iz + parity,2)
            do ix=ix1,ixmax,2

              if (iondensity .ne. 0.) then
                rhoe = iondensity*exp((phi(ix,iz)-plasmapotential)/
     &                                electrontemperature)
                rhoe = min(rhoe,electrondensitymaxscale*iondensity)
              else
                rhoe = 0.
              endif
              if (.true.) then
c               --- ZR
                if (ix == 0) then
                  denom = -4.*dxsqi - 2.*dzsqi - rhoe/eps0/electrontemperature
                  Lphi =  4.*(phi(1,iz) - phi(0,iz))*dxsqi
     &                  + (phi(0,iz-1) - 2.*phi(0,iz) + phi(0,iz+1))*dzsqi
     &                  - rhoe/eps0
                else
                  denom = -2.*dxsqi - 2.*dzsqi - rhoe/eps0/electrontemperature
                  Lphi =  ((ix-0.5)*phi(ix-1,iz) - 2*ix*phi(ix,iz) +
     &                     (ix+0.5)*phi(ix+1,iz))*dxsqi/ix
     &                  + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*dzsqi
     &                  - rhoe/eps0
                endif
              else
c               --- ZX
                denom = -2.*(dxsqi+dzsqi) - rhoe/eps0/electrontemperature
                Lphi =  (phi(ix-1,iz  )+phi(ix+1,iz  ))*dxsqi
     &               +  (phi(ix  ,iz-1)+phi(ix  ,iz+1))*dzsqi
     &                -  phi(ix,iz)*2.*(dxsqi+dzsqi) - rhoe/eps0
              endif
              phi(ix,iz) = phi(ix,iz) - mgparam*(Lphi - rho(ix,iz))/denom
            enddo
        enddo
!$OMP END DO

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
          if (parity == 0) then
           call condbndymgberz(conductors%evensubgrid,nx,nz,phi,rho,
     &                       dxsqi,dzsqi,mgparam,bounds,
     &                       mglevel,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
          endif
          if (parity == 1) then
           call condbndymgberz(conductors%oddsubgrid,nx,nz,phi,rho,
     &                       dxsqi,dzsqi,mgparam,bounds,
     &                       mglevel,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgberz(conductors%interior,nx,nz,phi,mglevel)

c       --- set phi in the guard planes
c       --- This must be done inside the loop over parities so that the
c       --- guard planes are updated with the most recent values.
        call applytransversebcberz(nx,nz,phi,bounds)
        call applylongitudinalbcberz(nx,nz,nzfull,phi,bounds)

c     --- end of loop over even and odd points
      enddo

      return
      end
c=============================================================================
      subroutine cond_potmgberz(interior,nx,nz,phi,mglevel)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nz+1)

c Set conductor points to the desired potential. The potential is used since
c at all levels, phi is being operated on directly.

      integer(ISZ):: ic,ix,iz

!$OMP DO
      do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
        ix = interior%indx(0,ic)
        iz = interior%indx(2,ic)
        phi(ix,iz) = interior%volt(ic)
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine cond_potmgbezerorz(interior,nx,nz,u,mglevel,delt,delz)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,nz,mglevel,delt,delz
      real(kind=8):: u(-delt:nx+delt,-delz:nz+delz)

c Set data at conductor points to zero.

      integer(ISZ):: ic,ix,iz

!$OMP DO
      do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
        ix = interior%indx(0,ic)
        iz = interior%indx(2,ic)
        u(ix,iz) = 0.
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine condbndymgberz(subgrid,nx,nz,phi,rho,dxsqi,dzsqi,
     &                          mgparam,bounds,mglevel,icndbndy,
     &                          iondensity,electrontemperature,plasmapotential,
     &                          electrondensitymaxscale)
      use Constant
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nz+1), rho(0:nx,0:nz)
      real(kind=8):: dxsqi,dzsqi,mgparam
      integer(ISZ):: bounds(0:5),icndbndy
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c
c Temporary variables pxm, pzm, pxp, and pzp hold
c phi(i-+1)-phi(i) at minus and plus one in each direction.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c The Cx and Cz hold the numerator of the coefficients of phi(i,j,k).
c The delx and delz hold the denominator of the coefficients of the
c full finite difference of phi.
c For icndbndy==1, these coefficients are just 1. For icndbndy==2, they
c include the dels.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: pik,pxm,pzm,pxp,pzp,denom,rhoe,Lphi
      real(kind=8):: delx,delz,Cx,Cz
      real(kind=8):: delxsqi,delzsqi
      integer(ISZ):: ic
      integer(ISZ):: ix,iz
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

        if (ix == 0  .and. bounds(0) < 1) cycle
        if (ix == nx .and. bounds(1) < 1) cycle

        if (iz == 0  .and. bounds(4) < 1) cycle
        if (iz == nz .and. bounds(5) < 1) cycle

c       --- Set temporaries with initial values.
        pik = subgrid%prevphi(ic)
        if (ix == 0) then
          pxm = 0.
          pxp = 4.*(phi(1,iz) - pik)
        else
          pxm = (1.-0.5/ix)*(phi(ix-1,iz) - pik)
          pxp = (1.+0.5/ix)*(phi(ix+1,iz) - pik)
        endif
        pzm = phi(ix  ,iz-1) - pik
        pzp = phi(ix  ,iz+1) - pik
        delx = 1.
        delz = 1.
        if (ix == 0) then
          Cx = 4.
        else
          Cx = 2.
        endif
        Cz = 2.
        dosubgrid = .false.

c       --- the point lower in x is inside the conductor
        if (dels(0,ic) < 1.) then
          if (ix > 0) then
            pxm = (1.-0.5/ix)*(volt(0,ic) - pik)/dels(0,ic)
            Cx = Cx - 1. + 1./dels(0,ic)
            if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(0,ic)
            dosubgrid = .true.
          endif
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          if (ix == 0) then
            pxp = 4.*(volt(1,ic) - pik)/dels(1,ic)
            Cx = 4./dels(1,ic)
          else
            pxp = (1.+0.5/ix)*(volt(1,ic) - pik)/dels(1,ic)
            Cx = Cx - 1. + 1./dels(1,ic)
          endif
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(1,ic)
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (dels(4,ic) < 1.) then
          pzm = (volt(4,ic) - pik)/dels(4,ic)
          Cz = Cz - 1. + 1./dels(4,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(4,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in z is inside the conductor
        if (dels(5,ic) < 1.) then
          pzp = (volt(5,ic) - pik)/dels(5,ic)
          Cz = Cz - 1. + 1./dels(5,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(5,ic)
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          delxsqi = dxsqi/delx
          delzsqi = dzsqi/delz
          if (iondensity .ne. 0.) then
            rhoe = iondensity*exp((pik-plasmapotential)/
     &                            electrontemperature)
            rhoe = min(rhoe,electrondensitymaxscale*iondensity)
          else
            rhoe = 0.
          endif
          denom = - Cx*delxsqi - Cz*delzsqi
     &            -rhoe/eps0/electrontemperature
          Lphi = (pxm+pxp)*delxsqi + (pzm+pzp)*delzsqi
     &           -rhoe/eps0
          phi(ix,iz) = pik - mgparam*(Lphi - rho(ix,iz))/denom
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine residualberz(nx,nz,nzfull,dxsqi,dzsqi,phi,res,
     &                        mglevel,bounds,mgparam,
     &                        lcndbndy,icndbndy,conductors,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nx,nz,nzfull
      real(kind=8):: dxsqi,dzsqi
      real(kind=8):: phi(-1:nx+1,-1:nz+1)
      real(kind=8):: res(-1:nx+1,-1:nz+1)
      integer(ISZ):: mglevel,bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Calculate the residual on the grid including the BE term.

      real(kind=8):: rhoe
      integer(ISZ):: ix,iz,ic
      integer(ISZ):: ixmin,ixmax,izmin,izmax

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      izmin = 0
      izmax = nz
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx-1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

      res = 0.
c     --- Calculate the residual.
!$OMP DO
      do iz=izmin,izmax

          do ix=ixmin,ixmax

            if (iondensity .ne. 0.) then
              rhoe = iondensity*exp((phi(ix,iz)-plasmapotential)/
     &                              electrontemperature)
              rhoe = min(rhoe,electrondensitymaxscale*iondensity)
            else
              rhoe = 0.
            endif

c           res(ix,iz) = (phi(ix-1,iz  ) + phi(ix+1,iz  ))*dxsqi
c    &                +  (phi(ix  ,iz-1) + phi(ix  ,iz+1))*dzsqi
c    &                -  phi(ix,iz)*2.*(dxsqi + dzsqi)
c    &                 - rhoe/eps0

            if (ix == 0) then
              res(ix,iz) = 4*(phi(1,iz) - phi(0,iz))*dxsqi
     &                   + (phi(0,iz-1) - 2.*phi(0,iz) + phi(0,iz+1))*dzsqi
     &                   - rhoe/eps0
            else
              res(ix,iz) = ((ix-0.5)*phi(ix-1,iz) - 2*ix*phi(ix,iz) +
     &                      (ix+0.5)*phi(ix+1,iz))*dxsqi/ix
     &                   + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*dzsqi
     &                   - rhoe/eps0
            endif

          enddo
      enddo
!$OMP END DO

c     --- Zero the residual inside conductors.
      call cond_potmgbezerorz(conductors%interior,nx,nz,res,mglevel,1,1)

      if (lcndbndy) then
c       --- Calculate the residual near the conductor.
        call condbndyresberz(conductors%evensubgrid,nx,nz,phi,res,
     &                       dxsqi,dzsqi,bounds,mglevel,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
        call condbndyresberz(conductors%oddsubgrid,nx,nz,phi,res,
     &                       dxsqi,dzsqi,bounds,mglevel,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
      endif

c     --- Transverse boundaries
      if (bounds(0) == 0) res(-1:0,:) = 0.
      if (bounds(1) == 0) res(nx:nx+1,:) = 0.
      if (bounds(0) == 1) res(-1,izmin:izmax) = res(1,izmin:izmax)
      if (bounds(1) == 1) res(nx+1,izmin:izmax) = res(nx-1,izmin:izmax)
      if (bounds(0) == 2) res(-1,izmin:izmax) = res(nx-1,izmin:izmax)
      if (bounds(1) == 2) res(nx:nx+1,izmin:izmax) = res(0:1,izmin:izmax)

c     --- Longitudinal boundaries
      if (bounds(4) == 0) res(:,-1:0) = 0.
      if (bounds(5) == 0) res(:,nz:nz+1) = 0.
      if (bounds(4) == 1) res(:,-1) = res(:,1)
      if (bounds(5) == 1) res(:,nz+1) = res(:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) res(:,-1) = res(:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) res(:,nz+1) = res(:,1)

      return
      end
c=============================================================================
      subroutine condbndyresberz(subgrid,nx,nz,phi,res,dxsqi,dzsqi,
     &                           bounds,mglevel,icndbndy,
     &                           iondensity,electrontemperature,plasmapotential,
     &                           electrondensitymaxscale)
      use Constant
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nz+1), res(-1:nx+1,-1:nz+1)
      real(kind=8):: dxsqi,dzsqi
      integer(ISZ):: bounds(0:5),icndbndy
      real(kind=8):: iondensity,electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c
c Temporary variables pxm, pzm, pxp, and pzp hold
c phi(i-+1)-phi(i) at minus and plus one in each direction.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c The Cx and Cz hold the numerator of the coefficients of phi(i,j,k).
c The delx and delz hold the denominator of the coefficients of the
c full finite difference of phi.
c For icndbndy==1, these coefficients are just 1. For icndbndy==2, they
c include the dels.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: pik,pxm,pzm,pxp,pzp,denom,rhoe,Lphi,ppp
      real(kind=8):: delx,delz
      real(kind=8):: delxsqi,delzsqi
      integer(ISZ):: ic
      integer(ISZ):: ix,iz
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

        if (ix == 0  .and. bounds(0) < 1) cycle
        if (ix == nx .and. bounds(1) < 1) cycle

        if (iz == 0  .and. bounds(4) < 1) cycle
        if (iz == nz .and. bounds(5) < 1) cycle

c       --- Set temporaries with initial values.
        pik = phi(ix,iz)
        if (ix == 0) then
          pxm = 0.
          pxp = 4.*(phi(1,iz) - pik)
        else
          pxm = (1.-0.5/ix)*(phi(ix-1,iz) - pik)
          pxp = (1.+0.5/ix)*(phi(ix+1,iz) - pik)
        endif
        pzm = phi(ix  ,iz-1) - pik
        pzp = phi(ix  ,iz+1) - pik
        delx = 1.
        delz = 1.
        ppp = 1.
        dosubgrid = .false.

c       --- the point lower in x is inside the conductor
        if (dels(0,ic) < 1.) then
          if (ix > 0) then
            pxm = (1.-0.5/ix)*(volt(0,ic) - pik)/dels(0,ic)
            if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(0,ic)
            ppp = min(ppp,dels(0,ic))
            dosubgrid = .true.
          endif
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          if (ix == 0) then
            pxp = 4.*(volt(1,ic) - pik)/dels(1,ic)
          else
            pxp = (1.+0.5/ix)*(volt(1,ic) - pik)/dels(1,ic)
          endif
          if (icndbndy == 2) delx = delx - 0.5 + 0.5*dels(1,ic)
          ppp = min(ppp,dels(1,ic))
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (dels(4,ic) < 1.) then
          pzm = (volt(4,ic) - pik)/dels(4,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(4,ic)
          ppp = min(ppp,dels(4,ic))
          dosubgrid = .true.
        endif
c       --- the point higher in z is inside the conductor
        if (dels(5,ic) < 1.) then
          pzp = (volt(5,ic) - pik)/dels(5,ic)
          if (icndbndy == 2) delz = delz - 0.5 + 0.5*dels(5,ic)
          ppp = min(ppp,dels(5,ic))
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          delxsqi = dxsqi/delx
          delzsqi = dzsqi/delz
          if (iondensity .ne. 0.) then
            rhoe = iondensity*exp((pik-plasmapotential)/
     &                            electrontemperature)
            rhoe = min(rhoe,electrondensitymaxscale*iondensity)
          else
            rhoe = 0.
          endif
          res(ix,iz) = (+ (pxm+pxp)*delxsqi
     &                  + (pzm+pzp)*delzsqi
     &                  - rhoe/eps0)*ppp
        endif
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine copyrhowithguardcellsrz(nx,nz,nzfull,rhoin,rhoout,bounds)
      integer(ISZ):: nx,nz,nzfull,bounds(0:5)
      real(kind=8):: rhoin(0:nx,0:nz)
      real(kind=8):: rhoout(-1:nx+1,-1:nz+1)

      rhoout(0:nx,0:nz) = rhoin - rhoout(0:nx,0:nz)

c     --- Transverse boundaries
      call applytransversebcberz(nx,nz,rhoout,bounds)
      call applylongitudinalbcberz(nx,nz,nzfull,rhoout,bounds)

      return
      end
c=============================================================================
      subroutine applytransversebcberz(nx,nz,u,bounds)
      integer(ISZ):: nx,nz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:nz+1)
      if (bounds(0) == 0) u(-1,:)   = u(0,:)
      if (bounds(1) == 0) u(nx+1,:) = u(nx,:)
      if (bounds(0) == 1) u(-1,:)   = u(1,:)
      if (bounds(1) == 1) u(nx+1,:) = u(nx-1,:)
      if (bounds(0) == 2) u(-1,:)   = u(nx-1,:)
      if (bounds(1) == 2) u(nx+1,:) = u(1,:)
      return
      end
c=============================================================================
      subroutine applylongitudinalbcberz(nx,nz,nzfull,u,bounds)
      integer(ISZ):: nx,nz,nzfull
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:nz+1)
      if (bounds(4) == 0) u(:,-1)   = u(:,0)
      if (bounds(5) == 0) u(:,nz+1) = u(:,nz)
      if (bounds(4) == 1) u(:,-1)   = u(:,1)
      if (bounds(5) == 1) u(:,nz+1) = u(:,nz-1)
      if (nz == nzfull) then
        if (bounds(4) == 2) u(:,-1)   = u(:,nz-1)
        if (bounds(5) == 2) u(:,nz+1) = u(:,1)
      endif
      return
      end
c=============================================================================
