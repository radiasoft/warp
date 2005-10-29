#include "top.h"
c=============================================================================
c@(#) File FRZ_MGRID_BE.F,
c version $Revision: 1.6 $, $Date: 2005/10/29 00:05:02 $
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
      use InGen3d, only:solvergeom,RZgeom
      use Parallel
c     use ifcore
      type(GRIDtype):: grid
      real(kind=8):: accuracy
      real(kind=8),allocatable:: iondensityarray(:,:)

      integer(ISZ):: ixbemin,ixbemax,izbemin,izbemax
      logical(ISZ):: lrz
c     call for_set_fpe(FPE_M_TRAP_INV)

c     --- copy boundary positions from bound0, boundnz, and boundxy
      bounds(0) = grid%ixlbnd
      bounds(1) = grid%ixrbnd
      bounds(4) = grid%izlbnd
      bounds(5) = grid%izrbnd

c     --- Set flag for whether to do cylindrical or Cartesian.
      lrz = (solvergeom == RZgeom)

c     --- Create array which holds the specified ion density data
      allocate(iondensityarray(0:grid%nr,0: grid%nz))
      iondensityarray = 0.
      ixbemin = max(0 ,int((xbemin - grid%rmin)/grid%dr))
      ixbemax = min(grid%nr,int((xbemax - grid%rmin)/grid%dr))
      izbemin = max(0 ,int((zbemin - grid%zmin)/grid%dz))
      izbemax = min(grid%nz,int((zbemax - grid%zmin)/grid%dz))
      if (ixbemax >= ixbemin .and. izbemax >= izbemin) then
        if (luseparticleldensity) then
          iondensityarray(ixbemin:ixbemax,izbemin:izbemax) =
     &        grid%rho(ixbemin:ixbemax,izbemin:izbemax)
        else
          iondensityarray(ixbemin:ixbemax,izbemin:izbemax) = iondensity
        endif
      endif

      call multigridberzsolve(grid%nr,grid%nz,grid%nz,grid%dr,grid%dz,
     &                        grid%phi,grid%rho,
     &                        bounds,
     &                        grid%rmin,grid%zmin,
     &                        grid%mgparam,grid%ncycles,grid%ncmax,
     &                        grid%nlevels,mgerror,accuracy,
     &                        grid%npre,grid%npost,
     &                        lrz,lcndbndy,laddconductor,icndbndy,
     &                        gridmode,grid%bndfirst,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        iondensityarray,electrontemperature,
     &                        plasmapotential,electrondensitymaxscale)

      deallocate(iondensityarray)

      return
      end
c=============================================================================
      subroutine multigridberzsolve(nx,nz,nzfull,dx,dz,phi,rho,
     &                              bounds,
     &                              xmmin,zmmin,
     &                              mgparam,mgiters,mgmaxiters,
     &                              mgmaxlevels,mgerror,mgtol,
     &                              downpasses,uppasses,
     &                              lrz,lcndbndy,laddconductor,icndbndy,
     &                              gridmode,bnd,
     &                              my_index,nslaves,izfsslave,nzfsslave,
     &                              iondensity,electrontemperature,
     &                              plasmapotential,electrondensitymaxscale)
      use Subtimersfrz
      use GRIDtypemodule
      use BNDtypemodule
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
      logical(ISZ):: lrz,lcndbndy,laddconductor
      integer(ISZ):: icndbndy,gridmode
      type(BNDtype):: bnd
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      real(kind=8):: iondensity(0:nx,0:nz)
      real(kind=8):: electrontemperature,plasmapotential
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
      real(kind=8):: xminodx
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lfrztimesubs) substarttime = wtime()

c     --- Make sure that the conductor data is consistent.
c     --- This is presumably already taken care of elsewhere.

!$OMP PARALLEL
!$OMP&PRIVATE(ii,i,k,ix,iz)

c     --- Prepare rho by dividing it by -eps0
      rho = -rho/eps0

c     --- Calculate grid min in grid cells.
c     --- Note that it needs to be a real number in case xmmin is not an
c     --- integer multiple of dx.
      xminodx = xmmin/dx

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
        call vcycleberz(0,nx,nz,nzfull,dx,dz,xminodx,phi,rho,
     &                  bounds,mgparam,mgmaxlevels,
     &                  downpasses,uppasses,lrz,lcndbndy,icndbndy,bnd,
     &                  my_index,nslaves,izfsslave,nzfsslave,
     &                  iondensity,electrontemperature,plasmapotential,
     &                  electrondensitymaxscale)

c       --- Calculate the change in phi (the error).

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
        call remark("MultigridBE-RZ: Maximum number of iterations reached")
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
      RECURSIVE subroutine vcycleberz(mglevel,nx,nz,nzfull,dx,dz,xminodx,
     &                                phi,rho,bounds,
     &                                mgparam,
     &                                mgmaxlevels,downpasses,uppasses,
     &                                lrz,lcndbndy,icndbndy,bnd,
     &                                my_index,nslaves,izfsslave,nzfsslave,
     &                                iondensity,electrontemperature,
     &                                plasmapotential,electrondensitymaxscale)
      use GRIDtypemodule
      use BNDtypemodule
      integer(ISZ):: mglevel
      integer(ISZ):: nx,nz,nzfull
      real(kind=8):: dx,dz,xminodx
      real(kind=8):: phi(-1:nx+1,-1:nz+1),rho(0:nx,0:nz)
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(BNDtype):: bnd
      logical(ISZ):: lrz,lcndbndy
      integer(ISZ):: icndbndy
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      real(kind=8):: iondensity(0:nx,0:nz)
      real(kind=8):: electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dzsqi
      logical(ISZ):: partialcoarsening
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: Lphi(:,:),Lphicoarse(:,:)
      real(kind=8),allocatable:: phicoarsesave(:,:)
      real(kind=8),allocatable:: iondensitycoarse(:,:)
      integer(ISZ):: i,iszone=1
      integer(ISZ):: nxcoarse,nzcoarse,nzfullcoarse
      real(kind=8):: dxcoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dzcoarsesqi
      real(kind=8):: xminodxcoarse
      integer(ISZ):: lzoffsetall(0:nslaves-1),rzoffsetall(0:nslaves-1)
      integer(ISZ):: lzoffset,rzoffset
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleftc(0:nslaves-1), izsendingleftc(0:nslaves-1)
      integer(ISZ):: whosendingrightc(0:nslaves-1),izsendingrightc(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      integer(ISZ):: allocerror
      integer(ISZ):: i1,i2,ic,ix,iz

      dxsqi = 1./dx**2
      dzsqi = 1./dz**2

c     --- Do initial relaxation iterations
      do i=1,downpasses
        call relaxberz(mglevel,nx,nz,nzfull,phi,rho,dxsqi,dzsqi,xminodx,
     &                 bounds,mgparam,lrz,lcndbndy,icndbndy,bnd,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright,
     &                 iondensity,electrontemperature,plasmapotential,
     &                 electrondensitymaxscale)
      enddo

c     --- Recurse until the bottom layer is reached.
c     --- This assumes that the parameters at the coarse levels are
c     --- already calculated.
      if (associated(bnd%next) .and. mglevel < mgmaxlevels) then

        nxcoarse = bnd%next%nr
        nzcoarse = bnd%next%nz
        nzfullcoarse = nzcoarse
        dxcoarse = bnd%next%dr
        dzcoarse = bnd%next%dz

        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2
        xminodxcoarse = xminodx*nxcoarse/nx

        lzoffset = 0
        rzoffset = 0

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

        allocate(Lphi(-1:nx+1,-1:nz+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate Lphi to shape ",nx,nz
          stop
        endif
        allocate(Lphicoarse(-1:nxcoarse+1,-1:nzcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate Lphicoarse to shape ",
     &           nxcoarse,nzcoarse
          stop
        endif
        allocate(iondensitycoarse(0:nxcoarse,0:nzcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate iondensitycoarse to shape ",
     &           nxcoarse,nzcoarse
          stop
        endif
 
c       --- Calculate the coarsened phi
        call restrictberz(nx,nz,nzfull,phi,1,1,
     &                    nxcoarse,nzcoarse,nzfullcoarse,phicoarse,1,1,
     &                    bounds,bounds,lzoffset)
        call cond_potmgberz(bnd%next,nxcoarse,nzcoarse,phicoarse,mglevel+1)
        call applytransversebcberz(nxcoarse,nzcoarse,phicoarse,bounds)
        call applylongitudinalbcberz(nxcoarse,nzcoarse,nzfullcoarse,
     &                               phicoarse,bounds)

c       --- Calculate the coarsened Lphi, putting it into rhocoarse
        call Lphiberz(nx,nz,nzfull,dxsqi,dzsqi,xminodx,phi,Lphi,
     &                mglevel,bounds,mgparam,
     &                lrz,lcndbndy,icndbndy,bnd,
     &                iondensity,electrontemperature,plasmapotential,
     &                electrondensitymaxscale)
        Lphi(0:nx,0:nz) = rho - Lphi(0:nx,0:nz)
        call applytransversebcberz(nx,nz,Lphi,bounds)
        call applylongitudinalbcberz(nx,nz,nzfull,Lphi,bounds)
        call restrictberz(nx,nz,nzfull,Lphi,1,1,
     &                    nxcoarse,nzcoarse,nzfullcoarse,rhocoarse,0,0,
     &                    bounds,bounds,lzoffset)

        deallocate(Lphi)

c       --- Calculate L(R phi), adding it into rhocoarse
        call Lphiberz(nxcoarse,nzcoarse,nzfullcoarse,
     &                dxcoarsesqi,dzcoarsesqi,xminodxcoarse,
     &                phicoarse,Lphicoarse,
     &                mglevel+1,bounds,mgparam,
     &                lrz,lcndbndy,icndbndy,bnd%next,
     &                iondensity,electrontemperature,plasmapotential,
     &                electrondensitymaxscale)
        rhocoarse = rhocoarse + Lphicoarse(0:nxcoarse,0:nzcoarse)
        call cond_potmgbezerorz(bnd%next,nxcoarse,nzcoarse,rhocoarse,
     &                          mglevel+1,0,0)

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

c       --- Get the iondensity on the coarse mesh.
        call restrictberz(nx,nz,nzfull,iondensity,0,0,
     &                    nxcoarse,nzcoarse,nzfullcoarse,iondensitycoarse,0,0,
     &                    bounds,bounds,lzoffset)

c       --- Continue at the next coarsest level.
        call vcycleberz(mglevel+iszone,nxcoarse,nzcoarse,nzfullcoarse,
     &                  dxcoarse,dzcoarse,xminodxcoarse,phicoarse,rhocoarse,
     &                  bounds,mgparam,mgmaxlevels,downpasses,uppasses,
     &                  lrz,lcndbndy,icndbndy,bnd%next,
     &                  my_index,nslaves,izfsslavec,nzfsslavec,
     &                  iondensitycoarse,electrontemperature,plasmapotential,
     &                  electrondensitymaxscale)

c       --- Add in the correction term.
        phicoarse = phicoarse - phicoarsesave
        call expandberz(nx,nz,nzfull,phi,
     &                  nxcoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                  bounds,lzoffset)

        deallocate(Lphicoarse)
        deallocate(phicoarsesave)
        deallocate(phicoarse)
        deallocate(rhocoarse)
        deallocate(iondensitycoarse)

      endif

c     --- Do final relaxation passes.
      do i=1,uppasses
        call relaxberz(mglevel,nx,nz,nzfull,phi,rho,dxsqi,dzsqi,xminodx,
     &                 bounds,mgparam,lrz,lcndbndy,icndbndy,bnd,
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
     &                        delcx,delcz,
     &                        bounds,boundscoarse,lzoffset)
      integer(ISZ):: nx,nz,nzfull,delx,delz,delcx,delcz
      integer(ISZ):: nxcoarse,nzcoarse,nzfullcoarse
      real(kind=8):: u(-delx:nx+delx,-delz:nz+delz)
      real(kind=8):: ucoarse(-delcx:nxcoarse+delcx,-delcz:nzcoarse+delcz)
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

c     --- Alternate definition of w. Other is more correct when points
c     --- are left out.
c     w = dx*dz

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
      if (delcx == 1) then
        call applytransversebcberz(nxcoarse,nzcoarse,ucoarse,boundscoarse)
      endif

      if (delcz == 1) then
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
     &           (1.-wx)*(1.-wz)*phicoarse(jx  ,jz  ) +
     &               wx *(1.-wz)*phicoarse(jx+1,jz  ) +
     &           (1.-wx)*    wz *phicoarse(jx  ,jz+1) +
     &               wx *    wz *phicoarse(jx+1,jz+1)
        enddo
      enddo

      call applytransversebcberz(nx,nz,phi,bounds)
      call applylongitudinalbcberz(nx,nz,nzfull,phi,bounds)

      return
      end
c=============================================================================
      subroutine relaxberz(mglevel,nx,nz,nzfull,phi,rho,
     &                     dxsqi,dzsqi,xminodx,bounds,
     &                     mgparam,lrz,lcndbndy,icndbndy,bnd,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright,
     &                     iondensity,electrontemperature,plasmapotential,
     &                     electrondensitymaxscale)
      use Constant
      use BNDtypemodule
      use CONDtypemodule
      integer(ISZ):: mglevel,nx,nz,nzfull
      real(kind=8):: phi(-1:nx+1,-1:nz+1),rho(0:nx,0:nz)
      real(kind=8):: dxsqi,dzsqi,xminodx
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lrz,lcndbndy
      integer(ISZ):: icndbndy
      type(BNDtype):: bnd
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      real(kind=8):: iondensity(0:nx,0:nz)
      real(kind=8):: electrontemperature,plasmapotential
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

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ix,iz
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ix1
      integer(ISZ):: ic,ii
      real(kind=8):: rr
      real(kind=8):: rhoe,denom,Lphi
      type(CONDtype):: cnd

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgberz(bnd,nx,nz,phi,mglevel)

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
          do ic = 1,bnd%nb_conductors
            if (ic == 1) cnd = bnd%cndfirst
            if (ic > 1) cnd = cnd%next
            do ii = 1,cnd%nbbnd
              ix = cnd%jj(ii) - 1
              iz = cnd%kk(ii) - 1
              cnd%phi0(ii) = phi(ix,iz)
            enddo
          enddo
        endif

        if (lrz) then
!$OMP DO
          do iz=izmin,izmax

            ix1 = ixmin + mod(ixmin + iz + parity,2)
            do ix=ix1,ixmax,2

              if (iondensity(ix,iz) == 0) then
                rhoe = 0.
              else
                rhoe = iondensity(ix,iz)*exp((phi(ix,iz)-plasmapotential)/
     &                                       electrontemperature)
                rhoe = min(rhoe,electrondensitymaxscale*iondensity(ix,iz))
              endif
              rr = ix + xminodx
              if (rr == 0.) then
                denom = -4.*dxsqi - 2.*dzsqi - rhoe/eps0/electrontemperature
                Lphi =  4.*(phi(1,iz) - phi(0,iz))*dxsqi
     &                   + (phi(0,iz-1) - 2.*phi(0,iz) + phi(0,iz+1))*dzsqi
     &                   - rhoe/eps0
              else
                denom = -2.*dxsqi - 2.*dzsqi - rhoe/eps0/electrontemperature
                Lphi =  ((rr-0.5)*phi(ix-1,iz) - 2*rr*phi(ix,iz) +
     &                   (rr+0.5)*phi(ix+1,iz))*dxsqi/rr
     &                + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*dzsqi
     &                - rhoe/eps0
              endif
              phi(ix,iz) = phi(ix,iz) - mgparam*(Lphi - rho(ix,iz))/denom
            enddo
          enddo
!$OMP END DO
        else
!$OMP DO
          do iz=izmin,izmax

            ix1 = ixmin + mod(ixmin + iz + parity,2)
            do ix=ix1,ixmax,2

              if (iondensity(ix,iz) == 0) then
                rhoe = 0.
              else
                rhoe = iondensity(ix,iz)*exp((phi(ix,iz)-plasmapotential)/
     &                                       electrontemperature)
                rhoe = min(rhoe,electrondensitymaxscale*iondensity(ix,iz))
              endif
              denom = -2.*(dxsqi+dzsqi) - rhoe/eps0/electrontemperature
              Lphi =  (phi(ix-1,iz  )+phi(ix+1,iz  ))*dxsqi
     &             +  (phi(ix  ,iz-1)+phi(ix  ,iz+1))*dzsqi
     &              -  phi(ix,iz)*2.*(dxsqi+dzsqi) - rhoe/eps0
              phi(ix,iz) = phi(ix,iz) - mgparam*(Lphi - rho(ix,iz))/denom
            enddo
          enddo
!$OMP END DO
        endif
           

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
         call condbndymgberz(bnd,parity,nx,nz,phi,rho,
     &                       dxsqi,dzsqi,xminodx,mgparam,bounds,
     &                       mglevel,lrz,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgberz(bnd,nx,nz,phi,mglevel)

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
      subroutine cond_potmgberz(bnd,nx,nz,phi,mglevel)
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: nx,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nz+1)

c Set conductor points to the desired potential. The potential is used since
c at all levels, phi is being operated on directly.

      type(CONDtype):: cnd
      integer(ISZ):: ic,ii,ix,iz

!$OMP DO
      do ic = 1,bnd%nb_conductors
        if (ic == 1) cnd = bnd%cndfirst
        if (ic > 1) cnd = cnd%next
        do ii = 1,cnd%ncond
          ix = cnd%jcond(ii) - 1
          iz = cnd%kcond(ii) - 1
          phi(ix,iz) = cnd%voltage(ii)
        enddo
!$OMP END DO
      enddo

      return
      end
c=============================================================================
      subroutine cond_potmgbezerorz(bnd,nx,nz,u,mglevel,delt,delz)
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: nx,nz,mglevel,delt,delz
      real(kind=8):: u(-delt:nx+delt,-delz:nz+delz)

c Set data at conductor points to zero.

      type(CONDtype):: cnd
      integer(ISZ):: ic,ii,ix,iz

!$OMP DO
      do ic = 1,bnd%nb_conductors
        if (ic == 1) cnd = bnd%cndfirst
        if (ic > 1) cnd = cnd%next
        do ii = 1,cnd%ncond
          ix = cnd%jcond(ii) - 1
          iz = cnd%kcond(ii) - 1
          u(ix,iz) = 0.
        enddo
!$OMP END DO
      enddo

      return
      end
c=============================================================================
      subroutine condbndymgberz(bnd,parity,nx,nz,phi,rho,dxsqi,dzsqi,xminodx,
     &                          mgparam,bounds,mglevel,lrz,icndbndy,
     &                          iondensity,electrontemperature,plasmapotential,
     &                          electrondensitymaxscale)
      use Constant
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: parity
      integer(ISZ):: nx,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nz+1), rho(0:nx,0:nz)
      real(kind=8):: dxsqi,dzsqi,xminodx,mgparam
      integer(ISZ):: bounds(0:5),icndbndy
      logical(ISZ):: lrz
      real(kind=8):: iondensity(0:nx,0:nz)
      real(kind=8):: electrontemperature,plasmapotential
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

c The ppp factor (the minimum of the dels at each grid point) is included
c here to be consistent with its inclusion in the Lphi calculation.
c Because of the structure of the solver, terms included in Lphi
c must nearly cancel the same terms here, and so must be scaled
c the same. Note that is it still unknown why the ppp factor is needed
c at all for convergence (the solver diverges with ppp=1). Also note that
c here, ppp=1 for the base level since its rho is the raw rho and does
c not include Lphi.

      type(CONDtype):: cnd
      real(kind=8):: pik,pxm,pzm,pxp,pzp,denom,rhoe,Lphi,ppp
      real(kind=8):: delx,delz,Cx,Cz,rr
      real(kind=8):: delxsqi,delzsqi
      integer(ISZ):: ic,ii,i1,i2
      integer(ISZ):: ix,iz
      logical(ISZ):: dosubgrid
      real(kind=8):: dxi,dzi,dxm,dxp,dzm,dzp

c     --- These should really be passed in
      dxi = sqrt(dxsqi)
      dzi = sqrt(dzsqi)

c     --- loop over points near surface of conductors
      do ic = 1,bnd%nb_conductors
        if (ic == 1) cnd = bnd%cndfirst
        if (ic > 1) cnd = cnd%next
        if (parity == 0) then
c         --- red is even
          i1 = 1
          i2 = cnd%nbbndred
        else
c         --- black is odd
          i1 = cnd%nbbndred + 1
          i2 = cnd%nbbnd
        endif

!$OMP DO
        do ii = i1,i2

          ix = cnd%jj(ii) - 1
          iz = cnd%kk(ii) - 1
          rr = ix + xminodx

          if (ix == 0  .and. bounds(0) < 1) cycle
          if (ix == nx .and. bounds(1) < 1) cycle

          if (iz == 0  .and. bounds(4) < 1) cycle
          if (iz == nz .and. bounds(5) < 1) cycle

c         --- Set temporaries with initial values.
          pik = cnd%phi0(ii)
          if (lrz) then
            if (rr == 0.) then
              pxm = 0.
              pxp = 4.*(phi(1,iz) - pik)
              Cx = 4.
            else
              pxm = (1.-0.5/rr)*(phi(ix-1,iz) - pik)
              pxp = (1.+0.5/rr)*(phi(ix+1,iz) - pik)
              Cx = 2. !(1.-0.5/rr) + (1.+0.5/rr)
            endif
          else
            pxm = phi(ix-1,iz) - pik
            pxp = phi(ix+1,iz) - pik
            Cx = 2.
          endif
          pzm = phi(ix,iz-1) - pik
          pzp = phi(ix,iz+1) - pik
          Cz = 2.
          delx = 1.
          delz = 1.
          ppp = 1.
          dosubgrid = .false.

          dxm = cnd%dxm(ii)*dxi
          dxp = cnd%dxp(ii)*dxi
          dzm = cnd%dzm(ii)*dzi
          dzp = cnd%dzp(ii)*dzi

c         --- the point lower in x is inside the conductor
          if (dxm < 1.) then
            if (lrz) then
              if (rr > 0.) then
                pxm = (1.-0.5/rr)*(cnd%volt0xm(ii) - pik)/dxm
                Cx = Cx + (1.-0.5/rr)*(1./dxm - 1.)
                if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
                ppp = min(ppp,dxm)
                dosubgrid = .true.
              endif
            else
              pxm = (cnd%volt0xm(ii) - pik)/dxm
              Cx = Cx - 1. + 1./dxm
              if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
              ppp = min(ppp,dxm)
              dosubgrid = .true.
            endif
          endif
c         --- the point higher in x is inside the conductor
          if (dxp < 1.) then
            if (lrz) then
              if (rr == 0.) then
                pxp = 4.*(cnd%volt0xp(ii) - pik)/dxp
                Cx = 4./dxp
              else
                pxp = (1.+0.5/rr)*(cnd%volt0xp(ii) - pik)/dxp
                Cx = Cx + (1.+0.5/rr)*(1./dxp - 1.)
              endif
            else
              pxp = (cnd%volt0xp(ii) - pik)/dxp
              Cx = Cx - 1. + 1./dxp
            endif
            if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxp
            ppp = min(ppp,dxp)
            dosubgrid = .true.
          endif
c         --- the point lower in z is inside the conductor
          if (dzm < 1.) then
            pzm = (cnd%volt0zm(ii) - pik)/dzm
            Cz = Cz - 1. + 1./dzm
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzm
            ppp = min(ppp,dzm)
            dosubgrid = .true.
          endif
c         --- the point higher in z is inside the conductor
          if (dzp < 1.) then
            pzp = (cnd%volt0zp(ii) - pik)/dzp
            Cz = Cz - 1. + 1./dzp
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzp
            ppp = min(ppp,dzp)
            dosubgrid = .true.
          endif
          if (mglevel == 0) ppp = 1.
c         --- calculate the new phi based on the boundary conditions
          if (dosubgrid) then
            delxsqi = dxsqi/delx
            delzsqi = dzsqi/delz
            if (iondensity(ix,iz) == 0) then
              rhoe = 0.
            else
              rhoe = iondensity(ix,iz)*exp((pik-plasmapotential)/
     &                                     electrontemperature)
              rhoe = min(rhoe,electrondensitymaxscale*iondensity(ix,iz))
            endif
            denom = - Cx*delxsqi - Cz*delzsqi - rhoe/eps0/electrontemperature
            Lphi = (pxm+pxp)*delxsqi + (pzm+pzp)*delzsqi - rhoe/eps0
            phi(ix,iz) = pik - mgparam*(Lphi*ppp - rho(ix,iz))/denom
          endif
        enddo
!$OMP END DO
      enddo

      return
      end
c=============================================================================
      subroutine Lphiberz(nx,nz,nzfull,dxsqi,dzsqi,xminodx,phi,lphi,
     &                        mglevel,bounds,mgparam,
     &                        lrz,lcndbndy,icndbndy,bnd,
     &                        iondensity,electrontemperature,plasmapotential,
     &                        electrondensitymaxscale)
      use Constant
      use BNDtypemodule
      integer(ISZ):: nx,nz,nzfull
      real(kind=8):: dxsqi,dzsqi,xminodx
      real(kind=8):: phi(-1:nx+1,-1:nz+1)
      real(kind=8):: lphi(-1:nx+1,-1:nz+1)
      integer(ISZ):: mglevel,bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lrz,lcndbndy
      integer(ISZ):: icndbndy
      type(BNDtype):: bnd
      real(kind=8):: iondensity(0:nx,0:nz)
      real(kind=8):: electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Calculate the Lphi on the grid including the BE term.
c Note that this is note quite the residual since it does not have the
c source term included, and so will no go to zero at convergence.

      real(kind=8):: rhoe
      real(kind=8):: rr
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

      lphi = 0.
c     --- Calculate the Lphi.
      if (lrz) then
c       --- RZ
!$OMP DO
        do iz=izmin,izmax

          do ix=ixmin,ixmax

            if (iondensity(ix,iz) == 0.) then
              rhoe = 0.
            else
              rhoe = iondensity(ix,iz)*exp((phi(ix,iz)-plasmapotential)/
     &                                     electrontemperature)
              rhoe = min(rhoe,electrondensitymaxscale*iondensity(ix,iz))
            endif

            rr = ix + xminodx
            if (rr == 0.) then
              lphi(ix,iz) = 4*(phi(1,iz) - phi(0,iz))*dxsqi
     &                   + (phi(0,iz-1) - 2.*phi(0,iz) + phi(0,iz+1))*dzsqi
     &                   - rhoe/eps0
            else
              lphi(ix,iz) = ((rr-0.5)*phi(ix-1,iz) - 2*rr*phi(ix,iz) +
     &                      (rr+0.5)*phi(ix+1,iz))*dxsqi/rr
     &                   + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*dzsqi
     &                   - rhoe/eps0
            endif

          enddo
        enddo
!$OMP END DO
      else
c       --- XZ
!$OMP DO
        do iz=izmin,izmax

          do ix=ixmin,ixmax

            if (iondensity(ix,iz) == 0.) then
              rhoe = 0.
            else
              rhoe = iondensity(ix,iz)*exp((phi(ix,iz)-plasmapotential)/
     &                                     electrontemperature)
              rhoe = min(rhoe,electrondensitymaxscale*iondensity(ix,iz))
            endif

            lphi(ix,iz) = (phi(ix-1,iz) - 2.*phi(ix,iz) + phi(ix+1,iz))*dxsqi
     &                 + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*dzsqi
     &                 - rhoe/eps0

          enddo
        enddo
!$OMP END DO
      endif

c     --- Zero the Lphi inside conductors.
      call cond_potmgbezerorz(bnd,nx,nz,lphi,mglevel,1,1)

      if (lcndbndy) then
c       --- Calculate the Lphi near the conductor.
        call condbndylphiberz(bnd,nx,nz,phi,lphi,dxsqi,dzsqi,xminodx,bounds,
     &                       mglevel,lrz,icndbndy,
     &                       iondensity,electrontemperature,plasmapotential,
     &                       electrondensitymaxscale)
      endif

c     --- Transverse boundaries
      if (bounds(0) == 0) lphi(-1:0,:) = 0.
      if (bounds(1) == 0) lphi(nx:nx+1,:) = 0.
      if (bounds(0) == 1) lphi(-1,izmin:izmax) = lphi(1,izmin:izmax)
      if (bounds(1) == 1) lphi(nx+1,izmin:izmax) = lphi(nx-1,izmin:izmax)
      if (bounds(0) == 2) lphi(-1,izmin:izmax) = lphi(nx-1,izmin:izmax)
      if (bounds(1) == 2) lphi(nx:nx+1,izmin:izmax) = lphi(0:1,izmin:izmax)

c     --- Longitudinal boundaries
      if (bounds(4) == 0) lphi(:,-1:0) = 0.
      if (bounds(5) == 0) lphi(:,nz:nz+1) = 0.
      if (bounds(4) == 1) lphi(:,-1) = lphi(:,1)
      if (bounds(5) == 1) lphi(:,nz+1) = lphi(:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) lphi(:,-1) = lphi(:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) lphi(:,nz+1) = lphi(:,1)

      return
      end
c=============================================================================
      subroutine condbndylphiberz(bnd,nx,nz,phi,lphi,dxsqi,dzsqi,xminodx,
     &                           bounds,mglevel,lrz,icndbndy,
     &                           iondensity,electrontemperature,plasmapotential,
     &                           electrondensitymaxscale)
      use Constant
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: nx,nz,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nz+1), lphi(-1:nx+1,-1:nz+1)
      real(kind=8):: dxsqi,dzsqi,xminodx
      integer(ISZ):: bounds(0:5),icndbndy
      logical(ISZ):: lrz
      real(kind=8):: iondensity(0:nx,0:nz)
      real(kind=8):: electrontemperature,plasmapotential
      real(kind=8):: electrondensitymaxscale

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c
c Temporary variables pxm, pzm, pxp, and pzp hold
c phi(i-+1)-phi(i) at minus and plus one in each direction.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c The delx and delz hold the denominator of the coefficients of the
c full finite difference of phi.
c For icndbndy==1, these coefficients are just 1. For icndbndy==2, they
c include the dels.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      type(CONDtype):: cnd
      real(kind=8):: pik,pxm,pzm,pxp,pzp,denom,rhoe,ppp
      real(kind=8):: delx,delz
      real(kind=8):: delxsqi,delzsqi
      real(kind=8):: rr
      integer(ISZ):: ic,ii
      integer(ISZ):: ix,iz
      logical(ISZ):: dosubgrid
      real(kind=8):: dxi,dzi,dxm,dxp,dzm,dzp

c     --- These should really be passed in
      dxi = sqrt(dxsqi)
      dzi = sqrt(dzsqi)

c     --- loop over points near surface of conductors
      do ic = 1,bnd%nb_conductors
        if (ic == 1) cnd = bnd%cndfirst
        if (ic > 1) cnd = cnd%next

!$OMP DO
        do ii = 1,cnd%nbbnd

          ix = cnd%jj(ii) - 1
          iz = cnd%kk(ii) - 1
          rr = ix + xminodx

          if (ix == 0  .and. bounds(0) < 1) cycle
          if (ix == nx .and. bounds(1) < 1) cycle

          if (iz == 0  .and. bounds(4) < 1) cycle
          if (iz == nz .and. bounds(5) < 1) cycle

c         --- Set temporaries with initial values.
          pik = phi(ix,iz)
          if (lrz) then
            if (rr == 0.) then
              pxm = 0.
              pxp = 4.*(phi(1,iz) - pik)
            else
              pxm = (1.-0.5/rr)*(phi(ix-1,iz) - pik)
              pxp = (1.+0.5/rr)*(phi(ix+1,iz) - pik)
            endif
          else
            pxm = phi(ix-1,iz) - pik
            pxp = phi(ix+1,iz) - pik
          endif
          pzm = phi(ix  ,iz-1) - pik
          pzp = phi(ix  ,iz+1) - pik
          delx = 1.
          delz = 1.
          dosubgrid = .false.
          ppp = 1.

          dxm = cnd%dxm(ii)*dxi
          dxp = cnd%dxp(ii)*dxi
          dzm = cnd%dzm(ii)*dzi
          dzp = cnd%dzp(ii)*dzi

c         --- the point lower in x is inside the conductor
          if (dxm < 1.) then
            if (lrz) then
              if (rr > 0.) then
                pxm = (1.-0.5/rr)*(cnd%volt0xm(ii) - pik)/dxm
                if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
                ppp = min(ppp,dxm)
                dosubgrid = .true.
              endif
            else
              pxm = (cnd%volt0xm(ii) - pik)/dxm
              if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
              ppp = min(ppp,dxm)
              dosubgrid = .true.
            endif
          endif
c         --- the point higher in x is inside the conductor
          if (dxp < 1.) then
            if (lrz) then
              if (rr == 0.) then
                pxp = 4.*(cnd%volt0xp(ii) - pik)/dxp
              else
                pxp = (1.+0.5/rr)*(cnd%volt0xp(ii) - pik)/dxp
              endif
            else
              pxp = (cnd%volt0xp(ii) - pik)/dxp
            endif
            if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxp
            ppp = min(ppp,dxp)
            dosubgrid = .true.
          endif
c         --- the point lower in z is inside the conductor
          if (dzm < 1.) then
            pzm = (cnd%volt0zm(ii) - pik)/dzm
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzm
            ppp = min(ppp,dzm)
            dosubgrid = .true.
          endif
c         --- the point higher in z is inside the conductor
          if (dzp < 1.) then
            pzp = (cnd%volt0zp(ii) - pik)/dzp
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzp
            ppp = min(ppp,dzp)
            dosubgrid = .true.
          endif
c         --- calculate the new phi based on the boundary conditions
          if (dosubgrid) then
            delxsqi = dxsqi/delx
            delzsqi = dzsqi/delz
            if (iondensity(ix,iz) == 0.) then
              rhoe = 0.
            else
              rhoe = iondensity(ix,iz)*exp((pik-plasmapotential)/
     &                                     electrontemperature)
              rhoe = min(rhoe,electrondensitymaxscale*iondensity(ix,iz))
            endif
            lphi(ix,iz) = (+ (pxm+pxp)*delxsqi
     &                     + (pzm+pzp)*delzsqi
     &                     - rhoe/eps0)*ppp
          endif
        enddo
!$OMP END DO
      enddo

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
