#include "top.h"
c=============================================================================
c@(#) File FRZ_MGRID_BE.F,
c version $Revision: 1.31 $, $Date: 2008/01/09 19:32:44 $
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
      use Constant, only:eps0
      use Parallel
      use DKInterp, only:interpdk,icalceps,iusing2deps
c     use ifcore
      type(GRIDtype):: grid
      real(kind=8):: accuracy,rmax,zmax
c      real(kind=8),allocatable:: iondensityarray(:,:)
c      real(kind=8),allocatable:: epsilon(:,:),epszfac(:,:)
      real(kind=8)aviondensity
c      real(kind=8),allocatable:: electrontemparray(:,:)
c     integer(4):: fff

      integer(ISZ):: ixbemin,ixbemax,izbemin,izbemax,id
      integer(ISZ):: i,j, ip1
      logical(ISZ):: lrz,lanyuseparticledensity
      integer(ISZ):: nlevels
      real(kind=8):: wr,wz
      logical(ISZ):: check,lzavparticledensity
c     call for_set_fpe(FPE_M_TRAP_INV)
c     fff = for_set_fpe(15_4)

c     --- copy boundary positions from bound0, boundnz, and boundxy
      bounds(0) = grid%ixlbnd
      bounds(1) = grid%ixrbnd
      bounds(4) = grid%izlbnd
      bounds(5) = grid%izrbnd

c     --- Check for luseparticledensity
      lanyuseparticledensity = .false.
      lzavparticledensity = .false.
      do id=1,nberegions
        if (luseparticledensity(id) > 0) then
          lanyuseparticledensity = .true.
          if (luseparticledensity(id) ==  2) lzavparticledensity = .true.
          liondensitygrid2d(id) = .true.
c         --- Find the extrema of all of the regions
          iondensitygrid2d%xmin = min(xbemin(id),iondensitygrid2d%xmin)
          iondensitygrid2d%xmax = max(xbemax(id),iondensitygrid2d%xmax)
          iondensitygrid2d%ymin = min(zbemin(id),iondensitygrid2d%ymin)
          iondensitygrid2d%ymax = max(zbemax(id),iondensitygrid2d%ymax)
c         --- Make sure that are within the grid bounds
          iondensitygrid2d%xmin = max(grid%rmin,iondensitygrid2d%xmin)
          iondensitygrid2d%xmax = min(grid%rmax,iondensitygrid2d%xmax)
          iondensitygrid2d%ymin = max(grid%zmin,iondensitygrid2d%ymin)
          iondensitygrid2d%ymax = min(grid%zmax,iondensitygrid2d%ymax)
        endif
        logelecdenmaxscale(id) = log(electrondensitymaxscale(id))
      enddo
c
      if (lanyuseparticledensity) then
c       --- Force the min's and max's to coincide with field grid values
c       --- Note that the min's and max's are relative to the field grid,
c       --- so for example, the location if xmin in the lab frame is
c       --- iondensitygrid2d%xmin + grid%rmin. This is done so that
c       --- grid%rmin etc are not needed when iondensitygrid2d is used.
        wr = floor((iondensitygrid2d%xmin - grid%rmin)/grid%dr)
        iondensitygrid2d%xmin = grid%dr*nint(wr)
        wz = floor((iondensitygrid2d%ymin - grid%zmin)/grid%dz)
        iondensitygrid2d%ymin = grid%dz*nint(wz)
        wr = ceiling((iondensitygrid2d%xmax - grid%rmin)/grid%dr)
        iondensitygrid2d%xmax = grid%dr*nint(wr)
        wz = ceiling((iondensitygrid2d%ymax - grid%zmin)/grid%dz)
        iondensitygrid2d%ymax = grid%dz*nint(wz)
c       --- Now dx and dz can be set.
        iondensitygrid2d%dx = grid%dr
        iondensitygrid2d%dy = grid%dz
c       --- This will calculate nx and nzlocal and make sure the input is OK.
        call setupgrid2dtype(iondensitygrid2d,check)
        if (.not. check) then
           call kaboom("ERROR: luseparticledensity is being used but the input data is inconsistent")
           return
        endif
c     NOTE there is no analagous coding for electrontemperaturegrid2d as the user should
c     set this up in a python script -- need to set all attributes of that object.
        
c       --- Now, copy the rho data into the iondensitygrid2d grid
        ixbemin = nint(iondensitygrid2d%xmin/grid%dr)
c        ixbemax = nint(iondensitygrid2d%xmax/grid%dr)
        ixbemax = ixbemin + iondensitygrid2d%nx
        izbemin = nint(iondensitygrid2d%ymin/grid%dz)
c        izbemax = nint(iondensitygrid2d%ymax/grid%dz)
        izbemax = izbemin + iondensitygrid2d%ny

c       Copy rho's into iondensitygrid2d.  If luseparticledensity = 2, then
c        use the z-averaged value of rho.
        if (lzavparticledensity) then
           do i = ixbemin,ixbemax
              ip1 = i+1
c           NOTE grid%rho is 1-based, not 0-based, hence 1 offsets
              aviondensity = 0
              do j=izbemin+2,izbemax
                 aviondensity=aviondensity+grid%rho(i+1,j)
              enddo
              aviondensity = aviondensity+0.5*(grid%rho(ip1,izbemin+1)
     &           + grid%rho(ip1,izbemax+1))
              aviondensity = aviondensity/max(1,izbemax-izbemin)
              iondensitygrid2d%grid(i,:)= aviondensity
           enddo
        else
           iondensitygrid2d%grid(:,:) = grid%rho(ixbemin+1:ixbemax+1,izbemin+1:izbemax+1)
        endif
        iondensitygrid2d%grid(:,:) = iondensitygrid2d%grid(:,:)/eps0
      endif
c

c     --- Set flag for whether to do cylindrical or Cartesian.
      lrz = (solvergeom == RZgeom)
      nlevels = grid%nlevels

c     --- Create arrays for polarization dielectric function.  These
c     will be set to 1 as default; if we are interpolating any species,
c     we calculate the dielectric in set_polarization, in w3d_interp.F
c     Note we need these whether or not there are electrons at a 
c     particular location, so this is not bounded by ixbemin,ixmbemax
c     or izbemin,izbemax
 
      iusing2deps = 0
      if (maxval(interpdk) > 0 .and. icalceps > 0) then
        call set_polarization(grid%rho,
     &    grid%nr,grid%nz,grid%dr,grid%dz,grid%rmin,grid%zmin)
        iusing2deps = 1
      endif
c      print*,"epsilon,epszfac", epsilon(0,0),epszfac(0,0)
c     The epsilon defined as above is normalized to eps0, and epszfac
c     is the ratio of the epsilon for z derivs to that for x derivs.  

      call multigridberzsolve(grid%nr,grid%nz,grid%nz,grid%dr,grid%dz,
     &                        grid%phi,grid%rho,
     &                        bounds,
     &                        grid%rmin,grid%zmin,
     &                        grid%mgparam,grid%ncycles,grid%ncmax,
     &                        nlevels,mgerror,accuracy,mgverbose,
     &                        grid%npre,grid%npost,
     &                        lrz,lcndbndy,laddconductor,icndbndy,
     &                        gridmode,grid%bndfirst,
     &                        my_index,nslaves,izfsslave,nzfsslave)


      return
      end
c=============================================================================
      subroutine multigridberzsolve(nx,nzlocal,nz,dx,dz,phi,rho,
     &                              bounds,
     &                              xmmin,zmmin,
     &                              mgparam,mgiters,mgmaxiters,
     &                              mgmaxlevels,mgerror,mgtol,mgverbose,
     &                              downpasses,uppasses,
     &                              lrz,lcndbndy,laddconductor,icndbndy,
     &                              gridmode,bnd,
     &                              my_index,nslaves,izfsslave,nzfsslave)
      use Subtimersfrz
      use GRIDtypemodule
      use BNDtypemodule
      use Constant
      use BoltzmannElectrons
      use DKInterp, only: iusing2deps,epsilon2d
      integer(ISZ):: nx,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nx,0:nzlocal)
      integer(ISZ),allocatable:: regionid(:,:)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,zmmin
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels,mgverbose
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lrz,lcndbndy,laddconductor
      integer(ISZ):: icndbndy,gridmode
      type(BNDtype):: bnd
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Use the multigrid method for solving Poisson's equation on a 3-D Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution and includes the Boltzmann electron term.
c
c When the grid cells are rectangular, semi-coarsening is done until the
c grid cell dimensions are roughly equal. Roughly equal means that
c   2/3 dx < dz < 4/3 dx
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3.

      integer(ISZ):: i,ii,k,ix,iz,j
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
      if (iusing2deps > 0) then
         rho = -rho/(eps0*epsilon2d%grid)
      else
         rho = -rho/eps0
      endif
      iondensity = iondensity/eps0

c     --- Calculate grid min in grid cells.
c     --- Note that it needs to be a real number in case xmmin is not an
c     --- integer multiple of dx.
      xminodx = xmmin/dx

      allocate(phisave(0:nx,0:nzlocal),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigridrzsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,nzlocal
        call kaboom("multigridrzsolve: allocation error")
        return
      endif

c     --- Create and setup array which holds the BE region ids
      allocate(regionid(0:nx,0:nzlocal))
      call setupregionidsbe2d(regionid,nx,nzlocal,dx,dz,xmmin,zmmin)

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- very few iterations are done.
      mgiters = 0
      mgerror = 2.*mgtol + 1.
      do while (mgerror > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phi(0:nx,0:nzlocal)

c       --- Do one vcycle.
        call vcycleberz(0,nx,nzlocal,nz,dx,dz,xminodx,phi,rho,regionid,
     &                  xmmin,zmmin,bounds,mgparam,mgmaxlevels,
     &                  downpasses,uppasses,lrz,lcndbndy,icndbndy,bnd,
     &                  my_index,nslaves,izfsslave,nzfsslave)


c       --- Calculate the change in phi (the error).

c       --- This line seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop below.
c       mgerror = maxval(abs(phisave - phi))

        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nzlocal
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
      if (bounds(5) == 0) phi(:,nzlocal+1) = 2.*phi(:,nzlocal) - phi(:,nzlocal-1)

c     --- Make a print out.
      if (mgverbose >= 1 .or. mgerror > mgtol) then
        if (mgerror > mgtol) then
          call remark("MultigridBE-RZ: Maximum number of iterations reached")
        endif
        write(errline,20) mgerror,mgiters
  20    format("MultigridBE-RZ: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif

      deallocate(phisave)
      deallocate(regionid)

c     --- Undo the change of rho
      if (iusing2deps > 0) then
         rho = -rho*eps0*epsilon2d%grid
      else
         rho = -rho*eps0
      endif
      iondensity = iondensity*eps0

!$OMP END PARALLEL

      if (lfrztimesubs) timemultigridberzsolve = timemultigridberzsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycleberz(mglevel,nx,nzlocal,nz,dx,dz,xminodx,
     &                                phi,rho,regionid,xmmin,zmmin,bounds,
     &                                mgparam,
     &                                mgmaxlevels,downpasses,uppasses,
     &                                lrz,lcndbndy,icndbndy,bnd,
     &                                my_index,nslaves,izfsslave,nzfsslave)
      use GRIDtypemodule
      use BNDtypemodule
      use BoltzmannElectrons
      use Multigrid3d_diagnostic
      integer(ISZ):: mglevel
      integer(ISZ):: nx,nzlocal,nz
      real(kind=8):: dx,dz,xminodx
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1),rho(0:nx,0:nzlocal)
      integer(ISZ):: regionid(0:nx,0:nzlocal)
      real(kind=8):: xmmin,zmmin
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(BNDtype):: bnd
      logical(ISZ):: lrz,lcndbndy
      integer(ISZ):: icndbndy
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dzsqi
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: Lphi(:,:),Lphicoarse(:,:)
      real(kind=8),allocatable:: phicoarsesave(:,:)
      integer(ISZ),allocatable:: regionidcoarse(:,:)
      integer(ISZ):: i,iszone=1
      integer(ISZ):: nxcoarse,nzlocalcoarse,nzcoarse
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
        call relaxberz(mglevel,nx,nzlocal,nz,phi,rho,regionid,
     &                 dxsqi,dzsqi,dx,dz,xminodx,
     &                 bounds,mgparam,lrz,lcndbndy,icndbndy,bnd,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo
c     --- Recurse until the bottom layer is reached.
c     --- This assumes that the parameters at the coarse levels are
c     --- already calculated.
      if (associated(bnd%next) .and. mglevel < mgmaxlevels) then

        nxcoarse = bnd%next%nr
        nzlocalcoarse = bnd%next%nz
        nzcoarse = nzlocalcoarse
        dxcoarse = bnd%next%dr
        dzcoarse = bnd%next%dz

        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2
        xminodxcoarse = xminodx*nxcoarse/nx

        lzoffset = 0
        rzoffset = 0

c       --- Alloate new work space
        allocate(phicoarse(-1:nxcoarse+1,-1:nzlocalcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcycleberz: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxcoarse,0:nzlocalcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcycleberz: allocation error")
          return
        endif

        allocate(Lphi(-1:nx+1,-1:nzlocal+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate Lphi to shape ",nx,nzlocal
          call kaboom("vcycleberz: allocation error")
          return
        endif
        allocate(Lphicoarse(-1:nxcoarse+1,-1:nzlocalcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate Lphicoarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcycleberz: allocation error")
          return
        endif
 
c     --- Create and setup array which holds the BE region ids
        allocate(regionidcoarse(0:nxcoarse,0:nzlocalcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate regionidcoarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcycleberz: allocation error")
          return
        endif
        call setupregionidsbe2d(regionidcoarse,nxcoarse,nzlocalcoarse,
     &                          dxcoarse,dzcoarse,xmmin,zmmin)

c       --- Calculate the coarsened phi
        phicoarse = 0.
        call restrictberz(nx,nzlocal,nz,phi,1,1,
     &                    nxcoarse,nzlocalcoarse,nzcoarse,phicoarse,1,1,
     &                    bounds,bounds,lzoffset,.false.)
        call cond_potmgberz(bnd%next,nxcoarse,nzlocalcoarse,phicoarse,mglevel+1)
        call applytransversebcberz(nxcoarse,nzlocalcoarse,phicoarse,bounds)
        call applylongitudinalbcberz(nxcoarse,nzlocalcoarse,nzcoarse,
     &                               phicoarse,bounds)

c       --- Calculate the coarsened Lphi, putting it into rhocoarse
        call Lphiberz(nx,nzlocal,nz,dxsqi,dzsqi,dx,dz,xminodx,phi,Lphi,regionid,
     &                mglevel,bounds,lrz,lcndbndy,icndbndy,bnd)
        Lphi(0:nx,0:nzlocal) = rho - Lphi(0:nx,0:nzlocal)
        call applytransversebcberz(nx,nzlocal,Lphi,bounds)
        call applylongitudinalbcberz(nx,nzlocal,nz,Lphi,bounds)
        call restrictberz(nx,nzlocal,nz,Lphi,1,1,
     &                    nxcoarse,nzlocalcoarse,nzcoarse,rhocoarse,0,0,
     &                    bounds,bounds,lzoffset,.false.)

        deallocate(Lphi)

c       --- Calculate L(R phi), adding it into rhocoarse
        call Lphiberz(nxcoarse,nzlocalcoarse,nzcoarse,
     &                dxcoarsesqi,dzcoarsesqi,dxcoarse,dzcoarse,xminodxcoarse,
     &                phicoarse,Lphicoarse,regionidcoarse,
     &                mglevel+1,bounds,
     &                lrz,lcndbndy,icndbndy,bnd%next)
        rhocoarse = rhocoarse + Lphicoarse(0:nxcoarse,0:nzlocalcoarse)
        call cond_potmgbezerorz(bnd%next,nxcoarse,nzlocalcoarse,rhocoarse,
     &                          mglevel+1,0,0)

c       --- Save the current coarsened phi since it is needed to apply
c       --- the corrections after relaxations.
        allocate(phicoarsesave(-1:nxcoarse+1,-1:nzlocalcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleberz: allocation error ",allocerror,
     &           ": could not allocate phicoarsesave to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcycleberz: allocation error")
          return
        endif
        phicoarsesave = phicoarse

c       --- Continue at the next coarsest level.
        call vcycleberz(mglevel+iszone,nxcoarse,nzlocalcoarse,nzcoarse,
     &                  dxcoarse,dzcoarse,xminodxcoarse,phicoarse,rhocoarse,
     &                  regionidcoarse,xmmin,zmmin,
     &                  bounds,mgparam,mgmaxlevels,downpasses,uppasses,
     &                  lrz,lcndbndy,icndbndy,bnd%next,
     &                  my_index,nslaves,izfsslavec,nzfsslavec)

c       --- Add in the correction term.
        phicoarse = phicoarse - phicoarsesave
        if (lprintmgphimaxchange) then
          print*,"Max change in phi = ",maxval(abs(phicoarse))," at MG level ",mglevel
        endif
        call expandberz(nx,nzlocal,nz,phi,
     &                  nxcoarse,nzlocalcoarse,nzcoarse,phicoarse,
     &                  bounds,lzoffset)
        call clampphitophimaxberz(nx,nzlocal,dx,dz,phi,rho,regionid)

        deallocate(Lphicoarse)
        deallocate(phicoarsesave)
        deallocate(phicoarse)
        deallocate(rhocoarse)
        deallocate(regionidcoarse)

      endif

c     --- Do final relaxation passes.
      do i=1,uppasses
        call relaxberz(mglevel,nx,nzlocal,nz,phi,rho,regionid,
     &                 dxsqi,dzsqi,dx,dz,xminodx,
     &                 bounds,mgparam,lrz,lcndbndy,icndbndy,bnd,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
        call clampphitophimaxberz(nx,nzlocal,dx,dz,phi,rho,regionid)
      enddo

      return
      end
c=============================================================================
      subroutine restrictberz(nx,nzlocal,nz,u,delx,delz,
     &                        nxcoarse,nzlocalcoarse,nzcoarse,ucoarse,
     &                        delcx,delcz,
     &                        bounds,boundscoarse,lzoffset,withoutzeros)
      integer(ISZ):: nx,nzlocal,nz,delx,delz,delcx,delcz
      integer(ISZ):: nxcoarse,nzlocalcoarse,nzcoarse
      real(kind=8):: u(-delx:nx+delx,-delz:nzlocal+delz)
      real(kind=8):: ucoarse(-delcx:nxcoarse+delcx,-delcz:nzlocalcoarse+delcz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: boundscoarse(0:5)
      integer(ISZ):: lzoffset
      logical(ISZ):: withoutzeros
      
c Restrict to a coarser grid.
c The withoutzeros option is for the restriction of the iondensity, so that
c for a constant value, it stays constant within the region of the plasma
c domain.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixcoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: r,w,dx,dz,dxi,dzi,wx(0:3),wz(0:3)

      dx = 1.*nx/nxcoarse
      dz = 1.*nz/nzcoarse
      dxi = 1.*nxcoarse/nx
      dzi = 1.*nzcoarse/nz

c     --- Do the loops.
!$OMP DO
      do izcoarse=0,nzlocalcoarse
        izmin = ((izcoarse-1)*nz - lzoffset + 4*nzcoarse)/nzcoarse-3
        izmax = ((izcoarse+1)*nz - lzoffset - 1)/nzcoarse
        if (izmin < -delz) izmin = -delz
        if (izmax > nzlocal+delz) izmax = nzlocal+delz

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzcoarse)*dzi)
        enddo

        if (izcoarse == 0 .and. bounds(4) == 0) then
          izmin = 0
          izmax = 0
          wz(0) = 2.
        else if (izcoarse == nzlocalcoarse .and. bounds(5) == 0) then
          izmin = nzlocal
          izmax = nzlocal
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
                if (withoutzeros .and. u(ix,iz) == 0.) cycle
                r = r + wx(ix-ixmin)*wz(iz-izmin)*u(ix,iz)
                w = w + wx(ix-ixmin)*wz(iz-izmin)
              enddo
          enddo
          if (w > 0.) then
            ucoarse(ixcoarse,izcoarse) = r/w
          else
            ucoarse(ixcoarse,izcoarse) = 0.
          endif

        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine expandberz(nx,nzlocal,nz,phi,
     &                      nxcoarse,nzlocalcoarse,nzcoarse,phicoarse,
     &                      bounds,lzoffset)
      integer(ISZ):: nx,nzlocal,nz
      integer(ISZ):: nxcoarse,nzlocalcoarse,nzcoarse
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1)
      real(kind=8):: phicoarse(-1:nxcoarse+1,-1:nzlocalcoarse+1)
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
      izmax = nzlocal
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nzlocal - 1

      dx = 1.*nxcoarse/nx
      dz = 1.*nzcoarse/nz

!$OMP DO
      do iz=izmin,izmax
        jz = int((iz*nzcoarse + lzoffset)/nz)
        wz =  1.*(iz*nzcoarse + lzoffset)/nz - jz
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

      call applytransversebcberz(nx,nzlocal,phi,bounds)
      call applylongitudinalbcberz(nx,nzlocal,nz,phi,bounds)

      return
      end
c=============================================================================
      subroutine relaxberz(mglevel,nx,nzlocal,nz,phi,rho,regionid,
     &                     dxsqi,dzsqi,dx,dz,xminodx,bounds,
     &                     mgparam,lrz,lcndbndy,icndbndy,bnd,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright)
      use BNDtypemodule
      use CONDtypemodule
      use BoltzmannElectrons
      integer(ISZ):: mglevel,nx,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1),rho(0:nx,0:nzlocal)
      integer(ISZ):: regionid(0:nx,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,dx,dz,xminodx
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lrz,lcndbndy
      integer(ISZ):: icndbndy
      type(BNDtype):: bnd
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

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ix,iz,ixx,izz
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ix1
      integer(ISZ):: ic,ii,id
      real(kind=8):: rr,epsilontmp,normfac,phiav
      real(kind=8):: rhoe(0:nx,0:nzlocal),rhoeote(0:nx,0:nzlocal),denom,Lphi,expo
      real(kind=8)::  boltzfac(0:nx,0:nzlocal)
      real(kind=8):: iondensitytmp,electemptmp,boltzfac1,epszfactmp
      type(CONDtype):: cnd
      integer,save:: printed=0
      real(kind=8):: phiold

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgberz(bnd,nx,nzlocal,phi,mglevel)
      call condbndymgintberz(bnd,nx,nzlocal,dx,dz,phi,bounds,mglevel)
c     print*, "IN RELAXBE, iondensity =", iondensity(0,0), iondensity(4,4)

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
      izmax = nzlocal
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx-1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nzlocal-1

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

c       precompute Boltzmann factor and some other arrays
        do ix = ixmin,ixmax
           do iz = izmin,izmax
              rhoe(ix,iz) = 0.
              rhoeote(ix,iz) = 0.
              if (regionid(ix,iz) > 0) then
                id = regionid(ix,iz)
                call getlocaltemp(ix,iz,dx,dz,id,electemptmp)
                call getlocaldens(ix,iz,dx,dz,id,iondensitytmp)
                call getlocalepsilon(ix,iz,dx,dz,epsilontmp)
                if (iondensitytmp .ne. 0. .and.
     &              electemptmp .ne. 0) then
                  expo=(phi(ix,iz)-plasmapotential(id))/electemptmp
                  if (.not. luselinboltz) then
                     boltzfac(ix,iz) = exp(min(expo,logelecdenmaxscale(id)))
                  else
                     boltzfac(ix,iz) = 1.+ min(expo,logelecdenmaxscale(id))
                  endif
c                  print*,ix,iz, epsilontmp
                  rhoe(ix,iz) = iondensitytmp*boltzfac(ix,iz)/epsilontmp
                  rhoeote(ix,iz) = rhoe(ix,iz)/(electemptmp)
                  ix1 = ixmin + mod(ixmin + iz + parity,2)
c                  if (ix == ix1+2 .and. iz==izmin+1 .and. printed .lt. 10) then
c                   endif
                endif
              endif
           enddo
c          Normalize to z-averaged Boltzmann factor if requestted
           if (lnormtoavboltzfac) then
              normfac=nzlocal/sum(boltzfac(ix,1:nzlocal))
              boltzfac(ix,izmin:izmax) = boltzfac(ix,izmin:izmax)*normfac
              rhoe(ix,izmin:izmax)=rhoe(ix,izmin:izmax)*normfac*
     &         (1.-boltzfac(ix,izmin:izmax)/nzlocal)
              rhoeote(ix,izmin:izmax)=rhoeote(ix,izmin:izmax)*normfac
c           Note this doesn't quite do the conventional thing for linearized
c           Boltzmann, but it saves doing interpolation of Te a second time
c           and should be approximately right whenver linearization is
c           valid
           endif
        enddo

        if (lrz) then
!$OMP DO
          do iz=izmin,izmax

            ix1 = ixmin + mod(ixmin + iz + parity,2)
            do ix=ix1,ixmax,2
              call getlocalepszfac(ix,iz,dx,dz,epszfactmp)
              rr = ix + xminodx
              if (rr == 0.) then
                denom = -4.*dxsqi - 2.*epszfactmp - rhoeote(ix,iz)
                Lphi =  4.*(phi(1,iz) - phi(0,iz))*dxsqi
     &                   + (phi(0,iz-1) - 2.*phi(0,iz) + phi(0,iz+1))*epszfactmp
     &                   - rhoe(ix,iz)
              else
                denom = -2.*dxsqi - 2.*epszfactmp - rhoeote(ix,iz)
                Lphi =  ((rr-0.5)*phi(ix-1,iz) - 2*rr*phi(ix,iz) +
     &                   (rr+0.5)*phi(ix+1,iz))*dxsqi/rr
     &                + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*epszfactmp
     &                - rhoe(ix,iz)
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

              call getlocalepszfac(ix,iz,dx,dz,epszfactmp)
              denom = -2.*(dxsqi+epszfactmp)-rhoeote(ix,iz)
              Lphi =  (phi(ix-1,iz  )+phi(ix+1,iz  ))*dxsqi
     &             +  (phi(ix  ,iz-1)+phi(ix  ,iz+1))*epszfactmp
     &              -  phi(ix,iz)*2.*(dxsqi+epszfactmp) - rhoe(ix,iz)
              phiold = phi(ix,iz)
              phi(ix,iz) = phi(ix,iz) - mgparam*(Lphi - rho(ix,iz))/denom
c              if (ix == ix1+2 .and. iz==izmin+1 .and. printed .lt. 1000) then
c               print*, "printed,denom,Lphi,phiold,phi ", printed,denom,Lphi,phiold,phi(ix,iz)
c               print*, "dxsqi,epszfacdz,rhoeote",dxsqi,epszfactmp,rhoeote(ix,iz)
               printed = printed + 1
c              endif
            enddo
          enddo
c          print*, "summed phi =", sum(phi)
!$OMP END DO
        endif
           

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
         call condbndymgberz(bnd,parity,nx,nzlocal,phi,rho,regionid,
     &                       dxsqi,dzsqi,dx,dz,xminodx,mgparam,bounds,
     &                       mglevel,lrz,icndbndy)
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgberz(bnd,nx,nzlocal,phi,mglevel)
        call condbndymgintberz(bnd,nx,nzlocal,dx,dz,phi,bounds,mglevel)

c       --- set phi in the guard planes
c       --- This must be done inside the loop over parities so that the
c       --- guard planes are updated with the most recent values.
        call applytransversebcberz(nx,nzlocal,phi,bounds)
        call applylongitudinalbcberz(nx,nzlocal,nz,phi,bounds)

c     --- end of loop over even and odd points
      enddo

c     If we are normalizing Boltzman factor to its average, then 
c     we need to fix the absolute value of phi.  For now do this by
c     making the average of phi = 0.  THIS ONLY MAKES SENSE FOR 
c     DOUBLY PERIODIC SYSTEM.
      if (lnormtoavboltzfac) then
         phiav = sum(phi(1:nx,1:nzlocal))/(nx*nzlocal)
         phi = phi - phiav
      endif
c      print*, "end of relaxberez, phi = ", phi(ixmax,izmax),phi(ixmax,izmax-1)
      return
      end
c=============================================================================
      subroutine cond_potmgberz(bnd,nx,nzlocal,phi,mglevel)
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: nx,nzlocal,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1)

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
      subroutine cond_potmgbezerorz(bnd,nx,nzlocal,u,mglevel,delt,delz)
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: nx,nzlocal,mglevel,delt,delz
      real(kind=8):: u(-delt:nx+delt,-delz:nzlocal+delz)

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
      subroutine condbndymgberz(bnd,parity,nx,nzlocal,phi,rho,regionid,
     &                          dxsqi,dzsqi,dx,dz,xminodx,
     &                          mgparam,bounds,mglevel,lrz,icndbndy)
      use BNDtypemodule
      use CONDtypemodule
      use BoltzmannElectrons
      type(BNDtype):: bnd
      integer(ISZ):: parity
      integer(ISZ):: nx,nzlocal,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1), rho(0:nx,0:nzlocal)
      integer(ISZ):: regionid(0:nx,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,dx,dz,xminodx,mgparam
      integer(ISZ):: bounds(0:5),icndbndy
      logical(ISZ):: lrz

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
c The ppp factor (the minimum of the dels at each grid point) is included
c here to be consistent with its inclusion in the Lphi calculation.
c Because of the structure of the solver, terms included in Lphi
c must nearly cancel the same terms here, and so must be scaled
c the same. Note that is it still unknown why the ppp factor is needed
c at all for convergence (the solver diverges with ppp=1). Also note that
c here, ppp=1 for the base level since its rho is the raw rho and does
c not include Lphi.

      type(CONDtype):: cnd
      real(kind=8):: pik,pxm,pzm,pxp,pzp,denom,rhoe,Lphi,ppp,expo
      real(kind=8):: delx,delz,Cx,Cz,rr
      real(kind=8):: delxsqi,delzsqi
      integer(ISZ):: ic,ii,i1,i2
      integer(ISZ):: ix,iz,id
      real(kind=8):: dxi,dzi,dxm,dxp,dzm,dzp
      real(kind=8):: iondensitytmp,electemptmp

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

          dxm = cnd%dxm(ii)*dxi
          dxp = cnd%dxp(ii)*dxi
          dzm = cnd%dzm(ii)*dzi
          dzp = cnd%dzp(ii)*dzi

c         --- the point lower in x is inside the conductor
          if (0. < dxm .and. dxm < 1.) then
            if (lrz) then
              if (rr > 0.) then
                pxm = (1.-0.5/rr)*(cnd%volt0xm(ii) - pik)/dxm
                Cx = Cx + (1.-0.5/rr)*(1./dxm - 1.)
                if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
                ppp = min(ppp,dxm)
              endif
            else
              pxm = (cnd%volt0xm(ii) - pik)/dxm
              Cx = Cx - 1. + 1./dxm
              if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
              ppp = min(ppp,dxm)
            endif
          else if (-1. < dxm .and. dxm <= 0.) then
            pxm = 0.
            Cx = Cx - 1.
            delx = delx - 0.5 + (-dxm)
            if (-dxm > 0.) then
              ppp = min(ppp,-dxm)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- the point higher in x is inside the conductor
          if (0. < dxp .and. dxp < 1.) then
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
          else if (-1. < dxp .and. dxp <= 0.) then
            pxp = 0.
            Cx = Cx - 1.
            delx = delx - 0.5 + (-dxp)
            if (-dxp > 0.) then
              ppp = min(ppp,-dxp)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- the point lower in z is inside the conductor
          if (0. < dzm .and. dzm < 1.) then
            pzm = (cnd%volt0zm(ii) - pik)/dzm
            Cz = Cz - 1. + 1./dzm
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzm
            ppp = min(ppp,dzm)
          else if (-1. < dzm .and. dzm <= 0.) then
            pzm = 0.
            Cz = Cz - 1.
            delz = delz - 0.5 + (-dzm)
            if (-dzm > 0.) then
              ppp = min(ppp,-dzm)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- the point higher in z is inside the conductor
          if (0. < dzp .and. dzp < 1.) then
            pzp = (cnd%volt0zp(ii) - pik)/dzp
            Cz = Cz - 1. + 1./dzp
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzp
            ppp = min(ppp,dzp)
          else if (-1. < dzp .and. dzp <= 0.) then
            pzp = 0.
            Cz = Cz - 1.
            delz = delz - 0.5 + (-dzp)
            if (-dzp > 0.) then
              ppp = min(ppp,-dzp)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- calculate the new phi based on the boundary conditions
          if (ppp < 1.) then
            if (mglevel == 0) ppp = 1.
            delxsqi = dxsqi/dvnz(delx)
            delzsqi = dzsqi/dvnz(delz)
            rhoe = 0.
            denom = - Cx*delxsqi - Cz*delzsqi
            if (regionid(ix,iz) > 0) then
              id = regionid(ix,iz)
              call getlocaldens(ix,iz,dx,dz,id,iondensitytmp)
              call getlocaltemp(ix,iz,dx,dz,id,electemptmp)
              if (iondensitytmp .ne. 0. .and.
     &            electemptmp .ne. 0) then
                expo = (pik - plasmapotential(id))/electemptmp
                expo = min(expo,log(electrondensitymaxscale(id)))
                rhoe = iondensitytmp*exp(expo)
                denom = denom - rhoe/electemptmp
              endif
            endif
            Lphi = (pxm+pxp)*delxsqi + (pzm+pzp)*delzsqi - rhoe
            phi(ix,iz) = pik - mgparam*(Lphi*ppp - rho(ix,iz))/denom
          endif
        enddo
!$OMP END DO
      enddo

      return
      end
c=============================================================================
      subroutine condbndymgintberz(bnd,nx,nzlocal,dx,dz,phi,bounds,mglevel)
      use BNDtypemodule
      use CONDtypemodule
      type(BNDtype):: bnd
      integer(ISZ):: nx,nzlocal,mglevel
      real(kind=8):: dx,dz
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1)
      integer(ISZ):: bounds(0:5)

c Set the potential at points just beyond an internal Neumann boundary.

      type(CONDtype):: cnd
      integer(ISZ):: ic,ii,i1,i2,ix,iz
      real(kind=8):: dxi,dzi,dxm,dxp,dzm,dzp

c     --- These should really be passed in
      dxi = 1./dx
      dzi = 1./dz

c     --- loop over points near surface of conductors
      do ic = 1,bnd%nb_conductors
        if (ic == 1) cnd = bnd%cndfirst
        if (ic > 1) cnd = cnd%next

!$OMP DO
        do ii = 1,cnd%nbbnd

          ix = cnd%jj(ii) - 1
          iz = cnd%kk(ii) - 1

          dxm = cnd%dxm(ii)*dxi
          dxp = cnd%dxp(ii)*dxi
          dzm = cnd%dzm(ii)*dzi
          dzp = cnd%dzp(ii)*dzi

          if (0. > dxm .and. dxm > -1.) phi(ix-1,iz) = phi(ix,iz)
          if (0. > dxp .and. dxp > -1.) phi(ix+1,iz) = phi(ix,iz)
          if (0. > dzm .and. dzm > -1.) phi(ix,iz-1) = phi(ix,iz)
          if (0. > dzp .and. dzp > -1.) phi(ix,iz+1) = phi(ix,iz)

          if (0. == dxm) phi(ix+1,iz) = phi(ix,iz)
          if (0. == dxp) phi(ix-1,iz) = phi(ix,iz)
          if (0. == dzm) phi(ix,iz+1) = phi(ix,iz)
          if (0. == dzp) phi(ix,iz-1) = phi(ix,iz)

        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine Lphiberz(nx,nzlocal,nz,dxsqi,dzsqi,dx,dz,xminodx,
     &                    phi,Lphi,regionid,
     &                    mglevel,bounds,
     &                    lrz,lcndbndy,icndbndy,bnd)
      use BNDtypemodule
      use BoltzmannElectrons
      integer(ISZ):: nx,nzlocal,nz
      real(kind=8):: dxsqi,dzsqi,dx,dz,xminodx
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1)
      real(kind=8):: Lphi(-1:nx+1,-1:nzlocal+1)
      real(kind=8):: boltzfac(0:nx,0:nzlocal)
      real(kind=8):: iondenstmp,normfac,epszfactmp
      integer(ISZ):: regionid(0:nx,0:nzlocal)
      integer(ISZ):: mglevel,bounds(0:5)
      logical(ISZ):: lrz,lcndbndy
      integer(ISZ):: icndbndy
      type(BNDtype):: bnd

c Calculate the Lphi on the grid including the BE term.
c Note that this is note quite the residual since it does not have the
c source term included, and so will no go to zero at convergence.

      real(kind=8):: rhoe(-1:nx+1,-1:nzlocal+1),expo,epsilontmp
      real(kind=8):: rr
      integer(ISZ):: ix,iz,ic,id
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: iondensitytmp,electemptmp

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      izmin = 0
      izmax = nzlocal
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx-1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nzlocal-1

c     --- These are set here to avoid complaints from valgrind about
c     --- uninitialized variables.
      iondensitytmp = 0.
      electemptmp = 0.

c     precompute Boltzmann factor and some other arrays
      do ix = ixmin,ixmax
         do iz = izmin,izmax
            rhoe(ix,iz) = 0.
            if (regionid(ix,iz) > 0) then
               id = regionid(ix,iz)
               call getlocaltemp(ix,iz,dx,dz,id,electemptmp)
               call getlocaldens(ix,iz,dx,dz,id,iondenstmp)
               call getlocalepsilon(ix,iz,dx,dz,epsilontmp)
               if (iondensitytmp .ne. 0. .and.
     &              electemptmp .ne. 0) then
                  expo=(phi(ix,iz)-plasmapotential(id))/electemptmp
                  if (.not. luselinboltz) then
                     boltzfac(ix,iz) = exp(min(expo,logelecdenmaxscale(id)))
                  else
                     boltzfac(ix,iz) = 1.+ min(expo,logelecdenmaxscale(id))
                  endif
                  rhoe(ix,iz) = iondensitytmp*boltzfac(ix,iz)/epsilontmp
               endif
            endif
         enddo
c          Normalize to z-averaged Boltzmann factor if requestted
         if (lnormtoavboltzfac) then
            normfac=nzlocal/sum(boltzfac(ix,1:nzlocal))
            boltzfac(ix,izmin:izmax) = boltzfac(ix,izmin:izmax)*normfac
            rhoe(ix,izmin:izmax)=rhoe(ix,izmin:izmax)*normfac*
     &         (1.-boltzfac(ix,izmin:izmax)/nzlocal)
c           Note this doesn't quite do the conventional thing for linearized
c           Boltzmann, but it saves doing interpolation of Te a second time
c           and should be approximately right whenver linearization is
c           valid
         endif
      enddo

      Lphi = 0.
c     --- Calculate the Lphi.
      if (lrz) then
c       --- RZ
!$OMP DO
        do iz=izmin,izmax

          do ix=ixmin,ixmax
            call getlocalepszfac(ix,iz,dx,dz,epszfactmp)
            rr = ix + xminodx
            if (rr == 0.) then
              Lphi(ix,iz) = 4*(phi(1,iz) - phi(0,iz))*dxsqi
     &                   + (phi(0,iz-1) - 2.*phi(0,iz) + phi(0,iz+1))*epszfactmp
     &                   - rhoe(ix,iz)
            else
              Lphi(ix,iz) = ((rr-0.5)*phi(ix-1,iz) - 2*rr*phi(ix,iz) +
     &                      (rr+0.5)*phi(ix+1,iz))*dxsqi/rr
     &                   + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*epszfactmp
     &                   - rhoe(ix,iz)
            endif

          enddo
        enddo
!$OMP END DO
      else
c       --- XZ
!$OMP DO
        do iz=izmin,izmax

          do ix=ixmin,ixmax

            call getlocalepszfac(ix,iz,dx,dz,epszfactmp)

            Lphi(ix,iz) = (phi(ix-1,iz) - 2.*phi(ix,iz) + phi(ix+1,iz))*dxsqi
     &                  + (phi(ix,iz-1) - 2.*phi(ix,iz) + phi(ix,iz+1))*epszfactmp
     &                  - rhoe(ix,iz)

          enddo
        enddo
!$OMP END DO
      endif

c     --- Zero the Lphi inside conductors.
      call cond_potmgbezerorz(bnd,nx,nzlocal,Lphi,mglevel,1,1)

      if (lcndbndy) then
c       --- Calculate the Lphi near the conductor.
        call condbndyLphiberz(bnd,nx,nzlocal,phi,Lphi,regionid,dxsqi,dzsqi,
     &                        dx,dz,xminodx,bounds,mglevel,lrz,icndbndy)
      endif

c     --- Transverse boundaries
      if (bounds(0) == 0) Lphi(-1:0,:) = 0.
      if (bounds(1) == 0) Lphi(nx:nx+1,:) = 0.
      if (bounds(0) == 1) Lphi(-1,izmin:izmax) = Lphi(1,izmin:izmax)
      if (bounds(1) == 1) Lphi(nx+1,izmin:izmax) = Lphi(nx-1,izmin:izmax)
      if (bounds(0) == 2) Lphi(-1,izmin:izmax) = Lphi(nx-1,izmin:izmax)
      if (bounds(1) == 2) Lphi(nx:nx+1,izmin:izmax) = Lphi(0:1,izmin:izmax)

c     --- Longitudinal boundaries
      if (bounds(4) == 0) Lphi(:,-1:0) = 0.
      if (bounds(5) == 0) Lphi(:,nzlocal:nzlocal+1) = 0.
      if (bounds(4) == 1) Lphi(:,-1) = Lphi(:,1)
      if (bounds(5) == 1) Lphi(:,nzlocal+1) = Lphi(:,nzlocal-1)
      if (bounds(4) == 2 .and. nzlocal == nz) Lphi(:,-1) = Lphi(:,nzlocal-1)
      if (bounds(5) == 2 .and. nzlocal == nz) Lphi(:,nzlocal+1) = Lphi(:,1)

      return
      end
c=============================================================================
      subroutine condbndyLphiberz(bnd,nx,nzlocal,phi,Lphi,regionid,dxsqi,dzsqi,
     &                            dx,dz,xminodx,bounds,mglevel,lrz,icndbndy)
      use BNDtypemodule
      use CONDtypemodule
      use BoltzmannElectrons
      type(BNDtype):: bnd
      integer(ISZ):: nx,nzlocal,mglevel
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1), Lphi(-1:nx+1,-1:nzlocal+1)
      integer(ISZ):: regionid(0:nx,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,dx,dz,xminodx
      integer(ISZ):: bounds(0:5),icndbndy
      logical(ISZ):: lrz

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

      type(CONDtype):: cnd
      real(kind=8):: pik,pxm,pzm,pxp,pzp,denom,rhoe,ppp,expo
      real(kind=8):: delx,delz
      real(kind=8):: delxsqi,delzsqi
      real(kind=8):: rr
      integer(ISZ):: ic,ii
      integer(ISZ):: ix,iz,id
      real(kind=8):: dxi,dzi,dxm,dxp,dzm,dzp
      real(kind=8):: iondensitytmp,electemptmp

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
          ppp = 1.

          dxm = cnd%dxm(ii)*dxi
          dxp = cnd%dxp(ii)*dxi
          dzm = cnd%dzm(ii)*dzi
          dzp = cnd%dzp(ii)*dzi

c         --- the point lower in x is inside the conductor
          if (0. < dxm .and. dxm < 1.) then
            if (lrz) then
              if (rr > 0.) then
                pxm = (1.-0.5/rr)*(cnd%volt0xm(ii) - pik)/dxm
                if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
                ppp = min(ppp,dxm)
              endif
            else
              pxm = (cnd%volt0xm(ii) - pik)/dxm
              if (icndbndy == 2) delx = delx - 0.5 + 0.5*dxm
              ppp = min(ppp,dxm)
            endif
          else if (-1. < dxm .and. dxm <= 0.) then
            pxm = 0.
            delx = delx - 0.5 + (-dxm)
            if (-dxm > 0.) then
              ppp = min(ppp,-dxm)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- the point higher in x is inside the conductor
          if (0. < dxp .and. dxp < 1.) then
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
          else if (-1. < dxp .and. dxp <= 0.) then
            pxp = 0.
            delx = delx - 0.5 + (-dxp)
            if (-dxp > 0.) then
              ppp = min(ppp,-dxp)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- the point lower in z is inside the conductor
          if (0. < dzm .and. dzm < 1.) then
            pzm = (cnd%volt0zm(ii) - pik)/dzm
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzm
            ppp = min(ppp,dzm)
          else if (-1. < dzm .and. dzm <= 0.) then
            pzm = 0.
            delz = delz - 0.5 + (-dzm)
            if (-dzm > 0.) then
              ppp = min(ppp,-dzm)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- the point higher in z is inside the conductor
          if (0. < dzp .and. dzp < 1.) then
            pzp = (cnd%volt0zp(ii) - pik)/dzp
            if (icndbndy == 2) delz = delz - 0.5 + 0.5*dzp
            ppp = min(ppp,dzp)
          else if (-1. < dzp .and. dzp <= 0.) then
            pzp = 0.
            delz = delz - 0.5 + (-dzp)
            if (-dzp > 0.) then
              ppp = min(ppp,-dzp)
            else
              ppp = min(ppp,1.-1.e-9)
            endif
          endif
c         --- calculate the new phi based on the boundary conditions
          if (ppp < 1.) then
            delxsqi = dxsqi/dvnz(delx)
            delzsqi = dzsqi/dvnz(delz)
            rhoe = 0.
            if (regionid(ix,iz) > 0) then
              id = regionid(ix,iz)
              call getlocaldens(ix,iz,dx,dz,id,iondensitytmp)
              call getlocaltemp(ix,iz,dx,dz,id,electemptmp)
              if (iondensitytmp .ne. 0. .and.
     &            electemptmp .ne. 0) then
                expo = (pik - plasmapotential(id))/electemptmp
                expo = min(expo,log(electrondensitymaxscale(id)))
                rhoe = iondensitytmp*exp(expo)
              endif
            endif
            Lphi(ix,iz) = (+ (pxm+pxp)*delxsqi
     &                     + (pzm+pzp)*delzsqi
     &                     - rhoe)*ppp

            if (0. > dxm .and. dxm >= -1.) Lphi(ix-1,iz  ) = 0.
            if (0. > dxp .and. dxp >= -1.) Lphi(ix+1,iz  ) = 0.
            if (0. > dzm .and. dzm >= -1.) Lphi(ix  ,iz-1) = 0.
            if (0. > dzp .and. dzp >= -1.) Lphi(ix  ,iz+1) = 0.

            if (0. == dxm) Lphi(ix+1,iz  ) = 0.
            if (0. == dxp) Lphi(ix-1,iz  ) = 0.
            if (0. == dzm) Lphi(ix  ,iz+1) = 0.
            if (0. == dzp) Lphi(ix  ,iz-1) = 0.

          endif
        enddo
!$OMP END DO
      enddo

      return
      end
c=============================================================================
      subroutine applytransversebcberz(nx,nzlocal,u,bounds)
      integer(ISZ):: nx,nzlocal
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:nzlocal+1)
      if (bounds(0) == 0) u(-1,:)   = 2.*u(0,:) - u(1,:)
      if (bounds(1) == 0) u(nx+1,:) = 2.*u(nx,:) - u(nx-1,:)
      if (bounds(0) == 1) u(-1,:)   = u(1,:)
      if (bounds(1) == 1) u(nx+1,:) = u(nx-1,:)
      if (bounds(0) == 2) u(-1,:)   = u(nx-1,:)
      if (bounds(1) == 2) u(nx+1,:) = u(1,:)
      return
      end
c=============================================================================
      subroutine applylongitudinalbcberz(nx,nzlocal,nz,u,bounds)
      integer(ISZ):: nx,nzlocal,nz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-1:nx+1,-1:nzlocal+1)
      if (bounds(4) == 0) u(:,-1)   = 2.*u(:,0) - u(:,1)
      if (bounds(5) == 0) u(:,nzlocal+1) = 2.*u(:,nzlocal) - u(:,nzlocal-1)
      if (bounds(4) == 1) u(:,-1)   = u(:,1)
      if (bounds(5) == 1) u(:,nzlocal+1) = u(:,nzlocal-1)
      if (nzlocal == nz) then
        if (bounds(4) == 2) u(:,-1)   = u(:,nzlocal-1)
        if (bounds(5) == 2) u(:,nzlocal+1) = u(:,1)
      endif
      return
      end
c=============================================================================
      subroutine clampphitophimaxberz(nx,nzlocal,dx,dz,phi,rho,regionid)
      use BoltzmannElectrons
      integer(ISZ):: nx,nzlocal
      real(kind=8):: dx,dz
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1),rho(0:nx,0:nzlocal)
      integer(ISZ):: regionid(0:nx,0:nzlocal)

      integer(ISZ):: ix,iz,id
      real(kind=8):: expo,phimax,rhomax
      real(kind=8):: iondensitytmp,electemptmp

c     --- Clamp phi to phimax. The phimax is calculated assuming that
c     --- the max electron density will be the max of the iondensity and 
c     --- and the ion particle density (rho).
      do iz=0,nzlocal
        do ix=0,nx
          if (regionid(ix,iz) > 0) then
            id = regionid(ix,iz)
            call getlocaldens(ix,iz,dx,dz,id,iondensitytmp)
            call getlocaltemp(ix,iz,dx,dz,id,electemptmp)
            if (iondensitytmp == 0.) cycle
            if (luseparticledensity(id) > 0) then
              rhomax = max(electrondensitymaxscale(id),-rho(ix,iz)/iondensitytmp)
            else
              rhomax = electrondensitymaxscale(id)
            endif
            phimax = plasmapotential(id)+electemptmp*log(rhomax)
            if (phi(ix,iz) > phimax) phi(ix,iz) = phimax
          endif
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine setupregionidsbe2d(regionid,nx,nzlocal,dx,dz,xmmin,zmmin)
      use BoltzmannElectrons
      integer(ISZ):: nx,nzlocal
      integer(ISZ):: regionid(0:nx,0:nzlocal)
      real(kind=8):: dx,dy,dz,xmmin,zmmin

c Sets id of regions where the Boltzmann-electron distribution is used.

      integer(ISZ):: id,ixmin,ixmax,izmin,izmax,ix,iz

      regionid = 0
      do id=1,nberegions
        ixmin = (max(0.,xbemin(id) - xmmin))/dx
        ixmax = min(1.0*nx,(xbemax(id) - xmmin)/dx)
        izmin = (max(0.,zbemin(id) - zmmin))/dz
        izmax = min(1.0*nzlocal,(zbemax(id) - zmmin)/dz)
        regionid(ixmin:ixmax,izmin:izmax) = id
      enddo

      return
      end
c=============================================================================
      subroutine getlocaldens(ix,iz,dx,dz,id,iondensitytmp)
c     generate interpolated local iondensity and electron temperature
      use BoltzmannElectrons
      integer(ISZ):: ix,iz,id
      real(kind=8):: dx,dz,reixx,reizz
      real(kind=8):: iondensitytmp

      real(kind=8):: xx,zz,wxx,wzz
      integer(ISZ):: ixx,izz

c      print*,"ix,iz,id,liondensitygrid2d ", ix,iz,id,liondensitygrid2d
      iondensitytmp = 0.
      if (liondensitygrid2d(id)) then
        xx = ix*dx
        zz = iz*dz
        reixx = (xx - iondensitygrid2d%xmin)/iondensitygrid2d%dx
        reizz = (zz - iondensitygrid2d%ymin)/iondensitygrid2d%dy
        ixx = reixx
        izz = reizz
        wxx = reixx-ixx
        wzz = reizz-izz
        if (ixx+1 > iondensitygrid2d%nx) then
          ixx = iondensitygrid2d%nx-1
          wxx = 1.
        endif
        if (izz+1 > iondensitygrid2d%ny) then
          izz = iondensitygrid2d%ny-1
          wzz = 1.
        endif
        iondensitytmp =
     &       iondensitygrid2d%grid(ixx  ,izz  )*(1.-wxx)*(1.-wzz) +
     &       iondensitygrid2d%grid(ixx+1,izz  )*    wxx *(1.-wzz) +
     &       iondensitygrid2d%grid(ixx  ,izz+1)*(1.-wxx)*    wzz  +
     &       iondensitygrid2d%grid(ixx+1,izz+1)*    wxx *    wzz
      else if (iondensity(id) .ne. 0.) then
        iondensitytmp = iondensity(id)
      endif

      return
      end
c=============================================================================
      subroutine getlocaltemp(ix,iz,dx,dz,id,electemptmp)
c     generate interpolated local electemp and electron temperature
      use BoltzmannElectrons
      integer(ISZ):: ix,iz,id
      real(kind=8):: dx,dz,reixx,reizz
      real(kind=8):: electemptmp

      real(kind=8):: xx,zz,wxx,wzz
      integer(ISZ):: ixx,izz

c      print*,"ix,iz,id,electemperature ", ix,iz,id,electrontemperature(id)
      electemptmp = 0.
      if (lelectrontemperaturegrid2d(id)) then
        xx = ix*dx
        zz = iz*dz
        reixx = (xx - electrontemperaturegrid2d%xmin)/electrontemperaturegrid2d%dx
        reizz = (zz - electrontemperaturegrid2d%ymin)/electrontemperaturegrid2d%dy
        ixx = reixx
        izz = reizz
        wxx = reixx-ixx
        wzz = reizz-izz
        if (ixx+1 > electrontemperaturegrid2d%nx) then
           ixx = electrontemperaturegrid2d%nx-1
           wxx = 1.
        endif
        if (izz+1 > electrontemperaturegrid2d%ny) then
           izz = electrontemperaturegrid2d%ny-1
           wzz = 1.
        endif
        electemptmp =
     &       electrontemperaturegrid2d%grid(ixx  ,izz  )*(1.-wxx)*(1.-wzz) +
     &       electrontemperaturegrid2d%grid(ixx+1,izz  )*    wxx *(1.-wzz) +
     &       electrontemperaturegrid2d%grid(ixx  ,izz+1)*(1.-wxx)*    wzz  +
     &       electrontemperaturegrid2d%grid(ixx+1,izz+1)*    wxx *    wzz

      else if (electrontemperature(id) .ne. 0.) then
        electemptmp = electrontemperature(id)
      endif

      return
      end
c=============================================================================
      subroutine getlocalepsilon(ix,iz,dx,dz,epsilontmp)
c     generate interpolated local epsilon2d and electron temperature
      use DKInterp, only:interpdk,icalceps,epsilon2d,iusing2deps
      use BoltzmannElectrons
      integer(ISZ):: ix,iz,id
      real(kind=8):: dx,dz,reixx,reizz
      real(kind=8):: epsilontmp

      real(kind=8):: xx,zz,wxx,wzz
      integer(ISZ):: ixx,izz

c      print*,"ix,iz,id,lepsilongrid2d ", ix,iz,id,lepsilon2d2d
      epsilontmp = 0.
      if (iusing2deps>0) then
        xx = ix*dx
        zz = iz*dz
        reixx = (xx - epsilon2d%xmin)/epsilon2d%dx
        reizz = (zz - epsilon2d%ymin)/epsilon2d%dy
        ixx = reixx
        izz = reizz
        wxx = reixx-ixx
        wzz = reizz-izz
        if (ixx+1 > epsilon2d%nx) then
           ixx = epsilon2d%nx-1
           wxx = 1.
        endif
        if (izz+1 > epsilon2d%ny) then
           izz = epsilon2d%ny-1
           wzz = 1.
        endif
        epsilontmp =
     &       epsilon2d%grid(ixx  ,izz  )*(1.-wxx)*(1.-wzz) +
     &       epsilon2d%grid(ixx+1,izz  )*    wxx *(1.-wzz) +
     &       epsilon2d%grid(ixx  ,izz+1)*(1.-wxx)*    wzz  +
     &       epsilon2d%grid(ixx+1,izz+1)*    wxx *    wzz
      else
        epsilontmp = 1.
      endif

      return
      end

c=============================================================================
      subroutine getlocalepszfac(ix,iz,dx,dz,epszfactmp)
c     generate interpolated local epsilon2d and electron temperature
      use DKInterp, only:interpdk,icalceps,epszfac2d,iusing2deps
      use BoltzmannElectrons
      integer(ISZ):: ix,iz,id
      real(kind=8):: dx,dz,reixx,reizz
      real(kind=8):: epszfactmp

      real(kind=8):: xx,zz,wxx,wzz
      integer(ISZ):: ixx,izz

c      print*,"ix,iz,id,iusing2deps ", ix,iz,id,lepsilon2d2d
      epszfactmp = 0.
      if (iusing2deps>0) then
        xx = ix*dx
        zz = iz*dz
        reixx = (xx - epszfac2d%xmin)/epszfac2d%dx
        reizz = (zz - epszfac2d%ymin)/epszfac2d%dy
        ixx = reixx
        izz = reizz
        wxx = reixx-ixx
        wzz = reizz-izz
        if (ixx+1 > epszfac2d%nx) then 
           ixx = epszfac2d%nx-1
           wxx = 1.
        endif
        if (izz+1 > epszfac2d%ny) then 
           izz = epszfac2d%ny-1
           wzz = 1.
        endif
        epszfactmp =
     &       epszfac2d%grid(ixx  ,izz  )*(1.-wxx)*(1.-wzz) +
     &       epszfac2d%grid(ixx+1,izz  )*    wxx *(1.-wzz) +
     &       epszfac2d%grid(ixx  ,izz+1)*(1.-wxx)*    wzz  +
     &       epszfac2d%grid(ixx+1,izz+1)*    wxx *    wzz
      else
        epszfactmp = 1.
      endif
      epszfactmp = epszfactmp/(dz*dz)

      return
      end
c=============================================================================
