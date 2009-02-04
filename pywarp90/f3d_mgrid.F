#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID.M, version $Revision: 3.185 $, $Date: 2009/02/04 00:42:10 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid field solver which is part of the F3D
c  package of WARP.
c  David P. Grote, LLNL, (925)423-7194, LBNL (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigrid3df(iwhich,nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                        dx,dy,dz,phi,rho,rstar,linbend,
     &                        bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                        xmmin,ymmin,zmmin,zbeam,zgrid)
      use GlobalVars
      use Conductor3d
      use Multigrid3d
      use Parallel
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rstar(-1:nzlocal+1)
      logical(ISZ):: linbend
      integer(ISZ):: bound0,boundnz,boundxy
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid

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

      call multigrid3dsolve(iwhich,nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                      phi,rho,rstar,linbend,bounds,
     &                      xmmin,ymmin,zmmin,zbeam,zgrid,
     &                      mgparam,mgform,mgiters,mgmaxiters,
     &                      mgmaxlevels,mgerror,mgtol,mgverbose,
     &                      downpasses,uppasses,
     &                      lcndbndy,laddconductor,icndbndy,
     &                      gridmode,conductors,lprecalccoeffs,fsdecomp)

      return
      end
c=============================================================================
      subroutine multigrid3dsolve(iwhich,nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                            dx,dy,dz,phi,rho,
     &                            rstar,linbend,bounds,
     &                            xmmin,ymmin,zmmin,zbeam,zgrid,
     &                            mgparam,mgform,mgiters,mgmaxiters,
     &                            mgmaxlevels,mgerror,mgtol,mgverbose,
     &                            downpasses,uppasses,
     &                            lcndbndy,laddconductor,icndbndy,
     &                            gridmode,conductors,lprecalccoeffs,fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      use Decompositionmodule
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nzlocal+1)
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgform,mgiters,mgmaxiters,mgmaxlevels,mgverbose
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor,lprecalccoeffs
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

c Use the multigrid method for solving Poisson's equation on a 3-D Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution.
c
c When the grid cells are rectangular, semi-coarsening is done until the
c grid cell dimensions are roughly equal. Roughly equal means that
c   2/3 dx < dz < 4/3 dx
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3. Currently, it is still
c assumed that dx ~ dy and that semi-coarsening is not needed transversely.
c
c The first call to vcycle can be done using one of two forms. When mgform
c is 1, the normal form is used and phi and rho are passed directly into
c vcycle. When mgform is 2, the error and the residual are passed in instead.
c The two produce nearly identical results and there is no effect on
c convergence. The second form, residual correction form, was put in to be
c consistent with the Chombo AMR/MG field solver.

c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
      real(kind=8):: maxres
#endif
      real(kind=8):: dxsqi,dysqi,dzsqi,reps0c,rdel
      integer(ISZ):: i,ii,k,ix,iy,iz
      real(kind=8):: rs,x,r
      real(kind=8),allocatable:: phisave(:,:,:)
      real(kind=8):: bendx(-1:nxlocal+1)
c     --- The following only used when mgform == 2
      integer(ISZ):: delx,dely,delz
      real(kind=8),allocatable:: rhosave(:,:,:),res(:,:,:)
      integer(ISZ):: localbounds(0:5)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Note that nx and ny must be even.
      if (mod(nx,2) == 1 .or. mod(ny,2) == 1) then
        call kaboom("3D MULTIGRID field solver error, nx and ny must be even")
        return
      endif

c     --- Initialize temporaries
      dxsqi  = 1./dx**2
      dysqi  = 1./dy**2
      dzsqi  = 1./dz**2
      reps0c = mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
      rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

c     --- Setup the boundary conditions for the local domain.
      localbounds = bounds
#ifdef MPIPARALLEL
      if (fsdecomp%ix(fsdecomp%ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(fsdecomp%ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iy(fsdecomp%iyproc) > 0)          localbounds(2) = -1
      if (fsdecomp%iy(fsdecomp%iyproc)+nylocal < ny) localbounds(3) = -1
      if (fsdecomp%iz(fsdecomp%izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(fsdecomp%izproc)+nzlocal < nz) localbounds(5) = -1
#endif

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
      call checkconductors(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                     conductors,fsdecomp)

      if (lprecalccoeffs) then
        call precalculatecoefficients3d(nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                                  conductors,mgparam,localbounds)
!$OMP DO
        do iz=0,nzlocal
          do iy=0,nylocal
            do ix=0,nxlocal
              if (conductors%coeffs%data(7,ix,iy,iz) == 0.) then
                phi(ix,iy,iz) = conductors%coeffs%data(6,ix,iy,iz)
              endif
            enddo
          enddo
        enddo
!$OMP END DO
      endif

!$OMP PARALLEL
!$OMP&PRIVATE(ii,i,k,rs,x,r,ix,iy,iz)

c     --- Preset rho to increase performance (reducing the number of
c     --- multiplies in the main SOR sweep loop).
      if (.not. linbend) then
!OMP DO
        rho = rho*reps0c
!OMP END DO
      else
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
!$OMP DO
        do iz=0,nzlocal
          rs = rstar(iz)
          do ix=0,nxlocal
            x  = xmmin + (ix + fsdecomp%ix(fsdecomp%ixproc))*dx
            r  = rs + x
c           --- rearranged to reduce divides
c           --- rho(ix,:,iz) = rho(ix,:,iz)*(rs/r)*reps0c/
c           ---             ( 1. + (x/r)*((x/r)-2.)*rdel )
            rho(ix,:,iz) = rho(ix,:,iz)*reps0c*rs*r/(r*r + x*(x-2.*r)*rdel)
          enddo
        enddo
!$OMP END DO
c       --- Fill scratch array with x values so it can be looked up
c       --- in the bent beam loop instead of calculated.
!$OMP DO
        do ix = -1,nxlocal+1
          bendx(ix) = xmmin + (ix + fsdecomp%ix(fsdecomp%ixproc))*dx
        enddo
!$OMP END DO
c       --- Change rstar if using Nuemann boundary conditions
        if (bounds(4) == 1) rstar(-1) = rstar(1)
        if (bounds(5) == 1) rstar(nzlocal+1) = rstar(nzlocal-1)
      endif

#ifdef MPIPARALLEL
c     --- These calls break the parallel field solver
c     call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,phi,
c    &                             1,1,1,-1,0,localbounds,fsdecomp)
c     call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,-1,fsdecomp)
#endif

c     --- Make sure guard planes have sensible values before beginning.
      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                               localbounds,.false.,.false.)

c     --- If using residual correction form, need to save the original rho.
c     --- Also setup parallel arrays.
      if (mgform == 2) then
!$OMP SINGLE
        allocate(rhosave(0:nxlocal,0:nylocal,0:nzlocal),stat=allocerror)
        if (allocerror /= 0) then
          print*,"multigrid3dsolve: allocation error ",allocerror,
     &           ": could not allocate rhosave to shape ",
     &           nxlocal,nylocal,nzlocal
          call kaboom("multigrid3dsolve: allocation error")
          return
        endif

        delx = 1
        dely = 1
        delz = 1
        if (nx > nxlocal) delx = 3
        if (ny > nylocal) dely = 3
        if (nz > nzlocal) delz = 3

        allocate(res(-delx:nxlocal+delx,
     &               -dely:nylocal+dely,
     &               -delz:nzlocal+delz),stat=allocerror)
        if (allocerror /= 0) then
          print*,"multigrid3dsolve: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nxlocal,nylocal,nzlocal
          call kaboom("multigrid3dsolve: allocation error")
          return
        endif
!$OMP END SINGLE
        rhosave = rho
      endif

      allocate(phisave(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nxlocal,nylocal,nzlocal
        call kaboom("multigrid3dsolve: allocation error")
        return
      endif

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- very few iterations are done.
      mgiters = 0
      mgerror = 2.*mgtol + 1.
      do while (mgerror > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phi

c       --- If using residual correction form, calculate the residual and
c       --- copy it into rhosave, zero phisave (the initial error).
c       --- In the calls to cond_potmg and residual, the last argument
c       --- is true, telling the routines to use the actual value of
c       --- voltages rather than zero as is done otherwise for residual
c       --- correction form since it is operating on the error.
        if (mgform == 2) then
          if (.not.lprecalccoeffs) then
            call cond_potmg(conductors%interior,nxlocal,nylocal,nzlocal,1,1,1,
     &                      phisave,0,mgform,.true.)
          endif
          call condbndymgint(conductors%evensubgrid,nxlocal,nylocal,nzlocal,
     &                       1,1,1,phisave,localbounds,0,icndbndy)
          call condbndymgint(conductors%oddsubgrid,nxlocal,nylocal,nzlocal,
     &                       1,1,1,phisave,localbounds,0,icndbndy)
          call residual3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                    dxsqi,dysqi,dzsqi,phisave,rhosave,res,
     &                    0,localbounds,mgparam,mgform,.true.,
     &                    lcndbndy,icndbndy,conductors,lprecalccoeffs)
#ifdef MPIPARALLEL
          call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,res,
     &                                 delx,dely,delz,-1,0,localbounds,fsdecomp)
          call mgexchange_phi(1,nxlocal,nylocal,nzlocal,res,
     &                        delx,dely,delz,-1,-1,fsdecomp)
#endif
          rho = res(0:nxlocal,0:nylocal,0:nzlocal)
          phi = 0.
        endif

c       --- Do one vcycle.
        call vcycle(0,1.,nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,phi,rho,
     &              rstar,linbend,bendx,bounds,mgparam,mgform,mgmaxlevels,
     &              downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &              lprecalccoeffs,fsdecomp)

c       --- If using residual correction form, add the resulting error to phi.
        if (mgform == 2) phi = phi + phisave

c       --- When using residual correction form, the other planes do need
c       --- to be set when using other than Dirichlet boundaries since
c       --- those planes are only set with the error of phi.
        if (mgform == 2) then
          call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                                   localbounds,.false.,.false.)
#ifdef MPIPARALLEL
          call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,0,fsdecomp)
#endif
        endif

c       --- Calculate the change in phi.
        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nzlocal
          do iy=0,nylocal
            do ix=0,nxlocal
              mgerror = max(mgerror,abs(phisave(ix,iy,iz) - phi(ix,iy,iz)))
            enddo
          enddo
        enddo
!$OMP END DO

c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
        maxres = 0.
        if (mgform == 2) then
!$OMP DO REDUCTION(MAX:maxres)
           do i=0,nzlocal
              maxres = max(maxres, maxval(abs(res(:,:,iz))))
           enddo
!$OMP END DO
        endif
#endif

#ifdef MPIPARALLEL
        if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
c         --- calculate global sorerror
          call parallelmaxrealarray(mgerror,1)
c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
          if (mgform == 2) then
             call parallelmaxrealarray(maxres,1)
          endif
#endif
        endif
#endif
c       print*,mgiters,mgerror

c       --- This line below seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop above.
c       mgerror = maxval(abs(phisave - phi))

c       added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
        if (mgverbose>=1) then
          print *, 'iteration ', mgiters, ' mgerror=', mgerror
c       unscaled residual
          print *, 'max(residual)=', maxres*2.0*(dxsqi+dysqi+dzsqi)
        endif
#endif

      enddo

c     --- Set boundary conditions. This is only really needed for the
c     --- Dirichlet boundaries, but this is convenient to call.
      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                               localbounds,.true.,.false.)

c     --- Make a print out.
      if (mgverbose>=1 .or. (mgverbose>=0 .and. mgerror > mgtol)) then
        if (mgerror > mgtol) then
          call remark("Multigrid: Maximum number of iterations reached")
        endif
        write(errline,20) mgerror,mgiters
  20    format("Multigrid: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif
      
c     --- If using residual correction form, restore saved rho
      if (mgform == 2) then
        rho = rhosave
        deallocate(rhosave,res)
      endif

      deallocate(phisave)

c     --- Restore rho
      reps0c = 1./reps0c
      if (.not. linbend) then
        rho = rho*reps0c
      else
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
        do iz=0,nzlocal
          rs = rstar(iz)
          do ix=0,nxlocal
            x  = xmmin + (ix + fsdecomp%ix(fsdecomp%ixproc))*dx
            r  = rs + x
            rho(ix,:,iz) = rho(ix,:,iz)/rs*reps0c*( r + x*((x/r)-2.)*rdel )
          enddo
        enddo
      endif

!$OMP END PARALLEL

      if (lf3dtimesubs) timemultigrid3dsolve = timemultigrid3dsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      module formggetarraysuminterface
      contains
      function mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,a,fsdecomp)
      use Decompositionmodule
      real(kind=8),dimension(2):: mggetarraysum
      integer(ISZ):: nxlocal,nylocal,nzlocal
      integer(ISZ):: delx,dely,delz
      real(kind=8):: a(-delx:nxlocal+delx,-dely:nylocal+dely,-delz:nzlocal+delz)
      type(Decomposition):: fsdecomp

#ifdef MPIPARALLEL

      integer(ISZ):: ix1,ix2
      integer(ISZ):: iy1,iy2
      integer(ISZ):: iz1,iz2
      real(kind=8):: sss(2)

      ix1 = -delx
      if (fsdecomp%ixproc < fsdecomp%nxprocs-1) then
        ix2 = fsdecomp%ix(fsdecomp%ixproc+1)
     &        - fsdecomp%ix(fsdecomp%ixproc) - 1 - delx
      else
        ix2 = nxlocal + delx
      endif
      iy1 = -dely
      if (fsdecomp%iyproc < fsdecomp%nyprocs-1) then
        iy2 = fsdecomp%iy(fsdecomp%iyproc+1)
     &        - fsdecomp%iy(fsdecomp%iyproc) - 1 - dely
      else
        iy2 = nylocal + dely
      endif
      iz1 = -delz
      if (fsdecomp%izproc < fsdecomp%nzprocs-1) then
        iz2 = fsdecomp%iz(fsdecomp%izproc+1)
     &        - fsdecomp%iz(fsdecomp%izproc) - 1 - delz
      else
        iz2 = nzlocal + delz
      endif

      if (ix2 >= ix1 .and. iy2 >= iy1 .and. iz2 >= iz1) then
        sss(1) = sum(sum(sum(abs(a(ix1:ix2,iy1:iy2,iz1:iz2)),1),1),1)
      else
        sss(1) = 0.
      endif

      if (fsdecomp%ixproc > 0) then
        ix1 = fsdecomp%ix(fsdecomp%ixproc-1)
     &        + fsdecomp%nx(fsdecomp%ixproc-1)
     &        - fsdecomp%ix(fsdecomp%ixproc) + 1 + delx
      else
        ix1 = 0 - delx
      endif
      ix2 = nxlocal + delx
      if (fsdecomp%iyproc > 0) then
        iy1 = fsdecomp%iy(fsdecomp%iyproc-1)
     &        + fsdecomp%ny(fsdecomp%iyproc-1)
     &        - fsdecomp%iy(fsdecomp%iyproc) + 1 + dely
      else
        iy1 = 0 - dely
      endif
      iy2 = nylocal + dely
      if (fsdecomp%izproc > 0) then
        iz1 = fsdecomp%iz(fsdecomp%izproc-1)
     &        + fsdecomp%nz(fsdecomp%izproc-1)
     &        - fsdecomp%iz(fsdecomp%izproc) + 1 + delz
      else
        iz1 = 0 - delz
      endif
      iz2 = nzlocal + delz
      if (ix2 >= ix1 .and. iy2 >= iy1 .and. iz2 >= iz1) then
        sss(2) = sum(sum(sum(abs(a(ix1:ix2,iy1:iy2,iz1:iz2)),1),1),1)
      else
        sss(2) = 0.
      endif

      if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
        call parallelsumrealarray(sss,2)
      endif
      mggetarraysum = sss

#else

      mggetarraysum(1) = sum(abs(a))
      mggetarraysum(2) = mggetarraysum(1)

#endif

      return
      end function mggetarraysum
      end module formggetarraysuminterface
c=============================================================================
      RECURSIVE subroutine vcycle(mglevel,mgscale,nx,ny,nz,
     &                            nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                            phi,rho,rstar,linbend,bendx,globalbounds,
     &                            mgparam,mgform,
     &                            mgmaxlevels,downpasses,uppasses,
     &                            lcndbndy,icndbndy,conductors,lprecalccoeffs,
     &                            fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use Decompositionmodule
      use formggetarraysuminterface
      integer(ISZ):: mglevel
      real(kind=8):: mgscale
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rstar(-1:nzlocal+1)
      real(kind=8):: bendx(-1:nxlocal+1)
      logical(ISZ):: linbend
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lprecalccoeffs
      integer(ISZ):: icndbndy
      type(Decomposition):: fsdecomp

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: mingridsize
      real(kind=8),allocatable:: phicoarse(:,:,:),rhocoarse(:,:,:)
      real(kind=8),allocatable:: res(:,:,:)
      integer(ISZ):: i,iszone=1
      real(kind=8):: ff
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      integer(ISZ):: nxlocalcoarse,nylocalcoarse,nzlocalcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dycoarsesqi,dzcoarsesqi
      real(kind=8):: mgscalecoarse
      integer(ISZ):: ixproc,iyproc,izproc
      integer(ISZ):: localbounds(0:5),localboundsc(0:5)
      integer(ISZ):: lxoffsetall(0:fsdecomp%nxprocs-1)
      integer(ISZ):: rxoffsetall(0:fsdecomp%nxprocs-1)
      integer(ISZ):: lyoffsetall(0:fsdecomp%nyprocs-1)
      integer(ISZ):: ryoffsetall(0:fsdecomp%nyprocs-1)
      integer(ISZ):: lzoffsetall(0:fsdecomp%nzprocs-1)
      integer(ISZ):: rzoffsetall(0:fsdecomp%nzprocs-1)
      integer(ISZ):: lxoffset,rxoffset
      integer(ISZ):: lyoffset,ryoffset
      integer(ISZ):: lzoffset,rzoffset
      type(Decomposition):: coarsedecomp
      integer(ISZ):: allocerror
      integer(ISZ):: delx,dely,delz
      real(kind=8):: sss(2)
      logical(ISZ):: lpe0
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Check that the same coarsening was done for the conductors as is
c     --- being done now. Note that only nzlocal is saved for the conductors so
c     --- only it can be checked.
      if (nxlocal .ne. conductors%levelnx(mglevel) .or.
     &    nylocal .ne. conductors%levelny(mglevel) .or.
     &    nzlocal .ne. conductors%levelnz(mglevel)) then
        call kaboom("vcycle: error: the coarsening level of the conductor is
     &inconsistent with that used by the solver")
        return
      endif

      dxsqi = 1./dx**2
      dysqi = 1./dy**2
      dzsqi = 1./dz**2

      localbounds = globalbounds

#ifdef MPIPARALLEL
      ixproc = fsdecomp%ixproc
      iyproc = fsdecomp%iyproc
      izproc = fsdecomp%izproc
      if (fsdecomp%ix(ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iy(iyproc) > 0)          localbounds(2) = -1
      if (fsdecomp%iy(iyproc)+nylocal < ny) localbounds(3) = -1
      if (fsdecomp%iz(izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(izproc)+nzlocal < nz) localbounds(5) = -1
#endif

c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
c      print *, 'vcycle on dimensions ', nxlocal, nylocal, nzlocal
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DOING ', downpasses,
     &     ' SORPASS3D CALLS ON ', nxlocal, nylocal, nzlocal, ' Warp'
#endif

      if (lprintmgarraysumdiagnostic) then
        lpe0=(fsdecomp%ixproc==0.and.fsdecomp%iyproc==0.and.fsdecomp%izproc==0)
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
        if (lpe0) print*,"V1 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V1 rho",mglevel,sss
      endif

c     --- Do initial SOR passes.
      do i=1,downpasses
        call sorpass3d(mglevel,nxlocal,nylocal,nzlocal,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 nx,ny,nz,globalbounds,localbounds,mgparam,mgform,
     &                 lcndbndy,icndbndy,conductors,lprecalccoeffs,fsdecomp)
      enddo

#ifdef WITHCHOMBO
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DID ', downpasses,
     &     ' SORPASS3D CALLS ON ', nxlocal, nylocal, nzlocal, ' Warp'
#endif

      if (lprintmgarraysumdiagnostic) then
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
        if (lpe0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V2 rho",mglevel,sss
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. ny >= 4 .and. nz >= 4 .and.
     &    mglevel < mgmaxlevels) then

        delx = 1
        dely = 1
        delz = 1
        if (nx > nxlocal) delx = 3
        if (ny > nylocal) dely = 3
        if (nz > nzlocal) delz = 3

        allocate(res(-delx:nxlocal+delx,
     &               -dely:nylocal+dely,
     &               -delz:nzlocal+delz),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nxlocal,nylocal,nzlocal
          call kaboom("vcycle: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residual3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                  dxsqi,dysqi,dzsqi,phi,rho,res,
     &                  mglevel,localbounds,mgparam,mgform,.false.,
     &                  lcndbndy,icndbndy,conductors,lprecalccoeffs)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,res,
     &                               delx,dely,delz,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,nylocal,nzlocal,res,
     &                      delx,dely,delz,-3,-1,fsdecomp)
#endif
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                        res,fsdecomp)
          if (lpe0) print*,"V3 res",mglevel,sss
        endif

        call getnextcoarselevel3d(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                            nxcoarse,nycoarse,nzcoarse,
     &                            dxcoarse,dycoarse,dzcoarse)

        dxcoarsesqi = 1./dxcoarse**2
        dycoarsesqi = 1./dycoarse**2
        dzcoarsesqi = 1./dzcoarse**2
c       --- This option is not supported
c       mgscalecoarse = mgscale*dxcoarse*dycoarse*dzcoarse/(dx*dy*dz)
        mgscalecoarse = 1.

        localboundsc = globalbounds

#ifdef MPIPARALLEL
        coarsedecomp%nxglobal = nxcoarse
        coarsedecomp%nyglobal = nycoarse
        coarsedecomp%nzglobal = nzcoarse
        coarsedecomp%mpi_comm_x = fsdecomp%mpi_comm_x
        coarsedecomp%mpi_comm_y = fsdecomp%mpi_comm_y
        coarsedecomp%mpi_comm_z = fsdecomp%mpi_comm_z
        coarsedecomp%ixproc = fsdecomp%ixproc
        coarsedecomp%iyproc = fsdecomp%iyproc
        coarsedecomp%izproc = fsdecomp%izproc
        coarsedecomp%nxprocs = fsdecomp%nxprocs
        coarsedecomp%nyprocs = fsdecomp%nyprocs
        coarsedecomp%nzprocs = fsdecomp%nzprocs
        allocate(coarsedecomp%ix(0:fsdecomp%nxprocs-1))
        allocate(coarsedecomp%nx(0:fsdecomp%nxprocs-1))
        allocate(coarsedecomp%iy(0:fsdecomp%nyprocs-1))
        allocate(coarsedecomp%ny(0:fsdecomp%nyprocs-1))
        allocate(coarsedecomp%iz(0:fsdecomp%nzprocs-1))
        allocate(coarsedecomp%nz(0:fsdecomp%nzprocs-1))
c       --- Find domains in coarser grid
        call mgdividenz(fsdecomp,coarsedecomp,nx,ny,nz,
     &                  nxcoarse,nycoarse,nzcoarse,mgscale)
c       --- Reset value to corrected one
        nxlocalcoarse = coarsedecomp%nx(ixproc)
        nylocalcoarse = coarsedecomp%ny(iyproc)
        nzlocalcoarse = coarsedecomp%nz(izproc)
c       --- Difference between starts and ends of coarse and fine grids.
c       --- Should only be in the range 0-2.
        lxoffsetall = (nxcoarse*fsdecomp%ix-nx*coarsedecomp%ix)
        rxoffsetall = (nx*(coarsedecomp%ix + coarsedecomp%nx) -
     &                 nxcoarse*(fsdecomp%ix + fsdecomp%nx))
        lyoffsetall = (nycoarse*fsdecomp%iy-ny*coarsedecomp%iy)
        ryoffsetall = (ny*(coarsedecomp%iy + coarsedecomp%ny) -
     &                 nycoarse*(fsdecomp%iy + fsdecomp%ny))
        lzoffsetall = (nzcoarse*fsdecomp%iz-nz*coarsedecomp%iz)
        rzoffsetall = (nz*(coarsedecomp%iz + coarsedecomp%nz) -
     &                 nzcoarse*(fsdecomp%iz + fsdecomp%nz))
c       --- Note that the lzoffsetall and rzoffsetall can only be used in
c       --- MPIPARALLEL sections since they will be unallocated in the
c       --- serial code. So, separate scalars are used in code which is
c       --- used in the serial version.
        lxoffset = lxoffsetall(ixproc)
        rxoffset = rxoffsetall(ixproc)
        lyoffset = lyoffsetall(iyproc)
        ryoffset = ryoffsetall(iyproc)
        lzoffset = lzoffsetall(izproc)
        rzoffset = rzoffsetall(izproc)
        if (coarsedecomp%ix(ixproc) > 0) localboundsc(0) = -1
        if (coarsedecomp%ix(ixproc)+nxlocalcoarse < nxcoarse) localboundsc(1) = -1
        if (coarsedecomp%iy(iyproc) > 0) localboundsc(2) = -1
        if (coarsedecomp%iy(iyproc)+nylocalcoarse < nycoarse) localboundsc(3) = -1
        if (coarsedecomp%iz(izproc) > 0) localboundsc(4) = -1
        if (coarsedecomp%iz(izproc)+nzlocalcoarse < nzcoarse) localboundsc(5) = -1
#else
        nxlocalcoarse = nxcoarse
        nylocalcoarse = nycoarse
        nzlocalcoarse = nzcoarse
        lxoffset = 0
        rxoffset = 0
        lyoffset = 0
        ryoffset = 0
        lzoffset = 0
        rzoffset = 0
#endif

c       --- Alloate new work space
        allocate(phicoarse(-1:nxlocalcoarse+1,
     &                     -1:nylocalcoarse+1,
     &                     -1:nzlocalcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxlocalcoarse,nylocalcoarse,nzlocalcoarse
          call kaboom("vcycle: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxlocalcoarse,0:nylocalcoarse,0:nzlocalcoarse),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxlocalcoarse,nylocalcoarse,nzlocalcoarse
          call kaboom("vcycle: allocation error")
          return
        endif

        if (lprintmgarraysumdiagnostic) then
          rhocoarse = 0. ! needed since not all elements set in restrict3d
        endif
c       rhocoarse = 0. ! not needed since all needed elements set in restrict3d
        phicoarse = 0.

c       --- Restriction
        ff = (dxsqi+dysqi+dzsqi)/(dxcoarsesqi + dycoarsesqi + dzcoarsesqi)
        call restrict3d(nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                  res,delx,dely,delz,
     &                  nxcoarse,nycoarse,nzcoarse,
     &                  nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                  rhocoarse,0,0,0,ff,localbounds,localboundsc,
     &                  lxoffset,lyoffset,lzoffset)
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,nylocalcoarse,nzlocalcoarse,0,0,0,
     &                        rhocoarse,coarsedecomp)
          if (lpe0) print*,"V3 rhocoarse",mglevel,sss
        endif

        if (lprecalccoeffs) conductors%coeffs => conductors%coeffs%coarser

c       --- Continue at the next coarsest level.
        call vcycle(mglevel+iszone,mgscalecoarse,
     &              nxcoarse,nycoarse,nzcoarse,
     &              nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &              dxcoarse,dycoarse,dzcoarse,phicoarse,rhocoarse,
     &              rstar,linbend,bendx,globalbounds,mgparam,mgform,
     &              mgmaxlevels,downpasses,uppasses,
     &              lcndbndy,icndbndy,conductors,lprecalccoeffs,coarsedecomp)

        if (lprecalccoeffs) conductors%coeffs => conductors%coeffs%finer

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,nylocalcoarse,nzlocalcoarse,1,1,1,
     7                        phicoarse,coarsedecomp)
          if (lpe0) print*,"V4 phicoarse",mglevel,sss
          sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
          if (lpe0) print*,"V4 phi",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand3d(nx,ny,nz,nxlocal,nylocal,nzlocal,1,1,1,phi,
     &                nxcoarse,nycoarse,nzcoarse,
     &                nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                phicoarse,localbounds,lxoffset,lyoffset,lzoffset,
     &                conductors,lprecalccoeffs)
        call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                                 localbounds,.false.,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,phi,
     &                               1,1,1,-1,-1,localbounds,fsdecomp)
#endif

        deallocate(phicoarse,rhocoarse)
        deallocate(res)

#ifdef MPIPARALLEL
        deallocate(coarsedecomp%ix)
        deallocate(coarsedecomp%nx)
        deallocate(coarsedecomp%iy)
        deallocate(coarsedecomp%ny)
        deallocate(coarsedecomp%iz)
        deallocate(coarsedecomp%nz)
#endif

      endif

      if (lprintmgarraysumdiagnostic) then
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,phi,
     &                               1,1,1,0,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,-1,fsdecomp)
#endif
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
        if (lpe0) print*,"V5 phi",mglevel,sss
      endif

c     --- Do final SOR passes.
#ifdef WITHCHOMBO
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DOING ', uppasses,
     &     ' SORPASS3D CALLS ON ', nxlocal, nylocal, nzlocal, ' Warp'
#endif

      do i=1,uppasses
        call sorpass3d(mglevel,nxlocal,nylocal,nzlocal,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 nx,ny,nz,globalbounds,localbounds,mgparam,mgform,
     &                 lcndbndy,icndbndy,conductors,lprecalccoeffs,fsdecomp)
      enddo

#ifdef WITHCHOMBO
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DID ', uppasses,
     &     ' SORPASS3D CALLS ON ', nxlocal, nylocal, nzlocal, ' Warp'
#endif

      if (lf3dtimesubs .and. mglevel==0) timevcycle = timevcycle + wtime() - substarttime

      return
      end
c=============================================================================
      subroutine restrict3d(nx,ny,nz,nxlocal,nylocal,nzlocal,res,delx,dely,delz,
     &                      nxcoarse,nycoarse,nzcoarse,
     &                      nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                      rhocoarse,delcx,delcy,delcz,
     &                      ff,localbounds,localboundscoarse,
     &                      lxoffset,lyoffset,lzoffset)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      integer(ISZ):: delx,dely,delz,delcx,delcy,delcz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      integer(ISZ):: nxlocalcoarse,nylocalcoarse,nzlocalcoarse
      real(kind=8):: res(-delx:nxlocal+delx,-dely:nylocal+dely,-delz:nzlocal+delz)
      real(kind=8):: rhocoarse(-delcx:nxlocalcoarse+delcx,
     &                         -delcy:nylocalcoarse+delcy,
     &                         -delcz:nzlocalcoarse+delcz)
      real(kind=8):: ff
      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: localboundscoarse(0:5)
      integer(ISZ):: lxoffset,lyoffset,lzoffset
      
c Restrict to a coarser grid.  The factor ff is needed since the
c residual will be used as the source term in the next coarser grid and
c it implicitly includes the factor 0.5/(dxsqi+dysqi+dzsqi) which is
c too small with the current values of dxsqi etc.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixcoarse,iycoarse,izcoarse
      integer(ISZ):: ixcoarsemin,ixcoarsemax,iycoarsemin,iycoarsemax
      integer(ISZ):: izcoarsemin,izcoarsemax
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ),allocatable:: ixmina(:),ixmaxa(:)
      integer(ISZ),allocatable:: iymina(:),iymaxa(:)
      integer(ISZ),allocatable:: izmina(:),izmaxa(:)
      real(kind=8):: r,w,dxi,dyi,dzi,wx(0:3),wy(0:3),wz(0:3)
      real(kind=8),allocatable:: wxa(:,:),wya(:,:),wza(:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      allocate(ixmina(0:nxlocalcoarse),ixmaxa(0:nxlocalcoarse))
      allocate(iymina(0:nylocalcoarse),iymaxa(0:nylocalcoarse))
      allocate(izmina(0:nzlocalcoarse),izmaxa(0:nzlocalcoarse))
      allocate(wxa(0:3,0:nxlocalcoarse))
      allocate(wya(0:3,0:nylocalcoarse))
      allocate(wza(0:3,0:nzlocalcoarse))
      wxa = 0.
      wya = 0.
      wza = 0.

c     --- Set the loop limits, including edges when appropriate.
      ixcoarsemin = 0
      ixcoarsemax = nxlocalcoarse
      iycoarsemin = 0
      iycoarsemax = nylocalcoarse
      izcoarsemin = 0
      izcoarsemax = nzlocalcoarse
      if (localboundscoarse(0) == 0) ixcoarsemin = 1
      if (localboundscoarse(1) == 0) ixcoarsemax = nxlocalcoarse - 1
      if (localboundscoarse(2) == 0) iycoarsemin = 1
      if (localboundscoarse(3) == 0) iycoarsemax = nylocalcoarse - 1
      if (localboundscoarse(4) == 0) izcoarsemin = 1
      if (localboundscoarse(5) == 0) izcoarsemax = nzlocalcoarse - 1

      dxi = 1.*nxcoarse/nx
      if (ny > 0) then
        dyi = 1.*nycoarse/ny
      else
        dyi = 1.
      endif
      dzi = 1.*nzcoarse/nz
c     w = ff*dxi*dyi*dzi

c     --- Precalculate the loop limits and weights. This saves a surprisingly
c     --- substantial amount of time.
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = ((izcoarse-1)*nz - lzoffset + 4*nzcoarse)/nzcoarse-3
        izmax = ((izcoarse+1)*nz - lzoffset - 1)/nzcoarse
        izmina(izcoarse) = max(izmin,-delz)
        izmaxa(izcoarse) = min(izmax,nzlocal+delz)
        do iz=izmin,izmax
          wza(iz-izmin,izcoarse) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzcoarse)*dzi)
        enddo
      enddo

      do iycoarse=iycoarsemin,iycoarsemax
        iymin = ((iycoarse-1)*ny - lyoffset + 4*nycoarse)/nycoarse-3
        iymax = ((iycoarse+1)*ny - lyoffset - 1)/nycoarse
        iymina(iycoarse) = max(iymin,-dely)
        iymaxa(iycoarse) = min(iymax,nylocal+dely)
        do iy=iymin,iymax
          wya(iy-iymin,iycoarse) = 1. - abs(iycoarse - (iy + 1.*lyoffset/nycoarse)*dyi)
        enddo
      enddo

      do ixcoarse=ixcoarsemin,ixcoarsemax
        ixmin = ((ixcoarse-1)*nx - lxoffset + 4*nxcoarse)/nxcoarse-3
        ixmax = ((ixcoarse+1)*nx - lxoffset - 1)/nxcoarse
        ixmina(ixcoarse) = max(ixmin,-delx)
        ixmaxa(ixcoarse) = min(ixmax,nxlocal+delx)
        do ix=ixmin,ixmax
          wxa(ix-ixmin,ixcoarse) = 1. - abs(ixcoarse - (ix + 1.*lxoffset/nxcoarse)*dxi)
        enddo
      enddo

c     --- Do the loops.

      if (localboundscoarse(4) == 0) rhocoarse(:,:,0) = 0.

!$OMP DO
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = izmina(izcoarse)
        izmax = izmaxa(izcoarse)
        if (izmax < izmin) continue
        wz = wza(:,izcoarse)

        if (localboundscoarse(2) == 0) rhocoarse(:,0,izcoarse) = 0.
        do iycoarse=iycoarsemin,iycoarsemax
          iymin = iymina(iycoarse)
          iymax = iymaxa(iycoarse)
          if (iymax < iymin) continue
          wy = wya(:,iycoarse)

          if (localboundscoarse(0) == 0) rhocoarse(0,iycoarse,izcoarse) = 0.
          do ixcoarse=ixcoarsemin,ixcoarsemax
            ixmin = ixmina(ixcoarse)
            ixmax = ixmaxa(ixcoarse)
            if (ixmax < ixmin) continue
            wx = wxa(:,ixcoarse)

            r = 0.
            w = 0.
            do iz=izmin,izmax
              do iy=iymin,iymax
                do ix=ixmin,ixmax
                  r = r + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)*res(ix,iy,iz)
                  w = w + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)
                enddo
              enddo
            enddo
            if (w > 0.) then
              r = r/w
              rhocoarse(ixcoarse,iycoarse,izcoarse) = r*ff
            else
              rhocoarse(ixcoarse,iycoarse,izcoarse) = 0.
            endif

          enddo
          if (localboundscoarse(1) == 0) rhocoarse(nxlocalcoarse,iycoarse,izcoarse) = 0.
        enddo
        if (localboundscoarse(3) == 0) rhocoarse(:,nylocalcoarse,izcoarse) = 0.
      enddo
!$OMP END DO

      if (localboundscoarse(5) == 0) rhocoarse(:,:,nzlocalcoarse) = 0.

c     --- Set appropriate boundary values
c     --- Not needed since rhocoarse at the boundary is never used
c     --- with Dirichlet.
c     if (localbounds(0) == 0) rhocoarse(0,:,:) = 0.
c     if (localbounds(1) == 0) rhocoarse(nxcoarse,:,:) = 0.
c     if (localbounds(2) == 0) rhocoarse(:,0,:) = 0.
c     if (localbounds(3) == 0) rhocoarse(:,nycoarse,:) = 0.
c     if (localboundscoarse(4) == 0) rhocoarse(:,:,0) = 0.
c     if (localboundscoarse(5) == 0) rhocoarse(:,:,nzlocalcoarse) = 0.

      deallocate(ixmina,ixmaxa)
      deallocate(iymina,iymaxa)
      deallocate(izmina,izmaxa)
      deallocate(wxa)
      deallocate(wya)
      deallocate(wza)

      if (lf3dtimesubs) timerestrict3d = timerestrict3d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine expand3d(nx,ny,nz,nxlocal,nylocal,nzlocal,delx,dely,delz,phi,
     &                    nxcoarse,nycoarse,nzcoarse,
     &                    nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                    phicoarse,bounds,lxoffset,lyoffset,lzoffset,
     &                    conductors,lprecalccoeffs)
      use ConductorTypemodule
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal,delx,dely,delz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      integer(ISZ):: nxlocalcoarse,nylocalcoarse,nzlocalcoarse
      real(kind=8):: phi(-delx:nxlocal+delx,-dely:nylocal+dely,-delz:nzlocal+delz)
      real(kind=8):: phicoarse(-delx:nxlocalcoarse+delx,-dely:nylocalcoarse+dely,-delz:nzlocalcoarse+delz)
      integer(ISZ):: lxoffset,lyoffset,lzoffset
      integer(ISZ):: bounds(0:5)
      logical(ISZ):: lprecalccoeffs
      type(ConductorType):: conductors

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: jx,jy,jz
      real(kind=8):: wx,wy,wz
      integer(ISZ),pointer:: jxa(:),jya(:),jza(:)
      real(kind=8),pointer:: wxa(:),wya(:),wza(:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      allocate(jxa(0:nxlocal),jya(0:nylocal),jza(0:nzlocal))
      allocate(wxa(0:nxlocal),wya(0:nylocal),wza(0:nzlocal))

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nxlocal
      iymin = 0
      iymax = nylocal
      izmin = 0
      izmax = nzlocal
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nxlocal - 1
      if (bounds(2) == 0) iymin = 1
      if (bounds(3) == 0) iymax = nylocal - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nzlocal - 1

      if (ny > 0) then

c       --- Precalculate the indices and weights. This saves a surprisingly
c       --- substantial amount of time.
        do iz=izmin,izmax
          jza(iz) = int((iz*nzcoarse + lzoffset)/nz)
          wza(iz) =  1.*(iz*nzcoarse + lzoffset)/nz - jza(iz)
        enddo
        do iy=iymin,iymax
          jya(iy) = int((iy*nycoarse + lyoffset)/ny)
          wya(iy) =  1.*(iy*nycoarse + lyoffset)/ny - jya(iy)
        enddo
        do ix=ixmin,ixmax
          jxa(ix) = int((ix*nxcoarse + lxoffset)/nx)
          wxa(ix) =  1.*(ix*nxcoarse + lxoffset)/nx - jxa(ix)
        enddo

        if (lprecalccoeffs) then

!$OMP DO
          do iz=izmin,izmax
            jz = jza(iz)
            wz = wza(iz)
            do iy=iymin,iymax
              jy = jya(iy)
              wy = wya(iy)
              do ix=ixmin,ixmax
                jx = jxa(ix)
                wx = wxa(ix)

c               --- Skip points inside of conductors
                if (conductors%coeffs%data(7,ix,iy,iz) == 0.) cycle

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
!$OMP END DO

        else

!$OMP DO
          do iz=izmin,izmax
            jz = jza(iz)
            wz = wza(iz)
            do iy=iymin,iymax
              jy = jya(iy)
              wy = wya(iy)
              do ix=ixmin,ixmax
                jx = jxa(ix)
                wx = wxa(ix)

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
!$OMP END DO

        endif

      else

        do iz=izmin,izmax
          jza(iz) = int((iz*nzcoarse + lzoffset)/nz)
          wza(iz) =  1.*(iz*nzcoarse + lzoffset)/nz - jza(iz)
        enddo
        do ix=ixmin,ixmax
          jxa(ix) = int((ix*nxcoarse + lxoffset)/nx)
          wxa(ix) =  1.*(ix*nxcoarse + lxoffset)/nx - jxa(ix)
        enddo

!$OMP DO
        iy = 0
        jy = 0
        do iz=izmin,izmax
          jz = jza(iz)
          wz = wza(iz)
            do ix=ixmin,ixmax
              jx = jxa(ix)
              wx = wxa(ix)

              phi(ix,iy,iz) = phi(ix,iy,iz) +
     &             (1.-wx)*(1.-wz)*phicoarse(jx  ,jy  ,jz  ) +
     &                 wx *(1.-wz)*phicoarse(jx+1,jy  ,jz  ) +
     &             (1.-wx)*    wz *phicoarse(jx  ,jy  ,jz+1) +
     &                 wx *    wz *phicoarse(jx+1,jy  ,jz+1)

            enddo
        enddo
!$OMP END DO

      endif

      deallocate(jxa,jya,jza)
      deallocate(wxa,wya,wza)

      if (lf3dtimesubs) timeexpand3d = timeexpand3d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine sorpass3d(mglevel,nxlocal,nylocal,nzlocal,
     &                     phi,rho,rstar,
     &                     dxsqi,dysqi,dzsqi,linbend,bendx,
     &                     nx,ny,nz,globalbounds,localbounds,
     &                     mgparam,mgform,lcndbndy,icndbndy,conductors,
     &                     lprecalccoeffs,fsdecomp)
      use Subtimersf3d
      use Constant
      use ConductorTypemodule
      use Decompositionmodule
      use Parallel,Only:my_index
      integer(ISZ):: mglevel
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rstar(-1:nzlocal+1)
      real(kind=8):: bendx(-1:nxlocal+1)
      real(kind=8):: dxsqi,dysqi,dzsqi
      logical(ISZ):: linbend
      integer(ISZ):: nx,ny,nz,globalbounds(0:5),localbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
      logical(ISZ):: lcndbndy,lprecalccoeffs
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

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
c rstar(-1) and rstar(nzlocal+1) are set based on the axial boundary conditions.

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ix,iy,iz
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Put desired potential onto conductors in phi array.
      if (.not. lprecalccoeffs) then
        call cond_potmg(conductors%interior,nxlocal,nylocal,nzlocal,1,1,1,
     &                  phi,mglevel,mgform,.false.)
      endif
      call condbndymgint(conductors%evensubgrid,nxlocal,nylocal,nzlocal,
     &                   1,1,1,phi,localbounds,mglevel,icndbndy)
      call condbndymgint(conductors%oddsubgrid ,nxlocal,nylocal,nzlocal,
     &                   1,1,1,phi,localbounds,mglevel,icndbndy)

c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      parity = + fsdecomp%ix(fsdecomp%ixproc)
     &         + fsdecomp%iy(fsdecomp%iyproc)
     &         + fsdecomp%iz(fsdecomp%izproc)
      s_parity = mod(parity,2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

        call sorhalfpass3d(parity,mglevel,nxlocal,nylocal,nzlocal,
     &                     phi,rho,rstar,dxsqi,dysqi,dzsqi,linbend,bendx,
     &                     localbounds,mgparam,mgform,
     &                     lcndbndy,icndbndy,conductors,lprecalccoeffs)

        call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                                 localbounds,.false.,.false.)

#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,phi,
     &                               1,1,1,-1,0,localbounds,fsdecomp)
        if (parity == s_parity) then
          call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,0,0,fsdecomp)
        else
          call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,0,fsdecomp)
        endif
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- This doesn't seem to be needed.
c     --- Exchange phi in the z guard planes
c     call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,-1,fsdecomp)
#endif

      if (lf3dtimesubs) timesorpass3d = timesorpass3d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine sorhalfpass3d(parity,mglevel,nxlocal,nylocal,nzlocal,
     &                         phi,rho,rstar,
     &                         dxsqi,dysqi,dzsqi,linbend,bendx,localbounds,
     &                         mgparam,mgform,
     &                         lcndbndy,icndbndy,conductors,lprecalccoeffs)
      use Subtimersf3d
      use Constant
      use ConductorTypemodule
      integer(ISZ):: parity,mglevel,nxlocal,nylocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rstar(-1:nzlocal+1)
      real(kind=8):: bendx(-1:nxlocal+1)
      real(kind=8):: dxsqi,dysqi,dzsqi
      logical(ISZ):: linbend
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
      logical(ISZ):: lcndbndy,lprecalccoeffs
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors

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
c rstar(-1) and rstar(nzlocal+1) are set based on the axial boundary conditions.

      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      real(kind=8):: rdel,const,dxsqic,dysqic,dzsqic,spm1,dx
      integer(ISZ):: nxy,iimx,iipx,iimy,iipy,iimz,iipz
      integer(ISZ):: ii,ix,iy,iz,ix1,ic,i1,i2
      integer(ISZ):: nphi,nrho,nrstar
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      nphi = (3+nxlocal)*(3+nylocal)*(3+nzlocal)
      nrho = (1+nxlocal)*(1+nylocal)*(1+nzlocal)
      nrstar = 3+nzlocal

c     --- Set temporary variables (these are used to increase performance)
      dx = 1./sqrt(dxsqi)
      rdel = dzsqi/(dxsqi + dysqi + dzsqi)
      const = mgparam*0.5/(dxsqi + dysqi + dzsqi)
      dxsqic = dxsqi*const
      dysqic = dysqi*const
      dzsqic = dzsqi*const
      spm1 = 1. - mgparam

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nxlocal
      iymin = 0
      iymax = nylocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) == 0) ixmin = 1
      if (localbounds(1) == 0) ixmax = nxlocal - 1
      if (localbounds(2) == 0) iymin = 1
      if (localbounds(3) == 0) iymax = nylocal - 1
      if (localbounds(4) == 0) izmin = 1
      if (localbounds(5) == 0) izmax = nzlocal - 1

c     --- Save values just outside conductor surfaces. Only save phi at the
c     --- subgrid points which are to be used at the current level of
c     --- grid refinement.
      if (lcndbndy .and. .not. lprecalccoeffs) then
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

c     --- guard planes in z are already set

      if (lprecalccoeffs) then

!$OMP DO
        do iz=izmin,izmax
          do iy=iymin,iymax
            ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
            do ix=ix1,ixmax,2
              if (conductors%coeffs%data(7,ix,iy,iz) == 0.) then
                phi(ix,iy,iz) = conductors%coeffs%data(6,ix,iy,iz)
              else
                phi(ix,iy,iz) = (rho(ix,iy,iz)
     &            + phi(ix-1,iy  ,iz  )*conductors%coeffs%data(0,ix,iy,iz)
     &            + phi(ix+1,iy  ,iz  )*conductors%coeffs%data(1,ix,iy,iz)
     &            + phi(ix  ,iy-1,iz  )*conductors%coeffs%data(2,ix,iy,iz)
     &            + phi(ix  ,iy+1,iz  )*conductors%coeffs%data(3,ix,iy,iz)
     &            + phi(ix  ,iy  ,iz-1)*conductors%coeffs%data(4,ix,iy,iz)
     &            + phi(ix  ,iy  ,iz+1)*conductors%coeffs%data(5,ix,iy,iz)
     &            + conductors%coeffs%data(6,ix,iy,iz))
     &            *conductors%coeffs%data(7,ix,iy,iz)
     &            + spm1*phi(ix,iy,iz)

              endif
            enddo
          enddo
        enddo
!$OMP END DO

      else

c       --- Loop over the array. Boundary points are calculated too.
        call mgsor_loop(ixmin,ixmax,iymin,iymax,izmin,izmax,parity,
     &                  phi,rho,nxlocal,nylocal,nzlocal,dxsqic,dysqic,dzsqic,
     &                  spm1,linbend,
     &                  rstar,dx,dxsqi,dysqi,dzsqi,rdel,const,bendx)

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
          if (parity == 0) then
           call condbndymg(conductors%evensubgrid,nxlocal,nylocal,nzlocal,1,1,1,
     &                     phi,rho,
     &                     dxsqic,dysqic,dzsqic,spm1,mgparam,localbounds,
     &                     mglevel,mgform,icndbndy)
          endif
          if (parity == 1) then
           call condbndymg(conductors%oddsubgrid,nxlocal,nylocal,nzlocal,1,1,1,
     &                     phi,rho,
     &                     dxsqic,dysqic,dzsqic,spm1,mgparam,localbounds,
     &                     mglevel,mgform,icndbndy)
          endif

        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmg(conductors%interior,nxlocal,nylocal,nzlocal,1,1,1,
     &                  phi,mglevel,mgform,.false.)

      endif

c     -- These again...
      call condbndymgint(conductors%evensubgrid,nxlocal,nylocal,nzlocal,
     &                   1,1,1,phi,localbounds,mglevel,icndbndy)
      call condbndymgint(conductors%oddsubgrid,nxlocal,nylocal,nzlocal,
     &                   1,1,1,phi,localbounds,mglevel,icndbndy)

      if (lf3dtimesubs) timesorhalfpass3d = timesorhalfpass3d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine mgsor_loop(ixmin,ixmax,iymin,iymax,izmin,izmax,parity,
     &                      phi,rho,nxlocal,nylocal,nzlocal,
     &                      dxsqic,dysqic,dzsqic,spm1,linbends,
     &                      rstar,dx,dxsqi,dysqi,dzsqi,rdel,const,bendx)
      use Subtimersf3d
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: parity,nxlocal,nylocal,nzlocal,ii1
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rstar(-1:nzlocal+1),bendx(-1:nxlocal+1)
      real(kind=8):: dxsqic,dysqic,dzsqic,spm1,dx,dxsqi,dysqi,dzsqi,rdel,const
      logical(ISZ):: linbends

c This routine provides the loops over the phi array which solves the 
c iterative equation.
c The case with and without bends are seperated since
c the bends case has extra terms.

      integer(ISZ):: ix,iy,iz,ix1
      real(kind=8):: x,rs,r,xfact,del2bndc,dh,dxi
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (.not. linbends) then 
!$OMP DO
        do iz=izmin,izmax
          do iy=iymin,iymax
            ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
            do ix=ix1,ixmax,2
              phi(ix,iy,iz) = rho(ix,iy,iz) +
     &           (phi(ix-1,iy  ,iz  ) + phi(ix+1,iy  ,iz  ))*dxsqic +
     &           (phi(ix  ,iy-1,iz  ) + phi(ix  ,iy+1,iz  ))*dysqic +
     &           (phi(ix  ,iy  ,iz-1) + phi(ix  ,iy  ,iz+1))*dzsqic +
     &           spm1*phi(ix,iy,iz)
            enddo
          enddo
        enddo
!$OMP END DO
      else
c       --- bends, loop over z slices.  Comment: timing tests indicate that 
c       --- use of a 1d array is faster than a 3d array.   
        dxi = 1./dx
!$OMP DO
        do iz=izmin,izmax
          rs  = rstar(iz+1)
          if (maxval(abs(rstar(iz-1:iz+1))) > LARGEPOS*1.e-6) then
c           --- current z-slice not in bend or bordering a bend, do not 
c           --- include curvature terms from Poisson's equation    
            do iy=iymin,iymax
              ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
              do ix=ix1,ixmax,2
                phi(ix,iy,iz) = rho(ix,iy,iz) +
     &             (phi(ix-1,iy  ,iz  ) + phi(ix+1,iy  ,iz  ))*dxsqic +
     &             (phi(ix  ,iy-1,iz  ) + phi(ix  ,iy+1,iz  ))*dysqic +
     &             (phi(ix  ,iy  ,iz-1) + phi(ix  ,iy  ,iz+1))*dzsqic +
     &             spm1*phi(ix,iy,iz)
              enddo  
            enddo  
          else  
c           --- current z-slice in a bend and near border of a bend, include 
c           --- curvature terms from Poisson's equation.  Note: if bordering 
c           --- a bend, the dh/dz "jump term" should be included.  The 
c           --- current method of treating this jump term may be inaccurate. 
c           --- r is set to 1 over (rs+x) to reduce the number of divides.
            dh = 1./rstar(iz+1) - 1./rstar(iz-1)
            if (abs(dh) > SMALLPOS .and. abs(rs) < LARGEPOS*1.e-6) then
              do iy=iymin,iymax
                ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
                do ix=ix1,ixmax,2
                  x = bendx(ix)
                  r = 1./(rs + x)
                  xfact = (x*r)*((x*r)-2.) 
                  del2bndc = const/( 1.+xfact*rdel )
                  phi(ix,iy,iz) = rho(ix,iy,iz) + del2bndc*( 
     &             (phi(ix-1,iy  ,iz  ) + phi(ix+1,iy  ,iz  ))*dxsqi +
     &             (phi(ix  ,iy-1,iz  ) + phi(ix  ,iy+1,iz  ))*dysqi +
     &             (phi(ix  ,iy  ,iz-1) + phi(ix  ,iy  ,iz+1))*dzsqi +
     &              .5*(phi(ix+1,iy,iz)-phi(ix-1,iy,iz))*r*dxi + 
     &              xfact*(phi(ix,iy,iz+1) + phi(ix,iy,iz-1))*dzsqi - 
     &              .25*(rs*r)**3*x*(phi(ix,iy,iz+1) - phi(ix,iy,iz-1))*dh*dzsqi) +
     &              spm1*phi(ix,iy,iz)
                enddo 
              enddo 

c           --- current z-slice is near border of a bend only
            elseif (abs(dh) > SMALLPOS) then
              do iy=iymin,iymax
                ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
                do ix=ix1,ixmax,2
                  x = bendx(ix)
                  phi(ix,iy,iz) = rho(ix,iy,iz) + const*(
     &             (phi(ix-1,iy  ,iz  ) + phi(ix+1,iy  ,iz  ))*dxsqi +
     &             (phi(ix  ,iy-1,iz  ) + phi(ix  ,iy+1,iz  ))*dysqi +
     &             (phi(ix  ,iy  ,iz-1) + phi(ix  ,iy  ,iz+1))*dzsqi
     &              - .25*x*(phi(ix,iy,iz+1) - phi(ix,iy,iz-1))*dh*dzsqi)
     &              + spm1*phi(ix,iy,iz)
                enddo
              enddo

c           --- current z-slice is in a bend only
            elseif (abs(rs) < LARGEPOS*1.e-6) then
              do iy=iymin,iymax
                ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
                do ix=ix1,ixmax,2
                  x = bendx(ix)
                  r = 1./(rs + x)
                  xfact = (x*r)*((x*r)-2.)
                  del2bndc = const/( 1.+xfact*rdel )
                  phi(ix,iy,iz) = rho(ix,iy,iz) + del2bndc*(
     &             (phi(ix-1,iy  ,iz  ) + phi(ix+1,iy  ,iz  ))*dxsqi +
     &             (phi(ix  ,iy-1,iz  ) + phi(ix  ,iy+1,iz  ))*dysqi +
     &             (phi(ix  ,iy  ,iz-1) + phi(ix  ,iy  ,iz+1))*dzsqi +
     &              .5*(phi(ix+1,iy,iz)-phi(ix-1,iy,iz))*r*dxi + 
     &              xfact*(phi(ix,iy,iz+1) + phi(ix,iy,iz-1))*dzsqi)
     &              + spm1*phi(ix,iy,iz)
                enddo
              enddo
            endif
          endif 
        enddo 
!$OMP END DO
      endif 

      if (lf3dtimesubs) timemgsor_loop = timemgsor_loop +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine cond_potmg(interior,nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                      phi,mglevel,mgform,mgform2init)
      use Subtimersf3d
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nxlocal,nylocal,nzlocal,delx,dely,delz,mglevel,mgform
      real(kind=8):: phi(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      logical(ISZ):: mgform2init

c Set conductor points to the desired potential.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- When at the finest level and not calculating the residual, set
c     --- phi to the voltage of the conductor, otherwise, set it to zero.
      if (mglevel == 0 .and. (mgform == 1 .or. mgform2init)) then
!$OMP DO
        do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
          ix = interior%indx(0,ic)
          iy = interior%indx(1,ic)
          iz = interior%indx(2,ic)
          phi(ix,iy,iz) = interior%volt(ic)
        enddo
!$OMP END DO
      else
!$OMP DO
        do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
          ix = interior%indx(0,ic)
          iy = interior%indx(1,ic)
          iz = interior%indx(2,ic)
          phi(ix,iy,iz) = 0.
        enddo
!$OMP END DO
      endif

      if (lf3dtimesubs) timecond_potmg = timecond_potmg +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine cond_potmgres(interior,nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                         res,mglevel,mgform,mgform2init)
      use Subtimersf3d
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nxlocal,nylocal,nzlocal,mglevel,mgform,delx,dely,delz
      real(kind=8):: res(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      logical(ISZ):: mgform2init

c Set conductor points to the desired potential.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

!$OMP DO
      do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        res(ix,iy,iz) = 0.
      enddo
!$OMP END DO

      if (lf3dtimesubs) timecond_potmgres = timecond_potmgres +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine cond_zerorhointerior(interior,nxlocal,nylocal,nzlocal,rho)
      use Subtimersf3d
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
 
c Set rho to zero inside conductor points.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

!$OMP DO
      do ic = interior%istart(0),interior%istart(1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        rho(ix,iy,iz) = 0.
      enddo
!$OMP END DO

      if (lf3dtimesubs) timecond_zerorhointerior = timecond_zerorhointerior +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      real(kind=8) function cond_sumrhointerior(interior,
     &                                          nxlocal,nylocal,nzlocal,rho,
     &                                  ixmin,ixmax,iymin,iymax,izmin,izmax)
      use Subtimersf3d
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
 
c Set rho to zero inside conductor points.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: rhosum
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      rhosum = 0.

!$OMP DO
      do ic = interior%istart(0),interior%istart(1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        if (ixmin <= ix .and. ix <= ixmax .and.
     &      iymin <= iy .and. iy <= iymax .and.
     &      izmin <= iz .and. iz <= izmax) then
          rhosum = rhosum + rho(ix,iy,iz)
        endif
      enddo
!$OMP END DO

      cond_sumrhointerior = rhosum
      if (lf3dtimesubs) timecond_sumrhointerior = timecond_sumrhointerior +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine subcond_sumrhointerior(rhosum,interior,
     &                                  nxlocal,nylocal,nzlocal,rho,
     &                                  ixmin,ixmax,iymin,iymax,izmin,izmax)
      use Subtimersf3d
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
 
c Set rho to zero inside conductor points.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: rhosum
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      rhosum = 0.

!$OMP DO
      do ic = interior%istart(0),interior%istart(1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        if (ixmin <= ix .and. ix <= ixmax .and.
     &      iymin <= iy .and. iy <= iymax .and.
     &      izmin <= iz .and. iz <= izmax) then
          rhosum = rhosum + rho(ix,iy,iz)
        endif
      enddo
!$OMP END DO

      if (lf3dtimesubs) timesubcond_sumrhointerior = timesubcond_sumrhointerior +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine handlesubgrid2(im,ip,pm,pp,dsqicos,ic,nc,dels,volt,voltfac,
     &                          denom,ppp)
      integer(ISZ):: im,ip
      real(kind=8):: pm,pp,dsqicos
      integer(ISZ):: ic,nc
      real(kind=8):: dels(0:5,nc),volt(0:5,nc),voltfac,denom
      real(kind=8):: ppp

c Handle the special coding needed for subgrid points, for icndbndy=2.
c
c For Dirichlet points, this method uses the first order, non-symmetric
c finite difference form of the discrete Poisson equation. In 1-d, with grid
c cell size D, for a conductor a distance d below point i, Poisson's equation
c is written as
c
c ((phiC - phi_i)/d - (phi_i - phi_i+1)/D)/((d+D)/2)
c
c where phiC is the voltage on the conductor. Each direction is treated
c independently. This is then solved for phi_i.
c
c For Neumann points, the equations can be derived in two ways. For the first,
c with the same situation as above, the first finite difference derivative is
c centered about d, and phiC becomes the potential at a virtual point 2*d below
c point i. The equation becomes
c
c ((phiC - phi_i)/(2*d) - (phi_i - phi_i+1)/D)/(d+D/2)
c
c By the definition of the Neumann boundary, the first finite difference term
c is zero, so the equation reduces to
c
c (- (phi_i - phi_i+1)/D)/(d+D/2)
c
c This expression can also be obtained by fitting a parabola through the
c three points, i-d, i, i+1, where phi_i and phi_i+1 are known, and the
c derivative at i-d is zero by definition. The value at phi_i-1 is obtained
c and plugged into the standard for of the finit difference equation. Doing
c the algebra produces the same result.
c
c Note that special cases are needed when there are Dirichlet and/or Neumann
c subgrid points on both sides of phi_i.

      real(kind=8):: adel

      if (0. < dels(im,ic) .and. dels(im,ic) < +1. .and.
     &    0. < dels(ip,ic) .and. dels(ip,ic) < +1.) then
c       --- Both terms are Dirichlet subgrid
        pm = voltfac*volt(im,ic)*2./(dels(im,ic)*(dels(ip,ic)+dels(im,ic)))
        pp = voltfac*volt(ip,ic)*2./(dels(ip,ic)*(dels(ip,ic)+dels(im,ic)))
        denom = denom +
     &        2.*(1.-dels(ip,ic)*dels(im,ic))/(dels(ip,ic)*dels(im,ic))*dsqicos
        ppp = min(ppp,dels(im,ic))
        ppp = min(ppp,dels(ip,ic))
      else if (0. >= dels(im,ic) .and. dels(im,ic) > -1. .and.
     &         0. >= dels(ip,ic) .and. dels(ip,ic) > -1.) then
c       --- Both terms are Neumann subgrid
        pm = 0.
        pp = 0.
        denom = denom - 2.*dsqicos
        if (abs(dels(im,ic)) > 0.) then
          ppp = min(ppp,abs(dels(im,ic)))
        else
          ppp = min(ppp,1.-1.e-9)
        endif
        if (abs(dels(ip,ic)) > 0.) then
          ppp = min(ppp,abs(dels(ip,ic)))
        else
          ppp = min(ppp,1.-1.e-9)
        endif
      else if (0. <  dels(im,ic) .and. dels(im,ic) < +1. .and.
     &         0. >= dels(ip,ic) .and. dels(ip,ic) > -1.) then
c       --- Minus term is Dirichlet and plus term is Neumann subgrid
        adel = abs(dels(ip,ic))
        pm = voltfac*volt(im,ic)/(dels(im,ic)*(adel+0.5*dels(im,ic)))
        pp = 0.
        denom = denom + (1.-2.*dels(im,ic)*adel-dels(im,ic))/
     &                  (dels(im,ic)*(adel+0.5*dels(im,ic)))*dsqicos
        ppp = min(ppp,dels(im,ic))
        if (adel > 0.) then
          ppp = min(ppp,adel)
        else
          ppp = min(ppp,1.-1.e-9)
        endif
      else if (0. >= dels(im,ic) .and. dels(im,ic) > -1. .and.
     &         0. <  dels(ip,ic) .and. dels(ip,ic) < +1.) then
c       --- Minus term is Neumann and plus term is Dirichlet subgrid
        adel = abs(dels(im,ic))
        pm = 0.
        pp = voltfac*volt(ip,ic)/(dels(ip,ic)*(adel+0.5*dels(ip,ic)))
        denom = denom + (1.-2.*dels(ip,ic)*adel-dels(ip,ic))/
     &                  (dels(ip,ic)*(adel+0.5*dels(ip,ic)))*dsqicos
        if (adel > 0.) then
          ppp = min(ppp,adel)
        else
          ppp = min(ppp,1.-1.e-9)
        endif
        ppp = min(ppp,dels(ip,ic))
      else
        if (0. < dels(im,ic) .and. dels(im,ic) < +1.) then
c         --- Minus term is Dirichlet subgrid and plus term is in bulk
          pm = voltfac*volt(im,ic)*2./(dels(im,ic)*(1.+dels(im,ic)))
          pp = pp*(2./(1.+dels(im,ic)))
          denom = denom + 2.*(1.-dels(im,ic))/dels(im,ic)*dsqicos
          ppp = min(ppp,dels(im,ic))
        else if (0. >= dels(im,ic) .and. dels(im,ic) > -1.) then
c         --- Minus term is Neumann subgrid and plus term is in bulk
          adel = abs(dels(im,ic))
          pm = 0.
          pp = pp/(0.5+adel)
          denom = denom - 2.*adel/(adel + 0.5)*dsqicos
          if (adel > 0.) then
            ppp = min(ppp,adel)
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
        if (0. < dels(ip,ic) .and. dels(ip,ic) < +1.) then
c         --- Minus term is in bulk and plus term is Dirichlet subgrid
          pp = voltfac*volt(ip,ic)*2./(dels(ip,ic)*(1.+dels(ip,ic)))
          pm = pm*(2./(1.+dels(ip,ic)))
          denom = denom + 2.*(1.-dels(ip,ic))/dels(ip,ic)*dsqicos
          ppp = min(ppp,dels(ip,ic))
        else if (0. >= dels(ip,ic) .and. dels(ip,ic) > -1.) then
c         --- Minus term is in bulk and plus term is Neumann subgrid
          adel = abs(dels(ip,ic))
          pp = 0.
          pm = pm/(0.5+adel)
          denom = denom - 2.*adel/(adel + 0.5)*dsqicos
          if (adel > 0.) then
            ppp = min(ppp,adel)
          else
            ppp = min(ppp,1.-1.e-9)
          endif
        endif
      endif

      return
      end
c=============================================================================
      subroutine condbndymg(subgrid,nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                      phi,rho,dxsqic,dysqic,dzsqic,spm1,
     &                      mgparam,localbounds,mglevel,mgform,icndbndy)
      use Subtimersf3d
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nylocal,nzlocal,delx,dely,delz,mglevel
      real(kind=8):: phi(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: dxsqic,dysqic,dzsqic,spm1,mgparam
      integer(ISZ):: localbounds(0:5),mgform,icndbndy

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c NOTE that dxsqicos and dysqicos are dxsqic and dysqic over mgparam.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: dxsqicos,dysqicos,dzsqicos,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: voltfac
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz
      real(kind=8):: ppp
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      dxsqicos = dxsqic/mgparam
      dysqicos = dysqic/mgparam
      dzsqicos = dzsqic/mgparam
      dels => subgrid%dels
      volt => subgrid%volt

      if (icndbndy == 1) then
c     --- Linear interpolation

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

        if (ix == 0  .and. localbounds(0) < 1) cycle
        if (ix == nxlocal .and. localbounds(1) < 1) cycle
        if (iy == 0  .and. localbounds(2) < 1) cycle
        if (iy == nylocal .and. localbounds(3) < 1) cycle
        if (iz == 0  .and. localbounds(4) < 1) cycle
        if (iz == nzlocal .and. localbounds(5) < 1) cycle

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1
        izp1 = iz + 1
        izm1 = iz - 1

        if (nylocal == 0) then
          iym1 = iy
          iyp1 = iy
        endif

c       --- Set temporaries with initial values.
        pxm = phi(ixm1,iy  ,iz  )
        pxp = phi(ixp1,iy  ,iz  )
        pym = phi(ix  ,iym1,iz  )
        pyp = phi(ix  ,iyp1,iz  )
        pzm = phi(ix  ,iy  ,izm1)
        pzp = phi(ix  ,iy  ,izp1)
        denom = 1.
        dosubgrid = .false.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0 .and. mgform == 1) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- the point lower in x is inside the conductor
        if (0 < dels(0,ic) .and. dels(0,ic) < 1.) then
          pxm = voltfac*volt(0,ic)/dels(0,ic)
          denom = denom + (1.-dels(0,ic))/dels(0,ic)*dxsqicos
          dosubgrid = .true.
        endif
c       --- the point higher in x is inside the conductor
        if (0 < dels(1,ic) .and. dels(1,ic) < 1.) then
          pxp = voltfac*volt(1,ic)/dels(1,ic)
          denom = denom + (1.-dels(1,ic))/dels(1,ic)*dxsqicos
          dosubgrid = .true.
        endif
c       --- the point lower in y is inside the conductor
        if (0 < dels(2,ic) .and. dels(2,ic) < 1. .and. nylocal > 0) then
          pym = voltfac*volt(2,ic)/dels(2,ic)
          denom = denom + (1.-dels(2,ic))/dels(2,ic)*dysqicos
          dosubgrid = .true.
        endif
c       --- the point higher in y is inside the conductor
        if (0 < dels(3,ic) .and. dels(3,ic) < 1. .and. nylocal > 0) then
          pyp = voltfac*volt(3,ic)/dels(3,ic)
          denom = denom + (1.-dels(3,ic))/dels(3,ic)*dysqicos
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (0 < dels(4,ic) .and. dels(4,ic) < 1.) then
          pzm = voltfac*volt(4,ic)/dels(4,ic)
          denom = denom + (1.-dels(4,ic))/dels(4,ic)*dzsqicos
          dosubgrid = .true.
        endif
c       --- the point higher in z is inside the conductor
        if (0 < dels(5,ic) .and. dels(5,ic) < 1.) then
          pzp = voltfac*volt(5,ic)/dels(5,ic)
          denom = denom + (1.-dels(5,ic))/dels(5,ic)*dzsqicos
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          phi(ix,iy,iz) = (rho(ix,iy,iz) +
     &      (pxm+pxp)*dxsqic + (pym+pyp)*dysqic + (pzm+pzp)*dzsqic)/denom +
     &      spm1*subgrid%prevphi(ic)
        endif
      enddo
!$OMP END DO

      elseif (icndbndy == 2) then
c     --- Quadratic interpolation

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

        if (ix == 0  .and. localbounds(0) < 1) cycle
        if (ix == nxlocal .and. localbounds(1) < 1) cycle
        if (iy == 0  .and. localbounds(2) < 1) cycle
        if (iy == nylocal .and. localbounds(3) < 1) cycle
        if (iz == 0  .and. localbounds(4) < 1) cycle
        if (iz == nzlocal .and. localbounds(5) < 1) cycle

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1
        izp1 = iz + 1
        izm1 = iz - 1

        if (nylocal == 0) then
          iym1 = iy
          iyp1 = iy
        endif

c       --- Set temporaries with initial values.
        pxm = phi(ixm1,iy   ,iz   )
        pxp = phi(ixp1,iy   ,iz   )
        pym = phi(ix   ,iym1,iz   )
        pyp = phi(ix   ,iyp1,iz   )
        pzm = phi(ix   ,iy   ,izm1)
        pzp = phi(ix   ,iy   ,izp1)
        denom = 1.
        dosubgrid = .false.
        ppp = 1.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0 .and. mgform == 1) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

        call handlesubgrid2(0,1,pxm,pxp,dxsqicos,ic,subgrid%nmax,dels,volt,voltfac,denom,ppp)
        if (nylocal > 0) call handlesubgrid2(2,3,pym,pyp,dysqicos,ic,subgrid%nmax,dels,volt,voltfac,denom,ppp)
        call handlesubgrid2(4,5,pzm,pzp,dzsqicos,ic,subgrid%nmax,dels,volt,voltfac,denom,ppp)

c       --- calculate the new phi based on the boundary conditions
        if (ppp < 1.) then
          phi(ix,iy,iz) = (rho(ix,iy,iz) +
     &      (pxm+pxp)*dxsqic + (pym+pyp)*dysqic + (pzm+pzp)*dzsqic)/denom +
     &      spm1*subgrid%prevphi(ic)
        endif

      enddo
!$OMP END DO

      endif

      if (lf3dtimesubs) timecondbndymg = timecondbndymg +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine condbndymgint(subgrid,nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                         phi,localbounds,mglevel,icndbndy)
      use Subtimersf3d
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nylocal,nzlocal,delx,dely,delz,mglevel
      real(kind=8):: phi(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      integer(ISZ):: localbounds(0:5),icndbndy

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c NOTE that dxsqicos and dysqicos are dxsqic and dysqic over mgparam.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz
      real(kind=8),pointer:: dels(:,:),volt(:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      dels => subgrid%dels
      volt => subgrid%volt

      if (icndbndy == 1) then
c     --- Linear interpolation
        return

      elseif (icndbndy == 2) then
c     --- Quadratic interpolation

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

c       --- Old code, kept in case there is a problem with the changes.
c       if (.not. ((ixm1 == 0 .and. localbounds(0) == 0) .or.
c    &             (ixm1 == -1 .and. localbounds(0) == 1) .or.
c    &             (ixm1 == -1 .and. localbounds(0) == 2) .or.
c    &             (ixp1 == nx .and. localbounds(1) == 0) .or.
c    &             (ixp1 == nx+1 .and. localbounds(1) == 1) .or.
c    &             (ixp1 == nx+1 .and. localbounds(1) == 2) .or.
c    &             (ixm1 == -1 .or. ixp1 == nx+1))) then
c         if (0. >= dels(0,ic) .and. dels(0,ic) > -1.) phi(ixm1,iy,iz)=phi(ix,iy,iz)
c         if (0. >= dels(1,ic) .and. dels(1,ic) > -1.) phi(ixp1,iy,iz)=phi(ix,iy,iz)
c       endif

c       if (.not. ((iym1 == 0 .and. localbounds(2) == 0) .or.
c    &             (iym1 == -1 .and. localbounds(2) == 1) .or.
c    &             (iym1 == -1 .and. localbounds(2) == 2) .or.
c    &             (iyp1 == ny .and. localbounds(3) == 0) .or.
c    &             (iyp1 == ny+1 .and. localbounds(3) == 1) .or.
c    &             (iyp1 == ny+1 .and. localbounds(3) == 2) .or.
c    &             (iym1 == -1 .or. iyp1 == ny+1))) then
c         if (0. >= dels(2,ic) .and. dels(2,ic) > -1.) phi(ix,iym1,iz)=phi(ix,iy,iz)
c         if (0. >= dels(3,ic) .and. dels(3,ic) > -1.) phi(ix,iyp1,iz)=phi(ix,iy,iz)
c       endif

        if (.not. ((ixm1 == 0 .and. localbounds(0) == 0) .or.
     &             (ixp1 == nxlocal .and. localbounds(1) == 0) .or.
     &             (ix == 0  .and. localbounds(0) < 1) .or.
     &             (ix == nxlocal .and. localbounds(1) < 1))) then
          if (0. >= dels(0,ic) .and. dels(0,ic) > -1.) phi(ixm1,iy,iz)=phi(ix,iy,iz)
          if (0. >= dels(1,ic) .and. dels(1,ic) > -1.) phi(ixp1,iy,iz)=phi(ix,iy,iz)
        endif

        if (.not. ((iym1 == 0 .and. localbounds(2) == 0) .or.
     &             (iyp1 == nylocal .and. localbounds(3) == 0) .or.
     &             (iy == 0  .and. localbounds(2) < 1) .or.
     &             (iy == nylocal .and. localbounds(3) < 1))) then
          if (0. >= dels(2,ic) .and. dels(2,ic) > -1.) phi(ix,iym1,iz)=phi(ix,iy,iz)
          if (0. >= dels(3,ic) .and. dels(3,ic) > -1.) phi(ix,iyp1,iz)=phi(ix,iy,iz)
        endif

        if (.not. ((izm1 == 0 .and. localbounds(4) == 0) .or.
     &             (izp1 == nzlocal .and. localbounds(5) == 0) .or.
     &             (iz == 0  .and. localbounds(4) < 1) .or.
     &             (iz == nzlocal .and. localbounds(5) < 1))) then
          if (0. >= dels(4,ic) .and. dels(4,ic) > -1.) phi(ix,iy,izm1)=phi(ix,iy,iz)
          if (0. >= dels(5,ic) .and. dels(5,ic) > -1.) phi(ix,iy,izp1)=phi(ix,iy,iz)
        endif

      enddo
!$OMP END DO

      endif

      if (lf3dtimesubs) timecondbndymgint = timecondbndymgint +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine condbndyres(subgrid,nxlocal,nylocal,nzlocal,
     &                       delx,dely,delz,
     &                       phi,rho,res,dxsqi,dysqi,dzsqi,
     &                       mgparam,localbounds,mglevel,
     &                       mgform,mgform2init,icndbndy)
      use Subtimersf3d
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nylocal,nzlocal,delx,dely,delz
      integer(ISZ):: mglevel
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: res(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      real(kind=8):: dxsqi,dysqi,dzsqi,mgparam
      integer(ISZ):: localbounds(0:5),mgform,icndbndy
      logical(ISZ):: mgform2init

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c
c mgform is used to specify what form of operator is being used for
c multigrid. This only effects the routine at the finest level. When 1, normal
c form is being used - the residual is calculated directly from phi at level 1.
c Hence the actual values of the voltages are used. When 2, residual
c correction form is used and so the residual is being calculated from the
c error. Use zero for the voltages (zero error).
c When the normal form is used, the result is scaled by the minimum of the
c deltas. This is done since the the correct term can get erroneously large
c as delta approaches zero which hinder convergence. With residual correction
c form, the opposite is true, when delta nears zero, the large residual is
c needed to allow rapid convergence.
c The logical mgform2init is true on the first call to residual when the
c residual correction form is being used. In that case, the actual voltages
c need to be used since the residual is operating on phi (and not the error).
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: const,dxsqic,dysqic,dzsqic,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: voltfac
      real(kind=8):: dxsqics,dysqics,dzsqics,ppp
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz
      real(kind=8),pointer:: dels(:,:),volt(:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      const = 0.5/(dxsqi+dysqi+dzsqi)
      dxsqic = dxsqi*const
      dysqic = dysqi*const
      dzsqic = dzsqi*const
      dxsqics = mgparam*dxsqi*const
      dysqics = mgparam*dysqi*const
      dzsqics = mgparam*dzsqi*const
      dels => subgrid%dels
      volt => subgrid%volt

      if (icndbndy == 1) then
c     --- Linear interpolation

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

        if (ix == 0  .and. localbounds(0) == 0) cycle
        if (ix == nxlocal .and. localbounds(1) == 0) cycle
        if (iy == 0  .and. localbounds(2) == 0) cycle
        if (iy == nylocal .and. localbounds(3) == 0) cycle
        if (iz == 0  .and. localbounds(4) == 0) cycle
        if (iz == nzlocal .and. localbounds(5) == 0) cycle

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1
        izp1 = iz + 1
        izm1 = iz - 1

        if (nylocal == 0) then
          iym1 = iy
          iyp1 = iy
        endif

c       --- set temporaries with initial values
        pxm = phi(ixm1,iy  ,iz  )
        pxp = phi(ixp1,iy  ,iz  )
        pym = phi(ix  ,iym1,iz  )
        pyp = phi(ix  ,iyp1,iz  )
        pzm = phi(ix  ,iy  ,izm1)
        pzp = phi(ix  ,iy  ,izp1)
        denom = 1.
        ppp = 1.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0 .and. (mgform == 1 .or. mgform2init)) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- the point lower in x is inside the conductor
        if (0 < dels(0,ic) .and. dels(0,ic) < 1.) then
          pxm = voltfac*volt(0,ic)/dels(0,ic)
          denom = denom + (1.-dels(0,ic))/dels(0,ic)*dxsqic
          ppp = min(ppp,dels(0,ic))
        endif
c       --- the point higher in x is inside the conductor
        if (0 < dels(1,ic) .and. dels(1,ic) < 1.) then
          pxp = voltfac*volt(1,ic)/dels(1,ic)
          denom = denom + (1.-dels(1,ic))/dels(1,ic)*dxsqic
          ppp = min(ppp,dels(1,ic))
        endif
c       --- the point lower in y is inside the conductor
        if (0 < dels(2,ic) .and. dels(2,ic) < 1. .and. nylocal > 0) then
          pym = voltfac*volt(2,ic)/dels(2,ic)
          denom = denom + (1.-dels(2,ic))/dels(2,ic)*dysqic
          ppp = min(ppp,dels(2,ic))
        endif
c       --- the point higher in y is inside the conductor
        if (0 < dels(3,ic) .and. dels(3,ic) < 1. .and. nylocal > 0) then
          pyp = voltfac*volt(3,ic)/dels(3,ic)
          denom = denom + (1.-dels(3,ic))/dels(3,ic)*dysqic
          ppp = min(ppp,dels(3,ic))
        endif
c       --- the point lower in z is inside the conductor
        if (0 < dels(4,ic) .and. dels(4,ic) < 1.) then
          pzm = voltfac*volt(4,ic)/dels(4,ic)
          denom = denom + (1.-dels(4,ic))/dels(4,ic)*dzsqic
          ppp = min(ppp,dels(4,ic))
        endif
c       --- the point higher in z is inside the conductor
        if (0 < dels(5,ic) .and. dels(5,ic) < 1.) then
          pzp = voltfac*volt(5,ic)/dels(5,ic)
          denom = denom + (1.-dels(5,ic))/dels(5,ic)*dzsqic
          ppp = min(ppp,dels(5,ic))
        endif
c       --- calculate the residual based on the boundary conditions
        if (ppp < 1.) then
          if (mgform2init) ppp = 1.
#ifdef WITHCHOMBO
          ppp = 1.
#endif
          res(ix,iy,iz) = ppp*(rho(ix,iy,iz)
     &           + (pxm+pxp)*dxsqics + (pym+pyp)*dysqics + (pzm+pzp)*dzsqics
     &           - phi(ix,iy,iz)*mgparam*denom)
        endif
      enddo
!$OMP END DO

      else if (icndbndy == 2) then
c     --- Quadratic interpolation

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

        if (ix == 0  .and. localbounds(0) == 0) cycle
        if (ix == nxlocal .and. localbounds(1) == 0) cycle
        if (iy == 0  .and. localbounds(2) == 0) cycle
        if (iy == nylocal .and. localbounds(3) == 0) cycle
        if (iz == 0  .and. localbounds(4) == 0) cycle
        if (iz == nzlocal .and. localbounds(5) == 0) cycle

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1
        izp1 = iz + 1
        izm1 = iz - 1

        if (nylocal == 0) then
          iym1 = iy
          iyp1 = iy
        endif

c       --- set temporaries with initial values
        pxm = phi(ixm1,iy  ,iz  )
        pxp = phi(ixp1,iy  ,iz  )
        pym = phi(ix  ,iym1,iz  )
        pyp = phi(ix  ,iyp1,iz  )
        pzm = phi(ix  ,iy  ,izm1)
        pzp = phi(ix  ,iy  ,izp1)
        denom = 1.
        ppp = 1.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0 .and. (mgform == 1 .or. mgform2init)) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

        call handlesubgrid2(0,1,pxm,pxp,dxsqic,ic,subgrid%nmax,dels,volt,voltfac,denom,ppp)
        if (nylocal > 0) call handlesubgrid2(2,3,pym,pyp,dysqic,ic,subgrid%nmax,dels,volt,voltfac,denom,ppp)
        call handlesubgrid2(4,5,pzm,pzp,dzsqic,ic,subgrid%nmax,dels,volt,voltfac,denom,ppp)

c       --- calculate the residual based on the boundary conditions
        if (ppp < 1.) then
          if (mgform2init) ppp = 1.
#ifdef WITHCHOMBO
          ppp = 1.
#endif
          res(ix,iy,iz) = ppp*(rho(ix,iy,iz)
     &           + (pxm+pxp)*dxsqics + (pym+pyp)*dysqics + (pzm+pzp)*dzsqics
     &           - phi(ix,iy,iz)*mgparam*denom)

          if (0. >= dels(0,ic) .and. dels(0,ic) >= -1.) res(ixm1,iy,iz) = 0.
          if (0. >= dels(1,ic) .and. dels(1,ic) >= -1.) res(ixp1,iy,iz) = 0.
          if (0. >= dels(2,ic) .and. dels(2,ic) >= -1. .and. nylocal > 0) res(ix,iym1,iz) = 0.
          if (0. >= dels(3,ic) .and. dels(3,ic) >= -1. .and. nylocal > 0) res(ix,iyp1,iz) = 0.
          if (0. >= dels(4,ic) .and. dels(4,ic) >= -1.) res(ix,iy,izm1) = 0.
          if (0. >= dels(5,ic) .and. dels(5,ic) >= -1.) res(ix,iy,izp1) = 0.

        endif
      enddo
!$OMP END DO

      endif

      if (lf3dtimesubs) timecondbndyres = timecondbndyres +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine residual3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                      dxsqi,dysqi,dzsqi,phi,rho,res,
     &                      mglevel,localbounds,mgparam,mgform,mgform2init,
     &                      lcndbndy,icndbndy,conductors,lprecalccoeffs)
      use Subtimersf3d
      use ConductorTypemodule
      integer(ISZ):: nxlocal,nylocal,nzlocal
      integer(ISZ):: delx,dely,delz
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: res(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      integer(ISZ):: mglevel,localbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
      logical(ISZ):: mgform2init
      logical(ISZ):: lcndbndy,lprecalccoeffs
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.
c taking into account the premultiplication of rho by
c   mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
c The resulting residual is also implicitly multiplied by the same constant.
c Note that then for restriction of the residual to a coarser grid, it must
c be scaled by the ratio old(dxsqi+dysqi+dzsqi)/new(dxsqi+dysqi+dzsqi).
c This is done in the restrict routine automatically.
c
c For internal conductors, the residual is set to zero inside and calculated
c using the modified form of the finite differenced Poisson's equation near
c the surface.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      real(kind=8):: const,dxsqic,dysqic,dzsqic
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      const = 0.5/(dxsqi+dysqi+dzsqi)
      dxsqic = dxsqi*mgparam*const
      dysqic = dysqi*mgparam*const
      dzsqic = dzsqi*mgparam*const

      ixmin = 0
      ixmax = nxlocal
      iymin = 0
      iymax = nylocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) == 0) ixmin = 1
      if (localbounds(1) == 0) ixmax = nxlocal - 1
      if (localbounds(2) == 0) iymin = 1
      if (localbounds(3) == 0) iymax = nylocal - 1
      if (localbounds(4) == 0) izmin = 1
      if (localbounds(5) == 0) izmax = nzlocal - 1

c     --- Calculate the residual.

      if (lprecalccoeffs) then

c       --- This zeroing out is done inside the loop so that the cache look
c       --- ups for the zeroing and the calculation are done at the same time.
        if (localbounds(4) == 0) then
c         --- At Dirichlet boundaries, zero out res at the boundary.
          res(:,:,-delz:0) = 0.
        else if (localbounds(4)  > 0 .and. delz > 1) then
c         --- Otherwise, only zero out the guard cells that won't be otherwise
c         --- filled in.
          res(:,:,-delz:-2) = 0
        endif

!$OMP DO
        do iz=izmin,izmax

          if (localbounds(2) == 0) then
            res(:,-dely:0,iz) = 0.
          else if (localbounds(2)  > 0 .and. dely > 1) then
            res(:,-dely:-2,iz) = 0.
          endif

          do iy=iymin,iymax

            if (localbounds(0) == 0) then
              res(-delx:0,iy,iz) = 0.
            else if (localbounds(0)  > 0 .and. delx > 1) then
              res(-delx:-2,iy,iz) = 0.
            endif

            do ix=ixmin,ixmax

              if (conductors%coeffs%data(7,ix,iy,iz) == 0.) then
                res(ix,iy,iz) = 0.
              else
                res(ix,iy,iz) = (rho(ix,iy,iz)
     &              + phi(ix-1,iy  ,iz  )*conductors%coeffs%data(0,ix,iy,iz)
     &              + phi(ix+1,iy  ,iz  )*conductors%coeffs%data(1,ix,iy,iz)
     &              + phi(ix  ,iy-1,iz  )*conductors%coeffs%data(2,ix,iy,iz)
     &              + phi(ix  ,iy+1,iz  )*conductors%coeffs%data(3,ix,iy,iz)
     &              + phi(ix  ,iy  ,iz-1)*conductors%coeffs%data(4,ix,iy,iz)
     &              + phi(ix  ,iy  ,iz+1)*conductors%coeffs%data(5,ix,iy,iz)
     &              + conductors%coeffs%data(6,ix,iy,iz)
     &              - phi(ix,iy,iz)*mgparam/conductors%coeffs%data(7,ix,iy,iz))
     &              *conductors%coeffs%data(8,ix,iy,iz)
              endif


            enddo

            if (localbounds(1) == 0) then
              res(nxlocal:nxlocal+delx,iy,iz) = 0.
            else if (localbounds(1)  > 0 .and. delx > 1) then
              res(nxlocal+2:nxlocal+delx,iy,iz) = 0.
            endif

          enddo

          if (localbounds(3) == 0) then
            res(:,nylocal:nylocal+dely,iz) = 0.
          else if (localbounds(3)  > 0 .and. dely > 1) then
            res(:,nylocal+2:nylocal+dely,iz) = 0.
          endif

        enddo
!$OMP END DO

        if (localbounds(5) == 0) then
          res(:,:,nzlocal:nzlocal+delz) = 0.
        else if (localbounds(5)  > 0 .and. delz > 1) then
          res(:,:,nzlocal+2:nzlocal+delz) = 0.
        endif

      else

c       --- This zeroing out is done inside the loop so that the cache look
c       --- ups for the zeroing and the calculation are done at the same time.
        if (localbounds(4) == 0) then
c         --- At Dirichlet boundaries, zero out res at the boundary.
          res(:,:,-delz:0) = 0.
        else if (localbounds(4)  > 0 .and. delz > 1) then
c         --- Otherwise, only zero out the guard cells that won't be otherwise
c         --- filled in.
          res(:,:,-delz:-2) = 0
        endif

!$OMP DO
        do iz=izmin,izmax

          if (localbounds(2) == 0) then
            res(:,-dely:0,iz) = 0.
          else if (localbounds(2)  > 0 .and. dely > 1) then
            res(:,-dely:-2,iz) = 0.
          endif

          do iy=iymin,iymax

            if (localbounds(0) == 0) then
              res(-delx:0,iy,iz) = 0.
            else if (localbounds(0)  > 0 .and. delx > 1) then
              res(-delx:-2,iy,iz) = 0.
            endif

            do ix=ixmin,ixmax

              res(ix,iy,iz) = rho(ix,iy,iz)
     &            +  (phi(ix-1,iy  ,iz  )+phi(ix+1,iy  ,iz  ))*dxsqic
     &            +  (phi(ix  ,iy-1,iz  )+phi(ix  ,iy+1,iz  ))*dysqic
     &            +  (phi(ix  ,iy  ,iz-1)+phi(ix  ,iy  ,iz+1))*dzsqic
     &            -  phi(ix,iy,iz)*mgparam

            enddo

            if (localbounds(1) == 0) then
              res(nxlocal:nxlocal+delx,iy,iz) = 0.
            else if (localbounds(1)  > 0 .and. delx > 1) then
              res(nxlocal+2:nxlocal+delx,iy,iz) = 0.
            endif

          enddo

          if (localbounds(3) == 0) then
            res(:,nylocal:nylocal+dely,iz) = 0.
          else if (localbounds(3)  > 0 .and. dely > 1) then
            res(:,nylocal+2:nylocal+dely,iz) = 0.
          endif

        enddo
!$OMP END DO

        if (localbounds(5) == 0) then
          res(:,:,nzlocal:nzlocal+delz) = 0.
        else if (localbounds(5)  > 0 .and. delz > 1) then
          res(:,:,nzlocal+2:nzlocal+delz) = 0.
        endif

c       --- Zero the residual inside conductors.
        call cond_potmgres(conductors%interior,nxlocal,nylocal,nzlocal,
     &                     delx,dely,delz,
     &                     res,mglevel,mgform,.false.)

        if (lcndbndy) then
c         --- Calculate the residual near the conductor.
          call condbndyres(conductors%evensubgrid,nxlocal,nylocal,nzlocal,
     &                     delx,dely,delz,
     &                     phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,localbounds,
     &                     mglevel,mgform,mgform2init,icndbndy)
          call condbndyres(conductors%oddsubgrid,nxlocal,nylocal,nzlocal,
     &                     delx,dely,delz,
     &                     phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,localbounds,
     &                     mglevel,mgform,mgform2init,icndbndy)
        endif

      endif

      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,
     &                               delx,dely,delz,
     &                               res,1,localbounds,.false.,.true.)

      if (lf3dtimesubs) timeresidual = timeresidual +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
c=============================================================================
      subroutine getnextcoarselevel3d(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                                nxcoarse,nycoarse,nzcoarse,
     &                                dxcoarse,dycoarse,dzcoarse)
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      integer(ISZ):: nxlocalcoarse,nylocalcoarse,nzlocalcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse

c     --- Calculate the size of the next coarsest grid.
c     --- The resulting number of coarse cells must be >= n/2.
c     --- Also, the number of cells in x and y must be even.
c     --- Note that nzlocalcoarse is only used if full-coarsening is done.

      real(kind=8):: mingridsize

      nxcoarse = (nx+1)/2
      if (mod(nxcoarse,2) == 1) nxcoarse = nxcoarse + 1

      nycoarse = (ny+1)/2
      if (mod(nycoarse,2) == 1) nycoarse = nycoarse + 1

      nzcoarse = (nz+1)/2
      if (mod(nzcoarse,2) == 1) nzcoarse = nzcoarse + 1

      mingridsize = min(dx,dy)
      mingridsize = min(mingridsize,dz)

      if (dx > 4./3.*mingridsize) then
        nxlocalcoarse = nxlocal
        nxcoarse = nx
      endif
      if (dy > 4./3.*mingridsize) then
        nylocalcoarse = nylocal
        nycoarse = ny
      endif
      if (dz > 4./3.*mingridsize) then
        nzlocalcoarse = nzlocal
        nzcoarse = nz
      endif

      dxcoarse = dx*nx/nxcoarse
      if (nycoarse == 0) then
        dycoarse = dxcoarse
      else
        dycoarse = dy*ny/nycoarse
      endif
      dzcoarse = dz*nz/nzcoarse

      return
      end
c=============================================================================
      subroutine getmglevels(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                       conductors,fsdecomp)
      use ConductorTypemodule
      use Decompositionmodule
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

      integer(ISZ):: iszzero = 0
      integer(ISZ):: mglevel
      real(kind=8):: mgscale
      type(Decomposition):: finedecomp
      type(Decomposition):: coarsedecomp
      integer(ISZ):: nxfine,nyfine,nzfine
      integer(ISZ):: nxlocalfine,nylocalfine,nzlocalfine
      real(kind=8):: dxfine,dyfine,dzfine
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse
      integer(ISZ):: nxlocalcoarse,nylocalcoarse,nzlocalcoarse
      real(kind=8):: lx,ly,lz
      real(kind=8):: dxcoarse,dycoarse,dzcoarse,cxf,cyf,czf

c     --- Only do this calculation if it is needed. This
c     --- is for parallel version. There, this is a global operation and
c     --- if it is not needed, it would be an unneccesary synchronization
c     --- point.  If any of the levels is greater than one, then this
c     --- operation has already been done for this set of conductors, so
c     --- skip it.
      if (maxval(conductors%levellx) > 1 .or.
     &    maxval(conductors%levelly) > 1 .or.
     &    maxval(conductors%levellz) > 1) return

      mglevel = 0
      mgscale = 1.

#ifdef MPIPARALLEL
c     --- Allocate temp space
      finedecomp%nxprocs = fsdecomp%nxprocs
      finedecomp%nyprocs = fsdecomp%nyprocs
      finedecomp%nzprocs = fsdecomp%nzprocs
      allocate(finedecomp%ix(0:fsdecomp%nxprocs-1))
      allocate(finedecomp%nx(0:fsdecomp%nxprocs-1))
      allocate(finedecomp%iy(0:fsdecomp%nyprocs-1))
      allocate(finedecomp%ny(0:fsdecomp%nyprocs-1))
      allocate(finedecomp%iz(0:fsdecomp%nzprocs-1))
      allocate(finedecomp%nz(0:fsdecomp%nzprocs-1))

      allocate(coarsedecomp%ix(0:fsdecomp%nxprocs-1))
      allocate(coarsedecomp%nx(0:fsdecomp%nxprocs-1))
      allocate(coarsedecomp%iy(0:fsdecomp%nyprocs-1))
      allocate(coarsedecomp%ny(0:fsdecomp%nyprocs-1))
      allocate(coarsedecomp%iz(0:fsdecomp%nzprocs-1))
      allocate(coarsedecomp%nz(0:fsdecomp%nzprocs-1))

      finedecomp%ix = fsdecomp%ix
      finedecomp%nx = fsdecomp%nx
      finedecomp%iy = fsdecomp%iy
      finedecomp%ny = fsdecomp%ny
      finedecomp%iz = fsdecomp%iz
      finedecomp%nz = fsdecomp%nz
#endif

      nxfine = nx
      nyfine = ny
      nzfine = nz
      nxlocalfine = nxlocal
      nylocalfine = nylocal
      nzlocalfine = nzlocal
      dxfine = dx
      dyfine = dy
      dzfine = dz
      lx = 1.
      ly = 1.
      lz = 1.

      do while (.true.)

#ifdef MPIPARALLEL
        conductors%levelix(mglevel) = finedecomp%ix(fsdecomp%ixproc)
        conductors%levelnx(mglevel) = finedecomp%nx(fsdecomp%ixproc)
        conductors%leveliy(mglevel) = finedecomp%iy(fsdecomp%iyproc)
        conductors%levelny(mglevel) = finedecomp%ny(fsdecomp%iyproc)
        conductors%leveliz(mglevel) = finedecomp%iz(fsdecomp%izproc)
        conductors%levelnz(mglevel) = finedecomp%nz(fsdecomp%izproc)
#else
        conductors%levelix(mglevel) = 0
        conductors%levelnx(mglevel) = nxlocalfine
        conductors%leveliy(mglevel) = 0
        conductors%levelny(mglevel) = nylocalfine
        conductors%leveliz(mglevel) = 0
        conductors%levelnz(mglevel) = nzlocalfine
#endif
        conductors%levellx(mglevel) = lx
        conductors%levelly(mglevel) = ly
        conductors%levellz(mglevel) = lz
        conductors%levels = mglevel + 1

        if (nxfine < 4 .or.
     &      (nyfine < 4 .and. nyfine > 0) .or.
     &      nzfine < 4 .or.
     &      mglevel == 100)
     &    exit

        call getnextcoarselevel3d(nxfine,nyfine,nzfine,
     &                            nxlocalfine,nylocalfine,nzlocalfine,
     &                            dxfine,dyfine,dzfine,
     &                            nxcoarse,nycoarse,nzcoarse,
     &                            dxcoarse,dycoarse,dzcoarse)

        cxf = dxcoarse/dxfine
        cyf = dycoarse/dyfine
        czf = dzcoarse/dzfine
c       --- This option is not supported
c       mgscale = mgscale*dxcoarse*dycoarse*dzcoarse/(dxfine*dyfine*dzfine)
        mgscale = 1.

#ifdef MPIPARALLEL
c       --- Find domains in coarser grid
        call mgdividenz(finedecomp,coarsedecomp,nxfine,nyfine,nzfine,
     &                  nxcoarse,nycoarse,nzcoarse,mgscale)
c       --- Set new value of nzlocal
        nxlocalcoarse = coarsedecomp%nx(fsdecomp%ixproc)
        nylocalcoarse = coarsedecomp%ny(fsdecomp%iyproc)
        nzlocalcoarse = coarsedecomp%nz(fsdecomp%izproc)
        finedecomp%ix = coarsedecomp%ix
        finedecomp%nx = coarsedecomp%nx
        finedecomp%iy = coarsedecomp%iy
        finedecomp%ny = coarsedecomp%ny
        finedecomp%iz = coarsedecomp%iz
        finedecomp%nz = coarsedecomp%nz
#else
        nxlocalcoarse = nxcoarse
        nylocalcoarse = nycoarse
        nzlocalcoarse = nzcoarse
#endif

        nxfine = nxcoarse
        nyfine = nycoarse
        nzfine = nzcoarse
        nxlocalfine = nxlocalcoarse
        nylocalfine = nylocalcoarse
        nzlocalfine = nzlocalcoarse
        dxfine = dxcoarse
        dyfine = dycoarse
        dzfine = dzcoarse
        mglevel = mglevel + 1
        lx = lx*cxf
        ly = ly*cyf
        lz = lz*czf

      enddo

#ifdef MPIPARALLEL
c     --- Deallocate temp space
      deallocate(finedecomp%ix)
      deallocate(finedecomp%nx)
      deallocate(finedecomp%iy)
      deallocate(finedecomp%ny)
      deallocate(finedecomp%iz)
      deallocate(finedecomp%nz)

      deallocate(coarsedecomp%ix)
      deallocate(coarsedecomp%nx)
      deallocate(coarsedecomp%iy)
      deallocate(coarsedecomp%ny)
      deallocate(coarsedecomp%iz)
      deallocate(coarsedecomp%nz)
#endif

      return
      end
c=============================================================================
      subroutine applyboundaryconditions3d(nx,ny,nz,delx,dely,delz,u,ncomp,
     &                                     bounds,lwithdirichlet,lzerodirichlet)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,delx,dely,delz,ncomp
      integer(ISZ):: bounds(0:5)
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz,ncomp)
      logical(ISZ):: lwithdirichlet,lzerodirichlet

      integer(ISZ):: ix,iy,iz
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (delx > 0) then
        if (lwithdirichlet) then
          if (bounds(0) == 0) then
            do ix=-1,-delx,-1
              u(ix,0:ny,0:nz,:)   = 2.*u(ix+1,0:ny,0:nz,:) - u(ix+2,0:ny,0:nz,:)
            enddo
          endif
          if (bounds(1) == 0) then
            do ix=nx+1,nx+delx
              u(ix,0:ny,0:nz,:) = 2.*u(ix-1,0:ny,0:nz,:) - u(ix-2,0:ny,0:nz,:)
            enddo
          endif
        else if (lzerodirichlet) then
          if (bounds(0) == 0) u(-delx:-1,:,:,:)   = 0.
          if (bounds(1) == 0) u(nx+1:nx+delx,:,:,:) = 0.
        endif
        if (bounds(0) == 1) u(-delx:-1,0:ny,0:nz,:) = u(delx:1:-1,0:ny,0:nz,:)
        if (bounds(1) == 1) u(nx+1:nx+delx,0:ny,0:nz,:) = u(nx-1:nx-delx:-1,0:ny,0:nz,:)
        if (bounds(0) == 2 .and. bounds(1) == 2) then
          u(-delx:-1,0:ny,0:nz,:) = u(nx-delx:nx-1,0:ny,0:nz,:)
          u(nx:nx+delx,0:ny,0:nz,:) = u(0:delx,0:ny,0:nz,:)
        endif
      endif

      if (dely > 0) then
        if (lwithdirichlet) then
          if (bounds(2) == 0) then
            do iy=-1,-dely,-1
              u(:,iy,0:nz,:)   = 2.*u(:,iy+1,0:nz,:) - u(:,iy+2,0:nz,:)
            enddo
          endif
          if (bounds(3) == 0) then
            do iy=ny+1,ny+dely
              u(:,iy,0:nz,:) = 2.*u(:,iy-1,0:nz,:) - u(:,iy-2,0:nz,:)
            enddo
          endif
        else if (lzerodirichlet) then
          if (bounds(2) == 0) u(:,-dely:-1,:,:)   = 0.
          if (bounds(3) == 0) u(:,ny+1:ny+dely,:,:) = 0.
        endif
        if (bounds(2) == 1) u(:,-dely:-1,0:nz,:)   = u(:,dely:1:-1,0:nz,:)
        if (bounds(3) == 1) u(:,ny+1:ny+dely,0:nz,:) = u(:,ny-1:ny-dely:-1,0:nz,:)
        if (bounds(2) == 2 .and. bounds(3) == 2) then
          u(:,-dely:-1,0:nz,:)   = u(:,ny-dely:ny-1,0:nz,:)
          u(:,ny:ny+dely,0:nz,:) = u(:,0:dely,0:nz,:)
        endif
      endif

      if (delz > 0) then
        if (lwithdirichlet) then
          if (bounds(4) == 0) then
            do iz=-1,-delz,-1
              u(:,:,iz,:)   = 2.*u(:,:,iz+1,:) - u(:,:,iz+2,:)
            enddo
          endif
          if (bounds(5) == 0) then
            do iz=nz+1,nz+delz
              u(:,:,iz,:) = 2.*u(:,:,iz-1,:) - u(:,:,iz-2,:)
            enddo
          endif
        else if (lzerodirichlet) then
          if (bounds(4) == 0) u(:,:,-delz:-1,:)   = 0.
          if (bounds(5) == 0) u(:,:,nz+1:nz+delz,:) = 0.
        endif
        if (bounds(4) == 1) u(:,:,-delz:-1,:)   = u(:,:,delz:1:-1,:)
        if (bounds(5) == 1) u(:,:,nz+1:nz+delz,:) = u(:,:,nz-1:nz-delz:-1,:)
        if (bounds(4) == 2 .and. bounds(5) == 2) then
          u(:,:,-delz:-1,:)   = u(:,:,nz-delz:nz-1,:)
          u(:,:,nz:nz+delz,:) = u(:,:,0:delz,:)
        endif
      endif

      if (lf3dtimesubs) timeapplyboundaryconditions3d =
     &                  timeapplyboundaryconditions3d + wtime() - substarttime

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      RECURSIVE subroutine checkconductors(nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                                     dx,dy,dz,conductors,fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Decompositionmodule
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

c Recursively calls the routine to generate the conductor data at the
c various mesh resolutions needed by the MG solver.

      integer(ISZ):: il,ic,ie,io,ix,iy,iz,nn,na,nmax
      integer(ISZ),allocatable:: isort(:)
      integer(ISZ):: nxlocal1,nylocal1,nzlocal1,levels
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Make sure the coarsening levels are setup.
      call getmglevels(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,conductors,
     &                 fsdecomp)

c     --- Call the work routine for each level
      do il=0,conductors%levels-1
        nxlocal1 = conductors%levelnx(il)
        nylocal1 = conductors%levelny(il)
        nzlocal1 = conductors%levelnz(il)
        call checkconductors_work(conductors%interior,
     &                            conductors%evensubgrid,conductors%oddsubgrid,
     &                            nxlocal1,nylocal1,nzlocal1,il)
      enddo

      conductors%interior%istart = 1
      conductors%evensubgrid%istart = 1
      conductors%oddsubgrid%istart = 1
      levels = conductors%levels

c     --- Sort the conductor data by level number
c     --- First, sort conductor points
      if (conductors%interior%n > 0) then
        allocate(isort(conductors%interior%n),stat=allocerror)
        if (allocerror /= 0) then
          print*,"checkconductors: allocation error ",allocerror,
     &           ": could not allocate isort to shape ",conductors%interior%n
          call kaboom("checkconductors: allocation error")
          return
        endif
        nn = conductors%interior%n
        call isortconductor(nn,conductors%interior%ilevel,isort,
     &                      conductors%interior%istart,levels)
        conductors%interior%n = conductors%interior%istart(levels) - 1
        na = conductors%interior%n
        call iswapconductor(na,isort,nn,3,conductors%interior%indx)
        call iswapconductor(na,isort,nn,1,conductors%interior%numb)
        call rswapconductor(na,isort,nn,1,conductors%interior%volt)
        call iswapconductor(na,isort,nn,1,conductors%interior%ilevel)
        deallocate(isort)
      endif

c     --- Sort even subgrid points
      if (conductors%evensubgrid%n > 0) then
        allocate(isort(conductors%evensubgrid%n),stat=allocerror)
        if (allocerror /= 0) then
          print*,"checkconductors: allocation error ",allocerror,
     &           ": could not allocate isort to shape ",conductors%evensubgrid%n
          call kaboom("checkconductors: allocation error")
          return
        endif
        nn = conductors%evensubgrid%n
        call isortconductor(nn,conductors%evensubgrid%ilevel,isort,
     &                      conductors%evensubgrid%istart,levels)
        conductors%evensubgrid%n = conductors%evensubgrid%istart(levels) - 1
        na = conductors%evensubgrid%n
        call iswapconductor(na,isort,nn,3,conductors%evensubgrid%indx)
        call rswapconductor(na,isort,nn,6,conductors%evensubgrid%dels)
        call rswapconductor(na,isort,nn,6,conductors%evensubgrid%volt)
        call iswapconductor(na,isort,nn,6,conductors%evensubgrid%numb)
        call iswapconductor(na,isort,nn,1,conductors%evensubgrid%ilevel)
        deallocate(isort)
      endif
    
c     --- Sort odd subgrid points
      if (conductors%oddsubgrid%n > 0) then
        allocate(isort(conductors%oddsubgrid%n),stat=allocerror)
        if (allocerror /= 0) then
          print*,"checkconductors: allocation error ",allocerror,
     &           ": could not allocate isort to shape ",conductors%oddsubgrid%n
          call kaboom("checkconductors: allocation error")
          return
        endif
        nn = conductors%oddsubgrid%n
        call isortconductor(nn,conductors%oddsubgrid%ilevel,isort,
     &                      conductors%oddsubgrid%istart,levels)
        conductors%oddsubgrid%n = conductors%oddsubgrid%istart(levels) - 1
        na = conductors%oddsubgrid%n
        call iswapconductor(na,isort,nn,3,conductors%oddsubgrid%indx)
        call rswapconductor(na,isort,nn,6,conductors%oddsubgrid%dels)
        call rswapconductor(na,isort,nn,6,conductors%oddsubgrid%volt)
        call iswapconductor(na,isort,nn,6,conductors%oddsubgrid%numb)
        call iswapconductor(na,isort,nn,1,conductors%oddsubgrid%ilevel)
        deallocate(isort)
      endif

      if (lf3dtimesubs) timecheckconductors = timecheckconductors +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine isortconductor(nc,condlevel,isort,istart,levels)
      integer(ISZ):: nc
      integer(ISZ):: condlevel(nc),isort(nc),istart(0:100)
      integer(ISZ):: levels

c     --- Sort the conductor data by level number
      integer(ISZ):: ilevel,nlevel,ii,ic

      ii = 1
      do ilevel=0,levels-1
        istart(ilevel) = ii
        do ic=1,nc
          if (condlevel(ic) == ilevel) then
            isort(ii) = ic
            ii = ii + 1
          endif
        enddo
      enddo
      istart(levels) = ii

      return
      end
c=============================================================================
      subroutine iswapconductor(nc,isort,nn,ni,icond)
      integer(ISZ):: nc,nn,ni,isort(nc),icond(ni,nn)
c Utility function for sort conductors
      integer(ISZ):: ic
      integer(ISZ),allocatable:: itemp(:,:)
      integer(ISZ):: allocerror
      allocate(itemp(ni,nn),stat=allocerror)
      if (allocerror /= 0) then
        print*,"iswapconductor: allocation error ",allocerror,
     &         ": could not allocate itemp to shape ",ni,nn
        call kaboom("iswapconductor: allocation error")
        return
      endif
      itemp = icond
      do ic=1,nc
        icond(:,ic) = itemp(:,isort(ic))
      enddo
      deallocate(itemp)
      return
      end
c=============================================================================
      subroutine rswapconductor(nc,isort,nn,ni,rcond)
      integer(ISZ):: nc,nn,ni,isort(nc)
      real(kind=8):: rcond(ni,nn)
c Utility function for sort conductors
      real(kind=8),allocatable:: rtemp(:,:)
      integer(ISZ):: ic
      integer(ISZ):: allocerror
      allocate(rtemp(ni,nn),stat=allocerror)
      if (allocerror /= 0) then
        print*,"rswapconductor: allocation error ",allocerror,
     &         ": could not allocate rtemp to shape ",ni,nn
        call kaboom("rswapconductor: allocation error")
        return
      endif
      rtemp = rcond
      do ic=1,nc
        rcond(:,ic) = rtemp(:,isort(ic))
      enddo
      deallocate(rtemp)
      return
      end
c=============================================================================
      subroutine checkconductors_work(interior,evensubgrid,oddsubgrid,
     &                                nxlocal,nylocal,nzlocal,mglevel)
      use ConductorInteriorTypemodule
      use ConductorSubGridTypemodule
      type(ConductorInteriorType):: interior
      type(ConductorSubGridType):: evensubgrid,oddsubgrid
      integer(ISZ):: nxlocal,nylocal,nzlocal,mglevel

c This checks the conductor dataset for consistency.
c  - removes any points outside of the mesh
c  - clean up data set, removing any subgrid points which may lie inside
c    of a conductor (those points are harmless to SOR but are damaging
c    to multigrid)
c  - removes any redundant subgrid points
c
c The notation for the 3D work grid is...
c   - all of the points inside of conductors are given a value larger
c     than the index of any subgrid point
c   - for all subgrid points, the index of that point is stored with a sign
c     attached - positive for even points, negative for odd points
c   - the value chosen for inside of conductors is large enough so that
c     is will not be the same as a subgrid point

      integer(ISZ),allocatable:: iii(:,:,:)
      integer(ISZ):: ic,i,ix,iy,iz,id,nmax
      integer(ISZ):: allocerror

      nmax = max(evensubgrid%nmax,oddsubgrid%nmax)

c     --- Set the conductor points.
      allocate(iii(0:nxlocal,0:nylocal,0:nzlocal),stat=allocerror)
      if (allocerror /= 0) then
        print*,"checkconductors_work: allocation error ",allocerror,
     &         ": could not allocate iii to shape ",nxlocal,nylocal,nzlocal
        call kaboom("checkconductors_work: allocation error")
        return
      endif
      iii = 0
      do ic=1,interior%n
        if (interior%ilevel(ic) /= mglevel) cycle
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        if (ix < 0 .or. nxlocal < ix .or.
     &      iy < 0 .or. nylocal < iy .or.
     &      iz < 0 .or. nzlocal < iz) then
          interior%ilevel(ic) = -1
          cycle
        endif
        iii(ix,iy,iz) = nmax + 1
      enddo

c     --- Scan through subgrid points:
c     ---   remove points which lie inside of a conductor
c     ---   register subgrid points in the work array iii
c     ---   check for redundant point (multiple points at grid location)
      do ic=1,evensubgrid%n

        if (evensubgrid%ilevel(ic) /= mglevel) cycle

        ix = evensubgrid%indx(0,ic)
        iy = evensubgrid%indx(1,ic)
        iz = evensubgrid%indx(2,ic)

        if (ix < 0 .or. nxlocal < ix .or.
     &      iy < 0 .or. nylocal < iy .or.
     &      iz < 0 .or. nzlocal < iz) then
          evensubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy,iz) == nmax+1) then
          evensubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If iii == 0, then this data point is outside any conductors and
c       --- is not redundant.
        if (iii(ix,iy,iz) == 0) then
          iii(ix,iy,iz) = ic
          cycle
        endif

c       iii(ix,iy,iz) < ncndmax+1
        i = iii(ix,iy,iz)

        if (i < 0) then
c         --- The point already there is odd so must be a different level
          iii(ix,iy,iz) = ic
          cycle
        endif

        if (evensubgrid%ilevel(ic) /= evensubgrid%ilevel(i)) then
c         --- The point already there is on a different level
          iii(ix,iy,iz) = ic
          cycle
        endif

c       --- There is another subgrid point here. Combine the data
c       --- of the two points.
c       --- For each direction, check if a conductor is nearer to this point
c       --- than the other point. If so, reset data for the other point.
        do id=0,5
          if (abs(evensubgrid%dels(id,ic)) < 1. .and.
     &        abs(evensubgrid%dels(id,ic)) < abs(evensubgrid%dels(id,i))) then
            evensubgrid%dels(id,i)  = evensubgrid%dels(id,ic)
            evensubgrid%volt(id,i) = evensubgrid%volt(id,ic)
            evensubgrid%numb(id,i) = evensubgrid%numb(id,ic)
          endif
        enddo
        evensubgrid%ilevel(ic) = -1

      enddo

c     --- Do the same for the odd conductor points.
c     --- Scan through subgrid points:
c     ---   remove points which lie inside of a conductor
c     ---   register subgrid points in the work array iii
c     ---   check for redundant point (multiple points at grid location)
      do ic=1,oddsubgrid%n

        if (oddsubgrid%ilevel(ic) /= mglevel) cycle

        ix = oddsubgrid%indx(0,ic)
        iy = oddsubgrid%indx(1,ic)
        iz = oddsubgrid%indx(2,ic)

        if (ix < 0 .or. nxlocal < ix .or.
     &      iy < 0 .or. nylocal < iy .or.
     &      iz < 0 .or. nzlocal < iz) then
          oddsubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy,iz) == nmax+1) then
          oddsubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If iii == 0, then this data point is outside any conductors and
c       --- is not redundant.
        if (iii(ix,iy,iz) == 0) then
          iii(ix,iy,iz) = ic
          cycle
        endif

c       iii(ix,iy,iz) < ncndmax+1
        i = iii(ix,iy,iz)

        if (i < 0) then
c         --- The point already there is odd so must be a different level
          iii(ix,iy,iz) = ic
          cycle
        endif

        if (oddsubgrid%ilevel(ic) /= oddsubgrid%ilevel(i)) then
c         --- The point already there is on a different level
          iii(ix,iy,iz) = ic
          cycle
        endif

c       --- There is another subgrid point here. Combine the data
c       --- of the two points.
c       --- For each direction, check if a conductor is nearer to this point
c       --- than the other point. If so, reset data for the other point.
        do id=0,5
          if (abs(oddsubgrid%dels(id,ic)) < 1. .and.
     &        abs(oddsubgrid%dels(id,ic)) < abs(oddsubgrid%dels(id,i))) then
            oddsubgrid%dels(id,i)  = oddsubgrid%dels(id,ic)
            oddsubgrid%volt(id,i) = oddsubgrid%volt(id,ic)
            oddsubgrid%numb(id,i) = oddsubgrid%numb(id,ic)
          endif
        enddo
        oddsubgrid%ilevel(ic) = -1

      enddo

      deallocate(iii)
      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine precalculatecoefficients3d(nxlocal,nylocal,nzlocal,dx,dy,dz,conductors,
     &                                      mgparam,localbounds)
      use ConductorTypemodule
      use MGCoefficientsmodule
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      type(ConductorType):: conductors
      real(kind=8):: mgparam
      integer(ISZ):: localbounds(0:5)

c Precalculate the finite difference coefficients. This includes the
c modified coefficients use with the subgrid scheme.

      integer(ISZ):: il,ic,ix,iy,iz
      real(kind=8):: dx1,dy1,dz1
      real(kind=8):: dxsqi,dysqi,dzsqi,const,dxsqic,dysqic,dzsqic
      real(kind=8):: dxsqicos,dysqicos,dzsqicos
      real(kind=8):: voltfac,denom,ppp,pxm,pxp,pym,pyp,pzm,pzp
      integer(ISZ):: nxlocal1,nylocal1,nzlocal1,nmax
      type(MGCoefficients),pointer:: coeffs,finercoeffs
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      do il=0,conductors%levels-1

c       --- Get the grid info for the coarse levels
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        dz1 = dz*conductors%levellz(il)
        nxlocal1 = conductors%levelnx(il)
        nylocal1 = conductors%levelny(il)
        nzlocal1 = conductors%levelnz(il)

c       --- Create the new instance which will hold the data. Setup the
c       --- appropriate pointers to the instance.
        coeffs => NewMGCoefficients()
        if (il == 0) then
          conductors%coeffs => coeffs
          NULLIFY(coeffs%finer)
        else
          finercoeffs%coarser => coeffs
          coeffs%finer => finercoeffs
        endif
        finercoeffs => coeffs

c       --- Allocate the data array.
        coeffs%nx = nxlocal1
        coeffs%ny = nylocal1
        coeffs%nz = nzlocal1
        call MGCoefficientschange(coeffs)

c       --- Calculate the finite difference coefficients
        dxsqi = 1./dx1**2
        dysqi = 1./dy1**2
        dzsqi = 1./dz1**2
        const = mgparam*0.5/(dxsqi + dysqi + dzsqi)
        dxsqic = dxsqi*const
        dysqic = dysqi*const
        dzsqic = dzsqi*const
        dxsqicos = dxsqic/mgparam
        dysqicos = dysqic/mgparam
        dzsqicos = dzsqic/mgparam
        coeffs%data(0:1,:,:,:) = dxsqic
        coeffs%data(2:3,:,:,:) = dysqic
        coeffs%data(4:5,:,:,:) = dzsqic
        coeffs%data(6,:,:,:) = 0. ! RHS
        coeffs%data(7,:,:,:) = 1. ! 1./denom
        coeffs%data(8,:,:,:) = 1. ! ppp

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (il == 0) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- For points interior to conductors, the coefficents are all zero
c       --- and the right hand side is the desired voltage or zero.
        do ic = conductors%interior%istart(il),
     &          conductors%interior%istart(il+1)-1
          ix = conductors%interior%indx(0,ic)
          iy = conductors%interior%indx(1,ic)
          iz = conductors%interior%indx(2,ic)
          coeffs%data(:,ix,iy,iz) = 0.
          if (il == 0) then
            coeffs%data(6,ix,iy,iz) = conductors%interior%volt(ic)
          endif
        enddo

        dels => conductors%evensubgrid%dels
        volt => conductors%evensubgrid%volt
        nmax = conductors%evensubgrid%nmax
        do ic = conductors%evensubgrid%istart(il),
     &          conductors%evensubgrid%istart(il+1)-1

          ix = conductors%evensubgrid%indx(0,ic)
          iy = conductors%evensubgrid%indx(1,ic)
          iz = conductors%evensubgrid%indx(2,ic)

          if (ix == 0       .and. localbounds(0) == 0) cycle
          if (ix == nxlocal .and. localbounds(1) == 0) cycle
          if (iy == 0       .and. localbounds(2) == 0) cycle
          if (iy == nylocal .and. localbounds(3) == 0) cycle
          if (iz == 0       .and. localbounds(4) == 0) cycle
          if (iz == nzlocal .and. localbounds(5) == 0) cycle

          pxm = 0.
          pxp = 0.
          pym = 0.
          pyp = 0.
          pzm = 0.
          pzp = 0.
          denom = 1.
          ppp = 1.

          call handlesubgrid2(0,1,
     &                        coeffs%data(0,ix,iy,iz),coeffs%data(1,ix,iy,iz),
     &                        dxsqicos,ic,nmax,
     &                        dels,volt,voltfac,denom,ppp)
          call handlesubgrid2(2,3,
     &                        coeffs%data(2,ix,iy,iz),coeffs%data(3,ix,iy,iz),
     &                        dysqicos,ic,nmax,
     &                        dels,volt,voltfac,denom,ppp)
          call handlesubgrid2(4,5,
     &                        coeffs%data(4,ix,iy,iz),coeffs%data(5,ix,iy,iz),
     &                        dzsqicos,ic,nmax,
     &                        dels,volt,voltfac,denom,ppp)

          if (ppp < 1.) then
            coeffs%data(7,ix,iy,iz) = 1./denom
            coeffs%data(8,ix,iy,iz) = ppp
            if (-1. < dels(0,ic) .and. dels(0,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(0,ix,iy,iz)*dxsqic
              coeffs%data(0,ix,iy,iz) = 0.
            endif
            if (-1. < dels(1,ic) .and. dels(1,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(1,ix,iy,iz)*dxsqic
              coeffs%data(1,ix,iy,iz) = 0.
            endif
            if (-1. < dels(2,ic) .and. dels(2,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(2,ix,iy,iz)*dysqic
              coeffs%data(2,ix,iy,iz) = 0.
            endif
            if (-1. < dels(3,ic) .and. dels(3,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(3,ix,iy,iz)*dysqic
              coeffs%data(3,ix,iy,iz) = 0.
            endif
            if (-1. < dels(4,ic) .and. dels(4,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(4,ix,iy,iz)*dzsqic
              coeffs%data(4,ix,iy,iz) = 0.
            endif
            if (-1. < dels(5,ic) .and. dels(5,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(5,ix,iy,iz)*dzsqic
              coeffs%data(5,ix,iy,iz) = 0.
            endif

          endif

        enddo

        dels => conductors%oddsubgrid%dels
        volt => conductors%oddsubgrid%volt
        nmax = conductors%oddsubgrid%nmax
        do ic = conductors%oddsubgrid%istart(il),
     &          conductors%oddsubgrid%istart(il+1)-1

          ix = conductors%oddsubgrid%indx(0,ic)
          iy = conductors%oddsubgrid%indx(1,ic)
          iz = conductors%oddsubgrid%indx(2,ic)

          if (ix == 0       .and. localbounds(0) == 0) cycle
          if (ix == nxlocal .and. localbounds(1) == 0) cycle
          if (iy == 0       .and. localbounds(2) == 0) cycle
          if (iy == nylocal .and. localbounds(3) == 0) cycle
          if (iz == 0       .and. localbounds(4) == 0) cycle
          if (iz == nzlocal .and. localbounds(5) == 0) cycle

          pxm = 0.
          pxp = 0.
          pym = 0.
          pyp = 0.
          pzm = 0.
          pzp = 0.
          denom = 1.
          ppp = 1.

          call handlesubgrid2(0,1,
     &                        coeffs%data(0,ix,iy,iz),coeffs%data(1,ix,iy,iz),
     &                        dxsqicos,ic,nmax,
     &                        dels,volt,voltfac,denom,ppp)
          call handlesubgrid2(2,3,
     &                        coeffs%data(2,ix,iy,iz),coeffs%data(3,ix,iy,iz),
     &                        dysqicos,ic,nmax,
     &                        dels,volt,voltfac,denom,ppp)
          call handlesubgrid2(4,5,
     &                        coeffs%data(4,ix,iy,iz),coeffs%data(5,ix,iy,iz),
     &                        dzsqicos,ic,nmax,
     &                        dels,volt,voltfac,denom,ppp)

          if (ppp < 1.) then
            coeffs%data(7,ix,iy,iz) = 1./denom
            coeffs%data(8,ix,iy,iz) = ppp
            if (-1. < dels(0,ic) .and. dels(0,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(0,ix,iy,iz)*dxsqic
              coeffs%data(0,ix,iy,iz) = 0.
            endif
            if (-1. < dels(1,ic) .and. dels(1,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(1,ix,iy,iz)*dxsqic
              coeffs%data(1,ix,iy,iz) = 0.
            endif
            if (-1. < dels(2,ic) .and. dels(2,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(2,ix,iy,iz)*dysqic
              coeffs%data(2,ix,iy,iz) = 0.
            endif
            if (-1. < dels(3,ic) .and. dels(3,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(3,ix,iy,iz)*dysqic
              coeffs%data(3,ix,iy,iz) = 0.
            endif
            if (-1. < dels(4,ic) .and. dels(4,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(4,ix,iy,iz)*dzsqic
              coeffs%data(4,ix,iy,iz) = 0.
            endif
            if (-1. < dels(5,ic) .and. dels(5,ic) < 1.) then
              coeffs%data(6,ix,iy,iz) = coeffs%data(6,ix,iy,iz) +
     &                                  coeffs%data(5,ix,iy,iz)*dzsqic
              coeffs%data(5,ix,iy,iz) = 0.
            endif

          endif

        enddo

      enddo

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine setupconductorfielddata(nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                                   dx,dy,dz,conductors,fsdecomp)
      use ConductorTypemodule
      use Decompositionmodule
      use Conductor3d,Only:lcorrectede
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

      integer(ISZ):: ix,iy,iz,nmax
      integer(ISZ):: ic,ie,io
      integer(ISZ):: allocerror
      logical(ISZ):: doalloc

      call checkconductors(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                     conductors,fsdecomp)

c     --- Initialize the data needed for the corrected E field calculation
      if ((lcorrectede .or. conductors%lcorrectede)) then

c       --- First, check if the icgrid array needs to be allocated or
c       --- reallocated.
        if (.not. ASSOCIATED(conductors%icgrid)) then
          doalloc = .true.
        else
          if (ALL((/nxlocal+1,nylocal+1,nzlocal+1/) == SHAPE(conductors%icgrid))) then
            doalloc = .false.
          else
            doalloc = .true.
            deallocate(conductors%icgrid)
          endif
        endif
        if (doalloc) then
          allocate(conductors%icgrid(0:nxlocal,0:nylocal,0:nzlocal),
     &             stat=allocerror)
          if (allocerror /= 0) then
            print*,"setupconductorfielddata: allocation error ",allocerror,
     &             ": could not allocate icgrid to shape ",1+nxlocal,1+nylocal,1+nzlocal
            call kaboom("setupconductorfielddata: allocation error")
            return
          endif
        endif

        nmax = max(conductors%evensubgrid%nmax,conductors%oddsubgrid%nmax)
        conductors%icgrid = 0
c       --- The E field at interior points is ignored for now
        do ic=conductors%interior%istart(0),conductors%interior%istart(1)-1
          ix = conductors%interior%indx(0,ic)
          iy = conductors%interior%indx(1,ic)
          iz = conductors%interior%indx(2,ic)
          conductors%icgrid(ix,iy,iz) = nmax + 1
        enddo
        do ie=conductors%evensubgrid%istart(0),conductors%evensubgrid%istart(1)-1
          ix = conductors%evensubgrid%indx(0,ie)
          iy = conductors%evensubgrid%indx(1,ie)
          iz = conductors%evensubgrid%indx(2,ie)
          conductors%icgrid(ix,iy,iz) = ie
        enddo
        do io=conductors%oddsubgrid%istart(0),conductors%oddsubgrid%istart(1)-1
          ix = conductors%oddsubgrid%indx(0,io)
          iy = conductors%oddsubgrid%indx(1,io)
          iz = conductors%oddsubgrid%indx(2,io)
          conductors%icgrid(ix,iy,iz) = -io
        enddo
      else

c       --- Since it is not needed, make sure that the memory is freed.
        if (ASSOCIATED(conductors%icgrid)) then
          deallocate(conductors%icgrid)
        endif

      endif

      return
      end
c=============================================================================
      subroutine getefieldatconductorsubgrid(conductors,phi,
     &                                       dx,dy,dz,nx,ny,nz,delx,dely,delz,
     &                                       bounds)
      use ConductorTypemodule
      use Subtimers3d
      use Conductor3d,Only:lcorrectede
      type(ConductorType):: conductors
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      integer(ISZ):: delx,dely,delz
      real(kind=8):: phi(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz)
      integer(ISZ):: bounds(0:5)

      integer(ISZ):: ic,ix,iy,iz,nmax
      integer(ISZ):: mx,my,mz,px,py,pz
      real(kind=8),pointer:: dels(:,:),volt(:,:)
      real(kind=8),pointer:: efield0(:,:),efieldd(:,:)

      if (.not. ((lcorrectede .or. conductors%lcorrectede) .and.
     &           conductors%evensubgrid%n > 0 .and.
     &           conductors%oddsubgrid%n > 0)) return

      nmax = max(conductors%evensubgrid%nmax,conductors%oddsubgrid%nmax)

      volt => conductors%evensubgrid%volt
      dels => conductors%evensubgrid%dels
      efield0 => conductors%evensubgrid%efield0
      efieldd => conductors%evensubgrid%efieldd
      do ic=conductors%evensubgrid%istart(0),conductors%evensubgrid%istart(1)-1
        ix = conductors%evensubgrid%indx(0,ic)
        iy = conductors%evensubgrid%indx(1,ic)
        iz = conductors%evensubgrid%indx(2,ic)
        call SETEE(0, -1, 0, 0, +1, 0, 0, 1./dx)
        call SETEE(1,  0,-1, 0,  0,+1, 0, 1./dy)
        call SETEE(2,  0, 0,-1,  0, 0,+1, 1./dz)
      enddo

      volt => conductors%oddsubgrid%volt
      dels => conductors%oddsubgrid%dels
      efield0 => conductors%oddsubgrid%efield0
      efieldd => conductors%oddsubgrid%efieldd
      do ic=conductors%oddsubgrid%istart(0),conductors%oddsubgrid%istart(1)-1
        ix = conductors%oddsubgrid%indx(0,ic)
        iy = conductors%oddsubgrid%indx(1,ic)
        iz = conductors%oddsubgrid%indx(2,ic)
        call SETEE(0, -1, 0, 0, +1, 0, 0, 1./dx)
        call SETEE(1,  0,-1, 0,  0,+1, 0, 1./dy)
        call SETEE(2,  0, 0,-1,  0, 0,+1, 1./dz)
      enddo

      return
      CONTAINS

c       ---------------------------------------------------------------------
        subroutine SETEE(id,mx,my,mz,px,py,pz,di)
        integer(ISZ):: id,mx,my,mz,px,py,pz
        real(kind=8):: di
        real(kind=8):: ppp,dd0,dd1,vv0,vv1
        integer(ISZ):: im,ip

        im = 2*id
        ip = 2*id+1

        if (0. < dels(im,ic) .and. dels(im,ic) <= 1.) then
          vv0 = volt(im,ic)
          dd0 = dels(im,ic)
        else
          vv0 = phi(ix+mx,iy+my,iz+mz)
          dd0 = 1.
        endif
        if (0. < dels(ip,ic) .and. dels(ip,ic) <= 1.) then
          vv1 = volt(ip,ic)
          dd1 = dels(ip,ic)
        else
          vv1 = phi(ix+px,iy+py,iz+pz)
          dd1 = 1.
        endif
        ppp = phi(ix,iy,iz)
        efield0(id,ic) = (dd1/(dd0 + dd1)*(vv0 - ppp)/dd0 +
     &                    dd0/(dd0 + dd1)*(ppp - vv1)/dd1)*di
        if (ix+mx >= 0 .and. iy+my >= 0 .and. iz+mz >= 0) then
          if (conductors%icgrid(ix+mx,iy+my,iz+mz) > nmax) then
            efieldd(im,ic) =
     &            ((vv0-ppp)/dd0 - (2.-dd0)/(2.+dd1)*(ppp-vv1)/dd1)*
     &            ((2.+dd1)/(dd1+dd0))*di
          endif
        endif
        if (ix+px <= nx .and. iy+py <= ny .and. iz+pz <= nz) then
          if (conductors%icgrid(ix+px,iy+py,iz+pz) > nmax) then
            efieldd(ip,ic) =
     &            ((ppp-vv1)/dd1 - (2.-dd1)/(2.+dd0)*(vv0-ppp)/dd0)*
     &            ((2.+dd0)/(dd0+dd1))*di
          endif
        endif

        return
        end subroutine SETEE

      end
c=============================================================================
      subroutine fixefieldatconductorpoints(conductors,field,dx,dy,dz,nx,ny,nz)
      use ConductorTypemodule
      use Subtimers3d
      use Conductor3d,Only:lcorrectede
      type(ConductorType):: conductors
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: field(0:2,0:nx,0:ny,0:nz)

c Doubles the field at conductor points. This is done to fix the finite
c difference of phi at conductor points.  The normal two point finite
c difference is E(i)=(phi(i-1)-phi(i+1))/2delta.  If the conductor covers
c cells i and i-1, the correct finite difference would be
c E=(phi(i) - phi(i+1))/delta.  Since in the conductor, phi(i-1) = phi(i),
c the E can be fixed by a simple multiplication by 2. Note that this will
c fail when a thin conductor only covers grid cell i and with cells i-1
c and i+1 in vacuum.  In that case, two seperate finite differences would
c need to be saved, one for each side of the conductor.  

      integer(ISZ):: ii,ix,iy,iz

      if (.not. ((lcorrectede .or. conductors%lcorrectede) .and.
     &           conductors%interior%n > 0)) return

      do ii=conductors%interior%istart(0),conductors%interior%istart(1)-1
        ix = conductors%interior%indx(0,ii)
        iy = conductors%interior%indx(1,ii)
        iz = conductors%interior%indx(2,ii)
        field(:,ix,iy,iz) = 2.*field(:,ix,iy,iz)
      enddo

      return
      end
c=========================================================================== 
      subroutine sete3dwithconductor(conductors,phi,selfe,np,xp,yp,zp,zgrid,
     &                  xmmin,ymmin,zmmin,
     &                  dx,dy,dz,nx,ny,nz,efetch,ex,ey,ez,l2symtry,l4symtry,
     &                  lcylindrical,delx,dely,delz)
      use ConductorTypemodule
      use Subtimers3d
      type(ConductorType):: conductors
      integer(ISZ):: np,nx,ny,nz
      integer(ISZ):: delx,dely,delz
      real(kind=8):: zgrid,xmmin,ymmin,zmmin,dx,dy,dz
      real(kind=8):: phi(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz)
      real(kind=8):: selfe(0:2,0:nx,0:ny,0:nz)
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      integer(ISZ):: efetch
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

c Gets self electric field for particles
c Note that the phi1d passed in is assumed to start at phi(0,0,-1).
c This is a change is behavior from what is was earlier. This change
c was made to avoid problems with bounds checking. phi needs to be accessed
c at the plane iz=-1, so when phi was passed in starting at phi(0,0,0), with
c bounds checking turned on, this would be caught (even though the code was
c technically correct).

c Algorithm notes: phi array is dimensioned (0:nx,0:ny,-1:nz+1) outside,
c but is made one dimensional in this routine
c so cell index into 1d phi array for vectorized deposition is:
c    i + j*(nx+1) + k*(nx+1)*(ny+1)
c The field is:
c    Ex = u0*v0*w0*ex(i  ,j  ,k  )
c       + u1*v0*w0*ex(i+1,j  ,k  )
c       + u0*v1*w0*ex(i  ,j+1,k  )
c       + ...

      integer(ISZ):: ip,i,j,k,inext,jnext,knext,nmax
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi
      real(kind=8):: u0,u1,v0,v1,w0,w1
      real(kind=8):: ext,eyt,ezt
      real(kind=8):: xi,yj,zk,x
      real(kind=8):: xinext,yjnext,zknext,xnext,ynext
      real(kind=8):: vv(0:1),dd(0:1)
      real(kind=8):: ee(0:2)
      real(kind=8):: eengp(0:2,0:1,0:1,0:1)
      real(kind=8):: wwsum(0:2)
      real(kind=8):: xsign,ysign,sx,sy
      real(kind=8),pointer:: dels(:),volt(:)

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- If the icgrid is not allocated, then the conductor data is to be
c     --- ignored, so call the base version of sete3d.
      if (.not. ASSOCIATED(conductors%icgrid)) then
        call sete3d(phi,selfe,np,xp,yp,zp,zgrid,xmmin,ymmin,zmmin,
     &              dx,dy,dz,nx,ny,nz,efetch,ex,ey,ez,l2symtry,l4symtry,
     &              lcylindrical,delx,dely,delz)
        return
      endif

c     --- Calculate some temporaries.
      tdxi = 1./(2.*dx)
      tdyi = 1./(2.*dy)
      tdzi = 1./(2.*dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      nmax = max(conductors%evensubgrid%nmax,conductors%oddsubgrid%nmax)

      if (efetch == 1 .or. efetch == 3 .or. efetch == 5 .or. efetch == 6) then

        if (lcylindrical) then
          xnext = sqrt(xp(1)**2 + yp(1)**2)
          ynext = ymmin
        else
          xnext = xp(1)
          ynext = yp(1)
        endif

        xinext = (xnext - xmmin)*dxi
        yjnext = (ynext - ymmin)*dyi
        zknext = (zp(1) - zgrid - zmmin)*dzi
        if (l4symtry) xinext = abs(xinext)
        if (l4symtry .or. l2symtry) yjnext = abs(yjnext)
        inext = xinext
        jnext = yjnext
        knext = zknext

        do ip = 1, np

          i = inext
          j = jnext
          k = knext
          x = xnext
          xi = xinext
          yj = yjnext
          zk = zknext
          if (ip < np) then
            if (lcylindrical) then
              xnext = sqrt(xp(ip+1)**2 + yp(ip+1)**2)
              ynext = ymmin
            else
              xnext = xp(ip+1)
              ynext = yp(ip+1)
            endif
            xinext = (xnext - xmmin)*dxi
            yjnext = (ynext - ymmin)*dyi
            zknext = (zp(ip+1) - zgrid - zmmin)*dzi
            if (l4symtry) xinext = abs(xinext)
            if (l4symtry .or. l2symtry) yjnext = abs(yjnext)
            inext = xinext
            jnext = yjnext
            knext = zknext
          endif

          if (xi < 0. .or. xi > nx .or.
     &        yj < 0. .or. yj > ny .or.
     &        zk < 0. .or. zk > nz) cycle

          u1 = xi - i
          v1 = yj - j
          w1 = zk - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

          if (conductors%icgrid(i  ,j  ,k  ) == 0 .and.
     &        conductors%icgrid(i+1,j  ,k  ) == 0 .and.
     &        conductors%icgrid(i  ,j+1,k  ) == 0 .and.
     &        conductors%icgrid(i+1,j+1,k  ) == 0 .and.
     &        conductors%icgrid(i  ,j  ,k+1) == 0 .and.
     &        conductors%icgrid(i+1,j  ,k+1) == 0 .and.
     &        conductors%icgrid(i  ,j+1,k+1) == 0 .and.
     &        conductors%icgrid(i+1,j+1,k+1) == 0) then

            if (efetch == 3) then

              ext = u0*v0*w0*selfe(0,i  ,j  ,k  )
     &            + u1*v0*w0*selfe(0,i+1,j  ,k  )
     &            + u0*v1*w0*selfe(0,i  ,j+1,k  )
     &            + u1*v1*w0*selfe(0,i+1,j+1,k  )
     &            + u0*v0*w1*selfe(0,i  ,j  ,k+1)
     &            + u1*v0*w1*selfe(0,i+1,j  ,k+1)
     &            + u0*v1*w1*selfe(0,i  ,j+1,k+1)
     &            + u1*v1*w1*selfe(0,i+1,j+1,k+1)

              eyt = u0*v0*w0*selfe(1,i  ,j  ,k  )
     &            + u1*v0*w0*selfe(1,i+1,j  ,k  )
     &            + u0*v1*w0*selfe(1,i  ,j+1,k  )
     &            + u1*v1*w0*selfe(1,i+1,j+1,k  )
     &            + u0*v0*w1*selfe(1,i  ,j  ,k+1)
     &            + u1*v0*w1*selfe(1,i+1,j  ,k+1)
     &            + u0*v1*w1*selfe(1,i  ,j+1,k+1)
     &            + u1*v1*w1*selfe(1,i+1,j+1,k+1)

              ezt = u0*v0*w0*selfe(2,i  ,j  ,k  )
     &            + u1*v0*w0*selfe(2,i+1,j  ,k  )
     &            + u0*v1*w0*selfe(2,i  ,j+1,k  )
     &            + u1*v1*w0*selfe(2,i+1,j+1,k  )
     &            + u0*v0*w1*selfe(2,i  ,j  ,k+1)
     &            + u1*v0*w1*selfe(2,i+1,j  ,k+1)
     &            + u0*v1*w1*selfe(2,i  ,j+1,k+1)
     &            + u1*v1*w1*selfe(2,i+1,j+1,k+1)

            else

              ext = u0*v0*w0*(phi(i-1,j  ,k  ) - phi(i+1,j  ,k  ))*tdxi +
     &              u1*v0*w0*(phi(i  ,j  ,k  ) - phi(i+2,j  ,k  ))*tdxi +
     &              u0*v1*w0*(phi(i-1,j+1,k  ) - phi(i+1,j+1,k  ))*tdxi +
     &              u1*v1*w0*(phi(i  ,j+1,k  ) - phi(i+2,j+1,k  ))*tdxi +
     &              u0*v0*w1*(phi(i-1,j  ,k+1) - phi(i+1,j  ,k+1))*tdxi +
     &              u1*v0*w1*(phi(i  ,j  ,k+1) - phi(i+2,j  ,k+1))*tdxi +
     &              u0*v1*w1*(phi(i-1,j+1,k+1) - phi(i+1,j+1,k+1))*tdxi +
     &              u1*v1*w1*(phi(i  ,j+1,k+1) - phi(i+2,j+1,k+1))*tdxi

              eyt = u0*v0*w0*(phi(i  ,j-1,k  ) - phi(i  ,j+1,k  ))*tdyi +
     &              u1*v0*w0*(phi(i+1,j-1,k  ) - phi(i+1,j+1,k  ))*tdyi +
     &              u0*v1*w0*(phi(i  ,j  ,k  ) - phi(i  ,j+2,k  ))*tdyi +
     &              u1*v1*w0*(phi(i+1,j  ,k  ) - phi(i+1,j+2,k  ))*tdyi +
     &              u0*v0*w1*(phi(i  ,j-1,k+1) - phi(i  ,j+1,k+1))*tdyi +
     &              u1*v0*w1*(phi(i+1,j-1,k+1) - phi(i+1,j+1,k+1))*tdyi +
     &              u0*v1*w1*(phi(i  ,j  ,k+1) - phi(i  ,j+2,k+1))*tdyi +
     &              u1*v1*w1*(phi(i+1,j  ,k+1) - phi(i+1,j+2,k+1))*tdyi 

              ezt = u0*v0*w0*(phi(i  ,j  ,k-1) - phi(i  ,j  ,k+1))*tdzi +
     &              u1*v0*w0*(phi(i+1,j  ,k-1) - phi(i+1,j  ,k+1))*tdzi +
     &              u0*v1*w0*(phi(i  ,j+1,k-1) - phi(i  ,j+1,k+1))*tdzi +
     &              u1*v1*w0*(phi(i+1,j+1,k-1) - phi(i+1,j+1,k+1))*tdzi +
     &              u0*v0*w1*(phi(i  ,j  ,k  ) - phi(i  ,j  ,k+2))*tdzi +
     &              u1*v0*w1*(phi(i+1,j  ,k  ) - phi(i+1,j  ,k+2))*tdzi +
     &              u0*v1*w1*(phi(i  ,j+1,k  ) - phi(i  ,j+1,k+2))*tdzi +
     &              u1*v1*w1*(phi(i+1,j+1,k  ) - phi(i+1,j+1,k+2))*tdzi
            endif

          else

            ee = 0.
            wwsum = 0.
            call GETEEFROME(0, 0,0,0, 1,0,0, dxi,u0*v0*w0,u1*v0*w0)
            call GETEEFROME(0, 0,1,0, 1,0,0, dxi,u0*v1*w0,u1*v1*w0)
            call GETEEFROME(0, 0,0,1, 1,0,0, dxi,u0*v0*w1,u1*v0*w1)
            call GETEEFROME(0, 0,1,1, 1,0,0, dxi,u0*v1*w1,u1*v1*w1)
            if (wwsum(0) .ne. 0.) wwsum(0) = 1./wwsum(0)

            call GETEEFROME(1, 0,0,0, 0,1,0, dyi,u0*v0*w0,u0*v1*w0)
            call GETEEFROME(1, 1,0,0, 0,1,0, dyi,u1*v0*w0,u1*v1*w0)
            call GETEEFROME(1, 0,0,1, 0,1,0, dyi,u0*v0*w1,u0*v1*w1)
            call GETEEFROME(1, 1,0,1, 0,1,0, dyi,u1*v0*w1,u1*v1*w1)
            if (wwsum(1) .ne. 0.) wwsum(1) = 1./wwsum(1)

            call GETEEFROME(2, 0,0,0, 0,0,1, dzi,u0*v0*w0,u0*v0*w1)
            call GETEEFROME(2, 1,0,0, 0,0,1, dzi,u1*v0*w0,u1*v0*w1)
            call GETEEFROME(2, 0,1,0, 0,0,1, dzi,u0*v1*w0,u0*v1*w1)
            call GETEEFROME(2, 1,1,0, 0,0,1, dzi,u1*v1*w0,u1*v1*w1)
            if (wwsum(2) .ne. 0.) wwsum(2) = 1./wwsum(2)

            ext = ee(0)*wwsum(0)
            eyt = ee(1)*wwsum(1)
            ezt = ee(2)*wwsum(2)

          endif

          if (xp(ip) < 0. .and. l4symtry) ext = -ext
          if (yp(ip) < 0. .and. (l4symtry .or. l2symtry)) ext = -eyt

          if (lcylindrical) then
            if (x > 0.) then
              eyt = ext*yp(ip)/x
              ext = ext*xp(ip)/x
            else
              ext = ext
              eyt = 0.
            endif
          endif

          ex(ip) = ex(ip) + ext
          ey(ip) = ey(ip) + eyt
          ez(ip) = ez(ip) + ezt

        enddo

      elseif (efetch == 4) then
c       --- Energy conserving

        if (lcylindrical) then
          xnext = sqrt(xp(1)**2 + yp(1)**2)
          ynext = ymmin
        else
          xnext = xp(1)
          ynext = yp(1)
        endif

        if (l4symtry) then
          sx = -1
        else
          sx = 1
        endif
        if (l4symtry .or. l2symtry) then
          sy = -1
        else
          sy = 1
        endif

        if (l4symtry) then
          xinext = (abs(xnext) - xmmin)*dxi
          yjnext = (abs(ynext) - ymmin)*dyi
        elseif (l2symtry) then
          xinext = (xnext - xmmin)*dxi
          yjnext = (abs(ynext) - ymmin)*dyi
        else
          xinext = (xnext - xmmin)*dxi
          yjnext = (ynext - ymmin)*dyi
        endif
        zknext = (zp(1) - zgrid - zmmin)*dzi
        inext = xinext
        jnext = yjnext
        knext = zknext

        do ip = 1, np

          i = inext
          j = jnext
          k = knext
          xi = xinext
          yj = yjnext
          zk = zknext
          if (ip < np) then
            if (lcylindrical) then
              xnext = sqrt(xp(ip+1)**2 + yp(ip+1)**2)
              ynext = ymmin
            else
              xnext = xp(ip+1)
              ynext = yp(ip+1)
            endif
            if (l4symtry) then
              xinext = (abs(xnext) - xmmin)*dxi
              yjnext = (abs(ynext) - ymmin)*dyi
            elseif (l2symtry) then
              xinext = (xnext - xmmin)*dxi
              yjnext = (abs(ynext) - ymmin)*dyi
            else
              xinext = (xnext - xmmin)*dxi
              yjnext = (ynext - ymmin)*dyi
            endif
            zknext = (zp(ip+1) - zgrid - zmmin)*dzi
            inext = xinext
            jnext = yjnext
            knext = zknext
          endif

          if (xi < 0. .or. xi > nx .or.
     &        yj < 0. .or. yj > ny .or.
     &        zk < 0. .or. zk > nz) cycle

          u1 = xi - i
          v1 = yj - j
          w1 = zk - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

          call GETEENGP(0,1,0,0,0,0,0,v0*w0,dxi)
          call GETEENGP(0,1,0,0,0,1,0,v1*w0,dxi)
          call GETEENGP(0,1,0,0,0,0,1,v0*w1,dxi)
          call GETEENGP(0,1,0,0,0,1,1,v1*w1,dxi)
          ext = (v0*w0*eengp(0,0,0,0) + v1*w0*eengp(0,0,1,0) +
     &           v0*w1*eengp(0,0,0,1) + v1*w1*eengp(0,0,1,1))
          if (xp(ip) < 0) ext = sx*ext

          call GETEENGP(1,0,1,0,0,0,0,u0*w0,dyi)
          call GETEENGP(1,0,1,0,1,0,0,u1*w0,dyi)
          call GETEENGP(1,0,1,0,0,0,1,u0*w1,dyi)
          call GETEENGP(1,0,1,0,1,0,1,u1*w1,dyi)
          eyt = (u0*w0*eengp(1,0,0,0) + u1*w0*eengp(1,1,0,0) +
     &           u0*w1*eengp(1,0,0,1) + u1*w1*eengp(1,1,0,1))
          if (yp(ip) < 0) eyt = sy*eyt
              
          call GETEENGP(2,0,0,1,0,0,0,u0*v0,dzi)
          call GETEENGP(2,0,0,1,1,0,0,u1*v0,dzi)
          call GETEENGP(2,0,0,1,0,1,0,u0*v1,dzi)
          call GETEENGP(2,0,0,1,1,1,0,u1*v1,dzi)
          ezt = (u0*v0*eengp(2,0,0,0) + u1*v0*eengp(2,1,0,0) +
     &           u0*v1*eengp(2,0,1,0) + u1*v1*eengp(2,1,1,0))

          if (lcylindrical) then
            if (x > 0.) then
              eyt = ext*yp(ip)/x
              ext = ext*xp(ip)/x
            else
              ext = ext
              eyt = 0.
            endif
          endif

          ex(ip) = ex(ip) + ext
          ey(ip) = ey(ip) + eyt
          ez(ip) = ez(ip) + ezt

        enddo

      endif

!$OMP MASTER
c     if (lf3dtimesubs) timesete3d = timesete3d + wtime() - substarttime
!$OMP END MASTER

      return
      CONTAINS

c       ---------------------------------------------------------------------
        subroutine GETEEFROME(id,ox,oy,oz,px,py,pz,di,wm,wp)
        integer(ISZ):: id,ox,oy,oz,px,py,pz
        real(kind=8):: di,wm,wp
        integer(ISZ):: ic
        real(kind=8):: ppp
        ic = conductors%icgrid(i+ox,j+oy,k+oz)
        if (ic == 0) then
          if (efetch == 3) then
            ee(id) = ee(id) + wm*selfe(id,i+ox,j+oy,k+oz)
          else
            ee(id) = ee(id) + wm*(  phi(i+ox-px,j+oy-py,k+oz-pz)
     &                            - phi(i+ox+px,j+oy+py,k+oz+pz))*0.5*di
          endif
          wwsum(id) = wwsum(id) + wm
        else if (ic < nmax+1) then
          if (ic < 0) then
            ee(id) = ee(id) + wm*conductors%oddsubgrid%efield0(id,-ic)
          else
            ee(id) = ee(id) + wm*conductors%evensubgrid%efield0(id,ic)
          endif
          wwsum(id) = wwsum(id) + wm
          if (conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz) > nmax) then
            if (ic < 0) then
              ee(id) = ee(id) + wp*conductors%oddsubgrid%efieldd(2*id+1,-ic)
            else
              ee(id) = ee(id) + wp*conductors%evensubgrid%efieldd(2*id+1,ic)
            endif
            wwsum(id) = wwsum(id) + wp
          endif
        endif
        ic = conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz)
        if (ic == 0) then
          if (efetch == 3) then
            ee(id) = ee(id) + wp*selfe(id,i+ox+px,j+oy+py,k+oz+pz)
          else
            ee(id) = ee(id) + wp*(  phi(i+ox     ,j+oy     ,k+oz)
     &                            - phi(i+ox+2*px,j+oy+2*py,k+oz+2*pz))*0.5*di
          endif
          wwsum(id) = wwsum(id) + wp
        else if (ic < nmax+1) then
          if (ic < 0) then
            ee(id) = ee(id) + wp*conductors%oddsubgrid%efield0(id,-ic)
          else
            ee(id) = ee(id) + wp*conductors%evensubgrid%efield0(id,ic)
          endif
          wwsum(id) = wwsum(id) + wp
          if (conductors%icgrid(i+ox,j+oy,k+oz) > nmax) then
            if (ic < 0) then
              ee(id) = ee(id) + wm*conductors%oddsubgrid%efieldd(2*id,-ic)
            else
              ee(id) = ee(id) + wm*conductors%evensubgrid%efieldd(2*id,ic)
            endif
            wwsum(id) = wwsum(id) + wm
          endif
        endif

        return
        end subroutine GETEEFROME

c       ---------------------------------------------------------------------
        subroutine GETEENGP(id,px,py,pz,ox,oy,oz,dw,dd)
        integer(ISZ):: id,px,py,pz,ox,oy,oz
        real(kind=8):: dw,dd
        real(kind=8):: vngp,dngp
        integer(ISZ):: ic
        if (conductors%icgrid(i+ox,j+oy,k+oz) == 0 .or.
     &      conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz) == 0) then
          eengp(id,ox,oy,oz) = (phi(i+ox,j+oy,k+oz) - phi(i+ox+px,j+oy+py,k+oz+pz))*dd
        elseif (conductors%icgrid(i+ox,j+oy,k+oz) < nmax+1) then
          ic = conductors%icgrid(i+ox,j+oy,k+oz)
          if (ic < 0) then
            dels => conductors%oddsubgrid%dels(:,-ic)
            volt => conductors%oddsubgrid%volt(:,-ic)
          else
            dels => conductors%evensubgrid%dels(:,ic)
            volt => conductors%evensubgrid%volt(:,ic)
          endif
          if (dels(2*id+1+1) <= 1.) then
            vngp = volt(2*id+1+1)
            dngp = dels(2*id+1+1)
          else
            vngp = phi(i+ox+px,j+oy+py,k+oz+pz)
            dngp = 1
          endif
          eengp(id,ox,oy,oz) = (phi(i+ox,j+oy,k+oz) - vngp)*dd/dngp
        elseif (conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz) < nmax+1) then
          ic = conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz)
          if (ic < 0) then
            dels => conductors%oddsubgrid%dels(:,-ic)
            volt => conductors%oddsubgrid%volt(:,-ic)
          else
            dels => conductors%evensubgrid%dels(:,ic)
            volt => conductors%evensubgrid%volt(:,ic)
          endif
          if (dels(2*id+1) <= 1.) then
            vngp = volt(2*id+1)
            dngp = dels(2*id+1)
          else
            vngp = phi(i+ox,j+oy,k+oz)
            dngp = 1.
          endif
          eengp(id,ox,oy,oz) = (vngp - phi(i+ox+px,j+oy+py,k+oz+pz))*dd/dngp
        else
          eengp(id,ox,oy,oz) = 0.
        endif

        end subroutine GETEENGP
      end
c=========================================================================== 
      subroutine sete3dongridwithconductor(conductors,phi,
     &                                     dx,dy,dz,
     &                                     nxlocal,nylocal,nzlocal,ex,ey,ez,
     &                                     delx,dely,delz,bounds)
      use GlobalVars,Only: dirichlet,neumann,periodic
      use ConductorTypemodule
      use Subtimers3d
      type(ConductorType):: conductors
      integer(ISZ):: nxlocal,nylocal,nzlocal
      integer(ISZ):: delx,dely,delz
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      real(kind=8):: ex(-1:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: ey(0:nxlocal,-1:nylocal,0:nzlocal)
      real(kind=8):: ez(0:nxlocal,0:nylocal,-1:nzlocal)
      integer(ISZ):: bounds(0:5)

c Gets self electric field normals for at dual cell grid points.
c The dual cell extends from i-1/2 to i+1/2 in each plane.
c For each cell (i,j,k), this gets
c  Ex(i+1/2,j,k) = (phi(i,j,k) - phi(i+1,j,k))/dx
c  Ey(i,j+1/2,k) = (phi(i,j,k) - phi(i,j+1,k))/dy
c  Ez(i,j,k+1/2) = (phi(i,j,k) - phi(i,j,k+1))/dz

      integer(ISZ):: i,j,k,nmax
      real(kind=8):: dxi,dyi,dzi
      real(kind=8):: eengp(0:2,0:1,0:1,0:1)
      real(kind=8),pointer:: dels(:),volt(:)

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- If the icgrid is not allocated, then the conductor data is to be
c     --- ignored, so call the base version of sete3d.
      if (.not. ASSOCIATED(conductors%icgrid)) then
        call kaboom("sete3dongridwithconductor: icgrid must be initialized")
        return
      endif

c     --- Calculate some temporaries.
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz
      nmax = max(conductors%evensubgrid%nmax,conductors%oddsubgrid%nmax)

c     --- Do the calculation in the bulk of the grid.
      do k=0,nzlocal-1
        do j=0,nylocal-1
          do i=0,nxlocal-1

            call GETEENGP(0,1,0,0,0,0,0,1.,dxi)
            ex(i,j,k) = eengp(0,0,0,0)

            call GETEENGP(1,0,1,0,0,0,0,1.,dyi)
            ey(i,j,k) = eengp(1,0,0,0)
              
            call GETEENGP(2,0,0,1,0,0,0,1.,dzi)
            ez(i,j,k) = eengp(2,0,0,0)

          enddo
        enddo
      enddo

c     --- Handle the values at nxlocal, nylocal, and nzlocal. These must be
c     --- done separately. For example, for i=nxlocal, ex cannot be calculated
c     --- since that would need the info at ix=nxlocal+1. phi is known there
c     --- but not icgrid.
      do k=0,nzlocal-1
        do j=0,nylocal-1
          i = nxlocal

            call GETEENGP(1,0,1,0,0,0,0,1.,dyi)
            ey(i,j,k) = eengp(1,0,0,0)
              
            call GETEENGP(2,0,0,1,0,0,0,1.,dzi)
            ez(i,j,k) = eengp(2,0,0,0)

        enddo
      enddo

      do k=0,nzlocal-1
        j = nylocal
          do i=0,nxlocal-1

            call GETEENGP(0,1,0,0,0,0,0,1.,dxi)
            ex(i,j,k) = eengp(0,0,0,0)

            call GETEENGP(2,0,0,1,0,0,0,1.,dzi)
            ez(i,j,k) = eengp(2,0,0,0)

          enddo
      enddo

      k = nzlocal
        do j=0,nylocal-1
          do i=0,nxlocal-1

            call GETEENGP(0,1,0,0,0,0,0,1.,dxi)
            ex(i,j,k) = eengp(0,0,0,0)

            call GETEENGP(1,0,1,0,0,0,0,1.,dyi)
            ey(i,j,k) = eengp(1,0,0,0)
              
          enddo
        enddo

c     --- Now, set the E fields in the guard cells by appropriate copying
c     --- or setting.
      select case (bounds(0))
        case (0) !(dirichlet)
          ex(-1,:,:) = 0
        case (1) !(neumann)
          ex(-1,:,:) = ex(0,:,:)
        case (2) !(periodic)
          ex(-1,:,:) = ex(nxlocal-1,:,:)
      end select
      select case (bounds(1))
        case (0) !(dirichlet)
          ex(nxlocal,:,:) = 0
        case (1) !(neumann)
          ex(nxlocal,:,:) = ex(nxlocal-1,:,:)
        case (2) !(periodic)
          ex(nxlocal,:,:) = ex(0,:,:)
      end select

      select case (bounds(2))
        case (0) !(dirichlet)
          ey(:,-1,:) = 0
        case (1) !(neumann)
          ey(:,-1,:) = ey(:,0,:)
        case (2) !(periodic)
          ey(:,-1,:) = ey(:,nylocal-1,:)
      end select
      select case (bounds(3))
        case (0) !(dirichlet)
          ey(:,nylocal,:) = 0
        case (1) !(neumann)
          ey(:,nylocal,:) = ey(:,nylocal-1,:)
        case (2) !(periodic)
          ey(:,nylocal,:) = ey(:,0,:)
      end select

      select case (bounds(4))
        case (0) !(dirichlet)
          ez(:,:,-1) = 0
        case (1) !(neumann)
          ez(:,:,-1) = ez(:,:,0)
        case (2) !(periodic)
          ez(:,:,-1) = ez(:,:,nzlocal-1)
      end select
      select case (bounds(5))
        case (0) !(dirichlet)
          ez(:,:,nzlocal) = 0
        case (1) !(neumann)
          ez(:,:,nzlocal) = ez(:,:,nzlocal-1)
        case (2) !(periodic)
          ez(:,:,nzlocal) = ez(:,:,0)
      end select

!$OMP MASTER
c     if (lf3dtimesubs) timesete3d = timesete3d + wtime() - substarttime
!$OMP END MASTER

      return
      CONTAINS

c       ---------------------------------------------------------------------
        subroutine GETEENGP(id,px,py,pz,ox,oy,oz,dw,dd)
        integer(ISZ):: id,px,py,pz,ox,oy,oz
        real(kind=8):: dw,dd
        real(kind=8):: vngp,dngp
        integer(ISZ):: ic
        if (conductors%icgrid(i+ox,j+oy,k+oz) == 0 .or.
     &      conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz) == 0) then
          eengp(id,ox,oy,oz) = (phi(i+ox,j+oy,k+oz) - phi(i+ox+px,j+oy+py,k+oz+pz))*dd
        elseif (conductors%icgrid(i+ox,j+oy,k+oz) < nmax+1) then
          ic = conductors%icgrid(i+ox,j+oy,k+oz)
          if (ic < 0) then
            dels => conductors%oddsubgrid%dels(:,-ic)
            volt => conductors%oddsubgrid%volt(:,-ic)
          else
            dels => conductors%evensubgrid%dels(:,ic)
            volt => conductors%evensubgrid%volt(:,ic)
          endif
          if (dels(2*id+1+1) <= 1.) then
            vngp = volt(2*id+1+1)
            dngp = dels(2*id+1+1)
          else
            vngp = phi(i+ox+px,j+oy+py,k+oz+pz)
            dngp = 1
          endif
          eengp(id,ox,oy,oz) = (phi(i+ox,j+oy,k+oz) - vngp)*dd/dngp
        elseif (conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz) < nmax+1) then
          ic = conductors%icgrid(i+ox+px,j+oy+py,k+oz+pz)
          if (ic < 0) then
            dels => conductors%oddsubgrid%dels(:,-ic)
            volt => conductors%oddsubgrid%volt(:,-ic)
          else
            dels => conductors%evensubgrid%dels(:,ic)
            volt => conductors%evensubgrid%volt(:,ic)
          endif
          if (dels(2*id+1) <= 1.) then
            vngp = volt(2*id+1)
            dngp = dels(2*id+1)
          else
            vngp = phi(i+ox,j+oy,k+oz)
            dngp = 1.
          endif
          eengp(id,ox,oy,oz) = (vngp - phi(i+ox+px,j+oy+py,k+oz+pz))*dd/dngp
        else
          eengp(id,ox,oy,oz) = 0.
        endif

        end subroutine GETEENGP
      end
c=========================================================================== 
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine gathersourcefromchild(source,nc,nn,childsource,cnn,
     &                                 l,u,fulllower,childlower,childupper,
     &                                 r,weights,radius,cradius,lcylindrical)
      use Subtimersf3d
      integer(ISZ):: nc,nn(0:2),nextra,cnn(0:2)
      integer(ISZ):: l(0:2),u(0:2),fulllower(0:2)
      integer(ISZ):: childlower(0:2),childupper(0:2)
      integer(ISZ):: r(0:2)
      real(kind=8):: source(0:nc-1,0:nn(0),0:nn(1),0:nn(2))
      real(kind=8):: childsource(0:nc-1,0:cnn(0),0:cnn(1),0:cnn(2))
      real(kind=8):: weights(-r(0)+1:r(0)-1,-r(1)+1:r(1)-1,-r(2)+1:r(2)-1)
      real(kind=8):: radius(0:nn(0)),cradius(0:cnn(0))
      logical(ISZ):: lcylindrical

      integer(ISZ):: ix,iy,iz,icx,icy,icz,iwx,iwy,iwz
      integer(ISZ):: ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2
      integer(ISZ):: iwx1,iwy1,iwz1,iwx2,iwy2,iwz2,ic
      real(kind=8):: dr,cdr,onethird,rr,crr
      real(kind=8):: wsource
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      dr = radius(1) - radius(0)
      cdr = cradius(1) - cradius(0)
      onethird = 1./3.

      do iz=l(2),u(2)
        iz0 = iz - fulllower(2)
        iz1 = max(iz*r(2) - r(2) + 1,childlower(2))
        iz2 = min(iz*r(2) + r(2) - 1,childupper(2))
        iwz1 = iz1 - iz*r(2)
        iwz2 = iz2 - iz*r(2)
        icz = iz*r(2) - childlower(2)
        do iy=l(1),u(1)
          iy0 = iy - fulllower(1)
          iy1 = max(iy*r(1) - r(1) + 1,childlower(1))
          iy2 = min(iy*r(1) + r(1) - 1,childupper(1))
          iwy1 = iy1 - iy*r(1)
          iwy2 = iy2 - iy*r(1)
          icy = iy*r(1) - childlower(1)
          do ix=l(0),u(0)
            ix0 = ix - fulllower(0)
            ix1 = max(ix*r(0) - r(0) + 1,childlower(0))
            ix2 = min(ix*r(0) + r(0) - 1,childupper(0))
            iwx1 = ix1 - ix*r(0)
            iwx2 = ix2 - ix*r(0)
            icx = ix*r(0) - childlower(0)

            if (lcylindrical) then
              if (radius(ix0) == 0.) then
                rr = 3./(dr)
              else
                rr = 1./(2.*radius(ix0))
              endif
            endif

            if (nc == 1) then
c             --- This case is separated out for optimization
              do iwz=iwz1,iwz2
                do iwy=iwy1,iwy2
                  do iwx=iwx1,iwx2
                    wsource=weights(iwx,iwy,iwz)*childsource(0,icx+iwx,icy+iwy,icz+iwz)
                    if (lcylindrical) then
c                     --- Calculate volumes. Note that the factors of pi and dz
c                     --- cancel out and so are left out. Also one
c                     --- factor of dr has beed factored out. Finally, note
c                     --- that crr is calculated as the volume to avoid a
c                     --- division.
                      if (cradius(icx+iwx) == 0.) then
                        crr = cdr*onethird
                      else
                        crr = 2.*cradius(icx+iwx)
                      endif
                      wsource = wsource*rr*crr
                    endif

                    source(0,ix0,iy0,iz0) = source(0,ix0,iy0,iz0) + wsource

                  enddo
                enddo
              enddo

            else

              do iwz=iwz1,iwz2
                do iwy=iwy1,iwy2
                  do iwx=iwx1,iwx2
                    do ic=0,nc-1
                      wsource=weights(iwx,iwy,iwz)*childsource(ic,icx+iwx,icy+iwy,icz+iwz)
                      if (lcylindrical) then
c                       --- Calculate volumes. Note that the factors of pi and
c                       --- dz cancel out and so are left out. Also one
c                       --- factor of dr has beed factored out. Finally, note
c                       --- that crr is calculated as the volume to avoid a
c                       --- division.
                        if (cradius(icx+iwx) == 0.) then
                          crr = cdr*onethird
                        else
                          crr = 2.*cradius(icx+iwx)
                        endif
                        wsource = wsource*r(2)*rr*crr
                      endif

                      source(ic,ix0,iy0,iz0) = source(ic,ix0,iy0,iz0) + wsource
                    enddo
                  enddo
                enddo
              enddo

            endif

          enddo
        enddo
      enddo

      if (lf3dtimesubs) timegathersourcefromchild = timegathersourcefromchild +
     &                                           wtime() - substarttime

      return
      end
c=============================================================================
      subroutine gatherpotentialfromparents(potential,nc,ng,nn,l,u,fulllower,
     &                                parentpotential,pnn,parentfulllower,r)
      use Subtimersf3d
      integer(ISZ):: nc,ng(0:2),nn(0:2),pnn(0:2)
      integer(ISZ):: r(0:2)
      integer(ISZ):: l(0:2),u(0:2),fulllower(0:2)
      integer(ISZ):: parentfulllower(0:2),parentupper(0:2)
      real(kind=8):: potential(0:nc-1,-ng(0):nn(0)+ng(0),-ng(1):nn(1)+ng(1),-ng(2):nn(2)+ng(2))
      real(kind=8):: parentpotential(0:nc-1,-ng(0):pnn(0)+ng(0),-ng(1):pnn(1)+ng(1),-ng(2):pnn(2)+ng(2))

      integer(ISZ):: ix0,iy0,iz0
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixp1,ixp2,iyp1,iyp2,izp1,izp2
      integer(ISZ),allocatable:: ixp1a(:),ixp2a(:)
      integer(ISZ),allocatable:: iyp1a(:),iyp2a(:)
      integer(ISZ),allocatable:: izp1a(:),izp2a(:)
      integer(ISZ):: ixp,iyp,izp,ic
      real(kind=8):: wx(0:1),wy(0:1),wz(0:1)
      real(kind=8),allocatable:: wxa(:,:),wya(:,:),wza(:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      allocate(ixp1a(0:nn(0)),ixp2a(0:nn(0)))
      allocate(iyp1a(0:nn(1)),iyp2a(0:nn(1)))
      allocate(izp1a(0:nn(2)),izp2a(0:nn(2)))
      allocate(wxa(0:1,0:nn(0)),wya(0:1,0:nn(1)),wza(0:1,0:nn(2)))

c     --- Precalculate the loop limits and weights.
      do iz0=l(2),u(2)
        iz = iz0 - fulllower(2)
        izp1a(iz) = int(iz0/r(2) - parentfulllower(2) + 1) - 1
        wza(1,iz) = 1.*mod(iz0,r(2))/r(2)
        wza(0,iz) = 1. - wza(1,iz)
        if (wza(1,iz) == 0.) then
          izp2a(iz) = izp1a(iz)
        else
          izp2a(iz) = izp1a(iz) + 1
        endif
      enddo

      do iy0=l(1),u(1)
        iy = iy0 - fulllower(1)
        iyp1a(iy) = int(iy0/r(1) - parentfulllower(1) + 1) - 1
        wya(1,iy) = 1.*mod(iy0,r(1))/r(1)
        wya(0,iy) = 1. - wya(1,iy)
        if (wya(1,iy) == 0.) then
          iyp2a(iy) = iyp1a(iy)
        else
          iyp2a(iy) = iyp1a(iy) + 1
        endif
      enddo

      do ix0=l(0),u(0)
        ix = ix0 - fulllower(0)
        ixp1a(ix) = int(ix0/r(0) - parentfulllower(0) + 1) - 1
        wxa(1,ix) = 1.*mod(ix0,r(0))/r(0)
        wxa(0,ix) = 1. - wxa(1,ix)
        if (wxa(1,ix) == 0.) then
          ixp2a(ix) = ixp1a(ix)
        else
          ixp2a(ix) = ixp1a(ix) + 1
        endif
      enddo

      do iz0=l(2),u(2)
        iz = iz0 - fulllower(2)
        izp1 = izp1a(iz)
        izp2 = izp2a(iz)
        wz(:) = wza(:,iz)

        do iy0=l(1),u(1)
          iy = iy0 - fulllower(1)
          iyp1 = iyp1a(iy)
          iyp2 = iyp2a(iy)
          wy(:) = wya(:,iy)

          do ix0=l(0),u(0)
            ix = ix0 - fulllower(0)
            ixp1 = ixp1a(ix)
            ixp2 = ixp2a(ix)
            wx(:) = wxa(:,ix)

            potential(:,ix,iy,iz) = 0.

            do izp=izp1,izp2
              do iyp=iyp1,iyp2
                do ixp=ixp1,ixp2

                  do ic=0,nc-1
                    potential(ic,ix,iy,iz) = potential(ic,ix,iy,iz) + 
     &                wx(ixp-ixp1)*wy(iyp-iyp1)*wz(izp-izp1)*
     &                parentpotential(ic,ixp,iyp,izp)
                  enddo

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      deallocate(ixp1a,ixp2a)
      deallocate(iyp1a,iyp2a)
      deallocate(izp1a,izp2a)
      deallocate(wxa,wya,wza)

      if (lf3dtimesubs) timegatherpotentialfromparents =
     &                                   timegatherpotentialfromparents +
     &                                   wtime() - substarttime

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigrid2dsolve(iwhich,nx,nz,nxlocal,nzlocal,dx,dz,
     &                            phi,rho,bounds,
     &                            xmminlocal,zbeam,zgrid,
     &                            mgparam,mgiters,mgmaxiters,
     &                            mgmaxlevels,mgerror,mgtol,mgverbose,
     &                            downpasses,uppasses,
     &                            lcndbndy,laddconductor,icndbndy,
     &                            gridmode,conductors,lrz,
     &                            fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      use Decompositionmodule
      integer(ISZ):: iwhich
      integer(ISZ):: nx,nz,nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nzlocal)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmminlocal,zbeam,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels,mgverbose
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      logical(ISZ):: lrz
      type(Decomposition):: fsdecomp

c Use the multigrid method for solving Poisson's equation on a 2D Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution.
c
c When the grid cells are rectangular, semi-coarsening is done until the
c grid cell dimensions are roughly equal. Roughly equal means that
c   2/3 dx < dz < 4/3 dx
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3. Currently, it is still
c assumed that dx ~ dy and that semi-coarsening is not needed transversely.

      integer(ISZ):: i,k,ix,iz
      real(kind=8),allocatable:: phisave(:,:)

      integer(ISZ):: localbounds(0:5)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Note that nx and nzlocal do not need to be even

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

c     --- Make sure that no decomposition in y is being done.
      if (fsdecomp%nyprocs > 1) then
        call kaboom("multigrid2dsolve: decomposition in y is not supported")
        return
      endif

      localbounds = bounds
#ifdef MPIPARALLEL
      if (fsdecomp%ix(fsdecomp%ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(fsdecomp%ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iz(fsdecomp%izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(fsdecomp%izproc)+nzlocal < nz) localbounds(5) = -1
#endif

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
      call checkconductors(nx,0,nz,nxlocal,0,nzlocal,dx,dx,dz,
     &                     conductors,fsdecomp)

!$OMP PARALLEL
!$OMP&PRIVATE(i,ix,iz)

c     --- Preset rho to increase performance (reducing the number of
c     --- multiplies in the main SOR sweep loop).
!OMP DO
      rho = rho/eps0
!OMP END DO

#ifdef MPIPARALLEL
c     --- These calls break the parallel field solver
c     call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
c    &                             1,0,1,-1,0,localbounds,fsdecomp)
c     call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif

c     --- Make sure guard planes have sensible values before beginning.
      call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,localbounds,
     &                               .false.,.false.)

      allocate(phisave(-1:nxlocal+1,-1:nzlocal+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid2dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nxlocal,nzlocal
        call kaboom("multigrid2dsolve: allocation error")
        return
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
        call vcycle2d(0,1.,nx,nz,nxlocal,nzlocal,dx,dz,phi,rho,
     &                bounds,mgparam,mgmaxlevels,
     &                downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &                xmminlocal,lrz,fsdecomp)

c       --- Calculate the change in phi.
        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nzlocal
          do ix=0,nxlocal
            mgerror = max(mgerror,abs(phisave(ix,iz) - phi(ix,iz)))
          enddo
        enddo
!$OMP END DO

#ifdef MPIPARALLEL
        if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
c         --- calculate global sorerror
          call parallelmaxrealarray(mgerror,1)
        endif
#endif

      enddo

c     --- Set boundary conditions. This is only really needed for the
c     --- Dirichlet boundaries, but this is convenient to call.
      call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,
     &                               localbounds,.true.,.false.)

      if (mgverbose>=1 .or. (mgverbose>=0 .and. mgerror > mgtol)) then
c       --- Make a print out.
        if (mgerror > mgtol) then
          call remark("Multigrid2d: Maximum number of iterations reached")
        endif
        write(errline,20) mgerror,mgiters
  20    format("Multigrid2d: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif

      deallocate(phisave)

c     --- Restore rho
      rho = rho*eps0

!$OMP END PARALLEL

      if (lf3dtimesubs) timemultigrid2dsolve = timemultigrid2dsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycle2d(mglevel,mgscale,nx,nz,nxlocal,nzlocal,dx,dz,
     &                              phi,rho,globalbounds,
     &                              mgparam,
     &                              mgmaxlevels,downpasses,uppasses,
     &                              lcndbndy,icndbndy,conductors,xmminlocal,lrz,
     &                              fsdecomp)
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use formggetarraysuminterface
      use Decompositionmodule
      integer(ISZ):: mglevel
      real(kind=8):: mgscale
      integer(ISZ):: nx,nz,nxlocal,nzlocal
      real(kind=8):: dx,dz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      real(kind=8):: xmminlocal
      logical(ISZ):: lcndbndy,lrz
      integer(ISZ):: icndbndy
      type(Decomposition):: fsdecomp

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: dxsqi,dzsqi
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: res(:,:)
      integer(ISZ):: i,iszone=1
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nylocal
      integer(ISZ):: nxlocalcoarse,nzlocalcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dzcoarsesqi
      real(kind=8):: xminodx,xmminlocalcoarse
      real(kind=8):: mgscalecoarse
      integer(ISZ):: ixproc,izproc
      integer(ISZ):: localbounds(0:5),localboundsc(0:5)
      integer(ISZ):: lxoffsetall(0:fsdecomp%nxprocs-1)
      integer(ISZ):: rxoffsetall(0:fsdecomp%nxprocs-1)
      integer(ISZ):: lzoffsetall(0:fsdecomp%nzprocs-1)
      integer(ISZ):: rzoffsetall(0:fsdecomp%nzprocs-1)
      integer(ISZ):: lxoffset,rxoffset
      integer(ISZ):: lzoffset,rzoffset
      type(Decomposition):: coarsedecomp
      integer(ISZ):: allocerror
      integer(ISZ):: delx,delz
      logical(ISZ):: lpe0
      real(kind=8):: sss(2)

      dxsqi = 1./dx**2
      dzsqi = 1./dz**2
      xminodx = xmminlocal/dx

      localbounds = globalbounds

#ifdef MPIPARALLEL
      ixproc = fsdecomp%ixproc
      izproc = fsdecomp%izproc
      if (fsdecomp%ix(ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iz(izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(izproc)+nzlocal < nz) localbounds(5) = -1
#endif

      if (lprintmgarraysumdiagnostic) then
#ifdef MPIPARALLEL
        lpe0=(fsdecomp%ixproc==0.and.fsdecomp%izproc==0)
#else
        lpe0 = .true.
#endif
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V1 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V1 rho",mglevel,sss
      endif

c     --- Do initial SOR passes.
      do i=1,downpasses
        call sorpass2d(mglevel,nxlocal,nzlocal,phi,rho,dxsqi,dzsqi,
     &                 xminodx,localbounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,lrz,fsdecomp)
      enddo

      if (lprintmgarraysumdiagnostic) then
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V2 rho",mglevel,sss
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. nz >= 4 .and.
     &    mglevel < mgmaxlevels) then

        delx = 1
        delz = 1
        if (nx > nxlocal) delx = 3
        if (nz > nzlocal) delz = 3

        allocate(res(-delx:nxlocal+delx,
     &               -delz:nzlocal+delz),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2d: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nxlocal,nzlocal
          call kaboom("vcycle2d: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residual2d(nxlocal,nzlocal,delx,delz,dxsqi,dzsqi,
     &                  xminodx,lrz,phi,rho,res,mglevel,localbounds,
     &                  lcndbndy,icndbndy,conductors)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,res,
     &                               delx,0,delz,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,0,nzlocal,res,delx,0,delz,-3,-1,fsdecomp)
#endif
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocal,0,nzlocal,delx,0,delz,res,fsdecomp)
          if (lpe0) print*,"V3 res",mglevel,sss
        endif

c       --- Note that some y quantities are included as dummies since the
c       --- routine will change the values.
        call getnextcoarselevel3d(nx,0,nz,nxlocal,nylocal,nzlocal,dx,dx,dz,
     &                            nxcoarse,nycoarse,nzcoarse,
     &                            dxcoarse,dycoarse,dzcoarse)

        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2
c       --- This option is not supported
c       mgscalecoarse = mgscale*dxcoarse*dzcoarse/(dx*dz)
        mgscalecoarse = 1.

        localboundsc = globalbounds

#ifdef MPIPARALLEL
        coarsedecomp%nxglobal = nxcoarse
        coarsedecomp%nyglobal = 0
        coarsedecomp%nzglobal = nzcoarse
        coarsedecomp%mpi_comm_x = fsdecomp%mpi_comm_x
        coarsedecomp%mpi_comm_z = fsdecomp%mpi_comm_z
        coarsedecomp%ixproc = fsdecomp%ixproc
        coarsedecomp%iyproc = fsdecomp%iyproc
        coarsedecomp%izproc = fsdecomp%izproc
        coarsedecomp%nxprocs = fsdecomp%nxprocs
        coarsedecomp%nyprocs = fsdecomp%nyprocs
        coarsedecomp%nzprocs = fsdecomp%nzprocs
        allocate(coarsedecomp%ix(0:fsdecomp%nxprocs-1))
        allocate(coarsedecomp%nx(0:fsdecomp%nxprocs-1))
        allocate(coarsedecomp%iy(0:fsdecomp%nyprocs-1))
        allocate(coarsedecomp%ny(0:fsdecomp%nyprocs-1))
        allocate(coarsedecomp%iz(0:fsdecomp%nzprocs-1))
        allocate(coarsedecomp%nz(0:fsdecomp%nzprocs-1))
c       --- Find domains in coarser grid
        call mgdividenz(fsdecomp,coarsedecomp,nx,0,nz,
     &                  nxcoarse,0,nzcoarse,mgscale)
c       --- Reset value to corrected one
        nxlocalcoarse = coarsedecomp%nx(ixproc)
        nzlocalcoarse = coarsedecomp%nz(izproc)
c       --- Difference between starts and ends of coarse and fine grids.
c       --- Should only be in the range 0-2.
        lxoffsetall = (nxcoarse*fsdecomp%ix-nx*coarsedecomp%ix)
        rxoffsetall = (nx*(coarsedecomp%ix + coarsedecomp%nx) -
     &                 nxcoarse*(fsdecomp%ix + fsdecomp%nx))
        lzoffsetall = (nzcoarse*fsdecomp%iz-nz*coarsedecomp%iz)
        rzoffsetall = (nz*(coarsedecomp%iz + coarsedecomp%nz) -
     &                 nzcoarse*(fsdecomp%iz + fsdecomp%nz))
c       --- Note that the lzoffsetall and rzoffsetall can only be used in
c       --- MPIPARALLEL sections since they will be unallocated in the
c       --- serial code. So, separate scalars are used in code which is
c       --- used in the serial version.
        lxoffset = lxoffsetall(ixproc)
        rxoffset = rxoffsetall(ixproc)
        lzoffset = lzoffsetall(izproc)
        rzoffset = rzoffsetall(izproc)
        if (coarsedecomp%ix(ixproc) > 0) localboundsc(0) = -1
        if (coarsedecomp%ix(ixproc)+nxlocalcoarse < nxcoarse) localboundsc(1) = -1
        if (coarsedecomp%iz(izproc) > 0) localboundsc(4) = -1
        if (coarsedecomp%iz(izproc)+nzlocalcoarse < nzcoarse) localboundsc(5) = -1
c       --- Calculate the xmminlocal of the coarse grid
        xmminlocalcoarse = xmminlocal
     &                          - fsdecomp%ix(fsdecomp%ixproc)*dx
     &                          + coarsedecomp%ix(fsdecomp%ixproc)*dxcoarse
#else
        nxlocalcoarse = nxcoarse
        nzlocalcoarse = nzcoarse
        lxoffset = 0
        rxoffset = 0
        lzoffset = 0
        rzoffset = 0
        xmminlocalcoarse = xmminlocal
#endif

c       --- Alloate new work space
        allocate(phicoarse(-1:nxlocalcoarse+1,-1:nzlocalcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2d: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcycle2d: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxlocalcoarse,0:nzlocalcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2d: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcycle2d: allocation error")
          return
        endif

        rhocoarse = 0.
        phicoarse = 0.

c       --- Restriction - note that no scaling factor is needed
        call restrict2d(nx,nz,nxlocal,nzlocal,res,delx,delz,
     &                  nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                  rhocoarse,0,0,
     &                  localbounds,localboundsc,lxoffset,lzoffset,.true.)

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,0,0,0,rhocoarse,coarsedecomp)
          if (lpe0) print*,"V3 rhocoarse",mglevel,sss
        endif

c       --- Continue at the next coarsest level.
        call vcycle2d(mglevel+iszone,mgscalecoarse,nxcoarse,nzcoarse,
     &                nxlocalcoarse,nzlocalcoarse,
     &                dxcoarse,dzcoarse,phicoarse,rhocoarse,
     &                globalbounds,mgparam,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors,xmminlocalcoarse,lrz,
     &                coarsedecomp)

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,1,0,1,phicoarse,coarsedecomp)
          if (lpe0) print*,"V4 phicoarse",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand2d(nx,nz,nxlocal,nzlocal,phi,
     &                nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,phicoarse,
     &                localbounds,lxoffset,lzoffset)
        call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,localbounds,
     &                                 .false.,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,-1,-1,localbounds,fsdecomp)
#endif
        deallocate(phicoarse,rhocoarse)
        deallocate(res)

#ifdef MPIPARALLEL
        deallocate(coarsedecomp%ix)
        deallocate(coarsedecomp%nx)
        deallocate(coarsedecomp%iy)
        deallocate(coarsedecomp%ny)
        deallocate(coarsedecomp%iz)
        deallocate(coarsedecomp%nz)
#endif

      endif

      if (lprintmgarraysumdiagnostic) then
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,0,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V5 phi",mglevel,sss
      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call sorpass2d(mglevel,nxlocal,nzlocal,phi,rho,
     &                 dxsqi,dzsqi,xminodx,
     &                 localbounds,mgparam,
     &                 lcndbndy,icndbndy,conductors,lrz,fsdecomp)
      enddo

      return
      end
c=============================================================================
      subroutine sorpass2d(mglevel,nxlocal,nzlocal,phi,rho,
     &                     dxsqi,dzsqi,xminodx,localbounds,
     &                     mgparam,lcndbndy,icndbndy,conductors,lrz,fsdecomp)
      use Constant
      use ConductorTypemodule
      use Decompositionmodule
      integer(ISZ):: mglevel,nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,xminodx
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      logical(ISZ):: lrz
      type(Decomposition):: fsdecomp

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.  It makes calls to the routines which specify internal
c conductors.
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

c     --- Put desired potential onto conductors in phi array.
      call cond_potmg(conductors%interior,nxlocal,0,nzlocal,1,0,1,
     &                phi,mglevel,1,.false.)
      call condbndymgint(conductors%evensubgrid,nxlocal,0,nzlocal,1,0,1,
     &                   phi,localbounds,mglevel,icndbndy)
      call condbndymgint(conductors%oddsubgrid ,nxlocal,0,nzlocal,1,0,1,
     &                   phi,localbounds,mglevel,icndbndy)

c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      parity = + fsdecomp%ix(fsdecomp%ixproc)
     &         + fsdecomp%iz(fsdecomp%izproc)
      s_parity = mod(parity,2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

        call sorhalfpass2d(parity,mglevel,nxlocal,nzlocal,phi,rho,
     &                     dxsqi,dzsqi,xminodx,
     &                     localbounds,mgparam,
     &                     lcndbndy,icndbndy,conductors,lrz)

        call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,localbounds,
     &                                 .false.,.false.)

#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,-1,0,localbounds,fsdecomp)
        if (parity == s_parity) then
          call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,0,0,fsdecomp)
        else
          call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,0,fsdecomp)
        endif
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- This doesn't seem to be needed.
c     --- Exchange phi in the z guard planes
c     call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif

      return
      end
c=============================================================================
      subroutine sorhalfpass2d(parity,mglevel,nxlocal,nzlocal,phi,rho,
     &                         dxsqi,dzsqi,xminodx,localbounds,
     &                         mgparam,lcndbndy,icndbndy,conductors,lrz)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: parity,mglevel,nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,xminodx
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: mgparam
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      logical(ISZ):: lrz

c This routine does one pass of point SOR with either even or odd
c parity. It makes calls to the routines which specify internal
c conductors.

      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ix,iz,ic,i1,i2

c     --- Set min and max indices for phi array.
c     --- If using Dirichlet boundary conditions, do not solve for the
c     --- potential on the grid edge.
      ixmin = 0
      ixmax = nxlocal
      izmin = 0
      izmax = nzlocal
      if (lrz .and. xminodx == 0.) ixmin = 0
      if (localbounds(0) < 1) ixmin = 1
      if (localbounds(1) < 1) ixmax = nxlocal - 1
      if (localbounds(4) < 1) izmin = 1
      if (localbounds(5) < 1) izmax = nzlocal - 1

c     --- Save values just outside conductor surfaces. Only save phi at the
c     --- subgrid points which are to be used at the current level of
c     --- grid refinement.
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

c     --- guard planes in z are already set

c     --- Loop over the bulk of the array, excluding boundaries.
      call mgsor_loop2d(ixmin,ixmax,izmin,izmax,parity,mgparam,
     &                  nxlocal,nzlocal,phi,rho,dxsqi,dzsqi,xminodx,lrz)

c     --- Apply altered difference equation to the points near the
c     --- surface of the conductor boundaries.
      if (lcndbndy) then
        if (parity == 0) then
         call condbndymg2d(conductors%evensubgrid,nxlocal,nzlocal,phi,rho,
     &                     dxsqi,dzsqi,xminodx,lrz,
     &                     mgparam,localbounds,mglevel,icndbndy)
        endif
        if (parity == 1) then
         call condbndymg2d(conductors%oddsubgrid,nxlocal,nzlocal,phi,rho,
     &                     dxsqi,dzsqi,xminodx,lrz,
     &                     mgparam,localbounds,mglevel,icndbndy)
        endif
      endif

c     --- Put desired potential onto conductors in phi array.
      call cond_potmg(conductors%interior,nxlocal,0,nzlocal,1,0,1,
     &                phi,mglevel,1,.false.)
      call condbndymgint(conductors%evensubgrid,nxlocal,0,nzlocal,1,0,1,
     &                   phi,localbounds,mglevel,icndbndy)
      call condbndymgint(conductors%oddsubgrid ,nxlocal,0,nzlocal,1,0,1,
     &                   phi,localbounds,mglevel,icndbndy)

      return
      end
c=============================================================================
      subroutine mgsor_loop2d(ixmin,ixmax,izmin,izmax,parity,mgparam,
     &                        nxlocal,nzlocal,phi,rho,
     &                        dxsqi,dzsqi,xminodx,lrz)
      integer(ISZ):: ixmin,ixmax,izmin,izmax,parity
      real(kind=8):: mgparam
      integer(ISZ):: nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,xminodx
      logical(ISZ):: lrz

      integer(ISZ):: ix,iz,ix1
      real(kind=8):: dxsqic0,dzsqic0,dxsqic,dzsqic
      real(kind=8):: const0,const,rr(0:nxlocal),rri(0:nxlocal)

c     --- Set temporary variables (these are used to increase performance)
      const = mgparam*0.5/(dxsqi + dzsqi)
      dxsqic = dxsqi*const
      dzsqic = dzsqi*const
      if (lrz) then
        const0 = mgparam*0.5/(2.*dxsqi + dzsqi)
        dxsqic0 = dxsqi*const0
        dzsqic0 = dzsqi*const0
      endif

      if (lrz) then

c       --- Precalculate the radius for efficiency
        do ix=0,nxlocal
          rr(ix) = ix + xminodx
          if (rr(ix) > 0.) rri(ix) = dxsqic/rr(ix)
        enddo

        do iz=izmin,izmax
          ix1 = ixmin + mod(ixmin + iz + parity,2)
          if (rr(ix1) == 0.) then
c           --- radius=0 is a special case
            phi(ix1,iz) = rho(ix1,iz)*const0 +
     &                   (4*phi(ix1+1,iz))*dxsqic0 +
     &                   (phi(ix1,iz-1) + phi(ix1,iz+1))*dzsqic0 +
     &                   (1. - mgparam)*phi(ix1,iz)
            ix1 = ix1 + 2
          endif
          do ix=ix1,ixmax,2
            phi(ix,iz) = rho(ix,iz)*const +
     &                   ((rr(ix)-0.5)*phi(ix-1,iz) + (rr(ix)+0.5)*phi(ix+1,iz))*rri(ix) +
     &                   (phi(ix,iz-1) + phi(ix,iz+1))*dzsqic +
     &                   (1. - mgparam)*phi(ix,iz)
          enddo
        enddo

      else

        do iz=izmin,izmax
          ix1 = ixmin + mod(ixmin + iz + parity,2)
          do ix=ix1,ixmax,2
            phi(ix,iz) = rho(ix,iz)*const +
     &                   (phi(ix-1,iz  ) + phi(ix+1,iz  ))*dxsqic +
     &                   (phi(ix  ,iz-1) + phi(ix  ,iz+1))*dzsqic +
     &                   (1. - mgparam)*phi(ix,iz)
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine condbndymg2d(subgrid,nxlocal,nzlocal,phi,rho,dxsqi,dzsqi,
     &                        xminodx,lrz,
     &                        mgparam,localbounds,mglevel,icndbndy)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nzlocal,mglevel
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1), rho(0:nxlocal,0:nzlocal)
      real(kind=8):: dxsqi,dzsqi,mgparam,xminodx
      logical(ISZ):: lrz
      integer(ISZ):: localbounds(0:5),icndbndy

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.

      integer(ISZ):: ic,ix,iz
      real(kind=8):: pxm,pzm,pxp,pzp
      real(kind=8):: dxm,dzm,dxp,dzp
      real(kind=8):: cxm,czm,cxp,czp
      real(kind=8):: voltfac,c0
      real(kind=8):: rr(0:nxlocal),rri(0:nxlocal)
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- Only use actual voltage on finest level. Set to zero for
c     --- coarser levels since solver for the residuals.
      if (mglevel == 0) then
        voltfac = 1.
      else
        voltfac = 0.
      endif

c     --- Note that in the loop below, for the Neumann boundary conditions,
c     --- the potentials are not set since they will be multiplied by zero
c     --- anyway. The code here just ensures that the variables have been
c     --- initialized so that debuggers and valgrind won't complain about
c     --- using uninitialized variables.
      pxm = 0.
      pxp = 0.
      pzm = 0.
      pzp = 0.

c     --- Precalculate the radius for efficiency
      if (lrz) then
        do ix=0,nxlocal
          rr(ix) = ix + xminodx
          if (rr(ix) > 0.) rri(ix) = 1./rr(ix)
        enddo
      endif

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet or parallel boundary
        if (ix == 0  .and. localbounds(0) < 1) cycle
        if (ix == nxlocal .and. localbounds(1) < 1) cycle
        if (iz == 0  .and. localbounds(4) < 1) cycle
        if (iz == nzlocal .and. localbounds(5) < 1) cycle

c       --- First, get the potential and effective grid cell sizes
c       --- Note that for the Neumann case, the potential is not
c       --- used and so is not set.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          pxm = voltfac*volt(0,ic)
          dxm = dels(0,ic)
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
        else
          pxm = phi(ix-1,iz  )
          dxm = 1.
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          pxp = voltfac*volt(1,ic)
          dxp = dels(1,ic)
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
        else
          pxp = phi(ix+1,iz  )
          dxp = 1.
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          pzm = voltfac*volt(4,ic)
          dzm = dels(4,ic)
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
        else
          pzm = phi(ix  ,iz-1)
          dzm = 1.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          dzp = dels(5,ic)
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
        else
          pzp = phi(ix  ,iz+1)
          dzp = 1.
        endif

c       --- Now construct the coefficients
        cxm = dxsqi/(dxm*(0.5*dxm + 0.5*dxp))
        cxp = dxsqi/(dxp*(0.5*dxm + 0.5*dxp))
        czm = dzsqi/(dzm*(0.5*dzm + 0.5*dzp))
        czp = dzsqi/(dzp*(0.5*dzm + 0.5*dzp))
        if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) cxm = 0.
        if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) cxp = 0.
        if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) czm = 0.
        if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) czp = 0.

c       --- Correct coefficients for axisymmetric case
        if (lrz) then
          if (rr(ix) > 0.) then
            cxm = cxm*(rr(ix) - 0.5*dxm)*rri(ix)
            cxp = cxp*(rr(ix) + 0.5*dxp)*rri(ix)
          else
            cxm = 0.
            cxp = 4.*cxp
          endif
        endif

c       --- Note that it is possible for c0 to be zero, but that is an
c       --- isolated case which doesn't make sense anyway, so let a NaN
c       --- happen rather than add an extra check on c0.

        c0 = cxm + cxp + czm + czp
        phi(ix,iz) = mgparam*(rho(ix,iz) + cxm*pxm + cxp*pxp + czm*pzm + czp*pzp)/c0
     &               + (1. - mgparam)*subgrid%prevphi(ic)

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine condbndymgres2d(subgrid,nxlocal,nzlocal,phi,rho,res,dxsqi,dzsqi,xminodx,lrz,
     &                           localbounds,mglevel,icndbndy,delx,delz)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nzlocal,mglevel
      integer(ISZ):: delx,delz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1), rho(0:nxlocal,0:nzlocal)
      real(kind=8):: res(-delx:nxlocal+delx,-delz:nzlocal+delz)
      real(kind=8):: dxsqi,dzsqi,xminodx
      logical(ISZ):: lrz
      integer(ISZ):: localbounds(0:5),icndbndy

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.

      integer(ISZ):: ic,ix,iz
      real(kind=8):: pxm,pzm,pxp,pzp
      real(kind=8):: dxm,dzm,dxp,dzp
      real(kind=8):: cxm,czm,cxp,czp
      real(kind=8):: voltfac,ppp,c0
      real(kind=8):: rr(0:nxlocal),rri(0:nxlocal)
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- Only use actual voltage on finest level. Set to zero for
c     --- coarser levels since solver for the residuals.
      if (mglevel == 0) then
        voltfac = 1.
      else
        voltfac = 0.
      endif

c     --- Note that in the loop below, for the Neumann boundary conditions,
c     --- the potentials are not set since they will be multiplied by zero
c     --- anyway. The code here just ensures that the variables have been
c     --- initialized so that debuggers and valgrind won't complain about
c     --- using uninitialized variables.
      pxm = 0.
      pxp = 0.
      pzm = 0.
      pzp = 0.

c     --- Precalculate the radius for efficiency
      if (lrz) then
        do ix=0,nxlocal
          rr(ix) = ix + xminodx
          if (rr(ix) > 0.) rri(ix) = 1./rr(ix)
        enddo
      endif

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet
        if (ix == 0  .and. localbounds(0) == 0) cycle
        if (ix == nxlocal .and. localbounds(1) == 0) cycle
        if (iz == 0  .and. localbounds(4) == 0) cycle
        if (iz == nzlocal .and. localbounds(5) == 0) cycle

c       --- First, get the potential and effective grid cell sizes
c       --- Note that for the Neumann case, the potential is not
c       --- used and so is not set.
        ppp = 1.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          pxm = voltfac*volt(0,ic)
          dxm = dels(0,ic)
          ppp = min(ppp,dels(0,ic))
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
          if (abs(dels(0,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(0,ic)))
          endif
        else
          pxm = phi(ix-1,iz  )
          dxm = 1.
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          pxp = voltfac*volt(1,ic)
          dxp = dels(1,ic)
          ppp = min(ppp,dels(1,ic))
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
          if (abs(dels(1,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(1,ic)))
          endif
        else
          pxp = phi(ix+1,iz  )
          dxp = 1.
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          pzm = voltfac*volt(4,ic)
          dzm = dels(4,ic)
          ppp = min(ppp,dels(4,ic))
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
          if (abs(dels(4,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(4,ic)))
          endif
        else
          pzm = phi(ix  ,iz-1)
          dzm = 1.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          dzp = dels(5,ic)
          ppp = min(ppp,dels(5,ic))
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
          if (abs(dels(5,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(5,ic)))
          endif
        else
          pzp = phi(ix  ,iz+1)
          dzp = 1.
        endif

c       --- Now construct the coefficients
        cxm = dxsqi/(dxm*(0.5*dxm + 0.5*dxp))
        cxp = dxsqi/(dxp*(0.5*dxm + 0.5*dxp))
        czm = dzsqi/(dzm*(0.5*dzm + 0.5*dzp))
        czp = dzsqi/(dzp*(0.5*dzm + 0.5*dzp))
        if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) cxm = 0.
        if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) cxp = 0.
        if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) czm = 0.
        if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) czp = 0.

c       --- Correct coefficients for axisymmetric case
        if (lrz) then
          if (rr(ix) > 0.) then
            cxm = cxm*(rr(ix) - 0.5*dxm)*rri(ix)
            cxp = cxp*(rr(ix) + 0.5*dxp)*rri(ix)
          else
            cxm = 0.
            cxp = 4.*cxp
          endif
        endif

        c0 = cxm + cxp + czm + czp
        res(ix,iz) = ppp*(rho(ix,iz) + cxm*pxm + cxp*pxp + czm*pzm + czp*pzp - c0*phi(ix,iz))

        if (0. >= dels(0,ic) .and. dels(0,ic) >= -1.) res(ix-1,iz) = 0.
        if (0. >= dels(1,ic) .and. dels(1,ic) >= -1.) res(ix+1,iz) = 0.
        if (0. >= dels(4,ic) .and. dels(4,ic) >= -1.) res(ix,iz-1) = 0.
        if (0. >= dels(5,ic) .and. dels(5,ic) >= -1.) res(ix,iz+1) = 0.

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine residual2d(nxlocal,nzlocal,delx,delz,dxsqi,dzsqi,
     &                      xminodx,lrz,phi,rho,res,mglevel,localbounds,
     &                      lcndbndy,icndbndy,conductors)
      use ConductorTypemodule
      integer(ISZ):: nxlocal,nzlocal
      integer(ISZ):: delx,delz
      real(kind=8):: dxsqi,dzsqi,xminodx
      logical(ISZ):: lrz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: res(-delx:nxlocal+delx,-delz:nzlocal+delz)
      integer(ISZ):: mglevel,localbounds(0:5)
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.
c
c For internal conductors, the residual is set to zero inside and calculated
c using the modified form of the finite differenced Poisson's equation near
c the surface.

      integer(ISZ):: ix,iz,ix1
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: const,const0
      real(kind=8):: rr(0:nxlocal),rri(0:nxlocal)

      const = 2.*(dxsqi+dzsqi)
      const0 = 2.*(2.*dxsqi + dzsqi)

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nxlocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) == 0) ixmin = 1
      if (localbounds(1) == 0) ixmax = nxlocal-1
      if (localbounds(4) == 0) izmin = 1
      if (localbounds(5) == 0) izmax = nzlocal-1

c     --- Precalculate the radius for efficiency
      do ix=0,nxlocal
        rr(ix) = ix + xminodx
        if (rr(ix) > 0.) rri(ix) = dxsqi/rr(ix)
      enddo

c     --- This zeroing out is done inside the loop so that the cache look
c     --- ups for the zeroing and the calculation are done at the same time.
      if (localbounds(4) == 0) then
c       --- At Dirichlet boundaries, zero out res at the boundary.
        res(:,-delz:0) = 0.
      else if (localbounds(4)  > 0 .and. delz > 1) then
c       --- Otherwise, only zero out the guard cells that won't be otherwise
c       --- filled in.
        res(:,-delz:-2) = 0
      endif

c     --- Calculate the residual.
      if (lrz) then

!$OMP DO
        do iz=izmin,izmax

          if (localbounds(0) == 0) then
            res(-delx:0,iz) = 0.
          else if (localbounds(0)  > 0 .and. delx > 1) then
            res(-delx:-2,iz) = 0.
          endif

          ix1 = ixmin
          if (rr(ix1) == 0) then
            res(ix1,iz) = rho(ix1,iz)
     &                    + (4*phi(ix1+1,iz  ))*dxsqi
     &                    + (phi(ix1,iz-1) + phi(ix1,iz+1))*dzsqi
     &                    - phi(ix1,iz)*const0
            ix1 = ix1 + 1
          endif
          do ix=ix1,ixmax
            res(ix,iz) = rho(ix,iz)
     &                   + ((rr(ix)-0.5)*phi(ix-1,iz  ) + (rr(ix)+0.5)*phi(ix+1,iz  ))*rri(ix)
     &                   + (phi(ix,iz-1) + phi(ix,iz+1))*dzsqi
     &                   - phi(ix,iz)*const
          enddo

          if (localbounds(1) == 0) then
            res(nxlocal:nxlocal+delx,iz) = 0.
          else if (localbounds(1)  > 0 .and. delx > 1) then
            res(nxlocal+2:nxlocal+delx,iz) = 0.
          endif

        enddo
!$OMP END DO

      else

!$OMP DO
        do iz=izmin,izmax

          if (localbounds(0) == 0) then
            res(-delx:0,iz) = 0.
          else if (localbounds(0)  > 0 .and. delx > 1) then
            res(-delx:-2,iz) = 0.
          endif

          do ix=ixmin,ixmax

            res(ix,iz) = rho(ix,iz)
     &          +  (phi(ix-1,iz  )+phi(ix+1,iz  ))*dxsqi
     &          +  (phi(ix  ,iz-1)+phi(ix  ,iz+1))*dzsqi
     &          -  phi(ix,iz)*const

          enddo

          if (localbounds(1) == 0) then
            res(nxlocal:nxlocal+delx,iz) = 0.
          else if (localbounds(1)  > 0 .and. delx > 1) then
            res(nxlocal+2:nxlocal+delx,iz) = 0.
          endif

        enddo
!$OMP END DO
      endif

      if (localbounds(5) == 0) then
        res(:,nzlocal:nzlocal+delz) = 0.
      else if (localbounds(5)  > 0 .and. delz > 1) then
        res(:,nzlocal+2:nzlocal+delz) = 0.
      endif

c     --- Zero the residual inside conductors.
      call cond_potmgres(conductors%interior,nxlocal,0,nzlocal,
     *                   delx,0,delz,res,mglevel,1,.false.)

      if (lcndbndy) then
c       --- Calculate the residual near the conductor.
        call condbndymgres2d(conductors%evensubgrid,nxlocal,nzlocal,phi,rho,res,
     &                       dxsqi,dzsqi,xminodx,lrz,localbounds,
     &                       mglevel,icndbndy,delx,delz)
        call condbndymgres2d(conductors%oddsubgrid,nxlocal,nzlocal,phi,rho,res,
     &                       dxsqi,dzsqi,xminodx,lrz,localbounds,
     &                       mglevel,icndbndy,delx,delz)
      endif

      call applyboundaryconditions3d(nxlocal,0,nzlocal,
     &                               delx,0,delz,
     &                               res,1,localbounds,.false.,.false.)
c     --- There is a problem with conductors and periodic boundaries. If a
c     --- conductor reaches the edge of the mesh, the residual at
c     --- the corresponding points on the opposite boundary will not be
c     --- zeroed out as it should be.

      return
      end
c=============================================================================
      subroutine restrict2d(nx,nz,nxlocal,nzlocal,u,delx,delz,
     &                      nxcoarse,nzcoarse,
     &                      nxlocalcoarse,nzlocalcoarse,ucoarse,
     &                      delcx,delcz,
     &                      localbounds,localboundscoarse,lxoffset,lzoffset,
     &                      lapplydirichlet)
      integer(ISZ):: nx,nz,nxlocal,nzlocal,delx,delz,delcx,delcz
      integer(ISZ):: nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse
      real(kind=8):: u(-delx:nxlocal+delx,-delz:nzlocal+delz)
      real(kind=8):: ucoarse(-delcx:nxlocalcoarse+delcx,
     &                       -delcz:nzlocalcoarse+delcz)
      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: localboundscoarse(0:5)
      integer(ISZ):: lxoffset,lzoffset
      logical(ISZ):: lapplydirichlet
      
c Restrict to a coarser grid.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixcoarse,izcoarse
      integer(ISZ):: ixcoarsemin,ixcoarsemax,izcoarsemin,izcoarsemax
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ),allocatable:: ixmina(:),ixmaxa(:)
      integer(ISZ),allocatable:: izmina(:),izmaxa(:)
      real(kind=8),allocatable:: wxa(:,:),wza(:,:)
      real(kind=8):: r,w,dx,dz,dxi,dzi,wx(0:3),wz(0:3)

      allocate(ixmina(0:nxlocalcoarse),ixmaxa(0:nxlocalcoarse))
      allocate(izmina(0:nzlocalcoarse),izmaxa(0:nzlocalcoarse))
      allocate(wxa(0:3,0:nxlocalcoarse))
      allocate(wza(0:3,0:nzlocalcoarse))
      wxa = 0.
      wza = 0.

c     --- Set the loop limits, including edges when appropriate.
      ixcoarsemin = 0
      ixcoarsemax = nxlocalcoarse
      izcoarsemin = 0
      izcoarsemax = nzlocalcoarse
      if (lapplydirichlet) then
        if (localboundscoarse(0) == 0) ixcoarsemin = 1
        if (localboundscoarse(1) == 0) ixcoarsemax = nxlocalcoarse - 1
        if (localboundscoarse(4) == 0) izcoarsemin = 1
        if (localboundscoarse(5) == 0) izcoarsemax = nzlocalcoarse - 1

        if (localboundscoarse(0) == 0) ucoarse(0,:) = 0.
        if (localboundscoarse(1) == 0) ucoarse(nxlocalcoarse,:) = 0.
        if (localboundscoarse(4) == 0) ucoarse(:,0) = 0.
        if (localboundscoarse(5) == 0) ucoarse(:,nzlocalcoarse) = 0.
      endif

      dx = 1.*nx/nxcoarse
      dz = 1.*nz/nzcoarse
      dxi = 1.*nxcoarse/nx
      dzi = 1.*nzcoarse/nz

c     --- Precalculate the loop limits and weights. This saves a surprisingly
c     --- substantial amount of time.
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = ((izcoarse-1)*nz - lzoffset + 4*nzcoarse)/nzcoarse-3
        izmax = ((izcoarse+1)*nz - lzoffset - 1)/nzcoarse
        if (izmin < -delz) izmin = -delz
        if (izmax > nzlocal+delz) izmax = nzlocal+delz
        if (izmax < izmin) continue

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzcoarse)*dzi)
        enddo

        if (izcoarse == 0 .and. localboundscoarse(4) == 0) then
          izmin = 0
          izmax = 0
          wz(0) = 2.
        else if (izcoarse == nzlocalcoarse .and. localboundscoarse(5) == 0) then
          izmin = nzlocal
          izmax = nzlocal
          wz(0) = 2.
        endif

        izmina(izcoarse) = izmin
        izmaxa(izcoarse) = izmax
        wza(:,izcoarse) = wz

      enddo

      do ixcoarse=ixcoarsemin,ixcoarsemax
        ixmin = ((ixcoarse-1)*nx - lxoffset + 4*nxcoarse)/nxcoarse-3
        ixmax = ((ixcoarse+1)*nx - lxoffset - 1)/nxcoarse
        if (ixmin < -delx) ixmin = -delx
        if (ixmax > nxlocal+delx) ixmax = nxlocal+delx
        if (ixmax < ixmin) continue

        do ix=ixmin,ixmax
          wx(ix-ixmin) = 1. - abs(ixcoarse - (ix + 1.*lxoffset/nxcoarse)*dxi)
        enddo

        if (ixcoarse == 0 .and. localboundscoarse(0) == 0) then
          ixmin = 0
          ixmax = 0
          wx(0) = 2.
        else if (ixcoarse == nxcoarse .and. localboundscoarse(1) == 0) then
          ixmin = nx
          ixmax = nx
          wx(0) = 2.
        endif

        ixmina(ixcoarse) = ixmin
        ixmaxa(ixcoarse) = ixmax
        wxa(:,ixcoarse) = wx

      enddo

c     --- Do the loops.
!$OMP DO
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = izmina(izcoarse)
        izmax = izmaxa(izcoarse)
        wz(:) = wza(:,izcoarse)

        do ixcoarse=ixcoarsemin,ixcoarsemax
          ixmin = ixmina(ixcoarse)
          ixmax = ixmaxa(ixcoarse)
          wx(:) = wxa(:,ixcoarse)

          r = 0.
          w = 0.
          do iz=izmin,izmax
              do ix=ixmin,ixmax
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

      deallocate(ixmina,ixmaxa)
      deallocate(izmina,izmaxa)
      deallocate(wxa,wza)

      return
      end
c=============================================================================
      subroutine expand2d(nx,nz,nxlocal,nzlocal,phi,
     &                    nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                    phicoarse,bounds,lxoffset,lzoffset)
      integer(ISZ):: nx,nz,nxlocal,nzlocal
      integer(ISZ):: nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1)
      real(kind=8):: phicoarse(-1:nxlocalcoarse+1,-1:nzlocalcoarse+1)
      integer(ISZ):: lxoffset,lzoffset
      integer(ISZ):: bounds(0:5)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ix,iz
      integer(ISZ):: jx,jz
      integer(ISZ),pointer:: jxa(:),jza(:)
      real(kind=8),pointer:: wxa(:),wza(:)
      real(kind=8):: dx,dz
      real(kind=8):: wx,wz

      allocate(jxa(0:nxlocal),jza(0:nzlocal))
      allocate(wxa(0:nxlocal),wza(0:nzlocal))

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nxlocal
      izmin = 0
      izmax = nzlocal
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nxlocal - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nzlocal - 1

      dx = 1.*nxcoarse/nx
      dz = 1.*nzcoarse/nz

c     --- Precalculate the indices and weights. This saves a surprisingly
c     --- substantial amount of time.
      do iz=izmin,izmax
        jza(iz) = int((iz*nzcoarse + lzoffset)/nz)
        wza(iz) =  1.*(iz*nzcoarse + lzoffset)/nz - jza(iz)
      enddo
      do ix=ixmin,ixmax
        jxa(ix) = int((ix*nxcoarse + lxoffset)/nx)
        wxa(ix) =  1.*(ix*nxcoarse + lxoffset)/nx - jxa(ix)
      enddo

!$OMP DO
      do iz=izmin,izmax
        jz = jza(iz)
        wz = wza(iz)
        do ix=ixmin,ixmax
          jx = jxa(ix)
          wx = wxa(ix)

          phi(ix,iz) = phi(ix,iz) +
     &           (1.-wx)*(1.-wz)*phicoarse(jx  ,jz  ) +
     &               wx *(1.-wz)*phicoarse(jx+1,jz  ) +
     &           (1.-wx)*    wz *phicoarse(jx  ,jz+1) +
     &               wx *    wz *phicoarse(jx+1,jz+1)
        enddo
      enddo
!$OMP ENDDO

      deallocate(jxa,jza)
      deallocate(wxa,wza)

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigrid2ddielectricsolve(iwhich,nx,nz,nxlocal,nzlocal,dx,dz,
     &                            phi,rho,epsilon,bounds,xmminlocal,
     &                            mgparam,mgiters,mgmaxiters,
     &                            mgmaxlevels,mgerror,mgtol,mgverbose,
     &                            downpasses,uppasses,
     &                            lcndbndy,laddconductor,
     &                            gridmode,conductors,lrz,
     &                            fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      use Decompositionmodule
      integer(ISZ):: iwhich
      integer(ISZ):: nx,nz,nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nzlocal)
      real(kind=8):: epsilon(0:nxlocal+1,0:nzlocal+1)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmminlocal
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels,mgverbose
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor
      integer(ISZ):: gridmode
      type(ConductorType):: conductors
      logical(ISZ):: lrz
      type(Decomposition):: fsdecomp

c Use the multigrid method for solving Poisson's equation on a 2D Cartesian
c mesh. The fieldsolver allows for a variable dielectric constant.
c
c When the grid cells are rectangular, semi-coarsening is done until the
c grid cell dimensions are roughly equal. Roughly equal means that
c   2/3 dx < dz < 4/3 dx
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3.

      integer(ISZ):: i,k,ix,iz
      real(kind=8),allocatable:: phisave(:,:)

      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: lxoffset(0:fsdecomp%nxprocs-1)
      integer(ISZ):: rxoffset(0:fsdecomp%nxprocs-1)
      integer(ISZ):: lzoffset(0:fsdecomp%nzprocs-1)
      integer(ISZ):: rzoffset(0:fsdecomp%nzprocs-1)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Note that nx and nzlocal do not need to be even

c     --- Only slab geometry is now supported.
      if (lrz) then
        call kaboom("multigrid2ddielectricsolve: only slab geometry supported")
        return
      endif

c     --- The parallel version does not yet work, due to the issue dealing
c     --- with epsilon, described below.
      if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
        print*,"multigrid2ddielectricsolve: does not yet work in parallel"
        call kaboom("multigrid2ddielectricsolve: does not yet work in parallel")
        return
      endif

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

      localbounds = bounds
#ifdef MPIPARALLEL
      if (fsdecomp%ix(fsdecomp%ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(fsdecomp%ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iz(fsdecomp%izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(fsdecomp%izproc)+nzlocal < nz) localbounds(5) = -1
#endif

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
      call checkconductors(nx,0,nz,nxlocal,0,nzlocal,dx,dx,dz,
     &                     conductors,fsdecomp)

!$OMP PARALLEL
!$OMP&PRIVATE(i,ix,iz)

#ifdef MPIPARALLEL
c     --- These calls break the parallel field solver
c     call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
c    &                             1,0,1,-1,0,localbounds,fsdecomp)
c     call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif

c     --- Make sure guard planes have sensible values before beginning.
      call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,localbounds,
     &                               .false.,.false.)

      allocate(phisave(-1:nxlocal+1,-1:nzlocal+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid2ddielectricsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nxlocal,nzlocal
        call kaboom("multigrid2ddielectricsolve: allocation error")
        return
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
        call vcycle2ddielectric(0,1.,nx,nz,nxlocal,nzlocal,dx,dz,
     &                phi,rho,epsilon,
     &                bounds,mgparam,mgmaxlevels,
     &                downpasses,uppasses,lcndbndy,conductors,
     &                xmminlocal,lrz,fsdecomp)

c       --- Calculate the change in phi.
        mgerror = 0.
!$OMP DO REDUCTION(MAX:mgerror)
        do iz=0,nzlocal
          do ix=0,nxlocal
            mgerror = max(mgerror,abs(phisave(ix,iz) - phi(ix,iz)))
          enddo
        enddo
!$OMP END DO

#ifdef MPIPARALLEL
        if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
c         --- calculate global sorerror
          call parallelmaxrealarray(mgerror,1)
        endif
#endif

      enddo

c     --- Set boundary conditions. This is only really needed for the
c     --- Dirichlet boundaries, but this is convenient to call.
      call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,
     &                               localbounds,.true.,.false.)

      if (mgverbose>=1 .or. (mgverbose>=0 .and. mgerror > mgtol)) then
c       --- Make a print out.
        if (mgerror > mgtol) then
          call remark("Multigrid2ddielectric: Maximum number of iterations reached")
        endif
        write(errline,20) mgerror,mgiters
  20    format("Multigrid2ddielectric: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif

      deallocate(phisave)

!$OMP END PARALLEL

      if (lf3dtimesubs) timemultigrid2dsolve = timemultigrid2dsolve +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycle2ddielectric(mglevel,mgscale,nx,nz,
     &                              nxlocal,nzlocal,dx,dz,
     &                              phi,rho,epsilon,globalbounds,mgparam,
     &                              mgmaxlevels,downpasses,uppasses,
     &                              lcndbndy,conductors,xmminlocal,lrz,
     &                              fsdecomp)
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use formggetarraysuminterface
      use Constant,Only: eps0
      use Decompositionmodule
      integer(ISZ):: mglevel
      real(kind=8):: mgscale
      integer(ISZ):: nx,nz,nxlocal,nzlocal
      real(kind=8):: dx,dz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: epsilon(0:nxlocal+1,0:nzlocal+1)
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      real(kind=8):: xmminlocal
      logical(ISZ):: lcndbndy,lrz
      type(Decomposition):: fsdecomp

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: epsiloncoarse(:,:)
      real(kind=8),allocatable:: res(:,:)
      integer(ISZ):: i,iszone=1
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nylocal
      integer(ISZ):: nxlocalcoarse,nzlocalcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dzcoarsesqi
      real(kind=8):: xminodx,xmminlocalcoarse
      real(kind=8):: mgscalecoarse
      integer(ISZ):: ixproc,izproc
      integer(ISZ):: localbounds(0:5),localboundsc(0:5)
      integer(ISZ):: lxoffsetall(0:fsdecomp%nxprocs-1)
      integer(ISZ):: rxoffsetall(0:fsdecomp%nxprocs-1)
      integer(ISZ):: lzoffsetall(0:fsdecomp%nzprocs-1)
      integer(ISZ):: rzoffsetall(0:fsdecomp%nzprocs-1)
      integer(ISZ):: lxoffset,rxoffset
      integer(ISZ):: lzoffset,rzoffset
      type(Decomposition):: coarsedecomp
      integer(ISZ):: allocerror
      integer(ISZ):: delx,delz
      logical(ISZ):: lpe0
      real(kind=8):: sss(2)

      localbounds = globalbounds
      xminodx = xmminlocal/dx

#ifdef MPIPARALLEL
      ixproc = fsdecomp%ixproc
      izproc = fsdecomp%izproc
      if (fsdecomp%ix(ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iz(izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(izproc)+nzlocal < nz) localbounds(5) = -1
#endif

      if (lprintmgarraysumdiagnostic) then
#ifdef MPIPARALLEL
        lpe0=(fsdecomp%ixproc==0.and.fsdecomp%izproc==0)
#else
        lpe0 = .true.
#endif
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V1 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V1 rho",mglevel,sss/eps0
      endif

c     --- Do initial SOR passes.
      do i=1,downpasses
        call relax2ddielectric(mglevel,nx,nz,nxlocal,nzlocal,phi,rho,epsilon,
     &                         dx,dz,xminodx,localbounds,mgparam,lcndbndy,
     &                         conductors,lrz,fsdecomp)
      enddo

      if (lprintmgarraysumdiagnostic) then
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V2 rho",mglevel,sss/eps0
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. nz >= 4 .and.
     &    mglevel < mgmaxlevels) then

        delx = 1
        delz = 1
        if (nx > nxlocal) delx = 3
        if (nz > nzlocal) delz = 3

        allocate(res(-delx:nxlocal+delx,-delz:nzlocal+delz),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2ddielectric: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nxlocal,nzlocal
          call kaboom("vcycle2ddielectric: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residual2ddielectric(nxlocal,nzlocal,delx,delz,
     &                            phi,rho,epsilon,res,dx,dz,
     &                            mglevel,localbounds,lcndbndy,conductors,
     &                            xminodx,lrz)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,res,
     &                               delx,0,delz,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,0,nzlocal,res,delx,0,delz,-3,-1,fsdecomp)
#endif
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocal,0,nzlocal,delx,0,delz,res,fsdecomp)
          if (lpe0) print*,"V3 res",mglevel,sss/eps0
        endif

c       --- Note that some y quantities are included as dummies since the
c       --- routine will change the values.
        call getnextcoarselevel3d(nx,0,nz,nxlocal,nylocal,nzlocal,dx,dx,dz,
     &                            nxcoarse,nycoarse,nzcoarse,
     &                            dxcoarse,dycoarse,dzcoarse)

        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2
c       --- This option is not supported
c       mgscalecoarse = mgscale*dxcoarse*dzcoarse/(dx*dz)
        mgscalecoarse = 1.

        localboundsc = globalbounds

#ifdef MPIPARALLEL
        coarsedecomp%nxglobal = nxcoarse
        coarsedecomp%nyglobal = 0
        coarsedecomp%nzglobal = nzcoarse
        coarsedecomp%mpi_comm_x = fsdecomp%mpi_comm_x
        coarsedecomp%mpi_comm_z = fsdecomp%mpi_comm_z
        coarsedecomp%ixproc = fsdecomp%ixproc
        coarsedecomp%iyproc = fsdecomp%iyproc
        coarsedecomp%izproc = fsdecomp%izproc
        coarsedecomp%nxprocs = fsdecomp%nxprocs
        coarsedecomp%nyprocs = fsdecomp%nyprocs
        coarsedecomp%nzprocs = fsdecomp%nzprocs
        allocate(coarsedecomp%ix(0:fsdecomp%nxprocs-1))
        allocate(coarsedecomp%nx(0:fsdecomp%nxprocs-1))
        allocate(coarsedecomp%iy(0:fsdecomp%nyprocs-1))
        allocate(coarsedecomp%ny(0:fsdecomp%nyprocs-1))
        allocate(coarsedecomp%iz(0:fsdecomp%nzprocs-1))
        allocate(coarsedecomp%nz(0:fsdecomp%nzprocs-1))
c       --- Find domains in coarser grid
        call mgdividenz(fsdecomp,coarsedecomp,nx,0,nz,
     &                  nxcoarse,0,nzcoarse,mgscale)
c       --- Reset value to corrected one
        nxlocalcoarse = coarsedecomp%nx(ixproc)
        nzlocalcoarse = coarsedecomp%nz(izproc)
c       --- Difference between starts and ends of coarse and fine grids.
c       --- Should only be in the range 0-2.
        lxoffsetall = (nxcoarse*fsdecomp%ix-nx*coarsedecomp%ix)
        rxoffsetall = (nx*(coarsedecomp%ix + coarsedecomp%nx) -
     &                 nxcoarse*(fsdecomp%ix + fsdecomp%nx))
        lzoffsetall = (nzcoarse*fsdecomp%iz-nz*coarsedecomp%iz)
        rzoffsetall = (nz*(coarsedecomp%iz + coarsedecomp%nz) -
     &                 nzcoarse*(fsdecomp%iz + fsdecomp%nz))
c       --- Note that the lzoffsetall and rzoffsetall can only be used in
c       --- MPIPARALLEL sections since they will be unallocated in the
c       --- serial code. So, separate scalars are used in code which is
c       --- used in the serial version.
        lxoffset = lxoffsetall(ixproc)
        rxoffset = rxoffsetall(ixproc)
        lzoffset = lzoffsetall(izproc)
        rzoffset = rzoffsetall(izproc)
        if (coarsedecomp%ix(ixproc) > 0) localboundsc(0) = -1
        if (coarsedecomp%ix(ixproc)+nxlocalcoarse < nxcoarse) localboundsc(1) = -1
        if (coarsedecomp%iz(izproc) > 0) localboundsc(4) = -1
        if (coarsedecomp%iz(izproc)+nzlocalcoarse < nzcoarse) localboundsc(5) = -1
c       --- Calculate the xmminlocal of the coarse grid
        xmminlocalcoarse = xmminlocal
     &                          - fsdecomp%ix(fsdecomp%ixproc)*dx
     &                          + coarsedecomp%ix(fsdecomp%ixproc)*dxcoarse
#else
        nxlocalcoarse = nxcoarse
        nzlocalcoarse = nzcoarse
        lxoffset = 0
        rxoffset = 0
        lzoffset = 0
        rzoffset = 0
        xmminlocalcoarse = xmminlocal
#endif

c       --- Alloate new work space
        allocate(phicoarse(-1:nxlocalcoarse+1,-1:nzlocalcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2ddielectric: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcycle2ddielectric: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxlocalcoarse,0:nzlocalcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2ddielectric: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcycle2ddielectric: allocation error")
          return
        endif
        allocate(epsiloncoarse(0:nxlocalcoarse+1,0:nzlocalcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycle2ddielectric: allocation error ",allocerror,
     &           ": could not allocate epsiloncoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcycle2ddielectric: allocation error")
          return
        endif

        rhocoarse = 0.
        phicoarse = 0.

c       --- Restriction - note that no scaling factor is needed
        call restrict2d(nx,nz,nxlocal,nzlocal,res,delx,delz,
     &                  nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                  rhocoarse,0,0,
     &                  localbounds,localboundsc,lxoffset,lzoffset,.true.)
        call restrict2dcellcentered(nx,nz,nxlocal,nzlocal,epsilon,
     &                    nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                    epsiloncoarse,
     &                    localbounds,localboundsc,lxoffset,lzoffset)
        call applyboundaryconditions3d(nxlocalcoarse-1,0,nzlocalcoarse-1,
     &                                 1,0,1,epsiloncoarse,1,localbounds,
     &                                 .true.,.false.)
c       --- XXX epsiloncoarse needs to be exchanged among neighboring
c       --- XXX processors. This code needs to be written, so this solver
c       --- XXX will not work in parallel. (This is the only issue -
c       --- XXX everything else here works OK in parallel.)

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,0,0,0,rhocoarse,coarsedecomp)
          if (lpe0) print*,"V3 rhocoarse",mglevel,sss/eps0
        endif

c       --- Continue at the next coarsest level.
        call vcycle2ddielectric(mglevel+iszone,mgscalecoarse,nxcoarse,nzcoarse,
     &                nxlocalcoarse,nzlocalcoarse,
     &                dxcoarse,dzcoarse,phicoarse,rhocoarse,epsiloncoarse,
     &                globalbounds,mgparam,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,conductors,xmminlocalcoarse,lrz,
     &                coarsedecomp)

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,1,0,1,phicoarse,coarsedecomp)
          if (lpe0) print*,"V4 phicoarse",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand2d(nx,nz,nxlocal,nzlocal,phi,
     &                nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,phicoarse,
     &                localbounds,lxoffset,lzoffset)
        call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,localbounds,
     &                                 .false.,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,-1,-1,localbounds,fsdecomp)
#endif

        deallocate(phicoarse,rhocoarse,epsiloncoarse)
        deallocate(res)

#ifdef MPIPARALLEL
        deallocate(coarsedecomp%ix)
        deallocate(coarsedecomp%nx)
        deallocate(coarsedecomp%iy)
        deallocate(coarsedecomp%ny)
        deallocate(coarsedecomp%iz)
        deallocate(coarsedecomp%nz)
#endif

      endif

      if (lprintmgarraysumdiagnostic) then
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,0,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V5 phi",mglevel,sss
      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call relax2ddielectric(mglevel,nx,nz,nxlocal,nzlocal,phi,rho,epsilon,
     &                         dx,dz,xminodx,localbounds,mgparam,lcndbndy,
     &                         conductors,lrz,fsdecomp)
      enddo

      return
      end
c=============================================================================
      subroutine relax2ddielectric(mglevel,nx,nz,nxlocal,nzlocal,phi,rho,
     &                             epsilon,dx,dz,xminodx,
     &                             localbounds,mgparam,lcndbndy,conductors,lrz,
     &                             fsdecomp)
      use Constant
      use ConductorTypemodule
      use Decompositionmodule
      integer(ISZ):: mglevel,nx,nz,nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: epsilon(0:nxlocal+1,0:nzlocal+1)
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: dx,dz,xminodx,mgparam
      logical(ISZ):: lcndbndy,lrz
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.
c phisave is needed because when there are B fields, the even and odd sweeps
c are no longer independent since the diagonal coefficients become non-zero.
c The simplest solution is to save the values of phi before a sweep and
c read those saved values for the update of phi.

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ixmin,ixmax,izmin,izmax,ix,iz,ix1
      integer(ISZ):: i1,i2,ic
      real(kind=8):: dxsqi,dzsqi

      dxsqi = 1./dx**2
      dzsqi = 1./dz**2

      call cond_potmg(conductors%interior,
     &                nxlocal,0,nzlocal,1,0,1,phi,mglevel,1,.false.)
      call condbndymgint(conductors%evensubgrid,nxlocal,0,nzlocal,1,0,1,phi,
     &                   localbounds,mglevel,2)
      call condbndymgint(conductors%oddsubgrid,nxlocal,0,nzlocal,1,0,1,phi,
     &                   localbounds,mglevel,2)

c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      parity = + fsdecomp%ix(fsdecomp%ixproc)
     &         + fsdecomp%iz(fsdecomp%izproc)
      s_parity = mod(parity,2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

c     --- Set min and max indices for phi array.
c     --- If using Dirichlet boundary conditions, do not solve for the
c     --- potential on the grid edge.
      ixmin = 0
      ixmax = nxlocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) < 1) ixmin = 1
      if (localbounds(1) < 1) ixmax = nxlocal - 1
      if (localbounds(4) < 1) izmin = 1
      if (localbounds(5) < 1) izmax = nzlocal - 1

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

        do iz=izmin,izmax
          ix1 = ixmin + mod(ixmin + iz + parity,2)
          do ix=ix1,ixmax,2

            phi(ix,iz) = mgparam*(rho(ix,iz) + 
     &        0.5*(epsilon(ix  ,iz)+epsilon(ix  ,iz+1))*phi(ix-1,iz  )*dxsqi +
     &        0.5*(epsilon(ix+1,iz)+epsilon(ix+1,iz+1))*phi(ix+1,iz  )*dxsqi +
     &        0.5*(epsilon(ix,iz  )+epsilon(ix+1,iz  ))*phi(ix  ,iz-1)*dzsqi +
     &        0.5*(epsilon(ix,iz+1)+epsilon(ix+1,iz+1))*phi(ix  ,iz+1)*dzsqi)/
     &                    (0.5*(epsilon(ix  ,iz)+epsilon(ix  ,iz+1) +
     &                          epsilon(ix+1,iz)+epsilon(ix+1,iz+1))*
     &                     (dxsqi + dzsqi)) +
     &                   (1.-mgparam)*phi(ix,iz)

          enddo
        enddo

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
          if (parity == 0) then
           call condbndymg2ddielectric(conductors%evensubgrid,nxlocal,nzlocal,
     &                                 phi,rho,epsilon,
     &                                 dxsqi,dzsqi,xminodx,lrz,
     &                                 mgparam,localbounds,mglevel)
          endif
          if (parity == 1) then
           call condbndymg2ddielectric(conductors%oddsubgrid,nxlocal,nzlocal,
     &                                 phi,rho,epsilon,
     &                                 dxsqi,dzsqi,xminodx,lrz,
     &                                 mgparam,localbounds,mglevel)
          endif
        endif

        call cond_potmg(conductors%interior,
     &                  nxlocal,0,nzlocal,1,0,1,phi,mglevel,1,.false.)
        call condbndymgint(conductors%evensubgrid,nxlocal,0,nzlocal,1,0,1,phi,
     &                     localbounds,mglevel,2)
        call condbndymgint(conductors%oddsubgrid ,nxlocal,0,nzlocal,1,0,1,phi,
     &                     localbounds,mglevel,2)

        call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,localbounds,
     &                                 .false.,.false.)

#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,0,fsdecomp)
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- Exchange phi in the z guard planes
      call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif

      return
      end
c=============================================================================
      subroutine residual2ddielectric(nxlocal,nzlocal,delx,delz,
     &                                phi,rho,epsilon,res,dx,dz,
     &                                mglevel,localbounds,lcndbndy,conductors,
     &                                xminodx,lrz)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nxlocal,nzlocal,delx,delz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: epsilon(0:nxlocal+1,0:nzlocal+1)
      real(kind=8):: res(-delx:nxlocal+delx,-delz:nzlocal+delz)
      real(kind=8):: dx,dz
      integer(ISZ):: mglevel,localbounds(0:5)
      logical(ISZ):: lcndbndy
      type(ConductorType):: conductors
      real(kind=8):: xminodx
      logical(ISZ):: lrz

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: dxsqi,dzsqi

      dxsqi = 1./dx**2
      dzsqi = 1./dz**2

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nxlocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) == 0) ixmin = 1
      if (localbounds(1) == 0) ixmax = nxlocal - 1
      if (localbounds(4) == 0) izmin = 1
      if (localbounds(5) == 0) izmax = nzlocal - 1

c     --- Calculate the residual.

c     --- This zeroing out is done inside the loop so that the cache look
c     --- ups for the zeroing and the calculation are done at the same time.
      if (localbounds(4) == 0) then
c       --- At Dirichlet boundaries, zero out res at the boundary.
        res(:,-delz:0) = 0.
      else if (localbounds(4)  > 0 .and. delz > 1) then
c       --- Otherwise, only zero out the guard cells that won't be otherwise
c       --- filled in.
        res(:,-delz:-2) = 0
      endif

!$OMP DO
      do iz=izmin,izmax

        if (localbounds(0) == 0) then
          res(-delx:0,iz) = 0.
        else if (localbounds(0)  > 0 .and. delx > 1) then
          res(-delx:-2,iz) = 0.
        endif

        do ix=ixmin,ixmax

          res(ix,iz) = rho(ix,iz)
     &      + 0.5*(epsilon(ix  ,iz)+epsilon(ix  ,iz+1))*phi(ix-1,iz  )*dxsqi
     &      + 0.5*(epsilon(ix+1,iz)+epsilon(ix+1,iz+1))*phi(ix+1,iz  )*dxsqi
     &      + 0.5*(epsilon(ix,iz  )+epsilon(ix+1,iz  ))*phi(ix  ,iz-1)*dzsqi
     &      + 0.5*(epsilon(ix,iz+1)+epsilon(ix+1,iz+1))*phi(ix  ,iz+1)*dzsqi
     &                  - 0.5*(epsilon(ix  ,iz)+epsilon(ix  ,iz+1) +
     &                         epsilon(ix+1,iz)+epsilon(ix+1,iz+1))*
     &                     (dxsqi + dzsqi)*phi(ix,iz)

        enddo

        if (localbounds(1) == 0) then
          res(nxlocal:nxlocal+delx,iz) = 0.
        else if (localbounds(1)  > 0 .and. delx > 1) then
          res(nxlocal+2:nxlocal+delx,iz) = 0.
        endif

      enddo
!$OMP END DO

      if (localbounds(5) == 0) then
        res(:,nzlocal:nzlocal+delz) = 0.
      else if (localbounds(5)  > 0 .and. delz > 1) then
        res(:,nzlocal+2:nzlocal+delz) = 0.
      endif

c     --- Zero the residual inside conductors.
      call cond_potmgres(conductors%interior,
     &                   nxlocal,0,nzlocal,delx,0,delz,
     &                   res,mglevel,1,.false.)

      if (lcndbndy) then
c       --- Calculate the residual near the conductor.
        call condbndymgres2ddielectric(conductors%evensubgrid,nxlocal,nzlocal,
     &                                 phi,rho,res,epsilon,
     &                                 dxsqi,dzsqi,xminodx,lrz,localbounds,
     &                                 mglevel,delx,delz)
        call condbndymgres2ddielectric(conductors%oddsubgrid,nxlocal,nzlocal,
     &                                 phi,rho,res,epsilon,
     &                                 dxsqi,dzsqi,xminodx,lrz,localbounds,
     &                                 mglevel,delx,delz)
      endif

      call applyboundaryconditions3d(nxlocal,0,nzlocal,
     &                               delx,0,delz,
     &                               res,1,localbounds,.false.,.false.)
c     --- There is a problem with conductors and periodic boundaries. If a
c     --- conductor reaches the edge of the mesh, the residual at
c     --- the corresponding points on the opposite boundary will not be
c     --- zeroed out as it should be.

      return
      end
c=============================================================================
      subroutine restrict2dcellcentered(nx,nz,nxlocal,nzlocal,u,
     &                                  nxcoarse,nzcoarse,
     &                                  nxlocalcoarse,nzlocalcoarse,ucoarse,
     &                                  localbounds,localboundscoarse,
     &                                  lxoffset,lzoffset)
      integer(ISZ):: nx,nz,nxlocal,nzlocal
      integer(ISZ):: nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse
      real(kind=8):: u(0:nxlocal+1,0:nzlocal+1)
      real(kind=8):: ucoarse(0:nxlocalcoarse+1,0:nzlocalcoarse+1)
      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: localboundscoarse(0:5)
      integer(ISZ):: lxoffset,lzoffset
      
c Restrict to a coarser grid.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixcoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ),allocatable:: ixmina(:),ixmaxa(:)
      integer(ISZ),allocatable:: izmina(:),izmaxa(:)
      real(kind=8),allocatable:: wxa(:,:),wza(:,:)
      real(kind=8):: r,w,dx,dz,dxi,dzi,wx(0:3),wz(0:3)

      allocate(ixmina(1:nxlocalcoarse),ixmaxa(1:nxlocalcoarse))
      allocate(izmina(1:nzlocalcoarse),izmaxa(1:nzlocalcoarse))
      allocate(wxa(0:3,1:nxlocalcoarse))
      allocate(wza(0:3,1:nzlocalcoarse))
      wxa = 0.
      wza = 0.

      dx = 1.*nx/nxcoarse
      dz = 1.*nz/nzcoarse
      dxi = 1.*nxcoarse/nx
      dzi = 1.*nzcoarse/nz

c     --- Precalculate the loop limits and weights. This saves a surprisingly
c     --- substantial amount of time.
      do izcoarse=1,nzlocalcoarse
        izmin = ((izcoarse-1)*nz - lzoffset + 4*nzcoarse)/nzcoarse-3
        izmax = ((izcoarse+1)*nz - lzoffset - 1)/nzcoarse
        if (izmin < 0) izmin = 0
        if (izmax > nzlocal+1) izmax = nzlocal+1
        if (izmax < izmin) continue

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzcoarse)*dzi)
        enddo
        izmina(izcoarse) = izmin
        izmaxa(izcoarse) = izmax
        wza(:,izcoarse) = wz
      enddo

      do ixcoarse=1,nxlocalcoarse
        ixmin = ((ixcoarse-1)*nx - lxoffset + 4*nxcoarse)/nxcoarse-3
        ixmax = ((ixcoarse+1)*nx - lxoffset - 1)/nxcoarse
        if (ixmin < 0) ixmin = 0
        if (ixmax > nx+1) ixmax = nx+1
        if (ixmax < ixmin) continue

        do ix=ixmin,ixmax
          wx(ix-ixmin) = 1. - abs(ixcoarse - ix*dxi)
        enddo
        ixmina(ixcoarse) = ixmin
        ixmaxa(ixcoarse) = ixmax
        wxa(:,ixcoarse) = wx
      enddo

c     --- Do the loops.
!$OMP DO
      do izcoarse=1,nzlocalcoarse
        izmin = izmina(izcoarse)
        izmax = izmaxa(izcoarse)
        wz(:) = wza(:,izcoarse)

        do ixcoarse=1,nxlocalcoarse
          ixmin = ixmina(ixcoarse)
          ixmax = ixmaxa(ixcoarse)
          wx(:) = wxa(:,ixcoarse)

          r = 0.
          w = 0.
          do iz=izmin,izmax
              do ix=ixmin,ixmax
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

      deallocate(ixmina,ixmaxa)
      deallocate(izmina,izmaxa)
      deallocate(wxa,wza)

      return
      end
c=============================================================================
      subroutine condbndymg2ddielectric(subgrid,nxlocal,nzlocal,phi,rho,epsilon,
     &                                  dxsqi,dzsqi,xminodx,lrz,
     &                                  mgparam,localbounds,mglevel)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nzlocal,mglevel
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1), rho(0:nxlocal,0:nzlocal)
      real(kind=8):: epsilon(0:nxlocal+1,0:nzlocal+1)
      real(kind=8):: dxsqi,dzsqi,mgparam,xminodx
      logical(ISZ):: lrz
      integer(ISZ):: localbounds(0:5)

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.

      integer(ISZ):: ic,ix,iz
      real(kind=8):: pxm,pzm,pxp,pzp
      real(kind=8):: dxm,dzm,dxp,dzp
      real(kind=8):: epxm,epzm,epxp,epzp
      real(kind=8):: cxm,czm,cxp,czp
      real(kind=8):: voltfac,c0
      real(kind=8):: rr(0:nxlocal),rri(0:nxlocal)
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- Only use actual voltage on finest level. Set to zero for
c     --- coarser levels since solver for the residuals.
      if (mglevel == 0) then
        voltfac = 1.
      else
        voltfac = 0.
      endif

c     --- Note that in the loop below, for the Neumann boundary conditions,
c     --- the potentials are not set since they will be multiplied by zero
c     --- anyway. The code here just ensures that the variables have been
c     --- initialized so that debuggers and valgrind won't complain about
c     --- using uninitialized variables.
      pxm = 0.
      pxp = 0.
      pzm = 0.
      pzp = 0.

c     --- Precalculate the radius for efficiency
      if (lrz) then
        do ix=0,nxlocal
          rr(ix) = ix + xminodx
          if (rr(ix) > 0.) rri(ix) = 1./rr(ix)
        enddo
      endif

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet or parallel boundary
        if (ix == 0  .and. localbounds(0) < 1) cycle
        if (ix == nxlocal .and. localbounds(1) < 1) cycle
        if (iz == 0  .and. localbounds(4) < 1) cycle
        if (iz == nzlocal .and. localbounds(5) < 1) cycle

c       --- First, get the potential and effective grid cell sizes
c       --- Note that for the Neumann case, the potential is not
c       --- used and so is not set.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          pxm = voltfac*volt(0,ic)
          dxm = dels(0,ic)
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
        else
          pxm = phi(ix-1,iz  )
          dxm = 1.
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          pxp = voltfac*volt(1,ic)
          dxp = dels(1,ic)
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
        else
          pxp = phi(ix+1,iz  )
          dxp = 1.
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          pzm = voltfac*volt(4,ic)
          dzm = dels(4,ic)
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
        else
          pzm = phi(ix  ,iz-1)
          dzm = 1.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          dzp = dels(5,ic)
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
        else
          pzp = phi(ix  ,iz+1)
          dzp = 1.
        endif

c       --- Setup the epsilon coefficients
        epxm = 0.5*(epsilon(ix  ,iz) + epsilon(ix  ,iz+1))
        epxp = 0.5*(epsilon(ix+1,iz) + epsilon(ix+1,iz+1))
        epzm = 0.5*(epsilon(ix,iz  ) + epsilon(ix+1,iz  ))
        epzp = 0.5*(epsilon(ix,iz+1) + epsilon(ix+1,iz+1))

c       --- Now construct the coefficients
        cxm = epxm*dxsqi/(dxm*(0.5*dxm + 0.5*dxp))
        cxp = epxp*dxsqi/(dxp*(0.5*dxm + 0.5*dxp))
        czm = epzm*dzsqi/(dzm*(0.5*dzm + 0.5*dzp))
        czp = epzp*dzsqi/(dzp*(0.5*dzm + 0.5*dzp))
        if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) cxm = 0.
        if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) cxp = 0.
        if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) czm = 0.
        if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) czp = 0.

c       --- Correct coefficients for axisymmetric case
        if (lrz) then
          if (rr(ix) > 0.) then
            cxm = cxm*(rr(ix) - 0.5*dxm)*rri(ix)
            cxp = cxp*(rr(ix) + 0.5*dxp)*rri(ix)
          else
            cxm = 0.
            cxp = 4.*cxp
          endif
        endif

c       --- Note that it is possible for c0 to be zero, but that is an
c       --- isolated case which doesn't make sense anyway, so let a NaN
c       --- happen rather than add an extra check on c0.

        c0 = cxm + cxp + czm + czp
        phi(ix,iz) = mgparam*(rho(ix,iz)
     &                        + cxm*pxm + cxp*pxp + czm*pzm + czp*pzp)/c0
     &               + (1. - mgparam)*subgrid%prevphi(ic)
c       print*,phi(ix,iz),rho(ix,iz),cxm,pxm,cxp,pxp,czm,pzm,czp,pzp,c0,
c    &               subgrid%prevphi(ic)

      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine condbndymgres2ddielectric(subgrid,nxlocal,nzlocal,phi,rho,res,
     &                                     epsilon,dxsqi,dzsqi,xminodx,lrz,
     &                                     localbounds,mglevel,
     &                                     resdelx,resdelz)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nxlocal,nzlocal,mglevel
      integer(ISZ):: resdelx,resdelz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1), rho(0:nxlocal,0:nzlocal)
      real(kind=8):: res(-resdelx:nxlocal+resdelx,-resdelz:nzlocal+resdelz)
      real(kind=8):: epsilon(0:nxlocal+1,0:nzlocal+1)
      real(kind=8):: dxsqi,dzsqi,xminodx
      logical(ISZ):: lrz
      integer(ISZ):: localbounds(0:5)

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.

      integer(ISZ):: ic,ix,iz
      real(kind=8):: pxm,pzm,pxp,pzp
      real(kind=8):: dxm,dzm,dxp,dzp
      real(kind=8):: cxm,czm,cxp,czp
      real(kind=8):: epxm,epzm,epxp,epzp
      real(kind=8):: voltfac,ppp,c0
      real(kind=8):: rr(0:nxlocal),rri(0:nxlocal)
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- Only use actual voltage on finest level. Set to zero for
c     --- coarser levels since solver for the residuals.
      if (mglevel == 0) then
        voltfac = 1.
      else
        voltfac = 0.
      endif

c     --- Note that in the loop below, for the Neumann boundary conditions,
c     --- the potentials are not set since they will be multiplied by zero
c     --- anyway. The code here just ensures that the variables have been
c     --- initialized so that debuggers and valgrind won't complain about
c     --- using uninitialized variables.
      pxm = 0.
      pxp = 0.
      pzm = 0.
      pzp = 0.

c     --- Precalculate the radius for efficiency
      if (lrz) then
        do ix=0,nxlocal
          rr(ix) = ix + xminodx
          if (rr(ix) > 0.) rri(ix) = 1./rr(ix)
        enddo
      endif

c     --- loop over points near surface of conductors
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet
        if (ix == 0  .and. localbounds(0) == 0) cycle
        if (ix == nxlocal .and. localbounds(1) == 0) cycle
        if (iz == 0  .and. localbounds(4) == 0) cycle
        if (iz == nzlocal .and. localbounds(5) == 0) cycle

c       --- First, get the potential and effective grid cell sizes
c       --- Note that for the Neumann case, the potential is not
c       --- used and so is not set.
        ppp = 1.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          pxm = voltfac*volt(0,ic)
          dxm = dels(0,ic)
          ppp = min(ppp,dels(0,ic))
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
          if (abs(dels(0,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(0,ic)))
          endif
        else
          pxm = phi(ix-1,iz  )
          dxm = 1.
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          pxp = voltfac*volt(1,ic)
          dxp = dels(1,ic)
          ppp = min(ppp,dels(1,ic))
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
          if (abs(dels(1,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(1,ic)))
          endif
        else
          pxp = phi(ix+1,iz  )
          dxp = 1.
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          pzm = voltfac*volt(4,ic)
          dzm = dels(4,ic)
          ppp = min(ppp,dels(4,ic))
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
          if (abs(dels(4,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(4,ic)))
          endif
        else
          pzm = phi(ix  ,iz-1)
          dzm = 1.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          dzp = dels(5,ic)
          ppp = min(ppp,dels(5,ic))
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
          if (abs(dels(5,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(5,ic)))
          endif
        else
          pzp = phi(ix  ,iz+1)
          dzp = 1.
        endif

c       --- Setup the epsilon coefficients
        epxm = 0.5*(epsilon(ix  ,iz) + epsilon(ix  ,iz+1))
        epxp = 0.5*(epsilon(ix+1,iz) + epsilon(ix+1,iz+1))
        epzm = 0.5*(epsilon(ix,iz  ) + epsilon(ix+1,iz  ))
        epzp = 0.5*(epsilon(ix,iz+1) + epsilon(ix+1,iz+1))

c       --- Now construct the coefficients
        cxm = epxm*dxsqi/(dxm*(0.5*dxm + 0.5*dxp))
        cxp = epxp*dxsqi/(dxp*(0.5*dxm + 0.5*dxp))
        czm = epzm*dzsqi/(dzm*(0.5*dzm + 0.5*dzp))
        czp = epzp*dzsqi/(dzp*(0.5*dzm + 0.5*dzp))
        if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) cxm = 0.
        if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) cxp = 0.
        if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) czm = 0.
        if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) czp = 0.

c       --- Correct coefficients for axisymmetric case
        if (lrz) then
          if (rr(ix) > 0.) then
            cxm = cxm*(rr(ix) - 0.5*dxm)*rri(ix)
            cxp = cxp*(rr(ix) + 0.5*dxp)*rri(ix)
          else
            cxm = 0.
            cxp = 4.*cxp
          endif
        endif

        c0 = cxm + cxp + czm + czp
        res(ix,iz) = ppp*(rho(ix,iz) + cxm*pxm + cxp*pxp + czm*pzm + czp*pzp
     &                    - c0*phi(ix,iz))

        if (0. >= dels(0,ic) .and. dels(0,ic) >= -1.) res(ix-1,iz) = 0.
        if (0. >= dels(1,ic) .and. dels(1,ic) >= -1.) res(ix+1,iz) = 0.
        if (0. >= dels(4,ic) .and. dels(4,ic) >= -1.) res(ix,iz-1) = 0.
        if (0. >= dels(5,ic) .and. dels(5,ic) >= -1.) res(ix,iz+1) = 0.

      enddo
!$OMP END DO

      return
      end
c=============================================================================
