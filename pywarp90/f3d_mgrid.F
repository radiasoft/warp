#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID.M, version $Revision: 3.87 $, $Date: 2004/12/22 22:03:12 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid field sovler which is part of the F3D
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigrid3df(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                        rstar,linbend,
     &                        bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                        xmmin,ymmin,zmmin,zbeam,zgrid)
      use GlobalVars
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use Parallel
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

      call multigrid3dsolve(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                            rstar,linbend,bounds,
     &                            xmmin,ymmin,zmmin,zbeam,zgrid,
     &                            mgparam,mgform,mgiters,mgmaxiters,
     &                            mgmaxlevels,mgerror,mgtol,
     &                            downpasses,uppasses,
     &                            lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                            gridmode,conductors,
     &                            my_index,nslaves,izfsslave,nzfsslave)

      return
      end
c=============================================================================
      subroutine multigrid3dsolve(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &                            rstar,linbend,bounds,
     &                            xmmin,ymmin,zmmin,zbeam,zgrid,
     &                            mgparam,mgform,mgiters,mgmaxiters,
     &                            mgmaxlevels,mgerror,mgtol,
     &                            downpasses,uppasses,
     &                            lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                            gridmode,conductors,
     &                            my_index,nslaves,izfsslave,nzfsslave)
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
      integer(ISZ):: mgform,mgiters,mgmaxiters,mgmaxlevels
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor,lbuildquads
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

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

      integer(ISZ):: nxy,nxyz
c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
      real(kind=8):: maxres
#endif
      real(kind=8):: dxsqi,dysqi,dzsqi,reps0c,rdel
      integer(ISZ):: i,ii,k,ix,iy,iz
      real(kind=8):: rs,x,r
      real(kind=8):: phisave(0:nx,0:ny,-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
c     --- The following only used when mgform == 2
      real(kind=8),allocatable:: rhosave(:,:,:),res(:,:,:)
      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: lzoffset(0:nslaves-1),rzoffset(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleftc(0:nslaves-1), izsendingleftc(0:nslaves-1)
      integer(ISZ):: whosendingrightc(0:nslaves-1),izsendingrightc(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      character(72):: errline
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Initialize temporaries
      nxy    = (nx+1)*(ny+1)
      nxyz   = (nx+1)*(ny+1)*(nz+1)
      dxsqi  = 1./dx**2
      dysqi  = 1./dy**2
      dzsqi  = 1./dz**2
      reps0c = mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
      rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

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
        do iz=0,nz
          rs = rstar(iz)
          do ix=0,nx
            x  = xmmin + ix*dx
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
        do ii = 1, nxy
          bendx(ii) = xmmin + mod(ii-1,nx+1)*dx
        enddo
!$OMP END DO
c       --- Change rstar if using Nuemann boundary conditions
#ifndef MPIPARALLEL
        if (bounds(4) == 1) rstar(-1) = rstar(1)
        if (bounds(5) == 1) rstar(nz+1) = rstar(nz-1)
#else
        if (bounds(4) == 1 .and. my_index == 0) rstar(-1) = rstar(1)
        if (bounds(5) == 1 .and. my_index == nslaves-1) rstar(nz+1)=rstar(nz-1)
#endif
      endif

c     --- If using residual correction form, need to save the original rho.
c     --- Also setup parallel arrays.
      if (mgform == 2) then
!$OMP SINGLE
        allocate(rhosave(0:nx,0:ny,0:nz),res(-1:nx+1,-1:ny+1,-3:nz+3))
!$OMP END SINGLE
        rhosave = rho
        localbounds = bounds
#ifdef MPIPARALLEL
        lzoffset = 0
        rzoffset = 0
        call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                        bounds,nzfull,
     &                        lzoffset,rzoffset,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
        if (izfsslave(my_index) > 0) localbounds(4) = -1
        if (izfsslave(my_index)+nz < nzfull) localbounds(5) = -1
#endif
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
          call cond_potmg(conductors%interior,
     &                    nx,ny,nz,phisave,0,mgform,.true.)
          call residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phisave,rhosave,res,
     &                  0,localbounds,mgparam,mgform,.true.,
     &                  lcndbndy,icndbndy,conductors)
#ifdef MPIPARALLEL
          call mgexchange_phi(nx,ny,nz,nzfull,res,localbounds,-1,1,3,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
          call mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localbounds,0,1,3,
     &                                my_index,nslaves,izfsslave,
     &                                whosendingleft,whosendingright)
#endif
          rho = res(0:nx,0:ny,0:nz)
          phi = 0.
        endif

c       --- Do one vcycle.
        call vcycle(0,nx,ny,nz,nzfull,dx,dy,dz,phi,rho,
     &              rstar,linbend,bendx,bounds,mgparam,mgform,mgmaxlevels,
     &              downpasses,uppasses,lcndbndy,icndbndy,conductors,
     &              my_index,nslaves,izfsslave,nzfsslave)

c       --- If using residual correction form, add the resulting error to phi.
        if (mgform == 2) phi = phi + phisave

c       --- When using residual correction form, the other planes do need
c       --- to be set when using other than Dirichlet boundaries since
c       --- those planes are only set with the error of phi.
        if (mgform == 2) then
          if (localbounds(4) == 1) phi(:,:,-1) = phi(:,:,1)
          if (localbounds(5) == 1) phi(:,:,nz+1) = phi(:,:,nz-1)
          if (localbounds(4) == 2 .and. nz==nzfull) phi(:,:,-1) = phi(:,:,nz-1)
          if (localbounds(5) == 2 .and. nz==nzfull) phi(:,:,nz+1) = phi(:,:,1)
#ifdef MPIPARALLEL
          call mgexchange_phi(nx,ny,nz,nzfull,phi,localbounds,0,0,1,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
          call mgexchange_phi(nx,ny,nz,nzfull,phi,localbounds,-1,0,1,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
#endif
      endif

c       --- Calculate the change in phi.
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

c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
        maxres = 0.
        if (mgform == 2) then
!$OMP DO REDUCTION(MAX:maxres)
           do i=0,nz
              maxres = max(maxres, maxval(abs(res(:,:,iz))))
           enddo
!$OMP END DO
        endif
#endif

#ifdef MPIPARALLEL
        if (nslaves > 1) then
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

c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
      print *, 'iteration ', mgiters, ' mgerror=', mgerror
c     unscaled residual
      print *, 'max(residual)=', maxres*2.0*(dxsqi+dysqi+dzsqi)
#endif

      enddo

c     --- For Dirichlet boundary conditions, copy data into guard planes
c     --- For other boundary conditions, the guard planes are used during
c     --- the solve are so are already set.
      if (bounds(4) == 0) phi(:,:,-1) = phi(:,:,0)
      if (bounds(5) == 0) phi(:,:,nz+1) = phi(:,:,nz)

c     --- Make a print out.
      if (mgerror > mgtol) then
        call remark("Multigrid: Maximum number of iterations reached")
      endif
      write(errline,20) mgerror,mgiters
  20  format("Multigrid: Error converged to ",1pe11.3," in ",i5," v-cycles")
      call remark(errline)

c     --- If using residual correction form, restore saved rho
      if (mgform == 2) then
        rho = rhosave
        deallocate(rhosave,res)
      endif

c     --- Restore rho
      reps0c = 1./reps0c
      if (.not. linbend) then
        rho = rho*reps0c
      else
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
        do iz=0,nz
          rs = rstar(iz)
          do ix=0,nx
            x  = xmmin + ix*dx
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
      function mggetarraysum(nx,ny,nz,delt,delz,a,
     &                       izfsslave,nzfsslave,npes,my_index)
      real(kind=8),dimension(2):: mggetarraysum
      integer(ISZ):: nx,ny,nz,delt,delz,npes,my_index
      integer(ISZ):: izfsslave(0:npes-1),nzfsslave(0:npes-1)
      real(kind=8):: a(-delt:nx+delt,-delt:ny+delt,-delz:nz+delz)

#ifdef MPIPARALLEL

      integer(ISZ):: i1,i2
      real(kind=8):: sss(2)

      i1 = -delz
      if (my_index < npes-1) then
        i2 = izfsslave(my_index+1) - izfsslave(my_index) - 1 - delz
      else
        i2 = nz + delz
      endif
      sss(1) = sum(sum(sum(a(:,:,i1:i2),1),1),1)

      if (my_index > 0) then
        i1 = izfsslave(my_index-1) + nzfsslave(my_index-1) -
     &       izfsslave(my_index) + 1 + delz
      else
        i1 = 0 - delz
      endif
      i2 = nz + delz
      sss(2) = sum(sum(sum(a(:,:,i1:i2),1),1),1)

      if (npes > 1) call parallelsumrealarray(sss,2)
      mggetarraysum = sss

#else

      mggetarraysum(1) = sum(sum(sum(a,1),1),1)
      mggetarraysum(2) = mggetarraysum(1)

#endif

      return
      end function mggetarraysum
      end module formggetarraysuminterface
c=============================================================================
      RECURSIVE subroutine vcycle(mglevel,nx,ny,nz,nzfull,dx,dy,dz,
     &                            phi,rho,rstar,linbend,bendx,globalbounds,
     &                            mgparam,mgform,
     &                            mgmaxlevels,downpasses,uppasses,
     &                            lcndbndy,icndbndy,conductors,
     &                            my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
c     use formggetarraysuminterface
      integer(ISZ):: mglevel
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      logical(ISZ):: linbend
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
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
      real(kind=8):: res(-1:nx+1,-1:ny+1,-3:nz+3)
      integer(ISZ):: i
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
c     real(kind=8):: sss(2)

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

c     added by petermc, 26 Sep 2002
#ifdef WITHCHOMBO
c      print *, 'vcycle on dimensions ', nx, ny, nz
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DOING ', downpasses,
     &     ' SORPASS3D CALLS ON ', nx, ny, nz, ' Warp'
#endif

c     sss = mggetarraysum(nx,ny,nz,0,1,phi,
c    &                    izfsslave,nzfsslave,nslaves,my_index)
c     if (my_index == 0) print*,"V1 phi",mglevel,sss
c     sss = mggetarraysum(nx,ny,nz,0,0,rho,
c    &                    izfsslave,nzfsslave,nslaves,my_index)
c     if (my_index == 0) print*,"V1 rho",mglevel,sss

c     --- Do initial SOR passes.
c     print *, 'vcycle on dimensions ', nx, ny, nzfull,my_index
c     call printarray3d(nx,ny,nz,0,0,izfsslave(my_index),rho,my_index,nslaves,
c    &                  "RHO",mglevel)
c     call printarray3d(nx,ny,nz,0,1,izfsslave(my_index),rho,my_index,nslaves,
c    &                  "RHO",mglevel)
      do i=1,downpasses
        call sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 localbounds,mgparam,mgform,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

#ifdef WITHCHOMBO
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DID ', downpasses,
     &     ' SORPASS3D CALLS ON ', nx, ny, nz, ' Warp'
#endif

c     sss = mggetarraysum(nx,ny,nz,0,1,phi,
c    &                    izfsslave,nzfsslave,nslaves,my_index)
c     if (my_index == 0) print*,"V2 phi",mglevel,sss
c     sss = mggetarraysum(nx,ny,nz,0,0,rho,
c    &                    izfsslave,nzfsslave,nslaves,my_index)
c     if (my_index == 0) print*,"V2 rho",mglevel,sss

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. ny >= 4 .and. nzfull >= 4 .and.
     &    mglevel < mgmaxlevels) then

c       --- Get the residual on the current grid.
        call residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
     &                mglevel,localbounds,mgparam,mgform,.false.,
     &                lcndbndy,icndbndy,conductors)
#ifdef MPIPARALLEL
          call mgexchange_phi(nx,ny,nz,nzfull,res,localbounds,-1,1,3,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
          call mgexchange_phi(nx,ny,nz,nzfull,res,localbounds,-2,1,3,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
          call mgexchange_phi(nx,ny,nz,nzfull,res,localbounds,-3,1,3,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
          call mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localbounds,0,1,3,
     &                                my_index,nslaves,izfsslave,
     &                                whosendingleft,whosendingright)
#endif
c         sss = mggetarraysum(nx,ny,nz,1,3,res,
c    &                        izfsslave,nzfsslave,nslaves,my_index)
c       if (my_index == 0) print*,"V3 res",mglevel,sss


c       --- Calculate the size of the next coarsest grid. If the current
c       --- size is twice and odd integer, the next coarsest size will be
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

        dxcoarse = dx*nx/nxcoarse
        dycoarse = dy*ny/nycoarse
        dzcoarse = dz*nzfull/nzfullcoarse
        dxcoarsesqi = 1./dxcoarse**2
        dycoarsesqi = 1./dycoarse**2
        dzcoarsesqi = 1./dzcoarse**2

c       --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
c       --- all axis.  This is the same check that is done in getmglevels.
c       --- dz > 4/3 dx <=> (9/16) / dx^2 < 1 / dz^2
        partialcoarsening = (dz > 4./3.*dx)
#ifdef MPIPARALLEL
c       --- This must be a global operation since, due to roundoff, each
c       --- processor can get a different value if dz == 4./3.*dx.
c       call parallellor(partialcoarsening)
#endif
        if (partialcoarsening) then

c         --- Alloate new work space
          allocate(phicoarse(0:nxcoarse,0:nycoarse,-1:nz+1),
     &             rhocoarse(0:nxcoarse,0:nycoarse,0:nz))
          phicoarse = 0.

c         --- Ratio of old to new constant needed to scale the residual for
c         --- the restriction.
          ff = (dxsqi+dysqi+dzsqi)/(dxcoarsesqi + dycoarsesqi + dzsqi)
          call restrict2d(nx,ny,nz,res,nxcoarse,nycoarse,rhocoarse,ff,
     &                    localbounds)

c         --- Continue at the next coarsest level.
          call vcycle(mglevel+1,nxcoarse,nycoarse,nz,nzfull,dxcoarse,dycoarse,
     &                dz,phicoarse,rhocoarse,
     &                rstar,linbend,bendx,globalbounds,mgparam,mgform,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors,
     &                my_index,nslaves,izfsslave,nzfsslave)

c         --- Add in resulting error.
          call expand2d(nx,ny,nz,phi,nxcoarse,nycoarse,phicoarse,localbounds)

          deallocate(phicoarse,rhocoarse)
        else

          localboundsc = globalbounds

#ifdef MPIPARALLEL
c         --- Find domains in coarser grid
          call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                    nzfull,nzfullcoarse)
c         --- Reset value to corrected one
          nzcoarse = nzfsslavec(my_index)
c         --- Difference between starts and ends of coarse and fine grids.
c         --- Should only be in the range 0-2.
          lzoffsetall = (nzfullcoarse*izfsslave-nzfull*izfsslavec)
          rzoffsetall = (nzfull*(izfsslavec + nzfsslavec) -
     &                   nzfullcoarse*(izfsslave + nzfsslave))
c         --- Note that the lzoffsetall and rzoffsetall can only be used in
c         --- MPIPARALLEL sections since they will be unallocated in the
c         --- serial code. So, separate scalars are used in code which is
c         --- used in the serial version.
          lzoffset = lzoffsetall(my_index)
          rzoffset = rzoffsetall(my_index)
c         --- Get processor with which to exchange data on coarse grid
          call mggetexchangepes(nslaves,izfsslavec,nzfsslavec,my_index,
     &                          globalbounds,nzfullcoarse,
     &                          lzoffsetall,rzoffsetall,
     &                          whosendingleftc,izsendingleftc,
     &                          whosendingrightc,izsendingrightc)
          if (izfsslavec(my_index) > 0) localboundsc(4) = -1
          if (izfsslavec(my_index) + nzcoarse < nzfullcoarse) localboundsc(5) = -1
#else
          lzoffset = 0
          rzoffset = 0
#endif

c         --- Alloate new work space
          allocate(phicoarse(0:nxcoarse,0:nycoarse,-1:nzcoarse+1),
     &             rhocoarse(0:nxcoarse,0:nycoarse,0:nzcoarse))
          rhocoarse = 0.
          phicoarse = 0.

c         --- Restriction - note that scaling factor for residual is always
c         --- 4 for full-coarsening and is compiled into the restriction
c         --- routine.
          ff = (dxsqi+dysqi+dzsqi)/(dxcoarsesqi + dycoarsesqi + dzcoarsesqi)
          call restrict3d(nx,ny,nz,nzfull,res,
     &                    nxcoarse,nycoarse,nzcoarse,nzfullcoarse,rhocoarse,ff,
     &                    localbounds,localboundsc,lzoffset)
c         sss = mggetarraysum(nxcoarse,nycoarse,nzcoarse,0,0,rhocoarse,
c    &                        izfsslavec,nzfsslavec,nslaves,my_index)
c         if (my_index == 0) print*,"V3 rhocoarse",mglevel,sss

c         --- Continue at the next coarsest level.
          call vcycle(mglevel+1,nxcoarse,nycoarse,nzcoarse,nzfullcoarse,
     &                dxcoarse,dycoarse,dzcoarse,phicoarse,rhocoarse,
     &                rstar,linbend,bendx,globalbounds,mgparam,mgform,
     &                mgmaxlevels,downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors,
     &                my_index,nslaves,izfsslavec,nzfsslavec)

c         sss = mggetarraysum(nxcoarse,nycoarse,nzcoarse,0,1,phicoarse,
c    &                        izfsslavec,nzfsslavec,nslaves,my_index)
c         if (my_index == 0) print*,"V4 phicoarse",mglevel,sss

c         --- Add in resulting error.
          call expand3d(nx,ny,nz,nzfull,phi,
     &                  nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                  localbounds,lzoffset)
#ifdef MPIPARALLEL
          call mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,localbounds,1,0,1,
     &                                my_index,nslaves,izfsslave,
     &                                whosendingleft,whosendingright)
#endif
c         sss = mggetarraysum(nx,ny,nz,0,1,phi,
c    &                        izfsslave,nzfsslave,nslaves,my_index)
c         if (my_index == 0) print*,"V5 phi",mglevel,sss

          deallocate(phicoarse,rhocoarse)
        endif

      endif

c     --- Do final SOR passes.
#ifdef WITHCHOMBO
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DOING ', uppasses,
     &     ' SORPASS3D CALLS ON ', nx, ny, nz, ' Warp'
#endif

      do i=1,uppasses
        call sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,bendx,
     &                 localbounds,mgparam,mgform,
     &                 lcndbndy,icndbndy,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

#ifdef WITHCHOMBO
c     message added by petermc, 1 Oct 2002
      if (.false.) print *, 'DID ', uppasses,
     &     ' SORPASS3D CALLS ON ', nx, ny, nz, ' Warp'
#endif

      return
      end
c=============================================================================
      subroutine restrict2d(nx,ny,nz,res,nxcoarse,nycoarse,rhocoarse,
     &                      ff,bounds)
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nxcoarse,nycoarse
      real(kind=8):: res(-1:nx+1,-1:ny+1,-3:nz+3)
      real(kind=8):: rhocoarse(0:nxcoarse,0:nycoarse,0:nz)
      real(kind=8):: ff
      integer(ISZ):: bounds(0:5)
      
c Restrict transversely to a coarser grid.  The factor of ff is needed
c since the residual will be used as the source term in the next coarser
c grid and it implicitly includes the factor 0.5/(dxsqi+dysqi+dzsqi)
c which is ff times too small with the current values of dxsqi etc.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixcoarse,iycoarse
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixcoarsemin,ixcoarsemax,iycoarsemin,iycoarsemax
      real(kind=8):: w,r,dx,dy,dxi,dyi,wx(0:3),wy(0:3)

c     --- Set the loop limits, including edges when appropriate.
      ixcoarsemin = 1
      ixcoarsemax = nxcoarse-1
      iycoarsemin = 1
      iycoarsemax = nycoarse-1
      izmin = 0
      izmax = nz
      if (bounds(0) >  0) ixcoarsemin = 0
      if (bounds(1) == 1) ixcoarsemax = nxcoarse
      if (bounds(2) >  0) iycoarsemin = 0
      if (bounds(3) == 1) iycoarsemax = nycoarse
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

      dx = 1.*nx/nxcoarse
      dy = 1.*ny/nycoarse
      dxi = 1.*nxcoarse/nx
      dyi = 1.*nycoarse/ny
      w = ff*dxi*dyi

c     --- Do the loops.
!$OMP DO
      do iz=izmin,izmax
        do iycoarse=iycoarsemin,iycoarsemax
c         --- The 1.e-10 fuzz is needed when dy is 2, in which case
c         --- (iycoarse+-1)*dy would be integers. The fuzz forces truncation
c         --- in the correct manner. The fuzz needs to be small compared to
c         --- the maximum number of grid cells ever used.
          iymin = int(ceiling((iycoarse-1)*dy + 1.e-10))
          iymax = int(floor((iycoarse+1)*dy - 1.e-10))

          do iy=iymin,iymax
            wy(iy-iymin) = 1. - abs(iycoarse - iy*dyi)
          enddo

          do ixcoarse=ixcoarsemin,ixcoarsemax
            ixmin = int(ceiling((ixcoarse-1)*dx + 1.e-10))
            ixmax = int(floor((ixcoarse+1)*dx - 1.e-10))

            do ix=ixmin,ixmax
              wx(ix-ixmin) = 1. - abs(ixcoarse - ix*dxi)
            enddo

            r = 0.
            do iy=iymin,iymax
              do ix=ixmin,ixmax
                r = r + wx(ix-ixmin)*wy(iy-iymin)*res(ix,iy,iz)
              enddo
            enddo
            rhocoarse(ixcoarse,iycoarse,iz) = w*r

          enddo
        enddo
      enddo
!$OMP END DO

c     --- Set appropriate boundary values
      if (bounds(0) == 0) rhocoarse(nxcoarse,:,:) = 0.
      if (bounds(1) == 0) rhocoarse(nxcoarse,:,:) = 0.
      if (bounds(1) == 2) rhocoarse(nxcoarse,:,:) = rhocoarse(0,:,:)
      if (bounds(2) == 0) rhocoarse(:,nycoarse,:) = 0.
      if (bounds(3) == 0) rhocoarse(:,nycoarse,:) = 0.
      if (bounds(3) == 2) rhocoarse(:,nycoarse,:) = rhocoarse(:,0,:)

      if (bounds(4) == 0) rhocoarse(:,:,0) = 0.
      if (bounds(5) == 0) rhocoarse(:,:,nz) = 0.

      return
      end
c=============================================================================
      subroutine restrict3d(nx,ny,nz,nzfull,res,
     &                      nxcoarse,nycoarse,nzcoarse,nzfullcoarse,rhocoarse,
     &                      ff,bounds,boundscoarse,lzoffset)
      integer(ISZ):: nx,ny,nz,nzfull
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: res(-1:nx+1,-1:ny+1,-3:nz+3)
      real(kind=8):: rhocoarse(0:nxcoarse,0:nycoarse,0:nzcoarse)
      real(kind=8):: ff
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: boundscoarse(0:5)
      integer(ISZ):: lzoffset
      
c Restrict to a coarser grid.  The factor ff is needed since the
c residual will be used as the source term in the next coarser grid and
c it implicitly includes the factor 0.5/(dxsqi+dysqi+dzsqi) which is
c too small with the current values of dxsqi etc.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixcoarse,iycoarse,izcoarse
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixcoarsemin,ixcoarsemax,iycoarsemin,iycoarsemax
      integer(ISZ):: izcoarsemin,izcoarsemax
      real(kind=8):: r,w,dx,dy,dz,dxi,dyi,dzi,wx(0:3),wy(0:3),wz(0:3)

c     --- Set the loop limits, including edges when appropriate.
      ixcoarsemin = 1
      ixcoarsemax = nxcoarse-1
      iycoarsemin = 1
      iycoarsemax = nycoarse-1
      izcoarsemin = 0
      izcoarsemax = nzcoarse
      if (bounds(0) >  0) ixcoarsemin = 0
      if (bounds(1) == 1) ixcoarsemax = nxcoarse
      if (bounds(2) >  0) iycoarsemin = 0
      if (bounds(3) == 1) iycoarsemax = nycoarse
      if (bounds(4) == 0) izcoarsemin = 1
      if (bounds(5) == 0) izcoarsemax = nzcoarse - 1

      dx = 1.*nx/nxcoarse
      dy = 1.*ny/nycoarse
      dz = 1.*nzfull/nzfullcoarse
      dxi = 1.*nxcoarse/nx
      dyi = 1.*nycoarse/ny
      dzi = 1.*nzfullcoarse/nzfull
      w = ff*dxi*dyi*dzi

c     --- Do the loops.
!$OMP DO
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = ((izcoarse-1)*nzfull - lzoffset + 4*nzfullcoarse)/nzfullcoarse-3
        izmax = ((izcoarse+1)*nzfull - lzoffset - 1)/nzfullcoarse

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzfullcoarse)*dzi)
        enddo

        do iycoarse=iycoarsemin,iycoarsemax
          iymin = int(ceiling((iycoarse-1)*dy + 1.e-10))
          iymax = int(floor((iycoarse+1)*dy - 1.e-10))

          do iy=iymin,iymax
            wy(iy-iymin) = 1. - abs(iycoarse - iy*dyi)
          enddo

          do ixcoarse=ixcoarsemin,ixcoarsemax
            ixmin = int(ceiling((ixcoarse-1)*dx + 1.e-10))
            ixmax = int(floor((ixcoarse+1)*dx - 1.e-10))

            do ix=ixmin,ixmax
              wx(ix-ixmin) = 1. - abs(ixcoarse - ix*dxi)
            enddo

            r = 0.
            do iz=izmin,izmax
              do iy=iymin,iymax
                do ix=ixmin,ixmax
                  r = r + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)*res(ix,iy,iz)
                enddo
              enddo
            enddo
            rhocoarse(ixcoarse,iycoarse,izcoarse) = w*r

          enddo
        enddo
      enddo
!$OMP END DO

c     --- Set appropriate boundary values
      if (bounds(0) == 0) rhocoarse(nxcoarse,:,:) = 0.
      if (bounds(1) == 0) rhocoarse(nxcoarse,:,:) = 0.
      if (bounds(1) == 2) rhocoarse(nxcoarse,:,:) = rhocoarse(0,:,:)
      if (bounds(2) == 0) rhocoarse(:,nycoarse,:) = 0.
      if (bounds(3) == 0) rhocoarse(:,nycoarse,:) = 0.
      if (bounds(3) == 2) rhocoarse(:,nycoarse,:) = rhocoarse(:,0,:)

      if (boundscoarse(4) == 0) rhocoarse(:,:,0) = 0.
      if (boundscoarse(5) == 0) rhocoarse(:,:,nzcoarse) = 0.

      return
      end
c=============================================================================
      subroutine expand2d(nx,ny,nz,phi,nxcoarse,nycoarse,phicoarse,bounds)
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nxcoarse,nycoarse
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: phicoarse(0:nxcoarse,0:nycoarse,-1:nz+1)
      integer(ISZ):: bounds(0:5)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: izmin,izmax
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: jx,jy,jz
      real(kind=8):: wx,wy,dx,dy

c     --- Set the loop limits, including edges when appropriate.
      izmin = -1
      izmax = nz+1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

      dx = 1.*nxcoarse/nx
      dy = 1.*nycoarse/ny

!$OMP DO
      do iz=izmin,izmax
        do iy=0,ny-1
          jy = int(iy*dy)
          wy =     iy*dy - jy
          do ix=0,nx-1
            jx = int(ix*dx)
            wx =     ix*dx - jx

            phi(ix,iy,iz) = phi(ix,iy,iz) +
     &             (1.-wx)*(1.-wy)*phicoarse(jx  ,jy  ,iz) +
     &                 wx *(1.-wy)*phicoarse(jx+1,jy  ,iz) +
     &             (1.-wx)*    wy *phicoarse(jx  ,jy+1,iz) +
     &                 wx *    wy *phicoarse(jx+1,jy+1,iz)
          enddo
        enddo
      enddo

      if (bounds(1) > 0 .or. bounds(3) > 0) then
        do iz=izmin,izmax
          phi(nx,ny,iz) = phi(nx,ny,iz) + phicoarse(nxcoarse,nycoarse,iz)
          if (bounds(1) > 0) then
            do iy=0,ny-1
              jy = int(iy*dy)
              wy =     iy*dy - jy
              phi(nx,iy,iz) = phi(nx,iy,iz) +
     &                       (1.-wy)*phicoarse(nxcoarse,jy  ,iz) +
     &                           wy *phicoarse(nxcoarse,jy+1,iz)
            enddo
          endif
          if (bounds(3) > 0) then
            do ix=0,nx-1
              jx = int(ix*dx)
              wx =     ix*dx - jx
              phi(ix,ny,iz) = phi(ix,ny,iz) +
     &                       (1.-wx)*phicoarse(jx  ,nycoarse,iz) +
     &                           wx *phicoarse(jx+1,nycoarse,iz)
            enddo
          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine expand3d(nx,ny,nz,nzfull,phi,
     &                    nxcoarse,nycoarse,nzcoarse,nzfullcoarse,phicoarse,
     &                    bounds,lzoffset)
      integer(ISZ):: nx,ny,nz,nzfull
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: phicoarse(0:nxcoarse,0:nycoarse,-1:nzcoarse+1)
      integer(ISZ):: lzoffset
      integer(ISZ):: bounds(0:5)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: izmin,izmax
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: jx,jy,jz
      real(kind=8):: dx,dy,dz
      real(kind=8):: wx,wy,wz

c     --- Set the loop limits, including edges when appropriate.
      izmin = 0
      izmax = nz

      dx = 1.*nxcoarse/nx
      dy = 1.*nycoarse/ny
      dz = 1.*nzfullcoarse/nzfull

!$OMP DO
      do iz=izmin,izmax
        jz = int((iz*nzfullcoarse + lzoffset)/nzfull)
        wz =  1.*(iz*nzfullcoarse + lzoffset)/nzfull - jz
        do iy=0,ny-1
          jy = int(iy*dy)
          wy =     iy*dy - jy
          do ix=0,nx-1
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
      subroutine sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                     rdx2,rdy2,rdz2,linbend,bendx,bounds,
     &                     mgparam,mgform,lcndbndy,icndbndy,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: mglevel,nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz),rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      real(kind=8):: rdx2,rdy2,rdz2
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
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

c     --- Put desired potential onto conductors in phi array.
      call cond_potmg(conductors%interior,
     &                nx,ny,nz,phi,mglevel,mgform,.false.)

c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      s_parity = mod(izfsslave(my_index),2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

        call sorhalfpass3d(parity,mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                     rdx2,rdy2,rdz2,linbend,bendx,
     &                     bounds,mgparam,mgform,
     &                     lcndbndy,icndbndy,conductors)

        if (bounds(4) == 2 .and. nz == nzfull) phi(:,:,-1) = phi(:,:,nz-1)
        if (bounds(5) == 2 .and. nz == nzfull) phi(:,:,nz:nz+1) = phi(:,:,0:1)
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
c         --- Exchange phi in the z guard planes
          call mgexchange_phi(nx,ny,nz,nzfull,phi,bounds,-1,0,1,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
#endif

      return
      end
c=============================================================================
      subroutine sorhalfpass3d(parity,mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                     rdx2,rdy2,rdz2,linbend,bendx,bounds,
     &                     mgparam,mgform,lcndbndy,icndbndy,conductors)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: parity,mglevel,nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz),rstar(-1:nz+1)
      real(kind=8):: bendx((nx+1)*(ny+1))
      real(kind=8):: rdx2,rdy2,rdz2
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
      logical(ISZ):: lcndbndy
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
c Note that loops over all directions assume that nx and ny are even.
c
c The arrangement of the loops was done to increase performance.  The entire
c grid is looped over as if it were a 1D array, ignoring boundaries.
c The boundaries are then reset, the previous value was destroyed.
c
c rstar(-1) and rstar(nz+1) are set based on the axial boundary conditions.

      real(kind=8):: rdel,const,rdx2c,rdy2c,rdz2c,spm1,dx
      real(kind=8):: boundsavex(0:nx,2,0:nz)
      real(kind=8):: boundsavey(0:ny,2,0:nz)
      integer(ISZ):: nxy,iimx,iipx,iimy,iipy,iimz,iipz,izl,izr
      integer(ISZ):: ii,ix,iy,iz,ic,i1,i2

c     --- Set temporary variables (these are used to increase performance)
      dx = 1./sqrt(rdx2)
      rdel = rdz2/(rdx2 + rdy2 + rdz2)
      const = mgparam*0.5/(rdx2 + rdy2 + rdz2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdz2c = rdz2*const
      spm1 = 1. - mgparam

c     --- Set indices for 1d arrays used in the seven point finite difference
c     --- form of Poisson's equation.
      nxy = (nx+1)*(ny+1)
      iimx = -1
      iipx = +1
      iimy = -nx-1
      iipy = +nx+1
      iimz = -nxy
      iipz = +nxy

c     --- Set longitudinal indices for 1d phi array.  '1' is added to izl and
c     --- izr since phi1d is passed as phi(0,0,-1).  If using Dirichlet
c     --- boundary conditions, do not solve for the potential on the end planes.
      izl = 1
      izr = nz + 1
      if (bounds(4) < 1) izl = izl + 1
      if (bounds(5) < 1) izr = izr - 1

c     --- Save values on the transverse boundaries.
      do iz=izl-1,izr-1
        do ix=mod(iz+parity,2),nx,2
          boundsavex(ix,1,iz) = phi(ix,0,iz)
          boundsavex(ix,2,iz) = phi(ix,ny,iz)
        enddo
        do iy=mod(iz+parity,2),ny,2
          boundsavey(iy,1,iz) = phi(0,iy,iz)
          boundsavey(iy,2,iz) = phi(nx,iy,iz)
        enddo
      enddo

c     --- Save values just outside conductor surfaces. Only save phi at the
c     --- subgrid points which are to be used at the current level of
c     --- grid refinement.
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

c     --- guard planes in z are already set

c     --- Loop over the rest of the array. Boundary points are calculated
c     --- too, even though the equation is wrong.  They are recalculated
c     --- later.  The case with and without bends are seperated since
c     --- the loops are done differently.  Without bends, the whole array
c     --- can be processed in one loop as a 1d array.  With bends, an
c     --- outer loop over z is needed.
      if (.not. linbend) then 
c       --- no bends 
        call mgsor_loop(parity,1,0,0,1,izl*nxy+2,(izr+1)*nxy,2,
     &                  iimx,iipx,iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      else 
c       --- bends
        call mgsor_loop(parity,1,izl,izr,1,2,nxy,2,
     &                  iimx,iipx,iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif 

c     --- Transverse boundaries     
c     --- Restore only even or odd boundary points to previous value
c     --- since only want to restore the values changed from the 1d
c     --- loop.  This automatically takes care of Dirichlet boundaries.
      do iz=izl-1,izr-1
        do ix=mod(iz+parity,2),nx,2
          phi(ix,0,iz) = boundsavex(ix,1,iz)
          phi(ix,ny,iz) = boundsavex(ix,2,iz)
        enddo
        do iy=mod(iz+parity,2),ny,2
          phi(0,iy,iz) = boundsavey(iy,1,iz)
          phi(nx,iy,iz) = boundsavey(iy,2,iz)
        enddo
      enddo

c     --- Neumann boundary conditions for the transverse planes
      if (bounds(0) == 1) then
        call mgsor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,2*(nx+1),
     &                  -iimx,iipx,iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(1) == 1) then
        call mgsor_loop(parity,nx+1,izl,izr,1,2*(nx+1),(nx+1)*ny,2*(nx+1),
     &                  iimx,-iipx,iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(2) == 1) then
        call mgsor_loop(parity,1,izl,izr,1,2,nx,2,
     &                  iimx,iipx,-iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(3) == 1) then
        call mgsor_loop(parity,1,izl,izr,1,(nx+1)*ny+2,(nx+1)*(ny+1)-1,2,
     &                  iimx,iipx,iimy,-iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif

c     --- lines at transverse edges
      if (bounds(0) == 1 .or. bounds(2) == 1) then
        call mgsor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,2*nxy,
     &                  -iimx,iipx,-iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(1) == 1 .or. bounds(2) == 1) then
        call mgsor_loop(parity,nxy,izl,izl,1,nx+1,nxy*(izr-izl)+nx+1,2*nxy,
     &                  iimx,-iipx,-iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(0) == 1 .or. bounds(3) == 1) then
        call mgsor_loop(parity,nxy,izl,izl,1,nxy-nx,nxy*(izr-izl)+nxy-nx,2*nxy,
     &                  -iimx,iipx,iimy,-iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(1) == 1 .or. bounds(3) == 1) then
        call mgsor_loop(parity,nxy,izl,izl,1,nxy,nxy*(izr-izl)+nxy,2*nxy,
     &                  iimx,-iipx,iimy,-iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif

c     --- Periodic boundary conditions for the transverse planes
      if (bounds(0) == 2) then
        call mgsor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,2*(nx+1),
     &                  nx-1,iipx,iimy,iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(2) == 2) then
        call mgsor_loop(parity,1,izl,izr,1,2,nx,2,
     &                  iimx,iipx,(nx+1)*(ny-1),iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
      endif
      if (bounds(1) == 2 .or. bounds(3) == 2) then
        do iz=izl-1,izr-1
          if (bounds(3) == 2) then
            do ix=1+mod(iz+parity+1,2),nx-1,2
              phi(ix,ny,iz) = phi(ix,0,iz)
            enddo
          endif
          if (bounds(1) == 2) then
            do iy=1+mod(iz+parity+1,2),ny-1,2
              phi(nx,iy,iz) = phi(0,iy,iz)
            enddo
          endif
        enddo
      endif
c     --- lines at transverse edges
      if (ANY(bounds(0:3) == 2)) then
        call mgsor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,2*nxy,
     &                  nx-1,iipx,(nx+1)*(ny-1),iipy,iimz,iipz,
     &                  phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend,
     &                  rstar,dx,rdx2,rdy2,rdz2,rdel,const,bendx)
        do iz=izl-1+mod(izl-1+parity,2),izr-1,2
          phi(nx,0,iz) = phi(0,0,iz)
          phi(0,ny,iz) = phi(0,0,iz)
          phi(nx,ny,iz) = phi(0,0,iz)
        enddo
      endif
c     --- end of transverse boundaries

c     --- Restore even or odd boundary points to previous value
c     --- along the transverse corners of the mesh when applying
c     --- Dirichlet boundary conditions there. Otherwise, the corners
c     --- will be left with Nuemann or periodic boundaries, which is not
c     --- quite correct.
      if (ANY(bounds(0:3)==0)) then
        do iz=izl-1,izr-1
          if (mod(iz+parity,2) == 0) then
            if (bounds(0) == 0) then
              phi(0,0,iz) = boundsavex(0,1,iz)
              phi(nx,0,iz) = boundsavex(nx,1,iz)
            endif
            if (bounds(1) == 0) then
              phi(0,ny,iz) = boundsavex(0,2,iz)
              phi(nx,ny,iz) = boundsavex(nx,2,iz)
            endif
          endif
          if (mod(iz+parity,2) == 0) then
            if (bounds(2) == 0) then
              phi(0,0,iz) = boundsavey(0,1,iz)
              phi(0,ny,iz) = boundsavey(ny,1,iz)
            endif
            if (bounds(3) == 0) then
              phi(nx,0,iz) = boundsavey(0,2,iz)
              phi(nx,ny,iz) = boundsavey(ny,2,iz)
            endif
          endif
        enddo
      endif

c     --- Apply altered difference equation to the points near the
c     --- surface of the conductor boundaries.
      if (lcndbndy) then
        if (parity == 0) then
         call condbndymg(conductors%evensubgrid,nx,ny,nz,phi,rho,
     &                   rdx2c,rdy2c,rdz2c,spm1,mgparam,bounds,
     &                   mglevel,mgform,icndbndy)
        endif
        if (parity == 1) then
         call condbndymg(conductors%oddsubgrid,nx,ny,nz,phi,rho,
     &                   rdx2c,rdy2c,rdz2c,spm1,mgparam,bounds,
     &                   mglevel,mgform,icndbndy)
        endif
      endif

c     --- Put desired potential onto conductors in phi array.
      call cond_potmg(conductors%interior,
     &                nx,ny,nz,phi,mglevel,mgform,.false.)

c     --- set phi in the z guard planes
      if (bounds(4) == 1) phi(:,:,izl-2) = phi(:,:,izl)
      if (bounds(5) == 1) phi(:,:,izr)   = phi(:,:,izr-2)
      if (bounds(4) == 2 .and. nz == nzfull) phi(:,:,izl-2) = phi(:,:,izr-2)
      if (bounds(5) == 2 .and. nz == nzfull) phi(:,:,nz:nz+1) = phi(:,:,0:1)

      return
      end
c=============================================================================
      subroutine mgsor_loop(i_parity,s_parity,izl,izr,izstep,
     &                      ii1,ii2,iistep,iimx,iipx,iimy,iipy,iimz,iipz,
     &                      phi1d,rho1d,nxy,rdx2c,rdy2c,rdz2c,spm1,linbends,
     &                      rstar,dx,rdx2,rdy2,rdz2,rdel,const,xx)
      real(kind=8):: phi1d(*),rho1d(*),rstar(*),xx(*)
      integer(ISZ):: i_parity,s_parity,izl,izr,izstep,nxy,ii1,ii2,iistep
      integer(ISZ):: iimx,iipx,iimy,iipy,iimz,iipz
      real(kind=8):: rdx2c,rdy2c,rdz2c,spm1,dx,rdx2,rdy2,rdz2,rdel,const
      logical(ISZ):: linbends

c This routine provides the loops over the phi array which solves the 
c iterative equation.  The loops are set up to do an outer loop over
c the axial dimension (the third dimension, or z) and an inner loop
c over the transverse dimensions.  By passing the proper loop limits when
c possible, the entire array can be looped over in the inner loop, increasing
c performance.  The structure of having two seperate loops makes this routine
c usable in the calculation of boundary conditions on the grid edges.  The
c structure is also amenable to the bent field solver.
c
c The variable parity is used to ensure that the loops maintain the parity as
c passed by i_parity.  The value of parity is not necessarily the same as
c i_parity.  When parity is calculated, the expression (iz-1) is used since
c phi1d is passed starting at iz=-1 and the array location phi(0,0,0) is
c considered even.  The expression (ii1-1) is used since the 1d array starts
c at an index of one, but the first element is even.
c
c The index of rho1d is (ii-nxy) since rho has dimensions rho(,,0:nz), whereas
c the value of ii is based on an array size of (,,-1:nz+1).
c
c Inputs:
c    i_parity  :input parity
c    s_parity  :parity multiplier
c    izl       :left, or starting z index
c    izr       :right, or ending z index
c    izstep    :step for z loop
c    ii1       :starting index for inner loop
c    ii2       :ending index for inner loop
c    iistep    :step for inner loop
c    iimx      :change in 1d index to array element phi(ix-1,iy  ,iz  )
c    iipx      :change in 1d index to array element phi(ix+1,iy  ,iz  )
c    iimy      :change in 1d index to array element phi(ix  ,iy-1,iz  )
c    iipy      :change in 1d index to array element phi(ix  ,iy+1,iz  )
c    iimz      :change in 1d index to array element phi(ix  ,iy  ,iz-1)
c    iipz      :change in 1d index to array element phi(ix  ,iy  ,iz+1)
c    xx        :value of x used with bends (see scrtch in psor3d)
c The rest of the inputs are the same as in psor3d.
c
c The if test (iistep == 2) is needed for the workstation version.  Without
c the loop stride explicitly given, there is about a 50% performance reduction
c on HIF, an HP 735 workstation.

      integer(ISZ):: iz,ii,parity
      real(kind=8):: x,rs,r,xfact,del2bndc,dh,dxi

      if (.not. linbends) then 
        if (izl == izr) then
          iz=izl
          parity = mod(i_parity*s_parity + (iz-1)*nxy + ii1-1,2)
          if (iistep == 2) then
!$OMP DO
            do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, 2
              phi1d(ii) = rho1d(ii-nxy) +
     &                  (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                  (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                  (phi1d(ii+iimz) + phi1d(ii+iipz))*rdz2c +
     &                  spm1*phi1d(ii)
            enddo
!$OMP END DO
          else
!$OMP DO
            do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
              phi1d(ii) = rho1d(ii-nxy) +
     &                  (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                  (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                  (phi1d(ii+iimz) + phi1d(ii+iipz))*rdz2c +
     &                  spm1*phi1d(ii)
            enddo
!$OMP END DO
          endif
        else
!$OMP DO
          do iz=izl,izr,izstep
            parity = mod(i_parity*s_parity + (iz-1)*nxy + ii1-1,2)
            if (iistep == 2) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, 2
                phi1d(ii) = rho1d(ii-nxy) +
     &                    (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                    (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                    (phi1d(ii+iimz) + phi1d(ii+iipz))*rdz2c +
     &                    spm1*phi1d(ii)
              enddo
            else
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                phi1d(ii) = rho1d(ii-nxy) +
     &                    (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                    (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                    (phi1d(ii+iimz) + phi1d(ii+iipz))*rdz2c +
     &                    spm1*phi1d(ii)
              enddo
            endif
          enddo
!$OMP END DO
        endif
      else
c       --- bends, loop over z slices.  Comment: timing tests indicate that 
c       --- use of a 1d array is faster than a 3d array.   
        dxi = 1./dx
!$OMP DO
        do iz=izl,izr,izstep
          parity = mod(i_parity*s_parity + (iz-1)*nxy + ii1-1,2)
          rs  = rstar(iz+1)
          if ( abs(rstar(iz  )) > LARGEPOS*1.e-6 .and. 
     &         abs(rstar(iz+2)) > LARGEPOS*1.e-6 ) then 
c           --- current z-slice not in bend or bordering a bend, do not 
c           --- include curvature terms from Poisson's equation    
            do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
              phi1d(ii) = rho1d(ii-nxy) +
     &            (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2c +
     &            (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2c +
     &            (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2c +
     &            spm1*phi1d(ii)
            enddo  
          else  
c           --- current z-slice in a bend and near border of a bend, include 
c           --- curvature terms from Poisson's equation.  Note: if bordering 
c           --- a bend, the dh/dz "jump term" should be included.  The 
c           --- current method of treating this jump term may be inaccurate. 
c           --- r is set to 1 over (rs+x) to reduce the number of divides.
            dh = 1./rstar(iz+2) - 1./rstar(iz)
            if (abs(dh) > SMALLPOS .and. abs(rs) < LARGEPOS*1.e-6) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                x = xx(ii-iz*nxy)
                r = 1./(rs + x)
                xfact    = (x*r)*((x*r)-2.) 
                del2bndc = const/( 1.+xfact*rdel )
                phi1d(ii) = rho1d(ii-nxy) + del2bndc*( 
     &              (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2 +
     &              (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2 +
     &              (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2+
     &              .5*(phi1d(ii+iipx)-phi1d(ii+iimx))*r*dxi + 
     &              xfact*(phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2 - 
     &              .25*(rs*r)**3*x*(phi1d(ii+iipz) - phi1d(ii+iimz))*dh*rdz2) +
     &              spm1*phi1d(ii)
              enddo 

c           --- current z-slice is near border of a bend only
            elseif (abs(dh) > SMALLPOS) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                x = xx(ii-iz*nxy)
                phi1d(ii) = rho1d(ii-nxy) + const*(
     &              (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2 +
     &              (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2 +
     &              (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2 -
     &              .25*x*(phi1d(ii+iipz) - phi1d(ii+iimz))*dh*rdz2) +
     &              spm1*phi1d(ii)
              enddo

c           --- current z-slice is in a bend only
            elseif (abs(rs) < LARGEPOS*1.e-6) then
              do ii = parity*s_parity + iz*nxy + ii1, iz*nxy + ii2, iistep
                x = xx(ii-iz*nxy)
                r = 1./(rs + x)
                xfact    = (x*r)*((x*r)-2.)
                del2bndc = const/( 1.+xfact*rdel )
                phi1d(ii) = rho1d(ii-nxy) + del2bndc*(
     &              (phi1d(ii+iipx) + phi1d(ii+iimx))*rdx2 +
     &              (phi1d(ii+iipy) + phi1d(ii+iimy))*rdy2 +
     &              (phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2+
     &              .5*(phi1d(ii+iipx) - phi1d(ii+iimx))*r*dxi +
     &              xfact*(phi1d(ii+iipz) + phi1d(ii+iimz))*rdz2) +
     &              spm1*phi1d(ii)
              enddo
            endif
          endif 
        enddo 
!$OMP END DO
      endif 

      return
      end
c=============================================================================
      subroutine cond_potmg(interior,nx,ny,nz,phi,mglevel,
     &                      mgform,mgform2init)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,ny,nz,mglevel,mgform
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      logical(ISZ):: mgform2init

c Set conductor points to the desired potential.

      integer(ISZ):: ic,ix,iy,iz

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

      return
      end
c=============================================================================
      subroutine cond_potmgres(interior,nx,ny,nz,res,mglevel,
     &                         mgform,mgform2init)
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,ny,nz,mglevel,mgform
      real(kind=8):: res(-1:nx+1,-1:ny+1,-3:nz+3)
      logical(ISZ):: mgform2init

c Set conductor points to the desired potential.

      integer(ISZ):: ic,ix,iy,iz

!$OMP DO
      do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        res(ix,iy,iz) = 0.
      enddo
!$OMP END DO

      return
      end
c=============================================================================
      subroutine condbndymg(subgrid,nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,srp,
     &                      bounds,mglevel,mgform,icndbndy)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), rho(0:nx,0:ny,0:nz)
      real(kind=8):: rdx2c,rdy2c,rdz2c,spm1,srp
      integer(ISZ):: bounds(0:5),mgform,icndbndy

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c NOTE that rdx2cos and rdy2cos are rdx2c and rdy2c over mgparam.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: rdx2cos,rdy2cos,rdz2cos,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: voltfac
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz,il
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      rdx2cos = rdx2c/srp
      rdy2cos = rdy2c/srp
      rdz2cos = rdz2c/srp
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
        pxm = phi(ixm1,iy   ,iz   )
        pxp = phi(ixp1,iy   ,iz   )
        pym = phi(ix   ,iym1,iz   )
        pyp = phi(ix   ,iyp1,iz   )
        pzm = phi(ix   ,iy   ,izm1)
        pzp = phi(ix   ,iy   ,izp1)
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
        if (dels(0,ic) < 1.) then
          pxm = voltfac*volt(0,ic)/dels(0,ic)
          denom = denom + (1.-dels(0,ic))/dels(0,ic)*rdx2cos
          dosubgrid = .true.
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          pxp = voltfac*volt(1,ic)/dels(1,ic)
          denom = denom + (1.-dels(1,ic))/dels(1,ic)*rdx2cos
          dosubgrid = .true.
        endif
c       --- the point lower in y is inside the conductor
        if (dels(2,ic) < 1.) then
          pym = voltfac*volt(2,ic)/dels(2,ic)
          denom = denom + (1.-dels(2,ic))/dels(2,ic)*rdy2cos
          dosubgrid = .true.
        endif
c       --- the point higher in y is inside the conductor
        if (dels(3,ic) < 1.) then
          pyp = voltfac*volt(3,ic)/dels(3,ic)
          denom = denom + (1.-dels(3,ic))/dels(3,ic)*rdy2cos
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (dels(4,ic) < 1.) then
          pzm = voltfac*volt(4,ic)/dels(4,ic)
          denom = denom + (1.-dels(4,ic))/dels(4,ic)*rdz2cos
          dosubgrid = .true.
        endif
c       --- the point higher in z is inside the conductor
        if (dels(5,ic) < 1.) then
          pzp = voltfac*volt(5,ic)/dels(5,ic)
          denom = denom + (1.-dels(5,ic))/dels(5,ic)*rdz2cos
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          phi(ix,iy,iz) = (rho(ix,iy,iz) +
     &      (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c + (pzm+pzp)*rdz2c)/denom +
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
        pxm = phi(ixm1,iy   ,iz   )
        pxp = phi(ixp1,iy   ,iz   )
        pym = phi(ix   ,iym1,iz   )
        pyp = phi(ix   ,iyp1,iz   )
        pzm = phi(ix   ,iy   ,izm1)
        pzp = phi(ix   ,iy   ,izp1)
        denom = 1.
        dosubgrid = .false.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0 .and. mgform == 1) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- If both points in x are inside the conductor, both terms change
        if (dels(0,ic) < 1. .and. dels(1,ic) < 1.) then
          pxm = voltfac*volt(0,ic)*2./(dels(0,ic)*(dels(1,ic)+dels(0,ic)))
          pxp = voltfac*volt(1,ic)*2./(dels(1,ic)*(dels(1,ic)+dels(0,ic)))
          denom = denom +
     &          2.*(1.-dels(1,ic)*dels(0,ic))/(dels(1,ic)*dels(0,ic))*rdx2cos
          dosubgrid = .true.
        else
c         --- the point lower in x is inside the conductor
          if (dels(0,ic) < 1.) then
            pxm = voltfac*volt(0,ic)*2./(dels(0,ic)*(1.+dels(0,ic)))
            pxp = pxp*(2./(1.+dels(0,ic)))
            denom = denom + 2.*(1.-dels(0,ic))/dels(0,ic)*rdx2cos
            dosubgrid = .true.
          endif
c         --- the point higher in x is inside the conductor
          if (dels(1,ic) < 1.) then
            pxp = voltfac*volt(1,ic)*2./(dels(1,ic)*(1.+dels(1,ic)))
            pxm = pxm*(2./(1.+dels(1,ic)))
            denom = denom + 2.*(1.-dels(1,ic))/dels(1,ic)*rdx2cos
            dosubgrid = .true.
          endif
        endif
c       --- If both points in y are inside the conductor, both terms change
        if (dels(2,ic) < 1. .and. dels(3,ic) < 1.) then
          pym = voltfac*volt(2,ic)*2./(dels(2,ic)*(dels(3,ic)+dels(2,ic)))
          pyp = voltfac*volt(3,ic)*2./(dels(3,ic)*(dels(3,ic)+dels(2,ic)))
          denom = denom +
     &          2.*(1.-dels(3,ic)*dels(2,ic))/(dels(3,ic)*dels(2,ic))*rdy2cos
          dosubgrid = .true.
        else
c         --- the point lower in y is inside the conductor
          if (dels(2,ic) < 1.) then
            pym = voltfac*volt(2,ic)*2./(dels(2,ic)*(1.+dels(2,ic)))
            pyp = pyp*(2./(1.+dels(2,ic)))
            denom = denom + 2.*(1.-dels(2,ic))/dels(2,ic)*rdy2cos
            dosubgrid = .true.
          endif
c         --- the point higher in y is inside the conductor
          if (dels(3,ic) < 1.) then
            pyp = voltfac*volt(3,ic)*2./(dels(3,ic)*(1.+dels(3,ic)))
            pym = pym*(2./(1.+dels(3,ic)))
            denom = denom + 2.*(1.-dels(3,ic))/dels(3,ic)*rdy2cos
            dosubgrid = .true.
          endif
        endif
c       --- If both points in z are inside the conductor, both terms change
        if (dels(4,ic) < 1. .and. dels(5,ic) < 1.) then
          pzm = voltfac*volt(4,ic)*2./(dels(4,ic)*(dels(5,ic)+dels(4,ic)))
          pzp = voltfac*volt(5,ic)*2./(dels(5,ic)*(dels(5,ic)+dels(4,ic)))
          denom = denom +
     &          2.*(1.-dels(5,ic)*dels(4,ic))/(dels(5,ic)*dels(4,ic))*rdz2cos
          dosubgrid = .true.
        else
c         --- the point lower in z is inside the conductor
          if (dels(4,ic) < 1.) then
            pzm = voltfac*volt(4,ic)*2./(dels(4,ic)*(1.+dels(4,ic)))
            pzp = pzp*(2./(1.+dels(4,ic)))
            denom = denom + 2.*(1.-dels(4,ic))/dels(4,ic)*rdz2cos
            dosubgrid = .true.
          endif
c         --- the point higher in z is inside the conductor
          if (dels(5,ic) < 1.) then
            pzp = voltfac*volt(5,ic)*2./(dels(5,ic)*(1.+dels(5,ic)))
            pzm = pzm*(2./(1.+dels(5,ic)))
            denom = denom + 2.*(1.-dels(5,ic))/dels(5,ic)*rdz2cos
            dosubgrid = .true.
          endif
        endif

c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          phi(ix,iy,iz) = (rho(ix,iy,iz) +
     &      (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c + (pzm+pzp)*rdz2c)/denom +
     &      spm1*subgrid%prevphi(ic)
        endif

      enddo
!$OMP END DO

      endif

      return
      end
c=============================================================================
      subroutine condbndyres(subgrid,nx,ny,nz,phi,rho,res,rdx2,rdy2,rdz2,
     &                       mgparam,bounds,mglevel,
     &                       mgform,mgform2init,icndbndy)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: res(-1:nx+1,-1:ny+1,-3:nz+3)
      real(kind=8):: rdx2,rdy2,rdz2,mgparam
      integer(ISZ):: bounds(0:5),cstart(0:100),mgform,icndbndy
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

      real(kind=8):: const,rdx2c,rdy2c,rdz2c,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: voltfac
      real(kind=8):: rdx2cs,rdy2cs,rdz2cs,ppp
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1,izp1,izm1
      integer(ISZ):: ix,iy,iz,il
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      const = 0.5/(rdx2+rdy2+rdz2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdz2c = rdz2*const
      rdx2cs = mgparam*rdx2*const
      rdy2cs = mgparam*rdy2*const
      rdz2cs = mgparam*rdz2*const
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

        if (iz == 0  .and. bounds(4) == 0) cycle
        if (iz == nz .and. bounds(5) == 0) cycle

c       --- set temporaries with initial values
        pxm = phi(ixm1,iy   ,iz   )
        pxp = phi(ixp1,iy   ,iz   )
        pym = phi(ix   ,iym1,iz   )
        pyp = phi(ix   ,iyp1,iz   )
        pzm = phi(ix   ,iy   ,izm1)
        pzp = phi(ix   ,iy   ,izp1)
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
        if (dels(0,ic) < 1.) then
          pxm = voltfac*volt(0,ic)/dels(0,ic)
          denom = denom + (1.-dels(0,ic))/dels(0,ic)*rdx2c
          ppp = min(ppp,dels(0,ic))
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          pxp = voltfac*volt(1,ic)/dels(1,ic)
          denom = denom + (1.-dels(1,ic))/dels(1,ic)*rdx2c
          ppp = min(ppp,dels(1,ic))
        endif
c       --- the point lower in y is inside the conductor
        if (dels(2,ic) < 1.) then
          pym = voltfac*volt(2,ic)/dels(2,ic)
          denom = denom + (1.-dels(2,ic))/dels(2,ic)*rdy2c
          ppp = min(ppp,dels(2,ic))
        endif
c       --- the point higher in y is inside the conductor
        if (dels(3,ic) < 1.) then
          pyp = voltfac*volt(3,ic)/dels(3,ic)
          denom = denom + (1.-dels(3,ic))/dels(3,ic)*rdy2c
          ppp = min(ppp,dels(3,ic))
        endif
c       --- the point lower in z is inside the conductor
        if (dels(4,ic) < 1.) then
          pzm = voltfac*volt(4,ic)/dels(4,ic)
          denom = denom + (1.-dels(4,ic))/dels(4,ic)*rdz2c
          ppp = min(ppp,dels(4,ic))
        endif
c       --- the point higher in z is inside the conductor
        if (dels(5,ic) < 1.) then
          pzp = voltfac*volt(5,ic)/dels(5,ic)
          denom = denom + (1.-dels(5,ic))/dels(5,ic)*rdz2c
          ppp = min(ppp,dels(5,ic))
        endif
c       --- calculate the residual based on the boundary conditions
        if (ppp < 1.) then
          if (mgform2init) ppp = 1.
#ifdef WITHCHOMBO
          ppp = 1.
#endif
          res(ix,iy,iz) = ppp*(rho(ix,iy,iz)
     &           + (pxm+pxp)*rdx2cs + (pym+pyp)*rdy2cs + (pzm+pzp)*rdz2cs
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

        if (iz == 0  .and. bounds(4) == 0) cycle
        if (iz == nz .and. bounds(5) == 0) cycle

c       --- set temporaries with initial values
        pxm = phi(ixm1,iy   ,iz   )
        pxp = phi(ixp1,iy   ,iz   )
        pym = phi(ix   ,iym1,iz   )
        pyp = phi(ix   ,iyp1,iz   )
        pzm = phi(ix   ,iy   ,izm1)
        pzp = phi(ix   ,iy   ,izp1)
        denom = 1.
        ppp = 1.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0 .and. (mgform == 1 .or. mgform2init)) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- If both points in x are inside the conductor, change both sides
        if (dels(0,ic) < 1. .and. dels(1,ic) < 1.) then
          pxm = voltfac*volt(0,ic)*2./(dels(0,ic)*(dels(1,ic)+dels(0,ic)))
          pxp = voltfac*volt(1,ic)*2./(dels(1,ic)*(dels(1,ic)+dels(0,ic)))
          denom = denom +
     &          2.*(1.-dels(1,ic)*dels(0,ic))/(dels(1,ic)*dels(0,ic))*rdx2c
          ppp = min(ppp,dels(0,ic))
          ppp = min(ppp,dels(1,ic))
        else
c         --- the point lower in x is inside the conductor
          if (dels(0,ic) < 1.) then
            pxm = voltfac*volt(0,ic)*2./(dels(0,ic)*(1.+dels(0,ic)))
            pxp = pxp*(2./(1.+dels(0,ic)))
            denom = denom + 2.*(1.-dels(0,ic))/dels(0,ic)*rdx2c
            ppp = min(ppp,dels(0,ic))
          endif
c         --- the point higher in x is inside the conductor
          if (dels(1,ic) < 1.) then
            pxp = voltfac*volt(1,ic)*2./(dels(1,ic)*(1.+dels(1,ic)))
            pxm = pxm*(2./(1.+dels(1,ic)))
            denom = denom + 2.*(1.-dels(1,ic))/dels(1,ic)*rdx2c
            ppp = min(ppp,dels(1,ic))
          endif
        endif
c       --- If both points in y are inside the conductor, change both sides
        if (dels(2,ic) < 1. .and. dels(3,ic) < 1.) then
          pym = voltfac*volt(2,ic)*2./(dels(2,ic)*(dels(3,ic)+dels(2,ic)))
          pyp = voltfac*volt(3,ic)*2./(dels(3,ic)*(dels(3,ic)+dels(2,ic)))
          denom = denom +
     &          2.*(1.-dels(3,ic)*dels(2,ic))/(dels(3,ic)*dels(2,ic))*rdy2c
          ppp = min(ppp,dels(2,ic))
          ppp = min(ppp,dels(3,ic))
        else
c         --- the point lower in y is inside the conductor
          if (dels(2,ic) < 1.) then
            pym = voltfac*volt(2,ic)*2./(dels(2,ic)*(1.+dels(2,ic)))
            pyp = pyp*(2./(1.+dels(2,ic)))
            denom = denom + 2.*(1.-dels(2,ic))/dels(2,ic)*rdy2c
            ppp = min(ppp,dels(2,ic))
          endif
c         --- the point higher in y is inside the conductor
          if (dels(3,ic) < 1.) then
            pyp = voltfac*volt(3,ic)*2./(dels(3,ic)*(1.+dels(3,ic)))
            pym = pym*(2./(1.+dels(3,ic)))
            denom = denom + 2.*(1.-dels(3,ic))/dels(3,ic)*rdy2c
            ppp = min(ppp,dels(3,ic))
          endif
        endif
c       --- If both points in z are inside the conductor, change both sides
        if (dels(4,ic) < 1. .and. dels(5,ic) < 1.) then
          pzm = voltfac*volt(4,ic)*2./(dels(4,ic)*(dels(5,ic)+dels(4,ic)))
          pzp = voltfac*volt(5,ic)*2./(dels(5,ic)*(dels(5,ic)+dels(4,ic)))
          denom = denom +
     &          2.*(1.-dels(5,ic)*dels(4,ic))/(dels(5,ic)*dels(4,ic))*rdz2c
          ppp = min(ppp,dels(4,ic))
          ppp = min(ppp,dels(5,ic))
        else
c         --- the point lower in z is inside the conductor
          if (dels(4,ic) < 1.) then
            pzm = voltfac*volt(4,ic)*2./(dels(4,ic)*(1.+dels(4,ic)))
            pzp = pzp*(2./(1.+dels(4,ic)))
            denom = denom + 2.*(1.-dels(4,ic))/dels(4,ic)*rdz2c
            ppp = min(ppp,dels(4,ic))
          endif
c         --- the point higher in z is inside the conductor
          if (dels(5,ic) < 1.) then
            pzp = voltfac*volt(5,ic)*2./(dels(5,ic)*(1.+dels(5,ic)))
            pzm = pzm*(2./(1.+dels(5,ic)))
            denom = denom + 2.*(1.-dels(5,ic))/dels(5,ic)*rdz2c
            ppp = min(ppp,dels(5,ic))
          endif
        endif

c       --- calculate the residual based on the boundary conditions
        if (ppp < 1.) then
          if (mgform2init) ppp = 1.
#ifdef WITHCHOMBO
          ppp = 1.
#endif
          res(ix,iy,iz) = ppp*(rho(ix,iy,iz)
     &           + (pxm+pxp)*rdx2cs + (pym+pyp)*rdy2cs + (pzm+pzp)*rdz2cs
     &           - phi(ix,iy,iz)*mgparam*denom)
        endif
      enddo
!$OMP END DO

      endif

      return
      end
c=============================================================================
      subroutine residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
     &                    mglevel,bounds,mgparam,mgform,mgform2init,
     &                    lcndbndy,icndbndy,conductors)
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      real(kind=8):: res(-1:nx+1,-1:ny+1,-3:nz+3)
      integer(ISZ):: mglevel,bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgform
      logical(ISZ):: mgform2init
      logical(ISZ):: lcndbndy
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
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1
      real(kind=8):: const,dxsqic,dysqic,dzsqic
      const = 0.5/(dxsqi+dysqi+dzsqi)
      dxsqic = dxsqi*mgparam*const
      dysqic = dysqi*mgparam*const
      dzsqic = dzsqi*mgparam*const

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

            res(ix,iy,iz) = rho(ix,iy,iz)
     &          +  (phi(ixm1,iy  ,iz  )+phi(ixp1,iy  ,iz  ))*dxsqic
     &          +  (phi(ix  ,iym1,iz  )+phi(ix  ,iyp1,iz  ))*dysqic
     &          +  (phi(ix  ,iy  ,izm1)+phi(ix  ,iy  ,izp1))*dzsqic
     &          -  phi(ix,iy,iz)*mgparam

          enddo
        enddo
      enddo
!$OMP END DO

c     --- Zero the residual inside conductors.
      call cond_potmgres(conductors%interior,
     &                nx,ny,nz,res,mglevel,mgform,.false.)

      if (lcndbndy) then
c       --- Calculate the residual near the conductor.
        call condbndyres(conductors%evensubgrid,nx,ny,nz,phi,rho,res,
     &                   dxsqi,dysqi,dzsqi,mgparam,bounds,
     &                   mglevel,mgform,mgform2init,icndbndy)
        call condbndyres(conductors%oddsubgrid,nx,ny,nz,phi,rho,res,
     &                   dxsqi,dysqi,dzsqi,mgparam,bounds,
     &                   mglevel,mgform,mgform2init,icndbndy)
      endif

c     --- Transverse boundaries
      if (bounds(0) == 0) res(-1:0,:,:) = 0.
      if (bounds(1) == 0) res(nx:nx+1,:,:) = 0.
      if (bounds(0) == 1) res(-1,:,:) = res(1,:,:)
      if (bounds(1) == 1) res(nx+1,:,:) = res(nx-1,:,:)
      if (bounds(0) == 2) res(-1,:,:) = res(nx-1,:,:)
      if (bounds(1) == 2) res(nx:nx+1,:,:) = res(0:1,:,:)

      if (bounds(2) == 0) res(:,-1:0,:) = 0.
      if (bounds(3) == 0) res(:,ny:ny+1,:) = 0.
      if (bounds(2) == 1) res(:,-1,:) = res(:,1,:)
      if (bounds(3) == 1) res(:,ny+1,:) = res(:,ny-1,:)
      if (bounds(2) == 2) res(:,-1,:) = res(:,ny-1,:)
      if (bounds(3) == 2) res(:,ny:ny+1,:) = res(:,0:1,:)

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
c=============================================================================
c=============================================================================
      subroutine getmglevels(nx,ny,nz,nzfull,dx,dy,dz,conductors,
     &                       my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c     --- Only do this calculation if it is needed. This
c     --- is for parallel version. There, this is a global operation and
c     --- if it is not needed, it would be an unneccesary synchronization
c     --- point.  If any of the levels is greater than one, then this
c     --- operation has already been done for this set of conductors, so
c     --- skip it.
      if (maxval(conductors%levellx) > 1 .or.
     &    maxval(conductors%levelly) > 1 .or.
     &    maxval(conductors%levellz) > 1) return

      call getmglevelsrecur(nx,ny,nz,nzfull,dx,dy,dz,0,1.,1.,1.,
     &                      conductors,
     &                      my_index,nslaves,izfsslave,nzfsslave)

      return
      end
c=============================================================================
      RECURSIVE subroutine getmglevelsrecur(nx,ny,nz,nzfull,dx,dy,dz,
     &                                      mglevel,lx,ly,lz,
     &                                      conductors,
     &                                     my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      integer(ISZ):: mglevel
      real(kind=8):: lx,ly,lz
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves),nzfsslave(0:nslaves)

c Calculate the levels of coarsening, saving the nx, ny, and nz for each level.

      integer(ISZ):: izfsslavec(0:nslaves),nzfsslavec(0:nslaves)
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse,cxf,cyf,czf
      logical(ISZ):: partialcoarsening

#ifdef MPIPARALLEL
      conductors%leveliz(mglevel) = izfsslave(my_index)
      conductors%levelnz(mglevel) = nzfsslave(my_index)
#else
      conductors%leveliz(mglevel) = 0
      conductors%levelnz(mglevel) = nz
#endif
      conductors%levellx(mglevel) = lx
      conductors%levelly(mglevel) = ly
      conductors%levellz(mglevel) = lz
      conductors%levels = mglevel + 1

      if (nx < 4 .or. ny < 4 .or. nzfull < 4 .or. mglevel == 100) return

c     --- Calculate the size of the next coarsest grid. If the current
c     --- size is twice and odd integer, the next coarsest size will be
c     --- n/2+1, gauranteeing that it is even.
c     --- Note that nzcoarse is only used if full-coarsening is done.
      nxcoarse = nx/2
      if (mod(nxcoarse,2) == 1) nxcoarse = nxcoarse + 1
      nycoarse = ny/2
      if (mod(nycoarse,2) == 1) nycoarse = nycoarse + 1
      nzcoarse = nz/2
      if (mod(nzcoarse,2) == 1) nzcoarse = nzcoarse + 1
      nzfullcoarse = nzfull/2
      if (mod(nzfullcoarse,2) == 1) nzfullcoarse = nzfullcoarse + 1

      dxcoarse = dx*nx/nxcoarse
      dycoarse = dy*ny/nycoarse
      dzcoarse = dz*nzfull/nzfullcoarse
      cxf = dxcoarse/dx
      cyf = dycoarse/dy
      czf = dzcoarse/dz

c     --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
c     --- all axis.
      partialcoarsening = (dz > 4./3.*dx)

#ifdef MPIPARALLEL
c     --- This must be a global operation since, due to roundoff, each
c     --- processor can get a different value if dz == 4./3.*dx.
c     call parallellor(partialcoarsening)
#endif

      if (partialcoarsening) then
        call getmglevelsrecur(nxcoarse,nycoarse,nz,nzfull,dxcoarse,dycoarse,dz,
     &                        mglevel+1,lx*cxf,ly*cyf,lz,conductors,
     &                        my_index,nslaves,izfsslave,nzfsslave)
      else

#ifdef MPIPARALLEL
c       --- Find domains in coarser grid
        call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                  nzfull,nzfullcoarse)
c       --- Set new value of nz
        nzcoarse = nzfsslavec(my_index)
#endif

        call getmglevelsrecur(nxcoarse,nycoarse,nzcoarse,nzfullcoarse,
     &                        dxcoarse,dycoarse,dzcoarse,
     &                        mglevel+1,lx*cxf,ly*cyf,lz*czf,conductors,
     &                        my_index,nslaves,izfsslavec,nzfsslavec)

      endif

      return
      end
c=============================================================================
c=============================================================================
c     --- These were only needed to get around a compiler bug that appears
c     --- to be fixed. They are left in incase they are in fact still needed.
c=============================================================================
      integer(ISZ) function getthestupidnumber(conductors)
      use ConductorTypemodule
      type(ConductorType):: conductors
      getthestupidnumber = conductors%interior%n
      return
      end
c=============================================================================
      integer(ISZ) function getthestupidevennumber(conductors)
      use ConductorTypemodule
      type(ConductorType):: conductors
      getthestupidevennumber = conductors%evensubgrid%n
      return
      end
c=============================================================================
      integer(ISZ) function getthestupidoddnumber(conductors)
      use ConductorTypemodule
      type(ConductorType):: conductors
      getthestupidoddnumber = conductors%oddsubgrid%n
      return
      end
c=============================================================================
      subroutine setcndtr3dmg(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                        bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                        conductors)
      use ConductorTypemodule
      use Parallel
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: bound0,boundnz,boundxy
      logical(ISZ):: l2symtry,l4symtry
      type(ConductorType):: conductors

c Makes call to the recursive routine which generates the conductors.

      integer(ISZ):: nzfull
      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: zmmin1,dx1,dy1,dz1
      integer(ISZ):: nzlocal,nx1,ny1,nz1,nzfull1
      integer(ISZ):: globalb0,globalbnz
      integer(ISZ):: localb0,localbnz
      integer(ISZ):: iiii
c     integer(ISZ):: getthestupidnumber
c     integer(ISZ):: getthestupidevennumber,getthestupidoddnumber

c     --- nzfull and zmminglobal are needed, but difficult to pass down from
c     --- above. But, they can be obtained from local info
#ifdef MPIPARALLEL
      nzfull = izfsslave(nslaves-1)+nzfsslave(nslaves-1)
      zmminglobal = zmslmin(0)
      nzlocal = nzfsslave(my_index)
#else
      nzfull = nz
      zmminglobal = zmmin
      nzlocal = nz
#endif
      globalb0 = bound0
      globalbnz = boundnz

c     --- Calculate levels of coarsening
      call getmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave)

      do il=0,conductors%levels-1

c       --- Set the grid parameters for this level
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        dz1 = dz*conductors%levellz(il)
        nx1 = nint(nx/conductors%levellx(il))
        ny1 = nint(ny/conductors%levelly(il))
        nz1 = conductors%levelnz(il)
        nzfull1 = nzfull/conductors%levellx(il)
        zmmin1 = zmminglobal + dz1*conductors%leveliz(il)
        localb0 = globalb0
        localbnz = globalbnz
#ifdef MPIPARALLEL
        if (conductors%leveliz(il) > 0) localb0 = -1
        if (conductors%leveliz(il)+nz1 < nzfull1) localbnz = -1
#endif

c       --- Get the conductors for this level
c       --- The stupidnumber things help get around a compiler bug.
        ic1 = conductors%interior%n + 1
c       iiii = getthestupidnumber(conductors)
        iiii = conductors%interior%n
        ie1 = conductors%evensubgrid%n + 1
c       iiii = getthestupidevennumber(conductors)
        iiii = conductors%evensubgrid%n
        io1 = conductors%oddsubgrid%n + 1
c       iiii = getthestupidoddnumber(conductors)
        iiii = conductors%oddsubgrid%n
        bound0 = localb0
        boundnz = localbnz
        call setcndtr3d(xmmin,ymmin,zmmin1,zbeam,zgrid,nx1,ny1,nz1,dx1,dy1,dz1,
     &                  bound0,boundnz,boundxy,l2symtry,l4symtry)
        ic2 = conductors%interior%n
        ie2 = conductors%evensubgrid%n
        io2 = conductors%oddsubgrid%n
        conductors%interior%ilevel(ic1:ic2) = il
        conductors%evensubgrid%ilevel(ie1:ie2) = il
        conductors%oddsubgrid%ilevel(io1:io2) = il

      enddo

      bound0 = globalb0
      boundnz = globalbnz

      return
      end
c=============================================================================
      subroutine srfrvout3dmg(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                        xmin,xmax,ymin,ymax,lshell,
     &                        zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                        ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,
     &                        condid,conductors)
      use ConductorTypemodule
      use Parallel
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmax
      logical(ISZ):: lfill,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: condid
      type(ConductorType):: conductors

c Makes calls to srfrvout3d for each level

      integer(ISZ):: nzfull
      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: zmmin1,zmmax1,dx1,dy1,dz1
      integer(ISZ):: nzlocal,nx1,ny1,nz1
      integer(ISZ):: ix_axis1,iy_axis1

c     --- nzfull and zmminglobal are needed, but difficult to pass down from
c     --- above. But, they can be obtained from local info
#ifdef MPIPARALLEL
      nzfull = izfsslave(nslaves-1)+nzfsslave(nslaves-1)
      zmminglobal = zmslmin(0)
      nzlocal = nzfsslave(my_index)
#else
      nzfull = nz
      zmminglobal = zmmin
      nzlocal = nz
#endif

c     --- Calculate levels of coarsening
      call getmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave)

      do il=0,conductors%levels-1

c       --- Set the grid parameters for this level
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        dz1 = dz*conductors%levellz(il)
        nx1 = nint(nx/conductors%levellx(il))
        ny1 = nint(ny/conductors%levelly(il))
        nz1 = conductors%levelnz(il)
        zmmin1 = zmminglobal + dz1*conductors%leveliz(il)
        zmmax1 = zmmin1 + dz1*nz1
        ix_axis1 = ix_axis/conductors%levellx(il)
        iy_axis1 = iy_axis/conductors%levelly(il)

c       --- Get the conductors for this level
        ic1 = conductors%interior%n + 1
        ie1 = conductors%evensubgrid%n + 1
        io1 = conductors%oddsubgrid%n + 1
        call srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                  xmin,xmax,ymin,ymax,lshell,
     &                  zmmin1,zmmax1,zbeam,dx1,dy1,dz1,nx1,ny1,nz1,
     &                  ix_axis1,iy_axis1,xmesh,ymesh,l2symtry,l4symtry,condid)
        ic2 = conductors%interior%n
        ie2 = conductors%evensubgrid%n
        io2 = conductors%oddsubgrid%n
        conductors%interior%ilevel(ic1:ic2) = il
        conductors%evensubgrid%ilevel(ie1:ie2) = il
        conductors%oddsubgrid%ilevel(io1:io2) = il

      enddo

      return
      end
c=============================================================================
      subroutine srfrvin3dmg(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                       xmin,xmax,ymin,ymax,lshell,
     &                       zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                       ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,
     &                       condid,conductors)
      use ConductorTypemodule
      use Parallel
      character(*):: rofzfunc
      real(kind=8):: volt,zmin,zmax,xcent,ycent,rmin
      logical(ISZ):: lfill,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: condid
      type(ConductorType):: conductors

c Makes call to the recursive routine which generates the conductors.

      integer(ISZ):: nzfull
      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: zmmin1,zmmax1,dx1,dy1,dz1
      integer(ISZ):: nzlocal,nx1,ny1,nz1
      integer(ISZ):: ix_axis1,iy_axis1

c     --- nzfull and zmminglobal are needed, but difficult to pass down from
c     --- above. But, they can be obtained from local info
#ifdef MPIPARALLEL
      nzfull = izfsslave(nslaves-1)+nzfsslave(nslaves-1)
      zmminglobal = zmslmin(0)
      nzlocal = nzfsslave(my_index)
#else
      nzfull = nz
      zmminglobal = zmmin
      nzlocal = nz
#endif

c     --- Calculate levels of coarsening
      call getmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave)

      do il=0,conductors%levels-1

c       --- Set the grid parameters for this level
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        dz1 = dz*conductors%levellz(il)
        nx1 = nint(nx/conductors%levellx(il))
        ny1 = nint(ny/conductors%levelly(il))
        nz1 = conductors%levelnz(il)
        zmmin1 = zmminglobal + dz1*conductors%leveliz(il)
        zmmax1 = zmmin1 + dz1*nz1
        ix_axis1 = ix_axis/conductors%levellx(il)
        iy_axis1 = iy_axis/conductors%levelly(il)

c       --- Get the conductors for this level
        ic1 = conductors%interior%n + 1
        ie1 = conductors%evensubgrid%n + 1
        io1 = conductors%oddsubgrid%n + 1
        call srfrvin3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                  xmin,xmax,ymin,ymax,lshell,
     &                  zmmin1,zmmax1,zbeam,dx1,dy1,dz1,nx1,ny1,nz1,
     &                  ix_axis1,iy_axis1,xmesh,ymesh,l2symtry,l4symtry,condid)
        ic2 = conductors%interior%n
        ie2 = conductors%evensubgrid%n
        io2 = conductors%oddsubgrid%n
        conductors%interior%ilevel(ic1:ic2) = il
        conductors%evensubgrid%ilevel(ie1:ie2) = il
        conductors%oddsubgrid%ilevel(io1:io2) = il

      enddo

      return
      end
c=============================================================================
      subroutine srfrvinout3dmg(rminofz,rmaxofz,volt,zmin,zmax,
     &                          xcent,ycent,lzend,xmin,xmax,ymin,ymax,lshell,
     &                          zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                          ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,
     &                          condid,conductors)
      use ConductorTypemodule
      use Parallel
      character(*):: rminofz,rmaxofz
      real(kind=8):: volt,zmin,zmax,xcent,ycent
      logical(ISZ):: lzend,lshell
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: zmmin,zmmax,zbeam,dx,dy,dz
      integer(ISZ):: nx,ny,nz,ix_axis,iy_axis
      real(kind=8):: xmesh(0:nx),ymesh(0:ny)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: condid
      type(ConductorType):: conductors

c Makes call to the recursive routine which generates the conductors.

      integer(ISZ):: nzfull
      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: zmmin1,zmmax1,dx1,dy1,dz1
      integer(ISZ):: nzlocal,nx1,ny1,nz1
      integer(ISZ):: ix_axis1,iy_axis1

c     --- nzfull and zmminglobal are needed, but difficult to pass down from
c     --- above. But, they can be obtained from local info
#ifdef MPIPARALLEL
      nzfull = izfsslave(nslaves-1)+nzfsslave(nslaves-1)
      zmminglobal = zmslmin(0)
      nzlocal = nzfsslave(my_index)
#else
      nzfull = nz
      zmminglobal = zmmin
      nzlocal = nz
#endif

c     --- Calculate levels of coarsening
      call getmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave)

      do il=0,conductors%levels-1

c       --- Set the grid parameters for this level
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        dz1 = dz*conductors%levellz(il)
        nx1 = nint(nx/conductors%levellx(il))
        ny1 = nint(ny/conductors%levelly(il))
        nz1 = conductors%levelnz(il)
        zmmin1 = zmminglobal + dz1*conductors%leveliz(il)
        zmmax1 = zmmin1 + dz1*nz1
        ix_axis1 = ix_axis/conductors%levellx(il)
        iy_axis1 = iy_axis/conductors%levelly(il)

c       --- Get the conductors for this level
        ic1 = conductors%interior%n + 1
        ie1 = conductors%evensubgrid%n + 1
        io1 = conductors%oddsubgrid%n + 1
        call srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,
     &                    xcent,ycent,lzend,xmin,xmax,ymin,ymax,lshell,
     &                    zmmin1,zmmax1,zbeam,dx1,dy1,dz1,nx1,ny1,nz1,
     &                    ix_axis1,iy_axis1,xmesh,ymesh,l2symtry,l4symtry,
     &                    condid)
        ic2 = conductors%interior%n
        ie2 = conductors%evensubgrid%n
        io2 = conductors%oddsubgrid%n
        conductors%interior%ilevel(ic1:ic2) = il
        conductors%evensubgrid%ilevel(ie1:ie2) = il
        conductors%oddsubgrid%ilevel(io1:io2) = il

      enddo

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      RECURSIVE subroutine checkconductors(nx,ny,nz,nzfull,dx,dy,dz,
     &                                     conductors,
     &                                     my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Recursively calls the routine to generate the conductor data at the
c various mesh resolutions needed by the MG solver.

      integer(ISZ):: ic,ix,iy,iz,ii,nn,na
      integer(ISZ),allocatable:: isort(:)

      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: dx1,dy1,dz1
      integer(ISZ):: nx1,ny1,nz1,levels
      integer(ISZ):: ix_axis1,iy_axis1

c     --- Make sure the coarsening levels are setup.
      call getmglevels(nx,ny,nz,nzfull,dx,dy,dz,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave)

c     --- Call the work routine for each level
      do il=0,conductors%levels-1
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        dz1 = dz*conductors%levellz(il)
        nx1 = nint(nx/conductors%levellx(il))
        ny1 = nint(ny/conductors%levelly(il))
        nz1 = conductors%levelnz(il)
        call checkconductors_work(conductors%interior,
     &                            conductors%evensubgrid,conductors%oddsubgrid,
     &                            nx1,ny1,nz1,dx1,dy1,dz1,il)
      enddo

      conductors%interior%istart = 1
      conductors%evensubgrid%istart = 1
      conductors%oddsubgrid%istart = 1
      levels = conductors%levels

c     --- Sort the conductor data by level number
c     --- First, sort conductor points
      if (conductors%interior%n > 0) then
        allocate(isort(conductors%interior%n))
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
        allocate(isort(conductors%evensubgrid%n))
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
        allocate(isort(conductors%oddsubgrid%n))
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
      integer(ISZ):: ic,itemp(ni,nn)
      itemp = icond
      do ic=1,nc
        icond(:,ic) = itemp(:,isort(ic))
      enddo
      return
      end
c=============================================================================
      subroutine rswapconductor(nc,isort,nn,ni,rcond)
      integer(ISZ):: nc,nn,ni,isort(nc)
      real(kind=8):: rcond(ni,nn)
c Utility function for sort conductors
      real(kind=8):: rtemp(ni,nn)
      integer(ISZ):: ic
      rtemp = rcond
      do ic=1,nc
        rcond(:,ic) = rtemp(:,isort(ic))
      enddo
      return
      end
c=============================================================================
      subroutine checkconductors_work(interior,evensubgrid,oddsubgrid,
     &                                nx,ny,nz,dx,dy,dz,mglevel)
      use ConductorInteriorTypemodule
      use ConductorSubGridTypemodule
      type(ConductorInteriorType):: interior
      type(ConductorSubGridType):: evensubgrid,oddsubgrid
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: dx,dy,dz

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
      integer(ISZ):: ic,i,ix,iy,iz,id

c     --- Set the conductor points.
      allocate(iii(0:nx,0:ny,0:nz))
      iii = 0
      do ic=1,interior%n
        if (interior%ilevel(ic) /= mglevel) cycle
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        iz = interior%indx(2,ic)
        if (ix < 0 .or. nx < ix .or.
     &      iy < 0 .or. ny < iy .or.
     &      iz < 0 .or. nz < iz) then
          interior%ilevel(ic) = -1
          cycle
        endif
        iii(ix,iy,iz) = interior%nmax + 1
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

        if (ix < 0 .or. nx < ix .or.
     &      iy < 0 .or. ny < iy .or.
     &      iz < 0 .or. nz < iz) then
          evensubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy,iz) == interior%nmax+1) then
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
          if (evensubgrid%dels(id,ic) < 1. .and.
     &        evensubgrid%dels(id,ic) < evensubgrid%dels(id,i)) then
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

        if (ix < 0 .or. nx < ix .or.
     &      iy < 0 .or. ny < iy .or.
     &      iz < 0 .or. nz < iz) then
          oddsubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy,iz) == interior%nmax+1) then
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
          if (oddsubgrid%dels(id,ic) < 1. .and.
     &        oddsubgrid%dels(id,ic) < oddsubgrid%dels(id,i)) then
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
      subroutine gatherrhofromchild(rho,nx,ny,nz,childrho,cnx,cny,cnz,
     &                              l,u,fulllower,childlower,childupper,
     &                              r,weights,owndomain,
     &                              dopbounds,pbounds,rootdims)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,cnx,cny,cnz
      integer(ISZ):: l(0:2),u(0:2),fulllower(0:2)
      integer(ISZ):: childlower(0:2),childupper(0:2)
      integer(ISZ):: r
      integer(ISZ):: dopbounds,pbounds(0:5),rootdims(0:2)
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      real(kind=8):: childrho(0:cnx,0:cny,0:cnz)
      real(kind=8):: weights(-r+1:r-1,-r+1:r-1,-r+1:r-1)
      integer(ISZ):: owndomain(0:cnx,0:cny,0:cnz)

      integer(ISZ):: ix,iy,iz,icx,icy,icz,iwx,iwy,iwz
      integer(ISZ):: ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2
      integer(ISZ):: iwx1,iwy1,iwz1,iwx2,iwy2,iwz2
      real(kind=8):: wrho
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      do iz=l(2),u(2)
        iz0 = iz - fulllower(2)
        iz1 = max(iz*r - r + 1,childlower(2))
        iz2 = min(iz*r + r - 1,childupper(2))
        iwz1 = iz1 - iz*r
        iwz2 = iz2 - iz*r
        icz = iz*r - childlower(2)
        do iy=l(1),u(1)
          iy0 = iy - fulllower(1)
          iy1 = max(iy*r - r + 1,childlower(1))
          iy2 = min(iy*r + r - 1,childupper(1))
          iwy1 = iy1 - iy*r
          iwy2 = iy2 - iy*r
          icy = iy*r - childlower(1)
          do ix=l(0),u(0)
            ix0 = ix - fulllower(0)
            ix1 = max(ix*r - r + 1,childlower(0))
            ix2 = min(ix*r + r - 1,childupper(0))
            iwx1 = ix1 - ix*r
            iwx2 = ix2 - ix*r
            icx = ix*r - childlower(0)

            do iwz=iwz1,iwz2
              do iwy=iwy1,iwy2
                do iwx=iwx1,iwx2
                  if (owndomain(icx+iwx,icy+iwy,icz+iwz) == 1) then
                    wrho=weights(iwx,iwy,iwz)*childrho(icx+iwx,icy+iwy,icz+iwz)

c                   --- Adjust for symmetries
                    if (dopbounds == 1) then
                      if (pbounds(0) == 1 .and. ix == 0 .and. iwx > 0)
     &                  wrho = 2*wrho
                      if (pbounds(1) == 1 .and. ix == rootdims(0) .and. iwx < 0)
     &                  wrho = 2*wrho
                      if (pbounds(2) == 1 .and. iy == 0 .and. iwy > 0)
     &                  wrho = 2*wrho
                      if (pbounds(3) == 1 .and. iy == rootdims(1) .and. iwy < 0)
     &                  wrho = 2*wrho
                      if (pbounds(4) == 1 .and. iz == 0 .and. iwz > 0)
     &                  wrho = 2*wrho
                      if (pbounds(5) == 1 .and. iz == rootdims(2) .and. iwz < 0)
     &                  wrho = 2*wrho
                    endif

                    rho(ix0,iy0,iz0) = rho(ix0,iy0,iz0) + wrho

                  endif

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      if (lf3dtimesubs) timegatherrhofromchild = timegatherrhofromchild +
     &                                           wtime() - substarttime

      return
      end
c=============================================================================
