#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID.M, version $Revision: 3.40 $, $Date: 2002/08/27 23:15:23 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid field sovler which is part of the F3D
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigrid3df(iwhich,nx,ny,nz,nzfull,dx,dy,dz,phi,rho1d,
     &                       rstar,linbend,l2symtry,l4symtry,scrtch,
     &                       xmmin,ymmin,zmmin,zbeam,zgrid)
      use Constant
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use Parallel
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho1d(*)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbend
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: scrtch(*)
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid

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
      real(kind=8):: maxerr
      real(kind=8):: dxsqi,dysqi,dzsqi,reps0c,rdel
      real(kind=8):: tdx,tdy,tdz
      integer(ISZ):: tnx,tny,tnz
      integer(ISZ):: i,ii,k,ix,iy,iz
      real(kind=8):: rs,x,r
      real(kind=8):: phisave(0:nx,0:ny,-1:nz+1)
c     --- The following only used when mgform == 2
      real(kind=8),allocatable:: rhosave(:),res(:)
      integer(ISZ):: localb0,localbnz
      integer(ISZ):: local2b0,local2bnz
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleft2(0:nslaves-1), izsendingleft2(0:nslaves-1)
      integer(ISZ):: whosendingright2(0:nslaves-1),izsendingright2(0:nslaves-1)
      integer(ISZ):: izfsslave2(0:nslaves-1),nzfsslave2(0:nslaves-1)
      character(72):: errline

c     --- Initialize temporaries
      nxy    = (nx+1)*(ny+1)
      nxyz   = (nx+1)*(ny+1)*(nz+1)
      dxsqi  = 1./dx**2
      dysqi  = 1./dy**2
      dzsqi  = 1./dz**2
      reps0c = mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
      rdel   = dzsqi/(dxsqi + dysqi + dzsqi)

c     --- Calculate levels of coarsening
      call setmglevels(nx,ny,nz,nzfull,dx,dy,dz)

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

c     --- Determine the points that make up the conductor.  This takes extra
c     --- time and so should not be done if the grid is not moving in the lab
c     --- frame.  Set gridmode to 1 to avoid this call. The data is then
c     --- converted and expanded for the multigrid solver.
      if (gridmode == 0 .or. iwhich == -2) then
        ncond = 0
        necndbdy = 0
        nocndbdy = 0
        call setcndtr(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                l2symtry,l4symtry)
      endif
      call checkconductors(nx,ny,nz,nzfull,dx,dy,dz,l2symtry,l4symtry)

c     --- Preset rho to increase performance (reducing the number of
c     --- multiplies in the main SOR sweep loop).
      if (.not. linbend) then
        rho1d(1:nxyz) = rho1d(1:nxyz)*reps0c
      else
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
        do ii=1,nxyz
          i = mod(ii-1,nx+1)
          k = (ii-1)/nxy
          rs = rstar(k)
          x  = xmmin + i*dx
          r  = rs + x
c         --- rearranged to reduce divides
c         --- rho1d(ii) = rho1d(ii)*(rs/r)*reps0c/
c         ---             ( 1. + (x/r)*((x/r)-2.)*rdel )
          rho1d(ii) = rho1d(ii)*reps0c*rs*r/(r*r + x*(x-2.*r)*rdel)
        enddo
c       --- Fill scratch array with x values so it can be looked up
c       --- in the bent beam loop instead of calculated.
        do ii = 1, nxy
          scrtch(ii) = xmmin + mod(ii-1,nx+1)*dx
        enddo
c       --- Change rstar if using Nuemann boundary conditions
#ifndef MPIPARALLEL
        if (bound0 == 1) rstar(-1) = rstar(1)
        if (boundnz == 1) rstar(nz+1) = rstar(nz-1)
#else
        if (bound0 == 1 .and. my_index == 0) rstar(-1) = rstar(1)
        if (boundnz == 1 .and. my_index == nslaves-1) rstar(nz+1) = rstar(nz-1)
#endif
      endif

c     --- If using residual correction form, need to save the original rho.
c     --- Also setup parallel arrays.
      if (mgform == 2) then
        allocate(rhosave(nxyz),res((nx+1)*(ny+1)*(nz+3)))
        rhosave = rho1d(1:nxyz)
#ifdef MPIPARALLEL
        call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                        bound0,boundnz,nzfull,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
        localb0 = bound0
        localbnz = boundnz
        if (izfsslave(my_index) > 0) localb0 = -1
        if (izfsslave(my_index)+nz < nzfull) localbnz = -1
#else
        localb0 = bound0
        localbnz = boundnz
#endif
      endif

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- very few iterations are done.
      mgiters = 0
      maxerr = 2.*mgtol + 1.
      do while (maxerr > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phi

c       --- If using residual correction form, calculate the residual and
c       --- copy it into rho_temp. Also, zero phi_temp (the initial error).
c       --- In the calls to cond_potmg and residual, the last argument
c       --- is true, telling the routines to use the actual value of
c       --- voltages rather than zero as is done otherwise for residual
c       --- correction form since it is operating on the error.
        if (mgform == 2) then
          call cond_potmg(nx,ny,nz,phisave,0,.false.,mgform,.true.)
          call residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phisave,rhosave,res,
     &                  0,localb0,localbnz,l2symtry,l4symtry,.true.,
     &                  my_index,nslaves,izfsslave,nzfsslave,
     &                  whosendingleft,izsendingleft,
     &                  whosendingright,izsendingright)
          rho1d(1:nxyz) = res(1+nxy:nxy+nxyz)
          phi = 0.
        endif

c       --- Do one vcycle.
        call vcycle(0,nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho1d,
     &              rstar,linbend,l2symtry,l4symtry,scrtch,
     &              boundxy,bound0,boundnz,my_index,nslaves,izfsslave,nzfsslave)

c       --- If using residual correction form, add the resulting error to phi.
        if (mgform == 2) phi = phi + phisave

c       --- When using residual correction form, the other planes do need
c       --- to be set when using other than Dirichlet boundaries since
c       --- those planes are only set with the error of phi.
        if (mgform == 2) then
          if (localb0  == 1) phi(:,:,-1) = phi(:,:,1)
          if (localbnz == 1) phi(:,:,nz+1) = phi(:,:,nz-1)
#ifndef MPIPARALLEL
          if (localb0  == 2) phi(:,:,-1) = phi(:,:,nz-1)
          if (localbnz == 2) phi(:,:,nz+1) = phi(:,:,1)
#else
          call mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,0,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
          call mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,-1,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
#endif
      endif

c       --- Calculate the change in phi.
        maxerr = 0.
        do iz=0,nz
          do iy=0,ny
            do ix=0,nx
              maxerr = max(maxerr,abs(phisave(ix,iy,iz) - phi(ix,iy,iz)))
            enddo
          enddo
        enddo
#ifdef MPIPARALLEL
c           --- calculate global sorerror
        call parallelmaxrealarray(maxerr,1)
#endif
c       print*,mgiters,maxerr

c       --- This line below seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop above.
c       maxerr = maxval(abs(phisave - phi))

      enddo

c     --- For Dirichlet boundary conditions, copy data into guard planes
c     --- For other boundary conditions, the guard planes are used during
c     --- the solve are so are already set.
      if (bound0 == 0) phi(:,:,-1) = phi(:,:,0)
      if (boundnz == 0) phi(:,:,nz+1) = phi(:,:,nz)

c     --- Make a print out.
      if (maxerr > mgtol) then
        call remark("Multigrid: Maximum number of iterations reached")
      endif
      write(errline,20) maxerr,mgiters
  20  format("Multigrid: Error converged to ",1pe11.3," in ",i4," v-cycles")
      call remark(errline)

c     --- If using residual correction form, restore saved rho
      if (mgform == 2) then
        rho1d(1:nxyz) = rhosave
        deallocate(rhosave,res)
      endif

c     --- Restore rho
      reps0c = 1./reps0c
      if (.not. linbend) then
        rho1d(1:nxyz) = rho1d(1:nxyz)*reps0c
      else
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
        do ii=1,nxyz
          i = mod(ii-1,nx+1)
          k = (ii-1)/nxy
          rs = rstar(k)
          x  = xmmin + i*dx
          r  = rs + x
          rho1d(ii) = rho1d(ii)/rs*reps0c*( r + x*((x/r)-2.)*rdel )
        enddo
      endif

      return
      end
c=============================================================================
      subroutine setmglevels(nx,ny,nz,nzfull,dx,dy,dz)
      use Parallel
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz

      call setmglevelsrecur(nx,ny,nz,nzfull,dx,dy,dz,0,1,1,1,
     &                      my_index,nslaves,izfsslave,nzfsslave)

      return
      end
c=============================================================================
      RECURSIVE subroutine setmglevelsrecur(nx,ny,nz,nzfull,dx,dy,dz,
     &                                      mglevel,lx,ly,lz,
     &                                     my_index,nslaves,izfsslave,nzfsslave)
      use Multigrid3d
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      integer(ISZ):: mglevel,lx,ly,lz
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves),nzfsslave(0:nslaves)

c Calculate the levels of coarsening, saving the nx, ny, and nz for each level.

      integer(ISZ):: izfsslave2(0:nslaves),nzfsslave2(0:nslaves)
      integer(ISZ):: nznew

      mglevelsnx(mglevel) = nx
      mglevelsny(mglevel) = ny
      mglevelsnzfull(mglevel) = nzfull
      mglevelsiz(mglevel) = izfsslave(my_index)
      mglevelsnz(mglevel) = nz
      mglevelslx(mglevel) = lx
      mglevelsly(mglevel) = ly
      mglevelslz(mglevel) = lz
      mglevels = mglevel + 1

      if (mod(nx,4) /= 0 .or.
     &    mod(ny,4) /= 0 .or.
     &    mod(nzfull,4) /= 0 .or.
     &    mglevel == mgmaxlevels) return

c     --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
c     --- all axis.
      if (dz > 4./3.*dx) then
        call setmglevelsrecur(nx/2,ny/2,nz,nzfull,dx*2,dy*2,dz,
     &                        mglevel+1,lx*2,ly*2,lz,
     &                        my_index,nslaves,izfsslave,nzfsslave)
      else

#ifdef MPIPARALLEL
c       --- Find domains in coarser grid
        call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslave2,nzfsslave2,
     &                  nzfull)
c       --- Set new value of nz
        nznew = nzfsslave2(my_index)
#else
        izfsslave2 = 0
        nznew = nz/2
#endif

        call setmglevelsrecur(nx/2,ny/2,nznew,nzfull/2,dx*2,dy*2,dz*2,
     &                        mglevel+1,lx*2,ly*2,lz*2,
     &                        my_index,nslaves,izfsslave2,nzfsslave2)
      endif

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycle(mglevel,nx,ny,nz,nzfull,
     &                            dxsqi,dysqi,dzsqi,
     &                            phi,rho,rstar,linbend,l2symtry,l4symtry,
     &                            scrtch,boundxy,globalb0,globalbnz,
     &                            my_index,nslaves,izfsslave,nzfsslave)
      use Multigrid3d
      integer(ISZ):: mglevel
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      real(kind=8):: rstar(-1:nz+1)
      real(kind=8):: scrtch(*)
      logical(ISZ):: linbend
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: boundxy,globalb0,globalbnz
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8),allocatable:: phi2(:,:,:),rho2(:,:,:)
      real(kind=8):: res(0:nx,0:ny,-1:nz+1)
      integer(ISZ):: i
      real(kind=8):: ff
      integer(ISZ):: nznew,lparity,rparity
      integer(ISZ):: localb0,localbnz
      integer(ISZ):: local2b0,local2bnz
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleft2(0:nslaves-1), izsendingleft2(0:nslaves-1)
      integer(ISZ):: whosendingright2(0:nslaves-1),izsendingright2(0:nslaves-1)
      integer(ISZ):: izfsslave2(0:nslaves-1),nzfsslave2(0:nslaves-1)

#ifdef MPIPARALLEL
      call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                      globalb0,globalbnz,nzfull,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
      localb0 = globalb0
      localbnz = globalbnz
      if (izfsslave(my_index) > 0) localb0 = -1
      if (izfsslave(my_index)+nz < nzfull) localbnz = -1
#else
      localb0 = globalb0
      localbnz = globalbnz
#endif

c     --- Do initial SOR passes.
c     print *, 'vcycle on dimensions ', nx, ny, nzfull,my_index
      do i=1,downpasses
        call sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,l2symtry,l4symtry,scrtch,
     &                 localb0,localbnz,my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

c     --- If grid is not at its coarsest level in any of the axis or and
c     --- all dimensions are even, continue the coarsening.
c     if (mod(nx,4) == 0 .and. mod(ny,4) == 0 .and. mod(nzfull,4) == 0) then
c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening.
      if (mglevel < mglevels-1) then

c       --- Get the residual on the current grid.
        call residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
     &                mglevel,localb0,localbnz,l2symtry,l4symtry,.false.,
     &                my_index,nslaves,izfsslave,nzfsslave,
     &                whosendingleft,izsendingleft,
     &                whosendingright,izsendingright)

c       --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
c       --- all axis.
c       --- dz > 4/3 dx <=> (9/16) / dx^2 < 1 / dz^2
        if (0.5625*dxsqi > dzsqi) then

c         --- Alloate new work space
          allocate(phi2(0:nx/2,0:ny/2,-1:nz+1),rho2(0:nx/2,0:ny/2,0:nz))
          phi2 = 0.
          rho2 = 0.

c         --- Ratio of old to new constant needed to scale the residual for
c         --- the restriction.
          ff = (dxsqi+dysqi+dzsqi)/(dxsqi*0.25 + dysqi*0.25 + dzsqi)
          call restrict2d(nx,ny,nz,nzfull,res,rho2,ff,boundxy,localb0,localbnz,
     &                    l2symtry,l4symtry)

c         --- Continue at the next coarsest level.
          call vcycle(mglevel+1,nx/2,ny/2,nz,nzfull,
     &                dxsqi*0.25,dysqi*0.25,dzsqi,phi2,rho2(0,0,0),
     &                rstar,linbend,l2symtry,l4symtry,scrtch,
     &                boundxy,globalb0,globalbnz,
     &                my_index,nslaves,izfsslave,nzfsslave)

c         --- Add in resulting error.
          call expand2d(nx/2,ny/2,nz,nzfull,phi2,phi,boundxy,localb0,localbnz)

          deallocate(phi2,rho2)
        else

#ifdef MPIPARALLEL
c         --- Find domains in coarser grid
          call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslave2,nzfsslave2,
     &                    nzfull)
c         --- Set new value of nz
          nznew = nzfsslave2(my_index)
c         --- Get processor with which to exchange data on coarse grid
          call mggetexchangepes(nslaves,izfsslave2,nzfsslave2,my_index,
     &                          globalb0,globalbnz,nzfull/2,
     &                          whosendingleft2,izsendingleft2,
     &                          whosendingright2,izsendingright2)
c         --- Difference between starts and ends of coarse and fine grids.
c         --- Should only be in the range 0-2.
          lparity = izfsslave(my_index) - 2*izfsslave2(my_index)
          rparity = 2*(izfsslave2(my_index)+nzfsslave2(my_index)) -
     &                (izfsslave(my_index)+nzfsslave(my_index))
          local2b0 = globalb0
          local2bnz = globalbnz
          if (izfsslave2(my_index) > 0) local2b0 = -1
          if (izfsslave2(my_index) + nznew < nzfull/2) local2bnz = -1
#else
          nznew = nz/2
          lparity = 0
          rparity = 0
          local2b0 = globalb0
          local2bnz = globalbnz
#endif

c         --- Alloate new work space
          allocate(phi2(0:nx/2,0:ny/2,-1:nznew+1),
     &             rho2(0:nx/2,0:ny/2,-1:nznew+1))
          phi2 = 0.
          rho2 = 0.

c         --- Restriction - note that scaling factor for residual is always
c         --- 4 for full-coarsening and is compiled into the restriction
c         --- routine.
          call restrict3d(nx,ny,nz,nznew,nzfull,res,rho2,boundxy,
     &                    local2b0,local2bnz,localb0,localbnz,
     &                    lparity,rparity,l2symtry,l4symtry,
     &                    my_index,nslaves,
     &                    izfsslave,nzfsslave,izfsslave2,nzfsslave2)

c         --- Continue at the next coarsest level.
          call vcycle(mglevel+1,nx/2,ny/2,nznew,nzfull/2,
     &                dxsqi*0.25,dysqi*0.25,dzsqi*0.25,phi2,rho2(0,0,0),
     &                rstar,linbend,l2symtry,l4symtry,scrtch,
     &                boundxy,globalb0,globalbnz,
     &                my_index,nslaves,izfsslave2,nzfsslave2)

c         --- Add in resulting error.
          call expand3d(nx/2,ny/2,nznew,nz,nzfull/2,phi2,phi,
     &                  boundxy,localb0,localbnz,lparity,rparity,
     &                  my_index,nslaves,izfsslave,
     &                  whosendingleft,whosendingright)

          deallocate(phi2,rho2)
        endif

      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbend,l2symtry,l4symtry,scrtch,
     &                 localb0,localbnz,my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

      return
      end
c=============================================================================
      subroutine restrict2d(nx,ny,nz,nzfull,res,rho2,ff,
     &                      boundxy,localb0,localbnz,l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: res(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho2(0:nx/2,0:ny/2,0:nz)
      real(kind=8):: ff
      integer(ISZ):: boundxy,localb0,localbnz
      logical(ISZ):: l2symtry,l4symtry

c Restrict transversely to a coarser grid.  The factor of ff is needed
c since the residual will be used as the source term in the next coarser
c grid and it implicitly includes the factor 0.5/(dxsqi+dysqi+dzsqi)
c which is ff times too small with the current values of dxsqi etc.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 2
      ixmax = nx-2
      iymin = 2
      iymax = ny-2
      izmin = 0
      izmax = nz
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny
      if (localb0  == 0) izmin = 1
      if (localbnz == 0) izmax = nz-1

c     --- Do the loops.
      do iz=izmin,izmax
        do iy=iymin,iymax,2
          iym1 = iy - 1
          if (iy==0 .and.  boundxy==2) iym1 = ny-1
          if (iy==0 .and. (boundxy==1 .or. l4symtry)) iym1 = 1
          iyp1 = iy + 1
          if (iy == ny .and. boundxy == 1)  iyp1 = ny-1
          if (iy == ny .and. boundxy == 2)  iyp1 = 1

          do ix=ixmin,ixmax,2
            ixm1 = ix - 1
            if (ix==0 .and.  boundxy==2) ixm1 = nx-1
            if (ix==0 .and. (boundxy==1 .or. l2symtry .or. l4symtry)) ixm1 = 1
            ixp1 = ix + 1
            if (ix == nx .and. boundxy == 1)  ixp1 = nx-1
            if (ix == nx .and. boundxy == 2)  ixp1 = 1

            rho2(ix/2,iy/2,iz) =
     &        ff*0.2500*res(ix  ,iy  ,iz) +
     &        ff*0.1250*res(ixm1,iy  ,iz) +
     &        ff*0.1250*res(ixp1,iy  ,iz) +
     &        ff*0.1250*res(ix  ,iym1,iz) +
     &        ff*0.1250*res(ix  ,iyp1,iz) +
     &        ff*0.0625*res(ixm1,iym1,iz) +
     &        ff*0.0625*res(ixp1,iym1,iz) +
     &        ff*0.0625*res(ixm1,iyp1,iz) +
     &        ff*0.0625*res(ixp1,iyp1,iz)

          enddo
        enddo
      enddo

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        rho2(:,ny/2,:) = rho2(:,0,:)
        rho2(nx/2,:,:) = rho2(0,:,:)
      endif

      if (localb0  == 0) rho2(:,:,0) = 0.
      if (localbnz == 0) rho2(:,:,nz) = 0.

      return
      end
c=============================================================================
      subroutine restrict3d(nx,ny,nz,nznew,nzfull,res,rho2,boundxy,
     &                      local2b0,local2bnz,localb0,localbnz,
     &                      lparity,rparity,l2symtry,l4symtry,
     &                      my_index,nslaves,
     &                      izfsslave,nzfsslave,izfsslave2,nzfsslave2)
      integer(ISZ):: nx,ny,nz,nznew,nzfull
      real(kind=8):: res(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho2(0:nx/2,0:ny/2,-1:nznew+1)
      integer(ISZ):: boundxy,local2b0,local2bnz,localb0,localbnz
      integer(ISZ):: lparity,rparity
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: izfsslave2(0:nslaves-1),nzfsslave2(0:nslaves-1)

c Restrict to a coarser grid.  The factor of 4 is needed since the
c residual will be used as the source term in the next coarser grid and
c it implicitly includes the factor 0.5/(dxsqi+dysqi+dzsqi) which is 4
c times too small with the current values of dxsqi etc.

      integer(ISZ):: ix,iy,iz,ic
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1,iz00

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 2
      ixmax = nx-2
      iymin = 2
      iymax = ny-2
      izmin = 0
      izmax = nznew
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny
      if (localb0 == 0 .or. lparity >= 1) izmin = 1
      if (localbnz == 0 .or. rparity >= 1) izmax = nznew - 1

c     --- Do the loops.
      do iz=izmin,izmax
        izm1 = 2*iz - 1 - lparity
        izp1 = 2*iz + 1 - lparity
        iz00 = 2*iz     - lparity

        do iy=iymin,iymax,2
          iym1 = iy - 1
          if (iy==0 .and.  boundxy==2) iym1 = ny-1
          if (iy==0 .and. (boundxy==1 .or. l4symtry)) iym1 = 1
          iyp1 = iy + 1
          if (iy == ny .and. boundxy == 1)  iyp1 = ny-1
          if (iy == ny .and. boundxy == 2)  iyp1 = 1

          do ix=ixmin,ixmax,2
            ixm1 = ix - 1
            if (ix==0 .and.  boundxy==2) ixm1 = nx-1
            if (ix==0 .and. (boundxy==1 .or. l2symtry .or. l4symtry)) ixm1 = 1
            ixp1 = ix + 1
            if (ix == nx .and. boundxy == 1)  ixp1 = nx-1
            if (ix == nx .and. boundxy == 2)  ixp1 = 1

            rho2(ix/2,iy/2,iz) =
     &       4.*0.125000*res(ix  ,iy  ,iz00) +
     &       4.*0.062500*res(ixm1,iy  ,iz00) + 4.*0.062500*res(ixp1,iy  ,iz00) +
     &       4.*0.062500*res(ix  ,iym1,iz00) + 4.*0.062500*res(ix  ,iyp1,iz00) +
     &       4.*0.062500*res(ix  ,iy  ,izm1) + 4.*0.062500*res(ix  ,iy  ,izp1) +
     &       4.*0.031250*res(ixm1,iym1,iz00) + 4.*0.031250*res(ixp1,iym1,iz00) +
     &       4.*0.031250*res(ixm1,iyp1,iz00) + 4.*0.031250*res(ixp1,iyp1,iz00) +
     &       4.*0.031250*res(ixm1,iy  ,izm1) + 4.*0.031250*res(ixp1,iy  ,izm1) +
     &       4.*0.031250*res(ixm1,iy  ,izp1) + 4.*0.031250*res(ixp1,iy  ,izp1) +
     &       4.*0.031250*res(ix  ,iym1,izm1) + 4.*0.031250*res(ix  ,iyp1,izm1) +
     &       4.*0.031250*res(ix  ,iym1,izp1) + 4.*0.031250*res(ix  ,iyp1,izp1) +
     &       4.*0.015625*res(ixm1,iym1,izm1) + 4.*0.015625*res(ixp1,iym1,izm1) +
     &       4.*0.015625*res(ixm1,iyp1,izm1) + 4.*0.015625*res(ixp1,iyp1,izm1) +
     &       4.*0.015625*res(ixm1,iym1,izp1) + 4.*0.015625*res(ixp1,iym1,izp1) +
     &       4.*0.015625*res(ixm1,iyp1,izp1) + 4.*0.015625*res(ixp1,iyp1,izp1)

          enddo
        enddo
      enddo

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        rho2(:,ny/2,:) = rho2(:,0,:)
        rho2(nx/2,:,:) = rho2(0,:,:)
      endif
      if (local2b0  == 0) rho2(:,:,0) = 0.
      if (local2bnz == 0) rho2(:,:,nznew) = 0.
#ifdef MPIPARALLEL
      call mgexchange_res(nx/2,ny/2,nznew,nzfull/2,rho2,my_index,nslaves,
     &                    izfsslave,nzfsslave,izfsslave2,nzfsslave2)
#endif

      return
      end
c=============================================================================
      subroutine expand2d(nx,ny,nz,nzfull,p2h,ph,boundxy,localb0,localbnz)
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: p2h(0:nx,0:ny,-1:nz+1)
      real(kind=8):: ph(0:nx*2,0:ny*2,-1:nz+1)
      integer(ISZ):: boundxy,localb0,localbnz

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: izmin,izmax

c     --- Set the loop limits, including edges when appropriate.
      izmin = -1
      izmax = nz+1
      if (localb0  == 0) izmin = 1
      if (localbnz == 0) izmax = nz-1

c     --- Do the loops.
      do iz=izmin,izmax
        do iy=0,ny-1
          do ix=0,nx-1

            ph(2*ix  ,2*iy  ,iz)=ph(2*ix  ,2*iy  ,iz)+      p2h(ix  ,iy  ,iz)
            ph(2*ix+1,2*iy  ,iz)=ph(2*ix+1,2*iy  ,iz)+0.50*(p2h(ix  ,iy  ,iz)+
     &                                                      p2h(ix+1,iy  ,iz))
            ph(2*ix  ,2*iy+1,iz)=ph(2*ix  ,2*iy+1,iz)+0.50*(p2h(ix  ,iy  ,iz)+
     &                                                      p2h(ix  ,iy+1,iz))
            ph(2*ix+1,2*iy+1,iz)=ph(2*ix+1,2*iy+1,iz)+0.25*(p2h(ix  ,iy  ,iz)+
     &                                                      p2h(ix+1,iy  ,iz)+
     &                                                      p2h(ix  ,iy+1,iz)+
     &                                                      p2h(ix+1,iy+1,iz))
          enddo
        enddo
      enddo

      if (boundxy > 0) then
        do iz=izmin,izmax

c         --- Expand ix=nx, iy=ny line.
          ph(2*nx,2*ny,iz) = ph(2*nx,2*ny,iz)+p2h(nx,ny,iz)

c         --- Expand ix=nx plane.
          do iy=0,ny-1
            ph(2*nx,2*iy  ,iz) = ph(2*nx,2*iy  ,iz) +      p2h(nx,iy  ,iz)
            ph(2*nx,2*iy+1,iz) = ph(2*nx,2*iy+1,iz) + 0.5*(p2h(nx,iy  ,iz) +
     &                                                     p2h(nx,iy+1,iz))
          enddo

c         --- Expand iy=ny plane.
          do ix=0,nx-1
            ph(2*ix  ,2*ny,iz) = ph(2*ix  ,2*ny,iz) +      p2h(ix  ,ny,iz)
            ph(2*ix+1,2*ny,iz) = ph(2*ix+1,2*ny,iz) + 0.5*(p2h(ix  ,ny,iz) +
     &                                                     p2h(ix+1,ny,iz))
          enddo

        enddo
      endif

      return
      end
c=============================================================================
      subroutine expand3d(nx,ny,nz,nz2,nzfull,p2h,ph,boundxy,localb0,localbnz,
     &                    lparity,rparity,my_index,nslaves,izfsslave,
     &                    whosendingleft,whosendingright)
      integer(ISZ):: nx,ny,nz,nz2,nzfull
      real(kind=8):: p2h(0:nx,0:ny,-1:nz+1)
      real(kind=8):: ph(0:nx*2,0:ny*2,-1:nz2+1)
      integer(ISZ):: boundxy,localb0,localbnz,lparity,rparity
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1),whosendingright(0:nslaves-1)

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is done along all axis.

      integer(ISZ):: ix,iy,iz,izmin,izmax
      integer(ISZ):: iz2

c     --- If lparity of even, start iz loop at 0 (since 0 lines up with
c     --- cell in coarse grid). If lparity is odd, start iz loop at -1 (since
c     --- -1 lines up with cell in coarse grid). It is almost certainly true
c     --- that lparity is only 0 or 1, and lparity is never negative.
c     --- Same for rparity.
      izmin = -mod(lparity,2)
      izmax = nz2 - mod(rparity,2)

c     --- Expand the bulk of the grid.
      do iz=izmin,izmax,2
        iz2 = (iz + lparity)/2
        do iy=0,ny-1
          do ix=0,nx-1

      ph(2*ix  ,2*iy  ,iz  )=ph(2*ix  ,2*iy  ,iz  ) +
     &                                                   p2h(ix  ,iy  ,iz2  )
      ph(2*ix+1,2*iy  ,iz  )=ph(2*ix+1,2*iy  ,iz  ) +
     &                                            0.500*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix+1,iy  ,iz2  ))
      ph(2*ix  ,2*iy+1,iz  )=ph(2*ix  ,2*iy+1,iz  ) +
     &                                            0.500*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix  ,iy+1,iz2  ))
      ph(2*ix  ,2*iy  ,iz+1)=ph(2*ix  ,2*iy  ,iz+1) +
     &                                            0.500*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix  ,iy  ,iz2+1))
      ph(2*ix+1,2*iy+1,iz  )=ph(2*ix+1,2*iy+1,iz  ) +
     &                                            0.250*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix+1,iy  ,iz2  ) +
     &                                                   p2h(ix  ,iy+1,iz2  ) +
     &                                                   p2h(ix+1,iy+1,iz2  ))
      ph(2*ix+1,2*iy  ,iz+1)=ph(2*ix+1,2*iy  ,iz+1) +
     &                                            0.250*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix+1,iy  ,iz2  ) +
     &                                                   p2h(ix  ,iy  ,iz2+1) +
     &                                                   p2h(ix+1,iy  ,iz2+1))
      ph(2*ix  ,2*iy+1,iz+1)=ph(2*ix  ,2*iy+1,iz+1) +
     &                                            0.250*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix  ,iy+1,iz2  ) +
     &                                                   p2h(ix  ,iy  ,iz2+1) +
     &                                                   p2h(ix  ,iy+1,iz2+1))
      ph(2*ix+1,2*iy+1,iz+1)=ph(2*ix+1,2*iy+1,iz+1) +
     &                                            0.125*(p2h(ix  ,iy  ,iz2  ) +
     &                                                   p2h(ix+1,iy  ,iz2  ) +
     &                                                   p2h(ix  ,iy+1,iz2  ) +
     &                                                   p2h(ix+1,iy+1,iz2  ) +
     &                                                   p2h(ix  ,iy  ,iz2+1) +
     &                                                   p2h(ix+1,iy  ,iz2+1) +
     &                                                   p2h(ix  ,iy+1,iz2+1) +
     &                                                   p2h(ix+1,iy+1,iz2+1))

          enddo
        enddo
      enddo

      if (boundxy > 0) then
        do iz=izmin,izmax,2
          iz2 = (iz + lparity)/2

c         --- Expand ix=nx, iy=ny line.
          ph(2*nx,2*ny,iz  )=ph(2*nx,2*ny,iz  ) +        p2h(nx,ny,iz2  )
          ph(2*nx,2*ny,iz+1)=ph(2*nx,2*ny,iz+1) + 0.500*(p2h(nx,ny,iz2  ) +
     &                                                   p2h(nx,ny,iz2+1))

c         --- Expand ix=nx plane.
          do iy=0,ny-1
            ph(2*nx,2*iy  ,iz  )=ph(2*nx,2*iy  ,iz  ) +
     &                                                      p2h(nx,iy  ,iz2  )
            ph(2*nx,2*iy+1,iz  )=ph(2*nx,2*iy+1,iz  )+
     &                                               0.500*(p2h(nx,iy  ,iz2  ) +
     &                                                      p2h(nx,iy+1,iz2  ))
            ph(2*nx,2*iy  ,iz+1)=ph(2*nx,2*iy  ,iz+1)+
     &                                               0.500*(p2h(nx,iy  ,iz2  ) +
     &                                                      p2h(nx,iy  ,iz2+1))
            ph(2*nx,2*iy+1,iz+1)=ph(2*nx,2*iy+1,iz+1)+
     &                                               0.250*(p2h(nx,iy  ,iz2  ) +
     &                                                      p2h(nx,iy+1,iz2  ) +
     &                                                      p2h(nx,iy  ,iz2+1) +
     &                                                      p2h(nx,iy+1,iz2+1))
          enddo

c         --- Expand iy=ny plane.
          do ix=0,nx-1
            ph(2*ix  ,2*ny,iz  ) = ph(2*ix  ,2*ny,iz  ) +
     &                                                      p2h(ix  ,ny,iz2  )
            ph(2*ix+1,2*ny,iz  ) = ph(2*ix+1,2*ny,iz  )+
     &                                               0.500*(p2h(ix  ,ny,iz2  ) +
     &                                                      p2h(ix+1,ny,iz2  ))
            ph(2*ix  ,2*ny,iz+1) = ph(2*ix  ,2*ny,iz+1)+
     &                                               0.500*(p2h(ix  ,ny,iz2  ) +
     &                                                      p2h(ix  ,ny,iz2+1))
            ph(2*ix+1,2*ny,iz+1) = ph(2*ix+1,2*ny,iz+1)+
     &                                               0.250*(p2h(ix  ,ny,iz2  ) +
     &                                                      p2h(ix+1,ny,iz2  ) +
     &                                                      p2h(ix  ,ny,iz2+1) +
     &                                                      p2h(ix+1,ny,iz2+1))
          enddo

        enddo
      endif

      if (localb0  == 1) ph(:,:,-1) = ph(:,:,1)
      if (localbnz == 1) ph(:,:,nz2+1) = ph(:,:,nz2-1)
#ifndef MPIPARALLEL
      if (localb0  == 2) ph(:,:,-1) = ph(:,:,nz2-1)
      if (localbnz == 2) ph(:,:,nz2+1) = ph(:,:,1)
#else
      call mgexchange_phiperiodic(2*nx,2*ny,nz2,2*nzfull,ph,localb0,localbnz,1,
     &                            my_index,nslaves,izfsslave,
     &                            whosendingleft,whosendingright)
#endif

      return
      end
c=============================================================================
c=============================================================================
      subroutine sorpass3d(mglevel,nx,ny,nz,nzfull,phi,rho,rstar,
     &                     rdx2,rdy2,rdz2,linbend,l2symtry,l4symtry,scrtch,
     &                     localb0,localbnz,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use Constant
      integer(ISZ):: mglevel,nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz),rstar(-1:nz+1)
      real(kind=8):: scrtch(*)
      real(kind=8):: rdx2,rdy2,rdz2
      logical(ISZ):: linbend
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: localb0,localbnz
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

      real(kind=8):: rdel,const,rdx2c,rdy2c,rdz2c,spm1,dx
      real(kind=8):: boundsavex(0:nx,2,0:nz)
      real(kind=8):: boundsavey(0:ny,2,0:nz)
      integer(ISZ):: nxy,iimx,iipx,iimy,iipy,iimz,iipz,izl,izr,parity
      integer(ISZ):: ii,ix,iy,iz,ic
      integer(ISZ):: s_parity,e_parity

c     --- Macro for common arguments to the psor_loop calls
#define LOOP_ARGS iimz,iipz,phi,rho,nxy,rdx2c,rdy2c,rdz2c,spm1,linbend, \
      rstar(-1),dx,rdx2,rdy2,rdz2,rdel,const,scrtch(1)

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
      if (localb0  < 1) izl = izl + 1
      if (localbnz < 1) izr = izr - 1

c     --- Put desired potential onto conductors in phi array.
      call cond_potmg(nx,ny,nz,phi,mglevel,.false.,mgform,.false.)

c     --- Save values on the transverse boundaries.
c     --- Both even and odd points are saved for all transverse boundaries.
      do iz=izl-1,izr-1
        do ix=0,nx
          boundsavex(ix,1,iz) = phi(ix,0,iz)
          boundsavex(ix,2,iz) = phi(ix,ny,iz)
        enddo
        do iy=0,ny
          boundsavey(iy,1,iz) = phi(0,iy,iz)
          boundsavey(iy,2,iz) = phi(nx,iy,iz)
        enddo
      enddo

c     --- Save values just outside conductor surfaces. Only save phi at the
c     --- subgrid points which are to be used at the current level of
c     --- grid refinement.
      if (lcndbndy) then
        do ic = ecstart(mglevel),ecstart(mglevel+1)-1
          ecndpvph(ic) = phi(iecndx(ic),iecndy(ic),iecndz(ic))
        enddo
        do ic = ocstart(mglevel),ocstart(mglevel+1)-1
          ocndpvph(ic) = phi(iocndx(ic),iocndy(ic),iocndz(ic))
        enddo
      endif

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

c       --- guard planes in z are already set

c       --- Loop over the rest of the array. Boundary points are calculated
c       --- too, even though the equation is wrong.  They are recalculated
c       --- later.  The case with and without bends are seperated since
c       --- the loops are done differently.  Without bends, the whole array
c       --- can be processed in one loop as a 1d array.  With bends, an
c       --- outer loop over z is needed.
        if (.not. linbend) then 
c         --- no bends 
          call psor_loop(parity,1,0,0,1,izl*nxy+2,(izr+1)*nxy,2,
     &                   iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
        else 
c         --- bends
          call psor_loop(parity,1,izl,izr,1,2,nxy,2,
     &                   iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
        endif 

c       --- Transverse boundaries     
c       --- Restore only even or odd boundary points to previous value
c       --- since only want to restore the values changed from the 1d
c       --- loop.  This automatically takes care of Dirichlet boundaries.
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

        if (boundxy == 1 .or. l2symtry .or. l4symtry) then
c         --- surfaces at ix=0, ix=nx, iy=0, and iy=ny
c         --- if only 2-fold, apply to surfaces at iy=0
c         --- if only 4-fold, apply to surfaces at ix=0 and iy=0
c         --- if also boundxy=1, then apply to all transverse surfaces
          call psor_loop(parity,1,izl,izr,1,
     &                   2,nx,2,iimx,iipx,-iimy,iipy,
     &LOOP_ARGS)
          if (boundxy == 1 .or. l4symtry) then
            call psor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,
     &                     2*(nx+1),-iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
c           --- line at transverse edge (ix=0, iy=0)
            call psor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,
     &                     2*nxy,-iimx,iipx,-iimy,iipy,
     &LOOP_ARGS)
          endif
c         --- now do planes at ix=nx and iy=ny
          if (boundxy == 1) then
            call psor_loop(parity,1,izl,izr,1,(nx+1)*ny+2,(nx+1)*(ny+1)-1,2,
     &                     iimx,iipx,iimy,-iipy,
     &LOOP_ARGS)
            call psor_loop(parity,nx+1,izl,izr,1,2*(nx+1),(nx+1)*ny,
     &                     2*(nx+1),iimx,-iipx,iimy,iipy,
     &LOOP_ARGS)
c           --- lines at other transverse edges
            call psor_loop(parity,nxy,izl,izl,1,nx+1,nxy*(izr-izl)+nx+1,
     &                     2*nxy,iimx,-iipx,-iimy,iipy,
     &LOOP_ARGS)
            call psor_loop(parity,nxy,izl,izl,1,nxy-nx,nxy*(izr-izl)+nxy-nx,
     &                     2*nxy,-iimx,iipx,iimy,-iipy,
     &LOOP_ARGS)
            call psor_loop(parity,nxy,izl,izl,1,nxy,nxy*(izr-izl)+nxy,
     &                     2*nxy,iimx,-iipx,iimy,-iipy,
     &LOOP_ARGS)
          endif
        else if (boundxy == 2) then
c         --- surfaces at ix=0, ix=nx, iy=0, and iy=ny
          call psor_loop(parity,1,izl,izr,1,
     &                   2,nx,2,iimx,iipx,(nx+1)*(ny-1),iipy,
     &LOOP_ARGS)
          call psor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,
     &                   2*(nx+1),nx-1,iipx,iimy,iipy,
     &LOOP_ARGS)
          do iz=izl-1,izr-1
            do ix=1+mod(iz+parity+1,2),nx-1,2
              phi(ix,ny,iz) = phi(ix,0,iz)
            enddo
            do iy=1+mod(iz+parity+1,2),ny-1,2
              phi(nx,iy,iz) = phi(0,iy,iz)
            enddo
          enddo
c         --- lines at transverse edges
          call psor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,
     &                   2*nxy,nx-1,iipx,(nx+1)*(ny-1),iipy,
     &LOOP_ARGS)
          do iz=izl-1+mod(izl-1+parity,2),izr-1,2
            phi(nx,0,iz) = phi(0,0,iz)
            phi(0,ny,iz) = phi(0,0,iz)
            phi(nx,ny,iz) = phi(0,0,iz)
          enddo
        endif
c       --- end of transverse boundaries

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
          if (parity == 0) then
           call condbndymg(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,mgparam,
     &                    ncndmax,necndbdy,ecdelmx,ecdelmy,ecdelmz,
     &                    ecdelpx,ecdelpy,ecdelpz,iecndx,iecndy,iecndz,
     &                    iecndlevel,
     &                    ecvoltmx,ecvoltpx,ecvoltmy,ecvoltpy,ecvoltmz,ecvoltpz,
     &                    ecndpvph,localb0,localbnz,boundxy,l2symtry,l4symtry,
     &                    mglevel,ecstart,mgform)
          endif
          if (parity == 1) then
           call condbndymg(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,mgparam,
     &                    ncndmax,nocndbdy,ocdelmx,ocdelmy,ocdelmz,
     &                    ocdelpx,ocdelpy,ocdelpz,iocndx,iocndy,iocndz,
     &                    iocndlevel,
     &                    ocvoltmx,ocvoltpx,ocvoltmy,ocvoltpy,ocvoltmz,ocvoltpz,
     &                    ocndpvph,localb0,localbnz,boundxy,l2symtry,l4symtry,
     &                    mglevel,ocstart,mgform)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmg(nx,ny,nz,phi,mglevel,.false.,mgform,.false.)

c       --- set phi in the z guard planes
        if (localb0  == 1) phi(:,:,izl-2) = phi(:,:,izl)
        if (localbnz == 1) phi(:,:,izr) = phi(:,:,izr-2)
#ifndef MPIPARALLEL
        if (localb0  == 2) phi(:,:,izl-2) = phi(:,:,izr-2)
        if (localbnz == 2) phi(:,:,nz:nz+1) = phi(:,:,0:1)
#else
        call mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,0,
     &                      my_index,nslaves,izfsslave,nzfsslave,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
        call mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,localb0,localbnz,1,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- Exchange phi in the z guard planes
      call mgexchange_phi(nx,ny,nz,nzfull,phi,localb0,localbnz,-1,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
#endif

      return
      end
c=============================================================================
      subroutine cond_potmg(nx,ny,nz,phi,mglevel,lresidual,
     &                      mgform,mgform2init)
      use PSOR3d
      use Conductor3d
      integer(ISZ):: nx,ny,nz,mglevel,mgform
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      logical(ISZ):: lresidual,mgform2init

c Set conductor points to the desired potential.

      integer(ISZ):: ic

c     --- When at the finest level and not calculating the residual, set
c     --- phi to the voltage of the conductor, otherwise, set it to zero.
      if (mglevel == 0 .and. .not. lresidual .and.
     &    (mgform == 1 .or. mgform2init)) then
        do ic = icstart(mglevel),icstart(mglevel+1)-1
          phi(ixcond(ic),iycond(ic),izcond(ic)) = condvolt(ic)
        enddo
      else
        do ic=icstart(mglevel),icstart(mglevel+1)-1
          phi(ixcond(ic),iycond(ic),izcond(ic)) = 0.
        enddo
      endif

      return
      end
c=============================================================================
      subroutine condbndymg(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,srp,
     &                      ncndmax,ncndbdy,cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,
     &                      cdelpz,icndx,icndy,icndz,icndlevel,
     &                      cvoltmx,cvoltpx,cvoltmy,cvoltpy,cvoltmz,cvoltpz,
     &                      cndpvph,localb0,localbnz,boundxy,
     &                      l2symtry,l4symtry,mglevel,cstart,mgform)
      integer(ISZ):: nx,ny,nz,ncndmax,ncndbdy,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), rho(0:nx,0:ny,0:nz)
      real(kind=8):: rdx2c,rdy2c,rdz2c,spm1,srp
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax), cdelmz(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax), cdelpz(ncndmax)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax), icndz(ncndmax)
      integer(ISZ):: icndlevel(ncndmax)
      real(kind=8):: cvoltmx(ncndmax),cvoltmy(ncndmax),cvoltmz(ncndmax)
      real(kind=8):: cvoltpx(ncndmax),cvoltpy(ncndmax),cvoltpz(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
      integer(ISZ):: localb0,localbnz,boundxy,cstart(0:100),mgform
      logical(ISZ):: l2symtry,l4symtry

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
c to hold icndx-1, icndx+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: rdx2cos,rdy2cos,rdz2cos,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: voltfac
      integer(ISZ):: ic,icndxp1,icndxm1,icndyp1,icndym1,icndzp1,icndzm1
      integer(ISZ):: icdx,icdy,icdz,il
      logical(ISZ):: dosubgrid

      rdx2cos = rdx2c/srp
      rdy2cos = rdy2c/srp
      rdz2cos = rdz2c/srp

c     --- loop over points near surface of conductors
      do ic = cstart(mglevel),cstart(mglevel+1)-1

        icdx = icndx(ic)
        icdy = icndy(ic)
        icdz = icndz(ic)

c       if (mod(icdx+icdy+icdz,2) /= parity) cycle

c       --- set temporaries for boundaries
        icndxp1 = icdx + 1
        icndxm1 = icdx - 1
        icndyp1 = icdy + 1
        icndym1 = icdy - 1
        icndzp1 = icdz + 1
        icndzm1 = icdz - 1

        if (icndxm1 == -1 .and. boundxy == 2) icndxm1 = nx-1
        if (icndxm1 == -1 .and.
     &     (boundxy == 1 .or. l2symtry .or. l4symtry)) icndxm1 = 1
        if (icndxp1 == nx+1 .and. boundxy == 1) icndxp1 = nx-1
        if (icndxp1 == nx+1 .and. boundxy == 2) icndxp1 = 1
        if (icndxm1 == -1 .or. icndxp1 == nx+1) cycle

        if (icndym1 == -1 .and. boundxy == 2) icndym1 = ny-1
        if (icndym1 == -1 .and.
     &     (boundxy == 1 .or. l4symtry)) icndym1 = 1
        if (icndyp1 == ny+1 .and. boundxy == 1) icndyp1 = ny-1
        if (icndyp1 == ny+1 .and. boundxy == 2) icndyp1 = 1
        if (icndym1 == -1 .or. icndyp1 == ny+1) cycle

        if (icdz == 0  .and. localb0  < 1) cycle
        if (icdz == nz .and. localbnz < 1) cycle

c       --- Set temporaries with initial values.
        pxm = phi(icndxm1,icdy   ,icdz   )
        pxp = phi(icndxp1,icdy   ,icdz   )
        pym = phi(icdx   ,icndym1,icdz   )
        pyp = phi(icdx   ,icndyp1,icdz   )
        pzm = phi(icdx   ,icdy   ,icndzm1)
        pzp = phi(icdx   ,icdy   ,icndzp1)
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
        if (cdelmx(ic) < 1.) then
          pxm = voltfac*cvoltmx(ic)/cdelmx(ic)
          denom = denom + (1.-cdelmx(ic))/cdelmx(ic)*rdx2cos
          dosubgrid = .true.
        endif
c       --- the point higher in x is inside the conductor
        if (cdelpx(ic) < 1.) then
          pxp = voltfac*cvoltpx(ic)/cdelpx(ic)
          denom = denom + (1.-cdelpx(ic))/cdelpx(ic)*rdx2cos
          dosubgrid = .true.
        endif
c       --- the point lower in y is inside the conductor
        if (cdelmy(ic) < 1.) then
          pym = voltfac*cvoltmy(ic)/cdelmy(ic)
          denom = denom + (1.-cdelmy(ic))/cdelmy(ic)*rdy2cos
          dosubgrid = .true.
        endif
c       --- the point higher in y is inside the conductor
        if (cdelpy(ic) < 1.) then
          pyp = voltfac*cvoltpy(ic)/cdelpy(ic)
          denom = denom + (1.-cdelpy(ic))/cdelpy(ic)*rdy2cos
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (cdelmz(ic) < 1.) then
          pzm = voltfac*cvoltmz(ic)/cdelmz(ic)
          denom = denom + (1.-cdelmz(ic))/cdelmz(ic)*rdz2cos
          dosubgrid = .true.
        endif
c       --- the point lower in z is inside the conductor
        if (cdelpz(ic) < 1.) then
          pzp = voltfac*cvoltpz(ic)/cdelpz(ic)
          denom = denom + (1.-cdelpz(ic))/cdelpz(ic)*rdz2cos
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
          phi(icdx,icdy,icdz) = (rho(icdx,icdy,icdz) +
     &      (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c + (pzm+pzp)*rdz2c)/denom +
     &      spm1*cndpvph(ic)
        endif
      enddo

      return
      end
c=============================================================================
      subroutine condbndyres(nx,ny,nz,phi,rho,res,rdx2,rdy2,rdz2,mgparam,
     &                       ncndmax,ncndbdy,cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,
     &                       cdelpz,icndx,icndy,icndz,icndlevel,
     &                       cvoltmx,cvoltpx,cvoltmy,cvoltpy,cvoltmz,cvoltpz,
     &                       cndpvph,localb0,localbnz,boundxy,
     &                       l2symtry,l4symtry,mglevel,cstart,
     &                       mgform,mgform2init)
      integer(ISZ):: nx,ny,nz,ncndmax,ncndbdy,mglevel
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rho(0:nx,0:ny,0:nz), res(0:nx,0:ny,-1:nz+1)
      real(kind=8):: rdx2,rdy2,rdz2,mgparam
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax), cdelmz(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax), cdelpz(ncndmax)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax), icndz(ncndmax)
      integer(ISZ):: icndlevel(ncndmax)
      real(kind=8):: cvoltmx(ncndmax),cvoltmy(ncndmax),cvoltmz(ncndmax)
      real(kind=8):: cvoltpx(ncndmax),cvoltpy(ncndmax),cvoltpz(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
      integer(ISZ):: localb0,localbnz,boundxy,cstart(0:100),mgform
      logical(ISZ):: l2symtry,l4symtry,mgform2init

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
c to hold icndx-1, icndx+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: const,rdx2c,rdy2c,rdz2c,pxm,pym,pzm,pxp,pyp,pzp,denom
      real(kind=8):: voltfac
      real(kind=8):: rdx2cs,rdy2cs,rdz2cs,ppp
      integer(ISZ):: ic,icndxp1,icndxm1,icndyp1,icndym1,icndzp1,icndzm1
      integer(ISZ):: icdx,icdy,icdz,il

      const = 0.5/(rdx2+rdy2+rdz2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdz2c = rdz2*const
      rdx2cs = mgparam*rdx2*const
      rdy2cs = mgparam*rdy2*const
      rdz2cs = mgparam*rdz2*const

c     --- loop over points near surface of conductors
      do ic = cstart(mglevel),cstart(mglevel+1)-1

        icdx = icndx(ic)
        icdy = icndy(ic)
        icdz = icndz(ic)

c       --- set temporaries for boundaries
        icndxp1 = icdx + 1
        icndxm1 = icdx - 1
        icndyp1 = icdy + 1
        icndym1 = icdy - 1
        icndzp1 = icdz + 1
        icndzm1 = icdz - 1

        if (icndxm1 == -1 .and. boundxy == 2) icndxm1 = nx-1
        if (icndxm1 == -1 .and.
     &     (boundxy == 1 .or. l2symtry .or. l4symtry)) icndxm1 = 1
        if (icndxp1 == nx+1 .and. boundxy == 1) icndxp1 = nx-1
        if (icndxp1 == nx+1 .and. boundxy == 2) icndxp1 = 1
        if (icndxm1 == -1 .or. icndxp1 == nx+1) cycle

        if (icndym1 == -1 .and. boundxy == 2) icndym1 = ny-1
        if (icndym1 == -1 .and.
     &     (boundxy == 1 .or. l4symtry)) icndym1 = 1
        if (icndyp1 == ny+1 .and. boundxy == 1) icndyp1 = ny-1
        if (icndyp1 == ny+1 .and. boundxy == 2) icndyp1 = 1
        if (icndym1 == -1 .or. icndyp1 == ny+1) cycle

        if (icdz == 0  .and. localb0  == 0) cycle
        if (icdz == nz .and. localbnz == 0) cycle

c       --- set temporaries with initial values
        pxm = phi(icndxm1,icdy   ,icdz   )
        pxp = phi(icndxp1,icdy   ,icdz   )
        pym = phi(icdx   ,icndym1,icdz   )
        pyp = phi(icdx   ,icndyp1,icdz   )
        pzm = phi(icdx   ,icdy   ,icndzm1)
        pzp = phi(icdx   ,icdy   ,icndzp1)
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
        if (cdelmx(ic) < 1.) then
          pxm = voltfac*cvoltmx(ic)/cdelmx(ic)
          denom = denom + (1.-cdelmx(ic))/cdelmx(ic)*rdx2c
          ppp = min(ppp,cdelmx(ic))
        endif
c       --- the point higher in x is inside the conductor
        if (cdelpx(ic) < 1.) then
          pxp = voltfac*cvoltpx(ic)/cdelpx(ic)
          denom = denom + (1.-cdelpx(ic))/cdelpx(ic)*rdx2c
          ppp = min(ppp,cdelpx(ic))
        endif
c       --- the point lower in y is inside the conductor
        if (cdelmy(ic) < 1.) then
          pym = voltfac*cvoltmy(ic)/cdelmy(ic)
          denom = denom + (1.-cdelmy(ic))/cdelmy(ic)*rdy2c
          ppp = min(ppp,cdelmy(ic))
        endif
c       --- the point higher in y is inside the conductor
        if (cdelpy(ic) < 1.) then
          pyp = voltfac*cvoltpy(ic)/cdelpy(ic)
          denom = denom + (1.-cdelpy(ic))/cdelpy(ic)*rdy2c
          ppp = min(ppp,cdelpy(ic))
        endif
c       --- the point lower in z is inside the conductor
        if (cdelmz(ic) < 1.) then
          pzm = voltfac*cvoltmz(ic)/cdelmz(ic)
          denom = denom + (1.-cdelmz(ic))/cdelmz(ic)*rdz2c
          ppp = min(ppp,cdelmz(ic))
        endif
c       --- the point higher in z is inside the conductor
        if (cdelpz(ic) < 1.) then
          pzp = voltfac*cvoltpz(ic)/cdelpz(ic)
          denom = denom + (1.-cdelpz(ic))/cdelpz(ic)*rdz2c
          ppp = min(ppp,cdelpz(ic))
        endif
c       --- calculate the residual based on the boundary conditions
        if (ppp < 1.) then
          if (mgform2init) ppp = 1.
          res(icdx,icdy,icdz) = ppp*(rho(icdx,icdy,icdz)
     &           + (pxm+pxp)*rdx2cs + (pym+pyp)*rdy2cs + (pzm+pzp)*rdz2cs
     &           - phi(icdx,icdy,icdz)*mgparam*denom)
        endif
      enddo

      return
      end
c=============================================================================
      subroutine residual(nx,ny,nz,nzfull,dxsqi,dysqi,dzsqi,phi,rho,res,
     &                    mglevel,localb0,localbnz,l2symtry,l4symtry,
     &                    mgform2init,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      real(kind=8):: res(0:nx,0:ny,-1:nz+1)
      integer(ISZ):: mglevel,localb0,localbnz
      logical(ISZ):: l2symtry,l4symtry,mgform2init
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

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
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny
      if (localb0 == 0) izmin = 1
      if (localbnz == 0) izmax = nz-1

c     --- Calculate the residual.
      do iz=izmin,izmax
        izm1 = iz - 1
        izp1 = iz + 1

        do iy=iymin,iymax
          iym1 = iy - 1
          if (iy==0 .and.  boundxy==2) iym1 = ny-1
          if (iy==0 .and. (boundxy==1 .or. l4symtry)) iym1 = 1
          iyp1 = iy + 1
          if (iy == ny .and. boundxy == 1)  iyp1 = ny-1
          if (iy == ny .and. boundxy == 2)  iyp1 = 1

          do ix=ixmin,ixmax
            ixm1 = ix - 1
            if (ix==0 .and.  boundxy==2) ixm1 = nx-1
            if (ix==0 .and. (boundxy==1 .or. l2symtry .or. l4symtry)) ixm1 = 1
            ixp1 = ix + 1
            if (ix == nx .and. boundxy == 1)  ixp1 = nx-1
            if (ix == nx .and. boundxy == 2)  ixp1 = 1

            res(ix,iy,iz) = rho(ix,iy,iz)
     &          +  (phi(ixm1,iy  ,iz  )+phi(ixp1,iy  ,iz  ))*dxsqic
     &          +  (phi(ix  ,iym1,iz  )+phi(ix  ,iyp1,iz  ))*dysqic
     &          +  (phi(ix  ,iy  ,izm1)+phi(ix  ,iy  ,izp1))*dzsqic
     &          -  phi(ix,iy,iz)*mgparam

          enddo
        enddo
      enddo


c     --- Zero the residual inside conductors.
      call cond_potmg(nx,ny,nz,res,mglevel,.true.,mgform,.false.)

      if (lcndbndy) then
c     --- Calculate the residual near the conductor.
      call condbndyres(nx,ny,nz,phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,
     &                ncndmax,necndbdy,ecdelmx,ecdelmy,ecdelmz,ecdelpx,ecdelpy,
     &                ecdelpz,iecndx,iecndy,iecndz,iecndlevel,
     &                ecvoltmx,ecvoltpx,ecvoltmy,ecvoltpy,ecvoltmz,ecvoltpz,
     &                ecndpvph,localb0,localbnz,boundxy,
     &                l2symtry,l4symtry,mglevel,ecstart,
     &                mgform,mgform2init)
      call condbndyres(nx,ny,nz,phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,
     &                ncndmax,nocndbdy,ocdelmx,ocdelmy,ocdelmz,ocdelpx,ocdelpy,
     &                ocdelpz,iocndx,iocndy,iocndz,iocndlevel,
     &                ocvoltmx,ocvoltpx,ocvoltmy,ocvoltpy,ocvoltmz,ocvoltpz,
     &                ocndpvph,localb0,localbnz,boundxy,
     &                l2symtry,l4symtry,mglevel,ocstart,
     &                mgform,mgform2init)

      endif

c     --- Transverse boundaries
      if (boundxy == 2) then
        res(:,ny,:) = res(:,0,:)
        res(nx,:,:) = res(0,:,:)
      endif

c     --- Longitudinal boundaries
      if (localb0  == 0) res(:,:,-1:0) = 0.
      if (localbnz == 0) res(:,:,nz:nz+1) = 0.
      if (localb0  == 1) res(:,:,-1) = res(:,:,1)
      if (localbnz == 1) res(:,:,nz+1) = res(:,:,nz-1)
#ifndef MPIPARALLEL
      if (localb0  == 2) res(:,:,-1) = res(:,:,nz-1)
      if (localbnz == 2) res(:,:,nz+1) = res(:,:,1)
#else
      call mgexchange_phi(nx,ny,nz,nzfull,res,localb0,localbnz,-1,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
      call mgexchange_phiperiodic(nx,ny,nz,nzfull,res,localb0,localbnz,0,
     &                            my_index,nslaves,izfsslave,
     &                            whosendingleft,whosendingright)
#endif

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine setcndtr3dmg(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                        l2symtry,l4symtry)
      use PSOR3d
      use Multigrid3d
      use Conductor3d
      use Parallel
      real(kind=8):: xmmin,ymmin,zmmin,zbeam,zgrid,dx,dy,dz
      integer(ISZ):: nx,ny,nz
      logical(ISZ):: l2symtry,l4symtry

c Makes call to the recursive routine which generates the conductors.

      integer(ISZ):: nzfull
      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: zmmin1,dx1,dy1,dz1
      integer(ISZ):: nzlocal,nx1,ny1,nz1,nzfull1
      integer(ISZ):: globalb0,globalbnz
      integer(ISZ):: localb0,localbnz

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
      call setmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz)

      do il=0,mglevels-1

c       --- Set the grid parameters for this level
        dx1 = dx*mglevelslx(il)
        dy1 = dy*mglevelsly(il)
        dz1 = dz*mglevelslz(il)
        nx1 = mglevelsnx(il)
        ny1 = mglevelsny(il)
        nz1 = mglevelsnz(il)
        nzfull1 = mglevelsnzfull(il)
        zmmin1 = zmminglobal + dz1*mglevelsiz(il)
        localb0 = globalb0
        localbnz = globalbnz
#ifdef MPIPARALLEL
        if (mglevelsiz(il) > 0) localb0 = -1
        if (mglevelsiz(il)+nz1 < nzfull1) localbnz = -1
#endif

c       --- Get the conductors for this level
        ic1 = ncond + 1
        ie1 = necndbdy + 1
        io1 = nocndbdy + 1
        bound0 = localb0
        boundnz = localbnz
        call setcndtr3d(xmmin,ymmin,zmmin1,zbeam,zgrid,nx1,ny1,nz1,dx1,dy1,dz1,
     &                  l2symtry,l4symtry)
        ic2 = ncond
        ie2 = necndbdy
        io2 = nocndbdy
        icondlevel(ic1:ic2) = il
        iecndlevel(ie1:ie2) = il
        iocndlevel(io1:io2) = il

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
     &                        condid)
      use Multigrid3d
      use Conductor3d
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
      call setmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz)

      do il=0,mglevels-1

c       --- Set the grid parameters for this level
        dx1 = dx*mglevelslx(il)
        dy1 = dy*mglevelsly(il)
        dz1 = dz*mglevelslz(il)
        nx1 = mglevelsnx(il)
        ny1 = mglevelsny(il)
        nz1 = mglevelsnz(il)
        zmmin1 = zmminglobal + dz1*mglevelsiz(il)
        zmmax1 = zmmin1 + dz1*nz1
        ix_axis1 = ix_axis/mglevelslx(il)
        iy_axis1 = iy_axis/mglevelsly(il)

c       --- Get the conductors for this level
        ic1 = ncond + 1
        ie1 = necndbdy + 1
        io1 = nocndbdy + 1
        call srfrvout3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmax,lfill,
     &                  xmin,xmax,ymin,ymax,lshell,
     &                  zmmin1,zmmax1,zbeam,dx1,dy1,dz1,nx1,ny1,nz1,
     &                  ix_axis1,iy_axis1,xmesh,ymesh,l2symtry,l4symtry,condid)
        ic2 = ncond
        ie2 = necndbdy
        io2 = nocndbdy
        icondlevel(ic1:ic2) = il
        iecndlevel(ie1:ie2) = il
        iocndlevel(io1:io2) = il

      enddo

      return
      end
c=============================================================================
      subroutine srfrvin3dmg(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                       xmin,xmax,ymin,ymax,lshell,
     &                       zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                       ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,
     &                       condid)
      use Multigrid3d
      use Conductor3d
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
      call setmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz)

      do il=0,mglevels-1

c       --- Set the grid parameters for this level
        dx1 = dx*mglevelslx(il)
        dy1 = dy*mglevelsly(il)
        dz1 = dz*mglevelslz(il)
        nx1 = mglevelsnx(il)
        ny1 = mglevelsny(il)
        nz1 = mglevelsnz(il)
        zmmin1 = zmminglobal + dz1*mglevelsiz(il)
        zmmax1 = zmmin1 + dz1*nz1
        ix_axis1 = ix_axis/mglevelslx(il)
        iy_axis1 = iy_axis/mglevelsly(il)

c       --- Get the conductors for this level
        ic1 = ncond + 1
        ie1 = necndbdy + 1
        io1 = nocndbdy + 1
        call srfrvin3d(rofzfunc,volt,zmin,zmax,xcent,ycent,rmin,lfill,
     &                  xmin,xmax,ymin,ymax,lshell,
     &                  zmmin1,zmmax1,zbeam,dx1,dy1,dz1,nx1,ny1,nz1,
     &                  ix_axis1,iy_axis1,xmesh,ymesh,l2symtry,l4symtry,condid)
        ic2 = ncond
        ie2 = necndbdy
        io2 = nocndbdy
        icondlevel(ic1:ic2) = il
        iecndlevel(ie1:ie2) = il
        iocndlevel(io1:io2) = il

      enddo

      return
      end
c=============================================================================
      subroutine srfrvinout3dmg(rminofz,rmaxofz,volt,zmin,zmax,
     &                          xcent,ycent,lzend,xmin,xmax,ymin,ymax,lshell,
     &                          zmmin,zmmax,zbeam,dx,dy,dz,nx,ny,nz,
     &                          ix_axis,iy_axis,xmesh,ymesh,l2symtry,l4symtry,
     &                          condid)
      use Multigrid3d
      use Conductor3d
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
      call setmglevels(nx,ny,nzlocal,nzfull,dx,dy,dz)

      do il=0,mglevels-1

c       --- Set the grid parameters for this level
        dx1 = dx*mglevelslx(il)
        dy1 = dy*mglevelsly(il)
        dz1 = dz*mglevelslz(il)
        nx1 = mglevelsnx(il)
        ny1 = mglevelsny(il)
        nz1 = mglevelsnz(il)
        zmmin1 = zmminglobal + dz1*mglevelsiz(il)
        zmmax1 = zmmin1 + dz1*nz1
        ix_axis1 = ix_axis/mglevelslx(il)
        iy_axis1 = iy_axis/mglevelsly(il)

c       --- Get the conductors for this level
        ic1 = ncond + 1
        ie1 = necndbdy + 1
        io1 = nocndbdy + 1
        call srfrvinout3d(rminofz,rmaxofz,volt,zmin,zmax,
     &                    xcent,ycent,lzend,xmin,xmax,ymin,ymax,lshell,
     &                    zmmin1,zmmax1,zbeam,dx1,dy1,dz1,nx1,ny1,nz1,
     &                    ix_axis1,iy_axis1,xmesh,ymesh,l2symtry,l4symtry,
     &                    condid)
        ic2 = ncond
        ie2 = necndbdy
        io2 = nocndbdy
        icondlevel(ic1:ic2) = il
        iecndlevel(ie1:ie2) = il
        iocndlevel(io1:io2) = il

      enddo

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      RECURSIVE subroutine checkconductors(nx,ny,nz,nzfull,dx,dy,dz,
     &                                     l2symtry,l4symtry)
      use Multigrid3d
      use Conductor3d
      use Parallel
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

c Recursively calls the routine to generate the conductor data at the
c various mesh resolutions needed by the MG solver.

      integer(ISZ):: ic,ix,iy,iz,ii,nn
      integer(ISZ),allocatable:: isort(:)

      real(kind=8):: zmminglobal
      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2
      real(kind=8):: dx1,dy1,dz1
      integer(ISZ):: nx1,ny1,nz1
      integer(ISZ):: ix_axis1,iy_axis1

c     --- Call the work routine for each level
      do il=0,mglevels-1
        dx1 = dx*mglevelslx(il)
        dy1 = dy*mglevelsly(il)
        dz1 = dz*mglevelslz(il)
        nx1 = mglevelsnx(il)
        ny1 = mglevelsny(il)
        nz1 = mglevelsnz(il)
        call checkconductors_work(nx1,ny1,nz1,dx1,dy1,dz1,
     &                            l2symtry,l4symtry,il)
      enddo

c     --- Sort the conductor data by level number
c     --- First, sort conductor points
      if (ncond > 0) then
        allocate(isort(ncond))
        call isortconductor(ncond,icondlevel,isort,icstart)
        nn = ncond
        ncond = maxval(icstart) - 1
        call iswapconductor(ncond,isort,nn,ixcond)
        call iswapconductor(ncond,isort,nn,iycond)
        call iswapconductor(ncond,isort,nn,izcond)
        call iswapconductor(ncond,isort,nn,condnumb)
        call rswapconductor(ncond,isort,nn,condvolt)
        call iswapconductor(ncond,isort,nn,icondlevel)
        deallocate(isort)
      endif

c     --- Sort even subgrid points
      if (necndbdy > 0) then
        allocate(isort(necndbdy))
        call isortconductor(necndbdy,iecndlevel,isort,ecstart)
        nn = necndbdy
        necndbdy = maxval(ecstart) - 1
        call iswapconductor(necndbdy,isort,nn,iecndx)
        call iswapconductor(necndbdy,isort,nn,iecndy)
        call iswapconductor(necndbdy,isort,nn,iecndz)
        call rswapconductor(necndbdy,isort,nn,ecdelmx)
        call rswapconductor(necndbdy,isort,nn,ecdelmy)
        call rswapconductor(necndbdy,isort,nn,ecdelmz)
        call rswapconductor(necndbdy,isort,nn,ecdelpx)
        call rswapconductor(necndbdy,isort,nn,ecdelpy)
        call rswapconductor(necndbdy,isort,nn,ecdelpz)
        call rswapconductor(necndbdy,isort,nn,ecvolt)
        call iswapconductor(necndbdy,isort,nn,ecnumb)
        call rswapconductor(necndbdy,isort,nn,ecvoltmx)
        call rswapconductor(necndbdy,isort,nn,ecvoltpx)
        call rswapconductor(necndbdy,isort,nn,ecvoltmy)
        call rswapconductor(necndbdy,isort,nn,ecvoltpy)
        call rswapconductor(necndbdy,isort,nn,ecvoltmz)
        call rswapconductor(necndbdy,isort,nn,ecvoltpz)
        call iswapconductor(necndbdy,isort,nn,ecnumbmx)
        call iswapconductor(necndbdy,isort,nn,ecnumbpx)
        call iswapconductor(necndbdy,isort,nn,ecnumbmy)
        call iswapconductor(necndbdy,isort,nn,ecnumbpy)
        call iswapconductor(necndbdy,isort,nn,ecnumbmz)
        call iswapconductor(necndbdy,isort,nn,ecnumbpz)
        call iswapconductor(necndbdy,isort,nn,iecndlevel)
        deallocate(isort)
      endif
    
c     --- Sort odd subgrid points
      if (nocndbdy > 0) then
        allocate(isort(nocndbdy))
        call isortconductor(nocndbdy,iocndlevel,isort,ocstart)
        nn = nocndbdy
        nocndbdy = maxval(ocstart) - 1
        call iswapconductor(nocndbdy,isort,nn,iocndx)
        call iswapconductor(nocndbdy,isort,nn,iocndy)
        call iswapconductor(nocndbdy,isort,nn,iocndz)
        call rswapconductor(nocndbdy,isort,nn,ocdelmx)
        call rswapconductor(nocndbdy,isort,nn,ocdelmy)
        call rswapconductor(nocndbdy,isort,nn,ocdelmz)
        call rswapconductor(nocndbdy,isort,nn,ocdelpx)
        call rswapconductor(nocndbdy,isort,nn,ocdelpy)
        call rswapconductor(nocndbdy,isort,nn,ocdelpz)
        call rswapconductor(nocndbdy,isort,nn,ocvolt)
        call iswapconductor(nocndbdy,isort,nn,ocnumb)
        call rswapconductor(nocndbdy,isort,nn,ocvoltmx)
        call rswapconductor(nocndbdy,isort,nn,ocvoltpx)
        call rswapconductor(nocndbdy,isort,nn,ocvoltmy)
        call rswapconductor(nocndbdy,isort,nn,ocvoltpy)
        call rswapconductor(nocndbdy,isort,nn,ocvoltmz)
        call rswapconductor(nocndbdy,isort,nn,ocvoltpz)
        call iswapconductor(nocndbdy,isort,nn,ocnumbmx)
        call iswapconductor(nocndbdy,isort,nn,ocnumbpx)
        call iswapconductor(nocndbdy,isort,nn,ocnumbmy)
        call iswapconductor(nocndbdy,isort,nn,ocnumbpy)
        call iswapconductor(nocndbdy,isort,nn,ocnumbmz)
        call iswapconductor(nocndbdy,isort,nn,ocnumbpz)
        call iswapconductor(nocndbdy,isort,nn,iocndlevel)
        deallocate(isort)
      endif

      return
      end
c=============================================================================
      subroutine isortconductor(nc,condlevel,isort,icstart)
      integer(ISZ):: nc
      integer(ISZ):: condlevel(nc),isort(nc),icstart(0:100)

c     --- Sort the conductor data by level number
      integer(ISZ):: ilevel,nlevel,ii,ic

      icstart = 0
      ilevel = 0
      nlevel = 1
      ii = 0
      do while (nlevel > 0)
        icstart(ilevel) = ii + 1
        nlevel = ii
        do ic=1,nc
          if (condlevel(ic) == ilevel) then
            ii = ii + 1
            isort(ii) = ic
          endif
        enddo
        nlevel = ii - nlevel
        ilevel = ilevel + 1
      enddo

      return
      end
c=============================================================================
      subroutine iswapconductor(nc,isort,nn,icond)
      integer(ISZ):: nc,nn,isort(nc),icond(nn)
c Utility function for sort conductors
      integer(ISZ):: ic,itemp(nn)
      itemp = icond
      do ic=1,nc
        icond(ic) = itemp(isort(ic))
      enddo
      return
      end
c=============================================================================
      subroutine rswapconductor(nc,isort,nn,rcond)
      integer(ISZ):: nc,nn,isort(nc)
      real(kind=8):: rcond(nn)
c Utility function for sort conductors
      real(kind=8):: rtemp(nn)
      integer(ISZ):: ic
      rtemp = rcond
      do ic=1,nc
        rcond(ic) = rtemp(isort(ic))
      enddo
      return
      end
c=============================================================================
      subroutine checkconductors_work(nx,ny,nz,dx,dy,dz,l2symtry,l4symtry,
     &                                mglevel)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny,nz,mglevel
      real(kind=8):: dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

c This checks the conductor dataset for consistency.
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

      integer(ISZ):: iii(0:nx,0:ny,0:nz)
      integer(ISZ):: ic,i,ix,iy,iz

c     --- Set the conductor points.
      iii = 0
      do ic=1,ncond
        if (icondlevel(ic) /= mglevel) cycle
        iii(ixcond(ic),iycond(ic),izcond(ic)) = ncndmax + 1
      enddo

c     --- Scan through subgrid points:
c     ---   remove points which lie inside of a conductor
c     ---   register subgrid points in the work array iii
c     ---   check for redundant point (multiple points at grid location)
      do ic=1,necndbdy

        if (iecndlevel(ic) /= mglevel) cycle

        ix = iecndx(ic)
        iy = iecndy(ic)
        iz = iecndz(ic)

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy,iz) == ncndmax+1) then
          iecndlevel(ic) = -1
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

        if (iecndlevel(ic) /= iecndlevel(i)) then
c         --- The point already there is on a different level
          iii(ix,iy,iz) = ic
          cycle
        endif

c       --- There is another subgrid point here. Combine the data
c       --- of the two points.
c       --- For each direction, check if a conductor is nearer to this point
c       --- than the other point. If so, reset data for the other point.
        if (ecdelmx(ic) < 1. .and. ecdelmx(ic) < ecdelmx(i)) then
          ecdelmx(i)  = ecdelmx(ic)
          ecvoltmx(i) = ecvoltmx(ic)
          ecnumbmx(i) = ecnumbmx(ic)
        endif
        if (ecdelpx(ic) < 1. .and. ecdelpx(ic) < ecdelpx(i)) then
          ecdelpx(i)  = ecdelpx(ic)
          ecvoltpx(i) = ecvoltpx(ic)
          ecnumbpx(i) = ecnumbpx(ic)
        endif
        if (ecdelmy(ic) < 1. .and. ecdelmy(ic) < ecdelmy(i)) then
          ecdelmy(i)  = ecdelmy(ic)
          ecvoltmy(i) = ecvoltmy(ic)
          ecnumbmy(i) = ecnumbmy(ic)
        endif
        if (ecdelpy(ic) < 1. .and. ecdelpy(ic) < ecdelpy(i)) then
          ecdelpy(i)  = ecdelpy(ic)
          ecvoltpy(i) = ecvoltpy(ic)
          ecnumbpy(i) = ecnumbpy(ic)
        endif
        if (ecdelmz(ic) < 1. .and. ecdelmz(ic) < ecdelmz(i)) then
          ecdelmz(i)  = ecdelmz(ic)
          ecvoltmz(i) = ecvoltmz(ic)
          ecnumbmz(i) = ecnumbmz(ic)
        endif
        if (ecdelpz(ic) < 1. .and. ecdelpz(ic) < ecdelpz(i)) then
          ecdelpz(i)  = ecdelpz(ic)
          ecvoltpz(i) = ecvoltpz(ic)
          ecnumbpz(i) = ecnumbpz(ic)
        endif
        iecndlevel(ic) = -1

      enddo

c     --- Do the same for the odd conductor points.
c     --- Scan through subgrid points:
c     ---   remove points which lie inside of a conductor
c     ---   register subgrid points in the work array iii
c     ---   check for redundant point (multiple points at grid location)
      do ic=1,nocndbdy

        if (iocndlevel(ic) /= mglevel) cycle

        ix = iocndx(ic)
        iy = iocndy(ic)
        iz = iocndz(ic)

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy,iz) == ncndmax+1) then
          iocndlevel(ic) = -1
          cycle
        endif

c       --- If iii == 0, then this data point is outside any conductors and
c       --- is not redundant.
        if (iii(ix,iy,iz) == 0) then
          iii(ix,iy,iz) = -ic
          cycle
        endif

c       iii(ix,iy,iz) < ncndmax+1
        i = -iii(ix,iy,iz)

        if (i < 0) then
c         --- The point already there is even so must be a different level
          iii(ix,iy,iz) = -ic
          cycle
        endif

        if (iocndlevel(ic) /= iocndlevel(i)) then
c         --- The point is on a different level
          iii(ix,iy,iz) = -ic
          cycle
        endif

c       --- There is another subgrid point here. Combine the data
c       --- of the two points.
c       --- For each direction, check if a conductor is nearer to this point
c       --- than the other point. If so, reset data for the other point.
        if (ocdelmx(ic) < 1. .and. ocdelmx(ic) < ocdelmx(i)) then
          ocdelmx(i)  = ocdelmx(ic)
          ocvoltmx(i) = ocvoltmx(ic)
          ocnumbmx(i) = ocnumbmx(ic)
        endif
        if (ocdelpx(ic) < 1. .and. ocdelpx(ic) < ocdelpx(i)) then
          ocdelpx(i)  = ocdelpx(ic)
          ocvoltpx(i) = ocvoltpx(ic)
          ocnumbpx(i) = ocnumbpx(ic)
        endif
        if (ocdelmy(ic) < 1. .and. ocdelmy(ic) < ocdelmy(i)) then
          ocdelmy(i)  = ocdelmy(ic)
          ocvoltmy(i) = ocvoltmy(ic)
          ocnumbmy(i) = ocnumbmy(ic)
        endif
        if (ocdelpy(ic) < 1. .and. ocdelpy(ic) < ocdelpy(i)) then
          ocdelpy(i)  = ocdelpy(ic)
          ocvoltpy(i) = ocvoltpy(ic)
          ocnumbpy(i) = ocnumbpy(ic)
        endif
        if (ocdelmz(ic) < 1. .and. ocdelmz(ic) < ocdelmz(i)) then
          ocdelmz(i)  = ocdelmz(ic)
          ocvoltmz(i) = ocvoltmz(ic)
          ocnumbmz(i) = ocnumbmz(ic)
        endif
        if (ocdelpz(ic) < 1. .and. ocdelpz(ic) < ocdelpz(i)) then
          ocdelpz(i)  = ocdelpz(ic)
          ocvoltpz(i) = ocvoltpz(ic)
          ocnumbpz(i) = ocnumbpz(ic)
        endif
        iocndlevel(ic) = -1
      enddo

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine subgrid_sor_to_mg(nx,ny,nz,dx,dy,dz,l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry

c This routine now is not needed.

      return
      end
c=============================================================================
