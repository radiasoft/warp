#include "top.h"
c=============================================================================
c@(#) File F3D_MGRID.M, version $Revision: 3.4 $, $Date: 2001/05/03 19:08:47 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid field sovler which is part of the F3D
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigrid3df(iwhich,nx,ny,nz,dx,dy,dz,phi1d,rho1d,
     &                       rstar,linbends,l2symtry,l4symtry,scrtch,
     &                       xmmin,ymmin,zmmin,zbeam,zgrid)
      use Constant
      use PSOR3d
      use Multigrid3d
      use Multigrid3d_res
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi1d(*)
      real(kind=8):: rho1d(*)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nz+1)
      logical(ISZ):: linbends
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
c All dimensions, nx, ny, and nz, must be powers of 2.



      integer(ISZ):: nxy,nxyz,nxyz3
      integer(ISZ):: log2nx,log2ny,log2nz
      real(kind=8):: maxerr
      real(kind=8):: dxsqi,dysqi,dzsqi,reps0c,rdel
      real(kind=8):: tdx,tdy,tdz
      integer(ISZ):: tnx,tny,tnz
      integer(ISZ):: i,ii,minlog,k
      real(kind=8):: rs,x,r

c     --- Make sure that all dimensions are powers of 2.
      log2nx = log(real(nx))/log(2.)+0.5
      log2ny = log(real(ny))/log(2.)+0.5
      log2nz = log(real(nz))/log(2.)+0.5
c     if (2**log2nx .ne. nx .or. 2**log2ny .ne. ny .or. 2**log2nz .ne. nz) then
c       print*,"Multigrid: nx, ny, and nz must be powers of 2"
c       return
c     endif

c     --- If iwhich equals zero or one, do initialization.
c     --- psorinit does the same things that mginit would do
c     --- except that multigrid does not use sorerrar.
      if (iwhich .eq. 0 .or. iwhich .eq. 1)
     &  call psorinit(nx,ny,nz,dx,dy,dz,l2symtry,l4symtry)

c     --- Initialize temporaries
      nxy  =(nx+1)*(ny+1)
      nxyz =(nx+1)*(ny+1)*(nz+1)
      nxyz3=(nx+1)*(ny+1)*(nz+3)
      dxsqi=1./dx**2
      dysqi=1./dy**2
      dzsqi=1./dz**2
      reps0c = mgparam/(eps0*2.*(dxsqi+dysqi+dzsqi))
      rdel = dzsqi/(dxsqi + dysqi + dzsqi)
      minlog = min(log2nx,log2ny,log2nz)

c     --- Allocate the res array.
      if (iwhich .eq. 0 .or. iwhich .eq. 1) then
        res_size = nxyz
        call gchange("Multigrid3d_res",0)
      endif

c     --- Calculate the size of the work arrays. The work arrays must be
c     --- big enough to hold the phi and the rho at each level of coarseness,
c     --- including the finest. The size calculation include semi-coarsening.
c     --- Only do the allocation during initialization.
c     --- This loops over the grid levels, adding in the size of the grid
c     --- at each level. Note that the size is calculated for the phi array,
c     --- which has extra planes, so that the size if a slight over-estimate
c     --- for rho.
      if (iwhich .eq. 0 .or. iwhich .eq. 1) then
        tnx = nx
        tny = ny
        tnz = nz
        tdx = dx
        tdy = dy
        tdz = dz
        tempsize = nxyz3
        do while (tnx > 2 .and. tny > 2 .and. tnz > 2 .and.
     &            mod(tnx,4) == 0 .and. mod(tny,4) == 0 .and. mod(tnz,4) == 0)
          tnx = tnx/2
          tny = tny/2
          tdx = tdx*2.
          tdy = tdy*2.
          if (1.5*tdz < tdx) then
            tnz = tnz/2
            tdz = tdz*2.
          endif
          tempsize = tempsize + (tnx+1)*(tny+1)*(tnz+3)
        end do
        call gallot("Multigrid3d",0)
      endif

c     --- If doing initialization only, then exit.
      if (iwhich .eq. 1) return

c     --- Determine the points that make up the conductor.  This takes extra
c     --- time and so should not be done if the grid is not moving in the lab
c     --- frame.  Set gridmode to 1 to avoid this call. The data is then
c     --- converted and expanded for the multigrid solver.
      if (gridmode .eq. 0 .or. iwhich .eq. -2) then
        call setcndtr(xmmin,ymmin,zmmin,zbeam,zgrid,nx,ny,nz,dx,dy,dz,
     &                l2symtry,l4symtry)
        call subgrid_sor_to_mg(nx,ny,nz,dx,dy,dz)
      endif

c     --- Copy phi to the temporary array and zero out rho_temp.
      call copyarry(phi1d,phi_temp,nxyz3)
c     call zeroarry(rho_temp(nxyz3+1),tempsize-nxyz3)
      rho_temp(nxyz3+1:) = 0.

c     --- Preset rho to increase performance (reducing the number of
c     --- multiplies in the main SOR sweep loop). Also, copy rho1d
c     --- into rho_temp.
      if (.not. linbends) then
        do ii=1,nxyz
          rho_temp(ii) = rho1d(ii)*reps0c
        enddo
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
c         --- rho_temp(ii) = rho1d(ii)*(rs/r)*reps0c/
c         ---             ( 1. + (x/r)*((x/r)-2.)*rdel )
          rho_temp(ii) = rho1d(ii)*reps0c*rs*r/(r*r + x*(x-2.*r)*rdel)
        enddo
c       --- Fill scratch array with x values so it can be looked up
c       --- in the bent beam loop instead of calculated.
        do ii = 1, nxy
          scrtch(ii) = xmmin + mod(ii-1,nx+1)*dx
        enddo
c       --- Change rstar if using Nuemann boundary conditions
        if (bound0 .eq. 1) rstar(-1) = rstar(1)
        if (boundnz .eq. 1) rstar(nz+1) = rstar(nz-1)
      endif

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- vary few iterations are done.
      mgiters = 0
      maxerr = 2.*mgtol + 1.
      do while (maxerr > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Zero out the coarse grid data area.
c       call zeroarry(phi_temp(nxyz3+1),tempsize-nxyz3)
        phi_temp(nxyz3+1:) = 0.

c       --- Do one vcycle.
        call vcycle(1,1,nx,ny,nz,dxsqi,dysqi,dzsqi,phi_temp,rho_temp,res,
     &              rstar,linbends,l2symtry,l4symtry,scrtch,
     &              boundxy,bound0,boundnz)

c       --- Calculate the change in phi.
        maxerr = 0.
        do i=1+nxy,nxy+nxyz
          maxerr = max(maxerr,abs(phi_temp(i) - phi1d(i)))
        enddo

c       --- Save the current phi. The original phi1d array is used as the
c       --- scratch space to hold the current value of phi which is used
c       --- after the iteration to check for convergence.
        call copyarry(phi_temp(1+nxy),phi1d(1+nxy),nxyz)

      end do

c     --- Set planes iz=-1 and nz+1 according to the axial boundary conditions
      if (bound0 .eq. 0) then
        call copyarry(phi1d(1+   nxy),phi1d(1),nxy)
      else if (bound0 .eq. 1) then
        call copyarry(phi1d(1+ 2*nxy),phi1d(1),nxy)
      else if (bound0 .eq. 2) then
        call copyarry(phi1d(1+nz*nxy),phi1d(1),nxy)
      endif
      if (boundnz .eq. 0) then
        call copyarry(phi1d(1+  nxyz),phi1d(1+nxy+nxyz),nxy)
      else if (boundnz .eq. 1) then
        call copyarry(phi1d(1+nz*nxy),phi1d(1+nxy+nxyz),nxy)
      else if (boundnz .eq. 2) then
        call copyarry(phi1d(1+ 2*nxy),phi1d(1+nxy+nxyz),nxy)
      endif

c     --- Make a print out.
      if (maxerr > mgtol) then
        print*,"Multigrid: Maximum number of iterations reached"
      endif
      print*,"Multigrid: Error converged to ",maxerr," in ",mgiters," v-cycles"

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycle(levelxy,levelz,nx,ny,nz,dxsqi,dysqi,dzsqi,
     &                            phi,rho,res,rstar,linbends,l2symtry,l4symtry,
     &                            scrtch,boundxy,bound0,boundnz)
      use Multigrid3d
      integer(ISZ):: levelxy,levelz
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(*),rho(*),res(*),rstar(*),scrtch(*)
      logical(ISZ):: linbends
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: boundxy,bound0,boundnz

c Routine that does the v-cycle for multigrid. Note that it is recursive.


      integer(ISZ):: i2h,i
      real(kind=8):: ff

c     --- Do initial SOR passes.
      do i=1,downpasses
        call sorpass3d(levelxy,levelz,nx,ny,nz,phi,rho,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbends,l2symtry,l4symtry,scrtch)
      enddo

c     --- If grid is not at its coarsest level in any of the axis or and
c     --- all dimensions are even, continue the coarsening.
      if (nx > 2 .and. ny > 2 .and. nz > 2 .and.
     &    mod(nx,4) == 0 .and. mod(ny,4) == 0 .and. mod(nz,4) == 0) then

c       --- Location of the next block of data.
        i2h = (nx+1)*(ny+1)*(nz+3) + 1

c       --- Get the residual on the current grid.
        call residual(nx,ny,nz,dxsqi,dysqi,dzsqi,phi,rho,res,levelxy,levelz,
     &                l2symtry,l4symtry)

c       --- If dz > 4/3 dx then only coarsen transversely, otherwise coarsen
c       --- all axis.
        if (0.5625*dxsqi > dzsqi) then

c         --- Ratio of old to new constant needed to scale the residual for
c         --- the restriction.
          ff = (dxsqi+dysqi+dzsqi)/(dxsqi*0.25 + dysqi*0.25 + dzsqi)
          call restrict2d(nx,ny,nz,res,rho(i2h),ff,boundxy,bound0,boundnz,
     &                    l2symtry,l4symtry)

c         --- Continue at the next coarsest level.
          call vcycle(levelxy*2,levelz,nx/2,ny/2,nz,
     &                dxsqi*0.25,dysqi*0.25,dzsqi,phi(i2h),rho(i2h),res,
     &                rstar,linbends,l2symtry,l4symtry,scrtch,
     &                boundxy,bound0,boundnz)

c         --- Add in resulting error.
          call expand2d(nx/2,ny/2,nz,phi(i2h),phi,boundxy,bound0,boundnz)

        else

c         --- Restiction - note that scaling factor for residual is always
c         --- 4 for full-coarsening and is compiled into the restriction
c         --- routine.
          call restrict3d(nx,ny,nz,res,rho(i2h),boundxy,bound0,boundnz,
     &                    l2symtry,l4symtry)

c         --- Continue at the next coarsest level.
          call vcycle(levelxy*2,levelz*2,nx/2,ny/2,nz/2,
     &                dxsqi*0.25,dysqi*0.25,dzsqi*0.25,phi(i2h),rho(i2h),res,
     &                rstar,linbends,l2symtry,l4symtry,scrtch,
     &                boundxy,bound0,boundnz)

c         --- Add in resulting error.
          call expand3d(nx/2,ny/2,nz/2,phi(i2h),phi,boundxy,boundnz)

        endif

      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call sorpass3d(levelxy,levelz,nx,ny,nz,phi,rho,phi,rho,rstar,
     &                 dxsqi,dysqi,dzsqi,linbends,l2symtry,l4symtry,scrtch)
      enddo

      return
      end
c=============================================================================
      subroutine restrict2d(nx,ny,nz,ph,p2h,ff,boundxy,bound0,boundnz,
     &                      l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: ph(0:nx,0:ny,0:nz)
      real(kind=8):: p2h(0:nx/2,0:ny/2,0:nz)
      real(kind=8):: ff
      integer(ISZ):: boundxy,bound0,boundnz
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
      izmin = 1
      izmax = nz-1
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny
      if (bound0 > 0) izmin = 0
      if (boundnz == 1) izmax = nz

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

            p2h(ix/2,iy/2,iz) =
     &        ff*0.2500*ph(ix  ,iy  ,iz) +
     &        ff*0.1250*ph(ixm1,iy  ,iz) +
     &        ff*0.1250*ph(ixp1,iy  ,iz) +
     &        ff*0.1250*ph(ix  ,iym1,iz) +
     &        ff*0.1250*ph(ix  ,iyp1,iz) +
     &        ff*0.0625*ph(ixm1,iym1,iz) +
     &        ff*0.0625*ph(ixp1,iym1,iz) +
     &        ff*0.0625*ph(ixm1,iyp1,iz) +
     &        ff*0.0625*ph(ixp1,iyp1,iz)

          enddo
        enddo
      enddo

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        do iz=izmin/2,izmax/2
          do ix=ixmin/2,ixmax/2
            p2h(ix,ny/2,iz) = p2h(ix,0,iz)
          enddo
          do iy=iymin/2,iymax/2
            p2h(nx/2,iy,iz) = p2h(0,iy,iz)
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine restrict3d(nx,ny,nz,ph,p2h,boundxy,bound0,boundnz,
     &                      l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: ph(0:nx,0:ny,0:nz)
      real(kind=8):: p2h(0:nx/2,0:ny/2,0:nz/2)
      integer(ISZ):: boundxy,bound0,boundnz
      logical(ISZ):: l2symtry,l4symtry

c Restrict to a coarser grid.  The factor of 4 is needed since the
c residual will be used as the source term in the next coarser grid and
c it implicitly includes the factor 0.5/(dxsqi+dysqi+dzsqi) which is 4
c times too small with the current values of dxsqi etc.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 2
      ixmax = nx-2
      iymin = 2
      iymax = ny-2
      izmin = 2
      izmax = nz-2
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny
      if (bound0 > 0) izmin = 0
      if (boundnz == 1) izmax = nz

c     --- Do the loops.
      do iz=izmin,izmax,2
        izm1 = iz - 1
        if (iz == 0 .and. bound0 == 1)  izm1 = 1
        if (iz == 0 .and. bound0 == 2)  izm1 = nz-1
        izp1 = iz + 1
        if (iz == nz .and. boundnz == 1)  izp1 = nz-1
        if (iz == nz .and. boundnz == 2)  izp1 = 1

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

            p2h(ix/2,iy/2,iz/2) =
     &        4.0*0.125000*ph(ix  ,iy  ,iz  ) +
     &        4.0*0.062500*ph(ixm1,iy  ,iz  ) +
     &        4.0*0.062500*ph(ixp1,iy  ,iz  ) +
     &        4.0*0.062500*ph(ix  ,iym1,iz  ) +
     &        4.0*0.062500*ph(ix  ,iyp1,iz  ) +
     &        4.0*0.062500*ph(ix  ,iy  ,izm1) +
     &        4.0*0.062500*ph(ix  ,iy  ,izp1) +
     &        4.0*0.031250*ph(ixm1,iym1,iz  ) +
     &        4.0*0.031250*ph(ixp1,iym1,iz  ) +
     &        4.0*0.031250*ph(ixm1,iyp1,iz  ) +
     &        4.0*0.031250*ph(ixp1,iyp1,iz  ) +
     &        4.0*0.031250*ph(ixm1,iy  ,izm1) +
     &        4.0*0.031250*ph(ixp1,iy  ,izm1) +
     &        4.0*0.031250*ph(ixm1,iy  ,izp1) +
     &        4.0*0.031250*ph(ixp1,iy  ,izp1) +
     &        4.0*0.031250*ph(ix  ,iym1,izm1) +
     &        4.0*0.031250*ph(ix  ,iyp1,izm1) +
     &        4.0*0.031250*ph(ix  ,iym1,izp1) +
     &        4.0*0.031250*ph(ix  ,iyp1,izp1) +
     &        4.0*0.015625*ph(ixm1,iym1,izm1) +
     &        4.0*0.015625*ph(ixp1,iym1,izm1) +
     &        4.0*0.015625*ph(ixm1,iyp1,izm1) +
     &        4.0*0.015625*ph(ixp1,iyp1,izm1) +
     &        4.0*0.015625*ph(ixm1,iym1,izp1) +
     &        4.0*0.015625*ph(ixp1,iym1,izp1) +
     &        4.0*0.015625*ph(ixm1,iyp1,izp1) +
     &        4.0*0.015625*ph(ixp1,iyp1,izp1)

          enddo
        enddo
      enddo

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        do iz=izmin/2,izmax/2
          do ix=ixmin/2,ixmax/2
            p2h(ix,ny/2,iz) = p2h(ix,0,iz)
          enddo
          do iy=iymin/2,iymax/2
            p2h(nx/2,iy,iz) = p2h(0,iy,iz)
          enddo
        enddo
      endif
      if (boundnz == 2) then
        do iy=iymin/2,iymax/2
          do ix=ixmin/2,ixmax/2
            p2h(ix,iy,nz/2) = p2h(ix,iy,0)
          enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine expand2d(nx,ny,nz,p2h,ph,boundxy,bound0,boundnz)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: p2h(0:nx,0:ny,-1:nz+1)
      real(kind=8):: ph(0:nx*2,0:ny*2,-1:nz+1)
      integer(ISZ):: boundxy,bound0,boundnz

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: izmin,izmax

c     --- Set the loop limits, including edges when appropriate.
      izmin = 1
      izmax = nz-1
      if (bound0 > 0) izmin = 0
      if (boundnz == 1) izmax = nz

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
      subroutine expand3d(nx,ny,nz,p2h,ph,boundxy,boundnz)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: p2h(0:nx,0:ny,-1:nz+1)
      real(kind=8):: ph(0:nx*2,0:ny*2,-1:nz*2+1)
      integer(ISZ):: boundxy,boundnz

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is done along all axis.

      integer(ISZ):: ix,iy,iz

c     --- Expand the bulk of the grid.
      do iz=0,nz-1
        do iy=0,ny-1
          do ix=0,nx-1

      ph(2*ix  ,2*iy  ,2*iz  )=ph(2*ix  ,2*iy  ,2*iz  ) +
     &                                                   p2h(ix  ,iy  ,iz  )
      ph(2*ix+1,2*iy  ,2*iz  )=ph(2*ix+1,2*iy  ,2*iz  ) +
     &                                            0.500*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix+1,iy  ,iz  ))
      ph(2*ix  ,2*iy+1,2*iz  )=ph(2*ix  ,2*iy+1,2*iz  ) +
     &                                            0.500*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix  ,iy+1,iz  ))
      ph(2*ix  ,2*iy  ,2*iz+1)=ph(2*ix  ,2*iy  ,2*iz+1) +
     &                                            0.500*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix  ,iy  ,iz+1))
      ph(2*ix+1,2*iy+1,2*iz  )=ph(2*ix+1,2*iy+1,2*iz  ) +
     &                                            0.250*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix+1,iy  ,iz  ) +
     &                                                   p2h(ix  ,iy+1,iz  ) +
     &                                                   p2h(ix+1,iy+1,iz  ))
      ph(2*ix+1,2*iy  ,2*iz+1)=ph(2*ix+1,2*iy  ,2*iz+1) +
     &                                            0.250*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix+1,iy  ,iz  ) +
     &                                                   p2h(ix  ,iy  ,iz+1) +
     &                                                   p2h(ix+1,iy  ,iz+1))
      ph(2*ix  ,2*iy+1,2*iz+1)=ph(2*ix  ,2*iy+1,2*iz+1) +
     &                                            0.250*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix  ,iy+1,iz  ) +
     &                                                   p2h(ix  ,iy  ,iz+1) +
     &                                                   p2h(ix  ,iy+1,iz+1))
      ph(2*ix+1,2*iy+1,2*iz+1)=ph(2*ix+1,2*iy+1,2*iz+1) +
     &                                            0.125*(p2h(ix  ,iy  ,iz  ) +
     &                                                   p2h(ix+1,iy  ,iz  ) +
     &                                                   p2h(ix  ,iy+1,iz  ) +
     &                                                   p2h(ix+1,iy+1,iz  ) +
     &                                                   p2h(ix  ,iy  ,iz+1) +
     &                                                   p2h(ix+1,iy  ,iz+1) +
     &                                                   p2h(ix  ,iy+1,iz+1) +
     &                                                   p2h(ix+1,iy+1,iz+1))

          enddo
        enddo
      enddo

      if (boundxy > 0) then
        do iz=0,nz-1

c         --- Expand ix=nx, iy=ny line.
          ph(2*nx,2*ny,2*iz  )=ph(2*nx,2*ny,2*iz  ) +        p2h(nx,ny,iz  )
          ph(2*nx,2*ny,2*iz+1)=ph(2*nx,2*ny,2*iz+1) + 0.500*(p2h(nx,ny,iz  ) +
     &                                                       p2h(nx,ny,iz+1))

c         --- Expand ix=nx plane.
          do iy=0,ny-1
            ph(2*nx,2*iy  ,2*iz  )=ph(2*nx,2*iy  ,2*iz  ) +
     &                                                      p2h(nx,iy  ,iz  )
            ph(2*nx,2*iy+1,2*iz  )=ph(2*nx,2*iy+1,2*iz  )+
     &                                               0.500*(p2h(nx,iy  ,iz  ) +
     &                                                      p2h(nx,iy+1,iz  ))
            ph(2*nx,2*iy  ,2*iz+1)=ph(2*nx,2*iy  ,2*iz+1)+
     &                                               0.500*(p2h(nx,iy  ,iz  ) +
     &                                                      p2h(nx,iy  ,iz+1))
            ph(2*nx,2*iy+1,2*iz+1)=ph(2*nx,2*iy+1,2*iz+1)+
     &                                               0.250*(p2h(nx,iy  ,iz  ) +
     &                                                      p2h(nx,iy+1,iz  ) +
     &                                                      p2h(nx,iy  ,iz+1) +
     &                                                      p2h(nx,iy+1,iz+1))
          enddo

c         --- Expand iy=ny plane.
          do ix=0,nx-1
            ph(2*ix  ,2*ny,2*iz  ) = ph(2*ix  ,2*ny,2*iz  ) +
     &                                                      p2h(ix  ,ny,iz  )
            ph(2*ix+1,2*ny,2*iz  ) = ph(2*ix+1,2*ny,2*iz  )+
     &                                               0.500*(p2h(ix  ,ny,iz  ) +
     &                                                      p2h(ix+1,ny,iz  ))
            ph(2*ix  ,2*ny,2*iz+1) = ph(2*ix  ,2*ny,2*iz+1)+
     &                                               0.500*(p2h(ix  ,ny,iz  ) +
     &                                                      p2h(ix  ,ny,iz+1))
            ph(2*ix+1,2*ny,2*iz+1) = ph(2*ix+1,2*ny,2*iz+1)+
     &                                               0.250*(p2h(ix  ,ny,iz  ) +
     &                                                      p2h(ix+1,ny,iz  ) +
     &                                                      p2h(ix  ,ny,iz+1) +
     &                                                      p2h(ix+1,ny,iz+1))
          enddo

        enddo
      endif

      if (boundnz > 0) then
c       --- Expand iz=nz plane.
        ph(2*nx,2*ny,2*nz) = ph(2*nx,2*ny,2*nz)+p2h(nx,ny,nz)
        do iy=0,ny-1
          ph(2*nx,2*iy  ,2*nz) = ph(2*nx,2*iy  ,2*nz)+     p2h(nx,iy  ,nz)
          ph(2*nx,2*iy+1,2*nz) = ph(2*nx,2*iy+1,2*nz)+0.5*(p2h(nx,iy  ,nz) +
     &                                                     p2h(nx,iy+1,nz))
        enddo
        do ix=0,nx-1
          ph(2*ix  ,2*ny,2*nz) = ph(2*ix  ,2*ny,2*nz) +    p2h(ix  ,ny,nz)
          ph(2*ix+1,2*ny,2*nz) = ph(2*ix+1,2*ny,2*nz)+0.5*(p2h(ix  ,ny,nz) +
     &                                                     p2h(ix+1,ny,nz))
        enddo
        do iy=0,ny-1
          do ix=0,nx-1

            ph(2*ix  ,2*iy  ,2*nz) = ph(2*ix  ,2*iy  ,2*nz) +
     &                                                      p2h(ix  ,iy  ,nz)
            ph(2*ix+1,2*iy  ,2*nz) = ph(2*ix+1,2*iy  ,2*nz)+
     &                                               0.500*(p2h(ix  ,iy  ,nz) +
     &                                                      p2h(ix+1,iy  ,nz))
            ph(2*ix  ,2*iy+1,2*nz) = ph(2*ix  ,2*iy+1,2*nz)+
     &                                               0.500*(p2h(ix  ,iy  ,nz) +
     &                                                      p2h(ix  ,iy+1,nz))
            ph(2*ix+1,2*iy+1,2*nz) = ph(2*ix+1,2*iy+1,2*nz)+
     &                                               0.250*(p2h(ix  ,iy  ,nz) +
     &                                                      p2h(ix+1,iy  ,nz) +
     &                                                      p2h(ix  ,iy+1,nz) +
     &                                                      p2h(ix+1,iy+1,nz))

          enddo
        enddo
      endif

      return
      end
c=============================================================================
c=============================================================================
      subroutine sorpass3d(levelxy,levelz,nx,ny,nz,phi,rho,phi1d,rho1d,rstar,
     &                     rdx2,rdy2,rdz2,linbends,l2symtry,l4symtry,scrtch)
      use PSOR3d
      use MultigridConductor3d
      use Multigrid3d
      use Constant
#ifdef PARALLEL
      use Parallel
#endif
      integer(ISZ):: levelxy,levelz,nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz),rstar(-1:nz+1)
      real(kind=8):: phi1d(*),rho1d(*),scrtch(*)
      real(kind=8):: rdx2,rdy2,rdz2
      logical(ISZ):: linbends
      logical(ISZ):: l2symtry,l4symtry

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


      real(kind=8):: rdel,reps0,const,rdx2c,rdy2c,rdz2c,spm1,dx
      integer(ISZ):: nxy,iimx,iipx,iimy,iipy,iimz,iipz,izl,izr,parity
      integer(ISZ):: ii,ix,iy,iz,ic
      integer(ISZ):: s_parity,e_parity

c     --- Macro for common arguments to the psor_loop calls
#define LOOP_ARGS iimz,iipz,phi1d,rho1d,nxy,rdx2c,rdy2c,rdz2c,spm1,linbends, \
         rstar,dx,rdx2,rdy2,rdz2,rdel,const,scrtch

c     --- Set temporary variables (these are used to increase performance)
      dx = 1./sqrt(rdx2)
      rdel = rdz2/(rdx2 + rdy2 + rdz2)
      reps0 = 1./eps0
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
      if (bound0 .eq. 0) izl = izl + 1
      if (boundnz.eq. 0) izr = izr - 1

c     --- Save values on the transverse boundaries.
c     --- Both even and odd points are saved for all transverse boundaries.
      do iz=izl-1,izr-1
        do ix=0,nx
          boundarr(ix,2,iz) = phi(ix,0,iz)
          boundarr(ix,4,iz) = phi(ix,ny,iz)
        enddo
        do iy=0,ny
          boundarr(iy,1,iz) = phi(0,iy,iz)
          boundarr(iy,3,iz) = phi(nx,iy,iz)
        enddo
      enddo

c     --- Save values just outside conductor surfaces. Only save phi at the
c     --- subgrid points which are to be used at the current level of
c     --- grid refinement.
      if (lcndbndy) then
        do ic=1,necndbdy
          if (iecndlxy(ic) >= levelxy .and. iecndlz(ic) >= levelz) then
            ecndpvph(ic) = phi(iecndx(ic)/levelxy,
     &                         iecndy(ic)/levelxy,
     &                         iecndz(ic)/levelz)
          endif
        enddo
        if (levelxy .eq. 1 .or. levelz .eq. 1) then
          do ic=1,nocndbdy
            if (iocndlxy(ic) >= levelxy .and. iocndlz(ic) >= levelz) then
              ocndpvph(ic) = phi(iocndx(ic)/levelxy,
     &                           iocndy(ic)/levelxy,
     &                           iocndz(ic)/levelz)
            endif
          enddo
        endif
      endif

c     --- Put desired potential onto conductors in phi array.
      call cond_potmg(nx,ny,nz,phi(0,0,0),levelxy,levelz,.false.)

c     --- Set starting and ending parity.
#ifdef PARALLEL
      s_parity = mod(izslave(my_index),2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 1
      e_parity = 0
#endif

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

#ifdef PARALLEL
c*********************Parallel Routine*********************
c       --- exchange phi on the boundary
        call exchange_phi(nx,ny,nz,phi(0,0,-1),bound0,boundnz,0)
c**********************************************************
#endif

c       --- Set guard planes in z appropriately for the boundary conditions
c       --- For periodic boundaries, this need not be done for the
c       --- parallel version since it automatically happens in exchange_phi

        if (bound0 .eq. 1) then
          do ii = (izl-1)*nxy+1+mod(parity+izl-1,2),(izl-1)*nxy+nxy
            phi1d(ii) = phi1d(ii + 2*nxy)
          enddo
#ifndef PARALLEL
        elseif (bound0 .eq. 2) then
          do ii = (izl-1)*nxy+1+mod(parity+izl-1,2),(izl-1)*nxy+nxy
            phi1d(ii) = phi1d(ii + nz*nxy)
          enddo
#endif
        endif
        if (boundnz .eq. 1) then
          do ii = (izr+1)*nxy+1+mod(parity+izr+1,2),(izr+1)*nxy+nxy
            phi1d(ii) = phi1d(ii - 2*nxy)
          enddo
#ifndef PARALLEL
        elseif (boundnz .eq. 2) then
          do ii = (izr+1)*nxy+1+mod(parity+izr+1,2),(izr+1)*nxy+nxy
            phi1d(ii) = phi1d(ii - nz*nxy)
          enddo
#endif
        endif

c       --- Loop over the rest of the array. Boundary points are calculated
c       --- too, even though the equation is wrong.  They are recalculated
c       --- later.  The case with and without bends are seperated since
c       --- the loops are done differently.  Without bends, the whole array
c       --- can be processed in one loop as a 1d array.  With bends, an
c       --- outer loop over z is needed.
        if (.not. linbends) then 
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
            phi(ix,0,iz) = boundarr(ix,2,iz)
            phi(ix,ny,iz) = boundarr(ix,4,iz)
          enddo
          do iy=mod(iz+parity,2),ny,2
            phi(0,iy,iz) = boundarr(iy,1,iz)
            phi(nx,iy,iz) = boundarr(iy,3,iz)
          enddo
        enddo

        if (boundxy .eq. 1 .or. l2symtry .or. l4symtry) then
c         --- surfaces at ix=0, ix=nx, iy=0, and iy=ny
c         --- if only 2-fold, apply to surfaces at iy=0
c         --- if only 4-fold, apply to surfaces at ix=0 and iy=0
c         --- if also boundxy=1, then apply to all transverse surfaces
          call psor_loop(parity,1,izl,izr,1,
     &                   2,nx,2,iimx,iipx,-iimy,iipy,
     &LOOP_ARGS)
          if (boundxy .eq. 1 .or. l4symtry) then
            call psor_loop(parity,iipy,izl,izr,1,nx+2,(nx+1)*ny-nx,
     &                     2*(nx+1),-iimx,iipx,iimy,iipy,
     &LOOP_ARGS)
c           --- line at transverse edge (ix=0, iy=0)
            call psor_loop(parity,nxy,izl,izl,1,1,nxy*(izr-izl)+1,
     &                     2*nxy,-iimx,iipx,-iimy,iipy,
     &LOOP_ARGS)
          endif
c         --- now do planes at ix=nx and iy=ny
          if (boundxy .eq. 1) then
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
        else if (boundxy .eq. 2) then
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
          if (parity == 0 .or. (levelxy > 1 .or. levelz > 1)) then
           call condbndymg(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,mgparam,
     &                    ncndmax,necndbdy,ecdelmx,ecdelmy,ecdelmz,
     &                    ecdelpx,ecdelpy,ecdelpz,iecndx,iecndy,iecndz,
     &                    iecndlxy,iecndlz,
     &                    ecvoltmx,ecvoltpx,ecvoltmy,ecvoltpy,ecvoltmz,ecvoltpz,
     &                    ecndpvph,bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                    levelxy,levelz,parity)
          endif
          if (parity == 1 .and. (levelxy == 1 .or. levelz == 1)) then
           call condbndymg(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,mgparam,
     &                    ncndmax,nocndbdy,ocdelmx,ocdelmy,ocdelmz,
     &                    ocdelpx,ocdelpy,ocdelpz,iocndx,iocndy,iocndz,
     &                    iocndlxy,iocndlz,
     &                    ocvoltmx,ocvoltpx,ocvoltmy,ocvoltpy,ocvoltmz,ocvoltpz,
     &                    ocndpvph,bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                    levelxy,levelz,parity)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmg(nx,ny,nz,phi(0,0,0),levelxy,levelz,.false.)

c     --- end of loop over even and odd points
      enddo

      return
      end
c=============================================================================
      subroutine cond_potmg(nx,ny,nz,phi,levelxy,levelz,lresidual)
      use PSOR3d
      use MultigridConductor3d
      integer(ISZ):: nx,ny,nz,levelxy,levelz
      real(kind=8):: phi(0:nx,0:ny,0:nz)
      logical(ISZ):: lresidual

c Set conductor points to the desired potential.


      integer(ISZ):: ic

c     --- When at the finest level and not calculating the residual, set
c     --- phi to the voltage of the conductor, otherwise, set it to zero.
      if (levelxy == 1 .and. levelz == 1 .and. .not. lresidual) then
        do ic = 1,ncond
          phi(ixcond(ic),iycond(ic),izcond(ic)) = condvolt(ic)
        enddo
      else
        do ic=1,ncond
c         --- Only set points which are on the coarse grid at the current
c         --- level of refinement.
          if (icondlxy(ic) >= levelxy .and. icondlz(ic) >= levelz) then
            phi(ixcond(ic)/levelxy,iycond(ic)/levelxy,izcond(ic)/levelz) = 0.
          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine condbndymg(nx,ny,nz,phi,rho,rdx2c,rdy2c,rdz2c,spm1,srp,
     &                      ncndmax,ncndbdy,cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,
     &                      cdelpz,icndx,icndy,icndz,icndlxy,icndlz,
     &                      cvoltmx,cvoltpx,cvoltmy,cvoltpy,cvoltmz,cvoltpz,
     &                      cndpvph,bound0,boundnz,boundxy,
     &                      l2symtry,l4symtry,levelxy,levelz,parity)
      integer(ISZ):: nx,ny,nz,ncndmax,ncndbdy,levelxy,levelz,parity
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), rho(0:nx,0:ny,0:nz)
      real(kind=8):: rdx2c,rdy2c,rdz2c,spm1,srp
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax), cdelmz(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax), cdelpz(ncndmax)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax), icndz(ncndmax)
      integer(ISZ):: icndlxy(ncndmax), icndlz(ncndmax)
      real(kind=8):: cvoltmx(ncndmax),cvoltmy(ncndmax),cvoltmz(ncndmax)
      real(kind=8):: cvoltpx(ncndmax),cvoltpy(ncndmax),cvoltpz(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
      integer(ISZ):: bound0,boundnz,boundxy
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
      integer(ISZ):: icdx,icdy,icdz
      logical(ISZ):: dosubgrid

      rdx2cos = rdx2c/srp
      rdy2cos = rdy2c/srp
      rdz2cos = rdz2c/srp

c     --- loop over points near surface of conductors
      do ic = 1,ncndbdy

c       --- Only use point if it is usable at this level or higher
c       --- and if the parity is correct.
        if (icndlxy(ic) >= levelxy .and. icndlz(ic) >= levelz) then

          icdx = icndx(ic)/levelxy
          icdy = icndy(ic)/levelxy
          icdz = icndz(ic)/levelz

          if (mod(icdx+icdy+icdz,2)==parity) then

c           --- set temporaries for boundaries
            icndxp1 = icdx + 1
            icndxm1 = icdx - 1
            icndyp1 = icdy + 1
            icndym1 = icdy - 1
            icndzp1 = icdz + 1
            icndzm1 = icdz - 1

            if (icndxm1 == -1 .and. boundxy == 2) icndxm1 = nx-1
            if (icndxm1 == -1 .and.
     &         (boundxy == 1 .or. l2symtry .or. l4symtry)) icndxm1 = 1
            if (icndxp1 == nx+1 .and. boundxy == 1) icndxp1 = nx-1
            if (icndxp1 == nx+1 .and. boundxy == 2) icndxp1 = 1

            if (icndym1 == -1 .and. boundxy == 2) icndym1 = ny-1
            if (icndym1 == -1 .and.
     &         (boundxy == 1 .or. l4symtry)) icndym1 = 1
            if (icndyp1 == ny+1 .and. boundxy == 1) icndyp1 = ny-1
            if (icndyp1 == ny+1 .and. boundxy == 2) icndyp1 = 1

            if (icndzm1 == -1 .and. bound0 == 1) icndzm1 = +1
            if (icndzm1 == -1 .and. bound0 == 2) icndzm1 = nz-1
            if (icndzp1 == nz+1 .and. boundnz == 1) icndzp1 = nz-1
            if (icndzp1 == nz+1 .and. boundnz == 2) icndzp1 = +1

c           --- Set temporaries with initial values.
            pxm = phi(icndxm1,icdy   ,icdz   )
            pxp = phi(icndxp1,icdy   ,icdz   )
            pym = phi(icdx   ,icndym1,icdz   )
            pyp = phi(icdx   ,icndyp1,icdz   )
            pzm = phi(icdx   ,icdy   ,icndzm1)
            pzp = phi(icdx   ,icdy   ,icndzp1)
            denom = 1.
            dosubgrid = .false.

c           --- Only use actual voltage on finest level. Set to zero for
c           --- coarser levels since solver for the residuals.
            if (levelxy == 1 .and. levelz == 1) then
              voltfac = 1.
            else
              voltfac = 0.
            endif

c           --- the point lower in x is inside the conductor
            if (cdelmx(ic) .lt. levelxy) then
              pxm = voltfac*cvoltmx(ic)/cdelmx(ic)
              denom = denom + (levelxy-cdelmx(ic))/cdelmx(ic)*rdx2cos
              dosubgrid = .true.
            endif
c           --- the point higher in x is inside the conductor
            if (cdelpx(ic) .lt. levelxy) then
              pxp = voltfac*cvoltpx(ic)/cdelpx(ic)
              denom = denom + (levelxy-cdelpx(ic))/cdelpx(ic)*rdx2cos
              dosubgrid = .true.
            endif
c           --- the point lower in y is inside the conductor
            if (cdelmy(ic) .lt. levelxy) then
              pym = voltfac*cvoltmy(ic)/cdelmy(ic)
              denom = denom + (levelxy-cdelmy(ic))/cdelmy(ic)*rdy2cos
              dosubgrid = .true.
            endif
c           --- the point higher in y is inside the conductor
            if (cdelpy(ic) .lt. levelxy) then
              pyp = voltfac*cvoltpy(ic)/cdelpy(ic)
              denom = denom + (levelxy-cdelpy(ic))/cdelpy(ic)*rdy2cos
              dosubgrid = .true.
            endif
c           --- the point lower in z is inside the conductor
            if (cdelmz(ic) .lt. levelz) then
              pzm = voltfac*cvoltmz(ic)/cdelmz(ic)
              denom = denom + (levelz-cdelmz(ic))/cdelmz(ic)*rdz2cos
              dosubgrid = .true.
            endif
c           --- the point lower in z is inside the conductor
            if (cdelpz(ic) .lt. levelz) then
              pzp = voltfac*cvoltpz(ic)/cdelpz(ic)
              denom = denom + (levelz-cdelpz(ic))/cdelpz(ic)*rdz2cos
              dosubgrid = .true.
            endif
c           --- calculate the new phi based on the boundary conditions
            if (dosubgrid)
     &        phi(icdx,icdy,icdz) = (rho(icdx,icdy,icdz) +
     &          (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c + (pzm+pzp)*rdz2c)/denom +
     &          spm1*cndpvph(ic)
          endif
        endif
      enddo

      return
      end
c=============================================================================
      subroutine condbndyres(nx,ny,nz,phi,rho,res,rdx2,rdy2,rdz2,mgparam,
     &                       ncndmax,ncndbdy,cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,
     &                       cdelpz,icndx,icndy,icndz,icndlxy,icndlz,
     &                       cvoltmx,cvoltpx,cvoltmy,cvoltpy,cvoltmz,cvoltpz,
     &                       cndpvph,bound0,boundnz,boundxy,
     &                       l2symtry,l4symtry,levelxy,levelz,parity)
      integer(ISZ):: nx,ny,nz,ncndmax,ncndbdy,levelxy,levelz,parity
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1), rho(0:nx,0:ny,0:nz), res(0:nx,0:ny,0:nz)
      real(kind=8):: rdx2,rdy2,rdz2,mgparam
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax), cdelmz(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax), cdelpz(ncndmax)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax), icndz(ncndmax)
      integer(ISZ):: icndlxy(ncndmax), icndlz(ncndmax)
      real(kind=8):: cvoltmx(ncndmax),cvoltmy(ncndmax),cvoltmz(ncndmax)
      real(kind=8):: cvoltpx(ncndmax),cvoltpy(ncndmax),cvoltpz(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
      integer(ISZ):: bound0,boundnz,boundxy
      logical(ISZ):: l2symtry,l4symtry

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c The result is scaled by the minimum of the deltas. This is done since the
c the correct term can get erroneously large as delta approaches zero.
c
c Temporary variables pxm, pym, pzm, pxp, pyp, and pzp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold icndx-1, icndx+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: const,rdx2c,rdy2c,rdz2c,pxm,pym,pzm,pxp,pyp,pzp,denom,voltfac
      real(kind=8):: rdx2cs,rdy2cs,rdz2cs,ppp
      integer(ISZ):: ic,icndxp1,icndxm1,icndyp1,icndym1,icndzp1,icndzm1
      integer(ISZ):: icdx,icdy,icdz
      real(kind=8):: levelxyi,levelzi

      const = 0.5/(rdx2+rdy2+rdz2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdz2c = rdz2*const
      rdx2cs = mgparam*rdx2*const
      rdy2cs = mgparam*rdy2*const
      rdz2cs = mgparam*rdz2*const
      levelxyi = 1./levelxy
      levelzi = 1./levelz

c     --- loop over points near surface of conductors
      do ic = 1,ncndbdy

c       --- Only use point if it is usable at this level or higher
c       --- and if the parity is correct.
        if (icndlxy(ic) >= levelxy .and. icndlz(ic) >= levelz) then

          icdx = icndx(ic)/levelxy
          icdy = icndy(ic)/levelxy
          icdz = icndz(ic)/levelz

          if (mod(icdx+icdy+icdz,2)==parity) then

c           --- set temporaries for boundaries
            icndxp1 = icdx + 1
            icndxm1 = icdx - 1
            icndyp1 = icdy + 1
            icndym1 = icdy - 1
            icndzp1 = icdz + 1
            icndzm1 = icdz - 1

            if (icndxm1 == -1 .and. boundxy == 2) icndxm1 = nx-1
            if (icndxm1 == -1 .and.
     &         (boundxy == 1 .or. l2symtry .or. l4symtry)) icndxm1 = 1
            if (icndxp1 == nx+1 .and. boundxy == 1) icndxp1 = nx-1
            if (icndxp1 == nx+1 .and. boundxy == 2) icndxp1 = 1

            if (icndym1 == -1 .and. boundxy == 2) icndym1 = ny-1
            if (icndym1 == -1 .and.
     &         (boundxy == 1 .or. l4symtry)) icndym1 = 1
            if (icndyp1 == ny+1 .and. boundxy == 1) icndyp1 = ny-1
            if (icndyp1 == ny+1 .and. boundxy == 2) icndyp1 = 1

            if (icndzm1 == -1 .and. bound0 == 1) icndzm1 = +1
            if (icndzm1 == -1 .and. bound0 == 2) icndzm1 = nz-1
            if (icndzp1 == nz+1 .and. boundnz == 1) icndzp1 = nz-1
            if (icndzp1 == nz+1 .and. boundnz == 2) icndzp1 = +1

c           --- set temporaries with initial values
            pxm = phi(icndxm1,icdy   ,icdz   )
            pxp = phi(icndxp1,icdy   ,icdz   )
            pym = phi(icdx   ,icndym1,icdz   )
            pyp = phi(icdx   ,icndyp1,icdz   )
            pzm = phi(icdx   ,icdy   ,icndzm1)
            pzp = phi(icdx   ,icdy   ,icndzp1)
            denom = 1.
            ppp = 1.

c           --- Only use actual voltage on finest level. Set to zero for
c           --- coarser levels since solver for the residuals.
            if (levelxy == 1 .and. levelz == 1) then
              voltfac = 1.
            else
              voltfac = 0.
            endif

c           --- the point lower in x is inside the conductor
            if (cdelmx(ic) .lt. levelxy) then
              pxm = voltfac*cvoltmx(ic)/cdelmx(ic)
              denom = denom + (levelxy-cdelmx(ic))/cdelmx(ic)*rdx2c
              ppp = min(ppp,cdelmx(ic)*levelxyi)
            endif
c           --- the point higher in x is inside the conductor
            if (cdelpx(ic) .lt. levelxy) then
              pxp = voltfac*cvoltpx(ic)/cdelpx(ic)
              denom = denom + (levelxy-cdelpx(ic))/cdelpx(ic)*rdx2c
              ppp = min(ppp,cdelpx(ic)*levelxyi)
            endif
c           --- the point lower in y is inside the conductor
            if (cdelmy(ic) .lt. levelxy) then
              pym = voltfac*cvoltmy(ic)/cdelmy(ic)
              denom = denom + (levelxy-cdelmy(ic))/cdelmy(ic)*rdy2c
              ppp = min(ppp,cdelmy(ic)*levelxyi)
            endif
c           --- the point higher in y is inside the conductor
            if (cdelpy(ic) .lt. levelxy) then
              pyp = voltfac*cvoltpy(ic)/cdelpy(ic)
              denom = denom + (levelxy-cdelpy(ic))/cdelpy(ic)*rdy2c
              ppp = min(ppp,cdelpy(ic)*levelxyi)
            endif
c           --- the point lower in z is inside the conductor
            if (cdelmz(ic) .lt. levelz) then
              pzm = voltfac*cvoltmz(ic)/cdelmz(ic)
              denom = denom + (levelz-cdelmz(ic))/cdelmz(ic)*rdz2c
              ppp = min(ppp,cdelmz(ic)*levelzi)
            endif
c           --- the point lower in z is inside the conductor
            if (cdelpz(ic) .lt. levelz) then
              pzp = voltfac*cvoltpz(ic)/cdelpz(ic)
              denom = denom + (levelz-cdelpz(ic))/cdelpz(ic)*rdz2c
              ppp = min(ppp,cdelpz(ic)*levelzi)
            endif
c           --- calculate the residual based on the boundary conditions
            if (ppp < 1.)
     &        res(icdx,icdy,icdz) = ppp*(rho(icdx,icdy,icdz)
     &               + (pxm+pxp)*rdx2cs + (pym+pyp)*rdy2cs + (pzm+pzp)*rdz2cs
     &               - phi(icdx,icdy,icdz)*mgparam*denom)
          endif
        endif
      enddo

      return
      end
c=============================================================================
      subroutine residual(nx,ny,nz,dxsqi,dysqi,dzsqi,phi,rho,res,levelxy,levelz,
     &                    l2symtry,l4symtry)
      use PSOR3d
      use MultigridConductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dxsqi,dysqi,dzsqi
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1),rho(0:nx,0:ny,0:nz)
      real(kind=8):: res(0:nx,0:ny,0:nz)
      integer(ISZ):: levelxy,levelz
      logical(ISZ):: l2symtry,l4symtry

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
      izmin = 1
      izmax = nz-1
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny
      if (bound0 > 0) izmin = 0
      if (boundnz == 1) izmax = nz

c     --- Calculate the residual.
      do iz=izmin,izmax
        izm1 = iz - 1
        if (iz == 0 .and. bound0 == 1)  izm1 = 1
        if (iz == 0 .and. bound0 == 2)  izm1 = nz-1
        izp1 = iz + 1
        if (iz == nz .and. boundnz == 1)  izp1 = nz-1
        if (iz == nz .and. boundnz == 2)  izp1 = 1

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

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        do iz=izmin,izmax
          do ix=ixmin,ixmax
            res(ix,ny,iz) = res(ix,0,iz)
          enddo
          do iy=iymin,iymax
            res(nx,iy,iz) = res(0,iy,iz)
          enddo
        enddo
      endif
      if (boundnz == 2) then
        do iy=iymin,iymax
          do ix=ixmin,ixmax
            res(ix,iy,nz) = res(ix,iy,0)
          enddo
        enddo
      endif

c     --- Zero the residual inside conductors.
      call cond_potmg(nx,ny,nz,res,levelxy,levelz,.true.)

c     --- Calculate the residual near the conductor. Note that the odd points
c     --- are only used at the finest level. Also, for higher levels, the even
c     --- points are passed down with an odd parity since they may be odd at
c     --- that level.
      if (levelxy == 1 .or. levelz == 1) then
        call condbndyres(nx,ny,nz,phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,
     &                ncndmax,nocndbdy,ocdelmx,ocdelmy,ocdelmz,ocdelpx,ocdelpy,
     &                ocdelpz,iocndx,iocndy,iocndz,iocndlxy,iocndlz,
     &                ocvoltmx,ocvoltpx,ocvoltmy,ocvoltpy,ocvoltmz,ocvoltpz,
     &                ocndpvph,bound0,boundnz,boundxy,
     &                l2symtry,l4symtry,levelxy,levelz,1)
      endif
      if (levelxy > 1 .or. levelz > 1) then
        call condbndyres(nx,ny,nz,phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,
     &                ncndmax,necndbdy,ecdelmx,ecdelmy,ecdelmz,ecdelpx,ecdelpy,
     &                ecdelpz,iecndx,iecndy,iecndz,iecndlxy,iecndlz,
     &                ecvoltmx,ecvoltpx,ecvoltmy,ecvoltpy,ecvoltmz,ecvoltpz,
     &                ecndpvph,bound0,boundnz,boundxy,
     &                l2symtry,l4symtry,levelxy,levelz,1)
      endif
      call condbndyres(nx,ny,nz,phi,rho,res,dxsqi,dysqi,dzsqi,mgparam,
     &                ncndmax,necndbdy,ecdelmx,ecdelmy,ecdelmz,ecdelpx,ecdelpy,
     &                ecdelpz,iecndx,iecndy,iecndz,iecndlxy,iecndlz,
     &                ecvoltmx,ecvoltpx,ecvoltmy,ecvoltpy,ecvoltmz,ecvoltpz,
     &                ecndpvph,bound0,boundnz,boundxy,
     &                l2symtry,l4symtry,levelxy,levelz,0)

      return
      end
c=============================================================================
      subroutine conductor_data_level(nx,ny,nz,dx,dy,dz)
      use PSOR3d
      use MultigridConductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz

c Set the coarsest level at which each conductor point is still on the grid.
c Also, remove subgrid data points which are not used - those points which are
c not near a conducting surface at any level of refinement.


      integer(ISZ):: ic,ne,no,levelxy,levelz
      integer(ISZ):: tnx,tny,tnz
      real(kind=8):: tdx,tdy,tdz

c     --- Zero arrays, assuming that no points are used.
c     call zeroarry(icondlxy,ncond)
c     call zeroarry(icondlz ,ncond)
c     call zeroarry(iecndlxy,necndbdy)
c     call zeroarry(iecndlz ,necndbdy)
c     call zeroarry(iocndlxy,nocndbdy)
c     call zeroarry(iocndlz ,nocndbdy)
      icondlxy(:ncond) = 0.
      icondlz(:ncond) = 0.
      iecndlxy(:necndbdy) = 0.
      iecndlz(:necndbdy) = 0.
      iocndlxy(:nocndbdy) = 0.
      iocndlz(:nocndbdy) = 0.

c     --- Initialize the temps
      tnx = nx
      tny = ny
      tnz = nz
      tdx = dx
      tdy = dy
      tdz = dz
      levelxy = 1
      levelz = 1

c     --- Loop until at least one of the axis is at the coarsest level.
      do while (tnx >= 2 .and. tny >= 2 .and. tnz >= 2 .and.
     &          mod(tnx,4) == 0 .and. mod(tny,4) == 0 .and. mod(tnz,4) == 0)

c       --- Interior points
        do ic=1,ncond
c         --- Check if point is on the current grid.
          if (mod(ixcond(ic),levelxy) == 0 .and.
     &        mod(iycond(ic),levelxy) == 0 .and.
     &        mod(izcond(ic),levelz ) == 0) then
            icondlxy(ic) = levelxy
            icondlz(ic)  = levelz
          endif
        enddo

c       --- Even points
        do ic=1,necndbdy
c         --- Check if point is near a conductor.
c         --- i.e. is within one grid cell of a conductor on the current grid.
          if (ecdelmx(ic)/levelxy < 1. .or. ecdelpx(ic)/levelxy < 1. .or.
     &        ecdelmy(ic)/levelxy < 1. .or. ecdelpy(ic)/levelxy < 1. .or.
     &        ecdelmz(ic)/levelz  < 1. .or. ecdelpz(ic)/levelz  < 1.) then
c           --- Now, check if point is on the current grid.
            if (mod(iecndx(ic),levelxy) == 0 .and.
     &          mod(iecndy(ic),levelxy) == 0 .and.
     &          mod(iecndz(ic),levelz ) == 0) then
              iecndlxy(ic) = levelxy
              iecndlz(ic)  = levelz
            endif
          endif
        enddo

c       --- Odd points
        do ic=1,nocndbdy
c         --- Check if point is near a conductor.
c         --- i.e. is within one grid cell of a conductor on the current grid.
          if (ocdelmx(ic)/levelxy < 1. .or. ocdelpx(ic)/levelxy < 1. .or.
     &        ocdelmy(ic)/levelxy < 1. .or. ocdelpy(ic)/levelxy < 1. .or.
     &        ocdelmz(ic)/levelz  < 1. .or. ocdelpz(ic)/levelz  < 1.) then
c           --- Now, check if point is on the current grid.
            if (mod(iocndx(ic),levelxy) == 0 .and.
     &          mod(iocndy(ic),levelxy) == 0 .and.
     &          mod(iocndz(ic),levelz ) == 0) then
              iocndlxy(ic) = levelxy
              iocndlz(ic)  = levelz
            endif
          endif
        enddo

c       --- Go to the next level of coarseness.
        tnx = tnx/2
        tny = tny/2
        tdx = tdx*2.
        tdy = tdy*2.
        levelxy = 2*levelxy
        if (1.5*tdz < tdx) then
          tnz = tnz/2
          tdz = tdz*2.
          levelz = 2*levelz
        endif

      end do

c     --- Now, clear out all of the points which were not flagged as usable.
      ne = 0
      do ic=1,necndbdy
        if (iecndlxy(ic) > 0) then
          ne = ne + 1
          if (ne < ic) then
            iecndx(ne)   = iecndx(ic)
            iecndy(ne)   = iecndy(ic)
            iecndz(ne)   = iecndz(ic)
            ecdelmx(ne)  = ecdelmx(ic)
            ecdelmy(ne)  = ecdelmy(ic)
            ecdelmz(ne)  = ecdelmz(ic)
            ecdelpx(ne)  = ecdelpx(ic)
            ecdelpy(ne)  = ecdelpy(ic)
            ecdelpz(ne)  = ecdelpz(ic)
            ecvolt(ne)   = ecvolt(ic)
            ecvoltmx(ne) = ecvoltmx(ic)
            ecvoltpx(ne) = ecvoltpx(ic)
            ecvoltmy(ne) = ecvoltmy(ic)
            ecvoltpy(ne) = ecvoltpy(ic)
            ecvoltmz(ne) = ecvoltmz(ic)
            ecvoltpz(ne) = ecvoltpz(ic)
            iecndlxy(ne) = iecndlxy(ic)
            iecndlz(ne)  = iecndlz(ic)
          endif
        endif
      enddo
      necndbdy = ne

      no = 0
      do ic=1,nocndbdy
        if (iocndlxy(ic) > 0) then
          no = no + 1
          if (no < ic) then
            iocndx(no)   = iocndx(ic)
            iocndy(no)   = iocndy(ic)
            iocndz(no)   = iocndz(ic)
            ocdelmx(no)  = ocdelmx(ic)
            ocdelmy(no)  = ocdelmy(ic)
            ocdelmz(no)  = ocdelmz(ic)
            ocdelpx(no)  = ocdelpx(ic)
            ocdelpy(no)  = ocdelpy(ic)
            ocdelpz(no)  = ocdelpz(ic)
            ocvolt(no)   = ocvolt(ic)
            ocvoltmx(no) = ocvoltmx(ic)
            ocvoltpx(no) = ocvoltpx(ic)
            ocvoltmy(no) = ocvoltmy(ic)
            ocvoltpy(no) = ocvoltpy(ic)
            ocvoltmz(no) = ocvoltmz(ic)
            ocvoltpz(no) = ocvoltpz(ic)
            iocndlxy(no) = iocndlxy(ic)
            iocndlz(no)  = iocndlz(ic)
          endif
        endif
      enddo
      nocndbdy = no

      return
      end
c=============================================================================
      subroutine subgrid_sor_to_mg(nx,ny,nz,dx,dy,dz)
      use PSOR3d
      use MultigridConductor3d
      use Multigrid3d_work
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz

c Convert subgrid data from the format used for the SOR field solver to the
c the format used for multigrid. This includes expanding the data set for
c use on the coarser grid levels. Note that after the conversion, the data
c set will still work for the SOR (though it will be slightly less
c efficient since there will be extra points in the data set that the SOR
c doesn't need).


      integer(ISZ):: ix,iy,iz

c     --- Make sure there is enough room for the extra data.
c     --- The factor of 3 is an over estimate since the number of points
c     --- should typically only double.
      if (ncndmax < 3*necndbdy) then
        ncndmax = 3*necndbdy
        call gchange("PSOR3d",0)
      endif
      call gchange("MultigridConductor3d",0)

c     --- Allocate the work space.
      wnx = nx
      wny = ny
      wnz = nz
      call gallot("Multigrid3d_work",0)
c     --- Explicit loop to zero iii
c     --- Note that zeroarry cannot be used since that routine assumes
c     --- that the input array is of type real(kind=8):: with 64 bit presicion.
c     --- iii is an integer which is not 64 bits long on all machines.
      do iz=0,wnz
        do iy=0,wny
          do ix=0,wnx
            iii(ix,iy,iz) = 0
          enddo
        enddo
      enddo

c     --- Do the conversion.
      call subgrid_sor_to_mg_work(nx,ny,nz,dx,dy,dz)

c     --- Now, get the refinement level for all of the points.
      call conductor_data_level(nx,ny,nz,dx,dy,dz)


      return
      end
c=============================================================================
      subroutine subgrid_sor_to_mg_work(nx,ny,nz,dx,dy,dz)
      use PSOR3d
      use MultigridConductor3d
      use Multigrid3d
      use Multigrid3d_work
      integer(ISZ):: nx,ny,nz
      real(kind=8):: dx,dy,dz

c This does the actual work of converting the subgrid data from the SOR format
c to the multigrid format. To do the conversion if does the following...
c  - register the conducting points on a 3D grid (this just makes it easy
c    to check for the presence of the conductor.
c  - clean up data set, removing any subgrid points which may lie inside
c    of a conductor (those points are harmless to SOR but are damaging
c    to multigrid)
c  - for existing data set, calculate fully the distances to conductors which
c    are greater than the grid cell size.
c  - register the subgrid points on the 3D grid
c  - finally, scan through the grid at each coarse level of resolution, and
c    save points which are within on grid cell of a conductor
c
c The notation for the 3D work grid is...
c   - all of the points inside of conductors are given the same large integer
c     value =((wnx+1)*(wny+1)*(wnz+1) + 1)
c   - for all subgrid points, the index of that point is stored with a sign
c     attached - positive for even points, negative for odd points
c   - the value chosen for inside of conductors is large enough so that
c     is will not be the same as a subgrid point


      integer(ISZ):: icond,ix,iy,iz
      integer(ISZ):: cds(3),maxcdel
      integer(ISZ):: ic,ne,no,levelxy,levelz,tnx,tny,tnz
      real(kind=8):: tdx,tdy,tdz

c     --- If there are not conductors, then immediately return.
      if (ncond == 0 .and. necndbdy == 0 .and. nocndbdy == 0) return

c     --- Set the known conductor points.
      icond = (wnx+1)*(wny+1)*(wnz+1) + 1
      do ic=1,ncond
        iii(ixcond(ic),iycond(ic),izcond(ic)) = icond + ic
      enddo

c     --- Scan through subgrid points and remove all points which lie inside
c     --- of a conductor.
      ne = 0
      do ic=1,necndbdy
        if (iii(iecndx(ic),iecndy(ic),iecndz(ic)) == 0) then
          ne = ne + 1
          if (ne < ic) then
            iecndx(ne)   = iecndx(ic)
            iecndy(ne)   = iecndy(ic)
            iecndz(ne)   = iecndz(ic)
            ecdelmx(ne)  = ecdelmx(ic)
            ecdelmy(ne)  = ecdelmy(ic)
            ecdelmz(ne)  = ecdelmz(ic)
            ecdelpx(ne)  = ecdelpx(ic)
            ecdelpy(ne)  = ecdelpy(ic)
            ecdelpz(ne)  = ecdelpz(ic)
            ecvolt(ne)   = ecvolt(ic)
          endif
        endif
      enddo
      necndbdy = ne

      no = 0
      do ic=1,nocndbdy
        if (iii(iocndx(ic),iocndy(ic),iocndz(ic)) == 0) then
          no = no + 1
          if (no < ic) then
            iocndx(no)   = iocndx(ic)
            iocndy(no)   = iocndy(ic)
            iocndz(no)   = iocndz(ic)
            ocdelmx(no)  = ocdelmx(ic)
            ocdelmy(no)  = ocdelmy(ic)
            ocdelmz(no)  = ocdelmz(ic)
            ocdelpx(no)  = ocdelpx(ic)
            ocdelpy(no)  = ocdelpy(ic)
            ocdelpz(no)  = ocdelpz(ic)
            ocvolt(no)   = ocvolt(ic)
          endif
        endif
      enddo
      nocndbdy = no

c     --- Reset distances which are greater than one to LARGEPOS and
c     --- register subgrid points in the work array iii.
      do ic=1,necndbdy
        iii(iecndx(ic),iecndy(ic),iecndz(ic)) = ic
        if (ecdelmx(ic) > 1.) ecdelmx(ic) = LARGEPOS
        if (ecdelpx(ic) > 1.) ecdelpx(ic) = LARGEPOS
        if (ecdelmy(ic) > 1.) ecdelmy(ic) = LARGEPOS
        if (ecdelpy(ic) > 1.) ecdelpy(ic) = LARGEPOS
        if (ecdelmz(ic) > 1.) ecdelmz(ic) = LARGEPOS
        if (ecdelpz(ic) > 1.) ecdelpz(ic) = LARGEPOS
      enddo
      do ic=1,nocndbdy
        iii(iocndx(ic),iocndy(ic),iocndz(ic)) = -ic
        if (ocdelmx(ic) > 1.) ocdelmx(ic) = LARGEPOS
        if (ocdelpx(ic) > 1.) ocdelpx(ic) = LARGEPOS
        if (ocdelmy(ic) > 1.) ocdelmy(ic) = LARGEPOS
        if (ocdelpy(ic) > 1.) ocdelpy(ic) = LARGEPOS
        if (ocdelmz(ic) > 1.) ocdelmz(ic) = LARGEPOS
        if (ocdelpz(ic) > 1.) ocdelpz(ic) = LARGEPOS
      enddo

c     --- Check through all of the SOR points to reset any distances which
c     --- should not be LARGEPOS. This also copies ecvolt to the voltage for
c     --- the appropriate direction.
      maxcdel = max(nx,ny,nz) + 1
      do ic=1,necndbdy
        cds(1) = iecndx(ic)
        cds(2) = iecndy(ic)
        cds(3) = iecndz(ic)
        call subgrid_scan_axis(1,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                         ecdelmx,ocdelmx,ecdelmx(ic),iecndx(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ecvoltmx(ic))
        call subgrid_scan_axis(1,+1,cds,nx+1,maxcdel,nx,ny,nz,iii,
     &                         ecdelpx,ocdelpx,ecdelpx(ic),iecndx(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ecvoltpx(ic))
        call subgrid_scan_axis(2,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                         ecdelmy,ocdelmy,ecdelmy(ic),iecndy(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ecvoltmy(ic))
        call subgrid_scan_axis(2,+1,cds,ny+1,maxcdel,nx,ny,nz,iii,
     &                         ecdelpy,ocdelpy,ecdelpy(ic),iecndy(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ecvoltpy(ic))
        call subgrid_scan_axis(3,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                         ecdelmz,ocdelmz,ecdelmz(ic),iecndz(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ecvoltmz(ic))
        call subgrid_scan_axis(3,+1,cds,nz+1,maxcdel,nx,ny,nz,iii,
     &                         ecdelpz,ocdelpz,ecdelpz(ic),iecndz(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ecvoltpz(ic))
      enddo
      do ic=1,nocndbdy
        cds(1) = iocndx(ic)
        cds(2) = iocndy(ic)
        cds(3) = iocndz(ic)
        call subgrid_scan_axis(1,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                         ecdelmx,ocdelmx,ocdelmx(ic),iocndx(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ocvoltmx(ic))
        call subgrid_scan_axis(1,+1,cds,nx+1,maxcdel,nx,ny,nz,iii,
     &                         ecdelpx,ocdelpx,ocdelpx(ic),iocndx(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ocvoltpx(ic))
        call subgrid_scan_axis(2,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                         ecdelmy,ocdelmy,ocdelmy(ic),iocndy(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ocvoltmy(ic))
        call subgrid_scan_axis(2,+1,cds,ny+1,maxcdel,nx,ny,nz,iii,
     &                         ecdelpy,ocdelpy,ocdelpy(ic),iocndy(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ocvoltpy(ic))
        call subgrid_scan_axis(3,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                         ecdelmz,ocdelmz,ocdelmz(ic),iocndz(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ocvoltmz(ic))
        call subgrid_scan_axis(3,+1,cds,nz+1,maxcdel,nx,ny,nz,iii,
     &                         ecdelpz,ocdelpz,ocdelpz(ic),iocndz(ic),icond,
     &                         condvolt,ecvolt,ocvolt,ocvoltpz(ic))
      enddo


c     --- Now, scan the grid at each level of refinement to gather points
c     --- which are within one grid cell of a conductor are the coarser
c     --- grid levels.

c     --- Initialize the temps
      tnx = nx
      tny = ny
      tnz = nz
      tdx = dx
      tdy = dy
      tdz = dz
      levelxy = 1
      levelz = 1

c     --- Increase size of conductor arrays if needed. Increase in size
c     --- is only an estimate.
      if (necndbdy == ncndmax .or. nocndbdy == ncndmax) then
        ncndmax = 1.5*max(necndbdy,nocndbdy) + 1
        call gchange("PSOR3d",0)
        call gchange("MultigridConductor3d",0)
      endif

c     --- Loop until at least one of the axis is at the coarsest level.
      ic = necndbdy + 1
      do while (tnx > 2 .and. tny > 2 .and. tnz > 2 .and.
     &          mod(tnx,4) == 0 .and. mod(tny,4) == 0 .and. mod(tnz,4) == 0)

c       --- Go to the next level of coarseness.
        tnx = tnx/2
        tny = tny/2
        tdx = tdx*2.
        tdy = tdy*2.
        levelxy = 2*levelxy
        if (1.5*tdz < tdx) then
          tnz = tnz/2
          tdz = tdz*2.
          levelz = 2*levelz
        endif

c       --- Loop over current grid and scan in each direction to find the
c       --- nearest conductor. If the resulting distances put it near a
c       --- conductor on the coarse grid, then add it to the list of even
c       --- points. (Odd points only appear on the finest level grid.)
        do iz=levelz,nz-levelz,levelz
          cds(3) = iz
          do iy=levelxy,ny-levelxy,levelxy
            cds(2) = iy
            do ix=levelxy,nx-levelxy,levelxy
              cds(1) = ix

c             --- Check if there is already a conductor or subgrid point
c             --- at this location.
              if (iii(ix,iy,iz) == 0) then

c               --- Scan along the 6 directions.
                call subgrid_scan_axis(1,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                                 ecdelmx,ocdelmx,ecdelmx(ic),ix,icond,
     &                                 condvolt,ecvoltmx,ocvoltmx,ecvoltmx(ic))
                call subgrid_scan_axis(1,+1,cds,nx+1,maxcdel,nx,ny,nz,iii,
     &                                 ecdelpx,ocdelpx,ecdelpx(ic),ix,icond,
     &                                 condvolt,ecvoltpx,ocvoltpx,ecvoltpx(ic))
                call subgrid_scan_axis(2,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                                 ecdelmy,ocdelmy,ecdelmy(ic),iy,icond,
     &                                 condvolt,ecvoltmy,ocvoltmy,ecvoltmy(ic))
                call subgrid_scan_axis(2,+1,cds,ny+1,maxcdel,nx,ny,nz,iii,
     &                                 ecdelpy,ocdelpy,ecdelpy(ic),iy,icond,
     &                                 condvolt,ecvoltpy,ocvoltpy,ecvoltpy(ic))
                call subgrid_scan_axis(3,-1,cds,  -1,maxcdel,nx,ny,nz,iii,
     &                                 ecdelmz,ocdelmz,ecdelmz(ic),iz,icond,
     &                                 condvolt,ecvoltmz,ocvoltmz,ecvoltmz(ic))
                call subgrid_scan_axis(3,+1,cds,nz+1,maxcdel,nx,ny,nz,iii,
     &                                 ecdelpz,ocdelpz,ecdelpz(ic),iz,icond,
     &                                 condvolt,ecvoltpz,ocvoltpz,ecvoltpz(ic))

c               --- If the point is near a conductor, save it.
                if (ecdelmx(ic)/levelxy < 1. .or.
     &              ecdelpx(ic)/levelxy < 1. .or.
     &              ecdelmy(ic)/levelxy < 1. .or.
     &              ecdelpy(ic)/levelxy < 1. .or.
     &              ecdelmz(ic)/levelz  < 1. .or.
     &              ecdelpz(ic)/levelz  < 1.) then
                  iecndx(ic) = ix
                  iecndy(ic) = iy
                  iecndz(ic) = iz
                  iii(ix,iy,iz) = ic
                  ic = ic + 1

c                 --- Increase size of conductor arrays if needed.
c                 --- Increase in size is only an estimate.
                  if (ic >= ncndmax) then
                    ncndmax = 1.1*max(ic,nocndbdy) + 1
                    call gchange("PSOR3d",0)
                    call gchange("MultigridConductor3d",0)
                  endif

                endif

              endif
            enddo
          enddo
        enddo

      end do

c     --- Reset necndbdy
      necndbdy = ic - 1
 
      return
      end
c=============================================================================
      subroutine subgrid_scan_axis(ia,sgn,cds,limit,maxcdel,nx,ny,nz,iii,
     &                             ecdel,ocdel,cdel,cnd,icond,
     &                             condvolt,ecvolt,ocvolt,volt)
      integer(ISZ):: ia,sgn,cds(3),limit,maxcdel,nx,ny,nz,icond
      integer(ISZ):: iii(0:nx,0:ny,0:nz)
      real(kind=8):: ecdel(*),ocdel(*),condvolt(*),ecvolt(*),ocvolt(*)
      real(kind=8):: cdel,volt
      integer(ISZ):: cnd

c Scan in the specified direction looking for a conductor or grid edge.
c If a conductor is found, save the distance to and the voltage on the
c conductor.

      do while (cds(ia) .ne. limit)

        if (iii(cds(1),cds(2),cds(3)) >= icond) then
c         --- A point inside a conductor was found. That means that the
c         --- conductor is aligned with the grid.
          cdel = sgn*(cds(ia) - cnd)
          volt = condvolt(iii(cds(1),cds(2),cds(3))-icond)
          cds(ia) = limit
        elseif (iii(cds(1),cds(2),cds(3)) > 0) then
c         --- An even subgrid point was found. Check if it is near a
c         --- conductor.
          if (ecdel(iii(cds(1),cds(2),cds(3))) < maxcdel) then
c           --- If it is near a conductor, save the info.
            cdel = ecdel(iii(cds(1),cds(2),cds(3))) + sgn*(cds(ia) - cnd)
            volt = ecvolt(iii(cds(1),cds(2),cds(3)))
            cds(ia) = limit
          else
c           --- If it is not near a conductor, then keep looking.
            cds(ia) = cds(ia) + sgn
            cdel = LARGEPOS
          endif
        elseif (iii(cds(1),cds(2),cds(3)) < 0) then
c         --- An odd subgrid point was found. Check if it is near a
c         --- conductor.
          if (ocdel(-iii(cds(1),cds(2),cds(3))) < maxcdel) then
c           --- If it is near a conductor, save the info.
            cdel = ocdel(-iii(cds(1),cds(2),cds(3))) + sgn*(cds(ia) - cnd)
            volt = ocvolt(-iii(cds(1),cds(2),cds(3)))
            cds(ia) = limit
          else
c           --- If it is not near a conductor, then keep looking.
            cds(ia) = cds(ia) + sgn
            cdel = LARGEPOS
          endif
        else
c         --- Otherwise, keep looking.
          cds(ia) = cds(ia) + sgn
          cdel = LARGEPOS
        endif
      end do

c     --- Reset the index.
      cds(ia) = cnd

      return
      end
c=============================================================================
      subroutine CoveredNodes(dx, mask, xlo, ylo, zlo, xhi, yhi, zhi)
      use PSOR3d
      use MultigridConductor3d
      integer(ISZ):: xlo, ylo, zlo, xhi, yhi, zhi
      integer(ISZ):: mask(xlo:xhi, ylo:yhi, zlo:zhi)
      real(kind=8):: dx
 
c Returns the mask, which is set to 1 for points which are inside a
c conductor. The input argument dx is used to determine the refinement level.
 
      integer(ISZ):: level,ic
 
c     --- Get the level of coarseness above the finest level.
      level = nint(dx/dxpsor)                                                    
c     --- Zero out the mask
      mask = 0

      do ic=1,ncond
c       --- First check if the conductor point is on the grid at this
c       --- level of refinement.
        if (icondlxy(ic) >= level .and. icondlz(ic) >= level) then
c         --- Then check if the point is within the requested region.
          if (xlo <= ixcond(ic)/level .and. ixcond(ic)/level <= xhi .and.
     &        ylo <= iycond(ic)/level .and. iycond(ic)/level <= yhi .and.
     &        zlo <= izcond(ic)/level .and. izcond(ic)/level <= zhi) then
c           --- If so, set mask to 1.
            mask(ixcond(ic)/level,iycond(ic)/level,izcond(ic)/level) = 1
          endif
        endif
      enddo
 
      return
      end
c=============================================================================
      subroutine NodalCoefficients(dx,coeffs,xlo,ylo,zlo,xhi,yhi,zhi,ncomp)
      use PSOR3d
      use MultigridConductor3d
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi,ncomp
      real(kind=8):: dx,coeffs(xlo:xhi,ylo:yhi,zlo:zhi,0:7)

c Fills the array coeffs with the coeeficients of the 7-point finite
c differenced Poisson's equation, taking into account the changed form of the
c equation near a boundary using the subgrid-scale technique.

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: level
      real(kind=8):: rdx2,coeff0

c     --- First put in the default values
      rdx2 = 1./dx**2
      coeff0 = 2.*(rdx2 + rdx2 + rdx2)
      coeffs(:,:,:,0) = coeff0
      coeffs(:,:,:,1:6) = rdx2
      coeffs(:,:,:,7) = 0.
      
c     --- Get level of coarseness
      level = nint(dx/dxpsor)

c     --- Get coefficients for even points.
      call NodalCoefficientswork(coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,
     &                        necndbdy,
     &                        ecdelmx,ecdelmy,ecdelmz,ecdelpx,ecdelpy,ecdelpz,
     &                        iecndx,iecndy,iecndz,iecndlxy,iecndlz,
     &                        ecvoltmx,ecvoltpx,ecvoltmy,
     &                        ecvoltpy,ecvoltmz,ecvoltpz,
     &                        level)

c     --- Get coefficients for odd points. Only needs to be done at the finest
c     --- level.
      if (level == 1) then
        call NodalCoefficientswork(coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,
     &                        nocndbdy,
     &                        ocdelmx,ocdelmy,ocdelmz,ocdelpx,ocdelpy,ocdelpz,
     &                        iocndx,iocndy,iocndz,iocndlxy,iocndlz,
     &                        ocvoltmx,ocvoltpx,ocvoltmy,
     &                        ocvoltpy,ocvoltmz,ocvoltpz,
     &                        level)
      endif

      return
      end
c=============================================================================
      subroutine NodalCoefficientswork(coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,
     &                        ncndbdy,
     &                        cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,
     &                        cdelpz,icndx,icndy,icndz,icndlxy,icndlz,
     &                        cvoltmx,cvoltpx,cvoltmy,cvoltpy,cvoltmz,cvoltpz,
     &                        level)
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi
      real(kind=8):: coeffs(xlo:xhi,ylo:yhi,zlo:zhi,0:7)
      real(kind=8):: rdx2
      integer(ISZ):: ncndbdy,level
      real(kind=8):: cdelmx(ncndbdy), cdelmy(ncndbdy), cdelmz(ncndbdy)
      real(kind=8):: cdelpx(ncndbdy), cdelpy(ncndbdy), cdelpz(ncndbdy)
      integer(ISZ):: icndx(ncndbdy), icndy(ncndbdy), icndz(ncndbdy)
      integer(ISZ):: icndlxy(ncndbdy), icndlz(ncndbdy)
      real(kind=8):: cvoltmx(ncndbdy),cvoltmy(ncndbdy),cvoltmz(ncndbdy)
      real(kind=8):: cvoltpx(ncndbdy),cvoltpy(ncndbdy),cvoltpz(ncndbdy)

c Does the work of setting the coefficients. The even or odd conductor
c points can be passed into this routine.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: coeff0(6)

c     --- loop over points near surface of conductors
      do ic = 1,ncndbdy
 
c       --- Only use point if it is usable at this level or higher
c       --- and if the point is within the requested domain
        if (icndlxy(ic) >= level .and. icndlz(ic) >= level) then
 
          ix = icndx(ic)/level
          iy = icndy(ic)/level
          iz = icndz(ic)/level
 
          if (xlo <= ix .and. ix <= xhi .and.
     &        ylo <= iy .and. iy <= yhi .and.
     &        zlo <= iz .and. iz <= zhi) then
 
c           --- Contributions of terms with no conductor boundary nearby.
            coeff0 = -rdx2

c           --- the point lower in x is inside the conductor
            if (cdelmx(ic) .lt. level) then
              coeff0(2) = 1./cdelmx(ic)
              coeffs(ix,iy,iz,2) = 0.
              coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7)+cvoltmx(ic)/cdelmx(ic)*rdx2
            endif
c           --- the point higher in x is inside the conductor
            if (cdelpx(ic) .lt. level) then
              coeff0(1) = 1./cdelpx(ic)
              coeffs(ix,iy,iz,1) = 0.
              coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7)+cvoltpx(ic)/cdelpx(ic)*rdx2
            endif
c           --- the point lower in y is inside the conductor
            if (cdelmy(ic) .lt. level) then
              coeff0(4) = 1./cdelmy(ic)
              coeffs(ix,iy,iz,4) = 0.
              coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7)+cvoltmy(ic)/cdelmy(ic)*rdx2
            endif
c           --- the point higher in y is inside the conductor
            if (cdelpy(ic) .lt. level) then
              coeff0(3) = 1./cdelpy(ic)
              coeffs(ix,iy,iz,3) = 0.
              coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7)+cvoltpy(ic)/cdelpy(ic)*rdx2
            endif
c           --- the point lower in z is inside the conductor
            if (cdelmz(ic) .lt. level) then
              coeff0(6) = 1./cdelmz(ic)
              coeffs(ix,iy,iz,6) = 0.
              coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7)+cvoltmz(ic)/cdelmz(ic)*rdx2
            endif
c           --- the point lower in z is inside the conductor
            if (cdelpz(ic) .lt. level) then
              coeff0(5) = 1./cdelpz(ic)
              coeffs(ix,iy,iz,5) = 0.
              coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7)+cvoltpz(ic)/cdelpz(ic)*rdx2
            endif
c           --- Sum contributions to center term.
            coeffs(ix,iy,iz,0) = sum(coeff0)

          endif
        endif
      enddo

      return
      end
c=============================================================================
