#include "top.h"
c=============================================================================
c@(#) File FXY_MGRID.M, version $Revision: 1.13 $, $Date: 2009/02/04 00:42:10 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 2D multigrid field sovler which is part of the FXY
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigridxyf(iwhich,nx,ny,dx,dy,phi,rho,
     &                       l2symtry,l4symtry,xmmin,ymmin)
      use Constant
      use GridBoundary3d
      use Multigrid3d
      use Conductor3d
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny
      real(kind=8):: phi(0:nx,0:ny)
      real(kind=8):: rho(0:nx,0:ny)
      real(kind=8):: dx,dy
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: xmmin,ymmin

c Use the multigrid method for solving Poisson's equation on a 2-D Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution.
c
c It is assumed that dx ~ dy and that semi-coarsening is not needed
c transversely.
c
c All dimensions, nx and ny, must be powers of 2.

      integer(ISZ):: nxy
      real(kind=8):: maxerr
      real(kind=8):: dxsqi,dysqi,reps0c,rdel
      real(kind=8):: tdx,tdy
      integer(ISZ):: tnx,tny
      integer(ISZ):: i,ii,k
      real(kind=8):: phisave(0:nx,0:ny)

c     --- Initialize temporaries
      nxy  =(nx+1)*(ny+1)
      dxsqi=1./dx**2
      dysqi=1./dy**2
      reps0c = mgparam/(eps0*2.*(dxsqi+dysqi))

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

      call checkconductorsxy(nx,ny,dx,dy,l2symtry,l4symtry,conductors)

c     --- Preset rho to increase performance (reducing the number of
c     --- multiplies in the main SOR sweep loop).
      rho = rho*reps0c

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- vary few iterations are done.
      mgiters = 0
      maxerr = 1.
      do while (maxerr > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Save current value of phi
        phisave = phi

c       --- Do one vcycle.
        call vcyclexy(0,nx,ny,dxsqi,dysqi,phi,rho,
     &                l2symtry,l4symtry,boundxy,mgmaxlevels,
     &                downpasses,uppasses,
     &                lcndbndy,icndbndy,conductors)

        maxerr = maxval(abs(phisave - phi))

      enddo

c     --- Make a print out.
      if (maxerr > mgtol) then
        print*,"Multigrid: Maximum number of iterations reached"
      endif
      print*,"Multigrid: Error converged to ",maxerr," in ",mgiters," v-cycles"

      reps0c = 1./reps0c
      rho = rho*reps0c

      return
      end
c=============================================================================
      subroutine getmglevelsxy(nx,ny,dx,dy,conductors)
      use Parallel
      use ConductorTypemodule
      integer(ISZ):: nx,ny
      real(kind=8):: dx,dy
      type(ConductorType):: conductors

      call getmglevelsrecurxy(nx,ny,dx,dy,0,1,1,conductors)

      return
      end
c=============================================================================
      RECURSIVE subroutine getmglevelsrecurxy(nx,ny,dx,dy,mglevel,lx,ly,
     &                                        conductors)
      use ConductorTypemodule
      integer(ISZ):: nx,ny
      real(kind=8):: dx,dy
      integer(ISZ):: mglevel,lx,ly
      type(ConductorType):: conductors

c Calculate the levels of coarsening, saving the nx, and ny for each level.

      conductors%levellx(mglevel) = lx
      conductors%levelly(mglevel) = ly
      conductors%levels = mglevel + 1

      if (mod(nx,4) /= 0 .or.
     &    mod(ny,4) /= 0) return

      call getmglevelsrecurxy(nx/2,ny/2,dx*2,dy*2,mglevel+1,lx*2,ly*2,
     &                        conductors)

      return
      end
c=============================================================================
      RECURSIVE subroutine vcyclexy(mglevel,nx,ny,dxsqi,dysqi,
     &                              phi,rho,l2symtry,l4symtry,boundxy,
     &                              mgmaxlevels,downpasses,uppasses,
     &                              lcndbndy,icndbndy,conductors)
      use ConductorTypemodule
      integer(ISZ):: mglevel
      integer(ISZ):: nx,ny
      real(kind=8):: dxsqi,dysqi
      real(kind=8):: phi(0:nx,0:ny),rho(0:nx,0:ny)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: boundxy,mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8),allocatable:: phi2(:,:),rho2(:,:)
      real(kind=8):: res(0:nx,0:ny)
      integer(ISZ):: i2h,i

c     --- Do initial SOR passes.
      do i=1,downpasses
        call sorpassxy(mglevel,nx,ny,phi,rho,phi,rho,
     &                 dxsqi,dysqi,l2symtry,l4symtry,
     &                 lcndbndy,icndbndy,conductors)
      enddo

      if (mglevel < mgmaxlevels-1) then

c       --- Get the residual on the current grid.
        call residualxy(nx,ny,dxsqi,dysqi,phi,rho,res,
     &                  mglevel,l2symtry,l4symtry,
     &                  lcndbndy,icndbndy,conductors)

        allocate(phi2(0:nx/2,0:ny/2),
     &           rho2(0:nx/2,0:ny/2))
        phi2 = 0.
        rho2 = 0.

        call restrictxy(nx,ny,res,rho2,boundxy,l2symtry,l4symtry)

c         --- Continue at the next coarsest level.
          call vcyclexy(mglevel+1,nx/2,ny/2,dxsqi*0.25,dysqi*0.25,phi2,rho2,
     &                  l2symtry,l4symtry,boundxy,mgmaxlevels,
     &                  downpasses,uppasses,
     &                  lcndbndy,icndbndy,conductors)

c         --- Add in resulting error.
          call expandxy(nx/2,ny/2,phi2,phi,boundxy)

          deallocate(phi2,rho2)
      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call sorpassxy(mglevel,nx,ny,phi,rho,phi,rho,
     &                 dxsqi,dysqi,l2symtry,l4symtry,
     &                 lcndbndy,icndbndy,conductors)
      enddo

      return
      end
c=============================================================================
      subroutine restrictxy(nx,ny,ph,p2h,boundxy,l2symtry,l4symtry)
      integer(ISZ):: nx,ny
      real(kind=8):: ph(0:nx,0:ny)
      real(kind=8):: p2h(0:nx/2,0:ny/2)
      integer(ISZ):: boundxy
      logical(ISZ):: l2symtry,l4symtry

      integer(ISZ):: ix,iy
      integer(ISZ):: ixmin,ixmax,iymin,iymax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 2
      ixmax = nx-2
      iymin = 2
      iymax = ny-2
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny

c     --- Do the loops.
      do iy=iymin,iymax,2
        iym1 = iy - 1
        if (iy==0 .and.  boundxy==2) iym1 = ny-1
        if (iy==0 .and. (boundxy==1 .or. l4symtry)) iym1 = 1
        iyp1 = iy + 1
        if (iy == ny .and. boundxy == 1) iyp1 = ny-1
        if (iy == ny .and. boundxy == 2) iyp1 = 1

        do ix=ixmin,ixmax,2
          ixm1 = ix - 1
          if (ix==0 .and.  boundxy==2) ixm1 = nx-1
          if (ix==0 .and. (boundxy==1 .or. l2symtry .or. l4symtry)) ixm1 = 1
          ixp1 = ix + 1
          if (ix == nx .and. boundxy == 1)  ixp1 = nx-1
          if (ix == nx .and. boundxy == 2)  ixp1 = 1

          p2h(ix/2,iy/2) =
     &      4.*0.2500*ph(ix  ,iy  ) +
     &      4.*0.1250*ph(ixm1,iy  ) +
     &      4.*0.1250*ph(ixp1,iy  ) +
     &      4.*0.1250*ph(ix  ,iym1) +
     &      4.*0.1250*ph(ix  ,iyp1) +
     &      4.*0.0625*ph(ixm1,iym1) +
     &      4.*0.0625*ph(ixp1,iym1) +
     &      4.*0.0625*ph(ixm1,iyp1) +
     &      4.*0.0625*ph(ixp1,iyp1)

        enddo
      enddo

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        p2h(ixmin/2:ixmax/2,ny/2) = p2h(ixmin/2:ixmax/2,0)
        p2h(nx/2,iymin/2:iymax/2) = p2h(0,iymin/2:iymax/2)
      endif

      return
      end
c=============================================================================
      subroutine expandxy(nx,ny,p2h,ph,boundxy)
      integer(ISZ):: nx,ny
      real(kind=8):: p2h(0:nx,0:ny)
      real(kind=8):: ph(0:nx*2,0:ny*2)
      integer(ISZ):: boundxy

c Add the error on the coarser grid to the current value on the finer grid.
c The expansion is only transverse.

      integer(ISZ):: ix,iy

c     --- Do the loops.
      do iy=0,ny-1
        do ix=0,nx-1

          ph(2*ix  ,2*iy  )=ph(2*ix  ,2*iy  )+      p2h(ix  ,iy  )
          ph(2*ix+1,2*iy  )=ph(2*ix+1,2*iy  )+0.50*(p2h(ix  ,iy  )+
     &                                              p2h(ix+1,iy  ))
          ph(2*ix  ,2*iy+1)=ph(2*ix  ,2*iy+1)+0.50*(p2h(ix  ,iy  )+
     &                                              p2h(ix  ,iy+1))
          ph(2*ix+1,2*iy+1)=ph(2*ix+1,2*iy+1)+0.25*(p2h(ix  ,iy  )+
     &                                              p2h(ix+1,iy  )+
     &                                              p2h(ix  ,iy+1)+
     &                                              p2h(ix+1,iy+1))
        enddo
      enddo

      if (boundxy > 0) then
c       --- Expand ix=nx, iy=ny point.
        ph(2*nx,2*ny) = ph(2*nx,2*ny)+p2h(nx,ny)

c       --- Expand ix=nx line.
        do iy=0,ny-1
          ph(2*nx,2*iy  ) = ph(2*nx,2*iy  ) +      p2h(nx,iy  )
          ph(2*nx,2*iy+1) = ph(2*nx,2*iy+1) + 0.5*(p2h(nx,iy  ) +
     &                                             p2h(nx,iy+1))
        enddo

c       --- Expand iy=ny line.
        do ix=0,nx-1
          ph(2*ix  ,2*ny) = ph(2*ix  ,2*ny) +      p2h(ix  ,ny)
          ph(2*ix+1,2*ny) = ph(2*ix+1,2*ny) + 0.5*(p2h(ix  ,ny) +
     &                                             p2h(ix+1,ny))
        enddo
      endif

      return
      end
c=============================================================================
c=============================================================================
      subroutine sorpassxy(mglevel,nx,ny,phi,rho,phi1d,rho1d,
     &                     rdx2,rdy2,l2symtry,l4symtry,
     &                     lcndbndy,icndbndy,conductors)
      use GridBoundary3d
      use Multigrid3d
      use Constant
      use ConductorTypemodule
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      integer(ISZ):: mglevel,nx,ny
      real(kind=8):: phi(0:nx,0:ny),rho(0:nx,0:ny)
      real(kind=8):: phi1d(*),rho1d(*)
      real(kind=8):: rdx2,rdy2
      logical(ISZ):: l2symtry,l4symtry
      type(ConductorType):: conductors

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.  It makes calls to the routines which specify internal
c conductors.
c
c The tranverse boundaries can either be held constant, have zero normal
c derivative, or be periodic.  When BOUNDXY is zero, the boundaries are held
c constant, when 1, they have zero normal derivative, and when 2, the
c boundaries are periodic.
c
c Note that loops over all directions assume that nx and ny are even.
c
c The arrangement of the loops was done to increase performance.  The entire
c grid is looped over as if it were a 1D array, ignoring boundaries.
c The boundaries are then reset, the previous value was destroyed.

      real(kind=8):: const,rdx2c,rdy2c,spm1,dx
      integer(ISZ):: nxy,iimx,iipx,iimy,iipy,parity
      integer(ISZ):: ii,ix,iy,ic,i1,i2
      integer(ISZ):: s_parity,e_parity
      real(kind=8):: boundarrx(0:nx,2),boundarry(0:ny,2)

c     --- Set temporary variables (these are used to increase performance)
      dx = 1./sqrt(rdx2)
      const = mgparam*0.5/(rdx2 + rdy2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      spm1 = 1. - mgparam

c     --- Set indices for 1d arrays used in the five point finite difference
c     --- form of Poisson's equation.
      nxy = (nx+1)*(ny+1)
      iimx = -1
      iipx = +1
      iimy = -nx-1
      iipy = +nx+1

c     --- Put desired potential onto conductors in phi array.
      call cond_potmgxy(conductors%interior,nx,ny,phi(0,0),mglevel,.false.)

c     --- Save values on the transverse boundaries.
c     --- Both even and odd points are saved for all transverse boundaries.
      do ix=0,nx
        boundarrx(ix,1) = phi(ix,0)
        boundarrx(ix,2) = phi(ix,ny)
      enddo
      do iy=0,ny
        boundarry(iy,1) = phi(0,iy)
        boundarry(iy,2) = phi(nx,iy)
      enddo

c     --- Save values just outside conductor surfaces. Only save phi at the
c     --- subgrid points which are to be used at the current level of
c     --- grid refinement.
      if (lcndbndy) then
        i1 = conductors%evensubgrid%istart(mglevel)
        i2 = conductors%evensubgrid%istart(mglevel+1)-1
        do ic = i1,i2
          ix = conductors%evensubgrid%indx(0,ic)
          iy = conductors%evensubgrid%indx(1,ic)
          conductors%evensubgrid%prevphi(ic) = phi(ix,iy)
        enddo
        i1 = conductors%oddsubgrid%istart(mglevel)
        i2 = conductors%oddsubgrid%istart(mglevel+1)-1
        do ic = i1,i2
          ix = conductors%oddsubgrid%indx(0,ic)
          iy = conductors%oddsubgrid%indx(1,ic)
          conductors%oddsubgrid%prevphi(ic) = phi(ix,iy)
        enddo
      endif

c     --- Set starting and ending parity.
      s_parity = 1
      e_parity = 0

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

c       --- Loop over the the array. Boundary points are calculated
c       --- too, even though the equation is wrong.  They are recalculated
c       --- later.
        do ii = parity+1,nxy, 2
          phi1d(ii) = rho1d(ii) +
     &                (phi1d(ii+iimx) + phi1d(ii+iipx))*rdx2c +
     &                (phi1d(ii+iimy) + phi1d(ii+iipy))*rdy2c +
     &                spm1*phi1d(ii)
        enddo

c       --- Transverse boundaries     
c       --- Restore only even or odd boundary points to previous value
c       --- since only want to restore the values changed from the 1d
c       --- loop.  This automatically takes care of Dirichlet boundaries.
        do ix=parity,nx,2
          phi(ix,0)  = boundarrx(ix,1)
          phi(ix,ny) = boundarrx(ix,2)
        enddo
        do iy=parity,ny,2
          phi(0,iy)  = boundarry(iy,1)
          phi(nx,iy) = boundarry(iy,2)
        enddo

        if (boundxy == 1 .or. l2symtry .or. l4symtry) then
c         --- lines at ix=0, ix=nx, iy=0, and iy=ny
c         --- if only 2-fold, apply to surfaces at iy=0
c         --- if only 4-fold, apply to surfaces at ix=0 and iy=0
c         --- if also boundxy=1, then apply to all transverse surfaces
          do ix=2-parity,nx-1,2
            phi(ix,0) = rho(ix,0) + (phi(ix-1,0) + phi(ix+1,0))*rdx2c +
     &                              (phi(ix  ,1) + phi(ix  ,1))*rdy2c +
     &                              spm1*phi(ix,0)
          enddo
          if (boundxy == 1 .or. l4symtry) then
            do iy=2-parity,ny-1,2
              phi(0,iy) = rho(0,iy) + (phi(1,iy  ) + phi(1,iy  ))*rdx2c +
     &                                (phi(0,iy-1) + phi(0,iy+1))*rdy2c +
     &                                spm1*phi(0,iy)
            enddo
c           --- point at transverse edge (ix=0, iy=0)
            if (parity == 0) then
              phi(0,0) = rho(0,0) + (phi(1,0) + phi(1,0))*rdx2c +
     &                              (phi(0,1) + phi(0,1))*rdy2c +
     &                              spm1*phi(0,0)
            endif
          endif
c         --- now do lines at ix=nx and iy=ny
          if (boundxy == 1) then
            do ix=2-parity,nx-1,2
              phi(ix,ny) = rho(ix,ny)+(phi(ix-1,ny  ) + phi(ix+1,ny  ))*rdx2c +
     &                                (phi(ix  ,ny-1) + phi(ix  ,ny-1))*rdy2c +
     &                                spm1*phi(ix,ny)
            enddo
            do iy=2-parity,ny-1,2
              phi(ny,iy) = rho(ny,iy)+(phi(nx-1,iy  ) + phi(nx-1,iy  ))*rdx2c +
     &                                (phi(ny  ,iy-1) + phi(ny  ,iy+1))*rdy2c +
     &                                spm1*phi(nx,iy)
            enddo
c           --- points at other transverse edges
            if (parity == 0) then
              phi(0,ny) = rho(0,ny) + (phi(1,ny  ) + phi(1,ny  ))*rdx2c +
     &                                (phi(0,ny-1) + phi(0,ny-1))*rdy2c +
     &                                spm1*phi(0,ny)
              phi(ny,0) = rho(ny,0) + (phi(nx-1,0) + phi(nx-1,0))*rdx2c +
     &                                (phi(ny  ,1) + phi(ny  ,1))*rdy2c +
     &                                spm1*phi(nx,0)
              phi(nx,ny) = rho(nx,ny)+(phi(nx-1,ny  ) + phi(nx-1,ny  ))*rdx2c +
     &                                (phi(nx  ,ny-1) + phi(nx  ,ny-1))*rdy2c +
     &                                spm1*phi(nx,ny)
            endif
          endif
        else if (boundxy == 2) then
c         --- lines surfaces at ix=0, ix=nx, iy=0, and iy=ny
          do ix=2-parity,nx-1,2
            phi(ix,0) = rho(ix,0) + (phi(ix-1,0   ) + phi(ix+1,0))*rdx2c +
     &                              (phi(ix  ,ny-1) + phi(ix  ,1))*rdy2c +
     &                              spm1*phi(ix,0)
          enddo
          do iy=2-parity,ny-1,2
            phi(0,iy) = rho(0,iy) + (phi(nx-1,iy  ) + phi(1,iy  ))*rdx2c +
     &                              (phi(0   ,iy-1) + phi(0,iy+1))*rdy2c +
     &                              spm1*phi(0,iy)
          enddo
          phi(2-parity:nx-1:2,ny) = phi(2-parity:nx-1:2,0)
          phi(nx,2-parity:ny-1:2) = phi(0,2-parity:ny-1:2)
c         --- points at transverse edges
          if (parity == 0) then
            phi(0,0) = rho(0,0) + (phi(nx-1,0   ) + phi(1,0))*rdx2c +
     &                            (phi(0   ,ny-1) + phi(0,1))*rdy2c +
     &                            spm1*phi(0,0)
            phi(nx,0) = phi(0,0)
            phi(0,ny) = phi(0,0)
            phi(nx,ny) = phi(0,0)
          endif
        endif
c       --- end of transverse boundaries

c       --- Apply altered difference equation to the points near the
c       --- surface of the conductor boundaries.
        if (lcndbndy) then
          if (parity == 0) then
            call condbndymgxy(conductors%evensubgrid,nx,ny,phi,rho,rdx2c,rdy2c,
     &                        spm1,mgparam,boundxy,l2symtry,l4symtry,mglevel)
          endif
          if (parity == 1) then
            call condbndymgxy(conductors%oddsubgrid,nx,ny,phi,rho,rdx2c,rdy2c,
     &                        spm1,mgparam,boundxy,l2symtry,l4symtry,mglevel)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgxy(conductors%interior,nx,ny,phi(0,0),mglevel,.false.)

c     --- end of loop over even and odd points
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine cond_potmgxy(interior,nx,ny,phi,mglevel,lresidual)
      use GridBoundary3d
      use ConductorInteriorTypemodule
      type(ConductorInteriorType):: interior
      integer(ISZ):: nx,ny,mglevel
      real(kind=8):: phi(0:nx,0:ny)
      logical(ISZ):: lresidual

c Set conductor points to the desired potential.

      integer(ISZ):: ic,ix,iy,iz

c     --- When at the finest level and not calculating the residual, set
c     --- phi to the voltage of the conductor, otherwise, set it to zero.
      if (mglevel == 0 .and. .not. lresidual) then
        do ic = 1,interior%n
          ix = interior%indx(0,ic)
          iy = interior%indx(1,ic)
          phi(ix,iy) = interior%volt(ic)
        enddo
      else
        do ic = interior%istart(mglevel),interior%istart(mglevel+1)-1
          ix = interior%indx(0,ic)
          iy = interior%indx(1,ic)
          phi(ix,iy) = 0.
        enddo
      endif

      return
      end
c=============================================================================
      subroutine condbndymgxy(subgrid,nx,ny,phi,rho,rdx2c,rdy2c,spm1,srp,
     &                        boundxy,l2symtry,l4symtry,mglevel)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,mglevel,parity
      real(kind=8):: phi(0:nx,0:ny), rho(0:nx,0:ny)
      real(kind=8):: rdx2c,rdy2c,spm1,srp
      integer(ISZ):: boundxy
      logical(ISZ):: l2symtry,l4symtry

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c NOTE that rdx2cos and rdy2cos are rdx2c and rdy2c over mgparam.
c
c Temporary variables pxm, pym, pxp, and pyp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: rdx2cos,rdy2cos,pxm,pym,pxp,pyp,denom
      real(kind=8):: dxm,dxp,dym,dyp
      real(kind=8):: voltfac,rhovolume
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1
      integer(ISZ):: ix,iy
      logical(ISZ):: dosubgrid
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      rdx2cos = rdx2c/srp
      rdy2cos = rdy2c/srp
      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1

        if (ixm1 == -1 .and. boundxy == 2) ixm1 = nx-1
        if (ixm1 == -1 .and.
     &     (boundxy == 1 .or. l2symtry .or. l4symtry)) ixm1 = 1
        if (ixp1 == nx+1 .and. boundxy == 1) ixp1 = nx-1
        if (ixp1 == nx+1 .and. boundxy == 2) ixp1 = 1

        if (iym1 == -1 .and. boundxy == 2) iym1 = ny-1
        if (iym1 == -1 .and.
     &     (boundxy == 1 .or. l4symtry)) iym1 = 1
        if (iyp1 == ny+1 .and. boundxy == 1) iyp1 = ny-1
        if (iyp1 == ny+1 .and. boundxy == 2) iyp1 = 1

c       --- Set temporaries with initial values.
        pxm = phi(ixm1,iy   )
        pxp = phi(ixp1,iy   )
        pym = phi(ix   ,iym1)
        pyp = phi(ix   ,iyp1)
        denom = 1.
        dosubgrid = .false.
        dxm = 1.
        dxp = 1.
        dym = 1.
        dyp = 1.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- the point lower in x is inside the conductor
        if (dels(0,ic) < 1.) then
          pxm = voltfac*volt(0,ic)/dels(0,ic)
          denom = denom + (1.-dels(0,ic))/dels(0,ic)*rdx2cos
          dxm = dels(0,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          pxp = voltfac*volt(1,ic)/dels(1,ic)
          denom = denom + (1.-dels(1,ic))/dels(1,ic)*rdx2cos
          dxp = dels(1,ic)
          dosubgrid = .true.
        endif
c       --- the point lower in y is inside the conductor
        if (dels(2,ic) < 1.) then
          pym = voltfac*volt(2,ic)/dels(2,ic)
          denom = denom + (1.-dels(2,ic))/dels(2,ic)*rdy2cos
          dym = dels(2,ic)
          dosubgrid = .true.
        endif
c       --- the point higher in y is inside the conductor
        if (dels(3,ic) < 1.) then
          pyp = voltfac*volt(3,ic)/dels(3,ic)
          denom = denom + (1.-dels(3,ic))/dels(3,ic)*rdy2cos
          dyp = dels(3,ic)
          dosubgrid = .true.
        endif
c       --- calculate the new phi based on the boundary conditions
        if (dosubgrid) then
c         rhovolume = (dxp+dxm)*(dyp+dym)*0.125
          rhovolume = 1.
          phi(ix,iy) = (rho(ix,iy)*rhovolume +
     &      (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c)/denom +
     &      spm1*subgrid%prevphi(ic)
        endif
      enddo

      return
      end
c=============================================================================
      subroutine condbndyresxy(subgrid,nx,ny,phi,rho,res,rdx2,rdy2,mgparam,
     &                         boundxy,l2symtry,l4symtry,mglevel)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: nx,ny,mglevel
      real(kind=8):: phi(0:nx,0:ny)
      real(kind=8):: rho(0:nx,0:ny), res(0:nx,0:ny)
      real(kind=8):: rdx2,rdy2,mgparam
      integer(ISZ):: boundxy
      logical(ISZ):: l2symtry,l4symtry

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.
c The result is scaled by the minimum of the deltas. This is done since the
c the correct term can get erroneously large as delta approaches zero.
c
c Temporary variables pxm, pym, pxp, and pyp hold phi at minus
c and plus one in each direction from the current point.
c These are changed when the finite difference in the appropriate direction
c includes the boundary condition.
c
c Note that care has been taken with the boundaries.  Temps are set up
c to hold ix-1, ix+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: const,rdx2c,rdy2c,pxm,pym,pxp,pyp,denom
      real(kind=8):: dxm,dxp,dym,dyp,rhovolume
      real(kind=8):: voltfac
      real(kind=8):: rdx2cs,rdy2cs,ppp
      integer(ISZ):: ic,ixp1,ixm1,iyp1,iym1
      integer(ISZ):: ix,iy
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      const = 0.5/(rdx2+rdy2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdx2cs = mgparam*rdx2*const
      rdy2cs = mgparam*rdy2*const
      dels => subgrid%dels
      volt => subgrid%volt

c     --- loop over points near surface of conductors
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)

c       --- set temporaries for boundaries
        ixp1 = ix + 1
        ixm1 = ix - 1
        iyp1 = iy + 1
        iym1 = iy - 1

        if (ixm1 == -1 .and. boundxy == 2) ixm1 = nx-1
        if (ixm1 == -1 .and.
     &     (boundxy == 1 .or. l2symtry .or. l4symtry)) ixm1 = 1
        if (ixp1 == nx+1 .and. boundxy == 1) ixp1 = nx-1
        if (ixp1 == nx+1 .and. boundxy == 2) ixp1 = 1

        if (iym1 == -1 .and. boundxy == 2) iym1 = ny-1
        if (iym1 == -1 .and.
     &     (boundxy == 1 .or. l4symtry)) iym1 = 1
        if (iyp1 == ny+1 .and. boundxy == 1) iyp1 = ny-1
        if (iyp1 == ny+1 .and. boundxy == 2) iyp1 = 1

c       --- set temporaries with initial values
        pxm = phi(ixm1,iy   )
        pxp = phi(ixp1,iy   )
        pym = phi(ix   ,iym1)
        pyp = phi(ix   ,iyp1)
        denom = 1.
        ppp = 1.
        dxm = 1.
        dxp = 1.
        dym = 1.
        dyp = 1.

c       --- Only use actual voltage on finest level. Set to zero for
c       --- coarser levels since solver for the residuals.
        if (mglevel == 0) then
          voltfac = 1.
        else
          voltfac = 0.
        endif

c       --- the point lower in x is inside the conductor
        if (dels(0,ic) < 1.) then
          pxm = voltfac*volt(0,ic)/dels(0,ic)
          denom = denom + (1.-dels(0,ic))/dels(0,ic)*rdx2c
          dxm = dels(0,ic)
          ppp = min(ppp,dels(0,ic))
        endif
c       --- the point higher in x is inside the conductor
        if (dels(1,ic) < 1.) then
          pxp = voltfac*volt(1,ic)/dels(1,ic)
          denom = denom + (1.-dels(1,ic))/dels(1,ic)*rdx2c
          dxp = dels(1,ic)
          ppp = min(ppp,dels(1,ic))
        endif
c       --- the point lower in y is inside the conductor
        if (dels(2,ic) < 1.) then
          pym = voltfac*volt(2,ic)/dels(2,ic)
          denom = denom + (1.-dels(2,ic))/dels(2,ic)*rdy2c
          dym = dels(2,ic)
          ppp = min(ppp,dels(2,ic))
        endif
c       --- the point higher in y is inside the conductor
        if (dels(3,ic) < 1.) then
          pyp = voltfac*volt(3,ic)/dels(3,ic)
          denom = denom + (1.-dels(3,ic))/dels(3,ic)*rdy2c
          dyp = dels(3,ic)
          ppp = min(ppp,dels(3,ic))
        endif
c       --- calculate the residual based on the boundary conditions
        if (ppp < 1.) then
c         rhovolume = (dxp+dxm)*(dyp+dym)*0.125
          rhovolume = 1.
          res(ix,iy) = ppp*(rho(ix,iy)*rhovolume
     &           + (pxm+pxp)*rdx2cs + (pym+pyp)*rdy2cs
     &           - phi(ix,iy)*mgparam*denom)
        endif
      enddo

      return
      end
c=============================================================================
      subroutine residualxy(nx,ny,dxsqi,dysqi,phi,rho,res,mglevel,
     &                      l2symtry,l4symtry,
     &                      lcndbndy,icndbndy,conductors)
      use GridBoundary3d
      use Multigrid3d
      use ConductorTypemodule
      integer(ISZ):: nx,ny
      real(kind=8):: dxsqi,dysqi
      real(kind=8):: phi(0:nx,0:ny),rho(0:nx,0:ny)
      real(kind=8):: res(0:nx,0:ny)
      integer(ISZ):: mglevel
      logical(ISZ):: l2symtry,l4symtry
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.
c taking into account the premultiplication of rho by
c   mgparam/(eps0*2.*(dxsqi+dysqi))
c The resulting residual is also implicitly multiplied by the same constant.
c Note that then for restriction of the residual to a coarser grid, it must
c be scaled by the ratio old(dxsqi+dysqi)/new(dxsqi+dysqi).
c This is done in the restrict routine automatically.
c
c For internal conductors, the residual is set to zero inside and calculated
c using the modified form of the finite differenced Poisson's equation near
c the surface.

      integer(ISZ):: ix,iy
      integer(ISZ):: ixmin,ixmax,iymin,iymax
      integer(ISZ):: ixm1,ixp1,iym1,iyp1
      real(kind=8):: const,dxsqic,dysqic
      const = 0.5/(dxsqi+dysqi)
      dxsqic = dxsqi*mgparam*const
      dysqic = dysqi*mgparam*const

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 1
      ixmax = nx-1
      iymin = 1
      iymax = ny-1
      if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
      if (boundxy == 1) ixmax = nx
      if (boundxy > 0 .or. l4symtry) iymin = 0
      if (boundxy == 1) iymax = ny

c     --- Calculate the residual.
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

          res(ix,iy) = rho(ix,iy)
     &        +  (phi(ixm1,iy  )+phi(ixp1,iy  ))*dxsqic
     &        +  (phi(ix  ,iym1)+phi(ix  ,iyp1))*dysqic
     &        -  phi(ix,iy)*mgparam

        enddo
      enddo

c     --- Make copies for the periodic boundaries.
      if (boundxy == 2) then
        res(ixmin:ixmax,ny) = res(ixmin:ixmax,0)
        res(nx,iymin:iymax) = res(0,iymin:iymax)
      endif

c     --- Zero the residual inside conductors.
      call cond_potmgxy(conductors%interior,nx,ny,res,mglevel,.true.)

c     --- Calculate the residual near the conductor.
      call condbndyresxy(conductors%evensubgrid,nx,ny,phi,rho,res,
     &                   dxsqi,dysqi,mgparam,boundxy,l2symtry,l4symtry,mglevel)
      call condbndyresxy(conductors%oddsubgrid,nx,ny,phi,rho,res,
     &                   dxsqi,dysqi,mgparam,boundxy,l2symtry,l4symtry,mglevel)

      return
      end
c=============================================================================
c=============================================================================
      RECURSIVE subroutine checkconductorsxy(nx,ny,dx,dy,l2symtry,l4symtry,
     &                                       conductors)
      use Parallel
      use ConductorTypemodule
      integer(ISZ):: nx,ny
      real(kind=8):: dx,dy
      logical(ISZ):: l2symtry,l4symtry
      type(ConductorType):: conductors

c Recursively calls the routine to generate the conductor data at the
c various mesh resolutions needed by the MG solver.

      integer(ISZ):: ic,ix,iy,ii,nn
      integer(ISZ),allocatable:: isort(:)

      integer(ISZ):: il,ic1,ie1,io1,ic2,ie2,io2,mglevels
      real(kind=8):: dx1,dy1
      integer(ISZ):: nx1,ny1
      integer(ISZ):: ix_axis1,iy_axis1

      mglevels = conductors%levels

c     --- Call the work routine for each level
      do il=0,conductors%levels-1
        dx1 = dx*conductors%levellx(il)
        dy1 = dy*conductors%levelly(il)
        nx1 = nx/conductors%levellx(il)
        ny1 = nx/conductors%levelly(il)
        call checkconductors_workxy(conductors%interior,
     &                            conductors%evensubgrid,conductors%oddsubgrid,
     &                            nx1,ny1,dx1,dy1,l2symtry,l4symtry,il)
      enddo

      conductors%interior%istart = 1
      conductors%evensubgrid%istart = 1
      conductors%oddsubgrid%istart = 1

c     --- Sort the conductor data by level number
c     --- First, sort conductor points
      if (conductors%interior%n > 0) then
        allocate(isort(conductors%interior%n))
        call isortconductor(conductors%interior%n,conductors%interior%ilevel,isort,conductors%interior%istart,
     &                      mglevels)
        nn = conductors%interior%n
        conductors%interior%n = conductors%interior%istart(mglevels) - 1
        call iswapconductor(conductors%interior%n,isort,nn,2,conductors%interior%indx)
        call iswapconductor(conductors%interior%n,isort,nn,1,conductors%interior%numb)
        call rswapconductor(conductors%interior%n,isort,nn,1,conductors%interior%volt)
        call iswapconductor(conductors%interior%n,isort,nn,1,conductors%interior%ilevel)
        deallocate(isort)
      endif

c     --- Sort even subgrid points
      if (conductors%evensubgrid%n > 0) then
        allocate(isort(conductors%evensubgrid%n))
        call isortconductor(conductors%evensubgrid%n,conductors%evensubgrid%ilevel,isort,
     &                      conductors%evensubgrid%istart,mglevels)
        nn = conductors%evensubgrid%n
        conductors%evensubgrid%n = conductors%evensubgrid%istart(mglevels) - 1
        call iswapconductor(conductors%evensubgrid%n,isort,nn,2,conductors%evensubgrid%indx)
        call rswapconductor(conductors%evensubgrid%n,isort,nn,4,conductors%evensubgrid%dels)
        call rswapconductor(conductors%evensubgrid%n,isort,nn,5,conductors%evensubgrid%volt)
        call iswapconductor(conductors%evensubgrid%n,isort,nn,5,conductors%evensubgrid%numb)
        call iswapconductor(conductors%evensubgrid%n,isort,nn,1,conductors%evensubgrid%ilevel)
        deallocate(isort)
      endif
    
c     --- Sort odd subgrid points
      if (conductors%oddsubgrid%n > 0) then
        allocate(isort(conductors%oddsubgrid%n))
        call isortconductor(conductors%oddsubgrid%n,conductors%oddsubgrid%ilevel,isort,
     &                      conductors%oddsubgrid%istart,mglevels)
        nn = conductors%oddsubgrid%n
        conductors%oddsubgrid%n = conductors%oddsubgrid%istart(mglevels) - 1
        call iswapconductor(conductors%oddsubgrid%n,isort,nn,2,conductors%oddsubgrid%indx)
        call rswapconductor(conductors%oddsubgrid%n,isort,nn,4,conductors%oddsubgrid%dels)
        call rswapconductor(conductors%oddsubgrid%n,isort,nn,5,conductors%oddsubgrid%volt)
        call iswapconductor(conductors%oddsubgrid%n,isort,nn,5,conductors%oddsubgrid%numb)
        call iswapconductor(conductors%oddsubgrid%n,isort,nn,1,conductors%oddsubgrid%ilevel)
        deallocate(isort)
      endif

      return
      end
c=============================================================================
      subroutine checkconductors_workxy(interior,evensubgrid,oddsubgrid,
     &                                  nx,ny,dx,dy,l2symtry,l4symtry,mglevel)
      use Multigrid3d
      use ConductorInteriorTypemodule
      use ConductorSubGridTypemodule
      type(ConductorInteriorType):: interior
      type(ConductorSubGridType):: evensubgrid,oddsubgrid
      integer(ISZ):: nx,ny,mglevel
      real(kind=8):: dx,dy
      logical(ISZ):: l2symtry,l4symtry

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

      integer(ISZ):: iii(0:nx,0:ny)
      integer(ISZ):: ic,i,ix,iy,id

c     --- Set the conductor points.
      iii = 0
      do ic=1,interior%n
        if (interior%ilevel(ic) /= mglevel) cycle
        ix = interior%indx(0,ic)
        iy = interior%indx(1,ic)
        if (ix < 0 .or. nx < ix .or.
     &      iy < 0 .or. ny < iy) then
          interior%ilevel(ic) = -1
          cycle
        endif
        iii(ix,iy) = interior%nmax + 1
      enddo

c     --- Scan through subgrid points:
c     ---   remove points which lie inside of a conductor
c     ---   register subgrid points in the work array iii
c     ---   check for redundant point (multiple points at grid location)
      do ic=1,evensubgrid%n

        if (evensubgrid%ilevel(ic) /= mglevel) cycle

        ix = evensubgrid%indx(0,ic)
        iy = evensubgrid%indx(1,ic)

        if (ix < 0 .or. nx < ix .or.
     &      iy < 0 .or. ny < iy) then
          evensubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy) == interior%nmax+1) then
          evensubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If iii == 0, then this data point is outside any conductors and
c       --- is not redundant.
        if (iii(ix,iy) == 0) then
          iii(ix,iy) = ic
          cycle
        endif

c       iii(ix,iy) < ncndmax+1
        i = iii(ix,iy)

        if (i < 0) then
c         --- The point already there is odd so must be a different level
          iii(ix,iy) = ic
          cycle
        endif

        if (evensubgrid%ilevel(ic) /= evensubgrid%ilevel(i)) then
c         --- The point already there is on a different level
          iii(ix,iy) = ic
          cycle
        endif

c       --- There is another subgrid point here. Combine the data
c       --- of the two points.
c       --- For each direction, check if a conductor is nearer to this point
c       --- than the other point. If so, reset data for the other point.
        do id=0,3
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

        if (ix < 0 .or. nx < ix .or.
     &      iy < 0 .or. ny < iy) then
          oddsubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If this point lies on a conductor point, kill it.
        if (iii(ix,iy) == interior%nmax+1) then
          oddsubgrid%ilevel(ic) = -1
          cycle
        endif

c       --- If iii == 0, then this data point is outside any conductors and
c       --- is not redundant.
        if (iii(ix,iy) == 0) then
          iii(ix,iy) = ic
          cycle
        endif

c       iii(ix,iy) < ncndmax+1
        i = iii(ix,iy)

        if (i < 0) then
c         --- The point already there is odd so must be a different level
          iii(ix,iy) = ic
          cycle
        endif

        if (oddsubgrid%ilevel(ic) /= oddsubgrid%ilevel(i)) then
c         --- The point already there is on a different level
          iii(ix,iy) = ic
          cycle
        endif

c       --- There is another subgrid point here. Combine the data
c       --- of the two points.
c       --- For each direction, check if a conductor is nearer to this point
c       --- than the other point. If so, reset data for the other point.
        do id=0,3
          if (oddsubgrid%dels(id,ic) < 1. .and.
     &        oddsubgrid%dels(id,ic) < oddsubgrid%dels(id,i)) then
            oddsubgrid%dels(id,i)  = oddsubgrid%dels(id,ic)
            oddsubgrid%volt(id,i) = oddsubgrid%volt(id,ic)
            oddsubgrid%numb(id,i) = oddsubgrid%numb(id,ic)
          endif
        enddo
        oddsubgrid%ilevel(ic) = -1

      enddo


      return
      end
c=============================================================================
