#include "top.h"
c=============================================================================
c@(#) File FXY_MGRID.M, version $Revision: 1.3 $, $Date: 2002/06/26 00:46:01 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 2D multigrid field sovler which is part of the FXY
c  package of WARP.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine multigridxyf(iwhich,nx,ny,dx,dy,phi1d,rho1d,
     &                       l2symtry,l4symtry,xmmin,ymmin)
      use Constant
      use PSOR3d
      use Multigrid3d
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny
      real(kind=8):: phi1d(*)
      real(kind=8):: rho1d(*)
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
      real(kind=8),allocatable:: res(:)

c     --- Initialize temporaries
      nxy  =(nx+1)*(ny+1)
      dxsqi=1./dx**2
      dysqi=1./dy**2
      reps0c = mgparam/(eps0*2.*(dxsqi+dysqi))

c     --- Allocate the res array.
      allocate(res(nxy))

c     --- Calculate the size of the work arrays. The work arrays must be
c     --- big enough to hold the phi and the rho at each level of coarseness,
c     --- including the finest. The size calculation include semi-coarsening.
c     --- Only do the allocation during initialization.
c     --- This loops over the grid levels, adding in the size of the grid
c     --- at each level.
      if (iwhich == 0 .or. iwhich == 1) then
        tnx = nx
        tny = ny
        tdx = dx
        tdy = dy
        tempsize = nxy
        do while (tnx > 2 .and. tny > 2 .and.
     &            mod(tnx,4) == 0 .and. mod(tny,4) == 0)
          tnx = tnx/2
          tny = tny/2
          tdx = tdx*2.
          tdy = tdy*2.
          tempsize = tempsize + (tnx+1)*(tny+1)
        enddo
        call gchange("Multigrid3d",0)
      endif

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

c     --- Determine the points that make up the conductor.  This takes extra
c     --- time and so should not be done if the grid is not moving in the lab
c     --- frame.  Set gridmode to 1 to avoid this call. The data is then
c     --- converted and expanded for the multigrid solver.
      if (gridmode == 0 .or. iwhich == -2) then
        call subgrid_sor_to_mgxy(nx,ny,dx,dy,l2symtry,l4symtry)
      endif

c     --- Copy phi to the temporary array and zero out rho_temp.
      call copyarry(phi1d(1),phi_temp(1),nxy)
      rho_temp(nxy+1:) = 0.

c     --- Preset rho to increase performance (reducing the number of
c     --- multiplies in the main SOR sweep loop). Also, copy rho1d
c     --- into rho_temp.
      rho_temp(1:nxy) = rho1d(1:nxy)*reps0c

c     --- Main multigrid v-cycle loop. Calculate error each iteration since
c     --- vary few iterations are done.
      mgiters = 0
      maxerr = 1.
      do while (maxerr > mgtol .and. mgiters < mgmaxiters)
        mgiters = mgiters + 1

c       --- Zero out the coarse grid data area.
        phi_temp(nxy+1:) = 0.

c       --- Do one vcycle.
        call vcyclexy(1,nx,ny,dxsqi,dysqi,phi_temp(1),rho_temp(1),res(1),
     &                l2symtry,l4symtry,boundxy)

        maxerr = maxval(abs(phi_temp(1:nxy) - phi1d(1:nxy)))

c       --- Save the current phi. The original phi1d array is used as the
c       --- scratch space to hold the current value of phi which is used
c       --- after the iteration to check for convergence.
        phi1d(1:nxy) = phi_temp(1:nxy)

      enddo

c     --- Make a print out.
      if (maxerr > mgtol) then
        print*,"Multigrid: Maximum number of iterations reached"
      endif
      print*,"Multigrid: Error converged to ",maxerr," in ",mgiters," v-cycles"

      return
      end
c=============================================================================
      RECURSIVE subroutine vcyclexy(levelxy,nx,ny,dxsqi,dysqi,
     &                            phi,rho,res,l2symtry,l4symtry,boundxy)
      use Multigrid3d
      integer(ISZ):: levelxy
      integer(ISZ):: nx,ny
      real(kind=8):: dxsqi,dysqi
      real(kind=8):: phi(*),rho(*),res(*)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: boundxy

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      integer(ISZ):: i2h,i
      real(kind=8):: ff

c     --- Do initial SOR passes.
      do i=1,downpasses
        call sorpassxy(levelxy,nx,ny,phi(1),rho(1),phi(1),rho(1),
     &                 dxsqi,dysqi,l2symtry,l4symtry)
      enddo

c     --- If grid is not at its coarsest level in any of the axis or and
c     --- all dimensions are even, continue the coarsening.
      if (nx > 2 .and. ny > 2 .and.
     &    mod(nx,4) == 0 .and. mod(ny,4) == 0) then

c       --- Location of the next block of data.
        i2h = (nx+1)*(ny+1) + 1

c       --- Get the residual on the current grid.
        call residualxy(nx,ny,dxsqi,dysqi,phi(1),rho(1),res(1),
     &                  levelxy,l2symtry,l4symtry)

c         --- Ratio of old to new constant needed to scale the residual for
c         --- the restriction.
c         ff = (dxsqi+dysqi)/(dxsqi*0.25 + dysqi*0.25)
          ff = 4.
          call restrictxy(nx,ny,res(1),rho(i2h),ff,boundxy,
     &                    l2symtry,l4symtry)

c         --- Continue at the next coarsest level.
          call vcyclexy(levelxy*2,nx/2,ny/2,
     &                  dxsqi*0.25,dysqi*0.25,phi(i2h),rho(i2h),res(1),
     &                  l2symtry,l4symtry,boundxy)

c         --- Add in resulting error.
          call expandxy(nx/2,ny/2,phi(i2h),phi(1),boundxy)

      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call sorpassxy(levelxy,nx,ny,phi(1),rho(1),phi(1),rho(1),
     &                 dxsqi,dysqi,l2symtry,l4symtry)
      enddo

      return
      end
c=============================================================================
      subroutine restrictxy(nx,ny,ph,p2h,ff,boundxy,l2symtry,l4symtry)
      integer(ISZ):: nx,ny
      real(kind=8):: ph(0:nx,0:ny)
      real(kind=8):: p2h(0:nx/2,0:ny/2)
      real(kind=8):: ff
      integer(ISZ):: boundxy
      logical(ISZ):: l2symtry,l4symtry

c Restrict transversely to a coarser grid.  The factor of ff is needed
c since the residual will be used as the source term in the next coarser
c grid and it implicitly includes the factor 0.5/(dxsqi+dysqi)
c which is ff times too small with the current values of dxsqi etc.

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
     &      ff*0.2500*ph(ix  ,iy  ) +
     &      ff*0.1250*ph(ixm1,iy  ) +
     &      ff*0.1250*ph(ixp1,iy  ) +
     &      ff*0.1250*ph(ix  ,iym1) +
     &      ff*0.1250*ph(ix  ,iyp1) +
     &      ff*0.0625*ph(ixm1,iym1) +
     &      ff*0.0625*ph(ixp1,iym1) +
     &      ff*0.0625*ph(ixm1,iyp1) +
     &      ff*0.0625*ph(ixp1,iyp1)

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
      subroutine sorpassxy(levelxy,nx,ny,phi,rho,phi1d,rho1d,
     &                     rdx2,rdy2,l2symtry,l4symtry)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      use Constant
      integer(ISZ):: levelxy,nx,ny
      real(kind=8):: phi(0:nx,0:ny),rho(0:nx,0:ny)
      real(kind=8):: phi1d(*),rho1d(*)
      real(kind=8):: rdx2,rdy2
      logical(ISZ):: l2symtry,l4symtry

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
      integer(ISZ):: ii,ix,iy,ic
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
      call cond_potmgxy(nx,ny,phi(0,0),levelxy,.false.)

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
        do ic=1,necndbdy
          if (iecndlevel(ic) >= levelxy) then
            ecndpvph(ic) = phi(iecndx(ic)/levelxy,
     &                         iecndy(ic)/levelxy)
          endif
        enddo
        if (levelxy == 1) then
          do ic=1,nocndbdy
            if (iocndlevel(ic) >= levelxy) then
              ocndpvph(ic) = phi(iocndx(ic)/levelxy,
     &                           iocndy(ic)/levelxy)
            endif
          enddo
        endif
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
        if (lcndbndy .and. ncndmax > 0) then
          if (parity == 0 .or. levelxy > 1) then
            call condbndymgxy(nx,ny,phi(0,0),rho(0,0),rdx2c,rdy2c,spm1,mgparam,
     &                        ncndmax,necndbdy,ecdelmx(1),ecdelmy(1),
     &                        ecdelpx(1),ecdelpy(1),iecndx(1),iecndy(1),
     &                        iecndlevel(1),
     &                        ecvoltmx(1),ecvoltpx(1),ecvoltmy(1),ecvoltpy(1),
     &                        ecndpvph(1),boundxy,l2symtry,l4symtry,
     &                        levelxy,parity)
          endif
          if (parity == 1 .and. levelxy == 1) then
            call condbndymgxy(nx,ny,phi(0,0),rho(0,0),rdx2c,rdy2c,spm1,mgparam,
     &                        ncndmax,nocndbdy,ocdelmx(1),ocdelmy(1),
     &                        ocdelpx(1),ocdelpy(1),iocndx(1),iocndy(1),
     &                        iocndlevel(1),
     &                        ocvoltmx(1),ocvoltpx(1),ocvoltmy(1),ocvoltpy(1),
     &                        ocndpvph(1),boundxy,l2symtry,l4symtry,
     &                        levelxy,parity)
          endif
        endif

c       --- Put desired potential onto conductors in phi array.
        call cond_potmgxy(nx,ny,phi(0,0),levelxy,.false.)

c     --- end of loop over even and odd points
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine cond_potmgxy(nx,ny,phi,levelxy,lresidual)
      use PSOR3d
      use Conductor3d
      integer(ISZ):: nx,ny,levelxy
      real(kind=8):: phi(0:nx,0:ny)
      logical(ISZ):: lresidual

c Set conductor points to the desired potential.

      integer(ISZ):: ic

c     --- When at the finest level and not calculating the residual, set
c     --- phi to the voltage of the conductor, otherwise, set it to zero.
      if (levelxy == 1 .and. .not. lresidual) then
        do ic = 1,ncond
          phi(ixcond(ic),iycond(ic)) = condvolt(ic)
        enddo
      else
        do ic=1,ncond
c         --- Only set points which are on the coarse grid at the current
c         --- level of refinement.
          if (icondlevel(ic) >= levelxy) then
            phi(ixcond(ic)/levelxy,iycond(ic)/levelxy) = 0.
          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine condbndymgxy(nx,ny,phi,rho,rdx2c,rdy2c,spm1,srp,
     &                      ncndmax,ncndbdy,cdelmx,cdelmy,cdelpx,cdelpy,
     &                      icndx,icndy,icndlevel,
     &                      cvoltmx,cvoltpx,cvoltmy,cvoltpy,
     &                      cndpvph,boundxy,
     &                      l2symtry,l4symtry,levelxy,parity)
      integer(ISZ):: nx,ny,ncndmax,ncndbdy,levelxy,parity
      real(kind=8):: phi(0:nx,0:ny), rho(0:nx,0:ny)
      real(kind=8):: rdx2c,rdy2c,spm1,srp
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax)
      integer(ISZ):: icndlevel(ncndmax)
      real(kind=8):: cvoltmx(ncndmax),cvoltmy(ncndmax)
      real(kind=8):: cvoltpx(ncndmax),cvoltpy(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
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
c to hold icndx-1, icndx+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: rdx2cos,rdy2cos,pxm,pym,pxp,pyp,denom
      real(kind=8):: dxm,dxp,dym,dyp
      real(kind=8):: voltfac,rhovolume
      real(kind=8):: levelxyi
      integer(ISZ):: ic,icndxp1,icndxm1,icndyp1,icndym1
      integer(ISZ):: icdx,icdy
      logical(ISZ):: dosubgrid

      rdx2cos = rdx2c/srp
      rdy2cos = rdy2c/srp
      levelxyi = 1./levelxy

c     --- loop over points near surface of conductors
      do ic = 1,ncndbdy

c       --- Only use point if it is usable at this level or higher
c       --- and if the parity is correct.
        if (icndlevel(ic) >= levelxy) then

          icdx = icndx(ic)/levelxy
          icdy = icndy(ic)/levelxy

          if (mod(icdx+icdy,2)==parity) then

c           --- set temporaries for boundaries
            icndxp1 = icdx + 1
            icndxm1 = icdx - 1
            icndyp1 = icdy + 1
            icndym1 = icdy - 1

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

c           --- Set temporaries with initial values.
            pxm = phi(icndxm1,icdy   )
            pxp = phi(icndxp1,icdy   )
            pym = phi(icdx   ,icndym1)
            pyp = phi(icdx   ,icndyp1)
            denom = 1.
            dosubgrid = .false.
            dxm = 1.
            dxp = 1.
            dym = 1.
            dyp = 1.

c           --- Only use actual voltage on finest level. Set to zero for
c           --- coarser levels since solver for the residuals.
            if (levelxy == 1) then
              voltfac = 1.
            else
              voltfac = 0.
            endif

c           --- the point lower in x is inside the conductor
            if (cdelmx(ic) < levelxy) then
              pxm = voltfac*cvoltmx(ic)/cdelmx(ic)
              denom = denom + (levelxy-cdelmx(ic))/cdelmx(ic)*rdx2cos
              dxm = cdelmx(ic)*levelxyi
              dosubgrid = .true.
            endif
c           --- the point higher in x is inside the conductor
            if (cdelpx(ic) < levelxy) then
              pxp = voltfac*cvoltpx(ic)/cdelpx(ic)
              denom = denom + (levelxy-cdelpx(ic))/cdelpx(ic)*rdx2cos
              dxp = cdelpx(ic)*levelxyi
              dosubgrid = .true.
            endif
c           --- the point lower in y is inside the conductor
            if (cdelmy(ic) < levelxy) then
              pym = voltfac*cvoltmy(ic)/cdelmy(ic)
              denom = denom + (levelxy-cdelmy(ic))/cdelmy(ic)*rdy2cos
              dym = cdelmy(ic)*levelxyi
              dosubgrid = .true.
            endif
c           --- the point higher in y is inside the conductor
            if (cdelpy(ic) < levelxy) then
              pyp = voltfac*cvoltpy(ic)/cdelpy(ic)
              denom = denom + (levelxy-cdelpy(ic))/cdelpy(ic)*rdy2cos
              dyp = cdelpy(ic)*levelxyi
              dosubgrid = .true.
            endif
c           --- calculate the new phi based on the boundary conditions
            if (dosubgrid) then
c             rhovolume = (dxp+dxm)*(dyp+dym)*0.125
              rhovolume = 1.
              phi(icdx,icdy) = (rho(icdx,icdy)*rhovolume +
     &          (pxm+pxp)*rdx2c + (pym+pyp)*rdy2c)/denom +
     &          spm1*cndpvph(ic)
            endif
          endif
        endif
      enddo

      return
      end
c=============================================================================
      subroutine condbndyresxy(nx,ny,phi,rho,res,rdx2,rdy2,mgparam,
     &                       ncndmax,ncndbdy,cdelmx,cdelmy,cdelpx,cdelpy,
     &                       icndx,icndy,icndlevel,
     &                       cvoltmx,cvoltpx,cvoltmy,cvoltpy,
     &                       cndpvph,boundxy,
     &                       l2symtry,l4symtry,levelxy,parity)
      integer(ISZ):: nx,ny,ncndmax,ncndbdy,levelxy,parity
      real(kind=8):: phi(0:nx,0:ny)
      real(kind=8):: rho(0:nx,0:ny), res(0:nx,0:ny)
      real(kind=8):: rdx2,rdy2,mgparam
      real(kind=8):: cdelmx(ncndmax), cdelmy(ncndmax)
      real(kind=8):: cdelpx(ncndmax), cdelpy(ncndmax)
      integer(ISZ):: icndx(ncndmax), icndy(ncndmax)
      integer(ISZ):: icndlevel(ncndmax)
      real(kind=8):: cvoltmx(ncndmax),cvoltmy(ncndmax)
      real(kind=8):: cvoltpx(ncndmax),cvoltpy(ncndmax)
      real(kind=8):: cndpvph(ncndmax)
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
c to hold icndx-1, icndx+1 etc which are are adjusted appopriately for
c points on the boundary.

      real(kind=8):: const,rdx2c,rdy2c,pxm,pym,pxp,pyp,denom
      real(kind=8):: dxm,dxp,dym,dyp,rhovolume
      real(kind=8):: voltfac
      real(kind=8):: rdx2cs,rdy2cs,ppp
      integer(ISZ):: ic,icndxp1,icndxm1,icndyp1,icndym1
      integer(ISZ):: icdx,icdy
      real(kind=8):: levelxyi

      const = 0.5/(rdx2+rdy2)
      rdx2c = rdx2*const
      rdy2c = rdy2*const
      rdx2cs = mgparam*rdx2*const
      rdy2cs = mgparam*rdy2*const
      levelxyi = 1./levelxy

c     --- loop over points near surface of conductors
      do ic = 1,ncndbdy

c       --- Only use point if it is usable at this level or higher
c       --- and if the parity is correct.
        if (icndlevel(ic) >= levelxy) then

          icdx = icndx(ic)/levelxy
          icdy = icndy(ic)/levelxy

          if (mod(icdx+icdy,2)==parity) then

c           --- set temporaries for boundaries
            icndxp1 = icdx + 1
            icndxm1 = icdx - 1
            icndyp1 = icdy + 1
            icndym1 = icdy - 1

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

c           --- set temporaries with initial values
            pxm = phi(icndxm1,icdy   )
            pxp = phi(icndxp1,icdy   )
            pym = phi(icdx   ,icndym1)
            pyp = phi(icdx   ,icndyp1)
            denom = 1.
            ppp = 1.
            dxm = 1.
            dxp = 1.
            dym = 1.
            dyp = 1.

c           --- Only use actual voltage on finest level. Set to zero for
c           --- coarser levels since solver for the residuals.
            if (levelxy == 1) then
              voltfac = 1.
            else
              voltfac = 0.
            endif

c           --- the point lower in x is inside the conductor
            if (cdelmx(ic) < levelxy) then
              pxm = voltfac*cvoltmx(ic)/cdelmx(ic)
              denom = denom + (levelxy-cdelmx(ic))/cdelmx(ic)*rdx2c
              dxm = cdelmx(ic)*levelxyi
              ppp = min(ppp,cdelmx(ic)*levelxyi)
            endif
c           --- the point higher in x is inside the conductor
            if (cdelpx(ic) < levelxy) then
              pxp = voltfac*cvoltpx(ic)/cdelpx(ic)
              denom = denom + (levelxy-cdelpx(ic))/cdelpx(ic)*rdx2c
              dxp = cdelpx(ic)*levelxyi
              ppp = min(ppp,cdelpx(ic)*levelxyi)
            endif
c           --- the point lower in y is inside the conductor
            if (cdelmy(ic) < levelxy) then
              pym = voltfac*cvoltmy(ic)/cdelmy(ic)
              denom = denom + (levelxy-cdelmy(ic))/cdelmy(ic)*rdy2c
              dym = cdelmy(ic)*levelxyi
              ppp = min(ppp,cdelmy(ic)*levelxyi)
            endif
c           --- the point higher in y is inside the conductor
            if (cdelpy(ic) < levelxy) then
              pyp = voltfac*cvoltpy(ic)/cdelpy(ic)
              denom = denom + (levelxy-cdelpy(ic))/cdelpy(ic)*rdy2c
              dyp = cdelpy(ic)*levelxyi
              ppp = min(ppp,cdelpy(ic)*levelxyi)
            endif
c           --- calculate the residual based on the boundary conditions
            if (ppp < 1.) then
c             rhovolume = (dxp+dxm)*(dyp+dym)*0.125
              rhovolume = 1.
              res(icdx,icdy) = ppp*(rho(icdx,icdy)*rhovolume
     &               + (pxm+pxp)*rdx2cs + (pym+pyp)*rdy2cs
     &               - phi(icdx,icdy)*mgparam*denom)
            endif
          endif
        endif
      enddo

      return
      end
c=============================================================================
      subroutine residualxy(nx,ny,dxsqi,dysqi,phi,rho,res,levelxy,
     &                      l2symtry,l4symtry)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny
      real(kind=8):: dxsqi,dysqi
      real(kind=8):: phi(0:nx,0:ny),rho(0:nx,0:ny)
      real(kind=8):: res(0:nx,0:ny)
      integer(ISZ):: levelxy
      logical(ISZ):: l2symtry,l4symtry

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
      call cond_potmgxy(nx,ny,res(0,0),levelxy,.true.)

      if (lcndbndy .and. ncndmax > 0) then
c     --- Calculate the residual near the conductor. Note that the odd points
c     --- are only used at the finest level. Also, for higher levels, the even
c     --- points are passed down with an odd parity since they may be odd at
c     --- that level.
      if (levelxy == 1) then
        call condbndyresxy(nx,ny,phi(0,0),rho(0,0),res(0,0),dxsqi,dysqi,mgparam,
     &             ncndmax,nocndbdy,ocdelmx(1),ocdelmy(1),ocdelpx(1),ocdelpy(1),
     &             iocndx(1),iocndy(1),iocndlevel(1),
     &             ocvoltmx(1),ocvoltpx(1),ocvoltmy(1),ocvoltpy(1),
     &             ocndpvph(1),boundxy,
     &             l2symtry,l4symtry,levelxy,1)
      endif
      if (levelxy > 1) then
        call condbndyresxy(nx,ny,phi(0,0),rho(0,0),res(0,0),dxsqi,dysqi,mgparam,
     &             ncndmax,necndbdy,ecdelmx(1),ecdelmy(1),ecdelpx(1),ecdelpy(1),
     &             iecndx(1),iecndy(1),iecndlevel(1),
     &             ecvoltmx(1),ecvoltpx(1),ecvoltmy(1),ecvoltpy(1),
     &             ecndpvph(1),boundxy,
     &             l2symtry,l4symtry,levelxy,1)
      endif
      call condbndyresxy(nx,ny,phi(0,0),rho(0,0),res(0,0),dxsqi,dysqi,mgparam,
     &             ncndmax,necndbdy,ecdelmx(1),ecdelmy(1),ecdelpx(1),ecdelpy(1),
     &             iecndx(1),iecndy(1),iecndlevel(1),
     &             ecvoltmx(1),ecvoltpx(1),ecvoltmy(1),ecvoltpy(1),
     &             ecndpvph(1),boundxy,
     &             l2symtry,l4symtry,levelxy,0)

      endif

      return
      end
c=============================================================================
      subroutine conductor_data_levelxy(nx,ny,dx,dy)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny
      real(kind=8):: dx,dy

c Set the coarsest level at which each conductor point is still on the grid.
c Also, remove subgrid data points which are not used - those points which are
c not near a conducting surface at any level of refinement.

      integer(ISZ):: ic,ne,no,levelxy
      integer(ISZ):: tnx,tny
      real(kind=8):: tdx,tdy

c     --- If there are not conductors, then immediately return.
      if (ncond == 0 .and. necndbdy == 0 .and. nocndbdy == 0) return

c     --- Zero arrays, assuming that no points are used.
      icondlevel(:ncond) = 0
      iecndlevel(:necndbdy) = 0
      iocndlevel(:nocndbdy) = 0

c     --- Initialize the temps
      tnx = nx
      tny = ny
      tdx = dx
      tdy = dy
      levelxy = 1

c     --- Loop until at least one of the axis is at the coarsest level.
      do while (tnx >= 2 .and. tny >= 2 .and.
     &          mod(tnx,4) == 0 .and. mod(tny,4) == 0)

c       --- Interior points
        do ic=1,ncond
c         --- Check if point is on the current grid.
          if (mod(ixcond(ic),levelxy) == 0 .and.
     &        mod(iycond(ic),levelxy) == 0) then
            icondlevel(ic) = levelxy
          endif
        enddo

c       --- Even points
        do ic=1,necndbdy
c         --- Check if point is near a conductor.
c         --- i.e. is within one grid cell of a conductor on the current grid.
          if (ecdelmx(ic)/levelxy < 1. .or. ecdelpx(ic)/levelxy < 1. .or.
     &        ecdelmy(ic)/levelxy < 1. .or. ecdelpy(ic)/levelxy < 1.) then
c           --- Now, check if point is on the current grid.
            if (mod(iecndx(ic),levelxy) == 0 .and.
     &          mod(iecndy(ic),levelxy) == 0) then
              iecndlevel(ic) = levelxy
            endif
          endif
        enddo

c       --- Odd points
        do ic=1,nocndbdy
c         --- Check if point is near a conductor.
c         --- i.e. is within one grid cell of a conductor on the current grid.
          if (ocdelmx(ic)/levelxy < 1. .or. ocdelpx(ic)/levelxy < 1. .or.
     &        ocdelmy(ic)/levelxy < 1. .or. ocdelpy(ic)/levelxy < 1.) then
c           --- Now, check if point is on the current grid.
            if (mod(iocndx(ic),levelxy) == 0 .and.
     &          mod(iocndy(ic),levelxy) == 0) then
              iocndlevel(ic) = levelxy
            endif
          endif
        enddo

c       --- Go to the next level of coarseness.
        tnx = tnx/2
        tny = tny/2
        tdx = tdx*2.
        tdy = tdy*2.
        levelxy = 2*levelxy

      enddo

c     --- Now, clear out all of the points which were not flagged as usable.
      ne = 0
      do ic=1,necndbdy
        if (iecndlevel(ic) > 0) then
          ne = ne + 1
          if (ne < ic) then
            iecndx(ne)   = iecndx(ic)
            iecndy(ne)   = iecndy(ic)
            ecdelmx(ne)  = ecdelmx(ic)
            ecdelmy(ne)  = ecdelmy(ic)
            ecdelpx(ne)  = ecdelpx(ic)
            ecdelpy(ne)  = ecdelpy(ic)
            ecvolt(ne)   = ecvolt(ic)
            ecvoltmx(ne) = ecvoltmx(ic)
            ecvoltpx(ne) = ecvoltpx(ic)
            ecvoltmy(ne) = ecvoltmy(ic)
            ecvoltpy(ne) = ecvoltpy(ic)
            iecndlevel(ne) = iecndlevel(ic)
          endif
        endif
      enddo
      necndbdy = ne

      no = 0
      do ic=1,nocndbdy
        if (iocndlevel(ic) > 0) then
          no = no + 1
          if (no < ic) then
            iocndx(no)   = iocndx(ic)
            iocndy(no)   = iocndy(ic)
            ocdelmx(no)  = ocdelmx(ic)
            ocdelmy(no)  = ocdelmy(ic)
            ocdelpx(no)  = ocdelpx(ic)
            ocdelpy(no)  = ocdelpy(ic)
            ocvolt(no)   = ocvolt(ic)
            ocvoltmx(no) = ocvoltmx(ic)
            ocvoltpx(no) = ocvoltpx(ic)
            ocvoltmy(no) = ocvoltmy(ic)
            ocvoltpy(no) = ocvoltpy(ic)
            iocndlevel(no) = iocndlevel(ic)
          endif
        endif
      enddo
      nocndbdy = no

      return
      end
c=============================================================================
      subroutine subgrid_sor_to_mgxy(nx,ny,dx,dy,l2symtry,l4symtry)
      use PSOR3d
      use Conductor3d
      integer(ISZ):: nx,ny
      real(kind=8):: dx,dy
      logical(ISZ):: l2symtry,l4symtry

c Convert subgrid data from the format used for the SOR field solver to the
c the format used for multigrid. This includes expanding the data set for
c use on the coarser grid levels. Note that after the conversion, the data
c set will still work for the SOR (though it will be slightly less
c efficient since there will be extra points in the data set that the SOR
c doesn't need).

      integer(ISZ):: ix,iy

c     --- If there are not conductors, then immediately return.
      if (ncond == 0 .and. necndbdy == 0 .and. nocndbdy == 0) return

c     --- Make sure there is enough room for the extra data.
c     --- The factor of 3 is an over estimate since the number of points
c     --- should typically only double.
      if (ncndmax < 3*necndbdy) then
        ncndmax = 3*necndbdy
        call gchange("Conductor3d",0)
      endif

c     --- Do the conversion.
      call subgrid_sor_to_mgxy_work(nx,ny,dx,dy,l2symtry,l4symtry)

c     --- Now, get the refinement level for all of the points.
      call conductor_data_levelxy(nx,ny,dx,dy)

      return
      end
c=============================================================================
      subroutine subgrid_sor_to_mgxy_work(nx,ny,dx,dy,l2symtry,l4symtry)
      use PSOR3d
      use Conductor3d
      use Multigrid3d
      integer(ISZ):: nx,ny
      real(kind=8):: dx,dy
      logical(ISZ):: l2symtry,l4symtry

c This does the actual work of converting the subgrid data from the SOR format
c to the multigrid format. To do the conversion if does the following...
c  - register the conducting points on a 2D grid (this just makes it easy
c    to check for the presence of the conductor.
c  - clean up data set, removing any subgrid points which may lie inside
c    of a conductor (those points are harmless to SOR but are damaging
c    to multigrid)
c  - for existing data set, calculate fully the distances to conductors which
c    are greater than the grid cell size.
c  - register the subgrid points on the 2D grid
c  - finally, scan through the grid at each coarse level of resolution, and
c    save points which are within on grid cell of a conductor
c
c The notation for the 2D work grid is...
c   - all of the points inside of conductors are given the same large integer
c     value =((nx+1)*(ny+1) + 1)
c   - for all subgrid points, the index of that point is stored with a sign
c     attached - positive for even points, negative for odd points
c   - the value chosen for inside of conductors is large enough so that
c     is will not be the same as a subgrid point

      integer(ISZ):: icond,ix,iy
      integer(ISZ):: cds(2),maxcdel
      integer(ISZ):: ic,ne,no,levelxy,tnx,tny
      integer(ISZ):: ixmin,ixmax,iymin,iymax
      integer(ISZ):: iii(0:nx,0:ny)
      real(kind=8):: tdx,tdy

c     --- If there are not conductors, then immediately return.
      if (ncond == 0 .and. necndbdy == 0 .and. nocndbdy == 0) return

c     --- Make sure that iii is zeroed out.
      iii = 0

c     --- Set the known conductor points.
      icond = (nx+1)*(ny+1) + 1
      do ic=1,ncond
        iii(ixcond(ic),iycond(ic)) = icond + ic
      enddo

c     --- Scan through subgrid points and remove all points which lie inside
c     --- of a conductor.
      ne = 0
      do ic=1,necndbdy
        if (iii(iecndx(ic),iecndy(ic)) == 0) then
          ne = ne + 1
          if (ne < ic) then
            iecndx(ne)   = iecndx(ic)
            iecndy(ne)   = iecndy(ic)
            ecdelmx(ne)  = ecdelmx(ic)
            ecdelmy(ne)  = ecdelmy(ic)
            ecdelpx(ne)  = ecdelpx(ic)
            ecdelpy(ne)  = ecdelpy(ic)
            ecvolt(ne)   = ecvolt(ic)
            ecnumb(ne)   = ecnumb(ic)
          endif
        endif
      enddo
      necndbdy = ne

      no = 0
      do ic=1,nocndbdy
        if (iii(iocndx(ic),iocndy(ic)) == 0) then
          no = no + 1
          if (no < ic) then
            iocndx(no)   = iocndx(ic)
            iocndy(no)   = iocndy(ic)
            ocdelmx(no)  = ocdelmx(ic)
            ocdelmy(no)  = ocdelmy(ic)
            ocdelpx(no)  = ocdelpx(ic)
            ocdelpy(no)  = ocdelpy(ic)
            ocvolt(no)   = ocvolt(ic)
            ocnumb(no)   = ocnumb(ic)
          endif
        endif
      enddo
      nocndbdy = no

c     --- Reset distances which are greater than one to LARGEPOS and
c     --- register subgrid points in the work array iii. Also copy the
c     --- voltage and conductor number into the appropriate arrays.
      do ic=1,necndbdy
        iii(iecndx(ic),iecndy(ic)) = ic
        if (ecdelmx(ic) > 1.) ecdelmx(ic) = LARGEPOS
        if (ecdelpx(ic) > 1.) ecdelpx(ic) = LARGEPOS
        if (ecdelmy(ic) > 1.) ecdelmy(ic) = LARGEPOS
        if (ecdelpy(ic) > 1.) ecdelpy(ic) = LARGEPOS
        if (ecdelmx(ic) < 1.) then
          ecvoltmx(ic) = ecvolt(ic)
          ecnumbmx(ic) = ecnumb(ic)
        endif
        if (ecdelpx(ic) < 1.) then
          ecvoltpx(ic) = ecvolt(ic)
          ecnumbpx(ic) = ecnumb(ic)
        endif
        if (ecdelmy(ic) < 1.) then
          ecvoltmy(ic) = ecvolt(ic)
          ecnumbmy(ic) = ecnumb(ic)
        endif
        if (ecdelpy(ic) < 1.) then
          ecvoltpy(ic) = ecvolt(ic)
          ecnumbpy(ic) = ecnumb(ic)
        endif
      enddo
      do ic=1,nocndbdy
        iii(iocndx(ic),iocndy(ic)) = -ic
        if (ocdelmx(ic) > 1.) ocdelmx(ic) = LARGEPOS
        if (ocdelpx(ic) > 1.) ocdelpx(ic) = LARGEPOS
        if (ocdelmy(ic) > 1.) ocdelmy(ic) = LARGEPOS
        if (ocdelpy(ic) > 1.) ocdelpy(ic) = LARGEPOS
        if (ocdelmx(ic) < 1.) then
          ocvoltmx(ic) = ocvolt(ic)
          ocnumbmx(ic) = ocnumb(ic)
        endif
        if (ocdelpx(ic) < 1.) then
          ocvoltpx(ic) = ocvolt(ic)
          ocnumbpx(ic) = ocnumb(ic)
        endif
        if (ocdelmy(ic) < 1.) then
          ocvoltmy(ic) = ocvolt(ic)
          ocnumbmy(ic) = ocnumb(ic)
        endif
        if (ocdelpy(ic) < 1.) then
          ocvoltpy(ic) = ocvolt(ic)
          ocnumbpy(ic) = ocnumb(ic)
        endif
      enddo

c     --- Check through all of the SOR points to reset any distances which
c     --- should not be LARGEPOS.
      maxcdel = max(nx,ny) + 1
      do ic=1,necndbdy
        cds(1) = iecndx(ic)
        cds(2) = iecndy(ic)
        call subgrid_scan_axisxy(1,-1,cds,  -1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelmx(1),ocdelmx(1),condvolt(1),ecvoltmx(1),ocvoltmx(1),
     &               condnumb(1),ecnumbmx(1),ocnumbmx(1),
     &               iecndx(ic),ecdelmx(ic),ecvoltmx(ic),ecnumbmx(ic))
        call subgrid_scan_axisxy(1,+1,cds,nx+1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelpx(1),ocdelpx(1),condvolt(1),ecvoltpx(1),ocvoltpx(1),
     &               condnumb(1),ecnumbpx(1),ocnumbpx(1),
     &               iecndx(ic),ecdelpx(ic),ecvoltpx(ic),ecnumbpx(ic))
        call subgrid_scan_axisxy(2,-1,cds,  -1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelmy(1),ocdelmy(1),condvolt(1),ecvoltmy(1),ocvoltmy(1),
     &               condnumb(1),ecnumbmy(1),ocnumbmy(1),
     &               iecndy(ic),ecdelmy(ic),ecvoltmy(ic),ecnumbmy(ic))
        call subgrid_scan_axisxy(2,+1,cds,ny+1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelpy(1),ocdelpy(1),condvolt(1),ecvoltpy(1),ocvoltpy(1),
     &               condnumb(1),ecnumbpy(1),ocnumbpy(1),
     &               iecndy(ic),ecdelpy(ic),ecvoltpy(ic),ecnumbpy(ic))
      enddo
      do ic=1,nocndbdy
        cds(1) = iocndx(ic)
        cds(2) = iocndy(ic)
        call subgrid_scan_axisxy(1,-1,cds,  -1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelmx(1),ocdelmx(1),condvolt(1),ecvoltmx(1),ocvoltmx(1),
     &               condnumb(1),ecnumbmx(1),ocnumbmx(1),
     &               iocndx(ic),ocdelmx(ic),ocvoltmx(ic),ocnumbmx(ic))
        call subgrid_scan_axisxy(1,+1,cds,nx+1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelpx(1),ocdelpx(1),condvolt(1),ecvoltpx(1),ocvoltpx(1),
     &               condnumb(1),ecnumbpx(1),ocnumbpx(1),
     &               iocndx(ic),ocdelpx(ic),ocvoltpx(ic),ocnumbpx(ic))
        call subgrid_scan_axisxy(2,-1,cds,  -1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelmy(1),ocdelmy(1),condvolt(1),ecvoltmy(1),ocvoltmy(1),
     &               condnumb(1),ecnumbmy(1),ocnumbmy(1),
     &               iocndy(ic),ocdelmy(ic),ocvoltmy(ic),ocnumbmy(ic))
        call subgrid_scan_axisxy(2,+1,cds,ny+1,maxcdel,icond,nx,ny,iii(0,0),
     &               ecdelpy(1),ocdelpy(1),condvolt(1),ecvoltpy(1),ocvoltpy(1),
     &               condnumb(1),ecnumbpy(1),ocnumbpy(1),
     &               iocndy(ic),ocdelpy(ic),ocvoltpy(ic),ocnumbpy(ic))
      enddo


c     --- Now, scan the grid at each level of refinement to gather points
c     --- which are within one grid cell of a conductor are the coarser
c     --- grid levels.

c     --- Initialize the temps
      tnx = nx
      tny = ny
      tdx = dx
      tdy = dy
      levelxy = 1

c     --- Increase size of conductor arrays if needed. Increase in size
c     --- is only an estimate.
      if (necndbdy == ncndmax .or. nocndbdy == ncndmax) then
        ncndmax = 1.5*max(necndbdy,nocndbdy) + 1
        call gchange("Conductor3d",0)
      endif

c     --- Loop until at least one of the axis is at the coarsest level.
      ic = necndbdy + 1
      do while (tnx > 2 .and. tny > 2 .and.
     &          mod(tnx,4) == 0 .and. mod(tny,4) == 0)

c       --- Go to the next level of coarseness.
        tnx = tnx/2
        tny = tny/2
        tdx = tdx*2.
        tdy = tdy*2.
        levelxy = 2*levelxy

c       --- Loop over current grid and scan in each direction to find the
c       --- nearest conductor. If the resulting distances put it near a
c       --- conductor on the coarse grid, then add it to the list of even
c       --- points. (Odd points only appear on the finest level grid.)
        ixmin = levelxy
        ixmax = nx-levelxy
        iymin = levelxy
        iymax = ny-levelxy
        if (boundxy > 0 .or. l2symtry .or. l4symtry) ixmin = 0
        if (boundxy == 1) ixmax = nx
        if (boundxy > 0 .or. l4symtry) iymin = 0
        if (boundxy == 1) iymax = ny
        do iy=iymin,iymax,levelxy
          cds(2) = iy
          do ix=ixmin,ixmax,levelxy
            cds(1) = ix

c           --- Check if there is already a conductor or subgrid point
c           --- at this location.
            if (iii(ix,iy) == 0) then
              ecdelmx(ic) = LARGEPOS
              ecdelpx(ic) = LARGEPOS
              ecdelmy(ic) = LARGEPOS
              ecdelpy(ic) = LARGEPOS

c             --- Scan along the 6 directions.
              call subgrid_scan_axisxy(1,-1,cds,  -1,maxcdel,icond,nx,ny,
     &                               iii(0,0),ecdelmx(1),ocdelmx(1),
     &                               condvolt(1),ecvoltmx(1),ocvoltmx(1),
     &                               condnumb(1),ecnumbmx(1),ocnumbmx(1),
     &                               ix,ecdelmx(ic),ecvoltmx(ic),ecnumbmx(ic))
              call subgrid_scan_axisxy(1,+1,cds,nx+1,maxcdel,icond,nx,ny,
     &                               iii(0,0),ecdelpx(1),ocdelpx(1),
     &                               condvolt(1),ecvoltpx(1),ocvoltpx(1),
     &                               condnumb(1),ecnumbpx(1),ocnumbpx(1),
     &                               ix,ecdelpx(ic),ecvoltpx(ic),ecnumbpx(ic))
              call subgrid_scan_axisxy(2,-1,cds,  -1,maxcdel,icond,nx,ny,
     &                               iii(0,0),ecdelmy(1),ocdelmy(1),
     &                               condvolt(1),ecvoltmy(1),ocvoltmy(1),
     &                               condnumb(1),ecnumbmy(1),ocnumbmy(1),
     &                               iy,ecdelmy(ic),ecvoltmy(ic),ecnumbmy(ic))
              call subgrid_scan_axisxy(2,+1,cds,ny+1,maxcdel,icond,nx,ny,
     &                               iii(0,0),ecdelpy(1),ocdelpy(1),
     &                               condvolt(1),ecvoltpy(1),ocvoltpy(1),
     &                               condnumb(1),ecnumbpy(1),ocnumbpy(1),
     &                               iy,ecdelpy(ic),ecvoltpy(ic),ecnumbpy(ic))

c             --- If the point is near a conductor, save it.
              if (ecdelmx(ic)/levelxy < 1. .or.
     &            ecdelpx(ic)/levelxy < 1. .or.
     &            ecdelmy(ic)/levelxy < 1. .or.
     &            ecdelpy(ic)/levelxy < 1.) then
                iecndx(ic) = ix
                iecndy(ic) = iy
                iii(ix,iy) = ic
                ic = ic + 1

c               --- Increase size of conductor arrays if needed.
c               --- Increase in size is only an estimate.
                if (ic >= ncndmax) then
                  ncndmax = 1.1*max(ic,nocndbdy) + 1
                  call gchange("Conductor3d",0)
                endif

              endif

            endif
          enddo
        enddo

      enddo

c     --- Reset necndbdy
      necndbdy = ic - 1
 
      return
      end
c=============================================================================
      subroutine subgrid_scan_axisxy(ia,sgn,cds,limit,maxcdel,icond,nx,ny,iii,
     &                             ecdel,ocdel,condvolt,ecvolt,ocvolt,
     &                             condnumb,ecnumb,ocnumb,
     &                             cnd,cdel,volt,numb)
      integer(ISZ):: ia,sgn,cds(2),limit,maxcdel,nx,ny,icond
      integer(ISZ):: iii(0:nx,0:ny)
      real(kind=8):: ecdel(*),ocdel(*),condvolt(*),ecvolt(*),ocvolt(*)
      integer(ISZ):: condnumb(*),ecnumb(*),ocnumb(*)
      real(kind=8):: cdel,volt
      integer(ISZ):: cnd,numb

c Scan in the specified direction looking for a conductor or grid edge.
c If a conductor is found, save the distance to and the voltage on the
c conductor.

      integer(ISZ):: ic(2),ival
      ic = cds

      do while (ic(ia) .ne. limit)
        ival = iii(ic(1),ic(2))

        if (ival >= icond) then
c         --- A point inside a conductor was found. That means that the
c         --- conductor is aligned with the grid.
          cdel = sgn*(ic(ia) - cnd)
          volt = condvolt(ival-icond)
          numb = condnumb(ival-icond)
          return
        elseif (ival > 0) then
c         --- An even subgrid point was found. Check if it is near a
c         --- conductor.
          if (ecdel(ival) < maxcdel) then
c           --- If it is near a conductor, save the info.
            cdel = ecdel(ival) + sgn*(ic(ia) - cnd)
            volt = ecvolt(ival)
            numb = ecnumb(ival)
            return
          endif
        elseif (ival < 0) then
c         --- An odd subgrid point was found. Check if it is near a
c         --- conductor.
          if (ocdel(-ival) < maxcdel) then
c           --- If it is near a conductor, save the info.
            cdel = ocdel(-ival) + sgn*(ic(ia) - cnd)
            volt = ocvolt(-ival)
            numb = ocnumb(-ival)
            return
          endif
        endif
c       --- Otherwise, keep looking.
        ic(ia) = ic(ia) + sgn
      enddo

      return
      end
c=============================================================================
c=============================================================================
