#include "top.h"
c=============================================================================
c@(#) File frz_ImplicitES.F, version $Revision: 1.10 $, $Date: 2007/04/26 18:12:03 $
c# Copyright (c) 2007-2007, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the Poisson solver for the electrostatic implicit scheme.
c  David P. Grote, LLNL, (925)423-7194 or (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      module coefficientsmodulees2d
      contains
      function getbfieldsongrides2d(nx,nz,nzfull,dx,dz,xmmin,zmmin,bounds) result(bongrid)
      real(kind=8),pointer:: bongrid(:,:,:)
      integer(ISZ):: nx,nz,nzfull
      real(kind=8):: dx,dz,xmmin,zmmin
      integer(ISZ):: bounds(0:5)

      real(kind=8),allocatable:: xx(:),yy(:),zz(:),uz(:),gaminv(:)
      real(kind=8),allocatable:: ex(:),ey(:),ez(:),bendres(:),bendradi(:)
      integer(ISZ)::ix,iz

c     --- Allocate space for B field
      allocate(bongrid(-1:nx+1,-3:nz+3,0:2))

c     --- Create temporary space to hold the coordinates of the grid
c     --- needed to look up the B fields.
      allocate(xx(0:nx),yy(0:nx),zz(0:nx),uz(0:nx),gaminv(0:nx))
      allocate(ex(0:nx),ey(0:nx),ez(0:nx))
      allocate(bendres(0:nx),bendradi(0:nx))

c     --- Fill the transverse position arrays.
      do ix=0,nx
        xx(ix) = xmmin + ix*dx
        yy(ix) = 0.
      enddo
      uz = -1.
      gaminv = 1.
      bendres = 0.
      bendradi = LARGEPOS
      ex = 0.
      ey = 0.
      ez = 0.
      bongrid = 0.

c     --- Get the B field, one iz line at a time. One iz is done at a
c     --- time to save memory.
      do iz=0,nz
        zz = zmmin + iz*dz
        call exteb3d(1+nx,xx,yy,zz,uz,gaminv,0.,0.,
     &               bongrid(0:nx,iz,0),bongrid(0:nx,iz,1),bongrid(0:nx,iz,2),
     &               ex,ey,ez,1.,1.,bendres,bendradi,1.,1.)
      enddo

c     --- Apply the boundary conditions, setting the values in the guard cells.
      call applyboundaryconditionses2d(bongrid,nx,0,nz,nzfull,1,0,3,3,bounds)

      deallocate(xx,yy,zz,uz,gaminv)
      deallocate(ex,ey,ez)
      deallocate(bendres,bendradi)

      return
      end function getbfieldsongrides2d
c=============================================================================
      function getcoefficientses2d(nx,nz,dx,dz,xmmin,ns,qomdt,chi0,bongrid,
     &                             lcndbndy,icndbndy,conductors,lrz)
     &result(coeffs)
      use Constant
      use ConductorTypemodule
      real(kind=8),pointer:: coeffs(:,:,:,:)
      integer(ISZ):: nx,nz,ns
      real(kind=8):: dx,dz,xmmin
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: chi0(-1:nx+1,-3:nz+3,0:ns-1)
      real(kind=8):: bongrid(-1:nx+1,-3:nz+3,0:2)
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lrz

      real(kind=8),pointer:: chi(:,:,:,:)
      integer(ISZ):: ix,iz,js
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: ox,oy,oz,oo
      real(kind=8):: rr,rrm,rrp

      allocate(coeffs(0:2,0:2,0:nx,0:nz))
      allocate(chi(0:1,0:1,-1:nx+1,-1:nz+1))

c     --- First, generate the chi tensor from chi0 and the B field.
      chi = 0.
      do js=0,ns-1
        do iz=-1,nz+1
          do ix=-1,nx+1

            ox = 0.5*qomdt(js)*bongrid(ix,iz,0)
            oy = 0.5*qomdt(js)*bongrid(ix,iz,1)
            oz = 0.5*qomdt(js)*bongrid(ix,iz,2)
            oo = 1./(1. + sqrt(ox**2 + oy**2 + oz**2))

            chi(0,0,ix,iz) = chi(0,0,ix,iz) + chi0(ix,iz,js)*oo*(1. + ox**2)
            chi(0,1,ix,iz) = chi(0,1,ix,iz) + chi0(ix,iz,js)*oo*(ox*oz + oy)
            chi(1,0,ix,iz) = chi(1,0,ix,iz) + chi0(ix,iz,js)*oo*(ox*oz - oy)
            chi(1,1,ix,iz) = chi(1,1,ix,iz) + chi0(ix,iz,js)*oo*(1. + oz**2)

          enddo
        enddo
      enddo

      if (lrz) then
c       --- RZ Axisymmetry

c       --- If xmmin == 0, then the case of r=0 is included and must be treated
c       --- specially.
        ixmin = 0
        if (xmmin == 0.) ixmin = 1

        do iz=0,nz
          if (xmmin == 0.) then
            ix = 0

c           --- The r=0 case
            coeffs(1,1,ix,iz) = -2.*(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
     &                          -2.*(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
     &                          -(1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
     &                          -(1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,1,ix,iz) = +2.*(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
            coeffs(0,1,ix,iz) = +2.*(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
            coeffs(1,2,ix,iz) = +(1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
            coeffs(1,0,ix,iz) = +(1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,2,ix,iz) = 0.
            coeffs(0,2,ix,iz) = 0.
            coeffs(2,0,ix,iz) = 0.
            coeffs(0,0,ix,iz) = 0.

          endif

c         --- The general case
          do ix=ixmin,nx
            rr = xmmin + ix*dx
            rrm = (rr - 0.5*dx)/rr
            rrp = (rr + 0.5*dx)/rr

            coeffs(1,1,ix,iz) = -rrp*(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
     &                          -rrm*(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
     &                          -    (1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
     &                          -    (1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,1,ix,iz) = +rrp*(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
            coeffs(0,1,ix,iz) = +rrm*(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
            coeffs(1,2,ix,iz) = +0.5*(rrp-rrm)*(chi(1,0,ix  ,iz  ))/(2.*dz*dx)
     &                          +    (1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
            coeffs(1,0,ix,iz) = -0.5*(rrp-rrm)*(chi(1,0,ix  ,iz  ))/(2.*dz*dx)
     &                          +    (1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,2,ix,iz) = +0.5*rrp*chi(1,0,ix+1,iz  )/(2.*dz*dx)
     &                          +0.5*    chi(0,1,ix  ,iz+1)/(2.*dx*dz)
            coeffs(0,2,ix,iz) = -0.5*rrm*chi(1,0,ix-1,iz  )/(2.*dz*dx)
     &                          -0.5*    chi(0,1,ix  ,iz+1)/(2.*dx*dz)
            coeffs(2,0,ix,iz) = -0.5*rrp*chi(1,0,ix+1,iz  )/(2.*dz*dx)
     &                          -0.5*    chi(0,1,ix  ,iz-1)/(2.*dx*dz)
            coeffs(0,0,ix,iz) = +0.5*rrm*chi(1,0,ix-1,iz  )/(2.*dz*dx)
     &                          +0.5*    chi(0,1,ix  ,iz-1)/(2.*dx*dz)

          enddo
        enddo

      else
c       --- XZ Planar

        do iz=0,nz
          do ix=0,nx

            coeffs(1,1,ix,iz) = -(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
     &                          -(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
     &                          -(1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
     &                          -(1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,1,ix,iz) = +(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
            coeffs(0,1,ix,iz) = +(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
            coeffs(1,2,ix,iz) = +(1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
            coeffs(1,0,ix,iz) = +(1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,2,ix,iz) = +0.5*chi(1,0,ix+1,iz  )/(2.*dz*dx)
     &                          +0.5*chi(0,1,ix  ,iz+1)/(2.*dx*dz)
            coeffs(0,2,ix,iz) = -0.5*chi(1,0,ix-1,iz  )/(2.*dz*dx)
     &                          -0.5*chi(0,1,ix  ,iz+1)/(2.*dx*dz)
            coeffs(2,0,ix,iz) = -0.5*chi(1,0,ix+1,iz  )/(2.*dz*dx)
     &                          -0.5*chi(0,1,ix  ,iz-1)/(2.*dx*dz)
            coeffs(0,0,ix,iz) = +0.5*chi(1,0,ix-1,iz  )/(2.*dz*dx)
     &                          +0.5*chi(0,1,ix  ,iz-1)/(2.*dx*dz)

          enddo
        enddo

      endif

      deallocate(chi)

      return
      end function getcoefficientses2d
      end module coefficientsmodulees2d
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine mgsolveimplicites2d(iwhich,nx,nz,nzfull,dx,dz,phi,rho,ns,qomdt,chi0,
     &                               bounds,xmmin,zmmin,zmminglobal,zbeam,zgrid,
     &                               mgparam,mgiters,mgmaxiters,
     &                               mgmaxlevels,mgerror,mgtol,mgverbose,
     &                               downpasses,uppasses,
     &                               lcndbndy,laddconductor,icndbndy,lbuildquads,
     &                               gridmode,conductors,lrz,
     &                               my_index,nslaves,izfsslave,nzfsslave)
      use Subtimersfrz
      use ConductorTypemodule
      use Constant
      use coefficientsmodulees2d
      integer(ISZ):: iwhich
      integer(ISZ):: nx,nz,nzfull,ns
      real(kind=8):: phi(0:nx,-1:nz+1)
      real(kind=8):: rho(0:nx,0:nz)
      real(kind=8):: chi0(0:nx,0:nz,0:ns-1)
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,zmmin,zmminglobal,zbeam,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: mgverbose,downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor,lbuildquads
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      logical(ISZ):: lrz
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Use the multigrid method for solving Poisson's equation on a RZ Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution.
c
c When the grid cells are rectangular, semi-coarsening is done until the
c grid cell dimensions are roughly equal. Roughly equal means that
c   2/3 dx < dz < 4/3 dx
c This keeps (max(dz,dx) - min(dz,dx))/dx < 1/3. Currently, it is still
c assumed that dx ~ dy and that semi-coarsening is not needed transversely.
c
c Note that chi0 is the (q/m/2*rho/eps0*dt**2) for each implicit species.

      integer(ISZ):: i,k,ix,iz
      real(kind=8),pointer:: bongrid(:,:,:)
      real(kind=8),pointer:: chi0temp(:,:,:)
      real(kind=8),allocatable:: phisave(:,:)
      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: lzoffset(0:nslaves-1),rzoffset(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lfrztimesubs) substarttime = wtime()

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

c     --- Fill the B field array
c     --- This will be passed to vcycle and coarsened appropriately.
      bongrid => getbfieldsongrides2d(nx,nz,nzfull,dx,dz,xmmin,zmmin,bounds)
c     --- Parallel B.C.s are set below.

c     --- Create copy of rho with extra guard cells, and fill it in.
      allocate(chi0temp(-1:nx+1,-3:nz+3,0:ns-1))
      chi0temp(0:nx,0:nz,:) = chi0
      call applyboundaryconditionses2d(chi0temp,nx,0,nz,nzfull,1,0,3,ns,bounds)
c     --- Parallel B.C.s are set below.

c     --- Determine the points that make up the conductor.  This takes extra
c     --- time and so should not be done if the grid is not moving in the lab
c     --- frame.  Set gridmode to 1 to avoid this call. The data is then
c     --- converted and expanded for the multigrid solver.
      if (gridmode == 0 .or. iwhich == -2) then
        conductors%interior%n = 0
        conductors%evensubgrid%n = 0
        conductors%oddsubgrid%n = 0
        if (lbuildquads) then
          call setcndtr(xmmin,0.,zmmin,zmminglobal,zbeam,zgrid,nx,0,nz,dx,dx,dz,
     &                  bounds(4),bounds(5),bounds(1),.false.,.false.)
        endif
        if (laddconductor) call execuser("calladdconductor")
      endif
      call checkconductors(nx,0,nz,nzfull,dx,dx,dz,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave)


!$OMP PARALLEL
!$OMP&PRIVATE(i,ix,iz)

c     --- Make sure guard planes have sensible values before beginning.
      if (bounds(4) == 1) phi(:,-1) = phi(:,1)
      if (bounds(5) == 1) phi(:,nz+1) = phi(:,nz-1)
      if (bounds(4) == 2 .and. nz == nzfull) phi(:,-1) = phi(:,nz-1)
      if (bounds(5) == 2 .and. nz == nzfull) phi(:,nz:nz+1) = phi(:,0:1)
#ifdef MPIPARALLEL
      lzoffset = 0
      rzoffset = 0
      if (izfsslave(my_index) > 0) localbounds(4) = -1
      if (izfsslave(my_index)+nz < nzfull) localbounds(5) = -1
      call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                      bounds,nzfull,
     &                      lzoffset,rzoffset,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
c     --- These calls break the parallel field solver
c     call mgexchange_phi(nx,0,nz,nzfull,phi,localbounds,-1,0,0,1,
c    &                    my_index,nslaves,izfsslave,nzfsslave,
c    &                    whosendingleft,izsendingleft,
c    &                    whosendingright,izsendingright)
c     call mgexchange_phiperiodic(nx,0,nz,nzfull,phi,localbounds,1,0,0,1,
c    &                            my_index,nslaves,izfsslave,
c    &                            whosendingleft,whosendingright)
      call applyparallelboundaryconditionses2d(bongrid,nx,0,nz,nzfull,
     &                                         1,0,3,3,1,bounds,
     &                                         my_index,nslaves,
     &                                         izfsslave,nzfsslave,
     &                                         whosendingleft,izsendingleft,
     &                                         whosendingright,izsendingright)
      call applyparallelboundaryconditionses2d(chi0temp,nx,0,nz,nzfull,
     &                                         1,0,3,ns,1,bounds,
     &                                         my_index,nslaves,
     &                                         izfsslave,nzfsslave,
     &                                         whosendingleft,izsendingleft,
     &                                         whosendingright,izsendingright)
#endif

      allocate(phisave(0:nx,-1:nz+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,nz
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

c       --- Do one vcycle.
        call vcyclees2d(0,nx,nz,nzfull,dx,dz,xmmin,phi,rho,ns,qomdt,chi0temp,bongrid,
     &                  bounds,mgparam,mgmaxlevels,
     &                  downpasses,uppasses,lcndbndy,icndbndy,conductors,lrz,
     &                  my_index,nslaves,izfsslave,nzfsslave)

c       --- Calculate the change in phi.
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
c       print*,mgiters,mgerror

c       --- This line below seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop above.
c       mgerror = maxval(abs(phisave - phi))

      enddo

c     --- For Dirichlet boundary conditions, copy data into guard planes
c     --- For other boundary conditions, the guard planes are used during
c     --- the solve are so are already set.
      if (bounds(4) == 0) phi(:,-1) = 2.*phi(:,0) - phi(:,1)
      if (bounds(5) == 0) phi(:,nz+1) = 2.*phi(:,nz) - phi(:,nz-1)

c     --- Make a print out.
      if (mgerror > mgtol) then
        call remark("ImplicitMultigridES2d: Maximum number of iterations reached")
      endif
      if (mgverbose >= 1) then
        write(errline,20) mgerror,mgiters
  20    format("ImplicitMultigridES2d: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif

      deallocate(phisave)
      deallocate(chi0temp)
      deallocate(bongrid)

!$OMP END PARALLEL

      if (lfrztimesubs) timemgsolveimplicites2d = timemgsolveimplicites2d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcyclees2d(mglevel,nx,nz,nzfull,dx,dz,xmmin,
     &                                phi,rho,ns,qomdt,chi0,bongrid,globalbounds,
     &                                mgparam,
     &                                mgmaxlevels,downpasses,uppasses,
     &                                lcndbndy,icndbndy,conductors,lrz,
     &                                my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use coefficientsmodulees2d
      use formggetarraysuminterface
      integer(ISZ):: mglevel
      integer(ISZ):: nx,nz,nzfull,ns
      real(kind=8):: dx,dz,xmmin
      real(kind=8):: phi(0:nx,-1:nz+1),rho(0:nx,0:nz)
      real(kind=8):: chi0(-1:nx+1,-3:nz+3,0:ns-1),bongrid(-1:nx+1,-3:nz+3,0:2)
      real(kind=8):: qomdt(0:ns-1)
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lrz
      integer(ISZ):: icndbndy
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: mingridsize
      real(kind=8),pointer:: coeffs(:,:,:,:)
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: chi0coarse(:,:,:),bongridcoarse(:,:,:)
      real(kind=8),allocatable:: res(:,:)
      integer(ISZ):: i,ic,iszone=1
      integer(ISZ):: nxcoarse,nzcoarse,nzfullcoarse
      real(kind=8):: dxcoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dzcoarsesqi
      integer(ISZ):: localbounds(0:5),localboundsc(0:5)
      integer(ISZ):: lzoffsetall(0:nslaves-1),rzoffsetall(0:nslaves-1)
      integer(ISZ):: lzoffset,rzoffset
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)
      integer(ISZ):: whosendingleftc(0:nslaves-1), izsendingleftc(0:nslaves-1)
      integer(ISZ):: whosendingrightc(0:nslaves-1),izsendingrightc(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      integer(ISZ):: allocerror
      real(kind=8):: sss(2)
      logical(ISZ):: larraysumdiagnostic=.false.

c     --- Get the finite difference coefficients for the Poisson operator at
c     --- the current level of refinement.
      coeffs => getcoefficientses2d(nx,nz,dx,dz,xmmin,ns,qomdt,chi0,bongrid,
     &                              lcndbndy,icndbndy,conductors,lrz)

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

      if (larraysumdiagnostic) then
        sss = mggetarraysum(nx,0,nz,0,0,1,phi,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V1 phi",mglevel,sss
        sss = mggetarraysum(nx,0,nz,0,0,0,rho,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V1 rho",mglevel,sss
        sss = mggetarraysum(nx,0,nz,1,0,3,chi0,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V1 chi0",mglevel,sss
      endif

c     --- Do initial SOR passes.
      do i=1,downpasses
        call relaxes2d(nx,nz,nzfull,phi,rho,coeffs,localbounds,mgparam,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

      if (larraysumdiagnostic) then
        sss = mggetarraysum(nx,0,nz,0,0,1,phi,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nx,0,nz,0,0,0,rho,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V2 rho",mglevel,sss
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. nzfull >= 4 .and. mglevel < mgmaxlevels) then

        allocate(res(-1:nx+1,-3:nz+3),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nx,nz
          call kaboom("vcyclees2d: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residuales2d(nx,nz,nzfull,phi,rho,coeffs,res,
     &                    mglevel,localbounds,mgparam)
#ifdef MPIPARALLEL
        call applyparallelboundaryconditionses2d(res,nx,0,nz,nzfull,1,0,3,1,3,
     &                                           localbounds,
     &                                           my_index,nslaves,
     &                                           izfsslave,nzfsslave,
     &                                           whosendingleft,izsendingleft,
     &                                           whosendingright,izsendingright)
#endif

        if (larraysumdiagnostic) then
          sss = mggetarraysum(nx,0,nz,1,0,3,res,
     &                        izfsslave,nzfsslave,nslaves,my_index)
          if (my_index == 0) print*,"V3 res",mglevel,sss
        endif

c       --- Calculate the size of the next coarsest grid. If the current
c       --- size is twice and odd integer, the next coarsest size will be
c       --- n/2+1, gauranteeing that it is even.
c       --- Note that nzcoarse is only used if full-coarsening is done.
        nxcoarse = nx/2
        if (mod(nxcoarse,2) == 1) nxcoarse = nxcoarse + 1
        nzcoarse = nz/2
        if (mod(nzcoarse,2) == 1) nzcoarse = nzcoarse + 1
        nzfullcoarse = nzfull/2
        if (mod(nzfullcoarse,2) == 1) nzfullcoarse = nzfullcoarse + 1

        mingridsize = min(dx,dy)
        mingridsize = min(mingridsize,dz)

        if (dx > 4./3.*mingridsize) then
          nxcoarse = nx
        endif
        if (dz > 4./3.*mingridsize) then
          nzcoarse = nz
          nzfullcoarse = nzfull
        endif

        dxcoarse = dx*nx/nxcoarse
        dzcoarse = dz*nzfull/nzfullcoarse
        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2

        localboundsc = globalbounds

#ifdef MPIPARALLEL
c       --- Find domains in coarser grid
        call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                  nzfull,nzfullcoarse)
c       --- Reset value to corrected one
        nzcoarse = nzfsslavec(my_index)
c       --- Difference between starts and ends of coarse and fine grids.
c       --- Should only be in the range 0-2.
        lzoffsetall = (nzfullcoarse*izfsslave-nzfull*izfsslavec)
        rzoffsetall = (nzfull*(izfsslavec + nzfsslavec) -
     &                 nzfullcoarse*(izfsslave + nzfsslave))
c       --- Note that the lzoffsetall and rzoffsetall can only be used in
c       --- MPIPARALLEL sections since they will be unallocated in the
c       --- serial code. So, separate scalars are used in code which is
c       --- used in the serial version.
        lzoffset = lzoffsetall(my_index)
        rzoffset = rzoffsetall(my_index)
c       --- Get processor with which to exchange data on coarse grid
        call mggetexchangepes(nslaves,izfsslavec,nzfsslavec,my_index,
     &                        globalbounds,nzfullcoarse,
     &                        lzoffsetall,rzoffsetall,
     &                        whosendingleftc,izsendingleftc,
     &                        whosendingrightc,izsendingrightc)
        if (izfsslavec(my_index) > 0) localboundsc(4) = -1
        if (izfsslavec(my_index) + nzcoarse < nzfullcoarse) localboundsc(5) = -1
#else
        lzoffset = 0
        rzoffset = 0
#endif

c       --- Alloate new work space
        allocate(phicoarse(0:nxcoarse,-1:nzcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxcoarse,nzcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxcoarse,0:nzcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxcoarse,nzcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(chi0coarse(-1:nxcoarse+1,-3:nzcoarse+3,0:ns-1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate chi0coarse to shape ",
     &           nxcoarse,nzcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(bongridcoarse(-1:nxcoarse+1,-3:nzcoarse+3,0:2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate bongridcoarse to shape ",
     &           nxcoarse,nzcoarse,3
          call kaboom("vcyclees2d: allocation error")
          return
        endif

        rhocoarse = 0.
        phicoarse = 0.
        chi0coarse = 0.
        bongridcoarse = 0.

c       --- Restriction
        call restrictes3d(nx,0,nz,nzfull,res,1,0,3,
     &                    nxcoarse,0,nzcoarse,nzfullcoarse,rhocoarse,0,0,0,
     &                    localbounds,localboundsc,lzoffset)
        if (larraysumdiagnostic) then
          sss = mggetarraysum(nxcoarse,0,nzcoarse,0,0,0,rhocoarse,
     &                        izfsslavec,nzfsslavec,nslaves,my_index)
          if (my_index == 0) print*,"V3 rhocoarse",mglevel,sss
          sss = mggetarraysum(nxcoarse,0,nzcoarse,1,0,3,chi0coarse,
     &                        izfsslavec,nzfsslavec,nslaves,my_index)
          if (my_index == 0) print*,"V3 chi0coarse",mglevel,sss
        endif

        do ic=0,ns-1
          call restrictes3d(nx,0,nz,nzfull,chi0(:,:,ic),1,0,3,
     &                      nxcoarse,0,nzcoarse,nzfullcoarse,chi0coarse(:,:,ic),1,0,3,
     &                      localbounds,localboundsc,lzoffset)
        enddo
        do ic=0,2
          call restrictes3d(nx,0,nz,nzfull,bongrid(:,:,ic),1,0,3,
     &                      nxcoarse,0,nzcoarse,nzfullcoarse,bongridcoarse(:,:,ic),1,0,3,
     &                      localbounds,localboundsc,lzoffset)
        enddo
        call applyboundaryconditionses2d(chi0coarse,nxcoarse,0,nzcoarse,nzfullcoarse,1,0,3,ns,localboundsc)
        call applyboundaryconditionses2d(bongridcoarse,nxcoarse,0,nzcoarse,nzfullcoarse,1,0,3,3,localboundsc)
        call applyparallelboundaryconditionses2d(chi0coarse,nxcoarse,0,nzcoarse,nzfullcoarse,
     &                                           1,0,3,ns,1,localboundsc,
     &                                           my_index,nslaves,
     &                                           izfsslavec,nzfsslavec,
     &                                           whosendingleftc,izsendingleftc,
     &                                           whosendingrightc,izsendingrightc)
        call applyparallelboundaryconditionses2d(bongridcoarse,nxcoarse,0,nzcoarse,nzfullcoarse,
     &                                           1,0,3,1,1,localboundsc,
     &                                           my_index,nslaves,
     &                                           izfsslavec,nzfsslavec,
     &                                           whosendingleftc,izsendingleftc,
     &                                           whosendingrightc,izsendingrightc)

c       --- Continue at the next coarsest level.
        call vcyclees2d(mglevel+iszone,nxcoarse,nzcoarse,nzfullcoarse,
     &                  dxcoarse,dzcoarse,xmmin,
     &                  phicoarse,rhocoarse,ns,qomdt,chi0coarse,bongridcoarse,
     &                  globalbounds,mgparam,
     &                  mgmaxlevels,downpasses,uppasses,
     &                  lcndbndy,icndbndy,conductors,lrz,
     &                  my_index,nslaves,izfsslavec,nzfsslavec)

        if (larraysumdiagnostic) then
          sss = mggetarraysum(nxcoarse,0,nzcoarse,0,0,1,phicoarse,
     &                        izfsslavec,nzfsslavec,nslaves,my_index)
          if (my_index == 0) print*,"V4 phicoarse",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand3d(nx,0,nz,nzfull,phi,
     &                nxcoarse,0,nzcoarse,nzfullcoarse,phicoarse,
     &                localbounds,lzoffset)
#ifdef MPIPARALLEL
        call mgexchange_phiperiodic(nx,0,nz,nzfull,phi,localbounds,1,0,0,1,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)
#endif
        if (larraysumdiagnostic) then
          sss = mggetarraysum(nx,0,nz,0,0,1,phi,
     &                        izfsslave,nzfsslave,nslaves,my_index)
          if (my_index == 0) print*,"V5 phi",mglevel,sss
        endif
        
        deallocate(res)
        deallocate(phicoarse,rhocoarse)
        deallocate(chi0coarse,bongridcoarse)

      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call relaxes2d(nx,nz,nzfull,phi,rho,coeffs,localbounds,mgparam,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

      deallocate(coeffs)

      return
      end
c=============================================================================
      subroutine relaxes2d(nx,nz,nzfull,phi,rho,coeffs,bounds,mgparam,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,-1:nz+1),rho(0:nx,0:nz),coeffs(0:2,0:2,0:nx,0:nz)
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.

      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ixmin,ixmax,izmin,izmax,ix,iz,ix1
      integer(ISZ):: ixm1,ixp1,izm1,izp1

c     --- Set starting and ending parity.
#ifdef MPIPARALLEL
      s_parity = mod(izfsslave(my_index),2)
      e_parity = mod(s_parity+1,2)
#else
      s_parity = 0
      e_parity = 1
#endif

c     --- Set min and max indices for phi array.
c     --- If using Dirichlet boundary conditions, do not solve for the
c     --- potential on the grid edge.
      ixmin = 1
      ixmax = nx-1
      izmin = 0
      izmax = nz
      if (bounds(0) > 0) ixmin = 0
      if (bounds(1) > 0) ixmax = nx
      if (bounds(4) < 1) izmin = 1
      if (bounds(5) < 1) izmax = nz - 1

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity

        do iz=izmin,izmax
          izm1 = iz - 1
          izp1 = iz + 1
          ix1 = ixmin + mod(ixmin + iz + parity,2)
          do ix=ix1,ixmax,2
            ixm1 = ix - 1
            if (ix==0 .and. bounds(0) == 1) ixm1 = 1
            if (ix==0 .and. bounds(0) == 2) ixm1 = nx-1
            ixp1 = ix + 1
            if (ix == nx .and. bounds(1) == 1)  ixp1 = nx-1
            if (ix == nx .and. bounds(1) == 2)  ixp1 = 1

            phi(ix,iz) = mgparam*(rho(ix,iz)/eps0 +
     &                            coeffs(0,0,ix,iz)*phi(ixm1,izm1) +
     &                            coeffs(1,0,ix,iz)*phi(ix  ,izm1) +
     &                            coeffs(2,0,ix,iz)*phi(ixp1,izm1) +
     &                            coeffs(0,1,ix,iz)*phi(ixm1,iz  ) +
     &                            coeffs(2,1,ix,iz)*phi(ixp1,iz  ) +
     &                            coeffs(0,2,ix,iz)*phi(ixm1,izp1) +
     &                            coeffs(1,2,ix,iz)*phi(ix  ,izp1) +
     &                            coeffs(2,2,ix,iz)*phi(ixp1,izp1))/
     &                           (-coeffs(1,1,ix,iz)) +
     &                   (1.-mgparam)*phi(ix,iz)

          enddo
        enddo

c       --- set phi in the z guard planes
        if (bounds(4) == 1) phi(:,-1) = phi(:,+1)
        if (bounds(5) == 1) phi(:,nz+1) = phi(:,nz-1)
        if (bounds(4) == 2 .and. nz == nzfull) phi(:,-1) = phi(:,nz-1)
        if (bounds(5) == 2 .and. nz == nzfull) phi(:,nz:nz+1) = phi(:,0:1)
#ifdef MPIPARALLEL
        call mgexchange_phi(nx,0,nz,nzfull,phi,bounds,0,0,0,1,
     &                      my_index,nslaves,izfsslave,nzfsslave,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
        call mgexchange_phiperiodic(nx,0,nz,nzfull,phi,bounds,1,0,0,1,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- Exchange phi in the z guard planes
      call mgexchange_phi(nx,0,nz,nzfull,phi,bounds,-1,0,0,1,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
#endif

      return
      end
c=============================================================================
      subroutine residuales2d(nx,nz,nzfull,phi,rho,coeffs,res,
     &                        mglevel,bounds,mgparam)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: phi(0:nx,-1:nz+1),rho(0:nx,0:nz)
      real(kind=8):: res(-1:nx+1,-3:nz+3),coeffs(0:2,0:2,0:nx,0:nz)
      integer(ISZ):: mglevel,bounds(0:5)
      real(kind=8):: mgparam

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      integer(ISZ):: ixm1,ixp1,izm1,izp1

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 1
      ixmax = nx-1
      izmin = 0
      izmax = nz
      if (bounds(0) >  0) ixmin = 0
      if (bounds(1) == 1) ixmax = nx
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nz-1

      res = 0.
c     --- Calculate the residual.
!$OMP DO
      do iz=izmin,izmax
        izm1 = iz - 1
        izp1 = iz + 1

        do ix=ixmin,ixmax
          ixm1 = ix - 1
          if (ix==0 .and. bounds(0) == 1) ixm1 = 1
          if (ix==0 .and. bounds(0) == 2) ixm1 = nx-1
          ixp1 = ix + 1
          if (ix == nx .and. bounds(1) == 1)  ixp1 = nx-1
          if (ix == nx .and. bounds(1) == 2)  ixp1 = 1

          res(ix,iz) = rho(ix,iz) + eps0*(
     &                  coeffs(0,0,ix,iz)*phi(ixm1,izm1) +
     &                  coeffs(1,0,ix,iz)*phi(ix  ,izm1) +
     &                  coeffs(2,0,ix,iz)*phi(ixp1,izm1) +
     &                  coeffs(0,1,ix,iz)*phi(ixm1,iz  ) +
     &                  coeffs(1,1,ix,iz)*phi(ix  ,iz  ) +
     &                  coeffs(2,1,ix,iz)*phi(ixp1,iz  ) +
     &                  coeffs(0,2,ix,iz)*phi(ixm1,izp1) +
     &                  coeffs(1,2,ix,iz)*phi(ix  ,izp1) +
     &                  coeffs(2,2,ix,iz)*phi(ixp1,izp1))

        enddo
      enddo
!$OMP END DO

c     --- Transverse boundaries
      if (bounds(0) == 0) res(-1:0,:) = 0.
      if (bounds(1) == 0) res(nx:nx+1,:) = 0.
      if (bounds(0) == 1) res(-1,:) = res(1,:)
      if (bounds(1) == 1) res(nx+1,:) = res(nx-1,:)
      if (bounds(0) == 2) res(-1,:) = res(nx-1,:)
      if (bounds(1) == 2) res(nx:nx+1,:) = res(0:1,:)

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
      subroutine restrictes3d(nx,ny,nz,nzfull,u,delx,dely,delz,
     &                        nxcoarse,nycoarse,nzcoarse,nzfullcoarse,ucoarse,
     &                        delcx,delcy,delcz,
     &                        bounds,boundscoarse,lzoffset)

      integer(ISZ):: nx,ny,nz,nzfull,delx,dely,delz,delcx,delcy,delcz
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nzfullcoarse
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz)
      real(kind=8):: ucoarse(-delcx:nxcoarse+delcx,-delcy:nycoarse+delcy,-delcz:nzcoarse+delcz)
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
      dxi = 1.*nxcoarse/nx
      if (ny > 0) then
        dy = 1.*ny/nycoarse
        dyi = 1.*nycoarse/ny
      else
        dy = 1.
        dyi = 1.
      endif
      dz = 1.*nzfull/nzfullcoarse
      dzi = 1.*nzfullcoarse/nzfull

c     --- Do the loops.
!$OMP DO
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = ((izcoarse-1)*nzfull - lzoffset + 4*nzfullcoarse)/nzfullcoarse-3
        izmax = ((izcoarse+1)*nzfull - lzoffset - 1)/nzfullcoarse
        if (izmin < -delz) izmin = -delz
        if (izmax > nz+delz) izmax = nz+delz

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzfullcoarse)*dzi)
        enddo

        do iycoarse=iycoarsemin,iycoarsemax
          iymin = int(ceiling((iycoarse-1)*dy + 1.e-10))
          iymax = int(floor((iycoarse+1)*dy - 1.e-10))
          if (iymin < -dely) iymin = -dely
          if (iymax > ny+dely) iymax = ny+dely

          do iy=iymin,iymax
            wy(iy-iymin) = 1. - abs(iycoarse - iy*dyi)
          enddo

          do ixcoarse=ixcoarsemin,ixcoarsemax
            ixmin = int(ceiling((ixcoarse-1)*dx + 1.e-10))
            ixmax = int(floor((ixcoarse+1)*dx - 1.e-10))
            if (ixmin < -delx) ixmin = -delx
            if (ixmax > nx+delx) ixmax = nx+delx

            do ix=ixmin,ixmax
              wx(ix-ixmin) = 1. - abs(ixcoarse - ix*dxi)
            enddo

            r = 0.
            w = 0.
            do iz=izmin,izmax
              do iy=iymin,iymax
                do ix=ixmin,ixmax
                  r = r + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)*u(ix,iy,iz)
                  w = w + wx(ix-ixmin)*wy(iy-iymin)*wz(iz-izmin)
                enddo
              enddo
            enddo
            if (w > 0.) then
              ucoarse(ixcoarse,iycoarse,izcoarse) = r/w
            else
              ucoarse(ixcoarse,iycoarse,izcoarse) = 0.
            endif

          enddo
        enddo
      enddo
!$OMP END DO

      return
      end
c=============================================================================
c=============================================================================
      subroutine applyboundaryconditionses2d(u,nx,ny,nz,nzfull,
     &                                       delx,dely,delz,ncomp,bounds)
      integer(ISZ):: nx,ny,nz,nzfull,delx,dely,delz,ncomp
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz,ncomp)
      integer(ISZ):: bounds(0:5)

c     --- Apply the boundary conditions, setting the values in the guard cells
      if (delx > 0) then
        if (bounds(0) == 0) u(-1,0:ny,0:nz,:) = 0.
        if (bounds(0) == 1) u(-1,0:ny,0:nz,:) = u(+1,0:ny,0:nz,:)
        if (bounds(0) == 2) u(-1,0:ny,0:nz,:) = u(nx-1,0:ny,0:nz,:)
        if (bounds(1) == 0) u(nx+1,0:ny,0:nz,:) = 0.
        if (bounds(1) == 1) u(nx+1,0:ny,0:nz,:) = u(nx-1,0:ny,0:nz,:)
        if (bounds(1) == 2) u(nx+1,0:ny,0:nz,:) = u(+1,0:ny,0:nz,:)
      endif

      if (dely > 0) then
        if (bounds(2) == 0) u(:,-1,0:nz,:) = 0.
        if (bounds(2) == 1) u(:,-1,0:nz,:) = u(:,+1,0:nz,:)
        if (bounds(2) == 2) u(:,-1,0:nz,:) = u(:,ny-1,0:nz,:)
        if (bounds(3) == 0) u(:,ny+1,0:nz,:) = 0.
        if (bounds(3) == 1) u(:,ny+1,0:nz,:) = u(:,ny-1,0:nz,:)
        if (bounds(3) == 2) u(:,ny+1,0:nz,:) = u(:,+1,0:nz,:)
      endif

      if (delz > 0) then
        if (bounds(4) == 0) u(:,:,-delz:-1,:) = 0.
        if (bounds(4) == 1) u(:,:,-delz:-1,:) = u(:,:,delz:+1:-1,:)
        if (bounds(4) == 2 .and. nz == nzfull) u(:,:,-delz:-1,:) = u(:,:,nz-delz:nz-1,:)
        if (bounds(5) == 0) u(:,:,nz+1:nz+delz,:) = 0.
        if (bounds(5) == 1) u(:,:,nz+1:nz+delz,:) = u(:,:,nz-1:nz-delz:-1,:)
        if (bounds(5) == 2 .and. nz == nzfull) u(:,:,nz+1:nz+delz,:) = u(:,:,+1:delz,:)
      endif

c     --- Need code to set z guard planes for periodic B.C.'s for parallel

      return
      end
c=============================================================================
      subroutine applyparallelboundaryconditionses2d(u,nx,ny,nz,nzfull,
     &                                               delx,dely,delz,ncomp,zsend,bounds,
     &                                               my_index,nslaves,
     &                                               izfsslave,nzfsslave,
     &                                               whosendingleft,izsendingleft,
     &                                               whosendingright,izsendingright)
      integer(ISZ):: nx,ny,nz,nzfull,delx,dely,delz,ncomp,zsend
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nz+delz,ncomp)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

      integer(ISZ):: iz,ic

#ifdef MPIPARALLEL

      do ic=1,ncomp
        do iz=1,zsend
          call mgexchange_phi(nx,ny,nz,nzfull,u(:,:,:,ic),bounds,-iz,delx,dely,delz,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
        enddo

        call mgexchange_phiperiodic(nx,ny,nz,nzfull,u(:,:,:,ic),bounds,0,delx,dely,delz,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)

      enddo

#endif

      return
      end
c=============================================================================
