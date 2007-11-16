#include "top.h"
c=============================================================================
c@(#) File frz_ImplicitES.F, version $Revision: 1.22 $, $Date: 2007/11/16 00:54:14 $
c# Copyright (c) 2007-2007, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the Poisson solver for the electrostatic implicit scheme.
c  David P. Grote, LLNL, (925)423-7194 or (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      module coefficientsmodulees2d
      contains
      function getbfieldsongrides2d(nx,nzlocal,nz,dx,dz,xmmin,zmminlocal,bounds) result(bongrid)
      real(kind=8),pointer:: bongrid(:,:,:)
      integer(ISZ):: nx,nzlocal,nz
      real(kind=8):: dx,dz,xmmin,zmminlocal
      integer(ISZ):: bounds(0:5)

      real(kind=8),allocatable:: xx(:),yy(:),zz(:),uz(:),gaminv(:)
      real(kind=8),allocatable:: ex(:),ey(:),ez(:),bendres(:),bendradi(:)
      integer(ISZ)::ix,iz

c     --- Allocate space for B field
      allocate(bongrid(-1:nx+1,-3:nzlocal+3,0:2))

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
      do iz=0,nzlocal
        zz = zmminlocal + iz*dz
        call exteb3d(1+nx,xx,yy,zz,uz,gaminv,0.,0.,
     &               bongrid(0:nx,iz,0),bongrid(0:nx,iz,1),bongrid(0:nx,iz,2),
     &               ex,ey,ez,1.,1.,bendres,bendradi,1.,1.)
      enddo

c     --- Apply the boundary conditions, setting the values in the guard cells.
      call applyboundaryconditionses2d(bongrid,nx,0,nzlocal,nz,1,0,3,3,bounds,.true.)

      deallocate(xx,yy,zz,uz,gaminv)
      deallocate(ex,ey,ez)
      deallocate(bendres,bendradi)

      return
      end function getbfieldsongrides2d
c=============================================================================
      subroutine getcoefficientses2d(nx,nzlocal,dx,dz,xmmin,ns,qomdt,rho,chi0,
     &                               bongrid,
     &                               mglevel,lcndbndy,icndbndy,conductors,lrz,
     &                               bounds,coeffs,rhs)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nx,nzlocal,ns
      real(kind=8):: dx,dz,xmmin
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: rho(0:nx,0:nzlocal)
      real(kind=8):: chi0(-1:nx+1,-3:nzlocal+3,0:ns-1)
      real(kind=8):: bongrid(-1:nx+1,-3:nzlocal+3,0:2)
      integer(ISZ):: mglevel
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lrz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: coeffs(0:2,0:2,0:nx,0:nzlocal)
      real(kind=8):: rhs(0:nx,0:nzlocal)

      real(kind=8),pointer:: chi(:,:,:,:)
      integer(ISZ):: ix,iz,js,ic
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: ox,oy,oz,oo
      real(kind=8):: rr,rrm,rrp

      allocate(chi(0:1,0:1,-1:nx+1,-1:nzlocal+1))
      chi = 0.

c     --- First, generate the chi tensor from chi0 and the B field.
      do js=0,ns-1
        do iz=-1,nzlocal+1
          do ix=-1,nx+1

            ox = 0.5*qomdt(js)*bongrid(ix,iz,0)
            oy = 0.5*qomdt(js)*bongrid(ix,iz,1)
            oz = 0.5*qomdt(js)*bongrid(ix,iz,2)
            oo = 1./(1. + ox**2 + oy**2 + oz**2)

            chi(0,0,ix,iz) = chi(0,0,ix,iz) + chi0(ix,iz,js)*oo*(1. + ox**2)
            chi(0,1,ix,iz) = chi(0,1,ix,iz) + chi0(ix,iz,js)*oo*(ox*oz + oy)
            chi(1,0,ix,iz) = chi(1,0,ix,iz) + chi0(ix,iz,js)*oo*(ox*oz - oy)
            chi(1,1,ix,iz) = chi(1,1,ix,iz) + chi0(ix,iz,js)*oo*(1. + oz**2)

          enddo
        enddo

c       --- Force chi to zero inside of any conductors. This is because
c       --- the charge density is always zero there.
c       --- This is probably not needed.
        do ic = conductors%interior%istart(mglevel),
     &          conductors%interior%istart(mglevel+1)-1
          ix = conductors%interior%indx(0,ic)
          iz = conductors%interior%indx(2,ic)
          chi(:,:,ix,iz) = 0.
        enddo

      enddo

      if (lrz) then
c       --- RZ Axisymmetry

c       --- If xmmin == 0, then the case of r=0 is included and must be treated
c       --- specially.
        ixmin = 0
        if (xmmin == 0.) ixmin = 1

        do iz=0,nzlocal
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

        do iz=0,nzlocal
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

      rhs = 0.
      call getcoefficientssubgrides2d(conductors%evensubgrid,nx,nzlocal,dx,dz,
     &                                xmmin,mglevel,lcndbndy,icndbndy,lrz,chi,
     &                                bounds,coeffs,rhs)
      call getcoefficientssubgrides2d(conductors%oddsubgrid,nx,nzlocal,dx,dz,
     &                                xmmin,mglevel,lcndbndy,icndbndy,lrz,chi,
     &                                bounds,coeffs,rhs)

      deallocate(chi)

      return
      end subroutine getcoefficientses2d
c=============================================================================
      subroutine getcoefficientssubgrides2d(subgrid,nx,nzlocal,dx,dz,xmmin,
     &                                      mglevel,lcndbndy,icndbndy,lrz,chi,
     &                                      bounds,coeffs,rhs)
      use Constant
      use ConductorSubgridTypemodule
      type(ConductorSubgridType):: subgrid
      integer(ISZ):: nx,nzlocal
      real(kind=8):: dx,dz,xmmin
      integer(ISZ):: mglevel
      integer(ISZ):: icndbndy
      logical(ISZ):: lcndbndy,lrz
      real(kind=8):: chi(0:1,0:1,-1:nx+1,-1:nzlocal+1)
      integer(ISZ):: bounds(0:5)
      real(kind=8):: coeffs(0:2,0:2,0:nx,0:nzlocal)
      real(kind=8):: rhs(0:nx,0:nzlocal)

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.

      integer(ISZ):: ic,ix,iz
      real(kind=8):: dxm,dzm,dxp,dzp
      real(kind=8):: pxm,pzm,pxp,pzp
      real(kind=8):: rrm,rrp
      real(kind=8):: chixm(0:1,0:1),chizm(0:1,0:1)
      real(kind=8):: chixp(0:1,0:1),chizp(0:1,0:1)
      real(kind=8):: cxm,czm,cxp,czp
      real(kind=8):: dx1,dx2,dz1,dz2
      real(kind=8):: voltfac,c0
      real(kind=8):: rr(0:nx),rri(0:nx)
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      if (.not. lcndbndy) return

      if (icndbndy == 1) then
        call kaboom("mgsolveimplicites2d: icndbndy == 1 not supported")
        return
      endif

      dels => subgrid%dels
      volt => subgrid%volt

c     --- Only use actual voltage on finest level. Set to zero for
c     --- coarser levels since solver for the residuals.
      if (mglevel == 0) then
        voltfac = 1.
      else
        voltfac = 0.
      endif

c     --- Precalculate the radius for efficiency
      if (lrz) then
        do ix=0,nx
          rr(ix) = ix*dx + xmmin
          if (rr(ix) > 0.) rri(ix) = 1./rr(ix)
        enddo
      endif

c     --- Loop over points near surface of conductors.
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet or parallel boundary
        if (ix == 0  .and. bounds(0) < 1) cycle
        if (ix == nx .and. bounds(1) < 1) cycle
        if (iz == 0  .and. bounds(4) < 1) cycle
        if (iz == nzlocal .and. bounds(5) < 1) cycle

c       --- First, get the effective grid cell sizes and the chi's.
c       --- The chi's inside the conductor are always zero because the
c       --- charge density is zero there.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          dxm = dels(0,ic)
          chixm = 0.
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
          chixm = 0.
        else
          dxm = 1.
          chixm = chi(:,:,ix-1,iz)
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          dxp = dels(1,ic)
          chixp = 0.
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
          chixp = 0.
        else
          dxp = 1.
          chixp = chi(:,:,ix+1,iz)
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          dzm = dels(4,ic)
          chizm = 0.
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
          chizm = 0.
        else
          dzm = 1.
          chizm = chi(:,:,ix,iz-1)
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          dzp = dels(5,ic)
          chizp = 0.
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
          chizp = 0.
        else
          dzp = 1.
          chizp = chi(:,:,ix,iz+1)
        endif

c       --- Now construct the coefficients
        cxm = 1./(dxm*(0.5*dxm + 0.5*dxp))
        cxp = 1./(dxp*(0.5*dxm + 0.5*dxp))
        czm = 1./(dzm*(0.5*dzm + 0.5*dzp))
        czp = 1./(dzp*(0.5*dzm + 0.5*dzp))
        if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) cxm = 0.
        if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) cxp = 0.
        if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) czm = 0.
        if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) czp = 0.

c       --- Correct coefficients for axisymmetric case
        if (lrz) then
          if (rr(ix) > 0.) then
            rrm = (rr(ix) - 0.5*dxm*dx)*rri(ix)
            rrp = (rr(ix) + 0.5*dxm*dx)*rri(ix)
            cxm = cxm*rrm
            cxp = cxp*rrp
          else
            rrm = 0.
            rrp = 0.
            cxm = 0.
            cxp = 4.*cxp
          endif
        endif

c       --- The dx1 and dz1 are for finite differences (i+1/2) - (i-1/2)
c       --- the dx2 and dz2 are for finite differences (i+1) - (i-1)
        dx1 = dx*(0.5*dxm + 0.5*dxp)
        dz1 = dz*(0.5*dzm + 0.5*dzp)
        dx2 = dx*(dxm + dxp)
        dz2 = dz*(dzm + dzp)

        coeffs(2,1,ix,iz) = +cxp*(1. + 0.5*(chi(0,0,ix,iz)+chixp(0,0)))/dx**2
        coeffs(0,1,ix,iz) = +cxm*(1. + 0.5*(chixm(0,0)+chi(0,0,ix,iz)))/dx**2
        coeffs(1,2,ix,iz) = +czp*(1. + 0.5*(chi(1,1,ix,iz)+chizp(1,1)))/dz**2
        coeffs(1,0,ix,iz) = +czm*(1. + 0.5*(chizm(1,1)+chi(1,1,ix,iz)))/dz**2
        coeffs(1,1,ix,iz) = -coeffs(2,1,ix,iz)-coeffs(0,1,ix,iz)
     &                      -coeffs(1,2,ix,iz)-coeffs(1,0,ix,iz)
        if (lrz) then
          coeffs(1,2,ix,iz) = coeffs(1,2,ix,iz)
     &                        + 0.5*(rrp-rrm)*(chi(1,0,ix,iz))/(dz2*dx1)
          coeffs(1,0,ix,iz) = coeffs(1,0,ix,iz)
     &                        - 0.5*(rrp-rrm)*(chi(1,0,ix,iz))/(dz2*dx1)
        endif

c       --- Now the right hand side and be generated and the associated
c       --- coefficients zeroed out.
c       --- Note that for the Neumann case, the potential is not
c       --- used and so there are not contributions to the rhs.
c       --- Changes are made to the corner coefficients for neighboring
c       --- cells. These corner terms appear because of finite differences
c       --- across the neighboring cells. For example, the template for
c       --- cell (i+1,k) includes the finite difference
c       --- (phi(i,k-1)-phi(i,k+1)). This finite difference is modified
c       --- by the subgrid info at cell (i,k).
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          pxm = voltfac*volt(0,ic)
          rhs(ix,iz) = rhs(ix,iz) + coeffs(0,1,ix,iz)*pxm
          coeffs(0,1,ix,iz) = 0.
          if (iz > 0) then
            rhs(ix,iz-1) = rhs(ix,iz-1) - 0.5*chi(0,1,ix,iz)/(dz*dx2)*pxm
            coeffs(0,2,ix,iz-1) = 0.
          endif
          if (iz < nzlocal) then
            rhs(ix,iz+1) = rhs(ix,iz+1) + 0.5*chi(0,1,ix,iz)/(dz*dx2)*pxm
            coeffs(0,0,ix,iz+1) = 0.
          endif
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          coeffs(0,1,ix,iz) = 0.
          if (iz > 0) coeffs(0,2,ix,iz-1) = 0.
          if (iz < nzlocal) coeffs(0,0,ix,iz+1) = 0.
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          pxp = voltfac*volt(1,ic)
          rhs(ix,iz) = rhs(ix,iz) + coeffs(2,1,ix,iz)*pxp
          coeffs(2,1,ix,iz) = 0.
          if (iz > 0) then
            rhs(ix,iz-1) = rhs(ix,iz-1) + 0.5*chi(0,1,ix,iz)/(dz*dx2)*pxp
            coeffs(2,2,ix,iz-1) = 0.
          endif
          if (iz < nzlocal) then
            rhs(ix,iz+1) = rhs(ix,iz+1) - 0.5*chi(0,1,ix,iz)/(dz*dx2)*pxp
            coeffs(2,0,ix,iz+1) = 0.
          endif
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          coeffs(2,1,ix,iz) = 0.
          if (iz > 0) coeffs(2,2,ix,iz-1) = 0.
          if (iz < nzlocal) coeffs(2,0,ix,iz+1) = 0.
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          pzm = voltfac*volt(4,ic)
          rhs(ix,iz) = rhs(ix,iz) + coeffs(1,0,ix,iz)*pzm
          coeffs(1,0,ix,iz) = 0.
          if (ix > 1 .or. (ix == 1 .and. xmmin > 0.)) then
            rhs(ix-1,iz) = rhs(ix-1,iz)
     &         - 0.5*(rr(ix)-0.5*dx)/(rr(ix)-dx)*chi(1,0,ix,iz)/(dx*dz2)*pzm
            coeffs(2,0,ix-1,iz) = 0.
          endif
          if (ix < nx) then
            rhs(ix+1,iz) = rhs(ix+1,iz)
     &         + 0.5*(rr(ix)+0.5*dx)/(rr(ix)+dx)*chi(1,0,ix,iz)/(dx*dz2)*pzm
            coeffs(0,0,ix+1,iz) = 0.
          endif
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          coeffs(1,0,ix,iz) = 0.
          if (ix > 0) coeffs(2,0,ix-1,iz) = 0.
          if (ix < nx) coeffs(0,0,ix+1,iz) = 0.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          rhs(ix,iz) = rhs(ix,iz) + coeffs(1,2,ix,iz)*pzp
          coeffs(1,2,ix,iz) = 0.
          if (ix > 1 .or. (ix == 1 .and. xmmin > 0.)) then
            rhs(ix-1,iz) = rhs(ix-1,iz)
     &         + 0.5*(rr(ix)-0.5*dx)/(rr(ix)-dx)*chi(1,0,ix,iz)/(dx*dz2)*pzp
            coeffs(2,2,ix-1,iz) = 0.
          endif
          if (ix < nx) then
            rhs(ix+1,iz) = rhs(ix+1,iz)
     &         - 0.5*(rr(ix)+0.5*dx)/(rr(ix)+dx)*chi(1,0,ix,iz)/(dx*dz2)*pzp
            coeffs(0,2,ix+1,iz) = 0.
          endif
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          coeffs(1,2,ix,iz) = 0.
          if (ix > 0) coeffs(2,2,ix-1,iz) = 0.
          if (ix < nx) coeffs(0,2,ix+1,iz) = 0.
        endif

      enddo
!$OMP END DO

      end subroutine getcoefficientssubgrides2d
      end module coefficientsmodulees2d
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine mgsolveimplicites2d(iwhich,nx,nzlocal,nz,dx,dz,phi,rho,ns,qomdt,chi0,
     &                               bounds,xmmin,zmminlocal,zmmin,zbeam,zgrid,
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
      integer(ISZ):: nx,nzlocal,nz,ns
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nx,0:nzlocal)
      real(kind=8):: chi0(0:nx,0:nzlocal,0:ns-1)
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmmin,zmminlocal,zmmin,zbeam,zgrid
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
      bongrid => getbfieldsongrides2d(nx,nzlocal,nz,dx,dz,xmmin,zmminlocal,bounds)
c     --- Parallel B.C.s are set below.

c     --- Create copy of rho with extra guard cells, and fill it in.
      allocate(chi0temp(-1:nx+1,-3:nzlocal+3,0:ns-1))
      chi0temp(0:nx,0:nzlocal,:) = chi0
      call applyboundaryconditionses2d(chi0temp,nx,0,nzlocal,nz,1,0,3,ns,bounds,.true.)
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
          call setcndtr(xmmin,0.,zmminlocal,zmmin,zbeam,zgrid,nx,0,nzlocal,dx,dx,dz,
     &                  bounds(4),bounds(5),bounds(1),.false.,.false.)
        endif
        if (laddconductor) call execuser("calladdconductor")
      endif
      call checkconductors(nx,0,nzlocal,nz,dx,dx,dz,conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave)


!$OMP PARALLEL
!$OMP&PRIVATE(i,ix,iz)

c     --- Make sure guard planes have sensible values before beginning.
      call applyboundaryconditionses2d(phi,nx,0,nzlocal,nz,1,0,1,1,bounds,.false.)

      localbounds = bounds

#ifdef MPIPARALLEL
      lzoffset = 0
      rzoffset = 0
      if (izfsslave(my_index) > 0) localbounds(4) = -1
      if (izfsslave(my_index)+nzlocal < nz) localbounds(5) = -1
      call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                      bounds,nz,
     &                      lzoffset,rzoffset,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
c     --- These calls break the parallel field solver
c     call mgexchange_phi(nx,0,nzlocal,nz,phi,localbounds,-1,0,0,1,
c    &                    my_index,nslaves,izfsslave,nzfsslave,
c    &                    whosendingleft,izsendingleft,
c    &                    whosendingright,izsendingright)
c     call mgexchange_phiperiodic(nx,0,nzlocal,nz,phi,localbounds,1,0,0,1,
c    &                            my_index,nslaves,izfsslave,
c    &                            whosendingleft,whosendingright)
      call applyparallelboundaryconditionses2d(bongrid,nx,0,nzlocal,nz,
     &                                         1,0,3,3,1,bounds,
     &                                         my_index,nslaves,
     &                                         izfsslave,nzfsslave,
     &                                         whosendingleft,izsendingleft,
     &                                         whosendingright,izsendingright)
      call applyparallelboundaryconditionses2d(chi0temp,nx,0,nzlocal,nz,
     &                                         1,0,3,ns,1,bounds,
     &                                         my_index,nslaves,
     &                                         izfsslave,nzfsslave,
     &                                         whosendingleft,izsendingleft,
     &                                         whosendingright,izsendingright)
#endif

      allocate(phisave(-1:nx+1,-1:nzlocal+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"multigrid3dsolve: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nx,nzlocal
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
        call vcyclees2d(0,nx,nzlocal,nz,dx,dz,xmmin,phi,rho,ns,qomdt,chi0temp,bongrid,
     &                  bounds,mgparam,mgmaxlevels,
     &                  downpasses,uppasses,lcndbndy,icndbndy,conductors,lrz,
     &                  my_index,nslaves,izfsslave,nzfsslave)

c       --- Calculate the change in phi.
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
c       print*,mgiters,mgerror

c       --- This line below seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop above.
c       mgerror = maxval(abs(phisave - phi))

      enddo

c     --- For Dirichlet boundary conditions, copy data into guard planes
c     --- For other boundary conditions, the guard planes are used during
c     --- the solve are so are already set.
      if (localbounds(0) == 0) phi(-1,:) = 2.*phi(0,:) - phi(1,:)
      if (localbounds(1) == 0) phi(nx+1,:) = 2.*phi(nx,:) - phi(nx-1,:)
      if (localbounds(4) == 0) phi(:,-1) = 2.*phi(:,0) - phi(:,1)
      if (localbounds(5) == 0) phi(:,nzlocal+1) = 2.*phi(:,nzlocal) - phi(:,nzlocal-1)

c     --- Make a print out.
      if (mgerror > mgtol) then
        call remark("ImplicitMultigridES2d: Maximum number of iterations reached")
      endif
      if (mgverbose >= 1 .or. mgerror > mgtol) then
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
      RECURSIVE subroutine vcyclees2d(mglevel,nx,nzlocal,nz,dx,dz,xmmin,
     &                                phi,rho,ns,qomdt,chi0,bongrid,globalbounds,
     &                                mgparam,
     &                                mgmaxlevels,downpasses,uppasses,
     &                                lcndbndy,icndbndy,conductors,lrz,
     &                                my_index,nslaves,izfsslave,nzfsslave)
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use coefficientsmodulees2d
      use formggetarraysuminterface
c     use ImplicitMG2D
      integer(ISZ):: mglevel
      integer(ISZ):: nx,nzlocal,nz,ns
      real(kind=8):: dx,dz,xmmin
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1),rho(0:nx,0:nzlocal)
      real(kind=8):: chi0(-1:nx+1,-3:nzlocal+3,0:ns-1)
      real(kind=8):: bongrid(-1:nx+1,-3:nzlocal+3,0:2)
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
      real(kind=8),pointer:: coeffs(:,:,:,:),rhs(:,:)
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: chi0coarse(:,:,:),bongridcoarse(:,:,:)
      real(kind=8),allocatable:: res(:,:)
      integer(ISZ):: i,ic,iszone=1
      integer(ISZ):: nxcoarse,nzlocalcoarse,nzcoarse
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

      localbounds = globalbounds

#ifdef MPIPARALLEL
      lzoffsetall = 0
      rzoffsetall = 0
      call mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                      globalbounds,nz,
     &                      lzoffsetall,rzoffsetall,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
      if (izfsslave(my_index) > 0) localbounds(4) = -1
      if (izfsslave(my_index)+nzlocal < nz) localbounds(5) = -1
#endif

c     --- Get the finite difference coefficients for the Poisson operator at
c     --- the current level of refinement.
      allocate(coeffs(0:2,0:2,0:nx,0:nzlocal))
      allocate(rhs(0:nx,0:nzlocal))
      call getcoefficientses2d(nx,nzlocal,dx,dz,xmmin,ns,qomdt,rho,chi0,bongrid,
     &                         mglevel,lcndbndy,icndbndy,conductors,lrz,
     &                         localbounds,coeffs,rhs)
c     if (mglevel == 0) then
c       allocate(coeffs1(0:2,0:2,0:nx,0:nzlocal))
c       allocate(rhs1(0:nx,0:nzlocal))
c       coeffs1 = coeffs
c       rhs1 = rhs
c     endif

      if (larraysumdiagnostic) then
        sss = mggetarraysum(nx,0,nzlocal,1,0,1,phi,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V1 phi",mglevel,sss
        sss = mggetarraysum(nx,0,nzlocal,0,0,0,rho,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V1 rho",mglevel,sss
        sss = mggetarraysum(nx,0,nzlocal,1,0,3,chi0,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V1 chi0",mglevel,sss
      endif

c     --- Do initial SOR passes.
      do i=1,downpasses
        call relaxes2d(mglevel,nx,nzlocal,nz,phi,rho,rhs,coeffs,
     &                 localbounds,mgparam,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

      if (larraysumdiagnostic) then
        sss = mggetarraysum(nx,0,nzlocal,1,0,1,phi,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nx,0,nzlocal,0,0,0,rho,
     &                      izfsslave,nzfsslave,nslaves,my_index)
        if (my_index == 0) print*,"V2 rho",mglevel,sss
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. nz >= 4 .and. mglevel < mgmaxlevels) then

        allocate(res(-1:nx+1,-3:nzlocal+3),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nx,nzlocal
          call kaboom("vcyclees2d: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residuales2d(nx,nzlocal,nz,phi,rho,coeffs,rhs,res,
     &                    mglevel,localbounds,mgparam,conductors)
#ifdef MPIPARALLEL
        call applyparallelboundaryconditionses2d(res,nx,0,nzlocal,nz,1,0,3,1,3,
     &                                           localbounds,
     &                                           my_index,nslaves,
     &                                           izfsslave,nzfsslave,
     &                                           whosendingleft,izsendingleft,
     &                                           whosendingright,izsendingright)
#endif

        if (larraysumdiagnostic) then
          sss = mggetarraysum(nx,0,nzlocal,1,0,3,res,
     &                        izfsslave,nzfsslave,nslaves,my_index)
          if (my_index == 0) print*,"V3 res",mglevel,sss
        endif

c       --- Calculate the size of the next coarsest grid. If the current
c       --- size is twice and odd integer, the next coarsest size will be
c       --- n/2+1, gauranteeing that it is even.
c       --- Note that nzlocalcoarse is only used if full-coarsening is done.
        nxcoarse = nx/2
        if (mod(nxcoarse,2) == 1) nxcoarse = nxcoarse + 1
        nzlocalcoarse = nzlocal/2
        if (mod(nzlocalcoarse,2) == 1) nzlocalcoarse = nzlocalcoarse + 1
        nzcoarse = nz/2
        if (mod(nzcoarse,2) == 1) nzcoarse = nzcoarse + 1

        mingridsize = min(dx,dz)

        if (dx > 4./3.*mingridsize) then
          nxcoarse = nx
        endif
        if (dz > 4./3.*mingridsize) then
          nzlocalcoarse = nzlocal
          nzcoarse = nz
        endif

        dxcoarse = dx*nx/nxcoarse
        dzcoarse = dz*nz/nzcoarse
        dxcoarsesqi = 1./dxcoarse**2
        dzcoarsesqi = 1./dzcoarse**2

        localboundsc = globalbounds

#ifdef MPIPARALLEL
c       --- Find domains in coarser grid
        call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                  nz,nzcoarse)
c       --- Reset value to corrected one
        nzlocalcoarse = nzfsslavec(my_index)
c       --- Difference between starts and ends of coarse and fine grids.
c       --- Should only be in the range 0-2.
        lzoffsetall = (nzcoarse*izfsslave-nz*izfsslavec)
        rzoffsetall = (nz*(izfsslavec + nzfsslavec) -
     &                 nzcoarse*(izfsslave + nzfsslave))
c       --- Note that the lzoffsetall and rzoffsetall can only be used in
c       --- MPIPARALLEL sections since they will be unallocated in the
c       --- serial code. So, separate scalars are used in code which is
c       --- used in the serial version.
        lzoffset = lzoffsetall(my_index)
        rzoffset = rzoffsetall(my_index)
c       --- Get processor with which to exchange data on coarse grid
        call mggetexchangepes(nslaves,izfsslavec,nzfsslavec,my_index,
     &                        globalbounds,nzcoarse,
     &                        lzoffsetall,rzoffsetall,
     &                        whosendingleftc,izsendingleftc,
     &                        whosendingrightc,izsendingrightc)
        if (izfsslavec(my_index) > 0) localboundsc(4) = -1
        if (izfsslavec(my_index) + nzlocalcoarse < nzcoarse) localboundsc(5) = -1
#else
        lzoffset = 0
        rzoffset = 0
#endif

c       --- Alloate new work space
        allocate(phicoarse(-1:nxcoarse+1,-1:nzlocalcoarse+1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxcoarse,0:nzlocalcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(chi0coarse(-1:nxcoarse+1,-3:nzlocalcoarse+3,0:ns-1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate chi0coarse to shape ",
     &           nxcoarse,nzlocalcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(bongridcoarse(-1:nxcoarse+1,-3:nzlocalcoarse+3,0:2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate bongridcoarse to shape ",
     &           nxcoarse,nzlocalcoarse,3
          call kaboom("vcyclees2d: allocation error")
          return
        endif

c       rhocoarse = 0.
        phicoarse = 0.
c       chi0coarse = 0.
c       bongridcoarse = 0.

c       --- Restriction
        call restrictes3d(nx,0,nzlocal,nz,res,1,0,3,
     &                    nxcoarse,0,nzlocalcoarse,nzcoarse,rhocoarse,0,0,0,
     &                    localbounds,localboundsc,lzoffset)
        if (larraysumdiagnostic) then
          sss = mggetarraysum(nxcoarse,0,nzlocalcoarse,0,0,0,rhocoarse,
     &                        izfsslavec,nzfsslavec,nslaves,my_index)
          if (my_index == 0) print*,"V3 rhocoarse",mglevel,sss
          sss = mggetarraysum(nxcoarse,0,nzlocalcoarse,1,0,3,chi0coarse,
     &                        izfsslavec,nzfsslavec,nslaves,my_index)
          if (my_index == 0) print*,"V3 chi0coarse",mglevel,sss
        endif

        do ic=0,ns-1
          call restrictes3d(nx,0,nzlocal,nz,chi0(:,:,ic),1,0,3,
     &                      nxcoarse,0,nzlocalcoarse,nzcoarse,chi0coarse(:,:,ic),1,0,3,
     &                      localbounds,localboundsc,lzoffset)
        enddo
        do ic=0,2
          call restrictes3d(nx,0,nzlocal,nz,bongrid(:,:,ic),1,0,3,
     &                      nxcoarse,0,nzlocalcoarse,nzcoarse,bongridcoarse(:,:,ic),1,0,3,
     &                      localbounds,localboundsc,lzoffset)
        enddo
        call applyboundaryconditionses2d(chi0coarse,nxcoarse,0,nzlocalcoarse,nzcoarse,1,0,3,ns,localboundsc,.true.)
        call applyboundaryconditionses2d(bongridcoarse,nxcoarse,0,nzlocalcoarse,nzcoarse,1,0,3,3,localboundsc,.true.)
        call applyparallelboundaryconditionses2d(chi0coarse,nxcoarse,0,nzlocalcoarse,nzcoarse,
     &                                           1,0,3,ns,1,localboundsc,
     &                                           my_index,nslaves,
     &                                           izfsslavec,nzfsslavec,
     &                                           whosendingleftc,izsendingleftc,
     &                                           whosendingrightc,izsendingrightc)
        call applyparallelboundaryconditionses2d(bongridcoarse,nxcoarse,0,nzlocalcoarse,nzcoarse,
     &                                           1,0,3,1,1,localboundsc,
     &                                           my_index,nslaves,
     &                                           izfsslavec,nzfsslavec,
     &                                           whosendingleftc,izsendingleftc,
     &                                           whosendingrightc,izsendingrightc)

c       --- Continue at the next coarsest level.
        call vcyclees2d(mglevel+iszone,nxcoarse,nzlocalcoarse,nzcoarse,
     &                  dxcoarse,dzcoarse,xmmin,
     &                  phicoarse,rhocoarse,ns,qomdt,chi0coarse,bongridcoarse,
     &                  globalbounds,mgparam,
     &                  mgmaxlevels,downpasses,uppasses,
     &                  lcndbndy,icndbndy,conductors,lrz,
     &                  my_index,nslaves,izfsslavec,nzfsslavec)

        if (larraysumdiagnostic) then
          sss = mggetarraysum(nxcoarse,0,nzlocalcoarse,1,0,1,phicoarse,
     &                        izfsslavec,nzfsslavec,nslaves,my_index)
          if (my_index == 0) print*,"V4 phicoarse",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand3d(nx,0,nzlocal,nz,1,0,1,phi,
     &                nxcoarse,0,nzlocalcoarse,nzcoarse,phicoarse,
     &                localbounds,lzoffset)
        call applyboundaryconditionses2d(phi,nx,0,nzlocal,nz,1,0,1,1,localboundsc,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phiperiodic(nx,0,nzlocal,nz,phi,localbounds,1,1,0,1,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)
#endif
        call averageperiodicphi2d(phi,nx,0,nzlocal,nz,1,0,1,localbounds)
        if (larraysumdiagnostic) then
          sss = mggetarraysum(nx,0,nzlocal,1,0,1,phi,
     &                        izfsslave,nzfsslave,nslaves,my_index)
          if (my_index == 0) print*,"V5 phi",mglevel,sss
        endif
        
        deallocate(res)
        deallocate(phicoarse,rhocoarse)
        deallocate(chi0coarse,bongridcoarse)

      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call relaxes2d(mglevel,nx,nzlocal,nz,phi,rho,rhs,coeffs,
     &                 localbounds,mgparam,conductors,
     &                 my_index,nslaves,izfsslave,nzfsslave,
     &                 whosendingleft,izsendingleft,
     &                 whosendingright,izsendingright)
      enddo

      deallocate(coeffs)
      deallocate(rhs)

      return
      end
c=============================================================================
      subroutine relaxes2d(mglevel,nx,nzlocal,nz,phi,rho,rhs,coeffs,bounds,mgparam,
     &                     conductors,
     &                     my_index,nslaves,izfsslave,nzfsslave,
     &                     whosendingleft,izsendingleft,
     &                     whosendingright,izsendingright)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: mglevel,nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1),rho(0:nx,0:nzlocal)
      real(kind=8):: rhs(0:nx,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:2,0:nx,0:nzlocal)
      integer(ISZ):: bounds(0:5)
      real(kind=8):: mgparam
      type(ConductorType):: conductors
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.
c phisave is needed because when there are B fields, the even and odd sweeps
c are no longer independent since the diagonal coefficients become non-zero.
c The simplest solution is to save the values of phi before a sweep and
c read those saved values for the update of phi.

      real(kind=8),allocatable:: phisave(:,:)
      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ixmin,ixmax,izmin,izmax,ix,iz,ix1

      allocate(phisave(-1:nx+1,-1:nzlocal+1))

      call cond_potmg(conductors%interior,
     &                nx,0,nzlocal,1,0,1,phi,mglevel,1,.false.)
      call condbndymgint(conductors%evensubgrid,nx,0,nzlocal,1,0,1,phi,
     &                   bounds,mglevel,2)
      call condbndymgint(conductors%oddsubgrid,nx,0,nzlocal,1,0,1,phi,
     &                   bounds,mglevel,2)

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
      ixmin = 0
      ixmax = nx
      izmin = 0
      izmax = nzlocal
      if (bounds(0) < 1) ixmin = 1
      if (bounds(1) < 1) ixmax = nx - 1
      if (bounds(4) < 1) izmin = 1
      if (bounds(5) < 1) izmax = nzlocal - 1

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity
        phisave = phi

        do iz=izmin,izmax
          ix1 = ixmin + mod(ixmin + iz + parity,2)
          do ix=ix1,ixmax,2

            phi(ix,iz) = mgparam*(rho(ix,iz)/eps0 + rhs(ix,iz) +
     &                            coeffs(0,0,ix,iz)*phisave(ix-1,iz-1) +
     &                            coeffs(1,0,ix,iz)*phisave(ix  ,iz-1) +
     &                            coeffs(2,0,ix,iz)*phisave(ix+1,iz-1) +
     &                            coeffs(0,1,ix,iz)*phisave(ix-1,iz  ) +
     &                            coeffs(2,1,ix,iz)*phisave(ix+1,iz  ) +
     &                            coeffs(0,2,ix,iz)*phisave(ix-1,iz+1) +
     &                            coeffs(1,2,ix,iz)*phisave(ix  ,iz+1) +
     &                            coeffs(2,2,ix,iz)*phisave(ix+1,iz+1))/
     &                           (-coeffs(1,1,ix,iz)) +
     &                   (1.-mgparam)*phi(ix,iz)

          enddo
        enddo

        call cond_potmg(conductors%interior,
     &                  nx,0,nzlocal,1,0,1,phi,mglevel,1,.false.)

c       --- set phi in the z guard planes
        call applyboundaryconditionses2d(phi,nx,0,nzlocal,nz,1,0,1,1,bounds,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi(nx,0,nzlocal,nz,phi,bounds,0,1,0,1,
     &                      my_index,nslaves,izfsslave,nzfsslave,
     &                      whosendingleft,izsendingleft,
     &                      whosendingright,izsendingright)
        call mgexchange_phiperiodic(nx,0,nzlocal,nz,phi,bounds,1,1,0,1,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- Exchange phi in the z guard planes
      call mgexchange_phi(nx,0,nzlocal,nz,phi,bounds,-1,1,0,1,
     &                    my_index,nslaves,izfsslave,nzfsslave,
     &                    whosendingleft,izsendingleft,
     &                    whosendingright,izsendingright)
#endif

      call averageperiodicphi2d(phi,nx,0,nzlocal,nz,1,0,1,bounds)

      deallocate(phisave)

      return
      end
c=============================================================================
      subroutine residuales2d(nx,nzlocal,nz,phi,rho,coeffs,rhs,res,
     &                        mglevel,bounds,mgparam,conductors)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nx,ny,nzlocal,nz
      real(kind=8):: phi(-1:nx+1,-1:nzlocal+1),rho(0:nx,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:2,0:nx,0:nzlocal)
      real(kind=8):: rhs(0:nx,0:nzlocal)
      real(kind=8):: res(-1:nx+1,-3:nzlocal+3)
      integer(ISZ):: mglevel,bounds(0:5)
      real(kind=8):: mgparam
      type(ConductorType):: conductors

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixmin,ixmax,izmin,izmax

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nx
      izmin = 0
      izmax = nzlocal
      if (bounds(0) == 0) ixmin = 1
      if (bounds(1) == 0) ixmax = nx - 1
      if (bounds(4) == 0) izmin = 1
      if (bounds(5) == 0) izmax = nzlocal - 1

      res = 0.
c     --- Calculate the residual.
!$OMP DO
      do iz=izmin,izmax

        do ix=ixmin,ixmax

          res(ix,iz) = rho(ix,iz) + eps0*(rhs(ix,iz) +
     &                  coeffs(0,0,ix,iz)*phi(ix-1,iz-1) +
     &                  coeffs(1,0,ix,iz)*phi(ix  ,iz-1) +
     &                  coeffs(2,0,ix,iz)*phi(ix+1,iz-1) +
     &                  coeffs(0,1,ix,iz)*phi(ix-1,iz  ) +
     &                  coeffs(1,1,ix,iz)*phi(ix  ,iz  ) +
     &                  coeffs(2,1,ix,iz)*phi(ix+1,iz  ) +
     &                  coeffs(0,2,ix,iz)*phi(ix-1,iz+1) +
     &                  coeffs(1,2,ix,iz)*phi(ix  ,iz+1) +
     &                  coeffs(2,2,ix,iz)*phi(ix+1,iz+1))

        enddo
      enddo
!$OMP END DO

      call cond_potmgres(conductors%interior,
     &                   nx,0,nzlocal,res,mglevel,1,.false.,1,0)
      call applyboundaryconditionses2d(res,nx,0,nzlocal,nz,1,0,3,1,bounds,
     &                                 .true.)

      return
      end
c=============================================================================
      subroutine restrictes3d(nx,ny,nzlocal,nz,u,delx,dely,delz,
     &                        nxcoarse,nycoarse,nzlocalcoarse,nzcoarse,ucoarse,
     &                        delcx,delcy,delcz,
     &                        bounds,boundscoarse,lzoffset)

      integer(ISZ):: nx,ny,nzlocal,nz,delx,dely,delz,delcx,delcy,delcz
      integer(ISZ):: nxcoarse,nycoarse,nzlocalcoarse,nzcoarse
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz)
      real(kind=8):: ucoarse(-delcx:nxcoarse+delcx,-delcy:nycoarse+delcy,-delcz:nzlocalcoarse+delcz)
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

c     --- Set the loop limits, skipping Dirichlet edges
      ixcoarsemin = 0
      ixcoarsemax = nxcoarse
      iycoarsemin = 0
      iycoarsemax = nycoarse
      izcoarsemin = 0
      izcoarsemax = nzlocalcoarse
      if (bounds(0) == 0) ixcoarsemin = 1
      if (bounds(1) == 0) ixcoarsemax = nxcoarse - 1
      if (bounds(2) == 0) iycoarsemin = 1
      if (bounds(3) == 0) iycoarsemax = nycoarse - 1
      if (bounds(4) == 0) izcoarsemin = 1
      if (bounds(5) == 0) izcoarsemax = nzlocalcoarse - 1

      dx = 1.*nx/nxcoarse
      dxi = 1.*nxcoarse/nx
      if (ny > 0) then
        dy = 1.*ny/nycoarse
        dyi = 1.*nycoarse/ny
      else
        dy = 1.
        dyi = 1.
      endif
      dz = 1.*nz/nzcoarse
      dzi = 1.*nzcoarse/nz

c     --- Do the loops.
!$OMP DO
      do izcoarse=izcoarsemin,izcoarsemax
        izmin = ((izcoarse-1)*nz - lzoffset + 4*nzcoarse)/nzcoarse-3
        izmax = ((izcoarse+1)*nz - lzoffset - 1)/nzcoarse
        if (izmin < -delz) izmin = -delz
        if (izmax > nzlocal+delz) izmax = nzlocal+delz

        do iz=izmin,izmax
          wz(iz-izmin) = 1. - abs(izcoarse - (iz + 1.*lzoffset/nzcoarse)*dzi)
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
      subroutine applyboundaryconditionses2d(u,nx,ny,nzlocal,nz,
     &                                       delx,dely,delz,ncomp,bounds,
     &                                       lzerodirichlet)
      integer(ISZ):: nx,ny,nzlocal,nz,delx,dely,delz,ncomp
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz,ncomp)
      integer(ISZ):: bounds(0:5)
      logical(ISZ):: lzerodirichlet

      integer(ISZ):: iz

c     --- Apply the boundary conditions, setting the values in the guard cells
      if (delx > 0) then
        if (bounds(0) == 0) then
          if (lzerodirichlet) then
            u(-1,0:ny,0:nzlocal,:) = 0.
          else
            u(-1,0:ny,0:nzlocal,:) = 2.*u(0,0:ny,0:nzlocal,:) - u(1,0:ny,0:nzlocal,:)
          endif
        endif
        if (bounds(1) == 0) then
          if (lzerodirichlet) then
            u(nx+1,0:ny,0:nzlocal,:) = 0.
          else
            u(nx+1,0:ny,0:nzlocal,:)=2.*u(nx,0:ny,0:nzlocal,:)-u(nx-1,0:ny,0:nzlocal,:)
          endif
        endif
        if (bounds(0) == 1) u(-1,0:ny,0:nzlocal,:) = u(+1,0:ny,0:nzlocal,:)
        if (bounds(1) == 1) u(nx+1,0:ny,0:nzlocal,:) = u(nx-1,0:ny,0:nzlocal,:)
        if (bounds(0) == 2) u(-1,0:ny,0:nzlocal,:) = u(nx-1,0:ny,0:nzlocal,:)
        if (bounds(1) == 2) u(nx+1,0:ny,0:nzlocal,:) = u(+1,0:ny,0:nzlocal,:)
      endif

      if (dely > 0) then
        if (bounds(2) == 0) then
          if (lzerodirichlet) then
            u(:,-1,0:nzlocal,:) = 2.*u(:,0,0:nzlocal,:) - u(:,1,0:nzlocal,:)
          else
            u(:,-1,0:nzlocal,:) = 0.
          endif
        endif
        if (bounds(3) == 0) then
          if (lzerodirichlet) then
            u(:,ny+1,0:nzlocal,:) = 0.
          else
            u(:,ny+1,0:nzlocal,:) = 2.*u(:,ny,0:nzlocal,:) - u(:,ny-1,0:nzlocal,:)
          endif
        endif
        if (bounds(2) == 1) u(:,-1,0:nzlocal,:) = u(:,+1,0:nzlocal,:)
        if (bounds(3) == 1) u(:,ny+1,0:nzlocal,:) = u(:,ny-1,0:nzlocal,:)
        if (bounds(2) == 2) u(:,-1,0:nzlocal,:) = u(:,ny-1,0:nzlocal,:)
        if (bounds(3) == 2) u(:,ny+1,0:nzlocal,:) = u(:,+1,0:nzlocal,:)
      endif

      if (delz > 0) then
        if (bounds(4) == 0) then
          if (lzerodirichlet) then
            u(:,:,-delz:-1,:) = 0.
          else
            do iz=-1,-delz,-1
              u(:,:,iz,:) = 2.*u(:,:,iz+1,:) - u(:,:,iz+2,:)
            enddo
          endif
        endif
        if (bounds(5) == 0) then
          if (lzerodirichlet) then
            u(:,:,nzlocal+1:nzlocal+delz,:) = 0.
          else
            do iz=nzlocal+1,nzlocal+delz
              u(:,:,iz,:) = 2.*u(:,:,iz-1,:) - u(:,:,iz-2,:)
            enddo
          endif
        endif
        if (bounds(4) == 1) u(:,:,-delz:-1,:) = u(:,:,delz:+1:-1,:)
        if (bounds(5) == 1)
     &    u(:,:,nzlocal+1:nzlocal+delz,:) = u(:,:,nzlocal-1:nzlocal-delz:-1,:)
        if (bounds(4) == 2 .and. nzlocal == nz)
     &    u(:,:,-delz:-1,:) = u(:,:,nzlocal-delz:nzlocal-1,:)
        if (bounds(5) == 2 .and. nzlocal == nz)
     &    u(:,:,nzlocal+1:nzlocal+delz,:) = u(:,:,+1:delz,:)
      endif

c     --- Need code to set z guard planes for periodic B.C.'s for parallel

      return
      end
c=============================================================================
      subroutine applyparallelboundaryconditionses2d(u,nx,ny,nzlocal,nz,
     &                                               delx,dely,delz,ncomp,zsend,bounds,
     &                                               my_index,nslaves,
     &                                               izfsslave,nzfsslave,
     &                                               whosendingleft,izsendingleft,
     &                                               whosendingright,izsendingright)
      integer(ISZ):: nx,ny,nzlocal,nz,delx,dely,delz,ncomp,zsend
      real(kind=8):: u(-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz,ncomp)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

      integer(ISZ):: iz,ic

#ifdef MPIPARALLEL

      do ic=1,ncomp
        do iz=1,zsend
          call mgexchange_phi(nx,ny,nzlocal,nz,u(:,:,:,ic),bounds,-iz,delx,dely,delz,
     &                        my_index,nslaves,izfsslave,nzfsslave,
     &                        whosendingleft,izsendingleft,
     &                        whosendingright,izsendingright)
        enddo

        call mgexchange_phiperiodic(nx,ny,nzlocal,nz,u(:,:,:,ic),bounds,0,delx,dely,delz,
     &                              my_index,nslaves,izfsslave,
     &                              whosendingleft,whosendingright)

      enddo

#endif

      return
      end
c=============================================================================
      subroutine averageperiodicphi2d(phi,nx,ny,nzlocal,nz,
     &                                delx,dely,delz,bounds)
      integer(ISZ):: nx,ny,nzlocal,nz,delx,dely,delz
      real(kind=8):: phi(-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz)
      integer(ISZ):: bounds(0:5)

c Subtracts out the average phi.
c Warning - does not work in parallel.

      real(kind=8):: phiave

      if (bounds(0) == 2 .and. bounds(1) == 2 .and.
     &    bounds(4) == 2 .and. bounds(5) == 2) then

        phiave = sum(phi(0:nx-1,0,0:nzlocal-1))/(nx*nz)
        phi = phi - phiave

      endif

      return
      end
c=============================================================================
