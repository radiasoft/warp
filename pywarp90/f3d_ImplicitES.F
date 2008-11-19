#include "top.h"
c=============================================================================
c@(#) Filef3d_ImplicitES.F, version $Revision: 1.8 $, $Date: 2008/11/19 18:29:49 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the 3D multigrid implicit electrostatic field solver.
c  David P. Grote, LLNL, (925)423-7194, LBNL (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      module coefficientsmodulees3d
      contains
      function getbfieldsongrid3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                            dx,dy,dz,
     &                            xmminlocal,ymminlocal,zmminlocal,zgrid,
     &                            localbounds)
     &  result(bongrid)
      use Subtimersf3d
      real(kind=8),pointer:: bongrid(:,:,:,:)
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz,xmminlocal,ymminlocal,zmminlocal,zgrid
      integer(ISZ):: localbounds(0:5)
      integer(ISZ):: delx,dely,delz

      integer(ISZ):: ix,iy,iz
      real(kind=8),allocatable:: xx(:),yy(:),zz(:),uz(:),gaminv(:)
      real(kind=8),allocatable:: ex(:),ey(:),ez(:),bendres(:),bendradi(:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Allocate space for B field
      allocate(bongrid(-delx:nxlocal+delx,
     &                 -dely:nylocal+dely,
     &                 -delz:nzlocal+delz,0:2))

c     --- Create temporary space to hold the coordinates of the grid
c     --- needed to look up the B fields.
      allocate(xx(0:nxlocal),yy(0:nxlocal),zz(0:nxlocal))
      allocate(uz(0:nxlocal),gaminv(0:nxlocal))
      allocate(ex(0:nxlocal),ey(0:nxlocal),ez(0:nxlocal))
      allocate(bendres(0:nxlocal),bendradi(0:nxlocal))

c     --- Fill the transverse position arrays.
      do ix=0,nxlocal
        xx(ix) = xmminlocal + ix*dx
      enddo
      uz = -1.
      gaminv = 1.
      bendres = 0.
      bendradi = LARGEPOS
      ex = 0.
      ey = 0.
      ez = 0.
      bongrid = 0.

c     --- Get the B field, one x line at a time. One x line is done at a
c     --- time to save memory.
      do iz=0,nzlocal
        zz = zmminlocal + iz*dz + zgrid
        do iy=0,nylocal
          yy = ymminlocal + iy*dy
          call exteb3d(1+nxlocal,xx,yy,zz,uz,gaminv,0.,0.,
     &                 bongrid(0:nxlocal,iy,iz,0),
     &                 bongrid(0:nxlocal,iy,iz,1),
     &                 bongrid(0:nxlocal,iy,iz,2),
     &                 ex,ey,ez,1.,1.,bendres,bendradi,1.,1.)
        enddo
      enddo

c     --- Apply the boundary conditions, setting the values in the guard cells.
      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                               bongrid,3,localbounds,.false.,.true.)

      deallocate(xx,yy,zz,uz,gaminv)
      deallocate(ex,ey,ez)
      deallocate(bendres,bendradi)

      if (lf3dtimesubs) timegetbfieldsongrid3d = timegetbfieldsongrid3d +
     &                                             wtime() - substarttime

      return
      end function getbfieldsongrid3d
c=============================================================================
      subroutine getcoefficientsimplicites3d(nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                                       ns,qomdt,chi0,bongrid,
     &                                       mglevel,lcndbndy,icndbndy,
     &                                       conductors,localbounds,
     &                                       coeffs0,coeffs,rhs,resscale,
     &                                       delx,dely,delz)
      use Subtimersf3d
      use Constant
      use ConductorTypemodule
      use formggetarraysuminterface
      integer(ISZ):: nxlocal,nylocal,nzlocal,ns
      real(kind=8):: dx,dy,dz
      integer(ISZ):: delx,dely,delz
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: chi0(-delx:nxlocal+delx,
     &                    -dely:nylocal+dely,
     &                    -delz:nzlocal+delz,0:ns-1)
      real(kind=8):: bongrid(-delx:nxlocal+delx,
     &                       -dely:nylocal+dely,
     &                       -delz:nzlocal+delz,0:2)
      integer(ISZ):: mglevel
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lrz
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: coeffs0(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:1,-1:1,0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: resscale(0:nxlocal,0:nylocal,0:nzlocal)

c Layout of coeffs0 and coeffs:
c  - coeffs0 is the coefficient of the central term in the template, phiijk
c  - coeffs, 1st dimension is direction, x(0), y(1), or z(2)
c            2nd dimension is offset in the direction, -1(0), +1(1)
c            3rd dimension is offset in the next direction, -1,0,+1
c    for example, coeffs(0,0,0,...) is for phi(i-1,j+0,k)
c                 coeffs(0,1,0,...) is for phi(i+1,j+0,k)
c                 coeffs(0,0,-1,...) is for phi(i-1,j-1,k)
c                 coeffs(2,0,0,...) is for phi(i+0,j,k-1)

      real(kind=8),pointer:: chi(:,:,:,:,:)
      integer(ISZ):: ix,iy,iz,js,ic
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
      real(kind=8):: ox,oy,oz,oo
      real(kind=8):: rr,rrm,rrp
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      allocate(chi(0:2,0:2,-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1))
      chi = 0.

c     --- First, generate the chi tensor from chi0 and the B field.
      do js=0,ns-1
        do iz=-1,nzlocal+1
          do iy=-1,nylocal+1
            do ix=-1,nxlocal+1

              ox = 0.5*qomdt(js)*bongrid(ix,iy,iz,0)
              oy = 0.5*qomdt(js)*bongrid(ix,iy,iz,1)
              oz = 0.5*qomdt(js)*bongrid(ix,iy,iz,2)
              oo = chi0(ix,iy,iz,js)/(1. + ox**2 + oy**2 + oz**2)

              chi(0,0,ix,iy,iz) = chi(0,0,ix,iy,iz) + oo*(1. + ox**2)
              chi(1,0,ix,iy,iz) = chi(1,0,ix,iy,iz) + oo*(ox*oy - oz)
              chi(2,0,ix,iy,iz) = chi(2,0,ix,iy,iz) + oo*(ox*oz + oy)

              chi(0,1,ix,iy,iz) = chi(0,1,ix,iy,iz) + oo*(ox*oy + oz)
              chi(1,1,ix,iy,iz) = chi(1,1,ix,iy,iz) + oo*(1. + oy**2)
              chi(2,1,ix,iy,iz) = chi(2,1,ix,iy,iz) + oo*(oy*oz - ox)

              chi(0,2,ix,iy,iz) = chi(0,2,ix,iy,iz) + oo*(ox*oz - oy)
              chi(1,2,ix,iy,iz) = chi(1,2,ix,iy,iz) + oo*(oy*oz + ox)
              chi(2,2,ix,iy,iz) = chi(2,2,ix,iy,iz) + oo*(1. + oz**2)

            enddo
          enddo
        enddo

c       --- Force chi to zero inside of any conductors. This is because
c       --- the charge density is always zero there.
c       --- This is probably not needed.
        do ic = conductors%interior%istart(mglevel),
     &          conductors%interior%istart(mglevel+1)-1
          ix = conductors%interior%indx(0,ic)
          iy = conductors%interior%indx(1,ic)
          iz = conductors%interior%indx(2,ic)
          chi(:,:,ix,iy,iz) = 0.
        enddo

      enddo

      do iz=0,nzlocal
        do iy=0,nylocal
          do ix=0,nxlocal

            coeffs(0,1,0,ix,iy,iz) = (1.+0.5*(chi(0,0,ix  ,iy  ,iz  ) +
     &                                        chi(0,0,ix+1,iy  ,iz  )))/dx**2
            coeffs(0,0,0,ix,iy,iz) = (1.+0.5*(chi(0,0,ix-1,iy  ,iz  ) +
     &                                        chi(0,0,ix  ,iy  ,iz  )))/dx**2

            coeffs(1,1,0,ix,iy,iz) = (1.+0.5*(chi(1,1,ix  ,iy  ,iz  ) +
     &                                        chi(1,1,ix  ,iy+1,iz  )))/dy**2
            coeffs(1,0,0,ix,iy,iz) = (1.+0.5*(chi(1,1,ix  ,iy-1,iz  ) +
     &                                        chi(1,1,ix  ,iy  ,iz  )))/dy**2

            coeffs(2,1,0,ix,iy,iz) = (1.+0.5*(chi(2,2,ix  ,iy  ,iz  ) +
     &                                        chi(2,2,ix  ,iy  ,iz+1)))/dz**2
            coeffs(2,0,0,ix,iy,iz) = (1.+0.5*(chi(2,2,ix  ,iy  ,iz-1) +
     &                                        chi(2,2,ix  ,iy  ,iz  )))/dz**2

            coeffs0(ix,iy,iz) = -coeffs(0,1,0,ix,iy,iz)
     &                          -coeffs(0,0,0,ix,iy,iz)
     &                          -coeffs(1,1,0,ix,iy,iz)
     &                          -coeffs(1,0,0,ix,iy,iz)
     &                          -coeffs(2,1,0,ix,iy,iz)
     &                          -coeffs(2,0,0,ix,iy,iz)


            coeffs(0,1,+1,ix,iy,iz) = +0.5*chi(0,1,ix+1,iy  ,iz  )/(2.*dy*dx)
     &                                +0.5*chi(1,0,ix  ,iy+1,iz  )/(2.*dx*dy)
            coeffs(0,0,+1,ix,iy,iz) = -0.5*chi(0,1,ix-1,iy  ,iz  )/(2.*dy*dx)
     &                                -0.5*chi(1,0,ix  ,iy+1,iz  )/(2.*dx*dy)
            coeffs(0,1,-1,ix,iy,iz) = -0.5*chi(0,1,ix+1,iy  ,iz  )/(2.*dy*dx)
     &                                -0.5*chi(1,0,ix  ,iy-1,iz  )/(2.*dx*dy)
            coeffs(0,0,-1,ix,iy,iz) = +0.5*chi(0,1,ix-1,iy  ,iz  )/(2.*dy*dx)
     &                                +0.5*chi(1,0,ix  ,iy-1,iz  )/(2.*dx*dy)

            coeffs(1,1,+1,ix,iy,iz) = +0.5*chi(1,2,ix  ,iy+1,iz  )/(2.*dz*dy)
     &                                +0.5*chi(2,1,ix  ,iy  ,iz+1)/(2.*dy*dz)
            coeffs(1,0,+1,ix,iy,iz) = -0.5*chi(1,2,ix  ,iy-1,iz  )/(2.*dz*dy)
     &                                -0.5*chi(2,1,ix  ,iy  ,iz+1)/(2.*dy*dz)
            coeffs(1,1,-1,ix,iy,iz) = -0.5*chi(1,2,ix  ,iy+1,iz  )/(2.*dz*dy)
     &                                -0.5*chi(2,1,ix  ,iy  ,iz-1)/(2.*dy*dz)
            coeffs(1,0,-1,ix,iy,iz) = +0.5*chi(1,2,ix  ,iy-1,iz  )/(2.*dz*dy)
     &                                +0.5*chi(2,1,ix  ,iy  ,iz-1)/(2.*dy*dz)

            coeffs(2,1,+1,ix,iy,iz) = +0.5*chi(2,0,ix  ,iy  ,iz+1)/(2.*dx*dz)
     &                                +0.5*chi(0,2,ix+1,iy  ,iz  )/(2.*dz*dx)
            coeffs(2,0,+1,ix,iy,iz) = -0.5*chi(2,0,ix  ,iy  ,iz-1)/(2.*dx*dz)
     &                                -0.5*chi(0,2,ix+1,iy  ,iz  )/(2.*dz*dx)
            coeffs(2,1,-1,ix,iy,iz) = -0.5*chi(2,0,ix  ,iy  ,iz+1)/(2.*dx*dz)
     &                                -0.5*chi(0,2,ix-1,iy  ,iz  )/(2.*dz*dx)
            coeffs(2,0,-1,ix,iy,iz) = +0.5*chi(2,0,ix  ,iy  ,iz-1)/(2.*dx*dz)
     &                                +0.5*chi(0,2,ix-1,iy  ,iz  )/(2.*dz*dx)

          enddo
        enddo
      enddo

      rhs = 0.
      resscale = 1.
      call getcoefficientssubgridimplicites3d(conductors%evensubgrid,
     &                                nxlocal,nylocal,nzlocal,
     &                                dx,dy,dz,mglevel,lcndbndy,icndbndy,chi,
     &                                localbounds,coeffs0,coeffs,rhs,resscale)
      call getcoefficientssubgridimplicites3d(conductors%oddsubgrid,
     &                                nxlocal,nylocal,nzlocal,
     &                                dx,dy,dz,mglevel,lcndbndy,icndbndy,chi,
     &                                localbounds,coeffs0,coeffs,rhs,resscale)

      deallocate(chi)

      if (lf3dtimesubs) timegetcoefficientsimplicites3d =
     &                  timegetcoefficientsimplicites3d +
     &                                            wtime() - substarttime

      return
      end subroutine getcoefficientsimplicites3d
c=============================================================================
      subroutine getcoefficientssubgridimplicites3d(subgrid,
     &                                nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                                mglevel,lcndbndy,icndbndy,chi,
     &                                localbounds,coeffs0,coeffs,rhs,resscale)
      use Subtimersf3d
      use Constant
      use ConductorSubgridTypemodule
      type(ConductorSubgridType):: subgrid
      integer(ISZ):: nxlocal,nylocal,nzlocal
      real(kind=8):: dx,dy,dz
      integer(ISZ):: mglevel
      integer(ISZ):: icndbndy
      logical(ISZ):: lcndbndy,lrz
      real(kind=8):: chi(0:2,0:2,-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: coeffs0(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:1,-1:1,0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: resscale(0:nxlocal,0:nylocal,0:nzlocal)

c Uses adjusted difference equation to enforce sub-grid level placement of 
c conductor boundaries for points near conductor surface.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: dxm,dym,dzm,dxp,dyp,dzp
      real(kind=8):: pxm,pym,pzm,pxp,pyp,pzp
      real(kind=8):: chixm(0:2,0:2),chiym(0:2,0:2),chizm(0:2,0:2)
      real(kind=8):: chixp(0:2,0:2),chiyp(0:2,0:2),chizp(0:2,0:2)
      real(kind=8):: cxm,cym,czm,cxp,cyp,czp
      real(kind=8):: dx1,dx2,dy1,dy2,dz1,dz2
      real(kind=8):: voltfac,c0,ppp
      real(kind=8),pointer:: dels(:,:),volt(:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (.not. lcndbndy) return

      if (icndbndy == 1) then
        call kaboom("mgsolveimplicites3d: icndbndy == 1 not supported")
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

c     --- Loop over points near surface of conductors.
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet.
c       --- Note that data in parallel boundaries (where localbounds==-1)
c       --- is handled since those coefficients are needed in the calculation
c       --- of the residual.
        if (ix == 0  .and. localbounds(0) == 0) cycle
        if (ix == nxlocal .and. localbounds(1) == 0) cycle
        if (iy == 0  .and. localbounds(2) == 0) cycle
        if (iy == nylocal .and. localbounds(3) == 0) cycle
        if (iz == 0  .and. localbounds(4) == 0) cycle
        if (iz == nzlocal .and. localbounds(5) == 0) cycle

        ppp = 1.

c       --- First, get the effective grid cell sizes and the chi's.
c       --- The chi's inside the conductor are always zero because the
c       --- charge density is zero there.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          dxm = dels(0,ic)
          chixm = 0.
          ppp = min(ppp,dxm)
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
          chixm = 0.
          ppp = min(ppp,1.-1.e-9)
        else
          dxm = 1.
          chixm = chi(:,:,ix-1,iy,iz)
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          dxp = dels(1,ic)
          chixp = 0.
          ppp = min(ppp,dxp)
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
          chixp = 0.
          ppp = min(ppp,1.-1.e-9)
        else
          dxp = 1.
          chixp = chi(:,:,ix+1,iy,iz)
        endif

        if (0. < dels(2,ic) .and. dels(2,ic) < +1.) then
          dym = dels(2,ic)
          chiym = 0.
          ppp = min(ppp,dym)
        elseif (-1. < dels(2,ic) .and. dels(2,ic) <= 0.) then
          dym = -2.*dels(2,ic)
          chiym = 0.
          ppp = min(ppp,1.-1.e-9)
        else
          dym = 1.
          chiym = chi(:,:,ix,iy-1,iz)
        endif

        if (0. < dels(3,ic) .and. dels(3,ic) < +1.) then
          dyp = dels(3,ic)
          chiyp = 0.
          ppp = min(ppp,dyp)
        elseif (-1. < dels(3,ic) .and. dels(3,ic) <= 0.) then
          dyp = -2.*dels(3,ic)
          chiyp = 0.
          ppp = min(ppp,1.-1.e-9)
        else
          dyp = 1.
          chiyp = chi(:,:,ix,iy+1,iz)
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          dzm = dels(4,ic)
          chizm = 0.
          ppp = min(ppp,dzm)
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
          chizm = 0.
          ppp = min(ppp,1.-1.e-9)
        else
          dzm = 1.
          chizm = chi(:,:,ix,iy,iz-1)
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          dzp = dels(5,ic)
          chizp = 0.
          ppp = min(ppp,dzp)
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
          chizp = 0.
          ppp = min(ppp,1.-1.e-9)
        else
          dzp = 1.
          chizp = chi(:,:,ix,iy,iz+1)
        endif

        resscale(ix,iy,iz) = ppp

c       --- Now construct the coefficients
        cxm = 1./(dxm*(0.5*dxm + 0.5*dxp))
        cxp = 1./(dxp*(0.5*dxm + 0.5*dxp))
        cym = 1./(dym*(0.5*dym + 0.5*dyp))
        cyp = 1./(dyp*(0.5*dym + 0.5*dyp))
        czm = 1./(dzm*(0.5*dzm + 0.5*dzp))
        czp = 1./(dzp*(0.5*dzm + 0.5*dzp))
        if (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) cxm = 0.
        if (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) cxp = 0.
        if (-1. < dels(2,ic) .and. dels(2,ic) <= 0.) cym = 0.
        if (-1. < dels(3,ic) .and. dels(3,ic) <= 0.) cyp = 0.
        if (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) czm = 0.
        if (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) czp = 0.

c       --- The dx1 etc are for finite differences (i+1/2) - (i-1/2)
c       --- the dx2 etc are for finite differences (i+1) - (i-1)
        dx1 = dx*(0.5*dxm + 0.5*dxp)
        dy1 = dy*(0.5*dym + 0.5*dyp)
        dz1 = dz*(0.5*dzm + 0.5*dzp)
        dx2 = dx*(dxm + dxp)
        dy2 = dy*(dym + dyp)
        dz2 = dz*(dzm + dzp)

        coeffs(0,1,0,ix,iy,iz)=cxp*(1.+0.5*(chi(0,0,ix,iy,iz)+chixp(0,0)))/dx**2
        coeffs(0,0,0,ix,iy,iz)=cxm*(1.+0.5*(chixm(0,0)+chi(0,0,ix,iy,iz)))/dx**2
        coeffs(1,1,0,ix,iy,iz)=cyp*(1.+0.5*(chi(1,1,ix,iy,iz)+chiyp(1,1)))/dy**2
        coeffs(1,0,0,ix,iy,iz)=cym*(1.+0.5*(chiym(1,1)+chi(1,1,ix,iy,iz)))/dy**2
        coeffs(2,1,0,ix,iy,iz)=czp*(1.+0.5*(chi(2,2,ix,iy,iz)+chizp(2,2)))/dz**2
        coeffs(2,0,0,ix,iy,iz)=czm*(1.+0.5*(chizm(2,2)+chi(2,2,ix,iy,iz)))/dz**2

        coeffs0(ix,iy,iz) = -coeffs(0,1,0,ix,iy,iz)-coeffs(0,0,0,ix,iy,iz)
     &                      -coeffs(1,1,0,ix,iy,iz)-coeffs(1,0,0,ix,iy,iz)
     &                      -coeffs(2,1,0,ix,iy,iz)-coeffs(2,0,0,ix,iy,iz)

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
          rhs(ix,iy,iz) = rhs(ix,iy,iz) + coeffs(0,0,0,ix,iy,iz)*pxm
          coeffs(0,0,0,ix,iy,iz) = 0.
          if (iy > 0) then
            rhs(ix,iy-1,iz)=rhs(ix,iy-1,iz) - 0.5*chi(1,0,ix,iy,iz)/(dy*dx2)*pxm
            coeffs(0,0,+1,ix,iy-1,iz) = 0.
          endif
          if (iy < nylocal) then
            rhs(ix,iy+1,iz)=rhs(ix,iy+1,iz) + 0.5*chi(1,0,ix,iy,iz)/(dy*dx2)*pxm
            coeffs(0,0,-1,ix,iy+1,iz) = 0.
          endif
          if (iz > 0) then
            rhs(ix,iy,iz-1)=rhs(ix,iy,iz-1) - 0.5*chi(2,0,ix,iy,iz)/(dz*dx2)*pxm
            coeffs(2,1,-1,ix,iy,iz-1) = 0.
          endif
          if (iz < nzlocal) then
            rhs(ix,iy,iz+1)=rhs(ix,iy,iz+1) + 0.5*chi(2,0,ix,iy,iz)/(dz*dx2)*pxm
            coeffs(2,0,-1,ix,iy,iz+1) = 0.
          endif
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          coeffs(0,0,0,ix,iy,iz) = 0.
          if (iy > 0) coeffs(0,0,+1,ix,iy-1,iz) = 0.
          if (iy < nylocal) coeffs(0,0,-1,ix,iy+1,iz) = 0.
          if (iz > 0) coeffs(2,1,-1,ix,iy,iz-1) = 0.
          if (iz < nzlocal) coeffs(2,0,-1,ix,iy,iz+1) = 0.
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          pxp = voltfac*volt(1,ic)
          rhs(ix,iy,iz) = rhs(ix,iy,iz) + coeffs(0,1,0,ix,iy,iz)*pxp
          coeffs(0,1,0,ix,iy,iz) = 0.
          if (iy > 0) then
            rhs(ix,iy-1,iz)=rhs(ix,iy-1,iz) + 0.5*chi(1,0,ix,iy,iz)/(dy*dx2)*pxp
            coeffs(0,1,+1,ix,iy-1,iz) = 0.
          endif
          if (iy < nylocal) then
            rhs(ix,iy+1,iz)=rhs(ix,iy+1,iz) - 0.5*chi(1,0,ix,iy,iz)/(dy*dx2)*pxp
            coeffs(0,1,-1,ix,iy+1,iz) = 0.
          endif
          if (iz > 0) then
            rhs(ix,iy,iz-1)=rhs(ix,iy,iz-1) + 0.5*chi(2,0,ix,iy,iz)/(dz*dx2)*pxp
            coeffs(2,1,+1,ix,iy,iz-1) = 0.
          endif
          if (iz < nzlocal) then
            rhs(ix,iy,iz+1)=rhs(ix,iy,iz+1) - 0.5*chi(2,0,ix,iy,iz)/(dz*dx2)*pxp
            coeffs(2,0,+1,ix,iy,iz+1) = 0.
          endif
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          coeffs(0,1,0,ix,iy,iz) = 0.
          if (iy > 0) coeffs(0,1,+1,ix,iy-1,iz) = 0.
          if (iy < nylocal) coeffs(0,1,-1,ix,iy+1,iz) = 0.
          if (iz > 0) coeffs(2,1,+1,ix,iy,iz-1) = 0.
          if (iz < nzlocal) coeffs(2,0,+1,ix,iy,iz+1) = 0.
        endif

        if (0. < dels(2,ic) .and. dels(2,ic) < +1.) then
          pym = voltfac*volt(0,ic)
          rhs(ix,iy,iz) = rhs(ix,iy,iz) + coeffs(1,0,0,ix,iy,iz)*pym
          coeffs(1,0,0,ix,iy,iz) = 0.
          if (iz > 0) then
            rhs(ix,iy,iz-1)=rhs(ix,iy,iz-1) - 0.5*chi(2,1,ix,iy,iz)/(dz*dy2)*pym
            coeffs(1,0,+1,ix,iy,iz-1) = 0.
          endif
          if (iz < nzlocal) then
            rhs(ix,iy,iz+1)=rhs(ix,iy,iz+1) + 0.5*chi(2,1,ix,iy,iz)/(dz*dy2)*pym
            coeffs(1,0,-1,ix,iy,iz+1) = 0.
          endif
          if (ix > 0) then
            rhs(ix-1,iy,iz)=rhs(ix-1,iy,iz) - 0.5*chi(0,1,ix,iy,iz)/(dx*dy2)*pym
            coeffs(0,1,-1,ix-1,iy,iz) = 0.
          endif
          if (ix < nxlocal) then
            rhs(ix+1,iy,iz)=rhs(ix+1,iy,iz) + 0.5*chi(0,1,ix,iy,iz)/(dx*dy2)*pym
            coeffs(0,0,-1,ix+1,iy,iz) = 0.
          endif
        elseif (-1. < dels(2,ic) .and. dels(2,ic) <= 0.) then
          coeffs(1,0,0,ix,iy,iz) = 0.
          if (iz > 0) coeffs(1,0,+1,ix,iy,iz-1) = 0.
          if (iz < nzlocal) coeffs(1,0,-1,ix,iy,iz+1) = 0.
          if (ix > 0) coeffs(0,1,-1,ix-1,iy,iz) = 0.
          if (ix < nxlocal) coeffs(0,0,-1,ix+1,iy,iz) = 0.
        endif

        if (0. < dels(3,ic) .and. dels(3,ic) < +1.) then
          pyp = voltfac*volt(1,ic)
          rhs(ix,iy,iz) = rhs(ix,iy,iz) + coeffs(1,1,0,ix,iy,iz)*pyp
          coeffs(1,1,0,ix,iy,iz) = 0.
          if (iz > 0) then
            rhs(ix,iy,iz-1) = rhs(ix,iy,iz-1) + 0.5*chi(2,1,ix,iy,iz)/(dz*dy2)*pyp
            coeffs(1,1,+1,ix,iy,iz-1) = 0.
          endif
          if (iz < nzlocal) then
            rhs(ix,iy,iz+1) = rhs(ix,iy,iz+1) - 0.5*chi(2,1,ix,iy,iz)/(dz*dy2)*pyp
            coeffs(1,1,-1,ix,iy,iz+1) = 0.
          endif
          if (ix > 0) then
            rhs(ix-1,iy,iz) = rhs(ix-1,iy,iz) + 0.5*chi(0,1,ix,iy,iz)/(dx*dy2)*pyp
            coeffs(0,1,+1,ix-1,iy,iz) = 0.
          endif
          if (ix < nxlocal) then
            rhs(ix+1,iy,iz) = rhs(ix+1,iy,iz) - 0.5*chi(0,1,ix,iy,iz)/(dx*dy2)*pyp
            coeffs(0,0,+1,ix+1,iy,iz) = 0.
          endif
        elseif (-1. < dels(3,ic) .and. dels(3,ic) <= 0.) then
          coeffs(1,1,0,ix,iy,iz) = 0.
          if (iz > 0) coeffs(1,1,+1,ix,iy,iz-1) = 0.
          if (iz < nzlocal) coeffs(1,1,-1,ix,iy,iz+1) = 0.
          if (ix > 0) coeffs(0,1,+1,ix-1,iy,iz) = 0.
          if (ix < nxlocal) coeffs(0,0,+1,ix+1,iy,iz) = 0.
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          pzm = voltfac*volt(4,ic)
          rhs(ix,iy,iz) = rhs(ix,iy,iz) + coeffs(2,0,0,ix,iy,iz)*pzm
          coeffs(2,0,0,ix,iy,iz) = 0.
          if (ix > 0) then
            rhs(ix-1,iy,iz) = rhs(ix-1,iy,iz) - 0.5*chi(0,2,ix,iy,iz)/(dx*dz2)*pzm
            coeffs(2,0,+1,ix-1,iy,iz) = 0.
          endif
          if (ix < nxlocal) then
            rhs(ix+1,iy,iz) = rhs(ix+1,iy,iz) + 0.5*chi(0,2,ix,iy,iz)/(dx*dz2)*pzm
            coeffs(2,0,-1,ix+1,iy,iz) = 0.
          endif
          if (iy > 0) then
            rhs(ix,iy-1,iz) = rhs(ix,iy-1,iz) - 0.5*chi(1,2,ix,iy,iz)/(dy*dz2)*pzm
            coeffs(1,1,-1,ix,iy-1,iz) = 0.
          endif
          if (iy < nylocal) then
            rhs(ix,iy+1,iz) = rhs(ix,iy+1,iz) + 0.5*chi(1,2,ix,iy,iz)/(dy*dz2)*pzm
            coeffs(1,0,-1,ix,iy+1,iz) = 0.
          endif
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          coeffs(2,0,0,ix,iy,iz) = 0.
          if (ix > 0) coeffs(2,0,+1,ix-1,iy,iz) = 0.
          if (ix < nxlocal) coeffs(2,0,-1,ix+1,iy,iz) = 0.
          if (iy > 0) coeffs(1,1,-1,ix,iy-1,iz) = 0.
          if (iy < nylocal) coeffs(1,0,-1,ix,iy+1,iz) = 0.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          rhs(ix,iy,iz) = rhs(ix,iy,iz) + coeffs(2,1,0,ix,iy,iz)*pzp
          coeffs(2,1,0,ix,iy,iz) = 0.
          if (ix > 0) then
            rhs(ix-1,iy,iz) = rhs(ix-1,iy,iz) + 0.5*chi(0,2,ix,iy,iz)/(dx*dz2)*pzp
            coeffs(2,1,+1,ix-1,iy,iz) = 0.
          endif
          if (ix < nxlocal) then
            rhs(ix+1,iy,iz) = rhs(ix+1,iy,iz) - 0.5*chi(0,2,ix,iy,iz)/(dx*dz2)*pzp
            coeffs(2,1,-1,ix+1,iy,iz) = 0.
          endif
          if (iy > 0) then
            rhs(ix,iy-1,iz) = rhs(ix,iy-1,iz) + 0.5*chi(1,2,ix,iy,iz)/(dy*dz2)*pzp
            coeffs(1,1,+1,ix,iy-1,iz) = 0.
          endif
          if (iy < nylocal) then
            rhs(ix,iy+1,iz) = rhs(ix,iy+1,iz) - 0.5*chi(1,2,ix,iy,iz)/(dy*dz2)*pzp
            coeffs(1,0,+1,ix,iy+1,iz) = 0.
          endif
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          coeffs(2,1,0,ix,iy,iz) = 0.
          if (ix > 0) coeffs(2,1,+1,ix-1,iy,iz) = 0.
          if (ix < nxlocal) coeffs(2,1,-1,ix+1,iy,iz) = 0.
          if (iy > 0) coeffs(1,1,+1,ix,iy-1,iz) = 0.
          if (iy < nylocal) coeffs(1,0,+1,ix,iy+1,iz) = 0.
        endif

      enddo
!$OMP END DO

      if (lf3dtimesubs) timegetcoefficientssubgridimplicites3d =
     &                  timegetcoefficientssubgridimplicites3d
     &                  + wtime() - substarttime

      end subroutine getcoefficientssubgridimplicites3d
      end module coefficientsmodulees3d
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine mgsolveimplicites3d(iwhich,nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                               dx,dy,dz,phi,rho,
     &                               ns,qomdt,chi0,
     &                               rstar,linbend,bounds,
     &                               xmminlocal,ymminlocal,zmminlocal,zgrid,
     &                               mgparam,mgiters,mgmaxiters,
     &                               mgmaxlevels,mgerror,mgtol,mgverbose,
     &                               downpasses,uppasses,
     &                               lcndbndy,laddconductor,icndbndy,
     &                               gridmode,conductors,fsdecomp)
      use Subtimersf3d
      use ConductorTypemodule
      use Constant
      use coefficientsmodulees3d
      use Decompositionmodule
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal,ns
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: chi0(0:nxlocal,0:nylocal,0:nzlocal,0:ns-1)
      real(kind=8):: dx,dy,dz
      real(kind=8):: rstar(-1:nzlocal+1)
      logical(ISZ):: linbend
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmminlocal,ymminlocal,zmminlocal,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels,mgverbose
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

c Use the multigrid method for solving Poisson's equation on a 3-D Cartesian
c mesh. The fieldsolver allows internal conductors with subgrid scale
c resolution.

      integer(ISZ):: i,ii,k,ix,iy,iz
      integer(ISZ):: delx,dely,delz
      real(kind=8):: rs,x,r
      real(kind=8):: dxsqi,dysqi,dzsqi,rdel
      real(kind=8),pointer:: bongrid(:,:,:,:)
      real(kind=8),pointer:: chi0temp(:,:,:,:)
      real(kind=8),allocatable:: phisave(:,:,:)
      real(kind=8):: bendx(-1:nxlocal+1)
      integer(ISZ):: localbounds(0:5)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

      delx = 1
      dely = 1
      delz = 1
      if (nx > nxlocal) delx = 3
      if (ny > nylocal) dely = 3
      if (nz > nzlocal) delz = 3

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

      bongrid => getbfieldsongrid3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                              dx,dy,dz,
     &                              xmminlocal,ymminlocal,zmminlocal,zgrid,
     &                              localbounds)
c     --- Parallel B.C.s are set below.

c     --- Create copy of chi0 with extra guard cells, and fill it in.
      allocate(chi0temp(-delx:nxlocal+delx,
     &                  -dely:nylocal+dely,
     &                  -delz:nzlocal+delz,0:ns-1))
      chi0temp(0:nxlocal,0:nylocal,0:nzlocal,:) = chi0
      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                               chi0temp,ns,localbounds,.false.,.true.)
c     --- Parallel B.C.s are set below.

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

!$OMP PARALLEL
!$OMP&PRIVATE(ii,i,k,rs,x,r,ix,iy,iz)

      if (linbend) then
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
        dxsqi  = 1./dx**2
        dysqi  = 1./dy**2
        dzsqi  = 1./dz**2
        rdel   = dzsqi/(dxsqi + dysqi + dzsqi)
!$OMP DO
        do iz=0,nzlocal
          rs = rstar(iz)
          do ix=0,nxlocal
            x  = xmminlocal + ix*dx
            r  = rs + x
c           --- rearranged to reduce divides
c           --- rho(ix,:,iz) = rho(ix,:,iz)*(rs/r)/
c           ---             ( 1. + (x/r)*((x/r)-2.)*rdel )
            rho(ix,:,iz) = rho(ix,:,iz)*rs*r/(r*r + x*(x-2.*r)*rdel)
          enddo
        enddo
!$OMP END DO
c       --- Fill scratch array with x values so it can be looked up
c       --- in the bent beam loop instead of calculated.
!$OMP DO
        do ix = -1,nxlocal+1
          bendx(ix) = xmminlocal + ix*dx
        enddo
!$OMP END DO
c       --- Change rstar if using Nuemann boundary conditions
        if (localbounds(4) == 1) rstar(-1) = rstar(1)
        if (localbounds(5) == 1) rstar(nzlocal+1) = rstar(nzlocal-1)
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

#ifdef MPIPARALLEL
      call applyparallelboundaryconditions3d(chi0temp,
     &               nxlocal,nylocal,nzlocal,
     &               delx,dely,delz,ns,3,localbounds,fsdecomp)
      call applyparallelboundaryconditions3d(bongrid,
     &               nxlocal,nylocal,nzlocal,
     &               delx,dely,delz,3,3,localbounds,fsdecomp)
#endif

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

c       --- Do one vcycle.
        call vcycleimplicites3d(0,1.,nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                          dx,dy,dz,phi,rho,
     &                          ns,qomdt,chi0temp,bongrid,
     &                          rstar,linbend,bendx,bounds,mgparam,mgmaxlevels,
     &                          downpasses,uppasses,
     &                          lcndbndy,icndbndy,conductors,
     &                          delx,dely,delz,fsdecomp)

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

#ifdef MPIPARALLEL
        if (fsdecomp%nxprocs*fsdecomp%nyprocs*fsdecomp%nzprocs > 1) then
c         --- calculate global sorerror
          call parallelmaxrealarray(mgerror,1)
c     added by petermc, 26 Sep 2002
        endif
#endif

      enddo

c     --- Set boundary conditions. This is only really needed for the
c     --- Dirichlet boundaries, but this is convenient to call.
      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                               localbounds,.true.,.false.)

c     --- Make a print out.
      if (mgverbose>=1 .or. mgerror > mgtol) then
        if (mgerror > mgtol) then
          call remark("Multigrid: Maximum number of iterations reached")
        endif
        write(errline,20) mgerror,mgiters
  20    format("Multigrid: Error converged to ",1pe11.3," in ",i5," v-cycles")
        call remark(errline)
      endif
      
      deallocate(phisave)
      deallocate(chi0temp)
      deallocate(bongrid)

c     --- Restore rho
      if (linbend) then
c       --- For bends, also include curvature corrections. Comment: Timing tests
c       --- show that the use of 1d array is slightly faster than a 3d array.
        do iz=0,nzlocal
          rs = rstar(iz)
          do ix=0,nxlocal
            x  = xmminlocal + ix*dx
            r  = rs + x
            rho(ix,:,iz) = rho(ix,:,iz)/rs*( r + x*((x/r)-2.)*rdel )
          enddo
        enddo
      endif

!$OMP END PARALLEL

      if (lf3dtimesubs) timemgsolveimplicites3d = timemgsolveimplicites3d +
     &                                            wtime() - substarttime

      return
      end
c=============================================================================
      RECURSIVE subroutine vcycleimplicites3d(mglevel,mgscale,nx,ny,nz,
     &                                   nxlocal,nylocal,nzlocal,
     &                                   dx,dy,dz,
     &                                   phi,rho,ns,qomdt,chi0,bongrid,
     &                                   rstar,linbend,bendx,globalbounds,
     &                                   mgparam,
     &                                   mgmaxlevels,downpasses,uppasses,
     &                                   lcndbndy,icndbndy,conductors,
     &                                   delx,dely,delz,fsdecomp)
      use Subtimersf3d
      use Constant
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use coefficientsmodulees3d
      use formggetarraysuminterface
      use Decompositionmodule
      integer(ISZ):: mglevel
      real(kind=8):: mgscale
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal,ns
      integer(ISZ):: delx,dely,delz
      real(kind=8):: dx,dy,dz
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: chi0(-delx:nxlocal+delx,
     &                    -dely:nylocal+dely,
     &                    -delz:nzlocal+delz,0:ns-1)
      real(kind=8):: bongrid(-delx:nxlocal+delx,
     &                       -dely:nylocal+dely,
     &                       -delz:nzlocal+delz,0:2)
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: rstar(-1:nzlocal+1)
      real(kind=8):: bendx(-1:nxlocal+1)
      logical(ISZ):: linbend
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy
      integer(ISZ):: icndbndy
      type(Decomposition):: fsdecomp

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8),pointer:: coeffs0(:,:,:),coeffs(:,:,:,:,:,:),rhs(:,:,:)
      real(kind=8),pointer:: resscale(:,:,:)
      real(kind=8),allocatable:: phicoarse(:,:,:),rhocoarse(:,:,:)
      real(kind=8),allocatable:: chi0coarse(:,:,:,:),bongridcoarse(:,:,:,:)
      real(kind=8),allocatable:: res(:,:,:)
      integer(ISZ):: i,ic,iszone=1
      integer(ISZ):: i1,i2,i3
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
      real(kind=8):: sss(2),reps0c
      logical(ISZ):: lpe0
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Check that the same coarsening was done for the conductors as is
c     --- being done now. Note that only nzlocal is saved for the conductors so
c     --- only it can be checked.
      if (nzlocal .ne. conductors%levelnz(mglevel)) then
        call kaboom("vcycleimplicites3d: error: the coarsening level of the
     &conductor is inconsistent with that used by the solver")
        return
      endif

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

c     --- Get the finite difference coefficients for the Poisson operator at
c     --- the current level of refinement.
      allocate(coeffs0(0:nxlocal,0:nylocal,0:nzlocal))
      allocate(coeffs(0:2,0:1,-1:1,0:nxlocal,0:nylocal,0:nzlocal))
      allocate(rhs(0:nxlocal,0:nylocal,0:nzlocal))
      allocate(resscale(0:nxlocal,0:nylocal,0:nzlocal))
      call getcoefficientsimplicites3d(nxlocal,nylocal,nzlocal,dx,dy,dz,ns,
     &                                 qomdt,chi0,bongrid,
     &                                 mglevel,lcndbndy,icndbndy,conductors,
     &                                 localbounds,coeffs0,coeffs,rhs,resscale,
     &                                 delx,dely,delz)

      if (lprintmgarraysumdiagnostic) then
#ifdef MPIPARALLEL
        lpe0=(fsdecomp%ixproc==0.and.fsdecomp%iyproc==0.and.fsdecomp%izproc==0)
#else
        lpe0 = .true.
#endif
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
        if (lpe0) print*,"V1 phi",mglevel,sss
        reps0c = mgparam/(eps0*2.*(1./dx**2+1./dy**2+1./dz**2))
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V1 rho",mglevel,sss*reps0c
        do i=0,ns-1
          sss = mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,chi0(:,:,:,i),fsdecomp)
          if (lpe0) print*,"V1 chi0",i,mglevel,sss
        enddo
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,coeffs0,fsdecomp)
        if (lpe0) print*,"V1 coeffs0",mglevel,sss
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,rhs,fsdecomp)
        if (lpe0) print*,"V1 rhs",mglevel,sss
        do i3=-1,1
          do i2=0,1
            do i1=0,2
              sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,coeffs(i1,i2,i3,:,:,:),fsdecomp)
              if (lpe0) print*,"V1 coeffs",i1,i2,i3,mglevel,sss
            enddo
          enddo
        enddo
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,bongrid(:,:,:,0),fsdecomp)
        if (lpe0) print*,"V1 bongrid 0",mglevel,sss
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,bongrid(:,:,:,1),fsdecomp)
        if (lpe0) print*,"V1 bongrid 1",mglevel,sss
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,bongrid(:,:,:,2),fsdecomp)
        if (lpe0) print*,"V1 bongrid 2",mglevel,sss
      endif

      do i=1,downpasses
        call relaximplicites3d(mglevel,nx,ny,nz,nxlocal,nylocal,nzlocal,phi,rho,
     &                         coeffs0,coeffs,rhs,
     &                         localbounds,mgparam,conductors,fsdecomp)
      enddo

      if (lprintmgarraysumdiagnostic) then
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
        if (lpe0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,nylocal,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V2 rho",mglevel,sss*reps0c
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. ny >= 4 .and. nz >= 4 .and.
     &    mglevel < mgmaxlevels) then

        allocate(res(-delx:nxlocal+delx,
     &               -dely:nylocal+dely,
     &               -delz:nzlocal+delz),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleimplicites3d: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nxlocal,nylocal,nzlocal
          call kaboom("vcycleimplicites3d: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residualimplicites3d(nxlocal,nylocal,nzlocal,phi,rho,
     &                            coeffs0,coeffs,rhs,resscale,res,
     &                            mglevel,localbounds,mgparam,conductors,
     &                            delx,dely,delz)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,res,
     &                               delx,dely,delz,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,nylocal,nzlocal,res,
     &                      delx,dely,delz,-3,-1,fsdecomp)
#endif
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocal,nylocal,nzlocal,delx,dely,delz,res,fsdecomp)
          if (lpe0) print*,"V3 res",mglevel,sss*reps0c
        endif

        call getnextcoarselevel3d(nx,ny,nz,nxlocal,nylocal,nzlocal,dx,dy,dz,
     &                            nxcoarse,nycoarse,nzcoarse,
     &                            dxcoarse,dycoarse,dzcoarse)

        dxcoarsesqi = 1./dxcoarse**2
        dycoarsesqi = 1./dycoarse**2
        dzcoarsesqi = 1./dzcoarse**2
        mgscalecoarse = mgscale*dxcoarse*dycoarse*dzcoarse/(dx*dy*dz)

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
        allocate(phicoarse(-1:nxlocalcoarse+1,-1:nylocalcoarse+1,-1:nzlocalcoarse+1),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleimplicites3d: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxlocalcoarse,nylocalcoarse,nzlocalcoarse
          call kaboom("vcycleimplicites3d: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxlocalcoarse,0:nylocalcoarse,0:nzlocalcoarse),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleimplicites3d: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxlocalcoarse,nylocalcoarse,nzlocalcoarse
          call kaboom("vcycleimplicites3d: allocation error")
          return
        endif
        allocate(chi0coarse(-delx:nxlocalcoarse+delx,
     &                      -dely:nylocalcoarse+dely,
     &                      -delz:nzlocalcoarse+delz,0:ns-1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleimplicites3d: allocation error ",allocerror,
     &           ": could not allocate chi0coarse to shape ",
     &           nxlocalcoarse,nylocalcoarse,nzlocalcoarse
          call kaboom("vcycleimplicites3d: allocation error")
          return
        endif
        allocate(bongridcoarse(-delx:nxlocalcoarse+delx,
     &                         -dely:nylocalcoarse+dely,
     &                         -delz:nzlocalcoarse+delz,0:2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcycleimplicites3d: allocation error ",allocerror,
     &           ": could not allocate bongridcoarse to shape ",
     &           nxlocalcoarse,nylocalcoarse,nzlocalcoarse,delz
          call kaboom("vcycleimplicites3d: allocation error")
          return
        endif

        rhocoarse = 0. ! not needed since all elements set in restrict3d
        phicoarse = 0.
        chi0coarse = 0.
        bongridcoarse = 0.

c       --- Restriction
        call restrict3d(nx,ny,nz,nxlocal,nylocal,nzlocal,res,delx,dely,delz,
     &                  nxcoarse,nycoarse,nzcoarse,
     &                  nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                  rhocoarse,0,0,0,1.,localbounds,localboundsc,
     &                  lxoffset,lyoffset,lzoffset)
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,nylocalcoarse,nzlocalcoarse,0,0,0,rhocoarse,coarsedecomp)
          reps0c = mgparam/(eps0*2.*(1./dxcoarse**2+1./dycoarse**2+1./dzcoarse**2))
          if (lpe0) print*,"V3 rhocoarse",mglevel,sss*reps0c
        endif

        do ic=0,ns-1
          call restrict3d(nx,ny,nz,nxlocal,nylocal,nzlocal,chi0(:,:,:,ic),
     &                    delx,dely,delz,
     &                    nxcoarse,nycoarse,nzcoarse,
     &                    nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                    chi0coarse(:,:,:,ic),delx,dely,delz,
     &                    1.,localbounds,localboundsc,
     &                    lxoffset,lyoffset,lzoffset)
        enddo
        do ic=0,2
          call restrict3d(nx,ny,nz,nxlocal,nylocal,nzlocal,bongrid(:,:,:,ic),
     &                    delx,dely,delz,
     &                    nxcoarse,nycoarse,nzcoarse,
     &                    nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                    bongridcoarse(:,:,:,ic),delx,1,delz,
     &                    1.,localbounds,localboundsc,
     &                    lxoffset,lyoffset,lzoffset)
        enddo
        call applyboundaryconditions3d(nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                                 delx,dely,delz,chi0coarse,ns,
     &                                 localboundsc,.false.,.true.)
        call applyboundaryconditions3d(nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                                 delx,dely,delz,bongridcoarse,3,
     &                                 localboundsc,.false.,.true.)
#ifdef MPIPARALLEL
        call applyparallelboundaryconditions3d(chi0coarse,
     &               nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &               delx,dely,delz,ns,3,localboundsc,coarsedecomp)
        call applyparallelboundaryconditions3d(bongridcoarse,
     &               nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &               delx,dely,delz,3,3,localboundsc,coarsedecomp)
#endif

c       --- Continue at the next coarsest level.
        call vcycleimplicites3d(mglevel+iszone,mgscalecoarse,
     &              nxcoarse,nycoarse,nzcoarse,
     &              nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &              dxcoarse,dycoarse,dzcoarse,phicoarse,rhocoarse,
     &              ns,qomdt,chi0coarse,bongridcoarse,
     &              rstar,linbend,bendx,globalbounds,mgparam,
     &              mgmaxlevels,downpasses,uppasses,
     &              lcndbndy,icndbndy,conductors,
     &              delx,dely,delz,coarsedecomp)

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,nylocalcoarse,nzlocalcoarse,1,1,1,phicoarse,coarsedecomp)
          if (lpe0) print*,"V4 phicoarse",mglevel,sss
          sss = mggetarraysum(nxlocal,nylocal,nzlocal,1,1,1,phi,fsdecomp)
          if (lpe0) print*,"V4 phi",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand3d(nx,ny,nz,nxlocal,nylocal,nzlocal,1,1,1,phi,
     &                nxcoarse,nycoarse,nzcoarse,
     &                nxlocalcoarse,nylocalcoarse,nzlocalcoarse,
     &                phicoarse,localbounds,lxoffset,lyoffset,lzoffset,
     &                conductors,.false.)
        call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                                 localbounds,.false.,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,phi,
     &                               1,1,1,-1,-1,localbounds,fsdecomp)
#endif

        deallocate(res)
        deallocate(phicoarse,rhocoarse)
        deallocate(chi0coarse,bongridcoarse)

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
      do i=1,uppasses
        call relaximplicites3d(mglevel,nx,ny,nz,nxlocal,nylocal,nzlocal,phi,rho,
     &                         coeffs0,coeffs,rhs,
     &                         localbounds,mgparam,conductors,fsdecomp)
      enddo

      deallocate(coeffs0,coeffs)
      deallocate(rhs,resscale)

      if (lf3dtimesubs) timevcycleimplicites3d = timevcycleimplicites3d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
      subroutine relaximplicites3d(mglevel,nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                             phi,rho,
     &                             coeffs0,coeffs,rhs,
     &                             localbounds,mgparam,conductors,fsdecomp)
      use Subtimersf3d
      use Constant
      use ConductorTypemodule
      use Decompositionmodule
      integer(ISZ):: mglevel,nx,ny,nz,nxlocal,nylocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: coeffs0(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:1,-1:1,0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nylocal,0:nzlocal)
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: mgparam
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.
c phisave is needed because when there are B fields, the even and odd sweeps
c are no longer independent since the diagonal coefficients become non-zero.
c The simplest solution is to save the values of phi before a sweep and
c read those saved values for the update of phi.

      real(kind=8),allocatable:: phisave(:,:,:)
      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax,ix,iy,iz,ix1
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      allocate(phisave(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1))

      call cond_potmg(conductors%interior,
     &                nxlocal,nylocal,nzlocal,1,1,1,phi,mglevel,1,.false.)
      call condbndymgint(conductors%evensubgrid,nxlocal,nylocal,nzlocal,
     &                   1,1,1,phi,localbounds,mglevel,2)
      call condbndymgint(conductors%oddsubgrid,nxlocal,nylocal,nzlocal,
     &                   1,1,1,phi,localbounds,mglevel,2)

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

c     --- Set min and max indices for phi array.
c     --- If using Dirichlet boundary conditions, do not solve for the
c     --- potential on the grid edge.
      ixmin = 0
      ixmax = nxlocal
      iymin = 0
      iymax = nylocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) < 1) ixmin = 1
      if (localbounds(1) < 1) ixmax = nxlocal - 1
      if (localbounds(2) < 1) iymin = 1
      if (localbounds(3) < 1) iymax = nylocal - 1
      if (localbounds(4) < 1) izmin = 1
      if (localbounds(5) < 1) izmax = nzlocal - 1

c     --- do loop to cover even and odd points
      do parity=s_parity,e_parity,e_parity-s_parity
        phisave = phi

        do iz=izmin,izmax
          do iy=iymin,iymax
            ix1 = ixmin + mod(ixmin + iy + iz + parity,2)
            do ix=ix1,ixmax,2

              phi(ix,iy,iz) = mgparam*(rho(ix,iy,iz)/eps0 + rhs(ix,iy,iz)

     &                       + coeffs(1,0,-1,ix,iy,iz)*phisave(ix  ,iy-1,iz-1)
     &                       + coeffs(2,0,-1,ix,iy,iz)*phisave(ix-1,iy  ,iz-1)
     &                       + coeffs(2,0, 0,ix,iy,iz)*phisave(ix  ,iy  ,iz-1)
     &                       + coeffs(2,0,+1,ix,iy,iz)*phisave(ix+1,iy  ,iz-1)
     &                       + coeffs(1,1,-1,ix,iy,iz)*phisave(ix  ,iy+1,iz-1)

     &                       + coeffs(0,0,-1,ix,iy,iz)*phisave(ix-1,iy-1,iz  )
     &                       + coeffs(1,0, 0,ix,iy,iz)*phisave(ix  ,iy-1,iz  )
     &                       + coeffs(0,1,-1,ix,iy,iz)*phisave(ix+1,iy-1,iz  )
     &                       + coeffs(0,0, 0,ix,iy,iz)*phisave(ix-1,iy  ,iz  )
     &                       + coeffs(0,1, 0,ix,iy,iz)*phisave(ix+1,iy  ,iz  )
     &                       + coeffs(0,0,+1,ix,iy,iz)*phisave(ix-1,iy+1,iz  )
     &                       + coeffs(1,1, 0,ix,iy,iz)*phisave(ix  ,iy+1,iz  )
     &                       + coeffs(0,1,+1,ix,iy,iz)*phisave(ix+1,iy+1,iz  )

     &                       + coeffs(1,0,+1,ix,iy,iz)*phisave(ix  ,iy-1,iz+1)
     &                       + coeffs(2,1,-1,ix,iy,iz)*phisave(ix-1,iy  ,iz+1)
     &                       + coeffs(2,1, 0,ix,iy,iz)*phisave(ix  ,iy  ,iz+1)
     &                       + coeffs(2,1,+1,ix,iy,iz)*phisave(ix+1,iy  ,iz+1)
     &                       + coeffs(1,1,+1,ix,iy,iz)*phisave(ix  ,iy+1,iz+1)

     &                        )/(-coeffs0(ix,iy,iz))
     &                      + (1.-mgparam)*phi(ix,iy,iz)

            enddo
          enddo
        enddo

        call cond_potmg(conductors%interior,
     &                  nxlocal,nylocal,nzlocal,1,1,1,phi,mglevel,1,.false.)

c       --- set phi in the z guard planes
        call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,1,1,1,phi,1,
     &                                 localbounds,.false.,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,phi,
     &                               1,1,1,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,0,fsdecomp)
#endif

c     --- end of loop over even and odd points
      enddo

#ifdef MPIPARALLEL
c     --- Exchange phi in the z guard planes
      call mgexchange_phi(1,nxlocal,nylocal,nzlocal,phi,1,1,1,-1,-1,fsdecomp)
#endif

      call averageperiodicphi3d(phi,nx,ny,nz,nxlocal,nylocal,nzlocal,1,1,1,
     &                          localbounds)

      deallocate(phisave)

      if (lf3dtimesubs) timerelaximplicites3d = timerelaximplicites3d +
     &                                        wtime() - substarttime

      return
      end
c=============================================================================
      subroutine residualimplicites3d(nxlocal,nylocal,nzlocal,phi,rho,
     &                                coeffs0,coeffs,rhs,resscale,res,
     &                                mglevel,localbounds,mgparam,conductors,
     &                                delx,dely,delz)
      use Subtimersf3d
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nxlocal,nylocal,nzlocal,delx,dely,delz
      real(kind=8):: phi(-1:nxlocal+1,-1:nylocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: coeffs0(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:1,-1:1,0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: resscale(0:nxlocal,0:nylocal,0:nzlocal)
      real(kind=8):: res(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      integer(ISZ):: mglevel,localbounds(0:5)
      real(kind=8):: mgparam
      type(ConductorType):: conductors

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.

      integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax,ix,iy,iz,ix1
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

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

      res = 0.
c     --- Calculate the residual.
!$OMP DO
      do iz=izmin,izmax
        do iy=iymin,iymax
          do ix=ixmin,ixmax

            res(ix,iy,iz) = (rho(ix,iy,iz) + eps0*(rhs(ix,iy,iz)

     &                       + coeffs(1,0,-1,ix,iy,iz)*phi(ix  ,iy-1,iz-1)
     &                       + coeffs(2,0,-1,ix,iy,iz)*phi(ix-1,iy  ,iz-1)
     &                       + coeffs(2,0, 0,ix,iy,iz)*phi(ix  ,iy  ,iz-1)
     &                       + coeffs(2,0,+1,ix,iy,iz)*phi(ix+1,iy  ,iz-1)
     &                       + coeffs(1,1,-1,ix,iy,iz)*phi(ix  ,iy+1,iz-1)

     &                       + coeffs(0,0,-1,ix,iy,iz)*phi(ix-1,iy-1,iz  )
     &                       + coeffs(1,0, 0,ix,iy,iz)*phi(ix  ,iy-1,iz  )
     &                       + coeffs(0,1,-1,ix,iy,iz)*phi(ix+1,iy-1,iz  )
     &                       + coeffs(0,0, 0,ix,iy,iz)*phi(ix-1,iy  ,iz  )
     &                       + coeffs0(ix,iy,iz)*phi(ix  ,iy  ,iz  )
     &                       + coeffs(0,1, 0,ix,iy,iz)*phi(ix+1,iy  ,iz  )
     &                       + coeffs(0,0,+1,ix,iy,iz)*phi(ix-1,iy+1,iz  )
     &                       + coeffs(1,1, 0,ix,iy,iz)*phi(ix  ,iy+1,iz  )
     &                       + coeffs(0,1,+1,ix,iy,iz)*phi(ix+1,iy+1,iz  )

     &                       + coeffs(1,0,+1,ix,iy,iz)*phi(ix  ,iy-1,iz+1)
     &                       + coeffs(2,1,-1,ix,iy,iz)*phi(ix-1,iy  ,iz+1)
     &                       + coeffs(2,1, 0,ix,iy,iz)*phi(ix  ,iy  ,iz+1)
     &                       + coeffs(2,1,+1,ix,iy,iz)*phi(ix+1,iy  ,iz+1)
     &                       + coeffs(1,1,+1,ix,iy,iz)*phi(ix  ,iy+1,iz+1)
     &                    ))*resscale(ix,iy,iz)

          enddo
        enddo
      enddo
!$OMP END DO

      call cond_potmgres(conductors%interior,
     &                   nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                   res,mglevel,1,.false.)
      call applyboundaryconditions3d(nxlocal,nylocal,nzlocal,delx,dely,delz,
     &                               res,1,localbounds,.true.,.false.)

      if (lf3dtimesubs) timeresidualimplicites3d = timeresidualimplicites3d +
     &                                           wtime() - substarttime

      return
      end
c=============================================================================
c=============================================================================
      subroutine averageperiodicphi3d(phi,nx,ny,nz,nxlocal,nylocal,nzlocal,
     &                                delx,dely,delz,localbounds)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,nxlocal,nylocal,nzlocal,delx,dely,delz
      real(kind=8):: phi(-delx:nxlocal+delx,
     &                   -dely:nylocal+dely,
     &                   -delz:nzlocal+delz)
      integer(ISZ):: localbounds(0:5)

c Subtracts out the average phi.
c Warning - does not work in parallel.

      real(kind=8):: phiave
      integer(ISZ):: ncells
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (localbounds(0) == 2 .and. localbounds(1) == 2 .and.
     &    localbounds(2) == 2 .and. localbounds(3) == 2 .and.
     &    localbounds(4) == 2 .and. localbounds(5) == 2) then

        ncells = 1
        if (nx > 0) ncells = ncells*nx
        if (ny > 0) ncells = ncells*ny
        if (nz > 0) ncells = ncells*nz

        phiave = sum(phi(0:nxlocal-1,0:nylocal-1,0:nzlocal-1))/ncells
c       --- This needs a parallel sum (handling overlapping boundaries
c       --- appropriately)
        phi = phi - phiave

      endif

      if (lf3dtimesubs) timeaverageperiodicphi3d = timeaverageperiodicphi3d +
     &                                             wtime() - substarttime

      return
      end
c=============================================================================
c=============================================================================
      subroutine applyparallelboundaryconditions3d(u,nxlocal,nylocal,nzlocal,
     &                 delx,dely,delz,ncomp,zsend,localbounds,fsdecomp)
      use Subtimersf3d
      use Decompositionmodule
      integer(ISZ):: nxlocal,nylocal,nzlocal,delx,dely,delz,ncomp,zsend
      real(kind=8):: u(-delx:nxlocal+delx,
     &                 -dely:nylocal+dely,
     &                 -delz:nzlocal+delz,ncomp)
      integer(ISZ):: localbounds(0:5)
      type(Decomposition):: fsdecomp

      integer(ISZ):: iz,ic
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

#ifdef MPIPARALLEL

      do ic=1,ncomp
        call mgexchange_phi_periodic(1,nxlocal,nylocal,nzlocal,u(:,:,:,ic),
     &                               delx,dely,delz,-1,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,nylocal,nzlocal,u(:,:,:,ic),
     &                      delx,dely,delz,-zsend,-1,fsdecomp)
      enddo

#endif

      if (lf3dtimesubs) timeapplyparallelboundaryconditions3d =
     &                  timeapplyparallelboundaryconditions3d +
     &                                         wtime() - substarttime

      return
      end
c=============================================================================
