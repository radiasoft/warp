#include "top.h"
c=============================================================================
c@(#) File frz_ImplicitES.F, version $Revision: 1.33 $, $Date: 2008/11/19 18:29:49 $
c# Copyright (c) 2007-2007, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the Poisson solver for the electrostatic implicit scheme.
c  David P. Grote, LLNL, (925)423-7194 or (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      module coefficientsmodulees2d
      contains
      function getbfieldsongrides2d(nxlocal,nzlocal,delx,delz,dx,dz,
     &                              xmminlocal,zmminlocal,zgrid,localbounds)
     &                              result(bongrid)
      real(kind=8),pointer:: bongrid(:,:,:)
      integer(ISZ):: nxlocal,nzlocal,delx,delz
      real(kind=8):: dx,dz,xmminlocal,zmminlocal,zgrid
      integer(ISZ):: localbounds(0:5)

      real(kind=8),allocatable:: xx(:),yy(:),zz(:),uz(:),gaminv(:)
      real(kind=8),allocatable:: ex(:),ey(:),ez(:),bendres(:),bendradi(:)
      integer(ISZ)::ix,iz

c     --- Allocate space for B field
      allocate(bongrid(-delx:nxlocal+delx,-delz:nzlocal+delz,0:2))

c     --- Create temporary space to hold the coordinates of the grid
c     --- needed to look up the B fields.
      allocate(xx(0:nxlocal),yy(0:nxlocal),zz(0:nxlocal))
      allocate(uz(0:nxlocal),gaminv(0:nxlocal))
      allocate(ex(0:nxlocal),ey(0:nxlocal),ez(0:nxlocal))
      allocate(bendres(0:nxlocal),bendradi(0:nxlocal))

c     --- Fill the transverse position arrays.
      do ix=0,nxlocal
        xx(ix) = xmminlocal + ix*dx
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
        zz = zmminlocal + iz*dz + zgrid
        call exteb3d(1+nxlocal,xx,yy,zz,uz,gaminv,0.,0.,
     &               bongrid(0:nxlocal,iz,0),
     &               bongrid(0:nxlocal,iz,1),
     &               bongrid(0:nxlocal,iz,2),
     &               ex,ey,ez,1.,1.,bendres,bendradi,1.,1.)
      enddo

c     --- Apply the boundary conditions, setting the values in the guard cells.
      call applyboundaryconditions3d(nxlocal,0,nzlocal,delx,0,delz,bongrid,3,
     &                               localbounds,.false.,.true.)

      deallocate(xx,yy,zz,uz,gaminv)
      deallocate(ex,ey,ez)
      deallocate(bendres,bendradi)

      return
      end function getbfieldsongrides2d
c=============================================================================
      subroutine getcoefficientses2d(nxlocal,nzlocal,delx,delz,dx,dz,
     &                               xmminlocal,ns,qomdt,rho,chi0,bongrid,
     &                               mglevel,lcndbndy,icndbndy,conductors,lrz,
     &                               localbounds,coeffs,rhs,fsdecomp,resdelfac)
      use Constant
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use DKInterp
      use Decompositionmodule
      use formggetarraysuminterface
      integer(ISZ):: nxlocal,nzlocal,delx,delz,ns
      real(kind=8):: dx,dz,xmminlocal
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: rho(0:nxlocal,0:nzlocal)
      real(kind=8):: chi0(-delx:nxlocal+delx,-delz:nzlocal+delz,0:ns-1)
      real(kind=8):: bongrid(-delx:nxlocal+delx,-delz:nzlocal+delz,0:2)
      integer(ISZ):: mglevel
      integer(ISZ):: icndbndy
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lrz
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: coeffs(0:2,0:2,0:nxlocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nzlocal)
      real(kind=8):: resdelfac(0:nxlocal,0:nzlocal)
      type(Decomposition):: fsdecomp

      real(kind=8),allocatable:: chi(:,:,:,:)
      integer(ISZ):: ix,iz,js,ic,js1
      integer(ISZ):: ixmin,ixmax,izmin,izmax
      real(kind=8):: ox,oy,oz,oo
      real(kind=8):: rr,rrm,rrp
      real(kind=8):: alpha,alphabar,oneminus
      real(kind=8):: osq,txx,txz,tzx,tzz,halfosqi,omcdti,bsqi,by
      real(kind=8):: bbxx,bbxz,bbzx,bbzz,iperpxx,iperpxz,iperpzx,iperpzz
      real(kind=8):: bcrossixz,bcrossizx
      real(kind=8):: sss(2)
      logical(ISZ):: lpe0

      oneminus = 1.0 - 1.e-12

      allocate(chi(0:1,0:1,-1:nxlocal+1,-1:nzlocal+1))
      chi = 0.

c     --- First, generate the chi tensor from chi0 and the B field.
      do js=0,ns-1
        js1=js+1
        do iz=-1,nzlocal+1
          do ix=-1,nxlocal+1

            ox = 0.5*qomdt(js)*bongrid(ix,iz,0)
            oy = 0.5*qomdt(js)*bongrid(ix,iz,1)
            oz = 0.5*qomdt(js)*bongrid(ix,iz,2)
c           -- (omega_c dt/2)**2 ... 
            osq = ox**2+oy**2+oz**2
            oo = 1./(1. + osq)
c           Birdsall-Langdon T matrix components..
            txx = oo*(1.+ox**2)
            txz = oo*(ox*oz-oy)
            tzx = oo*(ox*oz+oy)
            tzz = oo*(1.+oz**2)

            if (interpdk(js1) == 0 .or. oo > oneminus) then
c           -- Not interpolating this species
c           -- Reminder, first two indices of chi are the x-z tensor indices
c           -- Now fixed for standard tensoral notation,i.e. used as
c           -- Sum_ij (d/dx_i) chi_ij dphi/dx_j
               chi(0,0,ix,iz) = chi(0,0,ix,iz) + chi0(ix,iz,js)*txx
               chi(0,1,ix,iz) = chi(0,1,ix,iz) + chi0(ix,iz,js)*txz
               chi(1,0,ix,iz) = chi(1,0,ix,iz) + chi0(ix,iz,js)*tzx
               chi(1,1,ix,iz) = chi(1,1,ix,iz) + chi0(ix,iz,js)*tzz
            else
c           -- Interpolated species
               alpha = sqrt(oo)
               alphabar = 1.-alpha
c           -- store 2/(omega_c dt)**2 and 1/(omega_c dt)
               halfosqi=0.5/osq
               omcdti = sqrt(0.25/osq)
c           -- component of unit vector in direction of B
               bsqi = 1./(bongrid(ix,iz,0)**2+bongrid(ix,iz,1)**2 +
     &           bongrid(ix,iz,2)**2)
               by = bongrid(ix,iz,1)*sqrt(bsqi)
c           -- components of bb diad
               bbxx = bongrid(ix,iz,0)**2*bsqi
               bbxz = bongrid(ix,iz,0)*bongrid(ix,iz,2)*bsqi
               bbzx = bbxz
               bbzz = bongrid(ix,iz,2)**2*bsqi
c           -- components of identity tensor perpendicular to B
               iperpxx = 1.-bbxx
               iperpxz = -bbxz
               iperpzx = iperpxz
               iperpzz = 1.-bbzz
c           -- components of b cross I
c           --  note diagonal terms bcrossixx = bcrossizz = 0
               bcrossixz = by
               bcrossizx = - by

               chi(0,0,ix,iz) = chi(0,0,ix,iz) + chi0(ix,iz,js)*
     &           (alpha*txx+alphabar*(bbxx +alphabar*halfosqi*iperpxx))
               chi(0,1,ix,iz) = chi(0,1,ix,iz) + chi0(ix,iz,js)*
     &           (alpha*txz+alphabar*(bbxz +alphabar*halfosqi*iperpxz -
     &            omcdti*bcrossixz))
               chi(1,0,ix,iz) = chi(1,0,ix,iz) + chi0(ix,iz,js)*
     &           (alpha*tzx+alphabar*(bbzx +alphabar*halfosqi*iperpzx -
     &            omcdti*bcrossizx))
               chi(1,1,ix,iz) = chi(1,1,ix,iz) + chi0(ix,iz,js)*
     &           (alpha*tzz+alphabar*(bbzz +alphabar*halfosqi*iperpzz))
            endif
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

      if (lprintmgarraysumdiagnostic) then
        lpe0=(fsdecomp%ixproc==0.and.fsdecomp%izproc==0)
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,chi(0,0,:,:),fsdecomp)
        if (lpe0) print*,"V1 chi00",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,chi(1,0,:,:),fsdecomp)
        if (lpe0) print*,"V1 chi10",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,chi(0,1,:,:),fsdecomp)
        if (lpe0) print*,"V1 chi01",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,chi(1,1,:,:),fsdecomp)
        if (lpe0) print*,"V1 chi11",mglevel,sss
      endif

      if (lrz) then
c       --- RZ Axisymmetry

c       --- If xmminlocal == 0, then the case of r=0 is included and must be
c       --- treated specially.
        ixmin = 0
        if (xmminlocal == 0.) ixmin = 1

        do iz=0,nzlocal
          if (xmminlocal == 0.) then
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
          do ix=ixmin,nxlocal
            rr = xmminlocal + ix*dx
            rrm = (rr - 0.5*dx)/rr
            rrp = (rr + 0.5*dx)/rr

            coeffs(1,1,ix,iz) = -rrp*(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
     &                          -rrm*(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
     &                          -    (1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
     &                          -    (1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,1,ix,iz) = +rrp*(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
            coeffs(0,1,ix,iz) = +rrm*(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
            coeffs(1,2,ix,iz) = +0.5*(rrp-rrm)*(chi(0,1,ix  ,iz  ))/(2.*dz*dx)
     &                          +    (1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
            coeffs(1,0,ix,iz) = -0.5*(rrp-rrm)*(chi(0,1,ix  ,iz  ))/(2.*dz*dx)
     &                          +    (1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,2,ix,iz) = +0.5*rrp*chi(0,1,ix+1,iz  )/(2.*dz*dx)
     &                          +0.5*    chi(1,0,ix  ,iz+1)/(2.*dx*dz)
            coeffs(0,2,ix,iz) = -0.5*rrm*chi(0,1,ix-1,iz  )/(2.*dz*dx)
     &                          -0.5*    chi(1,0,ix  ,iz+1)/(2.*dx*dz)
            coeffs(2,0,ix,iz) = -0.5*rrp*chi(0,1,ix+1,iz  )/(2.*dz*dx)
     &                          -0.5*    chi(1,0,ix  ,iz-1)/(2.*dx*dz)
            coeffs(0,0,ix,iz) = +0.5*rrm*chi(0,1,ix-1,iz  )/(2.*dz*dx)
     &                          +0.5*    chi(1,0,ix  ,iz-1)/(2.*dx*dz)

          enddo
        enddo

      else
c       --- XZ Planar

        do iz=0,nzlocal
          do ix=0,nxlocal

            coeffs(1,1,ix,iz) = -(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
     &                          -(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
     &                          -(1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
     &                          -(1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,1,ix,iz) = +(1. + 0.5*(chi(0,0,ix  ,iz  )+chi(0,0,ix+1,iz  )))/dx**2
            coeffs(0,1,ix,iz) = +(1. + 0.5*(chi(0,0,ix-1,iz  )+chi(0,0,ix  ,iz  )))/dx**2
            coeffs(1,2,ix,iz) = +(1. + 0.5*(chi(1,1,ix  ,iz  )+chi(1,1,ix  ,iz+1)))/dz**2
            coeffs(1,0,ix,iz) = +(1. + 0.5*(chi(1,1,ix  ,iz-1)+chi(1,1,ix  ,iz  )))/dz**2
            coeffs(2,2,ix,iz) = +0.5*chi(0,1,ix+1,iz  )/(2.*dz*dx)
     &                          +0.5*chi(1,0,ix  ,iz+1)/(2.*dx*dz)
            coeffs(0,2,ix,iz) = -0.5*chi(0,1,ix-1,iz  )/(2.*dz*dx)
     &                          -0.5*chi(1,0,ix  ,iz+1)/(2.*dx*dz)
            coeffs(2,0,ix,iz) = -0.5*chi(0,1,ix+1,iz  )/(2.*dz*dx)
     &                          -0.5*chi(1,0,ix  ,iz-1)/(2.*dx*dz)
            coeffs(0,0,ix,iz) = +0.5*chi(0,1,ix-1,iz  )/(2.*dz*dx)
     &                          +0.5*chi(1,0,ix  ,iz-1)/(2.*dx*dz)

          enddo
        enddo

      endif

      rhs = 0.
      resdelfac = 1.
      call getcoefficientssubgrides2d(conductors%evensubgrid,nxlocal,nzlocal,
     &                                dx,dz,xmminlocal,mglevel,
     &                                lcndbndy,icndbndy,lrz,chi,
     &                                localbounds,coeffs,rhs,resdelfac)
      call getcoefficientssubgrides2d(conductors%oddsubgrid,nxlocal,nzlocal,
     &                                dx,dz,xmminlocal,mglevel,
     &                                lcndbndy,icndbndy,lrz,chi,
     &                                localbounds,coeffs,rhs,resdelfac)

      deallocate(chi)

      return
      end subroutine getcoefficientses2d
c=============================================================================
      subroutine getcoefficientssubgrides2d(subgrid,nxlocal,nzlocal,dx,dz,
     &                                      xmminlocal,
     &                                      mglevel,lcndbndy,icndbndy,lrz,chi,
     &                                      localbounds,coeffs,rhs,resdelfac)
      use Constant
      use ConductorSubgridTypemodule
      type(ConductorSubgridType):: subgrid
      integer(ISZ):: nxlocal,nzlocal
      real(kind=8):: dx,dz,xmminlocal
      integer(ISZ):: mglevel
      integer(ISZ):: icndbndy
      logical(ISZ):: lcndbndy,lrz
      real(kind=8):: chi(0:1,0:1,-1:nxlocal+1,-1:nzlocal+1)
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: coeffs(0:2,0:2,0:nxlocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nzlocal)
      real(kind=8):: resdelfac(0:nxlocal,0:nzlocal)

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
      real(kind=8):: voltfac,c0,ppp
      real(kind=8):: rr(0:nxlocal),rri(0:nxlocal)
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
        do ix=0,nxlocal
          rr(ix) = ix*dx + xmminlocal
          if (rr(ix) > 0.) rri(ix) = 1./rr(ix)
        enddo
      endif

c     --- Loop over points near surface of conductors.
!$OMP DO
      do ic = subgrid%istart(mglevel),subgrid%istart(mglevel+1)-1

        ix = subgrid%indx(0,ic)
        iz = subgrid%indx(2,ic)

c       --- Skip the data point if it is on a Dirichlet or parallel boundary
        if (ix == 0  .and. localbounds(0) < 1) cycle
        if (ix == nxlocal .and. localbounds(1) < 1) cycle
        if (iz == 0  .and. localbounds(4) < 1) cycle
        if (iz == nzlocal .and. localbounds(5) < 1) cycle

c       --- First, get the effective grid cell sizes and the chi's.
c       --- The chi's inside the conductor are always zero because the
c       --- charge density is zero there.
        ppp = 1.
        if (0. < dels(0,ic) .and. dels(0,ic) < +1.) then
          dxm = dels(0,ic)
          chixm = 0.
          ppp = min(ppp,dels(0,ic))
        elseif (-1. < dels(0,ic) .and. dels(0,ic) <= 0.) then
          dxm = -2.*dels(0,ic)
          chixm = 0.
          if (abs(dels(0,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(0,ic)))
          endif
        else
          dxm = 1.
          chixm = chi(:,:,ix-1,iz)
        endif

        if (0. < dels(1,ic) .and. dels(1,ic) < +1.) then
          dxp = dels(1,ic)
          chixp = 0.
          ppp = min(ppp,dels(1,ic))
        elseif (-1. < dels(1,ic) .and. dels(1,ic) <= 0.) then
          dxp = -2.*dels(1,ic)
          chixp = 0.
          if (abs(dels(1,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(1,ic)))
          endif
        else
          dxp = 1.
          chixp = chi(:,:,ix+1,iz)
        endif

        if (0. < dels(4,ic) .and. dels(4,ic) < +1.) then
          dzm = dels(4,ic)
          chizm = 0.
          ppp = min(ppp,dels(4,ic))
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          dzm = -2.*dels(4,ic)
          chizm = 0.
          if (abs(dels(4,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(4,ic)))
          endif
        else
          dzm = 1.
          chizm = chi(:,:,ix,iz-1)
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          dzp = dels(5,ic)
          chizp = 0.
          ppp = min(ppp,dels(5,ic))
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          dzp = -2.*dels(5,ic)
          chizp = 0.
          if (abs(dels(5,ic)) == 0.) then
            ppp = min(ppp,1.-1.e-9)
          else
            ppp = min(ppp,abs(dels(5,ic)))
          endif
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
            rrp = (rr(ix) + 0.5*dxp*dx)*rri(ix)
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
        resdelfac(ix,iz) = ppp

c       --- Now the right hand side can be generated and the associated
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
          if (lrz) then
            if (ix > 1 .or. (ix == 1 .and. xmminlocal > 0.)) then
              rhs(ix-1,iz) = rhs(ix-1,iz)
     &           - 0.5*(rr(ix)-0.5*dx)/(rr(ix)-dx)*chi(1,0,ix,iz)/(dx*dz2)*pzm
              coeffs(2,0,ix-1,iz) = 0.
            endif
            if (ix < nxlocal) then
              rhs(ix+1,iz) = rhs(ix+1,iz)
     &           + 0.5*(rr(ix)+0.5*dx)/(rr(ix)+dx)*chi(1,0,ix,iz)/(dx*dz2)*pzm
              coeffs(0,0,ix+1,iz) = 0.
            endif
          else
            if (ix > 0) then
              rhs(ix-1,iz) = rhs(ix-1,iz) - 0.5*chi(1,0,ix,iz)/(dx*dz2)*pzm
              coeffs(2,0,ix-1,iz) = 0.
            endif
            if (ix < nxlocal) then
              rhs(ix+1,iz) = rhs(ix+1,iz) + 0.5*chi(1,0,ix,iz)/(dx*dz2)*pzm
              coeffs(0,0,ix+1,iz) = 0.
            endif
          endif
        elseif (-1. < dels(4,ic) .and. dels(4,ic) <= 0.) then
          coeffs(1,0,ix,iz) = 0.
          if (ix > 0) coeffs(2,0,ix-1,iz) = 0.
          if (ix < nxlocal) coeffs(0,0,ix+1,iz) = 0.
        endif

        if (0. < dels(5,ic) .and. dels(5,ic) < +1.) then
          pzp = voltfac*volt(5,ic)
          rhs(ix,iz) = rhs(ix,iz) + coeffs(1,2,ix,iz)*pzp
          coeffs(1,2,ix,iz) = 0.
          if (lrz) then
            if (ix > 1 .or. (ix == 1 .and. xmminlocal > 0.)) then
              rhs(ix-1,iz) = rhs(ix-1,iz)
     &           + 0.5*(rr(ix)-0.5*dx)/(rr(ix)-dx)*chi(1,0,ix,iz)/(dx*dz2)*pzp
              coeffs(2,2,ix-1,iz) = 0.
            endif
            if (ix < nxlocal) then
              rhs(ix+1,iz) = rhs(ix+1,iz)
     &           - 0.5*(rr(ix)+0.5*dx)/(rr(ix)+dx)*chi(1,0,ix,iz)/(dx*dz2)*pzp
              coeffs(0,2,ix+1,iz) = 0.
            endif
          else
            if (ix > 0) then
              rhs(ix-1,iz) = rhs(ix-1,iz) + 0.5*chi(1,0,ix,iz)/(dx*dz2)*pzp
              coeffs(2,2,ix-1,iz) = 0.
            endif
            if (ix < nxlocal) then
              rhs(ix+1,iz) = rhs(ix+1,iz) - 0.5*chi(1,0,ix,iz)/(dx*dz2)*pzp
              coeffs(0,2,ix+1,iz) = 0.
            endif
          endif
        elseif (-1. < dels(5,ic) .and. dels(5,ic) <= 0.) then
          coeffs(1,2,ix,iz) = 0.
          if (ix > 0) coeffs(2,2,ix-1,iz) = 0.
          if (ix < nxlocal) coeffs(0,2,ix+1,iz) = 0.
        endif

      enddo
!$OMP END DO

      end subroutine getcoefficientssubgrides2d
      end module coefficientsmodulees2d
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine mgsolveimplicites2d(iwhich,nx,nz,nxlocal,nzlocal,dx,dz,
     &                               phi,rho,ns,qomdt,chi0,
     &                               bounds,xmminlocal,zmminlocal,zgrid,
     &                               mgparam,mgiters,mgmaxiters,
     &                               mgmaxlevels,mgerror,mgtol,mgverbose,
     &                               downpasses,uppasses,
     &                               lcndbndy,laddconductor,icndbndy,
     &                               gridmode,conductors,lrz,fsdecomp)
      use Subtimersfrz
      use ConductorTypemodule
      use Constant
      use coefficientsmodulees2d
      use Decompositionmodule
      integer(ISZ):: iwhich
      integer(ISZ):: nx,nz,nxlocal,nzlocal,ns
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1)
      real(kind=8):: rho(0:nxlocal,0:nzlocal)
      real(kind=8):: chi0(0:nxlocal,0:nzlocal,0:ns-1)
      real(kind=8):: qomdt(0:ns-1)
      real(kind=8):: dx,dz
      integer(ISZ):: bounds(0:5)
      real(kind=8):: xmminlocal,zmminlocal,zgrid
      real(kind=8):: mgparam
      integer(ISZ):: mgiters,mgmaxiters,mgmaxlevels
      real(kind=8):: mgerror,mgtol
      integer(ISZ):: mgverbose,downpasses,uppasses
      logical(ISZ):: lcndbndy,laddconductor
      integer(ISZ):: icndbndy,gridmode
      type(ConductorType):: conductors
      logical(ISZ):: lrz
      type(Decomposition):: fsdecomp

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

      integer(ISZ):: i,k,ix,iz,delx,delz
      real(kind=8),pointer:: bongrid(:,:,:)
      real(kind=8),pointer:: chi0temp(:,:,:)
      real(kind=8),allocatable:: phisave(:,:)
      integer(ISZ):: localbounds(0:5)
      character(72):: errline
      integer(ISZ):: allocerror
      real(kind=8):: substarttime,wtime
      if (lfrztimesubs) substarttime = wtime()

c     --- If doing initialization only, then exit.
      if (iwhich == 1) return

      delx = 1
      delz = 1
      if (nx > nxlocal) delx = 3
      if (nz > nzlocal) delz = 3

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

c     --- Fill the B field array
c     --- This will be passed to vcycle and coarsened appropriately.
      bongrid => getbfieldsongrides2d(nxlocal,nzlocal,delx,delz,dx,dz,
     &                                xmminlocal,zmminlocal,zgrid,localbounds)
c     --- Parallel B.C.s are set below.

c     --- Create copy of rho with extra guard cells, and fill it in.
      allocate(chi0temp(-delx:nxlocal+delx,-delz:nzlocal+delz,0:ns-1))
      chi0temp = 0.
      chi0temp(0:nxlocal,0:nzlocal,:) = chi0
      call applyboundaryconditions3d(nxlocal,0,nzlocal,delx,0,delz,chi0temp,ns,
     &                               localbounds,.false.,.true.)
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
      call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,
     &                               localbounds,.true.,.false.)

#ifdef MPIPARALLEL
      call applyparallelboundaryconditions3d(chi0temp,nxlocal,0,nzlocal,
     &                                       delx,0,delz,ns,3,localbounds,fsdecomp)
      call applyparallelboundaryconditions3d(bongrid,nxlocal,0,nzlocal,
     &                                       delx,0,delz,3,3,localbounds,fsdecomp)
#endif

      allocate(phisave(-1:nxlocal+1,-1:nzlocal+1),stat=allocerror)
      if (allocerror /= 0) then
        print*,"mgsolveimplicites2d: allocation error ",allocerror,
     &         ": could not allocate phisave to shape ",nxlocal,nzlocal
        call kaboom("mgsolveimplicites2d: allocation error")
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
        call vcyclees2d(0,1.,nx,nz,nxlocal,nzlocal,delx,delz,dx,dz,xmminlocal,
     &                  phi,rho,
     &                  ns,qomdt,chi0temp,bongrid,bounds,mgparam,mgmaxlevels,
     &                  downpasses,uppasses,lcndbndy,icndbndy,conductors,lrz,
     &                  fsdecomp)

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
c       print*,mgiters,mgerror

c       --- This line below seems to create a large temporary which can
c       --- cause problems when memory is close to full. So it was replaced
c       --- with the explicit loop above.
c       mgerror = maxval(abs(phisave - phi))

      enddo

c     --- Set boundary conditions. This is only really needed for the
c     --- Dirichlet boundaries, but this is convenient to call.
      call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,
     &                               localbounds,.true.,.false.)

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
      RECURSIVE subroutine vcyclees2d(mglevel,mgscale,nx,nz,nxlocal,nzlocal,
     &                                delx,delz,dx,dz,xmminlocal,
     &                                phi,rho,ns,qomdt,chi0,bongrid,
     &                                globalbounds,mgparam,
     &                                mgmaxlevels,downpasses,uppasses,
     &                                lcndbndy,icndbndy,conductors,lrz,
     &                                fsdecomp)
      use Constant,Only:eps0
      use ConductorTypemodule
      use Multigrid3d_diagnostic
      use coefficientsmodulees2d
      use formggetarraysuminterface
      use Decompositionmodule
c     use ImplicitMG2D
      integer(ISZ):: mglevel
      real(kind=8):: mgscale
      integer(ISZ):: nx,nz,nxlocal,nzlocal,delx,delz,ns
      real(kind=8):: dx,dz,xmminlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: chi0(-delx:nxlocal+delx,-delz:nzlocal+delz,0:ns-1)
      real(kind=8):: bongrid(-delx:nxlocal+delx,-delz:nzlocal+delz,0:2)
      real(kind=8):: qomdt(0:ns-1)
      integer(ISZ):: globalbounds(0:5)
      real(kind=8):: mgparam
      integer(ISZ):: mgmaxlevels,downpasses,uppasses
      type(ConductorType):: conductors
      logical(ISZ):: lcndbndy,lrz
      integer(ISZ):: icndbndy
      type(Decomposition):: fsdecomp

c Routine that does the v-cycle for multigrid. Note that it is recursive.

      real(kind=8):: mingridsize
      real(kind=8),pointer:: coeffs(:,:,:,:),rhs(:,:),resdelfac(:,:)
      real(kind=8),allocatable:: phicoarse(:,:),rhocoarse(:,:)
      real(kind=8),allocatable:: chi0coarse(:,:,:),bongridcoarse(:,:,:)
      real(kind=8),allocatable:: res(:,:)
      integer(ISZ):: i,ic,iszone=1
      integer(ISZ):: nxcoarse,nycoarse,nzcoarse,nylocal
      integer(ISZ):: nxlocalcoarse,nzlocalcoarse
      real(kind=8):: dxcoarse,dycoarse,dzcoarse
      real(kind=8):: dxcoarsesqi,dzcoarsesqi
      real(kind=8):: xmminlocalcoarse
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
      logical(ISZ):: lgoserial,lpe0
      integer(ISZ):: allocerror
      real(kind=8):: sss(2)

      localbounds = globalbounds

#ifdef MPIPARALLEL
      ixproc = fsdecomp%ixproc
      izproc = fsdecomp%izproc
      if (fsdecomp%ix(ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iz(izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(izproc)+nzlocal < nz) localbounds(5) = -1
#endif

c     --- Get the finite difference coefficients for the Poisson operator at
c     --- the current level of refinement.
      allocate(coeffs(0:2,0:2,0:nxlocal,0:nzlocal))
      allocate(rhs(0:nxlocal,0:nzlocal))
      allocate(resdelfac(0:nxlocal,0:nzlocal))
      call getcoefficientses2d(nxlocal,nzlocal,delx,delz,dx,dz,xmminlocal,
     &                         ns,qomdt,rho,chi0,bongrid,
     &                         mglevel,lcndbndy,icndbndy,conductors,lrz,
     &                         localbounds,coeffs,rhs,fsdecomp,resdelfac)

      if (lprintmgarraysumdiagnostic) then
        lpe0=(fsdecomp%ixproc==0.and.fsdecomp%izproc==0)
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V1 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V1 rho",mglevel,sss/eps0
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rhs,fsdecomp)
        if (lpe0) print*,"V1 rhs",mglevel,sss/eps0
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(0,0,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs00",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(1,0,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs10",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(2,0,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs20",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(0,1,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs01",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(1,1,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs11",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(2,1,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs21",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(0,2,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs02",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(1,2,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs12",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,coeffs(2,2,:,:),fsdecomp)
        if (lpe0) print*,"V1 coeffs22",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,delx,0,delz,chi0,fsdecomp)
        if (lpe0) print*,"V1 chi0",mglevel,sss/eps0
        sss = mggetarraysum(nxlocal,0,nzlocal,delx,0,delz,bongrid(:,:,0),fsdecomp)
        if (lpe0) print*,"V1 bongrid0",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,delx,0,delz,bongrid(:,:,1),fsdecomp)
        if (lpe0) print*,"V1 bongrid1",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,delx,0,delz,bongrid(:,:,2),fsdecomp)
        if (lpe0) print*,"V1 bongrid2",mglevel,sss
      endif

c     --- Do initial SOR passes.
      do i=1,downpasses
        call relaximplicites2d(mglevel,nx,nz,nxlocal,nzlocal,phi,rho,rhs,coeffs,
     &                         localbounds,mgparam,conductors,fsdecomp)
      enddo

      if (lprintmgarraysumdiagnostic) then
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V2 phi",mglevel,sss
        sss = mggetarraysum(nxlocal,0,nzlocal,0,0,0,rho,fsdecomp)
        if (lpe0) print*,"V2 rho",mglevel,sss/eps0
      endif

c     --- Check if this is the finest level. If so, then don't do any further
c     --- coarsening. This is the same check that is done in getmglevels.
      if (nx >= 4 .and. nz >= 4 .and. mglevel < mgmaxlevels) then

        allocate(res(-delx:nxlocal+delx,-delz:nzlocal+delz),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate res to shape ",nxlocal,nzlocal
          call kaboom("vcyclees2d: allocation error")
          return
        endif

c       --- Get the residual on the current grid.
        call residualimplicites2d(nxlocal,nzlocal,delx,delz,
     &                            phi,rho,coeffs,rhs,resdelfac,res,
     &                            mglevel,localbounds,mgparam,conductors)
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
        mgscalecoarse = mgscale*dxcoarse*dzcoarse/(dx*dz)

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
        lgoserial = (maxval(lxoffsetall) > 2*nxcoarse .or.
     &               maxval(rxoffsetall) > 2*nxcoarse .or.
     &               maxval(lzoffsetall) > 2*nzcoarse .or.
     &               maxval(rzoffsetall) > 2*nzcoarse)
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
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate phicoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(rhocoarse(0:nxlocalcoarse,0:nzlocalcoarse),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate rhocoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(chi0coarse(-delx:nxlocalcoarse+delx,
     &                      -delz:nzlocalcoarse+delz,0:ns-1),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate chi0coarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse
          call kaboom("vcyclees2d: allocation error")
          return
        endif
        allocate(bongridcoarse(-delx:nxlocalcoarse+delx,
     &                         -delz:nzlocalcoarse+delz,0:2),stat=allocerror)
        if (allocerror /= 0) then
          print*,"vcyclees2d: allocation error ",allocerror,
     &           ": could not allocate bongridcoarse to shape ",
     &           nxlocalcoarse,nzlocalcoarse,3
          call kaboom("vcyclees2d: allocation error")
          return
        endif

        rhocoarse = 0.
        phicoarse = 0.
        chi0coarse = 0.
        bongridcoarse = 0.

c       --- Restriction
        call restrict2d(nx,nz,nxlocal,nzlocal,res,delx,delz,
     &                  nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                  rhocoarse,0,0,
     &                  localbounds,localboundsc,lxoffset,lzoffset)

        do ic=0,ns-1
          call restrict2d(nx,nz,nxlocal,nzlocal,chi0(:,:,ic),delx,delz,
     &                    nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                    chi0coarse(:,:,ic),delx,delz,
     &                    localbounds,localboundsc,lxoffset,lzoffset)
        enddo
        do ic=0,2
          call restrict2d(nx,nz,nxlocal,nzlocal,bongrid(:,:,ic),delx,delz,
     &                    nxcoarse,nzcoarse,nxlocalcoarse,nzlocalcoarse,
     &                    bongridcoarse(:,:,ic),delx,delz,
     &                    localbounds,localboundsc,lxoffset,lzoffset)
        enddo
        call applyboundaryconditions3d(nxlocalcoarse,0,nzlocalcoarse,
     &                                 delx,0,delz,chi0coarse,ns,
     &                                 localboundsc,.false.,.true.)
        call applyboundaryconditions3d(nxlocalcoarse,0,nzlocalcoarse,
     &                                 delx,0,delz,bongridcoarse,3,
     &                                 localboundsc,.false.,.true.)
#ifdef MPIPARALLEL
        call applyparallelboundaryconditions3d(chi0coarse,
     &               nxlocalcoarse,0,nzlocalcoarse,
     &               delx,0,delz,ns,3,localboundsc,coarsedecomp)
        call applyparallelboundaryconditions3d(bongridcoarse,
     &               nxlocalcoarse,0,nzlocalcoarse,
     &               delx,0,delz,3,3,localboundsc,coarsedecomp)
#endif
        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,0,0,0,rhocoarse,
     &                        coarsedecomp)
          if (lpe0) print*,"V3 rhocoarse",mglevel,sss/eps0
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,delx,0,delz,
     &                        chi0coarse,coarsedecomp)
          if (lpe0) print*,"V3 chi0coarse",mglevel,sss/eps0
        endif

c       --- Continue at the next coarsest level.
        call vcyclees2d(mglevel+iszone,mgscalecoarse,nxcoarse,nzcoarse,
     &                  nxlocalcoarse,nzlocalcoarse,delx,delz,
     &                  dxcoarse,dzcoarse,xmminlocalcoarse,
     &                  phicoarse,rhocoarse,ns,qomdt,chi0coarse,bongridcoarse,
     &                  globalbounds,mgparam,
     &                  mgmaxlevels,downpasses,uppasses,
     &                  lcndbndy,icndbndy,conductors,lrz,coarsedecomp)

        if (lprintmgarraysumdiagnostic) then
          sss = mggetarraysum(nxlocalcoarse,0,nzlocalcoarse,1,0,1,phicoarse,
     &                        coarsedecomp)
          if (lpe0) print*,"V4 phicoarse",mglevel,sss
        endif

c       --- Add in resulting error.
        call expand3d(nx,0,nz,nxlocal,0,nzlocal,1,0,1,phi,
     &                nxcoarse,0,nzcoarse,nxlocalcoarse,0,nzlocalcoarse,
     &                phicoarse,localbounds,lxoffset,0,lzoffset)
        call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,
     &                                 localbounds,.false.,.false.)
#ifdef MPIPARALLEL
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,-1,-1,localbounds,fsdecomp)
#endif
        call averageperiodicphi3d(nx,0,nz,nxlocal,0,nzlocal,phi,1,0,1,
     &                            globalbounds)
        
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
        call mgexchange_phi_periodic(1,nxlocal,0,nzlocal,phi,
     &                               1,0,1,0,0,localbounds,fsdecomp)
        call mgexchange_phi(1,nxlocal,0,nzlocal,phi,1,0,1,-1,-1,fsdecomp)
#endif
        sss = mggetarraysum(nxlocal,0,nzlocal,1,0,1,phi,fsdecomp)
        if (lpe0) print*,"V5 phi",mglevel,sss
      endif

c     --- Do final SOR passes.
      do i=1,uppasses
        call relaximplicites2d(mglevel,nx,nz,nxlocal,nzlocal,phi,rho,rhs,coeffs,
     &                         localbounds,mgparam,conductors,fsdecomp)
      enddo

      deallocate(coeffs)
      deallocate(rhs)
      deallocate(resdelfac)

      return
      end
c=============================================================================
      subroutine relaximplicites2d(mglevel,nx,nz,nxlocal,nzlocal,
     &                             phi,rho,rhs,coeffs,localbounds,mgparam,
     &                             conductors,fsdecomp)
      use Constant
      use ConductorTypemodule
      use Decompositionmodule
      integer(ISZ):: mglevel,nx,nz,nxlocal,nzlocal
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:2,0:nxlocal,0:nzlocal)
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: mgparam,phiav
      type(ConductorType):: conductors
      type(Decomposition):: fsdecomp

c This routine does one pass of point SOR with even-odd (red-black)
c ordering.
c phisave is needed because when there are B fields, the even and odd sweeps
c are no longer independent since the diagonal coefficients become non-zero.
c The simplest solution is to save the values of phi before a sweep and
c read those saved values for the update of phi.

      real(kind=8),allocatable:: phisave(:,:)
      integer(ISZ):: parity,s_parity,e_parity
      integer(ISZ):: ixmin,ixmax,izmin,izmax,ix,iz,ix1

      allocate(phisave(-1:nxlocal+1,-1:nzlocal+1))

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
     &                  nxlocal,0,nzlocal,1,0,1,phi,mglevel,1,.false.)

c       --- set phi in the z guard planes
        call applyboundaryconditions3d(nxlocal,0,nzlocal,1,0,1,phi,1,
     &                                 localbounds,.false.,.false.)

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

      call averageperiodicphi3d(nx,0,nz,nxlocal,0,nzlocal,phi,1,0,1,localbounds)

      deallocate(phisave)
      return
      end
c=============================================================================
      subroutine residualimplicites2d(nxlocal,nzlocal,delx,delz,
     &                                phi,rho,coeffs,rhs,resdelfac,res,
     &                                mglevel,localbounds,mgparam,conductors)
      use Constant
      use ConductorTypemodule
      integer(ISZ):: nxlocal,nzlocal,delx,delz
      real(kind=8):: phi(-1:nxlocal+1,-1:nzlocal+1),rho(0:nxlocal,0:nzlocal)
      real(kind=8):: coeffs(0:2,0:2,0:nxlocal,0:nzlocal)
      real(kind=8):: rhs(0:nxlocal,0:nzlocal)
      real(kind=8):: resdelfac(0:nxlocal,0:nzlocal)
      real(kind=8):: res(-delx:nxlocal+delx,-delz:nzlocal+delz)
      integer(ISZ):: mglevel,localbounds(0:5)
      real(kind=8):: mgparam
      type(ConductorType):: conductors

c Calculate the residual on the grid. Residual = r.h.s. - l.h.s.

      integer(ISZ):: ix,iz
      integer(ISZ):: ixmin,ixmax,izmin,izmax

c     --- Set the loop limits, including edges when appropriate.
      ixmin = 0
      ixmax = nxlocal
      izmin = 0
      izmax = nzlocal
      if (localbounds(0) == 0) ixmin = 1
      if (localbounds(1) == 0) ixmax = nxlocal - 1
      if (localbounds(4) == 0) izmin = 1
      if (localbounds(5) == 0) izmax = nzlocal - 1

      res = 0.
c     --- Calculate the residual.
!$OMP DO
      do iz=izmin,izmax

        do ix=ixmin,ixmax

          res(ix,iz) = (rho(ix,iz) + eps0*(rhs(ix,iz) +
     &                 coeffs(0,0,ix,iz)*phi(ix-1,iz-1) +
     &                 coeffs(1,0,ix,iz)*phi(ix  ,iz-1) +
     &                 coeffs(2,0,ix,iz)*phi(ix+1,iz-1) +
     &                 coeffs(0,1,ix,iz)*phi(ix-1,iz  ) +
     &                 coeffs(1,1,ix,iz)*phi(ix  ,iz  ) +
     &                 coeffs(2,1,ix,iz)*phi(ix+1,iz  ) +
     &                 coeffs(0,2,ix,iz)*phi(ix-1,iz+1) +
     &                 coeffs(1,2,ix,iz)*phi(ix  ,iz+1) +
     &                 coeffs(2,2,ix,iz)*phi(ix+1,iz+1)))
     &                 *resdelfac(ix,iz)

        enddo
      enddo
!$OMP END DO

      call cond_potmgres(conductors%interior,
     &                   nxlocal,0,nzlocal,delx,0,delz,res,mglevel,1,.false.)
      call applyboundaryconditions3d(nxlocal,0,nzlocal,delx,0,delz,res,1,
     &                               localbounds,.false.,.false.)

      return
      end
c=============================================================================
c=============================================================================
