#include "top.h"
c=============================================================================
c@(#) File CHO.F, version $Revision: 1.2 $, $Date: 2002/08/27 23:11:36 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the interface to the Chombo/AMR fieldsolver.
c  David P. Grote, LLNL, (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine cho_solve3d(iwhich,nx,ny,nz,nzfull,dx,dy,dz,
     &                       l2symtry,l4symtry,xmmin,ymmin,zmmin,zmminglobal)
      use CHOHandle
      use CHOInput
      use Conductor3d
      use Parallel
      integer(ISZ):: iwhich
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      logical(ISZ):: l2symtry,l4symtry
      real(kind=8):: xmmin,ymmin,zmmin,zmminglobal

c Interface to the Chombo/AMR fieldsolver.

      real(kind=8):: x0(3)
      integer(ISZ):: domain(3,2)
      integer(ISZ):: subdomains(3,2,0:maxslaves-1)
      integer(ISZ):: refratios(cho_nlevels-1)

#ifdef CHOMBO

      if (iwhich == 0 .or. iwhich == 1) then

        call chp_create(cho_handle,cho_status)
        if (cho_status < 0) then
          print*,"Fatal error from chp_create"
          stop
        endif

        x0(1) = xmmin
        x0(2) = ymmin
        x0(3) = zmminglobal
        domain(:,1) = 0
        domain(1,2) = nx
        domain(2,2) = ny
        domain(3,2) = nzfull
        subdomains(1,1,:) = 0
        subdomains(1,2,:) = nx
        subdomains(2,1,:) = 0
        subdomains(2,2,:) = ny
        subdomains(3,1,:) = izpslave
        subdomains(3,2,:) = izpslave + nzpslave
        refratios = cho_refratio
        call chp_setgridparams(cho_handle,x0,dx,domain,subdomains,
     &                         cho_nlevels,refratios,cho_maxparticlespercell,
     &                         cho_status)
        if (cho_status < 0) then
          print*,"Fatal error from chombo"
          stop
        endif

        if (l2symtry) then
          cho_bcflags(2,1) = 1
        elseif (l4symtry) then
          cho_bcflags(1,1) = 1
          cho_bcflags(2,1) = 1
        endif
        call chp_solverparams(cho_handle,cho_tol,cho_bcflags,cho_status)
        if (cho_status < 0) then
          print*,"Fatal error from chp_solverparams"
          stop
        endif

      endif

      if (iwhich == 0 .or. iwhich == -1) then
        call chp_solve(cho_handle)
        if (cho_status < 0) then
          print*,"Fatal error from chp_solve"
          stop
        endif
      endif

#endif

      return
      end
c=============================================================================
      subroutine cho_setrho3d(np,xp,yp,zp,zgrid,sq,sw,js,ip)
      use CHOHandle
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: zgrid,sq,sw
      integer(ISZ):: js,ip

c Interface to the chp_putparticles routine.

#ifdef CHOMBO

c     --- Update bunchid. If this is the first batch of particles, then
c     --- set bunchid = 1, otherwise increment it. If js == -1, then
c     --- set bunchid = -1.
      if (js == 1 .and. ip == 0) then
        cho_bunchid = 1
      else if (js == -1) then
        cho_bunchid = -1
      else
        cho_bunchid = cho_bunchid + 1
      endif

c     --- This may not be the best way of shifting particles to the
c     --- grid frame. Though this does avoid a large temporary
c     --- array, there maybe a small problem with roundoff.
      zp = zp - zgrid

      call chp_putparticles(cho_handle,cho_bunchid,sq*sw,np,1,xp,yp,zp,
     &                      cho_status)
      if (cho_status < 0) then
        print*,"Fatal error from chp_putparticles"
        stop
      endif

c     --- Reset the particles z position.
      zp = zp + zgrid

#endif

      return
      end
c=============================================================================
      subroutine cho_gete3d(np,xp,yp,zp,zgrid,ex,ey,ez,js,ip)
      use CHOHandle
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      real(kind=8):: zgrid
      integer(ISZ):: js,ip

c Interface to the chp_getefield routine.

#ifdef CHOMBO

c     --- Update bunchid. If this is the first batch of particles, then
c     --- set bunchid = 1, otherwise increment it. If js == -1, then
c     --- set bunchid = -1.
      if (js == 1 .and. ip == 0) then
        cho_bunchid = 1
      else if (js == -1) then
        cho_bunchid = -1
      else
        cho_bunchid = cho_bunchid + 1
      endif

c     --- This may not be the best way of shifting particles to the
c     --- grid frame. Though this does avoid a large temporary
c     --- array, there maybe a small problem with roundoff.
      zp = zp - zgrid

      call chp_getefield(cho_handle,cho_bunchid,np,1,xp,yp,zp,ex,ey,ez,
     &                   cho_status)
      if (cho_status < 0) then
        print*,"Fatal error from chp_getefield"
        stop
      endif

c     --- Reset the particles z position.
      zp = zp + zgrid

#endif

      return
      end
c=============================================================================
      subroutine cho_getphi3d(np,xp,yp,zp,zgrid,p,js,ip)
      use CHOHandle
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: p(np)
      real(kind=8):: zgrid
      integer(ISZ):: js,ip

c Interface to the chp_getphi routine.

#ifdef CHOMBO

c     --- Update bunchid. If this is the first batch of particles, then
c     --- set bunchid = 1, otherwise increment it. If js == -1, then
c     --- set bunchid = -1.
      if (js == 1 .and. ip == 0) then
        cho_bunchid = 1
      else if (js == -1) then
        cho_bunchid = -1
      else
        cho_bunchid = cho_bunchid + 1
      endif

c     --- This may not be the best way of shifting particles to the
c     --- grid frame. Though this does avoid a large temporary
c     --- array, there maybe a small problem with roundoff.
      zp = zp - zgrid

      call chp_getphi(cho_handle,cho_bunchid,np,1,xp,yp,zp,p,cho_status)
      if (cho_status < 0) then
        print*,"Fatal error from chp_getphi"
        stop
      endif

c     --- Reset the particles z position.
      zp = zp + zgrid

#endif

      return
      end
c=============================================================================
      subroutine cho_getrho3d(np,xp,yp,zp,zgrid,r,js,ip)
      use CHOHandle
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: r(np)
      real(kind=8):: zgrid
      integer(ISZ):: js,ip

c Interface to the chp_getrho routine.

#ifdef CHOMBO

c     --- Update bunchid. If this is the first batch of particles, then
c     --- set bunchid = 1, otherwise increment it. If js == -1, then
c     --- set bunchid = -1.
      if (js == 1 .and. ip == 0) then
        cho_bunchid = 1
      else if (js == -1) then
        cho_bunchid = -1
      else
        cho_bunchid = cho_bunchid + 1
      endif

c     --- This may not be the best way of shifting particles to the
c     --- grid frame. Though this does avoid a large temporary
c     --- array, there maybe a small problem with roundoff.
      zp = zp - zgrid

      call chp_getrho(cho_handle,cho_bunchid,np,1,xp,yp,zp,r,cho_status)
      if (cho_status < 0) then
        print*,"Fatal error from chp_getrho"
        stop
      endif

c     --- Reset the particles z position.
      zp = zp + zgrid

#endif

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine ReachableNodes(dx, mask, xlo, ylo, zlo, xhi, yhi, zhi, ncomp)
      use PSOR3d
      use Conductor3d
      integer(ISZ):: xlo, ylo, zlo, xhi, yhi, zhi, ncomp
      integer(ISZ):: mask(xlo:xhi, ylo:yhi, zlo:zhi, 0:6)
      real(kind=8):: dx

c This subroutine sets array mask:
c
c mask(i, j, k, 0) = 1 if (i, j, k) is in the domain (uncovered); else 0
c               1    1 if (i+1, j  , k  ) is reachable from (i, j, k); else 0
c               2    1 if (i-1, j  , k  ) is reachable from (i, j, k); else 0
c               3    1 if (i  , j+1, k  ) is reachable from (i, j, k); else 0
c               4    1 if (i  , j-1, k  ) is reachable from (i, j, k); else 0
c               5    1 if (i  , j  , k+1) is reachable from (i, j, k); else 0
c               6    1 if (i  , j  , k-1) is reachable from (i, j, k); else 0

      integer(ISZ):: level
      integer(ISZ):: ic,ix,iy,iz,il

c     --- Get level of coarseness
      level = nint(dx/dxfine)
      il = int(log(1.*level)/log(2.) + 0.5)

c     --- By default, all points are within the domain and are reachable.
      mask = 1

c     --- Check for point in conductors.
      do ic=icstart(il),icstart(il+1)-1
        ix = ixcond(ic)
        iy = iycond(ic)
        iz = izcond(ic)
c       --- Then check if the point is within the requested region.
c       --- If so, set mask to 0.
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) mask(ix,iy,iz,0:6) = 0
c       --- Set mask to 0 for neighboring points
        if (xlo <= ix-1 .and. ix-1 <= xhi .and.
     &      ylo <= iy   .and. iy   <= yhi .and.
     &      zlo <= iz   .and. iz   <= zhi)  mask(ix-1,iy  ,iz  ,1) = 0
        if (xlo <= ix+1 .and. ix+1 <= xhi .and.
     &      ylo <= iy   .and. iy   <= yhi .and.
     &      zlo <= iz   .and. iz   <= zhi)  mask(ix+1,iy  ,iz  ,2) = 0
        if (xlo <= ix   .and. ix   <= xhi .and.
     &      ylo <= iy-1 .and. iy-1 <= yhi .and.
     &      zlo <= iz   .and. iz   <= zhi)  mask(ix  ,iy-1,iz  ,3) = 0
        if (xlo <= ix   .and. ix   <= xhi .and.
     &      ylo <= iy+1 .and. iy+1 <= yhi .and.
     &      zlo <= iz   .and. iz   <= zhi)  mask(ix  ,iy+1,iz  ,4) = 0
        if (xlo <= ix   .and. ix   <= xhi .and.
     &      ylo <= iy   .and. iy   <= yhi .and.
     &      zlo <= iz-1 .and. iz-1 <= zhi)  mask(ix  ,iy  ,iz-1,5) = 0
        if (xlo <= ix   .and. ix   <= xhi .and.
     &      ylo <= iy   .and. iy   <= yhi .and.
     &      zlo <= iz+1 .and. iz+1 <= zhi)  mask(ix  ,iy  ,iz+1,6) = 0
      enddo

c     --- Check points near boundary.
      call ReachableNodeswork(mask,xlo,ylo,zlo,xhi,yhi,zhi,
     &                        necndbdy,
     &                        ecdelmx,ecdelmy,ecdelmz,ecdelpx,ecdelpy,ecdelpz,
     &                        iecndx,iecndy,iecndz,iecndlevel,
     &                        level,il,ecstart)
      call ReachableNodeswork(mask,xlo,ylo,zlo,xhi,yhi,zhi,
     &                        nocndbdy,
     &                        ocdelmx,ocdelmy,ocdelmz,ocdelpx,ocdelpy,ocdelpz,
     &                        iocndx,iocndy,iocndz,iocndlevel,
     &                        level,il,ocstart)

      return
      end
c=============================================================================
      subroutine ReachableNodeswork(mask,xlo,ylo,zlo,xhi,yhi,zhi,
     &                              ncndbdy,
     &                              cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,
     &                              cdelpz,icndx,icndy,icndz,icndlevel,
     &                              level,il,cstart)
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi
      integer(ISZ):: mask(xlo:xhi, ylo:yhi, zlo:zhi, 0:6)
      integer(ISZ):: ncndbdy,level,il,cstart(0:100)
      real(kind=8):: cdelmx(ncndbdy), cdelmy(ncndbdy), cdelmz(ncndbdy)
      real(kind=8):: cdelpx(ncndbdy), cdelpy(ncndbdy), cdelpz(ncndbdy)
      integer(ISZ):: icndx(ncndbdy), icndy(ncndbdy), icndz(ncndbdy)
      integer(ISZ):: icndlevel(ncndbdy)

c Does the work of determining whether neighboring points are "reachable"
c from the central point. The even or odd conductor points can be passed
c into this routine.

      integer(ISZ):: ic,ix,iy,iz

c     --- loop over points near surface of conductors
      do ic = cstart(il),cstart(il+1)-1
 
        ix = icndx(ic)
        iy = icndy(ic)
        iz = icndz(ic)
 
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then

c         --- the point higher in x is inside the conductor
          if (cdelpx(ic) < 1.) mask(ix,iy,iz,1) = 0

c         --- the point lower in x is inside the conductor
          if (cdelmx(ic) < 1.) mask(ix,iy,iz,2) = 0

c         --- the point higher in y is inside the conductor
          if (cdelpy(ic) < 1.) mask(ix,iy,iz,3) = 0

c         --- the point lower in y is inside the conductor
          if (cdelmy(ic) < 1.) mask(ix,iy,iz,4) = 0

c         --- the point lower in z is inside the conductor
          if (cdelpz(ic) < 1.) mask(ix,iy,iz,5) = 0

c         --- the point lower in z is inside the conductor
          if (cdelmz(ic) < 1.) mask(ix,iy,iz,6) = 0

        endif
      enddo

      return
      end
c=============================================================================
c=============================================================================
      subroutine CoveredNodes(dx, mask, xlo, ylo, zlo, xhi, yhi, zhi)
      use PSOR3d
      use Conductor3d
      integer(ISZ):: xlo, ylo, zlo, xhi, yhi, zhi
      integer(ISZ):: mask(xlo:xhi, ylo:yhi, zlo:zhi)
      real(kind=8):: dx
 
c Returns the mask, which is set to 1 for points which are inside a
c conductor. The input argument dx is used to determine the refinement level.
 
      integer(ISZ):: level,ic,il
 
c     --- Get the level of coarseness above the finest level.
      level = nint(dx/dxfine)
      il = int(log(1.*level)/log(2.) + 0.5)

c     --- Zero out the mask
      mask = 0

      do ic=icstart(il),icstart(il+1)-1
c       --- Then check if the point is within the requested region.
        if (xlo <= ixcond(ic) .and. ixcond(ic) <= xhi .and.
     &      ylo <= iycond(ic) .and. iycond(ic) <= yhi .and.
     &      zlo <= izcond(ic) .and. izcond(ic) <= zhi) then
c         --- If so, set mask to 1.
          mask(ixcond(ic),iycond(ic),izcond(ic)) = 1
        endif
      enddo
 
      return
      end
c=============================================================================
      subroutine NodalCoefficients(dx,coeffs,xlo,ylo,zlo,xhi,yhi,zhi,ncomp)
      use PSOR3d
      use Conductor3d
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi,ncomp
      real(kind=8):: dx,coeffs(xlo:xhi,ylo:yhi,zlo:zhi,0:7)

c Fills the array coeffs with the coeeficients of the 7-point finite
c differenced Poisson's equation, taking into account the changed form of the
c equation near a boundary using the subgrid-scale technique.

      integer(ISZ):: ic,ix,iy,iz,il
      integer(ISZ):: level
      real(kind=8):: rdx2,coeff0

c     --- First put in the default values
      rdx2 = 1./dx**2
      coeff0 = -2.*(rdx2 + rdx2 + rdx2)
      coeffs(:,:,:,0) = coeff0
      coeffs(:,:,:,1:6) = rdx2
      coeffs(:,:,:,7) = 0.
      
c     --- Get level of coarseness
      level = nint(dx/dxfine)
      il = int(log(1.*level)/log(2.) + 0.5)

c     --- Check for point in conductors.
      do ic=icstart(il),icstart(il+1)-1
        ix = ixcond(ic)
        iy = iycond(ic)
        iz = izcond(ic)
c       --- Then check if the point is within the requested region.
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then
          coeffs(ix,iy,iz,0) = 1.
          coeffs(ix,iy,iz,1:6) = 0.
          coeffs(ix,iy,iz,7) = condvolt(ic)
        endif
      enddo

c     --- Get coefficients for even points.
      call NodalCoefficientswork(coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,
     &                        necndbdy,
     &                        ecdelmx,ecdelmy,ecdelmz,ecdelpx,ecdelpy,ecdelpz,
     &                        iecndx,iecndy,iecndz,iecndlevel,
     &                        ecvoltmx,ecvoltpx,ecvoltmy,
     &                        ecvoltpy,ecvoltmz,ecvoltpz,
     &                        level,il,ecstart)

c     --- Get coefficients for odd points.
      call NodalCoefficientswork(coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,
     &                        nocndbdy,
     &                        ocdelmx,ocdelmy,ocdelmz,ocdelpx,ocdelpy,ocdelpz,
     &                        iocndx,iocndy,iocndz,iocndlevel,
     &                        ocvoltmx,ocvoltpx,ocvoltmy,
     &                        ocvoltpy,ocvoltmz,ocvoltpz,
     &                        level,il,ecstart)

      return
      end
c=============================================================================
      subroutine NodalCoefficientswork(coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,
     &                        ncndbdy,
     &                        cdelmx,cdelmy,cdelmz,cdelpx,cdelpy,cdelpz,
     &                        icndx,icndy,icndz,icndlevel,
     &                        cvoltmx,cvoltpx,cvoltmy,cvoltpy,cvoltmz,cvoltpz,
     &                        level,il,cstart)
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi
      real(kind=8):: coeffs(xlo:xhi,ylo:yhi,zlo:zhi,0:7)
      real(kind=8):: rdx2
      integer(ISZ):: ncndbdy,level,il,cstart(0:100)
      real(kind=8):: cdelmx(ncndbdy), cdelmy(ncndbdy), cdelmz(ncndbdy)
      real(kind=8):: cdelpx(ncndbdy), cdelpy(ncndbdy), cdelpz(ncndbdy)
      integer(ISZ):: icndx(ncndbdy), icndy(ncndbdy), icndz(ncndbdy)
      integer(ISZ):: icndlevel(ncndbdy)
      real(kind=8):: cvoltmx(ncndbdy),cvoltmy(ncndbdy),cvoltmz(ncndbdy)
      real(kind=8):: cvoltpx(ncndbdy),cvoltpy(ncndbdy),cvoltpz(ncndbdy)

c Does the work of setting the coefficients. The even or odd conductor
c points can be passed into this routine.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: coeff0(6),voltswitch

c     --- Voltages are only applied at level == 1, the finest level.
c     --- Otherwise, the voltages are set to zero, meaning that the
c     --- residual is zero
c     if (level == 1) then
        voltswitch = 1.
c     else
c       voltswitch = 0.
c     endif

c     --- loop over points near surface of conductors
      do ic = cstart(il),cstart(il+1)-1
 
        ix = icndx(ic)
        iy = icndy(ic)
        iz = icndz(ic)
 
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then
 
c         --- Contributions of terms with no conductor boundary nearby.
          coeff0 = -rdx2

c         --- the point higher in x is inside the conductor
          if (cdelpx(ic) < 1.) then
            coeff0(1) = -1./cdelpx(ic)*rdx2
            coeffs(ix,iy,iz,1) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - cvoltpx(ic)/cdelpx(ic)*rdx2*voltswitch
          endif
c         --- the point lower in x is inside the conductor
          if (cdelmx(ic) < 1.) then
            coeff0(2) = -1./cdelmx(ic)*rdx2
            coeffs(ix,iy,iz,2) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - cvoltmx(ic)/cdelmx(ic)*rdx2*voltswitch
          endif
c         --- the point higher in y is inside the conductor
          if (cdelpy(ic) < 1.) then
            coeff0(3) = -1./cdelpy(ic)*rdx2
            coeffs(ix,iy,iz,3) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - cvoltpy(ic)/cdelpy(ic)*rdx2*voltswitch
          endif
c         --- the point lower in y is inside the conductor
          if (cdelmy(ic) < 1.) then
            coeff0(4) = -1./cdelmy(ic)*rdx2
            coeffs(ix,iy,iz,4) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - cvoltmy(ic)/cdelmy(ic)*rdx2*voltswitch
          endif
c         --- the point lower in z is inside the conductor
          if (cdelpz(ic) < 1.) then
            coeff0(5) = -1./cdelpz(ic)*rdx2
            coeffs(ix,iy,iz,5) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - cvoltpz(ic)/cdelpz(ic)*rdx2*voltswitch
          endif
c         --- the point lower in z is inside the conductor
          if (cdelmz(ic) < 1.) then
            coeff0(6) = -1./cdelmz(ic)*rdx2
            coeffs(ix,iy,iz,6) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - cvoltmz(ic)/cdelmz(ic)*rdx2*voltswitch
          endif

c         --- Sum contributions to center term.
          coeffs(ix,iy,iz,0) = sum(coeff0)

        endif
      enddo

      return
      end
c=============================================================================
