#include "top.h"
c=============================================================================
c@(#) File CHO.F, version $Revision: 1.16 $, $Date: 2007/04/24 21:37:34 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is the interface to the Chombo/AMR fieldsolver.
c  David P. Grote, LLNL, (510)495-2961
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine cho_setamrgrids(nx, ny, nz, dx, xmmin, ymmin, zmmin,
     &       numLevels, refratio, i, j, k, level, numtags,
     &       bcxlo, bcxhi, bcylo, bcyhi, bczlo, bczhi)
      integer(ISZ) nx, ny, nz, numlevels
      integer(ISZ) i(*), j(*), k(*), refratio(*)
      double precision dx, xmmin, ymmin, zmmin
      integer(ISZ) level(*), numtags
      integer(ISZ) bcxlo, bcxhi, bcylo, bcyhi, bczlo, bczhi
      
c     --- Warning: the following call will not work in parallel!
      call checkconductors(nx,ny,nz,nz,dx,dx,dx,.false.,.false.)

#ifdef WITHCHOMBO
      call setamrgrids(nx, ny, nz, dx, xmmin, ymmin, zmmin,
     &       numLevels, refratio, i, j, k, level, numtags,
     &       bcxlo, bcxhi, bcylo, bcyhi, bczlo, bczhi)
#endif

      return
      end

#ifndef WITHCHOMBO
      subroutine returnphi(phiout)
      double precision phiout(*)
      print*,"no chombo routines"
      return
      end

      subroutine returnphic(phiout)
      double precision phiout(*)
      print*,"no chombo routines"
      return
      end

      subroutine returnrho(rhoout)
      double precision rhoout(*)
      print*,"no chombo routines"
      return
      end
#endif
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
      integer(ISZ):: refratios(cho_nlevels)

#ifdef WITHCHOMBO
      call amrfieldsolve(iwhich)
#endif

#ifdef WITHCHP
      if (iwhich == 0 .or. iwhich == 1) then

c       if (cho_handle /= 0) then
c         --- If the handle is already set, destroy the old one.
c         call chp_destroy(cho_handle,cho_status)
c         if (cho_status > 0) then
c           print*,"Warning from chp_destroy, status = ",cho_status
c         else if (cho_status < 0) then
c           print*,"Fatal error from chp_destroy, status = ",cho_status
c           call kaboom("cho_solve3d: Fatal error")
c           return
c         endif
c       endif

c       --- Create a handle refering to an instance of the Chombo solver
        call chp_create(cho_handle,cho_status)
        if (cho_status > 0) then
          print*,"Warning from chp_create, status = ",cho_status
        else if (cho_status < 0) then
          print*,"Fatal error from chp_create, status = ",cho_status
          call kaboom("cho_solve3d: Fatal error")
          return
        endif

c       --- Initialize data describing the mesh and pass it to Chombo.
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
#ifdef MPIPARALLEL
        subdomains(3,1,:) = izpslave
        subdomains(3,2,:) = izpslave + nzpslave
#else
        subdomains(3,1,:) = 0
        subdomains(3,2,:) = nz
#endif
        refratios = cho_refratio
        call chp_setgridparams(cho_handle,x0,dx,domain,subdomains,
     &                         cho_nlevels,refratios,cho_maxparticlespercell,
     &                         cho_tagsbuffercells,cho_fillratio,
     &                         cho_status)
        if (cho_status > 0) then
          print*,"Warning from chp_setgridparams, status = ",cho_status
        else if (cho_status < 0) then
          print*,"Fatal error from chp_setgridparams, status = ",cho_status
          call kaboom("cho_solve3d: Fatal error")
          return
        endif

c       --- Setup arrays describing the boundary conditions and pass them to
c       --- Chombo.
        if (l2symtry) then
          cho_bcflags(2,1) = 1
          cho_bcvals(2,1) = 0.
        elseif (l4symtry) then
          cho_bcflags(1,1) = 1
          cho_bcflags(2,1) = 1
          cho_bcvals(1,1) = 0.
          cho_bcvals(2,1) = 0.
        endif
        call chp_setsolveparams(cho_handle,cho_tol,cho_bcflags,cho_bcvals,
     &                          cho_chargeinterptype,cho_debug,cho_status)
        if (cho_status > 0) then
          print*,"Warning from chp_solverparams, status = ",cho_status
        else if (cho_status < 0) then
          print*,"Fatal error from chp_solverparams, status = ",cho_status
          call kaboom("cho_solve3d: Fatal error")
          return
        endif

      endif

      if (iwhich == 0 .or. iwhich == -1) then
c       --- Make the call to Chombo to do the field solve.
        call chp_solve(cho_handle,cho_status)
        if (cho_status > 0) then
          print*,"Warning from chp_solve, status = ",cho_status
        else if (cho_status < 0) then
          print*,"Fatal error from chp_solve, status = ",cho_status
          call kaboom("cho_solve3d: Fatal error")
          return
        endif
      endif

#endif

      return
      end
c=============================================================================
      subroutine cho_setrho3d(np,xp,yp,zp,zgrid,sq,sw,js,ip)
      use Constant
      use CHOHandle
      integer(ISZ):: np
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: zgrid,sq,sw
      integer(ISZ):: js,ip

c Interface to the chp_putparticles routine.

#ifdef WITHCHP

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

      call chp_putparticles(cho_handle,cho_bunchid,sq*sw/eps0,np,1,xp,yp,zp,
     &                      cho_status)
      if (cho_status > 0) then
        print*,"Warning from chp_putparticles, status = ",cho_status
      else if (cho_status < 0) then
        print*,"Fatal error from chp_putparticles, status = ",cho_status
        call kaboom("chp_putparticles: Fatal error")
        return
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

#ifdef WITHCHP

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
      if (cho_status > 0) then
        print*,"Warning from chp_getefield, status = ",cho_status
      else if (cho_status < 0) then
        print*,"Fatal error from chp_getefield, status = ",cho_status
        call kaboom("cho_gete3d: Fatal error")
        return
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

#ifdef WITHCHP

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
      if (cho_status > 0) then
        print*,"Warning from chp_getphi, status = ",cho_status
      else if (cho_status < 0) then
        print*,"Fatal error from chp_getphi, status = ",cho_status
        call kaboom("cho_getphi: Fatal error")
        return
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

#ifdef WITHCHP

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

CCCCCCcall chp_getrho(cho_handle,cho_bunchid,np,1,xp,yp,zp,r,cho_status)
      print*,"Error: Not yet implemented!!!"
      if (cho_status > 0) then
        print*,"Warning from chp_getrho, status = ",cho_status
      else if (cho_status < 0) then
        print*,"Fatal error from chp_getrho, status = ",cho_status
        call kaboom("cho_getrho: Fatal error")
        return
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
      use CHOInput
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

c     --- If all nodes are reachable, then just return.
      if (lcho_allnodesreachable) return

c     --- Check for point in conductors.
      do ic=conductors%interior%istart(il),conductors%interior%istart(il+1)-1
        ix = conductors%interior%indx(0,ic)
        iy = conductors%interior%indx(1,ic)
        iz = conductors%interior%indx(2,ic)
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
      call ReachableNodeswork(conductors%evensubgrid,
     &                        mask,xlo,ylo,zlo,xhi,yhi,zhi,level,il)
      call ReachableNodeswork(conductors%oddsubgrid,
     &                        mask,xlo,ylo,zlo,xhi,yhi,zhi,level,il)

      return
      end
c=============================================================================
      subroutine ReachableNodeswork(subgrid,mask,xlo,ylo,zlo,xhi,yhi,zhi,
     &                              level,il)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi,il,level
      integer(ISZ):: mask(xlo:xhi, ylo:yhi, zlo:zhi, 0:6)

c Does the work of determining whether neighboring points are "reachable"
c from the central point. The even or odd conductor points can be passed
c into this routine.

      integer(ISZ):: ic,ix,iy,iz

c     --- loop over points near surface of conductors
      do ic = subgrid%istart(il),subgrid%istart(il+1)-1
 
        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)
 
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then

c         --- the point lower in x is inside the conductor
          if (subgrid%dels(0,ic) < 1.) mask(ix,iy,iz,2) = 0

c         --- the point higher in x is inside the conductor
          if (subgrid%dels(1,ic) < 1.) mask(ix,iy,iz,1) = 0

c         --- the point lower in y is inside the conductor
          if (subgrid%dels(2,ic) < 1.) mask(ix,iy,iz,4) = 0

c         --- the point higher in y is inside the conductor
          if (subgrid%dels(3,ic) < 1.) mask(ix,iy,iz,3) = 0

c         --- the point lower in z is inside the conductor
          if (subgrid%dels(4,ic) < 1.) mask(ix,iy,iz,6) = 0

c         --- the point higher in z is inside the conductor
          if (subgrid%dels(5,ic) < 1.) mask(ix,iy,iz,5) = 0

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
 
      integer(ISZ):: level,ic,il,ix,iy,iz
 
c     --- Get the level of coarseness above the finest level.
      level = nint(dx/dxfine)
      il = int(log(1.*level)/log(2.) + 0.5)

c     --- Zero out the mask
      mask = 0

      do ic=conductors%interior%istart(il),conductors%interior%istart(il+1)-1
        ix = conductors%interior%indx(0,ic)
        iy = conductors%interior%indx(1,ic)
        iz = conductors%interior%indx(2,ic)
c       --- Then check if the point is within the requested region.
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then
c         --- If so, set mask to 1.
          mask(ix,iy,iz) = 1
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
      do ic=conductors%interior%istart(il),conductors%interior%istart(il+1)-1
        ix = conductors%interior%indx(0,ic)
        iy = conductors%interior%indx(1,ic)
        iz = conductors%interior%indx(2,ic)
c       --- Then check if the point is within the requested region.
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then
          coeffs(ix,iy,iz,0) = 1.
          coeffs(ix,iy,iz,1:6) = 0.
          coeffs(ix,iy,iz,7) = conductors%interior%volt(ic)
        endif
      enddo

c     --- Get coefficients for even points.
      call NodalCoefficientswork(conductors%evensubgrid,
     &                           coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,level,il)

c     --- Get coefficients for odd points.
      call NodalCoefficientswork(conductors%oddsubgrid,
     &                           coeffs,xlo,ylo,zlo,xhi,yhi,zhi,rdx2,level,il)

      return
      end
c=============================================================================
      subroutine NodalCoefficientswork(subgrid,coeffs,xlo,ylo,zlo,xhi,yhi,zhi,
     &                                 rdx2,level,il)
      use ConductorSubGridTypemodule
      type(ConductorSubGridType):: subgrid
      integer(ISZ):: xlo,ylo,zlo,xhi,yhi,zhi
      real(kind=8):: coeffs(xlo:xhi,ylo:yhi,zlo:zhi,0:7)
      real(kind=8):: rdx2
      integer(ISZ):: level,il

c Does the work of setting the coefficients. The even or odd conductor
c points can be passed into this routine.

      integer(ISZ):: ic,ix,iy,iz
      real(kind=8):: coeff0(6),voltswitch
      real(kind=8),pointer:: dels(:,:),volt(:,:)

      dels => subgrid%dels
      volt => subgrid%volt

c     --- Voltages are only applied at level == 1, the finest level.
c     --- Otherwise, the voltages are set to zero, meaning that the
c     --- residual is zero
c     if (level == 1) then
        voltswitch = 1.
c     else
c       voltswitch = 0.
c     endif

c     --- loop over points near surface of conductors
      do ic = subgrid%istart(il),subgrid%istart(il+1)-1
 
        ix = subgrid%indx(0,ic)
        iy = subgrid%indx(1,ic)
        iz = subgrid%indx(2,ic)
 
        if (xlo <= ix .and. ix <= xhi .and.
     &      ylo <= iy .and. iy <= yhi .and.
     &      zlo <= iz .and. iz <= zhi) then
 
c         --- Contributions of terms with no conductor boundary nearby.
          coeff0 = -rdx2

c         --- the point higher in x is inside the conductor
          if (dels(1,ic) < 1.) then
            coeff0(1) = -1./dels(1,ic)*rdx2
            coeffs(ix,iy,iz,1) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - volt(1,ic)/dels(1,ic)*rdx2*voltswitch
          endif
c         --- the point lower in x is inside the conductor
          if (dels(0,ic) < 1.) then
            coeff0(2) = -1./dels(0,ic)*rdx2
            coeffs(ix,iy,iz,2) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - volt(0,ic)/dels(0,ic)*rdx2*voltswitch
          endif
c         --- the point higher in y is inside the conductor
          if (dels(3,ic) < 1.) then
            coeff0(3) = -1./dels(3,ic)*rdx2
            coeffs(ix,iy,iz,3) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - volt(3,ic)/dels(3,ic)*rdx2*voltswitch
          endif
c         --- the point lower in y is inside the conductor
          if (dels(2,ic) < 1.) then
            coeff0(4) = -1./dels(2,ic)*rdx2
            coeffs(ix,iy,iz,4) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - volt(2,ic)/dels(2,ic)*rdx2*voltswitch
          endif
c         --- the point higher in z is inside the conductor
          if (dels(5,ic) < 1.) then
            coeff0(5) = -1./dels(5,ic)*rdx2
            coeffs(ix,iy,iz,5) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - volt(5,ic)/dels(5,ic)*rdx2*voltswitch
          endif
c         --- the point lower in z is inside the conductor
          if (dels(4,ic) < 1.) then
            coeff0(6) = -1./dels(4,ic)*rdx2
            coeffs(ix,iy,iz,6) = 0.
            coeffs(ix,iy,iz,7)=coeffs(ix,iy,iz,7) - volt(4,ic)/dels(4,ic)*rdx2*voltswitch
          endif

c         --- Sum contributions to center term.
          coeffs(ix,iy,iz,0) = sum(coeff0)

        endif
      enddo

      return
      end
c=============================================================================
