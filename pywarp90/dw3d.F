#include "top.h"
c=============================================================================
c@(#) File DW3D.F, version $Revision: 3.155 $, $Date: 2006/10/03 23:44:07 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is part of the package W3D of the code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  This file contains the diagnostic and plotting routines.
c  Alex Friedman,  LLNL, (510) 422-0827
c  David P. Grote, LLNL, (510) 423-7194
c=============================================================================
      subroutine prntpa3d(lprntpara)
      use Beam_acc
      use InMesh3d
      use Fields3d
      use InGen
      use Lattice
      use Picglb3d
      use Io
      use Ch_var
      logical(ISZ):: lprntpara

c  Prints out various parameters to a plot frame and an output file or tty

      character(80):: textline
      real(kind=8):: xrbend,xbbend,xbendlen,xstralen,xz0bend

c Exit now if parameter are not to be printed
      if (.not. lprntpara) return

c     --- Call script version of this routine.
      call execuser("printparameters3d")

   20 format(1x,a,1pe11.4,a)
   30 format(1x,a,i8,a)

c  Gather up various lattice qtys for printout

      if (nbend >= 1) then
         xrbend = bendrc(1)
         if (ndipo >= 1) then
           xbbend = dipoby(1)
         else
           xbbend = 0.
         endif
         xbendlen = bendze(1) - bendzs(1)
         xstralen = bendzs(1) - bendze(0)
         xz0bend = bendzs(1)
      else
         xrbend = 0.
         xbbend = 0.
         xbendlen = 0.
         xstralen = 0.
         xz0bend = 0.
      endif

c  Write to plot frame
      write (textline,30) "Number of grid points in x = ",nx," "
      call remark(trim(textline))
      write (textline,30) "Number of grid points in y = ",ny," "
      call remark(trim(textline))
      write (textline,30) "Number of grid points in z = ",nzfull," "
      call remark(trim(textline))
      write (textline,20) "Grid spacing in x = ",dx," m"
      call remark(trim(textline))
      write (textline,20) "Grid spacing in y = ",dy," m"
      call remark(trim(textline))
      write (textline,20) "Grid spacing in z = ",dz," m"
      call remark(trim(textline))
      write (textline,20) "Bend radius = ",xrbend," m"
      call remark(trim(textline))
      write (textline,20) "Bending field = ",xbbend," T"
      call remark(trim(textline))
      write (textline,20) "Bend length = ",xbendlen," m"
      call remark(trim(textline))
      write (textline,20) "Straight section length = ",xstralen," m"
      call remark(trim(textline))
      write (textline,20) "Z at start of first bend = ",xz0bend," m"
      call remark(trim(textline))

c  Write to text output file
      if (warpout > -1) then
        call edit (warpout, "nx")
        call edit (warpout, "ny")
        call edit (warpout, "nz")
        call edit (warpout, "dx")
        call edit (warpout, "dy")
        call edit (warpout, "dz")
      endif

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine rhodia
      use InGen3d
      if(solvergeom==XYZgeom) then
        call rhodia3d
      end if
      return
      end
c=============================================================================
      subroutine gtlchg
      use InGen3d
      if(solvergeom==XYZgeom) then
        call gtlchg3d
      else if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
        call gtlchgrz
      end if
      return
      end
c=============================================================================
      subroutine srhoax
      use InGen3d
      if(solvergeom==XYZgeom) then
        call srhoax3d
      end if
      return
      end
c=============================================================================
      subroutine getese
      use InGen3d
      if(solvergeom==XYZgeom) then
        call getese3d
      end if
      return
      end
c=============================================================================
      subroutine sphiax
      use InGen3d
      if(solvergeom==XYZgeom) then
        call sphiax3d
      end if
      return
      end
c=============================================================================
      subroutine  sezax
      use InGen3d
      if(solvergeom==XYZgeom) then
        call sezax3d
      end if
      return
      end
c=============================================================================
      subroutine gtlchg3d
      use Subtimers3d
      use Picglb
      use InGen3d
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Z_arrays
      use InDiag3d

c  Calculates the line charge density from rho.

      real(kind=8):: zz,zg,wzg,dzi
      integer(ISZ):: ix,iy,iz,izg
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      if (.not. lgtlchg3d .or. .not. associated(rho)) return

      dzi = 1./dz

c     conversion factor to go from grid frame to beam frame
      zz = zgrid + zmmin - zzmin - zbeam

c First, sum up the charge density in each axial plane.

!$OMP PARALLEL DO
      do iz=0,nz

        if (l2symtry) then
c         --- 2-fold: account for only one side of rho being stored
          scrtch(iz,0) = sum(rho(:,0,iz)) + 2.*sum(rho(:,1:,iz))

        elseif (l4symtry) then
c         --- 4-fold: account for only one quadrant of rho being stored
          scrtch(iz,0) = rho(0,0,iz) + 4.*sum(rho(1:,1:,iz)) +
     &                   2.*sum(rho(1:,0,iz)) + 2.*sum(rho(0,1:,iz))

        else
c         --- no symmetries...
          scrtch(iz,0) = sum(rho(:,:,iz))

        endif

      enddo
!$OMP END PARALLEL DO

c Then, linearly interpolate the data into the beam frame array, linechg.

      do iz=0,nzzarr

        zg = (iz*dzz - zz)*dzi
        izg = int(zg + nz) - nz
        wzg = zg - izg

        if (0 <= izg .and. izg <= nz) then
          linechg(iz) = scrtch(izg,0)*dx*dy*(1. - wzg)
        else
          linechg(iz) = 0.
        endif
        if (0 <= izg+1 .and. izg+1 <= nz) then
          linechg(iz) = linechg(iz) + scrtch(izg+1,0)*dx*dy*wzg
        endif

      enddo

#ifdef MPIPARALLEL
      call parallelnonzerorealarray(linechg,nzzarr+1)
#endif

!$OMP MASTER
      if (lw3dtimesubs) timegtlchg3d = timegtlchg3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine rhodia3d
      use InGen3d
      use InDiag
      use InDiag3d
      use InMesh3d
      use Picglb
      use Picglb3d
      use Fields3d
      use Win_Moments
      use Z_Moments

      integer(ISZ):: iwin,ix,iy,iz,ii
      real(kind=8):: zwin,z0,z1,ww

c  Sets rho window diagnostics (they need data from the 3d mesh)

      if (.not. lrhodia3d .or. .not. associated(rho)) return

!$OMP PARALLEL DO PRIVATE(zwin,iz,z1,z0)
      do iwin = 0,NWINDOWS
        zwin = 0.5 * (zwindows(1,iwin) + zwindows(2,iwin))
        iz = (zwin + zbeam - zgrid - zmmin)/dz
        if (zwindows(1,iwin).ne.zwindows(2,iwin) .and.
     &      0 <= iz .and. iz < nz) then
          z1 = (zwin + zbeam - zgrid - zmmin)/dz - iz
          z0 = 1. - z1
          rhomid(iwin) = z0*rho(ix_axis,iy_axis,iz) +
     &                   z1*rho(ix_axis,iy_axis,iz+1)
          rhomax(iwin) = 0.
          do iy = 0,ny
            do ix = 0,nx
               rhomax(iwin) = max(rhomax(iwin),
     &                              z0*rho(ix,iy,iz) + z1*rho(ix,iy,iz+1))
            enddo
          enddo
        endif
      enddo
!$OMP END PARALLEL DO

c  Sets rho grid diagnostics (they need data from the 3d mesh)
      if (ifzmmnt > 0) then

c       --- rho on axis
        do iz = 0,nzmmnt
           ii = (zmntmesh(iz) + zbeam - zgrid - zmmin)/dz
           ww = (zmntmesh(iz) + zbeam - zgrid - zmmin)/dz - ii
           if (ii >= 0 .and. ii < nz)
     &       rhomidz(iz) = rho(ix_axis,iy_axis,ii)*(1. - ww) +
     &                     rho(ix_axis,iy_axis,ii+1)*ww
        enddo

c       --- max rho at each z (vectorized)
        do iz = 0,nzmmnt
           rhomaxz(iz) = -LARGEPOS
        enddo
!$OMP PARALLEL DO PRIVATE(ii,ww)
        do iz = 0,nzmmnt
          ii = (zmntmesh(iz) + zbeam - zgrid - zmmin)/dz
          ww = (zmntmesh(iz) + zbeam - zgrid - zmmin)/dz - ii
          if (ii >= 0 .and. ii < nz) then
            do iy = 0,ny
              do ix = 0,nx
                rhomaxz(iz) = max(rhomaxz(iz),rho(ix,iy,ii)*(1. - ww) +
     &                                        rho(ix,iy,ii+1)*ww)
              enddo
            enddo
          elseif (ii == nz) then
            do iy = 0,ny
              do ix = 0,nx
                rhomaxz(iz) = max(rhomaxz(iz),rho(ix,iy,nz))
              enddo
            enddo
          endif
        enddo
!$OMP END PARALLEL DO
      endif

      return
      end
c=============================================================================
      subroutine srhoax3d
c  Sets 1d array for the charge density on the axis
      use InGen3d
      use InMesh3d
      use InDiag3d
      use Picglb
      use Picglb3d
      use Fields3d
      use Z_arrays
#ifdef MPIPARALLEL
      use Parallel
#endif

      integer(ISZ):: ix,iy,iz,iz0,izg
      real(kind=8):: zz,wzg,dzi

      if(.not. lsrhoax3d .or. .not. associated(rho)) return

      dzi = 1./dz

      ix = ix_axis
      iy = iy_axis

      rhoax = 0.

#ifdef MPIPARALLEL
      iz0 = izfsslave(my_index)
#else
      iz0 = 0
#endif

      zz = zzmin + zbeam - zgrid - zmminglobal

      if (nz == nzzarr .and. abs(zz) <= dz*1.e-5) then
        do iz=0,nz
          rhoax(iz+iz0) = rho(ix,iy,iz)
        enddo
      else
        do iz=0,nzzarr
          izg = int((iz*dzz + zz)*dzi) - iz0
          if (0 <= izg .and. izg < nz) then
            wzg = (iz*dzz + zz)*dzi - iz0 - izg
            rhoax(iz) = rho(ix,iy,izg)*(1. - wzg) + rho(ix,iy,izg+1)*wzg
          else if (izg == nz) then
            rhoax(iz) = rho(ix,iy,izg)
          endif
        enddo
      endif

#ifdef MPIPARALLEL
      call parallelnonzerorealarray(rhoax,nzzarr+1)
#endif

      return
      end
c=============================================================================
      subroutine getese3d
      use Subtimers3d
      use Beam_acc
      use Picglb
      use Picglb3d
      use InMesh3d
      use Fields3d
      use Moments
      use InPart
      use InGen
      use InGen3d
      use Z_arrays
      use InDiag3d

c Calculate the electrostatic energy, sum of rho*phi*dx*dy*dz/2
c Includes external fields: Uniform focusing and eears for cigars
c This does it in a partially vectorized manner.

      integer(ISZ):: ix,iy,iz,izb
      real(kind=8):: zm,zz,phiextun,wzb
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      if (.not. lgetese3d .or.
     &    .not. associated(rho) .or. .not. associated(phi)) then
!$OMP MASTER
        if (lw3dtimesubs) timegetese3d = timegetese3d + wtime() - substarttime
!$OMP END MASTER
        return
      endif

      ese = 0.
      zm = (zimax - zimin)*(1. - straight)*0.5

c     conversion factor to go from grid frame to beam frame
      zz = zmmin + zgrid - zbeam

      do iy=0,ny
        do ix=0,nx
          scrtch(ix,iy) = 0.0
        enddo
      enddo

c     --- Loop order switched for the OpenMP parallel version - this
c     --- avoids requiring the accumulation into scrtch to be put
c     --- into a critical section.
!$OMP PARALLEL
!$OMP DO PRIVATE(phiextun,izb,wzb,ix,iy,iz)
#ifdef _OPENMP
      do iy=0,ny
#endif
      do iz=izfsmin,izfsmax-1
        phiextun = 0.
        if (ifeears == 1) then
          if (iz*dz+zz < -zm) phiextun=eears*.5*(iz*dz+zz+zm)**2
          if (iz*dz+zz >  zm) phiextun=eears*.5*(iz*dz+zz-zm)**2
        elseif (ifeears == 2) then
c         --- This is not quite right for the parallel version since
c         --- the longitudinal extent of eearsofz (the same as the extent
c         --- of the particles) will not necessarily overlap the extent of
c         --- the field solver.
          izb = (iz*dz+zz-zzmin)*dzzi
          wzb = (iz*dz+zz-zzmin)*dzzi - izb
          if (0 <= iz .and. iz < nzzarr)
     &      phiextun = eearsofz(izb)*(1. - wzb) + eearsofz(izb+1)*wzb
        endif
#ifndef _OPENMP
        do iy=0,ny
#endif
          do ix=0,nx
            scrtch(ix,iy) = scrtch(ix,iy) + rho(ix,iy,iz)*(phi(ix,iy,iz)
     &                       - dedr*.5*((ix*dx+xmmin)**2 + (iy*dy+ymmin)**2)
     &                       - dexdx*.5*(ix*dx+xmmin)**2 
     &                       - deydy*.5*(iy*dy+ymmin)**2
     &                       + phiextun)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

c     --- Apply appropriate multiplication factors for symmetries
      if (l4symtry) then
c       --- Center point is unchanged
c       --- Edges are doubled
c       --- Bulk is quadrupoled
        scrtch(0,1:ny) = 2.*scrtch(0,1:ny)
        scrtch(1:nx,0) = 2.*scrtch(1:nx,0)
        scrtch(1:nx,1:ny) = 4.*scrtch(1:nx,1:ny)
      else if (l2symtry) then
c       --- Edge is unchanged
c       --- Bulk is doubled
        scrtch(1:nx,0:ny) = 2.*scrtch(1:nx,0:ny)
      endif

      do iy=1,ny
        do ix=0,nx
          scrtch(ix,0) = scrtch(ix,0) + scrtch(ix,iy)
        enddo
      enddo

      do ix=0,nx
        ese = ese + scrtch(ix,0)
      enddo

      ese = dx*dy*dz*0.5*ese

!$OMP MASTER
      if (lw3dtimesubs) timegetese3d = timegetese3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine sphiax3d
c  Sets 1d array for the potential on the axis
      use InGen3d
      use InMesh3d
      use InDiag3d
      use Z_arrays
      use Picglb
      use Picglb3d
      use Fields3d
#ifdef MPIPARALLEL
      use Parallel
#endif

      integer(ISZ):: ix,iy,iz,iz0,izg
      real(kind=8):: zz,wzg,dzi

      if(.not. lsphiax3d .or. .not. associated(phi)) return

      dzi = 1./dz

      ix = ix_axis
      iy = iy_axis
      
      phiax = 0.

#ifdef MPIPARALLEL
      iz0 = izfsslave(my_index)
#else
      iz0 = 0
#endif

      zz = zzmin + zbeam - zgrid - zmminglobal

      if (nzfull == nzzarr .and. abs(zz) <= dz*1.e-5) then
        do iz=0,nz
          phiax(iz+iz0) = phi(ix,iy,iz)
        enddo
      else
        do iz=0,nzzarr
          izg = int((iz*dzz + zz)*dzi) - iz0
          if (0 <= izg .and. izg <= nz) then
            wzg = (iz*dzz + zz)*dzi - iz0  - izg
            phiax(iz) = phi(ix,iy,izg)*(1. - wzg) + phi(ix,iy,izg+1)*wzg
          endif
        enddo
      endif

#ifdef MPIPARALLEL
      call parallelnonzerorealarray(phiax,nzzarr+1)
#endif

      return
      end
c=============================================================================
      subroutine sezax3d
c  Sets 1d array for the space charge E field on the axis
      use InMesh3d
      use InGen3d
      use InDiag3d
      use Picglb
      use Picglb3d
      use Z_arrays
      use Fields3d
#ifdef MPIPARALLEL
      use Parallel
#endif

      integer(ISZ):: ix,iy,iz,iz0,izg
      real(kind=8):: wx,wy,zz,wzg,wx1,wy1,dzi

      if(.not. lsezax3d .or. .not. associated(phi)) return

      dzi = 1./dz

      ix = int(-xmmin/dx)
      wx = int(-xmmin/dx) - ix
      wx1 = 1. - wx
      iy = int(-ymmin/dy)
      wy = int(-ymmin/dy) - iy
      wy1 = 1. - wy

      ezax = 0.

#ifdef MPIPARALLEL
      iz0 = izfsslave(my_index)
#else
      iz0 = 0
#endif

      zz = zzmin + zbeam - zgrid - zmminglobal

      if ((xmmax == -xmmin .and.  ymmax == -ymmin) .or.
     &    (wx == 0. .and. wy == 0.)) then
        do iz=0,nzzarr
          izg = int((iz*dzz + zz)*dzi) - iz0
          if (0 <= izg .and. izg <= nz-1) then
            wzg = (iz*dzz + zz)*dzi - iz0 - izg
            ezax(iz) = ((phi(ix,iy,izg-1) - phi(ix,iy,izg+1))*(1.-wzg) +
     &                  (phi(ix,iy,izg  ) - phi(ix,iy,izg+2))*    wzg )*0.5*dzi
          else if (izg == nz) then
            ezax(iz) = (phi(ix,iy,izg-1) - phi(ix,iy,izg+1))*0.5*dzi
          endif
        enddo
      else
        do iz=0,nzzarr
          izg = int((iz*dzz + zz)*dzi) - iz0
          if (0 <= izg .and. izg <= nz) then
            wzg = (iz*dzz + zz)*dzi - iz0 - izg
        ezax(iz)=((phi(ix  ,iy  ,izg-1)-phi(ix  ,iy  ,izg+1))*wx1*wy1*(1.-wzg)+
     &            (phi(ix+1,iy  ,izg-1)-phi(ix+1,iy  ,izg+1))*wx *wy1*(1.-wzg)+
     &            (phi(ix  ,iy+1,izg-1)-phi(ix  ,iy+1,izg+1))*wx1*wy *(1.-wzg)+
     &            (phi(ix+1,iy+1,izg-1)-phi(ix+1,iy+1,izg+1))*wx *wy *(1.-wzg)+
     &            (phi(ix  ,iy  ,izg  )-phi(ix  ,iy  ,izg+2))*wx1*wy1*    wzg +
     &            (phi(ix+1,iy  ,izg  )-phi(ix+1,iy  ,izg+2))*wx *wy1*    wzg +
     &            (phi(ix  ,iy+1,izg  )-phi(ix  ,iy+1,izg+2))*wx1*wy *    wzg +
     &            (phi(ix+1,iy+1,izg  )-phi(ix+1,iy+1,izg+2))*wx *wy *    wzg )
     &            *0.5*dzi
          else if (izg == nz) then
            ezax(iz)=((phi(ix  ,iy  ,izg-1)-phi(ix  ,iy  ,izg+1))*wx1*wy1 +
     &                (phi(ix+1,iy  ,izg-1)-phi(ix+1,iy  ,izg+1))*wx *wy1 +
     &                (phi(ix  ,iy+1,izg-1)-phi(ix  ,iy+1,izg+1))*wx1*wy  +
     &                (phi(ix+1,iy+1,izg-1)-phi(ix+1,iy+1,izg+1))*wx *wy)
     &                *0.5*dzi
          endif
        enddo
      endif

#ifdef MPIPARALLEL
      call parallelnonzerorealarray(ezax,nzzarr+1)
#endif

      return
      end
c=============================================================================
      subroutine pltfld3d(fld,freqflag)
      use InDiag3d
      use InPltCtl3d
      use Timers
      character(3):: fld
      integer(ISZ):: freqflag

c  Master control for doing field plots


      real(kind=8):: timetemp,wtime
      integer(ISZ):: i

      timetemp = wtime()

      if (lpltfld3d) then
        if (freqflag == ALWAYS) call parsestr("pltfld3d(1,always)")
        if (freqflag == SELDOM) call parsestr("pltfld3d(2,seldom)")
        return
      endif

      plottime = plottime + (wtime() - timetemp)

      return
      end
c=============================================================================
c=============================================================================
c============== Additional phase space diagnostic routines ===================
c=============================================================================
      subroutine divxy(pgroup,iz,ndiv,divx,divy,divvx,divvx2,divvy,divvy2,
     &                 wnpx,wnpy,itask)
      use ParticleGroupmodule
      use InMesh3d
      use InDiag
      use InPart
      use Picglb
      use Picglb3d
      use Fields3d
      type(ParticleGroup):: pgroup
      integer(ISZ):: iz,ndiv,itask
      real(kind=8):: divx(0:ndiv), divy(0:ndiv)
      real(kind=8):: divvx(0:ndiv), divvx2(0:ndiv)
      real(kind=8):: divvy(0:ndiv), divvy2(0:ndiv)
      real(kind=8):: wnpx(0:ndiv), wnpy(0:ndiv)

c Calculates rms Vx versus X and rms Vy versus Y and returns the data
c in the arrays divx and divy.  The arrays divvx, divvx2, divvy, and divvy2
c are scratch arrays that hold average vx, vx^2, vy and vy^2 versus X and Y
c respectively.
c Wnpx and wnpy are also scratch arrays that hold the numbers of points versus
c X and Y.
c The arrays divvx, divvx2, divvy, and divvy2 set up so that this routine can
c be called several times to build up data.  When itask is non-zero, the
c arrays are not cleared out.

      integer(ISZ):: i,ip,ix,iy
      real(kind=8):: ddxi,ddyi,wx,wy

c Clear everything if requested
      if (itask == 0) then
        do i=0,ndiv
          divvx(i) = 0.
          divvy(i) = 0.
          divvx2(i) = 0.
          divvx2(i) = 0.
          wnpx(i) = 0.
          wnpy(i) = 0.
        enddo
      endif

c Set inverse grid cell lengths
      ddxi = ndiv/(xplmax - xplmin)
      ddyi = ndiv/(yplmax - yplmin)

c loop through particles
      do ip=pgroup%ins(1),pgroup%ins(1)+pgroup%nps(1)-1

c if particle within z grid cell and not lost, accumulate
        if (zmesh(iz)+zbeam < pgroup%zp(ip) .and.
     &                           pgroup%zp(ip) < zmesh(iz+1)+zbeam
     &      .and. pgroup%uzp(ip) /= 0.) then
          ix = (pgroup%xp(ip) - xplmin)*ddxi
          iy = (pgroup%yp(ip) - yplmin)*ddyi
c         if within x grid, accumulate
          if (ix >= 0 .and. ix < ndiv) then
            wx = (pgroup%xp(ip) - xplmin)*ddxi - ix
            wnpx(ix) = wnpx(ix) + 1. - wx
            wnpx(ix+1) = wnpx(ix+1) + wx
            divvx(ix+1)  = divvx(ix+1) +   pgroup%uxp(ip)/pgroup%uzp(ip)*wx
            divvx2(ix+1) = divvx2(ix+1) + (pgroup%uxp(ip)/pgroup%uzp(ip))**2*wx
            divvx(ix)  = divvx(ix) +   pgroup%uxp(ip)/pgroup%uzp(ip)*(1. - wx)
            divvx2(ix) = divvx2(ix) + (pgroup%uxp(ip)/pgroup%uzp(ip))**2*(1.-wx)
          endif
c         if within y grid, accumulate
          if (iy >= 0 .and. iy < ndiv) then
            wy = (pgroup%yp(ip) - yplmin)*ddyi - iy
            wnpy(iy) = wnpy(iy) + 1. - wy
            wnpy(iy+1) = wnpy(iy+1) + wy
            divvy(iy)  = divvy(iy) +   pgroup%uyp(ip)/pgroup%uzp(ip)*(1.-wy)
            divvy2(iy) = divvy2(iy) + (pgroup%uyp(ip)/pgroup%uzp(ip))**2*(1.-wy)
            divvy(iy+1)  = divvy(iy+1) +   pgroup%uyp(ip)/pgroup%uzp(ip)*wy
            divvy2(iy+1) = divvy2(iy+1) + (pgroup%uyp(ip)/pgroup%uzp(ip))**2*wy
          endif
        endif
      enddo

c finally, calculate RMS quantities.
      do i=0,ndiv
        divx(i) = sqrt(max(1.e-50,divvx2(i)/dvnz(wnpx(i)) -
     &                            divvx(i)**2/dvnz(wnpx(i))))
        divy(i) = sqrt(max(1.e-50,divvy2(i)/dvnz(wnpy(i)) -
     &                            divvy(i)**2/dvnz(wnpy(i))))
      enddo

      return
      end
c=============================================================================
      subroutine multpole(lmod,nlmod,irpowmx,lcosex,lsinex,aper,xcen,ycen,
     &                    nmult,nres,tol)
      use InMesh3d
      use Multipole
      integer(ISZ):: nlmod, irpowmx
      real(kind=8):: aper, xcen, ycen 
      integer(ISZ):: lmod(nlmod) 
      logical(ISZ):: lcosex, lsinex 
      integer(ISZ):: nmult,nres
      real(kind=8):: tol

c Calculates the multipole moments as an axial function of z from the 
c electrostatic potential using a linear least squares fit to extract the 
c radial dependence.  Both the multipole expansion and the numerical 
c method employed is described in detail in Appendix A and B of HIF note 
c number 94-3 by Steve Lund.    
c
c Inputs:
c
c    nlmod           number of azimuthal mode numbers to calculate moments 
c    lmod(1:nlmod)   specific  azimuthal mode numbers to calculate moments 
c    irpowmx         maximum radial power cutoff of moments 
c    lcosex          log. flag -- calc. multipoles, cos(theta) expansion terms
c    lsinex          log. flag -- calc. multipoles, sin(theta) expansion terms
c    aper            aperture of moments calculation (for no apeture set to 1.)
c    xcen            transverse x-coordinate moments are calculated about 
c    ycen            transverse y-coordinate moments are calculated about 
c    nmult           multipling factor on the number of radial points to fit
c    nres            scaling for minimum radius to fit
c    tol             tolerance for fitting
c
c Outputs (passed via common block) interpreter variables:
c [ All but nmom are actually calculated in the work routine multpwrk 
c   called in multpole. ]
c  
c    nmom                  number of terms in multipole expansion 
c    lazi(1:nmom)          vector of azimuthal mode numbers of moments 
c    irpow(1:nmom)         vector of radial powers of moments 
c    rmomcex(1:nmom,0:nz)  multipoles, cos(theta) expansion terms on axial grid
c    rmomsex(1:nmom,0:nz)  multipoles, sin(theta) expansion terms on axial grid
c 
c
c The multipole expansion is of the following schematic form 
c
c                                               
c
c                        
c                    /--  /--
c phi(r,theta,iz) =  \--  \--   (  rmomcex(iz)_l,k cos(l theta) + 
c                       l   k       
c                                  rmomsex(iz)_l,k sin(l theta)   )(r/aper)^k
c
c
c   where /-- denote constrained sums.  The azimuthal mode number sum l 
c         \--
c   is over the nlmod mode numbers specified by the vector lmod(1:nlmod), 
c   and the sum over radial powers k ranges from l to irpowmx in steps of 
c   2 (as a consequence of Maxwell's equations).  
c         


      integer(ISZ):: l, k

c     --- Set number of z grid points.
      nzmom = nz

c     --- find the total number of moments nmom to be calculated and 
c     --- resize dynamic arrays in the common block Multipole 
      nmom = 0 
      do l = 1,nlmod
        do k = lmod(l),irpowmx,2
          nmom = nmom + 1
        enddo 
      enddo
      call gchange("Multipole",0)

c     --- call work routine to do the actual multipole calculation.  This 
c     --- structure is necessary to avoid addressing problems associated 
c     --- with resizing multidimensional arrays
      call  multpwrk(lmod,nlmod,irpowmx,lcosex,lsinex,aper,xcen,ycen,
     &               nmult,nres,tol)

      return 
      end 
c=============================================================================
      subroutine multpwrk(lmod,nlmod,irpowmx,lcosex,lsinex,aper,xcen,ycen,
     &                    nmult,nres,tol)
      use Constant
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Fields3d
      use Multipole
      integer(ISZ):: nlmod, irpowmx
      real(kind=8):: aper, xcen, ycen 
      integer(ISZ):: lmod(nlmod) 
      logical(ISZ):: lcosex, lsinex 
      integer(ISZ):: nmult,nres
      real(kind=8):: tol

c Work routine for the calculation of multipole moments as an axial 
c function of z.  See comments in routine multpole for variable and 
c general method descriptions.   
c Standard values for nmult=2, nres=10, tol=1.e-5.

     
      integer(ISZ):: nptsmx, nfitmx 
      parameter(nfitmx = 20, nptsmx = 500)
     
      integer(ISZ):: k, l, iz, imom, lfund, nfit, npts, ngmin 
      integer(ISZ):: np, ith, ipow, ifit, nth, ix, iy    
      real(kind=8):: ds, xsymin, ysymin, dth, xc, yc, phiinter, chisq   
      real(kind=8):: rdat(nptsmx), basis(nptsmx,nfitmx)
      real(kind=8):: usvd(nptsmx,nfitmx), vsvd(nfitmx,nfitmx), wsvd(nfitmx)
      real(kind=8):: tmp(nptsmx) 
      real(kind=8):: azicosint(nptsmx), azisinint(nptsmx) 
      real(kind=8):: rmomcst(nfitmx), rmomsst(nfitmx)

c     --- define constants 
      ds = sqrt((dx*dx+dy*dy)/2.)

c     --- for each moment index i = 1 to nmom, specify azimuthal mode numbers 
c     --- lmod(i) and radial powers irpow(i) consistent with the mode numbers 
c     --- given in lmod(nlmod) and the maximum radial cutoff power irpowmx 
      imom = 0 
      do l = 1,nlmod
        do k = lmod(l),irpowmx,2
          imom = imom + 1
          lazi(imom)  = lmod(l)
          irpow(imom) = k 
        enddo 
      enddo 

c     --- start moments calculation 
      imom = 0
c     --- loop over each moment index with same azimuthal mode number l 
      do lfund = 1,nlmod 
c       --- azimuthal mode number l   
        l = lmod(lfund)
c       --- number of radial power terms nfit used in least-squares fit
        nfit = irpowmx/2 + 1 - (l+mod(irpowmx+1,2))/2
        if (nfit > nfitmx) then
          call remark("routine multpole -- too many radial multipoles, inc. nfitmx")
          call kaboom (0)  
        endif  
c       --- calculate number radial data points npts, min grid index ngmin 
        npts  = nmult*nfit
        if (npts > nptsmx) then 
          call remark("routine multpole -- too many radial data points, inc. nptsmx")
          call kaboom (0)  
        endif 
        ngmin = (nres-1)*l/(2.*pi)
        if (ngmin < 2) ngmin = 2
c       --- check for bound errors
        xsymin = xmmin 
        ysymin = ymmin 
        if (l2symtry) ysymin = -ymmax 
        if (l4symtry) then 
          xsymin = -xmmax 
          ysymin = -ymmax 
        endif  
        if ( (xcen + (npts-1+ngmin)*ds >= xmmax ) .or. 
     &       (xcen - (npts-1+ngmin)*ds <= xsymin)       ) then 
          call remark("routine multpole -- data out of x-grid bound")
          call kaboom (0)
        endif   
        if ( (ycen + (npts-1+ngmin)*ds >= ymmax ) .or. 
     &       (ycen - (npts-1+ngmin)*ds <= ysymin)       ) then 
          call remark("routine multpole -- data out of y-grid bound")
          call kaboom (0) 
        endif   
c       --- calculate npts radial data points and a npts x (kmax-l)/2 basis 
c       --- array of radial powers basis for use in radial least squares fit
        do np = 1,npts
c         --- assign radial data point 
          rdat(np) = (np-1+ngmin)*ds  
c         --- fill (kmax-l)/2 vector column of basis matrix 
          ifit = 0
          do ipow = l,irpowmx,2
            ifit = ifit + 1
            basis(np,ifit) = rdat(np)**ipow 
          enddo 
        enddo 
c       --- calculate singular value decomposition of basis matrix for 
c       --- use in the linear least squares fit, first copy basis to usvd  
        do np = 1,npts 
          do ifit = 1,nfit
            usvd(np,ifit) = basis(np,ifit)  
          enddo 
        enddo 
        call dvdcmp(usvd,npts,nfit,nptsmx,nfitmx,wsvd,vsvd,tmp)
c       --- loop over axial grid index iz, and calculate multipoles 
c       --- with azimuthal mode number l at each index 
        do iz = 0,nz
c         --- azimuthal integrals of the potential at each radial data point 
          do np = 1,npts
c           --- angular increments in integral 
            dth = ds/rdat(np)
            nth = int(2.*pi/dth)
            dth = 2.*pi/nth
c           --- azimuthal integral 
            azicosint(np) = 0.
            azisinint(np) = 0.
            do ith = 0,nth-1
c             --- transverse x and y coordinates to calculate potential  
              xc = xcen + rdat(np)*cos(dth*ith)
              yc = ycen + rdat(np)*sin(dth*ith)
c             --- symmetry operations for 2-fold (l2symtry) and 
c             --- 4-fold (l4symtry) symmetry about the y and x-y axes 
              if (l2symtry) then 
                if ((yc-ymmin) < 0.) yc = 2.*ymmin - yc 
              endif 
              if (l4symtry) then
                if ((xc-xmmin) < 0.) xc = 2.*xmmin - xc 
                if ((yc-ymmin) < 0.) yc = 2.*ymmin - yc 
              endif 
c             --- bottom left grid indices of transverse grid box 
c             --- bounding (xc,yc) 
              ix = int( ( xc - xmmin )/dx ) 
              iy = int( ( yc - ymmin )/dy ) 
c             --- area interpolated potential at (xc,yc) 
              phiinter = 
     &          ( phi(ix  ,iy  ,iz)*
     &              (xmesh(ix+1) - xc       )*(ymesh(iy+1) - yc       ) +
     &            phi(ix+1,iy  ,iz)*
     &              (xc          - xmesh(ix))*(ymesh(iy+1) - yc       ) +
     &            phi(ix+1,iy+1,iz)*
     &              (xc          - xmesh(ix))*(yc          - ymesh(iy)) +
     &            phi(ix  ,iy+1,iz)*
     &              (xmesh(ix+1) - xc       )*(yc          - ymesh(iy))    )/ 
     &          (dx*dy) 
c             --- add contributions for cosine and sine azimuthal integrals 
              azicosint(np) = azicosint(np) + phiinter*cos(l*dth*ith)
              azisinint(np) = azisinint(np) + phiinter*sin(l*dth*ith)
            enddo
c           --- set azimuthal integrals to zero if not calculated 
            if (.not. lcosex) azicosint(np) = 0.
            if (.not. lsinex) azisinint(np) = 0.   
c           --- rescale azimuthal integrals  
            azicosint(np) = (2./real(nth))*azicosint(np)
            azisinint(np) = (2./real(nth))*azisinint(np)
            if (l == 0) then 
              azicosint(np) = 0.5*azicosint(np)
              azisinint(np) = 0.
            endif  
          enddo  
c         --- calculate vectors of moments, rmomcst and rmomsst for cosine 
c         --- and sine expansion terms using linear least squares algorithm 
          call dvdfit(rdat,azicosint,npts,nptsmx,rmomcst,basis,nfit,nfitmx,
     &                usvd,wsvd,vsvd,tmp,tol,chisq)
          call dvdfit(rdat,azisinint,npts,nptsmx,rmomsst,basis,nfit,nfitmx,
     &                usvd,wsvd,vsvd,tmp,tol,chisq)
c         --- load moments calculated in axial grid indexed vectors of 
c         --- multipole moments 
          do ifit = 1,nfit
            rmomcex(imom+ifit,iz) = rmomcst(ifit)*aper**irpow(imom+ifit)
            rmomsex(imom+ifit,iz) = rmomsst(ifit)*aper**irpow(imom+ifit)  
          enddo
        enddo 
        imom = imom + nfit 
      enddo

      return 
      end 
c=============================================================================
