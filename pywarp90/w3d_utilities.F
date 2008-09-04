#include "top.h"
c=============================================================================
c@(#) File w3d_utilities.F, version $Revision: 1.23 $, $Date: 2008/09/04 21:06:06 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package W3D of code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
      subroutine sortparticlesbyindex1(n,indx,x,y,z,uz,nblocks,
     &                                 xout,yout,zout,uzout,pcounts)
      use Subtimers3d
      integer(ISZ):: n,nblocks
      integer(ISZ):: indx(n)
      real(kind=8):: x(n),y(n),z(n),uz(n)
      real(kind=8):: xout(n),yout(n),zout(n),uzout(n)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = indx(ip)
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = indx(ip)
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        uzout(ii) = uz(ip)
        tcounts(ib) = tcounts(ib) + 1
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesortpart_byindex = timesortpart_byindex +
     &                                             wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sortparticlesbyindex1w(n,indx,x,y,z,uz,wfact,nblocks,
     &                                 xout,yout,zout,uzout,wfactout,pcounts)
      use Subtimers3d
      integer(ISZ):: n,nblocks
      integer(ISZ):: indx(n)
      real(kind=8):: x(n),y(n),z(n),uz(n),wfact(n)
      real(kind=8):: xout(n),yout(n),zout(n),uzout(n),wfactout(n)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = indx(ip)
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = indx(ip)
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        uzout(ii) = uz(ip)
        wfactout(ii) = wfact(ip)
        tcounts(ib) = tcounts(ib) + 1
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesortpart_byindex = timesortpart_byindex +
     &                                             wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sortparticlesbyindex2(n,indx,x,y,z,ux,uy,uz,gaminv,nw,wght,
     &                                 nblocks,
     &                                 xout,yout,zout,uxout,uyout,uzout,
     &                                 gaminvout,wghtout,pcounts)
      use Subtimers3d
      integer(ISZ):: n,nw,nblocks
      integer(ISZ):: indx(n)
      real(kind=8):: x(n),y(n),z(n),ux(n),uy(n),uz(n),gaminv(n)
      real(kind=8):: wght(nw)
      real(kind=8):: xout(n),yout(n),zout(n),uxout(n),uyout(n),uzout(n)
      real(kind=8):: gaminvout(n)
      real(kind=8):: wghtout(nw)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = indx(ip)
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = indx(ip)
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        uxout(ii) = ux(ip)
        uyout(ii) = uy(ip)
        uzout(ii) = uz(ip)
        gaminvout(ii) = gaminv(ip)
        if (nw == n) wghtout(ii) = wght(ip)
        tcounts(ib) = tcounts(ib) + 1
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesortpart_byindex = timesortpart_byindex +
     &                                             wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sortparticlesbyindex2w(n,indx,x,y,z,ux,uy,uz,gaminv,wfact,nw,wght,
     &                                 nblocks,
     &                                 xout,yout,zout,uxout,uyout,uzout,
     &                                 gaminvout,wfactout,wghtout,pcounts)
      use Subtimers3d
      integer(ISZ):: n,nw,nblocks
      integer(ISZ):: indx(n)
      real(kind=8):: x(n),y(n),z(n),ux(n),uy(n),uz(n),gaminv(n),wfact(n)
      real(kind=8):: wght(nw)
      real(kind=8):: xout(n),yout(n),zout(n),uxout(n),uyout(n),uzout(n)
      real(kind=8):: gaminvout(n),wfactout(n)
      real(kind=8):: wghtout(nw)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
     
c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = indx(ip)
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = indx(ip)
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        uxout(ii) = ux(ip)
        uyout(ii) = uy(ip)
        uzout(ii) = uz(ip)
        gaminvout(ii) = gaminv(ip)
        wfactout(ii) = wfact(ip)
        if (nw == n) wghtout(ii) = wght(ip)
        tcounts(ib) = tcounts(ib) + 1
      enddo
!$OMP MASTER
      if (lw3dtimesubs) timesortpart_byindex = timesortpart_byindex +
     &                                             wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sortparticlesbyindexgetisort(n,indx,x,y,z,nblocks,
     &                                        xout,yout,zout,isort,pcounts)
      use Subtimers3d
      integer(ISZ):: n,nblocks
      integer(ISZ):: indx(n)
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: xout(n),yout(n),zout(n)
      integer(ISZ):: isort(n)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = indx(ip)
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = indx(ip)
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        isort(ii) = ip - 1
        tcounts(ib) = tcounts(ib) + 1
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesortpart_byindexgetisort =
     &            timesortpart_byindexgetisort + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine getichild(gridnumb,np,x,y,z,ichild,nx,ny,nz,grid,
     &                     xmin,xmax,ymin,ymax,zmin,zmax,zgrid,
     &                     l2symtry,l4symtry)
      use Subtimers3d
      integer(ISZ):: gridnumb,nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np)
      integer(ISZ):: ichild(np)
      integer(ISZ):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax,zgrid
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.
c Same as getgridngp3d, but only gathers the data if it is positive.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: xsymmetryplane,ysymmetryplane
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- These now default to zero, but could be input quantities.
      xsymmetryplane = 0.
      ysymmetryplane = 0.

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      if (ymin == ymax) then
        dyi = 1.
      else
        dyi = ny/(ymax - ymin)
      endif
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np
        if (ichild(ip) .ne. gridnumb) cycle

c       --- find location on grid, taking into account any symmetries.
        if (l4symtry) then
          gx = abs(x(ip) - xsymmetryplane) + xsymmetryplane
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else if (l2symtry) then
          gx = x(ip)
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else
          gx = x(ip)
          gy = y(ip)
        endif

        gz = z(ip) - zgrid

c       --- if not within grid, skip it
        if (gx < xmin .or. gx > xmax .or.
     &      gy < ymin .or. gy > ymax .or.
     &      gz < zmin .or. gz > zmax) cycle

        gx = (gx - xmin)*dxi
        gy = (gy - ymin)*dyi
        gz = (gz - zmin)*dzi

        ix = int(gx)
        iy = int(gy)
        iz = int(gz)
        ichild(ip) = abs(grid(ix,iy,iz))

      enddo

!$OMP MASTER
      if (lw3dtimesubs) timegetichild = timegetichild + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine getichildpositiveonly(gridnumb,np,x,y,z,ichild,nx,ny,nz,grid,
     &                                 xmin,xmax,ymin,ymax,zmin,zmax,zgrid,
     &                                 l2symtry,l4symtry)
      use Subtimers3d
      integer(ISZ):: gridnumb,nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np)
      integer(ISZ):: ichild(np)
      integer(ISZ):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax,zgrid
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.
c Same as getgridngp3d, but only gathers the data if it is positive.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: xsymmetryplane,ysymmetryplane
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- These now default to zero, but could be input quantities.
      xsymmetryplane = 0.
      ysymmetryplane = 0.

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      if (ymin == ymax) then
        dyi = 1.
      else
        dyi = ny/(ymax - ymin)
      endif
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np
        if (ichild(ip) .ne. gridnumb) cycle

c       --- find location on grid, taking into account any symmetries.
        if (l4symtry) then
          gx = abs(x(ip) - xsymmetryplane) + xsymmetryplane
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else if (l2symtry) then
          gx = x(ip)
          gy = abs(y(ip) - ysymmetryplane) + ysymmetryplane
        else
          gx = x(ip)
          gy = y(ip)
        endif

        gz = z(ip) - zgrid

c       --- if not within grid, skip it
        if (gx < xmin .or. gx > xmax .or.
     &      gy < ymin .or. gy > ymax .or.
     &      gz < zmin .or. gz > zmax) cycle

        gx = (gx - xmin)*dxi
        gy = (gy - ymin)*dyi
        gz = (gz - zmin)*dzi

        ix = int(gx)
        iy = int(gy)
        iz = int(gz)
        if (grid(ix,iy,iz) >= 0) ichild(ip) = grid(ix,iy,iz)

      enddo

!$OMP MASTER
      if (lw3dtimesubs) timegetichildpositiveonly = timegetichildpositiveonly +
     &                                              wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine addsortedefield(n,isort,tex,tey,tez,ex,ey,ez)
      use Subtimers3d
      integer(ISZ):: n
      integer(ISZ):: isort(0:n-1)
      real(kind=8),dimension(0:n-1):: tex,tey,tez,ex,ey,ez

      integer(ISZ):: i
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      do i=0,n-1
        ex(isort(i)) = ex(isort(i)) + tex(i)
        ey(isort(i)) = ey(isort(i)) + tey(i)
        ez(isort(i)) = ez(isort(i)) + tez(i)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timeaddsortedefield = timeaddsortedefield +
     &                                        wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine getabsgrad3d(nx,ny,nz,f,gr,dx,dy,dz)
      use Subtimers3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: f(0:nx,0:ny,0:nz)
      real(kind=8):: gr(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz

      integer(ISZ):: ix,iy,iz
      real(kind=8):: dxi,dyi,dzi
      integer(ISZ):: ixm1,ixp1,iym1,iyp1,izm1,izp1
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dxi = 1./dx**2
      dyi = 1./dy**2
      dzi = 1./dz**2

      do iz=0,nz
        do iy=0,ny
          do ix=0,nx
            ixm1 = ix - 1
            ixp1 = ix + 1
            if (ix == 0) ixm1 = 1
            if (ix == nx) ixp1 = nx-1
            iym1 = iy - 1
            iyp1 = iy + 1
            if (iy == 0) iym1 = 1
            if (iy == ny) iyp1 = ny-1
            izm1 = iz - 1
            izp1 = iz + 1
            if (iz == 0) izm1 = 1
            if (iz == nz) izp1 = nz-1

c            gr(ix,iy,iz) = (abs(f(ixp1,iy  ,iz  )-f(ix  ,iy  ,iz  )) +
c     &                      abs(f(ix  ,iy  ,iz  )-f(ixm1,iy  ,iz  )))*dxi +
c     &                     (abs(f(ix  ,iyp1,iz  )-f(ix  ,iy  ,iz  )) +
c     &                      abs(f(ix  ,iy  ,iz  )-f(ix  ,iym1,iz  )))*dyi +
c     &                     (abs(f(ix  ,iy  ,izp1)-f(ix  ,iy  ,iz  )) +
c     &                      abs(f(ix  ,iy  ,iz  )-f(ix  ,iy  ,izm1)))*dzi
            gr(ix,iy,iz) = sqrt( max( abs(f(ixp1,iy  ,iz  )-f(ix  ,iy  ,iz  )),
     &                                abs(f(ix  ,iy  ,iz  )-f(ixm1,iy  ,iz  )))**2*dxi +
     &                           max( abs(f(ix  ,iyp1,iz  )-f(ix  ,iy  ,iz  )),
     &                                abs(f(ix  ,iy  ,iz  )-f(ix  ,iym1,iz  )))**2*dyi +
     &                           max( abs(f(ix  ,iy  ,izp1)-f(ix  ,iy  ,iz  )),
     &                                abs(f(ix  ,iy  ,iz  )-f(ix  ,iy  ,izm1)))**2*dzi)
          enddo
        enddo
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timegetabsgrad = timegetabsgrad + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      logical(ISZ) function arecoordinatesconsist(nx,dx,xmin,xmax,dxi)
      integer(ISZ):: nx
      real(kind=8):: dx,xmin,xmax,dxi

c For the Gridndtype (or any grid coordinate specifications), makes sure
c that the input is consistent. At most one can be unspecified. If all are
c set, the values must be consistent. If one is not specified, it is
c calculated. The grid cell size inverses are always calculated

      integer(ISZ):: defaultcount

      defaultcount = 0
      if (nx == -1) defaultcount = defaultcount + 1
      if (dx == LARGEPOS) defaultcount = defaultcount + 1
      if (xmin == +LARGEPOS) defaultcount = defaultcount + 1
      if (xmax == -LARGEPOS) defaultcount = defaultcount + 1

      if (defaultcount > 1) then
        arecoordinatesconsist = .false.
        return
      else if (defaultcount == 0) then
        arecoordinatesconsist = (nx == nint((xmax - xmin)/dx))
        return
      else
        if (nx == -1) nx = nint((xmax - xmin)/dx)
        if (dx == LARGEPOS) dx = (xmax - xmin)/nx
        if (xmin == +LARGEPOS) xmin = xmax - nx*dx
        if (xmax == -LARGEPOS) xmax = xmin + nx*dx
        dxi = 1./dx
        arecoordinatesconsist = .true.
        return
      endif

      return
      end
c=============================================================================
      subroutine setupgrid3dtype(grid,check)
      use Grid3dtypemodule
      type(Grid3dtype):: grid
      logical(ISZ):: check

c This checks if the input coordinates are consistent, and if so, allocates
c array.

      logical(ISZ):: arecoordinatesconsist

      check = arecoordinatesconsist(grid%nx,grid%dx,grid%xmin,grid%xmax,
     &                              grid%dxi)
      if (.not. check) return
      check = arecoordinatesconsist(grid%ny,grid%dy,grid%ymin,grid%ymax,
     &                              grid%dyi)
      if (.not. check) return
      check = arecoordinatesconsist(grid%nz,grid%dz,grid%zmin,grid%zmax,
     &                              grid%dzi)
      if (.not. check) return

      call Grid3dtypechange(grid)

      return
      end
c=============================================================================
      subroutine setupgrid2dtype(grid,check)
      use Grid2dtypemodule
      type(Grid2dtype):: grid
      logical(ISZ):: check

c This checks if the input coordinates are consistent, and if so, allocates
c array.

      logical(ISZ):: arecoordinatesconsist

      check = arecoordinatesconsist(grid%nx,grid%dx,grid%xmin,grid%xmax,
     &                              grid%dxi)
      if (.not. check) return
      check = arecoordinatesconsist(grid%ny,grid%dy,grid%ymin,grid%ymax,
     &                              grid%dyi)
      if (.not. check) return

      call Grid2dtypechange(grid)

      return
      end
c=============================================================================
c=============================================================================
c Some routines used by the parallel code that need to be accessible
c from python
c=============================================================================
c=============================================================================
      subroutine domaindecomposefields(nz,nslaves,lfsautodecomp,
     &                                 izfsslave,nzfsslave,overlap)
      use InGen, Only: fstype
      integer(ISZ):: nz,nslaves
      logical(ISZ):: lfsautodecomp
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: overlap

      real(kind=8):: zperproc,avezpp,ztot
      integer(ISZ):: bestnz
      integer(ISZ):: i

c---------------------------------------------------------------------------
c     --- An overlap of one plane is needed for the FFT field solvers.
c     --- An overlap of two planes is needed for SOR field solver since
c     --- plane 0 of each processor must overlap with a plane which is
c     --- calculated in the neighboring processor.
c     --- As an example, for nz=8 divided among 2 PE's...
c     ---   for overlap=1        0 1 2 3 4 5 6 7 8
c     ---                       |_________|
c     ---                               |_________|
c     ---   for overlap=2        0 1 2 3 4 5 6 7 8
c     ---                       |___________|
c     ---                               |_________|
      if (overlap == 0) then
        if (fstype == 3 .or. fstype == 7 .or. fstype == 8 .or. fstype == 13) then
          overlap = 2
        else
          overlap = 1
        endif
      endif

c---------------------------------------------------------------------------
c     --- Do the domain decomposition for the field solver.
c     --- The domain decomposition by default is done so that each processor
c     --- has nearly the same number of z planes. If lautodecomp is false and
c     --- fstype == 3, then the decomposition is supplied by the user (input
c     --- through nzfsslave).
      if (lfsautodecomp .or.
     &    (fstype /= 3 .and. fstype /= 7 .and. fstype /= 12 .and.
     &     fstype /= 13)) then

c       --- Calculate average number of z planes per processor, including the
c       --- extra space for overlap.  For FFT, the minimum number of planes
c       --- allowable is 1 and for PSOR, the minimum number of planes allowable
c       --- is 2.
        zperproc = (nz + (nslaves - 1.)*(overlap - 1))/nslaves
        if (zperproc < overlap) zperproc = overlap

c       --- bestnz is zperproc rounded down.
        bestnz = int(zperproc)
     
c       --- First processor is easy
        nzfsslave(0) = bestnz
        izfsslave(0) = 0
        ztot = nzfsslave(0)

c       --- loop over processors until used all processors or have assigned
c       --- all of grid.
        i = 0
        do while (i < nslaves-1 .and.
     &            izfsslave(i)+nzfsslave(i)-overlap+1 < nz)
          i = i + 1

c         --- This processor starts at the end of the region covered
c         --- by the previous processor, overlapping it by the value of overlap.
          izfsslave(i) = izfsslave(i-1) + nzfsslave(i-1) + 1 - overlap

c         --- The number of z planes given to this processor is first assumed to
c         --- be bestnz.  If this gives an average number of z planes per
c         --- processor that is less then zperproc, it is increased by 1.
          nzfsslave(i) = bestnz
          ztot = ztot + nzfsslave(i)
          avezpp = ztot/(i+1)
          if (avezpp < zperproc) then
            nzfsslave(i) = nzfsslave(i) + 1
            ztot = ztot + 1
          endif

c         --- Check if this region extends past the end of the grid.  If so,
c         --- recalculate nzfsslave(i).
          if (izfsslave(i) + nzfsslave(i) > nz) then
            nzfsslave(i) = nz - izfsslave(i)
c           --- If nzfsslave(i) is less than 3 then skip this
c           --- processor and give remaining zones to previous processor.
c           if (nzfsslave(i) < 3) then
c             i = i - 1
c             nzfsslave(i) = nz - izfsslave(i)
c           endif
          endif

        enddo

c       --- Save the number of processors that have part of the grid assigned
c       --- to them.
        nslaves = i+1

      else

c       --- nzfsslave is assumed to be input by the user and is assumed to
c       --- not include the overlap.

c       --- First check to make sure that all values are > 0.
        do i=0,nslaves-1
          if (nzfsslave(i) == 0) then
            call kaboom("domaindecomposefields: ERROR: nz for the field solver
     & for each processor must be greater than zero")
            return
          endif
        enddo

c       --- Fill in the array izfsslave, based on the inputted nzfsslave and
c       --- add the overlap to nzfsslave.
        izfsslave(0) = 0
        do i=1,nslaves-1
          izfsslave(i) = izfsslave(i-1) + nzfsslave(i-1)
        enddo
c       --- Note that the last processor has no overlap
        do i=0,nslaves-2
          nzfsslave(i) = nzfsslave(i) + overlap - 1
        enddo
c       --- Get the new value of nz.
        nz = izfsslave(nslaves-1) + nzfsslave(nslaves-1)

      endif

      return
      end
c=============================================================================
      subroutine domaindecomposeparticles(nz,nslaves,izfsslave,nzfsslave,
     &                                    overlap,nzpguard,zmmin,zmmax,dz,zslave,
     &                                    lautodecomp,izpslave,nzpslave,
     &                                    zpslmin,zpslmax)
      use InGen3d,Only: solvergeom,AMRgeom
      integer(ISZ):: nz,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: overlap,nzpguard
      real(kind=8):: zmmin,zmmax,dz
      real(kind=8):: zslave(0:nslaves-1)
      logical(ISZ):: lautodecomp
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)
      real(kind=8):: zpslmin(0:nslaves-1),zpslmax(0:nslaves-1)

      integer(ISZ):: i
      real(kind=8):: sumzslave,zlast


c     --- Now, set the domain decomposition for the particles. This can
c     --- either be input from the user or set the same as the decompostion
c     --- for the field solver.
      if (lautodecomp) then

c       --- Set from field solver decompostion
        do i=0,nslaves-1
          izpslave(i) = izfsslave(i)
          nzpslave(i) = nzfsslave(i) - (overlap - 1)
          if (i == nslaves-1) nzpslave(i) = nzfsslave(i)
c         --- These shouldn't include the guard cells
          zpslmin(i) = izpslave(i)*dz + zmmin
          zpslmax(i) = (izpslave(i)+nzpslave(i))*dz + zmmin
c         --- Now include the guard cells
          izpslave(i) = izpslave(i) - nzpguard
          nzpslave(i) = nzpslave(i) + 2*nzpguard
c         --- Make sure that the processors doesn't have grid cells
c         --- sticking out the ends.
          if (izpslave(i) < 0) then
            nzpslave(i) = nzpslave(i) + izpslave(i)
            izpslave(i) = 0
          endif
          if (izpslave(i) + nzpslave(i) > nz) nzpslave(i) = nz - izpslave(i)
        enddo
c       --- Set zpslmax for the last process so that it is exactly equal to
c       --- zmmax (avoiding problems with roundoff errors).
        zpslmax(nslaves-1) = zmmax

      else
c       --- It is assumed that the user supplied decomposition is specified
c       --- in the array zslave, which is the fractional z-ranges
c       --- of the particles for each processor. It is assumed that zslave
c       --- has been properly normalized so that if the sum of the z-ranges
c       --- covers the whole grid, then the sum is one.

c       --- Get sum of zslave to allow proper scaling.
        sumzslave = 1.
c       sumzslave = zslave(0)
c       do i=1,nslaves-1
c         sumzslave = sumzslave + zslave(i)
c       enddo

c       --- All values of zslave must be > 0.
        do i=0,nslaves-1
c         if (zslave(i)/sumzslave*(zmmax-zmmin) <= 0.) then
          if (zslave(i) <= 0.) then
            call kaboom("domaindecomposeparticles: ERROR: The length of all
     & particle domains must be positive. Fix zslave appropriately.")
            return
          endif
        enddo

c       --- Set minimum z of each processor.
        zlast = zmmin
        do i=0,nslaves-1
          zpslmin(i) = zlast
          zpslmax(i) = zlast + zslave(i)/sumzslave*(zmmax - zmmin)
          zlast = zpslmax(i)
        enddo

c       --- When using solvergeom==AMRgeom, the particle decomposition must
c       --- be aligned with the grid.
        if (solvergeom==AMRgeom) then
         do i=0,nslaves-1
           zpslmin(i) = nint((zpslmin(i) - zmmin)/dz)*dz
           zpslmax(i) = nint((zpslmax(i) - zmmin)/dz)*dz
         enddo
        endif

c       --- This is only needed to avoid problems from round off in the
c       --- accumulation. From the loop above, zpslmax(nslaves-1) will
c       --- not be exactly the same as zmmax due to roundoff.
        if (zpslmax(nslaves-1) > zmmax) zpslmax(nslaves-1) = zmmax

c       --- Set iz and nz. This is done so that zmesh(izpslave) < zpslmin, and
c       --- zmesh(izpslave+nzpslave) > zpslmax.
        do i=0,nslaves-1
          izpslave(i) = int((zpslmin(i) - zmmin)/dz) - nzpguard
          nzpslave(i) = int((zpslmax(i) - zmmin)/dz) - izpslave(i) + 1 +
     &                  2*nzpguard
c         --- Make sure that the processors doesn't have grid cells
c         --- sticking out the end.
          if (izpslave(i) < 0) then
            nzpslave(i) = nzpslave(i) + izpslave(i)
            izpslave(i) = 0
          endif
          if (izpslave(i) + nzpslave(i) > nz) nzpslave(i) = nz - izpslave(i)
        enddo

      endif

      return
      end
c=============================================================================
c=============================================================================
