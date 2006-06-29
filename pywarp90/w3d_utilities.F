#include "top.h"
c=============================================================================
c@(#) File w3d_utilities.F, version $Revision: 1.10 $, $Date: 2006/06/29 16:52:43 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package W3D of code WARP
c  3d electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
c=============================================================================
      subroutine sortparticlesbyindex(n,indx,x,y,z,uz,nblocks,
     &                                xout,yout,zout,uzout,pcounts)
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
      subroutine sortparticlesbyindexwithcounts(n,indx,x,y,z,uz,nblocks,
     &                                          xout,yout,zout,uzout,pcounts)
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

c     --- pcounts is passed in.

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
      if (lw3dtimesubs) timesortpart_byindexwithcounts =
     &                         timesortpart_byindexwithcounts +
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
      dyi = ny/(ymax - ymin)
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
      subroutine getichildandcount(np,x,y,z,ichild,nblocks,nperchild,
     &                             nx,ny,nz,grid,
     &                             xmin,xmax,ymin,ymax,zmin,zmax,zgrid,
     &                             l2symtry,l4symtry)
      use Subtimers3d
      integer(ISZ):: nx,ny,nz,np,nblocks
      real(kind=8):: x(np), y(np), z(np)
      integer(ISZ):: ichild(np)
      integer(ISZ):: nperchild(0:nblocks-1)
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
      dyi = ny/(ymax - ymin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

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
c       if (gx < xmin .or. gx > xmax .or.
c    &      gy < ymin .or. gy > ymax .or.
c    &      gz < zmin .or. gz > zmax) cycle

        gx = (gx - xmin)*dxi
        gy = (gy - ymin)*dyi
        gz = (gz - zmin)*dzi

        ix = int(gx)
        iy = int(gy)
        iz = int(gz)
        ichild(ip) = abs(grid(ix,iy,iz))
        nperchild(ichild(ip)) = nperchild(ichild(ip)) + 1

      enddo

!$OMP MASTER
      if (lw3dtimesubs) timegetichildandcount = timegetichildandcount +
     &                                          wtime() - substarttime
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
      dyi = ny/(ymax - ymin)
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
      subroutine putsortedefield(n,isort,tex,tey,tez,ex,ey,ez)
      use Subtimers3d
      integer(ISZ):: n
      integer(ISZ):: isort(0:n-1)
      real(kind=8),dimension(0:n-1):: tex,tey,tez,ex,ey,ez

      integer(ISZ):: i
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      do i=0,n-1
        ex(isort(i)) = tex(i)
        ey(isort(i)) = tey(i)
        ez(isort(i)) = tez(i)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timeputsortedefield = timeputsortedefield +
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
c --- Routines used by the mesh refined B field solver
c=============================================================================
      subroutine sortposandvelbyindex(n,indx,x,y,z,ux,uy,uz,gaminv,nblocks,
     &                                xout,yout,zout,uxout,uyout,uzout,gout,
     &                                pcounts)
      use Subtimers3d
      integer(ISZ):: n,nblocks
      integer(ISZ):: indx(n)
      real(kind=8):: x(n),y(n),z(n),ux(n),uy(n),uz(n),gaminv(n)
      real(kind=8):: xout(n),yout(n),zout(n),uxout(n),uyout(n),uzout(n),gout(n)
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
        gout(ii) = gaminv(ip)
        tcounts(ib) = tcounts(ib) + 1
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesortpart_byindex = timesortpart_byindex +
     &                                             wtime() - substarttime
!$OMP END MASTER

      return
      end
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
      print*,nx,dx,xmin,xmax,dxi,defaultcount

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
      print*,"X ",check
      if (.not. check) return
      check = arecoordinatesconsist(grid%ny,grid%dy,grid%ymin,grid%ymax,
     &                              grid%dyi)
      print*,"Y ",check
      if (.not. check) return
      check = arecoordinatesconsist(grid%nz,grid%dz,grid%zmin,grid%zmax,
     &                              grid%dzi)
      print*,"Z ",check
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
      integer(ISZ) function getnsndtsforsubcycling()
      use InGen3d,Only: solvergeom,XYZgeom
      use Fields3dParticles,Only: ndtsaveraging,nsndts

      integer(ISZ):: tmpnsndts

      if (solvergeom == XYZgeom) then
        if (ndtsaveraging == 0) then
          tmpnsndts = 1
        elseif (ndtsaveraging == 1) then
          tmpnsndts = nsndts
        elseif (ndtsaveraging == 2) then
          tmpnsndts = nsndts
        endif
      else
        tmpnsndts = 1
      endif

      getnsndtsforsubcycling = tmpnsndts
      return
      end
c=============================================================================
