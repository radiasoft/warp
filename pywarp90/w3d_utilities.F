#include "top.h"
c=============================================================================
c@(#) File w3d_utilities.F, version $Revision: 1.2 $, $Date: 2004/07/29 17:26:30 $
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
      integer(ISZ):: n,nblocks
      real(kind=8):: indx(n),x(n),y(n),z(n),uz(n)
      real(kind=8):: xout(n),yout(n),zout(n),uzout(n)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)

c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = nint(indx(ip))
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = nint(indx(ip))
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        uzout(ii) = uz(ip)
        tcounts(ib) = tcounts(ib) + 1
      enddo

      return
      end
c=============================================================================
      subroutine sortparticlesbyindexgetisort(n,indx,x,y,z,nblocks,
     &                                        xout,yout,zout,isort,pcounts)
      integer(ISZ):: n,nblocks
      real(kind=8):: indx(n),x(n),y(n),z(n)
      real(kind=8):: xout(n),yout(n),zout(n)
      integer(ISZ):: isort(n)
      integer(ISZ):: pcounts(0:nblocks-1)

c Given particle data and an index array, the routine sorts the particles
c based on the index. It returns the sorted particles and a list of the
c number of particles for each index.

      integer(ISZ):: ip,ii,ib
      integer(ISZ):: tcounts(0:nblocks-1)

c     --- First, count how many in each block there are.
c     --- Note that the blocks numbers start at zero.
      pcounts = 0
      do ip=1,n
        ib = nint(indx(ip))
        pcounts(ib) = pcounts(ib) + 1
      enddo

c     --- Generate partial sum of pcounts
      tcounts(0) = 1
      do ib=0,nblocks-2
        tcounts(ib+1) = tcounts(ib) + pcounts(ib)
      enddo

c     --- Rearrange the particles
      do ip=1,n
        ib = nint(indx(ip))
        ii = tcounts(ib)
        xout(ii) = x(ip)
        yout(ii) = y(ip)
        zout(ii) = z(ip)
        isort(ii) = ip - 1
        tcounts(ib) = tcounts(ib) + 1
      enddo

      return
      end
c=============================================================================
      subroutine getgridngp3dpositiveonly(np,x,y,z,f,nx,ny,nz,grid,
     &                                    xmin,xmax,ymin,ymax,zmin,zmax,zgrid,
     &                                    l2symtry,l4symtry)
      integer(ISZ):: nx,ny,nz,np
      real(kind=8):: x(np), y(np), z(np), f(np)
      real(kind=8):: grid(0:nx,0:ny,0:nz)
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax,zgrid
      logical(ISZ):: l2symtry,l4symtry

c Gathers data from a 3-D mesh onto particle positions.
c Same as getgridngp3d, but only gathers the data if it is positive.

      integer(ISZ):: ip,ix,iy,iz
      real(kind=8):: gx,gy,gz
      real(kind=8):: dxi,dyi,dzi,wx,wy,wz
      real(kind=8):: xsymmetryplane,ysymmetryplane

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
        gx = (gx    - xmin)*dxi
        gy = (gy    - ymin)*dyi
        gz = (z(ip) - zmin - zgrid)*dzi

c       --- if within grid, gather
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gy .and. gy < ny .and.
     &      0. <= gz .and. gz < nz) then
          ix = int(gx)
          iy = int(gy)
          iz = int(gz)
          if (grid(ix,iy,iz) >= 0.) f(ip) = grid(ix,iy,iz)
        endif

      enddo

      return
      end
c=============================================================================
