#include "top.h"
c=============================================================================
c@(#) File w3d_utilities.F, version $Revision: 1.1 $, $Date: 2004/07/22 23:11:05 $
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
