#include "top.h"
c=======================================================================
c@(#) File F3D_CONDUCTORS.F, $Revision: 1.9 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  written by David P. Grote
c
c  It contains the routines for calculating the intersection of points
c  with various geometrical objects.
c=======================================================================
c=======================================================================
      subroutine setconductorparity(nn,ix,iy,iz,dels,parity,fuzz,fuzzsign,dfill)
      integer(ISZ):: nn
      integer(ISZ):: ix(nn),iy(nn),iz(nn),parity(nn)
      real(kind=8):: dels(0:5,nn)
      real(kind=8):: fuzz
      integer(ISZ):: fuzzsign
      real(kind=8):: dfill

c Set parity. For points inside, this is set to -1. For points near the surface,
c this is set to the parity of ix+iy+iz. Otherwise defaults to large integer.
c This assumes that the data has already been normalized with respect to the
c grid cell sizes.
c Format of dels array:
c  dels(0,:) -> delmx
c  dels(1,:) -> delpx
c  dels(2,:) -> delmy
c  dels(3,:) -> delpy
c  dels(4,:) -> delmz
c  dels(5,:) -> delpz

      integer(ISZ):: i,iparity
      integer(ISZ),parameter:: outside = 999
      
      do i=1,nn
        parity(i) = outside

        iparity = mod(ix(i)+iy(i)+iz(i),2)
        if (dels(0,i) < 1.+fuzzsign*fuzz) parity(i) = iparity
        if (dels(1,i) < 1.+fuzzsign*fuzz) parity(i) = iparity
        if (dels(2,i) < 1.+fuzzsign*fuzz) parity(i) = iparity
        if (dels(3,i) < 1.+fuzzsign*fuzz) parity(i) = iparity
        if (dels(4,i) < 1.+fuzzsign*fuzz) parity(i) = iparity
        if (dels(5,i) < 1.+fuzzsign*fuzz) parity(i) = iparity

c       --- If one direction is negative, all are negative. However,
c       --- because of the fuzz factor, each direction must be checked
c       --- since if one direction is between 0 and fuzz, the others
c       --- will all be positive.
        if (dels(0,i) < 0.+fuzz) parity(i) = -1
        if (dels(1,i) < 0.+fuzz) parity(i) = -1
        if (dels(2,i) < 0.+fuzz) parity(i) = -1
        if (dels(3,i) < 0.+fuzz) parity(i) = -1
        if (dels(4,i) < 0.+fuzz) parity(i) = -1
        if (dels(5,i) < 0.+fuzz) parity(i) = -1

c       --- To minimize the amount of conductor data, points deep inside the
c       --- conductor can be skipped. The depth is set by the input dfill.
c       --- The point must be at the depth in all directions to be skipped.
c       --- dfill is assumed to be passed in as a positive number.
        if (dels(0,i) <= -dfill .and. dels(1,i) <= -dfill .and.
     &      dels(2,i) <= -dfill .and. dels(3,i) <= -dfill .and.
     &      dels(4,i) <= -dfill .and. dels(5,i) <= -dfill) parity(i) = outside

      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine getconductorfacets(condid,nc,icnd,dels,gridn,griddd,gridmin)
      integer(ISZ):: condid
      integer(ISZ):: nc
      integer(ISZ):: icnd(0:2,nc)
      real(kind=8):: dels(0:5,nc)
      integer(ISZ):: gridn(0:2)
      real(kind=8):: griddd(0:2)
      real(kind=8):: gridmin(0:2)

      integer(ISZ):: iii(gridn(0)+1,gridn(1)+1,gridn(2)+1)
      integer(ISZ):: ic,ix,iy,iz,ox,oy,oz

      do ic=1,nc
        iii(icnd(0,ic),icnd(1,ic),icnd(2,ic)) = ic
      enddo

c     oos = [[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1]]
c     dds = [[+1,0,0],[-1,0,0],[+1,0,0],[-1,0,0],
c            [+1,0,0],[-1,0,0],[+1,0,0],[-1,0,0]]

c     do iz=0,gridn(2)-1
c       do iy=0,gridn(1)-1
c         do ix=0,gridn(0)-1

c         do oz=0,1
c         do oy=0,1
c         do ox=0,1

c       for oo,dd in map(None,oos,dds):
c         if iii[ix+oo[0],iy+oo[1],iz+oo[2]] > 0:
c           pp = newfacet(ix,iy,iz,oo,dd,+1,iii,dels,gridmin,griddd)
c           if len(pp) == 0:
c             pp = newfacet(ix,iy,iz,oo,dd,0,iii,dels,gridmin,griddd)
c           if len(pp) > 0:
c             model.FacetedPolygon(pp,color=color)
c       iii[ix:ix+2,iy:iy+2,iz:iz+2] = abs(iii[ix:ix+2,iy:iy+2,iz:iz+2])
c model.Display()
c return model.scene

cef nextdir(position,direction,parity):
c parity = (sum(array(position)) + parity + 1) % 2
c newdirection = [0,0,0]
c if direction[0] != 0:   newdirection[1+parity] = 1
c elif direction[1] != 0: newdirection[2-2*parity] = 1
c elif direction[2] != 0: newdirection[parity] = 1
c if position[0] == 1: newdirection[0] = -newdirection[0]
c if position[1] == 1: newdirection[1] = -newdirection[1]
c if position[2] == 1: newdirection[2] = -newdirection[2]
c return newdirection

c --- This would be a little less painful is lists were hashable.
cef newfacet(ix,iy,iz,oo,dd,parity,iii,dels,gridmin,griddd,
c dddict = {211:[+1,0,0],121:[0,+1,0],112:[0,0,+1],
c            11:[-1,0,0],101:[0,-1,0],110:[0,0,-1]},
c dirdict = {  0:[{211:0,112:1,121:2},{0:211,1:112,2:121}],
c            100:[{112:0, 11:1,121:2},{0:112,1: 11,2:121}],
c             10:[{211:0,101:1,112:2},{0:211,1:101,2:112}],
c            110:[{ 11:0,112:1,101:2},{0: 11,1:112,2:101}],
c              1:[{211:0,121:1,110:2},{0:211,1:121,2:110}],
c            101:[{110:0,121:1, 11:2},{0:110,1:121,2: 11}],
c             11:[{211:0,110:1,101:2},{0:211,1:110,2:101}],
c            111:[{ 11:0,101:1,110:2},{0: 11,1:101,2:110}]},
c delsdict = { 11:0,211:1,101:2,121:3,110:4,112:5},
c too = lambda oo : 100*oo[0] + 10*oo[1] + oo[2],
c tdd = lambda dd : 100*(1+dd[0]) + 10*(1+dd[1]) + (1+dd[2])):
c """
cork routine used by visualizeconductors.
c """
c pp = []
c done = 0
c currentoo = oo
c currentdd = dd
c ncorners = 1
c while not done:
c   ic = abs(iii[ix+currentoo[0],iy+currentoo[1],iz+currentoo[2]])
c   iii[ix+currentoo[0],iy+currentoo[1],iz+currentoo[2]] = -ic
c   if ic != 0 and dels[delsdict[tdd(currentdd)],ic-1] <= 1.000001:
c     newpoint = (array([ix+currentoo[0],iy+currentoo[1],iz+currentoo[2]]) +
c                       dels[delsdict[tdd(currentdd)],ic-1]*array(currentdd))
c     pp.append(newpoint*griddd+gridmin)
c     #dirs = dirdict[too(currentoo)]
c     #currentddnum = dirs[0][tdd(currentdd)]
c     #nextddnum = (currentddnum + 2*parity-1) % 3
c     #nextdd = dddict[dirs[1][nextddnum]]
c     nextdd = nextdir(currentoo,currentdd,parity)
c     nextoo = currentoo
c   else:
c     nextoo = list(array(currentoo) + array(currentdd))
c     #nextdirs = dirdict[too(nextoo)]
c     #nextddnum = (nextdirs[0][tdd(list(-array(currentdd)))] + 2*parity-1) % 3
c     #nextdd = dddict[nextdirs[1][nextddnum]]
c     nextdd = nextdir(nextoo,-array(currentdd),parity)
c     ncorners = ncorners + 1
c   currentoo = nextoo
c   currentdd = nextdd
c   if ncorners > 8: done = 1
c   if len(pp) > 6:
c     pp = []
c     done = 1
c   if currentoo == oo and currentdd == dd: done = 1
c return pp





      return
      end
c=======================================================================
c=======================================================================
      subroutine ZPlaneConductorF(zcent,zsign,
     &                            n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                            fuzz)
      real(kind=8):: zcent,zsign
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a z plane
c The notation is that a negative distance means that the point is inside.
c A positive zsign means that the conductor is the region above zcent.

      integer(ISZ):: i

      do i=1,n

        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

        if (zsign > 0.) then
          if (z(i) <= zcent) delpz(i) = zcent - z(i)
          if (z(i) >= zcent) delmz(i) = zcent - z(i)
        else if (zsign < 0.) then
          if (z(i) <= zcent) delpz(i) = z(i) - zcent
          if (z(i) >= zcent) delmz(i) = z(i) - zcent
        endif

      enddo

      return
      end
c=======================================================================
      subroutine BoxConductorF(xsize,ysize,zsize,xcent,ycent,zcent,
     &                            n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                            fuzz)
      real(kind=8):: xsize,ysize,zsize,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a z plane
c The notation is that a negative distance means that the point is inside.
c A positive zsign means that the conductor is the region above zcent.

      integer(ISZ):: i
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax

      xmin = xcent - xsize/2.
      xmax = xcent + xsize/2.
      ymin = ycent - ysize/2.
      ymax = ycent + ysize/2.
      zmin = zcent - zsize/2.
      zmax = zcent + zsize/2.

      do i=1,n

        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

        if (xmin - fuzz <= x(i) .and. x(i) <= xmax + fuzz .and.
     &      ymin - fuzz <= y(i) .and. y(i) <= ymax + fuzz) then

          if (z(i) < zmin) then
            delpz(i) = zmin - z(i)
          else if (zmin <= z(i) .and. z(i) <= zmax) then
            delmz(i) = zmin - z(i)
            delpz(i) = z(i) - zmax
          else if (zmax < z(i)) then
            delmz(i) = z(i) - zmax
          endif

        endif

        if (ymin - fuzz <= y(i) .and. y(i) <= ymax + fuzz .and.
     &      zmin - fuzz <= z(i) .and. z(i) <= zmax + fuzz) then

          if (x(i) < xmin) then
            delpx(i) = xmin - x(i)
          else if (xmin <= x(i) .and. x(i) <= xmax) then
            delmx(i) = xmin - x(i)
            delpx(i) = x(i) - xmax
          else if (xmax < x(i)) then
            delmx(i) = x(i) - xmax
          endif

        endif

        if (zmin - fuzz <= z(i) .and. z(i) <= zmax + fuzz .and.
     &      xmin - fuzz <= x(i) .and. x(i) <= xmax + fuzz) then

          if (y(i) < ymin) then
            delpy(i) = ymin - y(i)
          else if (ymin <= y(i) .and. y(i) <= ymax) then
            delmy(i) = ymin - y(i)
            delpy(i) = y(i) - ymax
          else if (ymax < y(i)) then
            delmy(i) = y(i) - ymax
          endif

        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZCylinderConductorF(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder aligned with
c the z-axis.
c The notation is that a negative distance means that the point is inside.

      real(kind=8):: zmin,zmax
      integer(ISZ):: i
      real(kind=8):: xx,yy,rsq

      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = (xx**2 + yy**2)

c       --- If point is not within the z extent of the cylinder
        if (z(i) < zmin - fuzz) then
          if (rsq <= rad**2) delpz(i) = zmin - z(i)
          cycle
        endif
        if (z(i) > zmax + fuzz) then
          if (rsq <= rad**2) delmz(i) = z(i) - zmax
          cycle
        endif

c       --- Point inside the cylinder
        if (rsq < rad**2) then
          delmx(i) = -(sqrt(rad**2 - yy**2) + xx)
          delpx(i) = -(sqrt(rad**2 - yy**2) - xx)
          delmy(i) = -(sqrt(rad**2 - xx**2) + yy)
          delpy(i) = -(sqrt(rad**2 - xx**2) - yy)
          delmz(i) = zmin - z(i)
          delpz(i) = z(i) - zmax
          cycle
        endif

c       --- Point is outside the cylinder
        if (rsq > rad**2) then
          if (abs(yy) <= rad) then
            if (xx > 0) delmx(i) = ( xx - sqrt(max(0.,rad**2 - yy**2)))
            if (xx < 0) delpx(i) = (-xx - sqrt(max(0.,rad**2 - yy**2)))
          endif
          if (abs(xx) <= rad) then
            if (yy > 0) delmy(i) = ( yy - sqrt(max(0.,rad**2 - xx**2)))
            if (yy < 0) delpy(i) = (-yy - sqrt(max(0.,rad**2 - xx**2)))
          endif
          cycle
        endif

c       --- The point is on the cylinder exactly.
        delmx(i) = 0.
        delpx(i) = 0.
        delmy(i) = 0.
        delpy(i) = 0.
        delmz(i) = 0.
        delpz(i) = 0.

      enddo

      return
      end
c=======================================================================
      subroutine ZRoundedCylinderConductorF(rad,length,rad2,
     &                              xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder aligned with
c the z-axis. The corners of the cylinder are rounded with radius rad2.
c The notation is that a negative distance means that the point is inside.

      real(kind=8):: zmin,zmax
      integer(ISZ):: i
      real(kind=8):: xx,yy,rsq,rr,r,zround

      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = (xx**2 + yy**2)
        rr = sqrt(rsq)
        if (rad-rad2 <= rr .and. rr <= rad) then
          zround = sqrt(max(0.,rad2**2 - (rr-rad+rad2)**2))
        endif

c       --- If point is not within the z extent of the cylinder
        if (z(i) < zmin - fuzz) then
          if (rsq <= (rad-rad2)**2) then
            delpz(i) = zmin - z(i)
          else if (rsq <= rad**2) then
            delpz(i) = zmin + rad2 - z(i) - zround
          endif
          cycle
        endif
        if (z(i) > zmax + fuzz) then
          if (rsq <= (rad-rad2)**2) then
            delmz(i) = z(i) - zmax
          else if (rsq <= rad**2) then
            delmz(i) = z(i) - zmax + rad2 - zround
          endif
          cycle
        endif

        if (rsq <= rad**2) then
          if (rsq > (rad-rad2)**2) then
            if (z(i) <= zmin + rad2 - zround) then
              delpz(i) = zmin + rad2 - z(i) - zround
            else if (z(i) < zmax - rad2 + zround) then
              delpz(i) = z(i) - zmax + rad2 - zround
              delmz(i) = zmin + rad2 - zround - z(i)
            else
              delmz(i) = z(i) - zmax + rad2 - zround
            endif
          else
            delmz(i) = zmin - z(i)
            delpz(i) = z(i) - zmax
          endif
        endif

c       --- Get radius of conductor.  Care is needed in the if's to avoid
c       --- sqrt's of negative numbers.
        if (z(i) <= zmin) then
          r = rad - rad2
        else if (z(i) < zmin + rad2) then
          r = rad - rad2 + sqrt(rad2**2 - (zmin+rad2-z(i))**2)
        else if (z(i) <= zmax - rad2) then
          r = rad
        else if (z(i) < zmax) then
          r = rad - rad2 + sqrt(rad2**2 - (z(i)-zmax+rad2)**2)
        else
          r = rad - rad2
        endif

c       --- Point inside the cylinder
        if (rsq < r**2) then
          delmx(i) = -(sqrt(r**2 - yy**2) + xx)
          delpx(i) = -(sqrt(r**2 - yy**2) - xx)
          delmy(i) = -(sqrt(r**2 - xx**2) + yy)
          delpy(i) = -(sqrt(r**2 - xx**2) - yy)
          cycle
        endif

c       --- Point is outside the cylinder
        if (rsq > r**2) then
          if (abs(yy) <= r) then
            if (xx > 0) delmx(i) = ( xx - sqrt(max(0.,r**2 - yy**2)))
            if (xx < 0) delpx(i) = (-xx - sqrt(max(0.,r**2 - yy**2)))
          endif
          if (abs(xx) <= r) then
            if (yy > 0) delmy(i) = ( yy - sqrt(max(0.,r**2 - xx**2)))
            if (yy < 0) delpy(i) = (-yy - sqrt(max(0.,r**2 - xx**2)))
          endif
          cycle
        endif

c       --- The point is on the cylinder exactly.
        delmx(i) = 0.
        delpx(i) = 0.
        delmy(i) = 0.
        delpy(i) = 0.
        delmz(i) = 0.
        delpz(i) = 0.

      enddo

      return
      end
c=======================================================================
      subroutine ZCylinderOutConductorF(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder aligned with
c the z-axis.
c The notation is that a negative distance means that the point is inside.

      real(kind=8):: zmin,zmax
      integer(ISZ):: i
      real(kind=8):: xx,yy,rsq

      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = (xx**2 + yy**2)

c       --- If point is not within the z extent of the cylinder
        if (z(i) < zmin - fuzz) then
          if (rsq >= rad**2) delpz(i) = zmin - z(i)
          cycle
        endif
        if (z(i) > zmax + fuzz) then
          if (rsq >= rad**2) delmz(i) = z(i) - zmax
          cycle
        endif

c       --- Point inside the cylinder
        if (rsq < rad**2) then
          delmx(i) = (sqrt(rad**2 - yy**2) + xx)
          delpx(i) = (sqrt(rad**2 - yy**2) - xx)
          delmy(i) = (sqrt(rad**2 - xx**2) + yy)
          delpy(i) = (sqrt(rad**2 - xx**2) - yy)
          cycle
        endif

c       --- Point is outside the cylinder
        if (rsq > rad**2) then
          delmx(i) = -LARGEPOS
          delpx(i) = -LARGEPOS
          delmy(i) = -LARGEPOS
          delpy(i) = -LARGEPOS
          if (abs(yy) <= rad) then
            if (xx > 0) delmx(i) = -( xx - sqrt(max(0.,rad**2 - yy**2)))
            if (xx < 0) delpx(i) = -(-xx - sqrt(max(0.,rad**2 - yy**2)))
          endif
          if (abs(xx) <= rad) then
            if (yy > 0) delmy(i) = -( yy - sqrt(max(0.,rad**2 - xx**2)))
            if (yy < 0) delpy(i) = -(-yy - sqrt(max(0.,rad**2 - xx**2)))
          endif
          delmz(i) = zmin - z(i)
          delpz(i) = z(i) - zmax
          cycle
        endif

c       --- The point is on the cylinder exactly.
        delmx(i) = 0.
        delpx(i) = 0.
        delmy(i) = 0.
        delpy(i) = 0.
        delmz(i) = 0.
        delpz(i) = 0.

      enddo

      return
      end
c=======================================================================
      subroutine ZRoundedCylinderOutConductorF(rad,length,rad2,
     &                              xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder aligned with
c the z-axis. The corners of the cylinder are rounded with radius rad2.
c The notation is that a negative distance means that the point is inside.

      real(kind=8):: zmin,zmax
      integer(ISZ):: i
      real(kind=8):: xx,yy,rsq,rr,r,zround

      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = (xx**2 + yy**2)
        rr = sqrt(rsq)
        if (rad <= rr .and. rr <= rad+rad2) then
          zround = sqrt(max(0.,rad2**2 - (rad+rad2-rr)**2))
        endif

c       --- If point is not within the z extent of the cylinder
        if (z(i) < zmin - fuzz) then
          if (rsq >= (rad+rad2)**2) then
            delpz(i) = zmin - z(i)
          else if (rsq >= rad**2) then
            delpz(i) = zmin + rad2 - z(i) - zround
          endif
          cycle
        endif
        if (z(i) > zmax + fuzz) then
          if (rsq >= (rad+rad2)**2) then
            delmz(i) = z(i) - zmax
          else if (rsq >= rad**2) then
            delmz(i) = z(i) - zmax + rad2 - zround
          endif
          cycle
        endif

        if (rsq >= rad**2) then
          if (rsq < (rad+rad2)**2) then
            if (z(i) <= zmin + rad2 - zround) then
              delpz(i) = zmin + rad2 - z(i) - zround
            else if (z(i) < zmax - rad2 + zround) then
              delpz(i) = z(i) - zmax + rad2 - zround
              delmz(i) = zmin + rad2 - zround - z(i)
            else
              delmz(i) = z(i) - zmax + rad2 - zround
            endif
          else
            delmz(i) = zmin - z(i)
            delpz(i) = z(i) - zmax
          endif
        endif

c       --- Get radius of conductor.  Care is needed in the if's to avoid
c       --- sqrt's of negative numbers.
        if (z(i) <= zmin) then
          r = rad + rad2
        else if (z(i) < zmin + rad2) then
          r = rad + rad2 - sqrt(rad2**2 - (zmin+rad2-z(i))**2)
        else if (z(i) <= zmax - rad2) then
          r = rad
        else if (z(i) < zmax) then
          r = rad + rad2 - sqrt(rad2**2 - (z(i)-zmax+rad2)**2)
        else
          r = rad + rad2
        endif

c       --- Point inside the cylinder
        if (rsq < r**2) then
          delmx(i) = (sqrt(r**2 - yy**2) + xx)
          delpx(i) = (sqrt(r**2 - yy**2) - xx)
          delmy(i) = (sqrt(r**2 - xx**2) + yy)
          delpy(i) = (sqrt(r**2 - xx**2) - yy)
          cycle
        endif

c       --- Point is outside the cylinder
        if (rsq > r**2) then
          delmx(i) = -LARGEPOS
          delpx(i) = -LARGEPOS
          delmy(i) = -LARGEPOS
          delpy(i) = -LARGEPOS
          if (abs(yy) <= r) then
            if (xx > 0) delmx(i) = -( xx - sqrt(max(0.,r**2 - yy**2)))
            if (xx < 0) delpx(i) = -(-xx - sqrt(max(0.,r**2 - yy**2)))
          endif
          if (abs(xx) <= r) then
            if (yy > 0) delmy(i) = -( yy - sqrt(max(0.,r**2 - xx**2)))
            if (yy < 0) delpy(i) = -(-yy - sqrt(max(0.,r**2 - xx**2)))
          endif
          cycle
        endif

c       --- The point is on the cylinder exactly.
        delmx(i) = 0.
        delpx(i) = 0.
        delmy(i) = 0.
        delpy(i) = 0.
        delmz(i) = 0.
        delpz(i) = 0.

      enddo

      return
      end
c=======================================================================
      subroutine XCylinderConductorF(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder aligned with
c the x-axis.
c The notation is that a negative distance means that the point is inside.

      call ZCylinderConductorF(rad,length,ycent,zcent,xcent,
     &                         n,y,z,x,delmy,delpy,delmz,delpz,delmx,delpx,fuzz)

      return
      end
c=======================================================================
      subroutine YCylinderConductorF(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder aligned with
c the y-axis.
c The notation is that a negative distance means that the point is inside.

      call ZCylinderConductorF(rad,length,zcent,xcent,ycent,
     &                         n,z,x,y,delmz,delpz,delmx,delpx,delmy,delpy,fuzz)

      return
      end
c=======================================================================
      subroutine SphereConductorF(rad,xcent,ycent,zcent,
     &                            n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                            fuzz)
      real(kind=8):: rad,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a sphere.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rsq

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rsq = xx**2 + yy**2 + zz**2

c       --- Point inside the sphere
        if (rsq < rad**2) then
          delmx(i) = -(sqrt(rad**2 - yy**2 - zz**2) + xx)
          delpx(i) = -(sqrt(rad**2 - yy**2 - zz**2) - xx)
          delmy(i) = -(sqrt(rad**2 - xx**2 - zz**2) + yy)
          delpy(i) = -(sqrt(rad**2 - xx**2 - zz**2) - yy)
          delmz(i) = -(sqrt(rad**2 - xx**2 - yy**2) + zz)
          delpz(i) = -(sqrt(rad**2 - xx**2 - yy**2) - zz)
        endif

c       --- Point outside the sphere
        if (rsq > rad**2) then
          if ((yy**2 + zz**2) <= rad**2) then
            if (xx > 0) delmx(i) = ( xx - sqrt(max(0.,rad**2 - yy**2 - zz**2)))
            if (xx < 0) delpx(i) = (-xx - sqrt(max(0.,rad**2 - yy**2 - zz**2)))
          endif
          if ((xx**2 + zz**2) <= rad**2) then
            if (yy > 0) delmy(i) = ( yy - sqrt(max(0.,rad**2 - xx**2 - zz**2)))
            if (yy < 0) delpy(i) = (-yy - sqrt(max(0.,rad**2 - xx**2 - zz**2)))
          endif
          if ((xx**2 + yy**2) <= rad**2) then
            if (zz > 0) delmz(i) = ( zz - sqrt(max(0.,rad**2 - xx**2 - yy**2)))
            if (zz < 0) delpz(i) = (-zz - sqrt(max(0.,rad**2 - xx**2 - yy**2)))
          endif
        endif

c       --- The point is on the sphere exactly.
        if (rsq == rad**2) then
          delmx(i) = 0.
          delpx(i) = 0.
          delmy(i) = 0.
          delpy(i) = 0.
          delmz(i) = 0.
          delpz(i) = 0.
        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZConeConductorF(r_zmin,r_zmax,length,xcent,ycent,zcent,
     &                           n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                           fuzz)
      real(kind=8):: r_zmin,r_zmax,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cone aligned with
c the z-axis.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: zmin,zmax
      real(kind=8):: xx,yy,rsq,r,rcone,zcone

c     --- First, find min and max extent of the cone.
      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Find radius of point relative to cone axis
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = xx**2 + yy**2
        r = sqrt(rsq)

c       --- Find z location where transverse line through point intersects
c       --- the cone.
        zcone = zmin + (zmax - zmin)*(r - r_zmin)/(r_zmax - r_zmin)
        if (r == r_zmax) zcone = zmax

c       --- Calculate z deltas
        if (r <= min(r_zmin,r_zmax)) then
c         --- point is inside of smaller radius
          if (z(i) >= zmin) delmz(i) = zmin - z(i)
          if (z(i) >= zmax) delmz(i) = z(i) - zmax
          if (z(i) <= zmax) delpz(i) = z(i) - zmax
          if (z(i) <= zmin) delpz(i) = zmin - z(i)
        else if (r <= max(r_zmin,r_zmax)) then
c         --- point is inside of larger radius
          if (r_zmin < r_zmax) then
            if (z(i) < zcone) delpz(i) = zcone - z(i)
            if (zcone <= z(i) .and. z(i) <= zmax) then
              delmz(i) = zcone - z(i)
              delpz(i) = z(i) - zmax
            endif
            if (z(i) > zmax) delmz(i) = z(i) - zmax
          endif
          if (r_zmin > r_zmax) then
            if (z(i) > zcone) delmz(i) = z(i) - zcone
            if (zcone >= z(i) .and. z(i) >= zmin) then
              delpz(i) = z(i) - zcone
              delmz(i) = zmin - z(i)
            endif
            if (z(i) < zmin) delpz(i) = zmin - z(i)
          endif
        endif

c       --- If point outside the z-extent of the cone, then skip the next
c       --- part.
        if (z(i) < zmin - fuzz .or. zmax + fuzz < z(i)) cycle

c       --- Calculate the radius of the cone at the z-location of the point
        rcone = r_zmin + (r_zmax - r_zmin)*(z(i) - zmin)/(zmax - zmin)
c       --- If z == zmax, the set rcone equal to r_zmax exactly since the
c       --- expression above will have roundoff problems.
        if (z(i) == zmax) rcone = r_zmax

c       --- Point inside the cone
        if (r < rcone) then
          delmx(i) = -(sqrt(rcone**2 - yy**2) + xx)
          delpx(i) = -(sqrt(rcone**2 - yy**2) - xx)
          delmy(i) = -(sqrt(rcone**2 - xx**2) + yy)
          delpy(i) = -(sqrt(rcone**2 - xx**2) - yy)

c       --- Point is outside the cone
        else if (r > rcone) then
          if (abs(yy) <= rcone) then
            if (xx > 0) delmx(i) = ( xx - sqrt(max(0.,rcone**2 - yy**2)))
            if (xx < 0) delpx(i) = (-xx - sqrt(max(0.,rcone**2 - yy**2)))
          endif
          if (abs(xx) <= rcone) then
            if (yy > 0) delmy(i) = ( yy - sqrt(max(0.,rcone**2 - xx**2)))
            if (yy < 0) delpy(i) = (-yy - sqrt(max(0.,rcone**2 - xx**2)))
          endif

        else
c         --- The point is on the cone exactly.
          delmx(i) = 0.
          delpx(i) = 0.
          delmy(i) = 0.
          delpy(i) = 0.
          delmz(i) = 0.
          delpz(i) = 0.

        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZConeOutConductorF(r_zmin,r_zmax,length,xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: r_zmin,r_zmax,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cone aligned with
c the z-axis.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: zmin,zmax
      real(kind=8):: xx,yy,rsq,r,rcone,zcone

c     --- First, find min and max extent of the cone.
      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Find radius of point relative to cone axis
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = xx**2 + yy**2
        r = sqrt(rsq)

c       --- Find z location where transverse line through point intersects
c       --- the cone.
        zcone = zmin + (zmax - zmin)*(r - r_zmin)/(r_zmax - r_zmin)
        if (r == r_zmax) zcone = zmax

c       --- Calculate z deltas
        if (r >= max(r_zmin,r_zmax)) then
c         --- point is outside of larger radius
          if (z(i) >= zmin) delmz(i) = zmin - z(i)
          if (z(i) >= zmax) delmz(i) = z(i) - zmax
          if (z(i) <= zmax) delpz(i) = z(i) - zmax
          if (z(i) <= zmin) delpz(i) = zmin - z(i)
        else if (r >= min(r_zmin,r_zmax)) then
c         --- point is outside of smaller radius
          if (r_zmin > r_zmax) then
            if (z(i) < zcone) delpz(i) = zcone - z(i)
            if (zcone <= z(i) .and. z(i) <= zmax) then
              delmz(i) = zcone - z(i)
              delpz(i) = z(i) - zmax
            endif
            if (z(i) > zmax) delmz(i) = z(i) - zmax
          endif
          if (r_zmin < r_zmax) then
            if (z(i) > zcone) delmz(i) = z(i) - zcone
            if (zcone >= z(i) .and. z(i) >= zmin) then
              delpz(i) = z(i) - zcone
              delmz(i) = zmin - z(i)
            endif
            if (z(i) < zmin) delpz(i) = zmin - z(i)
          endif
        endif

c       --- If point outside the z-extent of the cone, then skip the next
c       --- part.
        if (z(i) < zmin - fuzz .or. zmax + fuzz < z(i)) cycle

c       --- Calculate the radius of the cone at the z-location of the point
        rcone = r_zmin + (r_zmax - r_zmin)*(z(i) - zmin)/(zmax - zmin)
c       --- If z == zmax, the set rcone equal to r_zmax exactly since the
c       --- expression above will have roundoff problems.
        if (z(i) == zmax) rcone = r_zmax

c       --- Point inside the cone
        if (r < rcone) then
          delmx(i) = (sqrt(rcone**2 - yy**2) + xx)
          delpx(i) = (sqrt(rcone**2 - yy**2) - xx)
          delmy(i) = (sqrt(rcone**2 - xx**2) + yy)
          delpy(i) = (sqrt(rcone**2 - xx**2) - yy)

c       --- Point is outside the cone
        else if (r > rcone) then
          delmx(i) = -LARGEPOS
          delpx(i) = -LARGEPOS
          delmy(i) = -LARGEPOS
          delpy(i) = -LARGEPOS
          if (abs(yy) <= rcone) then
            if (xx > 0) delmx(i) = -( xx - sqrt(max(0.,rcone**2 - yy**2)))
            if (xx < 0) delpx(i) = -(-xx - sqrt(max(0.,rcone**2 - yy**2)))
          endif
          if (abs(xx) <= rcone) then
            if (yy > 0) delmy(i) = -( yy - sqrt(max(0.,rcone**2 - xx**2)))
            if (yy < 0) delpy(i) = -(-yy - sqrt(max(0.,rcone**2 - xx**2)))
          endif

        else
c         --- The point is on the cone exactly.
          delmx(i) = 0.
          delpx(i) = 0.
          delmy(i) = 0.
          delpy(i) = 0.
          delmz(i) = 0.
          delpz(i) = 0.

        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZTorusConductorF(r1,r2,xcent,ycent,zcent,
     &                             n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                            fuzz)
      real(kind=8):: r1,r2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a torus aligned with
c the z-axis.
c r1 is radius relative to z-axis
c r2 is radius relative to theta-axis
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: zl,zr,ri,ro
      real(kind=8):: xx,yy,rsq,r

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Find radius of point relative to torus axis
        xx = x(i) - xcent
        yy = y(i) - ycent
        rsq = xx**2 + yy**2
        r = sqrt(rsq)

c       --- Only need to set delmz and delpz for points with r within r1+-r2
        if (r1 - r2 <= r .and. r <= r1 + r2) then
          zl = zcent - sqrt(max(0.,r2**2 - (r - r1)**2))
          zr = zcent + sqrt(max(0.,r2**2 - (r - r1)**2))
          if (r == r1 - r2 .or. r == r1 + r2) then
            zl = zcent
            zr = zcent
          endif
          if (z(i) < zl) delpz(i) = zl - z(i)
          if (zl <= z(i) .and. z(i) <= zr) then
            delmz(i) = zl - z(i)
            delpz(i) = z(i) - zr
          endif
          if (z(i) > zr) delmz(i) = z(i) - zr
        endif

c       --- Only need to set transverse distances for points with
c       --- z within zcent+-r2.
        if (zcent-r2 <= z(i) .and. z(i) <= zcent+r2) then
          ri = r1 - sqrt(max(0.,r2**2 - (z(i) - zcent)**2))
          ro = r1 + sqrt(max(0.,r2**2 - (z(i) - zcent)**2))
          if (z(i) == zcent) then
            ri = r1 - r2
            ro = r1 + r2
          else if (z(i) == zcent - r2 .or. z(i) == zcent + r2) then
            ri = r1
            ro = r1
          endif
          if (r < ri) then
            delmx(i) = sqrt(ri**2 - yy**2) + xx
            delpx(i) = sqrt(ri**2 - yy**2) - xx
            delmy(i) = sqrt(ri**2 - xx**2) + yy
            delpy(i) = sqrt(ri**2 - xx**2) - yy
          else if (ri == r .or. r == ro) then
            delmx(i) = 0.
            delpx(i) = 0.
            delmy(i) = 0.
            delpy(i) = 0.
          else if (ri < r .and. r < ro) then
            delmx(i) = -(sqrt(ro**2 - yy**2) + xx)
            delpx(i) = -(sqrt(ro**2 - yy**2) - xx)
            delmy(i) = -(sqrt(ro**2 - xx**2) + yy)
            delpy(i) = -(sqrt(ro**2 - xx**2) - yy)
            if (abs(yy) < ri) then
              if (xx <= 0.) delpx(i) = xx + sqrt(max(0.,ri**2 - yy**2))
              if (xx >= 0.) delmx(i) = sqrt(max(0.,ri**2 - yy**2)) - xx
            else if (abs(yy) == ri) then
              if (xx <= 0.) delpx(i) = xx
              if (xx >= 0.) delmx(i) = -xx
            endif
            if (abs(xx) < ri) then
              if (yy <= 0.) delpy(i) = yy + sqrt(max(0.,ri**2 - xx**2))
              if (yy >= 0.) delmy(i) = sqrt(max(0.,ri**2 - xx**2)) - yy
            else if (abs(xx) == ri) then
              if (yy <= 0.) delpy(i) = yy
              if (yy >= 0.) delmy(i) = -yy
            endif
          else if (r > ro) then
            if (abs(yy) < ro) then
              if (xx <= 0.) delpx(i) = -sqrt(max(0.,ro**2 - yy**2)) - xx
              if (xx >= 0.) delmx(i) = xx - sqrt(max(0.,ro**2 - yy**2))
            endif
            if (abs(xx) < ro) then
              if (yy <= 0.) delpy(i) = -sqrt(max(0.,ro**2 - xx**2)) - yy
              if (yy >= 0.) delmy(i) = yy - sqrt(max(0.,ro**2 - xx**2))
            endif
          endif

        endif
      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine CylinderConductorF(rad,length,theta,phi,xcent,ycent,zcent,
     &                              n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      real(kind=8):: rad,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder at an arbitrary
c angle.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: ctheta,stheta,cphi,sphi,ss,cci,sci,si,lenhalf
      real(kind=8):: xx,yy,zz
      real(kind=8):: a,b,c
      real(kind=8):: z1,z2,zl,zu,zlend,zuend
      real(kind=8):: x1,x2,xl,xu,xlend,xuend
      real(kind=8):: y1,y2,yl,yu,ylend,yuend

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)
      cci = 1./(ctheta*cphi)
      sci = 1./(stheta*cphi)
      si  = 1./sphi
      lenhalf = length/2.

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
c       zp = zz*cos(theta)*cos(phi) + xx*sin(theta)*cos(phi) + yy*sin(phi)
c       xp = xx*cos(theta) - zz*sin(theta)
c       yp = yy*cos(phi) - zz*cos(theta)*sin(phi) - xx*sin(theta)*sin(phi)
c       rsq = (Ap**2 + yp**2)

c       --- delta z
c       --- Coefficients of the quadratic
        a = stheta**2 + (ctheta*sphi)**2
        b = -2*xx*ctheta*stheta - 2*(yy*cphi - xx*stheta*sphi)*ctheta*sphi
        c = (xx*ctheta)**2 + (yy*cphi - xx*stheta*sphi)**2 - rad**2
c       --- Check if z-line intersects the cylinder
        if (b**2 - 4*a*c >= 0) then
c         --- Get the two points of intersection
          if (a == 0) then
            z1 = -c/b
            z2 = z1
          else
            ss = sqrt(max(0.,b**2 - 4*a*c))
            z1 = (-b-ss)/(2*a)
            z2 = (-b+ss)/(2*a)
          endif
          zl = min(z1,z2)
          zu = max(z1,z2)
c         --- Get the point where the z-line intersects the end planes
          z1 = (-lenhalf - xx*stheta*cphi - yy*sphi)*cci
          z2 = (+lenhalf - xx*stheta*cphi - yy*sphi)*cci
          zlend = min(z1,z2)
          zuend = max(z1,z2)
c         --- Now, get distances to surface
          if (zz < zl) then
            if (zl <= zlend .and. zlend <= zu) delpz(i) = zlend - zz
            if (zlend <= zl .and. zl <= zuend) delpz(i) = zl - zz
          else if (zl <= zz .and. zz <= zu) then
            if (zl <= zlend .and. zlend <= zu .and. zz <= zuend) then
              if (zz <= zlend) then
                delpz(i) = zlend - zz
              else
                delpz(i) = zz - min(zu,zuend)
                delmz(i) = zlend - zz
              endif
            else if (zlend <= zl .and. zu <= zuend) then
              delmz(i) = zl - zz
              delpz(i) = zz - zu
            else if (zl <= zuend .and. zuend <= zu) then
              if (zz <= zuend) then
                delpz(i) = zz - zuend
                delmz(i) = max(zl,zlend) - zz
              else
                delmz(i) = zz - zuend
              endif
            endif
          else if (zu < zz) then
            if (zu >= zuend .and. zuend >= zl) delmz(i) = zz - zuend
            if (zuend >= zu .and. zu >= zlend) delmz(i) = zz - zu
          endif
        endif

c       --- delta x
        a = (ctheta)**2 + (stheta*sphi)**2
        b = -2*zz*stheta*ctheta - 2*(yy*cphi - zz*ctheta*sphi)*stheta*sphi
        c = (zz*stheta)**2 + (yy*cphi - zz*ctheta*sphi)**2 - rad**2
        if (b**2 - 4*a*c >= 0) then
          if (a == 0) then
            x1 = -c/b
            x2 = x1
          else
            ss = sqrt(max(0.,b**2 - 4*a*c))
            x1 = (-b-ss)/(2*a)
            x2 = (-b+ss)/(2*a)
          endif
          xl = min(x1,x2)
          xu = max(x1,x2)
c         --- Get the point where the x-line intersects the end planes
          x1 = (-lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          x2 = (+lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          xlend = min(x1,x2)
          xuend = max(x1,x2)
c         --- Now, get distances to surface
          if (xx < xl) then
            if (xl <= xlend .and. xlend <= xu) delpx(i) = xlend - xx
            if (xlend <= xl .and. xl <= xuend) delpx(i) = xl - xx
          else if (xl <= xx .and. xx <= xu) then
            if (xl <= xlend .and. xlend <= xu .and. xx <= xuend) then
              if (xx <= xlend) then
                delpx(i) = xlend - xx
              else
                delpx(i) = xx - min(xu,xuend)
                delmx(i) = xlend - xx
              endif
            else if (xlend <= xl .and. xu <= xuend) then
              delmx(i) = xl - xx
              delpx(i) = xx - xu
            else if (xl <= xuend .and. xuend <= xu) then
              if (xx <= xuend) then
                delpx(i) = xx - xuend
                delmx(i) = max(xl,xlend) - xx
              else
                delmx(i) = xx - xuend
              endif
            endif
          else if (xu < xx) then
            if (xu >= xuend .and. xuend >= xl) delmx(i) = xx - xuend
            if (xuend >= xu .and. xu >= xlend) delmx(i) = xx - xu
          endif
        endif

c       --- delta y
        a = (cphi)**2
        b = -2*(zz*ctheta*sphi + xx*stheta*sphi)*cphi
        c = (xx*ctheta - zz*stheta)**2 + (zz*ctheta*sphi + xx*stheta*sphi)**2
     &      -rad**2
        if (b**2 - 4*a*c >= 0) then
          if (a == 0) then
            y1 = -c/b
            y2 = y1
          else
            ss = sqrt(max(0.,b**2 - 4*a*c))
            y1 = (-b-ss)/(2*a)
            y2 = (-b+ss)/(2*a)
          endif
          yl = min(y1,y2)
          yu = max(y1,y2)
c         --- Get the point where the z-line intersects the end planes
          y1 = (-lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          y2 = (+lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          ylend = min(y1,y2)
          yuend = max(y1,y2)
c         --- Now, get distances to surface
          if (yy < yl) then
            if (yl <= ylend .and. ylend <= yu) delpy(i) = ylend - yy
            if (ylend <= yl .and. yl <= yuend) delpy(i) = yl - yy
          else if (yl <= yy .and. yy <= yu) then
            if (yl <= ylend .and. ylend <= yu .and. yy <= yuend) then
              if (yy <= ylend) then
                delpy(i) = ylend - yy
              else
                delpy(i) = yy - min(yu,yuend)
                delmy(i) = ylend - yy
              endif
            else if (ylend <= yl .and. yu <= yuend) then
              delmy(i) = yl - yy
              delpy(i) = yy - yu
            else if (yl <= yuend .and. yuend <= yu) then
              if (yy <= yuend) then
                delpy(i) = yy - yuend
                delmy(i) = max(yl,ylend) - yy
              else
                delmy(i) = yy - yuend
              endif
            endif
          else if (yu < yy) then
            if (yu >= yuend .and. yuend >= yl) delmy(i) = yy - yuend
            if (yuend >= yu .and. yu >= ylend) delmy(i) = yy - yu
          endif
        endif

      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine ConeConductorF(r_zmin,r_zmax,length,theta,phi,
     &                          xcent,ycent,zcent,
     &                          n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                          fuzz)
      real(kind=8):: r_zmin,r_zmax,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a cylinder at an arbitrary
c angle.
c The notation is that a negative distance means that the point is inside.
c
c Write equation of surface as xx^2 + yy^2 = (r1 + (zz-z1)/(z2-z1)*(r2-r1))^2
c                              xx^2 + yy^2 = (s1 + s2*zz)^2
c Here, xx, yy, zz are in frame relative to cone axis. Conversion from
c frame of mesh, xg, yg, zg
c       xx = xg*cos(theta) - zg*sin(theta)
c       yy = yg*cos(phi) - zg*cos(theta)*sin(phi) - xg*sin(theta)*sin(phi)
c       zz = zg*cos(theta)*cos(phi) + xg*sin(theta)*cos(phi) + yg*sin(phi)

      integer(ISZ):: i
      real(kind=8):: ctheta,stheta,cphi,sphi,ss,cci,sci,si,lenhalf
      real(kind=8):: xx,yy,zz,zcross
      real(kind=8):: s1,s2,t1,t2,t3
      real(kind=8):: a,b,c
      real(kind=8):: z1,z2,zgcross,zl,zu,zlend,zuend,zi1,zi2
      real(kind=8):: x1,x2,xgcross,xl,xu,xlend,xuend,xi1,xi2
      real(kind=8):: y1,y2,ygcross,yl,yu,ylend,yuend,yi1,yi2

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)
      cci = 1./dvnz(ctheta*cphi)
      sci = 1./dvnz(stheta*cphi)
      si  = 1./dvnz(sphi)
      lenhalf = length/2.

      s1 = r_zmin - (-lenhalf)*(r_zmax - r_zmin)/length
      s2 = (r_zmax - r_zmin)/length
      zcross = -s1/dvnz(s2)

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent

c       --- delta z
        t1 = xx*ctheta
        t2 = yy*cphi - xx*stheta*sphi
        t3 = xx*stheta*cphi + yy*sphi
c       --- Coefficients of the quadratic
        a = stheta**2 + (ctheta*sphi)**2 - (s2*ctheta*cphi)**2
        b = -2*t1*stheta - 2*t2*ctheta*sphi - 2*(s1+s2*t3)*s2*ctheta*cphi
        c = t1**2 + t2**2 - (s1+s2*t3)**2
c       --- Check if z-line intersects the cone
        if (b**2 - 4*a*c >= -fuzz) then
c         --- Get the two points of intersection
          if (a == 0) then
            z1 = -c/b
            z2 = z1
          else
            ss = sqrt(max(0.,b**2 - 4*a*c))
            z1 = (-b-ss)/(2*a)
            z2 = (-b+ss)/(2*a)
          endif
          zl = min(z1,z2)
          zu = max(z1,z2)
c         --- Get the point where the z-line intersects the end planes
          z1 = (-lenhalf - xx*stheta*cphi - yy*sphi)*cci
          z2 = (+lenhalf - xx*stheta*cphi - yy*sphi)*cci
          zlend = min(z1,z2)
          zuend = max(z1,z2)
          zgcross = (zcross - xx*stheta*cphi - yy*sphi)*cci
c         --- Now, get distances to surface
          zi1 = +LARGEPOS
          zi2 = -LARGEPOS
          if (zl <= zgcross .and. zgcross <= zu .and.
     &        (zlend <= zl .or. zu <= zuend)) then
            if (zl < zlend) then
              zi1 = max(zu,zlend)
              zi2 = zuend
            else if (zu > zuend) then
              zi1 = zlend
              zi2 = min(zl,zuend)
            else if (zz <= zl) then
              zi1 = zlend
              zi2 = zl
            else if (zz <= zu) then
              zi2 = zl
              zi1 = zu
            else if (zu <= zz) then
              zi1 = zu
              zi2 = zuend
            endif
          else if (((zl <= zgcross .and. zu <= zgcross) .or.
     &              (zl >= zgcross .and. zu >= zgcross)) .and.
     &             (zl <= zuend .and. zlend <= zu)) then
            zi1 = max(zl,zlend)
            zi2 = min(zu,zuend)
          endif
            
          if (zz <= zi1) delpz(i) = zi1 - zz
          if (zi1 <= zz .and. zz <= zi2) then
            delmz(i) = zi1 - zz
            delpz(i) = zz - zi2
          endif
          if (zi2 <= zz) delmz(i) = zz - zi2
        endif


c       --- delta x
        t1 = zz*stheta
        t2 = yy*cphi - zz*ctheta*sphi
        t3 = zz*ctheta*cphi + yy*sphi
        a = (ctheta)**2 + (stheta*sphi)**2 - (s2*stheta*cphi)**2
        b = -2*t1*ctheta - 2*t2*stheta*sphi - 2*(s1+s2*t3)*s2*stheta*cphi
        c = t1**2 + t2**2 - (s1 + s2*t3)**2
        if (b**2 - 4*a*c >= 0) then
          if (a == 0) then
            x1 = -c/b
            x2 = x1
          else
            ss = sqrt(b**2 - 4*a*c)
            x1 = (-b-ss)/(2*a)
            x2 = (-b+ss)/(2*a)
          endif
          xl = min(x1,x2)
          xu = max(x1,x2)
c         --- Get the point where the x-line intersects the end planes
          x1 = (-lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          x2 = (+lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          xlend = min(x1,x2)
          xuend = max(x1,x2)
          xgcross = (zcross - zz*ctheta*cphi - yy*sphi)*sci
c         --- Now, get distances to surface
          xi1 = +LARGEPOS
          xi2 = -LARGEPOS
          if (xl < xgcross .and. xgcross < xu .and.
     &        (xlend <= xl .or. xu <= xuend)) then
            if (xl < xlend) then
              xi1 = max(xu,xlend)
              xi2 = xuend
            else if (xu > xuend) then
              xi1 = xlend
              xi2 = min(xl,xuend)
            else if (xx <= xl) then
              xi1 = xlend
              xi2 = xl
            else if (xx <= xu) then
              xi2 = xl
              xi1 = xu
            else if (xu <= xx) then
              xi1 = xu
              xi2 = xuend
            endif
          else if (((xl <= xgcross .and. xu <= xgcross) .or.
     &              (xl >= xgcross .and. xu >= xgcross)) .and.
     &             (xl <= xuend .and. xlend <= xu)) then
            xi1 = max(xl,xlend)
            xi2 = min(xu,xuend)
          endif
            
          if (xx <= xi1) delpx(i) = xi1 - xx
          if (xi1 <= xx .and. xx <= xi2) then
            delmx(i) = xi1 - xx
            delpx(i) = xx - xi2
          endif
          if (xi2 <= xx) delmx(i) = xx - xi2
        endif

c       --- delta y
        t1 = xx*ctheta - zz*stheta
        t2 = zz*ctheta*sphi + xx*stheta*sphi
        t3 = zz*ctheta*cphi + xx*stheta*cphi
        a = (cphi)**2 - (s2*sphi)**2
        b = -2*t2*cphi - 2*(s1+s2*t3)*s2*sphi
        c = t1**2 + t2**2 - (s1+s2*t3)**2
        if (b**2 - 4*a*c >= 0) then
          if (a == 0) then
            y1 = -c/b
            y2 = y1
          else
            ss = sqrt(b**2 - 4*a*c)
            y1 = (-b-ss)/(2*a)
            y2 = (-b+ss)/(2*a)
          endif
          yl = min(y1,y2)
          yu = max(y1,y2)
c         --- Get the point where the z-line intersects the end planes
          y1 = (-lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          y2 = (+lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          ylend = min(y1,y2)
          yuend = max(y1,y2)
          ygcross = (zcross - zz*ctheta*cphi - xx*stheta*cphi)*si
c         --- Now, get distances to surface
          yi1 = +LARGEPOS
          yi2 = -LARGEPOS
          if (yl < ygcross .and. ygcross < yu .and.
     &        (ylend <= yl .or. yu <= yuend)) then
            if (yl < ylend) then
              yi1 = max(yu,ylend)
              yi2 = yuend
            else if (yu > yuend) then
              yi1 = ylend
              yi2 = min(yl,yuend)
            else if (yy <= yl) then
              yi1 = ylend
              yi2 = yl
            else if (yy <= yu) then
              yi2 = yl
              yi1 = yu
            else if (yu <= yy) then
              yi1 = yu
              yi2 = yuend
            endif
          else if (((yl <= ygcross .and. yu <= ygcross) .or.
     &              (yl >= ygcross .and. yu >= ygcross)) .and.
     &             (yl <= yuend .and. ylend <= yu)) then
            yi1 = max(yl,ylend)
            yi2 = min(yu,yuend)
          endif
            
          if (yy <= yi1) delpy(i) = yi1 - yy
          if (yi1 <= yy .and. yy <= yi2) then
            delmy(i) = yi1 - yy
            delpy(i) = yy - yi2
          endif
          if (yi2 <= yy) delmy(i) = yy - yi2
        endif

      enddo

      return
      end
c=======================================================================
      subroutine BeamletplateConductorF(za,zb,z0,thickness,xcent,ycent,zcent,
     &                            n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                            fuzz)
      real(kind=8):: za,zb,z0,thickness,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a plate for the beamlet
c source pre-accelerator.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz
      real(kind=8):: x0l,xl1,xl2,xlave
      real(kind=8):: x0r,xr1,xr2,xrave
      real(kind=8):: y0l,yl1,yl2,ylave
      real(kind=8):: y0r,yr1,yr2,yrave
      real(kind=8):: z0l,zl1,zl2,zlave
      real(kind=8):: z0r,zr1,zr2,zrave

      do i=1,n

c       --- Set defaults
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent

c       --- Get z-intersect with left and right hand side of the plate
        z0l = z0 - 0.5*thickness
        zl1 = zb - sqrt((zb-za+sqrt((za-z0l)**2-xx**2))**2 - yy**2)
        zl2 = za - sqrt((za-zb+sqrt((zb-z0l)**2-yy**2))**2 - xx**2)
        zlave = 0.5*(zl1 + zl2)

        z0r = z0 + 0.5*thickness
        zr1 = zb - sqrt((zb-za+sqrt((za-z0r)**2-xx**2))**2 - yy**2)
        zr2 = za - sqrt((za-zb+sqrt((zb-z0r)**2-yy**2))**2 - xx**2)
        zrave = 0.5*(zr1 + zr2)

c       --- Point far to the left of the plate
        if (zz < z0l) then
          delpz(i) = zlave - zz
          cycle
        endif

c       --- Point to the left of plate
        if (zz <= zlave) then
          xl1 = ((za-z0l)**2 - (sqrt(yy**2+(zb-zz)**2)+za-zb)**2)
          yl1 = ((zb-za+sqrt((za-z0l)**2-xx**2))**2 - (zb-zz)**2)
          xl2 = ((za-zb+sqrt((zb-z0l)**2-yy**2))**2 - (za-zz)**2)
          yl2 = ((zb-z0l)**2 - (sqrt(xx**2+(za-zz)**2)+zb-za)**2)
          xlave = 0.5*(xl1 + xl2)
          ylave = 0.5*(yl1 + yl2)
          if (xlave < 0.) then
            xlave = -LARGEPOS
          else
            xlave = sqrt(xlave)
          endif
          if (ylave < 0.) then
            ylave = -LARGEPOS
          else
            ylave = sqrt(ylave)
          endif
          delmx(i) = xx - xlave
          delmy(i) = yy - ylave
          delpz(i) = zlave - zz
          cycle
        endif

c       --- Point in the plate
        if (zlave < zz .and. zz <= zrave) then
          xl1 = sqrt((za-z0l)**2 - (sqrt(yy**2+(zb-zz)**2)+za-zb)**2)
          yl1 = sqrt((zb-za+sqrt((za-z0l)**2-xx**2))**2 - (zb-zz)**2)
          xl2 = sqrt((za-zb+sqrt((zb-z0l)**2-yy**2))**2 - (za-zz)**2)
          yl2 = sqrt((zb-z0l)**2 - (sqrt(xx**2+(za-zz)**2)+zb-za)**2)
          xlave = 0.5*(xl1 + xl2)
          ylave = 0.5*(yl1 + yl2)

          xr1 = ((za-z0r)**2 - (sqrt(yy**2+(zb-zz)**2)+za-zb)**2)
          yr1 = ((zb-za+sqrt((za-z0r)**2-xx**2))**2 - (zb-zz)**2)
          xr2 = ((za-zb+sqrt((zb-z0r)**2-yy**2))**2 - (za-zz)**2)
          yr2 = ((zb-z0r)**2 - (sqrt(xx**2+(za-zz)**2)+zb-za)**2)
          xrave = 0.5*(xr1 + xr2)
          yrave = 0.5*(yr1 + yr2)
          if (xrave < 0.) then
            xrave = -LARGEPOS
          else
            xrave = sqrt(xrave)
          endif
          if (yrave < 0.) then
            yrave = -LARGEPOS
          else
            yrave = sqrt(yrave)
          endif

          delmx(i) = xrave - xx
          delpx(i) = xx - xlave
          delmy(i) = yrave - yy
          delpy(i) = yy - ylave
          delmz(i) = zlave - zz
          delpz(i) = zz - zrave
          cycle
        endif

c       --- Point to the right of plate
        if (zrave < zz) then
          xr1 = sqrt((za-z0r)**2 - (sqrt(yy**2+(zb-zz)**2)+za-zb)**2)
          yr1 = sqrt((zb-za+sqrt((za-z0r)**2-xx**2))**2 - (zb-zz)**2)
          xr2 = sqrt((za-zb+sqrt((zb-z0r)**2-yy**2))**2 - (za-zz)**2)
          yr2 = sqrt((zb-z0r)**2 - (sqrt(xx**2+(za-zz)**2)+zb-za)**2)
          xrave = 0.5*(xr1 + xr2)
          yrave = 0.5*(yr1 + yr2)
          delpx(i) = xrave - xx
          delpy(i) = yrave - yy
          delmz(i) = zz - zrave
          cycle
        endif

      enddo

      return
      end
c=======================================================================
