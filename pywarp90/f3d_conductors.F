#include "top.h"
c=======================================================================
c@(#) File F3D_CONDUCTORS.F, $Revision: 1.25 $
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
      subroutine getconductorsnextdir(position,direction,parity,newdirection)
      integer(ISZ):: position(0:2)
      integer(ISZ):: direction(0:2)
      integer(ISZ):: parity
      integer(ISZ):: newdirection(0:2)
      integer(ISZ):: iparity

      iparity = mod(sum(position) + parity + 1, 2)
      newdirection = 0
      if (direction(0) /= 0) then
        newdirection(1+iparity) = 1
      else if (direction(1) /= 0) then
        newdirection(2-2*iparity) = 1
      else if (direction(2) /= 0) then
        newdirection(iparity) = 1
      endif
      if (position(0) == 1) newdirection(0) = -newdirection(0)
      if (position(1) == 1) newdirection(1) = -newdirection(1)
      if (position(2) == 1) newdirection(2) = -newdirection(2)

      return
      end
c=======================================================================
      subroutine getconductorfacets(nc,icnd,dels,gridn,griddd,gridmin)
      use ConductorGeometryVisualization
      integer(ISZ):: nc
      integer(ISZ):: icnd(0:2,nc)
      real(kind=8):: dels(0:5,nc)
      integer(ISZ):: gridn(0:2)
      real(kind=8):: griddd(0:2)
      real(kind=8):: gridmin(0:2)

c Generates a triangulated surface from the subgrid data.
c That surface can be viewed using some visualization tool, such as OpenDX.

      integer(ISZ):: iii(0:gridn(0),0:gridn(1),0:gridn(2))
      integer(ISZ):: ic,ix,iy,iz,ox,oy,oz,it,t1,t2
      integer(ISZ):: oo(0:2)
      integer(ISZ):: npp
      real(kind=8):: pp(0:2,0:5)
      real(kind=8):: v1(0:2),v2(0:2),tnorm(0:2),nm

c     --- First, register all of the conductors in a mesh. This list will
c     --- be used to give the conductor number at each grid point.
      iii = 0
      do ic=1,nc
        iii(icnd(0,ic),icnd(1,ic),icnd(2,ic)) = ic
      enddo

c     --- Loop over all grid points
      do iz=0,gridn(2)-1
        do iy=0,gridn(1)-1
          do ix=0,gridn(0)-1

c         --- For each grid point, loop over the eight corners. This is
c         --- done so that no surfaces are missed.
          do oz=0,1
          do oy=0,1
          do ox=0,1

c           --- oo holds the offset relative to the current grid point.
            oo(0) = ox
            oo(1) = oy
            oo(2) = oz

c         --- If there is a conductor at this grid point, try to generate
c         --- the surface.
          if (iii(ix+ox,iy+oy,iz+oz) > 0) then
c           --- Generate the surface use one search order.
            call getconductorsnewfacet(ix,iy,iz,oo,1,gridn,iii,nc,dels,
     &                                 gridmin,griddd,pp,npp)
            t1 = 1
            t2 = 2
            if (npp == 0) then
c             --- If a surface was not found, try the other search order.
              call getconductorsnewfacet(ix,iy,iz,oo,0,gridn,iii,nc,dels,
     &                                   gridmin,griddd,pp,npp)
              t1 = 2
              t2 = 1
            endif
            if (npp > 0) then
c             --- If a surface was found, triangulate it.
              if (ntriangles + (npp-2) > maxtriangles) then
                maxtriangles = maxtriangles + 999
                call gchange("ConductorGeometryVisualization",0)
              endif
c             --- Triangles are generated from the starting point and each
c             --- subsequent pair of points. Note that the order in which
c             --- the points are saved is important. The order is switched
c             --- depending on the seach order used to find the surface
c             --- (the order is set by t1 and t2).
              do it=0,npp-3
                ntriangles = ntriangles + 1
                triangles(:,0,ntriangles) = pp(:,0)
                triangles(:,1,ntriangles) = pp(:,it+t1)
                triangles(:,2,ntriangles) = pp(:,it+t2)
c               --- For each triangle, the normal to it is saved.
                v1 = pp(:,it+t1) - pp(:,0)
                v2 = pp(:,it+t2) - pp(:,0)
                tnorm(0) = v1(1)*v2(2) - v1(2)*v2(1)
                tnorm(1) = v1(2)*v2(0) - v1(0)*v2(2)
                tnorm(2) = v1(0)*v2(1) - v1(1)*v2(0)
                nm = sqrt(sum(tnorm**2))
                tnorm = tnorm/dvnz(nm)
                normals(:,0,ntriangles) = tnorm
                normals(:,1,ntriangles) = tnorm
                normals(:,2,ntriangles) = tnorm
              enddo
            endif
          endif

          enddo
          enddo
          enddo

c         --- The calls getconductorsnewfacet change the sign of the iii data
c         --- to flag grid points which where included in the surfaces.
c         --- This line resets the signs.
          iii(ix:ix+1,iy:iy+1,iz:iz+1) = abs(iii(ix:ix+1,iy:iy+1,iz:iz+1))

          enddo
        enddo
      enddo

c     --- This is only for convenience.
      maxtriangles = ntriangles
      call gchange("ConductorGeometryVisualization",0)

      return
      end
c=======================================================================
      subroutine getconductorsnewfacet(ix,iy,iz,oo,parity,gridn,iii,nc,dels,
     &                                 gridmin,griddd,pp,npp)
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: oo(0:2)
      integer(ISZ):: parity
      integer(ISZ):: gridn(0:2)
      integer(ISZ):: iii(0:gridn(0),0:gridn(1),0:gridn(2))
      integer(ISZ):: nc
      real(kind=8):: dels(0:5,nc)
      real(kind=8):: gridmin(0:2)
      real(kind=8):: griddd(0:2)
      real(kind=8):: pp(0:2,0:5)
      integer(ISZ):: npp

c Work routine used by visualizeconductors. Given a grid point and an offset
c relative to that point, it searches for a closed surface in the grid cell.
c The parity determines the search order. This is important since there are
c some cases where a search in one order will not properly find the surface.
c (i.e. when the surface slices the cell in half, with four grid points on one
c side. With one direction of search, the algorithm gets stuck in a loop
c cycling about those four points without ever using the subgrid info. The
c other search direction works OK.)

      logical(ISZ):: done
      integer(ISZ):: dd(0:2),currentoo(0:2),currentdd(0:2)
      real(kind=8):: newpoint(0:2)
      integer(ISZ):: nextoo(0:2),nextdd(0:2)
      integer(ISZ):: ncorners,ic,id

c     --- Pick the starting direction to be along the x axis.
      dd = 0
      if (oo(0) == 0) dd(0) = +1
      if (oo(0) == 1) dd(0) = -1

      npp = 0
      done = .false.
      currentoo = oo
      currentdd = dd
      ncorners = 1
      do while (.not. done)
        ic = abs(iii(ix+currentoo(0),iy+currentoo(1),iz+currentoo(2)))
c       --- Flag this grid point as already checked.
        iii(ix+currentoo(0),iy+currentoo(1),iz+currentoo(2)) = -ic
c       --- Depending on the current search direction, set which
c       --- subgrid data is to be used.
        if (currentdd(0) == -1) id = 0
        if (currentdd(0) == +1) id = 1
        if (currentdd(1) == -1) id = 2
        if (currentdd(1) == +1) id = 3
        if (currentdd(2) == -1) id = 4
        if (currentdd(2) == +1) id = 5
        if (ic /= 0 .and. dels(id,ic) <= 1.000001) then
c         --- This next  point is part of the surface so include its data.
          newpoint(0) = ix + currentoo(0) + dels(id,ic)*currentdd(0)
          newpoint(1) = iy + currentoo(1) + dels(id,ic)*currentdd(1)
          newpoint(2) = iz + currentoo(2) + dels(id,ic)*currentdd(2)
          if (npp <= 5) pp(:,npp) = newpoint*griddd + gridmin
          npp = npp + 1
c         --- Get the next direction to check (relative to the same grid point).
          call getconductorsnextdir(currentoo,currentdd,parity,nextdd)
          nextoo = currentoo
        else
c         --- This point is not on the surface so get the next grid point to
c         --- check.
          nextoo = currentoo + currentdd
          call getconductorsnextdir(nextoo,-currentdd,parity,nextdd)
          ncorners = ncorners + 1
        endif
        currentoo = nextoo
        currentdd = nextdd
c       --- If the number of corners counter is too large, then the algorithm
c       --- is stuck in an infinite loop, so just quit.
        if (ncorners > 8) return
c       --- If the current point is the same as the original, then all grid
c       --- points have been checked, so return.
        if (all(currentoo == oo) .and. all(currentdd == dd)) return
c       --- If more than six points have been found on the surface, then
c       --- there is an error. Clear the point and return. Note that this
c       --- should never happen, but is here as a leftover from code
c       --- debugging. If is left in case some unusual data set foils
c       --- the algorithm.
        if (npp > 6) then
          npp = 0
          return
        endif
      enddo

      return
      end
c=======================================================================
      subroutine conductorsmoothshading(tt,ii)
      use ConductorGeometryVisualization
      real(kind=8):: tt(3*ntriangles)
      integer(ISZ):: ii(3*ntriangles)

c Modify the normals of the data so that a smoother surface will be drawn.
c For each set of point which have the same coordinates, this averages
c the normals.
c A much, much better way of doing this is to sort the data first.

      integer(ISZ):: ip,iq,i1,j1,i2,j2

      maxpoints = 3*ntriangles
      call gchange("ConductorGeometryVisualization")
      npoints = 0
      ip = 1
      do while (ip < 3*ntriangles)
        i1 = ii(ip)/3 + 1
        j1 = ii(ip) - i1*3
        npoints = npoints + 1
        points(:,npoints) = triangles(:,j1,i1)
        pnormals(:,npoints) = normals(:,j1,j1)
        connections(j1,i1) = npoints
        iq = ip + 1
        do while (abs(tt(iq) - tt(ip)) < 1.e-10)
          i2 = ii(iq)/3 + 1
          j2 = ii(iq) - i2*3
          if (all(abs(triangles(:,j1,i1)-triangles(:,j2,i2)) < 1.e-10)) then
            pnormals(:,npoints) = pnormals(:,npoints) + normals(:,j2,i2)
            connections(j2,i2) = npoints
          endif
          iq = iq + 1
        end do
        pnormals(:,npoints) = pnormals(:,npoints)/
     &                        dvnz(sqrt(sum(pnormals(:,npoints)**2)))
        ip = ip + 1
        i2 = ii(ip)/3 + 1
        j2 = ii(ip) - i2*3
        do while (all(abs(triangles(:,j1,i1)-triangles(:,j2,i2)) < 1.e-10))
          ip = ip + 1
          i2 = ii(ip)/3 + 1
          j2 = ii(ip) - i2*3
        end do
      end do


c     integer(ISZ):: i1,i2,j1,j2
c     real(kind=8):: nmsum(0:2)
c     integer(ISZ):: nn
c     integer(ISZ):: done(3,ntriangles)
c     integer(ISZ):: nlist,ilist(100),jlist(100),il

c     --- The 'done' array is used to keep track of which points have had
c     --- their normals averaged.
c     done = 0

c     --- Loop over all points of all triangles.
c     do i1=1,ntriangles
c       do j1=0,2

c         --- Skip points which have already been treated.
c         if (done(j1,i1) == 1) cycle

c         --- Begin accumulating the normals.
c         nmsum = normals(:,j1,i1)
c         nlist = 0

c         --- Loop over the points of the rest of the triangles.
c         do i2=i1+1,ntriangles
c           do j2=0,2

c             --- If the coordinates are the same (i.e. within some tolerance)
c             --- then add this point to the list.
c             if (all(abs(triangles(:,j1,i1)-triangles(:,j2,i2)) < 1.e-10)) then
c               done(j2,i2) = 1
c               nmsum = nmsum + normals(:,j2,i2)
c               nlist = nlist + 1
c               ilist(nlist) = i2
c               jlist(nlist) = j2
c             endif
c           enddo
c         enddo

c         --- Renormalize the summed normals are reset the normals in the list.
c         if (nlist > 0) then
c           nmsum = nmsum/dvnz(sqrt(sum(nmsum**2)))
c           normals(:,j1,i1) = nmsum
c           do il=1,nlist
c             normals(:,jlist(il),ilist(il)) = nmsum
c           enddo
c         endif

c       enddo
c     enddo

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
      subroutine ZPlaneConductorD(zcent,zsign,n,x,y,z,distance)
      real(kind=8):: zcent,zsign
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances to surface
c The notation is that a negative distance means that the point is inside.
c A positive zsign means that the conductor is the region above zcent.

      integer(ISZ):: i

      distance = zsign*(z - zcent)

      return
      end
c=======================================================================
      subroutine PlaneConductorF(z0,zsign,theta,phi,xcent,ycent,zcent,
     &                            n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                            fuzz)
      real(kind=8):: z0,zsign,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a plane
c The notation is that a negative distance means that the point is inside.
c A positive zsign means that the conductor is the region above zcent.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz
      real(kind=8):: xi,yi,zi,zp,zs
      real(kind=8):: ctheta,stheta,cphi,sphi,cci,sci,si

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)
c     --- These ensure that divisions by zero do not happen.
      if (ctheta == 0) ctheta = SMALLPOS
      if (stheta == 0) stheta = SMALLPOS
      if (cphi == 0) cphi = SMALLPOS
      if (sphi == 0) sphi = SMALLPOS
      cci = 1./(ctheta*cphi)
      sci = 1./(stheta*cphi)
      si  = 1./(sphi)

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent

        xi = (z0 - zz*ctheta*cphi - yy*sphi)*sci
        yi = (z0 - xx*stheta*cphi - zz*ctheta*cphi)*si
        zi = (z0 - xx*stheta*cphi - yy*sphi)*cci
        zp = xx*stheta*cphi + yy*sphi + zz*ctheta*cphi
        if (zp >= z0) then
          zs = -zsign
        else
          zs = +zsign
        endif

        if (xx <= xi) then
          delmx(i) = LARGEPOS*zs
          delpx(i) = abs(xi - xx)*zs
        else
          delmx(i) = abs(xi - xx)*zs
          delpx(i) = LARGEPOS*zs
        endif
        if (yy <= yi) then
          delmy(i) = LARGEPOS*zs
          delpy(i) = abs(yi - yy)*zs
        else
          delmy(i) = abs(yi - yy)*zs
          delpy(i) = LARGEPOS*zs
        endif
        if (zz <= zi) then
          delmz(i) = LARGEPOS*zs
          delpz(i) = abs(zi - zz)*zs
        else
          delmz(i) = abs(zi - zz)*zs
          delpz(i) = LARGEPOS*zs
        endif

      enddo

      return
      end
c=======================================================================
      subroutine PlaneConductorD(z0,zsign,theta,phi,xcent,ycent,zcent,
     &                           n,x,y,z,distance)
      real(kind=8):: z0,zsign,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances to surface
c The notation is that a negative distance means that the point is inside.
c A positive zsign means that the conductor is the region above zcent.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz
      real(kind=8):: zp

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        zp = xx*sin(theta)*cos(phi) + yy*sin(phi) + zz*cos(theta)*cos(phi)
        distance(i) = (z0 - zp)*zsign

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

c Finds distances along the axis of the points from a box
c The notation is that a negative distance means that the point is inside.

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
      subroutine BoxConductorD(xsize,ysize,zsize,xcent,ycent,zcent,
     &                         n,x,y,z,distance)
      real(kind=8):: xsize,ysize,zsize,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a box.
c The notation is that a negative distance means that the point is inside.
c A positive zsign means that the conductor is the region above zcent.

      integer(ISZ):: i
      real(kind=8):: xd,yd,zd,dd

      do i=1,n
        xd = abs(x(i) - xcent) - xsize/2.
        yd = abs(y(i) - ycent) - ysize/2.
        zd = abs(z(i) - zcent) - zsize/2.

        dd = sqrt(max(0.,xd)**2 + max(0.,yd)**2 + max(0.,zd)**2)

        if (dd == 0.) dd = max(xd,yd,zd)

        distance(i) = dd
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
      subroutine ZCylinderConductorD(rad,length,xcent,ycent,zcent,
     &                               n,x,y,z,distance)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder aligned with the z-axis.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: rr,zz,dd

      do i=1,n
        rr = sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2) - rad
        zz = abs(z(i) - zcent) - length/2.

        dd = sqrt(max(0.,rr)**2 + max(0.,zz)**2)
        if (dd == 0.) dd = max(rr,zz)
        distance(i) = dd
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
      subroutine ZRoundedCylinderConductorD(rad,length,rad2,xcent,ycent,zcent,
     &                                      n,x,y,z,distance)
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder aligned with
c the z-axis. The corners of the cylinder are rounded with radius rad2.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: z2,rr,zz

      z2 = max(length/2.-rad2,0.)

      do i=1,n

        rr = sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2)
        zz = abs(z(i) - zcent)
        if (zz >= z2 .and. rr <= rad-rad2) then
          distance(i) = zz - length/2.
        else if (zz <= z2) then
          distance(i) = max(rr - rad,zz - length/2.)
        else
          distance(i) = sqrt((zz - z2)**2 + (rr - (rad - rad2))**2) - rad2
        endif
          
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
      subroutine ZCylinderOutConductorD(rad,length,xcent,ycent,zcent,
     &                                  n,x,y,z,distance)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder aligned with the z-axis.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: rr,zz,dd

      do i=1,n
        rr = rad - sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2)
        zz = abs(z(i) - zcent) - length/2.

        dd = sqrt(max(0.,rr)**2 + max(0.,zz)**2)
        if (dd == 0.) dd = max(rr,zz)
        distance(i) = dd
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
      subroutine ZRoundedCylinderOutConductorD(rad,length,rad2,
     &                                         xcent,ycent,zcent,
     &                                         n,x,y,z,distance)
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder aligned with
c the z-axis. The corners of the cylinder are rounded with radius rad2.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: z2,rr,zz

      z2 = max(length/2.-rad2,0.)

      do i=1,n

        rr = sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2)
        zz = abs(z(i) - zcent)
        if (zz >= z2 .and. rr >= rad+rad2) then
          distance(i) = zz - length/2.
        else if (zz <= z2) then
          distance(i) = max(rad - rr,zz - length/2.)
        else
          distance(i) = sqrt((zz - z2)**2 + (rr - (rad + rad2))**2) - rad2
        endif
          
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
      subroutine XCylinderConductorD(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,distance)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder aligned with the x-axis.
c The notation is that a negative distance means that the point is inside.

      call ZCylinderConductorD(rad,length,ycent,zcent,xcent,
     &                         n,y,z,x,distance)

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
      subroutine YCylinderConductorD(rad,length,xcent,ycent,zcent,
     &                               n,x,y,z,distance)
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder aligned with the y-axis.
c The notation is that a negative distance means that the point is inside.

      call ZCylinderConductorD(rad,length,zcent,xcent,ycent,
     &                         n,z,x,y,distance)

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
      subroutine SphereConductorD(rad,xcent,ycent,zcent,n,x,y,z,distance)
      real(kind=8):: rad,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a sphere.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2 + zz**2)
        
        distance(i) = rr - rad

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
      subroutine ZConeConductorD(r_zmin,r_zmax,length,xcent,ycent,zcent,
     &                           n,x,y,z,distance)
      use Constant
      real(kind=8):: r_zmin,r_zmax,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from a cone aligned with
c the z-axis.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr
      real(kind=8):: d1,ang_zmin,ang_zmax,dd,dz

      do i=1,n
        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        d1 = sqrt((rr - r_zmin)**2 + (zz + length/2.)**2)
        ang_zmin = atan2((rr - r_zmin),(zz + length/2.)) -
     &         atan2((r_zmax - r_zmin),length)
        if (ang_zmin > +pi) ang_zmin = ang_zmin - 2*pi
        if (ang_zmin < -pi) ang_zmin = ang_zmin + 2*pi
        ang_zmax = atan2((r_zmin - r_zmax),-length) -
     &          atan2((rr - r_zmax),(zz - length/2.))
        if (ang_zmax > +pi) ang_zmax = ang_zmax - 2*pi
        if (ang_zmax < -pi) ang_zmax = ang_zmax + 2*pi
        dd = d1*sin(ang_zmin)
        dz = abs(zz) - length/2.

        if (zz <= -length/2. .and. rr <= r_zmin) then
          distance(i) = dz
        else if (ang_zmin >= pi/2.) then
          distance(i) = sqrt((rr - r_zmin)**2 + (zz + length/2.)**2)
        else if (zz >= length/2. .and. rr <= r_zmax) then
          distance(i) = dz
        else if (ang_zmax >= pi/2.) then
          distance(i) = sqrt((rr - r_zmax)**2 + (zz - length/2.)**2)
        else if (dd >= 0.) then
          distance(i) = dd
        else
          distance(i) = max(dd,dz)
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
      subroutine ZConeOutConductorD(r_zmin,r_zmax,length,xcent,ycent,zcent,
     &                              n,x,y,z,distance)
      use Constant
      real(kind=8):: r_zmin,r_zmax,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from a cone aligned with
c the z-axis.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr
      real(kind=8):: d1,ang_zmin,ang_zmax,dd,dz

      do i=1,n
        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        d1 = sqrt((rr - r_zmin)**2 + (zz + length/2.)**2)
        ang_zmin = atan2((r_zmax - r_zmin),length) -
     &         atan2((rr - r_zmin),(zz + length/2.))
        if (ang_zmin > +pi) ang_zmin = ang_zmin - 2*pi
        if (ang_zmin < -pi) ang_zmin = ang_zmin + 2*pi
        ang_zmax = atan2((rr - r_zmax),(zz - length/2.)) -
     &          atan2((r_zmin - r_zmax),-length)
        if (ang_zmax > +pi) ang_zmax = ang_zmax - 2*pi
        if (ang_zmax < -pi) ang_zmax = ang_zmax + 2*pi
        dd = d1*sin(ang_zmin)
        dz = abs(zz) - length/2.

        if (zz <= -length/2. .and. rr >= r_zmin) then
          distance(i) = dz
        else if (ang_zmin >= pi/2.) then
          distance(i) = sqrt((rr - r_zmin)**2 + (zz + length/2.)**2)
        else if (zz >= length/2. .and. rr >= r_zmax) then
          distance(i) = dz
        else if (ang_zmax >= pi/2.) then
          distance(i) = sqrt((rr - r_zmax)**2 + (zz - length/2.)**2)
        else if (dd >= 0.) then
          distance(i) = dd
        else
          distance(i) = max(dd,dz)
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
      subroutine ZTorusConductorD(r1,r2,xcent,ycent,zcent,
     &                            n,x,y,z,distance)
      real(kind=8):: r1,r2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a torus aligned with the z-axis.
c r1 is radius relative to z-axis
c r2 is radius relative to theta-axis
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        distance(i) = sqrt(zz**2 + (r1 - rr)**2) - r2

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
      real(kind=8):: z1,z2,zl,zu,zlend,zuend,zi1,zi2
      real(kind=8):: x1,x2,xl,xu,xlend,xuend,xi1,xi2
      real(kind=8):: y1,y2,yl,yu,ylend,yuend,yi1,yi2

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)
c     --- In the limiting cases with the angles being exactly zero (or p/2)
c     --- the algorithm breaks down. The sign of expressions dependent on
c     --- the angles is important, and is lost when the angles are exactly
c     --- zero (or pi/2). So the limiting cases is used where the angle
c     --- approaches but does not equal to zero (or pi/2).
      if (ctheta == 0) ctheta = SMALLPOS
      if (stheta == 0) stheta = SMALLPOS
      if (cphi == 0) cphi = SMALLPOS
      if (sphi == 0) sphi = SMALLPOS
      cci = 1./(ctheta*cphi)
      sci = 1./(stheta*cphi)
      si  = 1./(sphi)
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
          zi1 = +LARGEPOS
          zi2 = -LARGEPOS
          if (zl <= zuend .and. zlend <= zu) then
            zi1 = max(zl,zlend)
            zi2 = min(zu,zuend)
          endif
          if (zz <= zi1) then
            delpz(i) = zi1 - zz
          else if (zi1 <= zz .and. zz <= zi2) then
            delmz(i) = zi1 - zz
            delpz(i) = zz - zi2
          else if (zi2 <= zz) then
            delmz(i) = zz - zi2
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
          xi1 = +LARGEPOS
          xi2 = -LARGEPOS
          if (xl <= xuend .and. xlend <= xu) then
            xi1 = max(xl,xlend)
            xi2 = min(xu,xuend)
          endif
          if (xx <= xi1) then
            delpx(i) = xi1 - xx
          else if (xi1 <= xx .and. xx <= xi2) then
            delmx(i) = xi1 - xx
            delpx(i) = xx - xi2
          else if (xi2 <= xx) then
            delmx(i) = xx - xi2
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
          yi1 = +LARGEPOS
          yi2 = -LARGEPOS
          if (yl <= yuend .and. ylend <= yu) then
            yi1 = max(yl,ylend)
            yi2 = min(yu,yuend)
          endif
          if (yy <= yi1) then
            delpy(i) = yi1 - yy
          else if (yi1 <= yy .and. yy <= yi2) then
            delmy(i) = yi1 - yy
            delpy(i) = yy - yi2
          else if (yi2 <= yy) then
            delmy(i) = yy - yi2
          endif
        endif

      enddo

      return
      end
c=======================================================================
      subroutine CylinderConductorD(rad,length,theta,phi,xcent,ycent,zcent,
     &                              n,x,y,z,distance)
      real(kind=8):: rad,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cylinder at an arbitrary angle.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: ctheta,stheta,cphi,sphi
      real(kind=8):: xx,yy,zz,xp,yp,zp,rr,dd

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        xp = xx*ctheta - zz*stheta
        yp = yy*cphi - zz*ctheta*sphi - xx*stheta*sphi
        zp = zz*ctheta*cphi + xx*stheta*cphi + yy*sphi

        rr = sqrt(xp**2 + yp**2) - rad
        zz = abs(zp) - length/2.

        dd = sqrt(max(0.,rr)**2 + max(0.,zz)**2)
        if (dd == 0.) dd = max(rr,zz)
        distance(i) = dd
      enddo

      return
      end
c=======================================================================
      subroutine CylindersConductorF(ncylinders,rad,length,theta,phi,
     &                           xcent,ycent,zcent,
     &                           n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                           fuzz)
      use Constant
      integer(ISZ):: ncylinders
      real(kind=8):: rad(ncylinders),length(ncylinders)
      real(kind=8):: theta(ncylinders),phi(ncylinders)
      real(kind=8):: xcent(ncylinders),ycent(ncylinders),zcent(ncylinders)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a list of cylinders at
c arbitrary angles.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: nx,ny,ic
      real(kind=8):: dx,dy
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: xc(1),yc(1),zc(1),zz
      real(kind=8):: mxc(1),pxc(1),myc(1),pyc(1),mzc(1),pzc(1)
      integer(ISZ):: nt,ixmin,ixmax,iymin,iymax,tix,tiy,ii,id
      real(kind=8):: xt(n),yt(n),zt(n)
      real(kind=8):: mxt(n),pxt(n),myt(n),pyt(n),mzt(n),pzt(n)

c     --- Very kludgy code to extract number of grid points and cell
c     --- sizes from the input coordinates.
      ymin = minval(y)
      ymax = maxval(y)
      dy = y(2) - y(1)
      if (dy /= 0.) then
        ny = maxval(int((y - ymin)/dy))
      else
        dy = 1.
        ny = 0
      endif
      
      dx = x(ny+2) - x(1)
      xmin = minval(x)
      xmax = maxval(x)
      nx = maxval(int((x - xmin)/dx))

      delmx = LARGEPOS
      delpx = LARGEPOS
      delmy = LARGEPOS
      delpy = LARGEPOS
      delmz = LARGEPOS
      delpz = LARGEPOS

c     --- This code assumes small angles
      if (maxval(abs(theta)) > pi/4. .or. maxval(abs(phi)) > pi/4.) then
        print*,"Error: The Cylinders object assumes small angles."
        print*,"Cylinders with large angles must be treated independently"
        stop
      endif

      do ic=1,ncylinders

c       --- First, get extent of cylinder by getting deltas relative to
c       --- center.
c       --- Get x and y at center of object at z = z(1)
        zz = z(1) - zcent(ic)
        xc(1) = zz*tan(theta(ic)) + xcent(ic)
        yc(1) = zz*(cos(theta(ic))*sin(phi(ic)) +
     &              tan(theta(ic))*sin(theta(ic))*sin(phi(ic)))/cos(phi(ic)) +
     &          ycent(ic)
        zc(1) = z(1)
        call CylinderConductorF(rad(ic),length(ic),theta(ic),phi(ic),
     &                          xcent(ic),ycent(ic),zcent(ic),
     &                          1,xc,yc,zc,mxc,pxc,myc,pyc,mzc,pzc,
     &                          fuzz)
c       --- Check if point outside of object.
        if (mxc(1) > 0.) then
c         --- Check distance in z to object. If the distance to left or right
c         --- is less than LARGEPOS, then set z to that location and find
c         --- x and y deltas there. Otherwise this conductor can be skipped.
          if (mzc(1) < LARGEPOS) then
            zz = z(1) - mzc(1)*1.000001
          else if (pzc(1) < LARGEPOS) then
            zz = z(1) + pzc(1)*1.000001
          else
            cycle
          endif
c         xc(1) = zz*tan(theta(ic)) + xcent(ic)
c         yc(1) = zz*(cos(theta(ic))*sin(phi(ic)) +
c    &              tan(theta(ic))*sin(theta(ic))*sin(phi(ic)))/cos(phi(ic)) +
c    &            ycent(ic)
          zc(1) = zz
          call CylinderConductorF(rad(ic),length(ic),theta(ic),phi(ic),
     &                        xcent(ic),ycent(ic),zcent(ic),
     &                        1,xc,yc,zc,mxc,pxc,myc,pyc,mzc,pzc,fuzz)
        endif
c       --- Get extent relative to grid
        ixmin = max(0 ,int((xcent(ic) - (-mxc(1)) - xmin)/dx) - 1)
        ixmax = min(nx,int((xcent(ic) + (-pxc(1)) - xmin)/dx) + 2)
        iymin = max(0 ,int((ycent(ic) - (-myc(1)) - ymin)/dy) - 1)
        iymax = min(ny,int((ycent(ic) + (-pxc(1)) - ymin)/dy) + 2)
c       --- Generate grid points within the circumscribing rectangle
        ii = 0
        do tix = ixmin,ixmax
          do tiy = iymin,iymax
            ii = ii + 1
            xt(ii) = xmin + tix*dx
            yt(ii) = ymin + tiy*dy
            zt(ii) = z(1)
          enddo
        enddo

c       --- Get deltas for points in that rectangle
        call CylinderConductorF(rad(ic),length(ic),theta(ic),phi(ic),
     &                          xcent(ic),ycent(ic),zcent(ic),
     &                          ii,xt,yt,zt,mxt,pxt,myt,pyt,mzt,pzt,
     &                          fuzz)

c       --- Fold the data into the arrays
        ii = 0
        do tix = ixmin,ixmax
          do tiy = iymin,iymax
            ii = ii + 1
            id = tix*(ny+1) + tiy + 1
            delmx(id) = min(delmx(id),mxt(ii))
            delpx(id) = min(delpx(id),pxt(ii))
            delmy(id) = min(delmy(id),myt(ii))
            delpy(id) = min(delpy(id),pyt(ii))
            delmz(id) = min(delmz(id),mzt(ii))
            delpz(id) = min(delpz(id),pzt(ii))
          enddo
        enddo

      enddo

      return
      end
c=======================================================================
      subroutine CylindersConductorD(ncylinders,rad,length,theta,phi,
     &                               xcent,ycent,zcent,n,x,y,z,distance)
      use Constant
      integer(ISZ):: ncylinders
      real(kind=8):: rad(ncylinders),length(ncylinders)
      real(kind=8):: theta(ncylinders),phi(ncylinders)
      real(kind=8):: xcent(ncylinders),ycent(ncylinders),zcent(ncylinders)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a list of cylinders at arbitrary angles.
c The notation is that a negative distance means that the point is inside.

c The very kludgy code used to calculate the deltas won't work here
c since the coordinates passed in may not be regularly spaced.
c This code assumes that the cylinders do not overlap transversely.
c Otherwise, the cylinders can be at any angle.

      real(kind=8):: dd
      integer(ISZ):: i,imin(1),ic

      distance = LARGEPOS
      do i=1,n
        imin = minloc(abs(x(i) - xcent) + abs(y(i) - ycent))
        ic = imin(1)
        call CylinderConductorD(rad(ic),length(ic),theta(ic),phi(ic),
     &                          xcent(ic),ycent(ic),zcent(ic),
     &                          1,x(i),y(i),z(i),dd)
        distance(i) = min(distance(i),dd)
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

c Finds distances along the axis of the points from a cone at an arbitrary
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
c     --- In the limiting cases with the angles being exactly zero (or p/2)
c     --- the algorithm breaks down. The sign of expressions dependent on
c     --- the angles is important, and is lost when the angles are exactly
c     --- zero (or pi/2). So the limiting cases is used where the angle
c     --- approaches but does not equal to zero (or pi/2).
      if (ctheta == 0) ctheta = SMALLPOS
      if (stheta == 0) stheta = SMALLPOS
      if (cphi == 0) cphi = SMALLPOS
      if (sphi == 0) sphi = SMALLPOS
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
c         --- Get the point where the z-line intersects the end planes
          z1 = (-lenhalf - xx*stheta*cphi - yy*sphi)*cci
          z2 = (+lenhalf - xx*stheta*cphi - yy*sphi)*cci
          zlend = min(z1,z2)
          zuend = max(z1,z2)
          zgcross = (zcross - xx*stheta*cphi - yy*sphi)*cci
c         --- Get the two points of intersection
          if (a == 0) then
            z1 = -c/b
            z2 = z1
          else
            ss = sqrt(max(0.,b**2 - 4*a*c))
            z1 = (-b-ss)/(2*a)
            z2 = (-b+ss)/(2*a)
            if (ss == 0. .and. abs(z1 - zgcross) < fuzz) then
              z1 = zgcross
              z2 = zgcross
            endif
          endif
          zl = min(z1,z2)
          zu = max(z1,z2)
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
c         --- Get the point where the x-line intersects the end planes
          x1 = (-lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          x2 = (+lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          xlend = min(x1,x2)
          xuend = max(x1,x2)
          xgcross = (zcross - zz*ctheta*cphi - yy*sphi)*sci
c         --- Get the two points of intersection
          if (a == 0) then
            x1 = -c/b
            x2 = x1
          else
            ss = sqrt(b**2 - 4*a*c)
            x1 = (-b-ss)/(2*a)
            x2 = (-b+ss)/(2*a)
            if (ss == 0. .and. abs(x1 - xgcross) < fuzz) then
              x1 = xgcross
              x2 = xgcross
            endif
          endif
          xl = min(x1,x2)
          xu = max(x1,x2)
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
c         --- Get the point where the z-line intersects the end planes
          y1 = (-lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          y2 = (+lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          ylend = min(y1,y2)
          yuend = max(y1,y2)
          ygcross = (zcross - zz*ctheta*cphi - xx*stheta*cphi)*si
c         --- Get the two points of intersection
          if (a == 0) then
            y1 = -c/b
            y2 = y1
          else
            ss = sqrt(b**2 - 4*a*c)
            y1 = (-b-ss)/(2*a)
            y2 = (-b+ss)/(2*a)
            if (ss == 0. .and. abs(y1 - ygcross) < fuzz) then
              y1 = ygcross
              y2 = ygcross
            endif
          endif
          yl = min(y1,y2)
          yu = max(y1,y2)
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
      subroutine ConeConductorD(r_zmin,r_zmax,length,theta,phi,
     &                          xcent,ycent,zcent,n,x,y,z,distance)
      use Constant
      real(kind=8):: r_zmin,r_zmax,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a cone at an arbitrary angle.
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
      real(kind=8):: ctheta,stheta,cphi,sphi
      real(kind=8):: xx,yy,zz,xp,yp,zp,rp
      real(kind=8):: d1,ang_zmin,ang_zmax,dd,dz

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        xp = xx*ctheta - zz*stheta
        yp = yy*cphi - zz*ctheta*sphi - xx*stheta*sphi
        zp = zz*ctheta*cphi + xx*stheta*cphi + yy*sphi
        rp = sqrt(xp**2 + yp**2)

        d1 = sqrt((rp - r_zmin)**2 + (zp + length/2.)**2)
        ang_zmin = atan2((rp - r_zmin),(zp + length/2.)) -
     &         atan2((r_zmax - r_zmin),length)
        if (ang_zmin > +pi) ang_zmin = ang_zmin - 2*pi
        if (ang_zmin < -pi) ang_zmin = ang_zmin + 2*pi
        ang_zmax = atan2((r_zmin - r_zmax),-length) -
     &          atan2((rp - r_zmax),(zp - length/2.))
        if (ang_zmax > +pi) ang_zmax = ang_zmax - 2*pi
        if (ang_zmax < -pi) ang_zmax = ang_zmax + 2*pi
        dd = d1*sin(ang_zmin)
        dz = abs(zp) - length/2.

        if (zp <= -length/2. .and. rp <= r_zmin) then
          distance(i) = dz
        else if (ang_zmin >= pi/2.) then
          distance(i) = sqrt((rp - r_zmin)**2 + (zp + length/2.)**2)
        else if (zp >= length/2. .and. rp <= r_zmax) then
          distance(i) = dz
        else if (ang_zmax >= pi/2.) then
          distance(i) = sqrt((rp - r_zmax)**2 + (zp - length/2.)**2)
        else if (dd >= 0.) then
          distance(i) = dd
        else
          distance(i) = max(dd,dz)
        endif
      enddo

      return
      end
c=======================================================================
      subroutine ConesConductorF(ncones,r_zmin,r_zmax,length,theta,phi,
     &                           xcent,ycent,zcent,
     &                           n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                           fuzz)
      use Constant
      integer(ISZ):: ncones
      real(kind=8):: r_zmin(ncones),r_zmax(ncones),length(ncones)
      real(kind=8):: theta(ncones),phi(ncones)
      real(kind=8):: xcent(ncones),ycent(ncones),zcent(ncones)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a list of cones at
c arbitrary angles.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: nx,ny,ic
      real(kind=8):: dx,dy
      real(kind=8):: xmin,xmax,ymin,ymax
      real(kind=8):: xc(1),yc(1),zc(1),zz
      real(kind=8):: mxc(1),pxc(1),myc(1),pyc(1),mzc(1),pzc(1)
      integer(ISZ):: nt,ixmin,ixmax,iymin,iymax,tix,tiy,ii,id
      real(kind=8):: xt(n),yt(n),zt(n)
      real(kind=8):: mxt(n),pxt(n),myt(n),pyt(n),mzt(n),pzt(n)

c     --- Very kludgy code to extract number of grid points and cell
c     --- sizes from the input coordinates.
      ymin = minval(y)
      ymax = maxval(y)
      dy = y(2) - y(1)
      if (dy /= 0.) then
        ny = maxval(int((y - ymin)/dy))
      else
        dy = 1.
        ny = 0
      endif
      
      xmin = minval(x)
      xmax = maxval(x)
      dx = x(ny+2) - x(1)
      nx = maxval(int((x - xmin)/dx))

      delmx = LARGEPOS
      delpx = LARGEPOS
      delmy = LARGEPOS
      delpy = LARGEPOS
      delmz = LARGEPOS
      delpz = LARGEPOS

c     --- This code assumes small angles
      if (maxval(abs(theta)) > pi/4. .or. maxval(abs(phi)) > pi/4.) then
        print*,"Error: The Cones object assumes small angles."
        print*,"Cones with large angles must be treated independently"
        stop
      endif

      do ic=1,ncones

c       --- First, get extent of cone by getting deltas relative to
c       --- center.
c       --- Get x and y at center of object at z = z(1)
        zz = z(1) - zcent(ic)
        xc(1) = zz*tan(theta(ic)) + xcent(ic)
        yc(1) = zz*(cos(theta(ic))*sin(phi(ic)) +
     &              tan(theta(ic))*sin(theta(ic))*sin(phi(ic)))/cos(phi(ic)) +
     &          ycent(ic)
        zc(1) = z(1)
        call ConeConductorF(r_zmin(ic),r_zmax(ic),length(ic),theta(ic),phi(ic),
     &                      xcent(ic),ycent(ic),zcent(ic),
     &                      1,xc,yc,zc,mxc,pxc,myc,pyc,mzc,pzc,
     &                      fuzz)
c       --- Check if point outside of object.
        if (mxc(1) > 0.) then
c         --- Check distance in z to object. If the distance to left or right
c         --- is less than LARGEPOS, then set z to that location and find
c         --- x and y deltas there. Otherwise this conductor can be skipped.
          if (mzc(1) < LARGEPOS) then
            zz = z(1) - mzc(1)*1.000001
          else if (pzc(1) < LARGEPOS) then
            zz = z(1) + pzc(1)*1.000001
          else
            cycle
          endif
c         xc(1) = zz*tan(theta(ic)) + xcent(ic)
c         yc(1) = zz*(cos(theta(ic))*sin(phi(ic)) +
c    &              tan(theta(ic))*sin(theta(ic))*sin(phi(ic)))/cos(phi(ic)) +
c    &            ycent(ic)
          zc(1) = zz
          call ConeConductorF(r_zmin(ic),r_zmax(ic),length(ic),
     &                        theta(ic),phi(ic),xcent(ic),ycent(ic),zcent(ic),
     &                        1,xc,yc,zc,mxc,pxc,myc,pyc,mzc,pzc,fuzz)
        endif
c       --- Get extent relative to grid
        ixmin = max(0 ,int((xcent(ic) - (-mxc(1)) - xmin)/dx) - 1)
        ixmax = min(nx,int((xcent(ic) + (-pxc(1)) - xmin)/dx) + 2)
        iymin = max(0 ,int((ycent(ic) - (-myc(1)) - ymin)/dy) - 1)
        iymax = min(ny,int((ycent(ic) + (-pxc(1)) - ymin)/dy) + 2)
c       --- Generate grid points within the circumscribing rectangle
        ii = 0
        do tix = ixmin,ixmax
          do tiy = iymin,iymax
            ii = ii + 1
            xt(ii) = xmin + tix*dx
            yt(ii) = ymin + tiy*dy
            zt(ii) = z(1)
          enddo
        enddo

c       --- Get deltas for points in that rectangle
        call ConeConductorF(r_zmin(ic),r_zmax(ic),length(ic),theta(ic),phi(ic),
     &                      xcent(ic),ycent(ic),zcent(ic),
     &                      ii,xt,yt,zt,mxt,pxt,myt,pyt,mzt,pzt,
     &                      fuzz)

c       --- Fold the data into the arrays
        ii = 0
        do tix = ixmin,ixmax
          do tiy = iymin,iymax
            ii = ii + 1
            id = tix*(ny+1) + tiy + 1
            delmx(id) = min(delmx(id),mxt(ii))
            delpx(id) = min(delpx(id),pxt(ii))
            delmy(id) = min(delmy(id),myt(ii))
            delpy(id) = min(delpy(id),pyt(ii))
            delmz(id) = min(delmz(id),mzt(ii))
            delpz(id) = min(delpz(id),pzt(ii))
          enddo
        enddo

      enddo

      return
      end
c=======================================================================
      subroutine ConesConductorD(ncones,r_zmin,r_zmax,length,theta,phi,
     &                           xcent,ycent,zcent,n,x,y,z,distance)
      use Constant
      integer(ISZ):: ncones
      real(kind=8):: r_zmin(ncones),r_zmax(ncones),length(ncones)
      real(kind=8):: theta(ncones),phi(ncones)
      real(kind=8):: xcent(ncones),ycent(ncones),zcent(ncones)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a list of cones at arbitrary angles.
c The notation is that a negative distance means that the point is inside.

c The very kludgy code used to calculate the deltas won't work here
c since the coordinates passed in may not be regularly spaced.
c This code assumes that the cones do not overlap transversely.
c Otherwise, the cones can be at any angle.

      real(kind=8):: dd
      integer(ISZ):: i,imin(1),ic

      distance = LARGEPOS
      do i=1,n
        imin = minloc(abs(x(i) - xcent) + abs(y(i) - ycent))
        ic = imin(1)
        call ConeConductorD(r_zmin(ic),r_zmax(ic),length(ic),theta(ic),phi(ic),
     &                      xcent(ic),ycent(ic),zcent(ic),
     &                      1,x(i),y(i),z(i),dd)
        distance(i) = min(distance(i),dd)
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
      real(kind=8):: xx,yy,zz,aa,bb
      real(kind=8):: x0l,xl1sq,xl2sq,xl1,xl2,xlave
      real(kind=8):: x0r,xr1sq,xr2sq,xr1,xr2,xrave
      real(kind=8):: y0l,yl1sq,yl2sq,yl1,yl2,ylave
      real(kind=8):: y0r,yr1sq,yr2sq,yr1,yr2,yrave
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
        aa = (2.*z0l*za - z0l**2 + xx**2)/(za+sqrt((za-z0l)**2-xx**2))
        zl1 = (2.*aa*zb - aa**2 + yy**2)/(zb+sqrt((zb-aa)**2-yy**2))
        bb = (2.*z0l*zb - z0l**2 + yy**2)/(zb+sqrt((zb-z0l)**2-yy**2))
        zl2 = (2.*bb*za - bb**2 + xx**2)/(za+sqrt((za-bb)**2-xx**2))
        zlave = 0.5*(zl1 + zl2)

        z0r = z0 + 0.5*thickness
        aa = (2.*z0r*za - z0r**2 + xx**2)/(za+sqrt((za-z0r)**2-xx**2))
        zr1 = (2.*aa*zb - aa**2 + yy**2)/(zb+sqrt((zb-aa)**2-yy**2))
        bb = (2.*z0r*zb - z0r**2 + yy**2)/(zb+sqrt((zb-z0r)**2-yy**2))
        zr2 = (2.*bb*za - bb**2 + xx**2)/(za+sqrt((za-bb)**2-xx**2))
        zrave = 0.5*(zr1 + zr2)

c       --- Point far to the left of the plate
        if (zz < z0l) then
          delpz(i) = zlave - zz
          cycle
        endif

c       xl1 = ((za-z0l)**2 - (sqrt(yy**2+(zb-zz)**2)+za-zb)**2)
c       yl1 = ((zb-za+sqrt((za-z0l)**2-xx**2))**2 - (zb-zz)**2)
c       xl2 = ((za-zb+sqrt((zb-z0l)**2-yy**2))**2 - (za-zz)**2)
c       yl2 = ((zb-z0l)**2 - (sqrt(xx**2+(za-zz)**2)+zb-za)**2)
        bb = (zz**2 - 2.*zz*zb + yy**2)/(zb+sqrt((zb-zz)**2+yy**2))
        xl1sq = z0l**2 - bb**2 - 2.*(z0l + bb)*za
        aa = (2.*z0l*za - z0l**2 + xx**2)/(za+sqrt((za-z0l)**2-xx**2))
        yl1sq = aa**2 - zz**2 - 2.*(aa - zz)*zb
        bb = (2.*z0l*zb - z0l**2 + yy**2)/(zb+sqrt((zb-z0l)**2-yy**2))
        xl2sq = bb**2 - zz**2 - 2.*(bb - zz)*za
        aa = (zz**2 - 2.*zz*za + xx**2)/(za+sqrt((za-zz)**2+xx**2))
        yl2sq = z0l**2 - aa**2 - 2.*(z0l + aa)*zb

c       xr1 = ((za-z0r)**2 - (sqrt(yy**2+(zb-zz)**2)+za-zb)**2)
c       yr1 = ((zb-za+sqrt((za-z0r)**2-xx**2))**2 - (zb-zz)**2)
c       xr2 = ((za-zb+sqrt((zb-z0r)**2-yy**2))**2 - (za-zz)**2)
c       yr2 = ((zb-z0r)**2 - (sqrt(xx**2+(za-zz)**2)+zb-za)**2)
        bb = (zz**2 - 2.*zz*zb + yy**2)/(zb+sqrt((zb-zz)**2+yy**2))
        xr1sq = z0r**2 - bb**2 - 2.*(z0r + bb)*za
        aa = (2.*z0r*za - z0r**2 + xx**2)/(za+sqrt((za-z0r)**2-xx**2))
        yr1sq = aa**2 - zz**2 - 2.*(aa - zz)*zb
        bb = (2.*z0r*zb - z0r**2 + yy**2)/(zb+sqrt((zb-z0r)**2-yy**2))
        xr2sq = bb**2 - zz**2 - 2.*(bb - zz)*za
        aa = (zz**2 - 2.*zz*za + xx**2)/(za+sqrt((za-zz)**2+xx**2))
        yr2sq = z0r**2 - aa**2 - 2.*(z0r + aa)*zb

c       --- Point to the left of plate
        if (zz <= zlave) then
          xlave = 0.5*(xl1sq + xl2sq)
          ylave = 0.5*(yl1sq + yl2sq)
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
          xl1 = sqrt(xl1sq)
          yl1 = sqrt(yl1sq)
          xl2 = sqrt(xl2sq)
          yl2 = sqrt(yl2sq)
          xlave = 0.5*(xl1 + xl2)
          ylave = 0.5*(yl1 + yl2)

          xrave = 0.5*(xr1sq + xr2sq)
          yrave = 0.5*(yr1sq + yr2sq)
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
          xr1 = sqrt(xr1sq)
          yr1 = sqrt(yr1sq)
          xr2 = sqrt(xr2sq)
          yr2 = sqrt(yr2sq)
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
      subroutine BeamletplateConductorD(za,zb,z0,thickness,xcent,ycent,zcent,
     &                                  n,x,y,z,distance)
      real(kind=8):: za,zb,z0,thickness,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from a plate for the beamlet
c source pre-accelerator.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i,ii
      real(kind=8):: xx,yy,zz,aa,bb
      real(kind=8):: x1,y1,z1,z2,zave,dd

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent

c       --- Iterate to find intersection of normal line going through point
c       --- and the surface
        x1 = xx
        y1 = yy
        if (za < 1.e10 .and. zb < 1.e10) then
          do ii=1,5
c           z1 = zb - sqrt((zb-za+sqrt((za-z0)**2-x1**2))**2 - y1**2)
c           z2 = za - sqrt((za-zb+sqrt((zb-z0)**2-y1**2))**2 - x1**2)
            aa = (2.*z0*za - z0**2 + x1**2)/(za+sqrt((za-z0)**2-x1**2))
            z1 = (2.*aa*zb - aa**2 + y1**2)/(zb+sqrt((zb-aa)**2-y1**2))
            bb = (2.*z0*zb - z0**2 + y1**2)/(zb+sqrt((zb-z0)**2-y1**2))
            z2 = (2.*bb*za - bb**2 + x1**2)/(za+sqrt((za-bb)**2-x1**2))
            zave = 0.5*(z1 + z2)
            x1 = xx*(za - zave)/(za-zz)
            y1 = yy*(zb - zave)/(zb-zz)
          enddo
        else
          zave = z0
        endif

        dd = sqrt((xx - x1)**2 + (yy - y1)**2 + (zz - zave)**2)
        distance(i) = dd - thickness/2.

      enddo

      return
      end
c=======================================================================
c=======================================================================
c=======================================================================
      subroutine ZSrfrvOutConductorF(rofzfunc,zmin,zmax,rmax,griddz,
     &                               xcent,ycent,zcent,n,x,y,z,
     &                               delmx,delpx,delmy,delpy,delmz,delpz,
     &                               fuzz)
      character(*):: rofzfunc
      real(kind=8):: zmin,zmax,rmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

      integer(ISZ):: i,iterrf,ii
      real(kind=8):: srfrv_f
      real(kind=8):: xx,yy,zz,rr,rsrf
      real(kind=8):: xisrf,yisrf,ximax,yimax
      real(kind=8):: zzz,zlower,zhigher,rlower,rhigher,rfuzz
      real(kind=8):: ddz,znew,rnew,zi(2),ri(2)

c     --- Initialize srfrv_f data
      call srfrv_finit

      do i=1,n
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rr = sqrt(xx**2 + yy**2)
        zz = z(i) - zcent

        if (zmin <= zz .and. zz <= zmax) then
c         --- Point inside the z-extent
          rsrf = srfrv_f(zz,rofzfunc,0,0)
          xisrf = sqrt(max(0.,rsrf**2 - yy**2))
          yisrf = sqrt(max(0.,rsrf**2 - xx**2))
          ximax = sqrt(max(0.,rmax**2 - yy**2))
          yimax = sqrt(max(0.,rmax**2 - xx**2))
          if (rr <= rsrf) then
            delmx(i) = xisrf + xx
            delpx(i) = xisrf - xx
            delmy(i) = yisrf + yy
            delpy(i) = yisrf - yy
          else if (rr < rmax) then
            if (yy > rsrf .or. -rsrf > yy .or. rsrf == 0.) then
              delmx(i) = - ximax - xx
              delpx(i) = + xx - ximax
            else
              if (xx > 0.) then
                delmx(i) = xisrf - xx
                delpx(i) = xx - ximax
              else
                delmx(i) = - ximax - xx
                delpx(i) = + xx + xisrf
              endif
            endif
            if (xx > rsrf .or. -rsrf > xx .or. rsrf == 0.) then
              delmy(i) = - yimax - yy
              delpy(i) = + yy - yimax
            else
              if (yy > 0.) then
                delmy(i) = yisrf - yy
                delpy(i) = yy - yimax
              else
                delmy(i) = - yimax - yy
                delpy(i) = yy + yisrf
              endif
            endif
          else
            if (abs(yy) <= rmax) then
              if (xx >= 0.) delmx(i) = + xx - ximax
              if (xx <= 0.) delpx(i) = - ximax - xx
            endif
            if (abs(xx) <= rmax) then
              if (yy >= 0.) delmy(i) = + yy - yimax
              if (yy <= 0.) delpy(i) = - yimax - yy
            endif
          endif
        endif

c Now get delmz delpz
c use griddz

c       --- Only find deltas if point within griddz of the z-extent of
c       --- conductor or within the max radius.
        if (zz >= zmin - griddz .and.
     &      zz <= zmax + griddz .and.
     &      rr <= rmax) then

          zzz = min(max(zz,zmin),zmax)
          zlower = max(zz-griddz,zmin)
          zhigher = min(zz+griddz,zmax)

          rsrf = srfrv_f(zzz,rofzfunc,0,0)
          rlower = srfrv_f(zlower,rofzfunc,0,0)
          rhigher = srfrv_f(zhigher,rofzfunc,0,0)
          rfuzz = max(rsrf,rlower,rhigher)*1.e-12

          if (zz >= zmax .and. rr >= rsrf) then
            delmz(i) = zz - zmax
          elseif (zz < zmin) then
            delmz(i) = LARGEPOS
          elseif (rr > max(rlower,rsrf)) then
            delmz(i) = zmin - zz
          elseif (min(rlower,rsrf) <= rr .and. rr <= max(rlower,rsrf)) then
c           --- If grid location is between rsrf and rlower find
c           --- intersection point.
c           --- First, try Regula-Falsi (secant method) iteration to
c           --- find distance in z to the surface. If that doesn't
c           --- converge rapidly, then go to the brute force binary search.
c           --- Empirically, it was found that stopping at 12 iterations
c           --- gives the best performance, with a broad minimum around
c           --- that point.
            zi(1) = zzz
            ri(1) = rsrf - rr
            zi(2) = zlower
            ri(2) = rlower - rr
            iterrf = 0
            if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
            ii = 1
            do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
              ii = 1
              if (abs(ri(2)) > abs(ri(1))) ii = 2
              if (ri(1) == ri(2)) exit
              zi(ii) = min((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmax)
              ri(ii) = srfrv_f(zi(ii),rofzfunc,0,0) - rr
              iterrf = iterrf + 1
            end do
            if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &          zlower <= zi(ii) .and. zi(ii) < zz) then
              delmz(i) = zz - zi(ii)
            else
c             --- A simple binary search is used to find distance in z
c             --- to the surface. This is used since it is more robust
c             --- than other root finding methods. The major problem
c             --- occurs when the radius changes like a step function.
c             --- The secant method (the original coding) converged
c             --- very slowly in that case, leaving errors as large as
c             --- 10% in the delta and slowing down the code.
              ddz = (min(zz,zmax) - zlower)*0.5
              znew = zlower
              do iterrf=1,32
                znew = znew + ddz
                rnew = srfrv_f(znew,rofzfunc,0,0)
                if (rnew > rr) then
                  ddz = -abs(ddz)*0.5
                else if (rnew <= rr) then
                  ddz = +abs(ddz)*0.5
                endif
              enddo
              delmz(i) = zz - znew
            endif
            if (rr >= rsrf) delmz(i) = -delmz(i)
          endif

c         --- Calculate distance to conductor at higher z.

          if (zz <= zmin .and. rr >= rsrf) then
            delpz(i) = zmin - zz
          elseif (zz > zmax) then
            delpz(i) = LARGEPOS
          elseif (rr > max(rhigher,rsrf)) then
            delpz(i) = zz - zmax
c         --- If grid location is between rrr and rhigher.
          elseif (min(rhigher,rsrf) <= rr .and. rr <= max(rhigher,rsrf)) then
c           --- First, try Regula-Falsi (secant method) iteration to
c           --- find distance in z to the surface. If that doesn't
c           --- converge rapidly, then go to the brute force binary search.
c           --- Empirically, it was found that stopping at 12 iterations
c           --- gives the best performance, with a broad minimum around
c           --- that point.
            zi(1) = zzz
            ri(1) = rsrf - rr
            zi(2) = zhigher
            ri(2) = rhigher - rr
            iterrf = 0
            if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
            ii = 1
            do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
              ii = 1
              if (abs(ri(2)) > abs(ri(1))) ii = 2
              if (ri(1) == ri(2)) exit
              zi(ii) = max((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmin)
              ri(ii) = srfrv_f(zi(ii),rofzfunc,0,0) - rr
              iterrf = iterrf + 1
            end do
            if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &          zz < zi(ii) .and. zi(ii) <= zhigher) then
              delpz(i) = zi(ii) - zz
            else
c             --- A simple binary search is used to find distance in z
c             --- to the surface. This is used since it is more robust
c             --- than other root finding methods. The major problem
c             --- occurs when the radius changes like a step function.
c             --- The secant method (the original coding) converged
c             --- very slowly in that case, leaving errors as large as
c             --- 10% in the delta and slowing down the code.
              ddz = (max(zz,zmin) - zhigher)*0.5
              znew = zhigher
              do iterrf=1,32
                znew = znew + ddz
                rnew = srfrv_f(znew,rofzfunc,0,0)
                if (rnew <= rr) then
                  ddz = -abs(ddz)*0.5
                else if (rnew > rr) then
                  ddz = +abs(ddz)*0.5
                endif
              enddo
              delpz(i) = znew - zz
            endif
            if (rr >= rsrf) delpz(i) = -delpz(i)
          endif

        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZSrfrvOutConductorD(rofzfunc,zmin,zmax,rmax,griddz,
     &                               xcent,ycent,zcent,n,x,y,z,distance)
      character(*):: rofzfunc
      real(kind=8):: zmin,zmax,rmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from the outside of
c a surface of revolution.
c The notation is that a negative distance means that the point is inside.
c NOTICE: This only computes distances radially and can only reliably used
c to determine whether a point is inside or not. The true distance calculation
c would have to account for distances longitudinally.
c A simple (but expensive) solution would be to generate a table of radii
c along z, and to find the one closest to each input data point. Since that
c would be an approximation, care would be needed to ensure that for particles
c near the surface, the distance has the correct sign.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr,rsrf,zd,rd
      real(kind=8):: srfrv_f

c     --- Initialize srfrv_f data
      call srfrv_finit

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        rr = sqrt(xx**2 + yy**2)
        zz = z(i) - zcent

        if (zmin <= zz .and. zz <= zmax) then
          rsrf = srfrv_f(zz,rofzfunc,0,0)
          if (rr <= rsrf) then 
            distance(i) = rsrf - rr
          elseif (rr < rmax) then
            zd = max(zmin-zz,zz-zmax)
            rd = max(rsrf-rr,rr-rmax)
            distance(i) = max(rd,zd)
          else
            distance(i) = rr - rmax
          endif

        elseif (zz < zmin) then
          rsrf = srfrv_f(zmin,rofzfunc,0,0)
          if (rr >= rsrf) then
            distance(i) = zmin - zz
          else
            distance(i) = sqrt((zmin-zz)**2 + (rr-rsrf)**2)
          endif

        else
          rsrf = srfrv_f(zmax,rofzfunc,0,0)
          if (rr >= rsrf) then
            distance(i) = zz - zmax
          else
            distance(i) = sqrt((zz-zmax)**2 + (rr-rsrf)**2)
          endif

        endif

      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine ZSrfrvInConductorF(rofzfunc,zmin,zmax,rmin,griddz,
     &                              xcent,ycent,zcent,n,x,y,z,
     &                              delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      character(*):: rofzfunc
      real(kind=8):: zmin,zmax,rmin,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

      integer(ISZ):: i,iterrf,ii
      real(kind=8):: srfrv_f
      real(kind=8):: xx,yy,zz,rr,rsrf
      real(kind=8):: xisrf,yisrf,ximin,yimin
      real(kind=8):: zzz,zlower,zhigher,rlower,rhigher,rfuzz
      real(kind=8):: ddz,znew,rnew,zi(2),ri(2)

c     --- Initialize srfrv_f data
      call srfrv_finit

      do i=1,n
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rr = sqrt(xx**2 + yy**2)
        zz = z(i) - zcent

        if (zmin <= zz .and. zz <= zmax) then
c         --- Point inside the z-extent
          rsrf = srfrv_f(zz,rofzfunc,0,0)
          xisrf = sqrt(max(0.,rsrf**2 - yy**2))
          yisrf = sqrt(max(0.,rsrf**2 - xx**2))
          ximin = sqrt(max(0.,rmin**2 - yy**2))
          yimin = sqrt(max(0.,rmin**2 - xx**2))
          if (rr <= rmin) then
            delmx(i) = ximin + xx
            delpx(i) = ximin - xx
            delmy(i) = yimin + yy
            delpy(i) = yimin - yy
          else if (rr < rsrf) then
            if (yy > rmin .or. -rmin > yy .or. rmin == 0.) then
              delmx(i) = - xisrf - xx
              delpx(i) = + xx - xisrf
            else
              if (xx > 0.) then
                delmx(i) = ximin - xx
                delpx(i) = xx - xisrf
              else
                delmx(i) = - xisrf - xx
                delpx(i) = + xx + ximin
              endif
            endif
            if (xx > rmin .or. -rmin > xx .or. rmin == 0.) then
              delmy(i) = - yisrf - yy
              delpy(i) = + yy - yisrf
            else
              if (yy > 0.) then
                delmy(i) = yimin - yy
                delpy(i) = yy - yisrf
              else
                delmy(i) = - yisrf - yy
                delpy(i) = yy + yimin
              endif
            endif
          else
            if (abs(yy) <= rsrf) then
              if (xx >= 0.) delmx(i) = + xx - xisrf
              if (xx <= 0.) delpx(i) = - xisrf - xx
            endif
            if (abs(xx) <= rsrf) then
              if (yy >= 0.) delmy(i) = + yy - yisrf
              if (yy <= 0.) delpy(i) = - yisrf - yy
            endif
          endif
        endif

c Now get delmz delpz
c use griddz

c       --- Only find deltas if point within griddz of the z-extent of
c       --- conductor or within the max radius.
        if (zz >= zmin - griddz .and.
     &      zz <= zmax + griddz .and.
     &      rr >= rmin) then

          zzz = min(max(zz,zmin),zmax)
          zlower = max(zz-griddz,zmin)
          zhigher = min(zz+griddz,zmax)

          rsrf = srfrv_f(zzz,rofzfunc,0,0)
          rlower = srfrv_f(zlower,rofzfunc,0,0)
          rhigher = srfrv_f(zhigher,rofzfunc,0,0)
          rfuzz = max(rsrf,rlower,rhigher)*1.e-12

          if (zz >= zmax .and. rr <= rsrf) then
            delmz(i) = zz - zmax
          elseif (zz < zmin) then
            delmz(i) = LARGEPOS
          elseif (rr < min(rlower,rsrf)) then
            delmz(i) = zmin - zz
          elseif (min(rlower,rsrf) <= rr .and. rr <= max(rlower,rsrf)) then
c           --- If grid location is between rsrf and rlower find
c           --- intersection point.
c           --- First, try Regula-Falsi (secant method) iteration to
c           --- find distance in z to the surface. If that doesn't
c           --- converge rapidly, then go to the brute force binary search.
c           --- Empirically, it was found that stopping at 12 iterations
c           --- gives the best performance, with a broad minimum around
c           --- that point.
            zi(1) = zzz
            ri(1) = rsrf - rr
            zi(2) = zlower
            ri(2) = rlower - rr
            iterrf = 0
            if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
            ii = 1
            do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
              ii = 1
              if (abs(ri(2)) > abs(ri(1))) ii = 2
              if (ri(1) == ri(2)) exit
              zi(ii) = min((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmax)
              ri(ii) = srfrv_f(zi(ii),rofzfunc,0,0) - rr
              iterrf = iterrf + 1
            end do
            if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &          zlower <= zi(ii) .and. zi(ii) < zz) then
              delmz(i) = zz - zi(ii)
            else
c             --- A simple binary search is used to find distance in z
c             --- to the surface. This is used since it is more robust
c             --- than other root finding methods. The major problem
c             --- occurs when the radius changes like a step function.
c             --- The secant method (the original coding) converged
c             --- very slowly in that case, leaving errors as large as
c             --- 10% in the delta and slowing down the code.
              ddz = (min(zz,zmax) - zlower)*0.5
              znew = zlower
              do iterrf=1,32
                znew = znew + ddz
                rnew = srfrv_f(znew,rofzfunc,0,0)
                if (rnew > rr) then
                  ddz = -abs(ddz)*0.5
                else if (rnew <= rr) then
                  ddz = +abs(ddz)*0.5
                endif
              enddo
              delmz(i) = zz - znew
            endif
            if (rr <= rsrf) delmz(i) = -delmz(i)
          endif

c         --- Calculate distance to conductor at higher z.

          if (zz <= zmin .and. rr <= rsrf) then
            delpz(i) = zmin - zz
          elseif (zz > zmax) then
            delpz(i) = LARGEPOS
          elseif (rr < min(rhigher,rsrf)) then
            delpz(i) = zz - zmax
c         --- If grid location is between rrr and rhigher.
          elseif (min(rhigher,rsrf) <= rr .and. rr <= max(rhigher,rsrf)) then
c           --- First, try Regula-Falsi (secant method) iteration to
c           --- find distance in z to the surface. If that doesn't
c           --- converge rapidly, then go to the brute force binary search.
c           --- Empirically, it was found that stopping at 12 iterations
c           --- gives the best performance, with a broad minimum around
c           --- that point.
            zi(1) = zzz
            ri(1) = rsrf - rr
            zi(2) = zhigher
            ri(2) = rhigher - rr
            iterrf = 0
            if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
            ii = 1
            do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
              ii = 1
              if (abs(ri(2)) > abs(ri(1))) ii = 2
              if (ri(1) == ri(2)) exit
              zi(ii) = max((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmin)
              ri(ii) = srfrv_f(zi(ii),rofzfunc,0,0) - rr
              iterrf = iterrf + 1
            end do
            if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &          zz < zi(ii) .and. zi(ii) <= zhigher) then
              delpz(i) = zi(ii) - zz
            else
c             --- A simple binary search is used to find distance in z
c             --- to the surface. This is used since it is more robust
c             --- than other root finding methods. The major problem
c             --- occurs when the radius changes like a step function.
c             --- The secant method (the original coding) converged
c             --- very slowly in that case, leaving errors as large as
c             --- 10% in the delta and slowing down the code.
              ddz = (max(zz,zmin) - zhigher)*0.5
              znew = zhigher
              do iterrf=1,32
                znew = znew + ddz
                rnew = srfrv_f(znew,rofzfunc,0,0)
                if (rnew <= rr) then
                  ddz = -abs(ddz)*0.5
                else if (rnew > rr) then
                  ddz = +abs(ddz)*0.5
                endif
              enddo
              delpz(i) = znew - zz
            endif
            if (rr <= rsrf) delpz(i) = -delpz(i)
          endif

        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZSrfrvInConductorD(rofzfunc,zmin,zmax,rmin,griddz,
     &                              xcent,ycent,zcent,n,x,y,z,distance)
      character(*):: rofzfunc
      real(kind=8):: zmin,zmax,rmin,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from the inside of
c a surface of revolution.
c The notation is that a negative distance means that the point is inside.
c NOTICE: This only computes distances radially and can only reliably used
c to determine whether a point is inside or not. The true distance calculation
c would have to account for distances longitudinally.
c A simple (but expensive) solution would be to generate a table of radii
c along z, and to find the one closest to each input data point. Since that
c would be an approximation, care would be needed to ensure that for particles
c near the surface, the distance has the correct sign.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr,rsrf,zd,rd
      real(kind=8):: srfrv_f

c     --- Initialize srfrv_f data
      call srfrv_finit

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        if (zmin <= zz .and. zz <= zmax) then
          rsrf = srfrv_f(zz,rofzfunc,0,0)
          if (rr >= rsrf) then 
            distance(i) = rr - rsrf
          elseif (rr > rmin) then
            zd = max(zmin-zz,zz-zmax)
            rd = max(rr-rsrf,rmin-rr)
            distance(i) = max(rd,zd)
          else
            distance(i) = rmin - rr
          endif

        elseif (zz < zmin) then
          rsrf = srfrv_f(zmin,rofzfunc,0,0)
          if (rr <= rsrf) then
            distance(i) = zmin - zz
          else
            distance(i) = sqrt((zmin-zz)**2 + (rr-rsrf)**2)
          endif

        else
          rsrf = srfrv_f(zmax,rofzfunc,0,0)
          if (rr <= rsrf) then
            distance(i) = zz - zmax
          else
            distance(i) = sqrt((zz-zmax)**2 + (rr-rsrf)**2)
          endif

        endif

      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine ZSrfrvInOutConductorF(rminofz,rmaxofz,zmin,zmax,griddz,
     &                                 xcent,ycent,zcent,n,x,y,z,
     &                                 delmx,delpx,delmy,delpy,delmz,delpz,fuzz)
      character(*):: rminofz,rmaxofz
      real(kind=8):: zmin,zmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

      integer(ISZ):: i,iterrf,ii
      real(kind=8):: srfrv_f
      real(kind=8):: xx,yy,zz,rr,rmin,rmax
      real(kind=8):: ximin,yimin,ximax,yimax
      real(kind=8):: zzz,zlower,zhigher,rfuzz
      real(kind=8):: rminlower,rminhigher,rmaxlower,rmaxhigher
      real(kind=8):: ddz,znew,rnew,zi(2),ri(2)
      real(kind=8):: delmz1,delpz1

c     --- Initialize srfrv_f data
      call srfrv_finit

      do i=1,n
        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

c       --- Set temps
        xx = x(i) - xcent
        yy = y(i) - ycent
        rr = sqrt(xx**2 + yy**2)
        zz = z(i) - zcent

        if (zmin <= zz .and. zz <= zmax) then
c         --- Point inside the z-extent
          rmin = srfrv_f(zz,rminofz,1,0)
          rmax = srfrv_f(zz,rmaxofz,2,0)
          ximin = sqrt(max(0.,rmin**2 - yy**2))
          yimin = sqrt(max(0.,rmin**2 - xx**2))
          ximax = sqrt(max(0.,rmax**2 - yy**2))
          yimax = sqrt(max(0.,rmax**2 - xx**2))
          if (rr < rmin) then
            delmx(i) = ximin + xx
            delpx(i) = ximin - xx
            delmy(i) = yimin + yy
            delpy(i) = yimin - yy
          else if (rr < rmax) then
            if (yy > rmin .or. -rmin > yy .or. rmin == 0.) then
              delmx(i) = - ximax - xx
              delpx(i) = + xx - ximax
            else
              if (xx > 0.) then
                delmx(i) = ximin - xx
                delpx(i) = xx - ximax
              else
                delmx(i) = - ximax - xx
                delpx(i) = + xx + ximin
              endif
            endif
            if (xx > rmin .or. -rmin > xx .or. rmin == 0.) then
              delmy(i) = - yimax - yy
              delpy(i) = + yy - yimax
            else
              if (yy > 0.) then
                delmy(i) = yimin - yy
                delpy(i) = yy - yimax
              else
                delmy(i) = - yimax - yy
                delpy(i) = yy + yimin
              endif
            endif
          else
            if (abs(yy) <= rmax) then
              if (xx >= 0.) delmx(i) = + xx - ximax
              if (xx <= 0.) delpx(i) = - ximax - xx
            endif
            if (abs(xx) <= rmax) then
              if (yy >= 0.) delmy(i) = + yy - yimax
              if (yy <= 0.) delpy(i) = - yimax - yy
            endif
          endif
        endif

c Now get delmz delpz
c use griddz

c       --- Only find deltas if point within griddz of the z-extent of
c       --- conductor or within the max radius.
        if (zz >= zmin - griddz .and.
     &      zz <= zmax + griddz) then

          zzz = min(max(zz,zmin),zmax)
          zlower = max(zz-griddz,zmin)
          zhigher = min(zz+griddz,zmax)

          rmin = srfrv_f(zzz,rminofz,1,0)
          rmax = srfrv_f(zzz,rmaxofz,2,0)
          rminlower = srfrv_f(zlower,rminofz,1,-1)
          rminhigher = srfrv_f(zhigher,rminofz,1,+1)
          rmaxlower = srfrv_f(zlower,rmaxofz,2,-1)
          rmaxhigher = srfrv_f(zhigher,rmaxofz,2,+1)
          rfuzz = max(rmin,rmax,rminlower,rminhigher,rmaxlower,rmaxhigher)*1.e-12

          if (zz < zmin .or. rr < min(rminlower,rmin) .or.
     &                       max(rmaxlower,rmax) < rr) then
            delmz(i) = LARGEPOS
          elseif (zz >= zmax .and. rmin <= rr .and. rr <= rmax) then
            delmz(i) = zz - zmax
          elseif (max(rminlower,rmin) < rr .and. rr < min(rmaxlower,rmax)) then
            delmz(i) = zmin - zz
          else
            delmz(i) = LARGEPOS
            if (min(rminlower,rmin) <= rr .and. rr <= max(rminlower,rmin)) then
c             --- If grid location is between rmin and rminlower find
c             --- intersection point.
c             --- First, try Regula-Falsi (secant method) iteration to
c             --- find distance in z to the surface. If that doesn't
c             --- converge rapidly, then go to the brute force binary search.
c             --- Empirically, it was found that stopping at 12 iterations
c             --- gives the best performance, with a broad minimum around
c             --- that point.
              zi(1) = zzz
              ri(1) = rmin - rr
              zi(2) = zlower
              ri(2) = rminlower - rr
              iterrf = 0
              if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
              ii = 1
              do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
                ii = 1
                if (abs(ri(2)) > abs(ri(1))) ii = 2
                if (ri(1) == ri(2)) exit
                zi(ii) = min((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmax)
                ri(ii) = srfrv_f(zi(ii),rminofz,1,0) - rr
                iterrf = iterrf + 1
              end do
              if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &            zlower <= zi(ii) .and. zi(ii) < zz) then
                delmz(i) = zz - zi(ii)
              else
c               --- A simple binary search is used to find distance in z
c               --- to the surface. This is used since it is more robust
c               --- than other root finding methods. The major problem
c               --- occurs when the radius changes like a step function.
c               --- The secant method (the original coding) converged
c               --- very slowly in that case, leaving errors as large as
c               --- 10% in the delta and slowing down the code.
                ddz = (min(zz,zmax) - zlower)*0.5
                znew = zlower
                do iterrf=1,32
                  znew = znew + ddz
                  rnew = srfrv_f(znew,rminofz,1,0)
                  if (rnew > rr) then
                    ddz = -abs(ddz)*0.5
                  else if (rnew <= rr) then
                    ddz = +abs(ddz)*0.5
                  endif
                enddo
                delmz(i) = zz - znew
              endif
              if (rr >= rmin) delmz(i) = -delmz(i)
            endif
            if (min(rmaxlower,rmax) <= rr .and. rr <= max(rmaxlower,rmax)) then
c             --- If grid location is between rmax and rmaxlower find
c             --- intersection point.
c             --- First, try Regula-Falsi (secant method) iteration to
c             --- find distance in z to the surface. If that doesn't
c             --- converge rapidly, then go to the brute force binary search.
c             --- Empirically, it was found that stopping at 12 iterations
c             --- gives the best performance, with a broad minimum around
c             --- that point.
              zi(1) = zzz
              ri(1) = rmax - rr
              zi(2) = zlower
              ri(2) = rmaxlower - rr
              iterrf = 0
              if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
              ii = 1
              do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
                ii = 1
                if (abs(ri(2)) > abs(ri(1))) ii = 2
                if (ri(1) == ri(2)) exit
                zi(ii) = min((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmax)
                ri(ii) = srfrv_f(zi(ii),rmaxofz,2,0) - rr
                iterrf = iterrf + 1
              end do
              if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &            zlower <= zi(ii) .and. zi(ii) < zz) then
                delmz1 = zz - zi(ii)
              else
c               --- A simple binary search is used to find distance in z
c               --- to the surface. This is used since it is more robust
c               --- than other root finding methods. The major problem
c               --- occurs when the radius changes like a step function.
c               --- The secant method (the original coding) converged
c               --- very slowly in that case, leaving errors as large as
c               --- 10% in the delta and slowing down the code.
                ddz = (min(zz,zmax) - zlower)*0.5
                znew = zlower
                do iterrf=1,32
                  znew = znew + ddz
                  rnew = srfrv_f(znew,rmaxofz,2,0)
                  if (rnew > rr) then
                    ddz = -abs(ddz)*0.5
                  else if (rnew <= rr) then
                    ddz = +abs(ddz)*0.5
                  endif
                enddo
                delmz1 = zz - znew
              endif
              if (delmz1 < abs(delmz(i))) then
                if (rr <= rmax) then
                  delmz(i) = -delmz1
                else
                  delmz(i) = +delmz1
                endif
              endif
            endif
          endif

c         --- Calculate distance to conductor at higher z.

          if (zz > zmax .or. rr < min(rminhigher,rmin) .or.
     &                       max(rmaxhigher,rmax) < rr) then
            delpz(i) = LARGEPOS
          elseif (zz <= zmin .and. rmin <= rr .and. rr <= rmax) then
            delpz(i) = zmin - zz
          elseif (max(rminhigher,rmin) < rr .and.
     &            rr < min(rmaxhigher,rmax)) then
            delpz(i) = zz - zmax
          else
            delpz(i) = LARGEPOS
            if (min(rminhigher,rmin) <= rr .and.
     &          rr <= max(rminhigher,rmin)) then
c             --- First, try Regula-Falsi (secant method) iteration to
c             --- find distance in z to the surface. If that doesn't
c             --- converge rapidly, then go to the brute force binary search.
c             --- Empirically, it was found that stopping at 12 iterations
c             --- gives the best performance, with a broad minimum around
c             --- that point.
              zi(1) = zzz
              ri(1) = rmin - rr
              zi(2) = zhigher
              ri(2) = rminhigher - rr
              iterrf = 0
              if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
              ii = 1
              do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
                ii = 1
                if (abs(ri(2)) > abs(ri(1))) ii = 2
                if (ri(1) == ri(2)) exit
                zi(ii) = max((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmin)
                ri(ii) = srfrv_f(zi(ii),rminofz,1,0) - rr
                iterrf = iterrf + 1
              end do
              if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &            zz < zi(ii) .and. zi(ii) <= zhigher) then
                delpz(i) = zi(ii) - zz
              else
c               --- A simple binary search is used to find distance in z
c               --- to the surface. This is used since it is more robust
c               --- than other root finding methods. The major problem
c               --- occurs when the radius changes like a step function.
c               --- The secant method (the original coding) converged
c               --- very slowly in that case, leaving errors as large as
c               --- 10% in the delta and slowing down the code.
                ddz = (max(zz,zmin) - zhigher)*0.5
                znew = zhigher
                do iterrf=1,32
                  znew = znew + ddz
                  rnew = srfrv_f(znew,rminofz,1,0)
                  if (rnew <= rr) then
                    ddz = -abs(ddz)*0.5
                  else if (rnew > rr) then
                    ddz = +abs(ddz)*0.5
                  endif
                enddo
                delpz(i) = znew - zz
              endif
              if (rr >= rmin) delpz(i) = -delpz(i)
            endif
            if (min(rmaxhigher,rmax) <= rr .and.
     &          rr <= max(rmaxhigher,rmax)) then
c             --- First, try Regula-Falsi (secant method) iteration to
c             --- find distance in z to the surface. If that doesn't
c             --- converge rapidly, then go to the brute force binary search.
c             --- Empirically, it was found that stopping at 12 iterations
c             --- gives the best performance, with a broad minimum around
c             --- that point.
              zi(1) = zzz
              ri(1) = rmax - rr
              zi(2) = zhigher
              ri(2) = rmaxhigher - rr
              iterrf = 0
              if (ri(1) == 0. .or. ri(2) == 0.) iterrf = 12
              ii = 1
              do while (abs(ri(2))+abs(ri(1)) > rfuzz .and. iterrf < 12)
                ii = 1
                if (abs(ri(2)) > abs(ri(1))) ii = 2
                if (ri(1) == ri(2)) exit
                zi(ii) = max((zi(1)*ri(2)-zi(2)*ri(1))/(ri(2)-ri(1)),zmin)
                ri(ii) = srfrv_f(zi(ii),rmaxofz,2,0) - rr
                iterrf = iterrf + 1
              end do
              if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &            zz < zi(ii) .and. zi(ii) <= zhigher) then
                delpz1 = zi(ii) - zz
              else
c               --- A simple binary search is used to find distance in z
c               --- to the surface. This is used since it is more robust
c               --- than other root finding methods. The major problem
c               --- occurs when the radius changes like a step function.
c               --- The secant method (the original coding) converged
c               --- very slowly in that case, leaving errors as large as
c               --- 10% in the delta and slowing down the code.
                ddz = (max(zz,zmin) - zhigher)*0.5
                znew = zhigher
                do iterrf=1,32
                  znew = znew + ddz
                  rnew = srfrv_f(znew,rmaxofz,2,0)
                  if (rnew <= rr) then
                    ddz = -abs(ddz)*0.5
                  else if (rnew > rr) then
                    ddz = +abs(ddz)*0.5
                  endif
                enddo
                delpz1 = znew - zz
              endif
              if (delpz1 < abs(delpz(i))) then
                if (rr <= rmax) then
                  delpz(i) = -delpz1
                else
                  delpz(i) = +delpz1
                endif
              endif
            endif
          endif

        endif

      enddo
      return
      end
c=======================================================================
      subroutine ZSrfrvInOutConductorD(rminofz,rmaxofz,zmin,zmax,griddz,
     &                                 xcent,ycent,zcent,n,x,y,z,distance)
      character(*):: rminofz,rmaxofz
      real(kind=8):: zmin,zmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from between
c surfaces of revolution.
c The notation is that a negative distance means that the point is inside.
c NOTICE: This only computes distances radially and can only reliably used
c to determine whether a point is inside or not. The true distance calculation
c would have to account for distances longitudinally.
c A simple (but expensive) solution would be to generate a table of radii
c along z, and to find the one closest to each input data point. Since that
c would be an approximation, care would be needed to ensure that for particles
c near the surface, the distance has the correct sign.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr,rmin,rmax,zd,rd
      real(kind=8):: srfrv_f

c     --- Initialize srfrv_f data
      call srfrv_finit

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        if (zmin <= zz .and. zz <= zmax) then
          rmin = srfrv_f(zz,rminofz,1,0)
          rmax = srfrv_f(zz,rmaxofz,2,0)
          if (rr >= rmax) then 
            distance(i) = rr - rmax
          elseif (rr > rmin) then
            zd = max(zmin-zz,zz-zmax)
            rd = max(rr-rmax,rmin-rr)
            distance(i) = max(rd,zd)
          else
            distance(i) = rmin - rr
          endif

        elseif (zz < zmin) then
          rmin = srfrv_f(zmin,rminofz,1,0)
          rmax = srfrv_f(zmin,rmaxofz,2,0)
          if (rr < rmin) then
            distance(i) = sqrt((zmin-zz)**2 + (rr-rmin)**2)
          else if (rr <= rmax) then
            distance(i) = zmin - zz
          else
            distance(i) = sqrt((zmin-zz)**2 + (rr-rmax)**2)
          endif

        else
          rmin = srfrv_f(zmax,rminofz,1,0)
          rmax = srfrv_f(zmax,rmaxofz,2,0)
          if (rr < rmin) then
            distance(i) = sqrt((zz-zmax)**2 + (rr-rmin)**2)
          else if (rr <= rmax) then
            distance(i) = zz - zmax
          else
            distance(i) = sqrt((zz-zmax)**2 + (rr-rmax)**2)
          endif

        endif
      enddo

      return
      end
c=======================================================================
