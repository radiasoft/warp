#include "top.h"
c=======================================================================
c@(#) File F3D_CONDUCTORS.F, $Revision: 1.4 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  written by David P. Grote
c
c  It contains the routines for calculating the intersection of points
c  with various geometrical objects.
c=======================================================================
c=======================================================================
      subroutine setconductorparity(nn,ix,iy,iz,dels,parity,fuzz)
      integer(ISZ):: nn
      integer(ISZ):: ix(nn),iy(nn),iz(nn),parity(nn)
      real(kind=8):: dels(0:5,nn)
      real(kind=8):: fuzz

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
      
      do i=1,nn
        parity(i) = 999

        iparity = mod(ix(i)+iy(i)+iz(i),2)
        if (dels(0,i) < 1.-fuzz) parity(i) = iparity
        if (dels(1,i) < 1.-fuzz) parity(i) = iparity
        if (dels(2,i) < 1.-fuzz) parity(i) = iparity
        if (dels(3,i) < 1.-fuzz) parity(i) = iparity
        if (dels(4,i) < 1.-fuzz) parity(i) = iparity
        if (dels(5,i) < 1.-fuzz) parity(i) = iparity

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

      enddo

      return
      end
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
          delmx(i) = -(sqrt(rad**2 - yy*yy) + xx)
          delpx(i) = -(sqrt(rad**2 - yy*yy) - xx)
          delmy(i) = -(sqrt(rad**2 - xx*xx) + yy)
          delpy(i) = -(sqrt(rad**2 - xx*xx) - yy)
          delmz(i) = zmin - z(i)
          delpz(i) = z(i) - zmax
          cycle
        endif

c       --- Point is outside the cylinder
        if (rsq > rad**2) then
          if (abs(yy) <= rad) then
            if (xx > 0) delmx(i) = ( xx - sqrt(rad**2 - yy*yy))
            if (xx < 0) delpx(i) = (-xx - sqrt(rad**2 - yy*yy))
          endif
          if (abs(xx) <= rad) then
            if (yy > 0) delmy(i) = ( yy - sqrt(rad**2 - xx*xx))
            if (yy < 0) delpy(i) = (-yy - sqrt(rad**2 - xx*xx))
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
          delmx(i) = (sqrt(rad**2 - yy*yy) + xx)
          delpx(i) = (sqrt(rad**2 - yy*yy) - xx)
          delmy(i) = (sqrt(rad**2 - xx*xx) + yy)
          delpy(i) = (sqrt(rad**2 - xx*xx) - yy)
          cycle
        endif

c       --- Point is outside the cylinder
        if (rsq > rad**2) then
          if (abs(yy) <= rad) then
            if (xx > 0) delmx(i) = -( xx - sqrt(rad**2 - yy*yy))
            if (xx < 0) delpx(i) = -(-xx - sqrt(rad**2 - yy*yy))
          endif
          if (abs(xx) <= rad) then
            if (yy > 0) delmy(i) = -( yy - sqrt(rad**2 - xx*xx))
            if (yy < 0) delpy(i) = -(-yy - sqrt(rad**2 - xx*xx))
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
          delmx(i) = -(sqrt(rad**2 - yy*yy - zz*zz) + xx)
          delpx(i) = -(sqrt(rad**2 - yy*yy - zz*zz) - xx)
          delmy(i) = -(sqrt(rad**2 - xx*xx - zz*zz) + yy)
          delpy(i) = -(sqrt(rad**2 - xx*xx - zz*zz) - yy)
          delmz(i) = -(sqrt(rad**2 - xx*xx - yy*yy) + zz)
          delpz(i) = -(sqrt(rad**2 - xx*xx - yy*yy) - zz)
        endif

c       --- Point outside the sphere
        if (rsq > rad**2) then
          if ((yy*yy + zz*zz) <= rad**2) then
            if (xx > 0) delmx(i) = ( xx - sqrt(rad**2 - yy*yy - zz*zz))
            if (xx < 0) delpx(i) = (-xx - sqrt(rad**2 - yy*yy - zz*zz))
          endif
          if ((xx*xx + zz*zz) <= rad**2) then
            if (yy > 0) delmy(i) = ( yy - sqrt(rad**2 - xx*xx - zz*zz))
            if (yy < 0) delpy(i) = (-yy - sqrt(rad**2 - xx*xx - zz*zz))
          endif
          if ((xx*xx + yy*yy) <= rad**2) then
            if (zz > 0) delmz(i) = ( zz - sqrt(rad**2 - xx*xx - yy*yy))
            if (zz < 0) delpz(i) = (-zz - sqrt(rad**2 - xx*xx - yy*yy))
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
          delmx(i) = -(sqrt(rcone**2 - yy*yy) + xx)
          delpx(i) = -(sqrt(rcone**2 - yy*yy) - xx)
          delmy(i) = -(sqrt(rcone**2 - xx*xx) + yy)
          delpy(i) = -(sqrt(rcone**2 - xx*xx) - yy)

c       --- Point is outside the cone
        else if (r > rcone) then
          if (abs(yy) <= rcone) then
            if (xx > 0) delmx(i) = ( xx - sqrt(rcone**2 - yy*yy))
            if (xx < 0) delpx(i) = (-xx - sqrt(rcone**2 - yy*yy))
          endif
          if (abs(xx) <= rcone) then
            if (yy > 0) delmy(i) = ( yy - sqrt(rcone**2 - xx*xx))
            if (yy < 0) delpy(i) = (-yy - sqrt(rcone**2 - xx*xx))
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
          delmx(i) = (sqrt(rcone**2 - yy*yy) + xx)
          delpx(i) = (sqrt(rcone**2 - yy*yy) - xx)
          delmy(i) = (sqrt(rcone**2 - xx*xx) + yy)
          delpy(i) = (sqrt(rcone**2 - xx*xx) - yy)

c       --- Point is outside the cone
        else if (r > rcone) then
          if (abs(yy) <= rcone) then
            if (xx > 0) delmx(i) = -( xx - sqrt(rcone**2 - yy*yy))
            if (xx < 0) delpx(i) = -(-xx - sqrt(rcone**2 - yy*yy))
          endif
          if (abs(xx) <= rcone) then
            if (yy > 0) delmy(i) = -( yy - sqrt(rcone**2 - xx*xx))
            if (yy < 0) delpy(i) = -(-yy - sqrt(rcone**2 - xx*xx))
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
          zl = zcent - sqrt(r2**2 - (r - r1)**2)
          zr = zcent + sqrt(r2**2 - (r - r1)**2)
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
          ri = r1 - sqrt(r2**2 - (z(i) - zcent)**2)
          ro = r1 + sqrt(r2**2 - (z(i) - zcent)**2)
          if (z(i) == zcent) then
            ri = r1 - r2
            ro = r1 + r2
          else if (z(i) == zcent - r2 .or. z(i) == zcent + r2) then
            ri = r1
            ro = r1
          endif
          if (r < ri) then
            delmx(i) = sqrt(ri**2 - yy*yy) + xx
            delpx(i) = sqrt(ri**2 - yy*yy) - xx
            delmy(i) = sqrt(ri**2 - xx*xx) + yy
            delpy(i) = sqrt(ri**2 - xx*xx) - yy
          else if (ri == r .or. r == ro) then
            delmx(i) = 0.
            delpx(i) = 0.
            delmy(i) = 0.
            delpy(i) = 0.
          else if (ri < r .and. r < ro) then
            delmx(i) = -(sqrt(ro**2 - yy*yy) + xx)
            delpx(i) = -(sqrt(ro**2 - yy*yy) - xx)
            delmy(i) = -(sqrt(ro**2 - xx*xx) + yy)
            delpy(i) = -(sqrt(ro**2 - xx*xx) - yy)
            if (abs(yy) < ri) then
              if (xx <= 0.) delpx(i) = xx + sqrt(ri**2 - yy*yy)
              if (xx >= 0.) delmx(i) = sqrt(ri**2 - yy*yy) - xx
            else if (abs(yy) == ri) then
              if (xx <= 0.) delpx(i) = xx
              if (xx >= 0.) delmx(i) = -xx
            endif
            if (abs(xx) < ri) then
              if (yy <= 0.) delpy(i) = yy + sqrt(ri**2 - xx*xx)
              if (yy >= 0.) delmy(i) = sqrt(ri**2 - xx*xx) - yy
            else if (abs(xx) == ri) then
              if (yy <= 0.) delpy(i) = yy
              if (yy >= 0.) delmy(i) = -yy
            endif
          else if (r > ro) then
            if (abs(yy) < ro) then
              if (xx <= 0.) delpx(i) = -sqrt(ro**2 - yy*yy) - xx
              if (xx >= 0.) delmx(i) = xx - sqrt(ro**2 - yy*yy)
            endif
            if (abs(xx) < ro) then
              if (yy <= 0.) delpy(i) = -sqrt(ro**2 - xx*xx) - yy
              if (yy >= 0.) delmy(i) = yy - sqrt(ro**2 - xx*xx)
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
            ss = sqrt(b**2 - 4*a*c)
            z1 = (-b-ss)/(2*a)
            z2 = (-b+ss)/(2*a)
          endif
          zl = min(z1,z2)
          zu = max(z1,z2)
c         --- Get the point where the z-line intersects the end planes
          z1 = (zcent-lenhalf - xx*stheta*cphi - yy*sphi)*cci
          z2 = (zcent+lenhalf - xx*stheta*cphi - yy*sphi)*cci
          zlend = min(z1,z2)
          zuend = max(z1,z2)
c         --- Now, get distances to surface
          if (z(i) < zl) then
            if (zl <= zlend .and. zlend <= zu) delpz(i) = zlend - z(i)
            if (zlend <= zl .and. zl <= zuend) delpz(i) = zl - z(i)
          else if (zl <= z(i) .and. z(i) <= zu) then
            if (zl <= zlend .and. zlend <= zu .and. z(i) <= zuend) then
              if (z(i) <= zlend) then
                delpz(i) = zlend - z(i)
              else
                delpz(i) = z(i) - min(zu,zuend)
                delmz(i) = zlend - z(i)
              endif
            else if (zlend <= zl .and. zu <= zuend) then
              delmz(i) = zl - z(i)
              delpz(i) = z(i) - zu
            else if (zl <= zuend .and. zuend <= zu) then
              if (z(i) <= zuend) then
                delpz(i) = z(i) - zuend
                delmz(i) = max(zl,zlend) - z(i)
              else
                delmz(i) = z(i) - zuend
              endif
            endif
          else if (zu < z(i)) then
            if (zu >= zuend .and. zuend >= zl) delmz(i) = z(i) - zuend
            if (zuend >= zu .and. zu >= zlend) delmz(i) = z(i) - zu
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
            ss = sqrt(b**2 - 4*a*c)
            x1 = (-b-ss)/(2*a)
            x2 = (-b+ss)/(2*a)
          endif
          xl = min(x1,x2)
          xu = max(x1,x2)
c         --- Get the point where the x-line intersects the end planes
          x1 = (zcent-lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          x2 = (zcent+lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          xlend = min(x1,x2)
          xuend = max(x1,x2)
c         --- Now, get distances to surface
          if (x(i) < xl) then
            if (xl <= xlend .and. xlend <= xu) delpx(i) = xlend - x(i)
            if (xlend <= xl .and. xl <= xuend) delpx(i) = xl - x(i)
          else if (xl <= x(i) .and. x(i) <= xu) then
            if (xl <= xlend .and. xlend <= xu .and. x(i) <= xuend) then
              if (x(i) <= xlend) then
                delpx(i) = xlend - x(i)
              else
                delpx(i) = x(i) - min(xu,xuend)
                delmx(i) = xlend - x(i)
              endif
            else if (xlend <= xl .and. xu <= xuend) then
              delmx(i) = xl - x(i)
              delpx(i) = x(i) - xu
            else if (xl <= xuend .and. xuend <= xu) then
              if (x(i) <= xuend) then
                delpx(i) = x(i) - xuend
                delmx(i) = max(xl,xlend) - x(i)
              else
                delmx(i) = x(i) - xuend
              endif
            endif
          else if (xu < x(i)) then
            if (xu >= xuend .and. xuend >= xl) delmx(i) = x(i) - xuend
            if (xuend >= xu .and. xu >= xlend) delmx(i) = x(i) - xu
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
            ss = sqrt(b**2 - 4*a*c)
            y1 = (-b-ss)/(2*a)
            y2 = (-b+ss)/(2*a)
          endif
          yl = min(y1,y2)
          yu = max(y1,y2)
c         --- Get the point where the z-line intersects the end planes
          y1 = (zcent-lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          y2 = (zcent+lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          ylend = min(y1,y2)
          yuend = max(y1,y2)
c         --- Now, get distances to surface
          if (y(i) < yl) then
            if (yl <= ylend .and. ylend <= yu) delpy(i) = ylend - y(i)
            if (ylend <= yl .and. yl <= yuend) delpy(i) = yl - y(i)
          else if (yl <= y(i) .and. y(i) <= yu) then
            if (yl <= ylend .and. ylend <= yu .and. y(i) <= yuend) then
              if (y(i) <= ylend) then
                delpy(i) = ylend - y(i)
              else
                delpy(i) = y(i) - min(yu,yuend)
                delmy(i) = ylend - y(i)
              endif
            else if (ylend <= yl .and. yu <= yuend) then
              delmy(i) = yl - y(i)
              delpy(i) = y(i) - yu
            else if (yl <= yuend .and. yuend <= yu) then
              if (y(i) <= yuend) then
                delpy(i) = y(i) - yuend
                delmy(i) = max(yl,ylend) - y(i)
              else
                delmy(i) = y(i) - yuend
              endif
            endif
          else if (yu < y(i)) then
            if (yu >= yuend .and. yuend >= yl) delmy(i) = y(i) - yuend
            if (yuend >= yu .and. yu >= ylend) delmy(i) = y(i) - yu
          endif
        endif

      enddo

      return
      end
c=======================================================================
