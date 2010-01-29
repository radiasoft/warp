#include "top.h"
c=======================================================================
c@(#) File F3D_CONDUCTORS.F, $Revision: 1.89 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  written by David P. Grote
c
c  It contains the routines for calculating the intersection of points
c  with various geometrical objects.
c=======================================================================
c=======================================================================
      subroutine conductordelfromintercepts(intercepts,conductors,
     &                                      dfill,fuzz)
      use ConductorInterceptTypeModule
      use ConductorTypeModule
      type(ConductorInterceptType):: intercepts
      type(ConductorType):: conductors
      real(kind=8):: dfill,fuzz

c Given information about the intercepts of grid lines with an object,
c contained in the intercepts instance, generate the subgrid information
c needed by the field solve, putting the data directly into the conductors
c instance.
c The data in intercepts is a list of interception points along the grid
c lines in each of the three dimensions. The interceptions need to be in
c increasing order and the first point must be at the lower edge. Note that
c if the lower edge extends to -infinity, the value of -LARGEPOS should be
c used. If there are no intercepts, the list should be filled with LARGEPOS.
c This routine assumes that the conductors object is empty. If there is
c any data in it, it is freed.

      integer(ISZ):: iouter

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: ixlocal,iylocal,izlocal
      integer(ISZ):: nxicpt,nyicpt,nzicpt
      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)
      real(kind=8),pointer:: xvoltages(:,:,:)
      real(kind=8),pointer:: yvoltages(:,:,:)
      real(kind=8),pointer:: zvoltages(:,:,:)
      integer(ISZ),pointer:: xcondids(:,:,:)
      integer(ISZ),pointer:: ycondids(:,:,:)
      integer(ISZ),pointer:: zcondids(:,:,:)

      integer(ISZ):: nc,ne,no
      integer(ISZ):: ncmax,nemax,nomax
      integer(ISZ):: ix,iy,iz,ii,idel
      integer(ISZ):: ix1,ix2,iy1,iy2,iz1,iz2
      real(kind=8):: dxi,dyi,dzi
      real(kind=8):: x1,x2,y1,y2,z1,z2
      real(kind=8):: xx,yy,zz
      real(kind=8):: delx,dely,delz
      real(kind=8):: xside,yside,zside
      integer(ISZ):: indx(0:2)
      real(kind=8):: dels(0:5),volt(0:5)
      integer(ISZ):: numb(0:5)

c     --- Make sure that there is no data in conductors
      if (conductors%interior%nmax > 0 .or.
     &    conductors%evensubgrid%nmax > 0 .or.
     &    conductors%oddsubgrid%nmax > 0) then
        conductors%interior%nmax = 0
        conductors%evensubgrid%nmax = 0
        conductors%oddsubgrid%nmax = 0
        conductors%interior%n = 0
        conductors%evensubgrid%n = 0
        conductors%oddsubgrid%n = 0
        call ConductorTypefree(conductors)
      endif

c     --- If there are no grid points, then just return
      if (intercepts%nx < 0 .or. intercepts%ny < 0 .or. intercepts%nz < 0) then
        conductors%interior%n = 0
        conductors%evensubgrid%n = 0
        conductors%oddsubgrid%n = 0
        return
      endif

c     --- Create locals to make the coding below more readable
      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz
      ixlocal = intercepts%ix
      iylocal = intercepts%iy
      izlocal = intercepts%iz
      nxicpt = intercepts%nxicpt
      nyicpt = intercepts%nyicpt
      nzicpt = intercepts%nzicpt

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts
      xvoltages => intercepts%xvoltages
      yvoltages => intercepts%yvoltages
      zvoltages => intercepts%zvoltages
      xcondids => intercepts%xcondids
      ycondids => intercepts%ycondids
      zcondids => intercepts%zcondids

c     --- These are handy to precalculate
      xmmax = xmmin + nx*dx
      ymmax = ymmin + ny*dy
      zmmax = zmmin + nz*dz
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

c     --- Create locals to hold the counters
      ncmax = 0
      nemax = 0
      nomax = 0
      nc = 0
      ne = 0
      no = 0

c     --- Do this calculation twice. The first time, only count the number
c     --- of data points accumulated. Then allocate the space all at once.
c     --- The second time through saves the data.
c     --- Interestingly, it takes about as much time to allocate the arrays
c     --- as it does to do the calculation. It is faster to do the calculation
c     --- twice and do the allocation all at once than it is to do the
c     --- allocations incrementally during the calculation (which would
c     --- involve a fair amount of data copying, moving the existing data
c     --- from old arrays to new arrays).
      do iouter=1,2

        if (iouter == 2) then
c         --- Allocate the space all at once before the 2nd time through
          conductors%interior%nmax = ncmax
          conductors%evensubgrid%nmax = nemax
          conductors%oddsubgrid%nmax = nomax
c         --- Calling ConductorTypeallot has extra overhead and is slower
c         call ConductorTypeallot(conductors)
c         --- Direct allocation is slightly faster. This assumes that
c         --- the arrays are not already associated.
          allocate(conductors%interior%indx(0:2,ncmax))
          allocate(conductors%interior%volt(ncmax))
          allocate(conductors%interior%numb(ncmax))
          allocate(conductors%interior%ilevel(ncmax))
c         allocate(conductors%evensubgrid%prevphi(nemax))
          allocate(conductors%evensubgrid%indx(0:2,nemax))
          allocate(conductors%evensubgrid%dels(0:5,nemax))
          allocate(conductors%evensubgrid%volt(0:5,nemax))
          allocate(conductors%evensubgrid%numb(0:5,nemax))
c         allocate(conductors%evensubgrid%efield0(0:2,nemax))
c         allocate(conductors%evensubgrid%efieldd(0:5,nemax))
          allocate(conductors%evensubgrid%ilevel(nemax))
c         allocate(conductors%oddsubgrid%prevphi(nomax))
          allocate(conductors%oddsubgrid%indx(0:2,nomax))
          allocate(conductors%oddsubgrid%dels(0:5,nomax))
          allocate(conductors%oddsubgrid%volt(0:5,nomax))
          allocate(conductors%oddsubgrid%numb(0:5,nomax))
c         allocate(conductors%oddsubgrid%efield0(0:2,nomax))
c         allocate(conductors%oddsubgrid%efieldd(0:5,nomax))
          allocate(conductors%oddsubgrid%ilevel(nomax))
        endif

c       --- Find all of the grid points within the object.
c       --- Only one of the x, y and zintercepts needs to be scanned.
        do iz=0,nz
          zz = zmmin + iz*dz
          do iy=0,ny
            yy = ymmin + iy*dy

c           --- Loop over xinetercept points, looping over all points,
c           --- or stopping if the intercept is LARGEPOS (which may mean
c           --- that there are no intercepts along the current grid line).
            ii = 1
            do while (ii <= nxicpt .and. xintercepts(ii,iy,iz) < LARGEPOS)

c             --- x1 and x2 are the intercept points, with the object
c             --- in between
              x1 = xintercepts(ii,iy,iz)
              ii = ii + 1
              x2 = xintercepts(ii,iy,iz)

c             --- Check if the object overlaps the grid.
              if (x2 >= xmmin .and. x1 <= xmmax) then

c               --- Get the range in x grid points that are inside the object
                ix1 = int((max(x1-fuzz,xmmin) - xmmin)*dxi)
                ix2 = int((min(x2+fuzz,xmmax) - xmmin)*dxi)
                if (xmmin + ix1*dx <= x1 .and. ix1 < nx) ix1 = ix1 + 1
                if (xmmin + (ix2+1)*dx <= x2 .and. ix2 < nx) ix2 = ix2 + 1

                do ix=ix1,ix2

c                 --- Get the closest distance from the object surface to the
c                 --- grid point and only include the point if it is within
c                 --- dfill of the surface.
c                 --- An optimization would skip this calculation if dfill
c                 --- is LARGEPOS.
                  if (dfill < LARGEPOS) then
                    xx = xmmin + ix*dx
                    delx = min(xx - x1,x2 - xx)
                    dely = LARGEPOS
                    delz = LARGEPOS
                    do iy1=1,nyicpt
                      dely = min(dely,abs(yy - yintercepts(iy1,ix,iz)))
                      if (yintercepts(iy1,ix,iz) > yy) exit
                    enddo
                    do iz1=1,nzicpt
                      delz = min(delz,abs(zz - zintercepts(iz1,ix,iy)))
                      if (zintercepts(iz1,ix,iy) > zz) exit
                    enddo
                  else
                    delx = 0.
                    dely = 0.
                    delz = 0.
                  endif

                  if (delx < dfill*dx .or.
     &                dely < dfill*dy .or.
     &                delz < dfill*dz) then
                    if (iouter == 1) then
c                     --- The first time through the loop, only count the
c                     --- number of points
                      ncmax = ncmax + 1
                    else
c                     --- The second time, save the data
                      nc = nc + 1
                      conductors%interior%indx(0,nc) = ix + ixlocal
                      conductors%interior%indx(1,nc) = iy + iylocal
                      conductors%interior%indx(2,nc) = iz + izlocal
                      conductors%interior%volt(nc) = xvoltages(ii/2,iy,iz)
                      conductors%interior%numb(nc) = xcondids(ii/2,iy,iz)
                      conductors%interior%ilevel(nc) = intercepts%mglevel
                    endif
                  endif

                enddo

              endif
              ii = ii + 1

            enddo ! ii
          enddo ! iy
        enddo ! iz

c       --- Find the deltas. All three of the x, y, and zintercepts need
c       --- to be scanned. x is done first. When scanning y, a check is made
c       --- so that points already gathered from the x scan are skipped, and
c       --- similarly for z.

        do iz=0,nz
          zz = zmmin + iz*dz
          do iy=0,ny
            yy = ymmin + iy*dy

            ii = 1
            xside = +1.
            do while (ii <= nxicpt .and. xintercepts(ii,iy,iz) < LARGEPOS)

c             --- Check that the intercept point is within the grid
              if (xmmin < xintercepts(ii,iy,iz) .and.
     &            xintercepts(ii,iy,iz) < xmmax) then

                ix = (xintercepts(ii,iy,iz) - xmmin)*dxi
                if (xside > 0.) then
c                 --- delpx
c                 --- Make sure that ix is outside the object
                  xx = xmmin + ix*dx
                  if (ii == 1 .or. xx > xintercepts(ii-1,iy,iz)) then
                    delx = (xintercepts(ii,iy,iz) - xmmin)*dxi - ix
                  else
                    delx = LARGEPOS
                  endif
                  idel = 1
                else
c                 --- delmx
c                 --- Make sure that ix is outside the object
                  ix = ix + 1
                  xx = xmmin + ix*dx
                  if (ii == nxicpt .or. xx < xintercepts(ii+1,iy,iz)) then
                    delx = ix - (xintercepts(ii,iy,iz) - xmmin)*dxi
                  else
                    delx = LARGEPOS
                  endif
                  idel = 0
                endif

c               --- Make sure that delx is not zero. Note that delx would
c               --- never be negative. Also, skip points that are close
c               --- to one.
                if (0. < delx .and. delx < 1.-fuzz) then

                  if (iouter == 1) then
                    if (mod(ix+ixlocal+iy+iylocal+iz+izlocal,2) == 0) then
                      nemax = nemax + 1
                    else
                      nomax = nomax + 1
                    endif
                  else

c                   --- Save the data in temps at first. This separates the
c                   --- generation of the data from saving it, which needs
c                   --- check on whether the point is even of odd.
c                   --- This simplifies the code at the expense of copies.
                    dels = LARGEPOS
                    volt = 0.
                    numb = 0
                    indx(0) = ix + ixlocal
                    indx(1) = iy + iylocal
                    indx(2) = iz + izlocal
                    dels(idel) = delx
                    volt(idel) = xvoltages((ii+1)/2,iy,iz)
                    numb(idel) = xcondids((ii+1)/2,iy,iz)

c                   --- Scan over the yintercepts to see if the grid point
c                   --- is near any.
                    iy1 = 1
                    yside = +1.
                    do while (iy1 <= nyicpt .and.
     &                        yintercepts(iy1,ix,iz) < LARGEPOS)
                      if (yside > 0.) then
c                       --- delpy
                        if (iy1 == 1 .or. yy > yintercepts(iy1-1,ix,iz)) then
                          dely = (yintercepts(iy1,ix,iz) - yy)*dyi
                        else
                          dely = LARGEPOS
                        endif
                        idel = 3
                      else
c                       --- delmy
                        if (iy1 == nyicpt .or.
     &                      yy < yintercepts(iy1+1,ix,iz)) then
                          dely = (yy - yintercepts(iy1,ix,iz))*dyi
                        else
                          dely = LARGEPOS
                        endif
                        idel = 2
                      endif
                      if (0. < dely .and. dely < 1.-fuzz .and.
     &                    dely < dels(idel)) then
                        dels(idel) = dely
                        volt(idel) = yvoltages((iy1+1)/2,ix,iz)
                        numb(idel) = ycondids((iy1+1)/2,ix,iz)
                      endif
                      iy1 = iy1 + 1
                      yside = -yside
                    enddo

c                   --- Scan over the zintercepts to see if the grid point
c                   --- is near any.
                    iz1 = 1
                    zside = +1.
                    do while (iz1 <= nzicpt .and.
     &                        zintercepts(iz1,ix,iy) < LARGEPOS)
                      if (zside > 0.) then
c                       --- delpz
                        if (iz1 == 1 .or. zz > zintercepts(iz1-1,ix,iy)) then
                          delz = (zintercepts(iz1,ix,iy) - zz)*dzi
                        else
                          delz = LARGEPOS
                        endif
                        idel = 5
                      else
c                       --- delmz
                        if (iz1 == nzicpt .or.
     &                      zz < zintercepts(iz1+1,ix,iy)) then
                          delz = (zz - zintercepts(iz1,ix,iy))*dzi
                        else
                          delz = LARGEPOS
                        endif
                        idel = 4
                      endif
                      if (0. < delz .and. delz < 1.-fuzz .and.
     &                    delz < dels(idel)) then
                        dels(idel) = delz
                        volt(idel) = zvoltages((iz1+1)/2,ix,iy)
                        numb(idel) = zcondids((iz1+1)/2,ix,iy)
                      endif
                      iz1 = iz1 + 1
                      zside = -zside
                    enddo

c                   --- Now save the data in the appropriate place
                    if (mod(ix+ixlocal+iy+iylocal+iz+izlocal,2) == 0) then
                      ne = ne + 1
                      conductors%evensubgrid%indx(:,ne) = indx
                      conductors%evensubgrid%dels(:,ne) = dels
                      conductors%evensubgrid%volt(:,ne) = volt
                      conductors%evensubgrid%numb(:,ne) = numb
                      conductors%evensubgrid%ilevel(ne) = intercepts%mglevel
                    else
                      no = no + 1
                      conductors%oddsubgrid%indx(:,no) = indx
                      conductors%oddsubgrid%dels(:,no) = dels
                      conductors%oddsubgrid%volt(:,no) = volt
                      conductors%oddsubgrid%numb(:,no) = numb
                      conductors%oddsubgrid%ilevel(no) = intercepts%mglevel
                    endif
                  endif

                endif
              endif

              ii = ii + 1
              xside = -xside

            enddo
          enddo
        enddo

c       --- Scan over the yintercepts
        do iz=0,nz
          zz = zmmin + iz*dz
          do ix=0,nx
            xx = xmmin + ix*dx

            ii = 1
            yside = +1.
            do while (ii <= nyicpt .and. yintercepts(ii,ix,iz) < LARGEPOS)

              if (ymmin < yintercepts(ii,ix,iz) .and.
     &            yintercepts(ii,ix,iz) < ymmax) then

                iy = (yintercepts(ii,ix,iz) - ymmin)*dyi
                if (yside > 0.) then
c                 --- delpy
                  yy = ymmin + iy*dy
                  if (ii == 1 .or. yy > yintercepts(ii-1,ix,iz)) then
                    dely = (yintercepts(ii,ix,iz) - ymmin)*dyi - iy
                  else
                    dely = LARGEPOS
                  endif
                  idel = 3
                else
c                 --- delmy
                  iy = iy + 1
                  yy = ymmin + iy*dy
                  if (ii == nyicpt .or. yy < yintercepts(ii+1,ix,iz)) then
                    dely = iy - (yintercepts(ii,ix,iz) - ymmin)*dyi
                  else
                    dely = LARGEPOS
                  endif
                  idel = 2
                endif

c               --- Check if the del in x would be less then 1. If so, then
c               --- this point would have already been gathered when examining
c               --- x intercepts.
                ix1 = 1
                xside = +1.
                delx = 2.
                do while (ix1 <= nxicpt .and.
     &                    xintercepts(ix1,iy,iz) < LARGEPOS .and. delx >= 1.)
                  if (xside > 0.) then
c                   --- delpx
                    if (xintercepts(ix1,iy,iz) >= xx) then
                      delx = (xintercepts(ix1,iy,iz) - xx)*dxi
                    endif
                  else
c                   --- delmx
                    if (xintercepts(ix1,iy,iz) <= xx) then
                      delx = (xx - xintercepts(ix1,iy,iz))*dxi
                    endif
                  endif
                  ix1 = ix1 + 1
                  xside = -xside
                enddo

                if (0. < dely .and. dely < 1.-fuzz .and. delx >= 1.-fuzz) then

                  if (iouter == 1) then
                    if (mod(ix+ixlocal+iy+iylocal+iz+izlocal,2) == 0) then
                      nemax = nemax + 1
                    else
                      nomax = nomax + 1
                    endif
                  else

                    dels = LARGEPOS
                    volt = 0.
                    numb = 0
                    indx(0) = ix + ixlocal
                    indx(1) = iy + iylocal
                    indx(2) = iz + izlocal
                    dels(idel) = dely
                    volt(idel) = yvoltages((ii+1)/2,ix,iz)
                    numb(idel) = ycondids((ii+1)/2,ix,iz)

                    iz1 = 1
                    zside = +1.
                    do while (iz1 <= nzicpt .and.
     &                        zintercepts(iz1,ix,iy) < LARGEPOS)
                      if (zside > 0.) then
c                       --- delpz
                        if (iz1 == 1 .or. zz > zintercepts(iz1-1,ix,iy)) then
                          delz = (zintercepts(iz1,ix,iy) - zz)*dzi
                        else
                          delz = LARGEPOS
                        endif
                        idel = 5
                      else
c                       --- delmz
                        if (iz1 == nzicpt .or.
     &                      zz < zintercepts(iz1+1,ix,iy)) then
                          delz = (zz - zintercepts(iz1,ix,iy))*dzi
                        else
                          delz = LARGEPOS
                        endif
                        idel = 4
                      endif
                      if (0. < delz .and. delz < 1.-fuzz .and.
     &                    delz < dels(idel)) then
                        dels(idel) = delz
                        volt(idel) = zvoltages((iz1+1)/2,ix,iy)
                        numb(idel) = zcondids((iz1+1)/2,ix,iy)
                      endif
                      iz1 = iz1 + 1
                      zside = -zside
                    enddo

                    if (mod(ix+ixlocal+iy+iylocal+iz+izlocal,2) == 0) then
                      ne = ne + 1
                      conductors%evensubgrid%indx(:,ne) = indx
                      conductors%evensubgrid%dels(:,ne) = dels
                      conductors%evensubgrid%volt(:,ne) = volt
                      conductors%evensubgrid%numb(:,ne) = numb
                      conductors%evensubgrid%ilevel(ne) = intercepts%mglevel
                    else
                      no = no + 1
                      conductors%oddsubgrid%indx(:,no) = indx
                      conductors%oddsubgrid%dels(:,no) = dels
                      conductors%oddsubgrid%volt(:,no) = volt
                      conductors%oddsubgrid%numb(:,no) = numb
                      conductors%oddsubgrid%ilevel(no) = intercepts%mglevel
                    endif
                  endif

                endif
              endif

              ii = ii + 1
              yside = -yside

            enddo
          enddo
        enddo

c       --- Scan over zintercepts
        do iy=0,ny
          yy = ymmin + iy*dy
          do ix=0,nx
            xx = xmmin + ix*dx

            ii = 1
            zside = +1.
            do while (ii <= nzicpt .and. zintercepts(ii,ix,iy) < LARGEPOS)

              if (zmmin < zintercepts(ii,ix,iy) .and.
     &            zintercepts(ii,ix,iy) < zmmax) then

                iz = (zintercepts(ii,ix,iy) - zmmin)*dzi
                if (zside > 0.) then
c                 --- delpz
                  zz = zmmin + iz*dz
                  if (ii == 1 .or. zz > zintercepts(ii-1,ix,iy)) then
                    delz = (zintercepts(ii,ix,iy) - zmmin)*dzi - iz
                  else
                    delz = LARGEPOS
                  endif
                  idel = 5
                else
c                 --- delmz
                  iz = iz + 1
                  zz = zmmin + iz*dz
                  if (ii == nzicpt .or. zz < zintercepts(ii+1,ix,iy)) then
                    delz = iz - (zintercepts(ii,ix,iy) - zmmin)*dzi
                  else
                    delz = LARGEPOS
                  endif
                  idel = 4
                endif

c               --- Check if the del in x would be less then 1. If so, then
c               --- this point would have already been gathered when examining
c               --- x intercepts.
                ix1 = 1
                xside = +1.
                delx = 2.
                do while (ix1 <= nxicpt .and.
     &                    xintercepts(ix1,iy,iz) < LARGEPOS .and. delx >= 1.)
                  if (xside > 0.) then
c                   --- delpx
                    if (xintercepts(ix1,iy,iz) >= xx) then
                      delx = (xintercepts(ix1,iy,iz) - xx)*dxi
                    endif
                  else
c                   --- delmx
                    if (xintercepts(ix1,iy,iz) <= xx) then
                      delx = (xx - xintercepts(ix1,iy,iz))*dxi
                    endif
                  endif
                  ix1 = ix1 + 1
                  xside = -xside
                enddo

c               --- Check if the del in y would be less then 1. If so, then
c               --- this point would have already been gathered when examining
c               --- y intercepts.
                iy1 = 1
                yside = +1.
                dely = 2.
                do while (iy1 <= nyicpt .and.
     &                    yintercepts(iy1,ix,iz) < LARGEPOS .and. dely >= 1.)
                  if (yside > 0.) then
c                   --- delpy
                    if (yintercepts(iy1,ix,iz) >= yy) then
                      dely = (yintercepts(iy1,ix,iz) - yy)*dyi
                    endif
                  else
c                   --- delmy
                    if (yintercepts(iy1,ix,iz) <= yy) then
                      dely = (yy - yintercepts(iy1,ix,iz))*dyi
                    endif
                  endif
                  iy1 = iy1 + 1
                  yside = -yside
                enddo

                if (0. < delz .and. delz < 1.-fuzz .and.
     &              delx >= 1.-fuzz .and. dely >= 1.-fuzz) then

                  if (iouter == 1) then
                    if (mod(ix+ixlocal+iy+iylocal+iz+izlocal,2) == 0) then
                      nemax = nemax + 1
                    else
                      nomax = nomax + 1
                    endif
                  else

                    dels = LARGEPOS
                    volt = 0.
                    numb = 0
                    indx(0) = ix + ixlocal
                    indx(1) = iy + iylocal
                    indx(2) = iz + izlocal
                    dels(idel) = delz
                    volt(idel) = zvoltages((ii+1)/2,ix,iy)
                    numb(idel) = zcondids((ii+1)/2,ix,iy)
 
                    if (mod(ix+ixlocal+iy+iylocal+iz+izlocal,2) == 0) then
                      ne = ne + 1
                      conductors%evensubgrid%indx(:,ne) = indx
                      conductors%evensubgrid%dels(:,ne) = dels
                      conductors%evensubgrid%volt(:,ne) = volt
                      conductors%evensubgrid%numb(:,ne) = numb
                      conductors%evensubgrid%ilevel(ne) = intercepts%mglevel
                    else
                      no = no + 1
                      conductors%oddsubgrid%indx(:,no) = indx
                      conductors%oddsubgrid%dels(:,no) = dels
                      conductors%oddsubgrid%volt(:,no) = volt
                      conductors%oddsubgrid%numb(:,no) = numb
                      conductors%oddsubgrid%ilevel(no) = intercepts%mglevel
                    endif
                  endif

                endif
              endif

              ii = ii + 1
              zside = -zside

            enddo
          enddo
        enddo

      enddo ! iouter

      conductors%interior%n = nc
      conductors%evensubgrid%n = ne
      conductors%oddsubgrid%n = no

      return
      end
c=======================================================================
      subroutine intercepts_or(ileft,iright,iresult)
      use ConductorInterceptTypeModule
      type(ConductorInterceptType):: ileft,iright,iresult

      integer(ISZ):: lnicpt,rnicpt,nicpt
      real(kind=8),pointer:: lintercepts(:,:,:),rintercepts(:,:,:)
      real(kind=8),pointer:: lvoltages(:,:,:),rvoltages(:,:,:)
      integer(ISZ),pointer:: lcondids(:,:,:),rcondids(:,:,:)
      real(kind=8),pointer:: intercepts(:,:,:)
      real(kind=8),pointer:: voltages(:,:,:)
      integer(ISZ),pointer:: condids(:,:,:)

c     --- Make sure that the intercepts have the same grid information
      if (ileft%mglevel .ne. iright%mglevel .or.
     &    ileft%xmmin .ne. iright%xmmin .or.
     &    ileft%ymmin .ne. iright%ymmin .or.
     &    ileft%zmmin .ne. iright%zmmin .or.
     &    ileft%dx .ne. iright%dx .or.
     &    ileft%dy .ne. iright%dy .or.
     &    ileft%dz .ne. iright%dz .or.
     &    ileft%nx .ne. iright%nx .or.
     &    ileft%ny .ne. iright%ny .or.
     &    ileft%nz .ne. iright%nz .or.
     &    ileft%ix .ne. iright%ix .or.
     &    ileft%iy .ne. iright%iy .or.
     &    ileft%iz .ne. iright%iz) then
        call kaboom("The left and right intercepts grids do not match")
        return
      endif

      iresult%mglevel = ileft%mglevel
      iresult%xmmin = ileft%xmmin
      iresult%ymmin = ileft%ymmin
      iresult%zmmin = ileft%zmmin
      iresult%dx = ileft%dx
      iresult%dy = ileft%dy
      iresult%dz = ileft%dz
      iresult%nx = ileft%nx
      iresult%ny = ileft%ny
      iresult%nz = ileft%nz
      iresult%ix = ileft%ix
      iresult%iy = ileft%iy
      iresult%iz = ileft%iz

      iresult%nxicpt = ileft%nxicpt + iright%nxicpt
      iresult%nyicpt = ileft%nyicpt + iright%nyicpt
      iresult%nzicpt = ileft%nzicpt + iright%nzicpt
      call ConductorInterceptTypechange(iresult)

      lnicpt = ileft%nxicpt
      lintercepts => ileft%xintercepts
      lvoltages => ileft%xvoltages
      lcondids => ileft%xcondids
      rnicpt = iright%nxicpt
      rintercepts => iright%xintercepts
      rvoltages => iright%xvoltages
      rcondids => iright%xcondids
      nicpt = iresult%nxicpt
      intercepts => iresult%xintercepts
      voltages => iresult%xvoltages
      condids => iresult%xcondids
      call interceptor(ileft%ny,ileft%nz)

      lnicpt = ileft%nyicpt
      lintercepts => ileft%yintercepts
      lvoltages => ileft%yvoltages
      lcondids => ileft%ycondids
      rnicpt = iright%nyicpt
      rintercepts => iright%yintercepts
      rvoltages => iright%yvoltages
      rcondids => iright%ycondids
      nicpt = iresult%nyicpt
      intercepts => iresult%yintercepts
      voltages => iresult%yvoltages
      condids => iresult%ycondids
      call interceptor(ileft%nx,ileft%nz)

      lnicpt = ileft%nzicpt
      lintercepts => ileft%zintercepts
      lvoltages => ileft%zvoltages
      lcondids => ileft%zcondids
      rnicpt = iright%nzicpt
      rintercepts => iright%zintercepts
      rvoltages => iright%zvoltages
      rcondids => iright%zcondids
      nicpt = iresult%nzicpt
      intercepts => iresult%zintercepts
      voltages => iresult%zvoltages
      condids => iresult%zcondids
      call interceptor(ileft%nx,ileft%ny)

      return
      CONTAINS

        subroutine interceptor(n1,n2)
        integer(ISZ):: n1,n2

        integer(ISZ):: i1,i2
        integer(ISZ):: i0left,i0right,i0result
        real(kind=8):: p1,p2,volt
        integer(ISZ):: condid
        logical(ISZ):: ldone

        do i2=0,n2
          do i1=0,n1
            i0left = 1
            i0right = 1
            i0result = 1
            do while (i0left < lnicpt .and. i0right < rnicpt)

              if (lintercepts(i0left,i1,i2) <= rintercepts(i0right,i1,i2)) then
                p1 = lintercepts(i0left,i1,i2)
                p2 = lintercepts(i0left+1,i1,i2)
                volt = lvoltages((i0left+1)/2,i1,i2)
                condid = lcondids((i0left+1)/2,i1,i2)
                i0left = i0left + 2
              else
                p1 = rintercepts(i0right,i1,i2)
                p2 = rintercepts(i0right+1,i1,i2)
                volt = rvoltages((i0right+1)/2,i1,i2)
                condid = rcondids((i0right+1)/2,i1,i2)
                i0right = i0right + 2
              endif

              ldone = .false.
              do while (.not. ldone)

                  do while (i0right < rnicpt .and.
     &                      rintercepts(i0right+1,i1,i2) <= p2)
                    i0right = i0right + 2
                  enddo
                  if (i0right < rnicpt) then
                    if (rintercepts(i0right,i1,i2) < p2) then
                      p2 = rintercepts(i0right+1,i1,i2)
                      i0right = i0right + 2
                      ldone = .false.
                    else
                      ldone = .true.
                    endif
                  else
                    ldone = .true.
                  endif

                  do while (i0left < lnicpt .and.
     &                      lintercepts(i0left+1,i1,i2) <= p2)
                    i0left = i0left + 2
                  enddo
                  if (i0left < lnicpt) then
                    if (lintercepts(i0left,i1,i2) < p2) then
                      p2 = lintercepts(i0left+1,i1,i2)
                      i0left = i0left + 2
                      ldone = .false.
                    else
                      ldone = .true.
                    endif
                  else
                    ldone = .true.
                  endif

              enddo

              intercepts(i0result,i1,i2) = p1
              intercepts(i0result+1,i1,i2) = p2
              voltages((i0result+1)/2,i1,i2) = volt
              condids((i0result+1)/2,i1,i2) = condid
              i0result = i0result + 2

            enddo
            do while (i0left < lnicpt)
              intercepts(i0result,i1,i2) = lintercepts(i0left,i1,i2)
              intercepts(i0result+1,i1,i2) = lintercepts(i0left+1,i1,i2)
              voltages((i0result+1)/2,i1,i2) = lvoltages((i0left+1)/2,i1,i2)
              condids((i0result+1)/2,i1,i2) = lcondids((i0left+1)/2,i1,i2)
              i0left = i0left + 2
              i0result = i0result + 2
            enddo
            do while (i0right < rnicpt)
              intercepts(i0result,i1,i2) = rintercepts(i0right,i1,i2)
              intercepts(i0result+1,i1,i2) = rintercepts(i0right+1,i1,i2)
              voltages((i0result+1)/2,i1,i2) = rvoltages((i0right+1)/2,i1,i2)
              condids((i0result+1)/2,i1,i2) = rcondids((i0right+1)/2,i1,i2)
              i0right = i0right + 2
              i0result = i0result + 2
            enddo
            do while (i0result < nicpt)
              intercepts(i0result,i1,i2) = LARGEPOS
              intercepts(i0result+1,i1,i2) = LARGEPOS
              i0result = i0result + 2
            enddo
            
          enddo
        enddo

        return
        end subroutine interceptor

      end
c=======================================================================
      subroutine intercepts_and(ileft,iright,iresult)
      use ConductorInterceptTypeModule
      type(ConductorInterceptType):: ileft,iright,iresult

      integer(ISZ):: lnicpt,rnicpt,nicpt
      real(kind=8),pointer:: lintercepts(:,:,:),rintercepts(:,:,:)
      real(kind=8),pointer:: lvoltages(:,:,:),rvoltages(:,:,:)
      integer(ISZ),pointer:: lcondids(:,:,:),rcondids(:,:,:)
      real(kind=8),pointer:: intercepts(:,:,:)
      real(kind=8),pointer:: voltages(:,:,:)
      integer(ISZ),pointer:: condids(:,:,:)

c     --- Make sure that the intercepts have the same grid information
      if (ileft%mglevel .ne. iright%mglevel .or.
     &    ileft%xmmin .ne. iright%xmmin .or.
     &    ileft%ymmin .ne. iright%ymmin .or.
     &    ileft%zmmin .ne. iright%zmmin .or.
     &    ileft%dx .ne. iright%dx .or.
     &    ileft%dy .ne. iright%dy .or.
     &    ileft%dz .ne. iright%dz .or.
     &    ileft%nx .ne. iright%nx .or.
     &    ileft%ny .ne. iright%ny .or.
     &    ileft%nz .ne. iright%nz .or.
     &    ileft%ix .ne. iright%ix .or.
     &    ileft%iy .ne. iright%iy .or.
     &    ileft%iz .ne. iright%iz) then
        call kaboom("The left and right intercepts grids do not match")
        return
      endif

      iresult%mglevel = ileft%mglevel
      iresult%xmmin = ileft%xmmin
      iresult%ymmin = ileft%ymmin
      iresult%zmmin = ileft%zmmin
      iresult%dx = ileft%dx
      iresult%dy = ileft%dy
      iresult%dz = ileft%dz
      iresult%nx = ileft%nx
      iresult%ny = ileft%ny
      iresult%nz = ileft%nz
      iresult%ix = ileft%ix
      iresult%iy = ileft%iy
      iresult%iz = ileft%iz

      iresult%nxicpt = max(ileft%nxicpt,iright%nxicpt)
      iresult%nyicpt = max(ileft%nyicpt,iright%nyicpt)
      iresult%nzicpt = max(ileft%nzicpt,iright%nzicpt)
      call ConductorInterceptTypechange(iresult)

      lnicpt = ileft%nxicpt
      lintercepts => ileft%xintercepts
      lvoltages => ileft%xvoltages
      lcondids => ileft%xcondids
      rnicpt = iright%nxicpt
      rintercepts => iright%xintercepts
      rvoltages => iright%xvoltages
      rcondids => iright%xcondids
      nicpt = iresult%nxicpt
      intercepts => iresult%xintercepts
      voltages => iresult%xvoltages
      condids => iresult%xcondids
      call interceptand(ileft%ny,ileft%nz)

      lnicpt = ileft%nyicpt
      lintercepts => ileft%yintercepts
      lvoltages => ileft%yvoltages
      lcondids => ileft%ycondids
      rnicpt = iright%nyicpt
      rintercepts => iright%yintercepts
      rvoltages => iright%yvoltages
      rcondids => iright%ycondids
      nicpt = iresult%nyicpt
      intercepts => iresult%yintercepts
      voltages => iresult%yvoltages
      condids => iresult%ycondids
      call interceptand(ileft%nx,ileft%nz)

      lnicpt = ileft%nzicpt
      lintercepts => ileft%zintercepts
      lvoltages => ileft%zvoltages
      lcondids => ileft%zcondids
      rnicpt = iright%nzicpt
      rintercepts => iright%zintercepts
      rvoltages => iright%zvoltages
      rcondids => iright%zcondids
      nicpt = iresult%nzicpt
      intercepts => iresult%zintercepts
      voltages => iresult%zvoltages
      condids => iresult%zcondids
      call interceptand(ileft%nx,ileft%ny)

      return
      CONTAINS

        subroutine interceptand(n1,n2)
        integer(ISZ):: n1,n2

        integer(ISZ):: i1,i2
        integer(ISZ):: i0left,i0right,i0result
        real(kind=8):: p1,p2

        do i2=0,n2
          do i1=0,n1
            i0left = 1
            i0right = 1
            i0result = 1
            do while (i0left < lnicpt .and. i0right < rnicpt)

              p1 = max(lintercepts(i0left,i1,i2),
     &                 rintercepts(i0right,i1,i2))
              p2 = min(lintercepts(i0left+1,i1,i2),
     &                 rintercepts(i0right+1,i1,i2))

              if (p2 > p1) then
                intercepts(i0result,i1,i2) = p1
                intercepts(i0result+1,i1,i2) = p2
                voltages((i0result+1)/2,i1,i2) = lvoltages((i0left+1)/2,i1,i2)
                condids((i0result+1)/2,i1,i2) = lcondids((i0left+1)/2,i1,i2)
                i0result = i0result + 2
              endif

              if (lintercepts(i0left+1,i1,i2) < rintercepts(i0right+1,i1,i2)) then
                i0left = i0left + 2
              else if (rintercepts(i0right+1,i1,i2) < lintercepts(i0left+1,i1,i2)) then
                i0right = i0right + 2
              else
                i0left = i0left + 2
                i0right = i0right + 2
              endif

            enddo
            do while (i0result < nicpt)
              intercepts(i0result,i1,i2) = LARGEPOS
              intercepts(i0result+1,i1,i2) = LARGEPOS
              i0result = i0result + 2
            enddo
            
          enddo
        enddo

        return
        end subroutine interceptand

      end
c=======================================================================
      subroutine intercepts_not(ileft,iresult)
      use ConductorInterceptTypeModule
      type(ConductorInterceptType):: ileft,iresult

      integer(ISZ):: lnicpt,rnicpt,nicpt
      real(kind=8),pointer:: lintercepts(:,:,:)
      real(kind=8),pointer:: lvoltages(:,:,:)
      integer(ISZ),pointer:: lcondids(:,:,:)
      real(kind=8),pointer:: intercepts(:,:,:)
      real(kind=8),pointer:: voltages(:,:,:)
      integer(ISZ),pointer:: condids(:,:,:)

      iresult%mglevel = ileft%mglevel
      iresult%xmmin = ileft%xmmin
      iresult%ymmin = ileft%ymmin
      iresult%zmmin = ileft%zmmin
      iresult%dx = ileft%dx
      iresult%dy = ileft%dy
      iresult%dz = ileft%dz
      iresult%nx = ileft%nx
      iresult%ny = ileft%ny
      iresult%nz = ileft%nz
      iresult%ix = ileft%ix
      iresult%iy = ileft%iy
      iresult%iz = ileft%iz

      iresult%nxicpt = ileft%nxicpt + 2
      iresult%nyicpt = ileft%nyicpt + 2
      iresult%nzicpt = ileft%nzicpt + 2
      call ConductorInterceptTypechange(iresult)

      lnicpt = ileft%nxicpt
      lintercepts => ileft%xintercepts
      lvoltages => ileft%xvoltages
      lcondids => ileft%xcondids
      nicpt = iresult%nxicpt
      intercepts => iresult%xintercepts
      voltages => iresult%xvoltages
      condids => iresult%xcondids
      call interceptnot(ileft%ny,ileft%nz)

      lnicpt = ileft%nyicpt
      lintercepts => ileft%yintercepts
      lvoltages => ileft%yvoltages
      lcondids => ileft%ycondids
      nicpt = iresult%nyicpt
      intercepts => iresult%yintercepts
      voltages => iresult%yvoltages
      condids => iresult%ycondids
      call interceptnot(ileft%nx,ileft%nz)

      lnicpt = ileft%nzicpt
      lintercepts => ileft%zintercepts
      lvoltages => ileft%zvoltages
      lcondids => ileft%zcondids
      nicpt = iresult%nzicpt
      intercepts => iresult%zintercepts
      voltages => iresult%zvoltages
      condids => iresult%zcondids
      call interceptnot(ileft%nx,ileft%ny)

      return
      CONTAINS

        subroutine interceptnot(n1,n2)
        integer(ISZ):: n1,n2

        integer(ISZ):: i1,i2
        integer(ISZ):: i0left,i0result

        do i2=0,n2
          do i1=0,n1

            if (lintercepts(1,i1,i2) == -LARGEPOS) then
c             --- skip the first intercept point (skipping -LARGEPOS)
              i0left = 2
              i0result = 1
            else
c             --- prepend -LARGEPOS to the list of intercepts
              i0left = 1
              i0result = 2
              intercepts(1,i1,i2) = -LARGEPOS
              voltages((1+1)/2,i1,i2) = lvoltages((i0left+1)/2,i1,i2)
              condids((1+1)/2,i1,i2) = lcondids((i0left+1)/2,i1,i2)
            endif

c           --- copy the rest of the data over
            do while (i0left <= lnicpt)
              intercepts(i0result,i1,i2) = lintercepts(i0left,i1,i2)
              voltages((i0result+1)/2,i1,i2) = lvoltages((i0left+1)/2,i1,i2)
              condids((i0result+1)/2,i1,i2) = lcondids((i0left+1)/2,i1,i2)
              i0left = i0left + 1
              i0result = i0result + 1
            enddo
c           --- fill in any remainder with LARGEPOS
            do while (i0result <= nicpt)
              intercepts(i0result,i1,i2) = LARGEPOS
              intercepts(i0result+1,i1,i2) = LARGEPOS
              i0result = i0result + 1
            enddo
            
          enddo
        enddo

        return
        end subroutine interceptnot

      end
c=======================================================================
c=======================================================================
      subroutine setconductorparity(nn,ix,iy,iz,dels,parity,fuzz0,fuzz1,fuzzsign,dfill)
      integer(ISZ):: nn
      integer(ISZ):: ix(nn),iy(nn),iz(nn),parity(nn)
      real(kind=8):: dels(0:5,nn)
      real(kind=8):: fuzz0,fuzz1
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
        if (dels(0,i) < 1.+fuzzsign*fuzz1) then
          parity(i) = iparity
        elseif (dels(1,i) < 1.+fuzzsign*fuzz1) then
          parity(i) = iparity
        elseif (dels(2,i) < 1.+fuzzsign*fuzz1) then
          parity(i) = iparity
        elseif (dels(3,i) < 1.+fuzzsign*fuzz1) then
          parity(i) = iparity
        elseif (dels(4,i) < 1.+fuzzsign*fuzz1) then
          parity(i) = iparity
        elseif (dels(5,i) < 1.+fuzzsign*fuzz1) then
          parity(i) = iparity
        endif

c       --- If one direction is negative, all are negative. However,
c       --- because of the fuzz factor, each direction must be checked
c       --- since if one direction is between 0 and fuzz, the others
c       --- will all be positive.
        if (fuzz0 > 0.) then
          if (dels(0,i) < 0.+fuzz0) then
            parity(i) = -1
          elseif (dels(1,i) < 0.+fuzz0) then
            parity(i) = -1
          elseif (dels(2,i) < 0.+fuzz0) then
            parity(i) = -1
          elseif (dels(3,i) < 0.+fuzz0) then
            parity(i) = -1
          elseif (dels(4,i) < 0.+fuzz0) then
            parity(i) = -1
          elseif (dels(5,i) < 0.+fuzz0) then
            parity(i) = -1
          endif
        endif

c       --- To minimize the amount of conductor data, points deep inside the
c       --- conductor can be skipped. The depth is set by the input dfill.
c       --- The point must be at the depth in all directions to be skipped.
c       --- dfill is assumed to be passed in as a positive number.
        if (dels(0,i) < -dfill .and. dels(1,i) < -dfill .and.
     &      dels(2,i) < -dfill .and. dels(3,i) < -dfill .and.
     &      dels(4,i) < -dfill .and. dels(5,i) < -dfill) parity(i) = outside

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
      call gchange("ConductorGeometryVisualization",0)
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
      subroutine ZPlaneConductorFnew(zcent,zsign,
     &                               intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: zcent,zsign
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.
c A positive zsign means that the conductor is the region above zcent.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      integer(ISZ):: ix,iy,iz
      real(kind=8):: zz

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2
      intercepts%nyicpt = 2
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      do iz=0,nz
        zz = zmmin + iz*dz
        do iy=0,ny
          if ((zsign > 0. .and. zz > zcent) .or.
     &        (zsign < 0. .and. zz < zcent))then
c           --- inside the object
            xintercepts(1,iy,iz) = -LARGEPOS
            xintercepts(2,iy,iz) = +LARGEPOS
          else
c           --- outside the object
            xintercepts(1,iy,iz) = +LARGEPOS
            xintercepts(2,iy,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz
        do ix=0,nx
          if ((zsign > 0. .and. zz > zcent) .or.
     &        (zsign < 0. .and. zz < zcent))then
c           --- inside the object
            yintercepts(1,ix,iz) = -LARGEPOS
            yintercepts(2,ix,iz) = +LARGEPOS
          else
c           --- outside the object
            yintercepts(1,ix,iz) = +LARGEPOS
            yintercepts(2,ix,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        do ix=0,nx
          if (zsign > 0.) then
            zintercepts(1,ix,iy) = zcent
            zintercepts(2,ix,iy) = +LARGEPOS
          else
            zintercepts(1,ix,iy) = -LARGEPOS
            zintercepts(2,ix,iy) = zcent
          endif
        enddo
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

      distance = zsign*(z - zcent)

      return
      end
c=======================================================================
      subroutine ZPlaneIntercept(zcent,zsign,n,x,y,z,vx,vy,vz,xi,yi,zi,
     &                           itheta,iphi)
      use Constant
      real(kind=8):: zcent,zsign
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i

      do i=1,n
        if (vz(i) .ne. 0.) then
          xi(i) = vx(i)/vz(i)*(zcent - z(i)) + x(i)
          yi(i) = vy(i)/vz(i)*(zcent - z(i)) + y(i)
          zi(i) = zcent
          if (zsign > 0) then
            itheta(i) = pi
          else
            itheta(i) = 0.
          endif
          iphi(i) = 0.
        else
          xi(i) = LARGEPOS
          yi(i) = LARGEPOS
          zi(i) = zcent
          itheta(i) = 0.
          iphi(i) = 0.
        endif
      enddo

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
      subroutine PlaneConductorFnew(z0,zsign,theta,phi,xcent,ycent,zcent,
     &                              intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: z0,zsign,theta,phi,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.
c A positive zsign means that the conductor is the region above zcent.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      integer(ISZ):: ix,iy,iz
      real(kind=8):: ctheta,stheta,cphi,sphi,cci,sci,si,zp
      real(kind=8):: xx,yy,zz
      real(kind=8):: xi,yi,zi

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2
      intercepts%nyicpt = 2
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)
      if (ctheta*cphi .ne. 0.) cci = 1./(ctheta*cphi)
      if (stheta*cphi .ne. 0.) sci = 1./(stheta*cphi)
      if (sphi .ne. 0.) si  = 1./(sphi)

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          zp = yy*sphi + zz*ctheta*cphi
          if (stheta*cphi == 0.) then
            if ((zsign > 0. .and. zp > z0) .or.
     &          (zsign < 0. .and. zp < z0))then
c             --- inside the object
              xintercepts(1,iy,iz) = -LARGEPOS
              xintercepts(2,iy,iz) = +LARGEPOS
            else
c             --- outside the object
              xintercepts(1,iy,iz) = +LARGEPOS
              xintercepts(2,iy,iz) = +LARGEPOS
            endif
          else
            xi = (z0 - zz*ctheta*cphi - yy*sphi)*sci
            if (zsign > 0.) then
              xintercepts(1,iy,iz) = xi
              xintercepts(2,iy,iz) = +LARGEPOS
            else
              xintercepts(1,iy,iz) = -LARGEPOS
              xintercepts(2,iy,iz) = xi
            endif
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          zp = xx*stheta*cphi + zz*ctheta*cphi
          if (sphi == 0.) then
            if ((zsign > 0. .and. zp > z0) .or.
     &          (zsign < 0. .and. zp < z0))then
c             --- inside the object
              yintercepts(1,ix,iz) = -LARGEPOS
              yintercepts(2,ix,iz) = +LARGEPOS
            else
c             --- outside the object
              yintercepts(1,ix,iz) = +LARGEPOS
              yintercepts(2,ix,iz) = +LARGEPOS
            endif
          else
            yi = (z0 - xx*stheta*cphi - zz*ctheta*cphi)*si
            if (zsign > 0.) then
              yintercepts(1,ix,iz) = yi
              yintercepts(2,ix,iz) = +LARGEPOS
            else
              yintercepts(1,ix,iz) = -LARGEPOS
              yintercepts(2,ix,iz) = yi
            endif
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          zp = xx*stheta*cphi + yy*sphi
          if (ctheta*cphi == 0.) then
            if ((zsign > 0. .and. zp > z0) .or.
     &          (zsign < 0. .and. zp < z0)) then
c             --- inside the object
              zintercepts(1,ix,iy) = -LARGEPOS
              zintercepts(2,ix,iy) = +LARGEPOS
            else
c             --- outside the object
              zintercepts(1,ix,iy) = +LARGEPOS
              zintercepts(2,ix,iy) = +LARGEPOS
            endif
          else
            zi = (z0 - xx*stheta*cphi - yy*sphi)*cci
            if (zsign > 0.) then
              zintercepts(1,ix,iy) = zcent
              zintercepts(2,ix,iy) = +LARGEPOS
            else
              zintercepts(1,ix,iy) = -LARGEPOS
              zintercepts(2,ix,iy) = zcent
            endif
          endif
        enddo
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
      subroutine PlaneIntercept(z0,zsign,theta,phi,xcent,ycent,zcent,
     &                          n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: z0,zsign,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,vsq,tt,zp0,vzp1
      real(kind=8):: ctheta,stheta,cphi,sphi

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)

      do i=1,n

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) then
          xi(i) = LARGEPOS
          yi(i) = LARGEPOS
          zi(i) = LARGEPOS
          itheta(i) = 0.
          iphi(i) = 0.

        else

          xx = x(i) - xcent
          yy = y(i) - ycent
          zz = z(i) - zcent
          zp0 = xx*stheta*cphi + yy*sphi + zz*ctheta*cphi - z0
          vzp1 = vx(i)*stheta*cphi + vy(i)*sphi + vz(i)*ctheta*cphi
          if (vzp1 == 0.) then
            tt = LARGEPOS
          else
            tt = -zp0/vzp1
          endif

          xi(i) = vx(i)*tt + x(i)
          yi(i) = vy(i)*tt + y(i)
          zi(i) = vz(i)*tt + z(i)

c         --- Note that the theta and phi input do not have the same meaning
c         --- as the itheta and iphi output. The input are a rotation about
c         --- the y-axis and a rotation about the new x-axis. The output are
c         --- a rotation about the z-axis and a rotation about the new y-axis.
c         --- It was not meant to be confusing, but that is how it turned out.
          itheta(i) = acos(-ctheta*cphi)
          iphi(i) = atan2(-zsign*sphi,-zsign*stheta*cphi)

        endif
      enddo

      return
      end
c=======================================================================
      subroutine BoxConductorF(xsize,ysize,zsize,xcent,ycent,zcent,
     &                         n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                         fuzz)
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

      xmin = xcent - xsize*0.5
      xmax = xcent + xsize*0.5
      ymin = ycent - ysize*0.5
      ymax = ycent + ysize*0.5
      zmin = zcent - zsize*0.5
      zmax = zcent + zsize*0.5

      do i=1,n

        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

        if (xmin - fuzz <= x(i) .and. x(i) <= xmax + fuzz .and.
     &      ymin - fuzz <= y(i) .and. y(i) <= ymax + fuzz) then

          if (z(i) < zmin - fuzz) then
            delpz(i) = zmin - z(i)
          else if (zmin - fuzz <= z(i) .and. z(i) <= zmax + fuzz) then
            delmz(i) = zmin - z(i)
            delpz(i) = z(i) - zmax
          else if (zmax + fuzz < z(i)) then
            delmz(i) = z(i) - zmax
          endif
          if (0 < delmz(i) .and. delmz(i) < fuzz) delmz(i) = 0.
          if (0 < delpz(i) .and. delpz(i) < fuzz) delpz(i) = 0.

        endif

        if (ymin - fuzz <= y(i) .and. y(i) <= ymax + fuzz .and.
     &      zmin - fuzz <= z(i) .and. z(i) <= zmax + fuzz) then

          if (x(i) < xmin - fuzz) then
            delpx(i) = xmin - x(i)
          else if (xmin - fuzz <= x(i) .and. x(i) <= xmax + fuzz) then
            delmx(i) = xmin - x(i)
            delpx(i) = x(i) - xmax
          else if (xmax + fuzz < x(i)) then
            delmx(i) = x(i) - xmax
          endif
          if (0 < delmx(i) .and. delmx(i) < fuzz) delmx(i) = 0.
          if (0 < delpx(i) .and. delpx(i) < fuzz) delpx(i) = 0.

        endif

        if (zmin - fuzz <= z(i) .and. z(i) <= zmax + fuzz .and.
     &      xmin - fuzz <= x(i) .and. x(i) <= xmax + fuzz) then

          if (y(i) < ymin - fuzz) then
            delpy(i) = ymin - y(i)
          else if (ymin - fuzz <= y(i) .and. y(i) <= ymax + fuzz) then
            delmy(i) = ymin - y(i)
            delpy(i) = y(i) - ymax
          else if (ymax + fuzz < y(i)) then
            delmy(i) = y(i) - ymax
          endif
          if (0 < delmy(i) .and. delmy(i) < fuzz) delmy(i) = 0.
          if (0 < delpy(i) .and. delpy(i) < fuzz) delpy(i) = 0.

        endif

      enddo

      return
      end
c=======================================================================
      subroutine BoxConductorFnew(xsize,ysize,zsize,xcent,ycent,zcent,
     &                            intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: xsize,ysize,zsize,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2
      intercepts%nyicpt = 2
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          if (abs(yy) <= ysize*0.5 .and. abs(zz) <= zsize*0.5) then
            xintercepts(1,iy,iz) = -xsize*0.5 + xcent
            xintercepts(2,iy,iz) = +xsize*0.5 + xcent
          else
            xintercepts(1,iy,iz) = +LARGEPOS
            xintercepts(2,iy,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          if (abs(xx) <= xsize*0.5 .and. abs(zz) <= zsize*0.5) then
            yintercepts(1,ix,iz) = -ysize*0.5 + ycent
            yintercepts(2,ix,iz) = +ysize*0.5 + ycent
          else
            yintercepts(1,ix,iz) = +LARGEPOS
            yintercepts(2,ix,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          if (abs(xx) <= xsize*0.5 .and. abs(yy) <= ysize*0.5) then
            zintercepts(1,ix,iy) = -zsize*0.5 + zcent
            zintercepts(2,ix,iy) = +zsize*0.5 + zcent
          else
            zintercepts(1,ix,iy) = +LARGEPOS
            zintercepts(2,ix,iy) = +LARGEPOS
          endif
        enddo
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

      integer(ISZ):: i
      real(kind=8):: xd,yd,zd,dd

      do i=1,n
        xd = abs(x(i) - xcent) - xsize*0.5
        yd = abs(y(i) - ycent) - ysize*0.5
        zd = abs(z(i) - zcent) - zsize*0.5

        dd = sqrt(max(0.,xd)**2 + max(0.,yd)**2 + max(0.,zd)**2)

        if (dd == 0.) dd = max(xd,yd,zd)

        distance(i) = dd
      enddo

      return
      end
c=======================================================================
      subroutine BoxIntercept(xsize,ysize,zsize,xcent,ycent,zcent,
     &                        n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: xsize,ysize,zsize,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: vsq
      real(kind=8):: distance
      real(kind=8):: tx,txm,txp,xix,yix,zix,dix
      real(kind=8):: ty,tym,typ,xiy,yiy,ziy,diy
      real(kind=8):: tz,tzm,tzp,xiz,yiz,ziz,diz
      real(kind=8):: ctheta,stheta,cphi,sphi
      real(kind=8):: xmin,xmax,ymin,ymax,zmin,zmax

      xmin = xcent - xsize*0.5
      xmax = xcent + xsize*0.5
      ymin = ycent - ysize*0.5
      ymax = ycent + ysize*0.5
      zmin = zcent - zsize*0.5
      zmax = zcent + zsize*0.5

      do i=1,n

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) then
          xi(i) = LARGEPOS
          yi(i) = LARGEPOS
          zi(i) = LARGEPOS
          itheta(i) = 0.
          iphi(i) = 0.
          cycle
        endif

        if (vx(i) == 0.) then
          tx = LARGEPOS
        else
          txm = -(xcent - xsize*0.5 - x(i))/vx(i)
          txp = -(xcent + xsize*0.5 - x(i))/vx(i)
          if (txm > 0. .and. txp > 0.) then
            tx = min(txm,txp)
          else
            tx = max(txm,txp)
          endif
        endif
        xix = -vx(i)*tx + x(i)
        yix = -vy(i)*tx + y(i)
        zix = -vz(i)*tx + z(i)
        if (tx > 0. .and.
     &      ymin <= yix .and. yix <= ymax .and.
     &      zmin <= zix .and. zix <= zmax) then
          dix = (xix - x(i))**2 + (yix - y(i))**2 + (zix - z(i))**2
        else
          dix = LARGEPOS
        endif

        if (vy(i) == 0.) then
          ty = LARGEPOS
        else
          tym = -(ycent - ysize*0.5 - y(i))/vy(i)
          typ = -(ycent + ysize*0.5 - y(i))/vy(i)
          if (tym > 0. .and. typ > 0.) then
            ty = min(tym,typ)
          else
            ty = max(tym,typ)
          endif
        endif
        xiy = -vx(i)*ty + x(i)
        yiy = -vy(i)*ty + y(i)
        ziy = -vz(i)*ty + z(i)
        if (ty > 0. .and.
     &      xmin <= xiy .and. xiy <= xmax .and.
     &      zmin <= ziy .and. ziy <= zmax) then
          diy = (xiy - x(i))**2 + (yiy - y(i))**2 + (ziy - z(i))**2
        else
          diy = LARGEPOS
        endif

        if (vz(i) == 0.) then
          tz = LARGEPOS
        else
          tzm = -(zcent - zsize*0.5 - z(i))/vz(i)
          tzp = -(zcent + zsize*0.5 - z(i))/vz(i)
          if (tzm > 0. .and. tzp > 0.) then
            tz = min(tzm,tzp)
          else
            tz = max(tzm,tzp)
          endif
        endif
        xiz = -vx(i)*tz + x(i)
        yiz = -vy(i)*tz + y(i)
        ziz = -vz(i)*tz + z(i)
        if (tz > 0. .and.
     &      xmin <= xiz .and. xiz <= xmax .and.
     &      ymin <= yiz .and. yiz <= ymax) then
          diz = (xiz - x(i))**2 + (yiz - y(i))**2 + (ziz - z(i))**2
        else
          diz = LARGEPOS
        endif

        if (dix < diy) then
          xi(i) = xix
          yi(i) = yix
          zi(i) = zix
          distance = dix
          if (xix > xcent) then
            itheta(i) = +pi*0.5
          else
            itheta(i) = -pi*0.5
          endif
          iphi(i) = 0.
        else if (dix > diy) then
          xi(i) = xiy
          yi(i) = yiy
          zi(i) = ziy
          distance = diy
          itheta(i) = pi*0.5
          if (yiy > ycent) then
            iphi(i) = +pi*0.5
          else
            iphi(i) = -pi*0.5
          endif
        else
          xi(i) = LARGEPOS
          yi(i) = LARGEPOS
          zi(i) = LARGEPOS
          itheta(i) = 0.
          iphi(i) = 0.
          distance = LARGEPOS
        endif
        if (diz < distance) then
          xi(i) = xiz
          yi(i) = yiz
          zi(i) = ziz
          if (ziz > zcent) then
            itheta(i) = 0.
          else
            itheta(i) = pi
          endif
          iphi(i) = 0.
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

      zmin = zcent - length*0.5
      zmax = zcent + length*0.5

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
          delmz(i) = min(0.,zmin - z(i))
          delpz(i) = min(0.,z(i) - zmax)
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
        if (xx > 0.) then
          delmx(i) = -(sqrt(max(0.,rad**2 - yy**2)) + xx)
          delpx(i) = 0.
        else
          delmx(i) = 0.
          delpx(i) = -(sqrt(max(0.,rad**2 - yy**2)) - xx)
        endif
        if (yy > 0.) then
          delmy(i) = -(sqrt(max(0.,rad**2 - xx**2)) + yy)
          delpy(i) = 0.
        else
          delmy(i) = 0.
          delpy(i) = -(sqrt(max(0.,rad**2 - xx**2)) - yy)
        endif
        delmz(i) = min(0.,zmin - z(i))
        delpz(i) = min(0.,z(i) - zmax)

      enddo

      return
      end
c=======================================================================
      subroutine ZCylinderConductorFnew(rad,length,xcent,ycent,zcent,
     &                                  intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: rad,length,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      real(kind=8):: zmin,zmax
      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz,rsq

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2
      intercepts%nyicpt = 2
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      zmin = -length*0.5
      zmax = +length*0.5

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz .and.
     &        abs(yy) <= rad) then
            xintercepts(1,iy,iz) = xcent - sqrt(max(0.,rad**2 - yy**2))
            xintercepts(2,iy,iz) = xcent + sqrt(max(0.,rad**2 - yy**2))
          else
            xintercepts(1:2,iy,iz) = LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz .and.
     &        abs(xx) <= rad) then
            yintercepts(1,ix,iz) = ycent - sqrt(max(0.,rad**2 - xx**2))
            yintercepts(2,ix,iz) = ycent + sqrt(max(0.,rad**2 - xx**2))
          else
            yintercepts(1:2,ix,iz) = LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rsq = xx**2 + yy**2
          if (rsq <= rad**2) then
            zintercepts(1,ix,iy) = zcent + zmin
            zintercepts(2,ix,iy) = zcent + zmax
          else
            zintercepts(1:2,ix,iy) = LARGEPOS
          endif
        enddo
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
        zz = abs(z(i) - zcent) - length*0.5

        dd = sqrt(max(0.,rr)**2 + max(0.,zz)**2)
        if (dd == 0.) dd = max(rr,zz)
        distance(i) = dd
      enddo

      return
      end
c=======================================================================
      subroutine ZCylinderIntercept(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: zmin,zmax,vsq,tt,theta,phi
      real(kind=8):: ttemp,xx,yy,zz,rr,rrsq,t1,t2,a,b,c

      zmin = zcent - length*0.5
      zmax = zcent + length*0.5

      do i=1,n

        xi(i) = LARGEPOS
        yi(i) = LARGEPOS
        zi(i) = LARGEPOS
        itheta(i) = 0.
        iphi(i) = 0.

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

        tt = LARGEPOS
        theta = 0.
        phi = 0.

c       --- First, check if trajectory hits conductor end.
        if (vz(i) .ne. 0.) then
          ttemp = -(zmin - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = (xx - xcent)**2 + (yy - ycent)**2
          if (rrsq <= rad**2) then
            if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &          (ttemp > tt .and. tt < 0.)) then
              tt = ttemp
              theta = pi
              phi = 0.
            endif
          endif

          ttemp = -(zmax - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = (xx - xcent)**2 + (yy - ycent)**2
          if (rrsq <= rad**2) then
            if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &          (ttemp > tt .and. tt < 0.)) then
              tt = ttemp
              theta = 0.
              phi = 0.
            endif
          endif
        endif

c       --- Check if it hits the cylinder.
        a = vx(i)**2 + vy(i)**2
        if (a > 0.) then
          b = -2*(vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
          c = (x(i)-xcent)**2 + (y(i)-ycent)**2 - rad**2
          if ((b**2 - 4*a*c) > 0.) then
            t1 = (-b-sqrt(b**2 - 4*a*c))/(2*a)
            t2 = (-b+sqrt(b**2 - 4*a*c))/(2*a)
            if (t1 > 0. .and. t2 > 0.) then
              ttemp = min(t1,t2)
            else
              ttemp = max(t1,t2)
            endif
            zz = -vz(i)*ttemp + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &            (ttemp > tt .and. tt < 0.)) then
                tt = ttemp
                xx = -vx(i)*ttemp + x(i) - xcent
                yy = -vy(i)*ttemp + y(i) - ycent
                theta = pi*0.5
                phi = atan2(yy,xx)
              endif
            endif
          endif
        endif

        if (tt == LARGEPOS) cycle

        xi(i) = -vx(i)*tt + x(i)
        yi(i) = -vy(i)*tt + y(i)
        zi(i) = -vz(i)*tt + z(i)
        itheta(i) = theta
        iphi(i) = phi

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

      zmin = zcent - length*0.5
      zmax = zcent + length*0.5

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
            delmz(i) = min(0.,zmin - z(i))
            delpz(i) = min(0.,z(i) - zmax)
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
        if (xx > 0.) then
          delmx(i) = -(sqrt(max(0.,r**2 - yy**2)) + xx)
          delpx(i) = 0.
        else
          delmx(i) = 0.
          delpx(i) = -(sqrt(max(0.,r**2 - yy**2)) - xx)
        endif
        if (yy > 0.) then
          delmy(i) = -(sqrt(max(0.,r**2 - xx**2)) + yy)
          delpy(i) = 0.
        else
          delmy(i) = 0.
          delpy(i) = -(sqrt(max(0.,r**2 - xx**2)) - yy)
        endif
        delmz(i) = min(0.,zmin - z(i))
        delpz(i) = min(0.,z(i) - zmax)

      enddo

      return
      end
c=======================================================================
      subroutine ZRoundedCylinderConductorFnew(rad,length,rad2,
     &                                         xcent,ycent,zcent,
     &                                         intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      real(kind=8):: zmin,zmax
      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz,rsq,rr,ri

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2
      intercepts%nyicpt = 2
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      zmin = -length*0.5
      zmax = +length*0.5

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz) then
            if (zmin-fuzz <= zz .and. zz < zmin+rad2) then
              rr = rad - rad2 + sqrt(max(0.,rad2**2 - (zz - (zmin+rad2))**2))
            else if (zmax-rad2 < zz .and. zz <= zmax+fuzz) then
              rr = rad - rad2 + sqrt(max(0.,rad2**2 - (zz - (zmax-rad2))**2))
            else
              rr = rad
            endif
            if (abs(yy) <= rr) then
              xintercepts(1,iy,iz) = xcent - sqrt(max(0.,rr**2 - yy**2))
              xintercepts(2,iy,iz) = xcent + sqrt(max(0.,rr**2 - yy**2))
            else
              xintercepts(1:2,iy,iz) = +LARGEPOS
            endif
          else
            xintercepts(1:2,iy,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz) then
            if (zmin-fuzz <= zz .and. zz < zmin+rad2) then
              rr = rad - rad2 + sqrt(max(0.,rad2**2 - (zz - (zmin+rad2))**2))
            else if (zmax-rad2 < zz .and. zz <= zmax+fuzz) then
              rr = rad - rad2 + sqrt(max(0.,rad2**2 - (zz - (zmax-rad2))**2))
            else
              rr = rad
            endif
            if (abs(xx) <= rr) then
              yintercepts(1,ix,iz) = ycent - sqrt(max(0.,rr**2 - xx**2))
              yintercepts(2,ix,iz) = ycent + sqrt(max(0.,rr**2 - xx**2))
            else
              yintercepts(1:2,ix,iz) = +LARGEPOS
            endif
          else
            yintercepts(1:2,ix,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rsq = xx**2 + yy**2
          if ((rad-rad2)**2 < rsq .and. rsq <= rad**2) then
            rr = sqrt(rsq)
            ri = rad2 - sqrt(max(0.,rad2**2 - (rad-rad2-rr)**2))
            zintercepts(1,ix,iy) = zcent + zmin + ri
            zintercepts(2,ix,iy) = zcent + zmax - ri
          else if (rsq <= (rad-rad2)**2) then
            zintercepts(1,ix,iy) = zcent + zmin
            zintercepts(2,ix,iy) = zcent + zmax
          else
            zintercepts(1:2,ix,iy) = LARGEPOS
          endif
        enddo
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

      z2 = max(length*0.5-rad2,0.)

      do i=1,n

        rr = sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2)
        zz = abs(z(i) - zcent)
        if (zz >= z2 .and. rr <= rad-rad2) then
          distance(i) = zz - length*0.5
        else if (zz <= z2) then
          distance(i) = max(rr - rad,zz - length*0.5)
        else
          distance(i) = sqrt((zz - z2)**2 + (rr - (rad - rad2))**2) - rad2
        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZRoundedCylinderIntercept(rad,length,rad2,xcent,ycent,zcent,
     &                                     n,x,y,z,vx,vy,vz,xi,yi,zi,
     &                                     itheta,iphi)
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: ri2,rr,zz,dr

c     --- first calculate for non-rounded cylinder and then readjust
c     --- the ones that are affected by rounded corners
      call ZCylinderIntercept(rad,length,xcent,ycent,zcent,
     &                              n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)

      do i=1,n
c       --- figure out which one are in the corners

c       --- between the rounded corners in Z => not affected
        if( zi(i)-zcent < 0.5*length-rad2 ) cycle

c       --- between the rounded corners in R => not affected
        ri2 =  (xi(i)-xcent)**2 + (yi(i)-ycent)**2
        if( ri2 < (rad-rad2)**2 ) cycle

c       --- these are the ones we need to change

c          R-Z diagram of top left corner:

c               particle
c           +------+------------------------------+
c           |      |-\        -----/.....|........|
c           |      |  \    --/...........|........|
c           |   dr |   \+-/..............|........|
c           |      |  -/+\...............|........|
c           |      |-/..|.\..............|........|
c           |      |....|..\.............|........|   ..   inside of a
c        rr |    -/|....|...-\...........|........|   ..   rounded cylinder
c           |   /..|....|.....\.rad2.....|........|
c           |  /...|....|......\.........|........|
c           | /....|....|.......-\.......|........|
c           | |....|....|.........\......|........|
c           |/.....|....|..........\...------.....|
c           ||.....|....|...........-\/..|...\....|
c           |......|....|.............\.itheta\...|
c           |......|.dz.|..............\.|.....\..|
c           |......|....|...............\|......|.|
c           +------+----+-------------------------\
c                   <-------- zz (<0) --->

        rr = sqrt(ri2) - (rad-rad2)
        zz = zi(i)-zcent - (0.5*length-rad2)

c       --- calculate angles
        iphi(i)   = atan2(yi(i)-ycent,xi(i)-xcent)
        itheta(i) = atan2(rr,zz)

c       --- adjust position
        dr    = rr - sin(itheta(i))*rad2
        xi(i) = xi(i) - dr*cos(iphi(i))
        yi(i) = yi(i) - dr*sin(iphi(i))
        zi(i) = zi(i) - (zz - cos(itheta(i))*rad2)
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

      zmin = zcent - length*0.5
      zmax = zcent + length*0.5

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
          delmz(i) = min(0.,zmin - z(i))
          delpz(i) = min(0.,z(i) - zmax)
          cycle
        endif

c       --- The point is on the cylinder exactly.
        if (xx > 0.) then
          delmx(i) = 0.
        else
          delpx(i) = 0.
        endif
        if (yy > 0.) then
          delmy(i) = 0.
        else
          delpy(i) = 0.
        endif
        delmz(i) = min(0.,zmin - z(i))
        delpz(i) = min(0.,z(i) - zmax)

      enddo

      return
      end
c=======================================================================
      subroutine ZCylinderOutConductorFnew(rad,length,xcent,ycent,zcent,
     &                                     intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: rad,length,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      real(kind=8):: zmin,zmax
      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz,rsq

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 4
      intercepts%nyicpt = 4
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      zmin = -length*0.5
      zmax = +length*0.5

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz) then
            if (abs(yy) <= rad) then
              xintercepts(1,iy,iz) = -LARGEPOS
              xintercepts(2,iy,iz) = xcent - sqrt(max(0.,rad**2 - yy**2))
              xintercepts(3,iy,iz) = xcent + sqrt(max(0.,rad**2 - yy**2))
              xintercepts(4,iy,iz) = +LARGEPOS
            else
              xintercepts(1,iy,iz) = -LARGEPOS
              xintercepts(2:4,iy,iz) = +LARGEPOS
            endif
          else
            xintercepts(1:4,iy,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz) then
            if (abs(xx) <= rad) then
              yintercepts(1,ix,iz) = -LARGEPOS
              yintercepts(2,ix,iz) = ycent - sqrt(max(0.,rad**2 - xx**2))
              yintercepts(3,ix,iz) = ycent + sqrt(max(0.,rad**2 - xx**2))
              yintercepts(4,ix,iz) = +LARGEPOS
            else
              yintercepts(1,ix,iz) = -LARGEPOS
              yintercepts(2:4,ix,iz) = +LARGEPOS
            endif
          else
            yintercepts(1:4,ix,iz) = LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rsq = xx**2 + yy**2
          if (rsq >= rad**2) then
            zintercepts(1,ix,iy) = zcent + zmin
            zintercepts(2,ix,iy) = zcent + zmax
          else
            zintercepts(1:2,ix,iy) = LARGEPOS
          endif
        enddo
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
        zz = abs(z(i) - zcent) - length*0.5

        dd = sqrt(max(0.,rr)**2 + max(0.,zz)**2)
        if (dd == 0.) dd = max(rr,zz)
        distance(i) = dd
      enddo

      return
      end
c=======================================================================
      subroutine ZCylinderOutIntercept(rad,length,xcent,ycent,zcent,
     &                                 n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: rad,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: zmin,zmax,vsq,tt,theta,phi
      real(kind=8):: ttemp,xx,yy,zz,rr,rrsq,t1,t2,a,b,c

      zmin = zcent - length*0.5
      zmax = zcent + length*0.5

      do i=1,n

        xi(i) = LARGEPOS
        yi(i) = LARGEPOS
        zi(i) = LARGEPOS
        itheta(i) = 0.
        iphi(i) = 0.

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

        tt = LARGEPOS
        theta = 0.
        phi = 0.

c       --- First, check if trajectory hits conductor end.
        if (vz(i) .ne. 0.) then
          ttemp = -(zmin - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = (xx - xcent)**2 + (yy - ycent)**2
          if (rrsq >= rad**2) then
            if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &          (ttemp > tt .and. tt < 0.)) then
              tt = ttemp
              theta = pi
              phi = 0.
            endif
          endif

          ttemp = -(zmax - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = (xx - xcent)**2 + (yy - ycent)**2
          if (rrsq >= rad**2) then
            if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &          (ttemp > tt .and. tt < 0.)) then
              tt = ttemp
              theta = 0.
              phi = 0.
            endif
          endif
        endif

c       --- Check if it hits the cylinder.
        a = vx(i)**2 + vy(i)**2
        if (a > 0.) then
          b = -2*(vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
          c = (x(i)-xcent)**2 + (y(i)-ycent)**2 - rad**2
          if ((b**2 - 4*a*c) > 0.) then
            t1 = (-b-sqrt(b**2 - 4*a*c))/(2*a)
            t2 = (-b+sqrt(b**2 - 4*a*c))/(2*a)
            if (t1 > 0. .and. t2 > 0.) then
              ttemp = min(t1,t2)
            else
              ttemp = max(t1,t2)
            endif
            zz = -vz(i)*ttemp + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &            (ttemp > tt .and. tt < 0.)) then
                tt = ttemp
                xx = -vx(i)*ttemp + x(i) - xcent
                yy = -vy(i)*ttemp + y(i) - ycent
                theta = pi*0.5
                phi = atan2(-yy,-xx)
              endif
            endif
          endif
        endif

        if (tt == LARGEPOS) cycle

        xi(i) = -vx(i)*tt + x(i)
        yi(i) = -vy(i)*tt + y(i)
        zi(i) = -vz(i)*tt + z(i)
        itheta(i) = theta
        iphi(i) = phi

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

      zmin = zcent - length*0.5
      zmax = zcent + length*0.5

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
            delmz(i) = min(0.,zmin - z(i))
            delpz(i) = min(0.,z(i) - zmax)
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
        if (xx > 0.) then
          delmx(i) = 0.
        else
          delpx(i) = 0.
        endif
        if (yy > 0.) then
          delmy(i) = 0.
        else
          delpy(i) = 0.
        endif
        delmz(i) = min(0.,zmin - z(i))
        delpz(i) = min(0.,z(i) - zmax)

      enddo

      return
      end
c=======================================================================
      subroutine ZRoundedCylinderOutConductorFnew(rad,length,rad2,
     &                                            xcent,ycent,zcent,
     &                                            intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      real(kind=8):: zmin,zmax
      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz,rsq,rr,ri

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 4
      intercepts%nyicpt = 4
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      zmin = -length*0.5
      zmax = +length*0.5

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz) then
            if (zmin-fuzz <= zz .and. zz < zmin+rad2) then
              rr = rad + rad2 - sqrt(max(0.,rad2**2 - (zz - (zmin+rad2))**2))
            else if (zmax-rad2 < zz .and. zz <= zmax+fuzz) then
              rr = rad + rad2 - sqrt(max(0.,rad2**2 - (zz - (zmax-rad2))**2))
            else
              rr = rad
            endif
            if (abs(yy) <= rr) then
              xintercepts(1,iy,iz) = -LARGEPOS
              xintercepts(2,iy,iz) = xcent - sqrt(max(0.,rr**2 - yy**2))
              xintercepts(3,iy,iz) = xcent + sqrt(max(0.,rr**2 - yy**2))
              xintercepts(4,iy,iz) = +LARGEPOS
            else
              xintercepts(1,iy,iz) = -LARGEPOS
              xintercepts(2:4,iy,iz) = +LARGEPOS
            endif
          else
            xintercepts(1:4,iy,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          if (zmin-fuzz <= zz .and. zz <= zmax+fuzz) then
            if (zmin-fuzz <= zz .and. zz < zmin+rad2) then
              rr = rad + rad2 - sqrt(max(0.,rad2**2 - (zz - (zmin+rad2))**2))
            else if (zmax-rad2 < zz .and. zz <= zmax+fuzz) then
              rr = rad + rad2 - sqrt(max(0.,rad2**2 - (zz - (zmax-rad2))**2))
            else
              rr = rad
            endif
            if (abs(xx) <= rr) then
              yintercepts(1,ix,iz) = -LARGEPOS
              yintercepts(2,ix,iz) = ycent - sqrt(max(0.,rr**2 - xx**2))
              yintercepts(3,ix,iz) = ycent + sqrt(max(0.,rr**2 - xx**2))
              yintercepts(4,ix,iz) = +LARGEPOS
            else
              yintercepts(1,ix,iz) = -LARGEPOS
              yintercepts(2:4,ix,iz) = +LARGEPOS
            endif
          else
            yintercepts(1:4,ix,iz) = LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rsq = xx**2 + yy**2
          if (rsq >= (rad+rad2)**2) then
            zintercepts(1,ix,iy) = zcent + zmin
            zintercepts(2,ix,iy) = zcent + zmax
          else if (rad**2 <= rsq .and. rsq < (rad+rad2)**2) then
            rr = sqrt(rsq)
            ri = rad2 - sqrt(max(0.,rad2**2 - (rad+rad2-rr)**2))
            zintercepts(1,ix,iy) = zcent + zmin + ri
            zintercepts(2,ix,iy) = zcent + zmax - ri
          else
            zintercepts(1:2,ix,iy) = LARGEPOS
          endif
        enddo
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

      z2 = max(length*0.5-rad2,0.)

      do i=1,n

        rr = sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2)
        zz = abs(z(i) - zcent)
        if (zz >= z2 .and. rr >= rad+rad2) then
          distance(i) = zz - length*0.5
        else if (zz <= z2) then
          distance(i) = max(rad - rr,zz - length*0.5)
        else
          distance(i) = sqrt((zz - z2)**2 + (rr - (rad + rad2))**2) - rad2
        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZRoundedCylinderOutIntercept(rad,length,rad2,xcent,ycent,zcent,
     &                                        n,x,y,z,vx,vy,vz,xi,yi,zi,
     &                                        itheta,iphi)
      use Constant
      real(kind=8):: rad,length,rad2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)
      real(kind=8):: ri2,rr,zz,dr,theta

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i

c     --- first calculate for non-rounded cylinder and then readjust
c     --- the ones that are affected by rounded corners
      call ZCylinderOutIntercept(rad,length,xcent,ycent,zcent,
     &                           n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)


      do i=1,n
c        --- figure out which one are in the corners

c        --- between the rounded corners in Z => not affected
         if( zi(i)-zcent < 0.5*length-rad2 ) cycle

c        --- outside the rounded corners in R => not affected
         ri2 =  (xi(i)-xcent)**2 + (yi(i)-ycent)**2
         if( ri2 > (rad+rad2)**2 ) cycle


         rr = (rad+rad2) - sqrt(ri2)
         zz = zi(i)-zcent - (0.5*length-rad2)

c        --- calculate angles
         iphi(i)   = atan2(yi(i)-ycent,xi(i)-xcent)
         theta     = atan2(rr,zz)
         itheta(i) = atan2(rr,zz) + pi

c        --- adjust position
         dr    = rr - sin(theta)*rad2
         xi(i) = xi(i) - dr*cos(iphi(i))
         yi(i) = yi(i) - dr*sin(iphi(i))
         zi(i) = zi(i) - (zz - cos(theta)*rad2)

      enddo

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
          if (xx > 0.) then
            delmx(i) = -(sqrt(max(0.,rad**2 - yy**2 - zz**2)) + xx)
            delpx(i) = 0.
          else
            delpx(i) = -(sqrt(max(0.,rad**2 - yy**2 - zz**2)) - xx)
            delmx(i) = 0.
          endif
          if (yy > 0.) then
            delmy(i) = -(sqrt(max(0.,rad**2 - xx**2 - zz**2)) + yy)
            delpy(i) = 0.
          else
            delmy(i) = 0.
            delpy(i) = -(sqrt(max(0.,rad**2 - xx**2 - zz**2)) - yy)
          endif
          if (zz > 0.) then
            delmz(i) = -(sqrt(max(0.,rad**2 - xx**2 - yy**2)) + zz)
            delpz(i) = 0.
          else
            delmz(i) = 0.
            delpz(i) = -(sqrt(max(0.,rad**2 - xx**2 - yy**2)) - zz)
          endif
        endif

      enddo

      return
      end
c=======================================================================
      subroutine SphereConductorFnew(rad,xcent,ycent,zcent,
     &                               intercepts,fuzz)
      use ConductorInterceptTypeModule
      real(kind=8):: rad,xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      integer(ISZ):: ix,iy,iz
      real(kind=8):: xx,yy,zz,rsq

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2
      intercepts%nyicpt = 2
      intercepts%nzicpt = 2
      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent
          rsq = yy**2 + zz**2
          if (rsq <= rad**2) then
            xintercepts(1,iy,iz) = xcent - sqrt(rad**2 - rsq)
            xintercepts(2,iy,iz) = xcent + sqrt(rad**2 - rsq)
          else
            xintercepts(1,iy,iz) = +LARGEPOS
            xintercepts(2,iy,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rsq = xx**2 + zz**2
          if (rsq <= rad**2) then
            yintercepts(1,ix,iz) = ycent - sqrt(rad**2 - rsq)
            yintercepts(2,ix,iz) = ycent + sqrt(rad**2 - rsq)
          else
            yintercepts(1,ix,iz) = +LARGEPOS
            yintercepts(2,ix,iz) = +LARGEPOS
          endif
        enddo
      enddo

      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rsq = xx**2 + yy**2
          if (rsq <= rad**2) then
            zintercepts(1,ix,iy) = zcent - sqrt(rad**2 - rsq)
            zintercepts(2,ix,iy) = zcent + sqrt(rad**2 - rsq)
          else
            zintercepts(1,ix,iy) = +LARGEPOS
            zintercepts(2,ix,iy) = +LARGEPOS
          endif
        enddo
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
      subroutine SphereIntercept(rad,xcent,ycent,zcent,
     &                           n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      real(kind=8):: rad,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: xx,yy,zz,rr
      real(kind=8):: a,b,c,t1,t2,tt

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        a = vx(i)**2 + vy(i)**2 + vz(i)**2
        b = -2.*(vx(i)*xx + vy(i)*yy + vz(i)*zz)
        c = xx**2 + yy**2 + zz**2 - rad**2
        if ((b**2 - 4*a*c) < 0. .or. a == 0.) then
          xi(i) = LARGEPOS
          yi(i) = LARGEPOS
          zi(i) = LARGEPOS
          itheta(i) = 0.
          iphi(i) = 0.
        else
          t1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
          t2 = (-b - sqrt(b**2 - 4*a*c))/(2*a)
          if (t1 > 0. .and. t2 > 0.) then
            tt = min(t1,t2)
          else
            tt = max(t1,t2)
          endif
          xi(i) = -vx(i)*tt + x(i)
          yi(i) = -vy(i)*tt + y(i)
          zi(i) = -vz(i)*tt + z(i)

          rr = sqrt(xx**2 + yy**2)
          itheta(i) = atan2(rr,zz)
          iphi(i) = atan2(yy,xx)
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
      subroutine ztorus_intercept(r1,r2,thmin,thmax,xcent,ycent,zcent,
     &                            n,x,y,z,vx,vy,vz,tt,itheta,iphi)
      use Constant
      real(kind=8):: r1,r2,thmin,thmax,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: tt(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i,nhit,ii
      real(kind=8):: zmin,zmax
      real(kind=8):: vsq,d1,d2
      real(kind=8):: xx,yy,zz,rr,th,dd,rrsq,t0,ttemp(4),ttfinal
      real(kind=8):: a,b,c
      real(kind=8):: v(2),f(2)
      integer(ISZ):: niter,imax,imin

      zmin = zcent - r2
      zmax = zcent + r2

      do i=1,n
        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

c       --- This check on whether a point is inside may be removed at a
c       --- future time.
c       rr = sqrt((x(i) - xcent)**2 + (y(i) - ycent)**2)
c       dd = sqrt((z(i)-zcent)**2 + (r1 - rr)**2) - r2
c       if (dd > 0.) cycle

        nhit = 0

c       --- First, check if trajectory hits the planes just outside the torus.
        if (vz(i) .ne. 0.) then
          t0 = -(zmin - z(i))/vz(i)
          xx = -vx(i)*t0 + x(i)
          yy = -vy(i)*t0 + y(i)
          rrsq = (xx - xcent)**2 + (yy - ycent)**2
          if (max(0.,(r1-r2))**2 <= rrsq .and. rrsq <= (r1+r2)**2) then
            nhit = nhit + 1
            ttemp(nhit) = t0
          endif

          t0 = -(zmax - z(i))/vz(i)
          xx = -vx(i)*t0 + x(i)
          yy = -vy(i)*t0 + y(i)
          rrsq = (xx - xcent)**2 + (yy - ycent)**2
          if (max(0.,(r1-r2))**2 <= rrsq .and. rrsq <= (r1+r2)**2) then
            nhit = nhit + 1
            ttemp(nhit) = t0
          endif
        endif

c       --- Check if it hits the cylinder just outside the torus.
        a = vx(i)**2 + vy(i)**2
        if (a > 0.) then
          b = -2*(vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
          c = (x(i)-xcent)**2 + (y(i)-ycent)**2 - (r1+r2)**2
          if ((b**2 - 4*a*c) > 0.) then
            t0 = (-b-sqrt(b**2 - 4*a*c))/(2*a)
            zz = -vz(i)*t0 + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              nhit = nhit + 1
              ttemp(nhit) = t0
            endif
            t0 = (-b+sqrt(b**2 - 4*a*c))/(2*a)
            zz = -vz(i)*t0 + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              nhit = nhit + 1
              ttemp(nhit) = t0
            endif
          endif
        endif

c       --- If none of those surfaces is intersected, then no reason to
c       --- continue (this saves the work below).
        if (nhit == 0) cycle

c       --- Check if it hits the cylinder just inside the torus.
        a = vx(i)**2 + vy(i)**2
        if (a > 0.) then
          b = -2*(vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
          c = (x(i)-xcent)**2 + (y(i)-ycent)**2 - max(0.,(r1-r2))**2
          if ((b**2 - 4*a*c) > 0.) then
            t0 = (-b-sqrt(b**2 - 4*a*c))/(2*a)
            zz = -vz(i)*t0 + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              ttemp(nhit+1) = t0
              nhit = nhit + 1
            endif
            t0 = (-b+sqrt(b**2 - 4*a*c))/(2*a)
            zz = -vz(i)*t0 + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              ttemp(nhit+1) = t0
              nhit = nhit + 1
            endif
          endif
        endif

c       --- Check each of the hits and refine the ones where ttemp is positive.
        do ii=1,nhit

c         --- Do a simple iteration to try to refine the intersection point.
          v(1) = ttemp(ii)*0.6
          xx = -vx(i)*v(1) + x(i) - xcent
          yy = -vy(i)*v(1) + y(i) - ycent
          zz = -vz(i)*v(1) + z(i) - zcent
          rr = sqrt(xx**2 + yy**2)
          f(1) = r2**2 - (rr - r1)**2 - zz**2
          v(2) = ttemp(ii)
          xx = -vx(i)*v(2) + x(i) - xcent
          yy = -vy(i)*v(2) + y(i) - ycent
          zz = -vz(i)*v(2) + z(i) - zcent
          rr = sqrt(xx**2 + yy**2)
          f(2) = r2**2 - (rr - r1)**2 - zz**2
          niter = 0
          if (f(1) == 0. .or. f(2) == 0.) niter = 12
          imax = 1
          do while (min(abs(f(2)),abs(f(1))) > 1.e-12 .and. niter < 12)
            imax = 1
            if (abs(f(2)) > abs(f(1))) imax = 2
            if (f(1) == f(2)) exit
            v(imax) = (v(1)*f(2)-v(2)*f(1))/(f(2)-f(1))
            xx = -vx(i)*v(imax) + x(i) - xcent
            yy = -vy(i)*v(imax) + y(i) - ycent
            zz = -vz(i)*v(imax) + z(i) - zcent
            rr = sqrt(xx**2 + yy**2)
            f(imax) = r2**2 - (rr - r1)**2 - zz**2
            niter = niter + 1
          end do
          if (min(abs(f(2)),abs(f(1))) < 1.e-12) then
            imin = 1
            if (abs(f(2)) < abs(f(1))) imin = 2
            ttemp(ii) = v(imin)
          else
            ttemp(ii) = LARGEPOS
          endif

        enddo

c       --- Now find the minimum ttemp, only checking points with in
c       --- thmin and thmax.
        ttfinal = LARGEPOS
        do ii=1,nhit
          if ((ttfinal > ttemp(ii) .and. ttemp(ii) > 0. .and.
     &         ttfinal > 0.) .or.
     &        (ttemp(ii) > ttfinal .and. ttfinal < 0.)) then

            xx = -vx(i)*ttemp(ii) + x(i) - xcent
            yy = -vy(i)*ttemp(ii) + y(i) - ycent
            zz = -vz(i)*ttemp(ii) + z(i) - zcent
            rr = sqrt(xx**2 + yy**2)
            th = atan2(rr-r1,zz)
            if (th < thmin) th = th + 2*pi
            if (thmin <= th .and. th <= thmax) then
              imin = ii
              ttfinal = ttemp(ii)
            endif

          endif
        enddo

        if (ttfinal == LARGEPOS) cycle

        xx = -vx(i)*ttemp(imin) + x(i) - xcent
        yy = -vy(i)*ttemp(imin) + y(i) - ycent
        zz = -vz(i)*ttemp(imin) + z(i) - zcent
        rr = sqrt(xx**2 + yy**2)
        itheta(i) = atan2(rr-r1,zz)
        iphi(i)   = atan2(yy,xx)
        tt(i) = ttemp(imin)

      enddo

      return
      end
c=======================================================================
      subroutine ZTorusIntercept(r1,r2,xcent,ycent,zcent,
     &                           n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: r1,r2,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: tt(n),theta(n),phi(n)

      xi = LARGEPOS
      yi = LARGEPOS
      zi = LARGEPOS
      tt = LARGEPOS
      theta = 0.
      phi = 0.

      call ztorus_intercept(r1,r2,0.,2.*pi,xcent,ycent,zcent,
     &                      n,x,y,z,vx,vy,vz,tt,theta,phi)

      do i=1,n

        if (tt(i) == LARGEPOS) cycle

        xi(i) = -vx(i)*tt(i) + x(i)
        yi(i) = -vy(i)*tt(i) + y(i)
        zi(i) = -vz(i)*tt(i) + z(i)

        itheta(i) = theta(i)
        iphi(i) = phi(i)
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
c     --- In the limiting cases with the angles being exactly zero (or pi/2)
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
      lenhalf = length*0.5

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
c       xp = +xx*cos(theta)                        - zz*sin(theta)
c       yp = -xx*sin(theta)*sin(phi) + yy*cos(phi) - zz*cos(theta)*sin(phi)
c       zp = +xx*sin(theta)*cos(phi) + yy*sin(phi) + zz*cos(theta)*cos(phi)
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
        zz = abs(zp) - length*0.5

        dd = sqrt(max(0.,rr)**2 + max(0.,zz)**2)
        if (dd == 0.) dd = max(rr,zz)
        distance(i) = dd
      enddo

      return
      end
c=======================================================================
      subroutine CylinderIntercept(rad,length,theta,phi,xcent,ycent,zcent,
     &                             n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      real(kind=8):: rad,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      real(kind=8):: ttheta(n),tphi(n)
      real(kind=8):: xx(n),yy(n),zz(n)
      real(kind=8):: xp(n),yp(n),zp(n),vxp(n),vyp(n),vzp(n)
      real(kind=8):: xn(n),yn(n),zn(n),rn(n)

c     --- Translate and rotate to the Cylinder frame
      xx = x - xcent
      yy = y - ycent
      zz = z - zcent
      xp = +xx*cos(theta)                        - zz*sin(theta)
      yp = -xx*sin(theta)*sin(phi) + yy*cos(phi) - zz*cos(theta)*sin(phi)
      zp = +xx*sin(theta)*cos(phi) + yy*sin(phi) + zz*cos(theta)*cos(phi)
      vxp = +vx*cos(theta)                        - vz*sin(theta)
      vyp = -vx*sin(theta)*sin(phi) + vy*cos(phi) - vz*cos(theta)*sin(phi)
      vzp = +vx*sin(theta)*cos(phi) + vy*sin(phi) + vz*cos(theta)*cos(phi)

c     --- Get data
      call ZCylinderIntercept(rad,length,0.,0.,0.,
     &                        n,xp,yp,zp,vxp,vyp,vzp,xn,yn,zn,ttheta,tphi)

c     --- Rotate and translate interception points to lab frame
      xi = +xn*cos(theta) - yn*sin(theta)*sin(phi) + zn*sin(theta)*cos(phi)
      yi =                + yn*cos(phi)            + zn*sin(phi)
      zi = -xn*sin(theta) - yn*cos(theta)*sin(phi) + zn*cos(theta)*cos(phi)
      xi = xi + xcent
      yi = yi + ycent
      zi = zi + zcent

c     --- Get points on surface normal at the unit sphere
      xp = sin(ttheta)*cos(tphi)
      yp = sin(ttheta)*sin(tphi)
      zp = cos(ttheta)
c     --- Rotate to the lab frame
      xn = +xp*cos(theta) - yp*sin(theta)*sin(phi) + zp*sin(theta)*cos(phi)
      yn =                + yp*cos(phi)            + zp*sin(phi)
      zn = -xp*sin(theta) - yp*cos(theta)*sin(phi) + zp*cos(theta)*cos(phi)
c     --- Calculate angles in lab frame
      rn = sqrt(xn**2 + yn**2)
      itheta = atan2(rn,zn)
      iphi = atan2(yn,xn)

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
      real(kind=8):: zz,rmax
      integer(ISZ):: nt,ixmin,ixmax,iymin,iymax,tix,tiy,ii,id
      real(kind=8),allocatable:: xt(:),yt(:),zt(:)
      real(kind=8),allocatable:: mxt(:),pxt(:),myt(:),pyt(:),mzt(:),pzt(:)

      allocate(xt(n),yt(n),zt(n))
      allocate(mxt(n),pxt(n),myt(n),pyt(n),mzt(n),pzt(n))

c     --- Very kludgy code to extract number of grid points and cell
c     --- sizes from the input coordinates.
      ymin = minval(y)
      ymax = maxval(y)
      dy = y(2) - y(1)
      if (dy /= 0.) then
        ny = maxval(nint((y - ymin)/dy))
      else
        dy = 1.
        ny = 0
      endif

      dx = x(ny+2) - x(1)
      xmin = minval(x)
      xmax = maxval(x)
      nx = maxval(nint((x - xmin)/dx))

      delmx = LARGEPOS
      delpx = LARGEPOS
      delmy = LARGEPOS
      delpy = LARGEPOS
      delmz = LARGEPOS
      delpz = LARGEPOS

      do ic=1,ncylinders

c       --- First, get maximum extent of the cylinder
c       --- This won't be very efficient for a long skinny cylinder, but it is
c       --- gauranteed to work.
        rmax = sqrt(rad(ic)**2 + (length(ic)/2.)**2)
        zz = z(1) - zcent(ic)
        if (abs(zz) > rmax) cycle
        rmax = sqrt(rmax**2 - zz**2)
        ixmin = max(0 ,int((xcent(ic) - rmax - xmin)/dx) - 1)
        ixmax = min(nx,int((xcent(ic) + rmax - xmin)/dx) + 2)
        iymin = max(0 ,int((ycent(ic) - rmax - ymin)/dy) - 1)
        iymax = min(ny,int((ycent(ic) + rmax - ymin)/dy) + 2)

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

      deallocate(xt,yt,zt)
      deallocate(mxt,pxt,myt,pyt,mzt,pzt)

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
      subroutine CylindersIntercept(ncylinders,rad,length,theta,phi,
     &                              xcent,ycent,zcent,
     &                              n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      integer(ISZ):: ncylinders
      real(kind=8):: rad(ncylinders),length(ncylinders)
      real(kind=8):: theta(ncylinders),phi(ncylinders)
      real(kind=8):: xcent(ncylinders),ycent(ncylinders),zcent(ncylinders)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: ic,i
      real(kind=8):: txi(n),tyi(n),tzi(n),ttheta(n),tphi(n),tt(n),ttemp

      xi = LARGEPOS
      yi = LARGEPOS
      zi = LARGEPOS
      tt = LARGEPOS
      itheta = 0.
      iphi = 0.

      do ic=1,ncylinders
        call CylinderIntercept(rad(ic),length(ic),theta(ic),phi(ic),
     &                         xcent(ic),ycent(ic),zcent(ic),
     &                         n,x,y,z,vx,vy,vz,txi,tyi,tzi,ttheta,tphi)
        do i=1,n
          if (vx(i)**2 + vy(i)**2 + vz(i)**2 == 0.) cycle
          if (vx(i) .ne. 0.) then
            ttemp = -(txi(i) - x(i))/vx(i)
          else if (vy(i) .ne. 0.) then
            ttemp = -(tyi(i) - y(i))/vy(i)
          else if (vz(i) .ne. 0.) then
            ttemp = -(tzi(i) - z(i))/vz(i)
          endif
          if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &        (ttemp > tt(i) .and. tt(i) < 0.)) then
            tt(i) = ttemp
            xi(i) = txi(i)
            yi(i) = tyi(i)
            zi(i) = tzi(i)
            itheta(i) = ttheta(i)
            iphi(i) = tphi(i)
          endif
        enddo

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
      lenhalf = length*0.5

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

        d1 = sqrt((rp - r_zmin)**2 + (zp + length*0.5)**2)
        ang_zmin = atan2((rp - r_zmin),(zp + length*0.5)) -
     &         atan2((r_zmax - r_zmin),length)
        if (ang_zmin > +pi) ang_zmin = ang_zmin - 2*pi
        if (ang_zmin < -pi) ang_zmin = ang_zmin + 2*pi
        ang_zmax = atan2((r_zmin - r_zmax),-length) -
     &          atan2((rp - r_zmax),(zp - length*0.5))
        if (ang_zmax > +pi) ang_zmax = ang_zmax - 2*pi
        if (ang_zmax < -pi) ang_zmax = ang_zmax + 2*pi
        dd = d1*sin(ang_zmin)
        dz = abs(zp) - length*0.5

        if (zp <= -length*0.5 .and. rp <= r_zmin) then
          distance(i) = dz
        else if (ang_zmin >= pi*0.5) then
          distance(i) = sqrt((rp - r_zmin)**2 + (zp + length*0.5)**2)
        else if (zp >= length*0.5 .and. rp <= r_zmax) then
          distance(i) = dz
        else if (ang_zmax >= pi*0.5) then
          distance(i) = sqrt((rp - r_zmax)**2 + (zp - length*0.5)**2)
        else if (dd >= 0.) then
          distance(i) = dd
        else
          distance(i) = max(dd,dz)
        endif
      enddo

      return
      end
c=======================================================================
      subroutine zconeintercept(r_zmin,r_zmax,length,xcent,ycent,zcent,
     &                          n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      real(kind=8):: r_zmin,r_zmax,length,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i

      print*,"The intercept calculation is not yet implemented for Cone"
      do i=1,n
        xi(i) = LARGEPOS
        yi(i) = LARGEPOS
        zi(i) = LARGEPOS
        itheta(i) = 0.
        iphi(i) = 0.
      enddo

      return
      end
c=======================================================================
      subroutine ConeIntercept(r_zmin,r_zmax,length,theta,phi,
     &                         xcent,ycent,zcent,
     &                         n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: r_zmin,r_zmax,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      real(kind=8):: ttheta(n),tphi(n)
      real(kind=8):: xx(n),yy(n),zz(n)
      real(kind=8):: xp(n),yp(n),zp(n),vxp(n),vyp(n),vzp(n)
      real(kind=8):: xn(n),yn(n),zn(n),rn(n)

c     --- Translate and rotate to frame of Cone
      xx = x - xcent
      yy = y - ycent
      zz = z - zcent
      xp = +xx*cos(theta)                        - zz*sin(theta)
      yp = -xx*sin(theta)*sin(phi) + yy*cos(phi) - zz*cos(theta)*sin(phi)
      zp = +xx*sin(theta)*cos(phi) + yy*sin(phi) + zz*cos(theta)*cos(phi)
      vxp = +vx*cos(theta)                        - vz*sin(theta)
      vyp = -vx*sin(theta)*sin(phi) + vy*cos(phi) - vz*cos(theta)*sin(phi)
      vzp = +vx*sin(theta)*cos(phi) + vy*sin(phi) + vz*cos(theta)*cos(phi)

c     --- Get data
      call zconeintercept(r_zmin,r_zmax,length,0.,0.,0.,
     &                    n,x,y,z,vx,vy,vz,xi,yi,zi,ttheta,tphi)

c     --- Rotate and translate interception points to lab frame
      xi = +xn*cos(theta) - yn*sin(theta)*sin(phi) + zn*sin(theta)*cos(phi)
      yi =                + yn*cos(phi)            + zn*sin(phi)
      zi = -xn*sin(theta) - yn*cos(theta)*sin(phi) + zn*cos(theta)*cos(phi)
      xi = xi + xcent
      yi = yi + ycent
      zi = zi + zcent

c     --- Find point along surface normal at the unit sphere
      xp = sin(ttheta)*cos(tphi)
      yp = sin(ttheta)*sin(tphi)
      zp = cos(ttheta)
c     --- Convert to lab frame coordinates
      xn = +xp*cos(theta) - yp*sin(theta)*sin(phi) + zp*sin(theta)*cos(phi)
      yn =                + yp*cos(phi)            + zp*sin(phi)
      zn = -xp*sin(theta) - yp*cos(theta)*sin(phi) + zp*cos(theta)*cos(phi)
c     --- Generate angles on lab frame
      rn = sqrt(xn**2 + yn**2)
      itheta = atan2(rn,zn)
      iphi = atan2(yn,xn)

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
      real(kind=8):: zz,rmax
      integer(ISZ):: nt,ixmin,ixmax,iymin,iymax,tix,tiy,ii,id
      real(kind=8),allocatable:: xt(:),yt(:),zt(:)
      real(kind=8),allocatable:: mxt(:),pxt(:),myt(:),pyt(:),mzt(:),pzt(:)

      allocate(xt(n),yt(n),zt(n))
      allocate(mxt(n),pxt(n),myt(n),pyt(n),mzt(n),pzt(n))

c     --- Very kludgy code to extract number of grid points and cell
c     --- sizes from the input coordinates.
      ymin = minval(y)
      ymax = maxval(y)
      dy = y(2) - y(1)
      if (dy /= 0.) then
        ny = maxval(nint((y - ymin)/dy))
      else
        dy = 1.
        ny = 0
      endif

      xmin = minval(x)
      xmax = maxval(x)
      dx = x(ny+2) - x(1)
      nx = maxval(nint((x - xmin)/dx))

      delmx = LARGEPOS
      delpx = LARGEPOS
      delmy = LARGEPOS
      delpy = LARGEPOS
      delmz = LARGEPOS
      delpz = LARGEPOS

      do ic=1,ncones

c       --- First, get maximum extent of cone
c       --- This won't be very efficient for a long skinny cone, but it is
c       --- gauranteed to work.
        rmax = sqrt(max(abs(r_zmin(ic)),abs(r_zmax(ic)))**2 +
     &                  (length(ic)/2.)**2)
        zz = z(1) - zcent(ic)
        if (abs(zz) > rmax) cycle
        rmax = sqrt(rmax**2 - zz**2)
        ixmin = max(0 ,int((xcent(ic) - rmax - xmin)/dx) - 1)
        ixmax = min(nx,int((xcent(ic) + rmax - xmin)/dx) + 2)
        iymin = max(0 ,int((ycent(ic) - rmax - ymin)/dy) - 1)
        iymax = min(ny,int((ycent(ic) + rmax - ymin)/dy) + 2)


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

      deallocate(xt,yt,zt)
      deallocate(mxt,pxt,myt,pyt,mzt,pzt)

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
      subroutine ConesIntercept(ncones,r_zmin,r_zmax,length,theta,phi,
     &                          xcent,ycent,zcent,
     &                          n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      integer(ISZ):: ncones
      real(kind=8):: r_zmin(ncones),r_zmax(ncones),length(ncones)
      real(kind=8):: theta(ncones),phi(ncones)
      real(kind=8):: xcent(ncones),ycent(ncones),zcent(ncones)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: ic,i
      real(kind=8):: txi(n),tyi(n),tzi(n),ttheta(n),tphi(n),tt(n),ttemp

      xi = LARGEPOS
      yi = LARGEPOS
      zi = LARGEPOS
      tt = LARGEPOS
      itheta = 0.
      iphi = 0.

      do ic=1,ncones
        call ConeIntercept(r_zmin(ic),r_zmax(ic),length(ic),theta(ic),phi(ic),
     &                     xcent(ic),ycent(ic),zcent(ic),
     &                     n,x,y,z,vx,vy,vz,txi,tyi,tzi,ttheta,tphi)
        do i=1,n
          if (vx(i)**2 + vy(i)**2 + vz(i)**2 == 0.) cycle
          if (vx(i) .ne. 0.) then
            ttemp = -(txi(i) - x(i))/vx(i)
          else if (vy(i) .ne. 0.) then
            ttemp = -(tyi(i) - y(i))/vy(i)
          else if (vz(i) .ne. 0.) then
            ttemp = -(tzi(i) - z(i))/vz(i)
          endif
          if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &        (ttemp > tt(i) .and. tt(i) < 0.)) then
            tt(i) = ttemp
            xi(i) = txi(i)
            yi(i) = tyi(i)
            zi(i) = tzi(i)
            itheta(i) = ttheta(i)
            iphi(i) = tphi(i)
          endif
        enddo

      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine alignedannulus(rmin,rmax,zmin,zmax,xx,yy,zz,
     &                          delmx,delpx,delmy,delpy,delmz,delpz,fuzz)
      real(kind=8):: rmin,rmax,zmin,zmax
      real(kind=8):: xx,yy,zz
      real(kind=8):: delmx,delpx,delmy,delpy,delmz,delpz,fuzz

      real(kind=8):: rsq

      rsq = xx**2 + yy**2

c     --- If point is not within the z extent of the cylinder
      if (zz < zmin - fuzz) then
        if (rmin**2 <= rsq .and. rsq <= rmax**2) delpz = zmin - zz
        return
      endif
      if (zz > zmax + fuzz) then
        if (rmin**2 <= rsq .and. rsq <= rmax**2) delmz = zz - zmax
        return
      endif

c     --- Point inside the inner cylinder
      if (rsq < rmin**2) then
        delmx = (sqrt(rmin**2 - yy**2) + xx)
        delpx = (sqrt(rmin**2 - yy**2) - xx)
        delmy = (sqrt(rmin**2 - xx**2) + yy)
        delpy = (sqrt(rmin**2 - xx**2) - yy)
        return
      endif

c     --- Point is outside the outer cylinder
      if (rsq > rmax**2) then
        if (abs(yy) <= rmax) then
          if (xx > 0) delmx = ( xx - sqrt(max(0.,rmax**2 - yy**2)))
          if (xx < 0) delpx = (-xx - sqrt(max(0.,rmax**2 - yy**2)))
        endif
        if (abs(xx) <= rmax) then
          if (yy > 0) delmy = ( yy - sqrt(max(0.,rmax**2 - xx**2)))
          if (yy < 0) delpy = (-yy - sqrt(max(0.,rmax**2 - xx**2)))
        endif
        return
      endif

c     --- Point is between the two cylinders
      if (rmin**2 < rsq .and. rsq < rmax**2) then
        delmz = min(0.,zmin - zz)
        delpz = min(0.,zz - zmax)
        delmx = -(sqrt(rmax**2 - yy**2) + xx)
        delpx = -(sqrt(rmax**2 - yy**2) - xx)
        delmy = -(sqrt(rmax**2 - xx**2) + yy)
        delpy = -(sqrt(rmax**2 - xx**2) - yy)
        if (abs(yy) <= rmin) then
          if (xx > 0) delmx = -( xx - sqrt(max(0.,rmin**2 - yy**2)))
          if (xx < 0) delpx = -(-xx - sqrt(max(0.,rmin**2 - yy**2)))
        endif
        if (abs(xx) <= rmin) then
          if (yy > 0) delmy = -( yy - sqrt(max(0.,rmin**2 - xx**2)))
          if (yy < 0) delpy = -(-yy - sqrt(max(0.,rmin**2 - xx**2)))
        endif
        return
      endif

c     --- The point is on the inner cylinder exactly.
      if (rsq < rmax**2) then
        if (xx > 0.) then
          delmx = -(sqrt(max(0.,rmin**2 - yy**2)) + xx)
          delpx = 0.
        else
          delmx = 0.
          delpx = -(sqrt(max(0.,rmin**2 - yy**2)) - xx)
        endif
        if (yy > 0.) then
          delmy = -(sqrt(max(0.,rmin**2 - xx**2)) + yy)
          delpy = 0.
        else
          delmy = 0.
          delpy = -(sqrt(max(0.,rmin**2 - xx**2)) - yy)
        endif
        delmz = min(0.,zmin - zz)
        delpz = min(0.,zz - zmax)
      else if (rsq > rmin**2) then
        if (xx > 0.) then
          delmx = -(sqrt(max(0.,rmax**2 - yy**2)) + xx)
          delpx = 0.
        else
          delmx = 0.
          delpx = -(sqrt(max(0.,rmax**2 - yy**2)) - xx)
        endif
        if (yy > 0.) then
          delmy = -(sqrt(max(0.,rmax**2 - xx**2)) + yy)
          delpy = 0.
        else
          delmy = 0.
          delpy = -(sqrt(max(0.,rmax**2 - xx**2)) - yy)
        endif
        delmz = min(0.,zmin - zz)
        delpz = min(0.,zz - zmax)
      endif

      return
      end
c=======================================================================
      subroutine nonalignedannulus(zlout,zlin,zuin,zuout,zlend,zuend,zz,
     &                             delmz,delpz)
      real(kind=8):: zlout,zlin,zuin,zuout,zlend,zuend
      real(kind=8):: zz
      real(kind=8):: delmz,delpz

c     --- Only deal with points where the z-line hits at least the outer
c     --- cylinder
      if (zlout <= zuout) then

c       --- Deal with points where the z-line intersects the both the outer
c       --- and inner cylinders.
        if (zlin <= zuin) then
          if (zz < zlout) then
            if (zlend <= zlout .and. zlout <= zuend) then
              delpz = zlout - zz
            else if (zlend <= zlin .and. zlin <= zuend) then
              delpz = zlend - zz
            else if (zlend <= zuin .and. zuin <= zuend) then
              delpz = zuin - zz
            endif
          else if (zlout <= zz .and. zz <= zlin) then
            if (zlend <= zz .and. zz <= zuend) then
              delmz = max(zlend,zlout) - zz
              delpz = zz - min(zuend,zlin)
            else if (zz < zlend) then
              if (zlend <= zlin) then
                delpz = zlend - zz
              else if (zlend <= zuin .and. zuin <= zuend) then
                delpz = zuin - zz
              else if (zuin <= zlend .and. zlend <= zuout) then
                delpz = zlend - zz
              endif
            else if (zz > zuend) then
              if (zuend >= zlout) delmz = zz - zuend
            endif
          else if (zlin < zz .and. zz < zuin) then
            if (zlend <= zlin .and. zlin <= zuend) delmz = zz - zlin
            if (zuend >= zuin .and. zuin >= zlend) delpz = zuin - zz
          else if (zuin <= zz .and. zz <= zuout) then
            if (zlend <= zz .and. zz <= zuend) then
              delmz = max(zlend,zuin) - zz
              delpz = zz - min(zuend,zuout)
            else if (zz > zuend) then
              if (zuin <= zuend) then
                delmz = zz - zuend
              else if (zlend <= zlin .and. zlin <= zuend) then
                delmz = zz - zlin
              else if (zlout <= zuend .and. zuend <= zlin) then
                delmz = zz - zuend
              endif
            else if (zz < zlend) then
              if (zlend <= zuout) delpz = zlend - zz
            endif
          else if (zz > zuout) then
            if (zuend >= zuout .and. zuout >= zlend) then
              delmz = zz - zuout
            else if (zuend >= zuin .and. zuin >= zlend) then
              delmz = zz - zuend
            else if (zuend >= zlin .and. zlin >= zlend) then
              delmz = zz - zlin
            endif
          endif

        else
c       --- Deal with point where the z-line only intersects the outer
c       --- cylinder
          if (zz < zlout) then
            if (zlend <= zlout .and. zlout <= zuend) delpz = zlout - zz
          else if (zlout <= zz .and. zz <= zuout) then
            if (zlend <= zz .and. zz <= zuend) then
              delmz = max(zlend,zlout) - zz
              delpz = zz - min(zuend,zuout)
            else if (zz < zlend) then
              if (zlend <= zuout) delpz = zlend - zz
            else if (zz > zuend) then
              if (zuend >= zlout) delmz = zz - zuend
            endif
          else if (zz > zuout) then
            if (zuend >= zuout .and. zuout >= zlend) delmz = zz - zuout
          endif
        endif
              
      endif

      return
      end
c=======================================================================
      subroutine AnnulusConductorF(rmin,rmax,length,theta,phi,xcent,ycent,zcent,
     &                             n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                             fuzz)
      real(kind=8):: rmin,rmax,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a annulus at an arbitrary
c angle.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: ctheta,stheta,cphi,sphi,ss,cci,sci,si,lenhalf
      real(kind=8):: xx,yy,zz
      real(kind=8):: a,b,cmin,cmax
      real(kind=8):: z1,z2,zlin,zuin,zlout,zuout,zlend,zuend,zi1,zi2
      real(kind=8):: x1,x2,xlin,xuin,xlout,xuout,xlend,xuend,xi1,xi2
      real(kind=8):: y1,y2,ylin,yuin,ylout,yuout,ylend,yuend,yi1,yi2

      ctheta = cos(theta)
      stheta = sin(theta)
      cphi = cos(phi)
      sphi = sin(phi)
c     --- In the limiting cases with the angles being exactly zero (or pi/2)
c     --- the algorithm breaks down. The sign of expressions dependent on
c     --- the angles is important, and is lost when the angles are exactly
c     --- zero (or pi/2). So the limiting cases is used where the angle
c     --- approaches but does not equal to zero (or pi/2).
c     if (ctheta == 0) ctheta = SMALLPOS
c     if (stheta == 0) stheta = SMALLPOS
c     if (cphi == 0) cphi = SMALLPOS
c     if (sphi == 0) sphi = SMALLPOS
      if (ctheta*cphi .ne. 0.) cci = 1./(ctheta*cphi)
      if (stheta*cphi .ne. 0.) sci = 1./(stheta*cphi)
      if (sphi .ne. 0.) si  = 1./(sphi)
      lenhalf = length*0.5

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

c       --- First, handle the special cases, where the annulus is aligned
c       --- with one of the grid axes.
        if (stheta**2 + (ctheta*sphi)**2 == 0.) then
c         --- Along the z-axis
          call alignedannulus(rmin,rmax,-lenhalf,+lenhalf,xx,yy,zz,
     &                        delmx(i),delpx(i),
     &                        delmy(i),delpy(i),
     &                        delmz(i),delpz(i),fuzz)
          cycle
        endif

        if ((ctheta)**2 + (stheta*sphi)**2 == 0.) then
c         --- Along the x-axis
          call alignedannulus(rmin,rmax,-lenhalf,+lenhalf,yy,zz,xx,
     &                        delmy(i),delpy(i),
     &                        delmz(i),delpz(i),
     &                        delmx(i),delpx(i),fuzz)
          cycle
        endif

        if ((cphi)**2 == 0.) then
c         --- Along the y-axis
          call alignedannulus(rmin,rmax,-lenhalf,+lenhalf,zz,xx,yy,
     &                        delmz(i),delpz(i),
     &                        delmx(i),delpx(i),
     &                        delmy(i),delpy(i),fuzz)
          cycle
        endif

c       --- delta z
c       --- Coefficients of the quadratic
        a = stheta**2 + (ctheta*sphi)**2
        b = -2*xx*ctheta*stheta - 2*(yy*cphi - xx*stheta*sphi)*ctheta*sphi
        cmin = (xx*ctheta)**2 + (yy*cphi - xx*stheta*sphi)**2 - rmin**2
        cmax = (xx*ctheta)**2 + (yy*cphi - xx*stheta*sphi)**2 - rmax**2

c       --- Get the point where the z-line intersects the end planes
        if (ctheta*cphi .ne. 0.) then
          z1 = (-lenhalf - xx*stheta*cphi - yy*sphi)*cci
          z2 = (+lenhalf - xx*stheta*cphi - yy*sphi)*cci
          zlend = min(z1,z2)
          zuend = max(z1,z2)
        else
          zlend = -LARGEPOS
          zuend = +LARGEPOS
c         --- The x ends need to be calculated to check if the point is
c         --- between the end planes
          x1 = (-lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          x2 = (+lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          xlend = min(x1,x2)
          xuend = max(x1,x2)
        endif

c       --- Check if z-line intersects the inner cylinder
        if (b**2 - 4*a*cmin >= 0) then
c         --- Get the two points of intersection
          ss = sqrt(b**2 - 4*a*cmin)
          z1 = (-b-ss)/(2*a)
          z2 = (-b+ss)/(2*a)
          zlin = min(z1,z2)
          zuin = max(z1,z2)
        else
          zlin = +LARGEPOS
          zuin = -LARGEPOS
        endif
c       --- Check if z-line intersects the outer cylinder
        if (b**2 - 4*a*cmax >= 0) then
c         --- Get the two points of intersection
          ss = sqrt(b**2 - 4*a*cmax)
          z1 = (-b-ss)/(2*a)
          z2 = (-b+ss)/(2*a)
          zlout = min(z1,z2)
          zuout = max(z1,z2)
        else
          zlout = +LARGEPOS
          zuout = -LARGEPOS
        endif

        if (ctheta*cphi .ne. 0. .or. (xlend <= xx .and. xx < xuend)) then
          call nonalignedannulus(zlout,zlin,zuin,zuout,zlend,zuend,zz,
     &                           delmz(i),delpz(i))
        endif

c       --- delta x
        a = (ctheta)**2 + (stheta*sphi)**2
        b = -2*zz*stheta*ctheta - 2*(yy*cphi - zz*ctheta*sphi)*stheta*sphi
        cmin = (zz*stheta)**2 + (yy*cphi - zz*ctheta*sphi)**2 - rmin**2
        cmax = (zz*stheta)**2 + (yy*cphi - zz*ctheta*sphi)**2 - rmax**2

c       --- Get the point where the x-line intersects the end planes
        if (stheta*cphi .ne. 0.) then
          x1 = (-lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          x2 = (+lenhalf - zz*ctheta*cphi - yy*sphi)*sci
          xlend = min(x1,x2)
          xuend = max(x1,x2)
        else
          xlend = -LARGEPOS
          xuend = +LARGEPOS
        endif

c       --- Check if x-line intersects the inner cylinder
        if (b**2 - 4*a*cmin >= 0) then
c         --- Get the two points of intersection
          ss = sqrt(b**2 - 4*a*cmin)
          x1 = (-b-ss)/(2*a)
          x2 = (-b+ss)/(2*a)
          xlin = min(x1,x2)
          xuin = max(x1,x2)
        else
          xlin = +LARGEPOS
          xuin = -LARGEPOS
        endif
c       --- Check if x-line intersects the outer cylinder
        if (b**2 - 4*a*cmax >= 0) then
c         --- Get the two points of intersection
          ss = sqrt(b**2 - 4*a*cmax)
          x1 = (-b-ss)/(2*a)
          x2 = (-b+ss)/(2*a)
          xlout = min(x1,x2)
          xuout = max(x1,x2)
        else
          xlout = +LARGEPOS
          xuout = -LARGEPOS
        endif

        if (stheta*cphi .ne. 0. .or. (zlend <= zz .and. zz <= zuend)) then
          call nonalignedannulus(xlout,xlin,xuin,xuout,xlend,xuend,xx,
     &                           delmx(i),delpx(i))
        endif

c       --- delta y
        a = (cphi)**2
        b = -2*(zz*ctheta*sphi + xx*stheta*sphi)*cphi
        cmin = (xx*ctheta - zz*stheta)**2 + (zz*ctheta*sphi + xx*stheta*sphi)**2
     &      -rmin**2
        cmax = (xx*ctheta - zz*stheta)**2 + (zz*ctheta*sphi + xx*stheta*sphi)**2
     &      -rmax**2

c       --- Get the point where the y-line intersects the end planes
        if (sphi .ne. 0.) then
          y1 = (-lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          y2 = (+lenhalf - zz*ctheta*cphi - xx*stheta*cphi)*si
          ylend = min(y1,y2)
          yuend = max(y1,y2)
        else
          ylend = -LARGEPOS
          yuend = +LARGEPOS
        endif

c       --- Check if y-line intersects the inner cylinder
        if (b**2 - 4*a*cmin >= 0) then
c         --- Get the two points of intersection
          ss = sqrt(b**2 - 4*a*cmin)
          y1 = (-b-ss)/(2*a)
          y2 = (-b+ss)/(2*a)
          ylin = min(y1,y2)
          yuin = max(y1,y2)
        else
          ylin = +LARGEPOS
          yuin = -LARGEPOS
        endif
c       --- Check if y-line intersects the outer cylinder
        if (b**2 - 4*a*cmax >= 0) then
c         --- Get the two points of intersection
          ss = sqrt(b**2 - 4*a*cmax)
          y1 = (-b-ss)/(2*a)
          y2 = (-b+ss)/(2*a)
          ylout = min(y1,y2)
          yuout = max(y1,y2)
        else
          ylout = +LARGEPOS
          yuout = -LARGEPOS
        endif

        if (sphi .ne. 0. .or. (zlend <= zz .and. zz <= zuend)) then
          call nonalignedannulus(ylout,ylin,yuin,yuout,ylend,yuend,yy,
     &                           delmy(i),delpy(i))
        endif

      enddo

      return
      end
c=======================================================================
      subroutine AnnulusConductorD(rmin,rmax,length,theta,phi,
     &                             xcent,ycent,zcent,
     &                             n,x,y,z,distance)
      real(kind=8):: rmin,rmax,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances from a annulus at an arbitrary angle.
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: ctheta,stheta,cphi,sphi
      real(kind=8):: xx,yy,zz,xp,yp,zp,rrin,rrout,dd

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

        rrin  = rmin - sqrt(xp**2 + yp**2)
        rrout = sqrt(xp**2 + yp**2) - rmax
        zz = abs(zp) - length*0.5

        if (rrout > 0.) then
c         --- Outside the outer cylinder
          dd = sqrt(rrout**2 + max(0.,zz)**2)
        else if (rrin > 0.) then
c         --- Inside the inner cylinder
          dd = sqrt(rrin**2 + max(0.,zz)**2)
        else if (zz > 0.) then
c         --- Within the annular radii, but beyond the z extent
          dd = zz
        else
c         --- Finally, within the annulus
          dd = max(max(rrin,rrout),zz)
        endif
        distance(i) = dd

      enddo

      return
      end
c=======================================================================
      subroutine AnnulusIntercept(rmin,rmax,length,theta,phi,xcent,ycent,zcent,
     &                            n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      real(kind=8):: rmin,rmax,length,theta,phi,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: ttheta(n),tphi(n)
      real(kind=8):: tthetaout(n),tphiout(n)
      real(kind=8):: xx(n),yy(n),zz(n)
      real(kind=8):: xp(n),yp(n),zp(n),vxp(n),vyp(n),vzp(n)
      real(kind=8):: xn(n),yn(n),zn(n),rn(n)
      real(kind=8):: xnout(n),ynout(n),znout(n)
      real(kind=8):: vsq,tin,tout

c     --- Translate and rotate to the Annulus frame
      xx = x - xcent
      yy = y - ycent
      zz = z - zcent
      xp = +xx*cos(theta)                        - zz*sin(theta)
      yp = -xx*sin(theta)*sin(phi) + yy*cos(phi) - zz*cos(theta)*sin(phi)
      zp = +xx*sin(theta)*cos(phi) + yy*sin(phi) + zz*cos(theta)*cos(phi)
      vxp = +vx*cos(theta)                        - vz*sin(theta)
      vyp = -vx*sin(theta)*sin(phi) + vy*cos(phi) - vz*cos(theta)*sin(phi)
      vzp = +vx*sin(theta)*cos(phi) + vy*sin(phi) + vz*cos(theta)*cos(phi)

c     --- Get data for both inner and out cylinders
      call ZCylinderIntercept(rmin,length,0.,0.,0.,
     &                        n,xp,yp,zp,vxp,vyp,vzp,xn,yn,zn,ttheta,tphi)
      call ZCylinderIntercept(rmax,length,0.,0.,0.,
     &                        n,xp,yp,zp,vxp,vyp,vzp,xnout,ynout,znout,
     &                        tthetaout,tphiout)

c     --- Find the most recent intercept time
      do i=1,n
        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

c       --- Find the intercept times to the inner and outer cylinders.
c       --- It shouldn't matter which direction is used, as long as the
c       --- velocity is not zero.
c       --- Note that t is time since the intercept.
        if (vxp(i) .ne. 0.) then
          tin = (xp(i) - xn(i))/vxp(i)
          tout = (xp(i) - xnout(i))/vxp(i)
        else if (vyp(i) .ne. 0.) then
          tin = (yp(i) - yn(i))/vyp(i)
          tout = (yp(i) - ynout(i))/vyp(i)
        else if (vzp(i) .ne. 0.) then
          tin = (zp(i) - zn(i))/vzp(i)
          tout = (zp(i) - znout(i))/vzp(i)
        endif

c       --- The output is placed in the xn etc arrays. If tout is the
c       --- smaller positive intercept time, then swap that data into
c       --- those output arrays.
        if ((tin > 0. .and. tout > 0. .and. tout < tin) .or.
     &      (tin < 0. .and. tout > 0.)) then
          xn(i) = xnout(i)
          yn(i) = ynout(i)
          zn(i) = znout(i)
          ttheta(i) = tthetaout(i)
          tphi(i) = tphiout(i)
        endif

      enddo

c     --- Rotate and translate interception points to lab frame
      xi = +xn*cos(theta) - yn*sin(theta)*sin(phi) + zn*sin(theta)*cos(phi)
      yi =                + yn*cos(phi)            + zn*sin(phi)
      zi = -xn*sin(theta) - yn*cos(theta)*sin(phi) + zn*cos(theta)*cos(phi)
      xi = xi + xcent
      yi = yi + ycent
      zi = zi + zcent

c     --- Get points on surface normal at the unit sphere
      xp = sin(ttheta)*cos(tphi)
      yp = sin(ttheta)*sin(tphi)
      zp = cos(ttheta)
c     --- Rotate to the lab frame
      xn = +xp*cos(theta) - yp*sin(theta)*sin(phi) + zp*sin(theta)*cos(phi)
      yn =                + yp*cos(phi)            + zp*sin(phi)
      zn = -xp*sin(theta) - yp*cos(theta)*sin(phi) + zp*cos(theta)*cos(phi)
c     --- Calculate angles in lab frame
      rn = sqrt(xn**2 + yn**2)
      itheta = atan2(rn,zn)
      iphi = atan2(yn,xn)

      return
      end
c=======================================================================
c=======================================================================
      subroutine ZGridConductorF(xcellsize,ycellsize,length,thickness,
     &                           xcent,ycent,zcent,
     &                           n,x,y,z,delmx,delpx,delmy,delpy,delmz,delpz,
     &                           fuzz)
      real(kind=8):: xcellsize,ycellsize,length,thickness
      real(kind=8):: xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

c Finds distances along the axis of the points from a grid
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i
      real(kind=8):: zmin,zmax

      zmin = zcent - length/2.
      zmax = zcent + length/2.

c     --- This cheats for now, putting in a solid conductor within the
c     --- z extent.
      do i=1,n

        delmx(i) = LARGEPOS
        delpx(i) = LARGEPOS
        delmy(i) = LARGEPOS
        delpy(i) = LARGEPOS
        delmz(i) = LARGEPOS
        delpz(i) = LARGEPOS

        if (z(i) <= zmin) then
          delpz(i) = zmin - z(i)
        else if (z(i) <= zmax) then
          delpz(i) = z(i) - zmax
        endif
        if (z(i) >= zmax) then
          delmz(i) = z(i) - zmax
        else if (z(i) >= zmin) then
          delmz(i) = zmin - z(i)
        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZGridConductorD(xcellsize,ycellsize,length,thickness,
     &                           xcent,ycent,zcent,
     &                           n,x,y,z,distance)
      real(kind=8):: xcellsize,ycellsize,length,thickness
      real(kind=8):: xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances to the grid
c The notation is that a negative distance means that the point is inside.

      integer(ISZ):: i,ix,iy
      real(kind=8):: zmin,zmax
      real(kind=8):: dx,dy,dxl,dxu,dyl,dyu
      real(kind=8):: dxmin,dymin

      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

        if (z(i) < zmin) then
          distance(i) = zmin - z(i)
        else if (z(i) > zmax) then
          distance(i) = z(i) - zmax
        else
          ix = int(abs(x(i) - xcent)/xcellsize)
          dx = abs(x(i) - xcent) - ix*xcellsize
          iy = int(abs(y(i) - ycent)/ycellsize)
          dy = abs(y(i) - ycent) - iy*ycellsize

          if (dx <= thickness/2.) then
            dxl = dx - thickness/2.
            dxu = LARGEPOS
          else if (dx >= xcellsize - thickness/2.) then
            dxl = LARGEPOS
            dxu = xcellsize - thickness/2. - dx
          else
            dxl = dx - thickness/2.
            dxu = xcellsize - thickness/2. - dx
          endif

          if (dy <= thickness/2.) then
            dyl = dy - thickness/2.
            dyu = LARGEPOS
          else if (dy >= ycellsize - thickness/2.) then
            dyl = LARGEPOS
            dyu = ycellsize - thickness/2. - dy
          else
            dyl = dy - thickness/2.
            dyu = ycellsize - thickness/2. - dy
          endif

          dxmin = min(dxl,dxu)
          dymin = min(dyl,dyu)

          if (dxmin < 0. .and. dymin < 0.) then
            distance(i) = -sqrt(dxmin**2 + dymin**2)
          else
            distance(i) = min(dxmin,dymin)
          endif

        endif

      enddo

      return
      end
c=======================================================================
      subroutine ZGridIntercept(xcellsize,ycellsize,length,thickness,
     &                          xcent,ycent,zcent,
     &                          n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      real(kind=8):: xcellsize,ycellsize,length,thickness
      real(kind=8):: xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i,ix,iy,signx,signy
      real(kind=8):: zmin,zmax
      real(kind=8):: dx,dy,dt
      real(kind=8):: dtx,dty,dtz
      real(kind=8):: ithetax,iphix
      real(kind=8):: ithetay,iphiy
      real(kind=8):: xiz,yiz,ziz,ithetaz,iphiz
      real(kind=8):: dd
      logical(ISZ):: linside

      zmin = zcent - length/2.
      zmax = zcent + length/2.

      do i=1,n

        if (x(i) < xcent) then
          signx = -1
        else
          signx = +1
        endif

        if (y(i) < ycent) then
          signy = -1
        else
          signy = +1
        endif

        ix = int(signx*(x(i) - xcent)/xcellsize)
        dx = signx*(x(i) - xcent) - ix*xcellsize
        iy = int(signy*(y(i) - ycent)/ycellsize)
        dy = signy*(y(i) - ycent) - iy*ycellsize

        if (signx*vx(i) < 0.) then
          if (dx < thickness/2.) then
            dtx = (thickness/2. - dx)/(signx*vx(i))
          else
            dtx = (xcellsize + thickness/2. - dx)/(signx*vx(i))
          endif
          ithetax = pi*0.5
          iphix = 0.
        else if (signx*vx(i) > 0.) then
          if (dx < xcellsize - thickness/2.) then
            dtx = (-thickness/2. - dx)/(signx*vx(i))
          else
            dtx = (xcellsize - thickness/2. - dx)/(signx*vx(i))
          endif
          ithetax = pi*0.5
          iphix = pi
        else
          dtx = -LARGEPOS
          ithetax = 0.
          iphix = 0.
        endif

        if (signy*vy(i) < 0.) then
          if (dy < thickness/2.) then
            dty = (thickness/2. - dy)/(signy*vy(i))
          else
            dty = (ycellsize + thickness/2. - dy)/(signy*vy(i))
          endif
          ithetay = pi*0.5
          iphiy = 0.
        else if (signy*vy(i) > 0.) then
          if (dy < ycellsize - thickness/2.) then
            dty = (-thickness/2. - dy)/(signy*vy(i))
          else
            dty = (ycellsize - thickness/2. - dy)/(signy*vy(i))
          endif
          ithetay = pi*0.5
          iphiy = pi
        else
          dty = -LARGEPOS
          ithetay = 0.
          iphiy = 0.
        endif

        if (vz(i) < 0.) then
          dtz = (zmax - z(i))/vz(i)
          xiz = vx(i)*dtz + x(i)
          yiz = vy(i)*dtz + y(i)
          ziz = zmax
          ithetaz = 0.
          iphiz = 0.
        else if (vz(i) > 0.) then
          dtz = (zmin - z(i))/vz(i)
          xiz = vx(i)*dtz + x(i)
          yiz = vy(i)*dtz + y(i)
          ziz = zmin
          ithetaz = pi
          iphiz = 0.
        else
          dtz = -LARGEPOS
          xiz = LARGEPOS
          yiz = LARGEPOS
          ziz = zcent
          ithetaz = 0.
          iphiz = 0.
        endif

        if (dtx == -LARGEPOS .or. dty == -LARGEPOS) then
          xi(i) = xiz
          yi(i) = yiz
          zi(i) = ziz
          itheta(i) = ithetaz
          iphi(i) = iphiz
          cycle
        endif

        linside = .true.
        do while (linside)

          if (-dtz < min(-dtx,-dty)) then
            xi(i) = xiz
            yi(i) = yiz
            zi(i) = ziz
            itheta(i) = ithetaz
            iphi(i) = iphiz
            exit
          endif

          if (-dtx < -dty) then
            xi(i) = vx(i)*dtx + x(i)
            yi(i) = vy(i)*dtx + y(i)
            zi(i) = vz(i)*dtx + z(i)
            itheta(i) = ithetax
            iphi(i) = iphix
          else
            xi(i) = vx(i)*dty + x(i)
            yi(i) = vy(i)*dty + y(i)
            zi(i) = vz(i)*dty + z(i)
            itheta(i) = ithetay
            iphi(i) = iphiy
          endif

          call ZGridConductorD(xcellsize,ycellsize,length,thickness,
     &                         xcent,ycent,zcent,
     &                         1,xi(i),yi(i),zi(i),dd)

c         linside = (dd < 0.)
          linside = (dd < -thickness/1000.)

          if (linside) then
            if (-dtx < -dty) then
              dtx = dtx - xcellsize/abs(vx(i))
            else
              dty = dty - ycellsize/abs(vy(i))
            endif
          endif

        enddo

      enddo

      return
      end

c=======================================================================
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
          xr1 = sqrt(max(0.,xr1sq))
          yr1 = sqrt(max(0.,yr1sq))
          xr2 = sqrt(max(0.,xr2sq))
          yr2 = sqrt(max(0.,yr2sq))
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
            aa = (2.*z0*za - z0**2 + x1**2)/(za+sqrt(max(0.,(za-z0)**2-x1**2)))
            z1 = (2.*aa*zb - aa**2 + y1**2)/(zb+sqrt(max(0.,(zb-aa)**2-y1**2)))
            bb = (2.*z0*zb - z0**2 + y1**2)/(zb+sqrt(max(0.,(zb-z0)**2-y1**2)))
            z2 = (2.*bb*za - bb**2 + x1**2)/(za+sqrt(max(0.,(za-bb)**2-x1**2)))
            zave = 0.5*(z1 + z2)
            x1 = xx*(za - zave)/(za-zz)
            y1 = yy*(zb - zave)/(zb-zz)
          enddo
        else
          zave = z0
        endif

        dd = sqrt((xx - x1)**2 + (yy - y1)**2 + (zz - zave)**2)
        distance(i) = dd - thickness*0.5

      enddo

      return
      end
c=======================================================================
      subroutine BeamletplateIntercept(za,zb,z0,thickness,xcent,ycent,zcent,
     &                                 n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      real(kind=8):: za,zb,z0,thickness,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i

      print*,"The intercept calculation is not yet implemented for Beamletplate"
      do i=1,n
        xi(i) = LARGEPOS
        yi(i) = LARGEPOS
        zi(i) = LARGEPOS
        itheta(i) = 0.
        iphi(i) = 0.
      enddo

      return
      end
c=======================================================================
c=======================================================================
c=======================================================================
      subroutine ZSrfrvConductorFnew(nn,rsrf,zsrf,rad,rc,zc,xcent,ycent,zcent,
     &                               intercepts,fuzz)
      use Constant
      use ConductorInterceptTypeModule
      integer(ISZ):: nn
      real(kind=8):: rsrf(nn),zsrf(nn),rad(nn-1),rc(nn-1),zc(nn-1)
      real(kind=8):: xcent,ycent,zcent
      type(ConductorInterceptType):: intercepts
      real(kind=8):: fuzz

c Finds the intercepts of the object with grid lines along the three axis.
c A positive zsign means that the conductor is the region above zcent.

      real(kind=8):: xmmin,ymmin,zmmin
      real(kind=8):: xmmax,ymmax,zmmax
      real(kind=8):: dx,dy,dz
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nxicpt,nyicpt,nzicpt

      real(kind=8),pointer:: xintercepts(:,:,:)
      real(kind=8),pointer:: yintercepts(:,:,:)
      real(kind=8),pointer:: zintercepts(:,:,:)

      real(kind=8),pointer:: rmin(:),rmax(:)
      real(kind=8),pointer:: zmin(:),zmax(:)
      real(kind=8),pointer:: thmin(:),thmax(:)
      real(kind=8),pointer:: th1(:),th2(:)
      integer(ISZ):: ix,iy,iz,ii,iim1,iip1
      real(kind=8):: xx,yy,zz
      real(kind=8):: rr,ww,rsq
      real(kind=8):: x1,x2,y1,y2,t1,t2,z1,z2
      logical(ISZ):: l1,l2,ll
      integer(ISZ):: nt,it
      real(kind=8):: rt(4),tt

      integer(ISZ):: nirflat,nizflat,iiflat,ic
      integer(ISZ),pointer:: irflat(:),izflat(:)

      xmmin = intercepts%xmmin
      ymmin = intercepts%ymmin
      zmmin = intercepts%zmmin
      dx = intercepts%dx
      dy = intercepts%dy
      dz = intercepts%dz
      nx = intercepts%nx
      ny = intercepts%ny
      nz = intercepts%nz

      intercepts%nxicpt = 2*(nn-1)
      intercepts%nyicpt = 2*(nn-1)
      intercepts%nzicpt = 2*(nn-1)

      allocate(rmin(nn-1),rmax(nn-1))
      allocate(zmin(nn-1),zmax(nn-1))
      allocate(thmin(nn-1),thmax(nn-1))
      allocate(th1(nn-1),th2(nn-1))

      do ii=1,nn-1

        rmin(ii) = min(rsrf(ii),rsrf(ii+1))
        rmax(ii) = max(rsrf(ii),rsrf(ii+1))
        zmin(ii) = min(zsrf(ii),zsrf(ii+1))
        zmax(ii) = max(zsrf(ii),zsrf(ii+1))

        if (rad(ii) .ne. LARGEPOS) then

c         --- The angle in the r-z plane of the start and end of the segment
          th1(ii) = atan2(rsrf(ii) - rc(ii),zsrf(ii) - zc(ii))
          th2(ii) = atan2(rsrf(ii+1) - rc(ii),zsrf(ii+1) - zc(ii))

c         --- If rad < 0, then the angle from start to end should be
c         --- increasing. Otherwise decreasing.
c         --- Note that both th1 and th2 will also be within the
c         --- range of -pi to 3*pi.
          if (th1(ii) > th2(ii) .and. rad(ii) < 0. .and. zsrf(ii) < zsrf(ii+1)) then
            th2(ii) = th2(ii) + 2*pi
          endif
          if (th1(ii) < th2(ii) .and. rad(ii) > 0. .and. zsrf(ii) < zsrf(ii+1)) then
            th1(ii) = th1(ii) + 2*pi
          endif

c         --- If the range of angles includes the extremum of the circle,
c         --- then the mins or maxs may need to be recalculated.
          thmin(ii) = min(th1(ii),th2(ii))
          thmax(ii) = max(th1(ii),th2(ii))
          if ((thmin(ii) < -pi*0.5 .and. thmax(ii) > -pi*0.5) .or.
     &        (thmin(ii) < pi*1.5 .and. thmax(ii) > pi*1.5)) then
            rmin(ii) = rc(ii) - abs(rad(ii))
            intercepts%nxicpt = intercepts%nxicpt + 2
            intercepts%nyicpt = intercepts%nyicpt + 2
          endif
          if ((thmin(ii) < pi*0.5 .and. thmax(ii) > pi*0.5) .or.
     &        (thmin(ii) < pi*2.5 .and. thmax(ii) > pi*2.5)) then
            rmax(ii) = rc(ii) + abs(rad(ii))
            intercepts%nxicpt = intercepts%nxicpt + 2
            intercepts%nyicpt = intercepts%nyicpt + 2
          endif
          if (thmin(ii) < pi .and. thmax(ii) > pi) then
            zmin(ii) = zc(ii) - abs(rad(ii))
            intercepts%nzicpt = intercepts%nzicpt + 2
          endif
          if ((thmin(ii) < 0. .and. thmax(ii) > 0.) .or.
     &        (thmin(ii) < 2*pi .and. thmax(ii) > 2*pi)) then
            zmax(ii) = zc(ii) + abs(rad(ii))
            intercepts%nzicpt = intercepts%nzicpt + 2
          endif

        endif

      enddo

      call ConductorInterceptTypeallot(intercepts)

      xintercepts => intercepts%xintercepts
      yintercepts => intercepts%yintercepts
      zintercepts => intercepts%zintercepts

c     --- irflat and izflat tag line segments that have the same r or z at
c     --- both ends and that align with an r or z grid cell. Some special
c     --- handling is needed to gaurantee that the line segment is considered
c     --- inside of the object. Note that the same array can be used for
c     --- both x and y.
      allocate(izflat(intercepts%nxicpt))
      allocate(irflat(intercepts%nzicpt))

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do iy=0,ny
          yy = ymmin + iy*dy - ycent

          nxicpt = 0
          nizflat = 0

c         --- Loop over r-z line segments
          do ii=1,nn-1

c           --- Check if the z of the grid is within the z range of the segment
c           --- ll will be true if the point (zz) is within the extent of
c           --- the line segment, including the starting point and
c           --- excluding the end point.
            ll = (zmin(ii) < zz .and. zz < zmax(ii)) .or. zz == zsrf(ii)
            if (zz == zsrf(ii+1)) then
c             --- In this end case, only include the point if the lines
c             --- before and after it are both either above or below it.
c             --- This will be a peak - two intercepts are gathered.
              iip1 = ii + 1
              if (ii == nn-1) iip1 = 1
              if ((zmin(ii) < zz .and. zmin(iip1) < zz) .or.
     &            (zmax(ii) > zz .and. zmax(iip1) > zz)) then
                ll = .true.
              endif
            endif
            if (ll) then
              if (rad(ii) == LARGEPOS) then
                if (zsrf(ii+1) - zsrf(ii) == 0.) then
c                 --- Based on the logic for ll above, this will only ever
c                 --- happen if it is also true that zz = zsrf(ii).
c                 --- This is tricky (and still not quite right).
c                 --- Only include the intercept if the lines before and
c                 --- after it are both either above or below it. In that
c                 --- case, this is a local min or max and is treated the
c                 --- same as a peak - the intercepts at the two ends of
c                 --- segment will be included. This does
c                 --- get the number of intercepts correct, but doesn't
c                 --- consistently make the line inside the object.
c                 --- Further processing is needed after scanning all of
c                 --- the line segments.
c                 --- Flag this line segment.
                  nizflat = nizflat + 1
                  izflat(nizflat) = ii
c                 --- Get the line segments before and after.
                  iim1 = ii - 1
                  if (iim1 == 0) iim1 = nn-1
                  iip1 = ii + 1
                  if (iip1 == nn) iip1 = 1
c                 --- Search for the next line that has a change in radius.
c                 --- This is done to skip sequential line segments that
c                 --- are all at the same radius.
                  do while (zmin(iip1) == zmax(iip1) .and. iip1 .ne. ii)
                    iip1 = iip1 + 1
                    if (iip1 == nn) iip1 = 1
                  enddo
c                 --- Note that since a search is not done for iim1, only
c                 --- the first a of series of line segments would ever
c                 --- be included.
                  if ((zmin(iim1) < zz .and. zmin(iip1) < zz) .or.
     &                (zmax(iim1) > zz .and. zmax(iip1) > zz)) then
                    rr = rsrf(ii)
                    if (yy**2 <= rr**2) then
                      x1 = -sqrt(max(0.,rr**2 - yy**2))
                      x2 = +sqrt(max(0.,rr**2 - yy**2))
                      call intserticpt(x1 + xcent,nxicpt,iy,iz,
     &                                 intercepts%nxicpt,ny,nz,xintercepts)
                      call intserticpt(x2 + xcent,nxicpt,iy,iz,
     &                                 intercepts%nxicpt,ny,nz,xintercepts)
c                     nxicpt = nxicpt + 1
c                     xintercepts(nxicpt,iy,iz) = x1 + xcent
c                     nxicpt = nxicpt + 1
c                     xintercepts(nxicpt,iy,iz) = x2 + xcent
                    endif
                  endif
                else
c                 --- A straight line segment is straightforward
                  ww = (zz - zsrf(ii))/(zsrf(ii+1) - zsrf(ii))
                  rr = rsrf(ii)*(1.-ww) + rsrf(ii+1)*ww
                  if (yy**2 <= rr**2) then
                    x1 = -sqrt(max(0.,rr**2 - yy**2))
                    x2 = +sqrt(max(0.,rr**2 - yy**2))
                    call intserticpt(x1 + xcent,nxicpt,iy,iz,
     &                               intercepts%nxicpt,ny,nz,xintercepts)
                    call intserticpt(x2 + xcent,nxicpt,iy,iz,
     &                               intercepts%nxicpt,ny,nz,xintercepts)
c                   nxicpt = nxicpt + 1
c                   xintercepts(nxicpt,iy,iz) = x1 + xcent
c                   nxicpt = nxicpt + 1
c                   xintercepts(nxicpt,iy,iz) = x2 + xcent
                  endif
                endif
              else
c               --- A circular arc is complicated. This is a torus.
c               --- There can be up to four intercept points.
c               --- Check each one, starting with the one that has the
c               --- smallest x intercept. Only the first two have to be
c               --- carefully checked. The second two would be the same,
c               --- but with the opposite sign of x. This checks if the
c               --- angle of the intercept point is within the angular
c               --- range of the arc. Note that the starting point, th1,
c               --- is included, whereas the ending point, th2, is not.
                rsq = rad(ii)**2 - (zz - zc(ii))**2
                nt = 0

                rr = sqrt(max(0.,rsq))
                tt = atan2(+rr,zz - zc(ii))
                if (tt < thmin(ii)) tt = tt + 2*pi
                if (thmin(ii) <= tt .and. tt <= thmax(ii)) then
                  nt = nt + 1
                  rt(nt) = -rc(ii) - rr
                endif
                tt = atan2(-rr,zz - zc(ii))
                if (tt < thmin(ii)) tt = tt + 2*pi
                if (thmin(ii) <= tt .and. tt <= thmax(ii)) then
                  nt = nt + 1
                  rt(nt) = -rc(ii) + rr
                endif

c               --- Get the possible upper points.
                if (nt == 1) then
                  nt = nt + 1
                  rt(2) = -rt(1)
                else if (nt == 2) then
                  nt = nt + 2
                  rt(3) = -rt(2)
                  rt(4) = -rt(1)
                endif

c               --- Check if any of those four points intercept the grid line
                do it=1,nt
                  if (yy**2 <= rt(it)**2) then
                    x1 = sqrt(max(0.,rt(it)**2 - yy**2))
                    if (rt(it) < 0.) x1 = -x1
                    call intserticpt(x1 + xcent,nxicpt,iy,iz,
     &                               intercepts%nxicpt,ny,nz,xintercepts)
c                   nxicpt = nxicpt + 1
c                   xintercepts(nxicpt,iy,iz) = x1 + xcent
                  endif
                enddo

              endif ! straight or curved segment
            endif ! within z range
          enddo ! segments

c         --- Do any clean up that is needed if there are line segments
c         --- that are flat in z. This needs to be done after the loop over
c         --- line segments since all of the intercepts are needed.
c         --- This is rather ugly code, but that is typical for code handling
c         --- the edge cases.
c         --- Loop over the flat line segments
          do iiflat=1,nizflat
            ii = izflat(iiflat)

c           --- If yy is outside the larger radius, then there are no
c           --- intercepts and nothing needs to be done.
            if (yy > max(rsrf(ii),rsrf(ii+1))) cycle

c           --- Note that x1 and x2 have different meanings that above.
c           --- Here they are the lower and upper ends of the line segment.
            x1 = sqrt(max(0.,min(rsrf(ii),rsrf(ii+1))**2 - yy**2))
            x2 = sqrt(max(0.,max(rsrf(ii),rsrf(ii+1))**2 - yy**2))

c           --- Use a while loop since nxicpt may change
            ic = 1
            do while (ic < nxicpt)

              if (xintercepts(ic,iy,iz) < -x2 + xcent .and.
     &            xintercepts(ic+1,iy,iz) == -x2 + xcent) then
c               --- Replace the lower point with the upper point, putting
c               --- edge inside of the object.
                xintercepts(ic+1,iy,iz) = -x1 + xcent

              else if (xintercepts(ic,iy,iz) == -x1 + xcent .and.
     &                 xintercepts(ic+1,iy,iz) > -x1 + xcent) then
c               --- Replace the upper point with the lower point, putting
c               --- edge inside of the object.
                xintercepts(ic,iy,iz) = -x2 + xcent
c               --- It can happen that, with the above change,
c               --- xintercepts(ic) will be less than xintercept(ic-1). This
c               --- can only happen if (ic-1) is the lower edge of the
c               --- line segment. The two points can be removed since they
c               --- are both on the segment and therefore inside.
                if (ic > 1) then
                  if (xintercepts(ic-1,iy,iz) == -x1 + xcent) then
                    xintercepts(ic-1:nxicpt-3,iy,iz) = xintercepts(ic+1:nxicpt-1,iy,iz)
                    nxicpt = nxicpt - 2
                    ic = ic - 2
                  endif
                endif

              else if (xintercepts(ic,iy,iz) > -x2 + xcent .and.
     &                 xintercepts(ic+1,iy,iz) < -x1 + xcent) then
c               --- Remove the two points so that the edge between
c               --- them is inside the object.
                xintercepts(ic:nxicpt-2,iy,iz) = xintercepts(ic+2:nxicpt,iy,iz)
                nxicpt = nxicpt - 2
                ic = ic - 2

              else if (xintercepts(ic,iy,iz) < x1 + xcent .and.
     &                 xintercepts(ic+1,iy,iz) == x1 + xcent) then
c               --- Replace the lower point with the upper point, putting
c               --- edge inside of the object.
                xintercepts(ic+1,iy,iz) = x2 + xcent

              else if (xintercepts(ic,iy,iz) == x2 + xcent .and.
     &                 xintercepts(ic+1,iy,iz) > x2 + xcent) then
c               --- Replace the upper point with the lower point, putting
c               --- edge inside of the object.
                xintercepts(ic,iy,iz) = x1 + xcent
c               --- It can happen that, with the above change,
c               --- xintercepts(ic) will be less than xintercept(ic-1). This
c               --- can only happen if (ic-1) is the lower edge of the
c               --- line segment. The two points can be removed since they
c               --- are both on the segment and therefore inside.
                if (ic > 1) then
                  if (xintercepts(ic-1,iy,iz) == x2 + xcent) then
                    xintercepts(ic-1:nxicpt-3,iy,iz) = xintercepts(ic+1:nxicpt-1,iy,iz)
                    nxicpt = nxicpt - 2
                    ic = ic - 2
                  endif
                endif

              else if (xintercepts(ic,iy,iz) > x1 + xcent .and.
     &                 xintercepts(ic+1,iy,iz) < x2 + xcent) then
c               --- Remove the two points so that the edge between
c               --- them is inside the object.
                xintercepts(ic:nxicpt-2,iy,iz) = xintercepts(ic+2:nxicpt,iy,iz)
                nxicpt = nxicpt - 2
                ic = ic - 2

              endif
              ic = ic + 2
            enddo
          enddo

          xintercepts(nxicpt+1:,iy,iz) = LARGEPOS

        enddo ! iy
      enddo ! iz

      do iz=0,nz
        zz = zmmin + iz*dz - zcent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent

          nyicpt = 0
          nizflat = 0

c         --- Loop over r-z line segments
          do ii=1,nn-1

c           --- Check if the z of the grid is within the z range of the segment
c           --- ll will be true if the point (zz) is within the extent of
c           --- the line segment, including the starting point and
c           --- excluding the end point.
            ll = (zmin(ii) < zz .and. zz < zmax(ii)) .or. zz == zsrf(ii)
            if (zz == zsrf(ii+1)) then
c             --- In this end case, only include the point if the lines
c             --- before and after it are both either above or below it.
c             --- This will be a peak - two intercepts are gathered.
              iip1 = ii + 1
              if (ii == nn-1) iip1 = 1
              if ((zmin(ii) < zz .and. zmin(iip1) < zz) .or.
     &            (zmax(ii) > zz .and. zmax(iip1) > zz)) then
                ll = .true.
              endif
            endif
            if (ll) then

              if (rad(ii) == LARGEPOS) then
c               --- A straight line segment is straightforward
                if (zsrf(ii+1) - zsrf(ii) == 0.) then
c                 --- Based on the logic for ll above, this will only ever
c                 --- happen if it is also true that zz = zsrf(ii).
c                 --- This is tricky (and still not quite right).
c                 --- Only include the intercept if the lines before and
c                 --- after it are both either above or below it. In that
c                 --- case, this is a local min or max and is treated the
c                 --- same as a peak - the intercepts at the two ends of
c                 --- segment will be included. This does
c                 --- get the number of intercepts correct, but doesn't
c                 --- consistently make the line inside the object.
c                 --- Flag this line segment.
                  nizflat = nizflat + 1
                  izflat(nizflat) = ii
c                 --- Get the line segments before and after.
                  iim1 = ii - 1
                  if (iim1 == 0) iim1 = nn-1
                  iip1 = ii + 1
                  if (iip1 == nn) iip1 = 1
c                 --- Search for the next line that has a change in radius.
c                 --- This is done to skip sequential line segments that
c                 --- are all at the same radius.
                  do while (zmin(iip1) == zmax(iip1) .and. iip1 .ne. ii)
                    iip1 = iip1 + 1
                    if (iip1 == nn) iip1 = 1
                  enddo
c                 --- Note that since a search is not done for iim1, only
c                 --- the first a of series of line segments would ever
c                 --- be included.
                  if ((zmin(iim1) < zz .and. zmin(iip1) < zz) .or.
     &                (zmax(iim1) > zz .and. zmax(iip1) > zz)) then
                    rr = rsrf(ii)
                    if (xx**2 <= rr**2) then
                      y1 = -sqrt(max(0.,rr**2 - xx**2))
                      y2 = +sqrt(max(0.,rr**2 - xx**2))
                      call intserticpt(y1 + ycent,nyicpt,ix,iz,
     &                                 intercepts%nyicpt,nx,nz,yintercepts)
                      call intserticpt(y2 + ycent,nyicpt,ix,iz,
     &                                 intercepts%nyicpt,nx,nz,yintercepts)
c                     nyicpt = nyicpt + 1
c                     yintercepts(nyicpt,ix,iz) = y1 + ycent
c                     nyicpt = nxicpt + 1
c                     yintercepts(nyicpt,ix,iz) = y2 + ycent
                    endif
                  endif
                else
                  ww = (zz - zsrf(ii))/(zsrf(ii+1) - zsrf(ii))
                  rr = rsrf(ii)*(1.-ww) + rsrf(ii+1)*ww
                  if (xx**2 <= rr**2) then
                    y1 = -sqrt(max(0.,rr**2 - xx**2))
                    y2 = +sqrt(max(0.,rr**2 - xx**2))
                    call intserticpt(y1 + ycent,nyicpt,ix,iz,
     &                               intercepts%nyicpt,nx,nz,yintercepts)
                    call intserticpt(y2 + ycent,nyicpt,ix,iz,
     &                               intercepts%nyicpt,nx,nz,yintercepts)
c                   nyicpt = nyicpt + 1
c                   yintercepts(nyicpt,ix,iz) = y1 + ycent
c                   nyicpt = nyicpt + 1
c                   yintercepts(nyicpt,ix,iz) = y2 + ycent
                  endif
                endif
              else
c               --- A circular arc is complicated. This is a torus.
c               --- There can be up to four intercept points.
c               --- Check each one, starting with the one that has the
c               --- smallest y intercept. Only the first two have to be
c               --- carefully checked. The second two would be the same,
c               --- but with the opposite sign of y. This checks if the
c               --- angle of the intercept point is within the angular
c               --- range of the arc. Note that the starting point, th1,
c               --- is included, whereas the ending point, th2, is not.
                rsq = rad(ii)**2 - (zz - zc(ii))**2
                nt = 0

                rr = sqrt(max(0.,rsq))
                tt = atan2(+rr,zz - zc(ii))
                if (tt < thmin(ii)) tt = tt + 2*pi
                if (thmin(ii) <= tt .and. tt <= thmax(ii)) then
                  nt = nt + 1
                  rt(nt) = -rc(ii) - rr
                endif
                tt = atan2(-rr,zz - zc(ii))
                if (tt < thmin(ii)) tt = tt + 2*pi
                if (thmin(ii) <= tt .and. tt <= thmax(ii)) then
                  nt = nt + 1
                  rt(nt) = -rc(ii) + rr
                endif

c               --- Get the possible upper points.
                if (nt == 1) then
                  nt = nt + 1
                  rt(2) = -rt(1)
                else if (nt == 2) then
                  nt = nt + 2
                  rt(3) = -rt(2)
                  rt(4) = -rt(1)
                endif

c               --- Check if any of those four points intercept the grid line
                do it=1,nt
                  if (xx**2 <= rt(it)**2) then
                    y1 = sqrt(max(0.,rt(it)**2 - xx**2))
                    if (rt(it) < 0.) y1 = -y1
                    call intserticpt(y1 + ycent,nyicpt,ix,iz,
     &                               intercepts%nyicpt,nx,nz,yintercepts)
c                   nyicpt = nyicpt + 1
c                   yintercepts(nyicpt,ix,iz) = y1 + ycent
                  endif
                enddo

              endif ! straight or curved segment
            endif ! within z range
          enddo ! segments

c         --- Do any clean up that is needed if there are line segments
c         --- that are flat in z. This needs to be done after the loop over
c         --- line segments since all of the intercepts are needed.
c         --- This is rather ugly code, but that is typical for code handling
c         --- the edge cases.
c         --- Loop over the flat line segments
          do iiflat=1,nizflat
            ii = izflat(iiflat)

c           --- If xx is outside the larger radius, then there are no
c           --- intercepts and nothing needs to be done.
            if (xx > max(rsrf(ii),rsrf(ii+1))) cycle

c           --- Note that y1 and y2 have different meanings that above.
c           --- Here they are the lower and upper ends of the line segment.
            y1 = sqrt(max(0.,min(rsrf(ii),rsrf(ii+1))**2 - xx**2))
            y2 = sqrt(max(0.,max(rsrf(ii),rsrf(ii+1))**2 - xx**2))

c           --- Use a while loop since nyicpt may change
            ic = 1
            do while (ic < nyicpt)

              if (yintercepts(ic,ix,iz) < -y2 + ycent .and.
     &            yintercepts(ic+1,ix,iz) == -y2 + ycent) then
c               --- Replace the lower point with the upper point, putting
c               --- edge inside of the object.
                yintercepts(ic+1,ix,iz) = -y1 + ycent

              else if (yintercepts(ic,ix,iz) == -y1 + ycent .and.
     &                 yintercepts(ic+1,ix,iz) > -y1 + ycent) then
c               --- Replace the upper point with the lower point, putting
c               --- edge inside of the object.
                yintercepts(ic,ix,iz) = -y2 + ycent
c               --- It can happen that, with the above change,
c               --- yintercepts(ic) will be less than yintercept(ic-1). This
c               --- can only happen if (ic-1) is the lower edge of the
c               --- line segment. The two points can be removed since they
c               --- are both on the segment and therefore inside.
                if (ic > 1) then
                  if (yintercepts(ic-1,ix,iz) == -y1 + ycent) then
                    yintercepts(ic-1:nyicpt-3,ix,iz) = yintercepts(ic+1:nyicpt-1,ix,iz)
                    nyicpt = nyicpt - 2
                    ic = ic - 2
                  endif
                endif

              else if (yintercepts(ic,ix,iz) > -y2 + ycent .and.
     &                 yintercepts(ic+1,ix,iz) < -y1 + ycent) then
c               --- Remove the two points so that the edge between
c               --- them is inside the object.
                yintercepts(ic:nyicpt-2,ix,iz) = yintercepts(ic+2:nyicpt,ix,iz)
                nyicpt = nyicpt - 2
                ic = ic - 2

              else if (yintercepts(ic,ix,iz) < y1 + ycent .and.
     &                 yintercepts(ic+1,ix,iz) == y1 + ycent) then
c               --- Replace the lower point with the upper point, putting
c               --- edge inside of the object.
                yintercepts(ic+1,ix,iz) = y2 + ycent

              else if (yintercepts(ic,ix,iz) == y2 + ycent .and.
     &                 yintercepts(ic+1,ix,iz) > y2 + ycent) then
c               --- Replace the upper point with the lower point, putting
c               --- edge inside of the object.
                yintercepts(ic,ix,iz) = y1 + ycent
c               --- It can happen that, with the above change,
c               --- yintercepts(ic) will be less than yintercept(ic-1). This
c               --- can only happen if (ic-1) is the lower edge of the
c               --- line segment. The two points can be removed since they
c               --- are both on the segment and therefore inside.
                if (ic > 1) then
                  if (yintercepts(ic-1,ix,iz) == y2 + ycent) then
                    yintercepts(ic-1:nyicpt-3,ix,iz) = yintercepts(ic+1:nyicpt-1,ix,iz)
                    nyicpt = nyicpt - 2
                    ic = ic - 2
                  endif
                endif

              else if (yintercepts(ic,ix,iz) > y1 + ycent .and.
     &                 yintercepts(ic+1,ix,iz) < y2 + ycent) then
c               --- Remove the two points so that the edge between
c               --- them is inside the object.
                yintercepts(ic:nyicpt-2,ix,iz) = yintercepts(ic+2:nyicpt,ix,iz)
                nyicpt = nyicpt - 2
                ic = ic - 2

              endif
              ic = ic + 2
            enddo
          enddo

          yintercepts(nyicpt+1:,ix,iz) = LARGEPOS

        enddo ! ix
      enddo ! iz


      do iy=0,ny
        yy = ymmin + iy*dy - ycent
        do ix=0,nx
          xx = xmmin + ix*dx - xcent
          rr = sqrt(xx**2 + yy**2)

          nzicpt = 0
          nirflat = 0

c         --- Loop over r-z line segments
          do ii=1,nn-1

c           --- ll will be true if the point (rr) is within the radius of
c           --- the line segment, including the starting point and
c           --- excluding the end point.
            ll = (rmin(ii) < rr .and. rr < rmax(ii)) .or. rr == rsrf(ii)
            if (rr == rsrf(ii+1)) then
c             --- In this end case, only include the point if the lines
c             --- before and after it are both either above or below it.
c             --- This will be a peak - two intercepts are gathered.
              iip1 = ii + 1
              if (ii == nn-1) iip1 = 1
              if ((rmin(ii) < rr .and. rmin(iip1) < rr) .or.
     &            (rmax(ii) > rr .and. rmax(iip1) > rr)) then
                ll = .true.
              endif
            endif
            if (ll) then
              if (rad(ii) == LARGEPOS) then
                if (rsrf(ii+1) - rsrf(ii) == 0.) then
c                 --- Based on the logic for ll above, this will only ever
c                 --- happen if it is also true that rr = rsrf(ii).
c                 --- This is tricky (and still not quite right).
c                 --- Only include the intercept if the lines before and
c                 --- after it are both either above or below it. In that
c                 --- case, this is a local min or max and is treated the
c                 --- same as a peak - the intercepts at the two ends of
c                 --- segment will be included. This does
c                 --- get the number of intercepts correct, but doesn't
c                 --- consistently make the line inside the object.
c                 --- Flag this line segment.
                  nirflat = nirflat + 1
                  irflat(nirflat) = ii
c                 --- Get the line segments before and after.
                  iim1 = ii - 1
                  if (iim1 == 0) iim1 = nn-1
                  iip1 = ii + 1
                  if (iip1 == nn) iip1 = 1
c                 --- Search for the next line that has a change in radius.
c                 --- This is done to skip sequential line segments that
c                 --- are all at the same radius.
                  do while (rmin(iip1) == rmax(iip1) .and. iip1 .ne. ii)
                    iip1 = iip1 + 1
                    if (iip1 == nn) iip1 = 1
                  enddo
c                 --- Note that since a search is not done for iim1, only
c                 --- the first a of series of line segments would ever
c                 --- be included.
                  if ((rmin(iim1) < rr .and. rmin(iip1) < rr) .or.
     &                (rmax(iim1) > rr .and. rmax(iip1) > rr)) then
                    zz = zsrf(ii)
                    call intserticpt(zz + zcent,nzicpt,ix,iy,
     &                               intercepts%nzicpt,nx,ny,zintercepts)
c                   nzicpt = nzicpt + 1
c                   zintercepts(nzicpt,ix,iy) = zz + zcent
                  endif
                else
                  ww = (rr - rsrf(ii))/(rsrf(ii+1) - rsrf(ii))
                  zz = zsrf(ii)*(1.-ww) + zsrf(ii+1)*ww
                  call intserticpt(zz + zcent,nzicpt,ix,iy,
     &                             intercepts%nzicpt,nx,ny,zintercepts)
c                 nzicpt = nzicpt + 1
c                 zintercepts(nzicpt,ix,iy) = zz + zcent
                endif
              else

c               --- A circular arc is complicated. This is a torus.
c               --- There can be up to two intercept points.
c               --- Check each one, starting with the one that has the
c               --- smallest z intercept. This checks if the
c               --- angle of the intercept point is within the angular
c               --- range of the arc. Note that the starting point, th1,
c               --- is included, whereas the ending point, th2, is not.
                z1 = sqrt(max(0.,rad(ii)**2 - (rr - rc(ii))**2))

                tt = atan2(rr - rc(ii),-z1)
                if (tt < thmin(ii)) tt = tt + 2*pi
                if (thmin(ii) < tt .and. tt < thmax(ii) .or. tt == th1(ii)) then
                  call intserticpt(zc(ii) - z1 + zcent,nzicpt,ix,iy,
     &                             intercepts%nzicpt,nx,ny,zintercepts)
c                 nzicpt = nzicpt + 1
c                 zintercepts(nzicpt,ix,iy) = zc(ii) - z1 + zcent
                endif
                
                tt = atan2(rr - rc(ii),+z1)
                if (tt < thmin(ii)) tt = tt + 2*pi
                if (thmin(ii) < tt .and. tt < thmax(ii) .or. tt == th1(ii)) then
                  call intserticpt(zc(ii) + z1 + zcent,nzicpt,ix,iy,
     &                             intercepts%nzicpt,nx,ny,zintercepts)
c                 nzicpt = nzicpt + 1
c                 zintercepts(nzicpt,ix,iy) = zc(ii) + z1 + zcent
                endif

              endif ! straight or curved segment
            endif ! within r range
          enddo ! segments

c         --- Do any clean up that is needed if there are line segments
c         --- that are flat in r. This needs to be done after the loop over
c         --- line segments since all of the intercepts are needed.
c         --- This is rather ugly code, but that is typical for code handling
c         --- the edge cases.
c         --- Loop over the flat line segments
          do iiflat=1,nirflat
            ii = irflat(iiflat)

c           --- Note that z1 and z2 have different meanings that above.
c           --- Here they are the lower and upper ends of the line segment.
            z1 = min(zsrf(ii),zsrf(ii+1))
            z2 = max(zsrf(ii),zsrf(ii+1))

c           --- Use a while loop since nzicpt may change
            ic = 1
            do while (ic < nzicpt)

              if (zintercepts(ic,ix,iy) < z1 + zcent .and.
     &            zintercepts(ic+1,ix,iy) == z1 + zcent) then
c               --- Replace the lower point with the upper point, putting
c               --- edge inside of the object.
                zintercepts(ic+1,ix,iy) = z2 + zcent

              else if (zintercepts(ic,ix,iy) == z2 + zcent .and.
     &                 zintercepts(ic+1,ix,iy) > z2 + zcent) then
c               --- Replace the upper point with the lower point, putting
c               --- edge inside of the object.
                zintercepts(ic,ix,iy) = z1 + zcent
c               --- It can happen that, with the above change,
c               --- zintercepts(ic) will be less than zintercept(ic-1). This
c               --- can only happen if (ic-1) is the lower edge of the
c               --- line segment. The two points can be removed since they
c               --- are both on the segment and therefore inside.
                if (ic > 1) then
                  if (zintercepts(ic-1,ix,iy) == z2 + zcent) then
                    zintercepts(ic-1:nzicpt-3,ix,iy) = zintercepts(ic+1:nzicpt-1,ix,iy)
                    nzicpt = nzicpt - 2
                    ic = ic - 2
                  endif
                endif

              else if (zintercepts(ic,ix,iy) > z1 + zcent .and.
     &                 zintercepts(ic+1,ix,iy) < z2 + zcent) then
c               --- Remove the two points so that the edge between
c               --- them is inside the object.
                zintercepts(ic:nzicpt-2,ix,iy) = zintercepts(ic+2:nzicpt,ix,iy)
                nzicpt = nzicpt - 2
                ic = ic - 2

              endif
              ic = ic + 2
            enddo
          enddo

          zintercepts(nzicpt+1:,ix,iy) = LARGEPOS

        enddo ! ix
      enddo ! iy

      deallocate(rmin,rmax)
      deallocate(zmin,zmax)
      deallocate(thmin,thmax)
      deallocate(th1,th2)
      deallocate(irflat,izflat)

      return
      CONTAINS

        subroutine intserticpt(icpt,nicpt,i1,i2,nicptmax,n1,n2,intercepts)
        real(kind=8):: icpt
        integer(ISZ):: nicpt,i1,i2
        integer(ISZ):: nicptmax,n1,n2
        real(kind=8):: intercepts(nicptmax,0:n1,0:n2)
        integer(ISZ):: ii

        nicpt = nicpt + 1

        if (nicpt == 1) then
          intercepts(nicpt,i1,i2) = icpt
          return
        endif

        do ii=nicpt-1,1,-1
          if (icpt >= intercepts(ii,i1,i2)) then
            intercepts(ii+1,i1,i2) = icpt
            exit
          else
            intercepts(ii+1,i1,i2) = intercepts(ii,i1,i2)
          endif
        enddo

        if (icpt < intercepts(1,i1,i2)) then
          intercepts(1,i1,i2) = icpt
        endif

        return
        end subroutine intserticpt

      end
c=======================================================================
      subroutine ZSrfrvConductorD(nn,rsrf,zsrf,rad,rc,zc,xcent,ycent,zcent,
     &                            n,x,y,z,distance)
      use Constant
      integer(ISZ):: nn
      real(kind=8):: rsrf(nn),zsrf(nn),rad(nn-1),rc(nn-1),zc(nn-1)
      real(kind=8):: xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: distance(n)

c Finds distances along the axis of the points from a surface of revolution.
c The notation is that a negative distance means that the point is inside.
c NOTICE: This only computes distances radially and can only reliably used
c to determine whether a point is inside or not. The true distance calculation
c would have to account for distances longitudinally.

      integer(ISZ):: nzicpt

      real(kind=8),pointer:: rmin(:),rmax(:)
      real(kind=8),pointer:: zmin(:),zmax(:)
      real(kind=8),pointer:: thmin(:),thmax(:)
      real(kind=8),pointer:: th1(:),th2(:)
      real(kind=8),pointer:: zintercepts(:)
      integer(ISZ):: ii,i,izicpt,iip1,iim1
      real(kind=8):: xx,yy,rr,zz,tt,rp,zp,tp,m
      real(kind=8):: ww,z1,dsq,dsqtemp
      logical(ISZ):: ll

      nzicpt = 2*(nn-1)

      allocate(rmin(nn-1),rmax(nn-1))
      allocate(zmin(nn-1),zmax(nn-1))
      allocate(thmin(nn-1),thmax(nn-1))
      allocate(th1(nn-1),th2(nn-1))

      do ii=1,nn-1

        rmin(ii) = min(rsrf(ii),rsrf(ii+1))
        rmax(ii) = max(rsrf(ii),rsrf(ii+1))
        zmin(ii) = min(zsrf(ii),zsrf(ii+1))
        zmax(ii) = max(zsrf(ii),zsrf(ii+1))

        if (rad(ii) .ne. LARGEPOS) then

c         --- The angle in the r-z plane of the start and end of the segment
          th1(ii) = atan2(rsrf(ii) - rc(ii),zsrf(ii) - zc(ii))
          th2(ii) = atan2(rsrf(ii+1) - rc(ii),zsrf(ii+1) - zc(ii))

c         --- If rad < 0, then the angle from start to end should be
c         --- increasing. Otherwise decreasing.
c         --- Note that both th1 and th2 will also be within the
c         --- range of -pi to 3*pi.
          if (th1(ii) > th2(ii) .and. rad(ii) < 0. .and. zsrf(ii) < zsrf(ii+1)) then
            th2(ii) = th2(ii) + 2*pi
          endif
          if (th1(ii) < th2(ii) .and. rad(ii) > 0. .and. zsrf(ii) < zsrf(ii+1)) then
            th1(ii) = th1(ii) + 2*pi
          endif

c         --- If the range of angles includes the extremum of the circle,
c         --- then the mins or maxs may need to be recalculated.
          thmin(ii) = min(th1(ii),th2(ii))
          thmax(ii) = max(th1(ii),th2(ii))
          if ((thmin(ii) < -pi*0.5 .and. thmax(ii) > -pi*0.5) .or.
     &        (thmin(ii) < pi*1.5 .and. thmax(ii) > pi*1.5)) then
            rmin(ii) = rc(ii) - abs(rad(ii))
          endif
          if ((thmin(ii) < pi*0.5 .and. thmax(ii) > pi*0.5) .or.
     &        (thmin(ii) < pi*2.5 .and. thmax(ii) > pi*2.5)) then
            rmax(ii) = rc(ii) + abs(rad(ii))
          endif
          if (thmin(ii) < pi .and. thmax(ii) > pi) then
            zmin(ii) = zc(ii) - abs(rad(ii))
            nzicpt = nzicpt + 2
          endif
          if ((thmin(ii) < 0. .and. thmax(ii) > 0.) .or.
     &        (thmin(ii) < 2*pi .and. thmax(ii) > 2*pi)) then
            zmax(ii) = zc(ii) + abs(rad(ii))
            nzicpt = nzicpt + 2
          endif

        endif

      enddo

      allocate(zintercepts(nzicpt))

      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        rr = sqrt(xx**2 + yy**2)
        zz = z(i) - zcent
        dsq = LARGEPOS
        izicpt = 0

        do ii=1,nn-1

c         --- Find the distance between the point and the line segment
          if (rad(ii) == LARGEPOS) then
c           --- The coordinates (rp,zp) will be the point on the line
c           --- closest to the point (rr,zz)
            if (zsrf(ii+1) == zsrf(ii)) then
              zp = zsrf(ii)
              rp = rr
            else
c             --- This expression was derived by finding the line
c             --- perpendicular to the line segment and that passes
c             --- through the point (rr,zz).
              m = (rsrf(ii+1) - rsrf(ii))/(zsrf(ii+1) - zsrf(ii))
              zp = (zz + m*rr + m**2*zsrf(ii) - m*rsrf(ii))/(1 + m**2)
              rp = m*(zp - zsrf(ii)) + rsrf(ii)
            endif
            if (rmin(ii) <= rp .and. rp <= rmax(ii) .and.
     &          zmin(ii) <= zp .and. zp <= zmax(ii)) then
c             --- The point (rp,zp) is within the line segment
              dsqtemp = (rr - rp)**2 + (zz - zp)**2
            else
c             --- The point (rp,zp) is not within the line segment, so
c             --- find the distance to the end of the line segment that
c             --- is closest to (rr,zz).
              dsqtemp = min((rsrf(ii  ) - rr)**2 + (zsrf(ii  ) - zz)**2,
     &                      (rsrf(ii+1) - rr)**2 + (zsrf(ii+1) - zz)**2)
            endif

          else

c           --- The angle of the line intercepting (rr,zz) and (rc,zc).
            tp = atan2(rr - rc(ii),zz - zc(ii))
            if (tp < thmin(ii)) tp = tp + 2*pi
            if (thmin(ii) <= tp .and. tp <= thmax(ii)) then
c             --- The angle is within the arc - get the point on the
c             --- circular arc at the angle.
              rp = abs(rad(ii))*sin(tp) + rc(ii)
              zp = abs(rad(ii))*cos(tp) + zc(ii)
              dsqtemp = (rr - rp)**2 + (zz - zp)**2
            else
c             --- The angle is not within the line segment, so
c             --- find the distance to the end of the line segment that
c             --- is closest to (rr,zz).
              dsqtemp = min((rsrf(ii  ) - rr)**2 + (zsrf(ii  ) - zz)**2,
     &                      (rsrf(ii+1) - rr)**2 + (zsrf(ii+1) - zz)**2)
            endif

          endif
          dsq = min(dsq,dsqtemp)

c         --- Find any intercepts along the z axis, so that it can be
c         --- determined whether (rr,zz) is inside or outside of the object.
c         --- ll will be true if the point (rr) is within the radius of
c         --- the line segment, including the starting point and
c         --- excluding the end point.
          ll = (rmin(ii) < rr .and. rr < rmax(ii)) .or. rr == rsrf(ii)
          if (rr == rsrf(ii+1)) then
c           --- In this end case, also include the point if the lines
c           --- before and after it are both either above or below it.
c           --- This will be a peak - two intercepts are gathered.
            iip1 = ii + 1
            if (ii == nn-1) iip1 = 1
            if ((rmin(ii) < rr .and. rmin(iip1) < rr) .or.
     &          (rmax(ii) > rr .and. rmax(iip1) > rr)) then
              ll = .true.
            endif
          endif
          if (ll) then
            if (rad(ii) == LARGEPOS) then
              if (rsrf(ii+1) == rsrf(ii)) then
c               --- Based on the logic for ll above, this will only ever
c               --- happen if it is also true that rr = rsrf(ii).
c               --- This is tricky (and still not quite right).
c               --- Only include the intercept if the lines before and
c               --- after it are both either above or below it. In that
c               --- case, this is a local min or max and is treated the
c               --- same as a peak - the intercepts at the two ends of
c               --- segment will be included. This does
c               --- get the number of intercepts correct, but doesn't
c               --- consistently make the line inside the object.
c               --- For the distance calculation though, this doesn't
c               --- matter since the distance would be zero if the
c               --- point is on the line segment. 
                iim1 = ii - 1
                if (iim1 == 0) iim1 = nn-1
                iip1 = ii + 1
                if (iip1 == nn) iip1 = 1
c               --- Search for the next line that has a change in radius.
c               --- This is done to skip sequential line segments that
c               --- are all at the same radius.
                do while (rmin(iip1) == rmax(iip1) .and. iip1 .ne. ii)
                  iip1 = iip1 + 1
                  if (iip1 == nn) iip1 = 1
                enddo
c               --- Note that since a search is not done for iim1, only
c               --- the first a of series of line segments would ever
c               --- be included.
                if ((rmin(iim1) < rr .and. rmin(iip1) < rr) .or.
     &              (rmax(iim1) > rr .and. rmax(iip1) > rr)) then
                  z1 = zsrf(ii)
                  call intsertzicpt(z1 + zcent)
                endif
              else
                ww = (rr - rsrf(ii))/(rsrf(ii+1) - rsrf(ii))
                z1 = zsrf(ii)*(1.-ww) + zsrf(ii+1)*ww
                call intsertzicpt(z1 + zcent)
              endif
            else

c             --- A circular arc is complicated. This is a torus.
c             --- There can be up to two intercept points.
c             --- Check each one, starting with the one that has the
c             --- smallest z intercept. This checks if the
c             --- angle of the intercept point is within the angular
c             --- range of the arc. Note that the starting point, th1,
c             --- is included, whereas the ending point, th2, is not.
              z1 = sqrt(max(0.,rad(ii)**2 - (rr - rc(ii))**2))
              if (rr == rsrf(ii)) z1 = abs(zsrf(ii) - zc(ii))
              if (rr == rsrf(ii+1)) z1 = abs(zsrf(ii+1) - zc(ii))

              tt = atan2(rr - rc(ii),-z1)
              if (tt < thmin(ii)) tt = tt + 2*pi
              if (thmin(ii) < tt .and. tt < thmax(ii) .or. tt == th1(ii)) then
                call intsertzicpt(zc(ii) - z1 + zcent)
              endif
                
              tt = atan2(rr - rc(ii),+z1)
              if (tt < thmin(ii)) tt = tt + 2*pi
              if (thmin(ii) < tt .and. tt < thmax(ii) .or. tt == th1(ii)) then
                call intsertzicpt(zc(ii) + z1 + zcent)
              endif

            endif ! straight or curved segment
          endif ! within r range
        enddo ! segments

        distance(i) = sqrt(max(0.,dsq))
        if (mod(izicpt,2) == 1) then
          print*,"ODD ",izicpt,zintercepts
        endif

        do ii=1,izicpt,2
          if (zintercepts(ii) <= zz .and. zz <= zintercepts(ii+1)) then
            distance(i) = -distance(i)
          endif
        enddo

      enddo

      deallocate(rmin,rmax)
      deallocate(zmin,zmax)
      deallocate(thmin,thmax)
      deallocate(th1,th2)
      deallocate(zintercepts)

      return
      CONTAINS

        subroutine intsertzicpt(zicpt)
        real(kind=8):: zicpt
        integer(ISZ):: ii

        izicpt = izicpt + 1

        if (izicpt == 1) then
          zintercepts(izicpt) = zicpt
          return
        endif

        do ii=izicpt-1,1,-1
          if (zicpt >= zintercepts(ii)) then
            zintercepts(ii+1) = zicpt
            exit
          else
            zintercepts(ii+1) = zintercepts(ii)
          endif
        enddo

        if (zicpt < zintercepts(1)) then
          zintercepts(1) = zicpt
        endif

        return
        end subroutine intsertzicpt

      end
c=======================================================================
      subroutine ZSrfrvIntercept(nn,rsrf,zsrf,rad,rc,zc,xcent,ycent,zcent,
     &                           n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      integer(ISZ):: nn
      real(kind=8):: zsrf(nn),rsrf(nn),rad(nn-1),zc(nn-1),rc(nn-1)
      real(kind=8):: xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      real(kind=8),pointer:: rmin(:),rmax(:)
      real(kind=8),pointer:: zmin(:),zmax(:)
      real(kind=8),pointer:: thmin(:),thmax(:)
      real(kind=8),pointer:: th1(:),th2(:)
      integer(ISZ):: i,ii
      real(kind=8):: vsq,xx,yy,zz,rr,r1,t1,t2,z1,z2,tt,ttemp,a1,a2,vr
      real(kind=8):: a,b,c,s,p
      real(kind=8):: vdotn,dd,th

      allocate(rmin(nn-1),rmax(nn-1))
      allocate(zmin(nn-1),zmax(nn-1))
      allocate(thmin(nn-1),thmax(nn-1))
      allocate(th1(nn-1),th2(nn-1))

      do ii=1,nn-1

        rmin(ii) = min(rsrf(ii),rsrf(ii+1))
        rmax(ii) = max(rsrf(ii),rsrf(ii+1))
        zmin(ii) = min(zsrf(ii),zsrf(ii+1))
        zmax(ii) = max(zsrf(ii),zsrf(ii+1))

        if (rad(ii) .ne. LARGEPOS) then

c         --- The angle in the r-z plane of the start and end of the segment
          th1(ii) = atan2(rsrf(ii) - rc(ii),zsrf(ii) - zc(ii))
          th2(ii) = atan2(rsrf(ii+1) - rc(ii),zsrf(ii+1) - zc(ii))

c         --- If rad < 0, then the angle from start to end should be
c         --- increasing. Otherwise decreasing.
c         --- Note that both th1 and th2 will also be within the
c         --- range of -pi to 3*pi.
          if (th1(ii) > th2(ii) .and. rad(ii) < 0. .and. zsrf(ii) < zsrf(ii+1)) then
            th2(ii) = th2(ii) + 2*pi
          endif
          if (th1(ii) < th2(ii) .and. rad(ii) > 0. .and. zsrf(ii) < zsrf(ii+1)) then
            th1(ii) = th1(ii) + 2*pi
          endif

c         --- If the range of angles includes the extremum of the circle,
c         --- then the mins or maxs may need to be recalculated.
          thmin(ii) = min(th1(ii),th2(ii))
          thmax(ii) = max(th1(ii),th2(ii))
          if ((thmin(ii) < -pi*0.5 .and. thmax(ii) > -pi*0.5) .or.
     &        (thmin(ii) < pi*1.5 .and. thmax(ii) > pi*1.5)) then
            rmin(ii) = rc(ii) - abs(rad(ii))
          endif
          if ((thmin(ii) < pi*0.5 .and. thmax(ii) > pi*0.5) .or.
     &        (thmin(ii) < pi*2.5 .and. thmax(ii) > pi*2.5)) then
            rmax(ii) = rc(ii) + abs(rad(ii))
          endif
          if (thmin(ii) < pi .and. thmax(ii) > pi) then
            zmin(ii) = zc(ii) - abs(rad(ii))
          endif
          if ((thmin(ii) < 0. .and. thmax(ii) > 0.) .or.
     &        (thmin(ii) < 2*pi .and. thmax(ii) > 2*pi)) then
            zmax(ii) = zc(ii) + abs(rad(ii))
          endif

        endif

      enddo

      do i=1,n

        tt = LARGEPOS
        xi(i) = LARGEPOS
        yi(i) = LARGEPOS
        zi(i) = LARGEPOS
        itheta(i) = 0.
        iphi(i) = 0.

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

c       --- Loop over the line segments and calculate the intersection of
c       --- each with the trajectory. Pick the one where the intersect
c       --- is within the segment and the intercept time is smallest.
        do ii=1,nn-1
          if (rad(ii) == LARGEPOS) then
            if (zsrf(ii+1) == zsrf(ii)) then
c             --- Special case of a z plane, where there is only one
c             --- possible intercept. If vz is zero, then there is no
c             --- intercept, so skip.
              if (vz(i) == 0.) cycle
              t1 = ((z(i)-zcent) - zsrf(ii))/vz(i)
              xx = -vx(i)*t1 + x(i) - xcent
              yy = -vy(i)*t1 + y(i) - ycent
              rr = sqrt(xx**2 + yy**2)
              if (rr >= min(rsrf(ii),rsrf(ii+1)) .and.
     &            rr <= max(rsrf(ii),rsrf(ii+1))) then
                if ((tt > t1 .and. t1 > 0. .and. tt > 0.) .or.
     &              (t1 > tt .and. tt < 0.)) then
                  tt = t1
                  itheta(i) = 0.
                  iphi(i) = 0.
                endif
              endif
            else
              s = (rsrf(ii+1) - rsrf(ii))/(zsrf(ii+1) - zsrf(ii))
              p = rsrf(ii) - s*zsrf(ii)
              a = (s*vz(i))**2 - vx(i)**2 - vy(i)**2
              b = 2*(-s*vz(i)*(s*(z(i)-zcent)+p) +
     &            vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
              c = (s*(z(i)-zcent)+p)**2 - (x(i)-xcent)**2 - (y(i)-ycent)**2
              if ((b**2 - 4*a*c) > 0.) then
                if (a .ne. 0.) then
                  t1 = (-b - sqrt(b**2 - 4*a*c))/(2*a)
                  t2 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
                else
                  t1 = -c/b
                  t2 = t1
                endif
                z1 = -vz(i)*t1 + z(i) - zcent
                z2 = -vz(i)*t2 + z(i) - zcent
              else
c               z1 = LARGEPOS
c               z2 = LARGEPOS
                cycle
              endif
              if (zmin(ii) <= z1 .and. z1 <= zmax(ii)) then
                if ((tt > t1 .and. t1 > 0. .and. tt > 0.) .or.
     &              (t1 > tt .and. tt < 0.)) then
                  tt = t1
                  xx = -vx(i)*t1 + x(i) - xcent
                  yy = -vy(i)*t1 + y(i) - ycent
                  itheta(i) = atan2((zsrf(ii+1) - zsrf(ii)),
     &                              (rsrf(ii) - rsrf(ii+1)))
                  iphi(i) = atan2(yy,xx)
                endif
              endif
              if (zmin(ii) <= z2 .and. z2 <= zmax(ii)) then
                if ((tt > t2 .and. t2 > 0. .and. tt > 0.) .or.
     &              (t2 > tt .and. tt < 0.)) then
                  tt = t2
                  xx = -vx(i)*t2 + x(i) - xcent
                  yy = -vy(i)*t2 + y(i) - ycent
                  itheta(i) = atan2((zsrf(ii+1) - zsrf(ii)),
     &                              (rsrf(ii) - rsrf(ii+1)))
                  iphi(i) = atan2(yy,xx)
                endif
              endif

            endif

          else ! --- Circular line segment

c           --- Use the routine already written for tori.
            ttemp = LARGEPOS
            call ztorus_intercept(rc(ii),abs(rad(ii)),thmin(ii),thmax(ii),
     &                            xcent,ycent,zcent+zc(ii),
     &                            1,x(i),y(i),z(i),vx(i),vy(i),vz(i),
     &                            ttemp,a1,a2)

c           --- If the point is the best so far, then save the data.
c           --- Note that ztorus_intercept checks if the intercept is
c           --- within thmin and thmax.
            if ((tt > ttemp .and. ttemp > 0. .and. tt > 0.) .or.
     &          (ttemp > tt .and. tt < 0.)) then
              tt = ttemp
              itheta(i) = a1
              iphi(i) = a2
            endif

          endif
        enddo

        if (tt .ne. LARGEPOS) then
          xi(i) = -vx(i)*tt + x(i)
          yi(i) = -vy(i)*tt + y(i)
          zi(i) = -vz(i)*tt + z(i)

c         --- Make sure that the normal is oriented to point outward.
c         --- This is a relatively robust, but somewhat expensive, method
c         --- and is independent of the direction of the velocity.
c         --- Take a point along the normal a little bit away
c         --- from the surface and check if it is inside or not.
c         --- A "little bit" is a small fraction of the total
c         --- size of the object. This method will still fail in
c         --- the case when the little bit crosses a small gap in
c         --- the object and the point ends up inside the object
c         --- even though the normal is pointing in the correct
c         --- direction. This won't be fixed since it is very unlikely.
          s = 1.e-10*sqrt((maxval(rmax)-minval(rmin))**2 + 
     &                    (maxval(zmax)-minval(zmin))**2)
          xx = xi(i) + s*sin(itheta(i))*cos(iphi(i))
          yy = yi(i) + s*sin(itheta(i))*sin(iphi(i))
          zz = zi(i) + s*cos(itheta(i))
          call ZSrfrvConductorD(nn,rsrf,zsrf,rad,rc,zc,
     &                          xcent,ycent,zcent,1,xx,yy,zz,dd)
          if (dd < 0.) then
c           --- dd is inside, so change itheta by pi
            if (itheta(i) > 0.) then
              itheta(i) = itheta(i) - pi
            else if (itheta(i) < 0.) then
              itheta(i) = itheta(i) + pi
            endif
          endif

        endif

      enddo

      deallocate(rmin,rmax)
      deallocate(zmin,zmax)
      deallocate(thmin,thmax)
      deallocate(th1,th2)

      return
      end
c=======================================================================
c=======================================================================
c=======================================================================
      subroutine ZSrfrvOutConductorF(lrofzfunc,zmin,zmax,rmax,griddz,
     &                               xcent,ycent,zcent,n,x,y,z,
     &                               delmx,delpx,delmy,delpy,delmz,delpz,
     &                               fuzz)
      logical(ISZ):: lrofzfunc
      real(kind=8):: zmin,zmax,rmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

      integer(ISZ):: i,iterrf,ii
      real(kind=8):: srfrv_f
      real(kind=8):: zzprev1,zzprev2
      real(kind=8):: xx,yy,zz,rr,rsrf
      real(kind=8):: risrf,xisrf,yisrf,ximax,yimax
      real(kind=8):: zzz,zlower,zhigher,rlower,rhigher,rfuzz
      real(kind=8):: ddz,znew,rnew,zi(2),ri(2)

c     --- This guarantees that zzprev is never equal to (z(1)-zcent).
      zzprev1 = 1.5*abs(z(1) - zcent) + 1.
      zzprev2 = 1.5*abs(z(1) - zcent) + 1.
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
          if (zz .ne. zzprev1) risrf = srfrv_f(zz,lrofzfunc,0,0)
          zzprev1 = zz
          xisrf = sqrt(max(0.,risrf**2 - yy**2))
          yisrf = sqrt(max(0.,risrf**2 - xx**2))
          ximax = sqrt(max(0.,rmax**2 - yy**2))
          yimax = sqrt(max(0.,rmax**2 - xx**2))
          if (rr <= risrf) then
            delmx(i) = xisrf + xx
            delpx(i) = xisrf - xx
            delmy(i) = yisrf + yy
            delpy(i) = yisrf - yy
          else if (rr < rmax) then
            if (yy > risrf .or. -risrf > yy .or. risrf == 0.) then
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
            if (xx > risrf .or. -risrf > xx .or. risrf == 0.) then
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

          if (zz .ne. zzprev2) then
            zzprev2 = zz
            zzz = min(max(zz,zmin),zmax)
            zlower = max(zz-griddz,zmin)
            zhigher = min(zz+griddz,zmax)

            rsrf = srfrv_f(zzz,lrofzfunc,0,0)
            rlower = srfrv_f(zlower,lrofzfunc,0,0)
            rhigher = srfrv_f(zhigher,lrofzfunc,0,0)
c           --- The commented code has not been tested.
c           if (.not. lrofzfunc .and. lsrfindrextremum) then
c             rlower = linrzrextremum(zlower,zzz,rsrf,npnts_sr,
c    &                                z_sr,r_sr,rad_sr,zc_sr,rc_sr,+1)
c             rhigher = linrzrextremum(zzz,zhigher,rsrf,npnts_sr,
c    &                                 z_sr,r_sr,rad_sr,zc_sr,rc_sr,+1)
c           else
              rlower = srfrv_f(zlower,lrofzfunc,0,0)
              rhigher = srfrv_f(zhigher,lrofzfunc,0,0)
c           endif
            rfuzz = max(rsrf,rlower,rhigher)*1.e-12
          endif

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
              ri(ii) = srfrv_f(zi(ii),lrofzfunc,0,0) - rr
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
                rnew = srfrv_f(znew,lrofzfunc,0,0)
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
              ri(ii) = srfrv_f(zi(ii),lrofzfunc,0,0) - rr
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
                rnew = srfrv_f(znew,lrofzfunc,0,0)
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
      subroutine ZSrfrvOutConductorD(lrofzfunc,zmin,zmax,rmax,griddz,
     &                               xcent,ycent,zcent,n,x,y,z,distance)
      logical(ISZ):: lrofzfunc
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
      real(kind=8):: zzprev,xx,yy,zz,rr,rsrf,zd,rd
      real(kind=8):: srfrv_f

c     --- This guarantees that zzprev is never equal to (z(1)-zcent).
      zzprev = 1.5*abs(z(1) - zcent) + 1.
      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        rr = sqrt(xx**2 + yy**2)
        zz = z(i) - zcent

        if (zmin <= zz .and. zz <= zmax) then
          if (zz .ne. zzprev) rsrf = srfrv_f(zz,lrofzfunc,0,0)
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
          if (zz .ne. zzprev) rsrf = srfrv_f(zmin,lrofzfunc,0,0)
          if (rr >= rsrf) then
            distance(i) = zmin - zz
          else
            distance(i) = sqrt((zmin-zz)**2 + (rr-rsrf)**2)
          endif

        else
          if (zz .ne. zzprev) rsrf = srfrv_f(zmax,lrofzfunc,0,0)
          if (rr >= rsrf) then
            distance(i) = zz - zmax
          else
            distance(i) = sqrt((zz-zmax)**2 + (rr-rsrf)**2)
          endif

        endif

        zzprev = zz
      enddo

      return
      end
c=======================================================================
      subroutine zsrfrv_intercept(lrofzfunc,zmin,zmax,xcent,ycent,zcent,inside,
     &                            npnts_sr,z_sr,r_sr,rad_sr,zc_sr,rc_sr,
     &                            n,x,y,z,vx,vy,vz,tt,itheta,iphi)
      use Constant
      logical(ISZ):: lrofzfunc
      real(kind=8):: zmin,zmax,xcent,ycent,zcent
      integer(ISZ):: inside
      integer(ISZ):: npnts_sr
      real(kind=8):: z_sr(npnts_sr),r_sr(npnts_sr)
      real(kind=8):: rad_sr(npnts_sr-1),zc_sr(npnts_sr-1),rc_sr(npnts_sr-1)
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: tt(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.
c When inside is == +1, conductor is inside the surface, -1 outside.
c This is only needed to get the correct orientation of the surface normal.

      integer(ISZ):: i,ii,il,niter
      real(kind=8):: vsq,vv,dd,xx,yy,zz,rr,rrsq,t1,t2,z1,z2,ttemp,a1,a2
      real(kind=8):: a,b,c,s,p
      real(kind=8):: r_zmin,r_zmax
      real(kind=8):: t(2),d(2)
      real(kind=8):: srfrv_f

      if (.not. lrofzfunc) then

        do i=1,n

          vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
          if (vsq == 0.) cycle

c         --- Loop over the line segments and calculate the intersection of
c         --- each with the trajectory. Pick the one where the intersect
c         --- is within the segment and the intercept time is smallest.
          do il=1,npnts_sr-1
            if (rad_sr(il) == LARGEPOS) then
              if (z_sr(il+1) == z_sr(il)) then
c               --- Special case of a z plane, where there is only one
c               --- possible intercept. If vz is zero, then there is no
c               --- intercept, so skip.
                if (vz(i) == 0.) cycle
                t1 = ((z(i)-zcent) - z_sr(il))/vz(i)
                xx = -vx(i)*t1 + x(i) - xcent
                yy = -vy(i)*t1 + y(i) - ycent
                rr = sqrt(xx**2 + yy**2)
                if (rr >= min(r_sr(il),r_sr(il+1)) .and.
     &              rr <= max(r_sr(il),r_sr(il+1))) then
                  if ((tt(i) > t1 .and. t1 > 0. .and. tt(i) > 0.) .or.
     &                (t1 > tt(i) .and. tt(i) < 0.)) then
                    tt(i) = t1
                    itheta(i) = atan2((z_sr(il+1) - z_sr(il))*inside,
     &                                (r_sr(il) - r_sr(il+1))*inside)
                    iphi(i) = atan2(yy*inside,xx*inside)
                  endif
                endif
              else
                s = (r_sr(il+1) - r_sr(il))/(z_sr(il+1) - z_sr(il))
                p = r_sr(il) - s*z_sr(il)
                a = (s*vz(i))**2 - vx(i)**2 - vy(i)**2
                b = 2*(-s*vz(i)*(s*(z(i)-zcent)+p) +
     &              vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
                c = (s*(z(i)-zcent)+p)**2 - (x(i)-xcent)**2 - (y(i)-ycent)**2
                if ((b**2 - 4*a*c) > 0.) then
                  if (a .ne. 0.) then
                    t1 = (-b - sqrt(b**2 - 4*a*c))/(2*a)
                    t2 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
                  else
                    t1 = -c/b
                    t2 = t1
                  endif
                  z1 = -vz(i)*t1 + z(i) - zcent
                  z2 = -vz(i)*t2 + z(i) - zcent
                else
c                 z1 = LARGEPOS
c                 z2 = LARGEPOS
                  cycle
                endif
                if (z_sr(il) <= z1 .and. z1 <= z_sr(il+1)) then
                  if ((tt(i) > t1 .and. t1 > 0. .and. tt(i) > 0.) .or.
     &                (t1 > tt(i) .and. tt(i) < 0.)) then
                    tt(i) = t1
                    xx = -vx(i)*t1 + x(i) - xcent
                    yy = -vy(i)*t1 + y(i) - ycent
                    rr = sqrt(xx**2 + yy**2)
c                   --- These angles need to be verified.
                    itheta(i) = atan2((z_sr(il+1) - z_sr(il))*inside,
     &                                (r_sr(il) - r_sr(il+1))*inside)
                    iphi(i) = atan2(yy*inside,xx*inside)
                  endif
                endif
                if (z_sr(il) <= z2 .and. z2 <= z_sr(il+1)) then
                  if ((tt(i) > t2 .and. t2 > 0. .and. tt(i) > 0.) .or.
     &                (t2 > tt(i) .and. tt(i) < 0.)) then
                    tt(i) = t2
                    xx = -vx(i)*t2 + x(i) - xcent
                    yy = -vy(i)*t2 + y(i) - ycent
                    rr = sqrt(xx**2 + yy**2)
c                   --- These angles need to be verified.
                    itheta(i) = atan2((z_sr(il+1) - z_sr(il))*inside,
     &                                (r_sr(il) - r_sr(il+1))*inside)
                    iphi(i) = atan2(yy*inside,xx*inside)
                  endif
                endif
              endif
            else

c             --- Use the routine already written for tori.
              ttemp = LARGEPOS
c             ---- Note that thmin and thmax should be calculated.
              call ztorus_intercept(rc_sr(il),abs(rad_sr(il)),0.,2.*pi,
     &                              xcent,ycent,zcent+zc_sr(il),
     &                              1,x(i),y(i),z(i),vx(i),vy(i),vz(i),
     &                              ttemp,a1,a2)

c             --- Check if the intersection point is within the line segment.
              xx = -vx(i)*ttemp + x(i) - xcent
              yy = -vy(i)*ttemp + y(i) - ycent
              zz = -vz(i)*ttemp + z(i) - ycent
              rr = sqrt(xx**2 + yy**2)
              if (z_sr(il) <= zz .and. zz <= z_sr(il+1) .and.
     &            ((rad_sr(il) > 0. .and. rr > rc_sr(il)) .or.
     &             (rad_sr(il) < 0. .and. rr < rc_sr(il)))) then
c               --- If it is, and the point is the best so far, then save
c               --- the data.
                if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &              (ttemp > tt(i) .and. tt(i) < 0.)) then
                  tt(i) = ttemp
                  if (inside*rad_sr(il) < 0.) a1 = a1 + pi
                  itheta(i) = a1
                  iphi(i) = a2
                endif
              endif

            endif
          enddo

        enddo

      else

c       --- WARNING!! This code is not yet complete and almost certainly
c       --- will not work!!!
        do i=1,n

c         --- Set parameter to initial guess.
          tt(i) = LARGEPOS
          itheta(i) = 0.
          iphi(i) = 0.

          vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
          if (vsq == 0.) cycle
          vv = sqrt(vsq)

c         --- Do a search to find the intersection.

c         --- Find the shortest distance of the point to the surface
c         call ZSrfrvOutConductorD(lrofzfunc,zmin,zmax,rmax,1.,
c    &                             xcent,ycent,zcent,1,x(i),y(i),z(i),dd)

          t(1) = 0.
          d(1) = dd

          t(2) = dd/vv
          xx = -vx(i)*t(2) + x(i)
          yy = -vy(i)*t(2) + y(i)
          zz = -vz(i)*t(2) + z(i)
c         call ZSrfrvOutConductorD(lrofzfunc,zmin,zmax,rmax,1.,
c    &                             xcent,ycent,zcent,1,xx,yy,zz,dd)
          d(2) = dd

          niter = 0
          if (d(1) == 0. .or. d(2) == 0) niter = 12
          ii = 1
          do while ((abs(d(1))+abs(d(2))) > 1.e-8 .and. niter < 12)
            ii = 1
            if (abs(d(2)) > abs(d(1))) ii = 2
            if (d(1) == d(2)) exit
            t(ii) = (t(1)*d(2)-t(2)*d(1))/(d(2)-d(1))
            xx = -vx(i)*t(ii) + x(i)
            yy = -vy(i)*t(ii) + y(i)
            zz = -vz(i)*t(ii) + z(i)
            zz = max(zmin,min(zmax,zz))
c           call ZSrfrvOutConductorD(lrofzfunc,zmin,zmax,rmax,griddz,
c    &                               xcent,ycent,zcent,1,xx,yy,zz,d(ii))
            niter = niter + 1
          enddo
          if (abs(d(2))+abs(d(1)) < 1.e-8) then
            tt(i) = t(ii)
c           --- !!!!!!! not correct !!!!!!
            itheta(i) = 0.
            iphi(i) = 0.
          endif

        enddo

      endif

      return
      end
c=======================================================================
      subroutine ZSrfrvOutIntercept(lrofzfunc,zmin,zmax,rmax,griddz,
     &                              xcent,ycent,zcent,
     &                              n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      use Surface_of_Rev
      logical(ISZ):: lrofzfunc
      real(kind=8):: zmin,zmax,rmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: tt(n)
      real(kind=8):: vsq,dd,xx,yy,zz,rr,rrsq,t1,t2,ttemp
      real(kind=8):: a,b,c,s,p
      real(kind=8):: r_zmin,r_zmax
      real(kind=8):: srfrv_f

      r_zmin = srfrv_f(zmin,lrofzfunc,0,0)
      r_zmax = srfrv_f(zmax,lrofzfunc,0,0)

      tt = LARGEPOS
      xi = LARGEPOS
      yi = LARGEPOS
      zi = LARGEPOS
      itheta = 0.
      iphi = 0.

      call zsrfrv_intercept(lrofzfunc,zmin,zmax,xcent,ycent,zcent,-1,
     &                      npnts_sr,z_sr,r_sr,rad_sr,zc_sr,rc_sr,
     &                      n,x,y,z,vx,vy,vz,tt,itheta,iphi)

      do i=1,n

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

c       --- First, check if trajectory hits conductor end.
        if (vz(i) .ne. 0.) then
          ttemp = -(zmin - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = xx**2 + yy**2
          if (r_zmin**2 <= rrsq .and. rrsq <= rmax**2) then
            if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &          (ttemp > tt(i) .and. tt(i) < 0.)) then
              tt(i) = ttemp
              itheta(i) = pi
              iphi(i) = 0.
            endif
          endif

          ttemp = -(zmax - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = xx**2 + yy**2
          if (r_zmax**2 <= rrsq .and. rrsq <= rmax**2) then
            if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &          (ttemp > tt(i) .and. tt(i) < 0.)) then
              tt(i) = ttemp
              itheta(i) = 0.
              iphi(i) = 0.
            endif
          endif
        endif

c       --- Check if it hits the outer cylinder.
        a = vx(i)**2 + vy(i)**2
        if (a > 0.) then
          b = -2*(vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
          c = (x(i)-xcent)**2 + (y(i)-ycent)**2 - rmax**2
          if ((b**2 - 4*a*c) > 0.) then
            t1 = (-b-sqrt(b**2 - 4*a*c))/(2*a)
            t2 = (-b+sqrt(b**2 - 4*a*c))/(2*a)
            if (t1 > 0. .and. t2 > 0.) then
              ttemp = min(t1,t2)
            else
              ttemp = max(t1,t2)
            endif
            zz = -vz(i)*ttemp + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &            (ttemp > tt(i) .and. tt(i) < 0.)) then
                tt(i) = ttemp
                xx = -vx(i)*ttemp + x(i) - xcent
                yy = -vy(i)*ttemp + y(i) - ycent
                rr = sqrt(xx**2 + yy**2)
                itheta(i) = pi*0.5
                iphi(i) = atan2(yy,xx)
              endif
            endif
          endif
        endif

        if (tt(i) .ne. LARGEPOS) then
          xi(i) = -vx(i)*tt(i) + x(i)
          yi(i) = -vy(i)*tt(i) + y(i)
          zi(i) = -vz(i)*tt(i) + z(i)
        endif

      enddo

      return
      end
c=======================================================================
c=======================================================================
      subroutine ZSrfrvInConductorF(lrofzfunc,zmin,zmax,rmin,griddz,
     &                              xcent,ycent,zcent,n,x,y,z,
     &                              delmx,delpx,delmy,delpy,delmz,delpz,
     &                              fuzz)
      logical(ISZ):: lrofzfunc
      real(kind=8):: zmin,zmax,rmin,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

      integer(ISZ):: i,iterrf,ii
      real(kind=8):: srfrv_f
      real(kind=8):: zzprev1,zzprev2
      real(kind=8):: xx,yy,zz,rr,rsrf
      real(kind=8):: risrf,xisrf,yisrf,ximin,yimin
      real(kind=8):: zzz,zlower,zhigher,rlower,rhigher,rfuzz
      real(kind=8):: ddz,znew,rnew,zi(2),ri(2)

c     --- This guarantees that zzprev is never equal to (z(1)-zcent).
      zzprev1 = 1.5*abs(z(1) - zcent) + 1.
      zzprev2 = 1.5*abs(z(1) - zcent) + 1.
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
          if (zz .ne. zzprev1) risrf = srfrv_f(zz,lrofzfunc,0,0)
          zzprev1 = zz
          xisrf = sqrt(max(0.,risrf**2 - yy**2))
          yisrf = sqrt(max(0.,risrf**2 - xx**2))
          ximin = sqrt(max(0.,rmin**2 - yy**2))
          yimin = sqrt(max(0.,rmin**2 - xx**2))
          if (rr <= rmin) then
            delmx(i) = ximin + xx
            delpx(i) = ximin - xx
            delmy(i) = yimin + yy
            delpy(i) = yimin - yy
          else if (rr < risrf) then
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
            if (abs(yy) <= risrf) then
              if (xx >= 0.) delmx(i) = + xx - xisrf
              if (xx <= 0.) delpx(i) = - xisrf - xx
            endif
            if (abs(xx) <= risrf) then
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

          if (zz .ne. zzprev2) then
            zzprev2 = zz
            zzz = min(max(zz,zmin),zmax)
            zlower = max(zz-griddz,zmin)
            zhigher = min(zz+griddz,zmax)

            rsrf = srfrv_f(zzz,lrofzfunc,0,0)
c           --- The commented code has not been tested.
c           if (.not. lrofzfunc .and. lsrfindrextremum) then
c             rlower = linrzrextremum(zlower,zzz,rsrf,npnts_sr,
c    &                                z_sr,r_sr,rad_sr,zc_sr,rc_sr,-1)
c             rhigher = linrzrextremum(zzz,zhigher,rsrf,npnts_sr,
c    &                                 z_sr,r_sr,rad_sr,zc_sr,rc_sr,-1)
c           else
              rlower = srfrv_f(zlower,lrofzfunc,0,0)
              rhigher = srfrv_f(zhigher,lrofzfunc,0,0)
c           endif
            rfuzz = max(rsrf,rlower,rhigher)*1.e-12
          endif

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
              ri(ii) = srfrv_f(zi(ii),lrofzfunc,0,0) - rr
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
                rnew = srfrv_f(znew,lrofzfunc,0,0)
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
              ri(ii) = srfrv_f(zi(ii),lrofzfunc,0,0) - rr
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
                rnew = srfrv_f(znew,lrofzfunc,0,0)
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
      subroutine ZSrfrvInConductorD(lrofzfunc,zmin,zmax,rmin,griddz,
     &                              xcent,ycent,zcent,n,x,y,z,distance)
      logical(ISZ):: lrofzfunc
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
      real(kind=8):: zzprev,xx,yy,zz,rr,rsrf,zd,rd
      real(kind=8):: srfrv_f

c     --- This guarantees that zzprev is never equal to (z(1)-zcent).
      zzprev = 1.5*abs(z(1) - zcent) + 1.
      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        if (zmin <= zz .and. zz <= zmax) then
          if (zz .ne. zzprev) rsrf = srfrv_f(zz,lrofzfunc,0,0)
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
          if (zz .ne. zzprev) rsrf = srfrv_f(zmin,lrofzfunc,0,0)
          if (rr <= rsrf) then
            distance(i) = zmin - zz
          else
            distance(i) = sqrt((zmin-zz)**2 + (rr-rsrf)**2)
          endif

        else
          if (zz .ne. zzprev) rsrf = srfrv_f(zmax,lrofzfunc,0,0)
          if (rr <= rsrf) then
            distance(i) = zz - zmax
          else
            distance(i) = sqrt((zz-zmax)**2 + (rr-rsrf)**2)
          endif

        endif

        zzprev = zz
      enddo

      return
      end
c=======================================================================
      subroutine ZSrfrvInIntercept(lrofzfunc,zmin,zmax,rmin,griddz,
     &                             xcent,ycent,zcent,
     &                             n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      use Surface_of_Rev
      logical(ISZ):: lrofzfunc
      real(kind=8):: zmin,zmax,rmin,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: tt(n)
      real(kind=8):: vsq,dd,xx,yy,zz,rr,rrsq,t1,t2,ttemp
      real(kind=8):: a,b,c,s,p
      real(kind=8):: r_zmin,r_zmax
      real(kind=8):: srfrv_f

      r_zmin = srfrv_f(zmin,lrofzfunc,0,0)
      r_zmax = srfrv_f(zmax,lrofzfunc,0,0)

      tt = LARGEPOS
      xi = LARGEPOS
      yi = LARGEPOS
      zi = LARGEPOS
      itheta = 0.
      iphi = 0.

      call zsrfrv_intercept(lrofzfunc,zmin,zmax,xcent,ycent,zcent,+1,
     &                      npnts_sr,z_sr,r_sr,rad_sr,zc_sr,rc_sr,
     &                      n,x,y,z,vx,vy,vz,tt,itheta,iphi)

      do i=1,n

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

c       --- First, check if trajectory hits conductor end.
        if (vz(i) .ne. 0.) then
          ttemp = -(zmin - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = xx**2 + yy**2
          if (r_zmin**2 >= rrsq .and. rrsq >= rmin**2) then
            if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &          (ttemp > tt(i) .and. tt(i) < 0.)) then
              tt(i) = ttemp
              itheta(i) = pi
              iphi(i) = 0.
            endif
          endif

          ttemp = -(zmax - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = xx**2 + yy**2
          if (r_zmax**2 >= rrsq .and. rrsq >= rmin**2) then
            if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &          (ttemp > tt(i) .and. tt(i) < 0.)) then
              tt(i) = ttemp
              itheta(i) = 0.
              iphi(i) = 0.
            endif
          endif
        endif

c       --- Check if it hits the inner cylinder.
        a = vx(i)**2 + vy(i)**2
        if (a > 0.) then
          b = -2*(vx(i)*(x(i)-xcent) + vy(i)*(y(i)-ycent))
          c = (x(i)-xcent)**2 + (y(i)-ycent)**2 - rmin**2
          if ((b**2 - 4*a*c) > 0.) then
            t1 = (-b-sqrt(b**2 - 4*a*c))/(2*a)
            t2 = (-b+sqrt(b**2 - 4*a*c))/(2*a)
            if (t1 > 0. .and. t2 > 0.) then
              ttemp = min(t1,t2)
            else
              ttemp = max(t1,t2)
            endif
            zz = -vz(i)*ttemp + z(i)
            if (zmin <= zz .and. zz <= zmax) then
              if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &            (ttemp > tt(i) .and. tt(i) < 0.)) then
                tt(i) = ttemp
                xx = -vx(i)*ttemp + x(i) - xcent
                yy = -vy(i)*ttemp + y(i) - ycent
                rr = sqrt(xx**2 + yy**2)
                itheta(i) = pi*0.5
                iphi(i) = atan2(-yy,-xx)
              endif
            endif
          endif
        endif

        if (tt(i) .ne. LARGEPOS) then
          xi(i) = -vx(i)*tt(i) + x(i)
          yi(i) = -vy(i)*tt(i) + y(i)
          zi(i) = -vz(i)*tt(i) + z(i)
        endif

      enddo

      return
      end
c=======================================================================
      real(kind=8) function linrzrextremum(z1,z2,rr,npnts_sr,z_sr,r_sr,
     &                                     rad_sr,zc_sr,rc_sr,rsign)
      real(kind=8):: z1,z2,rr
      integer(ISZ):: npnts_sr
      real(kind=8):: z_sr(npnts_sr),r_sr(npnts_sr)
      real(kind=8):: rad_sr(npnts_sr-1),zc_sr(npnts_sr-1),rc_sr(npnts_sr-1)
      integer(ISZ):: rsign

c This searches through the line segments to find the maximum radius
c within the range between z1 and z2. When rsign is positive, it finds
c the max, when negative the min.

      integer(ISZ):: ii
      real(kind=8):: rextremum,r1

c     --- Return the input value if the data is completely outside the range.
      if (z2 < z_sr(1) .or. z1 > z_sr(npnts_sr) .or. npnts_sr == 1 .or.
     &    z1 == z2) then
        linrzrextremum = rr
        return
      endif

c     --- Find the line segment which covers the range of z where the point is.
c     --- Find the line segment end just to the left of the data point.
c     --- If z1 if less that z_sr(1), then ii will be 1. If z1 greater than
c     --- the end of the last line segment, then ii will be npnts_sr-1.
      do ii=1,npnts_sr-2
        if (z1 < z_sr(ii+1)) exit
      enddo

c     --- Set a sensible starting value
      rextremum = rr

c     --- Search segment by segment until they go beyond the end of the region
c     --- or there are no more segments.
      do while (z_sr(ii) <= z2 .and. ii < npnts_sr)

c       --- Get the radius at the left end point of the segment and compare
c       --- it to the rextremum.
        if (z_sr(ii) < z1) then
c         --- The segment starts below z1, so find the radius at z1.
          if (rad_sr(ii) == LARGEPOS) then
c           --- straight line
            r1 = r_sr(ii) + (r_sr(ii+1) - r_sr(ii))*(z1 - z_sr(ii))/
     &                                              (z_sr(ii+1) - z_sr(ii))
          else
c           --- arc
            if (rad_sr(ii) > 0) then
              r1 = rc_sr(ii)
     &             + sqrt(max(0.,rad_sr(ii)**2 - (z1 - zc_sr(ii))**2))
            else
              r1 = rc_sr(ii)
     &             - sqrt(max(0.,rad_sr(ii)**2 - (z1 - zc_sr(ii))**2))
            endif
          endif
        else
c         --- The segment starts within the range so just use the radius
c         --- of the left end.
          r1 = r_sr(ii)
        endif
c       --- Compare the radius at the left end to rextremum.
        if (rsign > 0) rextremum = max(rextremum,r1)
        if (rsign < 0) rextremum = min(rextremum,r1)

c       --- If there is an arc, check if there is an extremum within the
c       --- part of the segment within the z range. This is only done if
c       --- the center of the circle is within the range and if the sign
c       --- of the radius is the same as the rsign flag.
        if (rad_sr(ii) < LARGEPOS .and.
     &      max(z1,z_sr(ii)) <= zc_sr(ii) .and.
     &      zc_sr(ii) <= min(z2,z_sr(ii+1)) .and.
     &      rsign*rad_sr(ii) > 0.) then
          r1 = rc_sr(ii) + rad_sr(ii)
          if (rsign > 0) rextremum = max(rextremum,r1)
          if (rsign < 0) rextremum = min(rextremum,r1)
        endif

c       --- Advance to the next segment
        ii = ii + 1

      enddo

c     --- Now, check the end point of the last segment
      if (z_sr(ii) >= z2) then
c       --- It extends beyond the range, so get the radius at the end of
c       --- the range.
        if (rad_sr(ii) == LARGEPOS) then
          r1 = r_sr(ii) + (r_sr(ii) - r_sr(ii-1))*(z_sr(ii) - z2)/
     &                                            (z_sr(ii) - z_sr(ii-1))
        else
          if (rad_sr(ii-1) > 0) then
            r1 = rc_sr(ii-1)
     &           + sqrt(max(0.,rad_sr(ii-1)**2 - (z2 - zc_sr(ii-1))**2))
          else
            r1 = rc_sr(ii-1)
     &           - sqrt(max(0.,rad_sr(ii-1)**2 - (z2 - zc_sr(ii-1))**2))
          endif
        endif
      else
c       --- Otherwise, use the radius at the end of the segment.
        r1 = r_sr(ii)
      endif
c     --- Compare the radius at the right end to rextremum.
      if (rsign > 0) rextremum = max(rextremum,r1)
      if (rsign < 0) rextremum = min(rextremum,r1)

      linrzrextremum = rextremum

      return
      end
c=======================================================================
      real(kind=8) function linrzintersection(zz,rr,npnts_sr,z_sr,r_sr,
     &                                        rad_sr,zc_sr,rc_sr,zend,istep)
      real(kind=8):: zz,rr
      integer(ISZ):: npnts_sr
      real(kind=8):: z_sr(npnts_sr),r_sr(npnts_sr)
      real(kind=8):: rad_sr(npnts_sr-1),zc_sr(npnts_sr-1),rc_sr(npnts_sr-1)
      real(kind=8):: zend
      integer(ISZ):: istep

c This searches through the line segments to find the one that the z-line
c will intersect. From that, the location of the z intersection can be
c calculated directly.

      integer(ISZ):: ii
      real(kind=8):: r1,r2,zi,z1
      logical(ISZ):: done

c     --- Find the line segment which covers the range of z where the point is.
c     --- The search can either be to the left or to the right. But in both
c     --- cases, find the line segment end just to the left of the data point.
c     --- If zz if less that z_sr(1), then ii will be 1. If zz greater than
c     --- the end of the last line segment, then ii will be npnts_sr-1.
      do ii=1,npnts_sr-2
        if (zz < z_sr(ii+1)) exit
      enddo

c     --- Search up (or down) in the line segments to find one which crosses
c     --- the r of the point.
      done = .false.
      do while (.not. done)
        r1 = min(r_sr(ii+1),r_sr(ii))
        r2 = max(r_sr(ii+1),r_sr(ii))
        if (rad_sr(ii) < LARGEPOS .and.
     &      (z_sr(ii) < zc_sr(ii) .and. zc_sr(ii) < z_sr(ii+1))) then
c         --- If the segment is an arc and if the z center of the arc is
c         --- within the z range, then the r max or min of the segment will
c         --- be somewhere along the arc rather than at an endpoint.
          if (rad_sr(ii) < 0.) then
            r1 = rc_sr(ii) + rad_sr(ii)
          else
            r2 = rc_sr(ii) + rad_sr(ii)
          endif
        endif
        if (r1 <= rr .and. rr <= r2) then
c         --- Directly calculate the z intersection for this segment.
          if (rad_sr(ii) == LARGEPOS) then
c           --- For a straight line
            if (r_sr(ii+1) - r_sr(ii) == 0.) then
c             --- A special case is needed if two adjacent radii are equal and
c             --- if the radius of the point is at that same radius.
              if (zz < z_sr(ii)) then
c               --- This will only ever be done if there is only one line
c               --- segment and the z of the left end is greater than the z of
c               --- the left end of the conductor (i.e. the data does not span
c               --- the conductor).
                zi = z_sr(ii)
              else if (zz > z_sr(ii+1)) then
c               --- This will only ever be done if there is only one line
c               --- segment and the z of the right end is less than the z of
c               --- the right end of the conductor (i.e. the data does not span
c               --- the conductor).
                zi = z_sr(ii+1)
              else
c               --- In all other cases, the z of the point will be within
c               --- the range of the line segment and so the point is on
c               --- the surface of the conductor.
                zi = zz
              endif
            else
              zi = (rr - r_sr(ii))*(z_sr(ii+1) - z_sr(ii))/
     &             (r_sr(ii+1) - r_sr(ii)) + z_sr(ii)
            endif
          else
c           --- For an arc
            z1 = sqrt(max(0.,rad_sr(ii)**2 - (rr-rc_sr(ii))**2))
c           --- The line along z intersects the arc at the two points zc+-z1.
c           --- Pick first the one seen first along the direction of search.
            zi = zc_sr(ii) + istep*z1
            if ((zi > zz .and. istep == 1).or.(zi < zz .and. istep == -1)) then
c             --- If that one is not in the correct direction relative to the
c             --- point, then pick the other.
              zi = zc_sr(ii) - istep*z1
            endif
          endif
c         --- Only accept this segment if the intersection is on the correct
c         --- side of the point. Otherwise, try the next segment.
          if ((zz > zi .and. istep == 1) .or. (zz < zi .and. istep == -1)) then
            done = .true.
            exit
          endif
        endif
        ii = ii - istep
c       --- Check if the search has gone beyond the extent of the data.
c       --- This will only occur if the line segments do not span the conductor.
c       --- If this does happen, return the end point of the conductor.
        if (ii == 0 .or. ii == npnts_sr) then
          zi = zend
          done = .true.
        endif
      enddo

      linrzintersection = zi

      return
      end
c=======================================================================
      subroutine ZSrfrvInOutConductorF(lrminofz,lrmaxofz,zmin,zmax,griddz,
     &                                 xcent,ycent,zcent,n,x,y,z,
     &                                 delmx,delpx,delmy,delpy,delmz,delpz,fuzz)
      use Surface_of_Rev
      logical(ISZ):: lrminofz,lrmaxofz
      real(kind=8):: zmin,zmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n)
      real(kind=8):: delmx(n),delpx(n)
      real(kind=8):: delmy(n),delpy(n)
      real(kind=8):: delmz(n),delpz(n)
      real(kind=8):: fuzz

      integer(ISZ):: i,iterrf,ii
      real(kind=8):: srfrv_f
      real(kind=8):: zzprev1,zzprev2
      real(kind=8):: xx,yy,zz,rr,rimin,rimax,rmin,rmax
      real(kind=8):: ximin,yimin,ximax,yimax
      real(kind=8):: zzz,zlower,zhigher,rfuzz
      real(kind=8):: rminlower,rminhigher,rmaxlower,rmaxhigher
      real(kind=8):: ddz,znew,rnew,zi(2),ri(2),r1,r2
      real(kind=8):: delmz1,delpz1
      real(kind=8):: linrzintersection,linrzrextremum

c     --- This guarantees that zzprev is never equal to (z(1)-zcent).
      zzprev1 = 1.5*abs(z(1) - zcent) + 1.
      zzprev2 = 1.5*abs(z(1) - zcent) + 1.
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
          if (zz .ne. zzprev1) then
            zzprev1 = zz
            rimin = min(srfrv_f(zz,lrminofz,1,-1),srfrv_f(zz,lrminofz,1,+1))
            rimax = max(srfrv_f(zz,lrmaxofz,2,-1),srfrv_f(zz,lrmaxofz,2,+1))
          endif
          ximin = sqrt(max(0.,rimin**2 - yy**2))
          yimin = sqrt(max(0.,rimin**2 - xx**2))
          ximax = sqrt(max(0.,rimax**2 - yy**2))
          yimax = sqrt(max(0.,rimax**2 - xx**2))
          if (rr < rimin) then
            delmx(i) = ximin + xx
            delpx(i) = ximin - xx
            delmy(i) = yimin + yy
            delpy(i) = yimin - yy
          else if (rr < rimax) then
            if (yy > rimin .or. -rimin > yy .or. rimin == 0.) then
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
            if (xx > rimin .or. -rimin > xx .or. rimin == 0.) then
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
            if (abs(yy) <= rimax) then
              if (xx >= 0.) delmx(i) = + xx - ximax
              if (xx <= 0.) delpx(i) = - ximax - xx
            endif
            if (abs(xx) <= rimax) then
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

          if (zz .ne. zzprev2) then
            zzprev2 = zz
            zzz = min(max(zz,zmin),zmax)
            zlower = max(zz-griddz,zmin)
            zhigher = min(zz+griddz,zmax)

            rmin = srfrv_f(zzz,lrminofz,1,0)
            rmax = srfrv_f(zzz,lrmaxofz,2,0)
            if (.not. lrminofz .and. lsrfindrextremum) then
c             --- XXX This code needs further debuggins. In the second call,
c             --- zzz+griddz should be zhigher, but when that change is made,
c             --- the code breaks and the conductors come out incorrectly - a
c             --- spike forms at the zhigher location. The code works as is for
c             --- the case of interest however and so is left as is. Also, note
c             --- that lsrfindrextremum defaults to false so this code isn't
c             --- used in general, since it will be slower.
              rminlower = linrzrextremum(zlower,zzz,rmin,npnts_srmin,
     &                           z_srmin,r_srmin,rad_srmin,zc_srmin,rc_srmin,-1)
              rminhigher = linrzrextremum(zzz,zzz+griddz,rmin,npnts_srmin,
     &                           z_srmin,r_srmin,rad_srmin,zc_srmin,rc_srmin,-1)
            else
              rminlower = srfrv_f(zlower,lrminofz,1,-1)
              rminhigher = srfrv_f(zhigher,lrminofz,1,+1)
            endif
            if (.not. lrmaxofz .and. lsrfindrextremum) then
              rmaxlower = linrzrextremum(zlower,zzz,rmax,npnts_srmax,
     &                           z_srmax,r_srmax,rad_srmax,zc_srmax,rc_srmax,+1)
              rmaxhigher = linrzrextremum(zzz,zzz+griddz,rmax,npnts_srmax,
     &                           z_srmax,r_srmax,rad_srmax,zc_srmax,rc_srmax,+1)
            else
              rmaxlower = srfrv_f(zlower,lrmaxofz,2,-1)
              rmaxhigher = srfrv_f(zhigher,lrmaxofz,2,+1)
            endif
            rfuzz = max(rmin,rmax,rminlower,rminhigher,rmaxlower,rmaxhigher)*1.e-12
          endif

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
              if (.not. lrminofz) then
                znew = linrzintersection(zz,rr,npnts_srmin,z_srmin,r_srmin,
     &                                   rad_srmin,zc_srmin,rc_srmin,zmin,1)
                delmz(i) = zz - znew
              else
c               --- If grid location is between rmin and rminlower find
c               --- intersection point.
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
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
                  ri(ii) = srfrv_f(zi(ii),lrminofz,1,0) - rr
                  iterrf = iterrf + 1
                end do
                if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &              zlower <= zi(ii) .and. zi(ii) < zz) then
                  delmz(i) = zz - zi(ii)
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (min(zz,zmax) - zlower)*0.5
                  znew = zlower
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,lrminofz,1,0)
                    if (rnew > rr) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew <= rr) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delmz(i) = zz - znew
                endif
              endif
              if (rr >= rmin) delmz(i) = -delmz(i)
            endif
            if (min(rmaxlower,rmax) <= rr .and. rr <= max(rmaxlower,rmax)) then
              if (.not. lrmaxofz) then
                znew = linrzintersection(zz,rr,npnts_srmax,z_srmax,r_srmax,
     &                                   rad_srmax,zc_srmax,rc_srmax,zmin,1)
                delmz1 = zz - znew
              else
c               --- If grid location is between rmax and rmaxlower find
c               --- intersection point.
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
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
                  ri(ii) = srfrv_f(zi(ii),lrmaxofz,2,0) - rr
                  iterrf = iterrf + 1
                end do
                if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &              zlower <= zi(ii) .and. zi(ii) < zz) then
                  delmz1 = zz - zi(ii)
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (min(zz,zmax) - zlower)*0.5
                  znew = zlower
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,lrmaxofz,2,0)
                    if (rnew > rr) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew <= rr) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delmz1 = zz - znew
                endif
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
              if (.not. lrminofz) then
                znew = linrzintersection(zz,rr,npnts_srmin,z_srmin,r_srmin,
     &                                   rad_srmin,zc_srmin,rc_srmin,zmax,-1)
                delpz(i) = znew - zz
              else
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
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
                  ri(ii) = srfrv_f(zi(ii),lrminofz,1,0) - rr
                  iterrf = iterrf + 1
                end do
                if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &              zz < zi(ii) .and. zi(ii) <= zhigher) then
                  delpz(i) = zi(ii) - zz
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (max(zz,zmin) - zhigher)*0.5
                  znew = zhigher
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,lrminofz,1,0)
                    if (rnew <= rr) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew > rr) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delpz(i) = znew - zz
                endif
              endif
              if (rr >= rmin) delpz(i) = -delpz(i)
            endif
            if (min(rmaxhigher,rmax) <= rr .and.
     &          rr <= max(rmaxhigher,rmax)) then
              if (.not. lrmaxofz) then
                znew = linrzintersection(zz,rr,npnts_srmax,z_srmax,r_srmax,
     &                                   rad_srmax,zc_srmax,rc_srmax,zmax,-1)
                delpz1 = znew - zz
              else
c               --- First, try Regula-Falsi (secant method) iteration to
c               --- find distance in z to the surface. If that doesn't
c               --- converge rapidly, then go to the brute force binary search.
c               --- Empirically, it was found that stopping at 12 iterations
c               --- gives the best performance, with a broad minimum around
c               --- that point.
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
                  ri(ii) = srfrv_f(zi(ii),lrmaxofz,2,0) - rr
                  iterrf = iterrf + 1
                end do
                if (abs(ri(2))+abs(ri(1)) < rfuzz .and.
     &              zz < zi(ii) .and. zi(ii) <= zhigher) then
                  delpz1 = zi(ii) - zz
                else
c                 --- A simple binary search is used to find distance in z
c                 --- to the surface. This is used since it is more robust
c                 --- than other root finding methods. The major problem
c                 --- occurs when the radius changes like a step function.
c                 --- The secant method (the original coding) converged
c                 --- very slowly in that case, leaving errors as large as
c                 --- 10% in the delta and slowing down the code.
                  ddz = (max(zz,zmin) - zhigher)*0.5
                  znew = zhigher
                  do iterrf=1,32
                    znew = znew + ddz
                    rnew = srfrv_f(znew,lrmaxofz,2,0)
                    if (rnew <= rr) then
                      ddz = -abs(ddz)*0.5
                    else if (rnew > rr) then
                      ddz = +abs(ddz)*0.5
                    endif
                  enddo
                  delpz1 = znew - zz
                endif
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
      subroutine ZSrfrvInOutConductorD(lrminofz,lrmaxofz,zmin,zmax,griddz,
     &                                 xcent,ycent,zcent,n,x,y,z,distance)
      logical(ISZ):: lrminofz,lrmaxofz
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
      real(kind=8):: zzprev,xx,yy,zz,rr,rmin,rmax,zd,rd
      real(kind=8):: srfrv_f

c     --- This guarantees that zzprev is never equal to (z(1)-zcent).
      zzprev = 1.5*abs(z(1) - zcent) + 1.
      do i=1,n

        xx = x(i) - xcent
        yy = y(i) - ycent
        zz = z(i) - zcent
        rr = sqrt(xx**2 + yy**2)

        if (zmin <= zz .and. zz <= zmax) then
          if (zz .ne. zzprev) then
            rmin = srfrv_f(zz,lrminofz,1,0)
            rmax = srfrv_f(zz,lrmaxofz,2,0)
          endif
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
          if (zz .ne. zzprev) then
            rmin = srfrv_f(zmin,lrminofz,1,0)
            rmax = srfrv_f(zmin,lrmaxofz,2,0)
          endif
          if (rr < rmin) then
            distance(i) = sqrt((zmin-zz)**2 + (rr-rmin)**2)
          else if (rr <= rmax) then
            distance(i) = zmin - zz
          else
            distance(i) = sqrt((zmin-zz)**2 + (rr-rmax)**2)
          endif

        else
          if (zz .ne. zzprev) then
            rmin = srfrv_f(zmax,lrminofz,1,0)
            rmax = srfrv_f(zmax,lrmaxofz,2,0)
          endif
          if (rr < rmin) then
            distance(i) = sqrt((zz-zmax)**2 + (rr-rmin)**2)
          else if (rr <= rmax) then
            distance(i) = zz - zmax
          else
            distance(i) = sqrt((zz-zmax)**2 + (rr-rmax)**2)
          endif

        endif
        zzprev = zz
      enddo

      return
      end
c=======================================================================
      subroutine ZSrfrvInOutIntercept(lrminofz,lrmaxofz,zmin,zmax,griddz,
     &                                xcent,ycent,zcent,
     &                                n,x,y,z,vx,vy,vz,xi,yi,zi,itheta,iphi)
      use Constant
      use Surface_of_Rev
      logical(ISZ):: lrminofz,lrmaxofz
      real(kind=8):: zmin,zmax,griddz,xcent,ycent,zcent
      integer(ISZ):: n
      real(kind=8):: x(n),y(n),z(n),vx(n),vy(n),vz(n)
      real(kind=8):: xi(n),yi(n),zi(n),itheta(n),iphi(n)

c Finds place where particle inside intercepted the surface and the
c normal angle of the surface at that point.

      integer(ISZ):: i
      real(kind=8):: tt(n)
      real(kind=8):: vsq,dd,xx,yy,zz,rr,rrsq,t1,t2,ttemp
      real(kind=8):: a,b,c,s,p
      real(kind=8):: rmin_zmin,rmin_zmax
      real(kind=8):: rmax_zmin,rmax_zmax
      real(kind=8):: srfrv_f

      rmin_zmin = srfrv_f(zmin,lrminofz,1,0)
      rmin_zmax = srfrv_f(zmax,lrminofz,1,0)
      rmax_zmin = srfrv_f(zmin,lrmaxofz,2,0)
      rmax_zmax = srfrv_f(zmax,lrmaxofz,2,0)

      tt = LARGEPOS
      xi = LARGEPOS
      yi = LARGEPOS
      zi = LARGEPOS
      itheta = 0.
      iphi = 0.

      call zsrfrv_intercept(lrminofz,zmin,zmax,xcent,ycent,zcent,-1,
     &                      npnts_srmin,z_srmin,r_srmin,
     &                      rad_srmin,zc_srmin,rc_srmin,
     &                      n,x,y,z,vx,vy,vz,tt,itheta,iphi)
      call zsrfrv_intercept(lrmaxofz,zmin,zmax,xcent,ycent,zcent,+1,
     &                      npnts_srmax,z_srmax,r_srmax,
     &                      rad_srmax,zc_srmax,rc_srmax,
     &                      n,x,y,z,vx,vy,vz,tt,itheta,iphi)

      do i=1,n

        vsq = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq == 0.) cycle

c       --- First, check if trajectory hits conductor end.
        if (vz(i) .ne. 0.) then
          ttemp = -(zmin - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = xx**2 + yy**2
          if (rmin_zmin**2 <= rrsq .and. rrsq <= rmax_zmin**2) then
            if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &          (ttemp > tt(i) .and. tt(i) < 0.)) then
              tt(i) = ttemp
              itheta(i) = pi
              iphi(i) = 0.
            endif
          endif

          ttemp = -(zmax - z(i))/vz(i)
          xx = -vx(i)*ttemp + x(i)
          yy = -vy(i)*ttemp + y(i)
          rrsq = xx**2 + yy**2
          if (rmin_zmax**2 <= rrsq .and. rrsq <= rmax_zmax**2) then
            if ((tt(i) > ttemp .and. ttemp > 0. .and. tt(i) > 0.) .or.
     &          (ttemp > tt(i) .and. tt(i) < 0.)) then
              tt(i) = ttemp
              itheta(i) = 0.
              iphi(i) = 0.
            endif
          endif
        endif

        if (tt(i) .ne. LARGEPOS) then
          xi(i) = -vx(i)*tt(i) + x(i)
          yi(i) = -vy(i)*tt(i) + y(i)
          zi(i) = -vz(i)*tt(i) + z(i)
        endif

      enddo

      return
      end
c========================================================================
      real(kind=8) function srfrv_f(zz,lrofzfunc,icase,izflag)
      use Machine_param
      use Surface_of_Rev
      real(kind=8):: zz
      logical(ISZ):: lrofzfunc
      integer(ISZ):: icase,izflag

c     This function returns the value of radius as a function of z.
c     That can be calculated either from a piecewise-linear curve, or from
c     the Python routine lrofzfunc.
c     The radius returned is always greater than or equal to zero.

      integer(ISZ):: i,ii,iz
      real(kind=8):: z_sriim1

      if (lrofzfunc) then

c       --- Get r(z) from the Python function.
        srfrv_z = zz
        if (icase == 0) then
          call callpythonfunc("rofzfunc","generateconductors")
        else if (icase == 1) then
          call callpythonfunc("rminofz","generateconductors")
        else if (icase == 2) then
          call callpythonfunc("rmaxofz","generateconductors")
        endif

      else
        if (icase == 0) then

c       --- Calculate r(z) from a piecewise-linear curve.

c       --- Find which piece z is in.
          ii = 0
          do i=1,npnts_sr
            if (zz >= z_sr(i)) ii = i
          enddo

c         --- Calculate r from z.
          if (ii == 0) then
            srfrv_r = r_sr(1)
          elseif (ii == npnts_sr) then
            srfrv_r = r_sr(npnts_sr)
          else
c           --- Get z_sr at ii-1 to avoid accessing z_sr(0)
            z_sriim1 = 0.
            if (ii > 1) z_sriim1 = z_sr(ii-1)
            if (ii > 1 .and.
     &          z_sriim1 == z_sr(ii) .and. zz == z_sr(ii) .and.
     &          izflag == 1) then
              if (ii-1 == 1) then
                srfrv_r = r_sr(1)
              else
                if (rad_sr(ii-2) == largepos) then
                  srfrv_r = r_sr(ii-2) + (zz - z_sr(ii-1))*
     &                         (r_sr(ii-2) - r_sr(ii-1))/(z_sr(ii-2) - z_sr(ii-1))
                else
                  if (rad_sr(ii-2) > 0) then
                    srfrv_r = rc_sr(ii-2)
     &                   + sqrt(max(0.,rad_sr(ii-2)**2 - (zz - zc_sr(ii-2))**2))
                  else
                    srfrv_r = rc_sr(ii-2)
     &                   - sqrt(max(0.,rad_sr(ii-2)**2 - (zz - zc_sr(ii-2))**2))
                  endif
                endif
              endif
            else
              if (rad_sr(ii) == largepos) then
                srfrv_r = r_sr(ii) + (zz - z_sr(ii))*
     &                       (r_sr(ii+1) - r_sr(ii))/(z_sr(ii+1) - z_sr(ii))
              else
                if (rad_sr(ii) > 0) then
                  srfrv_r = rc_sr(ii)
     &                   + sqrt(max(0.,rad_sr(ii)**2 - (zz - zc_sr(ii))**2))
                else
                  srfrv_r = rc_sr(ii)
     &                   - sqrt(max(0.,rad_sr(ii)**2 - (zz - zc_sr(ii))**2))
                endif
              endif
            endif
          endif

        elseif (icase == 1) then

c         --- Calculate rmin(z) from a piecewise-linear curve.
c         --- Depeding on which direction is being looked, the order in which
c         --- the points are searched changes. This only really matters when
c         --- the conductor description includes step function in r(z), i.e.
c         --- where two neighboring points have the same z value.

          if (izflag <= 0) then

c           --- Find which piece z is in.
            ii = 0
            do i=1,npnts_srmin
              if (zz >= z_srmin(i)) ii = i
            enddo

c           --- Calculate r from z.
            if (ii == 0) then
              srfrv_r = r_srmin(1)
            elseif (ii == npnts_srmin) then
              srfrv_r = r_srmin(npnts_srmin)
            else
              if (rad_srmin(ii) == largepos) then
                srfrv_r = r_srmin(ii) + (zz - z_srmin(ii))*
     &                (r_srmin(ii+1) - r_srmin(ii))/(z_srmin(ii+1) - z_srmin(ii))
              else
                if (rad_srmin(ii) > 0) then
                  srfrv_r = rc_srmin(ii)
     &                  + sqrt(max(0.,rad_srmin(ii)**2 - (zz - zc_srmin(ii))**2))
                else
                  srfrv_r = rc_srmin(ii)
     &                  - sqrt(max(0.,rad_srmin(ii)**2 - (zz - zc_srmin(ii))**2))
                endif
              endif
            endif

          else

c           --- Find which piece z is in.
            ii = npnts_srmin+1
            do i=npnts_srmin,1,-1
              if (zz <= z_srmin(i)) ii = i
            enddo

c           --- Calculate r from z.
            if (ii == npnts_srmin+1) then
              srfrv_r = r_srmin(npnts_srmin)
            elseif (ii == 1) then
              srfrv_r = r_srmin(1)
            else
              if (rad_srmin(ii-1) == largepos) then
                srfrv_r = r_srmin(ii) + (zz - z_srmin(ii))*
     &           (r_srmin(ii-1) - r_srmin(ii))/(z_srmin(ii-1) - z_srmin(ii))
              else
                if (rad_srmin(ii-1) > 0) then
                  srfrv_r = rc_srmin(ii-1)
     &              + sqrt(max(0.,rad_srmin(ii-1)**2 - (zz - zc_srmin(ii-1))**2))
                else
                  srfrv_r = rc_srmin(ii-1)
     &              - sqrt(max(0.,rad_srmin(ii-1)**2 - (zz - zc_srmin(ii-1))**2))
                endif
              endif
            endif

          endif

        elseif (icase == 2) then

c         --- Calculate rmax(z) from a piecewise-linear curve.
c         --- Depeding on which direction is being looked, the order in which
c         --- the points are searched changes. This only really matters when
c         --- the conductor description includes step function in r(z), i.e.
c         --- where two neighboring points have the same z value.

          if (izflag <= 0) then

c           --- Find which piece z is in.
            ii = 0
            do i=1,npnts_srmax
              if (zz >= z_srmax(i)) ii = i
            enddo

c           --- Calculate r from z.
            if (ii == 0) then
              srfrv_r = r_srmax(1)
            elseif (ii == npnts_srmax) then
              srfrv_r = r_srmax(npnts_srmax)
            else
              if (rad_srmax(ii) == largepos) then
                srfrv_r = r_srmax(ii) + (zz - z_srmax(ii))*
     &                (r_srmax(ii+1) - r_srmax(ii))/(z_srmax(ii+1) - z_srmax(ii))
              else
                if (rad_srmax(ii) > 0) then
                  srfrv_r = rc_srmax(ii)
     &                  + sqrt(max(0.,rad_srmax(ii)**2 - (zz - zc_srmax(ii))**2))
                else
                  srfrv_r = rc_srmax(ii)
     &                  - sqrt(max(0.,rad_srmax(ii)**2 - (zz - zc_srmax(ii))**2))
                endif
              endif
            endif

          else

c           --- Find which piece z is in.
            ii = npnts_srmax+1
            do i=npnts_srmax,1,-1
              if (zz <= z_srmax(i)) ii = i
            enddo

c           --- Calculate r from z.
            if (ii == npnts_srmax+1) then
              srfrv_r = r_srmax(npnts_srmax)
            elseif (ii == 1) then
              srfrv_r = r_srmax(1)
            else
              if (rad_srmax(ii-1) == largepos) then
                srfrv_r = r_srmax(ii) + (zz - z_srmax(ii))*
     &           (r_srmax(ii-1) - r_srmax(ii))/(z_srmax(ii-1) - z_srmax(ii))
              else
                if (rad_srmax(ii-1) > 0) then
                  srfrv_r = rc_srmax(ii-1)
     &              + sqrt(max(0.,rad_srmax(ii-1)**2 - (zz - zc_srmax(ii-1))**2))
                else
                  srfrv_r = rc_srmax(ii-1)
     &              - sqrt(max(0.,rad_srmax(ii-1)**2 - (zz - zc_srmax(ii-1))**2))
                endif
              endif
            endif

          endif

        endif

      endif

      srfrv_f = max(0.,srfrv_r)

      return
      end
c=======================================================================
      subroutine solvequartic(a0,a1,a2,a3,x1,x2,x3,x4)
      real(kind=8):: a0,a1,a2,a3
c
c     subroutine to solve for the double precision roots of a quartic
c     equation of the form x**4+a3*x**3+a2*x**2+a1*x+a0.  method is
c     from abromowitz and stegun p. 17.
c
c
      real(kind=8):: zero,one,two,b0,b1,b2,r1,r0,p1,p0
      real(kind=8):: sq1,sq0,arg,c1m,c1p,c0m,c0p
      integer(ISZ):: noreal,noroot
      complex*16 x1,x2,x3,x4,u1,u2,u3,ci
      one=1.d0
      zero=0.d0
      two=2.d0
      ci=cmplx(zero,one)
c     calculate the coefficients of the auxilliary cubic
      b0=-(a1**2+a0*a3**2-4.d0*a0*a2)
      b1=a1*a3-4.d0*a0
      b2=-a2
c      calculate the roots of the cubic
      call solvecubic(b0,b1,b2,u1,u2,u3,noreal)
c     calculate the coefficients of the auxilliary quadratic
      no root=0
   10 r1=(a3/two)**2+u1-a2
      r0=(u1/2)**2-a0
      no root=no root+1
c     if roots are all real seek the root which makes the coefficients real
      if(noreal.eq.1) go to 20
      p1=r1
      p0=r0
      if(p1.ge.zero.and.p0.ge.zero) go to 20
      u1=u2
      if(noroot.eq.1) go to 10
      u1=u3
      if(noroot.eq.2) go to 10
      write(6,60)
   60 format('   no real quadratics found')
   20 continue
      if(r0.lt.zero) go to 100
      if(r1.lt.zero) go to 200
   30 sq0=dsqrt(r0)
      sq1=dsqrt(r1)
      arg=a3*u1/two-a1
      if(arg.ne.zero) sq1=dsign(sq1,arg)
      c1m=a3/two-sq1
      c1p=a3/two+sq1
      c0m=u1/two-sq0
      c0p=u1/two+sq0
c      calculate the actual roots from the quadratics
      call solvequadratic(c0m,c1m,x1,x2)
      call solvequadratic(c0p,c1p,x3,x4)
      return
  100 write(6,65)
      r0=zero
   65 format('   r0 set to zero')
      go to 30
  200 write(6,70)
   70 format('   r1 set to zero')
      r1=zero
      go to 30
      end
      subroutine solvequadratic(a0,a1,x1,x2)
      real(kind=8):: a0,a1
      complex*16 x1,x2,rt,ci
      real(kind=8):: zero,two,r
      zero=0.d0
      two=2.d0
      ci=cmplx(0.d0,1.d0)
      r=a1**2-4.d0*a0
      if (r.ge.zero) rt=dsqrt(r)
      if(r.lt.zero) rt=ci*dsqrt(-r)
      x1=(-a1+rt)/two
      x2=(-a1-rt)/two
      return
      end
      subroutine solvecubic(a0,a1,a2,x1,x2,x3,noreal)
      real(kind=8):: a0,a1,a2
      integer(ISZ):: noreal
c
c     subroutine to find double precision roots of a cubic equation
c     of the form  x**3+a2*x**2+a1*x+a0=0.  the roots are returned in the
c     complex numbers x1,x2,x3.  from abromiwitz and stegun p. 17.
c     noreal is the number of real roots.
c
c
      real(kind=8):: one,two,three,third,root3,sixth,a2d3,q,r,t,roott,arg
      real(kind=8):: s1,s2,re,ai,phi,amp
      complex*16 x1,x2,x3,ci
      one=1.d0
      two=2.d0
      three=3.d0
      third=one/three
      root3=dsqrt(three)
      third=one/three
      sixth=one/6.d0
      ci=cmplx(0.d0,1.d0)
      a2d3=a2*third
      q=a1/3.d0-a2**2/9.d0
      r=(a1*a2-3.d0*a0)/6.d0-a2**3/27.d0
c
c     test to find character of roots
c
      t=q**3+r**2
      no real=1
      if(t.le.0.d0) no real=3
      if(t.lt.0.d0) go to 100
c
c     complex roots or degenerate real roots
c
      roott=dsqrt(t)
      arg=r+roott
      if(arg.ge.0.d0) s1=arg**third
      if(arg.lt.0.d0) s1=-((-arg)**third)
      arg=r-roott
      if(arg.ge.0.d0) s2=arg**third
      if(arg.lt.0.d0) s2=-((-arg)**third)
      re=s1+s2
      ai=root3/two*(s1-s2)
      x1=re-a2d3
      x2=-re/two-a2d3+ci*ai
      x3=-re/two-a2d3-ci*ai
      return
c
c     three real roots
c
  100 phi=datan2(dsqrt(-t),r)*third
      amp=(r**2-t)**sixth
      re=amp*dcos(phi)
      ai=amp*dsin(phi)
      x1=two*re-a2d3
      x2=-re-a2d3-root3*ai
      x3=-re-a2d3+root3*ai
      return
      end

