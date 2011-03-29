#include "top.h"
c=============================================================================
c Routines for handing the Forward semi-Lagrangian method
c=============================================================================
      subroutine applyboundaryconditionsonfgrid2d2v(nvx,nvz,nx,nz,fgrid,pbounds)
      use Subtimerstop
      integer(ISZ):: nvx,nvz,nx,nz
      real(kind=8):: fgrid(-1:nvx+1,-1:nvz+1,-1:nx+1,-1:nz+1)
      integer(ISZ):: pbounds(0:5)

      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Always zero out the data at the extrema in velocity
      fgrid(-1,:,:,:) = 0.
      fgrid(nvx+1,:,:,:) = 0.
      fgrid(:,-1,:,:) = 0.
      fgrid(:,nvz+1,:,:) = 0.

      if (pbounds(0) == 0) fgrid(:,:,-1,0:nz) = 0.
      if (pbounds(1) == 0) fgrid(:,:,nx+1,0:nz) = 0.
      if (pbounds(0) == 1) then
        fgrid(:,:,-1,0:nz) = fgrid(:,:,1,0:nz)
        fgrid(:,:,0,0:nz) = 2.*fgrid(:,:,0,0:nz)
      endif
      if (pbounds(1) == 1) then
        fgrid(:,:,nx+1,0:nz) = fgrid(:,:,nx-1,0:nz)
        fgrid(:,:,nx,0:nz) = 2.*fgrid(:,:,nx,0:nz)
      endif
      if (pbounds(0) == 2 .and. pbounds(1) == 2) then
        fgrid(:,:,nx,0:nz) = fgrid(:,:,nx,0:nz) + fgrid(:,:,0,0:nz)
        fgrid(:,:,-1:0,0:nz) = fgrid(:,:,nx-1:nx,0:nz)
        fgrid(:,:,nx+1,0:nz) = fgrid(:,:,1,0:nz)
      endif

      if (pbounds(4) == 0) fgrid(:,:,:,-1) = 0.
      if (pbounds(5) == 0) fgrid(:,:,:,nz+1) = 0.
      if (pbounds(4) == 1) then
        fgrid(:,:,:,-1) = fgrid(:,:,:,1)
        fgrid(:,:,:,0) = 2.*fgrid(:,:,:,0)
      endif
      if (pbounds(5) == 1) then
        fgrid(:,:,:,nz+1) = fgrid(:,:,:,nz-1)
        fgrid(:,:,:,nz) = 2.*fgrid(:,:,:,nz)
      endif
      if (pbounds(4) == 2 .and. pbounds(5) == 2) then
        fgrid(:,:,:,nz) = fgrid(:,:,:,nz) + fgrid(:,:,:,0)
        fgrid(:,:,:,-1:0) = fgrid(:,:,:,nz-1:nz)
        fgrid(:,:,:,nz+1) = fgrid(:,:,:,1)
      endif

!$OMP MASTER
      if (ltoptimesubs) timeapplyboundaryconditionsonfgrid2d2v =
     &                  timeapplyboundaryconditionsonfgrid2d2v +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine applyminmaxboundaryconditions2d(nx,nz,v,pbounds,lmax)
      use Subtimerstop
      integer(ISZ):: nx,nz
      real(kind=8):: v(-1:nx+1,-1:nz+1)
      integer(ISZ):: pbounds(0:5)
      logical(ISZ):: lmax

      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (pbounds(0) == 0 .or. pbounds(0) == 1) v(-1,0:nz) = v(1,0:nz)
      if (pbounds(1) == 0 .or. pbounds(1) == 1) v(nx+1,0:nz) = v(nx-1,0:nz)
      if (pbounds(0) == 2 .and. pbounds(1) == 2) then
        v(-1,0:nz) = v(nx-1,0:nz)
        if (lmax) then
          v(nx,0:nz) = max(v(nx,0:nz),v(0,0:nz))
        else
          v(nx,0:nz) = min(v(nx,0:nz),v(0,0:nz))
        endif
        v(0,0:nz) = v(nx,0:nz)
        v(nx+1,0:nz) = v(1,0:nz)
      endif

      if (pbounds(4) == 0 .or. pbounds(4) == 1) v(:,-1) = v(:,1)
      if (pbounds(5) == 0 .or. pbounds(5) == 1) v(:,nz+1) = v(:,nz-1)
      if (pbounds(4) == 2 .and. pbounds(5) == 2) then
        v(:,-1) = v(:,nz-1)
        if (lmax) then
          v(:,nz) = max(v(:,nz),v(:,0))
        else
          v(:,nz) = min(v(:,nz),v(:,0))
        endif
        v(:,0) = v(:,nz)
        v(:,nz+1) = v(:,1)
      endif

!$OMP MASTER
      if (ltoptimesubs) timeapplyminmaxboundaryconditions2d =
     &                  timeapplyminmaxboundaryconditions2d +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine createparticlesfromfgrid(pgroup,js,geometry,zbeam,vbeam,
     &                                    pbounds)
      use Subtimerstop
      use ForwardSemiLagrange
      use ParticleGroupmodule
      use Particles,Only: wpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: js,geometry
      real(kind=8):: zbeam,vbeam
      integer(ISZ):: pbounds(0:5)

      real(kind=8),pointer:: xp(:),yp(:),zp(:),uxp(:),uyp(:),uzp(:)
      real(kind=8):: vybeam,vzbeam
      integer(ISZ):: npnew,ip
      integer(ISZ):: ix,iz,ivx,ivz
      integer(ISZ):: ix0,iz0
      integer(ISZ):: ix1,iz1
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (wpid == 0) then
        call kaboom("createparticlesfromfgrid: wpid must be set")
        return
      endif

      if (pbounds(0) == 0) then
        ix0 = 1
      else
        ix0 = 0
      endif
      if (pbounds(1) == 0 .or. pbounds(1) == 2) then
        ix1 = nx_fsl - 1
      else
        ix1 = nx_fsl
      endif

      if (pbounds(4) == 0) then
        iz0 = 1
      else
        iz0 = 0
      endif
      if (pbounds(5) == 0 .or. pbounds(5) == 2) then
        iz1 = nz_fsl - 1
      else
        iz1 = nz_fsl
      endif

c     --- Count how many particles will be created, the number of fgrid
c     --- cells that are positive.
      npnew = 0
      do iz=iz0,iz1
        do ix=ix0,ix1
          do ivz=0,nvz_fsl
            do ivx=0,nvx_fsl
              if (fgrid_fsl(ivx,ivz,ix,iz) > 0.) npnew = npnew + 1
            enddo
          enddo
        enddo
      enddo
      if (npnew == 0) return

      if (npnew > pgroup%nps(js+1)) then
        pgroup%nps(js+1) = 0
        call checkparticlegroup(pgroup,js+1,0,npnew)
      endif

      if (geometry == 3) then
c       -- XZgeom
        xp => pgroup%xp
        yp => pgroup%yp
        zp => pgroup%zp
        uxp => pgroup%uxp
        uyp => pgroup%uyp
        uzp => pgroup%uzp
        vybeam = 0.
        vzbeam = vbeam
      else if (geometry == 4) then
c       -- XYgeom
        xp => pgroup%xp
        yp => pgroup%zp
        zp => pgroup%yp
        uxp => pgroup%uxp
        uyp => pgroup%uzp
        uzp => pgroup%uyp
        vybeam = vbeam
        vzbeam = 0.
      endif

c     --- Regenerate the weighted particles at the grid locations
      ip = pgroup%ins(js+1)
      do iz=iz0,iz1
        do ix=ix0,ix1
          do ivz=0,nvz_fsl
            do ivx=0,nvx_fsl
              if (fgrid_fsl(ivx,ivz,ix,iz) > 0.) then
                xp(ip) = xmin_fsl + ix*dx_fsl
                yp(ip) = 0.
                zp(ip) = zmin_fsl + iz*dz_fsl + zbeam
                uxp(ip) = vxmin_fsl(ix,iz) + ivx*dvx_fsl
                uyp(ip) = vybeam
                uzp(ip) = vzmin_fsl(ix,iz) + ivz*dvz_fsl + vzbeam
                pgroup%gaminv(ip) = 1.
                pgroup%pid(ip,wpid) = fgrid_fsl(ivx,ivz,ix,iz)
                ip = ip + 1
              endif
            enddo
          enddo
        enddo
      enddo
      pgroup%nps(js+1) = npnew

!$OMP MASTER
      if (ltoptimesubs) timecreateparticlesfromfgrid =
     &                  timecreateparticlesfromfgrid +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine enforcepositivity2d(nx,nz,fgrid,threshold)
      use Subtimerstop
      integer(ISZ):: nx,nz
      real(kind=8):: fgrid(-1:nx+1,-1:nz+1)
      real(kind=8):: threshold

c Given a distribution function, fgrid, this routine enforces positivity
c by taking any negative values and distributing it to nearby, already
c positive, locations. If ensures the the integral of fgrid remains fixed.
c The action is local, avoiding a global calculation.
c Much care was needed to handle the odd cases, for example if a
c negative value is surrounded by negative and zero values.

      integer(ISZ):: ix,iz,iparity,ix0,niter,dd,id
      real(kind=8):: sumcap,ferror,fdel
      integer(ISZ):: nitermax
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Repeat the operation until there are no negative locations.
      nitermax = min(nx/2,nz/2)

c     --- Loop over the array using even-odd ordering
      do iparity = 0,1
        do iz = 0,nz
          ix0 = mod(iz + iparity,2)
          do ix = ix0,nx,2

c           --- For each location, iterate until its value is nonzero
c           --- (of maximum number of iterations is reached).
c           --- For each iteration, if sumcap is zero, then increase
c           --- dd so that a larger region will be searched.
            niter = 0
            dd = 1
            do while (fgrid(ix,iz) < threshold .and.
     &                fgrid(ix,iz) .ne. 0. .and.
     &                niter < nitermax)
              niter = niter + 1

c             --- Create the sum of the surrounding, positive values.
c             --- this loops over the square of grid cells around the
c             --- center cell in question, including corners.
              sumcap = 0.
              do id=-dd,dd-1
                if (0 <= ix+id .and. ix+id <= nx .and.
     &              0 <= iz-dd .and. iz-dd <= nz .and.
     &              fgrid(ix+id,iz-dd) > threshold) then
                  sumcap = sumcap + fgrid(ix+id,iz-dd)
                endif
                if (0 <= ix+id+1 .and. ix+id+1 <= nx .and.
     &              0 <= iz+dd .and. iz+dd <= nz .and.
     &              fgrid(ix+id+1,iz+dd) > threshold) then
                  sumcap = sumcap + fgrid(ix+id+1,iz+dd)
                endif
                if (0 <= ix-dd .and. ix-dd <= nx .and.
     &              0 <= iz+id+1 .and. iz+id+1 <= nz .and.
     &              fgrid(ix-dd,iz+id+1) > threshold) then
                  sumcap = sumcap + fgrid(ix-dd,iz+id+1)
                endif
                if (0 <= ix+dd .and. ix+dd <= nx .and.
     &              0 <= iz+id .and. iz+id <= nz .and.
     &              fgrid(ix+dd,iz+id) > threshold) then
                  sumcap = sumcap + fgrid(ix+dd,iz+id)
                endif
              enddo
              if (sumcap > 0.) then
c               --- If the sum is > 0, then the negative f can be
c               --- distributed. The minimum value of the correction
c               --- factor is -1. This prevents a nearby value from
c               --- being forced negative. Becuase of this, the
c               --- negative f will not be fully corrected and another
c               --- iteration will be needed.
c               --- Add and subtract the error from each neighbor.
c               --- Writing the code this way seems to remove roundoff
c               --- problems, and allows the case when ferror is maxed
c               --- out and the correction is incomplete.
                fdel = max(-1.,fgrid(ix,iz)/sumcap)
                do id=-dd,dd-1
                  if (0 <= ix+id .and. ix+id <= nx .and.
     &                0 <= iz-dd .and. iz-dd <= nz .and.
     &                fgrid(ix+id,iz-dd) > threshold) then
                    ferror = fgrid(ix+id,iz-dd)*fdel
                    fgrid(ix+id,iz-dd) = fgrid(ix+id,iz-dd) + ferror
                    fgrid(ix  ,iz) = fgrid(ix  ,iz) - ferror
                  endif
                  if (0 <= ix+id+1 .and. ix+id+1 <= nx .and.
     &                0 <= iz+dd .and. iz+dd <= nz .and.
     &                fgrid(ix+id+1,iz+dd) > threshold) then
                    ferror = fgrid(ix+id+1,iz+dd)*fdel
                    fgrid(ix+id+1,iz+dd) = fgrid(ix+id+1,iz+dd) + ferror
                    fgrid(ix  ,iz) = fgrid(ix  ,iz) - ferror
                  endif
                  if (0 <= ix-dd .and. ix-dd <= nx .and.
     &                0 <= iz+id+1 .and. iz+id+1 <= nz .and.
     &                fgrid(ix-dd,iz+id+1) > threshold) then
                    ferror = fgrid(ix-dd,iz+id+1)*fdel
                    fgrid(ix-dd,iz+id+1) = fgrid(ix-dd,iz+id+1) + ferror
                    fgrid(ix,iz  ) = fgrid(ix,iz  ) - ferror
                  endif
                  if (0 <= ix+dd .and. ix+dd <= nx .and.
     &                0 <= iz+id .and. iz+id <= nz .and.
     &                fgrid(ix+dd,iz+id) > threshold) then
                    ferror = fgrid(ix+dd,iz+id)*fdel
                    fgrid(ix+dd,iz+id) = fgrid(ix+dd,iz+id) + ferror
                    fgrid(ix,iz  ) = fgrid(ix,iz  ) - ferror
                  endif
                enddo
                if (fdel > -1.) fgrid(ix,iz) = 0.
              endif
c             --- Increment the size of the area checked. If the
c             --- correction is incomplete, that means that there was
c             --- not enough positive f to distribute the error too,
c             --- so go to a bigger area.
              dd = dd + 1

            enddo

          enddo
        enddo
      enddo

!$OMP MASTER
      if (ltoptimesubs) timeenforcepositivity2d =
     &                  timeenforcepositivity2d +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine enforcepositivity2d2v(nvx,nvz,nx,nz,fgrid,threshold)
      use Subtimerstop
      integer(ISZ):: nvx,nvz,nx,nz
      real(kind=8):: fgrid(-1:nvx+1,-1:nvz+1,-1:nx+1,-1:nz+1)
      real(kind=8):: threshold

c Given a distribution function, fgrid, this routine enforces positivity
c by taking any negative values and distributing it to nearby, already
c positive, locations. If ensures the the integral of fgrid remains fixed.
c The action is local, avoiding a global calculation.
c Much care was needed to handle the odd cases, for example if a
c negative value is surrounded by negative and zero values.

      integer(ISZ):: ivx,ivz,ix,iz,iparity,ivx0,niter,dd,id
      integer(ISZ):: jvx,jvz,jx,jz
      integer(ISZ):: jvx1,jvz1,jx1,jz1
      integer(ISZ):: jvx2,jvz2,jx2,jz2
      real(kind=8):: sumcap,fdel
      logical(ISZ):: lcappedcorrection
      integer(ISZ):: nitermax
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Only one sweep over the grid is needed, since for each
c     --- location that is negative, the inner loop iterates until its
c     --- data is appropriately spread out.

      nitermax = min(min(nvx/2,nvz/2),min(nx/2,nz/2))

c     --- Loop over the array using even-odd ordering
      do iparity = 0,1
        do iz = 0,nz
          do ix = 0,nx
            do ivz = 0,nvz
              ivx0 = mod(iz + ix + ivz + iparity,2)
              do ivx = ivx0,nvx,2

                if (fgrid(ivx,ivz,ix,iz) < threshold .and.
     &              fgrid(ivx,ivz,ix,iz) .ne. 0.) then

c                 --- For each location, iterate until its value is zero
c                 --- (or maximum number of iterations is reached).
c                 --- For each iteration, if sumcap is zero or there is
c                 --- not enough positive f nearby to do the full
c                 --- correction, then increase dd so that a larger region
c                 --- will be searched.
                  niter = 0
                  dd = 1
                  do 
                    niter = niter + 1
                    if (niter >= nitermax) exit

c                   --- Create the sum of the surrounding, positive values,
c                   --- including all of the grid cells within a square,
c                   --- and including guard cells.
c                   --- Note that the center cell will be excluded
c                   --- since it less than the threshold.
                    jvx1 = max(-1,ivx - dd)
                    jvx2 = min(nvx+1,ivx + dd)
                    jvz1 = max(-1,ivz - dd)
                    jvz2 = min(nvz+1,ivz + dd)
                    jx1 = max(-1,ix - dd)
                    jx2 = min(nx+1,ix + dd)
                    jz1 = max(-1,iz - dd)
                    jz2 = min(nz+1,iz + dd)
                    sumcap = sum(fgrid(jvx1:jvx2,jvz1:jvz2,jx1:jx2,jz1:jz2),
     &               mask=fgrid(jvx1:jvx2,jvz1:jvz2,jx1:jx2,jz1:jz2)>=threshold)

                    if (sumcap > 0. .and. -fgrid(ivx,ivz,ix,iz) <= sumcap) then
c                     --- If the sum is > 0, and there is enough positive
c                     --- f nearby, then the distribution of the negative
c                     --- f can be done. If there is not enough, the loop
c                     --- will iterate and examine a bigger area.
                      fdel = fgrid(ivx,ivz,ix,iz)/sumcap
                      fgrid(ivx,ivz,ix,iz) = 0.
                      where (fgrid(jvx1:jvx2,jvz1:jvz2,jx1:jx2,jz1:jz2) >= threshold)
                        fgrid(jvx1:jvx2,jvz1:jvz2,jx1:jx2,jz1:jz2) =
     &                  fgrid(jvx1:jvx2,jvz1:jvz2,jx1:jx2,jz1:jz2) +
     &                  fgrid(jvx1:jvx2,jvz1:jvz2,jx1:jx2,jz1:jz2)*fdel
                      end where
c                     --- Exit the loop since f has been fully corrected.
                      exit
                    endif

c                   --- Increment the size of the area checked. If the
c                   --- correction is incomplete, that means that there was
c                   --- not enough positive f to distribute the error too,
c                   --- so go to a bigger area.
                    dd = dd + 1

                  enddo

                endif

              enddo
            enddo
          enddo
        enddo
      enddo

!$OMP MASTER
      if (ltoptimesubs) timeenforcepositivity2d2v =
     &                  timeenforcepositivity2d2v +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine findminmaxongrid2d2v(np,x,z,vx,vz,nx,nz,
     &                                xmin,xmax,zmin,zmax,
     &                                vxmin,vxmax,vzmin,vzmax)
      use Subtimerstop
      integer(ISZ):: np,nx,nz
      real(kind=8):: x(np), z(np), vx(np), vz(np)
      real(kind=8):: xmin,xmax,zmin,zmax
      real(kind=8):: vxmin(-1:nx+1,-1:nz+1),vxmax(-1:nx+1,-1:nz+1)
      real(kind=8):: vzmin(-1:nx+1,-1:nz+1),vzmax(-1:nx+1,-1:nz+1)

c Finds the min and max velocities as a function of position

      integer(ISZ):: ip,ix,iz
      real(kind=8):: gx,gz,dxi,dzi
      real(kind=8):: lastvalue1,lastvalue2
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (xmax == xmin .or. zmax == zmin) then
        call kaboom("findminmaxongrid2d2v: the grid extent is zero, xmin==xmax or zmin==ymax")
        return
      endif

      vxmin = +LARGEPOS
      vxmax = -LARGEPOS
      vzmin = +LARGEPOS
      vzmax = -LARGEPOS

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gz = (z(ip) - zmin)*dzi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gz .and. gz < nz) then
          ix = int(gx)
          iz = int(gz)
          vxmin(ix  ,iz  ) = min(vxmin(ix  ,iz  ),vx(ip))
          vxmin(ix+1,iz  ) = min(vxmin(ix+1,iz  ),vx(ip))
          vxmin(ix  ,iz+1) = min(vxmin(ix  ,iz+1),vx(ip))
          vxmin(ix+1,iz+1) = min(vxmin(ix+1,iz+1),vx(ip))
          vxmax(ix  ,iz  ) = max(vxmax(ix  ,iz  ),vx(ip))
          vxmax(ix+1,iz  ) = max(vxmax(ix+1,iz  ),vx(ip))
          vxmax(ix  ,iz+1) = max(vxmax(ix  ,iz+1),vx(ip))
          vxmax(ix+1,iz+1) = max(vxmax(ix+1,iz+1),vx(ip))
          vzmin(ix  ,iz  ) = min(vzmin(ix  ,iz  ),vz(ip))
          vzmin(ix+1,iz  ) = min(vzmin(ix+1,iz  ),vz(ip))
          vzmin(ix  ,iz+1) = min(vzmin(ix  ,iz+1),vz(ip))
          vzmin(ix+1,iz+1) = min(vzmin(ix+1,iz+1),vz(ip))
          vzmax(ix  ,iz  ) = max(vzmax(ix  ,iz  ),vz(ip))
          vzmax(ix+1,iz  ) = max(vzmax(ix+1,iz  ),vz(ip))
          vzmax(ix  ,iz+1) = max(vzmax(ix  ,iz+1),vz(ip))
          vzmax(ix+1,iz+1) = max(vzmax(ix+1,iz+1),vz(ip))

        endif

      enddo

      call cleanuplargepos(vxmin,+1)
      call cleanuplargepos(vxmax,-1)
      call cleanuplargepos(vzmin,+1)
      call cleanuplargepos(vzmax,-1)

!$OMP MASTER
      if (ltoptimesubs) timefindminmaxongrid2d2v =
     &                  timefindminmaxongrid2d2v +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      CONTAINS
        subroutine cleanuplargepos(v,lpsign)
        real(kind=8):: v(-1:nx+1,-1:nz+1)
        integer(ISZ):: lpsign
c       --- Remove the LARGEPOS values, replacing them with nearby values.
c       --- Two sweeps are done, one in each direction, and any LARGEPOS
c       --- are replaced by the previous value found during the sweep.
        do iz=-1,nz+1
          lastvalue1 = lpsign*LARGEPOS
          lastvalue2 = lpsign*LARGEPOS
          do ix=-1,nx+1
            if (v(ix,iz) == lpsign*LARGEPOS) then
              v(ix,iz) = lastvalue1
            else
              lastvalue1 = v(ix,iz)
            endif
            if (v(nx-ix,iz) == lpsign*LARGEPOS) then
              v(nx-ix,iz) = lastvalue2
            else
              lastvalue2 = v(nx-ix,iz)
            endif
          enddo
        enddo
        do ix=-1,nx+1
          lastvalue1 = lpsign*LARGEPOS
          lastvalue2 = lpsign*LARGEPOS
          do iz=-1,nz+1
            if (v(ix,iz) == lpsign*LARGEPOS) then
              v(ix,iz) = lastvalue1
            else
              lastvalue1 = v(ix,iz)
            endif
            if (v(ix,nz-iz) == lpsign*LARGEPOS) then
              v(ix,nz-iz) = lastvalue2
            else
              lastvalue2 = v(ix,nz-iz)
            endif
          enddo
        enddo

        end subroutine cleanuplargepos

      end
c=============================================================================
      subroutine loadfgridguassian(nz,density,radius,
     &                             vthermalperp,vthermalparallel)
      use Subtimerstop
      use Constant
      use ForwardSemiLagrange
      integer(ISZ):: nz
      real(kind=8):: density(0:nz),radius(0:nz)
      real(kind=8):: vthermalperp,vthermalparallel
      integer(ISZ):: ivx,ivz,ix,iz
      real(kind=8):: xx,vx,vz,fvz,fvx,sumf

      real(kind=8):: fvx0,fvz0
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (nz .ne. nz_fsl) then
        call kaboom("loadfgridguassian: nz must be equal to nz_fsl")
        return
      endif

      if (dvx_fsl == 0. .and. nvx_fsl > 0) then
        dvx_fsl = 6.*vthermalperp/nvx_fsl
      else if (nvx_fsl == 0 .and. dvx_fsl > 0.) then
        nvx_fsl = 6.*vthermalperp/dvx_fsl
      else
        call kaboom("loadfgridguassian: either dvx_fsl or nvx_fsl must be nonzero")
        return
      endif
      if (dvz_fsl == 0. .and. nvz_fsl > 0) then
        dvz_fsl = 6.*vthermalparallel/nvz_fsl
      else if (nvz_fsl == 0 .and. dvz_fsl > 0.) then
        nvz_fsl = 6.*vthermalparallel/dvz_fsl
      else
        call kaboom("loadfgridguassian: either dvz_fsl or nvz_fsl must be nonzero")
        return
      endif

      call gchange("ForwardSemiLagrange",0)

      vxmin_fsl = -dvx_fsl*nvx_fsl/2.
      vxmax_fsl = +dvx_fsl*nvx_fsl/2.
      vzmin_fsl = -dvz_fsl*nvz_fsl/2.
      vzmax_fsl = +dvz_fsl*nvz_fsl/2.

      fvz0 = sqrt(1./(2.*pi*vthermalparallel**2))
      fvx0 = sqrt(1./(2.*pi*vthermalperp**2))

      do iz=0,nz_fsl
        if (density(iz) == 0. .or. radius(iz) == 0.) cycle
        do ix=0,nx_fsl
          xx = xmin_fsl + ix*dx_fsl
          if (abs(xx) > radius(iz)) cycle
          do ivz=0,nvz_fsl
            vz = vzmin_fsl(ix,iz) + ivz*dvz_fsl
            fvz = fvz0*exp(-vz**2/(2.*vthermalparallel**2))
            do ivx=0,nvx_fsl
              vx = vxmin_fsl(ix,iz) + ivx*dvx_fsl
              fvx = fvx0*exp(-vx**2/(2.*vthermalperp**2))
              fgrid_fsl(ivx,ivz,ix,iz) = density(iz)*fvx*fvz
            enddo
          enddo
        enddo
      enddo

c     --- Normalize the distribution
c     sumf = sum(fgrid_fsl)
c     if (sumf > 0.) then
c       fgrid_fsl = fgrid_fsl/sumf
c     endif

!$OMP MASTER
      if (ltoptimesubs) timeloadfgridguassian = timeloadfgridguassian +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine remapparticles2d2v(pgroup,geometry,zbeam,l4symtry,l2symtry)
      use Subtimerstop
      use ForwardSemiLagrange
      use ParticleGroupmodule
      use Particles,Only: wpid
      use InGen,Only: pboundxy,pbound0,pboundnz,pbounds
      type(ParticleGroup):: pgroup
      integer(ISZ):: geometry
      real(kind=8):: zbeam
      logical(ISZ):: l4symtry,l2symtry

      real(kind=8),pointer:: xp(:),zp(:),uxp(:),uzp(:)
      integer(ISZ):: js,ip,i1,i2
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (geometry == 3) then
c       -- XZgeom
        xp => pgroup%xp
        zp => pgroup%zp
        uxp => pgroup%uxp
        uzp => pgroup%uzp
      else if (geometry == 4) then
c       -- XYgeom
        xp => pgroup%xp
        zp => pgroup%yp
        uxp => pgroup%uxp
        uzp => pgroup%uyp
      endif

      do js=0,pgroup%ns - 1

        i1 = pgroup%ins(js+1)
        i2 = i1 + pgroup%nps(js+1) - 1

c       --- Find the min and max velocity as a function of position
        call findminmaxongrid2d2v(pgroup%nps(js+1),xp(i1),zp(i1),
     &                            uxp(i1),uzp(i1),nx_fsl,nz_fsl,
     &                            xmin_fsl,xmax_fsl,
     &                            zmin_fsl+zbeam,zmax_fsl+zbeam,
     &                            vxmin_fsl,vxmax_fsl,vzmin_fsl,vzmax_fsl)

c       --- Do the parallel nastiness of sharing boundary data
c       --- This is a serial place holder.
        call setboundsfromflags(pbounds,pboundxy,pbound0,pboundnz,
     &                          l2symtry,l4symtry)
        call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,vxmin_fsl,
     &                                       pbounds,.false.)
        call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,vxmax_fsl,
     &                                       pbounds,.true.)
        call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,vzmin_fsl,
     &                                       pbounds,.false.)
        call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,vzmax_fsl,
     &                                       pbounds,.true.)

c       --- Round the velocity mins and maxs to the grid cells
        vxmin_fsl = floor(vxmin_fsl/dvx_fsl)*dvx_fsl
        vzmin_fsl = floor(vzmin_fsl/dvz_fsl)*dvz_fsl
        vxmax_fsl = ceiling(vxmax_fsl/dvx_fsl)*dvx_fsl
        vzmax_fsl = ceiling(vzmax_fsl/dvz_fsl)*dvz_fsl
c       --- Temporary until there are better diagnostics
        vxmin_fsl = minval(vxmin_fsl)
        vxmax_fsl = maxval(vxmax_fsl)
        vzmin_fsl = minval(vzmin_fsl)
        vzmax_fsl = maxval(vzmax_fsl)

c       --- Resize the grid if needed
        nvx_fsl = max(nvx_fsl,nint(maxval((vxmax_fsl - vxmin_fsl)/dvx_fsl)))
        nvz_fsl = max(nvz_fsl,nint(maxval((vzmax_fsl - vzmin_fsl)/dvz_fsl)))
        call gchange("ForwardSemiLagrange",0)

c       --- Deposit the particles into the distribution grid
c       fgrid_fsl(0:nvx_fsl,0:nvz_fsl,0:nx_fsl,0:nz_fsl) = 0.
        fgrid_fsl = 0.
        call setgrid2d2vmodbspline2w(pgroup%nps(js+1),
     &                               xp(i1),zp(i1),uxp(i1),uzp(i1),
     &                               pgroup%pid(i1,wpid),
     &                               nvx_fsl,nvz_fsl,
     &                               nx_fsl,nz_fsl,fgrid_fsl,
     &                               vxmin_fsl,dvx_fsl,vzmin_fsl,dvz_fsl,
     &                               xmin_fsl,dx_fsl,
     &                               zmin_fsl+zbeam,dz_fsl)

c       --- Do the parallel nastiness of sharing boundary data
c       --- This is a serial place holder.
        call applyboundaryconditionsonfgrid2d2v(nvx_fsl,nvz_fsl,nx_fsl,nz_fsl,
     &                                          fgrid_fsl,pbounds)

        call enforcepositivity2d2v(nvx_fsl,nvz_fsl,nx_fsl,nz_fsl,fgrid_fsl,
     &                             fthreshold_fsl)

c       --- Regenerate the particles from the distribution
        call createparticlesfromfgrid(pgroup,js,geometry,zbeam,0.,pbounds)
        i1 = pgroup%ins(js+1)
        i2 = i1 + pgroup%nps(js+1) - 1

      enddo

!$OMP MASTER
      if (ltoptimesubs) timeremapparticles2d2v = timeremapparticles2d2v +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine setgrid2dmodbspline2w(np,x,z,w,nx,nz,grid,xmin,xmax,zmin,zmax)
      use Subtimerstop
      integer(ISZ):: np,nx,nz
      real(kind=8):: x(np), z(np), w(np)
      real(kind=8):: grid(-1:nx+1,-1:nz+1)
      real(kind=8):: xmin,xmax,zmin,zmax

c Lays weighted particles down onto a node centered 2-D mesh using a
c modified second order B-spline.
c The shape function is as follows:
c S(x) = | 1 - 5/2 x**2 + 3/2 |x|**3  0 <= |x| <= 1
c        | 1/2 (2 - |x|)**2 (1 - |x|) 1 <= |x| <= 2
c        | 0                          otherwise
c This has an FWHM of 1 grid cell. Note that it is negative in the
c region 1 <= |x| <= 2.
c
c Arugments:
c   np: number of particles
c   x,z: arrays of particle positions
c   w: array of particle weights
c   nx,nz: number of grid cells
c   grid: array where particle number will be accumulated. Note that with
c         normal operation, grid needs to be set to zero before being passed in.
c   xmin,xmax,zmin,zmax: extent of the grid, i.e. grid(0,0) is at (xmin,zmin),
c                        grid(nx,nz) is at (xmax,zmax)

      integer(ISZ):: ip,ix,iz
      real(kind=8):: gx,gz,dxi,dzi,wx,wz
      real(kind=8):: sxm1,sx00,sxp1,sxp2,szm1,sz00,szp1,szp2
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (xmax == xmin .or. zmax == zmin) then
        call kaboom("setgrid2dmodbspline2w: the grid extent is zero, xmin==xmax or zmin==zmax")
        return
      endif

c     --- set grid cell inverse sizes
      dxi = nx/(xmax - xmin)
      dzi = nz/(zmax - zmin)

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gz = (z(ip) - zmin)*dzi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx < nx .and.
     &      0. <= gz .and. gz < nz) then
          ix = int(gx)
          iz = int(gz)
          wx = gx - ix
          wz = gz - iz

          sxm1 = -0.5*wx*(1. - wx)**2
          sx00 = 1. - 2.5*wx**2 + 1.5*wx**3
          sxp1 = 1. - 2.5*(1. - wx)**2 + 1.5*(1. - wx)**3
          sxp2 = 0.5*wx**2*(wx - 1.)
          szm1 = -0.5*wz*(1. - wz)**2
          sz00 = 1. - 2.5*wz**2 + 1.5*wz**3
          szp1 = 1. - 2.5*(1. - wz)**2 + 1.5*(1. - wz)**3
          szp2 = 0.5*wz**2*(wz - 1.)

          grid(ix-1,iz-1) = grid(ix-1,iz-1) + sxm1*szm1*w(ip)
          grid(ix  ,iz-1) = grid(ix  ,iz-1) + sx00*szm1*w(ip)
          grid(ix+1,iz-1) = grid(ix+1,iz-1) + sxp1*szm1*w(ip)
          grid(ix+2,iz-1) = grid(ix+2,iz-1) + sxp2*szm1*w(ip)
          grid(ix-1,iz  ) = grid(ix-1,iz  ) + sxm1*sz00*w(ip)
          grid(ix  ,iz  ) = grid(ix  ,iz  ) + sx00*sz00*w(ip)
          grid(ix+1,iz  ) = grid(ix+1,iz  ) + sxp1*sz00*w(ip)
          grid(ix+2,iz  ) = grid(ix+2,iz  ) + sxp2*sz00*w(ip)
          grid(ix-1,iz+1) = grid(ix-1,iz+1) + sxm1*szp1*w(ip)
          grid(ix  ,iz+1) = grid(ix  ,iz+1) + sx00*szp1*w(ip)
          grid(ix+1,iz+1) = grid(ix+1,iz+1) + sxp1*szp1*w(ip)
          grid(ix+2,iz+1) = grid(ix+2,iz+1) + sxp2*szp1*w(ip)
          grid(ix-1,iz+2) = grid(ix-1,iz+2) + sxm1*szp2*w(ip)
          grid(ix  ,iz+2) = grid(ix  ,iz+2) + sx00*szp2*w(ip)
          grid(ix+1,iz+2) = grid(ix+1,iz+2) + sxp1*szp2*w(ip)
          grid(ix+2,iz+2) = grid(ix+2,iz+2) + sxp2*szp2*w(ip)

        endif

      enddo

!$OMP MASTER
      if (ltoptimesubs) timesetgrid2dmodbspline2w =
     &                  timesetgrid2dmodbspline2w +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine setgrid2d2vmodbspline2w(np,x,z,vx,vz,w,nvx,nvz,nx,nz,grid,
     &                                   vxmin,dvx,vzmin,dvz,
     &                                   xmin,dx,zmin,dz)
      use Subtimerstop
      integer(ISZ):: np,nvx,nvz,nx,nz
      real(kind=8):: x(np), z(np), vx(np), vz(np), w(np)
      real(kind=8):: grid(-1:nvx+1,-1:nvz+1,-1:nx+1,-1:nz+1)
      real(kind=8):: vxmin(-1:nx+1,-1:nz+1),dvx
      real(kind=8):: vzmin(-1:nx+1,-1:nz+1),dvz
      real(kind=8):: xmin,dx,zmin,dz

c Lays weighted particles down onto a node centered 2-D-2-V mesh using a
c modified second order B-spline.
c The shape function is as follows:
c S(x) = | 1 - 5/2 x**2 + 3/2 |x|**3  0 <= |x| <= 1
c        | 1/2 (2 - |x|)**2 (1 - |x|) 1 <= |x| <= 2
c        | 0                          otherwise
c This has an FWHM of 1 grid cell. Note that it is negative in the
c region 1 <= |x| <= 2.
c
c Arugments:
c   np: number of particles
c   x,z,vx,vz: arrays of particle positions and velocities
c   w: array of particle weights
c   nvx,nvz,nx,nz: number of grid cells
c   grid: array where particle number will be accumulated. Note that with
c         normal operation, grid needs to be set to zero before being passed in.
c   vxmin,dvx,vzmin,dvz,xmin,dx,zmin,dz: extent of the grid,
c         i.e. grid(0,0,0,0) is at (vxmin,vzmin,xmin,zmin),
c         grid(nvx,nvz,nx,nz) is at (vxmax,vzmax,xmax,zmax)

      integer(ISZ):: ip,ix,iz,ivx,ivz
      integer(ISZ):: jx,jz,jvx,jvz
      integer(ISZ):: jx2,jz2,jvx2,jvz2
      real(kind=8):: gx,gz,gvx,gvz,dxi,dzi,dvxi,dvzi,wx,wz,wvx,wvz
      real(kind=8):: sx(-1:2),sz(-1:2)
      real(kind=8):: svx(-1:2),svz(-1:2)
      real(kind=8):: sxz
      real(kind=8):: sm1,s00,sp1,sp2
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- set grid cell inverse sizes
      dxi = 1./dx
      dzi = 1./dz
      dvxi = 1./dvx
      dvzi = 1./dvz

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin)*dxi
        gz = (z(ip) - zmin)*dzi
        ix = int(gx)
        iz = int(gz)

        gvx = (vx(ip) - vxmin(ix,iz))*dvxi
        gvz = (vz(ip) - vzmin(ix,iz))*dvzi

c       --- if within grid, accumulate
        if (0. <= gx .and. gx <= nx .and.
     &      0. <= gz .and. gz <= nz .and.
     &      0. <= gvx .and. gvx <= nvx .and.
     &      0. <= gvz .and. gvz <= nvz) then
          ivx = int(gvx)
          ivz = int(gvz)
          wx = gx - ix
          wz = gz - iz
          wvx = gvx - ivx
          wvz = gvz - ivz

c         --- Apply the weight here to reduce the number of multiplies
          sx(-1) = (-0.5*wx*(1. - wx)**2)*w(ip)
          sx( 0) = (1. - 2.5*wx**2 + 1.5*wx**3)*w(ip)
          sx(+1) = (1. - 2.5*(1. - wx)**2 + 1.5*(1. - wx)**3)*w(ip)
          sx(+2) = (0.5*wx**2*(wx - 1.))*w(ip)

          sz(-1) = -0.5*wz*(1. - wz)**2
          sz( 0) = 1. - 2.5*wz**2 + 1.5*wz**3
          sz(+1) = 1. - 2.5*(1. - wz)**2 + 1.5*(1. - wz)**3
          sz(+2) = 0.5*wz**2*(wz - 1.)

          svx(-1) = -0.5*wvx*(1. - wvx)**2
          svx( 0) = 1. - 2.5*wvx**2 + 1.5*wvx**3
          svx(+1) = 1. - 2.5*(1. - wvx)**2 + 1.5*(1. - wvx)**3
          svx(+2) = 0.5*wvx**2*(wvx - 1.)

          svz(-1) = -0.5*wvz*(1. - wvz)**2
          svz( 0) = 1. - 2.5*wvz**2 + 1.5*wvz**3
          svz(+1) = 1. - 2.5*(1. - wvz)**2 + 1.5*(1. - wvz)**3
          svz(+2) = 0.5*wvz**2*(wvz - 1.)

          if (ix == nx) then
            jx2 = 1
          else
            jx2 = 2
          endif
          if (iz == nz) then
            jz2 = 1
          else
            jz2 = 2
          endif
c         if (ivx == nvx) then
c           jvx2 = 1
c         else
c           jvx2 = 2
c         endif
c         if (ivz == nvz) then
c           jvz2 = 1
c         else
c           jvz2 = 2
c         endif

          if (ivx < nvx .and. ivz < nvz) then
            do jz=iz-1,iz+jz2
              do jx=ix-1,ix+jx2
c               do jvz=ivz-1,ivz+jvz2
c                 do jvx=ivx-1,ivx+jvx2

c                   --- Unrolling the loop over jvx made a big difference
c                   --- in the timing, cutting the time in half.
c                   --- Unrolling jvz only made a small change, a few percent.

c                   --- Removing common subexpressions helped by about 15%
                    sxz = sx(jx-ix)*sz(jz-iz)
                    sm1 = svz(-1)*sxz
                    s00 = svz( 0)*sxz
                    sp1 = svz(+1)*sxz
                    sp2 = svz(+2)*sxz

                    grid(ivx-1,ivz-1,jx,jz)=grid(ivx-1,ivz-1,jx,jz)+svx(-1)*sm1
                    grid(ivx  ,ivz-1,jx,jz)=grid(ivx  ,ivz-1,jx,jz)+svx(0)*sm1
                    grid(ivx+1,ivz-1,jx,jz)=grid(ivx+1,ivz-1,jx,jz)+svx(1)*sm1
                    grid(ivx+2,ivz-1,jx,jz)=grid(ivx+2,ivz-1,jx,jz)+svx(2)*sm1

                    grid(ivx-1,ivz  ,jx,jz)=grid(ivx-1,ivz  ,jx,jz)+svx(-1)*s00
                    grid(ivx  ,ivz  ,jx,jz)=grid(ivx  ,ivz  ,jx,jz)+svx(0)*s00
                    grid(ivx+1,ivz  ,jx,jz)=grid(ivx+1,ivz  ,jx,jz)+svx(1)*s00
                    grid(ivx+2,ivz  ,jx,jz)=grid(ivx+2,ivz  ,jx,jz)+svx(2)*s00

                    grid(ivx-1,ivz+1,jx,jz)=grid(ivx-1,ivz+1,jx,jz)+svx(-1)*sp1
                    grid(ivx  ,ivz+1,jx,jz)=grid(ivx  ,ivz+1,jx,jz)+svx(0)*sp1
                    grid(ivx+1,ivz+1,jx,jz)=grid(ivx+1,ivz+1,jx,jz)+svx(1)*sp1
                    grid(ivx+2,ivz+1,jx,jz)=grid(ivx+2,ivz+1,jx,jz)+svx(2)*sp1

                    grid(ivx-1,ivz+2,jx,jz)=grid(ivx-1,ivz+2,jx,jz)+svx(-1)*sp2
                    grid(ivx  ,ivz+2,jx,jz)=grid(ivx  ,ivz+2,jx,jz)+svx(0)*sp2
                    grid(ivx+1,ivz+2,jx,jz)=grid(ivx+1,ivz+2,jx,jz)+svx(1)*sp2
                    grid(ivx+2,ivz+2,jx,jz)=grid(ivx+2,ivz+2,jx,jz)+svx(2)*sp2

c                 enddo
c               enddo
              enddo
            enddo

          else

            do jz=iz-1,iz+jz2
              do jx=ix-1,ix+jx2
c               do jvz=ivz-1,ivz+jvz2
c                 do jvx=ivx-1,ivx+jvx2

c                   --- Unrolling the loop over jvx made a big difference
c                   --- in the timing, cutting the time in half.
c                   --- Unrolling jvz only made a small change, a few percent.

c                   --- Removing common subexpressions helped by about 15%
                    sxz = sx(jx-ix)*sz(jz-iz)
                    sm1 = svz(-1)*sxz
                    s00 = svz( 0)*sxz
                    sp1 = svz(+1)*sxz

                    grid(ivx-1,ivz-1,jx,jz)=grid(ivx-1,ivz-1,jx,jz)+svx(-1)*sm1
                    grid(ivx  ,ivz-1,jx,jz)=grid(ivx  ,ivz-1,jx,jz)+svx(0)*sm1
                    grid(ivx+1,ivz-1,jx,jz)=grid(ivx+1,ivz-1,jx,jz)+svx(1)*sm1

                    grid(ivx-1,ivz  ,jx,jz)=grid(ivx-1,ivz  ,jx,jz)+svx(-1)*s00
                    grid(ivx  ,ivz  ,jx,jz)=grid(ivx  ,ivz  ,jx,jz)+svx(0)*s00
                    grid(ivx+1,ivz  ,jx,jz)=grid(ivx+1,ivz  ,jx,jz)+svx(1)*s00

                    grid(ivx-1,ivz+1,jx,jz)=grid(ivx-1,ivz+1,jx,jz)+svx(-1)*sp1
                    grid(ivx  ,ivz+1,jx,jz)=grid(ivx  ,ivz+1,jx,jz)+svx(0)*sp1
                    grid(ivx+1,ivz+1,jx,jz)=grid(ivx+1,ivz+1,jx,jz)+svx(1)*sp1

                    if (ivx < nvx) then
                      grid(ivx+2,ivz-1,jx,jz)=grid(ivx+2,ivz-1,jx,jz)+svx(2)*sm1
                      grid(ivx+2,ivz  ,jx,jz)=grid(ivx+2,ivz  ,jx,jz)+svx(2)*s00
                      grid(ivx+2,ivz+1,jx,jz)=grid(ivx+2,ivz+1,jx,jz)+svx(2)*sp1
                    endif

                    if (ivz < nvz) then
                     sp2 = svz(+2)*sxz
                     grid(ivx-1,ivz+2,jx,jz)=grid(ivx-1,ivz+2,jx,jz)+svx(-1)*sp2
                     grid(ivx  ,ivz+2,jx,jz)=grid(ivx  ,ivz+2,jx,jz)+svx(0)*sp2
                     grid(ivx+1,ivz+2,jx,jz)=grid(ivx+1,ivz+2,jx,jz)+svx(1)*sp2
c                    if (ivx < nvx) then
c                     grid(ivx+2,ivz+2,jx,jz)=grid(ivx+2,ivz+2,jx,jz)+svx(2)*sp2
c                    endif
                    endif

c                 enddo
c               enddo
              enddo
            enddo

          endif

        endif

      enddo

!$OMP MASTER
      if (ltoptimesubs) timesetgrid2d2vmodbspline2w =
     &                  timesetgrid2d2vmodbspline2w +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
