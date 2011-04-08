#include "top.h"
c=============================================================================
c Routines for handing the Forward semi-Lagrangian method
c=============================================================================
      subroutine applyboundaryconditionsonfgrid2d2v(nvx,nvz,nx,nz,
     &                                              nxguard,nzguard,
     &                                              fgrid,
     &                                              localpbounds,
     &                                              flost,lcellcentered)
      use Subtimerstop
      integer(ISZ):: nvx,nvz,nx,nz,nxguard,nzguard
      real(kind=8):: fgrid(-1:nvx+1,-1:nvz+1,-nxguard:nx+1,-nzguard:nz+1)
      integer(ISZ):: localpbounds(0:5)
      real(kind=8):: flost
      logical(ISZ):: lcellcentered

      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- The the data from the guard cells on velocity and move it to
c     --- the edge of the domain. This avoids throwing data away, but
c     --- does make errors in energy and momentum.
c     --- Note that the velocity grid is always node centered.
      fgrid(0,:,:,:) = fgrid(0,:,:,:) + fgrid(-1,:,:,:)
      fgrid(-1,:,:,:) = 0.
      fgrid(nvx,:,:,:) = fgrid(nvx,:,:,:) + fgrid(nvx+1,:,:,:)
      fgrid(nvx+1,:,:,:) = 0.
      fgrid(:,0,:,:) = fgrid(:,0,:,:) + fgrid(:,-1,:,:)
      fgrid(:,-1,:,:) = 0.
      fgrid(:,nvz,:,:) = fgrid(:,nvz,:,:) + fgrid(:,nvz+1,:,:)
      fgrid(:,nvz+1,:,:) = 0.

      if (lcellcentered) then

c       --- The grid includes ix=0 through ix=nx-1. ix=-2,-1,nx,nx+1 are all
c       --- guard cells.

        if (localpbounds(0) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,-2:-1,:))
          fgrid(:,:,-2:-1,:) = 0.
        endif
        if (localpbounds(1) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,nx:nx+1,:))
          fgrid(:,:,nx:nx+1,:) = 0.
        endif
        if (localpbounds(0) == 1) then
          fgrid(:,:,0,:) = fgrid(:,:,0,:) + fgrid(:,:,-1,:)
          fgrid(:,:,1,:) = fgrid(:,:,1,:) + fgrid(:,:,-2,:)
          fgrid(:,:,-1,:) = fgrid(:,:,0,:)
          fgrid(:,:,-2,:) = fgrid(:,:,1,:)
        endif
        if (localpbounds(1) == 1) then
          fgrid(:,:,nx-1,:) = fgrid(:,:,nx-1,:) + fgrid(:,:,nx,:)
          fgrid(:,:,nx-2,:) = fgrid(:,:,nx-2,:) + fgrid(:,:,nx+1,:)
          fgrid(:,:,nx,:) = fgrid(:,:,nx-1,:)
          fgrid(:,:,nx+1,:) = fgrid(:,:,nx-2,:)
        endif
#ifndef MPIPARALLEL
        if (localpbounds(0) == 2 .and. localpbounds(1) == 2) then
          fgrid(:,:,-2:1,:) = fgrid(:,:,-2:1,:) + fgrid(:,:,nx-2:nx+1,:)
          fgrid(:,:,nx-2:nx+1,:) = fgrid(:,:,-2:1,:)
        endif
#endif

        if (localpbounds(4) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,:,-2:-1))
          fgrid(:,:,:,-2:-1) = 0.
        endif
        if (localpbounds(5) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,:,nz:nz+1))
          fgrid(:,:,:,nz:nz+1) = 0.
        endif
        if (localpbounds(4) == 1) then
          fgrid(:,:,:,0) = fgrid(:,:,:,0) + fgrid(:,:,:,-1)
          fgrid(:,:,:,1) = fgrid(:,:,:,1) + fgrid(:,:,:,-2)
          fgrid(:,:,:,-1) = fgrid(:,:,:,0)
          fgrid(:,:,:,-2) = fgrid(:,:,:,1)
        endif
        if (localpbounds(5) == 1) then
          fgrid(:,:,:,nz-1) = fgrid(:,:,:,nz-1) + fgrid(:,:,:,nz)
          fgrid(:,:,:,nz-2) = fgrid(:,:,:,nz-2) + fgrid(:,:,:,nz+1)
          fgrid(:,:,:,nz) = fgrid(:,:,:,nz-1)
          fgrid(:,:,:,nz+1) = fgrid(:,:,:,nz-2)
        endif
#ifndef MPIPARALLEL
        if (localpbounds(4) == 2 .and. localpbounds(5) == 2) then
          fgrid(:,:,:,-2:1) = fgrid(:,:,:,-2:1) + fgrid(:,:,:,nz-2:nz+1)
          fgrid(:,:,:,nz-2:nz+1) = fgrid(:,:,:,-2:1)
        endif
#endif

      else if (.not. lcellcentered) then

        if (localpbounds(0) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,-1:0,:))
          fgrid(:,:,-1:0,:) = 0.
        endif
        if (localpbounds(1) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,nx:nx+1,:))
          fgrid(:,:,nx:nx+1,:) = 0.
        endif
        if (localpbounds(0) == 1) then
          fgrid(:,:,1,:) = fgrid(:,:,1,:) + fgrid(:,:,-1,:)
          fgrid(:,:,-1,:) = fgrid(:,:,1,:)
          fgrid(:,:,0,:) = 2.*fgrid(:,:,0,:)
        endif
        if (localpbounds(1) == 1) then
          fgrid(:,:,nx-1,:) = fgrid(:,:,nx-1,:) + fgrid(:,:,nx+1,:)
          fgrid(:,:,nx+1,:) = fgrid(:,:,nx-1,:)
          fgrid(:,:,nx,:) = 2.*fgrid(:,:,nx,:)
        endif
#ifndef MPIPARALLEL
        if (localpbounds(0) == 2 .and. localpbounds(1) == 2) then
          fgrid(:,:,-1:1,:) = fgrid(:,:,-1:1,:) + fgrid(:,:,nx-1:nx+1,:)
          fgrid(:,:,nx-1:nx+1,:) = fgrid(:,:,-1:1,:)
        endif
#endif

        if (localpbounds(4) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,:,-1:0))
          fgrid(:,:,:,-1:0) = 0.
        endif
        if (localpbounds(5) == 0) then
c         --- Save the sum of the data in the guard cell which is about to
c         --- be thrown away.
          flost = flost + sum(fgrid(:,:,:,nz:nz+1))
          fgrid(:,:,:,nz:nz+1) = 0.
        endif
        if (localpbounds(4) == 1) then
          fgrid(:,:,:,1) = fgrid(:,:,:,1) + fgrid(:,:,:,-1)
          fgrid(:,:,:,-1) = fgrid(:,:,:,1)
          fgrid(:,:,:,0) = 2.*fgrid(:,:,:,0)
        endif
        if (localpbounds(5) == 1) then
          fgrid(:,:,:,nz-1) = fgrid(:,:,:,nz-1) + fgrid(:,:,:,nz+1)
          fgrid(:,:,:,nz+1) = fgrid(:,:,:,nz-1)
          fgrid(:,:,:,nz) = 2.*fgrid(:,:,:,nz)
        endif
#ifndef MPIPARALLEL
        if (localpbounds(4) == 2 .and. localpbounds(5) == 2) then
          fgrid(:,:,:,-1:1) = fgrid(:,:,:,-1:1) + fgrid(:,:,:,nz-1:nz+1)
          fgrid(:,:,:,nz-1:nz+1) = fgrid(:,:,:,-1:1)
        endif
#endif

      endif

#ifdef MPIPARALLEL
      call sumfondomainboundaries(nvx,nvz,nx,nz,nxguard,nzguard,
     &                            fgrid,localpbounds)
#endif

!$OMP MASTER
      if (ltoptimesubs) timeapplyboundaryconditionsonfgrid2d2v =
     &                  timeapplyboundaryconditionsonfgrid2d2v +
     &                                         wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine applyminmaxboundaryconditions2d(nx,nz,nxguard,nzguard,
     &                                           v,pbounds,lmax,lcellcentered)
      use Subtimerstop
      integer(ISZ):: nx,nz,nxguard,nzguard
      real(kind=8):: v(-nxguard:nx+1,-nzguard:nz+1)
      integer(ISZ):: pbounds(0:5)
      logical(ISZ):: lmax
      logical(ISZ):: lcellcentered

      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (lcellcentered) then

        if (pbounds(0) == 0 .or. pbounds(0) == 1) then
          v(-1,:) = v(0,:)
          v(-2,:) = v(0,:)
        endif
        if (pbounds(1) == 0 .or. pbounds(1) == 1) then
          v(nx,:) = v(nx-1,:)
          v(nx+1,:) = v(nx-1,:)
        endif
        if (pbounds(0) == 2 .and. pbounds(1) == 2) then
          if (lmax) then
            v(-2:1,:) = max(v(-2:1,:),v(nx-2:nx+1,:))
          else
            v(-2:1,:) = min(v(-2:1,:),v(nx-2:nx+1,:))
          endif
          v(nx:nx+1,:) = v(-2:-1,:)
        endif

        if (pbounds(4) == 0 .or. pbounds(4) == 1) then
          v(:,-1) = v(:,0)
          v(:,-2) = v(:,0)
        endif
        if (pbounds(5) == 0 .or. pbounds(5) == 1) then
          v(:,nz) = v(:,nz-1)
          v(:,nz+1) = v(:,nz-1)
        endif
        if (pbounds(4) == 2 .and. pbounds(5) == 2) then
          if (lmax) then
            v(:,-2:1) = max(v(:,-2:1),v(:,nz-2:nz+1))
          else
            v(:,-2:1) = min(v(:,-2:1),v(:,nz-2:nz+1))
          endif
          v(:,nz-2:nz+1) = v(:,-2:1)
        endif

      else

        if (pbounds(0) == 0 .or. pbounds(0) == 1) v(-1,:) = v(1,:)
        if (pbounds(1) == 0 .or. pbounds(1) == 1) v(nx+1,:) = v(nx-1,:)
        if (pbounds(0) == 2 .and. pbounds(1) == 2) then
          if (lmax) then
            v(-1:1,:) = max(v(-1:1,:),v(nx-1:nx+1,:))
          else
            v(-1:1,:) = min(v(-1:1,:),v(nx-1:nx+1,:))
          endif
          v(nx-1:nx+1,:) = v(-1:1,:)
        endif

        if (pbounds(4) == 0 .or. pbounds(4) == 1) v(:,-1) = v(:,1)
        if (pbounds(5) == 0 .or. pbounds(5) == 1) v(:,nz+1) = v(:,nz-1)
        if (pbounds(4) == 2 .and. pbounds(5) == 2) then
          if (lmax) then
            v(:,-1:1) = max(v(:,-1:1),v(:,nz-1:nz+1))
          else
            v(:,-1:1) = min(v(:,-1:1),v(:,nz-1:nz+1))
          endif
          v(:,nz-1:nz+1) = v(:,-1:1)
        endif

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
     &                                    localpbounds)
      use Subtimerstop
      use ForwardSemiLagrange
      use ParticleGroupmodule
      use Particles,Only: wpid
      type(ParticleGroup):: pgroup
      integer(ISZ):: js,geometry
      real(kind=8):: zbeam,vbeam
      integer(ISZ):: localpbounds(0:5)

      real(kind=8),pointer:: xp(:),yp(:),zp(:),uxp(:),uyp(:),uzp(:)
      real(kind=8):: vybeam,vzbeam
      integer(ISZ):: npnew,ip
      integer(ISZ):: ix,iz,ivx,ivz
      integer(ISZ):: ix0,iz0
      integer(ISZ):: ix1,iz1
      real(kind=8):: xoffset,zoffset
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (wpid == 0) then
        call kaboom("createparticlesfromfgrid: wpid must be set")
        return
      endif

      if (lcellcentered_fsl) then

c       --- When cell centered, these values are independent of the
c       --- boundary conditions.
        ix0 = 0
        ix1 = nx_fsl - 1
        iz0 = 0
        iz1 = nz_fsl - 1
        xoffset = dx_fsl/2.
        zoffset = dz_fsl/2.

      else

c       --- Be careful not to create particles at the edge of the mesh
c       --- with Dirichelt boundaries or the upper edge with peridoic.
        if (localpbounds(0) == 0) then
          ix0 = 1
        else
          ix0 = 0
        endif
        if (localpbounds(1) == 0 .or. localpbounds(1) == 2 .or.
     &      localpbounds(1) == -1) then
          ix1 = nx_fsl - 1
        else
          ix1 = nx_fsl
        endif
  
        if (localpbounds(4) == 0) then
          iz0 = 1
        else
          iz0 = 0
        endif
        if (localpbounds(5) == 0 .or. localpbounds(5) == 2 .or.
     &      localpbounds(5) == -1) then
          iz1 = nz_fsl - 1
        else
          iz1 = nz_fsl
        endif
        xoffset = 0.
        zoffset = 0.

      endif

c     --- Count how many particles will be created, the number of fgrid
c     --- cells that are positive.
c     --- Note that particles will be created everywhere that f is
c     --- nonzero (even if negative). This is needed for conservation
c     --- of f.
      npnew = 0
      do iz=iz0,iz1
        do ix=ix0,ix1
          do ivz=0,nvz_fsl
            do ivx=0,nvx_fsl
              if (fgrid_fsl(ivx,ivz,ix,iz) .ne. 0.) npnew = npnew + 1
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
              if (fgrid_fsl(ivx,ivz,ix,iz) .ne. 0.) then
                xp(ip) = xmin_fsl + ix*dx_fsl + xoffset
                yp(ip) = 0.
                zp(ip) = zmin_fsl + iz*dz_fsl + zoffset + zbeam
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

      integer(ISZ):: ix,iz,iparity,ix0,niterinner,dd,id
      real(kind=8):: sumcap,ferror,fdel
      integer(ISZ):: niterinnermax
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Repeat the operation until there are no negative locations.
      niterinnermax = min(nx/2,nz/2)

c     --- Loop over the array using even-odd ordering
      do iparity = 0,1
        do iz = 0,nz
          ix0 = mod(iz + iparity,2)
          do ix = ix0,nx,2

c           --- For each location, iterate until its value is nonzero
c           --- (of maximum number of iterations is reached).
c           --- For each iteration, if sumcap is zero, then increase
c           --- dd so that a larger region will be searched.
            niterinner = 0
            dd = 1
            do while (fgrid(ix,iz) < threshold .and.
     &                fgrid(ix,iz) .ne. 0. .and.
     &                niterinner < niterinnermax)
              niterinner = niterinner + 1

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
      subroutine enforcepositivity2d2v(nvx,nvz,nx,nz,nxguard,nzguard,
     &                                 fgrid,threshold,lcellcentered)
      use Subtimerstop
      integer(ISZ):: nvx,nvz,nx,nz,nxguard,nzguard
      real(kind=8):: fgrid(-1:nvx+1,-1:nvz+1,-nxguard:nx+1,-nzguard:nz+1)
      real(kind=8):: threshold
      logical(ISZ):: lcellcentered

c Given a distribution function, fgrid, this routine enforces positivity
c by taking any negative values and distributing it to nearby, already
c positive, locations. It ensures that the integral of fgrid remains fixed.
c The action is local, avoiding a global calculation.

c This sweeps from the edge of the grid inward. This ordering seems to
c work well and doesn't leave isolated cells. The distribution is over
c neighbors, not including diagonals (though this may not be important).
c A future version might allow the mid point to move, perhaps being the
c location of the max f for example.

      integer(ISZ):: kx1,kz1,k1
      integer(ISZ):: kx2,kz2
      integer(ISZ):: kvx,kvz,kx,kz
      integer(ISZ):: kvxmid,kvzmid,kxmid,kzmid
      integer(ISZ):: ivx,ivz,ix,iz
      integer(ISZ):: niter,nitermax
      integer(ISZ):: jvx1,jvz1,jx1,jz1
      integer(ISZ):: jvx2,jvz2,jx2,jz2
      real(kind=8):: sumcap,fdel,fminimum
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      nitermax = 5
      kvxmid = nvx/2
      kvzmid = nvz/2
      kxmid = nx/2
      kzmid = nz/2

      if (lcellcentered) then
        kx1 = 0
        kx2 = nx - 1
        kz1 = 0
        kz2 = nz - 1
        k1 = 1
      else
        kx1 = 0
        kx2 = nx
        kz1 = 0
        kz2 = nz
        k1 = 0
      endif

      niter = 0
      do
        niter = niter + 1
        if (niter > nitermax) exit

        fminimum = threshold

c       --- Loop over the array, sweeping in from the edges to the center
        do kz = kz1,kz2
          if (kz < kzmid) then
            iz = kz
          else
            iz = nz - kz + kzmid - k1
          endif

          do kx = kx1,kx2
            if (kx < kxmid) then
              ix = kx
            else
              ix = nx - kx + kxmid - k1
            endif

            do kvz = 0,nvz
              if (kvz < kvzmid) then
                ivz = kvz
              else
                ivz = nvz - kvz + kvzmid
              endif

              do kvx = 0,nvx
                if (kvx < kvxmid) then
                  ivx = kvx
                else
                  ivx = nvx - kvx + kvxmid
                endif

                if (fgrid(ivx,ivz,ix,iz) < threshold .and.
     &              fgrid(ivx,ivz,ix,iz) .ne. 0.) then

c                 --- Create the sum of the surrounding, positive values,
c                 --- including guard cells, but excluding corners.
                  sumcap = 0.
                  if (fgrid(ivx-1,ivz  ,ix  ,iz  ) >= threshold)
     &              sumcap = sumcap + fgrid(ivx-1,ivz  ,ix  ,iz  )
                  if (fgrid(ivx+1,ivz  ,ix  ,iz  ) >= threshold)
     &              sumcap = sumcap + fgrid(ivx+1,ivz  ,ix  ,iz  )
                  if (fgrid(ivx  ,ivz-1,ix  ,iz  ) >= threshold)
     &              sumcap = sumcap + fgrid(ivx  ,ivz-1,ix  ,iz  )
                  if (fgrid(ivx  ,ivz+1,ix  ,iz  ) >= threshold)
     &              sumcap = sumcap + fgrid(ivx  ,ivz+1,ix  ,iz  )
                  if (fgrid(ivx  ,ivz  ,ix-1,iz  ) >= threshold)
     &              sumcap = sumcap + fgrid(ivx  ,ivz  ,ix-1,iz  )
                  if (fgrid(ivx  ,ivz  ,ix+1,iz  ) >= threshold)
     &              sumcap = sumcap + fgrid(ivx  ,ivz  ,ix+1,iz  )
                  if (fgrid(ivx  ,ivz  ,ix  ,iz-1) >= threshold)
     &              sumcap = sumcap + fgrid(ivx  ,ivz  ,ix  ,iz-1)
                  if (fgrid(ivx  ,ivz  ,ix  ,iz+1) >= threshold)
     &              sumcap = sumcap + fgrid(ivx  ,ivz  ,ix  ,iz+1)

                  if (sumcap > 0.) then
c                   --- If the sum is > 0, then the distribution of the
c                   --- negative f can be done. If there is not enough,
c                   --- some areas will go negative, and will need to
c                   --- be corrected in the next iteration.
                    fdel = fgrid(ivx,ivz,ix,iz)/sumcap
                    fgrid(ivx,ivz,ix,iz) = 0.
                    if (fgrid(ivx-1,ivz  ,ix  ,iz  ) >= threshold) then
                      fgrid(ivx-1,ivz  ,ix  ,iz  ) = 
     &                fgrid(ivx-1,ivz  ,ix  ,iz  ) + 
     &                fgrid(ivx-1,ivz  ,ix  ,iz  )*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx-1,ivz  ,ix  ,iz  ))
                    endif
                    if (fgrid(ivx+1,ivz  ,ix  ,iz  ) >= threshold) then
                      fgrid(ivx+1,ivz  ,ix  ,iz  ) = 
     &                fgrid(ivx+1,ivz  ,ix  ,iz  ) + 
     &                fgrid(ivx+1,ivz  ,ix  ,iz  )*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx+1,ivz  ,ix  ,iz  ))
                    endif
                    if (fgrid(ivx  ,ivz-1,ix  ,iz  ) >= threshold) then
                      fgrid(ivx  ,ivz-1,ix  ,iz  ) = 
     &                fgrid(ivx  ,ivz-1,ix  ,iz  ) + 
     &                fgrid(ivx  ,ivz-1,ix  ,iz  )*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx  ,ivz-1,ix  ,iz  ))
                    endif
                    if (fgrid(ivx  ,ivz+1,ix  ,iz  ) >= threshold) then
                      fgrid(ivx  ,ivz+1,ix  ,iz  ) = 
     &                fgrid(ivx  ,ivz+1,ix  ,iz  ) + 
     &                fgrid(ivx  ,ivz+1,ix  ,iz  )*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx  ,ivz+1,ix  ,iz  ))
                    endif
                    if (fgrid(ivx  ,ivz  ,ix-1,iz  ) >= threshold) then
                      fgrid(ivx  ,ivz  ,ix-1,iz  ) = 
     &                fgrid(ivx  ,ivz  ,ix-1,iz  ) + 
     &                fgrid(ivx  ,ivz  ,ix-1,iz  )*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx  ,ivz  ,ix-1,iz  ))
                    endif
                    if (fgrid(ivx  ,ivz  ,ix+1,iz  ) >= threshold) then
                      fgrid(ivx  ,ivz  ,ix+1,iz  ) = 
     &                fgrid(ivx  ,ivz  ,ix+1,iz  ) + 
     &                fgrid(ivx  ,ivz  ,ix+1,iz  )*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx  ,ivz  ,ix+1,iz  ))
                    endif
                    if (fgrid(ivx  ,ivz  ,ix  ,iz-1) >= threshold) then
                      fgrid(ivx  ,ivz  ,ix  ,iz-1) = 
     &                fgrid(ivx  ,ivz  ,ix  ,iz-1) + 
     &                fgrid(ivx  ,ivz  ,ix  ,iz-1)*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx  ,ivz  ,ix  ,iz-1))
                    endif
                    if (fgrid(ivx  ,ivz  ,ix  ,iz+1) >= threshold) then
                      fgrid(ivx  ,ivz  ,ix  ,iz+1) = 
     &                fgrid(ivx  ,ivz  ,ix  ,iz+1) + 
     &                fgrid(ivx  ,ivz  ,ix  ,iz+1)*fdel
                      fminimum = min(fminimum,
     &                               fgrid(ivx  ,ivz  ,ix  ,iz+1))
                    endif

                  endif

                endif

              enddo
            enddo
          enddo
        enddo

        if (fminimum >= threshold .or. fminimum == 0.) exit

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
     &                                nxguard,nzguard,
     &                                xmin,dx,zmin,dz,
     &                                vxmin,vxmax,vzmin,vzmax,
     &                                lcellcentered)
      use Subtimerstop
      integer(ISZ):: np,nx,nz,nxguard,nzguard
      real(kind=8):: x(np), z(np), vx(np), vz(np)
      real(kind=8):: xmin,dx,zmin,dz
      real(kind=8):: vxmin(-nxguard:nx+1,-nzguard:nz+1)
      real(kind=8):: vxmax(-nxguard:nx+1,-nzguard:nz+1)
      real(kind=8):: vzmin(-nxguard:nx+1,-nzguard:nz+1)
      real(kind=8):: vzmax(-nxguard:nx+1,-nzguard:nz+1)
      logical(ISZ):: lcellcentered

c Finds the min and max velocities as a function of position

      real(kind=8):: xoffset,zoffset
      real(kind=8):: gxoffset,gzoffset
      integer(ISZ):: ip,ix,iz,i1
      real(kind=8):: gx,gz,dxi,dzi
      real(kind=8):: lastvalue1,lastvalue2
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- set grid cell inverse sizes
      dxi = 1./dx
      dzi = 1./dz

      if (lcellcentered) then
        xoffset = dx/2.
        zoffset = dz/2.
        gxoffset = 1.
        gzoffset = 1.
        i1 = 1
      else
        xoffset = 0.
        zoffset = 0.
        gxoffset = 0.
        gzoffset = 0.
        i1 = 0
      endif

      vxmin = +LARGEPOS
      vxmax = -LARGEPOS
      vzmin = +LARGEPOS
      vzmax = -LARGEPOS

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin - xoffset)*dxi
        gz = (z(ip) - zmin - zoffset)*dzi

c       --- if within grid, accumulate
        if (0.-gxoffset <= gx .and. gx < nx .and.
     &      0.-gzoffset <= gz .and. gz < nz) then
          ix = floor(gx)
          iz = floor(gz)
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
        real(kind=8):: v(-nxguard:nx+1,-nzguard:nz+1)
        integer(ISZ):: lpsign
c       --- Remove the LARGEPOS values, replacing them with nearby values.
c       --- Two sweeps are done, one in each direction, and any LARGEPOS
c       --- are replaced by the previous value found during the sweep.
        do iz=-nzguard,nz+1
          lastvalue1 = lpsign*LARGEPOS
          lastvalue2 = lpsign*LARGEPOS
          do ix=-nxguard,nx+1
            if (v(ix,iz) == lpsign*LARGEPOS) then
              v(ix,iz) = lastvalue1
            else
              lastvalue1 = v(ix,iz)
            endif
            if (v(nx-ix-i1,iz) == lpsign*LARGEPOS) then
              v(nx-ix-i1,iz) = lastvalue2
            else
              lastvalue2 = v(nx-ix-i1,iz)
            endif
          enddo
        enddo
        do ix=-nxguard,nx+1
          lastvalue1 = lpsign*LARGEPOS
          lastvalue2 = lpsign*LARGEPOS
          do iz=-nzguard,nz+1
            if (v(ix,iz) == lpsign*LARGEPOS) then
              v(ix,iz) = lastvalue1
            else
              lastvalue1 = v(ix,iz)
            endif
            if (v(ix,nz-iz-i1) == lpsign*LARGEPOS) then
              v(ix,nz-iz-i1) = lastvalue2
            else
              lastvalue2 = v(ix,nz-iz-i1)
            endif
          enddo
        enddo

        end subroutine cleanuplargepos

      end
c=============================================================================
      subroutine loadfgridguassian(nz,density,radius,
     &                             vthermalperp,vthermalparallel,
     &                             lcellcentered)
      use Subtimerstop
      use Constant
      use ForwardSemiLagrange
      integer(ISZ):: nz
      real(kind=8):: density(0:nz),radius(0:nz)
      real(kind=8):: vthermalperp,vthermalparallel
      logical(ISZ):: lcellcentered

      real(kind=8):: xoffset,zoffset
      integer(ISZ):: ivx,ivz,ix,iz,ix2,iz2
      real(kind=8):: xx,vx,vz,fvz,fvx,sumf
      real(kind=8):: fvx0,fvz0
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      lcellcentered_fsl = lcellcentered
      if (lcellcentered) then
        nxguard_fsl = 2
        nzguard_fsl = 2
        xoffset = dx_fsl/2.
        zoffset = dz_fsl/2.
        ix2 = nx_fsl - 1
        iz2 = nz_fsl - 1
      else
        nxguard_fsl = 1
        nzguard_fsl = 1
        xoffset = 0.
        zoffset = 0.
        ix2 = nx_fsl
        iz2 = nz_fsl
      endif

      if (nz .ne. nz_fsl) then
        call kaboom("loadfgridguassian: nz must be equal to nz_fsl")
        return
      endif

      if (dvx_fsl == 0. .and. nvx_fsl > 0) then
        dvx_fsl = 8.*vthermalperp/nvx_fsl
      else if (nvx_fsl == 0 .and. dvx_fsl > 0.) then
        nvx_fsl = 8.*vthermalperp/dvx_fsl
      else
        call kaboom("loadfgridguassian: either dvx_fsl or nvx_fsl must be nonzero")
        return
      endif
      if (dvz_fsl == 0. .and. nvz_fsl > 0) then
        dvz_fsl = 8.*vthermalparallel/nvz_fsl
      else if (nvz_fsl == 0 .and. dvz_fsl > 0.) then
        nvz_fsl = 8.*vthermalparallel/dvz_fsl
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

      do iz=0,iz2
        if (density(iz) == 0. .or. radius(iz) == 0.) cycle
        do ix=0,ix2
          xx = xmin_fsl + ix*dx_fsl + xoffset
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
      use Parallel
      type(ParticleGroup):: pgroup
      integer(ISZ):: geometry
      real(kind=8):: zbeam
      logical(ISZ):: l4symtry,l2symtry

      real(kind=8),pointer:: xp(:),zp(:),uxp(:),uzp(:)
      integer(ISZ):: js,ip,i1,i2
      integer(ISZ):: localpbounds(0:5)
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Setup the boundary conditions for the local domain.
      call setboundsfromflags(pbounds,pboundxy,pbound0,pboundnz,
     &                        l2symtry,l4symtry)
      localpbounds = pbounds
#ifdef MPIPARALLEL
      if (ppdecomp%ixproc > 0)                  localpbounds(0) = -1
      if (ppdecomp%ixproc < ppdecomp%nxprocs-1) localpbounds(1) = -1
      if (ppdecomp%izproc > 0)                  localpbounds(4) = -1
      if (ppdecomp%izproc < ppdecomp%nzprocs-1) localpbounds(5) = -1
#endif

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

c --- Skip the calculation of the min and max for now. This leads to
c --- an ever growing grid if nonzero f gets to the edge f the mesh.
c --- Future versions may pass the threshold into findminmaxongrid2d2v
c --- so the min and max of only those cells with f > threshold are
c --- counted.

c       --- Find the min and max velocity as a function of position
c       call findminmaxongrid2d2v(pgroup%nps(js+1),xp(i1),zp(i1),
c    &                            uxp(i1),uzp(i1),nx_fsl,nz_fsl,
c    &                            nxguard_fsl,nzguard_fsl,
c    &                            xmin_fsl,dx,zmin_fsl+zbeam,dz,
c    &                            vxmin_fsl,vxmax_fsl,vzmin_fsl,vzmax_fsl,
c    &                            lcellcentered_fsl)

c       --- Do the parallel nastiness of sharing boundary data
c       --- This is a serial place holder.
c       call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,
c    &                                       nxguard_fsl,nzguard_fsl,
c    &                                       vxmin_fsl,
c    &                                       localpbounds,.false.,
c    &                                       lcellcentered_fsl)
c       call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,
c    &                                       nxguard_fsl,nzguard_fsl,
c    &                                       vxmax_fsl,
c    &                                       localpbounds,.true.,
c    &                                       lcellcentered_fsl)
c       call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,
c    &                                       nxguard_fsl,nzguard_fsl,
c    &                                       vzmin_fsl,
c    &                                       localpbounds,.false.,
c    &                                       lcellcentered_fsl)
c       call applyminmaxboundaryconditions2d(nx_fsl,nz_fsl,
c    &                                       nxguard_fsl,nzguard_fsl,
c    &                                       vzmax_fsl,
c    &                                       localpbounds,.true.,
c    &                                       lcellcentered_fsl)

c       --- Round the velocity mins and maxs to the grid cells
c       vxmin_fsl = floor(vxmin_fsl/dvx_fsl)*dvx_fsl
c       vzmin_fsl = floor(vzmin_fsl/dvz_fsl)*dvz_fsl
c       vxmax_fsl = ceiling(vxmax_fsl/dvx_fsl)*dvx_fsl
c       vzmax_fsl = ceiling(vzmax_fsl/dvz_fsl)*dvz_fsl
c       --- Temporary until there are better diagnostics
c       vxmin_fsl = minval(vxmin_fsl)
c       vxmax_fsl = maxval(vxmax_fsl)
c       vzmin_fsl = minval(vzmin_fsl)
c       vzmax_fsl = maxval(vzmax_fsl)

c       --- Resize the grid if needed
c       nvx_fsl = max(nvx_fsl,nint(maxval((vxmax_fsl - vxmin_fsl)/dvx_fsl)))
c       nvz_fsl = max(nvz_fsl,nint(maxval((vzmax_fsl - vzmin_fsl)/dvz_fsl)))
c       call gchange("ForwardSemiLagrange",0)

c       --- Deposit the particles into the distribution grid
        fgrid_fsl = 0.
        call setgrid2d2vmodbspline2w(pgroup%nps(js+1),
     &                               xp(i1),zp(i1),uxp(i1),uzp(i1),
     &                               pgroup%pid(i1,wpid),
     &                               nvx_fsl,nvz_fsl,
     &                               nx_fsl,nz_fsl,nxguard_fsl,nzguard_fsl,
     &                               fgrid_fsl,
     &                               vxmin_fsl,dvx_fsl,vzmin_fsl,dvz_fsl,
     &                               xmin_fsl,dx_fsl,
     &                               zmin_fsl+zbeam,dz_fsl,lcellcentered_fsl)

c       --- Do the parallel nastiness of sharing boundary data
c       --- This is a serial place holder.
        call applyboundaryconditionsonfgrid2d2v(nvx_fsl,nvz_fsl,nx_fsl,nz_fsl,
     &                                          nxguard_fsl,nzguard_fsl,
     &                                          fgrid_fsl,localpbounds,
     &                                          flost_fsl,lcellcentered_fsl)

        call enforcepositivity2d2v(nvx_fsl,nvz_fsl,nx_fsl,nz_fsl,
     &                             nxguard_fsl,nzguard_fsl,
     &                             fgrid_fsl,fthreshold_fsl,lcellcentered_fsl)

c       --- Regenerate the particles from the distribution
        call createparticlesfromfgrid(pgroup,js,geometry,zbeam,0.,localpbounds)

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
      subroutine setgrid2d2vmodbspline2w(np,x,z,vx,vz,w,nvx,nvz,nx,nz,
     &                                   nxguard,nzguard,grid,
     &                                   vxmin,dvx,vzmin,dvz,
     &                                   xmin,dx,zmin,dz,lcellcentered)
      use Subtimerstop
      integer(ISZ):: np,nvx,nvz,nx,nz,nxguard,nzguard
      real(kind=8):: x(np), z(np), vx(np), vz(np), w(np)
      real(kind=8):: grid(-1:nvx+1,-1:nvz+1,-nxguard:nx+1,-nzguard:nz+1)
      real(kind=8):: vxmin(-nxguard:nx+1,-nzguard:nz+1),dvx
      real(kind=8):: vzmin(-nxguard:nx+1,-nzguard:nz+1),dvz
      real(kind=8):: xmin,dx,zmin,dz
      logical(ISZ):: lcellcentered

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

      real(kind=8):: xoffset,zoffset
      real(kind=8):: gxoffset,gzoffset
      integer(ISZ):: ip,ix,iz,ivx,ivz,i1
      integer(ISZ):: jx,jz,jvx,jvz
      integer(ISZ):: jx2,jz2,jvx2,jvz2
      real(kind=8):: gx,gz,gvx,gvz,dxi,dzi,dvxi,dvzi,wx,wz,wvx,wvz
      real(kind=8):: sx(-1:2),sz(-1:2)
      real(kind=8):: svx(-1:2),svz(-1:2)
      real(kind=8):: sxz
      real(kind=8):: sm1,s00,sp1,sp2
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      if (lcellcentered) then
        xoffset = dx/2.
        zoffset = dz/2.
        gxoffset = 1.
        gzoffset = 1.
        i1 = 1
      else
        xoffset = 0.
        zoffset = 0.
        gxoffset = 0.
        gzoffset = 0.
        i1 = 0
      endif

c     --- set grid cell inverse sizes
      dxi = 1./dx
      dzi = 1./dz
      dvxi = 1./dvx
      dvzi = 1./dvz

c     --- loop over particles
      do ip=1,np

c       --- find location on grid
        gx = (x(ip) - xmin - xoffset)*dxi
        gz = (z(ip) - zmin - zoffset)*dzi
        ix = floor(gx)
        iz = floor(gz)

        gvx = (vx(ip) - vxmin(ix,iz))*dvxi
        gvz = (vz(ip) - vzmin(ix,iz))*dvzi

c       --- Particles are always accumulated. Outliers are artificially
c       --- put at the edge of the mesh.
        ivx = int(gvx)
        ivz = int(gvz)
        wx = gx - ix
        wz = gz - iz
        wvx = gvx - ivx
        wvz = gvz - ivz

c       --- Apply the weight here to reduce the number of multiplies
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

        if (0.-gxoffset <= gx .and. gx < nx .and.
     &      0.-gzoffset <= gz .and. gz < nz .and.
     &      0. <= gvx .and. gvx < nvx .and.
     &      0. <= gvz .and. gvz < nvz) then
          do jz=iz-1,iz+jz2
            do jx=ix-1,ix+jx2
c             do jvz=ivz-1,ivz+jvz2
c               do jvx=ivx-1,ivx+jvx2

c                 --- Unrolling the loop over jvx made a big difference
c                 --- in the timing, cutting the time in half.
c                 --- Unrolling jvz only made a small change, a few percent.

c                 --- Removing common subexpressions helped by about 15%
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

c               enddo
c             enddo
            enddo
          enddo

        else

          if (gvx < 0.) then
            ivx = 0
            svx = (/0.,1.,0.,0./)
          endif
          if (gvz < 0.) then
            ivz = 0
            svz = (/0.,1.,0.,0./)
          endif
          if (ivx == nvx) svx = (/0.,1.,0.,0./)
          if (ivz == nvz) svz = (/0.,1.,0.,0./)

          do jz=iz-1,iz+jz2
            do jx=ix-1,ix+jx2
c             do jvz=ivz-1,ivz+jvz2
c               do jvx=ivx-1,ivx+jvx2

c                 --- Unrolling the loop over jvx made a big difference
c                 --- in the timing, cutting the time in half.
c                 --- Unrolling jvz only made a small change, a few percent.

c                 --- Removing common subexpressions helped by about 15%
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
                    if (ivx < nvx) then
                      grid(ivx+2,ivz+2,jx,jz)=grid(ivx+2,ivz+2,jx,jz)+svx(2)*sp2
                    endif
                  endif

c               enddo
c             enddo
            enddo
          enddo

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
      subroutine sumfondomainboundaries(nvx,nvz,nx,nz,nxguard,nzguard,
     &                                  fgrid,localpbounds)
      use Subtimerstop
      use Parallel
      integer(ISZ):: nvx,nvz,nx,nz,nxguard,nzguard
      real(kind=8):: fgrid(-1:nvx+1,-1:nvz+1,-nxguard:nx+1,-nzguard:nz+1)
      integer(ISZ):: localpbounds(0:5)

c Sum the distribution on domain boundaries.
c This assumes that each processor "owns" f from iz=0 to iz=nzlocal-1.
c Each is only responsible for sending out the source is owns.

#ifdef MPIPARALLEL

      real(kind=8),allocatable:: fgridrecv(:,:,:,:)
      integer(MPIISZ):: dims,nlocal(-2:1),starts(-2:1),sizes(-2:1)
      integer(MPIISZ):: iplower,ipupper,iparity
      integer(MPIISZ):: lowertype,uppertype
      integer(MPIISZ):: nn
      include "mpif.h"
      integer(MPIISZ):: mpi_comm
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE)
      integer(MPIISZ):: messid = 63
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      allocate(fgridrecv(-1:nvx+1,-1:nvz+1,-nxguard:nx+1,-nzguard:1))

      dims = 4
      nlocal = (/3+nvx,3+nvz,nxguard+2+nx,nzguard+2+nz/)
      starts = 0
      sizes = nlocal

      starts(1) = -1 + 1
      sizes(1) = 2 + nzguard
      call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                              MPI_ORDER_FORTRAN,
     &                              MPI_DOUBLE_PRECISION,
     &                              lowertype,mpierror)
      call MPI_TYPE_COMMIT(lowertype,mpierror)

      starts(1) = nz - 1 + 1
      sizes(1) = 2 + nzguard
      call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                              MPI_ORDER_FORTRAN,
     &                              MPI_DOUBLE_PRECISION,
     &                              uppertype,mpierror)
      call MPI_TYPE_COMMIT(uppertype,mpierror)

      mpi_comm = ppdecomp%mpi_comm_z

      iplower = ppdecomp%iprocgrid(2) - 1
      ipupper = ppdecomp%iprocgrid(2) + 1
      if (localpbounds(4) == 2 .and. iplower == -1) then
        iplower = ppdecomp%nprocgrid(2) - 1
      endif
      if (localpbounds(5) == 2 .and. ipupper == ppdecomp%nprocgrid(2)) then
        ipupper = 0
      endif

      do iparity=0,1
        if (mod(ppdecomp%iprocgrid(2),2) == iparity .and. iplower >= 0) then
          nn = 1
          call MPI_SEND(fgrid,nn,lowertype,
     &                  iplower,messid,mpi_comm,mpierror)
          call MPI_RECV(fgrid,nn,lowertype,
     &                  iplower,messid,mpi_comm,mpistatus,mpierror)
        elseif (mod(ppdecomp%iprocgrid(2),2) == 1-iparity .and. 
     &          ipupper < ppdecomp%nprocgrid(2)) then
          nn = (3+nvx)*(3+nvz)*(3+nx)*(3)
          call MPI_RECV(fgridrecv,nn,MPI_DOUBLE_PRECISION,
     &                  ipupper,messid,mpi_comm,mpistatus,mpierror)
          fgrid(:,:,:,nz-nzguard:nz+1) = fgrid(:,:,:,nz-nzguard:nz+1)+fgridrecv
          nn = 1
          call MPI_SEND(fgrid,nn,uppertype,
     &                  ipupper,messid,mpi_comm,mpierror)
        endif
      enddo

      call MPI_TYPE_FREE(lowertype,mpierror)
      call MPI_TYPE_FREE(uppertype,mpierror)

      deallocate(fgridrecv)

!$OMP MASTER
      if (ltoptimesubs) timesumfondomainboundaries = timesumfondomainboundaries + wtime() - substarttime
!$OMP END MASTER

#endif

      return
      end
c=============================================================================
