#include "top.h"
c=============================================================================
c@(#) File F3D_BFIELD.F, version $Revision: 1.38 $, $Date: 2007/06/06 23:10:26 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This contains routines to handle the B field, calculated from the current
c  density.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      module f3d_bfield_interfaces
      interface

      subroutine setj3d(j,j1d,np,xp,yp,zp,zgrid,uxp,uyp,uzp,gaminv,
     &                  q,wght,nw,wghtp,depos,nx,ny,nzlocal,dx,dy,dz,xmmin,ymmin,zmmin,
     &                  l2symtry,l4symtry,lcylindrical)
      integer(ISZ):: np,nw
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: zgrid,q,wght
      real(kind=8):: j(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: j1d(0:3*(1+nx)*(1+ny)*(1+nzlocal)-1)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8),target:: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8),target:: wghtp(nw)
      character(8):: depos
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical
      end subroutine setj3d

      subroutine setb3d(b,np,xp,yp,zp,zgrid,bx,by,bz,
     &                  nx,ny,nzlocal,dx,dy,dz,xmmin,ymmin,zmmin,
     &                  l2symtry,l4symtry,lcylindrical)
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: zgrid
      real(kind=8):: b(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8):: bx(np), by(np), bz(np)
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical
      end subroutine setb3d

      subroutine fetchafrompositions3d(a,np,xp,yp,zp,zgrid,ap,
     &                                 nx,ny,nzlocal,dx,dy,dz,xmmin,ymmin,zmmin,
     &                                 l2symtry,l4symtry,lcylindrical)
      use Subtimersf3d
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: a(0:2,-1:nx+1,-1:ny+1,-1:nzlocal+1)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8):: ap(0:2,np)
      real(kind=8):: zgrid
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical
      end subroutine fetchafrompositions3d

      end interface
      end module f3d_bfield_interfaces
c=============================================================================
      subroutine init_bfieldsolver(bfstype)
      use GlobalVars
      use BFieldGrid
      use InMesh3d
      use InGen3d, only: l2symtry,l4symtry,solvergeom,RZgeom
      use GridBoundary3d
      use Multigrid3d
      use Conductor3d, only: icndbndy

      integer(ISZ):: bfstype

      if (bfstype == -1) return
      if (bfstype == 12) return

      if (bfield%xmmin == 0) bfield%xmmin = xmmin
      if (bfield%xmmax == 0) bfield%xmmax = xmmax
      if (bfield%ymmin == 0) bfield%ymmin = ymmin
      if (bfield%ymmax == 0) bfield%ymmax = ymmax
      if (bfield%zmminlocal == 0) bfield%zmminlocal = zmminlocal
      if (bfield%zmmaxlocal == 0) bfield%zmmaxlocal = zmmaxlocal
      if (bfield%zmmin == 0) bfield%zmmin = zmmin
      if (bfield%zmmax == 0) bfield%zmmax = zmmax
      if (bfield%nx == 0) bfield%nx = nx
      if (bfield%ny == 0) bfield%ny = ny
      if (bfield%nzlocal == 0) bfield%nzlocal = nzlocal
      if (bfield%nz == 0) bfield%nz = nz
      if (bfield%dx == 0) bfield%dx = (bfield%xmmax-bfield%xmmin)/bfield%nx
      if (bfield%dy == 0 .and. bfield%ny > 0) bfield%dy = (bfield%ymmax-bfield%ymmin)/bfield%ny
      if (bfield%dz == 0) bfield%dz = (bfield%zmmax-bfield%zmmin)/bfield%nz
      call BFieldGridTypechange(bfield)

      bfield%icndbndy = icndbndy
      bfield%mgparam = mgparam
      bfield%mgmaxiters = mgmaxiters
      bfield%mgmaxlevels = mgmaxlevels
      if (bfield%mgtol(0) == 0) bfield%mgtol(0) = mgtol
      if (bfield%mgtol(1) == 0) bfield%mgtol(1) = mgtol
      if (bfield%mgtol(2) == 0) bfield%mgtol(2) = mgtol
      bfield%mgform = mgform
      bfield%downpasses = downpasses
      bfield%uppasses = uppasses

      bfield%bounds(0) = boundxy
      bfield%bounds(1) = boundxy
      bfield%bounds(2) = boundxy
      bfield%bounds(3) = boundxy
      bfield%bounds(4) = bound0
      bfield%bounds(5) = boundnz
      if (l2symtry) then
        bfield%bounds(2) = neumann
        if (boundxy == 2) bfield%bounds(3) = neumann
      else if (l4symtry) then
        bfield%bounds(0) = neumann
        bfield%bounds(2) = neumann
        if (boundxy == 2) bfield%bounds(1) = neumann
        if (boundxy == 2) bfield%bounds(3) = neumann
      endif

      bfield%lcylindrical = (solvergeom == RZgeom)
      if (bfield%lcylindrical) then
        bfield%dy = bfield%dx
        if (l4symtry .or. (bfield%lcylindrical .and. bfield%xmmin==0.)) then
          bfield%bounds(0) = neumann
        endif
        bfield%bounds(2) = dirichlet
        bfield%bounds(3) = dirichlet
      endif

      call bvp3d(1,bfstype)

      return
      end
c=============================================================================
      subroutine loadj3d(pgroup,ins_i,nps_i,is_i,lzero)
      use ParticleGroupmodule
      use GlobalVars
      use Subtimersf3d
      use InGen
      use InGen3d
      use Picglb
      use Picglb3d
      use Particles,Only: wpid
      use BFieldGrid
      use GridBoundary3d
      use FieldSolveAPI, Only: lzerorhofsapi
      use f3d_bfield_interfaces
      type(ParticleGroup):: pgroup
      integer(ISZ):: ins_i,nps_i,is_i
      logical(ISZ):: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- jp array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.

      integer(ISZ):: ins_u,nps_u
      integer(ISZ):: is1,is2
      integer(ISZ):: ip,ipmin,is
      real(kind=8):: swtmp,wptmp(1)
      real(kind=8):: substarttime,wtime
      real(kind=8),dimension(:,:,:,:),pointer::jpcopy
      real(kind=8),pointer:: pidtmp(:)
      integer(ISZ):: allocerror

      if (lf3dtimesubs) substarttime = wtime()

#ifdef MPIPARALLEL
c     --- Set the flag signalling that j has been changed and that
c     --- the parallel sums on boundaries will need to be done again.
      ljsummedonboundaries = .false.
#endif

      if (depos == 'none') return

      if (bfstype == 12) then
        lzerorhofsapi = lzero
        call execuser("bloadjregistered")
      elseif (fstype == 12 .and. bfstype < 0) then
        lzerorhofsapi = lzero
        call execuser("loadjregistered")
      elseif (bfstype >= 0) then

c       --- Ensure that the jp array is setup properly
        if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
          call setupbfieldsforparticles3d(pgroup%ns,pgroup%ndts,it,bfield,bfieldp)
        endif

c       --- zero jp if requested
        if (lzero) then
          if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
              call zeroarry(bfieldp%j,
     &                    3*(bfieldp%nx+1)*(bfieldp%ny+1)*(bfieldp%nzlocal+1))
          endif
        end if

c       --- set limits on loop over species
        if (is_i == -1) then
          is1 = 1
          is2 = pgroup%ns
        else
          is1 = is_i
          is2 = is_i
        endif

c       --- set initial limits from input
c       --- (will be changed if necessary in the loop)
        ins_u = ins_i
        nps_u = nps_i

c       --- loop over species
        do is=is1,is2

c         --- get loop limits for particles if needed
          if (ins_i == -1) ins_u = pgroup%ins(is)
          if (nps_i == -1) nps_u = pgroup%nps(is)

c         --- Scale the weight, sw, by the time step scale size. This only
c         --- makes sense for steady-state and slice modes. In time-dependent
c         --- mode, it is assumed that dtscale has not been changed from 1.
          swtmp = pgroup%sw(is)*pgroup%dtscale(is)

          if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
            ipmin = ins_u
            if (nps_u > 0) then
              if (pgroup%ndts(is-1)>1) then
c               --- This code would be better if it were in setj3d, but since
c               --- jp is passed into it, this must be done here.
                if (mod(it+1,pgroup%ndts(is-1))==0) then
c                  --- If this species is being advanced this step, then
c                  --- deposit its j into jptmp.
                   jpcopy => bfieldp%j
                   bfieldp%j => bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
                   bfieldp%j = 0.
                else
c                  --- If this species is not being advanced this step, then
c                  --- just add its saved jptmp into jp.
                   bfieldp%j = bfieldp%j + bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
                   cycle
                 end if
              end if
              if(wpid==0) then
                call setj3d(bfieldp%j,bfieldp%j,nps_u,
     &                      pgroup%xp(ipmin:ipmin+nps_u-1),
     &                      pgroup%yp(ipmin:ipmin+nps_u-1),
     &                      pgroup%zp(ipmin:ipmin+nps_u-1),zgrid,
     &                      pgroup%uxp(ipmin:ipmin+nps_u-1),
     &                      pgroup%uyp(ipmin:ipmin+nps_u-1),
     &                      pgroup%uzp(ipmin:ipmin+nps_u-1),
     &                      pgroup%gaminv(ipmin:ipmin+nps_u-1),
     &                      pgroup%sq(is),swtmp,0,wptmp,depos,
     &                      bfield%nx,bfield%ny,bfield%nzlocal,
     &                      bfield%dx,bfield%dy,bfield%dz,
     &                      bfield%xmmin,bfield%ymmin,bfield%zmminlocal,
     &                      l2symtry,l4symtry,bfield%lcylindrical)
              else
c               --- Due to compiler bug, the pid array needs to be temporarily
c               --- pointed to in order to be properly passed into setj3d.
c               --- When pid is passed directly into setj3d, it is just all
c               --- zeros.
                pidtmp => pgroup%pid(ipmin:ipmin+nps_u-1,wpid)
                call setj3d(bfieldp%j,bfieldp%j,nps_u,
     &                      pgroup%xp(ipmin:ipmin+nps_u-1),
     &                      pgroup%yp(ipmin:ipmin+nps_u-1),
     &                      pgroup%zp(ipmin:ipmin+nps_u-1),zgrid,
     &                      pgroup%uxp(ipmin:ipmin+nps_u-1),
     &                      pgroup%uyp(ipmin:ipmin+nps_u-1),
     &                      pgroup%uzp(ipmin:ipmin+nps_u-1),
     &                      pgroup%gaminv(ipmin:ipmin+nps_u-1),
     &                      pgroup%sq(is),swtmp,nps_u,pidtmp,depos,
     &                      bfield%nx,bfield%ny,bfield%nzlocal,
     &                      bfield%dx,bfield%dy,bfield%dz,
     &                      bfield%xmmin,bfield%ymmin,bfield%zmminlocal,
     &                      l2symtry,l4symtry,bfield%lcylindrical)
              endif
              if (pgroup%ndts(is-1)>1 .and. mod(it+1,pgroup%ndts(is-1))==0) then
c               --- If this species is being advanced this step, then restore
c               --- jp and copy its j into jp.
                bfieldp%j => jpcopy
                bfieldp%j = bfieldp%j + bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
              end if
            endif
          elseif(solvergeom==XZgeom) then
            print*,"Error: Self magnetic field not support with solvergeom=XZgeom"
          elseif(solvergeom==Zgeom) then
            print*,"Error: Self magnetic field not support with solvergeom=Zgeom"
          elseif(solvergeom==Rgeom) then
            print*,"Error: Self magnetic field not support with solvergeom=Rgeom"
          elseif(solvergeom==AMRgeom) then
            print*,"Error: Self magnetic field not support with solvergeom=AMRgeom"
          end if
        enddo


c       --- For parallel version, each processor sends j to neighboring
c       --- processors whose field solve region overlap its particle region.
        if (lzero) call getjforfieldsolve()

c       --- enforce periodicity
c       --- This is done on the j array, rather than jp, since j will
c       --- be used for the field solve. Also, in some cases, the jp array
c       --- will not cover the full axial extent of the system, so 
c       --- periodicity can only be enforced on j which does cover the full
c       --- system.
        if (lzero) then
          call perj3d(bfield,bound0,boundxy)
        endif

      endif

!$OMP MASTER
      if (lf3dtimesubs) timeloadj3d = timeloadj3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setaboundaries3d(a,nx,ny,nzlocal,
     &                            zmminlocal,zmmaxlocal,zmmin,zmmax,
     &                            bounds,lcylindrical,llongitudinal)
      use Subtimersf3d
      use GlobalVars
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: a(0:2,-1:nx+1,-1:ny+1,-1:nzlocal+1)
      real(kind=8):: zmminlocal,zmmaxlocal,zmmin,zmmax
      integer(ISZ):: bounds(0:5)
      logical(ISZ):: lcylindrical,llongitudinal

c  Sets the slices on the exterior of A, depending on the boundary conditions.
c  When llongitudinal is true, set the longitudinal boundaries, otherwise
c  skip that code block.

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- With the cylindrical solver, the guard cells are already set.
      if (.not. lcylindrical) then
        if (bounds(0) == periodic .or. bounds(1) == periodic) then
          a(:,-1,:,:)   = a(:,nx-1,:,:)
          a(:,nx,:,:)   = a(:,0,:,:)
          a(:,nx+1,:,:) = a(:,1,:,:)
        else
          if (bounds(0) == neumann) then
            a(:,-1,:,:) = a(:,1,:,:)
          elseif (bounds(0) == dirichlet) then
            a(:,-1,:,:) = 2.*a(:,0,:,:) - a(:,1,:,:)
          endif
          if (bounds(1) == neumann) then
            a(:,nx+1,:,:) = a(:,nx-1,:,:)
          elseif (bounds(1) == dirichlet) then
            a(:,nx+1,:,:) = 2.*a(:,nx,:,:) - a(:,nx-1,:,:)
          endif
        endif
      endif

      if (lcylindrical) then
c       --- In this case, actually ny == 0. This ensures that d/dtheta is zero.
        a(:,:,-1,:) = a(:,:,0,:)
        a(:,:,ny+1,:) = a(:,:,ny,:)
      else
        if (bounds(2) == periodic .or. bounds(3) == periodic) then
          a(:,:,-1,:)   = a(:,:,ny-1,:)
          a(:,:,ny,:)   = a(:,:,0,:)
          a(:,:,ny+1,:) = a(:,:,1,:)
        else
          if (bounds(2) == neumann) then
            a(:,:,-1,:) = a(:,:,1,:)
          elseif (bounds(2) == dirichlet) then
            a(:,:,-1,:) = 2.*a(:,:,0,:) - a(:,:,1,:)
          endif
          if (bounds(3) == neumann) then
            a(:,:,ny+1,:) = a(:,:,ny-1,:)
          elseif (bounds(3) == dirichlet) then
            a(:,:,ny+1,:) = 2.*a(:,:,ny,:) - a(:,:,ny-1,:)
          endif
        endif
      endif

c     --- With the cylindrical solver, the guard cells are already set.
      if (llongitudinal .and. .not. lcylindrical) then
        if (bounds(4) == periodic .or. bounds(5) == periodic) then
#ifdef MPIPARALLEL
          call perpot3d_slave(a,3,nx,ny,nzlocal,1,1)
#else
          a(:,:,:,-1)   = a(:,:,:,nzlocal-1)
          a(:,:,:,nzlocal)   = a(:,:,:,0)
          a(:,:,:,nzlocal+1) = a(:,:,:,1)
#endif
        else
          if (zmminlocal == zmmin) then
            if (bounds(4) == neumann) then
              a(:,:,:,-1) = a(:,:,:,1)
            elseif (bounds(4) == dirichlet) then
              a(:,:,:,-1) = 2.*a(:,:,:,0) - a(:,:,:,1)
            endif
          endif
          if (zmmaxlocal == zmmax) then
            if (bounds(5) == neumann) then
              a(:,:,:,nzlocal+1) = a(:,:,:,nzlocal-1)
            elseif (bounds(5) == dirichlet) then
              a(:,:,:,nzlocal+1) = 2.*a(:,:,:,nzlocal) - a(:,:,:,nzlocal-1)
            endif
          endif
        endif
      endif

!$OMP MASTER
      if (lf3dtimesubs) timepera3d = timepera3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine perj3d(bfield,bound0,boundxy)
      use BFieldGridTypemodule
      use GlobalVars
      use Subtimersf3d
      use InGen3d
      type(BFieldGridType):: bfield
      integer(ISZ):: bound0,boundxy

c  Sums the first and last slices of j for periodicity
c  and puts the result into both slices.

      real(kind=8),pointer:: j(:,:,:,:)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (solvergeom==AMRgeom) return

      j => bfield%j

c  Enforce transverse periodicity if requested.
      if (solvergeom==XYZgeom) then
        if (boundxy==periodic) then
          j(:,:,0,:) = j(:,:,0,:) + j(:,:,bfield%ny,:)
          j(:,0,:,:) = j(:,0,:,:) + j(:,bfield%nx,:,:)
          j(:,:,bfield%ny,:) = j(:,:,0,:)
          j(:,bfield%nx,:,:) = j(:,0,:,:)
        endif
        if (boundxy==neumann) then
          if (.not. (l2symtry .or. l4symtry)) j(:,:,0,:) = 2.*j(:,:,0,:)
          if (.not. l4symtry) j(:,0,:,:) = 2.*j(:,0,:,:)
          j(:,:,bfield%ny,:) = 2.*j(:,:,bfield%ny,:)
          j(:,bfield%nx,:,:) = 2.*j(:,bfield%nx,:,:)
        endif
      endif

      if (solvergeom==XYZgeom .and. bound0==periodic) then
#ifdef MPIPARALLEL
        call persource3d_slave(j,3,bfield%nx,bfield%ny,bfield%nzlocal)
#else
        j(:,:,:,0)  = j(:,:,:,0) + j(:,:,:,bfield%nzlocal)
        j(:,:,:,bfield%nzlocal) = j(:,:,:,0)
#endif
      endif

!$OMP MASTER
      if (lf3dtimesubs) timeperj3d = timeperj3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setb3d(b,np,xp,yp,zp,zgrid,bx,by,bz,
     &                  nx,ny,nzlocal,dx,dy,dz,xmmin,ymmin,zmmin,
     &                  l2symtry,l4symtry,lcylindrical)
      use Subtimersf3d
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: zgrid
      real(kind=8):: b(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8):: bx(np), by(np), bz(np)
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

c  Sets magnetic field for particles

c     --- Temp arrays to hold particle data
c     --- These are needed when lcylindrical is true, in which case x = r and
c     --- y = 0.
      real(kind=8),pointer:: x(:),y(:),z(:),bxtemp(:)

      integer(ISZ):: ip,i,j,k
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      real(kind=8):: sx,sy

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (lcylindrical) then
c       --- Allocate and set the temporary particle arrays which will hold
c       --- the radius.
        allocate(x(np),y(np),bxtemp(np))
        x = sqrt(xp**2 + yp**2)
        y = 0.
      else
        x => xp
        y => yp
      endif
      z => zp

c  Evaluation of B, vectorized over particles
      tdxi = 1./(2.*dx)
      tdyi = 1./(2.*dy)
      tdzi = 1./(2.*dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      if ((.not. (l2symtry .or. l4symtry)) .or. lcylindrical) then
        do ip = 1, np

          i = (x(ip) - xmmin) * dxi
          j = (y(ip) - ymmin) * dyi
          k = (z(ip) - zgrid - zmmin) * dzi

          u1 = (x(ip) - xmmin) * dxi - i
          v1 = (y(ip) - ymmin) * dyi - j
          w1 = (z(ip) - zgrid - zmmin) * dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

          bx(ip) = u0*v0*w0*b(0,i  ,j  ,k  )
     &           + u1*v0*w0*b(0,i+1,j  ,k  )
     &           + u0*v1*w0*b(0,i  ,j+1,k  )
     &           + u1*v1*w0*b(0,i+1,j+1,k  )
     &           + u0*v0*w1*b(0,i  ,j  ,k+1)
     &           + u1*v0*w1*b(0,i+1,j  ,k+1)
     &           + u0*v1*w1*b(0,i  ,j+1,k+1)
     &           + u1*v1*w1*b(0,i+1,j+1,k+1)

          by(ip) = u0*v0*w0*b(1,i  ,j  ,k  )
     &           + u1*v0*w0*b(1,i+1,j  ,k  )
     &           + u0*v1*w0*b(1,i  ,j+1,k  )
     &           + u1*v1*w0*b(1,i+1,j+1,k  )
     &           + u0*v0*w1*b(1,i  ,j  ,k+1)
     &           + u1*v0*w1*b(1,i+1,j  ,k+1)
     &           + u0*v1*w1*b(1,i  ,j+1,k+1)
     &           + u1*v1*w1*b(1,i+1,j+1,k+1)

          bz(ip) = u0*v0*w0*b(2,i  ,j  ,k  )
     &           + u1*v0*w0*b(2,i+1,j  ,k  )
     &           + u0*v1*w0*b(2,i  ,j+1,k  )
     &           + u1*v1*w0*b(2,i+1,j+1,k  )
     &           + u0*v0*w1*b(2,i  ,j  ,k+1)
     &           + u1*v0*w1*b(2,i+1,j  ,k+1)
     &           + u0*v1*w1*b(2,i  ,j+1,k+1)
     &           + u1*v1*w1*b(2,i+1,j+1,k+1)

        enddo

      else

c       --- Set the signs of the B field for particles on negative side of
c       --- the axis of symmetry.
        sy = -1.
        sx = 1.
        if (l4symtry) sx = -1.

c       --- special loop symmetry is used
        do ip = 1, np

          i = (abs(x(ip)) - xmmin)*dxi
          j = (abs(y(ip)) - ymmin)*dyi
          k = (z(ip) - zgrid - zmmin)*dzi

          u1 = (abs(x(ip)) - xmmin)*dxi - i
          v1 = (abs(y(ip)) - ymmin)*dyi - j
          w1 = (z(ip) - zgrid - zmmin)*dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Adjust sign of B field for appropriate quadrant.
          xsign = +1.
          ysign = +1.
          if (x(ip) < 0.) xsign = sx
          if (y(ip) < 0.) ysign = sy

          bx(ip) = xsign*(u0*v0*w0*b(0,i  ,j  ,k  )
     &                  + u1*v0*w0*b(0,i+1,j  ,k  )
     &                  + u0*v1*w0*b(0,i  ,j+1,k  )
     &                  + u1*v1*w0*b(0,i+1,j+1,k  )
     &                  + u0*v0*w1*b(0,i  ,j  ,k+1)
     &                  + u1*v0*w1*b(0,i+1,j  ,k+1)
     &                  + u0*v1*w1*b(0,i  ,j+1,k+1)
     &                  + u1*v1*w1*b(0,i+1,j+1,k+1))

          by(ip) = ysign*(u0*v0*w0*b(1,i  ,j  ,k  )
     &                  + u1*v0*w0*b(1,i+1,j  ,k  )
     &                  + u0*v1*w0*b(1,i  ,j+1,k  )
     &                  + u1*v1*w0*b(1,i+1,j+1,k  )
     &                  + u0*v0*w1*b(1,i  ,j  ,k+1)
     &                  + u1*v0*w1*b(1,i+1,j  ,k+1)
     &                  + u0*v1*w1*b(1,i  ,j+1,k+1)
     &                  + u1*v1*w1*b(1,i+1,j+1,k+1))

          bz(ip) =        u0*v0*w0*b(2,i  ,j  ,k  )
     &                  + u1*v0*w0*b(2,i+1,j  ,k  )
     &                  + u0*v1*w0*b(2,i  ,j+1,k  )
     &                  + u1*v1*w0*b(2,i+1,j+1,k  )
     &                  + u0*v0*w1*b(2,i  ,j  ,k+1)
     &                  + u1*v0*w1*b(2,i+1,j  ,k+1)
     &                  + u0*v1*w1*b(2,i  ,j+1,k+1)
     &                  + u1*v1*w1*b(2,i+1,j+1,k+1)

        enddo

      endif

      if (lcylindrical) then
c       --- Transform Br and Btheta into Bx and By
        y = atan2(yp,xp)
        bxtemp = bx
        bx = bxtemp*cos(y) - by*sin(y)
        by = bxtemp*sin(y) + by*cos(y)
        deallocate(x,y,bxtemp)
      endif

!$OMP MASTER
      if (lf3dtimesubs) timesetb3d = timesetb3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fetchafrompositions3d(a,np,xp,yp,zp,zgrid,ap,
     &                                 nx,ny,nzlocal,dx,dy,dz,xmmin,ymmin,zmmin,
     &                                 l2symtry,l4symtry,lcylindrical)
      use Subtimersf3d
      integer(ISZ):: np
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: a(0:2,-1:nx+1,-1:ny+1,-1:nzlocal+1)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8):: ap(0:2,np)
      real(kind=8):: zgrid
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

c  Gets the magnetic vector potential

c     --- Temp arrays to hold particle data
c     --- These are needed when lcylindrical is true, in which case x = r and
c     --- y = 0.
      real(kind=8),pointer:: x(:),y(:),z(:),axtemp(:)

      integer(ISZ):: ip,i,j,k
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      real(kind=8):: sx,sy

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (lcylindrical) then
c       --- Allocate and set the temporary particle arrays which will hold
c       --- the radius.
        allocate(x(np),y(np),axtemp(np))
        x = sqrt(xp**2 + yp**2)
        y = 0.
      else
        x => xp
        y => yp
      endif
      z => zp

c  Evaluation of A, vectorized over particles
      tdxi = 1./(2.*dx)
      tdyi = 1./(2.*dy)
      tdzi = 1./(2.*dz)
      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      if ((.not. (l2symtry .or. l4symtry)) .or. lcylindrical) then
        do ip = 1, np

          i = (x(ip) - xmmin) * dxi
          j = (y(ip) - ymmin) * dyi
          k = (z(ip) - zgrid - zmmin) * dzi

          u1 = (x(ip) - xmmin) * dxi - i
          v1 = (y(ip) - ymmin) * dyi - j
          w1 = (z(ip) - zgrid - zmmin) * dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

          ap(0,ip) = u0*v0*w0*a(0,i  ,j  ,k  )
     &             + u1*v0*w0*a(0,i+1,j  ,k  )
     &             + u0*v1*w0*a(0,i  ,j+1,k  )
     &             + u1*v1*w0*a(0,i+1,j+1,k  )
     &             + u0*v0*w1*a(0,i  ,j  ,k+1)
     &             + u1*v0*w1*a(0,i+1,j  ,k+1)
     &             + u0*v1*w1*a(0,i  ,j+1,k+1)
     &             + u1*v1*w1*a(0,i+1,j+1,k+1)

          ap(1,ip) = u0*v0*w0*a(1,i  ,j  ,k  )
     &             + u1*v0*w0*a(1,i+1,j  ,k  )
     &             + u0*v1*w0*a(1,i  ,j+1,k  )
     &             + u1*v1*w0*a(1,i+1,j+1,k  )
     &             + u0*v0*w1*a(1,i  ,j  ,k+1)
     &             + u1*v0*w1*a(1,i+1,j  ,k+1)
     &             + u0*v1*w1*a(1,i  ,j+1,k+1)
     &             + u1*v1*w1*a(1,i+1,j+1,k+1)

          ap(2,ip) = u0*v0*w0*a(2,i  ,j  ,k  )
     &             + u1*v0*w0*a(2,i+1,j  ,k  )
     &             + u0*v1*w0*a(2,i  ,j+1,k  )
     &             + u1*v1*w0*a(2,i+1,j+1,k  )
     &             + u0*v0*w1*a(2,i  ,j  ,k+1)
     &             + u1*v0*w1*a(2,i+1,j  ,k+1)
     &             + u0*v1*w1*a(2,i  ,j+1,k+1)
     &             + u1*v1*w1*a(2,i+1,j+1,k+1)

        enddo

      else

c       --- Set the signs of the B field for particles on negative side of
c       --- the axis of symmetry.
        sy = -1.
        sx = 1.
        if (l4symtry) sx = -1.

c       --- special loop symmetry is used
        do ip = 1, np

          i = (abs(x(ip)) - xmmin)*dxi
          j = (abs(y(ip)) - ymmin)*dyi
          k = (z(ip) - zgrid - zmmin)*dzi

          u1 = (abs(x(ip)) - xmmin)*dxi - i
          v1 = (abs(y(ip)) - ymmin)*dyi - j
          w1 = (z(ip) - zgrid - zmmin)*dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Adjust sign of B field for appropriate quadrant.
          xsign = +1.
          ysign = +1.
          if (x(ip) < 0.) xsign = sx
          if (y(ip) < 0.) ysign = sy

          ap(0,ip) = xsign*(u0*v0*w0*a(0,i  ,j  ,k  )
     &                    + u1*v0*w0*a(0,i+1,j  ,k  )
     &                    + u0*v1*w0*a(0,i  ,j+1,k  )
     &                    + u1*v1*w0*a(0,i+1,j+1,k  )
     &                    + u0*v0*w1*a(0,i  ,j  ,k+1)
     &                    + u1*v0*w1*a(0,i+1,j  ,k+1)
     &                    + u0*v1*w1*a(0,i  ,j+1,k+1)
     &                    + u1*v1*w1*a(0,i+1,j+1,k+1))

          ap(1,ip) = ysign*(u0*v0*w0*a(1,i  ,j  ,k  )
     &                    + u1*v0*w0*a(1,i+1,j  ,k  )
     &                    + u0*v1*w0*a(1,i  ,j+1,k  )
     &                    + u1*v1*w0*a(1,i+1,j+1,k  )
     &                    + u0*v0*w1*a(1,i  ,j  ,k+1)
     &                    + u1*v0*w1*a(1,i+1,j  ,k+1)
     &                    + u0*v1*w1*a(1,i  ,j+1,k+1)
     &                    + u1*v1*w1*a(1,i+1,j+1,k+1))

          ap(2,ip) =        u0*v0*w0*a(2,i  ,j  ,k  )
     &                    + u1*v0*w0*a(2,i+1,j  ,k  )
     &                    + u0*v1*w0*a(2,i  ,j+1,k  )
     &                    + u1*v1*w0*a(2,i+1,j+1,k  )
     &                    + u0*v0*w1*a(2,i  ,j  ,k+1)
     &                    + u1*v0*w1*a(2,i+1,j  ,k+1)
     &                    + u0*v1*w1*a(2,i  ,j+1,k+1)
     &                    + u1*v1*w1*a(2,i+1,j+1,k+1)

        enddo

      endif

      if (lcylindrical) then
c       --- Transform Br and Btheta into Bx and By
        y = atan2(yp,xp)
        axtemp = ap(0,:)
        ap(0,:) = axtemp*cos(y) - ap(1,:)*sin(y)
        ap(1,:) = axtemp*sin(y) + ap(1,:)*cos(y)
        deallocate(x,y,axtemp)
      endif

!$OMP MASTER
      if (lf3dtimesubs) timefetchafrompositions3d = timefetchafrompositions3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=========================================================================== 
      subroutine getbfroma3d(a,b,nx,ny,nzlocal,dx,dy,dz,xmmin,
     &                       lcylindrical,lusevectorpotential)
      use GlobalVars
      use Subtimersf3d
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: a(0:2,-1:nx+1,-1:ny+1,-1:nzlocal+1)
      real(kind=8):: b(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin
      logical(ISZ):: lcylindrical,lusevectorpotential

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (lusevectorpotential) then

c       --- Calculate B via finite differences of A.
        call curl3d(a,b,nx,ny,nzlocal,
     &              dx,dy,dz,
     &              xmmin,lcylindrical,
     &              1,1,0,0)

      else

c       --- The solver calculated B directly, so just do a copy.
        b = a(:,0:nx,0:ny,0:nzlocal)

      endif

!$OMP MASTER
      if (lf3dtimesubs) timegetbfroma3d = timegetbfroma3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=========================================================================== 
      subroutine curl3d(a,b,nx,ny,nzlocal,dx,dy,dz,xmmin,lcylindrical,
     &                  adelx,adelz,bdelx,bdelz)
      integer(ISZ):: nx,ny,nzlocal,adelx,adelz,bdelx,bdelz
      real(kind=8):: a(0:2,-adelx:nx+adelx,-adelx:ny+adelx,-adelz:nzlocal+adelz)
      real(kind=8):: b(0:2,-bdelx:nx+bdelx,-bdelx:ny+bdelx,-bdelz:nzlocal+bdelz)
      real(kind=8):: dx,dy,dz,xmmin
      logical(ISZ):: lcylindrical

c Calculate the curl of A using finite differences.

      integer(ISZ):: ix,iy,iz
      real(kind=8):: tdxi,tdyi,tdzi,r,ri

      tdxi = 0.5/dx
      tdyi = 0.5/dy
      tdzi = 0.5/dz

c     --- Do the calculation. In this case, guard cells have been added to
c     --- A since having special cases for each B component on each
c     --- boundary for each boundary condition begins to get very complicated.
      if (ny > 0) then
        do iz=-adelz+1,nzlocal+adelz-1
          do iy=-adelx+1,ny+adelx-1
            do ix=-adelx+1,nx+adelx-1
              b(0,ix,iy,iz) = + tdyi*(a(2,ix  ,iy+1,iz  ) - a(2,ix  ,iy-1,iz  ))
     &                        - tdzi*(a(1,ix  ,iy  ,iz+1) - a(1,ix  ,iy  ,iz-1))
              b(1,ix,iy,iz) = + tdzi*(a(0,ix  ,iy  ,iz+1) - a(0,ix  ,iy  ,iz-1))
     &                        - tdxi*(a(2,ix+1,iy  ,iz  ) - a(2,ix-1,iy  ,iz  ))
              b(2,ix,iy,iz) = + tdxi*(a(1,ix+1,iy  ,iz  ) - a(1,ix-1,iy  ,iz  ))
     &                        - tdyi*(a(0,ix  ,iy+1,iz  ) - a(0,ix  ,iy-1,iz  ))
            enddo
          enddo
        enddo
      else
        iy = 0
        do iz=-adelz+1,nzlocal+adelz-1
          do ix=-adelx+1,nx+adelx-1
            b(0,ix,iy,iz) = + 0.
     &                      - tdzi*(a(1,ix  ,iy  ,iz+1) - a(1,ix  ,iy  ,iz-1))
            b(1,ix,iy,iz) = + tdzi*(a(0,ix  ,iy  ,iz+1) - a(0,ix  ,iy  ,iz-1))
     &                      - tdxi*(a(2,ix+1,iy  ,iz  ) - a(2,ix-1,iy  ,iz  ))
            b(2,ix,iy,iz) = + tdxi*(a(1,ix+1,iy  ,iz  ) - a(1,ix-1,iy  ,iz  ))
     &                      - 0.
          enddo
        enddo
      endif

c     --- When lcylindrical is true, an extra term is needed for Bz
c     ---      1 d               d Atheta   Atheta
c     --- Bz = - -- (r Atheta) = -------- + ------
c     ---      r dr                 dr        r
c     --- The derivative is calculated above. Note that the d/dtheta terms
c     --- which appear in Br and Bz are still ignored, assuming axisymmetry.
c     --- Because of the way the guard cells are set for the y plane,
c     --- the derivatives there calculated above are zero.
      if (lcylindrical) then
        do ix=-adelx+1,nx+adelx-1
          r = (xmmin + dx*ix)
          if (r > 0.) then
            ri = 1./r
            b(2,ix,:,:)=b(2,ix,:,:)+a(1,ix,-bdelx:ny+bdelx,-bdelz:nzlocal+bdelz)*ri
          else
c           --- This can be justified by writing out the del cross A in
c           --- Cartesian coordinates, and using
c           --- a(1,1,0,iz) = a(0,0,-1,iz) = -a(1,-1,0,iz) = -a(0,0,+1,iz)
            b(2,ix,:,:) = 2.*b(2,ix,:,:)
          endif
        enddo
      endif

      return
      end
c=========================================================================== 
      subroutine getefroma3d(e,nx,ny,nzlocal,dz,vframe)
      use GlobalVars
      use BFieldGrid
      use Subtimersf3d
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: e(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: dz,vframe

      integer(ISZ):: ix,iy,iz
      real(kind=8):: substarttime,wtime,tdzi
      if (lf3dtimesubs) substarttime = wtime()

      tdzi = 0.5*vframe/dz

      do iz=0,nzlocal
        do iy=0,ny
          do ix=0,nx
            e(0,ix,iy,iz) = e(0,ix,iy,iz)+tdzi*(bfield%a(0,ix ,iy ,iz+1) - bfield%a(0,ix ,iy ,iz-1))
            e(1,ix,iy,iz) = e(1,ix,iy,iz)+tdzi*(bfield%a(1,ix ,iy ,iz+1) - bfield%a(1,ix ,iy ,iz-1))
            e(2,ix,iy,iz) = e(2,ix,iy,iz)+tdzi*(bfield%a(2,ix ,iy ,iz+1) - bfield%a(2,ix ,iy ,iz-1))
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine getanalyticbtheta(b,j,nx,ny,nzlocal,dx,xmmin)
      use Constant
      use BFieldGridTypemodule
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: b(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: j(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: dx,xmmin

      integer(ISZ):: ir
      real(kind=8):: r
      real(kind=8):: izi(0:nzlocal)

c     --- If xmmin > 0, this assumes that there is no current below xmmin.
c     --- Note that a future version could extract a current below xmmin
c     --- by inspecting a nonzero Btheta at xmmin.
      izi = 0.
      b(1,0,0,:) = 0.

      do ir=0,nx-1
c       --- The total current is calculated assuming linear interpolation
c       --- between grid points, i.e. the trapezoidal rule.
        izi = izi + dx**2*(j(2,ir  ,0,:)*(ir/2.+1./6.) +
     &                            j(2,ir+1,0,:)*(ir/2.+1./3.))
        r = xmmin + (ir+1)*dx
        b(1,ir+1,0,:) = mu0*izi/r
      enddo

      return
      end
c=============================================================================
      subroutine setj3d(j,j1d,np,xp,yp,zp,zgrid,uxp,uyp,uzp,gaminv,
     &                  q,w,nw,wghtp,depos,nx,ny,nzlocal,dx,dy,dz,xmmin,ymmin,zmmin,
     &                  l2symtry,l4symtry,lcylindrical)
      use GlobalVars
      use Constant
      use Subtimersf3d
      integer(ISZ):: np,nw
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: zgrid,q,w
      real(kind=8):: j(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: j1d(0:3*(1+nx)*(1+ny)*(1+nzlocal)-1)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8),target:: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8),target:: wghtp(nw)
      character(8):: depos
      real(kind=8):: dx,dy,dz
      real(kind=8):: xmmin,ymmin,zmmin
      logical(ISZ):: l2symtry,l4symtry,lcylindrical

c  Sets current density

c  Algorithm notes: j array is dimensioned (3,0:nx,0:ny,0:nzlocal) outside,
c  but is made one dimensional in this routine
c  so cell index into 1d j array for vectorized deposition is:
c     3*(i + j*(nx+1) + k*(nx+1)*(ny+1))
c  In each case,
c     j(d,i  ,j  ,k  ) = j(d,i  ,j  ,k  ) + u0 * v0 * w0 * g
c     j(d,i+1,j  ,k  ) = j(d,i+1,j  ,k  ) + u1 * v0 * w0 * g
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  j is not zeroed here (to allow handling of blocks of particles
c  at a time)

c     --- For vectorized algorithm
      integer(ISZ):: moff(0:7)
      integer(ISZ),allocatable:: indx(:,:)
      real(kind=8),allocatable:: s(:,:),v(:,:)
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ),allocatable:: ii(:), jj(:), kk(:)
c     --- For both algorithms
c     --- Work array holding q/cell volume, the charge density per
c     --- real particle. This is primarily needed for the RZ version
c     --- since the cell volume there has radial dependence.
      real(kind=8):: cdens(0:nx)
c     --- Work array holding q*w/cell volume, the charge density per
c     --- simulation particle.
      real(kind=8):: cdensw(0:nx)
c     --- Temp arrays to hold particle data
c     --- These are needed when lcylindrical is true, in which case x=r, y=0,
c     --- and vx = vr, and vy = vtheta.
      real(kind=8),pointer:: x(:),y(:),z(:),ux(:),uy(:),uz(:),gi(:),wght(:)

      integer(ISZ):: ipmin,nptmp,ix,iy,iz,ind0,m,ir
      real(kind=8):: g,gw,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxf,gyf
      real(kind=8):: vv(0:2)
      real(kind=8):: xsign,ysign
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c  Set up offset array for vectorized deposition:

      moff(0) = 0
      moff(1) = 1
      moff(2) = nx+1
      moff(3) = nx+2
      moff(4) = (nx+1)*(ny+1)
      moff(5) = (nx+1)*(ny+1)+1
      moff(6) = (nx+1)*(ny+1)+nx+1
      moff(7) = (nx+1)*(ny+1)+nx+2

      dxi = 1./dx
      dyi = 1./dy
      dzi = 1./dz

      if (lcylindrical) then
        if (xmmin == 0.) then
c         --- The factor 0.75 corrects for overdeposition due to linear
c         --- weighting (for uniform distribution)
c         --- see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
c         --- and Verboncoeur, J. of Comp. Phys.,
          cdens(0) = 0.75/(pi*(0.5*0.5*dx*dx*dz))
        else
          cdens(0) = 1./(2.*pi*(xmmin)*dx*dz)
        endif
        do ix = 1,nx
          cdens(ix) = 1./(2.*pi*(ix*dx+xmmin)*dx*dz)
        enddo
        cdens = cdens*q*w
      else
        g = w*q/(dx*dy*dz)
        if (l2symtry) then
c         --- The particle weight is reduced by a factor of 2 except near the
c         --- transverse boundaries.
          g = g*0.5
        elseif (l4symtry) then
c         --- The particle weight is reduced by a factor of 4 except near the
c         --- transverse boundaries.
          g = g*0.25
        endif
        cdens = g
      endif

      if (nw == 0) then
        gw = g
        cdensw = cdens
      endif

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ipmin,nptmp,i,j,k,u1,u0,v1,v0,w1,w0,ir,ip,ind0,indx,
!$OMP&                 gyf,gxf,s,m,ii,jj,kk,x,y,z,ux,uy,uz)
!$OMP&FIRSTPRIVATE(np)

      if (lcylindrical) then
c       --- Allocate temporary particle arrays which will hold
c       --- the radius.
c       --- This is done in the parallel block since each thread must allocate
c       --- its own work arrays.
        allocate(x(nparpgrp),y(nparpgrp))
        allocate(ux(nparpgrp),uy(nparpgrp))
      endif

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)

c       --- Setup temporary particle arrays
        if (lcylindrical) then
          x(1:nptmp) = sqrt(xp(ipmin:ipmin+nptmp-1)**2 +
     &                      yp(ipmin:ipmin+nptmp-1)**2)
          y(1:nptmp) = atan2(yp(ipmin:ipmin+nptmp-1),
     &                       xp(ipmin:ipmin+nptmp-1))
          ux(1:nptmp) = uxp(ipmin:ipmin+nptmp-1)*cos(y(1:nptmp)) +
     &                  uyp(ipmin:ipmin+nptmp-1)*sin(y(1:nptmp))
          uy(1:nptmp) = -uxp(ipmin:ipmin+nptmp-1)*sin(y(1:nptmp)) +
     &                   uyp(ipmin:ipmin+nptmp-1)*cos(y(1:nptmp))
          y(1:nptmp) = 0.
        else
          x => xp(ipmin:ipmin+nptmp-1)
          y => yp(ipmin:ipmin+nptmp-1)
          ux => uxp(ipmin:ipmin+nptmp-1)
          uy => uyp(ipmin:ipmin+nptmp-1)
        endif
        z => zp(ipmin:ipmin+nptmp-1)
        uz => uzp(ipmin:ipmin+nptmp-1)
        gi => gaminv(ipmin:ipmin+nptmp-1)
        if (nw == np) wght => wghtp(ipmin:ipmin+nptmp-1)

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

      allocate(indx(0:7,0:nparpgrp-1))
      allocate(s(0:7,0:nparpgrp-1))
      allocate(v(0:7,0:nparpgrp-1))

c     --- vectorized loop to compute indices, weights
      if (l2symtry .and. .not. lcylindrical) then
c       --- special loop for 2-fold symmetry
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        do ir = 1,nptmp
           ix  = int((x(ir) - xmmin) * dxi)
           u1 =      (x(ir) - xmmin) * dxi - ix
           u0 = 1. - u1
           iy  = int((abs(y(ir)) - ymmin)*dyi)
           v1 =      (abs(y(ir)) - ymmin)*dyi - iy
           v0 = 1. - v1
           iz  = int((z(ir) - zgrid - zmmin) * dzi)
           w1 =      (z(ir) - zgrid - zmmin) * dzi - iz
           w0 = 1. - w1
           ind0 = ix + iy*(nx+1) + iz*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gyf = 1.
           if (iy == 0 .and. ymmin == 0.) gyf = 2.
           if (nw == np) gw = g*wght(ir)
           s(0,ir) = u0 * v0 * w0 * gw*gyf
           s(1,ir) = u1 * v0 * w0 * gw*gyf
           s(2,ir) = u0 * v1 * w0 * gw
           s(3,ir) = u1 * v1 * w0 * gw
           s(4,ir) = u0 * v0 * w1 * gw*gyf
           s(5,ir) = u1 * v0 * w1 * gw*gyf
           s(6,ir) = u0 * v1 * w1 * gw
           s(7,ir) = u1 * v1 * w1 * gw
           ysign = 1.
           if (y(ir) < 0.) ysign = -1.
           v(0,ir) = ux(ir)*gi(ir)*ysign
           v(1,ir) = uy(ir)*gi(ir)
           v(2,ir) = uz(ir)*gi(ir)
        enddo
      elseif (l4symtry .and. .not. lcylindrical) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        do ir = 1,nptmp
           ix  = int((abs(x(ir)) - xmmin)*dxi)
           u1 =      (abs(x(ir)) - xmmin)*dxi - ix
           u0 = 1. - u1
           iy  = int((abs(y(ir)) - ymmin)*dyi)
           v1 =      (abs(y(ir)) - ymmin)*dyi - iy
           v0 = 1. - v1
           iz  = int((z(ir) - zgrid - zmmin) * dzi)
           w1 =      (z(ir) - zgrid - zmmin) * dzi - iz
           w0 = 1. - w1
           ind0 = ix + iy*(nx+1) + iz*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gxf = 1.
           gyf = 1.
           if (ix == 0 .and. xmmin == 0.) gxf = 2.
           if (iy == 0 .and. ymmin == 0.) gyf = 2.
           if (nw == np) gw = g*wght(ir)
           s(0,ir) = u0 * v0 * w0 * gw*gxf*gyf
           s(1,ir) = u1 * v0 * w0 * gw*gyf
           s(2,ir) = u0 * v1 * w0 * gw*gxf
           s(3,ir) = u1 * v1 * w0 * gw
           s(4,ir) = u0 * v0 * w1 * gw*gxf*gyf
           s(5,ir) = u1 * v0 * w1 * gw*gyf
           s(6,ir) = u0 * v1 * w1 * gw*gxf
           s(7,ir) = u1 * v1 * w1 * gw
           xsign = 1.
           ysign = 1.
           if (x(ir) < 0.) xsign = -1.
           if (y(ir) < 0.) ysign = -1.
           v(0,ir) = ux(ir)*gi(ir)*ysign
           v(1,ir) = uy(ir)*gi(ir)*xsign
           v(2,ir) = uz(ir)*gi(ir)
        enddo
      else
c       --- normal loop
        do ir = 1,nptmp
           ix  = int((x(ir) - xmmin) * dxi)
           u1 =      (x(ir) - xmmin) * dxi - ix
           u0 = 1. - u1
           iy  = int((y(ir) - ymmin) * dyi)
           v1 =      (y(ir) - ymmin) * dyi - iy
           v0 = 1. - v1
           iz  = int((z(ir) - zgrid - zmmin) * dzi)
           w1 =      (z(ir) - zgrid - zmmin) * dzi - iz
           w0 = 1. - w1
           ind0 = ix + iy*(nx+1) + iz*(nx+1)*(ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           if (nw == np) cdensw(ix:ix+1) = cdens(ix:ix+1)*wght(ir)
           s(0,ir) = u0 * v0 * w0 * cdensw(ix  )
           s(1,ir) = u1 * v0 * w0 * cdensw(ix+1)
           s(2,ir) = u0 * v1 * w0 * cdensw(ix  )
           s(3,ir) = u1 * v1 * w0 * cdensw(ix+1)
           s(4,ir) = u0 * v0 * w1 * cdensw(ix  )
           s(5,ir) = u1 * v0 * w1 * cdensw(ix+1)
           s(6,ir) = u0 * v1 * w1 * cdensw(ix  )
           s(7,ir) = u1 * v1 * w1 * cdensw(ix+1)
           v(0,ir) = ux(ir)*gi(ir)
           v(1,ir) = uy(ir)*gi(ir)
           v(2,ir) = uz(ir)*gi(ir)
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETJ3D1)
      do ir = 1,nptmp
        if (uz(ir) /= 0.) then
          do m = 0, 7
             ind0 = 3*indx(m,ir)
             j1d(ind0  ) = j1d(ind0  ) + s(m,ir)*v(0,ir)
             j1d(ind0+1) = j1d(ind0+1) + s(m,ir)*v(1,ir)
             j1d(ind0+2) = j1d(ind0+2) + s(m,ir)*v(2,ir)
          enddo
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETJ3D1)

      deallocate(indx,s,v)

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species
      allocate(ii(0:nparpgrp-1), jj(0:nparpgrp-1), kk(0:nparpgrp-1))
      allocate(s(0:7,0:nparpgrp-1))


c     --- vectorized loop to compute indices, weights
      if (l2symtry .and. .not. lcylindrical) then
        do ir = 1,nptmp
          ii(ir) = int((x(ir) - xmmin) * dxi)
          u1     =     (x(ir) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(y(ir)) - ymmin)*dyi)
          v1     =     (abs(y(ir)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((z(ir) - zgrid - zmmin) * dzi)
          w1     =     (z(ir) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gyf = 1.
          if (jj(ir) == 0 .and. ymmin == 0.) gyf = 2.
          if (nw == np) gw = g*wght(ir)
          s(0,ir) = u0 * v0 * w0 * gw*gyf
          s(1,ir) = u1 * v0 * w0 * gw*gyf
          s(2,ir) = u0 * v1 * w0 * gw
          s(3,ir) = u1 * v1 * w0 * gw
          s(4,ir) = u0 * v0 * w1 * gw*gyf
          s(5,ir) = u1 * v0 * w1 * gw*gyf
          s(6,ir) = u0 * v1 * w1 * gw
          s(7,ir) = u1 * v1 * w1 * gw
          ysign = 1.
          if (y(ir) < 0.) ysign = -1.
          v(0,ir) = ux(ir)*gi(ir)*ysign
          v(1,ir) = uy(ir)*gi(ir)
          v(2,ir) = uz(ir)*gi(ir)
        enddo
      elseif (l4symtry .and. .not. lcylindrical) then
        do ir = 1,nptmp
          ii(ir) = int((abs(x(ir)) - xmmin)*dxi)
          u1     =     (abs(x(ir)) - xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(y(ir)) - ymmin)*dyi)
          v1     =     (abs(y(ir)) - ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((z(ir) - zgrid - zmmin) * dzi)
          w1     =     (z(ir) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gxf = 1.
          gyf = 1.
          if (ii(ir) == 0 .and. xmmin == 0.) gxf = 2.
          if (jj(ir) == 0 .and. ymmin == 0.) gyf = 2.
          if (nw == np) gw = g*wght(ir)
          s(0,ir) = u0 * v0 * w0 * gw*gxf*gyf
          s(1,ir) = u1 * v0 * w0 * gw*gyf
          s(2,ir) = u0 * v1 * w0 * gw*gxf
          s(3,ir) = u1 * v1 * w0 * gw
          s(4,ir) = u0 * v0 * w1 * gw*gxf*gyf
          s(5,ir) = u1 * v0 * w1 * gw*gyf
          s(6,ir) = u0 * v1 * w1 * gw*gxf
          s(7,ir) = u1 * v1 * w1 * gw
          xsign = 1.
          ysign = 1.
          if (x(ir) < 0.) xsign = -1.
          if (y(ir) < 0.) ysign = -1.
          v(0,ir) = ux(ir)*gi(ir)*ysign
          v(1,ir) = uy(ir)*gi(ir)*xsign
          v(2,ir) = uz(ir)*gi(ir)
        enddo
      else
c       --- normal loop
        do ir = 1,nptmp
          ii(ir) = int((x(ir) - xmmin) * dxi)
          u1     =     (x(ir) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((y(ir) - ymmin) * dyi)
          v1     =     (y(ir) - ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((z(ir) - zgrid - zmmin) * dzi)
          w1     =     (z(ir) - zgrid - zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          if (nw == np) cdensw(ii(ir):ii(ir)+1)=cdens(ii(ir):ii(ir)+1)*wght(ir)
          s(0,ir) = u0 * v0 * w0 * cdensw(ii(ir)  )
          s(1,ir) = u1 * v0 * w0 * cdensw(ii(ir)+1)
          s(2,ir) = u0 * v1 * w0 * cdensw(ii(ir)  )
          s(3,ir) = u1 * v1 * w0 * cdensw(ii(ir)+1)
          s(4,ir) = u0 * v0 * w1 * cdensw(ii(ir)  )
          s(5,ir) = u1 * v0 * w1 * cdensw(ii(ir)+1)
          s(6,ir) = u0 * v1 * w1 * cdensw(ii(ir)  )
          s(7,ir) = u1 * v1 * w1 * cdensw(ii(ir)+1)
          v(0,ir) = ux(ir)*gi(ir)
          v(1,ir) = uy(ir)*gi(ir)
          v(2,ir) = uz(ir)*gi(ir)
        enddo
      endif
c     --- scalar loop does the actual deposition
!$OMP CRITICAL (CRITICAL_SETJ3D2)
      do ir = 1, nptmp
         if (uz(ir) /= 0) then
         ix = ii(ir)
         iy = jj(ir)
         iz = kk(ir)
         j(0,ix  ,iy  ,iz  ) = j(0,ix  ,iy  ,iz  ) + s(0,ir)*v(0,ir)
         j(0,ix+1,iy  ,iz  ) = j(0,ix+1,iy  ,iz  ) + s(1,ir)*v(0,ir)
         j(0,ix  ,iy+1,iz  ) = j(0,ix  ,iy+1,iz  ) + s(2,ir)*v(0,ir)
         j(0,ix+1,iy+1,iz  ) = j(0,ix+1,iy+1,iz  ) + s(3,ir)*v(0,ir)
         j(0,ix  ,iy  ,iz+1) = j(0,ix  ,iy  ,iz+1) + s(4,ir)*v(0,ir)
         j(0,ix+1,iy  ,iz+1) = j(0,ix+1,iy  ,iz+1) + s(5,ir)*v(0,ir)
         j(0,ix  ,iy+1,iz+1) = j(0,ix  ,iy+1,iz+1) + s(6,ir)*v(0,ir)
         j(0,ix+1,iy+1,iz+1) = j(0,ix+1,iy+1,iz+1) + s(7,ir)*v(0,ir)
         j(1,ix  ,iy  ,iz  ) = j(1,ix  ,iy  ,iz  ) + s(0,ir)*v(1,ir)
         j(1,ix+1,iy  ,iz  ) = j(1,ix+1,iy  ,iz  ) + s(1,ir)*v(1,ir)
         j(1,ix  ,iy+1,iz  ) = j(1,ix  ,iy+1,iz  ) + s(2,ir)*v(1,ir)
         j(1,ix+1,iy+1,iz  ) = j(1,ix+1,iy+1,iz  ) + s(3,ir)*v(1,ir)
         j(1,ix  ,iy  ,iz+1) = j(1,ix  ,iy  ,iz+1) + s(4,ir)*v(1,ir)
         j(1,ix+1,iy  ,iz+1) = j(1,ix+1,iy  ,iz+1) + s(5,ir)*v(1,ir)
         j(1,ix  ,iy+1,iz+1) = j(1,ix  ,iy+1,iz+1) + s(6,ir)*v(1,ir)
         j(1,ix+1,iy+1,iz+1) = j(1,ix+1,iy+1,iz+1) + s(7,ir)*v(1,ir)
         j(2,ix  ,iy  ,iz  ) = j(2,ix  ,iy  ,iz  ) + s(0,ir)*v(2,ir)
         j(2,ix+1,iy  ,iz  ) = j(2,ix+1,iy  ,iz  ) + s(1,ir)*v(2,ir)
         j(2,ix  ,iy+1,iz  ) = j(2,ix  ,iy+1,iz  ) + s(2,ir)*v(2,ir)
         j(2,ix+1,iy+1,iz  ) = j(2,ix+1,iy+1,iz  ) + s(3,ir)*v(2,ir)
         j(2,ix  ,iy  ,iz+1) = j(2,ix  ,iy  ,iz+1) + s(4,ir)*v(2,ir)
         j(2,ix+1,iy  ,iz+1) = j(2,ix+1,iy  ,iz+1) + s(5,ir)*v(2,ir)
         j(2,ix  ,iy+1,iz+1) = j(2,ix  ,iy+1,iz+1) + s(6,ir)*v(2,ir)
         j(2,ix+1,iy+1,iz+1) = j(2,ix+1,iy+1,iz+1) + s(7,ir)*v(2,ir)
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETJ3D2)

      deallocate(ii, jj, kk)
      deallocate(s)

c--------------------------------------
c  Begin direct deposition loop
c--------------------------------------
      elseif (depos == "direct1") then

c  Begin main loop over species

c     --- vectorized loop to compute indices, weights
      if (l2symtry .and. .not. lcylindrical) then
        do ir = 1,nptmp
          ix = int((x(ir) - xmmin) * dxi)
          iy = int((abs(y(ir)) - ymmin)*dyi)
          iz = int((z(ir) - zgrid - zmmin) * dzi)
          u1 =     (x(ir) - xmmin) * dxi - ix
          u0 = 1. - u1
          v1 =     (abs(y(ir)) - ymmin)*dyi - iy
          v0 = 1. - v1
          w1 =     (z(ir) - zgrid - zmmin) * dzi - iz
          w0 = 1. - w1
          gyf = 1.
          if (iy == 0 .and. ymmin == 0.) gyf = 2.
          if (nw == np) gw = g*wght(ir)
          ysign = 1.
          if (y(ir) < 0.) ysign = -1.
          vv(0) = ux(ir)*gi(ir)*ysign
          vv(1) = uy(ir)*gi(ir)
          vv(2) = uz(ir)*gi(ir)
          j(0,ix  ,iy  ,iz  ) = j(0,ix  ,iy  ,iz  ) + u0*v0*w0*gw*vv(0)*gyf
          j(0,ix+1,iy  ,iz  ) = j(0,ix+1,iy  ,iz  ) + u1*v0*w0*gw*vv(0)*gyf
          j(0,ix  ,iy+1,iz  ) = j(0,ix  ,iy+1,iz  ) + u0*v1*w0*gw*vv(0)
          j(0,ix+1,iy+1,iz  ) = j(0,ix+1,iy+1,iz  ) + u1*v1*w0*gw*vv(0)
          j(0,ix  ,iy  ,iz+1) = j(0,ix  ,iy  ,iz+1) + u0*v0*w1*gw*vv(0)*gyf
          j(0,ix+1,iy  ,iz+1) = j(0,ix+1,iy  ,iz+1) + u1*v0*w1*gw*vv(0)*gyf
          j(0,ix  ,iy+1,iz+1) = j(0,ix  ,iy+1,iz+1) + u0*v1*w1*gw*vv(0)
          j(0,ix+1,iy+1,iz+1) = j(0,ix+1,iy+1,iz+1) + u1*v1*w1*gw*vv(0)
          j(1,ix  ,iy  ,iz  ) = j(1,ix  ,iy  ,iz  ) + u0*v0*w0*gw*vv(1)*gyf
          j(1,ix+1,iy  ,iz  ) = j(1,ix+1,iy  ,iz  ) + u1*v0*w0*gw*vv(1)*gyf
          j(1,ix  ,iy+1,iz  ) = j(1,ix  ,iy+1,iz  ) + u0*v1*w0*gw*vv(1)
          j(1,ix+1,iy+1,iz  ) = j(1,ix+1,iy+1,iz  ) + u1*v1*w0*gw*vv(1)
          j(1,ix  ,iy  ,iz+1) = j(1,ix  ,iy  ,iz+1) + u0*v0*w1*gw*vv(1)*gyf
          j(1,ix+1,iy  ,iz+1) = j(1,ix+1,iy  ,iz+1) + u1*v0*w1*gw*vv(1)*gyf
          j(1,ix  ,iy+1,iz+1) = j(1,ix  ,iy+1,iz+1) + u0*v1*w1*gw*vv(1)
          j(1,ix+1,iy+1,iz+1) = j(1,ix+1,iy+1,iz+1) + u1*v1*w1*gw*vv(1)
          j(2,ix  ,iy  ,iz  ) = j(2,ix  ,iy  ,iz  ) + u0*v0*w0*gw*vv(2)*gyf
          j(2,ix+1,iy  ,iz  ) = j(2,ix+1,iy  ,iz  ) + u1*v0*w0*gw*vv(2)*gyf
          j(2,ix  ,iy+1,iz  ) = j(2,ix  ,iy+1,iz  ) + u0*v1*w0*gw*vv(2)
          j(2,ix+1,iy+1,iz  ) = j(2,ix+1,iy+1,iz  ) + u1*v1*w0*gw*vv(2)
          j(2,ix  ,iy  ,iz+1) = j(2,ix  ,iy  ,iz+1) + u0*v0*w1*gw*vv(2)*gyf
          j(2,ix+1,iy  ,iz+1) = j(2,ix+1,iy  ,iz+1) + u1*v0*w1*gw*vv(2)*gyf
          j(2,ix  ,iy+1,iz+1) = j(2,ix  ,iy+1,iz+1) + u0*v1*w1*gw*vv(2)
          j(2,ix+1,iy+1,iz+1) = j(2,ix+1,iy+1,iz+1) + u1*v1*w1*gw*vv(2)
        enddo
      elseif (l4symtry .and. .not. lcylindrical) then
        do ir = 1,nptmp
          ix = int((abs(x(ir)) - xmmin)*dxi)
          iy = int((abs(y(ir)) - ymmin)*dyi)
          iz = int((z(ir) - zgrid - zmmin) * dzi)
          u1 =     (abs(x(ir)) - xmmin)*dxi - ix
          u0 = 1. - u1
          v1 =     (abs(y(ir)) - ymmin)*dyi - iy
          v0 = 1. - v1
          w1     =     (z(ir) - zgrid - zmmin) * dzi - iz
          w0     = 1. - w1
          gxf = 1.
          gyf = 1.
          if (ix == 0 .and. xmmin == 0.) gxf = 2.
          if (iy == 0 .and. ymmin == 0.) gyf = 2.
          if (nw == np) gw = g*wght(ir)
          xsign = 1.
          ysign = 1.
          if (x(ir) < 0.) xsign = -1.
          if (y(ir) < 0.) ysign = -1.
          vv(0) = ux(ir)*gi(ir)*ysign
          vv(1) = uy(ir)*gi(ir)*xsign
          vv(2) = uz(ir)*gi(ir)
          j(0,ix  ,iy  ,iz  ) = j(0,ix  ,iy  ,iz  ) + u0*v0*w0*gw*vv(0)*gxf*gyf
          j(0,ix+1,iy  ,iz  ) = j(0,ix+1,iy  ,iz  ) + u1*v0*w0*gw*vv(0)*gyf
          j(0,ix  ,iy+1,iz  ) = j(0,ix  ,iy+1,iz  ) + u0*v1*w0*gw*vv(0)*gxf
          j(0,ix+1,iy+1,iz  ) = j(0,ix+1,iy+1,iz  ) + u1*v1*w0*gw*vv(0)
          j(0,ix  ,iy  ,iz+1) = j(0,ix  ,iy  ,iz+1) + u0*v0*w1*gw*vv(0)*gxf*gyf
          j(0,ix+1,iy  ,iz+1) = j(0,ix+1,iy  ,iz+1) + u1*v0*w1*gw*vv(0)*gyf
          j(0,ix  ,iy+1,iz+1) = j(0,ix  ,iy+1,iz+1) + u0*v1*w1*gw*vv(0)*gxf
          j(0,ix+1,iy+1,iz+1) = j(0,ix+1,iy+1,iz+1) + u1*v1*w1*gw*vv(0)
          j(1,ix  ,iy  ,iz  ) = j(1,ix  ,iy  ,iz  ) + u0*v0*w0*gw*vv(1)*gxf*gyf
          j(1,ix+1,iy  ,iz  ) = j(1,ix+1,iy  ,iz  ) + u1*v0*w0*gw*vv(1)*gyf
          j(1,ix  ,iy+1,iz  ) = j(1,ix  ,iy+1,iz  ) + u0*v1*w0*gw*vv(1)*gxf
          j(1,ix+1,iy+1,iz  ) = j(1,ix+1,iy+1,iz  ) + u1*v1*w0*gw*vv(1)
          j(1,ix  ,iy  ,iz+1) = j(1,ix  ,iy  ,iz+1) + u0*v0*w1*gw*vv(1)*gxf*gyf
          j(1,ix+1,iy  ,iz+1) = j(1,ix+1,iy  ,iz+1) + u1*v0*w1*gw*vv(1)*gyf
          j(1,ix  ,iy+1,iz+1) = j(1,ix  ,iy+1,iz+1) + u0*v1*w1*gw*vv(1)*gxf
          j(1,ix+1,iy+1,iz+1) = j(1,ix+1,iy+1,iz+1) + u1*v1*w1*gw*vv(1)
          j(2,ix  ,iy  ,iz  ) = j(2,ix  ,iy  ,iz  ) + u0*v0*w0*gw*vv(2)*gxf*gyf
          j(2,ix+1,iy  ,iz  ) = j(2,ix+1,iy  ,iz  ) + u1*v0*w0*gw*vv(2)*gyf
          j(2,ix  ,iy+1,iz  ) = j(2,ix  ,iy+1,iz  ) + u0*v1*w0*gw*vv(2)*gxf
          j(2,ix+1,iy+1,iz  ) = j(2,ix+1,iy+1,iz  ) + u1*v1*w0*gw*vv(2)
          j(2,ix  ,iy  ,iz+1) = j(2,ix  ,iy  ,iz+1) + u0*v0*w1*gw*vv(2)*gxf*gyf
          j(2,ix+1,iy  ,iz+1) = j(2,ix+1,iy  ,iz+1) + u1*v0*w1*gw*vv(2)*gyf
          j(2,ix  ,iy+1,iz+1) = j(2,ix  ,iy+1,iz+1) + u0*v1*w1*gw*vv(2)*gxf
          j(2,ix+1,iy+1,iz+1) = j(2,ix+1,iy+1,iz+1) + u1*v1*w1*gw*vv(2)
        enddo
      else
c       --- normal loop
        do ir = 1,nptmp
          ix = int((x(ir) - xmmin) * dxi)
          iy = int((y(ir) - ymmin) * dyi)
          iz = int((z(ir) - zgrid - zmmin) * dzi)
          u1 =     (x(ir) - xmmin) * dxi - ix
          u0 = 1. - u1
          v1 =     (y(ir) - ymmin) * dyi - iy
          v0 = 1. - v1
          w1 =     (z(ir) - zgrid - zmmin) * dzi - iz
          w0 = 1. - w1
          if (nw == np) cdensw(ix:ix+1)=cdens(ix:ix+1)*wght(ir)
          u0 = u0*cdensw(ix  )
          u1 = u1*cdensw(ix+1)
          vv(0) = ux(ir)*gi(ir)
          vv(1) = uy(ir)*gi(ir)
          vv(2) = uz(ir)*gi(ir)
          j(0,ix  ,iy  ,iz  ) = j(0,ix  ,iy  ,iz  ) + u0*v0*w0*vv(0)
          j(0,ix+1,iy  ,iz  ) = j(0,ix+1,iy  ,iz  ) + u1*v0*w0*vv(0)
          j(0,ix  ,iy+1,iz  ) = j(0,ix  ,iy+1,iz  ) + u0*v1*w0*vv(0)
          j(0,ix+1,iy+1,iz  ) = j(0,ix+1,iy+1,iz  ) + u1*v1*w0*vv(0)
          j(0,ix  ,iy  ,iz+1) = j(0,ix  ,iy  ,iz+1) + u0*v0*w1*vv(0)
          j(0,ix+1,iy  ,iz+1) = j(0,ix+1,iy  ,iz+1) + u1*v0*w1*vv(0)
          j(0,ix  ,iy+1,iz+1) = j(0,ix  ,iy+1,iz+1) + u0*v1*w1*vv(0)
          j(0,ix+1,iy+1,iz+1) = j(0,ix+1,iy+1,iz+1) + u1*v1*w1*vv(0)
          j(1,ix  ,iy  ,iz  ) = j(1,ix  ,iy  ,iz  ) + u0*v0*w0*vv(1)
          j(1,ix+1,iy  ,iz  ) = j(1,ix+1,iy  ,iz  ) + u1*v0*w0*vv(1)
          j(1,ix  ,iy+1,iz  ) = j(1,ix  ,iy+1,iz  ) + u0*v1*w0*vv(1)
          j(1,ix+1,iy+1,iz  ) = j(1,ix+1,iy+1,iz  ) + u1*v1*w0*vv(1)
          j(1,ix  ,iy  ,iz+1) = j(1,ix  ,iy  ,iz+1) + u0*v0*w1*vv(1)
          j(1,ix+1,iy  ,iz+1) = j(1,ix+1,iy  ,iz+1) + u1*v0*w1*vv(1)
          j(1,ix  ,iy+1,iz+1) = j(1,ix  ,iy+1,iz+1) + u0*v1*w1*vv(1)
          j(1,ix+1,iy+1,iz+1) = j(1,ix+1,iy+1,iz+1) + u1*v1*w1*vv(1)
          j(2,ix  ,iy  ,iz  ) = j(2,ix  ,iy  ,iz  ) + u0*v0*w0*vv(2)
          j(2,ix+1,iy  ,iz  ) = j(2,ix+1,iy  ,iz  ) + u1*v0*w0*vv(2)
          j(2,ix  ,iy+1,iz  ) = j(2,ix  ,iy+1,iz  ) + u0*v1*w0*vv(2)
          j(2,ix+1,iy+1,iz  ) = j(2,ix+1,iy+1,iz  ) + u1*v1*w0*vv(2)
          j(2,ix  ,iy  ,iz+1) = j(2,ix  ,iy  ,iz+1) + u0*v0*w1*vv(2)
          j(2,ix+1,iy  ,iz+1) = j(2,ix+1,iy  ,iz+1) + u1*v0*w1*vv(2)
          j(2,ix  ,iy+1,iz+1) = j(2,ix  ,iy+1,iz+1) + u0*v1*w1*vv(2)
          j(2,ix+1,iy+1,iz+1) = j(2,ix+1,iy+1,iz+1) + u1*v1*w1*vv(2)
        enddo
      endif

      endif

      enddo
!$OMP END DO

      if (lcylindrical) then
        deallocate(x,y)
        deallocate(ux,uy)
      endif

!$OMP END PARALLEL

      if (lf3dtimesubs) timesetj3d = timesetj3d + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine getjforfieldsolve()
      use InGen
      use InGen3d
      use BFieldGrid

#ifndef MPIPARALLEL

      if (solvergeom==XYZgeom) then
c       --- If jp is not associated with j, then copy the data.
        if (.not. associated(bfield%j,bfieldp%j)) then
          if (bfieldp%nx == bfield%nx .and.
     &        bfieldp%ny == bfield%ny .and.
     &        bfieldp%nzlocal == bfield%nzlocal) then
            bfield%j = bfieldp%j
          else
            call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            call remark("ERROR!! jp and j are not the same shape!   ")
            call remark("     The current density will not be properly used!")
            call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          endif
        endif
      endif

#else

      use Parallel
      use Fields3dParticles,Only: nzpguard

c  For parallel version, each processor sends j to neighboring processors
c  whose field solve region overlap its particle region.
      if(solvergeom==XZgeom) then
c       print*,"Error: Self magnetic field not supported with solvergeom=XZgeom"
c       call getjforfieldsolverz(nx,nzlocal,j)
c     else if(solvergeom==Zgeom) then
c       print*,"Error: Self magnetic field not supported with solvergeom=Zgeom"
c       call getjforfieldsolvez(nzlocal,j)
      elseif (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        call setjforfieldsolve3d(bfield%nx,bfield%ny,bfield%nzlocal,bfield%j,
     &                           bfieldp%nx,bfieldp%ny,bfieldp%nzlocal,bfieldp%j,
     &                           nzpguard,
     &                           my_index,nslaves,izpslave,nzpslave,
     &                           izfsslave,nzfsslave)
      end if

#endif

      return
      end
c=============================================================================
      subroutine setjforfieldsolve3d(nx,ny,nzlocal,j,nxp,nyp,nzp,jp,nzpguard,
     &                               my_index,nslaves,izpslave,nzpslave,
     &                               izfsslave,nzfsslave)

      use BFieldGrid,Only: ljsummedonboundaries
      integer(ISZ):: nx,ny,nzlocal,nxp,nyp,nzp,nzpguard
      real(kind=8):: j(0:2,0:nx,0:ny,0:nzlocal)
      real(kind=8):: jp(0:2,0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

#ifdef MPIPARALLEL

      if (.not. ljsummedonboundaries) then
        call sumsourcepondomainboundaries(jp,3,nxp,nyp,nzp,nzpguard,
     &                                    my_index,nslaves,izpslave,nzpslave)
      endif
      call setsourceforfieldsolve3d_parallel(3,nx,ny,nzlocal,j,nxp,nyp,nzp,jp,
     &                                       my_index,nslaves,izpslave,nzpslave,
     &                                       izfsslave,nzfsslave)

#endif

      return
      end
c=============================================================================
      subroutine setupbfieldsforparticles3d(ns,ndts,it,bfield,bfieldp)
      use Picglb3d
      use BFieldGridTypemodule
#ifdef MPIPARALLEL
      use Parallel
#endif

      integer(ISZ):: ns,ndts(0:ns-1),it
      type(BFieldGridType):: bfield,bfieldp

c Ensures that the bp and jp arrays are setup properly.
c Also, setup jptmp for the case when species have different time step sizes.
c The setup of jptmp is skipped is ns == -1.

      integer(ISZ):: js,i,oldnsndtsj
      integer(ISZ),allocatable:: lijtmp(:)

      bfieldp%nx = bfield%nx
      bfieldp%ny = bfield%ny
      bfieldp%nz = bfield%nz
#ifndef MPIPARALLEL
      bfieldp%nzlocal = bfield%nzlocal
#else
      bfieldp%nzlocal = nzpslave(my_index)
#endif
      bfieldp%dx = bfield%dx
      bfieldp%dy = bfield%dy
      bfieldp%dz = bfield%dz
      bfieldp%xmmin = bfield%xmmin
      bfieldp%xmmax = bfield%xmmax
      bfieldp%ymmin = bfield%ymmin
      bfieldp%ymmax = bfield%ymmax
      bfieldp%zmminlocal = bfield%zmminlocal
      bfieldp%zmmaxlocal = bfield%zmmaxlocal
      bfieldp%zmmin = bfield%zmmin
      bfieldp%zmmax = bfield%zmmax
      bfieldp%bounds = bfield%bounds
      bfieldp%lcylindrical = bfield%lcylindrical

c     --- If either jp of bp are not associated or either does not have
c     --- the correct dimensions as given by nxj, etc, something 
c     --- needs to be done for fix them.
      if (.not. associated(bfieldp%j) .or. .not. associated(bfieldp%b) .or.
     &  ANY(ubound(bfieldp%b).ne.(/2,bfieldp%nx,bfieldp%ny,bfieldp%nzlocal/)) .or.
     &  ANY(ubound(bfieldp%j).ne.(/2,bfieldp%nx,bfieldp%ny,bfieldp%nzlocal/))) then

c       --- Free whatever the jp and bp arrays are now refering to.
        call BFieldGridTypefree(bfieldp)

        if (ALL( (/ bfieldp%nx,bfieldp%ny,bfieldp%nzlocal /) ==
     &           (/ bfield%nx, bfield%ny, bfield%nz /) )) then
c         --- If the sizes of jp and j, and bp and b,
c         --- are the same,
c         --- then associate jp with j and bp with b
          bfieldp%j => bfield%j
          bfieldp%b => bfield%b
        else
c         --- If the sizes are not the same, then allocate the group
          call BFieldGridTypeallot(bfieldp)
        endif
      endif

c     --- If -1 was passed in for ns, then skip the setup of jptmp
      if (ns == -1) return

c     --- Setup jptmp

c     --- First, count how many species take larger time steps and increase
c     --- the size of the arrays appropriately.
      oldnsndtsj = bfieldp%nsndtsj
      bfieldp%nsndtsj = 0
      do js = 0, ns-1
        if (ndts(js) > 1) bfieldp%nsndtsj = bfieldp%nsndtsj + 1
      enddo

      if (bfieldp%nsndtsj > 0) then

c       --- Set local copy of the number of species, used to size jsjptmp
        bfieldp%nsjtmp = ns

c       --- Only change size if it has increased. If it is decreased, do the
c       --- gchange afterward so that data that may still be needed at the end
c       --- of the array is not thrown away.
        if (bfieldp%nsndtsj > oldnsndtsj) then
          call BFieldGridTypechange(bfieldp)
        endif
c       --- Now set the species index into the jptmp array's last dimension.
c       --- Clear out any that have have ndts set back to 1 and find which
c       --- indices are already taken. If ndts is changed and a value is
c       --- changed from or to 1, the species that already have ndts > 1
c       --- must keep the existing data in jptmp. This code fills in any
c       --- empty spaces with new species or ones at the end of the list.
c       --- The lijtmp array flags which indices have already been claimed.
c       --- Also, jsjptmp defaults to -1 for species which have ndts == 1.
        allocate(lijtmp(0:bfieldp%nsndtsj-1))
        lijtmp = 0
        do js = 0, ns-1
          if (ndts(js) == 1) bfieldp%jsjtmp(js) = -1
          if (ndts(js) > 0 .and. bfieldp%jsjtmp(js) >= 0)
     &       lijtmp(bfieldp%jsjtmp(js)) = 1
        enddo

        i = 0
        do js = 0, ns-1
          if (ndts(js) > 1) then
c           --- If ndts is now set, or if this species index falls beyond the
c           --- new size of the array, find a spot for it.
            if (bfieldp%jsjtmp(js) == -1 .or. bfieldp%jsjtmp(js) >= bfieldp%nsndtsj) then
c             --- Find the next available index
              do while (lijtmp(i) == 1)
                i = i + 1
              enddo
c             --- If the species was beyond the end of the array, copy the
c             --- data into the new place, unless it will be recalculated
c             --- anyway this step.
              if (bfieldp%jsjtmp(js) >= bfieldp%nsndtsj .and.
     &            mod(it+1,ndts(js)) /= 0) then
                bfieldp%jtmp(:,:,:,:,i)=bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(js))
              endif
c             --- Finally, do the assignment
              bfieldp%jsjtmp(js) = i
              lijtmp(i) = 1
            endif
          endif
        enddo
        deallocate(lijtmp)
      endif

c     --- If the new size is smaller, the change can now be done since the
c     --- data that needed to be saved is copied to its new place.
c     --- Note that this is done outside the if-block above to catch the case
c     --- when nsndtsj==0 and oldnsndtsj>0.
      if (bfieldp%nsndtsj < oldnsndtsj) then
        call BFieldGridTypechange(bfieldp)
      endif

      return
      end
c=============================================================================
      subroutine fetchb3dfrompositions(jsid,ndts,n,x,y,z,bx,by,bz)
      use GlobalVars
      use Subtimersf3d
      use Picglb
      use Picglb3d
      use InGen
      use InGen3d
      use BFieldGrid
      use FieldSolveAPI
      use f3d_bfield_interfaces
      integer(ISZ):: jsid,ndts,n
      real(kind=8),target:: x(n),y(n),z(n),bx(n),by(n),bz(n)

c     --- Obtain the magnetic field

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if ((fstype == 12 .and. bfstype < 0) .or. bfstype == 12) then
#if defined(XLF) || defined(XLF2)
c       --- The pointer assignment breaks with the xlf compilers, so instead
c       --- new arrays are allocated and those are used as the work space.
c       --- This is less efficient since there is overhead in the allocation
c       --- and the data needs to be copied around.
        allocate(xfsapi(n),yfsapi(n),zfsapi(n))
        allocate(bxfsapi(n),byfsapi(n),bzfsapi(n))
        xfsapi = x
        yfsapi = y
        zfsapi = z
        bxfsapi = 0.
        byfsapi = 0.
        bzfsapi = 0.
#else
        xfsapi => x
        yfsapi => y
        zfsapi => z
        bxfsapi => bx
        byfsapi => by
        bzfsapi => bz
#endif
        npfsapi = n
        jsfsapi = jsid
        ndtsfsapi = ndts
        if (bfstype == 12) then
          call execuser("bfetchbregistered")
        else
          call execuser("fetchbregistered")
        endif
        npfsapi = 0
        jsfsapi = -1
#if defined(XLF) || defined(XLF2)
        bx(1:n) = bxfsapi
        by(1:n) = byfsapi
        bz(1:n) = bzfsapi
        deallocate(xfsapi,yfsapi,zfsapi)
        deallocate(bxfsapi,byfsapi,bzfsapi)
#endif
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(bxfsapi)
        nullify(byfsapi)
        nullify(bzfsapi)
      elseif (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        if (bfstype >= 0) then
          call setb3d(bfieldp%b,n,x,y,z,zgridprv,bx,by,bz,
     &                bfield%nx,bfield%ny,bfield%nzlocal,
     &                bfield%dx,bfield%dy,bfield%dz,
     &                bfield%xmmin,bfield%ymmin,bfield%zmminlocal,
     &                l2symtry,l4symtry,bfield%lcylindrical)
        endif
      else
        bx = 0.
        by = 0.
        bz = 0.
      endif

!$OMP MASTER
      if (lf3dtimesubs) timefetchb3dfrompositions = timefetchb3dfrompositions + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fetcha(n,x,y,z,a)
      use BFieldGrid
      use InGen,Only: fstype,bfstype
      use InGen3d
      use Picglb
      use FieldSolveAPI
      use f3d_bfield_interfaces
      integer(ISZ):: n
      real(kind=8),target:: x(n),y(n),z(n),a(0:2,n)

c     --- Call the appropriate routine to get the A
      if ((fstype == 12 .and. bfstype < 0) .or. bfstype == 12) then
        xfsapi => x
        yfsapi => y
        zfsapi => z
        afsapi => a
        npfsapi = n
        jsfsapi = -1
        ndtsfsapi = 1
        if (bfstype == 12) then
          call execuser("bfetcharegistered")
        else
          call execuser("fetcharegistered")
        endif
        npfsapi = 0
        jsfsapi = -1
        ndtsfsapi = 0
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(afsapi)
      else if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        call fetchafrompositions3d(bfield%a,n,x,y,z,zgrid,a,
     &                             bfield%nx,bfield%ny,bfield%nzlocal,
     &                             bfield%dx,bfield%dy,bfield%dz,
     &                             bfield%xmmin,bfield%ymmin,bfield%zmminlocal,
     &                             l2symtry,l4symtry,bfield%lcylindrical)
c     elseif(solvergeom==XZgeom) then
c       call setarz(n,x,y,z,p,zgrid)
c     elseif(solvergeom==Zgeom) then
c       call setaz(n,z,p,zgrid)
c     elseif(solvergeom==AMRgeom) then
c       call cho_geta3d(n,x,y,z,0.,p,-1,-1)
      endif

      return
      end
c=============================================================================
      subroutine getbforparticles()
      use InGen3d
      use InMesh3d
      use BFieldGrid

c     --- Ensure that the Bp array is setup properly
      if(solvergeom==XYZgeom) then
c       --- Note the -1 is passed in for ns as a flag to skip the
c       --- setting up of jptmp.
        call setupbfieldsforparticles3d(-1,0,0,bfield,bfieldp)
      endif

#ifndef MPIPARALLEL

      if(solvergeom==XYZgeom .or.
     &   solvergeom==RZgeom) then
c       --- If Bp is not associated with B, then copy the data.
        if (.not. associated(bfield%b,bfieldp%b)) then
          if (bfieldp%nx == bfield%nx .and.
     &        bfieldp%ny == bfield%ny .and.
     &        bfieldp%nzlocal == bfield%nzlocal) then
            bfieldp%b = bfield%b
          else
            call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            call remark("ERROR!! bp and b are not the same shape!        ")
            call remark("        The self field will not be properly applied!")
            call remark("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          endif
        endif
      endif

#else

c     --- Distribute B among the processes so each has B in its
c     --- particle domain.
      if(solvergeom==XZgeom .or. solvergeom==Zgeom) then
c       call getbforparticlesrz()
c      else if(solvergeom==Zgeom) then
c        call getbforparticlesz()
      elseif (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        call getbforparticles3d(bfield,bfieldp)
      end if
#endif

      return
      end
c===========================================================================
      subroutine bfieldsol3d(iwhich)
      use GlobalVars
      use Subtimersf3d
      use InGen
      use LatticeInternal
      use Fields3d
      use Timers
      use Picglb, only: zgrid
      use GridBoundary3d
      use BFieldGrid
      integer(ISZ):: iwhich

c  Field solver for 3d warped Cartesian geometry. 
c  Enter with charge density in J array, old potential in B array.
c  Exit with new potential in B array, and J unchanged.
c  
c  For field solve in a bend there are two possible cases:
c      1) Call to an SOR field solver.  
c      2) Iteration loop for a perturbative calculation of fields. 
c         The loop is in bendfieldsol3d

      real(kind=8):: timetemp
      real(kind=8):: substarttime,wtime
      integer(ISZ):: izmin,izmax,tnz
      if (lf3dtimesubs) substarttime = wtime()
      timetemp = wtime()

c     --- If no field solve, return 
      if ( bfstype <= -1) return 

c     --- Calculate rstar here since it is only used by the field solver and
c     --- to ensure that it is set when the field solver is called.
      if (linbend) then
        call setrstar(bfield%rstar(-1),bfield%nzlocal,bfield%dz,bfield%zmminlocal,zgrid)
      endif

      if (bfstype == 7 .or. bfstype == 12) then
c       --- multigrid field solvers - deal with bends directly and set
c       --- axial boundary conditions
        call bvp3d(iwhich,bfstype)

      else if (.not. linbend) then
c       --- If not in bend, call VP3D
        call bvp3d(iwhich,bfstype)

      else
c       --- Call perturbative bent beam field solver
        call bendfieldsol3d

      endif

#ifdef MPIPARALLEL
c     --- Set the flag signalling that j has been parallel summed
c     --- on the boundaries. This is done here only after j has been
c     --- setup for all subcycling and selfb groups.
      ljsummedonboundaries = .true.
#endif

c     --- Get B from A, either taking the finite difference curl of A to get B,
c     --- or just copy it when B is solved for directly.
      call getbfroma3d(bfield%a,bfield%b,bfield%nx,bfield%ny,bfield%nzlocal,
     &                 bfield%dx,bfield%dy,bfield%dz,bfield%xmmin,
     &                 bfield%lcylindrical,bfield%lusevectorpotential)

c     --- If using the analytic form of Btheta, calculate it here.
      if (bfield%lanalyticbtheta) then
        call getanalyticbtheta(bfield%b,bfield%j,
     &                         bfield%nx,bfield%ny,bfield%nzlocal,
     &                         bfield%dx,bfield%xmmin)
      endif

c     --- Distribute B among the processes so each has B in its
c     --- particle domain.
      call getbforparticles()

!$OMP MASTER
      if (lf3dtimesubs) timebfieldsol3d = timebfieldsol3d + wtime() - substarttime
!$OMP END MASTER
      fstime = fstime + (wtime() - timetemp)
      return
      end
c=============================================================================
      subroutine bvp3d(iwhich,bfstype)
      use Constant
      use BFieldGrid
      use InGen3d, only: filt,solvergeom,RZgeom
      integer(ISZ):: iwhich,bfstype

      integer(ISZ):: id,idmax
      real(kind=8),pointer:: atemp(:,:,:),jtemp(:,:,:)
      real(kind=8),pointer:: j(:,:,:,:),curlj(:,:,:,:)

      if (bfstype < 0) return

      if (iwhich==0 .or. iwhich==1) call execuser('initbfieldsolver')

      if (bfstype == 12) then
        if (iwhich <= 0) then
          call execuser('bfieldsolregistered')
          return
        endif
      endif

      if (solvergeom == RZgeom) then
        bfstype = 10
        call init_bworkgrid(bfield%nx,bfield%nzlocal,bfield%dx,bfield%dz,
     &                      bfield%xmmin,bfield%zmminlocal,bfield%bounds,
     &                      bfield%a,bfield%j)
      endif

      if (bfield%lusevectorpotential) then
        j => bfield%j
      else
        allocate(curlj(0:2,0:bfield%nx,0:bfield%ny,0:bfield%nzlocal))
        curlj = 0.
        call curl3d(bfield%j,curlj,bfield%nx,bfield%ny,bfield%nzlocal,
     &              bfield%dx,bfield%dy,bfield%dz,
     &              bfield%xmmin,bfield%lcylindrical,
     &              0,0,0,0)
        j => curlj
      endif

      allocate(atemp(0:bfield%nx,0:bfield%ny,-1:bfield%nzlocal+1))

      if (bfstype <= 4) then
c       --- For the FFT field solvers, only one temp array is needed since
c       --- j is not passed and the field solve happens in place.
c       --- jtemp is pointed to atemp just to give it a non-null value.
        jtemp => atemp
      else
        allocate(jtemp(0:bfield%nx,0:bfield%ny,0:bfield%nzlocal))
      endif

      j = j*mu0*eps0

c     --- Note that the arrays being passed in are not contiguous, which means
c     --- that copies are being done.
c     --- If only initialization is being done (iwhich==1) then the bvp3d_work
c     --- routine only needs to be called once. Proper arrays are still passed
c     --- though they should never be needed during initialization.
      idmax = 2
      if (iwhich == 1) idmax = 0
      do id=0,idmax
        if (bfield%lanalyticbtheta .and.
     &      ((bfield%lusevectorpotential .and. (id == 0 .or. id == 2)) .or.
     &      (.not. bfield%lusevectorpotential .and. id == 1))) cycle
        if (bfstype <= 4) then
          atemp(:,:,0:bfield%nzlocal) = j(id,:,:,:)
        else
          atemp = bfield%a(id,0:bfield%nx,0:bfield%ny,:)
          jtemp = j(id,:,:,:)
        endif
        call bvp3d_work(iwhich,bfstype,atemp,jtemp,id,bfield)
        if (bfstype .ne. 10) bfield%a(id,0:bfield%nx,0:bfield%ny,:) = atemp
      enddo

c     --- Rescale, replacing 1/eps0 with mu0
c     bfield%a = bfield%a*mu0*eps0

#ifdef MPIPARALLEL
      if (bfstype <= 4) then
        call getaforfields3d(bfield%a,bfield%nx,bfield%ny,bfield%nzlocal)
      endif
#endif

c     --- The last argument, llongitudinal, is only true for the FFT based
c     --- solvers since those do not set the z guard cells.
      call setaboundaries3d(bfield%a,bfield%nx,bfield%ny,bfield%nzlocal,
     &                      bfield%zmminlocal,bfield%zmmaxlocal,
     &                      bfield%zmmin,bfield%zmmax,
     &                      bfield%bounds,bfield%lcylindrical,
     &                      bfstype <= 4)

      deallocate(atemp)
      if (bfstype > 4) deallocate(jtemp)

      if (bfield%lusevectorpotential) then
c       --- Unscale the current density
        j = j/(mu0*eps0)
      else
c       --- The unscaling is not needed since the data is thrown away.
        deallocate(curlj)
      endif

      return
      end
c=============================================================================
      subroutine bvp3d_work(iwhich,bfstype,a,j,iaxis,bfield)
      use Subtimersf3d
      use BFieldGridTypemodule
      use LatticeInternal, only: linbend
      use Picglb, only: zbeam,zgrid
      use InGen3d, only: filt,solvergeom,l2symtry,l4symtry
      use GridBoundary3d, only: bound0,boundnz,boundxy
      use PSOR3d, only: gridmode
      use Parallel, only: my_index,nslaves,izfsslave,nzfsslave
      integer(ISZ):: iwhich,bfstype,iaxis
      type(BFieldGridType):: bfield
      real(kind=8):: a(0:bfield%nx,0:bfield%ny,-1:bfield%nzlocal+1)
      real(kind=8):: j(0:bfield%nx,0:bfield%ny,0:bfield%nzlocal)

c Interface to field solver for magnetic vector potential
c Note that bfield is passed in to have access to the data and work arrays
c kxsq, attx, scrtch etc.

      integer(ISZ):: nx,ny,nzlocal,nz
      real(kind=8):: dx,dy,dz
      real(kind=8):: xlen,ylen,zlen
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c Note that the select case statement was giving the f90 compiler on HPUX
c fits and so was replaced with if statements.

      nx = bfield%nx
      ny = bfield%ny
      nzlocal = bfield%nzlocal
      nz = bfield%nz
      dx = bfield%dx
      dy = bfield%dy
      dz = bfield%dz
      xlen = bfield%xmmax - bfield%xmmin
      ylen = bfield%ymmax - bfield%ymmin
      zlen = bfield%zmmax - bfield%zmmin

      if (bfstype == 0) then
        call vpois3d(iwhich,a(0,0,0),a(0,0,0),
     &               bfield%kxsq,bfield%kysq,bfield%kzsq,
     &               bfield%attx,bfield%atty,bfield%attz,
     &               filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &               bfield%scrtch,bfield%xywork,bfield%zwork,
     &               0,l2symtry,l4symtry,bound0,boundnz,boundxy)
      elseif (bfstype == 4) then
        if (iwhich == 1 .or. iwhich == 0) then
          call vpois3d(1,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
        endif
        if (iwhich == -1 .or. iwhich == 0) then
          call vpois3d(12,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
          call vpois3d(14,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
          call vpois3d(13,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nzlocal,nz,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
        endif

      elseif (bfstype == 7) then

        call multigrid3dsolve(iwhich,nx,ny,nzlocal,nz,dx,dy,dz,a,j,
     &                      bfield%rstar,linbend,bfield%bounds,
     &                      bfield%xmmin,bfield%ymmin,bfield%zmminlocal,
     &                      bfield%zmmin,zbeam,zgrid,
     &                      bfield%mgparam(iaxis),bfield%mgform(iaxis),
     &                      bfield%mgiters(iaxis),bfield%mgmaxiters(iaxis),
     &                      bfield%mgmaxlevels(iaxis),bfield%mgerror(iaxis),
     &                      bfield%mgtol(iaxis),bfield%mgverbose(iaxis),
     &                      bfield%downpasses(iaxis),bfield%uppasses(iaxis),
     &                      bfield%lcndbndy,bfield%laddconductor,
     &                      bfield%icndbndy,.false.,
     &                      gridmode,bfield%conductors,
     &                      my_index,nslaves,izfsslave,nzfsslave)


      elseif (bfstype == 10) then
        call multigridrzb(iwhich,iaxis,bfield%a(iaxis,:,0,:),
     &                    j(:,0,:),bfield%nx,bfield%nzlocal,bfield%mgtol(iaxis))

      endif

!$OMP MASTER
      if (lf3dtimesubs) timebvp3d = timebvp3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
