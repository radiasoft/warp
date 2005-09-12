#include "top.h"
c=============================================================================
c@(#) File F3D_BFIELD.F, version $Revision: 1.2 $, $Date: 2005/09/12 07:22:53 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This contains routines to handle the B field, calculated from the current
c  density.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      module f3d_bfield_interfaces
      interface

      subroutine setj3d(bfield,j1d,np,xp,yp,zp,zgrid,uxp,uyp,uzp,gaminv,
     &                  q,wght,depos,l2symtry,l4symtry)
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield
      integer(ISZ):: np
      real(kind=8):: zgrid,q,wght
      real(kind=8):: j1d(0:3*(1+bfield%nx)*(1+bfield%ny)*(1+bfield%nz)-1)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8),target:: uxp(np), uyp(np), uzp(np), gaminv(np)
      character(8):: depos
      logical(ISZ):: l2symtry,l4symtry
      end subroutine setj3d

      subroutine setb3d(bfield,np,xp,yp,zp,zgrid,bx,by,bz,l2symtry,l4symtry)
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield
      integer(ISZ):: np
      real(kind=8):: zgrid
      real(kind=8),target:: xp(np),yp(np),zp(np)
      real(kind=8):: bx(np),by(np),bz(np)
      logical(ISZ):: l2symtry,l4symtry
      end subroutine setb3d

      subroutine fetchafrompositions3d(np,xp,yp,zp,zgrid,bfield,
     &                                 l2symtry,l4symtry)
      use BFieldGridTypemodule
      integer(ISZ):: np
      real(kind=8):: zgrid
      real(kind=8),target:: xp(np),yp(np),zp(np)
      real(kind=8):: a(0:2,np)
      type(BFieldGridType):: bfield
      logical(ISZ):: l2symtry,l4symtry
      end subroutine fetchafrompositions3d

      end interface
      end module f3d_bfield_interfaces
c=============================================================================
      subroutine init_bfieldsolver(bfstype)
      use GlobalVars
      use BFieldGrid
      use InMesh3d
      use InGen3d, only: l2symtry,l4symtry,solvergeom,RZgeom
      use Fields3d, only: nzfull
      use GridBoundary3d
      use Multigrid3d

      integer(ISZ):: bfstype

      if (bfield%xmmin == 0) bfield%xmmin = xmmin
      if (bfield%xmmax == 0) bfield%xmmax = xmmax
      if (bfield%ymmin == 0) bfield%ymmin = ymmin
      if (bfield%ymmax == 0) bfield%ymmax = ymmax
      if (bfield%zmmin == 0) bfield%zmmin = zmmin
      if (bfield%zmmax == 0) bfield%zmmax = zmmax
      if (bfield%zmminglobal == 0) bfield%zmminglobal = zmminglobal
      if (bfield%zmmaxglobal == 0) bfield%zmmaxglobal = zmmaxglobal
      if (bfield%nx == 0) bfield%nx = nx
      if (bfield%ny == 0) bfield%ny = ny
      if (bfield%nz == 0) bfield%nz = nz
      if (bfield%nzfull == 0) bfield%nzfull = nzfull
      if (bfield%dx == 0) bfield%dx = (bfield%xmmax-bfield%xmmin)/bfield%nx
      if (bfield%dy == 0 .and. bfield%ny > 0) bfield%dy = (bfield%ymmax-bfield%ymmin)/bfield%ny
      if (bfield%dz == 0) bfield%dz = (bfield%zmmax-bfield%zmmin)/bfield%nz
      call BFieldGridTypechange(bfield)

      bfield%icndbndy = icndbndy
      bfield%mgparam = mgparam
      bfield%mgmaxiters = mgmaxiters
      bfield%mgmaxlevels = mgmaxlevels
      bfield%mgtol = mgtol
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

      call init_bworkgrid(bfield%nx,bfield%nz,bfield%dx,bfield%dz,
     &                    bfield%xmmin,bfield%zmmin,bfield%bounds,
     &                    bfield%a,bfield%j)
      call bvp3d(1,bfstype)

      return
      end
c=============================================================================
      subroutine loadj3d(ins_i,nps_i,is_i,lzero)
      use GlobalVars
      use Subtimersf3d
      use InGen
      use InGen3d
      use InPart
      use Picglb
      use Picglb3d
      use Particles
      use BFieldGrid
      use GridBoundary3d
      use f3d_bfield_interfaces
      integer(ISZ):: ins_i,nps_i,is_i
      logical(ISZ):: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- jp array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.

      integer(ISZ):: ins_u,nps_u
      integer(ISZ):: is1,is2
      integer(ISZ):: ip,ipmin,is
      real(kind=8):: swtmp
      real(kind=8):: substarttime,wtime
      real(kind=8),dimension(:,:,:,:),pointer::jpcopy
      integer(ISZ):: allocerror

      if (lf3dtimesubs) substarttime = wtime()

c     --- Ensure that the jp array is setup properly
      if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        call setupbfieldsforparticles3d(ns,ndts,it,bfield,bfieldp)
      endif

c     --- zero jp if requested
      if (lzero) then
        if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
          call zeroarry(bfieldp%j,
     &                  3*(bfieldp%nx+1)*(bfieldp%ny+1)*(bfieldp%nz+1))
        endif
      end if

c     --- set limits on loop over species
      if (is_i == -1) then
        is1 = 1
        is2 = ns
      else
        is1 = is_i
        is2 = is_i
      endif

      if(depos /= 'none') then

c     --- set initial limits from input
c     --- (will be changed if necessary in the loop)
      ins_u = ins_i
      nps_u = nps_i

c     --- loop over species
      do is=is1,is2

c       --- For the AMR/Chombo version, all of the lost particles must be
c       --- cleared out.
c       if (solvergeom == AMRgeom) call clearpart(is,1)

c       --- get loop limits for particles if needed
        if (ins_i == -1) ins_u = ins(is)
        if (nps_i == -1) nps_u = nps(is)

c       --- Scale the weight, sw, by the time step scale size. This only
c       --- makes sense for steady-state and slice modes. In time-dependent
c       --- mode, it is assumed that dtscale has not been changed from 1.
        swtmp = sw(is)*dtscale(is)

        if (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
          ipmin = ins_u
          if (nps_u > 0) then
            if (ndts(is)>1) then
c             --- This code would be better if it were in setj3d, but since
c             --- jp is passed into it, this must be done here.
              if (mod(it+1,ndts(is))==0) then
c                --- If this species is being advanced this step, then
c                --- deposit its j into jptmp.
                 jpcopy => bfieldp%j
                 bfieldp%j => bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
                 bfieldp%j = 0.
              else
c                --- If this species is not being advanced this step, then
c                --- just add its saved jptmp into jp.
                 bfieldp%j = bfieldp%j + bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
                 cycle
               end if
            end if
            call setj3d(bfieldp,bfieldp%j,nps_u,
     &                  xp(ipmin:ipmin+nps_u-1),yp(ipmin:ipmin+nps_u-1),
     &                  zp(ipmin:ipmin+nps_u-1),zgrid,
     &                  uxp(ipmin:ipmin+nps_u-1),uyp(ipmin:ipmin+nps_u-1),
     &                  uzp(ipmin:ipmin+nps_u-1),gaminv(ipmin:ipmin+nps_u-1),
     &                  sq(is),swtmp,depos,l2symtry,l4symtry)
            if (ndts(is)>1 .and. mod(it+1,ndts(is))==0) then
c             --- If this species is being advanced this step, then restore
c             --- jp and copy its j into jp.
              bfieldp%j => jpcopy
              bfieldp%j = bfieldp%j + bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
            end if
          endif
        elseif(solvergeom==XZgeom) then
          print*,"Error: Self magnetic field not support with solvergeom=RZgeom or XZgeom"
        elseif(solvergeom==Zgeom) then
          print*,"Error: Self magnetic field not support with solvergeom=Zgeom"
        elseif(solvergeom==Rgeom) then
          print*,"Error: Self magnetic field not support with solvergeom=Rgeom"
        elseif(solvergeom==AMRgeom) then
          print*,"Error: Self magnetic field not support with solvergeom=AMRgeom"
        end if
      enddo

      if (solvergeom == XYZgeomMR) then
        call execuser("loadjMR")
      endif

      endif

c  For parallel version, each processor sends j to neighboring processors
c  whose field solve region overlap its particle region.
      if (lzero) call getjforfieldsolve()

c     --- enforce periodicity
c     --- This is done on the j array, rather than jp, since j will
c     --- be used for the field solve. Also, in some cases, the jp array
c     --- will not cover the full axial extent of the system, so 
c     --- periodicity can only be enforced on j which does cover the full
c     --- system.
      if (lzero) then
        call perj3d(bfield,bound0,boundxy)
      endif

!$OMP MASTER
      if (lf3dtimesubs) timeloadj3d = timeloadj3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setaboundaries3d(bfield)
      use Subtimersf3d
      use GlobalVars
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield

c  Sets the slices on the exterior of A, depending on the boundary conditions.

      real(kind=8),pointer,dimension(:,:,:,:):: a
      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c       call perarz()
c       return
c     end if

      a => bfield%a

c     --- With the cylindrical solver, the guard cells are already set.
      if (.not. bfield%lcylindrical) then
        if (bfield%bounds(0) == periodic .or. bfield%bounds(1) == periodic) then
          a(:,-1,:,:)   = a(:,bfield%nx-1,:,:)
          a(:,bfield%nx,:,:)   = a(:,0,:,:)
          a(:,bfield%nx+1,:,:) = a(:,1,:,:)
        else
          if (bfield%bounds(0) == neumann) then
            a(:,-1,:,:) = a(:,1,:,:)
          elseif (bfield%bounds(0) == dirichlet) then
            a(:,-1,:,:) = 2.*a(:,0,:,:) - a(:,1,:,:)
          endif
          if (bfield%bounds(1) == neumann) then
            a(:,bfield%nx+1,:,:) = a(:,bfield%nx-1,:,:)
          elseif (bfield%bounds(1) == dirichlet) then
            a(:,bfield%nx+1,:,:) = 2.*a(:,bfield%nx,:,:) - a(:,bfield%nx-1,:,:)
          endif
        endif
      endif

      if (bfield%lcylindrical) then
c       --- In this case, actually ny == 0. This ensures that d/dtheta is zero.
        a(:,:,-1,:) = a(:,:,0,:)
        a(:,:,bfield%ny+1,:) = a(:,:,bfield%ny,:)
      else
        if (bfield%bounds(2) == periodic .or. bfield%bounds(3) == periodic) then
          a(:,:,-1,:)   = a(:,:,bfield%ny-1,:)
          a(:,:,bfield%ny,:)   = a(:,:,0,:)
          a(:,:,bfield%ny+1,:) = a(:,:,1,:)
        else
          if (bfield%bounds(2) == neumann) then
            a(:,:,-1,:) = a(:,:,1,:)
          elseif (bfield%bounds(2) == dirichlet) then
            a(:,:,-1,:) = 2.*a(:,:,0,:) - a(:,:,1,:)
          endif
          if (bfield%bounds(3) == neumann) then
            a(:,:,bfield%ny+1,:) = a(:,:,bfield%ny-1,:)
          elseif (bfield%bounds(3) == dirichlet) then
            a(:,:,bfield%ny+1,:) = 2.*a(:,:,bfield%ny,:) - a(:,:,bfield%ny-1,:)
          endif
        endif
      endif

c     --- With the cylindrical solver, the guard cells are already set.
      if (.not. bfield%lcylindrical) then
        if (bfield%bounds(4) == periodic .or. bfield%bounds(5) == periodic) then
#ifdef MPIPARALLEL
          call getaforfields3d(bfield)
          call pera3d_slave(bfield)
#else
          a(:,:,:,-1)   = a(:,:,:,bfield%nz-1)
          a(:,:,:,bfield%nz)   = a(:,:,:,0)
          a(:,:,:,bfield%nz+1) = a(:,:,:,1)
#endif
        else
          if (bfield%bounds(4) == neumann) then
            a(:,:,:,-1) = a(:,:,:,1)
          elseif (bfield%bounds(4) == dirichlet) then
            a(:,:,:,-1) = 2.*a(:,:,:,0) - a(:,:,:,1)
          endif
          if (bfield%bounds(5) == neumann) then
            a(:,:,:,bfield%nz+1) = a(:,:,:,bfield%nz-1)
          elseif (bfield%bounds(5) == dirichlet) then
            a(:,:,:,bfield%nz+1) = 2.*a(:,:,:,bfield%nz) - a(:,:,:,bfield%nz-1)
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

c  Distribute j for 2d solver
c     if(solvergeom==RZgeom .or. solvergeom==XZgeom .or.
c    &   solvergeom==XYgeom .or. solvergeom==Zgeom) then
c       call distribute_j_rz()
c     end if

c  Enforce longitudinal periodicity.
c     if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==XYgeom) then
c       call jbndrz()
cifdef MPIPARALLEL
c       if(bound0==periodic) call perjrz()
cendif
c     end if

c  Copy charge density from frz.basegrid to f3d.j
c     if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c       call get_j_rz(j,nx,nz,1,0)
c     else if(solvergeom==XYgeom) then
c       call get_j_rz(j,nx,ny,1,0)
c     else if(solvergeom==Zgeom) then
c       call get_j_z(j,nz,1,0)
c     else if(solvergeom==Rgeom) then
c       call get_j_r(j,nx,1,0)
c     end if

      if (solvergeom==XYZgeom .and. bound0==periodic) then
#ifdef MPIPARALLEL
        call perj3d_slave(bfield)
#else
        j(:,:,:,0)  = j(:,:,:,0) + j(:,:,:,bfield%nz)
        j(:,:,:,bfield%nz) = j(:,:,:,0)
#endif
      endif

!$OMP MASTER
      if (lf3dtimesubs) timeperj3d = timeperj3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setb3d(bfield,np,xp,yp,zp,zgrid,bx,by,bz,l2symtry,l4symtry)
      use Subtimersf3d
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield
      integer(ISZ):: np
      real(kind=8):: zgrid
      real(kind=8),target:: xp(np),yp(np),zp(np)
      real(kind=8):: bx(np),by(np),bz(np)
      logical(ISZ):: l2symtry,l4symtry

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

      if (bfield%lcylindrical) then
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
      tdxi = 1./(2.*bfield%dx)
      tdyi = 1./(2.*bfield%dy)
      tdzi = 1./(2.*bfield%dz)
      dxi = 1./bfield%dx
      dyi = 1./bfield%dy
      dzi = 1./bfield%dz

      if (.not. (l2symtry .or. l4symtry)) then
        do ip = 1, np

          i = (x(ip) - bfield%xmmin) * dxi
          j = (y(ip) - bfield%ymmin) * dyi
          k = (z(ip) - zgrid - bfield%zmmin) * dzi

          u1 = (x(ip) - bfield%xmmin) * dxi - i
          v1 = (y(ip) - bfield%ymmin) * dyi - j
          w1 = (z(ip) - zgrid - bfield%zmmin) * dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

          bx(ip) = u0*v0*w0*bfield%b(0,i  ,j  ,k  )
     &           + u1*v0*w0*bfield%b(0,i+1,j  ,k  )
     &           + u0*v1*w0*bfield%b(0,i  ,j+1,k  )
     &           + u1*v1*w0*bfield%b(0,i+1,j+1,k  )
     &           + u0*v0*w1*bfield%b(0,i  ,j  ,k+1)
     &           + u1*v0*w1*bfield%b(0,i+1,j  ,k+1)
     &           + u0*v1*w1*bfield%b(0,i  ,j+1,k+1)
     &           + u1*v1*w1*bfield%b(0,i+1,j+1,k+1)

          by(ip) = u0*v0*w0*bfield%b(1,i  ,j  ,k  )
     &           + u1*v0*w0*bfield%b(1,i+1,j  ,k  )
     &           + u0*v1*w0*bfield%b(1,i  ,j+1,k  )
     &           + u1*v1*w0*bfield%b(1,i+1,j+1,k  )
     &           + u0*v0*w1*bfield%b(1,i  ,j  ,k+1)
     &           + u1*v0*w1*bfield%b(1,i+1,j  ,k+1)
     &           + u0*v1*w1*bfield%b(1,i  ,j+1,k+1)
     &           + u1*v1*w1*bfield%b(1,i+1,j+1,k+1)

          bz(ip) = u0*v0*w0*bfield%b(2,i  ,j  ,k  )
     &           + u1*v0*w0*bfield%b(2,i+1,j  ,k  )
     &           + u0*v1*w0*bfield%b(2,i  ,j+1,k  )
     &           + u1*v1*w0*bfield%b(2,i+1,j+1,k  )
     &           + u0*v0*w1*bfield%b(2,i  ,j  ,k+1)
     &           + u1*v0*w1*bfield%b(2,i+1,j  ,k+1)
     &           + u0*v1*w1*bfield%b(2,i  ,j+1,k+1)
     &           + u1*v1*w1*bfield%b(2,i+1,j+1,k+1)

        enddo

      else

c       --- Set the signs of the B field for particles on negative side of
c       --- the axis of symmetry.
        sy = -1.
        sx = 1.
        if (l4symtry) sx = -1.

c       --- special loop symmetry is used
        do ip = 1, np

          i = (abs(x(ip)) - bfield%xmmin)*dxi
          j = (abs(y(ip)) - bfield%ymmin)*dyi
          k = (z(ip) - zgrid - bfield%zmmin)*dzi

          u1 = (abs(x(ip)) - bfield%xmmin)*dxi - i
          v1 = (abs(y(ip)) - bfield%ymmin)*dyi - j
          w1 = (z(ip) - zgrid - bfield%zmmin)*dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Adjust sign of B field for appropriate quadrant.
          xsign = +1.
          ysign = +1.
          if (x(ip) < 0.) xsign = sx
          if (y(ip) < 0.) ysign = sy

          bx(ip) = xsign*(u0*v0*w0*bfield%b(0,i  ,j  ,k  )
     &                  + u1*v0*w0*bfield%b(0,i+1,j  ,k  )
     &                  + u0*v1*w0*bfield%b(0,i  ,j+1,k  )
     &                  + u1*v1*w0*bfield%b(0,i+1,j+1,k  )
     &                  + u0*v0*w1*bfield%b(0,i  ,j  ,k+1)
     &                  + u1*v0*w1*bfield%b(0,i+1,j  ,k+1)
     &                  + u0*v1*w1*bfield%b(0,i  ,j+1,k+1)
     &                  + u1*v1*w1*bfield%b(0,i+1,j+1,k+1))

          by(ip) = ysign*(u0*v0*w0*bfield%b(1,i  ,j  ,k  )
     &                  + u1*v0*w0*bfield%b(1,i+1,j  ,k  )
     &                  + u0*v1*w0*bfield%b(1,i  ,j+1,k  )
     &                  + u1*v1*w0*bfield%b(1,i+1,j+1,k  )
     &                  + u0*v0*w1*bfield%b(1,i  ,j  ,k+1)
     &                  + u1*v0*w1*bfield%b(1,i+1,j  ,k+1)
     &                  + u0*v1*w1*bfield%b(1,i  ,j+1,k+1)
     &                  + u1*v1*w1*bfield%b(1,i+1,j+1,k+1))

          bz(ip) =        u0*v0*w0*bfield%b(2,i  ,j  ,k  )
     &                  + u1*v0*w0*bfield%b(2,i+1,j  ,k  )
     &                  + u0*v1*w0*bfield%b(2,i  ,j+1,k  )
     &                  + u1*v1*w0*bfield%b(2,i+1,j+1,k  )
     &                  + u0*v0*w1*bfield%b(2,i  ,j  ,k+1)
     &                  + u1*v0*w1*bfield%b(2,i+1,j  ,k+1)
     &                  + u0*v1*w1*bfield%b(2,i  ,j+1,k+1)
     &                  + u1*v1*w1*bfield%b(2,i+1,j+1,k+1)

        enddo

      endif

      if (bfield%lcylindrical) then
c       --- Transform Br and Btheta into Bx and By
        y = atan2(yp,xp)
        bxtemp = bx
        bx = bx*cos(y) + by*sin(y)
        by = by*cos(y) - bxtemp*sin(y)
        deallocate(x,y,bxtemp)
      endif

!$OMP MASTER
      if (lf3dtimesubs) timesetb3d = timesetb3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fetchafrompositions3d(np,xp,yp,zp,zgrid,bfield,
     &                                 l2symtry,l4symtry)
      use Subtimersf3d
      use BFieldGridTypemodule
      integer(ISZ):: np
      real(kind=8):: zgrid
      real(kind=8),target:: xp(np),yp(np),zp(np)
      real(kind=8):: a(0:2,np)
      type(BFieldGridType):: bfield
      logical(ISZ):: l2symtry,l4symtry

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

      if (bfield%lcylindrical) then
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
      tdxi = 1./(2.*bfield%dx)
      tdyi = 1./(2.*bfield%dy)
      tdzi = 1./(2.*bfield%dz)
      dxi = 1./bfield%dx
      dyi = 1./bfield%dy
      dzi = 1./bfield%dz

      if (.not. (l2symtry .or. l4symtry)) then
        do ip = 1, np

          i = (x(ip) - bfield%xmmin) * dxi
          j = (y(ip) - bfield%ymmin) * dyi
          k = (z(ip) - zgrid - bfield%zmmin) * dzi

          u1 = (x(ip) - bfield%xmmin) * dxi - i
          v1 = (y(ip) - bfield%ymmin) * dyi - j
          w1 = (z(ip) - zgrid - bfield%zmmin) * dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

          a(0,ip) = u0*v0*w0*bfield%a(0,i  ,j  ,k  )
     &            + u1*v0*w0*bfield%a(0,i+1,j  ,k  )
     &            + u0*v1*w0*bfield%a(0,i  ,j+1,k  )
     &            + u1*v1*w0*bfield%a(0,i+1,j+1,k  )
     &            + u0*v0*w1*bfield%a(0,i  ,j  ,k+1)
     &            + u1*v0*w1*bfield%a(0,i+1,j  ,k+1)
     &            + u0*v1*w1*bfield%a(0,i  ,j+1,k+1)
     &            + u1*v1*w1*bfield%a(0,i+1,j+1,k+1)

          a(1,ip) = u0*v0*w0*bfield%a(1,i  ,j  ,k  )
     &            + u1*v0*w0*bfield%a(1,i+1,j  ,k  )
     &            + u0*v1*w0*bfield%a(1,i  ,j+1,k  )
     &            + u1*v1*w0*bfield%a(1,i+1,j+1,k  )
     &            + u0*v0*w1*bfield%a(1,i  ,j  ,k+1)
     &            + u1*v0*w1*bfield%a(1,i+1,j  ,k+1)
     &            + u0*v1*w1*bfield%a(1,i  ,j+1,k+1)
     &            + u1*v1*w1*bfield%a(1,i+1,j+1,k+1)

          a(2,ip) = u0*v0*w0*bfield%a(2,i  ,j  ,k  )
     &            + u1*v0*w0*bfield%a(2,i+1,j  ,k  )
     &            + u0*v1*w0*bfield%a(2,i  ,j+1,k  )
     &            + u1*v1*w0*bfield%a(2,i+1,j+1,k  )
     &            + u0*v0*w1*bfield%a(2,i  ,j  ,k+1)
     &            + u1*v0*w1*bfield%a(2,i+1,j  ,k+1)
     &            + u0*v1*w1*bfield%a(2,i  ,j+1,k+1)
     &            + u1*v1*w1*bfield%a(2,i+1,j+1,k+1)

        enddo

      else

c       --- Set the signs of the B field for particles on negative side of
c       --- the axis of symmetry.
        sy = -1.
        sx = 1.
        if (l4symtry) sx = -1.

c       --- special loop symmetry is used
        do ip = 1, np

          i = (abs(x(ip)) - bfield%xmmin)*dxi
          j = (abs(y(ip)) - bfield%ymmin)*dyi
          k = (z(ip) - zgrid - bfield%zmmin)*dzi

          u1 = (abs(x(ip)) - bfield%xmmin)*dxi - i
          v1 = (abs(y(ip)) - bfield%ymmin)*dyi - j
          w1 = (z(ip) - zgrid - bfield%zmmin)*dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Adjust sign of B field for appropriate quadrant.
          xsign = +1.
          ysign = +1.
          if (x(ip) < 0.) xsign = sx
          if (y(ip) < 0.) ysign = sy

          a(0,ip) = xsign*(u0*v0*w0*bfield%a(0,i  ,j  ,k  )
     &                   + u1*v0*w0*bfield%a(0,i+1,j  ,k  )
     &                   + u0*v1*w0*bfield%a(0,i  ,j+1,k  )
     &                   + u1*v1*w0*bfield%a(0,i+1,j+1,k  )
     &                   + u0*v0*w1*bfield%a(0,i  ,j  ,k+1)
     &                   + u1*v0*w1*bfield%a(0,i+1,j  ,k+1)
     &                   + u0*v1*w1*bfield%a(0,i  ,j+1,k+1)
     &                   + u1*v1*w1*bfield%a(0,i+1,j+1,k+1))

          a(1,ip) = ysign*(u0*v0*w0*bfield%a(1,i  ,j  ,k  )
     &                   + u1*v0*w0*bfield%a(1,i+1,j  ,k  )
     &                   + u0*v1*w0*bfield%a(1,i  ,j+1,k  )
     &                   + u1*v1*w0*bfield%a(1,i+1,j+1,k  )
     &                   + u0*v0*w1*bfield%a(1,i  ,j  ,k+1)
     &                   + u1*v0*w1*bfield%a(1,i+1,j  ,k+1)
     &                   + u0*v1*w1*bfield%a(1,i  ,j+1,k+1)
     &                   + u1*v1*w1*bfield%a(1,i+1,j+1,k+1))

          a(2,ip) =        u0*v0*w0*bfield%a(2,i  ,j  ,k  )
     &                   + u1*v0*w0*bfield%a(2,i+1,j  ,k  )
     &                   + u0*v1*w0*bfield%a(2,i  ,j+1,k  )
     &                   + u1*v1*w0*bfield%a(2,i+1,j+1,k  )
     &                   + u0*v0*w1*bfield%a(2,i  ,j  ,k+1)
     &                   + u1*v0*w1*bfield%a(2,i+1,j  ,k+1)
     &                   + u0*v1*w1*bfield%a(2,i  ,j+1,k+1)
     &                   + u1*v1*w1*bfield%a(2,i+1,j+1,k+1)

        enddo

      endif

      if (bfield%lcylindrical) then
c       --- Transform Br and Btheta into Bx and By
        y = atan2(yp,xp)
        axtemp = a(0,:)
        a(0,:) = a(0,:)*cos(y) + a(1,:)*sin(y)
        a(1,:) = a(1,:)*cos(y) - axtemp*sin(y)
        deallocate(x,y,axtemp)
      endif

!$OMP MASTER
      if (lf3dtimesubs) timefetchafrompositions3d = timefetchafrompositions3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=========================================================================== 
      subroutine getbfroma3d(bfield)
      use GlobalVars
      use Subtimersf3d
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield

c Calculate the B via finite differences of A.

      integer(ISZ):: ix,iy,iz
      real(kind=8):: tdxi,tdyi,tdzi,r,ri
      real(kind=8),pointer,dimension(:,:,:,:):: a,b
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      tdxi = 0.5/bfield%dx
      tdyi = 0.5/bfield%dy
      tdzi = 0.5/bfield%dz

      a => bfield%a
      b => bfield%b

c     --- Do the calculation. In this case, guard cells have been added to
c     --- A since having special cases for each B component on each
c     --- boundary for each boundary condition begins to get very complicated.
      do iz=0,bfield%nz
        do iy=0,bfield%ny
          do ix=0,bfield%nx
            b(0,ix,iy,iz) = + tdyi*(a(2,ix  ,iy+1,iz  ) - a(2,ix  ,iy-1,iz  ))
     &                      - tdzi*(a(1,ix  ,iy  ,iz+1) - a(1,ix  ,iy  ,iz-1))
            b(1,ix,iy,iz) = + tdzi*(a(0,ix  ,iy  ,iz+1) - a(0,ix  ,iy  ,iz-1))
     &                      - tdxi*(a(2,ix+1,iy  ,iz  ) - a(2,ix-1,iy  ,iz  ))
            b(2,ix,iy,iz) = + tdxi*(a(1,ix+1,iy  ,iz  ) - a(1,ix-1,iy  ,iz  ))
     &                      - tdyi*(a(0,ix  ,iy+1,iz  ) - a(0,ix  ,iy-1,iz  ))
          enddo
        enddo
      enddo

c     --- When lcylindrical is true, an extra term is needed for Bz
c     ---      1 d               d Atheta   Atheta
c     --- Bz = - -- (r Atheta) = -------- + ------
c     ---      r dr                 dr        r
c     --- The derivative is calculated above. Note that the d/dtheta terms
c     --- which appear in Br and Bz are still ignored, assuming axisymmetry.
c     --- Because of the way the guard cells are set for the y plane,
c     --- the derivatives there calculated above are zero.
      if (bfield%lcylindrical) then
        do ix=0,bfield%nx
          r = (bfield%xmmin + bfield%dx*ix)
          if (r > 0.) then
            ri = 1./r
            b(2,ix,:,:) = b(2,ix,:,:) + a(1,ix,:,:)*ri
          else
c           --- This can be justified by writing out the del cross A in
c           --- Cartesian coordinates, and using
c           --- a(1,1,0,iz) = a(0,0,-1,iz) = -a(1,-1,0,iz) = -a(0,0,+1,iz)
            b(2,ix,:,:) = 2.*b(2,ix,:,:)
          endif
        enddo
      endif

!$OMP MASTER
      if (lf3dtimesubs) timegetbfroma3d = timegetbfroma3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setj3d(bfield,j1d,np,xp,yp,zp,zgrid,uxp,uyp,uzp,gaminv,
     &                  q,wght,depos,l2symtry,l4symtry)
      use GlobalVars
      use Constant
      use Subtimersf3d
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield
      integer(ISZ):: np
      real(kind=8):: zgrid,q,wght
      real(kind=8):: j1d(0:3*(1+bfield%nx)*(1+bfield%ny)*(1+bfield%nz)-1)
      real(kind=8),target:: xp(np), yp(np), zp(np)
      real(kind=8),target:: uxp(np), uyp(np), uzp(np), gaminv(np)
      character(8):: depos
      logical(ISZ):: l2symtry,l4symtry

c  Sets current density

c  Algorithm notes: j array is dimensioned (3,0:nx,0:ny,0:nz) outside,
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
      integer(ISZ):: moff(0:7), indx(0:7,nparpgrp)
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(nparpgrp), jj(nparpgrp), kk(nparpgrp)
c     --- For both algorithms
      real(kind=8):: s(0:7,nparpgrp),v(0:2,nparpgrp)
c     --- Work array holding q*w/cell volume, the charge density per
c     --- simulation particle. This is primarily needed for the RZ version
c     --- since the cell volume there has radial dependence.
      real(kind=8):: cdens(0:bfield%nx)
c     --- Temp arrays to hold particle data
c     --- These are needed when lcylindrical is true, in which case x=r, y=0,
c     --- and vx = vr, and vy = vtheta.
      real(kind=8),pointer:: x(:),y(:),z(:),ux(:),uy(:),uz(:),gi(:)

      real(kind=8),pointer,dimension(:,:,:,:):: j
      integer(ISZ):: ipmin,nptmp,ix,iy,iz,ind0,m,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: xsign,ysign
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c  Set up offset array for vectorized deposition:

      moff(0) = 0
      moff(1) = 1
      moff(2) = bfield%nx+1
      moff(3) = bfield%nx+2
      moff(4) = (bfield%nx+1)*(bfield%ny+1)
      moff(5) = (bfield%nx+1)*(bfield%ny+1)+1
      moff(6) = (bfield%nx+1)*(bfield%ny+1)+bfield%nx+1
      moff(7) = (bfield%nx+1)*(bfield%ny+1)+bfield%nx+2

      dxi = 1./bfield%dx
      dyi = 1./bfield%dy
      dzi = 1./bfield%dz

      if (bfield%lcylindrical) then
        if (bfield%xmmin == 0.) then
c         --- The factor 0.75 corrects for overdeposition due to linear
c         --- weighting (for uniform distribution)
c         --- see Larson et al., Comp. Phys. Comm., 90:260-266, 1995
c         --- and Verboncoeur, J. of Comp. Phys.,
          cdens(0) = 0.75/(pi*(0.5*0.5*bfield%dx*bfield%dx*bfield%dz))
        else
          cdens(0) = 1./(2.*pi*(bfield%xmmin)*bfield%dx*bfield%dz)
        endif
        do ix = 1,bfield%nx
          cdens(ix) = 1./(2.*pi*(ix*bfield%dx+bfield%xmmin)*bfield%dx*bfield%dz)
        enddo
        cdens = cdens*wght*q
      else
        g = wght * q / (bfield%dx * bfield%dy * bfield%dz)
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

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ipmin,nptmp,i,j,k,u1,u0,v1,v0,w1,w0,ir,ip,ind0,indx,
!$OMP&                 gyfact,gxfact,s,m,ii,jj,kk,x,y,z,ux,uy,uz)
!$OMP&FIRSTPRIVATE(np)

      if (bfield%lcylindrical) then
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
        if (bfield%lcylindrical) then
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

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
c       --- special loop for 2-fold symmetry
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        do ir = 1,nptmp
           ix  = int((x(ir) - bfield%xmmin) * dxi)
           u1 =      (x(ir) - bfield%xmmin) * dxi - ix
           u0 = 1. - u1
           iy  = int((abs(y(ir)) - bfield%ymmin)*dyi)
           v1 =      (abs(y(ir)) - bfield%ymmin)*dyi - iy
           v0 = 1. - v1
           iz  = int((z(ir) - zgrid - bfield%zmmin) * dzi)
           w1 =      (z(ir) - zgrid - bfield%zmmin) * dzi - iz
           w0 = 1. - w1
           ind0 = ix + iy*(bfield%nx+1) + iz*(bfield%nx+1)*(bfield%ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gyfact = 1.
           if (iy == 0 .and. bfield%ymmin == 0.) gyfact = 2.
           g = cdens(ix)
           s(0,ir) = u0 * v0 * w0 * g*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
           ysign = 1.
           if (y(ir) < 0.) ysign = -1.
           v(0,ir) = ux(ir)*gi(ir)*ysign
           v(1,ir) = uy(ir)*gi(ir)
           v(2,ir) = uz(ir)*gi(ir)
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        do ir = 1,nptmp
           ix  = int((abs(x(ir)) - bfield%xmmin)*dxi)
           u1 =      (abs(x(ir)) - bfield%xmmin)*dxi - ix
           u0 = 1. - u1
           iy  = int((abs(y(ir)) - bfield%ymmin)*dyi)
           v1 =      (abs(y(ir)) - bfield%ymmin)*dyi - iy
           v0 = 1. - v1
           iz  = int((z(ir) - zgrid - bfield%zmmin) * dzi)
           w1 =      (z(ir) - zgrid - bfield%zmmin) * dzi - iz
           w0 = 1. - w1
           ind0 = ix + iy*(bfield%nx+1) + iz*(bfield%nx+1)*(bfield%ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           gxfact = 1.
           gyfact = 1.
           if (ix == 0 .and. bfield%xmmin == 0.) gxfact = 2.
           if (iy == 0 .and. bfield%ymmin == 0.) gyfact = 2.
           g = cdens(ix)
           s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g*gxfact
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g*gxfact
           s(7,ir) = u1 * v1 * w1 * g
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
           ix  = int((x(ir) - bfield%xmmin) * dxi)
           u1 =      (x(ir) - bfield%xmmin) * dxi - ix
           u0 = 1. - u1
           iy  = int((y(ir) - bfield%ymmin) * dyi)
           v1 =      (y(ir) - bfield%ymmin) * dyi - iy
           v0 = 1. - v1
           iz  = int((z(ir) - zgrid - bfield%zmmin) * dzi)
           w1 =      (z(ir) - zgrid - bfield%zmmin) * dzi - iz
           w0 = 1. - w1
           ind0 = ix + iy*(bfield%nx+1) + iz*(bfield%nx+1)*(bfield%ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           g = cdens(ix)
           s(0,ir) = u0 * v0 * w0 * g
           s(1,ir) = u1 * v0 * w0 * g
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g
           s(5,ir) = u1 * v0 * w1 * g
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
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

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        do ir = 1,nptmp
          ii(ir) = int((x(ir) - bfield%xmmin) * dxi)
          u1     =     (x(ir) - bfield%xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(y(ir)) - bfield%ymmin)*dyi)
          v1     =     (abs(y(ir)) - bfield%ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((z(ir) - zgrid - bfield%zmmin) * dzi)
          w1     =     (z(ir) - zgrid - bfield%zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gyfact = 1.
          if (jj(ir) == 0 .and. bfield%ymmin == 0.) gyfact = 2.
          g = cdens(ii(ir))
          s(0,ir) = u0 * v0 * w0 * g*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
          ysign = 1.
          if (y(ir) < 0.) ysign = -1.
          v(0,ir) = ux(ir)*gi(ir)*ysign
          v(1,ir) = uy(ir)*gi(ir)
          v(2,ir) = uz(ir)*gi(ir)
        enddo
      elseif (l4symtry) then
        do ir = 1,nptmp
          ii(ir) = int((abs(x(ir)) - bfield%xmmin)*dxi)
          u1     =     (abs(x(ir)) - bfield%xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(y(ir)) - bfield%ymmin)*dyi)
          v1     =     (abs(y(ir)) - bfield%ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((z(ir) - zgrid - bfield%zmmin) * dzi)
          w1     =     (z(ir) - zgrid - bfield%zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0 .and. bfield%xmmin == 0.) gxfact = 2.
          if (jj(ir) == 0 .and. bfield%ymmin == 0.) gyfact = 2.
          g = cdens(ii(ir))
          s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g*gxfact
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g*gxfact
          s(7,ir) = u1 * v1 * w1 * g
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
          ii(ir) = int((x(ir) - bfield%xmmin) * dxi)
          u1     =     (x(ir) - bfield%xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((y(ir) - bfield%ymmin) * dyi)
          v1     =     (y(ir) - bfield%ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((z(ir) - zgrid - bfield%zmmin) * dzi)
          w1     =     (z(ir) - zgrid - bfield%zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          g = cdens(ii(ir))
          s(0,ir) = u0 * v0 * w0 * g
          s(1,ir) = u1 * v0 * w0 * g
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g
          s(5,ir) = u1 * v0 * w1 * g
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
          v(0,ir) = ux(ir)*gi(ir)
          v(1,ir) = uy(ir)*gi(ir)
          v(2,ir) = uz(ir)*gi(ir)
        enddo
      endif
c     --- scalar loop does the actual deposition
      j => bfield%j
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

      endif

      enddo
!$OMP END DO

      if (bfield%lcylindrical) then
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
     &        bfieldp%nz == bfield%nz) then
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

c  For parallel version, each processor sends j to neighboring processors
c  whose field solve region overlap its particle region.
      if(solvergeom==XZgeom) then
c       print*,"Error: Self magnetic field not supported with solvergeom=RZgeom or XZgeom"
c       call getjforfieldsolverz(nx,nz,j)
c     else if(solvergeom==Zgeom) then
c       print*,"Error: Self magnetic field not supported with solvergeom=Zgeom"
c       call getjforfieldsolvez(nz,j)
      elseif (solvergeom==XYZgeom .or. solvergeom==RZgeom) then
        call getjforfieldsolve3d(bfield,bfieldp)
      end if

#endif

      return
      end
c=============================================================================
      subroutine getjforfieldsolve3d(bfield,bfieldp)
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield,bfieldp

#ifdef MPIPARALLEL

      call sumjondomainboundaries(bfieldp)
      call getjforfieldsolve3d_parallel(bfield,bfieldp)

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
      bfieldp%nzfull = bfield%nzfull
#ifndef MPIPARALLEL
      bfieldp%nz = bfield%nz
#else
      bfieldp%nz = nzpslave(my_index)
#endif
      bfieldp%dx = bfield%dx
      bfieldp%dy = bfield%dy
      bfieldp%dz = bfield%dz
      bfieldp%xmmin = bfield%xmmin
      bfieldp%xmmax = bfield%xmmax
      bfieldp%ymmin = bfield%ymmin
      bfieldp%ymmax = bfield%ymmax
      bfieldp%zmmin = bfield%zmmin
      bfieldp%zmmax = bfield%zmmax
      bfieldp%zmminglobal = bfield%zmminglobal
      bfieldp%zmmaxglobal = bfield%zmmaxglobal
      bfieldp%bounds = bfield%bounds
      bfieldp%lcylindrical = bfield%lcylindrical

c     --- If either jp of bp are not associated or either does not have
c     --- the correct dimensions as given by nxj, etc, something 
c     --- needs to be done for fix them.
      if (.not. associated(bfieldp%j) .or. .not. associated(bfieldp%b) .or.
     &  ANY(ubound(bfieldp%b).ne.(/3,bfieldp%nx,bfieldp%ny,bfieldp%nz+1/)) .or.
     &  ANY(ubound(bfieldp%j).ne.(/3,bfieldp%nx,bfieldp%ny,bfieldp%nz  /))) then

c       --- Free whatever the jp and bp arrays are now refering to.
        call BFieldGridTypefree(bfieldp)

        if (ALL( (/ bfieldp%nx,bfieldp%ny,bfieldp%nz /) ==
     &           (/ bfield%nx, bfield%ny, bfield%nzfull /) )) then
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
      subroutine fetchb3dfrompositions(is,n,x,y,z,bx,by,bz)
      use GlobalVars
      use Subtimersf3d
      use Picglb
      use Picglb3d
      use InGen
      use InGen3d
      use BFieldGrid
      use FieldSolveAPI
      use f3d_bfield_interfaces
      integer(ISZ):: is,n
      real(kind=8),target:: x(n),y(n),z(n),bx(n),by(n),bz(n)

c     --- Obtain the magnetic field

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (.not. lcallfetchb) then
        bx = 0.
        by = 0.
        bz = 0.
        return
      endif

      if(solvergeom==XYZgeom .and. bfstype >= 0) then
        call setb3d(bfieldp,n,x,y,z,zgrid,bx,by,bz,l2symtry,l4symtry)
      elseif(solvergeom==XYZgeomMR) then
        isfsapi = is
        xfsapi => x
        yfsapi => y
        zfsapi => z
        bxfsapi => bx
        byfsapi => by
        bzfsapi => bz
        call execuser("fetchbMR")
        isfsapi = -1
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(bxfsapi)
        nullify(byfsapi)
        nullify(bzfsapi)
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
      use InGen3d
      use Picglb
      use FieldSolveAPI
      use f3d_bfield_interfaces
      integer(ISZ):: n
      real(kind=8),target:: x(n),y(n),z(n),a(0:2,n)

c     --- Call the appropriate routine to get the A
      if(solvergeom==XYZgeom) then
        call fetchafrompositions3d(n,x,y,z,zgrid,bfield,l2symtry,l4symtry)
c     elseif(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c       call setarz(n,x,y,z,p,zgrid)
c     elseif(solvergeom==Zgeom) then
c       call setaz(n,z,p,zgrid)
c     elseif(solvergeom==AMRgeom) then
c       call cho_geta3d(n,x,y,z,0.,p,-1,-1)
      elseif(solvergeom==XYZgeomMR) then
        xfsapi => x
        yfsapi => y
        zfsapi => z
        afsapi => a
        call execuser("fetchaMR")
        nullify(xfsapi)
        nullify(yfsapi)
        nullify(zfsapi)
        nullify(afsapi)
      endif

      return
      end
c=============================================================================
      subroutine getbforparticles()
      use InGen3d
      use InMesh3d
      use BFieldGrid

c     --- Ensure that the Bp array is setup properly
      if(solvergeom==XYZgeom .or. solvergeom==XYZgeomMR) then
c       --- Note the -1 is passed in for ns as a flag to skip the
c       --- setting up of jptmp.
        call setupbfieldsforparticles3d(-1,0,0,bfield,bfieldp)
      endif

#ifndef MPIPARALLEL

      if(solvergeom==XYZgeom .or. solvergeom==XYZgeomMR .or.
     &   solvergeom==RZgeom) then
c       --- If Bp is not associated with B, then copy the data.
        if (.not. associated(bfield%b,bfieldp%b)) then
          if (bfieldp%nx == bfield%nx .and.
     &        bfieldp%ny == bfield%ny .and.
     &        bfieldp%nz == bfield%nz) then
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
        call setrstar(bfield%rstar(-1),bfield%nz,bfield%dz,bfield%zmmin,zgrid)
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

c     --- Take the finite differences of A to get B.
      call getbfroma3d(bfield)

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
      integer(ISZ):: iwhich,bfstype

      integer(ISZ):: id,idmax
      real(kind=8),pointer:: atemp(:,:,:),jtemp(:,:,:)

      allocate(atemp(0:bfield%nx,0:bfield%ny,-1:bfield%nz+1))

      if (bfstype <= 4) then
c       --- For the FFT field solvers, only one temp array is needed since
c       --- rho is not passed and the field solve happens in place.
c       --- jtemp is pointed to atemp just to give it a non-null value.
        jtemp => atemp
      else
        allocate(jtemp(0:bfield%nx,0:bfield%ny,0:bfield%nz))
      endif

      bfield%j = bfield%j*mu0*eps0

c     --- Note that the arrays being passed in are not contiguous, which means
c     --- that copies are being done.
c     --- If only initialization is being done (iwhich==1) then the bvp3d_work
c     --- routine only needs to be called once. Proper arrays are still passed
c     --- though they should never be needed during initialization.
      idmax = 2
      if (iwhich == 1) idmax = 0
      do id=0,idmax
        if (bfstype <= 4) then
          atemp(:,:,0:bfield%nz) = bfield%j(id,:,:,:)
        else
          atemp = bfield%a(id,0:bfield%nx,0:bfield%ny,:)
          jtemp = bfield%j(id,:,:,:)
        endif
        call bvp3d_work(iwhich,bfstype,atemp,jtemp,id,bfield)
        if (bfstype .ne. 10) bfield%a(id,0:bfield%nx,0:bfield%ny,:) = atemp
      enddo

c     --- Rescale, replacing 1/eps0 with mu0
c     bfield%a = bfield%a*mu0*eps0

c     --- This is slightly inefficient in some cases, since for example, the
c     --- MG solver already takes care of the longitudinal BC's.
      call setaboundaries3d(bfield)

      deallocate(atemp)
      if (bfstype > 4) deallocate(jtemp)

      bfield%j = bfield%j/(mu0*eps0)

      return
      end
c=============================================================================
      subroutine bvp3d_work(iwhich,bfstype,a,j,iaxis,bfield)
      use Subtimersf3d
      use BFieldGridTypemodule
      use LatticeInternal, only: linbend
      use Picglb, only: zbeam,zgrid
      use InGen3d, only: filt,solvergeom,XYZgeomMR,RZgeom,l2symtry,l4symtry
      use GridBoundary3d, only: bound0,boundnz,boundxy
      use PSOR3d, only: gridmode
      use Parallel, only: my_index,nslaves,izfsslave,nzfsslave
      integer(ISZ):: iwhich,bfstype,iaxis
      type(BFieldGridType):: bfield
      real(kind=8):: a(0:bfield%nx,0:bfield%ny,-1:bfield%nz+1)
      real(kind=8):: j(0:bfield%nx,0:bfield%ny,0:bfield%nz)

c Interface to field solver for magnetic vector potential
c Note that bfield is passed in to have access to the data and work arrays
c kxsq, attx, scrtch etc.

      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: dx,dy,dz
      real(kind=8):: xlen,ylen,zlen
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (bfstype >= 0) then
        if (solvergeom == RZgeom) bfstype = 10
        if (solvergeom == XYZgeomMR) bfstype = 12
      endif

      if (iwhich==0 .or. iwhich==1) call execuser('initfieldsolver')

c Note that the select case statement was giving the f90 compiler on HPUX
c fits and so was replaced with if statements.

      nx = bfield%nx
      ny = bfield%ny
      nz = bfield%nz
      nzfull = bfield%nzfull
      dx = bfield%dx
      dy = bfield%dy
      dz = bfield%dz
      xlen = bfield%xmmax - bfield%xmmin
      ylen = bfield%ymmax - bfield%ymmin
      zlen = bfield%zmmaxglobal - bfield%zmminglobal

      if (bfstype == 0) then
        call vpois3d(iwhich,a(0,0,0),a(0,0,0),
     &               bfield%kxsq,bfield%kysq,bfield%kzsq,
     &               bfield%attx,bfield%atty,bfield%attz,
     &               filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &               bfield%scrtch,bfield%xywork,bfield%zwork,
     &               0,l2symtry,l4symtry,bound0,boundnz,boundxy)
      elseif (bfstype == 4) then
        if (iwhich == 1 .or. iwhich == 0) then
          call vpois3d(1,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
        endif
        if (iwhich == -1 .or. iwhich == 0) then
          call vpois3d(12,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
          call vpois3d(14,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
          call vpois3d(13,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
        endif

      elseif (bfstype == 7) then

        call multigrid3dsolve(iwhich,nx,ny,nz,nzfull,dx,dy,dz,a,j,
     &                      bfield%rstar,linbend,bfield%bounds,
     &                      bfield%xmmin,bfield%ymmin,bfield%zmmin,zbeam,zgrid,
     &                      bfield%mgparam(iaxis),bfield%mgform(iaxis),
     &                      bfield%mgiters(iaxis),bfield%mgmaxiters(iaxis),
     &                      bfield%mgmaxlevels(iaxis),bfield%mgerror(iaxis),
     &                      bfield%mgtol(iaxis),
     &                      bfield%downpasses(iaxis),bfield%uppasses(iaxis),
     &                      bfield%lcndbndy,bfield%laddconductor,
     &                      bfield%icndbndy,.false.,
     &                      gridmode,bfield%conductors,
     &                      my_index,nslaves,izfsslave,nzfsslave)


      elseif (bfstype == 10) then
        call multigridrzb(iwhich,iaxis,bfield%a(iaxis,:,0,:),bfield%j(iaxis,:,0,:),
     &                    bfield%nx,bfield%nz)
      elseif (bfstype == 12) then
        if (iwhich <= 0) then
          call execuser('fieldsolMR')
        endif

      endif

!$OMP MASTER
      if (lf3dtimesubs) timebvp3d = timebvp3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
