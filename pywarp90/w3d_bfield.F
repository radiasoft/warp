#include "top.h"
c=============================================================================
c@(#) File W3D_SELFB.F, version $Revision: 1.3 $, $Date: 2005/08/19 16:28:23 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This contains routines to handle the B field, calculated from the current
c  density.
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine init_bfieldsolver(bfstype)
      use BFieldGrid
      use InMesh3d
      use Fields3d, only: nzfull

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
      if (bfield%dy == 0) bfield%dy = (bfield%ymmax-bfield%ymmin)/bfield%ny
      if (bfield%dz == 0) bfield%dz = (bfield%zmmax-bfield%zmmin)/bfield%nz
      call BFieldGridTypechange(bfield)

      call bvp3d(1,bfstype)

      return
      end
c=============================================================================
      subroutine loadj3d(ins_i,nps_i,is_i,lzero)
      use GlobalVars
      use Subtimers3d
      use InGen
      use InGen3d
      use InPart
      use Picglb
      use Picglb3d
      use Particles
      use BFieldGrid
      use GridBoundary3d
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

      if (lw3dtimesubs) substarttime = wtime()

c     --- Ensure that the jp array is setup properly
      if(solvergeom==XYZgeom) then
        call setupbfieldsforparticles3d(ns,ndts,it,bfield,bfieldp)
      endif

c     --- zero jp if requested
      if (lzero) then
        if (solvergeom == XYZgeom) then
          call zeroarry(bfieldp%j,
     &                  3*(bfieldp%nx+1)*(bfieldp%ny+1)*(bfieldp%nz+1))
        endif
c       if(solvergeom==RZgeom .or. solvergeom==XZgeom .or.
c    &     solvergeom==Rgeom  .or. solvergeom==Zgeom) call reset_rzmgrid_j()
      end if

c     --- change AMR grid if necessary
c     if(solvergeom==RZgeom .or. solvergeom==XZgeom .or.
c    &   solvergeom==Rgeom  .or. solvergeom==Zgeom) call change_loc_part()

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

c        --- For the AMR/Chombo version, all of the lost particles must be
c        --- cleared out.
c        if (solvergeom == AMRgeom) call clearpart(is,1)

c        --- get loop limits for particles if needed
         if (ins_i == -1) ins_u = ins(is)
         if (nps_i == -1) nps_u = nps(is)

c        --- Scale the weight, sw, by the time step scale size. This only
c        --- makes sense for steady-state and slice modes. In time-dependent
c        --- mode, it is assumed that dtscale has not been changed from 1.
         swtmp = sw(is)*dtscale(is)

         if(solvergeom==XYZgeom) then
           ipmin = ins_u
           if (nps_u > 0) then
             if (ndts(is)>1) then
c              --- This code would be better if it were in setj3d, but since
c              --- jp is passed into it, this must be done here.
               if (mod(it+1,ndts(is))==0) then
c                 --- If this species is being advanced this step, then
c                 --- deposit its j into jptmp.
                  jpcopy => bfieldp%j
                  bfieldp%j => bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
                  bfieldp%j = 0.
               else
c                 --- If this species is not being advanced this step, then
c                 --- just add its saved jptmp into jp.
                  bfieldp%j = bfieldp%j + bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
                  cycle
                end if
             end if
             call setj3d(bfieldp,bfieldp%j,
     &                   nps_u,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
     &                   uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),sq(is),
     &                   swtmp,depos,l2symtry,l4symtry)
             if (ndts(is)>1 .and. mod(it+1,ndts(is))==0) then
c               --- If this species is being advanced this step, then restore
c               --- jp and copy its j into jp.
                bfieldp%j => jpcopy
                bfieldp%j = bfieldp%j + bfieldp%jtmp(:,:,:,:,bfieldp%jsjtmp(is-1))
             end if
           endif
         elseif(solvergeom==RZgeom .or. solvergeom==XZgeom) then
           print*,"Error: Self magnetic field not support with solvergeom=RZgeom or XZgeom"
c          --- loop over particle blocks
c          do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
c            ip = min(nparpgrp, ins_u+nps_u-ipmin)
c            if(wpid==0) then
c              call jweightrz(xp(ipmin),yp(ipmin),zp(ipmin),ip,
c    &                          sq(is)*swtmp,nx,nz,dx,dz,xmmin,zgrid)
c            else
c              call jweightrz_weights(xp(ipmin),yp(ipmin),zp(ipmin),pid(ipmin,wpid),ip,
c    &                          sq(is)*swtmp,nx,nz,dx,dz,xmmin,zgrid)
c            end if
c          enddo
         elseif(solvergeom==Zgeom) then
           print*,"Error: Self magnetic field not support with solvergeom=Zgeom"
c          --- loop over particle blocks
c          do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
c             ip = min(nparpgrp, ins_u+nps_u-ipmin)
c             call jweightz(zp(ipmin),ip,sq(is)*swtmp,nz,dz,zgrid)
c          enddo
         elseif(solvergeom==Rgeom) then
           print*,"Error: Self magnetic field not support with solvergeom=Rgeom"
c          --- loop over particle blocks
c          do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
c             ip = min(nparpgrp, ins_u+nps_u-ipmin)
c             if(wpid==0) then
c               call jweightr(xp(ipmin),yp(ipmin),ip,sq(is)*swtmp,nx,dx,xmmin)
c             else
c               call jweightr_weights(xp(ipmin),yp(ipmin),pid(ipmin,wpid),ip,
c    &                                  sq(is)*swtmp,nx,dx,xmmin)
c             end if
c          enddo
         elseif(solvergeom==AMRgeom) then
           print*,"Error: Self magnetic field not support with solvergeom=AMRgeom"
c          --- loop over particle blocks
c          do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
c             ip = min(nparpgrp, ins_u+nps_u-ipmin)
c             call cho_setj3d(ip,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
c    &                        sq(is),swtmp,is,(ipmin-ins(is))/nparpgrp)
c          enddo
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
      if (lw3dtimesubs) timeloadj3d = timeloadj3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setaboundaries3d(bfield)
      use Subtimers3d
      use GlobalVars
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield

c  Sets the slices on the exterior of A, depending on the boundary conditions.

      real(kind=8),pointer,dimension(:,:,:,:):: a
      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c       call perarz()
c       return
c     end if

      a => bfield%a

      if (bfield%bounds(0) == periodic .or. bfield%bounds(1) == periodic) then
        a(:,-1,:,:)   = a(:,bfield%nx-1,:,:)
        a(:,bfield%nx,:,:)   = a(:,0,:,:)
        a(:,bfield%nx+1,:,:) = a(:,1,:,:)
      else
        if (bfield%bounds(0) == neumann) then
          a(:,-1,:,:) = a(:,1,:,:)
        elseif (bfield%bounds(0) == dirichlet) then
          a(:,-1,:,:) = a(:,0,:,:)
        endif
        if (bfield%bounds(1) == neumann) then
          a(:,bfield%nx+1,:,:) = a(:,bfield%nx-1,:,:)
        elseif (bfield%bounds(1) == dirichlet) then
          a(:,bfield%nx+1,:,:) = a(:,bfield%nx,:,:)
        endif
      endif

      if (bfield%bounds(2) == periodic .or. bfield%bounds(3) == periodic) then
        a(:,:,-1,:)   = a(:,:,bfield%ny-1,:)
        a(:,:,bfield%ny,:)   = a(:,:,0,:)
        a(:,:,bfield%ny+1,:) = a(:,:,1,:)
      else
        if (bfield%bounds(2) == neumann) then
          a(:,:,-1,:) = a(:,:,1,:)
        elseif (bfield%bounds(2) == dirichlet) then
          a(:,:,-1,:) = a(:,:,0,:)
        endif
        if (bfield%bounds(3) == neumann) then
          a(:,:,bfield%ny+1,:) = a(:,:,bfield%ny-1,:)
        elseif (bfield%bounds(3) == dirichlet) then
          a(:,:,bfield%ny+1,:) = a(:,:,bfield%ny,:)
        endif
      endif

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
          a(:,:,:,-1) = a(:,:,:,0)
        endif
        if (bfield%bounds(5) == neumann) then
          a(:,:,:,bfield%nz+1) = a(:,:,:,bfield%nz-1)
        elseif (bfield%bounds(5) == dirichlet) then
          a(:,:,:,bfield%nz+1) = a(:,:,:,bfield%nz)
        endif
      endif

!$OMP MASTER
      if (lw3dtimesubs) timepera3d = timepera3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine perj3d(bfield,bound0,boundxy)
      use BFieldGridTypemodule
      use GlobalVars
      use Subtimers3d
      use InGen3d
      type(BFieldGridType):: bfield
      integer(ISZ):: bound0,boundxy

c  Sums the first and last slices of j for periodicity
c  and puts the result into both slices.

      real(kind=8),pointer:: j(:,:,:,:)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

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

c  Copy charge density from frz.basegrid to w3d.j
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
      if (lw3dtimesubs) timeperj3d = timeperj3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setb3d(bfield,np,xp,yp,zp,zgrid,bx,by,bz,l2symtry,l4symtry)
      use Subtimers3d
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield
      integer(ISZ):: np
      real(kind=8):: zgrid
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: bx(np),by(np),bz(np)
      logical(ISZ):: l2symtry,l4symtry

c  Sets magnetic field for particles

      integer(ISZ):: ip,i,j,k
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      real(kind=8):: sx,sy

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Evaluation of B, vectorized over particles
      tdxi = 1./(2.*bfield%dx)
      tdyi = 1./(2.*bfield%dy)
      tdzi = 1./(2.*bfield%dz)
      dxi = 1./bfield%dx
      dyi = 1./bfield%dy
      dzi = 1./bfield%dz

      if (.not. (l2symtry .or. l4symtry)) then
        do ip = 1, np

          i = (xp(ip) - bfield%xmmin) * dxi
          j = (yp(ip) - bfield%ymmin) * dyi
          k = (zp(ip) - zgrid - bfield%zmmin) * dzi

          u1 = (xp(ip) - bfield%xmmin) * dxi - i
          v1 = (yp(ip) - bfield%ymmin) * dyi - j
          w1 = (zp(ip) - zgrid - bfield%zmmin) * dzi - k

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

          i = (abs(xp(ip)) - bfield%xmmin)*dxi
          j = (abs(yp(ip)) - bfield%ymmin)*dyi
          k = (zp(ip) - zgrid - bfield%zmmin)*dzi

          u1 = (abs(xp(ip)) - bfield%xmmin)*dxi - i
          v1 = (abs(yp(ip)) - bfield%ymmin)*dyi - j
          w1 = (zp(ip) - zgrid - bfield%zmmin)*dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Adjust sign of B field for appropriate quadrant.
          xsign = +1.
          ysign = +1.
          if (xp(ip) < 0.) xsign = sx
          if (yp(ip) < 0.) ysign = sy

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

!$OMP MASTER
      if (lw3dtimesubs) timesetb3d = timesetb3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fetchafrompositions3d(np,xp,yp,zp,zgrid,bfield,
     &                                 l2symtry,l4symtry)
      use Subtimers3d
      use BFieldGridTypemodule
      integer(ISZ):: np
      real(kind=8):: zgrid
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: a(0:2,np)
      type(BFieldGridType):: bfield
      logical(ISZ):: l2symtry,l4symtry

c  Gets the magnetic vector potential

      integer(ISZ):: ip,i,j,k
      real(kind=8):: dxi,dyi,dzi,tdxi,tdyi,tdzi,u0,u1,v0,v1,w0,w1,ysign,xsign
      real(kind=8):: sx,sy

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Evaluation of A, vectorized over particles
      tdxi = 1./(2.*bfield%dx)
      tdyi = 1./(2.*bfield%dy)
      tdzi = 1./(2.*bfield%dz)
      dxi = 1./bfield%dx
      dyi = 1./bfield%dy
      dzi = 1./bfield%dz

      if (.not. (l2symtry .or. l4symtry)) then
        do ip = 1, np

          i = (xp(ip) - bfield%xmmin) * dxi
          j = (yp(ip) - bfield%ymmin) * dyi
          k = (zp(ip) - zgrid - bfield%zmmin) * dzi

          u1 = (xp(ip) - bfield%xmmin) * dxi - i
          v1 = (yp(ip) - bfield%ymmin) * dyi - j
          w1 = (zp(ip) - zgrid - bfield%zmmin) * dzi - k

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

          i = (abs(xp(ip)) - bfield%xmmin)*dxi
          j = (abs(yp(ip)) - bfield%ymmin)*dyi
          k = (zp(ip) - zgrid - bfield%zmmin)*dzi

          u1 = (abs(xp(ip)) - bfield%xmmin)*dxi - i
          v1 = (abs(yp(ip)) - bfield%ymmin)*dyi - j
          w1 = (zp(ip) - zgrid - bfield%zmmin)*dzi - k

          u0 = 1. - u1
          v0 = 1. - v1
          w0 = 1. - w1

c         --- Adjust sign of B field for appropriate quadrant.
          xsign = +1.
          ysign = +1.
          if (xp(ip) < 0.) xsign = sx
          if (yp(ip) < 0.) ysign = sy

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

!$OMP MASTER
      if (lw3dtimesubs) timefetchafrompositions3d = timefetchafrompositions3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=========================================================================== 
      subroutine getbfroma3d(bfield)
      use GlobalVars
      use Subtimers3d
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield

c Calculate the B via finite differences of A.

      integer(ISZ):: ix,iy,iz
      real(kind=8):: tdxi,tdyi,tdzi
      real(kind=8),pointer,dimension(:,:,:,:):: a,b
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

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
            b(0,ix,iy,iz) = + tdyi*(a(2,ix  ,iy-1,iz  ) - a(2,ix  ,iy+1,iz  ))
     &                      - tdzi*(a(1,ix  ,iy  ,iz-1) - a(1,ix  ,iy  ,iz+1))
            b(1,ix,iy,iz) = + tdzi*(a(0,ix  ,iy  ,iz-1) - a(0,ix  ,iy  ,iz+1))
     &                      - tdxi*(a(2,ix-1,iy  ,iz  ) - a(2,ix+1,iy  ,iz  ))
            b(2,ix,iy,iz) = + tdxi*(a(1,ix-1,iy  ,iz  ) - a(1,ix+1,iy  ,iz  ))
     &                      - tdyi*(a(0,ix  ,iy-1,iz  ) - a(0,ix  ,iy+1,iz  ))
          enddo
        enddo
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timegetbfroma3d = timegetbfroma3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine setj3d(bfield,j1d,np,xp,yp,zp,zgrid,uxp,uyp,uzp,gaminv,
     &                  q,wght,depos,l2symtry,l4symtry)
      use GlobalVars
      use Subtimers3d
      use BFieldGridTypemodule
      type(BFieldGridType):: bfield
      integer(ISZ):: np
      real(kind=8):: zgrid,q,wght
      real(kind=8):: j1d(0:3*(1+bfield%nx)*(1+bfield%ny)*(1+bfield%nz)-1)
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
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
      integer(ISZ):: moff(0:7), indx(0:7,0:nparpgrp-1)
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(0:nparpgrp-1), jj(0:nparpgrp-1), kk(0:nparpgrp-1)
c     --- For both algorithms
      real(kind=8):: s(0:7,0:nparpgrp-1)

      real(kind=8),pointer,dimension(:,:,:,:):: j
      integer(ISZ):: ipmin,nptmp,ip,ix,iy,iz,ind0,m,ir
      real(kind=8):: g,dxi,dyi,dzi,u0,u1,v0,v1,w0,w1,gxfact,gyfact
      real(kind=8):: vx,vy,vz
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c  Set up offset array for vectorized deposition:

      moff(0) = 0
      moff(1) = 1
      moff(2) = bfield%nx+1
      moff(3) = bfield%nx+2
      moff(4) = (bfield%nx+1)*(bfield%ny+1)
      moff(5) = (bfield%nx+1)*(bfield%ny+1)+1
      moff(6) = (bfield%nx+1)*(bfield%ny+1)+bfield%nx+1
      moff(7) = (bfield%nx+1)*(bfield%ny+1)+bfield%nx+2

      g = wght * q / (bfield%dx * bfield%dy * bfield%dz)
      dxi = 1./bfield%dx
      dyi = 1./bfield%dy
      dzi = 1./bfield%dz
      if (l2symtry) then
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        g = g*0.5
      elseif (l4symtry) then
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        g = g*0.25
      endif

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ipmin,nptmp,i,j,k,u1,u0,v1,v0,w1,w0,ir,ip,ind0,indx,
!$OMP&                 gyfact,gxfact,s,m,ii,jj,kk)
!$OMP&FIRSTPRIVATE(np)

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
c       --- special loop for 2-fold symmetry
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        do ip = ipmin,ipmin+nptmp-1
           ix  = int((xp(ip) - bfield%xmmin) * dxi)
           u1 =      (xp(ip) - bfield%xmmin) * dxi - ix
           u0 = 1. - u1
           iy  = int((abs(yp(ip)) - bfield%ymmin)*dyi)
           v1 =      (abs(yp(ip)) - bfield%ymmin)*dyi - iy
           v0 = 1. - v1
           iz  = int((zp(ip) - zgrid - bfield%zmmin) * dzi)
           w1 =      (zp(ip) - zgrid - bfield%zmmin) * dzi - iz
           w0 = 1. - w1
           ir = ip - ipmin
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
           s(0,ir) = u0 * v0 * w0 * g*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        do ip = ipmin,ipmin+nptmp-1
           ix  = int((abs(xp(ip)) - bfield%xmmin)*dxi)
           u1 =      (abs(xp(ip)) - bfield%xmmin)*dxi - ix
           u0 = 1. - u1
           iy  = int((abs(yp(ip)) - bfield%ymmin)*dyi)
           v1 =      (abs(yp(ip)) - bfield%ymmin)*dyi - iy
           v0 = 1. - v1
           iz  = int((zp(ip) - zgrid - bfield%zmmin) * dzi)
           w1 =      (zp(ip) - zgrid - bfield%zmmin) * dzi - iz
           w0 = 1. - w1
           ir = ip - ipmin
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
           s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
           s(1,ir) = u1 * v0 * w0 * g*gyfact
           s(2,ir) = u0 * v1 * w0 * g*gxfact
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
           s(5,ir) = u1 * v0 * w1 * g*gyfact
           s(6,ir) = u0 * v1 * w1 * g*gxfact
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
           ix  = int((xp(ip) - bfield%xmmin) * dxi)
           u1 =      (xp(ip) - bfield%xmmin) * dxi - ix
           u0 = 1. - u1
           iy  = int((yp(ip) - bfield%ymmin) * dyi)
           v1 =      (yp(ip) - bfield%ymmin) * dyi - iy
           v0 = 1. - v1
           iz  = int((zp(ip) - zgrid - bfield%zmmin) * dzi)
           w1 =      (zp(ip) - zgrid - bfield%zmmin) * dzi - iz
           w0 = 1. - w1
           ir = ip - ipmin
           ind0 = ix + iy*(bfield%nx+1) + iz*(bfield%nx+1)*(bfield%ny+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           indx(4,ir) = ind0 + moff(4)
           indx(5,ir) = ind0 + moff(5)
           indx(6,ir) = ind0 + moff(6)
           indx(7,ir) = ind0 + moff(7)
           s(0,ir) = u0 * v0 * w0 * g
           s(1,ir) = u1 * v0 * w0 * g
           s(2,ir) = u0 * v1 * w0 * g
           s(3,ir) = u1 * v1 * w0 * g
           s(4,ir) = u0 * v0 * w1 * g
           s(5,ir) = u1 * v0 * w1 * g
           s(6,ir) = u0 * v1 * w1 * g
           s(7,ir) = u1 * v1 * w1 * g
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETJ3D1)
      do ir = 0,nptmp-1
        if (uzp(ir+ipmin) /= 0.) then
          do m = 0, 7
             ind0 = 3.*indx(m,ir)
             j1d(ind0  ) = j1d(ind0  ) + s(m,ir)*uxp(ir+ipmin)*gaminv(ir+ipmin)
             j1d(ind0+1) = j1d(ind0+1) + s(m,ir)*uyp(ir+ipmin)*gaminv(ir+ipmin)
             j1d(ind0+2) = j1d(ind0+2) + s(m,ir)*uzp(ir+ipmin)*gaminv(ir+ipmin)
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
        do ip = ipmin,ipmin+nptmp-1
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - bfield%xmmin) * dxi)
          u1     =     (xp(ip) - bfield%xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - bfield%ymmin)*dyi)
          v1     =     (abs(yp(ip)) - bfield%ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - bfield%zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - bfield%zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gyfact = 1.
          if (jj(ir) == 0 .and. bfield%ymmin == 0.) gyfact = 2.
          s(0,ir) = u0 * v0 * w0 * g*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      elseif (l4symtry) then
        do ip = ipmin,ipmin+nptmp-1
          ir = ip - ipmin
          ii(ir) = int((abs(xp(ip)) - bfield%xmmin)*dxi)
          u1     =     (abs(xp(ip)) - bfield%xmmin)*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((abs(yp(ip)) - bfield%ymmin)*dyi)
          v1     =     (abs(yp(ip)) - bfield%ymmin)*dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - bfield%zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - bfield%zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0 .and. bfield%xmmin == 0.) gxfact = 2.
          if (jj(ir) == 0 .and. bfield%ymmin == 0.) gyfact = 2.
          s(0,ir) = u0 * v0 * w0 * g*gxfact*gyfact
          s(1,ir) = u1 * v0 * w0 * g*gyfact
          s(2,ir) = u0 * v1 * w0 * g*gxfact
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g*gxfact*gyfact
          s(5,ir) = u1 * v0 * w1 * g*gyfact
          s(6,ir) = u0 * v1 * w1 * g*gxfact
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          ir = ip - ipmin
          ii(ir) = int((xp(ip) - bfield%xmmin) * dxi)
          u1     =     (xp(ip) - bfield%xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = int((yp(ip) - bfield%ymmin) * dyi)
          v1     =     (yp(ip) - bfield%ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          kk(ir) = int((zp(ip) - zgrid - bfield%zmmin) * dzi)
          w1     =     (zp(ip) - zgrid - bfield%zmmin) * dzi - kk(ir)
          w0     = 1. - w1
          s(0,ir) = u0 * v0 * w0 * g
          s(1,ir) = u1 * v0 * w0 * g
          s(2,ir) = u0 * v1 * w0 * g
          s(3,ir) = u1 * v1 * w0 * g
          s(4,ir) = u0 * v0 * w1 * g
          s(5,ir) = u1 * v0 * w1 * g
          s(6,ir) = u0 * v1 * w1 * g
          s(7,ir) = u1 * v1 * w1 * g
        enddo
      endif
c     --- scalar loop does the actual deposition
      j => bfield%j
!$OMP CRITICAL (CRITICAL_SETJ3D2)
      do ir = 0, nptmp-1
         if (uzp(ir+ipmin) /= 0) then
         ix = ii(ir)
         iy = jj(ir)
         iz = kk(ir)
         vx = uxp(ir+ipmin)*gaminv(ir+ipmin)
         vy = uyp(ir+ipmin)*gaminv(ir+ipmin)
         vz = uzp(ir+ipmin)*gaminv(ir+ipmin)
         j(0,ix  ,iy  ,iz  ) = j(0,ix  ,iy  ,iz  ) + s(0,ir)*vx
         j(0,ix+1,iy  ,iz  ) = j(0,ix+1,iy  ,iz  ) + s(1,ir)*vx
         j(0,ix  ,iy+1,iz  ) = j(0,ix  ,iy+1,iz  ) + s(2,ir)*vx
         j(0,ix+1,iy+1,iz  ) = j(0,ix+1,iy+1,iz  ) + s(3,ir)*vx
         j(0,ix  ,iy  ,iz+1) = j(0,ix  ,iy  ,iz+1) + s(4,ir)*vx
         j(0,ix+1,iy  ,iz+1) = j(0,ix+1,iy  ,iz+1) + s(5,ir)*vx
         j(0,ix  ,iy+1,iz+1) = j(0,ix  ,iy+1,iz+1) + s(6,ir)*vx
         j(0,ix+1,iy+1,iz+1) = j(0,ix+1,iy+1,iz+1) + s(7,ir)*vx
         j(1,ix  ,iy  ,iz  ) = j(1,ix  ,iy  ,iz  ) + s(0,ir)*vy
         j(1,ix+1,iy  ,iz  ) = j(1,ix+1,iy  ,iz  ) + s(1,ir)*vy
         j(1,ix  ,iy+1,iz  ) = j(1,ix  ,iy+1,iz  ) + s(2,ir)*vy
         j(1,ix+1,iy+1,iz  ) = j(1,ix+1,iy+1,iz  ) + s(3,ir)*vy
         j(1,ix  ,iy  ,iz+1) = j(1,ix  ,iy  ,iz+1) + s(4,ir)*vy
         j(1,ix+1,iy  ,iz+1) = j(1,ix+1,iy  ,iz+1) + s(5,ir)*vy
         j(1,ix  ,iy+1,iz+1) = j(1,ix  ,iy+1,iz+1) + s(6,ir)*vy
         j(1,ix+1,iy+1,iz+1) = j(1,ix+1,iy+1,iz+1) + s(7,ir)*vy
         j(2,ix  ,iy  ,iz  ) = j(2,ix  ,iy  ,iz  ) + s(0,ir)*vz
         j(2,ix+1,iy  ,iz  ) = j(2,ix+1,iy  ,iz  ) + s(1,ir)*vz
         j(2,ix  ,iy+1,iz  ) = j(2,ix  ,iy+1,iz  ) + s(2,ir)*vz
         j(2,ix+1,iy+1,iz  ) = j(2,ix+1,iy+1,iz  ) + s(3,ir)*vz
         j(2,ix  ,iy  ,iz+1) = j(2,ix  ,iy  ,iz+1) + s(4,ir)*vz
         j(2,ix+1,iy  ,iz+1) = j(2,ix+1,iy  ,iz+1) + s(5,ir)*vz
         j(2,ix  ,iy+1,iz+1) = j(2,ix  ,iy+1,iz+1) + s(6,ir)*vz
         j(2,ix+1,iy+1,iz+1) = j(2,ix+1,iy+1,iz+1) + s(7,ir)*vz
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETJ3D2)

      endif

      enddo
!$OMP END DO

!$OMP END PARALLEL

      if (lw3dtimesubs) timesetj3d = timesetj3d + wtime() - substarttime
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
      if(solvergeom==RZgeom .or. solvergeom==XZgeom) then
c       print*,"Error: Self magnetic field not supported with solvergeom=RZgeom or XZgeom"
c       call getjforfieldsolverz(nx,nz,j)
c     else if(solvergeom==Zgeom) then
c       print*,"Error: Self magnetic field not supported with solvergeom=Zgeom"
c       call getjforfieldsolvez(nz,j)
      elseif (solvergeom==XYZgeom) then
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
      bfieldp%bounds = bfield%bounds

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
      use Subtimers3d
      use Picglb
      use Picglb3d
      use InGen
      use InGen3d
      use BFieldGrid
      use FieldSolveAPI
      integer(ISZ):: is,n
      real(kind=8),target:: x(n),y(n),z(n),bx(n),by(n),bz(n)

c     --- Obtain the magnetic field

      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

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
      if (lw3dtimesubs) timefetchb3dfrompositions = timefetchb3dfrompositions + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
      subroutine fetcha(n,x,y,z,a)
      use BFieldGrid
      use InGen3d
      use Picglb
      use FieldSolveAPI
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

      if(solvergeom==XYZgeom .or. solvergeom==XYZgeomMR) then
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
      if(solvergeom==RZgeom .or. solvergeom==XZgeom .or. solvergeom==Zgeom) then
c       call getbforparticlesrz()
c      else if(solvergeom==Zgeom) then
c        call getbforparticlesz()
      elseif (solvergeom==XYZgeom) then
        call getbforparticles3d(bfield,bfieldp)
      end if
#endif

      return
      end
c===========================================================================
      subroutine bfieldsol3d(iwhich)
      use GlobalVars
      use Subtimers3d
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
      if (lw3dtimesubs) substarttime = wtime()
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
      if (lw3dtimesubs) timebfieldsol3d = timebfieldsol3d + wtime() - substarttime
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
        call bvp3d_work(iwhich,bfstype,
     &                  atemp,jtemp,
     &                  bfield%nx,bfield%ny,bfield%nz,bfield%nzfull,
     &                  bfield%dx,bfield%dy,bfield%dz,
     &                  bfield%xmmin,bfield%xmmax,bfield%ymmin,bfield%ymmax,
     &                  bfield%zmmin,bfield%zmmax,
     &                  bfield%zmminglobal,bfield%zmmaxglobal,
     &                  bfield)
         bfield%a(id,0:bfield%nx,0:bfield%ny,:) = atemp
      enddo

c     --- Rescale, replacing 1/eps0 with mu0
      bfield%a = bfield%a*mu0*eps0

c     --- This is slightly inefficient in some cases, since for example, the
c     --- MG solver already takes care of the longitudinal BC's.
      call setaboundaries3d(bfield)

      deallocate(atemp)
      if (bfstype > 4) deallocate(jtemp)

      return
      end
c=============================================================================
      subroutine bvp3d_work(iwhich,bfstype,a,j,nx,ny,nz,nzfull,dx,dy,dz,
     &                      xmmin,xmmax,ymmin,ymmax,zmmin,zmmax,
     &                      zmminglobal,zmmaxglobal,
     &                      bfield)
      use Subtimers3d
      use BFieldGridTypemodule
      use LatticeInternal, only: linbend
      use Picglb, only: zbeam,zgrid
      use InGen3d, only: filt,solvergeom,XYZgeomMR,l2symtry,l4symtry
      use GridBoundary3d, only: bound0,boundnz,boundxy
      integer(ISZ):: iwhich,bfstype
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: a(0:nx,0:ny,-1:nz+1),j(0:nx,0:ny,0:nz)
      real(kind=8):: dx,dy,dz,xmmin,xmmax,ymmin,ymmax,zmmin,zmmax
      real(kind=8):: zmminglobal,zmmaxglobal
      type(BFieldGridType):: bfield

c Interface to field solver for magnetic vector potential
c Note that bfield is passed in to have access to the data and work arrays
c kxsq, attx, scrtch etc.

      real(kind=8):: xlen,ylen,zlen
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()


c LatticeInternal is included to have access to the variable linbend.

      if (solvergeom == XYZgeomMR) bfstype = 12

      if (iwhich==0 .or. iwhich==1) call execuser('initfieldsolver')

c Note that the select case statement was giving the f90 compiler on HPUX
c fits and so was replaced with if statements.

      xlen = xmmax - xmmin
      ylen = ymmax - ymmin
      zlen = zmmaxglobal - zmminglobal

      if (bfstype == 0) then
          call vpois3d(iwhich,a(0,0,0),a(0,0,0),
     &                 bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                 bfield%attx,bfield%atty,bfield%attz,
     &                 filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                 bfield%scrtch,bfield%xywork,bfield%zwork,
     &                 0,l2symtry,l4symtry,bound0,boundnz,boundxy)
      elseif (bfstype == 4) then
          if (iwhich == 1 .or. iwhich == 0) then
            call vpois3d(1,a(0,0,0),a(0,0,0),
     &                   bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                   bfield%attx,bfield%atty,bfield%attz,
     &                   filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                   bfield%scrtch,bfield%xywork,bfield%zwork,
     &                   0,l2symtry,l4symtry,bound0,boundnz,boundxy)
          endif
          if (iwhich == -1 .or. iwhich == 0) then
            call vpois3d(12,a(0,0,0),a(0,0,0),
     &                   bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                   bfield%attx,bfield%atty,bfield%attz,
     &                   filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                   bfield%scrtch,bfield%xywork,bfield%zwork,
     &                   0,l2symtry,l4symtry,bound0,boundnz,boundxy)
            call vpois3d(14,a(0,0,0),a(0,0,0),
     &                   bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                   bfield%attx,bfield%atty,bfield%attz,
     &                   filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                   bfield%scrtch,bfield%xywork,bfield%zwork,
     &                   0,l2symtry,l4symtry,bound0,boundnz,boundxy)
            call vpois3d(13,a(0,0,0),a(0,0,0),
     &                   bfield%kxsq,bfield%kysq,bfield%kzsq,
     &                   bfield%attx,bfield%atty,bfield%attz,
     &                   filt,xlen,ylen,zlen,nx,ny,nz,nzfull,
     &                   bfield%scrtch,bfield%xywork,bfield%zwork,
     &                   0,l2symtry,l4symtry,bound0,boundnz,boundxy)
          endif

      elseif (bfstype == 7) then
          call multigrid3df(iwhich,nx,ny,nz,nzfull,dx,dy,dz,a,j,bfield%rstar,
     &                      linbend,bound0,boundnz,boundxy,l2symtry,l4symtry,
     &                      xmmin,ymmin,zmmin,zbeam,zgrid)

      elseif (bfstype == 12) then
        if (iwhich <= 0) then
          call execuser('fieldsolMR')
        endif

      endif

!$OMP MASTER
      if (lw3dtimesubs) timebvp3d = timebvp3d + wtime() - substarttime
!$OMP END MASTER
      return
      end
c=============================================================================
