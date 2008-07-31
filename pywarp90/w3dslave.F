#include "top.h"
c===========================================================================
c W3DSLAVE.M, version $Revision: 1.48 $, $Date: 2008/07/31 00:45:19 $
c Slave routines associated with the W3D package.
c===========================================================================
c===========================================================================
      subroutine init_w3d_parallel
      use Subtimers3d
      use Parallel
      use InGen
      use InGen3d
      use InDiag
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      use Fields3d
      use Fields3dParticles
      use Io
      use Z_arrays
      use LatticeInternal
      use Z_Moments
      use InjectVars
      use InjectVars3d

c This is the routine which divies up the work among the slaves.

      logical(MPIISZ):: w3dinitialized
      data w3dinitialized/.false./
      save w3dinitialized
      logical(MPIISZ):: initialized
      integer(MPIISZ):: mpierror,nslavestmp,my_indextmp
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c This routine should only be called once.
      if (w3dinitialized) return
      w3dinitialized = .true.

c get number of slaves and my_pe
      call MPI_INITIALIZED(initialized,mpierror)
      if (.not. initialized) call MPI_INIT(mpierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nslavestmp,mpierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_indextmp,mpierror)
      nslaves = nslavestmp
      my_index = my_indextmp
      maxslaves = nslaves
      call gchange("Parallel",0)

      if (my_index > 0) then
        verbosity = 0
        lprntpara = .false.
      endif

c---------------------------------------------------------------------------
c     --- Calculate how the work is to be arranged among the processors.

      call domaindecomposefields(nz,nslaves,lfsautodecomp,
     &                           izfsslave,nzfsslave,grid_overlap)

c     --- Set the grid cell size, which is needed below. It must be set here
c     --- after the field solve domain decomposition since nzlocal may be
c     --- determined from the user supplied decomposition.
      dz = (zmmax - zmmin)/nz

      call domaindecomposeparticles(nz,nslaves,izfsslave,nzfsslave,
     &                              grid_overlap,nzpguard,zmmin,zmmax,
     &                              dz,zslave,lautodecomp,izpslave,nzpslave,
     &                              zpslmin,zpslmax)

c---------------------------------------------------------------------------
c     --- Reset local values
      if (injctspc > 0) then
        injctspc = max(1,int(injctspc/nslaves*1.1))
      endif
      nzlocal = nzfsslave(my_index)
      nzp = nzpslave(my_index)
      zmminp = zmmin + izpslave(my_index)*dz
      zmmaxp = zmmin + (izpslave(my_index) + nzpslave(my_index))*dz
      izfsmin = 0
      izfsmax = nzfsslave(my_index)
      zmminlocal = zmmin + izfsslave(my_index)*dz
      zmmaxlocal = zmmin + (izfsslave(my_index) + nzfsslave(my_index))*dz

      xpminlocal = xmmin
      xpmaxlocal = xmmax
      ypminlocal = ymmin
      ypmaxlocal = ymmax
      zpminlocal = zpslmin(my_index)
      zpmaxlocal = zpslmax(my_index)

!$OMP MASTER
      if (lw3dtimesubs) timeinit_w3d_parallel = timeinit_w3d_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
      subroutine sw_globalsum(ns,sw)
      use Subtimers3d
      integer(ISZ):: ns
      real(kind=8):: sw(ns)
c As calculated in stptl3d, sw depends on 1 over the number of particles
c in each processor. The new value of sw is 1 over the sum of the
c reciprocals of the sw of each processor, and so depends only on the
c total number of particles.
      real(kind=8):: swtemp(ns)
      integer(ISZ):: is
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo
      swtemp = sw
      call MPI_ALLREDUCE(swtemp,sw,int(ns,MPIISZ),MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesw_globalsum = timesw_globalsum + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
c===========================================================================
      subroutine sumsourcepondomainboundaries(sourcep,nc,nxp,nyp,nzp,nzpguard,
     &                                        my_index,nslaves,izpslave,nzpslave)
      use Subtimers3d
      integer(ISZ):: nc,nxp,nyp,nzp,nzpguard
      real(kind=8):: sourcep(0:nc-1,0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)

      integer(ISZ):: allocerror
      integer(MPIISZ):: i,nn,i1,i2,iz,nzbuff
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      real(kind=8),allocatable:: sendbuffer(:,:,:,:),recvbuffer(:,:,:,:)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpierror
      integer(MPIISZ):: mpirequest(nslaves)
      integer(MPIISZ):: w
      integer(MPIISZ):: messid = 50
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent.
      do i=0,nslaves-1
        izsend(i) = max(izpslave(my_index),izpslave(i)) - izpslave(my_index)
        nzsend(i) = min(izpslave(my_index)+nzpslave(my_index),
     &                  izpslave(i)+nzpslave(i))
     &             - max(izpslave(my_index),izpslave(i)) + 1
        nzsend(i) = max(int(0,MPIISZ),nzsend(i))
      enddo
      nzsend(my_index) = 0

c     --- Create buffer for MPI work, copying only the first two and last
c     --- two planes of sourcep.
      nzbuff = min(nzp+1,4+4*nzpguard)
      allocate(sendbuffer(nc,0:nxp,0:nyp,nzbuff))
      if (nzbuff == nzp+1) then
        sendbuffer = sourcep
      else
        nn = min(1+2*nzpguard,nzp)
        sendbuffer(:,:,:,1:nn+1) = sourcep(:,:,:,0:nn)
        nn = max(nzp-1-2*nzpguard,nn+1)
        sendbuffer(:,:,:,nzbuff-nzp+nn:nzbuff) = sourcep(:,:,:,nn:nzp)
      endif

c     --- Send the data out to processors that need it.
c     --- Copy any data locally as well.
      w = 0
      do i=0,nslaves-1
        if (i == my_index) cycle
        if (nzsend(i) > 0) then
          w = w + 1
          if (i < my_index) iz = izsend(i) + 1
          if (i > my_index) iz = nzbuff - (nzp - izsend(i))
          nn = nzsend(i)*(nxp+1)*(nyp+1)*nc
          call MPI_ISEND(sendbuffer(:,:,:,iz),nn,MPI_DOUBLE_PRECISION,
     &                   i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      enddo

c     --- Data to be received. Do this calculation after the send so that
c     --- they get posted as soon as possible.
      do i=0,nslaves-1
        izrecv(i) = max(izpslave(my_index),izpslave(i)) - izpslave(my_index)
        nzrecv(i) = min(izpslave(my_index)+nzpslave(my_index),
     &                  izpslave(i)+nzpslave(i))
     &             - max(izpslave(my_index),izpslave(i)) + 1
        nzrecv(i) = max(int(0,MPIISZ),nzrecv(i))
      enddo
      nzrecv(my_index) = 0

c     --- Then, gather up the data sent to this processor.
      allocate(recvbuffer(nc,0:nxp,0:nyp,maxval(nzrecv)),stat=allocerror)
      do i=0,nslaves-1
        if (i == my_index) cycle
        if (nzrecv(i) > 0) then
          nn = nzrecv(i)*(nxp+1)*(nyp+1)*nc
          call MPI_RECV(recvbuffer,nn,MPI_DOUBLE_PRECISION,
     &                  i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          i1 = izrecv(i)
          i2 = izrecv(i) + nzrecv(i) - 1
          sourcep(:,:,:,i1:i2) = sourcep(:,:,:,i1:i2) + recvbuffer(:,:,:,1:nzrecv(i))
        endif
      enddo
      deallocate(recvbuffer)

c     --- Wait for all sends to complete before deleting the sendbuffer.
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
      deallocate(sendbuffer)

!$OMP MASTER
      if (lw3dtimesubs) timesumsourcepondomainboundaries = timesumsourcepondomainboundaries + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine persource3d_slave(source,nc,nx,ny,nzlocal)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nc,nx,ny,nzlocal
      real(kind=8):: source(0:nc-1,0:nx,0:ny,0:nzlocal)

c  Sets the slices on the exterior of source for periodicity
c  sets slice at -1 equal to the slice at nzlocal-1
c  sets slice at nzlocal+1 equal to the slice at 1
c  Only first and last processors do anything.

      real(kind=8),allocatable:: sourcetemp(:,:,:)
      integer(MPIISZ):: nn,pe
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
      integer(MPIISZ):: messid = 70
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      allocate(sourcetemp(0:nc-1,0:nx,0:ny))
 
      nn = nc*(nx+1)*(ny+1)
      if (my_index == nslaves-1) then
        pe = 0
        call MPI_ISEND(source(:,:,:,nzlocal),nn,MPI_DOUBLE_PRECISION,
     &                 pe,messid,MPI_COMM_WORLD,mpirequest,mpierror)
        call MPI_RECV(source(:,:,:,nzlocal),nn,MPI_DOUBLE_PRECISION,
     &                pe,messid,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == 0) then
        pe = nslaves-1
        call MPI_RECV(sourcetemp,nn,MPI_DOUBLE_PRECISION,
     &                pe,messid,MPI_COMM_WORLD,mpistatus,mpierror)
        source(:,:,:,0) = source(:,:,:,0) + sourcetemp
        call MPI_ISEND(source(:,:,:,0),nn,MPI_DOUBLE_PRECISION,
     &                 pe,messid,MPI_COMM_WORLD,mpirequest,mpierror)
      endif

      deallocate(sourcetemp)
 
      if (my_index == 0 .or. my_index == nslaves-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timepersource3d_slave = timepersource3d_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine setsourceforfieldsolve3d_parallel(nc,nx,ny,nzlocal,source,nxp,nyp,nzp,sourcep,
     &                                             my_index,nslaves,izpslave,nzpslave,
     &                                             izfsslave,nzfsslave)
      use Subtimers3d
      integer(ISZ):: nc,nx,ny,nzlocal
      real(kind=8):: source(0:nc-1,0:nx,0:ny,0:nzlocal)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: sourcep(0:nc-1,0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" source from
c iz=0 to either iz=nzlocal-1 or nzlocal-2 depending on the overlap. Each is
c only responsible for sending out the source is owns.

      integer(MPIISZ):: i,right_nz,onz,nn
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpierror
      integer(MPIISZ):: mpirequest(nslaves)
      integer(MPIISZ):: w
      integer(MPIISZ):: messid = 60
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent and received.
c     --- Data to be sent
      right_nz = 1
      if (my_index < nslaves-1)
     &  right_nz=izpslave(my_index)+nzpslave(my_index)-izpslave(my_index+1)+1
      do i=0,nslaves-1
        izsend(i) = max(izpslave(my_index),izfsslave(i)) - izpslave(my_index)
        nzsend(i) = min(izpslave(my_index)+nzpslave(my_index)-right_nz,
     &                 izfsslave(i)+nzfsslave(i))
     &             - max(izpslave(my_index),izfsslave(i)) + 1
c       --- For all to all
c       if (nzsend(i) <= 0) then
c         nzsend(i) = 0
c         izsend(i) = 0
c       endif
      enddo

c     --- Data to be received
      do i=0,nslaves-1
        right_nz = 1
        if (i < nslaves-1)
     &    right_nz=izpslave(i)+nzpslave(i)-izpslave(i+1)+1
        izrecv(i) = max(izfsslave(my_index),izpslave(i)) - izfsslave(my_index)
        nzrecv(i) = min(izfsslave(my_index)+nzfsslave(my_index),
     &                 izpslave(i)+nzpslave(i)-right_nz)
     &             - max(izfsslave(my_index),izpslave(i)) + 1
c       --- For all to all
c       if (nzrecv(i) <= 0) then
c         nzrecv(i) = 0
c         izrecv(i) = 0
c       endif
      enddo

c     --- For all to all
c     nzsend = nzsend*(nxp+1)*(nyp+1)
c     izsend = izsend*(nxp+1)*(nyp+1)
c     nzrecv = nzrecv*(nx +1)*(ny +1)
c     izrecv = izrecv*(nx +1)*(ny +1)

c     call MPI_ALLTOALLV(sourcep,nzsend,izsend,MPI_DOUBLE_PRECISION,
c    &                   source ,nzrecv,izrecv,MPI_DOUBLE_PRECISION,
c    &                   MPI_COMM_WORLD,mpierror)
c     --- Everything after this point should be delete for all to all

c     --- Send the data out to processors that need it.
c     --- Copy any data locally as well.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if (i /= my_index) then
            w = w + 1
            nn = nzsend(i)*(nxp+1)*(nyp+1)*nc
            call MPI_ISEND(sourcep(0,0,0,izsend(i)),nn,MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          else
            if (loc(source) .ne. loc(sourcep)) then
c             --- should be izrecv(i)+nzrecv(i)-1?
c             --- should be fixed so no unnecessary data is sent!!!
              onz = 1
              if (my_index == nslaves-1) onz = 0
              source(:,:,:,izrecv(i):izrecv(i)+nzrecv(i) - onz) =
     &          sourcep(:,:,:,izsend(i):izsend(i)+nzsend(i) - onz)
            endif
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if (i /= my_index) then
          if (nzrecv(i) > 0) then
            nn = nzrecv(i)*(nx+1)*(ny+1)*nc
            call MPI_RECV(source(0,0,0,izrecv(i)),nn,MPI_DOUBLE_PRECISION,
     &                    i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          endif
        endif
      enddo

c     --- This is needed since the sends are done without buffering. This
c     --- only matters though for the SOR field solver which modifies source
c     --- for optimization purposes.
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
c     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (lw3dtimesubs) timesetsourceforfieldsolve3d = timesetsourceforfieldsolve3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine perpot3d_slave(pot,nc,nx,ny,nzlocal,delx,dely)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nc,nx,ny,nzlocal,delx,dely
      real(kind=8):: pot(0:nc-1,-delx:nx+delx,-dely:ny+dely,-1:nzlocal+1)
c  Sets the slices on the exterior of a potential for periodicity
c  sets slice at -1 equal to the slice at nzlocal-1
c  sets slice at nzlocal+1 equal to the slice at 1
c  Only first and last processors do anything.
      integer(MPIISZ):: nn1,nn2,pe0,pens
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
      integer(MPIISZ):: messid = 70
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
 
      nn1 = 1*nc*(nx+1+2*delx)*(ny+1+2*dely)
      nn2 = 2*nc*(nx+1+2*delx)*(ny+1+2*dely)
      pe0 = 0
      pens = nslaves-1
      if (my_index == 0) then
        call MPI_ISEND(pot(:,:,:,0:1),nn2,MPI_DOUBLE_PRECISION,
     &                 pens,messid,MPI_COMM_WORLD,mpirequest,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_ISEND(pot(:,:,:,nzlocal-1),nn1,MPI_DOUBLE_PRECISION,
     &                 pe0,messid,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
 
      if (my_index == 0) then
        call MPI_RECV(pot(:,:,:,-1),nn1,MPI_DOUBLE_PRECISION,
     &                pens,messid,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_RECV(pot(:,:,:,nzlocal:nzlocal+1),nn2,MPI_DOUBLE_PRECISION,
     &                pe0,messid,MPI_COMM_WORLD,mpistatus,mpierror)
      endif

      if (my_index == 0 .or. my_index == nslaves-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timeperpot3d_slave = timeperpot3d_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine getphipforparticles3d_parallel(nc,nx,ny,nzlocal,phi,
     &                                          nxp,nyp,nzp,phip,delx,dely,delz,
     &                                          my_index,nslaves,izpslave,nzpslave,
     &                                          izfsslave,nzfsslave)
      use Subtimers3d
      integer(ISZ):: nc,nx,ny,nzlocal,delx,dely,delz
      real(kind=8):: phi(0:nc-1,-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: phip(0:nc-1,-delx:nxp+delx,-dely:nyp+dely,-delz:nzp+delz)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

c Get the phi for the extent where the particles are. This gets phi from
c neighboring processors, and at the very least gets phi in the guard planes,
c iz=-1 and +1.

      integer(MPIISZ):: i,nn
      integer(MPIISZ):: izglobal,izmaxp,izmaxfs
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpirequest(nslaves)
      integer(MPIISZ):: mpierror,w
      integer(MPIISZ):: messid = 80
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent and received.
c     --- Data to be sent
      do i=0,nslaves-1
        izglobal = max(izfsslave(my_index)-delz,izpslave(i)-delz)
        izmaxp   = izpslave(i)+nzpslave(i)+delz
        izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)+delz

        izsend(i) = izglobal - izfsslave(my_index)
        nzsend(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        izglobal = max(izpslave(my_index)-delz,izfsslave(i)-delz)
        izmaxp   = izpslave(my_index)+nzpslave(my_index)+delz
        izmaxfs  = izfsslave(i)+nzfsslave(i)+delz

        izrecv(i) = izglobal - izpslave(my_index)
        nzrecv(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if ( i /= my_index) then
            w = w + 1
            nn = nc*(1+nx+2*delx)*(1+ny+2*dely)*nzsend(i)
            call MPI_ISEND(phi(0,-delx,-dely,izsend(i)),nn,MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          else
            phip(:,:,:,izrecv(i):izrecv(i)+nzrecv(i)-1) =
     &        phi(:,:,:,izsend(i):izsend(i)+nzsend(i)-1)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then
          if (nzrecv(i) > 0) then
            nn = nc*(1+nxp+2*delx)*(1+nyp+2*dely)*nzrecv(i)
            call MPI_RECV(phip(0,-delx,-dely,izrecv(i)),nn,MPI_DOUBLE_PRECISION,
     &                    i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          endif
        endif
      enddo
 
c     --- This is needed since the sends are done without buffering.
c     --- No problems have been seem without it, but this just for robustness.
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
c     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (lw3dtimesubs) timegetphiforparticles3d = timegetphiforparticles3d + wtime() - substarttime
!$OMP END MASTER

      return
      end                                                                       
c===========================================================================
      subroutine getphiforfields3d(nx,ny,nzlocal,phi)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nzlocal+1)

c Get the phi for the full extent where the fields are. This gets phi from
c neighboring processors, and at the very least gets phi in the guard planes,
c iz=-1 and +1.

      integer(MPIISZ):: i,nn
      integer(MPIISZ):: izglobal,izmax,izmaxfs
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpirequest(nslaves)
      integer(MPIISZ):: mpierror,w
      integer(MPIISZ):: messid = 81
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent and received.
c     --- Note that special cases are needed for the first and last processors
c     --- so that the guard cells are filled properly.
c     --- Data to be sent
      do i=0,nslaves-1
        if (my_index == 0) then
          izglobal = max(izfsslave(my_index)-1,izfsslave(i)-1)
        else
          izglobal = max(izfsslave(my_index),izfsslave(i)-1)
        endif
        izmax   = izfsslave(i)+nzfsslave(i)+1
        if (my_index == nslaves-1) then
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)+1
        else
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)-1
        endif

        izsend(i) = izglobal - izfsslave(my_index)
        nzsend(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if (i == 0) then
          izglobal = max(izfsslave(my_index)-1,izfsslave(i)-1)
        else
          izglobal = max(izfsslave(my_index)-1,izfsslave(i))
        endif
        izmax   = izfsslave(my_index)+nzfsslave(my_index)+1
        if (i == nslaves-1) then
          izmaxfs  = izfsslave(i)+nzfsslave(i)+1
        else
          izmaxfs  = izfsslave(i)+nzfsslave(i)-1
        endif

        izrecv(i) = izglobal - izfsslave(my_index)
        nzrecv(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if ( i /= my_index) then
            w = w + 1
            nn = nzsend(i)*(nx+3)*(ny+3)
            call MPI_ISEND(phi(-1,-1,izsend(i)),nn,
     &                     MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then
          if (nzrecv(i) > 0) then
            nn = nzrecv(i)*(nx+3)*(ny+3)
            call MPI_RECV(phi(-1,-1,izrecv(i)),nn,
     &                    MPI_DOUBLE_PRECISION,
     &                    i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          endif
        endif
      enddo
 
c     --- This is needed since the sends are done without buffering.
c     --- No problems have been seem without it, but this just for robustness.
c     if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
c     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (lw3dtimesubs) timegetphiforfields3d = timegetphiforfields3d + wtime() - substarttime
!$OMP END MASTER

      return
      end                                                                       
c===========================================================================
