#include "top.h"
c===========================================================================
c W3DSLAVE.M, version $Revision: 1.49 $, $Date: 2008/11/19 18:29:50 $
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
      integer(MPIISZ):: mpierror,nprocstmp,my_indextmp
      integer(ISZ):: convertindextoproc
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c This routine should only be called once.
      if (w3dinitialized) return
      w3dinitialized = .true.

c get number of slaves and my_pe
      call MPI_INITIALIZED(initialized,mpierror)
      if (.not. initialized) call MPI_INIT(mpierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocstmp,mpierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_indextmp,mpierror)
      nprocs = nprocstmp
      my_index = my_indextmp

!$OMP MASTER
      if (lw3dtimesubs) timeinit_w3d_parallel = timeinit_w3d_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
      subroutine initializedecomp_work(decomp)
      use Decompositionmodule
      type(Decomposition):: decomp

      integer(MPIISZ):: mpierror
      include "mpif.h"

c     --- Create communicator groups for processors along each axis.
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,
     &                    decomp%iyproc+decomp%izproc*decomp%nyprocs,
     &                    decomp%my_index,decomp%mpi_comm_x,mpierror)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,
     &                    decomp%izproc+decomp%ixproc*decomp%nzprocs,
     &                    decomp%my_index,decomp%mpi_comm_y,mpierror)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,
     &                    decomp%ixproc+decomp%iyproc*decomp%nxprocs,
     &                    decomp%my_index,decomp%mpi_comm_z,mpierror)

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
      subroutine sumsourcepondomainboundaries(nc,nxp,nyp,nzp,sourcep,ppdecomp)
      use Subtimers3d
      use Decompositionmodule
      use Parallel,only:my_index
      integer(ISZ):: nc,nxp,nyp,nzp
      real(kind=8):: sourcep(0:nc-1,0:nxp,0:nyp,0:nzp)
      type(Decomposition):: ppdecomp

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" source from
c iz=0 to either iz=nzlocal-1 or nzlocal-2 depending on the overlap. Each is
c only responsible for sending out the source is owns.

      real(kind=8),allocatable:: sourcepsend(:,:,:,:)
      real(kind=8),allocatable:: sourceprecv(:,:,:,:)
      integer(ISZ):: npx,npy,npz
      integer(ISZ):: ixglobal,ixmaxp
      integer(ISZ):: iyglobal,iymaxp
      integer(ISZ):: izglobal,izmaxp
      integer(ISZ):: ix,iy,iz,ii(0:2),axis
      integer(ISZ),allocatable:: isend(:,:,:,:),nsend(:,:,:,:)
      integer(ISZ):: my_ixpp,my_iypp,my_izpp
      integer(ISZ):: my_nxpp,my_nypp,my_nzpp
      integer(ISZ):: i1(0:2),i2(0:2)
      integer(MPIISZ):: dims,nlocal(-1:2),starts(-1:2),sizes(-1:2)
      integer(MPIISZ):: sendtype,recvtype
      integer(MPIISZ):: itemp,ineighbor
      include "mpif.h"
      integer(ISZ):: mpi_comm
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE)
      integer(MPIISZ):: w
      integer(MPIISZ):: messid = 60
      logical(ISZ):: lsend
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Note that mpi_comm is type ISZ instead of MPIISZ since the
c     --- mpi_comm_x in ppdecomp are also ISZ, but it is handled consistently
c     --- so the correct data is passed around.

      allocate(sourcepsend(0:nc-1,0:nxp,0:nyp,0:nzp))
      allocate(sourceprecv(0:nc-1,0:nxp,0:nyp,0:nzp))

      npx = ppdecomp%nxprocs
      npy = ppdecomp%nyprocs
      npz = ppdecomp%nzprocs
      allocate(isend(0:2,0:npx-1,0:npy-1,0:npz-1))
      allocate(nsend(0:2,0:npx-1,0:npy-1,0:npz-1))

      dims = 4
      nlocal = (/nc,1+nxp,1+nyp,1+nzp/)
      starts = 0
      sizes = nlocal

c     --- Calculate what data needs to be sent and received.
c     --- Data is exchanged with other processors, so the data slice sent
c     --- is the same slice that is received.
      my_ixpp = ppdecomp%ix(ppdecomp%ixproc)
      my_nxpp = ppdecomp%nx(ppdecomp%ixproc)
      my_iypp = ppdecomp%iy(ppdecomp%iyproc)
      my_nypp = ppdecomp%ny(ppdecomp%iyproc)
      my_izpp = ppdecomp%iz(ppdecomp%izproc)
      my_nzpp = ppdecomp%nz(ppdecomp%izproc)

      do iz=0,ppdecomp%nzprocs-1
        do iy=0,ppdecomp%nyprocs-1
          do ix=0,ppdecomp%nxprocs-1

            ixglobal = max(my_ixpp,ppdecomp%ix(ix))
            ixmaxp   = min(ppdecomp%ix(ix) + ppdecomp%nx(ix),
     &                     my_ixpp + my_nxpp)
            isend(0,ix,iy,iz) = ixglobal - my_ixpp
            nsend(0,ix,iy,iz) = max(0,ixmaxp - ixglobal + 1)

            iyglobal = max(my_iypp,ppdecomp%iy(iy))
            iymaxp   = min(ppdecomp%iy(iy) + ppdecomp%ny(iy),
     &                     my_iypp + my_nypp)
            isend(1,ix,iy,iz) = iyglobal - my_iypp
            nsend(1,ix,iy,iz) = max(0,iymaxp - iyglobal + 1)

            izglobal = max(my_izpp,ppdecomp%iz(iz))
            izmaxp   = min(ppdecomp%iz(iz) + ppdecomp%nz(iz),
     &                     my_izpp + my_nzpp)
            isend(2,ix,iy,iz) = izglobal - my_izpp
            nsend(2,ix,iy,iz) = max(0,izmaxp - izglobal + 1)

          enddo
        enddo
      enddo
      nsend(:,ppdecomp%ixproc,ppdecomp%iyproc,ppdecomp%izproc) = 0

c     --- Send the data out to processors that need it.
c     --- Copy any data locally as well.

      do axis=0,2
c       --- Skip the non-domain-decomposed axis.
        if (npx == 1 .and. axis == 0) cycle
        if (npy == 1 .and. axis == 1) cycle
        if (npz == 1 .and. axis == 2) cycle

        ii = ppdecomp%iprocgrid
        if (axis == 0) mpi_comm = ppdecomp%mpi_comm_x
        if (axis == 1) mpi_comm = ppdecomp%mpi_comm_y
        if (axis == 2) mpi_comm = ppdecomp%mpi_comm_z

c       --- Only copy the sourcep that will actually be sent.
c       sourcepsend = sourcep
        do itemp=0,ppdecomp%nprocgrid(axis)-1
          ii(axis) = itemp
          if (product(nsend(:,ii(0),ii(1),ii(2))) > 0) then
            starts(0:2) = isend(:,ii(0),ii(1),ii(2))
            sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
            i1 = starts(0:2)
            i2 = i1 + sizes(0:2) - 1
            sourcepsend(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) =
     &          sourcep(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2))
          endif
        enddo

c       --- The loop ordering guarantees that there are no lock ups, that is
c       --- no two processors attempt to send data to each other at the same
c       --- time.
c       --- The loop starts with nearest neighbors and fans outwards, sending
c       --- data to further and further processors. This is the ineighbor
c       --- loop index.
c       --- The inner loops exchange data with the two processors downward
c       --- and upward by ineighbor spots (skipping indices out of range).
c       --- The first send is done by the processors that have
c       --- (iproc/ineighbor) even, and the second send odd.

        do ineighbor=1,ppdecomp%nprocgrid(axis)-1

          if (mod(ppdecomp%iprocgrid(axis)/ineighbor,2) == 0) then
            do itemp=ppdecomp%iprocgrid(axis)-ineighbor,
     &               ppdecomp%iprocgrid(axis)+ineighbor,
     &               2*ineighbor
              if (itemp < 0 .or. itemp > ppdecomp%nprocgrid(axis)-1) cycle
              ii(axis) = itemp
              if (product(nsend(:,ii(0),ii(1),ii(2))) > 0) then
                starts(0:2) = isend(:,ii(0),ii(1),ii(2))
                sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
                call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                        MPI_ORDER_FORTRAN,
     &                                        MPI_DOUBLE_PRECISION,
     &                                        sendtype,mpierror)
                call MPI_TYPE_COMMIT(sendtype,mpierror)
                call MPI_SEND(sourcepsend,1,sendtype,
     &                        itemp,messid,mpi_comm,mpierror)
                call MPI_TYPE_FREE(sendtype,mpierror)
              endif
            enddo
          endif

c         --- Gather up the data sent to this processor.
          do itemp=ppdecomp%iprocgrid(axis)-ineighbor,
     &             ppdecomp%iprocgrid(axis)+ineighbor,
     &             2*ineighbor
            if (itemp < 0 .or. itemp > ppdecomp%nprocgrid(axis)-1) cycle
            ii(axis) = itemp
            if (product(nsend(:,ii(0),ii(1),ii(2))) > 0) then
              starts(0:2) = isend(:,ii(0),ii(1),ii(2))
              sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
              call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      recvtype,mpierror)
              call MPI_TYPE_COMMIT(recvtype,mpierror)
              call MPI_RECV(sourceprecv,1,recvtype,
     &                      itemp,messid,mpi_comm,mpistatus,mpierror)
              call MPI_TYPE_FREE(recvtype,mpierror)
              i1 = starts(0:2)
              i2 = i1 + sizes(0:2) - 1
              sourcep(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) =
     &              sourcep(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) +
     &          sourceprecv(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2))
            endif
          enddo

          if (mod(ppdecomp%iprocgrid(axis)/ineighbor,2) == 1) then
            do itemp=ppdecomp%iprocgrid(axis)-ineighbor,
     &               ppdecomp%iprocgrid(axis)+ineighbor,
     &               2*ineighbor
              if (itemp < 0 .or. itemp > ppdecomp%nprocgrid(axis)-1) cycle
              ii(axis) = itemp
              if (product(nsend(:,ii(0),ii(1),ii(2))) > 0) then
                starts(0:2) = isend(:,ii(0),ii(1),ii(2))
                sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
                call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                        MPI_ORDER_FORTRAN,
     &                                        MPI_DOUBLE_PRECISION,
     &                                        sendtype,mpierror)
                call MPI_TYPE_COMMIT(sendtype,mpierror)
                call MPI_SEND(sourcepsend,1,sendtype,
     &                        itemp,messid,mpi_comm,mpierror)
                call MPI_TYPE_FREE(sendtype,mpierror)
              endif
            enddo
          endif
        enddo
      enddo

      deallocate(isend)
      deallocate(nsend)
      deallocate(sourcepsend)
      deallocate(sourceprecv)

!$OMP MASTER
      if (lw3dtimesubs) timesumsourcepondomainboundaries = timesumsourcepondomainboundaries + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine makesourceperiodic_slave_work(axis,source,nc,
     &                                         nxlocal,nylocal,nzlocal,
     &                                         nprocs,iproc,mpi_comm)
      use Subtimers3d
      integer(ISZ):: axis,nc,nxlocal,nylocal,nzlocal
      real(kind=8):: source(0:nc-1,0:nxlocal,0:nylocal,0:nzlocal)
      integer(ISZ):: nprocs,iproc,mpi_comm

c  Sets the slices on the exterior of source for periodicity
c  sets slice at -1 equal to the slice at nzlocal-1
c  sets slice at nzlocal+1 equal to the slice at 1
c  Only first and last processors do anything.

      integer(ISZ):: i1(0:2),i2(0:2)
      real(kind=8),allocatable:: sourcetemp(:,:,:,:)
      integer(MPIISZ):: nn,pe
      integer(ISZ):: convertindextoproc
      include "mpif.h"
      integer(MPIISZ):: dims,nlocal(-1:2),starts(-1:2),sizes(-1:2)
      integer(MPIISZ):: sendtype
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpierror
      integer(MPIISZ):: messid = 70
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      dims = 4
      nlocal = (/nc,1+nxlocal,1+nylocal,1+nzlocal/)
      starts = 0
      sizes = nlocal

      if (iproc == nprocs - 1) then
        
c       --- Select the chunk of data to send.
        starts(axis) = nlocal(axis) - 1
        sizes(axis) = 1
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtype,mpierror)
        call MPI_TYPE_COMMIT(sendtype,mpierror)
        call MPI_SEND(source,1,sendtype,0,messid,mpi_comm,mpierror)

c       --- receive back the data sent plus the data from the other processor
        call MPI_RECV(source,1,sendtype,0,messid,mpi_comm,mpistatus,mpierror)
        call MPI_TYPE_FREE(sendtype,mpierror)

      else if (iproc == 0) then

c       --- Select the chunk of data to exchange.
        i1 = 0
        i2 = (/nxlocal,nylocal,nzlocal/)
        i2(axis) = 0
        nn = nc*(i2(0)-i1(0)+1)*(i2(1)-i1(1)+1)*(i2(2)-i1(2)+1)

c       --- allocate temporary space to receive the incoming data into
        allocate(sourcetemp(0:nc-1,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)))

        call MPI_RECV(sourcetemp,nn,MPI_DOUBLE_PRECISION,
     &                nprocs-1,messid,mpi_comm,mpistatus,mpierror)

c       --- do the sum
        source(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) =
     &  source(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) + sourcetemp

c       --- send the summed data back
        starts(axis) = 0
        sizes(axis) = 1
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtype,mpierror)
        call MPI_TYPE_COMMIT(sendtype,mpierror)
        call MPI_SEND(source,1,sendtype,nprocs-1,messid,mpi_comm,mpierror)
        call MPI_TYPE_FREE(sendtype,mpierror)

c       --- The temp array can be deallocated now that the send is complete.
        deallocate(sourcetemp)

      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timemakesourceperiodic_slave_work = timemakesourceperiodic_slave_work + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine setsourceforfieldsolve3d_parallel(nc,
     &                                         nxlocal,nylocal,nzlocal,source,
     &                                         nxp,nyp,nzp,sourcep,
     &                                         fsdecomp,ppdecomp)
      use Subtimers3d
      use Decompositionmodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal
      real(kind=8):: source(0:nc-1,0:nxlocal,0:nylocal,0:nzlocal)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: sourcep(0:nc-1,0:nxp,0:nyp,0:nzp)
      type(Decomposition):: fsdecomp,ppdecomp

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" source from
c iz=0 to either iz=nzlocal-1 or nzlocal-2 depending on the overlap. Each is
c only responsible for sending out the source is owns.

      integer(ISZ):: npx,npy,npz
      integer(ISZ):: right_nx,right_ny,right_nz
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: my_ixpp,my_iypp,my_izpp
      integer(ISZ):: my_nxpp,my_nypp,my_nzpp
      integer(ISZ):: my_ixfs,my_iyfs,my_izfs
      integer(ISZ):: my_nxfs,my_nyfs,my_nzfs
      include "mpif.h"
      integer(MPIISZ):: dims,nlocal(-1:2),nplocal(-1:2),starts(-1:2),sizes(-1:2)
      integer(MPIISZ),allocatable:: recvtypes(:,:,:),sendtypes(:,:,:)
      integer(MPIISZ),allocatable:: sendcounts(:,:,:),recvcounts(:,:,:)
      integer(MPIISZ),allocatable:: displs(:,:,:)
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      npx = fsdecomp%nxprocs
      npy = fsdecomp%nyprocs
      npz = fsdecomp%nzprocs
      allocate(sendtypes(0:npx-1,0:npy-1,0:npz-1))
      allocate(recvtypes(0:npx-1,0:npy-1,0:npz-1))
      allocate(sendcounts(0:npx-1,0:npy-1,0:npz-1))
      allocate(recvcounts(0:npx-1,0:npy-1,0:npz-1))
      allocate(displs(0:npx-1,0:npy-1,0:npz-1))

      dims = 4
      nlocal = (/nc,1+nxlocal,1+nylocal,1+nzlocal/)
      nplocal = (/nc,1+nxp,1+nyp,1+nzp/)
      displs = 0

c     --- First, calculate what data needs to be sent and received.

c     --- Send and recv all of the first dimension.
      starts(-1) = 0
      sizes(-1) = nc

c     --- Data to be sent
      my_ixpp = ppdecomp%ix(ppdecomp%ixproc)
      my_nxpp = ppdecomp%nx(ppdecomp%ixproc)
      my_iypp = ppdecomp%iy(ppdecomp%iyproc)
      my_nypp = ppdecomp%ny(ppdecomp%iyproc)
      my_izpp = ppdecomp%iz(ppdecomp%izproc)
      my_nzpp = ppdecomp%nz(ppdecomp%izproc)

      right_nx = 0 ! was 1 XXX
      right_ny = 0 ! was 1 XXX
      right_nz = 0 ! was 1 XXX
      if (fsdecomp%ixproc < fsdecomp%nxprocs-1) then
        right_nx = my_ixpp + my_nxpp - ppdecomp%ix(ppdecomp%ixproc+1) + 1
      endif
      if (fsdecomp%iyproc < fsdecomp%nyprocs-1) then
        right_ny = my_iypp + my_nypp - ppdecomp%iy(ppdecomp%iyproc+1) + 1
      endif
      if (fsdecomp%izproc < fsdecomp%nzprocs-1) then
        right_nz = my_izpp + my_nzpp - ppdecomp%iz(ppdecomp%izproc+1) + 1
      endif

      do iz=0,fsdecomp%nzprocs-1
        do iy=0,fsdecomp%nyprocs-1
          do ix=0,fsdecomp%nxprocs-1

            starts(0) = max(my_ixpp,fsdecomp%ix(ix)) - my_ixpp
            sizes(0) = min(my_ixpp+my_nxpp-right_nx,
     &                     fsdecomp%ix(ix)+fsdecomp%nx(ix))
     &                 - max(my_ixpp,fsdecomp%ix(ix)) + 1

            starts(1) = max(my_iypp,fsdecomp%iy(iy)) - my_iypp
            sizes(1) = min(my_iypp+my_nypp-right_ny,
     &                     fsdecomp%iy(iy)+fsdecomp%ny(iy))
     &                 - max(my_iypp,fsdecomp%iy(iy)) + 1

            starts(2) = max(my_izpp,fsdecomp%iz(iz)) - my_izpp
            sizes(2) = min(my_izpp+my_nzpp-right_nz,
     &                     fsdecomp%iz(iz)+fsdecomp%nz(iz))
     &                 - max(my_izpp,fsdecomp%iz(iz)) + 1

            if (ALL(sizes > 0)) then
              call MPI_TYPE_CREATE_SUBARRAY(dims,nplocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      sendtypes(ix,iy,iz),mpierror)
              call MPI_TYPE_COMMIT(sendtypes(ix,iy,iz),mpierror)
              sendcounts(ix,iy,iz) = 1
            else
              sendtypes(ix,iy,iz) = MPI_DOUBLE_PRECISION
              sendcounts(ix,iy,iz) = 0
            endif

          enddo
        enddo
      enddo

      my_ixfs = fsdecomp%ix(fsdecomp%ixproc)
      my_nxfs = fsdecomp%nx(fsdecomp%ixproc)
      my_iyfs = fsdecomp%iy(fsdecomp%iyproc)
      my_nyfs = fsdecomp%ny(fsdecomp%iyproc)
      my_izfs = fsdecomp%iz(fsdecomp%izproc)
      my_nzfs = fsdecomp%nz(fsdecomp%izproc)
      do iz=0,fsdecomp%nzprocs-1
        right_nz = 0 ! was 1 XXX
        if (iz < fsdecomp%nzprocs-1) then
          right_nz=ppdecomp%iz(iz)+ppdecomp%nz(iz)-ppdecomp%iz(iz+1)+1
        endif

        do iy=0,fsdecomp%nyprocs-1
          right_ny = 0 ! was 1 XXX
          if (iy < fsdecomp%nyprocs-1) then
            right_ny=ppdecomp%iy(iy)+ppdecomp%ny(iy)-ppdecomp%iy(iy+1)+1
          endif

          do ix=0,fsdecomp%nxprocs-1
            right_nx = 0 ! was 1 XXX
            if (ix < fsdecomp%nxprocs-1) then
              right_nx=ppdecomp%ix(ix)+ppdecomp%nx(ix)-ppdecomp%ix(ix+1)+1
            endif

            starts(0) = max(my_ixfs,ppdecomp%ix(ix)) - my_ixfs
            sizes(0) = min(my_ixfs+my_nxfs,
     &                     ppdecomp%ix(ix)+ppdecomp%nx(ix)-right_nx)
     &                - max(my_ixfs,ppdecomp%ix(ix)) + 1

            starts(1) = max(my_iyfs,ppdecomp%iy(iy)) - my_iyfs
            sizes(1) = min(my_iyfs+my_nyfs,
     &                     ppdecomp%iy(iy)+ppdecomp%ny(iy)-right_ny)
     &                - max(my_iyfs,ppdecomp%iy(iy)) + 1

            starts(2) = max(my_izfs,ppdecomp%iz(iz)) - my_izfs
            sizes(2) = min(my_izfs+my_nzfs,
     &                     ppdecomp%iz(iz)+ppdecomp%nz(iz)-right_nz)
     &                - max(my_izfs,ppdecomp%iz(iz)) + 1

            if (ALL(sizes > 0)) then
              call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      recvtypes(ix,iy,iz),mpierror)
              call MPI_TYPE_COMMIT(recvtypes(ix,iy,iz),mpierror)
              recvcounts(ix,iy,iz) = 1
            else
              recvtypes(ix,iy,iz) = MPI_DOUBLE_PRECISION
              recvcounts(ix,iy,iz) = 0
            endif

          enddo
        enddo
      enddo

c     --- Now, do all of the message passing at once.
      call MPI_ALLTOALLW(sourcep,sendcounts,displs,sendtypes,
     &                   source ,recvcounts,displs,recvtypes,
     &                   MPI_COMM_WORLD,mpierror)

c     --- Free all of the types
      do iz=0,fsdecomp%nzprocs-1
        do iy=0,fsdecomp%nyprocs-1
          do ix=0,fsdecomp%nxprocs-1
            if (sendcounts(ix,iy,iz) > 0) then
              call MPI_TYPE_FREE(sendtypes(ix,iy,iz),mpierror)
            endif
            if (recvcounts(ix,iy,iz) > 0) then
              call MPI_TYPE_FREE(recvtypes(ix,iy,iz),mpierror)
            endif
          enddo
        enddo
      enddo

      deallocate(recvtypes)
      deallocate(sendtypes)
      deallocate(sendcounts)
      deallocate(recvcounts)
      deallocate(displs)

!$OMP MASTER
      if (lw3dtimesubs) timesetsourceforfieldsolve3d_parallel = timesetsourceforfieldsolve3d_parallel + wtime() - substarttime
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
      pens = nzprocs-1
      if (my_index == 0) then
        call MPI_ISEND(pot(:,:,:,0:1),nn2,MPI_DOUBLE_PRECISION,
     &                 pens,messid,MPI_COMM_WORLD,mpirequest,mpierror)
      else if (my_index == nzprocs-1) then
        call MPI_ISEND(pot(:,:,:,nzlocal-1),nn1,MPI_DOUBLE_PRECISION,
     &                 pe0,messid,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
 
      if (my_index == 0) then
        call MPI_RECV(pot(:,:,:,-1),nn1,MPI_DOUBLE_PRECISION,
     &                pens,messid,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == nzprocs-1) then
        call MPI_RECV(pot(:,:,:,nzlocal:nzlocal+1),nn2,MPI_DOUBLE_PRECISION,
     &                pe0,messid,MPI_COMM_WORLD,mpistatus,mpierror)
      endif

      if (my_index == 0 .or. my_index == nzprocs-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timeperpot3d_slave = timeperpot3d_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine getphipforparticles3d_parallel(nc,
     &                                          nxlocal,nylocal,nzlocal,phi,
     &                                          nxp,nyp,nzp,phip,delx,dely,delz,
     &                                          fsdecomp,ppdecomp)
      use Subtimers3d
      use Decompositionmodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal
      integer(ISZ):: delx,dely,delz
      real(kind=8):: phi(0:nc-1,-delx:nxlocal+delx,-dely:nylocal+dely,-delz:nzlocal+delz)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: phip(0:nc-1,-delx:nxp+delx,-dely:nyp+dely,-delz:nzp+delz)
      type(Decomposition):: fsdecomp,ppdecomp

c Gather the potential in the region where the particles are.
c This gets phi from neighboring processors, and at the very least gets
c phi in the guard planes, iz=-1 and +1.

      integer(ISZ):: npx,npy,npz
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: my_ixpp,my_iypp,my_izpp
      integer(ISZ):: my_nxpp,my_nypp,my_nzpp
      integer(ISZ):: my_ixfs,my_iyfs,my_izfs
      integer(ISZ):: my_nxfs,my_nyfs,my_nzfs
      integer(ISZ):: ixglobal,ixmaxp,ixmaxfs
      integer(ISZ):: iyglobal,iymaxp,iymaxfs
      integer(ISZ):: izglobal,izmaxp,izmaxfs
      integer(MPIISZ):: iproc,nn
      include "mpif.h"
      integer(MPIISZ):: dims,nlocal(-1:2),nplocal(-1:2),starts(-1:2),sizes(-1:2)
      integer(MPIISZ),allocatable:: recvtypes(:,:,:),sendtypes(:,:,:)
      integer(MPIISZ),allocatable:: sendcounts(:,:,:),recvcounts(:,:,:)
      integer(MPIISZ),allocatable:: displs(:,:,:)
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      npx = fsdecomp%nxprocs
      npy = fsdecomp%nyprocs
      npz = fsdecomp%nzprocs
      allocate(sendtypes(0:npx-1,0:npy-1,0:npz-1))
      allocate(recvtypes(0:npx-1,0:npy-1,0:npz-1))
      allocate(sendcounts(0:npx-1,0:npy-1,0:npz-1))
      allocate(recvcounts(0:npx-1,0:npy-1,0:npz-1))
      allocate(displs(0:npx-1,0:npy-1,0:npz-1))

      dims = 4
      nlocal = (/nc,1+nxlocal+2*delx,1+nylocal+2*dely,1+nzlocal+2*delz/)
      nplocal = (/nc,1+nxp+2*delx,1+nyp+2*dely,1+nzp+2*delz/)
      displs = 0

c     --- First, calculate what data needs to be sent and received.

c     --- Send and recv all of the first dimension.
      starts(-1) = 0
      sizes(-1) = nc

c     --- Data to be sent
      my_ixfs = fsdecomp%ix(fsdecomp%ixproc)
      my_nxfs = fsdecomp%nx(fsdecomp%ixproc)
      my_iyfs = fsdecomp%iy(fsdecomp%iyproc)
      my_nyfs = fsdecomp%ny(fsdecomp%iyproc)
      my_izfs = fsdecomp%iz(fsdecomp%izproc)
      my_nzfs = fsdecomp%nz(fsdecomp%izproc)

      do iz=0,fsdecomp%nzprocs-1
        do iy=0,fsdecomp%nyprocs-1
          do ix=0,fsdecomp%nxprocs-1

            ixglobal = max(my_ixfs - delx,ppdecomp%ix(ix) - delx)
            ixmaxp   = ppdecomp%ix(ix) + ppdecomp%nx(ix) + delx
            ixmaxfs  = my_ixfs + my_nxfs + delx
            starts(0) = ixglobal - my_ixfs + delx
            sizes(0) = min(ixmaxp,ixmaxfs) - ixglobal + 1

            iyglobal = max(my_iyfs - dely,ppdecomp%iy(iy) - dely)
            iymaxp   = ppdecomp%iy(iy) + ppdecomp%ny(iy) + dely
            iymaxfs  = my_iyfs + my_nyfs + dely
            starts(1) = iyglobal - my_iyfs + dely
            sizes(1) = min(iymaxp,iymaxfs) - iyglobal + 1

            izglobal = max(my_izfs - delz,ppdecomp%iz(iz) - delz)
            izmaxp   = ppdecomp%iz(iz) + ppdecomp%nz(iz) + delz
            izmaxfs  = my_izfs + my_nzfs + delz
            starts(2) = izglobal - my_izfs + delz
            sizes(2) = min(izmaxp,izmaxfs) - izglobal + 1

            if (ALL(sizes > 0)) then
              call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      sendtypes(ix,iy,iz),mpierror)
              call MPI_TYPE_COMMIT(sendtypes(ix,iy,iz),mpierror)
              sendcounts(ix,iy,iz) = 1
            else
              sendtypes(ix,iy,iz) = MPI_DOUBLE_PRECISION
              sendcounts(ix,iy,iz) = 0
            endif

          enddo
        enddo
      enddo

      my_ixpp = ppdecomp%ix(ppdecomp%ixproc)
      my_nxpp = ppdecomp%nx(ppdecomp%ixproc)
      my_iypp = ppdecomp%iy(ppdecomp%iyproc)
      my_nypp = ppdecomp%ny(ppdecomp%iyproc)
      my_izpp = ppdecomp%iz(ppdecomp%izproc)
      my_nzpp = ppdecomp%nz(ppdecomp%izproc)

      do iz=0,fsdecomp%nzprocs-1
        do iy=0,fsdecomp%nyprocs-1
          do ix=0,fsdecomp%nxprocs-1

            ixglobal = max(my_ixpp - delx,fsdecomp%ix(ix) - delx)
            ixmaxp   = my_ixpp + my_nxpp + delx
            ixmaxfs  = fsdecomp%ix(ix) + fsdecomp%nx(ix) + delx
            starts(0) = ixglobal - my_ixpp + delx
            sizes(0) = min(ixmaxp,ixmaxfs) - ixglobal + 1

            iyglobal = max(my_iypp - dely,fsdecomp%iy(iy) - dely)
            iymaxp   = my_iypp + my_nypp + dely
            iymaxfs  = fsdecomp%iy(iy) + fsdecomp%ny(iy) + dely
            starts(1) = iyglobal - my_iypp + dely
            sizes(1) = min(iymaxp,iymaxfs) - iyglobal + 1

            izglobal = max(my_izpp - delz,fsdecomp%iz(iz) - delz)
            izmaxp   = my_izpp + my_nzpp + delz
            izmaxfs  = fsdecomp%iz(iz) + fsdecomp%nz(iz) + delz
            starts(2) = izglobal - my_izpp + delz
            sizes(2) = min(izmaxp,izmaxfs) - izglobal + 1

            if (ALL(sizes > 0)) then
              call MPI_TYPE_CREATE_SUBARRAY(dims,nplocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      recvtypes(ix,iy,iz),mpierror)
              call MPI_TYPE_COMMIT(recvtypes(ix,iy,iz),mpierror)
              recvcounts(ix,iy,iz) = 1
            else
              recvtypes(ix,iy,iz) = MPI_DOUBLE_PRECISION
              recvcounts(ix,iy,iz) = 0
            endif

          enddo
        enddo
      enddo

c     --- Now, do all of the message passing at once.
      call MPI_ALLTOALLW(phi ,sendcounts,displs,sendtypes,
     &                   phip,recvcounts,displs,recvtypes,
     &                   MPI_COMM_WORLD,mpierror)

c     --- Free all of the types
      do iz=0,fsdecomp%nzprocs-1
        do iy=0,fsdecomp%nyprocs-1
          do ix=0,fsdecomp%nxprocs-1
            if (sendcounts(ix,iy,iz) > 0) then
              call MPI_TYPE_FREE(sendtypes(ix,iy,iz),mpierror)
            endif
            if (recvcounts(ix,iy,iz) > 0) then
              call MPI_TYPE_FREE(recvtypes(ix,iy,iz),mpierror)
            endif
          enddo
        enddo
      enddo

      deallocate(recvtypes)
      deallocate(sendtypes)
      deallocate(sendcounts)
      deallocate(recvcounts)
      deallocate(displs)

!$OMP MASTER
      if (lw3dtimesubs) timegetphipforparticles3d_parallel = timegetphipforparticles3d_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
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
      integer(MPIISZ):: izsend(0:nzprocs-1),nzsend(0:nzprocs-1)
      integer(MPIISZ):: izrecv(0:nzprocs-1),nzrecv(0:nzprocs-1)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nzprocs),mpirequest(nzprocs)
      integer(MPIISZ):: mpierror,w
      integer(MPIISZ):: messid = 81
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent and received.
c     --- Note that special cases are needed for the first and last processors
c     --- so that the guard cells are filled properly.
c     --- Data to be sent
      do i=0,nzprocs-1
        if (my_index == 0) then
          izglobal = max(izfsslave(my_index)-1,izfsslave(i)-1)
        else
          izglobal = max(izfsslave(my_index),izfsslave(i)-1)
        endif
        izmax   = izfsslave(i)+nzfsslave(i)+1
        if (my_index == nzprocs-1) then
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)+1
        else
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)-1
        endif

        izsend(i) = izglobal - izfsslave(my_index)
        nzsend(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nzprocs-1
        if (i == 0) then
          izglobal = max(izfsslave(my_index)-1,izfsslave(i)-1)
        else
          izglobal = max(izfsslave(my_index)-1,izfsslave(i))
        endif
        izmax   = izfsslave(my_index)+nzfsslave(my_index)+1
        if (i == nzprocs-1) then
          izmaxfs  = izfsslave(i)+nzfsslave(i)+1
        else
          izmaxfs  = izfsslave(i)+nzfsslave(i)-1
        endif

        izrecv(i) = izglobal - izfsslave(my_index)
        nzrecv(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nzprocs-1
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
      do i=0,nzprocs-1
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
