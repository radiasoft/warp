#include "top.h"
c===========================================================================
c W3DSLAVE.M, version $Revision: 1.60 $, $Date: 2011/08/27 00:35:51 $
c Slave routines associated with the W3D package.
c===========================================================================
c===========================================================================
      subroutine init_w3d_parallel
      use Subtimersw3d
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
      integer(MPIISZ):: mpierror,comm_world_mpiisz,nprocstmp,my_indextmp
      integer(ISZ):: convertindextoproc
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      comm_world_mpiisz = comm_world

c This routine should only be called once.
      if (w3dinitialized) return
      w3dinitialized = .true.

c get number of slaves and my_pe
      call MPI_INITIALIZED(initialized,mpierror)
      if (.not. initialized) call MPI_INIT(mpierror)
      if (.not. lcomm_world_initted) then
        comm_world = MPI_COMM_WORLD
        lcomm_world_initted = .true.
      endif
      comm_world_mpiisz = comm_world
      call MPI_COMM_SIZE(comm_world_mpiisz,nprocstmp,mpierror)
      call MPI_COMM_RANK(comm_world_mpiisz,my_indextmp,mpierror)
      nprocs = nprocstmp
      my_index = my_indextmp

!$OMP MASTER
      if (lw3dtimesubs) timeinit_w3d_parallel = timeinit_w3d_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
      subroutine initializedecomp_work(decomp)
      use Parallel,Only: comm_world
      use Decompositionmodule
      type(Decomposition):: decomp

      integer(MPIISZ):: mpierror,comm_world_mpiisz
      include "mpif.h"
      comm_world_mpiisz = comm_world

      call MPI_COMM_DUP(comm_world_mpiisz,decomp%mpi_comm,mpierror)

c     --- Create communicator groups for processors along each axis.
      call MPI_COMM_SPLIT(decomp%mpi_comm,
     &                    decomp%iyproc+decomp%izproc*decomp%nyprocs,
     &                    decomp%my_index,decomp%mpi_comm_x,mpierror)
      call MPI_COMM_SPLIT(decomp%mpi_comm,
     &                    decomp%izproc+decomp%ixproc*decomp%nzprocs,
     &                    decomp%my_index,decomp%mpi_comm_y,mpierror)
      call MPI_COMM_SPLIT(decomp%mpi_comm,
     &                    decomp%ixproc+decomp%iyproc*decomp%nxprocs,
     &                    decomp%my_index,decomp%mpi_comm_z,mpierror)

      return
      end
c===========================================================================
      subroutine sw_globalsum(ns,sw)
      use Subtimersw3d
      use Parallel,Only: comm_world
      integer(ISZ):: ns
      real(kind=8):: sw(ns)
c As calculated in stptl3d, sw depends on 1 over the number of particles
c in each processor. The new value of sw is 1 over the sum of the
c reciprocals of the sw of each processor, and so depends only on the
c total number of particles.
      real(kind=8):: swtemp(ns)
      integer(ISZ):: is
      integer(MPIISZ):: mpierror,comm_world_mpiisz
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      comm_world_mpiisz = comm_world
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo
      swtemp = sw
      call MPI_ALLREDUCE(swtemp,sw,int(ns,MPIISZ),MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   comm_world_mpiisz,mpierror)
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo

!$OMP MASTER
      if (lw3dtimesubs) timesw_globalsum = timesw_globalsum + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
      subroutine sumsourcepondomainboundaries(nc,nxp,nyp,nzp,
     &                                        nxguardrho,nyguardrho,nzguardrho,
     &                                        sourcep,ppdecomp)
      use Subtimersw3d
      use Decompositionmodule
      integer(ISZ):: nc,nxp,nyp,nzp,nxguardrho,nyguardrho,nzguardrho
      real(kind=8):: sourcep(0:nc-1,-nxguardrho:nxp+nxguardrho,
     &                              -nyguardrho:nyp+nyguardrho,
     &                              -nzguardrho:nzp+nzguardrho)
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
      integer(ISZ):: ix,iy,iz,ii(0:2),axis,nguards(0:2)
      integer(ISZ),allocatable:: isend(:,:,:,:),nsend(:,:,:,:)
      integer(ISZ):: my_ixpp,my_iypp,my_izpp
      integer(ISZ):: my_nxpp,my_nypp,my_nzpp
      integer(ISZ):: i1(0:2),i2(0:2)
      integer(MPIISZ):: dims,nlocal(-1:2),starts(-1:2),sizes(-1:2)
      integer(MPIISZ):: sendtype,recvtype
      integer(MPIISZ):: itemp,ineighbor,nn
      include "mpif.h"
      integer(MPIISZ):: mpi_comm
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE)
      integer(MPIISZ):: w
      integer(MPIISZ):: messid = 60
      logical(ISZ):: lsend
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      allocate(sourcepsend(0:nc-1,-nxguardrho:nxp+nxguardrho,
     &                            -nyguardrho:nyp+nyguardrho,
     &                            -nzguardrho:nzp+nzguardrho))
      allocate(sourceprecv(0:nc-1,-nxguardrho:nxp+nxguardrho,
     &                            -nyguardrho:nyp+nyguardrho,
     &                            -nzguardrho:nzp+nzguardrho))

      npx = ppdecomp%nxprocs
      npy = ppdecomp%nyprocs
      npz = ppdecomp%nzprocs
      allocate(isend(0:2,0:npx-1,0:npy-1,0:npz-1))
      allocate(nsend(0:2,0:npx-1,0:npy-1,0:npz-1))

      dims = 4
      nguards = (/nxguardrho,nyguardrho,nzguardrho/)
      nlocal = (/nc,1+nxp+2*nxguardrho,
     &              1+nyp+2*nyguardrho,
     &              1+nzp+2*nzguardrho/)
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

            ixglobal = max(my_ixpp,ppdecomp%ix(ix)) - nxguardrho
            ixmaxp   = min(ppdecomp%ix(ix) + ppdecomp%nx(ix),
     &                     my_ixpp + my_nxpp) + nxguardrho
            isend(0,ix,iy,iz) = ixglobal - my_ixpp
            nsend(0,ix,iy,iz) = max(0,ixmaxp - ixglobal + 1)

            iyglobal = max(my_iypp,ppdecomp%iy(iy)) - nyguardrho
            iymaxp   = min(ppdecomp%iy(iy) + ppdecomp%ny(iy),
     &                     my_iypp + my_nypp) + nyguardrho
            isend(1,ix,iy,iz) = iyglobal - my_iypp
            nsend(1,ix,iy,iz) = max(0,iymaxp - iyglobal + 1)

            izglobal = max(my_izpp,ppdecomp%iz(iz)) - nzguardrho
            izmaxp   = min(ppdecomp%iz(iz) + ppdecomp%nz(iz),
     &                     my_izpp + my_nzpp) + nzguardrho
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
                starts(0:2) = isend(:,ii(0),ii(1),ii(2)) + nguards
                sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
                call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                        MPI_ORDER_FORTRAN,
     &                                        MPI_DOUBLE_PRECISION,
     &                                        sendtype,mpierror)
                call MPI_TYPE_COMMIT(sendtype,mpierror)
                nn = 1
                call MPI_SEND(sourcepsend,nn,sendtype,
     &                        itemp,messid,mpi_comm,mpierror)
                call MPI_TYPE_FREE(sendtype,mpierror)
              endif
            enddo
          endif

c         --- Gather up the data sent to this processor.
          do itemp=ppdecomp%iprocgrid(axis)+ineighbor,
     &             ppdecomp%iprocgrid(axis)-ineighbor,
     &             -2*ineighbor
            if (itemp < 0 .or. itemp > ppdecomp%nprocgrid(axis)-1) cycle
            ii(axis) = itemp
            if (product(nsend(:,ii(0),ii(1),ii(2))) > 0) then
              starts(0:2) = isend(:,ii(0),ii(1),ii(2)) + nguards
              sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
              call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      recvtype,mpierror)
              call MPI_TYPE_COMMIT(recvtype,mpierror)
              nn = 1
              call MPI_RECV(sourceprecv,nn,recvtype,
     &                      itemp,messid,mpi_comm,mpistatus,mpierror)
              call MPI_TYPE_FREE(recvtype,mpierror)
              i1 = starts(0:2) - nguards
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
                starts(0:2) = isend(:,ii(0),ii(1),ii(2)) + nguards
                sizes(0:2) = nsend(:,ii(0),ii(1),ii(2))
                call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                        MPI_ORDER_FORTRAN,
     &                                        MPI_DOUBLE_PRECISION,
     &                                        sendtype,mpierror)
                call MPI_TYPE_COMMIT(sendtype,mpierror)
                nn = 1
                call MPI_SEND(sourcepsend,nn,sendtype,
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
     &                                         nxguardrho,nyguardrho,nzguardrho,
     &                                         nprocs,iproc,mpi_comm)
      use Subtimersw3d
      integer(ISZ):: axis,nc,nxlocal,nylocal,nzlocal
      integer(ISZ):: nxguardrho,nyguardrho,nzguardrho
      real(kind=8):: source(0:nc-1,-nxguardrho:nxlocal+nxguardrho,
     &                             -nyguardrho:nylocal+nyguardrho,
     &                             -nzguardrho:nzlocal+nzguardrho)
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
      integer(MPIISZ):: nguards(0:2)
      integer(MPIISZ):: sendtype
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpierror
      integer(MPIISZ):: mpi_comm_mpiisz
      integer(MPIISZ):: isend,irecv,messid = 71
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

      mpi_comm_mpiisz = mpi_comm
      dims = 4
      nlocal = (/nc,1+nxlocal+2*nxguardrho,
     &              1+nylocal+2*nyguardrho,
     &              1+nzlocal+2*nzguardrho/)
      starts = 0
      sizes = nlocal
      nguards = (/nxguardrho,nyguardrho,nzguardrho/)

      if (iproc == nprocs - 1) then
        
c       --- Select the chunk of data to send.
        starts(axis) = nlocal(axis) - 1 - 2*nguards(axis)
        sizes(axis) = 1 + 2*nguards(axis)
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtype,mpierror)
        call MPI_TYPE_COMMIT(sendtype,mpierror)
        nn = 1
        isend = 0
        call MPI_SEND(source,nn,sendtype,isend,messid,mpi_comm_mpiisz,mpierror)

c       --- receive back the data sent plus the data from the other processor
        nn = 1
        irecv = 0
        call MPI_RECV(source,nn,sendtype,irecv,messid,mpi_comm_mpiisz,mpistatus,mpierror)
        call MPI_TYPE_FREE(sendtype,mpierror)

      else if (iproc == 0) then

c       --- Select the chunk of data to exchange.
        i1 = -nguards
        i2 = (/nxlocal,nylocal,nzlocal/) + nguards
        i2(axis) = nguards(axis)
        nn = nc*(i2(0)-i1(0)+1)*(i2(1)-i1(1)+1)*(i2(2)-i1(2)+1)

c       --- allocate temporary space to receive the incoming data into
        allocate(sourcetemp(0:nc-1,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)))

        irecv = nprocs-1
        call MPI_RECV(sourcetemp,nn,MPI_DOUBLE_PRECISION,
     &                irecv,messid,mpi_comm_mpiisz,mpistatus,mpierror)

c       --- do the sum
        source(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) =
     &  source(:,i1(0):i2(0),i1(1):i2(1),i1(2):i2(2)) + sourcetemp

c       --- send the summed data back
        starts(axis) = 0
        sizes(axis) = 1 + 2*nguards(axis)
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtype,mpierror)
        call MPI_TYPE_COMMIT(sendtype,mpierror)
        nn = 1
        isend = nprocs-1
        call MPI_SEND(source,nn,sendtype,isend,messid,mpi_comm_mpiisz,mpierror)
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
     &                                         nxguardrho,nyguardrho,nzguardrho,
     &                                         fsdecomp,ppdecomp)
      use Subtimersw3d
      use Decompositionmodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal
      integer(ISZ):: nxguardrho,nyguardrho,nzguardrho
      real(kind=8):: source(0:nc-1,-nxguardrho:nxlocal+nxguardrho,
     &                             -nyguardrho:nylocal+nyguardrho,
     &                             -nzguardrho:nzlocal+nzguardrho)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: sourcep(0:nc-1,-nxguardrho:nxp+nxguardrho,
     &                              -nyguardrho:nyp+nyguardrho,
     &                              -nzguardrho:nzp+nzguardrho)
      type(Decomposition):: fsdecomp,ppdecomp

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" source from
c iz=0 to either iz=nzlocal-1 or nzlocal-2 depending on the overlap. Each is
c only responsible for sending out the source it owns.

      integer(ISZ):: npx,npy,npz
      integer(ISZ):: right_nx,right_ny,right_nz
      integer(ISZ):: ix,iy,iz
      integer(ISZ):: my_ixpp,my_iypp,my_izpp
      integer(ISZ):: my_nxpp,my_nypp,my_nzpp
      integer(ISZ):: my_ixfs,my_iyfs,my_izfs
      integer(ISZ):: my_nxfs,my_nyfs,my_nzfs
      include "mpif.h"
      integer(MPIISZ):: dims,nlocal(-1:2),nplocal(-1:2),starts(-1:2),sizes(-1:2)
      integer(MPIISZ):: nguards(0:2)
      integer(MPIISZ),allocatable:: recvtypes(:,:,:),sendtypes(:,:,:)
      integer(MPIISZ),allocatable:: sendcounts(:,:,:),recvcounts(:,:,:)
      integer(MPIISZ),allocatable:: displs(:,:,:)
      integer(MPIISZ):: mpierror,mpi_comm_mpiisz
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
      nguards = (/nxguardrho,nyguardrho,nzguardrho/)
      nlocal = (/nc,1+nxlocal+2*nxguardrho,
     &              1+nylocal+2*nyguardrho,
     &              1+nzlocal+2*nzguardrho/)
      nplocal = (/nc,1+nxp+2*nxguardrho,
     &               1+nyp+2*nyguardrho,
     &               1+nzp+2*nzguardrho/)
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
      if (ppdecomp%ixproc < ppdecomp%nxprocs-1) then
        right_nx = my_ixpp + my_nxpp - ppdecomp%ix(ppdecomp%ixproc+1) + 1
      endif
      if (ppdecomp%iyproc < ppdecomp%nyprocs-1) then
        right_ny = my_iypp + my_nypp - ppdecomp%iy(ppdecomp%iyproc+1) + 1
      endif
      if (ppdecomp%izproc < ppdecomp%nzprocs-1) then
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
              starts(0:2) = starts(0:2) + nguards
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
      do iz=0,ppdecomp%nzprocs-1
        right_nz = 0 ! was 1 XXX
        if (iz < ppdecomp%nzprocs-1) then
          right_nz=ppdecomp%iz(iz)+ppdecomp%nz(iz)-ppdecomp%iz(iz+1)+1
        endif

        do iy=0,ppdecomp%nyprocs-1
          right_ny = 0 ! was 1 XXX
          if (iy < ppdecomp%nyprocs-1) then
            right_ny=ppdecomp%iy(iy)+ppdecomp%ny(iy)-ppdecomp%iy(iy+1)+1
          endif

          do ix=0,ppdecomp%nxprocs-1
            right_nx = 0 ! was 1 XXX
            if (ix < ppdecomp%nxprocs-1) then
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
              starts(0:2) = starts(0:2) + nguards
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
      mpi_comm_mpiisz = fsdecomp%mpi_comm
      call MPI_ALLTOALLW(sourcep,sendcounts,displs,sendtypes,
     &                   source ,recvcounts,displs,recvtypes,
     &                   mpi_comm_mpiisz,mpierror)

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
      use Subtimersw3d
      use Parallel
      integer(ISZ):: nc,nx,ny,nzlocal,delx,dely
      real(kind=8):: pot(0:nc-1,-delx:nx+delx,-dely:ny+dely,-1:nzlocal+1)
c  Sets the slices on the exterior of a potential for periodicity
c  sets slice at -1 equal to the slice at nzlocal-1
c  sets slice at nzlocal+1 equal to the slice at 1
c  Only first and last processors do anything.
      integer(MPIISZ):: nn1,nn2,pe0,pens
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest
      integer(MPIISZ):: mpierror,comm_world_mpiisz
      integer(MPIISZ):: messid = 70
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      comm_world_mpiisz = comm_world
 
      nn1 = 1*nc*(nx+1+2*delx)*(ny+1+2*dely)
      nn2 = 2*nc*(nx+1+2*delx)*(ny+1+2*dely)
      pe0 = 0
      pens = nzprocs-1
      if (izproc == 0) then
        call MPI_ISEND(pot(:,:,:,0:1),nn2,MPI_DOUBLE_PRECISION,
     &                 pens,messid,comm_world_mpiisz,mpirequest,mpierror)
      else if (izproc == nzprocs-1) then
        call MPI_ISEND(pot(:,:,:,nzlocal-1),nn1,MPI_DOUBLE_PRECISION,
     &                 pe0,messid,comm_world_mpiisz,mpirequest,mpierror)
      endif
 
      if (izproc == 0) then
        call MPI_RECV(pot(:,:,:,-1),nn1,MPI_DOUBLE_PRECISION,
     &                pens,messid,comm_world_mpiisz,mpistatus,mpierror)
      else if (izproc == nzprocs-1) then
        call MPI_RECV(pot(:,:,:,nzlocal:nzlocal+1),nn2,MPI_DOUBLE_PRECISION,
     &                pe0,messid,comm_world_mpiisz,mpistatus,mpierror)
      endif

      if (izproc == 0 .or. izproc == nzprocs-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timeperpot3d_slave = timeperpot3d_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine getphipforparticles3d_parallel(nc,nxlocal,nylocal,nzlocal,
     &                                         nxguardphi,nyguardphi,nzguardphi,
     &                                         phi,
     &                                         nxp,nyp,nzp,phip,
     &                                         fsdecomp,ppdecomp)
      use Subtimersw3d
      use Decompositionmodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal
      integer(ISZ):: nxguardphi,nyguardphi,nzguardphi
      real(kind=8):: phi(0:nc-1,-nxguardphi:nxlocal+nxguardphi,
     &                          -nyguardphi:nylocal+nyguardphi,
     &                          -nzguardphi:nzlocal+nzguardphi)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: phip(0:nc-1,-nxguardphi:nxp+nxguardphi,
     &                           -nyguardphi:nyp+nyguardphi,
     &                           -nzguardphi:nzp+nzguardphi)
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
      integer(MPIISZ):: mpierror,mpi_comm_mpiisz
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
      nlocal = (/nc,1+nxlocal+2*nxguardphi,
     &              1+nylocal+2*nyguardphi,
     &              1+nzlocal+2*nzguardphi/)
      nplocal = (/nc,1+nxp+2*nxguardphi,
     &               1+nyp+2*nyguardphi,
     &               1+nzp+2*nzguardphi/)
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

            ixglobal = max(my_ixfs - nxguardphi,ppdecomp%ix(ix) - nxguardphi)
            ixmaxp   = ppdecomp%ix(ix) + ppdecomp%nx(ix) + nxguardphi
            ixmaxfs  = my_ixfs + my_nxfs + nxguardphi
            starts(0) = ixglobal - my_ixfs + nxguardphi
            sizes(0) = min(ixmaxp,ixmaxfs) - ixglobal + 1

            iyglobal = max(my_iyfs - nyguardphi,ppdecomp%iy(iy) - nyguardphi)
            iymaxp   = ppdecomp%iy(iy) + ppdecomp%ny(iy) + nyguardphi
            iymaxfs  = my_iyfs + my_nyfs + nyguardphi
            starts(1) = iyglobal - my_iyfs + nyguardphi
            sizes(1) = min(iymaxp,iymaxfs) - iyglobal + 1

            izglobal = max(my_izfs - nzguardphi,ppdecomp%iz(iz) - nzguardphi)
            izmaxp   = ppdecomp%iz(iz) + ppdecomp%nz(iz) + nzguardphi
            izmaxfs  = my_izfs + my_nzfs + nzguardphi
            starts(2) = izglobal - my_izfs + nzguardphi
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

            ixglobal = max(my_ixpp - nxguardphi,fsdecomp%ix(ix) - nxguardphi)
            ixmaxp   = my_ixpp + my_nxpp + nxguardphi
            ixmaxfs  = fsdecomp%ix(ix) + fsdecomp%nx(ix) + nxguardphi
            starts(0) = ixglobal - my_ixpp + nxguardphi
            sizes(0) = min(ixmaxp,ixmaxfs) - ixglobal + 1

            iyglobal = max(my_iypp - nyguardphi,fsdecomp%iy(iy) - nyguardphi)
            iymaxp   = my_iypp + my_nypp + nyguardphi
            iymaxfs  = fsdecomp%iy(iy) + fsdecomp%ny(iy) + nyguardphi
            starts(1) = iyglobal - my_iypp + nyguardphi
            sizes(1) = min(iymaxp,iymaxfs) - iyglobal + 1

            izglobal = max(my_izpp - nzguardphi,fsdecomp%iz(iz) - nzguardphi)
            izmaxp   = my_izpp + my_nzpp + nzguardphi
            izmaxfs  = fsdecomp%iz(iz) + fsdecomp%nz(iz) + nzguardphi
            starts(2) = izglobal - my_izpp + nzguardphi
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
      mpi_comm_mpiisz = ppdecomp%mpi_comm
      call MPI_ALLTOALLW(phi ,sendcounts,displs,sendtypes,
     &                   phip,recvcounts,displs,recvtypes,
     &                   mpi_comm_mpiisz,mpierror)

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
      subroutine getphiforfields3d(nx,ny,nzlocal,
     &                             nxguardphi,nyguardphi,nzguardphi,
     &                             phi)
      use Subtimersw3d
      use Parallel
      integer(ISZ):: nx,ny,nzlocal
      integer(ISZ):: nxguardphi,nyguardphi,nzguardphi
      real(kind=8):: phi(-nxguardphi:nx+nxguardphi,
     &                   -nyguardphi:ny+nyguardphi,
     &                   -nzguardphi:nzlocal+nzguardphi)

c THIS IS OBSOLETE, AND ONLY INCLUDES THE DECOMPOSITION IN Z.
c DO NOT USE!

c Get the phi for the full extent where the fields are. This gets phi from
c neighboring processors, and at the very least gets phi in the guard planes,
c iz=-1 and +1.

      integer(ISZ):: i
      integer(ISZ):: izglobal,izmax,izmaxfs
      integer(MPIISZ):: nn,iproc,nnxy
      integer(MPIISZ):: izsend(0:nzprocs-1),nzsend(0:nzprocs-1)
      integer(MPIISZ):: izrecv(0:nzprocs-1),nzrecv(0:nzprocs-1)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nzprocs),mpirequest(nzprocs)
      integer(MPIISZ):: mpierror,w,comm_world_mpiisz
      integer(MPIISZ):: messid = 81
      integer(ISZ):: convertindextoproc
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      comm_world_mpiisz = comm_world

c     --- First, calculate what data needs to be sent and received.
c     --- Note that special cases are needed for the first and last processors
c     --- so that the guard cells are filled properly.
c     --- Data to be sent
      do i=0,nzprocs-1
        if (izproc == 0) then
          izglobal = max(izfsslave(izproc)-1,izfsslave(i)-nzguardphi)
        else
          izglobal = max(izfsslave(izproc),izfsslave(i)-nzguardphi)
        endif
        izmax = izfsslave(i) + nzfsslave(i) + nzguardphi
        if (izproc == nzprocs-1) then
          izmaxfs  = izfsslave(izproc)+nzfsslave(izproc)+nzguardphi
        else
          izmaxfs  = izfsslave(izproc)+nzfsslave(izproc)-1
        endif

        izsend(i) = izglobal - izfsslave(izproc)
        nzsend(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nzprocs-1
        if (i == 0) then
          izglobal = max(izfsslave(izproc)-nzguardphi,izfsslave(i)-1)
        else
          izglobal = max(izfsslave(izproc)-nzguardphi,izfsslave(i))
        endif
        izmax   = izfsslave(izproc)+nzfsslave(izproc)+nzguardphi
        if (i == nzprocs-1) then
          izmaxfs  = izfsslave(i)+nzfsslave(i)+nzguardphi
        else
          izmaxfs  = izfsslave(i)+nzfsslave(i)-1
        endif

        izrecv(i) = izglobal - izfsslave(izproc)
        nzrecv(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
c     --- Only even processors first to avoid a lock-up.
      nnxy = (1 + nx + 2*nxguardphi)*(1 + ny + 2*nyguardphi)
      if (mod(izproc,2) == 0) then
        w = 0
        do i=0,nzprocs-1
          if (nzsend(i) > 0) then
            if ( i /= izproc) then
              w = w + 1
              nn = nzsend(i)*nnxy
              iproc = convertindextoproc(ixproc,iyproc,i,
     &                                   nxprocs,nyprocs,nzprocs)
              call MPI_ISEND(phi(:,:,izsend(i):izsend(i)+nzsend(i)-1),nn,
     &                       MPI_DOUBLE_PRECISION,
     &                       iproc,messid,comm_world_mpiisz,mpirequest(w),mpierror)
            endif
          endif
        enddo
      endif

c     --- Then, gather up the data sent to this processor.
      do i=0,nzprocs-1
        if ( i /= izproc) then
          if (nzrecv(i) > 0) then
            nn = nzrecv(i)*nnxy
            iproc = convertindextoproc(ixproc,iyproc,i,
     &                                 nxprocs,nyprocs,nzprocs)
            call MPI_RECV(phi(:,:,izrecv(i):izrecv(i)+nzrecv(i)-1),nn,
     &                    MPI_DOUBLE_PRECISION,
     &                    iproc,messid,comm_world_mpiisz,mpistatus,mpierror)
          endif
        endif
      enddo
 
c     --- Now the odd processors
      if (mod(izproc,2) == 1) then
        w = 0
        do i=0,nzprocs-1
          if (nzsend(i) > 0) then
            if ( i /= izproc) then
              w = w + 1
              nn = nzsend(i)*nnxy
              iproc = convertindextoproc(ixproc,iyproc,i,
     &                                   nxprocs,nyprocs,nzprocs)
              call MPI_ISEND(phi(:,:,izsend(i):izsend(i)+nzsend(i)-1),nn,
     &                       MPI_DOUBLE_PRECISION,
     &                       iproc,messid,comm_world_mpiisz,mpirequest(w),mpierror)
            endif
          endif
        enddo
      endif

c     --- This is needed since the sends are done without buffering.
c     --- No problems have been seem without it, but this just for robustness.
c     if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
c     call MPI_BARRIER(comm_world_mpiisz,mpierror)

!$OMP MASTER
      if (lw3dtimesubs) timegetphiforfields3d = timegetphiforfields3d + wtime() - substarttime
!$OMP END MASTER

      return
      end                                                                       
c===========================================================================










c=============================================================================
c=============================================================================
c=============================================================================
      module sourcetransfermessagedatatypemodule
c       --- This module holds a data type which can store all of the data
c       --- needed for transferring sourcep to source. This data is cached
c       --- since it is expensive to calculate and to setup the MPI
c       --- type definitions.
        SAVE
        type messagedatatype
          integer(ISZ),pointer:: ipp1(:),ipp2(:)
          integer(ISZ),pointer:: ifs1(:,:),ifs2(:,:)
          integer(MPIISZ),pointer:: recvtypes(:),sendtypes(:)
          integer(MPIISZ),pointer:: sendcounts(:),recvcounts(:)
          integer(MPIISZ),pointer:: displs(:)
          integer(ISZ):: nc_cache = -1
          integer(ISZ):: nxguardrho_cache
          integer(ISZ):: nyguardrho_cache
          integer(ISZ):: nzguardrho_cache
          integer(ISZ),pointer:: fsix_cache(:)
          integer(ISZ),pointer:: fsiy_cache(:)
          integer(ISZ),pointer:: fsiz_cache(:)
          integer(ISZ),pointer:: fsnx_cache(:)
          integer(ISZ),pointer:: fsny_cache(:)
          integer(ISZ),pointer:: fsnz_cache(:)
          integer(ISZ),pointer:: ppix_cache(:)
          integer(ISZ),pointer:: ppiy_cache(:)
          integer(ISZ),pointer:: ppiz_cache(:)
          integer(ISZ),pointer:: ppnx_cache(:)
          integer(ISZ),pointer:: ppny_cache(:)
          integer(ISZ),pointer:: ppnz_cache(:)
        end type messagedatatype
        integer:: icache = 0
        type(messagedatatype):: messagedatacache(200)
      contains
        subroutine getsourcetransfermessagedata(nc,
     &                                        nxguardrho,nyguardrho,nzguardrho,
     &                                        ppdecomp,fsdecomp,
     &                                        messagedata)
        use Subtimersw3d
        use Decompositionmodule
        use Parallel,Only:nprocs,my_index
        integer(ISZ):: nc
        integer(ISZ):: nxguardrho,nyguardrho,nzguardrho
        type(Decomposition):: fsdecomp,ppdecomp
        type(messagedatatype):: messagedata

c Calculate the data needed to do the transfer of sourcep to source.
c The data is cached, and if already calculated, the cached data is
c returned, giving a substantial time savings.

        integer(ISZ):: ic
        integer(ISZ):: pe
        integer(ISZ):: me_fsix,me_fsiy,me_fsiz
        integer(ISZ):: me_ppix,me_ppiy,me_ppiz
        integer(ISZ):: pe_fsix,pe_fsiy,pe_fsiz
        integer(ISZ):: pe_ppix,pe_ppiy,pe_ppiz
        integer(ISZ):: ixmin,ixmax,iymin,iymax,izmin,izmax
        include "mpif.h"
        integer(MPIISZ):: dims,nlocal(-1:2),nlocalp(-1:2)
        integer(MPIISZ):: starts(-1:2),sizes(-1:2)
        integer(MPIISZ):: mpierror
        real(kind=8):: substarttime,wtime
        if (lw3dtimesubs) substarttime = wtime()

c       --- Search through the cached data to see if this data set has
c       --- already been calculated.
c       --- Note, to turn caching off, only this loop needs to be commented out
        do ic = 1,SIZE(messagedatacache)

          if (messagedatacache(ic)%nc_cache < 0) exit

          if (messagedatacache(ic)%nc_cache == nc .and.
     &        messagedatacache(ic)%nxguardrho_cache == nxguardrho .and.
     &        messagedatacache(ic)%nyguardrho_cache == nyguardrho .and.
     &        messagedatacache(ic)%nzguardrho_cache == nzguardrho .and.
     &        size(messagedatacache(ic)%fsix_cache) == fsdecomp%nxprocs .and.
     &        size(messagedatacache(ic)%fsiy_cache) == fsdecomp%nyprocs .and.
     &        size(messagedatacache(ic)%fsiz_cache) == fsdecomp%nzprocs .and.
     &        size(messagedatacache(ic)%ppix_cache) == ppdecomp%nxprocs .and.
     &        size(messagedatacache(ic)%ppiy_cache) == ppdecomp%nyprocs .and.
     &        size(messagedatacache(ic)%ppiz_cache) == ppdecomp%nzprocs) then

            if (all(messagedatacache(ic)%fsix_cache == fsdecomp%ix) .and.
     &          all(messagedatacache(ic)%fsiy_cache == fsdecomp%iy) .and.
     &          all(messagedatacache(ic)%fsiz_cache == fsdecomp%iz) .and.
     &          all(messagedatacache(ic)%fsnx_cache == fsdecomp%nx) .and.
     &          all(messagedatacache(ic)%fsny_cache == fsdecomp%ny) .and.
     &          all(messagedatacache(ic)%fsnz_cache == fsdecomp%nz) .and.
     &          all(messagedatacache(ic)%ppix_cache == ppdecomp%ix) .and.
     &          all(messagedatacache(ic)%ppiy_cache == ppdecomp%iy) .and.
     &          all(messagedatacache(ic)%ppiz_cache == ppdecomp%iz) .and.
     &          all(messagedatacache(ic)%ppnx_cache == ppdecomp%nx) .and.
     &          all(messagedatacache(ic)%ppny_cache == ppdecomp%ny) .and.
     &          all(messagedatacache(ic)%ppnz_cache == ppdecomp%nz)) then

c             --- The cached data is the same, then there's no need
c             --- to recalculate.

c             --- Setup the array pointers.
              messagedata%ipp1 => messagedatacache(ic)%ipp1
              messagedata%ipp2 => messagedatacache(ic)%ipp2
              messagedata%ifs1 => messagedatacache(ic)%ifs1
              messagedata%ifs2 => messagedatacache(ic)%ifs2
              messagedata%sendtypes => messagedatacache(ic)%sendtypes
              messagedata%recvtypes => messagedatacache(ic)%recvtypes
              messagedata%sendcounts => messagedatacache(ic)%sendcounts
              messagedata%recvcounts => messagedatacache(ic)%recvcounts
              messagedata%displs => messagedatacache(ic)%displs

              if (lw3dtimesubs) timegetsourcetransfermessagedata =
     &               timegetsourcetransfermessagedata + wtime() - substarttime
              return
            endif
          endif

        enddo

c       --- A new data set is needed. Increment the counter and put the data
c       --- in the next available spot.
        icache = mod(icache,SIZE(messagedatacache)) + 1

c       --- If nc_cache has been set, then the previous data in the
c       --- cache needs to be freed.
        if (messagedatacache(icache)%nc_cache > 0) then
          call freesourcetransfermessagedata(icache)
        endif

        messagedatacache(icache)%nc_cache = nc
        messagedatacache(icache)%nxguardrho_cache = nxguardrho
        messagedatacache(icache)%nyguardrho_cache = nyguardrho
        messagedatacache(icache)%nzguardrho_cache = nzguardrho

        allocate(messagedatacache(icache)%fsix_cache(fsdecomp%nxprocs))
        allocate(messagedatacache(icache)%fsiy_cache(fsdecomp%nyprocs))
        allocate(messagedatacache(icache)%fsiz_cache(fsdecomp%nzprocs))
        allocate(messagedatacache(icache)%fsnx_cache(fsdecomp%nxprocs))
        allocate(messagedatacache(icache)%fsny_cache(fsdecomp%nyprocs))
        allocate(messagedatacache(icache)%fsnz_cache(fsdecomp%nzprocs))
        allocate(messagedatacache(icache)%ppix_cache(ppdecomp%nxprocs))
        allocate(messagedatacache(icache)%ppiy_cache(ppdecomp%nyprocs))
        allocate(messagedatacache(icache)%ppiz_cache(ppdecomp%nzprocs))
        allocate(messagedatacache(icache)%ppnx_cache(ppdecomp%nxprocs))
        allocate(messagedatacache(icache)%ppny_cache(ppdecomp%nyprocs))
        allocate(messagedatacache(icache)%ppnz_cache(ppdecomp%nzprocs))

        messagedatacache(icache)%fsix_cache = fsdecomp%ix
        messagedatacache(icache)%fsiy_cache = fsdecomp%iy
        messagedatacache(icache)%fsiz_cache = fsdecomp%iz
        messagedatacache(icache)%fsnx_cache = fsdecomp%nx
        messagedatacache(icache)%fsny_cache = fsdecomp%ny
        messagedatacache(icache)%fsnz_cache = fsdecomp%nz
        messagedatacache(icache)%ppix_cache = ppdecomp%ix
        messagedatacache(icache)%ppiy_cache = ppdecomp%iy
        messagedatacache(icache)%ppiz_cache = ppdecomp%iz
        messagedatacache(icache)%ppnx_cache = ppdecomp%nx
        messagedatacache(icache)%ppny_cache = ppdecomp%ny
        messagedatacache(icache)%ppnz_cache = ppdecomp%nz

        allocate(messagedatacache(icache)%ipp1(0:2))
        allocate(messagedatacache(icache)%ipp2(0:2))
        allocate(messagedatacache(icache)%ifs1(0:2,0:nprocs-1))
        allocate(messagedatacache(icache)%ifs2(0:2,0:nprocs-1))
        allocate(messagedatacache(icache)%sendtypes(0:nprocs-1))
        allocate(messagedatacache(icache)%recvtypes(0:nprocs-1))
        allocate(messagedatacache(icache)%sendcounts(0:nprocs-1))
        allocate(messagedatacache(icache)%recvcounts(0:nprocs-1))
        allocate(messagedatacache(icache)%displs(0:nprocs-1))

        messagedata%ipp1 => messagedatacache(icache)%ipp1
        messagedata%ipp2 => messagedatacache(icache)%ipp2
        messagedata%ifs1 => messagedatacache(icache)%ifs1
        messagedata%ifs2 => messagedatacache(icache)%ifs2
        messagedata%sendtypes => messagedatacache(icache)%sendtypes
        messagedata%recvtypes => messagedatacache(icache)%recvtypes
        messagedata%sendcounts => messagedatacache(icache)%sendcounts
        messagedata%recvcounts => messagedatacache(icache)%recvcounts
        messagedata%displs => messagedatacache(icache)%displs

c       --- Now, do the actual calculation.

        dims = 4

        me_fsix = fsdecomp%ixproc
        me_fsiy = fsdecomp%iyproc
        me_fsiz = fsdecomp%izproc
        nlocal = (/nc,1+fsdecomp%nx(me_fsix)+2*nxguardrho,
     &                1+fsdecomp%ny(me_fsiy)+2*nyguardrho,
     &                1+fsdecomp%nz(me_fsiz)+2*nzguardrho/)

        me_ppix = ppdecomp%ixproc
        me_ppiy = ppdecomp%iyproc
        me_ppiz = ppdecomp%izproc
        nlocalp = (/nc,1+ppdecomp%nx(me_ppix)+2*nxguardrho,
     &                 1+ppdecomp%ny(me_ppiy)+2*nyguardrho,
     &                 1+ppdecomp%nz(me_ppiz)+2*nzguardrho/)

        messagedata%displs = 0

c       --- Send and recv all of the first dimension.
        starts(-1) = 0
        sizes(-1) = nc

c       --- loop over all processors
        do pe=0,nprocs-1

          pe_ppix = mod(pe,ppdecomp%nxprocs)
          pe_ppiy = mod(pe,ppdecomp%nxprocs*ppdecomp%nyprocs)/ppdecomp%nxprocs
          pe_ppiz = pe/(ppdecomp%nxprocs*ppdecomp%nyprocs)

          pe_fsix = mod(pe,fsdecomp%nxprocs)
          pe_fsiy = mod(pe,fsdecomp%nxprocs*fsdecomp%nyprocs)/fsdecomp%nxprocs
          pe_fsiz = pe/(fsdecomp%nxprocs*fsdecomp%nyprocs)

          ixmin = max(ppdecomp%ix(me_ppix),fsdecomp%ix(pe_fsix)) - nxguardrho
          iymin = max(ppdecomp%iy(me_ppiy),fsdecomp%iy(pe_fsiy)) - nyguardrho
          izmin = max(ppdecomp%iz(me_ppiz),fsdecomp%iz(pe_fsiz)) - nzguardrho

          ixmax = min(ppdecomp%ix(me_ppix)+ppdecomp%nx(me_ppix),
     &                fsdecomp%ix(pe_fsix)+fsdecomp%nx(pe_fsix)) + nxguardrho
          iymax = min(ppdecomp%iy(me_ppiy)+ppdecomp%ny(me_ppiy),
     &                fsdecomp%iy(pe_fsiy)+fsdecomp%ny(pe_fsiy)) + nyguardrho
          izmax = min(ppdecomp%iz(me_ppiz)+ppdecomp%nz(me_ppiz),
     &                fsdecomp%iz(pe_fsiz)+fsdecomp%nz(pe_fsiz)) + nzguardrho

c         --- Note that the starts are zero based.
          starts(0) = ixmin + nxguardrho - ppdecomp%ix(me_ppix)
          starts(1) = iymin + nyguardrho - ppdecomp%iy(me_ppiy)
          starts(2) = izmin + nzguardrho - ppdecomp%iz(me_ppiz)
          sizes(0) = max(0,ixmax - ixmin + 1)
          sizes(1) = max(0,iymax - iymin + 1)
          sizes(2) = max(0,izmax - izmin + 1)

          if (ALL(sizes > 0)) then
            messagedata%sendcounts(pe) = 1
            if (pe == my_index) then
              messagedata%ipp1(0) = ixmin - ppdecomp%ix(me_ppix)
              messagedata%ipp1(1) = iymin - ppdecomp%iy(me_ppiy)
              messagedata%ipp1(2) = izmin - ppdecomp%iz(me_ppiz)
              messagedata%ipp2(0) = ixmax - ppdecomp%ix(me_ppix)
              messagedata%ipp2(1) = iymax - ppdecomp%iy(me_ppiy)
              messagedata%ipp2(2) = izmax - ppdecomp%iz(me_ppiz)
            else
              call MPI_TYPE_CREATE_SUBARRAY(dims,nlocalp,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      messagedata%sendtypes(pe),mpierror)
              call MPI_TYPE_COMMIT(messagedata%sendtypes(pe),mpierror)
            endif
          else
            messagedata%sendcounts(pe) = 0
            messagedata%sendtypes(pe) = MPI_DOUBLE_PRECISION
          endif

          ixmin = max(ppdecomp%ix(pe_ppix),fsdecomp%ix(me_fsix)) - nxguardrho
          iymin = max(ppdecomp%iy(pe_ppiy),fsdecomp%iy(me_fsiy)) - nyguardrho
          izmin = max(ppdecomp%iz(pe_ppiz),fsdecomp%iz(me_fsiz)) - nzguardrho

          ixmax = min(ppdecomp%ix(pe_ppix)+ppdecomp%nx(pe_ppix),
     &                fsdecomp%ix(me_fsix)+fsdecomp%nx(me_fsix)) + nxguardrho
          iymax = min(ppdecomp%iy(pe_ppiy)+ppdecomp%ny(pe_ppiy),
     &                fsdecomp%iy(me_fsiy)+fsdecomp%ny(me_fsiy)) + nyguardrho
          izmax = min(ppdecomp%iz(pe_ppiz)+ppdecomp%nz(pe_ppiz),
     &                fsdecomp%iz(me_fsiz)+fsdecomp%nz(me_fsiz)) + nzguardrho

          starts(0) = ixmin + nxguardrho - fsdecomp%ix(me_fsix)
          starts(1) = iymin + nyguardrho - fsdecomp%iy(me_fsiy)
          starts(2) = izmin + nzguardrho - fsdecomp%iz(me_fsiz)
          sizes(0) = max(0,ixmax - ixmin + 1)
          sizes(1) = max(0,iymax - iymin + 1)
          sizes(2) = max(0,izmax - izmin + 1)

          if (ALL(sizes > 0)) then
            messagedata%recvcounts(pe) = 1
            if (pe .ne. my_index) then
              call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                      MPI_ORDER_FORTRAN,
     &                                      MPI_DOUBLE_PRECISION,
     &                                      messagedata%recvtypes(pe),mpierror)
              call MPI_TYPE_COMMIT(messagedata%recvtypes(pe),mpierror)
            endif
            messagedata%ifs1(0,pe) = ixmin - fsdecomp%ix(me_fsix)
            messagedata%ifs1(1,pe) = iymin - fsdecomp%iy(me_fsiy)
            messagedata%ifs1(2,pe) = izmin - fsdecomp%iz(me_fsiz)
            messagedata%ifs2(0,pe) = ixmax - fsdecomp%ix(me_fsix)
            messagedata%ifs2(1,pe) = iymax - fsdecomp%iy(me_fsiy)
            messagedata%ifs2(2,pe) = izmax - fsdecomp%iz(me_fsiz)
          else
            messagedata%recvcounts(pe) = 0
            messagedata%recvtypes(pe) = MPI_DOUBLE_PRECISION
          endif

        enddo

!$OMP MASTER
        if (lw3dtimesubs) timegetsourcetransfermessagedata =
     &              timegetsourcetransfermessagedata + wtime() - substarttime
!$OMP END MASTER

        return
        end subroutine getsourcetransfermessagedata
c       ======================================================================
        subroutine freesourcetransfermessagedata(ic)
        use Decompositionmodule
        use Parallel,Only:my_index
        integer(ISZ):: ic

        integer(ISZ):: pe
        integer(MPIISZ):: mpierror

        deallocate(messagedatacache(ic)%fsix_cache)
        deallocate(messagedatacache(ic)%fsiy_cache)
        deallocate(messagedatacache(ic)%fsiz_cache)
        deallocate(messagedatacache(ic)%fsnx_cache)
        deallocate(messagedatacache(ic)%fsny_cache)
        deallocate(messagedatacache(ic)%fsnz_cache)
        deallocate(messagedatacache(ic)%ppix_cache)
        deallocate(messagedatacache(ic)%ppiy_cache)
        deallocate(messagedatacache(ic)%ppiz_cache)
        deallocate(messagedatacache(ic)%ppnx_cache)
        deallocate(messagedatacache(ic)%ppny_cache)
        deallocate(messagedatacache(ic)%ppnz_cache)

c       --- Free all of the types
        do pe=0,SIZE(messagedatacache(ic)%sendcounts)-1
          if (messagedatacache(ic)%sendcounts(pe) > 0 .and.
     &        pe .ne. my_index) then
            call MPI_TYPE_FREE(messagedatacache(ic)%sendtypes(pe),mpierror)
          endif
          if (messagedatacache(ic)%recvcounts(pe) > 0 .and.
     &        pe .ne. my_index) then
            call MPI_TYPE_FREE(messagedatacache(ic)%recvtypes(pe),mpierror)
          endif
        enddo

        deallocate(messagedatacache(ic)%ipp1)
        deallocate(messagedatacache(ic)%ipp2)
        deallocate(messagedatacache(ic)%ifs1)
        deallocate(messagedatacache(ic)%ifs2)
        deallocate(messagedatacache(ic)%sendtypes)
        deallocate(messagedatacache(ic)%recvtypes)
        deallocate(messagedatacache(ic)%recvcounts)
        deallocate(messagedatacache(ic)%sendcounts)
        deallocate(messagedatacache(ic)%displs)

        return
        end subroutine freesourcetransfermessagedata
      end module sourcetransfermessagedatatypemodule
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine transfersourceptosource3d(nc,nxp,nyp,nzp,sourcep,
     &                                     nxlocal,nylocal,nzlocal,source,
     &                                     nxguardrho,nyguardrho,nzguardrho,
     &                                     ppdecomp,fsdecomp)
      use Subtimersw3d
      use Decompositionmodule
      use Parallel,Only:nprocs,my_index,comm_world
      use sourcetransfermessagedatatypemodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal,nxp,nyp,nzp
      integer(ISZ):: nxguardrho,nyguardrho,nzguardrho
      real(kind=8):: sourcep(0:nc-1,-nxguardrho:nxp+nxguardrho,
     &                              -nyguardrho:nyp+nyguardrho,
     &                              -nzguardrho:nzp+nzguardrho)
      real(kind=8):: source(0:nc-1,-nxguardrho:nxlocal+nxguardrho,
     &                             -nyguardrho:nylocal+nyguardrho,
     &                             -nzguardrho:nzlocal+nzguardrho)
      type(Decomposition):: ppdecomp,fsdecomp

c Exchange the potential on the boundary and/or in the gaurd cells.

      type(messagedatatype):: messagedata
      save messagedata

      integer(ISZ):: pe,ip
      real(kind=8),pointer:: sourcetemp(:,:,:,:)
      integer(ISZ):: ifs1(0:2),ifs2(0:2)
      integer(ISZ):: ipp1(0:2),ipp2(0:2)

      include "mpif.h"
      integer(MPIISZ):: ipisz
      integer(MPIISZ):: messid = 62
      integer(MPIISZ):: mpierror,comm_world_mpiisz
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      comm_world_mpiisz = comm_world

c     --- Generate the MPI data types needed for the communication.
      call getsourcetransfermessagedata(nc,
     &                                  nxguardrho,nyguardrho,nzguardrho,
     &                                  ppdecomp,fsdecomp,
     &                                  messagedata)

      allocate(sourcetemp(0:nc-1,-nxguardrho:nxlocal+nxguardrho,
     &                           -nyguardrho:nylocal+nyguardrho,
     &                           -nzguardrho:nzlocal+nzguardrho))

      source = 0.

      do pe=0,nprocs-1
        ip = nprocs - 1 - my_index - pe
        ip = mod(ip+nprocs,nprocs)
        ipisz = ip
        if (ip .ne. my_index) then
          if (messagedata%sendcounts(ip) > 0 .or.
     &        messagedata%recvcounts(ip) > 0) then
            call MPI_SENDRECV(sourcep,messagedata%sendcounts(ip),
     &                        messagedata%sendtypes(ip),ipisz,messid,
     &                        sourcetemp,messagedata%recvcounts(ip),
     &                        messagedata%recvtypes(ip),ipisz,messid,
     &                        comm_world_mpiisz,mpistatus,mpierror)
            if (messagedata%recvcounts(ip) > 0) then
              ifs1 = messagedata%ifs1(:,ip)
              ifs2 = messagedata%ifs2(:,ip)
              source(:,ifs1(0):ifs2(0),ifs1(1):ifs2(1),ifs1(2):ifs2(2)) =
     &        source(:,ifs1(0):ifs2(0),ifs1(1):ifs2(1),ifs1(2):ifs2(2)) +
     &        sourcetemp(:,ifs1(0):ifs2(0),ifs1(1):ifs2(1),ifs1(2):ifs2(2))
            endif
          endif
        else
          if (messagedata%sendcounts(my_index) > 0) then
            ifs1 = messagedata%ifs1(:,my_index)
            ifs2 = messagedata%ifs2(:,my_index)
            ipp1 = messagedata%ipp1
            ipp2 = messagedata%ipp2
            source(:,ifs1(0):ifs2(0),ifs1(1):ifs2(1),ifs1(2):ifs2(2)) =
     &      source(:,ifs1(0):ifs2(0),ifs1(1):ifs2(1),ifs1(2):ifs2(2)) +
     &     sourcep(:,ipp1(0):ipp2(0),ipp1(1):ipp2(1),ipp1(2):ipp2(2))
          endif
        endif
      enddo

      deallocate(sourcetemp)

!$OMP MASTER
      if (lw3dtimesubs) timetransfersourceptosource3d =
     &                  timetransfersourceptosource3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
