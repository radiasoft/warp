#include "top.h"
c===========================================================================
c W3DSLAVE.M, version $Revision: 1.37 $, $Date: 2006/11/03 13:46:24 $
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

      integer(ISZ):: i,is,overlap,bestnz
      integer(ISZ):: nzinj
      real(kind=8):: zinjmax(ninject)
      real(kind=8):: zperproc,avezpp,ztot,sumzslave,zlast
      real(kind=8):: zpmin,zpmax
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

c Calculate how the work is to be arranged among the processors.

c---------------------------------------------------------------------------
c     --- An overlap of one plane is needed for the FFT field solvers.
c     --- An overlap of two planes is needed for SOR field solver since
c     --- plane 0 of each processor must overlap with a plane which is
c     --- calculated in the neighboring processor.
c     --- As an example, for nz=8 divided among 2 PE's...
c     ---   for overlap=1        0 1 2 3 4 5 6 7 8
c     ---                       |_________|
c     ---                               |_________|
c     ---   for overlap=2        0 1 2 3 4 5 6 7 8
c     ---                       |___________|
c     ---                               |_________|
      if (fstype == 3 .or. fstype == 7 .or. fstype == 8 .or. fstype == 13) then
        overlap = 2
      else
        overlap = 1
      endif
      grid_overlap = overlap

c---------------------------------------------------------------------------
c     --- Do the domain decomposition for the field solver.
c     --- The domain decomposition by default is done so that each processor
c     --- has nearly the same number of z planes. If lautodecomp is false and
c     --- fstype == 3, then the decomposition is supplied by the user (input
c     --- through nzfsslave).
      if (lfsautodecomp .or.
     &    (fstype /= 3 .and. fstype /= 7 .and. fstype /= 13)) then

c       --- Calculate average number of z planes per processor, including the
c       --- extra space for overlap.  For FFT, the minimum number of planes
c       --- allowable is 1 and for PSOR, the minimum number of planes allowable
c       --- is 2.
        zperproc = (nz + (nslaves - 1.)*(overlap - 1))/nslaves
        if (zperproc < overlap) zperproc = overlap

c       --- bestnz is zperproc rounded down.
        bestnz = int(zperproc)
     
c       --- First processor is easy
        nzfsslave(0) = bestnz
        izfsslave(0) = 0
        ztot = nzfsslave(0)

c       --- loop over processors until used all processors or have assigned
c       --- all of grid.
        i = 0
        do while (i < nslaves-1 .and.
     &            izfsslave(i)+nzfsslave(i)-overlap+1 < nz)
          i = i + 1

c         --- This processor starts at the end of the region covered
c         --- by the previous processor, overlapping it by the value of overlap.
          izfsslave(i) = izfsslave(i-1) + nzfsslave(i-1) + 1 - overlap

c         --- The number of z planes given to this processor is first assumed to
c         --- be bestnz.  If this gives an average number of z planes per
c         --- processor that is less then zperproc, it is increased by 1.
          nzfsslave(i) = bestnz
          ztot = ztot + nzfsslave(i)
          avezpp = ztot/(i+1)
          if (avezpp < zperproc) then
            nzfsslave(i) = nzfsslave(i) + 1
            ztot = ztot + 1
          endif

c         --- Check if this region extends past the end of the grid.  If so,
c         --- recalculate nzfsslave(i).
          if (izfsslave(i) + nzfsslave(i) > nz) then
            nzfsslave(i) = nz - izfsslave(i)
c           --- If nzfsslave(i) is less than 3 then skip this
c           --- processor and give remaining zones to previous processor.
c           if (nzfsslave(i) < 3) then
c             i = i - 1
c             nzfsslave(i) = nz - izfsslave(i)
c           endif
          endif

        enddo

c       --- Save the number of processors that have part of the grid assigned
c       --- to them.
        nslaves = i+1

      else

c       --- nzfsslave is assumed to be input by the user and is assumed to
c       --- not include the overlap.

c       --- First check to make sure that all values are > 0.
        do i=0,nslaves-1
          if (nzfsslave(i) == 0) then
            call remark("ERROR: nz for the field solver for each processor")
            call remark("       must be greater than zero")
            call kaboom(1)
          endif
        enddo

c       --- Fill in the array izfsslave, based on the inputted nzfsslave and
c       --- add the overlap to nzfsslave.
        izfsslave(0) = 0
        do i=1,nslaves-1
          izfsslave(i) = izfsslave(i-1) + nzfsslave(i-1)
        enddo
c       --- Note that the last processor has no overlap
        do i=0,nslaves-2
          nzfsslave(i) = nzfsslave(i) + overlap - 1
        enddo
c       --- Get the new value of nz.
        nz = izfsslave(nslaves-1) + nzfsslave(nslaves-1)

      endif

c---------------------------------------------------------------------------
c     --- Set the grid cell size, which is needed below. It must be set here
c     --- after the field solve domain decomposition since nz may be
c     --- determined from the user supplied decomposition.
      dz = (zmmax - zmmin)/nz

c---------------------------------------------------------------------------
c     --- Now, set the domain decomposition for the particles. This can
c     --- either be input from the user or set the same as the decompostion
c     --- for the field solver.
      if (lautodecomp) then

c       --- Set from field solver decompostion
        do i=0,nslaves-1
          izpslave(i) = izfsslave(i)
          nzpslave(i) = nzfsslave(i) - (overlap - 1)
c         --- These shouldn't include the guard cells
          zpslmin(i) = izpslave(i)*dz + zmmin
          zpslmax(i) = (izpslave(i)+nzpslave(i))*dz + zmmin
c         --- Now include the guard cells
          izpslave(i) = izpslave(i) - nzpguard
          nzpslave(i) = nzpslave(i) + 2*nzpguard
c         --- Make sure that the processors doesn't have grid cells
c         --- sticking out the ends.
          if (izpslave(i) < 0) then
            nzpslave(i) = nzpslave(i) + izpslave(i)
            izpslave(i) = 0
          endif
          if (izpslave(i) + nzpslave(i) > nz) nzpslave(i) = nz - izpslave(i)
        enddo
        zpslmax(nslaves-1) = zmmax

      else
c       --- It is assumed that the user supplied decomposition is specified
c       --- in the array zslave, which is the fractional z-ranges
c       --- of the particles for each processor. It is assumed that zslave
c       --- has been properly normalized so that if the sum of the z-ranges
c       --- covers the whole grid, then the sum is one.

c       --- Get sum of zslave to allow proper scaling.
        sumzslave = 1.
c       sumzslave = zslave(0)
c       do i=1,nslaves-1
c         sumzslave = sumzslave + zslave(i)
c       enddo

c       --- All values of zslave must be > 0.
        do i=0,nslaves-1
c         if (zslave(i)/sumzslave*(zmmax-zmmin) <= 0.) then
          if (zslave(i) <= 0.) then
            call remark("ERROR: The length of all particle domains must be")
            call remark("       positive. Fix zslave appropriately.")
            call kaboom(1)
          endif
        enddo

c       --- Set minimum z of each processor.
        zlast = zmmin
        do i=0,nslaves-1
          zpslmin(i) = zlast
          zpslmax(i) = zlast + zslave(i)/sumzslave*(zmmax - zmmin)
          zlast = zpslmax(i)
        enddo

c       --- When using solvergeom==AMRgeom, the particle decomposition must
c       --- be aligned with the grid.
        if (solvergeom==AMRgeom) then
         do i=0,nslaves-1
           zpslmin(i) = nint((zpslmin(i) - zmmin)/dz)*dz
           zpslmax(i) = nint((zpslmax(i) - zmmin)/dz)*dz
         enddo
        endif

c       --- This is only needed to avoid problems from round off in the
c       --- accumulation. From the loop above, zpslmax(nslaves-1) will
c       --- not be exactly the same as zmmax due to roundoff.
        if (zpslmax(nslaves-1) > zmmax) zpslmax(nslaves-1) = zmmax

c       --- Set iz and nz. This is done so that zmesh(izpslave) < zpslmin, and
c       --- zmesh(izpslave+nzpslave) > zpslmax.
        do i=0,nslaves-1
          izpslave(i) = int((zpslmin(i) - zmmin)/dz) - nzpguard
          nzpslave(i) = int((zpslmax(i) - zmmin)/dz) - izpslave(i) + 1 +
     &                  2*nzpguard
c         --- Make sure that the processors doesn't have grid cells
c         --- sticking out the end.
          if (izpslave(i) < 0) then
            nzpslave(i) = nzpslave(i) + izpslave(i)
            izpslave(i) = 0
          endif
          if (izpslave(i) + nzpslave(i) > nz) nzpslave(i) = nz - izpslave(i)
        enddo

      endif

c---------------------------------------------------------------------------
c     --- Now set the axial extent of each slaves domain, which includes
c     --- both the particle and field solve domain.
      do i=0,nslaves-1
        izslave(i) = izfsslave(i)
        nzslave(i) = nzfsslave(i)
        zmslmin(i) = izfsslave(i)*dz + zmmin
        zmslmax(i) = (izfsslave(i) + nzfsslave(i))*dz + zmmin
      enddo

c---------------------------------------------------------------------------
c     --- Save global array sizes
      nzfull = nz
      zplmin = zmmin
      zplmax = zmmax

c---------------------------------------------------------------------------
c     --- Reset local values
      if (injctspc > 0) then
        injctspc = max(1,int(injctspc/nslaves*1.1))
      endif
      nz = nzfsslave(my_index)
      nzp = nzpslave(my_index)
      zmminp = zmminglobal + izpslave(my_index)*dz
      zmmaxp = zmminp + nzp*dz
      zpmin = zmmin + izpslave(my_index)*dz
      zpmax = (izpslave(my_index)+nzpslave(my_index))*dz + zmmin
      izfsmin = 0
      izfsmax = nzfsslave(my_index)
      zmmin = zmslmin(my_index)
      zmmax = zmslmax(my_index)

c     --- Set value of izextra, which is the amount of extra space at the
c     --- end of phi. For the case when (nx+ny) < nslaves, phi needs
c     --- to be made larger to hold the transposed array. izextra could
c     --- actualy be one less then the value calculated here, but I don't
c     --- subtract one so there is some extra slop.
      if(solvergeom==XYZgeom .and. 0 <= fstype .and. fstype < 7) then
        if (nx+ny < nslaves) then
          izextra = max(1,(nzfull-nzfull/nslaves*(nx+ny))/(nx*ny))
        endif
      endif

!$OMP MASTER
      if (lw3dtimesubs) timeinit_w3d_parallel = timeinit_w3d_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
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
      subroutine sumrhoondomainboundariesold(rhop,nxp,nyp,nzp)
      use Subtimers3d
      use GlobalVars
      use Parallel
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)

c This routine sums rhop in the overlapping boundaries with the neighboring
c processors. Note that this does not send the summed rhop back to the left
c since that is done in getrhoforfieldsolve3d.
c This allows for cases where more than two processors have deposited rhop
c into the same plane. A check is made if the data being sent to the right
c overlaps with the data being received from the left. If so, the process
c waits to receive the data from the left and then adds it to its own data.
c That data is then sent to the right. This way, the data is accumulated as it
c is passed from process to process moving toward the right. If there is no
c overlap, then the data is sent first to the right and then received from
c the left. This is done so that processes to the right do not have to wait
c for processes to the left. Note that process 0 always sends it data first
c so that there won't be a lock up.
c Since no two particles domains ever overlap, at most only two z planes
c will ever have to be exchanged. The number of grid cells exchanged,
c left_nz and right_nz, are limited to 2 or less. This allows the value
c of nzpslave to overlap further than the 2 cells - currently needed only
c for special cases of injection. Nonetheless, having the limits makes the
c more robust to future changes.

      integer(MPIISZ):: left_pe,right_pe
      integer(MPIISZ):: left_iz,right_iz,left_nz,right_nz,nn
      integer(MPIISZ):: ix,iy,iz
      real(kind=8),allocatable:: buffer(:,:,:)
      include "mpif.h"
      integer(MPIISZ):: messid = 50
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpierror,mpirequest
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Location where rhop planes are to be exchanged, and number of
c     --- planes to exchange.
      if (my_index > 0) then
        left_iz = 0
        left_nz=izpslave(my_index-1)+nzpslave(my_index-1)-izpslave(my_index)+1
      else
        left_iz = 0
        left_nz = 1
      endif
      if (my_index < nslaves-1) then
        right_iz = izpslave(my_index+1) - izpslave(my_index)
        right_nz = izpslave(my_index)+nzpslave(my_index)-izpslave(my_index+1)+1
      else
        right_iz = nzp
        right_nz = 1
      endif

c     --- Check if data being sent overlaps with data being received.
c     --- If so, wait for the incoming data first.
      if (my_index > 0 .and. right_iz <= left_iz+left_nz) then
        nn = left_nz*(nxp+1)*(nyp+1)
        allocate(buffer(0:nxp,0:nyp,left_nz))
        call MPI_RECV(buffer,nn,
     &                MPI_DOUBLE_PRECISION,
     &                left_pe,messid,MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Sum the two planes. The utility call sumarry is used since
c       --- with the gchange above, the local dimension of buffer may not
c       --- be correct, so explicitly writing out the loops may not work.
        call sumarry(buffer,rhop(0,0,left_iz),left_nz*(nxp+1)*(nyp+1))
        deallocate(buffer)

      endif

c     --- Send rhop to the process to the right
      if (my_index < nslaves-1) then
        nn = right_nz*(nxp+1)*(nyp+1)
        call MPI_ISEND(rhop(0,0,right_iz),nn,
     &                 MPI_DOUBLE_PRECISION,
     &                 right_pe,messid,MPI_COMM_WORLD,mpirequest,mpierror)
      endif

c     --- If there is no overlap, wait for the data after sending so
c     --- that processors to the right do not have to also wait.
      if (my_index > 0 .and. right_iz > left_iz+left_nz) then
        nn = left_nz*(nxp+1)*(nyp+1)
        allocate(buffer(0:nxp,0:nyp,left_nz))
        call MPI_RECV(buffer,nn,
     &                MPI_DOUBLE_PRECISION,
     &                left_pe,messid,MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Sum the two planes. The utility call sumarry is used since
c       --- with the gchange above, the local dimension of buffer may not
c       --- be correct, so explicitly writing out the loops may not work.
        call sumarry(buffer,rhop(0,0,left_iz),left_nz*(nxp+1)*(nyp+1))
        deallocate(buffer)

      endif

c     --- Now, wait just to make sure that the send was recv'ed
      if (my_index < nslaves-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif

!$OMP MASTER
      if (lw3dtimesubs) timesumrhoondomainboundaries = timesumrhoondomainboundaries + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sumrhoondomainboundariestry1(rhop,nxp,nyp,nzp)
      use Subtimers3d
      use GlobalVars
      use Parallel
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)

c This routine sums rhop in the overlapping boundaries with the neighboring
c processors. Note that this does not send the summed rhop back to the left
c since that is done in getrhoforfieldsolve3d.
c This allows for cases where more than two processors have deposited rhop
c into the same plane. A check is made if the data being sent to the right
c overlaps with the data being received from the left. If so, the process
c waits to receive the data from the left and then adds it to its own data.
c That data is then sent to the right. This way, the data is accumulated as it
c is passed from process to process moving toward the right. If there is no
c overlap, then the data is sent first to the right and then received from
c the left. This is done so that processes to the right do not have to wait
c for processes to the left. Note that process 0 always sends it data first
c so that there won't be a lock up.
c Since no two particles domains ever overlap, at most only two z planes
c will ever have to be exchanged. The number of grid cells exchanged,
c left_nz and right_nz, are limited to 2 or less. This allows the value
c of nzpslave to overlap further than the 2 cells - currently needed only
c for special cases of injection. Nonetheless, having the limits makes the
c more robust to future changes.

      integer(MPIISZ):: left_pe,right_pe
      integer(MPIISZ):: left_iz,right_iz,left_nz,right_nz
      integer(MPIISZ):: nzpfull,nn,nz,izl,nzl
      integer(MPIISZ),allocatable:: izoverlap(:)
      integer(MPIISZ):: iz_comm
      integer(MPIISZ):: iz,ip,iz1,iz2,iw
      integer(MPIISZ):: color
      real(kind=8),allocatable:: buffer(:,:,:),buffer2d(:,:)
      include "mpif.h"
      integer(MPIISZ):: messid = 50
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,2),mpierror,mpirequest(2)
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
c     print*,"sumrhoondomainboundaries ",my_index
c     if (my_index == 0) print*,'iz ',izpslave
c     if (my_index == 0) print*,'nz ',nzpslave

c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Find how many processors overlap each grid cell. Where there are
c     --- three or more, special handling is needed.
      nzpfull = izpslave(nslaves-1)+nzpslave(nslaves-1)
      allocate(izoverlap(0:nzpfull))
      izoverlap = 0
      do ip=0,nslaves-1
        iz1 = izpslave(ip)
        iz2 = izpslave(ip) + nzpslave(ip)
        izoverlap(iz1:iz2) = izoverlap(iz1:iz2) + 1
      enddo
c     if (my_index == 0) print*,'ov ',nzpfull,izoverlap

c     --- Location where rhop planes are to be exchanged, and number of
c     --- planes to exchange.
      if (my_index > 0) then
        left_iz = 0
        left_nz=izpslave(my_index-1)+nzpslave(my_index-1)-izpslave(my_index)+1
      endif
      if (my_index < nslaves-1) then
        right_iz = izpslave(my_index+1) - izpslave(my_index)
        right_nz = izpslave(my_index)+nzpslave(my_index)-izpslave(my_index+1)+1
      endif
      iw = 0

c     --- Send rhop to the process to the right if the planes are only
c     --- shared by the two processors.
c     --- Figure out whether one or both planes can be sent (since one
c     --- plane may have an overlap > 2 and the other doesn't).
      if (my_index < nslaves-1) then
        if (izoverlap(izpslave(my_index)+right_iz) == 2) then
          nz = 1
          iz = 0
          if (right_nz == 2 .and.
     &        izoverlap(izpslave(my_index)+right_iz+1) == 2) then
             nz = 2
          endif
        else if (right_nz == 2 .and.
     &           izoverlap(izpslave(my_index)+right_iz+1) == 2) then
          nz = 1
          iz = 1
        else
          nz = 0
        endif
        if (nz > 0) then
          nn = nz*(nxp+1)*(nyp+1)
          iw = iw + 1
          call MPI_ISEND(rhop(:,:,right_iz+iz),nn,MPI_DOUBLE_PRECISION,
     &                   right_pe,messid,MPI_COMM_WORLD,mpirequest(iw),mpierror)
c         print*,"isend ",my_index,izpslave(my_index)+right_iz+iz,nz
        endif
      endif

c     --- Recieve rhop from the process to the left if the planes are only
c     --- shared by the two processors.
      nzl = 0
      if (my_index > 0) then
        if (izoverlap(izpslave(my_index)+left_iz) == 2) then
          nzl = 1
          izl = left_iz
          if (left_nz == 2 .and.
     &        izoverlap(izpslave(my_index)+left_iz+1) == 2) then
             nzl = 2
          endif
        else if (left_nz == 2 .and.
     &           izoverlap(izpslave(my_index)+left_iz+1) == 2) then
          nzl = 1
          izl = left_iz + 1
        else
          nzl = 0
        endif
        if (nzl > 0) then
          nn = nzl*(nxp+1)*(nyp+1)
          iw = iw + 1
          allocate(buffer(0:nxp,0:nyp,nzl))
          buffer = 3.
          call MPI_IRECV(buffer,nn,MPI_DOUBLE_PRECISION,
     &                   left_pe,messid,MPI_COMM_WORLD,mpirequest(iw),mpierror)
c         --- Note - don't change izl or nzl since they are used below.
c         --- The buffer is added into rhop at the end after the WAIT.
c       print*,"irecv ",my_index,izpslave(my_index)+left_iz+izl,nzl,left_nz

        endif
      endif

c     print*,"sumrho recv done ",my_index,maxval(rhop),maxval(izoverlap)
      if (maxval(izoverlap) >= 3) then
c       --- Now do REDUCEs on planes which have overlaps of three or more
c       --- processors.
        nn = (1+nxp)*(1+nyp)
        allocate(buffer2d(0:nxp,0:nyp))
        do ip=0,1
          do iz=ip,nzpfull,2
            if (izoverlap(iz) >= 3) then
c             --- First setup the communicator
              if (izpslave(my_index) <= iz .and.
     &            iz <= izpslave(my_index)+nzpslave(my_index)) then
                color = 1
              else
                color = MPI_UNDEFINED
              endif
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,int(0,MPIISZ),iz_comm,mpierror)
              if (color == 1) then
                call MPI_ALLREDUCE(rhop(:,:,iz-izpslave(my_index)),buffer2d,nn,
     &                             MPI_DOUBLE_PRECISION,MPI_SUM,iz_comm,
     &                             mpierror)
                call MPI_COMM_FREE(iz_comm,mpierror)
                rhop(:,:,iz-izpslave(my_index)) = buffer2d
              endif
            endif
          enddo
        enddo
        deallocate(buffer2d)
      endif

      deallocate(izoverlap)

c     print*,"sumrho reduce done ",my_index,maxval(rhop),iw
c     --- Now, wait for any outstanding messages to clear.
      if (iw > 0) then
        call MPI_WAITALL(iw,mpirequest,mpistatus,mpierror)
      endif

c     --- Now, add in rhop from the left.
      if (nzl > 0) then
c       print*,"buffer ",my_index,maxval(buffer)
        rhop(:,:,izl:izl+nzl-1) = rhop(:,:,izl:izl+nzl-1) + buffer
        deallocate(buffer)
      endif
c     print*,"sumrho all done ",my_index,maxval(rhop)

!$OMP MASTER
      if (lw3dtimesubs) timesumrhoondomainboundaries = timesumrhoondomainboundaries + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sumrhoondomainboundariestry2(rhop,nxp,nyp,nzp)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)

      integer(ISZ):: allocerror
      integer(MPIISZ):: i,nn,i1,i2,iz,nzbuff
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      real(kind=8),allocatable:: sendbuffer(:,:,:),recvbuffer(:,:,:)
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
c     --- two planes of rhop.
      nzbuff = min(nzp+1,4)
      allocate(sendbuffer(0:nxp,0:nyp,nzbuff))
      sendbuffer(:,:,1:2) = rhop(:,:,0:1)
      if (nzbuff == 3) then
        sendbuffer(:,:,3) = rhop(:,:,nzp)
      else if (nzbuff >= 4) then
        sendbuffer(:,:,3:4) = rhop(:,:,nzp-1:nzp)
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
          nn = nzsend(i)*(nxp+1)*(nyp+1)
          call MPI_ISEND(sendbuffer(:,:,iz),nn,MPI_DOUBLE_PRECISION,
     &                    i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
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
      allocate(recvbuffer(0:nxp,0:nyp,maxval(nzrecv)),stat=allocerror)
      do i=0,nslaves-1
        if (i == my_index) cycle
        if (nzrecv(i) > 0) then
          nn = nzrecv(i)*(nxp+1)*(nyp+1)
          call MPI_RECV(recvbuffer,nn,MPI_DOUBLE_PRECISION,
     &                  i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          i1 = izrecv(i)
          i2 = izrecv(i) + nzrecv(i) - 1
          rhop(:,:,i1:i2) = rhop(:,:,i1:i2) + recvbuffer(:,:,1:nzrecv(i))
        endif
      enddo
      deallocate(recvbuffer)

c     --- Wait for all sends to complete before deleting the sendbuffer.
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
      deallocate(sendbuffer)

!$OMP MASTER
      if (lw3dtimesubs) timesumrhoondomainboundaries = timesumrhoondomainboundaries + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine sumrhoondomainboundaries(rhop,nxp,nyp,nzp,nzpguard)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nxp,nyp,nzp,nzpguard
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)

      integer(ISZ):: allocerror
      integer(MPIISZ):: i,nn,i1,i2,iz,nzbuff
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      real(kind=8),allocatable:: sendbuffer(:,:,:),recvbuffer(:,:,:)
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
c     --- two planes of rhop.
      nzbuff = min(nzp+1,4+4*nzpguard)
      allocate(sendbuffer(0:nxp,0:nyp,nzbuff))
      if (nzbuff == nzp+1) then
        sendbuffer = rhop
      else
        nn = min(1+2*nzpguard,nzp)
        sendbuffer(:,:,1:nn+1) = rhop(:,:,0:nn)
        nn = max(nzp-1-2*nzpguard,nn+1)
        sendbuffer(:,:,nzbuff-nzp+nn:nzbuff) = rhop(:,:,nn:nzp)
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
          nn = nzsend(i)*(nxp+1)*(nyp+1)
          call MPI_ISEND(sendbuffer(:,:,iz),nn,MPI_DOUBLE_PRECISION,
     &                    i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
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
      allocate(recvbuffer(0:nxp,0:nyp,maxval(nzrecv)),stat=allocerror)
      do i=0,nslaves-1
        if (i == my_index) cycle
        if (nzrecv(i) > 0) then
          nn = nzrecv(i)*(nxp+1)*(nyp+1)
          call MPI_RECV(recvbuffer,nn,MPI_DOUBLE_PRECISION,
     &                  i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          i1 = izrecv(i)
          i2 = izrecv(i) + nzrecv(i) - 1
          rhop(:,:,i1:i2) = rhop(:,:,i1:i2) + recvbuffer(:,:,1:nzrecv(i))
        endif
      enddo
      deallocate(recvbuffer)

c     --- Wait for all sends to complete before deleting the sendbuffer.
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
      deallocate(sendbuffer)

!$OMP MASTER
      if (lw3dtimesubs) timesumrhoondomainboundaries = timesumrhoondomainboundaries + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine perrho3d_slave(rho,nx,ny,nz)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
c  Sets the slices on the exterior of rho for periodicity
c  sets slice at -1 equal to the slice at nz-1
c  sets slice at nz+1 equal to the slice at 1
c  Only first and last processors do anything.
      real(kind=8):: rhotemp(0:nx,0:ny)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
      integer(MPIISZ):: messid = 70
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
 
      if (my_index == nslaves-1) then
        call MPI_ISEND(rho(0,0,nz),int((nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &                 int(0,MPIISZ),messid,MPI_COMM_WORLD,mpirequest,mpierror)
        call MPI_RECV(rho(0,0,nz),int((nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &                int(0,MPIISZ),messid,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == 0) then
        call MPI_RECV(rhotemp,int((nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &            int(nslaves-1,MPIISZ),messid,MPI_COMM_WORLD,mpistatus,mpierror)
        rho(:,:,0) = rho(:,:,0) + rhotemp
        call MPI_ISEND(rho(0,0,0),int((nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &             int(nslaves-1,MPIISZ),messid,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
 
      if (my_index == 0 .or. my_index == nslaves-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timeperrho3d_slave = timeperrho3d_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine setrhoforfieldsolve3d_parallel(nx,ny,nz,rho,nxp,nyp,nzp,rhop)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" rho from
c iz=0 to either iz=nz-1 or nz-2 depending on the overlap. Each is
c only responsible for sending out the rho is owns.

      integer(MPIISZ):: i,right_nz,onz
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

c     call MPI_ALLTOALLV(rhop,nzsend,izsend,MPI_DOUBLE_PRECISION,
c    &                   rho ,nzrecv,izrecv,MPI_DOUBLE_PRECISION,
c    &                   MPI_COMM_WORLD,mpierror)
c     --- Everything after this point should be delete for all to all

c     --- Send the data out to processors that need it.
c     --- Copy any data locally as well.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if (i /= my_index) then
            w = w + 1
            call MPI_ISEND(rhop(0,0,izsend(i)),int(nzsend(i)*(nxp+1)*(nyp+1),MPIISZ),
     &                     MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          else
            if (loc(rho) .ne. loc(rhop)) then
c             --- should be izrecv(i)+nzrecv(i)-1?
c             --- should be fixed so no unnecessary data is sent!!!
              onz = 1
              if (my_index == nslaves-1) onz = 0
              rho(:,:,izrecv(i):izrecv(i)+nzrecv(i) - onz) =
     &          rhop(:,:,izsend(i):izsend(i)+nzsend(i) - onz)
            endif
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if (i /= my_index) then
          if (nzrecv(i) > 0) then
            call MPI_RECV(rho(0,0,izrecv(i)),int(nzrecv(i)*(nx+1)*(ny+1),MPIISZ),
     &                    MPI_DOUBLE_PRECISION,
     &                    i,messid,MPI_COMM_WORLD,mpistatus,mpierror)
          endif
        endif
      enddo

c     --- This is needed since the sends are done without buffering. This
c     --- only matters though for the SOR field solver which modifies rho
c     --- for optimization purposes.
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)
c     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (lw3dtimesubs) timegetrhoforfieldsolve3d = timegetrhoforfieldsolve3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine perphi3d_slave(phi,nx,ny,nz)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
c  Sets the slices on the exterior of phi for periodicity
c  sets slice at -1 equal to the slice at nz-1
c  sets slice at nz+1 equal to the slice at 1
c  Only first and last processors do anything.
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
      integer(MPIISZ):: messid = 70
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
 
      if (my_index == 0) then
        call MPI_ISEND(phi(0,0,0),int(2*(nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &             int(nslaves-1,MPIISZ),messid,MPI_COMM_WORLD,mpirequest,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_ISEND(phi(0,0,nz-1),int((nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &                 int(0,MPIISZ),messid,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
 
      if (my_index == 0) then
        call MPI_RECV(phi(0,0,-1),int((nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &                int(nslaves-1,MPIISZ),messid,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_RECV(phi(0,0,nz),int(2*(nx+1)*(ny+1),MPIISZ),MPI_DOUBLE_PRECISION,
     &                int(0,MPIISZ),messid,MPI_COMM_WORLD,mpistatus,mpierror)
      endif

      if (my_index == 0 .or. my_index == nslaves-1) then
        call MPI_WAIT(mpirequest,mpistatus,mpierror)
      endif
 
!$OMP MASTER
      if (lw3dtimesubs) timeperphi3d_slave = timeperphi3d_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine getphiforparticles3d(nx,ny,nz,phi,nxp,nyp,nzp,phip)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: phip(0:nxp,0:nyp,-1:nzp+1)

c Get the phi for the extent where the particles are. This gets phi from
c neighboring processors, and at the very least gets phi in the guard planes,
c iz=-1 and +1.

      integer(MPIISZ):: i
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
c     --- Note that special cases are needed for the first and last processors
c     --- so that the guard cells are filled properly.
c     --- Data to be sent
      do i=0,nslaves-1
        izglobal = max(izfsslave(my_index)-1,izpslave(i)-1)
        izmaxp   = izpslave(i)+nzpslave(i)+1
        izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)+1

        izsend(i) = izglobal - izfsslave(my_index)
        nzsend(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        izglobal = max(izpslave(my_index)-1,izfsslave(i)-1)
        izmaxp   = izpslave(my_index)+nzpslave(my_index)+1
        izmaxfs  = izfsslave(i)+nzfsslave(i)+1

        izrecv(i) = izglobal - izpslave(my_index)
        nzrecv(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if ( i /= my_index) then
            w = w + 1
            call MPI_ISEND(phi(0,0,izsend(i)),int(nzsend(i)*(nx+1)*(ny+1),MPIISZ),
     &                     MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          else
            phip(:,:,izrecv(i):izrecv(i)+nzrecv(i)-1) =
     &        phi(:,:,izsend(i):izsend(i)+nzsend(i)-1)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then
          if (nzrecv(i) > 0) then
            call MPI_RECV(phip(0,0,izrecv(i)),int(nzrecv(i)*(nxp+1)*(nyp+1),MPIISZ),
     &                    MPI_DOUBLE_PRECISION,
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
      subroutine getphiforfields3d(nx,ny,nz,phi)
      use Subtimers3d
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)

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
          izglobal = max(izfsslave(my_index)-1,izslave(i)-1)
        else
          izglobal = max(izfsslave(my_index),izslave(i)-1)
        endif
        izmax   = izslave(i)+nzslave(i)+1
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
          izglobal = max(izslave(my_index)-1,izfsslave(i)-1)
        else
          izglobal = max(izslave(my_index)-1,izfsslave(i))
        endif
        izmax   = izslave(my_index)+nzslave(my_index)+1
        if (i == nslaves-1) then
          izmaxfs  = izfsslave(i)+nzfsslave(i)+1
        else
          izmaxfs  = izfsslave(i)+nzfsslave(i)-1
        endif

        izrecv(i) = izglobal - izslave(my_index)
        nzrecv(i) = min(izmax,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if ( i /= my_index) then
            w = w + 1
            nn = nzsend(i)*(nx+1)*(ny+1)
            call MPI_ISEND(phi(0,0,izsend(i)),nn,
     &                     MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then
          if (nzrecv(i) > 0) then
            nn = nzrecv(i)*(nx+1)*(ny+1)
            call MPI_RECV(phi(0,0,izrecv(i)),nn,
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
