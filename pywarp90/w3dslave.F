#include "top.h"
c===========================================================================
c W3DSLAVE.M, version $Revision: 1.25 $, $Date: 2005/06/15 21:09:34 $
c Slave routines associated with the W3D package.
c===========================================================================
c===========================================================================
      subroutine init_w3d_parallel
      use Subtimers3d
      use Parallel
      use Databuffers
      use InGen
      use InGen3d
      use InDiag
      use InPart
      use InMesh3d
      use Picglb
      use Picglb3d
      use Fields3d
      use Fields3dParticles
      use Particles
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
      logical(ISZ):: w3dinitialized
      data w3dinitialized/.false./
      save w3dinitialized
      logical(ISZ):: initialized
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c This routine should only be called once.
      if (w3dinitialized) return
      w3dinitialized = .true.

c get number of slaves and my_pe
      call MPI_INITIALIZED(initialized,mpierror)
      if (.not. initialized) call MPI_INIT(mpierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nslaves,mpierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_index,mpierror)
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
          zpslmin(i) = izpslave(i)*dz + zmmin
          zpslmax(i) = (izpslave(i)+nzpslave(i))*dz + zmmin
        enddo
        nzpslave(nslaves-1) = nzfsslave(nslaves-1)
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
          izpslave(i) = int((zpslmin(i) - zmmin)/dz)
          nzpslave(i) = int((zpslmax(i) - zmmin)/dz) - izpslave(i) + 1
        enddo

c       --- Make sure that the last processors doesn't have grid cells
c       --- sticking out the end.
        if (izpslave(nslaves-1)+nzpslave(nslaves-1) > nz) then
          nzpslave(nslaves-1) = nz - izpslave(nslaves-1)
        endif

      endif

c---------------------------------------------------------------------------
c     --- Now set the axial extent of each slaves domain, which includes
c     --- both the particle and field solve domain.
      if(solvergeom==XYZgeom) then
       do i=0,nslaves-1
        izslave(i) = izfsslave(i)
        nzslave(i) = nzfsslave(i)
        zmslmin(i) = izfsslave(i)*dz + zmmin
        zmslmax(i) = (izfsslave(i) + nzfsslave(i))*dz + zmmin
       enddo
      else
       do i=0,nslaves-1
        izslave(i) = izfsslave(i)
        nzslave(i) = nzfsslave(i)
        zmslmin(i) = izfsslave(i)*dz + zmmin
        zmslmax(i) = (izfsslave(i) + nzfsslave(i))*dz + zmmin
       end do
      end if

c---------------------------------------------------------------------------
c     --- Save global array sizes
      nzfull = nz
      slavenp = npmax
      np = npmax
      zplmin = zmmin
      zplmax = zmmax

c---------------------------------------------------------------------------
c     --- Reset local values
      npmax = npmax/nslaves*1.1
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

c---------------------------------------------------------------------------
c     --- Set np_s here using slavenp, which is the total number of
c     --- particles, instead of using the value set in alotpart (based
c     --- on npmax, which is only a local number of particles).
c     --- This is done so that in stptcl3d, the total number of particles is
c     --- known and can be divied up among the appropriate processors.
      do is=1,ns
        if (np_s(is) == 0) np_s(is) = slavenp*sp_fract(is)
      enddo

c     --- Set npmax_s here using the local number of particles per
c     --- processor. The estimated scaling of np_s here is the same
c     --- as that used to scale npmax above.
      npmax_s(0) = 0
      do is=1,ns
        if (npmax_s(is) == 0) npmax_s(is) = npmax_s(is-1)+np_s(is)/nslaves*1.1
      enddo

c     --- Set value of izextra, which is the amount of extra space at the
c     --- end of phi. For the case when (nx+ny) < nslaves, phi needs
c     --- to be made larger to hold the transposed array. izextra could
c     --- actualy be one less then the value calculated here, but I don't
c     --- subtract one so there is some extra slop.
      if(solvergeom==XYZgeom) then
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
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo
      swtemp = sw
      call MPI_ALLREDUCE(swtemp,sw,ns,MPI_DOUBLE_PRECISION,MPI_SUM,
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
      subroutine sumrhoondomainboundaries(rhop,nxp,nyp,nzp,
     &                                    my_index,nslaves,izpslave,nzpslave)
      use Subtimers3d
      use GlobalVars
      use Databuffers
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)

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

      integer(ISZ):: left_pe,right_pe
      integer(ISZ):: left_iz,right_iz,left_nz,right_nz
      integer(ISZ):: ix,iy,iz
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpierror,mpirequest
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- Allocate some buffer space.
      if (b2d1 < nxp .or. b2d2 < nyp) then
        b2d1 = nxp
        b2d2 = nyp
        call gchange("Databuffers",0)
      endif

c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Location where rhop planes are to be exchanged, and number of
c     --- planes to exchange.
      if (my_index > 0) then
        left_iz = 0
        left_nz=izpslave(my_index-1)+nzpslave(my_index-1)-izpslave(my_index)+1
        left_nz = min(left_nz, 2)
      else
        left_iz = 0
        left_nz = 1
      endif
      if (my_index < nslaves-1) then
        right_iz = izpslave(my_index+1) - izpslave(my_index)
        right_nz=izpslave(my_index)+nzpslave(my_index)-izpslave(my_index+1)+1
        right_nz = min(right_nz, 2)
      else
        right_iz = nzp
        right_nz = 1
      endif

c     --- Check if data being sent overlaps with data being received.
c     --- If so, wait for the incoming data first.
      if (my_index > 0 .and. right_iz <= left_iz+left_nz) then
        call MPI_RECV(buffer2d,left_nz*(nxp+1)*(nyp+1),MPI_DOUBLE_PRECISION,
     &                left_pe,50,MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Sum the two planes. The utility call sumarry is used since
c       --- with the gchange above, the local dimension of buffer2d may not
c       --- be correct, so explicitly writing out the loops may not work.
        call sumarry(buffer2d,rhop(0,0,left_iz),left_nz*(nxp+1)*(nyp+1))

      endif

c     --- Send rhop to the process to the right
      if (my_index < nslaves-1) then
        call MPI_ISEND(rhop(0,0,right_iz),right_nz*(nxp+1)*(nyp+1),
     &                 MPI_DOUBLE_PRECISION,
     &                 right_pe,50,MPI_COMM_WORLD,mpirequest,mpierror)
      endif

c     --- If there is no overlap, wait for the data after sending so
c     --- that processors to the right do not have to also wait.
      if (my_index > 0 .and. right_iz > left_iz+left_nz) then
        call MPI_RECV(buffer2d,left_nz*(nxp+1)*(nyp+1),MPI_DOUBLE_PRECISION,
     &                left_pe,50,MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Sum the two planes. The utility call sumarry is used since
c       --- with the gchange above, the local dimension of buffer2d may not
c       --- be correct, so explicitly writing out the loops may not work.
        call sumarry(buffer2d,rhop(0,0,left_iz),left_nz*(nxp+1)*(nyp+1))

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
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
 
      if (my_index == nslaves-1) then
        call MPI_ISEND(rho(0,0,nz),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 0,70,MPI_COMM_WORLD,mpirequest,mpierror)
        call MPI_RECV(rho(0,0,nz),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                0,70,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == 0) then
        call MPI_RECV(rhotemp,(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                nslaves-1,70,MPI_COMM_WORLD,mpistatus,mpierror)
        rho(:,:,0) = rho(:,:,0) + rhotemp
        call MPI_ISEND(rho(0,0,0),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 nslaves-1,70,MPI_COMM_WORLD,mpirequest,mpierror)
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
      subroutine getrhoforfieldsolve3d_parallel(nx,ny,nz,rho,nxp,nyp,nzp,rhop,
     &                                          my_index,nslaves,
     &                                          izfsslave,nzfsslave,
     &                                          izpslave,nzpslave)
      use Subtimers3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      integer(ISZ):: nxp,nyp,nzp
      real(kind=8):: rhop(0:nxp,0:nyp,0:nzp)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: izpslave(0:nslaves-1),nzpslave(0:nslaves-1)

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" rho from
c iz=0 to either iz=nz-1 or nz-2 depending on the overlap. Each is
c only responsible for sending out the rho is owns.

      integer(ISZ):: i,right_nz
      integer(ISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(ISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpierror
      integer(ISZ):: mpirequest(nslaves)
      integer(ISZ):: w
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
      enddo

c     --- Send the data out to processors that need it.
c     --- Copy any data locally as well.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if (i /= my_index) then
            w = w + 1
            call MPI_ISEND(rhop(0,0,izsend(i)),nzsend(i)*(nxp+1)*(nyp+1),
     &                     MPI_DOUBLE_PRECISION,
     &                     i,60,MPI_COMM_WORLD,mpirequest(w),mpierror)
          else
            rho(:,:,izrecv(i):izrecv(i)+nzrecv(i)-1) =
     &        rhop(:,:,izsend(i):izsend(i)+nzsend(i)-1)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if (i /= my_index) then
          if (nzrecv(i) > 0) then
            call MPI_RECV(rho(0,0,izrecv(i)),nzrecv(i)*(nx+1)*(ny+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    i,60,MPI_COMM_WORLD,mpistatus,mpierror)
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
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()
 
      if (my_index == 0) then
        call MPI_ISEND(phi(0,0,0),2*(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 nslaves-1,70,MPI_COMM_WORLD,mpirequest,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_ISEND(phi(0,0,nz-1),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 0,70,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
 
      if (my_index == 0) then
        call MPI_RECV(phi(0,0,-1),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                nslaves-1,70,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_RECV(phi(0,0,nz),2*(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                0,70,MPI_COMM_WORLD,mpistatus,mpierror)
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

      integer(ISZ):: i
      integer(ISZ):: izglobal,izmaxp,izmaxfs
      integer(ISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(ISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpirequest(nslaves)
      integer(ISZ):: mpierror,w
      real(kind=8):: substarttime,wtime
      if (lw3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent and received.
c     --- Note that special cases are needed for the first and last processors
c     --- so that the guard cells are filled properly.
c     --- Data to be sent
      do i=0,nslaves-1
        if (my_index == 0) then
          izglobal = max(izfsslave(my_index)-1,izpslave(i)-1)
        else
          izglobal = max(izfsslave(my_index),izpslave(i)-1)
        endif
        izmaxp   = izpslave(i)+nzpslave(i)+1
        if (my_index == nslaves-1) then
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)+1
        else
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)-grid_overlap
        endif

        izsend(i) = izglobal - izfsslave(my_index)
        nzsend(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if (i == 0) then
          izglobal = max(izpslave(my_index)-1,izfsslave(i)-1)
        else
          izglobal = max(izpslave(my_index)-1,izfsslave(i))
        endif
        izmaxp   = izpslave(my_index)+nzpslave(my_index)+1
        if (i == nslaves-1) then
          izmaxfs  = izfsslave(i)+nzfsslave(i)+1
        else
          izmaxfs  = izfsslave(i)+nzfsslave(i)-grid_overlap
        endif

        izrecv(i) = izglobal - izpslave(my_index)
        nzrecv(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if ( i /= my_index) then
            w = w + 1
            call MPI_ISEND(phi(0,0,izsend(i)),nzsend(i)*(nx+1)*(ny+1),
     &                     MPI_DOUBLE_PRECISION,
     &                     i,80,MPI_COMM_WORLD,mpirequest(w),mpierror)
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
            call MPI_RECV(phip(0,0,izrecv(i)),nzrecv(i)*(nxp+1)*(nyp+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    i,80,MPI_COMM_WORLD,mpistatus,mpierror)
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
