#include "top.h"
c===========================================================================
c W3DSLAVE.M, version $Revision: 1.13 $, $Date: 2002/08/29 18:23:46 $
c Slave routines associated with the W3D package.
c===========================================================================
c===========================================================================
      subroutine init_w3d_parallel
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
      if (fstype == 3 .or. fstype == 7 .or. fstype == 8) then
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
c     --- through nzslave).
      if (lfsautodecomp .or. (fstype /= 3 .and. fstype /= 7)) then

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
        do while (i < nslaves-1 .and. izfsslave(i)+nzfsslave(i) < nz)
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
            if (nzfsslave(i) < 3) then
              i = i - 1
              nzfsslave(i) = nz - izfsslave(i)
            endif
          endif

        enddo

c       --- Save the number of processors that have part of the grid assigned
c       --- to them.
        nslaves = i+1

      else

c       --- nzfsslaves is assumed to be input by the user and is assumed to
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
c       --- in the array zslave, which is an unscaled weighting of the z-ranges
c       --- of the particles for each processor.

c       --- Get sum of zslave to allow proper scaling.
        sumzslave = zslave(0)
        do i=1,nslaves-1
          sumzslave = sumzslave + zslave(i)
        enddo

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
        zpslmax(nslaves-1) = zmmax

c       --- Set iz and nz. This is done so that zmesh(izpslave) < zpslmin, and
c       --- zmesh(izpslave+nzpslave) > zpslmax.
        do i=0,nslaves-1
          izpslave(i) = int((zpslmin(i) - zmmin)/dz)
          nzpslave(i) = int((zpslmax(i) - zmmin)/dz) - izpslave(i) + 1
        enddo

c       --- Make sure that the last processors doesn't have grid cells
c       --- sticking out the end.
        nzpslave(nslaves-1) = nz - izpslave(nslaves-1)

      endif

c---------------------------------------------------------------------------
c     --- Now set the axial extent of each slaves domain, which includes
c     --- both the particle and field solve domain.
      if(solvergeom==XYZgeom) then
       do i=0,nslaves-1
        izslave(i) = min(izpslave(i),izfsslave(i))
        nzslave(i) = max(izpslave(i) + nzpslave(i),izfsslave(i) + nzfsslave(i))
     &               - izslave(i)
        zmslmin(i) = izslave(i)*dz + zmmin
        zmslmax(i) = (izslave(i) + nzslave(i))*dz + zmmin
       enddo
      else
       do i=0,nslaves-1
        izslave(i) = izfsslave(i)
        nzslave(i) = izfsslave(i) + nzfsslave(i)
     &               - izslave(i)
        zmslmin(i) = izslave(i)*dz + zmmin
        zmslmax(i) = (izslave(i) + nzslave(i))*dz + zmmin
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
      nz = nzslave(my_index)
      nzzarr = nzpslave(my_index)
      nzlmax = nzpslave(my_index)
      nzmmnt = nzpslave(my_index)
      zpmin = zmmin + izpslave(my_index)*dz
      zpmax = (izpslave(my_index)+nzpslave(my_index))*dz + zmmin
      if (zzmin == 0.) zzmin = zpmin
      if (zzmax == 0.) zzmax = zpmax
      if (zlmin == 0.) zlmin = zpmin
      if (zlmax == 0.) zlmax = zpmax
      if (zmmntmin == 0.) zmmntmin = zpmin
      if (zmmntmax == 0.) zmmntmax = zpmax
      izfsmin = izfsslave(my_index) - izslave(my_index)
      izfsmax = izfsmin + nzfsslave(my_index)
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

      return
      end
c===========================================================================
c===========================================================================
      subroutine reorg_particles
      use Parallel
      use Picglb
      use Picglb3d
      use InMesh3d
      use InPart
      use Particles

c Reorganize particles:  distributes particles to the appropriate slaves.
c Particles are shifted multiple times until they are in the proper region.

      integer(ISZ):: checkzpartbnd
      integer(ISZ):: nout1,nout2
      integer(ISZ):: mpierror
      include "mpif.h"

c First, call zpartbnd to redistribute particles.  This checks if there are
c any particles outside of the mesh and if so, loops until there aren't any.
c This loop is necessary since in the initial load, the particles may be
c loaded more than one region away from where the will be in the simulation,
c but the routine zpartbnd only moves particles to the nearest neighbors.
      nout1 = checkzpartbnd(zgrid)
      call MPI_ALLREDUCE(nout1,nout2,1,MPI_INTEGER,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      do while (nout2 > 0)
        call zpartbnd(zmmax,zmmin,dz,zgrid)
        nout1 = checkzpartbnd(zgrid)
        call MPI_ALLREDUCE(nout1,nout2,1,MPI_INTEGER,MPI_MAX,
     &                     MPI_COMM_WORLD,mpierror)
      enddo

      return
      end
c===========================================================================
      subroutine sw_globalsum(ns,sw)
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
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo
      swtemp = sw
      call MPI_ALLREDUCE(swtemp,sw,ns,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)
      do is=1,ns
        if (sw(is) /= 0.) sw(is) = 1./sw(is)
      enddo

      return
      end
c===========================================================================
      subroutine perrho3d_slave(rho,nx,ny,nz,periinz)
      use Parallel
      use Databuffers
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      logical(ISZ):: periinz

c This routine sums rho in the overlapping boundaries with the neighboring
c processors. Note that this does not send the summed rho back to the left
c since that is done in getrhoforfieldsolve.
c This allows for cases where more than two processors have deposited rho
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

c     --- Allocate some buffer space.
      if (b2d1 < nx .or. b2d2 < ny) then
        b2d1 = nx
        b2d2 = ny
        call gchange("Databuffers",0)
      endif

c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Location where rho planes are to be exchanged, and number of
c     --- planes to exchange.
      if (my_index > 0) then
        left_iz = izpslave(my_index) - izslave(my_index)
        left_nz=izpslave(my_index-1)+nzpslave(my_index-1)-izpslave(my_index)+1
        left_nz = min(left_nz, 2)
      else
        left_iz = 0
        left_nz = 1
      endif
      if (my_index < nslaves-1) then
        right_iz = izpslave(my_index+1) - izslave(my_index)
        right_nz=izpslave(my_index)+nzpslave(my_index)-izpslave(my_index+1)+1
        right_nz = min(right_nz, 2)
      else
        right_iz = nz
        right_nz = 1
      endif

c     --- Check if data being sent overlaps with data being received.
c     --- If so, wait for the incoming data first.
      if (my_index > 0 .and. right_iz <= left_iz+left_nz) then
        call MPI_RECV(buffer2d,left_nz*(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                left_pe,0,MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Sum the two planes. The utility call sumarry is used since
c       --- with the gchange above, the local dimension of buffer2d may not
c       --- be correct, so explicitly writing out the loops may not work.
        call sumarry(buffer2d,rho(0,0,left_iz),left_nz*(nx+1)*(ny+1))

      endif

c     --- Send rho to the process to the right
      if (my_index < nslaves-1 .or. periinz) then
        call MPI_ISEND(rho(0,0,right_iz),right_nz*(nx+1)*(ny+1),
     &                 MPI_DOUBLE_PRECISION,
     &                 right_pe,0,MPI_COMM_WORLD,mpirequest,mpierror)
      endif

c     --- If there is no overlap, wait for the data after sending so
c     --- that processors to the right do not have to also wait.
      if ((my_index > 0 .and. right_iz > left_iz+left_nz) .or.
     &    (my_index == 0 .and. periinz)) then
        call MPI_RECV(buffer2d,left_nz*(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                left_pe,0,MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Sum the two planes. The utility call sumarry is used since
c       --- with the gchange above, the local dimension of buffer2d may not
c       --- be correct, so explicitly writing out the loops may not work.
        call sumarry(buffer2d,rho(0,0,left_iz),left_nz*(nx+1)*(ny+1))

      endif

      return
      end
c=============================================================================
      subroutine getrhoforfieldsolve(nx,ny,nz,rho)
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)

c Gather the charge density in the region where the field solve
c will be done. This assumes that each processor "owns" rho from
c iz=0 to either iz=nz-1 or nz-2 depending on the overlap. Each is
c only responsible for sending out the rho is owns.

      integer(ISZ):: i,right_nz
      integer(ISZ):: izsend,nzsend,izrecv,nzrecv
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpierror,mpirequest

c     --- First, send the data out to processors that need it.
      right_nz = 1
      if (my_index < nslaves-1)
     &  right_nz=izpslave(my_index)+nzpslave(my_index)-izpslave(my_index+1)+1
      do i=0,nslaves-1
        if (i /= my_index) then
          izsend = max(izpslave(my_index),izfsslave(i)) - izslave(my_index)
          nzsend = min(izpslave(my_index)+nzpslave(my_index)-right_nz,
     &                 izfsslave(i)+nzfsslave(i))
     &             - max(izpslave(my_index),izfsslave(i)) + 1
          if (nzsend > 0) then
            call MPI_ISEND(rho(0,0,izsend),nzsend*(nx+1)*(ny+1),
     &                     MPI_DOUBLE_PRECISION,
     &                     i,0,MPI_COMM_WORLD,mpirequest,mpierror)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if (i /= my_index) then
          right_nz = 1
          if (i < nslaves-1)
     &      right_nz=izpslave(i)+nzpslave(i)-izpslave(i+1)+1
          izrecv = max(izfsslave(my_index),izpslave(i)) - izslave(my_index)
          nzrecv = min(izfsslave(my_index)+nzfsslave(my_index),
     &                 izpslave(i)+nzpslave(i)-right_nz)
     &             - max(izfsslave(my_index),izpslave(i)) + 1
          if (nzrecv > 0) then
            call MPI_RECV(rho(0,0,izrecv),nzrecv*(nx+1)*(ny+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    i,0,MPI_COMM_WORLD,mpistatus,mpierror)
          endif
        endif
      enddo

c     --- This is needed since the sends are done without buffering. This
c     --- only matters though for the SOR field solver which modifies rho
c     --- for optimization purposes.
      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

      return
      end
c=============================================================================
      subroutine perphi3d_slave(phi,nx,ny,nz)
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)
c  Sets the slices on the exterior of phi for periodicity
c  sets slice at -1 equal to the slice at nz-1
c  sets slice at nz+1 equal to the slice at 1
c  Only first and last processors do anything.
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror
 
      if (my_index == 0) then
        call MPI_ISEND(phi(0,0,0),2*(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 nslaves-1,0,MPI_COMM_WORLD,mpirequest,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_ISEND(phi(0,0,nz-1),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 0,0,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
 
      if (my_index == 0) then
        call MPI_RECV(phi(0,0,-1),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                nslaves-1,0,MPI_COMM_WORLD,mpistatus,mpierror)
      else if (my_index == nslaves-1) then
        call MPI_RECV(phi(0,0,nz),2*(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                0,0,MPI_COMM_WORLD,mpistatus,mpierror)
      endif
 
      return
      end
c=============================================================================
      subroutine getphiforparticles(nx,ny,nz,phi)
      use Parallel
      use InjectVars3d
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)

c Get the phi for the extent where the particles are. This gets phi from
c neighboring processors, and at the very least gets phi in the guard planes,
c iz=-1 and +1.

      integer(ISZ):: i
      integer(ISZ):: izglobal,izmaxp,izmaxfs
      integer(ISZ):: izsend,nzsend,izrecv,nzrecv
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpirequest,mpierror

c     --- First, send the data out to processors that need it.
      do i=0,nslaves-1
        if ( i /= my_index) then

c         izsend = max(izfsslave(my_index),izpslave(i)-1) - izslave(my_index)
c         nzsend = min(izfsslave(my_index)+nzfsslave(my_index)-grid_overlap,
c    &                 izpslave(i)+nzpslave(i)+1)
c    &             - max(izfsslave(my_index),izpslave(i)-1) + 1

          izglobal = max(izfsslave(my_index),izpslave(i)-1)
          izmaxp   = izpslave(i)+nzpslave(i)+1
          izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)-grid_overlap

          izsend   = izglobal - izslave(my_index)
          nzsend   = min(izmaxp,izmaxfs) - izglobal + 1

          if (nzsend > 0) then
            call MPI_ISEND(phi(0,0,izsend),nzsend*(nx+1)*(ny+1),
     &                     MPI_DOUBLE_PRECISION,
     &                     i,0,MPI_COMM_WORLD,mpirequest,mpierror)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then

c         izrecv = max(izpslave(my_index)-1,izfsslave(i)) - izslave(my_index)
c         nzrecv = min(izpslave(my_index)+nzpslave(my_index)+1,
c    &                 izfsslave(i)+nzfsslave(i)-grid_overlap)
c    &             - max(izpslave(my_index)-1,izfsslave(i)) + 1

          izglobal = max(izpslave(my_index)-1,izfsslave(i))
          izmaxp   = izpslave(my_index)+nzpslave(my_index)+1
          izmaxfs  = izfsslave(i)+nzfsslave(i)-grid_overlap

          izrecv   = izglobal - izslave(my_index)
          nzrecv   = min(izmaxp,izmaxfs) - izglobal + 1

          if (nzrecv > 0) then
            call MPI_RECV(phi(0,0,izrecv),nzrecv*(nx+1)*(ny+1),
     &                    MPI_DOUBLE_PRECISION,
     &                    i,0,MPI_COMM_WORLD,mpistatus,mpierror)
          endif
        endif
      enddo
 
c     --- This is needed since the sends are done without buffering.
c     --- No problems have been seem without it, but this just for robustness.
      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

      return
      end                                                                       
c===========================================================================
      subroutine seteears_slave(workspace)
      use Beam_acc
      use InGen
      use InPart
      use InMesh3d
      use Picglb3d
      use Constant
      use Particles
      use Z_arrays
      use Fields3d
      use Parallel
      use Databuffers
      real(kind=8):: workspace(0:*)

c Gather ezax from all slaves into eearsofz array, modify eearsofz
c as appropriate, and redistribute eearsofz.
c The first processor does the gathering and redistributing.
c This assumes that nzzarr = nz
c
c Note: this is really not necessary since each processor has enough
c information to do this calculation on its own. But since it is working
c and only done at the beginning, there is little incentive to fix this.

      integer(ISZ):: nzla,nzlb,nzma,nzmb,iz
      real(kind=8):: zm,zs,zl,eearsprs
      integer(ISZ):: left_pe,right_pe,i
      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE),mpierror

c     --- Gather data from all processes
      if (b1size < 2*nslaves) then
        b1size = 2*nslaves
        call gchange("Databuffers",0)
      endif
c     --- Get number of data points sent from each processor
      call MPI_GATHER(nzzarr,1,MPI_INTEGER,buffer1(1),1,MPI_INTEGER,0,
     &                MPI_COMM_WORLD,mpierror)
c     --- Accumulate the displacments
      buffer1(nslaves+1) = 0
      do i=2,nslaves
        buffer1(nslaves+i) = buffer1(nslaves+i) + buffer1(i)
      enddo
c     --- Gather the data.
      call MPI_GATHERV(ezax,nzzarr,MPI_DOUBLE_PRECISION,
     &                 workspace,buffer1(1),buffer1(nslaves+1),
     &                 MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierror)

c     --- Processor 0 does all of the work
      if (my_index == 0) then

c       --- Now adjust Eears appropriately. (Any changes made in the
c       --- seteears subroutine should also be made here and vice versa.)
        zl = (zimax - zimin)*0.5
        nzla = nzfull/2-int(zl*dzzi)
        nzlb = nzfull/2+int(zl*dzzi)
        zs = (zimax - zimin)*straight*0.5
        nzma = nzfull/2-int(zs*dzzi)
        nzmb = nzfull/2+int(zs*dzzi)
        do iz=0,nzla
          workspace(iz) = - workspace(nzla)
        enddo
        do iz=nzla+1,nzma+5
          workspace(iz) = - workspace(iz)
        enddo
        do iz=nzma+6,nzmb-6
          workspace(iz) = 0.
        enddo
        do iz=nzmb-5,nzlb
          workspace(iz) = - workspace(iz)
        enddo
c       --- Note that workspace(nzlb) has already beed negated.
        do iz=nzlb+1,nzfull
          workspace(iz) =   workspace(nzlb)
        enddo
c       --- Add on linear pressure term
        zm = (zimax - zimin)*(1. - straight)*0.5
        emitlong = zm*2.*vthz/vbeam
        if (emitlong /= 0.) then
          eearsprs = - sm(1)*vbeam**2*emitlong**2/(sq(1)*zm**4)
          do iz=0,nzma
            workspace(iz) = workspace(iz) + eearsprs*(zzmin + iz*dzz + zs)
          enddo
          do iz=nzmb,nzfull
            workspace(iz) = workspace(iz) + eearsprs*(zzmin + iz*dzz - zs)
          enddo
        endif

      endif

c     --- Send data back to the processors
      call MPI_SCATTERV(workspace,buffer1(1),buffer1(nslaves+1),
     &                  MPI_DOUBLE_PRECISION,
     &                  eearsofz,nzzarr,MPI_DOUBLE_PRECISION,
     &                  0,MPI_COMM_WORLD,mpierror)


c     --- Now get overlapped data
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)
      call MPI_SENDRECV(eearsofz(0),1,MPI_DOUBLE_PRECISION,left_pe,0,
     &                  eearsofz(nz),1,MPI_DOUBLE_PRECISION,right_pe,0,
     &                  MPI_COMM_WORLD,mpistatus,mpierror)

      return
      end
c=============================================================================
