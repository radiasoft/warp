#include "top.h"
c=============================================================================
c F3DSLAVE.M, version $Revision: 1.23 $, $Date: 2004/11/25 19:54:38 $
c Slave routines related to F3D package.
c=============================================================================
c=============================================================================
      subroutine exchange_phi(nx,ny,nz,phi,bound0,boundnz,zsend)
      use Subtimersf3d
      use Parallel
      integer(ISZ):: nx,ny,nz,bound0,boundnz,zsend
      real(kind=8):: phi(0:nx,0:ny,-1:nz+1)

c This routine sends out and receives boundary data for the SOR field solver.

      integer(ISZ):: left_pe,right_pe
      integer(ISZ):: left_recv,right_recv,left_send,right_send
      integer(ISZ):: waitcount,waitstart

      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE,4),mpirequest(4),mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- calculate index of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Location where incoming phi plane is to be loaded.  Normally, it
c     --- is put at the edge of the grid, either at iz=0 or iz=nz.  For
c     --- periodic boundaries, it is put just outside the grid, at iz=-1
c     --- or iz=nz+1.  For the above, zsend=0.  At the end of the field solve,
c     --- phi at the outside of the grid is passed, so phi is put at iz=-1 and
c     --- iz=nz+1 so zsend needs to be -1.
      left_recv = zsend
      right_recv = zsend
      if (bound0 == 2) left_recv = -1
      if (boundnz == 2) right_recv = -1

c     --- Location from where incoming phi is obtained.  The same cases as
c     --- above are used.
      left_send = 1 - zsend
      right_send = 1 - zsend
      if (bound0 == 2) left_send = 1
      if (boundnz == 2) right_send = 1

      waitcount = 4
      waitstart = 1
      if (my_index == 0 .and. bound0 /= 2) then
        waitcount = 2
        waitstart = 3
      else if (my_index == nslaves-1 .and. boundnz /= 2) then
        waitcount = 2
        waitstart = 1
      endif

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.
      if (my_index > 0 .or. bound0 == 2) then
        call MPI_ISEND(phi(0,0,left_send),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 left_pe,nx+ny,MPI_COMM_WORLD,mpirequest(1),mpierror)
        call MPI_IRECV(phi(0,0,left_recv),(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                 left_pe,nx+ny,MPI_COMM_WORLD,mpirequest(2),mpierror)
      endif

      if (my_index < nslaves-1 .or. boundnz == 2) then
        call MPI_ISEND(phi(0,0,nz-right_send),(nx+1)*(ny+1),
     &                 MPI_DOUBLE_PRECISION,right_pe,nx+ny,
     &                 MPI_COMM_WORLD,mpirequest(3),mpierror)
        call MPI_IRECV(phi(0,0,nz-right_recv),(nx+1)*(ny+1),
     &                 MPI_DOUBLE_PRECISION,right_pe,nx+ny,
     &                 MPI_COMM_WORLD,mpirequest(4),mpierror)
      endif

c     --- Now wait for everything to finish
      call MPI_WAITALL(waitcount,mpirequest(waitstart),mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timeexchange_phi = timeexchange_phi + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      subroutine transpose(nx,ny,nz,phi,nx_tran,ny_tran)
      use Subtimersf3d
      use Parallel
      use Transpose_work_space
      integer(ISZ):: nx,ny,nz,nx_tran,ny_tran
      real(kind=8):: phi(0:nx,0:ny,0:nz)

c Transpose 3-D array.  Only the volume (0:nx-1,0:ny-1,0:nz-1) is transposed.
c Transpose is NOT done in place but uses scratch array.
c
c NOTICE: In cases where (nx+ny) < nslaves, the size of the transposed
c phi is LARGER than the size of phi. Extra padding must be added to the
c end of phi of size
c izextra = (nzfull - nzfull/nslaves*(nx+ny)/(nx*ny))
c and with phi dimensions as
c phi(0:nx,0:ny,0:nz+izextra)
c
c Algorithm assumes that all dimensions are powers of two and the number
c of processors is a power of two.

      integer(ISZ):: nn
      integer(ISZ):: is,ix,iy,j,left_pe,right_pe
      include "mpif.h"
      integer(ISZ):: mpinewtype,mpierror
      integer(ISZ):: mpinewtype2,mpidoublesize
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE)
      integer(ISZ):: blen(2),bdisp(2),btype(2)
      integer(ISZ):: sdispls(0:nslaves-1),rdispls(0:nslaves-1)
      integer(ISZ):: sendcounts(0:nslaves-1),recvcounts(0:nslaves-1)
c     integer(ISZ),parameter:: M_TRANSPOSE=10
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Calculate what the dimensions of the grid will be after the
c     --- transpose.  The dimensions will be (0:nx_tran,ny_tran).  The
c     --- full x dimension is used since the z tridiag solver uses the
c     --- elements ix=nx for scratch space.
c     --- The work is divided up first along the y direction.  Since
c     --- all dimensions and number of processors are powers of 2, ny/nslaves
c     --- will always also be a power of 2 (or zero).
c     --- If ny >= nslaves, set first dimension of array after transpose to
c     --- be the same as nx, otherwise work has to also be divided along x.
      if (ny .ge. nslaves) then
        nx_tran = nx
        ny_tran = ny/nslaves
      else
        nx_tran = nx*ny/nslaves
        ny_tran = 1
      endif
      nn = (nx_tran+1)*ny_tran

c     --- Allocate the temporary space to copy the transpose phi into.
c     --- At the end, the data is copy back into the phi array.
      if (phi_trnspssize < nn*nz*nslaves) then
        phi_trnspssize = nn*nz*nslaves
        call gallot("Transpose_work_space",0)
      endif

c     --- Create an MPI type to refer to the distributed pieces of data that
c     --- are sent to each of the processors.
      call MPI_TYPE_VECTOR(nz,nn,(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                     mpinewtype,mpierror)
      call MPI_TYPE_COMMIT(mpinewtype,mpierror)


c     --- The transpose is much less complicated when ny > nslaves.
      if (ny .ge. nslaves) then

c       On the T3E, the loop over slaves using SENDRECV is 20 to 50% faster
c       than the ALLTOALL. Go figure.
c       On the SP, the SENDRECV is only as fast as the call to ALLTOALL.
c       That is now used since it has the potential of improving if the
c       MPI implementation is further optimized.

c       --- Create a second type so that the upper bound (and therefore) the
c       --- extent of the new type can be set. The extent needs to be modified
c       --- since the data being sent with the alltoall is interleaved.
c       --- The extent is set to be equal to the byte length of the contiguous
c       --- chunks of data, i.e. nn times number of bytes per real.
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,mpidoublesize,mpierror)
        blen(1) = 1
        blen(2) = 1
        bdisp(1) = 0
        bdisp(2) = nn*mpidoublesize
        btype(1) = mpinewtype
        btype(2) = MPI_UB
        call MPI_TYPE_STRUCT(2,blen,bdisp,btype,mpinewtype2,mpierror)
        call MPI_TYPE_COMMIT(mpinewtype2,mpierror)

c       --- Do all of the communication at once.
        call MPI_ALLTOALL(phi(0,0,0),1,mpinewtype2,
     &                    phi_trnsps,nn*nz,MPI_DOUBLE_PRECISION,
     &                    MPI_COMM_WORLD,mpierror)

c       --- Free the data type
        call MPI_TYPE_FREE(mpinewtype2,mpierror)

      else
c       --- Things are very different with ny < nslaves
c       --- Note though that the same ALLTOALL code would actually work just
c       --- as well for ny >= nslaves.

c       --- Create a copy of the new data type, but defined with an upper
c       --- bound so that the extent is equal to the size a one double.
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,mpidoublesize,mpierror)
        blen(1) = 1
        blen(2) = 1
        bdisp(1) = 0
        bdisp(2) = mpidoublesize
        btype(1) = mpinewtype
        btype(2) = MPI_UB
        call MPI_TYPE_STRUCT(2,blen,bdisp,btype,mpinewtype2,mpierror)
        call MPI_TYPE_COMMIT(mpinewtype2,mpierror)

c       --- Setup displacement arrays.
        do is = 0,nslaves-1
          ix = mod(is*nx_tran,nx)
          iy = is*ny/nslaves
          sdispls(is) = ix + iy*(ny+1)
          sendcounts(is) = 1
          rdispls(is) = is*nn*nz
          recvcounts(is) = nn*nz
        enddo

c       --- Do all of the communication at once.
        call MPI_ALLTOALLV(phi(0,0,0),sendcounts,sdispls,mpinewtype2,
     &                     phi_trnsps,recvcounts,rdispls,MPI_DOUBLE_PRECISION,
     &                     MPI_COMM_WORLD,mpierror)

      endif

c     --- Free the data type
      call MPI_TYPE_FREE(mpinewtype,mpierror)

c     --- Now the array phi_trnsps contains 3-D data with dimensions
c     --- (0:nx_tran,ny_tran,nz*nslaves).

c     --- Copy phi_trnsps into phi (don't want field solving routine to use
c     --- the temporary array phi_trnsps).  This is somewhat inefficient, but
c     --- hopefully the transpose will eventually be done in place, avoiding
c     --- the temp space.
      call copyarry(phi_trnsps,phi,nn*nz*nslaves)

c DONE!
!$OMP MASTER
      if (lf3dtimesubs) timetranspose = timetranspose + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine transposei(nx,ny,nz,phi,nx_tran,ny_tran)
      use Subtimersf3d
      use Parallel
      use Transpose_work_space
      integer(ISZ):: nx,ny,nz,nx_tran,ny_tran
      real(kind=8):: phi(0:nx,0:ny,0:nz)

c Undo above transpose of 3-D array.  Only the volume (0:nx-1,0:ny-1,0:nz-1)
c is transposed back.  Now, nx_tran and ny_tran are used as input.
c Transpose is NOT done in place but uses scratch array.
c
c See the NOTICE in transpose above.
c
c Algorithm assumes that all dimensions are powers of two and the number
c of processors is a power of two.

      integer(ISZ):: nn
      integer(ISZ):: is,ix,iy,j,left_pe,right_pe
      include "mpif.h"
      integer(ISZ):: mpinewtype,mpierror
      integer(ISZ):: mpinewtype2,mpidoublesize
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE)
      integer(ISZ):: blen(2),bdisp(2),btype(2)
      integer(ISZ):: sdispls(0:nslaves-1),rdispls(0:nslaves-1)
      integer(ISZ):: sendcounts(0:nslaves-1),recvcounts(0:nslaves-1)
c     integer(ISZ),parameter:: M_TRANSPOSE=10
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Allocate the work space.
      if (phi_trnspssize < (nx_tran+1)*ny_tran*nz*nslaves) then
        phi_trnspssize = (nx_tran+1)*ny_tran*nz*nslaves
        call gallot("Transpose_work_space",0)
      endif

c Copy phi into phi_trnsps array  (since field solving routine doesn't use
c the temporary array phi_trnsps).  This is somewhat inefficient, but
c hopefully the transpose will eventually be done in place, avoiding
c the temp space.
      nn = (nx_tran+1)*ny_tran
      call copyarry(phi,phi_trnsps,nn*nz*nslaves)


      call MPI_TYPE_VECTOR(nz,nn,(nx+1)*(ny+1),MPI_DOUBLE_PRECISION,
     &                     mpinewtype,mpierror)
      call MPI_TYPE_COMMIT(mpinewtype,mpierror)

      if (ny .ge. nslaves) then

        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,mpidoublesize,mpierror)
        blen(1) = 1
        blen(2) = 1
        bdisp(1) = 0
        bdisp(2) = nn*mpidoublesize
        btype(1) = mpinewtype
        btype(2) = MPI_UB
        call MPI_TYPE_STRUCT(2,blen,bdisp,btype,mpinewtype2,mpierror)
        call MPI_TYPE_COMMIT(mpinewtype2,mpierror)

c       --- Do all of the communication at once.
        call MPI_ALLTOALL(phi_trnsps,nn*nz,MPI_DOUBLE_PRECISION,
     &                    phi(0,0,0),1,mpinewtype2,
     &                    MPI_COMM_WORLD,mpierror)

c       --- Free the data type
        call MPI_TYPE_FREE(mpinewtype2,mpierror)

      else

c       --- Create a copy of the new data type, but defined with an upper
c       --- bound so that the extent is equal to the size a one double.
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,mpidoublesize,mpierror)
        blen(1) = 1
        blen(2) = 1
        bdisp(1) = 0
        bdisp(2) = mpidoublesize
        btype(1) = mpinewtype
        btype(2) = MPI_UB
        call MPI_TYPE_STRUCT(2,blen,bdisp,btype,mpinewtype2,mpierror)
        call MPI_TYPE_COMMIT(mpinewtype2,mpierror)

c       --- Setup displacement arrays.
        do is = 0,nslaves-1
          sdispls(is) = is*nn*nz
          sendcounts(is) = nn*nz
          ix = mod(is*nx_tran,nx)
          iy = is*ny/nslaves
          rdispls(is) = ix + iy*(ny+1)
          recvcounts(is) = 1
        enddo

c       --- Do all of the communication at once.
        call MPI_ALLTOALLV(phi_trnsps,sendcounts,sdispls,MPI_DOUBLE_PRECISION,
     &                     phi(0,0,0),recvcounts,rdispls,mpinewtype2,
     &                     MPI_COMM_WORLD,mpierror)

      endif

c     --- Free the data type
      call MPI_TYPE_FREE(mpinewtype,mpierror)

c     --- Now the array phi is returned to normal.

c DONE!
!$OMP MASTER
      if (lf3dtimesubs) timetransposei = timetransposei + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine lantzsolver(iwhich,a,kxsq,kysq,kzsq,attx,atty,attz,
     &                       filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,
     &                       l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Subtimersf3d
      use Constant
      use Parallel
      use LantzSolverTemp
      integer(ISZ):: iwhich,nx,ny,nz
      real(kind=8):: a(0:nx,0:ny,-1:nz)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz-1)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz-1)
      real(kind=8):: filt(5,3)
      real(kind=8):: lx,ly,lz,scrtch(*),xywork(*),zwork(2,0:nx,0:nz)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c Poisson solver based on the method developed by Stephen Lantz.
c First, FFT's are applied to the transverse planes, reducing Poisson's
c equation to a tridiagonal matrix. Then, row reduction is done on the
c matrix in each process, changing the form of the matrix into one resembing
c the capital letter 'N'. In this form, the first and last rows of the matrix
c on each processor can be seperated out to form an independent tridiagonal
c matrix of smaller size than the original. This smaller matrix is transposed,
c solved and transposed back. The remaining interior elements can then be
c directly calculated from the solution of the smaller tridiag system.
c Finally, inverse FFT's are applied transversely.

      integer(ISZ):: ikxmin,ikymin,ikxm,ikym,jx,jy
      real(kind=8):: norm,t
      integer(ISZ):: ix,iy,iz,nx_tran,ny_tran
      logical(ISZ):: lalloted
      data lalloted/.false./
      save lalloted
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Allocate space if needed or requested
      if (iwhich == 0 .or. iwhich == 1 .or. .not. lalloted) then
        nxlan = nx
        nylan = ny
        nzlan = nz
        nzfulllan = 2*nslaves
        if (ny .ge. nslaves) then
          nxtranlan = nx
          nytranlan = ny/nslaves
        else
          nxtranlan = nx*ny/nslaves
          nytranlan = 1
        endif
        call gchange("LantzSolverTemp",0)
        lalloted = .true.
c       --- Initialize the k's
        call vpois3d(1,a(0,0,0),a(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry,bound0,boundnz,boundxy)
      endif

      if (iwhich == 1) return

      print*,0,a(4,5,4)
c     --- Apply FFT's to transverse planes.
      call vpois3d(12,a(0,0,0),a(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

      print*,1,a(4,5,4)
c     --- Minimum index of kxsq and kysq that is used
      ikxmin = 1
      if (l4symtry) ikxmin = 0
      ikymin = 1
      if (l2symtry .or. l4symtry) ikymin = 0

c     --- Normalization factor
      norm = (lz/nz)**2/eps0

c     --- Initial setup for the solve
      do iz=0,nz-1
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
c           --- the RHS is the transverse FFTs of rho times dz**2/eps0
c           --- This multiply could actually be done within the tridiag
c           --- solving loops below, making them somewhat more complicated.
c           --- That saves a few percent in the run time.
            a(ix,iy,iz) = a(ix,iy,iz)*norm
          enddo
        enddo
      enddo
c     --- set the end points using Dirichlet boundary conditions.
      if (my_index == 0) then
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            a(ix,iy,0) = a(ix,iy,0) + a(ix,iy,-1)
          enddo
        enddo
      endif
      if (my_index == nslaves-1) then
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            a(ix,iy,nz-1) = a(ix,iy,nz-1) + a(ix,iy,nz)
          enddo
        enddo
      endif
      do iy=ikymin,ny-1
        do ix=ikxmin,nx-1
c         --- set some initial values for the matrix diagonals
          blan(ix,iy,0) = 2. + (kxsq(ix)+kysq(iy))*norm
          blan(ix,iy,1) = 2. + (kxsq(ix)+kysq(iy))*norm
          alan(ix,iy,1) = -1.
          clan(ix,iy,nz-2) = -1.
        enddo
      enddo
c     --- Fill in the zeros since the generaltridiag
c     --- solver requires all b's to be /= 0.
c     do iy=0,ny
c       blan(nx,iy,0) = 1.
c       blan(nx,iy,nz-1) = 1.
c     enddo
c     do ix=0,nx
c       blan(ix,ny,0) = 1.
c       blan(ix,ny,nz-1) = 1.
c     enddo
c     if (ikxmin == 1) then
c       do iy=0,ny
c         blan(0,iy,0) = 1.
c         blan(0,iy,nz-1) = 1.
c       enddo
c     endif
c     if (ikymin == 1) then
c       do ix=0,nx
c         blan(ix,0,0) = 1.
c         blan(ix,0,nz-1) = 1.
c       enddo
c     endif

c     --- Downward row reduction
      do iz=2,nz-1
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            t = (-1.)/blan(ix,iy,iz-1)
            alan(ix,iy,iz) = -alan(ix,iy,iz-1)*t
            blan(ix,iy,iz) = blan(ix,iy,0) - (-1.)*t
            a(ix,iy,iz) = a(ix,iy,iz) - a(ix,iy,iz-1)*t
          enddo
        enddo
      enddo
      print*,2,a(4,5,4)

c     --- Upward row reduction
      do iz=nz-1-2,0,-1
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            t = (-1.)/blan(ix,iy,iz+1)
            alan(ix,iy,iz) = alan(ix,iy,iz) - alan(ix,iy,iz+1)*t
            clan(ix,iy,iz) = -clan(ix,iy,iz+1)*t
            a(ix,iy,iz) = a(ix,iy,iz) - a(ix,iy,iz+1)*t
          enddo
        enddo
      enddo
      print*,3,a(4,5,4)

c     --- Assemble reduced rhs
      do iy=ikymin,ny-1
        do ix=ikxmin,nx-1
          dtranlan(ix,iy,0) = a(ix,iy,0)
          dtranlan(ix,iy,1) = a(ix,iy,nz-1)
        enddo
      enddo

c     --- Transpose reduced tridiagonal matrix
      call transpose(nx,ny,2,dtranlan,nx_tran,ny_tran)

c     --- Find location of transverse data relative to original dimensions.
      jx = mod(my_index*nx_tran,nx)
      jy = int(my_index*ny/real(nslaves))
c     --- Check if work includes boundary points.  If so, use precalculated
c     --- min's.
      if (jx == 0) then
        ikxm = ikxmin
      else
        ikxm = 0
      endif
      if (jy == 0) then
        ikym = ikymin
      else
        ikym = 0
      endif

c     --- Assemble reduced matrix.
      do iz=0,2*nslaves-1,2
        do iy=ikym,ny_tran-1
          do ix=ikxm,nx_tran-1
            atranlan(ix,iy,0+iz) = -1.
            atranlan(ix,iy,1+iz) = alan(jx+ix,jy+iy,nz-1)
            btranlan(ix,iy,0+iz) = blan(jx+ix,jy+iy,0)
            btranlan(ix,iy,1+iz) = blan(jx+ix,jy+iy,nz-1)
            ctranlan(ix,iy,0+iz) = clan(jx+ix,jy+iy,0)
            ctranlan(ix,iy,1+iz) = -1.
          enddo
        enddo
      enddo

c     --- Do tridiag solve
      call generaltridiag(nx_tran,ny_tran-1,2*nslaves,atranlan,btranlan,
     &                    ctranlan,dtranlan,ikxmin,ikymin,0)

c     --- Transpose back
      call transposei(nx,ny,2,dtranlan,nx_tran,ny_tran)

c     --- Load solved values back into 'a' array.
      do iy=ikymin,ny-1
        do ix=ikxmin,nx-1
          a(ix,iy,0) = dtranlan(ix,iy,0)
          a(ix,iy,nz-1) = dtranlan(ix,iy,1)
        enddo
      enddo
      print*,4,a(4,5,4)

c     --- Calculate remaining values
      do iz=1,nz-2
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            a(ix,iy,iz) = (a(ix,iy,iz) - alan(ix,iy,iz)*a(ix,iy,0) -
     &                                   clan(ix,iy,iz)*a(ix,iy,nz-1))/
     &                     blan(ix,iy,iz)
          enddo
        enddo
      enddo
      print*,5,a(4,5,4)

c     --- Apply inverse FFT's to transverse planes.
      call vpois3d(13,a(0,0,0),a(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

      print*,6,a(4,5,4)
!$OMP MASTER
      if (lf3dtimesubs) timelantzsolver = timelantzsolver + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine generaltridiag(nx,ny,nz,a,b,c,d,ikxmin,ikymin,jey)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,ikxmin,ikymin,jey
      real(kind=8):: a(0:nx,0:ny,nz)
      real(kind=8):: b(0:nx,0:ny,nz)
      real(kind=8):: c(0:nx,0:ny,nz)
      real(kind=8):: d(0:nx,0:ny,nz)

      integer(ISZ):: ix,iy,iz
      real(kind=8):: t
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      do iz = 2,nz
        do iy=ikymin,ny-jey
          do ix=ikxmin,nx-1
            t = a(ix,iy,iz)/b(ix,iy,iz-1)
            b(ix,iy,iz) = b(ix,iy,iz) - c(ix,iy,iz-1)*t
            d(ix,iy,iz) = d(ix,iy,iz) - d(ix,iy,iz-1)*t
          enddo
        enddo
      enddo

      do iy=ikymin,ny-jey
        do ix=ikxmin,nx-1
          d(ix,iy,nz) = d(ix,iy,nz)/b(ix,iy,nz)
        enddo
      enddo
      do iz = nz-1,1,-1
        do iy=ikymin,ny-jey
          do ix=ikxmin,nx-1
            d(ix,iy,iz) = (d(ix,iy,iz) - c(ix,iy,iz)*d(ix,iy,iz+1))/b(ix,iy,iz)
          enddo
        enddo
      enddo

!$OMP MASTER
      if (lf3dtimesubs) timegeneraltridiag = timegeneraltridiag + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine paralleltridiag(a,kxsq,kysq,kzsq,attx,atty,attz,filt,lx,ly,lz,
     &                     nx,ny,nz,scrtch,xywork,zwork,l2symtry,l4symtry,
     &                     bound0,boundnz,boundxy)
      use Subtimersf3d
      use Constant
      integer(ISZ):: nx,ny,nz
      real(kind=8):: a(0:nx,0:ny,-1:nz+1)
      real(kind=8):: kxsq(0:nx-1),kysq(0:ny-1),kzsq(0:nz-1)
      real(kind=8):: attx(0:nx-1),atty(0:ny-1),attz(0:nz-1)
      real(kind=8):: filt(5,3)
      real(kind=8):: lx,ly,lz,scrtch(*),xywork(0:nx,0:ny),zwork(2,0:nx,0:nz)
      logical(ISZ):: l2symtry,l4symtry
      integer(ISZ):: bound0,boundnz,boundxy

c Experimental parallel Poisson solver. The tridiag solver is used only
c locally. Iteration is done, exchange boundary information. The hope
c is that few enough iterations can be used so that the extra time spent
c on the iterations and communication time will be less than the time of
c a global transpose of the matrix.
c
c The iterations continue until the change in the left hand plane is less
c than the tolerance.

      integer(ISZ):: ikxmin,ikymin
      real(kind=8):: norm
      real(kind=8):: achange,err
      integer(ISZ):: ix,iy
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Apply FFT's to transverse planes.
      call vpois3d(12,a(0,0,0),a(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &                filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,0,
     &                l2symtry,l4symtry,bound0,boundnz,boundxy)

c     --- Minimum index of kxsq and kysq that is used
      ikxmin = 1
      if (l4symtry) ikxmin = 0
      ikymin = 1
      if (l2symtry .or. l4symtry) ikymin = 0

c     --- Normalization factor
      norm = (lz/nz)**2/eps0

c     --- Iteration loop
      achange = LARGEPOS
      do while (achange > 1.e-6)

c       --- Save sample values
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            xywork(ix,iy) = a(ix,iy,1)
          enddo
        enddo

c       --- Exchange phi at the boundaries
        call exchange_phi(nx,ny,nz,a,0,0,0)

c       --- Do the tridiag solve
        call tridiag(nx,ny,nz,nz/2,a(0,0,0),norm,kxsq,kysq,kzsq,
     &               ikxmin,ikymin,1,zwork)

c       --- Get error
        achange = 0.
        do iy=ikymin,ny-1
          do ix=ikxmin,nx-1
            err = abs(a(ix,iy,1) - xywork(ix,iy))
            if (err > achange) achange = err
          enddo
        enddo
        call parallelmaxrealarray(achange,1)

      enddo

c     --- Apply inverse FFT's to transverse planes.
      call vpois3d(13,a(0,0,0),a(0,0,0),kxsq,kysq,kzsq,attx,atty,attz,
     &             filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,0,
     &             l2symtry,l4symtry,bound0,boundnz,boundxy)

!$OMP MASTER
      if (lf3dtimesubs) timeparalleltridiag = timeparalleltridiag + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine parallelgatherall(data,ndata,nproc,nstep)
      use Subtimersf3d
      use Parallel
      integer(ISZ):: ndata,nproc,nstep
      real(kind=8):: data(ndata,nproc,nstep)
c Gathers data in subgroups of processors
      include "mpif.h"
      integer(ISZ):: isend
      integer:: datatype,subcomm,ierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- First, create subgroups of processors
      if (nproc < nslaves) then
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,my_index/nproc,0,subcomm,ierror)
      else
        call MPI_COMM_DUP(MPI_COMM_WORLD,subcomm,ierror)
      endif

c     --- First nstep > 1, create new data type
      if (nstep > 1) then
c       --- NOTE: thist coding is not complete
        call MPI_TYPE_VECTOR(nproc,ndata,ndata*nproc,MPI_DOUBLE_PRECISION,
     &                       datatype,ierror)
        call MPI_TYPE_COMMIT(datatype)
      else
        datatype = MPI_DOUBLE_PRECISION
      endif

c     --- Gather the data
      isend = mod(my_index,nproc) + 1
      call MPI_ALLGATHER(data(1:ndata,isend,1),ndata,MPI_DOUBLE_PRECISION,
     &                   data,ndata,MPI_DOUBLE_PRECISION,subcomm,ierror)

c     --- Delete the subgroup communicator
      call MPI_COMM_FREE(subcomm,ierror)

!$OMP MASTER
      if (lf3dtimesubs) timeparallelgatherall = timeparallelgatherall + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c==== Parallel routines for the multigrid fieldsolver
c=============================================================================
      subroutine mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                      nzfull,nzfullcoarse)
      use Subtimersf3d
      integer(ISZ):: nslaves,nzfull,nzfullcoarse
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)

c Divides nz for each processor by two. For odd numbers, the starting end
c is rounded down, and the upper end up. For each processor, if nz ends up
c being equal to 1, then change it to 2. For the last processor, if nz is 1,
c also, subtract 1 from izfsslave.

      integer(ISZ):: is,nzmin = 2
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      do is=0,nslaves-1

c       --- First divide nz by 2. izfsslave is rounded down, and nzfsslave
c       --- is rounded up.
        izfsslavec(is) = int(floor(1.*izfsslave(is)*nzfullcoarse/nzfull))
        nzfsslavec(is) = int(ceiling(1.*(izfsslave(is)+nzfsslave(is))*
     &                                nzfullcoarse/nzfull - izfsslavec(is)))

c       --- Then check that no nz are less than nzmin.
c       --- Note that if nzfsslave is already < nzmin, leave it as is.
c       --- When nzfsslave < nzmin, in some cases setting nzfsslavec to nzmin
c       --- will cause the coarsened grid to extend too far beyond the finer
c       --- grid's gaurd cell, so it would need data from a neighboring cell,
c       --- where it would normally be able to use data from the gaurd cells.
c       --- That is a rare occurance so it is simpler to avoid the problem
c       --- then to write the code to exchange the data for that case.
        if (nzfsslavec(is) < nzmin .and. nzfsslavec(is) < nzfsslave(is)) then
          nzfsslavec(is) = nzfsslavec(is) + 1
          if (izfsslavec(is)+nzfsslavec(is) > nzfullcoarse) then
            izfsslavec(is) = izfsslavec(is) - 1
            nzfsslavec(is) = nzfullcoarse - izfsslavec(is)
          endif
        endif

      enddo

!$OMP MASTER
      if (lf3dtimesubs) timemgdividenz = timemgdividenz + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                            bounds,nzfull,lzoffset,rzoffset,
     &                            whosendingleft,izsendingleft,
     &                            whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nslaves,my_index
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: bounds(0:5),nzfull
      integer(ISZ):: lzoffset(0:nslaves-1),rzoffset(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1)
      integer(ISZ):: izsendingleft(0:nslaves-1)
      integer(ISZ):: izsendingright(0:nslaves-1)

c Returns list of processors to the left and right to exchange data with.
c The algorithm tries to share the message passing among the
c processors holding the same data. For each processor, it finds the other
c processors which have the data it needs. It then picks the one which
c is sending out the fewest messages so far.

      integer(ISZ):: numsendingleft(0:nslaves-1)
      integer(ISZ):: numsendingright(0:nslaves-1)
      integer(ISZ):: idx,js,ie,mm(1)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- By default, no processors send data (flagged by -1)
      whosendingleft = -1
      numsendingleft = 0
      izsendingleft = 0
      whosendingright = -1
      numsendingright = 0
      izsendingright = 0

c     --- For each processor, find a processor on the right which has the
c     --- data needed.
c     --- Loop over all slaves since this must be a global process.
      do idx=0,nslaves-1

c       --- If this process covers the full extent of the mesh on this and the
c       --- next finer mesh, then it doesn't need any data sent to it.
        if (izfsslave(idx) == 0 .and. izfsslave(idx)+nzfsslave(idx) == nzfull
     &      .and. rzoffset(idx) == 0.)
     &    cycle

c       --- Skip processors which have the same right mesh extent (and
c       --- therefore don't have the data needed for this processor.
        js = idx + 1
        do while (js < nslaves .and.
     &            izfsslave(js)+nzfsslave(js)==izfsslave(idx)+nzfsslave(idx))
          js = js + 1
        enddo

c       --- Processors whose extent extends to the right edge are treated
c       --- differently. They are sent data from the left with periodic
c       --- boundary conditions. Otherwise, all receive data from the last
c       --- processor.
c       --- Do all but that one first.
        if (js < nslaves) then

c         --- Find all processors who data extent overlaps the right
c         --- end of this processor.
          ie = js
          do while (ie < nslaves .and. 
     &              izfsslave(ie) < izfsslave(idx)+nzfsslave(idx))
            ie = ie + 1
          enddo
          ie = ie - 1

        elseif (bounds(5) == 2) then
c         --- Now, treat the right most processor.

c         --- Find all processors who data extent overlaps the right
c         --- end of this processor, wrapping around to the left edge.
          js = 0
          ie = js
          do while (ie < nslaves .and. izfsslave(ie) == 0 .and.
     &              lzoffset(ie) == 0)
            ie = ie + 1
          enddo
          ie = ie - 1
        else

c         --- Some of the processors which extend to the right hand edge will
c         --- need data there if they did not extend to the right edge
c         --- at a the next finer level. These processors are known by
c         --- having rzoffset > 0.
          if (idx < nslaves-1 .and. rzoffset(idx) > 0.) then
            ie = nslaves - 1
            js = ie - 1
            do while (js > idx .and. rzoffset(js) == 0. .and.
     &            izfsslave(js)+nzfsslave(js)==izfsslave(idx)+nzfsslave(idx))
              js = js - 1
            enddo
            js = js + 1
          else
            cycle
          endif
          
        endif

c       --- Get the data from the processor which is doing the least amount
c       --- of sending so far. Note that minloc actually returns an array
c       --- rather than a scalar.
        mm = minloc(numsendingleft(js:ie))
        whosendingleft(idx) = js + mm(1) - 1
        if (js > idx) then
c         --- Difference between data location to be sent and the
c         --- starting value of iz for the sending processor.
          izsendingleft(idx) = izfsslave(idx) + nzfsslave(idx) -
     &                         izfsslave(whosendingleft(idx))
        else
c         --- For periodic boundary conditions, this always has the same value.
          izsendingleft(idx) = 1
        endif

c       --- Increment the number of sends that processor makes.
        numsendingleft(whosendingleft(idx)) = numsendingleft(whosendingleft(idx)) + 1
      enddo

c     --- For each processor, find a processor on the left which has the
c     --- data needed.
c     --- Loop over all slaves since this most be a global process.
      do idx=0,nslaves-1

c       --- If this process covers the full extent on this and the
c       --- next finer mesh, then it doesn't need any data sent to it.
        if (izfsslave(idx) == 0 .and. izfsslave(idx)+nzfsslave(idx) == nzfull
     &      .and. lzoffset(idx) == 0)
     &    cycle

c       --- Skip processors which have the same left mesh extent (and
c       --- therefore don't have the data needed for this processor.
        js = idx - 1
        do while (js >= 0 .and. izfsslave(js) == izfsslave(idx))
          js = js - 1
        enddo

c       --- The left most processor is treated differently (and is only
c       --- sent data from the right with periodic boundary conditions).
c       --- Do all but that one first.
        if (js >= 0) then

c         --- Find all processors who data extent overlaps the left
c         --- end of this processor.
          ie = js
          do while (ie >= 0 .and.
     &              izfsslave(idx) < izfsslave(ie)+nzfsslave(ie))
            ie = ie - 1
          enddo
          ie = ie + 1

        elseif (bounds(4) == 2) then
c         --- Now, treat the left most processor.

c         --- Find all processors who data extent overlaps the left
c         --- end of this processor, wrapping around to the left edge.
          ie = nslaves-1
          js = ie
          do while (js >= 0 .and. izfsslave(js)+nzfsslave(js) == nzfull .and.
     &              rzoffset(js) == 0.)
            js = js - 1
          enddo
          js = js + 1
        else

c         --- Some of the processors which extend to the left hand edge will
c         --- need data there if they did not extend to the left edge
c         --- at a the next finer level. These processors are known by
c         --- having lzoffset > 0.
          if (idx > 0 .and. lzoffset(idx) > 0) then
            js = 0
            ie = js + 1
            do while (ie < idx .and. lzoffset(ie) == 0 .and. izfsslave(ie) == 0)
              ie = ie + 1
            enddo
            ie = ie - 1
          else
            cycle
          endif

        endif

c       --- Get the data from the processor which is doing the least amount
c       --- of sending so far. Note that minloc actually returns an array
c       --- rather than a scalar.
        mm = minloc(numsendingright(js:ie))
        whosendingright(idx) = js + mm(1) - 1
        if (js < idx) then
c         --- Difference between data location to be sent and the
c         --- starting value of iz for the sending processor.
          izsendingright(idx)=izfsslave(idx)-izfsslave(whosendingright(idx))

        else
c         --- For periodic boundary conditions, this always has the same value.
          izsendingright(idx) = nzfsslave(whosendingright(idx)) - 1
        endif

c       --- Increment the number of sends that processor makes.
        numsendingright(whosendingright(idx)) = numsendingright(whosendingright(idx)) + 1
      enddo

!$OMP MASTER
      if (lf3dtimesubs) timemggetexchangepes = timemggetexchangepes + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_phi(nx,ny,nz,nzfull,phi,bounds,zsend,delt,delz,
     &                          my_index,nslaves,izfsslave,nzfsslave,
     &                          whosendingleft,izsendingleft,
     &                          whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,nzfull,zsend,delt,delz
      real(kind=8):: phi(-delt:nx+delt,-delt:ny+delt,-delz:nz+delz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c This routine sends out and receives boundary data for the MG field solver.

      integer(ISZ):: leftdelta,rightdelta
      integer(ISZ):: lr,rr,ls,rs
      integer(ISZ):: ip,w

      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE,2*nslaves+2)
      integer(ISZ):: mpirequest(2*nslaves+2),mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Location where incoming phi plane is to be loaded.  Normally, it
c     --- is put at the edge of the grid, either at iz=0 or iz=nz.  For
c     --- periodic boundaries, it is put just outside the grid, at iz=-1
c     --- or iz=nz+1.  For the above, zsend=0.  At the end of the field solve,
c     --- phi at the outside of the grid is passed, so phi is put at iz=-1 and
c     --- iz=nz+1 so zsend needs to be -1.
      lr = zsend
      rr = nz - zsend
      if (bounds(4) == 2) lr = -1
      if (bounds(5) == 2) rr = nz + 1

c     --- The actual value of zsend used depends on the boundary conditions.
c     --- For periodic boundary conditions, processors which are on the ends
c     --- of the mess, always send the same data, independent of what zsend is.
      leftdelta = zsend
      rightdelta = zsend
      if (izfsslave(my_index) == 0 .and. bounds(4) == 2) leftdelta = 0
      if (izfsslave(my_index)+nzfsslave(my_index) == nzfull .and. bounds(5) == 2)
     &  rightdelta = 0

c     --- Set count to zero
      w = 0

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.
      do ip=0,nslaves-1
        if (my_index == whosendingleft(ip)) then
          w = w + 1
          ls = izsendingleft(ip) - leftdelta
          call MPI_ISEND(phi(-delt,-delt,ls),(nx+1+2*delt)*(ny+1+2*delt),
     &                   MPI_DOUBLE_PRECISION,ip,nzfull,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      enddo
      if (whosendingright(my_index) /= -1) then
        w = w + 1
        call MPI_IRECV(phi(-delt,-delt,lr),(nx+1+2*delt)*(ny+1+2*delt),
     &                 MPI_DOUBLE_PRECISION,whosendingright(my_index),nzfull,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif

      do ip=0,nslaves-1
        if (my_index == whosendingright(ip)) then
          w = w + 1
          rs = izsendingright(ip) + rightdelta
          call MPI_ISEND(phi(-delt,-delt,rs),(nx+1+2*delt)*(ny+1+2*delt),
     &                   MPI_DOUBLE_PRECISION,ip,nzfull,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      enddo
      if (whosendingleft(my_index) /= -1) then
        w = w + 1
        call MPI_IRECV(phi(-delt,-delt,rr),(nx+1+2*delt)*(ny+1+2*delt),
     &                 MPI_DOUBLE_PRECISION,whosendingleft(my_index),nzfull,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif

c     --- Now wait for everything to finish
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phi = timemgexchange_phi + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_phiperiodic(nx,ny,nz,nzfull,phi,bounds,zsend,delt,
     &                                  delz,my_index,nslaves,izfsslave,
     &                                  whosendingleft,whosendingright)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,nzfull,zsend,delt,delz
      real(kind=8):: phi(-delt:nx+delt,-delt:ny+delt,-delz:nz+delz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1),whosendingright(0:nslaves-1)

c This routine sends out and receives boundary data for the MG field solver.
c This only sends plane iz=0 when there are periodic boundary conditions.

      integer(ISZ):: ip,w

      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE,nslaves)
      integer(ISZ):: mpirequest(nslaves),mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Does nothing if boundaries are not periodic.
      if (.not. (bounds(4) == 2 .or. bounds(5) == 2)) return

c     --- Set count to zero
      w = 0

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.
      if (izfsslave(my_index) == 0) then
c       --- The loop only covers processors to the right since any processors
c       --- to the left would not need phi on the left boundary. In some cases,
c       --- this processor my be sending other than boundary data to the left.
c       --- This avoids the send, which would cause this process to lock up
c       --- on the waitall since that message would never be received.
        do ip=my_index+1,nslaves-1
          if (my_index == whosendingleft(ip)) then
            w = w + 1
            call MPI_ISEND(phi(-delt,-delt,zsend),(nx+1+2*delt)*(ny+1+2*delt),
     &                     MPI_DOUBLE_PRECISION,ip,10,
     &                     MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif
      if (izfsslave(my_index)+nz == nzfull .and. zsend /= 0) then
c       --- See comments on loop above.
        do ip=0,my_index-1
          if (my_index == whosendingright(ip)) then
            w = w + 1
            call MPI_ISEND(phi(-delt,-delt,nz-zsend),(nx+1+2*delt)*(ny+1+2*delt),
     &                     MPI_DOUBLE_PRECISION,ip,10,
     &                     MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif

      if (izfsslave(my_index)+nz == nzfull) then
        if (whosendingleft(my_index) /= -1) then
          w = w + 1
          call MPI_IRECV(phi(-delt,-delt,nz+zsend),(nx+1+2*delt)*(ny+1+2*delt),
     &                   MPI_DOUBLE_PRECISION,whosendingleft(my_index),10,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      endif
      if (izfsslave(my_index) == 0 .and. zsend /= 0) then
        if (whosendingright(my_index) /= -1) then
          w = w + 1
          call MPI_IRECV(phi(-delt,-delt,-zsend),(nx+1+2*delt)*(ny+1+2*delt),
     &                   MPI_DOUBLE_PRECISION,whosendingright(my_index),10,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      endif

c     --- Now wait for everything to finish
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phiperiodic = timemgexchange_phiperiodic + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_rho(nx,ny,nz,nzfull,rho,
     &                          my_index,nslaves,izfsslave,nzfsslave,
     &                          lzoffset,rzoffset,
     &                          whosendingleft,izsendingleft,
     &                          whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,nzfull
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: lzoffset(0:nslaves-1),rzoffset(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c This routine sends out and receives boundary data for the MG field solver.

      integer(ISZ):: ls,rs
      integer(ISZ):: ip,w

      include "mpif.h"
      integer(ISZ):: mpistatus(MPI_STATUS_SIZE,2*nslaves+2)
      integer(ISZ):: mpirequest(2*nslaves+2),mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Set count to zero
      w = 0

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.

c     --- The last processor is skipped since it doesn't need to have rho
c     --- on the right end sent to it.
      do ip=0,nslaves-2

c       --- Also skip processors that extend to the right edge of this and the
c       --- next finer level.
        if (izfsslave(ip)+nzfsslave(ip)==nzfull .and. rzoffset(ip)==0.) cycle

        if (my_index == whosendingleft(ip)) then
          w = w + 1
          ls = izsendingleft(ip)

c         --- When my_index < ip, it means that the data being sent is periodic.
c         --- In this case, the izsendingleft(ip) is set to send the guard cell,
c         --- so override it.
          if (my_index < ip) ls = 0

          call MPI_ISEND(rho(0,0,ls),(nx+1)*(ny+1),
     &                   MPI_DOUBLE_PRECISION,ip,nzfull,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)

        endif
      enddo

c     --- Any processor that extends to the left edge of this and the next
c     --- finer level don't need a rho sent to it.
      if (whosendingright(my_index) /= -1 .and. 
     &    .not. (izfsslave(my_index)==0 .and. lzoffset(my_index)==0)) then
        w = w + 1
        call MPI_IRECV(rho(0,0,0),(nx+1)*(ny+1),
     &                 MPI_DOUBLE_PRECISION,whosendingright(my_index),nzfull,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif

c     --- The first processor is skipped since it doesn't need to have rho
c     --- on the left end sent to it.
      do ip=1,nslaves-1

c       --- Also skip processors that extend to the left edge of this and the
c       --- next finer level.
        if (izfsslave(ip)==0 .and. lzoffset(ip)==0) cycle

        if (my_index == whosendingright(ip)) then
          w = w + 1
          rs = izsendingright(ip)

c         --- When my_index > ip, it means that the data being sent is periodic.
c         --- In this case, the izsendingright(ip) is set to send the guard cell,
c         --- so override it.
          if (my_index > ip) rs = nz

          call MPI_ISEND(rho(0,0,rs),(nx+1)*(ny+1),
     &                   MPI_DOUBLE_PRECISION,ip,nzfull,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)

        endif
      enddo

c     --- Any processor that extends to the right edge of this and the next
c     --- finer level don't need a rho sent to it.
      if (whosendingleft(my_index) /= -1 .and.
     &    .not. (izfsslave(my_index)+nzfsslave(my_index)==nzfull .and.
     &           rzoffset(my_index)==0.)) then
        w = w + 1
        call MPI_IRECV(rho(0,0,nz),(nx+1)*(ny+1),
     &                 MPI_DOUBLE_PRECISION,whosendingleft(my_index),nzfull,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif

c     --- Now wait for everything to finish
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_rho = timemgexchange_rho + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine printarray3d(nx,ny,nz,iz0,nz0,izs,a,my_index,nslaves,text,mglevel)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,iz0,nz0,izs,my_index,nslaves,mglevel
      real(kind=8):: a(0:nx,0:ny,iz0:nz)
      character(*):: text

c Specialized routine for printing out an array in parallel, keeping the same
c ordering as would be in serial. Primarily used for debugging purposes.
c It is maintained in case it is needed for future use since there is some subtlty
c to getting it correct.

      integer(ISZ):: ix,iy,iz,lastz,firstz
      include "mpif.h"
      integer(ISZ):: mpierror,mpistatus(MPI_STATUS_SIZE)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

c The value of lastz sets which processor's data is printed out in
c overlapping regions.
c   when lastz == nz, the lower number processor's data is printed
c   when lastz == nz-1 (or -2), the higher number processor's data is printed
      lastz = nz
c     lastz = nz - 2
      lastz = nz - nz0
      if (my_index == nslaves-1) lastz = nz

c Each processor prints out its data and then sends a message to the next so it
c can print its data.
      if (my_index > 0) then
        call MPI_RECV(firstz,1,MPI_INTEGER,my_index-1,999,MPI_COMM_WORLD,mpistatus,mpierror)
        firstz = firstz - izs
      else
        firstz = 0
      endif

c The file needs to be reopened each time since there is no coordination in how
c the processors deal with files.
      if (nslaves == 16) then
        open(22,action="write",position="append",file="out16")
      else
        open(22,action="write",position="append",file="out32")
      endif
      do iz=firstz,lastz
        do iy=0,ny
          do ix=0,nx
            write(22,'(I3,x,I3,x,I3,1pE15.5,x,A,I3,I3)') ix,iy,izs+iz,a(ix,iy,iz),text,mglevel,nz0
          enddo
        enddo
      enddo
      close(22)

      if (my_index < nslaves - 1) then
        call MPI_SEND(izs+lastz+1,1,MPI_INTEGER,my_index+1,999,MPI_COMM_WORLD,mpierror)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timeprintarray3d = timeprintarray3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
