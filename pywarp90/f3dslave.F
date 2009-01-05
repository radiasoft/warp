#include "top.h"
c=============================================================================
c F3DSLAVE.M, version $Revision: 1.40 $, $Date: 2009/01/05 23:29:38 $
c Slave routines related to F3D package.
c=============================================================================
c=============================================================================
      subroutine exchange_phi(nx,ny,nz,phi,bound0,boundnz,zsend)
      use Subtimersf3d
      use Parallel
      integer(ISZ):: nx,ny,nz,bound0,boundnz,zsend
      real(kind=8):: phi(-1:nx+1,-1:ny+1,-1:nz+1)

c This routine sends out and receives boundary data for the SOR field solver.

      integer(MPIISZ):: left_pe,right_pe
      integer(MPIISZ):: left_recv,right_recv,left_send,right_send
      integer(MPIISZ):: waitcount,waitstart

      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,4),mpirequest(4),mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- calculate index of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nzprocs+my_index-1,nzprocs)
      right_pe = mod(my_index+1,nzprocs)

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
      else if (my_index == nzprocs-1 .and. boundnz /= 2) then
        waitcount = 2
        waitstart = 1
      endif

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.
      if (my_index > 0 .or. bound0 == 2) then
        call MPI_ISEND(phi(-1,-1,left_send),int((nx+3)*(ny+3),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,left_pe,int(nx+ny,MPIISZ),
     &                 MPI_COMM_WORLD,mpirequest(1),mpierror)
        call MPI_IRECV(phi(-1,-1,left_recv),int((nx+3)*(ny+3),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,left_pe,int(nx+ny,MPIISZ),
     &                 MPI_COMM_WORLD,mpirequest(2),mpierror)
      endif

      if (my_index < nzprocs-1 .or. boundnz == 2) then
        call MPI_ISEND(phi(-1,-1,nz-right_send),int((nx+3)*(ny+3),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,right_pe,int(nx+ny,MPIISZ),
     &                 MPI_COMM_WORLD,mpirequest(3),mpierror)
        call MPI_IRECV(phi(-1,-1,nz-right_recv),int((nx+3)*(ny+3),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,right_pe,int(nx+ny,MPIISZ),
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
c=============================================================================
c=============================================================================
      subroutine copy_tophitranspose(nx,ny,nz,phi,phi_trnsps)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nz)
      real(kind=8):: phi_trnsps(0:nx,0:ny,0:nz)

      integer(ISZ):: ix,iy,iz

      do iz=0,nz
        do iy=0,ny
          do ix=0,nx
            phi_trnsps(ix,iy,iz) = phi(ix,iy,iz)
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine copy_fromphitranspose(nx,ny,nz,phi,phi_trnsps)
      integer(ISZ):: nx,ny,nz
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nz)
      real(kind=8):: phi_trnsps(0:nx,0:ny,0:nz)

      integer(ISZ):: ix,iy,iz

      do iz=nz,0,-1
        do iy=ny,0,-1
          do ix=nx,0,-1
            phi(ix,iy,iz) = phi_trnsps(ix,iy,iz)
          enddo
        enddo
      enddo

      return
      end
c=============================================================================
      subroutine printsumphiarray(nn,phi,s)
      use Parallel
      integer(ISZ):: nn
      real(kind=8):: phi(nn)
      character*(*):: s
      print*,my_index,s,nn,sum(phi)
      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine warptranspose(nx,ny,nzlocal,phi,nx_tran,ny_tran)
      use Subtimersf3d
      use Parallel
      use Transpose_work_space
      integer(ISZ):: nx,ny,nzlocal,nx_tran,ny_tran
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)

c 3-D transpose. The volume (0:nx-1,0:ny-1,0:nzlocal-1) is transposed.
c Transpose is NOT done in place but uses a separate work array.
c
c Algorithm assumes that all dimensions are powers of two and the number
c of processors is a power of two.

      integer(ISZ):: ip
      integer(MPIISZ):: dims,sizes(0:2),starts(0:2)
      integer(MPIISZ):: nlocal(0:2),nlocal_trnsps(0:2)
      include "mpif.h"
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: sendtypes(0:nzprocs-1),recvtypes(0:nzprocs-1)
      integer(MPIISZ):: sdispls(0:nzprocs-1),rdispls(0:nzprocs-1)
      integer(MPIISZ):: sendcounts(0:nzprocs-1),recvcounts(0:nzprocs-1)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Calculate what the dimensions of the grid will be after the
c     --- transpose, (0:nx_tran,0:ny_tran).
c     --- The work is divided up first along the y direction. Since
c     --- all dimensions and number of processors are powers of 2, ny/nzprocs
c     --- will always also be a power of 2 (or zero).
c     --- If ny > nzprocs, each processor is given ny/nzprocs lines along x.
c     --- If nx*ny > nzprocs, then the x lines are divided up as well.
c     --- Note that then the transposed data is stored as 2 'y' lines, each
c     --- nx*ny/nzprocs/2 long in x. This makes it easier to use the 2 at a
c     --- time operation of vpftz.
      if (ny > nzprocs) then
        nx_tran = nx
        ny_tran = ny/nzprocs
      else if (nx*ny > nzprocs) then
        nx_tran = nx*ny/nzprocs/2
        ny_tran = 2
      else
        call kaboom("warptranspose: the number of processors must be less than nx*ny")
        return
      endif

c     --- Create the work space that phi is transposed into.
      phi_trnspsnx = nx_tran
      phi_trnspsny = ny_tran
      phi_trnspsnz = nzlocal*nzprocs
      call gchange("Transpose_work_space",0)

c     --- Create an MPI type to refer to the distributed pieces of data that
c     --- are sent to each of the processors.
      dims = 3
      nlocal = (/nx+3,ny+3,nzlocal+1/)
      nlocal_trnsps = (/nx_tran+3,ny_tran+3,nzlocal*nzprocs+1/)
      sdispls = 0
      rdispls = 0
      sendcounts = 1
      recvcounts = 1

      do ip = 0,nzprocs-1

        if (ny < nzprocs) then
          starts(0) = 1 + (2*nx_tran)*mod(ip,nzprocs/ny)
          starts(1) = 1 + int(ip*ny/nzprocs)
          starts(2) = 0
          sizes = (/2*nx_tran,1,nzlocal/)
        else
          starts(0) = 1
          starts(1) = 1 + ny_tran*ip
          starts(2) = 0
          sizes = (/nx_tran,ny_tran,nzlocal/)
        endif
        if (izproc == nzprocs-1) sizes(2) = nzlocal + 1
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtypes(ip),mpierror)
        call MPI_TYPE_COMMIT(sendtypes(ip),mpierror)

        starts(0) = 1
        starts(1) = 1
        starts(2) = ip*nzlocal
        sizes = (/nx_tran,ny_tran,nzlocal/)
        if (ip == nzprocs-1) sizes(2) = nzlocal + 1
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal_trnsps,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                recvtypes(ip),mpierror)
        call MPI_TYPE_COMMIT(recvtypes(ip),mpierror)

      enddo

c     --- Do all of the communication at once.
      call MPI_ALLTOALLW(phi,sendcounts,sdispls,sendtypes,
     &                   phi_trnsps,recvcounts,rdispls,recvtypes,
     &                   MPI_COMM_WORLD,mpierror)

c     --- Free the data types
      do ip = 0,nzprocs-1
        call MPI_TYPE_FREE(sendtypes(ip),mpierror)
        call MPI_TYPE_FREE(recvtypes(ip),mpierror)
      enddo

c     --- Now the array phi_trnsps contains 3-D data with dimensions
c     --- (0:nx_tran,0:ny_tran,nzlocal*nzprocs).

c DONE!
!$OMP MASTER
      if (lf3dtimesubs) timetranspose = timetranspose + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine warptransposei(nx,ny,nzlocal,phi,nx_tran,ny_tran)
      use Subtimersf3d
      use Parallel
      use Transpose_work_space
      integer(ISZ):: nx,ny,nzlocal,nx_tran,ny_tran
      real(kind=8):: phi(-1:nx+1,-1:ny+1,0:nzlocal)

c Inverse 3-D transpose. The volume (0:nx-1,0:ny-1,0:nzlocal-1) is transposed.
c Transpose is NOT done in place but uses a work array.
c
c Algorithm assumes that all dimensions are powers of two and the number
c of processors is a power of two.

      integer(ISZ):: ip
      integer(MPIISZ):: dims,sizes(0:2),starts(0:2)
      integer(MPIISZ):: nlocal(0:2),nlocal_trnsps(0:2)
      include "mpif.h"
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: sendtypes(0:nzprocs-1),recvtypes(0:nzprocs-1)
      integer(MPIISZ):: sdispls(0:nzprocs-1),rdispls(0:nzprocs-1)
      integer(MPIISZ):: sendcounts(0:nzprocs-1),recvcounts(0:nzprocs-1)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Create an MPI type to refer to the distributed pieces of data that
c     --- are sent to each of the processors.
      dims = 3
      nlocal = (/nx+3,ny+3,nzlocal/)
      nlocal_trnsps = (/nx_tran+3,ny_tran+3,nzlocal*nzprocs/)
      sdispls = 0
      rdispls = 0
      sendcounts = 1
      recvcounts = 1

      do ip = 0,nzprocs-1

        starts(0) = 1
        starts(1) = 1
        starts(2) = ip*nzlocal
        sizes = (/nx_tran,ny_tran,nzlocal/)
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal_trnsps,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtypes(ip),mpierror)
        call MPI_TYPE_COMMIT(sendtypes(ip),mpierror)

        if (ny < nzprocs) then
          starts(0) = 1 + (2*nx_tran)*mod(ip,nzprocs/ny)
          starts(1) = 1 + int(ip*ny/nzprocs)
          starts(2) = 0
          sizes = (/2*nx_tran,1,nzlocal/)
        else
          starts(0) = 1
          starts(1) = 1 + ny_tran*ip
          starts(2) = 0
          sizes = (/nx_tran,ny_tran,nzlocal/)
        endif
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                recvtypes(ip),mpierror)
        call MPI_TYPE_COMMIT(recvtypes(ip),mpierror)

      enddo

c     --- Do all of the communication at once.
      call MPI_ALLTOALLW(phi_trnsps,sendcounts,sdispls,sendtypes,
     &                   phi,recvcounts,rdispls,recvtypes,
     &                   MPI_COMM_WORLD,mpierror)


c     --- Free the data types
      do ip = 0,nzprocs-1
        call MPI_TYPE_FREE(sendtypes(ip),mpierror)
        call MPI_TYPE_FREE(recvtypes(ip),mpierror)
      enddo

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
c=============================================================================
      subroutine lantzsolver(iwhich,a,kxsq,kysq,kzsq,attx,atty,attz,
     &                       filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,
     &                       l2symtry,l4symtry,bound0,boundnz,boundxy)
      use Subtimersf3d
      use Constant
      use Parallel
      use LantzSolverTemp
      integer(ISZ):: iwhich,nx,ny,nz
      real(kind=8):: a(-1:nx+1,-1:ny+1,-1:nz)
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
        nzlan = 2*nzprocs
        nzlocallan = nz
        if (ny >= nzprocs) then
          nxtranlan = nx
          nytranlan = ny/nzprocs
        else
          nxtranlan = nx*ny/nzprocs
          nytranlan = 1
        endif
        call gchange("LantzSolverTemp",0)
        lalloted = .true.
c       --- Initialize the k's
        call vpois3d(1,a(-1,-1,0),a(-1,-1,0),kxsq,kysq,kzsq,attx,atty,attz,
     &               filt,lx,ly,lz,nx,ny,nz,scrtch,xywork,zwork,0,
     &               l2symtry,l4symtry,bound0,boundnz,boundxy)
      endif

      if (iwhich == 1) return

      print*,0,a(4,5,4)
c     --- Apply FFT's to transverse planes.
      call vpois3d(12,a(-1,-1,0),a(-1,-1,0),kxsq,kysq,kzsq,attx,atty,attz,
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
      if (my_index == nzprocs-1) then
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
      call warptranspose(nx,ny,2,dtranlan,nx_tran,ny_tran)

c     --- Find location of transverse data relative to original dimensions.
      jx = mod(my_index*nx_tran,nx)
      jy = int(my_index*ny/real(nzprocs))
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
      do iz=0,2*nzprocs-1,2
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
      call generaltridiag(nx_tran,ny_tran-1,2*nzprocs,atranlan,btranlan,
     &                    ctranlan,dtranlan,ikxmin,ikymin,0)

c     --- Transpose back
      call warptransposei(nx,ny,2,dtranlan,nx_tran,ny_tran)

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
      call vpois3d(13,a(-1,-1,0),a(-1,-1,0),kxsq,kysq,kzsq,attx,atty,attz,
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
      real(kind=8):: a(-1:nx+1,-1:ny+1,-1:nz+1)
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
      call vpois3d(12,a(-1,-1,0),a(-1,-1,0),kxsq,kysq,kzsq,attx,atty,attz,
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
        call tridiag(nx,ny,nz,nz/2,a(-1,-1,0),norm,kxsq,kysq,kzsq,
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
      call vpois3d(13,a(-1,-1,0),a(-1,-1,0),kxsq,kysq,kzsq,attx,atty,attz,
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
      integer(MPIISZ):: isend
      integer(MPIISZ):: datatype,subcomm,ierror
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- First, create subgroups of processors
      if (nproc < nzprocs) then
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,int(my_index/nproc,MPIISZ),0,subcomm,
     &                      ierror)
      else
        call MPI_COMM_DUP(MPI_COMM_WORLD,subcomm,ierror)
      endif

c     --- First nstep > 1, create new data type
      if (nstep > 1) then
c       --- NOTE: thist coding is not complete
        call MPI_TYPE_VECTOR(int(nproc,MPIISZ),int(ndata,MPIISZ),int(ndata*nproc,MPIISZ),
     &                       MPI_DOUBLE_PRECISION,
     &                       datatype,ierror)
        call MPI_TYPE_COMMIT(datatype,mpierror)
      else
        datatype = MPI_DOUBLE_PRECISION
      endif

c     --- Gather the data
      isend = mod(my_index,nproc) + 1
      call MPI_ALLGATHER(data(1:ndata,isend,1),int(ndata,MPIISZ),
     &                   MPI_DOUBLE_PRECISION,
     &                   data,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,subcomm,ierror)

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
      subroutine mgdividenz(finedecomp,coarsedecomp,
     &                      nx,ny,nz,nxcoarse,nycoarse,nzcoarse,mgscale)
      use Subtimersf3d
      use Multigrid3d,Only:mgscaleserial
      use Decompositionmodule
      type(Decomposition):: finedecomp,coarsedecomp
      integer(ISZ):: nx,ny,nz,nxcoarse,nycoarse,nzcoarse
      real(kind=8):: mgscale

      call mgdividenz_work(finedecomp%nxprocs,finedecomp%ix,finedecomp%nx,
     &                     coarsedecomp%ix,coarsedecomp%nx,nx,nxcoarse,mgscale)
      call mgdividenz_work(finedecomp%nyprocs,finedecomp%iy,finedecomp%ny,
     &                     coarsedecomp%iy,coarsedecomp%ny,ny,nycoarse,mgscale)
      call mgdividenz_work(finedecomp%nzprocs,finedecomp%iz,finedecomp%nz,
     &                     coarsedecomp%iz,coarsedecomp%nz,nz,nzcoarse,mgscale)

      return
      end
c=============================================================================
      subroutine mgdividenz_work(nzprocs,izfsdecomp,nzfsdecomp,
     &                           izfsdecompc,nzfsdecompc,nz,nzcoarse,mgscale)
      use Subtimersf3d
      use Multigrid3d,Only:mgscaleserial
      integer(ISZ):: nzprocs,nz,nzcoarse
      integer(ISZ):: izfsdecomp(0:nzprocs-1),nzfsdecomp(0:nzprocs-1)
      integer(ISZ):: izfsdecompc(0:nzprocs-1),nzfsdecompc(0:nzprocs-1)
      real(kind=8):: mgscale

c Divides nz for each processor by two. For odd numbers, the starting end
c is rounded down, and the upper end up. For each processor, if nz ends up
c being equal to 1, then change it to 2. For the last processor, if nz is 1,
c also, subtract 1 from izfsdecomp.

      integer(ISZ):: is,nzmin = 2
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (mgscale < mgscaleserial .and. nz > 0) then

        do is=0,nzprocs-1

c         --- First divide nz by 2. izfsdecomp is rounded down, and nzfsdecomp
c         --- is rounded up.
          izfsdecompc(is) = int(floor(1.*izfsdecomp(is)*nzcoarse/nz))
          nzfsdecompc(is) = int(ceiling(1.*(izfsdecomp(is)+nzfsdecomp(is))*
     &                                  nzcoarse/nz - izfsdecompc(is)))

c         --- Then check that no nz are less than nzmin.
c         --- Note that if nzfsdecomp is already < nzmin, leave it as is.
c         --- When nzfsdecomp < nzmin, in some cases setting nzfsdecompc to
c         --- nzmin will cause the coarsened grid to extend too far beyond the
c         --- finer grid's gaurd cell, so it would need data from a neighboring
c         --- cell, where it would normally be able to use data from the gaurd
c         --- cells. That is a rare occurance so it is simpler to avoid the
c         --- problem then to write the code to exchange the data for that case.
          if (nzfsdecompc(is) < nzmin .and. nzfsdecompc(is) < nzfsdecomp(is)) then
            nzfsdecompc(is) = nzfsdecompc(is) + 1
            if (izfsdecompc(is)+nzfsdecompc(is) > nzcoarse) then
              izfsdecompc(is) = izfsdecompc(is) - 1
              nzfsdecompc(is) = nzcoarse - izfsdecompc(is)
            endif
          endif

        enddo

      else

        if (mgscale > mgscaleserial) then
          call kaboom("Error: mgscaleserial is not supported, don't set it!")
          return
        endif

c       --- When the volume of the coarsened grid gets small compared to
c       --- the finest level, then expand the domain of each processor to
c       --- cover the full system, to avoid passing around lots of small
c       --- messages.
        izfsdecompc = 0
        nzfsdecompc = nzcoarse

      endif

!$OMP MASTER
      if (lf3dtimesubs) timemgdividenz = timemgdividenz + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mggetexchangepes(fsdecomp,bounds,nx,ny,nz,
     &                            lxoffset,rxoffset,
     &                            lyoffset,ryoffset,
     &                            lzoffset,rzoffset,
     &                            whosendingdown,whosendingup)
      use Decompositionmodule
      type(Decomposition):: fsdecomp
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: nx,ny,nz
      integer(ISZ):: lxoffset(fsdecomp%nxprocs),rxoffset(fsdecomp%nxprocs)
      integer(ISZ):: lyoffset(fsdecomp%nyprocs),ryoffset(fsdecomp%nyprocs)
      integer(ISZ):: lzoffset(fsdecomp%nzprocs),rzoffset(fsdecomp%nzprocs)
      type(Decomposition):: whosendingdown
      type(Decomposition):: whosendingup

      integer(ISZ):: ix,iy,iz
      integer(ISZ):: convertindextoproc

      whosendingdown%nxprocs = fsdecomp%nxprocs
      whosendingdown%nyprocs = fsdecomp%nyprocs
      whosendingdown%nzprocs = fsdecomp%nzprocs
      allocate(whosendingdown%nx(0:fsdecomp%nxprocs-1))
      allocate(whosendingdown%ny(0:fsdecomp%nyprocs-1))
      allocate(whosendingdown%nz(0:fsdecomp%nzprocs-1))
      allocate(whosendingdown%ix(0:fsdecomp%nxprocs-1))
      allocate(whosendingdown%iy(0:fsdecomp%nyprocs-1))
      allocate(whosendingdown%iz(0:fsdecomp%nzprocs-1))

      whosendingup%nxprocs = fsdecomp%nxprocs
      whosendingup%nyprocs = fsdecomp%nyprocs
      whosendingup%nzprocs = fsdecomp%nzprocs
      allocate(whosendingup%nx(0:fsdecomp%nxprocs-1))
      allocate(whosendingup%ny(0:fsdecomp%nyprocs-1))
      allocate(whosendingup%nz(0:fsdecomp%nzprocs-1))
      allocate(whosendingup%ix(0:fsdecomp%nxprocs-1))
      allocate(whosendingup%iy(0:fsdecomp%nyprocs-1))
      allocate(whosendingup%iz(0:fsdecomp%nzprocs-1))

c     --- The Decomposition type is used as a convenient place to put the
c     --- information about where data is to be sent. The nx, ny and nz will
c     --- hold the processor ID that data needs to be sent to, and the ix, iy
c     --- and iz will hold the local starting position of the data. Using
c     --- Decomposition avoids having six separate arrays floating around.

c     --- First, get the process IDs and ix's relative to the ix line.
      call mggetexchangepes_work(fsdecomp%nxprocs,fsdecomp%ix,fsdecomp%nx,
     &                           bounds(0:1),nx,lxoffset,rxoffset,
     &                           whosendingdown%nx,whosendingdown%ix,
     &                           whosendingup%nx,whosendingup%ix)

c     --- Convert the process IDs to be global.
c     iy = fsdecomp%iyproc
c     iz = fsdecomp%izproc
c     do ix=0,fsdecomp%nxprocs-1
c       if (whosendingdown%nx(ix) > -1) then
c         ix = whosendingdown%nx(ix)
c         whosendingdown%nx(ix) = convertindextoproc(ix,iy,iz,
c    &                                               fsdecomp%nxprocs,
c    &                                               fsdecomp%nyprocs,
c    &                                               fsdecomp%nzprocs)
c       endif
c       if (whosendingup%nx(ix) > -1) then
c         ix = whosendingup%nx(ix)
c         whosendingup%nx(ix) = convertindextoproc(ix,iy,iz,
c    &                                             fsdecomp%nxprocs,
c    &                                             fsdecomp%nyprocs,
c    &                                             fsdecomp%nzprocs)
c       endif
c     enddo

c     --- First, get the process IDs and iy's relative to the iy line.
      call mggetexchangepes_work(fsdecomp%nyprocs,fsdecomp%iy,fsdecomp%ny,
     &                           bounds(2:3),ny,lyoffset,ryoffset,
     &                           whosendingdown%ny,whosendingdown%iy,
     &                           whosendingup%ny,whosendingup%iy)

c     --- Convert the process IDs to be global.
c     ix = fsdecomp%ixproc
c     iz = fsdecomp%izproc
c     do iy=0,fsdecomp%nyprocs-1
c       if (whosendingdown%ny(iy) > -1) then
c         iy = whosendingdown%ny(iy)
c         whosendingdown%ny(iy) = convertindextoproc(ix,iy,iz,
c    &                                               fsdecomp%nxprocs,
c    &                                               fsdecomp%nyprocs,
c    &                                               fsdecomp%nzprocs)
c       endif
c       if (whosendingup%nx(ix) > -1) then
c         iy = whosendingup%ny(iy)
c         whosendingup%ny(iy) = convertindextoproc(ix,iy,iz,
c    &                                             fsdecomp%nxprocs,
c    &                                             fsdecomp%nyprocs,
c    &                                             fsdecomp%nzprocs)
c       endif
c     enddo

c     --- First, get the process IDs and iz's relative to the iz line.
      call mggetexchangepes_work(fsdecomp%nzprocs,fsdecomp%iz,fsdecomp%nz,
     &                           bounds(4:5),nz,lzoffset,rzoffset,
     &                           whosendingdown%nz,whosendingdown%iz,
     &                           whosendingup%nz,whosendingup%iz)

c     --- Convert the process IDs to be global.
c     ix = fsdecomp%ixproc
c     iy = fsdecomp%iyproc
c     do iz=0,fsdecomp%nzprocs-1
c       if (whosendingdown%nz(iz) > -1) then
c         iz = whosendingdown%nz(iz)
c         whosendingdown%nz(iz) = convertindextoproc(ix,iy,iz,
c    &                                               fsdecomp%nxprocs,
c    &                                               fsdecomp%nyprocs,
c    &                                               fsdecomp%nzprocs)
c       endif
c       if (whosendingup%nx(ix) > -1) then
c         iz = whosendingup%nz(iz)
c         whosendingup%nz(iz) = convertindextoproc(ix,iy,iz,
c    &                                             fsdecomp%nxprocs,
c    &                                             fsdecomp%nyprocs,
c    &                                             fsdecomp%nzprocs)
c       endif
c     enddo

      return
      end
c=============================================================================
      subroutine mggetexchangepes_work(nzprocs,izfsdecomp,nzfsdecomp,
     &                                 bounds,nz,lzoffset,rzoffset,
     &                                 whosendingleft,izsendingleft,
     &                                 whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nzprocs
      integer(ISZ):: izfsdecomp(0:nzprocs-1),nzfsdecomp(0:nzprocs-1)
      integer(ISZ):: bounds(0:1),nz
      integer(ISZ):: lzoffset(0:nzprocs-1),rzoffset(0:nzprocs-1)
      integer(ISZ):: whosendingleft(0:nzprocs-1)
      integer(ISZ):: whosendingright(0:nzprocs-1)
      integer(ISZ):: izsendingleft(0:nzprocs-1)
      integer(ISZ):: izsendingright(0:nzprocs-1)

c Returns list of processors to the left and right to exchange data with.
c The algorithm tries to share the message passing among the
c processors holding the same data. For each processor, it finds the other
c processors which have the data it needs. It then picks the one which
c is sending out the fewest messages so far.

      integer(ISZ):: numsendingleft(0:nzprocs-1)
      integer(ISZ):: numsendingright(0:nzprocs-1)
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
c     --- Loop over all processors since this must be a global process.
      do idx=0,nzprocs-1

c       --- If this process covers the full extent of the mesh on this and the
c       --- next finer mesh, then it doesn't need any data sent to it.
        if (izfsdecomp(idx) == 0 .and. izfsdecomp(idx)+nzfsdecomp(idx) == nz
     &      .and. rzoffset(idx) == 0.)
     &    cycle

c       --- Skip processors which have the same right mesh extent (and
c       --- therefore don't have the data needed for this processor.
        js = idx + 1
        do while (js < nzprocs .and.
     &         izfsdecomp(js)+nzfsdecomp(js)==izfsdecomp(idx)+nzfsdecomp(idx))
          js = js + 1
        enddo

c       --- Processors whose extent extends to the right edge are treated
c       --- differently. They are sent data from the left with periodic
c       --- boundary conditions. Otherwise, all receive data from the last
c       --- processor.
c       --- Do all but that one first.
        if (js < nzprocs) then

c         --- Find all processors whose data extent overlaps the right
c         --- end of this processor.
          ie = js
          do while (ie < nzprocs .and.
     &              izfsdecomp(ie) < izfsdecomp(idx)+nzfsdecomp(idx))
            ie = ie + 1
          enddo
          ie = ie - 1

        elseif (bounds(1) == 2) then
c         --- Now, treat the right most processor.

c         --- Find all processors who data extent overlaps the right
c         --- end of this processor, wrapping around to the left edge.
          js = 0
          ie = js
          do while (ie < nzprocs .and. izfsdecomp(ie) == 0 .and.
     &              lzoffset(ie) == 0)
            ie = ie + 1
          enddo
          ie = ie - 1
        else

c         --- Some of the processors which extend to the right hand edge will
c         --- need data there if they did not extend to the right edge
c         --- at a the next finer level. These processors are known by
c         --- having rzoffset > 0.
          if (idx < nzprocs-1 .and. rzoffset(idx) > 0.) then
            ie = nzprocs - 1
            js = ie - 1
            do while (js > idx .and. rzoffset(js) == 0. .and.
     &          izfsdecomp(js)+nzfsdecomp(js)==izfsdecomp(idx)+nzfsdecomp(idx))
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
          izsendingleft(idx) = izfsdecomp(idx) + nzfsdecomp(idx) -
     &                         izfsdecomp(whosendingleft(idx))
        else
c         --- For periodic boundary conditions, this always has the same value.
          izsendingleft(idx) = 1
        endif

c       --- Increment the number of sends that processor makes.
        numsendingleft(whosendingleft(idx)) =
     &                       numsendingleft(whosendingleft(idx)) + 1
      enddo

c     --- For each processor, find a processor on the left which has the
c     --- data needed.
c     --- Loop over all processors since this most be a global process.
      do idx=0,nzprocs-1

c       --- If this process covers the full extent on this and the
c       --- next finer mesh, then it doesn't need any data sent to it.
        if (izfsdecomp(idx) == 0 .and. izfsdecomp(idx)+nzfsdecomp(idx) == nz
     &      .and. lzoffset(idx) == 0)
     &    cycle

c       --- Skip processors which have the same left mesh extent (and
c       --- therefore don't have the data needed for this processor.
        js = idx - 1
        do while (js >= 0 .and. izfsdecomp(js) == izfsdecomp(idx))
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
     &              izfsdecomp(idx) < izfsdecomp(ie)+nzfsdecomp(ie))
            ie = ie - 1
          enddo
          ie = ie + 1

        elseif (bounds(0) == 2) then
c         --- Now, treat the left most processor.

c         --- Find all processors who data extent overlaps the left
c         --- end of this processor, wrapping around to the left edge.
          ie = nzprocs-1
          js = ie
          do while (js >= 0 .and. izfsdecomp(js)+nzfsdecomp(js) == nz .and.
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
            do while (ie < idx .and. lzoffset(ie) == 0 .and.
     &                izfsdecomp(ie) == 0)
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
          izsendingright(idx)=izfsdecomp(idx)-izfsdecomp(whosendingright(idx))

        else
c         --- For periodic boundary conditions, this always has the same value.
          izsendingright(idx) = nzfsdecomp(whosendingright(idx)) - 1
        endif

c       --- Increment the number of sends that processor makes.
        numsendingright(whosendingright(idx)) =
     &              numsendingright(whosendingright(idx)) + 1
      enddo


!$OMP MASTER
      if (lf3dtimesubs) timemggetexchangepes = timemggetexchangepes + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c===========================================================================
      subroutine printarray3d(nx,ny,nz,iz0,nz0,izs,a,my_index,nzprocs,text,
     &                        mglevel)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,iz0,nz0,izs,my_index,nzprocs,mglevel
      real(kind=8):: a(-1:nx+1,-1:ny+1,iz0:nz)
      character(*):: text

c Specialized routine for printing out an array in parallel, keeping the same
c ordering as would be in serial. Primarily used for debugging purposes.
c It is maintained in case it is needed for future use since there is some subtlty
c to getting it correct.

      integer(MPIISZ):: ix,iy,iz,lastz,firstz
      include "mpif.h"
      integer(MPIISZ):: mpierror,mpistatus(MPI_STATUS_SIZE)
      integer(MPIISZ):: messid = 999
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
      if (my_index == nzprocs-1) lastz = nz

c Each processor prints out its data and then sends a message to the next so it
c can print its data.
      if (my_index > 0) then
        call MPI_RECV(firstz,int(1,MPIISZ),MPI_INTEGER,int(my_index-1,MPIISZ),messid,
     &                MPI_COMM_WORLD,mpistatus,mpierror)
        firstz = firstz - izs
      else
        firstz = 0
      endif

c The file needs to be reopened each time since there is no coordination in how
c the processors deal with files.
      if (nzprocs == 16) then
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

      if (my_index < nzprocs - 1) then
        call MPI_SEND(izs+lastz+1,int(1,MPIISZ),MPI_INTEGER,int(my_index+1,MPIISZ),messid,
     &                MPI_COMM_WORLD,mpierror)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timeprintarray3d = timeprintarray3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================

c===========================================================================
c===========================================================================
c===========================================================================
c===========================================================================
      subroutine getbforparticles3d(bfield,bfieldp)
      use BFieldGridTypemodule
      use Subtimersf3d
      use Parallel
      type(BFieldGridType):: bfield,bfieldp

c Get the B for the extent where the particles are. This gets B from
c neighboring processors, and at the very least gets B in the guard planes,
c iz=-1 and +1.

      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      call getphipforparticles3d_parallel(3,
     &             bfield%nxlocal,bfield%nylocal,bfield%nzlocal,bfield%b,
     &             bfieldp%nxlocal,bfieldp%nylocal,bfieldp%nzlocal,bfieldp%b,
     &             0,0,0,fsdecomp,ppdecomp)

!$OMP MASTER
      if (lf3dtimesubs) timegetbforparticles3d = timegetbforparticles3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
      subroutine getaforfields3d(bfield)
      use BFieldGridTypemodule
      use Decompositionmodule
      use Subtimersf3d
      use Parallel
      type(BFieldGridType):: bfield

c Get the A for the full extent where the fields are. This gets A from
c neighboring processors, and at the very least gets A in the guard planes,
c iz=-1 and +1.

      integer(ISZ):: nx,ny,nz
      integer(ISZ):: nxlocal,nylocal,nzlocal
      integer(ISZ):: localbounds(0:5)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      nx = bfield%nx
      ny = bfield%ny
      nz = bfield%nz
      nxlocal = bfield%nxlocal
      nylocal = bfield%nylocal
      nzlocal = bfield%nzlocal

      localbounds = bfield%bounds
      if (fsdecomp%ix(fsdecomp%ixproc) > 0)          localbounds(0) = -1
      if (fsdecomp%ix(fsdecomp%ixproc)+nxlocal < nx) localbounds(1) = -1
      if (fsdecomp%iy(fsdecomp%iyproc) > 0)          localbounds(2) = -1
      if (fsdecomp%iy(fsdecomp%iyproc)+nylocal < ny) localbounds(3) = -1
      if (fsdecomp%iz(fsdecomp%izproc) > 0)          localbounds(4) = -1
      if (fsdecomp%iz(fsdecomp%izproc)+nzlocal < nz) localbounds(5) = -1

      call mgexchange_phi_periodic(3,nxlocal,nylocal,nzlocal,bfield%a,
     &                             1,1,1,0,0,localbounds,fsdecomp)
      call mgexchange_phi(3,nxlocal,nylocal,nzlocal,bfield%a,1,1,1,
     &                    -1,0,fsdecomp)

!$OMP MASTER
      if (lf3dtimesubs) timegetaforfields3d = timegetaforfields3d + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c=============================================================================
      module messagedatatypemodule
c       --- This module holds a data type which can store all of the data
c       --- needed for exchange phi at the boundaries. This data is cached
c       --- since it is expensive to calculate and to setup the MPI
c       --- type definitions. mgexcange_phi will be called with four
c       --- different sets of values for each level of coarsening. The size
c       --- of the cache will then allow 50 levels, which should be sufficient.
        SAVE
        type messagedatatype
          integer(ISZ):: id = 0
          integer(MPIISZ),pointer:: xrecvtypes(:),xsendtypes(:)
          integer(MPIISZ),pointer:: xsendcounts(:),xrecvcounts(:)
          integer(MPIISZ),pointer:: xdispls(:)
          integer(MPIISZ),pointer:: yrecvtypes(:),ysendtypes(:)
          integer(MPIISZ),pointer:: ysendcounts(:),yrecvcounts(:)
          integer(MPIISZ),pointer:: ydispls(:)
          integer(MPIISZ),pointer:: zrecvtypes(:),zsendtypes(:)
          integer(MPIISZ),pointer:: zsendcounts(:),zrecvcounts(:)
          integer(MPIISZ),pointer:: zdispls(:)
        end type messagedatatype
        integer:: icache = 0
        type(messagedatatype):: messagedatacache(200)
      contains
        subroutine getmessagedata(nc,nxlocal,nylocal,nzlocal,
     &                            delx,dely,delz,ils,ile,ius,iue,fsdecomp,
     &                            messagedata)
        use Subtimersf3d
        use Decompositionmodule
        integer(ISZ):: nc,nxlocal,nylocal,nzlocal,ils,ile,ius,iue
        integer(ISZ):: delx,dely,delz
        type(Decomposition):: fsdecomp
        type(messagedatatype):: messagedata

c Calculate the data needed to do the exchange of potential on the boundary
c and/or in the gaurd cells. The data is cached, and if already calculated,
c the cached data is returned, giving a substantial time savings.

c ils, ile are the start and end of the range of slices at the lower edges
c          of the grid that need to be obtained from neighbors. The values are
c          relative to 0.
c ius, iue are the start and end of the range of slices at the upper edges
c          of the grid that need to be obtained from neighbors. The values are
c          relative to the upper index, nx, ny or nz.

        integer(ISZ):: ic,id
        integer(ISZ):: npx,npy,npz
        integer(ISZ):: px,py,pz
        integer(ISZ):: me_ix,me_iy,me_iz
        integer(ISZ):: me_nx,me_ny,me_nz
        integer(ISZ):: pe_ix,pe_iy,pe_iz
        integer(ISZ):: pe_nx,pe_ny,pe_nz
        integer(ISZ):: ixls,iyls,izls,ixle,iyle,izle
        integer(ISZ):: ixus,iyus,izus,ixue,iyue,izue
        integer(ISZ):: me_ixlgood,me_ixugood
        integer(ISZ):: me_iylgood,me_iyugood
        integer(ISZ):: me_izlgood,me_izugood
        integer(ISZ):: pe_ixlgood,pe_ixugood
        integer(ISZ):: pe_iylgood,pe_iyugood
        integer(ISZ):: pe_izlgood,pe_izugood
        integer(ISZ):: i1global,i2global
        integer(MPIISZ):: iproc,nn
        include "mpif.h"
        integer(MPIISZ):: dims,nlocal(-1:2),starts(-1:2),sizes(-1:2)
        integer(MPIISZ):: mpierror
        real(kind=8):: substarttime,wtime,s1,s2
        if (lf3dtimesubs) substarttime = wtime()
        s1 = wtime()

c       --- Compute a hash of the input parameters.
c       --- Note that this expression is ad hoc and makes no gaurantee
c       --- of uniqueness.
        id = + (delx + 3)
     &       + (dely + 3)*10
     &       + (delz + 3)*100
     &       + (ils + 5)*1000
     &       + (ile + 5)*10000
     &       + (ius + 5)*100000
     &       + (iue + 5)*1000000
     &       + sum(fsdecomp%ix**1)
     &       + sum(fsdecomp%iy**2)
     &       + sum(fsdecomp%iz**3)
     &       + sum(fsdecomp%nx**3)
     &       + sum(fsdecomp%ny**2)
     &       + sum(fsdecomp%nz**1)

c       --- Search through the cached data to see if this data set has
c       --- already been calculated.
c       --- Note, to turn caching off, only this loop needs to be commented out
        do ic = 1,SIZE(messagedatacache)

          if (messagedatacache(ic)%id == id) then

c           --- The id is the same, then there's no need to recalculate.

c           --- Setup the array pointers.
            messagedata%xsendtypes => messagedatacache(ic)%xsendtypes
            messagedata%xrecvtypes => messagedatacache(ic)%xrecvtypes
            messagedata%xsendcounts => messagedatacache(ic)%xsendcounts
            messagedata%xrecvcounts => messagedatacache(ic)%xrecvcounts
            messagedata%xdispls => messagedatacache(ic)%xdispls
            messagedata%ysendtypes => messagedatacache(ic)%ysendtypes
            messagedata%yrecvtypes => messagedatacache(ic)%yrecvtypes
            messagedata%ysendcounts => messagedatacache(ic)%ysendcounts
            messagedata%yrecvcounts => messagedatacache(ic)%yrecvcounts
            messagedata%ydispls => messagedatacache(ic)%ydispls
            messagedata%zsendtypes => messagedatacache(ic)%zsendtypes
            messagedata%zrecvtypes => messagedatacache(ic)%zrecvtypes
            messagedata%zsendcounts => messagedatacache(ic)%zsendcounts
            messagedata%zrecvcounts => messagedatacache(ic)%zrecvcounts
            messagedata%zdispls => messagedatacache(ic)%zdispls

            return
          endif

        enddo

c       --- A new data set is needed. Increment the counter and put the data
c       --- in the next available spot.
        icache = mod(icache,SIZE(messagedatacache)) + 1

c       --- If the id is non zero, then the previous data in the cache needs
c       --- to be cleaned up.
        if (messagedatacache(icache)%id .ne. 0) then
          call freemessagedata(messagedatacache(icache))
        endif

        messagedatacache(icache)%id = id

        npx = fsdecomp%nxprocs
        npy = fsdecomp%nyprocs
        npz = fsdecomp%nzprocs
        allocate(messagedatacache(icache)%xsendtypes(0:npx-1))
        allocate(messagedatacache(icache)%xrecvtypes(0:npx-1))
        allocate(messagedatacache(icache)%xsendcounts(0:npx-1))
        allocate(messagedatacache(icache)%xrecvcounts(0:npx-1))
        allocate(messagedatacache(icache)%xdispls(0:npx-1))
        allocate(messagedatacache(icache)%ysendtypes(0:npy-1))
        allocate(messagedatacache(icache)%yrecvtypes(0:npy-1))
        allocate(messagedatacache(icache)%ysendcounts(0:npy-1))
        allocate(messagedatacache(icache)%yrecvcounts(0:npy-1))
        allocate(messagedatacache(icache)%ydispls(0:npy-1))
        allocate(messagedatacache(icache)%zsendtypes(0:npz-1))
        allocate(messagedatacache(icache)%zrecvtypes(0:npz-1))
        allocate(messagedatacache(icache)%zsendcounts(0:npz-1))
        allocate(messagedatacache(icache)%zrecvcounts(0:npz-1))
        allocate(messagedatacache(icache)%zdispls(0:npz-1))

        messagedata%xsendtypes => messagedatacache(icache)%xsendtypes
        messagedata%xrecvtypes => messagedatacache(icache)%xrecvtypes
        messagedata%xsendcounts => messagedatacache(icache)%xsendcounts
        messagedata%xrecvcounts => messagedatacache(icache)%xrecvcounts
        messagedata%xdispls => messagedatacache(icache)%xdispls
        messagedata%ysendtypes => messagedatacache(icache)%ysendtypes
        messagedata%yrecvtypes => messagedatacache(icache)%yrecvtypes
        messagedata%ysendcounts => messagedatacache(icache)%ysendcounts
        messagedata%yrecvcounts => messagedatacache(icache)%yrecvcounts
        messagedata%ydispls => messagedatacache(icache)%ydispls
        messagedata%zsendtypes => messagedatacache(icache)%zsendtypes
        messagedata%zrecvtypes => messagedatacache(icache)%zrecvtypes
        messagedata%zsendcounts => messagedatacache(icache)%zsendcounts
        messagedata%zrecvcounts => messagedatacache(icache)%zrecvcounts
        messagedata%zdispls => messagedatacache(icache)%zdispls

c       --- Now, do the actual calculation.

        dims = 4
        nlocal = (/nc,1+nxlocal+2*delx,1+nylocal+2*dely,1+nzlocal+2*delz/)
        messagedata%xdispls = 0
        messagedata%ydispls = 0
        messagedata%zdispls = 0

c     --- First, calculate what data needs to be sent and received.

c     --- Send and recv all of the first dimension.
        starts(-1) = 0
        sizes(-1) = nc

c     --- Data to be sent
        me_ix = fsdecomp%ix(fsdecomp%ixproc)
        me_nx = fsdecomp%nx(fsdecomp%ixproc)
        me_iy = fsdecomp%iy(fsdecomp%iyproc)
        me_ny = fsdecomp%ny(fsdecomp%iyproc)
        me_iz = fsdecomp%iz(fsdecomp%izproc)
        me_nz = fsdecomp%nz(fsdecomp%izproc)

        ixls = -max(-delx,ils)
        ixle = -ile
        iyls = -max(-dely,ils)
        iyle = -ile
        izls = -max(-delz,ils)
        izle = -ile

        ixus = min(delx,ius)
        ixue = iue
        iyus = min(dely,ius)
        iyue = iue
        izus = min(delz,ius)
        izue = iue

c     --- Get the lower and upper range of the good data. Good data means
c     --- data that was calculated on the processor and is available to send
c     --- to other processors that need it. The default value is the range of
c     --- data out to the data that is needed, set by the ixle etc input
c     --- arguments. If the processor sits on the edge of the grid, then
c     --- include all of the data in the guard cells as well.
        me_ixlgood = me_ix + 1 - ixle
        me_ixugood = me_ix + me_nx - 1 + ixus
        me_iylgood = me_iy + 1 - iyle
        me_iyugood = me_iy + me_ny - 1 + iyus
        me_izlgood = me_iz + 1 - izle
        me_izugood = me_iz + me_nz - 1 + izus
        if (me_ix == 0) me_ixlgood = -delx
        if (me_ix + me_nx == fsdecomp%nxglobal) me_ixugood = me_ix + me_nx + delx
        if (me_iy == 0) me_iylgood = -dely
        if (me_iy + me_ny == fsdecomp%nyglobal) me_iyugood = me_iy + me_ny + dely
        if (me_iz == 0) me_izlgood = -delz
        if (me_iz + me_nz == fsdecomp%nzglobal) me_izugood = me_iz + me_nz + delz

        do px=0,fsdecomp%nxprocs-1

          starts(1) = 0
          sizes(1) = me_ny + 2*dely + 1
          starts(2) = 0
          sizes(2) = me_nz + 2*delz + 1

          pe_ix = fsdecomp%ix(px)
          pe_nx = fsdecomp%nx(px)

          if (px > fsdecomp%ixproc) then
c         --- Sending to the lower guard cells on the right
            i1global = max(pe_ix - ixls,me_ixlgood)
            i2global = min(pe_ix - ixle,me_ixugood)
            if (px > fsdecomp%ixproc+1) then
              i2global = fsdecomp%ix(fsdecomp%ixproc+1) - ixle
            endif
          else if (px < fsdecomp%ixproc) then
c         --- Sending to the upper guard cells on the left
            i1global = max(pe_ix + pe_nx + ixus,me_ixlgood)
            i2global = min(pe_ix + pe_nx + ixue,me_ixugood)
            if (px < fsdecomp%ixproc-1) then
              i1global = fsdecomp%ix(fsdecomp%ixproc-1) +
     &                   fsdecomp%nx(fsdecomp%ixproc-1) + ixus
            endif
          else
c         --- Send nothing to self
            i1global = 0
            i2global = -1
          endif
          starts(0) = i1global - me_ix + delx
          sizes(0) = i2global - i1global + 1

c       --- Create the data types for data that is to be sent
          if (ALL(sizes > 0)) then
            call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                  MPI_ORDER_FORTRAN,
     &                                  MPI_DOUBLE_PRECISION,
     &                                  messagedata%xsendtypes(px),mpierror)
            call MPI_TYPE_COMMIT(messagedata%xsendtypes(px),mpierror)
            messagedata%xsendcounts(px) = 1
          else
            messagedata%xsendtypes(px) = MPI_DOUBLE_PRECISION
            messagedata%xsendcounts(px) = 0
          endif

          pe_ixlgood = pe_ix + 1 - ixle
          pe_ixugood = pe_ix + pe_nx - 1 + ixus
          if (pe_ix == 0) pe_ixlgood = -delx
          if (pe_ix+pe_nx == fsdecomp%nxglobal) pe_ixugood = pe_ix + pe_nx + delx

          if (px < fsdecomp%ixproc) then
c         --- Receiving to the lower guard cells from the left
            i1global = max(me_ix - ixls,pe_ixlgood)
            i2global = min(me_ix - ixle,pe_ixugood)
            if (px < fsdecomp%ixproc-1) then
              i2global = fsdecomp%ix(px+1) - ixle
            endif
          else if (px > fsdecomp%ixproc) then
c         --- Receiving to the upper guard cells from the right
            i1global = max(me_ix + me_nx + ixus,pe_ixlgood)
            i2global = min(me_ix + me_nx + ixue,pe_ixugood)
            if (px > fsdecomp%ixproc+1) then
              i1global = fsdecomp%ix(px-1) + fsdecomp%nx(px-1) + ixus
            endif
          else
c         --- Receive nothing from self
            i1global = 0
            i2global = -1
          endif
          starts(0) = i1global - me_ix + delx
          sizes(0) = i2global - i1global + 1

c       --- Create the data types for data that is to be received
          if (ALL(sizes > 0)) then
            call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                  MPI_ORDER_FORTRAN,
     &                                  MPI_DOUBLE_PRECISION,
     &                                  messagedata%xrecvtypes(px),mpierror)
            call MPI_TYPE_COMMIT(messagedata%xrecvtypes(px),mpierror)
            messagedata%xrecvcounts(px) = 1
          else
            messagedata%xrecvtypes(px) = MPI_DOUBLE_PRECISION
            messagedata%xrecvcounts(px) = 0
          endif

        enddo

        do py=0,fsdecomp%nyprocs-1

          starts(0) = 0
          sizes(0) = me_nx + 2*delx + 1
          starts(2) = 0
          sizes(2) = me_nz + 2*delz + 1

          pe_iy = fsdecomp%iy(py)
          pe_ny = fsdecomp%ny(py)

          if (py > fsdecomp%iyproc) then
c         --- Sending to the lower guard cells on the right
            i1global = max(pe_iy - iyls,me_iylgood)
            i2global = min(pe_iy - iyle,me_iyugood)
            if (py > fsdecomp%iyproc+1) then
              i2global = fsdecomp%iy(fsdecomp%iyproc+1) - iyle
            endif
          else if (py < fsdecomp%iyproc) then
c         --- Sending to the upper guard cells on the left
            i1global = max(pe_iy + pe_ny + iyus,me_iylgood)
            i2global = min(pe_iy + pe_ny + iyue,me_iyugood)
            if (py < fsdecomp%iyproc-1) then
              i1global = fsdecomp%iy(fsdecomp%iyproc-1) +
     &                   fsdecomp%ny(fsdecomp%iyproc-1) + iyus
            endif
          else
c           --- Send nothing to self
            i1global = 0
            i2global = -1
          endif
          starts(1) = i1global - me_iy + dely
          sizes(1) = i2global - i1global + 1

c         --- Create the data types for data that is to be sent
          if (ALL(sizes > 0)) then
            call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                    MPI_ORDER_FORTRAN,
     &                                    MPI_DOUBLE_PRECISION,
     &                                    messagedata%ysendtypes(py),mpierror)
            call MPI_TYPE_COMMIT(messagedata%ysendtypes(py),mpierror)
            messagedata%ysendcounts(py) = 1
          else
            messagedata%ysendtypes(py) = MPI_DOUBLE_PRECISION
            messagedata%ysendcounts(py) = 0
          endif

          pe_iylgood = pe_iy + 1 - iyle
          pe_iyugood = pe_iy + pe_ny - 1 + iyus
          if (pe_iy == 0) pe_iylgood = -dely
          if (pe_iy+pe_ny == fsdecomp%nyglobal) pe_iyugood = pe_iy + pe_ny + dely

          if (py < fsdecomp%iyproc) then
c           --- Receiving to the lower guard cells from the left
            i1global = max(me_iy - iyls,pe_iylgood)
            i2global = min(me_iy - iyle,pe_iyugood)
            if (py < fsdecomp%iyproc-1) then
              i2global = fsdecomp%iy(py+1) - iyle
            endif
          else if (py > fsdecomp%iyproc) then
c           --- Receiving to the upper guard cells from the right
            i1global = max(me_iy + me_ny + iyus,pe_iylgood)
            i2global = min(me_iy + me_ny + iyue,pe_iyugood)
            if (py > fsdecomp%iyproc+1) then
              i1global = fsdecomp%iy(py-1) + fsdecomp%ny(py-1) + iyus
            endif
          else
c           --- Receive nothing from self
            i1global = 0
            i2global = -1
          endif
          starts(1) = i1global - me_iy + dely
          sizes(1) = i2global - i1global + 1

c         --- Create the data types for data that is to be received
          if (ALL(sizes > 0)) then
            call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                    MPI_ORDER_FORTRAN,
     &                                    MPI_DOUBLE_PRECISION,
     &                                    messagedata%yrecvtypes(py),mpierror)
            call MPI_TYPE_COMMIT(messagedata%yrecvtypes(py),mpierror)
            messagedata%yrecvcounts(py) = 1
          else
            messagedata%yrecvtypes(py) = MPI_DOUBLE_PRECISION
            messagedata%yrecvcounts(py) = 0
          endif

        enddo

        do pz=0,fsdecomp%nzprocs-1

          starts(0) = 0
          sizes(0) = me_nx + 2*delx + 1
          starts(1) = 0
          sizes(1) = me_ny + 2*dely + 1

          pe_iz = fsdecomp%iz(pz)
          pe_nz = fsdecomp%nz(pz)

          if (pz > fsdecomp%izproc) then
c           --- Sending to the lower guard cells on the right
            i1global = max(pe_iz - izls,me_izlgood)
            i2global = min(pe_iz - izle,me_izugood)
            if (pz > fsdecomp%izproc+1) then
              i2global = fsdecomp%iz(fsdecomp%izproc+1) - izle
            endif
          else if (pz < fsdecomp%izproc) then
c           --- Sending to the upper guard cells on the left
            i1global = max(pe_iz + pe_nz + izus,me_izlgood)
            i2global = min(pe_iz + pe_nz + izue,me_izugood)
            if (pz < fsdecomp%izproc-1) then
              i1global = fsdecomp%iz(fsdecomp%izproc-1) +
     &                   fsdecomp%nz(fsdecomp%izproc-1) + izus
            endif
          else
c           --- Send nothing to self
            i1global = 0
            i2global = -1
          endif
          starts(2) = i1global - me_iz + delz
          sizes(2) = i2global - i1global + 1

c         --- Create the data types for data that is to be sent
          if (ALL(sizes > 0)) then
            call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                    MPI_ORDER_FORTRAN,
     &                                    MPI_DOUBLE_PRECISION,
     &                                    messagedata%zsendtypes(pz),mpierror)
            call MPI_TYPE_COMMIT(messagedata%zsendtypes(pz),mpierror)
            messagedata%zsendcounts(pz) = 1
          else
            messagedata%zsendtypes(pz) = MPI_DOUBLE_PRECISION
            messagedata%zsendcounts(pz) = 0
          endif

          pe_izlgood = pe_iz + 1 - izle
          pe_izugood = pe_iz + pe_nz - 1 + izus
          if (pe_iz == 0) pe_izlgood = -delz
          if (pe_iz+pe_nz == fsdecomp%nzglobal) pe_izugood = pe_iz + pe_nz + delz

          if (pz < fsdecomp%izproc) then
c           --- Receiving to the lower guard cells from the left
            i1global = max(me_iz - izls,pe_izlgood)
            i2global = min(me_iz - izle,pe_izugood)
            if (pz < fsdecomp%izproc-1) then
              i2global = fsdecomp%iz(pz+1) - izle
            endif
          else if (pz > fsdecomp%izproc) then
c           --- Receiving to the upper guard cells from the right
            i1global = max(me_iz + me_nz + izus,pe_izlgood)
            i2global = min(me_iz + me_nz + izue,pe_izugood)
            if (pz > fsdecomp%izproc+1) then
              i1global = fsdecomp%iz(pz-1) + fsdecomp%nz(pz-1) + izus
            endif
          else
c           --- Receive nothing from self
            i1global = 0
            i2global = -1
          endif
          starts(2) = i1global - me_iz + delz
          sizes(2) = i2global - i1global + 1

c         --- Create the data types for data that is to be received
          if (ALL(sizes > 0)) then
            call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal,sizes,starts,
     &                                    MPI_ORDER_FORTRAN,
     &                                    MPI_DOUBLE_PRECISION,
     &                                    messagedata%zrecvtypes(pz),mpierror)
            call MPI_TYPE_COMMIT(messagedata%zrecvtypes(pz),mpierror)
            messagedata%zrecvcounts(pz) = 1
          else
            messagedata%zrecvtypes(pz) = MPI_DOUBLE_PRECISION
            messagedata%zrecvcounts(pz) = 0
          endif

        enddo

        return
        end subroutine getmessagedata
c       ======================================================================
        subroutine freemessagedata(messagedata)
        use Decompositionmodule
        type(messagedatatype):: messagedata

        integer(ISZ):: px,py,pz
        integer(MPIISZ):: mpierror

c       --- Note that fsdecomp%nxprocs etc are not used since they may
c       --- have changed

c       --- Free all of the types
        do px=0,SIZE(messagedata%xsendcounts)-1
          if (messagedata%xsendcounts(px) > 0) then
            call MPI_TYPE_FREE(messagedata%xsendtypes(px),mpierror)
          endif
          if (messagedata%xrecvcounts(px) > 0) then
            call MPI_TYPE_FREE(messagedata%xrecvtypes(px),mpierror)
          endif
        enddo
        do py=0,SIZE(messagedata%ysendcounts)-1
          if (messagedata%ysendcounts(py) > 0) then
            call MPI_TYPE_FREE(messagedata%ysendtypes(py),mpierror)
          endif
          if (messagedata%yrecvcounts(py) > 0) then
            call MPI_TYPE_FREE(messagedata%yrecvtypes(py),mpierror)
          endif
        enddo
        do pz=0,SIZE(messagedata%zsendcounts)-1
          if (messagedata%zsendcounts(pz) > 0) then
            call MPI_TYPE_FREE(messagedata%zsendtypes(pz),mpierror)
          endif
          if (messagedata%zrecvcounts(pz) > 0) then
            call MPI_TYPE_FREE(messagedata%zrecvtypes(pz),mpierror)
          endif
        enddo

        deallocate(messagedata%xrecvtypes)
        deallocate(messagedata%xsendtypes)
        deallocate(messagedata%xsendcounts)
        deallocate(messagedata%xrecvcounts)
        deallocate(messagedata%xdispls)
        deallocate(messagedata%yrecvtypes)
        deallocate(messagedata%ysendtypes)
        deallocate(messagedata%ysendcounts)
        deallocate(messagedata%yrecvcounts)
        deallocate(messagedata%ydispls)
        deallocate(messagedata%zrecvtypes)
        deallocate(messagedata%zsendtypes)
        deallocate(messagedata%zsendcounts)
        deallocate(messagedata%zrecvcounts)
        deallocate(messagedata%zdispls)

        return
        end subroutine freemessagedata
      end module messagedatatypemodule
c=============================================================================
c=============================================================================
c=============================================================================
      subroutine mgexchange_phi2(nc,nxlocal,nylocal,nzlocal,phi,
     &                           delx,dely,delz,ils,ile,ius,iue,fsdecomp)
      use Subtimersf3d
      use Decompositionmodule
      use messagedatatypemodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal,ils,ile,ius,iue
      integer(ISZ):: delx,dely,delz
      real(kind=8):: phi(0:nc-1,-delx:nxlocal+delx,
     &                          -dely:nylocal+dely,
     &                          -delz:nzlocal+delz)
      type(Decomposition):: fsdecomp

c Exchange the potential on the boundary and/or in the gaurd cells.

      type(messagedatatype):: messagedata
      save messagedata

      include "mpif.h"
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime,s1,s2
      if (lf3dtimesubs) substarttime = wtime()
      s1 = wtime()

c     --- Generate the MPI data types needed for the communication.
      call getmessagedata(nc,nxlocal,nylocal,nzlocal,
     &                    delx,dely,delz,ils,ile,ius,iue,fsdecomp,
     &                    messagedata)

c     --- Now, do all of the message passing at once.
c     --- Note that there should be no problem with having phi being both
c     --- the send and receive buffer since there is no overlap between
c     --- out going and in coming data.
      if (fsdecomp%nxprocs > 1) then
        call MPI_ALLTOALLW(phi,messagedata%xsendcounts,messagedata%xdispls,
     &                     messagedata%xsendtypes,
     &                     phi,messagedata%xrecvcounts,messagedata%xdispls,
     &                     messagedata%xrecvtypes,
     &                     fsdecomp%mpi_comm_x,mpierror)
      endif
      if (fsdecomp%nyprocs > 1) then
        call MPI_ALLTOALLW(phi,messagedata%ysendcounts,messagedata%ydispls,
     &                     messagedata%ysendtypes,
     &                     phi,messagedata%yrecvcounts,messagedata%ydispls,
     &                     messagedata%yrecvtypes,
     &                     fsdecomp%mpi_comm_y,mpierror)
      endif
      if (fsdecomp%nzprocs > 1) then
        call MPI_ALLTOALLW(phi,messagedata%zsendcounts,messagedata%zdispls,
     &                     messagedata%zsendtypes,
     &                     phi,messagedata%zrecvcounts,messagedata%zdispls,
     &                     messagedata%zrecvtypes,
     &                     fsdecomp%mpi_comm_z,mpierror)
      endif

c     --- Note that the MPI types generated are cached and do not need to be
c     --- released.

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phi = timemgexchange_phi + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_phi(nc,nxlocal,nylocal,nzlocal,phi,
     &                          delx,dely,delz,i1,i2,fsdecomp)
      use Subtimersf3d
      use Decompositionmodule
      use messagedatatypemodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal,i1,i2
      integer(ISZ):: delx,dely,delz
      real(kind=8):: phi(0:nc-1,-delx:nxlocal+delx,
     &                          -dely:nylocal+dely,
     &                          -delz:nzlocal+delz)
      type(Decomposition):: fsdecomp

c Exchange the potential on the boundary and/or in the gaurd cells.

      type(messagedatatype):: messagedata
      save messagedata

      include "mpif.h"
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime,s1,s2
      if (lf3dtimesubs) substarttime = wtime()
      s1 = wtime()

c     --- Generate the MPI data types needed for the communication.
      call getmessagedata(nc,nxlocal,nylocal,nzlocal,
     &                    delx,dely,delz,i1,i2,-i2,-i1,fsdecomp,
     &                    messagedata)

c     --- Now, do all of the message passing at once.
c     --- Note that there should be no problem with having phi being both
c     --- the send and receive buffer since there is no overlap between
c     --- out going and in coming data.
      if (fsdecomp%nxprocs > 1) then
        call MPI_ALLTOALLW(phi,messagedata%xsendcounts,messagedata%xdispls,
     &                     messagedata%xsendtypes,
     &                     phi,messagedata%xrecvcounts,messagedata%xdispls,
     &                     messagedata%xrecvtypes,
     &                     fsdecomp%mpi_comm_x,mpierror)
      endif
      if (fsdecomp%nyprocs > 1) then
        call MPI_ALLTOALLW(phi,messagedata%ysendcounts,messagedata%ydispls,
     &                     messagedata%ysendtypes,
     &                     phi,messagedata%yrecvcounts,messagedata%ydispls,
     &                     messagedata%yrecvtypes,
     &                     fsdecomp%mpi_comm_y,mpierror)
      endif
      if (fsdecomp%nzprocs > 1) then
        call MPI_ALLTOALLW(phi,messagedata%zsendcounts,messagedata%zdispls,
     &                     messagedata%zsendtypes,
     &                     phi,messagedata%zrecvcounts,messagedata%zdispls,
     &                     messagedata%zrecvtypes,
     &                     fsdecomp%mpi_comm_z,mpierror)
      endif

c     --- Note that the MPI types generated are cached and do not need to be
c     --- released.

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phi = timemgexchange_phi + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      subroutine mgexchange_phi_periodic(nc,nxlocal,nylocal,nzlocal,phi,
     &                                   delx,dely,delz,i1,i2,
     &                                   localbounds,fsdecomp)
      use Subtimersf3d
      use Decompositionmodule
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal,nx,ny,nz
      integer(ISZ):: delx,dely,delz,i1,i2
      real(kind=8):: phi(0:nc-1,-delx:nxlocal+delx,
     &                          -dely:nylocal+dely,
     &                          -delz:nzlocal+delz)
      integer(ISZ):: localbounds(0:5)
      type(Decomposition):: fsdecomp
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if ((localbounds(0) == 2 .or. localbounds(1) == 2) .and.
     &    fsdecomp%nxprocs > 0) then
        call mgexchange_phiperiodic_work(0,nc,nxlocal,nylocal,nzlocal,phi,
     &                                   i1,i2,delx,dely,delz,
     &                                   fsdecomp%mpi_comm_x,
     &                                   fsdecomp%ixproc,fsdecomp%nxprocs,
     &                                   fsdecomp%nxglobal,
     &                                   fsdecomp%ix,fsdecomp%nx)
      endif
      if ((localbounds(2) == 2 .or. localbounds(3) == 2) .and.
     &    fsdecomp%nyprocs > 0) then
        call mgexchange_phiperiodic_work(1,nc,nxlocal,nylocal,nzlocal,phi,
     &                                   i1,i2,delx,dely,delz,
     &                                   fsdecomp%mpi_comm_y,
     &                                   fsdecomp%iyproc,fsdecomp%nyprocs,
     &                                   fsdecomp%nyglobal,
     &                                   fsdecomp%iy,fsdecomp%ny)
      endif
      if ((localbounds(4) == 2 .or. localbounds(5) == 2) .and.
     &    fsdecomp%nzprocs > 0) then
        call mgexchange_phiperiodic_work(2,nc,nxlocal,nylocal,nzlocal,phi,
     &                                   i1,i2,delx,dely,delz,
     &                                   fsdecomp%mpi_comm_z,
     &                                   fsdecomp%izproc,fsdecomp%nzprocs,
     &                                   fsdecomp%nzglobal,
     &                                   fsdecomp%iz,fsdecomp%nz)
      endif

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phi_periodic = timemgexchange_phi_periodic + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_phiperiodic_work(axis,nc,nxlocal,nylocal,nzlocal,phi,
     &                                       i1,i2,delx,dely,delz,
     &                                       mpi_comm,iproc,nprocs,nz,
     &                                       pe_iz,pe_nz)
      use Subtimersf3d
      integer(ISZ):: axis
      integer(ISZ):: nc,nxlocal,nylocal,nzlocal,i1,i2,delx,dely,delz
      real(kind=8):: phi(0:nc-1,-delx:nxlocal+delx,
     &                          -dely:nylocal+dely,
     &                          -delz:nzlocal+delz)
      integer(ISZ):: mpi_comm,iproc,nprocs,nz
      integer(ISZ):: pe_iz(0:nprocs-1),pe_nz(0:nprocs-1)

c This routine sends out and receives boundary data for the MG field solver.
c This only sends plane iz=0 when there are periodic boundary conditions.
c Note that this operation needs to be separate from mgexchange_phi since
c there are cases when the same two processors exchange data on the left
c and right, and this cannot be done (easily) with the alltoall.

      real(kind=8),pointer:: buffer(:,:,:,:)
      integer(MPIISZ):: ip
      integer(MPIISZ):: ix1,ix2,iy1,iy2,iz1,iz2

      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE)
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: messid = 10
      integer(MPIISZ):: nn

c     --- By default, select for full dimension of each axis. The down
c     --- select is done below.
      ix1 = -delx
      ix2 = nxlocal + delx
      iy1 = -dely
      iy2 = nylocal + dely
      iz1 = -delz
      iz2 = nzlocal + delz

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.
      if (iproc == 0) then
c       --- The loop only covers processors to the right since any processors
c       --- to the left would not need phi on the left boundary. In some cases,
c       --- this processor my be sending other than boundary data to the left.
c       --- This avoids the send, which would cause this process to lock up
c       --- on the waitall since that message would never be received.
        do ip=1,nprocs-1
          if (pe_iz(ip) + pe_nz(ip) == nz) then
            if (axis == 0) then
              ix1 = -i2
              ix2 = -i1
            elseif (axis == 1) then
              iy1 = -i2
              iy2 = -i1
            elseif (axis == 2) then
              iz1 = -i2
              iz2 = -i1
            endif
            allocate(buffer(0:nc-1,ix1:ix2,iy1:iy2,iz1:iz2))
            buffer = phi(:,ix1:ix2,iy1:iy2,iz1:iz2)
            nn = nc*(ix2-ix1+1)*(iy2-iy1+1)*(iz2-iz1+1)
            call MPI_SEND(buffer,nn,
     &                    MPI_DOUBLE_PRECISION,ip,messid,
     &                    mpi_comm,mpierror)
            deallocate(buffer)
          endif
        enddo
      endif

      if (iproc > 0 .and. pe_iz(iproc) + pe_nz(iproc) == nz) then
        if (axis == 0) then
          ix1 = nxlocal - i2
          ix2 = nxlocal - i1
        elseif (axis == 1) then
          iy1 = nylocal - i2
          iy2 = nylocal - i1
        elseif (axis == 2) then
          iz1 = nzlocal - i2
          iz2 = nzlocal - i1
        endif
        allocate(buffer(0:nc-1,ix1:ix2,iy1:iy2,iz1:iz2))
        nn = nc*(ix2-ix1+1)*(iy2-iy1+1)*(iz2-iz1+1)
        call MPI_RECV(buffer,nn,
     &                MPI_DOUBLE_PRECISION,0,messid,
     &                mpi_comm,mpistatus,mpierror)
        phi(:,ix1:ix2,iy1:iy2,iz1:iz2) = buffer
        deallocate(buffer)
      endif

      if (iproc == nprocs - 1 .and. i1 < 0) then
c       --- See comments on loop above.
        do ip=0,nprocs - 2
          if (pe_iz(ip) == 0) then
            if (axis == 0) then
              ix1 = nxlocal + i1
              ix2 = nxlocal + i1
            elseif (axis == 1) then
              iy1 = nylocal + i1
              iy2 = nylocal + i1
            elseif (axis == 2) then
              iz1 = nzlocal + i1
              iz2 = nzlocal + i1
            endif
            allocate(buffer(0:nc-1,ix1:ix2,iy1:iy2,iz1:iz2))
            buffer = phi(:,ix1:ix2,iy1:iy2,iz1:iz2)
            nn = nc*(ix2-ix1+1)*(iy2-iy1+1)*(iz2-iz1+1)
            call MPI_SEND(buffer,nn,
     &                    MPI_DOUBLE_PRECISION,ip,messid,
     &                    mpi_comm,mpierror)
            deallocate(buffer)
          endif
        enddo
      endif

      if (iproc < nprocs - 1 .and. pe_iz(iproc) == 0 .and. i1 < 0) then
        if (axis == 0) then
          ix1 = i1
          ix2 = i1
        elseif (axis == 1) then
          iy1 = i1
          iy2 = i1
        elseif (axis == 2) then
          iz1 = i1
          iz2 = i1
        endif
        nn = nc*(ix2-ix1+1)*(iy2-iy1+1)*(iz2-iz1+1)
        allocate(buffer(0:nc-1,ix1:ix2,iy1:iy2,iz1:iz2))
        call MPI_RECV(buffer,nn,
     &                MPI_DOUBLE_PRECISION,nprocs-1,messid,
     &                mpi_comm,mpistatus,mpierror)
        phi(:,ix1:ix2,iy1:iy2,iz1:iz2) = buffer
        deallocate(buffer)
      endif

      return
      end
c=============================================================================
      subroutine MY_MPI_ALLTOALLW(ph1,sendcounts,senddispls,sendtypes,
     &                            ph2,recvcounts,recvdispls,recvtypes,
     &                            mpi_comm,mpierror)
      real(kind=8):: ph1(*),ph2(*)
      integer(4),dimension(:):: sendcounts(0:*),senddispls(0:*),sendtypes(0:*)
      integer(4),dimension(:):: recvcounts(0:*),recvdispls(0:*),recvtypes(0:*)
      integer(4):: mpi_comm,mpierror

      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE)

      integer(4):: nprocs,iproc,ip

c Temporary replacement of MPI_ALLTOALLW until the memory leak is fixed.
c This ignores the displacements and assumes that ph1 is the same as ph2,
c (or that the processor sends nothing to itself).

      call MPI_COMM_SIZE(mpi_comm,nprocs,mpierror)
      call MPI_COMM_RANK(mpi_comm,iproc,mpierror)

      do ip=0,nprocs-1
        if (iproc < ip) then
          if (sendcounts(ip) > 0) then
            call MPI_SEND(ph1,sendcounts(ip),sendtypes(ip),ip,7777,
     &                    mpi_comm,mpierror)
          endif
          if (recvcounts(ip) > 0) then
            call MPI_RECV(ph2,recvcounts(ip),recvtypes(ip),ip,7777,
     &                    mpi_comm,mpistatus,mpierror)
          endif
        else if (iproc > ip) then
          if (recvcounts(ip) > 0) then
            call MPI_RECV(ph2,recvcounts(ip),recvtypes(ip),ip,7777,
     &                    mpi_comm,mpistatus,mpierror)
          endif
          if (sendcounts(ip) > 0) then
            call MPI_SEND(ph1,sendcounts(ip),sendtypes(ip),ip,7777,
     &                    mpi_comm,mpierror)
          endif
        else
          if (sendcounts(ip) > 0) then
            call MPI_SENDRECV(ph1,sendcounts(ip),sendtypes(ip),ip,7777,
     &                        ph2,recvcounts(ip),recvtypes(ip),ip,7777,
     &                        mpi_comm,mpistatus,mpierror)
          endif
        endif
      enddo

      return
      end
c=============================================================================
