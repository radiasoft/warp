#include "top.h"
c=============================================================================
c F3DSLAVE.M, version $Revision: 1.38 $, $Date: 2008/07/31 21:01:04 $
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
        call MPI_ISEND(phi(-1,-1,left_send),int((nx+3)*(ny+3),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,left_pe,int(nx+ny,MPIISZ),
     &                 MPI_COMM_WORLD,mpirequest(1),mpierror)
        call MPI_IRECV(phi(-1,-1,left_recv),int((nx+3)*(ny+3),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,left_pe,int(nx+ny,MPIISZ),
     &                 MPI_COMM_WORLD,mpirequest(2),mpierror)
      endif

      if (my_index < nslaves-1 .or. boundnz == 2) then
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
      integer(MPIISZ):: sendtypes(0:nslaves-1),recvtypes(0:nslaves-1)
      integer(MPIISZ):: sdispls(0:nslaves-1),rdispls(0:nslaves-1)
      integer(MPIISZ):: sendcounts(0:nslaves-1),recvcounts(0:nslaves-1)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Calculate what the dimensions of the grid will be after the
c     --- transpose, (0:nx_tran,0:ny_tran).
c     --- The work is divided up first along the y direction. Since
c     --- all dimensions and number of processors are powers of 2, ny/nslaves
c     --- will always also be a power of 2 (or zero).
c     --- If ny > nslaves, each processor is given ny/nslaves lines along x.
c     --- If nx*ny > nslaves, then the x lines are divided up as well.
c     --- Note that then the transposed data is stored as 2 'y' lines, each
c     --- nx*ny/nslaves/2 long in x. This makes it easier to use the 2 at a
c     --- time operation of vpftz.
      if (ny > nslaves) then
        nx_tran = nx
        ny_tran = ny/nslaves
      else if (nx*ny > nslaves) then
        nx_tran = nx*ny/nslaves/2
        ny_tran = 2
      else
        call kaboom("warptranspose: the number of processors must be less than nx*ny")
        return
      endif

c     --- Create the work space that phi is transposed into.
      phi_trnspsnx = nx_tran
      phi_trnspsny = ny_tran
      phi_trnspsnz = nzlocal*nslaves
      call gchange("Transpose_work_space",0)

c     --- Create an MPI type to refer to the distributed pieces of data that
c     --- are sent to each of the processors.
      dims = 3
      nlocal = (/nx+3,ny+3,nzlocal/)
      nlocal_trnsps = (/nx_tran+3,ny_tran+3,nzlocal*nslaves/)
      sdispls = 0
      rdispls = 0
      sendcounts = 1
      recvcounts = 1

      do ip = 0,nslaves-1

        if (ny < nslaves) then
          starts(0) = 1 + (2*nx_tran)*mod(ip,nslaves/ny)
          starts(1) = 1 + int(ip*ny/nslaves)
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
     &                                sendtypes(ip),mpierror)
        call MPI_TYPE_COMMIT(sendtypes(ip),mpierror)

        starts(0) = 1
        starts(1) = 1
        starts(2) = ip*nzlocal
        sizes = (/nx_tran,ny_tran,nzlocal/)
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
      do ip = 0,nslaves-1
        call MPI_TYPE_FREE(sendtypes(ip),mpierror)
        call MPI_TYPE_FREE(recvtypes(ip),mpierror)
      enddo

c     --- Now the array phi_trnsps contains 3-D data with dimensions
c     --- (0:nx_tran,0:ny_tran,nzlocal*nslaves).

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
      integer(MPIISZ):: sendtypes(0:nslaves-1),recvtypes(0:nslaves-1)
      integer(MPIISZ):: sdispls(0:nslaves-1),rdispls(0:nslaves-1)
      integer(MPIISZ):: sendcounts(0:nslaves-1),recvcounts(0:nslaves-1)
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- Create an MPI type to refer to the distributed pieces of data that
c     --- are sent to each of the processors.
      dims = 3
      nlocal = (/nx+3,ny+3,nzlocal/)
      nlocal_trnsps = (/nx_tran+3,ny_tran+3,nzlocal*nslaves/)
      sdispls = 0
      rdispls = 0
      sendcounts = 1
      recvcounts = 1

      do ip = 0,nslaves-1

        starts(0) = 1
        starts(1) = 1
        starts(2) = ip*nzlocal
        sizes = (/nx_tran,ny_tran,nzlocal/)
        call MPI_TYPE_CREATE_SUBARRAY(dims,nlocal_trnsps,sizes,starts,
     &                                MPI_ORDER_FORTRAN,
     &                                MPI_DOUBLE_PRECISION,
     &                                sendtypes(ip),mpierror)
        call MPI_TYPE_COMMIT(sendtypes(ip),mpierror)

        if (ny < nslaves) then
          starts(0) = 1 + (2*nx_tran)*mod(ip,nslaves/ny)
          starts(1) = 1 + int(ip*ny/nslaves)
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
      do ip = 0,nslaves-1
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
        nzlan = 2*nslaves
        nzlocallan = nz
        if (ny >= nslaves) then
          nxtranlan = nx
          nytranlan = ny/nslaves
        else
          nxtranlan = nx*ny/nslaves
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
      call warptranspose(nx,ny,2,dtranlan,nx_tran,ny_tran)

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
      if (nproc < nslaves) then
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
      subroutine mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                      nz,nzcoarse,mgscale)
      use Subtimersf3d
      use Multigrid3d,Only:mgscaleserial
      integer(ISZ):: nslaves,nz,nzcoarse
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      real(kind=8):: mgscale

c Divides nz for each processor by two. For odd numbers, the starting end
c is rounded down, and the upper end up. For each processor, if nz ends up
c being equal to 1, then change it to 2. For the last processor, if nz is 1,
c also, subtract 1 from izfsslave.

      integer(ISZ):: is,nzmin = 2
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      if (mgscale < mgscaleserial) then

        do is=0,nslaves-1

c         --- First divide nz by 2. izfsslave is rounded down, and nzfsslave
c         --- is rounded up.
          izfsslavec(is) = int(floor(1.*izfsslave(is)*nzcoarse/nz))
          nzfsslavec(is) = int(ceiling(1.*(izfsslave(is)+nzfsslave(is))*
     &                                  nzcoarse/nz - izfsslavec(is)))

c         --- Then check that no nz are less than nzmin.
c         --- Note that if nzfsslave is already < nzmin, leave it as is.
c         --- When nzfsslave < nzmin, in some cases setting nzfsslavec to nzmin
c         --- will cause the coarsened grid to extend too far beyond the finer
c         --- grid's gaurd cell, so it would need data from a neighboring cell,
c         --- where it would normally be able to use data from the gaurd cells.
c         --- That is a rare occurance so it is simpler to avoid the problem
c         --- then to write the code to exchange the data for that case.
          if (nzfsslavec(is) < nzmin .and. nzfsslavec(is) < nzfsslave(is)) then
            nzfsslavec(is) = nzfsslavec(is) + 1
            if (izfsslavec(is)+nzfsslavec(is) > nzcoarse) then
              izfsslavec(is) = izfsslavec(is) - 1
              nzfsslavec(is) = nzcoarse - izfsslavec(is)
            endif
          endif

        enddo

      else

c       --- When the volume of the coarsened grid gets small compared to
c       --- the finest level, then expand the domain of each processor to
c       --- cover the full system, to avoid passing around lots of small
c       --- messages.
        izfsslavec = 0
        nzfsslavec = nzcoarse

      endif

!$OMP MASTER
      if (lf3dtimesubs) timemgdividenz = timemgdividenz + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mggetexchangepes(nslaves,izfsslave,nzfsslave,my_index,
     &                            bounds,nz,lzoffset,rzoffset,
     &                            whosendingleft,izsendingleft,
     &                            whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nslaves,my_index
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: bounds(0:5),nz
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
        if (izfsslave(idx) == 0 .and. izfsslave(idx)+nzfsslave(idx) == nz
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

c         --- Find all processors whose data extent overlaps the right
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
        if (izfsslave(idx) == 0 .and. izfsslave(idx)+nzfsslave(idx) == nz
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
          do while (js >= 0 .and. izfsslave(js)+nzfsslave(js) == nz .and.
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
      subroutine mgexchange_phi(nx,ny,nzlocal,nz,phi,bounds,zsend,delx,dely,delz,
     &                          my_index,nslaves,izfsslave,nzfsslave,
     &                          whosendingleft,izsendingleft,
     &                          whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nzlocal,nz,zsend,delx,dely,delz
      real(kind=8):: phi(-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c This routine sends out and receives boundary data for the MG field solver.

      integer(MPIISZ):: leftdelta,rightdelta
      integer(MPIISZ):: lr,rr,ls,rs
      integer(MPIISZ):: ip,w

      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,2*nslaves+2)
      integer(MPIISZ):: mpirequest(2*nslaves+2),mpierror
      integer(MPIISZ):: messid
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()
      messid = nz

c     --- Location where incoming phi plane is to be loaded.  Normally, it
c     --- is put at the edge of the grid, either at iz=0 or iz=nzlocal.  For
c     --- periodic boundaries, it is put just outside the grid, at iz=-1
c     --- or iz=nzlocal+1.  For the above, zsend=0.  At the end of the field solve,
c     --- phi at the outside of the grid is passed, so phi is put at iz=-1 and
c     --- iz=nzlocal+1 so zsend needs to be -1.
      lr = zsend
      rr = nzlocal - zsend
      if (bounds(4) == 2) lr = -1
      if (bounds(5) == 2) rr = nzlocal + 1

c     --- The actual value of zsend used depends on the boundary conditions.
c     --- For periodic boundary conditions, processors which are on the ends
c     --- of the mess, always send the same data, independent of what zsend is.
      leftdelta = zsend
      rightdelta = zsend
      if (izfsslave(my_index) == 0 .and. bounds(4) == 2) leftdelta = 0
      if (izfsslave(my_index)+nzfsslave(my_index) == nz .and. bounds(5) == 2)
     &  rightdelta = 0

c     --- Set count to zero
      w = 0

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.
      if (mod(my_index,2) == 0) then
        do ip=0,nslaves-1
          if (my_index == whosendingleft(ip)) then
            w = w + 1
            ls = izsendingleft(ip) - leftdelta
            call MPI_SEND(phi(-delx,-dely,ls),int((nx+1+2*delx)*(ny+1+2*dely),MPIISZ),
     &                     MPI_DOUBLE_PRECISION,ip,messid,
     &                     MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif
      if (whosendingleft(my_index) /= -1) then
        w = w + 1
        call MPI_RECV(phi(-delx,-dely,rr),int((nx+1+2*delx)*(ny+1+2*dely),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,int(whosendingleft(my_index),MPIISZ),
     &                 messid,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif
      if (mod(my_index,2) == 1) then
        do ip=0,nslaves-1
          if (my_index == whosendingleft(ip)) then
            w = w + 1
            ls = izsendingleft(ip) - leftdelta
            call MPI_SEND(phi(-delx,-dely,ls),int((nx+1+2*delx)*(ny+1+2*dely),MPIISZ),
     &                     MPI_DOUBLE_PRECISION,ip,messid,
     &                     MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif

      if (mod(my_index,2) == 0) then
        do ip=0,nslaves-1
          if (my_index == whosendingright(ip)) then
            w = w + 1
            rs = izsendingright(ip) + rightdelta
            call MPI_SEND(phi(-delx,-dely,rs),int((nx+1+2*delx)*(ny+1+2*dely),MPIISZ),
     &                     MPI_DOUBLE_PRECISION,ip,messid,
     &                     MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif
      if (whosendingright(my_index) /= -1) then
        w = w + 1
        call MPI_RECV(phi(-delx,-dely,lr),int((nx+1+2*delx)*(ny+1+2*dely),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,int(whosendingright(my_index),MPIISZ),
     &                 messid,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif
      if (mod(my_index,2) == 1) then
        do ip=0,nslaves-1
          if (my_index == whosendingright(ip)) then
            w = w + 1
            rs = izsendingright(ip) + rightdelta
            call MPI_SEND(phi(-delx,-dely,rs),int((nx+1+2*delx)*(ny+1+2*dely),MPIISZ),
     &                     MPI_DOUBLE_PRECISION,ip,messid,
     &                     MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif

c     --- Now wait for everything to finish
c     if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phi = timemgexchange_phi + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_phiperiodic(nx,ny,nzlocal,nz,phi,bounds,zsend,
     &                                  delx,dely,delz,my_index,nslaves,izfsslave,
     &                                  whosendingleft,whosendingright)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nzlocal,nz,zsend,delx,dely,delz
      real(kind=8):: phi(-delx:nx+delx,-dely:ny+dely,-delz:nzlocal+delz)
      integer(ISZ):: bounds(0:5)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1),whosendingright(0:nslaves-1)

c This routine sends out and receives boundary data for the MG field solver.
c This only sends plane iz=0 when there are periodic boundary conditions.

      integer(MPIISZ):: ip,w

      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nslaves)
      integer(MPIISZ):: mpirequest(nslaves),mpierror
      integer(MPIISZ):: messid = 10
      integer(MPIISZ):: nn
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
            nn = (nx+1+2*delx)*(ny+1+2*dely)
            call MPI_SEND(phi(-delx,-dely,zsend),nn,
     &                    MPI_DOUBLE_PRECISION,ip,messid,
     &                    MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif

      if (izfsslave(my_index)+nzlocal == nz) then
        if (whosendingleft(my_index) /= -1) then
          w = w + 1
          nn = (nx+1+2*delx)*(ny+1+2*dely)
          ip = whosendingleft(my_index)
          call MPI_RECV(phi(-delx,-dely,nzlocal+zsend),nn,
     &                  MPI_DOUBLE_PRECISION,ip,messid,
     &                  MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      endif

      if (izfsslave(my_index)+nzlocal == nz .and. zsend /= 0) then
c       --- See comments on loop above.
        do ip=0,my_index-1
          if (my_index == whosendingright(ip)) then
            w = w + 1
            nn = (nx+1+2*delx)*(ny+1+2*dely)
            call MPI_SEND(phi(-delx,-dely,nzlocal-zsend),nn,
     &                    MPI_DOUBLE_PRECISION,ip,messid,
     &                    MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        enddo
      endif

      if (izfsslave(my_index) == 0 .and. zsend /= 0) then
        if (whosendingright(my_index) /= -1) then
          w = w + 1
          nn = (nx+1+2*delx)*(ny+1+2*dely)
          ip = whosendingright(my_index)
          call MPI_RECV(phi(-delx,-dely,-zsend),nn,
     &                  MPI_DOUBLE_PRECISION,ip,messid,
     &                  MPI_COMM_WORLD,mpirequest(w),mpierror)
        endif
      endif

c     --- Now wait for everything to finish
c     if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_phiperiodic = timemgexchange_phiperiodic + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine mgexchange_rho(nx,ny,nzlocal,nz,rho,
     &                          my_index,nslaves,izfsslave,nzfsslave,
     &                          lzoffset,rzoffset,
     &                          whosendingleft,izsendingleft,
     &                          whosendingright,izsendingright)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nzlocal,nz
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)
      integer(ISZ):: lzoffset(0:nslaves-1),rzoffset(0:nslaves-1)
      integer(ISZ):: whosendingleft(0:nslaves-1), izsendingleft(0:nslaves-1)
      integer(ISZ):: whosendingright(0:nslaves-1),izsendingright(0:nslaves-1)

c This routine sends out and receives boundary data for the MG field solver.

      integer(MPIISZ):: ls,rs
      integer(MPIISZ):: ip,w

      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,2*nslaves+2)
      integer(MPIISZ):: mpirequest(2*nslaves+2),mpierror
      integer(MPIISZ):: messid
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()
      messid = nz

c     --- Set count to zero
      w = 0

c     --- Sends and receives can all be done synchronously since none
c     --- of the data overlaps.

c     --- The last processor is skipped since it doesn't need to have rho
c     --- on the right end sent to it.
      do ip=0,nslaves-2

c       --- Also skip processors that extend to the right edge of this and the
c       --- next finer level.
        if (izfsslave(ip)+nzfsslave(ip)==nz .and. rzoffset(ip)==0.) cycle

        if (my_index == whosendingleft(ip)) then
          w = w + 1
          ls = izsendingleft(ip)

c         --- When my_index < ip, it means that the data being sent is periodic.
c         --- In this case, the izsendingleft(ip) is set to send the guard cell,
c         --- so override it.
          if (my_index < ip) ls = 0

          call MPI_ISEND(rho(0,0,ls),int((nx+1)*(ny+1),MPIISZ),
     &                   MPI_DOUBLE_PRECISION,ip,messid,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)

        endif
      enddo

c     --- Any processor that extends to the left edge of this and the next
c     --- finer level don't need a rho sent to it.
      if (whosendingright(my_index) /= -1 .and. 
     &    .not. (izfsslave(my_index)==0 .and. lzoffset(my_index)==0)) then
        w = w + 1
        call MPI_IRECV(rho(0,0,0),int((nx+1)*(ny+1),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,int(whosendingright(my_index),MPIISZ),
     &                 messid,
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
          if (my_index > ip) rs = nzlocal

          call MPI_ISEND(rho(0,0,rs),int((nx+1)*(ny+1),MPIISZ),
     &                   MPI_DOUBLE_PRECISION,ip,messid,
     &                   MPI_COMM_WORLD,mpirequest(w),mpierror)

        endif
      enddo

c     --- Any processor that extends to the right edge of this and the next
c     --- finer level don't need a rho sent to it.
      if (whosendingleft(my_index) /= -1 .and.
     &    .not. (izfsslave(my_index)+nzfsslave(my_index)==nz .and.
     &           rzoffset(my_index)==0.)) then
        w = w + 1
        call MPI_IRECV(rho(0,0,nzlocal),int((nx+1)*(ny+1),MPIISZ),
     &                 MPI_DOUBLE_PRECISION,int(whosendingleft(my_index),MPIISZ),
     &                 messid,
     &                 MPI_COMM_WORLD,mpirequest(w),mpierror)
      endif

c     --- Now wait for everything to finish
      if (w > 0) call MPI_WAITALL(w,mpirequest,mpistatus,mpierror)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchange_rho = timemgexchange_rho + wtime() - substarttime
!$OMP END MASTER

      return
      end
c===========================================================================
      subroutine mgexchangeallrhocoarse(nc,nx,ny,nz,source,nzfine,
     &                                  my_index,nslaves,izfsslave,nzfsslave)
      use Subtimersf3d
      integer(ISZ):: nc,nx,ny,nz,nzfine
      real(kind=8):: source(0:nc-1,0:nx,0:ny,0:nz)
      integer(ISZ):: my_index,nslaves
      integer(ISZ):: izfsslave(0:nslaves-1),nzfsslave(0:nslaves-1)

      real(kind=8),allocatable:: source2(:,:,:,:)
      integer(MPIISZ):: i
      integer(ISZ):: izfsslavec(0:nslaves-1),nzfsslavec(0:nslaves-1)
      integer(MPIISZ):: ii(0:nslaves-1),nn(0:nslaves-1)
      integer(MPIISZ):: is(0:nslaves-1),ns(0:nslaves-1)
      include "mpif.h"
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

      call mgdividenz(nslaves,izfsslave,nzfsslave,izfsslavec,nzfsslavec,
     &                nzfine,nz,1.)

      do i=0,nslaves-2
        ii(i) = izfsslavec(i)
        nn(i) = max(0,izfsslavec(i+1) - izfsslavec(i))
      enddo
      ii(nslaves-1) = izfsslavec(nslaves-1)
      nn(nslaves-1) = nzfsslavec(nslaves-1) + 1

      ii = ii*nc*(nx+1)*(ny+1)
      nn = nn*nc*(nx+1)*(ny+1)

      is = ii(my_index)
      ns = nn(my_index)

c     allocate(source2(0:nc-1,0:nx,0:ny,0:nz))
c     source2 = source
      call MPI_ALLTOALLV(source,ns,is,MPI_DOUBLE_PRECISION,
     &                   source,nn,ii,MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD,mpierror)
c     deallocate(source2)

!$OMP MASTER
      if (lf3dtimesubs) timemgexchangeallrhocoarse = timemgexchangeallrhocoarse + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      subroutine printarray3d(nx,ny,nz,iz0,nz0,izs,a,my_index,nslaves,text,mglevel)
      use Subtimersf3d
      integer(ISZ):: nx,ny,nz,iz0,nz0,izs,my_index,nslaves,mglevel
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
      if (my_index == nslaves-1) lastz = nz

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

      integer(MPIISZ):: i,nn
      integer(MPIISZ):: izglobal,izmaxp,izmaxfs
      integer(MPIISZ):: izsend(0:nslaves-1),nzsend(0:nslaves-1)
      integer(MPIISZ):: izrecv(0:nslaves-1),nzrecv(0:nslaves-1)
      include "mpif.h"
      integer(MPIISZ):: mpistatus(MPI_STATUS_SIZE,nslaves),mpirequest(nslaves)
      integer(MPIISZ):: mpierror,w
      integer(MPIISZ):: messid = 80
      real(kind=8):: substarttime,wtime
      if (lf3dtimesubs) substarttime = wtime()

c     --- First, calculate what data needs to be sent and received.
c     --- Note that special cases are needed for the first and last processors
c     --- so that the guard cells are filled properly.
c     --- Data to be sent
      do i=0,nslaves-1
        izglobal = max(izfsslave(my_index),izpslave(i))
        izmaxp   = izpslave(i)+nzpslave(i)
        izmaxfs  = izfsslave(my_index)+nzfsslave(my_index)

        izsend(i) = izglobal - izfsslave(my_index)
        nzsend(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        izglobal = max(izpslave(my_index),izfsslave(i))
        izmaxp   = izpslave(my_index)+nzpslave(my_index)
        izmaxfs  = izfsslave(i)+nzfsslave(i)

        izrecv(i) = izglobal - izpslave(my_index)
        nzrecv(i) = min(izmaxp,izmaxfs) - izglobal + 1
      enddo
 
c     --- Send the data out to processors that need it.
      w = 0
      do i=0,nslaves-1
        if (nzsend(i) > 0) then
          if ( i /= my_index) then
            w = w + 1
            nn = int(3*nzsend(i)*(bfield%nx+1)*(bfield%ny+1),MPIISZ)
            call MPI_ISEND(bfield%b(:,:,:,izsend(i)),nn,
     &                     MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          else
            bfieldp%b(:,:,:,izrecv(i):izrecv(i)+nzrecv(i)-1) =
     &        bfield%b(:,:,:,izsend(i):izsend(i)+nzsend(i)-1)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then
          if (nzrecv(i) > 0) then
            nn = int(3*nzrecv(i)*(bfieldp%nx+1)*(bfieldp%ny+1),MPIISZ)
            call MPI_RECV(bfieldp%b(:,:,:,izrecv(i)),nn,
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
      if (lf3dtimesubs) timegetbforparticles3d = timegetbforparticles3d + wtime() - substarttime
!$OMP END MASTER

      return
      end                                                                       
c===========================================================================
      subroutine getaforfields3d(a,nx,ny,nz)
      use Subtimersf3d
      use Parallel
      integer(ISZ):: nx,ny,nz
      real(kind=8):: a(0:2,-1:nx+1,-1:ny+1,-1:nz+1)

c Get the A for the full extent where the fields are. This gets A from
c neighboring processors, and at the very least gets A in the guard planes,
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
      if (lf3dtimesubs) substarttime = wtime()
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
            nn = 3*nzsend(i)*(nx+3)*(ny+3)
            call MPI_ISEND(a(:,:,:,izsend(i)),nn,
     &                     MPI_DOUBLE_PRECISION,
     &                     i,messid,MPI_COMM_WORLD,mpirequest(w),mpierror)
          endif
        endif
      enddo

c     --- Then, gather up the data sent to this processor.
      do i=0,nslaves-1
        if ( i /= my_index) then
          if (nzrecv(i) > 0) then
            nn = 3*nzrecv(i)*(nx+3)*(ny+3)
            call MPI_RECV(a(:,:,:,izrecv(i)),nn,
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
      if (lf3dtimesubs) timegetaforfields3d = timegetaforfields3d + wtime() - substarttime
!$OMP END MASTER

      return
      end                                                                       
c===========================================================================
