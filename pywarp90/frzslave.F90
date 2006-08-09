!     Last change:  JLV  11 Jun 2002    4:35 pm
#ifdef MPIPARALLEL
#include "top.h"
module mpirz

INCLUDE 'mpif.h'
INTEGER(MPIISZ) :: ierr
INTEGER(MPIISZ) :: mpi_status(mpi_status_size)

contains

  SUBROUTINE mpi_send_real_array(r, tid, tag)
    IMPLICIT NONE
    INTEGER(ISZ), INTENT(IN) :: TID, tag
    REAL(8), DIMENSION(:) :: r
!    WRITE(0,*) 'send to ',tid
      call mpi_send(r,int(SIZE(r),MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,ierr)
!    WRITE(0,*) 'sent to ',tid
  END SUBROUTINE mpi_send_real_array

  function mpi_recv_real_array(isize,tid,tag)
    implicit none
    INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
    REAL(8), DIMENSION(isize) :: mpi_recv_real_array
!    WRITE(0,*) 'recv from ',tid
      call mpi_recv(mpi_recv_real_array,int(isize,MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world,mpi_status,ierr)
!    WRITE(0,*) 'recvd from ',tid
    return
  end function mpi_recv_real_array

  SUBROUTINE mpi_isend_real_array(r, tid, tag, rqu)
    IMPLICIT NONE
    INTEGER(ISZ), INTENT(IN) :: TID, tag
    INTEGER(MPIISZ), INTENT(IN OUT) :: rqu
    REAL(8), DIMENSION(:) :: r

!    WRITE(0,*) 'send to ',tid
      call mpi_isend(r,int(SIZE(r),MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,rqu,ierr)
!    WRITE(0,*) 'sent to ',tid
  END SUBROUTINE mpi_isend_real_array

  function mpi_irecv_real_array(isize,tid,tag,rqu)
    implicit none
    INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
    INTEGER(MPIISZ), INTENT(IN OUT) :: rqu
    REAL(8), DIMENSION(isize) :: mpi_irecv_real_array
!    WRITE(0,*) 'recv from ',tid
      call mpi_irecv(mpi_irecv_real_array,int(isize,MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world,rqu,ierr)
!    WRITE(0,*) 'recvd from ',tid
    return
  end function mpi_irecv_real_array

  function mpi_global_compute_real(DATA,op)
    implicit none
    REAL(8) :: DATA, mpi_global_compute_real
    INTEGER(MPIISZ) :: op
!    WRITE(0,*) 'enter mpi_allreduce '
      call mpi_allreduce(data,mpi_global_compute_real,int(1,MPIISZ),mpi_double_precision,int(op,MPIISZ),mpi_comm_world,ierr)
!    WRITE(0,*) 'exit mpi_allreduce '
    return
  end function mpi_global_compute_real

end module mpirz
#endif
