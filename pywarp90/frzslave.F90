!     Last change:  JLV   7 Sep 2001    9:59 am
#ifdef PARALLEL
#include "top.h"
module mpirz

INCLUDE 'mpif.h'
INTEGER(ISZ) :: ierr
INTEGER(ISZ) :: mpi_status(mpi_status_size)

contains

  SUBROUTINE mpi_send_real_array(r, tid, tag)
    IMPLICIT NONE
    INTEGER(ISZ), INTENT(IN) :: TID, tag
    REAL(8), DIMENSION(:) :: r
      call mpi_send(r,SIZE(r),mpi_real8,tid,tag,MPI_COMM_WORLD,ierr)
  END SUBROUTINE mpi_send_real_array

  function mpi_recv_real_array(isize,tid,tag)
    implicit none
    INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
    REAL(8), DIMENSION(isize) :: mpi_recv_real_array
      call mpi_recv(mpi_recv_real_array,isize,mpi_real8,tid,tag,mpi_comm_world,mpi_status,ierr)
    return
  end function mpi_recv_real_array

  function mpi_global_compute_real(DATA,op)
    implicit none
    REAL(8) :: DATA, mpi_global_compute_real
    INTEGER(ISZ) :: op
      call mpi_allreduce(data,mpi_global_compute_real,1,mpi_real8,op,mpi_comm_world,ierr)
    return
  end function mpi_global_compute_real

end module mpirz
#endif
