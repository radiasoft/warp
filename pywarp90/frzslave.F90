!     Last change:  JLV  11 Jun 2002    4:35 pm
#ifdef MPIPARALLEL
#include "top.h"
module mpirz
use Parallel
INCLUDE 'mpif.h'
INTEGER(MPIISZ) :: ierr, size_packbuffer, pack_pos
INTEGER(MPIISZ) :: mpi_status(mpi_status_size)
 integer*8, allocatable :: packbuffer(:)
 logical, parameter :: l_mpiverbose=.false.

! INTERFACE mpi_send_int
!   MODULE PROCEDURE mpi_send_int_scalar
!   MODULE PROCEDURE mpi_send_int_array
! END interface
! INTERFACE mpi_send_real
!   MODULE PROCEDURE mpi_send_real_scalar
!   MODULE PROCEDURE mpi_send_real_array
! END interface
 INTERFACE mympi_send
   MODULE PROCEDURE mpi_send_int_scalar
   MODULE PROCEDURE mpi_send_int_array
   MODULE PROCEDURE mpi_send_real_scalar
   MODULE PROCEDURE mpi_send_real_array
 END interface
 INTERFACE mympi_isend
   MODULE PROCEDURE mpi_isend_real_array
 END interface
 INTERFACE mpi_recv_int
   MODULE PROCEDURE mpi_recv_int_scalar
   MODULE PROCEDURE mpi_recv_int_array
 END interface
 INTERFACE mpi_recv_real
   MODULE PROCEDURE mpi_recv_real_scalar
   MODULE PROCEDURE mpi_recv_real_array
 END interface
 INTERFACE mpi_irecv_real
   MODULE PROCEDURE mpi_irecv_real_array
 END interface
 INTERFACE mympi_pack
   MODULE PROCEDURE mpi_pack_int_scalar
   MODULE PROCEDURE mpi_pack_int_array
   MODULE PROCEDURE mpi_pack_real_scalar
   MODULE PROCEDURE mpi_pack_real_1darray
   MODULE PROCEDURE mpi_pack_real_2darray
   MODULE PROCEDURE mpi_pack_real_3darray
 END interface
 INTERFACE mpi_unpack_int
   MODULE PROCEDURE mpi_unpack_int_scalar
   MODULE PROCEDURE mpi_unpack_int_array
 END interface
 INTERFACE mpi_unpack_real
   MODULE PROCEDURE mpi_unpack_real_scalar
   MODULE PROCEDURE mpi_unpack_real_array
 END interface

contains

 SUBROUTINE mpi_send_int_scalar(i, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, i, tag
  call mpi_send(i,1,mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,ierr)
 END SUBROUTINE mpi_send_int_scalar

 SUBROUTINE mpi_send_int_array(i, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, i(:), tag
  call mpi_send(i,SIZE(i),mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,ierr)
 END SUBROUTINE mpi_send_int_array

 function mpi_recv_int_scalar(tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   INTEGER(ISZ) :: mpi_recv_int_scalar
     call mpi_recv(mpi_recv_int_scalar,1,mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world,mpi_status_ignore,ierr)
   return
 end function mpi_recv_int_scalar

 function mpi_recv_int_array(isize,tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
   INTEGER(ISZ), DIMENSION(isize) :: mpi_recv_int_array
     call mpi_recv(mpi_recv_int_array,int(isize,MPIISZ),mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, &
                   mpi_status_ignore,ierr)
   return
 end function mpi_recv_int_array

 SUBROUTINE mpi_send_real_scalar(r, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, tag
  REAL(8) :: r
  call mpi_send(r,1,mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,ierr)
 END SUBROUTINE mpi_send_real_scalar

 SUBROUTINE mpi_send_real_array(r, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, tag
  REAL(8), DIMENSION(:) :: r
  if (l_mpiverbose) WRITE(0,*) my_index,'send to ',tid
 call mpi_send(r,int(SIZE(r),MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,ierr)
  if (l_mpiverbose) WRITE(0,*) my_index,'sent to ',tid
 END SUBROUTINE mpi_send_real_array

 function mpi_recv_real_scalar(tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   REAL(8) :: mpi_recv_real_scalar
     call mpi_recv(mpi_recv_real_scalar,1,mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, &
                   mpi_status_ignore,ierr)
   return
 end function mpi_recv_real_scalar

 function mpi_recv_real_array(isize,tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
   REAL(8), DIMENSION(isize) :: mpi_recv_real_array
     if (l_mpiverbose) WRITE(0,*) my_index,'recv from ',tid,isize
     call mpi_recv(mpi_recv_real_array,int(isize,MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, &
                   mpi_status_ignore,ierr)
    if (l_mpiverbose) WRITE(0,*) my_index,'recvd from ',tid
   return
 end function mpi_recv_real_array

  SUBROUTINE mpi_isend_real_array(r, tid, tag, rqu)
    IMPLICIT NONE
    INTEGER(ISZ), INTENT(IN) :: TID, tag
    INTEGER(MPIISZ), INTENT(IN OUT) :: rqu
    REAL(8), DIMENSION(:) :: r

!    WRITE(0,*) 'isend to ',tid,size(r)
      call mpi_isend(r,int(SIZE(r),MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),MPI_COMM_WORLD,rqu,ierr)
!    WRITE(0,*) 'isent to ',tid
  END SUBROUTINE mpi_isend_real_array

  function mpi_irecv_real_array(isize,tid,tag,rqu)
    implicit none
    INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
    INTEGER(MPIISZ), INTENT(IN OUT) :: rqu
    REAL(8), DIMENSION(isize) :: mpi_irecv_real_array
!    WRITE(0,*) 'irecv from ',tid
      call mpi_irecv(mpi_irecv_real_array,int(isize,MPIISZ),mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ), &
                     mpi_comm_world,rqu,ierr)
!    WRITE(0,*) 'irecvd from ',tid
    return
  end function mpi_irecv_real_array

  function mpi_global_compute_real(DATA,op)
    implicit none
    REAL(8) :: DATA, mpi_global_compute_real
    INTEGER(MPIISZ) :: op
      call mpi_allreduce(data,mpi_global_compute_real,int(1,MPIISZ),mpi_double_precision,int(op,MPIISZ),mpi_comm_world,ierr)
    return
  end function mpi_global_compute_real

 subroutine mpi_packbuffer_init(isize)
 implicit none
 INTEGER(ISZ), INTENT(IN) :: isize
   pack_pos=0
   size_packbuffer = 8*isize
   IF(ALLOCATED(packbuffer)) DEALLOCATE(packbuffer)
   ALLOCATE(packbuffer(size_packbuffer))
 return
 END subroutine mpi_packbuffer_init

 subroutine mpi_pack_int_scalar(a)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: a
     call mpi_pack(a, 1, mpi_integer, packbuffer, size_packbuffer, pack_pos, mpi_comm_world, ierr)
   return
 end subroutine mpi_pack_int_scalar

 subroutine mpi_pack_int_array(a)
   implicit none
   INTEGER(ISZ), DIMENSION(:), INTENT(IN) :: a
     call mpi_pack(a, SIZE(a), mpi_integer, packbuffer, size_packbuffer, pack_pos, mpi_comm_world, ierr)
   return
 end subroutine mpi_pack_int_array

 subroutine mpi_pack_real_scalar(a)
   implicit none
   REAL(8), INTENT(IN) :: a
     call mpi_pack(a, 1, mpi_double_precision, packbuffer, size_packbuffer, pack_pos, mpi_comm_world, ierr)
   return
 end subroutine mpi_pack_real_scalar

 subroutine mpi_pack_real_1darray(a)
   implicit none
   REAL(8), DIMENSION(:), INTENT(IN) :: a
     call mpi_pack(a, SIZE(a), mpi_double_precision, packbuffer, size_packbuffer, pack_pos, mpi_comm_world, ierr)
   return
 end subroutine mpi_pack_real_1darray

 subroutine mpi_pack_real_2darray(a)
   implicit none
   REAL(8), DIMENSION(:,:), INTENT(IN) :: a
     call mpi_pack(a(1,1), SIZE(a), mpi_double_precision, packbuffer, size_packbuffer, pack_pos, mpi_comm_world, ierr)
   return
 end subroutine mpi_pack_real_2darray

 subroutine mpi_pack_real_3darray(a)
   implicit none
   REAL(8), DIMENSION(:,:,:), INTENT(IN) :: a
     call mpi_pack(a(1,1,1), SIZE(a), mpi_double_precision, packbuffer, size_packbuffer, pack_pos, mpi_comm_world, ierr)
   return
 end subroutine mpi_pack_real_3darray

 subroutine mpi_send_pack(tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag
     call mpi_send(packbuffer,pack_pos,mpi_packed,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, ierr)
   return
 end subroutine mpi_send_pack

 subroutine mpi_isend_pack(tid,tag,rqu)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag,rqu
     call mpi_isend(packbuffer,pack_pos,mpi_packed,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, rqu, ierr)
   return
 end subroutine mpi_isend_pack

 subroutine mpi_recv_pack(tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag
     call mpi_recv(packbuffer,size_packbuffer,mpi_packed,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, &
                   mpi_status_ignore, ierr)
   return
 end subroutine mpi_recv_pack

 subroutine mpi_irecv_pack(tid,tag,rqu)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag,rqu
     call mpi_irecv(packbuffer,size_packbuffer,mpi_packed,int(tid,MPIISZ),int(tag,MPIISZ),mpi_comm_world, &
                    mpi_status_ignore, rqu, ierr)
   return
 end subroutine mpi_irecv_pack

 function mpi_unpack_int_scalar()
   implicit none
   INTEGER(ISZ) :: mpi_unpack_int_scalar
     call mpi_unpack(packbuffer,size_packbuffer,pack_pos,mpi_unpack_int_scalar,1,mpi_integer,mpi_comm_world,ierr)
   return
 end function mpi_unpack_int_scalar

 function mpi_unpack_int_array(isize)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: isize
   INTEGER(ISZ), DIMENSION(isize) :: mpi_unpack_int_array
     call mpi_unpack(packbuffer,size_packbuffer,pack_pos,mpi_unpack_int_array,int(isize,MPIISZ),mpi_integer,mpi_comm_world,ierr)
   return
 end function mpi_unpack_int_array

 function mpi_unpack_real_scalar()
   implicit none
   REAL(8) :: mpi_unpack_real_scalar
     call mpi_unpack(packbuffer,size_packbuffer,pack_pos,mpi_unpack_real_scalar,1,mpi_double_precision,mpi_comm_world,ierr)
   return
 end function mpi_unpack_real_scalar

 function mpi_unpack_real_array(isize)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: isize
   REAL(8), DIMENSION(isize) :: mpi_unpack_real_array
     call mpi_unpack(packbuffer,size_packbuffer,pack_pos,mpi_unpack_real_array,int(isize,MPIISZ), &
                     mpi_double_precision,mpi_comm_world,ierr)
   return
 end function mpi_unpack_real_array

end module mpirz
#endif
