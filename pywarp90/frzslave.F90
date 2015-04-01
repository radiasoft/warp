!     Last change:  JLV  11 Jun 2002    4:35 pm
#ifdef MPIPARALLEL
#include "top.h"
module mpirz
use Parallel
INCLUDE 'mpif.h'
INTEGER(MPIISZ) :: ierr, length, temp
CHARACTER( LEN = MPI_MAX_ERROR_STRING ) :: message
INTEGER(MPIISZ) :: mpi_status(mpi_status_size),mpirequests(1000),mpireqpnt=0
INTEGER(MPIISZ) :: comm_world_mpiisz
type mpibuffertype
  integer(MPIISZ) :: pack_pos=0
  integer(MPIISZ) :: pack_size=0
  integer*8, allocatable :: buffer(:)
  logical(ISZ) :: l_allocated=.false.
end type mpibuffertype
type(mpibuffertype), dimension(1000), save :: mpibuffers
!  integer*8, allocatable :: packbuffer

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
   MODULE PROCEDURE mpi_pack_logical_3darray
   MODULE PROCEDURE mpi_pack_complex_1darray
   MODULE PROCEDURE mpi_pack_complex_2darray
   MODULE PROCEDURE mpi_pack_complex_3darray
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
  comm_world_mpiisz = comm_world
  call mpi_send(i,1,mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ),comm_world_mpiisz,ierr)
  if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_send_int_scalar: ', message(1:length)
   endif
 END SUBROUTINE mpi_send_int_scalar

 SUBROUTINE mpi_send_int_array(i, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, i(:), tag
  comm_world_mpiisz = comm_world
  call mpi_send(i,SIZE(i),mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ),comm_world_mpiisz,ierr)
  if (ierr /= MPI_SUCCESS) then
     call MPI_Error_string( ierr, message, length, temp )
     write(STDOUT,*) '***** Error in mpi_send_int_array: ', message(1:length)
  endif
 END SUBROUTINE mpi_send_int_array

 function mpi_recv_int_scalar(tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   INTEGER(ISZ) :: mpi_recv_int_scalar
   comm_world_mpiisz = comm_world
   call mpi_recv(mpi_recv_int_scalar,1,mpi_integer,int(tid,MPIISZ),int(tag,MPIISZ), &
        comm_world_mpiisz,mpi_status_ignore,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_recv_int_scalar: ', message(1:length)
   endif
   return
 end function mpi_recv_int_scalar

 function mpi_recv_int_array(isize,tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
   INTEGER(ISZ), DIMENSION(isize) :: mpi_recv_int_array
   comm_world_mpiisz = comm_world
   call mpi_recv(mpi_recv_int_array,int(isize,MPIISZ),mpi_integer,int(tid,MPIISZ), &
        int(tag,MPIISZ),comm_world_mpiisz, mpi_status_ignore,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_recv_int_array: ', message(1:length)
   endif
   return
 end function mpi_recv_int_array

 SUBROUTINE mpi_send_real_scalar(r, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, tag
  REAL(8) :: r
  comm_world_mpiisz = comm_world
  call mpi_send(r,1,mpi_double_precision,int(tid,MPIISZ),int(tag,MPIISZ),comm_world_mpiisz,ierr)
  if (ierr /= MPI_SUCCESS) then
     call MPI_Error_string( ierr, message, length, temp )
     write(STDOUT,*) '***** Error in mpi_send_real_scalar: ', message(1:length)
  endif
 END SUBROUTINE mpi_send_real_scalar

 SUBROUTINE mpi_send_real_array(r, tid, tag)
  IMPLICIT NONE
  INTEGER(ISZ), INTENT(IN) :: TID, tag
  REAL(8), DIMENSION(:) :: r
  if (l_mpiverbose) WRITE(0,*) my_index,'send to ',tid
  comm_world_mpiisz = comm_world
  call mpi_send(r,int(SIZE(r),MPIISZ),mpi_double_precision,int(tid,MPIISZ), &
       int(tag,MPIISZ),comm_world_mpiisz,ierr)
  if (ierr /= MPI_SUCCESS) then
     call MPI_Error_string( ierr, message, length, temp )
     write(STDOUT,*) '***** Error in mpi_send_real_array: ', message(1:length)
  endif
  if (l_mpiverbose) WRITE(0,*) my_index,'sent to ',tid
 END SUBROUTINE mpi_send_real_array

 function mpi_recv_real_scalar(tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   REAL(8) :: mpi_recv_real_scalar
   comm_world_mpiisz = comm_world
   call mpi_recv(mpi_recv_real_scalar,1,mpi_double_precision,int(tid,MPIISZ), &
        int(tag,MPIISZ),comm_world_mpiisz, mpi_status_ignore,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_recv_real_scalar: ', message(1:length)
   endif
   return
 end function mpi_recv_real_scalar

 function mpi_recv_real_array(isize,tid,tag)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
   REAL(8), DIMENSION(isize) :: mpi_recv_real_array
     if (l_mpiverbose) WRITE(0,*) my_index,'recv from ',tid,isize
     comm_world_mpiisz = comm_world
     call mpi_recv(mpi_recv_real_array,int(isize,MPIISZ),mpi_double_precision,int(tid,MPIISZ), &
          int(tag,MPIISZ), comm_world_mpiisz, mpi_status_ignore,ierr)
     if (ierr /= MPI_SUCCESS) then
        call MPI_Error_string( ierr, message, length, temp )
        write(STDOUT,*) '***** Error in mpi_recv_real_array: ', message(1:length)
     endif
     if (l_mpiverbose) WRITE(0,*) my_index,'recvd from ',tid
   return
 end function mpi_recv_real_array

  SUBROUTINE mpi_isend_real_array(r, tid, tag)
    IMPLICIT NONE
    INTEGER(ISZ), INTENT(IN) :: TID, tag
    REAL(8), DIMENSION(:) :: r

!    WRITE(0,*) 'isend to ',tid,size(r)
      comm_world_mpiisz = comm_world
      call mpi_isend(r,int(SIZE(r),MPIISZ),mpi_double_precision,int(tid,MPIISZ), &
           int(tag,MPIISZ),comm_world_mpiisz,mpirequests(mpireqpnt+1),ierr)
      mpireqpnt=mpireqpnt+1
      if (ierr /= MPI_SUCCESS) then
        call MPI_Error_string( ierr, message, length, temp )
        write(STDOUT,*) '***** Error in mpi_isend_real_array: ', message(1:length)
     endif
!    WRITE(0,*) 'isent to ',tid
  END SUBROUTINE mpi_isend_real_array

  function mpi_irecv_real_array(isize,tid,tag)
    implicit none
    INTEGER(ISZ), INTENT(IN) :: isize,tid,tag
    REAL(8), DIMENSION(isize) :: mpi_irecv_real_array
      comm_world_mpiisz = comm_world
!    WRITE(0,*) 'irecv from ',tid
      call mpi_irecv(mpi_irecv_real_array,int(isize,MPIISZ),mpi_double_precision, &
           int(tid,MPIISZ),int(tag,MPIISZ), comm_world_mpiisz,mpirequests(mpireqpnt+1),ierr)
      mpireqpnt=mpireqpnt+1
      if (ierr /= MPI_SUCCESS) then
         call MPI_Error_string( ierr, message, length, temp )
         write(STDOUT,*) '***** Error in mpi_irecv_real_array: ', message(1:length)
      endif
      !    WRITE(0,*) 'irecvd from ',tid
    return
  end function mpi_irecv_real_array

  function mpi_global_compute_real(DATA,op)
    implicit none
    REAL(8) :: DATA, mpi_global_compute_real
    INTEGER(MPIISZ) :: op
    comm_world_mpiisz = comm_world
    call mpi_allreduce(data,mpi_global_compute_real,int(1,MPIISZ),mpi_double_precision, &
         int(op,MPIISZ),comm_world_mpiisz,ierr)
    if (ierr /= MPI_SUCCESS) then
       call MPI_Error_string( ierr, message, length, temp )
       write(STDOUT,*) '***** Error in mpi_global_compute_real: ', message(1:length)
    endif
    return
  end function mpi_global_compute_real

  subroutine mpi_packbuffer_init(isize,ibuf)
    implicit none
    INTEGER(ISZ), INTENT(IN) :: isize,ibuf
    mpibuffers(ibuf)%pack_pos=0
    !   IF(ALLOCATED(mpibuffers(ibuf)%buffer)) then
    IF(mpibuffers(ibuf)%l_allocated) then
       if (mpibuffers(ibuf)%pack_size==isize) return
       DEALLOCATE(mpibuffers(ibuf)%buffer)
    end if
    ALLOCATE(mpibuffers(ibuf)%buffer(8*isize))
    mpibuffers(ibuf)%l_allocated=.true.
    mpibuffers(ibuf)%pack_size=isize
    return
  END subroutine mpi_packbuffer_init

 subroutine mpi_pack_int_scalar(a,ibuf)
   implicit none
   INTEGER(ISZ), INTENT(IN) :: a,ibuf
   comm_world_mpiisz = comm_world
   call mpi_pack(a, 1, mpi_integer, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer),MPIISZ), &
        mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_int_scalar: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_int_scalar

 subroutine mpi_pack_int_array(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), DIMENSION(:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_integer, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer),MPIISZ), &
        mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_int_array: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_int_array

 subroutine mpi_pack_real_scalar(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   REAL(8), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, 1, mpi_double_precision, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer), &
        MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_real_scalar: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_real_scalar

 subroutine mpi_pack_real_1darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   REAL(8), DIMENSION(:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_double_precision, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer), &
        MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_real_1darray: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_real_1darray

 subroutine mpi_pack_real_2darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   REAL(8), DIMENSION(:,:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_double_precision, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer), &
        MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_real_2array: ', message(1:length)
   endif
   if (l_mpiverbose) &
        write(STDOUT,*) 'pack_real', mpibuffers(ibuf)%pack_pos, int(size(mpibuffers(ibuf)%buffer),MPIISZ)
   return
 end subroutine mpi_pack_real_2darray

 subroutine mpi_pack_real_3darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   REAL(8), DIMENSION(:,:,:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_double_precision, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer) &
        ,MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_real_3darray: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_real_3darray

  subroutine mpi_pack_complex_1darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   COMPLEX(8), DIMENSION(:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_double_complex, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer), &
        MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_complex_1darray: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_complex_1darray

 subroutine mpi_pack_complex_2darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   COMPLEX(8), DIMENSION(:,:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_double_complex, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer), &
        MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_complex_2darray: ', message(1:length)
   endif
   if (l_mpiverbose) &
        write(STDOUT,*) 'pack_cmpx', mpibuffers(ibuf)%pack_pos, int(size(mpibuffers(ibuf)%buffer),MPIISZ)     
   return
 end subroutine mpi_pack_complex_2darray

 subroutine mpi_pack_complex_3darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   COMPLEX(8), DIMENSION(:,:,:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
   call mpi_pack(a, SIZE(a), mpi_double_complex, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer), &
        MPIISZ), mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_complex_3darray: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_complex_3darray

 subroutine mpi_pack_logical_3darray(a,ibuf)
   implicit none
   integer(ISZ)::ibuf
   logical(ISZ), DIMENSION(:,:,:), INTENT(IN) :: a
   comm_world_mpiisz = comm_world
     call mpi_pack(a, SIZE(a), mpi_double_precision, mpibuffers(ibuf)%buffer, int(size(mpibuffers(ibuf)%buffer),MPIISZ), &
          mpibuffers(ibuf)%pack_pos, comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_logical_3darray: ', message(1:length)
   endif
   return
 end subroutine mpi_pack_logical_3darray

 subroutine mpi_send_pack(tid,tag,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   comm_world_mpiisz = comm_world
   call mpi_send(mpibuffers(ibuf)%buffer,mpibuffers(ibuf)%pack_pos,mpi_packed,int(tid,MPIISZ),int(tag,MPIISZ), &
        comm_world_mpiisz, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_send_pack: ', message(1:length)
   endif
   return
 end subroutine mpi_send_pack

 subroutine mpi_isend_pack(tid,tag,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   comm_world_mpiisz = comm_world
   call mpi_isend(mpibuffers(ibuf)%buffer,mpibuffers(ibuf)%pack_pos,mpi_packed,int(tid,MPIISZ),int(tag,MPIISZ), &
        comm_world_mpiisz, mpirequests(mpireqpnt+1), ierr)
   mpireqpnt=mpireqpnt+1
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_isend_pack: ', message(1:length)
   endif
   return
 end subroutine mpi_isend_pack

 subroutine mpi_recv_pack(tid,tag,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   comm_world_mpiisz = comm_world
   call mpi_recv(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpi_packed,int(tid,MPIISZ), &
        int(tag,MPIISZ), comm_world_mpiisz, mpi_status_ignore, ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_recv_pack: ', message(1:length)
   endif
   return
 end subroutine mpi_recv_pack

 subroutine mpi_irecv_pack(tid,tag,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: tid,tag
   comm_world_mpiisz = comm_world
   call mpi_irecv(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpi_packed,int(tid,MPIISZ), &
        int(tag,MPIISZ), comm_world_mpiisz, mpi_status_ignore, mpirequests(mpireqpnt+1), ierr)
   mpireqpnt=mpireqpnt+1
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_irecv_pack: ', message(1:length)
   endif
   return
 end subroutine mpi_irecv_pack

 function mpi_unpack_int_scalar(ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ) :: mpi_unpack_int_scalar
   comm_world_mpiisz = comm_world
   call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
        mpi_unpack_int_scalar, 1,mpi_integer,comm_world_mpiisz,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_unpack_int_scalar: ', message(1:length)
   endif
   return
 end function mpi_unpack_int_scalar

 function mpi_unpack_int_array(isize,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: isize
   INTEGER(ISZ), DIMENSION(isize) :: mpi_unpack_int_array
   comm_world_mpiisz = comm_world
   call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
        mpi_unpack_int_array, int(isize,MPIISZ),mpi_integer,comm_world_mpiisz,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_unpack_int_array: ', message(1:length)
   endif
   return
 end function mpi_unpack_int_array

 function mpi_unpack_real_scalar(ibuf)
   implicit none
   integer(ISZ)::ibuf
   REAL(8) :: mpi_unpack_real_scalar
   comm_world_mpiisz = comm_world
   call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
        mpi_unpack_real_scalar, 1,mpi_double_precision,comm_world_mpiisz,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_unpack_real_scalar: ', message(1:length)
   endif
   return
 end function mpi_unpack_real_scalar

 function mpi_unpack_real_array(isize,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: isize
   REAL(8), DIMENSION(isize) :: mpi_unpack_real_array
   comm_world_mpiisz = comm_world
   call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
        mpi_unpack_real_array, int(isize,MPIISZ), mpi_double_precision,comm_world_mpiisz,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_pack_real_array: ', message(1:length)
   endif
   if (l_mpiverbose) &
        write(STDOUT,*) 'unpack_real', mpibuffers(ibuf)%pack_pos, int(size(mpibuffers(ibuf)%buffer),MPIISZ)
   return
 end function mpi_unpack_real_array
 
 function mpi_unpack_complex_array(isize,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: isize
   COMPLEX(8), DIMENSION(isize) :: mpi_unpack_complex_array
   comm_world_mpiisz = comm_world
   call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
        mpi_unpack_complex_array, int(isize,MPIISZ), mpi_double_complex,comm_world_mpiisz,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_unpack_complex_array: ', message(1:length)
   endif
   if (l_mpiverbose) &
        write(STDOUT,*) 'unpack_cmpx', mpibuffers(ibuf)%pack_pos, int(size(mpibuffers(ibuf)%buffer),MPIISZ)
   return
 end function mpi_unpack_complex_array

 function mpi_unpack_logical_array(isize,ibuf)
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: isize
   logical(ISZ), DIMENSION(isize) :: mpi_unpack_logical_array
   comm_world_mpiisz = comm_world
     call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
          mpi_unpack_logical_array, int(isize,MPIISZ), mpi_double_precision,comm_world_mpiisz,ierr)
   return
 end function mpi_unpack_logical_array

 subroutine mpi_waitall_requests()
   implicit none
   call mpi_waitall(mpireqpnt,mpirequests,mpi_status,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in mpi_waitall_requests: ', message(1:length)
   endif
   mpireqpnt=0
 end subroutine mpi_waitall_requests

end module mpirz

 subroutine submpi_unpack_real_array(a,isize,ibuf)
   use mpirz
   implicit none
   integer(ISZ)::ibuf
   INTEGER(ISZ), INTENT(IN) :: isize
   REAL(8), DIMENSION(isize) :: a
   comm_world_mpiisz = comm_world
   call mpi_unpack(mpibuffers(ibuf)%buffer,int(size(mpibuffers(ibuf)%buffer),MPIISZ),mpibuffers(ibuf)%pack_pos, &
        a, int(isize,MPIISZ), mpi_double_precision,comm_world_mpiisz,ierr)
   if (ierr /= MPI_SUCCESS) then
      call MPI_Error_string( ierr, message, length, temp )
      write(STDOUT,*) '***** Error in submpi_unpack_real_array: ', message(1:length)
   endif
   return
 end subroutine submpi_unpack_real_array

#endif
