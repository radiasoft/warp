#include "top.h"
c=============================================================================
c TOPSLAVE.M, version $Revision: 1.5 $, $Date: 2001/10/26 00:39:49 $
c Slave routines related to top package.
c D. P. Grote
c=============================================================================
c===========================================================================
      subroutine zpartbnd_slave(zmmax,zmmin,dz,zgrid)
      use InGen
      use InPart
      use Particles
      use Parallel
      use Databuffers
      real(kind=8):: zmmax,zmmin,dz,zgrid

c  Impose boundary conditions on zp.  Puts particles that exit to the left
c  lower in the arrays and particles that exit to the right higher.  Keeps
c  track of how many left and sends to the appropriate processor.

      real(kind=8):: temp,syslen
      integer(ISZ):: is,ip,ins_init,nps_init,idest
      integer(ISZ):: nsendleft,nsendright,nrecvleft,nrecvright
      integer(ISZ):: left_pe,right_pe
      integer(ISZ):: nquant

c     integer(ISZ):: ii
c     integer(ISZ):: blocklengths(7),displacements(7)
      include "mpif.h"
      integer(ISZ):: mpierror,mpistatus(MPI_STATUS_SIZE)
c     integer(ISZ):: mpinewtype,mpirequest


c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Loop over species
      do is=1,ns

c       --- initialize counter and save particle's start and number
        ip = ins(is)
        ins_init = ins(is)
        nps_init = nps(is)

c       --- Loop over particles and pick out the ones which have ventured
c       --- into the regions of neighboring processors.
        do while (ip < ins(is) + nps(is))

c         --- If particle exited to the left, switch with lowest particle
c         --- in the arrays.
          if (zp(ip)-zgrid < zpslmin(my_index)) then
            idest = ins(is)
            ins(is) = ins(is) + 1
            nps(is) = nps(is) - 1

            temp = xp(ip)
            xp(ip) = xp(idest)
            xp(idest) = temp
            temp = yp(ip)
            yp(ip) = yp(idest)
            yp(idest) = temp
            temp = zp(ip)
            zp(ip) = zp(idest)
            zp(idest) = temp
            temp = uxp(ip)
            uxp(ip) = uxp(idest)
            uxp(idest) = temp
            temp = uyp(ip)
            uyp(ip) = uyp(idest)
            uyp(idest) = temp
            temp = uzp(ip)
            uzp(ip) = uzp(idest)
            uzp(idest) = temp
            temp = gaminv(ip)
            gaminv(ip) = gaminv(idest)
            gaminv(idest) = temp
	    if (npmaxi > 0) then
              temp = pid(ip)
              pid(ip) = pid(idest)
              pid(idest) = temp
	    endif
c         endif

c         --- If particle exited to the right, switch with highest particle
c         --- in the arrays.
          elseif (zp(ip)-zgrid >= zpslmax(my_index)) then
            idest = ins(is) + nps(is) - 1
            nps(is) = nps(is) - 1

            temp = xp(ip)
            xp(ip) = xp(idest)
            xp(idest) = temp
            temp = yp(ip)
            yp(ip) = yp(idest)
            yp(idest) = temp
            temp = zp(ip)
            zp(ip) = zp(idest)
            zp(idest) = temp
            temp = uxp(ip)
            uxp(ip) = uxp(idest)
            uxp(idest) = temp
            temp = uyp(ip)
            uyp(ip) = uyp(idest)
            uyp(idest) = temp
            temp = uzp(ip)
            uzp(ip) = uzp(idest)
            uzp(idest) = temp
            temp = gaminv(ip)
            gaminv(ip) = gaminv(idest)
            gaminv(idest) = temp
	    if (npmaxi > 0) then
              temp = pid(ip)
              pid(ip) = pid(idest)
              pid(idest) = temp
	    endif
            ip = ip - 1
          endif

c         --- advance ip
          ip = ip + 1
        enddo

c       --- Apply appropriate boundary conditions on particles at the axial
c       --- ends of the full mesh.

c       --- Periodic boundary conditions: all particles picked out have
c       --- their z adjusted.
        if (my_index == 0 .and. periinz) then
          syslen = zpslmax(nslaves-1) - zpslmin(0)
          do ip=ins_init,ins(is)-1
            zp(ip) = zp(ip) + syslen
          enddo
        elseif (my_index == nslaves-1 .and. periinz) then
          syslen = zpslmax(nslaves-1) - zpslmin(0)
          do ip=ins(is)+nps(is),ins_init+nps_init-1
            zp(ip) = zp(ip) - syslen
          enddo
c       --- Sticky boundary condition: particles picked out are not passed.
        elseif (my_index == 0 .and. stickyz) then
          nps_init = nps_init - (ins(is) - ins_init)
          ins_init = ins(is)
        elseif (my_index == nslaves-1 .and. stickyz) then
          nps_init = nps(is) + (ins(is) - ins_init)
        endif

c Here is what happens...
c  First, particle data to be sent is copied into temporary buffers.
c  Second, exchange number of particles to be sent in each direction.
c  Third, make sure there is room in particle arrays for incoming data.
c  Fourth, For each direction, create an mpi type which which has the
c          explicit addresses of the places where the particle data goes.
c  Fifth, exchange the data.
c
c Step four is done so that the incoming data does not have to be buffered,
c but can be received directly into the correct memory locations.

c       --- Make sure these are zero in case they are not set below.
        nsendleft = 0
        nsendright = 0

c       --- Number of quantities to be exchanged. Normally it is seven,
c       --- but when npmaxi > 0, then it is 8 since pid is also exchanged.
	nquant = 7
	if (npmaxi > 0) nquant = 8

c       --- Copy particles being sent to the left to buffer1.
        if (periinz .or. my_index > 0) then
          nsendleft = ins(is) - ins_init
          if (nquant*nsendleft > b1size) then
            b1size = nquant*nsendleft
            call gchange("Databuffers",0)
          endif
          if (nsendleft > 0) then
            call copyarry(xp(ins_init),buffer1(1),nsendleft)
            call copyarry(yp(ins_init),buffer1(1+nsendleft),nsendleft)
            call copyarry(zp(ins_init),buffer1(1+2*nsendleft),nsendleft)
            call copyarry(uxp(ins_init),buffer1(1+3*nsendleft),nsendleft)
            call copyarry(uyp(ins_init),buffer1(1+4*nsendleft),nsendleft)
            call copyarry(uzp(ins_init),buffer1(1+5*nsendleft),nsendleft)
            call copyarry(gaminv(ins_init),buffer1(1+6*nsendleft),nsendleft)
	    if (npmaxi > 0)
     &        call copyarry(pid(ins_init),buffer1(1+7*nsendleft),nsendleft)
          endif
        endif

c       --- Copy particles being sent to the right to buffer2.
        if (periinz .or. my_index < nslaves-1) then
          nsendright = ins_init + nps_init - ins(is) - nps(is)
          if (nquant*nsendright > b2size) then
            b2size = nquant*nsendright
            call gchange("Databuffers",0)
          endif
          if (nsendright > 0) then
            ip = ins(is) + nps(is)
            call copyarry(xp(ip),buffer2(1),nsendright)
            call copyarry(yp(ip),buffer2(1+nsendright),nsendright)
            call copyarry(zp(ip),buffer2(1+2*nsendright),nsendright)
            call copyarry(uxp(ip),buffer2(1+3*nsendright),nsendright)
            call copyarry(uyp(ip),buffer2(1+4*nsendright),nsendright)
            call copyarry(uzp(ip),buffer2(1+5*nsendright),nsendright)
            call copyarry(gaminv(ip),buffer2(1+6*nsendright),nsendright)
	    if (npmaxi > 0)
     &        call copyarry(pid(ip),buffer2(1+7*nsendright),nsendright)
          endif
        endif

c       --- Now send the data.

c       --- First send the number of particles to be sent
        call MPI_SENDRECV(nsendleft,1,MPI_INTEGER,left_pe,0,
     &                    nrecvright,1,MPI_INTEGER,right_pe,0,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)
        call MPI_SENDRECV(nsendright,1,MPI_INTEGER,right_pe,0,
     &                    nrecvleft,1,MPI_INTEGER,left_pe,0,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)

c       --- Make sure that there is space for the incoming data.
        call chckpart(is,nrecvleft,nrecvright,.false.)

        if (nquant*nrecvleft > b3size) then
          b3size = nquant*nrecvleft
          call gchange("Databuffers",0)
        endif
        if (nquant*nrecvright > b4size) then
          b4size = nquant*nrecvright
          call gchange("Databuffers",0)
        endif

        call MPI_SENDRECV(buffer2,nquant*nsendright,MPI_DOUBLE_PRECISION,right_pe,0,
     &                    buffer3,nquant*nrecvleft,MPI_DOUBLE_PRECISION,left_pe,0,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)
        call MPI_SENDRECV(buffer1,nquant*nsendleft,MPI_DOUBLE_PRECISION,left_pe,0,
     &                    buffer4,nquant*nrecvright,MPI_DOUBLE_PRECISION,right_pe,0,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)

        ip = ins(is) - nrecvleft
        call copyarry(buffer3(1),xp(ip),nrecvleft)
        call copyarry(buffer3(1+nrecvleft),yp(ip),nrecvleft)
        call copyarry(buffer3(1+2*nrecvleft),zp(ip),nrecvleft)
        call copyarry(buffer3(1+3*nrecvleft),uxp(ip),nrecvleft)
        call copyarry(buffer3(1+4*nrecvleft),uyp(ip),nrecvleft)
        call copyarry(buffer3(1+5*nrecvleft),uzp(ip),nrecvleft)
        call copyarry(buffer3(1+6*nrecvleft),gaminv(ip),nrecvleft)
        if (npmaxi > 0)
     &    call copyarry(buffer3(1+7*nrecvleft),pid(ip),nrecvleft)
        ins(is) = ins(is) - nrecvleft
        nps(is) = nps(is) + nrecvleft

        ip = ins(is) + nps(is)
        call copyarry(buffer4(1),xp(ip),nrecvright)
        call copyarry(buffer4(1+nrecvright),yp(ip),nrecvright)
        call copyarry(buffer4(1+2*nrecvright),zp(ip),nrecvright)
        call copyarry(buffer4(1+3*nrecvright),uxp(ip),nrecvright)
        call copyarry(buffer4(1+4*nrecvright),uyp(ip),nrecvright)
        call copyarry(buffer4(1+5*nrecvright),uzp(ip),nrecvright)
        call copyarry(buffer4(1+6*nrecvright),gaminv(ip),nrecvright)
        if (npmaxi > 0)
     &    call copyarry(buffer4(1+7*nrecvright),pid(ip),nrecvright)
        nps(is) = nps(is) + nrecvright

c What follows (but is commented out) is the preferred version, but
c is causes a core dump for some unknown reason.
c       --- Create new type for receiving particles from left
c       do ii=1,nquant
c         blocklengths(ii) = nrecvleft
c       enddo
c       ip = ins(is) - nrecvleft
c       call MPI_ADDRESS(xp(ip),displacements(1),mpierror)
c       call MPI_ADDRESS(yp(ip),displacements(2),mpierror)
c       call MPI_ADDRESS(zp(ip),displacements(3),mpierror)
c       call MPI_ADDRESS(uxp(ip),displacements(4),mpierror)
c       call MPI_ADDRESS(uyp(ip),displacements(5),mpierror)
c       call MPI_ADDRESS(uzp(ip),displacements(6),mpierror)
c       call MPI_ADDRESS(gaminv(ip),displacements(7),mpierror)
c       if (npmaxi > 0)
c    &    call MPI_ADDRESS(pid(ip),displacements(8),mpierror)
c       call MPI_TYPE_HINDEXED(nquant,blocklengths,displacements,
c    &                         MPI_DOUBLE_PRECISION,mpinewtype,mpierror)
c       call MPI_TYPE_COMMIT(mpinewtype,mpierror)

c       --- Pass particle data to the right
c       call MPI_SENDRECV(buffer2,nquant*nsendright,MPI_DOUBLE_PRECISION,right_pe,0,
c    &                    MPI_BOTTOM,1,mpinewtype,left_pe,0,
c    &                    MPI_COMM_WORLD,mpistatus,mpierror)
c       ins(is) = ins(is) - nrecvleft
c       nps(is) = nps(is) + nrecvleft
c       call MPI_TYPE_FREE(mpinewtype,mpierror)

c       --- Create new type for receiving particles from right
c       do ii=1,nquant
c         blocklengths(ii) = nrecvright
c       enddo
c       ip = ins(is) + nps(is)
c       call MPI_ADDRESS(xp(ip),displacements(1),mpierror)
c       call MPI_ADDRESS(yp(ip),displacements(2),mpierror)
c       call MPI_ADDRESS(zp(ip),displacements(3),mpierror)
c       call MPI_ADDRESS(uxp(ip),displacements(4),mpierror)
c       call MPI_ADDRESS(uyp(ip),displacements(5),mpierror)
c       call MPI_ADDRESS(uzp(ip),displacements(6),mpierror)
c       call MPI_ADDRESS(gaminv(ip),displacements(7),mpierror)
c       if (npmaxi > 0)
c    &    call MPI_ADDRESS(pid(ip),displacements(8),mpierror)
c       call MPI_TYPE_HINDEXED(nquant,blocklengths,displacements,
c    &                         MPI_DOUBLE_PRECISION,mpinewtype,mpierror)
c       call MPI_TYPE_COMMIT(mpinewtype,mpierror)

c       --- Pass particle data to the left
c       call MPI_SENDRECV(buffer1,nquant*nsendleft,MPI_DOUBLE_PRECISION,left_pe,0,
c    &                    MPI_BOTTOM,1,mpinewtype,right_pe,0,
c    &                    MPI_COMM_WORLD,mpistatus,mpierror)
c       nps(is) = nps(is) + nrecvright
c       call MPI_TYPE_FREE(mpinewtype,mpierror)

      enddo

      return
      end
c=============================================================================
      integer(ISZ) function checkzpartbnd(zgrid)
      use InPart
      use Particles
      use Parallel
      real(kind=8):: zgrid
c Check if all particles are within the mesh.  Returns number of particles
c outside the mesh.
      integer(ISZ):: is,ip,nout
      nout = 0
      do is=1,ns
        do ip=ins(is),ins(is)+nps(is)-1
          if (zp(ip)-zgrid < zpslmin(my_index) .or.
     &        zp(ip)-zgrid >= zpslmax(my_index)) then
            nout = nout + 1
          endif
        enddo
      enddo
      checkzpartbnd = nout
      return
      end
c=============================================================================
c=============================================================================
      subroutine sum_mmnts
      use Parallel
      use Moments
      use Z_Moments

c Use reduction routines to sum whole beam moments across all
c of the processors.  It also shares z moment data at PE boundaries.

c     --- temporary for z moments
      real(kind=8):: temp0(NUMZMMNT),temp1(NUMZMMNT,0:1),temp2(NUMZMMNT,0:1)
      integer(ISZ):: im,iz
      integer(ISZ):: left_pe,right_pe
      integer(ISZ):: left_iz,right_iz
      integer(ISZ):: left_nz,right_nz
      real(kind=8):: maxes1(6),mines1(6)
      real(kind=8):: maxes2(6),mines2(6)
      include "mpif.h"
      integer(ISZ):: mpierror,mpirequest,mpistatus(MPI_STATUS_SIZE)

c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Do reduction on whole beam moments.
      temp0 = zmmnts0
      call MPI_ALLREDUCE(temp0,zmmnts0,NUMZMMNT,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)

c     --- Get number of planes to exchange. The default values of zero are
c     --- only applied to prevent kaboom from being called on the first
c     --- and last processors.
      left_nz = 0
      right_nz = 0
      if (my_index > 0) then
        left_iz = 0
        left_nz = izpslave(my_index-1)+nzpslave(my_index-1)-izpslave(my_index)
	left_nz = min(left_nz, 1)
      endif
      if (my_index < nslaves-1) then
        right_iz = izpslave(my_index+1) - izpslave(my_index)
        right_nz = nzpslave(my_index) - right_iz
	right_nz = min(right_nz, 1)
      endif

c     --- Share z moments data at overlapping PE boundaries.
c     --- Send data to the right.
      if (my_index < nslaves-1) then
        do iz=right_iz,right_iz+right_nz
          do im=1,NUMZMMNT
            temp1(im,iz-right_iz) = zmmnts(iz,im)
          enddo
        enddo
        call MPI_ISEND(temp1,(right_nz+1)*NUMZMMNT,MPI_DOUBLE_PRECISION,
     &                 right_pe,0,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
      if (my_index > 0) then
        call MPI_RECV(temp2,(left_nz+1)*NUMZMMNT,MPI_DOUBLE_PRECISION,
     &                left_pe,0,MPI_COMM_WORLD,mpistatus,mpierror)
        do iz=left_iz,left_iz+left_nz
          do im=1,NUMZMMNT
            zmmnts(iz,im) = zmmnts(iz,im) + temp2(im,iz-left_iz)
            temp2(im,iz-left_iz) = zmmnts(iz,im)
          enddo
        enddo
c       --- Send summed data back to the left.
        call MPI_ISEND(temp2,(left_nz+1)*NUMZMMNT,MPI_DOUBLE_PRECISION,
     &                 left_pe,0,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
      if (my_index < nslaves-1) then
        call MPI_RECV(temp1,(right_nz+1)*NUMZMMNT,MPI_DOUBLE_PRECISION,
     &                right_pe,0,MPI_COMM_WORLD,mpistatus,mpierror)
        do iz=right_iz,right_iz+right_nz
          do im=1,NUMZMMNT
            zmmnts(iz,im) = temp1(im,iz-right_iz)
          enddo
        enddo
      endif

c     --- Find global max's and min's
c     --- This is done minimizing the amount of interprocessor communication
      maxes1(1) = xmaxp
      maxes1(2) = ymaxp
      maxes1(3) = zmaxp
      maxes1(4) = vxmaxp
      maxes1(5) = vymaxp
      maxes1(6) = vzmaxp
      mines1(1) = xminp
      mines1(2) = yminp
      mines1(3) = zminp
      mines1(4) = vxminp
      mines1(5) = vyminp
      mines1(6) = vzminp
      call MPI_ALLREDUCE(maxes1,maxes2,6,MPI_DOUBLE_PRECISION,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      call MPI_ALLREDUCE(mines1,mines2,6,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                   MPI_COMM_WORLD,mpierror)
      xmaxp = maxes2(1)
      ymaxp = maxes2(2)
      zmaxp = maxes2(3)
      vxmaxp = maxes2(4)
      vymaxp = maxes2(5)
      vzmaxp = maxes2(6)
      xminp = mines2(1)
      yminp = mines2(2)
      zminp = mines2(3)
      vxminp = mines2(4)
      vyminp = mines2(5)
      vzminp = mines2(6)

      return
      end
c=============================================================================
      subroutine percurr_slave(curr,nzzarr,periinz)
      use Parallel
      integer(ISZ):: nzzarr
      real(kind=8):: curr(0:nzzarr)
      logical(ISZ):: periinz
c Make current periodic and share current at PE boundaries.
c Assumes a maximum grid overlap of 11.

      integer(ISZ):: left_pe,right_pe
      integer(ISZ):: left_iz,right_iz
      integer(ISZ):: left_nz,right_nz
      integer(ISZ):: iz
      real(kind=8):: temp(0:10)
      include "mpif.h"
      integer(ISZ):: mpierror,mpirequest,mpistatus(MPI_STATUS_SIZE)

c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

c     --- Get number of planes to exchange.
      left_iz = 0
      left_nz = 1
      right_iz = nzpslave(my_index)
      right_nz = 1
      if (my_index > 0) then
        left_iz = 0
        left_nz = izpslave(my_index-1)+nzpslave(my_index-1)-izpslave(my_index)
	left_nz = min(left_nz, 1)
      endif
      if (my_index < nslaves-1) then
        right_iz = izpslave(my_index+1) - izpslave(my_index)
        right_nz = nzpslave(my_index) - right_iz
	right_nz = min(right_nz, 1)
      endif

c     --- This should never happen, but make sure that no more than 11
c     --- elements of data are to be exchanged.
      if (left_nz > 11 .or. right_nz > 11) then
        call remark("ERROR: in percurr_slave more than 11 values are to be")
        call remark("       exchanged. This just doesn't work. Call Dave.")
        call kaboom(1)
      endif

c     --- Send the data to the right.
      if (my_index < nslaves-1 .or. periinz) then
        call MPI_ISEND(curr(right_iz),right_nz+1,MPI_DOUBLE_PRECISION,
     &                 right_pe,0,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
      if (my_index > 0 .or. periinz) then
        call MPI_RECV(temp,left_nz+1,MPI_DOUBLE_PRECISION,
     &                left_pe,0,MPI_COMM_WORLD,mpistatus,mpierror)
        do iz=left_iz,left_iz+left_nz
          curr(iz) = curr(iz) + temp(iz-left_iz)
        enddo
        call MPI_ISEND(curr(left_iz),left_nz+1,MPI_DOUBLE_PRECISION,
     &                 left_pe,0,MPI_COMM_WORLD,mpirequest,mpierror)
      endif
      if (my_index < nslaves-1 .or. periinz) then
        call MPI_RECV(curr(right_iz),right_nz+1,MPI_DOUBLE_PRECISION,
     &                right_pe,0,MPI_COMM_WORLD,mpistatus,mpierror)
      endif

      return
      end
c=============================================================================
c=============================================================================
c  Some parallel utility routines
c=============================================================================
c=============================================================================
      subroutine parallelsum(data,ndata)
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Does a sum over all of the processors.  Resulting data is broadcasted
c to all processors.

      real(kind=8):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"

      call MPI_ALLREDUCE(data,data2,ndata,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

      return
      end
c==========================================================================
      subroutine parallelmax(data,ndata)
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Find the max of data over all processors.  Resulting data is broadcasted
c to all processors.

      real(kind=8):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"

      call MPI_ALLREDUCE(data,data2,ndata,MPI_DOUBLE_PRECISION,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

      return
      end
c==========================================================================
      subroutine parallelmin(data,ndata)
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Find the min of data over all processors.  Resulting data is broadcasted
c to all processors.

      real(kind=8):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"

      call MPI_ALLREDUCE(data,data2,ndata,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

      return
      end
c==========================================================================
c==========================================================================
      subroutine parallelbroadcast(data,ndata,root)
      integer(ISZ):: ndata,root
      real(kind=8):: data(ndata)
c Broadcast the data to all processors.

      integer(ISZ):: mpierror
      include "mpif.h"

      call MPI_BCAST(data,ndata,MPI_DOUBLE_PRECISION,
     &               root,MPI_COMM_WORLD,mpierror)

      return
      end
c==========================================================================
c==========================================================================
