#include "top.h"
c=============================================================================
c TOPSLAVE.M, version $Revision: 1.28 $, $Date: 2005/01/12 17:04:58 $
c Slave routines related to top package.
c D. P. Grote
c=============================================================================
c===========================================================================
      subroutine zpartbnd_slave(zmmax,zmmin,dz,zgrid)
      use Subtimerstop
      use GlobalVars
      use InGen
      use InPart
      use Particles
      use Parallel
      use Databuffers
      real(kind=8):: zmmax,zmmin,dz,zgrid

c  Impose boundary conditions on zp.  Puts particles that exit to the left
c  lower in the arrays and particles that exit to the right higher.  Keeps
c  track of how many left and sends to the appropriate processor.

      logical(ISZ):: ldoagain
      real(kind=8):: temp,pidtemp(npidmax),syslen
      integer(ISZ):: is,ip,idest,ipid,ii,nn
      integer(ISZ):: ins_init(ns),nps_init(ns)
      integer(ISZ):: nsendleft(ns),nsendright(ns),nrecvleft(ns),nrecvright(ns)
      integer(ISZ):: nsendleftsum,nsendrightsum,nrecvleftsum,nrecvrightsum
      integer(ISZ):: left_pe,right_pe
      integer(ISZ):: nquant,ierr

c     integer(ISZ):: blocklengths(100),displacements(100)
      include "mpif.h"
      integer(ISZ):: mpierror,mpistatus(MPI_STATUS_SIZE)
c     integer(ISZ):: mpinewtype,mpirequest
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

#ifdef VAMPIR
      call VTSYMDEF(10, "ZBND", "SORT", IERR)
      call VTSYMDEF(20, "ZBND", "CPY1", IERR)
      call VTSYMDEF(30, "ZBND", "SNDN", IERR)
      call VTSYMDEF(40, "ZBND", "CHCK", IERR)
      call VTSYMDEF(50, "ZBND", "SNDP", IERR)
      call VTSYMDEF(60, "ZBND", "CPY2", IERR)
      call VTSYMDEF(70, "ZBND", "PLOR", IERR)
      call VTBEGIN(10,IERR)
#endif

c     --- Set ldoagain flag to false. Later on, particles received are checked
c     --- if they are within the bounds and if any are not, ldoagain is set
c     --- to true.
      ldoagain = .false.

c     --- make sure that particles marked as lost in this routine will be
c     --- cleaned by clearpart.
      if(clearlostpart==0) clearlostpart = 1
      
c     --- calculate numbers of the processors to the left and to the right,
c     --- wrapping around at the ends.
      left_pe = mod(nslaves+my_index-1,nslaves)
      right_pe = mod(my_index+1,nslaves)

      ldoagain = .true.
      do while (ldoagain)
        ldoagain = .false.

c       --- Loop over species
        do is=1,ns

c         --- initialize counter and save particle's start and number
          ip = ins(is)
          ins_init(is) = ins(is)
          nps_init(is) = nps(is)

c         --- Loop over particles and pick out the ones which have ventured
c         --- into the regions of neighboring processors.
          do while (ip < ins(is) + nps(is))

c           --- If particle exited to the left, switch with lowest particle
c           --- in the arrays.
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
              if (npid > 0) then
                pidtemp = pid(ip,:)
                pid(ip,:) = pid(idest,:)
                pid(idest,:) = pidtemp
              endif
c           endif

c           --- If particle exited to the right, switch with highest particle
c           --- in the arrays.
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
              if (npid > 0) then
                pidtemp = pid(ip,:)
                pid(ip,:) = pid(idest,:)
                pid(idest,:) = pidtemp
              endif
              ip = ip - 1
            endif

c           --- advance ip
            ip = ip + 1
          enddo

c         --- Apply appropriate boundary conditions on particles at the axial
c         --- ends of the full mesh.

c         --- Periodic boundary conditions: all particles picked out have
c         --- their z adjusted.
          if (my_index == 0 .and. pbound0==periodic) then
            syslen = zpslmax(nslaves-1) - zpslmin(0)
!$OMP PARALLEL DO
            do ip=ins_init(is),ins(is)-1
              zp(ip) = zp(ip) + syslen
            enddo
!$OMP END PARALLEL DO
          elseif (my_index == nslaves-1 .and. pboundnz==periodic) then
            syslen = zpslmax(nslaves-1) - zpslmin(0)
!$OMP PARALLEL DO
            do ip=ins(is)+nps(is),ins_init(is)+nps_init(is)-1
              zp(ip) = zp(ip) - syslen
            enddo
!$OMP END PARALLEL DO
c         --- Sticky boundary condition: particles picked out are not passed.
c         --- Instead, they are marked as lost particles and ins and nps are
c         --- reset to include them.
          elseif (my_index == 0 .and. pbound0==absorb) then
            gaminv(ins_init(is):ins(is)-1) = 0.
            nps(is) = nps(is) + (ins(is) - ins_init(is))
            ins(is) = ins_init(is)
          elseif (my_index == nslaves-1 .and. pboundnz==absorb) then
            gaminv(ins(is)+nps(is):ins_init(is)+nps_init(is)-1) = 0.
            nps(is) = ins_init(is) + nps_init(is) - ins(is)
          endif
c       --- End of loop over species
        enddo
#ifdef VAMPIR
      call VTEND(10,IERR)
      call VTBEGIN(20,IERR)
#endif

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
        nsendleftsum = 0
        nsendrightsum = 0

c       --- Number of quantities to be exchanged. Normally it is seven,
c       --- but when npid>0, then it is more since pid is also exchanged.
        nquant = 7
        if (npid > 0) nquant = nquant + npid

c       --- Copy particles being sent to the left to buffer1.
        if (pbound0==periodic .or. my_index > 0) then
          nsendleft = ins - ins_init
          nsendleftsum = sum(nsendleft)
          if (nquant*nsendleftsum > b1size) then
            b1size = nquant*nsendleftsum
            call gchange("Databuffers",0)
          endif
          if (nsendleftsum > 0) then
            ii = 1
            do is=1,ns
              ip = ins_init(is)
              nn = nsendleft(is)
              if (nn > 0) then
                buffer1(ii     :ii+  nn-1) = xp(ip:ip+nn-1)
                buffer1(ii+  nn:ii+2*nn-1) = yp(ip:ip+nn-1)
                buffer1(ii+2*nn:ii+3*nn-1) = zp(ip:ip+nn-1)
                buffer1(ii+3*nn:ii+4*nn-1) = uxp(ip:ip+nn-1)
                buffer1(ii+4*nn:ii+5*nn-1) = uyp(ip:ip+nn-1)
                buffer1(ii+5*nn:ii+6*nn-1) = uzp(ip:ip+nn-1)
                buffer1(ii+6*nn:ii+7*nn-1) = gaminv(ip:ip+nn-1)
                do ipid=1,npid
                  buffer1(ii+(6+ipid)*nn:ii+(7+ipid)*nn-1) = pid(ip:ip+nn-1,ipid)
                enddo
                ii = ii + nn*nquant
              endif
            enddo
          endif
        endif

c       --- Copy particles being sent to the right to buffer2.
        if (pboundnz==periodic .or. my_index < nslaves-1) then
          nsendright = ins_init + nps_init - ins - nps
          nsendrightsum = sum(nsendright)
          if (nquant*nsendrightsum > b2size) then
            b2size = nquant*nsendrightsum
            call gchange("Databuffers",0)
          endif
          if (nsendrightsum > 0) then
            ii = 1
            do is=1,ns
              ip = ins(is) + nps(is)
              nn = nsendright(is)
              if (nn > 0) then
                buffer2(ii     :ii+  nn-1) = xp(ip:ip+nn-1)
                buffer2(ii+  nn:ii+2*nn-1) = yp(ip:ip+nn-1)
                buffer2(ii+2*nn:ii+3*nn-1) = zp(ip:ip+nn-1)
                buffer2(ii+3*nn:ii+4*nn-1) = uxp(ip:ip+nn-1)
                buffer2(ii+4*nn:ii+5*nn-1) = uyp(ip:ip+nn-1)
                buffer2(ii+5*nn:ii+6*nn-1) = uzp(ip:ip+nn-1)
                buffer2(ii+6*nn:ii+7*nn-1) = gaminv(ip:ip+nn-1)
                do ipid=1,npid
                  buffer2(ii+(6+ipid)*nn:ii+(7+ipid)*nn-1) = pid(ip:ip+nn-1,ipid)
                enddo
                ii = ii + nn*nquant
              endif
            enddo
          endif
        endif

#ifdef VAMPIR
      call VTEND(20,IERR)
      call VTBEGIN(30,IERR)
#endif
c       --- Now send the data.

c       --- First send the number of particles to be sent
        call MPI_SENDRECV(nsendleft,ns,MPI_INTEGER,left_pe,20,
     &                    nrecvright,ns,MPI_INTEGER,right_pe,20,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)
        call MPI_SENDRECV(nsendright,ns,MPI_INTEGER,right_pe,20,
     &                    nrecvleft,ns,MPI_INTEGER,left_pe,20,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)

#ifdef VAMPIR
      call VTEND(30,IERR)
      call VTBEGIN(40,IERR)
#endif
c       --- Make sure that there is space for the incoming data.
        do is=1,ns
          call chckpart(is,nrecvleft(is),nrecvright(is),.false.)
        enddo
#ifdef VAMPIR
      call VTEND(40,IERR)
      call VTBEGIN(50,IERR)
#endif

        nrecvleftsum = sum(nrecvleft)
        nrecvrightsum = sum(nrecvright)
        if (nquant*nrecvleftsum > b3size) then
          b3size = nquant*nrecvleftsum
          call gchange("Databuffers",0)
        endif
        if (nquant*nrecvrightsum > b4size) then
          b4size = nquant*nrecvrightsum
          call gchange("Databuffers",0)
        endif

        call MPI_SENDRECV(buffer2,nquant*nsendrightsum,MPI_DOUBLE_PRECISION,right_pe,20,
     &                    buffer3,nquant*nrecvleftsum,MPI_DOUBLE_PRECISION,left_pe,20,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)
        call MPI_SENDRECV(buffer1,nquant*nsendleftsum,MPI_DOUBLE_PRECISION,left_pe,20,
     &                    buffer4,nquant*nrecvrightsum,MPI_DOUBLE_PRECISION,right_pe,20,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)

#ifdef VAMPIR
      call VTEND(50,IERR)
      call VTBEGIN(60,IERR)
#endif
        if (nrecvleftsum > 0) then
          ii = 1
          do is=1,ns
            nn = nrecvleft(is)
            if (nn > 0) then
              ip = ins(is) - nn
              xp(ip:ip+nn-1)  = buffer3(ii     :ii+  nn-1)
              yp(ip:ip+nn-1)  = buffer3(ii+  nn:ii+2*nn-1)
              zp(ip:ip+nn-1)  = buffer3(ii+2*nn:ii+3*nn-1)
              uxp(ip:ip+nn-1) = buffer3(ii+3*nn:ii+4*nn-1)
              uyp(ip:ip+nn-1) = buffer3(ii+4*nn:ii+5*nn-1)
              uzp(ip:ip+nn-1) = buffer3(ii+5*nn:ii+6*nn-1)
              gaminv(ip:ip+nn-1) = buffer3(ii+6*nn:ii+7*nn-1)
              do ipid=1,npid
                pid(ip:ip+nn-1,ipid) = buffer3(ii+(6+ipid)*nn:ii+(7+ipid)*nn-1)
              enddo
              ii = ii + nquant*nn
              ins(is) = ins(is) - nn
              nps(is) = nps(is) + nn
              if (maxval(zp(ip:ip+nn-1))-zgrid > zpslmax(my_index)) then
                ldoagain = .true.
              endif
            endif
          enddo
        endif

        if (nrecvrightsum > 0) then
          ii = 1
          do is=1,ns
            nn = nrecvright(is)
            if (nn > 0) then
              ip = ins(is) + nps(is)
              xp(ip:ip+nn-1)  = buffer4(ii     :ii+1*nn-1)
              yp(ip:ip+nn-1)  = buffer4(ii+  nn:ii+2*nn-1)
              zp(ip:ip+nn-1)  = buffer4(ii+2*nn:ii+3*nn-1)
              uxp(ip:ip+nn-1) = buffer4(ii+3*nn:ii+4*nn-1)
              uyp(ip:ip+nn-1) = buffer4(ii+4*nn:ii+5*nn-1)
              uzp(ip:ip+nn-1) = buffer4(ii+5*nn:ii+6*nn-1)
              gaminv(ip:ip+nn-1) = buffer4(ii+6*nn:ii+7*nn-1)
              do ipid=1,npid
                pid(ip:ip+nn-1,ipid) = buffer4(ii+(6+ipid)*nn:ii+(7+ipid)*nn-1)
              enddo
              ii = ii + nquant*nn
              nps(is) = nps(is) + nn
              if (minval(zp(ip:ip+nn-1))-zgrid < zpslmin(my_index)) then
                ldoagain = .true.
              endif
            endif
          enddo
        endif

#ifdef VAMPIR
      call VTEND(60,IERR)
      call VTBEGIN(70,IERR)
#endif
        call parallellor(ldoagain)
#ifdef VAMPIR
      call VTEND(70,IERR)
#endif

      enddo

c What follows (but is commented out) is the preferred version, but
c is causes a core dump for some unknown reason.
c Now it needs to be adjusted to send all species at once.
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
c       do ipid=1,npid
c         call MPI_ADDRESS(pid(ip,ipid),displacements(7+ipid),mpierror)
c       enddo
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
c       do ipid=1,npid
c         call MPI_ADDRESS(pid(ip,ipid),displacements(7+ipid),mpierror)
c       enddo
c       call MPI_TYPE_HINDEXED(nquant,blocklengths,displacements,
c    &                         MPI_DOUBLE_PRECISION,mpinewtype,mpierror)
c       call MPI_TYPE_COMMIT(mpinewtype,mpierror)

c       --- Pass particle data to the left
c       call MPI_SENDRECV(buffer1,nquant*nsendleft,MPI_DOUBLE_PRECISION,left_pe,0,
c    &                    MPI_BOTTOM,1,mpinewtype,right_pe,0,
c    &                    MPI_COMM_WORLD,mpistatus,mpierror)
c       nps(is) = nps(is) + nrecvright
c       call MPI_TYPE_FREE(mpinewtype,mpierror)

!$OMP MASTER
      if (ltoptimesubs) timezpartbnd_slave = timezpartbnd_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      subroutine reorgparticles_parallel()
      use Subtimerstop
      use InPart
      use Picglb
      use Particles
      use Parallel

c Do all to all communication to pass particles to where they belong. This
c sorts through all of the particles first, then does all of the communication
c at once.

      integer(ISZ),allocatable:: zprocs(:)
      integer(ISZ):: nzp,izp,ip,izp1,izp2
      real(kind=8):: wzp,dzp,zpmin,zpmax
      real(kind=8),allocatable:: ptemp(:,:),precv(:,:)
      integer(ISZ),allocatable:: pprocs(:)
      integer(ISZ):: pins(0:ns-1,0:nslaves-1),pnps(0:ns-1,0:nslaves-1)
      integer(ISZ):: piii(0:ns-1,0:nslaves-1)
      integer(ISZ):: inrecv(0:ns-1,0:nslaves-1),nprecv(0:ns-1,0:nslaves-1)
      integer(ISZ):: pnp,pnpmax,ncoords,ii,is,pip,ip1,ip2,iitemp,ii1,ii2
      real(kind=8):: zz
      integer(ISZ):: mpierror
      integer(ISZ):: sendcounts(0:nslaves-1),sdispls(0:nslaves-1)
      integer(ISZ):: recvcounts(0:nslaves-1),rdispls(0:nslaves-1)
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- First, create z-array that specifies which processors owns which
c     --- location. This is used to make the sorting of particles to
c     --- processors efficient. Using this array avoids having to loop
c     --- over every processor to find which one the particle is in.
      zpmin = zpslmin(0)
      zpmax = zpslmax(nslaves-1)
      dzp = minval(zpslmax-zpslmin)/10.
      nzp = min(100000,int((zpmax-zpmin)/dzp))
      dzp = (zpmax - zpmin)/nzp
      allocate(zprocs(0:nzp))
      do ip=0,nslaves-1
        izp1 = nzp - int((zpmax - zpslmin(ip))/dzp)
        izp2 = int((zpslmax(ip) - zpmin)/dzp)
        zprocs(izp1:izp2) = ip
      enddo
      zprocs(0) = 0
      zprocs(nzp) = nslaves-1

c     --- Next, create temporary arrays to hold the sorted particles.
      pnp = sum(nps)
      ncoords = 7 + npid
      pnpmax = max(1,pnp)
      allocate(ptemp(0:ncoords-1,pnpmax),pprocs(pnpmax))

c     --- Loop through the particles, finding to which processor each particle
c     --- belongs.
      pip = 0
      pnps = 0
      do is=0,ns-1
        do ii=ins(is+1),ins(is+1)+nps(is+1)-1
          zz = zp(ii) - zbeam
          izp = int((zz - zpmin)/dzp)
          ip1 = zprocs(izp)
          ip2 = zprocs(izp+1)
          pip = pip + 1
          do ip=ip1,ip2
            if (zpslmin(ip) <= zz .and. zz < zpslmax(ip)) then
              pprocs(pip) = ip
              pnps(is,ip) = pnps(is,ip) + 1
              exit
            endif
          enddo
        enddo
      enddo

c     --- Calculate starting index of particle groups for each species and
c     --- processor
      pins(0,0) = 1
      do ip=0,nslaves-1
        do is=0,ns-1
          if (is == 0) then
            if (ip == 0) then
              pins(is,ip) = 1
            else
              pins(is,ip) = pins(ns-1,ip-1) + pnps(ns-1,ip-1)
            endif
          else
            pins(is,ip) = pins(is-1,ip) + pnps(is-1,ip)
          endif
        enddo
      enddo

c     --- Copy sorted particle data into the temporary array
      pip = 0
      piii = pins
      do is=0,ns-1
        do ii=ins(is+1),ins(is+1)+nps(is+1)-1
          pip = pip + 1
          ip = pprocs(pip)
          iitemp = piii(is,ip)
          piii(is,ip) = piii(is,ip) + 1
          ptemp(0,iitemp) = xp(ii)
          ptemp(1,iitemp) = yp(ii)
          ptemp(2,iitemp) = zp(ii)
          ptemp(3,iitemp) = uxp(ii)
          ptemp(4,iitemp) = uyp(ii)
          ptemp(5,iitemp) = uzp(ii)
          ptemp(6,iitemp) = gaminv(ii)
          if (npid > 0) then
            ptemp(7:,iitemp) = pid(ii,:)
          endif
        enddo
      enddo

c     --- Calculate pins for data set sent to each processor relative to 0
      do ip=0,nslaves-1
        pins(:,ip) = pins(:,ip) - pins(0,ip)
      enddo

c     --- Exchange the number of particles which are to be sent to each
c     --- processor
      call MPI_ALLTOALL(pins,ns,MPI_INTEGER,inrecv,ns,MPI_INTEGER,
     &                  MPI_COMM_WORLD,mpierror)
      call MPI_ALLTOALL(pnps,ns,MPI_INTEGER,nprecv,ns,MPI_INTEGER,
     &                  MPI_COMM_WORLD,mpierror)

c     --- Resize particle arrays if neccesary. If there is more room available
c     --- than incoming particles, then spread particles out by giving each
c     --- species a fraction equal to its number of particles over the total.
c     --- If there is not enough room, than add extra space and spread it out
c     --- among the species the same way. Ensure that npmax == sum(np_s) so
c     --- that if more space is not needed, an unneccesary realloction is not
c     --- done.
      pnp = sum(nprecv)
      np_s = sum(nprecv,2)
      if (pnp > npmax) npmax = int(1.1*pnp)
      if (pnp > 0) then
        np_s = int(np_s*npmax/pnp)
      else
        np_s = npmax/ns
      endif
      np_s(ns) = npmax - sum(np_s(1:ns-1))
      call alotpart()
      
c     --- Create space for incoming data
      pnpmax = max(1,pnp)
      allocate(precv(0:ncoords-1,pnpmax))

c     --- Do the communication
      sendcounts = sum(pnps,1)*ncoords
      recvcounts = sum(nprecv,1)*ncoords
      sdispls(0) = 0
      rdispls(0) = 0
      do ip=1,nslaves-1
        sdispls(ip) = sdispls(ip-1) + sendcounts(ip-1)
        rdispls(ip) = rdispls(ip-1) + recvcounts(ip-1)
      enddo
      call MPI_ALLTOALLV(ptemp,sendcounts,sdispls,MPI_DOUBLE_PRECISION,
     &                   precv,recvcounts,rdispls,MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD,mpierror)

c     --- Calculate inrecv for data set sent to each processor relative
c     --- to rdispls
      do ip=0,nslaves-1
        inrecv(:,ip) = inrecv(:,ip) + rdispls(ip)/ncoords + 1
      enddo

c     --- Now copy the data into the particle arrays
      nps = 0
      do is=0,ns-1
        do ip=0,nslaves-1
          if (nprecv(is,ip) > 0) then
            ip1 = inrecv(is,ip)
            ip2 = ip1 + nprecv(is,ip) - 1
            ii1 = ins(is+1) + nps(is+1)
            ii2 = ii1 + nprecv(is,ip) - 1
            nps(is+1) = nps(is+1) + nprecv(is,ip)
            xp(ii1:ii2) = precv(0,ip1:ip2)
            yp(ii1:ii2) = precv(1,ip1:ip2)
            zp(ii1:ii2) = precv(2,ip1:ip2)
            uxp(ii1:ii2) = precv(3,ip1:ip2)
            uyp(ii1:ii2) = precv(4,ip1:ip2)
            uzp(ii1:ii2) = precv(5,ip1:ip2)
            gaminv(ii1:ii2) = precv(6,ip1:ip2)
            do ii=1,npid
              pid(ii1:ii2,ii) = precv(6+ii,ip1:ip2)
            enddo
          endif
        enddo
      enddo

c     --- Free work space
      deallocate(zprocs,ptemp,precv,pprocs)

c     --- Done!

!$OMP MASTER
      if (ltoptimesubs) timereorgparticles_parallel = timereorgparticles_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      integer(ISZ) function checkzpartbnd(zgrid)
      use Subtimerstop
      use InPart
      use Particles
      use Parallel
      real(kind=8):: zgrid
c Check if all particles are within the mesh.  Returns number of particles
c outside the mesh.
      integer(ISZ):: is,ip,nout
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      nout = 0
      do is=1,ns
!$OMP PARALLEL DO REDUCTION(+:nout)
        do ip=ins(is),ins(is)+nps(is)-1
          if (zp(ip)-zgrid < zpslmin(my_index) .or.
     &        zp(ip)-zgrid >= zpslmax(my_index)) then
            nout = nout + 1
          endif
        enddo
      enddo
!$OMP END PARALLEL DO
      checkzpartbnd = nout

!$OMP MASTER
      if (ltoptimesubs) timecheckzpartbnd = timecheckzpartbnd + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      subroutine parallel_sum_mmnts(zmmnts0,zmmnts)
      use Subtimerstop
      use Parallel
      use Moments
      use Z_Moments
      real(kind=8):: zmmnts0(NUMZMMNT,0:nszmmnt)
      real(kind=8):: zmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt)

c Use reduction routines to sum whole beam moments across all
c of the processors.  It also shares z moment data at PE boundaries.

c     --- temporary for z moments
      real(kind=8):: ztemp0(NUMZMMNT,0:nszmmnt)
      real(kind=8):: ztemp(0:nzmmnt,NUMZMMNT,0:nszmmnt)
      real(kind=8):: maxes1(6,0:nszmmnt),mines1(6,0:nszmmnt)
      real(kind=8):: maxes2(6,0:nszmmnt),mines2(6,0:nszmmnt)
      include "mpif.h"
      integer(ISZ):: mpierror
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- Do reduction on whole beam moments.
      ztemp0 = zmmnts0
      call MPI_ALLREDUCE(ztemp0,zmmnts0,NUMZMMNT*(1+nszmmnt),
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)

c     --- Do reduction on beam z moments.
      ztemp = zmmnts
      call MPI_ALLREDUCE(ztemp,zmmnts,(1+nzmmnt)*NUMZMMNT*(1+nszmmnt),
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)

c     --- Find global max's and min's
c     --- This is done minimizing the amount of interprocessor communication
      maxes1(1,:) = xmaxp
      maxes1(2,:) = ymaxp
      maxes1(3,:) = zmaxp
      maxes1(4,:) = vxmaxp
      maxes1(5,:) = vymaxp
      maxes1(6,:) = vzmaxp
      mines1(1,:) = xminp
      mines1(2,:) = yminp
      mines1(3,:) = zminp
      mines1(4,:) = vxminp
      mines1(5,:) = vyminp
      mines1(6,:) = vzminp
      call MPI_ALLREDUCE(maxes1,maxes2,6*(1+nszmmnt),MPI_DOUBLE_PRECISION,
     &                   MPI_MAX,MPI_COMM_WORLD,mpierror)
      call MPI_ALLREDUCE(mines1,mines2,6*(1+nszmmnt),MPI_DOUBLE_PRECISION,
     &                   MPI_MIN,MPI_COMM_WORLD,mpierror)
      xmaxp = maxes2(1,:)
      ymaxp = maxes2(2,:)
      zmaxp = maxes2(3,:)
      vxmaxp = maxes2(4,:)
      vymaxp = maxes2(5,:)
      vzmaxp = maxes2(6,:)
      xminp = mines2(1,:)
      yminp = mines2(2,:)
      zminp = mines2(3,:)
      vxminp = mines2(4,:)
      vyminp = mines2(5,:)
      vzminp = mines2(6,:)

!$OMP MASTER
      if (ltoptimesubs) timeparallel_sum_mmnts = timeparallel_sum_mmnts + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine parallel_sum_temperature()
      use Subtimerstop
      use Temperatures
      integer(ISZ):: nn
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c Do a global sum of the data used to calculate local temperatures.

      nn = (1+nxtslices)*(1+nytslices)*nztslices
      call parallelsumrealarray(pnumt,nn)
      call parallelsumrealarray(pnumtw,nn)
      call parallelsumrealarray(vxbart,nn)
      call parallelsumrealarray(vybart,nn)
      call parallelsumrealarray(vzbart,nn)
      call parallelsumrealarray(vxsqbart,nn)
      call parallelsumrealarray(vysqbart,nn)
      call parallelsumrealarray(vzsqbart,nn)

      if (l_temp_rmcorrelations) then

        nn = (1+nxtslicesc)*(1+nytslicesc)*nztslicesc
        call parallelsumrealarray(xbart,nn)
        call parallelsumrealarray(ybart,nn)
        call parallelsumrealarray(zbart,nn)
        call parallelsumrealarray(xsqbart,nn)
        call parallelsumrealarray(ysqbart,nn)
        call parallelsumrealarray(zsqbart,nn)
        call parallelsumrealarray(xvxbart,nn)
        call parallelsumrealarray(yvybart,nn)
        call parallelsumrealarray(zvzbart,nn)

      endif

!$OMP MASTER
      if (ltoptimesubs) timeparallel_sum_temperature = timeparallel_sum_temperature + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
c  Some parallel utility routines
c=============================================================================
c=============================================================================
      subroutine parallelsumrealarray(data,ndata)
      use Subtimerstop
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Does a sum over all of the processors.  Resulting data is broadcast
c to all processors.

      real(kind=8):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,data2,ndata,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

!$OMP MASTER
      if (ltoptimesubs) timeparallelsumrealarray = timeparallelsumrealarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
      subroutine parallelsumintegerarray(data,ndata)
      use Subtimerstop
      integer(ISZ):: ndata
      integer(ISZ):: data(ndata)
c Does a sum over all of the processors.  Resulting data is broadcast
c to all processors.

      integer(ISZ):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,data2,ndata,MPI_INTEGER,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

!$OMP MASTER
      if (ltoptimesubs) timeparallelsumintegerarray = timeparallelsumintegerarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallelmaxrealarray(data,ndata)
      use Subtimerstop
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Find the max of data over all processors.  Resulting data is broadcast
c to all processors.

      real(kind=8):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,data2,ndata,MPI_DOUBLE_PRECISION,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

!$OMP MASTER
      if (ltoptimesubs) timeparallelmaxrealarray = timeparallelmaxrealarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallelminrealarray(data,ndata)
      use Subtimerstop
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Find the min of data over all processors.  Resulting data is broadcast
c to all processors.

      real(kind=8):: data2(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,data2,ndata,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

!$OMP MASTER
      if (ltoptimesubs) timeparallelminrealarray = timeparallelminrealarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallellor(data)
      use Subtimerstop
      logical(ISZ):: data
c Logical OR of data over all processors.  Resulting data is broadcast
c to all processors.

      logical(ISZ):: data2
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,data2,1,MPI_LOGICAL,MPI_LOR,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2

!$OMP MASTER
      if (ltoptimesubs) timeparallellor = timeparallellor + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallelbroadcastrealarray(data,ndata,root)
      use Subtimerstop
      integer(ISZ):: ndata,root
      real(kind=8):: data(ndata)
c Broadcast the data to all processors.

      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_BCAST(data,ndata,MPI_DOUBLE_PRECISION,
     &               root,MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (ltoptimesubs) timeparallelbroadcastrealarray = timeparallelbroadcastrealarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallelbarrier()
      use Subtimerstop
c Calls MPI barrier
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (ltoptimesubs) timeparallelbarrier = timeparallelbarrier + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallelnonzerorealarray(data,ndata)
      use Subtimerstop
      integer(ISZ):: ndata
      real(kind=8):: data(ndata)
c Find the nonzero value of data over all processors. This assumes that the
c non-zero values for each index are the same for all processors.
c Resulting data is broadcast to all processors.

      real(kind=8):: dmax(ndata),dmin(ndata)
      integer(ISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,dmax,ndata,MPI_DOUBLE_PRECISION,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      call MPI_ALLREDUCE(data,dmin,ndata,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                   MPI_COMM_WORLD,mpierror)
      where (dmax /= 0)
        data = dmax
      elsewhere
        data = dmin
      endwhere

!$OMP MASTER
      if (ltoptimesubs) timeparallelnonzerorealarray = timeparallelnonzerorealarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
