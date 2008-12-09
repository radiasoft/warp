#include "top.h"
c=============================================================================
c TOPSLAVE.M, version $Revision: 1.59 $, $Date: 2008/12/09 20:02:37 $
c Slave routines related to top package.
c D. P. Grote
c=============================================================================
c===========================================================================
c===========================================================================
      subroutine particleboundaries_parallel(axis,pgroup,js1,js2,
     &                                       zmmax,zmmin,
     &                                       zpmaxlocal,zpminlocal,
     &                                       pboundleft,pboundright,
     &                                       lcountaslost,labs,lrz)
      use ParticleGroupmodule
      use Subtimerstop
      use GlobalVars
      use InGen,Only: clearlostpart
c     use InPart
      use Parallel
      use Databuffers
      use Subcycling, Only: ndtstorho,zgridndts
      integer(ISZ):: axis
      type(ParticleGroup):: pgroup
      integer(ISZ):: js1,js2
      real(kind=8):: zmmax,zmmin
      real(kind=8):: zpmaxlocal,zpminlocal
      integer(ISZ):: pboundleft,pboundright
      logical(ISZ):: lcountaslost,labs,lrz

c Impose boundary conditions baesd on zz. zz can be any coordinate.
c Note that zz must be on of the arrays in pgroup so that the rearranging
c below will affect it.
c axis is the axis of zz. 0 is x, 1 is y, 2 is z.
c Puts particles that exit to the left lower in the arrays and particles
c that exit to the right higher.  Keeps track of how many left and sends
c to the appropriate processor.

      logical(ISZ):: ldoagain
      logical(ISZ):: isleftmostproc,isrightmostproc
      real(kind=8),pointer:: zz(:),uz(:)
      real(kind=8):: ztemp,temp,pidtemp(pgroup%npid),syslen,zgrid,rr,rrin
      integer(MPIISZ):: is1,is2,is,ip,idest,ipid,ii,il,ir,nn,indts,i1,i2
      integer(MPIISZ):: ins_init(pgroup%ns),nps_init(pgroup%ns)
      integer(MPIISZ):: nsendleft(pgroup%ns),nsendright(pgroup%ns)
      integer(MPIISZ):: nrecvleft(pgroup%ns),nrecvright(pgroup%ns)
      integer(MPIISZ):: nsendleftsum,nsendrightsum,nrecvleftsum,nrecvrightsum
      integer(MPIISZ):: left_pe,right_pe
      integer(MPIISZ):: nquant,ierr

      include "mpif.h"
      integer(MPIISZ):: mpierror,mpistatus(MPI_STATUS_SIZE)
c     integer(MPIISZ),allocatable:: blocklengths(:),displacements(:)
c     integer(MPIISZ):: mpinewtype,mpirequest
      integer(MPIISZ):: messid = 20

      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- If this dimension has zero length, then skip the boundary
c     --- conditions since this will be an ignorable dimension.
      if (zmmin == zmmax) return
      if (zpminlocal == zpmaxlocal) return

#ifdef VAMPIR
c     call VTSYMDEF(10, "ZBND", "SORT", IERR)
c     call VTSYMDEF(20, "ZBND", "CPY1", IERR)
c     call VTSYMDEF(30, "ZBND", "SNDN", IERR)
c     call VTSYMDEF(40, "ZBND", "CHCK", IERR)
c     call VTSYMDEF(50, "ZBND", "SNDP", IERR)
c     call VTSYMDEF(60, "ZBND", "CPY2", IERR)
c     call VTSYMDEF(70, "ZBND", "PLOR", IERR)
c     call VTBEGIN(10,IERR)
#endif

      is1 = js1 + 1
      is2 = js2 + 1

c     --- Get the neighboring processors
      left_pe = procneighbors(0,axis)
      right_pe = procneighbors(1,axis)
      isleftmostproc = (iprocgrid(axis) == 0)
      isrightmostproc = (iprocgrid(axis) == nprocgrid(axis) - 1)

c     --- Setup pointer to correct data coordinates
      if (axis == 0) then
        zz => pgroup%xp
        uz => pgroup%uxp
      else if (axis == 1) then
        zz => pgroup%yp
        uz => pgroup%uyp
      else if (axis == 2) then
        zz => pgroup%zp
        uz => pgroup%uzp
      endif

      ldoagain = .true.
      do while (ldoagain)
c       --- Set ldoagain flag to false. Later on, particles received are checked
c       --- if they are within the bounds and if any are not, ldoagain is set
c       --- to true.
        ldoagain = .false.

c       --- Loop over species
        do is=is1,is2
          if (axis == 2) then
            indts = ndtstorho(pgroup%ndts(is-1))
            zgrid = zgridndts(indts)
          else
            zgrid = 0.
          endif

c         --- initialize counter and save particle's start and number
          ip = pgroup%ins(is)
          ins_init(is) = pgroup%ins(is)
          nps_init(is) = pgroup%nps(is)

c         --- Loop over particles and pick out the ones which have ventured
c         --- into the regions of neighboring processors.
          do while (ip < pgroup%ins(is) + pgroup%nps(is))

            if (lrz) then
              ztemp = sqrt(pgroup%xp(ip)**2 + pgroup%yp(ip)**2)
            else if (labs) then
              ztemp = abs(zz(ip)-zgrid)
            else
              ztemp = zz(ip)-zgrid
            endif

c           --- If particle exited to the left, switch with lowest particle
c           --- in the arrays.
            if ((ztemp <= zpminlocal .and. pboundleft == absorb) .or.
     &           ztemp < zpminlocal) then
              idest = pgroup%ins(is)
              pgroup%ins(is) = pgroup%ins(is) + 1
              pgroup%nps(is) = pgroup%nps(is) - 1

              temp = pgroup%xp(ip); pgroup%xp(ip) = pgroup%xp(idest); pgroup%xp(idest) = temp
              temp = pgroup%yp(ip); pgroup%yp(ip) = pgroup%yp(idest); pgroup%yp(idest) = temp
              temp = pgroup%zp(ip); pgroup%zp(ip) = pgroup%zp(idest); pgroup%zp(idest) = temp

              temp = pgroup%uxp(ip); pgroup%uxp(ip) = pgroup%uxp(idest); pgroup%uxp(idest) = temp
              temp = pgroup%uyp(ip); pgroup%uyp(ip) = pgroup%uyp(idest); pgroup%uyp(idest) = temp
              temp = pgroup%uzp(ip); pgroup%uzp(ip) = pgroup%uzp(idest); pgroup%uzp(idest) = temp

              temp = pgroup%gaminv(ip); pgroup%gaminv(ip) = pgroup%gaminv(idest); pgroup%gaminv(idest) = temp

              temp = pgroup%ex(ip); pgroup%ex(ip) = pgroup%ex(idest); pgroup%ex(idest) = temp
              temp = pgroup%ey(ip); pgroup%ey(ip) = pgroup%ey(idest); pgroup%ey(idest) = temp
              temp = pgroup%ez(ip); pgroup%ez(ip) = pgroup%ez(idest); pgroup%ez(idest) = temp

              temp = pgroup%bx(ip); pgroup%bx(ip) = pgroup%bx(idest); pgroup%bx(idest) = temp
              temp = pgroup%by(ip); pgroup%by(ip) = pgroup%by(idest); pgroup%by(idest) = temp
              temp = pgroup%bz(ip); pgroup%bz(ip) = pgroup%bz(idest); pgroup%bz(idest) = temp

              if (pgroup%npid > 0) then
                pidtemp = pgroup%pid(ip,:)
                pgroup%pid(ip,:) = pgroup%pid(idest,:)
                pgroup%pid(idest,:) = pidtemp
              endif

c           --- If particle exited to the right, switch with highest particle
c           --- in the arrays.
            elseif (ztemp >= zpmaxlocal) then
              idest = pgroup%ins(is) + pgroup%nps(is) - 1
              pgroup%nps(is) = pgroup%nps(is) - 1

              temp = pgroup%xp(ip); pgroup%xp(ip) = pgroup%xp(idest); pgroup%xp(idest) = temp
              temp = pgroup%yp(ip); pgroup%yp(ip) = pgroup%yp(idest); pgroup%yp(idest) = temp
              temp = pgroup%zp(ip); pgroup%zp(ip) = pgroup%zp(idest); pgroup%zp(idest) = temp

              temp = pgroup%uxp(ip); pgroup%uxp(ip) = pgroup%uxp(idest); pgroup%uxp(idest) = temp
              temp = pgroup%uyp(ip); pgroup%uyp(ip) = pgroup%uyp(idest); pgroup%uyp(idest) = temp
              temp = pgroup%uzp(ip); pgroup%uzp(ip) = pgroup%uzp(idest); pgroup%uzp(idest) = temp

              temp = pgroup%gaminv(ip); pgroup%gaminv(ip) = pgroup%gaminv(idest); pgroup%gaminv(idest) = temp

              temp = pgroup%ex(ip); pgroup%ex(ip) = pgroup%ex(idest); pgroup%ex(idest) = temp
              temp = pgroup%ey(ip); pgroup%ey(ip) = pgroup%ey(idest); pgroup%ey(idest) = temp
              temp = pgroup%ez(ip); pgroup%ez(ip) = pgroup%ez(idest); pgroup%ez(idest) = temp

              temp = pgroup%bx(ip); pgroup%bx(ip) = pgroup%bx(idest); pgroup%bx(idest) = temp
              temp = pgroup%by(ip); pgroup%by(ip) = pgroup%by(idest); pgroup%by(idest) = temp
              temp = pgroup%bz(ip); pgroup%bz(ip) = pgroup%bz(idest); pgroup%bz(idest) = temp

              if (pgroup%npid > 0) then
                pidtemp = pgroup%pid(ip,:)
                pgroup%pid(ip,:) = pgroup%pid(idest,:)
                pgroup%pid(idest,:) = pidtemp
              endif
              ip = ip - 1
            endif

c           --- advance ip
            ip = ip + 1
          enddo

c         --- Apply appropriate boundary conditions on particles at the axial
c         --- ends of the full mesh.

          i1 = ins_init(is)
          i2 = pgroup%ins(is) - 1
          if (isleftmostproc .and. i2 >= i1) then

            if (pboundleft==periodic) then
c             --- The position is adjusted, and the data will still be passed
c             --- to the appropriate neighbor.
              syslen = zmmax - zmmin
              do ip=i1,i2
                zz(ip) = zz(ip) + syslen
              enddo

            elseif (pboundleft==absorb) then
c             --- The gaminv is flagged as a lost particle. The ins and nps
c             --- counts are reset to include the lost particles so that they
c             --- will be handled properly.
              if (lcountaslost) then
                pgroup%gaminv(i1:i2) = 0.
              else
                pgroup%gaminv(i1:i2) = -1.
              endif
              pgroup%nps(is) = pgroup%nps(is) + (pgroup%ins(is)-ins_init(is))
              pgroup%ins(is) = ins_init(is)

            elseif (pboundleft==reflect) then
c             --- The position is adjusted, but the data is not sent.
c             --- The ins and nps are reset to include them.
c             --- Note that the lrz case is not needed, since r will never
c             --- be less than zero.
              do ip=i1,i2
                zz(ip) = zmmin + (zmmin - zz(ip))
                uz(ip) = -uz(ip)
              enddo
              pgroup%nps(is) = pgroup%nps(is) + (pgroup%ins(is) - ins_init(is))
              pgroup%ins(is) = ins_init(is)
            endif

          endif

          i1 = pgroup%ins(is) + pgroup%nps(is)
          i2 = ins_init(is) + nps_init(is) - 1
          if (isrightmostproc .and. i2 >= i1) then

            if (pboundright==periodic) then
c             --- The position is adjusted, and the data will still be passed
c             --- to the appropriate neighbor.
              syslen = zmmax - zmmin
              do ip=i1,i2
                zz(ip) = zz(ip) - syslen
              enddo

            elseif (pboundright==absorb) then
c             --- The gaminv is flagged as a lost particle. The ins and nps
c             --- counts are reset to include the lost particles so that they
c             --- will be handled properly. Check if there are any flagged
c             --- particles.
              if (lcountaslost) then
                pgroup%gaminv(i1:i2) = 0.
              else
                pgroup%gaminv(i1:i2) = -1.
              endif
              pgroup%nps(is) = ins_init(is) + nps_init(is) - pgroup%ins(is)

            elseif (pboundright==reflect) then
c             --- The position is adjusted, but the data is not sent.
c             --- The ins and nps are reset to include them.
              if (lrz) then
                do ip=i1,i2
                  rr = sqrt(pgroup%xp(ip)**2 + pgroup%yp(ip)**2)
                  if (rr > zmmax) then
                    pgroup%xp(ip) = (2.*zmmax - rr)*(pgroup%xp(ip)/rr)
                    pgroup%yp(ip) = (2.*zmmax - rr)*(pgroup%yp(ip)/rr)
                  else ! rr == zmmax
c                   --- Shift the particle in bounds by a number that is small
c                   --- compared to the dimensions of the grid.
                    rrin = zmmax - (zmmax - zmmin)*1.e-12
                    pgroup%xp(ip) = rrin*(pgroup%xp(ip)/rr)
                    pgroup%yp(ip) = rrin*(pgroup%yp(ip)/rr)
                  endif
                  pgroup%uxp(ip) = -pgroup%uxp(ip)
                  pgroup%uyp(ip) = -pgroup%uyp(ip)
                enddo
              else
                do ip=i1,i2
                  if (zz(ip) > zmmax) then
                    zz(ip) = zmmax - (zz(ip) - zmmax)
                  else ! zz == zmmax
c                   --- Shift the particle in bounds by a number that is small
c                   --- compared to the dimensions of the grid.
                    zz(ip) = zmmax - (zmmax - zmmin)*1.e-12
                  endif
                  uz(ip) = -uz(ip)
                enddo
              endif
              pgroup%nps(is) = ins_init(is) + nps_init(is) - pgroup%ins(is)

            endif
          endif

c       --- End of loop over species
        enddo

#ifdef VAMPIR
c     call VTEND(10,IERR)
c     call VTBEGIN(20,IERR)
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

c       --- Number of quantities to be exchanged. Normally it is thirteen,
c       --- the positions, velocities, gaminv, and the fields
c       --- but when pgroup%npid>0, then it is more since pid is also exchanged.
        nquant = 13
        if (pgroup%npid > 0) nquant = nquant + pgroup%npid

c       --- Copy particles being sent to the left to buffer1.
        if (pboundleft==periodic .or. .not. isleftmostproc) then
          nsendleft(is1:is2) = pgroup%ins(is1:is2) - ins_init(is1:is2)
          nsendleftsum = sum(nsendleft(is1:is2))
          if (nquant*nsendleftsum > b1size) then
            b1size = nquant*nsendleftsum
            call gchange("Databuffers",0)
          endif
          if (nsendleftsum > 0) then
            ii = 1
            do is=is1,is2
              ip = ins_init(is)
              nn = nsendleft(is)
              if (nn > 0) then
                buffer1(ii      :ii+   nn-1) = pgroup%xp(ip:ip+nn-1)
                buffer1(ii+   nn:ii+ 2*nn-1) = pgroup%yp(ip:ip+nn-1)
                buffer1(ii+ 2*nn:ii+ 3*nn-1) = pgroup%zp(ip:ip+nn-1)
                buffer1(ii+ 3*nn:ii+ 4*nn-1) = pgroup%uxp(ip:ip+nn-1)
                buffer1(ii+ 4*nn:ii+ 5*nn-1) = pgroup%uyp(ip:ip+nn-1)
                buffer1(ii+ 5*nn:ii+ 6*nn-1) = pgroup%uzp(ip:ip+nn-1)
                buffer1(ii+ 6*nn:ii+ 7*nn-1) = pgroup%gaminv(ip:ip+nn-1)
                buffer1(ii+ 7*nn:ii+ 8*nn-1) = pgroup%ex(ip:ip+nn-1)
                buffer1(ii+ 8*nn:ii+ 9*nn-1) = pgroup%ey(ip:ip+nn-1)
                buffer1(ii+ 9*nn:ii+10*nn-1) = pgroup%ez(ip:ip+nn-1)
                buffer1(ii+10*nn:ii+11*nn-1) = pgroup%bx(ip:ip+nn-1)
                buffer1(ii+11*nn:ii+12*nn-1) = pgroup%by(ip:ip+nn-1)
                buffer1(ii+12*nn:ii+13*nn-1) = pgroup%bz(ip:ip+nn-1)
                do ipid=1,pgroup%npid
                  buffer1(ii+(12+ipid)*nn:ii+(13+ipid)*nn-1) = pgroup%pid(ip:ip+nn-1,ipid)
                enddo
                ii = ii + nn*nquant
              endif
            enddo
          endif
        endif

c       --- Copy particles being sent to the right to buffer2.
        if (pboundright==periodic .or. .not. isrightmostproc) then
          nsendright(is1:is2) = ins_init(is1:is2) + nps_init(is1:is2)
     &                          - pgroup%ins(is1:is2) - pgroup%nps(is1:is2)
          nsendrightsum = sum(nsendright(is1:is2))
          if (nquant*nsendrightsum > b2size) then
            b2size = nquant*nsendrightsum
            call gchange("Databuffers",0)
          endif
          if (nsendrightsum > 0) then
            ii = 1
            do is=is1,is2
              ip = pgroup%ins(is) + pgroup%nps(is)
              nn = nsendright(is)
              if (nn > 0) then
                buffer2(ii      :ii+   nn-1) = pgroup%xp(ip:ip+nn-1)
                buffer2(ii+   nn:ii+ 2*nn-1) = pgroup%yp(ip:ip+nn-1)
                buffer2(ii+ 2*nn:ii+ 3*nn-1) = pgroup%zp(ip:ip+nn-1)
                buffer2(ii+ 3*nn:ii+ 4*nn-1) = pgroup%uxp(ip:ip+nn-1)
                buffer2(ii+ 4*nn:ii+ 5*nn-1) = pgroup%uyp(ip:ip+nn-1)
                buffer2(ii+ 5*nn:ii+ 6*nn-1) = pgroup%uzp(ip:ip+nn-1)
                buffer2(ii+ 6*nn:ii+ 7*nn-1) = pgroup%gaminv(ip:ip+nn-1)
                buffer2(ii+ 7*nn:ii+ 8*nn-1) = pgroup%ex(ip:ip+nn-1)
                buffer2(ii+ 8*nn:ii+ 9*nn-1) = pgroup%ey(ip:ip+nn-1)
                buffer2(ii+ 9*nn:ii+10*nn-1) = pgroup%ez(ip:ip+nn-1)
                buffer2(ii+10*nn:ii+11*nn-1) = pgroup%bx(ip:ip+nn-1)
                buffer2(ii+11*nn:ii+12*nn-1) = pgroup%by(ip:ip+nn-1)
                buffer2(ii+12*nn:ii+13*nn-1) = pgroup%bz(ip:ip+nn-1)
                do ipid=1,pgroup%npid
                  buffer2(ii+(12+ipid)*nn:ii+(13+ipid)*nn-1) = pgroup%pid(ip:ip+nn-1,ipid)
                enddo
                ii = ii + nn*nquant
              endif
            enddo
          endif
        endif

#ifdef VAMPIR
c     call VTEND(20,IERR)
c     call VTBEGIN(30,IERR)
#endif
c       --- Now send the data.

c       --- First send the number of particles to be sent
        call MPI_SENDRECV(nsendleft(is1:is2),int(is2-is1+1,MPIISZ),
     &                    MPI_INTEGER,left_pe,messid,
     &                    nrecvright(is1:is2),int(is2-is1+1,MPIISZ),
     &                    MPI_INTEGER,right_pe,messid,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)
        call MPI_SENDRECV(nsendright(is1:is2),int(is2-is1+1,MPIISZ),
     &                    MPI_INTEGER,right_pe,messid,
     &                    nrecvleft(is1:is2),int(is2-is1+1,MPIISZ),
     &                    MPI_INTEGER,left_pe,messid,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)

#ifdef VAMPIR
c     call VTEND(30,IERR)
c     call VTBEGIN(40,IERR)
#endif
#ifdef VAMPIR
c     call VTEND(40,IERR)
c     call VTBEGIN(50,IERR)
#endif

        nrecvleftsum = sum(nrecvleft(is1:is2))
        nrecvrightsum = sum(nrecvright(is1:is2))
        if (nquant*nrecvleftsum > b3size) then
          b3size = nquant*nrecvleftsum
          call gchange("Databuffers",0)
        endif
        if (nquant*nrecvrightsum > b4size) then
          b4size = nquant*nrecvrightsum
          call gchange("Databuffers",0)
        endif

        call MPI_SENDRECV(buffer2,nquant*nsendrightsum,MPI_DOUBLE_PRECISION,
     &                    right_pe,messid,
     &                    buffer3,nquant*nrecvleftsum,MPI_DOUBLE_PRECISION,
     &                    left_pe,messid,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)
        call MPI_SENDRECV(buffer1,nquant*nsendleftsum,MPI_DOUBLE_PRECISION,
     &                    left_pe,messid,
     &                    buffer4,nquant*nrecvrightsum,MPI_DOUBLE_PRECISION,
     &                    right_pe,messid,
     &                    MPI_COMM_WORLD,mpistatus,mpierror)

#ifdef VAMPIR
c     call VTEND(50,IERR)
c     call VTBEGIN(60,IERR)
#endif

        il = 1
        ir = 1
        do is=is1,is2

c         --- Make sure that there is space for the incoming data.
c         --- Note that chckpart must be called in this loop - if it is
c         --- called for all species at once, it could double count the
c         --- available space between species. If for example there where
c         --- 10 empty elements between species and each needed 8 elements,
c         --- chckpart wouldn't shift anything since it would see the 10
c         --- would be enough for the 8 for each species, not knowing that
c         --- 16 spaces is really needed.
          call chckpart(pgroup,int(is,ISZ),
     &                  int(nrecvleft(is),ISZ),int(nrecvright(is),ISZ))

c         --- Reset pointer to correct data coordinates.
c         --- This is needed in case the arrays are reallocated in chckpart.
          if (axis == 0) then
            zz => pgroup%xp
          else if (axis == 1) then
            zz => pgroup%yp
          else if (axis == 2) then
            zz => pgroup%zp
          endif

          if (nrecvleft(is) > 0) then
            nn = nrecvleft(is)
            if (nn > 0) then
              ip = pgroup%ins(is) - nn
              pgroup%xp(ip:ip+nn-1)     = buffer3(il      :il+   nn-1)
              pgroup%yp(ip:ip+nn-1)     = buffer3(il+   nn:il+ 2*nn-1)
              pgroup%zp(ip:ip+nn-1)     = buffer3(il+ 2*nn:il+ 3*nn-1)
              pgroup%uxp(ip:ip+nn-1)    = buffer3(il+ 3*nn:il+ 4*nn-1)
              pgroup%uyp(ip:ip+nn-1)    = buffer3(il+ 4*nn:il+ 5*nn-1)
              pgroup%uzp(ip:ip+nn-1)    = buffer3(il+ 5*nn:il+ 6*nn-1)
              pgroup%gaminv(ip:ip+nn-1) = buffer3(il+ 6*nn:il+ 7*nn-1)
              pgroup%ex(ip:ip+nn-1)     = buffer3(il+ 7*nn:il+ 8*nn-1)
              pgroup%ey(ip:ip+nn-1)     = buffer3(il+ 8*nn:il+ 9*nn-1)
              pgroup%ez(ip:ip+nn-1)     = buffer3(il+ 9*nn:il+10*nn-1)
              pgroup%bx(ip:ip+nn-1)     = buffer3(il+10*nn:il+11*nn-1)
              pgroup%by(ip:ip+nn-1)     = buffer3(il+11*nn:il+12*nn-1)
              pgroup%bz(ip:ip+nn-1)     = buffer3(il+12*nn:il+13*nn-1)
              do ipid=1,pgroup%npid
                pgroup%pid(ip:ip+nn-1,ipid) = buffer3(il+(12+ipid)*nn:il+(13+ipid)*nn-1)
              enddo
              il = il + nquant*nn
              pgroup%ins(is) = pgroup%ins(is) - nn
              pgroup%nps(is) = pgroup%nps(is) + nn
              if (lrz) then
                ztemp = maxval(sqrt(pgroup%xp(ip:ip+nn-1)**2 +
     &                              pgroup%yp(ip:ip+nn-1)**2))
              elseif (labs) then
                ztemp = maxval(abs(zz(ip:ip+nn-1)))-zgrid
              else
                ztemp = maxval(zz(ip:ip+nn-1))-zgrid
              endif
              if (ztemp > zpmaxlocal) then
                ldoagain = .true.
              endif
            endif
          endif

          if (nrecvright(is) > 0) then
            nn = nrecvright(is)
            if (nn > 0) then
              ip = pgroup%ins(is) + pgroup%nps(is)
              pgroup%xp(ip:ip+nn-1)     = buffer4(ir      :ir+ 1*nn-1)
              pgroup%yp(ip:ip+nn-1)     = buffer4(ir+   nn:ir+ 2*nn-1)
              pgroup%zp(ip:ip+nn-1)     = buffer4(ir+ 2*nn:ir+ 3*nn-1)
              pgroup%uxp(ip:ip+nn-1)    = buffer4(ir+ 3*nn:ir+ 4*nn-1)
              pgroup%uyp(ip:ip+nn-1)    = buffer4(ir+ 4*nn:ir+ 5*nn-1)
              pgroup%uzp(ip:ip+nn-1)    = buffer4(ir+ 5*nn:ir+ 6*nn-1)
              pgroup%gaminv(ip:ip+nn-1) = buffer4(ir+ 6*nn:ir+ 7*nn-1)
              pgroup%ex(ip:ip+nn-1)     = buffer4(ir+ 7*nn:ir+ 8*nn-1)
              pgroup%ey(ip:ip+nn-1)     = buffer4(ir+ 8*nn:ir+ 9*nn-1)
              pgroup%ez(ip:ip+nn-1)     = buffer4(ir+ 9*nn:ir+10*nn-1)
              pgroup%bx(ip:ip+nn-1)     = buffer4(ir+10*nn:ir+11*nn-1)
              pgroup%by(ip:ip+nn-1)     = buffer4(ir+11*nn:ir+12*nn-1)
              pgroup%bz(ip:ip+nn-1)     = buffer4(ir+12*nn:ir+13*nn-1)
              do ipid=1,pgroup%npid
                pgroup%pid(ip:ip+nn-1,ipid) = buffer4(ir+(12+ipid)*nn:ir+(13+ipid)*nn-1)
              enddo
              ir = ir + nquant*nn
              pgroup%nps(is) = pgroup%nps(is) + nn
              if (lrz) then
                ztemp = maxval(sqrt(pgroup%xp(ip:ip+nn-1)**2 +
     &                              pgroup%yp(ip:ip+nn-1)**2))
              elseif (labs) then
                ztemp = minval(abs(zz(ip:ip+nn-1)))-zgrid
              else
                ztemp = minval(zz(ip:ip+nn-1))-zgrid
              endif
              if (ztemp < zpminlocal) then
                ldoagain = .true.
              endif
            endif
          endif
        enddo

#ifdef VAMPIR
c     call VTEND(60,IERR)
c     call VTBEGIN(70,IERR)
#endif
        call parallellor(ldoagain)
#ifdef VAMPIR
c     call VTEND(70,IERR)
#endif

      enddo

c What follows (but is commented out) is the preferred version, but
c it causes a core dump for some unknown reason.
c Now it needs to be adjusted to send all species at once.
c       allocate(blocklengths(nquant),displacements(nquant))
c       --- Create new type for receiving particles from left
c       do ii=1,nquant
c         blocklengths(ii) = nrecvleft
c       enddo
c       ip = pgroup%ins(is) - nrecvleft
c       call MPI_ADDRESS(pgroup%xp(ip),displacements(1),mpierror)
c       call MPI_ADDRESS(pgroup%yp(ip),displacements(2),mpierror)
c       call MPI_ADDRESS(pgroup%zp(ip),displacements(3),mpierror)
c       call MPI_ADDRESS(pgroup%uxp(ip),displacements(4),mpierror)
c       call MPI_ADDRESS(pgroup%uyp(ip),displacements(5),mpierror)
c       call MPI_ADDRESS(pgroup%uzp(ip),displacements(6),mpierror)
c       call MPI_ADDRESS(pgroup%gaminv(ip),displacements(7),mpierror)
c       call MPI_ADDRESS(pgroup%ex(ip),displacements(8),mpierror)
c       call MPI_ADDRESS(pgroup%ey(ip),displacements(9),mpierror)
c       call MPI_ADDRESS(pgroup%ez(ip),displacements(10),mpierror)
c       call MPI_ADDRESS(pgroup%bx(ip),displacements(11),mpierror)
c       call MPI_ADDRESS(pgroup%by(ip),displacements(12),mpierror)
c       call MPI_ADDRESS(pgroup%bz(ip),displacements(13),mpierror)
c       do ipid=1,pgroup%npid
c         call MPI_ADDRESS(pgroup%pid(ip,ipid),displacements(nquant+ipid),mpierror)
c       enddo
c       call MPI_TYPE_HINDEXED(nquant,blocklengths,displacements,
c    &                         MPI_DOUBLE_PRECISION,mpinewtype,mpierror)
c       call MPI_TYPE_COMMIT(mpinewtype,mpierror)

c       --- Pass particle data to the right
c       call MPI_SENDRECV(buffer2,nquant*nsendright,MPI_DOUBLE_PRECISION,right_pe,0,
c    &                    MPI_BOTTOM,1,mpinewtype,left_pe,0,
c    &                    MPI_COMM_WORLD,mpistatus,mpierror)
c       pgroup%ins(is) = pgroup%ins(is) - nrecvleft
c       pgroup%nps(is) = pgroup%nps(is) + nrecvleft
c       call MPI_TYPE_FREE(mpinewtype,mpierror)

c       --- Create new type for receiving particles from right
c       do ii=1,nquant
c         blocklengths(ii) = nrecvright
c       enddo
c       ip = pgroup%ins(is) + pgroup%nps(is)
c       call MPI_ADDRESS(pgroup%xp(ip),displacements(1),mpierror)
c       call MPI_ADDRESS(pgroup%yp(ip),displacements(2),mpierror)
c       call MPI_ADDRESS(pgroup%zp(ip),displacements(3),mpierror)
c       call MPI_ADDRESS(pgroup%uxp(ip),displacements(4),mpierror)
c       call MPI_ADDRESS(pgroup%uyp(ip),displacements(5),mpierror)
c       call MPI_ADDRESS(pgroup%uzp(ip),displacements(6),mpierror)
c       call MPI_ADDRESS(pgroup%gaminv(ip),displacements(7),mpierror)
c       call MPI_ADDRESS(pgroup%ex(ip),displacements(8),mpierror)
c       call MPI_ADDRESS(pgroup%ey(ip),displacements(9),mpierror)
c       call MPI_ADDRESS(pgroup%ez(ip),displacements(10),mpierror)
c       call MPI_ADDRESS(pgroup%bx(ip),displacements(11),mpierror)
c       call MPI_ADDRESS(pgroup%by(ip),displacements(12),mpierror)
c       call MPI_ADDRESS(pgroup%bz(ip),displacements(13),mpierror)
c       do ipid=1,pgroup%npid
c         call MPI_ADDRESS(pgroup%pid(ip,ipid),displacements(nquant+ipid),mpierror)
c       enddo
c       call MPI_TYPE_HINDEXED(nquant,blocklengths,displacements,
c    &                         MPI_DOUBLE_PRECISION,mpinewtype,mpierror)
c       call MPI_TYPE_COMMIT(mpinewtype,mpierror)

c       --- Pass particle data to the left
c       call MPI_SENDRECV(buffer1,nquant*nsendleft,MPI_DOUBLE_PRECISION,left_pe,0,
c    &                    MPI_BOTTOM,1,mpinewtype,right_pe,0,
c    &                    MPI_COMM_WORLD,mpistatus,mpierror)
c       pgroup%nps(is) = pgroup%nps(is) + nrecvright
c       call MPI_TYPE_FREE(mpinewtype,mpierror)
c       deallocate(blocklengths,displacements)

!$OMP MASTER
      if (ltoptimesubs) timezpartbnd_slave = timezpartbnd_slave + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      subroutine reorgparticles_parallel(pgroup,l4symtry,l2symtry,lrz)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      use Picglb
      use Parallel
      type(ParticleGroup):: pgroup
      logical(ISZ):: l4symtry,l2symtry,lrz

c Do all to all communication to pass particles to where they belong. This
c sorts through all of the particles first, then does all of the communication
c at once.

      integer(MPIISZ),allocatable:: xprocs(:),yprocs(:),zprocs(:)
      integer(MPIISZ):: nxp,nyp,nzp
      real(kind=8):: dxp,dyp,dzp
      integer(MPIISZ):: ix,iy,iz,ip
      real(kind=8):: xx,yy,zz
      integer(MPIISZ):: ipx1,ipx2,ipy1,ipy2,ipz1,ipz2,ipx,ipy,ipz
      real(kind=8),allocatable:: ptemp(:,:),precv(:,:)
      integer(MPIISZ),allocatable:: pprocs(:)
      integer(MPIISZ):: pins(0:pgroup%ns-1,0:nprocs-1)
      integer(MPIISZ):: pnps(0:pgroup%ns-1,0:nprocs-1)
      integer(MPIISZ):: piii(0:pgroup%ns-1,0:nprocs-1)
      integer(MPIISZ):: inrecv(0:pgroup%ns-1,0:nprocs-1)
      integer(MPIISZ):: nprecv(0:pgroup%ns-1,0:nprocs-1)
      integer(MPIISZ):: pnp,pnpmax,ncoords,ii,js,pip,ip1,ip2,iitemp,ii1,ii2
      integer(MPIISZ):: mpierror
      integer(MPIISZ):: sendcounts(0:nprocs-1),sdispls(0:nprocs-1)
      integer(MPIISZ):: recvcounts(0:nprocs-1),rdispls(0:nprocs-1)
      integer(MPIISZ):: allocerror
      logical(ISZ):: skipx,skipy,skipz
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

c     --- First, create arrays that specify which processors owns which
c     --- location. This is used to make the sorting of particles to
c     --- processors efficient. Using this array avoids having to loop
c     --- over every processor to find which one the particle is in.
      if (fsdecomp%nxglobal > 0) then
        dxp = minval(ppdecomp%xmax-ppdecomp%xmin)/10.
        nxp = min(100000,int((xpmax-xpmin)/dxp))
        dxp = (xpmax - xpmin)/nxp
        skipx = .false.
      else
        nxp = 1
        dxp = 1.
        skipx = .true.
      endif
      if (fsdecomp%nyglobal > 0) then
        dyp = minval(ppdecomp%ymax-ppdecomp%ymin)/10.
        nyp = min(100000,int((ypmax-ypmin)/dyp))
        dyp = (ypmax - ypmin)/nyp
        skipy = .false.
      else
c       --- Set nyp to 1, so that ipy2 is OK.
        nyp = 1
        dyp = 1.
        skipy = .true.
      endif
      if (fsdecomp%nzglobal > 0) then
        dzp = minval(ppdecomp%zmax-ppdecomp%zmin)/10.
        nzp = min(100000,int((zpmax-zpmin)/dzp))
        dzp = (zpmax - zpmin)/nzp
        skipz = .false.
      else
        nzp = 1
        dzp = 1.
        skipz = .true.
      endif

      allocate(xprocs(0:nxp),yprocs(0:nyp),zprocs(0:nzp),stat=allocerror)
      if (allocerror /= 0) then
        print*,"reorgparticles_parallel: allocation error ",allocerror,
     &         ": could not allocate x,y,zprocs to shape ",nxp,nyp,nzp
        call kaboom("reorgparticles_parallel: allocation error")
        return
      endif

      ip = 0
      do ix=0,nxp
        xx = xpmin + ix*dxp
        do while (xx > ppdecomp%xmax(ip) .and. ip < nxprocs-1)
          ip = ip + 1
        enddo
        xprocs(ix) = ip
      enddo

      ip = 0
      do iy=0,nyp
        yy = ypmin + iy*dyp
        do while (yy > ppdecomp%ymax(ip) .and. ip < nyprocs-1)
          ip = ip + 1
        enddo
        yprocs(iy) = ip
      enddo

      ip = 0
      do iz=0,nzp
        zz = zpmin + iz*dzp
        do while (zz > ppdecomp%zmax(ip) .and. ip < nzprocs-1)
          ip = ip + 1
        enddo
        zprocs(iz) = ip
      enddo

c     --- Next, create temporary arrays to hold the sorted particles.
      pnp = sum(pgroup%nps)
      ncoords = 13 + pgroup%npid
      pnpmax = max(int(1,MPIISZ),pnp)
      allocate(ptemp(0:ncoords-1,pnpmax),pprocs(pnpmax),stat=allocerror)
      if (allocerror /= 0) then
        print*,"reorgparticles_parallel: allocation error ",allocerror,
     &         ": could not allocate ptemp and pprocs to shape ",ncoords,pnpmax
        call kaboom("reorgparticles_parallel: allocation error")
        return
      endif

c     --- Loop through the particles, finding to which processor each particle
c     --- belongs.
      pip = 0
      pnps = 0
      pprocs = 1000000
      do js=0,pgroup%ns-1
        PLOOP: do ii=pgroup%ins(js+1),pgroup%ins(js+1)+pgroup%nps(js+1)-1
          xx = pgroup%xp(ii)
          yy = pgroup%yp(ii)
          zz = pgroup%zp(ii) - zbeam
          if (l4symtry) xx = abs(xx)
          if (l4symtry .or. l2symtry) yy = abs(yy)
          if (lrz) then
            xx = sqrt(xx**2 + yy**2)
            yy = ypmin
          else if (skipy) then
            yy = ypmin
          endif
          ix = int((xx - xpmin)/dxp)
          iy = int((yy - ypmin)/dyp)
          iz = int((zz - zpmin)/dzp)
          ipx1 = xprocs(ix)
          ipx2 = xprocs(ix+1)
          ipy1 = yprocs(iy)
          ipy2 = yprocs(iy+1)
          ipz1 = zprocs(iz)
          ipz2 = zprocs(iz+1)
          pip = pip + 1
          do ipz=ipz1,ipz2
            do ipy=ipy1,ipy2
              do ipx=ipx1,ipx2
                if ((ppdecomp%xmin(ipx) <= xx .and. xx < ppdecomp%xmax(ipx)
     &               .or. skipx) .and.
     &              (ppdecomp%ymin(ipy) <= yy .and. yy < ppdecomp%ymax(ipy)
     &               .or. skipy) .and.
     &              (ppdecomp%zmin(ipz) <= zz .and. zz < ppdecomp%zmax(ipz)
     &               .or. skipz)) then
c                 ip = convertindextoproc(ipx,ipy,ipz,
c    &                                    ppdecomp%nxprocs,
c    &                                    ppdecomp%nyprocs,
c    &                                    ppdecomp%nzprocs)
c                 --- This is the same expression in convertindextoproc, but
c                 --- is written out here for optimization.
                  ip = ipx + ipy*ppdecomp%nxprocs +
     &                       ipz*ppdecomp%nxprocs*ppdecomp%nyprocs
                  pprocs(pip) = ip
                  pnps(js,ip) = pnps(js,ip) + 1
                  cycle PLOOP
                endif
              enddo
            enddo
          enddo
c         print*,"BADx ",js,ii,xx,ppdecomp%xmin(ipx1),ppdecomp%xmax(ipx2),ipx1,ipx2
c         print*,"BADy ",js,ii,yy,ppdecomp%ymin(ipy1),ppdecomp%ymax(ipy2),ipy1,ipy2
c         print*,"BADz ",js,ii,zz,ppdecomp%zmin(ipz1),ppdecomp%zmax(ipz2),ipz1,ipz2
        enddo PLOOP
      enddo

c     --- Calculate starting index of particle groups for each species and
c     --- processor
      pins(0,0) = 1
      do ip=0,nprocs-1
        do js=0,pgroup%ns-1
          if (js == 0) then
            if (ip == 0) then
              pins(js,ip) = 1
            else
              pins(js,ip) = pins(ns-1,ip-1) + pnps(ns-1,ip-1)
            endif
          else
            pins(js,ip) = pins(js-1,ip) + pnps(js-1,ip)
          endif
        enddo
      enddo

c     --- Copy sorted particle data into the temporary array
      pip = 0
      piii = pins
      do js=0,pgroup%ns-1
        do ii=pgroup%ins(js+1),pgroup%ins(js+1)+pgroup%nps(js+1)-1
          pip = pip + 1
          ip = pprocs(pip)
c         if (ip == 1000000) print*,"PIP ",js,ii,pgroup%xp(ii),pgroup%yp(ii),pgroup%zp(ii)
          iitemp = piii(js,ip)
          piii(js,ip) = piii(js,ip) + 1
          ptemp( 0,iitemp) = pgroup%xp(ii)
          ptemp( 1,iitemp) = pgroup%yp(ii)
          ptemp( 2,iitemp) = pgroup%zp(ii)
          ptemp( 3,iitemp) = pgroup%uxp(ii)
          ptemp( 4,iitemp) = pgroup%uyp(ii)
          ptemp( 5,iitemp) = pgroup%uzp(ii)
          ptemp( 6,iitemp) = pgroup%gaminv(ii)
          ptemp( 7,iitemp) = pgroup%ex(ii)
          ptemp( 8,iitemp) = pgroup%ey(ii)
          ptemp( 9,iitemp) = pgroup%ez(ii)
          ptemp(10,iitemp) = pgroup%bx(ii)
          ptemp(11,iitemp) = pgroup%by(ii)
          ptemp(12,iitemp) = pgroup%bz(ii)
          if (pgroup%npid > 0) then
            ptemp(13:,iitemp) = pgroup%pid(ii,:)
          endif
        enddo
      enddo

c     --- Calculate pins for data set sent to each processor relative to 0
c     do ip=0,nprocs-1
c       pins(:,ip) = pins(:,ip) - pins(0,ip)
c     enddo

c     --- Exchange the number of particles which are to be sent to each
c     --- processor
c     call MPI_ALLTOALL(pins,int(pgroup%ns,MPIISZ),MPI_INTEGER,
c    &                  inrecv,int(pgroup%ns,MPIISZ),MPI_INTEGER,
c    &                  MPI_COMM_WORLD,mpierror)
      call MPI_ALLTOALL(pnps,int(pgroup%ns,MPIISZ),MPI_INTEGER,
     &                  nprecv,int(pgroup%ns,MPIISZ),MPI_INTEGER,
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
      if (pnp > pgroup%npmax) pgroup%npmax = int(1.1*pnp)
      if (pnp > 0) then
        np_s = int(np_s*pgroup%npmax/pnp)
      else
        np_s = pgroup%npmax/ns
      endif
      np_s(ns) = pgroup%npmax - sum(np_s(1:ns-1))
      call alotpart(pgroup)
      
c     --- Create space for incoming data
      pnpmax = max(int(1,MPIISZ),pnp)
      allocate(precv(0:ncoords-1,pnpmax),stat=allocerror)
      if (allocerror /= 0) then
        print*,"reorgparticles_parallel: allocation error ",allocerror,
     &         ": could not allocate precv to shape ",ncoords,pnpmax
        call kaboom("reorgparticles_parallel: allocation error")
        return
      endif

c     --- Do the communication
      sendcounts = sum(pnps,1)*ncoords
      recvcounts = sum(nprecv,1)*ncoords
      sdispls(0) = 0
      rdispls(0) = 0
      do ip=1,nprocs-1
        sdispls(ip) = sdispls(ip-1) + sendcounts(ip-1)
        rdispls(ip) = rdispls(ip-1) + recvcounts(ip-1)
      enddo
      call MPI_ALLTOALLV(ptemp,sendcounts,sdispls,MPI_DOUBLE_PRECISION,
     &                   precv,recvcounts,rdispls,MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD,mpierror)

c     --- Calculate inrecv for data set sent to each processor relative
c     --- to rdispls
c     do ip=0,nprocs-1
c       inrecv(:,ip) = inrecv(:,ip) + rdispls(ip)/ncoords + 1
c     enddo

c     --- Calculate starting index of particle groups for each species and
c     --- processor
      inrecv(0,0) = 1
      do ip=0,nprocs-1
        do js=0,pgroup%ns-1
          if (js == 0) then
            if (ip == 0) then
              inrecv(js,ip) = 1
            else
              inrecv(js,ip) = inrecv(ns-1,ip-1) + nprecv(ns-1,ip-1)
            endif
          else
            inrecv(js,ip) = inrecv(js-1,ip) + nprecv(js-1,ip)
          endif
        enddo
      enddo

c     --- Now copy the data into the particle arrays
      pgroup%nps = 0
      do js=0,pgroup%ns-1
        do ip=0,nprocs-1
          if (nprecv(js,ip) > 0) then
            ip1 = inrecv(js,ip)
            ip2 = ip1 + nprecv(js,ip) - 1
            ii1 = pgroup%ins(js+1) + pgroup%nps(js+1)
            ii2 = ii1 + nprecv(js,ip) - 1
            pgroup%nps(js+1) = pgroup%nps(js+1) + nprecv(js,ip)
            pgroup%xp(ii1:ii2)     = precv( 0,ip1:ip2)
            pgroup%yp(ii1:ii2)     = precv( 1,ip1:ip2)
            pgroup%zp(ii1:ii2)     = precv( 2,ip1:ip2)
            pgroup%uxp(ii1:ii2)    = precv( 3,ip1:ip2)
            pgroup%uyp(ii1:ii2)    = precv( 4,ip1:ip2)
            pgroup%uzp(ii1:ii2)    = precv( 5,ip1:ip2)
            pgroup%gaminv(ii1:ii2) = precv( 6,ip1:ip2)
            pgroup%ex(ii1:ii2)     = precv( 7,ip1:ip2)
            pgroup%ey(ii1:ii2)     = precv( 8,ip1:ip2)
            pgroup%ez(ii1:ii2)     = precv( 9,ip1:ip2)
            pgroup%bx(ii1:ii2)     = precv(10,ip1:ip2)
            pgroup%by(ii1:ii2)     = precv(11,ip1:ip2)
            pgroup%bz(ii1:ii2)     = precv(12,ip1:ip2)
            do ii=1,pgroup%npid
              pgroup%pid(ii1:ii2,ii) = precv(12+ii,ip1:ip2)
            enddo
          endif
        enddo
      enddo

c     --- Free work space
      deallocate(xprocs,yprocs,zprocs,ptemp,precv,pprocs)

c     --- Done!

!$OMP MASTER
      if (ltoptimesubs) timereorgparticles_parallel = timereorgparticles_parallel + wtime() - substarttime
!$OMP END MASTER

      return
      end
c=============================================================================
c=============================================================================
      integer(ISZ) function checkzpartbnd(pgroup)
      use ParticleGroupmodule
      use Subtimerstop
      use InPart
      use Picglb,Only: zpminlocal,zpmaxlocal
      use Parallel
      use Subcycling, Only: ndtstorho,zgridndts
      type(ParticleGroup):: pgroup
c Check if all particles are within the mesh.  Returns number of particles
c outside the mesh.
      integer(ISZ):: is,ip,nout,indts
      real(kind=8):: zgrid
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      nout = 0
      do is=1,pgroup%ns
        indts = ndtstorho(pgroup%ndts(is-1))
        zgrid = zgridndts(indts)
!$OMP PARALLEL DO REDUCTION(+:nout)
        do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
          if (pgroup%zp(ip)-zgrid < zpminlocal .or.
     &        pgroup%zp(ip)-zgrid >= zpmaxlocal) then
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
      real(kind=8),allocatable:: ztemp0(:,:)
      real(kind=8),allocatable:: ztemp(:,:,:)
      real(kind=8),allocatable:: maxes1(:,:),mines1(:,:)
      real(kind=8),allocatable:: maxes2(:,:),mines2(:,:)
      integer(MPIISZ):: nn
      include "mpif.h"
      integer(MPIISZ):: mpierror
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      allocate(ztemp0(NUMZMMNT,0:nszmmnt))
      allocate(ztemp(0:nzmmnt,NUMZMMNT,0:nszmmnt))
      allocate(maxes1(6,0:nszmmnt),mines1(6,0:nszmmnt))
      allocate(maxes2(6,0:nszmmnt),mines2(6,0:nszmmnt))

c     --- Do reduction on whole beam moments.
      ztemp0 = zmmnts0
      nn = NUMZMMNT*(1+nszmmnt)
      call MPI_ALLREDUCE(ztemp0,zmmnts0,nn,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)

c     --- Do reduction on beam z moments.
      ztemp = zmmnts
      nn = (1+nzmmnt)*NUMZMMNT*(1+nszmmnt)
c     print*,"PSM1 ",my_index,nn
      call MPI_ALLREDUCE(ztemp,zmmnts,nn,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
c     call MPI_REDUCE(ztemp,zmmnts,nn,
c    &                MPI_DOUBLE_PRECISION,MPI_SUM,int(0,MPIISZ),MPI_COMM_WORLD,mpierror)
c     call MPI_BCAST(zmmnts,nn,MPI_DOUBLE_PRECISION,int(0,MPIISZ),MPI_COMM_WORLD,mpierror)

c     print*,"PSM2 ",my_index
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
      call MPI_ALLREDUCE(maxes1,maxes2,int(6*(1+nszmmnt),MPIISZ),
     &                   MPI_DOUBLE_PRECISION,
     &                   MPI_MAX,MPI_COMM_WORLD,mpierror)
      call MPI_ALLREDUCE(mines1,mines2,int(6*(1+nszmmnt),MPIISZ),
     &                   MPI_DOUBLE_PRECISION,
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

      deallocate(ztemp0)
      deallocate(ztemp)
      deallocate(maxes1,mines1)
      deallocate(maxes2,mines2)

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

      real(kind=8),allocatable:: data2(:)
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()
      allocate(data2(ndata))

      call MPI_ALLREDUCE(data,data2,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,MPI_SUM,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2
      deallocate(data2)

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

      integer(MPIISZ),allocatable:: data2(:)
      integer(MPIISZ),allocatable:: data3(:)
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()
      allocate(data2(ndata))

      if (ISZ == 4) then
        call MPI_ALLREDUCE(data,data2,ndata,MPI_INTEGER,MPI_SUM,
     &                     MPI_COMM_WORLD,mpierror)
      else
        allocate(data3(ndata))
        data3 = data
        call MPI_ALLREDUCE(data3,data2,int(ndata,MPIISZ),MPI_INTEGER,MPI_SUM,
     &                     MPI_COMM_WORLD,mpierror)
        deallocate(data3)
      endif
      data = data2
      deallocate(data2)

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

      real(kind=8),allocatable:: data2(:)
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()
      allocate(data2(ndata))

      call MPI_ALLREDUCE(data,data2,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2
      deallocate(data2)

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

      real(kind=8),allocatable:: data2(:)
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()
      allocate(data2(ndata))

      call MPI_ALLREDUCE(data,data2,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,MPI_MIN,
     &                   MPI_COMM_WORLD,mpierror)
      data = data2
      deallocate(data2)

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

      logical(MPIISZ):: data1,data2
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      data1 = data
      call MPI_ALLREDUCE(data1,data2,int(1,MPIISZ),MPI_LOGICAL,MPI_LOR,
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

      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_BCAST(data,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,
     &               int(root,MPIISZ),MPI_COMM_WORLD,mpierror)

!$OMP MASTER
      if (ltoptimesubs) timeparallelbroadcastrealarray = timeparallelbroadcastrealarray + wtime() - substarttime
!$OMP END MASTER

      return
      end
c==========================================================================
      subroutine parallelbarrier()
      use Subtimerstop
c Calls MPI barrier
      integer(MPIISZ):: mpierror
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
      integer(MPIISZ):: mpierror
      include "mpif.h"
      real(kind=8):: substarttime,wtime
      if (ltoptimesubs) substarttime = wtime()

      call MPI_ALLREDUCE(data,dmax,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,MPI_MAX,
     &                   MPI_COMM_WORLD,mpierror)
      call MPI_ALLREDUCE(data,dmin,int(ndata,MPIISZ),MPI_DOUBLE_PRECISION,MPI_MIN,
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
