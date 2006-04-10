#include "top.h"
c=============================================================================
c@(#) File WRZ.M, version $Revision: 3.19 $, $Date: 2006/04/10 18:23:20 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package WRZ of code WARP
c  RZ electrostatic PIC code, cylindrical geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c  Debbie Callahan, LLNL, (510)423-5926
c=============================================================================
      subroutine wrzinit
 
c  Called at first reference to package (not nec. a "run" etc.).
 
      call wrzvers (STDOUT)
 
      return
      end
c=============================================================================
      subroutine wrzvers (iout)
      use WRZversion
      integer(ISZ):: iout
c  Echoes code version, etc. to output files as they're created
      call printpkgversion(iout,"Particle package WRZ",verswrz)
      return
      end
c=============================================================================
      subroutine wrzgen()
      use GlobalVars
      use Beam_acc
      use Ch_var
      use Constant
      use InGen
      use InGenrz
      use InDiag
      use InPart
      use InPartrz
      use InMeshrz
      use Picglb
      use Picglbrz
      use Fieldsrz
      use Io
      use Particles,Only: pgroup,npmax,wpid,nplive
      use OutParams
      use Z_arrays
      use Z_Moments
      use Win_Moments
      use Moments
      use Hist
      use FFTDiagrz
      use Sortrz
      use LinearBeam
 
c  Invoked by the GENERATE command, it sets up the problem
 

c  declare scratch space for getzmmnt
       
       integer(ISZ):: i,k,ia,izb,is,ip,ipmin,iwin,iz,num,nzma,ismax
       integer(ISZ):: nzmb,ieearmax,ieearmin,ir
       real(kind=8):: time1,time2,er,et,ez,frext,zm,zs,eearsprs,eearmax,eearmin
       real(kind=8):: uxpo(nparpgrp), uypo(nparpgrp), uzpo(nparpgrp)
       real(kind=8):: wtimeoff

c    Heap check on NLTSS       
c      ierr = mzparam("check",1)
c      write(STDOUT,*)"Ierr = ",ierr
c    Heap check on NERSC
c      integer(ISZ):: list(7)
c      do ilist = 1,6
c      list(ilist)=-1
c      enddo
c      list(7) = 1
c      ierr = mzparam(list,7)
 
c  Announce that we're starting up
 
      write (STDOUT,10)
   10 format(" ***  particle simulation package RZ generating  ***")

c  Estimate wall radius, needed for g-factor calc
      
      rwall = rwallfac * rmmax

c  Calculate derived quantities

      call derivqty 
 
c  Initialize the cycle counter, time, etc.
 
      if (nrestart .eq. " ") then
         it = 0
         time = 0.
      endif
      call stepid (it, time, zgrid)
 
c  Print values of input variables, other interesting things to text file
      if (warpout > -1) then 
        call edit (warpout, "runid")
        call edit (warpout, "it")
        call edit (warpout, "time")
        call edit (warpout, "InGen")
        call edit (warpout, "InDiag")
        call edit (warpout, "InPartrz")
        call edit (warpout, "InMeshrz")
      endif
 
c  Create the dynamic arrays for fields, contour plot workspace
      nmrz  = max(nr, nz)
      nrz = (nr+1) * (nz+1)
      call gallot ("Fieldsrz", 0)
 
c  Calculate mesh dimensioning quantities
      dr = (rmmax - rmmin) / nr
      dz = (zmmax - zmmin) / nz
      do i = 0, nr
         rmesh(i) = i * dr + rmmin
      enddo

c  Allocate the 1-d arrays
      nzzarr = nz
      call gallot ("Z_arrays", 0)

      do k = 0, nz
         zmesh(k) = k * dz + zmmin
      enddo

      do k = 0, nz
         zplmesh(k) = zmesh(k)
      enddo
 
c  Create the dynamic arrays for particles (set npmax to an estimated 
c  length for now, for those loading schemes that don't actually
c  use a user-set npmax directly)
 
      if (xrandom .eq. "grid") npmax = nxstripe*nystripe*nzstripe
      if (xrandom .eq. "fibonacc") npmax = nfibgrps*fibg1
      if (shell) npmax = nz
      if (linbeam) then
          ia = int(a0/dr)
          izb = int((zimax - zimin)/dz)
c          npmax = 8*ia/2*(ia+1)*izb
          npmax = (4*(ia+1)*(ia+2)+1)*izb
      endif
      call gchangeparticles()

c  Load the particles, calculate the charge density and current
 
      if (shell) then
         call ldshell
      elseif (linbeam) then
         call ldlinbm
      else
         call stptclrz
      endif
c     if (.not. shell) then
c        call stptclrz
c     else
c        call ldshell
c     endif

      call setgamma(lrelativ)
c     call zeroarry(rho,nrz)
c     call zeroarry(curr,nz+1)
      rho = 0.
      curr = 0.

      call wtimeon 
      do is=1,ns
         do ipmin = pgroup%ins(is), pgroup%ins(is) + pgroup%nps(is) - 1, nparpgrp
            ip = min(nparpgrp, pgroup%ins(is)+pgroup%nps(is)-ipmin)
            call setrhorz (rho,ip,pgroup%xp(ipmin),pgroup%zp(ipmin),zbeam,
     &                     pgroup%uzp(ipmin),pgroup%sq(is),pgroup%sw(is),depos)
            call setcurr(curr,ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),pgroup%gaminv(ipmin),
     &                   pgroup%sq(is),pgroup%sw(is),zbeam,dz,zmmin)
         enddo
      enddo
      time1 = wtimeoff()
      if (periinz) call perrhorz(rho,nr,nz)
      call fixrhorz

      ignmax = 1
      npsort = npmax
      call gallot ("Sortrz", 0)

      if (vecrho) then
c     call zeroarry(rhov,nrz)
c     call zeroarry(currv,nz+1)
      rhov = 0.
      currv = 0.
      call wtimeon 
        do is = 1,ns
          call psort(pgroup%nps(is),zbeam,dz)
          call vstrhorz(rhov,zbeam,pgroup%sq(is),pgroup%sw(is),rmesh)
          call vstcurr (currv,zbeam,pgroup%sq(is),pgroup%sw(is))
        enddo
      time2 =  wtimeoff()
      call perrhorz(rhov,nr,nz)
      endif

      write(STDOUT,*)"Scalar time",time1
      write(STDOUT,*)"Vector time",time2
      if (frcetype .eq. "analytic") call getforce
      if (frcetype .eq. "dedr") call getforce

c  Initialize the surface charge if we are including resistivity
      call setsc
 
c  Change the size of the dynamic arrays for particles to true length
 
      call gchangeparticles()

c  Create the dynamic arrays for particle moments; calculate moments
 
      zwindows(1,0) = zmmin
      zwindows(2,0) = zmmax
      nzwind = 0
      do iwin = 1, NWINDOWS
         if (zwindows(1,iwin) .ne. zwindows(2,iwin)) nzwind = nzwind + 1
      enddo
      call gallot ("Win_Moments", 0)
      call gallot ("Moments", 0)
c     call getemits(1,xp,yp,pgroup%zp,pgroup%uxp,pgroup%uyp,pgroup%uzp,pgroup%gaminv,pgroup%sm(1),pgroup%sw(1),1,nplive)
c     do is=1,ns
c        do ipmin = pgroup%ins(is), pgroup%ins(is) + pgroup%nps(is) - 1, nparpgrp
c           ip = min(nparpgrp, pgroup%ins(is)+pgroup%nps(is)-ipmin)
c           nplive = ip
c           call getemits(ip,xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
c    &                    pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
c    &                    pgroup%gaminv(ipmin),pgroup%sm(1),pgroup%sw(1),2,nplive)
c        enddo
c        call getemits(1,xp,yp,pgroup%zp,pgroup%uxp,pgroup%uyp,pgroup%uzp,pgroup%gaminv,pgroup%sm(1),pgroup%sw(1),3,nplive)
 
c     enddo
c  Create the dynamic arrays for z moments; calculate moments

      if (ifzmmnt .gt. 0) then
        nzmmnt = nz
        zmmntmax = zmmax
        zmmntmin = zmmin
        call gallot ("Z_Moments", 0)
        do k=0,nzmmnt
           zmntmesh(k) = zmmntmin + k*(zmmntmax-zmmntmin)/nzmmnt
        enddo
        ismax = maxval(pgroup%sid)+1
        call getzmmnt(1,0.,0.,0.,0.,0.,
     &                0.,0.,0.,0.,
     &                0.,0.,1,
     &                nplive,0.,0.,0.,1,1,ismax,
     &                tempmaxp,tempminp,tempzmmnts0,tempzmmnts)
        do is=1,ns
           do ipmin = pgroup%ins(is), pgroup%ins(is) + pgroup%nps(is) - 1, nparpgrp
              num = min(nparpgrp, pgroup%ins(is)+pgroup%nps(is)-ipmin)
c             --- Use t=0 u data also as "old" data for extraps
              call copyarry (pgroup%uxp(ipmin),uxpo,num)
              call copyarry (pgroup%uyp(ipmin),uypo,num)
              call copyarry (pgroup%uzp(ipmin),uzpo,num)
              if(wpid==0) then
                call getzmmnt(num,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                        pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     & pgroup%gaminv(ipmin),pgroup%sq(1),pgroup%sm(1),pgroup%sw(1),dt/2.,2,nplive,
     &                        uxpo,uypo,uzpo,is,ns)
              else
                call getzmmnt_weights(num,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                        pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     & pgroup%gaminv(ipmin),pgroup%pid(ipmin,wpid),pgroup%sq(1),pgroup%sm(1),pgroup%sw(1),dt/2.,2,nplive,
     &                        uxpo,uypo,uzpo,is,ns)
              endif
           enddo
           call getzmmnt(1,pgroup%xp,pgroup%yp,pgroup%zp,pgroup%uxp,
     &                   pgroup%uyp,pgroup%uzp,
     &                   pgroup%gaminv,pgroup%sq(1),pgroup%sm(1),pgroup%sw(1),dt/2.,
     &                   3,nplive,uxpo,uypo,uzpo,1,ns)
         enddo
      endif
 
c  Create the scratch arrays for phase space plots (permanent, for now)
c  and set limits for plots
 
c     npsplt = min(nparpgrp, npmax)
      if (xplmin .eq. 0.) xplmin = rmmin
      if (xplmax .eq. 0.) xplmax = rmmax
      if (yplmin .eq. 0.) yplmin = rmmin
      if (yplmax .eq. 0.) yplmax = rmmax
      if (zplmin .eq. 0.) zplmin = zmmin
      if (zplmax .eq. 0.) zplmax = zmmax
 
c  Create the dynamic arrays for history data; set pointer into them
 
      lenhist = min ( nt/nhist + 1 , 100) 
      call gallot ("Hist", 0)
      jhist = -1

c  Allocate arrays which are used to find the fourier transform of
c  a quantity as a diagnostic
      
      nzfft = nzzarr
      call gallot ("FFTDiagrz", 0)
 
c  Print interesting things to plot file and teletype
      call prntpara(dr,dr,dz,lprntpara)
      call prntparz(lprntpara)
 
c  Initial call to fieldsolver in order to initialize attz, kzsq, etc.
      call vprz (1)

c  Call the field solver so  we can plot phi on axis at t=0

      call copyarry (rho, phi(0,0), nrz)
      call vprz (-1)
      if (periinz) call perphirz(phi(0,-1),nr,nz)

c  Call seterz to fill E_z and E_r arrays
         call seterz (-1,phi,force,1,pgroup%xp(1),pgroup%zp(1),zbeam,
     &                rmmin,zmmin,zmmax,dr,dz,nr,nz,
     &                efetch,er,et,ez,frext,erfld,ezfld,eoffrz)
         call setexte(it,dt,dz,vbeam,nr,nz)

c  If we are looking at a mode which is initially excited, find the index
c  corresponding to that mode
      if (fftdiag) call getindx(kzsq,kzpert,indxk,nz,dr)
 
c  Initial fieldsolve, diagnostics
      call steprz ("wrzgen  ")

c  Set up Eears of z
      zm = (zimax - zimin)*(1. - straight)/2.
      emitlong = zm*2.*vthz/vbeam
      if (ifeears .eq. 1) then
c       --- set to be linear in z (eears is actually slope of Eears)
        eears = eears*(- 2*ibeam/(vbeam*zm**2)*gfactor/(2.*Pi*eps0)
     &                  - pgroup%sm(1)*vbeam**2*emitlong**2/(pgroup%sq(1)*zm**4))
      elseif (ifeears .eq. 2) then
c       --- set to initial field on axis (with the center part zero)
c       --- the '5' is just a guestimate, it should be nz dependent
        zs = (zimax - zimin)*straight/2.
        nzma = nz/2-int(zs/dz)
        nzmb = nz/2+int(zs/dz)
        do iz=0,nzma+5
          eearsofz(iz) = - ezax(iz)
        enddo
        do iz=nzmb-5,nzzarr
          eearsofz(iz) = - ezax(iz)
        enddo
c       --- Add on linear pressure term
        if (emitlong .ne. 0.) then
          eearsprs = - pgroup%sm(1)*vbeam**2*emitlong**2/(pgroup%sq(1)*zm**4)
          do iz=0,nzma
            eearsofz(iz) = eearsofz(iz) + eearsprs*(zmesh(iz) + zs)
          enddo
          do iz=nzmb,nzzarr
            eearsofz(iz) = eearsofz(iz) + eearsprs*(zmesh(iz) - zs)
          enddo
        endif
c  Make the fields outside the same as the maximum value
        eearmax = -LARGEPOS
        eearmin =  LARGEPOS
        do iz = 0,nz/2
          if (eearsofz(iz) .gt. eearmax) then
              eearmax = eearsofz(iz)
              ieearmax = iz
          endif
        enddo
        do iz = 0,ieearmax
          eearsofz(iz) = eearmax
        enddo
        do iz = nz/2,nz
          if (eearsofz(iz) .lt. eearmin) then
             eearmin = eearsofz(iz)
             ieearmin = iz
          endif
        enddo
        do iz = ieearmin,nz
          eearsofz(iz) = eearmin
        enddo
      elseif (ifeears .eq. 3) then
        do iz = 0,nz
        do ir = 0,nr
           eearsrz(ir,iz) = - ezfld(ir,iz)
        enddo
        enddo
      endif

      return
      end
c=============================================================================
      subroutine wrzexe()
      use GlobalVars
      use InGen
      use InPart
      use InGenrz
      use Picglb
      use Ctl_to_pic
      use Beam_acc
      use InMeshrz
      use Particles,Only: pgroup,nplive
      use Picglbrz
      use Fieldsrz
 
c  Takes a time step

      integer(ISZ):: is,ipmin,ip
 
c  Announce that we're running
 
      if (it.eq.0) write (STDOUT,10)
   10 format(" ***  particle simulation package RZ running  ***")
 
c  If it=0 calculate rho.  This way we can play with the initial
c  conditions in Python after we generate.
 
      if (it.eq.0) then
c        call zeroarry(rho,nrz)
         rho = 0.
         do is=1,ns
            do ipmin = pgroup%ins(is), pgroup%ins(is) + pgroup%nps(is) - 1, nparpgrp
               ip = min(nparpgrp, pgroup%ins(is)+pgroup%nps(is)-ipmin)
               call setrhorz (rho,ip,pgroup%xp(ipmin),pgroup%zp(ipmin),zbeam,
     &                        pgroup%uzp(ipmin),pgroup%sq(is),pgroup%sw(is),depos)
            enddo
         enddo
         if (periinz) call perrhorz(rho,nr,nz)
         call fixrhorz

         if (vecrho) then
c           call zeroarry(rhov,nrz)
            rhov = 0.
            do is = 1,ns
              call psort(pgroup%nps(is),zbeam,dz)
              call vstrhorz(rhov,zbeam,pgroup%sq(is),pgroup%sw(is),rmesh)
            enddo
            call perrhorz(rhov,nr,nz)
         endif
         call copyarry (rho, phi(0,0), nrz)
         call vprz (-1)
         if (periinz) call perphirz(phi(0,-1),nr,nz)
      endif
 
c  set timestep counter, time
 
      it = it + 1
      time = it * dt
      zgrid = zgrid + vbeam*dt
      call stepid (it, time, zgrid)
 
c  set logicals
 
      lfirst = .false.
      if (ncall .eq. 1) lfirst = .true.
      llast = .false.
      if (ncall .eq. maxcalls) llast = .true.
 
c  call the routine that does the actual work
 
      call steprz ("wrzexe  ")
c  Have we reached the end of the run ?
 
      if ( lfinishd ) then
         call remark("wrzexe: problem completed.")
         return
      elseif (nplive .le. 0) then
         write (STDOUT,20) nplive
   20    format(" *** WRZEXE: stopping, nplive =",i6)
         return
      else
         return
      endif
 
      end
c=============================================================================
      subroutine wrzfin()
      use InGen
      use InGenrz
      use InDiag
      use InPart
      use InPartrz
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      use Win_Moments
      use Z_arrays
      use Z_Moments
      use Hist
      use Io
      use Particles,Only: pgroup
      use Picglb
 
c  Finish up at end of RUN, or on receipt of FIN
 
 
c  print final edits
 
c  perform diagnostics (unless we just did)
 
c  Create history plots
 
c  Make a restart dump (unless we just did, or the user inhibits it)
 
c  create final printouts
 
c  release storage
 
      call gfree ("Fieldsrz")
      call gfree ("Hist")
      call gfree ("Win_Moments")
      call gfree ("Moments")
      call gfree ("Particles")
      call gfree ("Z_arrays")
 
      return
      end
c=============================================================================
      subroutine steprz (caller)
      use Beam_acc
      use Constant
      use InGen
      use InGenrz
      use InDiag
      use InPart
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      use Io
      use Particles,Only: pgroup
      use Picglb
      use Moments
      use Z_arrays
      use FFTDiagrz
      character(*):: caller
      logical(ISZ):: thisstep,thiszbeam
 
c  Advances the system forward in time, when called by WRZEXE
c  Takes a step of zero size, to compute fields, diagnostics,
c  at start of run when called by WRZGEN.

      real(kind=8):: zbeaml,zbeamr,er,et,ez,frext
 
c  Set logicals
 
      lspecial = (lspecial .or. lfirst)
 
c  Main particle advance: x to t.l. it; v to t.l. it-1/2
c  Half-step in v from t.l. it-1 if last step was "special"
c  Full-step in v from t.l. it-3/2 if last step not "special"
c  No step at all if only generating.
 
      if (caller .eq. "wrzexe") then
         if (lspecial) then
            call padvncrz ("halfv   ")
         else
            call padvncrz ("fullv   ")
         endif
      endif

c  The next two variables are the left and right ends of the range centered
c  about the end of the current time step plus/minus one half a step.
c  The range is used is determining whether diagnostics are done which
c  are based on the z location of the beam frame.  The diagnostics are done
c  on the time step which ends closest to the value given in the controlling
c  arrays.
      zbeaml   = zbeam - 0.5*vbeam*dt
      zbeamr   = zbeam + 0.5*vbeam*dt

c  "Frequent" or "infreq" phase space plots or restart dump, or done,
c  at end of this step? Set flags
 
      lfinishd = (it .ge. nt) .or. (time .ge. tstop-fuz)
      lalways  = thisstep (it           ,itplalways,NCONTROL) .or.
     &           thiszbeam(zbeaml,zbeamr,zzplalways,NCONTROL) .or.
     &           thisstep (it           ,itplfreq,  NCONTROL) .or.
     &           thiszbeam(zbeaml,zbeamr,zzplfreq,  NCONTROL)
      lseldom  = thisstep (it           ,itplseldom,NCONTROL) .or.
     &           thiszbeam(zbeaml,zbeamr,zzplseldom,NCONTROL) .or.
     &           thisstep (it           ,itplps,    NCONTROL) .or.
     &           thiszbeam(zbeaml,zbeamr,zzplps,    NCONTROL)
      lmoments = thisstep (it           ,itmomnts,NCONTROL) .or.
     &           thiszbeam(zbeaml,zbeamr,zzmomnts,NCONTROL)
      ldump = mod(it, itdump) .eq. 0
      lspecial = (lfinishd .or. lalways .or. lseldom .or. ldump .or. lmoments
     &            .or. llast .or. allspecl)

c  Compute line charge density, vz-of-z, "gap" e field on 1d meshes.
c  Note that vzofz is computed from curr, linechg; also curr was computed
c  using "old" (half-step-behind) vz's.  We might compute linechg less often.
c  Also compute on-axis quantities.

      if (lspecial) then
      call gtlchgrz
      call srhoaxrz
      endif

c  Charge density and line charge plots
      if (lalways .or. lseldom) then
         call pltfldrz("rho",ALWAYS)
      endif
      if (lseldom) then
         call pltfldrz("rho",SELDOM)
      endif
 
c  Field-solve for potential (copy rho into phi first, compute e.s.e. after)
      if (lbeforefs) call execuser("beforefs")
      call copyarry (rho, phi(0,0), nrz)
      call vprz (-1)
      if (periinz) call perphirz(phi(0,-1),nr,nz)
      if (lafterfs) call execuser("afterfs")
      call geteserz

      if ((caller .eq. "wrzgen").and.(frcetype .eq. "regress")) then
         call getforce
      endif
      if ((caller .eq. "wrzgen").and.(ifeears .eq. 3)) then
c  Call seterz to fill E_z and E_r arrays
         call seterz (-1,phi,force,1,pgroup%xp(1),pgroup%zp(1),zbeam,
     &                rmmin,zmmin,zmmax,dr,dz,nr,nz,
     &                efetch,er,et,ez,frext,erfld,ezfld,eoffrz)
         call getforce
         call setexte(it,dt,dz,vbeam,nr,nz)
      endif
 
c  Electrostatic potential plots
      if (lalways .or. lseldom) call pltfldrz("phi",ALWAYS)
      if (lseldom) call pltfldrz("phi",SELDOM)

c  If a flag was set making this a "special" step,
c  we'll do a half-advance now to bring v to t.l. it
      if (caller .eq. "wrzexe" .and. lspecial) then
         call padvncrz ("synchv  ")
      endif

      if (lspecial) then
      call getvzofz
      call setegap
      call sphiaxrz
      call sezaxrz
      endif

c  1d array plots.

      if (lalways .or. lseldom) call onedplts(ALWAYS)
      if (lseldom) call onedplts(SELDOM)

c  Phi on axis and Vz vs Z plot in one frame for videos
c     if (lalways .or. lseldom) call ppphivz(ALWAYS)
c     if (lseldom) call ppphivz(SELDOM)

c  Phi on axis and Vz vs Z plot in one frame in color for videos
c     if (lalways .or. lseldom) call ppphivzc(ALWAYS)
c     if (lseldom) call ppphivzc(SELDOM)
 
c  Phase space diagnostics
      if (lalways .or. lseldom) call psplots (ALWAYS)
      if (lseldom) call psplots (SELDOM)
 
c  Finally, final diagnostics and history stuff
c  (put here so frame index looks nicer)

      if (caller .eq. "wrzexe") then
         if (lspecial) then
            if (fftdiag) call fftdrz 
            call minidiag (it,time,lmoments)
         endif
      else
         if (fftdiag) call fftdrz 
         call minidiag (it,time,lmoments)
      endif
 
      return
      end
c=============================================================================
      subroutine extbrz (np,rp,tp,zp,uzp,gaminv,br,bt,bz)
      use Beam_acc
      use InGen
      integer(ISZ):: np
      real(kind=8):: rp(np), tp(np), zp(np)
      real(kind=8):: uzp(np), gaminv(np)
      real(kind=8):: br(np), bt(np), bz(np)
 

      integer(ISZ):: ip
 
 
c -- Set Bz to axial field
      do ip = 1,np
 
        br(ip) = 0.
        bt(ip) = 0.
        bz(ip) = bz0
 
      enddo
 
      return
      end
c=============================================================================
      subroutine exterz (np,zp,uzp,zbeam,rp,ez0,er,et,ez,dedr)
      integer(ISZ):: np
      real(kind=8):: zbeam,ez0,dedr
      real(kind=8):: zp(np), uzp(np), rp(np)
      real(kind=8):: er(np), et(np), ez(np)
 
c  Set the electric fields from the quads etc.
c  Also zeros the electric field if particle is lost (uzp(ip)=0).

      integer(ISZ):: ip
 
      do ip=1,np

        if (dedr .ne. 0.) then
          er(ip) = er(ip) + dedr*rp(ip)
        endif

        if (uzp(ip) == 0) then
          er(ip) = 0.
          et(ip) = 0.
          ez(ip) = 0.
        endif

      enddo
 
      return
      end
c=============================================================================
      subroutine fixrhorz
 
c  Corrects charge density for curvature
 
      return
      end
c=============================================================================
      subroutine othererz (np,xp,yp,zp,zbeam,zimax,zimin,straight,
     &                   ifeears,eears,eearsofz,dr,dz,nr,nz,zmmin,
     &                   ex,ey,ez,eearsrz,rmmin,ezext)
      integer(ISZ):: np,ifeears,nr,nz
      real(kind=8):: zbeam,zimax,zimin,straight,eears,dr,dz,zmmin,rmmin
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: ex(np), ey(np), ez(np)
      real(kind=8):: eearsofz(0:nz),eearsrz(0:nr,0:nz),ezext(0:nz)

c  Set the electric fields from external sources

      integer(ISZ):: ip,iz,i,k
      real(kind=8):: zs,dzi,wz,r1,w1,r2,w2,dri

c  average axial confining ears for finite beam

      if (ifeears .eq. 1) then
        zs = (zimax - zimin)*straight/2.
        do ip=1,np
          if ((zp(ip)-zbeam) .gt. zs) then
            ez(ip) = ez(ip) + eears*(zp(ip) - zbeam - zs)
          elseif ((zp(ip)-zbeam) .lt. -zs) then
            ez(ip) = ez(ip) + eears*(zp(ip) - zbeam + zs)
          endif
        enddo
      endif

c  axial confining ears for finite beam as a function of z

      if (ifeears .eq. 2) then
        dzi = 1./dz
        do ip=1,np
          iz = (zp(ip) - zbeam - zmmin)*dzi
          wz = (zp(ip) - zbeam - zmmin)*dzi - iz
          ez(ip) = ez(ip) + eearsofz(iz)*(1.-wz) + eearsofz(iz+1)*wz
        enddo
      endif

c  Add external E_z that keeps beam from slowing down with resistive/
c  capacitive wall
      dzi = 1./dz
      do ip=1,np
        iz = (zp(ip) - zbeam - zmmin)*dzi
        wz = (zp(ip) - zbeam - zmmin)*dzi - iz
        ez(ip) = ez(ip) + ezext(iz)*(1.-wz) + ezext(iz+1)*wz
      enddo

c  axial confining ears for finite beam as a function of z and r

      if (ifeears .eq. 3) then
        dzi = 1./dz
        dri = 1./dr
        do ip = 1, np
         i =  (xp(ip) - rmmin) * dri
         k =  (zp(ip) - zbeam - zmmin) * dzi
         r1 = (xp(ip) - rmmin)  - i*dr
         w1 = (zp(ip) - zbeam - zmmin)  - k*dz
         r2 = dr - r1
         w2 = dz - w1
         w1 = w1 * dzi
         w2 = w2 * dzi
         r1 = r1 * dri
         r2 = r2 * dri
         ez(ip) = ez(ip) + (r2 *  w2 * eearsrz(i,k)
     &                   +  r1 *  w2 * eearsrz(i+1,k)
     &                   +  r2 *  w1 * eearsrz(i,k+1)
     &                   +  r1 *  w1 * eearsrz(i+1,k+1))
        enddo
      endif
      return
      end
c=============================================================================
      subroutine geteserz
      use Constant
      use Beam_acc
      use Picglbrz
      use InMeshrz
      use Fieldsrz
      use Moments
      use InPart
      use InGen
      use InGenrz
      use Z_arrays
 
c  Calculate the electrostatic energy, sum of rho*phi*r*pi
c  Includes external fields: Uniform focusing and eears for cigars
c  This does it in a partially vectorized manner.
 

      real(kind=8):: zm,phiextun
      integer(ISZ):: ir,iz
 
      ese = 0.
      zm = (zimax - zimin)*(1. - straight)/2.

c  Use scrtch(1,ir) to hold the cell volumes
      scrtch(1,0) = 0.25*pi*dr*dr*dz
      do ir=1,nr
         scrtch(1,ir) = 2.*pi*dr*dz*rmesh(ir)
      enddo

      do ir=0,nr
        phiextun = 0.
        if (ifeears .eq. 1) then
          if (zmmin .lt. -zm) phiextun = eears*.5*(zmmin + zm)**2
          if (zmmin .gt.  zm) phiextun = eears*.5*(zmmin - zm)**2
        endif
        if (ifeears .eq. 2) then
          phiextun = eearsofz(0)
        endif
        scrtch(0,ir) = rho(ir,0)*(0.5 * phi(ir,0) + phiextun
     &                       - dedr * 0.5 * rmesh(ir)**2)*scrtch(1,ir)
      enddo

      do iz=1,nz
        phiextun = 0.
        if (ifeears .eq. 1) then
          if (iz*dz+zmmin .lt. -zm) phiextun = eears*.5*(iz*dz+zmmin + zm)**2
          if (iz*dz+zmmin .gt.  zm) phiextun = eears*.5*(iz*dz+zmmin - zm)**2
        endif
        if (ifeears .eq. 2) then
          phiextun = eearsofz(iz)
        endif
        do ir=0,nr
          scrtch(0,ir) = scrtch(0,ir) + rho(ir,iz)*(0.5*phi(ir,iz) +
     &          phiextun - dedr * 0.5 * rmesh(ir)**2)*scrtch(1,ir)
        enddo
      enddo
 
      do ir=0,nr
        ese = ese + scrtch(0,ir)
      enddo
 
      return
      end
c=============================================================================
      subroutine getforce
      use Constant
      use InPart
      use Beam_acc
      use Particles,Only: pgroup
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      use InGen
      use InGenrz
      use LinearBeam
 

      real(kind=8):: density,sumr,sumr2,sume,sumer,sume2,denom,slope,yint
      integer(ISZ):: icount,ir,iz,irsave,irz,nbeam
 
c  Calculate the linear external force which balances the space-charge
c  force of the beam particles.
 
c  Method 1 -- calculate average density and then use analytic formula
c  for finding force
 
      if (frcetype .eq. "analytic") then
        if (.not. linbeam) then
         density = 0.0
         icount = 0
         do ir = 1,nr
           if (rmesh(ir).le. a0) then
              density = density + rho(ir,nz/2)
              icount = icount + 1
           endif
         enddo
 
         density = density/float(icount)
        else
         density = rho(0,nz/2)
        endif
 
         dedr = -pgroup%sq(1)*density*0.5/eps0
         do iz = 0,nz
         do ir = 0,nr
            force(ir,iz) = dedr * rmesh(ir)
         enddo
         enddo

      endif
 
C  Method 2 -- Do a linear regression of the electric field and
c  use this to calculate the force
 
      if (frcetype .eq. "regress") then
         dr = rmesh(3) - rmesh(2)
         do ir = 0,nr-1
            rhalf(ir) = 0.5*(rmesh(ir+1)+rmesh(ir))
         enddo
       do iz = 0,nz
         sumr = 0.
         sumr2= 0.
         do ir = 0,nr-1
c           if (rhalf(ir).lt.a0) then
            if (rho(ir,iz) .gt. 0.) then
               irsave = ir
               sumr = sumr + rhalf(ir)
               sumr2 = sumr2 + rhalf(ir)**2
            endif
         enddo
 
 
         sume = 0.0
         sumer = 0.0
         sume2 = 0.
         do ir = 0,irsave
            if (rho(ir,iz) .gt. 0) irz = ir
         enddo

         do ir = 0,irz
            force(ir,iz) = pgroup%sq(1)*(phi(ir,iz) - phi(ir+1,iz))/dr
         enddo
 
         do ir = 0,irz
            sumer = sumer + force(ir,iz)*rhalf(ir)
            sume  = sume  + force(ir,iz)
            sume2 = sume2 + force(ir,iz)**2
         enddo
 
         nbeam = irz + 1
         denom = nbeam*sumr2 - sumr**2
         if (denom .ne. 0.) then
         slope = (nbeam*sumer - sumr*sume)/denom
         yint = (sume*sumr2 - sumer*sumr)/denom
         if (iz .eq. nz/2) dedr = -slope
         else
         slope = 0.
         yint = 0.
         endif
 
         do ir = 0,nr
            force(ir,iz) = -(slope*rmesh(ir) + yint)
         enddo

 
       enddo
      endif

c  If forcetype is dedr, then external force is added via exterz
      if (frcetype .eq. "dedr") then
         do iz = 0,nz
         do ir = 0,nr
            force(ir,iz) = 0.0
         enddo
         enddo
      endif

      if (ifeears .eq. 3) then
         do iz = 0,nz
         do ir = 0,nr
            force(ir,iz) = -pgroup%sq(1)*erfld(ir,iz)
         enddo
         enddo
      endif
 
      return
      end
 
c=============================================================================
      subroutine gtlchgrzold
      use Constant
      use Picglbrz
      use InMeshrz
      use InGenrz
      use Fieldsrz
      use Z_arrays
 
c  Calculates the line charge density from rho.


      integer(ISZ):: iz,ir
 
      do iz=0,nz
        linechg(iz) = 0.0
        do ir=0,nr
          linechg(iz) = linechg(iz) + rho(ir,iz)*rmesh(ir)
        enddo
 
        linechg(iz) = linechg(iz)*2.*pi*dr
 
      enddo
 
      return
      end
c=============================================================================
      subroutine padvncrz(center)
      use GlobalVars
      use Beam_acc
      use InGen
      use InPart
      use InGenrz
      use InMeshrz
      use InGaps
      use Particles,Only: pgroup,nplive,wpid
      use Picglbrz
      use Picglb
      use Fieldsrz
      use Z_arrays
      use LinearBeam
      character(*):: center
 
c  Advances the particles position and velocity according to CENTER,
c  and also advances RHO.
 
 
      integer(ISZ):: is,ipmin,ip,i
      real(kind=8):: er(nparpgrp), et(nparpgrp), ez(nparpgrp)
      real(kind=8):: br(nparpgrp), bt(nparpgrp), bz(nparpgrp)
      real(kind=8):: frext(nparpgrp)
      real(kind=8):: uxpo(nparpgrp), uypo(nparpgrp), uzpo(nparpgrp)
 
c  Zero RHO if particles are advanced
      if (center.eq."fullv" .or. center.eq."halfv") then
c        call zeroarry(rho,nrz)
         rho = 0.
      endif
c     call zeroarry(curr,nz+1)
      curr = 0.
 
c  Zero the moments if center is synchv
      if (center .eq. "synchv") then
c       call getemits(1,xp,yp,pgroup%zp,pgroup%uxp,pgroup%uyp,pgroup%uzp,pgroup%gaminv,pgroup%sm(1),pgroup%sw(1),1,nplive)
        call getzmmnt(1,pgroup%xp,pgroup%yp,pgroup%zp,pgroup%uxp,pgroup%uyp,
     &                pgroup%uzp,pgroup%gaminv,pgroup%sq(1),pgroup%sm(1),
     &                pgroup%sw(1),dt/2.,1,
     &                nplive,uxpo,uypo,uzpo,1,ns)
      endif
 
c  Loop over species
      do is=1,ns
c  Call seterz to fill E_z and E_r arrays
         if (it .le. 1) then
         call seterz (-1,phi,force,1,pgroup%xp(1),pgroup%zp(1),zbeam,
     &                rmmin,zmmin,zmmax,dr,dz,nr,nz,
     &                efetch,er,et,ez,frext,erfld,ezfld,eoffrz)
         else
         call seterz (0,phi,force,1,pgroup%xp(1),pgroup%zp(1),zbeam,
     &                rmmin,zmmin,zmmax,dr,dz,nr,nz,
     &                efetch,er,et,ez,frext,erfld,ezfld,eoffrz)
         endif
         call setexte(it,dt,dz,vbeam,nr,nz)
 
c  Loop over particle block, and move each block seperately
         do ipmin = pgroup%ins(is), pgroup%ins(is) + pgroup%nps(is) - 1, nparpgrp
            ip = min(nparpgrp, pgroup%ins(is)+pgroup%nps(is)-ipmin)
            call seterz (1,phi,force,ip,pgroup%xp(ipmin),pgroup%zp(ipmin),zbeam,
     &                   rmmin,zmmin,zmmax,dr,dz,nr,nz,
     &                   efetch,er,et,ez,frext,erfld,ezfld,eoffrz)
            call othererz(ip,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                    zbeam,zimax,zimin,straight,ifeears,eears,
     &                    eearsofz,dr,dz,nr,nz,zmmin,er,et,ez,
     &                    eearsrz,rmmin,ezext)
            call exterz (ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),zbeam,pgroup%xp(ipmin),
     &                   ez0,er,et,ez,dedr)
            call extbrz (ip,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                   pgroup%uzp(ipmin),pgroup%gaminv(ipmin),br,bt,bz)
            call pushrz (center, ip, pgroup%xp(ipmin), pgroup%yp(ipmin), pgroup%zp(ipmin),
     &                   pgroup%uxp(ipmin), pgroup%uyp(ipmin), pgroup%uzp(ipmin), pgroup%gaminv(ipmin),
     &                   er, et, ez, br, bt, bz, frext, pgroup%sq(is), pgroup%sm(is),
     &                   dt, gamadv)
            if (center .eq. "synchv") then
c             --- copy 'old' velocity into temp array uxpo,uypo, and uzpo for
c             --- z moments
               do i=1,ip
                 uxpo(i) = pgroup%uxp(ipmin+i-1)
                 uypo(i) = pgroup%uyp(ipmin+i-1)
                 uzpo(i) = pgroup%uzp(ipmin+i-1)
               enddo
c              call getemits(ip,xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
c    &                       pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
c    &                       pgroup%gaminv(ipmin),pgroup%sm(is),pgroup%sw(is),2,nplive)
               if(wpid==0) then
                 call getzmmnt(ip,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                         pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     & pgroup%gaminv(ipmin),pgroup%sq(1),pgroup%sm(is),pgroup%sw(is),dt/2.,2,nplive,
     &                         uxpo,uypo,uzpo,is,ns)
               else
                 call getzmmnt_weights(ip,pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                         pgroup%uxp(ipmin),pgroup%uyp(ipmin),pgroup%uzp(ipmin),
     &                         pgroup%gaminv(ipmin),pgroup%pid(ipmin,wpid),
     &                         pgroup%sq(1),pgroup%sm(is),pgroup%sw(is),dt/2.,2,nplive,
     &                         uxpo,uypo,uzpo,is,ns)
               endif
               call setcurr(curr,ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),pgroup%gaminv(ipmin),
     &                      pgroup%sq(is),pgroup%sw(is),zgrid,dz,zmmin)
            endif
         enddo
         if (center .eq. "fullv" .or. center .eq. "halfv") then
            if (periinz) call periz(pgroup%nps(is),pgroup%zp(pgroup%ins(is)),zgrid,zmmax,zmmin)
            if (stickyz)
     &         call stckyz(pgroup%nps(is),pgroup%zp(pgroup%ins(is)),zmmax,zmmin,dz,pgroup%uxp(pgroup%ins(is)),
     & pgroup%uyp(pgroup%ins(is)),pgroup%uzp(pgroup%ins(is)),pgroup%gaminv(pgroup%ins(is)),zgrid)
            if (stickyxy) call stckyrz(pgroup%nps(is),
     &                                pgroup%xp(pgroup%ins(is)),rmmax,rmmin,dr,
     &                                 pgroup%zp(pgroup%ins(is)),zmmin,dz,
     & pgroup%uxp(pgroup%ins(is)),pgroup%uyp(pgroup%ins(is)),pgroup%uzp(pgroup%ins(is)),
     &                                 pgroup%gaminv(pgroup%ins(is)),
     &                                 zgrid,pgroup%ins(is),it)
            call processlostpart(is,clearlostpart,time,zbeam)
            do ipmin = pgroup%ins(is), pgroup%ins(is) + pgroup%nps(is) - 1, nparpgrp
               ip = min(nparpgrp, pgroup%ins(is)+pgroup%nps(is)-ipmin)
               call setrhorz (rho, ip, pgroup%xp(ipmin), pgroup%zp(ipmin), zgrid,
     &                        pgroup%uzp(ipmin), pgroup%sq(is), pgroup%sw(is), depos)
              if (ifgap) 
     &         call setcurr(curr,ip,pgroup%zp(ipmin),pgroup%uzp(ipmin),pgroup%gaminv(ipmin),
     &                      pgroup%sq(is),pgroup%sw(is),zgrid,dz,zmmin)
            enddo
         endif

         if (vecrho) then
c          call zeroarry(rhov,nrz)
           rhov = 0.
           call psort(pgroup%nps(is),zgrid,dz)
           call vstrhorz(rhov,zgrid,pgroup%sq(is),pgroup%sw(is),rmesh)
           call perrhorz(rhov,nr,nz)
         endif

c        ---  Do final stuff for moments calculation
         if (center .eq. "synchv") then
c         call getemits(1,xp,yp,pgroup%zp,pgroup%uxp,pgroup%uyp,pgroup%uzp,pgroup%gaminv,pgroup%sm(is),pgroup%sw(is),3,nplive)
          call getzmmnt(1,pgroup%xp,pgroup%yp,pgroup%zp,pgroup%uxp,pgroup%uyp,
     &                  pgroup%uzp,pgroup%gaminv,pgroup%sq(1),pgroup%sm(is),
     &                  pgroup%sw(is),dt/2.,
     &                  3,nplive,uxpo,uypo,uzpo,1,ns)
         endif
 
c        --- Make rho periodic if periodic in z
         if ((center.eq."fullv" .or. center.eq."halfv") .and. periinz)
     &         call perrhorz(rho,nr,nz)
 
c        --- Set center of beam, if particles were advanced
         if (center.eq."fullv" .or. center.eq."halfv") zbeam = zgrid
 
      enddo
c        --- calculate rho diagnostics
         if (center .eq. "synchv") call rhodiarz
 
      return
      end
c=============================================================================
      subroutine perphirzold(phi,nr,nz)
      integer(ISZ):: nr,nz
      real(kind=8):: phi(0:nr,-1:nz+1)
 
c  Sets the first and last slices of phi for periodicity

      integer(ISZ):: ir
 
      do ir=0,nr
          phi(ir,-1)   = phi(ir,nz-1)
          phi(ir,nz+1) = phi(ir,1)
      enddo
 
      return
      end
c=============================================================================
      subroutine perrhorzold(rho,nr,nz)
      integer(ISZ):: nr,nz
      real(kind=8):: rho(0:nr,0:nz)
 
c  Adds the first and last slices of rho for periodicity

      integer(ISZ):: ir
 
      do ir=0,nr
       rho(ir,0)  = (rho(ir,0) + rho(ir,nz))
       rho(ir,nz) = rho(ir,0)
      enddo
 
      return
      end
c=============================================================================
      subroutine pushrz (center, np, rp, tp, zp, urp, utp, uzp, gaminv,
     &                 er, et, ez, br, bt, bz, frext, q, m, dt, gamadv)
      use GlobalVars
      use Constant
      use Beam_acc
      use InMeshrz
      character(*):: center ,gamadv
      integer(ISZ):: np
      real(kind=8):: rp(np),tp(np),zp(np),urp(np),utp(np),uzp(np),gaminv(np)
      real(kind=8):: er(np),et(np),ez(np),br(np),bt(np),bz(np),frext(np)
      real(kind=8):: q,m,dt
      real(kind=8):: tmp(np)
 
c  Advances the particles
 

      integer(ISZ):: ip
      real(kind=8):: tr,tt,tz,tsq,sr,st,sz,urppr,utppr,uzppr,urppl,utppl,uzppl
      real(kind=8):: urpmi,utpmi,uzpmi
      real(kind=8):: x2p,y2p
 
c  -------  Half Velocity Advance  -------
c  from it-1 to it-1/2
 
      if (center .eq. "halfv") then
 
        do ip=1,np
          tr = gaminv(ip)*q*br(ip)*dt/(4.*m)
          tt = gaminv(ip)*q*bt(ip)*dt/(4.*m)
          tz = gaminv(ip)*q*bz(ip)*dt/(4.*m)
          tsq = tr**2 + tt**2 + tz**2
          sr = 2.*tr/(1. + tsq)
          st = 2.*tt/(1. + tsq)
          sz = 2.*tz/(1. + tsq)
 
          urppr = urp(ip) + utp(ip)*tz - uzp(ip)*tt
          utppr = utp(ip) + uzp(ip)*tr - urp(ip)*tz
          uzppr = uzp(ip) + urp(ip)*tt - utp(ip)*tr
 
          urppl = urp(ip) + utppr*sz - uzppr*st
          utppl = utp(ip) + uzppr*sr - urppr*sz
          uzppl = uzp(ip) + urppr*st - utppr*sr
 
          urp(ip) = urppl + (q*er(ip)+frext(ip))*dt/(2.*m)
          utp(ip) = utppl + q*et(ip)*dt/(2.*m)
          uzp(ip) = uzppl + q*ez(ip)*dt/(2.*m)
        enddo
 
      endif
 
c  -------  Synchronize Velocity With Position  -------
c  from it-1/2 to it
 
      if (center .eq. "synchv") then
 
        do ip=1,np
          urpmi = urp(ip) + (q*er(ip)+frext(ip))*dt/(2.*m)
          utpmi = utp(ip) + q*et(ip)*dt/(2.*m)
          uzpmi = uzp(ip) + q*ez(ip)*dt/(2.*m)
 
          tr = gaminv(ip)*q*br(ip)*dt/(4.*m)
          tt = gaminv(ip)*q*bt(ip)*dt/(4.*m)
          tz = gaminv(ip)*q*bz(ip)*dt/(4.*m)
          tsq = tr**2 + tt**2 + tz**2
          sr = 2.*tr/(1. + tsq)
          st = 2.*tt/(1. + tsq)
          sz = 2.*tz/(1. + tsq)
 
          urppr = urpmi + utpmi*tz - uzpmi*tt
          utppr = utpmi + uzpmi*tr - urpmi*tz
          uzppr = uzpmi + urpmi*tt - utpmi*tr
 
          urp(ip) = urpmi + utppr*sz - uzppr*st
          utp(ip) = utpmi + uzppr*sr - urppr*sz
          uzp(ip) = uzpmi + urppr*st - utppr*sr
        enddo
 
      endif
 
c  -------  Full Velocity Advance  -------
c  from it-3/2 to it-1/2
 
      if (center .eq. "fullv") then
 
        do ip=1,np
          urpmi = urp(ip) + (q*er(ip)+frext(ip))*dt/(2.*m)
          utpmi = utp(ip) + q*et(ip)*dt/(2.*m)
          uzpmi = uzp(ip) + q*ez(ip)*dt/(2.*m)
 
          tr = gaminv(ip)*q*br(ip)*dt/(2.*m)
          tt = gaminv(ip)*q*bt(ip)*dt/(2.*m)
          tz = gaminv(ip)*q*bz(ip)*dt/(2.*m)
          tsq = tr**2 + tt**2 + tz**2
          sr = 2.*tr/(1. + tsq)
          st = 2.*tt/(1. + tsq)
          sz = 2.*tz/(1. + tsq)
 
          urppr = urpmi + utpmi*tz - uzpmi*tt
          utppr = utpmi + uzpmi*tr - urpmi*tz
          uzppr = uzpmi + urpmi*tt - utpmi*tr
 
          urppl = urpmi + utppr*sz - uzppr*st
          utppl = utpmi + uzppr*sr - urppr*sz
          uzppl = uzpmi + urppr*st - utppr*sr
 
          urp(ip) = urppl + (q*er(ip)+frext(ip))*dt/(2.*m)
          utp(ip) = utppl + q*et(ip)*dt/(2.*m)
          uzp(ip) = uzppl + q*ez(ip)*dt/(2.*m)
        enddo
 
      endif
 
c  -------  Gamma Advance     -------
      call gammaadv(np,gaminv,urp,utp,uzp,gamadv,lrelativ)
c  -------  Advance Position  -------
c  from it-1 to it
 
      if (center .eq. "halfv" .or. center .eq. "fullv") then
 
c  Update the positions using the particle advance given in Birdsall and
c  Langdon p 338, then rotate velocities to match coordinates of the
c  particle -- See diagram p 339 of B&L
        do ip=1,np
          x2p = rp(ip) + urp(ip)*dt*gaminv(ip)
          y2p = utp(ip)*dt*gaminv(ip)
          rp(ip) = sqrt(x2p**2 + y2p**2)
          zp(ip) = zp(ip) + uzp(ip)*dt*gaminv(ip)
          tp(ip) = 0.0
          if (rp(ip) .ne. 0.0) then
             tmp(ip) =  urp(ip)*x2p/rp(ip) + utp(ip)*y2p/rp(ip)
             utp(ip) = -urp(ip)*y2p/rp(ip) + utp(ip)*x2p/rp(ip)
c            tp(ip) = tp(ip) + asin(y2p/rp(ip))
          endif
        enddo

        do ip=1,np
           if (rp(ip) .ne. 0.0) then
              urp(ip) = tmp(ip)
           endif
        enddo

        do ip = 1,np
        if (rp(ip) .lt. 0.) then
             rp(ip) = -rp(ip)
        endif
        enddo
 
      endif
 
      return
      end
c=============================================================================
      subroutine seterz (itask,phi,force,np,rp,zp,zbeam,rmmin,zmmin,zmmax,
     &                   dr,dz,nr,nz,efetch,er,et,ez,frext,
     &                   erfld,ezfld,eoffrz)
      use Constant
      use LinearBeam
      use InPart
      integer(ISZ):: np,itask,nr,nz
      real(kind=8):: zbeam,rmmin,zmmin,zmmax,dr,dz,eoffrz
      real(kind=8):: phi(0:nr,-1:nz+1), force(0:nr,0:nz), rp(np), zp(np)
      real(kind=8):: er(np), et(np), ez(np), frext(np)
      real(kind=8):: erfld(0:nr,0:nz), ezfld(0:nr,0:nz)
      integer(ISZ):: efetch
 
      integer(ISZ):: ip,iz,ir,i,k
      real(kind=8):: tdri,tdzi,r1,w1,r2,w2

c  Sets electric field and external force for a group of up to 64 particles
 
 
c  Algorithm notes: phi array is dimensioned (0:nr,0:nz),
c  so cell index into 1d phi array for vectorized deposition is:
c     i + k*(nr+1)
c  The field is:
c     Er = u2 *  w2 * er(i  ,k  )
c        + u1 *  w2 * er(i+1,k  )
c        + u2 *  w1 * er(i  ,k+1)
c        + u1 *  w1 * er(i+1,k+1)
c     ...
 
c  Evaluate E field on the grid
 
      tdri = 1. / (2. * dr)
      tdzi = 1. / (2. * dz)

c  itask = 0; fill the field arrays 
      if (itask .le. 0) then

c  Do the boundaries separately -- phi(-1,iz) = phi(1,iz) since the
c  system is cylindrically symmetric.
c     if ((.not. shell) .and. (.not. linbeam)) then 
      if (.not. shell) then 
      do iz = 0,nz
         erfld(0,iz) = tdri*(phi(1,iz) - phi(1,iz))
         ezfld(0,iz) = tdzi*(phi(0,iz-1) - phi(0,iz+1))
         erfld(nr,iz) = (1./dr)*(phi(nr-1,iz) - phi(nr,iz))
         ezfld(nr,iz) = tdzi*(phi(nr,iz-1) - phi(nr,iz+1))
      enddo
 
      do ir = 1,nr-1
      do iz = 1,nz
         erfld(ir,iz) = tdri*(phi(ir-1,iz) - phi(ir+1,iz))
         ezfld(ir,iz) = tdzi*(phi(ir,iz-1) - phi(ir,iz+1))
      enddo
      enddo

      do ir = 1,nr-1
         ezfld(ir,0) = tdzi*(phi(ir,nz-1) - phi(ir,1))
         ezfld(ir,nz) = ezfld(ir,0)
         erfld(ir,0) = tdri*(phi(ir-1,0) - phi(ir+1,0))
         erfld(ir,nz) = erfld(ir,0)
      enddo


c  If we've loaded a shell, use backward differencing to get the
c  radial electric field right
      else

      do iz = 0,nz
c        erfld(0,iz) = tdri*(phi(1,iz) - phi(1,iz))
         erfld(0,iz) = 0.0
         ezfld(0,iz) = tdzi*(phi(0,iz-1) - phi(0,iz+1))
c        erfld(nr,iz) = (1./dr)*(phi(nr-1,iz) - phi(nr,iz))
         erfld(nr,iz) = 0.0
         ezfld(nr,iz) = tdzi*(phi(nr,iz-1) - phi(nr,iz+1))
      enddo
 
      do ir = 1,nr-1
      do iz = 0,nz
c        erfld(ir,iz) = (phi(ir-1,iz) - phi(ir,iz))/dr
         erfld(ir,iz) = 0.0
         ezfld(ir,iz) = tdzi*(phi(ir,iz-1) - phi(ir,iz+1))
      enddo
      enddo
      endif

      endif

      if ((itask .eq. -1).and.(zmmax.ne.zimax).and.(zmmin.ne.zimin))
     &               eoffrz = ezfld(0,1)

      if (itask .le. 0) then
      do ir = 0,nr
      do iz = 0,nz
         ezfld(ir,iz) = ezfld(ir,iz) - eoffrz
      enddo
      enddo
      endif
 
c  Interpolate the electric field and the external force to the particles ...
c  We use volume weighting scheme in cylindrical coordinates.
 
      if (itask .eq. 1) then

      do ip = 1, np
 
         i =  (rp(ip) - rmmin) / dr
         k =  (zp(ip) - zbeam - zmmin) / dz
 
         r1 = (rp(ip) - rmmin)  - i*dr
         w1 = (zp(ip) - zbeam - zmmin)  - k*dz
 
         r2 = dr - r1
         w2 = dz - w1
 
         w1 = w1 / dz
         w2 = w2 / dz
 
         r1 = r1 / dr
         r2 = r2 / dr
 
         er(ip) = (r2 *  w2 * erfld(i,k)
     &           + r1 *  w2 * erfld(i+1,k)
     &           + r2 *  w1 * erfld(i,k+1)
     &           + r1 *  w1 * erfld(i+1,k+1))
 
         frext(ip) = (r2 *  w2 * force(i,k)
     &              + r1 *  w2 * force(i+1,k)
     &              + r2 *  w1 * force(i,k+1)
     &              + r1 *  w1 * force(i+1,k+1))
 
c   Set e_theta to zero ...
         et(ip) = 0.0
 
         ez(ip) = (r2 *  w2 * ezfld(i,k)
     &           + r1 *  w2 * ezfld(i+1,k)
     &           + r2 *  w1 * ezfld(i,k+1)
     &           + r1 *  w1 * ezfld(i+1,k+1))
 
      enddo

      endif
 
      return
      end
c=============================================================================
      subroutine setsc
      use Constant
      use InMeshrz
      use InGenrz
      use Fieldsrz
      use Picglbrz
c   Sets up the initial surface charge on the walls, using the average
c   charge density at each z point.


      integer(ISZ):: iz,ir
      real(kind=8):: vol
 
         vol = 2.*dz*dr*dr*pi
         do iz = 0,nz-1
c  the 1/8 is because the volume at r=0 is 1/8 *(2.*pi*dr*dr*dz)
           schrg(iz) = rho(0,iz)/8.
           do ir = 1,nr
            schrg(iz) = schrg(iz) + rho(ir,iz)*ir 
           enddo
         enddo
         do iz = 0,nz-1
            schrg(iz) = -schrg(iz)*vol
         enddo

c  Make the surface charge periodic
      schrg(nz) = schrg(0)
 
      return
      end
 
c=============================================================================
      subroutine stptclrz
      use GlobalVars
      use Beam_acc
      use Constant
      use InPart
      use InPartrz
      use InMeshrz
      use InGenrz
      use Picglbrz
      use Fieldsrz
      use Particles,Only: pgroup,npmax,np
      use Setpworkrz
      integer(ISZ):: envxport
      character(72):: errline
 
c  Loads particles; sets npmax, np, etc.
c  Note that species 1 is the main (beam) ion species, and obtains
c  its charge, mass, etc. from global quantities such as aion, zion.
 
      integer(ISZ):: npm,ip,np_main,is,k_outer,i1_outer,k,i1,i,izstripe,icheck
      integer(ISZ):: nfib2,nfib3,nfib4,j,ii,ipart,ipmin,ipt
      real(kind=8):: rnpmi,strgt,zmid,zlen,zleft,vtz,rr,th,r0,rp,rpp,vtx,vty
      real(kind=8):: phi1,phi2,r,vr,vt,gfact,z0
      real(kind=8):: uxpert(0:nparpgrp),uzpert(0:nparpgrp)
      real(kind=8):: wranf,rm,rnrev
 
c  Set npgrp so Setpworkrz arrays are correct size
      if (xrandom .eq. "fibonacc") npgrp = nparpgrp
      if (xrandom .eq. "digitrev") npgrp = nparpgrp
      if (xrandom .eq. "pseudo"  ) npgrp = nparpgrp
      if (xrandom .eq. "grid"    ) npgrp = nxstripe*nystripe
      if (npgrp .eq. 0) then
        write (errline,'("ERROR: xrandom has an improper value = ",a8)')
     &         xrandom
        call remark (errline)
        call kaboom (0)
      endif
 
      call gallot ("Setpworkrz", 0)
 
c  Set npm to loop through all possible particles
      if (xrandom .eq. "fibonacc") npm = nfibgrps*fibg1
      if (xrandom .eq. "digitrev") npm = npmax
      if (xrandom .eq. "pseudo"  ) npm = npmax
      if (xrandom .eq. "grid"    ) npm = nxstripe*nystripe*nzstripe
 
      ip = 0
      np = 0
      np_main=npm-np_pert
      rnpmi = 1./np_main
      strgt=straight
 
c----------------------------------
c  Begin main loop over species
c----------------------------------
 
      do 1000 is = 1, ns
 
      zmid = .5 * (zimax + zimin)
      zlen = zimax - zimin
      zleft=zimin
c  ADDED LOOP AROUND MAIN LOOP TO ALLOW FOR A SECOND GROUP OF PARTICLES
c  CODE CURRENTLY IMPLEMENT FOR semi-gauss, cigar LOADING
      do 1001 k_outer = 1, npm, np_main
c
      i1_outer = min(npm, k_outer+np_main-1)
c 
c  MAIN INITIALIZATION LOOP
c 
      do k = k_outer, i1_outer, npgrp
         i1 = min(i1_outer, k+npgrp-1)
c
c  Load longitudinal stuff
c        --- load normalized z (0 < z < 1)
         if (     xrandom .eq. "fibonacc"
     &       .or. xrandom .eq. "digitrev"
     &       .or. xrandom .eq. "pseudo") then
           do i=k,i1
              zt(i-k+1) = (i-k_outer+.5)*rnpmi
           enddo
         elseif (xrandom .eq. "grid") then
           izstripe = (k-1)/npgrp + 1
           do i=k,i1
              zt(i-k+1) = (izstripe-.5)/nzstripe + zjig*(wranf()-.5)/nzstripe
           enddo
         endif

c        --- zero vz if vthz is zero only the first time through
         if (vthz .eq. 0. .and. k .eq. 1) then
           do i=1,i1-k+1
              uzt(i) = 0.
           enddo
         endif

         if (cigarld) then
c          --- set vtz to vthz*2 since vtz is vz_max - vz_bar for cigar
           vtz = 2.*vthz
c          --- load uzt with linear distribution (-.5 < uzt < .5)
           if (vtz .ne. 0.) then
             if (vzrandom .eq. "pseudo") then
               do i=1,i1-k+1
                 uzt(i) = wranf() - .5
               enddo
             elseif (vzrandom .eq. "digitrev") then
               do i=1,i1-k+1
                 uzt(i) = rnrev(i+k-1+randoffset,dig5) - .5
               enddo
             else
               write (errline,'("ERROR: vzrandom has an improper value = ",a8)')
     &               vzrandom
               call remark (errline)
               call kaboom (0)
             endif
           endif
c          --- use cigar to adjust the z's
c          --- at,bt,apt,bpt passed as scratch arrays
           call cigar(i1-k+1,zt,uzt,zt,uzt,perpscal,strgt,at,bt,apt,bpt)
c          --- readjust the perpendicular scale in the perturbation to
c              match the scale at the center of the main distribution.
           if(k .gt. np_main) then
             do i=1,i1-k+1
               perpscal(i) = 1.
             enddo
           endif

c          --- Loading gaussian longitudinal velocity distribution.
           if (distr_l .eq. "gaussian") then
             do i=1,i1-k+1
               uzt(i) = 0.5*perpscal(i)*rm()
             enddo
           endif

         else
c          --- non-cigar load
c          --- set vtz to vthz
           vtz = vthz
c          --- load uzt with gaussian distribution (mean 0, variance 1)
           if (vtz .ne. 0.) then
             if (vzrandom .eq. "pseudo") then
               do i=1,i1-k+1
                 uzt(i) = rm()
               enddo
             elseif (vzrandom .eq. "digitrev") then
               call rnormdig(k+randoffset,i1-k+1,dig5,dig4,0.,uzt)
             else
               write (errline,'("ERROR: vzrandom has an improper value = ",a8)')
     &               vzrandom
               call remark (errline)
               call kaboom (0)
             endif
           endif
c          --- set perpscal to one first time through
           if (k .eq. 1) then
             do i=1,i1-k+1
               perpscal(i) = 1.0
             enddo
           endif
         endif

c        --- unormalize z to fetch envelope
         do i=1,i1-k+1
           zt(i) = zleft + zlen*zt(i)
         enddo

c        --- fetch envelope
         if (.not. cylinder) then
           icheck = envxport(i1-k+1,zt,at,apt,bt,bpt,xct,xpct,yct,ypct,
     &                       vzt,epsxt,epsyt,ibeamt)
           if (icheck .eq. 1) then
              call remark ("ERROR: out-of-range z sent to ENVXPORT")
              call kaboom (0)
           endif
         else
c          --- set only first time through
           if (k .eq. 1) then
             do i=1,i1-k+1
               at(i)  = a0
               bt(i)  = b0
               apt(i) = 0.0
               bpt(i) = 0.0
             enddo
           endif
         endif
           

c  Load tranverse stuff
c  semi-gaussian distribution
         if (distrbtn .eq. "semigaus") then

c          --- load normalized tranverse space variables
           if (xrandom .eq. "pseudo") then
             do i=1,i1-k+1
               xt(i) = 2.*wranf() - 1.
               yt(i) = 2.*wranf() - 1.
             enddo
           elseif (xrandom .eq. "fibonacc") then
             nfib2 = nfibgrps*fibg2
             nfib3 = nfibgrps*fibg3
             do i=1,i1-k+1
               xt(i) = 2.*mod((nfib2*(i+k-2) + 0.5)*rnpmi+randoffset,1.0) - 1.
               yt(i) = 2.*mod((nfib3*(i+k-2) + 0.5)*rnpmi+randoffset,1.0) - 1.
             enddo
           elseif (xrandom .eq. "digitrev") then
             do i=1,i1-k+1
               xt(i) = 2.*rnrev(i+k-1+randoffset,dig1) - 1.
               yt(i) = 2.*rnrev(i+k-1+randoffset,dig2) - 1.
             enddo
           elseif (xrandom .eq. "grid") then
             do j=1,nystripe
               do i=1,nxstripe
                 xt(k+i-1+(j-1)*nxstripe) = 2.*(i-.5)/nxstripe - 1.
                 yt(k+i-1+(j-1)*nxstripe) = 2.*(j-.5)/nystripe - 1.
               enddo
             enddo
           endif

c          --- load normalized transverse velocity variables
           if (vtrandom .eq. "pseudo") then
             do i=1,i1-k+1
               uxt(i) = rm()
               uyt(i) = rm()
             enddo
           elseif (vtrandom .eq. "digitrev") then
             call rnormdig(k+randoffset,i1-k+1,dig3,dig4,0.,uxt)
             call rnormdig(k+randoffset,i1-k+1,dig4,dig3,0.,uyt)
           else
             write (errline,'("ERROR: vtrandom has an improper value = ",a8)')
     &           vtrandom
             call remark (errline)
             call kaboom (0)
           endif

c          --- use random numbers to load particles in polar coordinates
           if (ldprfile .eq. "polar") then
             do i=1,i1-k+1
               rr = 0.5*(xt(i) + 1.)
               if (hollow .eq. 2) then
                 rr = (1 + hollow_h)*rr/
     &                (hollow_h + sqrt(hollow_h**2 + (1. - hollow_h**2)*rr))
               endif
               rr = sqrt(rr)
               th = yt(i)*pi
               xt(i) = rr*cos(th)
               yt(i) = rr*sin(th)
               indx(i) = i
             enddo
             j = i1-k+1
           elseif (ldprfile .eq. "streamls") then
c          --- carve into cylinder
             j=0
             do i=1,i1-k+1
               if (xt(i)**2 + yt(i)**2 .lt. 1.) then
                 j=j+1
                 indx(j) = i
               endif
             enddo
           elseif (ldprfile .eq. "stripes") then
c          --- carve into normalized envelope
             j=0
             r0 = max(a0,b0)
             do i=1,i1-k+1
               if ((xt(i)*bt(i))**2 + (yt(i)*at(i))**2
     &                                .lt. (at(i)*bt(i)/r0)**2) then
                 j=j+1
                 indx(j) = i
                 xt(i) = xt(i)*r0/at(i)
                 yt(i) = yt(i)*r0/bt(i)
               endif
             enddo
           else
             write (errline,'("ERROR: ldprfile has an improper value = ",a8)')
     &           ldprfile
             call remark (errline)
             call kaboom (0)
           endif

c          --- transform to hollow beam of type one
c          --- f(r^2)  =    f0 * (r/rmax)^2                r < rmax/2
c          ---              f0 * (1 - (r/rmax)^2)/3        r > rmax/2
c          --- x and y are multiplied by 1.08 to keep rbar and rrms roughly
c          --- the same as in the uniform beam
           if (hollow .eq. 1 .and. ldprfile .eq. "streamls")  then
             do i=1,j
               ii = indx(i)
               rp = sqrt(xt(ii)**2 + yt(ii)**2)
               rpp = sqrt(0.5*sqrt(rp**2))
               if (rp .gt. 0.5) rpp = sqrt(1. - 0.5*sqrt(3. - 3.*rp**2))
               xt(ii) = rpp*xt(ii)/rp*1.08
               yt(ii) = rpp*yt(ii)/rp*1.08
             enddo
           endif

c          --- unnormalize everything and load into particle arrays
           do i=1,j
             ii = indx(i)
             pgroup%xp(ip+i) = at(ii)*xt(ii)*perpscal(ii)
             pgroup%yp(ip+i) = bt(ii)*yt(ii)*perpscal(ii)
             pgroup%zp(ip+i) = zt(ii)
             vtx = .5*vbeam*emit/at(ii)*perpscal(ii)**2
             vty = .5*vbeam*emit/bt(ii)*perpscal(ii)**2
             pgroup%uxp(ip+i) = vbeam*apt(ii)*xt(ii)*perpscal(ii) + vtx*uxt(ii)
             pgroup%uyp(ip+i) = vbeam*bpt(ii)*yt(ii)*perpscal(ii) + vty*uyt(ii)
             pgroup%uzp(ip+i) = vbeam*(1. + vtilt*(zmid - zt(ii))/zlen) + vtz*uzt(ii)
           enddo

c          --- add sinusoidal perturbation to pgroup%uzp
           if (vzperamp .ne. 0.) then
             do i=ip+1,ip+j
               pgroup%uzp(i) = pgroup%uzp(i) + vzperamp*cos(2.*pi*pgroup%zp(i)/vzperlam + vzperphs)
             enddo
           endif

c          --- increment number of particles by size of current group
           ip = ip + j

c  K-V distribution
         elseif (distrbtn .eq. "K-V") then

c          --- fetch random numbers and put into xt,yt, & uxt  temporarily
           if (xrandom .eq. "pseudo") then
             do i=1,i1-k+1
               xt(i)  = wranf()
               yt(i)  = wranf()
               uxt(i) = wranf()
             enddo
           elseif (xrandom .eq. "fibonacc") then
             nfib2 = nfibgrps*fibg2
             nfib3 = nfibgrps*fibg3
             nfib4 = nfibgrps*fibg4
             do i=1,i1-k+1
               xt(i)  = mod((nfib2*(i+k-2) + 0.5)*rnpmi+randoffset,1.0)
               yt(i)  = mod((nfib3*(i+k-2) + 0.5)*rnpmi+randoffset,1.0)
               uxt(i) = mod((nfib4*(i+k-2) + 0.5)*rnpmi+randoffset,1.0)
             enddo
           elseif (xrandom .eq. "digitrev") then
             do i=1,i1-k+1
               xt(i)  = rnrev(i+k-1+randoffset,dig1)
               yt(i)  = rnrev(i+k-1+randoffset,dig2)
               uxt(i) = rnrev(i+k-1+randoffset,dig3)
             enddo
           else
             call remark ("ERROR: xrandom cannot be grid for K-V load.")
             call kaboom (0)
           endif

c          --- load x,y,ux, and uy evenly onto a 4-D ellipsoid
           do i=1,i1-k+1
             rr = sqrt(xt(i))
             phi1 = 2.*Pi*yt(i)
             phi2 = 2.*Pi*uxt(i)
             pgroup%xp(ip+i) = rr*cos(phi1)*at(i)*perpscal(i)
             pgroup%yp(ip+i) = rr*sin(phi1)*bt(i)*perpscal(i)
             pgroup%zp(ip+i) = zt(i)
             rr = sqrt(1.-rr*rr)
             pgroup%uxp(ip+i) = vbeam*(emit*rr*cos(phi2)*perpscal(i)**2 +
     &                          pgroup%xp(ip+i)*apt(i))/at(i)
             pgroup%uyp(ip+i) = vbeam*(emit*rr*sin(phi2)*perpscal(i)**2 +
     &                          pgroup%yp(ip+i)*bpt(i))/bt(i)
             pgroup%uzp(ip+i) = vbeam*(1. + vtilt*(zmid - zt(i))/zlen) + vtz*uzt(i)
           enddo

c          --- add sinusoidal perturbation to pgroup%uzp
           if (vzperamp .ne. 0.) then
             do i=ip+1,ip+i1-k+1
               pgroup%uzp(i) = pgroup%uzp(i) + vzperamp*cos(2.*pi*pgroup%zp(i)/vzperlam + vzperphs)
             enddo
           endif

c          --- increment number of particles by size of current group
           ip = ip + i1-k+1

         else
           write (errline,'("ERROR: distrbtn has an improper value = ",a8)')
     &           distrbtn
           call remark (errline)
           call kaboom (0)
         endif

      enddo
c
c  Reset some of the parameters for the perturbed particles.
      if(np_pert.eq.0) then
        rnpmi = 0.
      else
        rnpmi = 1./np_pert
      endif
      zmid = pertz_ctr
      zlen = pertz_len
      zleft=zmid-zlen/2
      strgt = pertz_strgt
c
c  END OF LOOP TO SPLIT SPECIES INTO MAIN AND PERTURBATION
c
1001  continue
c
c  Set particle number, indices, etc.
 
      pgroup%ins(is) = np + 1
      pgroup%nps(is) = ip - np
      np = ip
      npmax = np
      if (is .eq. 1) then
         pgroup%sq(is) = zion * echarge
         pgroup%sm(is) = aion * amu
         if (vbeam .ne. 0.)
     &    pgroup%sw(is) = ibeam * (zimax-zimin) / (vbeam * echarge * zion * np)
      endif
 
c     Adjust weighting of particles since cigar() makes beam more dense,
c     also adjust for removal of particles used in perturbation.
      if (cigarld) then
        pgroup%sw(is) = pgroup%sw(is)*(straight + (1. - straight)*2./3.)
     .      *float(np)/float(np_main)
      endif

      do ipart = 1,np
         r = sqrt(pgroup%xp(ipart)**2 + pgroup%yp(ipart)**2)
         if (r .ne. 0.) then
             vr = ( pgroup%uxp(ipart)*pgroup%xp(ipart) + pgroup%uyp(ipart)*pgroup%yp(ipart))/r
             vt = (-pgroup%uxp(ipart)*pgroup%yp(ipart) + pgroup%uyp(ipart)*pgroup%xp(ipart))/r
             pgroup%uxp(ipart) = vr
             pgroup%uyp(ipart) = vt
         endif
c        pgroup%yp(ipart) = atan(pgroup%yp(ipart)/pgroup%xp(ipart))
         pgroup%yp(ipart) = 0.
         pgroup%xp(ipart) = r
      enddo
c   Add perturbation in velocity 
      if (vzpert .gt. 0.) then
         gfact = log(rmmax/a0)/(2.*pi*eps0)
         vphase = sqrt(gfact*pgroup%sq(1)*ibeam/(vbeam*pgroup%sm(1)))
         do ipmin = 1,np,nparpgrp
            ipt = min(nparpgrp,np - ipmin + 1)
            do ipart=0,ipt-1
            uxpert(ipart) = 0.5*pgroup%xp(ipart+ipmin)*vzpert*kzpert*
     &                   sin(kzpert*pgroup%zp(ipart+ipmin))
            uzpert(ipart) = vzpert*cos(kzpert*pgroup%zp(ipart+ipmin))
            enddo

            do ipart=0,ipt-1
               pgroup%uxp(ipart+ipmin) = pgroup%uxp(ipart+ipmin) + uxpert(ipart)
               pgroup%uzp(ipart+ipmin) = pgroup%uzp(ipart+ipmin) + uzpert(ipart)
            enddo
         enddo
      endif

      if (vzpert .lt. 0.) then

         gfact = log(rmmax/a0)/(2.*pi*eps0)
         vphase = sqrt(gfact*pgroup%sq(1)*ibeam/(vbeam*pgroup%sm(1)))
         z0 = 1./kzpert
         
         do ipart = 1,np
            pgroup%uzp(ipart) = pgroup%uzp(ipart) + (vzpert*vphase)/
     &                   (1.+(pgroup%zp(ipart)/z0)**2)
         enddo
      endif
            

c  If particle has left the region, put it back using periodic bc's

        do ipart = 1,np
          if (pgroup%zp(ipart) .gt. zmmax) then
              pgroup%zp(ipart) = pgroup%zp(ipart) + zmmin - zmmax
          endif

          if (pgroup%zp(ipart) .lt. zmmin) then
              pgroup%zp(ipart) = pgroup%zp(ipart) + zmmax - zmmin
          endif

          if (pgroup%xp(ipart) .lt. 0.) then
             pgroup%xp(ipart) = -pgroup%xp(ipart)
          endif
        enddo


 1000 continue
 
c--------------------------------
c  End main loop over species
c--------------------------------
 
      call gfree ("Setpworkrz")
 
      return
      end
c=============================================================================
      subroutine setrhorz (rho1d, np, rp, zp, zbeam, uzp, q, wght, depos)
      use Constant
      use InMeshrz
      use Picglbrz
      use InGenrz
      use Fieldsrz
      integer(ISZ):: np
      real(kind=8):: zbeam,q,wght
      real(kind=8):: rho1d (0:*)
      real(kind=8):: rp(np), zp(np), uzp(np)
      character(*):: depos
 
c  Sets charge density

 
c  Algorithm notes: rho array is dimensioned (0:nr,0:nz),
c  so cell index into 1d rho array for vectorized deposition is:
c     i + k*(nr+1)
c  In each case,
c     rho(i  ,k  ) = rho(i  ,k  ) + u0 *  w0 * g
c     rho(i+1,k  ) = rho(i+1,k  ) + u1 *  w0 * g
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  RHO must be zeroed in ZERORHO since it is not zeroed here (to allow
c  handling of blocks of particles at a time)
 
      integer(ISZ):: ip,ir,i,k,ind0,m
      real(kind=8):: g,r1,r2,w1,w2

c     --- For vectorized algorithm
      integer(ISZ):: moff(0:3), indx(0:3,0:256-1)
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(0:256-1), kk(0:256-1)
c     --- For both algorithms
      real(kind=8):: s(0:3,0:256-1)
 
c  Set rmesh(0) to be dr/8 so that loop will vectorize and then
c  set it back at the end.
      rmesh(0) =  dr / 8.
 
c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos .eq. "vector") then
 
c  Set up offset array for vectorized deposition:
 
      moff(0) = 0
      moff(1) = 1
      moff(2) = nr+1
      moff(3) = nr+2
 
c  Begin main loop over species, groups of 256 particles
 
      g = wght * q / (pi * (dr * dz)**2)
c     --- set all weights equal to zero
      do ip = 0, np-1
         s(0,ip) = 0.0
         s(2,ip) = 0.0
         s(1,ip) = 0.0
         s(3,ip) = 0.0
      enddo

c     --- vectorized loop to compute indices, weights
      do ip = 1, np
         ir = ip - 1
         i =  (rp(ip) - rmmin) / dr
         r1 = (rp(ip) - rmmin)      - i*dr
         r2 = dr - r1
         k =  (zp(ip) - zbeam - zmmin) / dz
         w1 = (zp(ip) - zbeam - zmmin)     - k*dz
         w2 = dz - w1
         ind0 = i + k*(nr+1)
         indx(0,ir) = ind0 + moff(0)
         indx(1,ir) = ind0 + moff(1)
         indx(2,ir) = ind0 + moff(2)
         indx(3,ir) = ind0 + moff(3)
c  The case of i=0 has been handled by setting rmesh(0) such that it
c  comes out right.  Otherwise this loop won't vectorize.

c  Set the weight ne to zero if the particle isn't "lost"
        if (uzp(ip) /= 0.0) then
         s(0,ir) = g * r2 * w2 / (2.*rmesh(i))
         s(2,ir) = g * r2 * w1 / (2.*rmesh(i))
         s(1,ir) = g * r1 * w2 / (2.*rmesh(i+1))
         s(3,ir) = g * r1 * w1 / (2.*rmesh(i+1))
        endif
      enddo
c     --- vectorized deposition over the 4 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
      do m = 0, 3
         do ir = 0, np-1
            rho1d(indx(m,ir)) = rho1d(indx(m,ir)) + s(m,ir)
         enddo
      enddo
 
c    --- density at r=0 is only 1/2 of what it should be since particles
c        in the region r=-dr to r=0 would contribute the same amount of
c        density as those in r=0 to r=dr
c     do i = 0,nz
c        rho1d(i*(nr+1)) = 2.*rho1d(i*(nr+1))
c     enddo
 
c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos .eq. "scalar") then
 
c  Begin main loop over species, groups of 64 particles
 
      g = wght * q / (pi * (dr * dz)**2)
c     --- vectorized loop to compute indices, weights
      do ip = 1, np
         ir = ip - 1
         ii(ir) =  (rp(ip) - rmmin) / dr
         r1 = (rp(ip) - rmmin)         - ii(ir)*dr
         r2 = dr - r1
         kk(ir) =  (zp(ip) - zbeam - zmmin) / dz
         w1 = (zp(ip) - zbeam - zmmin)      - kk(ir)*dz
         w2 = dz - w1
         s(0,ir) = g * r2 * w2 / (2.*rmesh(ii(ir)))
         s(1,ir) = g * r1 * w2 / (2.*rmesh(ii(ir)+1))
         s(2,ir) = g * r2 * w1 / (2.*rmesh(ii(ir)))
         s(3,ir) = g * r1 * w1 / (2.*rmesh(ii(ir)+1))
      enddo
c     --- scalar loop does the actual deposition
      do ir = 0, np-1
       if (uzp(ir+1) /= 0) then
         rho(ii(ir)  ,kk(ir)  )
     &  =rho(ii(ir)  ,kk(ir)  ) + s(0,ir)
         rho(ii(ir)+1,kk(ir)  )
     &  =rho(ii(ir)+1,kk(ir)  ) + s(1,ir)
         rho(ii(ir)  ,kk(ir)+1)
     &  =rho(ii(ir)  ,kk(ir)+1) + s(2,ir)
         rho(ii(ir)+1,kk(ir)+1)
     &  =rho(ii(ir)+1,kk(ir)+1) + s(3,ir)
       endif
      enddo
 
c    --- density at r=0 is only 1/2 of what it should be since particles
c        in the region r=-dr to r=0 would contribute the same amount of
c        density as those in r=0 to r=dr
c     do i = 0,nz
c        rho1d(i*(nr+1)) = 2.*rho1d(i*(nr+1))
c     enddo
 
      endif

c  Put rmesh(0) back
      rmesh(0) = rmmin
 
      return
      end
c=============================================================================
      subroutine stckyrz(np,rp,rmmax,rmmin,dr,
     &                   zp,zmmin,dz,urp,utp,uzp,gaminv,zgrid,ipmin,it)
      integer(ISZ):: np,ipmin,it
      real(kind=8):: rmmax,rmmin,dr,zmmin,dz,zgrid
      real(kind=8):: rp(np), zp(np)
      real(kind=8):: urp(np), utp(np), uzp(np), gaminv(np)
 
c  Enforces sticky b.c.'s on the r wall

      integer(ISZ):: ip
 
      do ip=1,np
        if (rp(ip).ge.rmmax-dr) gaminv(ip) = 0.
c       if (rp(ip).ge.rmmax-dr) then
c         rp(ip) = rmmax - .1*dr
c         zp(ip) = zmmin - .9999*dz + zgrid
c         urp(ip) = 0.0
c         utp(ip) = 0.0
c         uzp(ip) = 0.0
c       endif
      enddo
 
      return
      end
c=============================================================================
      subroutine fieldsolrz(iwhich)
      use InGen
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      integer(ISZ):: iwhich

c Does the complete fieldsolver, including copying rho into phi and setting
c boundary conditions.

      call copyarry (rho, phi(0,0), nrz)
      call vprz (-1)
      if (periinz) call perphirz(phi(0,-1),nr,nz)

      return
      end
c=============================================================================
      subroutine vprz (iwhich)
      use Picglbrz
      use InGenrz
      use InMeshrz
      use Beam_acc
      use InGen
      use Fieldsrz
      integer(ISZ):: iwhich
 
c  Interface to VPOISRZ using variables from database of package RZ
 
 
      call vpoisrz (iwhich, phi(0,0), phi(0,0), kzsq, schrg, eta,
     & phikold, taurc,
     & attz, filt, dt, vbeam, rmmax-rmmin, zmmax-zmmin,
     & nr, nz, rfsmat, rwwork, rwwork2, ibc)
      return
      end
c=============================================================================
      subroutine getindx(kzsq,kzpert,indxk,nz,dr)
      integer(ISZ):: nz,indxk
      real(kind=8):: dr
      real(kind=8):: kzsq(0:nz)
      real(kind=8):: kzpert

      integer(ISZ):: nz2,ikz
      real(kind=8):: xkm2,diff,diffmax

      nz2 = nz/2
      xkm2 = (kzpert*dr)**2
      diffmax = LARGEPOS

      do ikz = 0,nz2
        diff = abs(xkm2 - kzsq(ikz))
        if (diff .lt. diffmax) then
           indxk = ikz
           diffmax = diff
        endif
      enddo

      return
      end

c=============================================================================
      subroutine fftdrz
      use InGenrz
      use FFTDiagrz
      use Z_arrays

      integer(ISZ):: iz
      
      do iz = 0,nzfft
         xreal(iz) = linechg(iz)
         ximag(iz) = 0.
      enddo

c  Shift the k space origin to the middle
      do iz = 1, nzfft-1,2
         xreal(iz) = -xreal(iz)
      enddo

c  Do the complex fourier transform
      call vcpft(xreal,ximag,nzfft,1,-1,1,1)

c  Save the mode of interest
c     lamkreal(kind=8):: = xreal(indxk)
c     lamkimag = ximag(indxk)
      do iz = 0, nzfft
         lamkreal(iz) = xreal(iz)
         lamkimag(iz) = ximag(iz)
      enddo

      return
      end
c=============================================================================
       subroutine psort(npart,zbeam,dz)
      use InPart
      use Particles,Only: pgroup
      use InMeshrz
      use Sortrz
       integer(ISZ):: npart
       real(kind=8):: zbeam,dz
       integer(ISZ):: n,ic,isgn,ivmin,ivmax,ignme,ignmo,ignold
       integer(ISZ):: ign

c      call zeroarry(npic,nz+1)
       npic = 0
c   Put each particle in a group.  The particles in odd cells go into
c   ``positive'' groups while the particles in even cells go into
c   ``negative'' groups.  nig is number of particles per cell here.
       do n=1,npart
          ic = (pgroup%zp(n) - zbeam - zmmin)/dz
          npic(ic) = npic(ic) + 1
          isgn = mod(ic,2)
          isgn = 2*isgn - 1
          igrp(n) = isgn * npic(ic)
       enddo
c   Find the number of groups, odd and even
c      ignme = abs(min(igrp)) + 1
c      ignmo = max(igrp)

c      call minmx(igrp,1,npart,1,ivmin,ivmax)
       ivmin = 100000
       ivmax = -100000
       do n = 1,npart
           if (igrp(n) .lt. ivmin) ivmin = igrp(n)
           if (igrp(n) .gt. ivmax) ivmax = igrp(n)
       enddo
       ignme = abs(ivmin) + 1
       ignmo = ivmax
c   ignmax is the total number of groups
       ignold = ignmax
       ignmax = ignme + ignmo
c   Make all the group numbers positive
       do n = 1,npart
          igrp(n) = igrp(n) + ignme
       enddo

c   Allocate space if ignmax is larger than it was the last time 
       if (ignmax .gt. ignold) call gchange("Sortrz",0)

c      call zeroarry(nig,ignmax)
       nig = 0

c   Get the number of particles in each group
       do n = 1,npart
          ign = igrp(n)
          nig(ign) = nig(ign) + 1
       enddo

       ipnt(1) = 1
       ipindex(1) = 1
c   Put the particle numbers into an array and use pointers to 
c   keep track of where the groups begin and end
       do ign = 2,ignmax
          ipnt(ign) = ipnt(ign-1) + nig(ign-1)
          ipindex(ign) = ipnt(ign)
       enddo
       ipnt(ignmax+1) = ipnt(ignmax) + nig(ignmax)

       do n = 1, npart
          ipgrp(ipindex(igrp(n))) = n
          ipindex(igrp(n)) = ipindex(igrp(n)) + 1
       enddo

       return
       end
c=============================================================================
      subroutine vstcurr (curr, zbeam, q, wght)
      use InPart
      use InMeshrz
      use Picglbrz
      use Particles,Only: pgroup
      use Sortrz
      real(kind=8):: zbeam,q,wght
      real(kind=8):: curr(0:*)
c  Sets charge density
 
c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------

      integer(ISZ):: ign,m,ip,ioffset,npart
      real(kind=8):: g,dzi,w0,w1
 
      g = wght * q / dz
      dzi = 1. / dz

c     --- vectorized loop to compute indices, weights
      do ign = 1, ignmax
       npart = ipnt(ign+1) - ipnt(ign)
       ioffset = ipnt(ign) - 1
       do m = 1,npart
         ip = ipgrp(m+ioffset)
         index(0,m) =  (pgroup%zp(ip) - zbeam - zmmin) * dzi
         w1 = (pgroup%zp(ip) - zbeam - zmmin) * dzi  - index(0,m)
         w0 = 1. - w1
         s(0,m) = g * w0 * pgroup%gaminv(m) * pgroup%uzp(m)
         s(1,m) = g * w1 * pgroup%gaminv(m) * pgroup%uzp(m)
       enddo

       do m = 1,npart
         if (pgroup%uzp(ipgrp(m+ioffset)) /= 0.) then
         curr(index(0,m))  =curr(index(0,m))   + s(0,m)
         curr(index(0,m)+1)=curr(index(0,m)+1) + s(1,m)
         endif
       enddo
      enddo
 
      return
      end

c=============================================================================
      subroutine vstrhorz (rho1d, zbeam, q, wght, rmesh)
      use Constant
      use InPart
      use InMeshrz
      use Picglbrz
      use Particles,Only: pgroup
      use Sortrz
      real(kind=8):: zbeam,q,wght
      real(kind=8):: rho1d(0:*)
      real(kind=8):: moff(0:3)
c  Sets charge density
 
      real(kind=8):: rmesh(0:nr)
      integer(ISZ):: ign,npart,ioffset,ip,i,k,ind0,m
      real(kind=8):: g,r1,r2,w1,w2

c   Set rmesh(0) to dr/8 so that the first loop vectorizes
c   set it back at the end of the routine
      rmesh(0) = dr/8.
 
c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
 
      g = wght * q / (pi * (dr * dz)**2)
      moff(0) = 0
      moff(1) = 1
      moff(2) = nr+1
      moff(3) = nr+2

c     --- vectorized loop to compute indices, weights
      do ign = 1, ignmax
       npart = ipnt(ign+1) - ipnt(ign)
       ioffset = ipnt(ign) - 1
       do m = 1,npart
         ip = ipgrp(m+ioffset)
         i =  (pgroup%xp(ip) - rmmin) / dr
         r1 = (pgroup%xp(ip) - rmmin)      - i*dr
         r2 = dr - r1
         k =  (pgroup%zp(ip) - zbeam - zmmin) / dz
         w1 = (pgroup%zp(ip) - zbeam - zmmin)     - k*dz
         w2 = dz - w1
         ind0 = i + k*(nr+1)
         index(0,m) = ind0 + moff(0)
         index(1,m) = ind0 + moff(1)
         index(2,m) = ind0 + moff(2)
         index(3,m) = ind0 + moff(3)
         s(0,m) = g * r2 * w2 / (2.*rmesh(i))
         s(2,m) = g * r2 * w1 / (2.*rmesh(i))
         s(1,m) = g * r1 * w2 / (2.*rmesh(i+1)) 
         s(3,m) = g * r1 * w1 / (2.*rmesh(i+1)) 
       enddo
         
       do m = 1,npart
         if (pgroup%uzp(ipgrp(m+ioffset)) /= 0.) then
         rho1d(index(0,m))=rho1d(index(0,m)) + s(0,m)
         rho1d(index(1,m))=rho1d(index(1,m)) + s(1,m)
         rho1d(index(2,m))=rho1d(index(2,m)) + s(2,m)
         rho1d(index(3,m))=rho1d(index(3,m)) + s(3,m)
         endif
       enddo
      enddo
 
c  Return rmesh(0) to the correct value
      rmesh(0) = rmmin

      return
      end

c=============================================================================
      subroutine ldshell
      use Beam_acc
      use Constant
      use InPart
      use InPartrz
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      use InGenrz
      use Particles,Only: pgroup,np

       integer(ISZ):: is,iz
       real(kind=8):: zlen,gfact,z0

       zlen = zimax - zimin
       np = pgroup%npmax
       do is = 1,pgroup%ns
          pgroup%sq(is) = zion * echarge
          pgroup%sm(is) = aion * amu
          if (vbeam .ne. 0.)
     &     pgroup%sw(is) = ibeam * zlen / (vbeam * echarge * zion * np)
       enddo
       gfact = log(rmmax/(a0+dr))/(2.*pi*eps0)
       vphase = sqrt(gfact*pgroup%sq(1)*ibeam/(vbeam*pgroup%sm(1)))

c  Load particles in a shell, one particle per cell
       
       if (vzpert .lt. 0.) then
       z0 = 1./kzpert
       do iz = 1,nz
          pgroup%zp(iz) = zimin + dz*(iz-1) + dz/2.
          pgroup%xp(iz) = a0 + dr/2. 
          pgroup%uzp(iz) = vbeam + vzpert*vphase/(1. + (pgroup%zp(iz)/z0)**2)
          pgroup%uxp(iz) = 0.
          pgroup%uyp(iz) = 0.
          pgroup%yp(iz) = 0.
       enddo
       elseif (vzpert .gt. 0.) then
       do iz = 1,nz
          pgroup%zp(iz) = zimin + dz*(iz-1) + dz/2.
          pgroup%xp(iz) = a0 + dr/2. 
          pgroup%uzp(iz) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(iz))
          pgroup%uxp(iz) = 0.
          pgroup%uyp(iz) = 0.
          pgroup%yp(iz) = 0.
       enddo

       do iz = 1,nz
          pgroup%zp(iz) = pgroup%zp(iz)  + vzpert * sin(kzpert*pgroup%zp(iz))/(kzpert)
       enddo
       else
       do iz = 1,nz
          pgroup%zp(iz) = zimin + dz*(iz-1) + dz/2.
          pgroup%xp(iz) = a0 + dr/2. 
          pgroup%uzp(iz) = vbeam 
          pgroup%uxp(iz) = 0.
          pgroup%uyp(iz) = 0.
          pgroup%yp(iz) = 0.
       enddo
       endif


       pgroup%nps(1) = nz

       return
       end
c=============================================================================
      subroutine ldlinbm
      use Beam_acc
      use Constant
      use InPart
      use InPartrz
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      use InGenrz
      use Particles,Only: pgroup,np


       integer(ISZ):: izso2,izendo2,ia,ipart,iz,isame,ir,ip,is
       real(kind=8):: zlen,zs,zends,gfact,z0,z02,zmid,zt,aend

       zlen = zimax - zimin

c  Load particles 

       zs = zlen * straight
       zends = zlen*(1.-straight)
       izso2 = int(zs/(2.*dz)) -1
       izendo2 = int(zends/(2.*dz))
       ia = int(a0/dr)
       ipart = 0
       gfact = log(rmmax/((ia+1)*dr))/(2.*pi*eps0)
       if (vbeam .ne. 0.) then
       vphase = sqrt(zion*echarge*gfact*ibeam/(vbeam*aion*amu))
       endif
       if (kzpert .ne. 0.) then
       z0 = 1./kzpert
c      z02 = (1./kzpert)**2
       endif
c      zmid = 0.5*(zimax + zimin)

c  Load the particles in the straight section
       do iz = nz/2, nz/2 + izso2
          do isame = 1,5
            ipart = ipart + 1
            pgroup%zp(ipart) = (iz + 0.5)*dz + zmmin
c           pgroup%xp(ipart) = 0.
            pgroup%xp(ipart) = .8 * dr
            pgroup%yp(ipart) = 0.
            pgroup%uxp(ipart) = 0.
            pgroup%uyp(ipart) = 0.
            if (vzpert .lt. 0.) then
              pgroup%uzp(ipart) = vbeam + vzpert*vphase/(1. + (pgroup%zp(ipart)/z0)**2)
c             pgroup%uzp(ipart) = vbeam + vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)
c             pgroup%uxp(ipart) = vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)*
c    &                     (pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
              pgroup%uxp(ipart) = vzpert*vphase*pgroup%xp(ipart)*pgroup%zp(ipart)/
     &                     (z0 *(1. + (pgroup%zp(ipart)/z0)**2))**2
            else
            pgroup%uzp(ipart) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(ipart))
            pgroup%uxp(ipart) = pgroup%xp(ipart)*vzpert*vphase*kzpert*
     &                   0.5*sin(kzpert*pgroup%zp(ipart))
            endif
          enddo
       do ir = 1,ia-1
c        do isame = 1,8*ir
         do isame = 1,(8*ir+4)
           ipart = ipart + 1
           pgroup%zp(ipart) = (iz + 0.5)*dz + zmmin
c          pgroup%xp(ipart) = ir*dr 
           pgroup%xp(ipart) = dr * (ir + 0.5*(1. + (1./(2.*ir+1.)))) 
           pgroup%yp(ipart) = 0.
           pgroup%uxp(ipart) = 0.
           pgroup%uyp(ipart) = 0.
           if (vzpert .lt. 0.) then
             pgroup%uzp(ipart) = vbeam + vzpert*vphase/(1. + (pgroup%zp(ipart)/z0)**2)
c          pgroup%uzp(ipart) = vbeam + vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)
c          pgroup%uxp(ipart) = vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)*
c    &                  (pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
           pgroup%uxp(ipart) = vzpert*vphase*pgroup%xp(ipart)*pgroup%zp(ipart)/
     &                  (z0 *(1. + (pgroup%zp(ipart)/z0)**2))**2
           else
           pgroup%uzp(ipart) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(ipart))
           pgroup%uxp(ipart) = pgroup%xp(ipart)*vzpert*vphase*kzpert*
     &                  0.5*sin(kzpert*pgroup%zp(ipart))
           endif
         enddo
       enddo
         do isame = 1,12*ia+8
           ipart = ipart + 1
           pgroup%zp(ipart) = (iz + 0.5)*dz + zmmin
c          pgroup%xp(ipart) = ir*dr 
           pgroup%xp(ipart) = dr * (ia + (8.*ia+8.)/(12.*ia+8.))
           pgroup%yp(ipart) = 0.
           pgroup%uxp(ipart) = 0.
           pgroup%uyp(ipart) = 0.
           if (vzpert .lt. 0.) then
             pgroup%uzp(ipart) = vbeam + vzpert*vphase/(1. + (pgroup%zp(ipart)/z0)**2)
c          pgroup%uzp(ipart) = vbeam + vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)
c          pgroup%uxp(ipart) = vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)*
c    &                  (pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
           pgroup%uxp(ipart) = vzpert*vphase*pgroup%xp(ipart)*pgroup%zp(ipart)/
     &                  (z0 *(1. + (pgroup%zp(ipart)/z0)**2))**2
           else
           pgroup%uzp(ipart) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(ipart))
           pgroup%uxp(ipart) = pgroup%xp(ipart)*vzpert*vphase*kzpert*
     &                  0.5*sin(kzpert*pgroup%zp(ipart))
           endif
         enddo

       enddo
c  Load the particles in the tapered ends
      if (straight .ne. 1.) then
       do iz = nz/2+izso2+1 , nz/2+izso2+izendo2
         zt = (iz + 0.5)*dz + zmmin
         aend = a0*(1. - (zt - 0.5*zs)**2/(0.5*zlen - 0.5*zs)**2)
         ia = int(aend/dr)
         do isame = 1,5
           ipart = ipart + 1
           pgroup%zp(ipart) = zt
c          pgroup%xp(ipart) = 0.
           pgroup%xp(ipart) = .8 * dr
           pgroup%yp(ipart) = 0.
           pgroup%uxp(ipart) = 0.
           pgroup%uyp(ipart) = 0.
           if (vzpert .lt. 0.) then
             pgroup%uzp(ipart) = vbeam + vzpert*vphase/(1. + (pgroup%zp(ipart)/z0)**2)
c            pgroup%uzp(ipart) = vbeam + vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)
c            pgroup%uxp(ipart) = pgroup%uzp(ipart)*(pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
c            pgroup%uxp(ipart) = vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)*
c    &                    (pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
             pgroup%uxp(ipart) = vzpert*vphase*pgroup%xp(ipart)*pgroup%zp(ipart)/
     &                    (z0 *(1. + (pgroup%zp(ipart)/z0)**2))**2
           else
             pgroup%uzp(ipart) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(ipart))
             pgroup%uxp(ipart) = pgroup%xp(ipart)*vzpert*vphase*kzpert*
     &                    0.5*sin(kzpert*pgroup%zp(ipart))
           endif
         enddo
         do ir = 1, ia-1
c          do isame = 1,8*ir
           do isame = 1,(8*ir+4)
             ipart = ipart + 1
             pgroup%zp(ipart) = zt
c            pgroup%xp(ipart) = ir*dr 
             pgroup%xp(ipart) = dr * (ir + 0.5*(1. + (1./(2.*ir+1.)))) 
             pgroup%yp(ipart) = 0.
             pgroup%uxp(ipart) = 0.
             pgroup%uyp(ipart) = 0.
             if (vzpert .lt. 0.) then
               pgroup%uzp(ipart) = vbeam + vzpert*vphase/(1. + (pgroup%zp(ipart)/z0)**2)
c              pgroup%uzp(ipart) = vbeam + vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)
c              pgroup%uxp(ipart) = pgroup%uzp(ipart)*(pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
c              pgroup%uxp(ipart) = vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)*
c    &                      (pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
               pgroup%uxp(ipart) = vzpert*vphase*pgroup%xp(ipart)*pgroup%zp(ipart)/
     &                      (z0 *(1. + (pgroup%zp(ipart)/z0)**2))**2
             else
               pgroup%uzp(ipart) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(ipart))
             endif
           enddo
         enddo
         do isame = 1,12*ia+8
           ipart = ipart + 1
           pgroup%zp(ipart) = (iz + 0.5)*dz + zmmin
c          pgroup%xp(ipart) = ir*dr 
           pgroup%xp(ipart) = dr * (ir + (8.*ia+8.)/(12.*ia+8.))
           pgroup%yp(ipart) = 0.
           pgroup%uxp(ipart) = 0.
           pgroup%uyp(ipart) = 0.
           if (vzpert .lt. 0.) then
             pgroup%uzp(ipart) = vbeam + vzpert*vphase/(1. + (pgroup%zp(ipart)/z0)**2)
c            pgroup%uzp(ipart) = vbeam + vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)
c            pgroup%uxp(ipart) = pgroup%uzp(ipart)*(pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
c            pgroup%uxp(ipart) = vzpert*vphase*exp(-(pgroup%zp(ipart)-zmid)**2/z02)*
c    &                    (pgroup%zp(ipart)-zmid)/z02 * xp(ipart)
             pgroup%uxp(ipart) = vzpert*vphase*pgroup%xp(ipart)*pgroup%zp(ipart)/
     &                    (z0 *(1. + (pgroup%zp(ipart)/z0)**2))**2
           else
             pgroup%uzp(ipart) = vbeam + vzpert*vphase*cos(kzpert*pgroup%zp(ipart))
c            pgroup%uxp(ipart) = xp(ipart)*vzpert*vphase*kzpert*
c    &                    0.5*sin(kzpert*pgroup%zp(ipart))
             pgroup%uxp(ipart) = 0.
           endif
         enddo
       enddo
      endif

      do ip = ipart+1, 2*ipart
        pgroup%zp(ip) = -pgroup%zp(ip-ipart)
        pgroup%xp(ip) = pgroup%xp(ip-ipart)
        pgroup%yp(ip) = 0.
        pgroup%uxp(ip) = -pgroup%uxp(ip-ipart)
c       pgroup%uxp(ip) = 0.
        pgroup%uyp(ip) = 0.
        pgroup%uzp(ip) = pgroup%uzp(ip-ipart)
      enddo
      pgroup%npmax = 2*ipart
      np = pgroup%npmax
      pgroup%nps(1) = pgroup%npmax

      do is = 1,ns
         pgroup%sq(is) = zion * echarge
         pgroup%sm(is) = aion * amu
         if (vbeam .ne. 0.)
     &     pgroup%sw(is) = ibeam * zlen / (vbeam * echarge * zion * np)
      enddo
      return
      end
c--------------------------------------------------------------------
      subroutine fldrz(nstep)
      use InPart
      use InMeshrz
      use InFluidrz
      use Fluidsrz
      use Picglbrz
      use Fieldsrz
      use InGen
      use InGenrz
      use Beam_acc
      use Particles,Only: pgroup
      integer(ISZ):: nstep
      call fluidrz(dt,nstep,vbeam,ibeam,pgroup%sq(1),pgroup%sm(1),a0)
      nfrstrz = nstep + nfrstrz
      return
      end
c--------------------------------------------------------------------
      subroutine fluidrz(dt,nstep,vbeam,ibeam,q,m,a0)
      use Constant
      use InMeshrz
      use Picglbrz
      use Fieldsrz
      use InGenrz
      use InFluidrz
      use Fluidsrz
      integer(ISZ):: nstep
      real(kind=8):: dt,vbeam,ibeam,q,m,a0
      integer(ISZ):: iz,istp,ikz
      real(kind=8):: lamb0,gfact,eoverm,taumdt,taupdt,tauodt,growthrate
      
c   This routine solves the 1-d fluid equations with resistance and
c   capacitance.  It uses the equations found in the article
c   "Asymptotic Analysis of the Longitudinal Instability of a Heavy
c   Ion Induction Linac", E.P. Lee and L. Smith, Proceedings of the 1990
c   Linear Accelerator Conference,  Albuquerque, NM Sept 1990
c   With the addition of the self force (-g d lambda/dz).

c   Fourier transform the initial lambda, vz, and Ez in space

      do iz = 0,nz
         arryr(iz) = real(lamfld(iz))
         arryi(iz) = aimag(lamfld(iz))
      enddo
      call forward(arryr,arryi,nz,dz)
      do iz = 0,nz
         lamfld(iz) = cmplx(arryr(iz),arryi(iz))
      enddo

      do iz = 0,nz
         arryr(iz) = real(vzfld(iz))
         arryi(iz) = aimag(vzfld(iz))
      enddo
      call forward(arryr,arryi,nz,dz)
      do iz = 0,nz
         vzfld(iz) = cmplx(arryr(iz),arryi(iz))
      enddo

      do iz = 0,nz
         arryr(iz) = real(e1fld(iz))
         arryi(iz) = aimag(e1fld(iz))
      enddo
      call forward(arryr,arryi,nz,dz)
      do iz = 0,nz
         e1fld(iz) = cmplx(arryr(iz),arryi(iz))
      enddo
c  set up the constants that are needed
      gfact = log(rmmax/a0)/(2.*pi*eps0)
      lamb0 = ibeam/vbeam
      eoverm = q/m
c     taumdt = 2.*taurc - dt
c     taupdt = 2.*taurc + dt
      tauodt = taurc/dt
      vphase = sqrt(eoverm*gfact*lamb0)

c  Set this up if this is the first time through
      if (nfrstrz .eq. 0) then
      do iz=0,nz-1
         ikdt(iz) = cmplx(0.,sqrt(kzsq(iz)/(dr*dr))*dt)
      enddo
      growthrate = eta*vbeam*vphase/(2.*gfact)
      write(STDOUT,*)"Linear growth rate = ",growthrate

      thistrz(0) = 0.
      lamrhist(0) = real(lamfld(ikmode))
      lamihist(0) = aimag(lamfld(ikmode))
      endif

c  step forward in time
      do 999 istp = nfrstrz, nfrstrz + nstep

c  save lambda at the old time
         do ikz = 0,nz-1
           lamold(ikz) = lamfld(ikz)
         enddo
c  update lambda
         do ikz = 0,nz-1
           lamfld(ikz) = lamfld(ikz) - ikdt(ikz)*lamb0*vzfld(ikz)
         enddo
c  update the electric field 
c         do ikz = 0,nz-1
c           e1fld(ikz) = e1fld(ikz) * taumdt/taupdt - (2.*eta*dt)/taupdt
c    &                * (lamb0 * vzfld(ikz) + 0.5 * vbeam *
c    &                (lamfld(ikz) + lamold(ikz)))
c        enddo
         do ikz = 0,nz-1
           e1fld(ikz) = (e1fld(ikz)*tauodt - 
     &                  eta*(vbeam*0.5*(lamfld(ikz) + lamold(ikz))
     &                 + lamb0*vzfld(ikz)))/
     &                  (1.+tauodt)
         enddo
c        do ikz = 0,nz-1
c          e1fld(ikz) = -eta*(vbeam*lamfld(ikz) + lamb0*vzfld(ikz))
c        enddo

c  update the velocity

         do ikz = 0,nz-1
           vzfld(ikz) = vzfld(ikz) - eoverm*gfact*ikdt(ikz)*lamfld(ikz)
     &               + dt * eoverm * e1fld(ikz)
         enddo


         if (mod(istp,5) .eq. 0.) then
            ihist = ihist + 1
            thistrz(ihist) = dt*(istp + nfrstrz)
            lamrhist(ihist) = real(lamfld(ikmode))
            lamihist(ihist) = aimag(lamfld(ikmode))
         endif

 999  continue

      do iz = 0,nz
         arryr(iz) = real(lamfld(iz))
         arryi(iz) = aimag(lamfld(iz))
      enddo
      call backward(arryr,arryi,nz,dz,zmmax,zmmin)
      do iz = 0,nz
         lamfld(iz) = cmplx(arryr(iz),arryi(iz))
      enddo

      do iz = 0,nz
         arryr(iz) = real(vzfld(iz))
         arryi(iz) = aimag(vzfld(iz))
      enddo
      call backward(arryr,arryi,nz,dz,zmmax,zmmin)
      do iz = 0,nz
         vzfld(iz) = cmplx(arryr(iz),arryi(iz))
      enddo

      do iz = 0,nz
         arryr(iz) = real(e1fld(iz))
         arryi(iz) = aimag(e1fld(iz))
      enddo
      call backward(arryr,arryi,nz,dz,zmmax,zmmin)
      do iz = 0,nz
         e1fld(iz) = cmplx(arryr(iz),arryi(iz))
      enddo

      return
      end
c----------------------------------------------------------------------------
      subroutine forward(arryr,arryi,nz,dz)
      integer(ISZ):: nz
      real(kind=8):: dz
      real(kind=8):: arryr(0:nz),arryi(0:nz)
c ----------------------------------------------------------------------------
c  Do the forward transform
c ----------------------------------------------------------------------------
      integer(ISZ):: iz
      real(kind=8):: znorm
c        --- shift k-space origin to middle.
         do 600 iz = 1, nz-1, 2
            arryr(iz) = -arryr(iz)
            arryi(iz) = -arryi(iz)
  600    continue
c        --- Normalize the array
         znorm = dz/2.
         do 601 iz = 0,nz
            arryr(iz) = arryr(iz)*znorm
            arryi(iz) = arryi(iz)*znorm
  601    continue
c        --- do the transforms
         call vcpft (arryr(0), arryi(0), nz, 1,-1,1,1)

      return
      end

      subroutine backward(arryr,arryi,nz,dz,zmmax,zmmin)
      integer(ISZ):: nz
      real(kind=8):: dz,zmmax,zmmin
      real(kind=8):: arryr(0:nz),arryi(0:nz)
      integer(ISZ):: iz
      real(kind=8):: znorm
    
c        --- do the back transforms
         call vcpft  (arryr(0), arryi(0), nz, 1,+1,1,1)
c        --- shift k-space origin back.
         do 900 iz = 1, nz-1, 2
            arryr(iz) = -arryr(iz)
            arryi(iz) = -arryi(iz)
  900 continue
         znorm = 2./(zmmax - zmmin)
         do 901 iz = 0,nz-1
            arryr(iz) = arryr(iz)*znorm
            arryi(iz) = arryi(iz)*znorm
  901    continue
         arryr(nz) = arryr(0)
         arryi(nz) = arryi(0)

      return
      end
c---------------------------------------------------------------
      subroutine setexte(it,dt,dz,vbeam,nr,nz)
      use Fieldsrz
      use InGenrz
      use Z_arrays
      integer(ISZ):: it,nr,nz
      real(kind=8):: dt,dz,vbeam
      integer(ISZ):: iz,nl

      if (it .eq. 0) then
         do iz = 0,nz
         curr0(iz) = curr(iz,nszarr)
         ezext(iz) = eta*curr0(iz)
         enddo
c  Also set up the integral of the current for phi on axis plots
         phiresist(0) = 0.
         do iz = 1,nz
           phiresist(iz) = phiresist(iz-1) + curr(iz,nszarr)*dz*eta
         enddo

      elseif (taurc .ne. 0.) then
c  Advect the external electric field 
      nl = int(vbeam*dt/dz)

      do iz = 0, (nz-nl)
         scrtch(iz,0) = ezext(iz + nl)
      enddo

      do iz = (nz-nl+1),nz
         scrtch(iz,0) = ezext(iz-nz+nl-1)
      enddo

      do iz = 0,nz
         ezext(iz) = scrtch(iz,0)
      enddo

c  Update the external electric field using 
c  taurc/dt * (E^(n+1) - E^n) + E^n = eta*Curr0
         do iz = 0,nz
            ezext(iz) = dt/taurc*(ezext(iz)*(taurc/dt - 1.)
     &                           + eta*curr0(iz))
         enddo

         phiresist(0) = 0.
         do iz = 1,nz
           phiresist(iz) = phiresist(iz-1) + ezext(iz)*dz
         enddo
      endif
      return
      end
