#include "top.h"
c=============================================================================
c@(#) File WXY.M, version $Revision: 3.1 $, $Date: 2001/04/19 23:43:17 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package WXY of code WARP
c  XY electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine wxyinit

c  Called at first reference to package (not nec. a "run" etc.).

      call wxyvers (STDOUT)

      return
      end
c=============================================================================
      subroutine wxyvers (iout)
      use WXYversion
      integer     :: iout
c  Echoes code version, etc. to output files as they're created
      call printpkgversion(iout,"Particle package WXY",verswxy)
      return
      end
c=============================================================================
      subroutine wxygen()
      use Ch_var
      use Constant
      use InGen
      use InGen3d
      use InGenxy
      use InDiag
      use InPart
      use InPart3d
      use InMesh3d
      use Fields3d
      use Efields3d
      use Io
      use Lattice
      use LatticeInternal
      use Particles
      use Particlesxy
      use Picglb
      use Picglb3d
      use OutParams
      use Beam_acc
      use Z_arrays
      use Win_Moments
      use Z_Moments
      use Lab_Moments
      use Moments
      use Damped_eom
      use Hist
      use Pspwork
#ifdef PARALLEL
      use Parallel
#endif

c  Invoked by the GENERATE command, it sets up the problem
c  This routine allots all of the neccesary dynamic arrays, calls the
c  particle loader and does the initial load onto the charge density
c  mesh, initializes arrays for the field solver and sets the mesh arrays,
c  does the initial field solve, and sets up other arrays that are needed.

      integer     :: i,j,k,is,ipmin,ip,iwin

c  Announce that we're starting up

      call remark(" ***  particle simulation package WXY generating")

c  Estimate wall radius, needed for g-factor calc 
c  (rwallfac = 1 is probably NOT a good guess)

      rwall = rwallfac * sqrt( xmmax**2 + xmmin**2 )

c  Calculate derived quantities and species related arrays (both set from 
c  derivqty and an internal call to species).  

      call derivqty

c  Calculate step size if not set by the user.
      if (ds == 0.) ds = vbeam*dt
      if (dt == 0.) dt = ds/vbeam

c  Set default values for the axial grid if using thin slice model.
c  Zmmin and zmmax are chosen so that dz = 1.
      if (.not. lthick) then
        if (zmmin == 0.0 .and. zmmax == 0.0) then
          zmmin = -0.5
          zmmax = +0.5
        endif
        if (zimin == 0.0 .and. zimax == 0.0) then
          zimin = zmmin
          zimax = zmmax
        endif
        if (nz == 2) nz = 1
      endif

#ifdef PARALLEL
c     --- Use routine from w3d to initialize the MPI and divide up
c     --- the problem. nz is set to nslaves so that the grid is divided
c     --- up properly among the processors. zmmin and zmmax are then set so
c     --- that dz = 1.
      if (.not. lthick) then
        nz = nslaves
c       zmmin = -nslaves/2.
c       zmmax = +nslaves/2.
c       zimin = zmmin
c       zimax = zmmax
      endif
      call init_w3d_parallel()
#endif

c  Call the species routine again so that the values of zimin_s and
c  zimax_s are set appropriately.  This resolves a circle where zimin is set
c  equal to zmmin which depends on vbeam which is calculated in derivqty
c  which requires things from species which sets zimin_s equal to zimin
c  which at that point will not be set (if the user didn't set it).

      call species

c  When lvzchang is not true, then change number of iterations that are
c  done in the calculation of dt to 1 (so no iteration is done).
      if (.not. lvzchang) niter_dt = 1

c  When using thick slice model, don't do any iterations
      if (lthick) niter_dt = 1

c  Initialize the cycle counter, time, etc.

      if (nrestart == " ") then
         it = 0
         time = 0.
      endif
      call stepid (it, time, zbeam)

c Print values of input variables, other interesting things to text file
      if (warpout > -1) then
        call edit (warpout, "runid")
        call edit (warpout, "it")
        call edit (warpout, "time")
        call edit (warpout, "InGen")
        call edit (warpout, "InDiag")
        call edit (warpout, "InPart3d")
        call edit (warpout, "InMesh3d")
      endif

c  Create the dynamic arrays in Z_arrays; set the z mesh

      if (nzzarr == 0) nzzarr = nz
      call gchange("Z_arrays", 0)
      if (zzmax == 0.) zzmax = zmmax
      if (zzmin == 0.) zzmin = zmmin
      dzz = (zzmax - zzmin)/nzzarr
      dzzi = 1./dzz
      do k = 0, nzzarr
         zplmesh(k) = zzmin + k*dzz
      enddo

c  Re-size the dynamic arrays for the lattice (scan for true length, first).
c  Also allocate internal lattice arrays.  Note that size is set to zero
c  since all particles will be at the same element.

      call remark(" ---  Resetting lattice array sizes")
      call resetlat
      nzl = 0

c  Create the dynamic arrays for fields, contour plot workspace
c  Also set vbeamfrm, adjust x and ymmin for symmetries and set
c  so that by default izfsmax = nz (for full field solve).

      if (nzfull == 0) nzfull = nz
      nmxy  = max(nx,ny)
      nmxyz = max(nx,ny,nzfull)
      call gallot("Fields3d", 0)
      if (efetch == 3) call gallot("Efields3d",0)
      vbeamfrm = vbeam
      if (l2symtry) then
        ymmin = 0.
      elseif (l4symtry) then
        xmmin = 0.
        ymmin = 0.
      endif
      if (izfsmax == 0) izfsmax = nz

c  Set the value of prwall, radius at which particles are lost
      if (prwall == 0.) prwall = xmmax - dx
      do k=0,nzzarr
        prwallz(k) = prwall
        if (prwallxz(k) == 0.) prwallxz(k) = prwallx
        if (prwallyz(k) == 0.) prwallyz(k) = prwally
        if (prwelips(k) == 0.) prwelips(k) = 1.
      enddo
        
c  Calculate mesh dimensioning quantities 

      dx = (xmmax - xmmin) / nx
      dy = (ymmax - ymmin) / ny
      dz = (zmmax - zmmin) / nz
      do i = 0, nx
         xmesh(i) = i * dx + xmmin
      enddo
      do j = 0, ny
         ymesh(j) = j * dy + ymmin
      enddo
      do k = 0, nz
         zmesh(k) = k * dz + zmmin
      enddo

c Calculate location of axis in mesh, the term dx*1.e-5 acts as fuzz
      ix_axis = int(-xmmin/dx + dx*1.e-5)
      iy_axis = int(-ymmin/dy + dy*1.e-5)
      iz_axis = int(-zmmin/dz + dz*1.e-5)

c  Create the dynamic arrays for particles (set npmax to an estimated 
c  length for now, for those loading schemes that don't actually
c  use a user-set npmax directly)

      call remark(" ---  Allocating space for particles")
      if (xrandom == "grid") npmax = nxstripe*nystripe*nzstripe
      if (xrandom == "fibonacc") npmax = nfibgrps*fibg1

      call alotpart
      npmaxxy = npmax
      call gallot("Particlesxy",0)

c  Load the particles, calculate the charge density

      call remark(" ---  Loading particles")
      cigarld = .false.
      cylinder = .true.
      call stptcl3d
      call setgamma(lrelativ)
      call remark(" ---  Setting charge density")
c     call zeroarry(rho,(nx+1)*(ny+1)*(nz+1))
c     call zeroarry(curr,nzzarr+1)
      rho = 0.
      curr = 0.
      do is=1,ns
c        --- Calculate initial guess of dt for each particle.  This will
c        --- not change if vz does not change.
         if (lthick) then
           do ip=ins(is),ins(is)+nps(is)-1
             dtp(ip) = dt
           enddo
         else
           do ip=ins(is),ins(is)+nps(is)-1
             dtp(ip) = dt*vbeam/uzp(ip)
           enddo
         endif
c        --- Set charge density
         do ipmin = ins(is), ins(is) + nps(is) - 1, NPARPGRP
            ip = min0 (NPARPGRP, ins(is)+nps(is)-ipmin)
            call setrhoxy(rho(0,0,0),ip,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
     &                    uzp(ipmin),gaminv(ipmin),sq(is),sw(is))
         enddo
      enddo
      call fixrhoxy(rho,nx,ny,nz,periinz,lthick)

c  Create the dynamic arrays for the partcle qtys needed for the
c  damped mover.  They are always allocated, but with length 1 if not used.

      npdamp = 1
      if (eomdamp .ne. 0.) npdamp = npmax
      call gallot("Damped_eom", 0)

c  Create the dynamic arrays for "window" moments

      call remark(" ---  Allocating Win_Moments")
      zwindows(1,0) = zmmin
      zwindows(2,0) = zmmax
      nzwind = 0
      do iwin = 1, NWINDOWS
         if (zwindows(1,iwin) .ne. zwindows(2,iwin)) nzwind = nzwind + 1
      enddo
      call gallot("Win_Moments", 0)

c  Create the dynamic arrays for z moments
c  The slice code should only calculate the global particle moments (and
c  not the z moments) since the distribution in z of the particles is
c  artificial, the particles are physically at the same z location.  The
c  thick slice model, though, should still calculate the zmoments.

      call remark(" ---  Allocating Z_Moments")
      if (.not. lthick .and. ifzmmnt > 0) ifzmmnt = 1
      if (nzmmnt == 0) nzmmnt = nz
      call gallot("Z_Moments", 0)
      if (zmmntmax == 0.) zmmntmax = zmmax
      if (zmmntmin == 0.) zmmntmin = zmmin
      dzm = (zmmntmax - zmmntmin)/nzmmnt
      dzmi = 1./dzm
      do k = 0, nzmmnt
         zmntmesh(k) = zmmntmin + k*dzm
      enddo

c  Create the dynamic arrays for lab frame moments

      call remark(" ---  Allocating Lab_Moments")
      nlabwn = 0
      do iwin = 1, MAXNUMLW
         if (zlw(iwin) .ne. 0.) nlabwn = nlabwn + 1
      enddo
      if (nlabwn == 0) iflabwn = 0
      if (itlabwn == 0) itlabwn = max(1,nhist)
      ntlabwn = (zmmax-zmmin)/(vbeam*dt*itlabwn) +.5
      if (iflabwn .ne. 0) call gallot("Lab_Moments", 0)

c  Create the scratch arrays for phase space plots (permanent, for now)
c  and set limits for plots

      call remark(" ---  Allocating scratch space for plots")
      npsplt = min0 (NPARPGRP-1, npmax)
      if (npsplt == 0) npsplt = NPARPGRP-1
      call gallot("Pspwork", 0)
      if (xplmin == 0.) xplmin = -xmmax
      if (xplmax == 0.) xplmax =  xmmax
      if (yplmin == 0.) yplmin = -ymmax
      if (yplmax == 0.) yplmax =  ymmax
      if (zplmin == 0.) zplmin =  zmmin
      if (zplmax == 0.) zplmax =  zmmax

c  Choose which particles will be plotted in "subset" plots

      call psubsets

c  Setup history mechanism

      if (nhist > 0) then
c       --- create the dynamic arrays for history data; set pointer into them
        call remark(" ---  Allocating history arrays")
        if (lenhist == 0) lenhist = min ( nt/nhist + 1, 100)
        call gallot("Hist", 0)
        jhist = -1
      elseif (nhist < 0) then
c       --- call interpreter routine to setup hst package
c       --- setup_hst is in bas.wrp
        call execuser("setup_hst")
      endif

c  Print interesting things to plot file and teletype

      call prntpara(dx,dy,dz,lprntpara)
      call prntpa3d(lprntpara)

c  Initial call to fieldsolver in order to initialize attx, kxsq, etc.

      call vpxy (1)

c  Initial fieldsolve, diagnostics

      call stepxy ("wxygen")

      return
      end
c=============================================================================
      subroutine wxyexe()
      use Picglb
      use InPart
      use InGen
      use InGenxy
      use Picglb3d
      use Ctl_to_pic
      use Particles
      use Particlesxy

c  Takes a time step, returns status OK, ERR, or DONE
c  This routine advances the mesh in the lab frame, sets the logicals
c  which control how this next step is to be done, and then calls
c  the routine STEPXY to do the step.

      real(kind=8):: vbeamfrm0,zbeam0
      real(kind=8):: ds0,zz,zcorrection
      integer     :: is,ip,i
      logical     :: lbendend

c  Announce that we're running

      if (it == 0) call remark(" ***  particle simulation package WXY running")

c  Get distance to get to start or end of bend, which ever is closer. If
c  that distance is less than ds, then change ds and dt and to take a
c  fractional timestep that will put the slice at exactly the start or
c  end of the bend.  The timestep size for each particle is also scaled.
c  The factor of (1.-1.e-9) multiplying 'ds' is there to prevent too small
c  of a step.  This also helps eliminate problems with round off when zz
c  is equal to ds.  The 1.e-9 probably should be a variable and changeable,
c  but that value wouldn't be changed much.

      call nextbend(zbeam,zz)
      if (zz < ds*(1.-1.e-9)) then
        lbendend = .true.
      else
        lbendend = .false.
      endif
      if (lbendend) then
        ds0 = ds
        ds = zz
        dt = dt*zz/ds0
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            dtp(ip) = dtp(ip)*dt/(ds0/vbeamfrm)
          enddo
        enddo
        ldiag = .false.
      endif

c  accelerate grid frame and rescale the time step size.
c  The iteration is done to get the correct dt.
c  Note that dt is only changed if the grid frame velocity was changed.
      zbeam0 = zbeam
      vbeamfrm0 = vbeamfrm
      do i=1,niter_dt
        zbeam = zbeam0
        vbeamfrm = vbeamfrm0
        call acclbfrm(zcorrection)
        dt = dt*ds/(zbeam + vbeamfrm*dt + zcorrection - zbeam0)
      enddo
      zbeam = zbeam0

c  set timestep counter, time, and advance grid frame

      it = it + 1
      time = time + dt
      if (lgridqnt) then
c       --- zgrid is integer number of dz's
        zgridprv = int(zbeam/dz + .5)*dz
        zgrid = int((zbeam + dt*vbeamfrm + zcorrection)/dz + .5)*dz
      else
c       --- zgrid is the same as the beam frame
        zgridprv = zbeam
        zgrid = zbeam + dt*vbeamfrm + zcorrection
      endif
      call stepid (it, time, zgrid)

c  set logicals

      lfirst = .false.
      if (ncall == 1) lfirst = .true.
      llast = .false.
      if (ncall == maxcalls) llast = .true.

c  Take the first fractional timestep to put the slice exactly at the edge
c  of the bend.  The next step will then be another fractional timestep so
c  that the sum of the two will push the beam a full ds.  The timestep
c  size is rescaled appropriately for that next step.

      if (lbendend) then
        call stepxy ("wxyexe")

        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            dtp(ip) = dtp(ip)*(ds0/vbeamfrm - dt)/dt
          enddo
        enddo
        ds = ds0 - ds
        dt = ds0/vbeamfrm - dt
        ldiag = .true.

c       --- accelerate grid frame and rescale the time step size.
c       ---  The iteration is done to get the correct dt.
c       --- Note that dt is only changed if the grid frame velocity was changed.
        zbeam0 = zbeam
        vbeamfrm0 = vbeamfrm
        do i=1,niter_dt
          zbeam = zbeam0
          vbeamfrm = vbeamfrm0
          call acclbfrm(zcorrection)
          dt = dt*ds/(zbeam + vbeamfrm*dt + zcorrection - zbeam0)
        enddo
        zbeam = zbeam0

c       --- set timestep counter, time, and advance grid frame
        time = time + dt
        if (lgridqnt) then
c         --- zgrid is integer number of dz's
          zgrid = int((zbeam + dt*vbeamfrm + zcorrection)/dz + .5)*dz
        else
c         --- zgrid is the same as the beam frame
          zgrid = zbeam + dt*vbeamfrm + zcorrection
        endif
        call stepid (it, time, zgrid)
      endif

c  call the routine that does the actual work

      call stepxy ("wxyexe")

c  Reset the step size if extra substep was taken.

      if (lbendend) then
        do is=1,ns
          do ip=ins(is),ins(is)+nps(is)-1
            dtp(ip) = dtp(ip)*ds0/vbeamfrm/dt
          enddo
        enddo
        ds = ds0
        dt = ds0/vbeamfrm
      endif

c  Have we reached the end of the run ?

      if ( lfinishd ) then
         call remark("wxyexe: problem completed.")
         return
      elseif (nplive <= 0) then
           call remark(" *** WXYEXE: stopping, nplive = 0")
           return
      else
         return
      endif

      end
c=============================================================================
      subroutine wxyfin()
      use InGen
      use InGen3d
      use InDiag
      use InPart
      use InPart3d
      use InMesh3d
      use Fields3d
      use Io
      use Lattice
      use LatticeInternal
      use Particles
      use Particlesxy
      use Picglb
      use Picglb3d
      use Win_Moments
      use Z_Moments
      use Z_arrays
      use Hist
      use Pspwork

c  Finish up at end of RUN, or on receipt of FIN
c  This routine is never called, at present; history plots are
c  made using a BASIS interpreter script (histplot), and we just
c  end the run.  If we wanted to chain runs so that an output qty
c  might be plotted vs a parameter, this routine might be useful.
c  It would be needed for a non-BASIS version of WARP.
c  For now it serves as a place-holder.

c  print final edits

c  perform diagnostics (unless we just did)

c  Create history plots

c  Make a restart dump (unless we just did, or the user inhibits it)

c  create final printouts

c  release storage

      call gfree ("Fields3d")
      call gfree ("Hist")
      call gfree ("Win_Moments")
      call gfree ("Z_Moments")
      call gfree ("Lab_Moments")
      call gfree ("Moments")
      call gfree ("Lattice")
      call gfree ("LatticeInternal")
      call gfree ("Particles")
      call gfree ("Particlesxy")
      call gfree ("Pspwork")
      call gfree ("Z_arrays")

      return
      end
c=============================================================================
      subroutine stepxy (caller)
      use Constant
      use InGen
      use InGen3d
      use InGenxy
      use InDiag
      use InPart
      use InMesh3d
      use Fields3d
      use Efields3d
      use Io
      use Particles
      use Picglb
      use Picglb3d
      use LatticeInternal
      use Timers
      character*(*):: caller
      logical     :: thisstep,thiszbeam,dolabwn

c  When called by WXYEXE, stepxy advances the system forward in time one 
c  timestep and gathers diagnostics.  When called by WXYGEN, stepxy takes 
c  a step of zero size, to compute fields, and gather diagnostics at start 
c  of run.


      real(kind=8):: zbeaml,zbeamr,timetemp,wtime

c --- Set the internal lattice variables. This is not generally necessary at
c --- this point (it is redundant most of the time, the next call to
c --- setlatt in this subroutine is sufficient). There are cases where
c --- this is required for consistency. Since it is cheap (time wise),
c --- it is better to make sure the data is consistent than to save a
c --- little bit of time. The value of nzl must be checked since other
c --- packages (like W3D) may have reset it. For example, if the
c --- W3D package is generated after the WXY package, nzl will be set to
c --- zero. Switching back to WXY and running step, the internal lattice
c --- would still be setup for the W3D package and so the step would produce
c --- erroneaous results.
      nzl = 0
      call setlatt

c  Main particle advance: x to t.l. it; v to t.l. it-1/2
c  Half-step in v from t.l. it-1   if last step was "special"
c  Full-step in v from t.l. it-3/2 if last step not "special"
c  No step at all if generating.

      if (caller == "wxyexe") then
        call padvncxy ("halfv")
      endif

c  The next two variables are the left and right ends of the range centered
c  about the end of the current time step plus/minus one half a step.
c  The range is used is determining whether diagnostics are done which
c  are based on the z location of the beam frame.  The diagnostics are done
c  on the time step which ends closest to the value given in the controlling
c  arrays.
      zbeaml   = zbeam - 0.5*vbeamfrm*dt
      zbeamr   = zbeam + 0.5*vbeamfrm*dt

c  Set logical flags to determine if "always" or "seldom" phase space 
c  plots, restart dumps, final timesteps, and moment accumulations should 
c  be done at the end of this step.

      lfinishd = (it .ge. nt) .or. (time .ge. tstop*(1.-MACHEPS)) .or.
     &                             (zbeam .ge. zstop)
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
      ldump    = mod(it, itdump) .eq. 0
      llabwn   = (mod(it,itlabwn) .eq. 0) .and. dolabwn()
      lspecial = (lfinishd .or. lalways .or. lseldom .or. ldump .or. lmoments
     &            .or. llabwn .or. llast .or. (it .eq. 0) .or. allspecl)

c  Set the "gap" electric field. 

      call setegap

c  Gather moments used in diagnostics at "special" timesteps only. 
c  Compute line charge density (gtlchg3d) and the axial line charge 
c  (srhoax3d) on 1-d meshes.  Note -- these moments accumulations are
c  done at this phase of the particle advance to allow for the eventual 
c  use of a single array for rho and phi.   

      if (lspecial .and. ldiag) then 
        call gtlchg3d
        call srhoax3d
      endif 

c  Set lattice; this is done just before field solve, and so is
c  relative to ZBEAM in the same way that self-fields are.

      call setlatt
      call setrstar(rstar,nz,dz,zmmin,zgrid)

c  Field-solve for potential 

      timetemp = wtime()
      if (lbeforefs) call execuser("beforefs")
      call fieldsolxy(-1)
      if (lafterfs) call execuser("afterfs")
      fstime = fstime + (wtime() - timetemp)

c  Pre-calculate the self-E if it is needed for sete3d. This is done after
c  the call to afterfs in case some manipulation is done to phi.
      if (efetch == 3)
     &  call getselfe3d(phi,selfe,nx,ny,nz,dx,dy,dz)

c  If a flag was set making this a "special" step,
c  do a half-advance to bring v to t.l. it 

      if (caller == "wxyexe") then
         call padvncxy ("synchv")
      elseif (caller == "wxygen") then
         call padvncxy ("gen")
      endif

c  Gather moments used in diagnostics at "special" timesteps only. 
c  Compute mean beam z velocity from current and line charge density 
c  on a 1-d mesh.  Also, calculate the electrostatic energy (getese3d), 
c  electrostatic potential on axis (sphiax3d), and the axial electric 
c  field (sezax3d).  

      if (lspecial .and. ldiag) then 
        call getvzofz
        call getese3d 
        call sphiax3d 
        call  sezax3d
      endif 

c  1d array plot diagnostics.

      if (ldiag .and. (lalways .or. lseldom)) call onedplts(ALWAYS)
      if (ldiag .and. (lseldom))              call onedplts(SELDOM)

c  Phase space diagnostics

      if (ldiag .and. (lalways .or. lseldom)) call psplots (ALWAYS)
      if (ldiag .and. (lseldom))              call psplots (SELDOM)

c  Finally, moment diagnostic printout and history storage

      if (ldiag .and. (caller == "wxygen" .or. lspecial))
     &  call minidiag (it,time,lmoments)

      return
      end
c=============================================================================
      subroutine extebxy(np,xp,yp,zp,uzp,gaminv,dtl,dtr,
     &           dodec,bz0,bx,by,bz,ex,ey,ez,m,q,bendres,bendradi,lexbend,
     &           gammabar,zbeam,vbeam,dt,time)
      use Lattice
      use LatticeInternal
      use Mult_data
      use BGRDdata
      use PGRDdata
      use Timers
      integer     :: np
      real(kind=8):: dtl,dtr,dodec,bz0,m,q,gammabar,zbeam,vbeam,dt,time
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uzp(np), gaminv(np)
      real(kind=8):: bx(np), by(np), bz(np),bendres,bendradi
      real(kind=8):: ex(np), ey(np), ez(np)
      logical     :: lexbend


c  Calculates "external" E, B fields
c  Calculates electric or magnetic AG focusing fields, bending and dipole
c  fields, and accelerating fields.
c  Includes back-rotation associated with coordinate transformation into By. 

c  NOTE: When we (someday) set B_self from a Lorentz transformation
c  on E_self, we'll have to carefully work out a sequence of calls,
c  since this routine is called more than once in the PADVNCXY loop on
c  a single step at present.  Perhaps we will also have to compute
c  B_self more than once.

      integer     :: j,ip,in,iele,ii 
      real(kind=8):: vz,zl,zr,fl,fr,frac
      real(kind=8):: xpmqoff,ypmqoff,qi
      real(kind=8):: vzi,dti
      real(kind=8):: n,v,f,fz,rpow,cosnt,sinnt,coeff,alpha,coeffp,alphap
      integer     :: izm(NPARPGRP),imm(NPARPGRP)
      real(kind=8):: xxx(NPARPGRP),yyy(NPARPGRP),zzz(NPARPGRP),ttt(NPARPGRP)
      real(kind=8):: wzm(NPARPGRP),rr(NPARPGRP),tt(NPARPGRP)
      real(kind=8):: mltph(NPARPGRP),mltsf(NPARPGRP)
      real(kind=8):: tex(NPARPGRP),tey(NPARPGRP),tez(NPARPGRP)
      real(kind=8):: tbx(NPARPGRP),tby(NPARPGRP),tbz(NPARPGRP)
      real(kind=8):: txpmqoff(NPARPGRP),typmqoff(NPARPGRP)
      real(kind=8):: trpmqoff(NPARPGRP),ttpmqoff(NPARPGRP)
      real(kind=8):: timetemp,wtime

      timetemp = wtime()

      dti = 1./(dtr-dtl)

c     --- For most elements, all particles are considered to have the same
c     --- axial velocity and z position.  The is equivalent to having
c     --- particles with different axial velocity and different timestep
c     --- sizes so that each travel the same axial distance on a time
c     --- step.  The one exception is the accelerating elements, where
c     --- the residence correction depends on the time spent inside the
c     --- element (and not the distance, which would be the same thing
c     --- IF the axial velocity did not change over a time step).
      vz = vbeam
      vzi = 1./dvnz(vz)
c     --- All particles lie at the same z-cell
      j = 0
c     --- "left" end of velocity advance step
      zl = zbeam + vz * dtl
c     --- "right" end of velocity advance step
      zr = zbeam + vz * dtr

c     --- handle quads
      if (quads) then
c        --- "left" end of velocity advance step
         fl = 0.
         if (zl >= cquadzs(j) .and. zl < cquadze(j)) fl = 1.
c        --- "right" end of velocity advance step
         fr = 0.
         if (zr >= cquadzs(j) .and. zr < cquadze(j)) fr = 1.
c        --- residence fraction
         frac = fl
         if (fl > fr) frac = (cquadze(j)-zl) * vzi*dti
         if (fr > fl) frac = (zr-cquadzs(j)) * vzi*dti
c        --- set the field if inside a quad
         if (frac > 0.) then
           do ip = 1, np
              xpmqoff  = xp(ip) - cqoffx(j)
              ypmqoff  = yp(ip) - cqoffy(j)
              bx(ip) = bx(ip) + cquaddb(j) * frac * ypmqoff
              by(ip) = by(ip) + cquaddb(j) * frac * xpmqoff
              ex(ip) = ex(ip) + cquadde(j) * frac * xpmqoff
              ey(ip) = ey(ip) - cquadde(j) * frac * ypmqoff
              if (dodec .ne. 0.) then
                ex(ip) = ex(ip) + dodec*cquadde(j)*frac*
     &                   (6.*xpmqoff**5 - 60.*xpmqoff**3*ypmqoff**2 + 
     &                   30.*xpmqoff*ypmqoff**4)
                ey(ip) = ey(ip) - dodec*cquadde(j)*frac*
     &                   (6.*ypmqoff**5 - 60.*ypmqoff**3*xpmqoff**2 + 
     &                   30.*ypmqoff*xpmqoff**4)
              endif
           enddo
         endif
      endif

c     --- handle dipos 
      if (dipos) then
c        --- find points of entry/exit for "nearest" dipole
c        --- "left" end of velocity advance step
         fl = 0.
         if (zl >= cdipozs(j) .and. zl < cdipoze(j)) fl = 1.
c        --- "right" end of velocity advance step
         fr = 0.
         if (zr >= cdipozs(j) .and. zr < cdipoze(j)) fr = 1.
c        --- set fraction of v-advance step actually within dipole
         frac = fl
         if (fl > fr) frac = (cdipoze(j)-zl) * vzi*dti
         if (fr > fl) frac = (zr-cdipozs(j)) * vzi*dti
c        --- set the field
         if (frac > 0.) then
           do ip = 1, np
              by(ip) = by(ip) + cdipoby(j) * frac 
              bx(ip) = bx(ip) + cdipobx(j) * frac 
              ex(ip) = ex(ip) + cdipoex(j) * frac 
              ey(ip) = ey(ip) + cdipoey(j) * frac 
           enddo
         endif
      endif

c     --- handle sexts
      if (sexts) then
c        --- "left" end of velocity advance step
         fl = 0.
         if (zl >= csextzs(j) .and. zl < csextze(j)) fl = 1.
c        --- "right" end of velocity advance step
         fr = 0.
         if (zr >= csextzs(j) .and. zr < csextze(j)) fr = 1.
c        --- residence fraction
         frac = fl
         if (fl > fr) frac = (csextze(j)-zl) * vzi*dti
         if (fr > fl) frac = (zr-csextzs(j)) * vzi*dti
c        --- set the field
         if (frac > 0.) then
           do ip = 1, np
              xpmqoff  = xp(ip)
              ypmqoff  = yp(ip)
              bx(ip) = bx(ip) + csextdb(j)*frac*3.*(xpmqoff**2 - ypmqoff**2)
              by(ip) = by(ip) + csextdb(j)*frac*(-6.)*xpmqoff*ypmqoff
              ex(ip) = ex(ip) + csextde(j)*frac*3.*(xpmqoff**2 - ypmqoff**2)
              ey(ip) = ey(ip) + csextde(j)*frac*(-6.)*xpmqoff*ypmqoff
           enddo
         endif
      endif

c     --- handle hard-edge electric and magnetic multipoles
      if (heles) then
c       --- calculate coordinates, element indices, and residence
c       --- fractions over particle block and store in temp arrays
c       --- These factors will be used repeatedly when the fields
c       --- are accumulated at each particle.
c       --- find element index
        iele = cheleid(j)
c       --- "left" end of velocity advance step
        fl = 0.
        if (zl >= chelezs(j) .and. zl < cheleze(j)) fl = 1.
c       --- "right" end of velocity advance step
        fr = 0.
        if (zr >= chelezs(j) .and. zr < cheleze(j)) fr = 1.
c       --- residence fraction
        frac = fl
        if (fl > fr) frac = (cheleze(j)-zl) * vzi*dti
        if (fr > fl) frac = (zr-chelezs(j)) * vzi*dti
        if (frac > 0.) then
          do ip = 1, np
c           --- x,y coordinates about multipole center
            xpmqoff = xp(ip) - heleox(iele)
            ypmqoff = yp(ip) - heleoy(iele)
            txpmqoff(ip) = xpmqoff
            typmqoff(ip) = ypmqoff
c           --- r,theta coordinates about multipole center
            trpmqoff(ip) = sqrt(xpmqoff**2 + ypmqoff**2)
            ttpmqoff(ip) = atan2(ypmqoff,dvnz(xpmqoff))
          enddo
c         --- accumulate field contributions at each particle in
c         --- block looping (outermost) over each multipole component
c         --- for vectorization
          do ii = 1, nhmlt
c           --- loop over particle block using stored coordinates, element
c           --- indices, and residence fractions
            do ip = 1,np
              n = hele_n(ii,iele)
              v = hele_v(ii,iele)
c             --- electric multipoles
              if (ii <= helene(iele)) then
                if (n == 0) then
c                 --- accelerating field
c                 --- Note that this really should not be used since the
c                 --- accl element does a better job.  The transverse
c                 --- fields should only be kicks at the entrance and exit.
                  rpow = dvnz(trpmqoff(ip))**(2*v)
                  f = 1./(2.*v+2.)
                  ex(ip) = ex(ip) - frac*heleep(ii,iele)*rpow*f*txpmqoff(ip)
                  ey(ip) = ey(ip) - frac*heleep(ii,iele)*rpow*f*typmqoff(ip)
                  ez(ip) = ez(ip) + frac*heleae(ii,iele)*rpow
                else
c                 --- rest of components
                  cosnt = cos(n*(ttpmqoff(ip)-helepe(ii,iele)))
                  sinnt = sin(n*(ttpmqoff(ip)-helepe(ii,iele)))
                  rpow = dvnz(trpmqoff(ip))**(n-2+2*v)
                  f = 1. + 2.*v/n
                  ex(ip) = ex(ip) + frac*heleae(ii,iele)*rpow*
     &                           (f*txpmqoff(ip)*cosnt + typmqoff(ip)*sinnt)
                  ey(ip) = ey(ip) + frac*heleae(ii,iele)*rpow*
     &                           (f*typmqoff(ip)*cosnt - txpmqoff(ip)*sinnt)
                endif
              endif
c             --- magnetic multipoles
              if (ii <= helenm(iele)) then
                if (n == 0) then
c                 --- Solenoid field
c                 --- The transverse fields should only be kicks at the entrance
c                 --- and exit.
                  rpow = dvnz(trpmqoff(ip))**(2*v)
                  f = 1./(2.*v+2.)
                  bx(ip) = bx(ip) - frac*helemp(ii,iele)*rpow*f*txpmqoff(ip)
                  by(ip) = by(ip) - frac*helemp(ii,iele)*rpow*f*typmqoff(ip)
                  bz(ip) = bz(ip) + frac*heleam(ii,iele)*rpow
                else
c                 --- rest of components
                  cosnt = cos(n*(ttpmqoff(ip)-helepm(ii,iele)))
                  sinnt = sin(n*(ttpmqoff(ip)-helepm(ii,iele)))
                  rpow = dvnz(trpmqoff(ip))**(n-2+2*v)
                  f = 1. + 2.*v/n
                  bx(ip) = bx(ip) + frac*heleam(ii,iele)*rpow*
     &                           (f*txpmqoff(ip)*sinnt - typmqoff(ip)*cosnt)
                  by(ip) = by(ip) + frac*heleam(ii,iele)*rpow*
     &                           (f*typmqoff(ip)*sinnt + txpmqoff(ip)*cosnt)
                endif
              endif
            enddo
          enddo
        endif
      endif

c     --- fold in coordinate transformation associated with bends
      if (bends .and. linbend .and. .not. lexbend) then
         qi = 1./q
         do ip = 1,np
            by(ip) = by(ip) - bendres*(m*qi)*uzp(ip)/(bendradi+xp(ip))
         enddo
      endif


c  ------------------------------------------------------------
c  End of quads, dipos, and sexts
c  ------------------------------------------------------------

c Apply axial field

      if (bz0 .ne. 0.) then
         do ip = 1,np
            bz(ip) = bz(ip) + bz0
         enddo
      endif

c  ------------------------------------------------------------
c  Apply electrostatic multipole components
c  ------------------------------------------------------------
c  For the XY code, the only change should be in mltlocatxy, so that it
c  uses the same axial position (zbeam) for all particles.  Note that
c  for fields near a bend, the axial position relative to the center of
c  an element will vary if the slice is at an angle (in a bend) relative
c  to the element axis.

      if (emlts .and. linemlt) then
c  Apply electrostatic multipoles.
c  See HIF-note 96-10, equation D17, with E replacing B and
c  -n*psi replacing psi.
c  Field is of the following form, where n and v are integers, n>0, v>=0
c    Er = - sum_v  [ E'0v/(2*v+2)*r**(2*v+1) ] +
c         sum_nv [ Env*(1+2*v/n)*r**(n-1+2*v)*cos(n*(theta-psi)) ]
c    Et = - sum_nv [ Env*r**(n-1+2*v)*sin(n*(theta-psi)) ]
c    Ez = sum_v  [ E0v*r**(2*v) ] +
c         sum_nv [ E'nv/n*r**(2*v+n)*cos(n*(theta-psi)) ]
c
c    Ex = Er*cos(theta) - Et*sin(theta)
c    Ey = Er*sin(theta) + Et*cos(theta)
c
c  Note that the total phase angle has two parts, one associated with the
c  multipole data, and one associated with the lattice element.  This allows
c  an element with multipoles with different phases to be physically rotated
c  by changing only one variable, the phase associated with the element.
c  This allows the same set of data to be used by lattice elements with
c  different errors in the angle.
c  An additional scale factor is included in the field, allowing different
c  elements to use the same data set scaled to different values.  The scale
c  factor is added to one (in mltlocatxy) and multiplies the applied field.
c  The addition is done in mltlocatxy so that it is only done once per particle
c  and not once per particle per multipole component.

c       --- Get location of particle relative to moment data.
        call mltlocatxy(np,xp,yp,zp,nemltsets,dzemlt,nzemlt,
     &                  nzl,cemltzs,cemltze,cemltph,cemltsf,cemltsc,
     &                  cemltim,cemltox,cemltoy,
     &                  tex,tey,tez,imm,izm,wzm,xxx,yyy,zzz,ttt,rr,tt,
     &                  mltph,mltsf,
     &                  zlmin,zbeam,dzli,bends,linbend,cbendzs,cbendze,cbendrc)

c       --- accumulate the E field of the electrostatic multipoles
        do in=1,nesmult
          n = emlt_n(in)
          v = emlt_v(in)
          if (n == 0) then
c           --- Apply accelerating field and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                     (esemltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = dvnz(rr(ip))**(2*v)
                tex(ip) = tex(ip) - coeffp*rpow*f*xxx(ip)
                tey(ip) = tey(ip) - coeffp*rpow*f*yyy(ip)
                tez(ip) = tez(ip) + coeff*rpow
              endif
            enddo
          else
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      esemltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      esemltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (esemltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (esemlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      esemlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (esemltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      esemltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = dvnz(rr(ip))**(n-2+2*v)
                tex(ip)=tex(ip) + coeff*rpow*(f*xxx(ip)*cosnt + yyy(ip)*sinnt)
                tey(ip)=tey(ip) + coeff*rpow*(f*yyy(ip)*cosnt - xxx(ip)*sinnt)
                tez(ip)=tez(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*cosnt - alphap*coeff*sinnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            ex(ip) = ex(ip) + tex(ip)*cos(ttt(ip)) - tez(ip)*sin(ttt(ip))
            ey(ip) = ey(ip) + tey(ip)
            ez(ip) = ez(ip) + tex(ip)*sin(ttt(ip)) + tez(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            ex(ip) = ex(ip) + tex(ip)
            ey(ip) = ey(ip) + tey(ip)
            ez(ip) = ez(ip) + tez(ip)
          enddo
        endif
      endif

c  ------------------------------------------------------------
c  Apply magnetostatic multipole components
c  ------------------------------------------------------------
c  (See comments above for emlt elements in XY code.)

      if (mmlts .and. linmmlt) then
c  Apply magnetostatic multipoles.
c  See HIF-note 96-10, equation D17, with -n*psi-pi/2 replacing psi.
c  Field is of the following form, where n and v are integers, n>0, v>=0
c    Br = - sum_v  [ B'0v/(2*v+2)*r**(2*v+1) ] +
c         sum_nv [ Bnv*(1+2*v/n)*r**(n-1+2*v)*sin(n*(theta-psi)) ]
c    Bt = sum_nv [ Bnv*r**(n-1+2*v)*cos(n*(theta-psi)) ]
c    Bz = sum_v  [ B0v*r**(2*v) ] +
c         sum_nv [ B'nv/n*r**(2*v)*sin(n*(theta-psi)) ]
c
c    Bx = Br*cos(theta) - Bt*sin(theta)
c    By = Br*sin(theta) + Bt*cos(theta)
c
c  Note that the total phase angle has two parts, one associated with the
c  multipole data, and one associated with the lattice element.  This allows
c  an element with multipoles with different phases to be physically rotated
c  by changing only one variable, the phase associated with the element.
c  This allows the same set of data to be used by lattice elements with
c  different errors in the angle.
c  An additional scale factor is included in the field, allowing different
c  elements to use the same data set scaled to different values.  The scale
c  factor is added to one (in mltlocatxy) and multiplies the applied field.
c  The addition is done in mltlocatxy so that it is only done once per particle
c  and not once per particle per multipole component.

c       --- Get location of particle relative to moment data.
        call mltlocatxy(np,xp,yp,zp,nmmltsets,dzmmlt,nzmmlt,
     &                  nzl,cmmltzs,cmmltze,cmmltph,cmmltsf,cmmltsc,
     &                  cmmltim,cmmltox,cmmltoy,
     &                  tbx,tby,tbz,imm,izm,wzm,xxx,yyy,zzz,ttt,rr,tt,
     &                  mltph,mltsf,
     &                  zlmin,zbeam,dzli,bends,linbend,cbendzs,cbendze,cbendrc)

c       --- accumulate the B field of the magnetostatic multipoles
        do in=1,nmsmult
          n = mmlt_n(in)
          v = mmlt_v(in)
          if (n == 0) then
c           --- Apply solenoidal field component and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                    (msmmltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                     msmmltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                    (msmmlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                     msmmlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = dvnz(rr(ip))**(2*v)
                tbx(ip) = tbx(ip) - coeffp*rpow*f*xxx(ip)
                tby(ip) = tby(ip) - coeffp*rpow*f*yyy(ip)
                tbz(ip) = tbz(ip) + coeff*rpow
              endif
            enddo
          else
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      msmmltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      msmmltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (msmmltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (msmmlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      msmmlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (msmmltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = dvnz(rr(ip))**(n-2+2*v)
                tbx(ip)=tbx(ip) + coeff*rpow*(f*xxx(ip)*sinnt - yyy(ip)*cosnt)
                tby(ip)=tby(ip) + coeff*rpow*(f*yyy(ip)*sinnt + xxx(ip)*cosnt)
                tbz(ip)=tbz(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*sinnt + alphap*coeff*cosnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)*cos(ttt(ip)) - tbz(ip)*sin(ttt(ip))
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbx(ip)*sin(ttt(ip)) + tbz(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbz(ip)
          enddo
        endif

      endif

      if (mmlt2s .and. linmmlt2) then
c  Apply 2nd element magnetostatic multipoles for overlapping elements 

c       --- Get location of particle relative to moment data.
        call mltlocatxy(np,xp,yp,zp,nmmltsets,dzmmlt,nzmmlt,
     &                  nzl,cmmlt2zs,cmmlt2ze,cmmlt2ph,cmmlt2sf,cmmlt2sc,
     &                  cmmlt2im,cmmlt2ox,cmmlt2oy,
     &                  tbx,tby,tbz,imm,izm,wzm,xxx,yyy,zzz,ttt,rr,tt,
     &                  mltph,mltsf,
     &                  zlmin,zbeam,dzli,bends,linbend,cbendzs,cbendze,cbendrc)

c       --- accumulate the B field of the magnetostatic multipoles
        do in=1,nmsmult
          n = mmlt_n(in)
          v = mmlt_v(in)
          if (n == 0) then
c           --- Apply solenoidal field component and it's pseudomultipoles.
            f = 1./(2.*v+2.)
            do ip=1,np
              if (imm(ip) > 0) then
                coeffp = mltsf(ip)*
     &                    (msmmltp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                     msmmltp(izm(ip)+1,in,imm(ip))*wzm(ip))
                coeff = mltsf(ip)*
     &                    (msmmlt(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                     msmmlt(izm(ip)+1,in,imm(ip))*wzm(ip))
                rpow = dvnz(rr(ip))**(2*v)
                tbx(ip) = tbx(ip) - coeffp*rpow*f*xxx(ip)
                tby(ip) = tby(ip) - coeffp*rpow*f*yyy(ip)
                tbz(ip) = tbz(ip) + coeff*rpow
              endif
            enddo
          else
c           --- Apply the rest of the multipoles.
            f = 1. + 2.*v/n
            fz = 1./n
            do ip=1,np
              if (imm(ip) > 0) then
                alpha = tt(ip) - mltph(ip) -
     &                      msmmltph (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) -
     &                      msmmltph (izm(ip)+1,in,imm(ip))*      wzm(ip)
                alphap =   (msmmltphp(izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltphp(izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeff = mltsf(ip)*
     &                     (msmmlt   (izm(ip)  ,in,imm(ip))*(1. - wzm(ip)) +
     &                      msmmlt   (izm(ip)+1,in,imm(ip))*      wzm(ip))
                coeffp = mltsf(ip)*
     &                     (msmmltp  (izm(ip)  ,in,imm(ip))*(1. - wzm(ip))+
     &                      msmmltp  (izm(ip)+1,in,imm(ip))*      wzm(ip))
                cosnt = cos(n*alpha)
                sinnt = sin(n*alpha)
                rpow = dvnz(rr(ip))**(n-2+2*v)
                tbx(ip)=tbx(ip) + coeff*rpow*(f*xxx(ip)*sinnt - yyy(ip)*cosnt)
                tby(ip)=tby(ip) + coeff*rpow*(f*yyy(ip)*sinnt + xxx(ip)*cosnt)
                tbz(ip)=tbz(ip) + fz*dvnz(rr(ip))**(n+2*v)*
     &                            (coeffp*sinnt + alphap*coeff*cosnt)
              endif
            enddo
          endif
        enddo

c       --- do coordinate transform on fields back to warped coordinates
        if (bends .and. linbend) then
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)*cos(ttt(ip)) - tbz(ip)*sin(ttt(ip))
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbx(ip)*sin(ttt(ip)) + tbz(ip)*cos(ttt(ip))
          enddo
        else
          do ip=1,np
            bx(ip) = bx(ip) + tbx(ip)
            by(ip) = by(ip) + tby(ip)
            bz(ip) = bz(ip) + tbz(ip)
          enddo
        endif

      endif

c  ------------------------------------------------------------
c  End of electrostatic and magnetostatic multipole components
c  ------------------------------------------------------------

c  ------------------------------------------------------------
c  Magnetic fields from 3-D grid
c  ------------------------------------------------------------
      if (bgrds .and. linbgrd) then
        in = (1 + bgrdnx)*(1 + bgrdny)*(1 + bgrdnz)
        call setbgrdxy(in,bgrdns,bgrdbx,bgrdby,bgrdbz,np,xp,yp,zbeam,
     &               nbgrd,bgrdxs,bgrdys,cbgrdzs,bgrdid,bgrdsf,bgrdsc,
     &               bgrdsy,bgrdox,bgrdoy,bgrdph,bgrdsp,bgrdcp,  
     &               bgrddxi,bgrddyi,bgrddzi,
     &               bgrdnx,bgrdny,bgrdnz,bx,by,bz,
     &               cbgrdid,nzl,zlmin,zbeam,dzli)
      endif

      if (bgrd2s .and. linbgrd2) then
        in = (1 + bgrdnx)*(1 + bgrdny)*(1 + bgrdnz)
        call setbgrdxy(in,bgrdns,bgrdbx,bgrdby,bgrdbz,np,xp,yp,zbeam,
     &               nbgrd2,bgrd2xs,bgrd2ys,cbgrd2zs,bgrd2id,bgrd2sf,bgrd2sc,
     &               bgrd2sy,bgrd2ox,bgrd2oy,bgrd2ph,bgrd2sp,bgrd2cp,  
     &               bgrddxi,bgrddyi,bgrddzi,
     &               bgrdnx,bgrdny,bgrdnz,bx,by,bz,
     &               cbgrd2id,nzl,zlmin,zbeam,dzli)
      endif

c  ------------------------------------------------------------
c  Electrostatic potential from 3-D grid
c  ------------------------------------------------------------
      if (pgrds .and. linpgrd) then
        in = (1 + pgrdnx)*(1 + pgrdny)*(3 + pgrdnz)
        call setpgrdxy(in,pgrdns,pgrd,np,xp,yp,zbeam,
     &               npgrd,pgrdxs,pgrdys,cpgrdzs,pgrdid,pgrdsf,pgrdsc,
     &               pgrdox,pgrdoy,pgrdph,pgrdsp,pgrdcp, 
     &               pgrddxi,pgrddyi,pgrddzi,
     &               pgrdnx,pgrdny,pgrdnz,ex,ey,ez,
     &               cpgrdid,nzl,zlmin,zbeam,dzli)
      endif

c  Set E to zero for dead particles

      do ip=1,np
        if (uzp(ip) == 0) then
          ex(ip) = 0.
          ey(ip) = 0.
          ez(ip) = 0.
        endif
      enddo

c  Accumulate time for applying fields from the lattice
      latticetime = latticetime + (wtime() - timetemp)

      return
      end
c=============================================================================
      subroutine mltlocatxy(np,xp,yp,zp,nmltsets,dzmlt,nzmlt,
     &                    nzl,cmltzs,cmltze,cmltph,cmltsf,cmltsc,
     &                    cmltid,cmltox,cmltoy,
     &                    tx,ty,tz,imm,izm,wzm,xxx,yyy,zzz,ttt,rr,tt,
     &                    mltph,mltsf,
     &                    zlmin,zbeam,dzli,bends,linbend,
     &                    cbendzs,cbendze,cbendrc)
      integer     :: np,nmltsets,nzl
      real(kind=8):: xp(np),yp(np),zp(np)
      real(kind=8):: dzmlt(nmltsets)
      integer     :: nzmlt(nmltsets)
      real(kind=8):: cmltzs(0:nzl),cmltze(0:nzl),cmltph(0:nzl)
      real(kind=8):: cmltsf(0:nzl),cmltsc(0:nzl)
      integer     :: cmltid(0:nzl)
      real(kind=8):: cmltox(0:nzl),cmltoy(0:nzl)
      real(kind=8):: tx(np),ty(np),tz(np)
      integer     :: imm(np),izm(np)
      real(kind=8):: wzm(np),xxx(np),yyy(np),zzz(np)
      real(kind=8):: ttt(np),rr(np),tt(np),mltph(np),mltsf(np)
      logical     :: bends,linbend
      real(kind=8):: cbendzs(0:nzl),cbendze(0:nzl),cbendrc(0:nzl)
      real(kind=8):: zlmin,zbeam,dzli

c     --- Calculate grid location and polar coordinates of particles
c     --- for mlt elements.  Includes change of coordinates out of
c     --- warped coordinates (when a bend is between the
c     --- particle and the nearest mult).
c  Assumes that at most one bend will be between particle and the mult center.

      integer     :: j,ip
      real(kind=8):: cent,dzi,rrr,bzs,bze

c     --- All particles lie at the same z-cell
      j = 0

      if (.not. bends .or. .not. linbend) then
        do ip=1,np
          imm(ip) = cmltid(j)
          mltph(ip) = cmltph(j)
          mltsf(ip) = cmltsc(j) + cmltsf(j)
          dzi = 1./dzmlt(imm(ip))
c         --- zero temporary field arrays
          tx(ip) = 0.
          ty(ip) = 0.
          tz(ip) = 0.
c         --- find z location relative to multipole data
          izm(ip) = int((zbeam-cmltzs(j))*dzi)
          wzm(ip) =     (zbeam-cmltzs(j))*dzi - izm(ip)
          if (izm(ip) < 1) then
            izm(ip) = 1
            wzm(ip) = 0
            imm(ip) = 0
          elseif (izm(ip) > nzmlt(imm(ip))-1) then
            izm(ip) = nzmlt(imm(ip))-1
            wzm(ip) = 0
            imm(ip) = 0
          endif
c         --- Set temporary coordinates and add offset if particle within
c         --- an element.
          if (imm(ip) > 0) then
            xxx(ip) = xp(ip) - cmltox(j)
            yyy(ip) = yp(ip) - cmltoy(j)
            zzz(ip) = zbeam
            rr(ip) = sqrt(xxx(ip)**2 + yyy(ip)**2)
            tt(ip) = atan2(yyy(ip),dvnz(xxx(ip)))
          endif
        enddo
      else
        do ip=1,np
          cent = 0.5*(cmltzs(j) + cmltze(j))
c         --- set temporary x and z
          xxx(ip) = xp(ip)
          zzz(ip) = zbeam
          imm(ip) = cmltid(j)
          dzi = 1./dzmlt(imm(ip))
c         --- zero temporary field arrays
          tx(ip) = 0.
          ty(ip) = 0.
          tz(ip) = 0.
          ttt(ip) = 0.
          if (cmltzs(j) <= zbeam .and. zbeam <= cmltze(j)) then
            mltph(ip) = cmltph(j)
            mltsf(ip) = cmltsc(j) + cmltsf(j)
c           --- apply coordinate change left of mult, using index of (j)
            if (zbeam <= cent) then
c             --- find end of bend: if bend extends beyond center of mult, use
c             --- that as the end of the bend since moments are in coordinate
c             --- system at the center of the mult
              if (cbendze(j) < cent) then
                bze = cbendze(j)
              else
                bze = cent
              endif
c             --- if particle is not in bend and full bend is between particle
c             --- and mult
              rrr = xp(ip) + cbendrc(j)
              if (zbeam < cbendzs(j) .and. cbendzs(j) < cent) then
                ttt(ip) = (bze - cbendzs(j))/cbendrc(j)
                zzz(ip) = bze - rrr*sin(ttt(ip)) -
     &                    (cbendzs(j) - zbeam)*cos(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(j) -
     &                    (cbendzs(j) - zbeam)*sin(ttt(ip))
c             --- if particle is in bend
              elseif (zbeam < bze .and. cbendzs(j) < cent) then
                ttt(ip) = (bze - zbeam)/cbendrc(j)
                zzz(ip) = bze - rrr*sin(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(j)
              endif
c           --- apply coordinate change right of mult, using index of (j)
            elseif (zbeam >= cent) then
c             --- find start of bend: if bend extends beyond center of mult, use
c             --- that as the start of the bend since moments are in coordinate
c             --- system at the center of the mult
              if (cbendzs(j) < cent) then
                bzs = cent
              else
                bzs = cbendzs(j)
              endif
c             --- if particle is not in bend and full bend is between particle
c             --- and mult
              rrr = xp(ip) + cbendrc(j)
              if (zbeam > cbendze(j) .and. cbendze(j) > cent) then
                ttt(ip) =  - (cbendze(j) - bzs)/cbendrc(j)
                zzz(ip) = bzs - rrr*sin(ttt(ip)) +
     &                    (zbeam - cbendze(j))*cos(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(j) +
     &                    (zbeam - cbendze(j))*sin(ttt(ip))
c             --- if particle is in bend
              elseif (zbeam > bzs .and. cbendze(j) > cent) then
                ttt(ip) = - (zbeam - bzs)/cbendrc(j)
                zzz(ip) = bzs - rrr*sin(ttt(ip))
                xxx(ip) = rrr*cos(ttt(ip)) - cbendrc(j)
              endif
            endif
          endif
c         --- end of coordinate transformation

          izm(ip) = int((zzz(ip)-cmltzs(j))*dzi)
          wzm(ip) =     (zzz(ip)-cmltzs(j))*dzi - izm(ip)
          if (izm(ip) < 1) then
            izm(ip) = 1
            wzm(ip) = 0
            imm(ip) = 0
          elseif (izm(ip) > nzmlt(imm(ip))-1) then
            izm(ip) = nzmlt(imm(ip))-1
            wzm(ip) = 0
            imm(ip) = 0
          endif
c         --- Add offset
          if (imm(ip) > 0) then
            xxx(ip) = xxx(ip) - cmltox(j)
            yyy(ip) = yp(ip) -  cmltoy(j)
            rr(ip) = sqrt(xxx(ip)**2 + yyy(ip)**2)
            tt(ip) = atan2(yyy(ip),dvnz(xxx(ip)))
          endif
        enddo
      endif

      return
      end
c=============================================================================
      subroutine acclxy(np,xp,zp,uzp,gaminv,dtp,dtl,dtr,ez,m,q,zbeam,dt)
      use Lattice
      use LatticeInternal
      integer     :: np
      real(kind=8):: dtl,dtr,m,q,zbeam,dt
      real(kind=8):: xp(np), zp(np)
      real(kind=8):: uzp(np), gaminv(np), dtp(np)
      real(kind=8):: ez(np)

c Apply acceleration for finite length gaps.
c This is done seperately from the rest of the elements (which are done
c in extebxy) since the residence corrections depends on dtp which
c changes during the iteration dealing with changing Vz.  It would be much
c less efficient to put the entire extebxy within that iteration loop so
c the accl elements was seperated out.


      real(kind=8):: dti,oneodt,qoverm,moverq,gapez,z1,z2,vn,frac,cacclz,zl,zr
      real(kind=8):: dtl_a,dtr_a,dti_a,vz_a
      integer     :: j,ip

c  --- Note that the accl elements in the slice code are slightly different
c  --- than other hard edged elements since the residence correction is
c  --- based on fraction of time spent inside element rather than fraction of
c  --- distance.  So, the particles individual velocity and dt are used.

      if (.not. accls) return

      dti = 1./(dtr-dtl)
      oneodt = 1./dt
      qoverm = q/m

c     --- All particles lie at the same z-cell
      j = 0

      if (.not. lacclzl) then

c       --- For velocity correction, first calculate the velocity at
c       --- time level n.
        do ip=1,np
          dtl_a = dtl*oneodt*dtp(ip)
          dtr_a = dtr*oneodt*dtp(ip)
          dti_a = dti*dt/dtp(ip)
          vz_a = uzp(ip)*gaminv(ip)
          if (vz_a == 0.) vz_a = LARGEPOS
          gapez = cacclez(j) + cacclxw(j)*xp(ip)
          z1 = zbeam + dtl_a*vz_a - 0.5*gapez*qoverm*dtl_a**2
          z2 = zbeam + dtl_a*vz_a
          if (zbeam <= cacclzs(j)) then
            vn = vz_a
          elseif (z1 < cacclzs(j)) then
            vn = sqrt(vz_a**2 + 2.*gapez*qoverm*(zbeam - cacclzs(j)))
          elseif (zbeam <= cacclze(j)) then
            vn = vz_a - gapez*qoverm*dtl_a
          elseif (z2 < cacclze(j)) then
            vn = 0.5*(vz_a + 0.5*gapez*qoverm*dtp(ip) + sqrt((vz_a +
     &           0.5*gapez*qoverm*dtp(ip))**2 -
     &           4.*gapez*qoverm*(zbeam - cacclze(j))))
          else
            vn = vz_a
          endif
c         --- Calculate the fraction of time in the gap.  Cases inside
c         --- and outside are included implicitly in the max and min calls.
c         --- Note that the max's inside the sqrt are for idiot proofing.
          if (zbeam <= cacclzs(j)) then
            frac = max((dtr_a + (zbeam - cacclzs(j))/vn)*dti_a , 0.)
          elseif (zbeam <= (cacclzs(j) + cacclze(j))*0.5) then
            frac = min((dtr_a + 2.*(zbeam - cacclzs(j))/
     &             (sqrt(max(0.,vn**2 - 2.*gapez*qoverm*
     &             (zbeam - cacclzs(j)))) + vn))*dti_a, 1.)
          elseif (zbeam <= cacclze(j)) then
            frac = min((-dtl_a + 2.*(cacclze(j) - zbeam)/
     &             (sqrt(max(0.,vn**2 + 2.*gapez*qoverm*
     &             (cacclze(j) - zbeam))) + vn))*dti_a, 1.)
          else
            frac = max((-dtl_a - (zbeam - cacclze(j))/vn)*dti_a , 0.)
          endif
c         --- add acceleration field to Ez field
          ez(ip) = ez(ip) + gapez*frac
        enddo

      else

c       --- Zero-length gaps.

        moverq = m/q
        cacclz = 0.5*(cacclzs(j) + cacclze(j))
c       --- Add acceleration to velocity.
c       --- Calculate the change in velocity and convert that into
c       --- and Ez - the expression reduces to q*V/m for small V.
        do ip = 1, np
          dtl_a = dtl*oneodt*dtp(ip)
          dtr_a = dtr*oneodt*dtp(ip)
          dti_a = dti*dt/dtp(ip)
          vz_a = uzp(ip)*gaminv(ip)
          zl = zbeam + dtl_a*vz_a
          zr = zbeam + dtr_a*vz_a
          if (zl <= cacclz .and. cacclz < zr) then
            gapez = cacclez(j) + cacclxw(j)*xp(ip)
            ez(ip) = ez(ip) + (sqrt(vz_a**2 +
     &          2.*qoverm*gapez*(cacclze(j)-cacclzs(j))) - vz_a)*moverq*dti_a
          endif
        enddo

      endif

c  ------------------------------------------------------------
c  End of acceleration
c  ------------------------------------------------------------

      return
      end
c=============================================================================
      subroutine otherexy (np,xp,yp,dedr,dexdx,deydy,dbdr,ex,ey,ez,bx,by,bz)
      integer     :: np
      real(kind=8):: dedr,dexdx,deydy,dbdr
      real(kind=8):: xp(np), yp(np)
      real(kind=8):: ex(np), ey(np), ez(np)
      real(kind=8):: bx(np), by(np), bz(np)

      integer     :: ip

c  Set the electric fields from external sources,
c  inculding uniform focusing fields.

c  uniform focusing force

c     --- radial electric field
      if (dedr .ne. 0.) then
        do ip=1,np
          ex(ip) = ex(ip) + dedr*xp(ip)
          ey(ip) = ey(ip) + dedr*yp(ip)
        enddo
      endif
c     --- x- and y-electric fields
      if ((dexdx .ne. 0.) .or. (deydy .ne. 0.)) then
        do ip=1,np
          ex(ip) = ex(ip) + dexdx*xp(ip)
          ey(ip) = ey(ip) + deydy*yp(ip)
        enddo
      endif
c     --- azimuthal magnetic field
      if (dbdr .ne. 0.) then
        do ip=1,np
          bx(ip) = bx(ip) - dbdr*yp(ip)
          by(ip) = by(ip) + dbdr*xp(ip)
        enddo
      endif

      return
      end
c=============================================================================
      subroutine padvncxy(center)
      use InMesh3d
      use InGen
      use InGen3d
      use InGenxy
      use InPart
      use InGaps
      use LatticeInternal
      use Particles
      use Particlesxy
      use Fields3d
      use Efields3d
      use Picglb
      use Picglb3d
      use Beam_acc
      use Z_arrays
      use Damped_eom
      use Lattice
      character*(*):: center

c  Advances the particles position and velocity according to CENTER,
c  and also loads RHO at the new time level.


      integer     :: is,ip,ipmin,i
      real(kind=8):: uxpadv,uypadv,uzpadv
      real(kind=8):: ex(NPARPGRP), ey(NPARPGRP), ez(NPARPGRP), ezo(NPARPGRP)
      real(kind=8):: bx(NPARPGRP), by(NPARPGRP), bz(NPARPGRP)
      real(kind=8):: xpo(NPARPGRP), ypo(NPARPGRP), zpo(NPARPGRP)
      real(kind=8):: uxpo(NPARPGRP), uypo(NPARPGRP), uzpo(NPARPGRP)
      real(kind=8):: gaminvo(NPARPGRP)
      real(kind=8):: bendres, bendradi
      integer     :: iter_dt

c  Zero the bend radius and residence fraction
      bendradi = 0.
      bendres = 0.

c  Zero rho only if particles are advanced
c     if (center == "halfv") call zeroarry(rho,(nx+1)*(ny+1)*(nz+1))
      if (center == "halfv") rho = 0.
      
c  Zero curr if center is synchv or if using gaps
c     if (center == "synchv" .or. ifgap) call zeroarry(curr,nzzarr+1)
      if (center == "synchv" .or. ifgap) curr = 0.

c  Zero the moments if center is synchv
      if (ldiag .and. (center == "synchv" .or. center == "gen")) then
        call getzmmnt(1,xp,yp,zp,uxp,uyp,uzp,gaminv,sq(1),sm(1),sw(1),dt*0.5,1,
     &                nplive,uxpo,uypo,uzpo,1,ns)
      endif

c  Loop over species
      do is=1,ns

c  Loop over particle blocks; move each block separately
        do ipmin = ins(is), ins(is) + nps(is) - 1, NPARPGRP
          ip = min0 (NPARPGRP, ins(is)+nps(is)-ipmin)

c         --- Obtain the self-field from the electrostatic potential
          call sete3d(phi(0,0,-1),selfe,ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                zgridprv,xmmin,ymmin,zmmin,dx,dy,dz,nx,ny,nz,
     &                efetch,ex,ey,ez,l2symtry,l4symtry)
c         --- Scale the self E-field to get the lowest order relativistic
c         --- correction.
          if (relativity == 1)
     &      call sete3d_relativity(ip,ex,ey,vbeam)

c         --- Zero out B field arrays
          do i = 1,ip
             bx(i) = 0.
             by(i) = 0.
             bz(i) = 0.
          enddo
c         --- Compute lag average for experimental damping algorithm
          if (eomdamp .ne. 0.)
     &      call edamp (eomdamp,it,itdamp,center,ip,ex,ey,ez,
     &                  exold(ipmin),eyold(ipmin),ezold(ipmin),
     &                  exlag(ipmin),eylag(ipmin),ezlag(ipmin))

c         --- HALFV
          if (center == "halfv") then
c           --- Obtain bend radii and residence factors
            call getbendxy(zbeam,vbeam,bendres,bendradi,0.,dt)
c           --- Correct Ez_self for warped mesh effect
            call bendezxy (ip,xp(ipmin),zp(ipmin),ez,
     &                     bendres,bendradi,bends,bnezflag,linbend)
c           --- Add in Ez from axially-smoothed gaps 
            call gapfield (ip,zp(ipmin),ez,zbeam,zzmin,egap,dzz)
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex,ey,ez,bx,by,bz)
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),0.,dt*0.5,dodec,bz0,bx,by,bz,
     &                    ex,ey,ez,sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)

c           --- Save current position and velocity and ez
            if (lvzchang) then
              do i=1,ip
                xpo(i) = xp(ipmin+i-1)
                ypo(i) = yp(ipmin+i-1)
                zpo(i) = zp(ipmin+i-1)
                uxpo(i) = uxp(ipmin+i-1)
                uypo(i) = uyp(ipmin+i-1)
                uzpo(i) = uzp(ipmin+i-1)
                gaminvo(i) = gaminv(ipmin+i-1)
                ezo(i) = ez(i)
              enddo
            endif

c           --- Iterate over advance to calculate dt
c           --- (Only done when lvzchang is true.)
            do iter_dt=1,niter_dt
c             --- Get field from accelerating elements (depends on dtp)
              call acclxy (ip,xp(ipmin),zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                     dtp(ipmin),0.,dt*0.5,ez,sm(is),sq(is),zbeam,dt)
c             --- Correction to z on entry/exit to accelerator gap
              call zgapcorrxy(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                        dtp(ipmin),zbeam,0., dt*0.5, dt,
     &                        sm(1), sq(1), time)
c             --- Magnetic field increment to momenta
              call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      bx, by, bz, sq(is), sm(is), dtp(ipmin),0.5,ibpush)
c             --- Final half-electric field increment to momenta
              call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                      ex, ey, ez, sq(is), sm(is), dtp(ipmin), 0.5)
c             --- Advance relativistic Gamma factor
              call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                      gamadv,lrelativ)
c             --- Position advance
              call xpushxy (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                      uyp(ipmin),uzp(ipmin),gaminv(ipmin),dtp(ipmin))
c             --- Simplified translation of position for warped mesh effect
              call bendcorxy(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),dtp(ipmin),zpo,ds,bendres,bendradi,
     &                       lexbend)

c             --- Calculate dt for particles.
              call setdtp(ip,dtp(ipmin),xp(ipmin),zp(ipmin),zpo,uzp(ipmin),
     &                    bendres,bendradi)
c             --- Restore position and velocity if not finished with iteration
              if (lvzchang .and. iter_dt < niter_dt) then
                do i=1,ip
                  xp(ipmin+i-1) = xpo(i)
                  yp(ipmin+i-1) = ypo(i)
                  zp(ipmin+i-1) = zpo(i)
                  uxp(ipmin+i-1) = uxpo(i)
                  uyp(ipmin+i-1) = uypo(i)
                  uzp(ipmin+i-1) = uzpo(i)
                  gaminv(ipmin+i-1) = gaminvo(i)
                  ez(i) = ezo(i)
                enddo
              endif
            enddo

c           --- Exact translation of position and velocity for warped
c           --- mesh effect
            call exbendcorxy(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       zpo,ds,bendres,bendradi,lexbend)
c           --- Correct position advance for slanted dipole entry/exit
            call sledgcorxy(ip,xp(ipmin),yp(ipmin),zbeam+vbeam*dt,zbeam,
     &                      uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      sm(is),sq(is),dt,zbeam,vbeam,dipos) 
c         --- PUSHV
          elseif (center == "pushv") then
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex,ey,ez,bx,by,bz)
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),-SMALLPOS,SMALLPOS,dodec,bz0,bx,by,bz,
     &                    ex,ey,ez,sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx, by, bz, sq(is), sm(is), dtp(ipmin),0.5,ibpush)
c           --- Final half-electric field increment to momenta
            call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                    ex, ey, ez, sq(is), sm(is), dtp(ipmin), 1.0)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx, by, bz, sq(is), sm(is), dtp(ipmin),0.5,ibpush)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)

c         --- PUSHV
          elseif (center == "lpushv") then
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex,ey,ez,bx,by,bz)
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),-dt*0.5,0.,dodec,bz0,bx,by,bz,
     &                    ex,ey,ez,sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)
c           --- Final half-electric field increment to momenta
            call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                    ex, ey, ez, sq(is), sm(is), dtp(ipmin), 0.5)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx, by, bz, sq(is), sm(is), dtp(ipmin),0.5,ibpush)

c         --- PUSHV
          elseif (center == "rpushv") then
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex,ey,ez,bx,by,bz)
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),0.,dt*0.5,dodec,bz0,bx,by,bz,
     &                    ex,ey,ez,sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx, by, bz, sq(is), sm(is), dtp(ipmin),0.5,ibpush)
c           --- Final half-electric field increment to momenta
            call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                    ex, ey, ez, sq(is), sm(is), dtp(ipmin), 0.5)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)

c         --- PUSHX
          elseif (center == "pushx") then
c             --- Position advance
              call xpushxy (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                      uyp(ipmin),uzp(ipmin),gaminv(ipmin),dtp(ipmin))

c         --- SYNCHV or GEN
          elseif (center == "synchv" .or. center == "gen") then
c            --- Copy 'old' velocity into uxpo, uypo, and uzpo
             do i=1,ip
                uxpo(i) = uxp(ipmin+i-1)
                uypo(i) = uyp(ipmin+i-1)
                uzpo(i) = uzp(ipmin+i-1)
             enddo
c            --- Obtain bend radii and residence factors
             call getbendxy(zbeam,vbeam,bendres,bendradi,-dt,0.)
c            --- Correct Ez_self for warped mesh effect
             call bendezxy (ip,xp(ipmin),zp(ipmin),ez,
     &                      bendres,bendradi,bends,bnezflag,linbend)
c            --- Add in Ez from axially-smoothed gaps 
             call gapfield (ip,zp(ipmin),ez,zbeam,zzmin,egap,dzz)
c            --- Add in ears and uniform focusing E field pieces
             call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                      ex,ey,ez,bx,by,bz)
c            --- Set quad, dipole E and B; All: Bz
             call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),-dt*0.5,0.,dodec,
     &                     bz0,bx,by,bz,ex,ey,ez,sm(is),sq(is),
     &                     bendres,bendradi,lexbend,
     &                     gammabar,zbeam,vbeam,dt,time)
c            --- Get field from accelerating elements (depends on dtp)
             call acclxy (ip,xp(ipmin),zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    dtp(ipmin),-dt*0.5,0.,ez,sm(is),sq(is),zbeam,dt)
c            --- Half electric field increment to momenta
             call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                     ex, ey, ez, sq(is), sm(is), dtp(ipmin), 0.5)
c            --- Half magnetic field increment to momenta
             call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                     bx, by, bz, sq(is), sm(is), dtp(ipmin),0.5,ibpush)
             if (center == "gen") then
c              --- Reset uxp to uxpo, set uxpo to half step backward
c              --- for interpolation in moments calculation
               do i=1,ip
                 uxpadv = uxp(ipmin+i-1)
                 uxp(ipmin+i-1) = uxpo(i)
                 uxpo(i) = uxp(ipmin+i-1) - (uxpadv - uxp(ipmin+i-1))
                 uypadv = uyp(ipmin+i-1)
                 uyp(ipmin+i-1) = uypo(i)
                 uypo(i) = uyp(ipmin+i-1) - (uypadv - uyp(ipmin+i-1))
                 uzpadv = uzp(ipmin+i-1)
                 uzp(ipmin+i-1) = uzpo(i)
                 uzpo(i) = uzp(ipmin+i-1) - (uzpadv - uzp(ipmin+i-1))
               enddo
             endif
c            --- Calculate moments over particles, now that we're sync'd
             if (ldiag)
     &         call getzmmnt(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                       uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),sq(1),sm(is),sw(is),dt*0.5,2,
     &                       nplive,uxpo,uypo,uzpo,is,ns)
c            --- Calculate current
             call setcurrxy(curr(0),ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      sq(is),sw(is),zbeam,dzz,zzmin,dz,lthick)
          endif

          if (center == "halfv") then
c            --- Apply particle boundary conditions
             if (periinz)
     &          call periz(ip,zp(ipmin),zgrid,zmmax,zmmin)
             if (stickyz)
     &          call stckyz(ip,zp(ipmin),zmmax,zmmin,dz,uxp(ipmin),
     &                      uyp(ipmin),uzp(ipmin),zgrid)
             if (stickyxy)
     &          call stckxy3d(ip,xp(ipmin),xmmax,xmmin,dx,yp(ipmin),ymmax,
     &                        ymmin,dy,zp(ipmin),zmmin,dz,uxp(ipmin),
     &                        uyp(ipmin),uzp(ipmin),zgrid,zbeam,
     &                        nzzarr,prwallz,prwallxz,prwallyz,prwelips,
     &                        l2symtry,l4symtry,lostpars,zzmin,dzzi)

c            --- Collect charge density (xy) and current (1d)
             call setrhoxy(rho(0,0,0),ip,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
     &                     uzp(ipmin),gaminv(ipmin),sq(is),sw(is))

             if (ifgap) 
     &         call setcurrxy(curr(0),ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                        sq(is),sw(is),zbeam,dzz,zzmin,dz,lthick)
          endif
        enddo

      enddo
c     --- End loop over species

c     ---  Do final stuff for moments calculation and fix current
      if (center == "synchv" .or. center == "gen") then

        if (ldiag) then
          call getzmmnt(1,xp,yp,zp,uxp,uyp,uzp,gaminv,sq(1),sm(1),sw(1),
     &                  dt*0.5,3,nplive,uxpo,uypo,uzpo,1,ns)
          call getlabwn()
c         --- load rho diagnostic qtys
          call rhodia3d()
        endif

c       --- Fix the current
        call fixcurrxy(curr,nzzarr,periinz,lthick)

      endif

      if (center == "halfv") then
c       --- Make rho periodic if periodic in z
        call fixrhoxy(rho,nx,ny,nz,periinz,lthick)

c       --- Fix the current if it was calculated
        if (ifgap) call fixcurrxy(curr,nzzarr,periinz,lthick)

        if (lbeamcom) then
c         --- Set zbeam so that it follows the center of mass of the beam.
          call zbeamcom(zbeam)
        else
c         --- Advance beam frame location using the nominal beam frame velocity.
          zbeam = zgrid
        endif

c       --- zgridprv needs to be updated for the "synchv" step
c       --- Note that zgridprv is also set at the beginning of wxyexe
        zgridprv = zgrid

      endif

      return
      end
c=============================================================================
      subroutine fixrhoxy(rho,nx,ny,nz,periinz,lthick)
      integer     :: nx,ny,nz
      real(kind=8):: rho(0:nx,0:ny,0:nz)
      logical     :: periinz,lthick

c  Sums the first and last slices of rho for periodicity
c  and puts the result into both slices.
c  For thin slice model, average out rho along the axis.

      integer     :: ix,iy,iz
      real(kind=8):: nzi

      if (periinz) then

        do ix=0,nx
          do iy=0,ny
            rho(ix,iy,0)  = rho(ix,iy,0) + rho(ix,iy,nz)
            rho(ix,iy,nz) = rho(ix,iy,0)
          enddo
        enddo

      endif

#ifdef PARALLEL
c     --- Sum up rho over the processors.
      if (.not. lthick) then
        call parallelsum(rho(0,0,0),(1+nx)*(1+ny))
      endif
#endif

      if (.not. lthick) then
c       --- In the setrhoxy routine, all of the charge was put at iz=0.

        if (nz > 1) then
          nzi = 1./real(nz)
          do iy=0,ny
            do ix=0,nx
              rho(ix,iy,0) = rho(ix,iy,0)*nzi
            enddo
          enddo
        endif

c       --- Distribute rho to all z planes.
        do iz=1,nz
          do iy=0,ny
            do ix=0,nx
              rho(ix,iy,iz) = rho(ix,iy,0)
            enddo
          enddo
        enddo

      endif

      return
      end
c=============================================================================
      subroutine fixcurrxy(curr,nzzarr,periinz,lthick)
      integer     :: nzzarr
      real(kind=8):: curr(0:nzzarr)
      logical     :: periinz,lthick

c  Sums the first and last slices of curr for periodicity
c  and puts the result into both slices.
c  For thin slice model, distribute the current along the axis.

      integer     :: iz

      if (periinz) then

        curr(0)  = curr(0) + curr(nzzarr)
        curr(nzzarr) = curr(0)

      endif

      if (.not. lthick .and. nzzarr > 1) then
c       --- In the setcurrxy routine, all of the current was put at iz=0.

c       --- Distribute the current.
        do iz=1,nzzarr
          curr(iz) = curr(0)
        enddo

      endif

      return
      end
c=============================================================================
      subroutine epushxy (np, uxp, uyp, uzp, ex, ey, ez, q, m, dtp, fdt)
      integer     :: np
      real(kind=8):: uxp(np), uyp(np), uzp(np)
      real(kind=8):: ex(np), ey(np), ez(np)
      real(kind=8):: q, m, dtp(np), fdt

c  Push the particle velocity with E field

      integer     :: ip
      real(kind=8):: const
      const = q*fdt/m

      do ip=1,np
        uxp(ip) = uxp(ip) + ex(ip)*const*dtp(ip)
        uyp(ip) = uyp(ip) + ey(ip)*const*dtp(ip)
        uzp(ip) = uzp(ip) + ez(ip)*const*dtp(ip)
      enddo

      return
      end
c============================================================================= 
      subroutine bpushxy (np, uxp, uyp, uzp, gaminv, bx, by, bz, q, m, dtp,
     &                    fdt, ibpush)
      integer     :: np,ibpush
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np)
      real(kind=8):: bx(np), by(np), bz(np)
      real(kind=8):: q, m, dtp(np), fdt

c  Push the particle velocity with B field

      integer     :: ip
      real(kind=8):: tx,ty,tz,tsqi,sx,sy,sz,uxppr,uyppr,uzppr
      real(kind=8):: const
      const = q*fdt*0.5/m

      if (ibpush == 1) then
c        --- fast b-field rotation algorithm
         do ip=1,np
            tx = gaminv(ip)*bx(ip)*const*dtp(ip)
            ty = gaminv(ip)*by(ip)*const*dtp(ip)
            tz = gaminv(ip)*bz(ip)*const*dtp(ip)
            tsqi = 2./(1. + tx**2 + ty**2 + tz**2)
            sx = tx*tsqi
            sy = ty*tsqi
            sz = tz*tsqi
            uxppr = uxp(ip) + uyp(ip)*tz - uzp(ip)*ty
            uyppr = uyp(ip) + uzp(ip)*tx - uxp(ip)*tz
            uzppr = uzp(ip) + uxp(ip)*ty - uyp(ip)*tx
            uxp(ip) = uxp(ip) + uyppr*sz - uzppr*sy
            uyp(ip) = uyp(ip) + uzppr*sx - uxppr*sz
            uzp(ip) = uzp(ip) + uxppr*sy - uyppr*sx
         enddo
      elseif (ibpush == 2) then
c        --- tan(alpha) / alpha algorithm
         do ip=1,np
            tx = tan (gaminv(ip)*bx(ip)*const*dtp(ip))
            ty = tan (gaminv(ip)*by(ip)*const*dtp(ip))
            tz = tan (gaminv(ip)*bz(ip)*const*dtp(ip))
            tsqi = 2./(1. + tx**2 + ty**2 + tz**2)
            sx = tx*tsqi
            sy = ty*tsqi
            sz = tz*tsqi
            uxppr = uxp(ip) + uyp(ip)*tz - uzp(ip)*ty
            uyppr = uyp(ip) + uzp(ip)*tx - uxp(ip)*tz
            uzppr = uzp(ip) + uxp(ip)*ty - uyp(ip)*tx
            uxp(ip) = uxp(ip) + uyppr*sz - uzppr*sy
            uyp(ip) = uyp(ip) + uzppr*sx - uxppr*sz
            uzp(ip) = uzp(ip) + uxppr*sy - uyppr*sx
         enddo
      endif

      return
      end
c=========================================================================== 
      subroutine xpushxy (np, xp, yp, zp, uxp, uyp, uzp, gaminv, dtp)
      integer     :: np
      real(kind=8):: dtp(np)
      real(kind=8):: xp(np), yp(np), zp(np), uxp(np), uyp(np), uzp(np), gaminv(np)

c  Advance particle positions

      integer     :: ip

      do ip=1,np
        xp(ip) = xp(ip) + uxp(ip)*dtp(ip)*gaminv(ip)
        yp(ip) = yp(ip) + uyp(ip)*dtp(ip)*gaminv(ip)
        zp(ip) = zp(ip) + uzp(ip)*dtp(ip)*gaminv(ip)
      enddo

      return
      end      
c============================================================================ 
      subroutine setdtp(np,dtp,xp,zp,zpo,uzp,bendres,bendradi)
      use InGenxy
      integer     :: np
      real(kind=8):: dtp(np),xp(np),zp(np),zpo(np),uzp(np)
      real(kind=8):: bendres,bendradi

c Set dt for each of the particles.  With changing vz, dt needs to be
c recalculated for each particle each timestep.


      integer     :: ip
      real(kind=8):: dtheta0,dtheta,rr

      if (lvzchang) then

c       --- This algorithm scales dt by the ratio of the desired ds and
c       --- distance actually travelled.  Note that in a bend, the distance
c       --- traveled is calculated as an angle around the bend.

        if (bendres > 0. .and. lexbend) then

          dtheta0 = ds/abs(bendradi)

          do ip=1,np
            if (uzp(ip) .ne. 0.) then
              rr = sqrt((bendradi+xp(ip))**2 + (zp(ip) - zpo(ip))**2)
              dtheta = asin((zp(ip) - zpo(ip))/rr)
              dtp(ip) = dtp(ip)*dtheta0/dtheta
            endif
          enddo

        else

c         --- with no bends...
          do ip=1,np
            if (uzp(ip) .ne. 0.) then
              dtp(ip) = dtp(ip)*ds/(zp(ip) - zpo(ip))
            endif
          enddo

        endif

      endif

      return
      end      
c============================================================================ 
      subroutine nextbend(zbeam,zz)
      use InGen
      use Lattice
      use LatticeInternal
      real(kind=8):: zbeam,zz

c Get distance to get to start or end of bend, which ever is closer.


      integer     :: j

      if (bends) then

c       --- All particles lie at the same z-cell
        j = 0

        if (zbeam < cbendzs(j)) then
c         --- distance to start of bend
          zz = cbendzs(j) - zbeam
        elseif (zbeam < cbendze(j)) then
c         --- distance to end of bend
          zz = cbendze(j) - zbeam
        else
c         --- slice is completely past bend, return a big number
          zz = LARGEPOS
        endif

      else

c       --- If no bends, return a big number.
        zz = LARGEPOS

      endif

      return
      end
c============================================================================ 
      subroutine getbendxy(zbeam,vbeam,bendres,bendradi,dtl,dtr)
      use InGen
      use Lattice
      use LatticeInternal
      real(kind=8):: dtl,dtr,zbeam,vbeam
      real(kind=8):: bendres, bendradi


c  Gets residence factor and radius for bends.
c  For periodic runs, assumes mesh period length = lattice period length
c  Note that for the slice code, since the steps are arranged to always
c  land at the edge of a bend, the residence fraction will always be
c  either zero or one.

      integer     :: j
      real(kind=8):: vz,zl,fl,zr,fr,dti,vzi,br

      if (.not. bends .or. .not. linbend) return

      dti = 1./(dtr-dtl)

c  Extract the local bend properties off the general lattice

      vz = vbeam
      vzi = 1./dvnz(vz)
c     --- All particles lie at the same z-cell
      j = 0
c     --- "left" end of velocity advance step
      zl = zbeam + vz * dtl
      fl = 0.
      if (zl >= cbendzs(j) .and. zl < cbendze(j)) fl = 1.
c     --- "right" end of velocity advance step
      zr = zbeam + vz * dtr
      fr = 0.
      if (zr >= cbendzs(j) .and. zr < cbendze(j)) fr = 1.
c     --- residence fraction; set to all-in or all-out, possibly override
      if (fr == fl) then
c        --- all in, or all out
         br = fl
      elseif (fr > fl) then
c        --- entered bend this v-advance step
         br = (zr - cbendzs(j)) * vzi*dti
      else
c        --- exited dipole this v-advance step
         br = (cbendze(j) - zl) * vzi*dti
      endif
c     --- residence fraction
      if (br < 0.5) then
        bendres = 0.
      else
        bendres = 1.
      endif
c     --- bend radius 
      bendradi = cbendrc(j)

      return
      end
c=========================================================================== 
      subroutine bendezxy (np, xp, zp, ez, bendres, bendradi, bends, bnezflag,
     &                     linbend)
      integer     :: np
      real(kind=8):: xp(np), zp(np), ez(np), bendres, bendradi
      logical     :: bends, bnezflag, linbend

c  Corrects axial electric field at particle position for warped geometry
c  via multiplying by r_star/r = 1 - x/r, in a residence-corrected way;
c  at smaller radii, zones are closer together, so field is larger.

      integer     :: ip

      if (.not. (bends.and.bnezflag) .or. .not. linbend) return

      do ip=1,np
         ez(ip) = ez(ip) * ( 1. 
     &    - bendres * xp(ip) / (bendradi + xp(ip)) )
      enddo

      return
      end      
c=========================================================================== 
      subroutine sledgcorxy(np,xp,yp,zp,zpo,uxp,uyp,uzp,gaminv,m,q,dt,zbeam,
     &                      vbeam,dipos) 
      use LatticeInternal
      integer     :: np
      real(kind=8):: xp(np), yp(np), zp, zpo 
      real(kind=8):: uxp(np), uyp(np), uzp(np), gaminv(np) 
      real(kind=8):: q, m, dt, zbeam, vbeam
      logical     :: dipos 


c Applies position and velocity jump corrections on entry/exit of dipoles 
c with slanted faces.  (Must be called after all z advancement is done this 
c timestep, so that entry/exit of the bend is detected exactly once.)    
c Accounts for the fact that the jumped velocities act on the particles for 
c less than the entire andvance step on entry/exit of the dipole.  

      integer     :: ip,j
      real(kind=8):: vz,vzi,qoverm

      qoverm = q/m

      if (.not.dipos) return

      vz = vbeam
      vzi = 1./dvnz(vz)
c     --- All particles lie at the same z-cell
      j = 0
      if ((zpo < cdipozs(j)) .and. (zp >= cdipozs(j))) then 
c       --- entered dipole
        do ip = 1,np
          uxp(ip) = uxp(ip) + xp(ip)*cdipoby(j)*cdipota(j)*qoverm
          uyp(ip) = uyp(ip) - yp(ip)*cdipoby(j)*cdipota(j)*qoverm
          xp(ip)  = xp(ip) + (zp-cdipozs(j))*vzi*
     &                       xp(ip)*cdipoby(j)*cdipota(j)*qoverm
          yp(ip)  = yp(ip) - (zp-cdipozs(j))*vzi*
     &                       yp(ip)*cdipoby(j)*cdipota(j)*qoverm
        enddo
      elseif ((zpo < cdipoze(j)).and.(zp >= cdipoze(j))) then  
c       --- exited dipole
        do ip = 1,np
          uxp(ip) = uxp(ip) + xp(ip)*cdipoby(j)*cdipotb(j)*qoverm
          uyp(ip) = uyp(ip) - yp(ip)*cdipoby(j)*cdipotb(j)*qoverm
          xp(ip)  = xp(ip) + (zp-cdipoze(j))*vzi*
     &                        xp(ip)*cdipoby(j)*cdipotb(j)*qoverm
          yp(ip)  = yp(ip) - (zp-cdipoze(j))*vzi*
     &                       yp(ip)*cdipoby(j)*cdipotb(j)*qoverm
        enddo
      endif 

      return
      end
c=========================================================================== 
      subroutine exbendcorxy(np,xp,zp,uxp,uzp,zpo,ds,bendres,bendradi,lexbend)
      integer     :: np
      real(kind=8):: xp(np), zp(np), zpo(np), uxp(np), uzp(np)
      real(kind=8):: ds,bendres,bendradi
      logical     :: lexbend

c Exact bend translation:
c Translate the particle position and velocity to the new rotated
c frame (when in a bend).
c The expression for translating xp is written in a way that is correct
c for either sign of bendradi.

      integer     :: ip
      real(kind=8):: rr
      real(kind=8):: dtheta0,ty,sy,uxppr,uzppr

      if (bendres > 0. .and. lexbend) then

c       --- constants for all the particles
        dtheta0 = ds/bendradi
        ty = -tan(dtheta0/2.)
        sy = 2.*ty/(1. + ty*ty)

        do ip=1,np

c         --- translate position
          rr = sqrt((bendradi+xp(ip))**2 + (zp(ip) - zpo(ip))**2)
          xp(ip) = (rr - abs(bendradi))*sign(1.,bendradi)
          zp(ip) = zpo(ip) + ds

c         --- translate velocity
          uxppr = uxp(ip) - uzp(ip)*ty
          uzppr = uzp(ip) + uxp(ip)*ty
          uxp(ip) = uxp(ip) - uzppr*sy
          uzp(ip) = uzp(ip) + uxppr*sy

        enddo

      endif

      return
      end
c=========================================================================== 
      subroutine bendcorxy(np,xp,zp,uxp,uzp,gaminv,dtp,zpo,ds,
     &                     bendres,bendradi,lexbend)
      integer     :: np
      real(kind=8):: xp(np), zp(np), zpo(np), uxp(np), uzp(np), gaminv(np), dtp(np)
      real(kind=8):: ds,bendres,bendradi
      logical     :: lexbend

c Simple translation, same model as in 3-D code:
c Translate the particle position and velocity to the new rotated
c frame (when in a bend).

      integer     :: ip
      real(kind=8):: xprv,xc

      if (bendres > 0. .and. .not. lexbend) then

        do ip=1,np
          xprv = xp(ip) - dtp(ip)*uxp(ip)*gaminv(ip)
          xc = 0.5*(xp(ip) + xprv)
          zp(ip) = zp(ip) + dtp(ip)*uzp(ip)*gaminv(ip)*bendres*
     &                      (bendradi/(bendradi+xc) - 1.)
        enddo

      endif

      return
      end      
c=========================================================================== 
      subroutine zgapcorrxy(np,zp,xp,uzp,gaminv,dtp,zbeam,dtl,dtr,dt,m,q,time)
      use Lattice
      use LatticeInternal
      integer     :: np
      real(kind=8):: zp(np), xp(np), uzp(np), gaminv(np), dtp(np)
      real(kind=8):: zbeam,dtl,dtr,dt,m,q,time

c Add residence correction to particle position when acceleration is done


      integer     :: ip,j
      real(kind=8):: gapez,zp2,delta,qoverm,vzi
      real(kind=8):: vz_a,dtl_a,dtr_a
      real(kind=8):: cacclz,zr,vz,deltav

      if (.not. accls) return

      qoverm = q/m

c     --- All particles lie at the same z-cell
      j = 0

c     --- Finite length gap
      if (.not. lacclzl) then

        do ip=1,np
c         --- Scale time step size by particles time step.
          dtl_a = dtl/dt*dtp(ip)
          dtr_a = dtr/dt*dtp(ip)
c         --- Particle axial velocity and inverse
          vz_a = uzp(ip)*gaminv(ip)
          vzi = 1./dvnz(vz_a)
c         --- Fetch accelerating voltage
          gapez = cacclez(j) + cacclxw(j)*xp(ip)
c         --- Approximate position of particle on next step
c         --- with simple residence correction for gaps.
c         --- Calculate time of entrance or exit, delta.
c         --- Space charge forces are ignored.
          zp2 = zbeam + vz_a*dtp(ip)
          delta = 0
          if (zbeam <= cacclzs(j) .and. cacclzs(j) <= zp2) then
            delta = (cacclzs(j) - zbeam)*vzi
          elseif (zbeam <= cacclze(j) .and.
     &         cacclze(j) <= zp2 + gapez*qoverm*dtp(ip)*(dtr_a-dtl_a)) then
            delta = 2.*(cacclze(j) - zbeam)/
     &              (sqrt((vz_a - gapez*qoverm*dtl_a)**2 +
     &              2.*(cacclze(j) - zbeam)*gapez*qoverm) +
     &              vz_a - gapez*qoverm*dtl_a)
            zp2 = zp2 + gapez*qoverm*dtp(ip)*(dtr_a-dtl_a)
          endif
c         --- calculate correction on particle position
          if (delta > 0.5*dtp(ip)) delta = dtp(ip) - delta
          if (zbeam <= cacclzs(j) .and. cacclzs(j) <= zp2) then
            zp(ip) = zp(ip) + 0.5*gapez*qoverm*delta**2
          elseif (zbeam <= cacclze(j) .and. cacclze(j) <= zp2) then
            zp(ip) = zp(ip) - 0.5*gapez*qoverm*delta**2
          endif
        enddo

      else

c       --- Zero length gap

        cacclz = 0.5*(cacclzs(j) + cacclze(j))
        do ip = 1, np
          vz = uzp(ip)*gaminv(ip)
c         --- "right" end of velocity advance step
          zr = zbeam + vz*dtr
c         --- calculate correction on position
          if (zbeam < cacclz .and. cacclz < zr) then
            vzi = 1./dvnz(vz)
            gapez = cacclez(j) + cacclxw(j)*xp(ip)
            deltav = sqrt(vz**2 + 2.*qoverm*gapez*(cacclze(j)-cacclzs(j))) - vz
            zp(ip) = zp(ip) - (cacclz - zbeam)*vzi*deltav
          endif
          if (zr < cacclz .and. cacclz < zbeam + vz*dtp(ip)) then
            vzi = 1./dvnz(vz)
            gapez = cacclez(j) + cacclxw(j)*xp(ip)
            deltav = sqrt(vz**2 + 2.*qoverm*gapez*(cacclze(j)-cacclzs(j))) - vz
            zp(ip) = zp(ip) + (dtp(ip) - (cacclz - zbeam)*vzi)*deltav
          endif
        enddo

      endif

      return
      end
c=========================================================================== 
      subroutine setbgrdxy(nxyz,bgrdns,bgrdbx,bgrdby,bgrdbz,np,xp,yp,zp,
     &                     nbgrd,bgrdxs,bgrdys,cbgrdzs,bgrdid,bgrdsf,bgrdsc,
     &                     bgrdsy,bgrdox,bgrdoy,bgrdph,bgrdsp,bgrdcp,  
     &                     bgrddxi,bgrddyi,bgrddzi,bgrdnx,bgrdny,bgrdnz,
     &                     bx,by,bz,id,nzl,zlmin,zbeam,dzli)
      integer     :: nxyz,bgrdns,np,nbgrd,bgrdnx,bgrdny,bgrdnz,nzl
      real(kind=8):: bgrdbx(0:nxyz-1,bgrdns)
      real(kind=8):: bgrdby(0:nxyz-1,bgrdns)
      real(kind=8):: bgrdbz(0:nxyz-1,bgrdns)
      real(kind=8):: xp(np),yp(np),zp
      real(kind=8):: bx(np),by(np),bz(np)
      real(kind=8):: bgrdxs(0:nbgrd),bgrdys(0:nbgrd),cbgrdzs(0:nzl)
      real(kind=8):: bgrdsf(0:nbgrd),bgrdsc(0:nbgrd)
      real(kind=8):: bgrdox(0:nbgrd),bgrdoy(0:nbgrd)
      real(kind=8):: bgrdph(0:nbgrd),bgrdsp(0:nbgrd),bgrdcp(0:nbgrd) 
      integer     :: bgrdid(0:nbgrd),bgrdsy(0:nbgrd)
      real(kind=8):: bgrddxi(bgrdns),bgrddyi(bgrdns),bgrddzi(bgrdns)
      integer     :: id(0:nzl)
      real(kind=8):: zlmin,zbeam,dzli

c  Sets magnetic field for a group of up to NPARPGRP particles from data sets
c  containing Bx, By, and Bz on 3-D grids.
c  Calling arguments:
c     nxyz           Product (1+nx)*(1+ny)*(1+nz)
c     bgrdns         number of data sets
c     bgrdbx         |arrays containing B field. In calling routine, they are
c     bgrdby         |dimensioned (0:nx,0:ny,0:nz,bgrdns), bgrdns different
c     bgrdbz         |sets of fields
c                    |Here, the first three dimensions are reduced to one.
c     np             number of particles
c     xp,yp,zp       coordinates of particles
c     nbgrd          number of sets of B field arrays
c     bgrdxs,bgrdys  mins of meshes of each set of B field arrays
c     cbgrdzs        z mins of data sets stored in the style of internal lattice
c                    arrays
c     bgrdsf         scale factor of each set of B field array
c     bgrdsc         scale factor of each set of B field array
c     bgrdsy         level of symmetry of B field data
c     bgrdox,oy      transverse offset of element
c     bgrdph         rotation angle of element 
c     bgrdsp,cp      sine and cosine of bgrdph 
c     bgrdid         index of array for set of B field arrays
c     bgrddxi,y,z    one over the grid cell sizes of each set of B field arrays
c     bgrdnx,y,z     Dimensions of arrays in calling routine 
c     bx,by,bz       returned B fields for the particles
c     id             index to data sets stored in the style of internal lattice
c                    arrays
c     nzl            size of internal lattice arrays
c     zlmin          minimum location of internal lattice arrays relative to
c                    beam frame
c     zbeam          location of beam frame
c     dzli           one over internal lattice arrays grid cell size

c  Algorithm notes: bgrdbx, bgrdby, and bgrdbz arrays are dimensioned
c  (0:bgrdnx,0:bgrdny,0:bgrdnz,bgrdns)
c  outside, but are made two dimensional in this routine
c  so cell index into 1D array for vectorized deposition is:
c     i + j*(bgrdnx+1) + k*(bgrdnx+1)*(bgrdny+1)
c  The field is:
c     Bx = u0 * v0 * w0 * bgrdbx(i  ,j  ,k  ,bgrdid)
c        + u1 * v0 * w0 * bgrdbx(i+1,j  ,k  ,bgrdid)
c        + u0 * v1 * w0 * bgrdbx(i  ,j+1,k  ,bgrdid)
c        + ...

      integer     :: nnx,nnxy,ip,i,j,k,ind0,iz,ii,ib
      real(kind=8):: u0,u1,v0,v1,w0,w1,bfac,xsign,ysign,zsign,txp,typ,temp
      real(kind=8):: tbx,tby
      integer     :: noff(8)
      save noff

      nnx  = bgrdnx + 1
      nnxy = (bgrdnx + 1)*(bgrdny + 1)

c  Evaluation of offsets occurs on the first call only provided
c  local arrays are SAVE'd.

      if (noff(8) .ne. 1 + nnx + nnxy) then
         noff(1) = 0
         noff(2) = 1
         noff(3) =     nnx
         noff(4) = 1 + nnx
         noff(5) =           nnxy
         noff(6) = 1       + nnxy
         noff(7) =     nnx + nnxy
         noff(8) = 1 + nnx + nnxy
      endif

c  Evaluation of B, vectorized over particles

c     --- All particles lie in the same z-cell
      iz = 0
      ii = id(iz)
      ib = bgrdid(ii)

c     --- Particle z-location in field grid and z linear weight factor 
      k  = (zp - cbgrdzs(iz))*bgrddzi(ib)
      w1 = (zp - cbgrdzs(iz))*bgrddzi(ib) - k
      w0 = 1. - w1

c     --- Exit field accumulation if particles are outside of axial grid
      if (k < 0 .or. k >= bgrdnz) return 

c     --- Field scale factors 
      bfac = bgrdsc(ii) + bgrdsf(ii)

      do ip = 1, np

c        --- find transverse particle coordinate in frame of gridded field 
c            --- transverse offsets
         txp = xp(ip) - bgrdox(ii)
         typ = yp(ip) - bgrdoy(ii)
c            --- transverse rotation to take into account an active rotation 
c                of the field element.  Particles are rotated in that oposite 
c                sense of the element. Later the field components accumulated 
c                must be rotated back. 
         if ( bgrdph(ii) .ne. 0. ) then
           temp = txp  
           txp =  temp*bgrdcp(ii) + typ*bgrdsp(ii) 
           typ = -temp*bgrdsp(ii) + typ*bgrdcp(ii) 
         endif 

c        --- Shift coordinates to measure from the edge of the field grid
         txp = txp - bgrdxs(ii)
         typ = typ - bgrdys(ii)

c        --- Set default sign of B field
         xsign = 1.
         ysign = 1.
         zsign = 1.

c        --- If B is quadrupolar symmetric, make transformations.
c        --- When the particle is in one of the even quadrants (either
c        --- x<0 or y<0 but not both), the transformation is done by
c        --- swapping x and y, and by swapping Bx and By (done at the end
c        --- of the loop). The sign of Bz is also changed. In the third
c        --- quadrant, the signs of both Bx and By are changed.
c                Quadrupole symmetries on field grid
c          
c            Quadrant      B_x            B_y             B_z 
c            --------------------------------------------------------
c               I          B_x( x, y,z)   B_y( x, y,z)    B_z( x, y,z) 
c               II         B_y( y,-x,z)  -B_x( y,-x,z)   -B_z( y,-x,z)  
c               III       -B_x(-x,-y,z)  -B_y(-x,-y,z)    B_z(-x,-y,z) 
c               IV        -B_y(-y, x,z)   B_x(-y, x,z)   -B_z(-y, x,z) 
c
         if (bgrdsy(ii) == 2) then
c          --- Get quadrant that the particle is in.
           if (txp < 0.) then
             xsign = -1.
             txp = -txp
           endif
           if (typ < 0.) then
             ysign = -1.
             typ = -typ
           endif
c          --- If in even quadrant...
           if (xsign*ysign < 0.) then
c            --- Switch sign of Bz.
             zsign = -1.
c            --- Swap x and y
             temp = txp
             txp = typ
             typ = temp
           endif
         endif

c        --- find location of particle in B field grid
         i = txp*bgrddxi(ib)
         j = typ*bgrddyi(ib)
         ind0 = i + j*nnx + k*nnxy

c        --- Calculate linear weights
         u1 = txp*bgrddxi(ib) - i
         v1 = typ*bgrddyi(ib) - j
         u0 = 1. - u1
         v0 = 1. - v1

c        --- Only calculate for particles inside the B field grid
         if (i >= 0 .and. i < bgrdnx .and.
     &       j >= 0 .and. j < bgrdny) then

           tbx = xsign*bfac*(u0*v0*w0*bgrdbx(noff(1)+ind0,ib) +
     &                       u1*v0*w0*bgrdbx(noff(2)+ind0,ib) +
     &                       u0*v1*w0*bgrdbx(noff(3)+ind0,ib) +
     &                       u1*v1*w0*bgrdbx(noff(4)+ind0,ib) +
     &                       u0*v0*w1*bgrdbx(noff(5)+ind0,ib) +
     &                       u1*v0*w1*bgrdbx(noff(6)+ind0,ib) +
     &                       u0*v1*w1*bgrdbx(noff(7)+ind0,ib) +
     &                       u1*v1*w1*bgrdbx(noff(8)+ind0,ib))

           tby = ysign*bfac*(u0*v0*w0*bgrdby(noff(1)+ind0,ib) +
     &                       u1*v0*w0*bgrdby(noff(2)+ind0,ib) +
     &                       u0*v1*w0*bgrdby(noff(3)+ind0,ib) +
     &                       u1*v1*w0*bgrdby(noff(4)+ind0,ib) +
     &                       u0*v0*w1*bgrdby(noff(5)+ind0,ib) +
     &                       u1*v0*w1*bgrdby(noff(6)+ind0,ib) +
     &                       u0*v1*w1*bgrdby(noff(7)+ind0,ib) +
     &                       u1*v1*w1*bgrdby(noff(8)+ind0,ib))

           if (bgrdsy(ii) == 2 .and. xsign*ysign < 0.) then
             temp = tbx 
             tbx = tby
             tby = temp
           endif 

c          --- rotate transverse field componets for correct lab frame 
c              orientation if the element is rotated.  
           if ( bgrdph(ii) .ne. 0. ) then
             temp = tbx 
             tbx = temp*bgrdcp(ii) - tby*bgrdsp(ii) 
             tby = temp*bgrdsp(ii) + tby*bgrdcp(ii)         
           endif 

c          --- tranverse field 
           bx(ip) = bx(ip) + tbx
           by(ip) = by(ip) + tby

           bz(ip) = bz(ip) + zsign*bfac*(u0*v0*w0*bgrdbz(noff(1)+ind0,ib) +
     &                                   u1*v0*w0*bgrdbz(noff(2)+ind0,ib) +
     &                                   u0*v1*w0*bgrdbz(noff(3)+ind0,ib) +
     &                                   u1*v1*w0*bgrdbz(noff(4)+ind0,ib) +
     &                                   u0*v0*w1*bgrdbz(noff(5)+ind0,ib) +
     &                                   u1*v0*w1*bgrdbz(noff(6)+ind0,ib) +
     &                                   u0*v1*w1*bgrdbz(noff(7)+ind0,ib) +
     &                                   u1*v1*w1*bgrdbz(noff(8)+ind0,ib))

         endif

      enddo

      return
      end
c=========================================================================== 
      subroutine setpgrdxy(nxyz,pgrdns,pgrdp,np,xp,yp,zp,
     &                     npgrd,pgrdxs,pgrdys,cpgrdzs,pgrdid,pgrdsf,pgrdsc,
     &                     pgrdox,pgrdoy,pgrdph,pgrdsp,pgrdcp, 
     &                     pgrddxi,pgrddyi,pgrddzi,pgrdnx,pgrdny,pgrdnz,ex,ey,ez,
     &                     id,nzl,zlmin,zbeam,dzli)
      integer     :: nxyz,pgrdns,np,npgrd,pgrdnx,pgrdny,pgrdnz,nzl
      real(kind=8):: pgrdp(0:nxyz-1,pgrdns)
      real(kind=8):: xp(np),yp(np),zp
      real(kind=8):: ex(np),ey(np),ez(np)
      real(kind=8):: pgrdxs(0:npgrd),pgrdys(0:npgrd),cpgrdzs(0:nzl)
      real(kind=8):: pgrdsf(0:npgrd),pgrdsc(0:npgrd)
      real(kind=8):: pgrdox(0:npgrd),pgrdoy(0:npgrd) 
      real(kind=8):: pgrdph(0:npgrd),pgrdsp(0:npgrd),pgrdcp(0:npgrd)
      integer     :: pgrdid(0:npgrd)
      real(kind=8):: pgrddxi(pgrdns),pgrddyi(pgrdns),pgrddzi(pgrdns)
      integer     :: id(0:nzl)
      real(kind=8):: zlmin,zbeam,dzli

c  Sets electric field for a group of up to NPARPGRP particles from data
c  sets containing the potential on a 3-D grid.
c  Calling arguments:
c     nxyz         Product (1+nx)*(1+ny)*(1+nz)
c     pgrdns       number of data sets
c     pgrdp        arrays containing potential. In calling routine, it is
c                  dimensioned (0:nx,0:ny,0:nz,pgrdns), pgrdns different sets
c                  of fields
c                  Here, the first three dimensions are reduced to one.
c     np           number of particles
c     xp,yp,zp     coordinates of particles
c     npgrd        number of sets of potential arrays
c     pgrdxs,ys    mins of meshes of each set of potential arrays
c     cpgrdzs      z mins of data sets stored in the style of internal lattice
c                  arrays
c     pgrdsf       scale factor of each set of potential array
c     pgrdsc       scale factor of each set of potential array
c     pgrdox,oy    transverse offset of element
c     pgrdph       rotation angle of element 
c     pgrdsp,cp    sine and cosine of pgrdph 
c     pgrdid       index of array for set of potential arrays
c     pgrddxi,yz   one over the grid cell sizes of each set of potential arrays
c     pgrdnx,ny,nz Dimensions of arrays in calling routine 
c     ex,ey,ez     returned potentials for the particles
c     id           index to data sets stored in the style of internal lattice
c                  arrays
c     nzl          size of internal lattice arrays
c     zlmin        minimum location of internal lattice arrays relative to beam
c                  frame
c     zbeam        location of beam frame
c     dzli         one over internal lattice arrays grid cell size

c  Algorithm notes: potential array is dimensioned
c  (0:pgrdnx,0:pgrdny,0:pgrdnz,pgrdns)
c  outside, but is made two dimensional in this routine
c  so cell index into 1D array for vectorized deposition is:
c     i + j*(pgrdnx+1) + k*(pgrdnx+1)*(pgrdny+1)
c  Note that nnxy is added to the offset since the axial dimension is
c  (-1:pgrdnz+1)
c  The field is:
c     Ex = u0 * v0 * w0 * ex(i  ,j  ,k  ,pgrdid)
c        + u1 * v0 * w0 * ex(i+1,j  ,k  ,pgrdid)
c        + u0 * v1 * w0 * ex(i  ,j+1,k  ,pgrdid)
c        + ...
c
c  Note that this routine is very similar to sete3d routine.  It might be
c  a good idea to combine them to reduce code complexity, but that would make
c  sete3d more complex and possibly slower.  Since the pgrd routine will not
c  be used very often and sete3d is used always, the author felt it was better
c  to have a seperate routine for pgrd.

      integer     :: nnx,nnxy,ip,i,j,k,ind0,iz,ii,ib
      real(kind=8):: u0,u1,v0,v1,w0,w1,xfac,yfac,zfac,txp,typ,temp
      integer     :: indx05,indx06,indx07,indx11
      integer     :: indx17,indx18,indx19,indx23
      integer     :: noff(32)
      save noff

      nnx  = pgrdnx + 1
      nnxy = (pgrdnx + 1)*(pgrdny + 1)

c  Evaluation of offsets occurs on the first call only provided
c  local arrays are SAVE'd.

      if (noff(1) .ne. -nnxy) then
         noff(1)  = -nnxy               + nnxy
         noff(2)  = -nnxy   + 1         + nnxy
         noff(3)  = -nnxy   + nnx       + nnxy
         noff(4)  = -nnxy   + nnx   + 1 + nnxy
         noff(5)  =         - nnx       + nnxy
         noff(6)  =         - nnx   + 1 + nnxy
         noff(7)  =                 - 1 + nnxy
         noff(8)  =                   0 + nnxy
         noff(9)  =                   1 + nnxy
         noff(10) =                   2 + nnxy
         noff(11) =           nnx   - 1 + nnxy
         noff(12) =           nnx       + nnxy
         noff(13) =           nnx   + 1 + nnxy
         noff(14) =           nnx   + 2 + nnxy
         noff(15) =           2*nnx     + nnxy
         noff(16) =           2*nnx + 1 + nnxy
         noff(17) =  nnxy   - nnx       + nnxy
         noff(18) =  nnxy   - nnx   + 1 + nnxy
         noff(19) =  nnxy           - 1 + nnxy
         noff(20) =  nnxy               + nnxy
         noff(21) =  nnxy           + 1 + nnxy
         noff(22) =  nnxy           + 2 + nnxy
         noff(23) =  nnxy   + nnx   - 1 + nnxy
         noff(24) =  nnxy   + nnx       + nnxy
         noff(25) =  nnxy   + nnx   + 1 + nnxy
         noff(26) =  nnxy   + nnx   + 2 + nnxy
         noff(27) =  nnxy   + 2*nnx     + nnxy
         noff(28) =  nnxy   + 2*nnx + 1 + nnxy
         noff(29) =  2*nnxy             + nnxy
         noff(30) =  2*nnxy         + 1 + nnxy
         noff(31) =  2*nnxy + nnx       + nnxy
         noff(32) =  2*nnxy + nnx   + 1 + nnxy
      endif

c  Evaluation of E, vectorized over particles

c     --- All particles lie in the same z-cell
      iz = 0
      ii = id(iz)
      ib = pgrdid(ii)

c     --- Particle z-location in field grid and z linear weight factor 
      k  = (zp - cpgrdzs(iz))*pgrddzi(ib)
      w1 = (zp - cpgrdzs(iz))*pgrddzi(ib) - k
      w0 = 1. - w1

c     --- Exit field accumulation if particles are outside of axial grid
      if (k < 0 .or. k >= pgrdnz) return 

c     --- Scale factors 
      xfac = pgrdsc(ii) + pgrdsf(ii)
      yfac = pgrdsc(ii) + pgrdsf(ii)
      zfac = pgrdsc(ii) + pgrdsf(ii)

      do ip = 1, np

c        --- find transverse particle coordinate in frame of gridded field 
c            --- transverse offsets
         txp = xp(ip) - pgrdox(ii)
         typ = yp(ip) - pgrdoy(ii)
c            --- transverse rotation to take into account an active rotation 
c                of the field element.  Particles are rotated in that oposite 
c                sense of the element. Later the field components accumulated 
c                must be rotated back. 
         if ( pgrdph(ii) .ne. 0. ) then
           temp = txp  
           txp =  temp*pgrdcp(ii) + typ*pgrdsp(ii) 
           typ = -temp*pgrdsp(ii) + typ*pgrdcp(ii) 
         endif 

c        --- find location of particle in potential grid
         i =  abs(txp - pgrdxs(ii))*pgrddxi(ib)
         j =  abs(typ - pgrdys(ii))*pgrddyi(ib)
         ind0 = i + j*nnx + k*nnxy

c        --- Calculate linear weights
         u1 =  abs(txp - pgrdxs(ii))*pgrddxi(ib) - i
         v1 =  abs(typ - pgrdys(ii))*pgrddyi(ib) - j
         u0 = 1. - u1
         v0 = 1. - v1

c        --- Only calculate for particles inside the B field grid
         if (i >= 0 .and. i < pgrdnx .and.
     &       j >= 0 .and. j < pgrdny) then

           indx05 = ind0 + noff(5)
           indx06 = ind0 + noff(6)
           indx07 = ind0 + noff(7)
           indx11 = ind0 + noff(11)
           indx17 = ind0 + noff(17)
           indx18 = ind0 + noff(18)
           indx19 = ind0 + noff(19)
           indx23 = ind0 + noff(23)
c          --- make sure all indices refer to first quadrant
           if (i == 0) then
             indx07 = ind0 + 1
             indx11 = ind0 + nnx + 1
             indx19 = ind0 + nnxy + 1
             indx23 = ind0 + nnxy + nnx + 1
           endif
           if (j == 0) then
             indx05 = ind0 + nnx
             indx06 = ind0 + nnx + 1
             indx17 = ind0 + nnxy + nnx
             indx18 = ind0 + nnxy + nnx + 1
           endif
c          --- adjust sign of E field for approiate quadrant
           xfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddxi(ib)*0.5
           yfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddyi(ib)*0.5
           zfac = (pgrdsc(ii) + pgrdsf(ii))*pgrddzi(ib)*0.5
           if (xp(ip) < pgrdxs(ii)) xfac = -xfac
           if (yp(ip) < pgrdys(ii)) yfac = -yfac

           ex(ip)=ex(ip)+xfac*
     &           (u0*v0*w0*(pgrdp(indx07,ib)        - pgrdp(ind0+noff(09),ib))
     &          + u1*v0*w0*(pgrdp(ind0+noff(08),ib) - pgrdp(ind0+noff(10),ib))
     &          + u0*v1*w0*(pgrdp(indx11,ib)        - pgrdp(ind0+noff(13),ib))
     &          + u1*v1*w0*(pgrdp(ind0+noff(12),ib) - pgrdp(ind0+noff(14),ib))
     &          + u0*v0*w1*(pgrdp(indx19,ib)        - pgrdp(ind0+noff(21),ib))
     &          + u1*v0*w1*(pgrdp(ind0+noff(20),ib) - pgrdp(ind0+noff(22),ib))
     &          + u0*v1*w1*(pgrdp(indx23,ib)        - pgrdp(ind0+noff(25),ib))
     &          + u1*v1*w1*(pgrdp(ind0+noff(24),ib) - pgrdp(ind0+noff(26),ib)))

           ey(ip)=ey(ip)+yfac*
     &           (u0*v0*w0*(pgrdp(indx05,ib)        - pgrdp(ind0+noff(12),ib))
     &          + u1*v0*w0*(pgrdp(indx06,ib)        - pgrdp(ind0+noff(13),ib))
     &          + u0*v1*w0*(pgrdp(ind0+noff(08),ib) - pgrdp(ind0+noff(15),ib))
     &          + u1*v1*w0*(pgrdp(ind0+noff(09),ib) - pgrdp(ind0+noff(16),ib))
     &          + u0*v0*w1*(pgrdp(indx17,ib)        - pgrdp(ind0+noff(24),ib))
     &          + u1*v0*w1*(pgrdp(indx18,ib)        - pgrdp(ind0+noff(25),ib))
     &          + u0*v1*w1*(pgrdp(ind0+noff(20),ib) - pgrdp(ind0+noff(27),ib))
     &          + u1*v1*w1*(pgrdp(ind0+noff(21),ib) - pgrdp(ind0+noff(28),ib)))

           ez(ip)=ez(ip)+zfac*
     &           (u0*v0*w0*(pgrdp(ind0+noff(01),ib) - pgrdp(ind0+noff(20),ib))
     &          + u1*v0*w0*(pgrdp(ind0+noff(02),ib) - pgrdp(ind0+noff(21),ib))
     &          + u0*v1*w0*(pgrdp(ind0+noff(03),ib) - pgrdp(ind0+noff(24),ib))
     &          + u1*v1*w0*(pgrdp(ind0+noff(04),ib) - pgrdp(ind0+noff(25),ib))
     &          + u0*v0*w1*(pgrdp(ind0+noff(08),ib) - pgrdp(ind0+noff(29),ib))
     &          + u1*v0*w1*(pgrdp(ind0+noff(09),ib) - pgrdp(ind0+noff(30),ib))
     &          + u0*v1*w1*(pgrdp(ind0+noff(12),ib) - pgrdp(ind0+noff(31),ib))
     &          + u1*v1*w1*(pgrdp(ind0+noff(13),ib) - pgrdp(ind0+noff(32),ib)))

         endif

      enddo

      return
      end
c=========================================================================== 
      subroutine setrhoxy(rho1d,np,xp,yp,zp,zgrid,uzp,gaminv,q,wght)
      use InGen
      use InGenxy
      use InGen3d
      use InMesh3d
      use Fields3d
      use Picglb3d
      integer     :: np
      real(kind=8):: q,wght,zgrid
      real(kind=8):: rho1d(0:*)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np), gaminv(np)

c  Sets charge density.
c  When using thick slice version, makes call to 3d version of setrho.
c  When using slice version, puts all charge at iz=0.  The call to
c  the fixrhoxy routine is then used to distribute the charge density
c  throughout the rho array after all of the particles have been deposited.


c  Algorithm notes: rho array is dimensioned (0:nx,0:ny,0:nz) outside,
c  but is made one dimensional in this routine
c  so cell index into 1d rho array for vectorized deposition is:
c     i + j*(nx+1)
c  In each case,
c     rho(i  ,j  ,0  ) = rho(i  ,j  ,0  ) + u0*v0*g*vbeamfrm/(uzp*gaminv)
c     rho(i+1,j  ,0  ) = rho(i+1,j  ,0  ) + u1*v0*g*vbeamfrm/(uzp*gaminv)
c  Note that many changes are possible; for example, we might define
c  ind0(ir) and not use indx; this saves some store operations but
c  leads to a more complicated indirect address for the vectorized
c  gather-add-scatter loop.  It seems about 3% slower than the present way.
c  RHO must be zeroed in ZERORHO since it is not zeroed here (to allow
c  handling of blocks of particles at a time)

c     --- For vectorized algorithm
      integer     :: moff(0:3), indx(0:3,NPARPGRP)
c     --- For "scalar" (actually partly vectorized) algorithm
      integer     :: ii(0:NPARPGRP-1), jj(0:NPARPGRP-1)
c     --- For both algorithms
      real(kind=8):: s(0:3,NPARPGRP)

      integer     :: ip,i,j,ind0,m,ir
      real(kind=8):: g,dxi,dyi,u0,u1,v0,v1,gxfact,gyfact,f

      if (lthick) then
c       --- Call 3d version.
        call setrho3d(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,depos)
      else

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

c  Set up offset array for vectorized deposition:

      moff(0) = 0
      moff(1) = 1
      moff(2) = nx+1
      moff(3) = nx+2

c  Begin main loop over species, groups of NPARPGRP particles

      g = wght * q / (dx * dy)
      dxi = 1./dx
      dyi = 1./dy

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
c       --- special loop for 2-fold symmetry
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        g = g*0.5
        do ip = 1, np
           f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
           i  = (xp(ip) - xmmin) * dxi
           u1 = (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = abs(yp(ip))*dyi
           v1 = abs(yp(ip))*dyi - j
           v0 = 1. - v1
           ind0 = i + j*(nx+1)
           indx(0,ip) = ind0 + moff(0)
           indx(1,ip) = ind0 + moff(1)
           indx(2,ip) = ind0 + moff(2)
           indx(3,ip) = ind0 + moff(3)
           gyfact = 1.
           if (j == 0) gyfact = 2.
           s(0,ip) = u0 * v0 * f*gyfact
           s(1,ip) = u1 * v0 * f*gyfact
           s(2,ip) = u0 * v1 * f
           s(3,ip) = u1 * v1 * f
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        g = g*0.25
        do ip = 1, np
           f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
           i  = abs(xp(ip))*dxi
           u1 = abs(xp(ip))*dxi - i
           u0 = 1. - u1
           j  = abs(yp(ip))*dyi
           v1 = abs(yp(ip))*dyi - j
           v0 = 1. - v1
           ind0 = i + j*(nx+1)
           indx(0,ip) = ind0 + moff(0)
           indx(1,ip) = ind0 + moff(1)
           indx(2,ip) = ind0 + moff(2)
           indx(3,ip) = ind0 + moff(3)
           gxfact = 1.
           gyfact = 1.
           if (i == 0) gxfact = 2.
           if (j == 0) gyfact = 2.
           s(0,ip) = u0 * v0 * f*gxfact*gyfact
           s(1,ip) = u1 * v0 * f*gyfact
           s(2,ip) = u0 * v1 * f*gxfact
           s(3,ip) = u1 * v1 * f
        enddo
      else
c       --- normal loop
        do ip = 1, np
           f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
           i  = (xp(ip) - xmmin) * dxi
           u1 = (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = (yp(ip) - ymmin) * dyi
           v1 = (yp(ip) - ymmin) * dyi - j
           v0 = 1. - v1
           ind0 = i + j*(nx+1)
           indx(0,ip) = ind0 + moff(0)
           indx(1,ip) = ind0 + moff(1)
           indx(2,ip) = ind0 + moff(2)
           indx(3,ip) = ind0 + moff(3)
           s(0,ip) = u0 * v0 * f
           s(1,ip) = u1 * v0 * f
           s(2,ip) = u0 * v1 * f
           s(3,ip) = u1 * v1 * f
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
      do ip = 1, np
         if (uzp(ip) > 0) then
         do m = 0, 3
            rho1d(indx(m,ip)) = rho1d(indx(m,ip)) + s(m,ip)
         enddo
         endif
      enddo

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species, groups of NPARPGRP particles

      g = wght * q / (dx * dy)
      dxi = 1./dx
      dyi = 1./dy

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        g = g*0.5
        do ip = 1, np
          f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          ir = ip - 1
          ii(ir) = (xp(ip) - xmmin) * dxi
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = abs(yp(ip))*dyi
          v1     = abs(yp(ip))*dyi - jj(ir)
          v0     = 1. - v1
          gyfact = 1.
          if (jj(ir) == 0) gyfact = 2.
          s(0,ip) = u0 * v0 * f*gyfact
          s(1,ip) = u1 * v0 * f*gyfact
          s(2,ip) = u0 * v1 * f
          s(3,ip) = u1 * v1 * f
        enddo
      elseif (l4symtry) then
        g = g*0.25
        do ip = 1, np
          f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          ir = ip - 1
          ii(ir) = abs(xp(ip))*dxi
          u1     = abs(xp(ip))*dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = abs(yp(ip))*dyi
          v1     = abs(yp(ip))*dyi - jj(ir)
          v0     = 1. - v1
          gxfact = 1.
          gyfact = 1.
          if (ii(ir) == 0) gxfact = 2.
          if (jj(ir) == 0) gyfact = 2.
          s(0,ip) = u0 * v0 * f*gxfact*gyfact
          s(1,ip) = u1 * v0 * f*gyfact
          s(2,ip) = u0 * v1 * f*gxfact
          s(3,ip) = u1 * v1 * f
        enddo
      else
c       --- normal loop
        do ip = 1, np
          f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          ir = ip - 1
          ii(ir) = (xp(ip) - xmmin) * dxi
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = (yp(ip) - ymmin) * dyi
          v1     = (yp(ip) - ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          s(0,ip) = u0 * v0 * f
          s(1,ip) = u1 * v0 * f
          s(2,ip) = u0 * v1 * f
          s(3,ip) = u1 * v1 * f
        enddo
      endif
c     --- scalar loop does the actual deposition
      do ir = 0, np-1
         if (uzp(ir+1) > 0) then
         rho(ii(ir)  ,jj(ir)  ,0) = rho(ii(ir)  ,jj(ir)  ,0) + s(0,ir+1)
         rho(ii(ir)+1,jj(ir)  ,0) = rho(ii(ir)+1,jj(ir)  ,0) + s(1,ir+1)
         rho(ii(ir)  ,jj(ir)+1,0) = rho(ii(ir)  ,jj(ir)+1,0) + s(2,ir+1)
         rho(ii(ir)+1,jj(ir)+1,0) = rho(ii(ir)+1,jj(ir)+1,0) + s(3,ir+1)
        endif
      enddo

      endif

c     --- End of slice version
      endif
      return
      end
c=============================================================================
      subroutine loadrhoxy(ins_i,nps_i,is_i,lzero)
      use InGen
      use InGenxy
      use InPart
      use Picglb
      use Picglb3d
      use Particles
      use Particlesxy
      use InMesh3d
      use Fields3d
      integer     :: ins_i,nps_i,is_i
      logical     :: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- rho array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.

      integer     :: ins_u,nps_u
      integer     :: is1,is2
      integer     :: ip,ipmin,is

c     --- zero rho if requested
c     if (lzero) call zeroarry(rho,(nx+1)*(ny+1)*(nz+1))
      if (lzero) rho = 0.

c     --- set limits on loop over species
      if (is_i == -1) then
        is1 = 1
        is2 = ns
      else
        is1 = is_i
        is2 = is_i
      endif

c     --- set initial limits from input
c     --- (will be changed if necessary in the loop)
      ins_u = ins_i
      nps_u = nps_i

c     --- loop over species
      do is=is1,is2

c        --- get loop limits for particles if needed
         if (ins_i == -1) ins_u = ins(is)
         if (nps_i == -1) nps_u = nps(is)

c        --- loop over particle blocks
         do ipmin = ins_u, ins_u + nps_u - 1, NPARPGRP
            ip = min(NPARPGRP, ins_u+nps_u-ipmin)
            call setrhoxy(rho(0,0,0),ip,xp(ipmin),yp(ipmin),zp(ipmin),zgrid,
     &                    uzp(ipmin),gaminv(ipmin),sq(is),sw(is))
         enddo
      enddo

c     --- enforce axial periodicity if rho was zeroed
      if (lzero) call fixrhoxy(rho,nx,ny,nz,periinz,lthick)

      return
      end
c=============================================================================
      subroutine setcurrxy(curr,np,zp,uzp,gaminv,q,wght,zbeam,dzz,zzmin,dz,
     &                     lthick)
      integer     :: np
      real(kind=8):: q,wght,zbeam,dzz,zzmin,dz
      real(kind=8):: curr(0:*), zp(np), uzp(np), gaminv(np)
      logical     :: lthick

c  Sets 1d beam current directly from particle data.
c  If running slice version, puts all data at iz=0, else if running thick
c  slice version, makes a call to the 3d version of setcurr.

      real(kind=8):: csum,g
      integer     :: ip

      if (lthick) then
c       --- Call 3d version.
        call setcurr(curr,np,zp,uzp,gaminv,q,wght,zbeam,dzz,zzmin)

      else

c       --- Sum the charge times velocity over all particles.
        g = wght*q/dz
        csum = 0.
        do ip = 1,np
           csum = csum + gaminv(ip)*uzp(ip)
        enddo
        curr(0) = curr(0) + csum*g

      endif

      return
      end
c=============================================================================
      subroutine fieldsolxy(iwhich)
      use Constant
      use InGen
      use InGenxy
      use InGen3d
      use Picglb
      use Picglb3d
      use InMesh3d
      use Fields3d
      use LatticeInternal
      integer     :: iwhich

c If using the 3D field solver, calls the fieldsol routine in package w3d,
c otherwise only does a FFT field solve in 2D.


      integer     :: i,j,iz
      real(kind=8):: r,x,dxi,phiref
      character*(120):: outstr

c     --- Return if not doing any field solves.
      if ( fstype == -1) return

c     --- Call appropriate field solver
      if (lvp3d) then
c       --- 3D field solver
        call fieldsol3d(iwhich)
      else

c       --- 2D fieldsolver

        if (.not. linbend .or. iwhich > 0) then

c         --- if not in a bend, then make direct call to vpxy
          call copyarry(rho, phi(0,0,0), (nx+1)*(ny+1)*(nz+1))
          call vpxy(iwhich)

        else

c         --- otherwise, include bend correction terms

          dxi = 1./dx
          bndfit = 0
          bndferr = 2.*bndftol
          do while ((bndferr > bndftol) .and. (bndfit < bndfitmx))
            bndfit = bndfit + 1

c           --- save phi for error measure, can save it into scrtch
c           --- since that isn't used for much
            call copyarry(phi(0,0,0),scrtch(0,0),(nx+1)*(ny+1))

c           --- add curvature terms to the source (rho) as a correction
c           --- (remember that the previous phi is stored in scrtch)
            do i = 1, nx-1
              x = xmmin + i*dx
              r = 1./(rstar(0) + x)
              do j = 0, ny
                phi(i,j,0) = rho(i,j,0)*rstar(0)*r
     &                   + eps0*(scrtch(i+1,j) - scrtch(i-1,j))*0.5*dxi*r
              enddo
            enddo

c           --- call Cartesian field solver
            call vpxy (-1)

c           --- compute error
            bndferr = 0.
            phiref = 0.
            do i = 0, nx
               do j = 0, ny
                  bndferr = max( bndferr, abs(phi(i,j,0)-scrtch(i,j)) )
                  phiref = max( phiref, phi(i,j,0) )
               enddo
            enddo

c           --- for debug, print out the error
            if (bnprflag) then
              write (outstr,9985) it, bndfit, bndferr
 9985         format ("It =",i7," Bent field iteration",i3," Rel Change = ",
     &                1pe12.4)
              call remark(outstr)
            endif

          end do

c         --- if failure to converge, report the bad news to user
          if (bndferr > bndftol) then
            write (outstr,9995) bndferr, bndfit
 9995       format ("*** NONCONVERGENCE in bent field iteration",
     &              /,"Relative change = ",1pe12.4," after ",i3," iterations.")
            call remark(outstr)
          endif

        endif
        
c       --- and fill rest of phi array
        call copyarry(phi(0,0,0),phi(0,0,-1),(nx+1)*(ny+1))
        do iz=1,nz+1
          call copyarry(phi(0,0,0),phi(0,0,iz),(nx+1)*(ny+1))
        enddo

      endif

      return
      end
c=============================================================================
      subroutine vpxy(iwhich)
      use Constant
      use InGen
      use InGenxy
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Picglb
      use Fields3d

      integer     :: iwhich

c  Interface to VPOIS2D and VPOIS3D using variables from database of package 3D

      if (lvp3d) then
c       --- Use 3d field solver
        call vp3d(iwhich)
      else

c       --- Otherwise, use 2d field solver

        if (fstype == 0) then
          call vpois2d (iwhich, phi(0,0,0), phi(0,0,0), kxsq, kysq,
     &       attx, atty, filt, xmmax-xmmin, ymmax-ymmin,
     &       nx, ny, phi(0,0,1), 0, l2symtry,l4symtry)

        elseif (fstype == 1 .or. fstype == 2) then
c         --- Note that phi(,,-1) and phi(,,1) are passed in as
c         --- scratch space.
          call capmatxyf(iwhich,phi(0,0,0),kxsq,kysq,attx,atty,
     &             filt,xmmax-xmmin,ymmax-ymmin,nx,ny,dx,dy,xmmin,ymmin,
     &             phi(0,0,-1),phi(0,0,1),l2symtry,l4symtry)
        endif

      endif

      return
      end
c=============================================================================

