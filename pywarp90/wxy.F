#include "top.h"
c=============================================================================
c@(#) File WXY.M, version $Revision: 3.85 $, $Date: 2007/08/20 17:51:41 $
c# Copyright (c) 1990-1998, The Regents of the University of California.
c# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.
c  This is main file of package WXY of code WARP
c  XY electrostatic PIC code, Cartesian geometry, for beam problems
c  Alex Friedman, LLNL, (510)422-0827
c  David P. Grote, LLNL, (510)423-7194
c=============================================================================
      subroutine wxyinit
      use Subtimersxy
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c  Called at first reference to package (not nec. a "run" etc.).

      call wxyvers (STDOUT)

      if (lwxytimesubs) timewxyinit = timewxyinit + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine wxyvers (iout)
      use Subtimersxy
      use WXYversion
      integer(ISZ):: iout
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()
c  Echoes code version, etc. to output files as they're created
      call printpkgversion(iout,"Particle package WXY",verswxy)
      if (lwxytimesubs) timewxyvers = timewxyvers + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine wxygen()
      use Subtimersxy
      use GlobalVars
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
      use Fields3dParticles
      use Efields3d
      use Io
      use Lattice
      use LatticeInternal
      use Particles, Only: pgroup,npmax
      use Particlesxy
      use Picglb
      use Picglb3d
      use OutParams
      use Beam_acc
      use Z_arrays
      use Win_Moments
      use Z_Moments
      use Moments
      use Damped_eom
      use Hist
      use AMR
#ifdef MPIPARALLEL
      use Parallel
#endif

c  Invoked by the GENERATE command, it sets up the problem
c  This routine allots all of the neccesary dynamic arrays, calls the
c  particle loader and does the initial load onto the charge density
c  mesh, initializes arrays for the field solver and sets the mesh arrays,
c  does the initial field solve, and sets up other arrays that are needed.

      integer(ISZ):: i,j,k,is,ipmin,ip,iwin,nl
      integer(ISZ):: nextpid
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

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
      if (dt == 0. .and. vbeam .ne. 0.) dt = ds/vbeam

c  Set default values for the axial grid if using thin slice model.
c  Zmmin and zmmax are chosen so that dz = 1.
      if (.not. lthick) then
        if (zmminlocal == 0.0 .and. zmmaxlocal == 0.0) then
          zmminlocal = -0.5
          zmmaxlocal = -zmminlocal
        endif
        if (zimin == 0.0 .and. zimax == 0.0) then
          zimin = -0.5
          zimax = -zimin
        endif
        zmmin = zmminlocal
        zmmax = zmmaxlocal
        zmminglobal = zmmin
        zmmaxglobal = zmmax
        if (nz == 2) nz = 1
      endif

#ifdef MPIPARALLEL
c     --- Use routine from w3d to initialize the MPI and divide up
c     --- the problem. nz is set to nslaves so that the grid is divided
c     --- up properly among the processors. Note that then init_w3d_parallel
c     --- then divides nz by nslaves again for the decomposition.
      if (.not. lthick) then
        nz = nslaves
      endif
      nzfull = nz
      call init_w3d_parallel()
#else
      nzlocal = nz
      nzfull = nz
#endif

c  Call the species routine again so that the values of zimin_s and
c  zimax_s are set appropriately.  This resolves a circle where zimin is set
c  equal to zmminlocal which depends on vbeam which is calculated in derivqty
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
      call stepid (it,time,zbeam)

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
      zlmin = 0.
      zlmax = 0.

c  Set vbeamfrm, adjust x and ymmin for symmetries and set
c  so that by default izfsmax = nz (for full field solve).

      vbeamfrm = vbeam
      if (l2symtry) then
        ymmin = 0.
      elseif (l4symtry) then
        xmmin = 0.
        ymmin = 0.
      endif
      if (izfsmax == 0) izfsmax = nz

c  Setup arrays for potential and charge density for the particles.

      nxp = nx
      nyp = ny
#ifndef MPIPARALLEL
      nzp = nz
      zmminp = zmminlocal
      zmmaxp = zmmaxlocal
#endif

c  Create the dynamic arrays for fields

      nmxy  = max(nx,ny)
      nmxyz = max(nx,ny,nz)
      call gallot("Fields3d", 0)
      call setupSubcycling(pgroup)
      call setupSelfB(pgroup)
      call setupFields3dParticles()

c  This is a kludge for now to get the arrays setup.

      call assignrhoandphi(rhopndts(:,:,:,0,0),phipndts(:,:,:,0))
      call assignrhopandphip(rhopndts(:,:,:,0,0),phipndts(:,:,:,0))

c  Setup arrays for option which uses the selfe.

      if (ANY(efetch == 3) .and. fstype < 12) then
        nx_selfe = nxp
        ny_selfe = nyp
        nz_selfe = nzp
        call gallot("Efields3d",0)
      endif

c  Set the value of prwall, radius at which particles are lost
      if (prwall == 0.) prwall = xmmax - dx
      do k=0,nzzarr
        prwallz(k) = prwall
        if (prwallxz(k) == 0.) prwallxz(k) = prwallx
        if (prwallyz(k) == 0.) prwallyz(k) = prwally
        if (prwelipz(k) == 1.) prwelipz(k) = prwelip
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
      do k = 0, nzlocal
         zmeshlocal(k) = k * dz + zmminlocal
      enddo

c Calculate location of axis in mesh, the term dx*1.e-5 acts as fuzz
      ix_axis = nint(-xmmin/dx)
      iy_axis = nint(-ymmin/dy)
      iz_axis = nint(-zmminlocal/dz)

c Initialize grid when solvergeom=XYgeom
      if(AMRlevels>0) solvergeom=XYgeom
      if(solvergeom==XYgeom) then
        nzlocal = 0
        nz = nzlocal
        fstype = 10
        call init_base(nx,ny,dx,dy,xmmin,ymmin,.false.)
      end if

c  Create the dynamic arrays for particles (set npmax to an estimated 
c  length for now, for those loading schemes that don't actually
c  use a user-set npmax directly)

      call remark(" ---  Allocating space for particles")
      if (xrandom == "grid") npmax = nxstripe*nystripe*nzstripe
      if (xrandom == "fibonacc") npmax = nfibgrps*fibg1
      call alotlostpart

c  Create a spot in the pid array for the particle time step size

      if (dtpid == 0) dtpid = nextpid()

c  Load the particles, calculate the charge density

      call remark(" ---  Loading particles")
      cigarld = .false.
      cylinder = .true.
      call stptcl3d(pgroup)
      call setgamma(pgroup,lrelativ)
      call initdtp(pgroup)
      call initzpxy(pgroup)
      call remark(" ---  Setting charge density")
c     call zeroarry(rho,(nx+1)*(ny+1)*(nz+1))
c     call zeroarry(curr,nzzarr+1)
      call loadrhoxy(pgroup,-1,-1,-1,.true.)

c  Create the dynamic arrays for the partcle qtys needed for the
c  damped mover.  They are always allocated, but with length 1 if not used.

      if (eomdamp /= 0.) then
        if (exeomoldpid == 0) exeomoldpid = nextpid()
        if (eyeomoldpid == 0) eyeomoldpid = nextpid()
        if (ezeomoldpid == 0) ezeomoldpid = nextpid()
        if (exeomlagpid == 0) exeomlagpid = nextpid()
        if (eyeomlagpid == 0) eyeomlagpid = nextpid()
        if (ezeomlagpid == 0) ezeomlagpid = nextpid()
      endif

c  Create the dynamic arrays for "window" moments

      call remark(" ---  Allocating Win_Moments")
      zwindows(1,0) = zmminlocal
      zwindows(2,0) = zmmaxlocal
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
      if (zmmntmax == 0.) zmmntmax = zmmaxlocal
      if (zmmntmin == 0.) zmmntmin = zmminlocal
      dzm = (zmmntmax - zmmntmin)/nzmmnt
      dzmi = 1./dzm
      do k = 0, nzmmnt
         zmntmesh(k) = zmmntmin + k*dzm
      enddo

c  Create the dynamic arrays for lab frame moments

      call remark(" ---  Allocating Lab_Moments")
      call initlabwn(1)

c  Create the scratch arrays for phase space plots (permanent, for now)
c  and set limits for plots

c     call remark(" ---  Allocating scratch space for plots")
c     npsplt = min(nparpgrp-1, npmax)
c     if (npsplt == 0) npsplt = nparpgrp-1
      if (xplmin == 0.) xplmin = -xmmax
      if (xplmax == 0.) xplmax =  xmmax
      if (yplmin == 0.) yplmin = -ymmax
      if (yplmax == 0.) yplmax =  ymmax
      if (zplmin == 0.) zplmin =  zmminlocal
      if (zplmax == 0.) zplmax =  zmmaxlocal

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

      call prntpara(dx,dy,dz,lprntpara,pgroup)
      call prntpa3d(lprntpara)
c  Initial call to fieldsolver in order to initialize attx, kxsq, etc.
      call vpxy (1)

c  Initial fieldsolve, diagnostics
      call stepxy ("wxygen")

      if (lwxytimesubs) timewxygen = timewxygen + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine wxyexe()
      use Subtimersxy
      use Picglb
      use InGen
      use InGenxy
      use Picglb3d
      use Ctl_to_pic
      use Particles,Only: pgroup
      use Particlesxy

c  Takes a time step
c  This routine advances the mesh in the lab frame, sets the logicals
c  which control how this next step is to be done, and then calls
c  the routine STEPXY to do the step.

      real(kind=8):: vbeamfrm0,zbeam0
      real(kind=8):: ds0,zz,zcorrection
      integer(ISZ):: is,ip,i
      logical(ISZ):: lbendend
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

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
        do is=1,pgroup%ns
          do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
            pgroup%pid(ip,dtpid) = pgroup%pid(ip,dtpid)*dt/(ds0/vbeamfrm)
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
        if ((zbeam + vbeamfrm*dt + zcorrection - zbeam0) .ne. 0) then
          dt = dt*ds/(zbeam + vbeamfrm*dt + zcorrection - zbeam0)
        endif
      enddo
      zbeam = zbeam0

c  set timestep counter, time, and advance grid frame
c  Note that at the end of padvncxy, zbeam was set equal to zgrid.
c  With zgridprv=zbeam, the user only needs to set zbeam. zgrid and zgridprv
c  are then set from that.

      it = it + 1
      time = time + dt
      if (lgridqnt) then
c       --- zgrid is integer number of dz's
        zgridprv = zbeam
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

        do is=1,pgroup%ns
          do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
            pgroup%pid(ip,dtpid) = pgroup%pid(ip,dtpid)*(ds0/vbeamfrm - dt)/dt
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
        do is=1,pgroup%ns
          do ip=pgroup%ins(is),pgroup%ins(is)+pgroup%nps(is)-1
            pgroup%pid(ip,dtpid) = pgroup%pid(ip,dtpid)*ds0/vbeamfrm/dt
          enddo
        enddo
        ds = ds0
        dt = ds0/vbeamfrm
      endif

c  Have we reached the end of the run ?

c     if ( lfinishd ) then
c        call remark("wxyexe: problem completed.")
c        if (lwxytimesubs) timewxyexe = timewxyexe + wtime() - substarttime
c        return
c     elseif (nplive <= 0) then
c          call remark(" *** WXYEXE: stopping, nplive = 0")
c        if (lwxytimesubs) timewxyexe = timewxyexe + wtime() - substarttime
c          return
c     else
c        if (lwxytimesubs) timewxyexe = timewxyexe + wtime() - substarttime
c        return
c     endif

      if (lwxytimesubs) timewxyexe = timewxyexe + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine wxyfin()
      use Subtimersxy
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
      use Picglb
      use Picglb3d
      use Win_Moments
      use Z_Moments
      use Z_arrays
      use Hist

c  Finish up at end of RUN, or on receipt of FIN
c  This routine is never called, at present; history plots are
c  made using a Python interpreter script (histplot), and we just
c  end the run.  If we wanted to chain runs so that an output qty
c  might be plotted vs a parameter, this routine might be useful.
c  It would be needed for a non-Python version of WARP.
c  For now it serves as a place-holder.

      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()
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
      call gfree ("Z_arrays")

      if (lwxytimesubs) timewxyfin = timewxyfin + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine stepxy (caller)
      use Subtimersxy
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
      use Particles,Only: pgroup,nplive
      use Picglb
      use Picglb3d
      use Lattice
      use LatticeInternal
      use Timers
      use Z_Moments,Only: tempmaxp,tempminp,tempzmmnts0,tempzmmnts
      character(*):: caller
      logical(ISZ):: thisstep,thiszbeam,dolabwn

c  When called by WXYEXE, stepxy advances the system forward in time one 
c  timestep and gathers diagnostics.  When called by WXYGEN, stepxy takes 
c  a step of zero size, to compute fields, and gather diagnostics at start 
c  of run.


      real(kind=8):: zbeaml,zbeamr,timetemp,wtime
      real(kind=8):: substarttime
      if (lwxytimesubs) substarttime = wtime()

c --- Set the internal lattice variables. This is not generally necessary at
c --- this point (it is redundant most of the time, the next call to
c --- setlatt in this subroutine is sufficient). There are cases where
c --- this is required for consistency. Since it is cheap (time wise),
c --- it is better to make sure the data is consistent than to save a
c --- little bit of time. The value of nzl must be checked since other
c --- packages (like W3D) may have reset it. For example, if the
c --- W3D package is generated after the WXY package, nzl will be set to
c --- non-zero. Switching back to WXY and running step, the internal lattice
c --- would still be setup for the W3D package and so the step would produce
c --- erroneaous results.
      nzl = 0
      zlmin = 0.
      zlmax = 0.
      zlatbuffer = ds
      call setlatt

c  Main particle advance: x to t.l. it; v to t.l. it-1/2
c  Half-step in v from t.l. it-1   if last step was "special"
c  Full-step in v from t.l. it-3/2 if last step not "special"
c  No step at all if generating.

      if (caller == "wxyexe") then
        call padvncxy ("halfv",pgroup)
      endif

c  The next two variables are the left and right ends of the range centered
c  about the end of the current time step plus/minus one half a step.
c  The range is used is determining whether diagnostics are done which
c  are based on the z location of the beam frame.  The diagnostics are done
c  on the time step which ends closest to the value given in the controlling
c  arrays.
c  The absolute values are taken so that if dt < 0 or vbeamfrm < 0, then
c  it will still be true that zbeaml < zbeamr.
      zbeaml   = zbeam - abs(0.5*vbeamfrm*dt)
      zbeamr   = zbeam + abs(0.5*vbeamfrm*dt)

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
      lhist    = mod(it, nhist) .eq. 0
      ldump    = mod(it, itdump) .eq. 0
      llabwn   = dolabwn()
      lspecial = (lfinishd .or. lalways .or. lseldom .or. ldump .or.
     &            lmoments .or. lhist .or. llabwn .or. llast .or.
     &            (it .eq. 0) .or. allspecl)

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

c  Field-solve for potential 

      if (lbeforefs) call execuser("beforefs")
      call fieldsolxy(-1)
      if (lafterfs) call execuser("afterfs")

c  Pre-calculate the self-E if it is needed for sete3d. This is done
c  after the call to afterfs in case some manipulation is done to phi.
      if (ANY(efetch == 3))
     &  call getselfe3d(phi(0,0,-1),nx,ny,nzlocal,selfe(1,0,0,0),
     &                  nx_selfe,ny_selfe,nz_selfe,dx,dy,dz,pboundxy,.true.)

c  Initialize the moments arrays which are calculated during the synchv and
c  gen phases.
c  0. is passed in as a dummy for all of the particles coordinates
c  which are not used at this time.
      call getzmmnt(1,0.,0.,0.,0.,0.,0.,0.,
     &              0.,0.,0.,0.,0.,1,
     &              nplive,0.,0.,0.,1,-1,ns,
     &              tempmaxp,tempminp,tempzmmnts0,tempzmmnts)

c  If a flag was set making this a "special" step,
c  do a half-advance to bring v to t.l. it 

      if (caller == "wxyexe") then
         call padvncxy ("synchv",pgroup)
      elseif (caller == "wxygen") then
         call padvncxy ("gen",pgroup)
      endif

c  Finalize the moments calculation
      call getzmmnt(1,0.,0.,0.,0.,0.,0.,0.,
     &              0.,0.,0.,0.,0.,3,nplive,0.,0.,0.,
     &              1,0,ns,tempmaxp,tempminp,tempzmmnts0,tempzmmnts)

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

      if (lwxytimesubs) timestepxy = timestepxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine extebxy(np,xp,yp,zp,uzp,gaminv,dtl,dtr,
     &                   bx,by,bz,ex,ey,ez,m,q,bendres,bendradi,lexbend,
     &                   gammabar,zbeam,vbeam,dt,time)
      use Subtimersxy
      use Timers
      integer(ISZ):: np
      real(kind=8):: dtl,dtr,m,q,gammabar,zbeam,vbeam,dt,time
      real(kind=8):: xp(np), yp(np), zp(np)
      real(kind=8):: uzp(np), gaminv(np)
      real(kind=8):: bx(np), by(np), bz(np),bendres,bendradi
      real(kind=8):: ex(np), ey(np), ez(np)
      logical(ISZ):: lexbend


c  Calculates "external" E, B fields
c  Calculates electric or magnetic AG focusing fields, bending and dipole
c  fields, and accelerating fields.
c  Includes back-rotation associated with coordinate transformation into By. 

c  NOTE: When we (someday) set B_self from a Lorentz transformation
c  on E_self, we'll have to carefully work out a sequence of calls,
c  since this routine is called more than once in the PADVNCXY loop on
c  a single step at present.  Perhaps we will also have to compute
c  B_self more than once.

      integer(ISZ):: j,ip,in,iele,ii 
      real(kind=8):: timetemp,wtime
      real(kind=8):: substarttime
      if (lwxytimesubs) substarttime = wtime()

      timetemp = wtime()

c     --- handle quads
      call applyquad(np,xp,yp,1,zbeam,vbeam,1.,dtl,dtr,dt,.true.,ex,ey,bx,by)

c     --- handle dipos 
      call applydipo(np,1,zbeam,vbeam,1.,dtl,dtr,dt,.true.,ex,ey,bx,by)

c     --- handle sexts
      call applysext(np,xp,yp,1,zbeam,vbeam,1.,dtl,dtr,dt,.true.,ex,ey,bx,by)

c     --- handle hard-edge electric and magnetic multipoles
      call applyhele(np,xp,yp,1,zbeam,vbeam,1.,dtl,dtr,dt,.true.,
     &               ex,ey,ez,bx,by,bz)

c     --- fold in coordinate transformation associated with bends
      if (.not. lexbend) then
        call applybend(np,xp,uzp,1,bendres,bendradi,m,q,.true.,by)
      endif

c     --- uniform fields
      call applyuniformfields(np,ex,ey,ez,bx,by,bz)

c     --- electrostatic multipole components
      call applyemlt(np,xp,yp,1,zbeam,dtl,dtr,dt,.true.,ex,ey,ez)

c     --- magnetostatic multipole components
      call applymmlt(np,xp,yp,1,zbeam,dtl,dtr,dt,.true.,bx,by,bz)

c     --- handle electric fields from 3-D grid
      call applyegrd(np,xp,yp,1,zbeam,.true.,ex,ey,ez)

c     --- handle magnetic fields from 3-D grid
      call applybgrd(np,xp,yp,1,zbeam,.true.,bx,by,bz)

c     --- handle electrostatic potential from 3-D grid
      call applypgrd(np,xp,yp,1,zbeam,.true.,ex,ey,ez)

c     --- Set E to zero for lost particles
      do ip=1,np
        if (uzp(ip) == 0) then
          ex(ip) = 0.
          ey(ip) = 0.
          ez(ip) = 0.
        endif
      enddo

c     --- Accumulate time for applying fields from the lattice
      latticetime = latticetime + (wtime() - timetemp)

      if (lwxytimesubs) timeextebxy = timeextebxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine otherexy (np,xp,yp,dedr,dexdx,deydy,dbdr,ex,ey,ez,bx,by,bz)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: dedr,dexdx,deydy,dbdr
      real(kind=8):: xp(np), yp(np)
      real(kind=8):: ex(np), ey(np), ez(np)
      real(kind=8):: bx(np), by(np), bz(np)

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

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

      if (lwxytimesubs) timeotherexy = timeotherexy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine fetchexy(pgroup,ipmin,ip,is,ex,ey,ez)
      use ParticleGroupmodule
      use GlobalVars
      use Subtimersxy
      use Picglb
      use Picglb3d
      use InPart,Only: efetch
      use InGenxy
      use InGen3d
      use InMesh3d
      use Fields3d
      use Efields3d
      use FieldSolveAPI
      type(ParticleGroup):: pgroup
      integer(ISZ):: ipmin,ip,is
      real(kind=8), dimension(ip) :: ex,ey,ez,zpo

      if (lresetparticlee) then
        ex = 0.
        ey = 0.
        ez = 0.
      endif

      if(solvergeom/=XYgeom) then
c       --- Obtain the self-field from the electrostatic potential
        if (lwithez) then
c         --- Note that zp, zgrid, and zmminlocal are all passed in as zero
c         --- to avoid any possible funny business with rounding.
c         --- ds is passed in for dz since that is the distance between
c         --- the locations where phi is known.
          zpo = 0.
          call sete3d(phi,selfe,ip,
     &                pgroup%xp(ipmin),pgroup%yp(ipmin),zpo,
     &                0.,xmmin,ymmin,0.,dx,dy,ds,nx,ny,nzlocal,
     &                efetch(is),ex,ey,ez,l2symtry,l4symtry,.false.,0,0,1)
c         --- The zp is preserved so that the small correction done
c         --- in setdtp is still used.
          pgroup%zp(ipmin:ipmin+ip-1) = zpo
        else
          call sete3d(phi,selfe,ip,
     &                pgroup%xp(ipmin),pgroup%yp(ipmin),pgroup%zp(ipmin),
     &                zgridprv,xmmin,ymmin,zmminlocal,dx,dy,dz,nx,ny,nzlocal,
     &                efetch(is),ex,ey,ez,l2symtry,l4symtry,.false.,0,0,1)
        endif
      elseif(solvergeom==XYgeom) then
        call fieldweightxz(pgroup%xp(ipmin),pgroup%yp(ipmin),ex,ey,ip,zgridprv,
     &                     efetch(is))
      end if
      return
      end
c=============================================================================
      subroutine padvncxy(center,pgroup)
      use ParticleGroupmodule
      use Subtimersxy
      use GlobalVars
      use InMesh3d
      use InGen
      use InGen3d
      use InGenxy
      use InPart
      use InGaps
      use InDiag
      use Lattice
      use LatticeInternal
      use Particles,Only: nplive,wpid
      use Particlesxy
      use Fields3d
      use Efields3d
      use Picglb
      use Picglb3d
      use Beam_acc
      use Z_arrays
      use Z_Moments, only: nzmmnt,nszmmnt,tempmaxp,tempminp,tempzmmnts0,tempzmmnts
      use Damped_eom
      character(*):: center
      type(ParticleGroup):: pgroup

c  Advances the particles position and velocity according to CENTER,
c  and also loads RHO at the new time level.

c     --- Create local pointers to the arrays in pgroup.
      real(kind=8),pointer:: xp(:),yp(:),zp(:),uxp(:),uyp(:),uzp(:)
      real(kind=8),pointer:: ex(:),ey(:),ez(:)
      real(kind=8),pointer:: bx(:),by(:),bz(:)
      real(kind=8),pointer:: gaminv(:),pid(:,:)
      real(kind=8),pointer:: sm(:),sq(:),sw(:),dtscale(:)
      integer(ISZ),pointer:: ins(:),nps(:)

      integer(ISZ):: isid,is,ismax,ip,ipmin,i
      real(kind=8):: uxpadv,uypadv,uzpadv
      real(kind=8):: ezo(nparpgrp)
      real(kind=8):: xpo(nparpgrp), ypo(nparpgrp), zpo(nparpgrp)
      real(kind=8):: uxpo(nparpgrp), uypo(nparpgrp), uzpo(nparpgrp)
      real(kind=8):: gaminvo(nparpgrp)
      real(kind=8):: bendres, bendradi
#ifdef _OPENMP
      real(kind=8),allocatable:: threadmaxp(:,:,:),threadminp(:,:,:)
      real(kind=8),allocatable:: threadzmmnts0(:,:,:),threadzmmnts(:,:,:,:)
      integer(ISZ):: ithread,omp_get_thread_num
      integer(ISZ):: nthread,omp_get_num_threads
      integer(ISZ):: allocerror
#endif
      real(kind=8),pointer:: maxp(:,:),minp(:,:)
      real(kind=8),pointer:: zmmnts0(:,:),zmmnts(:,:,:)
      integer(ISZ):: iter_dt,niter_dt_tmp
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c     --- Create local pointers to the arrays in pgroup.
      xp => pgroup%xp
      yp => pgroup%yp
      zp => pgroup%zp
      uxp => pgroup%uxp
      uyp => pgroup%uyp
      uzp => pgroup%uzp
      gaminv => pgroup%gaminv
      ex => pgroup%ex
      ey => pgroup%ey
      ez => pgroup%ez
      bx => pgroup%bx
      by => pgroup%by
      bz => pgroup%bz
      pid => pgroup%pid

      sm => pgroup%sm
      sq => pgroup%sq
      sw => pgroup%sw
      ins => pgroup%ins
      nps => pgroup%nps
      dtscale => pgroup%dtscale

      ismax = maxval(pgroup%sid)+1
      if (ismax == 0) return

c  Zero curr if center is synchv or if using gaps
c     if (center == "synchv" .or. ifgap) call zeroarry(curr,nzzarr+1)
      if (center == "synchv" .or. ifgap) then
c       --- First set nszarr and update array allocation if needed
        if (lspeciesmoments) then
c         --- Check if the moments are to be calculated separately for
c         --- each species. If so, check if nszarr already has been set
c         --- appropriately. If not, set it and allocate the arrays.
c         --- If only one species, then don't have separate species data.
          if (nszarr /= ns .and. ns > 1) then
            nszarr = ns
            call gchange("Z_arrays",0)
          endif
        else
          if (nszarr /= 0) then
            nszarr = 0
            call gchange("Z_arrays",0)
          endif
        endif
        curr = 0.
      endif

c  Zero the moments if center is synchv
      if (ldiag .and. (center == "synchv" .or. center == "gen")) then
#ifdef _OPENMP
        allocate(maxp(6,0:nszmmnt),minp(6,0:nszmmnt),
     &           zmmnts0(NUMZMMNT,0:nszmmnt),
     &           zmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt),stat=allocerror)
        if (allocerror /= 0) then
          print*,"padvncxy: allocation error ",allocerror,
     &           ": could not allocate temp arrays to shape ",nszmmnt
          call kaboom("padvncxy: allocation error")
          return
        endif
        maxp = -LARGEPOS
        minp = +LARGEPOS
        zmmnts0 = 0.
        zmmnts = 0.
#else
        maxp => tempmaxp
        minp => tempminp
        zmmnts0 => tempzmmnts0
        zmmnts => tempzmmnts
#endif
      endif

!$OMP PARALLEL
!$OMP&PRIVATE(ip,ipmin,i)
!$OMP&PRIVATE(xpo,ypo,zpo,uxpo,uypo,uzpo,gaminvo,ezo,iter_dt)
!$OMP&PRIVATE(uxpadv,uypadv,uzpadv,bendres,bendradi)
!$OMP&FIRSTPRIVATE(maxp,minp,zmmnts0,zmmnts)

c  Zero the bend radius and residence fraction
      bendradi = 0.
      bendres = 0.

      if(solvergeom==XYgeom) ez = 0.

c  Loop over species
      do is=1,pgroup%ns
        isid = pgroup%sid(is-1) + 1

c  Loop over particle blocks; move each block separately
!$OMP DO
        do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
          ip = min(nparpgrp, ins(is)+nps(is)-ipmin)

          call fetchexy(pgroup,ipmin,ip,is,ex(ipmin),ey(ipmin),ez(ipmin))

c         --- Scale the self E-field to get the lowest order relativistic
c         --- correction.
          if (relativity == 1)
     &      call sete3d_relativity(ip,ex(ipmin),ey(ipmin),vbeam)

c         --- Zero out B field arrays
          bx(ipmin:ipmin+ip-1) = 0.
          by(ipmin:ipmin+ip-1) = 0.
          bz(ipmin:ipmin+ip-1) = 0.

c         --- Compute lag average for experimental damping algorithm
          if (eomdamp .ne. 0.)
     &      call edamp(eomdamp,it,itdamp,center,ip,
     &                 pid(ipmin,exeomoldpid),
     &                 pid(ipmin,eyeomoldpid),
     &                 pid(ipmin,ezeomoldpid),
     &                 pid(ipmin,exeomlagpid),
     &                 pid(ipmin,eyeomlagpid),
     &                 pid(ipmin,ezeomlagpid))

c         --- HALFV
          if (center == "halfv") then
c           --- Obtain bend radii and residence factors
            call getbend(1,1,zbeam,vbeam,1.,bendres,bendradi,0.,dt,.true.)
c           --- Correct Ez_self for warped mesh effect
            call bendezxy (ip,xp(ipmin),zp(ipmin),ez(ipmin),
     &                     bendres,bendradi,bends,bnezflag)
c           --- Add in Ez from axially-smoothed gaps 
            call gapfield (ip,zp(ipmin),ez(ipmin),zbeam,zzmin,egap,dzz)
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin))
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),0.,dt*0.5,bx(ipmin),by(ipmin),bz(ipmin),
     &                    ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)

            if (lvzchang) then
              if (lexactds .and. .not. linbend .and. .not. linaccl(0) .and.
     &            (ibpush == 0 .or.
     &             (maxval(abs(bx(ipmin:ipmin+ip-1))) == 0. .and.
     &              maxval(abs(by(ipmin:ipmin+ip-1))) == 0.))) then
                niter_dt_tmp = 1
c               --- Calculate new step size directly.
                call getnewdtpwithe(ip,pid(ipmin,dtpid),
     &                              uzp(ipmin),gaminv(ipmin),
     &                              ez(ipmin),sq(is),sm(is),ds)
              else
                niter_dt_tmp = niter_dt
c               --- Save current position and velocity and ez
                do i=1,ip
                  xpo(i) = xp(ipmin+i-1)
                  ypo(i) = yp(ipmin+i-1)
                  zpo(i) = zp(ipmin+i-1)
                  uxpo(i) = uxp(ipmin+i-1)
                  uypo(i) = uyp(ipmin+i-1)
                  uzpo(i) = uzp(ipmin+i-1)
                  gaminvo(i) = gaminv(ipmin+i-1)
                  ezo(i) = ez(ipmin+i-1)
                enddo
              endif
            endif


c           --- Iterate over advance to calculate dt
c           --- (Only done when lvzchang is true.)
            do iter_dt=1,niter_dt_tmp
c             --- Get field from accelerating elements (depends on dtp)
              call applyacclxy(ip,xp(ipmin),zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                         pid(ipmin,dtpid),0.,dt*0.5,
     &                         sm(is),sq(is),dt,
     &                         .true.,ez(ipmin))
c             --- Correction to z on entry/exit to accelerator gap
              call zgapcorrxy(ip,zp(ipmin),xp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                        pid(ipmin,dtpid),0., dt*0.5, dt,
     &                        sm(1), sq(1), time)
c             --- Magnetic field increment to momenta
              call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is),
     &                      pid(ipmin,dtpid),0.5,ibpush)
c             --- Final half-electric field increment to momenta
              call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                      ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is),
     &                      pid(ipmin,dtpid), 0.5)
c             --- Advance relativistic Gamma factor
              call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                      gamadv,lrelativ)
c             --- Position advance
              call xpushxy (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                      uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      pid(ipmin,dtpid))
c             --- Simplified translation of position for warped mesh effect
              call bendcorxy(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       gaminv(ipmin),pid(ipmin,dtpid),
     &                       zpo,ds,bendres,bendradi,
     &                       lexbend)

              if (lvzchang .and. iter_dt < niter_dt_tmp) then
c               --- Calculate dt for particles.
                call setdtp(ip,pid(ipmin,dtpid),
     &                      xp(ipmin),zp(ipmin),zpo,uzp(ipmin),
     &                      bendres,bendradi,zgridprv)
c               --- Restore position and velocity if not finished with iteration
                do i=1,ip
                  xp(ipmin+i-1) = xpo(i)
                  yp(ipmin+i-1) = ypo(i)
                  zp(ipmin+i-1) = zpo(i)
                  uxp(ipmin+i-1) = uxpo(i)
                  uyp(ipmin+i-1) = uypo(i)
                  uzp(ipmin+i-1) = uzpo(i)
                  gaminv(ipmin+i-1) = gaminvo(i)
                  ez(ipmin+i-1) = ezo(i)
                enddo
              endif
            enddo

c           --- If using a common z position, enforce it to remove small
c           --- differences. The iteration done
c           --- Exact translation of position and velocity for warped
c           --- mesh effect
            call exbendcorxy(ip,xp(ipmin),zp(ipmin),uxp(ipmin),uzp(ipmin),
     &                       zpo,ds,bendres,bendradi,lexbend)
c           --- Correct position advance for slanted dipole entry/exit
            call sledgcor(pgroup,ip,ipmin,zpo,zbeam+vbeam*dt,zbeam,vbeam,
     &                    pgroup%sm(is),pgroup%sq(is),.true.)
c         --- PUSHV
          elseif (center == "pushv") then
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin))
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),-SMALLPOS,SMALLPOS,bx(ipmin),by(ipmin),bz(ipmin),
     &                    ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid),0.5,ibpush)
c           --- Final half-electric field increment to momenta
            call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                    ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid), 1.0)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid),0.5,ibpush)

c         --- PUSHV
          elseif (center == "lpushv") then
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin))
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),-dt*0.5,0.,bx(ipmin),by(ipmin),bz(ipmin),
     &                    ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)
c           --- Final half-electric field increment to momenta
            call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                    ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid), 0.5)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid),0.5,ibpush)

c         --- PUSHV
          elseif (center == "rpushv") then
c           --- Add in ears and uniform focusing E field pieces
            call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                     ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin))
c           --- Set quad, dipole E and B;  All: Bz
            call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                    gaminv(ipmin),0.,dt*0.5,bx(ipmin),by(ipmin),bz(ipmin),
     &                    ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),bendres,bendradi,lexbend,
     &                    gammabar,zbeam,vbeam,dt,time)
c           --- Magnetic field increment to momenta
            call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                    bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid),0.5,ibpush)
c           --- Final half-electric field increment to momenta
            call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                    ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is),
     &                    pid(ipmin,dtpid), 0.5)
c           --- Advance relativistic Gamma factor
            call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                    gamadv,lrelativ)

c         --- PUSHX
          elseif (center == "pushx") then
c             --- Position advance
              call xpushxy (ip, xp(ipmin), yp(ipmin), zp(ipmin), uxp(ipmin),
     &                      uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                      pid(ipmin,dtpid))

c         --- SYNCHV or GEN
          elseif (center == "synchv" .or. center == "gen") then
c            --- Copy 'old' velocity into uxpo, uypo, and uzpo
             do i=1,ip
                uxpo(i) = uxp(ipmin+i-1)
                uypo(i) = uyp(ipmin+i-1)
                uzpo(i) = uzp(ipmin+i-1)
             enddo
c            --- Obtain bend radii and residence factors
             call getbend(1,1,zbeam,vbeam,1.,bendres,bendradi,-dt,0.,.true.)
c            --- Correct Ez_self for warped mesh effect
             call bendezxy (ip,xp(ipmin),zp(ipmin),ez(ipmin),
     &                      bendres,bendradi,bends,bnezflag)
c            --- Add in Ez from axially-smoothed gaps 
             call gapfield (ip,zp(ipmin),ez(ipmin),zbeam,zzmin,egap,dzz)
c            --- Add in ears and uniform focusing E field pieces
             call otherexy (ip,xp(ipmin),yp(ipmin),dedr,dexdx,deydy,dbdr,
     &                      ex(ipmin),ey(ipmin),ez(ipmin),bx(ipmin),by(ipmin),bz(ipmin))
c            --- Set quad, dipole E and B; All: Bz
             call extebxy (ip,xp(ipmin),yp(ipmin),zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),-dt*0.5,0.,
     &                     bx(ipmin),by(ipmin),bz(ipmin),ex(ipmin),ey(ipmin),ez(ipmin),sm(is),sq(is),
     &                     bendres,bendradi,lexbend,
     &                     gammabar,zbeam,vbeam,dt,time)
c            --- Get field from accelerating elements (depends on dtp)
             call applyacclxy(ip,xp(ipmin),zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                        pid(ipmin,dtpid),-dt*0.5,0.,
     &                        sm(is),sq(is),dt,
     &                        .true.,ez(ipmin))
c            --- Half electric field increment to momenta
             call epushxy (ip, uxp(ipmin), uyp(ipmin), uzp(ipmin),
     &                     ex(ipmin), ey(ipmin), ez(ipmin), sq(is), sm(is),
     &                     pid(ipmin,dtpid), 0.5)
c            --- Advance relativistic Gamma factor
             call gammaadv(ip,gaminv(ipmin),uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                     gamadv,lrelativ)
c            --- Half magnetic field increment to momenta
             call bpushxy (ip,uxp(ipmin),uyp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                     bx(ipmin), by(ipmin), bz(ipmin), sq(is), sm(is),
     &                     pid(ipmin,dtpid),0.5,ibpush)
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
             if (ldiag) then
               if(wpid==0) then
                 call getzmmnt(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                         uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                         gaminv(ipmin),sq(is),sm(is),sw(is),
     &                         dt*0.5,dtscale(is),2,nplive,uxpo,uypo,uzpo,
     &                         is,isid,ismax,
     &                         maxp,minp,zmmnts0,zmmnts)
               else
                 call getzmmnt_weights(ip,xp(ipmin),yp(ipmin),zp(ipmin),
     &                         uxp(ipmin),uyp(ipmin),uzp(ipmin),
     &                         gaminv(ipmin),pid(ipmin,wpid),
     &                         sq(1),sm(is),sw(is),dt*0.5,dtscale(is),
     &                         2,nplive,uxpo,uypo,uzpo,
     &                         is,isid,ismax,
     &                         maxp,minp,zmmnts0,zmmnts)
               endif
             endif
c            --- Calculate current
             call setcurrxy(nzzarr,nszarr,curr,ip,zp(ipmin),uzp(ipmin),
     &                      gaminv(ipmin),sq(is),sw(is),is,zbeam,dzz,zzmin,
     &                      ns,lspeciesmoments,dz,lthick)
          endif
        enddo
!$OMP END DO

      enddo
c     --- End loop over species

#ifdef _OPENMP
      if (center == "synchv" .or. center == "gen") then
        ithread = omp_get_thread_num() + 1
        nthread = omp_get_num_threads()
!$OMP SINGLE
        allocate(threadmaxp(6,0:nszmmnt,nthread),
     &           threadminp(6,0:nszmmnt,nthread),
     &           threadzmmnts0(NUMZMMNT,0:nszmmnt,nthread),
     &           threadzmmnts(0:nzmmnt,NUMZMMNT,0:nszmmnt,nthread),
     &           stat=allocerror)
        if (allocerror /= 0) then
          print*,"padvncxy: allocation error ",allocerror,
     &           ": could not allocate temp arrays to shape ",nszmmnt,nthread
          call kaboom("padvncxy: allocation error")
          return
        endif
!$OMP END SINGLE
        threadmaxp(:,:,ithread) = maxp
        threadminp(:,:,ithread) = minp
        threadzmmnts0(:,:,ithread) = zmmnts0
        threadzmmnts(:,:,:,ithread) = zmmnts
        deallocate(maxp,minp,zmmnts0,zmmnts)
      endif
#endif
!$OMP END PARALLEL

      if (center == "halfv") then
        do is=1,ns

c         --- Apply particle boundary conditions

          if (lthick) then
            if (periinz)
     &         call periz(nps(is),zp(ins(is)),zgrid,zmmaxlocal,zmminlocal)
            if (stickyz)
     &         call stckyz(nps(is),zp(ins(is)),zmmaxlocal,zmminlocal,dz,uxp(ins(is)),
     &                     uyp(ins(is)),uzp(ins(is)),zgrid)
          endif

          if (stickyxy)
     &       call stckxy3d(nps(is),xp(ins(is)),xmmax,xmmin,dx,
     &                     yp(ins(is)),ymmax,ymmin,dy,
     &                     zp(ins(is)),zmminlocal,dz,uxp(ins(is)),
     &                     uyp(ins(is)),uzp(ins(is)),gaminv(ins(is)),
     &                     zgrid,zbeam,l2symtry,l4symtry,pboundxy,.true.)
        enddo

        if (lcallscraper) call execuser("callscraper")

        do is=1,ns
          call processlostpart(pgroup,is,clearlostpart,time,zbeam)
        enddo

c       --- Collect charge density (xy) and current (1d)
        call loadrhoxy(pgroup,-1,-1,-1,.true.)

        if (ifgap) then
          do ipmin = ins(is), ins(is) + nps(is) - 1, nparpgrp
            ip = min(nparpgrp, ins(is)+nps(is)-ipmin)
            call setcurrxy(nzzarr,nszarr,curr,ip,zp(ipmin),uzp(ipmin),
     &                     gaminv(ipmin),sq(is),sw(is),is,zbeam,dzz,zzmin,
     &                     ns,lspeciesmoments,dz,lthick)
          enddo
          call fixcurrxy(curr,nzzarr,nszarr,periinz,lthick)
        endif

c       --- Advance beam frame location using the nominal beam frame velocity.
        zbeam = zgrid

c       --- zgridprv needs to be updated for the "synchv" step
c       --- Note that zgridprv is also set at the beginning of wxyexe
        zgridprv = zbeam

      endif

c     ---  Do final stuff for moments calculation and fix current
      if (center == "synchv" .or. center == "gen") then

        if (ldiag) then
#ifdef _OPENMP
          tempmaxp = max(tempmaxp,maxval(threadmaxp(:,1:nthread),3))
          tempminp = min(tempminp,minval(threadminp(:,1:nthread),3))
          tempzmmnts0 = tempzmmnts0 + sum(threadzmmnts0(:,1:nthread),3)
          tempzmmnts = tempzmmnts + sum(threadzmmnts(:,:,1:nthread),4)
          deallocate(threadmaxp,threadminp,threadzmmnts0,threadzmmnts)
#endif

          call getlabwn()
c         --- load rho diagnostic qtys
          call rhodia3d()
        endif

c       --- Fix the current
        call fixcurrxy(curr,nzzarr,nszarr,periinz,lthick)

      endif

      if (lwxytimesubs) timepadvncxy = timepadvncxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine fixrhoxy(rho,nx,ny,nzlocal,periinz,lthick)
      use Subtimersxy
      integer(ISZ):: nx,ny,nzlocal
      real(kind=8):: rho(0:nx,0:ny,0:nzlocal)
      logical(ISZ):: periinz,lthick

c  Sums the first and last slices of rho for periodicity
c  and puts the result into both slices.
c  For thin slice model, average out rho along the axis.

      integer(ISZ):: ix,iy,iz
      real(kind=8):: nzi
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      if(nzlocal==0) return

      if (periinz) then

!$OMP DO
        do ix=0,nx
          do iy=0,ny
            rho(ix,iy,0)  = rho(ix,iy,0) + rho(ix,iy,nzlocal)
            rho(ix,iy,nzlocal) = rho(ix,iy,0)
          enddo
        enddo
!$OMP END DO

      endif

#ifdef MPIPARALLEL
c     --- Sum up rho over the processors.
      if (.not. lthick) then
        call parallelsumrealarray(rho(0,0,0),(1+nx)*(1+ny))
      endif
#endif

      if (.not. lthick) then
c       --- In the setrhoxy routine, all of the charge was put at iz=0.

        if (nzlocal > 1) then
          nzi = 1./real(nzlocal)
!$OMP DO
          do iy=0,ny
            do ix=0,nx
              rho(ix,iy,0) = rho(ix,iy,0)*nzi
            enddo
          enddo
!$OMP END DO
        endif

c       --- Distribute rho to all z planes.
!$OMP DO
        do iz=1,nzlocal
          do iy=0,ny
            do ix=0,nx
              rho(ix,iy,iz) = rho(ix,iy,0)
            enddo
          enddo
        enddo
!$OMP END DO

      endif

      if (lwxytimesubs) timefixrhoxy = timefixrhoxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine fixcurrxy(curr,nzzarr,nszarr,periinz,lthick)
      use Subtimersxy
      integer(ISZ):: nzzarr,nszarr
      real(kind=8):: curr(0:nzzarr,0:nszarr)
      logical(ISZ):: periinz,lthick

c  Sums the first and last slices of curr for periodicity
c  and puts the result into both slices.
c  For thin slice model, distribute the current along the axis.

      integer(ISZ):: iz
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      if (nszarr > 0) then
c       --- Sum up the current from all of the speices.
        curr(:,nszarr) = sum(curr(:,0:nszarr-1),2)
      endif

      if (periinz) then

        curr(0,:)  = curr(0,:) + curr(nzzarr,:)
        curr(nzzarr,:) = curr(0,:)

      endif

      if (.not. lthick .and. nzzarr > 1) then
c       --- In the setcurrxy routine, all of the current was put at iz=0.

c       --- Distribute the current.
        do iz=1,nzzarr
          curr(iz,:) = curr(0,:)
        enddo

      endif

      if (lwxytimesubs) timefixcurrxy = timefixcurrxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine epushxy(np,uxp,uyp,uzp,ex,ey,ez,q,m,dtp,fdt)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: uxp(np),uyp(np),uzp(np)
      real(kind=8):: ex(np),ey(np),ez(np)
      real(kind=8):: q,m,dtp(np),fdt

c  Push the particle velocity with E field

      integer(ISZ):: ip
      real(kind=8):: const
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()
      const = q*fdt/m

      do ip=1,np
        uxp(ip) = uxp(ip) + ex(ip)*const*dtp(ip)
        uyp(ip) = uyp(ip) + ey(ip)*const*dtp(ip)
        uzp(ip) = uzp(ip) + ez(ip)*const*dtp(ip)
      enddo

      if (lwxytimesubs) timeepushxy = timeepushxy + wtime() - substarttime
      return
      end
c============================================================================= 
      subroutine bpushxy(np,uxp,uyp,uzp,gaminv,bx,by,bz,q,m,dtp,fdt,ibpush)
      use Subtimersxy
      integer(ISZ):: np,ibpush
      real(kind=8):: uxp(np),uyp(np),uzp(np),gaminv(np)
      real(kind=8):: bx(np),by(np),bz(np)
      real(kind=8):: q,m,dtp(np),fdt

c  Push the particle velocity with B field

      integer(ISZ):: ip
      real(kind=8):: btot,btotinv,tanalpha
      real(kind=8):: tx,ty,tz,tsqi,sx,sy,sz,uxppr,uyppr,uzppr
      real(kind=8):: const
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()
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
            btot = sqrt(bx(ip)**2 + by(ip)**2 + bz(ip)**2)
            if (btot == 0.) cycle
            btotinv = 1./btot
            tanalpha = tan(gaminv(ip)*btot*const*dtp(ip))
            tx = bx(ip)*tanalpha*btotinv
            ty = by(ip)*tanalpha*btotinv
            tz = bz(ip)*tanalpha*btotinv
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

      if (lwxytimesubs) timebpushxy = timebpushxy + wtime() - substarttime
      return
      end
c=========================================================================== 
      subroutine xpushxy(np,xp,yp,zp,uxp,uyp,uzp,gaminv,dtp)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: dtp(np)
      real(kind=8):: xp(np),yp(np),zp(np),uxp(np),uyp(np),uzp(np),gaminv(np)

c  Advance particle positions

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      do ip=1,np
        xp(ip) = xp(ip) + uxp(ip)*dtp(ip)*gaminv(ip)
        yp(ip) = yp(ip) + uyp(ip)*dtp(ip)*gaminv(ip)
        zp(ip) = zp(ip) + uzp(ip)*dtp(ip)*gaminv(ip)
      enddo

      if (lwxytimesubs) timexpushxy = timexpushxy + wtime() - substarttime
      return
      end      
c============================================================================ 
      subroutine initzpxy(pgroup)
      use ParticleGroupmodule
      use InPart
      use Picglb
      use InGenxy
      type(ParticleGroup):: pgroup

c If requested to, this routine sets the z value of all of the particles to
c zbeam.

      integer(ISZ):: is

      if (lcommonz .or. lwithez) then
        do is=1,pgroup%ns
          pgroup%zp(pgroup%ins(is):pgroup%ins(is)+pgroup%nps(is)-1) = zbeam
        enddo
      endif

      return
      end
c============================================================================ 
      subroutine initdtp(pgroup)
      use ParticleGroupmodule
      use Subtimersxy
      use Beam_acc
      use InGen
      use InGenxy
      use Particlesxy
      type(ParticleGroup):: pgroup
c Initializes dtp
      integer(ISZ):: is,i1,i2
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c     --- Give an initial value to dtp if it has not already been set or
c     --- if any values are zero.
!$OMP DO PRIVATE(i1,i2)
      do is=1,pgroup%ns
        i1 = pgroup%ins(is)
        i2 = pgroup%ins(is) + pgroup%nps(is) - 1
        if (minval(abs(pgroup%pid(i1:i2,dtpid))) == 0.) then
          if (lthick .or. vbeam == 0.) then
            pgroup%pid(i1:i2,dtpid) = dt
          else
            where (pgroup%uzp(i1:i2) .ne. 0.)
              pgroup%pid(i1:i2,dtpid) = dt*vbeam/
     &                               (pgroup%uzp(i1:i2)*pgroup%gaminv(i1:i2))
            elsewhere
              pgroup%pid(i1:i2,dtpid) = 0.
            endwhere
          endif
        endif
      enddo
!$OMP END DO

      if (lwxytimesubs) timeinitdtp = timeinitdtp + wtime() - substarttime
      return
      end
c============================================================================ 
      subroutine setdtp(np,dtp,xp,zp,zpo,uzp,bendres,bendradi,zgridprv)
      use Subtimersxy
      use InGenxy
      integer(ISZ):: np
      real(kind=8):: dtp(np),xp(np),zp(np),zpo(np),uzp(np)
      real(kind=8):: bendres,bendradi,zgridprv

c Set dt for each of the particles.  With changing vz, dt needs to be
c recalculated for each particle each timestep.

      integer(ISZ):: ip
      real(kind=8):: dsp(np)
      real(kind=8):: dtheta0,dtheta,rr
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      if (lvzchang) then

c       --- This algorithm scales dt by the ratio of the desired ds and
c       --- distance actually travelled.  Note that in a bend, the distance
c       --- traveled is calculated as an angle around the bend.

c       --- bendres is compared to 0.5 since the only possible values are 0
c       --- and 1 - this avoids potential problems with roundoff errors.
        if (bendres > 0.5 .and. lexbend) then

          dtheta0 = ds/abs(bendradi)

          do ip=1,np
            if (uzp(ip) /= 0.) then
              rr = sqrt((bendradi+xp(ip))**2 + (zp(ip) - zpo(ip))**2)
              dtheta = asin((zp(ip) - zpo(ip))/rr)
              dtp(ip) = dtp(ip)*dtheta0/dtheta
            endif
          enddo

        else

c         --- with no bends...
          if ((lcommonz .or. lwithez) .and. lzstepcorrection) then
c           --- If the particles share a common z, a correction can be done
c           --- when calculating the step size to account for errors on
c           --- previous steps. If the particle position on the previous
c           --- step is short, then the dt will be larger with this correction.
c           --- The correction tries to move the particle to the new zgrid
c           --- rather than just move the particle ds.
            dsp = ds + (zgridprv - zpo)
          else
            dsp = ds
          endif

          do ip=1,np
            if (uzp(ip) /= 0. .and. ds .ne. 0.) then
              dtp(ip) = dtp(ip)*dsp(ip)/(zp(ip) - zpo(ip))
            endif
          enddo

        endif

      endif

      if (lwxytimesubs) timesetdtp = timesetdtp + wtime() - substarttime
      return
      end      
c============================================================================ 
      subroutine getnewdtpwithe(np,dtp,uzp,gaminv,ez,q,m,ds)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: dtp(np),uzp(np),gaminv(np),ez(np)
      real(kind=8):: q,m,ds
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      dtp = 2.*ds/(sqrt((uzp*gaminv)**2 + 2.*ez*q/m*ds*gaminv) + uzp*gaminv)

      if (lwxytimesubs) timegetnewdtpwithe = timegetnewdtpwithe + wtime() - substarttime
      return
      end
c============================================================================ 
      subroutine nextbend(zbeam,zz)
      use Subtimersxy
      use InGen
      use Lattice
      use LatticeInternal
      real(kind=8):: zbeam,zz

c Get distance to get to start or end of bend, which ever is closer.

      integer(ISZ):: j
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

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

      if (lwxytimesubs) timenextbend = timenextbend + wtime() - substarttime
      return
      end
c============================================================================ 
      subroutine bendezxy(np,xp,zp,ez,bendres,bendradi,bends,bnezflag)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: xp(np),zp(np),ez(np),bendres,bendradi
      logical(ISZ):: bends,bnezflag

c  Corrects axial electric field at particle position for warped geometry
c  via multiplying by r_star/r = 1 - x/r, in a residence-corrected way;
c  at smaller radii, zones are closer together, so field is larger.

      integer(ISZ):: ip
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      if (.not. (bends.and.bnezflag) .or. bendres < 0.5) return

      do ip=1,np
         ez(ip) = ez(ip)*( 1. - bendres*xp(ip)/(bendradi + xp(ip)))
      enddo

      if (lwxytimesubs) timebendezxy = timebendezxy + wtime() - substarttime
      return
      end      
c=========================================================================== 
      subroutine exbendcorxy(np,xp,zp,uxp,uzp,zpo,ds,bendres,bendradi,lexbend)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: xp(np), zp(np), zpo(np), uxp(np), uzp(np)
      real(kind=8):: ds,bendres,bendradi
      logical(ISZ):: lexbend

c Exact bend translation:
c Translate the particle position and velocity to the new rotated
c frame (when in a bend).
c The expression for translating xp is written in a way that is correct
c for either sign of bendradi.

      integer(ISZ):: ip
      real(kind=8):: rr
      real(kind=8):: dtheta0,ty,sy,uxppr,uzppr
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c       --- bendres is compared to 0.5 since the only possible values are 0
c       --- and 1 - this avoids potential problems with roundoff errors.
      if (bendres > 0.5 .and. lexbend) then

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

      if (lwxytimesubs) timeexbendcorxy = timeexbendcorxy + wtime() - substarttime
      return
      end
c=========================================================================== 
      subroutine bendcorxy(np,xp,zp,uxp,uzp,gaminv,dtp,zpo,ds,
     &                     bendres,bendradi,lexbend)
      use Subtimersxy
      integer(ISZ):: np
      real(kind=8):: xp(np),zp(np),zpo(np),uxp(np),uzp(np),gaminv(np),dtp(np)
      real(kind=8):: ds,bendres,bendradi
      logical(ISZ):: lexbend

c Simple translation, same model as in 3-D code:
c Translate the particle position and velocity to the new rotated
c frame (when in a bend).

      integer(ISZ):: ip
      real(kind=8):: xprv,xc
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c       --- bendres is compared to 0.5 since the only possible values are 0
c       --- and 1 - this avoids potential problems with roundoff errors.
      if (bendres > 0.5 .and. .not. lexbend) then

        do ip=1,np
          xprv = xp(ip) - dtp(ip)*uxp(ip)*gaminv(ip)
          xc = 0.5*(xp(ip) + xprv)
          zp(ip) = zp(ip) + dtp(ip)*uzp(ip)*gaminv(ip)*bendres*
     &                      (bendradi/(bendradi+xc) - 1.)
        enddo

      endif

      if (lwxytimesubs) timebendcorxy = timebendcorxy + wtime() - substarttime
      return
      end      
c=========================================================================== 
      subroutine setrhoxy(rho1d,np,xp,yp,zp,zgrid,uzp,gaminv,q,wght)
      use Subtimersxy
      use GlobalVars
      use InGen
      use InGenxy
      use InGen3d
      use InMesh3d
      use Fields3d
      use Picglb3d
      integer(ISZ):: np
      real(kind=8):: q,wght,zgrid
      real(kind=8):: rho1d(0:*)
      real(kind=8):: xp(np), yp(np), zp(np), uzp(np), gaminv(np)

c  Sets charge density.
c  When using thick slice version, makes call to 3d version of setrho.
c  When using slice version, puts all charge at iz=0.  The call to
c  the fixrhoxy routine is then used to distribute the charge density
c  throughout the rho array after all of the particles have been deposited.


c  Algorithm notes: rho array is dimensioned (0:nx,0:ny,0:nzlocal) outside,
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
      integer(ISZ):: moff(0:3), indx(0:3,0:nparpgrp-1)
c     --- For "scalar" (actually partly vectorized) algorithm
      integer(ISZ):: ii(0:nparpgrp-1), jj(0:nparpgrp-1)
c     --- For both algorithms
      real(kind=8):: s(0:3,0:nparpgrp-1)
      integer(ISZ):: ipmin,nptmp

      integer(ISZ):: ip,i,j,ind0,m,ir
      real(kind=8):: g,dxi,dyi,u0,u1,v0,v1,gxfact,gyfact,f
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      if (lthick) then
c       --- Call 3d version.
        call setrho3d(rho1d,np,xp,yp,zp,zgrid,uzp,q,wght,depos)
      else

! np was made FIRSTPRIVATE to get around a bug when the expression
! np+1-ipmin was evaluating to 1-ipmin (as if np was zero).
! I don't know why it works, but it does.
!$OMP PARALLEL PRIVATE(ipmin,nptmp,f,i,j,u1,u0,v1,v0,ir,ip,ind0,indx,gyfact,gxfact,s,m,ii,jj),FIRSTPRIVATE(np)


c  Set up offset array for vectorized deposition:

      moff(0) = 0
      moff(1) = 1
      moff(2) = nx+1
      moff(3) = nx+2

c  Begin main loop over species, groups of nparpgrp particles

      g = wght * q / (dx * dy)
      dxi = 1./dx
      dyi = 1./dy
      if (l2symtry) then
c       --- The particle weight is reduced by a factor of 2 except near the
c       --- transverse boundaries.
        g = g*0.5
      elseif (l4symtry) then
c       --- The particle weight is reduced by a factor of 4 except near the
c       --- transverse boundaries.
        g = g*0.25
      endif

!$OMP DO
      do ipmin = 1,np,nparpgrp
        nptmp = min(nparpgrp, np+1-ipmin)

c--------------------------------------
c  Begin vectorized deposition loop
c--------------------------------------
      if (depos == "vector") then

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
c       --- special loop for 2-fold symmetry
        do ip = ipmin,ipmin+nptmp-1
          if (vbeamfrm .ne. 0.) then
             f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
           else
             f = 1.
           endif
           i  = (xp(ip) - xmmin) * dxi
           u1 = (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = abs(yp(ip))*dyi
           v1 = abs(yp(ip))*dyi - j
           v0 = 1. - v1
           ir = ip - ipmin
           ind0 = i + j*(nx+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           gyfact = 1.
           if (j == 0) gyfact = 2.
           s(0,ir) = u0 * v0 * f*gyfact
           s(1,ir) = u1 * v0 * f*gyfact
           s(2,ir) = u0 * v1 * f
           s(3,ir) = u1 * v1 * f
        enddo
      elseif (l4symtry) then
c       --- special loop for 4-fold symmetry
        do ip = ipmin,ipmin+nptmp-1
          if (vbeamfrm .ne. 0.) then
             f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
           else
             f = 1.
           endif
           i  = abs(xp(ip))*dxi
           u1 = abs(xp(ip))*dxi - i
           u0 = 1. - u1
           j  = abs(yp(ip))*dyi
           v1 = abs(yp(ip))*dyi - j
           v0 = 1. - v1
           ir = ip - ipmin
           ind0 = i + j*(nx+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           gxfact = 1.
           gyfact = 1.
           if (i == 0) gxfact = 2.
           if (j == 0) gyfact = 2.
           s(0,ir) = u0 * v0 * f*gxfact*gyfact
           s(1,ir) = u1 * v0 * f*gyfact
           s(2,ir) = u0 * v1 * f*gxfact
           s(3,ir) = u1 * v1 * f
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          if (vbeamfrm .ne. 0.) then
             f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
           else
             f = 1.
           endif
           i  = (xp(ip) - xmmin) * dxi
           u1 = (xp(ip) - xmmin) * dxi - i
           u0 = 1. - u1
           j  = (yp(ip) - ymmin) * dyi
           v1 = (yp(ip) - ymmin) * dyi - j
           v0 = 1. - v1
           ir = ip - ipmin
           ind0 = i + j*(nx+1)
           indx(0,ir) = ind0 + moff(0)
           indx(1,ir) = ind0 + moff(1)
           indx(2,ir) = ind0 + moff(2)
           indx(3,ir) = ind0 + moff(3)
           s(0,ir) = u0 * v0 * f
           s(1,ir) = u1 * v0 * f
           s(2,ir) = u0 * v1 * f
           s(3,ir) = u1 * v1 * f
        enddo
      endif
c     --- vectorized deposition over the 8 cells touched;
c     --- there'd be a hazard if we interchanged the loops.
!$OMP CRITICAL (CRITICAL_SETRHOXY1)
      do ir = 0,nptmp-1
        if (uzp(ir+ipmin) /= 0) then
          do m = 0, 3
            rho1d(indx(m,ir)) = rho1d(indx(m,ir)) + s(m,ir)
          enddo
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETRHOXY1)

c--------------------------------------
c  Begin scalar deposition loop
c--------------------------------------
      elseif (depos == "scalar") then

c  Begin main loop over species, groups of nparpgrp particles

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        do ip = ipmin,ipmin+nptmp-1
          if (vbeamfrm .ne. 0.) then
            f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          else
            f = 1.
          endif
          ir = ip - ipmin
          ii(ir) = (xp(ip) - xmmin) * dxi
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = abs(yp(ip))*dyi
          v1     = abs(yp(ip))*dyi - jj(ir)
          v0     = 1. - v1
          gyfact = 1.
          if (jj(ir) == 0) gyfact = 2.
          s(0,ir) = u0 * v0 * f*gyfact
          s(1,ir) = u1 * v0 * f*gyfact
          s(2,ir) = u0 * v1 * f
          s(3,ir) = u1 * v1 * f
        enddo
      elseif (l4symtry) then
        do ip = ipmin,ipmin+nptmp-1
          if (vbeamfrm .ne. 0.) then
            f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          else
            f = 1.
          endif
          ir = ip - ipmin
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
          s(0,ir) = u0 * v0 * f*gxfact*gyfact
          s(1,ir) = u1 * v0 * f*gyfact
          s(2,ir) = u0 * v1 * f*gxfact
          s(3,ir) = u1 * v1 * f
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          if (vbeamfrm .ne. 0.) then
            f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          else
            f = 1.
          endif
          ir = ip - ipmin
          ii(ir) = (xp(ip) - xmmin) * dxi
          u1     = (xp(ip) - xmmin) * dxi - ii(ir)
          u0     = 1. - u1
          jj(ir) = (yp(ip) - ymmin) * dyi
          v1     = (yp(ip) - ymmin) * dyi - jj(ir)
          v0     = 1. - v1
          s(0,ir) = u0 * v0 * f
          s(1,ir) = u1 * v0 * f
          s(2,ir) = u0 * v1 * f
          s(3,ir) = u1 * v1 * f
        enddo
      endif
c     --- scalar loop does the actual deposition
!$OMP CRITICAL (CRITICAL_SETRHOXY2)
      do ir = 0, nptmp-1
         if (uzp(ir+ipmin) /= 0) then
         rho(ii(ir)  ,jj(ir)  ,0) = rho(ii(ir)  ,jj(ir)  ,0) + s(0,ir)
         rho(ii(ir)+1,jj(ir)  ,0) = rho(ii(ir)+1,jj(ir)  ,0) + s(1,ir)
         rho(ii(ir)  ,jj(ir)+1,0) = rho(ii(ir)  ,jj(ir)+1,0) + s(2,ir)
         rho(ii(ir)+1,jj(ir)+1,0) = rho(ii(ir)+1,jj(ir)+1,0) + s(3,ir)
        endif
      enddo
!$OMP END CRITICAL (CRITICAL_SETRHOXY2)

c--------------------------------------
c  Direct deposition loop with precalculated integer conversions
c--------------------------------------
      elseif (depos == "direct1") then

c  Begin main loop over species, groups of nparpgrp particles

c     --- vectorized loop to compute indices, weights
      if (l2symtry) then
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0) cycle
          if (vbeamfrm .ne. 0.) then
            f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          else
            f = 1.
          endif
          i      = (xp(ip) - xmmin) * dxi
          u1     = (xp(ip) - xmmin) * dxi - i
          u0     = 1. - u1
          j      = abs(yp(ip))*dyi
          v1     = abs(yp(ip))*dyi - j
          v0     = 1. - v1
          gyfact = 1.
          if (j == 0) gyfact = 2.
          rho(i  ,j  ,0) = rho(i  ,j  ,0) + u0 * v0 * f*gyfact
          rho(i+1,j  ,0) = rho(i+1,j  ,0) + u1 * v0 * f*gyfact
          rho(i  ,j+1,0) = rho(i  ,j+1,0) + u0 * v1 * f
          rho(i+1,j+1,0) = rho(i+1,j+1,0) + u1 * v1 * f
        enddo
      elseif (l4symtry) then
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0) cycle
          if (vbeamfrm .ne. 0.) then
            f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          else
            f = 1.
          endif
          i      = abs(xp(ip))*dxi
          u1     = abs(xp(ip))*dxi - i
          u0     = 1. - u1
          j      = abs(yp(ip))*dyi
          v1     = abs(yp(ip))*dyi - j
          v0     = 1. - v1
          gxfact = 1.
          gyfact = 1.
          if (i == 0) gxfact = 2.
          if (j == 0) gyfact = 2.
          rho(i  ,j  ,0) = rho(i  ,j  ,0) + u0 * v0 * f*gxfact*gyfact
          rho(i+1,j  ,0) = rho(i+1,j  ,0) + u1 * v0 * f*gyfact
          rho(i  ,j+1,0) = rho(i  ,j+1,0) + u0 * v1 * f*gxfact
          rho(i+1,j+1,0) = rho(i+1,j+1,0) + u1 * v1 * f
        enddo
      else
c       --- normal loop
        do ip = ipmin,ipmin+nptmp-1
          if (uzp(ip) == 0) cycle
          if (vbeamfrm .ne. 0.) then
            f = g*vbeamfrm/dvnz(uzp(ip)*gaminv(ip))
          else
            f = 1.
          endif
          i      = (xp(ip) - xmmin) * dxi
          u1     = (xp(ip) - xmmin) * dxi - i
          u0     = 1. - u1
          j      = (yp(ip) - ymmin) * dyi
          v1     = (yp(ip) - ymmin) * dyi - j
          v0     = 1. - v1
          rho(i  ,j  ,0) = rho(i  ,j  ,0) + u0 * v0 * f
          rho(i+1,j  ,0) = rho(i+1,j  ,0) + u1 * v0 * f
          rho(i  ,j+1,0) = rho(i  ,j+1,0) + u0 * v1 * f
          rho(i+1,j+1,0) = rho(i+1,j+1,0) + u1 * v1 * f
        enddo
      endif

      endif

      enddo
!$OMP END PARALLEL

c     --- End of slice version
      endif
      if (lwxytimesubs) timesetrhoxy = timesetrhoxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine loadrhoxy(pgroup,ins_i,nps_i,is_i,lzero)
      use ParticleGroupmodule
      use Subtimersxy
      use GlobalVars
      use InGen
      use InGenxy
      use InGen3d
      use Picglb
      use Picglb3d
      use Particles,Only: wpid
      use Particlesxy
      use InMesh3d
      use Fields3d
      type(ParticleGroup):: pgroup
      integer(ISZ):: ins_i,nps_i,is_i
      logical(ISZ):: lzero

c --- This routine provides a simple call from the interpreter to load the
c --- rho array.  The value '-1' is used as a flag in the input to use
c --- all of the particles, otherwise the specified particles are loaded.

      integer(ISZ):: ins_u,nps_u
      integer(ISZ):: is1,is2
      integer(ISZ):: ip,ipmin,is
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c     --- zero rho if requested
c     if (lzero) call zeroarry(rho,(nx+1)*(ny+1)*(nzlocal+1))
      if(lzero) then
        if(solvergeom==XYgeom) then
          call reset_rzmgrid_rho()
        else
          rho = 0.
        end if
      end if

c     --- set limits on loop over species
      if (is_i == -1) then
        is1 = 1
        is2 = pgroup%ns
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
         if (ins_i == -1) ins_u = pgroup%ins(is)
         if (nps_i == -1) nps_u = pgroup%nps(is)

c        --- loop over particle blocks
         if(solvergeom==XYgeom) then
           do ipmin = ins_u, ins_u + nps_u - 1, nparpgrp
             ip = min(nparpgrp, ins_u+nps_u-ipmin)
             if(wpid==0) then
               call rhoweightrz(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                          pgroup%yp(ipmin),ip,
     &                          pgroup%sq(is)*pgroup%sw(is),
     &                          nx,ny,dx,dy,xmmin,ymmin)
             else
               call rhoweightrz_weights(pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                                  pgroup%yp(ipmin),pgroup%pid(ipmin,wpid),
     &                                  ip,pgroup%sq(is)*pgroup%sw(is),
     &                                  nx,ny,dx,dy,xmmin,ymmin)
             end if
           enddo
         else
           ip = nps_u
           ipmin = ins_u
           call setrhoxy(rho,ip,pgroup%xp(ipmin),pgroup%yp(ipmin),
     &                   pgroup%zp(ipmin),zgrid,
     &                   pgroup%uzp(ipmin),pgroup%gaminv(ipmin),
     &                   pgroup%sq(is),pgroup%sw(is))
         end if
      enddo

c     --- enforce axial periodicity if rho was zeroed
      if (lzero) call fixrhoxy(rho(0,0,0),nx,ny,nzlocal,periinz,lthick)

      if(solvergeom==XYgeom) then
c  Distribute rho for AMR solver
        call distribute_rho_rz()
c  Enforce boundary conditions
        call rhobndrz()
c  Copy charge density from frz.basegrid to w3d.rho
        call get_rho_rz(rho(0,0,0),nx,ny,1,0)
      end if

      if (lwxytimesubs) timeloadrhoxy = timeloadrhoxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine setcurrxy(nzzarr,nszarr,curr,np,zp,uzp,gaminv,q,wght,is,
     &                     zbeam,dzz,zzmin,
     &                     ns,lspeciesmoments,dz,lthick)
      use GlobalVars
      use Subtimersxy
      integer(ISZ):: nzzarr,nszarr,np,is,ns
      real(kind=8):: q,wght,zbeam,dzz,zzmin,dz
      real(kind=8):: curr(0:nzzarr,0:nszarr), zp(np), uzp(np), gaminv(np)
      logical(ISZ):: lspeciesmoments,lthick

c  Sets 1d beam current directly from particle data.
c  If running slice version, puts all data at iz=0, else if running thick
c  slice version, makes a call to the 3d version of setcurr.

      real(kind=8):: csum,g
      integer(ISZ):: ip,ipmin,js
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

      if (lthick) then
c       --- Call 3d version.
!$OMP PARALLEL PRIVATE(ip)
!$OMP DO
        do ipmin = 1,np,nparpgrp
          ip = min(nparpgrp, np+1-ipmin)
          call setcurr(nzzarr,nszarr,curr,ip,zp(ipmin),uzp(ipmin),gaminv(ipmin),
     &                 q,wght,is,zbeam,dzz,zzmin,ns,lspeciesmoments)
        enddo
!$OMP END DO
!$OMP END PARALLEL

      else

        if (lspeciesmoments) then
c         --- Check if the moments are to be calculated separately for
c         --- each species. If so, check if nszarr already has been set
c         --- appropriately. If not, set it and allocate the arrays.
c         --- If only one species, then don't have separate species data.
          if (nszarr /= ns .and. ns > 1) then
            nszarr = ns
            call gchange("Z_arrays",0)
          endif
          js = is - 1
        else
          if (nszarr /= 0) then
            nszarr = 0
            call gchange("Z_arrays",0)
          endif
          js = 0
        endif

c       --- Sum the charge times velocity over all particles.
        g = wght*q/dz
        csum = 0.
!$OMP PARALLEL DO REDUCTION(+:csum)
        do ip = 1,np
           csum = csum + gaminv(ip)*uzp(ip)
        enddo
!$OMP END PARALLEL DO
        curr(0,js) = curr(0,js) + csum*g
        if (js > 0) curr(0,0) = curr(0,0) + csum*g

      endif

      if (lwxytimesubs) timesetcurrxy = timesetcurrxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine bendfieldsolxy(bendradi)
      use Subtimersxy
      use Constant
      use InGen3d
      use Picglb3d
      use InMesh3d
      use Fields3d
      real(kind=8):: bendradi

c FFT/bend iterative field solver.

      integer(ISZ):: i,j
      real(kind=8):: r,x,dxi,phiref
      real(kind=8):: phisave(0:nx,0:ny)
      character(120):: outstr
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

!$OMP PARALLEL PRIVATE(x,r)
      dxi = 1./dx
      bndfit = 0
      bndferr = 2.*bndftol
      do while ((bndferr > bndftol) .and. (bndfit < bndfitmx))
        bndfit = bndfit + 1

c       --- save phi for error measure
        call copyarry(phi(0,0,0),phisave(0,0),(nx+1)*(ny+1))

c       --- add curvature terms to the source (rho) as a correction
!$OMP DO
        do i = 1, nx-1
          x = xmmin + i*dx
          r = 1./(bendradi + x)
          do j = 0, ny
            phi(i,j,0) = rho(i,j,0)*bendradi*r
     &               + eps0*(phisave(i+1,j) - phisave(i-1,j))*0.5*dxi*r
          enddo
        enddo
!$OMP END DO

c       --- call Cartesian field solver
        call vpxy (-1)

c       --- compute error
        bndferr = 0.
        phiref = 0.
!$OMP DO REDUCTION(max:bndferr,phiref)
        do i = 0, nx
           do j = 0, ny
              bndferr = max( bndferr, abs(phi(i,j,0)-phisave(i,j)) )
              phiref = max( phiref, phi(i,j,0) )
           enddo
        enddo
!$OMP END DO

c       --- for debug, print out the error
!$OMP SINGLE
        if (bnprflag) then
          write (outstr,9985) bndfit, bndferr
 9985     format ("Bent field iteration",i3," Rel Change = ",1pe12.4)
          call remark(outstr)
        endif
!$OMP END SINGLE NOWAIT

      end do

!$OMP END PARALLEL

c     --- if failure to converge, report the bad news to user
      if (bndferr > bndftol) then
        write (outstr,9995) bndferr, bndfit
 9995   format ("*** NONCONVERGENCE in bent field iteration",
     &          /,"Relative change = ",1pe12.4," after ",i3," iterations.")
        call remark(outstr)
      endif

      if (lwxytimesubs) timebendfieldsolxy = timebendfieldsolxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine fieldsolxy(iwhich)
      use Subtimersxy
      use Timers
      use InGen
      use InGenxy
      use InGen3d
      use InMesh3d
      use Picglb
      use Picglb3d
      use Fields3d
      use Efields3d
      use LatticeInternal
      integer(ISZ):: iwhich

c If using the 3D field solver, calls the fieldsol routine in package w3d,
c otherwise only does a FFT field solve in 2D.

      integer(ISZ):: iz
      real(kind=8):: bendres,bendradi
      real(kind=8):: timetemp,wtime
      real(kind=8):: substarttime
      if (lwxytimesubs) substarttime = wtime()
      timetemp = wtime()

c     --- Return if not doing any field solves.
      if ( fstype == -1) return

c     --- Call appropriate field solver
      if (lvp3d) then
c       --- 3D field solver
        call fieldsol3d(iwhich)
      else

c       --- 2D fieldsolver

c       --- If using approximate Ez, copy old phi into the plane iz=-1.
        if (lwithez) then
          phi(:,:,-1) = phi(:,:,0)
        endif

c       --- Get information about the closest bend.
        call getbend(1,1,zbeam,0.,1.,bendres,bendradi,0.,0.,.true.)

        if (bendres < 0.5 .or. iwhich > 0) then
c         --- if not in a bend, then make direct call to vpxy
          if (fstype == 0 .or. fstype == 1 .or. fstype == 2) then
            call copyarry(rho(0,0,0), phi(0,0,0), (nx+1)*(ny+1)*(nzlocal+1))
          endif
          call vpxy(iwhich)

        else
c         --- otherwise, include bend correction terms
          call bendfieldsolxy(bendradi)

        endif

c       --- fill rest of phi array
        if (lwithez .and. it > 0) then

c         --- Use old phi and new phi to extrapolate a value to the plate iz=+1.
c         --- This is skipped on the first step since there is no old phi.
          phi(:,:,1) = 2.*phi(:,:,0) - phi(:,:,-1)

        else

c         --- Fill rest of phi array with copy of new phi
          call copyarry(phi(0,0,0),phi(0,0,-1),(nx+1)*(ny+1))
          do iz=1,nzlocal+1
            call copyarry(phi(0,0,0),phi(0,0,iz),(nx+1)*(ny+1))
          enddo

        endif

      endif

      fstime = fstime + (wtime() - timetemp)
      if (lwxytimesubs) timefieldsolxy = timefieldsolxy + wtime() - substarttime
      return
      end
c=============================================================================
      subroutine vpxy(iwhich)
      use Subtimersxy
      use Constant
      use InGen
      use InGenxy
      use InGen3d
      use InMesh3d
      use Picglb3d
      use Picglb
      use Fields3d

      integer(ISZ):: iwhich
      real(kind=8):: substarttime,wtime
      if (lwxytimesubs) substarttime = wtime()

c  Interface to VPOIS2D and VPOIS3D using variables from database of package 3D

      if (lvp3d) then
c       --- Use 3d field solver
        call vp3d(iwhich)
      else

c       --- Otherwise, use 2d field solver

        if (iwhich==0 .or. iwhich==1) call execuser('initfieldsolver')

        if (fstype == 0) then
          call vpois2d (iwhich, phi(0,0,0), phi(0,0,0), kxsq, kysq,
     &       attx, atty, filt, xmmax-xmmin, ymmax-ymmin,
     &       nx, ny, scrtch, phi(0,0,1), 0, l2symtry,l4symtry)

        elseif (fstype == 1 .or. fstype == 2) then
c         --- Note that phi(,,-1) and phi(,,1) are passed in as
c         --- scratch space.
          call capmatxyf(iwhich,phi(0,0,0),kxsq,kysq,attx,atty,
     &             filt,xmmax-xmmin,ymmax-ymmin,nx,ny,dx,dy,xmmin,ymmin,
     &             scrtch,phi(0,0,2),phi(0,0,1),l2symtry,l4symtry)
        elseif (fstype == 7) then
c         --- Multigrid solver
          call multigridxyf(iwhich,nx,ny,dx,dy,phi(0,0,0),rho(0,0,0),
     &                      l2symtry,l4symtry,xmmin,ymmin)
c         --- AMR multigrid solver
        elseif (fstype == 10) then
          call multigridxyf2(iwhich,phi(:,:,0),rho(:,:,0),nx,ny)
        endif

      endif

      if (lwxytimesubs) timevpxy = timevpxy + wtime() - substarttime
      return
      end
c=============================================================================

